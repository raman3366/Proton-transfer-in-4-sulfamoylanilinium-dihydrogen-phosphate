#!/usr/bin/env python3
"""Gate G-Z of md/PLAN.md for one zinc set, read from the validation run X0 of that set.

The plan (stage S2, entries of 5 Oct 2026, 9:41 PM IST, and 6 Oct 2026, 3:02 AM IST):
"after minimisation and 100 ps", (a) Zn-N(sulfonamide) within 0.10 A of 1.93 A, (b) Zn-N(His)
within 0.10 A of 2.01-2.02 A, (c) each mean N-Zn-N angle within 10 deg of its value in 6RL9,
(d) an oxygen (ligand sulfonyl, water or protein) within 2.6 A of the zinc in at most 5 % of the
frames; the distance from the zinc to Glu106 OE1 is reported by name.  The four parts and their
numbers are those of GV3 of md/analysis/SPEC.md, section 11 (version 1.3), taken as printed there.

"After minimisation and 100 ps" is the first equilibration step of the run (eq01: 100 ps, heating
from 50 K to 300 K, heavy atoms of the solute held with k = 1000 kJ mol-1 nm-2).  The gate is read
on the frames of that step: 10 solute frames (10 to 100 ps) and 1 all-atom frame (100 ps), from the
series the analysis writes for the equilibration (timeseries_equil_solute.csv and
timeseries_equil_full.csv, rows with stage 1).  The means are rounded to six decimals before they
are compared, as the specification does for gate G-V.

Because the heavy atoms are held in that step, the same four parts are also given, FOR INFORMATION
and without a verdict, for the minimised structure (state_min.xml), for the frame at 100 ps alone,
for each later equilibration stage and for the last stage, which has no restraint; the zinc geometry
of the free 5 ns is part GV3 of gate G-V and is copied here from summary.json.

Usage:  gate_GZ.py RUN_DIR --out FILE.json [--analysis analysis] [--second analysis_recheck]
Nothing is written into the run folder.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path

# constants of md/analysis/SPEC.md, section 11, as printed there
WIN_ZN_N2 = (1.83, 2.03)
WIN_ZN_NHIS = (1.91, 2.12)
ANGLES = [  # column, reference in 6RL9 (deg)
    ("ang_N2_Zn_H94NE2_deg", 108.9),
    ("ang_N2_Zn_H96NE2_deg", 114.3),
    ("ang_N2_Zn_H119ND1_deg", 116.8),
    ("ang_H94NE2_Zn_H96NE2_deg", 103.8),
    ("ang_H94NE2_Zn_H119ND1_deg", 113.6),
    ("ang_H96NE2_Zn_H119ND1_deg", 98.4),
]
ANGLE_TOL = 10.0
D_HIS = ["d_Zn_H94NE2_A", "d_Zn_H96NE2_A", "d_Zn_H119ND1_A"]
FIFTH_A = 2.6
USED = ["d_Zn_N2_A"] + D_HIS + [a for a, _ in ANGLES] + [
    "fifth_sulfonylO", "d_Zn_Osulfonyl_min_A", "d_Zn_E106OE1_A", "d_Zn_E106OE2_A"]
USED_FULL = ["fifth_any", "d_Zn_O_any_min_A", "d_Zn_Owater_min_A", "d_Zn_Oprotein_min_A"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def now():
    u = datetime.now(timezone.utc)
    i = u.astimezone(timezone(timedelta(hours=5, minutes=30)))
    return {"utc": u.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "ist": i.strftime("%-d %b %Y, %-I:%M:%S %p IST")}


def rows_of(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def mean(v):
    return round(sum(v) / len(v), 6)


def sd(v):
    if len(v) < 2:
        return None
    m = sum(v) / len(v)
    return round(math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)), 6)


def parts(sol, full):
    """Parts a to d on a set of solute rows and all-atom rows (dicts of the two series)."""
    out = {"n_solute_frames": len(sol), "n_full_frames": len(full)}
    if not sol:
        out["passed"] = None
        return out
    col = lambda c: [float(r[c]) for r in sol]
    m = mean(col("d_Zn_N2_A"))
    a = {"mean_A": m, "sd_A": sd(col("d_Zn_N2_A")), "min_A": min(col("d_Zn_N2_A")),
         "max_A": max(col("d_Zn_N2_A")), "window_A": list(WIN_ZN_N2),
         "passed": WIN_ZN_N2[0] <= m <= WIN_ZN_N2[1]}
    b = {}
    for c in D_HIS:
        m = mean(col(c))
        b[c] = {"mean_A": m, "sd_A": sd(col(c)), "min_A": min(col(c)), "max_A": max(col(c)),
                "window_A": list(WIN_ZN_NHIS), "passed": WIN_ZN_NHIS[0] <= m <= WIN_ZN_NHIS[1]}
    cpart = {}
    for c, ref in ANGLES:
        m = mean(col(c))
        lo, hi = round(ref - ANGLE_TOL, 1), round(ref + ANGLE_TOL, 1)
        cpart[c] = {"mean_deg": m, "sd_deg": sd(col(c)), "min_deg": min(col(c)), "max_deg": max(col(c)),
                    "reference_6RL9_deg": ref, "deviation_deg": round(m - ref, 6),
                    "window_deg": [lo, hi], "passed": lo <= m <= hi}
    cs = sum(int(r["fifth_sulfonylO"]) for r in sol)
    d = {"limit_max_fraction": 0.05,
         "solute_frames": {"flag": "fifth_sulfonylO", "count": cs, "n": len(sol),
                           "fraction": round(cs / len(sol), 6),
                           "nearest_sulfonyl_O_min_A": min(col("d_Zn_Osulfonyl_min_A")),
                           "passed": 20 * cs <= len(sol)}}
    if full:
        cf = sum(int(r["fifth_any"]) for r in full)
        d["full_frames"] = {"flag": "fifth_any", "count": cf, "n": len(full),
                            "fraction": round(cf / len(full), 6),
                            "nearest_O_of_any_kind_min_A": min(float(r["d_Zn_O_any_min_A"]) for r in full),
                            "nearest_water_O_min_A": min(float(r["d_Zn_Owater_min_A"]) for r in full),
                            "nearest_protein_O_min_A": min(float(r["d_Zn_Oprotein_min_A"]) for r in full),
                            "passed": 20 * cf <= len(full)}
    else:
        d["full_frames"] = {"flag": "fifth_any", "n": 0, "passed": False}
    d["passed"] = d["solute_frames"]["passed"] and d["full_frames"]["passed"]
    glu = {}
    for c in ("d_Zn_E106OE1_A", "d_Zn_E106OE2_A"):
        v = col(c)
        n_in = sum(1 for x in v if x <= FIFTH_A)
        glu[c] = {"mean_A": mean(v), "min_A": min(v), "count_within": n_in, "n": len(v),
                  "fraction_within": round(n_in / len(v), 6)}
    out.update({"a_Zn_N2": a, "b_Zn_NHis": b, "c_N_Zn_N": cpart, "d_fifth_ligand": d,
                "reported_zn_glu106": dict(glu, note="reported by name; no criterion", distance_max_A=FIFTH_A)})
    out["passed"] = bool(a["passed"] and all(x["passed"] for x in b.values())
                         and all(x["passed"] for x in cpart.values()) and d["passed"])
    out["failed_parts"] = ([] if a["passed"] else ["a: d_Zn_N2_A"]) \
        + ["b: " + k for k, x in b.items() if not x["passed"]] \
        + ["c: " + k for k, x in cpart.items() if not x["passed"]] \
        + ([] if d["passed"] else ["d: fifth ligand"])
    return out


def minimised(run, sel):
    """Zinc geometry of the minimised structure, from state_min.xml (positions in nm)."""
    p = run / "state_min.xml"
    if not p.exists():
        return {"available": False}
    need = {"zn": sel["zn"], "N2": sel["lig_heavy"]["N2"], "O1": sel["lig_heavy"]["O1"],
            "O2": sel["lig_heavy"]["O2"], "H94NE2": sel["his94_NE2"], "H96NE2": sel["his96_NE2"],
            "H119ND1": sel["his119_ND1"], "E106OE1": sel["glu106_OE1"], "E106OE2": sel["glu106_OE2"]}
    want = {v: k for k, v in need.items()}
    pos = {}
    pat = re.compile(r'<Position x="([^"]+)" y="([^"]+)" z="([^"]+)"')
    i = 0
    with open(p) as f:
        for line in f:
            mt = pat.search(line)
            if mt:
                if i in want:
                    pos[want[i]] = tuple(10.0 * float(g) for g in mt.groups())
                i += 1
            elif "</Positions>" in line:
                break
    dist = lambda x, y: math.sqrt(sum((pos[x][k] - pos[y][k]) ** 2 for k in range(3)))

    def ang(x, c, y):
        u = [pos[x][k] - pos[c][k] for k in range(3)]
        v = [pos[y][k] - pos[c][k] for k in range(3)]
        cosv = sum(p * q for p, q in zip(u, v)) / math.sqrt(sum(p * p for p in u) * sum(q * q for q in v))
        return math.degrees(math.acos(max(-1.0, min(1.0, cosv))))
    r = {"available": True, "file": str(p), "sha256": sha256(p), "n_positions_read": i,
         "note": "solute atoms are never wrapped; plain distances",
         "d_Zn_N2_A": round(dist("zn", "N2"), 4)}
    for c, k in zip(D_HIS, ("H94NE2", "H96NE2", "H119ND1")):
        r[c] = round(dist("zn", k), 4)
    names = {"N2": "N2", "H94NE2": "H94NE2", "H96NE2": "H96NE2", "H119ND1": "H119ND1"}
    for c, ref in ANGLES:
        x, y = c[len("ang_"):-len("_deg")].split("_Zn_")
        r[c] = round(ang(names[x], "zn", names[y]), 3)
        r[c + "_minus_6RL9"] = round(r[c] - ref, 3)
    r["d_Zn_O1_A"] = round(dist("zn", "O1"), 4)
    r["d_Zn_O2_A"] = round(dist("zn", "O2"), 4)
    r["d_Zn_E106OE1_A"] = round(dist("zn", "E106OE1"), 4)
    r["d_Zn_E106OE2_A"] = round(dist("zn", "E106OE2"), 4)
    r["each_value_in_its_window"] = bool(
        WIN_ZN_N2[0] <= r["d_Zn_N2_A"] <= WIN_ZN_N2[1]
        and all(WIN_ZN_NHIS[0] <= r[c] <= WIN_ZN_NHIS[1] for c in D_HIS)
        and all(abs(r[c] - ref) <= ANGLE_TOL for c, ref in ANGLES)
        and min(r["d_Zn_O1_A"], r["d_Zn_O2_A"]) > FIFTH_A)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--analysis", default="analysis", help="folder of the primary analysis in the run folder")
    ap.add_argument("--second", default="analysis_recheck", help="folder of the second program")
    ap.add_argument("--test", action="store_true",
                    help="accept a run that is not X0 (a test of this program; the record says so)")
    a = ap.parse_args()
    run = a.run.resolve()
    man = json.load(open(run / "manifest.json"))
    if man.get("run_id") != "X0" and not a.test:
        sys.exit("gate G-Z is read from a validation run X0; this is run %r" % man.get("run_id"))
    sysdir = Path(man["system_dir"])
    sel = json.load(open(sysdir / "selections.json"))
    es, ef = run / a.analysis / "timeseries_equil_solute.csv", run / a.analysis / "timeseries_equil_full.csv"
    sol, full = rows_of(es), rows_of(ef)
    rec = {"schema": "cqc-sfa-md gate G-Z, version 1", "written": now(),
           "test_of_the_program_only": bool(a.test and man.get("run_id") != "X0"),
           "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
           "run": {"dir": str(run), "run_id": man["run_id"], "status": man.get("status"),
                   "system_dir": man["system_dir"], "model": man.get("build_json", {}).get("model"),
                   "parameter_set": man.get("build_json", {}).get("parameter_set"),
                   "system_prmtop_sha256": man["inputs"]["system.prmtop"]["sha256"],
                   "run_seed": man["seed"]["run_seed"]},
           "inputs": {str(es): sha256(es), str(ef): sha256(ef)},
           "criteria": {"source": "md/analysis/SPEC.md, section 11, GV3 parts a to d, as printed",
                        "a_window_A": list(WIN_ZN_N2), "b_window_A": list(WIN_ZN_NHIS),
                        "c_references_6RL9_deg": dict(ANGLES), "c_tolerance_deg": ANGLE_TOL,
                        "d": "an oxygen within 2.6 A of the zinc in at most 5 % of the frames, in the solute "
                             "frames (sulfonyl O) and in the all-atom frames (any O), and at least one all-atom frame"}}
    # the two programs on the frames used
    s2, f2 = run / a.second / "timeseries_equil_solute.csv", run / a.second / "timeseries_equil_full.csv"
    if s2.exists() and f2.exists():
        sol2, full2 = rows_of(s2), rows_of(f2)
        worst, n = 0.0, 0
        same_rows = len(sol2) == len(sol) and len(full2) == len(full)
        if same_rows:
            for r1, r2 in zip(sol, sol2):
                if r1["stage"] != "1":
                    continue
                for c in USED:
                    worst = max(worst, abs(float(r1[c]) - float(r2[c]))); n += 1
            for r1, r2 in zip(full, full2):
                if r1["stage"] != "1":
                    continue
                for c in USED_FULL:
                    worst = max(worst, abs(float(r1[c]) - float(r2[c]))); n += 1
        rec["second_program"] = {"files": {str(s2): sha256(s2), str(f2): sha256(f2)},
                                 "same_number_of_rows": same_rows, "values_compared_on_the_gate_frames": n,
                                 "largest_difference": worst,
                                 "files_identical_byte_for_byte": sha256(s2) == sha256(es) and sha256(f2) == sha256(ef),
                                 "agree": bool(same_rows and n > 0 and worst == 0.0)}
    else:
        rec["second_program"] = {"available": False, "agree": None}
    st1s = [r for r in sol if r["stage"] == "1" and r["extension"] == "0"]
    st1f = [r for r in full if r["stage"] == "1" and r["extension"] == "0"]
    e1 = json.load(open(run / "rec_eq01.json"))
    rmin = json.load(open(run / "rec_min.json"))
    rec["frames_of_the_gate"] = {
        "step": "eq01", "settings": e1["settings"],
        "time_after_minimisation_ps": [e1["time_start_ps"], e1["time_end_ps"]],
        "solute_frames": {"n": len(st1s), "frame_numbers_of_the_equilibration_series": [int(r["frame"]) for r in st1s],
                          "time_ns_of_the_series": [float(r["time_ns"]) for r in st1s],
                          "interval_ps": 10.0, "file": "solute_eq01.dcd",
                          "sha256": e1["files"]["solute_eq01.dcd"]["sha256"]},
        "full_frames": {"n": len(st1f), "frame_numbers_of_the_equilibration_series": [int(r["frame"]) for r in st1f],
                        "time_ns_of_the_series": [float(r["time_ns"]) for r in st1f],
                        "interval_ps": 100.0, "file": "full_eq01.dcd",
                        "sha256": e1["files"]["full_eq01.dcd"]["sha256"]},
        "temperature_at_the_end_K": e1["final"]["temperature_K"],
        "note": "time zero of the series is the end of the equilibration, so the first 100 ps after the "
                "minimisation are -1.99 to -1.90 ns; the heavy atoms of the solute (protein, zinc, ligand) are "
                "held in this step to the coordinates of system.inpcrd, for start X the crystal positions "
                "(U = 0.5 k d^2, k = 1000 kJ mol-1 nm-2 = 2.39 kcal mol-1 A-2): far weaker than a zinc bond term, "
                "but not negligible beside an angle term at the zinc, so the angles of this step lean towards "
                "the crystal; the stages from the fourth on hold no atom of the zinc site"}
    gate = parts(st1s, st1f)
    rec["gate_GZ"] = gate
    rec["gate_GZ"]["verdict"] = "PASSED" if gate["passed"] else "FAILED"
    info = {"note": "for information; no verdict of the gate is taken from these"}
    info["minimisation"] = {"lbfgs_iterations_per_pass": rmin["lbfgs_iterations_per_pass"],
                            "restrained_atoms_displacement_A": rmin["restrained_atoms_displacement_A"],
                            "monitor_of_the_driver": {k: v for k, v in rmin["monitor"].items() if k.startswith("zn_")},
                            "structure": minimised(run, sel)}
    info["frame_at_100_ps_alone"] = parts(st1s[-1:], st1f[-1:])
    stages = {}
    for st in sorted({r["stage"] for r in sol}, key=int):
        ss = [r for r in sol if r["stage"] == st]
        ff = [r for r in full if r["stage"] == st]
        p = parts(ss, ff)
        try:
            p["settings"] = json.load(open(run / ("rec_eq%02d.json" % int(st))))["settings"]
        except Exception:
            pass
        p["includes_extension_steps"] = any(r["extension"] != "0" for r in ss)
        stages["stage_%s" % st] = p
    info["by_equilibration_stage"] = stages
    summ = run / a.analysis / "summary.json"
    if summ.exists():
        g = json.load(open(summ)).get("gate_GV", {})
        info["production_GV3_of_gate_GV"] = {"file": str(summ), "sha256": sha256(summ), "GV3_zinc": g.get("GV3_zinc")}
    rec["for_information"] = info
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(rec, f, indent=1)
        f.write("\n")
    # table
    print("gate G-Z, %s, run %s (%s)" % (rec["run"]["model"], man["run_id"], run))
    print("frames of the gate: eq01, %d solute frames (10 to 100 ps after the minimisation), %d all-atom frame(s)"
          % (len(st1s), len(st1f)))
    print("%-28s %10s %10s %10s %10s  %-17s %s" % ("quantity", "minimised", "mean eq01", "at 100 ps", "last stage", "window", ""))
    ms = info["minimisation"]["structure"]
    last = stages[sorted(stages, key=lambda k: int(k.split("_")[1]))[-1]]
    f100 = info["frame_at_100_ps_alone"]

    def line(name, key, grp, unit):
        g1 = gate[grp] if key is None else gate[grp][key]
        l1 = last[grp] if key is None else last[grp][key]
        h1 = f100[grp] if key is None else f100[grp][key]
        col = key or "d_Zn_N2_A"
        mv = ms.get(col, float("nan")) if ms.get("available") else float("nan")
        w = g1["window_" + unit]
        print("%-28s %10.3f %10.3f %10.3f %10.3f  %7.2f to %-7.2f %s" % (
            name, mv, g1["mean_" + unit], h1["mean_" + unit], l1["mean_" + unit], w[0], w[1],
            "ok" if g1["passed"] else "OUTSIDE"))
    line("Zn-N2 (A)", None, "a_Zn_N2", "A")
    for c in D_HIS:
        line(c, c, "b_Zn_NHis", "A")
    for c, _ in ANGLES:
        line(c, c, "c_N_Zn_N", "deg")
    d = gate["d_fifth_ligand"]
    print("fifth ligand: sulfonyl O within 2.6 A in %d of %d solute frames (nearest %.3f A); any O in %s of %s all-atom frames: %s"
          % (d["solute_frames"]["count"], d["solute_frames"]["n"], d["solute_frames"]["nearest_sulfonyl_O_min_A"],
             d["full_frames"].get("count"), d["full_frames"]["n"], "ok" if d["passed"] else "NOT MET"))
    gl = gate["reported_zn_glu106"]["d_Zn_E106OE1_A"]
    print("Zn...Glu106 OE1: mean %.3f A, smallest %.3f A, within 2.6 A in %d of %d frames (reported, no criterion)"
          % (gl["mean_A"], gl["min_A"], gl["count_within"], gl["n"]))
    print("the two programs on these frames: %s" % rec["second_program"])
    print("GATE G-Z: %s%s" % (rec["gate_GZ"]["verdict"],
                              "" if gate["passed"] else "  (" + "; ".join(gate["failed_parts"]) + ")"))
    print("written:", a.out)
    return 0 if gate["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
