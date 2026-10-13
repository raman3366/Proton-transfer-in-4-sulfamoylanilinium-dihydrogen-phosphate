#!/usr/bin/env python3
"""Brings the result files of the campaign into md/results/ (so that the numbers are in git) and derives the
few numbers that md/RESULTS.md needs and that no result file holds as such.

1. Copies, byte for byte, with the sha256 of source and copy compared:
     campaign tables and records of both programs, the comparison report, the log, the times and the exit
     status of the gate of stage S6, and for each of the 13 analysed runs (9 production, 2 sensitivity,
     2 validation) the summary.json of both programs and the small series of the first program
     (equilibration, all-atom frames, RMSF). The solute series of the production (about 3.4 MB per run and
     program) and the per-run figures of the first program are NOT copied; they are listed with their
     checksums and stay in /root/cqc-sfa-md/runs.
   -> md/results/SOURCES.json
2. Health of the runs, from the driver's manifest and its state logs (log_<step>.csv):
   -> md/results/run_health.json
3. What happened in each run around the release of the ligand and afterwards, counted from the time series
   of both programs (the two must give the same counts, or the script stops):
   -> md/results/run_events.json

Nothing is computed from coordinates and no result file under /root/cqc-sfa-md/runs is changed.

    /root/cqc-sfa-tools/md/bin/python md/production/collect_results.py
"""
import csv
import datetime
import hashlib
import json
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MD = os.path.dirname(HERE)
OUT = os.path.join(REPO_MD, "results")
BULK = "/root/cqc-sfa-md"
RUNS = os.path.join(BULK, "runs")
LOGS = os.path.join(BULK, "analysis_logs")
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))

PRODUCTION = ["bonded_2e/" + r for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")]
SENSITIVITY = ["bonded_2s/R1", "bonded_2e_n2fit/R1"]
VALIDATION = ["bonded_2e/X0", "bonded_2s/X0"]
ALL = PRODUCTION + SENSITIVITY + VALIDATION
SMALL_SERIES = ["timeseries_equil_solute.csv", "timeseries_equil_full.csv", "timeseries_full.csv", "rmsf_ca.csv"]
DWELL = 100


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def J(path):
    with open(path) as fh:
        return json.load(fh)


def rows_of(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def copy(src, rel, record):
    dst = os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    a, b = sha256_file(src), sha256_file(dst)
    if a != b:
        raise SystemExit("copy of %s differs from its source" % src)
    record.append({"copy": "md/results/" + rel, "source": src, "bytes": os.path.getsize(dst), "sha256": a,
                   "source_written": datetime.datetime.fromtimestamp(os.path.getmtime(src), IST)
                   .strftime("%-d %b %Y, %-I:%M:%S %p IST")})


def mean_sd(v):
    n = len(v)
    m = sum(v) / n
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / n)


# ------------------------------------------------------------------------------------------------ 1. copies
def collect_copies():
    copied, left_out = [], []
    for f in ("campaign_table.md", "campaign_table.csv", "campaign.json"):
        copy(os.path.join(RUNS, "analysis_campaign", f), "campaign/first_program/" + f, copied)
    for f in ("campaign_table.csv", "campaign.json", "recheck_report.md", "recheck_report.json"):
        copy(os.path.join(RUNS, "analysis_campaign_recheck", f), "campaign/second_program/" + f, copied)
    for f in sorted(os.listdir(os.path.join(RUNS, "analysis_campaign"))):
        if f.endswith(".png"):
            p = os.path.join(RUNS, "analysis_campaign", f)
            left_out.append({"file": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                             "why": "figure of the first program; md/results/fig_*.png show the same series"})
    for f in ("gate_S6.log", "gate_S6.times", "gate_S6.exit_status", "run_gate_S6.sh"):
        copy(os.path.join(LOGS, f), "gate/" + f, copied)
    for lab in ALL:
        d = os.path.join(RUNS, lab)
        copy(os.path.join(d, "analysis", "summary.json"), "runs/%s/first_program_summary.json" % lab, copied)
        copy(os.path.join(d, "analysis_recheck", "summary.json"), "runs/%s/second_program_summary.json" % lab, copied)
        for f in SMALL_SERIES:
            copy(os.path.join(d, "analysis", f), "runs/%s/first_program_%s" % (lab, f), copied)
            p = os.path.join(d, "analysis_recheck", f)
            left_out.append({"file": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                             "why": "second program's copy of a small series; the same numbers within the "
                                    "tolerance of the gate (md/results/campaign/second_program/recheck_report.json)"})
        for sub in ("analysis", "analysis_recheck"):
            p = os.path.join(d, sub, "timeseries_solute.csv")
            left_out.append({"file": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                             "why": "solute series of the production, one line per 10 ps; larger than 1 MB"
                             if os.path.getsize(p) > 1e6 else
                             "solute series of the production (validation run, 500 lines); left out like those of the long runs"})
        for f in sorted(os.listdir(os.path.join(d, "analysis"))):
            if f.endswith(".png"):
                p = os.path.join(d, "analysis", f)
                left_out.append({"file": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                                 "why": "per-run figure of the first program"})
        for f in ("recheck_aux.json",):
            p = os.path.join(d, "analysis_recheck", f)
            if os.path.isfile(p):
                left_out.append({"file": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                                 "why": "auxiliary record of the second program (timing and file list)"})
    rec = {"what": "Result files of the MD campaign copied into md/results/, byte for byte, and those left where "
                   "they are. Written by md/production/collect_results.py.",
           "written": datetime.datetime.now(IST).strftime("%-d %b %Y, %-I:%M %p IST"),
           "n_copied": len(copied), "bytes_copied": sum(c["bytes"] for c in copied),
           "n_left_out": len(left_out), "bytes_left_out": sum(c["bytes"] for c in left_out),
           "copied": copied, "left_out": left_out}
    with open(os.path.join(OUT, "SOURCES.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
        fh.write("\n")
    return rec


# ------------------------------------------------------------------------------------------------ 2. health
def log_stats(path):
    rows = rows_of(path)
    t = [float(r["temperature_K"]) for r in rows]
    d = [float(r["density_g_per_cm3"]) for r in rows]
    e = [float(r["potential_energy_kj_per_mol"]) for r in rows]
    bad = sum(1 for v in t + d + e if not math.isfinite(v))
    return {"rows": len(rows), "T": t, "rho": d, "non_finite": bad,
            "thermostat_K": sorted({float(r["thermostat_temperature_K"]) for r in rows})}


def health_of(lab):
    d = os.path.join(RUNS, lab)
    man = J(os.path.join(d, "manifest.json"))
    summ = J(os.path.join(d, "analysis", "summary.json"))
    out = {"run": lab, "manifest": os.path.join(d, "manifest.json"), "manifest_sha256": sha256_file(os.path.join(d, "manifest.json"))}
    inv = man["invocations"]
    out["driver"] = {
        "status": man["status"], "all_steps_complete": man["all_steps_complete"], "test_only": man["test_only"],
        "equil_scale": man["equil_scale"], "run_id": man["run_id"], "run_seed": man["seed"]["run_seed"],
        "system_dir": man["system_dir"], "n_atoms": man["system"]["n_atoms"],
        "n_water_molecules": man["system"]["n_water_molecules"],
        "crystal_water_O_listed_in_selections": man["system"].get("crystal_water_O_listed_in_selections"),
        "script_sha256": man["script"]["sha256"],
        "protocol_sha256": man["inputs"]["protocol"].get("sha256"),
        "inputs_sha256": {k: man["inputs"][k]["sha256"] for k in ("system.prmtop", "system.inpcrd")},
        "n_invocations": len(inv),
        "invocations": [{"started_ist": i["started"]["ist"], "ended_ist": (i.get("ended") or {}).get("ist"),
                         "exit_code": i.get("exit_code"), "outcome": i.get("outcome"),
                         "host": (i.get("host") or {}).get("hostname") if isinstance(i.get("host"), dict) else i.get("host")}
                        for i in inv],
        "platform": man["minimisation"]["platform"]["name"],
        "gpu": man["minimisation"]["platform"]["properties_actual"].get("DeviceName"),
        "precision": man["minimisation"]["platform"]["properties_actual"].get("Precision"),
        "openmm": man["minimisation"]["openmm"],
        "minimisation_energy_change_kj_per_mol": man["minimisation"]["energy_change_kj_per_mol"],
        "equilibration": {k: man["equilibration"][k] for k in ("complete", "n_stages", "n_stages_complete",
                                                              "n_extension_steps", "n_extension_steps_complete",
                                                              "completed_ps")},
        "production": {k: man["production"][k] for k in ("target_ns", "segment_ns", "n_segments_target",
                                                        "n_segments_complete", "completed_ns", "complete",
                                                        "frames_solute", "frames_full", "wall_s", "ns_per_day")},
    }
    b = man["build_json"]
    out["build"] = {"model": b["model"], "model_kind": b["model_kind"], "parameter_set": b["parameter_set"],
                    "start": b["start"], "status_note": b["status_note"],
                    "terms_at_N2": dict({k: b["terms_at_N2"].get(k) for k in ("source", "given_file", "first_line_of_the_file",
                                                                            "validated_by_gate_GR", "sha256")},
                                        class_of_this_set=b["terms_at_N2"]["holds"]["class_of_this_set"],
                                        entries_beyond_the_terms_as_built=b["terms_at_N2"]["holds"]["entries_beyond_the_terms_as_built"]),
                    "composition": b["composition"],
                    "geometry_as_built": b["geometry_as_built"]["values"],
                    "script_sha256": b["script_sha256"]}
    hc = man["hydration_check"]
    out["hydration_check"] = {"requested": hc["requested"], "applies": hc.get("applies"), "skipped": hc.get("skipped"),
                              "settings": {k: hc["settings"][k] for k in ("after_step", "site_atoms", "water_cutoff_A",
                                                                         "applies_if_ligand_mean_distance_above_A",
                                                                         "max_extensions")},
                              "reference_guard": hc.get("reference_guard"),
                              "counts": hc.get("counts"), "n_extensions": hc.get("n_extensions"),
                              "final_decision": hc.get("final_decision"),
                              "site_without_water_when_the_solute_was_released":
                                  hc.get("site_without_water_when_the_solute_was_released")}
    stages = []
    for s in man["equilibration"]["stages"]:
        st = s["settings"]
        ls = log_stats(os.path.join(d, "log_%s.csv" % s["tag"]))
        mt, _ = mean_sd(ls["T"])
        mr, _ = mean_sd(ls["rho"])
        stages.append({"tag": s["tag"], "name": st["name"], "ensemble": st["ensemble"], "timestep_fs": st["timestep_fs"],
                       "length_ps": st["length_ps"], "time_start_ps": s["time_start_ps"], "time_end_ps": s["time_end_ps"],
                       "restraint_atoms": st["restraint_atoms"], "restraint_k_kj_per_mol_nm2": st["restraint_k_kj_per_mol_nm2"],
                       "n_restrained_atoms": st["n_restrained_atoms"], "log_rows": ls["rows"],
                       "log_mean_temperature_K": mt, "log_mean_density_g_per_cm3": mr,
                       "log_non_finite": ls["non_finite"]})
    out["equilibration_stages"] = stages
    # production: every row of every state log
    T, RHO, seg_T, seg_rho, rows, bad, thermo, steps = [], [], [], [], 0, 0, set(), set()
    for s in man["production"]["segments"]:
        ls = log_stats(os.path.join(d, "log_%s.csv" % s["tag"]))
        T += ls["T"]
        RHO += ls["rho"]
        rows += ls["rows"]
        bad += ls["non_finite"]
        thermo |= set(ls["thermostat_K"])
        steps.add(s["settings"]["timestep_fs"])
        seg_T.append(mean_sd(ls["T"])[0])
        seg_rho.append(mean_sd(ls["rho"])[0])
        # the driver's own mean of the same log
        if abs(s["mean_over_state_log"]["temperature_K"] - seg_T[-1]) > 0.006:
            raise SystemExit("%s %s: mean temperature of the log %.3f, driver's record %.2f"
                             % (lab, s["tag"], seg_T[-1], s["mean_over_state_log"]["temperature_K"]))
    mt, st_ = mean_sd(T)
    mr, sr = mean_sd(RHO)
    half = len(T) // 2
    out["production_state_log"] = {
        "files": "log_segNNNN.csv of the run folder, one row per 10 ps", "rows": rows, "non_finite_values": bad,
        "timestep_fs": sorted(steps), "thermostat_K": sorted(thermo),
        "precision": "means and standard deviations are written as computed (no rounding here: a number rounded "
                     "here and again when it is printed can come out one unit of the last digit off)",
        "temperature_K": {"mean": mt, "sd": st_, "min": min(T), "max": max(T),
                          "mean_first_half": mean_sd(T[:half])[0], "mean_last_half": mean_sd(T[half:])[0],
                          "segment_means_min": min(seg_T), "segment_means_max": max(seg_T)},
        "density_g_per_cm3": {"mean": mr, "sd": sr, "min": min(RHO), "max": max(RHO),
                              "mean_first_half": mean_sd(RHO[:half])[0],
                              "mean_last_half": mean_sd(RHO[half:])[0],
                              "segment_means_min": min(seg_rho), "segment_means_max": max(seg_rho)}}
    # checksums: what the first program hashed when it read the trajectories against the driver's own record
    drv = {}
    for s in man["equilibration"]["stages"] + man["production"]["segments"]:
        for f, v in s["files"].items():
            drv[f] = v["sha256"]
    ins = summ["inputs_sha256"]
    dcd = sorted(f for f in ins if f.endswith(".dcd"))
    same = [f for f in dcd if drv.get(f) == ins[f]]
    out["trajectory_checksums"] = {"what": "sha256 of every DCD file as the first program read it "
                                           "(summary.json, inputs_sha256) against the driver's record in the manifest",
                                   "n_dcd_files_read": len(dcd), "n_equal_to_the_driver_record": len(same),
                                   "not_equal_or_not_recorded": sorted(set(dcd) - set(same))}
    # the driver's console log
    with open(os.path.join(d, "run_md.log")) as fh:
        lines = fh.read().splitlines()
    words = ("warn", "error", "fail", "traceback", "nan ", "not finite")
    out["driver_log"] = {"file": os.path.join(d, "run_md.log"), "lines": len(lines),
                         "lines_with_warning_error_or_failure": [ln for ln in lines if any(w in ln.lower() for w in words)],
                         "last_line": lines[-1] if lines else None}
    marker = os.path.join(d, "SENSITIVITY_RUN.txt")
    out["sensitivity_marker"] = open(marker).read().strip() if os.path.isfile(marker) else None
    return out


# ------------------------------------------------------------------------------------------------ 3. events
def events_of(lab, conv_rec):
    d = os.path.join(RUNS, lab)
    per = {}
    for sub in ("analysis", "analysis_recheck"):
        summ = J(os.path.join(d, sub, "summary.json"))
        start_state = summ["orientation"]["start_state"]
        opp = summ["orientation"]["opposite_state"]
        sign0 = summ["n2_configuration"]["sign_of_start"]
        eq = rows_of(os.path.join(d, sub, "timeseries_equil_solute.csv"))
        pr = rows_of(os.path.join(d, sub, "timeseries_solute.csv"))
        eqf = rows_of(os.path.join(d, sub, "timeseries_equil_full.csv"))
        prf = rows_of(os.path.join(d, sub, "timeseries_full.csv"))
        held = [r for r in eq if int(r["stage"]) < 4]
        rel = [r for r in eq if int(r["stage"]) >= 4 and int(r["extension"]) == 0]
        t_rel = float(rel[0]["time_ns"]) - 0.01

        def from_release_ps(r):
            return round((float(r["time_ns"]) - t_rel) * 1000.0, 3)

        def counts(rs):
            return {"n": len(rs), "C": sum(1 for r in rs if r["state"] == "C"), "R": sum(1 for r in rs if r["state"] == "R"),
                    "I": sum(1 for r in rs if r["state"] == "I"),
                    "n2_positive": sum(1 for r in rs if int(r["n2_sign"]) > 0),
                    "n2_negative": sum(1 for r in rs if int(r["n2_sign"]) < 0),
                    "hb_N2H_T199OG1": sum(int(r["hb_N2H_T199OG1"]) for r in rs)}

        def first(rs, cond):
            for k, r in enumerate(rs):
                if cond(r):
                    return {"released_frame": k, "time_from_release_ps": from_release_ps(r), "stage": int(r["stage"]),
                            "state": r["state"], "n2_sign": int(r["n2_sign"]), "n2_volume": float(r["n2_volume"]),
                            "hb_N2H_T199OG1": int(r["hb_N2H_T199OG1"]), "lig_rmsd_A": float(r["lig_rmsd_A"])}
            return None

        e = {"state_of_start": start_state, "opposite_state": opp, "n2_sign_of_start": sign0,
             "held_part": counts(held), "released_equilibration": counts(rel),
             "released_part_ps": round(-t_rel * 1000.0, 3),
             "first_released_frame": first(rel, lambda r: True),
             "first_released_frame_n2_sign_not_of_start": first(rel, lambda r: int(r["n2_sign"]) != sign0),
             "first_released_frame_not_in_state_of_start": first(rel, lambda r: r["state"] != start_state),
             "first_released_frame_in_opposite_state": first(rel, lambda r: r["state"] == opp)}
        ks = [k for k, r in enumerate(rel) if r["state"] == start_state]
        e["last_released_frame_in_state_of_start"] = (
            {"released_frame": ks[-1], "time_from_release_ps": from_release_ps(rel[ks[-1]]),
             "stage": int(rel[ks[-1]]["stage"]), "n2_sign": int(rel[ks[-1]]["n2_sign"])} if ks else None)
        e["released_frames_states"] = "".join(r["state"] for r in rel)
        e["released_frames_n2_sign"] = "".join("+" if int(r["n2_sign"]) > 0 else "-" if int(r["n2_sign"]) < 0 else "0"
                                               for r in rel)
        e["released_frames_hb_N2H_T199OG1"] = "".join(r["hb_N2H_T199OG1"] for r in rel)
        conv = conv_rec["programs"][sub]["conversion"]
        if conv and conv["in_equilibration"]:
            k = conv["joined_index"]
            e["conversion_released_frame"] = k
            e["released_before_the_conversion"] = counts(rel[:k])
            e["released_from_the_conversion_on"] = counts(rel[k:])
            k0 = e["first_released_frame_n2_sign_not_of_start"]["released_frame"]
            win = rel[max(0, min(k0, k) - 2): k + 3]
            e["frames_around_the_conversion"] = [
                {"time_from_release_ps": from_release_ps(r), "state": r["state"], "n2_volume": float(r["n2_volume"]),
                 "hb_N2H_T199OG1": int(r["hb_N2H_T199OG1"]), "d_N2_T199OG1_A": float(r["d_N2_T199OG1_A"]),
                 "ang_N2_HN2_T199OG1_deg": float(r["ang_N2_HN2_T199OG1_deg"]), "lig_rmsd_A": float(r["lig_rmsd_A"]),
                 "d1_A": float(r["d1_N1_F131CZ_A"]), "d2_A": float(r["d2_N1_N67OD1_A"])} for r in win]
            a, b = sorted((k0, k))
            e["hb_from_first_sign_change_to_the_conversion_frame"] = counts(rel[a: b + 1])
            e["stage_of_the_conversion"] = conv["stage"]
        e["production"] = counts(pr)
        # added after the audits of 8 Oct 2026: the state (d1, d2) against the sign at N2, the signed volume,
        # the changes of sign, and the ligand by state
        joined = rel + pr

        def state_and_sign(rs):
            return {"n": len(rs),
                    "R": sum(1 for r in rs if r["state"] == "R"),
                    "R_with_sign_positive": sum(1 for r in rs if r["state"] == "R" and int(r["n2_sign"]) > 0),
                    "R_with_sign_negative": sum(1 for r in rs if r["state"] == "R" and int(r["n2_sign"]) < 0),
                    "C": sum(1 for r in rs if r["state"] == "C"),
                    "C_with_sign_negative": sum(1 for r in rs if r["state"] == "C" and int(r["n2_sign"]) < 0),
                    "I": sum(1 for r in rs if r["state"] == "I"),
                    "I_with_sign_negative": sum(1 for r in rs if r["state"] == "I" and int(r["n2_sign"]) < 0)}

        def volume(rs):
            v = [float(r["n2_volume"]) for r in rs]
            sg = [int(r["n2_sign"]) for r in rs]
            return {"n": len(v), "min": min(v), "max": max(v), "n_abs_below_0.2": sum(1 for x in v if abs(x) < 0.2),
                    "n_changes_of_sign": sum(1 for a_, b_ in zip(sg, sg[1:]) if a_ != b_)}

        e["state_and_n2_sign"] = {"production": state_and_sign(pr), "from_the_release": state_and_sign(joined)}
        e["n2_volume"] = {"production": volume(pr), "from_the_release": volume(joined)}
        by = {}
        for stt in ("C", "R", "I"):
            rs = [r for r in pr if r["state"] == stt]
            by[stt] = {"n": len(rs),
                       "lig_rmsd_mean_A": sum(float(r["lig_rmsd_A"]) for r in rs) / len(rs) if rs else None,
                       "tau_fold_mean_deg": sum(float(r["tau_fold_deg"]) for r in rs) / len(rs) if rs else None}
        e["production_by_state"] = by
        # all-atom frames of the equilibration with the state of the solute frame of the same time
        smap = {round(float(r["time_ns"]), 6): r for r in eq}
        e["water_in_the_place_of_the_ring_equilibration"] = [
            {"time_ns": float(r["time_ns"]), "stage": int(r["stage"]), "extension": int(r["extension"]),
             "n_water_ring_place": int(r["n_water_ring_place"]),
             "d_water_ring_place_min_A": float(r["d_water_ring_place_min_A"]),
             "state_of_the_solute_frame": smap[round(float(r["time_ns"]), 6)]["state"]} for r in eqf]
        n = [int(r["n_water_ring_place"]) for r in prf]
        e["water_in_the_place_of_the_ring_production"] = {
            "n_full_frames": len(n), "frames_with_one_or_more": sum(1 for v in n if v >= 1), "largest_count": max(n),
            "frames_with_one_or_more_while_state_C": sum(1 for r in prf if int(r["n_water_ring_place"]) >= 1 and
                                                         pr[10 * (int(r["frame"]) + 1) - 1]["state"] == "C"),
            "frames_with_one_or_more_while_state_R": sum(1 for r in prf if int(r["n_water_ring_place"]) >= 1 and
                                                         pr[10 * (int(r["frame"]) + 1) - 1]["state"] == "R"),
            "frames_with_one_or_more_while_state_I": sum(1 for r in prf if int(r["n_water_ring_place"]) >= 1 and
                                                         pr[10 * (int(r["frame"]) + 1) - 1]["state"] == "I")}
        # visits of the production by assigned state (the column state_assigned of the programs)
        a = [r["state_assigned"] for r in pr]
        vis = []
        i = 0
        while i < len(a):
            j = i
            while j + 1 < len(a) and a[j + 1] == a[i]:
                j += 1
            vis.append((a[i], i, j - i + 1))
            i = j + 1
        byst = {}
        for st in ("C", "R"):
            v = [x for x in vis if x[0] == st]
            longest = max(v, key=lambda x: x[2]) if v else None
            byst[st] = {"n_visits": len(v), "n_long_visits": sum(1 for x in v if x[2] >= DWELL),
                        "frames_assigned": sum(x[2] for x in v),
                        "longest_visit_frames": longest[2] if longest else 0,
                        "longest_visit_first_frame": longest[1] if longest else None,
                        "longest_visit_starts_ns": round((longest[1] + 1) * 0.01, 2) if longest else None}
        e["production_visits_by_assigned_state"] = byst
        blocks = []
        for b0 in range(0, len(pr), 2000):
            blk = pr[b0: b0 + 2000]
            c = counts(blk)
            c["from_ns"], c["to_ns"] = b0 * 0.01, (b0 + len(blk)) * 0.01
            c["lig_rmsd_mean_A"] = round(sum(float(r["lig_rmsd_A"]) for r in blk) / len(blk), 4)
            blocks.append(c)
        e["production_blocks_of_20_ns"] = blocks
        per[sub] = e

    def strip(x):  # values that two correct programs may write differently in the last decimal
        if isinstance(x, dict):
            return {k: strip(v) for k, v in x.items()}
        if isinstance(x, list):
            return [strip(v) for v in x]
        if isinstance(x, float):
            return round(x, 2)
        return x
    if strip(per["analysis"]) != strip(per["analysis_recheck"]):
        raise SystemExit("%s: the series of the two programs give different counts" % lab)
    out = per["analysis"]
    out["run"] = lab
    out["counted_from"] = "analysis/ (first program); the same counts from analysis_recheck/ (second program)"
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    src = collect_copies()
    print("copied %d files (%.2f MB); left out %d files (%.1f MB), listed in md/results/SOURCES.json"
          % (src["n_copied"], src["bytes_copied"] / 1e6, src["n_left_out"], src["bytes_left_out"] / 1e6))
    # the guard's record of the end of the work (read only; nothing under status/ is written)
    fv_path = os.path.join(BULK, "status", "final_verification.json")
    fv = J(fv_path)
    final = {"file": fv_path, "sha256": sha256_file(fv_path), "ist": fv["ist"], "at": fv["at"], "why": fv["why"],
             "verified": fv["verified"], "not_verified_because": fv["not_verified_because"], "checked": fv["checked"],
             "tries": fv["tries"],
             "driver_checksums": {k: fv["driver_checksums"][k] for k in ("runs", "steps", "files_ok", "bytes_hashed",
                                                                        "missing", "bad", "ok")},
             "pi_against_machine": {k: fv["pi_against_machine"][k] for k in ("files", "same", "differ", "missing",
                                                                            "not_hashed", "live_state_not_compared", "ok")},
             "left_on_machine": fv["left_on_machine"]}
    health = {"what": "Health of the runs from the driver's manifest and state logs. Written by "
                      "md/production/collect_results.py.",
              "written": datetime.datetime.now(IST).strftime("%-d %b %Y, %-I:%M %p IST"),
              "final_verification": final,
              "runs": [health_of(lab) for lab in ALL]}
    with open(os.path.join(OUT, "run_health.json"), "w") as fh:
        json.dump(health, fh, indent=1)
        fh.write("\n")
    ct = J(os.path.join(OUT, "conversion_times.json"))
    cmap = {r["label"]: r for r in ct["runs"]}
    events = {"what": "Counts around the release of the ligand and in the production, from the time series of "
                      "both programs. Written by md/production/collect_results.py.",
              "written": datetime.datetime.now(IST).strftime("%-d %b %Y, %-I:%M %p IST"),
              "runs": [events_of(lab, cmap[lab]) for lab in ALL]}
    with open(os.path.join(OUT, "run_events.json"), "w") as fh:
        json.dump(events, fh, indent=1)
        fh.write("\n")
    print("written: SOURCES.json, run_health.json, run_events.json in %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
