#!/usr/bin/env python
"""pi_test_look.py: every snapshot a test run left on disk, from the built start to the last frame.

TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper (md/test/PI_TEST.md).

  python pi_test_look.py RUN_DIR --out-prefix md/test/results/X1 [--no-figure]

md/analysis/analyse_md.py reads the production segments only. This script also reads the built start
(system.inpcrd), the minimised state, the equilibration frames and the saved state at the end of every
step, so that the minimisation and the release of the restraints can be looked at. It is written
independently of analyse_md.py (own readers, own formulas after md/analysis/SPEC.md), and its production
rows are compared with timeseries_solute.csv of that program if the file exists.

Snapshots:  built      system.inpcrd                                   all atoms
            min        state_min.xml                                   all atoms
            eqNN       solute_eqNN.dcd frames (solute only), state_eqNN.xml at the end (all atoms)
            segNNNN    solute_segNNNN.dcd frames (solute only), full_segNNNN.dcd frames (all atoms)
Water can be seen in the all-atom snapshots only.

Output:  <prefix>_look.csv   one line per snapshot
         <prefix>_log.csv    the state logs of all steps in one table
         <prefix>_look.json  summary: minimisation, gate quantities, events at the zinc
         <prefix>_look.png   the time series (unless --no-figure)

Units: Angstrom, degrees, ps, K, g/cm3. Time zero is the end of the equilibration (start of production).
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_test_common as C  # noqa: E402

# N-Zn-N angles of 6RL9 (md/analysis/SPEC.md, section 11, third column)
ANGLE_REF = {"ang_N2_Zn_H94NE2_deg": 108.920, "ang_N2_Zn_H96NE2_deg": 114.280, "ang_N2_Zn_H119ND1_deg": 116.833,
             "ang_H94NE2_Zn_H96NE2_deg": 103.819, "ang_H94NE2_Zn_H119ND1_deg": 113.591,
             "ang_H96NE2_Zn_H119ND1_deg": 98.409}
ZN_N_REF = {"d_Zn_N2_A": 1.9314, "d_Zn_H94NE2_A": 2.0099, "d_Zn_H96NE2_A": 2.0206, "d_Zn_H119ND1_A": 2.0183}
GATE = {"d_Zn_N2_A": (1.83, 2.03), "d_Zn_His_A": (1.91, 2.12), "angle_tolerance_deg": 10.0, "fifth_ligand_A": 2.6,
        "fifth_ligand_max_fraction": 0.05, "hbond_d_A": 3.5, "hbond_angle_deg": 135.0, "hbond_min_fraction": 0.8,
        "lig_rmsd_max_A": 1.0, "ang_C4_S_N2_deg": (105.8, 111.8), "tau_fold_deg": (112.1, 152.1),
        "midpoint_d1_A": 6.7, "midpoint_d2_A": 5.2}
NEAR_ZN_LIST_A = 3.0


class System:
    def __init__(self, system_dir, run_dir):
        self.sel = C.read_json(os.path.join(system_dir, "selections.json"))
        s = self.sel
        self.solute = np.array(s["solute"])
        atoms = C.read_pdb_atoms(os.path.join(run_dir, "full.pdb"))
        self.names = [a[1] for a in atoms]
        self.resnames = [a[2] for a in atoms]
        self.n_atoms = len(atoms)
        # residue index (tleap) of every atom, from the residue changes in file order
        res, last, k = [], None, -1
        for a in atoms:
            key = (a[2], a[3])
            if key != last:
                k, last = k + 1, key
            res.append(k)
        self.res_index = np.array(res)
        self.dep = {}
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "..", "system", "residue_map.csv")) as f:
            for row in csv.DictReader(f):
                self.dep[int(row["tleap_index"]) - 1] = int(row["deposited_number"])
        self.crystal_water = dict(zip(s.get("crystal_water_O", []), s.get("crystal_water_numbers_deposited", [])))
        self.water_O = np.array(s["water_O"])
        self.heavy_solute = np.array(sorted(set(s["protein_heavy"]) | set(s["lig_heavy"].values())))
        ref = C.read_pdb_atoms(os.path.join(system_dir, "reference.pdb"))
        self.ref_ca = np.array([a[4] for a in ref if a[1] == "CA"])
        self.ref_lig = {a[1]: np.array(a[4]) for a in ref if a[1] in C.LIG_NAMES and a[2] == "SAN"}
        self.ca = np.array(s["protein_CA"])
        if len(self.ref_ca) != len(self.ca) or len(self.ref_lig) != 11:
            raise SystemExit("reference.pdb does not hold the CA atoms and the 11 ligand atoms expected")

    def label(self, i):
        i = int(i)
        if i in self.crystal_water:
            return f"water {self.crystal_water[i]} O (crystal)"
        if self.resnames[i] in ("WAT", "HOH"):
            return f"water O (added by tleap, atom index {i})"
        if i == self.sel["zn"]:
            return "ZN"
        if i in self.sel["lig_heavy"].values() or self.resnames[i] == "SAN":
            return f"SAN {self.names[i]}"
        r = int(self.res_index[i])
        return f"{self.resnames[i]}{self.dep.get(r, '?')} {self.names[i]}"


def measure(S, pos, cell, all_atoms, built_ca, built_lig):
    """All quantities of one snapshot. pos: (n_atoms, 3) Angstrom, NaN for atoms the snapshot does not hold."""
    s, L = S.sel, S.sel["lig_heavy"]
    zn = pos[s["zn"]]
    P = {k: pos[v] for k, v in L.items()}
    out = {}
    N4 = {"N2": P["N2"], "H94NE2": pos[s["his94_NE2"]], "H96NE2": pos[s["his96_NE2"]], "H119ND1": pos[s["his119_ND1"]]}
    for k, v in N4.items():
        out[f"d_Zn_{k}_A"] = C.dist(zn, v)
    keys = list(N4)
    for a in range(4):
        for b in range(a + 1, 4):
            out[f"ang_{keys[a]}_Zn_{keys[b]}_deg"] = C.angle(N4[keys[a]], zn, N4[keys[b]])
    out["d_Zn_O1_A"], out["d_Zn_O2_A"] = C.dist(zn, P["O1"]), C.dist(zn, P["O2"])
    out["d_Zn_Osulfonyl_min_A"] = min(out["d_Zn_O1_A"], out["d_Zn_O2_A"])
    out["d_Zn_E106OE1_A"], out["d_Zn_E106OE2_A"] = C.dist(zn, pos[s["glu106_OE1"]]), C.dist(zn, pos[s["glu106_OE2"]])
    out["d_Zn_T199OG1_A"] = C.dist(zn, pos[s["thr199_OG1"]])
    # everything heavy of protein and ligand near the zinc (the solute is written whole, no image search)
    d = np.linalg.norm(pos[S.heavy_solute] - zn, axis=1)
    near = [(float(d[k]), int(S.heavy_solute[k])) for k in np.nonzero(d <= NEAR_ZN_LIST_A)[0]]
    out["n_heavy_within_2p6A_of_Zn"] = sum(1 for dd, _i in near if dd <= GATE["fifth_ligand_A"])   # protein and ligand only
    prot_o = [(dd, i) for dd, i in zip(d, S.heavy_solute) if S.names[i].startswith("O") and S.resnames[i] != "SAN"]
    dpo, ipo = min(prot_o)
    out["d_Zn_Oprotein_min_A"], out["protein_O_nearest_Zn"] = float(dpo), S.label(ipo)
    # hydrogen bonds
    og1, hn2 = pos[s["thr199_OG1"]], pos[s["lig_HN2"]]
    out["d_N2_T199OG1_A"] = C.dist(P["N2"], og1)
    out["ang_N2_HN2_T199OG1_deg"] = C.angle(P["N2"], hn2, og1)
    out["hb_N2H_T199OG1"] = int(out["d_N2_T199OG1_A"] <= GATE["hbond_d_A"] and out["ang_N2_HN2_T199OG1_deg"] >= GATE["hbond_angle_deg"])
    hg1 = pos[s["thr199_HG1"]]
    de = [C.dist(og1, pos[s[k]]) for k in ("glu106_OE1", "glu106_OE2")]
    ae = [C.angle(og1, hg1, pos[s[k]]) for k in ("glu106_OE1", "glu106_OE2")]
    out["d_T199OG1_E106OE_min_A"] = min(de)
    out["ang_T199OG1_HG1_E106OE_best_deg"] = max(ae)
    out["hb_T199OG1H_E106OE"] = int(any(x <= 3.5 and y >= 135.0 for x, y in zip(de, ae)))
    tn, th = pos[s["thr199_N"]], pos[s["thr199_H"]]
    dn = [C.dist(tn, P[k]) for k in ("O1", "O2")]
    an = [C.angle(tn, th, P[k]) for k in ("O1", "O2")]
    out["d_T199N_Osulfonyl_min_A"] = min(dn)
    out["hb_T199NH_Osulfonyl"] = int(any(x <= 3.5 and y >= 135.0 for x, y in zip(dn, an)))
    # orientation and the ligand itself
    out["d1_N1_F131CZ_A"] = C.dist(P["N1"], pos[s["phe131_CZ"]])
    out["d2_N1_N67OD1_A"] = C.dist(P["N1"], pos[s["asn67_OD1"]])
    if s.get("asn67_ND2") is not None:
        out["d_N1_N67ND2_A"] = C.dist(P["N1"], pos[s["asn67_ND2"]])
    out["state"] = ("C" if out["d1_N1_F131CZ_A"] < GATE["midpoint_d1_A"] and out["d2_N1_N67OD1_A"] > GATE["midpoint_d2_A"]
                    else "R" if out["d1_N1_F131CZ_A"] > GATE["midpoint_d1_A"] and out["d2_N1_N67OD1_A"] < GATE["midpoint_d2_A"] else "I")
    out["tors_Zn_N2_S_C4_deg"] = C.torsion(zn, P["N2"], P["S"], P["C4"])
    t3, t5 = C.torsion(P["C3"], P["C4"], P["S"], P["N2"]), C.torsion(P["C5"], P["C4"], P["S"], P["N2"])
    out["tau_fold_deg"] = max(abs(t3), abs(t5))
    out["ang_C4_S_N2_deg"] = C.angle(P["C4"], P["S"], P["N2"])
    out["ang_S_N2_Zn_deg"] = C.angle(P["S"], P["N2"], zn)
    out["d_S_N2_A"] = C.dist(P["S"], P["N2"])
    # fits on the CA atoms: to the crystal structure, and to the start as built
    ca = pos[S.ca]
    rot, t = C.kabsch(ca, S.ref_ca)
    out["ca_rmsd_to_crystal_A"] = C.rmsd(ca @ rot.T + t, S.ref_ca)
    v, lab = C.ligand_rmsd_symmetric({k: P[k] @ rot.T + t for k in C.LIG_NAMES}, S.ref_lig)
    out["lig_rmsd_to_crystal_A"], out["lig_rmsd_labelling"] = v, lab
    rot, t = C.kabsch(ca, built_ca)
    out["lig_rmsd_to_built_start_A"] = C.rmsd([P[k] @ rot.T + t for k in C.LIG_NAMES], [built_lig[k] for k in C.LIG_NAMES])
    out["lig_rmsd_to_built_start_no_fit_A"] = C.rmsd([P[k] for k in C.LIG_NAMES], [built_lig[k] for k in C.LIG_NAMES])
    # water, in the all-atom snapshots only
    if all_atoms:
        dw = np.linalg.norm(C.min_image(pos[S.water_O] - zn, cell), axis=1)
        k = int(np.argmin(dw))
        out["d_Zn_Owater_min_A"], out["water_nearest_Zn"] = float(dw[k]), S.label(S.water_O[k])
        out["n_water_O_within_5A_of_Zn"] = int((dw <= 5.0).sum())
        out["n_water_O_within_2p6A_of_Zn"] = int((dw <= GATE["fifth_ligand_A"]).sum())
        near += [(float(dw[j]), int(S.water_O[j])) for j in np.nonzero(dw <= NEAR_ZN_LIST_A)[0]]
        centre = np.mean([P[k2] for k2 in C.LIG_NAMES], axis=0)
        dc = C.min_image(pos[S.water_O] - centre, cell)
        close = dc[np.einsum("ij,ij->i", dc, dc) < 12.0 ** 2] + centre
        ligxyz = np.array([P[k2] for k2 in C.LIG_NAMES])
        dmin = np.sqrt(((close[:, None, :] - ligxyz[None, :, :]) ** 2).sum(axis=2)).min(axis=1) if len(close) else np.zeros(0)
        out["n_water_O_within_3p5A_of_ligand"] = int((dmin <= 3.5).sum())
        out["n_water_O_within_3p5A_of_N1"] = int((np.linalg.norm(close - P["N1"], axis=1) <= 3.5).sum()) if len(close) else 0
        out["d_Zn_O_any_min_A"] = min(out["d_Zn_Owater_min_A"], out["d_Zn_Oprotein_min_A"], out["d_Zn_Osulfonyl_min_A"])
    near.sort()
    four = {L["N2"], s["his94_NE2"], s["his96_NE2"], s["his119_ND1"]}
    out["extra_within_2p6A_of_Zn"] = "; ".join(f"{S.label(i)} {dd:.2f}" for dd, i in near
                                               if dd <= GATE["fifth_ligand_A"] and i not in four)
    out["within_3A_of_Zn"] = "; ".join(f"{S.label(i)} {dd:.2f}" for dd, i in near)
    out["four_N_among_the_nearest_four"] = int({i for _d, i in near[:4]} == four) if len(near) >= 4 else 0
    return out


def collect(run_dir):
    man = C.read_json(os.path.join(run_dir, "manifest.json"))
    system_dir = man["system_dir"]
    S = System(system_dir, run_dir)
    plan = C.expected_steps(man)
    eq_total = sum(s["length_ps"] for s in plan if s["kind"] == "equilibration")
    rows, logs = [], []

    built, built_cell = C.read_inpcrd(os.path.join(system_dir, "system.inpcrd"))
    built_ca = built[S.ca]
    built_lig = {k: built[v] for k, v in S.sel["lig_heavy"].items()}

    def add(phase, source, t_prod_ps, pos, cell, all_atoms):
        m = measure(S, pos, cell, all_atoms, built_ca, built_lig)
        rows.append(dict({"phase": phase, "source": source, "time_ps": None if t_prod_ps is None else round(t_prod_ps, 6),
                          "all_atoms": int(all_atoms)}, **m))

    add("built", "system.inpcrd", None, built, built_cell, True)
    for st in plan:
        tag = st["tag"]
        if not os.path.isfile(os.path.join(run_dir, f"rec_{tag}.json")):
            # a step that was cut off: the whole frames it wrote are read, nothing else of it
            sp = os.path.join(run_dir, f"solute_{tag}.dcd")
            if st["kind"] != "minimisation" and os.path.isfile(sp) and os.path.getsize(sp) > C.DCD_HEADER_BYTES:
                offset = st["time_start_ps"] - (eq_total if st["kind"] == "equilibration" else 0.0)
                _hs, xs, cs = C.read_dcd(sp)
                for k in range(len(xs)):
                    pos = np.full((S.n_atoms, 3), np.nan)
                    pos[S.solute] = xs[k]
                    add(tag, f"solute_{tag}.dcd (step cut off)",
                        offset + (k + 1) * st["every_steps"]["solute_frame_interval_ps"] * st["timestep_fs"] / 1000.0, pos, cs[k], False)
                lp = os.path.join(run_dir, f"log_{tag}.csv")
                if os.path.isfile(lp) and os.path.getsize(lp) > 0:
                    header, log = C.read_log_csv(lp)
                    for r in log:
                        d = dict(zip(header, r))
                        d["phase"] = tag
                        d["time_ps"] = round(d["time_ps"] - (eq_total if st["kind"] == "equilibration" else 0.0), 6)
                        logs.append(d)
            break
        if st["kind"] == "minimisation":
            pos, _v, cell = C.read_state_xml(os.path.join(run_dir, "state_min.xml"))
            add("min", "state_min.xml", None, pos, cell, True)
            continue
        offset = st["time_start_ps"] - (eq_total if st["kind"] == "equilibration" else 0.0)
        dt_ps = st["timestep_fs"] / 1000.0
        hs, xs, cs = C.read_dcd(os.path.join(run_dir, f"solute_{tag}.dcd"))
        hf, xf, cf = C.read_dcd(os.path.join(run_dir, f"full_{tag}.dcd"))
        full_times = {round((j + 1) * st["every_steps"]["full_frame_interval_ps"] * dt_ps, 6): j for j in range(len(xf))}
        for k in range(len(xs)):
            t_in = round((k + 1) * st["every_steps"]["solute_frame_interval_ps"] * dt_ps, 6)
            if t_in in full_times:
                j = full_times[t_in]
                add(tag, f"full_{tag}.dcd", offset + t_in, xf[j].astype(np.float64), cf[j], True)
            else:
                pos = np.full((S.n_atoms, 3), np.nan)
                pos[S.solute] = xs[k]
                add(tag, f"solute_{tag}.dcd", offset + t_in, pos, cs[k], False)
        t_end = offset + st["length_ps"]
        if not (rows and rows[-1]["all_atoms"] and abs((rows[-1]["time_ps"] or 1e9) - t_end) < 1e-6):
            pos, _v, cell = C.read_state_xml(os.path.join(run_dir, f"state_{tag}.xml"))
            if rows and rows[-1]["time_ps"] is not None and abs(rows[-1]["time_ps"] - t_end) < 1e-6:
                rows.pop()                                   # the solute frame of the same instant
            add(tag, f"state_{tag}.xml", t_end, pos, cell, True)
        header, log = C.read_log_csv(os.path.join(run_dir, f"log_{tag}.csv"))
        for r in log:
            d = dict(zip(header, r))
            d["phase"] = tag
            d["time_ps"] = round(d["time_ps"] - (eq_total if st["kind"] == "equilibration" else 0.0), 6)
            logs.append(d)
    return man, S, plan, eq_total, rows, logs


def stats(values):
    v = np.array([x for x in values if x is not None], dtype=np.float64)
    if not len(v):
        return None
    return {"n": int(len(v)), "mean": round(float(v.mean()), 4), "sd": round(float(v.std()), 4),
            "min": round(float(v.min()), 4), "max": round(float(v.max()), 4)}


def summarise(run_dir, man, S, plan, eq_total, rows, logs):
    by = {r["phase"]: r for r in rows if r["phase"] in ("built", "min")}
    prod = [r for r in rows if r["phase"].startswith("seg")]
    eq = [r for r in rows if r["phase"].startswith("eq")]
    out = {"note": C.NOTE, "run_dir": run_dir, "run_id": man["run_id"], "system_dir": man["system_dir"],
           "written": C.now_ist(), "software": C.software(),
           "script_sha256": C.sha256_file(os.path.abspath(__file__)),
           "inputs_sha256": {"manifest.json": C.sha256_file(os.path.join(run_dir, "manifest.json")),
                             "selections.json": C.sha256_file(os.path.join(man["system_dir"], "selections.json")),
                             "reference.pdb": C.sha256_file(os.path.join(man["system_dir"], "reference.pdb")),
                             "system.inpcrd": C.sha256_file(os.path.join(man["system_dir"], "system.inpcrd"))},
           "units": "Angstrom, degrees, ps, K, g/cm3; time zero is the end of the equilibration",
           "equilibration_ps": round(eq_total, 6),
           "production_ps_seen": round(max([r["time_ps"] for r in prod] or [0.0]), 6),
           "production_segments_complete": sum(1 for s in plan if s["kind"] == "production"
                                               and os.path.isfile(os.path.join(run_dir, f"rec_{s['tag']}.json"))),
           "n_snapshots": {"all": len(rows), "equilibration": len(eq), "production": len(prod),
                           "production_all_atoms": sum(r["all_atoms"] for r in prod)}}

    if "min" in by:
        keys = ["d_Zn_N2_A", "d_Zn_H94NE2_A", "d_Zn_H96NE2_A", "d_Zn_H119ND1_A", "d_S_N2_A", "d_Zn_O1_A", "d_Zn_O2_A",
                "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "state", "lig_rmsd_to_crystal_A", "tors_Zn_N2_S_C4_deg",
                "ang_C4_S_N2_deg", "ang_S_N2_Zn_deg", "tau_fold_deg", "d_N2_T199OG1_A", "ang_N2_HN2_T199OG1_deg",
                "hb_N2H_T199OG1", "n_heavy_within_2p6A_of_Zn", "extra_within_2p6A_of_Zn", "d_Zn_Owater_min_A",
                "n_water_O_within_3p5A_of_ligand"] + list(ANGLE_REF)
        rec = C.read_json(os.path.join(run_dir, "rec_min.json"))
        built_xyz, _c = C.read_inpcrd(os.path.join(man["system_dir"], "system.inpcrd"))
        mpos, _v, _c = C.read_state_xml(os.path.join(run_dir, "state_min.xml"))
        moved = {k: round(C.dist(mpos[i], built_xyz[i]), 4) for k, i in S.sel["lig_heavy"].items()}
        moved["ZN"] = round(C.dist(mpos[S.sel["zn"]], built_xyz[S.sel["zn"]]), 4)
        out["minimisation"] = {
            "before_built_start": {k: by["built"].get(k) for k in keys},
            "after_minimisation": {k: by["min"].get(k) for k in keys},
            "ligand_heavy_atom_rmsd_minimised_against_built_start_A": {
                "after_fit_on_CA": round(by["min"]["lig_rmsd_to_built_start_A"], 4),
                "no_fit": round(by["min"]["lig_rmsd_to_built_start_no_fit_A"], 4)},
            "displacement_of_each_atom_no_fit_A": moved,
            "largest_displacement_no_fit_A": max(moved.items(), key=lambda kv: kv[1]),
            "energy_before_kj_per_mol": rec["energies_before"]["potential_energy_kj_per_mol"],
            "energy_after_kj_per_mol": rec["energies_after"]["potential_energy_kj_per_mol"],
            "energies_after": rec["energies_after"],
            "lbfgs_iterations_per_pass": rec["lbfgs_iterations_per_pass"],
            "max_iterations_per_pass": rec["settings"]["max_iterations"],
            "passes_stopped_by_the_iteration_limit": rec["passes_stopped_by_the_iteration_limit"],
            "last_pass_stopped_by_the_iteration_limit": rec["last_pass_stopped_by_the_iteration_limit"],
            "minimiser_last_report": rec.get("minimiser_last_report"),
            "restrained_atoms_displacement_A": rec["restrained_atoms_displacement_A"],
            "wall_s": rec["wall_s"], "threads": rec["platform"]["properties_set"].get("Threads")}

    def block(sub, label):
        if not sub:
            return None
        b = {"what": label, "n_snapshots": len(sub), "n_all_atom_snapshots": sum(r["all_atoms"] for r in sub)}
        for k in list(ZN_N_REF) + list(ANGLE_REF) + ["d_Zn_O1_A", "d_Zn_O2_A", "d_Zn_Osulfonyl_min_A", "d_Zn_E106OE1_A",
                                                     "d_Zn_E106OE2_A", "d_Zn_Oprotein_min_A", "d_N2_T199OG1_A",
                                                     "ang_N2_HN2_T199OG1_deg", "d_T199OG1_E106OE_min_A",
                                                     "lig_rmsd_to_crystal_A", "lig_rmsd_to_built_start_A",
                                                     "ca_rmsd_to_crystal_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A",
                                                     "tors_Zn_N2_S_C4_deg", "ang_C4_S_N2_deg", "tau_fold_deg", "d_S_N2_A",
                                                     "ang_S_N2_Zn_deg", "n_heavy_within_2p6A_of_Zn"]:
            b[k] = stats([r.get(k) for r in sub])
        for k in ("d_Zn_Owater_min_A", "n_water_O_within_5A_of_Zn", "n_water_O_within_2p6A_of_Zn",
                  "n_water_O_within_3p5A_of_ligand", "n_water_O_within_3p5A_of_N1", "d_Zn_O_any_min_A"):
            b[k] = stats([r.get(k) for r in sub if r["all_atoms"]])
        for k in ("hb_N2H_T199OG1", "hb_T199OG1H_E106OE", "hb_T199NH_Osulfonyl", "four_N_among_the_nearest_four"):
            b[k] = {"count": sum(r[k] for r in sub), "n": len(sub)}
        b["fifth_sulfonylO"] = {"count": sum(1 for r in sub if r["d_Zn_Osulfonyl_min_A"] <= GATE["fifth_ligand_A"]), "n": len(sub)}
        b["fifth_proteinO"] = {"count": sum(1 for r in sub if r["d_Zn_Oprotein_min_A"] <= GATE["fifth_ligand_A"]), "n": len(sub),
                               "which": sorted({r["protein_O_nearest_Zn"] for r in sub if r["d_Zn_Oprotein_min_A"] <= GATE["fifth_ligand_A"]})}
        aa = [r for r in sub if r["all_atoms"]]
        b["fifth_water"] = {"count": sum(1 for r in aa if r["d_Zn_Owater_min_A"] <= GATE["fifth_ligand_A"]), "n": len(aa),
                            "which": sorted({r["water_nearest_Zn"] for r in aa if r["d_Zn_Owater_min_A"] <= GATE["fifth_ligand_A"]})}
        b["fifth_any_all_atom_snapshots"] = {"count": sum(1 for r in aa if r["d_Zn_O_any_min_A"] <= GATE["fifth_ligand_A"]), "n": len(aa)}
        b["states"] = {k: sum(1 for r in sub if r["state"] == k) for k in "CRI"}
        return b

    out["equilibration"] = block(eq, "all snapshots of the equilibration (frames and end states)")
    out["production"] = block(prod, "all production frames")

    # gate quantities over the production frames (what little this run can say)
    if prod:
        p = out["production"]
        g = {"frames": len(prod), "all_atom_frames": p["n_all_atom_snapshots"]}
        g["zn_distances"] = {k: {"mean_A": p[k]["mean"], "sd_A": p[k]["sd"], "r0_A": ZN_N_REF[k],
                                 "window_A": list(GATE["d_Zn_N2_A"] if k == "d_Zn_N2_A" else GATE["d_Zn_His_A"]),
                                 "inside": (GATE["d_Zn_N2_A"] if k == "d_Zn_N2_A" else GATE["d_Zn_His_A"])[0] <= p[k]["mean"]
                                 <= (GATE["d_Zn_N2_A"] if k == "d_Zn_N2_A" else GATE["d_Zn_His_A"])[1]} for k in ZN_N_REF}
        g["angles"] = {k: {"mean_deg": p[k]["mean"], "sd_deg": p[k]["sd"], "in_6RL9_deg": ANGLE_REF[k],
                           "deviation_deg": round(p[k]["mean"] - ANGLE_REF[k], 3),
                           "inside": abs(p[k]["mean"] - round(ANGLE_REF[k], 1)) <= GATE["angle_tolerance_deg"]} for k in ANGLE_REF}
        g["fifth_ligand"] = {"sulfonyl_O_within_2p6A_solute_frames": p["fifth_sulfonylO"],
                             "protein_O_within_2p6A_all_frames": p["fifth_proteinO"],
                             "water_O_within_2p6A_all_atom_frames": p["fifth_water"],
                             "any_O_within_2p6A_all_atom_frames": p["fifth_any_all_atom_snapshots"]}
        g["hb_N2H_T199OG1"] = dict(p["hb_N2H_T199OG1"], fraction=round(p["hb_N2H_T199OG1"]["count"] / len(prod), 4))
        g["lig_rmsd_to_crystal_A"] = p["lig_rmsd_to_crystal_A"]
        g["ang_C4_S_N2_deg"], g["tau_fold_deg"] = p["ang_C4_S_N2_deg"], p["tau_fold_deg"]
        out["gate_quantities_production"] = g

    # temperature and density from the state logs: every step, and the last 100 ps of the production
    per = {}
    for tag in dict.fromkeys(d["phase"] for d in logs):
        sub = [d for d in logs if d["phase"] == tag]
        per[tag] = {"rows": len(sub), "temperature_K": stats([d["temperature_K"] for d in sub]),
                    "density_g_per_cm3": stats([d["density_g_per_cm3"] for d in sub]),
                    "volume_nm3": stats([d["volume_nm3"] for d in sub]),
                    "potential_energy_kj_per_mol": stats([d["potential_energy_kj_per_mol"] for d in sub]),
                    "extra_bond_energy_kj_per_mol": stats([d["extra_bond_energy_kj_per_mol"] for d in sub])}
    out["state_log_by_step"] = per
    plog = [d for d in logs if d["phase"].startswith("seg")]
    if plog:
        t_end = max(d["time_ps"] for d in plog)
        last = [d for d in plog if d["time_ps"] > t_end - 100.0 + 1e-9]       # all rows if the production is shorter
        dof = man["system"]["degrees_of_freedom"]
        dens = np.array([d["density_g_per_cm3"] for d in last])
        tt = np.array([d["time_ps"] for d in last])
        out["last_100_ps"] = {"from_ps": round(max(t_end - 100.0, 0.0), 3), "to_ps": round(t_end, 3), "rows": len(last),
                              "note": "the last 100 ps of the production written, or all of it if it is shorter",
                              "temperature_K": stats([d["temperature_K"] for d in last]),
                              "temperature_expected_sd_K_one_sample": round(300.0 * (2.0 / dof) ** 0.5, 3),
                              "density_g_per_cm3": stats(dens), "volume_nm3": stats([d["volume_nm3"] for d in last]),
                              "density_slope_g_per_cm3_per_100_ps": (round(float(np.polyfit(tt, dens, 1)[0]) * 100.0, 5) if len(last) > 2 else None),
                              "potential_energy_kj_per_mol": stats([d["potential_energy_kj_per_mol"] for d in last]),
                              "all_energies_finite": bool(np.isfinite([d["potential_energy_kj_per_mol"] for d in plog]).all())}

    # events at the zinc, in order of time: protein and ligand atoms in every snapshot, water in the all-atom ones
    def changes(sub, keep):
        events, seen = [], ()
        for r in sub:
            items = [e for e in r["extra_within_2p6A_of_Zn"].split("; ") if e and keep(e)]
            names = tuple(sorted(e.rsplit(" ", 1)[0] for e in items))
            if names != seen:
                events.append({"phase": r["phase"], "time_ps": r["time_ps"], "source": r["source"],
                               "within_2p6A_of_Zn": "; ".join(items) or "none"})
                seen = names
        return events
    out["zinc_site_events"] = {
        "what": "each change in the set of atoms within 2.6 A of the zinc beyond the four nitrogens",
        "protein_and_ligand_atoms_every_snapshot": changes(rows, lambda e: not e.startswith("water")),
        "water_all_atom_snapshots_only": changes([r for r in rows if r["all_atoms"]], lambda e: e.startswith("water"))}
    out["end_states"] = [{"phase": r["phase"], "time_ps": r["time_ps"], "within_3A_of_Zn": r["within_3A_of_Zn"],
                          "d_Zn_Owater_min_A": round(r["d_Zn_Owater_min_A"], 3), "water_nearest_Zn": r["water_nearest_Zn"],
                          "n_water_O_within_3p5A_of_ligand": r["n_water_O_within_3p5A_of_ligand"],
                          "n_water_O_within_3p5A_of_N1": r["n_water_O_within_3p5A_of_N1"],
                          "n_water_O_within_5A_of_Zn": r["n_water_O_within_5A_of_Zn"]}
                         for r in rows if r["all_atoms"]]

    # comparison with analyse_md.py, where it has run
    ts = os.path.join(run_dir, "analysis", "timeseries_solute.csv")
    if os.path.isfile(ts) and prod:
        with open(ts) as f:
            theirs = list(csv.DictReader(f))
        pairs = {"lig_rmsd_A": "lig_rmsd_to_crystal_A", "d1_N1_F131CZ_A": "d1_N1_F131CZ_A", "d2_N1_N67OD1_A": "d2_N1_N67OD1_A",
                 "d_N2_T199OG1_A": "d_N2_T199OG1_A", "ang_N2_HN2_T199OG1_deg": "ang_N2_HN2_T199OG1_deg",
                 "hb_N2H_T199OG1": "hb_N2H_T199OG1", "d_Zn_N2_A": "d_Zn_N2_A", "d_Zn_H94NE2_A": "d_Zn_H94NE2_A",
                 "d_Zn_H96NE2_A": "d_Zn_H96NE2_A", "d_Zn_H119ND1_A": "d_Zn_H119ND1_A", "d_Zn_O1_A": "d_Zn_O1_A",
                 "d_Zn_O2_A": "d_Zn_O2_A", "tau_fold_deg": "tau_fold_deg", "ang_C4_S_N2_deg": "ang_C4_S_N2_deg",
                 "d_S_N2_A": "d_S_N2_A", "ca_rmsd_A": "ca_rmsd_to_crystal_A", "ang_S_N2_Zn_deg": "ang_S_N2_Zn_deg",
                 "hb_T199OG1H_E106OE": "hb_T199OG1H_E106OE", "hb_T199NH_Osulfonyl": "hb_T199NH_Osulfonyl"}
        pairs.update({k: k for k in ANGLE_REF})
        cmp = {"file": ts, "sha256": C.sha256_file(ts), "frames_theirs": len(theirs), "frames_here": len(prod),
               "largest_abs_difference": {}}
        if len(theirs) == len(prod):
            for a, b in pairs.items():
                cmp["largest_abs_difference"][a] = round(max(abs(float(t[a]) - float(r[b])) for t, r in zip(theirs, prod)), 5)
            cmp["states_identical"] = all(t["state"] == r["state"] for t, r in zip(theirs, prod))
            dmax = max(v for k, v in cmp["largest_abs_difference"].items() if k.endswith("_A"))
            amax = max(v for k, v in cmp["largest_abs_difference"].items() if k.endswith("_deg"))
            cmp["largest_distance_difference_A"], cmp["largest_angle_difference_deg"] = dmax, amax
            cmp["agree_within_0p001_A_and_0p01_deg"] = bool(dmax <= 0.001 and amax <= 0.01 and cmp["states_identical"]
                                                             and all(cmp["largest_abs_difference"][k] == 0 for k in pairs if k.startswith("hb_")))
        tf = os.path.join(run_dir, "analysis", "timeseries_full.csv")
        if os.path.isfile(tf):
            with open(tf) as f:
                full = list(csv.DictReader(f))
            mine = [r for r in prod if r["source"].startswith("full_")]
            if len(full) == len(mine) and len(full) > 0:
                cmp["full_frames"] = len(full)
                cmp["largest_abs_difference_full"] = {
                    k: round(max(abs(float(t[k]) - float(r[k])) for t, r in zip(full, mine)), 5)
                    for k in ("d_Zn_Owater_min_A", "d_Zn_Oprotein_min_A", "d_Zn_O_any_min_A")}
                cmp["largest_abs_difference_full"]["n_water_N1_3p5A"] = max(
                    abs(int(t["n_water_N1_3p5A"]) - r["n_water_O_within_3p5A_of_N1"]) for t, r in zip(full, mine))
        out["comparison_with_analyse_md"] = cmp
    return out


def figure(path, man, plan, eq_total, rows, logs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    # categorical slots of the project's chart guidance, in their fixed order (validated: results/palette_check.txt)
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
    ink, ink2, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2,
                         "ytick.color": ink2, "text.color": ink, "axes.facecolor": surface, "figure.facecolor": surface,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    gap = max(5.0, 0.04 * (eq_total + sum(s["length_ps"] for s in plan if s["kind"] == "production")))
    x_of = {"built": -eq_total - 2 * gap, "min": -eq_total - gap}
    xs = [x_of.get(r["phase"], r["time_ps"]) for r in rows]
    fig, axes = plt.subplots(5, 2, figsize=(11, 13), sharex=True)

    def series(ax, key, colour, label, only_all=False, marker_only=False, ref=None):
        x = [a for a, r in zip(xs, rows) if (r["all_atoms"] or not only_all) and r.get(key) is not None]
        y = [(r[key] - (ref or 0.0)) for r in rows if (r["all_atoms"] or not only_all) and r.get(key) is not None]
        ax.plot(x, y, color=colour, lw=0 if marker_only else 1.6, marker="o", ms=4 if marker_only else 3, label=label)

    def dress(ax, ylabel, title, lines=()):
        for y, text in lines:
            ax.axhline(y, color=ink2, lw=0.8, ls=(0, (4, 3)))
            ax.annotate(text, xy=(1.0, y), xycoords=("axes fraction", "data"), ha="right", va="bottom", fontsize=7, color=ink2)
        ax.axvline(0.0, color=ink2, lw=0.8)
        ax.axvline(-eq_total, color=ink2, lw=0.8)
        for s in plan:
            if s["kind"] == "equilibration" and s["time_start_ps"] > 0:
                ax.axvline(s["time_start_ps"] - eq_total, color=grid, lw=0.8)
        ax.grid(axis="y", color=grid, lw=0.6)
        ax.set_ylabel(ylabel)
        handles, _labels = ax.get_legend_handles_labels()
        if len(handles) > 1:                                  # the legend sits between the title and the plot
            ax.legend(fontsize=7, ncol=3, loc="lower left", bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2,
                      handlelength=1.6, columnspacing=1.2)
        ax.set_title(title, loc="left", fontsize=9, color=ink, pad=30 if len(handles) > 3 else 18 if len(handles) > 1 else 6)

    ax = axes[0, 0]
    for k, (key, lab) in enumerate((("d_Zn_N2_A", "N2 (ligand)"), ("d_Zn_H94NE2_A", "NE2 His94"),
                                    ("d_Zn_H96NE2_A", "NE2 His96"), ("d_Zn_H119ND1_A", "ND1 His119"))):
        series(ax, key, pal[k], lab)
    dress(ax, "distance (A)", "Zn-N distances (restrained to 1.93, 2.01, 2.02, 2.02 A)")
    ax = axes[0, 1]
    series(ax, "d_Zn_O1_A", pal[0], "sulfonyl O1")
    series(ax, "d_Zn_O2_A", pal[1], "sulfonyl O2")
    series(ax, "d_Zn_Oprotein_min_A", pal[2], "nearest protein O")
    series(ax, "d_Zn_Owater_min_A", pal[3], "nearest water O (all-atom snapshots)", only_all=True, marker_only=True)
    dress(ax, "distance (A)", "Oxygen atoms near the zinc", [(2.6, "fifth ligand: 2.6 A")])
    ax = axes[1, 0]
    labs = {"ang_N2_Zn_H94NE2_deg": "N2-Zn-His94", "ang_N2_Zn_H96NE2_deg": "N2-Zn-His96", "ang_N2_Zn_H119ND1_deg": "N2-Zn-His119",
            "ang_H94NE2_Zn_H96NE2_deg": "His94-Zn-His96", "ang_H94NE2_Zn_H119ND1_deg": "His94-Zn-His119",
            "ang_H96NE2_Zn_H119ND1_deg": "His96-Zn-His119"}
    for k, key in enumerate(ANGLE_REF):
        series(ax, key, pal[k], labs[key], ref=ANGLE_REF[key])
    dress(ax, "angle minus its value in 6RL9 (deg)", "N-Zn-N angles, as deviation from the crystal structure",
          [(10.0, "+10"), (-10.0, "-10")])
    ax = axes[1, 1]
    series(ax, "n_heavy_within_2p6A_of_Zn", pal[0], "protein and ligand atoms (every snapshot)")
    series(ax, "n_water_O_within_2p6A_of_Zn", pal[1], "water O (all-atom snapshots)", only_all=True, marker_only=True)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    dress(ax, "atoms within 2.6 A of Zn", "Atoms within 2.6 A of the zinc (protein and ligand: 4 = the four N)")
    ax = axes[2, 0]
    series(ax, "d_N2_T199OG1_A", pal[0], "N2 ... OG1 Thr199")
    series(ax, "d_T199OG1_E106OE_min_A", pal[1], "OG1 Thr199 ... OE Glu106")
    series(ax, "d_T199N_Osulfonyl_min_A", pal[2], "N Thr199 ... O sulfonyl")
    dress(ax, "distance (A)", "Hydrogen-bond distances at Thr199", [(3.5, "3.5 A")])
    ax = axes[2, 1]
    series(ax, "ang_N2_HN2_T199OG1_deg", pal[0], "N2-H ... OG1 Thr199")
    series(ax, "ang_T199OG1_HG1_E106OE_best_deg", pal[1], "OG1-H ... OE Glu106")
    dress(ax, "angle at the hydrogen (deg)", "Hydrogen-bond angles at Thr199", [(135.0, "135 deg")])
    ax = axes[3, 0]
    series(ax, "lig_rmsd_to_crystal_A", pal[0], "ligand, to the crystal pose")
    series(ax, "lig_rmsd_to_built_start_A", pal[1], "ligand, to the start as built")
    series(ax, "ca_rmsd_to_crystal_A", pal[2], "protein CA, to the crystal")
    dress(ax, "RMSD after the fit on CA (A)", "RMSD", [(1.0, "1.0 A"), (2.0, "2.0 A")])
    ax = axes[3, 1]
    series(ax, "d1_N1_F131CZ_A", pal[0], "d1: N1 to CZ Phe131")
    series(ax, "d2_N1_N67OD1_A", pal[1], "d2: N1 to OD1 Asn67")
    dress(ax, "distance (A)", "Orientation of the ring", [(6.7, "midpoint d1 6.7 A"), (5.2, "midpoint d2 5.2 A")])
    lx = [d["time_ps"] for d in logs]
    ax = axes[4, 0]
    ax.plot(lx, [d["density_g_per_cm3"] for d in logs], color=pal[0], lw=1.6, marker="o", ms=3)
    dress(ax, "density (g/cm3)", "Density, from the state logs")
    ax.set_xlabel("time (ps); 0 = end of the equilibration")
    ax = axes[4, 1]
    ax.plot(lx, [d["temperature_K"] for d in logs], color=pal[0], lw=1.6, marker="o", ms=3)
    dress(ax, "temperature (K)", "Temperature, from the state logs", [(300.0, "300 K")])
    ax.set_xlabel("time (ps); 0 = end of the equilibration")
    fig.suptitle(f"Pipeline test on the Raspberry Pi, run {man['run_id']} of the restrained-zinc model: "
                 f"NOT a result for the paper\nThe two points at the far left are the start as built and the minimised "
                 f"start; vertical lines: start of the equilibration, its stages, start of the production (0).",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(path, dpi=110, metadata={"Software": None})
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    run_dir = os.path.abspath(args.run_dir)
    man, S, plan, eq_total, rows, logs = collect(run_dir)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    cols = []
    for r in rows:
        cols.extend(k for k in r if k not in cols)
    with open(args.out_prefix + "_look.csv", "w", newline="") as f:
        f.write(f"# {C.NOTE}. Units: Angstrom, degrees, ps (time zero = end of the equilibration).\n")
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else ("" if v is None else v)) for k, v in r.items()})
    if logs:
        lcols = ["phase"] + [k for k in logs[0] if k != "phase"]
        with open(args.out_prefix + "_log.csv", "w", newline="") as f:
            f.write(f"# {C.NOTE}. State logs of all steps; time_ps counts from the end of the equilibration.\n")
            w = csv.DictWriter(f, fieldnames=lcols)
            w.writeheader()
            for d in logs:
                w.writerow(d)
    summary = summarise(run_dir, man, S, plan, eq_total, rows, logs)
    C.write_json(args.out_prefix + "_look.json", summary)
    if not args.no_figure:
        figure(args.out_prefix + "_look.png", man, plan, eq_total, rows, logs)
    print(f"{len(rows)} snapshots ({summary['n_snapshots']}); written {args.out_prefix}_look.csv, _log.csv, _look.json"
          + ("" if args.no_figure else ", _look.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
