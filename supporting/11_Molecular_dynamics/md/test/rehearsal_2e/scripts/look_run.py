#!/usr/bin/env python
"""look_run.py RUN_DIR OUT.json: the numbers of one REHEARSAL run for md/test/REHEARSAL_2e.md, gathered from the
files of the driver (manifest, records, state logs) and of the primary analysis (analysis/summary.json and its
CSV series); nothing is computed from coordinates here.  A rehearsal: not a result for the paper."""
import csv
import json
import os
import sys

import numpy as np

run, out_path = sys.argv[1], sys.argv[2]
man = json.load(open(os.path.join(run, "manifest.json")))
sm = json.load(open(os.path.join(run, "analysis", "summary.json")))
rk = json.load(open(os.path.join(run, "analysis_recheck", "recheck_report.json"))) if os.path.isfile(
    os.path.join(run, "analysis_recheck", "recheck_report.json")) else None


def table(name):
    with open(os.path.join(run, "analysis", name)) as f:
        rows = list(csv.DictReader(f))
    return rows


def col(rows, key, typ=float):
    return np.array([typ(r[key]) for r in rows])


def stat(x):
    x = np.asarray(x, float)
    return {"n": int(len(x)), "mean": round(float(x.mean()), 4), "sd": round(float(x.std()), 4),
            "min": round(float(x.min()), 4), "max": round(float(x.max()), 4)} if len(x) else None


res = {"what": "REHEARSAL run (not a result for the paper): numbers gathered from the driver's files and from "
               "md/analysis/analyse_md.py's output", "run_dir": run, "run_id": man["run_id"],
       "status": man["status"], "test_only": man["test_only"], "system_dir": man["system_dir"],
       "build_json": {k: man["build_json"].get(k) for k in ("model", "model_kind", "parameter_set", "start", "status_note")},
       "inputs_sha256": {k: (v or {}).get("sha256") for k, v in man["inputs"].items()},
       "script_sha256": man["script"]["sha256"],
       "platform": man["invocations"][-1]["platform"]["name"], "gpu": man["invocations"][-1]["platform"]["gpu_name"],
       "precision": man["invocations"][-1]["platform"]["properties_actual"].get("Precision"),
       "deterministic_forces": man["invocations"][-1]["platform"]["properties_actual"].get("DeterministicForces"),
       "invocations": [{k: i.get(k) for k in ("started", "ended", "exit_code", "outcome")} for i in man["invocations"]]}

# ---- the driver's own records
m = man["minimisation"]
res["minimisation"] = {"energy_before_kj_per_mol": m["energies_before"]["potential_energy_kj_per_mol"],
                       "energy_after_kj_per_mol": m["energies_after"]["potential_energy_kj_per_mol"],
                       "lbfgs_iterations_per_pass": m["lbfgs_iterations_per_pass"],
                       "passes_stopped_by_the_limit": m["passes_stopped_by_the_iteration_limit"],
                       "last_pass_stopped_by_the_limit": m["last_pass_stopped_by_the_iteration_limit"],
                       "wall_s": m["wall_s"], "restrained_atoms_displacement_A": m["restrained_atoms_displacement_A"],
                       "monitor": m["monitor"]}
steps = []
for r in man["equilibration"]["stages"] + man["production"]["segments"]:
    log = list(csv.DictReader(open(os.path.join(run, f"log_{r['tag']}.csv"))))
    T, rho, V = col(log, "temperature_K"), col(log, "density_g_per_cm3"), col(log, "volume_nm3")
    half = len(log) // 2
    steps.append({"tag": r["tag"], "name": r["settings"]["name"], "ensemble": r["settings"]["ensemble"],
                  "timestep_fs": r["settings"]["timestep_fs"], "n_steps": r["settings"]["n_steps"],
                  "length_ps": r["settings"]["length_ps"], "restraints": r["settings"]["restraint_atoms"],
                  "restraint_k": r["settings"]["restraint_k_kj_per_mol_nm2"],
                  "wall_s": r["wall_s"], "ns_per_day": r["ns_per_day"],
                  "T_final_K": r["final"]["temperature_K"], "T_mean_K": round(float(T.mean()), 2), "T_sd_K": round(float(T.std()), 2),
                  "T_min_K": round(float(T.min()), 2), "T_max_K": round(float(T.max()), 2),
                  "T_mean_second_half_K": round(float(T[half:].mean()), 2),
                  "density_final": r["final"]["density_g_per_cm3"], "density_mean_second_half": round(float(rho[half:].mean()), 5),
                  "density_sd_second_half": round(float(rho[half:].std()), 5),
                  "volume_first_row_nm3": round(float(V[0]), 3), "volume_last_row_nm3": round(float(V[-1]), 3),
                  "energy_finite": bool(np.isfinite(col(log, "potential_energy_kj_per_mol")).all()),
                  "frames": r["frames"], "monitor_at_end": r["monitor"], "hydration_check": r.get("hydration_check")})
res["steps"] = steps
seg = list(csv.DictReader(open(os.path.join(run, "log_seg0001.csv"))))
t, V, rho, T = col(seg, "time_ps") / 1000.0, col(seg, "volume_nm3"), col(seg, "density_g_per_cm3"), col(seg, "temperature_K")
slope = np.polyfit(t, V, 1)[0]
res["box_in_production"] = {"volume_mean_nm3": round(float(V.mean()), 3), "volume_sd_nm3": round(float(V.std()), 3),
                            "volume_slope_nm3_per_ns": round(float(slope), 4),
                            "volume_mean_first_ns_nm3": round(float(V[t <= 1.0].mean()), 3),
                            "volume_mean_last_ns_nm3": round(float(V[t > t.max() - 1.0].mean()), 3),
                            "density_mean_g_per_cm3": round(float(rho.mean()), 5), "density_sd": round(float(rho.std()), 5),
                            "density_mean_first_ns": round(float(rho[t <= 1.0].mean()), 5),
                            "density_mean_last_ns": round(float(rho[t > t.max() - 1.0].mean()), 5),
                            "temperature_mean_K": round(float(T.mean()), 2), "temperature_sd_K": round(float(T.std()), 2),
                            "expected_temperature_sd_K": round(300.0 * (2.0 / man["system"]["degrees_of_freedom"]) ** 0.5, 2)}
res["hydration_check"] = {k: man["hydration_check"].get(k) for k in ("status", "applies", "skipped", "counts", "n_extensions",
                                                                    "extension_ps", "final_decision",
                                                                    "site_without_water_when_the_solute_was_released",
                                                                    "settings", "reference_guard")}
res["speed"] = {"equilibration_ns_per_day": man["equilibration"]["ns_per_day"], "production_ns_per_day": man["production"]["ns_per_day"],
                "wall_s_total_of_the_steps": round(m["wall_s"] + man["equilibration"]["wall_s"] + man["production"]["wall_s"], 1)}

# ---- the analysis
g = sm["gate_GV"]
res["analysis"] = {"spec": sm["spec"], "script": sm["script"], "run": {k: sm["run"][k] for k in (
    "label", "model", "start", "replica", "role", "planned_ns", "seen_ns", "complete", "n_frames_solute", "n_frames_full", "sources")},
    "warnings": sm["warnings"],
    "guards": {k: sm["run"].get(k) for k in ("header_check_solute", "header_check_full", "segment_overlap_check",
                                             "manifest_segment_check", "system_files_match_manifest", "imaging_solute",
                                             "largest_neighbour_CA_distance_A", "consistency_full_vs_solute")},
    "gate_GV": g, "orientation": sm["orientation"], "n2_configuration": sm["n2_configuration"],
    "hbond_first_passage": sm["hbond_first_passage"], "start_of_production": sm["start_of_production"],
    "validation": sm["validation"], "lig_rmsd_labelling_counts": sm["lig_rmsd_labelling_counts"],
    "fractions": {k: {x: v[x] for x in ("n", "count", "fraction")} for k, v in sm["fractions"].items()},
    "rmsf": sm["rmsf"]}
keys = ["lig_rmsd_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "d_N1_N67ND2_A", "d_N2_T199OG1_A", "ang_N2_HN2_T199OG1_deg",
        "d_Zn_N2_A", "d_Zn_H94NE2_A", "d_Zn_H96NE2_A", "d_Zn_H119ND1_A", "ang_N2_Zn_H94NE2_deg", "ang_N2_Zn_H96NE2_deg",
        "ang_N2_Zn_H119ND1_deg", "ang_H94NE2_Zn_H96NE2_deg", "ang_H94NE2_Zn_H119ND1_deg", "ang_H96NE2_Zn_H119ND1_deg",
        "ang_S_N2_Zn_deg", "d_Zn_O1_A", "d_Zn_O2_A", "d_Zn_Osulfonyl_min_A", "d_Zn_E106OE1_A", "d_Zn_E106OE2_A",
        "d_Zn_Owater_min_A", "d_Zn_Oprotein_min_A", "d_Zn_O_any_min_A", "tau_fold_deg", "tors_C3_C4_S_N2_deg",
        "tors_C5_C4_S_N2_deg", "ang_C4_S_N2_deg", "d_S_N2_A", "ca_rmsd_A", "n2_volume", "n_water_N1_3p5A",
        "n_water_ring_site_3p5A", "d_water_ring_site_min_A", "d_T199N_Osulfonyl_min_A", "d_T199OG1_E106OE_min_A"]
res["analysis"]["stats"] = {k: {x: sm["stats"][k].get(x) for x in ("n", "mean", "sd", "min", "max", "mean_last_half")}
                            for k in keys if k in sm["stats"]}

# ---- the series, looked at stage by stage (equilibration) and ns by ns (production)
eq, pr = table("timeseries_equil_solute.csv"), table("timeseries_solute.csv")
eqf, prf = table("timeseries_equil_full.csv"), table("timeseries_full.csv")
want = ["lig_rmsd_A", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "n2_volume", "d_Zn_N2_A", "d_Zn_H94NE2_A", "d_Zn_H96NE2_A",
        "d_Zn_H119ND1_A", "d_Zn_Osulfonyl_min_A", "d_Zn_E106OE1_A", "d_N2_T199OG1_A", "ang_N2_HN2_T199OG1_deg",
        "tau_fold_deg", "tors_C3_C4_S_N2_deg", "ang_C4_S_N2_deg", "d_S_N2_A", "ang_S_N2_Zn_deg", "ca_rmsd_A"]


def block(rows, frows):
    b = {"n_frames": len(rows), "state_counts": {s: sum(1 for r in rows if r["state"] == s) for s in "CRI"},
         "n2_sign_counts": {s: sum(1 for r in rows if r["n2_sign"] == s) for s in ("1", "-1", "0")},
         "hb_N2H_T199OG1_frames": sum(int(r["hb_N2H_T199OG1"]) for r in rows),
         "hb_T199NH_Osulfonyl_frames": sum(int(r["hb_T199NH_Osulfonyl"]) for r in rows),
         "hb_T199OG1H_E106OE_frames": sum(int(r["hb_T199OG1H_E106OE"]) for r in rows),
         "fifth_sulfonylO_frames": sum(int(r["fifth_sulfonylO"]) for r in rows),
         "fifth_E106OE_frames": sum(int(r["fifth_E106OE"]) for r in rows),
         "labelling_counts": {s: sum(1 for r in rows if r["lig_rmsd_labelling"] == s) for s in "0123"}}
    for k in want:
        b[k] = stat(col(rows, k))
    if frows:
        b["full_frames"] = {"n": len(frows), "d_Zn_Owater_min_A": stat(col(frows, "d_Zn_Owater_min_A")),
                            "d_Zn_Oprotein_min_A": stat(col(frows, "d_Zn_Oprotein_min_A")),
                            "fifth_any_frames": sum(int(r["fifth_any"]) for r in frows),
                            "fifth_water_frames": sum(int(r["fifth_water"]) for r in frows),
                            "n_water_ring_site_3p5A": stat(col(frows, "n_water_ring_site_3p5A")),
                            "d_water_ring_site_min_A": stat(col(frows, "d_water_ring_site_min_A")),
                            "n_water_N1_3p5A": stat(col(frows, "n_water_N1_3p5A")),
                            "nearest_water_to_zinc_indices": sorted({r["water_O_index_nearest_Zn"] for r in frows})[:12],
                            "nearest_protein_O_to_zinc_indices": sorted({r["protein_O_index_nearest_Zn"] for r in frows})[:12]}
    return b


res["series"] = {"equilibration_by_step": {}, "production_by_ns": {}}
for st in sorted({(int(r["stage"]), int(r["extension"])) for r in eq}):
    rows = [r for r in eq if (int(r["stage"]), int(r["extension"])) == st]
    frows = [r for r in eqf if (int(r["stage"]), int(r["extension"])) == st]
    res["series"]["equilibration_by_step"]["eq%02d" % st[0] + ("x%d" % st[1] if st[1] else "")] = block(rows, frows)
tn = col(pr, "time_ns")
tf = col(prf, "time_ns")
for k in range(int(np.ceil(tn.max() - 1e-9))):
    rows = [r for r, x in zip(pr, tn) if k < x <= k + 1 + 1e-9]
    frows = [r for r, x in zip(prf, tf) if k < x <= k + 1 + 1e-9]
    res["series"]["production_by_ns"][f"{k}-{k + 1} ns"] = block(rows, frows)
res["series"]["production_all"] = block(pr, prf)
# when does anything change: first frames in another state or with another sign at N2, equilibration and production
ev = []
for name, rows in (("equilibration", eq), ("production", pr)):
    last_state, last_sign = None, None
    for r in rows:
        if r["state"] != last_state or r["n2_sign"] != last_sign:
            ev.append({"part": name, "time_ns": float(r["time_ns"]), "stage": r.get("stage"), "state": r["state"],
                       "n2_sign": r["n2_sign"], "d1": float(r["d1_N1_F131CZ_A"]), "d2": float(r["d2_N1_N67OD1_A"]),
                       "n2_volume": float(r["n2_volume"]), "lig_rmsd_A": float(r["lig_rmsd_A"])})
        last_state, last_sign = r["state"], r["n2_sign"]
res["series"]["changes_of_state_or_of_sign_at_N2"] = {"n": len(ev), "first_60": ev[:60]}
res["recheck"] = None if rk is None else {k: rk.get(k) for k in rk if k in ("result", "runs", "n_runs", "agree", "summary", "programs")}
json.dump(res, open(out_path, "w"), indent=1)
print("written", out_path)
