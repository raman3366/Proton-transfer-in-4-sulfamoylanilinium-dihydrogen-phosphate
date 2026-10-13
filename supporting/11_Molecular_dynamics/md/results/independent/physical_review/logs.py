#!/usr/bin/env python
"""Physical review: temperature, density, potential energy and volume per segment from log_*.csv and rec_*.json.
Reads only; writes logs_summary.json and logs_table.txt beside itself."""
import os, glob, json, csv, re
import numpy as np
RUNS = "/root/cqc-sfa-md/runs"
OUT = os.path.dirname(os.path.abspath(__file__))
LIST = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
res = {}
lines = []
for model, run in LIST:
    rd = f"{RUNS}/{model}/{run}"
    segs = sorted(glob.glob(f"{rd}/log_seg*.csv"))
    eqs = sorted(glob.glob(f"{rd}/log_eq*.csv"))
    per = []
    allrows = {k: [] for k in ("step", "time_ps", "T", "rho", "U", "K", "V", "speed")}
    prev_last = None
    cont = []
    dts = set()
    for p in segs:
        tag = re.search(r"log_(seg\d+)\.csv", p).group(1)
        rows = list(csv.DictReader(open(p)))
        step = np.array([int(r["step"]) for r in rows]); t = np.array([float(r["time_ps"]) for r in rows])
        T = np.array([float(r["temperature_K"]) for r in rows]); rho = np.array([float(r["density_g_per_cm3"]) for r in rows])
        U = np.array([float(r["potential_energy_kj_per_mol"]) for r in rows]); K = np.array([float(r["kinetic_energy_kj_per_mol"]) for r in rows])
        V = np.array([float(r["volume_nm3"]) for r in rows]); sp = np.array([float(r["speed_ns_per_day"]) for r in rows])
        rec = json.load(open(f"{rd}/rec_{tag}.json"))
        dts.add(rec["settings"]["timestep_fs"])
        # time per step inside the log
        dt_fs = np.diff(t) / np.diff(step) * 1000.0
        per.append(dict(tag=tag, n=len(rows), T=float(T.mean()), T_sd=float(T.std()), rho=float(rho.mean()), rho_sd=float(rho.std()),
                        U=float(U.mean()), U_sd=float(U.std()), V=float(V.mean()), first_step=int(step[0]), last_step=int(step[-1]),
                        first_t=float(t[0]), last_t=float(t[-1]), dt_fs_min=float(dt_fs.min()), dt_fs_max=float(dt_fs.max()),
                        rec_dt_fs=rec["settings"]["timestep_fs"], rec_nsteps=rec["settings"]["n_steps"],
                        rec_first=rec.get("first_step"), rec_last=rec.get("last_step"), wall_s=rec.get("wall_s"),
                        U_first=float(U[0]), U_last=float(U[-1]), rho_first=float(rho[0]), rho_last=float(rho[-1]), T_first=float(T[0]), T_last=float(T[-1])))
        if prev_last is not None:
            cont.append((int(step[0]) - prev_last[0], float(t[0]) - prev_last[1], float(U[0]) - prev_last[2], float(rho[0]) - prev_last[3], float(T[0]) - prev_last[4]))
        prev_last = (int(step[-1]), float(t[-1]), float(U[-1]), float(rho[-1]), float(T[-1]))
        for k, a in zip(allrows, (step, t, T, rho, U, K, V, sp)): allrows[k].extend(a.tolist())
    a = {k: np.array(v) for k, v in allrows.items()}
    segT = np.array([p["T"] for p in per]); segR = np.array([p["rho"] for p in per]); segU = np.array([p["U"] for p in per])
    # within-log row-to-row differences, to judge the differences across segment boundaries
    dU = np.abs(np.diff(a["U"])); dR = np.abs(np.diff(a["rho"])); dT = np.abs(np.diff(a["T"]))
    cont = np.array(cont)
    # linear drift of U and density over the run
    tt = a["time_ps"] / 1000.0
    sU = np.polyfit(tt, a["U"], 1)[0]; sR = np.polyfit(tt, a["rho"], 1)[0]; sT = np.polyfit(tt, a["T"], 1)[0]
    eq = []
    for p in eqs:
        tag = re.search(r"log_(eq\w+)\.csv", p).group(1)
        rows = list(csv.DictReader(open(p)))
        rec = json.load(open(f"{rd}/rec_{tag}.json"))
        T = np.array([float(r["temperature_K"]) for r in rows]); rho = np.array([float(r["density_g_per_cm3"]) for r in rows])
        eq.append(dict(tag=tag, n=len(rows), dt_fs=rec["settings"]["timestep_fs"], length_ps=rec["settings"]["length_ps"], ensemble=rec["settings"]["ensemble"],
                       restraint=rec["settings"]["restraint_atoms"], k=rec["settings"]["restraint_k_kj_per_mol_nm2"], T=float(T.mean()), rho_last=float(rho[-1]), rho=float(rho.mean())))
    res[f"{model}/{run}"] = dict(
        n_segments=len(per), rows=int(len(a["T"])), timestep_fs=sorted(dts),
        T_mean=float(a["T"].mean()), T_sd=float(a["T"].std()), T_seg_min=float(segT.min()), T_seg_max=float(segT.max()),
        rho_mean=float(a["rho"].mean()), rho_sd=float(a["rho"].std()), rho_seg_min=float(segR.min()), rho_seg_max=float(segR.max()),
        U_mean=float(a["U"].mean()), U_sd=float(a["U"].std()), U_seg_min=float(segU.min()), U_seg_max=float(segU.max()),
        V_mean=float(a["V"].mean()), speed_mean=float(a["speed"].mean()),
        T_max=float(a["T"].max()), T_min=float(a["T"].min()), rho_min=float(a["rho"].min()), rho_max=float(a["rho"].max()),
        U_slope_kJ_per_ns=float(sU), rho_slope_per_ns=float(sR), T_slope_K_per_ns=float(sT),
        boundary_step_gaps=sorted(set(int(c) for c in cont[:, 0])), boundary_time_gaps_ps=sorted(set(round(float(c), 3) for c in cont[:, 1])),
        boundary_dU_abs_max=float(np.abs(cont[:, 2]).max()), boundary_dU_abs_mean=float(np.abs(cont[:, 2]).mean()),
        within_dU_abs_mean=float(dU.mean()), within_dU_abs_p99=float(np.percentile(dU, 99)), within_dU_abs_max=float(dU.max()),
        boundary_drho_abs_max=float(np.abs(cont[:, 3]).max()), within_drho_abs_mean=float(dR.mean()), within_drho_abs_p99=float(np.percentile(dR, 99)),
        boundary_dT_abs_max=float(np.abs(cont[:, 4]).max()), within_dT_abs_mean=float(dT.mean()), within_dT_abs_p99=float(np.percentile(dT, 99)),
        log_dt_fs_range=[min(p["dt_fs_min"] for p in per), max(p["dt_fs_max"] for p in per)],
        last_step=per[-1]["last_step"], last_time_ps=per[-1]["last_t"], segments=per, equilibration=eq)
    r = res[f"{model}/{run}"]
    lines.append(f"{model}/{run}: {r['n_segments']} segments, {r['rows']} rows, step {r['timestep_fs']} fs (log: {r['log_dt_fs_range'][0]:.4f}-{r['log_dt_fs_range'][1]:.4f}), last step {r['last_step']}, last time {r['last_time_ps']} ps")
    lines.append(f"   T {r['T_mean']:.2f} +- {r['T_sd']:.2f} K (segment means {r['T_seg_min']:.2f}-{r['T_seg_max']:.2f}; rows {r['T_min']:.1f}-{r['T_max']:.1f}); slope {r['T_slope_K_per_ns']:+.4f} K/ns")
    lines.append(f"   density {r['rho_mean']:.5f} +- {r['rho_sd']:.5f} g/cm3 (segment means {r['rho_seg_min']:.5f}-{r['rho_seg_max']:.5f}; rows {r['rho_min']:.4f}-{r['rho_max']:.4f}); slope {r['rho_slope_per_ns']:+.2e} /ns")
    lines.append(f"   U {r['U_mean']:.0f} +- {r['U_sd']:.0f} kJ/mol (segment means {r['U_seg_min']:.0f} to {r['U_seg_max']:.0f}); slope {r['U_slope_kJ_per_ns']:+.2f} kJ/mol/ns")
    lines.append(f"   boundaries: step gaps {r['boundary_step_gaps']}, time gaps {r['boundary_time_gaps_ps']} ps; |dU| across a boundary max {r['boundary_dU_abs_max']:.0f}, mean {r['boundary_dU_abs_mean']:.0f} (rows inside a segment: mean {r['within_dU_abs_mean']:.0f}, 99th pct {r['within_dU_abs_p99']:.0f}, max {r['within_dU_abs_max']:.0f}); |d rho| max {r['boundary_drho_abs_max']:.5f} (inside: mean {r['within_drho_abs_mean']:.5f}, 99th pct {r['within_drho_abs_p99']:.5f}); |dT| max {r['boundary_dT_abs_max']:.2f} (inside: mean {r['within_dT_abs_mean']:.2f}, 99th pct {r['within_dT_abs_p99']:.2f})")
    lines.append("   equilibration: " + "; ".join(f"{e['tag']} {e['ensemble']} {e['length_ps']:.0f} ps at {e['dt_fs']} fs, {e['restraint']} k={e['k']:.0f}, T {e['T']:.1f}, rho end {e['rho_last']:.4f}" for e in eq))
json.dump(res, open(f"{OUT}/logs_summary.json", "w"), indent=1)
open(f"{OUT}/logs_table.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
