#!/usr/bin/env python
"""Physical review, further questions.
 A (from the npz files): water and protein oxygens near the zinc, frame by frame; ligand RMSD above 1 A and water in
   the place of the crystal ring; chloride at the aniline N; blocks of the candidate-terms run; His64 against W516.
 B (one more pass, solute frames every 100 ps): side-chain torsions of Asn62, Asn67, Gln92, Thr200 and His64, to see
   whether the runs from start R keep a memory of the rotated start in the residues the rotated ring touched.
Writes followup2.txt beside itself."""
import os, glob, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from common import *
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
L = []
def P(*a): L.append(" ".join(str(x) for x in a))
runs = [Run(m, r) for m, r in RUNLIST]
sn = runs[0].meta["site_names"]; wi = [sn.index(f"W{w}") for w in (423, 516, 553, 568)]
P("== A1. All-atom frames with a water oxygen within 3.0 A of the zinc (production): run | n | times (ns), distance (A), index of the water oxygen ==")
for r in runs:
    f = r.f; m = np.where(r.prodf & (f["d_Zn_Ow_min"] <= 3.0))[0]
    P(r.name, "|", len(m), "|", "; ".join("%.1f ns %.2f A (atom %d)" % (r.tf[k] / 1000, f["d_Zn_Ow_min"][k], f["i_Zn_Ow_min"][k]) for k in m))
P("")
P("== A2. Ligand RMSD to the crystal pose: % of production frames above 1.0, 1.5, 2.0 A; all-atom frames with a water O within 2.5 A of a crystal ring-carbon position: n, mean RMSD in them, mean RMSD in the other frames; mean tau_fold in them / in the others ==")
for r in runs:
    s = r.s; p = r.prod; rm = s["lig_rmsd_X"][p]
    f = r.f; pf = r.prodf; w = f["n_w_Xring_2p5"][pf] > 0
    rmf = s["lig_rmsd_X"][r.f2s][pf]; tf = np.maximum(np.abs(s["t_C3_C4_S_N2"]), np.abs(s["t_C5_C4_S_N2"]))[r.f2s][pf]
    P(f"{r.name} | {100*(rm > 1.0).mean():.2f}, {100*(rm > 1.5).mean():.2f}, {100*(rm > 2.0).mean():.2f} | {int(w.sum())}, {rmf[w].mean() if w.any() else float('nan'):.2f}, {rmf[~w].mean():.2f} | {tf[w].mean() if w.any() else float('nan'):.0f} / {tf[~w].mean():.0f}")
P("")
P("== A3. Chloride or sodium within 4.0 A of the aniline N1 (all-atom production frames): run | Cl- frames | Na+ frames ==")
for r in runs:
    f = r.f; pf = r.prodf
    P(f"{r.name} | {int((f['d_N1_Cl_min'][pf] <= 4.0).sum())} | {int((f['d_N1_Na_min'][pf] <= 4.0).sum())}")
P("")
P("== A4. Candidate-terms run and the Seminario run in blocks of 10 ns: block | % C / R / I | % V<0 | % |V|<0.2 | ligand RMSD | N2-H...OG1 % | Thr199 N-H...O % | W553 occupied % | water O in the crystal-ring place (mean n) | water within 3.5 A of N2 (mean n) | His64 g+ % ==")
for r in runs:
    if r.name not in ("n2fit/R1", "2s/R1"): continue
    s = r.s; f = r.f
    rows = []
    for b in range(10):
        m = (r.t > b * 10000) & (r.t <= (b + 1) * 10000); mf = (r.tf > b * 10000) & (r.tf <= (b + 1) * 10000)
        st = s["state"][m]; V = s["n2_vol"][m]
        row = [100 * (st == 0).mean(), 100 * (st == 1).mean(), 100 * (st == 2).mean(), 100 * (V < 0).mean(), 100 * (np.abs(V) < 0.2).mean(), s["lig_rmsd_X"][m].mean(),
               100 * hb(s["d_N2_T199OG1"][m], s["a_N2H_T199OG1"][m]).mean(), 100 * (hb(s["d_T199N_O1"][m], s["a_T199NH_O1"][m]) | hb(s["d_T199N_O2"][m], s["a_T199NH_O2"][m])).mean(),
               100 * (f["site_dmin"][mf][:, wi[2]] <= 1.4).mean(), f["n_w_Xring_2p5"][mf].mean(), f["n_w_N2_3p5"][mf].mean(), 100 * ((s["h64_chi1"][m] > 0) & (s["h64_chi1"][m] <= 120)).mean()]
        rows.append(row)
        P(f"{r.name} {b*10}-{(b+1)*10} ns | {row[0]:.1f} / {row[1]:.1f} / {row[2]:.1f} | {row[3]:.1f} | {row[4]:.1f} | {row[5]:.2f} | {row[6]:.1f} | {row[7]:.1f} | {row[8]:.0f} | {row[9]:.2f} | {row[10]:.2f} | {row[11]:.0f}")
    a = np.array(rows)
    P(f"   {r.name}: correlation over the ten blocks of % not-C with: water in the crystal-ring place {np.corrcoef(100 - a[:, 0], a[:, 9])[0, 1]:+.2f}; W553 occupancy {np.corrcoef(100 - a[:, 0], a[:, 8])[0, 1]:+.2f}; % V<0 {np.corrcoef(100 - a[:, 0], a[:, 3])[0, 1]:+.2f}")
P("")
P("== A5. Candidate-terms run, all-atom frames by state of the solute frame: state | n | W553 occupied % | water within 3.5 A of N2 (mean n) | water O in the crystal-ring place: mean n, % of frames with one ==")
for r in runs:
    if r.name != "n2fit/R1": continue
    f = r.f
    for k, nm in ((0, "C"), (1, "R"), (2, "I")):
        m = r.prodf & (r.statef == k)
        P(f"{r.name} {nm} | {int(m.sum())} | {100*(f['site_dmin'][m][:, wi[2]] <= 1.4).mean():.0f} | {f['n_w_N2_3p5'][m].mean():.2f} | {f['n_w_Xring_2p5'][m].mean():.2f}, {100*(f['n_w_Xring_2p5'][m] > 0).mean():.0f}")
    for sg, nm in ((1, "V>0"), (-1, "V<0")):
        m = r.prodf & (np.sign(r.s["n2_vol"][r.f2s]) == sg)
        P(f"{r.name} {nm} | {int(m.sum())} | {100*(f['site_dmin'][m][:, wi[2]] <= 1.4).mean():.0f} | {f['n_w_N2_3p5'][m].mean():.2f} | {f['n_w_Xring_2p5'][m].mean():.2f}, {100*(f['n_w_Xring_2p5'][m] > 0).mean():.0f}")
P("")
P("== A6. Water 516 against His64 (all-atom production frames, all eleven runs pooled and per run): W516 occupied % when His64 chi1 is g+ and Zn...NE2 <= 7.0 A ('deep in') | g+ otherwise | g- ==")
tot = {k: [0, 0] for k in ("deep", "gp", "gm")}
for r in runs:
    s = r.s; f = r.f; pf = r.prodf
    c1 = s["h64_chi1"][r.f2s][pf]; dz = s["d_Zn_H64NE2"][r.f2s][pf]; occ = f["site_dmin"][pf][:, wi[1]] <= 1.4
    gp = (c1 > 0) & (c1 <= 120); deep = gp & (dz <= 7.0); other = gp & ~deep; gm = ~gp
    row = []
    for k, m in (("deep", deep), ("gp", other), ("gm", gm)):
        tot[k][0] += int(occ[m].sum()); tot[k][1] += int(m.sum())
        row.append("%s (%d frames)" % ("%.0f" % (100 * occ[m].mean()) if m.sum() else "none", m.sum()))
    P(r.name, "|", " | ".join(row), "| deep-in, percent of production frames (10 ps series): %.1f; in the last 80 ns: %.1f" % (100 * (((s["h64_chi1"][r.prod] > 0) & (s["h64_chi1"][r.prod] <= 120)) & (s["d_Zn_H64NE2"][r.prod] <= 7.0)).mean(),
      100 * (((s["h64_chi1"][r.t > 20000] > 0) & (s["h64_chi1"][r.t > 20000] <= 120)) & (s["d_Zn_H64NE2"][r.t > 20000] <= 7.0)).mean()))
P("pooled | " + " | ".join("%.0f (%d frames)" % (100 * tot[k][0] / max(1, tot[k][1]), tot[k][1]) for k in ("deep", "gp", "gm")))
P("")
P("== A7. The R state before the conversion (frames from release with state R and V<0) against the held part (stages 2 and 3, ligand heavy atoms held): run | held: V mean, sum of angles, Zn-N2-S, N2-H...OG1 %, Thr199 N-H...O %, ligand RMSD to the rotated start | free: the same, and n ==")
for r in runs:
    if r.start != "R": continue
    s = r.s
    held = (r.t <= r.trel) & (np.array([str(x) for x in s["stage"]]) != "eq01")
    free = r.free & (s["state"] == 1) & (s["n2_vol"] < 0)
    def row(m):
        return "%+.2f, %.1f, %.1f, %.0f, %.0f, %.2f" % (s["n2_vol"][m].mean(), (s["a_Zn_N2_S"][m] + s["a_Zn_N2_H"][m] + s["a_S_N2_H"][m]).mean(), s["a_Zn_N2_S"][m].mean(), 100 * hb(s["d_N2_T199OG1"][m], s["a_N2H_T199OG1"][m]).mean(),
                                               100 * (hb(s["d_T199N_O1"][m], s["a_T199NH_O1"][m]) | hb(s["d_T199N_O2"][m], s["a_T199NH_O2"][m])).mean(), s["lig_rmsd_R"][m].mean())
    P(f"{r.name} | {row(held)} (n {int(held.sum())}) | {row(free) if free.sum() else 'none'} (n {int(free.sum())})")
# ---------------------------------------------------------------- B
P("")
P("== B. Side-chain torsions (solute frames every 100 ps of the production; circular mean +- sd, deg) in two windows: 0-10 ns | 20-100 ns ==")
RUNS_ = "/root/cqc-sfa-md/runs"
def dep_res(dep): return dep - 1 if dep <= 125 else dep - 2
TORS = {"Asn62 chi1": (62, ("N", "CA", "CB", "CG")), "Asn62 chi2": (62, ("CA", "CB", "CG", "OD1")), "Asn67 chi1": (67, ("N", "CA", "CB", "CG")), "Asn67 chi2": (67, ("CA", "CB", "CG", "OD1")),
        "Gln92 chi1": (92, ("N", "CA", "CB", "CG")), "Gln92 chi2": (92, ("CA", "CB", "CG", "CD")), "Gln92 chi3": (92, ("CB", "CG", "CD", "OE1")), "Thr200 chi1": (200, ("N", "CA", "CB", "OG1")),
        "His64 chi1": (64, ("N", "CA", "CB", "CG")), "His64 chi2": (64, ("CA", "CB", "CG", "ND1"))}
res = {}
for model, run in RUNLIST:
    rd = f"{RUNS_}/{model}/{run}"
    u = mda.Universe(f"{rd}/solute.pdb", sorted(glob.glob(f"{rd}/solute_seg*.dcd")))
    idx = {k: [u.select_atoms(f"resid {dep_res(d)} and name {n}")[0].index for n in nm] for k, (d, nm) in TORS.items()}
    vals = {k: [] for k in TORS}
    for ts in u.trajectory[9::10]:
        x = ts.positions
        for k, ii in idx.items(): vals[k].append(float(np.degrees(calc_dihedrals(x[ii[0]], x[ii[1]], x[ii[2]], x[ii[3]]))))
    res[label(model, run)] = {k: np.array(v) for k, v in vals.items()}
for k in TORS:
    P(k + ": " + " ; ".join(f"{n} {circ_mean(v[k][:100]):.0f}+-{circ_sd(v[k][:100]):.0f} | {circ_mean(v[k][200:]):.0f}+-{circ_sd(v[k][200:]):.0f}" for n, v in res.items()))
np.savez_compressed(f"{HERE}/followup2_chi.npz", **{f"{n}|{k}": v[k] for n, v in res.items() for k in v})
open(f"{HERE}/followup2.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
