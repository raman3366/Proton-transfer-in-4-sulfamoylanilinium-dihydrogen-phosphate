#!/usr/bin/env python
"""Physical review: orientation state against the configuration at N2, frame by frame, and what each change of the
orientation state goes with. Writes n2state.txt and n2state.json beside itself."""
import os, json
import numpy as np
from common import *
L = []; J = {}
def P(*a): L.append(" ".join(str(x) for x in a))
ST = "CRI"
runs = [Run(m, r) for m, r in RUNLIST if os.path.exists(f"{HERE}/data_{m}_{r}.npz")]
P("== A. Orientation state (d1, d2) against the sign of the signed volume V at N2; frames from the release of the ligand ==")
P("run | frames | C,V>0 | C,V<0 | R,V>0 | R,V<0 | I,V>0 | I,V<0 | in each cell: number of frames (mean ligand RMSD to the crystal pose / to the rotated start, A)")
for r in runs:
    s = r.s; fr = r.free; V = s["n2_vol"][fr]; st = s["state"][fr]
    row = []; J[r.name] = {"cells": {}}
    for k in (0, 1, 2):
        for sg, nm in ((1, "V>0"), (-1, "V<0")):
            m = (st == k) & (np.sign(V) == sg)
            if m.sum():
                row.append("%d (%.2f / %.2f)" % (m.sum(), s["lig_rmsd_X"][fr][m].mean(), s["lig_rmsd_R"][fr][m].mean()))
                J[r.name]["cells"][f"{ST[k]},{nm}"] = dict(n=int(m.sum()), rmsd_X=float(s["lig_rmsd_X"][fr][m].mean()), rmsd_R=float(s["lig_rmsd_R"][fr][m].mean()),
                    V_mean=float(V[m].mean()), sumang=float((s["a_Zn_N2_S"][fr][m] + s["a_Zn_N2_H"][fr][m] + s["a_S_N2_H"][fr][m]).mean()),
                    t_ZNSC=circ_mean(s["t_Zn_N2_S_C4"][fr][m]), t_ZNSC_sd=circ_sd(s["t_Zn_N2_S_C4"][fr][m]), phiS=circ_mean(s["t_H94_Zn_N2_S"][fr][m]), phiH=circ_mean(s["t_H94_Zn_N2_H"][fr][m]),
                    hb=float(hb(s["d_N2_T199OG1"][fr][m], s["a_N2H_T199OG1"][fr][m]).mean()), hb_T199N=float((hb(s["d_T199N_O1"][fr][m], s["a_T199NH_O1"][fr][m]) | hb(s["d_T199N_O2"][fr][m], s["a_T199NH_O2"][fr][m])).mean()),
                    a_ZnN2S=float(s["a_Zn_N2_S"][fr][m].mean()), a_ZnN2H=float(s["a_Zn_N2_H"][fr][m].mean()), a_SN2H=float(s["a_S_N2_H"][fr][m].mean()),
                    d1=float(s["d1"][fr][m].mean()), d2=float(s["d2"][fr][m].mean()), axis_to_X=float(s["axis_to_X"][fr][m].mean()), axis_to_R=float(s["axis_to_R"][fr][m].mean()),
                    Osulf=float(np.minimum(s["d_Zn_O1"][fr][m], s["d_Zn_O2"][fr][m]).mean()), tfold=float(np.maximum(np.abs(s["t_C3_C4_S_N2"][fr][m]), np.abs(s["t_C5_C4_S_N2"][fr][m])).mean()),
                    angH94=float(s["angZn_0"][fr][m].mean()), angH96=float(s["angZn_1"][fr][m].mean()), angH119=float(s["angZn_2"][fr][m].mean()))
            else: row.append("0")
    P(r.name, "|", int(fr.sum()), "|", " | ".join(row))
P("")
P("== B. Geometry of each occupied cell (frames from release): run cell | n | V mean | sum of the angles at N2 | Zn-N2-S, Zn-N2-H, S-N2-H | Zn-N2-S-C4 circ. mean +- sd | His94NE2-Zn-N2-S (phiS), -H (phiH) | N2-Zn-N(His94, His96, His119) | N2-H...OG1 % | Thr199 N-H...O % | d1, d2 | ring axis to crystal / to rotated start (deg) | nearer sulfonyl O to Zn | tau_fold")
for r in runs:
    for c, j in J[r.name]["cells"].items():
        if j["n"] < 5: continue
        P(f"{r.name} {c} | {j['n']} | {j['V_mean']:+.2f} | {j['sumang']:.1f} | {j['a_ZnN2S']:.1f}, {j['a_ZnN2H']:.1f}, {j['a_SN2H']:.1f} | {j['t_ZNSC']:.1f} +- {j['t_ZNSC_sd']:.1f} | {j['phiS']:.1f}, {j['phiH']:.1f} | {j['angH94']:.1f}, {j['angH96']:.1f}, {j['angH119']:.1f} | {100*j['hb']:.1f} | {100*j['hb_T199N']:.1f} | {j['d1']:.2f}, {j['d2']:.2f} | {j['axis_to_X']:.1f} / {j['axis_to_R']:.1f} | {j['Osulf']:.2f} | {j['tfold']:.1f}")
P("")
P("== C. Distribution of V (production frames): run | mean in frames with V>0 | mean in frames with V<0 | histogram of V in bins of 0.2 from -1 to 1 (% of frames) ==")
edges = np.linspace(-1, 1, 11)
for r in runs:
    V = r.s["n2_vol"][r.prod]
    h, _ = np.histogram(V, edges)
    J[r.name]["V_hist"] = [float(x) for x in 100 * h / len(V)]
    P(r.name, "| %s | %s |" % ("%.3f" % V[V > 0].mean() if (V > 0).any() else "none", "%.3f" % V[V < 0].mean() if (V < 0).any() else "none"), " ".join("%5.2f" % x for x in 100 * h / len(V)))
P("")
P("== D. Sum of the three angles at N2 (360 = planar): production; run | mean +- sd | % of frames above 350 | above 355 | largest ==")
for r in runs:
    s = r.s; p = r.prod; sa = s["a_Zn_N2_S"][p] + s["a_Zn_N2_H"][p] + s["a_S_N2_H"][p]
    J[r.name]["sumang"] = dict(mean=float(sa.mean()), sd=float(sa.std()), above350=float((sa > 350).mean()), above355=float((sa > 355).mean()), max=float(sa.max()))
    P(f"{r.name} | {sa.mean():.1f} +- {sa.std():.1f} | {100*(sa > 350).mean():.2f} | {100*(sa > 355).mean():.2f} | {sa.max():.1f}")
P("")
P("== E. Every change of the orientation state C <-> R (frames of state I belong to the state before them), from release: what the sign of V does within 5 frames (50 ps) before or after, and whether N2-H...OG1 holds in the 11 frames around it ==")
P("run | changes C->R: n; with a change of the sign of V within 50 ps; V<0 on arrival | changes R->C: n; with a sign change within 50 ps; V>0 on arrival | N2-H...OG1 present in all 11 frames around the change: % of changes | dwell times in R (frames of 10 ps): n, median, longest")
for r in runs:
    s = r.s; fr = r.free; V = s["n2_vol"][fr]; st = s["state"][fr]
    start_state = 1 if r.start == "R" else 0
    a = assigned(st, start_state)
    idx = np.where(a[1:] != a[:-1])[0] + 1
    hbf = hb(s["d_N2_T199OG1"][fr], s["a_N2H_T199OG1"][fr])
    sg = np.sign(V)
    cr = [i for i in idx if a[i] == 1]; rc = [i for i in idx if a[i] == 0]
    def sign_change_near(i, w=5):
        lo, hi = max(1, i - w), min(len(sg) - 1, i + w)
        return bool((sg[lo:hi + 1] != sg[lo - 1:hi]).any())
    def hb_all(i, w=5):
        return bool(hbf[max(0, i - w):min(len(hbf), i + w + 1)].all())
    dw = [n for v, i0, n in visits(a) if v == 1]
    J[r.name]["changes"] = dict(n_CR=len(cr), CR_with_signchange=sum(sign_change_near(i) for i in cr), CR_Vneg_on_arrival=sum(1 for i in cr if V[i] < 0),
                                n_RC=len(rc), RC_with_signchange=sum(sign_change_near(i) for i in rc), RC_Vpos_on_arrival=sum(1 for i in rc if V[i] > 0),
                                hb_all_frac=(float(np.mean([hb_all(i) for i in idx])) if len(idx) else None), R_dwells=dw)
    j = J[r.name]["changes"]
    P(f"{r.name} | {j['n_CR']}; {j['CR_with_signchange']}; {j['CR_Vneg_on_arrival']} | {j['n_RC']}; {j['RC_with_signchange']}; {j['RC_Vpos_on_arrival']} | {'none' if j['hb_all_frac'] is None else '%.0f' % (100*j['hb_all_frac'])} | {len(dw)}, {np.median(dw) if dw else 'none'}, {max(dw) if dw else 'none'}")
P("")
P("== F. Dwell of the sign of V (from release): run | stretches with V<0: n, total frames, median, longest (frames) | stretches with V>0: n, median, longest ==")
for r in runs:
    V = r.s["n2_vol"][r.free]; v = visits((V > 0).astype(int))
    neg = [n for a_, i0, n in v if a_ == 0]; pos = [n for a_, i0, n in v if a_ == 1]
    J[r.name]["V_dwell"] = dict(neg_n=len(neg), neg_total=int(sum(neg)), neg_median=float(np.median(neg)) if neg else None, neg_max=int(max(neg)) if neg else 0, pos_n=len(pos), pos_median=float(np.median(pos)) if pos else None, pos_max=int(max(pos)) if pos else 0)
    j = J[r.name]["V_dwell"]
    P(f"{r.name} | {j['neg_n']}, {j['neg_total']}, {j['neg_median']}, {j['neg_max']} | {j['pos_n']}, {j['pos_median']}, {j['pos_max']}")
json.dump(J, open(f"{HERE}/n2state.json", "w"), indent=1)
open(f"{HERE}/n2state.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
