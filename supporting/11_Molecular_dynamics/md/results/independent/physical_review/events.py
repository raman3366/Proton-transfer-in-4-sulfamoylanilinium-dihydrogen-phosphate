#!/usr/bin/env python
"""Physical review: how a change of orientation happens. For every run: the changes of sign of the signed volume at
N2 and the visits of the orientation state from the release of the ligand; for every change, the frames around it
(10 ps) and the water in the two ring places (100 ps). Writes events.txt and events.json beside itself."""
import os, json
import numpy as np
from common import *
L = []; J = {}
def P(*a): L.append(" ".join(str(x) for x in a))
ST = "CRI"
runs = [Run(m, r) for m, r in RUNLIST if os.path.exists(f"{HERE}/data_{m}_{r}.npz")]
sn = runs[0].meta["site_names"]; wi = [sn.index(f"W{w}") for w in (423, 516, 553, 568)]
def frame_line(r, i):
    s = r.s
    sang = s["a_Zn_N2_S"][i] + s["a_Zn_N2_H"][i] + s["a_S_N2_H"][i]
    return ("%7.0f %-7s %s d1 %5.2f d2 %5.2f | V %+.2f sum %5.1f ZnN2S %5.1f ZnN2H %5.1f SN2H %5.1f | phiS %7.1f phiH %7.1f ZnN2SC4 %7.1f C3C4SN2 %7.1f | N2..OG1 %.2f ang %5.1f hb %d | T199N..O1 %.2f O2 %.2f | ZnO1 %.2f ZnO2 %.2f ZnN2 %.2f | rmsdX %.2f rmsdR %.2f axis %4.1f" % (
        r.t[i] - r.trel, s["stage"][i], ST[s["state"][i]], s["d1"][i], s["d2"][i], s["n2_vol"][i], sang, s["a_Zn_N2_S"][i], s["a_Zn_N2_H"][i], s["a_S_N2_H"][i],
        s["t_H94_Zn_N2_S"][i], s["t_H94_Zn_N2_H"][i], s["t_Zn_N2_S_C4"][i], s["t_C3_C4_S_N2"][i], s["d_N2_T199OG1"][i], s["a_N2H_T199OG1"][i],
        int(hb(s["d_N2_T199OG1"][i], s["a_N2H_T199OG1"][i])), s["d_T199N_O1"][i], s["d_T199N_O2"][i], s["d_Zn_O1"][i], s["d_Zn_O2"][i], s["d_Zn_N2"][i], s["lig_rmsd_X"][i], s["lig_rmsd_R"][i], s["axis_to_X"][i]))
for r in runs:
    s = r.s; fr = np.where(r.free)[0]; t = r.t
    start_state = 1 if r.start == "R" else 0; opp = 1 - start_state
    V = s["n2_vol"][fr]; sg = np.sign(V)
    ch = np.where(sg[1:] != sg[:-1])[0] + 1          # index (within fr) of the first frame with the new sign
    a = assigned(s["state"][fr], start_state)
    vis = visits(a)
    long_vis = [(ST[v], float(t[fr][i0] - r.trel), n) for v, i0, n in vis if n >= 100]
    short_opp = [(float(t[fr][i0] - r.trel), n) for v, i0, n in vis if v == opp and n < 100]
    nR = int((s["state"][fr] == 1).sum()); nC = int((s["state"][fr] == 0).sum()); nI = int((s["state"][fr] == 2).sum())
    P("=" * 150)
    P(f"{r.name}  start {r.start}; stages {r.meta['eq_tags']}; release at {r.trel:.0f} ps (production time); frames from release {len(fr)}: C {nC}, R {nR}, I {nI}")
    P(f"  sign changes of V between consecutive frames from release: {len(ch)}; at (ps from release, new sign): " + ", ".join("%.0f %+d" % (t[fr][i] - r.trel, sg[i]) for i in ch[:60]) + (" ..." if len(ch) > 60 else ""))
    P(f"  long visits (>= 100 frames) of the assigned state: {long_vis}")
    P(f"  short visits in the opposite state: {len(short_opp)}; total frames {sum(n for _, n in short_opp)}; the longest: {sorted(short_opp, key=lambda z: -z[1])[:8]}")
    # frames whose orientation state and N2 configuration disagree
    opp_state_frames = s["state"][fr] == opp
    P(f"  frames in the opposite orientation state: {int(opp_state_frames.sum())}; of them with V of the start's sign: {int((opp_state_frames & (np.sign(V) == (1 if start_state == 0 else -1))).sum())}")
    cst = s["state"][fr] == 0; rst = s["state"][fr] == 1
    P(f"  frames in state C with V < 0: {int((cst & (V < 0)).sum())} of {int(cst.sum())}; frames in state R with V > 0: {int((rst & (V > 0)).sum())} of {int(rst.sum())}")
    hbf = hb(s["d_N2_T199OG1"][fr], s["a_N2H_T199OG1"][fr])
    P(f"  N2-H...OG1 from release: {100*hbf.mean():.2f} %; frames without it: {int((~hbf).sum())}; longest stretch without it: {max([n for v, i0, n in visits(hbf.astype(int)) if v == 0] or [0])} frames")
    J[r.name] = dict(sign_changes=[(float(t[fr][i] - r.trel), int(sg[i])) for i in ch], long_visits=long_vis, n_short_opp=len(short_opp), frames_short_opp=int(sum(n for _, n in short_opp)),
                     C_with_Vneg=int((cst & (V < 0)).sum()), R_with_Vpos=int((rst & (V > 0)).sum()), nC=nC, nR=nR, nI=nI, hb_frames_without=int((~hbf).sum()))
    # detail around the first few sign changes and around the first state change
    shown = set()
    centres = [int(i) for i in ch[:6]]
    first_opp = np.where(a == opp)[0]
    if len(first_opp): centres.append(int(first_opp[0]))
    for c in sorted(set(centres)):
        lo, hi = max(0, c - 8), min(len(fr), c + 12)
        if any(abs(c - x) < 6 for x in shown): continue
        shown.add(c)
        P(f"  --- frames around {t[fr][c] - r.trel:.0f} ps from release (time from release, stage, state, ...; phiS = His94NE2-Zn-N2-S, phiH = His94NE2-Zn-N2-H)")
        for k in range(lo, hi): P("   " + frame_line(r, fr[k]))
    # water at 100 ps, from the held part to 3 ns after the release when something happens early; else only a summary
    if r.start == "R" or len(ch) > 0 or len(short_opp) > 0:
        f = r.f
        P("  --- all-atom frames (100 ps): ps from release, stage, state of the solute frame | water O within 2.5 A of the crystal-ring carbons: n, nearest A | of the rotated-ring carbons: n, nearest A | nearest water O to W423, W516, W553, W568 (A)")
        sel = np.where((r.tf <= r.trel + 3000))[0]
        for k in sel:
            d = f["site_dmin"][k][wi]
            P("   %7.0f %-7s %s | %d %5.2f | %d %5.2f | %s" % (r.tf[k] - r.trel, f["stage"][k], ST[r.statef[k]], f["n_w_Xring_2p5"][k], f["d_w_Xring_min"][k], f["n_w_Rring_2p5"][k], f["d_w_Rring_min"][k], " ".join("%5.2f" % x for x in d)))
    # relaxation of the pose after the change: ligand RMSD in windows from the release
    if len(long_vis) and long_vis[0][0] == ST[opp] or (len(long_vis) > 1):
        P("  --- ligand RMSD to the crystal pose, d1, d2, tilt of the ring axis in windows after the release (mean): window ps | RMSD | d1 | d2 | axis to crystal | Zn-N2-S-C4 circ. mean | V mean")
        edges = [0, 100, 200, 300, 400, 500, 600, 800, 1000, 1300, 2000, 3000, 5000, 10000, 20000, 50000, 101300]
        tr = t[fr] - r.trel
        for a0, a1 in zip(edges[:-1], edges[1:]):
            m = (tr > a0) & (tr <= a1)
            if m.sum() == 0: continue
            P("   %6d-%6d | %.2f | %.2f | %.2f | %.1f | %.1f | %+.2f" % (a0, a1, s["lig_rmsd_X"][fr][m].mean(), s["d1"][fr][m].mean(), s["d2"][fr][m].mean(), s["axis_to_X"][fr][m].mean(), circ_mean(s["t_Zn_N2_S_C4"][fr][m]), V[m].mean()))
json.dump(J, open(f"{HERE}/events.json", "w"), indent=1)
open(f"{HERE}/events.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
