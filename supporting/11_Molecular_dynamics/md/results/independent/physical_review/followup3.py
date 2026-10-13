#!/usr/bin/env python
"""Physical review: hydration of the site before and after His64 moves in; reads the npz files. Writes followup3.txt."""
import numpy as np
from common import *
L = []
def P(*a): L.append(" ".join(str(x) for x in a))
runs = [Run(m, r) for m, r in RUNLIST]
sn = runs[0].meta["site_names"]; wi = [sn.index(f"W{w}") for w in (423, 516, 553, 568)]
P("== water O within 8 A of the zinc, and in the place of the rotated ring (within 2.5 A of its carbons), in windows; His64 'deep in' (chi1 g+ and Zn...NE2 <= 7.0 A) % in the same windows ==")
P("run | 0-5 ns: n8, nR, deep % | 5-10 ns | 10-20 ns | 20-50 ns | 50-100 ns | first time His64 is deep in for 1 ns without a break (ns)")
for r in runs:
    s = r.s; f = r.f
    deep_s = (s["h64_chi1"] > 0) & (s["h64_chi1"] <= 120) & (s["d_Zn_H64NE2"] <= 7.0)
    row = []
    for a, b in ((0, 5), (5, 10), (10, 20), (20, 50), (50, 100)):
        mf = (r.tf > a * 1000) & (r.tf <= b * 1000); m = (r.t > a * 1000) & (r.t <= b * 1000)
        row.append("%.2f, %.2f, %.0f" % (f["n_Ow_Zn_8p0"][mf].mean(), f["n_w_Rring_2p5"][mf].mean(), 100 * deep_s[m].mean()))
    first = None
    for v, i0, n in visits(deep_s[r.prod].astype(int)):
        if v == 1 and n >= 100: first = r.t[r.prod][i0] / 1000; break
    P(r.name, "|", " | ".join(row), "|", "none" if first is None else "%.2f" % first)
P("")
P("== the same two water counts by His64 class, production all-atom frames pooled over the nine runs of the primary set: class | frames | water O within 8 A of Zn | water O in the place of the rotated ring | W423, W516, W553, W568 occupied % ==")
acc = {"deep in": [], "g+ other": [], "g-": []}
for r in runs:
    if not r.name.startswith("2e/"): continue
    s = r.s; f = r.f; pf = r.prodf
    c1 = s["h64_chi1"][r.f2s][pf]; dz = s["d_Zn_H64NE2"][r.f2s][pf]
    gp = (c1 > 0) & (c1 <= 120); deep = gp & (dz <= 7.0)
    for k, m in (("deep in", deep), ("g+ other", gp & ~deep), ("g-", ~gp)):
        acc[k].append(np.column_stack([f["n_Ow_Zn_8p0"][pf][m], f["n_w_Rring_2p5"][pf][m], (f["site_dmin"][pf][m][:, wi] <= 1.4).astype(float)]))
for k, v in acc.items():
    a = np.vstack(v)
    P(f"{k} | {len(a)} | {a[:, 0].mean():.2f} | {a[:, 1].mean():.2f} | " + ", ".join("%.0f" % (100 * x) for x in a[:, 2:].mean(0)))
P("")
P("== within the 'deep in' class and within the g- class separately: runs from R against runs from X and N (primary set): class | start | frames | water O within 8 A of Zn | in the place of the rotated ring | ligand RMSD ==")
for cls in ("deep in", "g-"):
    for grp, sel in (("R", ("2e/R1", "2e/R2", "2e/R3")), ("X and N", ("2e/X1", "2e/X2", "2e/X3", "2e/N1", "2e/N2", "2e/N3"))):
        a8, aR, rm = [], [], []
        for r in runs:
            if r.name not in sel: continue
            s = r.s; f = r.f; pf = r.prodf & (r.tf > 20000)
            c1 = s["h64_chi1"][r.f2s][pf]; dz = s["d_Zn_H64NE2"][r.f2s][pf]
            gp = (c1 > 0) & (c1 <= 120); m = (gp & (dz <= 7.0)) if cls == "deep in" else ~gp
            a8.extend(f["n_Ow_Zn_8p0"][pf][m]); aR.extend(f["n_w_Rring_2p5"][pf][m]); rm.extend(s["lig_rmsd_X"][r.f2s][pf][m])
        P(f"{cls} | {grp} | {len(a8)} | {np.mean(a8) if a8 else float('nan'):.2f} | {np.mean(aR) if aR else float('nan'):.2f} | {np.mean(rm) if rm else float('nan'):.2f}")
open(f"{HERE}/followup3.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
