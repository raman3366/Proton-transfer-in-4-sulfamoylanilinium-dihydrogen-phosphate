#!/usr/bin/env python
"""Physical review: after a change of orientation, is a run from start R distinguishable from the crystal-start runs?
Window: production 20 to 100 ns of every run, cut into eight blocks of 10 ns. For each quantity the 24 block means
of X1, X2, X3 give a mean and a standard deviation; each other run is given with its own mean over its eight blocks
and that mean's distance from the X mean in units of the X block standard deviation (z).
Writes compare.txt and compare.json beside itself."""
import os, json
import numpy as np
from common import *
L = []; J = {}
def P(*a): L.append(" ".join(str(x) for x in a))
runs = [Run(m, r) for m, r in RUNLIST if os.path.exists(f"{HERE}/data_{m}_{r}.npz")]
sn = runs[0].meta["site_names"]; wi = [sn.index(f"W{w}") for w in (423, 516, 553, 568)]
def sol(r, fn):
    return lambda lo, hi: fn(r, (r.t > lo) & (r.t <= hi))
Q = {
 "ligand RMSD to crystal pose (A)": ("s", lambda r, m: r.s["lig_rmsd_X"][m].mean()),
 "d1 N1...Phe131 CZ (A)": ("s", lambda r, m: r.s["d1"][m].mean()),
 "d2 N1...Asn67 OD1 (A)": ("s", lambda r, m: r.s["d2"][m].mean()),
 "frames in state C (%)": ("s", lambda r, m: 100 * (r.s["state"][m] == 0).mean()),
 "V at N2": ("s", lambda r, m: r.s["n2_vol"][m].mean()),
 "sum of angles at N2 (deg)": ("s", lambda r, m: (r.s["a_Zn_N2_S"][m] + r.s["a_Zn_N2_H"][m] + r.s["a_S_N2_H"][m]).mean()),
 "Zn-N2-S (deg)": ("s", lambda r, m: r.s["a_Zn_N2_S"][m].mean()),
 "Zn-N2-S-C4 (deg, circular)": ("s", lambda r, m: circ_mean(r.s["t_Zn_N2_S_C4"][m])),
 "tau_fold (deg)": ("s", lambda r, m: np.maximum(np.abs(r.s["t_C3_C4_S_N2"][m]), np.abs(r.s["t_C5_C4_S_N2"][m])).mean()),
 "C4-S-N2 (deg)": ("s", lambda r, m: r.s["ang_C4_S_N2"][m].mean()),
 "ring axis against the crystal axis (deg)": ("s", lambda r, m: r.s["axis_to_X"][m].mean()),
 "ring normal against the crystal normal (deg)": ("s", lambda r, m: r.s["normal_to_X"][m].mean()),
 "N2-H...OG1 Thr199 (%)": ("s", lambda r, m: 100 * hb(r.s["d_N2_T199OG1"][m], r.s["a_N2H_T199OG1"][m]).mean()),
 "N2...OG1 (A)": ("s", lambda r, m: r.s["d_N2_T199OG1"][m].mean()),
 "Thr199 N-H...O sulfonyl (%)": ("s", lambda r, m: 100 * (hb(r.s["d_T199N_O1"][m], r.s["a_T199NH_O1"][m]) | hb(r.s["d_T199N_O2"][m], r.s["a_T199NH_O2"][m])).mean()),
 "Thr199 N...nearer sulfonyl O (A)": ("s", lambda r, m: np.minimum(r.s["d_T199N_O1"][m], r.s["d_T199N_O2"][m]).mean()),
 "Zn...nearer sulfonyl O (A)": ("s", lambda r, m: np.minimum(r.s["d_Zn_O1"][m], r.s["d_Zn_O2"][m]).mean()),
 "Zn-N2 (A)": ("s", lambda r, m: r.s["d_Zn_N2"][m].mean()),
 "Thr199 OG1-H...Glu106 OE (%)": ("s", lambda r, m: 100 * (hb(r.s["d_T199OG1_E106OE1"][m], r.s["a_T199OG1H_E106OE1"][m]) | hb(r.s["d_T199OG1_E106OE2"][m], r.s["a_T199OG1H_E106OE2"][m])).mean()),
 "CA RMSD (A)": ("s", lambda r, m: r.s["ca_rmsd"][m].mean()),
 "His64 chi1 in g+ (%)": ("s", lambda r, m: 100 * ((r.s["h64_chi1"][m] > 0) & (r.s["h64_chi1"][m] <= 120)).mean()),
 "ring to Phe131 (A)": ("s", lambda r, m: r.s["contacts"][m][:, r.meta["contact_res"].index(131)].mean()),
 "ring to Leu198 (A)": ("s", lambda r, m: r.s["contacts"][m][:, r.meta["contact_res"].index(198)].mean()),
 "ring to Val121 (A)": ("s", lambda r, m: r.s["contacts"][m][:, r.meta["contact_res"].index(121)].mean()),
 "ring to Thr200 (A)": ("s", lambda r, m: r.s["contacts"][m][:, r.meta["contact_res"].index(200)].mean()),
 "ring to Gln92 (A)": ("s", lambda r, m: r.s["contacts"][m][:, r.meta["contact_res"].index(92)].mean()),
 "W423 occupied (%)": ("f", lambda r, m: 100 * (r.f["site_dmin"][m][:, wi[0]] <= 1.4).mean()),
 "W516 occupied (%)": ("f", lambda r, m: 100 * (r.f["site_dmin"][m][:, wi[1]] <= 1.4).mean()),
 "W553 occupied (%)": ("f", lambda r, m: 100 * (r.f["site_dmin"][m][:, wi[2]] <= 1.4).mean()),
 "W568 occupied (%)": ("f", lambda r, m: 100 * (r.f["site_dmin"][m][:, wi[3]] <= 1.4).mean()),
 "water O in the place of the rotated ring (n within 2.5 A)": ("f", lambda r, m: r.f["n_w_Rring_2p5"][m].mean()),
 "water O in the place of the crystal ring (n within 2.5 A)": ("f", lambda r, m: r.f["n_w_Xring_2p5"][m].mean()),
 "water O within 3.5 A of N1 (n)": ("f", lambda r, m: r.f["n_w_N1_3p5"][m].mean()),
 "water O within 3.5 A of O1 (n)": ("f", lambda r, m: r.f["n_w_O1_3p5"][m].mean()),
 "water O within 3.5 A of O2 (n)": ("f", lambda r, m: r.f["n_w_O2_3p5"][m].mean()),
 "water O within 5 A of Zn (n)": ("f", lambda r, m: r.f["n_Ow_Zn_5p0"][m].mean()),
 "water O within 8 A of Zn (n)": ("f", lambda r, m: r.f["n_Ow_Zn_8p0"][m].mean()),
}
blocks = [(20000 + 10000 * i, 30000 + 10000 * i) for i in range(8)]
def block_means(r, kind, fn):
    out = []
    for lo, hi in blocks:
        t = r.t if kind == "s" else r.tf
        m = (t > lo) & (t <= hi)
        out.append(float(fn(r, m)))
    return np.array(out)
names = [r.name for r in runs]
P("window: production 20-100 ns, eight blocks of 10 ns per run. X block mean and sd over the 24 blocks of X1, X2, X3.")
P("quantity | X runs: mean +- sd of blocks (range of the three run means) | then for each other run: mean (z)")
P("runs: " + ", ".join(n for n in names if not n.startswith("2e/X")))
flag = {}
for q, (kind, fn) in Q.items():
    bm = {r.name: block_means(r, kind, fn) for r in runs}
    xs = np.concatenate([bm[n] for n in names if n.startswith("2e/X")])
    xm, xsd = xs.mean(), xs.std(ddof=1)
    xr = [bm[n].mean() for n in names if n.startswith("2e/X")]
    row = []
    J[q] = dict(X_mean=float(xm), X_block_sd=float(xsd), X_run_means=[float(x) for x in xr], runs={})
    for n in names:
        if n.startswith("2e/X"): continue
        m = bm[n].mean(); z = (m - xm) / xsd if xsd > 0 else float("nan")
        J[q]["runs"][n] = dict(mean=float(m), z=float(z), block_sd=float(bm[n].std(ddof=1)))
        row.append("%s %.2f (%+.1f)" % (n, m, z))
        if abs(z) > 2: flag.setdefault(n, []).append((q, float(m), float(z)))
    P(f"{q} | {xm:.2f} +- {xsd:.2f} ({min(xr):.2f} to {max(xr):.2f}) | " + " | ".join(row))
P("")
P("quantities with |z| above 2, per run:")
for n in names:
    if n.startswith("2e/X"): continue
    P(" ", n, ":", "; ".join("%s %.2f (z %+.1f)" % f for f in flag.get(n, [])) or "none")
J["_flags"] = flag
json.dump(J, open(f"{HERE}/compare.json", "w"), indent=1)
open(f"{HERE}/compare.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
