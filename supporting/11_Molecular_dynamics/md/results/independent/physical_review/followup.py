#!/usr/bin/env python
"""Physical review, second pass over the trajectories for questions raised by the first:
 (1) all-atom frames: raw (not imaged) distances Zn-NE2(His94) and N2-Zn, to see whether zinc or ligand are ever
     written in another periodic image than the protein;
 (2) the water at the position of crystal water 553: its distance from N2 and from the nearest protein N/O, and the
     angle O-H...N2 (does it donate to the sulfonamide nitrogen?);
 (3) His64: where its ring nitrogens lie against the positions of crystal waters 516 and 423 (CA fit), and the
     nearest protein N/O or ligand atom to NE2 and to ND1.
One frame per 100 ps (the all-atom frames). Writes followup_<model>_<run>.npz and followup.txt beside itself."""
import os, glob, json, warnings, sys
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
RUNS = "/root/cqc-sfa-md/runs"; SYSTEMS = "/root/cqc-sfa-md/systems"; OUT = os.path.dirname(os.path.abspath(__file__))
REPO = "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate"
LIST = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
SH = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
def mi(D, H, Hinv):
    f = D @ Hinv; f -= np.round(f); D0 = f @ H
    c = D0[:, None, :] + (SH @ H)[None]; k = (c * c).sum(-1).argmin(1)
    return c[np.arange(len(D)), k]
def ang(a, b, c):
    u = a - b; v = c - b
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), float(u @ v))))
def kabsch(X, R):
    gx = X.mean(0); gr = R.mean(0)
    U, S, Vt = np.linalg.svd((X - gx).T @ (R - gr)); d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1.0, 1.0, d]) @ Vt, gx, gr
def dep_of(resid): return resid + 1 if resid <= 124 else resid + 2
wat = {}
for l in open(f"{REPO}/md/system/waters_all.pdb"):
    if l.startswith("HETATM"): wat[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
lines = []
for model, run in LIST:
    rd = f"{RUNS}/{model}/{run}"
    ref = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(f"{SYSTEMS}/{model}/{run[0]}/reference.pdb") if l.startswith("ATOM") and l[12:16].strip() == "CA"])
    u = mda.Universe(f"{rd}/full.pdb", sorted(glob.glob(f"{rd}/full_seg*.dcd")))
    names, resnames, resids = u.atoms.names, u.atoms.resnames, u.atoms.resids
    nsol = int(np.where(resnames == "SN1")[0][-1]) + 1
    def A(dep, name):
        r = dep - 1 if dep <= 125 else dep - 2
        m = np.where((resids[:nsol] == r) & (names[:nsol] == name))[0]; assert len(m) == 1
        return int(m[0])
    def Lg(name):
        m = np.where((resnames == "SN1") & (names == name))[0]; assert len(m) == 1
        return int(m[0])
    ZN = int(np.where(resnames == "ZN1")[0][0]); N2 = Lg("N2"); NE2_94 = A(94, "NE2")
    CA = np.where((names[:nsol] == "CA") & (resnames[:nsol] != "SN1") & (resnames[:nsol] != "ZN1"))[0]
    WO = np.where((resnames == "WAT") & (names == "O"))[0]
    el = np.array([n[0] for n in names[:nsol]])
    pol = np.where(((el == "N") | (el == "O")) & (resnames[:nsol] != "ZN1"))[0]
    pol_lab = np.array([f"{resnames[i]}{dep_of(resids[i]) if resnames[i] not in ('SN1',) else ''}.{names[i]}" for i in pol])
    h64 = {n: A(64, n) for n in ("NE2", "ND1", "HE2")}
    not64 = np.array([resids[i] != resids[h64["NE2"]] for i in pol])
    prot_pol = pol[(resnames[pol] != "SN1")]
    out = {k: [] for k in ("raw_Zn_H94", "raw_N2_Zn", "w553_d", "w553_dN2", "w553_aOHN2", "w553_dprot", "w553_prot", "NE2_W516", "NE2_W423", "ND1_W516", "ND1_W423", "NE2_part_d", "NE2_part", "ND1_part_d", "ND1_part", "NE2_wat_d", "ND1_wat_d", "chi1")}
    from MDAnalysis.lib.distances import calc_dihedrals
    chi = [A(64, n) for n in ("N", "CA", "CB", "CG")]
    for ts in u.trajectory:
        x = ts.positions.astype(float); H = triclinic_vectors(ts.dimensions).astype(float); Hinv = np.linalg.inv(H)
        out["raw_Zn_H94"].append(float(np.linalg.norm(x[ZN] - x[NE2_94]))); out["raw_N2_Zn"].append(float(np.linalg.norm(x[N2] - x[ZN])))
        zn = x[ZN]
        Q, gx, gr = kabsch(x[CA], ref)
        site = lambda r: (r - gr) @ Q.T + gx
        w = zn + mi(x[WO] - zn, H, Hinv)
        sol = x[:nsol].copy()
        # water at W553
        s553 = site(wat[553]); dw = np.sqrt(((w - s553) ** 2).sum(1)); k = int(dw.argmin())
        o = w[k]; oi = WO[k]
        h1 = o + mi((x[oi + 1] - o)[None], H, Hinv)[0]; h2 = o + mi((x[oi + 2] - o)[None], H, Hinv)[0]
        n2 = zn + mi((x[N2] - zn)[None], H, Hinv)[0]
        out["w553_d"].append(float(dw[k])); out["w553_dN2"].append(float(np.linalg.norm(o - n2))); out["w553_aOHN2"].append(max(ang(o, h1, n2), ang(o, h2, n2)))
        dp = np.sqrt(((sol[prot_pol] - o) ** 2).sum(1)); kp = int(dp.argmin())
        out["w553_dprot"].append(float(dp[kp])); out["w553_prot"].append(int(prot_pol[kp]))
        # His64
        for nm in ("NE2", "ND1"):
            p = sol[h64[nm]]
            out[f"{nm}_W516"].append(float(np.linalg.norm(p - site(wat[516])))); out[f"{nm}_W423"].append(float(np.linalg.norm(p - site(wat[423]))))
            d = np.sqrt(((sol[pol] - p) ** 2).sum(1)); d[~not64] = 99.0; kk = int(d.argmin())
            out[f"{nm}_part_d"].append(float(d[kk])); out[f"{nm}_part"].append(int(pol[kk]))
            out[f"{nm}_wat_d"].append(float(np.sqrt(((w - p) ** 2).sum(1)).min()))
        out["chi1"].append(float(np.degrees(calc_dihedrals(*[sol[i] for i in chi]))))
    o_ = {k: np.array(v) for k, v in out.items()}
    np.savez_compressed(f"{OUT}/followup_{model}_{run}.npz", **o_, pol_names=np.array([f"{resnames[i]}{dep_of(resids[i]) if resnames[i] != 'SN1' else ''}.{names[i]}" for i in range(nsol)]))
    n = len(o_["chi1"]); labs = np.array([f"{resnames[i]}{dep_of(resids[i]) if resnames[i] != 'SN1' else ''}.{names[i]}" for i in range(nsol)])
    def top(idx, m, k=3):
        ids, c = np.unique(idx[m], return_counts=True)
        return ", ".join(f"{labs[i]} {100*cc/m.sum():.0f}%" for i, cc in sorted(zip(ids, c), key=lambda z: -z[1])[:k])
    gp = (o_["chi1"] > 0) & (o_["chi1"] <= 120); gm = ~gp
    occ = o_["w553_d"] <= 1.4
    lines.append(f"{model}/{run}: {n} all-atom frames; raw Zn-NE2(His94) max {o_['raw_Zn_H94'].max():.2f} A, raw N2-Zn max {o_['raw_N2_Zn'].max():.2f} A")
    lines.append(f"   water at W553 (within 1.4 A in {100*occ.mean():.1f} % of frames): O...N2 mean {o_['w553_dN2'][occ].mean():.2f} A; O-H...N2 angle >= 135 and O...N2 <= 3.5 in {100*((o_['w553_aOHN2'][occ] >= 135) & (o_['w553_dN2'][occ] <= 3.5)).mean():.1f} % of those frames; nearest protein N/O: {top(o_['w553_prot'], occ)} (mean distance {o_['w553_dprot'][occ].mean():.2f} A)")
    for nm, m in (("g+", gp), ("g-", gm)):
        if m.sum() < 5: continue
        lines.append(f"   His64 {nm} ({m.sum()} frames): NE2 to W516 site {o_['NE2_W516'][m].mean():.2f} A (within 1.5 A in {100*(o_['NE2_W516'][m] <= 1.5).mean():.0f} %), to W423 site {o_['NE2_W423'][m].mean():.2f} (within 1.5 A in {100*(o_['NE2_W423'][m] <= 1.5).mean():.0f} %); ND1 to W516 {o_['ND1_W516'][m].mean():.2f}, to W423 {o_['ND1_W423'][m].mean():.2f}; "
                     f"nearest N/O to NE2: {top(o_['NE2_part'], m)} (mean {o_['NE2_part_d'][m].mean():.2f} A; <= 3.2 A in {100*(o_['NE2_part_d'][m] <= 3.2).mean():.0f} %); to ND1: {top(o_['ND1_part'], m)} (mean {o_['ND1_part_d'][m].mean():.2f} A); nearest water O to NE2 {o_['NE2_wat_d'][m].mean():.2f}, to ND1 {o_['ND1_wat_d'][m].mean():.2f} A")
open(f"{OUT}/followup.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
