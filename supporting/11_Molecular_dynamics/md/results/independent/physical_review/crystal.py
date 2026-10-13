#!/usr/bin/env python
"""Physical review: the same site quantities in the starts as built (start X = heavy atoms of 6RL9 with built hydrogens).
Reads solute.pdb of one run of each start; writes crystal_values.json beside itself."""
import os, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
OUT = os.path.dirname(os.path.abspath(__file__))
def norm(v): return float(np.sqrt((v * v).sum()))
def dist(a, b): return norm(a - b)
def ang(a, b, c):
    u = a - b; v = c - b
    return float(np.degrees(np.arctan2(norm(np.cross(u, v)), float(u @ v))))
def tors(a, b, c, d):
    b1 = c - b; b1 = b1 / norm(b1)
    v = (a - b) - ((a - b) @ b1) * b1; w = (d - c) - ((d - c) @ b1) * b1
    return float(np.degrees(np.arctan2(float(np.cross(b1, v) @ w), float(v @ w))))
def plane(P):
    c = P.mean(0); u, s, vt = np.linalg.svd(P - c); return c, vt[2]
res = {}
for start, path in (("X", "/root/cqc-sfa-md/runs/bonded_2e/X1/solute.pdb"), ("N", "/root/cqc-sfa-md/runs/bonded_2e/N1/solute.pdb"), ("R", "/root/cqc-sfa-md/runs/bonded_2e/R1/solute.pdb")):
    u = mda.Universe(path); x = u.atoms.positions.astype(float)
    names, resids, resnames = u.atoms.names, u.atoms.resids, u.atoms.resnames
    def A(dep, name):
        r = dep - 1 if dep <= 125 else dep - 2
        m = np.where((resids == r) & (names == name))[0]; assert len(m) == 1
        return x[m[0]]
    def L(name):
        m = np.where((resnames == "SN1") & (names == name))[0]; assert len(m) == 1
        return x[m[0]]
    zn = x[np.where(resnames == "ZN1")[0][0]]
    N4 = [L("N2"), A(94, "NE2"), A(96, "NE2"), A(119, "ND1")]
    r = {}
    for nm, p in zip(("N2", "H94", "H96", "H119"), N4): r["d_Zn_" + nm] = dist(zn, p)
    k = 0; lab = ["N2-H94", "N2-H96", "N2-H119", "H94-H96", "H94-H119", "H96-H119"]
    for i in range(4):
        for j in range(i + 1, 4):
            r["angZn_%d_%s" % (k, lab[k])] = ang(N4[i], zn, N4[j]); k += 1
    for h in (94, 96, 119):
        pts = np.array([A(h, n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")])
        c, n = plane(pts); ref = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        if n @ ref < 0: n = -n
        r[f"zn_oop_{h}"] = float((zn - c) @ n)
    r["d_Zn_O1"] = dist(zn, L("O1")); r["d_Zn_O2"] = dist(zn, L("O2"))
    r["d_Zn_E106OE1"] = dist(zn, A(106, "OE1")); r["d_Zn_E106OE2"] = dist(zn, A(106, "OE2"))
    r["d_T199OG1_E106OE1"] = dist(A(199, "OG1"), A(106, "OE1")); r["a_T199OG1H_E106OE1"] = ang(A(199, "OG1"), A(199, "HG1"), A(106, "OE1"))
    r["d_E106OE2_R246N"] = dist(A(106, "OE2"), A(246, "N"))
    r["d_H94ND1_Q92OE1"] = dist(A(94, "ND1"), A(92, "OE1")); r["d_H96ND1_N244O"] = dist(A(96, "ND1"), A(244, "O"))
    r["d_H119NE2_E117OE"] = min(dist(A(119, "NE2"), A(117, "OE1")), dist(A(119, "NE2"), A(117, "OE2")))
    r["h64_chi1"] = tors(A(64, "N"), A(64, "CA"), A(64, "CB"), A(64, "CG")); r["h64_chi2"] = tors(A(64, "CA"), A(64, "CB"), A(64, "CG"), A(64, "ND1"))
    r["d_Zn_H64NE2"] = dist(zn, A(64, "NE2")); r["d_Zn_H64ND1"] = dist(zn, A(64, "ND1"))
    r["omega_S29_P30"] = tors(A(29, "CA"), A(29, "C"), A(30, "N"), A(30, "CA")); r["omega_P201_P202"] = tors(A(201, "CA"), A(201, "C"), A(202, "N"), A(202, "CA"))
    r["d_N2_T199OG1"] = dist(L("N2"), A(199, "OG1")); r["a_N2H_T199OG1"] = ang(L("N2"), L("HN21"), A(199, "OG1"))
    r["d_T199N_O1"] = dist(A(199, "N"), L("O1")); r["d_T199N_O2"] = dist(A(199, "N"), L("O2"))
    r["a_T199NH_O1"] = ang(A(199, "N"), A(199, "H"), L("O1")); r["a_T199NH_O2"] = ang(A(199, "N"), A(199, "H"), L("O2"))
    eZ = (zn - L("N2")) / dist(zn, L("N2")); eH = (L("HN21") - L("N2")) / dist(L("HN21"), L("N2")); eS = (L("S") - L("N2")) / dist(L("S"), L("N2"))
    r["n2_vol"] = float(eZ @ np.cross(eH, eS))
    r["a_Zn_N2_S"] = ang(zn, L("N2"), L("S")); r["a_Zn_N2_H"] = ang(zn, L("N2"), L("HN21")); r["a_S_N2_H"] = ang(L("S"), L("N2"), L("HN21"))
    r["t_Zn_N2_S_C4"] = tors(zn, L("N2"), L("S"), L("C4")); r["t_H_N2_S_C4"] = tors(L("HN21"), L("N2"), L("S"), L("C4"))
    r["t_C3_C4_S_N2"] = tors(L("C3"), L("C4"), L("S"), L("N2")); r["t_C5_C4_S_N2"] = tors(L("C5"), L("C4"), L("S"), L("N2"))
    r["t_H94_Zn_N2_S"] = tors(A(94, "NE2"), zn, L("N2"), L("S")); r["t_H94_Zn_N2_H"] = tors(A(94, "NE2"), zn, L("N2"), L("HN21"))
    r["d_S_N2"] = dist(L("S"), L("N2")); r["d_S_C4"] = dist(L("S"), L("C4")); r["ang_C4_S_N2"] = ang(L("C4"), L("S"), L("N2"))
    r["d1"] = dist(L("N1"), A(131, "CZ")); r["d2"] = dist(L("N1"), A(67, "OD1"))
    r["a_ZnN2_axis"] = float(np.degrees(np.arccos(((L("N2") - zn) / dist(L("N2"), zn)) @ ((L("C1") - L("C4")) / dist(L("C1"), L("C4"))))))
    ca = x[(names == "CA") & (resnames != "SN1")]
    r["rg_ca"] = float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()))
    ph = x[(resnames != "SN1") & (resnames != "ZN1") & (np.array([n[0] for n in names]) != "H")]
    r["rg_heavy"] = float(np.sqrt(((ph - ph.mean(0)) ** 2).sum(1).mean()))
    res[start] = r
json.dump(res, open(f"{OUT}/crystal_values.json", "w"), indent=1)
for k in res["X"]:
    print("%-28s X %9.3f   N %9.3f   R %9.3f" % (k, res["X"][k], res["N"][k], res["R"][k]))
