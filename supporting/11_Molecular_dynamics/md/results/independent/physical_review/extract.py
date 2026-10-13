#!/usr/bin/env python
"""Physical review of the MD campaign: per-frame quantities from the raw trajectories.

Own code of the physical reviewer; reads solute_*.dcd, full_*.dcd, solute.pdb, full.pdb and rec_*.json of a
run folder and writes ONE npz file into the scratch folder. Nothing is written into a run folder.
Atoms are chosen by residue number and atom name from the PDB files of the run (not from selections.json);
selections.json is read only to compare a few indices (an independent check of that file).

usage: extract.py <model> <run>      e.g. extract.py bonded_2e R1
"""
import sys, os, glob, json, re, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
from scipy.spatial import cKDTree

RUNS = "/root/cqc-sfa-md/runs"
SYSTEMS = "/root/cqc-sfa-md/systems"
REPO = "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate"
OUT = os.path.dirname(os.path.abspath(__file__))

model, run = sys.argv[1], sys.argv[2]
rd = f"{RUNS}/{model}/{run}"
start = run[0]
sysdir = f"{SYSTEMS}/{model}/{start}"

# ---------------------------------------------------------------- geometry helpers
def norm(v): return float(np.sqrt((v * v).sum()))
def dist(a, b): return norm(a - b)
def ang(a, b, c):
    u = a - b; v = c - b
    return float(np.degrees(np.arctan2(norm(np.cross(u, v)), float(u @ v))))
def tors(a, b, c, d):
    b1 = c - b; b1 = b1 / norm(b1)
    v = (a - b) - ((a - b) @ b1) * b1
    w = (d - c) - ((d - c) @ b1) * b1
    return float(np.degrees(np.arctan2(float(np.cross(b1, v) @ w), float(v @ w))))
def kabsch(X, R):
    """rotation Q and centroids so that (X-gx) @ Q + gr ~ R (row vectors), no reflection; returns Q, gx, gr, rmsd"""
    gx = X.mean(0); gr = R.mean(0)
    A = (X - gx).T @ (R - gr)
    U, S, Vt = np.linalg.svd(A)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    Q = U @ D @ Vt
    diff = (X - gx) @ Q - (R - gr)
    return Q, gx, gr, float(np.sqrt((diff * diff).sum() / len(X)))
def plane(P):
    c = P.mean(0)
    u, s, vt = np.linalg.svd(P - c)
    return c, vt[2], float(np.sqrt((((P - c) @ vt[2]) ** 2).mean()))

SHIFTS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
def min_image(D, H, Hinv):
    """minimum-image version of displacement vectors D (n,3) in the cell H (rows = lattice vectors).
    fractional rounding first, then the 27 neighbouring cells are searched."""
    f = D @ Hinv
    f -= np.round(f)
    D0 = f @ H
    cand = D0[:, None, :] + (SHIFTS @ H)[None, :, :]
    d2 = (cand * cand).sum(-1)
    k = d2.argmin(1)
    return cand[np.arange(len(D)), k], k
def min_image_vec(d, H, Hinv):
    v, k = min_image(d[None, :], H, Hinv)
    return v[0], bool(norm(v[0] - d) > 1e-6)

# ---------------------------------------------------------------- topology by residue number and name
us = mda.Universe(f"{rd}/solute.pdb")
NAT = us.atoms.n_atoms
names = us.atoms.names
resids = us.atoms.resids
resnames = us.atoms.resnames
def dep_of(resid):  # deposited residue number of 6RL9 (126 does not exist)
    return resid + 1 if resid <= 124 else resid + 2
def resid_of(dep):
    return dep - 1 if dep <= 125 else dep - 2
EXPECT = {29: "SER", 30: "PRO", 62: "ASN", 64: "HI", 67: "ASN", 92: "GLN", 94: "HD1", 96: "HD2", 106: "GLU",
          117: "GLU", 119: "HE1", 121: "VAL", 131: "PHE", 143: "VAL", 198: "LEU", 199: "THR", 200: "THR",
          201: "PRO", 202: "PRO", 209: "TRP", 244: "ASN", 246: "ARG"}
for d, rn in EXPECT.items():
    got = set(resnames[resids == resid_of(d)])
    assert len(got) == 1 and list(got)[0].startswith(rn), (d, rn, got)
def A(dep, name):
    m = np.where((resids == resid_of(dep)) & (names == name))[0]
    assert len(m) == 1, (dep, name, m)
    return int(m[0])
def L(name):
    m = np.where((resnames == "SN1") & (names == name))[0]
    assert len(m) == 1, name
    return int(m[0])
ZN = int(np.where(resnames == "ZN1")[0][0])
lig = {n: L(n) for n in ("C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S", "HN11", "HN12", "HN21",
                         "H2", "H3", "H5", "H6")}
HEAVY11 = ["C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S"]
RING = ["C1", "C2", "C3", "C4", "C5", "C6"]
lig_idx_all = np.where(resnames == "SN1")[0]
NE2_94, NE2_96, ND1_119 = A(94, "NE2"), A(96, "NE2"), A(119, "ND1")
IM = {94: [A(94, n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
      96: [A(96, n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
      119: [A(119, n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")]}
T199 = {n: A(199, n) for n in ("N", "H", "OG1", "HG1", "CA", "CB", "O")}
E106 = {n: A(106, n) for n in ("OE1", "OE2", "CD")}
F131CZ = A(131, "CZ"); N67OD1 = A(67, "OD1"); N67ND2 = A(67, "ND2")
H64 = {n: A(64, n) for n in ("N", "CA", "CB", "CG", "ND1", "NE2", "CD2", "CE1")}
Q92OE1, Q92NE2 = A(92, "OE1"), A(92, "NE2")
H94ND1, H94HD1 = A(94, "ND1"), A(94, "HD1")
H96ND1, H96HD1 = A(96, "ND1"), A(96, "HD1")
H119NE2, H119HE2 = A(119, "NE2"), A(119, "HE2")
N244O = A(244, "O")
E117 = [A(117, "OE1"), A(117, "OE2")]
R246N, R246H = A(246, "N"), A(246, "H")
OMEGA = {"S29_P30": [A(29, "CA"), A(29, "C"), A(30, "N"), A(30, "CA")],
         "P201_P202": [A(201, "CA"), A(201, "C"), A(202, "N"), A(202, "CA")]}
CA = np.where((names == "CA") & (resnames != "SN1") & (resnames != "ZN1"))[0]
assert len(CA) == 258
ca_dep = np.array([dep_of(r) for r in resids[CA]])
prot_mask = (resnames != "SN1") & (resnames != "ZN1")
elem = np.array([n[0] for n in names])
prot_heavy = np.where(prot_mask & (elem != "H"))[0]
polar = np.where(prot_mask & ((elem == "N") | (elem == "O")))[0]
polar_label = [f"{resnames[i]}{dep_of(resids[i])}.{names[i]}" if resnames[i] != "ACE" else f"ACE.{names[i]}" for i in polar]
prot_O = np.where(prot_mask & (elem == "O"))[0]
prot_O_label = [f"{resnames[i]}{dep_of(resids[i])}.{names[i]}" for i in prot_O]
# hydrogens bonded to each polar atom (from the start coordinates, same residue, within 1.25 A)
x0 = us.atoms.positions.astype(float)
polar_H = []
for i in polar:
    m = np.where((resids == resids[i]) & (elem == "H"))[0]
    dd = np.sqrt(((x0[m] - x0[i]) ** 2).sum(1))
    polar_H.append(m[dd < 1.25])
CONTACT_RES = [62, 64, 65, 67, 69, 91, 92, 94, 96, 119, 121, 131, 135, 141, 143, 198, 199, 200, 202, 204, 209]
contact_idx = [np.where((resids == resid_of(d)) & (elem != "H"))[0] for d in CONTACT_RES]
ringN1 = [lig[n] for n in RING + ["N1"]]

# reference (6RL9 frame): CA and the crystal ligand from reference.pdb; own reading of the file
ref_ca, ref_lig = [], {}
for line in open(f"{sysdir}/reference.pdb"):
    if line.startswith(("ATOM", "HETATM")):
        nm = line[12:16].strip(); xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if line[17:20] == "SN1": ref_lig[nm] = np.array(xyz)
        elif nm == "CA": ref_ca.append(xyz)
ref_ca = np.array(ref_ca); assert len(ref_ca) == 258
def read_pdb_atoms(path):
    out = {}
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out[(line[12:16].strip(), int(line[22:26]))] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out
ligR = {k[0]: v for k, v in read_pdb_atoms(f"{REPO}/md/system/ligand_R.pdb").items()}
ligX = {k[0]: v for k, v in read_pdb_atoms(f"{REPO}/md/system/ligand_X.pdb").items()}
for n in HEAVY11: assert dist(ligX[n], ref_lig[n]) < 1e-3
wat = read_pdb_atoms(f"{REPO}/md/system/waters_all.pdb")
WSITES = {w: wat[("O", w)] for w in (423, 516, 553, 568)}
zn_x = list(read_pdb_atoms(f"{REPO}/md/system/zinc.pdb").values())[0]
PERM = [dict(), {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3"}, {"O1": "O2", "O2": "O1"},
        {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3", "O1": "O2", "O2": "O1"}]
def lig_rmsd(P, ref):
    best = 1e9
    for pm in PERM:
        s = 0.0
        for n in HEAVY11:
            s += ((P[n] - ref[pm.get(n, n)]) ** 2).sum()
        best = min(best, np.sqrt(s / 11.0))
    return float(best)
ref_axis = ref_lig["C1"] - ref_lig["C4"]; ref_axis /= norm(ref_axis)
refR_axis = ligR["C1"] - ligR["C4"]; refR_axis /= norm(refR_axis)
_, ref_nrm, _ = plane(np.array([ref_lig[n] for n in RING]))
_, refR_nrm, _ = plane(np.array([ligR[n] for n in RING]))

# compare a few indices with selections.json (check of that file, not a source)
sel = json.load(open(f"{sysdir}/selections.json"))
chk = {"zn": ZN, "his94_NE2": NE2_94, "his96_NE2": NE2_96, "his119_ND1": ND1_119, "thr199_OG1": T199["OG1"],
       "thr199_HG1": T199["HG1"], "thr199_N": T199["N"], "thr199_H": T199["H"], "glu106_OE1": E106["OE1"],
       "glu106_OE2": E106["OE2"], "phe131_CZ": F131CZ, "asn67_OD1": N67OD1, "lig_HN2": lig["HN21"]}
sel_mismatch = [k for k, v in chk.items() if sel[k] != v]
sel_mismatch += [n for n in HEAVY11 if sel["lig_heavy"][n] != lig[n]]
if list(sel["protein_CA"]) != list(CA): sel_mismatch.append("protein_CA")
if list(sel["solute"]) != list(range(NAT)): sel_mismatch.append("solute")

# ---------------------------------------------------------------- stage list and time axis
def stage_files(kind):
    eq = sorted(glob.glob(f"{rd}/{kind}_eq*.dcd"))
    seg = sorted(glob.glob(f"{rd}/{kind}_seg*.dcd"))
    return eq, seg
eq_s, seg_s = stage_files("solute")
eq_f, seg_f = stage_files("full")
def tag_of(p): return re.search(r"_(eq\w+|seg\d+)\.dcd$", p).group(1)
eq_len = {}
for p in eq_s:
    t = tag_of(p)
    rec = json.load(open(f"{rd}/rec_{t}.json"))
    eq_len[t] = float(rec["settings"]["length_ps"])
eq_tags = [tag_of(p) for p in eq_s]
tot_eq = sum(eq_len[t] for t in eq_tags)
eq_start = {}
acc = -tot_eq
for t in eq_tags:
    eq_start[t] = acc; acc += eq_len[t]
assert abs(acc) < 1e-6
T_RELEASE = eq_start["eq04"]

# ---------------------------------------------------------------- solute frames
S = {}
def put(k, v): S.setdefault(k, []).append(v)
partner_idx = []; partner_code = []
def solute_frame(x, dims, t_ps, stage):
    x = x.astype(float)
    H = triclinic_vectors(dims).astype(float); Hinv = np.linalg.inv(H)
    put("t_ps", t_ps); put("stage", stage); put("box_a", float(dims[0]))
    # raw distances before any imaging
    put("raw_Zn_H94NE2", dist(x[ZN], x[NE2_94])); put("raw_Zn_N2", dist(x[ZN], x[lig["N2"]]))
    v, sh = min_image_vec(x[ZN] - x[NE2_94], H, Hinv)
    put("zn_shifted", int(sh)); x[ZN] = x[NE2_94] + v
    v, sh = min_image_vec(x[lig["N2"]] - x[ZN], H, Hinv)
    put("lig_shifted", int(sh)); x[lig_idx_all] += (x[ZN] + v) - x[lig["N2"]]
    # protein whole?
    ca = x[CA]
    dca = np.sqrt(((ca[1:] - ca[:-1]) ** 2).sum(1))
    put("ca_step_max", float(dca.max()))
    Q, gx, gr, rms = kabsch(ca, ref_ca)
    put("ca_rmsd", rms)
    put("rg_ca", float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean())))
    ph = x[prot_heavy]
    put("rg_heavy", float(np.sqrt(((ph - ph.mean(0)) ** 2).sum(1).mean())))
    P = {n: x[i] for n, i in lig.items()}
    zn = x[ZN]
    # ligand internal
    put("d_S_N2", dist(P["S"], P["N2"])); put("d_S_C4", dist(P["S"], P["C4"]))
    put("d_S_O1", dist(P["S"], P["O1"])); put("d_S_O2", dist(P["S"], P["O2"]))
    put("d_C1_N1", dist(P["C1"], P["N1"])); put("d_N2_H", dist(P["N2"], P["HN21"]))
    cc = [dist(P[RING[i]], P[RING[(i + 1) % 6]]) for i in range(6)]
    put("ring_cc_min", min(cc)); put("ring_cc_max", max(cc))
    rc, rn, rrms = plane(np.array([P[n] for n in RING]))
    put("ring_rms_plane", rrms)
    put("oop_S", abs(float((P["S"] - rc) @ rn))); put("oop_N1", abs(float((P["N1"] - rc) @ rn)))
    put("ang_C4_S_N2", ang(P["C4"], P["S"], P["N2"])); put("ang_O1_S_O2", ang(P["O1"], P["S"], P["O2"]))
    put("ang_O1_S_C4", ang(P["O1"], P["S"], P["C4"])); put("ang_O2_S_C4", ang(P["O2"], P["S"], P["C4"]))
    # zinc site
    N4 = [P["N2"], x[NE2_94], x[NE2_96], x[ND1_119]]
    for nm, p in zip(("N2", "H94", "H96", "H119"), N4): put("d_Zn_" + nm, dist(zn, p))
    k = 0
    for i in range(4):
        for j in range(i + 1, 4):
            put(f"angZn_{k}", ang(N4[i], zn, N4[j])); k += 1
    for r, idx in IM.items():
        pts = x[idx]
        c, n, _ = plane(pts)
        ref_n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        if n @ ref_n < 0: n = -n
        put(f"zn_oop_{r}", float((zn - c) @ n))
    put("d_Zn_O1", dist(zn, P["O1"])); put("d_Zn_O2", dist(zn, P["O2"]))
    put("d_Zn_E106OE1", dist(zn, x[E106["OE1"]])); put("d_Zn_E106OE2", dist(zn, x[E106["OE2"]]))
    # Thr199 OG1-HG1 ... Glu106 OE
    og, hg = x[T199["OG1"]], x[T199["HG1"]]
    for oe in ("OE1", "OE2"):
        put(f"d_T199OG1_E106{oe}", dist(og, x[E106[oe]])); put(f"a_T199OG1H_E106{oe}", ang(og, hg, x[E106[oe]]))
    put("d_T199HG1_N2", dist(hg, P["N2"])); put("d_T199HG1_O1", dist(hg, P["O1"])); put("d_T199HG1_O2", dist(hg, P["O2"]))
    put("d_E106OE2_R246N", dist(x[E106["OE2"]], x[R246N])); put("d_E106OE1_R246N", dist(x[E106["OE1"]], x[R246N]))
    put("a_R246NH_E106OE2", ang(x[R246N], x[R246H], x[E106["OE2"]])); put("a_R246NH_E106OE1", ang(x[R246N], x[R246H], x[E106["OE1"]]))
    # second shell of the zinc histidines
    put("d_H94ND1_Q92OE1", dist(x[H94ND1], x[Q92OE1])); put("a_H94ND1H_Q92OE1", ang(x[H94ND1], x[H94HD1], x[Q92OE1]))
    put("d_H94ND1_Q92NE2", dist(x[H94ND1], x[Q92NE2]))
    put("d_H96ND1_N244O", dist(x[H96ND1], x[N244O])); put("a_H96ND1H_N244O", ang(x[H96ND1], x[H96HD1], x[N244O]))
    dd = [dist(x[H119NE2], x[e]) for e in E117]
    aa = [ang(x[H119NE2], x[H119HE2], x[e]) for e in E117]
    kk = int(np.argmin(dd)); put("d_H119NE2_E117OE", dd[kk]); put("a_H119NE2H_E117OE", aa[kk])
    # His64
    put("h64_chi1", tors(x[H64["N"]], x[H64["CA"]], x[H64["CB"]], x[H64["CG"]]))
    put("h64_chi2", tors(x[H64["CA"]], x[H64["CB"]], x[H64["CG"]], x[H64["ND1"]]))
    put("d_Zn_H64NE2", dist(zn, x[H64["NE2"]])); put("d_Zn_H64ND1", dist(zn, x[H64["ND1"]]))
    put("d_N1_H64ring", min(dist(P["N1"], x[H64[n]]) for n in ("ND1", "NE2", "CD2", "CE1", "CG")))
    for nm, idx in OMEGA.items(): put("omega_" + nm, tors(*[x[i] for i in idx]))
    # ligand to protein
    put("d_N2_T199OG1", dist(P["N2"], og)); put("a_N2H_T199OG1", ang(P["N2"], P["HN21"], og))
    tn, th = x[T199["N"]], x[T199["H"]]
    put("d_T199N_O1", dist(tn, P["O1"])); put("a_T199NH_O1", ang(tn, th, P["O1"]))
    put("d_T199N_O2", dist(tn, P["O2"])); put("a_T199NH_O2", ang(tn, th, P["O2"]))
    put("d_T199N_N2", dist(tn, P["N2"])); put("a_T199NH_N2", ang(tn, th, P["N2"]))
    # configuration at N2
    eZ = (zn - P["N2"]) / dist(zn, P["N2"]); eH = (P["HN21"] - P["N2"]) / dist(P["HN21"], P["N2"]); eS = (P["S"] - P["N2"]) / dist(P["S"], P["N2"])
    put("n2_vol", float(eZ @ np.cross(eH, eS)))
    put("a_Zn_N2_S", ang(zn, P["N2"], P["S"])); put("a_Zn_N2_H", ang(zn, P["N2"], P["HN21"])); put("a_S_N2_H", ang(P["S"], P["N2"], P["HN21"]))
    put("t_Zn_N2_S_C4", tors(zn, P["N2"], P["S"], P["C4"])); put("t_H_N2_S_C4", tors(P["HN21"], P["N2"], P["S"], P["C4"]))
    put("t_C3_C4_S_N2", tors(P["C3"], P["C4"], P["S"], P["N2"])); put("t_C5_C4_S_N2", tors(P["C5"], P["C4"], P["S"], P["N2"]))
    put("t_H94_Zn_N2_S", tors(x[NE2_94], zn, P["N2"], P["S"])); put("t_H94_Zn_N2_H", tors(x[NE2_94], zn, P["N2"], P["HN21"]))
    put("t_H96_Zn_N2_S", tors(x[NE2_96], zn, P["N2"], P["S"])); put("t_H119_Zn_N2_S", tors(x[ND1_119], zn, P["N2"], P["S"]))
    # orientation
    d1 = dist(P["N1"], x[F131CZ]); d2 = dist(P["N1"], x[N67OD1])
    put("d1", d1); put("d2", d2); put("d_N1_N67ND2", dist(P["N1"], x[N67ND2]))
    put("state", 0 if (d1 < 6.7 and d2 > 5.2) else (1 if (d1 > 6.7 and d2 < 5.2) else 2))  # 0 C, 1 R, 2 I
    # fitted frame quantities
    Pf = {n: (P[n] - gx) @ Q + gr for n in HEAVY11}
    put("lig_rmsd_X", lig_rmsd(Pf, ref_lig)); put("lig_rmsd_R", lig_rmsd(Pf, ligR))
    ax = Pf["C1"] - Pf["C4"]; ax /= norm(ax)
    put("axis_to_X", float(np.degrees(np.arccos(np.clip(ax @ ref_axis, -1, 1)))))
    put("axis_to_R", float(np.degrees(np.arccos(np.clip(ax @ refR_axis, -1, 1)))))
    _, nf, _ = plane(np.array([Pf[n] for n in RING]))
    put("normal_to_X", float(np.degrees(np.arccos(min(1.0, abs(float(nf @ ref_nrm)))))))
    put("normal_to_R", float(np.degrees(np.arccos(min(1.0, abs(float(nf @ refR_nrm)))))))
    put("a_ZnN2_axis", float(np.degrees(np.arccos(np.clip(((P["N2"] - zn) / dist(P["N2"], zn)) @ ((P["C1"] - P["C4"]) / dist(P["C1"], P["C4"])), -1, 1)))))
    znf = (zn - gx) @ Q + gr
    put("zn_disp_fit", dist(znf, zn_x))
    put("N1_fit", (P["N1"] - gx) @ Q + gr)
    # contacts of the ring (C1..C6, N1) with residues
    rp = x[ringN1]
    put("contacts", [float(np.sqrt(((rp[:, None, :] - x[ci][None, :, :]) ** 2).sum(-1)).min()) for ci in contact_idx])
    # partners of the aniline N among protein N and O
    dN = np.sqrt(((x[polar] - P["N1"]) ** 2).sum(1))
    order = np.argsort(dN)[:4]
    pi = [-1] * 4; pc = [0] * 4
    for q, o in enumerate(order):
        if dN[o] > 3.5: break
        pi[q] = int(o)
        code = 0
        a1 = max(ang(P["N1"], P["HN11"], x[polar[o]]), ang(P["N1"], P["HN12"], x[polar[o]]))
        if a1 >= 135.0: code |= 1          # N1 donates
        for h in polar_H[o]:
            if ang(x[polar[o]], x[h], P["N1"]) >= 135.0: code |= 2   # N1 accepts
        pc[q] = code
    partner_idx.append(pi); partner_code.append(pc)
    put("d_ligcen_prot_min", float(np.sqrt(((x[prot_heavy] - np.mean([P[n] for n in HEAVY11], 0)) ** 2).sum(1)).min()))

def run_solute(files, is_eq):
    for p in files:
        tag = tag_of(p)
        u = mda.Universe(f"{rd}/solute.pdb", p)
        assert u.atoms.n_atoms == NAT
        n = u.trajectory.n_frames
        for i, ts in enumerate(u.trajectory):
            if is_eq: t = eq_start[tag] + (i + 1) * 10.0
            else: t = (int(tag[3:]) - 1) * 5000.0 + (i + 1) * 10.0
            solute_frame(ts.positions, ts.dimensions, t, tag)
        nfr.setdefault("solute", {})[tag] = n
nfr = {}
run_solute(eq_s, True); run_solute(seg_s, False)

# ---------------------------------------------------------------- full frames
uf0 = mda.Universe(f"{rd}/full.pdb")
fn, frn, frs = uf0.atoms.names, uf0.atoms.resnames, uf0.atoms.resids
assert list(fn[:NAT]) == list(names) and list(frn[:NAT]) == list(resnames)
WO = np.where((frn == "WAT") & (fn == "O"))[0]
assert list(fn[WO + 1]) == ["H1"] * len(WO) and list(fn[WO + 2]) == ["H2"] * len(WO)
NA = np.where((frn == "Na+") | (fn == "Na+") | (fn == "NA"))[0]
CL = np.where((frn == "Cl-") | (fn == "Cl-") | (fn == "CL"))[0]
if list(sel["water_O"]) != list(WO): sel_mismatch.append("water_O")
SITES = {}
for w, p in WSITES.items(): SITES[f"W{w}"] = p
for n in RING: SITES[f"X_{n}"] = ref_lig[n]
SITES["X_N1"] = ref_lig["N1"]
for n in RING: SITES[f"R_{n}"] = ligR[n]
SITES["R_N1"] = ligR["N1"]
site_names = list(SITES); site_ref = np.array([SITES[k] for k in site_names])
F = {}
def putf(k, v): F.setdefault(k, []).append(v)
def full_frame(x, dims, t_ps, stage, do_image):
    x = x.astype(float)
    H = triclinic_vectors(dims).astype(float); Hinv = np.linalg.inv(H)
    putf("t_ps", t_ps); putf("stage", stage)
    zn = x[ZN]
    ca = x[CA]
    putf("ca_step_max", float(np.sqrt(((ca[1:] - ca[:-1]) ** 2).sum(1)).max()))
    Q, gx, gr, rms = kabsch(ca, ref_ca)
    putf("ca_rmsd", rms)
    sites = (site_ref - gr) @ Q.T + gx
    # everything is placed at its image nearest to the zinc; the zinc itself is taken as it is, and the
    # solute atoms used here are imaged to the zinc as well (minimum image)
    def near_zn(idx):
        v, _ = min_image(x[idx] - zn, H, Hinv)
        return zn + v
    w = near_zn(WO)
    dzn = np.sqrt(((w - zn) ** 2).sum(1))
    kz = int(dzn.argmin())
    putf("d_Zn_Ow_min", float(dzn[kz])); putf("i_Zn_Ow_min", int(WO[kz]))
    putf("n_Ow_Zn_2p6", int((dzn <= 2.6).sum())); putf("n_Ow_Zn_3p0", int((dzn <= 3.0).sum()))
    putf("n_Ow_Zn_3p5", int((dzn <= 3.5).sum())); putf("n_Ow_Zn_5p0", int((dzn <= 5.0).sum())); putf("n_Ow_Zn_8p0", int((dzn <= 8.0).sum()))
    po = near_zn(prot_O)
    dpo = np.sqrt(((po - zn) ** 2).sum(1)); kp = int(dpo.argmin())
    putf("d_Zn_Oprot_min", float(dpo[kp])); putf("i_Zn_Oprot_min", kp)
    putf("n_Oprot_Zn_2p6", int((dpo <= 2.6).sum())); putf("n_Oprot_Zn_3p0", int((dpo <= 3.0).sum()))
    lg = {n: near_zn(np.array([i]))[0] for n, i in lig.items()}
    putf("d_Zn_Osulf_min", min(dist(zn, lg["O1"]), dist(zn, lg["O2"])))
    putf("d_Zn_N2_full", dist(zn, lg["N2"]))
    # sites of the crystal
    sites_i = zn + min_image(sites - zn, H, Hinv)[0]
    near = w[dzn < 25.0]; near_idx = WO[dzn < 25.0]
    D = np.sqrt(((near[:, None, :] - sites_i[None, :, :]) ** 2).sum(-1))
    kmin = D.argmin(0)
    putf("site_dmin", D.min(0)); putf("site_imin", near_idx[kmin])
    xr = [site_names.index(f"X_{n}") for n in RING]; rr = [site_names.index(f"R_{n}") for n in RING]
    dX = D[:, xr].min(1); dR = D[:, rr].min(1)
    putf("n_w_Xring_2p5", int((dX <= 2.5).sum())); putf("n_w_Rring_2p5", int((dR <= 2.5).sum()))
    putf("d_w_Xring_min", float(dX.min())); putf("d_w_Rring_min", float(dR.min()))
    dXw = D[:, xr + [site_names.index("X_N1")]].min(1)
    putf("n_w_Xring_wide_3p5", int((dXw <= 3.5).sum()))
    # water at the ligand's polar atoms
    for an in ("N1", "N2", "O1", "O2"):
        dl = np.sqrt(((near - lg[an]) ** 2).sum(1))
        putf(f"n_w_{an}_3p5", int((dl <= 3.5).sum()))
        if an == "N1":
            nd = na = 0
            for q in np.where(dl <= 3.5)[0]:
                o = near[q]; oi = near_idx[q]
                h1 = o + min_image((x[oi + 1] - o)[None], H, Hinv)[0][0]; h2 = o + min_image((x[oi + 2] - o)[None], H, Hinv)[0][0]
                if max(ang(lg["N1"], lg["HN11"], o), ang(lg["N1"], lg["HN12"], o)) >= 135.0: nd += 1
                if max(ang(o, h1, lg["N1"]), ang(o, h2, lg["N1"])) >= 135.0: na += 1
            putf("n_w_N1_donated", nd); putf("n_w_N1_accepted", na)
    # water bridging: within 3.2 of both N2 and anything? kept simple: water within 3.5 of Thr200 OG1 and of the ring
    # ions
    for nm, idx in (("Na", NA), ("Cl", CL)):
        if len(idx):
            pi_ = near_zn(idx)
            putf(f"d_Zn_{nm}_min", float(np.sqrt(((pi_ - zn) ** 2).sum(1)).min()))
            putf(f"d_N1_{nm}_min", float(np.sqrt(((pi_ - lg["N1"]) ** 2).sum(1)).min()))
    # distance of the protein from its own images (heavy atoms), on every 10th production frame
    if do_image:
        ph = x[prot_heavy]
        tree = cKDTree(ph)
        best = 1e9
        for s in SHIFTS:
            if not s.any(): continue
            dq, _ = tree.query(ph + s @ H, k=1)
            best = min(best, float(dq.min()))
        putf("prot_image_min", best)
    else:
        putf("prot_image_min", np.nan)

def run_full(files, is_eq):
    cnt = 0
    for p in files:
        tag = tag_of(p)
        u = mda.Universe(f"{rd}/full.pdb", p)
        n = u.trajectory.n_frames
        for i, ts in enumerate(u.trajectory):
            if is_eq: t = eq_start[tag] + (i + 1) * 100.0
            else: t = (int(tag[3:]) - 1) * 5000.0 + (i + 1) * 100.0
            full_frame(ts.positions, ts.dimensions, t, tag, (not is_eq) and (cnt % 10 == 9))
            cnt += 1
        nfr.setdefault("full", {})[tag] = n
run_full(eq_f, True); run_full(seg_f, False)

out = {}
for k, v in S.items(): out["s_" + k] = np.array(v)
for k, v in F.items(): out["f_" + k] = np.array(v)
out["s_partner_idx"] = np.array(partner_idx, dtype=np.int16); out["s_partner_code"] = np.array(partner_code, dtype=np.int8)
meta = dict(model=model, run=run, start=start, t_release_ps=T_RELEASE, eq_tags=eq_tags, eq_len=eq_len, frames=nfr,
            polar_label=polar_label, prot_O_label=prot_O_label, contact_res=CONTACT_RES, site_names=site_names,
            sel_mismatch=sel_mismatch, n_water=int(len(WO)), n_Na=int(len(NA)), n_Cl=int(len(CL)), natoms_full=int(uf0.atoms.n_atoms))
out["meta"] = np.array(json.dumps(meta))
np.savez_compressed(f"{OUT}/data_{model}_{run}.npz", **out)
print(model, run, "solute frames", len(S["t_ps"]), "full frames", len(F["t_ps"]), "release", T_RELEASE, "sel_mismatch", sel_mismatch)
