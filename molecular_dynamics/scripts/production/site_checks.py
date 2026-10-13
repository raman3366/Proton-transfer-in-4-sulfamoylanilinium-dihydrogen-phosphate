#!/usr/bin/env python3
"""Checks on the raw frames for md/RESULTS.md, made after the two audits and the physical review of 8 Oct 2026.

Own code (MDAnalysis reads the DCD files; atoms are chosen by residue number and name through
md/system/residue_map.csv). No code and no output of the two analysis programs is used to form a number; the
series of the two programs are read only to set this reading beside theirs (the largest difference is written
down, and the script stops if a frame is classed differently). Nothing under /root/cqc-sfa-md/runs is changed.

    --part solute   solute frames (10 ps) of the 13 runs, equilibration and production:
                      frames within the tolerance of a threshold (md/analysis/SPEC.md, section 14);
                    and for the eleven long runs:
                      His64 (Zn to NE2, chi1), the zinc out of the plane of the His119 ring, S-C4, O1-S-C4,
                      Zn-N2-H, the sum of the three angles at N2, the torsions of S and of H about Zn-N2,
                      the ligand RMSD to the docked rotated pose (md/system/ligand_R.pdb)
                    -> md/results/site_checks_solute.json
    --part full     all-atom frames (100 ps) of the 13 runs:
                      frames within the tolerance of a threshold;
                    and for the eleven long runs:
                      water at the places of crystal waters 423, 516, 553, 568 (md/system/waters_all.pdb),
                      a water that donates a hydrogen bond to N2, water oxygens within 8 A of the zinc
                    -> md/results/site_checks_full.json

Definitions that are not those of the specification are those of the physical review
(md/results/independent/physical_review/REVIEW.md, section 0), so that its numbers can be set beside:
a crystal-water place counts as occupied when a water oxygen lies within 1.4 A of it after the fit of the CA
atoms on reference.pdb; "His64 deep in" is chi1 between 0 and 120 degrees with Zn to NE2 at most 7.0 A (a limit
that the review chose after looking at the data; 6.5 and 7.5 A are given beside it).

    /root/cqc-sfa-tools/md/bin/python md/production/site_checks.py --part solute
    /root/cqc-sfa-tools/md/bin/python md/production/site_checks.py --part full
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import MDAnalysis as mda  # noqa: E402
from MDAnalysis.lib.distances import distance_array, minimize_vectors  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MD = os.path.dirname(HERE)
RUNS = "/root/cqc-sfa-md/runs"
OUT = os.path.join(REPO_MD, "results")
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))
NINE = ["bonded_2e/" + r for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")]
SENS = ["bonded_2s/R1", "bonded_2e_n2fit/R1"]
VAL = ["bonded_2e/X0", "bonded_2s/X0"]
ELEVEN = NINE + SENS
ALL = ELEVEN + VAL
RING = ["C1", "C2", "C3", "C4", "C5", "C6"]
HEAVY = RING + ["N1", "N2", "O1", "O2", "S"]
LABELLINGS = [{}, {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3"}, {"O1": "O2", "O2": "O1"},
              {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3", "O1": "O2", "O2": "O1"}]
# thresholds and tolerances of md/analysis/SPEC.md, section 14
M1, M2, HB_D, HB_A, ZN_O, RING_CUT, WAT_CUT = 6.7, 5.2, 3.5, 135.0, 2.6, 2.5, 3.5
TOL_D, TOL_A, TOL_V = 0.001, 0.01, 0.001
WATERS = (423, 516, 553, 568)
OCC = 1.4                     # A; a crystal-water place is occupied (review, section 0)
DEEP = (6.5, 7.0, 7.5)        # A; Zn to NE2(His64)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def J(path):
    with open(path) as fh:
        return json.load(fh)


def rows_of(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def tleap_index():
    """deposited residue number -> residue number of the built system"""
    return {int(r["deposited_number"]): int(r["tleap_index"]) for r in rows_of(os.path.join(REPO_MD, "system/residue_map.csv"))}


def pdb_atoms(path):
    out = []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            out.append((ln[12:16].strip(), ln[17:20].strip(), int(ln[22:26]), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return out


def norm(v):
    return np.linalg.norm(v, axis=-1)


def angle(a, b, c):
    u, v = a - b, c - b
    cs = np.einsum("fi,fi->f", u, v) / norm(u) / norm(v)
    return np.degrees(np.arccos(np.clip(cs, -1.0, 1.0)))


def torsion(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / norm(b1)[:, None]
    v = b0 - np.einsum("fi,fi->f", b0, b1n)[:, None] * b1n
    w = b2 - np.einsum("fi,fi->f", b2, b1n)[:, None] * b1n
    return np.degrees(np.arctan2(np.einsum("fi,fi->f", np.cross(b1n, v), w), np.einsum("fi,fi->f", v, w)))


def circ_mean(a):
    a = np.radians(np.asarray(a, float))
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())))


def kabsch(X, ref):
    """Rotations Q and centroids so that (x - gx) @ Q + gr lies on ref (row vectors), without reflection.
    X: (F, M, 3); ref: (M, 3)."""
    gx = X.mean(axis=1, keepdims=True)
    gr = ref.mean(axis=0)
    H = np.einsum("fmi,mj->fij", X - gx, ref - gr)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    D = np.zeros((X.shape[0], 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    return np.einsum("fij,fjk,fkl->fil", U, D, Vt), gx, gr


def reference(system_dir):
    at = pdb_atoms(os.path.join(system_dir, "reference.pdb"))
    ca = np.array([a[3] for a in at if a[0] == "CA"])
    lig = {a[0]: a[3] for a in at if a[1] == "SN1"}
    assert ca.shape == (258, 3) and set(lig) == set(HEAVY)
    return ca, lig


def crystal_sites():
    """Places in the frame of 6RL9: the four waters removed for start R, the docked rotated pose."""
    wat = {a[2]: a[3] for a in pdb_atoms(os.path.join(REPO_MD, "system/waters_all.pdb")) if a[0] == "O"}
    ligR = {a[0]: a[3] for a in pdb_atoms(os.path.join(REPO_MD, "system/ligand_R.pdb")) if a[0] in HEAVY}
    ligX = {a[0]: a[3] for a in pdb_atoms(os.path.join(REPO_MD, "system/ligand_X.pdb")) if a[0] in HEAVY}
    assert set(ligR) == set(HEAVY) and set(ligX) == set(HEAVY)
    return {w: wat[w] for w in WATERS}, ligR, ligX


def steps_of(man):
    """(tag, kind, stage, extension, solute frames, all-atom frames) of every step in the order of the manifest"""
    out = []
    for s in man["equilibration"]["stages"]:
        tag = s["tag"]
        stage = int(tag[2:4])
        ext = int(tag.split("x")[1]) if "x" in tag else 0
        out.append((tag, "equilibration", stage, ext, s["frames"]["solute"], s["frames"]["full"]))
    for s in man["production"]["segments"]:
        out.append((s["tag"], "production", None, 0, s["frames"]["solute"], s["frames"]["full"]))
    return out


def asis(value, decimals=None):
    """Numbers are written as computed. They are rounded once, where they are printed (a number rounded here and
    again in the record can come out one unit of the last digit off)."""
    return float(value)


def near(values, threshold, tol):
    return np.abs(np.asarray(values) - threshold) <= tol


# ================================================================================================ solute
def solute_run(lab, idx_of):
    d = os.path.join(RUNS, lab)
    man = J(os.path.join(d, "manifest.json"))
    steps = steps_of(man)
    u = mda.Universe(os.path.join(d, "solute.pdb"), [os.path.join(d, "solute_%s.dcd" % s[0]) for s in steps])
    n_expected = sum(s[4] for s in steps)
    assert len(u.trajectory) == n_expected, (lab, len(u.trajectory), n_expected)

    def one(resname, resid, name):
        sel = u.select_atoms("resname %s and resid %d and name %s" % (resname, resid, name))
        assert len(sel) == 1, (lab, resname, resid, name, len(sel))
        return int(sel.indices[0])

    lig_res = int(u.select_atoms("resname SN1").resids[0])
    zn_res = int(u.select_atoms("resname ZN1").resids[0])
    L = {n: one("SN1", lig_res, n) for n in HEAVY + ["HN21"]}
    A = {"ZN": one("ZN1", zn_res, "ZN"),
         "H94NE2": one("HD1", idx_of[94], "NE2"), "H96NE2": one("HD2", idx_of[96], "NE2"),
         "F131CZ": one("PHE", idx_of[131], "CZ"), "N67OD1": one("ASN", idx_of[67], "OD1"),
         "T199OG1": one("THR", idx_of[199], "OG1")}
    H119 = [one("HE1", idx_of[119], n) for n in ("CG", "ND1", "CD2", "CE1", "NE2")]
    H64 = {n: one("HIE", idx_of[64], n) for n in ("N", "CA", "CB", "CG", "NE2")}
    CA = np.array(sorted(set(int(i) for i in u.select_atoms("name CA").indices)))
    assert len(CA) == 258, (lab, len(CA))
    keep = sorted(set(list(L.values()) + list(A.values()) + H119 + list(H64.values()) + list(CA)))
    pos_of = {i: k for k, i in enumerate(keep)}
    X = np.empty((n_expected, len(keep), 3), dtype=np.float64)
    for k, ts in enumerate(u.trajectory):
        X[k] = ts.positions[keep]
    g = lambda i: X[:, pos_of[i], :]
    zn = g(A["ZN"])
    P = {n: g(i) for n, i in L.items()}
    # the molecule is written whole: guards
    ca = X[:, [pos_of[i] for i in CA], :]
    guard = {"largest_neighbour_CA_distance_A": float(norm(ca[:, 1:] - ca[:, :-1]).max()),
             "largest_raw_Zn_N2_A": float(norm(zn - P["N2"]).max()),
             "largest_raw_Zn_H94NE2_A": float(norm(zn - g(A["H94NE2"])).max())}
    assert guard["largest_neighbour_CA_distance_A"] < 4.5 and guard["largest_raw_Zn_N2_A"] < 2.6 and guard["largest_raw_Zn_H94NE2_A"] < 2.6, (lab, guard)

    part = np.array(["e"] * sum(s[4] for s in steps if s[1] == "equilibration") + ["p"] * sum(s[4] for s in steps if s[1] == "production"))
    stage = np.concatenate([np.full(s[4], -1 if s[2] is None else s[2]) for s in steps])
    ext = np.concatenate([np.full(s[4], s[3]) for s in steps])
    prod = part == "p"
    released = (part == "e") & (stage >= 4) & (ext == 0)
    n_prod = int(prod.sum())

    d1 = norm(P["N1"] - g(A["F131CZ"]))
    d2 = norm(P["N1"] - g(A["N67OD1"]))
    hbd = norm(P["N2"] - g(A["T199OG1"]))
    hba = angle(P["N2"], P["HN21"], g(A["T199OG1"]))
    zno = np.minimum(norm(zn - P["O1"]), norm(zn - P["O2"]))
    unit = lambda v: v / norm(v)[:, None]
    V = np.einsum("fi,fi->f", unit(zn - P["N2"]), np.cross(unit(P["HN21"] - P["N2"]), unit(P["S"] - P["N2"])))
    state = np.where((d1 < M1) & (d2 > M2), "C", np.where((d1 > M1) & (d2 < M2), "R", "I"))
    hb = (hbd <= HB_D) & (hba >= HB_A)
    fifth = zno <= ZN_O
    sign = np.sign(V).astype(int)

    # this reading beside the series of the two programs
    beside = {}
    for sub in ("analysis", "analysis_recheck"):
        eq = rows_of(os.path.join(d, sub, "timeseries_equil_solute.csv"))
        pr = rows_of(os.path.join(d, sub, "timeseries_solute.csv"))
        rows = eq + pr
        assert len(rows) == n_expected, (lab, sub, len(rows))
        col = lambda name: np.array([float(r[name]) for r in rows])
        worst = {"d1_A": float(np.abs(col("d1_N1_F131CZ_A") - d1).max()), "d2_A": float(np.abs(col("d2_N1_N67OD1_A") - d2).max()),
                 "d_N2_T199OG1_A": float(np.abs(col("d_N2_T199OG1_A") - hbd).max()),
                 "ang_N2_HN2_T199OG1_deg": float(np.abs(col("ang_N2_HN2_T199OG1_deg") - hba).max()),
                 "d_Zn_Osulfonyl_min_A": float(np.abs(col("d_Zn_Osulfonyl_min_A") - zno).max()),
                 "n2_volume": float(np.abs(col("n2_volume") - V).max())}
        ndiff = {"state": int((np.array([r["state"] for r in rows]) != state).sum()),
                 "hb_N2H_T199OG1": int((np.array([int(r["hb_N2H_T199OG1"]) for r in rows]) != hb.astype(int)).sum()),
                 "fifth_sulfonylO": int((np.array([int(r["fifth_sulfonylO"]) for r in rows]) != fifth.astype(int)).sum()),
                 "n2_sign": int((np.array([int(r["n2_sign"]) for r in rows]) != sign).sum())}
        if any(ndiff.values()) or max(worst[k] for k in worst if k != "ang_N2_HN2_T199OG1_deg") > 1.0e-4 or worst["ang_N2_HN2_T199OG1_deg"] > 1.0e-3:
            raise SystemExit("%s: this reading and %s differ: %s %s" % (lab, sub, worst, ndiff))
        beside[sub] = {"largest_difference": {k: asis(v, 6) for k, v in worst.items()}, "frames_classed_differently": ndiff}

    def frames_near(mask):
        return [{"part": "production" if part[k] == "p" else "equilibration",
                 "stage": None if stage[k] < 0 else int(stage[k]), "frame_of_the_joined_reading": int(k)} for k in np.where(mask)[0]]

    nd1, nd2, nhd, nha, nzo, nv = near(d1, M1, TOL_D), near(d2, M2, TOL_D), near(hbd, HB_D, TOL_D), near(hba, HB_A, TOL_A), near(zno, ZN_O, TOL_D), np.abs(V) < TOL_V
    any_near = nd1 | nd2 | nhd | nha | nzo | nv
    rec = {"run": lab, "n_solute_frames": int(n_expected), "n_production_frames": n_prod,
           "n_equilibration_frames": int((part == "e").sum()), "n_released_equilibration_frames": int(released.sum()),
           "guard": {k: asis(v, 3) for k, v in guard.items()},
           "beside_the_two_programs": beside,
           "threshold": {
               "what": "solute frames, equilibration and production, whose value lies within the tolerance of a threshold "
                       "(SPEC.md, section 14): d1 6.7 A, d2 5.2 A, N2...OG1(Thr199) 3.5 A (0.001 A each), N2-H...OG1 135 deg "
                       "(0.01 deg), nearer sulfonyl oxygen to the zinc 2.6 A (0.001 A), |signed volume at N2| below 0.001",
               "n_d1": int(nd1.sum()), "n_d2": int(nd2.sum()), "n_hb_distance": int(nhd.sum()), "n_hb_angle": int(nha.sum()),
               "n_sulfonyl_O": int(nzo.sum()), "n_n2_sign": int(nv.sum()), "n_frames_any": int(any_near.sum()),
               "n_frames_any_production": int((any_near & prod).sum()),
               "nearest_margin": {"d1_A": float(np.abs(d1 - M1).min()), "d2_A": float(np.abs(d2 - M2).min()),
                                  "hb_distance_A": float(np.abs(hbd - HB_D).min()), "hb_angle_deg": float(np.abs(hba - HB_A).min()),
                                  "sulfonyl_O_A": float(np.abs(zno - ZN_O).min()), "abs_n2_volume": float(np.abs(V).min())},
               "frames": frames_near(any_near)}}
    if lab in VAL:
        return rec

    # ---------------------------------------------------------------- the eleven long runs
    last80 = prod.copy()
    last80[np.where(prod)[0][:2000]] = False            # production from 20 to 100 ns
    # His64
    dz64 = norm(zn - g(H64["NE2"]))
    chi1 = torsion(g(H64["N"]), g(H64["CA"]), g(H64["CB"]), g(H64["CG"]))
    gplus = (chi1 > 0.0) & (chi1 < 120.0)
    rec["his64"] = {"what": "Zn to NE2(His64) and chi1 (N-CA-CB-CG); 'deep in' = chi1 between 0 and 120 deg with the distance "
                            "at most the limit; shares in per cent of the solute frames of the production from 20 to 100 ns",
                    "n_frames_20_to_100_ns": int(last80.sum()),
                    "deep_in_share_20_to_100_ns": {"%.1f" % c: asis(100.0 * float(((dz64 <= c) & gplus)[last80].mean()), 1) for c in DEEP},
                    "chi1_gauche_plus_share_production": asis(100.0 * float(gplus[prod].mean()), 1),
                    "Zn_NE2_mean_first_10_ns_A": asis(float(dz64[np.where(prod)[0][:1000]].mean()), 2),
                    "Zn_NE2_mean_20_to_100_ns_A": asis(float(dz64[last80].mean()), 2),
                    "Zn_NE2_mean_while_deep_in_7.0_A": (asis(float(dz64[last80 & (dz64 <= 7.0) & gplus].mean()), 2)
                                                        if (last80 & (dz64 <= 7.0) & gplus).any() else None)}
    # steady departures that no gate tests (plan, 7 Oct 2026, 12:10 AM IST, point 8)
    ring = X[:, [pos_of[i] for i in H119], :]
    c = ring.mean(axis=1)
    _, _, vt = np.linalg.svd(ring - c[:, None, :])
    nrm = vt[:, 2, :]
    # one sense of the normal in every frame: that of (ND1 - CG) x (CD2 - CG)
    sense = np.sign(np.einsum("fi,fi->f", nrm, np.cross(ring[:, 1] - ring[:, 0], ring[:, 2] - ring[:, 0])))
    oop = np.einsum("fi,fi->f", zn - c, nrm * sense[:, None])
    rec["departures"] = {"what": "means over the solute frames of the production; the zinc out of the plane of the His119 ring "
                                 "(least-squares plane of CG, ND1, CD2, CE1, NE2) is the mean of the signed distance, given "
                                 "without its sign, and the mean of the unsigned distance beside it",
                         "zinc_out_of_the_plane_of_the_His119_ring_A": asis(abs(float(oop[prod].mean())), 3),
                         "zinc_out_of_the_plane_of_the_His119_ring_mean_unsigned_A": asis(float(np.abs(oop[prod]).mean()), 3),
                         "d_S_C4_A": asis(float(norm(P["S"] - P["C4"])[prod].mean()), 4),
                         "ang_O1_S_C4_deg": asis(float(angle(P["O1"], P["S"], P["C4"])[prod].mean()), 2),
                         "ang_S_N2_Zn_deg": asis(float(angle(P["S"], P["N2"], zn)[prod].mean()), 2),
                         "d_S_N2_A": asis(float(norm(P["S"] - P["N2"])[prod].mean()), 4)}
    # the nitrogen
    a_znh = angle(zn, P["N2"], P["HN21"])
    sum_ang = angle(zn, P["N2"], P["S"]) + a_znh + angle(P["S"], P["N2"], P["HN21"])
    from_release = released | prod
    rec["nitrogen"] = {"what": "signed volume V at N2 (Zn, H, S) and the sum of the three angles at N2",
                       "production": {"V_min": asis(float(V[prod].min()), 4), "V_max": asis(float(V[prod].max()), 4),
                                      "n_abs_V_below_0.2": int((np.abs(V[prod]) < 0.2).sum()),
                                      "n_sum_of_angles_above_355_deg": int((sum_ang[prod] > 355.0).sum()),
                                      "sum_of_angles_mean_deg": asis(float(sum_ang[prod].mean()), 2),
                                      "ang_Zn_N2_H_mean_deg": asis(float(a_znh[prod].mean()), 2)},
                       "from_the_release": {"n_frames": int(from_release.sum()), "V_min": asis(float(V[from_release].min()), 4),
                                            "V_max": asis(float(V[from_release].max()), 4),
                                            "n_changes_of_sign": int((sign[from_release][1:] != sign[from_release][:-1]).sum())}}
    # the ligand against the docked rotated pose
    refca, reflig = reference(man["system_dir"])
    _, ligR, ligX = crystal_sites()
    assert max(float(np.linalg.norm(ligX[n] - reflig[n])) for n in HEAVY) < 1.0e-3
    Q, gx, gr = kabsch(ca, refca)
    fit = lambda arr: np.einsum("fi,fij->fj", arr - gx[:, 0, :], Q) + gr

    def lig_rmsd(ref):
        best = None
        for pm in LABELLINGS:
            s = sum(((fit(P[n]) - ref[pm.get(n, n)]) ** 2).sum(axis=1) for n in HEAVY)
            r = np.sqrt(s / len(HEAVY))
            best = r if best is None else np.minimum(best, r)
        return best

    rX, rR = lig_rmsd(reflig), lig_rmsd(ligR)
    held = (part == "e") & ~released
    rec["ligand_against_the_docked_rotated_pose"] = {
        "what": "RMSD of the 11 heavy atoms of the ligand to md/system/ligand_R.pdb after the fit of the CA atoms on "
                "reference.pdb, smallest of the four equivalent labellings",
        "production_smallest_A": asis(float(rR[prod].min()), 3), "production_mean_A": asis(float(rR[prod].mean()), 3),
        "production_frames_within_1.0_A": int((rR[prod] <= 1.0).sum()),
        "held_part_mean_A": asis(float(rR[held].mean()), 3),
        "rmsd_to_the_crystal_pose_production_mean_A": asis(float(rX[prod].mean()), 4)}
    # the rotated-like frames: which configuration at N2, and how far from the docked pose
    rs = (state == "R")
    rec["rotated_like_frames"] = {
        "from_the_release": {"n": int((rs & from_release).sum()), "n_with_V_positive": int((rs & from_release & (V > 0)).sum()),
                             "n_within_1.0_A_of_the_docked_rotated_pose": int((rs & from_release & (rR <= 1.0)).sum()),
                             "d1_mean_A": asis(float(d1[rs & from_release].mean()), 2) if (rs & from_release).any() else None,
                             "d2_mean_A": asis(float(d2[rs & from_release].mean()), 2) if (rs & from_release).any() else None},
        "production": {"n": int((rs & prod).sum()), "n_with_V_positive": int((rs & prod & (V > 0)).sum()),
                       "n_crystal_like_with_V_negative": int(((state == "C") & prod & (V < 0)).sum())}}
    # the change of the configuration at N2 in a run from start R
    if lab.endswith("/R1") or lab.endswith("/R2") or lab.endswith("/R3"):
        tS = torsion(g(A["H94NE2"]), zn, P["N2"], P["S"])
        tH = torsion(g(A["H94NE2"]), zn, P["N2"], P["HN21"])
        rel_idx = np.where(released)[0]
        k_first = int(np.where(V[rel_idx] > 0)[0][0])
        before = rel_idx[:k_first]
        assert (V[before] < 0).all()
        ev = {"what": "torsions NE2(His94)-Zn-N2-S and NE2(His94)-Zn-N2-H (circular means): the released frames before the "
                      "first frame with V positive, and the production",
              "first_frame_with_V_positive_ps_from_release": (k_first + 1) * 10,
              "n_frames_before": int(len(before)),
              "torsion_S_before_deg": asis(circ_mean(tS[before]), 1), "torsion_S_production_deg": asis(circ_mean(tS[prod]), 1),
              "torsion_H_before_deg": asis(circ_mean(tH[before]), 1), "torsion_H_production_deg": asis(circ_mean(tH[prod]), 1),
              "sum_of_angles_before_mean_deg": asis(float(sum_ang[before].mean()), 1),
              "ang_S_N2_Zn_before_mean_deg": asis(float(angle(P["S"], P["N2"], zn)[before].mean()), 1)}
        if lab == "bonded_2e/R1":
            k = rel_idx[52]                                # 530 ps after the release
            ev["frame_at_530_ps"] = {"d1_A": asis(float(d1[k]), 2), "d2_A": asis(float(d2[k]), 2), "V": asis(float(V[k]), 2),
                                     "sum_of_angles_deg": asis(float(sum_ang[k]), 1), "state": str(state[k])}
        rec["change_at_N2"] = ev
    return rec


# ================================================================================================ full
def full_run(lab, idx_of):
    d = os.path.join(RUNS, lab)
    man = J(os.path.join(d, "manifest.json"))
    steps = steps_of(man)
    u = mda.Universe(os.path.join(d, "full.pdb"), [os.path.join(d, "full_%s.dcd" % s[0]) for s in steps])
    n_expected = sum(s[5] for s in steps)
    assert len(u.trajectory) == n_expected, (lab, len(u.trajectory), n_expected)
    names, resn = u.atoms.names, u.atoms.resnames
    WO = np.where((resn == "WAT") & (names == "O"))[0]
    assert len(WO) == man["system"]["n_water_molecules"] and list(names[WO + 1]) == ["H1"] * len(WO) and list(names[WO + 2]) == ["H2"] * len(WO)
    not_protein = {"WAT", "Na+", "Cl-", "ZN1", "SN1"}
    PO = np.array([i for i in range(len(names)) if resn[i] not in not_protein and names[i].startswith("O")])
    CA = np.array([i for i in range(len(names)) if names[i] == "CA" and resn[i] not in not_protein])
    assert len(CA) == 258
    zn_i = int(np.where(resn == "ZN1")[0][0])
    L = {str(names[i]): int(i) for i in np.where(resn == "SN1")[0]}
    refca, reflig = reference(man["system_dir"])
    wsites, ligR, ligX = crystal_sites()
    # the files of the crystal frame belong together: water 553 lies 2.93 A from N2 of the crystal pose
    assert abs(float(np.linalg.norm(wsites[553] - reflig["N2"])) - 2.93) < 0.01
    ring_ref = np.array([reflig[n] for n in RING])
    wide_ref = np.array([reflig[n] for n in RING + ["N1"]])
    wat_ref = np.array([wsites[w] for w in WATERS])

    part = np.array(["e"] * sum(s[5] for s in steps if s[1] == "equilibration") + ["p"] * sum(s[5] for s in steps if s[1] == "production"))
    prod = part == "p"
    F = {k: [] for k in ("zn_ow", "zn_op", "zn_os", "n_ring", "d_ring", "n_n1", "n_wide", "near_ring_cut", "near_n1_cut", "near_wide_cut",
                         "occ", "d_occ", "donor", "d_donor", "n8")}
    for ts in u.trajectory:
        x = ts.positions.astype(np.float64)
        box = ts.dimensions
        zn = x[zn_i][None, :]
        W = x[WO]
        dzn = distance_array(zn, W, box=box)[0]
        F["zn_ow"].append(dzn.min())
        F["zn_op"].append(distance_array(zn, x[PO], box=box)[0].min())
        F["zn_os"].append(distance_array(zn, x[[L["O1"], L["O2"]]], box=box)[0].min())
        F["n8"].append(int((dzn <= 8.0).sum()))
        Q, gx, gr = kabsch(x[CA][None], refca)
        back = lambda pts: (pts - gr) @ Q[0].T + gx[0, 0]           # places of the crystal carried into the frame
        dr = distance_array(back(ring_ref), W, box=box).min(axis=0)
        F["n_ring"].append(int((dr <= RING_CUT).sum()))
        F["d_ring"].append(dr.min())
        F["near_ring_cut"].append(int((np.abs(dr - RING_CUT) <= TOL_D).sum()))
        dn1 = distance_array(x[L["N1"]][None, :], W, box=box)[0]
        F["n_n1"].append(int((dn1 <= WAT_CUT).sum()))
        F["near_n1_cut"].append(int((np.abs(dn1 - WAT_CUT) <= TOL_D).sum()))
        dw = distance_array(back(wide_ref), W, box=box).min(axis=0)
        F["n_wide"].append(int((dw <= WAT_CUT).sum()))
        F["near_wide_cut"].append(int((np.abs(dw - WAT_CUT) <= TOL_D).sum()))
        do = distance_array(back(wat_ref), W, box=box).min(axis=1)
        F["occ"].append(do <= OCC)
        F["d_occ"].append(do)
        # a water that donates to N2: O within 3.5 A of N2 and O-H...N2 at least 135 deg
        n2 = x[L["N2"]]
        dn2 = distance_array(n2[None, :], W, box=box)[0]
        best = None
        for k in np.where(dn2 <= HB_D)[0]:
            o = x[WO[k]]
            for h in (WO[k] + 1, WO[k] + 2):
                hp = x[h]
                # minimum image of both vectors from the hydrogen
                vo = minimize_vectors(np.array([o - hp]), box)[0]
                vn = minimize_vectors(np.array([n2 - hp]), box)[0]
                a = np.degrees(np.arccos(np.clip(float(vo @ vn) / np.linalg.norm(vo) / np.linalg.norm(vn), -1.0, 1.0)))
                if a >= HB_A and (best is None or dn2[k] < best):
                    best = float(dn2[k])
        F["donor"].append(best is not None)
        F["d_donor"].append(np.nan if best is None else best)
    A = {k: np.array(v) for k, v in F.items()}
    any_o = np.minimum(np.minimum(A["zn_ow"], A["zn_op"]), A["zn_os"])

    # this reading beside the all-atom series of the two programs
    beside = {}
    for sub in ("analysis", "analysis_recheck"):
        rows = rows_of(os.path.join(d, sub, "timeseries_equil_full.csv")) + rows_of(os.path.join(d, sub, "timeseries_full.csv"))
        assert len(rows) == n_expected, (lab, sub, len(rows))
        col = lambda name: np.array([float(r[name]) for r in rows])
        worst = {"d_Zn_Owater_min_A": float(np.abs(col("d_Zn_Owater_min_A") - A["zn_ow"]).max()),
                 "d_Zn_Oprotein_min_A": float(np.abs(col("d_Zn_Oprotein_min_A") - A["zn_op"]).max()),
                 "d_Zn_Osulfonyl_min_full_A": float(np.abs(col("d_Zn_Osulfonyl_min_full_A") - A["zn_os"]).max()),
                 "d_water_ring_place_min_A": float(np.abs(col("d_water_ring_place_min_A") - A["d_ring"]).max())}
        ndiff = {"n_water_ring_place": int((col("n_water_ring_place").astype(int) != A["n_ring"]).sum()),
                 "n_water_N1_3p5A": int((col("n_water_N1_3p5A").astype(int) != A["n_n1"]).sum()),
                 "n_water_ring_wide_3p5A": int((col("n_water_ring_wide_3p5A").astype(int) != A["n_wide"]).sum()),
                 "fifth_any": int((col("fifth_any").astype(int) != (any_o <= ZN_O).astype(int)).sum()),
                 "fifth_water": int((col("fifth_water").astype(int) != (A["zn_ow"] <= ZN_O).astype(int)).sum()),
                 "fifth_proteinO": int((col("fifth_proteinO").astype(int) != (A["zn_op"] <= ZN_O).astype(int)).sum())}
        if any(ndiff.values()) or max(worst.values()) > 2.0e-4:
            raise SystemExit("%s: this reading and %s differ in the all-atom frames: %s %s" % (lab, sub, worst, ndiff))
        beside[sub] = {"largest_difference": {k: asis(v, 6) for k, v in worst.items()}, "frames_classed_differently": ndiff}

    nz = near(A["zn_ow"], ZN_O, TOL_D) | near(A["zn_op"], ZN_O, TOL_D) | near(A["zn_os"], ZN_O, TOL_D)
    nr, nn, nw = A["near_ring_cut"] > 0, A["near_n1_cut"] > 0, A["near_wide_cut"] > 0
    any_near = nz | nr | nn | nw
    rec = {"run": lab, "n_all_atom_frames": int(n_expected), "n_production_frames": int(prod.sum()),
           "n_water_O": int(len(WO)), "n_protein_O": int(len(PO)),
           "beside_the_two_programs": beside,
           "threshold": {
               "what": "all-atom frames, equilibration and production, with a value within 0.001 A of a threshold (SPEC.md, "
                       "section 14): the nearest water, protein or sulfonyl oxygen 2.6 A from the zinc; a water oxygen at the "
                       "cut-off of n_water_ring_place (2.5 A), of n_water_N1_3p5A or of n_water_ring_wide_3p5A (3.5 A)",
               "n_zinc_oxygen": int(nz.sum()), "n_ring_place": int(nr.sum()), "n_water_N1": int(nn.sum()), "n_ring_wide": int(nw.sum()),
               "n_frames_any": int(any_near.sum()), "n_frames_any_production": int((any_near & prod).sum()),
               "n_frames_zinc_oxygen_or_ring_place": int((nz | nr).sum()),
               "n_frames_zinc_oxygen_or_ring_place_production": int(((nz | nr) & prod).sum()),
               "nearest_margin_zinc_oxygen_A": float(np.minimum(np.minimum(np.abs(A["zn_ow"] - ZN_O), np.abs(A["zn_op"] - ZN_O)), np.abs(A["zn_os"] - ZN_O)).min())}}
    if lab in VAL:
        return rec
    p = np.where(prod)[0]
    late = p[200:]                                            # production from 20 to 100 ns
    occ = A["occ"]
    rec["crystal_water_places"] = {
        "what": "a water oxygen within %.1f A of the place of the crystal water after the fit of the CA atoms; per cent of the "
                "all-atom frames" % OCC,
        "production": {str(w): asis(100.0 * float(occ[p, k].mean()), 1) for k, w in enumerate(WATERS)},
        "production_20_to_100_ns": {str(w): asis(100.0 * float(occ[late, k].mean()), 1) for k, w in enumerate(WATERS)},
        "n_frames_production": int(len(p)), "n_frames_20_to_100_ns": int(len(late))}
    don = A["donor"]
    rec["water_donating_to_N2"] = {
        "what": "a water with O within 3.5 A of N2 and O-H...N2 at least 135 deg (no crystal position is used); per cent of the "
                "all-atom frames of the production",
        "share_production": asis(100.0 * float(don[p].mean()), 1),
        "mean_O_N2_of_the_donor_A": asis(float(np.nanmean(A["d_donor"][p])), 2) if don[p].any() else None,
        "share_while_place_553_occupied": asis(100.0 * float(don[p][occ[p, WATERS.index(553)]].mean()), 1) if occ[p, WATERS.index(553)].any() else None}
    rec["water_within_8_A_of_the_zinc"] = {"mean_production_20_to_100_ns": asis(float(A["n8"][late].mean()), 2),
                                            "mean_first_10_ns": asis(float(A["n8"][p[:100]].mean()), 2)}
    # by the state of the solute frame of the same time, from the solute series of the first program (state only)
    st = [r["state"] for r in rows_of(os.path.join(d, "analysis", "timeseries_solute.csv"))]
    stf = np.array([st[10 * (k + 1) - 1] for k in range(len(p))])
    rec["by_state_of_the_solute_frame"] = {
        s: {"n_all_atom_frames": int((stf == s).sum()),
            "place_553_occupied_share": asis(100.0 * float(occ[p, WATERS.index(553)][stf == s].mean()), 1) if (stf == s).any() else None,
            "water_in_the_place_of_the_crystal_ring_share": asis(100.0 * float((A["n_ring"][p][stf == s] >= 1).mean()), 1) if (stf == s).any() else None}
        for s in "CRI"}
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--part", choices=("solute", "full"), required=True)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--runs", nargs="*", default=ALL)
    args = ap.parse_args(argv)
    idx_of = tleap_index()
    fn = solute_run if args.part == "solute" else full_run
    recs = []
    for lab in args.runs:
        recs.append(fn(lab, idx_of))
        print("%s %s done" % (args.part, lab), flush=True)
    out = {"what": "Checks on the raw frames for md/RESULTS.md (%s frames). Written by md/production/site_checks.py; own "
                   "reading of the DCD files, set beside the series of the two analysis programs." % ("solute" if args.part == "solute" else "all-atom"),
           "script": {"name": "md/production/site_checks.py", "sha256": sha256_file(os.path.abspath(__file__))},
           "written": datetime.datetime.now(IST).strftime("%-d %b %Y, %-I:%M %p IST"),
           "tolerances": {"distance_A": TOL_D, "angle_deg": TOL_A, "signed_volume": TOL_V},
           "runs": recs}
    path = os.path.join(args.out, "site_checks_%s.json" % args.part)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print("written to %s" % path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
