#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""analyse_md.py -- primary implementation of md/analysis/SPEC.md (version 1.3, fixed 6 Oct 2026).

Analysis of the molecular-dynamics runs of the hCA II / sulfanilamide complex (md/PLAN.md, stage S6).
Every definition, threshold, file name and JSON key used here is fixed in md/analysis/SPEC.md, which was
written before any production run existed; this script implements it and must not be changed to fit results.

Usage (python of the MD environment; MDAnalysis, NumPy, Matplotlib):

  analyse_md.py RUN_DIR [--system DIR] [--start X|N|R] [--planned-ns 100] [--out DIR] [--no-figures]
      one run directory (solute_segNNNN.dcd, full_segNNNN.dcd, solute.pdb, full.pdb, manifest.json);
      works on a run that is not finished and says how many ns it saw.  Run X0 of a zinc model is the
      validation run: planned length 5 ns, gate G-V with a verdict, never part of the campaign table.
      The equilibration frames (solute_eqNN.dcd, full_eqNN.dcd), if present, are read into a separate time
      series with negative times; nothing of them enters a statistic, a verdict or a gate.
  analyse_md.py --campaign DIR [--out DIR] [--reuse]
      every run directory below DIR, then the campaign table (CSV, Markdown), JSON and figures.
  analyse_md.py --check-starts [--out DIR]
      the definitions applied to the three starting poses of md/system/, against md/system/starts.csv.

Units: coordinates and distances in Angstrom (as stored in DCD and PDB files), angles in degrees, times in ns
(intervals in ps).  No random numbers are used; a second run on the same input gives the same files.
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')          # one thread: the machine is shared, and nothing here needs more

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np

SPEC_VERSION = '1.3 (6 Oct 2026)'

# ----------------------------------------------------------------------------------------------------------
# Constants of the SPEC (section numbers refer to md/analysis/SPEC.md)
# ----------------------------------------------------------------------------------------------------------
LIG_HEAVY = ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S')            # section 3
_RING = {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3'}
_OXY = {'O1': 'O2', 'O2': 'O1'}
LABELLINGS = (('identity', {}), ('ring_flip', _RING), ('oxygen_swap', _OXY),
              ('ring_flip+oxygen_swap', {**_RING, **_OXY}))

CONST = {
    'solute_interval_ps_default': 10.0,
    'full_interval_ps_default': 100.0,
    'planned_ns_default': 100.0,
    'midpoint_d1_N1_F131CZ_A': 6.7,          # section 4
    'midpoint_d2_N1_N67OD1_A': 5.2,
    'verdict_window_ns': 20.0,
    'verdict_fraction': [9, 10],             # at least 9/10 of the frames of the window
    'dwell_ns': 1.0,
    'hbond_distance_max_A': 3.5,             # donor ... acceptor, section 5
    'hbond_angle_min_deg': 135.0,            # donor - H ... acceptor, at the hydrogen
    'fifth_ligand_distance_max_A': 2.6,      # section 6
    'water_near_N1_distance_max_A': 3.5,     # section 7
    'rule1_rmsd_A': 2.0,                     # plan, section 5, rule 1
    'block_ns': 10.0,                        # section 10
    'validation_run': {'start': 'X', 'replica': 0, 'planned_ns': 5.0},      # run X0 of a zinc model, section 1
    'gate_GV': {                             # plan, section 4, gates G-V and G-Z (section 11)
        'lig_rmsd_mean_max_A': 1.0,
        'hb_N2H_T199OG1_fraction_min': [4, 5],
        'd_Zn_N2_mean_window_A': [1.83, 2.03],
        'd_Zn_NHis_mean_window_A': [1.91, 2.12],
        'ang_N_Zn_N_reference_6RL9_deg': {   # the six angles of the crystal start (6RL9), one decimal
            'ang_N2_Zn_H94NE2_deg': 108.9, 'ang_N2_Zn_H96NE2_deg': 114.3, 'ang_N2_Zn_H119ND1_deg': 116.8,
            'ang_H94NE2_Zn_H96NE2_deg': 103.8, 'ang_H94NE2_Zn_H119ND1_deg': 113.6,
            'ang_H96NE2_Zn_H119ND1_deg': 98.4},
        'ang_N_Zn_N_tolerance_deg': 10.0,
        'fifth_ligand_fraction_max': [1, 20],        # at most 5 % of the frames, in each of the two trajectories
        'ang_C4_S_N2_reference_deg': 108.8,
        'ang_C4_S_N2_tolerance_deg': 3.0,
        'tau_fold_reference_deg': 132.1,
        'tau_fold_tolerance_deg': 20.0,
    },
    'state_of_start': {'X': 'C', 'N': 'C', 'R': 'R'},
    # ---- version 1.2: reported quantities only (sections 9A to 9E); no verdict and no gate uses them
    'n2_sign_of_start': {'X': 1, 'N': 1, 'R': -1},        # section 9A: sign of the signed volume at N2 in the built starts
    'n2_dwell_ns': 1.0,                                   # a sign change counts if the new sign lasts this long
    'glu106_distance_max_A': 2.6,                         # section 9C: the limit of the fifth ligand
    # ---- version 1.3, section 9D: water in the place of the crystal ring.  The count the driver acts on in its
    # hydration check: site atoms and cut-off are those of the run (manifest hydration_check.settings) where the
    # manifest gives them; otherwise these two, the values of md/protocol.json after md/PLAN.md, 6 Oct 2026,
    # 8:01 AM point 7 and 8:56 AM point 6
    'ring_place_atoms_default': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
    'ring_place_water_cutoff_default_A': 2.5,
    # the wider count of version 1.2 (then named ..._ring_site_...): seven positions, 3.5 A; constants
    'ring_wide_atoms': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'],
    'ring_wide_water_distance_max_A': 3.5,
    'equilibration_stages_protocol_ps': {                 # section 9E: start and end of the six stages of md/protocol.json
        '01': [0.0, 100.0], '02': [100.0, 400.0], '03': [400.0, 700.0], '04': [700.0, 1100.0],
        '05': [1100.0, 1500.0], '06': [1500.0, 2000.0]},
}

# selections.json key -> atom name expected in solute.pdb / full.pdb (guard against a wrong mapping)
SITE = {'zn': 'ZN', 'his94_NE2': 'NE2', 'his96_NE2': 'NE2', 'his119_ND1': 'ND1', 'thr199_OG1': 'OG1',
        'thr199_HG1': 'HG1', 'thr199_N': 'N', 'thr199_H': 'H', 'glu106_OE1': 'OE1', 'glu106_OE2': 'OE2',
        'phe131_CZ': 'CZ', 'asn67_OD1': 'OD1'}
# optional key of selections.json -> atom name; read if the key is there and not null (section 4)
SITE_OPTIONAL = {'asn67_ND2': 'ND2'}

# column, kind.  kinds: A distance (Angstrom), deg angle, tors signed torsion (circular statistics),
# flag 0/1, int integer, state C/R/I
SOLUTE_COLUMNS = (
    ('lig_rmsd_A', 'A'), ('lig_rmsd_labelling', 'int'), ('lig_rmsd_identity_A', 'A'),
    ('d1_N1_F131CZ_A', 'A'), ('d2_N1_N67OD1_A', 'A'), ('state', 'state'), ('state_assigned', 'state'),
    ('d_N2_T199OG1_A', 'A'), ('ang_N2_HN2_T199OG1_deg', 'deg'), ('hb_N2H_T199OG1', 'flag'),
    ('d_Zn_N2_A', 'A'), ('d_Zn_H94NE2_A', 'A'), ('d_Zn_H96NE2_A', 'A'), ('d_Zn_H119ND1_A', 'A'),
    ('ang_N2_Zn_H94NE2_deg', 'deg'), ('ang_N2_Zn_H96NE2_deg', 'deg'), ('ang_N2_Zn_H119ND1_deg', 'deg'),
    ('ang_H94NE2_Zn_H96NE2_deg', 'deg'), ('ang_H94NE2_Zn_H119ND1_deg', 'deg'),
    ('ang_H96NE2_Zn_H119ND1_deg', 'deg'), ('ang_S_N2_Zn_deg', 'deg'),
    ('d_Zn_O1_A', 'A'), ('d_Zn_O2_A', 'A'), ('d_Zn_Osulfonyl_min_A', 'A'), ('fifth_sulfonylO', 'flag'),
    ('d_T199N_O1_A', 'A'), ('d_T199N_O2_A', 'A'), ('ang_T199N_H_O1_deg', 'deg'), ('ang_T199N_H_O2_deg', 'deg'),
    ('d_T199N_Osulfonyl_min_A', 'A'), ('hb_T199NH_Osulfonyl', 'flag'),
    ('d_T199OG1_E106OE1_A', 'A'), ('d_T199OG1_E106OE2_A', 'A'), ('ang_T199OG1_HG1_E106OE1_deg', 'deg'),
    ('ang_T199OG1_HG1_E106OE2_deg', 'deg'), ('d_T199OG1_E106OE_min_A', 'A'), ('hb_T199OG1H_E106OE', 'flag'),
    ('tors_C3_C4_S_N2_deg', 'tors'), ('tors_C5_C4_S_N2_deg', 'tors'), ('tau_fold_deg', 'deg'),
    ('ang_C4_S_N2_deg', 'deg'), ('d_S_N2_A', 'A'), ('ca_rmsd_A', 'A'),
)
FULL_COLUMNS = (
    ('d_Zn_Owater_min_A', 'A'), ('water_O_index_nearest_Zn', 'int'), ('d_Zn_Oprotein_min_A', 'A'),
    ('protein_O_index_nearest_Zn', 'int'), ('d_Zn_Osulfonyl_min_full_A', 'A'), ('d_Zn_O_any_min_A', 'A'),
    ('fifth_water', 'flag'), ('fifth_proteinO', 'flag'), ('fifth_any', 'flag'), ('n_water_N1_3p5A', 'int'),
    ('d_Zn_N2_full_A', 'A'),
)
# version 1.2: columns that follow ca_rmsd_A in timeseries_solute.csv (sections 9A and 9C), and the two columns that
# follow d_Zn_N2_full_A in timeseries_full.csv (section 9D; named ..._ring_site_... in version 1.2, ..._ring_wide_...
# since version 1.3).  vol: the signed volume at N2, no unit
V12_SOLUTE_COLUMNS = (
    ('n2_volume', 'vol'), ('n2_sign', 'int'), ('n2_sign_assigned', 'int'),
    ('d_Zn_E106OE1_A', 'A'), ('d_Zn_E106OE2_A', 'A'), ('d_Zn_E106OE_min_A', 'A'),
    ('fifth_E106OE1', 'flag'), ('fifth_E106OE2', 'flag'), ('fifth_E106OE', 'flag'),
)
V12_FULL_COLUMNS = (('n_water_ring_wide_3p5A', 'int'), ('d_water_ring_wide_min_A', 'A'))
# version 1.3: the two columns that follow them (section 9D): the count the driver acts on, and its nearest distance
V13_FULL_COLUMNS = (('n_water_ring_place', 'int'), ('d_water_ring_place_min_A', 'A'))
ALL_FULL_COLUMNS = FULL_COLUMNS + V12_FULL_COLUMNS + V13_FULL_COLUMNS
FULL_COUNT_COLUMNS = ('n_water_N1_3p5A', 'n_water_ring_wide_3p5A', 'n_water_ring_place')    # integer columns that get statistics
# columns of the production series that depend on the frames before them; not written for the equilibration (section 9E)
HISTORY_COLUMNS = ('state_assigned', 'n2_sign_assigned')
# columns appended to timeseries_solute.csv when selections.json has the optional key asn67_ND2 (section 4)
ASN67_COLUMNS = (('d_N1_N67ND2_A', 'A'), ('asn67_ND2_nearer_N1', 'flag'))
ZN_ANGLES = ('ang_N2_Zn_H94NE2_deg', 'ang_N2_Zn_H96NE2_deg', 'ang_N2_Zn_H119ND1_deg',
             'ang_H94NE2_Zn_H96NE2_deg', 'ang_H94NE2_Zn_H119ND1_deg', 'ang_H96NE2_Zn_H119ND1_deg')
ZN_HIS = ('d_Zn_H94NE2_A', 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A')
FMT = {'A': '%.4f', 'deg': '%.3f', 'tors': '%.3f', 'flag': '%d', 'int': '%d', 'state': '%s', 'vol': '%.4f'}


class AnalysisError(Exception):
    pass


# ----------------------------------------------------------------------------------------------------------
# Geometry (vectorised over frames; inputs [..., 3] in Angstrom)
# ----------------------------------------------------------------------------------------------------------
def dist(a, b):
    return np.linalg.norm(np.asarray(a, float) - np.asarray(b, float), axis=-1)


def angle(a, b, c):
    """Angle a-b-c at b, degrees, as atan2(|u x v|, u . v) with u = a - b, v = c - b."""
    u, v = np.asarray(a, float) - np.asarray(b, float), np.asarray(c, float) - np.asarray(b, float)
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v), axis=-1), (u * v).sum(-1)))


def dihedral(a, b, c, d):
    """Torsion a-b-c-d in (-180, 180], IUPAC sign (clockwise looking from b to c is positive)."""
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    b0, b1, b2 = a - b, c - b, d - c
    e = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * e).sum(-1, keepdims=True) * e
    w = b2 - (b2 * e).sum(-1, keepdims=True) * e
    return np.degrees(np.arctan2((np.cross(e, v) * w).sum(-1), (v * w).sum(-1)))


def unit_vector(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def kabsch(X, R):
    """Least-squares superposition (no weights) of X [n, m, 3] onto R [m, 3].
    Returns Q [n, 3, 3], cx [n, 1, 3], cr [3]; a point y of frame n maps to (y - cx[n]) @ Q[n] + cr."""
    X, R = np.asarray(X, float), np.asarray(R, float)
    cx, cr = X.mean(axis=1, keepdims=True), R.mean(axis=0)
    H = np.einsum('nmi,mj->nij', X - cx, R - cr)
    U, _, Vt = np.linalg.svd(H)
    D = np.zeros_like(H)
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))      # proper rotation only
    return U @ D @ Vt, cx, cr


def transform(Y, Q, cx, cr):
    Y = np.asarray(Y, float)
    if Y.ndim == 2:                                  # [n, 3]
        return np.einsum('ni,nij->nj', Y - cx[:, 0, :], Q) + cr
    return np.einsum('nki,nij->nkj', Y - cx, Q) + cr


def box_matrix(dims):
    """Lattice vectors (rows, Angstrom) from a, b, c, alpha, beta, gamma (degrees); a along x, b in xy."""
    a, b, c, al, be, ga = (float(x) for x in dims)
    al, be, ga = np.radians([al, be, ga])
    bx, by = b * np.cos(ga), b * np.sin(ga)
    cx_ = c * np.cos(be)
    cy_ = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz_ = math.sqrt(max(c * c - cx_ * cx_ - cy_ * cy_, 0.0))
    return np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx_, cy_, cz_]])


def min_image_shift(v, B):
    """Lattice vector s (a combination of the rows of B with integer coefficients) that makes |v + s| smallest.
    Search: coefficients round(-f) + {-1, 0, 1} in each direction, f the fractional coordinates of v."""
    f = np.linalg.solve(B.T, np.asarray(v, float))
    k0 = -np.rint(f)
    best, bestd = None, None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                s = (k0 + (i, j, k)) @ B
                d = float(np.linalg.norm(v + s))
                if bestd is None or d < bestd - 1e-9:
                    best, bestd = s, d
    return best


# ----------------------------------------------------------------------------------------------------------
# Small I/O helpers
# ----------------------------------------------------------------------------------------------------------
def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_pdb(path):
    """Atoms of the first model of a PDB file, by fixed columns: name, residue name, residue number, xyz."""
    names, resn, resi, xyz = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            if line.startswith(('ATOM  ', 'HETATM')):
                names.append(line[12:16].strip())
                resn.append(line[17:21].strip())
                resi.append(line[22:27].strip())
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return names, resn, resi, np.array(xyz, float).reshape(-1, 3)


def r6(v):
    if v is None:
        return None
    v = float(v)
    return round(v, 6) if math.isfinite(v) else None


def segment_files(run_dir, kind):
    out = []
    for p in sorted(Path(run_dir).glob(f'{kind}_seg*.dcd')):
        m = re.fullmatch(rf'{kind}_seg(\d+)\.dcd', p.name)
        if m:
            out.append((int(m.group(1)), p))
    return sorted(out)


EQ_STEP = r'eq(\d+)(?:x([1-9]\d*))?'           # eqNN: stage NN; eqNNxM: extension step M of stage NN (section 9E)


def step_code(nn, mm=0):
    """One whole number per equilibration step, in the order of the run: stage NN, then its extension steps 1, 2, ..."""
    return 1000 * int(nn) + int(mm or 0)


def step_tag(code):
    nn, mm = divmod(int(code), 1000)
    return f'eq{nn:02d}' + (f'x{mm}' if mm else '')


def equilibration_files(run_dir, kind):
    """Files <kind>_eqNN.dcd of the equilibration stages and <kind>_eqNNxM.dcd of the extension steps of a stage
    (section 9E), as (step code, path), in the order of the run."""
    out = []
    for p in sorted(Path(run_dir).glob(f'{kind}_eq*.dcd')):
        m = re.fullmatch(rf'{kind}_{EQ_STEP}\.dcd', p.name)
        if m:
            out.append((step_code(m.group(1), m.group(2)), p))
    return sorted(out)


def software_versions():
    import MDAnalysis
    import matplotlib
    return {'python': sys.version.split()[0], 'numpy': np.__version__, 'MDAnalysis': MDAnalysis.__version__,
            'matplotlib': matplotlib.__version__}


# ----------------------------------------------------------------------------------------------------------
# Quantities of the solute trajectory (SPEC sections 2 to 9)
# ----------------------------------------------------------------------------------------------------------
def solute_quantities(A, ca, ref_ca, ref_lig):
    """A: dict of [n, 3] arrays in the frame of the trajectory after imaging: the keys of SITE, the eleven
    names of LIG_HEAVY and 'HN2'; 'asn67_ND2' if selections.json has it.  ca [n, m, 3]; ref_ca [m, 3];
    ref_lig {name: [3]}.  Returns (columns {name: [n]}, fitted CA [n, m, 3])."""
    Q, cx, cr = kabsch(ca, ref_ca)                                                   # section 2
    caf = transform(ca, Q, cx, cr)
    q = {}
    q['ca_rmsd_A'] = np.sqrt(((caf - ref_ca) ** 2).sum(-1).mean(-1))                 # section 9

    lig = np.stack([transform(A[n], Q, cx, cr) for n in LIG_HEAVY], axis=1)          # [n, 11, 3], fitted frame
    rms = []
    for _, perm in LABELLINGS:                                                       # section 3
        ref = np.array([ref_lig[perm.get(n, n)] for n in LIG_HEAVY], float)
        rms.append(np.sqrt(((lig - ref) ** 2).sum(-1).mean(-1)))
    rms = np.array(rms)                                                              # [4, n]
    lab = np.argmin(rms, axis=0)                                                     # ties: lowest number
    q['lig_rmsd_A'] = rms[lab, np.arange(rms.shape[1])]
    q['lig_rmsd_labelling'] = lab
    q['lig_rmsd_identity_A'] = rms[0]

    d1 = dist(A['N1'], A['phe131_CZ'])                                               # section 4
    d2 = dist(A['N1'], A['asn67_OD1'])
    m1, m2 = CONST['midpoint_d1_N1_F131CZ_A'], CONST['midpoint_d2_N1_N67OD1_A']
    state = np.full(d1.shape, 'I', dtype='U1')
    state[(d1 < m1) & (d2 > m2)] = 'C'
    state[(d1 > m1) & (d2 < m2)] = 'R'
    q['d1_N1_F131CZ_A'], q['d2_N1_N67OD1_A'], q['state'] = d1, d2, state
    if 'asn67_ND2' in A:                     # optional: reported beside d2, never used for the state
        q['d_N1_N67ND2_A'] = dist(A['N1'], A['asn67_ND2'])
        q['asn67_ND2_nearer_N1'] = (q['d_N1_N67ND2_A'] < d2).astype(int)

    hd, ha = CONST['hbond_distance_max_A'], CONST['hbond_angle_min_deg']             # section 5
    q['d_N2_T199OG1_A'] = dist(A['N2'], A['thr199_OG1'])
    q['ang_N2_HN2_T199OG1_deg'] = angle(A['N2'], A['HN2'], A['thr199_OG1'])
    q['hb_N2H_T199OG1'] = ((q['d_N2_T199OG1_A'] <= hd) & (q['ang_N2_HN2_T199OG1_deg'] >= ha)).astype(int)

    zn = A['zn']                                                                     # section 6
    lig4 = (('N2', A['N2']), ('H94NE2', A['his94_NE2']), ('H96NE2', A['his96_NE2']), ('H119ND1', A['his119_ND1']))
    for name, x in lig4:
        q[f'd_Zn_{name}_A'] = dist(zn, x)
    for i in range(4):
        for j in range(i + 1, 4):
            q[f'ang_{lig4[i][0]}_Zn_{lig4[j][0]}_deg'] = angle(lig4[i][1], zn, lig4[j][1])
    q['ang_S_N2_Zn_deg'] = angle(A['S'], A['N2'], zn)
    q['d_Zn_O1_A'], q['d_Zn_O2_A'] = dist(zn, A['O1']), dist(zn, A['O2'])
    q['d_Zn_Osulfonyl_min_A'] = np.minimum(q['d_Zn_O1_A'], q['d_Zn_O2_A'])
    q['fifth_sulfonylO'] = (q['d_Zn_Osulfonyl_min_A'] <= CONST['fifth_ligand_distance_max_A']).astype(int)

    tn, th = A['thr199_N'], A['thr199_H']                                            # section 7
    q['d_T199N_O1_A'], q['d_T199N_O2_A'] = dist(tn, A['O1']), dist(tn, A['O2'])
    q['ang_T199N_H_O1_deg'], q['ang_T199N_H_O2_deg'] = angle(tn, th, A['O1']), angle(tn, th, A['O2'])
    q['d_T199N_Osulfonyl_min_A'] = np.minimum(q['d_T199N_O1_A'], q['d_T199N_O2_A'])
    q['hb_T199NH_Osulfonyl'] = (((q['d_T199N_O1_A'] <= hd) & (q['ang_T199N_H_O1_deg'] >= ha)) |
                                ((q['d_T199N_O2_A'] <= hd) & (q['ang_T199N_H_O2_deg'] >= ha))).astype(int)
    og, hg = A['thr199_OG1'], A['thr199_HG1']
    q['d_T199OG1_E106OE1_A'], q['d_T199OG1_E106OE2_A'] = dist(og, A['glu106_OE1']), dist(og, A['glu106_OE2'])
    q['ang_T199OG1_HG1_E106OE1_deg'] = angle(og, hg, A['glu106_OE1'])
    q['ang_T199OG1_HG1_E106OE2_deg'] = angle(og, hg, A['glu106_OE2'])
    q['d_T199OG1_E106OE_min_A'] = np.minimum(q['d_T199OG1_E106OE1_A'], q['d_T199OG1_E106OE2_A'])
    q['hb_T199OG1H_E106OE'] = (((q['d_T199OG1_E106OE1_A'] <= hd) & (q['ang_T199OG1_HG1_E106OE1_deg'] >= ha)) |
                               ((q['d_T199OG1_E106OE2_A'] <= hd) & (q['ang_T199OG1_HG1_E106OE2_deg'] >= ha))
                               ).astype(int)

    q['tors_C3_C4_S_N2_deg'] = dihedral(A['C3'], A['C4'], A['S'], A['N2'])           # section 8
    q['tors_C5_C4_S_N2_deg'] = dihedral(A['C5'], A['C4'], A['S'], A['N2'])
    q['tau_fold_deg'] = np.maximum(np.abs(q['tors_C3_C4_S_N2_deg']), np.abs(q['tors_C5_C4_S_N2_deg']))
    q['ang_C4_S_N2_deg'] = angle(A['C4'], A['S'], A['N2'])
    q['d_S_N2_A'] = dist(A['S'], A['N2'])

    # section 9A: signed volume at N2 = e(Zn) . (e(HN21) x e(S)), unit vectors from N2, in this order
    e_zn, e_h, e_s = (unit_vector(A[k] - A['N2']) for k in ('zn', 'HN2', 'S'))
    q['n2_volume'] = (e_zn * np.cross(e_h, e_s)).sum(-1)
    q['n2_sign'] = np.sign(q['n2_volume']).astype(int)
    # section 9C: the zinc to the two carboxylate oxygens of Glu106, by name
    cg = CONST['glu106_distance_max_A']
    q['d_Zn_E106OE1_A'], q['d_Zn_E106OE2_A'] = dist(zn, A['glu106_OE1']), dist(zn, A['glu106_OE2'])
    q['d_Zn_E106OE_min_A'] = np.minimum(q['d_Zn_E106OE1_A'], q['d_Zn_E106OE2_A'])
    q['fifth_E106OE1'] = (q['d_Zn_E106OE1_A'] <= cg).astype(int)
    q['fifth_E106OE2'] = (q['d_Zn_E106OE2_A'] <= cg).astype(int)
    q['fifth_E106OE'] = (q['d_Zn_E106OE_min_A'] <= cg).astype(int)
    return q, caf


# ----------------------------------------------------------------------------------------------------------
# Orientation: states, dwells, verdict (SPEC section 4)
# ----------------------------------------------------------------------------------------------------------
def assign_states(strict, start_state):
    """Frames in I take the last definite state (C or R) before them; before the first definite frame they
    take the state of the start."""
    out, cur = [], start_state
    for s in strict:
        if s in ('C', 'R'):
            cur = s
        out.append(cur)
    return np.array(out, dtype='U1')


def orientation(strict, start_state, dt_ns, n_planned_frames):
    n = len(strict)
    opp = 'R' if start_state == 'C' else 'C'
    assigned = assign_states(strict, start_state)
    dwell = max(1, int(round(CONST['dwell_ns'] / dt_ns)))
    visits, i = [], 0
    while i < n:                                             # maximal stretches of one assigned state
        j = i
        while j + 1 < n and assigned[j + 1] == assigned[i]:
            j += 1
        seg = strict[i:j + 1]
        visits.append({'state': str(assigned[i]), 'first_frame': i, 'n_frames': j - i + 1,
                       't_first_ns': r6((i + 1) * dt_ns), 't_last_ns': r6((j + 1) * dt_ns),
                       'long': (j - i + 1) >= dwell,
                       'fraction_strict': r6(float(np.mean(seg == assigned[i])))})
        i = j + 1
    longv = [v for v in visits if v['long']]
    seq = [start_state] + [v['state'] for v in longv]
    n_trans = sum(1 for a, b in zip(seq[:-1], seq[1:]) if a != b)
    first = next((v for v in longv if v['state'] == opp), None)

    def counts(sl):
        s = strict[sl]
        return {k: int(np.sum(s == k)) for k in ('C', 'R', 'I')}

    w = min(n, int(round(CONST['verdict_window_ns'] / dt_ns)))
    cw = counts(slice(n - w, n))
    num, den = CONST['verdict_fraction']
    if den * cw[opp] >= num * w:
        verdict = 'CONVERTED'
    elif den * cw[start_state] >= num * w:
        verdict = 'STAYED'
    else:
        verdict = 'MIXED'
    h = n // 2
    call, clast = counts(slice(0, n)), counts(slice(n - h, n))
    return assigned, {
        'start_state': start_state, 'opposite_state': opp,
        'midpoint_d1_A': CONST['midpoint_d1_N1_F131CZ_A'], 'midpoint_d2_A': CONST['midpoint_d2_N1_N67OD1_A'],
        'window_frames': w, 'window_ns': r6(w * dt_ns),
        'window_complete': w == int(round(CONST['verdict_window_ns'] / dt_ns)),
        'count_C_window': cw['C'], 'count_R_window': cw['R'], 'count_I_window': cw['I'],
        'fraction_start_state_window': r6(cw[start_state] / w), 'fraction_opposite_state_window': r6(cw[opp] / w),
        'fraction_intermediate_window': r6(cw['I'] / w),
        'verdict': verdict, 'provisional': n < n_planned_frames,
        'dwell_frames': dwell,
        'first_passage_ns': None if first is None else first['t_first_ns'],
        'first_passage_frame': None if first is None else first['first_frame'],
        'n_transitions': n_trans,
        'n_visits': len(visits), 'n_long_visits': len(longv), 'long_visits': longv,
        'state_first_frame': str(strict[0]),
        'count_C': call['C'], 'count_R': call['R'], 'count_I': call['I'],
        'fraction_C': r6(call['C'] / n), 'fraction_R': r6(call['R'] / n), 'fraction_I': r6(call['I'] / n),
        'n_last_half': h,
        'count_C_last_half': clast['C'], 'count_R_last_half': clast['R'], 'count_I_last_half': clast['I'],
    }


# ----------------------------------------------------------------------------------------------------------
# Configuration at N2 and the hydrogen bond around the first passage (SPEC sections 9A and 9B; reported only)
# ----------------------------------------------------------------------------------------------------------
def n2_configuration(volume, strict, start_sign, dt_ns):
    """Section 9A.  volume [n]: the signed volume at N2 of every frame; strict [n]: the column `state`.
    Returns (assigned sign [n], dictionary for summary.json)."""
    sign = np.sign(np.asarray(volume, float)).astype(int)
    n = len(sign)
    assigned, cur = np.zeros(n, int), int(start_sign)
    for i in range(n):                                       # a frame with volume zero keeps the sign before it
        if sign[i] != 0:
            cur = int(sign[i])
        assigned[i] = cur
    dwell = max(1, int(round(CONST['n2_dwell_ns'] / dt_ns)))
    stretches, a = [], 0
    while a < n:                                             # maximal stretches of one assigned sign
        b = a
        while b + 1 < n and assigned[b + 1] == assigned[a]:
            b += 1
        stretches.append({'sign': int(assigned[a]), 'first_frame': a, 'n_frames': b - a + 1,
                          't_first_ns': r6((a + 1) * dt_ns), 't_last_ns': r6((b + 1) * dt_ns),
                          'long': (b - a + 1) >= dwell})
        a = b + 1
    longs = [x for x in stretches if x['long']]
    changes, prev = [], int(start_sign)
    for x in longs:                                          # the sign of the start, then the long stretches in order
        if x['sign'] != prev:
            changes.append({'frame': x['first_frame'], 'time_ns': x['t_first_ns'], 'from_sign': prev, 'to_sign': x['sign']})
        prev = x['sign']

    def counts(v):
        return {'n': int(len(v)), 'count_positive': int(np.sum(v > 0)), 'count_negative': int(np.sum(v < 0)),
                'count_zero': int(np.sum(v == 0))}
    w = min(n, int(round(CONST['verdict_window_ns'] / dt_ns)))
    ca, cw = counts(sign), counts(sign[n - w:])
    strict = np.asarray(strict)

    def by_state(k):                                         # frames whose `state` (not the assigned state) is k
        c = counts(sign[strict == k])
        c['fraction_positive'] = r6(c['count_positive'] / c['n']) if c['n'] else None
        c['fraction_negative'] = r6(c['count_negative'] / c['n']) if c['n'] else None
        return c
    return assigned, {
        'order': 'e(Zn) . (e(HN21) x e(S)), unit vectors from N2',
        'sign_of_start': int(start_sign),
        'n': ca['n'], 'count_positive': ca['count_positive'], 'count_negative': ca['count_negative'],
        'count_zero': ca['count_zero'],
        'fraction_positive': r6(ca['count_positive'] / n), 'fraction_negative': r6(ca['count_negative'] / n),
        'window_frames': w, 'count_positive_window': cw['count_positive'], 'count_negative_window': cw['count_negative'],
        'count_zero_window': cw['count_zero'],
        'fraction_positive_window': r6(cw['count_positive'] / w), 'fraction_negative_window': r6(cw['count_negative'] / w),
        'by_state': {k: by_state(k) for k in ('C', 'R', 'I')},
        'sign_first_frame': int(sign[0]), 'sign_last_frame': int(sign[-1]),
        'dwell_frames': dwell, 'n_stretches': len(stretches), 'n_long_stretches': len(longs), 'long_stretches': longs,
        'n_lasting_sign_changes': len(changes), 'lasting_sign_changes': changes,
        'first_lasting_sign_change_ns': changes[0]['time_ns'] if changes else None,
    }


def hbond_first_passage(hb, orient):
    """Section 9B: the hydrogen bond N2-H...OG1(Thr199) before the first passage (frames 0 .. p - 1) and from it
    on (frames p .. N - 1); p is the first frame of the first long visit in the opposite state."""
    p = orient['first_passage_frame']
    out = {'available': p is not None, 'flag': 'hb_N2H_T199OG1', 'verdict': orient['verdict'],
           'first_passage_frame': p, 'first_passage_ns': orient['first_passage_ns'],
           'n_before': None, 'count_before': None, 'fraction_before': None,
           'n_after': None, 'count_after': None, 'fraction_after': None}
    if p is None:
        return out
    hb = np.asarray(hb, int)
    before, after = hb[:p], hb[p:]
    out.update(n_before=int(len(before)), count_before=int(before.sum()),
               fraction_before=r6(before.mean()) if len(before) else None,
               n_after=int(len(after)), count_after=int(after.sum()), fraction_after=r6(after.mean()))
    return out


# ----------------------------------------------------------------------------------------------------------
# Statistics (SPEC section 10)
# ----------------------------------------------------------------------------------------------------------
def _circ(x):
    a = np.radians(np.asarray(x, float))
    s, c = np.sin(a).mean(), np.cos(a).mean()
    R = min(1.0, math.hypot(s, c))
    return math.degrees(math.atan2(s, c)), (math.degrees(math.sqrt(-2.0 * math.log(R))) if R > 0 else None)


def summarise(x, kind, block):
    x = np.asarray(x, float)
    n = len(x)
    out = {'kind': kind, 'n': int(n), 'mean': None, 'sd': None, 'min': None, 'max': None, 'n_last_half': int(n // 2),
           'mean_last_half': None, 'sd_last_half': None, 'block_frames': int(block), 'block_means': [],
           'block_sem': None}
    if n == 0:
        return out
    last = x[n - n // 2:]
    blocks = [x[k * block:(k + 1) * block] for k in range(n // block)] if block > 0 else []
    if kind == 'tors':
        out['mean'], out['sd'] = (r6(v) for v in _circ(x))
        if len(last):
            out['mean_last_half'], out['sd_last_half'] = (r6(v) for v in _circ(last))
        out['block_means'] = [r6(_circ(b)[0]) for b in blocks]
    else:
        out['mean'], out['sd'] = r6(x.mean()), r6(x.std())           # population standard deviation
        if len(last):
            out['mean_last_half'], out['sd_last_half'] = r6(last.mean()), r6(last.std())
        bm = [float(b.mean()) for b in blocks]
        out['block_means'] = [r6(v) for v in bm]
        if len(bm) >= 2:
            out['block_sem'] = r6(np.std(bm, ddof=1) / math.sqrt(len(bm)))
    out['min'], out['max'] = r6(x.min()), r6(x.max())
    return out


def fraction(flag, block):
    f = np.asarray(flag, int)
    n = len(f)
    h = n // 2
    out = {'n': int(n), 'count': int(f.sum()), 'fraction': r6(f.mean()) if n else None,
           'n_last_half': int(h), 'count_last_half': int(f[n - h:].sum()),
           'fraction_last_half': r6(f[n - h:].mean()) if h else None, 'block_frames': int(block),
           'block_fractions': [r6(f[k * block:(k + 1) * block].mean()) for k in range(n // block)] if block > 0 else []}
    return out


# ----------------------------------------------------------------------------------------------------------
# One run
# ----------------------------------------------------------------------------------------------------------
def _dig(d, *path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _first(*pairs):
    """First value that is not None among (value, source) pairs; returns (value, source)."""
    for v, src in pairs:
        if v is not None:
            return v, src
    return None, None


def identify_run(run_dir, args):
    """Who the run is and where its system is (SPEC section 1).  Order: command line, manifest.json (the plain
    keys, then the keys of the manifest of md/run_md.py), folder name, default.  Start X with replica 0 (run X0)
    is the validation run: its planned length is 5 ns unless the command line gives another."""
    run_dir = Path(run_dir).resolve()
    man = {}
    mp = run_dir / 'manifest.json'
    if mp.is_file():
        try:
            man = json.loads(mp.read_text())
        except Exception as e:                                # an unreadable manifest is reported, not fatal
            man = {'_error': f'{type(e).__name__}: {e}'}
    pat = r'([XNR])(\d+)'
    mf = re.fullmatch(pat, run_dir.name)
    mi = re.fullmatch(pat, man['run_id']) if isinstance(man.get('run_id'), str) else None
    bs = _dig(man, 'build_json', 'start')
    src = {}
    start, src['start'] = _first((args.start, 'command line'), (man.get('start'), 'manifest start'),
                                 (bs if bs in ('X', 'N', 'R') else None, 'manifest build_json.start'),
                                 (mi.group(1) if mi else None, 'manifest run_id'),
                                 (mf.group(1) if mf else None, 'folder name'))
    if start not in ('X', 'N', 'R'):
        raise AnalysisError(f'{run_dir}: the start (X, N or R) is not given by --start, the manifest or the folder name')
    replica, src['replica'] = _first((man.get('replica'), 'manifest replica'),
                                     (int(mi.group(2)) if mi else None, 'manifest run_id'),
                                     (int(mf.group(2)) if mf else None, 'folder name'))
    model, src['model'] = _first((args.model, 'command line'), (man.get('model'), 'manifest model'),
                                 (_dig(man, 'build_json', 'model'), 'manifest build_json.model'),
                                 (run_dir.parent.name, 'parent folder name'))
    # system folder: command line; else the manifest's, if it is on this machine; else the folder of the conventions
    conv = run_dir.parent.parent.parent / 'systems' / str(model) / start
    if args.system:
        sysdir, src['system_dir'] = Path(args.system), 'command line'
    else:
        cands = []
        if isinstance(man.get('system_dir'), str):
            q = Path(man['system_dir'])
            cands.append((q if q.is_absolute() else run_dir / q, 'manifest system_dir'))
        cands.append((conv, 'convention'))
        sysdir, src['system_dir'] = next(((q, w) for q, w in cands if (q / 'selections.json').is_file()), cands[0])
        if src['system_dir'] == 'convention' and len(cands) == 2:
            src['system_dir'] = 'convention (the system_dir of the manifest is not on this machine)'
    try:
        rep_int = int(replica)
    except (TypeError, ValueError):
        rep_int = None
    vr = CONST['validation_run']
    role = 'validation' if (start == vr['start'] and rep_int == vr['replica']) else 'production'
    role_by_name = None                                       # what the folder name alone would say
    if mf:
        role_by_name = 'validation' if (mf.group(1) == vr['start'] and int(mf.group(2)) == vr['replica']) else 'production'
    man_planned, man_planned_src = _first((man.get('planned_ns'), 'manifest planned_ns'),
                                          (_dig(man, 'production', 'target_ns'), 'manifest production.target_ns'))
    if role == 'validation':                                  # fixed by the specification, not by the manifest
        planned, src['planned_ns'] = _first((args.planned_ns, 'command line'),
                                            (vr['planned_ns'], 'validation run X0: fixed by the specification'))
    else:
        planned, src['planned_ns'] = _first(
            (args.planned_ns, 'command line'), (man_planned, man_planned_src),
            (_dig(man, 'protocol', 'production', 'length_ns'), 'manifest protocol.production.length_ns'),
            (CONST['planned_ns_default'], 'default'))
    dts, src['solute_interval_ps'] = _first(
        (args.solute_interval_ps, 'command line'), (man.get('solute_interval_ps'), 'manifest solute_interval_ps'),
        (_dig(man, 'protocol', 'output', 'solute_frame_interval_ps'), 'manifest protocol.output.solute_frame_interval_ps'),
        (CONST['solute_interval_ps_default'], 'default'))
    dtf, src['full_interval_ps'] = _first(
        (args.full_interval_ps, 'command line'), (man.get('full_interval_ps'), 'manifest full_interval_ps'),
        (_dig(man, 'protocol', 'output', 'full_frame_interval_ps'), 'manifest protocol.output.full_frame_interval_ps'),
        (CONST['full_interval_ps_default'], 'default'))
    return {'run_dir': run_dir, 'manifest': man, 'manifest_path': mp if mp.is_file() else None, 'start': start,
            'replica': replica, 'model': str(model), 'system_dir': sysdir.resolve(), 'label': f'{model}/{run_dir.name}',
            'planned_ns': float(planned), 'solute_interval_ps': float(dts), 'full_interval_ps': float(dtf),
            'sources': src, 'start_by_folder_name': mf.group(1) if mf else None, 'role': role,
            'role_by_folder_name': role_by_name,
            'manifest_planned_ns': man_planned, 'manifest_planned_ns_source': man_planned_src}


def ring_place_settings(man, ref_lig, warn):
    """SPEC section 9D: the site atoms and the cut-off of the count that the driver acts on in its hydration check.
    Both come from the manifest of md/run_md.py (hydration_check.settings: site_atoms, water_cutoff_A) if both are
    usable there; otherwise both are the defaults of the specification.  A run whose settings are not the
    defaults is reported with its own settings, and a warning says so."""
    atoms0 = list(CONST['ring_place_atoms_default'])
    cut0 = float(CONST['ring_place_water_cutoff_default_A'])
    out = {'note': 'section 9D: n_water_ring_place counts the water oxygens within water_cutoff_A of the crystal '
                   'positions of site_atoms, the quantity of the hydration check of md/run_md.py',
           'site_atoms': atoms0, 'water_cutoff_A': cut0, 'source': 'specification (defaults)',
           'settings_are_the_defaults': True}
    st = _dig(man, 'hydration_check', 'settings')
    if st is None:
        return out
    atoms, cut = (st.get('site_atoms'), st.get('water_cutoff_A')) if isinstance(st, dict) else (None, None)
    usable = (isinstance(atoms, list) and len(atoms) > 0 and all(isinstance(a, str) and a in ref_lig for a in atoms)
              and len(set(atoms)) == len(atoms) and isinstance(cut, (int, float)) and not isinstance(cut, bool)
              and math.isfinite(cut) and cut > 0)
    if not usable:
        warn('the manifest holds hydration_check.settings without a usable site_atoms (names of ligand heavy atoms, '
             'each once) and water_cutoff_A (a positive number); water in the place of the crystal ring is counted '
             f'with the defaults of the specification ({", ".join(atoms0)}; {cut0:g} A)')
        return out
    same = set(atoms) == set(atoms0) and abs(float(cut) - cut0) <= 1e-9
    out.update(site_atoms=[str(a) for a in atoms], water_cutoff_A=float(cut), source='manifest hydration_check.settings',
               settings_are_the_defaults=bool(same))
    if not same:
        warn(f'the hydration check of this run was made with the site atoms {", ".join(atoms)} and the cut-off '
             f'{float(cut):g} A (manifest), not with those of the plan ({", ".join(atoms0)}; {cut0:g} A); '
             'n_water_ring_place follows the run')
    return out


def load_system(sysdir, solute_pdb):
    sel = json.loads((sysdir / 'selections.json').read_text())
    sol = [int(i) for i in sel['solute']]
    pos = {p: k for k, p in enumerate(sol)}                   # prmtop index -> position in the solute trajectory
    names, resn, _, _ = read_pdb(solute_pdb)
    if len(names) != len(sol):
        raise AnalysisError(f'{solute_pdb}: {len(names)} atoms, selections.json "solute" has {len(sol)}')
    need = {k: sel[k] for k in SITE}
    need.update({n: sel['lig_heavy'][n] for n in LIG_HEAVY})
    need['HN2'] = sel['lig_HN2']
    optional = [k for k in SITE_OPTIONAL if sel.get(k) is not None]
    need.update({k: sel[k] for k in optional})
    missing = [k for k, v in need.items() if v not in pos]
    if missing:
        raise AnalysisError(f'atoms not in the "solute" selection: {missing}')
    bad = []
    for k, v in need.items():                                 # guard: the names at the mapped positions
        exp = SITE.get(k) or SITE_OPTIONAL.get(k) or k
        got = names[pos[v]]
        ok = got.upper().startswith('H') if k == 'HN2' else got.upper() == exp.upper()
        if not ok:
            bad.append(f'{k}: expected {exp}, solute.pdb has {got}')
    ca = [int(i) for i in sel['protein_CA']]
    bad += [f'protein_CA index {i}: solute.pdb has {names[pos[i]]}' for i in ca if names[pos[i]] != 'CA'][:5]
    resnum = [int(i) for i in sel['residue_numbers_deposited_for_CA']]
    if len(resnum) != len(ca):
        raise AnalysisError('residue_numbers_deposited_for_CA and protein_CA differ in length')

    rn, rr, _, rx = read_pdb(sysdir / 'reference.pdb')        # reference: CA atoms in file order, ligand by name
    ica = [i for i, n in enumerate(rn) if n == 'CA']
    if len(ica) != len(ca):
        raise AnalysisError(f'reference.pdb has {len(ica)} CA atoms, selections.json has {len(ca)}')
    ref_lig = {}
    for n in LIG_HEAVY:
        hit = [i for i, m in enumerate(rn) if m == n]
        if len(hit) != 1:
            raise AnalysisError(f'reference.pdb: {len(hit)} atoms named {n} (one expected)')
        ref_lig[n] = rx[hit[0]]
    ca_resn = [resn[pos[i]] for i in ca]
    bad += [f'CA {k}: residue {ca_resn[k]} in solute.pdb, {rr[ica[k]]} in reference.pdb'
            for k in range(len(ca)) if ca_resn[k] != rr[ica[k]]][:5]
    return {'sel': sel, 'pos': pos, 'need': need, 'ca': ca, 'ca_resnum': resnum, 'ca_resname': ca_resn,
            'ref_ca': rx[ica], 'ref_lig': ref_lig, 'name_mismatches': bad, 'n_solute': len(sol),
            'lig_all': [int(i) for i in sel['lig_all']], 'has_asn67_ND2': 'asn67_ND2' in optional}


def open_segment(path, warn):
    """Header and reader of one DCD segment; (header, None) when the file holds no complete frame yet."""
    from MDAnalysis.coordinates.DCD import DCDReader
    from MDAnalysis.lib.formats.libdcd import DCDFile
    try:
        with DCDFile(str(path)) as f:
            hdr, nfr = dict(f.header), len(f)
    except (OSError, IOError, ValueError) as e:
        warn(f'{path.name}: not readable as a DCD file ({e}); taken as holding no frame')
        return {}, None
    if nfr == 0:
        return hdr, None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')          # MDAnalysis announces a change of its version 3
        return hdr, DCDReader(str(path))


def segment_info(num, path, hdr, rd, n_frames):
    return {'index': num, 'file': path.name, 'n_frames': n_frames, 'bytes': path.stat().st_size,
            'header_first_step': int(hdr.get('istart', -1)), 'header_step_interval': int(hdr.get('nsavc', -1)),
            'header_interval_ps': r6(rd.dt) if rd is not None else None,
            'header_periodic': bool(hdr.get('is_periodic', 0))}


def read_solute(segs, S, warn, allow_empty=False):
    pos = S['pos']
    keys = list(S['need'])
    atom_idx = np.array([pos[S['need'][k]] for k in keys] + [pos[i] for i in S['ca']] + [pos[i] for i in S['lig_all']])
    nk, nca = len(keys), len(S['ca'])
    P, D, segid, info = [], [], [], []
    for num, path in segs:
        hdr, rd = open_segment(path, warn)
        if rd is None:
            info.append(segment_info(num, path, hdr, None, 0))
            continue
        if rd.n_atoms != S['n_solute']:
            raise AnalysisError(f'{path}: {rd.n_atoms} atoms per frame, the solute has {S["n_solute"]}')
        n0 = len(P)
        for ts in rd:
            P.append(ts.positions[atom_idx].astype(float))
            D.append(np.array(ts.dimensions, float) if ts.dimensions is not None else np.full(6, np.nan))
            segid.append(num)
        info.append(segment_info(num, path, hdr, rd, len(P) - n0))
        rd.close()
    if not P:
        if not allow_empty:
            raise AnalysisError('no complete frame in the solute trajectory yet')
        P, D = np.zeros((0, len(atom_idx), 3)), np.zeros((0, 6))
    P, D = np.array(P), np.array(D)
    A = {k: P[:, i, :].copy() for i, k in enumerate(keys)}
    return A, P[:, nk:nk + nca, :].copy(), P[:, nk + nca:, :].copy(), D, np.array(segid), info


def image_solute(A, lig_all_xyz, D):
    """SPEC section 2, step 1: the zinc and the ligand are brought to the periodic image next to the protein."""
    n = len(D)
    lig_keys = list(LIG_HEAVY) + ['HN2']
    n_zn = n_lig = n_box = 0
    for f in range(n):
        d = D[f]
        if not (np.all(np.isfinite(d)) and np.all(d[:3] > 0)):
            continue
        n_box += 1
        B = box_matrix(d)
        s = min_image_shift(A['zn'][f] - A['his94_NE2'][f], B)
        if np.linalg.norm(s) > 1e-6:
            A['zn'][f] += s
            n_zn += 1
        s = min_image_shift(A['N2'][f] - A['zn'][f], B)
        if np.linalg.norm(s) > 1e-6:
            for k in lig_keys:
                A[k][f] += s
            lig_all_xyz[f] += s
            n_lig += 1
    return {'frames_with_box': n_box, 'frames_without_box': n - n_box, 'frames_zinc_moved': n_zn,
            'frames_ligand_moved': n_lig}


def read_full(segs, S, full_pdb, warn, place, what='full trajectory'):
    """Quantities of the full-system trajectory (SPEC sections 6, 7 and 9D), minimum-image distances."""
    from MDAnalysis.lib.distances import distance_array
    sel = S['sel']
    ica = np.array(S['ca'], int)                             # section 9D: the CA atoms, for the fit of the frame
    ring_counts = (                                          # site atoms, cut-off, the two columns
        (list(place['site_atoms']), float(place['water_cutoff_A']), 'n_water_ring_place', 'd_water_ring_place_min_A'),
        (list(CONST['ring_wide_atoms']), float(CONST['ring_wide_water_distance_max_A']),
         'n_water_ring_wide_3p5A', 'd_water_ring_wide_min_A'))
    names, _, _, _ = read_pdb(full_pdb)
    wat = np.array(sel['water_O'], int)
    prot_o = np.array([i for i in sel['protein_heavy'] if names[i].upper().startswith('O')], int)
    bad = [f'water_O index {i}: full.pdb has {names[i]}' for i in wat if not names[i].upper().startswith('O')][:5]
    izn, in1, in2 = sel['zn'], sel['lig_heavy']['N1'], sel['lig_heavy']['N2']
    io = np.array([sel['lig_heavy']['O1'], sel['lig_heavy']['O2']], int)
    cols = {c: [] for c, _ in ALL_FULL_COLUMNS}
    segid, info, n_nobox = [], [], 0
    cut5, cutw = CONST['fifth_ligand_distance_max_A'], CONST['water_near_N1_distance_max_A']
    for num, path in segs:
        hdr, rd = open_segment(path, warn)
        if rd is None:
            info.append(segment_info(num, path, hdr, None, 0))
            continue
        if rd.n_atoms != len(names):
            raise AnalysisError(f'{path}: {rd.n_atoms} atoms per frame, full.pdb has {len(names)}')
        n0 = len(segid)
        for ts in rd:
            x = ts.positions
            box = ts.dimensions
            if box is None or not np.all(np.asarray(box)[:3] > 0):
                box = None
                n_nobox += 1
            zn = x[izn][None, :]

            def mind(idx):
                if len(idx) == 0:
                    return float('nan'), -1
                d = distance_array(zn, x[idx], box=box)[0]
                k = int(np.argmin(d))
                return float(d[k]), int(idx[k])
            dw, iw = mind(wat)
            dp, ip = mind(prot_o)
            ds, _ = mind(io)
            dn = distance_array(x[in1][None, :], x[wat], box=box)[0] if len(wat) else np.zeros(0)
            cols['d_Zn_Owater_min_A'].append(dw)
            cols['water_O_index_nearest_Zn'].append(iw)
            cols['d_Zn_Oprotein_min_A'].append(dp)
            cols['protein_O_index_nearest_Zn'].append(ip)
            cols['d_Zn_Osulfonyl_min_full_A'].append(ds)
            cols['d_Zn_O_any_min_A'].append(float(np.nanmin([dw, dp, ds])))
            cols['fifth_water'].append(int(dw <= cut5))
            cols['fifth_proteinO'].append(int(dp <= cut5))
            cols['fifth_any'].append(int(np.nanmin([dw, dp, ds]) <= cut5))
            cols['n_water_N1_3p5A'].append(int(np.sum(dn <= cutw)))
            cols['d_Zn_N2_full_A'].append(float(distance_array(zn, x[in2][None, :], box=box)[0, 0]))
            # section 9D: crystal positions (reference.pdb) carried into the frame by the inverse of the CA fit of
            # section 2; a water oxygen counts once if it is within the cut-off of any of them (minimum image).
            # ring_place: site atoms and cut-off of the driver's hydration check; ring_wide: seven positions, 3.5 A
            Qf, cxf, crf = kabsch(x[ica][None, :, :], S['ref_ca'])
            for site_atoms, cut_site, col_n, col_d in ring_counts:
                ref_site = np.array([S['ref_lig'][a] for a in site_atoms], float)
                site = (ref_site - crf) @ Qf[0].T + cxf[0, 0]
                if len(wat):
                    dsite = distance_array(site.astype(np.float32), x[wat], box=box).min(axis=0)
                    cols[col_n].append(int(np.sum(dsite <= cut_site)))
                    cols[col_d].append(float(dsite.min()))
                else:
                    cols[col_n].append(0)
                    cols[col_d].append(float('nan'))
            segid.append(num)
        info.append(segment_info(num, path, hdr, rd, len(segid) - n0))
        rd.close()
    if n_nobox:
        warn(f'{what}: {n_nobox} frames without a periodic box; distances there are not minimum-image')
    return ({c: np.array(v) for c, v in cols.items()}, np.array(segid, int), info,
            {'n_water_O': int(len(wat)), 'n_protein_O': int(len(prot_o)), 'frames_without_box': n_nobox,
             'name_mismatches': bad})


def header_check(info, interval_ps, warn, what):
    """Guards of SPEC section 1 from the DCD headers: the interval, and, where the header counts the steps of
    the whole run (md/run_md.py does; OpenMM's own reporter starts every file at one interval and tells nothing),
    that each segment starts one interval after the last frame of the segment before it."""
    out = {'interval_matches': True, 'contiguous': None, 'details': []}
    segs = [a for a in info if a['n_frames'] > 0]
    for a in segs:
        if a['header_interval_ps'] is None or abs(a['header_interval_ps'] - interval_ps) > 1e-3 * interval_ps:
            out['interval_matches'] = False
            out['details'].append(f'segment {a["index"]}: header interval {a["header_interval_ps"]} ps, '
                                  f'analysis uses {interval_ps} ps')
    seen = []
    for a, b in zip(segs[:-1], segs[1:]):
        exp = a['header_first_step'] + a['n_frames'] * a['header_step_interval']
        if b['header_first_step'] == exp:
            seen.append(True)
        elif b['header_first_step'] == b['header_step_interval']:
            seen.append(None)                                 # a file that counts from its own start: no information
        else:
            seen.append(False)
            out['details'].append(f'segment {b["index"]}: first step {b["header_first_step"]} in the header, '
                                  f'{exp} expected after segment {a["index"]}')
    if False in seen:
        out['contiguous'] = False
    elif seen and all(v is True for v in seen):
        out['contiguous'] = True
    if not out['interval_matches']:
        warn(f'{what} trajectory: the interval in the DCD header differs from the interval used ('
             + '; '.join(d for d in out['details'] if 'interval' in d) + ')')
    if out['contiguous'] is False:
        warn(f'{what} trajectory: by their DCD headers the segments leave a gap or overlap ('
             + '; '.join(d for d in out['details'] if 'first step' in d) + '): the time axis is wrong')
    return out


def manifest_segment_check(man, sinfo, finfo, dt_ps, warn):
    """Guard of SPEC section 1 from the manifest of md/run_md.py (production.segments): the frames it records
    for a segment against the frames read, and its start time against the frames before it."""
    recs = _dig(man, 'production', 'segments')
    out = {'available': isinstance(recs, list), 'consistent': None, 'details': []}
    if not out['available']:
        return out
    by = {}
    for r in recs:
        m = re.fullmatch(r'seg(\d+)', str(r.get('tag', ''))) if isinstance(r, dict) else None
        if m:
            by[int(m.group(1))] = r
    for what, info in (('solute', sinfo), ('full', finfo)):
        before = 0
        for a in info:
            r = by.get(a['index'])
            if r is not None:
                nrec = _dig(r, 'frames', what)
                if nrec is not None and int(nrec) != a['n_frames']:
                    out['details'].append(f'{a["file"]}: {a["n_frames"]} frames read, the manifest records {nrec}')
                t0 = r.get('time_start_ps')
                if what == 'solute' and t0 is not None and abs(float(t0) - before * dt_ps) > 1e-6 * max(1.0, before * dt_ps):
                    out['details'].append(f'{a["file"]}: the manifest starts it at {t0} ps, the frames before it end '
                                          f'at {before * dt_ps:g} ps')
            before += a['n_frames']
    have = {a['index'] for a in sinfo}
    out['details'] += [f'segment {k} is recorded in the manifest and has no solute file' for k in sorted(by) if k not in have]
    out['consistent'] = not out['details']
    if not out['consistent']:
        warn('the manifest and the trajectory files disagree (' + '; '.join(out['details']) + ')')
    return out


def overlap_check(ca_raw, segid, warn):
    """Guard (SPEC section 1): a segment that repeats frames of the segment before it (a restart from an older
    checkpoint).  The first frame of each segment is compared with every frame of the previous segment; CA
    coordinates equal within 0.001 A mean the same frame was written twice."""
    out = {'suspected': False, 'boundaries': []}
    ids = list(dict.fromkeys(int(i) for i in segid))
    for a, b in zip(ids[:-1], ids[1:]):
        fa, fb = np.where(segid == a)[0], np.where(segid == b)[0]
        d = np.abs(ca_raw[fa] - ca_raw[fb[0]]).max(axis=(1, 2))
        k = int(np.argmin(d))
        rec = {'segments': [a, b], 'smallest_max_abs_diff_A': r6(d[k]), 'frame_of_previous_segment': int(fa[k]),
               'repeated': bool(d[k] < 1e-3)}
        out['boundaries'].append(rec)
        if rec['repeated']:
            out['suspected'] = True
            warn(f'segment {b} starts with a frame that segment {a} already holds (frame {int(fa[k])}): the segments '
                 f'overlap, the time axis is wrong')
    return out


def _window(ref, tol):
    """Closed interval [ref - tol, ref + tol], rounded to six decimals so that a limit is the number of the SPEC."""
    return [round(ref - tol, 6), round(ref + tol, 6)]


def gate_gv(stats, fr, role, complete):
    """Gate G-V (SPEC section 11).  The criteria are evaluated for every run; they are a gate, with a verdict,
    for the validation run (X0) only."""
    g = CONST['gate_GV']
    applies = role == 'validation'
    out = {'applies_to_this_run': applies,
           'note': 'gate G-V of md/PLAN.md is the gate of the validation run X0 (5 ns from the crystal pose); for '
                   'every other run the same numbers are given without a verdict'}
    m = stats['lig_rmsd_A']['mean']
    out['GV1_lig_rmsd'] = {'mean_A': m, 'limit_max_A': g['lig_rmsd_mean_max_A'], 'passed': m <= g['lig_rmsd_mean_max_A']}
    h = fr['hb_N2H_T199OG1']
    num, den = g['hb_N2H_T199OG1_fraction_min']
    out['GV2_hbond'] = {'fraction': h['fraction'], 'count': h['count'], 'n': h['n'], 'limit_min': num / den,
                        'passed': den * h['count'] >= num * h['n']}
    lo, hi = g['d_Zn_N2_mean_window_A']
    z = {'d_Zn_N2_A': {'mean_A': stats['d_Zn_N2_A']['mean'], 'window_A': [lo, hi],
                       'passed': lo <= stats['d_Zn_N2_A']['mean'] <= hi}}
    lo, hi = g['d_Zn_NHis_mean_window_A']
    for c in ZN_HIS:
        z[c] = {'mean_A': stats[c]['mean'], 'window_A': [lo, hi], 'passed': lo <= stats[c]['mean'] <= hi}
    tol = g['ang_N_Zn_N_tolerance_deg']
    for c in ZN_ANGLES:                                       # each mean angle against its own value in 6RL9
        ref = g['ang_N_Zn_N_reference_6RL9_deg'][c]
        lo, hi = _window(ref, tol)
        z[c] = {'mean_deg': stats[c]['mean'], 'reference_6RL9_deg': ref, 'tolerance_deg': tol,
                'deviation_deg': r6(stats[c]['mean'] - ref), 'window_deg': [lo, hi], 'passed': lo <= stats[c]['mean'] <= hi}
    # fifth ligand: two fractions, one per trajectory, each at most 5 % (integers: 20 x count <= frames)
    num, den = g['fifth_ligand_fraction_max']
    fs, ff = fr['fifth_sulfonylO'], fr['fifth_any']
    have_full = ff['n'] > 0
    ok_s = den * fs['count'] <= num * fs['n']
    ok_f = bool(have_full and den * ff['count'] <= num * ff['n'])
    z['fifth_ligand'] = {
        'limit_max_fraction': num / den,
        'solute_trajectory': {'flag': 'fifth_sulfonylO', 'count': fs['count'], 'n': fs['n'], 'fraction': fs['fraction'],
                              'passed': ok_s},
        'full_trajectory': {'flag': 'fifth_any', 'available': have_full, 'count': ff['count'] if have_full else None,
                            'n': ff['n'], 'fraction': ff['fraction'], 'passed': ok_f},
        'passed': bool(ok_s and ok_f)}
    out['GV3_zinc'] = {'criteria': z, 'passed': all(v['passed'] for v in z.values())}
    k = {}
    for c, name in (('ang_C4_S_N2_deg', 'ang_C4_S_N2'), ('tau_fold_deg', 'tau_fold')):
        ref, tol = g[f'{name}_reference_deg'], g[f'{name}_tolerance_deg']
        lo, hi = _window(ref, tol)
        k[c] = {'mean_deg': stats[c]['mean'], 'sd_deg': stats[c]['sd'], 'reference_deg': ref, 'tolerance_deg': tol,
                'deviation_deg': r6(stats[c]['mean'] - ref), 'window_deg': [lo, hi], 'passed': lo <= stats[c]['mean'] <= hi}
    out['GV4_ligand_internal'] = {'criteria': k, 'passed': all(v['passed'] for v in k.values())}
    out['passed'] = bool(out['GV1_lig_rmsd']['passed'] and out['GV2_hbond']['passed'] and out['GV3_zinc']['passed']
                         and out['GV4_ligand_internal']['passed'])
    out['verdict'] = ('PASSED' if out['passed'] else 'FAILED') if applies else None
    out['provisional'] = not complete
    return out


def write_csv(path, head, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)


# ----------------------------------------------------------------------------------------------------------
# Equilibration frames (SPEC section 9E): a separate time series; nothing of it enters a statistic or a rule
# ----------------------------------------------------------------------------------------------------------
def equilibration_stage_table(man):
    """Start and end (ps, counted from the start of the equilibration) of every step: from the manifest of
    md/run_md.py if it lists the steps (the stages eqNN and the extension steps eqNNxM), else the six stages of
    md/protocol.json.  Returns ({step code: (start, end)}, source)."""
    recs = _dig(man, 'equilibration', 'stages')
    if isinstance(recs, list):
        tab = {}
        for r in recs:
            m = re.fullmatch(EQ_STEP, str(r.get('tag', ''))) if isinstance(r, dict) else None
            try:
                if m:
                    tab[step_code(m.group(1), m.group(2))] = (float(r['time_start_ps']), float(r['time_end_ps']))
            except (KeyError, TypeError, ValueError):
                pass
        return tab, 'manifest equilibration.stages'
    return ({step_code(k): (float(a), float(b)) for k, (a, b) in CONST['equilibration_stages_protocol_ps'].items()},
            'protocol (constants of the specification)')


def driver_hydration_checks(man):
    """What the driver recorded of its own count of water in the place of the crystal ring at the end of a step
    (manifest equilibration.stages[].hydration_check, md/run_md.py of 6 Oct 2026): {step code: {...}}."""
    out = {}
    recs = _dig(man, 'equilibration', 'stages')
    for r in recs if isinstance(recs, list) else []:
        m = re.fullmatch(EQ_STEP, str(r.get('tag', ''))) if isinstance(r, dict) else None
        h = r.get('hydration_check') if m else None
        if isinstance(h, dict) and isinstance(h.get('n_water_O'), int) and not isinstance(h.get('n_water_O'), bool):
            out[step_code(m.group(1), m.group(2))] = {
                'n_water_O': h['n_water_O'], 'nearest_water_O_A': _dig(h, 'nearest_water_O', 'distance_A'),
                'applies': h.get('applies'), 'decision': h.get('decision')}
    return out


def expected_frames(start_ps, end_ps, interval_ps):
    """Frames of a stage: the largest m with m x interval <= length (the driver writes at interval, 2 x interval, ...)."""
    return int(math.floor((end_ps - start_ps) / interval_ps * (1.0 + 1e-9) + 1e-9))


def analyse_equilibration(run_dir, R, S, out_dir, solute_columns, start_state, start_sign, warn, place):
    names = ('timeseries_equil_solute.csv', 'timeseries_equil_full.csv')
    for name in names:                                       # no file of an earlier analysis is left behind
        if (out_dir / name).is_file():
            (out_dir / name).unlink()
    sfiles, ffiles = equilibration_files(run_dir, 'solute'), equilibration_files(run_dir, 'full')
    out = {'available': False,
           'note': 'equilibration frames (solute_eqNN.dcd, full_eqNN.dcd): reported in a time series of their own, '
                   'times in ns counted from the start of the production (negative); no statistic, verdict or gate '
                   'of the run uses them',
           'files_solute': [p.name for _, p in sfiles], 'files_full': [p.name for _, p in ffiles]}
    if not sfiles and not ffiles:
        return out, None
    table, source = equilibration_stage_table(R['manifest'])
    out['stage_table_source'] = source
    if not table:
        warn('equilibration: the manifest lists no stage with a start and an end; the equilibration files are not read')
        return out, None
    t_eq = max(e for _, e in table.values())
    out['total_ps'] = r6(t_eq)
    for what, files in (('solute', sfiles), ('full', ffiles)):
        left = [p.name for nn, p in files if nn not in table]
        if left:
            warn(f'equilibration: {", ".join(left)} belong to no step of the table used ({source}); not read')
    suse = [(nn, p) for nn, p in sfiles if nn in table]
    fuse = [(nn, p) for nn, p in ffiles if nn in table]
    have_full = bool(fuse) and (run_dir / 'full.pdb').is_file()

    def times(info, interval_ps, what):
        t, stage = [], []
        for a in info:
            code = a.pop('index')
            a.update(tag=step_tag(code), stage=code // 1000, extension=code % 1000)
            s0, s1 = table[code]
            t += [(s0 + (k + 1) * interval_ps - t_eq) / 1000.0 for k in range(a['n_frames'])]
            stage += [code] * a['n_frames']
            a['n_frames_expected'] = expected_frames(s0, s1, interval_ps)
            if a['n_frames'] != a['n_frames_expected']:
                warn(f'equilibration: {a["file"]} holds {a["n_frames"]} frames, {a["n_frames_expected"]} expected for a '
                     f'stage from {s0:g} to {s1:g} ps with {interval_ps:g} ps between frames')
            if a['n_frames'] and a['header_interval_ps'] is not None and a['header_step_interval'] > 0:
                if abs(a['header_interval_ps'] - interval_ps) > 1e-3 * interval_ps:
                    warn(f'equilibration: {a["file"]}: header interval {a["header_interval_ps"]} ps, the analysis uses '
                         f'{interval_ps} ps')
                th = a['header_first_step'] * a['header_interval_ps'] / a['header_step_interval']
                a['header_first_frame_ps'] = r6(th)
                tol = 1e-3 * interval_ps                      # a header that counts from the start of its own file says nothing
                if abs(th - (s0 + interval_ps)) > tol and abs(th - interval_ps) > tol:
                    warn(f'equilibration: by its DCD header the first frame of {a["file"]} is at {th:g} ps, the stage '
                         f'table puts it at {s0 + interval_ps:g} ps: the time axis of the equilibration is wrong')
        return np.array(t, float), np.array(stage, int)

    A, ca, lig_all_xyz, D, _, sinfo = read_solute(suse, S, warn, allow_empty=True)
    ts, sstage = times(sinfo, R['solute_interval_ps'], 'solute')
    ns = len(ts)
    cols = [(c, k) for c, k in solute_columns if c not in HISTORY_COLUMNS]
    q = None
    if ns:
        imaging = image_solute(A, lig_all_xyz, D)
        q, _ = solute_quantities(A, ca, S['ref_ca'], S['ref_lig'])
        out['imaging_solute'] = imaging
    if suse:
        write_csv(out_dir / names[0], ['frame', 'stage', 'extension', 'time_ns'] + [c for c, _ in cols],
                  ([i, int(sstage[i]) // 1000, int(sstage[i]) % 1000, '%.6f' % ts[i]] + [FMT[k] % q[c][i] for c, k in cols]
                   for i in range(ns)))
    fq, tf, fstage, finfo, nf = None, np.zeros(0), np.zeros(0, int), [], 0
    if have_full:
        fq, _, finfo, _ = read_full(fuse, S, run_dir / 'full.pdb', warn, place, what='equilibration, full trajectory')
        tf, fstage = times(finfo, R['full_interval_ps'], 'full')
        nf = len(tf)
        fcols = ALL_FULL_COLUMNS
        write_csv(out_dir / names[1], ['frame', 'stage', 'extension', 'time_ns'] + [c for c, _ in fcols],
                  ([i, int(fstage[i]) // 1000, int(fstage[i]) % 1000, '%.6f' % tf[i]] + [FMT[k] % fq[c][i] for c, k in fcols]
                   for i in range(nf)))
    elif fuse:
        warn('equilibration: full_eqNN.dcd without full.pdb; the full frames of the equilibration are not read')

    opp = 'R' if start_state == 'C' else 'C'
    stages = []
    by_s = {step_code(a['stage'], a['extension']): a for a in sinfo}
    by_f = {step_code(a['stage'], a['extension']): a for a in finfo}
    driver = driver_hydration_checks(R['manifest'])
    for nn in sorted(set(by_s) | set(by_f)):
        s0, s1 = table[nn]
        rec = {'stage': nn // 1000, 'extension': nn % 1000, 'tag': step_tag(nn), 'time_start_ps': r6(s0), 'time_end_ps': r6(s1),
               'time_start_ns': r6((s0 - t_eq) / 1000.0), 'time_end_ns': r6((s1 - t_eq) / 1000.0),
               'solute': by_s.get(nn), 'full': by_f.get(nn),
               'count_C': None, 'count_R': None, 'count_I': None, 'count_n2_positive': None, 'count_n2_negative': None,
               'last_solute_frame': None, 'last_full_frame': None}
        k = np.where(sstage == nn)[0]
        if nn in by_s:
            st, sg = (q['state'][k], q['n2_sign'][k]) if len(k) else (np.zeros(0, 'U1'), np.zeros(0, int))
            rec.update(count_C=int(np.sum(st == 'C')), count_R=int(np.sum(st == 'R')), count_I=int(np.sum(st == 'I')),
                       count_n2_positive=int(np.sum(sg > 0)), count_n2_negative=int(np.sum(sg < 0)))
        if len(k):
            i = int(k[-1])
            rec['last_solute_frame'] = {'frame': i, 'time_ns': r6(ts[i]), 'state': str(q['state'][i]),
                                        'lig_rmsd_A': r6(q['lig_rmsd_A'][i]), 'n2_volume': r6(q['n2_volume'][i]),
                                        'n2_sign': int(q['n2_sign'][i])}
        k = np.where(fstage == nn)[0]
        if len(k):
            i = int(k[-1])
            rec['last_full_frame'] = {'frame': i, 'time_ns': r6(tf[i]),
                                      'n_water_ring_place': int(fq['n_water_ring_place'][i]),
                                      'd_water_ring_place_min_A': r6(fq['d_water_ring_place_min_A'][i]),
                                      'n_water_ring_wide_3p5A': int(fq['n_water_ring_wide_3p5A'][i]),
                                      'd_water_ring_wide_min_A': r6(fq['d_water_ring_wide_min_A'][i])}
        # the driver's own count at the end of the step (its hydration check) beside n_water_ring_place of section 9D
        # in the full frame written at that time; compared only if the step ends on a full frame
        rec['driver_hydration_check'] = None
        if nn in driver:
            d = dict(driver[nn], compared=False, count_agrees=None, nearest_agrees=None)
            lf = rec['last_full_frame']
            m_full = by_f[nn]['n_frames'] if nn in by_f else 0        # the step ends on a full frame: S + m x interval = E
            if lf is not None and m_full >= 1 and abs(s0 + m_full * R['full_interval_ps'] - s1) <= 1e-6:
                d.update(compared=True, count_agrees=bool(lf['n_water_ring_place'] == d['n_water_O']))
                if (isinstance(d['nearest_water_O_A'], (int, float)) and not isinstance(d['nearest_water_O_A'], bool)
                        and lf['d_water_ring_place_min_A'] is not None):
                    d['nearest_agrees'] = bool(abs(lf['d_water_ring_place_min_A'] - float(d['nearest_water_O_A'])) <= 0.002)
                if d['count_agrees'] is False or d['nearest_agrees'] is False:
                    warn(f'equilibration: at the end of {rec["tag"]} the driver counted {d["n_water_O"]} water oxygens in the '
                         f'place of the crystal ring (nearest {d["nearest_water_O_A"]} A), the analysis '
                         f'{lf["n_water_ring_place"]} (nearest {lf["d_water_ring_place_min_A"]} A) in the full frame of that '
                         f'time, with the site atoms {", ".join(place["site_atoms"])} and the cut-off '
                         f'{place["water_cutoff_A"]:g} A ({place["source"]})')
            rec['driver_hydration_check'] = d
        stages.append(rec)

    def first(mask):
        k = np.where(mask)[0]
        if not len(k):
            return None
        i = int(k[0])
        return {'frame': i, 'stage': int(sstage[i]) // 1000, 'extension': int(sstage[i]) % 1000, 'time_ns': r6(ts[i]),
                'state': str(q['state'][i]), 'n2_sign': int(q['n2_sign'][i])}
    out.update(
        available=bool(ns or nf), n_frames_solute=int(ns), n_frames_full=int(nf),
        files_written=[n for n, ok in zip(names, (bool(suse), have_full)) if ok],
        state_of_start=start_state, n2_sign_of_start=int(start_sign), stages=stages,
        count_C=int(np.sum(q['state'] == 'C')) if ns else 0, count_R=int(np.sum(q['state'] == 'R')) if ns else 0,
        count_I=int(np.sum(q['state'] == 'I')) if ns else 0,
        count_n2_positive=int(np.sum(q['n2_sign'] > 0)) if ns else 0,
        count_n2_negative=int(np.sum(q['n2_sign'] < 0)) if ns else 0,
        first_frame_not_in_state_of_start=first(q['state'] != start_state) if ns else None,
        first_frame_in_opposite_state=first(q['state'] == opp) if ns else None,
        first_frame_n2_sign_not_of_start=first(q['n2_sign'] != start_sign) if ns else None,
        state_last_frame=str(q['state'][-1]) if ns else None, n2_sign_last_frame=int(q['n2_sign'][-1]) if ns else None,
        time_last_solute_frame_ns=r6(ts[-1]) if ns else None)
    series = {'t': ts, 'q': q, 'tf': tf, 'fq': fq} if (ns or nf) else None
    return out, series


def analyse_run(run_dir, args):
    warnings_list = []

    def warn(msg):
        warnings_list.append(msg)
        print(f'  WARNING: {msg}', file=sys.stderr)

    R = identify_run(run_dir, args)
    run_dir, sysdir = R['run_dir'], R['system_dir']
    out_dir = Path(args.out).resolve() if (args.out and not args.campaign) else run_dir / 'analysis'
    out_dir.mkdir(parents=True, exist_ok=True)
    for f in ('selections.json', 'reference.pdb'):
        if not (sysdir / f).is_file():
            raise AnalysisError(f'{sysdir / f} is missing (system folder: --system, manifest "system_dir" or the convention)')
    if not (run_dir / 'solute.pdb').is_file():
        raise AnalysisError(f'{run_dir / "solute.pdb"} is missing')
    if '_error' in R['manifest']:
        warn(f'manifest.json could not be read ({R["manifest"]["_error"]})')
    S = load_system(sysdir, run_dir / 'solute.pdb')
    # guard: the system folder is the one the run was made from (checksums recorded in the manifest of md/run_md.py)
    syscheck = {}
    for f in ('selections.json', 'reference.pdb'):
        rec = _dig(R['manifest'], 'inputs', f, 'sha256')
        syscheck[f] = None if rec is None else (rec == sha256(sysdir / f))
        if syscheck[f] is False:
            warn(f'{f} of the system folder used ({sysdir}) is not the file the run was made with (checksum in the '
                 f'manifest): the atom indices may not belong to this trajectory')
    if S['name_mismatches']:
        msg = 'atom names at the mapped positions differ: ' + '; '.join(S['name_mismatches'])
        if not args.skip_name_check:
            raise AnalysisError(msg + ' (use --skip-name-check to go on)')
        warn(msg)

    # ---- solute trajectory -----------------------------------------------------------------------------
    ssegs = segment_files(run_dir, 'solute')
    if not ssegs:
        raise AnalysisError(f'{run_dir}: no solute_segNNNN.dcd')
    A, ca, lig_all_xyz, D, segid, sinfo = read_solute(ssegs, S, warn)
    n = len(segid)
    overlap = overlap_check(ca, segid, warn)
    dt_ns = R['solute_interval_ps'] / 1000.0
    time_ns = (np.arange(n) + 1) * dt_ns
    imaging = image_solute(A, lig_all_xyz, D)
    if imaging['frames_without_box']:
        warn(f'solute trajectory: {imaging["frames_without_box"]} frames without a periodic box; no imaging there')
    q, caf = solute_quantities(A, ca, S['ref_ca'], S['ref_lig'])
    solute_columns = SOLUTE_COLUMNS + V12_SOLUTE_COLUMNS + (ASN67_COLUMNS if S['has_asn67_ND2'] else ())
    full_columns = ALL_FULL_COLUMNS
    start_state = CONST['state_of_start'][R['start']]
    start_sign = CONST['n2_sign_of_start'][R['start']]
    n_planned = int(round(R['planned_ns'] / dt_ns))
    if R['role_by_folder_name'] and R['role_by_folder_name'] != R['role']:
        warn(f'by its folder name ({run_dir.name}) this is a {R["role_by_folder_name"]} run, by start {R["start"]} '
             f'({R["sources"]["start"]}) and replica {R["replica"]} ({R["sources"]["replica"]}) it is a {R["role"]} run; '
             f'it is analysed as a {R["role"]} run')
    if R['role'] == 'validation':                             # run X0: 5 ns by the specification (section 1)
        mp = R['manifest_planned_ns']
        if R['sources']['planned_ns'] != 'command line' and mp is not None and abs(float(mp) - R['planned_ns']) > 1e-9:
            warn(f'validation run: the manifest plans {float(mp):g} ns ({R["manifest_planned_ns_source"]}), the '
                 f'specification fixes {R["planned_ns"]:g} ns; the analysis uses {R["planned_ns"]:g} ns')
        if n > n_planned:
            warn(f'validation run: {n} solute frames, more than the {n_planned} of the planned {R["planned_ns"]:g} ns; '
                 f'all frames are used, and gate G-V is defined for a run of the planned length')
    q['state_assigned'], orient = orientation(q['state'], start_state, dt_ns, n_planned)
    q['n2_sign_assigned'], n2 = n2_configuration(q['n2_volume'], q['state'], start_sign, dt_ns)     # section 9A
    hb_split = hbond_first_passage(q['hb_N2H_T199OG1'], orient)                                     # section 9B
    shead = header_check(sinfo, R['solute_interval_ps'], warn, 'solute')
    if R['start_by_folder_name'] and R['start_by_folder_name'] != R['start']:
        warn(f'the folder name says start {R["start_by_folder_name"]}, the analysis uses start {R["start"]} '
             f'({R["sources"]["start"]})')
    # the protein must be whole: neighbouring CA atoms (deposited numbers i, i+1; 125 and 127 are joined)
    rn = S['ca_resnum']
    pairs = [k for k in range(len(rn) - 1) if rn[k + 1] - rn[k] == 1 or (rn[k], rn[k + 1]) == (125, 127)]
    ca_gap = float(dist(ca[:, pairs, :], ca[:, [k + 1 for k in pairs], :]).max()) if pairs else None
    if ca_gap is not None and ca_gap > 4.5:
        warn(f'neighbouring CA atoms up to {ca_gap:.2f} A apart: the protein is not whole in some frame')

    # ---- full trajectory -------------------------------------------------------------------------------
    place = ring_place_settings(R['manifest'], S['ref_lig'], warn)                                  # section 9D
    fsegs = segment_files(run_dir, 'full')
    fq, fsegid, finfo, fmeta, fhead = {c: np.zeros(0) for c, _ in full_columns}, np.zeros(0, int), [], None, None
    dtf_ns = R['full_interval_ps'] / 1000.0
    if fsegs and (run_dir / 'full.pdb').is_file():
        fq, fsegid, finfo, fmeta = read_full(fsegs, S, run_dir / 'full.pdb', warn, place)
        if fmeta['name_mismatches']:
            msg = 'full.pdb: ' + '; '.join(fmeta['name_mismatches'])
            if not args.skip_name_check:
                raise AnalysisError(msg + ' (use --skip-name-check to go on)')
            warn(msg)
        fhead = header_check(finfo, R['full_interval_ps'], warn, 'full')
    else:
        warn('no full-system trajectory (full_segNNNN.dcd with full.pdb): water quantities are missing')
    nf = len(fsegid)
    ftime_ns = (np.arange(nf) + 1) * dtf_ns
    mcheck = manifest_segment_check(R['manifest'], sinfo, finfo, R['solute_interval_ps'], warn)
    # consistency of the two trajectories: Zn-N2 of full frame j against solute frame ratio*(j+1)-1
    ratio = R['full_interval_ps'] / R['solute_interval_ps']
    cons = {'ratio_full_to_solute_interval': r6(ratio), 'n_compared': 0, 'max_abs_diff_d_Zn_N2_A': None}
    if nf and abs(ratio - round(ratio)) < 1e-9:
        k = int(round(ratio)) * (np.arange(nf) + 1) - 1
        ok = k < n
        if ok.any():
            dmax = float(np.abs(fq['d_Zn_N2_full_A'][ok] - q['d_Zn_N2_A'][k[ok]]).max())
            cons.update(n_compared=int(ok.sum()), max_abs_diff_d_Zn_N2_A=r6(dmax))
            if dmax > 0.01:
                warn(f'the full and the solute trajectory disagree (Zn-N2 differs by up to {dmax:.3f} A at matching '
                     f'times): the segments are not aligned')
        if (~ok).any():
            warn(f'{int((~ok).sum())} frames of the full trajectory lie beyond the end of the solute trajectory')

    # ---- statistics ------------------------------------------------------------------------------------
    blk, blkf = int(round(CONST['block_ns'] / dt_ns)), int(round(CONST['block_ns'] / dtf_ns))
    stats, fracs = {}, {}
    for c, kind in solute_columns:
        if kind in ('A', 'deg', 'tors', 'vol'):
            stats[c] = summarise(q[c], kind, blk)
        elif kind == 'flag':
            fracs[c] = fraction(q[c], blk)
    for c, kind in full_columns:
        if kind == 'A' or c in FULL_COUNT_COLUMNS:
            stats[c] = summarise(fq[c], 'count' if kind == 'int' else kind, blkf)
        elif kind == 'flag':
            fracs[c] = fraction(fq[c], blkf)
    for c in list(stats):
        stats[c]['trajectory'] = 'full' if c in dict(full_columns) else 'solute'
    labels = {name: int(np.sum(q['lig_rmsd_labelling'] == i)) for i, (name, _) in enumerate(LABELLINGS)}
    asn = {'available': S['has_asn67_ND2'],
           'note': 'optional key asn67_ND2 of selections.json; the flag is 1 where ND2 is nearer to the aniline N than '
                   'OD1; the orientation states use OD1 only'}
    if S['has_asn67_ND2']:
        fl = np.asarray(q['asn67_ND2_nearer_N1'], int)
        asn.update(n=int(n), count_ND2_nearer=int(fl.sum()), fraction_ND2_nearer=r6(fl.mean()),
                   by_state={k: {'n': int(np.sum(q['state'] == k)), 'count_ND2_nearer': int(fl[q['state'] == k].sum())}
                             for k in ('C', 'R', 'I')})

    mean_ca = caf.mean(axis=0)                                                        # section 9
    rmsf = np.sqrt(((caf - mean_ca) ** 2).sum(-1).mean(0))
    kmax = int(np.argmax(rmsf))

    seen_ns = n * dt_ns
    complete = n >= n_planned
    val = {'lig_rmsd_mean_A': stats['lig_rmsd_A']['mean'], 'lig_rmsd_sd_A': stats['lig_rmsd_A']['sd'],
           'lig_rmsd_mean_last_half_A': stats['lig_rmsd_A']['mean_last_half'],
           'rule1_threshold_A': CONST['rule1_rmsd_A'],
           'rule1_left_crystal': bool(stats['lig_rmsd_A']['mean_last_half'] is not None
                                      and stats['lig_rmsd_A']['mean_last_half'] > CONST['rule1_rmsd_A']),
           'rule1_applies_to_this_run': R['start'] == 'X' and R['role'] == 'production',
           'hb_N2H_T199OG1_fraction': fracs['hb_N2H_T199OG1']['fraction']}
    gv = gate_gv(stats, fracs, R['role'], complete)
    # plan, change of 6 Oct 2026, point 4: the zinc to Glu106 by name beside the gate; reported, not a criterion
    gv['reported_zn_glu106'] = {
        'note': 'reported beside gate G-V (section 9C); no criterion, and not part of "passed" or of the verdict',
        'distance_max_A': CONST['glu106_distance_max_A'],
        **{c: {'mean_A': stats[c]['mean'], 'min_A': stats[c]['min'], 'count_within': fracs[f]['count'],
               'n': fracs[f]['n'], 'fraction_within': fracs[f]['fraction']}
           for c, f in (('d_Zn_E106OE1_A', 'fifth_E106OE1'), ('d_Zn_E106OE2_A', 'fifth_E106OE2'))}}
    # section 9E: the equilibration in a series of its own, and the first production frame against the built start
    equil, eq_series = analyse_equilibration(run_dir, R, S, out_dir, solute_columns, start_state, start_sign, warn, place)
    start_prod = {
        'note': 'the first frame of the production (solute frame 0) against the start as built; reported only',
        'time_ns': r6(time_ns[0]), 'state_of_start': start_state, 'state_first_frame': str(q['state'][0]),
        'state_differs': bool(q['state'][0] != start_state),
        'state_first_frame_is_opposite': bool(q['state'][0] == orient['opposite_state']),
        'n2_sign_of_start': int(start_sign), 'n2_sign_first_frame': int(q['n2_sign'][0]),
        'n2_sign_differs': bool(q['n2_sign'][0] != start_sign), 'n2_volume_first_frame': r6(q['n2_volume'][0]),
        'lig_rmsd_first_frame_A': r6(q['lig_rmsd_A'][0])}

    ang_means = [stats[c]['mean'] for c in ZN_ANGLES]
    ang_ref = CONST['gate_GV']['ang_N_Zn_N_reference_6RL9_deg']
    row = {
        'run': R['label'], 'model': R['model'], 'start': R['start'], 'replica': R['replica'],
        'seen_ns': r6(seen_ns), 'planned_ns': r6(R['planned_ns']), 'complete': complete,
        'verdict': orient['verdict'], 'provisional': orient['provisional'],
        'frac_start_state_last20ns': orient['fraction_start_state_window'],
        'frac_opposite_state_last20ns': orient['fraction_opposite_state_window'],
        'frac_intermediate_last20ns': orient['fraction_intermediate_window'],
        'first_passage_ns': orient['first_passage_ns'], 'n_transitions': orient['n_transitions'],
        'd1_N1_F131CZ_mean_last_half_A': stats['d1_N1_F131CZ_A']['mean_last_half'],
        'd2_N1_N67OD1_mean_last_half_A': stats['d2_N1_N67OD1_A']['mean_last_half'],
        'd_N1_N67ND2_mean_last_half_A': stats['d_N1_N67ND2_A']['mean_last_half'] if S['has_asn67_ND2'] else None,
        'asn67_ND2_nearer_N1_fraction': fracs['asn67_ND2_nearer_N1']['fraction'] if S['has_asn67_ND2'] else None,
        'lig_rmsd_mean_A': stats['lig_rmsd_A']['mean'], 'lig_rmsd_sd_A': stats['lig_rmsd_A']['sd'],
        'lig_rmsd_mean_last_half_A': stats['lig_rmsd_A']['mean_last_half'],
        'lig_rmsd_sd_last_half_A': stats['lig_rmsd_A']['sd_last_half'],
        'rule1_left_crystal': val['rule1_left_crystal'],
        'hb_N2H_T199OG1_fraction': fracs['hb_N2H_T199OG1']['fraction'],
        'hb_T199NH_Osulfonyl_fraction': fracs['hb_T199NH_Osulfonyl']['fraction'],
        'hb_T199OG1H_E106OE_fraction': fracs['hb_T199OG1H_E106OE']['fraction'],
        'd_Zn_N2_mean_A': stats['d_Zn_N2_A']['mean'], 'd_Zn_N2_sd_A': stats['d_Zn_N2_A']['sd'],
        'd_Zn_H94NE2_mean_A': stats['d_Zn_H94NE2_A']['mean'], 'd_Zn_H96NE2_mean_A': stats['d_Zn_H96NE2_A']['mean'],
        'd_Zn_H119ND1_mean_A': stats['d_Zn_H119ND1_A']['mean'],
        'ang_N_Zn_N_smallest_mean_deg': r6(min(ang_means)), 'ang_N_Zn_N_largest_mean_deg': r6(max(ang_means)),
        'ang_N_Zn_N_largest_abs_deviation_deg': r6(max(abs(stats[c]['mean'] - ang_ref[c]) for c in ZN_ANGLES)),
        'd_Zn_Osulfonyl_min_mean_A': stats['d_Zn_Osulfonyl_min_A']['mean'],
        'fifth_sulfonylO_fraction': fracs['fifth_sulfonylO']['fraction'],
        'd_Zn_Owater_min_mean_A': stats['d_Zn_Owater_min_A']['mean'],
        'd_Zn_Owater_min_min_A': stats['d_Zn_Owater_min_A']['min'],
        'fifth_any_fraction': fracs['fifth_any']['fraction'],
        'n_water_N1_3p5A_mean': stats['n_water_N1_3p5A']['mean'],
        'tau_fold_mean_deg': stats['tau_fold_deg']['mean'], 'tau_fold_sd_deg': stats['tau_fold_deg']['sd'],
        'ang_C4_S_N2_mean_deg': stats['ang_C4_S_N2_deg']['mean'], 'ang_C4_S_N2_sd_deg': stats['ang_C4_S_N2_deg']['sd'],
        'd_S_N2_mean_A': stats['d_S_N2_A']['mean'], 'd_S_N2_sd_A': stats['d_S_N2_A']['sd'],
        'ca_rmsd_mean_A': stats['ca_rmsd_A']['mean'], 'ca_rmsd_mean_last_half_A': stats['ca_rmsd_A']['mean_last_half'],
        'rmsf_ca_max_A': r6(rmsf[kmax]), 'rmsf_ca_max_residue': int(S['ca_resnum'][kmax]),
        'gate_GV_passed': gv['passed'], 'n_warnings': None,
        # ---- version 1.2 (sections 9A to 9E); reported only
        'n2_fraction_positive': n2['fraction_positive'], 'n2_fraction_positive_last20ns': n2['fraction_positive_window'],
        'n2_n_lasting_sign_changes': n2['n_lasting_sign_changes'],
        'n2_first_lasting_sign_change_ns': n2['first_lasting_sign_change_ns'],
        'hb_N2H_T199OG1_fraction_before_first_passage': hb_split['fraction_before'],
        'hb_N2H_T199OG1_fraction_after_first_passage': hb_split['fraction_after'],
        'd_Zn_E106OE1_mean_A': stats['d_Zn_E106OE1_A']['mean'], 'd_Zn_E106OE1_min_A': stats['d_Zn_E106OE1_A']['min'],
        'fifth_E106OE1_fraction': fracs['fifth_E106OE1']['fraction'],
        'd_Zn_E106OE2_mean_A': stats['d_Zn_E106OE2_A']['mean'], 'd_Zn_E106OE2_min_A': stats['d_Zn_E106OE2_A']['min'],
        'fifth_E106OE2_fraction': fracs['fifth_E106OE2']['fraction'],
        'n_water_ring_wide_3p5A_mean': stats['n_water_ring_wide_3p5A']['mean'],
        'd_water_ring_wide_min_min_A': stats['d_water_ring_wide_min_A']['min'],
        'state_first_frame': start_prod['state_first_frame'], 'state_first_frame_differs_from_start': start_prod['state_differs'],
        'n2_sign_first_frame': start_prod['n2_sign_first_frame'],
        'equilibration_frames_solute': equil.get('n_frames_solute', 0),
        # ---- version 1.3 (section 9D): the count the driver acts on; reported only
        'n_water_ring_place_mean': stats['n_water_ring_place']['mean'],
        'd_water_ring_place_min_min_A': stats['d_water_ring_place_min_A']['min'],
    }

    # ---- files -----------------------------------------------------------------------------------------
    scols = [c for c, _ in solute_columns]
    sfmt = [FMT[k] for _, k in solute_columns]
    write_csv(out_dir / 'timeseries_solute.csv', ['frame', 'segment', 'time_ns'] + scols,
              ([i, int(segid[i]), '%.3f' % time_ns[i]] + [f % q[c][i] for c, f in zip(scols, sfmt)] for i in range(n)))
    fcols = [c for c, _ in full_columns]
    ffmt = [FMT[k] for _, k in full_columns]
    write_csv(out_dir / 'timeseries_full.csv', ['frame', 'segment', 'time_ns'] + fcols,
              ([i, int(fsegid[i]), '%.3f' % ftime_ns[i]] + [f % fq[c][i] for c, f in zip(fcols, ffmt)] for i in range(nf)))
    write_csv(out_dir / 'rmsf_ca.csv', ['residue_number_deposited', 'residue_name', 'rmsf_A'],
              ([S['ca_resnum'][k], S['ca_resname'][k], '%.4f' % rmsf[k]] for k in range(len(rmsf))))

    inputs = {'selections.json': sha256(sysdir / 'selections.json'), 'reference.pdb': sha256(sysdir / 'reference.pdb'),
              'solute.pdb': sha256(run_dir / 'solute.pdb')}
    if (run_dir / 'full.pdb').is_file():
        inputs['full.pdb'] = sha256(run_dir / 'full.pdb')
    if R['manifest_path']:
        inputs['manifest.json'] = sha256(R['manifest_path'])
    if not args.no_hash_trajectories:
        for _, p in ssegs + fsegs + equilibration_files(run_dir, 'solute') + equilibration_files(run_dir, 'full'):
            inputs[p.name] = sha256(p)
    row['n_warnings'] = len(warnings_list)
    summary = {
        'spec': f'md/analysis/SPEC.md, version {SPEC_VERSION}',
        'script': {'name': Path(__file__).name, 'sha256': sha256(__file__)},
        'software': software_versions(),
        'random_numbers': 'none',
        'units': 'distances in Angstrom, angles in degrees, times in ns, intervals in ps',
        'run': {'run_dir': str(run_dir), 'label': R['label'], 'model': R['model'], 'start': R['start'],
                'replica': R['replica'], 'role': R['role'], 'system_dir': str(sysdir), 'planned_ns': r6(R['planned_ns']),
                'seen_ns': r6(seen_ns), 'complete': complete,
                'solute_interval_ps': r6(R['solute_interval_ps']), 'full_interval_ps': r6(R['full_interval_ps']),
                'n_frames_solute': int(n), 'n_frames_full': int(nf), 'seen_full_ns': r6(nf * dtf_ns),
                'n_atoms_solute': S['n_solute'], 'n_CA': len(S['ca']),
                'segments_solute': sinfo, 'segments_full': finfo,
                'header_check_solute': shead, 'header_check_full': fhead, 'segment_overlap_check': overlap,
                'manifest_segment_check': mcheck, 'sources': R['sources'],
                'system_files_match_manifest': syscheck,
                'imaging_solute': imaging, 'full_trajectory': fmeta,
                'largest_neighbour_CA_distance_A': r6(ca_gap), 'consistency_full_vs_solute': cons},
        'inputs_sha256': inputs,
        'constants': CONST,
        'stats': stats,
        'fractions': fracs,
        'lig_rmsd_labelling_counts': labels,
        'orientation': orient,
        'asn67_ND2': asn,
        'n2_configuration': n2,
        'hbond_first_passage': hb_split,
        'water_ring_place': place,
        'start_of_production': start_prod,
        'equilibration': equil,
        'validation': val,
        'gate_GV': gv,
        'rmsf': {'file': 'rmsf_ca.csv', 'n_residues': int(len(rmsf)), 'mean_A': r6(rmsf.mean()),
                 'max_A': r6(rmsf[kmax]), 'max_residue_number_deposited': int(S['ca_resnum'][kmax])},
        'warnings': warnings_list,
        'table_row': row,
    }
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if not args.no_figures:
        figures_run(out_dir, R, time_ns, q, ftime_ns, fq, S['ca_resnum'], rmsf)
        figures_v12(out_dir, R, time_ns, q, ftime_ns, fq, eq_series, place)
    print(f'{R["label"]}: {seen_ns:.6g} ns seen of {R["planned_ns"]:g} planned ({n} solute frames, {nf} full '
          f'frames){"" if complete else " -- NOT FINISHED"}; verdict {orient["verdict"]}'
          f'{" (provisional)" if orient["provisional"] else ""}; ligand RMSD {stats["lig_rmsd_A"]["mean"]:.3f} A; '
          + (f'VALIDATION RUN, gate G-V {gv["verdict"]}{" (provisional)" if gv["provisional"] else ""}; '
             if R['role'] == 'validation' else '')
          + f'{len(warnings_list)} warning(s); written to {out_dir}')
    return summary


# ----------------------------------------------------------------------------------------------------------
# Figures (Matplotlib, Agg).  One y axis per panel; colours in a fixed order, checked for colour-blind
# separation; identity is also given by the legend, and every plotted number is in the CSV files.
# ----------------------------------------------------------------------------------------------------------
PAL = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300')
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'text.color': INK,
                         'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6, 'axes.grid': True,
                         'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                         'axes.spines.right': False, 'legend.frameon': False, 'legend.fontsize': 7,
                         'axes.titlesize': 8.5, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
                         'figure.dpi': 150, 'savefig.dpi': 150, 'figure.facecolor': 'white',
                         'axes.facecolor': 'white', 'savefig.facecolor': 'white'})
    return plt


def _trace(ax, t, y, colour, label=None, win=100):
    """Every frame as a thin line; the running mean over `win` frames (1 ns) as the 2-point line."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    if len(y) == 0:
        return
    if len(y) >= 2 * win:
        ax.plot(t, y, color=colour, lw=0.4, alpha=0.35)
        k = np.ones(win) / win
        ax.plot(np.convolve(t, k, 'valid'), np.convolve(y, k, 'valid'), color=colour, lw=1.6, label=label)
    else:
        ax.plot(t, y, color=colour, lw=1.2, label=label)


def _ref(ax, y, text):
    ax.axhline(y, color=INK2, lw=0.7, zorder=1)
    ax.annotate(text, xy=(1.0, y), xycoords=('axes fraction', 'data'), xytext=(-2, 2), textcoords='offset points',
                ha='right', va='bottom', fontsize=6.5, color=INK2)


def _legend(ax, **kw):
    """Legend to the right of the panel, outside the data."""
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper left', bbox_to_anchor=(1.005, 1.0), borderaxespad=0.0, handlelength=1.4, **kw)


def _save(fig, path):
    fig.savefig(path, metadata={'Software': None})
    import matplotlib.pyplot as plt
    plt.close(fig)


def figures_run(out, R, t, q, tf, fq, resnum, rmsf):
    plt = _plt()
    head = f'{R["label"]} (start {R["start"]}), {len(t) * R["solute_interval_ps"] / 1000:.2f} ns'
    wf = max(1, int(round(1000.0 / R['full_interval_ps'])))

    fig, ax = plt.subplots(2, 1, figsize=(8.6, 4.6), sharex=True, gridspec_kw={'height_ratios': [3, 1]},
                           constrained_layout=True)
    _trace(ax[0], t, q['d1_N1_F131CZ_A'], PAL[0], 'd1: aniline N to Phe131 CZ')
    _trace(ax[0], t, q['d2_N1_N67OD1_A'], PAL[1], 'd2: aniline N to Asn67 OD1')
    if 'd_N1_N67ND2_A' in q:
        _trace(ax[0], t, q['d_N1_N67ND2_A'], PAL[2], 'aniline N to Asn67 ND2 (not used for the state)')
    _ref(ax[0], CONST['midpoint_d1_N1_F131CZ_A'], 'midpoint of d1, 6.7 Å')
    _ref(ax[0], CONST['midpoint_d2_N1_N67OD1_A'], 'midpoint of d2, 5.2 Å')
    ax[0].set_ylabel('distance (Å)')
    ax[0].set_title(f'Orientation of the ring: {head}')
    _legend(ax[0])
    y = np.select([q['state'] == 'C', q['state'] == 'I'], [2, 1], 0)
    ax[1].plot(t, y, ls='none', marker='|', ms=5, mew=0.5, color=INK)
    ax[1].set_yticks([0, 1, 2], ['rotated-like', 'intermediate', 'crystal-like'])
    ax[1].set_ylim(-0.5, 2.5)
    ax[1].set_xlabel('time (ns)')
    ax[1].set_title('State of each frame')
    _save(fig, out / 'fig_orientation.png')

    fig, ax = plt.subplots(1, 1, figsize=(7.2, 2.8), constrained_layout=True)
    _trace(ax, t, q['lig_rmsd_A'], PAL[0])
    _ref(ax, CONST['rule1_rmsd_A'], 'rule 1, 2.0 Å')
    _ref(ax, CONST['gate_GV']['lig_rmsd_mean_max_A'], 'gate G-V, 1.0 Å')
    ax.set_ylim(bottom=0)
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title(f'Ligand heavy-atom RMSD to the crystal pose (fit on CA, symmetry-corrected): {head}')
    _save(fig, out / 'fig_ligand.png')

    fig, ax = plt.subplots(3, 1, figsize=(8.8, 7.2), sharex=True, constrained_layout=True)
    for i, (c, lab) in enumerate((('d_Zn_N2_A', 'N2 (sulfonamide)'), ('d_Zn_H94NE2_A', 'His94 NE2'),
                                  ('d_Zn_H96NE2_A', 'His96 NE2'), ('d_Zn_H119ND1_A', 'His119 ND1'))):
        _trace(ax[0], t, q[c], PAL[i], lab)
    ax[0].set_ylabel('distance (Å)')
    ax[0].set_title(f'Zinc to its four nitrogen atoms: {head}')
    _legend(ax[0])
    nice = {'N2': 'N2', 'H94NE2': 'His94 NE2', 'H96NE2': 'His96 NE2', 'H119ND1': 'His119 ND1'}
    for i, c in enumerate(ZN_ANGLES):
        a, b = c[4:-4].split('_Zn_')
        _trace(ax[1], t, q[c], PAL[i], f'{nice[a]}–Zn–{nice[b]}')
    ax[1].set_ylabel('angle (°)')
    ax[1].set_title('The six N–Zn–N angles')
    _legend(ax[1])
    _trace(ax[2], t, q['d_Zn_Osulfonyl_min_A'], PAL[0], 'sulfonyl O (every 10 ps)')
    if len(tf):
        _trace(ax[2], tf, fq['d_Zn_Owater_min_A'], PAL[1], 'water O (every 100 ps)', win=wf)
        _trace(ax[2], tf, fq['d_Zn_Oprotein_min_A'], PAL[2], 'protein O (every 100 ps)', win=wf)
    _ref(ax[2], CONST['fifth_ligand_distance_max_A'], 'fifth ligand, 2.6 Å')
    ax[2].set_ylabel('distance (Å)')
    ax[2].set_xlabel('time (ns)')
    ax[2].set_title('Nearest oxygen atom to the zinc')
    _legend(ax[2])
    _save(fig, out / 'fig_zinc.png')

    fig, ax = plt.subplots(3, 1, figsize=(8.8, 7.2), sharex=True, constrained_layout=True)
    trio = (('d_N2_T199OG1_A', 'hb_N2H_T199OG1', 'N2–H···OG1 Thr199'),
            ('d_T199N_Osulfonyl_min_A', 'hb_T199NH_Osulfonyl', 'Thr199 N–H···O sulfonyl'),
            ('d_T199OG1_E106OE_min_A', 'hb_T199OG1H_E106OE', 'Thr199 OG1–H···OE Glu106'))
    for i, (c, _, lab) in enumerate(trio):
        _trace(ax[0], t, q[c], PAL[i], lab)
    _ref(ax[0], CONST['hbond_distance_max_A'], 'criterion, 3.5 Å')
    ax[0].set_ylabel('donor···acceptor (Å)')
    ax[0].set_title(f'Hydrogen bonds at the site: {head}')
    _legend(ax[0])
    win = 100
    for i, (_, c, lab) in enumerate(trio):
        f = np.asarray(q[c], float)
        if len(f) >= win:
            k = np.ones(win) / win
            ax[1].plot(np.convolve(t, k, 'valid'), 100 * np.convolve(f, k, 'valid'), color=PAL[i], lw=1.6, label=lab)
    ax[1].set_ylim(-3, 103)
    ax[1].set_ylabel('present (% of frames)')
    ax[1].set_title('Hydrogen bond present (distance and angle), running mean over 1 ns')
    _legend(ax[1])
    if len(tf):
        ax[2].plot(tf, fq['n_water_N1_3p5A'], color=PAL[0], lw=0.8, drawstyle='steps-mid')
    ax[2].set_ylabel('number')
    ax[2].set_xlabel('time (ns)')
    ax[2].set_title('Water oxygen atoms within 3.5 Å of the aniline N (every 100 ps)')
    _save(fig, out / 'fig_contacts.png')

    fig, ax = plt.subplots(4, 1, figsize=(8.2, 8.6), sharex=True, constrained_layout=True)
    _trace(ax[0], t, q['tau_fold_deg'], PAL[0])
    _ref(ax[0], CONST['gate_GV']['tau_fold_reference_deg'], 'crystal, 132.1°')
    ax[0].set_ylabel('τ_fold (°)')
    ax[0].set_title(f'τ_fold = max(|C3–C4–S–N2|, |C5–C4–S–N2|): {head}')
    ax[1].plot(t, q['tors_C3_C4_S_N2_deg'], ls='none', marker='.', ms=1.5, color=PAL[0], label='C3–C4–S–N2')
    ax[1].plot(t, q['tors_C5_C4_S_N2_deg'], ls='none', marker='.', ms=1.5, color=PAL[1], label='C5–C4–S–N2')
    ax[1].set_ylim(-180, 180)
    ax[1].set_yticks([-180, -90, 0, 90, 180])
    ax[1].set_ylabel('torsion (°)')
    ax[1].set_title('The two signed torsions')
    _legend(ax[1], markerscale=5)
    _trace(ax[2], t, q['ang_C4_S_N2_deg'], PAL[0])
    _ref(ax[2], CONST['gate_GV']['ang_C4_S_N2_reference_deg'], 'crystal, 108.8°')
    ax[2].set_ylabel('angle (°)')
    ax[2].set_title('Angle C4–S–N2')
    _trace(ax[3], t, q['d_S_N2_A'], PAL[0])
    ax[3].set_ylabel('bond (Å)')
    ax[3].set_xlabel('time (ns)')
    ax[3].set_title('Bond S–N2')
    _save(fig, out / 'fig_internal.png')

    fig, ax = plt.subplots(2, 1, figsize=(7.2, 4.8), constrained_layout=True)
    _trace(ax[0], t, q['ca_rmsd_A'], PAL[0])
    ax[0].set_ylim(bottom=0)
    ax[0].set_xlabel('time (ns)')
    ax[0].set_ylabel('RMSD (Å)')
    ax[0].set_title(f'CA RMSD to the crystal structure after the fit: {head}')
    ax[1].plot(resnum, rmsf, color=PAL[0], lw=1.0)
    ax[1].set_ylim(bottom=0)
    ax[1].set_xlabel('residue number (6RL9)')
    ax[1].set_ylabel('RMSF (Å)')
    ax[1].set_title('CA RMSF about the mean structure of the run')
    _save(fig, out / 'fig_protein.png')


def figures_v12(out, R, t, q, tf, fq, eq, place):
    """Versions 1.2 and 1.3: the configuration at N2, the zinc to Glu106, water in the place of the crystal ring
    (fig_n2_site.png); the equilibration frames, if there are any (fig_equilibration.png)."""
    plt = _plt()
    head = f'{R["label"]} (start {R["start"]}), {len(t) * R["solute_interval_ps"] / 1000:.2f} ns'
    wf = max(1, int(round(1000.0 / R['full_interval_ps'])))
    fig, ax = plt.subplots(4, 1, figsize=(8.8, 8.8), sharex=True, constrained_layout=True)
    _trace(ax[0], t, q['n2_volume'], PAL[0])
    _ref(ax[0], 0.0, 'planar nitrogen, 0')
    ax[0].set_ylim(-1.05, 1.05)
    ax[0].set_ylabel('signed volume')
    ax[0].set_title(f'Configuration at N2, e(Zn)·(e(H)×e(S)): {head}')
    _trace(ax[1], t, q['d_Zn_E106OE1_A'], PAL[0], 'Glu106 OE1')
    _trace(ax[1], t, q['d_Zn_E106OE2_A'], PAL[1], 'Glu106 OE2')
    _ref(ax[1], CONST['glu106_distance_max_A'], 'fifth ligand, 2.6 Å')
    ax[1].set_ylabel('distance (Å)')
    ax[1].set_title('Zinc to the carboxylate oxygen atoms of Glu106 (every 10 ps)')
    _legend(ax[1])
    cut, names = place['water_cutoff_A'], ', '.join(place['site_atoms'])
    if len(tf):
        ax[2].plot(tf, fq['n_water_ring_place'], color=PAL[0], lw=1.0, drawstyle='steps-mid',
                   label=f'in the place of the ring: within {cut:g} Å of the crystal positions of {names}')
        ax[2].plot(tf, fq['n_water_ring_wide_3p5A'], color=PAL[1], lw=0.6, drawstyle='steps-mid',
                   label='wider count: within 3.5 Å of the crystal positions of C1 to C6 and N1')
        _trace(ax[3], tf, fq['d_water_ring_place_min_A'], PAL[0], win=wf)
        _legend(ax[2])
    ax[2].set_ylabel('number')
    ax[2].set_title('Water oxygen atoms at the crystal positions of the ring (every 100 ps)')
    _ref(ax[3], cut, f'cut-off of the count, {cut:g} Å')
    ax[3].set_ylabel('distance (Å)')
    ax[3].set_xlabel('time (ns)')
    ax[3].set_title(f'Nearest water oxygen atom to the crystal positions of {names}')
    _save(fig, out / 'fig_n2_site.png')

    path = out / 'fig_equilibration.png'
    if eq is None:
        if path.is_file():
            path.unlink()
        return
    te, qe, tfe, fqe = eq['t'], eq['q'], eq['tf'], eq['fq']
    fig, ax = plt.subplots(4, 1, figsize=(8.8, 8.8), sharex=True, constrained_layout=True)
    if len(te):
        ax[0].plot(te, qe['d1_N1_F131CZ_A'], color=PAL[0], lw=1.2, label='d1: aniline N to Phe131 CZ')
        ax[0].plot(te, qe['d2_N1_N67OD1_A'], color=PAL[1], lw=1.2, label='d2: aniline N to Asn67 OD1')
        ax[1].plot(te, qe['lig_rmsd_A'], color=PAL[0], lw=1.2)
        ax[2].plot(te, qe['n2_volume'], color=PAL[0], lw=1.2)
    _ref(ax[0], CONST['midpoint_d1_N1_F131CZ_A'], 'midpoint of d1, 6.7 Å')
    _ref(ax[0], CONST['midpoint_d2_N1_N67OD1_A'], 'midpoint of d2, 5.2 Å')
    ax[0].set_ylabel('distance (Å)')
    ax[0].set_title(f'Equilibration (not part of any statistic), time from the start of the production: {R["label"]}')
    _legend(ax[0])
    ax[1].set_ylim(bottom=0)
    ax[1].set_ylabel('RMSD (Å)')
    ax[1].set_title('Ligand heavy-atom RMSD to the crystal pose')
    _ref(ax[2], 0.0, 'planar nitrogen, 0')
    ax[2].set_ylim(-1.05, 1.05)
    ax[2].set_ylabel('signed volume')
    ax[2].set_title('Configuration at N2')
    if fqe is not None and len(tfe):
        ax[3].plot(tfe, fqe['n_water_ring_place'], color=PAL[0], lw=1.0, marker='o', ms=3,
                   label=f'in the place of the ring ({place["water_cutoff_A"]:g} Å, {", ".join(place["site_atoms"])})')
        ax[3].plot(tfe, fqe['n_water_ring_wide_3p5A'], color=PAL[1], lw=0.6, marker='o', ms=2,
                   label='wider count (3.5 Å, C1 to C6 and N1)')
        _legend(ax[3])
    ax[3].set_ylabel('number')
    ax[3].set_xlabel('time (ns)')
    ax[3].set_title('Water oxygen atoms at the crystal positions of the ring')
    _save(fig, path)


def figures_campaign(out, runs):
    """One panel per run (rows: starts, columns: runs of the start); nothing is averaged over runs."""
    plt = _plt()
    by = {}
    for r in runs:
        by.setdefault((r['row']['model'], r['row']['start']), []).append(r)
    models = sorted({k[0] for k in by})
    for model in models:
        starts = [s for s in ('X', 'N', 'R') if (model, s) in by]
        ncol = max(len(by[(model, s)]) for s in starts)
        for name, what in (('campaign_orientation', 'o'), ('campaign_ligand_rmsd', 'r')):
            fig, axes = plt.subplots(len(starts), ncol, figsize=(2.6 * ncol + 0.6, 1.9 * len(starts) + 0.7),
                                     sharex=True, sharey=True, squeeze=False, constrained_layout=True)
            for i, s in enumerate(starts):
                for j in range(ncol):
                    ax = axes[i][j]
                    if j >= len(by[(model, s)]):
                        ax.set_visible(False)
                        continue
                    r = by[(model, s)][j]
                    ts = r['timeseries']
                    if what == 'o':
                        _trace(ax, ts['time_ns'], ts['d1_N1_F131CZ_A'], PAL[0], 'd1: aniline N to Phe131 CZ')
                        _trace(ax, ts['time_ns'], ts['d2_N1_N67OD1_A'], PAL[1], 'd2: aniline N to Asn67 OD1')
                        ax.axhline(CONST['midpoint_d1_N1_F131CZ_A'], color=INK2, lw=0.6)
                        ax.axhline(CONST['midpoint_d2_N1_N67OD1_A'], color=INK2, lw=0.6)
                    else:
                        _trace(ax, ts['time_ns'], ts['lig_rmsd_A'], PAL[0])
                        ax.axhline(CONST['rule1_rmsd_A'], color=INK2, lw=0.6)
                    v = r['row']['verdict'] + (' (provisional)' if r['row']['provisional'] else '')
                    ax.set_title(f'{r["row"]["run"]}: {v}', fontsize=7.5)
                    if i == len(starts) - 1:
                        ax.set_xlabel('time (ns)')
                    if j == 0:
                        ax.set_ylabel('distance (Å)' if what == 'o' else 'RMSD (Å)')
            if what == 'o':
                h, lab = axes[0][0].get_legend_handles_labels()
                fig.legend(h, lab, loc='outside lower center', ncol=2)
                fig.suptitle('Orientation distances of every run; lines at the midpoints 6.7 Å (d1) and 5.2 Å (d2)',
                             fontsize=8.5, fontweight='bold', x=0.01, ha='left')
            else:
                fig.suptitle('Ligand heavy-atom RMSD to the crystal pose of every run; line at 2.0 Å (rule 1)',
                             fontsize=8.5, fontweight='bold', x=0.01, ha='left')
            suffix = '' if len(models) == 1 else f'_{model}'
            _save(fig, out / f'{name}{suffix}.png')


# ----------------------------------------------------------------------------------------------------------
# Campaign (SPEC section 12)
# ----------------------------------------------------------------------------------------------------------
def find_runs(root):
    root = Path(root).resolve()
    dirs = sorted({p.parent for p in root.glob('**/solute_seg*.dcd') if re.fullmatch(r'solute_seg\d+\.dcd', p.name)})
    return dirs


def rules(rows):
    """Rules 1, 3 and 4 of section 5 of md/PLAN.md for the runs of one zinc model (SPEC section 12)."""
    by = {s: [r for r in rows if r['start'] == s] for s in ('X', 'N', 'R')}
    all_nine = all(len(by[s]) == 3 for s in by) and all(r['complete'] for r in rows)
    x_left = sum(1 for r in by['X'] if r['rule1_left_crystal'])
    rule1_failed = x_left >= 2
    conv = {s: sum(1 for r in by[s] if r['verdict'] == 'CONVERTED') for s in by}
    stay = {s: sum(1 for r in by[s] if r['verdict'] == 'STAYED') for s in by}
    xn = by['X'] + by['N']
    if rule1_failed:
        outcome = 'validation_failed'
    elif conv['R'] >= 2 and xn and all(r['verdict'] == 'STAYED' for r in xn):
        outcome = 'rotated_converts'
    elif rows and all(r['verdict'] == 'STAYED' for r in rows):
        outcome = 'both_persist'
    else:
        outcome = 'mixed'
    return {'n_runs': len(rows), 'n_complete': sum(1 for r in rows if r['complete']),
            'runs_per_start': {s: len(by[s]) for s in by}, 'all_nine_complete': all_nine,
            'final': all_nine,
            'rule1': {'n_X_runs': len(by['X']), 'n_X_left_crystal': x_left, 'failed': rule1_failed},
            'n_converted': conv, 'n_stayed': stay,
            'n_mixed': {s: len(by[s]) - conv[s] - stay[s] for s in by},
            'outcome': outcome,
            'rule4_sensitivity_run_required': outcome == 'rotated_converts'}


OUTCOME_TEXT = {
    'validation_failed': 'rule 1: the crystal-start runs left the crystal pose in two or three runs; the simulation '
                         'is a failed validation and supports nothing',
    'rotated_converts': 'rule 3: the rotated starts converted in two or three runs and the crystal and near-native '
                        'starts stayed; the rotated orientation is not stable outside the scoring function',
    'both_persist': 'rule 3: no run left its starting orientation; both orientations persisted and the simulation '
                    'does not decide',
    'mixed': 'rule 3: the runs disagree; they are reported run by run and nothing is averaged over them',
}


def campaign(root, args):
    root = Path(root).resolve()
    out = Path(args.out).resolve() if args.out else root / 'analysis_campaign'
    out.mkdir(parents=True, exist_ok=True)
    runs, failed, skipped = [], [], []
    spec_line = f'md/analysis/SPEC.md, version {SPEC_VERSION}'
    for d in find_runs(root):
        if not re.fullmatch(r'[XNR][1-9]\d*', d.name):         # production runs are <start><replica>, replica from 1;
            skipped.append(str(d))                            # the validation run X0 is analysed on its own
            print(f'{d}: not a production run folder (name), left out of the campaign')
            continue
        try:
            sj = d / 'analysis' / 'summary.json'
            s = None
            if args.reuse and sj.is_file():
                s = json.loads(sj.read_text())
                if s.get('spec') == spec_line:
                    print(f'{s["run"]["label"]}: summary.json reused')
                else:                                         # written after another version of the specification
                    s = None
            if s is None:
                s = analyse_run(d, args)
        except AnalysisError as e:
            failed.append({'run_dir': str(d), 'error': str(e)})
            print(f'{d}: FAILED: {e}', file=sys.stderr)
            continue
        if s['run']['role'] != 'production':                   # a validation run by its manifest under another name
            skipped.append(str(d))
            print(f'{d}: {s["run"]["role"]} run by its manifest, left out of the campaign')
            continue
        ts = {}
        with open(d / 'analysis' / 'timeseries_solute.csv') as f:
            rd = csv.DictReader(f)
            keep = ('time_ns', 'd1_N1_F131CZ_A', 'd2_N1_N67OD1_A', 'lig_rmsd_A')
            cols = {k: [] for k in keep}
            for line in rd:
                for k in keep:
                    cols[k].append(float(line[k]))
            ts = {k: np.array(v) for k, v in cols.items()}
        runs.append({'row': s['table_row'], 'summary_file': str(sj), 'summary_sha256': sha256(sj), 'timeseries': ts})
    order = {'X': 0, 'N': 1, 'R': 2}
    runs.sort(key=lambda r: (r['row']['model'], order[r['row']['start']], str(r['row']['replica']), r['row']['run']))
    rows = [r['row'] for r in runs]
    keys = list(rows[0].keys()) if rows else []
    write_csv(out / 'campaign_table.csv', keys, ([('' if r[k] is None else r[k]) for k in keys] for r in rows))
    per_model = {m: rules([r for r in rows if r['model'] == m]) for m in sorted({r['model'] for r in rows})}

    def f(v, nd):
        return '–' if v is None else f'{v:.{nd}f}'

    md = ['# MD campaign: one line per run', '',
          f'Written by `md/analysis/analyse_md.py` after `md/analysis/SPEC.md`, version {SPEC_VERSION}. Nothing is '
          'averaged over runs. All columns with full precision: `campaign_table.csv`; every number: `campaign.json` '
          'and the `summary.json` of each run.', '']
    for m, ru in per_model.items():
        md += [f'**Zinc model `{m}`**: {ru["n_runs"]} runs found, {ru["n_complete"]} complete '
               f'(X {ru["runs_per_start"]["X"]}, N {ru["runs_per_start"]["N"]}, R {ru["runs_per_start"]["R"]}). '
               f'Rule 1: {ru["rule1"]["n_X_left_crystal"]} of {ru["rule1"]["n_X_runs"]} crystal-start runs left the '
               f'crystal pose ({"FAILED" if ru["rule1"]["failed"] else "not failed"}). '
               f'Converted: X {ru["n_converted"]["X"]}, N {ru["n_converted"]["N"]}, R {ru["n_converted"]["R"]}. '
               f'Outcome{"" if ru["final"] else " (NOT FINAL: the nine runs are not all present and complete)"}: '
               f'`{ru["outcome"]}` – {OUTCOME_TEXT[ru["outcome"]]}.', '']
    md += ['| Run | ns seen / planned | Verdict | Start state, last 20 ns (%) | Opposite state, last 20 ns (%) | '
           'First passage (ns) | Transitions | Ligand RMSD, mean ± SD (Å) | Ligand RMSD, last half (Å) | '
           'Left crystal (rule 1) | N2–H···OG1 Thr199 (% of frames) | Zn–N2 (Å) | τ_fold, mean ± SD (°) | '
           'C4–S–N2, mean ± SD (°) | Fifth ligand (% of frames, 100 ps) | CA RMSD (Å) |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        md.append(
            f'| {r["run"]} | {f(r["seen_ns"], 2)} / {r["planned_ns"]:g} | '
            f'{r["verdict"]}{" (provisional)" if r["provisional"] else ""} | '
            f'{f(100 * r["frac_start_state_last20ns"], 1)} | {f(100 * r["frac_opposite_state_last20ns"], 1)} | '
            f'{f(r["first_passage_ns"], 2)} | {r["n_transitions"]} | {f(r["lig_rmsd_mean_A"], 2)} ± {f(r["lig_rmsd_sd_A"], 2)} | '
            f'{f(r["lig_rmsd_mean_last_half_A"], 2)} | {"yes" if r["rule1_left_crystal"] else "no"} | '
            f'{f(100 * r["hb_N2H_T199OG1_fraction"], 1)} | {f(r["d_Zn_N2_mean_A"], 3)} | '
            f'{f(r["tau_fold_mean_deg"], 1)} ± {f(r["tau_fold_sd_deg"], 1)} | '
            f'{f(r["ang_C4_S_N2_mean_deg"], 1)} ± {f(r["ang_C4_S_N2_sd_deg"], 1)} | '
            f'{f(None if r["fifth_any_fraction"] is None else 100 * r["fifth_any_fraction"], 1)} | '
            f'{f(r["ca_rmsd_mean_A"], 2)} |')
    if failed:
        md += ['', 'Run folders that could not be analysed:', ''] + [f'- `{x["run_dir"]}`: {x["error"]}' for x in failed]
    if skipped:
        md += ['', 'Folders with a trajectory that are not production runs (by their name, or a validation run by '
               'its manifest); left out, the validation run X0 is analysed on its own:', '']
        md += [f'- `{x}`' for x in skipped]
    (out / 'campaign_table.md').write_text('\n'.join(md) + '\n')
    cj = {'spec': spec_line,
          'script': {'name': Path(__file__).name, 'sha256': sha256(__file__)},
          'software': software_versions(), 'campaign_dir': str(root),
          'rules_per_model': per_model, 'outcome_text': OUTCOME_TEXT,
          'runs': [{'table_row': r['row'], 'summary_file': r['summary_file'], 'summary_sha256': r['summary_sha256']}
                   for r in runs],
          'failed': failed, 'skipped': skipped}
    (out / 'campaign.json').write_text(json.dumps(cj, indent=1) + '\n')
    if runs and not args.no_figures:
        figures_campaign(out, runs)
    print(f'campaign: {len(rows)} runs, {len(failed)} failed, {len(skipped)} left out; written to {out}')
    for m, ru in per_model.items():
        print(f'  model {m}: outcome {ru["outcome"]}{"" if ru["final"] else " (not final)"}')
    return cj


# ----------------------------------------------------------------------------------------------------------
# The definitions applied to the three starts (SPEC section 3, check against md/system/starts.csv)
# ----------------------------------------------------------------------------------------------------------
def check_starts(args):
    sysd = Path(__file__).resolve().parent.parent / 'system'
    files = {'protein': sysd / 'tleap' / 'protein_H.pdb', 'protein_heavy': sysd / 'protein_md.pdb',
             'zinc': sysd / 'zinc.pdb', 'starts': sysd / 'starts.csv',
             **{f'ligand_{s}': sysd / f'ligand_{s}.pdb' for s in 'XNR'},
             **{f'waters_{s}': sysd / f'waters_{s}.pdb' for s in 'XNR'}}
    pn, _, pi, px = read_pdb(files['protein'])
    hn, _, hi, hx = read_pdb(files['protein_heavy'])

    def P(num, name):
        hit = [k for k in range(len(pn)) if pn[k] == name and pi[k] == str(num)]
        assert len(hit) == 1, (num, name)
        return px[hit[0]]
    ca = np.array([hx[k] for k in range(len(hn)) if hn[k] == 'CA'])
    site = {'his94_NE2': P(94, 'NE2'), 'his96_NE2': P(96, 'NE2'), 'his119_ND1': P(119, 'ND1'),
            'thr199_OG1': P(199, 'OG1'), 'thr199_HG1': P(199, 'HG1'), 'thr199_N': P(199, 'N'),
            'thr199_H': P(199, 'H'), 'glu106_OE1': P(106, 'OE1'), 'glu106_OE2': P(106, 'OE2'),
            'phe131_CZ': P(131, 'CZ'), 'asn67_OD1': P(67, 'OD1'), 'asn67_ND2': P(67, 'ND2'),
            'zn': read_pdb(files['zinc'])[3][0]}
    lig = {}
    for s in 'XNR':
        n, _, _, x = read_pdb(files[f'ligand_{s}'])
        lig[s] = dict(zip(n, x))
    # a fixed rotation and translation of each start (protein, zinc and ligand together): the fit must remove it
    th = np.radians([35.0, -110.0, 200.0])
    rot = []
    for k, ax in enumerate(((1.0, 2.0, 3.0), (-2.0, 1.0, 0.5), (0.3, -1.0, 2.0))):
        u = np.array(ax) / np.linalg.norm(ax)
        K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
        rot.append(np.eye(3) + np.sin(th[k]) * K + (1 - np.cos(th[k])) * K @ K)
    shift = np.array([[12.0, -7.0, 30.0], [-40.0, 5.0, 2.0], [3.0, 60.0, -15.0]])
    A = {k: np.array([v @ rot[i].T + shift[i] for i in range(3)]) for k, v in site.items()}
    for n in LIG_HEAVY:
        A[n] = np.array([lig[s][n] @ rot[i].T + shift[i] for i, s in enumerate('XNR')])
    A['HN2'] = np.array([lig[s]['HN21'] @ rot[i].T + shift[i] for i, s in enumerate('XNR')])
    caf = np.array([ca @ rot[i].T + shift[i] for i in range(3)])
    q, _ = solute_quantities(A, caf, ca, {n: lig['X'][n] for n in LIG_HEAVY})
    with open(files['starts']) as f:
        ref = {r['start']: r for r in csv.DictReader(f)}
    pairs = (('rmsd_to_crystal_A', 'lig_rmsd_A', 'A'), ('Zn_N2', 'd_Zn_N2_A', 'A'), ('Zn_O_min', 'd_Zn_Osulfonyl_min_A', 'A'),
             ('Zn_O1', 'd_Zn_O1_A', 'A'), ('Zn_O2', 'd_Zn_O2_A', 'A'), ('N2_OG1', 'd_N2_T199OG1_A', 'A'),
             ('N2_H_OG1', 'ang_N2_HN2_T199OG1_deg', 'deg'), ('T199N_O_min', 'd_T199N_Osulfonyl_min_A', 'A'),
             ('C3_C4_S_N2', 'tors_C3_C4_S_N2_deg', 'deg'), ('C5_C4_S_N2', 'tors_C5_C4_S_N2_deg', 'deg'),
             ('tau_fold', 'tau_fold_deg', 'deg'), ('C4_S_N2', 'ang_C4_S_N2_deg', 'deg'), ('S_N2', 'd_S_N2_A', 'A'),
             ('S_N2_Zn', 'ang_S_N2_Zn_deg', 'deg'), ('N1_Phe131CZ', 'd1_N1_F131CZ_A', 'A'),
             ('N1_Asn67OD1', 'd2_N1_N67OD1_A', 'A'), ('N1_Asn67ND2', 'd_N1_N67ND2_A', 'A'))
    tol = {'A': 0.0015, 'deg': 0.06}        # starts.csv is rounded to 3 decimals; coordinates have 3 decimals
    rows, worst = [], {'A': 0.0, 'deg': 0.0}
    for i, s in enumerate('XNR'):
        for a, b, u in pairs:
            mine, theirs = float(q[b][i]), float(ref[s][a])
            d = abs(mine - theirs)
            worst[u] = max(worst[u], d)
            rows.append({'start': s, 'starts_csv_column': a, 'starts_csv': theirs, 'spec_column': b,
                         'spec_value': round(mine, 4), 'abs_diff': round(d, 4), 'unit': u, 'agrees': d <= tol[u]})
    extra = {s: {'lig_rmsd_labelling': LABELLINGS[int(q['lig_rmsd_labelling'][i])][0],
                 'lig_rmsd_identity_A': r6(q['lig_rmsd_identity_A'][i]),
                 'lig_rmsd_all_four_labellings_A': None,
                 'state': str(q['state'][i]), 'ca_rmsd_after_fit_A': r6(q['ca_rmsd_A'][i]),
                 'asn67_ND2_nearer_N1': int(q['asn67_ND2_nearer_N1'][i]),
                 'hb_N2H_T199OG1': int(q['hb_N2H_T199OG1'][i]), 'hb_T199NH_Osulfonyl': int(q['hb_T199NH_Osulfonyl'][i]),
                 'hb_T199OG1H_E106OE': int(q['hb_T199OG1H_E106OE'][i]),
                 'd_Zn_NHis_A': {c: r6(q[c][i]) for c in ZN_HIS},
                 'ang_N_Zn_N_deg': {c: r6(q[c][i]) for c in ZN_ANGLES},
                 'ang_T199N_H_O1_deg': r6(q['ang_T199N_H_O1_deg'][i]), 'd_T199N_O1_A': r6(q['d_T199N_O1_A'][i]),
                 'd_T199OG1_E106OE1_A': r6(q['d_T199OG1_E106OE1_A'][i]),
                 'ang_T199OG1_HG1_E106OE1_deg': r6(q['ang_T199OG1_HG1_E106OE1_deg'][i])}
             for i, s in enumerate('XNR')}
    for i, s in enumerate('XNR'):                            # the four labellings, in the 6RL9 frame, no fit
        vals = []
        for _, perm in LABELLINGS:
            vals.append(r6(math.sqrt(np.mean([((lig[s][n] - lig['X'][perm.get(n, n)]) ** 2).sum() for n in LIG_HEAVY]))))
        extra[s]['lig_rmsd_all_four_labellings_A'] = dict(zip([n for n, _ in LABELLINGS], vals))
    # version 1.2, section 9A: the signed volume at N2 of the three starts against the numbers of the independent review
    # of 6 Oct 2026 (paper/audit/research/md_setup_review_2026-10-06.json, finding 1, as md/PLAN.md quotes them)
    review = {'X': 0.72, 'N': 0.79, 'R': -0.71}
    n2rows = []
    for i, s in enumerate('XNR'):
        v = float(q['n2_volume'][i])
        n2rows.append({'start': s, 'n2_volume': round(v, 4), 'n2_sign': int(q['n2_sign'][i]),
                       'review_value': review[s], 'abs_diff': round(abs(v - review[s]), 4),
                       'agrees_with_the_review': abs(v - review[s]) <= 0.005,
                       'sign_of_start_in_the_specification': CONST['n2_sign_of_start'][s],
                       'sign_agrees': int(q['n2_sign'][i]) == CONST['n2_sign_of_start'][s]})
        extra[s].update(n2_volume=r6(v), n2_sign=int(q['n2_sign'][i]),
                        d_Zn_E106OE1_A=r6(q['d_Zn_E106OE1_A'][i]), d_Zn_E106OE2_A=r6(q['d_Zn_E106OE2_A'][i]))
    # section 9D: crystal water of each start (md/system/waters_<start>.pdb, oxygen atoms, 6RL9 frame, no periodic
    # cell) against the crystal positions of C1 to C6 (2.5 A: the defaults of the count the driver acts on) and of
    # C1 to C6 and N1 (3.5 A: the wider count of version 1.2)
    ring = {}
    for s in 'XNR':
        wn, _, wi, wx = read_pdb(files[f'waters_{s}'])
        ring[s] = {'n_crystal_water': int(len(wx))}
        for key, atoms, cut in (('place', CONST['ring_place_atoms_default'], CONST['ring_place_water_cutoff_default_A']),
                                ('wide', CONST['ring_wide_atoms'], CONST['ring_wide_water_distance_max_A'])):
            site = np.array([lig['X'][a] for a in atoms])
            d = np.sqrt(((wx[:, None, :] - site[None, :, :]) ** 2).sum(-1)).min(axis=1)
            k = np.argsort(d, kind='stable')[:3]
            ring[s].update({f'n_water_ring_{key}' + ('_3p5A' if key == 'wide' else ''): int(np.sum(d <= cut)),
                            f'd_water_ring_{key}_min_A': r6(d.min()),
                            f'nearest_three_{key}': [{'water_number_6RL9': wi[j], 'distance_A': round(float(d[j]), 3)} for j in k]})
    # the references of gate G-V (section 11) against the crystal start: each is the crystal value to one decimal
    g, gref = CONST['gate_GV'], []
    for c in ZN_ANGLES:
        gref.append({'column': c, 'crystal_start_deg': round(float(q[c][0]), 4),
                     'reference_in_the_specification_deg': g['ang_N_Zn_N_reference_6RL9_deg'][c]})
    gref.append({'column': 'ang_C4_S_N2_deg', 'crystal_start_deg': round(float(q['ang_C4_S_N2_deg'][0]), 4),
                 'reference_in_the_specification_deg': g['ang_C4_S_N2_reference_deg']})
    gref.append({'column': 'tau_fold_deg', 'crystal_start_deg': round(float(q['tau_fold_deg'][0]), 4),
                 'reference_in_the_specification_deg': g['tau_fold_reference_deg']})
    for r in gref:
        r['abs_diff_deg'] = round(abs(r['crystal_start_deg'] - r['reference_in_the_specification_deg']), 4)
        r['agrees'] = r['abs_diff_deg'] <= 0.05
    # the geometric criteria of the gate on the crystal start itself (one frame; no hydrogen bond or fifth-ligand count)
    one = {c: {'mean': r6(q[c][0]), 'sd': 0.0} for c in ('lig_rmsd_A', 'd_Zn_N2_A', 'ang_C4_S_N2_deg', 'tau_fold_deg')
           + ZN_HIS + ZN_ANGLES}
    f1 = {'hb_N2H_T199OG1': fraction([int(q['hb_N2H_T199OG1'][0])], 0),
          'fifth_sulfonylO': fraction([int(q['fifth_sulfonylO'][0])], 0), 'fifth_any': fraction([], 0)}
    gx = gate_gv(one, f1, 'production', True)
    crystal = {'d_Zn_N2_A': gx['GV3_zinc']['criteria']['d_Zn_N2_A'],
               **{c: gx['GV3_zinc']['criteria'][c] for c in ZN_HIS + ZN_ANGLES},
               **gx['GV4_ligand_internal']['criteria']}
    res = {'spec': f'md/analysis/SPEC.md, version {SPEC_VERSION}',
           'script': {'name': Path(__file__).name, 'sha256': sha256(__file__)},
           'inputs_sha256': {str(p.relative_to(sysd.parent)): sha256(p) for p in files.values()},
           'gate_references_against_the_crystal_start': {'tolerance_deg': 0.05, 'all_agree': all(r['agrees'] for r in gref),
                                                         'values': gref},
           'crystal_start_against_the_geometric_criteria_of_gate_GV': {
               'note': 'one frame, the crystal start X as built: the windows on the distances and angles of GV3 and GV4',
               'all_met': all(v['passed'] for v in crystal.values()), 'criteria': crystal},
           'note': 'each start was rotated and translated as a whole before the analysis; the fit on the CA atoms '
                   'removes this (ca_rmsd_after_fit_A)',
           'tolerance': tol, 'largest_abs_diff': {k: round(v, 4) for k, v in worst.items()},
           'n_compared': len(rows), 'n_agree': sum(r['agrees'] for r in rows),
           'all_agree': all(r['agrees'] for r in rows), 'comparison': rows, 'other_values_of_the_starts': extra,
           'n2_configuration_of_the_starts': {
               'order': 'e(Zn) . (e(HN21) x e(S)), unit vectors from N2', 'tolerance': 0.005,
               'review': 'paper/audit/research/md_setup_review_2026-10-06.json, finding 1 (two decimals, as in md/PLAN.md)',
               'all_agree': all(r['agrees_with_the_review'] and r['sign_agrees'] for r in n2rows), 'values': n2rows},
           'crystal_water_in_the_place_of_the_ring_of_the_starts': {
               'note': 'crystal water kept in each start before solvation (md/system/waters_<start>.pdb); in X and N '
                       'the ring itself lies at these positions.  place: within 2.5 A of the crystal positions of C1 to '
                       'C6 (the defaults of section 9D); wide: within 3.5 A of those of C1 to C6 and N1', 'values': ring},
           'midpoints_from_starts_A': {
               'd1': r6((float(q['d1_N1_F131CZ_A'][0]) + float(q['d1_N1_F131CZ_A'][2])) / 2),
               'd2': r6((float(q['d2_N1_N67OD1_A'][0]) + float(q['d2_N1_N67OD1_A'][2])) / 2)}}
    text = json.dumps(res, indent=1) + '\n'
    if args.out:
        Path(args.out).mkdir(parents=True, exist_ok=True)
        (Path(args.out) / 'check_starts.json').write_text(text)
    for s in 'XNR':
        i = 'XNR'.index(s)
        print(f'start {s}: ligand RMSD {q["lig_rmsd_A"][i]:.4f} A (labelling {extra[s]["lig_rmsd_labelling"]}; '
              f'starts.csv {ref[s]["rmsd_to_crystal_A"]}), d1 {q["d1_N1_F131CZ_A"][i]:.3f}, d2 {q["d2_N1_N67OD1_A"][i]:.3f}, '
              f'state {q["state"][i]}, CA RMSD after the fit {q["ca_rmsd_A"][i]:.1e} A')
    print(f'{res["n_agree"]} of {res["n_compared"]} values agree with md/system/starts.csv '
          f'(largest difference {worst["A"]:.4f} A, {worst["deg"]:.3f} deg)')
    for r in n2rows:
        print(f'start {r["start"]}: signed volume at N2 {r["n2_volume"]:+.4f} (review {r["review_value"]:+.2f}), sign '
              f'{r["n2_sign"]:+d}; Zn to Glu106 OE1 {extra[r["start"]]["d_Zn_E106OE1_A"]:.3f}, OE2 '
              f'{extra[r["start"]]["d_Zn_E106OE2_A"]:.3f} A; crystal water in the place of the ring (2.5 A, C1 to C6) '
              f'{ring[r["start"]]["n_water_ring_place"]} (nearest {ring[r["start"]]["d_water_ring_place_min_A"]:.3f} A), '
              f'wider count (3.5 A, C1 to C6 and N1) {ring[r["start"]]["n_water_ring_wide_3p5A"]} '
              f'(nearest {ring[r["start"]]["d_water_ring_wide_min_A"]:.3f} A)')
    print(f'signed volume at N2: the three starts agree with the review within 0.005 and with the signs of the '
          f'specification: {res["n2_configuration_of_the_starts"]["all_agree"]}')
    print('references of gate G-V against the crystal start: '
          + ', '.join(f'{r["reference_in_the_specification_deg"]:g} ({r["crystal_start_deg"]:.3f})' for r in gref)
          + f'; all agree within 0.05 deg: {all(r["agrees"] for r in gref)}; the crystal start meets the geometric '
          f'criteria: {res["crystal_start_against_the_geometric_criteria_of_gate_GV"]["all_met"]}')
    return 0 if (res['all_agree'] and all(r['agrees'] for r in gref)
                 and res['n2_configuration_of_the_starts']['all_agree']) else 1


# ----------------------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run_dir', nargs='?', help='run directory (solute_segNNNN.dcd, solute.pdb, ...)')
    ap.add_argument('--campaign', metavar='DIR', help='analyse every run directory below DIR and write the campaign table')
    ap.add_argument('--check-starts', action='store_true', help='apply the definitions to the three starts of md/system/')
    ap.add_argument('--system', metavar='DIR', help='system folder with selections.json and reference.pdb')
    ap.add_argument('--start', choices=('X', 'N', 'R'), help='start of the run (default: manifest, then folder name)')
    ap.add_argument('--model', help='zinc model (default: manifest, then the parent folder name)')
    ap.add_argument('--planned-ns', type=float, help='planned length of the run (default: manifest, then 100; '
                                                     'for the validation run X0: 5)')
    ap.add_argument('--solute-interval-ps', type=float, help='interval of the solute trajectory (default: manifest, then 10)')
    ap.add_argument('--full-interval-ps', type=float, help='interval of the full trajectory (default: manifest, then 100)')
    ap.add_argument('--out', metavar='DIR', help='output folder (default: RUN_DIR/analysis, or DIR/analysis_campaign)')
    ap.add_argument('--reuse', action='store_true', help='campaign: use an existing analysis/summary.json of a run')
    ap.add_argument('--no-figures', action='store_true')
    ap.add_argument('--no-hash-trajectories', action='store_true', help='do not compute the sha256 of the DCD files')
    ap.add_argument('--skip-name-check', action='store_true', help='go on when atom names at mapped positions differ')
    args = ap.parse_args(argv)
    warnings.filterwarnings('ignore', module='MDAnalysis')
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    try:
        if args.check_starts:
            return check_starts(args)
        if args.campaign:
            if args.start or args.system or args.run_dir:
                ap.error('--campaign takes neither --start nor --system nor a run directory')
            campaign(args.campaign, args)
            return 0
        if not args.run_dir:
            ap.error('give a run directory, --campaign DIR or --check-starts')
        analyse_run(args.run_dir, args)
        return 0
    except AnalysisError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
