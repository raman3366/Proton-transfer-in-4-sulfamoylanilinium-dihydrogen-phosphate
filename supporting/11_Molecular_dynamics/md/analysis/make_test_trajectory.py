#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""make_test_trajectory.py -- synthetic test of md/analysis/analyse_md.py with planted values.

  make_test_trajectory.py --out DIR        write the synthetic systems, runs and the planted values below DIR
  make_test_trajectory.py --compare DIR    compare what analyse_md.py wrote below DIR with the planted values;
                                           writes DIR/test_report.json and DIR/test_report.md; exit code 1 if
                                           any value is outside its tolerance

Between the two calls analyse_md.py is run on DIR/runs (see the end of md/analysis/SPEC.md).  This script
never imports analyse_md.py: the planted values come from (a) the construction itself, (b) md/system/starts.csv,
written by the system builder with its own code, and (c) a second set of geometry formulas in this file,
applied to the coordinates before they are rotated, wrapped and stored in single precision.

The synthetic system (410 atoms) follows the conventions of the campaign: the 258 CA atoms of 6RL9 and the named
site atoms as the "protein", the zinc, the 18-atom ligand, 40 water molecules and 2 sodium ions; five of the
water molecules lie between the protein and the zinc in the topology, so that an index of the full system is
not its position in the solute trajectory.  Files are written with OpenMM (DCDFile, PDBFile), as the production
driver will write them.  Coordinates in Angstrom, angles in degrees, times in ns.

What is planted (run by run in DIR/planted/):
  - the ligand in the crystal pose, then moved rigidly to the rotated pose, with a short excursion, a stretch of
    intermediate frames and a return (first passage, transitions, 90 % rule at its boundary);
  - atom labels swapped between equivalent atoms (ring, sulfonyl oxygens): symmetry-corrected RMSD unchanged;
  - the sulfonamide hydrogen turned away in every fifth frame: hydrogen-bond fraction 0.8;
  - the zinc moved 0.300 A away from N2 in every tenth frame: Zn-N2 known;
  - a water at 2.400 A from the zinc in every tenth frame of the full trajectory, 3.690 A otherwise; 0 to 3 water
    molecules at 3.000 A from the aniline N;
  - CA atoms scaled about their centroid by 1 +/- 0.01: CA RMSD and RMSF known, the fit unchanged;
  - every frame rotated at random and shifted; water wrapped into a truncated-octahedral cell whose size changes
    with time; the ligand and the zinc moved by lattice vectors in some frames;
  - runs in one or several segments, one run unfinished with a broken last frame, one run with repeated frames;
    manifests and DCD headers as md/run_md.py writes them (steps counted through the run) in X1, R1, N1 and X9,
    a manifest with plain keys in X2 and X3, DCD headers as OpenMM's own reporter writes them in X3.
Added with version 1.1 of the specification (gate G-V, the validation run X0, the optional key asn67_ND2):
  - runs/test/X0: a validation run (5 ns) in the crystal geometry beside the production runs of the campaign;
    gate G-V passes, and the campaign must leave the run out;
  - runs_gv/<case>/X0: validation runs of a second synthetic system (411 atoms, with ND2 of Asn67 and the key
    asn67_ND2): `edge` with every new criterion just inside its limit (gate passes), `angle`, `sulfonyl`, `water`,
    `csn`, `tau` with exactly one new criterion just outside (gate fails), `long` (6 ns, manifest target 6 ns:
    two warnings), `short` (3 ns, no manifest: provisional);
  - runs_gv/gate/R1: a production run of that system with the Asn67 amide turned in every fourth frame (OD1 and
    ND2 exchanged), in the rotated and in the crystal pose.
  The poses are bent for these runs by rigid rotations: ND1 of His119 about the zinc (one N-Zn-N angle), the ring
  about an axis through S (C4-S-N2) and about the C4-S bond (tau_fold); the zinc is put 2.500 A from the nearer
  sulfonyl oxygen in chosen frames.  The gate is evaluated a second time here, with the generator's own code
  (g_gate), from the planted values of every frame.
Added with version 1.2 of the specification (reported quantities of sections 9A to 9E):
  - every run: the signed volume at N2 and its sign, the zinc to OE1 and OE2 of Glu106, the water oxygens within
    3.5 A of the crystal positions of C1 to C6 and N1 (all by second formulas in this file); the stretches of one
    sign, the sign changes that last, the hydrogen bond before and after the first passage (own code here);
  - runs_v12/conv/R1 (22 ns, with the equilibration files of the full protocol and the manifest of the driver):
    rotated pose, then crystal pose; the hydrogen on N2 mirrored through the plane Zn-N2-S for 99 frames (too short
    to count) and for exactly 100 frames; OG1 of Thr199 moved away in chosen frames (hydrogen bond 301 of 600 before
    the first passage, 800 of 1600 after it); OE1 of Glu106 2.500 A and 2.700 A, OE2 2.550 A from the zinc in chosen
    frames; water at the centre of the crystal ring, 3.400 A from the crystal N1 and 3.600 A from the crystal C3;
  - runs_v12/flip/X1 (3 ns, no manifest, equilibration files with headers as OpenMM's own reporter writes them): the
    pose turns during the equilibration, the first production frame is intermediate with the other sign at N2;
  - runs_v12/scaled/N1 (2 ns, the manifest of a shortened equilibration: stages of 5, 15, 15, 20, 20 and 25 ps, so
    that files hold no frame, one frame or two, and the last frame lies 5 ps before the production);
  - runs_v12/glu/X0: a validation run with OE1 of Glu106 within 2.6 A in 6 % of the solute frames and in no full
    frame; gate G-V must pass, because Glu106 is reported and is no criterion of the solute fraction;
  - runs_v12/ext/R1 (2 ns): an equilibration whose held part was extended twice, as md/run_md.py of 6 Oct 2026
    does it (steps eq03x1 and eq03x2 between eq03 and eq04, 3000 ps in all); the manifest carries the driver's
    own count of water at the crystal ring positions at the end of eq03, eq03x1 and eq03x2 (0, 0, 1).
Added with version 1.3 of the specification (section 9D: the water count that the driver acts on):
  - every run: the water oxygens within 2.5 A of the crystal positions of C1 to C6 (n_water_ring_place) and the
    distance of the nearest one, beside the wider count of version 1.2 (3.5 A, C1 to C6 and N1), which is now
    named n_water_ring_wide_3p5A; both by second formulas in this file;
  - runs_v12/ext/R1: the same frames as in version 1.2; its manifest now carries the settings of the hydration
    check as the driver recorded them before the plan's lines of 6 Oct 2026, 8:01 and 8:56 AM (seven site atoms,
    3.5 A), so the analysis must follow the manifest, agree with the driver's counts (0, 0, 1) and warn once;
  - runs_v13/place/R1 (2 ns, new): an equilibration whose held part was extended three times, with the records and
    the settings of the hydration check as md/run_md.py writes them now (six ring carbons, 2.5 A); water planted
    2.400 A from the crystal C5 (inside the place of the ring) and 2.600 A from the crystal C2 (outside it, inside
    the wider count), and 2.300 A from the crystal N1 (in the wider count only); the driver's counts at the end of
    eq03, eq03x1, eq03x2, eq03x3 are 0, 0, 0, 1;
  - runs_v13/default/R1 (2 ns, new): the same production frames with a manifest of plain keys, which says nothing
    of the hydration check: the analysis must use the defaults of the specification.
Seed of the random rotations: 20261005 (equilibration frames: the same seed with a third number, 12).
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

SEED = 20261005
HERE = Path(__file__).resolve().parent
SYSD = HERE.parent / 'system'
LIG_HEAVY = ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S')
LIG_ALL = LIG_HEAVY + ('H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21')
KEEP_ND2 = {67: ('CA', 'OD1', 'ND2')}       # the second synthetic system also holds ND2 of Asn67
RING_ATOMS = ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'H2', 'H3', 'H5', 'H6', 'HN11', 'HN12')
KEEP = {67: ('CA', 'OD1'), 94: ('CA', 'NE2'), 96: ('CA', 'NE2'), 106: ('CA', 'OE1', 'OE2'), 119: ('CA', 'ND1'),
        131: ('CA', 'CZ'), 199: ('N', 'H', 'CA', 'OG1', 'HG1')}
RING_SWAP = {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3', 'H2': 'H6', 'H6': 'H2', 'H3': 'H5', 'H5': 'H3'}
OXY_SWAP = {'O1': 'O2', 'O2': 'O1'}
EPS = 0.01                      # scaling of the CA atoms about their centroid
DZN = 0.300                     # A, zinc moved away from N2 in every tenth frame
DW_NEAR, DW_FAR = 2.400, 3.690  # A, the water next to the zinc
DS_NEAR = 2.500                 # A, zinc to the nearer sulfonyl oxygen in the frames chosen for it
RHO_IN, RHO_OUT = 3.000, 6.000  # A, the water molecules at the aniline N
BOX_L, BOX_AMP = 90.0, 0.2      # A, edge of the truncated-octahedral cell and its variation
STEP_PS, SOL_STEPS, FULL_STEPS = 0.004, 2500, 25000
TOL = {'A': 0.001, 'deg': 0.01, 'vol': 0.001}
TOL_STARTS = {'A': 0.0015, 'deg': 0.06}     # starts.csv is rounded to three decimals
# gate G-V as SPEC.md, section 11, prints it (typed from the specification, not taken from analyse_md.py)
GV_ANGLE_REF = {'ang_N2_Zn_H94NE2_deg': 108.9, 'ang_N2_Zn_H96NE2_deg': 114.3, 'ang_N2_Zn_H119ND1_deg': 116.8,
                'ang_H94NE2_Zn_H96NE2_deg': 103.8, 'ang_H94NE2_Zn_H119ND1_deg': 113.6, 'ang_H96NE2_Zn_H119ND1_deg': 98.4}
GV_ANGLE_TOL, GV_CSN, GV_TAU = 10.0, (108.8, 3.0), (132.1, 20.0)
GV_ZN_N2, GV_ZN_HIS = (1.83, 2.03), (1.91, 2.12)
# version 1.2 (typed from SPEC.md, sections 9A to 9E)
# version 1.3 (typed from SPEC.md, section 9D): the count the driver acts on (defaults) and the wider count of 1.2
RING_PLACE, CUT_PLACE = ('C1', 'C2', 'C3', 'C4', 'C5', 'C6'), 2.5
RING_WIDE, CUT_WIDE = ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'), 3.5
D_RING_D, D_RING_E = 2.400, 2.600                           # A, planted water above the crystal C5 and below the crystal C2
D_RING_F = 2.300                                            # A, planted water from the crystal N1, away from the ring
SIGN_OF_START = {'X': 1, 'N': 1, 'R': -1}                   # sign of the signed volume at N2 in the built starts
STATE_OF_START = {'X': 'C', 'N': 'C', 'R': 'R'}
PROTOCOL_STAGES = ((0.0, 100.0), (100.0, 400.0), (400.0, 700.0), (700.0, 1100.0), (1100.0, 1500.0), (1500.0, 2000.0))
D_OE1_NEAR, D_OE1_FAR, D_OE2_NEAR = 2.500, 2.700, 2.550     # A, zinc to the oxygens of Glu106 in the frames chosen
D_OG1_AWAY = 2.000                                          # A, OG1 of Thr199 moved away from N2 (no hydrogen bond)
D_RING_B, D_RING_C = 3.400, 3.600                           # A, planted water from the crystal N1 and the crystal C3


# ----------------------------------------------------------------------------------------------------------
# second set of geometry formulas (not those of analyse_md.py)
# ----------------------------------------------------------------------------------------------------------
def g_dist(a, b):
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))


def g_angle(a, b, c):
    u, v = np.asarray(a, float) - b, np.asarray(c, float) - b
    return math.degrees(math.acos(max(-1.0, min(1.0, float(u @ v) / (math.sqrt(u @ u) * math.sqrt(v @ v))))))


def g_torsion(a, b, c, d):
    """Torsion from the normals of the two planes (IUPAC sign)."""
    b1, b2, b3 = np.asarray(b, float) - a, np.asarray(c, float) - b, np.asarray(d, float) - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return math.degrees(math.atan2(math.sqrt(b2 @ b2) * float(b1 @ n2), float(n1 @ n2)))


def g_volume(zn, n2, h, s):
    """Signed volume at N2: determinant of the unit vectors from N2 to Zn, to HN21 and to S (rows, in this order),
    written out by the rule of Sarrus."""
    rows = []
    for a in (zn, h, s):
        v = [float(a[i]) - float(n2[i]) for i in range(3)]
        m = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        rows.append([c / m for c in v])
    (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = rows
    return a1 * b2 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1 - a2 * b1 * c3 - a1 * b3 * c2


def g_sign(v):
    return 1 if v > 0 else -1 if v < 0 else 0


def g_stretches(signs, start_sign, dwell=100):
    """Section 9A by the generator's own code: stretches of one assigned sign as [sign, first frame, frames];
    the long ones; the changes of sign in (sign of the start, long stretches...)."""
    import itertools
    assigned, cur = [], start_sign
    for v in signs:
        cur = v if v != 0 else cur
        assigned.append(cur)
    out, k = [], 0
    for sg, grp in itertools.groupby(assigned):
        m = len(list(grp))
        out.append([sg, k, m])
        k += m
    longs = [x for x in out if x[2] >= dwell]
    changes, prev = [], start_sign
    for sg, first, m in longs:
        if sg != prev:
            changes.append({'frame': first, 'time_ns': round((first + 1) * STEP_PS * SOL_STEPS / 1000.0, 6), 'from_sign': prev, 'to_sign': sg})
        prev = sg
    return assigned, out, longs, changes


def g_first_passage(states, start_state, dwell=100):
    """First frame of the first visit of at least `dwell` frames in the opposite state (section 4), own code:
    intermediate frames take the state before them, before the first frame the state of the start."""
    opp = 'R' if start_state == 'C' else 'C'
    cur, run, first = start_state, 0, None
    for n, st in enumerate(states):
        a = st if st in 'CR' else cur
        if a != cur or n == 0:
            run, first = 0, n
        cur = a
        run += 1
        if cur == opp and run >= dwell:
            return first
    return None


def g_rmsd(lig, ref, perm):
    return math.sqrt(sum(g_dist(lig[n], ref[perm.get(n, n)]) ** 2 for n in LIG_HEAVY) / len(LIG_HEAVY))


def g_rotate(x, origin, axis, deg):
    """Point x turned about the line through `origin` along `axis` by `deg` (right-handed; Rodrigues)."""
    k = np.asarray(axis, float) / math.sqrt(float(np.dot(axis, axis)))
    v = np.asarray(x, float) - origin
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return origin + v * c + np.cross(k, v) * s + k * float(k @ v) * (1.0 - c)


def g_gate(srows, frows):
    """Gate G-V from the planted values of every frame, after SPEC.md section 11 (second implementation).
    Returns the means used, the deviations from the references, and the criteria that are met."""
    n, nf = len(srows), len(frows)

    def mean(c):
        return sum(r[c] for r in srows) / n
    m = {c: mean(c) for c in ('lig_rmsd_A', 'd_Zn_N2_A', 'd_Zn_H94NE2_A', 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A',
                              'ang_C4_S_N2_deg', 'tau_fold_deg') + tuple(GV_ANGLE_REF)}
    hb = sum(r['hb_N2H_T199OG1'] for r in srows)
    cs = sum(r['fifth_sulfonylO'] for r in srows)
    cf = sum(r['fifth_any'] for r in frows)
    met = {'GV1 ligand RMSD': m['lig_rmsd_A'] <= 1.0,
           'GV2 hydrogen bond': 5 * hb >= 4 * n,
           'GV3 Zn-N2': GV_ZN_N2[0] <= m['d_Zn_N2_A'] <= GV_ZN_N2[1]}
    for c in ('d_Zn_H94NE2_A', 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A'):
        met[f'GV3 {c}'] = GV_ZN_HIS[0] <= m[c] <= GV_ZN_HIS[1]
    dev = {c: m[c] - ref for c, ref in GV_ANGLE_REF.items()}
    for c in GV_ANGLE_REF:
        met[f'GV3 {c}'] = abs(dev[c]) <= GV_ANGLE_TOL
    met['GV3 fifth ligand, solute trajectory'] = 20 * cs <= n
    met['GV3 fifth ligand, full trajectory'] = nf > 0 and 20 * cf <= nf
    dev['ang_C4_S_N2_deg'], dev['tau_fold_deg'] = m['ang_C4_S_N2_deg'] - GV_CSN[0], m['tau_fold_deg'] - GV_TAU[0]
    met['GV4 ang_C4_S_N2_deg'] = abs(dev['ang_C4_S_N2_deg']) <= GV_CSN[1]
    met['GV4 tau_fold_deg'] = abs(dev['tau_fold_deg']) <= GV_TAU[1]
    return {'means': m, 'deviations': dev, 'met': met, 'not_met': sorted(k for k, v in met.items() if not v),
            'passed': all(met.values()), 'hb_count': hb, 'fifth_solute_count': cs, 'fifth_full_count': cf,
            'n_solute': n, 'n_full': nf}


LABELLINGS = ({}, {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3'}, {'O1': 'O2', 'O2': 'O1'},
              {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3', 'O1': 'O2', 'O2': 'O1'})


def read_pdb(path):
    out = []
    for line in open(path):
        if line.startswith(('ATOM  ', 'HETATM')):
            out.append((line[12:16].strip(), line[17:21].strip(), int(line[22:26]),
                        np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rigid_onto(P, Q):
    """Rotation and translation that put the points P onto Q by least squares (used to build the moved pose)."""
    cp, cq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    return lambda x: (np.asarray(x, float) - cp) @ R + cq


def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def cell(L):
    """Truncated octahedron as OpenMM and Amber store it: a = b = c = L, angles 109.4712 degrees (rows, A)."""
    return np.array([[L, 0.0, 0.0], [-L / 3, 2 * math.sqrt(2) * L / 3, 0.0],
                     [-L / 3, -math.sqrt(2) * L / 3, math.sqrt(6) * L / 3]])


# ----------------------------------------------------------------------------------------------------------
# the synthetic system
# ----------------------------------------------------------------------------------------------------------
class System:
    def __init__(self, with_nd2=False):
        prot = read_pdb(SYSD / 'tleap' / 'protein_H.pdb')
        self.with_nd2 = with_nd2
        keep = {**KEEP, **KEEP_ND2} if with_nd2 else KEEP
        self.atoms = []                                       # (name, resname, resnum, kind, xyz in the 6RL9 frame)
        for name, resn, resi, x in prot:
            if name in keep.get(resi, ('CA',)) and resn != 'ACE':
                self.atoms.append([name, resn, resi, 'protein', x])
        n_prot = len(self.atoms)
        for w in range(5):
            self.atoms += [[a, 'WAT', 401 + w, 'water', np.zeros(3)] for a in ('O', 'H1', 'H2')]
        zn = read_pdb(SYSD / 'zinc.pdb')[0][3]
        self.atoms.append(['ZN', 'ZN', 301, 'zinc', zn])
        self.lig = {s: {n: x for n, _, _, x in read_pdb(SYSD / f'ligand_{s}.pdb')} for s in 'XNR'}
        self.atoms += [[n, 'SAN', 304, 'ligand', self.lig['X'][n]] for n in LIG_ALL]
        self.atoms += [['Na+', 'Na+', 901 + i, 'ion', np.zeros(3)] for i in range(2)]
        for w in range(5, 40):
            self.atoms += [[a, 'WAT', 401 + w, 'water', np.zeros(3)] for a in ('O', 'H1', 'H2')]
        self.n = len(self.atoms)
        self.kind = [a[3] for a in self.atoms]
        self.name = [a[0] for a in self.atoms]
        self.x0 = np.array([a[4] for a in self.atoms])
        self.solute = [i for i in range(self.n) if self.kind[i] in ('protein', 'zinc', 'ligand')]
        self.izn = self.kind.index('zinc')
        self.ilig = {n: next(i for i in range(self.n) if self.kind[i] == 'ligand' and self.name[i] == n) for n in LIG_ALL}
        self.water_o = [i for i in range(self.n) if self.kind[i] == 'water' and self.name[i] == 'O']
        self.ions = [i for i in range(self.n) if self.kind[i] == 'ion']
        self.ca = [i for i in range(n_prot) if self.name[i] == 'CA']
        self.ca_resnum = [self.atoms[i][2] for i in self.ca]

        def site(resi, name):
            return next(i for i in range(n_prot) if self.atoms[i][2] == resi and self.name[i] == name)
        self.site = {'his94_NE2': site(94, 'NE2'), 'his96_NE2': site(96, 'NE2'), 'his119_ND1': site(119, 'ND1'),
                     'thr199_OG1': site(199, 'OG1'), 'thr199_HG1': site(199, 'HG1'), 'thr199_N': site(199, 'N'),
                     'thr199_H': site(199, 'H'), 'glu106_OE1': site(106, 'OE1'), 'glu106_OE2': site(106, 'OE2'),
                     'phe131_CZ': site(131, 'CZ'), 'asn67_OD1': site(67, 'OD1')}
        if with_nd2:
            self.site['asn67_ND2'] = site(67, 'ND2')
        self.protein_all = list(range(n_prot))
        self.protein_heavy = [i for i in range(n_prot) if not self.name[i].startswith('H')]
        self.protein_o = [i for i in self.protein_heavy if self.name[i].startswith('O')]
        self.ca_ref = self.x0[self.ca].copy()
        self.ca_centroid = self.ca_ref.mean(0)
        self.zn0 = zn.copy()
        self.S = {k: self.x0[i].copy() for k, i in self.site.items()}

        # poses of the ligand ({name: xyz}, 6RL9 frame)
        X = self.lig['X']
        fit = rigid_onto(np.array([X[n] for n in LIG_HEAVY]), np.array([self.lig['R'][n] for n in LIG_HEAVY]))
        u = (X['N1'] - self.S['phe131_CZ'])
        u /= np.linalg.norm(u)
        self.pose = {'X': X, 'T': {n: fit(X[n]) for n in LIG_ALL}, 'I': {n: X[n] + 3.0 * u for n in LIG_ALL},
                     'R': self.lig['R'], 'N': self.lig['N']}
        # directions for the planted water (fixed in the 6RL9 frame)
        self.u_zn = (zn - X['N1']) / np.linalg.norm(zn - X['N1'])
        self._bent = {}

    def bent_pose(self, pose, tau_deg, csn_deg):
        """The pose with the ring (and the aniline group) turned rigidly: by tau_deg about the C4-S bond, so that
        the torsion C3-C4-S-N2 grows by tau_deg, then by csn_deg about the axis through S normal to the plane
        C4-S-N2, so that the angle C4-S-N2 grows by csn_deg.  S, N2, O1, O2 and HN21 stay."""
        key = (pose, tau_deg, csn_deg)
        if key not in self._bent:
            P = {k: np.array(v, float) for k, v in self.pose[pose].items()}
            if tau_deg:
                t0 = g_torsion(P['C3'], P['C4'], P['S'], P['N2'])
                for sign in (1.0, -1.0):
                    Q = dict(P)
                    for k in RING_ATOMS:
                        Q[k] = g_rotate(P[k], P['C4'], P['S'] - P['C4'], sign * tau_deg)
                    if abs((g_torsion(Q['C3'], Q['C4'], Q['S'], Q['N2']) - t0 + 180.0) % 360.0 - 180.0 - tau_deg) < 1e-6:
                        break
                else:
                    raise AssertionError('torsion not turned as designed')
                P = Q
            if csn_deg:
                a0 = g_angle(P['C4'], P['S'], P['N2'])
                axis = np.cross(P['N2'] - P['S'], P['C4'] - P['S'])
                Q = dict(P)
                for k in RING_ATOMS:
                    Q[k] = g_rotate(P[k], P['S'], axis, csn_deg)
                assert abs(g_angle(Q['C4'], Q['S'], Q['N2']) - a0 - csn_deg) < 1e-6
                P = Q
            self._bent[key] = P
        return self._bent[key]

    def nd1_turned(self, deg):
        """ND1 of His119 turned about the zinc (crystal position), in the plane NE2(His96)-Zn-ND1 and away from
        NE2 of His96: the angle NE2(His96)-Zn-ND1(His119) grows by deg, Zn-ND1 keeps its length."""
        u, v = self.S['his96_NE2'] - self.zn0, self.S['his119_ND1'] - self.zn0
        new = g_rotate(self.S['his119_ND1'], self.zn0, np.cross(u, v), deg)
        assert abs(g_angle(self.S['his96_NE2'], self.zn0, new) - g_angle(self.S['his96_NE2'], self.zn0, self.S['his119_ND1']) - deg) < 1e-6
        return new

    def selections(self):
        return {'zn': self.izn, 'lig_heavy': {n: self.ilig[n] for n in LIG_HEAVY}, 'lig_HN2': self.ilig['HN21'],
                'lig_all': [self.ilig[n] for n in LIG_ALL], **self.site,
                'protein_CA': self.ca, 'protein_heavy': self.protein_heavy, 'protein_all': self.protein_all,
                'water_O': self.water_o, 'ions': self.ions, 'solute': self.solute,
                'residue_numbers_deposited_for_CA': self.ca_resnum}

    def topology(self, idx, L):
        import openmm.app as app
        import openmm.unit as u
        from openmm import Vec3
        top = app.Topology()
        group = {'protein': 'A', 'water': 'W', 'zinc': 'Z', 'ligand': 'Z', 'ion': 'I'}
        el = {'C': app.element.carbon, 'N': app.element.nitrogen, 'O': app.element.oxygen, 'H': app.element.hydrogen,
              'S': app.element.sulfur}
        ch, last_group, last, res = None, None, None, None
        for i in idx:                                         # a new chain whenever the kind of molecule changes, so
            name, resn, resi, kind, _ = self.atoms[i]         # that the PDB file keeps the atom order of the trajectory
            if group[kind] != last_group:
                ch = top.addChain(group[kind])
                last_group = group[kind]
            if (kind, resi) != last:
                res = top.addResidue(resn, ch, id=str(resi))
                last = (kind, resi)
            e = app.element.zinc if kind == 'zinc' else app.element.sodium if kind == 'ion' else el[name[0]]
            top.addAtom(name, e, res)
        B = cell(L) / 10.0
        top.setPeriodicBoxVectors(u.Quantity([Vec3(*B[k]) for k in range(3)], u.nanometer))
        return top


# ----------------------------------------------------------------------------------------------------------
# schedules of the runs
# ----------------------------------------------------------------------------------------------------------
def variant(n):
    return ('', 'ring', 'oxy', 'ring+oxy')[n % 4]


def schedule(run):
    """List of (pose, labels swapped) for every solute frame, and the description of the run."""
    def xv(n):
        return ('X', variant(n))

    def tv(n):
        return ('T', 'ring' if n % 2 else '')
    if run in ('X1', 'X3'):
        fr = []
        for n in range(3000):
            if 300 <= n < 350 or 900 <= n < 950 or 970 <= n < 1400 or 1600 <= n:
                fr.append(tv(n))
            elif 950 <= n < 970:
                fr.append(('I', ''))
            else:
                fr.append(xv(n))
        if run == 'X3':
            fr[2999] = xv(2999)
        return fr, {'start': 'X', 'planned_ns': 30.0, 'segments': [1000, 1500, 500] if run == 'X1' else [2000, 1000]}
    if run == 'X2':
        return [xv(n) for n in range(3000)], {'start': 'X', 'planned_ns': 30.0, 'segments': [3000]}
    if run == 'R1':
        return ([('R', '')] * 800 + [xv(n) for n in range(800, 3000)],
                {'start': 'R', 'planned_ns': 30.0, 'segments': [1200, 1800]})
    if run == 'N1':
        return [('N', '')] * 1234, {'start': 'N', 'planned_ns': 100.0, 'segments': [1000, 234], 'broken_last_frame': True}
    if run == 'X9':                                           # frames 150..199 written twice (a restart that overlaps)
        return [xv(n) for n in range(200)], {'start': 'X', 'planned_ns': 2.5, 'segments': [200, 50], 'repeat_from': 150}
    if run in GV_CASES:                                       # version 1.1: validation runs and the run with ND2
        c = GV_CASES[run]
        if run == 'gate/R1':
            fr = [('R', '')] * 250 + [xv(n) for n in range(250, 400)]
        else:
            fr = [xv(n) for n in range(c['frames'])]
        return fr, {'start': run.split('/')[1][0], 'planned_ns': c['planned_ns'], 'segments': c['segments'],
                    'manifest': c['manifest'], 'target_ns': c.get('target_ns'), 'opt': c.get('opt', {})}
    if run in V12_CASES:                                      # version 1.2
        c = V12_CASES[run]
        return ([c['pose'](n) for n in range(c['frames'])],
                {'start': run.split('/')[1][0], 'planned_ns': c['planned_ns'], 'segments': c['segments'],
                 'manifest': c['manifest'], 'target_ns': c.get('target_ns'), 'opt': c.get('opt', {})})
    raise KeyError(run)


def _xv(n):
    return ('X', variant(n))


def _s20(k):
    return [3 + 20 * i for i in range(k)]


# The runs added with version 1.1.  Key: <model>/<run>.  planned_ns is what the analysis must find (5 ns for X0:
# from the specification, never from the manifest).  opt: what is bent or planted in the frames, see
# frame_coordinates() and place_water().  Every X0 here has 10 ps and 100 ps between frames, like every run.
GV_CASES = {
    # crystal geometry, beside the production runs of the campaign: gate passes
    'test/X0': {'folder': 'runs', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'driver',
                'target_ns': 5.0, 'opt': {'water_near': [7, 31]}},
    # every new criterion just inside its limit: gate passes.  25 of 500 solute frames (5.0 %) with a sulfonyl
    # oxygen at 2.500 A, the last of them in a full frame; with one water frame, 2 of 50 full frames (4 %)
    'edge/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [300, 200], 'planned_ns': 5.0, 'manifest': 'plain-bare',
                'opt': {'nd1': 9.0, 'csn': 2.5, 'tau': 19.0, 'sulfonyl': _s20(24) + [499], 'water_near': [7],
                        'amide_turned': list(range(7, 500, 25))}},
    # one new criterion just outside its limit, everything else as in the crystal
    'angle/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'driver',
                 'target_ns': 5.0, 'opt': {'nd1': 11.0, 'water_near': [7, 31]}},
    'sulfonyl/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'plain',
                    'opt': {'sulfonyl': _s20(25) + [13], 'water_near': [7, 31]}},          # 26 of 500 = 5.2 %
    'water/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [250, 250], 'planned_ns': 5.0, 'manifest': 'driver',
                 'target_ns': 5.0, 'opt': {'sulfonyl': [3, 103, 203, 303, 499], 'water_near': [7, 31]}},  # 3 of 50 full
    'csn/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'plain-bare',
               'opt': {'csn': 3.5, 'water_near': [7, 31]}},
    'tau/X0': {'folder': 'runs_gv', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'driver',
               'target_ns': 5.0, 'opt': {'tau': 21.0, 'water_near': [7, 31]}},
    # the driver was told 6 ns: the analysis keeps 5 ns, warns twice (manifest, more frames than planned).
    # 3 of 60 full frames with a water at 2.400 A: exactly 5 %, the limit of the full-trajectory fraction
    'long/X0': {'folder': 'runs_gv', 'frames': 600, 'segments': [500, 100], 'planned_ns': 5.0, 'manifest': 'driver',
                'target_ns': 6.0, 'opt': {'water_near': [7, 31, 55]}},
    # no manifest at all, 3 ns of 5: role and length from the folder name and the specification; provisional
    'short/X0': {'folder': 'runs_gv', 'frames': 300, 'segments': [300], 'planned_ns': 5.0, 'manifest': 'none',
                 'opt': {'water_near': [7]}},
    # a production run of the system with ND2: rotated pose, then crystal pose; the amide turned in every 4th frame
    'gate/R1': {'folder': 'runs_gv', 'frames': 400, 'segments': [400], 'planned_ns': 4.0, 'manifest': 'driver',
                'target_ns': 4.0, 'opt': {'water_near': [7], 'amide_turned': list(range(1, 400, 4))}},
}


# The runs added with version 1.2.  opt['v12']: the hydrogen on N2 stays in place unless the frame is listed in
# `mirror` (hydrogen mirrored through the plane Zn-N2-S: the signed volume changes its sign and nothing else);
# `og1_away`: OG1 of Thr199 moved 2.000 A away from N2 (no hydrogen bond, signed volume unchanged); `oe1_near`,
# `oe1_far`, `oe2_near`: OE1 or OE2 of Glu106 put 2.500, 2.700 or 2.550 A from the zinc; `ring`: full frame ->
# planted water (A: centre of the crystal ring, B: 3.400 A from the crystal N1, C: 3.600 A from the crystal C3;
# version 1.3: D: 2.400 A from the crystal C5, E: 2.600 A from the crystal C2, both along the normal of the ring;
# F: 2.300 A from the crystal N1, away from the ring: near the aniline N only, so in the wider count alone).
# hydration (version 1.3): the settings of the driver's hydration check that the manifest carries.
# equil: the equilibration files; stages (start, end) in ps, the pose of every equilibration frame, how the DCD
# headers count, whether the manifest lists the stages.
_CONV_RING = {**{j: 'C' for j in range(20)}, **{j: 'BC' for j in range(20, 40)}, **{j: 'ABC' for j in range(40, 60)},
              **{j: 'C' for j in range(60, 220)}}
# version 1.3, full frames of the two new runs: E alone (outside the place of the ring), D (inside), D with E and F,
# A and D
_PLACE_RING = {**{j: 'CE' for j in range(5)}, **{j: 'CD' for j in range(5, 10)}, **{j: 'CDEF' for j in range(10, 15)},
               **{j: 'ACD' for j in range(15, 20)}}
V12_CASES = {
    'conv/R1': {'folder': 'runs_v12', 'frames': 2200, 'segments': [1200, 1000], 'planned_ns': 22.0, 'manifest': 'driver',
                'target_ns': 22.0, 'as_built_frame': ['R', 2],
                'pose': lambda n: ('R', '') if n < 600 else _xv(n),
                'opt': {'v12': True, 'water_near': [7],
                        'mirror': list(range(200, 299)) + list(range(400, 500)),
                        'og1_away': [n for n in range(600) if n % 4 == 1] + [n for n in range(600, 2200) if n % 2 == 0],
                        'oe1_near': [n for n in range(1000) if n % 50 == 9] + [n for n in range(1000, 2200) if n % 50 == 13],
                        'oe1_far': [n for n in range(2200) if n % 100 == 33],
                        'oe2_near': [n for n in range(2200) if n % 100 == 57], 'ring': _CONV_RING},
                'equil': {'stages': PROTOCOL_STAGES, 'in_manifest': True, 'header': 'driver', 'pose': lambda e: ('R', ''),
                          'opt': {'v12': True, 'water_near': [],
                                  'ring': {**{j: 'C' for j in range(4)}, **{j: 'BC' for j in range(4, 7)},
                                           **{j: 'ABC' for j in range(7, 20)}}}}},
    'flip/X1': {'folder': 'runs_v12', 'frames': 300, 'segments': [300], 'planned_ns': 100.0, 'manifest': 'none',
                'as_built_frame': ['X', 152],
                'pose': lambda n: ('I', '') if n < 5 else (('T', 'ring' if n % 2 else '') if n < 150 else _xv(n)),
                'opt': {'v12': True, 'water_near': [7], 'mirror': list(range(5)),
                        'og1_away': [n for n in range(150, 300) if n % 3 == 0], 'ring': {j: 'C' for j in range(30)}},
                'equil': {'stages': PROTOCOL_STAGES, 'in_manifest': False, 'header': 'reporter',
                          'pose': lambda e: _xv(e) if e < 120 else ('T', ''),
                          'opt': {'v12': True, 'water_near': [], 'ring': {j: 'C' for j in range(20)}}}},
    'scaled/N1': {'folder': 'runs_v12', 'frames': 200, 'segments': [200], 'planned_ns': 2.0, 'manifest': 'driver',
                  'target_ns': 2.0, 'as_built_frame': ['N', 1], 'pose': lambda n: ('N', ''),
                  'opt': {'v12': True, 'water_near': [7], 'ring': {j: 'C' for j in range(20)}},
                  'equil': {'stages': ((0.0, 5.0), (5.0, 20.0), (20.0, 35.0), (35.0, 55.0), (55.0, 75.0), (75.0, 100.0)),
                            'in_manifest': True, 'header': 'driver', 'pose': lambda e: ('N', ''),
                            'opt': {'v12': True, 'water_near': [], 'ring': {}}}},
    'ext/R1': {'folder': 'runs_v12', 'frames': 200, 'segments': [200], 'planned_ns': 2.0, 'manifest': 'driver',
               'target_ns': 2.0, 'as_built_frame': ['R', 2], 'pose': lambda n: ('R', ''),
               'opt': {'v12': True, 'water_near': [7], 'ring': {j: 'BC' for j in range(20)}},
               # version 1.3: the manifest names the settings of the driver before 8:01 AM of 6 Oct 2026 (seven atoms, 3.5 A)
               'hydration': {'site_atoms': list(RING_WIDE), 'water_cutoff_A': 3.5, 'form': 'settings only'},
               'equil': {'steps': (('eq01', 0.0, 100.0), ('eq02', 100.0, 400.0), ('eq03', 400.0, 700.0), ('eq03x1', 700.0, 1200.0),
                                   ('eq03x2', 1200.0, 1700.0), ('eq04', 1700.0, 2100.0), ('eq05', 2100.0, 2500.0), ('eq06', 2500.0, 3000.0)),
                         'in_manifest': True, 'header': 'driver', 'pose': lambda e: ('R', ''),
                         # what the driver records at the end of a step: water oxygens at the ring positions, nearest (A), decision
                         'driver_counts': {'eq03': [0, 3.6, 'extend'], 'eq03x1': [0, 3.6, 'extend'], 'eq03x2': [1, 3.4, 'go on']},
                         'opt': {'v12': True, 'water_near': [],
                                 'ring': {**{j: 'C' for j in range(15)}, **{j: 'BC' for j in range(15, 30)}}}}},
    # version 1.3: the count the driver acts on, with the records of the driver as it writes them now
    'place/R1': {'folder': 'runs_v13', 'frames': 200, 'segments': [200], 'planned_ns': 2.0, 'manifest': 'driver',
                 'target_ns': 2.0, 'as_built_frame': ['R', 2], 'pose': lambda n: ('R', ''),
                 'hydration': {'site_atoms': list(RING_PLACE), 'water_cutoff_A': 2.5, 'form': 'driver'},
                 'opt': {'v12': True, 'water_near': [7], 'ring': _PLACE_RING},
                 'equil': {'steps': (('eq01', 0.0, 100.0), ('eq02', 100.0, 400.0), ('eq03', 400.0, 700.0), ('eq03x1', 700.0, 1200.0),
                                     ('eq03x2', 1200.0, 1700.0), ('eq03x3', 1700.0, 2200.0), ('eq04', 2200.0, 2600.0),
                                     ('eq05', 2600.0, 3000.0), ('eq06', 3000.0, 3500.0)),
                           'in_manifest': True, 'header': 'driver', 'pose': lambda e: ('R', ''),
                           # the driver's count at the end of a step: water oxygens within 2.5 A of the six ring carbons, nearest (A)
                           'driver_counts': {'eq03': [0, 2.6, 'extend'], 'eq03x1': [0, 2.6, 'extend'], 'eq03x2': [0, 2.6, 'extend'],
                                             'eq03x3': [1, 2.4, 'go on']},
                           # full frames 0 to 20: only the water outside the place; frame 21, the last of eq03x3, and on: one
                           # inside it, and one nearer to the aniline N that the wider count sees and the count of the driver does not
                           'opt': {'v12': True, 'water_near': [],
                                   'ring': {**{j: 'CE' for j in range(21)}, **{j: 'CDF' for j in range(21, 35)}}}}},
    # version 1.3: the same production frames with a manifest that says nothing of the hydration check (plain keys):
    # the analysis must use the defaults of the specification (six ring carbons, 2.5 A)
    'default/R1': {'folder': 'runs_v13', 'frames': 200, 'segments': [200], 'planned_ns': 2.0, 'manifest': 'plain',
                   'as_built_frame': ['R', 2], 'pose': lambda n: ('R', ''),
                   'opt': {'v12': True, 'water_near': [7], 'ring': _PLACE_RING}},
    'glu/X0': {'folder': 'runs_v12', 'frames': 500, 'segments': [500], 'planned_ns': 5.0, 'manifest': 'plain',
               'as_built_frame': ['X', 8], 'pose': _xv,
               'opt': {'v12': True, 'water_near': [7, 31], 'og1_away': [n for n in range(500) if n % 5 == 4],
                       'oe1_near': [2 + 10 * k for k in range(30)], 'ring': {j: 'C' for j in range(50)}}},
}
# What the design of the new runs fixes by hand (typed, not computed; checked against the frames as built when the
# runs are written, and against the analysis in --compare).  Times in ns, frames counted from 0.
V12_DESIGN = {
    'conv/R1': {
        'n2': {'count_positive': 1799, 'count_negative': 401, 'count_zero': 0, 'window_frames': 2000,
               'count_positive_window': 1799, 'count_negative_window': 201,
               'by_state': {'C': [1600, 1600, 0], 'R': [600, 199, 401], 'I': [0, 0, 0]},      # frames, positive, negative
               'n_stretches': 6, 'n_long_stretches': 5, 'change_frames': [400, 500, 600], 'change_times_ns': [4.01, 5.01, 6.01],
               'change_to_sign': [1, -1, 1], 'sign_first_frame': -1},
        'hb': {'first_passage_frame': 600, 'n_before': 600, 'count_before': 301, 'n_after': 1600, 'count_after': 800},
        'glu': {'count_OE1': 44, 'count_OE2': 22, 'count_either': 66, 'min_OE1': 2.5, 'min_OE2': 2.55, 'full_frames_protein_O': 20},
        'ring': {'counts': [0] * 20 + [1] * 20 + [2] * 20 + [j % 4 for j in range(60, 220)], 'mean': 300 / 220,
                 'min_distance': 1.382},
        # version 1.3: only water A (the centre of the ring) is within 2.5 A of a ring carbon
        'place': {'counts': [0] * 40 + [1] * 20 + [0] * 160, 'mean': 20 / 220, 'min_distance': 1.382},
        'start': {'state_first_frame': 'R', 'state_differs': False, 'n2_sign_first_frame': -1, 'n2_sign_differs': False},
        'equil': {'available': True, 'source': 'manifest equilibration.stages', 'n_solute': 200, 'n_full': 20, 'total_ps': 2000.0,
                  'frames_per_stage': [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]],
                  'time_first_ns': -1.99, 'time_last_ns': 0.0, 'time_frame_10_ns': -1.89,
                  'ring_counts': [0] * 4 + [1] * 3 + [2] * 13, 'ring_last_frame_of_stage_3': 1,
                  'place_counts': [0] * 7 + [1] * 13, 'place_last_frame_of_stage_3': 0,
                  'first_not_start': None, 'first_opposite': None, 'first_sign': None, 'state_last': 'R', 'sign_last': -1},
        'warnings': 0},
    'flip/X1': {
        'n2': {'count_positive': 150, 'count_negative': 150, 'count_zero': 0, 'window_frames': 300,
               'count_positive_window': 150, 'count_negative_window': 150,
               'by_state': {'C': [150, 150, 0], 'R': [145, 0, 145], 'I': [5, 0, 5]},
               'n_stretches': 2, 'n_long_stretches': 2, 'change_frames': [0, 150], 'change_times_ns': [0.01, 1.51],
               'change_to_sign': [-1, 1], 'sign_first_frame': -1},
        'hb': {'first_passage_frame': 5, 'n_before': 5, 'count_before': 0, 'n_after': 295, 'count_after': 245},
        'glu': {'count_OE1': 0, 'count_OE2': 0, 'count_either': 0, 'min_OE1': None, 'min_OE2': None, 'full_frames_protein_O': 0},
        'ring': None,
        'start': {'state_first_frame': 'I', 'state_differs': True, 'n2_sign_first_frame': -1, 'n2_sign_differs': True},
        'equil': {'available': True, 'source': 'protocol (constants of the specification)', 'n_solute': 200, 'n_full': 20,
                  'total_ps': 2000.0, 'frames_per_stage': [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]],
                  'time_first_ns': -1.99, 'time_last_ns': 0.0, 'time_frame_10_ns': -1.89, 'ring_counts': None,
                  'ring_last_frame_of_stage_3': None,
                  'first_not_start': [120, 5, 0, -0.79], 'first_opposite': [120, 5, 0, -0.79], 'first_sign': [120, 5, 0, -0.79],
                  'state_last': 'R', 'sign_last': -1},
        'warnings': 0},
    'scaled/N1': {
        'n2': {'count_positive': 200, 'count_negative': 0, 'count_zero': 0, 'window_frames': 200,
               'count_positive_window': 200, 'count_negative_window': 0,
               'by_state': {'C': [200, 200, 0], 'R': [0, 0, 0], 'I': [0, 0, 0]},
               'n_stretches': 1, 'n_long_stretches': 1, 'change_frames': [], 'change_times_ns': [], 'change_to_sign': [],
               'sign_first_frame': 1},
        'hb': {'first_passage_frame': None, 'n_before': None, 'count_before': None, 'n_after': None, 'count_after': None},
        'glu': {'count_OE1': 0, 'count_OE2': 0, 'count_either': 0, 'min_OE1': None, 'min_OE2': None, 'full_frames_protein_O': 0},
        'ring': None,
        'start': {'state_first_frame': 'C', 'state_differs': False, 'n2_sign_first_frame': 1, 'n2_sign_differs': False},
        'equil': {'available': True, 'source': 'manifest equilibration.stages', 'n_solute': 8, 'n_full': 0, 'total_ps': 100.0,
                  'frames_per_stage': [[0, 0], [1, 0], [1, 0], [2, 0], [2, 0], [2, 0]],
                  'times_ns': [-0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005],
                  'time_first_ns': -0.085, 'time_last_ns': -0.005, 'ring_counts': [], 'ring_last_frame_of_stage_3': None,
                  'first_not_start': None, 'first_opposite': None, 'first_sign': None, 'state_last': 'C', 'sign_last': 1},
        'warnings': 0},
    'ext/R1': {
        'n2': {'count_positive': 0, 'count_negative': 200, 'count_zero': 0, 'window_frames': 200,
               'count_positive_window': 0, 'count_negative_window': 200,
               'by_state': {'C': [0, 0, 0], 'R': [200, 0, 200], 'I': [0, 0, 0]},
               'n_stretches': 1, 'n_long_stretches': 1, 'change_frames': [], 'change_times_ns': [], 'change_to_sign': [],
               'sign_first_frame': -1},
        'hb': {'first_passage_frame': None, 'n_before': None, 'count_before': None, 'n_after': None, 'count_after': None},
        'glu': {'count_OE1': 0, 'count_OE2': 0, 'count_either': 0, 'min_OE1': None, 'min_OE2': None, 'full_frames_protein_O': 0},
        'ring': {'counts': [1] * 20, 'mean': 1.0, 'min_distance': 3.4},
        # version 1.3: with the settings of its manifest (seven atoms, 3.5 A) the count of the driver is the wider count
        'place': {'counts': [1] * 20, 'mean': 1.0, 'min_distance': 3.4},
        'settings': {'site_atoms': list(RING_WIDE), 'water_cutoff_A': 3.5, 'source': 'manifest hydration_check.settings',
                     'settings_are_the_defaults': False},
        'start': {'state_first_frame': 'R', 'state_differs': False, 'n2_sign_first_frame': -1, 'n2_sign_differs': False},
        'equil': {'available': True, 'source': 'manifest equilibration.stages', 'n_solute': 300, 'n_full': 30, 'total_ps': 3000.0,
                  'steps': ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq05', 'eq06'],
                  'frames_per_stage': [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]],
                  'time_first_ns': -2.99, 'time_last_ns': 0.0, 'time_frame_10_ns': -2.89,
                  'frame_70': [3, 1, -2.29], 'frame_170': [4, 0, -1.29],            # stage, extension, time of solute frames 70 and 170
                  'ring_counts': [0] * 15 + [1] * 15, 'ring_last_frame_of_stage_3': 0,
                  'place_counts': [0] * 15 + [1] * 15, 'place_last_frame_of_stage_3': 0,
                  'ring_at_the_end_of': {'eq03': 0, 'eq03x1': 0, 'eq03x2': 1}, 'driver_counts_agree': [True, True, True],
                  'first_not_start': None, 'first_opposite': None, 'first_sign': None, 'state_last': 'R', 'sign_last': -1},
        'warnings': 1},            # version 1.3: the settings of the manifest are not those of the plan (one warning)
    'place/R1': {
        'n2': {'count_positive': 0, 'count_negative': 200, 'count_zero': 0, 'window_frames': 200,
               'count_positive_window': 0, 'count_negative_window': 200,
               'by_state': {'C': [0, 0, 0], 'R': [200, 0, 200], 'I': [0, 0, 0]},
               'n_stretches': 1, 'n_long_stretches': 1, 'change_frames': [], 'change_times_ns': [], 'change_to_sign': [],
               'sign_first_frame': -1},
        'hb': {'first_passage_frame': None, 'n_before': None, 'count_before': None, 'n_after': None, 'count_after': None},
        'glu': {'count_OE1': 0, 'count_OE2': 0, 'count_either': 0, 'min_OE1': None, 'min_OE2': None, 'full_frames_protein_O': 0},
        # the wider count (3.5 A, seven positions): E; D; D, E and F; A and D.  The count the driver acts on (2.5 A, six
        # ring carbons): none (nearest 2.600 A); D (2.400 A); D; A and D (A at the centre of the ring).  F (2.300 A from the
        # crystal N1) is in the wider count only
        'ring': {'counts': [1] * 10 + [3] * 5 + [2] * 5, 'mean': 1.75, 'min_distance': 1.382},
        'place': {'counts': [0] * 5 + [1] * 10 + [2] * 5, 'mean': 1.0, 'min_distance': 1.382,
                  'nearest': [2.6] * 5 + [2.4] * 10 + [1.382] * 5},
        'settings': {'site_atoms': list(RING_PLACE), 'water_cutoff_A': 2.5, 'source': 'manifest hydration_check.settings',
                     'settings_are_the_defaults': True},
        'start': {'state_first_frame': 'R', 'state_differs': False, 'n2_sign_first_frame': -1, 'n2_sign_differs': False},
        'equil': {'available': True, 'source': 'manifest equilibration.stages', 'n_solute': 350, 'n_full': 35, 'total_ps': 3500.0,
                  'steps': ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq03x3', 'eq04', 'eq05', 'eq06'],
                  'frames_per_stage': [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]],
                  'time_first_ns': -3.49, 'time_last_ns': 0.0, 'time_frame_10_ns': -3.39,
                  'frame_70': [3, 1, -2.79], 'frame_170': [3, 3, -1.79], 'frame_170_step': 'eq03x3',
                  'ring_counts': [1] * 21 + [2] * 14, 'ring_last_frame_of_stage_3': 1,
                  'place_counts': [0] * 21 + [1] * 14, 'place_last_frame_of_stage_3': 0,
                  'place_nearest': [2.6] * 21 + [2.4] * 14,
                  'ring_at_the_end_of': {'eq03': 0, 'eq03x1': 0, 'eq03x2': 0, 'eq03x3': 1},
                  'driver_counts_agree': [True, True, True, True],
                  'first_not_start': None, 'first_opposite': None, 'first_sign': None, 'state_last': 'R', 'sign_last': -1},
        'warnings': 0},
    'default/R1': {
        'n2': {'count_positive': 0, 'count_negative': 200, 'count_zero': 0, 'window_frames': 200,
               'count_positive_window': 0, 'count_negative_window': 200,
               'by_state': {'C': [0, 0, 0], 'R': [200, 0, 200], 'I': [0, 0, 0]},
               'n_stretches': 1, 'n_long_stretches': 1, 'change_frames': [], 'change_times_ns': [], 'change_to_sign': [],
               'sign_first_frame': -1},
        'hb': {'first_passage_frame': None, 'n_before': None, 'count_before': None, 'n_after': None, 'count_after': None},
        'glu': {'count_OE1': 0, 'count_OE2': 0, 'count_either': 0, 'min_OE1': None, 'min_OE2': None, 'full_frames_protein_O': 0},
        'ring': {'counts': [1] * 10 + [3] * 5 + [2] * 5, 'mean': 1.75, 'min_distance': 1.382},
        'place': {'counts': [0] * 5 + [1] * 10 + [2] * 5, 'mean': 1.0, 'min_distance': 1.382,
                  'nearest': [2.6] * 5 + [2.4] * 10 + [1.382] * 5},
        'settings': {'site_atoms': list(RING_PLACE), 'water_cutoff_A': 2.5, 'source': 'specification (defaults)',
                     'settings_are_the_defaults': True},
        'start': {'state_first_frame': 'R', 'state_differs': False, 'n2_sign_first_frame': -1, 'n2_sign_differs': False},
        'equil': {'available': False},
        'warnings': 0},
    'glu/X0': {
        'n2': {'count_positive': 500, 'count_negative': 0, 'count_zero': 0, 'window_frames': 500,
               'count_positive_window': 500, 'count_negative_window': 0,
               'by_state': {'C': [500, 500, 0], 'R': [0, 0, 0], 'I': [0, 0, 0]},
               'n_stretches': 1, 'n_long_stretches': 1, 'change_frames': [], 'change_times_ns': [], 'change_to_sign': [],
               'sign_first_frame': 1},
        'hb': {'first_passage_frame': None, 'n_before': None, 'count_before': None, 'n_after': None, 'count_after': None},
        'glu': {'count_OE1': 30, 'count_OE2': 0, 'count_either': 30, 'min_OE1': 2.5, 'min_OE2': None, 'full_frames_protein_O': 0},
        'ring': {'counts': [j % 4 for j in range(50)], 'mean': 1.46, 'min_distance': 3.0},
        'start': {'state_first_frame': 'C', 'state_differs': False, 'n2_sign_first_frame': 1, 'n2_sign_differs': False},
        'equil': {'available': False},
        'warnings': 0},
}

# what the design of each run fixes by hand (not computed): section "planted by construction" of the report
DESIGN = {
    'X1': {'verdict': 'CONVERTED', 'provisional': False, 'complete': True, 'first_passage_ns': 9.01, 'n_transitions': 3,
           'count_C_window': 200, 'count_R_window': 1800, 'count_I_window': 0, 'count_C': 1050, 'count_R': 1930,
           'count_I': 20, 'window_frames': 2000, 'n_long_visits': 5, 'rule1_left_crystal': True},
    'X2': {'verdict': 'STAYED', 'provisional': False, 'complete': True, 'first_passage_ns': None, 'n_transitions': 0,
           'count_C_window': 2000, 'count_R_window': 0, 'count_I_window': 0, 'count_C': 3000, 'count_R': 0,
           'count_I': 0, 'window_frames': 2000, 'n_long_visits': 1, 'rule1_left_crystal': False},
    'X3': {'verdict': 'MIXED', 'provisional': False, 'complete': True, 'first_passage_ns': 9.01, 'n_transitions': 3,
           'count_C_window': 201, 'count_R_window': 1799, 'count_I_window': 0, 'count_C': 1051, 'count_R': 1929,
           'count_I': 20, 'window_frames': 2000, 'n_long_visits': 5, 'rule1_left_crystal': True},
    'R1': {'verdict': 'CONVERTED', 'provisional': False, 'complete': True, 'first_passage_ns': 8.01, 'n_transitions': 1,
           'count_C_window': 2000, 'count_R_window': 0, 'count_I_window': 0, 'count_C': 2200, 'count_R': 800,
           'count_I': 0, 'window_frames': 2000, 'n_long_visits': 2, 'rule1_left_crystal': False},
    'N1': {'verdict': 'STAYED', 'provisional': True, 'complete': False, 'first_passage_ns': None, 'n_transitions': 0,
           'count_C_window': 1234, 'count_R_window': 0, 'count_I_window': 0, 'count_C': 1234, 'count_R': 0,
           'count_I': 0, 'window_frames': 1234, 'n_long_visits': 1, 'rule1_left_crystal': False},
}
def _stayed(n, complete=True):
    return {'verdict': 'STAYED', 'provisional': not complete, 'complete': complete, 'first_passage_ns': None,
            'n_transitions': 0, 'count_C_window': n, 'count_R_window': 0, 'count_I_window': 0, 'count_C': n, 'count_R': 0,
            'count_I': 0, 'window_frames': n, 'n_long_visits': 1, 'rule1_left_crystal': False}


DESIGN.update({k: _stayed(GV_CASES[k]['frames']) for k in ('test/X0', 'edge/X0', 'angle/X0', 'sulfonyl/X0', 'water/X0',
                                                           'csn/X0', 'tau/X0', 'long/X0')})
DESIGN['short/X0'] = _stayed(300, complete=False)
DESIGN['gate/R1'] = {'verdict': 'MIXED', 'provisional': False, 'complete': True, 'first_passage_ns': 2.51, 'n_transitions': 1,
                     'count_C_window': 150, 'count_R_window': 250, 'count_I_window': 0, 'count_C': 150, 'count_R': 250,
                     'count_I': 0, 'window_frames': 400, 'n_long_visits': 2, 'rule1_left_crystal': False}
# version 1.2
DESIGN['conv/R1'] = {'verdict': 'MIXED', 'provisional': False, 'complete': True, 'first_passage_ns': 6.01, 'n_transitions': 1,
                     'count_C_window': 1600, 'count_R_window': 400, 'count_I_window': 0, 'count_C': 1600, 'count_R': 600,
                     'count_I': 0, 'window_frames': 2000, 'n_long_visits': 2, 'rule1_left_crystal': False}
DESIGN['flip/X1'] = {'verdict': 'MIXED', 'provisional': True, 'complete': False, 'first_passage_ns': 0.06, 'n_transitions': 2,
                     'count_C_window': 150, 'count_R_window': 145, 'count_I_window': 5, 'count_C': 150, 'count_R': 145,
                     'count_I': 5, 'window_frames': 300, 'n_long_visits': 2, 'rule1_left_crystal': False}
DESIGN['scaled/N1'] = _stayed(200)
DESIGN['ext/R1'] = {'verdict': 'STAYED', 'provisional': False, 'complete': True, 'first_passage_ns': None, 'n_transitions': 0,
                    'count_C_window': 0, 'count_R_window': 200, 'count_I_window': 0, 'count_C': 0, 'count_R': 200,
                    'count_I': 0, 'window_frames': 200, 'n_long_visits': 1, 'rule1_left_crystal': True}
DESIGN['glu/X0'] = _stayed(500)
DESIGN['place/R1'] = dict(DESIGN['ext/R1'])                 # version 1.3: the rotated start throughout, as ext/R1
DESIGN['default/R1'] = dict(DESIGN['ext/R1'])
# gate G-V by design (typed by hand): the verdict, and the criteria that are not met (names of g_gate).
# fifth: frames with an oxygen within 2.6 A, [solute count, full count]; turned: [frames with ND2 nearer in state C,
# frames of state C, the same for state R]
GV_DESIGN = {
    'test/X0': {'role': 'validation', 'verdict': 'PASSED', 'not_met': [], 'provisional': False, 'fifth': [0, 2], 'warnings': 0},
    'edge/X0': {'role': 'validation', 'verdict': 'PASSED', 'not_met': [], 'provisional': False, 'fifth': [25, 2], 'warnings': 0,
                'nd2': [480, 500, 0, 0]},
    'angle/X0': {'role': 'validation', 'verdict': 'FAILED', 'not_met': ['GV3 ang_H96NE2_Zn_H119ND1_deg'], 'provisional': False,
                 'fifth': [0, 2], 'warnings': 0, 'nd2': [500, 500, 0, 0]},
    'sulfonyl/X0': {'role': 'validation', 'verdict': 'FAILED', 'not_met': ['GV3 fifth ligand, solute trajectory'],
                    'provisional': False, 'fifth': [26, 2], 'warnings': 0, 'nd2': [500, 500, 0, 0]},
    'water/X0': {'role': 'validation', 'verdict': 'FAILED', 'not_met': ['GV3 fifth ligand, full trajectory'],
                 'provisional': False, 'fifth': [5, 3], 'warnings': 0, 'nd2': [500, 500, 0, 0]},
    'csn/X0': {'role': 'validation', 'verdict': 'FAILED', 'not_met': ['GV4 ang_C4_S_N2_deg'], 'provisional': False,
               'fifth': [0, 2], 'warnings': 0, 'nd2': [500, 500, 0, 0]},
    'tau/X0': {'role': 'validation', 'verdict': 'FAILED', 'not_met': ['GV4 tau_fold_deg'], 'provisional': False,
               'fifth': [0, 2], 'warnings': 0, 'nd2': [500, 500, 0, 0]},
    'long/X0': {'role': 'validation', 'verdict': 'PASSED', 'not_met': [], 'provisional': False, 'fifth': [0, 3], 'warnings': 2,
                'nd2': [600, 600, 0, 0]},
    'short/X0': {'role': 'validation', 'verdict': 'PASSED', 'not_met': [], 'provisional': True, 'fifth': [0, 1], 'warnings': 0,
                 'nd2': [300, 300, 0, 0]},
    # a production run: no verdict; which criteria its numbers meet is left to g_gate (not typed by hand)
    'gate/R1': {'role': 'production', 'verdict': None, 'not_met': None, 'provisional': False, 'fifth': [0, 1],
                'warnings': 0, 'nd2': [113, 150, 63, 250]},
    # version 1.2: OE1 of Glu106 within 2.6 A in 30 of 500 solute frames (6 %), in no full frame: the gate passes
    'glu/X0': {'role': 'validation', 'verdict': 'PASSED', 'not_met': [], 'provisional': False, 'fifth': [0, 2], 'warnings': 0,
               'nd2': [500, 500, 0, 0]},
}
# the limits that the bent runs are meant to sit beside: (run, column, lowest and highest deviation allowed by design)
GV_BESIDE = (('edge/X0', 'ang_H96NE2_Zn_H119ND1_deg', 8.0, 10.0), ('edge/X0', 'ang_C4_S_N2_deg', 2.0, 3.0),
             ('edge/X0', 'tau_fold_deg', 15.0, 20.0), ('angle/X0', 'ang_H96NE2_Zn_H119ND1_deg', 10.0, 12.0),
             ('csn/X0', 'ang_C4_S_N2_deg', 3.0, 4.0), ('tau/X0', 'tau_fold_deg', 20.0, 25.0))
CAMPAIGN_DESIGN = {'n_runs': 5, 'n_complete': 4, 'all_nine_complete': False, 'final': False,
                   'rule1.n_X_runs': 3, 'rule1.n_X_left_crystal': 2, 'rule1.failed': True,
                   'n_converted.X': 1, 'n_converted.N': 0, 'n_converted.R': 1, 'n_stayed.X': 1, 'n_stayed.N': 1,
                   'n_stayed.R': 0, 'n_mixed.X': 1, 'outcome': 'validation_failed',
                   'rule4_sensitivity_run_required': False}

ASN67_COLS = (('d_N1_N67ND2_A', 'A'), ('asn67_ND2_nearer_N1', 'int'))        # only with the key asn67_ND2
SOLUTE_COLS = (
    ('lig_rmsd_A', 'A'), ('lig_rmsd_labelling', 'int'), ('lig_rmsd_identity_A', 'A'),
    ('d1_N1_F131CZ_A', 'A'), ('d2_N1_N67OD1_A', 'A'), ('state', 'state'),
    ('d_N2_T199OG1_A', 'A'), ('ang_N2_HN2_T199OG1_deg', 'deg'), ('hb_N2H_T199OG1', 'int'),
    ('d_Zn_N2_A', 'A'), ('d_Zn_H94NE2_A', 'A'), ('d_Zn_H96NE2_A', 'A'), ('d_Zn_H119ND1_A', 'A'),
    ('ang_N2_Zn_H94NE2_deg', 'deg'), ('ang_N2_Zn_H96NE2_deg', 'deg'), ('ang_N2_Zn_H119ND1_deg', 'deg'),
    ('ang_H94NE2_Zn_H96NE2_deg', 'deg'), ('ang_H94NE2_Zn_H119ND1_deg', 'deg'), ('ang_H96NE2_Zn_H119ND1_deg', 'deg'),
    ('ang_S_N2_Zn_deg', 'deg'), ('d_Zn_O1_A', 'A'), ('d_Zn_O2_A', 'A'), ('d_Zn_Osulfonyl_min_A', 'A'),
    ('fifth_sulfonylO', 'int'), ('d_T199N_O1_A', 'A'), ('d_T199N_O2_A', 'A'), ('ang_T199N_H_O1_deg', 'deg'),
    ('ang_T199N_H_O2_deg', 'deg'), ('d_T199N_Osulfonyl_min_A', 'A'), ('hb_T199NH_Osulfonyl', 'int'),
    ('d_T199OG1_E106OE1_A', 'A'), ('d_T199OG1_E106OE2_A', 'A'), ('ang_T199OG1_HG1_E106OE1_deg', 'deg'),
    ('ang_T199OG1_HG1_E106OE2_deg', 'deg'), ('d_T199OG1_E106OE_min_A', 'A'), ('hb_T199OG1H_E106OE', 'int'),
    ('tors_C3_C4_S_N2_deg', 'tors'), ('tors_C5_C4_S_N2_deg', 'tors'), ('tau_fold_deg', 'deg'),
    ('ang_C4_S_N2_deg', 'deg'), ('d_S_N2_A', 'A'), ('ca_rmsd_A', 'A'),
    # version 1.2 (n2_sign_assigned is planted too, by g_stretches; it is not written for the equilibration)
    ('n2_volume', 'vol'), ('n2_sign', 'int'), ('d_Zn_E106OE1_A', 'A'), ('d_Zn_E106OE2_A', 'A'),
    ('d_Zn_E106OE_min_A', 'A'), ('fifth_E106OE1', 'int'), ('fifth_E106OE2', 'int'), ('fifth_E106OE', 'int'))
FULL_COLS = (('d_Zn_Owater_min_A', 'A'), ('water_O_index_nearest_Zn', 'int'), ('d_Zn_Oprotein_min_A', 'A'),
             ('protein_O_index_nearest_Zn', 'int'), ('d_Zn_Osulfonyl_min_full_A', 'A'), ('d_Zn_O_any_min_A', 'A'),
             ('fifth_water', 'int'), ('fifth_proteinO', 'int'), ('fifth_any', 'int'), ('n_water_N1_3p5A', 'int'),
             ('d_Zn_N2_full_A', 'A'), ('n_water_ring_wide_3p5A', 'int'), ('d_water_ring_wide_min_A', 'A'),
             # version 1.3: the count the driver acts on (the two columns before it were named ..._ring_site_... in 1.2)
             ('n_water_ring_place', 'int'), ('d_water_ring_place_min_A', 'A'))
STARTS_MAP = (('rmsd_to_crystal_A', 'lig_rmsd_A', 'A'), ('Zn_N2', 'd_Zn_N2_A', 'A'),
              ('Zn_O_min', 'd_Zn_Osulfonyl_min_A', 'A'), ('N2_OG1', 'd_N2_T199OG1_A', 'A'),
              ('N2_H_OG1', 'ang_N2_HN2_T199OG1_deg', 'deg'), ('T199N_O_min', 'd_T199N_Osulfonyl_min_A', 'A'),
              ('C3_C4_S_N2', 'tors_C3_C4_S_N2_deg', 'deg'), ('C5_C4_S_N2', 'tors_C5_C4_S_N2_deg', 'deg'),
              ('tau_fold', 'tau_fold_deg', 'deg'), ('C4_S_N2', 'ang_C4_S_N2_deg', 'deg'), ('S_N2', 'd_S_N2_A', 'A'),
              ('S_N2_Zn', 'ang_S_N2_Zn_deg', 'deg'), ('N1_Phe131CZ', 'd1_N1_F131CZ_A', 'A'),
              ('N1_Asn67OD1', 'd2_N1_N67OD1_A', 'A'))


# ----------------------------------------------------------------------------------------------------------
# one frame: coordinates in the 6RL9 frame and the planted values
# ----------------------------------------------------------------------------------------------------------
def frame_coordinates(sy, n, pose, swap, opt=None):
    """Coordinates [n_atoms, 3] of solute frame n before rotation, and the dictionaries used for the values.
    opt (runs of version 1.1): tau, csn (ring turned, degrees), nd1 (ND1 of His119 turned about the zinc, degrees),
    sulfonyl (frames with the zinc 2.500 A from the nearer sulfonyl oxygen), amide_turned (frames with OD1 and ND2
    of Asn67 exchanged)."""
    opt = opt or {}
    x = sy.x0.copy()
    s = 1.0 if n % 2 == 0 else -1.0
    x[sy.ca] = sy.ca_centroid + (1.0 + EPS * s) * (sy.ca_ref - sy.ca_centroid)
    base = sy.pose[pose]
    if opt.get('tau') or opt.get('csn'):
        base = sy.bent_pose(pose, opt.get('tau', 0.0), opt.get('csn', 0.0))
    L = dict(base)
    if 'ring' in swap:
        L = {k: base[RING_SWAP.get(k, k)] for k in L}
    if 'oxy' in swap:
        L = {k: L[OXY_SWAP.get(k, k)] for k in L}
    L = {k: np.array(v, float) for k, v in L.items()}
    S = dict(sy.S)
    if opt.get('nd1'):
        S['his119_ND1'] = sy.nd1_turned(opt['nd1'])
    if n in opt.get('amide_turned', ()):
        S['asn67_OD1'], S['asn67_ND2'] = sy.S['asn67_ND2'], sy.S['asn67_OD1']
    for k, i in sy.site.items():
        x[i] = S[k]
    og = S['thr199_OG1']
    if n % 5 == 4 and not opt.get('v12'):                     # hydrogen turned away from Thr199 OG1
        d = L['N2'] - og
        L['HN21'] = L['N2'] + 1.021 * d / np.linalg.norm(d)
    if n in opt.get('mirror', ()):                            # version 1.2: hydrogen mirrored through the plane Zn-N2-S
        m = np.cross(sy.zn0 - L['N2'], L['S'] - L['N2'])
        m /= np.linalg.norm(m)
        L['HN21'] = L['HN21'] - 2.0 * float((L['HN21'] - L['N2']) @ m) * m
    if n in opt.get('og1_away', ()):                          # version 1.2: OG1 of Thr199 moved away from N2
        d = og - L['N2']
        S['thr199_OG1'] = og + D_OG1_AWAY * d / np.linalg.norm(d)
        x[sy.site['thr199_OG1']] = S['thr199_OG1']
    zn = sy.zn0.copy()
    if n % 10 == 0:                                           # zinc moved away from N2
        d = zn - L['N2']
        zn = zn + DZN * d / np.linalg.norm(d)
    if n in opt.get('sulfonyl', ()):                          # zinc put 2.500 A from the nearer sulfonyl oxygen
        o = min((L['O1'], L['O2']), key=lambda a: g_dist(zn, a))
        d = zn - o
        zn = o + DS_NEAR * d / np.linalg.norm(d)
    for key, atom, dz in (('oe1_near', 'glu106_OE1', D_OE1_NEAR), ('oe1_far', 'glu106_OE1', D_OE1_FAR),
                          ('oe2_near', 'glu106_OE2', D_OE2_NEAR)):       # version 1.2: Glu106 at a known distance
        if n in opt.get(key, ()):
            d = sy.S[atom] - zn
            S[atom] = zn + dz * d / np.linalg.norm(d)
            x[sy.site[atom]] = S[atom]
    x[sy.izn] = zn
    for k in LIG_ALL:
        x[sy.ilig[k]] = L[k]
    return x, L, zn, S


def solute_values(sy, n, L, zn, S=None):
    S, X = (sy.S if S is None else S), sy.pose['X']
    v = {}
    r = [g_rmsd(L, X, p) for p in LABELLINGS]
    v['lig_rmsd_labelling'] = int(min(range(4), key=lambda i: (r[i], i)))
    v['lig_rmsd_A'], v['lig_rmsd_identity_A'] = r[v['lig_rmsd_labelling']], r[0]
    d1, d2 = g_dist(L['N1'], S['phe131_CZ']), g_dist(L['N1'], S['asn67_OD1'])
    v['d1_N1_F131CZ_A'], v['d2_N1_N67OD1_A'] = d1, d2
    v['state'] = 'C' if (d1 < 6.7 and d2 > 5.2) else 'R' if (d1 > 6.7 and d2 < 5.2) else 'I'
    if 'asn67_ND2' in S:
        v['d_N1_N67ND2_A'] = g_dist(L['N1'], S['asn67_ND2'])
        v['asn67_ND2_nearer_N1'] = int(v['d_N1_N67ND2_A'] < d2)
    v['d_N2_T199OG1_A'] = g_dist(L['N2'], S['thr199_OG1'])
    v['ang_N2_HN2_T199OG1_deg'] = g_angle(L['N2'], L['HN21'], S['thr199_OG1'])
    v['hb_N2H_T199OG1'] = int(v['d_N2_T199OG1_A'] <= 3.5 and v['ang_N2_HN2_T199OG1_deg'] >= 135.0)
    four = (('N2', L['N2']), ('H94NE2', S['his94_NE2']), ('H96NE2', S['his96_NE2']), ('H119ND1', S['his119_ND1']))
    for a, xa in four:
        v[f'd_Zn_{a}_A'] = g_dist(zn, xa)
    for i in range(4):
        for j in range(i + 1, 4):
            v[f'ang_{four[i][0]}_Zn_{four[j][0]}_deg'] = g_angle(four[i][1], zn, four[j][1])
    v['ang_S_N2_Zn_deg'] = g_angle(L['S'], L['N2'], zn)
    v['d_Zn_O1_A'], v['d_Zn_O2_A'] = g_dist(zn, L['O1']), g_dist(zn, L['O2'])
    v['d_Zn_Osulfonyl_min_A'] = min(v['d_Zn_O1_A'], v['d_Zn_O2_A'])
    v['fifth_sulfonylO'] = int(v['d_Zn_Osulfonyl_min_A'] <= 2.6)
    for o in ('O1', 'O2'):
        v[f'd_T199N_{o}_A'] = g_dist(S['thr199_N'], L[o])
        v[f'ang_T199N_H_{o}_deg'] = g_angle(S['thr199_N'], S['thr199_H'], L[o])
    v['d_T199N_Osulfonyl_min_A'] = min(v['d_T199N_O1_A'], v['d_T199N_O2_A'])
    v['hb_T199NH_Osulfonyl'] = int(any(v[f'd_T199N_{o}_A'] <= 3.5 and v[f'ang_T199N_H_{o}_deg'] >= 135.0 for o in ('O1', 'O2')))
    for o in ('OE1', 'OE2'):
        v[f'd_T199OG1_E106{o}_A'] = g_dist(S['thr199_OG1'], S[f'glu106_{o}'])
        v[f'ang_T199OG1_HG1_E106{o}_deg'] = g_angle(S['thr199_OG1'], S['thr199_HG1'], S[f'glu106_{o}'])
    v['d_T199OG1_E106OE_min_A'] = min(v['d_T199OG1_E106OE1_A'], v['d_T199OG1_E106OE2_A'])
    v['hb_T199OG1H_E106OE'] = int(any(v[f'd_T199OG1_E106{o}_A'] <= 3.5 and v[f'ang_T199OG1_HG1_E106{o}_deg'] >= 135.0
                                      for o in ('OE1', 'OE2')))
    v['tors_C3_C4_S_N2_deg'] = g_torsion(L['C3'], L['C4'], L['S'], L['N2'])
    v['tors_C5_C4_S_N2_deg'] = g_torsion(L['C5'], L['C4'], L['S'], L['N2'])
    v['tau_fold_deg'] = max(abs(v['tors_C3_C4_S_N2_deg']), abs(v['tors_C5_C4_S_N2_deg']))
    v['ang_C4_S_N2_deg'] = g_angle(L['C4'], L['S'], L['N2'])
    v['d_S_N2_A'] = g_dist(L['S'], L['N2'])
    rg = math.sqrt(float(((sy.ca_ref - sy.ca_centroid) ** 2).sum(1).mean()))
    v['ca_rmsd_A'] = EPS * rg                                 # by construction
    # version 1.2
    v['n2_volume'] = g_volume(zn, L['N2'], L['HN21'], L['S'])
    v['n2_sign'] = g_sign(v['n2_volume'])
    v['d_Zn_E106OE1_A'], v['d_Zn_E106OE2_A'] = g_dist(zn, S['glu106_OE1']), g_dist(zn, S['glu106_OE2'])
    v['d_Zn_E106OE_min_A'] = min(v['d_Zn_E106OE1_A'], v['d_Zn_E106OE2_A'])
    v['fifth_E106OE1'], v['fifth_E106OE2'] = int(v['d_Zn_E106OE1_A'] <= 2.6), int(v['d_Zn_E106OE2_A'] <= 2.6)
    v['fifth_E106OE'] = int(v['d_Zn_E106OE_min_A'] <= 2.6)
    return v


def place_water(sy, x, j, L, zn, near=None, ring=None, place=None):
    """Water and ions of full frame j (6RL9 frame).  Returns the planted values of the frame.
    near: the full frames with the water 2.400 A from the zinc (default: every tenth, j = 3, 13, ...).
    ring (version 1.2): letters of the water planted at the crystal ring in this frame (A, B, C; D, E, F since version
    1.3; see V12_CASES).  place (version 1.3): (site atoms, cut-off) of the count the driver acts on; default: the six
    ring carbons and 2.5 A."""
    n1 = L['N1']
    d0 = (n1 - zn) / np.linalg.norm(n1 - zn)
    p = np.cross(d0, [0.0, 0.0, 1.0])
    p /= np.linalg.norm(p)
    q = np.cross(d0, p)
    c, s = math.cos(math.radians(50.0)), math.sin(math.radians(50.0))
    dirs = (d0, c * d0 + s * p, c * d0 + s * q)
    dw = DW_NEAR if (j % 10 == 3 if near is None else j in near) else DW_FAR
    k_in = j % 4
    pos = [zn + dw * sy.u_zn] + [n1 + (RHO_IN if k < k_in else RHO_OUT) * dirs[k] for k in range(3)]
    m = 36                                                    # the other water molecules: 15 A from the zinc
    for k in range(m):
        z = 1.0 - 2.0 * (k + 0.5) / m
        r = math.sqrt(1.0 - z * z)
        ph = k * math.pi * (3.0 - math.sqrt(5.0))
        pos.append(sy.zn0 + 15.0 * np.array([r * math.cos(ph), r * math.sin(ph), z]))
    X = sy.pose['X']
    if ring:                                                  # version 1.2: three of the far molecules are moved
        cen = sum(X[a] for a in RING_PLACE) / 6.0
        nrm = np.cross(X['C1'] - cen, X['C2'] - cen)
        nrm /= np.linalg.norm(nrm)
        out = (X['N1'] - cen) / np.linalg.norm(X['N1'] - cen)
        spot = {'A': cen, 'B': X['N1'] + D_RING_B * out, 'C': X['C3'] + D_RING_C * nrm,
                'D': X['C5'] + D_RING_D * nrm, 'E': X['C2'] - D_RING_E * nrm, 'F': X['N1'] + D_RING_F * out}
        for k, letter in enumerate('ABCDEF'):
            if letter in ring:
                pos[4 + k] = spot[letter]
    for w, io in enumerate(sy.water_o):
        x[io] = pos[w]
        x[io + 1] = pos[w] + (0.757, 0.586, 0.0)
        x[io + 2] = pos[w] + (-0.757, 0.586, 0.0)
    for k, i in enumerate(sy.ions):
        x[i] = sy.zn0 + (25.0, 5.0 * (k + 1), -20.0)
    # values, by brute force over all atoms concerned (second formulas, no periodic images needed here)
    dz = [g_dist(zn, x[i]) for i in sy.water_o]
    kw = int(np.argmin(dz))
    dp = [g_dist(zn, x[i]) for i in sy.protein_o]
    kp = int(np.argmin(dp))
    ds = min(g_dist(zn, L['O1']), g_dist(zn, L['O2']))
    nn = sum(1 for i in sy.water_o if g_dist(n1, x[i]) <= 3.5)
    assert abs(dz[kw] - dw) < 1e-9 and kw == 0, 'the planted water is not the nearest to the zinc'
    assert nn == k_in, f'planted {k_in} water molecules at the aniline N, found {nn}'
    # version 1.2: water oxygens within 3.5 A of a crystal position of C1 to C6 or N1, each molecule once (the wider
    # count); version 1.3: those within the cut-off of a crystal position of the site atoms (the count of the driver)
    dr = [min(g_dist(x[i], X[a]) for a in RING_WIDE) for i in sy.water_o]
    atoms, cut = place or (RING_PLACE, CUT_PLACE)
    dpl = [min(g_dist(x[i], X[a]) for a in atoms) for i in sy.water_o]
    v12 = {'n_water_ring_wide_3p5A': sum(1 for d in dr if d <= CUT_WIDE), 'd_water_ring_wide_min_A': min(dr),
           'n_water_ring_place': sum(1 for d in dpl if d <= cut), 'd_water_ring_place_min_A': min(dpl)}
    return {**v12, 'd_Zn_Owater_min_A': dw, 'water_O_index_nearest_Zn': sy.water_o[kw], 'd_Zn_Oprotein_min_A': dp[kp],
            'protein_O_index_nearest_Zn': sy.protein_o[kp], 'd_Zn_Osulfonyl_min_full_A': ds,
            'd_Zn_O_any_min_A': min(dw, dp[kp], ds), 'fifth_water': int(dw <= 2.6), 'fifth_proteinO': int(dp[kp] <= 2.6),
            'fifth_any': int(min(dw, dp[kp], ds) <= 2.6), 'n_water_N1_3p5A': k_in, 'd_Zn_N2_full_A': g_dist(zn, L['N2'])}


def to_lab(sy, x, n, rng, B, with_solvent):
    """Rotation about the zinc, shift, lattice moves and wrapping: what the trajectory file holds."""
    Q = random_rotation(rng)
    t = rng.uniform(-8.0, 8.0, size=3)
    y = (x - sy.zn0) @ Q.T + t
    if n % 7 == 3:                                            # the ligand in another periodic image
        y[[sy.ilig[k] for k in LIG_ALL]] += B[0] - B[2]
    if n % 11 == 5:                                           # the zinc in another periodic image
        y[sy.izn] += B[1]
    if with_solvent:
        for io in sy.water_o:                                 # molecules wrapped into the cell by their oxygen
            f = np.linalg.solve(B.T, y[io])
            y[io:io + 3] -= np.floor(f) @ B
        for i in sy.ions:
            f = np.linalg.solve(B.T, y[i])
            y[i] -= np.floor(f) @ B
    return y


# ----------------------------------------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------------------------------------
def write_pdb(path, top, xyz):
    import openmm.app as app
    import openmm.unit as u
    with open(path, 'w') as f:
        app.PDBFile.writeFile(top, xyz * u.angstrom, f, keepIds=True)


def stats(x):
    x = np.asarray(x, float)
    n = len(x)
    last = x[n - n // 2:]
    return {'mean': float(x.mean()), 'sd': float(x.std()), 'mean_last_half': float(last.mean()),
            'sd_last_half': float(last.std()), 'block_means': [float(x[k * 1000:(k + 1) * 1000].mean()) for k in range(n // 1000)]}


def equil_steps(eq):
    """The steps of an equilibration as (tag, start, end in ps): the stages eq01.., or the list `steps` of a case
    with extension steps (eq03x1, ..)."""
    return list(eq['steps']) if 'steps' in eq else [(f'eq{k:02d}', a, b) for k, (a, b) in enumerate(eq['stages'], 1)]


def generate_equilibration(sy, key, rdir, eq, top_s, top_f, pdir, tag, start, place=None):
    """Version 1.2: the files solute_eqNN.dcd and full_eqNN.dcd of a run and their planted values.  A stage from
    a to b ps holds the frames at a + 10, a + 20, ... ps (solute) and a + 100, a + 200, ... ps (full) that are not
    later than b; the time written is counted from the end of the last stage, in ns."""
    import openmm.app as app
    import openmm.unit as u
    from openmm import Vec3
    rng = np.random.default_rng([SEED, sum(ord(c) for c in key), 12])
    steps, opt = equil_steps(eq), eq['opt']
    t_end = steps[-1][2]
    dts, dtf = STEP_PS * SOL_STEPS, STEP_PS * FULL_STEPS
    srows, frows, per, e, j = [], [], [], 0, 0
    for tag_k, a, b in steps:
        k, ext = (int(v) for v in (tag_k[2:].split('x') + ['0'])[:2])           # stage number, extension number (0: the stage)
        ns, nfull = int(math.floor((b - a) / dts + 1e-9)), int(math.floor((b - a) / dtf + 1e-9))
        lab_s, lab_f, e0, j0 = [], [], e, j
        for i in range(ns):
            pose, swap = eq['pose'](e)
            x, L, zn, S = frame_coordinates(sy, e, pose, swap, opt)
            v = solute_values(sy, e, L, zn, S)
            v.update(frame=e, stage=k, extension=ext, time_ns=(a + (i + 1) * dts - t_end) / 1000.0, pose=pose)
            srows.append(v)
            B = cell(BOX_L + BOX_AMP * math.sin(2 * math.pi * e / 500.0))
            full = (i + 1) % 10 == 0 and (i + 1) // 10 <= nfull
            if full:
                w = place_water(sy, x, j, L, zn, opt.get('water_near'), opt.get('ring', {}).get(j), place)
                w.update(frame=j, stage=k, extension=ext, time_ns=(a + ((i + 1) // 10) * dtf - t_end) / 1000.0)
                frows.append(w)
                j += 1
            y = to_lab(sy, x, e, rng, B, full)
            lab_s.append((y[sy.solute], B))
            if full:
                lab_f.append((y, B))
            e += 1
        for name, top, items, interval in ((f'solute_{tag_k}.dcd', top_s, lab_s, SOL_STEPS), (f'full_{tag_k}.dcd', top_f, lab_f, FULL_STEPS)):
            # header: steps since the start of the equilibration (md/run_md.py), or one interval in every file
            first = interval if eq['header'] == 'reporter' else int(round(a / STEP_PS)) + interval
            with open(rdir / name, 'wb') as f:
                d = app.DCDFile(f, top, STEP_PS * u.picosecond, first, interval)
                for y, B in items:
                    d.writeModel(y * u.angstrom,
                                 periodicBoxVectors=u.Quantity([Vec3(*(B[m] / 10.0)) for m in range(3)], u.nanometer))
        rs, rf = srows[e0:e], frows[j0:j]
        dc = eq.get('driver_counts', {}).get(tag_k)
        if dc is not None:                                    # what the manifest says the driver counted, against the frame as built
            assert rf and rf[-1]['n_water_ring_place'] == dc[0] and abs(rf[-1]['d_water_ring_place_min_A'] - dc[1]) < 5e-4, (key, tag_k, rf[-1])
        per.append({'stage': k, 'extension': ext, 'tag': tag_k, 'driver_count': None if dc is None else dc[0],
                    'time_start_ps': a, 'time_end_ps': b, 'n_frames_solute': ns, 'n_frames_full': nfull,
                    'count_C': sum(r['state'] == 'C' for r in rs), 'count_R': sum(r['state'] == 'R' for r in rs),
                    'count_I': sum(r['state'] == 'I' for r in rs), 'count_n2_positive': sum(r['n2_sign'] > 0 for r in rs),
                    'count_n2_negative': sum(r['n2_sign'] < 0 for r in rs),
                    'last_solute_frame': None if not rs else {'frame': rs[-1]['frame'], 'time_ns': rs[-1]['time_ns'], 'state': rs[-1]['state'],
                                                              'lig_rmsd_A': rs[-1]['lig_rmsd_A'], 'n2_volume': rs[-1]['n2_volume'], 'n2_sign': rs[-1]['n2_sign']},
                    'last_full_frame': None if not rf else {'frame': rf[-1]['frame'], 'time_ns': rf[-1]['time_ns'],
                                                            'n_water_ring_wide_3p5A': rf[-1]['n_water_ring_wide_3p5A'],
                                                            'd_water_ring_wide_min_A': rf[-1]['d_water_ring_wide_min_A'],
                                                            'n_water_ring_place': rf[-1]['n_water_ring_place'],
                                                            'd_water_ring_place_min_A': rf[-1]['d_water_ring_place_min_A']}})
    scols = SOLUTE_COLS + (ASN67_COLS if sy.with_nd2 else ())
    with open(pdir / f'{tag}__equil_solute.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['frame', 'stage', 'extension', 'time_ns', 'pose'] + [c for c, _ in scols])
        for r in srows:
            w.writerow([r['frame'], r['stage'], r['extension'], '%.6f' % r['time_ns'], r['pose']]
                       + [r[c] if isinstance(r[c], (int, str)) else '%.6f' % r[c] for c, _ in scols])
    with open(pdir / f'{tag}__equil_full.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['frame', 'stage', 'extension', 'time_ns'] + [c for c, _ in FULL_COLS])
        for r in frows:
            w.writerow([r['frame'], r['stage'], r['extension'], '%.6f' % r['time_ns']] + [r[c] if isinstance(r[c], int) else '%.6f' % r[c] for c, _ in FULL_COLS])
    s0, g0 = STATE_OF_START[start], SIGN_OF_START[start]
    opp = 'R' if s0 == 'C' else 'C'

    def first(test):
        r = next((r for r in srows if test(r)), None)
        return None if r is None else {'frame': r['frame'], 'stage': r['stage'], 'extension': r['extension'], 'time_ns': r['time_ns']}
    return {'n_frames_solute': len(srows), 'n_frames_full': len(frows), 'total_ps': t_end, 'stages': per,
            'frames_70_170': [[srows[n]['stage'], srows[n]['extension'], srows[n]['time_ns']] for n in (70, 170)] if len(srows) > 170 else None,
            'in_manifest': eq['in_manifest'], 'header': eq['header'],
            'frames_solute_moved': {'ligand': sum(1 for r in srows if r['frame'] % 7 == 3), 'zinc': sum(1 for r in srows if r['frame'] % 11 == 5)},
            'count_C': sum(r['state'] == 'C' for r in srows), 'count_R': sum(r['state'] == 'R' for r in srows),
            'count_I': sum(r['state'] == 'I' for r in srows),
            'count_n2_positive': sum(r['n2_sign'] > 0 for r in srows), 'count_n2_negative': sum(r['n2_sign'] < 0 for r in srows),
            'first_frame_not_in_state_of_start': first(lambda r: r['state'] != s0),
            'first_frame_in_opposite_state': first(lambda r: r['state'] == opp),
            'first_frame_n2_sign_not_of_start': first(lambda r: r['n2_sign'] != g0),
            'state_last_frame': srows[-1]['state'] if srows else None, 'n2_sign_last_frame': srows[-1]['n2_sign'] if srows else None,
            'time_first_ns': srows[0]['time_ns'] if srows else None, 'time_last_ns': srows[-1]['time_ns'] if srows else None,
            'times_ns': [r['time_ns'] for r in srows] if len(srows) <= 10 else None,
            'ring_counts': [r['n_water_ring_wide_3p5A'] for r in frows],
            'place_counts': [r['n_water_ring_place'] for r in frows],
            'place_nearest': [r['d_water_ring_place_min_A'] for r in frows]}


def generate_run(sy, run, root, model='test', runs_folder='runs'):
    import openmm.app as app
    import openmm.unit as u
    from openmm import Vec3
    key = f'{model}/{run}'                                    # the runs of versions 1.1 and 1.2 are named <model>/<run>
    new = key in GV_CASES or key in V12_CASES
    case12 = V12_CASES.get(key)
    frames, info = schedule(key if new else run)
    opt = info.get('opt', {})
    near = opt.get('water_near')
    ring = opt.get('ring', {})
    hyd = case12.get('hydration') if case12 else None        # version 1.3: the settings of the hydration check in the manifest
    place = (tuple(hyd['site_atoms']), hyd['water_cutoff_A']) if hyd else None
    nfr = len(frames)
    rdir = root / runs_folder / model / run
    rdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([SEED, sum(ord(c) for c in (key if new else run))])
    top_s, top_f = sy.topology(sy.solute, BOX_L), sy.topology(range(sy.n), BOX_L)
    x_first, L_first, zn_first, _ = frame_coordinates(sy, 0, *frames[0], opt)
    if new:
        place_water(sy, x_first, 0, L_first, zn_first, near, ring.get(0), place)
    else:                                                     # as in version 1.0, so that the old runs keep their files
        place_water(sy, x_first, 0, sy.pose[frames[0][0]], sy.zn0)
    write_pdb(rdir / 'solute.pdb', top_s, x_first[sy.solute])
    write_pdb(rdir / 'full.pdb', top_f, x_first)
    seg = info['segments']
    rep = info.get('repeat_from')
    kind = info.get('manifest', 'plain' if run in ('X2', 'X3') else 'driver')    # the forms of manifest that are read
    plain = kind != 'driver'
    if kind == 'none':                                        # no manifest: folder name, convention, defaults
        man = None
    elif kind == 'plain-bare':                                # neither replica nor length nor system folder
        man = {'model': model, 'start': info['start']}
    elif plain:                                               # every plain key
        man = {'model': model, 'start': info['start'], 'replica': int(run[1:]), 'planned_ns': info['planned_ns'],
               'system_dir': f'../../../systems/{model}/{info["start"]}', 'solute_interval_ps': STEP_PS * SOL_STEPS,
               'full_interval_ps': STEP_PS * FULL_STEPS}
    else:                                                     # the keys of the manifest of md/run_md.py that are read
        recs, t0 = [], 0.0
        done = seg[:-1] if info.get('broken_last_frame') else seg      # a segment in progress has no record yet
        for k, m in enumerate(done, 1):
            if rep is not None and k == 2:
                t0 = rep * STEP_PS * SOL_STEPS                # the restart that overlaps
            recs.append({'tag': f'seg{k:04d}', 'kind': 'production', 'time_start_ps': t0,
                         'time_end_ps': t0 + m * STEP_PS * SOL_STEPS, 'frames': {'solute': m, 'full': m // 10}})
            t0 += m * STEP_PS * SOL_STEPS
        sysd = (root / 'systems' / model / info['start']).resolve()
        hashes = {f: {'sha256': sha256(sysd / f)} for f in ('selections.json', 'reference.pdb')}
        if rep is not None:                                   # the run that is wrong in every way: another system
            hashes['selections.json']['sha256'] = '0' * 64
        man = {'schema': 'cqc-sfa-md run manifest, version 1 (the keys read by analyse_md.py only)', 'run_id': run,
               # R1: as after the transfer from another machine, the folder named in the manifest is not here
               'system_dir': '/not/on/this/machine/systems/test/R' if (run == 'R1' and not new) else str(sysd),
               'inputs': hashes,
               'build_json': {'model': model, 'start': info['start']},
               'protocol': {'production': {'length_ns': 100.0},
                            'output': {'solute_frame_interval_ps': STEP_PS * SOL_STEPS,
                                       'full_frame_interval_ps': STEP_PS * FULL_STEPS}},
               'production': {'target_ns': info.get('target_ns') or info['planned_ns'], 'segments': recs}}
        if case12 and case12.get('equil', {}).get('in_manifest'):          # version 1.2: the stages as the driver lists them
            dc = case12['equil'].get('driver_counts', {})
            full_form = bool(hyd) and hyd['form'] == 'driver'  # version 1.3: the record as md/run_md.py writes it now

            def hyd_record(t, n_done):
                h = {'n_water_O': dc[t][0], 'nearest_water_O': {'distance_A': dc[t][1]}, 'applies': True, 'decision': dc[t][2]}
                if full_form:
                    nxt = f'eq03x{n_done + 1}' if dc[t][2] == 'extend' else 'eq04'
                    h = {'made_on': f'the coordinates at the end of this step (state_{t}.xml)', 'site_atoms': list(hyd['site_atoms']),
                         'ligand_mean_distance_from_site_A': 9.9, 'applies_if_ligand_mean_distance_above_A': 2.0,
                         'applies': True, 'skipped': False, 'water_cutoff_A': hyd['water_cutoff_A'], 'n_water_O': dc[t][0],
                         'water_O': [{'atom_index': 0, 'atom': 'WAT:O', 'distance_A': dc[t][1]}] * dc[t][0],
                         'nearest_water_O': {'atom_index': 0, 'atom': 'WAT:O', 'distance_A': dc[t][1]},
                         'n_extensions_done': n_done, 'max_extensions': 4, 'decision': dc[t][2], 'next_step': nxt}
                return h
            man['equilibration'] = {'stages': [
                {'tag': t, 'kind': 'equilibration', 'time_start_ps': a, 'time_end_ps': b,
                 'frames': {'solute': int(math.floor((b - a) / (STEP_PS * SOL_STEPS) + 1e-9)),
                            'full': int(math.floor((b - a) / (STEP_PS * FULL_STEPS) + 1e-9))},
                 **({'hydration_check': hyd_record(t, int((t.split('x') + ['0'])[1]))} if t in dc else {})}
                for t, a, b in equil_steps(case12['equil'])]}
        if hyd:                                               # version 1.3: the key hydration_check of the manifest
            dc = case12.get('equil', {}).get('driver_counts', {})
            ends = {t: b for t, a, b in equil_steps(case12['equil'])} if case12.get('equil') else {}
            man['hydration_check'] = {
                'requested': True,
                'settings': {'after_step': 'eq03', 'site_atoms': list(hyd['site_atoms']), 'water_cutoff_A': hyd['water_cutoff_A'],
                             'applies_if_ligand_mean_distance_above_A': 2.0, 'max_extensions': 4,
                             'extension_length_ps_of_the_protocol': 500.0, 'extension_length_ps_in_this_run': 500.0},
                'status': 'decided', 'applies': True, 'skipped': False,
                'counts': [{'step': t, 'time_end_ps': ends[t], 'n_water_O': v[0], 'nearest_water_O_A': v[1], 'applies': True,
                            'decision': v[2]} for t, v in dc.items()],
                'n_extensions': sum(1 for t in dc if 'x' in t), 'extension_ps': 500.0 * sum(1 for t in dc if 'x' in t)}
    if man is not None:
        man['note'] = 'synthetic test run of md/analysis/make_test_trajectory.py'
        (rdir / 'manifest.json').write_text(json.dumps(man, indent=1) + '\n')

    srows, frows, lab_s, lab_f = [], [], [], []
    for n in range(nfr):
        x, L, zn, S = frame_coordinates(sy, n, *frames[n], opt)
        srows.append(solute_values(sy, n, L, zn, S))
        Lbox = BOX_L + BOX_AMP * math.sin(2 * math.pi * n / 500.0)
        B = cell(Lbox)
        full = (n + 1) % 10 == 0
        if full:
            frows.append(place_water(sy, x, (n + 1) // 10 - 1, L, zn, near, ring.get((n + 1) // 10 - 1), place))
        y = to_lab(sy, x, n, rng, B, full)
        lab_s.append((y[sy.solute], B))
        if full:
            lab_f.append((y, B))

    def write_dcd(path, top, items, interval, frames_before):
        # first step in the header: counted through the run (md/run_md.py), or one interval in every file (the
        # reporter of OpenMM), which tells nothing about the place of the segment in the run
        first = interval if run == 'X3' else (frames_before + 1) * interval
        with open(path, 'wb') as f:
            d = app.DCDFile(f, top, STEP_PS * u.picosecond, first, interval)
            for y, B in items:
                d.writeModel(y * u.angstrom,
                             periodicBoxVectors=u.Quantity([Vec3(*(B[k] / 10.0)) for k in range(3)], u.nanometer))

    a = 0
    if rep is None:
        for k, m in enumerate(seg, 1):                        # segments are numbered from 1, as md/run_md.py does
            write_dcd(rdir / f'solute_seg{k:04d}.dcd', top_s, lab_s[a:a + m], SOL_STEPS, a)
            write_dcd(rdir / f'full_seg{k:04d}.dcd', top_f, lab_f[a // 10:(a + m) // 10], FULL_STEPS, a // 10)
            a += m
    else:                                                     # second segment repeats the frames from `rep` on
        write_dcd(rdir / 'solute_seg0001.dcd', top_s, lab_s, SOL_STEPS, 0)
        write_dcd(rdir / 'solute_seg0002.dcd', top_s, lab_s[rep:], SOL_STEPS, rep)
        write_dcd(rdir / 'full_seg0001.dcd', top_f, lab_f, FULL_STEPS, 0)
    if info.get('broken_last_frame'):                         # a frame cut off while it was written
        with open(rdir / f'solute_seg{len(seg):04d}.dcd', 'ab') as f:
            f.write(b'\x30\x00\x00\x00' + b'\x00' * 137)

    # ---- planted values ------------------------------------------------------------------------------------
    pdir = root / 'planted'
    pdir.mkdir(exist_ok=True)
    tag = f'{runs_folder}__{model}__{run}'
    # version 1.2: the sign at N2 carried over frames, the stretches, the hydrogen bond around the first passage
    signs = [r['n2_sign'] for r in srows]
    start_sign = SIGN_OF_START[info['start']]
    assigned, stretches, longs, changes = g_stretches(signs, start_sign)
    for r, a in zip(srows, assigned):
        r['n2_sign_assigned'] = a
    scols = SOLUTE_COLS + (('n2_sign_assigned', 'int'),) + (ASN67_COLS if sy.with_nd2 else ())
    equil = generate_equilibration(sy, key, rdir, case12['equil'], top_s, top_f, pdir, tag, info['start'], place) \
        if (case12 and case12.get('equil')) else None
    with open(pdir / f'{tag}__solute.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['frame', 'pose', 'labels_swapped'] + [c for c, _ in scols])
        for n, r in enumerate(srows):
            w.writerow([n, frames[n][0], frames[n][1] or '-'] + [r[c] if isinstance(r[c], (int, str)) else '%.6f' % r[c]
                                                               for c, _ in scols])
    with open(pdir / f'{tag}__full.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['frame'] + [c for c, _ in FULL_COLS])
        for j, r in enumerate(frows):
            w.writerow([j] + [r[c] if isinstance(r[c], int) else '%.6f' % r[c] for c, _ in FULL_COLS])
    with open(pdir / f'{tag}__rmsf.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['residue_number_deposited', 'rmsf_A'])
        for k, i in enumerate(sy.ca):
            w.writerow([sy.ca_resnum[k], '%.6f' % (EPS * g_dist(sy.ca_ref[k], sy.ca_centroid))])
    st = [r['state'] for r in srows]
    gate = g_gate(srows, frows)                               # gate G-V by the generator's own code
    gd = GV_DESIGN.get(key)
    if gd is not None:                                        # the design, typed by hand, against the frames as built
        if gd['not_met'] is not None:
            assert gate['not_met'] == sorted(gd['not_met']), (key, gate['not_met'])
            assert gate['passed'] == (gd['verdict'] == 'PASSED'), (key, gate['passed'])
        assert [gate['fifth_solute_count'], gate['fifth_full_count']] == gd['fifth'], (key, gate['fifth_solute_count'], gate['fifth_full_count'])
        assert 5 * gate['hb_count'] == 4 * nfr, (key, gate['hb_count'])        # the hydrogen bond at its limit, 80 %
    for k_, c_, lo_, hi_ in GV_BESIDE:
        assert k_ != key or lo_ < abs(gate['deviations'][c_]) <= hi_, (key, c_, gate['deviations'][c_])
    asn = None
    if sy.with_nd2:
        asn = {'count_ND2_nearer': sum(r['asn67_ND2_nearer_N1'] for r in srows),
               'by_state': {k: {'n': st.count(k), 'count_ND2_nearer': sum(r['asn67_ND2_nearer_N1'] for r in srows if r['state'] == k)}
                            for k in 'CRI'}}
        if gd is not None:
            b = asn['by_state']
            assert [b['C']['count_ND2_nearer'], b['C']['n'], b['R']['count_ND2_nearer'], b['R']['n']] == gd['nd2'], (key, b)
    # a frame with the pose as built (labels as named, hydrogen in place, zinc in place, nothing bent or turned)
    bent = any(opt.get(k) for k in ('tau', 'csn', 'nd1'))
    as_built = None if bent else [info['start'], 8 if info['start'] == 'X' else (2 if key == 'gate/R1' else 1)]
    if case12:
        as_built = case12['as_built_frame']
    # ---- version 1.2: what the frames as built hold, counted here
    nwin = min(nfr, 2000)

    def sgn(rows):
        return [len(rows), sum(1 for r in rows if r['n2_sign'] > 0), sum(1 for r in rows if r['n2_sign'] < 0),
                sum(1 for r in rows if r['n2_sign'] == 0)]
    n2p = {'sign_of_start': start_sign, 'all': sgn(srows), 'window_frames': nwin, 'window': sgn(srows[nfr - nwin:]),
           'by_state': {k: sgn([r for r in srows if r['state'] == k]) for k in 'CRI'},
           'n_stretches': len(stretches), 'n_long_stretches': len(longs),
           'long_stretches': [{'sign': a, 'first_frame': b, 'n_frames': c} for a, b, c in longs], 'changes': changes,
           'sign_first_frame': signs[0], 'sign_last_frame': signs[-1]}
    pf = g_first_passage(st, STATE_OF_START[info['start']])
    hbs = {'first_passage_frame': pf}
    if pf is not None:
        hbs.update(n_before=pf, count_before=sum(r['hb_N2H_T199OG1'] for r in srows[:pf]), n_after=nfr - pf,
                   count_after=sum(r['hb_N2H_T199OG1'] for r in srows[pf:]))
    glu = {'count_OE1': sum(r['fifth_E106OE1'] for r in srows), 'count_OE2': sum(r['fifth_E106OE2'] for r in srows),
           'count_either': sum(r['fifth_E106OE'] for r in srows),
           'min_OE1': min(r['d_Zn_E106OE1_A'] for r in srows), 'min_OE2': min(r['d_Zn_E106OE2_A'] for r in srows),
           'full_frames_protein_O': sum(r['fifth_proteinO'] for r in frows)}
    ringp = {'counts': [r['n_water_ring_wide_3p5A'] for r in frows],
             'mean': float(np.mean([r['n_water_ring_wide_3p5A'] for r in frows])),
             'min_distance': min(r['d_water_ring_wide_min_A'] for r in frows),
             'mean_distance': float(np.mean([r['d_water_ring_wide_min_A'] for r in frows]))}
    # version 1.3: the count the driver acts on, and the settings it is made with (the manifest's, else the defaults)
    placep = {'counts': [r['n_water_ring_place'] for r in frows],
              'nearest': [r['d_water_ring_place_min_A'] for r in frows],
              'mean': float(np.mean([r['n_water_ring_place'] for r in frows])),
              'min_distance': min(r['d_water_ring_place_min_A'] for r in frows),
              'mean_distance': float(np.mean([r['d_water_ring_place_min_A'] for r in frows])),
              'settings': {'site_atoms': list(place[0]) if place else list(RING_PLACE),
                           'water_cutoff_A': place[1] if place else CUT_PLACE,
                           'source': 'manifest hydration_check.settings' if place else 'specification (defaults)',
                           'settings_are_the_defaults': (not place) or (set(place[0]) == set(RING_PLACE) and place[1] == CUT_PLACE)}}
    startp = {'state_first_frame': st[0], 'state_differs': st[0] != STATE_OF_START[info['start']],
              'state_first_frame_is_opposite': st[0] in 'CR' and st[0] != STATE_OF_START[info['start']],
              'n2_sign_first_frame': signs[0], 'n2_sign_differs': signs[0] != start_sign,
              'n2_volume_first_frame': srows[0]['n2_volume'], 'lig_rmsd_first_frame_A': srows[0]['lig_rmsd_A']}
    d12 = V12_DESIGN.get(key)
    if d12 is not None:                                       # the design typed by hand against the frames as built
        dn = d12['n2']
        assert n2p['all'][1:] == [dn['count_positive'], dn['count_negative'], dn['count_zero']], (key, n2p['all'])
        assert n2p['window'][1:3] == [dn['count_positive_window'], dn['count_negative_window']], (key, n2p['window'])
        assert all(n2p['by_state'][k][:3] == dn['by_state'][k] for k in 'CRI'), (key, n2p['by_state'])
        assert [n2p['n_stretches'], n2p['n_long_stretches']] == [dn['n_stretches'], dn['n_long_stretches']], (key, n2p['n_stretches'], n2p['n_long_stretches'])
        assert [c['frame'] for c in changes] == dn['change_frames'] and [c['to_sign'] for c in changes] == dn['change_to_sign'], (key, changes)
        dh = d12['hb']
        assert [hbs.get(k) for k in ('first_passage_frame', 'n_before', 'count_before', 'n_after', 'count_after')] == \
               [dh[k] for k in ('first_passage_frame', 'n_before', 'count_before', 'n_after', 'count_after')], (key, hbs)
        dg = d12['glu']
        assert [glu[k] for k in ('count_OE1', 'count_OE2', 'count_either', 'full_frames_protein_O')] == \
               [dg[k] for k in ('count_OE1', 'count_OE2', 'count_either', 'full_frames_protein_O')], (key, glu)
        for k in ('min_OE1', 'min_OE2'):
            assert dg[k] is None or abs(glu[k] - dg[k]) < 1e-9, (key, k, glu[k])
        if d12['ring'] is not None:
            assert ringp['counts'] == d12['ring']['counts'] and abs(ringp['mean'] - d12['ring']['mean']) < 1e-9, (key, ringp['counts'])
            assert abs(ringp['min_distance'] - d12['ring']['min_distance']) < 5e-4, (key, ringp['min_distance'])
        if d12.get('place') is not None:                      # version 1.3
            dpl = d12['place']
            assert placep['counts'] == dpl['counts'] and abs(placep['mean'] - dpl['mean']) < 1e-9, (key, placep['counts'])
            assert abs(placep['min_distance'] - dpl['min_distance']) < 5e-4, (key, placep['min_distance'])
            assert 'nearest' not in dpl or all(abs(a - b) < 5e-4 for a, b in zip(placep['nearest'], dpl['nearest'])), (key, placep['nearest'])
        if d12.get('settings') is not None:
            assert placep['settings'] == d12['settings'], (key, placep['settings'])
        assert all(startp[k] == v for k, v in d12['start'].items()), (key, startp)
        de = d12['equil']
        assert (equil is not None) == de['available'], key
        if equil is not None:
            assert [equil['n_frames_solute'], equil['n_frames_full'], equil['total_ps']] == [de['n_solute'], de['n_full'], de['total_ps']], (key, equil['n_frames_solute'], equil['n_frames_full'])
            assert [[a['n_frames_solute'], a['n_frames_full']] for a in equil['stages']] == de['frames_per_stage'], key
            for k, got in (('first_not_start', equil['first_frame_not_in_state_of_start']), ('first_opposite', equil['first_frame_in_opposite_state']),
                           ('first_sign', equil['first_frame_n2_sign_not_of_start'])):
                assert (got is None and de[k] is None) or (got is not None and [got['frame'], got['stage'], got['extension']] == de[k][:3]
                                                           and abs(got['time_ns'] - de[k][3]) < 1e-9), (key, k, got)
            if 'steps' in de:
                assert [a['tag'] for a in equil['stages']] == de['steps'], key
                assert [[a, b, round(c, 9)] for a, b, c in equil['frames_70_170']] == [de['frame_70'], de['frame_170']], (key, equil['frames_70_170'])
                assert {a['tag']: a['last_full_frame']['n_water_ring_place'] for a in equil['stages'] if a['tag'] in de['ring_at_the_end_of']} == de['ring_at_the_end_of'], key
            assert [equil['state_last_frame'], equil['n2_sign_last_frame']] == [de['state_last'], de['sign_last']], key
            if de.get('ring_counts') is not None:
                assert equil['ring_counts'] == de['ring_counts'], (key, equil['ring_counts'])
            if de.get('place_counts') is not None:            # version 1.3
                assert equil['place_counts'] == de['place_counts'], (key, equil['place_counts'])
            if de.get('place_nearest') is not None:
                assert all(abs(a - b) < 5e-4 for a, b in zip(equil['place_nearest'], de['place_nearest'])), (key, equil['place_nearest'])
    sysdir_src = ('convention (the system_dir of the manifest is not on this machine)' if (run == 'R1' and not new) else
                  'convention' if kind in ('none', 'plain-bare') else 'manifest system_dir')
    P = {'run': f'{model}/{run}', 'folder': str(rdir.relative_to(root)), 'start': info['start'],
         'n_frames_solute': nfr, 'n_frames_full': len(frows), 'seen_ns': nfr * STEP_PS * SOL_STEPS / 1000.0,
         'planned_ns': info['planned_ns'], 'segments': seg,
         'frames_ligand_moved': sum(1 for n in range(nfr) if n % 7 == 3),
         'frames_zinc_moved': sum(1 for n in range(nfr) if n % 11 == 5),
         'design': DESIGN.get(key if new else run),
         'role': 'validation' if (info['start'] == 'X' and run[1:] == '0') else 'production',
         'with_nd2': sy.with_nd2, 'opt': opt, 'as_built_frame': as_built,
         'n_warnings': gd['warnings'] if gd is not None else (d12['warnings'] if d12 is not None else 0),
         'gate': gate, 'gate_design': gd, 'asn67': asn,
         'n2': n2p, 'hb_first_passage': hbs, 'glu106': glu, 'ring_wide': ringp, 'ring_place': placep,
         'start_of_production': startp,
         'equilibration': equil, 'design_v12': d12,
         'files': {'manifest': {'driver': 'keys of md/run_md.py', 'none': 'none'}.get(kind, 'plain keys'),
                   'header_contiguous': None if (len(seg) == 1 or run == 'X3') else (rep is None),
                   'manifest_segments_available': not plain,
                   'manifest_segments_consistent': None if plain else (rep is None),
                   'system_dir_source': sysdir_src,
                   'selections_match_manifest': None if plain else (rep is None),
                   'reference_match_manifest': None if plain else True},
         'counted': {'count_C': st.count('C'), 'count_R': st.count('R'), 'count_I': st.count('I'),
                     'hb_N2H_T199OG1_count': sum(r['hb_N2H_T199OG1'] for r in srows),
                     'fifth_water_count': sum(r['fifth_water'] for r in frows),
                     'labelling_counts': [sum(1 for r in srows if r['lig_rmsd_labelling'] == i) for i in range(4)]},
         'stats': {c: stats([r[c] for r in srows]) for c in ('lig_rmsd_A', 'd_Zn_N2_A', 'tau_fold_deg', 'ang_C4_S_N2_deg',
                                                            'd_S_N2_A', 'ca_rmsd_A', 'd1_N1_F131CZ_A', 'd2_N1_N67OD1_A',
                                                            'n2_volume', 'd_Zn_E106OE1_A', 'd_Zn_E106OE2_A', 'd_Zn_E106OE_min_A')},
         'stats_full': {c: {'mean': float(np.mean([r[c] for r in frows]))} for c in ('d_Zn_Owater_min_A', 'n_water_N1_3p5A')},
         'ca_rmsd_by_construction_A': srows[0]['ca_rmsd_A']}
    if rep is not None:
        P.update(n_frames_solute=nfr + (nfr - rep), overlap_expected=True, overlap_frame_of_previous_segment=rep)
    (pdir / f'{tag}.json').write_text(json.dumps(P, indent=1) + '\n')
    return P


def generate(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    sy = System()
    # checks of the construction against md/system/starts.csv (second formulas against the system builder)
    with open(SYSD / 'starts.csv') as f:
        starts = {r['start']: r for r in csv.DictReader(f)}
    for s in 'XNR':
        v = solute_values(sy, 1, sy.pose[s], sy.zn0)
        for a, b, u in STARTS_MAP:
            assert abs(v[b] - float(starts[s][a])) <= TOL_STARTS[u], (s, a, v[b], starts[s][a])
    vt, vi = solute_values(sy, 1, sy.pose['T'], sy.zn0), solute_values(sy, 1, sy.pose['I'], sy.zn0)
    assert vt['state'] == 'R' and vi['state'] == 'I', (vt['state'], vi['state'])
    assert abs(vi['lig_rmsd_A'] - 3.0) < 1e-9 and abs(vi['d1_N1_F131CZ_A'] - (float(g_dist(sy.pose['X']['N1'], sy.S['phe131_CZ'])) + 3.0)) < 1e-9

    # version 1.1: the second system (with ND2 of Asn67), the references of gate G-V, the bent poses
    sy2 = System(with_nd2=True)
    assert sy2.n == sy.n + 1 and len(sy2.solute) == len(sy.solute) + 1
    for s in 'XNR':
        v = solute_values(sy2, 1, sy2.pose[s], sy2.zn0)
        assert abs(v['d_N1_N67ND2_A'] - float(starts[s]['N1_Asn67ND2'])) <= TOL_STARTS['A'], (s, v['d_N1_N67ND2_A'])
        for a, b, u in STARTS_MAP:
            assert abs(v[b] - float(starts[s][a])) <= TOL_STARTS[u], (s, a, v[b], starts[s][a])
    vx = solute_values(sy, 1, sy.pose['X'], sy.zn0)
    crystal_refs = {c: vx[c] for c in tuple(GV_ANGLE_REF) + ('ang_C4_S_N2_deg', 'tau_fold_deg')}
    for c, ref in {**GV_ANGLE_REF, 'ang_C4_S_N2_deg': GV_CSN[0], 'tau_fold_deg': GV_TAU[0]}.items():
        assert abs(vx[c] - ref) <= 0.05, (c, vx[c], ref)       # each reference is the crystal value to one decimal

    sel = sy.selections()
    top_ref_atoms = [(sy.name[i], sy.atoms[i][1], sy.atoms[i][2], sy.x0[i]) for i in sy.ca] + \
                    [(n, 'SAN', 304, sy.pose['X'][n]) for n in LIG_HEAVY]
    folders = [('test', s, sel) for s in 'XNR']
    folders += [(k.split('/')[0], k.split('/')[1][0], sy2.selections()) for k in GV_CASES if not k.startswith('test/')]
    folders += [(k.split('/')[0], k.split('/')[1][0], sy2.selections()) for k in V12_CASES]
    for model_, s, sel_ in folders:
        d = root / 'systems' / model_ / s
        d.mkdir(parents=True, exist_ok=True)
        (d / 'selections.json').write_text(json.dumps(sel_) + '\n')
        with open(d / 'reference.pdb', 'w') as f:
            f.write('REMARK   synthetic test system: crystal coordinates (6RL9 frame, Angstrom) of the CA atoms and of SAN 304\n')
            for k, (name, resn, resi, x) in enumerate(top_ref_atoms, 1):
                rec = 'ATOM  ' if resn != 'SAN' else 'HETATM'
                f.write(f'{rec}{k:5d} {name:<4s} {resn:>3s} A{resi:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00\n')
            f.write('END\n')
        (d / 'extra_forces.json').write_text(json.dumps({'harmonic_bonds': []}) + '\n')
        (d / 'build.json').write_text(json.dumps({'note': 'synthetic test system, no topology', 'seed': SEED}) + '\n')
    planted = [generate_run(sy, r, root) for r in ('X1', 'X2', 'X3', 'N1', 'R1')]
    planted.append(generate_run(sy, 'X9', root, runs_folder='runs_overlap'))
    for k, c in GV_CASES.items():                             # version 1.1
        model_, run_ = k.split('/')
        planted.append(generate_run(sy if model_ == 'test' else sy2, run_, root, model=model_, runs_folder=c['folder']))
    for k, c in V12_CASES.items():                            # version 1.2
        model_, run_ = k.split('/')
        planted.append(generate_run(sy2, run_, root, model=model_, runs_folder=c['folder']))
    gv = {p['run']: p for p in planted if p['run'] in GV_CASES}
    meta = {'written_by': 'md/analysis/make_test_trajectory.py', 'sha256_of_this_script': sha256(__file__), 'seed': SEED,
            'inputs_sha256': {f'md/system/{p}': sha256(SYSD / p) for p in
                              ('tleap/protein_H.pdb', 'zinc.pdb', 'ligand_X.pdb', 'ligand_N.pdb', 'ligand_R.pdb', 'starts.csv')},
            'software': software(), 'n_atoms': sy.n, 'n_solute': len(sy.solute), 'n_CA': len(sy.ca),
            'second_system_with_ND2': {'n_atoms': sy2.n, 'n_solute': len(sy2.solute)},
            'gate_GV': {'references_typed_from_the_specification': {**GV_ANGLE_REF, 'ang_C4_S_N2_deg': GV_CSN[0],
                                                                    'tau_fold_deg': GV_TAU[0]},
                        'crystal_start_by_second_formulas': crystal_refs,
                        'runs': {k: {'verdict_by_design': GV_DESIGN[k]['verdict'], 'not_met': p['gate']['not_met'],
                                     'deviations_deg': p['gate']['deviations'],
                                     'means': {c: p['gate']['means'][c] for c in ('lig_rmsd_A', 'd_Zn_N2_A', 'd_Zn_H94NE2_A',
                                                                                 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A')},
                                     'fifth_solute': [p['gate']['fifth_solute_count'], p['gate']['n_solute']],
                                     'fifth_full': [p['gate']['fifth_full_count'], p['gate']['n_full']]}
                                 for k, p in gv.items()}},
            'poses': {'T_rigid_move': {k: solute_values(sy, 1, sy.pose['T'], sy.zn0)[k] for k in
                                       ('lig_rmsd_A', 'lig_rmsd_labelling', 'd1_N1_F131CZ_A', 'd2_N1_N67OD1_A', 'state',
                                        'tau_fold_deg', 'd_S_N2_A', 'd_Zn_N2_A')},
                      'I_translation_3A': {k: vi[k] for k in ('lig_rmsd_A', 'd1_N1_F131CZ_A', 'd2_N1_N67OD1_A', 'state')}},
            'constants': {'EPS': EPS, 'DZN_A': DZN, 'DW_NEAR_A': DW_NEAR, 'DW_FAR_A': DW_FAR, 'RHO_IN_A': RHO_IN,
                          'RHO_OUT_A': RHO_OUT, 'BOX_L_A': BOX_L, 'BOX_AMP_A': BOX_AMP},
            'campaign_design': CAMPAIGN_DESIGN, 'runs': [p['run'] for p in planted]}
    (root / 'planted' / 'planted.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(f'synthetic test written to {root}: {sy.n} atoms ({len(sy.solute)} in the solute, {len(sy.ca)} CA), '
          f'{len(planted)} runs')
    return 0


def software():
    import openmm
    return {'python': sys.version.split()[0], 'numpy': np.__version__, 'openmm': openmm.__version__}


# ----------------------------------------------------------------------------------------------------------
# comparison of the analysis with the planted values
# ----------------------------------------------------------------------------------------------------------
def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def compare(root):
    root = Path(root).resolve()
    meta = json.loads((root / 'planted' / 'planted.json').read_text())
    with open(SYSD / 'starts.csv') as f:
        starts = {r['start']: r for r in csv.DictReader(f)}
    rep = []                                                  # rows of the report

    def add(run, what, planted, returned, tol, source, unit=''):
        if isinstance(planted, bool) or isinstance(planted, str) or planted is None or tol is None:
            ok, diff = planted == returned, None
        elif returned is None:
            ok, diff = False, None
        else:
            diff = abs(float(planted) - float(returned))
            ok = diff <= tol
        rep.append({'run': run, 'quantity': what, 'unit': unit, 'planted': planted, 'returned': returned,
                    'abs_diff': None if diff is None else round(diff, 6), 'tolerance': tol, 'source': source, 'passed': bool(ok)})

    for tag in sorted(p.stem for p in (root / 'planted').glob('runs*__*__*.json')):
        P = json.loads((root / 'planted' / f'{tag}.json').read_text())
        run = P['run'] + (' (overlap)' if tag.startswith('runs_overlap__') else '')
        adir = root / P['folder'] / 'analysis'
        S = json.loads((adir / 'summary.json').read_text())
        add(run, 'solute frames seen', P['n_frames_solute'], S['run']['n_frames_solute'], 0, 'construction')
        F = P['files']
        add(run, 'DCD headers: segments contiguous (none = the headers tell nothing)', F['header_contiguous'],
            S['run']['header_check_solute']['contiguous'], None, 'construction')
        add(run, 'manifest lists the segments', F['manifest_segments_available'],
            S['run']['manifest_segment_check']['available'], None, 'construction')
        add(run, 'manifest and files agree (none = no list)', F['manifest_segments_consistent'],
            S['run']['manifest_segment_check']['consistent'], None, 'construction')
        add(run, 'interval in the DCD headers as used', True, S['run']['header_check_solute']['interval_matches'], None, 'construction')
        add(run, 'where the system folder was found', F['system_dir_source'], S['run']['sources']['system_dir'], None, 'construction')
        add(run, 'selections.json is the file of the manifest (none = no checksum there)', F['selections_match_manifest'],
            S['run']['system_files_match_manifest']['selections.json'], None, 'construction')
        add(run, 'reference.pdb is the file of the manifest (none = no checksum there)', F['reference_match_manifest'],
            S['run']['system_files_match_manifest']['reference.pdb'], None, 'construction')
        if P.get('overlap_expected'):
            oc = S['run']['segment_overlap_check']
            add(run, 'repeated frames found at the segment boundary', True, oc['suspected'], None, 'construction')
            add(run, 'frame of the first segment that is repeated', P['overlap_frame_of_previous_segment'],
                oc['boundaries'][0]['frame_of_previous_segment'], 0, 'construction')
            add(run, 'warnings written (headers, repeated frame, manifest, system file)', 4, len(S['warnings']), 0, 'construction')
            continue
        add(run, 'full frames seen', P['n_frames_full'], S['run']['n_frames_full'], 0, 'construction')
        add(run, 'time seen', P['seen_ns'], S['run']['seen_ns'], 1e-9, 'construction', 'ns')
        add(run, 'repeated frames found at a segment boundary', False, S['run']['segment_overlap_check']['suspected'], None, 'construction')
        add(run, 'frames in which the ligand was in another periodic image', P['frames_ligand_moved'],
            S['run']['imaging_solute']['frames_ligand_moved'], 0, 'construction')
        add(run, 'frames in which the zinc was in another periodic image', P['frames_zinc_moved'],
            S['run']['imaging_solute']['frames_zinc_moved'], 0, 'construction')
        add(run, 'warnings', P['n_warnings'], len(S['warnings']), 0, 'construction')

        # every frame, every column
        mine, theirs = read_csv(root / 'planted' / f'{tag}__solute.csv'), read_csv(adir / 'timeseries_solute.csv')
        add(run, 'lines of timeseries_solute.csv', len(mine), len(theirs), 0, 'construction')
        add(run, 'columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2)', P['with_nd2'],
            'd_N1_N67ND2_A' in theirs[0] and 'asn67_ND2_nearer_N1' in theirs[0], None, 'construction')
        for c, kind in SOLUTE_COLS + (('n2_sign_assigned', 'int'),) + (ASN67_COLS if P['with_nd2'] else ()):
            if kind in ('A', 'deg', 'tors', 'vol'):
                d = np.array([abs(float(a[c]) - float(b[c])) for a, b in zip(mine, theirs)])
                if kind == 'tors':
                    d = np.minimum(d, 360.0 - d)
                u = {'A': 'A', 'vol': 'vol'}.get(kind, 'deg')
                add(run, f'every frame: {c}, largest difference', 0.0, float(d.max()), TOL[u], 'second formulas', '' if u == 'vol' else u)
            else:
                bad = sum(1 for a, b in zip(mine, theirs) if str(a[c]) != str(b[c]))
                add(run, f'every frame: {c}, frames that differ', 0, bad, 0, 'second formulas')
        mine, theirs = read_csv(root / 'planted' / f'{tag}__full.csv'), read_csv(adir / 'timeseries_full.csv')
        add(run, 'lines of timeseries_full.csv', len(mine), len(theirs), 0, 'construction')
        for c, kind in FULL_COLS:
            if kind == 'A':
                d = np.array([abs(float(a[c]) - float(b[c])) for a, b in zip(mine, theirs)])
                add(run, f'every full frame: {c}, largest difference', 0.0, float(d.max()), TOL['A'],
                    'construction' if c == 'd_Zn_Owater_min_A' else 'second formulas', 'A')
            else:
                bad = sum(1 for a, b in zip(mine, theirs) if str(a[c]) != str(b[c]))
                add(run, f'every full frame: {c}, frames that differ', 0, bad, 0,
                    'construction' if c in ('n_water_N1_3p5A', 'fifth_water', 'water_O_index_nearest_Zn') else 'second formulas')
        mine, theirs = read_csv(root / 'planted' / f'{tag}__rmsf.csv'), read_csv(adir / 'rmsf_ca.csv')
        d = np.array([abs(float(a['rmsf_A']) - float(b['rmsf_A'])) for a, b in zip(mine, theirs)])
        add(run, 'CA RMSF of 258 residues, largest difference', 0.0, float(d.max()), TOL['A'], 'construction', 'A')
        add(run, 'residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261)', 0,
            sum(1 for a, b in zip(mine, theirs) if a['residue_number_deposited'] != b['residue_number_deposited'])
            + abs(len(mine) - len(theirs)), 0, 'construction')

        # the numbers fixed by the design of the run
        D, O, V = P['design'], S['orientation'], S['validation']
        add(run, 'verdict', D['verdict'], O['verdict'], None, 'construction')
        add(run, 'verdict provisional (run not finished)', D['provisional'], O['provisional'], None, 'construction')
        add(run, 'run complete', D['complete'], S['run']['complete'], None, 'construction')
        add(run, 'first passage into the opposite state', D['first_passage_ns'], O['first_passage_ns'],
            None if D['first_passage_ns'] is None else 1e-9, 'construction', 'ns')
        add(run, 'transitions', D['n_transitions'], O['n_transitions'], 0, 'construction')
        add(run, 'long visits (at least 1 ns)', D['n_long_visits'], O['n_long_visits'], 0, 'construction')
        add(run, 'frames of the verdict window', D['window_frames'], O['window_frames'], 0, 'construction')
        for k in ('count_C_window', 'count_R_window', 'count_I_window', 'count_C', 'count_R', 'count_I'):
            add(run, k, D[k], O[k], 0, 'construction')
        add(run, 'rule 1: left the crystal pose', D['rule1_left_crystal'], V['rule1_left_crystal'], None, 'construction')
        for k in ('count_C', 'count_R', 'count_I'):
            add(run, f'{k} (design against the planted frames)', D[k], P['counted'][k], 0, 'construction')
        add(run, 'frames with the hydrogen bond N2-H...OG1', P['counted']['hb_N2H_T199OG1_count'],
            S['fractions']['hb_N2H_T199OG1']['count'], 0, 'second formulas')
        add(run, 'full frames with a water within 2.6 A of the zinc', P['counted']['fifth_water_count'],
            S['fractions']['fifth_water']['count'], 0, 'construction')
        for i, name in enumerate(('identity', 'ring_flip', 'oxygen_swap', 'ring_flip+oxygen_swap')):
            add(run, f'frames with labelling {name}', P['counted']['labelling_counts'][i],
                S['lig_rmsd_labelling_counts'][name], 0, 'second formulas')
        for c, st in P['stats'].items():
            u = 'A' if c.endswith('_A') else 'vol' if c == 'n2_volume' else 'deg'
            for k in ('mean', 'sd', 'mean_last_half', 'sd_last_half'):
                add(run, f'{c}: {k}', round(st[k], 6), S['stats'][c][k], TOL[u], 'second formulas', '' if u == 'vol' else u)
            dmax = max([abs(a - b) for a, b in zip(st['block_means'], S['stats'][c]['block_means'])] or [0.0])
            add(run, f'{c}: number of 10 ns blocks', len(st['block_means']), len(S['stats'][c]['block_means']), 0, 'construction')
            add(run, f'{c}: block means, largest difference', 0.0, dmax, TOL[u], 'second formulas', '' if u == 'vol' else u)
        add(run, 'CA RMSD, mean (0.01 x radius of gyration of the CA atoms)', round(P['ca_rmsd_by_construction_A'], 6),
            S['stats']['ca_rmsd_A']['mean'], TOL['A'], 'construction', 'A')
        add(run, 'nearest water to the zinc, mean', round(P['stats_full']['d_Zn_Owater_min_A']['mean'], 6),
            S['stats']['d_Zn_Owater_min_A']['mean'], TOL['A'], 'construction', 'A')
        add(run, 'water molecules within 3.5 A of the aniline N, mean', round(P['stats_full']['n_water_N1_3p5A']['mean'], 6),
            S['stats']['n_water_N1_3p5A']['mean'], 1e-6, 'construction')

        # one frame with the pose as built (no zinc move, hydrogen in place, labels as named) against starts.csv;
        # not for the runs with a bent pose
        if P['as_built_frame'] is not None:
            srow, k1 = P['as_built_frame']
            f1 = read_csv(adir / 'timeseries_solute.csv')[k1]
            for a, b, u in STARTS_MAP + ((('N1_Asn67ND2', 'd_N1_N67ND2_A', 'A'),) if P['with_nd2'] else ()):
                add(run, f'frame {f1["frame"]}, {b} against starts.csv ({srow}, {a})', float(starts[srow][a]), float(f1[b]),
                    TOL_STARTS[u], 'starts.csv', u)

        # version 1.1: the role of the run, gate G-V (the generator's own evaluation, g_gate, and the design typed
        # by hand), the optional Asn67 ND2
        G, PG, GD = S['gate_GV'], P['gate'], P['gate_design']
        add(run, 'role of the run (X0 = validation)', P['role'], S['run']['role'], None, 'construction')
        add(run, 'planned length (X0: 5 ns by the specification)', P['planned_ns'], S['run']['planned_ns'], 1e-9, 'construction', 'ns')
        add(run, 'gate G-V applies to this run (only to X0)', P['role'] == 'validation', G['applies_to_this_run'], None, 'construction')
        add(run, 'rule 1 applies to this run (production runs of start X)', P['role'] == 'production' and P['start'] == 'X',
            V['rule1_applies_to_this_run'], None, 'construction')
        zc, kc = G['GV3_zinc']['criteria'], G['GV4_ligand_internal']['criteria']
        got = {'GV1 ligand RMSD': G['GV1_lig_rmsd']['passed'], 'GV2 hydrogen bond': G['GV2_hbond']['passed'],
               'GV3 Zn-N2': zc['d_Zn_N2_A']['passed'],
               **{f'GV3 {c}': zc[c]['passed'] for c in ('d_Zn_H94NE2_A', 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A') + tuple(GV_ANGLE_REF)},
               'GV3 fifth ligand, solute trajectory': zc['fifth_ligand']['solute_trajectory']['passed'],
               'GV3 fifth ligand, full trajectory': zc['fifth_ligand']['full_trajectory']['passed'],
               'GV4 ang_C4_S_N2_deg': kc['ang_C4_S_N2_deg']['passed'], 'GV4 tau_fold_deg': kc['tau_fold_deg']['passed']}
        add(run, 'gate G-V: the criteria that are not met', ', '.join(PG['not_met']) or 'none',
            ', '.join(sorted(k for k, v in got.items() if not v)) or 'none', None, 'second formulas')
        add(run, 'gate G-V: all criteria met', PG['passed'], G['passed'], None, 'second formulas')
        add(run, 'gate G-V: GV3 as a whole', all(v for k, v in PG['met'].items() if k.startswith('GV3')), G['GV3_zinc']['passed'], None, 'second formulas')
        add(run, 'gate G-V: GV4 as a whole', all(v for k, v in PG['met'].items() if k.startswith('GV4')), G['GV4_ligand_internal']['passed'], None, 'second formulas')
        add(run, 'gate G-V: fifth-ligand criterion as a whole (both fractions)',
            PG['met']['GV3 fifth ligand, solute trajectory'] and PG['met']['GV3 fifth ligand, full trajectory'],
            zc['fifth_ligand']['passed'], None, 'second formulas')
        add(run, 'gate G-V: verdict (none for a run that is not X0)',
            None if P['role'] != 'validation' else ('PASSED' if PG['passed'] else 'FAILED'), G['verdict'], None, 'second formulas')
        add(run, 'gate G-V: provisional (run not finished)', not D['complete'], G['provisional'], None, 'construction')
        for c in GV_ANGLE_REF:
            add(run, f'gate G-V: {c}, mean minus its value in 6RL9', round(PG['deviations'][c], 6), zc[c]['deviation_deg'], TOL['deg'], 'second formulas', 'deg')
        for c in ('ang_C4_S_N2_deg', 'tau_fold_deg'):
            add(run, f'gate G-V: {c}, mean minus the reference', round(PG['deviations'][c], 6), kc[c]['deviation_deg'], TOL['deg'], 'second formulas', 'deg')
        fl = zc['fifth_ligand']
        add(run, 'gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc', PG['fifth_solute_count'],
            fl['solute_trajectory']['count'], 0, 'second formulas')
        add(run, 'gate G-V: solute frames counted for that fraction', PG['n_solute'], fl['solute_trajectory']['n'], 0, 'construction')
        add(run, 'gate G-V: full frames with any oxygen within 2.6 A of the zinc', PG['fifth_full_count'],
            fl['full_trajectory']['count'], 0, 'second formulas')
        add(run, 'gate G-V: full frames counted for that fraction', PG['n_full'], fl['full_trajectory']['n'], 0, 'construction')
        add(run, 'campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9',
            round(max(abs(PG['deviations'][c]) for c in GV_ANGLE_REF), 6), S['table_row']['ang_N_Zn_N_largest_abs_deviation_deg'],
            TOL['deg'], 'second formulas', 'deg')
        if GD is not None:                                    # typed by hand in GV_DESIGN
            if GD['not_met'] is not None:
                add(run, 'gate G-V by design: verdict', GD['verdict'], G['verdict'], None, 'construction')
                add(run, 'gate G-V by design: the criteria that are not met', ', '.join(sorted(GD['not_met'])) or 'none',
                    ', '.join(sorted(k for k, v in got.items() if not v)) or 'none', None, 'construction')
            add(run, 'gate G-V by design: provisional', GD['provisional'], G['provisional'], None, 'construction')
            add(run, 'by design: solute frames with a sulfonyl oxygen within 2.6 A', GD['fifth'][0], S['fractions']['fifth_sulfonylO']['count'], 0, 'construction')
            add(run, 'by design: full frames with any oxygen within 2.6 A', GD['fifth'][1], S['fractions']['fifth_any']['count'], 0, 'construction')
        AS = S['asn67_ND2']
        add(run, 'Asn67 ND2 reported (the key asn67_ND2 is in selections.json)', P['with_nd2'], AS['available'], None, 'construction')
        if P['with_nd2']:
            add(run, 'frames with ND2 nearer to the aniline N than OD1', P['asn67']['count_ND2_nearer'], AS['count_ND2_nearer'], 0, 'second formulas')
            for k in 'CRI':
                add(run, f'state {k}: frames, and frames with ND2 nearer',
                    [P['asn67']['by_state'][k]['n'], P['asn67']['by_state'][k]['count_ND2_nearer']],
                    [AS['by_state'][k]['n'], AS['by_state'][k]['count_ND2_nearer']], None, 'second formulas')
            if GD is not None:
                b = AS['by_state']
                add(run, 'by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R', GD['nd2'],
                    [b['C']['count_ND2_nearer'], b['C']['n'], b['R']['count_ND2_nearer'], b['R']['n']], None, 'construction')
            add(run, 'campaign-table value: fraction of frames with ND2 nearer', round(P['asn67']['count_ND2_nearer'] / P['n_frames_solute'], 6),
                S['table_row']['asn67_ND2_nearer_N1_fraction'], 1e-6, 'second formulas')
        else:
            add(run, 'campaign-table value: fraction of frames with ND2 nearer (none without the key)', None,
                S['table_row']['asn67_ND2_nearer_N1_fraction'], None, 'construction')

        compare_v12(add, run, P, S, adir, root, tag)

    # values planted by hand for single runs (not derived from the planted frames)
    def S_of(r):
        return json.loads((root / 'runs' / 'test' / r / 'analysis' / 'summary.json').read_text())
    s2 = S_of('X2')
    add('test/X2', 'ligand RMSD, mean (labels swapped in 3 of 4 frames)', 0.0, s2['stats']['lig_rmsd_A']['mean'], TOL['A'], 'construction', 'A')
    add('test/X2', 'ligand RMSD, largest value', 0.0, s2['stats']['lig_rmsd_A']['max'], TOL['A'], 'construction', 'A')
    add('test/X2', 'ligand RMSD without the symmetry correction, largest value (must be large)', True,
        s2['stats']['lig_rmsd_identity_A']['max'] > 1.0, None, 'construction')
    add('test/X2', 'hydrogen bond N2-H...OG1, fraction of frames', 0.8, s2['fractions']['hb_N2H_T199OG1']['fraction'], 1e-9, 'construction')
    add('test/X2', 'hydrogen bond N2-H...OG1, fraction in the last half', 0.8, s2['fractions']['hb_N2H_T199OG1']['fraction_last_half'], 1e-9, 'construction')
    add('test/X2', 'hydrogen bond N2-H...OG1, the three block fractions', [0.8, 0.8, 0.8], s2['fractions']['hb_N2H_T199OG1']['block_fractions'], None, 'construction')
    add('test/X2', 'Zn-N2, mean (starts.csv 1.931 + 0.300 in every tenth frame)', round(float(starts['X']['Zn_N2']) + 0.03, 6),
        s2['stats']['d_Zn_N2_A']['mean'], 0.0015, 'starts.csv', 'A')
    add('test/X2', 'Zn-N2, largest value (starts.csv 1.931 + 0.300)', round(float(starts['X']['Zn_N2']) + 0.3, 6),
        s2['stats']['d_Zn_N2_A']['max'], 0.0015, 'starts.csv', 'A')
    add('test/X2', 'nearest water to the zinc, mean (0.1 x 2.400 + 0.9 x 3.690)', 3.561, s2['stats']['d_Zn_Owater_min_A']['mean'], TOL['A'], 'construction', 'A')
    add('test/X2', 'nearest water to the zinc, smallest', 2.4, s2['stats']['d_Zn_Owater_min_A']['min'], TOL['A'], 'construction', 'A')
    add('test/X2', 'fraction of full frames with a fifth ligand (water)', 0.1, s2['fractions']['fifth_water']['fraction'], 1e-9, 'construction')
    add('test/X2', 'water molecules within 3.5 A of the aniline N, mean', 1.5, s2['stats']['n_water_N1_3p5A']['mean'], 1e-9, 'construction')
    add('test/X2', 'tau_fold, mean (starts.csv, unchanged by the label swaps)', float(starts['X']['tau_fold']), s2['stats']['tau_fold_deg']['mean'], 0.06, 'starts.csv', 'deg')
    add('test/X2', 'tau_fold, standard deviation', 0.0, s2['stats']['tau_fold_deg']['sd'], TOL['deg'], 'construction', 'deg')
    g = s2['gate_GV']
    add('test/X2', 'gate G-V, criterion 1 (mean ligand RMSD at most 1.0 A)', True, g['GV1_lig_rmsd']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 2 (hydrogen bond in at least 80 % of frames; here exactly 80 %)', True, g['GV2_hbond']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 3, no fifth ligand (planted: water at 2.400 A in 10 % of the full frames)', False, g['GV3_zinc']['criteria']['fifth_ligand']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 3, fifth ligand in the solute trajectory alone (no sulfonyl oxygen near)', True,
        g['GV3_zinc']['criteria']['fifth_ligand']['solute_trajectory']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 3, full frames with an oxygen within 2.6 A, of 300', 30,
        g['GV3_zinc']['criteria']['fifth_ligand']['full_trajectory']['count'], 0, 'construction')
    add('test/X2', 'gate G-V, criterion 3, Zn-N2 window', True, g['GV3_zinc']['criteria']['d_Zn_N2_A']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 3, angle His96 NE2-Zn-His119 ND1 within 10 deg of 98.4 (version 1.0: in 100-120 deg, not met)', True,
        g['GV3_zinc']['criteria']['ang_H96NE2_Zn_H119ND1_deg']['passed'], None, 'construction')
    add('test/X2', 'gate G-V, criterion 4 (C4-S-N2 within 3 deg of 108.8, tau_fold within 20 deg of 132.1)', True,
        g['GV4_ligand_internal']['passed'], None, 'construction')
    add('test/X2', 'gate G-V passed', False, g['passed'], None, 'construction')
    add('test/X2', 'gate G-V verdict (a production run has none)', None, g['verdict'], None, 'construction')
    s1 = S_of('X1')
    add('test/X1', 'fraction of the last 20 ns in the opposite state (boundary of the 90 % rule)', 0.9,
        s1['orientation']['fraction_opposite_state_window'], 1e-9, 'construction')
    add('test/X1', 'tau_fold, mean (rigid moves and label swaps keep the crystal value of starts.csv)', float(starts['X']['tau_fold']),
        s1['stats']['tau_fold_deg']['mean'], 0.06, 'starts.csv', 'deg')
    add('test/X1', 'S-N2, standard deviation (rigid moves)', 0.0, s1['stats']['d_S_N2_A']['sd'], TOL['A'], 'construction', 'A')
    ts1 = read_csv(root / 'runs' / 'test' / 'X1' / 'analysis' / 'timeseries_solute.csv')
    add('test/X1', 'frame 960 (ligand translated by 3.000 A): ligand RMSD', 3.0, float(ts1[960]['lig_rmsd_A']), TOL['A'], 'construction', 'A')
    add('test/X1', 'frame 960: d1 (starts.csv 3.897 + 3.000)', round(float(starts['X']['N1_Phe131CZ']) + 3.0, 6),
        float(ts1[960]['d1_N1_F131CZ_A']), 0.0015, 'starts.csv', 'A')
    add('test/X1', 'frame 960: state', 'I', ts1[960]['state'], None, 'construction')
    add('test/X1', 'frame 960: assigned state (carried from the visit it lies in)', 'R', ts1[960]['state_assigned'], None, 'construction')
    add('test/X1', 'frames 1602 and 1603 (labels swapped in one): same ligand RMSD', 0.0,
        abs(float(ts1[1602]['lig_rmsd_A']) - float(ts1[1603]['lig_rmsd_A'])), TOL['A'], 'construction', 'A')
    add('test/X1', 'frames 1602 and 1603: labelling differs', True, ts1[1602]['lig_rmsd_labelling'] != ts1[1603]['lig_rmsd_labelling'], None, 'construction')
    add('test/X1', 'frame 320 (short excursion, 0.5 ns): assigned state', 'R', ts1[320]['state_assigned'], None, 'construction')
    sn = S_of('N1')
    add('test/N1', 'verdict window complete (20 ns)', False, sn['orientation']['window_complete'], None, 'construction')
    add('test/N1', 'ligand RMSD, mean (starts.csv, near-native pose)', float(starts['N']['rmsd_to_crystal_A']), sn['stats']['lig_rmsd_A']['mean'], 0.0015, 'starts.csv', 'A')

    # campaign
    cj = root / 'runs' / 'analysis_campaign' / 'campaign.json'
    if cj.is_file():
        C = json.loads(cj.read_text())['rules_per_model']['test']
        for k, v in meta['campaign_design'].items():
            got = C
            for part in k.split('.'):
                got = got[part]
            add('campaign', k, v, got, None if isinstance(v, (bool, str)) else 0, 'construction')
        Call = json.loads(cj.read_text())
        add('campaign', 'folders left out of the campaign (the validation run X0)', 'test/X0',
            ' '.join('/'.join(Path(x).parts[-2:]) for x in Call['skipped']), None, 'construction')
        add('campaign', 'runs that could not be analysed', 0, len(Call['failed']), 0, 'construction')
        rows = read_csv(root / 'runs' / 'analysis_campaign' / 'campaign_table.csv')
        add('campaign', 'lines of campaign_table.csv', 5, len(rows), 0, 'construction')
        add('campaign', 'order of the runs', 'test/X1 test/X2 test/X3 test/N1 test/R1', ' '.join(r['run'] for r in rows), None, 'construction')
    else:
        add('campaign', 'campaign.json exists', True, False, None, 'construction')

    n_pass = sum(r['passed'] for r in rep)
    out = {'written_by': 'md/analysis/make_test_trajectory.py --compare', 'sha256_of_this_script': sha256(__file__),
           'sha256_of_analyse_md.py': sha256(HERE / 'analyse_md.py'), 'planted': meta,
           'tolerances': {'distance_A': TOL['A'], 'angle_deg': TOL['deg'], 'against_starts_csv_A': TOL_STARTS['A'],
                          'against_starts_csv_deg': TOL_STARTS['deg'], 'counts': 'identical'},
           'n_values': len(rep), 'n_passed': n_pass, 'n_failed': len(rep) - n_pass, 'all_passed': n_pass == len(rep),
           'values': rep}
    (root / 'test_report.json').write_text(json.dumps(out, indent=1) + '\n')

    def fmt(v):
        if v is None:
            return 'none'
        if isinstance(v, bool):
            return 'yes' if v else 'no'
        if isinstance(v, float):
            return f'{v:.6g}'
        return str(v)
    md = ['# Synthetic test of analyse_md.py: planted against returned', '',
          f'{n_pass} of {len(rep)} values within tolerance ({len(rep) - n_pass} failed). Tolerances: 0.001 Å, 0.01°, '
          'identical counts; against `md/system/starts.csv` (three decimals) 0.0015 Å and 0.06°.', '',
          '| Run | Quantity | Planted | Returned | Difference | Tolerance | Source of the planted value | Result |',
          '|---|---|---|---|---|---|---|---|']
    for r in rep:
        md.append(f'| {r["run"]} | {r["quantity"]}{" (" + r["unit"] + ")" if r["unit"] else ""} | {fmt(r["planted"])} | '
                  f'{fmt(r["returned"])} | {fmt(r["abs_diff"]) if r["abs_diff"] is not None else "–"} | '
                  f'{fmt(r["tolerance"]) if r["tolerance"] is not None else "equal"} | {r["source"]} | '
                  f'{"passed" if r["passed"] else "FAILED"} |')
    (root / 'test_report.md').write_text('\n'.join(md) + '\n')
    for r in rep:
        if not r['passed']:
            print(f'FAILED  {r["run"]}: {r["quantity"]}: planted {fmt(r["planted"])}, returned {fmt(r["returned"])}')
    print(f'{n_pass} of {len(rep)} planted values returned within tolerance; report: {root / "test_report.md"}')
    return 0 if n_pass == len(rep) else 1


def compare_v12(add, run, P, S, adir, root, tag):
    """Version 1.2: the reported quantities of sections 9A to 9E of one run against what the generator counted in
    its own frames (second formulas) and against the design typed by hand (V12_DESIGN, construction)."""
    def l6(v):
        return None if v is None else round(v, 6)
    N2, PN = S['n2_configuration'], P['n2']
    nall = PN['all'][0]
    add(run, 'N2: sign of the start as built (X and N +1, R -1)', PN['sign_of_start'], N2['sign_of_start'], 0, 'construction')
    add(run, 'N2: frames with positive, negative and zero signed volume', PN['all'][1:],
        [N2['count_positive'], N2['count_negative'], N2['count_zero']], None, 'second formulas')
    add(run, 'N2: fraction of frames with positive signed volume', l6(PN['all'][1] / nall), N2['fraction_positive'], 1e-6, 'second formulas')
    add(run, 'N2: fraction of frames with negative signed volume', l6(PN['all'][2] / nall), N2['fraction_negative'], 1e-6, 'second formulas')
    add(run, 'N2: frames of the window (last 20 ns)', PN['window_frames'], N2['window_frames'], 0, 'construction')
    add(run, 'N2: positive, negative and zero in the last 20 ns', PN['window'][1:],
        [N2['count_positive_window'], N2['count_negative_window'], N2['count_zero_window']], None, 'second formulas')
    add(run, 'N2: fraction positive in the last 20 ns', l6(PN['window'][1] / PN['window'][0]), N2['fraction_positive_window'], 1e-6, 'second formulas')
    for k in 'CRI':
        b = N2['by_state'][k]
        add(run, f'N2, frames of state {k}: number, positive, negative, zero', PN['by_state'][k],
            [b['n'], b['count_positive'], b['count_negative'], b['count_zero']], None, 'second formulas')
        pk = PN['by_state'][k]
        add(run, f'N2, frames of state {k}: fraction positive (none without a frame of that state)',
            l6(pk[1] / pk[0]) if pk[0] else None, b['fraction_positive'], 1e-6 if pk[0] else None, 'second formulas')
    add(run, 'N2: stretches of one sign, and those of at least 1 ns', [PN['n_stretches'], PN['n_long_stretches']],
        [N2['n_stretches'], N2['n_long_stretches']], None, 'second formulas')
    add(run, 'N2: the stretches of at least 1 ns (sign, first frame, frames)',
        [[x['sign'], x['first_frame'], x['n_frames']] for x in PN['long_stretches']],
        [[x['sign'], x['first_frame'], x['n_frames']] for x in N2['long_stretches']], None, 'second formulas')
    add(run, 'N2: number of sign changes that last at least 1 ns', len(PN['changes']), N2['n_lasting_sign_changes'], 0, 'second formulas')
    add(run, 'N2: these sign changes (frame, from, to)', [[c['frame'], c['from_sign'], c['to_sign']] for c in PN['changes']],
        [[c['frame'], c['from_sign'], c['to_sign']] for c in N2['lasting_sign_changes']], None, 'second formulas')
    add(run, 'N2: times of these sign changes', [c['time_ns'] for c in PN['changes']],
        [c['time_ns'] for c in N2['lasting_sign_changes']], None, 'second formulas', 'ns')
    add(run, 'N2: time of the first lasting sign change', PN['changes'][0]['time_ns'] if PN['changes'] else None,
        N2['first_lasting_sign_change_ns'], 1e-9 if PN['changes'] else None, 'second formulas', 'ns')
    add(run, 'N2: sign of the first and of the last frame', [PN['sign_first_frame'], PN['sign_last_frame']],
        [N2['sign_first_frame'], N2['sign_last_frame']], None, 'second formulas')

    H, PH = S['hbond_first_passage'], P['hb_first_passage']
    pf = PH['first_passage_frame']
    add(run, 'hydrogen bond around the first passage: reported (a first passage exists)', pf is not None, H['available'], None, 'second formulas')
    add(run, 'hydrogen bond around the first passage: frame of the first passage', pf, H['first_passage_frame'], None, 'second formulas')
    if pf is not None:
        add(run, 'hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond', [PH['n_before'], PH['count_before']],
            [H['n_before'], H['count_before']], None, 'second formulas')
        add(run, 'hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond', [PH['n_after'], PH['count_after']],
            [H['n_after'], H['count_after']], None, 'second formulas')
        add(run, 'hydrogen bond N2-H...OG1, fraction before the first passage',
            l6(PH['count_before'] / PH['n_before']) if PH['n_before'] else None, H['fraction_before'], 1e-6 if PH['n_before'] else None, 'second formulas')
        add(run, 'hydrogen bond N2-H...OG1, fraction from the first passage on', l6(PH['count_after'] / PH['n_after']), H['fraction_after'], 1e-6, 'second formulas')
    else:
        add(run, 'hydrogen bond around the first passage: the six numbers (none without a first passage)', [None] * 6,
            [H[k] for k in ('n_before', 'count_before', 'fraction_before', 'n_after', 'count_after', 'fraction_after')], None, 'construction')
    add(run, 'campaign-table values: hydrogen-bond fraction before and after the first passage',
        [H['fraction_before'], H['fraction_after']],
        [S['table_row']['hb_N2H_T199OG1_fraction_before_first_passage'], S['table_row']['hb_N2H_T199OG1_fraction_after_first_passage']], None, 'construction')

    G6, FR, ST = P['glu106'], S['fractions'], S['stats']
    add(run, 'Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc', [G6['count_OE1'], G6['count_OE2'], G6['count_either']],
        [FR['fifth_E106OE1']['count'], FR['fifth_E106OE2']['count'], FR['fifth_E106OE']['count']], None, 'second formulas')
    add(run, 'Glu106: smallest distance zinc to OE1', l6(G6['min_OE1']), ST['d_Zn_E106OE1_A']['min'], TOL['A'], 'second formulas', 'A')
    add(run, 'Glu106: smallest distance zinc to OE2', l6(G6['min_OE2']), ST['d_Zn_E106OE2_A']['min'], TOL['A'], 'second formulas', 'A')
    add(run, 'Glu106: full frames with a protein oxygen within 2.6 A of the zinc', G6['full_frames_protein_O'], FR['fifth_proteinO']['count'], 0, 'second formulas')
    rg = S['gate_GV']['reported_zn_glu106']
    for o, c, f in (('OE1', 'd_Zn_E106OE1_A', 'count_OE1'), ('OE2', 'd_Zn_E106OE2_A', 'count_OE2')):
        add(run, f'Glu106 beside gate G-V, {o}: mean distance', l6(P['stats'][c]['mean']), rg[c]['mean_A'], TOL['A'], 'second formulas', 'A')
        add(run, f'Glu106 beside gate G-V, {o}: smallest distance', l6(G6['min_' + o]), rg[c]['min_A'], TOL['A'], 'second formulas', 'A')
        add(run, f'Glu106 beside gate G-V, {o}: frames within 2.6 A, frames counted', [G6[f], nall], [rg[c]['count_within'], rg[c]['n']], None, 'second formulas')
    add(run, 'campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A',
        [ST['d_Zn_E106OE1_A']['mean'], ST['d_Zn_E106OE1_A']['min'], l6(G6['count_OE1'] / nall)],
        [S['table_row']['d_Zn_E106OE1_mean_A'], S['table_row']['d_Zn_E106OE1_min_A'], S['table_row']['fifth_E106OE1_fraction']], None, 'second formulas')

    PR = P['ring_wide']
    add(run, 'water at the crystal ring positions: mean number within 3.5 A', l6(PR['mean']), ST['n_water_ring_wide_3p5A']['mean'], 1e-6, 'second formulas')
    add(run, 'water at the crystal ring positions: smallest distance of a water oxygen', l6(PR['min_distance']), ST['d_water_ring_wide_min_A']['min'], TOL['A'], 'second formulas', 'A')
    add(run, 'water at the crystal ring positions: mean distance of the nearest water oxygen', l6(PR['mean_distance']), ST['d_water_ring_wide_min_A']['mean'], TOL['A'], 'second formulas', 'A')
    add(run, 'campaign-table value: mean number of water oxygens at the crystal ring positions', l6(PR['mean']), S['table_row']['n_water_ring_wide_3p5A_mean'], 1e-6, 'second formulas')
    # version 1.3: the count the driver acts on
    PP, WP = P['ring_place'], S['water_ring_place']
    add(run, 'water in the place of the ring: site atoms and cut-off used', [PP['settings']['site_atoms'], PP['settings']['water_cutoff_A']],
        [WP['site_atoms'], WP['water_cutoff_A']], None, 'construction', 'A')
    add(run, 'water in the place of the ring: where the settings came from, and whether they are the defaults',
        [PP['settings']['source'], PP['settings']['settings_are_the_defaults']], [WP['source'], WP['settings_are_the_defaults']], None, 'construction')
    add(run, 'water in the place of the ring: mean number', l6(PP['mean']), ST['n_water_ring_place']['mean'], 1e-6, 'second formulas')
    add(run, 'water in the place of the ring: largest number in a full frame', max(PP['counts']), ST['n_water_ring_place']['max'], 1e-6, 'second formulas')
    add(run, 'water in the place of the ring: smallest distance of a water oxygen', l6(PP['min_distance']), ST['d_water_ring_place_min_A']['min'], TOL['A'], 'second formulas', 'A')
    add(run, 'water in the place of the ring: mean distance of the nearest water oxygen', l6(PP['mean_distance']), ST['d_water_ring_place_min_A']['mean'], TOL['A'], 'second formulas', 'A')
    add(run, 'campaign-table values: water in the place of the ring, mean number and smallest distance',
        [ST['n_water_ring_place']['mean'], ST['d_water_ring_place_min_A']['min']],
        [S['table_row']['n_water_ring_place_mean'], S['table_row']['d_water_ring_place_min_min_A']], None, 'second formulas')
    add(run, 'campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count',
        ST['d_water_ring_wide_min_A']['min'], S['table_row']['d_water_ring_wide_min_min_A'], None, 'second formulas')

    PS, SP = P['start_of_production'], S['start_of_production']
    add(run, 'first production frame: state, and whether it differs from the state of the start', [PS['state_first_frame'], PS['state_differs']],
        [SP['state_first_frame'], SP['state_differs']], None, 'second formulas')
    add(run, 'first production frame: in the opposite state', PS['state_first_frame_is_opposite'], SP['state_first_frame_is_opposite'], None, 'second formulas')
    add(run, 'first production frame: sign at N2, and whether it differs from the sign of the start', [PS['n2_sign_first_frame'], PS['n2_sign_differs']],
        [SP['n2_sign_first_frame'], SP['n2_sign_differs']], None, 'second formulas')
    add(run, 'first production frame: signed volume at N2', l6(PS['n2_volume_first_frame']), SP['n2_volume_first_frame'], TOL['vol'], 'second formulas')
    add(run, 'first production frame: ligand RMSD', l6(PS['lig_rmsd_first_frame_A']), SP['lig_rmsd_first_frame_A'], TOL['A'], 'second formulas', 'A')
    add(run, 'campaign-table values: state of the first frame, differs from the start, sign at N2',
        [PS['state_first_frame'], PS['state_differs'], PS['n2_sign_first_frame']],
        [S['table_row']['state_first_frame'], S['table_row']['state_first_frame_differs_from_start'], S['table_row']['n2_sign_first_frame']], None, 'second formulas')
    add(run, 'campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes',
        [N2['fraction_positive'], N2['fraction_positive_window'], len(PN['changes'])],
        [S['table_row']['n2_fraction_positive'], S['table_row']['n2_fraction_positive_last20ns'], S['table_row']['n2_n_lasting_sign_changes']], None, 'second formulas')

    # ---- equilibration (section 9E)
    PE, E = P['equilibration'], S['equilibration']
    fs, ff = adir / 'timeseries_equil_solute.csv', adir / 'timeseries_equil_full.csv'
    add(run, 'equilibration: frames read (files solute_eqNN.dcd present)', PE is not None, E['available'], None, 'construction')
    add(run, 'equilibration: time series written (only with equilibration files)', [PE is not None, PE is not None],
        [fs.is_file(), ff.is_file()], None, 'construction')
    add(run, 'campaign-table value: solute frames of the equilibration', 0 if PE is None else PE['n_frames_solute'],
        S['table_row']['equilibration_frames_solute'], 0, 'construction')
    if PE is not None:
        add(run, 'equilibration: solute frames and full frames read', [PE['n_frames_solute'], PE['n_frames_full']],
            [E['n_frames_solute'], E['n_frames_full']], None, 'construction')
        add(run, 'equilibration: length of the equilibration', PE['total_ps'], E['total_ps'], 1e-9, 'construction', 'ps')
        add(run, 'equilibration: where the stage table came from', 'manifest equilibration.stages' if PE['in_manifest'] else 'protocol (constants of the specification)',
            E['stage_table_source'], None, 'construction')
        mine, theirs = read_csv(root / 'planted' / f'{tag}__equil_solute.csv'), read_csv(fs)
        add(run, 'equilibration: lines of timeseries_equil_solute.csv', len(mine), len(theirs), 0, 'construction')
        with open(fs) as f:
            head = f.readline().strip().split(',')
        add(run, 'equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned)', True,
            'state_assigned' not in head and 'n2_sign_assigned' not in head and 'n2_sign' in head, None, 'construction')
        if mine and len(mine) == len(theirs):
            add(run, 'equilibration, every frame: stage and extension step, frames that differ', 0,
                sum(1 for a, b in zip(mine, theirs) if (a['stage'], a['extension']) != (b['stage'], b['extension'])), 0, 'construction')
            add(run, 'equilibration, every frame: time, largest difference', 0.0,
                max(abs(float(a['time_ns']) - float(b['time_ns'])) for a, b in zip(mine, theirs)), 1e-9, 'construction', 'ns')
            add(run, 'equilibration: no time above zero', True, all(float(b['time_ns']) <= 0.0 for b in theirs), None, 'construction')
            for c, kind in SOLUTE_COLS + (ASN67_COLS if P['with_nd2'] else ()):
                if kind in ('A', 'deg', 'tors', 'vol'):
                    d = np.array([abs(float(a[c]) - float(b[c])) for a, b in zip(mine, theirs)])
                    if kind == 'tors':
                        d = np.minimum(d, 360.0 - d)
                    u = {'A': 'A', 'vol': 'vol'}.get(kind, 'deg')
                    add(run, f'equilibration, every frame: {c}, largest difference', 0.0, float(d.max()), TOL[u], 'second formulas', '' if u == 'vol' else u)
                else:
                    add(run, f'equilibration, every frame: {c}, frames that differ', 0,
                        sum(1 for a, b in zip(mine, theirs) if str(a[c]) != str(b[c])), 0, 'second formulas')
        mine, theirs = read_csv(root / 'planted' / f'{tag}__equil_full.csv'), read_csv(ff)
        add(run, 'equilibration: lines of timeseries_equil_full.csv', len(mine), len(theirs), 0, 'construction')
        if mine and len(mine) == len(theirs):
            add(run, 'equilibration, every full frame: stage and extension step, frames that differ', 0,
                sum(1 for a, b in zip(mine, theirs) if (a['stage'], a['extension']) != (b['stage'], b['extension'])), 0, 'construction')
            add(run, 'equilibration, every full frame: time, largest difference', 0.0,
                max(abs(float(a['time_ns']) - float(b['time_ns'])) for a, b in zip(mine, theirs)), 1e-9, 'construction', 'ns')
            for c, kind in FULL_COLS:
                if kind == 'A':
                    d = np.array([abs(float(a[c]) - float(b[c])) for a, b in zip(mine, theirs)])
                    add(run, f'equilibration, every full frame: {c}, largest difference', 0.0, float(d.max()), TOL['A'], 'second formulas', 'A')
                else:
                    add(run, f'equilibration, every full frame: {c}, frames that differ', 0,
                        sum(1 for a, b in zip(mine, theirs) if str(a[c]) != str(b[c])), 0, 'second formulas')
        es = {(a['stage'], a['extension']): a for a in E['stages']}
        add(run, 'equilibration: the steps with a file, in the order of the run', [a['tag'] for a in PE['stages']], [a['tag'] for a in E['stages']], None, 'construction')
        for a in PE['stages']:
            b = es.get((a['stage'], a['extension']))
            if b is None:
                continue
            k = f'equilibration, stage {a["stage"]}' + (f', extension step {a["extension"]}' if a['extension'] else '')
            dh = b['driver_hydration_check']
            add(run, f'{k}: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it',
                [a['driver_count'], None if a['driver_count'] is None else True],
                [None if dh is None else dh['n_water_O'], None if dh is None else dh['count_agrees']], None, 'construction')
            add(run, f'{k}: start and end', [a['time_start_ps'], a['time_end_ps']], [b['time_start_ps'], b['time_end_ps']], None, 'construction', 'ps')
            add(run, f'{k}: solute and full frames read, and expected',
                [a['n_frames_solute'], a['n_frames_solute'], a['n_frames_full'], a['n_frames_full']],
                [b['solute']['n_frames'], b['solute']['n_frames_expected'], b['full']['n_frames'], b['full']['n_frames_expected']], None, 'construction')
            add(run, f'{k}: frames of state C, R, I; frames with positive and negative sign at N2',
                [a['count_C'], a['count_R'], a['count_I'], a['count_n2_positive'], a['count_n2_negative']],
                [b['count_C'], b['count_R'], b['count_I'], b['count_n2_positive'], b['count_n2_negative']], None, 'second formulas')
            la, lb = a['last_solute_frame'], b['last_solute_frame']
            add(run, f'{k}: last solute frame (frame, time, state, sign at N2)', None if la is None else [la['frame'], l6(la['time_ns']), la['state'], la['n2_sign']],
                None if lb is None else [lb['frame'], lb['time_ns'], lb['state'], lb['n2_sign']], None, 'second formulas')
            if la is not None and lb is not None:
                add(run, f'{k}: last solute frame, ligand RMSD', l6(la['lig_rmsd_A']), lb['lig_rmsd_A'], TOL['A'], 'second formulas', 'A')
                add(run, f'{k}: last solute frame, signed volume at N2', l6(la['n2_volume']), lb['n2_volume'], TOL['vol'], 'second formulas')
            la, lb = a['last_full_frame'], b['last_full_frame']
            add(run, f'{k}: last full frame (frame, time, water oxygens at the crystal ring positions)',
                None if la is None else [la['frame'], l6(la['time_ns']), la['n_water_ring_wide_3p5A']],
                None if lb is None else [lb['frame'], lb['time_ns'], lb['n_water_ring_wide_3p5A']], None, 'second formulas')
            add(run, f'{k}: last full frame, water oxygens in the place of the ring',          # version 1.3
                None if la is None else la['n_water_ring_place'], None if lb is None else lb['n_water_ring_place'], None, 'second formulas')
            if la is not None and lb is not None:
                add(run, f'{k}: last full frame, nearest water oxygen to the place of the ring', l6(la['d_water_ring_place_min_A']),
                    lb['d_water_ring_place_min_A'], TOL['A'], 'second formulas', 'A')
                add(run, f'{k}: last full frame, nearest water oxygen to the seven positions of the wider count', l6(la['d_water_ring_wide_min_A']),
                    lb['d_water_ring_wide_min_A'], TOL['A'], 'second formulas', 'A')
        add(run, 'equilibration: frames of state C, R, I; positive and negative sign at N2',
            [PE['count_C'], PE['count_R'], PE['count_I'], PE['count_n2_positive'], PE['count_n2_negative']],
            [E['count_C'], E['count_R'], E['count_I'], E['count_n2_positive'], E['count_n2_negative']], None, 'second formulas')
        for k in ('first_frame_not_in_state_of_start', 'first_frame_in_opposite_state', 'first_frame_n2_sign_not_of_start'):
            a, b = PE[k], E[k]
            add(run, f'equilibration: {k} (frame, stage, extension step, time)', None if a is None else [a['frame'], a['stage'], a['extension'], l6(a['time_ns'])],
                None if b is None else [b['frame'], b['stage'], b['extension'], b['time_ns']], None, 'second formulas')
        add(run, 'equilibration: state and sign at N2 of the last frame', [PE['state_last_frame'], PE['n2_sign_last_frame']],
            [E['state_last_frame'], E['n2_sign_last_frame']], None, 'second formulas')
        add(run, 'equilibration: time of the last solute frame', l6(PE['time_last_ns']), E['time_last_solute_frame_ns'], 1e-9 if PE['time_last_ns'] is not None else None, 'construction', 'ns')
        if PE['n_frames_solute']:
            add(run, 'equilibration: frames in which the ligand and the zinc were in another periodic image',
                [PE['frames_solute_moved']['ligand'], PE['frames_solute_moved']['zinc']],
                [E['imaging_solute']['frames_ligand_moved'], E['imaging_solute']['frames_zinc_moved']], None, 'construction')

    # ---- the design typed by hand (V12_DESIGN)
    D = P['design_v12']
    if D is None:
        return
    dn = D['n2']
    add(run, 'by design: frames with positive and with negative signed volume at N2', [dn['count_positive'], dn['count_negative']],
        [N2['count_positive'], N2['count_negative']], None, 'construction')
    add(run, 'by design: positive and negative in the last 20 ns', [dn['count_positive_window'], dn['count_negative_window']],
        [N2['count_positive_window'], N2['count_negative_window']], None, 'construction')
    for k in 'CRI':
        b = N2['by_state'][k]
        add(run, f'by design, frames of state {k}: number, positive, negative', dn['by_state'][k], [b['n'], b['count_positive'], b['count_negative']], None, 'construction')
    add(run, 'by design: stretches of one sign, and those of at least 1 ns', [dn['n_stretches'], dn['n_long_stretches']],
        [N2['n_stretches'], N2['n_long_stretches']], None, 'construction')
    add(run, 'by design: times of the sign changes that last at least 1 ns', dn['change_times_ns'],
        [c['time_ns'] for c in N2['lasting_sign_changes']], None, 'construction', 'ns')
    add(run, 'by design: sign after each of these changes', dn['change_to_sign'], [c['to_sign'] for c in N2['lasting_sign_changes']], None, 'construction')
    dh = D['hb']
    add(run, 'by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond)',
        [dh['n_before'], dh['count_before'], dh['n_after'], dh['count_after']], [H['n_before'], H['count_before'], H['n_after'], H['count_after']], None, 'construction')
    dg = D['glu']
    add(run, 'by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A', [dg['count_OE1'], dg['count_OE2'], dg['count_either']],
        [FR['fifth_E106OE1']['count'], FR['fifth_E106OE2']['count'], FR['fifth_E106OE']['count']], None, 'construction')
    if dg['min_OE1'] is not None:
        add(run, 'by design: smallest distance zinc to OE1 of Glu106', dg['min_OE1'], ST['d_Zn_E106OE1_A']['min'], TOL['A'], 'construction', 'A')
    if dg['min_OE2'] is not None:
        add(run, 'by design: smallest distance zinc to OE2 of Glu106', dg['min_OE2'], ST['d_Zn_E106OE2_A']['min'], TOL['A'], 'construction', 'A')
    if D['ring'] is not None:
        full = read_csv(adir / 'timeseries_full.csv')
        add(run, 'by design: water oxygens at the crystal ring positions, full frames that differ', 0,
            sum(1 for a, b in zip(D['ring']['counts'], full) if a != int(b['n_water_ring_wide_3p5A'])) + abs(len(full) - len(D['ring']['counts'])), 0, 'construction')
        add(run, 'by design: mean number of water oxygens at the crystal ring positions', l6(D['ring']['mean']), ST['n_water_ring_wide_3p5A']['mean'], 1e-6, 'construction')
        add(run, 'by design: smallest distance of a water oxygen to the crystal ring positions', D['ring']['min_distance'], ST['d_water_ring_wide_min_A']['min'], TOL['A'], 'construction', 'A')
    if D.get('place') is not None:                            # version 1.3: the count the driver acts on
        full = read_csv(adir / 'timeseries_full.csv')
        dpl = D['place']
        add(run, 'by design: water oxygens in the place of the ring in every full frame', dpl['counts'],
            [int(b['n_water_ring_place']) for b in full], None, 'construction')
        add(run, 'by design: mean number of water oxygens in the place of the ring', l6(dpl['mean']), ST['n_water_ring_place']['mean'], 1e-6, 'construction')
        add(run, 'by design: smallest distance of a water oxygen to the place of the ring', dpl['min_distance'], ST['d_water_ring_place_min_A']['min'], TOL['A'], 'construction', 'A')
        if 'nearest' in dpl:
            add(run, 'by design: nearest water oxygen to the place of the ring in every full frame, largest difference', 0.0,
                max(abs(a - float(b['d_water_ring_place_min_A'])) for a, b in zip(dpl['nearest'], full)), TOL['A'], 'construction', 'A')
    if D.get('settings') is not None:
        dst, WP = D['settings'], S['water_ring_place']
        add(run, 'by design: site atoms and cut-off of the count that the driver acts on', [dst['site_atoms'], dst['water_cutoff_A']],
            [WP['site_atoms'], WP['water_cutoff_A']], None, 'construction', 'A')
        add(run, 'by design: where these settings came from, and whether they are the defaults of the specification',
            [dst['source'], dst['settings_are_the_defaults']], [WP['source'], WP['settings_are_the_defaults']], None, 'construction')
    ds = D['start']
    add(run, 'by design: state of the first production frame, differs from the start; sign at N2, differs from the start',
        [ds['state_first_frame'], ds['state_differs'], ds['n2_sign_first_frame'], ds['n2_sign_differs']],
        [SP['state_first_frame'], SP['state_differs'], SP['n2_sign_first_frame'], SP['n2_sign_differs']], None, 'construction')
    de = D['equil']
    add(run, 'by design: equilibration frames read', de['available'], E['available'], None, 'construction')
    if de['available']:
        add(run, 'by design: equilibration, solute frames, full frames, length (ps)', [de['n_solute'], de['n_full'], de['total_ps']],
            [E['n_frames_solute'], E['n_frames_full'], E['total_ps']], None, 'construction')
        add(run, 'by design: equilibration, where the stage table came from', de['source'], E['stage_table_source'], None, 'construction')
        add(run, 'by design: equilibration, solute and full frames of each stage', de['frames_per_stage'],
            [[a['solute']['n_frames'], a['full']['n_frames']] for a in E['stages']], None, 'construction')
        eqs = read_csv(fs)
        add(run, 'by design: equilibration, time of the first and of the last solute frame', [de['time_first_ns'], de['time_last_ns']],
            [float(eqs[0]['time_ns']), float(eqs[-1]['time_ns'])], None, 'construction', 'ns')
        if 'times_ns' in de:
            add(run, 'by design: equilibration, the times of all solute frames', de['times_ns'], [float(r['time_ns']) for r in eqs], None, 'construction', 'ns')
        if 'time_frame_10_ns' in de:
            add(run, 'by design: equilibration, time of solute frame 10 (the first of stage 2)', [de['time_frame_10_ns'], '2'],
                [float(eqs[10]['time_ns']), eqs[10]['stage']], None, 'construction', 'ns')
        if de.get('ring_counts') is not None:
            add(run, 'by design: equilibration, water oxygens at the crystal ring positions in every full frame', de['ring_counts'],
                [int(r['n_water_ring_wide_3p5A']) for r in read_csv(ff)], None, 'construction')
        if de.get('ring_last_frame_of_stage_3') is not None:
            add(run, 'by design: equilibration, water oxygens at the crystal ring positions at the end of stage 3 (the held part)',
                de['ring_last_frame_of_stage_3'], E['stages'][2]['last_full_frame']['n_water_ring_wide_3p5A'], 0, 'construction')
        if de.get('place_counts') is not None:                # version 1.3
            add(run, 'by design: equilibration, water oxygens in the place of the ring in every full frame', de['place_counts'],
                [int(r['n_water_ring_place']) for r in read_csv(ff)], None, 'construction')
        if de.get('place_nearest') is not None:
            add(run, 'by design: equilibration, nearest water oxygen to the place of the ring in every full frame, largest difference', 0.0,
                max(abs(a - float(r['d_water_ring_place_min_A'])) for a, r in zip(de['place_nearest'], read_csv(ff))), TOL['A'], 'construction', 'A')
        if de.get('place_last_frame_of_stage_3') is not None:
            add(run, 'by design: equilibration, water oxygens in the place of the ring at the end of stage 3',
                de['place_last_frame_of_stage_3'], E['stages'][2]['last_full_frame']['n_water_ring_place'], 0, 'construction')
        for k, key in (('first_not_start', 'first_frame_not_in_state_of_start'), ('first_opposite', 'first_frame_in_opposite_state'),
                       ('first_sign', 'first_frame_n2_sign_not_of_start')):
            b = E[key]
            add(run, f'by design: equilibration, {key} (frame, stage, extension step, time)', de[k],
                None if b is None else [b['frame'], b['stage'], b['extension'], b['time_ns']], None, 'construction')
        add(run, 'by design: equilibration, state and sign at N2 of the last frame', [de['state_last'], de['sign_last']],
            [E['state_last_frame'], E['n2_sign_last_frame']], None, 'construction')
        if 'steps' in de:                                     # a run whose held part was extended
            add(run, 'by design: equilibration, the steps in the order of the run', de['steps'], [a['tag'] for a in E['stages']], None, 'construction')
            add(run, 'by design: equilibration, stage, extension step and time of solute frame 70 (the first of eq03x1)', de['frame_70'],
                [int(eqs[70]['stage']), int(eqs[70]['extension']), float(eqs[70]['time_ns'])], None, 'construction')
            add(run, f'by design: equilibration, stage, extension step and time of solute frame 170 (the first of {de.get("frame_170_step", "eq04")})', de['frame_170'],
                [int(eqs[170]['stage']), int(eqs[170]['extension']), float(eqs[170]['time_ns'])], None, 'construction')
            add(run, 'by design: equilibration, water oxygens in the place of the ring (the count the driver acts on) at the end of '
                + ', '.join(de['ring_at_the_end_of']), list(de['ring_at_the_end_of'].values()),
                [a['last_full_frame']['n_water_ring_place'] for a in E['stages'] if a['tag'] in de['ring_at_the_end_of']], None, 'construction')
            add(run, 'by design: equilibration, these steps were compared with the driver (their last full frame is at the end of the step)',
                [True] * len(de['ring_at_the_end_of']),
                [a['driver_hydration_check']['compared'] for a in E['stages'] if a['tag'] in de['ring_at_the_end_of']], None, 'construction')
            add(run, 'by design: equilibration, the count of the driver at the end of these steps agrees with the analysis', de['driver_counts_agree'],
                [a['driver_hydration_check']['count_agrees'] for a in E['stages'] if a['tag'] in de['ring_at_the_end_of']], None, 'construction')
            add(run, 'by design: equilibration, the nearest water oxygen of the driver at the end of these steps agrees with the analysis', de['driver_counts_agree'],
                [a['driver_hydration_check']['nearest_agrees'] for a in E['stages'] if a['tag'] in de['ring_at_the_end_of']], None, 'construction')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--out', metavar='DIR', help='write the synthetic test below DIR')
    g.add_argument('--compare', metavar='DIR', help='compare the analysis below DIR with the planted values')
    a = ap.parse_args()
    return generate(a.out) if a.out else compare(a.compare)


if __name__ == '__main__':
    sys.exit(main())
