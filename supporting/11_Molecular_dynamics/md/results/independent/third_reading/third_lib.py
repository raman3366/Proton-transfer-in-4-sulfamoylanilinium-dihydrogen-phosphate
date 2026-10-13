#!/usr/bin/env python
"""Blind third reading of the MD campaign: library.

Own code throughout: own PDB parsing (fixed columns), own DCD reader (numpy),
own imaging, Kabsch fit, states, visits.  Nothing is imported from
md/analysis/analyse_md.py or md/analysis/recheck_md.py, and none of their
outputs is read.
"""
import csv
import itertools
import json
import os
import struct

import numpy as np

REPO = '/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate'
RUNS = '/root/cqc-sfa-md/runs'
SYS = os.path.join(REPO, 'md', 'system')

M1 = 6.7          # midpoint of d1 = |N1 - CZ(Phe131)|, plan entry 5 Oct 8:40 PM
M2 = 5.2          # midpoint of d2 = |N1 - OD1(Asn67)|
WINDOW = 2000     # last 20 ns, a frame every 10 ps
LONG = 100        # a long visit: at least 100 frames (1 ns)
DT_PS = 10.0      # interval of the solute frames

LIG_HEAVY = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S']
LABELLINGS = [
    ('identity', {}),
    ('ring_flip', {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3'}),
    ('oxygen_swap', {'O1': 'O2', 'O2': 'O1'}),
    ('ring_flip+oxygen_swap', {'C2': 'C6', 'C6': 'C2', 'C3': 'C5', 'C5': 'C3',
                               'O1': 'O2', 'O2': 'O1'}),
]

RUN_LIST = [
    ('bonded_2e', 'X1', 'X', 'production'), ('bonded_2e', 'X2', 'X', 'production'),
    ('bonded_2e', 'X3', 'X', 'production'), ('bonded_2e', 'N1', 'N', 'production'),
    ('bonded_2e', 'N2', 'N', 'production'), ('bonded_2e', 'N3', 'N', 'production'),
    ('bonded_2e', 'R1', 'R', 'production'), ('bonded_2e', 'R2', 'R', 'production'),
    ('bonded_2e', 'R3', 'R', 'production'),
    ('bonded_2s', 'R1', 'R', 'sensitivity'), ('bonded_2e_n2fit', 'R1', 'R', 'sensitivity'),
]


# ---------------------------------------------------------------- PDB
def read_pdb(fn):
    """Fixed-column reader.  Returns list of dicts and CRYST1 (or None)."""
    atoms = []
    cryst = None
    with open(fn) as f:
        for line in f:
            rec = line[:6]
            if rec == 'CRYST1':
                cryst = [float(line[6:15]), float(line[15:24]), float(line[24:33]),
                         float(line[33:40]), float(line[40:47]), float(line[47:54])]
            elif rec in ('ATOM  ', 'HETATM'):
                atoms.append(dict(
                    name=line[12:16].strip(), resname=line[17:20].strip(),
                    chain=line[21], resseq=int(line[22:26]),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return atoms, cryst


def residue_map():
    """deposited residue number -> tleap index (md/system/residue_map.csv)."""
    out = {}
    with open(os.path.join(SYS, 'residue_map.csv')) as f:
        for row in csv.DictReader(f):
            out[int(row['deposited_number'])] = (int(row['tleap_index']), row['deposited_name'],
                                                 row['name_md'])
    return out


def build_reference():
    """Reference of the fit and of the ligand RMSD, built by this reader from
    md/system/protein_md.pdb (CA atoms, 6RL9 frame) and md/system/ligand_X.pdb
    (crystal pose, SAN 304 of 6RL9)."""
    prot, _ = read_pdb(os.path.join(SYS, 'protein_md.pdb'))
    ca = [(a['resseq'], a['resname'], a['xyz']) for a in prot if a['name'] == 'CA']
    lig, _ = read_pdb(os.path.join(SYS, 'ligand_X.pdb'))
    ligd = {a['name']: a['xyz'] for a in lig}
    return dict(ca_resseq=[c[0] for c in ca], ca_resname=[c[1] for c in ca],
                ca=np.array([c[2] for c in ca], float),
                lig={n: np.array(ligd[n], float) for n in LIG_HEAVY})


def select(solute_pdb):
    """Atom indices (0-based, order of solute.pdb = order of solute_*.dcd)."""
    atoms, cryst = read_pdb(solute_pdb)
    rmap = residue_map()

    def one(dep, name, allowed):
        tl = rmap[dep][0]
        hits = [i for i, a in enumerate(atoms) if a['resseq'] == tl and a['name'] == name]
        if len(hits) != 1:
            raise RuntimeError('%s of residue %d: %d hits' % (name, dep, len(hits)))
        if atoms[hits[0]]['resname'] not in allowed:
            raise RuntimeError('%s of residue %d: residue name %s' %
                               (name, dep, atoms[hits[0]]['resname']))
        return hits[0]

    zn = [i for i, a in enumerate(atoms) if a['name'] == 'ZN']
    if len(zn) != 1:
        raise RuntimeError('zinc: %d hits' % len(zn))
    zres = atoms[zn[0]]['resseq']
    lig_idx = [i for i, a in enumerate(atoms) if a['resname'] == 'SN1']
    lig_res = {atoms[i]['resseq'] for i in lig_idx}
    if len(lig_res) != 1 or len(lig_idx) != 18:
        raise RuntimeError('ligand: %d atoms in residues %s' % (len(lig_idx), lig_res))
    lig = {atoms[i]['name']: i for i in lig_idx}
    for n in LIG_HEAVY + ['HN21']:
        if n not in lig:
            raise RuntimeError('ligand atom %s missing' % n)
    lres = lig_res.pop()
    ca = [i for i, a in enumerate(atoms)
          if a['name'] == 'CA' and a['resseq'] not in (zres, lres)]
    sel = dict(
        n_atoms=len(atoms), cryst=cryst,
        zn=zn[0], lig_all=np.array(sorted(lig_idx)), lig=lig, ca=np.array(ca),
        ca_resseq=[atoms[i]['resseq'] for i in ca],
        ca_resname=[atoms[i]['resname'] for i in ca],
        h94_NE2=one(94, 'NE2', ('HD1', 'HID', 'HIS')),
        h96_NE2=one(96, 'NE2', ('HD2', 'HID', 'HIS')),
        h119_ND1=one(119, 'ND1', ('HE1', 'HIE', 'HIS')),
        t199_OG1=one(199, 'OG1', ('THR',)),
        f131_CZ=one(131, 'CZ', ('PHE',)),
        n67_OD1=one(67, 'OD1', ('ASN',)),
        xyz=np.array([a['xyz'] for a in atoms], float),
    )
    return sel


# ---------------------------------------------------------------- DCD
def read_dcd(fn, n_atoms):
    """Own reader of the DCD files of the driver (little endian, CHARMM layout
    with a unit-cell record of six doubles per frame).  Complete frames only.
    Returns coordinates (F, n_atoms, 3) float64 in Angstrom, lattice vectors
    (F, 3, 3) or None, and header facts."""
    with open(fn, 'rb') as f:
        raw = f.read()
    if len(raw) < 92 or struct.unpack('<i', raw[:4])[0] != 84 or raw[4:8] != b'CORD':
        raise RuntimeError('not a little-endian CORD DCD: ' + fn)
    ic = struct.unpack('<9i', raw[8:44])
    delta = struct.unpack('<f', raw[44:48])[0]
    ic2 = struct.unpack('<10i', raw[48:88])
    has_cell = ic2[0] == 1
    off = 92
    n = struct.unpack('<i', raw[off:off + 4])[0]
    title = raw[off + 8:off + 4 + n].decode('ascii', 'replace')
    off += 4 + n + 4
    m0, nat, m1 = struct.unpack('<3i', raw[off:off + 12])
    if (m0, m1) != (4, 4) or nat != n_atoms:
        raise RuntimeError('atom count %d in %s, expected %d' % (nat, fn, n_atoms))
    off += 12
    fields = []
    if has_cell:
        fields += [('c0', '<i4'), ('cell', '<f8', (6,)), ('c1', '<i4')]
    for ax in 'xyz':
        fields += [(ax + '0', '<i4'), (ax, '<f4', (nat,)), (ax + '1', '<i4')]
    dt = np.dtype(fields)
    nfr = (len(raw) - off) // dt.itemsize
    left = (len(raw) - off) - nfr * dt.itemsize
    arr = np.frombuffer(raw, dtype=dt, count=nfr, offset=off)
    for ax in 'xyz':
        if nfr and not (np.all(arr[ax + '0'] == 4 * nat) and np.all(arr[ax + '1'] == 4 * nat)):
            raise RuntimeError('record markers wrong in ' + fn)
    xyz = np.empty((nfr, nat, 3), float)
    xyz[:, :, 0] = arr['x']
    xyz[:, :, 1] = arr['y']
    xyz[:, :, 2] = arr['z']
    lat = None
    if has_cell:
        if nfr and not (np.all(arr['c0'] == 48) and np.all(arr['c1'] == 48)):
            raise RuntimeError('cell record markers wrong in ' + fn)
        cell = arr['cell']                       # a, cos(gamma), b, cos(beta), cos(alpha), c
        a, cg, b, cb, ca_, c = [cell[:, i].astype(float) for i in range(6)]
        if nfr and np.any(np.abs(np.stack([cg, cb, ca_])) > 1.0):
            # angles in degrees instead of cosines
            cg, cb, ca_ = [np.cos(np.radians(v)) for v in (cg, cb, ca_)]
        sg = np.sqrt(1.0 - cg ** 2)
        lat = np.zeros((nfr, 3, 3))
        lat[:, 0, 0] = a
        lat[:, 1, 0] = b * cg
        lat[:, 1, 1] = b * sg
        lat[:, 2, 0] = c * cb
        lat[:, 2, 1] = c * (ca_ - cb * cg) / sg
        lat[:, 2, 2] = np.sqrt(np.maximum(c ** 2 - lat[:, 2, 0] ** 2 - lat[:, 2, 1] ** 2, 0.0))
    head = dict(nset_header=ic[0], first_step=ic[1], interval_steps=ic[2], last_step=ic[3],
                delta_akma=delta, has_cell=has_cell, n_frames_read=int(nfr),
                bytes_left_over=int(left), title=title.strip())
    return xyz, lat, head


# ---------------------------------------------------------------- geometry
_COMBOS = np.array(list(itertools.product(range(-2, 3), repeat=3)), float)   # 125 cells


def best_lattice_shift(diff, lat):
    """For each frame the lattice vector t (searched over 125 cells) that makes
    |diff + t| smallest.  diff (F,3), lat (F,3,3) rows a, b, c."""
    shifts = np.einsum('ci,fij->fcj', _COMBOS, lat)          # (F,125,3)
    d2 = ((diff[:, None, :] + shifts) ** 2).sum(-1)
    k = d2.argmin(1)
    f = np.arange(len(diff))
    return shifts[f, k], _COMBOS[k]


def dist(a, b):
    return np.sqrt(((a - b) ** 2).sum(-1))


def angle_deg(a, b, c):
    """Angle at b."""
    u = a - b
    v = c - b
    cs = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    return np.degrees(np.arccos(np.clip(cs, -1.0, 1.0)))


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def kabsch(P, Q):
    """Rotation R (row-vector convention, no reflection) and means such that
    (P - mean(P)) @ R + mean(Q) is the least-squares superposition of P on Q.
    P (F,n,3), Q (n,3).  Equal weights."""
    pm = P.mean(1, keepdims=True)
    qm = Q.mean(0)
    Pc = P - pm
    Qc = Q - qm
    H = np.einsum('fni,nj->fij', Pc, Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.zeros_like(H)
    D[:, 0, 0] = 1.0
    D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    R = U @ D @ Vt
    fitted = Pc @ R
    rmsd = np.sqrt(((fitted - Qc) ** 2).sum(-1).mean(-1))
    return R, pm, qm, rmsd


def frame_quantities(xyz, lat, sel, ref):
    """All per-frame numbers of a block of solute frames."""
    F = xyz.shape[0]
    xyz = xyz.copy()
    out = {}
    # 1. images: zinc to NE2(His94), then the whole ligand to the zinc
    if lat is not None:
        t, k = best_lattice_shift(xyz[:, sel['zn']] - xyz[:, sel['h94_NE2']], lat)
        xyz[:, sel['zn']] += t
        out['zn_moved'] = np.any(k != 0, axis=1)
        t2, k2 = best_lattice_shift(xyz[:, sel['lig']['N2']] - xyz[:, sel['zn']], lat)
        xyz[:, sel['lig_all']] += t2[:, None, :]
        out['lig_moved'] = np.any(k2 != 0, axis=1)
    else:
        out['zn_moved'] = np.zeros(F, bool)
        out['lig_moved'] = np.zeros(F, bool)
    L = {n: xyz[:, i] for n, i in sel['lig'].items()}
    zn = xyz[:, sel['zn']]
    # 2. orientation coordinate
    out['d1'] = dist(L['N1'], xyz[:, sel['f131_CZ']])
    out['d2'] = dist(L['N1'], xyz[:, sel['n67_OD1']])
    if lat is not None:
        # diagnostic only: the same two distances as minimum-image distances
        for key, at in (('d1_mi', 'f131_CZ'), ('d2_mi', 'n67_OD1')):
            dd = L['N1'] - xyz[:, sel[at]]
            t3, _ = best_lattice_shift(dd, lat)
            out[key] = np.linalg.norm(dd + t3, axis=1)
    # 3. zinc site
    out['zn_n2'] = dist(zn, L['N2'])
    out['zn_h94'] = dist(zn, xyz[:, sel['h94_NE2']])
    out['zn_h96'] = dist(zn, xyz[:, sel['h96_NE2']])
    out['zn_h119'] = dist(zn, xyz[:, sel['h119_ND1']])
    out['zn_o1'] = dist(zn, L['O1'])
    out['zn_o2'] = dist(zn, L['O2'])
    out['zn_omin'] = np.minimum(out['zn_o1'], out['zn_o2'])
    # 4. hydrogen bond N2-H...OG1(Thr199)
    og1 = xyz[:, sel['t199_OG1']]
    out['n2_og1'] = dist(L['N2'], og1)
    out['ang_n2_h_og1'] = angle_deg(L['N2'], L['HN21'], og1)
    out['hb'] = (out['n2_og1'] <= 3.5) & (out['ang_n2_h_og1'] >= 135.0)
    # 5. configuration at N2: signed volume of the unit vectors N2->Zn, N2->H, N2->S
    u = unit(zn - L['N2'])
    v = unit(L['HN21'] - L['N2'])
    w = unit(L['S'] - L['N2'])
    out['n2_vol'] = (u * np.cross(v, w)).sum(-1)
    # 6. fit of the protein CA atoms on the reference, ligand RMSD to the crystal pose
    R, pm, qm, ca_rmsd = kabsch(xyz[:, sel['ca']], ref['ca'])
    out['ca_rmsd'] = ca_rmsd
    ligfit = {n: ((L[n][:, None, :] - pm) @ R)[:, 0, :] + qm for n in LIG_HEAVY}
    rm = np.empty((F, 4))
    for j, (_, perm) in enumerate(LABELLINGS):
        s = np.zeros(F)
        for n in LIG_HEAVY:
            s += ((ligfit[n] - ref['lig'][perm.get(n, n)]) ** 2).sum(-1)
        rm[:, j] = np.sqrt(s / len(LIG_HEAVY))
    out['lig_rmsd_all'] = rm
    out['lig_rmsd_labelling'] = rm.argmin(1)        # the lowest number if two are equal
    out['lig_rmsd'] = rm.min(1)
    out['lig_rmsd_identity'] = rm[:, 0]
    return out


# ---------------------------------------------------------------- states
def states_of(d1, d2):
    """'C' crystal-like, 'R' rotated-like, 'I' everything else (strict inequalities)."""
    st = np.full(len(d1), 'I', dtype='<U1')
    st[(d1 < M1) & (d2 > M2)] = 'C'
    st[(d1 > M1) & (d2 < M2)] = 'R'
    return st


def assign(seq, start, neutral):
    """Value of the frame if it is not the neutral one, else the value before;
    before the first frame the value of the start."""
    out = []
    last = start
    for s in seq:
        if s != neutral:
            last = s
        out.append(last)
    return out


def visits(assigned):
    """Maximal stretches of equal assigned value: list of (value, first, length)."""
    v = []
    i = 0
    n = len(assigned)
    while i < n:
        j = i
        while j + 1 < n and assigned[j + 1] == assigned[i]:
            j += 1
        v.append((assigned[i], i, j - i + 1))
        i = j + 1
    return v


def passage(seq, start, opposite, neutral):
    """Assigned values, visits, long visits, first long visit of the opposite
    value, transitions in (start, long visits...)."""
    a = assign(list(seq), start, neutral)
    vis = visits(a)
    longv = [x for x in vis if x[2] >= LONG]
    first = None
    for x in longv:
        if x[0] == opposite:
            first = x
            break
    chain = [start] + [x[0] for x in longv]
    ntr = sum(1 for i in range(1, len(chain)) if chain[i] != chain[i - 1])
    return a, vis, longv, first, ntr


def counts(st):
    n = len(st)
    c = {s: int(np.sum(st == s)) for s in 'CRI'}
    c['n'] = int(n)
    c['fraction'] = {s: (c[s] / n if n else None) for s in 'CRI'}
    return c


def verdict(st, start_state):
    opp = 'R' if start_state == 'C' else 'C'
    n = len(st)
    w = min(n, WINDOW)
    win = st[n - w:]
    c_opp = int(np.sum(win == opp))
    c_start = int(np.sum(win == start_state))
    if 10 * c_opp >= 9 * w:
        v = 'CONVERTED'
    elif 10 * c_start >= 9 * w:
        v = 'STAYED'
    else:
        v = 'MIXED'
    return v, w, c_start, c_opp


def segment_files(folder, kind='solute'):
    segs = sorted(f for f in os.listdir(folder)
                  if f.startswith(kind + '_seg') and f.endswith('.dcd'))
    return segs


def load_manifest(folder):
    with open(os.path.join(folder, 'manifest.json')) as f:
        return json.load(f)
