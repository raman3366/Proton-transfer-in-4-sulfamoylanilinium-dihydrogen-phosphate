"""Crystal environment of SAPH from the deposited CIF (gemmi), for fig01_scheme.py.

    ions(window)      -> list of Ion: every image of the cation and the anion (space-group operations of Pc x lattice
                         translations -2..2) whose centroid lies in the fractional window
    hbonds(ions)      -> list of HBondInst: every instance, between the ions given, of the nine hydrogen bonds of the
                         deposited table (D, H, A labels), identified by distance: D...A equal to the deposited value
                         (facts xray.hb.<D>_<H>_<A>.DA) within 0.002 A, H...A within 0.02 A of the deposited value,
                         angle D-H...A > 120 deg.  Nothing about which image accepts is hard-coded; the symmetry
                         operation found for each bond is checked against the SHELX code of the CIF in main().
    HB                -> the nine bonds in CIF order (numbered 1-9 in the figure, as in Table 4)
Coordinates: orthogonal Angstrom from gemmi's orthogonalisation matrix (x || a, y || b, z || c*).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import gemmi

import _saph as S

HB = [('N1', 'H1A', 'O5'), ('N1', 'H1B', 'O6'), ('N1', 'H1C', 'O6'), ('N2', 'H2A', 'O3'), ('N2', 'H2B', 'O2'),
      ('O3', 'H3A', 'O1'), ('O4', 'H4A', 'O5'), ('C3', 'H3', 'O6'), ('C6', 'H6', 'O4')]
# SHELX symmetry codes of the acceptor in the CIF hydrogen-bond table (operation number, translation + 555)
CIF_CODE = {HB[0]: '1_554', HB[1]: '2_554', HB[2]: '2_654', HB[3]: '2_564', HB[4]: '1_655', HB[5]: '1_455',
            HB[6]: '1_455', HB[7]: '1_655', HB[8]: '1_554'}
IN_PAIR = {HB[0], HB[8]}          # the two bonds kept in the 27-atom ion pair (they close the R2,2(8) ring)


@dataclass
class Ion:
    kind: str                     # 'cation' / 'anion'
    op: int                       # 1 = x, y, z ; 2 = x, -y, z + 1/2
    tr: tuple                     # lattice translation
    labels: list
    symbols: list
    F: np.ndarray                 # fractional
    X: np.ndarray                 # orthogonal, Angstrom
    key: tuple = field(default=None)

    def at(self, label):
        return self.X[self.labels.index(label)]


@dataclass
class HBondInst:
    n: int                        # 1..9, CIF order
    d: tuple                      # (ion index, label)
    h: tuple
    a: tuple
    da: float
    ha: float
    ang: float


@lru_cache(None)
def structure():
    st = gemmi.read_small_structure(str(S.CIF))
    M = np.array(st.cell.orth.mat.tolist())
    labels = [s.label for s in st.sites]
    sym = [s.type_symbol.strip().capitalize() for s in st.sites]
    F = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in st.sites])
    sg = gemmi.find_spacegroup_by_name(st.spacegroup_hm)
    ops = []
    for op in sg.operations():
        R = np.array(op.rot, float) / op.DEN
        t = np.array(op.tran, float) / op.DEN
        ops.append((op.triplet(), R, t))
    assert [o[0] for o in ops] == ['x,y,z', 'x,-y,z+1/2'], [o[0] for o in ops]
    return dict(st=st, M=M, labels=labels, sym=sym, F=F, ops=ops, cell=st.cell)


def _ion(kind, opn, tr):
    s = structure()
    labs = S.CATION if kind == 'cation' else S.ANION
    L = {l: i for i, l in enumerate(s['labels'])}
    idx = [L[l] for l in labs]
    _, R, t = s['ops'][opn - 1]
    Fm = s['F'][idx] @ R.T + t + np.array(tr, float)
    return Ion(kind, opn, tuple(tr), list(labs), [s['sym'][i] for i in idx], Fm, Fm @ s['M'].T,
               key=(kind, opn, tuple(tr)))


@lru_cache(None)
def all_ions(n=2):
    out = []
    for kind in ('cation', 'anion'):
        for opn in (1, 2):
            for tr in itertools.product(range(-n, n + 1), repeat=3):
                out.append(_ion(kind, opn, tr))
    return out


def ions(window):
    """window: ((x0, x1), (y0, y1), (z0, z1)) fractional limits for the ion centroid."""
    out = []
    for ion in all_ions():
        c = ion.F.mean(0)
        if all(lo <= v < hi for v, (lo, hi) in zip(c, window)):
            out.append(ion)
    return out


def ion(kind, opn=1, tr=(0, 0, 0)):
    return _ion(kind, opn, tr)


def _angle(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def hbonds(ion_list, only=None, donors=None):
    """Instances of the nine deposited hydrogen bonds between ions of ion_list (see module docstring).
    donors: optional indices of the ions allowed to donate (default all)."""
    out = []
    donors = range(len(ion_list)) if donors is None else donors
    for n, (d, h, a) in enumerate(HB, start=1):
        if only is not None and n not in only:
            continue
        DA = S.F(f'xray.hb.{d}_{h}_{a}.DA')
        HA = S.F(f'xray.hb.{d}_{h}_{a}.HA')
        acc = [(j, x.at(a)) for j, x in enumerate(ion_list) if a in x.labels]
        if not acc:
            continue
        J = np.array([j for j, _ in acc])
        XA = np.array([p for _, p in acc])
        for i in donors:
            di = ion_list[i]
            if d not in di.labels:
                continue
            Xd, Xh = di.at(d), di.at(h)
            dd = np.linalg.norm(XA - Xd, axis=1)
            for k in np.nonzero(np.abs(dd - DA) <= 0.002)[0]:
                j = int(J[k])
                if j == i:
                    continue
                Xa = XA[k]
                ha = float(np.linalg.norm(Xh - Xa))
                ang = _angle(Xd, Xh, Xa)
                if abs(ha - HA) < 0.02 and ang > 120.0:
                    out.append(HBondInst(n, (i, d), (i, h), (j, a), float(dd[k]), ha, ang))
    return out


def shelx_code(donor_ion, acceptor_ion):
    """SHELX code n_klm of the acceptor ion relative to the donor ion (the donor must be an identity image)."""
    assert donor_ion.op == 1
    t = np.array(acceptor_ion.tr) - np.array(donor_ion.tr)
    return f'{acceptor_ion.op}_{t[0] + 5}{t[1] + 5}{t[2] + 5}'


def check_codes():
    """Every deposited bond is found from the asymmetric-unit donor, with the CIF's symmetry code, exactly once."""
    pool = all_ions(2)
    for n, (d, h, a) in enumerate(HB, start=1):
        donor = ion('cation' if d in S.CATION else 'anion')
        lst = [donor] + [p for p in pool if p.key != donor.key]
        found = [x for x in hbonds(lst, only=[n], donors=[0])]
        codes = sorted(shelx_code(donor, lst[x.a[0]]) for x in found)
        assert codes == [CIF_CODE[(d, h, a)]], (n, d, h, a, codes)
    return True
