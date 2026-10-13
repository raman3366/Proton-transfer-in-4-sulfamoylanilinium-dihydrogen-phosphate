"""Shared helpers for the independent re-derivation scripts (paper/main/data/derive/d*.py).

These scripts recompute the headline numbers from the raw result files with their own code (own geometry,
own unit conversions, own parsers); they never import dft/ or dock/ code.  Each writes one JSON into
paper/main/data/derived/.  Inputs are read only; nothing under molecules/, dock/ or the archive is modified.
"""
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
RES = REPO / 'molecules/sfa/results'
ARCH = Path('/root/cqc-sfa-archive/vm/cqc-sfa')
ARES = ARCH / 'molecules/sfa/results'
DOCK = REPO / 'dock'
OUT = REPO / 'paper/main/data/derived'
CIF = REPO / 'e-70-0o997-sup1.cif'

# CODATA 2018-based conversions (the same values the project tools use)
HA_KCAL = 627.5095
HA_KJ = 2625.4996
HA_EV = 27.211386
RY_EV = 13.605693
EV_KJ = 96.485332
BOHR_A = 0.529177210903
CM1_KCAL = 2.859144e-3        # 1 cm-1 in kcal/mol (h c N_A)
AU_DEBYE = 2.541746473

CATION = ["C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S1",
          "H2", "H3", "H5", "H6", "H1A", "H1B", "H1C", "H2A", "H2B"]
ANION = ["P1", "O3", "O4", "O5", "O6", "H3A", "H4A"]
PAIR_LABELS = CATION + ANION
IDX = {l: i for i, l in enumerate(PAIR_LABELS)}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.array(xyz), lines[1]


def read_extxyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    head = lines[1]
    lat = head.split('Lattice="')[1].split('"')[0].split()
    cell = np.array([float(x) for x in lat]).reshape(3, 3)
    syms, xyz, labels = [], [], []
    for l in lines[2:2 + n]:
        p = l.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
        labels.append(p[4] if len(p) > 4 else p[0])
    return syms, np.array(xyz), labels, cell, head


def dist(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def angle(a, b, c):
    u, v = np.asarray(a) - np.asarray(b), np.asarray(c) - np.asarray(b)
    return float(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def dihedral(a, b, c, d):
    b0, b1, b2 = np.asarray(a) - np.asarray(b), np.asarray(c) - np.asarray(b), np.asarray(d) - np.asarray(c)
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1n) * b1n
    w = b2 - (b2 @ b1n) * b1n
    x = v @ w
    y = np.cross(b1n, v) @ w
    return float(np.degrees(np.arctan2(y, x)))


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal superposition (both N x 3)."""
    P = np.asarray(P) - np.mean(P, axis=0)
    Q = np.asarray(Q) - np.mean(Q, axis=0)
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P2 = (R @ P.T).T
    return float(np.sqrt(np.mean(np.sum((P2 - Q) ** 2, axis=1))))


def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    p.write_text(json.dumps(obj, indent=1, default=float))
    print(f'wrote {p.relative_to(REPO)}')


def rel(p):
    p = Path(p)
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)
