"""Seminario's projection of a Cartesian Hessian on bonds and angles, written independently of MCPB.py.

J. M. Seminario, Int. J. Quantum Chem. 1996, 60, 1271 (doi:10.1002/(SICI)1097-461X(1996)60:7<1271::AID-QUA8>3.0.CO;2-W).

Used (1) to cross-check what MCPB.py reads from the fchk file (same formula, own code, the Hessian taken from
memory and not from the file) and (2) for terms MCPB.py does not derive, such as the angles at sulfur.

Conventions: Hessian in Hartree/Bohr^2, coordinates in Bohr, atom indices from 0.  Results in the AMBER
form E = k (q - q0)^2: kcal/mol/A^2 for bonds, kcal/mol/rad^2 for angles (hence the factor 1/2).
The interatomic block taken is [a, b] (rows of a, columns of b), as MCPB.py takes it with bondfc_avg 0 and
anglefc_avg 0; the block is not symmetric and its eigenvalues can be complex, the real part is used.

The angle constant of this method is not exact even for a model with one harmonic angle: the projection mixes
the two arms through the cosine of the angle.  It is an estimate, as in every use of the method.
"""
import numpy as np

BOHR = 0.52917721092
HARTREE_KCAL = 627.5094740631


def _projected(h, a, b, direction):
    blk = -np.asarray(h)[3 * a:3 * a + 3, 3 * b:3 * b + 3]
    w, v = np.linalg.eig(blk)
    return float(np.real(sum(w[i] * abs(np.dot(direction, v[:, i])) for i in range(3))))


def bond(h, xyz_bohr, a, b):
    """(length in Angstrom, k in kcal/mol/A^2)."""
    x = np.asarray(xyz_bohr, dtype=float)
    d = x[b] - x[a]
    r = np.linalg.norm(d)
    k = _projected(h, a, b, d / r)
    return r * BOHR, 0.5 * k * HARTREE_KCAL / BOHR ** 2


def angle(h, xyz_bohr, a, b, c):
    """a-b-c with b in the middle: (angle in degrees, k in kcal/mol/rad^2)."""
    x = np.asarray(xyz_bohr, dtype=float)
    uab = x[b] - x[a]
    ucb = x[b] - x[c]
    rab, rcb = np.linalg.norm(uab), np.linalg.norm(ucb)
    uab, ucb = uab / rab, ucb / rcb
    un = np.cross(ucb, uab)
    un /= np.linalg.norm(un)
    upa = np.cross(un, uab)
    upc = np.cross(ucb, un)
    kab = _projected(h, a, b, upa)
    kcb = _projected(h, c, b, upc)
    k = 1.0 / (1.0 / (rab * rab * kab) + 1.0 / (rcb * rcb * kcb))
    theta = np.degrees(np.arccos(np.clip(np.dot(uab, ucb), -1, 1)))
    return float(theta), 0.5 * k * HARTREE_KCAL


def frequencies(h, masses_amu, xyz_bohr=None):
    """Harmonic wavenumbers in cm^-1 (negative = imaginary) of a full Cartesian Hessian in Hartree/Bohr^2.

    With coordinates the six translations and rotations are projected out; otherwise the 3N raw values.
    """
    amu = 1822.888486209
    m = np.repeat(np.asarray(masses_amu, dtype=float) * amu, 3)
    hm = np.asarray(h) / np.sqrt(np.outer(m, m))
    if xyz_bohr is not None:
        x = np.asarray(xyz_bohr, dtype=float)
        n = len(x)
        mm = np.asarray(masses_amu, dtype=float)
        com = (x * mm[:, None]).sum(0) / mm.sum()
        y = x - com
        basis = []
        for k in range(3):
            v = np.zeros((n, 3))
            v[:, k] = np.sqrt(mm)
            basis.append(v.ravel())
        for k in range(3):
            e = np.zeros(3)
            e[k] = 1.0
            v = np.cross(e, y) * np.sqrt(mm)[:, None]
            basis.append(v.ravel())
        q, _ = np.linalg.qr(np.array(basis).T)
        p = np.eye(3 * n) - q @ q.T
        hm = p @ hm @ p
    w = np.linalg.eigvalsh(0.5 * (hm + hm.T))
    au_to_cm = 219474.6313632
    return np.sign(w) * np.sqrt(np.abs(w)) * au_to_cm
