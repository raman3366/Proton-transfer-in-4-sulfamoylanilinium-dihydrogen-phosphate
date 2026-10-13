"""Adapters (b) and (c): files in Gaussian's formats that MCPB.py reads back, and the Merz-Kollman points.

Pure NumPy, so that it runs in the quantum-chemistry environment and in the MD environment alike.

(b) write_fchk          formatted checkpoint file with the geometry and the Cartesian force constants.
                        MCPB.py step 2 (pymsmt/mol/gauio.py: get_crds_from_fchk, get_matrix_from_fchk) reads
                          'Current cartesian coordinates'   3N numbers, Bohr, atom by atom (x1 y1 z1 x2 ...)
                          'Cartesian Force Constants'       3N(3N+1)/2 numbers, Hartree/Bohr^2, the LOWER
                                                            triangle row by row: H11, H21, H22, H31, H32, H33, ...
                        five numbers per line (Fortran 5E16.8).
(c) write_mk_log        the part of a Gaussian log that MCPB.py step 3 (get_esp_from_gau) and espgen read:
                          'Atomic Center  i is at  x y z'   ANGSTROM, F10.6 in columns 33-62
                          'ESP Fit Center i is at  x y z'   ANGSTROM, same columns
                          'i Atom  V' and 'i Fit  V'        potential in atomic units (Hartree per e), F14.6
                        The readers divide the coordinates by 0.529177249 to get Bohr for the .esp file.
    mk_points           fit points of the Merz-Singh-Kollman scheme as Gaussian's Pop=MK makes them: for each
                        atom and each layer (1.4, 1.6, 1.8 and 2.0 times the radius) points on a sphere,
                        `density` points per square Angstrom (IOp(6/42); the tutorial input has 6), kept when
                        they lie outside the scaled sphere of every other atom.  Radii: Gaussian's MK set
                        (H 1.20, C 1.50, N 1.50, O 1.40, S 1.75 Angstrom) and for the metal the value MCPB.py
                        writes after the geometry (ReadRadii; Zn 1.456).  The points on a sphere are placed
                        by the latitude-circle rule of Connolly's surface program, which the MK code uses.
                        Gaussian itself is not here to compare with: the construction is the published
                        one, the individual points need not coincide with Gaussian's.

Units are explicit in every signature: `_ang` Angstrom, `_bohr` Bohr.
"""
import math

import numpy as np

BOHR = 0.52917721092          # Angstrom per Bohr (PySCF 2.14, CODATA 2010)
MK_RADII = {'H': 1.20, 'C': 1.50, 'N': 1.50, 'O': 1.40, 'F': 1.35, 'P': 1.80, 'S': 1.75, 'Cl': 1.70}
MK_LAYERS = (1.4, 1.6, 1.8, 2.0)
Z_OF = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Zn': 30}
MASS = {1: 1.00782504, 6: 12.0, 7: 14.0030740, 8: 15.9949146, 16: 31.9720718, 30: 63.9291454}


# ---------------------------------------------------------------------------------------------------- (b) fchk
def _ints(values):
    out = []
    for i in range(0, len(values), 6):
        out.append(''.join('%12d' % v for v in values[i:i + 6]))
    return out


def _reals(values):
    out = []
    for i in range(0, len(values), 5):
        out.append(''.join('%16.8E' % v for v in values[i:i + 5]))
    return out


def lower_triangle(h):
    """H11, H21, H22, H31, ... of a symmetric matrix."""
    n = h.shape[0]
    return h[np.tril_indices(n)]


def write_fchk(path, atomic_numbers, coords_bohr, hessian=None, *, charge, multiplicity, energy=None,
               gradient=None, title='', job='Freq', method='RB3LYP', basis='6-31G(d)', nbasis=None,
               symmetrise_tolerance=1e-6):
    """coords_bohr (N,3) in Bohr; hessian (3N,3N) in Hartree/Bohr^2; gradient (N,3) in Hartree/Bohr."""
    z = [int(x) for x in atomic_numbers]
    n = len(z)
    xyz = np.asarray(coords_bohr, dtype=float).reshape(n, 3)
    lines = [title[:72], '%-10s%-30s%30s' % (job, method, basis),
             '%-40s   I     %12d' % ('Number of atoms', n),
             '%-40s   I     %12d' % ('Charge', charge),
             '%-40s   I     %12d' % ('Multiplicity', multiplicity),
             '%-40s   I     %12d' % ('Number of electrons', sum(z) - charge)]
    if nbasis is not None:
        lines.append('%-40s   I     %12d' % ('Number of basis functions', nbasis))
    lines.append('%-40s   I   N=%12d' % ('Atomic numbers', n))
    lines += _ints(z)
    lines.append('%-40s   R   N=%12d' % ('Nuclear charges', n))
    lines += _reals([float(v) for v in z])
    lines.append('%-40s   R   N=%12d' % ('Current cartesian coordinates', 3 * n))
    lines += _reals(xyz.ravel())
    lines.append('%-40s   R   N=%12d' % ('Real atomic weights', n))
    lines += _reals([MASS.get(v, 0.0) for v in z])
    if energy is not None:
        lines.append('%-40s   R     %22.15E' % ('Total Energy', energy))
    if gradient is not None:
        lines.append('%-40s   R   N=%12d' % ('Cartesian Gradient', 3 * n))
        lines += _reals(np.asarray(gradient, dtype=float).ravel())
    if hessian is not None:
        h = np.asarray(hessian, dtype=float)
        if h.shape != (3 * n, 3 * n):
            raise ValueError(f'Hessian has shape {h.shape}, expected {(3 * n, 3 * n)}')
        asym = float(np.abs(h - h.T).max())
        if asym > symmetrise_tolerance:
            raise ValueError(f'Hessian is not symmetric: largest |H - H^T| = {asym:.3e} Hartree/Bohr^2')
        h = 0.5 * (h + h.T)
        tri = lower_triangle(h)
        lines.append('%-40s   R   N=%12d' % ('Cartesian Force Constants', tri.size))
        lines += _reals(tri)
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def read_fchk_array(path, label):
    """Independent reader for the round-trip test (not the reader of MCPB.py)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    for i, l in enumerate(lines):
        if l.startswith(label):
            n = int(l.split('N=')[1])
            vals = []
            k = i + 1
            while len(vals) < n:
                vals += [float(lines[k][j:j + 16]) for j in range(0, len(lines[k]), 16)]
                k += 1
            return np.array(vals)
    raise KeyError(label)


# ------------------------------------------------------------------------------------------------ MK points
def unit_sphere(n):
    """About n points on the unit sphere on latitude circles (Connolly's 'genun')."""
    nequat = int(math.sqrt(math.pi * n))
    nvert = max(nequat // 2, 1)
    pts = []
    for i in range(nvert + 1):
        fi = math.pi * i / nvert
        z, xy = math.cos(fi), math.sin(fi)
        nhor = max(int(nequat * xy), 1)
        for j in range(nhor):
            fj = 2.0 * math.pi * j / nhor
            pts.append((math.cos(fj) * xy, math.sin(fj) * xy, z))
    return np.array(pts)


def mk_radii(symbols, extra=None):
    r = dict(MK_RADII)
    r.update(extra or {})
    missing = sorted({s for s in symbols if s not in r})
    if missing:
        raise KeyError(f'no Merz-Kollman radius for {missing}')
    return np.array([r[s] for s in symbols])


def mk_points(symbols, coords_ang, radii_extra=None, density=6.0, layers=MK_LAYERS):
    """Fit points in ANGSTROM, shape (M, 3), and for each point the layer and the atom it was built on."""
    xyz = np.asarray(coords_ang, dtype=float)
    rad = mk_radii(symbols, radii_extra)
    pts, lay, own = [], [], []
    for f in layers:
        rs = rad * f
        for i in range(len(symbols)):
            n = int(density * 4.0 * math.pi * rs[i] ** 2)
            p = xyz[i] + rs[i] * unit_sphere(n)
            d = np.linalg.norm(p[:, None, :] - xyz[None, :, :], axis=2)
            d[:, i] = np.inf
            keep = (d >= rs[None, :]).all(axis=1)
            pts.append(p[keep])
            lay += [f] * int(keep.sum())
            own += [i] * int(keep.sum())
    return np.vstack(pts), np.array(lay), np.array(own)


# ---------------------------------------------------------------------------------------------- (c) MK log
def write_mk_log(path, symbols, coords_ang, points_ang, esp_atoms, esp_points, *, header_lines=()):
    """coords_ang (N,3) and points_ang (M,3) in ANGSTROM; potentials in atomic units."""
    xyz = np.asarray(coords_ang, dtype=float)
    pts = np.asarray(points_ang, dtype=float)
    va = np.asarray(esp_atoms, dtype=float)
    vp = np.asarray(esp_points, dtype=float)
    if len(va) != len(xyz) or len(vp) != len(pts):
        raise ValueError('number of potentials does not match the number of centres')
    if len(pts) > 99999:
        # the .esp file that espgen and MCPB.py write has the header (2I5); resp cannot read more
        raise ValueError(f'{len(pts)} fit points: the .esp format of resp holds at most 99,999')
    if np.abs(xyz).max() >= 100.0 or np.abs(pts).max() >= 100.0:
        # F10.6 holds -99.999999; the readers take fixed columns
        raise ValueError('a coordinate does not fit the F10.6 field of the Gaussian format')
    star = ' ' + '*' * 70
    out = [' Log file in the format of Gaussian, written by md/qm/gauout.py from a PySCF calculation.',
           ' It holds the part that espgen and MCPB.py read: centres in Angstrom, potentials in atomic units.']
    out += [' ' + h for h in header_lines]
    out += ['', star, '', '            Electrostatic Properties Using The SCF Density', '', star, '']

    def idx(i, w):
        s = '%d' % i
        return s.rjust(w) if len(s) <= w else '*' * w

    for i, r in enumerate(xyz, 1):
        out.append('       Atomic Center %s is at %10.6f%10.6f%10.6f' % (idx(i, 4), r[0], r[1], r[2]))
    for i, r in enumerate(pts, 1):
        out.append('      ESP Fit Center %s is at %10.6f%10.6f%10.6f' % (idx(i, 4), r[0], r[1], r[2]))
    out += [' %5d points will be used for fitting atomic charges' % len(pts), '',
            '              Electrostatic Properties (Atomic Units)', '',
            ' -----------------------------------------------------------------',
            '    Center     Electric         -------- Electric Field --------',
            '               Potential          X             Y             Z',
            ' -----------------------------------------------------------------']
    k = 0
    for v in va:
        k += 1
        out.append('%s Atom%14.6f' % (idx(k, 5), v))
    for v in vp:
        k += 1
        out.append('%s Fit %14.6f' % (idx(k, 5), v))
    out += [' -----------------------------------------------------------------', '',
            ' Normal termination of the adapter (no Gaussian run).']
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')


def read_esp_file(path):
    """The .esp file of espgen / MCPB.py: Bohr and atomic units.  Returns atoms, points, potentials."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    nat, npt = int(lines[0][:5]), int(lines[0][5:10])
    atoms = np.array([[float(x) for x in l.split()] for l in lines[1:1 + nat]])
    rest = np.array([[float(x) for x in l.split()] for l in lines[1 + nat:1 + nat + npt]])
    return atoms, rest[:, 1:], rest[:, 0]
