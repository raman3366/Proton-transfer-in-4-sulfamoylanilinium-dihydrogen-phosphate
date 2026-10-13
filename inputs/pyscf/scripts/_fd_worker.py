"""
Worker for the parallel finite-difference Hessian (see dft/frequencies.py).
Computes the analytic gradient and the dipole moment at displaced geometries for a slice of Cartesian
coordinates and writes the Hessian rows and dipole derivatives (keys 'dip<col>') to an .npz file.

Usage: python dft/_fd_worker.py <case> <col_start> <col_end> <out.npz> <step_bohr>
Run one process per worker so each gets a clean OpenMP thread pool.
"""
import sys, os, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *

case = load_case(sys.argv[1]); d = case['dir']
c0, c1, out, h = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], float(sys.argv[5])

mol = mol_for_case(case, d / 'opt.xyz', verbose=0)
x0 = mol.atom_coords().copy()
natm = mol.natm
rows = {}
for col in range(c0, c1):
    a, k = divmod(col, 3)
    gs, dips = [], []
    for sgn in (+1, -1):
        x = x0.copy(); x[a, k] += sgn * h
        m2 = mol.copy(); m2.set_geom_(x, unit='Bohr')
        mf2 = make_mf(m2, conv_tol=1e-10, solvent=case.get('solvent'))
        mf2.verbose = 0
        # start from the converged density of the undisplaced geometry: a 0.005 Bohr shift is tiny,
        # so the SCF needs only a handful of cycles
        run_scf(mf2, init_chk=d / 'final.chk')
        gs.append(mf2.nuc_grad_method().kernel())
        # dipole (a.u.) at the displaced geometry: its central difference gives the dipole derivative
        # that turns the harmonic frequencies into IR intensities at no extra SCF cost
        dips.append(np.asarray(mf2.dip_moment(unit='AU', verbose=0)))
    rows[str(col)] = (gs[0] - gs[1]) / (2 * h)
    rows[f'dip{col}'] = (dips[0] - dips[1]) / (2 * h)
    print(f'[fd-worker {c0}-{c1}] column {col + 1}/{3 * natm} done', flush=True)
np.savez_compressed(out, **rows)
