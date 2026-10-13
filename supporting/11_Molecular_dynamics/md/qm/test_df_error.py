"""What density fitting and the grid cost in accuracy, measured on the small test complex.

    OMP_NUM_THREADS=3 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/test_df_error.py

System: [Zn(imidazole)(NH3)(OH)]+ at the hand-built geometry of md/qm/test_small_chain.py (16 atoms),
B3LYP/6-31G* with pure functions, SCF to 1e-10 Hartree.

 D1  density fitting (def2-universal-JKFIT) against exact four-centre integrals, grid level 2:
     total energy, analytic gradient, electrostatic potential at the Merz-Kollman points, and the charges
     of an unrestrained fit to that potential.
 D2  the gradient without and with the derivative of the grid weights (PySCF: grid_response), grid level 2
     and 3: sum of the forces over the atoms (zero for an exact gradient) and the largest change.
Result: md/qm/test_df_error.json.
"""
import json
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gauout  # noqa: E402
import pyscf_calc as pc  # noqa: E402
import run_qm  # noqa: E402
from test_small_chain import ATOMS  # noqa: E402
from pyscf import dft  # noqa: E402

sym, xyz = [a[0] for a in ATOMS], np.array([a[1:] for a in ATOMS])
s = pc.Settings(verbose=0, max_memory_mb=4000)
mol = pc.build_mol(sym, xyz, 1, 1, s)
pts, _, _ = gauout.mk_points(sym, xyz, {'Zn': 1.456}, density=6.0)
pb = pts / gauout.BOHR
out = dict(system='[Zn(imidazole)(NH3)(OH)]+, 16 atoms, B3LYP/6-31G* (5D, 7F)', nao=int(mol.nao), n_points=int(len(pts)),
           versions=pc.versions())


def scf(level, density_fit):
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.grids.level = level
    mf.conv_tol = 1e-10
    if density_fit:
        mf = mf.density_fit(auxbasis='def2-universal-jkfit')
    t = time.time()
    mf.kernel()
    return mf, time.time() - t


r = {}
for name, dfit in (('df', True), ('exact', False)):
    mf, t = scf(2, dfit)
    g = mf.nuc_grad_method().kernel()
    v = pc.esp_at_points(mol, mf.make_rdm1(), pb)
    q, rrms = run_qm.esp_fit_charges(mol.atom_coords(), pb, v, 1)
    r[name] = dict(e=float(mf.e_tot), g=g, v=v, q=q, rrms=rrms, seconds=round(t, 1))
    if dfit:
        mf_df = mf
out['D1_density_fitting'] = dict(
    energy_df_minus_exact_hartree=r['df']['e'] - r['exact']['e'],
    gradient_max_abs_difference=float(np.abs(r['df']['g'] - r['exact']['g']).max()),
    gradient_largest_element=float(np.abs(r['exact']['g']).max()),
    potential_max_abs_difference=float(np.abs(r['df']['v'] - r['exact']['v']).max()),
    potential_rms_difference=float(np.sqrt(((r['df']['v'] - r['exact']['v']) ** 2).mean())),
    potential_rms=float(np.sqrt((r['exact']['v'] ** 2).mean())),
    fitted_charges_max_abs_difference=float(np.abs(r['df']['q'] - r['exact']['q']).max()),
    zinc_charge_df=float(r['df']['q'][sym.index('Zn')]), zinc_charge_exact=float(r['exact']['q'][sym.index('Zn')]),
    scf_seconds_df=r['df']['seconds'], scf_seconds_exact=r['exact']['seconds'], units='atomic units; charges in e')
(HERE / 'test_df_error.json').write_text(json.dumps(out, indent=1) + '\n')
rows = []
for level in (2, 3):
    mf = mf_df if level == 2 else scf(3, True)[0]
    res = {}
    for gr in (False, True):
        gobj = mf.nuc_grad_method()
        gobj.grid_response = gr
        t = time.time()
        g = gobj.kernel()
        res[gr] = (g, time.time() - t)
    rows.append(dict(grid_level=level, grid_points=int(mf.grids.weights.size),
                     sum_of_forces_without=[float(x) for x in res[False][0].sum(axis=0)],
                     sum_of_forces_with=[float(x) for x in res[True][0].sum(axis=0)],
                     largest_change_of_a_gradient_element=float(np.abs(res[True][0] - res[False][0]).max()),
                     seconds_without=round(res[False][1], 1), seconds_with=round(res[True][1], 1)))
out['D2_grid_response'] = dict(rows=rows, units='Hartree/Bohr')
(HERE / 'test_df_error.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
