"""Which start makes the B3LYP SCF of a zinc cation converge in fewer cycles?  (small test complex)

    OMP_NUM_THREADS=1 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/test_scf_start.py

One SCF cycle of the real models costs minutes, so a better start would pay.  System: the hand-built
geometry of [Zn(imidazole)(NH3)(OH)]+ of md/qm/test_small_chain.py, B3LYP/6-31G*, density fitting, grid
level 2, convergence 1e-7 Hartree.  Run on 5 Oct 2026 with Cartesian d functions (the basis was changed to
pure functions afterwards; the count of cycles is what matters here).  Starts compared: PySCF's guesses
minao (default), atom, huckel, vsap, and a PBE calculation converged first.
Result: md/qm/test_scf_start.json.  The entry `small_basis_start` of that file was added by hand from the
logs of a run of the driver with `guess_basis: sto-3g` (job timing) on the same system.
Outcome: nothing beats the default by enough to matter (11 cycles; PBE first saves 4 expensive cycles for 9
cheap ones but goes through a Coulomb-only code path that solves with the ill-conditioned metric).
"""
import json
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pyscf_calc as pc  # noqa: E402
from test_small_chain import ATOMS  # noqa: E402
from pyscf import dft  # noqa: E402

s = pc.Settings(verbose=0, max_memory_mb=3000, cart=True)
mol = pc.build_mol([a[0] for a in ATOMS], [a[1:] for a in ATOMS], 1, 1, s)
out = {}


def new(xc='b3lyp'):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = 2
    mf.conv_tol = 1e-7
    return mf.density_fit(auxbasis='def2-universal-jkfit')


def run(mf, dm0=None):
    es = []
    mf.callback = lambda envs: es.append(float(envs['e_tot']))
    t = time.time()
    mf.kernel(dm0=dm0)
    return dict(cycles=len(es), converged=bool(mf.converged), e=float(mf.e_tot), seconds=round(time.time() - t, 1))


for guess in ('minao', 'atom', 'huckel', 'vsap'):
    mf = new()
    mf.init_guess = guess
    out[guess] = run(mf)
    print(guess, out[guess], flush=True)
mf = new('pbe')
mf.conv_tol = 1e-5
r1 = run(mf)
mf.xc = 'b3lyp'
mf.conv_tol = 1e-7
out['pbe first, then b3lyp'] = dict(pbe=r1, b3lyp=run(mf, mf.make_rdm1()))
old = json.loads((HERE / 'test_scf_start.json').read_text()) if (HERE / 'test_scf_start.json').exists() else {}
if 'small_basis_start' in old:
    out['small_basis_start'] = old['small_basis_start']
(HERE / 'test_scf_start.json').write_text(json.dumps(out, indent=1) + '\n')
