"""
Shared plumbing for the property and energetics scripts of burst 1b: dft/nlo.py (finite-field mu, alpha, beta),
dft/dscf.py (vertical IP/EA by Delta-SCF, conceptual-DFT descriptors), dft/tddft.py (UV-Vis), dft/fsp.py
(functional check single points) and dft/interaction.py (counterpoise interaction and deformation energies).

Why a module of its own:
  * These scripts are launched directly, not through ops/vm/run_case.sh, so importing this module does what
    that script does, BEFORE numpy loads: it pins every BLAS/OpenMP pool to DFT_THREADS (or --threads), defaults
    the PySCF memory to 6000 MB, and puts PySCF's scratch on the NVMe volume (/mnt/scratch) in a per-process
    directory removed on exit.  Unpinned BLAS oversubscribes the machine by an order of magnitude (run_case.sh).
  * They all take the same geometry / level-of-theory switches, so a code test at a tiny basis and grid on the
    Pi and the production run on the VM go through one code path:
        --xyz FILE        geometry to use (default: results/<case>/opt.xyz, the burst-1a optimum)
        --geom-from CASE  take results/<CASE>/opt.xyz instead (e.g. gas-phase response at the PCM geometry)
        --basis B         override the case basis (tests: sto-3g); --grid N  override the DFT grid level
        --outdir DIR      write results there instead of results/<case>/ (tests write to a scratch dir)
        --tag T           suffix for the output file names
        --threads N       thread count (else DFT_THREADS; else common.py's default)
  * Every result file carries the same provenance block (level of theory, geometry path + sha256, PySCF
    version, threads, host, wall time), so a number in the paper can be traced to the exact run.

IMPORTANT: import this module BEFORE `common` (i.e. before numpy / PySCF load).
"""
import os, sys, atexit, shutil, tempfile, pathlib


def _pin_environment():
    t = None
    if '--threads' in sys.argv:
        t = sys.argv[sys.argv.index('--threads') + 1]
        os.environ['DFT_THREADS'] = t
    t = t or os.environ.get('DFT_THREADS')
    if t:
        for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                  'VECLIB_MAXIMUM_THREADS'):
            os.environ[v] = t
    os.environ.setdefault('DFT_MAX_MEMORY_MB', '6000')
    if 'numpy' in sys.modules:
        print('[props] WARNING: numpy was imported before props_common; thread pinning may not apply',
              file=sys.stderr)
    if 'PYSCF_TMPDIR' not in os.environ and os.path.isdir('/mnt/scratch') and os.access('/mnt/scratch', os.W_OK):
        tmp = tempfile.mkdtemp(prefix=f'{pathlib.Path(sys.argv[0]).stem}_', dir='/mnt/scratch')
        os.environ['PYSCF_TMPDIR'] = os.environ['TMPDIR'] = tmp
        atexit.register(shutil.rmtree, tmp, True)


_pin_environment()

import argparse, contextlib, datetime, hashlib, json, socket, time        # noqa: E402
import numpy as np                                                         # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *                                                       # noqa: E402,F401,F403
import common as _common                                                   # noqa: E402
import pyscf                                                               # noqa: E402

KCAL_PER_HA = 627.509474          # CODATA-2018 hartree in kcal/mol
KJ_PER_HA = 2625.499639           # CODATA-2018 hartree in kJ/mol
AU_DEBYE = 2.541746473            # 1 e a0 in Debye


def parser(doc, geometry=True):
    """argparse with the switches every property script shares (see the module docstring)."""
    ap = argparse.ArgumentParser(description=doc.strip().split('\n')[0])
    ap.add_argument('case')
    if geometry:
        ap.add_argument('--xyz', help='geometry file (default results/<case>/opt.xyz)')
        ap.add_argument('--geom-from', help='use results/<CASE>/opt.xyz of another case')
    ap.add_argument('--basis', help='override the basis (code tests: sto-3g)')
    ap.add_argument('--grid', type=int, help='override the DFT grid level (code tests: 1)')
    ap.add_argument('--outdir', help='output directory (default results/<case>/)')
    ap.add_argument('--tag', help='suffix for the output file names')
    ap.add_argument('--threads', help='threads (else DFT_THREADS)')
    ap.add_argument('--verbose', type=int, default=3, help='PySCF verbosity (3 = SCF summaries, 4 = cycles)')
    return ap


def sha256(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def setup(args, need_geometry=True):
    """Load the case, apply the command-line overrides and resolve geometry, output directory and file tag."""
    case = load_case(args.case)
    if args.basis:
        case['basis'] = args.basis
    xyz, tag = None, args.tag
    if need_geometry:
        if getattr(args, 'xyz', None):
            xyz = pathlib.Path(args.xyz).resolve()
            tag = tag or f'xyz-{xyz.stem}'
        elif getattr(args, 'geom_from', None):
            xyz = load_case(args.geom_from)['dir'] / 'opt.xyz'
            # load_case() of the other case re-exported its DFT_XC / DFT_DISP; restore this case's definition
            load_case(args.case)
            tag = tag or f'geom-{args.geom_from}'
        else:
            xyz = case['dir'] / 'opt.xyz'
        if not xyz.exists():
            raise SystemExit(f'[props] geometry {xyz} not found -- run the optimisation first (burst 1a), '
                             f'or pass --xyz / --geom-from')
        syms, _ = read_xyz(xyz); syms0, _ = read_xyz(ROOT / case['xyz'])
        if [s.capitalize() for s in syms] != [s.capitalize() for s in syms0]:
            raise SystemExit(f'[props] {xyz}: element sequence differs from the case input {case["xyz"]}; '
                             f'refusing to use it with case {case["name"]}')
    outdir = pathlib.Path(args.outdir).resolve() if args.outdir else case['dir']
    outdir.mkdir(parents=True, exist_ok=True)
    if args.outdir:
        # a shared output directory (tests) holds several cases: keep their files apart
        tag = case['name'] + (f'_{tag}' if tag else '')
    return case, xyz, outdir, tag


def out_name(stem, tag, ext='.json'):
    return f'{stem}_{tag}{ext}' if tag else f'{stem}{ext}'


def grid_level(args):
    return args.grid if getattr(args, 'grid', None) else _common.GRID_LEVEL


@contextlib.contextmanager
def level(xc=None, disp=None):
    """Temporarily change the functional / dispersion that common.make_mf() reads from the environment.

    make_mf() does not receive the case; load_case() exports DFT_XC and DFT_DISP instead.  Scripts that need a
    different functional for one step (TD-DFT with CAM-B3LYP, the functional check) go through here so the
    case's own definition is back in place afterwards."""
    old = {k: os.environ.get(k) for k in ('DFT_XC', 'DFT_DISP')}
    try:
        if xc is not None:
            os.environ['DFT_XC'] = xc
        if disp is not None:
            os.environ['DFT_DISP'] = disp or 'none'
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def scf_for(case, mol, args=None, xc=None, disp=None, solvent='case', conv_tol=1e-9, chkfile=None):
    """make_mf() with this run's grid override; solvent='case' means the case's own environment."""
    solv = case.get('solvent') if solvent == 'case' else solvent
    with level(xc, disp):
        mf = make_mf(mol, conv_tol=conv_tol, grid_level=grid_level(args) if args else None, chkfile=chkfile,
                     solvent=solv)
    return mf


def environment_label(case, solvent='case'):
    solv = case.get('solvent') if solvent == 'case' else solvent
    if not solv:
        return 'gas phase'
    return (f"IEF-PCM ({solv}, eps = {SOLVENT_EPS[solv.lower()]}, Lebedev order "
            f"{os.environ.get('PCM_LEBEDEV_ORDER', '29')})")


def provenance(case, args, xyz=None, t0=None, extra=None):
    """The block every result file carries: enough to reproduce the number and to catch a mismatched run."""
    p = {'script': pathlib.Path(sys.argv[0]).name, 'argv': sys.argv[1:], 'case': case['name'],
         'case_file': str(pathlib.Path(case['case_file']).relative_to(ROOT)) if str(case['case_file']).startswith(str(ROOT))
         else str(case['case_file']),
         'charge': case['charge'], 'spin': case['spin'],
         'basis': basis_description(case), 'grid_level': grid_level(args) if args else _common.GRID_LEVEL,
         'environment': environment_label(case), 'case_xc': case.get('xc') or XC, 'case_disp': case.get('disp') or 'none',
         'density_fitting': True, 'pyscf_version': pyscf.__version__, 'threads': _common.NTHREADS,
         'host': socket.gethostname(), 'finished': datetime.datetime.now().isoformat(timespec='seconds')}
    if xyz is not None:
        p['geometry'] = str(xyz.relative_to(ROOT)) if str(xyz).startswith(str(ROOT)) else str(xyz)
        p['geometry_sha256'] = sha256(xyz)
        if case.get('proton'):
            # which side holds the proton in THIS geometry: the gas-phase M1 optimum is the neutral pair
            # (proton on O5), so its numbers describe sulfanilamide...H3PO4, not the salt
            _, x = read_xyz(xyz)
            pr = case['proton']
            dn = float(np.linalg.norm(x[pr['H']] - x[pr['donor_N']]))
            do = float(np.linalg.norm(x[pr['H']] - x[pr['acceptor_O']]))
            p['proton'] = {'d_H_N_A': dn, 'd_H_O_A': do, 'on': 'N (ionic pair)' if dn < do else 'O (neutral pair)'}
    if t0 is not None:
        p['wall_min'] = round((time.time() - t0) / 60, 2)
    if extra:
        p.update(extra)
    return p


def converge(mf, dm0=None, what='SCF'):
    """Run an SCF; if it does not converge, retry once with damping + level shift (ions and field-perturbed
    systems occasionally oscillate), then give up loudly -- a non-converged energy must never reach a table."""
    mf.kernel(dm0=dm0)
    if not mf.converged:
        print(f'[props] {what}: not converged in {mf.max_cycle} cycles; retrying with damping 0.3 / level shift 0.2',
              flush=True)
        dm = mf.make_rdm1()
        mf.damp, mf.level_shift, mf.max_cycle = 0.3, 0.2, 300
        mf.kernel(dm0=dm)
        mf.damp, mf.level_shift = 0, 0
        if mf.converged:                     # polish without the shift so orbital energies are unbiased
            mf.kernel(dm0=mf.make_rdm1())
    if not mf.converged:
        raise RuntimeError(f'{what} did not converge')
    return mf


def speed_up_pcm():
    """Faster, numerically identical PCM surface-potential contractions (PySCF 2.14 solvent/pcm.py).

    Why: PCM._get_v / _get_vmat recompute the (AO pair x surface point) integrals on every call and contract them
    with a single-threaded numpy.einsum, once per density.  Every PCM SCF cycle calls both, and a TD-DFT Davidson
    iteration calls both for all ~40 trial densities: for sulfanilamide (306 AOs, 2333 points) that was 16 s of
    integrals + 99 s of einsum per call on the Pi, against 14 s as one BLAS product (max difference 5e-16).
    Here the integrals are kept in memory when they fit PCM_CACHE_GB (default 8 GB; M1 needs 5.5 GB: 445 AOs, 3473
    surface points at Lebedev 29) and the
    contractions are matrix products.  Only the order of floating-point summation changes.  PCM_FAST=0 restores
    PySCF's own code.  Energies only: the PCM gradient code (not used by these scripts) is untouched.
    Checked on the VM (27 Sep 2026, sulfanilamide, CAM-B3LYP/6-31+G(d), IEF-PCM water, 5 TD-DFT roots): SCF energy
    identical to 2e-12 Ha, excitation energies identical to 1e-8 eV; SCF 46 s vs 155 s, TD-DFT 1653 s vs 2085 s."""
    if os.environ.get('PCM_FAST', '1') == '0':
        return
    from pyscf.solvent import pcm as _pcm
    from pyscf import df as _df, gto as _gto, lib as _lib
    cache = {}
    limit = float(os.environ.get('PCM_CACHE_GB', '8')) * 1e9

    def blocks(self):
        mol, s = self.mol, self.surface
        g, ex = s['grid_coords'], s['charge_exp']
        nao, ng = mol.nao, g.shape[0]
        key = (nao, ng, hashlib.sha1(g.tobytes() + ex.tobytes() + mol._env.tobytes() + mol._bas.tobytes()
                                     + mol._atm.tobytes()).hexdigest())
        if cache.get('key') == key:
            yield 0, ng, cache['v']
            return
        intor = mol._add_suffix('int3c2e')
        opt = _gto.moleintor.make_cintopt(mol._atm, mol._bas, mol._env, intor)

        def ints(p0, p1):
            fake = _gto.fakemol_for_charges(g[p0:p1], expnt=ex[p0:p1] ** 2)
            fake.cart = mol.cart
            # aux_e2 returns a Fortran-ordered (nao, nao, ngrid) array; a C-order reshape would copy it, doubling the
            # peak (11 GB instead of 5.5 GB for M1 at Lebedev 29).  The Fortran-order reshape is a view whose row
            # index is i + nao*j instead of i*nao + j, which changes nothing because (ij|g) = (ji|g).
            # [review 28 Sep 2026]
            return _df.incore.aux_e2(mol, fake, intor=intor, aosym='s1', cintopt=opt).reshape(nao * nao, p1 - p0,
                                                                                              order='F')
        if nao * nao * ng * 8 <= limit:
            cache.clear()
            cache.update(key=key, v=ints(0, ng))
            yield 0, ng, cache['v']
            return
        blk = int(max((self.max_memory - _lib.current_memory()[0]) * .9e6 / 8 / nao ** 2, 400))
        for p0, p1 in _lib.prange(0, ng, blk):
            yield p0, p1, ints(p0, p1)

    def _get_v(self, dms):
        nao = dms.shape[-1]
        d = dms.reshape(-1, nao * nao)
        out = np.empty((d.shape[0], self.surface['grid_coords'].shape[0]))
        for p0, p1, v in blocks(self):
            out[:, p0:p1] = d @ v
        return out

    def _get_vmat(self, q):
        nao = self.mol.nao
        ng = self.surface['grid_coords'].shape[0]
        q = q.reshape(-1, ng)
        out = np.zeros((q.shape[0], nao * nao))
        for p0, p1, v in blocks(self):
            out -= q[:, p0:p1] @ v.T
        return out.reshape(-1, nao, nao)

    _pcm.PCM._get_v = _get_v
    _pcm.PCM._get_vmat = _get_vmat


speed_up_pcm()


def md_table(rows, header):
    out = '| ' + ' | '.join(header) + ' |\n|' + '---|' * len(header) + '\n'
    for r in rows:
        out += '| ' + ' | '.join(str(v) for v in r) + ' |\n'
    return out
