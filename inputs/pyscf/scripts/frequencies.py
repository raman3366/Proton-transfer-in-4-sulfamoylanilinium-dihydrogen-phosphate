"""
Step 2 -- harmonic frequencies at the optimised geometry (confirms a true minimum: no imaginary modes).
Usage: python dft/frequencies.py <case>
Writes results/<case>/freq.json  (all frequencies in cm-1, imaginary ones negative; thermochemistry)

Hessian strategy
  * Analytic Hessians (DF-RKS, with or without PCM) build all derivative integrals at once and need far more
    memory than the basis-set size suggests: 3ATB at 400 basis functions already exceeded 8 GB.  So anything
    past ~300 basis functions goes to finite differences of analytic gradients.
  * Those 6N gradient evaluations are independent, and PySCF's own threading is strongly sublinear past about
    four cores.  They are therefore farmed out to FD_WORKERS separate processes of FD_WORKER_THREADS threads
    each (default 6 x 2), which is 2.5-3.5x faster than one big serial loop on the same machine.
  * A machine-wide lock (tools/freq.lock) keeps two molecules from running frequencies at the same time, so
    whichever case gets there first has all the cores to itself.
Environment: FREQ_ANALYTIC=1 (force the analytic Hessian), FREQ_FD=1 (force finite differences),
FD_WORKERS (upper bound; the pool is sized from idle cores), FD_WORKER_THREADS, FD_RESERVE_CORES,
FREQ_NOLOCK=1 (skip the machine-wide lock), FD_STEP (displacement in Bohr, default 0.005).
"""
import sys, os, time, fcntl, subprocess, tempfile, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *
from pyscf import lib
from pyscf.hessian import thermo

case = load_case(sys.argv[1]); d = case['dir']
WTHREADS = int(os.environ.get('FD_WORKER_THREADS', '2'))
# Displacement for the finite-difference Hessian.  0.005 Bohr is the default and is what every
# completed case used.  A larger step (FD_STEP=0.01) averages over more of the potential and is the
# remedy when a near-zero mode comes back marginally imaginary: at 0.005 the second difference of the
# gradients is the difference of two nearly equal numbers, so a -2.5 cm-1 mode is inside the noise.
STEP = float(os.environ.get('FD_STEP', '0.005'))
RESERVE = int(os.environ.get('FD_RESERVE_CORES', '2'))       # leave this many cores for other users of the box


def size_worker_pool():
    """How many FD workers fit right now.

    The machine-wide lock below only serialises frequency runs against each other -- other cases may still
    be optimising.  Sizing the pool from a fixed FD_WORKERS would then oversubscribe the machine (12 FD
    threads on top of a 6-thread optimisation on 16 cores).  So take the cores that are actually idle,
    measured once the lock is held, with FD_WORKERS as an upper bound.

    FD_MIN_WORKERS is a floor: a frequency run that starts while the box is full would otherwise be pinned
    to a single worker for its whole duration, which is far worse than mild oversubscription.
    """
    ncpu = os.cpu_count() or 4
    busy = os.getloadavg()[0]
    free = max(1.0, ncpu - busy - RESERVE)
    n = max(int(os.environ.get('FD_MIN_WORKERS', '2')), int(free // WTHREADS))
    cap = os.environ.get('FD_WORKERS')          # treated as an upper bound, not an exact count
    return min(n, int(cap)) if cap else n

def use_fd_planned():
    """Will this run use finite differences?  Decided before the molecule is built so the thread count can
    be set correctly: the nao > 300 rule is irrelevant when FREQ_ANALYTIC forces the analytic path."""
    if os.environ.get('FREQ_ANALYTIC') == '1':
        return False
    return bool(case.get('solvent')) or os.environ.get('FREQ_FD') == '1'


lock_path = ROOT / 'tools' / 'freq.lock'
lock_path.parent.mkdir(parents=True, exist_ok=True)
lock_fh = open(lock_path, 'w')
if os.environ.get('FREQ_NOLOCK') != '1':
    t_wait = time.time()
    try:
        fcntl.flock(lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        stage_note(case, 'frequencies: waiting for the machine-wide frequency lock (another case is running)')
        print('[freq] another case holds the frequency lock; waiting for the whole machine ...', flush=True)
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
    if time.time() - t_wait > 5:
        print(f'[freq] lock acquired after {(time.time()-t_wait)/60:.1f} min', flush=True)
    stage_note(case, 'frequencies: holding the machine-wide lock')

WORKERS = size_worker_pool()
print(f'[freq] {WORKERS} finite-difference workers x {WTHREADS} threads '
      f'(load {os.getloadavg()[0]:.1f} on {os.cpu_count()} cores, reserving {RESERVE})', flush=True)

# This stage holds the machine-wide lock, so it owns the box: an analytic Hessian is a single process and
# must be given every physical core, not the 4 threads the rest of the pipeline runs at.  (Finite differences
# get their parallelism from the worker pool instead, so they keep the per-worker thread count.)
if not use_fd_planned():
    # Size from cores that are actually idle, exactly as the finite-difference pool does.  The machine-wide
    # lock only serialises frequency runs against each other; other cases may still be optimising, so a fixed
    # "take every core" would oversubscribe the box.
    def analytic_threads():
        if os.environ.get('FREQ_THREADS'):
            return int(os.environ['FREQ_THREADS'])
        phys = max(1, (os.cpu_count() or 4) // 2)
        busy = os.getloadavg()[0] / 2.0          # load counts logical threads; halve for physical cores
        return max(4, int(phys - busy - RESERVE))
    nthreads = analytic_threads()
    lib.num_threads(nthreads)
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        os.environ[v] = str(nthreads)
    print(f'[freq] analytic Hessian holds the machine lock -> {nthreads} threads '
          f'(load {os.getloadavg()[0]:.1f} on {os.cpu_count()} logical)', flush=True)

mol = mol_for_case(case, d / 'opt.xyz')
mf = make_mf(mol, conv_tol=1e-10, chkfile=d / 'final.chk', solvent=case.get('solvent'))
run_scf(mf, init_chk=d / 'opt.chk')
print(f"[freq] SCF E = {mf.e_tot:.8f} Ha ({mf._wall:.0f} s); {mol.nao} basis functions, {mol.natm} atoms", flush=True)
t0 = time.time()
method = 'analytic'
# FREQ_ANALYTIC=1 forces the analytic Hessian: it has no finite-difference truncation error, so it is the
# better choice whenever the machine has the memory for it (roughly 8.4 GB gas / 11 GB PCM at 400 basis
# functions on the 3ATB run, scaling about as nao^3).
FORCE_ANALYTIC = os.environ.get('FREQ_ANALYTIC') == '1'
use_fd = not FORCE_ANALYTIC and (
    bool(case.get('solvent')) or os.environ.get('FREQ_FD') == '1' or mol.nao > 300)
hess = None
hess_is_fd = False
if not use_fd:
    try:
        print(f'[freq] analytic Hessian ({mol.nao} basis functions, max_memory {MAX_MEMORY} MB) ...', flush=True)
        stage_note(case, f'analytic Hessian attempted ({mol.nao} basis functions)')
        hess = mf.Hessian().kernel()
        method = 'analytic (no finite-difference truncation error)'
        print(f'[freq] analytic Hessian done in {(time.time()-t0)/60:.1f} min', flush=True)
    except Exception as err:
        print(f'[freq] analytic Hessian failed ({type(err).__name__}: {err}); falling back to finite differences', flush=True)
        stage_note(case, f'analytic Hessian failed ({type(err).__name__}); using finite differences')
        hess = None
        use_fd = True

if hess is None:
    ncol = 3 * mol.natm
    method = f'finite-difference (gradients, step {STEP} Bohr, {WORKERS} workers x {WTHREADS} threads)'
    stage_note(case, f'Hessian by parallel finite differences: {2*ncol} gradients, {WORKERS} workers x {WTHREADS} threads')
    print(f'[freq] parallel finite differences: {ncol} columns ({2*ncol} gradients) '
          f'over {WORKERS} workers x {WTHREADS} threads', flush=True)
    tmp = pathlib.Path(tempfile.mkdtemp(prefix='fd_', dir=str(d)))
    edges = [round(i * ncol / WORKERS) for i in range(WORKERS + 1)]
    procs, outs = [], []
    env = dict(os.environ, DFT_THREADS=str(WTHREADS), OMP_NUM_THREADS=str(WTHREADS),
               MKL_NUM_THREADS=str(WTHREADS), OPENBLAS_NUM_THREADS=str(WTHREADS))
    logs = d / 'logs'; logs.mkdir(exist_ok=True)
    # Remove worker logs left by earlier runs of this case.  A run with fewer workers than its
    # predecessor only truncates the logs it reuses, so the extra ones survive with their old column
    # lines in them -- and the status reporter, which counts column lines across fd_worker*.log,
    # then reports a brand-new run as already 24 columns in.  Observed on btpc_pcm, 18 Sep 2026.
    for stale in logs.glob('fd_worker*.log'):
        stale.unlink()
    for w in range(WORKERS):
        if edges[w] == edges[w + 1]: continue
        o = tmp / f'w{w}.npz'; outs.append(o)
        lf = open(logs / f'fd_worker{w}.log', 'w')
        procs.append(subprocess.Popen(
            [sys.executable, str(ROOT / 'dft' / '_fd_worker.py'), case['name'],
             str(edges[w]), str(edges[w + 1]), str(o), str(STEP)],
            stdout=lf, stderr=subprocess.STDOUT, cwd=str(ROOT), env=env))
    done = 0
    for pr in procs:
        rc = pr.wait(); done += 1
        print(f'[freq] worker {done}/{len(procs)} finished (exit {rc}, {(time.time()-t0)/60:.1f} min elapsed)', flush=True)
        if rc != 0:
            raise RuntimeError(f'FD worker failed with exit code {rc}; see {d}/logs/fd_worker*.log')
    hess = np.zeros((mol.natm, mol.natm, 3, 3))
    dipd = np.zeros((ncol, 3))
    got = 0
    for o in outs:
        z = np.load(o)
        for key in z.files:
            if key.startswith('dip'):
                dipd[int(key[3:])] = z[key]
                continue
            col = int(key); a, k = divmod(col, 3)
            hess[a, :, k, :] = z[key]; got += 1
    if got != ncol:
        raise RuntimeError(f'expected {ncol} Hessian columns, collected {got}')
    asym = float(np.abs(hess - hess.transpose(1, 0, 3, 2)).max())
    print(f'[freq] finite-difference Hessian asymmetry before symmetrising: {asym:.2e} Ha/bohr^2', flush=True)
    hess = 0.5 * (hess + hess.transpose(1, 0, 3, 2))
    hess_is_fd = True
    for o in outs: o.unlink(missing_ok=True)
    tmp.rmdir()

res = thermo.harmonic_analysis(mol, hess, imaginary_freq=True)
pyscf_freqs = np.sort(np.array(res['freq_wavenumber']))

# Authoritative frequencies come from our own projection -- see dft/vibrations.py for why PySCF's mode
# selection cannot be trusted on floppy clusters.
from vibrations import frequencies as project_frequencies
freqs, n_imag_proj, tr_res = project_frequencies(mol, hess)
if len(freqs) != len(pyscf_freqs) or not np.allclose(np.sort(freqs), pyscf_freqs, atol=1.0):
    print(f'[freq] NOTE: projected frequencies differ from PySCF. Using the projected set.', flush=True)
    print(f'[freq]   PySCF lowest 6   : {np.round(pyscf_freqs[:6], 1)}', flush=True)
    print(f'[freq]   projected lowest 6: {np.round(np.sort(freqs)[:6], 1)}', flush=True)
    stage_note(case, f'frequency projection disagreed with PySCF; using projected set '
                     f'({n_imag_proj} imaginary vs {int((pyscf_freqs < 0).sum())} reported)')
print(f'[freq] translation/rotation residual: {tr_res:.2f} cm-1 (should be << 1)', flush=True)
th = thermo.thermo(mf, res['freq_au'], 298.15, 101325)
n_imag = n_imag_proj
out = {'case': case['name'], 'E_Ha': float(mf.e_tot), 'n_imaginary': n_imag, 'hessian_method': method,
       'n_basis': int(mol.nao), 'n_atoms': int(mol.natm),
       'frequencies_cm-1': [float(x) for x in np.sort(freqs)], 'hessian_wall_min': (time.time() - t0) / 60,
       'frequencies_pyscf_cm-1': [float(x) for x in pyscf_freqs],
       'n_imaginary_pyscf': int((pyscf_freqs < 0).sum()), 'tr_residual_cm-1': tr_res,
       'fd_hessian_asymmetry': (asym if hess_is_fd else None),
       'ZPE_Ha': float(th['ZPE'][0]), 'H298_Ha': float(th['H_tot'][0]), 'G298_Ha': float(th['G_tot'][0])}
if hess_is_fd:
    from vibrations import ir_intensities
    f_ir, i_ir = ir_intensities(mol, hess, dipd)
    out['ir'] = [{'freq_cm-1': float(a), 'intensity_km_mol': float(b)} for a, b in zip(f_ir, i_ir)]
    np.save(d / 'dipole_derivatives.npy', dipd)
save_json(d / 'freq.json', out)
np.save(d / 'hessian.npy', hess)
stage_note(case, f"frequencies ({method}): {n_imag} imaginary, lowest {np.sort(freqs)[0]:.1f} cm-1, {(time.time()-t0)/60:.0f} min")
print(f"[freq] {len(freqs)} modes, lowest: {np.round(np.sort(freqs)[:6],1)} cm-1; imaginary = {n_imag} "
      f"({'TRUE MINIMUM' if n_imag == 0 else 'NOT a minimum -- re-optimise'}); Hessian took {(time.time()-t0)/60:.1f} min")
