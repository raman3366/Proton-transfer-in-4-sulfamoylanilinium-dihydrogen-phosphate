"""
Shared helpers for the proton-position studies of the M1 ion pair: dft/scan.py (relaxed N1-H1A scan),
dft/eps_ladder.py (dielectric ladder) and dft/zpe_barrier.py (zero-point energy against the barrier).

Why a module of its own: all three must measure the proton in exactly the same way -- which atoms are the
donor, the proton and the acceptor, delta = d(N1-H1A) - d(H1A...O5), and the P-O bond lengths that tell
H2PO4- (two short, two long P-O) from H3PO4 (one short, three long) -- and the scan and the ladder share one
resumable geomeTRIC driver.  One code path means the three studies cannot disagree about what "the proton
moved" means.

IMPORTANT: import this module BEFORE `common` (i.e. before numpy / PySCF load).  These scripts do not go through
ops/vm/run_case.sh, so importing this module does what that script does: it pins every BLAS/OpenMP pool to
DFT_THREADS (unpinned BLAS oversubscribes the machine by an order of magnitude, see run_case.sh), defaults the
PySCF memory to run_case.sh's 6000 MB, and puts PySCF's scratch on the NVMe volume (/mnt/scratch) in a
per-process directory that is removed on exit -- per process, because the two scan chains and the ladder
points of one case run concurrently.
"""
import os, sys, atexit, shutil, signal, tempfile, pathlib


def _pin_environment():
    t = os.environ.get('DFT_THREADS')
    if t:
        for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                  'VECLIB_MAXIMUM_THREADS'):
            os.environ[v] = t
    os.environ.setdefault('DFT_MAX_MEMORY_MB', '6000')
    if 'numpy' in sys.modules:
        print('[proton_tools] WARNING: numpy was imported before proton_tools; thread pinning may not apply',
              file=sys.stderr)
    if 'PYSCF_TMPDIR' not in os.environ and os.path.isdir('/mnt/scratch') and os.access('/mnt/scratch', os.W_OK):
        tmp = tempfile.mkdtemp(prefix=f'{pathlib.Path(sys.argv[0]).stem}_', dir='/mnt/scratch')
        os.environ['PYSCF_TMPDIR'] = os.environ['TMPDIR'] = tmp
        atexit.register(shutil.rmtree, tmp, True)
        # a killed task (SIGTERM) should still remove its scratch: turn the signal into a normal exit
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))


_pin_environment()

import json, time, datetime, fcntl, socket
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *                                     # noqa: E402  (after the pinning above, on purpose)
from pyscf.geomopt import geometric_solver                # noqa: E402

KCAL_PER_HA = 627.509474
# geomeTRIC switches to exact constraint enforcement once the violation is below this (a.u.).  0.1 is the value
# geomeTRIC's own constrained-optimisation tests use; without it a $set distance is only approached to within
# convergence_cmax (default 0.01 A), which is 10% of the 0.10 A scan step.
ENFORCE = 0.1
# geomeTRIC's "GAU_LOOSE" set -- only for code tests at tiny basis (--conv loose); production uses the defaults,
# which equal Gaussian's standard criteria and are what every burst-1a optimisation used.
LOOSE = dict(convergence_energy=1e-6, convergence_grms=1.7e-3, convergence_gmax=2.5e-3,
             convergence_drms=6.7e-3, convergence_dmax=1.0e-2)


# ------------------------------------------------------------------------------------------------ geometry --
def short_label(lab):
    """'O5#1_554' -> 'O5': the symmetry suffix only records which image of the anion the pair uses."""
    return lab.split('#')[0]


def proton_atoms(case):
    """Donor N, proton H, acceptor O (case['proton']) plus the phosphorus bonded to the acceptor and the list of
    its oxygens, all as 0-based indices into the case geometry (atom order = CIF labels)."""
    p = case['proton']
    syms, x0 = read_xyz(ROOT / case['xyz'])
    O = p['acceptor_O']
    P_list = [i for i, s in enumerate(syms) if s == 'P']
    P = min(P_list, key=lambda i: np.linalg.norm(x0[i] - x0[O]))
    PO = [i for i, s in enumerate(syms) if s == 'O' and np.linalg.norm(x0[i] - x0[P]) < 1.9]
    lab = [short_label(l) for l in (case.get('cif_labels') or atom_labels(syms))]
    return {'N': p['donor_N'], 'H': p['H'], 'O': O, 'P': P, 'PO': PO, 'labels': lab, 'syms': syms}


def proton_metrics(atoms, xyz):
    """Everything the proton studies report for one geometry (Angstrom in, Angstrom/degrees out).

    delta = d(N1-H1A) - d(H1A...O5): negative = proton on the anilinium N (salt), positive = on the phosphate O.
    The P-O bond lengths are the independent check: P-O5 lengthens from ~1.50-1.52 A (P=O/P-O-, H2PO4-) to
    ~1.56-1.60 A (P-OH) when O5 is protonated."""
    N, H, O, P = atoms['N'], atoms['H'], atoms['O'], atoms['P']
    lab = atoms['labels']
    dNH = float(np.linalg.norm(xyz[H] - xyz[N])); dHO = float(np.linalg.norm(xyz[H] - xyz[O]))
    v1, v2 = xyz[N] - xyz[H], xyz[O] - xyz[H]
    ang = float(np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))
    out = {'d_N1_H1A_A': dNH, 'd_H1A_O5_A': dHO, 'delta_A': dNH - dHO,
           'd_N1_O5_A': float(np.linalg.norm(xyz[N] - xyz[O])), 'angle_N1_H1A_O5_deg': ang,
           'd_P_O5_A': float(np.linalg.norm(xyz[P] - xyz[O])),
           'P_O_A': {f'{lab[P]}-{lab[i]}': float(np.linalg.norm(xyz[P] - xyz[i])) for i in atoms['PO']},
           'proton_on': 'N' if dNH < dHO else 'O'}
    return out


def o_protonated_guess(atoms, xyz, d_oh=1.00):
    """Starting geometry for the O-protonated form: H1A moved onto the O5->N1 line, d_oh from O5.

    Everything else is left at the input geometry; the optimiser relaxes the phosphate (P-O5 lengthens) and the
    now-neutral amine.  A straight N...H-O line is the natural guess for a proton handed across an N-H...O bond
    that is ~174 deg in the crystal."""
    x = np.array(xyz, dtype=float).copy()
    u = x[atoms['N']] - x[atoms['O']]; u /= np.linalg.norm(u)
    x[atoms['H']] = x[atoms['O']] + d_oh * u
    return x


def n_protonated_guess(atoms, xyz, d_nh=1.03):
    """The mirror image of o_protonated_guess: H1A put back on the N1->O5 line, d_nh from N1.  Used when a
    geometry meant to be the salt (e.g. an optimum in which the proton moved to O5) is not N-protonated."""
    x = np.array(xyz, dtype=float).copy()
    u = x[atoms['O']] - x[atoms['N']]; u /= np.linalg.norm(u)
    x[atoms['H']] = x[atoms['N']] + d_nh * u
    return x


def xyz_block(syms, xyz, comment):
    return f'{len(syms)}\n{comment}\n' + ''.join(f'{s:2s} {a:14.8f} {b:14.8f} {c:14.8f}\n'
                                                 for s, (a, b, c) in zip(syms, xyz))


# -------------------------------------------------------------------------------------- level of theory ----
def level_description(case, basis=None, eps=None):
    solv = case.get('solvent')
    if eps is not None:
        env = 'gas phase' if eps <= 1.0 else f'IEF-PCM eps={eps:g}'
    else:
        env = f'IEF-PCM {solv} (eps={SOLVENT_EPS[solv.lower()]})' if solv else 'gas phase'
    return {'xc': os.environ.get('DFT_XC', XC), 'disp': os.environ.get('DFT_DISP', 'none'),
            'basis': basis or basis_description(case), 'environment': env,
            'grid_level': GRID_LEVEL,
            'pcm_lebedev_order': (int(os.environ.get('PCM_LEBEDEV_ORDER', '29')) if env != 'gas phase' else None),
            'density_fitting': True}


def make_scanner(case, xyz_path, basis=None, eps=None, chkfile=None, verbose=3):
    """Gradient scanner for one chain / ladder point.  One scanner is kept for a whole chain so every new
    point starts its SCF from the converged density of the previous geometry.

    eps: None -> the case's own environment (solvent name or gas); <= 1 -> gas phase; otherwise IEF-PCM with
    that dielectric constant.  The numeric-eps path needs no change to common.make_mf: make_mf builds the
    case's PCM object and only the dielectric constant is replaced before anything is computed (PySCF's PCM
    reads .eps when it builds its K/R matrices, on the first SCF), so cavity, Lebedev order and every other
    setting stay identical to the water runs."""
    if basis:
        mol = build_mol(xyz_path, case['charge'], case['spin'], basis=basis, verbose=verbose,
                        basis_overrides=case.get('basis_overrides'))
    else:
        mol = mol_for_case(case, xyz_path, verbose=verbose)
    solvent = case.get('solvent')
    if eps is not None:
        solvent = None if eps <= 1.0 else (solvent or 'water')
    conv_tol = 1e-8 if solvent else 1e-9            # as dft/optimize.py
    mf = make_mf(mol, conv_tol=conv_tol, solvent=solvent)
    if eps is not None and eps > 1.0:
        mf.with_solvent.eps = float(eps)
    if chkfile:
        mf.chkfile = str(chkfile)
        # resume: take the first SCF guess from the chain's last converged density
        if os.path.exists(chkfile):
            mf.init_guess = 'chkfile'
    return mf.nuc_grad_method().as_scanner(), mol


# ------------------------------------------------------------------------------------- optimisation -------
class StopOptimisation(Exception):
    """Raised from the step callback to end an optimisation whose outcome is already known."""


class NegligibleForces(StopOptimisation):
    """Raised when an optimisation meets one of the two extra acceptance rules of optimise() (see ACCEPTANCE
    below): it is accepted as converged.  .rule says which ('negligible forces' or 'forces + energy')."""
    def __init__(self, msg, rule):
        super().__init__(msg)
        self.rule = rule


# Extra acceptance rules, each required on NEGLIGIBLE_STEPS consecutive gradients (constraint on target, SCF
# converged).  Forces are measured exactly as geomeTRIC measures them (per-atom gradient norms, rms and max over
# atoms, constraint component projected out), so they are the "Grad = rms/max" numbers of its log and of
# dft/finalize_opt.py.
#   "negligible forces": 10x inside the rms/max force thresholds (Gaussian's "optimization completed on the basis
#       of negligible forces" is the model).
#   "forces + energy":   inside the force thresholds and |dE| < convergence_energy (1e-6 Ha) between consecutive
#       gradients: the displacement criterion is the only one left unmet.  geomeTRIC offers the same idea as its
#       Q-Chem-style option (qccnv: grms and either drms or energy, checked once); this rule is stricter (rms AND
#       max force, three gradients in a row).
# Why the second rule (added in review, 28 Sep): the PCM ion pair has nearly free inter-ion rotations on which
# geomeTRIC crawls.  Burst 1a's sfa_pair_pcm_d3 optimisation spent steps 51-64 at rms 3.8e-5..1.9e-4 / max
# 6.9e-5..6.8e-4 Ha/bohr (inside the thresholds from step 55 on), |dE| 5e-7..3.4e-6 Ha per step and dmax
# 3-8e-3 A (> 1.8e-3), and was killed at step 64 (sfa_pair_pcm_b3lyp likewise at 68).  The 10x rule never fires
# on that trajectory (smallest rms 3.8e-5 > 3.0e-5), so every PCM scan point would run into its step cap; the
# second rule would have accepted it at step 59.  The whole 51->64 crawl gained 1.6e-5 Ha (0.01 kcal/mol), far
# below the scan's 0.1 kcal/mol reporting precision.
NEGLIGIBLE_STEPS = 3


def proton_left(atoms, expect, patience=3):
    """stop_if for an unconstrained anchor: True once the proton has sat on the other side for `patience`
    consecutive steps.  An energy minimiser that carries the proton across has shown that no minimum holds it
    on the expected side along that descent, and following it all the way into the other well costs another
    20-40 gradients for a geometry the other chain produces anyway."""
    n = [0]

    def stop(x):
        n[0] = n[0] + 1 if proton_metrics(atoms, x)['proton_on'] != expect else 0
        return n[0] >= patience
    return stop


def optimise(gs, syms, start_xyz, *, constraint=None, cfile=None, partial=None, traj=None, maxsteps=200,
             loose=False, tag='', stop_if=None):
    """One geomeTRIC optimisation with the gradient scanner gs, from start_xyz (Angstrom).

    constraint: None (unconstrained) or (i, j, value_A) with 0-based atom indices -> written to cfile as
        $set
        distance i+1 j+1 value          (geomeTRIC: 1-based atom numbers, Angstrom)
    partial: the latest geometry is rewritten here after every gradient, so an interrupted point restarts from
    where it stopped; traj: every step is appended here.
    stop_if: optional callable(xyz) -> bool checked after every gradient; True ends the run early (result has
    converged=False and 'stopped' set, and is not a stationary point).
    A run that meets one of the two extra acceptance rules above (NEGLIGIBLE_STEPS consecutive gradients) is
    accepted as converged ('accepted_on' = 'negligible forces' or 'forces + energy').
    Returns dict(converged, xyz, E_Ha, steps, wall_s, history, stopped, accepted_on).  An SCF that fails to converge is retried
    once from the last geometry with a fresh guess and a level shift (which changes the path to convergence,
    not the converged energy)."""
    kw = dict(LOOSE) if loose else {}
    g_rms_thr = kw.get('convergence_grms', 3.0e-4)            # geomeTRIC / Gaussian defaults: 3.0e-4, 4.5e-4
    g_max_thr = kw.get('convergence_gmax', 4.5e-4)
    e_thr = kw.get('convergence_energy', 1.0e-6)
    small, flat = [0], [0]
    if constraint is not None:
        i, j, val = constraint
        kw['enforce'] = ENFORCE
        pathlib.Path(cfile).write_text(f'$set\ndistance {i + 1} {j + 1} {val:.6f}\n')
    hist, t0 = [], time.time()
    tf = open(traj, 'a') if traj else None

    def cb(env):
        m, e, g = env['mol'], float(env['energy']), env['gradients']
        # the callback runs before PySCF's own SCF-convergence assertion: never count an unconverged gradient
        scf_ok = bool(getattr(env.get('g_scanner'), 'converged', True))
        x = m.atom_coords() * BOHR_A
        hist.append({'step': len(hist), 'E_Ha': e, 'gmax': float(abs(g).max()),
                     'grms': float(np.sqrt((g ** 2).mean())), 't_s': round(time.time() - t0, 1)})
        if not scf_ok:
            hist[-1]['scf_converged'] = False
        if partial:
            tmp = str(partial) + '.tmp'
            write_xyz(tmp, syms, x, f'{tag} step {len(hist) - 1} E = {e:.10f}')
            os.replace(tmp, partial)
        if tf:
            tf.write(xyz_block(syms, x, f'{tag} step {len(hist) - 1} E = {e:.10f}')); tf.flush()
        extra, gp = '', g
        if constraint is not None:
            # the raw gradient keeps the constraint force along d(i-j) and never goes to zero at a constrained
            # minimum; the gradient with that component projected out is the one to watch
            ci, cj = constraint[0], constraint[1]
            e_ = (x[cj] - x[ci]) / np.linalg.norm(x[cj] - x[ci])
            b = np.zeros_like(g); b[cj] = e_; b[ci] = -e_; b /= np.linalg.norm(b)
            gp = g - (g * b).sum() * b
            extra = (f'  max|g_perp| = {abs(gp).max():.2e}  d = {np.linalg.norm(x[ci] - x[cj]):.5f} A')
        print(f'[{tag}] step {len(hist) - 1:3d}  E = {e:.8f}  max|g| = {abs(g).max():.2e}{extra}', flush=True)
        last['x'] = x
        if stop_if is not None and stop_if(x):
            raise StopOptimisation(f'stopped at step {len(hist) - 1}: stop condition met')
        on_target = constraint is None or abs(np.linalg.norm(x[constraint[0]] - x[constraint[1]]) - constraint[2]) < 1e-3
        # geomeTRIC's force measure: per-atom norms of the (constraint-projected) gradient, rms and max over atoms
        ga = np.linalg.norm(np.asarray(gp).reshape(-1, 3), axis=1)
        c_rms, c_max = float(np.sqrt((ga ** 2).mean())), float(ga.max())
        hist[-1].update({'grms_conv': c_rms, 'gmax_conv': c_max})
        prev_ok = len(hist) >= 2 and hist[-2].get('scf_converged', True)
        dE = (e - hist[-2]['E_Ha']) if prev_ok else None
        ok = scf_ok and on_target
        small[0] = small[0] + 1 if (ok and c_rms < g_rms_thr / 10 and c_max < g_max_thr / 10) else 0
        flat[0] = flat[0] + 1 if (ok and c_rms < g_rms_thr and c_max < g_max_thr and dE is not None
                                  and abs(dE) < e_thr) else 0
        if small[0] >= NEGLIGIBLE_STEPS:
            raise NegligibleForces(f'accepted at step {len(hist) - 1} on negligible forces (rms < {g_rms_thr / 10:.1e}, '
                                   f'max < {g_max_thr / 10:.1e} for {NEGLIGIBLE_STEPS} gradients)', 'negligible forces')
        if flat[0] >= NEGLIGIBLE_STEPS:
            raise NegligibleForces(f'accepted at step {len(hist) - 1} on forces + energy (rms < {g_rms_thr:.1e}, max < '
                                   f'{g_max_thr:.1e}, |dE| < {e_thr:.0e} Ha for {NEGLIGIBLE_STEPS} gradients; '
                                   f'displacement criterion not met)', 'forces + energy')

    last = {}
    x_start = np.asarray(start_xyz, dtype=float)
    stopped = accepted = None
    for attempt in (1, 2):
        gs.reset(gs.mol.set_geom_(x_start, unit='Angstrom', inplace=False))
        try:
            conv, mol_f = geometric_solver.kernel(gs, constraints=(str(cfile) if constraint is not None else None),
                                                  callback=cb, maxsteps=maxsteps, **kw)
            x = mol_f.atom_coords() * BOHR_A
            break
        except NegligibleForces as err:
            print(f'[{tag}] {err}', flush=True)
            conv, accepted, x = True, err.rule, last['x']
            break
        except StopOptimisation as err:
            print(f'[{tag}] {err}', flush=True)
            conv, stopped, x = False, str(err), last['x']
            break
        except RuntimeError as err:
            if attempt == 2 or 'not converged' not in str(err):
                raise
            print(f'[{tag}] SCF failed ({err}); retrying once from the last geometry with a fresh guess', flush=True)
            if partial and os.path.exists(partial):
                x_start = read_xyz(partial)[1]
            gs.base.mo_coeff = None
            gs.base.init_guess = 'minao'
            gs.base.level_shift = 0.2
            gs.base.max_cycle = 300
    gs.base.level_shift = 0.0
    if tf:
        tf.close()
    return {'converged': bool(conv), 'xyz': x, 'E_Ha': hist[-1]['E_Ha'], 'steps': len(hist),
            'wall_s': time.time() - t0, 'history': hist, 'stopped': stopped, 'accepted_on': accepted}


# ------------------------------------------------------------------------------------------ bookkeeping ----
def stamp():
    return {'host': socket.gethostname(), 'time': datetime.datetime.now().isoformat(timespec='seconds'),
            'threads': int(os.environ.get('DFT_THREADS', '0') or 0)}


def atomic_json(path, obj):
    tmp = pathlib.Path(str(path) + '.tmp')
    save_json(tmp, obj)
    os.replace(tmp, path)


class FileLock:
    """Advisory lock so two concurrent chains never write the combined report at the same moment."""
    def __init__(self, path):
        self.fh = open(path, 'w')

    def __enter__(self):
        fcntl.flock(self.fh, fcntl.LOCK_EX); return self

    def __exit__(self, *a):
        fcntl.flock(self.fh, fcntl.LOCK_UN); self.fh.close()


def load_points(scan_dir, profile_only=True):
    """Finished scan points (per-point JSON written by dft/scan.py), both chains.  profile_only drops anchors
    that were stopped early (proton_left): they are not stationary points and do not belong on E(d)."""
    pts = []
    for p in sorted(pathlib.Path(scan_dir).glob('[AB]_[0-9][0-9]_*.json')):
        r = json.load(open(p))
        if r.get('status') == 'done' and not (profile_only and r.get('stopped')):
            pts.append(r)
    return pts


def lower_envelope(pts, tol=2e-3):
    """Lowest-energy point at each value of d(N1-H1A) (chains A and B share the grid; anchors add their own
    d).  A relaxed scan can follow different valleys from its two ends, so the lower envelope is the best
    available estimate of the minimum-energy profile along d.  Returns rows sorted by d."""
    rows = sorted(pts, key=lambda r: (r['metrics']['d_N1_H1A_A'], r['E_Ha']))
    env = []
    for r in rows:
        d = r['metrics']['d_N1_H1A_A']
        if env and abs(env[-1]['metrics']['d_N1_H1A_A'] - d) < tol:
            if r['E_Ha'] < env[-1]['E_Ha']:
                env[-1] = r
        else:
            env.append(r)
    return env
