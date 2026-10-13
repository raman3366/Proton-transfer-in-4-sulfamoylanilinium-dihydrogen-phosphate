"""
M9 (referee ref_chem.md M9; ref_dock.md M5, M6): relaxed B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM (water) scans of the
sulfonamide torsion of sulfanilamide (neutral) and of the sulfonamidate anion.  They replace the rigid single points
of the docked conformers (burst 2, results/<case>/dock_strain.json, SI Table S27), whose common 3-4.5 kcal/mol offset
came from MMFF94s-placed N-H hydrogens on a rigid DFT frame, by constrained optimisations: the torsion is held at
each value, every other coordinate relaxes, and E(tau) - E(minimum) is the torsional strain.

Usage (run from the repository root; ~/cqc-sfa on the VM)
  python dft/torsion_scan.py plan [--write] [--fragment PATH]   # Pi: requested torsions -> the plan (+ fragment)
  python dft/torsion_scan.py check                              # gate: plan, minima and every start geometry
  python dft/torsion_scan.py ref <case> [--max-min M]           # unconstrained re-optimisation of opt.xyz
  python dft/torsion_scan.py point <case> --tau T [--branch C5] [--max-min M]
  python dft/torsion_scan.py report [<case> ...]                # per-case and combined tables (cheap)
Cases: sulfanilamide_pcm_d3 (docking species sfa_neutral), sulfonamidate_pcm_d3 (sfa_anion).
Test-only options: --basis B, --conv loose, --out DIR (every output under DIR; nothing is written to results/).
Threads: DFT_THREADS (dft/proton_tools pins BLAS/OpenMP to it and puts the PySCF scratch on /mnt/scratch).

The coordinate
  tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|), exactly as dock/RESULTS.md section 6, Fig. 11 and Table S27 define
  it (ring-flip invariant; 90 deg = S1-N2 perpendicular to the ring plane, 180 deg = in it).  Atom indices are looked
  up from the case's cif_labels (0-based, the same in both cases): C3 2, C4 3, C5 4, N2 7, S1 10.
  The constrained coordinate is C3-C4-S1-N2 itself, held at +tau with geomeTRIC's "$set / dihedral 3 4 11 8 tau"
  (exact enforcement once within 0.1 rad, proton_tools.ENFORCE).  Both minima have C3-C4-S1-N2 > 0 (neutral +87.16,
  whose tau_fold 90.44 comes from C5; anion +101.70), so +tau for tau in [90, 180] is one continuous branch from
  each minimum.  For tau >= ~93 deg tau_fold = |C3-C4-S1-N2| (|C5-C4-S1-N2| = 180 - tau within the 2-3 deg
  non-planarity at C4); the measured tau_fold of every point is reported beside the target.
  One branch suffices: the other ring-flip branch (C5-C4-S1-N2 = -tau) is the mirror image of this one through the
  plane that contains C1...C4 and is perpendicular to the ring (it swaps C3/C5, O1/O2, H1B/H1C), so the relaxed
  energies of the two branches are equal by symmetry.  `point --branch C5` runs the mirror branch as a check
  (optional task m9_sfa_b150).
Starts
  Every point is independent, so all can run at once: the species' unconstrained PCM minimum results/<case>/opt.xyz
  (burst 1a) with the S1 side of the C4-S1 bond (S1, O1, O2, N2 and its H) rotated rigidly about C4-S1 to the target.
  No chain: no point inherits another's conformer or error.  The rest of the molecule relaxes from the minimum's
  conformation (amine and sulfonamide NH2 orientations included), so each point is the constrained minimum in that
  basin; a point that falls into another basin shows up as an outlier of the smooth profile, and one that ends below
  the reference is flagged (below_reference).
Reference
  `ref` re-optimises opt.xyz without constraint, with the same driver and settings (a few gradients: it is already
  the minimum).  strain(tau) = E(tau) - E(ref).  Until ref exists the report falls back to opt_summary.json's
  E_final_Ha (the same level: grid 4, Lebedev 29, DF; burst 2 reproduced it to 5e-9 Ha) and says which it used.
Level of theory
  The case's: B3LYP-D3(BJ)/6-311++G(d,p), density fitting, grid level 4, IEF-PCM water with Lebedev order 29 (set
  GRID_LEVEL / PCM_LEBEDEV_ORDER in the task env, as run_case.sh does).
Convergence
  geomeTRIC defaults (= Gaussian's) plus dft/scan.py's two extra acceptance rules (proton_tools: NEGLIGIBLE_STEPS
  gradients in a row with forces 10x inside the thresholds, or inside them with |dE| < 1e-6 Ha), measured on the
  gradient with the constraint direction (the dihedral's Wilson B vector) projected out, so a point cannot crawl on
  a soft NH2 rotation.  Step cap --maxsteps (60): a capped point is kept and flagged unconverged.
  --max-min: wall-clock limit of this process; at the limit the point stops cleanly with status 'time_limit' and
  its .partial.xyz / .chk (exit status 3, so the burst shows it as not done), and rerunning the same command
  resumes it.
Requested torsions (the plan: molecules/sfa/input/torsion_scan_plan.json, written by `plan --write` on the Pi)
  grid 90..180 every 10 deg; the co-crystal ligands measured from dock/T*/ref_*.sdf (heavy atoms): 6RL9
  sulfanilamide 132.1, 1AJ0 sulfanilamide 111.3, 3TZF sulfamethoxazole 105.6 (the referee's values, re-derived
  here); the docked tau_fold of every production job of the species (facts bulk.docktors.*.tau_fold.mean) and the
  Table S27 representative pose of each job (dock/results/poses_for_dft/manifest.json); the SAPH crystal (102.66).
  Each entry is 'exact' (a task computes it; tier critical or optional) or 'interp' (read off a cubic spline through
  every converged point of the species, x = the constrained C3-C4-S1-N2 value, the minimum included at its own
  C3-C4-S1-N2).  A requested value within 0.5 deg of an exact point is interpolated, not recomputed.  The report
  checks the interpolation: a spline through the grid points and the minimum alone predicts every exact off-grid
  point, and the largest miss is printed (interp_check_max_abs_kcal).
Outputs
  results/<case>/torsion_scan/: tau<100*tau>[_C5].{json,xyz,constraints}, *_traj.xyz, ref.{json,xyz},
    torsion_scan.csv, torsion_scan_summary.json.
  molecules/sfa/results/torsion_scan/: torsion_strain.{json,csv,md} (both species: every requested torsion with its
    relaxed strain, exact or interpolated, beside the MMFF94s and rigid-DFT strain of the same pose where one exists)
    and torsion_profile.png (a working plot).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from proton_tools import *     # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
from proton_tools import StopOptimisation, NegligibleForces, NEGLIGIBLE_STEPS, ENFORCE, LOOSE
import argparse, math

CASES = {'sulfanilamide_pcm_d3': 'sfa_neutral', 'sulfonamidate_pcm_d3': 'sfa_anion'}
SPECIES_LABEL = {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate'}
PLAN = ROOT / 'molecules' / 'sfa' / 'input' / 'torsion_scan_plan.json'
COMBINED = ROOT / 'molecules' / 'sfa' / 'results' / 'torsion_scan'
GRID = [90.0 + 10.0 * k for k in range(10)]
MERGE_DEG = 0.5          # a requested value this close to an exact point is interpolated, not recomputed
EXACT_TOL = 0.06         # a requested value matches a computed target within this (targets are rounded to 0.1 deg)
ON_TARGET_DEG = 0.05     # the constraint counts as met (extra acceptance rules) within this
CO_CRYSTAL = [           # (target dir, sdf, PDB, protein, ligand, the referee's value from ref_chem.md M9c)
    ('T3_6RL9', 'ref_SAN.sdf', '6RL9', 'hCA II', 'sulfanilamide', 132.1),
    ('T2_1AJ0', 'ref_SAN.sdf', '1AJ0', 'E. coli DHPS', 'sulfanilamide', 111.3),
    ('T1_3TZF', 'ref_08D.sdf', '3TZF', 'Y. pestis DHPS', 'sulfamethoxazole', 105.6)]
TARGET_NAME = {'T1_3TZF': '3TZF (DHPS)', 'T3_6RL9': '6RL9 (hCA II)', 'T4_3HS4': '3HS4 (hCA II)',
               'T5_7S2K': '7S2K (Sul2)', 'T2_1AJ0': '1AJ0 (DHPS)'}
SCORING_NAME = {'vina': 'Vina', 'vinardo': 'Vinardo', 'ad4zn': 'AD4Zn'}
# Which docked values the plan computes exactly (tier critical): the ones ref_chem.md M9c names -- 6RL9 AutoDock4Zn
# mean 147.0, 6RL9 Vina ~160, DHPS 107.7 (3TZF Vina representative pose) -- and, for the anion, the DHPS and 6RL9
# Vina means.  Every other docked job mean is an optional task (or interpolated when within MERGE_DEG of an exact
# point); the representative poses are interpolated.
CRITICAL_DOCKED = {'sfa_neutral': [('mean', 'T3_6RL9', 'ad4zn'), ('mean', 'T3_6RL9', 'vina'),
                                   ('rep', 'T1_3TZF', 'vina')],
                   'sfa_anion': [('mean', 'T1_3TZF', 'vina'), ('mean', 'T3_6RL9', 'vina')]}
# Cost model for the fragment's est_min (see the fragment notes): minutes per gradient at 4 threads on the shared
# 48-vCPU VM, and gradients per point as a function of the rigid rotation |delta| from the minimum.  Pi test
# (29 Sep 2026, B3LYP-D3(BJ)/6-31G*, grid 3, PCM Lebedev 17, default convergence): 13 gradients at |delta| 44.9
# (132.1 deg), 18 at 92.8 (180 deg), i.e. 8.3 + 0.104 |delta|; one gradient more is allowed for the production
# basis.  The fragment multiplies by 1.3.
MIN_PER_GRAD = 2.4
def est_steps(delta_deg):
    return 9 + 0.11 * abs(delta_deg)


# ------------------------------------------------------------------------------------------------ geometry --
def dihedral(x, i, j, k, l):
    """Dihedral i-j-k-l in degrees, (-180, 180] (IUPAC sign)."""
    b0, b1, b2 = x[i] - x[j], x[k] - x[j], x[l] - x[k]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1n) * b1n
    w = b2 - (b2 @ b1n) * b1n
    return float(np.degrees(np.arctan2(np.cross(b1n, v) @ w, v @ w)))


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def dihedral_bvec(x, idx, h=1e-4):
    """d(dihedral)/d(x) (rad per Angstrom, flattened), central differences: the constraint direction."""
    b = np.zeros(x.size)
    for a in idx:
        for c in range(3):
            xp, xm = x.copy(), x.copy()
            xp[a, c] += h; xm[a, c] -= h
            b[3 * a + c] = math.radians(wrap(dihedral(xp, *idx) - dihedral(xm, *idx))) / (2 * h)
    return b


def angle(x, i, j, k):
    u, v = x[i] - x[j], x[k] - x[j]
    return float(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def neighbours(syms, x, i):
    return [j for j in range(len(syms)) if j != i and
            np.linalg.norm(x[i] - x[j]) < 1.25 * (COV_R.get(syms[i], 0.8) + COV_R.get(syms[j], 0.8))]


def side_atoms(syms, x, a, b):
    """Atoms on b's side of the a-b bond (b included), by the covalent-radius connectivity of x."""
    seen, stack = {b}, [b]
    while stack:
        i = stack.pop()
        for j in neighbours(syms, x, i):
            if j != a and j not in seen:
                seen.add(j); stack.append(j)
    assert a not in seen, 'the a-b bond is in a ring'
    return sorted(seen)


def rotate(x, a, b, idx, deg):
    """Rigid rotation of atoms idx about the axis a->b by deg (Rodrigues)."""
    k = (x[b] - x[a]) / np.linalg.norm(x[b] - x[a])
    th = math.radians(deg)
    p = x[idx] - x[b]
    p = p * math.cos(th) + np.cross(k, p) * math.sin(th) + np.outer(p @ k, k) * (1 - math.cos(th))
    y = x.copy(); y[idx] = p + x[b]
    return y


class Atoms:
    """Named atoms of one case (0-based indices from the case's cif_labels)."""
    def __init__(self, case):
        lab = [short_label(l) for l in case['cif_labels']]
        self.syms, x0 = read_xyz(ROOT / case['xyz'])
        assert len(lab) == len(self.syms)
        ix = {l: i for i, l in enumerate(lab)}
        self.C1, self.C3, self.C4, self.C5 = ix['C1'], ix['C3'], ix['C4'], ix['C5']
        self.S1, self.N1, self.N2 = ix['S1'], ix['N1'], ix['N2']
        self.labels = lab
        self.T3 = (self.C3, self.C4, self.S1, self.N2)
        self.T5 = (self.C5, self.C4, self.S1, self.N2)

    def constrained(self, branch):
        return self.T3 if branch == 'C3' else self.T5

    def metrics(self, x):
        t3, t5 = dihedral(x, *self.T3), dihedral(x, *self.T5)
        m = {'tau_C3_deg': t3, 'tau_C5_deg': t5, 'tau_fold_deg': max(abs(t3), abs(t5)),
             'd_C4_S1_A': float(np.linalg.norm(x[self.C4] - x[self.S1])),
             'd_S1_N2_A': float(np.linalg.norm(x[self.S1] - x[self.N2])),
             'a_C4_S1_N2_deg': angle(x, self.C4, self.S1, self.N2),
             # out-of-plane angle of S1 at C4 (0 = S1 in the plane of C3, C4, C5)
             'S1_out_of_plane_deg': self._oop(x, self.S1, self.C4, self.C3, self.C5)}
        for nm, n in (('N1', self.N1), ('N2', self.N2)):
            nb = neighbours(self.syms, x, n)
            m[f'{nm}_angle_sum_deg'] = (sum(angle(x, p, n, q) for ii, p in enumerate(nb) for q in nb[ii + 1:])
                                        if len(nb) == 3 else None)
        # orientation of the sulfonamide N-H bonds: C4-S1-N2-H for each H on N2
        m['tau_C4_S1_N2_H_deg'] = {self.labels[h]: dihedral(x, self.C4, self.S1, self.N2, h)
                                   for h in neighbours(self.syms, x, self.N2) if self.syms[h] == 'H'}
        return m

    @staticmethod
    def _oop(x, a, c, n1, n2):
        nrm = np.cross(x[n1] - x[c], x[n2] - x[c]); nrm /= np.linalg.norm(nrm)
        v = x[a] - x[c]
        return float(np.degrees(np.arcsin(np.clip(v @ nrm / np.linalg.norm(v), -1, 1))))


def rotated_start(atoms, x0, tau, branch='C3'):
    """x0 with the S1 side of C4-S1 rotated rigidly so that the constrained dihedral equals its target:
    C3-C4-S1-N2 = +tau (branch C3) or C5-C4-S1-N2 = -tau (branch C5, the mirror branch)."""
    idx = atoms.constrained(branch)
    target = tau if branch == 'C3' else -tau
    side = side_atoms(atoms.syms, x0, atoms.C4, atoms.S1)
    delta = wrap(target - dihedral(x0, *idx))
    for sgn in (1, -1):
        x = rotate(x0, atoms.C4, atoms.S1, side, sgn * delta)
        if abs(wrap(dihedral(x, *idx) - target)) < 1e-6:
            break
    else:
        raise RuntimeError('rigid rotation did not reach the target dihedral')
    # no clash: the rotated group keeps every non-bonded contact above 1.5 A (H...H included)
    rest = [i for i in range(len(x)) if i not in side and i != atoms.C4]
    dmin = min(np.linalg.norm(x[i] - x[j]) for i in side for j in rest)
    assert dmin > 1.5, f'rotated start has a {dmin:.2f} A contact'
    return x, delta, float(dmin)


# ------------------------------------------------------------------------------------------------ optimiser --
class TimeLimit(StopOptimisation):
    """Raised from the step callback at the wall-clock limit (--max-min): the point is resumable."""


def optimise_torsion(gs, syms, start_xyz, *, dih=None, cfile=None, partial=None, traj=None, maxsteps=60,
                     loose=False, tag='', deadline=None):
    """One geomeTRIC optimisation with the gradient scanner gs from start_xyz (Angstrom): proton_tools.optimise
    with a dihedral instead of a distance constraint.

    dih: None (unconstrained) or (i, j, k, l, value_deg), 0-based -> cfile "$set / dihedral i+1 j+1 k+1 l+1 value".
    The extra acceptance rules of proton_tools.optimise apply to the gradient with the constraint direction (the
    dihedral's B vector) projected out; they need the dihedral within ON_TARGET_DEG of its value.
    Returns dict(converged, xyz, E_Ha, steps, wall_s, history, stopped, accepted_on, time_limit)."""
    kw = dict(LOOSE) if loose else {}
    g_rms_thr = kw.get('convergence_grms', 3.0e-4)
    g_max_thr = kw.get('convergence_gmax', 4.5e-4)
    e_thr = kw.get('convergence_energy', 1.0e-6)
    small, flat = [0], [0]
    if dih is not None:
        i, j, k, l, val = dih
        kw['enforce'] = ENFORCE
        pathlib.Path(cfile).write_text(f'$set\ndihedral {i + 1} {j + 1} {k + 1} {l + 1} {val:.6f}\n')
    hist, t0, last = [], time.time(), {}
    tf = open(traj, 'a') if traj else None

    def cb(env):
        m, e, g = env['mol'], float(env['energy']), np.asarray(env['gradients'])
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
        gp, extra, on_target = g.reshape(-1), '', True
        if dih is not None:
            b = dihedral_bvec(x, dih[:4]); b /= np.linalg.norm(b)
            gp = gp - (gp @ b) * b
            cur = dihedral(x, *dih[:4])
            on_target = abs(wrap(cur - dih[4])) < ON_TARGET_DEG
            extra = f'  max|g_perp| = {abs(gp).max():.2e}  dihedral = {cur:.3f}'
            hist[-1]['dihedral_deg'] = cur
        print(f'[{tag}] step {len(hist) - 1:3d}  E = {e:.8f}  max|g| = {abs(g).max():.2e}{extra}', flush=True)
        last['x'] = x
        ga = np.linalg.norm(gp.reshape(-1, 3), axis=1)
        c_rms, c_max = float(np.sqrt((ga ** 2).mean())), float(ga.max())
        hist[-1].update({'grms_conv': c_rms, 'gmax_conv': c_max})
        prev_ok = len(hist) >= 2 and hist[-2].get('scf_converged', True)
        dE = (e - hist[-2]['E_Ha']) if prev_ok else None
        ok = scf_ok and on_target
        small[0] = small[0] + 1 if (ok and c_rms < g_rms_thr / 10 and c_max < g_max_thr / 10) else 0
        flat[0] = flat[0] + 1 if (ok and c_rms < g_rms_thr and c_max < g_max_thr and dE is not None
                                  and abs(dE) < e_thr) else 0
        if small[0] >= NEGLIGIBLE_STEPS:
            raise NegligibleForces(f'accepted at step {len(hist) - 1} on negligible forces', 'negligible forces')
        if flat[0] >= NEGLIGIBLE_STEPS:
            raise NegligibleForces(f'accepted at step {len(hist) - 1} on forces + energy (displacement criterion '
                                   f'not met)', 'forces + energy')
        if deadline is not None and time.time() > deadline:
            raise TimeLimit(f'wall-clock limit reached after step {len(hist) - 1}')

    x_start = np.asarray(start_xyz, dtype=float)
    stopped = accepted = None
    tlim = False
    for attempt in (1, 2):
        gs.reset(gs.mol.set_geom_(x_start, unit='Angstrom', inplace=False))
        try:
            conv, mol_f = geometric_solver.kernel(gs, constraints=(str(cfile) if dih is not None else None),
                                                  callback=cb, maxsteps=maxsteps, **kw)
            x = mol_f.atom_coords() * BOHR_A
            break
        except NegligibleForces as err:
            print(f'[{tag}] {err}', flush=True)
            conv, accepted, x = True, err.rule, last['x']
            break
        except TimeLimit as err:
            print(f'[{tag}] {err}', flush=True)
            conv, stopped, x, tlim = False, str(err), last['x'], True
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
            'wall_s': time.time() - t0, 'history': hist, 'stopped': stopped, 'accepted_on': accepted,
            'time_limit': tlim}


# ------------------------------------------------------------------------------------------------- the plan --
def sdf_tau_fold(path):
    """tau_fold of a sulfonamide ligand from a V2000 SDF (heavy atoms suffice): S bonded to an aromatic C (C4) and
    an N (N2); C3/C5 = the ring neighbours of C4."""
    L = pathlib.Path(path).read_text().splitlines()
    na, nb = int(L[3][:3]), int(L[3][3:6])
    sy = [l[31:34].strip() for l in L[4:4 + na]]
    x = np.array([[float(l[0:10]), float(l[10:20]), float(l[20:30])] for l in L[4:4 + na]])
    nbr = {i: set() for i in range(na)}
    for l in L[4 + na:4 + na + nb]:
        a, b = int(l[:3]) - 1, int(l[3:6]) - 1
        nbr[a].add(b); nbr[b].add(a)
    (s,) = [i for i in range(na) if sy[i] == 'S']
    (c4,) = [j for j in nbr[s] if sy[j] == 'C']
    (n2,) = [j for j in nbr[s] if sy[j] == 'N']
    c3, c5 = sorted(j for j in nbr[c4] if sy[j] == 'C')
    t3, t5 = dihedral(x, c3, c4, s, n2), dihedral(x, c5, c4, s, n2)
    return max(abs(t3), abs(t5)), t3, t5


def build_plan():
    facts = json.loads((ROOT / 'paper' / 'main' / 'data' / 'facts.json').read_text())
    man = json.loads((ROOT / 'dock' / 'results' / 'poses_for_dft' / 'manifest.json').read_text())
    xtal = []
    for tdir, sdf, pdb, prot, lig, ref_val in CO_CRYSTAL:
        tf, t3, t5 = sdf_tau_fold(ROOT / 'dock' / tdir / sdf)
        assert abs(tf - ref_val) < 0.06, (pdb, tf, ref_val)
        xtal.append({'tau': round(tf, 1), 'tau_full': tf, 'kind': 'co-crystal', 'label': f'{pdb} {prot}, bound {lig}',
                     'source': f'dock/{tdir}/{sdf} (heavy atoms; C3-C4-S1-N2 {t3:.2f}, C5-C4-S1-N2 {t5:.2f})'})
    saph = facts['dock.tors.crystal']['value']
    plan = {'purpose': 'M9 requested torsions (dft/torsion_scan.py); written by `torsion_scan.py plan --write`',
            'generated': datetime.datetime.now().isoformat(timespec='seconds'),
            'definition': 'tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|); computed as C3-C4-S1-N2 = +tau',
            'merge_deg': MERGE_DEG, 'cases': {}}
    for case, sp in CASES.items():
        req = [{'tau': g, 'tau_full': g, 'kind': 'grid', 'label': f'grid {g:.0f}', 'source': 'grid'} for g in GRID]
        req += [dict(r) for r in xtal]
        req.append({'tau': round(saph, 1), 'tau_full': saph, 'kind': 'SAPH crystal',
                    'label': 'SAPH crystal (cation, CIF)', 'source': 'facts dock.tors.crystal'})
        for k, v in sorted(facts.items()):
            pre = 'bulk.docktors.'
            if k.startswith(pre) and k.endswith('.tau_fold.mean') and f'.{sp}__' in k:
                tgt, job = k[len(pre):].split('.')[:2]
                scoring = job.split('__')[1]
                req.append({'tau': round(v['value'], 1), 'tau_full': v['value'], 'kind': 'docked mean',
                            'label': f'{TARGET_NAME[tgt]}, {SCORING_NAME[scoring]}, mean of the top poses',
                            'target': tgt, 'scoring': scoring, 'source': f'facts {k}'})
        for p in man['poses']:
            if p['species'] != sp:
                continue
            req.append({'tau': round(p['tau_fold'], 1), 'tau_full': p['tau_fold'], 'kind': 'docked pose',
                        'label': f"{TARGET_NAME[p['target']]}, {SCORING_NAME[p['scoring']]}, representative "
                                 f"(seed {p['seed']}, Table S27)",
                        'target': p['target'], 'scoring': p['scoring'], 'pose': p['file'][:-4],
                        'strain_mmff_relaxed': p['strain_mmff_relaxed'], 'strain_mmff_rigid': p['strain_mmff_rigid'],
                        'source': 'dock/results/poses_for_dft/manifest.json'})
        crit = CRITICAL_DOCKED[sp]
        # decide what is computed: grid and co-crystal always (critical), then the named docked values (critical),
        # then the other job means (optional) unless within MERGE_DEG of something already exact
        exact = []

        def near(r):              # the unrounded value decides (110.45 is 0.45 deg from 110, not 0.5)
            return any(abs(r['tau_full'] - e) < MERGE_DEG for e in exact)
        for r in req:
            if r['kind'] in ('grid', 'co-crystal'):
                r['tier'] = 'critical'
                r['compute'] = 'interp' if near(r) else 'exact'
                if r['compute'] == 'exact':
                    exact.append(r['tau_full'])
        for r in req:
            key = ('mean' if r['kind'] == 'docked mean' else 'rep' if r['kind'] == 'docked pose' else None,
                   r.get('target'), r.get('scoring'))
            if key in crit:
                r['tier'] = 'critical'
                r['compute'] = 'interp' if near(r) else 'exact'
                if r['compute'] == 'exact':
                    exact.append(r['tau_full'])
        for r in req:
            if 'compute' in r:
                continue
            if r['kind'] == 'docked mean' and not near(r):
                r['compute'], r['tier'] = 'exact', 'optional'
                exact.append(r['tau_full'])
            else:
                r['compute'], r['tier'] = 'interp', '-'
        # tasks: one per distinct exact tau (the first entry that asked for it names the tier)
        tasks = {}
        for r in req:
            if r['compute'] == 'exact':
                t = tasks.setdefault(r['tau'], {'tau': r['tau'], 'tier': r['tier'], 'labels': [], 'kinds': []})
                if r['tier'] == 'critical':
                    t['tier'] = 'critical'
                t['labels'].append(r['label'])
                t['kinds'].append(r['kind'])
        plan['cases'][case] = {'species': sp, 'requested': sorted(req, key=lambda r: (r['tau'], r['kind'])),
                               'points': sorted(tasks.values(), key=lambda t: t['tau'])}
    return plan


def load_plan():
    return json.loads(PLAN.read_text()) if PLAN.exists() else None


def plan_labels(case_name, tau):
    p = load_plan()
    if not p:
        return []
    for t in p['cases'].get(case_name, {}).get('points', []):
        if abs(t['tau'] - tau) < EXACT_TOL:
            return t['labels']
    return []


# ----------------------------------------------------------------------------------------------- fragment --
def fragment(plan):
    """ops/vm/bursts/fragments/m9.json: the burst3b task objects of this plan (see the fragment notes)."""
    env = {'DFT_THREADS': 4, 'GRID_LEVEL': 4, 'PCM_LEBEDEV_ORDER': 29, 'DFT_MAX_MEMORY_MB': 6000}
    py = 'tools/env/bin/python dft/torsion_scan.py'
    # time guard (as m2): the gate stamps status/m9_t0; points stop cleanly (resumable) 2 h 50 min after it
    left = ('E=$(( $(date +%s) - $(cat status/m9_t0 2>/dev/null || date +%s) )); L=$(( (10200 - E) / 60 )); ')
    short = {'sulfanilamide_pcm_d3': 'sfa', 'sulfonamidate_pcm_d3': 'anion'}
    tasks = [{'id': 'm9_gate',
              'label': 'M9 gate: dft/torsion_scan.py compiles; plan, both PCM minima (opt.xyz, opt_summary.json) '
                       'present; every start geometry builds (rigid rotation exact, no contact < 1.5 A); stamps '
                       'status/m9_t0 (clock of the m9 time guards).  MUST depend on the burst3a waiter.',
              'threads': 1, 'est_min': 1, 'critical': True, 'deps': [],
              'cmd': f'tools/env/bin/python -m py_compile dft/torsion_scan.py && {py} check && mkdir -p status '
                     f'&& date +%s > status/m9_t0', 'env': dict(env, DFT_THREADS=1)}]
    crit_ids, opt_ids = [], []
    for case, c in plan['cases'].items():
        s = short[case]
        tasks.append({'id': f'm9_{s}_ref', 'threads': 4, 'est_min': 15, 'critical': True, 'deps': ['m9_gate'],
                      'tier': 'core',
                      'label': f'M9 {SPECIES_LABEL[c["species"]]}: reference = unconstrained re-optimisation of '
                               f'results/{case}/opt.xyz (B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM water, grid 4, '
                               f'Lebedev 29), 4 threads',
                      'cmd': f'{left}[ $L -ge 15 ] || L=15; timeout -s TERM $(( L * 60 + 900 )) {py} ref {case} '
                             f'--max-min $L', 'env': env})
        crit_ids.append(f'm9_{s}_ref')
    x0 = {case: read_xyz(load_case(case)['dir'] / 'opt.xyz')[1] for case in plan['cases']}
    rows = []
    for case, c in plan['cases'].items():
        atoms = Atoms(load_case(case))
        for p in c['points']:
            delta = abs(wrap(p['tau'] - dihedral(x0[case], *atoms.T3)))
            est = math.ceil(1.3 * MIN_PER_GRAD * est_steps(delta) / 5) * 5
            rows.append((p['tier'] != 'critical', case, p, est, delta))
    rows.sort(key=lambda r: (r[0], -r[3]))            # critical first, longest first within a tier
    for opt, case, p, est, delta in rows:
        s, c = short[case], plan['cases'][case]
        # trimming aid for the lead (burst.py ignores the key): core = grid + co-crystal torsions, named_docked = the
        # docked values ref_chem.md M9c names (interpolable from core, which the report checks), optional = the rest
        tier = ('optional' if opt else 'core' if set(p['kinds']) & {'grid', 'co-crystal'} else 'named_docked')
        tid = f"m9_{s}_{int(round(p['tau'] * 10)):04d}"
        what = '; '.join(p['labels'])
        if not opt:
            guard = f'{left}[ $L -ge 15 ] || L=15; '
        else:
            guard = (f'{left}[ $L -ge {max(30, est)} ] || {{ echo "m9: only $L min left: {tid} not started"; '
                     f'exit 1; }}; ')
        tasks.append({'id': tid, 'threads': 4, 'est_min': est, 'critical': not opt, 'deps': ['m9_gate'], 'tier': tier,
                      'label': f"M9 {SPECIES_LABEL[c['species']]}: C3-C4-S1-N2 held at {p['tau']:.1f} deg, all else "
                               f"relaxed, from opt.xyz rotated by {delta:.1f} deg ({what}); 4 threads; resumable",
                      'cmd': f"{guard}timeout -s TERM $(( L * 60 + 900 )) {py} point {case} --tau {p['tau']:.1f} "
                             f"--max-min $L", 'env': env})
        (opt_ids if opt else crit_ids).append(tid)
    # the mirror-branch check (neutral, 150 deg): equal energies by symmetry
    d150 = abs(wrap(-150.0 - dihedral(x0['sulfanilamide_pcm_d3'], *Atoms(load_case('sulfanilamide_pcm_d3')).T5)))
    est = math.ceil(1.3 * MIN_PER_GRAD * est_steps(d150) / 5) * 5
    tasks.append({'id': 'm9_sfa_b150', 'threads': 4, 'est_min': est, 'critical': False, 'deps': ['m9_gate'],
                  'tier': 'optional',
                  'label': 'M9 check (optional): sulfanilamide, mirror branch C5-C4-S1-N2 = -150 deg; should equal the '
                           '150 deg point of the C3 branch (symmetry argument of torsion_scan.py)',
                  'cmd': f'{left}[ $L -ge {max(30, est)} ] || {{ echo "m9: only $L min left: m9_sfa_b150 not '
                         f'started"; exit 1; }}; timeout -s TERM $(( L * 60 + 900 )) {py} point sulfanilamide_pcm_d3 '
                         f'--tau 150 --branch C5 --max-min $L', 'env': env})
    opt_ids.append('m9_sfa_b150')
    tasks.append({'id': 'm9_report', 'threads': 0, 'est_min': 1, 'critical': True, 'deps': crit_ids,
                  'label': 'M9 report: per-case torsion_scan.csv / _summary.json and the combined '
                           'molecules/sfa/results/torsion_scan/torsion_strain.{json,csv,md} + torsion_profile.png '
                           '(also refreshed by every point task)',
                  'cmd': f'{py} report', 'env': {'DFT_THREADS': 1}})
    return tasks, crit_ids, opt_ids


# ------------------------------------------------------------------------------------------------ run modes --
def run_one(args, kind):
    case = load_case(args.case)
    assert args.case in CASES, f'{args.case}: not one of {list(CASES)}'
    atoms = Atoms(case)
    OUT = pathlib.Path(args.out) / args.case if args.out else case['dir'] / 'torsion_scan'
    OUT.mkdir(parents=True, exist_ok=True)
    src = case['dir'] / 'opt.xyz'
    s, x0 = read_xyz(src)
    assert s == atoms.syms, 'opt.xyz atom order differs from the case'
    if kind == 'ref':
        stem, dih, delta, how = OUT / 'ref', None, 0.0, f'{src} (unconstrained)'
        tag = f'{args.case} ref'
        xs = x0
    else:
        branch = args.branch
        tau = round(float(args.tau), 1)
        assert 90.0 <= tau <= 180.0, 'tau_fold must lie in [90, 180]'
        idx = atoms.constrained(branch)
        target = tau if branch == 'C3' else -tau
        xs, delta, dmin = rotated_start(atoms, x0, tau, branch)
        stem = OUT / (f'tau{int(round(tau * 100)):05d}' + ('' if branch == 'C3' else '_C5'))
        dih = (*idx, target)
        how = f'{src} with the S1 side of C4-S1 rotated rigidly by {delta:+.2f} deg (closest contact {dmin:.2f} A)'
        tag = f'{args.case} tau {tau:.1f}{"" if branch == "C3" else " C5"}'
    js, xf, part, cf, chk, trj = (pathlib.Path(f'{stem}{e}') for e in
                                  ('.json', '.xyz', '.partial.xyz', '.constraints', '.chk', '_traj.xyz'))
    if js.exists():
        r = json.load(open(js))
        if r.get('status') == 'done' and (r['converged'] or not args.retry_unconverged):
            print(f'[torsion] {stem.name}: done earlier (E = {r["E_Ha"]:.8f}), skipping', flush=True)
            return 0
    x_first = xs
    if part.exists():
        x_first, how = read_xyz(part)[1], f'resumed from {part.name} (first start: {how})'
    elif js.exists() and xf.exists():
        x_first, how = read_xyz(xf)[1], f'retry of unconverged {xf.name}'
    lvl = level_description(case, basis=args.basis)
    print(f'[torsion] {tag}: {lvl}\n[torsion] start: {how}', flush=True)
    tmpxyz = OUT / f'.{stem.name}_mol.xyz'
    write_xyz(tmpxyz, atoms.syms, x_first, 'scanner template')
    gs, mol = make_scanner(case, tmpxyz, basis=args.basis, chkfile=chk)
    tmpxyz.unlink(missing_ok=True)
    lvl['n_basis'] = int(mol.nao)
    deadline = (time.time() + 60 * args.max_min) if args.max_min else None
    res = optimise_torsion(gs, atoms.syms, x_first, dih=dih, cfile=cf, partial=part, traj=trj,
                           maxsteps=(args.maxsteps if kind == 'point' else args.ref_maxsteps),
                           loose=(args.conv == 'loose'), tag=tag, deadline=deadline)
    m = atoms.metrics(res['xyz'])
    r = {'case': args.case, 'species': CASES[args.case], 'kind': 'reference' if kind == 'ref' else 'constrained',
         'target_tau_fold_deg': None if kind == 'ref' else tau,
         'branch': None if kind == 'ref' else args.branch,
         'constraint': None if dih is None else {'dihedral': [atoms.labels[i] for i in dih[:4]],
                                                  'atoms_1based': [i + 1 for i in dih[:4]], 'value_deg': dih[4]},
         'constraint_dev_deg': None if dih is None else abs(wrap(dihedral(res['xyz'], *dih[:4]) - dih[4])),
         'labels': [] if kind == 'ref' else plan_labels(args.case, tau),
         'rigid_rotation_deg': delta,
         'status': 'time_limit' if res['time_limit'] else 'done',
         'converged': res['converged'], 'accepted_on': res['accepted_on'], 'stopped': res['stopped'],
         'steps': res['steps'], 'wall_min': round(res['wall_s'] / 60, 2), 'E_Ha': res['E_Ha'], 'metrics': m,
         'started_from': how, 'level': lvl, 'history': res['history'], **stamp()}
    if js.exists() and part.exists():                      # a resumed run: keep the earlier steps' count
        try:
            old = json.load(open(js))
            r['steps_total'] = old.get('steps_total', old.get('steps', 0)) + r['steps']
            r['wall_min_total'] = round(old.get('wall_min_total', old.get('wall_min', 0)) + r['wall_min'], 2)
        except (OSError, ValueError):
            pass
    write_xyz(xf, atoms.syms, res['xyz'], f"{tag}  E = {res['E_Ha']:.10f} Ha  tau_fold = {m['tau_fold_deg']:.3f}  "
                                          f"(C3-C4-S1-N2 {m['tau_C3_deg']:.3f}, C5-C4-S1-N2 {m['tau_C5_deg']:.3f})")
    atomic_json(js, r)
    if not res['time_limit']:
        part.unlink(missing_ok=True); chk.unlink(missing_ok=True)
    state = ('TIME LIMIT (resumable)' if res['time_limit'] else 'converged' if res['converged'] else 'NOT CONVERGED')
    print(f"[torsion] {stem.name}: {state} in {res['steps']} steps, {res['wall_s'] / 60:.1f} min; "
          f"E = {res['E_Ha']:.8f}; tau_fold {m['tau_fold_deg']:.3f} (C3 {m['tau_C3_deg']:.3f}, C5 "
          f"{m['tau_C5_deg']:.3f}); C4-S1-N2 {m['a_C4_S1_N2_deg']:.2f}", flush=True)
    if not args.out:
        stage_note(case, f"torsion scan ({stem.name}): {state}, {res['steps']} steps, E = {res['E_Ha']:.8f}")
    try:
        report_all(args, [args.case])
    except Exception as err:                                # the report must never fail a finished point
        print(f'[torsion] report failed ({err!r}); the point itself is saved', flush=True)
    # a point stopped by its wall-clock limit is not finished: exit 3 so the burst shows it (rerun resumes it);
    # a point that hit --maxsteps is finished but flagged unconverged in the tables (as dft/scan.py)
    return 3 if res['time_limit'] else 0


# --------------------------------------------------------------------------------------------------- report --
def load_case_points(OUT):
    pts, ref = [], None
    for p in sorted(OUT.glob('*.json')):
        if p.name.startswith('torsion_scan'):
            continue
        try:
            r = json.load(open(p))
        except ValueError:
            continue
        if r.get('kind') == 'reference':
            ref = r
        elif r.get('kind') == 'constrained':
            pts.append(r)
    return pts, ref


def spline(xs, ys):
    from scipy.interpolate import CubicSpline
    o = np.argsort(xs)
    xs, ys = np.asarray(xs, float)[o], np.asarray(ys, float)[o]
    keep = np.concatenate([[True], np.diff(xs) > 1e-3])
    return CubicSpline(xs[keep], ys[keep], bc_type='not-a-knot') if keep.sum() >= 4 else None


def case_report(args, case_name):
    case = load_case(case_name)
    atoms = Atoms(case)
    OUT = pathlib.Path(args.out) / case_name if args.out else case['dir'] / 'torsion_scan'
    if not OUT.exists():
        return None
    pts, ref = load_case_points(OUT)
    x_opt = read_xyz(case['dir'] / 'opt.xyz')[1]
    if ref and ref['status'] == 'done':
        E0, ref_src, ref_m = ref['E_Ha'], f"ref.json (unconstrained re-optimisation, converged={ref['converged']})", \
            ref['metrics']
    else:
        E0 = json.load(open(case['dir'] / 'opt_summary.json'))['E_final_Ha']
        ref_src, ref_m = 'opt_summary.json E_final_Ha (ref task not finished)', atoms.metrics(x_opt)
        if args.basis:
            E0 = None                          # a test at another basis has no fallback reference
    rows = []
    for r in sorted(pts, key=lambda r: (r['branch'], r['target_tau_fold_deg'])):
        m = r['metrics']
        e = None if E0 is None else (r['E_Ha'] - E0) * KCAL_PER_HA
        rows.append({'target_tau_fold_deg': r['target_tau_fold_deg'], 'branch': r['branch'],
                     'tau_fold_deg': m['tau_fold_deg'], 'tau_C3_deg': m['tau_C3_deg'], 'tau_C5_deg': m['tau_C5_deg'],
                     'strain_kcal': e, 'strain_kJ': None if e is None else e * 4.184, 'E_Ha': r['E_Ha'],
                     'status': r['status'], 'converged': r['converged'], 'accepted_on': r['accepted_on'],
                     'steps': r.get('steps_total', r['steps']), 'wall_min': r.get('wall_min_total', r['wall_min']),
                     'rigid_rotation_deg': r['rigid_rotation_deg'], 'constraint_dev_deg': r['constraint_dev_deg'],
                     'a_C4_S1_N2_deg': m['a_C4_S1_N2_deg'], 'd_S1_N2_A': m['d_S1_N2_A'], 'd_C4_S1_A': m['d_C4_S1_A'],
                     'S1_out_of_plane_deg': m['S1_out_of_plane_deg'], 'N1_angle_sum_deg': m['N1_angle_sum_deg'],
                     'N2_angle_sum_deg': m['N2_angle_sum_deg'],
                     'below_reference': bool(e is not None and e < -0.01), 'labels': r.get('labels', [])})
    # the profile: converged C3-branch points + the minimum (x = C3-C4-S1-N2, the constrained coordinate)
    good = [w for w in rows if w['branch'] == 'C3' and w['status'] == 'done' and w['converged']
            and w['strain_kcal'] is not None]
    xs = [w['target_tau_fold_deg'] for w in good] + [ref_m['tau_C3_deg']]
    ys = [w['strain_kcal'] for w in good] + [0.0]
    cs = spline(xs, ys) if E0 is not None else None
    # interpolation check: grid + minimum alone predicts every exact off-grid point
    gx = [w['target_tau_fold_deg'] for w in good if w['target_tau_fold_deg'] in GRID] + [ref_m['tau_C3_deg']]
    gy = [w['strain_kcal'] for w in good if w['target_tau_fold_deg'] in GRID] + [0.0]
    cg = spline(gx, gy) if E0 is not None else None
    chk = [{'tau': w['target_tau_fold_deg'], 'exact_kcal': w['strain_kcal'],
            'grid_spline_kcal': float(cg(w['target_tau_fold_deg']))}
           for w in good if cg is not None and w['target_tau_fold_deg'] not in GRID
           and min(gx) <= w['target_tau_fold_deg'] <= max(gx)]
    for c in chk:
        c['diff_kcal'] = c['grid_spline_kcal'] - c['exact_kcal']
    mirror = None
    b5 = [w for w in rows if w['branch'] == 'C5' and w['status'] == 'done' and w['strain_kcal'] is not None]
    if b5:
        mirror = []
        for w in b5:
            twin = next((v for v in good if abs(v['target_tau_fold_deg'] - w['target_tau_fold_deg']) < EXACT_TOL), None)
            mirror.append({'tau': w['target_tau_fold_deg'], 'C5_branch_kcal': w['strain_kcal'],
                           'C3_branch_kcal': None if twin is None else twin['strain_kcal'],
                           'diff_kcal': None if twin is None else w['strain_kcal'] - twin['strain_kcal']})
    summ = {'case': case_name, 'species': CASES[case_name], 'scan_dir': str(OUT),
            'reference': {'E_Ha': E0, 'source': ref_src, 'tau_fold_deg': ref_m['tau_fold_deg'],
                          'tau_C3_deg': ref_m['tau_C3_deg']},
            'level': (pts[0]['level'] if pts else (ref or {}).get('level')),
            'n_points': len(rows), 'n_converged': sum(bool(w['converged']) for w in rows),
            'n_time_limit': sum(w['status'] == 'time_limit' for w in rows),
            'points': rows, 'interp_check': chk,
            'interp_check_max_abs_kcal': max((abs(c['diff_kcal']) for c in chk), default=None),
            'mirror_branch_check': mirror,
            'note': 'strain = E(constrained, all else relaxed) - E(unconstrained minimum), electronic energies, '
                    'kcal/mol; x of the profile = the constrained C3-C4-S1-N2; tau_fold measured per point.'}
    cols = ['target_tau_fold_deg', 'branch', 'tau_fold_deg', 'tau_C3_deg', 'tau_C5_deg', 'strain_kcal', 'strain_kJ',
            'status', 'converged', 'accepted_on', 'steps', 'wall_min', 'rigid_rotation_deg', 'constraint_dev_deg',
            'a_C4_S1_N2_deg', 'd_S1_N2_A', 'S1_out_of_plane_deg', 'N1_angle_sum_deg', 'N2_angle_sum_deg',
            'below_reference']
    write_csv(OUT / 'torsion_scan.csv', cols, rows)
    atomic_json(OUT / 'torsion_scan_summary.json', summ)
    return {'summary': summ, 'spline': cs, 'good': good, 'ref_m': ref_m}


def fmt(v, nd=2):
    if v is None:
        return ''
    if isinstance(v, float):
        return f'{v:.{nd}f}'
    return str(v)


def write_csv(path, cols, rows):
    tmp = pathlib.Path(str(path) + '.tmp')
    with open(tmp, 'w') as f:
        f.write(','.join(cols) + '\n')
        for w in rows:
            f.write(','.join('"' + str(w.get(c, '')).replace('"', "'") + '"' if isinstance(w.get(c), str)
                             and (',' in w.get(c) or ' ' in w.get(c)) else fmt(w.get(c), 4) for c in cols) + '\n')
    os.replace(tmp, path)


def report_all(args, cases=None):
    cases = cases or list(CASES)
    COMB = pathlib.Path(args.out) / 'combined' if args.out else COMBINED
    COMB.mkdir(parents=True, exist_ok=True)
    with FileLock(COMB / '.report.lock'):
        res = {}
        for c in CASES:                          # always rebuild both, so the combined table is complete
            res[c] = case_report(args, c)
        plan = load_plan()
        table = []
        for c, rr in res.items():
            if rr is None:
                continue
            sp = CASES[c]
            good, cs = rr['good'], rr['spline']
            ex = {round(w['target_tau_fold_deg'], 1): w for w in good}
            req = (plan['cases'][c]['requested'] if plan and c in plan['cases']
                   else [{'tau': t, 'tau_full': t, 'kind': 'grid', 'label': f'grid {t:.0f}'} for t in GRID])
            dstr = {}
            ds = load_case(c)['dir'] / 'dock_strain.json'
            if ds.exists():
                try:
                    dstr = json.load(open(ds))['relative']['b3lyp-d3bj']['vs_first_kcal']
                except (KeyError, ValueError):
                    dstr = {}
            for q in req:
                w = next((v for t, v in ex.items() if abs(t - q['tau']) < EXACT_TOL), None)
                if w is not None:
                    val, how = w['strain_kcal'], 'exact'
                elif cs is not None and min(cs.x) <= q['tau_full'] <= max(cs.x):
                    val, how = float(cs(q['tau_full'])), 'interpolated'
                else:
                    val, how = None, 'pending'
                pose = q.get('pose')
                table.append({'species': SPECIES_LABEL[sp], 'case': c, 'kind': q['kind'], 'label': q['label'],
                              'tau_fold_deg': q['tau_full'], 'strain_kcal': val,
                              'strain_kJ': None if val is None else val * 4.184, 'how': how,
                              'tau_fold_measured_deg': None if w is None else w['tau_fold_deg'],
                              'converged': None if w is None else w['converged'],
                              'mmff_relaxed_kcal': q.get('strain_mmff_relaxed'),
                              'dft_rigid_old_kcal': dstr.get(f'{pose}__Hrot') if pose else None,
                              'plan': f"{q.get('compute', '')}/{q.get('tier', '')}", 'source': q.get('source', '')})
        cols = ['species', 'kind', 'label', 'tau_fold_deg', 'strain_kcal', 'strain_kJ', 'how',
                'tau_fold_measured_deg', 'converged', 'mmff_relaxed_kcal', 'dft_rigid_old_kcal', 'plan', 'source']
        write_csv(COMB / 'torsion_strain.csv', cols, table)
        out = {'what': 'M9 relaxed torsional strain of sulfanilamide and the sulfonamidate anion, B3LYP-D3(BJ)/'
                       '6-311++G(d,p) IEF-PCM water (dft/torsion_scan.py)',
               'definition': 'tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|); strain = E(C3-C4-S1-N2 held at tau, '
                             'everything else relaxed) - E(unconstrained PCM minimum), electronic, kcal/mol',
               'columns': {'how': 'exact = a constrained optimisation at this tau (within 0.06 deg); interpolated = '
                                  'cubic spline through the converged points and the minimum; pending = not yet',
                           'mmff_relaxed_kcal': 'MMFF94s relaxed strain of the same pose (manifest.json)',
                           'dft_rigid_old_kcal': 'burst-2 rigid single point of the same pose, MMFF94s-placed N-H '
                                                 '(dock_strain.json, Table S27 first number): superseded'},
               'generated': datetime.datetime.now().isoformat(timespec='seconds'),
               'per_case': {c: (None if rr is None else {k: rr['summary'][k] for k in
                                ('reference', 'n_points', 'n_converged', 'n_time_limit', 'interp_check',
                                 'interp_check_max_abs_kcal', 'mirror_branch_check', 'level')})
                            for c, rr in res.items()},
               'profiles': {c: (None if rr is None else [{'tau_C3_deg': w['target_tau_fold_deg'],
                                                           'tau_fold_deg': w['tau_fold_deg'],
                                                           'strain_kcal': w['strain_kcal']} for w in
                                                          sorted(rr['good'], key=lambda v: v['target_tau_fold_deg'])])
                            for c, rr in res.items()},
               'table': table}
        atomic_json(COMB / 'torsion_strain.json', out)
        write_md(COMB / 'torsion_strain.md', out)
        try:
            plot(COMB / 'torsion_profile.png', res, table)
        except Exception as err:
            print(f'[torsion] plot skipped ({err!r})', flush=True)
    for c, rr in res.items():
        if rr is None:
            continue
        s = rr['summary']
        print(f"[torsion] {c}: {s['n_converged']}/{s['n_points']} points converged; reference {s['reference']['source']}"
              + ('' if s['interp_check_max_abs_kcal'] is None else
                 f"; grid-spline check max |diff| {s['interp_check_max_abs_kcal']:.3f} kcal/mol"), flush=True)
        for w in sorted(rr['good'], key=lambda v: v['target_tau_fold_deg']):
            print(f"   tau {w['target_tau_fold_deg']:6.1f}  (fold {w['tau_fold_deg']:6.2f})  strain "
                  f"{w['strain_kcal']:6.2f} kcal/mol  {w['steps']} steps", flush=True)
    return res


def write_md(path, out):
    L = [f"# M9 relaxed torsional strain ({out['generated']})", '', out['definition'] + '.', '']
    for c, s in out['per_case'].items():
        if s is None:
            continue
        L.append(f"- {c}: {s['n_converged']}/{s['n_points']} points converged, reference: {s['reference']['source']}"
                 + ('' if s['interp_check_max_abs_kcal'] is None else
                    f"; interpolation check (grid spline vs exact off-grid points) max |diff| "
                    f"{s['interp_check_max_abs_kcal']:.3f} kcal/mol")
                 + ('' if not s['mirror_branch_check'] else
                    '; mirror branch: ' + ', '.join(f"{m['tau']:.0f} deg {fmt(m['diff_kcal'], 3)} kcal/mol"
                                                     for m in s['mirror_branch_check'])))
    L += ['', '| Species | Torsion | tau_fold (deg) | Relaxed strain (kcal/mol) | How | MMFF94s relaxed | '
              'Rigid DFT, burst 2 (superseded) |', '|---|---|---:|---:|---|---:|---:|']
    for r in out['table']:
        L.append(f"| {r['species']} | {r['label']} | {r['tau_fold_deg']:.1f} | {fmt(r['strain_kcal'])} | {r['how']} | "
                 f"{fmt(r['mmff_relaxed_kcal'])} | {fmt(r['dft_rigid_old_kcal'])} |")
    tmp = pathlib.Path(str(path) + '.tmp')
    tmp.write_text('\n'.join(L) + '\n')
    os.replace(tmp, path)


def plot(path, res, table):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for (c, rr), col in zip(res.items(), ('C0', 'C3')):
        if rr is None or rr['spline'] is None:
            continue
        g = sorted(rr['good'], key=lambda v: v['target_tau_fold_deg'])
        xx = np.linspace(min(rr['spline'].x), max(rr['spline'].x), 300)
        ax.plot(xx, rr['spline'](xx), '-', color=col, lw=1, label=f'{SPECIES_LABEL[CASES[c]]} (spline)')
        ax.plot([w['target_tau_fold_deg'] for w in g], [w['strain_kcal'] for w in g], 'o', color=col, mfc='none')
        for q in table:
            if q['case'] == c and q['kind'] == 'co-crystal' and q['strain_kcal'] is not None:
                ax.plot(q['tau_fold_deg'], q['strain_kcal'], 'k*', ms=8)
    for q in {(q['tau_fold_deg'], q['label']) for q in table if q['kind'] == 'co-crystal'}:
        ax.axvline(q[0], color='0.6', lw=0.6, ls=':')
        ax.text(q[0], ax.get_ylim()[1], q[1].split()[0], rotation=90, va='top', ha='right', fontsize=7)
    ax.set_xlabel('tau_fold = C3-C4-S1-N2 held (deg)'); ax.set_ylabel('relaxed strain (kcal/mol)')
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


# ------------------------------------------------------------------------------------------------ check / plan --
def check(args):
    plan = load_plan()
    assert plan, f'{PLAN} missing (run `torsion_scan.py plan --write` on the Pi and sync it)'
    n = 0
    for c, cp in plan['cases'].items():
        case = load_case(c)
        atoms = Atoms(case)
        assert [atoms.C3, atoms.C4, atoms.C5, atoms.N2, atoms.S1] == [2, 3, 4, 7, 10], 'unexpected atom order'
        s, x0 = read_xyz(case['dir'] / 'opt.xyz')
        assert s == atoms.syms, f'{c}: opt.xyz atom order differs from the case'
        E = json.load(open(case['dir'] / 'opt_summary.json'))['E_final_Ha']
        m = atoms.metrics(x0)
        pts = [(p['tau'], 'C3') for p in cp['points']] + ([(150.0, 'C5')] if c == 'sulfanilamide_pcm_d3' else [])
        for tau, br in pts:
            x, d, dmin = rotated_start(atoms, x0, tau, br)
            mm = atoms.metrics(x)
            assert abs(mm['tau_fold_deg'] - tau) < 3.5, (c, tau, mm['tau_fold_deg'])
            n += 1
        print(f"[check] {c}: opt.xyz E = {E:.8f}, tau_fold {m['tau_fold_deg']:.2f} (C3 {m['tau_C3_deg']:.2f}); "
              f"{len(pts)} starts build", flush=True)
    print(f'[check] OK: {n} start geometries', flush=True)


def do_plan(args):
    plan = build_plan()
    for c, cp in plan['cases'].items():
        print(f'== {c}')
        for q in cp['requested']:
            print(f"  {q['tau']:6.1f}  {q['compute']:6s} {q['tier']:8s} {q['kind']:12s} {q['label']}")
        print(f"  -> {len(cp['points'])} exact points: " +
              ', '.join(f"{p['tau']:.1f}{'' if p['tier'] == 'critical' else '*'}" for p in cp['points']))
    if args.write:
        atomic_json(PLAN, plan)
        print(f'[plan] written {PLAN}')
    if args.fragment:
        tasks, crit, opt = fragment(plan)
        tot_c = sum(t['est_min'] for t in tasks if t['id'] in crit)
        tot_o = sum(t['est_min'] for t in tasks if t['id'] in opt)
        tier_min = {k: sum(t['est_min'] for t in tasks if t.get('tier') == k) for k in ('core', 'named_docked')}
        tier_n = {k: sum(t.get('tier') == k for t in tasks) for k in ('core', 'named_docked')}
        frag = {'fragment': 'm9', 'burst': 'burst3b', 'thread_budget_assumed': 40,
                'notes': FRAGMENT_NOTES.format(n_crit=len(crit), n_opt=len(opt), tot_c=tot_c, tot_o=tot_o,
                                               wall_c10=math.ceil(tot_c / 10), wall_c8=math.ceil(tot_c / 8),
                                               mpg=MIN_PER_GRAD, n_core=tier_n['core'], min_core=tier_min['core'],
                                               n_named=tier_n['named_docked'], min_named=tier_min['named_docked']),
                'tasks': tasks}
        pathlib.Path(args.fragment).write_text(json.dumps(frag, indent=1) + '\n')
        print(f'[plan] fragment {args.fragment}: {len(tasks)} tasks, critical {len(crit)} ({tot_c} task-min), '
              f'optional {len(opt)} ({tot_o} task-min)')


FRAGMENT_NOTES = (
    "M9 fragment (referee ref_chem.md M9; ref_dock.md M5/M6), generated by `dft/torsion_scan.py plan --write "
    "--fragment ...` from molecules/sfa/input/torsion_scan_plan.json.  Relaxed torsional strain of sulfanilamide "
    "(sulfanilamide_pcm_d3) and the sulfonamidate anion (sulfonamidate_pcm_d3) at B3LYP-D3(BJ)/6-311++G(d,p) DF, "
    "grid 4, IEF-PCM water Lebedev 29: C3-C4-S1-N2 held (geomeTRIC $set dihedral, exact enforcement) at each "
    "tau, all else relaxed; tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|) as Fig. 11/Table S27; one branch by "
    "mirror symmetry (optional m9_sfa_b150 checks it).  POINTS: grid 90-180 every 10, the co-crystal torsions "
    "(6RL9 132.1, 1AJ0 111.3, 3TZF 105.6, re-derived from dock/T*/ref_*.sdf), the referee-named docked values "
    "(neutral: 6RL9 AD4Zn mean 147.0, 3TZF Vina representative 107.7, and 6RL9 Vina mean 160.36, which lies "
    "0.36 deg from the 160 grid point and is read off the spline; anion: 3TZF Vina mean 115.05 (run at 115.0), "
    "6RL9 Vina mean 155.8) = critical; the other docked job means = optional (a value within 0.5 "
    "deg of an exact point is interpolated instead); every other requested value (representative poses, SAPH "
    "crystal) is read off the spline and the report checks the spline against the exact off-grid points.  "
    "STARTS: independent, each = opt.xyz with the S1 side of C4-S1 rotated rigidly to the target, so every point "
    "can run at once (no chain).  REFERENCE: m9_*_ref (unconstrained re-optimisation of opt.xyz, same driver).  "
    "THREADS: every task 4 (DFT_THREADS 4, BLAS pinned by proton_tools), memory 6 GB each; 10 concurrent = 40.  "
    "TIME GUARD (as m2): m9_gate stamps status/m9_t0 (it must depend on the burst3a waiter); every point stops "
    "cleanly at 2 h 50 min after the gate (status time_limit, resumable by rerunning the same command); optional "
    "points do not start with less than their est_min (>= 30 min) left.  est_min: {mpg} min per gradient at 4 "
    "threads on the shared VM (burst 1a: this molecule optimised at 1.8 min per gradient on 6 threads beside 38 "
    "other threads; burst 2: single points 38-67 s at 12 threads) x gradients per point 9 + 0.11 x |rotation| "
    "(deg) (Pi test at 6-31G*: 13 gradients for a 44.9 deg rotation, 18 for 92.8 deg; +1 for the production basis) x 1.3 margin.  Totals: critical {n_crit} tasks = {tot_c} "
    "task-min (~{wall_c10} min wall at 10 concurrent, ~{wall_c8} at 8), of which tier core (refs, grid, co-crystal) "
    "{n_core} tasks = {min_core} task-min and tier named_docked {n_named} tasks = {min_named} task-min (droppable "
    "first: the report interpolates them from core and checks the spline against the exact co-crystal points); "
    "optional {n_opt} tasks = {tot_o} task-min.  Every point is independent, so any subset can be run and a "
    "dropped point can be added in a later burst with the same command.  m9_report (0 threads) rebuilds molecules/sfa/results/torsion_scan/torsion_strain.{{json,csv,md}} "
    "and torsion_profile.png; every point task refreshes them too.  Outputs: molecules/sfa/results/"
    "{{sulfanilamide,sulfonamidate}}_pcm_d3/torsion_scan/ (per-point json/xyz/traj, ref.json, torsion_scan.csv, "
    "torsion_scan_summary.json).  Nothing existing is overwritten.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('mode', choices=['plan', 'check', 'ref', 'point', 'report'])
    ap.add_argument('case', nargs='*')
    ap.add_argument('--tau', type=float)
    ap.add_argument('--branch', choices=['C3', 'C5'], default='C3')
    ap.add_argument('--max-min', type=float, help='wall-clock limit of this run in minutes (resumable stop)')
    ap.add_argument('--maxsteps', type=int, default=60)
    ap.add_argument('--ref-maxsteps', type=int, default=40)
    ap.add_argument('--retry-unconverged', action='store_true')
    ap.add_argument('--write', action='store_true', help='plan: write molecules/sfa/input/torsion_scan_plan.json')
    ap.add_argument('--fragment', help='plan: also write the burst fragment here')
    ap.add_argument('--out'); ap.add_argument('--basis')
    ap.add_argument('--conv', choices=['default', 'loose'], default='default')
    args = ap.parse_args()
    if args.basis and not args.out:
        ap.error('--basis is for tests only: give --out as well')
    if args.mode == 'plan':
        return do_plan(args)
    if args.mode == 'check':
        return check(args)
    if args.mode == 'report':
        report_all(args, args.case or None)
        return 0
    if len(args.case) != 1:
        ap.error(f'{args.mode} takes one case')
    args.case = args.case[0]
    if args.mode == 'point' and args.tau is None:
        ap.error('point needs --tau')
    return run_one(args, args.mode)


if __name__ == '__main__':
    sys.exit(main() or 0)
