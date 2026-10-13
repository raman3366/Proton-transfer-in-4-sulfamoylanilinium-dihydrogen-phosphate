"""
Step 1 -- geometry optimisation at B3LYP/6-311++G(d,p).
Usage: python dft/optimize.py <case> [--restart]
Writes results/<case>/opt.xyz, opt_traj.xyz, opt_summary.json, geometry_vs_start.csv

Two-stage protocol (PREOPT_BASIS, default 6-31G(d)): the deposited structures have poor hydrogen positions
-- PZC's hydrogens were never refined and include a C-H of 1.14 A -- so most optimisation steps are spent
relaxing X-H bonds.  Doing that at a small basis first is ~15x cheaper per step and leaves only a short
refinement at the requested level.  The reported geometry is always converged at B3LYP/6-311++G(d,p) with
Gaussian's default criteria; the pre-optimisation only supplies a better starting point.
Set PREOPT_BASIS=none to disable.
"""
import sys, os, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from common import *
from pyscf.geomopt.geometric_solver import optimize

case = load_case(sys.argv[1]); d = case['dir']
start = d / 'opt.xyz' if ('--restart' in sys.argv and (d / 'opt.xyz').exists()) else ROOT / case['xyz']
print(f"[opt] case={case['name']} start geometry={start}")

# ---- stage 0: cheap pre-optimisation to fix the crystallographic hydrogen positions ----
PREOPT = os.environ.get('PREOPT_BASIS', '6-31g(d)')
pre_steps = 0
# Loose criteria for the pre-optimisation (geomeTRIC's GAU_LOOSE).  These clusters have very flat
# intermolecular modes -- PZC reached an rms gradient of 4.7e-5, six times below the target, and still kept
# failing the displacement criterion while the trust radius collapsed, crawling along a rotating-water mode.
# The pre-optimisation only has to supply a starting geometry, so chasing that is pure waste; the main
# optimisation below keeps Gaussian's default criteria unchanged.
PRE_CONV = dict(convergence_energy=1e-6, convergence_grms=1.7e-3, convergence_gmax=2.5e-3,
                convergence_drms=6.7e-3, convergence_dmax=1.0e-2)
if (d / 'preopt.xyz').exists() and '--restart' not in sys.argv:
    start = d / 'preopt.xyz'
    print(f'[opt] reusing existing {start.name}; skipping the pre-optimisation', flush=True)
    stage_note(case, 'reused existing preopt.xyz (pre-optimisation skipped)')
elif PREOPT.lower() not in ('none', '', '0') and '--restart' not in sys.argv:
    t_pre = time.time()
    pre_mol = build_mol(start, case['charge'], case['spin'], basis=PREOPT,
                        basis_overrides=case.get('basis_overrides'))
    pre_mf = make_mf(pre_mol, conv_tol=1e-8, solvent=case.get('solvent'))
    print(f'[opt] pre-optimisation at B3LYP/{PREOPT} ({pre_mol.nao} basis functions) ...', flush=True)
    stage_note(case, f'pre-optimisation at B3LYP/{PREOPT} ({pre_mol.nao} basis functions)')
    pre_hist = []
    pre_opt = optimize(pre_mf, maxsteps=300,
                       callback=lambda env: pre_hist.append(float(env['energy'])), **PRE_CONV)
    pre_steps = len(pre_hist)
    start = d / 'preopt.xyz'
    write_xyz(start, [a[0] for a in pre_opt._atom], pre_opt.atom_coords() * BOHR_A,
              f'pre-optimised at B3LYP/{PREOPT}')
    print(f'[opt] pre-optimisation done: {pre_steps} steps, {(time.time()-t_pre)/60:.1f} min -> {start.name}', flush=True)
    stage_note(case, f'pre-optimisation converged in {pre_steps} steps, {(time.time()-t_pre)/60:.1f} min')

mol = mol_for_case(case, start)
print(f'[opt] main optimisation at B3LYP/{basis_description(case)} ({mol.nao} basis functions)', flush=True)
mf = make_mf(mol, conv_tol=1e-9 if not case.get('solvent') else 1e-8, chkfile=d / 'opt.chk', solvent=case.get('solvent'))
stage_note(case, f"optimisation started from {start.name}; solvent={case.get('solvent') or 'none (gas phase)'}")

traj = open(d / 'opt_traj.xyz', 'w'); hist = []
def cb(env):
    m = env['mol']; e = env['energy']; g = env['gradients']
    hist.append({'step': len(hist), 'E_Ha': float(e), 'gmax': float(abs(g).max()), 'grms': float(np.sqrt((g**2).mean()))})
    syms = [a[0] for a in m._atom]; xyz = m.atom_coords() * BOHR_A
    traj.write(f'{len(syms)}\nstep {len(hist)-1}  E = {e:.8f} Ha\n' + ''.join(f'{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n' for s, (x, y, z) in zip(syms, xyz))); traj.flush()
    print(f"[opt] step {len(hist)-1:3d}  E = {e:.8f} Ha  max|g| = {abs(g).max():.2e}  rms g = {np.sqrt((g**2).mean()):.2e}", flush=True)

t0 = time.time()
# geomeTRIC defaults == Gaussian's default convergence criteria (max force 4.5e-4, RMS force 3e-4 Ha/Bohr, ...)
mol_opt = optimize(mf, maxsteps=300, callback=cb)
wall = time.time() - t0
syms = [a[0] for a in mol_opt._atom]; xyz = mol_opt.atom_coords() * BOHR_A
write_xyz(d / 'opt.xyz', syms, xyz, f"{case['title']} -- B3LYP/{basis_description(case)} optimised, E = {hist[-1]['E_Ha']:.8f} Ha")

# geometry comparison against the starting (X-ray) structure
s0, x0 = read_xyz(ROOT / case['xyz']); lab = case.get('cif_labels') or atom_labels(syms)
rows = []
for i, j, d0 in bonds_from_geometry(s0, x0):
    d1 = np.linalg.norm(xyz[i] - xyz[j]); rows.append((f'{lab[i]}-{lab[j]}', d0, d1))
with open(d / 'geometry_vs_start.csv', 'w') as f:
    f.write('bond,start_A,dft_A,diff_A\n')
    for n, a, b in rows: f.write(f'{n},{a:.4f},{b:.4f},{b-a:+.4f}\n')
diffs = np.array([b - a for _, a, b in rows])
summary = {'case': case['name'], 'E_final_Ha': hist[-1]['E_Ha'], 'steps': len(hist), 'wall_min': wall / 60,
           'preopt_basis': (PREOPT if pre_steps else None), 'preopt_steps': pre_steps, 'n_basis': int(mol.nao),
           'bond_rmsd_A': float(np.sqrt((diffs**2).mean())), 'bond_maxdev_A': float(abs(diffs).max()),
           'history': hist}
# extra distances requested in the case file (e.g. hydrogen-bond N...O distances)
for name, (i, j) in case.get('distances', {}).items():
    summary.setdefault('distances_A', {})[name] = {'start': float(np.linalg.norm(x0[i]-x0[j])), 'dft': float(np.linalg.norm(xyz[i]-xyz[j]))}
for name, (D, H, A) in case.get('hbonds', {}).items():
    summary.setdefault('hbonds', {})[name] = {'start': hbond_geometry(x0, D, H, A), 'dft': hbond_geometry(xyz, D, H, A)}
save_json(d / 'opt_summary.json', summary)
stage_note(case, f"optimisation converged in {len(hist)} steps, {wall/60:.1f} min, E = {hist[-1]['E_Ha']:.8f} Ha")
print(f"[opt] done: {len(hist)} steps, {wall/60:.1f} min, E = {hist[-1]['E_Ha']:.8f} Ha, bond RMSD vs start {summary['bond_rmsd_A']:.3f} A")
