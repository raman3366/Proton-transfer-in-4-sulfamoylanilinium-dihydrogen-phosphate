"""
Accept an optimisation that is converged on forces but still chasing the displacement criterion.

These hydrogen-bonded clusters have nearly-free inter-fragment rotations.  geomeTRIC can sit on such a mode
for hours, gaining under a microhartree per step, while the forces are already an order of magnitude inside
Gaussian's threshold (BTP's PCM model reached rms 1.8e-05 against a 3.0e-04 target).  Gaussian itself has an
escape for this -- "optimization completed on the basis of negligible forces" -- and this is the equivalent.

The structure is still proved to be a minimum, by the frequency calculation that follows: zero imaginary
modes is the real test, not the optimiser's displacement bookkeeping.

Usage: python dft/finalize_opt.py <case>
"""
import sys, json, time, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *

case = load_case(sys.argv[1]); d = case['dir']
# Sort by MODIFICATION TIME, not by path.  A case that was restarted leaves several scratch directories
# with random names; sorting by name picked btp_pcm's pre-OOM trajectory (10 frames) over the live one
# (16 frames), which would have fed a stale geometry into a 14-hour Hessian.
traj = sorted(list((ROOT / 'scratch').glob(f'{case["name"]}/*/*_optim.xyz')) +
              list(pathlib.Path('/mnt/scratch').glob(f'{case["name"]}/*/*_optim.xyz')),
              key=lambda p: p.stat().st_mtime)
if not traj:
    raise SystemExit(f'no geomeTRIC trajectory found for {case["name"]}')
T = traj[-1]
lines = T.read_text().splitlines()
n = int(lines[0].split()[0])
frames = len(lines) // (n + 2)
last = lines[(frames - 1) * (n + 2):frames * (n + 2)]
syms = [l.split()[0] for l in last[2:2 + n]]
xyz = np.array([[float(v) for v in l.split()[1:4]] for l in last[2:2 + n]])

# forces from the optimiser log, to record what we accepted on
import re
log = (d / 'logs' / 'optimize.log').read_text(errors='replace')
g = re.findall(r'Grad = ([\d.e+-]+)/([\d.e+-]+)', re.sub(r'\x1b\[[0-9;]*m', '', log))
grms, gmax = (float(g[-1][0]), float(g[-1][1])) if g else (float('nan'),) * 2
E = float(re.findall(r'Energy = (-?\d+\.\d+)', last[1])[0]) if re.findall(r'Energy = (-?\d+\.\d+)', last[1]) \
    else float(re.findall(r'-?\d+\.\d+', last[1])[-1])

write_xyz(d / 'opt.xyz', syms, xyz,
          f"{case['title']} -- B3LYP/{basis_description(case)} grid {GRID_LEVEL}, accepted on negligible forces "
          f"(rms {grms:.2e}, max {gmax:.2e}), E = {E:.8f} Ha")

s0, x0 = read_xyz(ROOT / case['xyz']); lab = case.get('cif_labels') or atom_labels(syms)
rows = [(f'{lab[i]}-{lab[j]}', d0, float(np.linalg.norm(xyz[i] - xyz[j])))
        for i, j, d0 in bonds_from_geometry(s0, x0)]
with open(d / 'geometry_vs_start.csv', 'w') as f:
    f.write('bond,start_A,dft_A,diff_A\n')
    for nm, a, b in rows: f.write(f'{nm},{a:.4f},{b:.4f},{b-a:+.4f}\n')
diffs = np.array([b - a for _, a, b in rows])
summary = {'case': case['name'], 'E_final_Ha': E, 'steps': frames - 1, 'wall_min': float('nan'),
           'n_basis': None, 'bond_rmsd_A': float(np.sqrt((diffs ** 2).mean())),
           'bond_maxdev_A': float(abs(diffs).max()),
           'converged_on': f'negligible forces (rms {grms:.2e}, max {gmax:.2e}; '
                           f'Gaussian thresholds 3.0e-4 / 4.5e-4)'}
for nm, (i, j) in case.get('distances', {}).items():
    summary.setdefault('distances_A', {})[nm] = {'start': float(np.linalg.norm(x0[i]-x0[j])),
                                                 'dft': float(np.linalg.norm(xyz[i]-xyz[j]))}
for nm, (D, H, A) in case.get('hbonds', {}).items():
    summary.setdefault('hbonds', {})[nm] = {'start': hbond_geometry(x0, D, H, A), 'dft': hbond_geometry(xyz, D, H, A)}
save_json(d / 'opt_summary.json', summary)
stage_note(case, f'optimisation accepted on negligible forces (rms {grms:.2e}, max {gmax:.2e}) '
                 f'after {frames-1} steps; minimum to be confirmed by frequencies')
print(f'[finalize] {case["name"]}: frame {frames-1}, E = {E:.8f} Ha, rms grad {grms:.2e}, max {gmax:.2e}')
print(f'[finalize] wrote opt.xyz and opt_summary.json -- run the pipeline from proton_check')
