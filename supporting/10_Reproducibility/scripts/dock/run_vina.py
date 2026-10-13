"""
Vina driver: one (target, ligand, scoring function) job per call, ten independent seeds.

    python dock/run_vina.py --target T3_6RL9 --ligand sfa_neutral
    python dock/run_vina.py --target T1_3TZF --ligand smx_anion --scoring vinardo --seeds 5
    python dock/run_vina.py --target T3_6RL9 --ligand sfa_anion --from-crystal   # control only

Why it is built this way:

*   **Ten seeds, each with its own starting conformer.**  A single Vina run reports one number and
    hides the variance of a stochastic search.  Ten seeds give a mean +- SD for the score and, for
    redocking, a success *fraction* rather than a single lucky pose.  Every seed also gets a fresh
    ETKDG conformer generated from that seed, so the internal geometry the search starts from is
    randomised too: Vina holds bond lengths and angles fixed, so starting all seeds from the crystal
    conformer would quietly hand the answer to the search.  `--from-crystal` exists only to quantify
    that effect, never to produce a reported redocking number.

*   **Politeness on the Pi.**  The Pi has 4 cores and runs other projects' cron jobs, so the default
    is `--cpu 3` and `nice`, and the jobs are run one after another rather than in parallel.

*   **Everything on disk.**  Each seed keeps its input conformer, its Vina log and all poses, so any
    number in the paper can be traced back to a file and recomputed.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from prep_ligands import SPECIES, input_mol_for_seed, write_pdbqt  # noqa: E402

DOCK = pathlib.Path(__file__).resolve().parent
VINA = 'vina'


def conformer_for_seed(species, seed, out_pdbqt, crystal_sdf=None, geometry='etkdg',
                       xyz_override=None):
    """A fresh 3D conformer for one seed, or the crystal conformer when --from-crystal is used, or
    the DFT geometry (same for every seed) with --geometry dft (production 2026-09-28)."""
    if crystal_sdf is not None:
        mol = Chem.AddHs(Chem.MolFromMolFile(str(crystal_sdf)), addCoords=True)
        write_pdbqt(mol, out_pdbqt, f'{species}_crystal')
        return mol, dict(source='crystal', sdf=str(crystal_sdf))
    mol3d, prov = input_mol_for_seed(species, seed, geometry, xyz_override)
    write_pdbqt(mol3d, out_pdbqt, species)
    return mol3d, prov


def parse_vina_out(path):
    """Scores and pose coordinates from a Vina output PDBQT (model by model)."""
    models, cur, score = [], None, None
    for line in path.read_text().splitlines():
        if line.startswith('MODEL'):
            cur, score = [], None
        elif line.startswith('REMARK VINA RESULT'):
            score = float(line.split()[3])
        elif line.startswith(('ATOM', 'HETATM')):
            if line[77:79].strip() != 'HD' and not line[12:16].strip().startswith('H'):
                cur.append((line[12:16].strip(), float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
        elif line.startswith('ENDMDL'):
            models.append(dict(score=score, atoms=cur))
    return models


def main():
    ap = argparse.ArgumentParser()
    # The choices are the prepared receptor *directories*, not the five target keys, because every
    # receptor variant (T1_3TZF_seed12345, T2_1AJ0_noSO4, T3_6RL9_his64HID, ...) is docked exactly like
    # a target.  A directory only counts once it holds both files a run needs, so a half-prepared or
    # leftover directory cannot be selected by mistake.  'C' directories are the control receptor
    # (C1_4DUH and its variants, 2026-10-04): docked by this same driver, analysed by control_gyrb.py.
    ap.add_argument('--target', required=True, choices=sorted(
        d.name for d in DOCK.iterdir()
        if d.is_dir() and d.name.startswith(('T', 'C')) and (d / 'box.json').exists()
        and (d / 'receptor.pdbqt').exists()))
    ap.add_argument('--ligand', required=True, choices=sorted(SPECIES))
    ap.add_argument('--scoring', default='vina', choices=['vina', 'vinardo'])
    ap.add_argument('--seeds', type=int, default=10)
    ap.add_argument('--exhaustiveness', type=int, default=32)
    ap.add_argument('--num-modes', type=int, default=20)
    ap.add_argument('--energy-range', type=float, default=3.0)
    ap.add_argument('--cpu', type=int, default=3)
    ap.add_argument('--box-size', type=float, default=None,
                    help='override the 20 A box edge (diagnostic: a 20 A box on a 280-residue '
                         'protein can contain a competing pocket, which is what happened in 1AJ0)')
    ap.add_argument('--box-json', default=None,
                    help='take the box centre and size from this JSON file (keys center, size) instead '
                         'of the target\'s box.json; used for the whole-fragment search of the control '
                         'receptor, whose box is not centred on the crystal ligand')
    ap.add_argument('--from-crystal', action='store_true',
                    help='start every seed from the crystal conformer (control, not a result)')
    ap.add_argument('--geometry', default='etkdg', choices=['etkdg', 'dft'],
                    help='etkdg: a fresh ETKDG conformer per seed (validation); dft: the DFT '
                         'geometry in dock/ligands/<species>.sdf for every seed (production)')
    ap.add_argument('--geometry-xyz', default=None,
                    help='dock this xyz instead of the species geometry in dock/ligands/ '
                         '(diagnostic; implies --geometry dft)')
    ap.add_argument('--tag-suffix', default='',
                    help='appended to the run directory name, e.g. "__repro", so that a repeat or a '
                         'variant never overwrites the run it is compared with')
    a = ap.parse_args()

    tdir = DOCK / a.target
    box = json.loads((tdir / 'box.json').read_text())
    if a.box_json:
        other = json.loads(pathlib.Path(a.box_json).read_text())
        box = dict(other, replaces=dict(center=box['center'], size=box['size']),
                   box_json=str(pathlib.Path(a.box_json).resolve().relative_to(DOCK)))
    if a.box_size:
        box = dict(box, size=[a.box_size] * 3, size_overridden_from=box['size'])
    receptor = tdir / 'receptor.pdbqt'
    if not receptor.exists():
        raise SystemExit(f'{receptor} missing; run prep_receptor.py first')
    tag = f'{a.ligand}__{a.scoring}' + ('__fromxtal' if a.from_crystal else '') \
        + (f'__box{a.box_size:g}' if a.box_size else '') + a.tag_suffix
    rundir = tdir / 'runs' / tag
    rundir.mkdir(parents=True, exist_ok=True)
    crystal = None
    if a.from_crystal:
        crystal = next(tdir.glob('ref_*.sdf'))

    results = dict(target=a.target, ligand=a.ligand, scoring=a.scoring, box=box, tag=tag,
                   exhaustiveness=a.exhaustiveness, num_modes=a.num_modes,
                   energy_range=a.energy_range, cpu=a.cpu, seeds=[],
                   from_crystal=bool(a.from_crystal),
                   geometry='dft' if a.geometry_xyz else a.geometry,
                   receptor=str(receptor.relative_to(DOCK)),
                   vina_version=subprocess.run([VINA, '--version'], capture_output=True,
                                               text=True).stdout.strip())
    for seed in range(1, a.seeds + 1):
        sdir = rundir / f'seed{seed:02d}'
        sdir.mkdir(exist_ok=True)
        lig = sdir / 'ligand.pdbqt'
        _, prov = conformer_for_seed(a.ligand, seed, lig, crystal_sdf=crystal,
                                     geometry='dft' if a.geometry_xyz else a.geometry,
                                     xyz_override=a.geometry_xyz)
        results.setdefault('input_geometry', prov)
        out = sdir / 'out.pdbqt'
        cmd = [VINA, '--receptor', str(receptor), '--ligand', str(lig),
               '--center_x', str(box['center'][0]), '--center_y', str(box['center'][1]),
               '--center_z', str(box['center'][2]),
               '--size_x', str(box['size'][0]), '--size_y', str(box['size'][1]),
               '--size_z', str(box['size'][2]),
               '--exhaustiveness', str(a.exhaustiveness), '--num_modes', str(a.num_modes),
               '--energy_range', str(a.energy_range), '--seed', str(seed),
               '--cpu', str(a.cpu), '--out', str(out)]
        if a.scoring != 'vina':
            cmd += ['--scoring', a.scoring]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.time() - t0
        (sdir / 'vina.log').write_text(r.stdout + '\n' + r.stderr)
        if r.returncode != 0 or not out.exists():
            results['seeds'].append(dict(seed=seed, ok=False, wall_s=round(wall, 1),
                                         stderr=r.stderr[-1500:]))
            print(f'  seed {seed}: FAILED ({r.stderr.strip().splitlines()[-1:]})')
            continue
        models = parse_vina_out(out)
        scores = [m['score'] for m in models]
        results['seeds'].append(dict(seed=seed, ok=True, wall_s=round(wall, 1),
                                     n_modes=len(models), scores=scores,
                                     best=scores[0] if scores else None,
                                     out=str(out.relative_to(DOCK))))
        print(f'  seed {seed:2d}: best {scores[0]:7.3f} kcal/mol, {len(models)} modes, {wall:.0f} s')
    best = [s['best'] for s in results['seeds'] if s.get('ok')]
    if best:
        results['summary'] = dict(n_ok=len(best), best_mean=round(float(np.mean(best)), 3),
                                  best_sd=round(float(np.std(best, ddof=1)), 3) if len(best) > 1
                                  else 0.0,
                                  best_min=round(float(np.min(best)), 3),
                                  best_max=round(float(np.max(best)), 3),
                                  wall_s_total=round(sum(s['wall_s'] for s in results['seeds']), 1))
    (rundir / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    if best:
        s = results['summary']
        print(f"{a.target} {tag}: {s['best_mean']:.3f} +- {s['best_sd']:.3f} kcal/mol "
              f"(n={s['n_ok']}), total {s['wall_s_total']:.0f} s")
    return 0


if __name__ == '__main__':
    sys.exit(main())
