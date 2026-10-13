"""
AutoDock4Zn docking for the two hCA II targets (Santos-Martins 2014, 10.1021/ci500209e).

    python dock/run_ad4zn.py --target T3_6RL9 --stage maps           # TZ receptor, GPF, maps
    python dock/run_ad4zn.py --target T3_6RL9 --ligand sfa_neutral   # dock (maps must exist)
    python dock/run_ad4zn.py --target T3_6RL9 --stage maps --vm      # make the maps on the VM

Why AD4Zn at all.  The Vina scoring function has no metal-coordination term, so the redocking of
sulfanilamide into 6RL9 with plain Vina put the sulfonamide nitrogen 2.83 A from the zinc where the
crystal has 1.93 A, and lost the Thr199 OG1 hydrogen bond.  AD4Zn adds a tetrahedral zinc pseudo atom
(type TZ) marking the free coordination position and a set of pairwise potentials, and was calibrated
on 292 zinc complexes.  A check worth recording: for 6RL9 the TZ atom placed by the official
zinc_pseudo.py falls 0.17 A from the crystallographic sulfonamide nitrogen and 2.00 A from the zinc,
so the pseudo atom really does mark where the ligand nitrogen belongs.

Why the maps are made elsewhere.  The AD4Zn potentials reach autogrid through `nbp_r_eps` records, and
the AutoDock-Vina documentation warns that this keyword changed behaviour between autogrid 4.2.6 and
4.2.7.  conda-forge and bioconda have no linux-aarch64 autogrid at all and their linux-64 build is
4.2.6, the wrong one, so the Pi cannot make these maps correctly.  The IBM VM has autogrid 4.2.9
(reporting "AutoGrid 4.2.7.x") in ~/cqc-sfa/tools/dock, so `--vm` ships the receptor and GPF to
~/sandbox/docking there, runs autogrid single-threaded, and brings the maps back.  Docking itself then
runs on the Pi from the maps: `vina --scoring ad4 --maps <prefix>` needs no autogrid.

The GPF is written here rather than by prepare_gpf4zn.py, which needs MGLTools and Python 2.  The
keyword order and the six nbp_r_eps records follow the documented example and the upstream script
byte for byte (dock/vendor/prepare_gpf4zn.py lines 157-162).
"""
import argparse
import json
import pathlib
import shlex
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from prep_ligands import SPECIES  # noqa: E402

DOCK = pathlib.Path(__file__).resolve().parent
VENDOR = DOCK / 'vendor'
SPACING = 0.375                  # A, the AutoDock default
VM = 'cqc-dft'
VM_SANDBOX = '~/sandbox/docking'
VM_AUTOGRID = '~/cqc-sfa/tools/dock/bin/autogrid4'

# The AutoDock4Zn pairwise potentials, copied verbatim from the upstream prepare_gpf4zn.py.
NBP_R_EPS = [
    'nbp_r_eps 0.25 23.2135 12 6 NA TZ',
    'nbp_r_eps 2.1   3.8453 12 6 OA Zn',
    'nbp_r_eps 2.25  7.5914 12 6 SA Zn',
    'nbp_r_eps 1.0   0.0    12 6 HD Zn',
    'nbp_r_eps 2.0   0.0060 12 6 NA Zn',
    'nbp_r_eps 2.0   0.2966 12 6  N Zn',
]


def pdbqt_types(path):
    return [ln[77:79].strip() for ln in path.read_text().splitlines() if ln.startswith('ATOM')]


def add_tz(receptor, out):
    """Place the tetrahedral zinc pseudo atoms with the official script (pure Python 3)."""
    r = subprocess.run([sys.executable, str(VENDOR / 'zinc_pseudo.py'), '-r', str(receptor),
                        '-o', str(out)], capture_output=True, text=True)
    if r.returncode != 0 or not out.exists():
        raise SystemExit(f'zinc_pseudo.py failed: {r.stdout}\n{r.stderr}')
    n = sum(1 for t in pdbqt_types(out) if t == 'TZ')
    return n, r.stdout.strip()


def write_gpf(gpf, prefix, receptor_tz, box, ligand_types, receptor_types):
    """The grid parameter file.

    npts must be even, so the requested 20.0 A box becomes 54 x 0.375 = 20.25 A: the AD4Zn search
    volume is 0.25 A wider per axis (0.125 A per face) than the Vina/Vinardo box for the same target.
    Both numbers are recorded in maps_report.json, and the difference is a confound whenever AD4Zn and
    Vina results are compared, so it was measured rather than argued away: redocking sulfanilamide into
    6RL9 with Vina in a 20.25 A box instead of 20.0 A is the control
    dock/T3_6RL9/runs/sfa_neutral__vina__box20.25, and dock/VALIDATION.md reports what it changed.
    Matching the two exactly would mean either an odd npts (autogrid refuses) or a 20.25 A Vina box for
    every target, which would invalidate the existing Vina runs for a 1 % volume difference."""
    npts = [int(2 * round(box['size'][i] / SPACING / 2)) for i in range(3)]
    lines = [f'npts {npts[0]} {npts[1]} {npts[2]}',
             'parameter_file AD4Zn.dat',
             f'gridfld {prefix}.maps.fld',
             f'spacing {SPACING}',
             'receptor_types ' + ' '.join(receptor_types),
             'ligand_types ' + ' '.join(ligand_types),
             f'receptor {receptor_tz.name}',
             'gridcenter ' + ' '.join(f'{v:.3f}' for v in box['center']),
             'smooth 0.5']
    lines += [f'map {prefix}.{t}.map' for t in ligand_types]
    lines += [f'elecmap {prefix}.e.map', f'dsolvmap {prefix}.d.map', 'dielectric -0.1465']
    lines += NBP_R_EPS
    gpf.write_text('\n'.join(lines) + '\n')
    return npts, [round(npts[i] * SPACING, 3) for i in range(3)]


def run_autogrid_local(workdir, gpf, glg):
    from shutil import which
    exe = which('autogrid4')
    if not exe:
        return None, 'autogrid4 not on PATH'
    r = subprocess.run(['nice', '-n', '10', exe, '-p', gpf.name, '-l', glg.name],
                       cwd=workdir, capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr)[-3000:]


def run_autogrid_vm(workdir, gpf, glg, prefix, report, remote_name=None):
    """Ship receptor + GPF + AD4Zn.dat to the VM sandbox, run autogrid there, bring the maps back.

    Single-threaded and confined to ~/sandbox/docking, so it does not disturb the DFT burst that owns
    ~/cqc-sfa.  Nothing is written outside the sandbox."""
    # one remote directory per target: workdir.name is always "ad4zn", so two targets prepared in
    # the same session would otherwise overwrite each other's maps on the VM
    remote = f'{VM_SANDBOX}/{remote_name or workdir.parent.name}'
    files = [workdir / gpf.name, workdir / f'{prefix}.pdbqt', VENDOR / 'AD4Zn.dat']
    cmds = [
        ['ssh', VM, f'mkdir -p {remote} && rm -f {remote}/*.map {remote}/*.fld'],
        ['scp', '-q', *[str(f) for f in files], f'{VM}:{remote}/'],
        ['ssh', VM, f'cd {remote} && nice -n 10 {VM_AUTOGRID} -p {gpf.name} -l {glg.name} '
                    f'&& ls -1 *.map *.fld | wc -l'],
    ]
    log = []
    for c in cmds:
        r = subprocess.run(c, capture_output=True, text=True, timeout=1800)
        log.append(dict(cmd=' '.join(shlex.quote(x) for x in c), rc=r.returncode,
                        out=r.stdout[-1500:], err=r.stderr[-1500:]))
        if r.returncode != 0:
            report['autogrid_vm'] = log
            raise SystemExit(f'VM autogrid step failed: {c}\n{r.stdout}\n{r.stderr}')
    r = subprocess.run(['bash', '-lc',
                        f'scp -q "{VM}:{remote}/*.map" "{VM}:{remote}/*.fld" '
                        f'"{VM}:{remote}/{glg.name}" {shlex.quote(str(workdir))}/'],
                       capture_output=True, text=True, timeout=1800)
    log.append(dict(cmd='scp maps back', rc=r.returncode, err=r.stderr[-1500:]))
    report['autogrid_vm'] = log
    if r.returncode != 0:
        raise SystemExit(f'could not copy the maps back: {r.stderr}')
    return len(list(workdir.glob('*.map')))


def stage_maps(target, use_vm, extra_types):
    tdir = DOCK / target
    box = json.loads((tdir / 'box.json').read_text())
    work = tdir / 'ad4zn'
    work.mkdir(exist_ok=True)
    prefix = 'receptor_tz'
    rec_tz = work / f'{prefix}.pdbqt'
    n_tz, msg = add_tz(tdir / 'receptor.pdbqt', rec_tz)
    if n_tz == 0:
        raise SystemExit(f'{target}: zinc_pseudo.py placed no TZ atoms -- is there a zinc?')
    rtypes = sorted(set(pdbqt_types(rec_tz)))
    # every ligand type we will ever dock against these maps must have a map
    ltypes = set(extra_types)
    for name in SPECIES:
        p = DOCK / 'ligands' / f'{name}.pdbqt'
        if p.exists():
            ltypes |= set(pdbqt_types(p))
    ltypes = sorted(ltypes)
    gpf, glg = work / f'{prefix}.gpf', work / f'{prefix}.glg'
    npts, size = write_gpf(gpf, prefix, rec_tz, box, ltypes, rtypes)
    report = dict(target=target, n_tz_atoms=n_tz, zinc_pseudo_message=msg,
                  receptor_types=rtypes, ligand_types=ltypes, npts=npts, spacing=SPACING,
                  box_size_from_npts=size, box_size_requested=box['size'],
                  gridcenter=box['center'], nbp_r_eps=NBP_R_EPS)
    # a check with real content: the TZ atom should mark where the crystal ligand's sulfonamide N is
    try:
        from rdkit import Chem
        ref = Chem.MolFromMolFile(str(next(tdir.glob('ref_*.sdf'))))
        conf = ref.GetConformer()
        tz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                       for ln in rec_tz.read_text().splitlines()
                       if ln.startswith('ATOM') and ln[77:79].strip() == 'TZ'])
        ns = [np.array(list(conf.GetAtomPosition(a.GetIdx()))) for a in ref.GetAtoms()
              if a.GetSymbol() == 'N' and any(n.GetSymbol() == 'S' for n in a.GetNeighbors())]
        if len(ns):
            report['tz_to_crystal_sulfonamide_N_A'] = round(
                float(min(np.linalg.norm(a - b) for a in tz for b in ns)), 2)
    except Exception as exc:
        report['tz_check_error'] = f'{type(exc).__name__}: {exc}'

    if not use_vm:
        rc, log = run_autogrid_local(work, gpf, glg)
        report['autogrid_local'] = dict(returncode=rc, log=log)
        if rc != 0:
            report['maps_ready'] = False
            report['note'] = ('no usable local autogrid (>=4.2.7 required for nbp_r_eps); '
                              'rerun with --vm, or run autogrid4 by hand in ' + str(work))
    else:
        (VENDOR / 'AD4Zn.dat').exists() or sys.exit('vendor/AD4Zn.dat missing')
        n_maps = run_autogrid_vm(work, gpf, glg, prefix, report, remote_name=target)
        report['n_maps'] = n_maps
    report['maps_ready'] = bool(list(work.glob(f'{prefix}.*.map')))
    report['map_files'] = sorted(p.name for p in work.glob(f'{prefix}.*.map'))
    (work / 'maps_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f"{target}: {n_tz} TZ atom(s), TZ->crystal sulfonamide N "
          f"{report.get('tz_to_crystal_sulfonamide_N_A')} A, npts {npts} "
          f"({size} A), maps_ready={report['maps_ready']} ({len(report['map_files'])} maps)")
    return report


def stage_dock(target, ligand, seeds, exhaustiveness, cpu, tag_suffix='', geometry='etkdg',
               xyz_override=None):
    from prep_ligands import input_mol_for_seed, write_pdbqt
    tdir = DOCK / target
    work = tdir / 'ad4zn'
    prefix = work / 'receptor_tz'
    if not list(work.glob('receptor_tz.*.map')):
        raise SystemExit(f'{target}: no maps in {work}; run --stage maps first')
    tag = f'{ligand}__ad4zn{tag_suffix}'
    rundir = tdir / 'runs' / tag
    rundir.mkdir(parents=True, exist_ok=True)
    box = json.loads((tdir / 'box.json').read_text())
    results = dict(target=target, ligand=ligand, scoring='ad4zn', box=box, cpu=cpu, tag=tag,
                   exhaustiveness=exhaustiveness, seeds=[], geometry=geometry,
                   maps=str(prefix.relative_to(DOCK)),
                   receptor=str((work / 'receptor_tz.pdbqt').relative_to(DOCK)))
    for seed in range(1, seeds + 1):
        sdir = rundir / f'seed{seed:02d}'
        sdir.mkdir(exist_ok=True)
        lig = sdir / 'ligand.pdbqt'
        # production 2026-09-28: --geometry dft docks the DFT structure (see prep_ligands.input_mol_for_seed)
        mol3d, prov = input_mol_for_seed(ligand, seed, geometry, xyz_override)
        results.setdefault('input_geometry', prov)
        write_pdbqt(mol3d, lig, ligand)
        out = sdir / 'out.pdbqt'
        cmd = ['nice', '-n', '10', 'vina', '--ligand', str(lig), '--maps', str(prefix),
               '--scoring', 'ad4', '--exhaustiveness', str(exhaustiveness),
               '--num_modes', '20', '--energy_range', '3', '--seed', str(seed),
               '--cpu', str(cpu), '--out', str(out)]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.time() - t0
        (sdir / 'vina.log').write_text(r.stdout + '\n' + r.stderr)
        if r.returncode != 0 or not out.exists():
            results['seeds'].append(dict(seed=seed, ok=False, wall_s=round(wall, 1),
                                         stderr=r.stderr[-1500:]))
            print(f'  seed {seed}: FAILED {r.stderr.strip()[-200:]}')
            continue
        scores = [float(ln.split()[3]) for ln in out.read_text().splitlines()
                  if ln.startswith('REMARK VINA RESULT')]
        results['seeds'].append(dict(seed=seed, ok=True, wall_s=round(wall, 1),
                                     n_modes=len(scores), scores=scores,
                                     best=scores[0] if scores else None,
                                     out=str(out.relative_to(DOCK))))
        print(f'  seed {seed:2d}: best {scores[0]:7.3f} kcal/mol, {len(scores)} modes, {wall:.0f} s')
    best = [s['best'] for s in results['seeds'] if s.get('ok')]
    if best:
        results['summary'] = dict(n_ok=len(best), best_mean=round(float(np.mean(best)), 3),
                                  best_sd=round(float(np.std(best, ddof=1)), 3) if len(best) > 1
                                  else 0.0,
                                  wall_s_total=round(sum(s['wall_s'] for s in results['seeds']), 1))
        print(f"{target} {ligand} ad4zn: {results['summary']['best_mean']:.3f} "
              f"+- {results['summary']['best_sd']:.3f} kcal/mol")
    (rundir / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', required=True)
    ap.add_argument('--stage', default='dock', choices=['maps', 'dock', 'both'])
    ap.add_argument('--ligand', default=None, choices=sorted(SPECIES))
    ap.add_argument('--vm', action='store_true', help='make the maps on the IBM VM')
    ap.add_argument('--seeds', type=int, default=10)
    ap.add_argument('--exhaustiveness', type=int, default=32)
    ap.add_argument('--cpu', type=int, default=3)
    ap.add_argument('--tag-suffix', default='',
                    help='appended to the run directory name (e.g. "__repro", "__exh64") so that a '
                         'repeat or a variant never overwrites the run it is compared with')
    ap.add_argument('--geometry', default='etkdg', choices=['etkdg', 'dft'],
                    help='etkdg: fresh ETKDG conformer per seed (validation); dft: the DFT geometry '
                         'from dock/ligands/<species>.sdf for every seed (production)')
    ap.add_argument('--geometry-xyz', default=None,
                    help='dock this xyz instead of the species geometry in dock/ligands/ (diagnostic; '
                         'implies --geometry dft)')
    ap.add_argument('--extra-types', nargs='*', default=[],
                    help='extra AutoDock ligand types to build maps for')
    a = ap.parse_args()
    if a.stage in ('maps', 'both'):
        stage_maps(a.target, a.vm, a.extra_types)
    if a.stage in ('dock', 'both'):
        if not a.ligand:
            raise SystemExit('--ligand is required for the docking stage')
        stage_dock(a.target, a.ligand, a.seeds, a.exhaustiveness, a.cpu, a.tag_suffix,
                   'dft' if a.geometry_xyz else a.geometry, a.geometry_xyz)
    return 0


if __name__ == '__main__':
    sys.exit(main())
