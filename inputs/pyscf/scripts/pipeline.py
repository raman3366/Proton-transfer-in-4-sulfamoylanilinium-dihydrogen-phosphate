"""
Orchestrator: runs every stage for a case, records status.json + STAGES.md + per-stage logs, and handles the
gas-phase proton-transfer check with an automatic IEF-PCM (water) rerun.
Usage: python dft/pipeline.py <case> [--from STAGE]
Stages: optimize, proton_check, frequencies, fmo, mwfn, plots, render, report
"""
import sys, os, subprocess, json, time, datetime, pathlib, shutil
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *

PY = sys.executable
PVBATCH = os.environ.get('PVBATCH', str(ROOT / 'tools' / 'bin' / 'pvbatch'))
if not os.path.exists(PVBATCH):
    PVBATCH = '/Applications/ParaView-6.1.1.app/Contents/bin/pvbatch'   # macOS fallback
STAGES = ['optimize', 'proton_check', 'frequencies', 'fmo', 'mwfn', 'plots', 'render', 'report']
# Frequencies are the dominant cost (6N displaced gradients).  With SKIP_GAS_FREQ=1 they are computed only
# for the solvent (primary) model; the gas-phase model still gets the optimisation, the proton-transfer
# check and all five wavefunction analyses, so every deliverable is present for both.
SKIP_GAS_FREQ = os.environ.get('SKIP_GAS_FREQ') == '1'
# With ALWAYS_PCM=1 the IEF-PCM rerun is started whether or not the gas-phase proton moved: for an ionic
# salt the continuum model is the one that represents the crystal, so it is always wanted.
ALWAYS_PCM = os.environ.get('ALWAYS_PCM') == '1'
# NO_AUTO_PCM=1 suppresses the automatic solvent rerun entirely -- used when every case (gas and PCM) is
# defined and scheduled explicitly, so the pipeline must not invent a duplicate.
NO_AUTO_PCM = os.environ.get('NO_AUTO_PCM') == '1'
# Stages whose failure must halt the pipeline.  'mwfn' belongs here: when PZC's Multiwfn stage failed the
# pipeline carried on and produced a RESULTS.md and 14 figures from partial data, reporting 7/8 rather than
# failing.  A report assembled from a crashed analysis is worse than no report.
STOP_ON_FAIL = {'optimize', 'fmo', 'mwfn'}
now = lambda: datetime.datetime.now().isoformat(timespec='seconds')

def load_status(case):
    p = case['dir'] / 'status.json'
    return json.load(open(p)) if p.exists() else {'case': case['name'], 'stages': {}}

def save_status(case, st):
    save_json(case['dir'] / 'status.json', st)

def proton_check(case):
    """Did the watched proton move from its donor to its acceptor during the optimisation? Also flags any change
    in covalent connectivity between the start and the optimised structure."""
    s0, x0 = read_xyz(ROOT / case['xyz']); s1, x1 = read_xyz(case['dir'] / 'opt.xyz')
    b0 = {(i, j) for i, j, _ in bonds_from_geometry(s0, x0)}; b1 = {(i, j) for i, j, _ in bonds_from_geometry(s1, x1)}
    lab = case.get('cif_labels') or atom_labels(s1)
    res = {'bonds_broken': [f'{lab[i]}-{lab[j]}' for i, j in sorted(b0 - b1)], 'bonds_formed': [f'{lab[i]}-{lab[j]}' for i, j in sorted(b1 - b0)]}
    if case.get('proton'):
        H, N, O = case['proton']['H'], case['proton']['donor_N'], case['proton']['acceptor_O']
        dHN, dHO = float(np.linalg.norm(x1[H] - x1[N])), float(np.linalg.norm(x1[H] - x1[O]))
        res.update({'d_H_N_A': dHN, 'd_H_O_A': dHO, 'proton_transferred': dHO < dHN})
    res['connectivity_changed'] = bool(res['bonds_broken'] or res['bonds_formed'])
    save_json(case['dir'] / 'proton_check.json', res)
    msg = (f"proton check: H-N {res.get('d_H_N_A', float('nan')):.3f} A, H-O {res.get('d_H_O_A', float('nan')):.3f} A -> "
           f"{'PROTON TRANSFERRED to the acceptor' if res.get('proton_transferred') else 'proton retained on the donor'}; "
           f"bonds broken {res['bonds_broken'] or 'none'}, formed {res['bonds_formed'] or 'none'}")
    stage_note(case, msg); print('[pipeline]', msg, flush=True)
    return res

def make_pcm_case(case):
    """Create the IEF-PCM(water) twin of a gas-phase case, starting again from the crystal geometry.

    `pzc_gas` becomes `pzc_pcm`, not `pzc_gas_pcm`, and the title says PCM instead of gas phase -- these
    names end up in file names, report headings and the hand-over package.
    """
    src = find_case_file(case['name'])
    new = json.load(open(src))
    new['solvent'] = 'water'
    title = new.get('title', case['name'])
    new['title'] = (title.replace('gas phase', 'IEF-PCM (water)') if 'gas phase' in title
                    else title + ' -- IEF-PCM (water)')
    why = ('because the gas-phase optimisation transferred the proton'
           if (case['dir'] / 'proton_check.json').exists()
           and json.load(open(case['dir'] / 'proton_check.json')).get('proton_transferred')
           else 'because the continuum model is the primary result for an ionic salt')
    new['notes'] = (new.get('notes', '') + f' | Rerun with an IEF-PCM (water) continuum {why}.').strip()
    base = case['name'][:-4] if case['name'].endswith('_gas') else case['name']
    name = f'{base}_pcm'
    json.dump(new, open(src.parent / f'{name}.json', 'w'), indent=2)
    return name

def run_stage(case, stage):
    st = load_status(case); logs = case['dir'] / 'logs'; logs.mkdir(exist_ok=True)
    if stage == 'frequencies' and SKIP_GAS_FREQ and not case.get('solvent'):
        st['stages'][stage] = {'status': 'skipped', 'note': 'frequencies are computed on the solvent (primary) model only'}
        save_status(case, st)
        stage_note(case, 'frequencies skipped: computed on the IEF-PCM (primary) model instead')
        print('[pipeline] frequencies: skipped (primary model only)', flush=True)
        return True
    st['stages'][stage] = {'status': 'running', 'start': now()}; save_status(case, st)
    t0 = time.time(); ok = True; note = ''
    if stage == 'proton_check':
        try:
            res = proton_check(case); note = 'proton transferred' if res.get('proton_transferred') else 'ionic form retained'
        except Exception as e:
            ok = False; note = f'{type(e).__name__}: {e}'
    else:
        script = ROOT / 'dft' / f'{stage}.py'
        if not script.exists():
            st['stages'][stage] = {'status': 'skipped', 'note': 'script not written yet'}; save_status(case, st)
            stage_note(case, f'{stage}: skipped (dft/{stage}.py not present yet)'); return True
        cmd = [PVBATCH, '--force-offscreen-rendering', str(script), case['name']] if stage == 'render' else [PY, str(script), case['name']]
        with open(logs / f'{stage}.log', 'w') as f:
            rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT).returncode
        ok = rc == 0; note = '' if ok else f'exit code {rc}, see logs/{stage}.log'
    st = load_status(case)
    st['stages'][stage] = {'status': 'done' if ok else 'failed', 'start': st['stages'][stage]['start'], 'end': now(), 'wall_min': round((time.time() - t0) / 60, 2), 'note': note}
    save_status(case, st)
    if not ok: stage_note(case, f'{stage} FAILED: {note}')
    print(f"[pipeline] {stage}: {'done' if ok else 'FAILED'} ({(time.time()-t0)/60:.1f} min) {note}", flush=True)
    return ok

def main(name, start_from='optimize'):
    case = load_case(name); (case['dir'] / 'logs').mkdir(exist_ok=True)
    stage_note(case, f"=== pipeline started for {name} (from stage '{start_from}') ===")
    st = load_status(case); st['pipeline_start'] = st.get('pipeline_start') or now(); save_status(case, st)
    for stage in STAGES[STAGES.index(start_from):]:
        ok = run_stage(case, stage)
        if not ok and stage in STOP_ON_FAIL:
            stage_note(case, f'pipeline stopped at {stage}'); return False
        # gas-phase proton transfer -> run the IEF-PCM(water) case FIRST (it is the salt form the paper is about),
        # then come back and finish the gas-phase stages for comparison
        pc = case['dir'] / 'proton_check.json'
        if stage == 'proton_check' and not NO_AUTO_PCM and not case.get('solvent') and pc.exists() and (
                ALWAYS_PCM or json.load(open(pc)).get('proton_transferred')):
            pcm = make_pcm_case(case)
            st = load_status(case); st['pcm_rerun'] = pcm; save_status(case, st)
            why = ('proton transferred in gas phase' if json.load(open(pc)).get('proton_transferred')
                   else 'ALWAYS_PCM: the continuum model is the primary result for an ionic salt')
            stage_note(case, f'{why} -> running IEF-PCM(water) case {pcm} first (primary result); gas-phase analyses follow')
            main(pcm)
            stage_note(case, f'PCM case {pcm} finished; resuming gas-phase stages')
    st = load_status(case); st['pipeline_end'] = now(); save_status(case, st)
    stage_note(case, f'=== pipeline finished for {name} ===')
    return True

if __name__ == '__main__':
    name = sys.argv[1]
    if '--stage' in sys.argv:                      # run one stage only (records status/notes like the full pipeline)
        sys.exit(0 if run_stage(load_case(name), sys.argv[sys.argv.index('--stage') + 1]) else 1)
    start = sys.argv[sys.argv.index('--from') + 1] if '--from' in sys.argv else 'optimize'
    sys.exit(0 if main(name, start) else 1)
