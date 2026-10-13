"""
Status table for the six requested analyses: progress %, elapsed, ETA. Usage: python dft/status.py <case>
Writes STATUS.md (repo root) and results/<case>/STATUS.md; prints the table.
"""
import sys, json, re, datetime, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *

EXPECTED_MIN = {'frequencies': 60, 'frequencies_fd': 230, 'fmo': 6, 'mwfn': 25, 'plots': 2, 'render': 6, 'report': 1}
EXPECTED_OPT_STEPS = 60

def fmt(m):
    if m is None: return '-'
    return f'{int(m)//60} h {int(m)%60:02d} min' if m >= 60 else f'{m:.0f} min'

def opt_progress(case):
    log = case['dir'] / 'logs' / 'optimize.log'
    if not log.exists(): return 0, None, None, ''
    txt = log.read_text(); steps = len(re.findall(r'^\[opt\] step', txt, re.M))
    m = re.findall(r'^\[opt\] step\s+(\d+)\s+E = (\S+) Ha\s+max\|g\| = (\S+)', txt, re.M)
    start = datetime.datetime.fromisoformat(json.load(open(case['dir'] / 'status.json'))['stages']['optimize']['start'])
    elapsed = (datetime.datetime.now() - start).total_seconds() / 60
    if steps == 0: return 0, elapsed, None, 'SCF of first step running'
    per = elapsed / steps; exp = max(EXPECTED_OPT_STEPS, steps + 8)
    eta = (exp - steps) * per; gmax = float(m[-1][2]) if m else float('nan')
    return steps / exp, elapsed, eta, f'step {steps}, max|g| {gmax:.1e} (target 4.5e-4), {per:.1f} min/step'

def table(case, offset_min=0.0, order_note=''):
    st = json.load(open(case['dir'] / 'status.json')) if (case['dir'] / 'status.json').exists() else {'stages': {}}
    S = st['stages']; status = lambda k: S.get(k, {}).get('status', 'pending')
    def elapsed_of(k):
        s = S.get(k, {})
        if 'wall_min' in s: return s['wall_min']
        if 'start' in s: return (datetime.datetime.now() - datetime.datetime.fromisoformat(s['start'])).total_seconds() / 60
        return None
    rows = []; remaining = 0.0
    # 1. geometry optimisation + frequencies
    p_opt, el_opt, eta_opt, note = opt_progress(case)
    if status('optimize') == 'done': p_opt, eta_opt = 1.0, 0.0
    fk = 'frequencies'; fd = (case['dir'] / 'logs' / 'frequencies.log').exists() and 'finite differences' in (case['dir'] / 'logs' / 'frequencies.log').read_text()
    exp_f = EXPECTED_MIN['frequencies_fd' if fd else 'frequencies']
    if status(fk) == 'done':
        p_f, eta_f = 1.0, 0.0; fj = json.load(open(case['dir'] / 'freq.json')); note = f"{fj['n_imaginary']} imaginary frequencies ({'true minimum' if fj['n_imaginary']==0 else 'NOT a minimum'}), lowest {fj['frequencies_cm-1'][0]:.0f} cm-1"
    elif status(fk) == 'running':
        p_f = min(0.95, elapsed_of(fk) / exp_f); eta_f = max(1, exp_f - elapsed_of(fk)); note = f"Hessian running ({'finite differences' if fd else 'analytic'})"
        cols = re.findall(r'FD column (\d+)/(\d+) done \(([\d.]+) min\)', (case['dir'] / 'logs' / 'frequencies.log').read_text()) if fd else []
        if cols:
            done, tot, tmin = int(cols[-1][0]), int(cols[-1][1]), float(cols[-1][2]); per = tmin / done
            p_f = done / tot; eta_f = (tot - done) * per; note = f'finite differences: {done}/{tot} displacement columns, {per:.1f} min each'
    else: p_f, eta_f = 0.0, (exp_f if status(fk) != 'failed' else 0)
    if status('optimize') == 'running': eta1 = (eta_opt or 0) + exp_f
    elif status(fk) == 'running': eta1 = eta_f
    elif status(fk) in ('done', 'failed'): eta1 = 0
    else: eta1 = (eta_opt or 0) + exp_f
    prog1 = 0.7 * p_opt + 0.3 * p_f
    stage1 = {'pending': 'waiting', 'running': 'optimising', 'done': 'optimised', 'failed': 'FAILED'}[status('optimize')]
    if status('optimize') == 'done': stage1 = {'pending': 'optimised, frequencies queued', 'running': 'frequencies running', 'done': 'optimised + frequencies done', 'failed': 'frequencies FAILED'}[status(fk)]
    rows.append(('1', 'Geometry optimisation + frequency check', stage1, prog1, (el_opt or 0) + (elapsed_of(fk) or 0), eta1, note)); remaining += eta1
    # 2. FMO
    st2 = status('fmo'); el2 = elapsed_of('fmo') or 0
    p2 = {'pending': 0, 'running': 0.5, 'done': 1, 'failed': 0, 'skipped': 0}[st2]; eta2 = {'pending': EXPECTED_MIN['fmo'], 'running': max(1, EXPECTED_MIN['fmo'] - el2)}.get(st2, 0)
    note2 = ''
    if st2 == 'done':
        f = json.load(open(case['dir'] / 'fmo.json')); note2 = f"HOMO {f['HOMO_eV']:.2f} eV, LUMO {f['LUMO_eV']:.2f} eV, gap {f['gap_eV']:.2f} eV"
    rows.append(('2', 'Frontier orbitals (HOMO/LUMO) + reactivity descriptors', {'pending': 'waiting for stage 1', 'running': 'single point running', 'done': 'done', 'failed': 'FAILED', 'skipped': 'skipped'}[st2], p2, el2, remaining + eta2, note2)); remaining += eta2
    # 3-6 depend on mwfn (Multiwfn) then plots/render
    stm, stp, str_ = status('mwfn'), status('plots'), status('render')
    elm = elapsed_of('mwfn') or 0; etam = {'pending': EXPECTED_MIN['mwfn'], 'running': max(1, EXPECTED_MIN['mwfn'] - elm)}.get(stm, 0)
    etapr = sum(EXPECTED_MIN[k] for k in ('plots', 'render', 'report') if status(k) in ('pending', 'running'))
    base = {'pending': 0, 'running': 0.35, 'done': 0.7, 'failed': 0, 'skipped': 0}[stm]
    fig_done = 0.15 * (stp == 'done') + 0.15 * (str_ == 'done')
    def stage36(figs_needed):
        if stm == 'pending': return 'waiting for stage 2'
        if stm == 'running': return 'Multiwfn running'
        if stm in ('failed', 'skipped'): return f'Multiwfn {stm}'
        done = [f for f in figs_needed if (case['dir'] / 'fig' / f).exists()]
        return 'done' if len(done) == len(figs_needed) else f'figures {len(done)}/{len(figs_needed)}'
    for num, name, figs, extra in [('3', 'Molecular electrostatic potential (MEP) surface', ['mep.png'], 'mep_extrema.csv'),
                                   ('4', 'RDG / non-covalent interaction (NCI) analysis', ['rdg_scatter.png', 'nci.png'], 'nci_interactions.md'),
                                   ('5', 'Electron localisation function (ELF) maps',
                                    [f'elf_{pl}.png' for pl in (case.get('planes') or {'plane': 0})], ''),
                                   ('6', 'Localised orbital locator (LOL) maps',
                                    [f'lol_{pl}.png' for pl in (case.get('planes') or {'plane': 0})], '')]:
        nf = sum((case['dir'] / 'fig' / f).exists() for f in figs)
        p = base + 0.3 * nf / len(figs) if stm == 'done' else base
        rows.append((num, name, stage36(figs), p, elm, remaining + etam + etapr, ''))
    hdr = f"### Status of {case['name']} — {datetime.datetime.now():%Y-%m-%d %H:%M}\n\n" + (order_note + '\n\n' if order_note else '')
    rows = [(num, name, stage, p, el, 0 if p >= 1 else ((eta + offset_min) if eta else eta), note) for num, name, stage, p, el, eta, note in rows]
    out = hdr + '| # | Analysis | Stage | Progress | Elapsed | ETA (from now) | Notes |\n|---|---|---|---|---|---|---|\n'
    for num, name, stage, p, el, eta, note in rows:
        out += f'| {num} | {name} | {stage} | {100*p:.0f}% | {fmt(el) if el else "-"} | {fmt(eta) if eta else ("done" if p >= 1 else "-")} | {note} |\n'
    pcs = case['dir'] / 'proton_check.json'
    if pcs.exists():
        pc = json.load(open(pcs)); L = case.get('cif_labels') or []
        pr = case.get('proton') or {}
        dn = L[pr['donor_N']] if pr and L else 'donor'
        ac = L[pr['acceptor_O']] if pr and L else 'acceptor'
        out += (f"\nProton check: H···{dn} {pc.get('d_H_N_A', 0):.3f} Å, H···{ac} {pc.get('d_H_O_A', 0):.3f} Å → "
                f"**{f'proton transferred to {ac}' if pc.get('proton_transferred') else f'proton retained on {dn}'}**\n")
    return out, max([r[5] or 0 for r in rows] + [0])

if __name__ == '__main__':
    case = load_case(sys.argv[1]); st0 = json.load(open(case['dir'] / 'status.json')) if (case['dir'] / 'status.json').exists() else {}
    pcm = st0.get('pcm_rerun'); out = ''
    if pcm:
        try:
            pcase = load_case(pcm); pst = json.load(open(pcase['dir'] / 'status.json')) if (pcase['dir'] / 'status.json').exists() else {}
            ptab, prem = table(pcase, 0.0, 'Primary result (salt form, IEF-PCM water). Runs first because the gas-phase optimisation transferred the proton.')
            prem = 0 if pst.get('pipeline_end') else prem
            gtab, _ = table(case, prem, f'Secondary result (gas phase, neutral acid–base complex after proton transfer). Its remaining stages start after `{pcm}` finishes' + (f' (about {fmt(prem)} from now).' if prem else '.'))
            out = ptab + '\n' + gtab
        except Exception as e:
            out, _ = table(case); out += f'\n(PCM case status not available yet: {e})\n'
    else:
        out, _ = table(case)
    open(ROOT / 'STATUS.md', 'w').write(out); open(case['dir'] / 'STATUS.md', 'w').write(out); print(out)
