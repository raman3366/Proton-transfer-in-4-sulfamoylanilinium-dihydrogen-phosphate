"""
Dielectric ladder for the M1 ion pair: optimise at eps = 4, 8 and 20 (IEF-PCM) between the two environments
burst 1a already covers (eps 1 = gas, sfa_pair_gas_d3; eps 78.36 = water, sfa_pair_pcm_d3), and report where
the N-protonated (salt) minimum disappears as the dielectric screening is turned down.

Why: the referee critique separates two things that hold the proton on N1 -- bulk dielectric screening of the
ion pair and specific hydrogen bonding.  A continuum with a variable eps isolates the first.  The crystal's
effective dielectric constant is not known, so the ladder asks at which screening the salt form first survives
rather than assuming one value.

Usage
  python dft/eps_ladder.py <case> --eps 4 [--eps 8 ...] [--forms N,O]    # optimise (one process per eps is best)
  python dft/eps_ladder.py <case> --report [--gas-case sfa_pair_gas_d3]   # table / plot / crossover only
  python dft/eps_ladder.py <case> --refine 2                              # + 2 bisection points in the bracket
Options: --start XYZ (default results/<case>/opt.xyz, the burst-1a PCM-water optimum), --maxsteps N,
  --out DIR (default results/<case>/eps_ladder), and for code tests only --basis B, --conv loose.
Threads: DFT_THREADS (dft/proton_tools pins BLAS/OpenMP to it).

Method
  * Same level as the case (B3LYP-D3(BJ)/6-311++G(d,p), DF, grid 4, IEF-PCM Lebedev 29); only the dielectric
    constant changes (proton_tools.make_scanner sets with_solvent.eps on the case's PCM object, so cavity and
    discretisation are identical to the water run; common.make_mf is unchanged).  eps <= 1 means gas phase.
  * Form N (default) starts from the PCM-water optimum, i.e. the salt.  If the N-protonated minimum still
    exists at that eps the optimisation stays there; if it has disappeared the proton slides to O5.  The
    crossover eps* is bracketed by the largest eps whose optimum is O-protonated and the smallest eps above
    it that stays N-protonated; --refine bisects that bracket geometrically (sqrt(lo*hi)).
  * Form O (optional, --forms N,O) starts from the O-protonated guess (H1A 1.00 A from O5 on the O5->N1 line)
    and gives, where both forms are minima, dE(O-N) vs eps.  Its sign change (linear interpolation in the
    IEF-PCM factor (eps-1)/(eps+1)) is the thermodynamic crossover, as opposed to the point where the
    N-protonated minimum ceases to exist.
  * Convergence as in dft/scan.py: geomeTRIC defaults plus proton_tools.optimise's two acceptance rules
    (negligible forces; forces + energy); --maxsteps 150.
  * Resumable: result.json per eps and form; the latest geometry of a running optimisation is kept in
    partial.xyz and a restart continues from it.
  * Caveat printed in the report: the eps = 1 row comes from burst 1a's gas optimisation, which started from
    the X-ray geometry, not from the PCM-water optimum.  Run --eps 1 here if a same-start gas row is wanted.
Outputs: <out>/eps_<eps>_<form>/{opt.xyz, opt_traj.xyz, result.json}; <out>/eps_ladder.{json,csv,png}.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from proton_tools import *     # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
import argparse

ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
ap.add_argument('case')
ap.add_argument('--eps', action='append', default=[], help='dielectric constant(s); repeat or comma-separate')
ap.add_argument('--forms', default='N', help='N, O or N,O')
ap.add_argument('--report', action='store_true')
ap.add_argument('--refine', type=int, default=0, help='extra bisection points (form N) inside the bracket')
ap.add_argument('--gas-case', default='sfa_pair_gas_d3')
ap.add_argument('--start'); ap.add_argument('--maxsteps', type=int, default=150)
ap.add_argument('--out'); ap.add_argument('--basis')
ap.add_argument('--conv', choices=['default', 'loose'], default='default')
args = ap.parse_args()
EPS = [float(v) for e in args.eps for v in e.split(',') if v]
FORMS = [f.strip().upper() for f in args.forms.split(',') if f.strip()]
if not (EPS or args.report or args.refine):
    ap.error('give --eps, --report or --refine')

case = load_case(args.case)
if not case.get('solvent'):
    raise SystemExit('[eps] use the PCM case (e.g. sfa_pair_pcm_d3): its PCM settings are the template')
atoms = proton_atoms(case)
SYMS = atoms['syms']
OUT = pathlib.Path(args.out) if args.out else case['dir'] / 'eps_ladder'
OUT.mkdir(parents=True, exist_ok=True)
EPS_WATER = SOLVENT_EPS[case['solvent'].lower()]


def note(text):
    """Stage log: the case's STAGES.md for production runs; a test run (--out) keeps its notes in its own dir."""
    if args.out:
        with open(OUT / 'STAGES.md', 'a') as f:
            f.write(f"- {datetime.datetime.now():%Y-%m-%d %H:%M}  {text}\n")
    else:
        stage_note(case, text)


def epsdir(eps, form):
    return OUT / f'eps_{eps:g}_{form}'


def start_geometry():
    p = pathlib.Path(args.start) if args.start else case['dir'] / 'opt.xyz'
    if not p.exists():
        raise SystemExit(f'[eps] no start geometry: {p} is missing (burst-1a PCM-water optimum) -- give --start')
    s, x = read_xyz(p)
    assert s == SYMS, 'start geometry atom order differs from the case'
    if proton_metrics(atoms, x)['proton_on'] != 'N':
        # the water optimum is not a salt: form N must still start N-protonated to test for an N minimum
        print(f'[eps] {p} is O-protonated; form N starts from it with H1A put back 1.03 A from N1', flush=True)
        return f'{p} with H1A put back 1.03 A from N1', n_protonated_guess(atoms, x)
    return p, x


def run_point(eps, form):
    d = epsdir(eps, form); d.mkdir(exist_ok=True)
    if (d / 'result.json').exists():
        r = json.load(open(d / 'result.json'))
        print(f"[eps {eps:g} {form}] done earlier: delta {r['metrics']['delta_A']:+.4f} A, skipping", flush=True)
        return r
    sp, xs = start_geometry()
    x0, src = (xs, str(sp)) if form == 'N' else (o_protonated_guess(atoms, xs), f'{sp} with H1A moved 1.00 A from O5')
    if not (d / 'start.xyz').exists():               # the record of where this optimisation began
        write_xyz(d / 'start.xyz', SYMS, x0, f'start for eps={eps:g} form {form}: {src}')
    if (d / 'partial.xyz').exists():
        x0, src = read_xyz(d / 'partial.xyz')[1], 'resumed from partial.xyz'
    lvl = level_description(case, basis=args.basis, eps=eps)
    print(f'[eps {eps:g} {form}] {lvl}; start: {src}', flush=True)
    note(f'eps ladder: eps={eps:g} form {form} started ({lvl["environment"]})')
    gs, mol = make_scanner(case, d / 'start.xyz', basis=args.basis, eps=eps, chkfile=d / 'opt.chk')
    lvl['n_basis'] = int(mol.nao)
    res = optimise(gs, SYMS, x0, partial=d / 'partial.xyz', traj=d / 'opt_traj.xyz', maxsteps=args.maxsteps,
                   loose=(args.conv == 'loose'), tag=f'eps {eps:g} {form}')
    m = proton_metrics(atoms, res['xyz'])
    write_xyz(d / 'opt.xyz', SYMS, res['xyz'], f"{case['name']} eps={eps:g} form {form} E = {res['E_Ha']:.10f} Ha")
    r = {'case': case['name'], 'eps': eps, 'form': form, 'source': 'eps_ladder.py', 'started_from': src,
         'converged': res['converged'], 'accepted_on': res['accepted_on'], 'steps': res['steps'],
         'wall_min': round(res['wall_s'] / 60, 2),
         'E_Ha': res['E_Ha'], 'metrics': m, 'level': lvl, 'history': res['history'], **stamp()}
    atomic_json(d / 'result.json', r)
    (d / 'partial.xyz').unlink(missing_ok=True)
    note(f"eps ladder: eps={eps:g} form {form}: proton on {m['proton_on']}, delta {m['delta_A']:+.3f} A, "
                     f"{res['steps']} steps")
    print(f"[eps {eps:g} {form}] {'converged' if res['converged'] else 'NOT CONVERGED'} in {res['steps']} steps "
          f"({res['wall_s'] / 60:.1f} min): proton on {m['proton_on']}, delta {m['delta_A']:+.4f} A, N...O "
          f"{m['d_N1_O5_A']:.4f} A, P-O5 {m['d_P_O5_A']:.4f} A", flush=True)
    return r


def burst1a_row(case_name, eps):
    """eps = 1 (gas) and eps = 78.36 (water) come from the burst-1a optimisations of the full cases."""
    # load_case() exports the loaded case's DFT_DISP / DFT_XC to the environment, which make_mf() reads: without
    # restoring them, the --refine optimisations that follow a report would silently run at the gas case's level
    # (e.g. with D3 for the plain-B3LYP case sfa_pair_pcm_b3lyp).  Fixed in review, 28 Sep.
    keep = {k: os.environ.get(k) for k in ('DFT_DISP', 'DFT_XC')}
    try:
        c = load_case(case_name)
    except FileNotFoundError:
        return None
    finally:
        for k, v in keep.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    # a row at another level of theory (dispersion, functional or basis) does not belong in this ladder
    lvl = lambda cc: (str(cc.get('disp') or 'none').lower(), str(cc.get('xc') or XC).lower(),
                      str(cc.get('basis', BASIS)).lower())
    if lvl(c) != lvl(case) or (args.basis and args.basis.lower() != lvl(case)[2]):
        print(f'[eps] {case_name} is at another level of theory {lvl(c)} than this ladder '
              f'({args.basis or lvl(case)}): its eps = {eps:g} row is left out (run --eps {eps:g} here instead)',
              flush=True)
        return None
    p = c['dir'] / 'opt.xyz'
    if not p.exists():
        return None
    s = json.load(open(c['dir'] / 'opt_summary.json')) if (c['dir'] / 'opt_summary.json').exists() else {}
    x = read_xyz(p)[1]
    return {'case': case_name, 'eps': eps, 'form': 'N', 'source': f'burst 1a ({case_name}/opt.xyz)',
            'started_from': 'X-ray geometry (burst 1a pipeline)', 'converged': True, 'steps': s.get('steps'),
            'accepted_on': s.get('converged_on'),        # e.g. 1a's negligible-forces acceptance (finalize_opt.py)
            'E_Ha': s.get('E_final_Ha'), 'metrics': proton_metrics(atoms, x)}


def collect():
    rows = [json.load(open(p)) for p in sorted(OUT.glob('eps_*_*/result.json'))]
    have = {(round(r['eps'], 6), r['form']) for r in rows}
    for cname, e in ((args.gas_case, 1.0), (case['name'], EPS_WATER)):
        if (round(e, 6), 'N') not in have:
            r = burst1a_row(cname, e)
            if r:
                rows.append(r)
    return sorted(rows, key=lambda r: (r['eps'], r['form']))


def bracket(rows):
    """(eps_lo, eps_hi): largest eps whose form-N optimum is O-protonated, smallest eps above it that is
    N-protonated.  None when the proton never moves (or always moves) in the rows available."""
    n = sorted([r for r in rows if r['form'] == 'N'], key=lambda r: r['eps'])
    lo = [r['eps'] for r in n if r['metrics']['proton_on'] == 'O']
    if not lo:
        return None
    hi = [r['eps'] for r in n if r['eps'] > max(lo) and r['metrics']['proton_on'] == 'N']
    return (max(lo), min(hi)) if hi else None


def report():
    rows = collect()
    if not rows:
        print('[eps] nothing to report yet'); return
    br = bracket(rows)
    fe = lambda e: (e - 1.0) / (e + 1.0)                   # IEF-PCM's own dielectric factor
    # thermodynamic crossover: dE(O - N) at eps values where both forms converged to their own side
    dE = []
    for e in sorted({r['eps'] for r in rows}):
        rn = next((r for r in rows if r['eps'] == e and r['form'] == 'N' and r['metrics']['proton_on'] == 'N'), None)
        ro = next((r for r in rows if r['eps'] == e and r['form'] == 'O' and r['metrics']['proton_on'] == 'O'), None)
        if rn and ro and rn.get('E_Ha') is not None and ro.get('E_Ha') is not None:
            dE.append((e, (ro['E_Ha'] - rn['E_Ha']) * KCAL_PER_HA))
    thermo = None
    for (e1, v1), (e2, v2) in zip(dE, dE[1:]):
        if v1 * v2 < 0:
            f = fe(e1) + (fe(e2) - fe(e1)) * v1 / (v1 - v2)
            thermo = {'eps_star': (1 + f) / (1 - f), 'between': [e1, e2], 'interpolated_in': '(eps-1)/(eps+1)'}
    if br:
        cross = (f'N-protonated minimum disappears between eps = {br[0]:g} (proton on O5) and eps = {br[1]:g} '
                 f'(proton kept on N1)')
    elif all(r['metrics']['proton_on'] == 'N' for r in rows if r['form'] == 'N'):
        cross = (f"N-protonated minimum survives at every eps computed (down to eps = "
                 f"{min(r['eps'] for r in rows if r['form'] == 'N'):g})")
    elif all(r['metrics']['proton_on'] == 'O' for r in rows if r['form'] == 'N'):
        cross = (f"no N-protonated minimum at any eps computed (up to eps = "
                 f"{max(r['eps'] for r in rows if r['form'] == 'N'):g}): the proton goes to O5 everywhere")
    else:
        cross = 'no clean bracket (proton position not monotonic in eps, or missing rows) -- inspect the table'
    out = {'case': case['name'], 'crossover': cross, 'bracket_eps': br, 'thermodynamic_crossover': thermo,
           'dE_O_minus_N_kcal': [{'eps': e, 'dE_kcal': v} for e, v in dE],
           'rows': [{k: r.get(k) for k in ('eps', 'form', 'source', 'started_from', 'converged', 'accepted_on',
                                          'steps', 'E_Ha', 'metrics')} for r in rows],
           'caveat': 'The eps = 1 row (if from burst 1a) started from the X-ray geometry, the others from the '
                     'PCM-water optimum.'}
    atomic_json(OUT / 'eps_ladder.json', out)
    po_keys = sorted({k for r in rows for k in r['metrics']['P_O_A']})
    with open(OUT / 'eps_ladder.csv', 'w') as f:
        f.write('eps,form,source,converged,steps,E_Ha,proton_on,delta_A,d_N1_H1A_A,d_H1A_O5_A,d_N1_O5_A,d_P_O5_A,'
                + ','.join(po_keys) + '\n')
        for r in rows:
            m = r['metrics']
            f.write(f"{r['eps']:g},{r['form']},{r['source']},{r['converged']},{r.get('steps')},"
                    f"{'' if r.get('E_Ha') is None else format(r['E_Ha'], '.10f')},{m['proton_on']},"
                    f"{m['delta_A']:+.4f},{m['d_N1_H1A_A']:.4f},{m['d_H1A_O5_A']:.4f},{m['d_N1_O5_A']:.4f},"
                    f"{m['d_P_O5_A']:.4f}," + ','.join(f"{m['P_O_A'].get(k, float('nan')):.4f}" for k in po_keys) + '\n')
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for form, mk in (('N', 'o'), ('O', 's')):
        rr = [r for r in rows if r['form'] == form]
        if not rr:
            continue
        e = [r['eps'] for r in rr]
        ax.plot(e, [r['metrics']['delta_A'] for r in rr], mk + '-', label=f'delta, start {form}')
        ax.plot(e, [r['metrics']['d_N1_O5_A'] for r in rr], mk + '--', label=f'N1...O5, start {form}')
        ax.plot(e, [r['metrics']['d_P_O5_A'] for r in rr], mk + ':', label=f'P-O5, start {form}')
    ax.axhline(0, color='k', lw=0.6); ax.set_xscale('log'); ax.set_xlabel('eps (IEF-PCM; 1 = gas)')
    ax.set_ylabel('A'); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2)
    import textwrap
    ax.set_title('\n'.join(textwrap.wrap(f"{case['name']}: dielectric ladder -- {cross}", 90)), fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / 'eps_ladder.png', dpi=150); plt.close(fig)
    print(f'[eps] {cross}')
    for r in rows:
        m = r['metrics']
        print(f"   eps {r['eps']:7.3f} {r['form']}  proton on {m['proton_on']}  delta {m['delta_A']:+.3f}  "
              f"N...O {m['d_N1_O5_A']:.3f}  P-O5 {m['d_P_O5_A']:.3f}  [{r['source']}]")
    if thermo:
        print(f"   thermodynamic crossover (dE(O-N) = 0) at eps* ~ {thermo['eps_star']:.2f}")
    return br


for form in FORMS:
    for e in EPS:
        run_point(e, form)
br = None
with FileLock(OUT / '.report.lock'):
    br = report()
for _ in range(args.refine):
    if not br:
        print('[eps] --refine: no bracket to refine'); break
    e_new = float(f'{np.sqrt(br[0] * br[1]):.3g}')
    if e_new in (br[0], br[1]):
        break
    run_point(e_new, 'N')
    with FileLock(OUT / '.report.lock'):
        br = report()
