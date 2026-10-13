"""
Relaxed scan of the N1-H1A distance in the M1 ion pair: the proton H1A between its donor N1 (anilinium) and
the acceptor O5 of H2PO4-.  This is the central question of the paper -- at delta-pKa ~ 0, what holds the
proton on N1? -- asked as "what does it cost to hand the proton to O5, and is there a second minimum there?".

Usage
  python dft/scan.py <case> --chain A        # N side -> O side (outward)
  python dft/scan.py <case> --chain B        # O side -> N side (inward); independent of A, run in parallel
  python dft/scan.py <case> --report         # rebuild CSV / PNG / summary from the per-point files only
Options: --grid LO:HI:STEP (default 1.05:1.75:0.10), --points d1,d2,... (override the grid for this chain),
  --start XYZ (N-side start; default results/<case>/opt.xyz from burst 1a), --start-b XYZ (O-side start),
  --no-anchor, --b-min D, --maxsteps N (60), --anchor-maxsteps N (100), --retry-unconverged, --out DIR,
  and for code tests only: --basis B, --conv loose (plus GRID_LEVEL / PCM_LEBEDEV_ORDER in the environment).
Threads: DFT_THREADS (dft/proton_tools pins BLAS/OpenMP to it).

Method
  * Level of theory = the case's (B3LYP-D3(BJ)/6-311++G(d,p), DF, grid 4; IEF-PCM water, Lebedev 29, for
    sfa_pair_pcm_d3; gas for sfa_pair_gas_d3).  Only d(N1-H1A) is constrained, through geomeTRIC's
    "$set / distance 7 16 <value>" (1-based atom numbers; N1 = atom 7, H1A = atom 16 in the CIF order), with
    exact enforcement switched on near the target (proton_tools.ENFORCE); every other coordinate relaxes.
  * Grid: one absolute grid for both chains, anchored at 1.75 A and stepping down by 0.10 A (1.75, 1.65, ...,
    1.05), so the two chains evaluate the SAME d values and hysteresis is read off directly.
  * Chain A starts at the N-protonated minimum: an unconstrained re-optimisation of the burst-1a optimum (the
    "anchor": 1a accepted that geometry on forces after geomeTRIC crawled on a soft inter-ion mode, so expect a
    few gradients until the forces + energy rule below accepts it again) fixes d0, then every grid point above d0 in
    turn, each started from the previous converged geometry.  If the anchor loses the proton to O5 (expected
    in the gas phase) there is no N-side minimum; the anchor is stopped as for chain B and the chain covers
    the whole grid from its start geometry.
  * Chain B starts on the O side: H1A is put on the O5->N1 line 1.00 A from O5 and the structure optimised
    without constraints (anchor B).  If that stays O-protonated it is the O-side minimum and the chain walks
    down from it; if the proton returns to N1 there is no O-side minimum (itself a result: the anchor is
    stopped three steps after the proton has crossed, proton_tools.proton_left, and kept out of the profile)
    and the chain walks down the whole grid from the O-protonated guess.  It stops above d0 of the start.
  * Why two chains: in a relaxed scan the other coordinates (N...O distance, phosphate and NH3+ rotations)
    follow the valley they came from, so a single chain can overshoot a crossing and give a profile that
    depends on the direction.  The lower envelope of both chains is the minimum-energy profile estimate; the
    A-B difference at shared points is the hysteresis measure.
  * Convergence: geomeTRIC's defaults (= Gaussian's) at every point, plus two acceptance rules for 3 gradients in
    a row (proton_tools.optimise): forces 10x inside the thresholds, or forces inside the thresholds with
    |dE| < 1e-6 Ha (only the displacement criterion unmet), so a point cannot spend dozens of steps on a free
    inter-ion rotation (burst 1a's PCM pair optimisations crawled like that to steps 64 and 68).  Step caps (--maxsteps 60 per point, --anchor-maxsteps 100) bound a chain's wall
    time; a point that hits its cap is kept, flagged unconverged in the CSV, and the chain continues from it.
  * Resumable: each point writes <chain>_<kk>_<tag>.json + .xyz when done, the latest geometry of a running
    point goes to <...>.partial.xyz after every gradient, and the chain's SCF checkpoint seeds the first SCF
    after a restart.  Rerunning the same command continues where it stopped.
Outputs (results/<case>/scan/): per-point JSON/XYZ, <chain>_traj.xyz, scan_profile.csv, scan_summary.json,
scan_profile.png (a working plot, not publication styling).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from proton_tools import *     # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
import argparse

ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
ap.add_argument('case')
ap.add_argument('--chain', choices=['A', 'B'])
ap.add_argument('--report', action='store_true')
ap.add_argument('--grid', default='1.05:1.75:0.10')
ap.add_argument('--points', help='comma-separated d(N1-H1A) values (A) for this chain, overriding the grid')
ap.add_argument('--start'); ap.add_argument('--start-b')
ap.add_argument('--no-anchor', action='store_true')
ap.add_argument('--b-min', type=float, help='chain B stops above this d (default: d0 of the start + 0.03 A)')
# Step caps bound the wall time of a chain (~3.5 min per PCM gradient on the loaded VM): a constrained point
# started from its converged neighbour needs 10-30 gradients; one that hits the cap is kept, flagged unconverged,
# and the chain moves on from its last geometry (--retry-unconverged revisits it).
ap.add_argument('--maxsteps', type=int, default=60, help='per constrained point (default 60)')
ap.add_argument('--anchor-maxsteps', type=int, default=100, help='per unconstrained anchor (default 100)')
ap.add_argument('--retry-unconverged', action='store_true')
ap.add_argument('--out'); ap.add_argument('--basis'); ap.add_argument('--conv', choices=['default', 'loose'],
                                                                      default='default')
args = ap.parse_args()
if not (args.chain or args.report):
    ap.error('give --chain A|B or --report')

case = load_case(args.case)
atoms = proton_atoms(case)
SYMS, N, H, O = atoms['syms'], atoms['N'], atoms['H'], atoms['O']
OUT = pathlib.Path(args.out) if args.out else case['dir'] / 'scan'
OUT.mkdir(parents=True, exist_ok=True)
GAP = 0.03          # a grid point closer than this to an anchor's d adds nothing: skip it


def note(text):
    """Stage log: the case's STAGES.md for production runs; a test run (--out) keeps its notes in its own dir."""
    if args.out:
        with open(OUT / 'STAGES.md', 'a') as f:
            f.write(f"- {datetime.datetime.now():%Y-%m-%d %H:%M}  {text}\n")
    else:
        stage_note(case, text)


def grid():
    lo, hi, step = (float(v) for v in args.grid.split(':'))
    g = [round(hi - k * step, 4) for k in range(int(round((hi - lo) / step)) + 1)]
    return sorted(v for v in g if v >= lo - 1e-9)


def n_side_start():
    p = pathlib.Path(args.start) if args.start else case['dir'] / 'opt.xyz'
    if not p.exists():
        raise SystemExit(f'[scan] no N-side start geometry: {p} is missing (burst-1a optimum) -- give --start')
    s, x = read_xyz(p)
    assert s == SYMS, 'start geometry atom order differs from the case'
    if proton_metrics(atoms, x)['proton_on'] != 'N':
        # e.g. an optimum in which the proton already moved to O5: chain A needs a salt geometry to start from
        print(f'[scan] {p} is O-protonated; the N side starts from it with H1A put back 1.03 A from N1', flush=True)
        return f'{p} with H1A put back 1.03 A from N1', n_protonated_guess(atoms, x)
    return p, x


# ----------------------------------------------------------------------------------------------- chain -----
def run_chain(chain):
    start_path, xA = n_side_start()
    d0_start = proton_metrics(atoms, xA)['d_N1_H1A_A']
    if chain == 'A':
        x_first, first_src = xA, str(start_path)
    else:
        if args.start_b:
            x_first, first_src = read_xyz(args.start_b)[1], args.start_b
        else:
            x_first, first_src = o_protonated_guess(atoms, xA), f'{start_path} with H1A moved 1.00 A from O5'
    lvl = level_description(case, basis=args.basis)
    print(f"[scan {chain}] case {case['name']}  {lvl}\n[scan {chain}] start: {first_src}", flush=True)
    note(f"scan chain {chain} started ({lvl['environment']}, out {OUT})")
    gs = None
    traj = OUT / f'{chain}_traj.xyz'

    def scanner():
        nonlocal gs
        if gs is None:
            tmpxyz = OUT / f'.{chain}_mol.xyz'
            write_xyz(tmpxyz, SYMS, x_first, 'scanner template')
            gs, mol = make_scanner(case, tmpxyz, basis=args.basis, chkfile=OUT / f'{chain}.chk')
            tmpxyz.unlink(missing_ok=True)
            lvl['n_basis'] = int(mol.nao)
        return gs

    def used_indices():
        """{index: tag} of every point of this chain on disk (finished or partial)."""
        out = {}
        for p in OUT.glob(f'{chain}_[0-9][0-9]_*'):
            if p.name.endswith(('.json', '.partial.xyz')):
                out[int(p.name[2:4])] = p.name[5:].replace('.partial.xyz', '').replace('.json', '')
        return out

    def index_for(k, tag):
        """(chain, index) must stay unique: fsp.py finds a point's geometry by <chain>_<index>_*.xyz and the report
        keys the lower envelope on it.  A point already on disk keeps its index; a new point takes k unless another
        point holds it (e.g. an extra --points run, such as zpe_barrier.py's refine suggestion, after the chain),
        then the next free index.  (Fixed in review, 28 Sep: extra points used to reuse index 1.)"""
        used = used_indices()
        mine = [i for i, t in used.items() if t == tag]
        if mine:
            return min(mine)
        return k if k not in used else max(used) + 1

    def do_point(k, kind, target, x_prev, src):
        tag = 'anchor' if kind == 'anchor' else f'd{target:.3f}'
        k = index_for(k, tag)
        stem = OUT / f'{chain}_{k:02d}_{tag}'
        # (not with_suffix: the stem itself contains the decimal point of the distance)
        js, xf, part, cf = (pathlib.Path(f'{stem}{ext}') for ext in ('.json', '.xyz', '.partial.xyz', '.constraints'))
        if js.exists():
            r = json.load(open(js))
            if r.get('status') == 'done' and (r['converged'] or r.get('stopped') or not args.retry_unconverged):
                print(f"[scan {chain}] {stem.name}: done earlier (E = {r['E_Ha']:.8f}), skipping", flush=True)
                return r, read_xyz(xf)[1]
        x0, how = x_prev, src
        if part.exists():
            x0, how = read_xyz(part)[1], f'resumed from {part.name}'
        elif js.exists() and xf.exists():
            x0, how = read_xyz(xf)[1], f'retry of unconverged {xf.name}'
        print(f"[scan {chain}] {stem.name}: {kind} {'' if target is None else f'd(N1-H1A) = {target:.3f} A '}"
              f"from {how}", flush=True)
        cons = None if target is None else (N, H, target)
        res = optimise(scanner(), SYMS, x0, constraint=cons, cfile=cf,
                       partial=part, traj=traj, loose=(args.conv == 'loose'),
                       maxsteps=(args.anchor_maxsteps if kind == 'anchor' else args.maxsteps),
                       tag=f'scan {chain} {tag}',
                       stop_if=(proton_left(atoms, 'N' if chain == 'A' else 'O') if kind == 'anchor' else None))
        m = proton_metrics(atoms, res['xyz'])
        dev = None if target is None else abs(m['d_N1_H1A_A'] - target)
        r = {'case': case['name'], 'chain': chain, 'index': k, 'kind': kind, 'target_d_A': target,
             'status': 'done', 'converged': res['converged'], 'steps': res['steps'],
             'wall_min': round(res['wall_s'] / 60, 2), 'E_Ha': res['E_Ha'], 'metrics': m,
             'constraint_dev_A': dev, 'started_from': how, 'stopped': res['stopped'],
             'accepted_on': res['accepted_on'], 'level': lvl,
             'history': res['history'], **stamp()}
        write_xyz(xf, SYMS, res['xyz'], f"{case['name']} scan {chain} {tag}  E = {res['E_Ha']:.10f} Ha  "
                                        f"d(N1-H1A) = {m['d_N1_H1A_A']:.4f} A")
        atomic_json(js, r)
        part.unlink(missing_ok=True)
        print(f"[scan {chain}] {stem.name}: {'converged' if res['converged'] else ('STOPPED' if res['stopped'] else 'NOT CONVERGED')} in "
              f"{res['steps']} steps, {res['wall_s'] / 60:.1f} min; E = {res['E_Ha']:.8f}; d(N-H) "
              f"{m['d_N1_H1A_A']:.4f}, d(H...O) {m['d_H1A_O5_A']:.4f}, delta {m['delta_A']:+.4f}, P-O5 "
              f"{m['d_P_O5_A']:.4f} A" + ('' if dev is None else f'; |d - target| = {dev:.1e} A'), flush=True)
        return r, res['xyz']

    k, x = 0, x_first
    expect = 'N' if chain == 'A' else 'O'
    if not args.no_anchor:
        r, xa = do_point(0, 'anchor', None, x_first, first_src)
        k = 1
        if r['metrics']['proton_on'] == expect:
            x = xa
            d_ref = r['metrics']['d_N1_H1A_A']
        else:
            print(f"[scan {chain}] the unconstrained anchor ended with the proton on "
                  f"{r['metrics']['proton_on']}: there is no {expect}-side minimum at this level; the chain "
                  f"continues from its start geometry", flush=True)
            note(f'scan chain {chain}: no {expect}-protonated minimum (anchor went to '
                             f"{r['metrics']['proton_on']})")
            d_ref = None
    else:
        d_ref = proton_metrics(atoms, x_first)['d_N1_H1A_A'] if chain == 'A' else None
    if args.points:
        pts = [float(v) for v in args.points.split(',')]
    elif chain == 'A':
        # no N-side minimum (anchor went to O5): the constrained chain has to cover the whole grid
        lo = d_ref + GAP if d_ref is not None else -1.0
        pts = [g for g in grid() if g > lo]
    else:
        top = d_ref - GAP if d_ref is not None else float('inf')
        bot = args.b_min if args.b_min is not None else d0_start + GAP
        pts = [g for g in grid()[::-1] if bot < g < top]
    print(f"[scan {chain}] constrained points: {', '.join(f'{p:.3f}' for p in pts)}", flush=True)
    src = first_src if (args.no_anchor or d_ref is None) else f'anchor {chain}_00'

    def nearest_done(target):
        """Geometry of this chain's finished point nearest in d(N1-H1A) (anchors stopped early excluded)."""
        best = None
        for p in OUT.glob(f'{chain}_[0-9][0-9]_*.json'):
            r = json.load(open(p))
            if r.get('status') != 'done' or r.get('stopped'):
                continue
            dist = abs(r['metrics']['d_N1_H1A_A'] - target)
            if best is None or dist < best[0]:
                best = (dist, p)
        if best is None:
            return None
        return read_xyz(best[1].with_suffix('.xyz'))[1], f'nearest finished point {best[1].stem}'

    for p in pts:
        x0, s0 = x, src
        if args.points:
            # extra points (e.g. a refinement at the barrier top) start from the chain's nearest finished point
            near = nearest_done(p)
            if near is not None:
                x0, s0 = near
        r, x = do_point(k, 'constrained', p, x0, s0)
        k, src = k + 1, 'previous point'
    note(f'scan chain {chain} finished: {k} points')


# ---------------------------------------------------------------------------------------------- report -----
def report():
    pts = load_points(OUT)                                   # stationary / constrained points: the profile
    allp = load_points(OUT, profile_only=False)              # + anchors stopped early, for the record
    if not pts:
        print('[scan] no finished points yet'); return
    env = lower_envelope(pts)
    ids_env = {(r['chain'], r['index']) for r in env}
    # reference: the N-side anchor of chain A if it is a genuine N-protonated minimum, else the lowest point
    ref = next((r for r in pts if r['chain'] == 'A' and r['kind'] == 'anchor' and r['metrics']['proton_on'] == 'N'),
               min(pts, key=lambda r: r['E_Ha']))
    E0 = ref['E_Ha']
    cols = ['chain', 'index', 'kind', 'target_d_A', 'd_N1_H1A_A', 'd_H1A_O5_A', 'delta_A', 'd_N1_O5_A',
            'angle_N1_H1A_O5_deg', 'd_P_O5_A', 'E_Ha', 'E_rel_kcal', 'converged', 'steps', 'wall_min',
            'constraint_dev_A', 'proton_on', 'on_lower_envelope', 'stopped_early']
    rows = []
    for r in sorted(allp, key=lambda r: (r['chain'], r['index'])):
        m = r['metrics']
        rows.append([r['chain'], r['index'], r['kind'], '' if r['target_d_A'] is None else f"{r['target_d_A']:.3f}",
                     f"{m['d_N1_H1A_A']:.4f}", f"{m['d_H1A_O5_A']:.4f}", f"{m['delta_A']:+.4f}",
                     f"{m['d_N1_O5_A']:.4f}", f"{m['angle_N1_H1A_O5_deg']:.1f}", f"{m['d_P_O5_A']:.4f}",
                     f"{r['E_Ha']:.10f}", f"{(r['E_Ha'] - E0) * KCAL_PER_HA:.3f}", r['converged'], r['steps'],
                     r['wall_min'], '' if r['constraint_dev_A'] is None else f"{r['constraint_dev_A']:.1e}",
                     m['proton_on'], (r['chain'], r['index']) in ids_env, bool(r.get('stopped'))])
    tmp = OUT / 'scan_profile.csv.tmp'
    with open(tmp, 'w') as f:
        f.write(','.join(cols) + '\n')
        for row in rows:
            f.write(','.join(str(v) for v in row) + '\n')
    os.replace(tmp, OUT / 'scan_profile.csv')

    # hysteresis at grid points evaluated by both chains
    byd = {}
    for r in pts:
        if r['target_d_A'] is not None:
            byd.setdefault(round(r['target_d_A'], 3), {})[r['chain']] = r['E_Ha']
    hyst = {f'{d:.3f}': (v['A'] - v['B']) * KCAL_PER_HA for d, v in sorted(byd.items()) if len(v) == 2}
    # stationary points on the lower envelope (discrete: a point lower / higher than both neighbours)
    E = np.array([r['E_Ha'] for r in env]); dd = np.array([r['metrics']['d_N1_H1A_A'] for r in env])
    mins = [i for i in range(len(E)) if (i == 0 or E[i] < E[i - 1]) and (i == len(E) - 1 or E[i] < E[i + 1])]
    maxs = [i for i in range(1, len(E) - 1) if E[i] > E[i - 1] and E[i] > E[i + 1]]
    summ = {'case': case['name'], 'scan_dir': str(OUT), 'n_points': len(pts),
            'n_unconverged': sum(not r['converged'] for r in pts),
            'reference': {'chain': ref['chain'], 'index': ref['index'], 'kind': ref['kind'], 'E_Ha': E0},
            'level': pts[0]['level'],
            'lower_envelope': [{'d_N1_H1A_A': float(dd[i]), 'E_rel_kcal': float((E[i] - E0) * KCAL_PER_HA),
                                'delta_A': env[i]['metrics']['delta_A'], 'd_P_O5_A': env[i]['metrics']['d_P_O5_A'],
                                'chain': env[i]['chain'], 'index': env[i]['index']} for i in range(len(env))],
            # a minimum at either end of the envelope that is not an unconstrained anchor is only the edge of the
            # scanned range (the energy is still falling there), not a stationary point
            'envelope_minima': [{'d_N1_H1A_A': float(dd[i]), 'E_rel_kcal': float((E[i] - E0) * KCAL_PER_HA),
                                 'proton_on': env[i]['metrics']['proton_on'],
                                 'at_scan_edge': bool(i in (0, len(E) - 1) and env[i]['kind'] != 'anchor')}
                                for i in mins],
            'envelope_maxima': [{'d_N1_H1A_A': float(dd[i]), 'E_rel_kcal': float((E[i] - E0) * KCAL_PER_HA)}
                                for i in maxs],
            'hysteresis_A_minus_B_kcal': hyst,
            'max_abs_hysteresis_kcal': (max(abs(v) for v in hyst.values()) if hyst else None),
            'anchors': {r['chain']: {'proton_on': r['metrics']['proton_on'], 'd_N1_H1A_A': r['metrics']['d_N1_H1A_A'],
                                     'E_rel_kcal': (r['E_Ha'] - E0) * KCAL_PER_HA, 'converged': r['converged'],
                                     'stopped_early': r.get('stopped'),
                                     'minimum_on_expected_side': bool(r['converged'] and r['metrics']['proton_on']
                                                                      == ('N' if r['chain'] == 'A' else 'O'))}
                        for r in allp if r['kind'] == 'anchor'},
            'note': ('Energies are electronic (no ZPE); see dft/zpe_barrier.py.  envelope_minima/maxima are the '
                     'discrete stationary points of the lower envelope at the 0.10 A grid resolution.')}
    atomic_json(OUT / 'scan_summary.json', summ)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 7.2), sharex=True)
    for ch, c, mk in (('A', 'C0', 'o'), ('B', 'C3', 's')):
        rr = sorted([r for r in pts if r['chain'] == ch], key=lambda r: r['metrics']['d_N1_H1A_A'])
        if not rr:
            continue
        x = [r['metrics']['d_N1_H1A_A'] for r in rr]
        a1.plot(x, [(r['E_Ha'] - E0) * KCAL_PER_HA for r in rr], mk + '-', color=c, mfc='none',
                label=f"chain {ch} ({'N->O' if ch == 'A' else 'O->N'})")
        a2.plot(x, [r['metrics']['delta_A'] for r in rr], mk + '-', color=c, mfc='none', label=f'delta, {ch}')
        a2.plot(x, [r['metrics']['d_P_O5_A'] for r in rr], mk + ':', color=c, label=f'P-O5, {ch}')
        for r in rr:
            if not r['converged']:
                a1.plot(r['metrics']['d_N1_H1A_A'], (r['E_Ha'] - E0) * KCAL_PER_HA, 'x', color='k', ms=10)
    a1.plot(dd, (E - E0) * KCAL_PER_HA, '-', color='0.6', lw=4, alpha=0.4, zorder=0, label='lower envelope')
    a1.set_ylabel('E - E(ref) (kcal/mol)'); a1.legend(fontsize=8); a1.grid(alpha=0.3)
    a1.set_title(f"{case['name']}: relaxed N1-H1A scan ({pts[0]['level']['basis']}, "
                 f"{pts[0]['level']['environment']})", fontsize=9)
    a2.axhline(0, color='k', lw=0.6); a2.set_xlabel('d(N1-H1A) (A)'); a2.set_ylabel('A'); a2.grid(alpha=0.3)
    a2.legend(fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(OUT / 'scan_profile.png', dpi=150); plt.close(fig)
    print(f"[scan] report: {len(pts)} points -> {OUT / 'scan_profile.csv'}, scan_summary.json, scan_profile.png")
    for r in env:
        m = r['metrics']
        print(f"   {r['chain']}{r['index']:02d} d(N-H) {m['d_N1_H1A_A']:.3f}  delta {m['delta_A']:+.3f}  "
              f"P-O5 {m['d_P_O5_A']:.3f}  E_rel {(r['E_Ha'] - E0) * KCAL_PER_HA:8.3f} kcal/mol"
              f"{'' if r['converged'] else '  (NOT converged)'}")
    if hyst:
        print('   hysteresis A-B (kcal/mol): ' + ', '.join(f'{d}: {v:+.3f}' for d, v in hyst.items()))


if args.chain:
    run_chain(args.chain)
with FileLock(OUT / '.report.lock'):
    report()
