"""d14: numbers added when the last open audit items were closed (29 Sep 2026, night); no new calculation.
Reads only; writes derived/d14_completion.json.

(1) Measured binding of the two hCA II controls (acetazolamide AZM, sulfanilamide SFM) from the Protein-Ligand Binding
    Database (PLBD; Linge et al. 2023), snapshot paper/audit/research/plbd_hca2_azm_sfm.json (rows used, retrieved
    29 Sep 2026 through the PLBD REST interface).  Observed dissociation constants by fluorescent thermal shift assay
    (table kd_obs_tsa_all_cas, column ca02 = CA II, 37 C) are re-derived here as the mean K_d = 1/K_b,obs of the 37 C,
    pH 7.0 experiments that PLBD attributes to their primary publication (Linkuviene et al. 2018); the intrinsic
    constants (binding-linked protonation of ligand and protein subtracted) come from table intrinsic_parameters
    (K_b,int and dG_int), re-derived as K_d,int = 1/K_b,int and dG = -RT ln K_b,int.
(2) One definition of "near-native" for the AutoDock4Zn orientation gaps (numbers audit round 2, m5): a mode is
    near-native when its symmetry-corrected heavy-atom RMSD to the crystal ligand is at most 1.0 A (the gate of the
    text).  For the 6RL9 validation run (ETKDG input, 10 seeds) and for the eight AutoDock4Zn runs of sulfanilamide
    (production and diagnostic, 6RL9 and 3HS4) every mode of every seed is read from out.pdbqt and its RMSD computed
    here (own graph-isomorphism RMSD, no superposition, receptor frame), then per seed
      gap = score(best-ranked near-native mode) - score(top mode)       (0 when the top mode is near-native), and
      sep = |score(best near-native mode) - score(best mode beyond 2.0 A)|   (the score separation of the two
            orientations, where a seed has both).
    Check: the RMSDs of the validation run against dock/analysis/summary.json (_rmsds, the modes it kept).
(3) The sulfonamide torsion tau_fold of every docked top pose, measured on the raw pose coordinates of out.pdbqt
    (numbers audit burst 3, round 2, m5): Table S31 printed the value measured on the DFT molecule rebuilt onto each
    pose (production_summary.json, one decimal), which differs from the raw pose by up to 0.1 deg.  The relaxed DFT
    strain at the raw torsions is evaluated with the spline of d13 (C3 branch, converged points, and the minimum).
"""
import itertools

import networkx as nx
from networkx.algorithms import isomorphism
from scipy.interpolate import CubicSpline

from common import *

R_KJ = 8.314462618e-3            # kJ/(mol K)
PLBD = REPO / 'paper/audit/research/plbd_hca2_azm_sfm.json'

# ------------------------------------------------------------------ (1) PLBD
P = json.loads(PLBD.read_text())
T = P['tables']
cid = {c['id']: c['name'] for c in T['compound']}                    # 1783 AZM, 2255 SFM
PUB = [p['id'] for p in T['publication'] if p['label'].strip() == 'Linkuviene2018b'][0]   # doi:10.1017/S0033583518000082
assert cid == {1783: 'AZM', 2255: 'SFM'}, cid
plbd = dict(retrieved=P['retrieved'], database_state=P['database_state'], snapshot=rel(PLBD))
for name in ('AZM', 'SFM'):
    row = {int(r['id'].split('-')[1]): r for r in T['kd_obs_tsa_all_cas'] if r['name'] == name}
    rint = {int(r['id'].split('-')[1]): r for r in T['kd_int_tsa_all_cas'] if r['name'] == name}
    code = [k for k, v in cid.items() if v == name][0]
    exps = [e for e in T['experiment_tsa'] if e['compound_name'] == name]
    e37 = [e for e in exps if e['temperature'] == 37 and e['ph'] == 7 and e.get('primary_publication_id') == PUB]
    e37_all = [e for e in exps if e['temperature'] == 37 and e['ph'] == 7]
    ip = {r['temperature']: r for r in T['intrinsic_parameters'] if r['compound_id'] == code and r['protein_id'] == 2}
    kb37 = ip[37]['kb_int']
    plbd[name] = dict(
        kd_obs_37_M=row[37]['ca02'], kd_obs_25_M=row[25]['ca02'],
        kd_obs_37_rederived_M=float(np.mean([1 / e['kb_obs'] for e in e37])), n_exp_37=len(e37),
        exp_37=[dict(id=e['id'], label=e['label'], kb_obs=e['kb_obs'], ph=e['ph'], date=e.get('date')) for e in e37],
        kd_obs_37_all_exp_M=float(np.mean([1 / e['kb_obs'] for e in e37_all])), n_exp_37_all=len(e37_all),
        kd_int_tsa_37_M=rint[37]['ca02'],
        kb_int_37=kb37, kd_int_37_M=1 / kb37, dG_int_37_kJ=ip[37]['dbg_int'],
        dG_int_37_rederived_kJ=-R_KJ * 310.15 * math.log(kb37), n_int_37=ip[37]['n'],
        int_37_tsa_ids=ip[37]['best_intrinsic_tsa_ids'], int_37_itc_ids=ip[37]['best_intrinsic_itc_ids'])
plbd['ratio_obs_37'] = plbd['SFM']['kd_obs_37_M'] / plbd['AZM']['kd_obs_37_M']
plbd['ratio_int_37'] = plbd['SFM']['kd_int_37_M'] / plbd['AZM']['kd_int_37_M']

# ------------------------------------------------------------------ (2) orientation gaps
AD2EL = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O', 'O': 'O', 'S': 'S', 'SA': 'S',
         'H': 'H', 'HD': 'H', 'HS': 'H', 'P': 'P'}
RCOV = {'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07, 'H': 0.31}


def pdbqt_models(path):
    """[(score, [(element, xyz), ...]), ...] in the order of the file (mode 1 first)."""
    models, cur, score = [], None, None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = []
        elif ln.startswith('REMARK VINA RESULT'):
            score = float(ln.split()[3])
        elif ln.startswith(('ATOM', 'HETATM')) and cur is not None:
            cur.append((AD2EL[ln[77:79].strip()], np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
        elif ln.startswith('ENDMDL'):
            models.append((score, cur))
            cur = None
    return models


def heavy_graph(atoms):
    G = nx.Graph()
    heavy = [(e, x) for e, x in atoms if e != 'H']
    for i, (e, x) in enumerate(heavy):
        G.add_node(i, el=e, x=x)
    for i, j in itertools.combinations(range(len(heavy)), 2):
        if dist(heavy[i][1], heavy[j][1]) < RCOV[heavy[i][0]] + RCOV[heavy[j][0]] + 0.45:
            G.add_edge(i, j)
    return G


def sdf_graph(path):
    lines = Path(path).read_text().splitlines()
    na, nb = int(lines[3][:3]), int(lines[3][3:6])
    atoms = [(l[31:34].strip(), np.array([float(l[0:10]), float(l[10:20]), float(l[20:30])])) for l in lines[4:4 + na]]
    keep = [i for i, (e, _) in enumerate(atoms) if e != 'H']
    G = nx.Graph()
    for k, i in enumerate(keep):
        G.add_node(k, el=atoms[i][0], x=atoms[i][1])
    re_ = {i: k for k, i in enumerate(keep)}
    for l in lines[4 + na:4 + na + nb]:
        i, j = int(l[:3]) - 1, int(l[3:6]) - 1
        if i in re_ and j in re_:
            G.add_edge(re_[i], re_[j])
    return G


_maps = {}


def rmsd_sym(Gref, Gp):
    k = (id(Gref), Gp.number_of_nodes())
    if k not in _maps:
        _maps[k] = list(isomorphism.GraphMatcher(Gref, Gp, node_match=lambda a, b: a['el'] == b['el']).isomorphisms_iter())
        assert _maps[k], 'no isomorphism between reference and pose'
    R = np.array([Gref.nodes[i]['x'] for i in Gref.nodes])
    return min(float(np.sqrt(((R - np.array([Gp.nodes[m[i]]['x'] for i in Gref.nodes])) ** 2).sum(1).mean()))
               for m in _maps[k])


def orientation(run, ref):
    Gr = sdf_graph(DOCK / ref)
    seeds = []
    for f in sorted((DOCK / run).glob('seed*/out.pdbqt')):
        ms = pdbqt_models(f)
        rr = [rmsd_sym(Gr, heavy_graph(at)) for _, at in ms]
        sc = [s for s, _ in ms]
        near = [i for i, r in enumerate(rr) if r <= 1.0]
        far = [i for i, r in enumerate(rr) if r > 2.0]
        seeds.append(dict(seed=f.parent.name, n_modes=len(ms), rmsd=rr, score=sc,
                          rank_near=(near[0] + 1) if near else None,
                          gap=(sc[near[0]] - sc[0]) if near else None,
                          sep=abs(min(sc[i] for i in near) - min(sc[i] for i in far)) if near and far else None,
                          rmsd_top=rr[0]))
    gaps = [s['gap'] for s in seeds if s['gap'] is not None]
    seps = [s['sep'] for s in seeds if s['sep'] is not None]
    return dict(run=run, ref=ref, n_seeds=len(seeds), n_with_near=len(gaps), n_with_both=len(seps),
                gap_min=min(gaps), gap_max=max(gaps), sep_min=min(seps) if seps else None,
                sep_max=max(seps) if seps else None, worst_rank=max(s['rank_near'] for s in seeds if s['rank_near']),
                seeds=seeds)


val = orientation('T3_6RL9/runs/sfa_neutral__ad4zn', 'T3_6RL9/ref_SAN.sdf')
# check against the validation analysis (dock/analysis/summary.json keeps the RMSDs of the modes it read)
VS = json.loads((DOCK / 'analysis' / 'summary.json').read_text())['T3_6RL9/sfa_neutral__ad4zn']['per_seed']
val['check_rmsd_vs_summary'] = max(abs(a - b) for s, v in zip(val['seeds'], VS) for a, b in zip(s['rmsd'], v['_rmsds']))
val['check_top_score_vs_summary'] = max(abs(s['score'][0] - v['score_top']) for s, v in zip(val['seeds'], VS))
eight = {}
for tg, ref in (('T3_6RL9', 'ref_SAN.sdf'), ('T4_3HS4', 'ref_cross_T3_6RL9_SAN.sdf')):
    for suf in ('prod', 'csn102', 'csn109', 'etkdg'):
        o = orientation(f'{tg}/runs/sfa_neutral__ad4zn__{suf}', f'{tg}/{ref}')
        o.pop('seeds')
        eight[f'{tg}/{suf}'] = o
eight_all = dict(gap_min=min(o['gap_min'] for o in eight.values()), gap_max=max(o['gap_max'] for o in eight.values()),
                 n_runs=len(eight), n_seeds=sum(o['n_seeds'] for o in eight.values()),
                 all_seeds_near=all(o['n_with_near'] == o['n_seeds'] for o in eight.values()),
                 worst_rank=max(o['worst_rank'] for o in eight.values()))

# ------------------------------------------------------------------ (3) raw-pose torsions of the docked top poses
CUT = {('C', 'C'): 1.65, ('C', 'S'): 1.9, ('N', 'S'): 1.8}


def tau_fold(atoms):
    heavy = [(e, x) for e, x in atoms if e != 'H']
    s = [a for a in heavy if a[0] == 'S']
    assert len(s) == 1
    s = s[0]
    bond = lambda a, b: dist(a[1], b[1]) < CUT.get(tuple(sorted((a[0], b[0]))), 0)
    c4 = [a for a in heavy if a[0] == 'C' and bond(a, s)]
    n = [a for a in heavy if a[0] == 'N' and bond(a, s)]
    assert len(c4) == 1 and len(n) == 1
    ring = [a for a in heavy if a[0] == 'C' and a is not c4[0] and bond(a, c4[0])]
    assert len(ring) == 2
    return max(abs(dihedral(r[1], c4[0][1], s[1], n[0][1])) for r in ring)


D13 = json.loads((OUT / 'd13_referee_calcs.json').read_text())['m9']
splines = {}
for sp in ('sfa_neutral', 'sfa_anion'):
    m = D13[sp]
    c3 = sorted([p for p in m['points'] if p['branch'] == 'C3' and p['converged']], key=lambda p: p['tau_fold'])
    xs = np.array([m['tau_ref']] + [p['tau_fold'] for p in c3])
    ys = np.array([0.0] + [p['strain'] for p in c3])
    o = np.argsort(xs)
    splines[sp] = CubicSpline(xs[o], ys[o])
PS = json.loads((DOCK / 'results' / 'production_summary.json').read_text())
raw = {}
for key, jb in PS['jobs'].items():
    if not jb['ligand'].startswith('sfa_') or jb.get('diagnostic') or not key.endswith('__prod') or not jb.get('torsion'):
        continue
    tg, job = key.split('/')
    taus = [tau_fold(pdbqt_models(f)[0][1]) for f in sorted((DOCK / tg / 'runs' / job).glob('seed*/out.pdbqt'))]
    rebuilt = [q['tau_fold'] for q in jb['torsion']['per_seed']]
    r_ = dict(n=len(taus), tau_min=min(taus), tau_max=max(taus), tau_mean=float(np.mean(taus)),
              max_diff_vs_rebuilt=max(abs(a - b) for a, b in zip(sorted(taus), sorted(rebuilt))))
    if jb['ligand'] in splines:
        st = [float(splines[jb['ligand']](t)) for t in taus]
        r_.update(strain_min=min(st), strain_max=max(st))
    raw[key] = r_
tau_all = [v for r_ in raw.values() for v in (r_['tau_min'], r_['tau_max'])]
raw_all = dict(n_poses=sum(r_['n'] for r_ in raw.values()), tau_min=min(tau_all), tau_max=max(tau_all),
               max_diff_vs_rebuilt=max(r_['max_diff_vs_rebuilt'] for r_ in raw.values()))

dump('d14_completion.json', dict(plbd=plbd, orientation_validation=val, orientation_eight=eight,
                                 orientation_eight_all=eight_all, tau_raw=raw, tau_raw_all=raw_all))
print('PLBD 37 C: AZM Kd,obs %.1f nM (re-derived %.1f), SFM %.2f uM (re-derived %.2f); Kd,int %.2f / %.2f nM' % (
    plbd['AZM']['kd_obs_37_M'] * 1e9, plbd['AZM']['kd_obs_37_rederived_M'] * 1e9, plbd['SFM']['kd_obs_37_M'] * 1e6,
    plbd['SFM']['kd_obs_37_rederived_M'] * 1e6, plbd['AZM']['kd_int_37_M'] * 1e9, plbd['SFM']['kd_int_37_M'] * 1e9))
print('validation 6RL9 AD4Zn: gap %.3f-%.3f, separation %.3f-%.3f, RMSD check %.4f, score check %.4f' % (
    val['gap_min'], val['gap_max'], val['sep_min'], val['sep_max'], val['check_rmsd_vs_summary'], val['check_top_score_vs_summary']))
print('eight AD4Zn runs: gap %.3f-%.3f, worst rank %d, near-native in every seed %s' % (
    eight_all['gap_min'], eight_all['gap_max'], eight_all['worst_rank'], eight_all['all_seeds_near']))
print('raw-pose tau: %d poses, %.2f-%.2f, max diff vs rebuilt %.3f' % (raw_all['n_poses'], raw_all['tau_min'], raw_all['tau_max'],
                                                                      raw_all['max_diff_vs_rebuilt']))
