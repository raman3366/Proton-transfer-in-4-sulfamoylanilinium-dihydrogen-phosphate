"""Gate G-R2, step 1: the new geometries and the split into fitting and test set (md/N2.md, section 3).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_points.py [--check]
    (NumPy, SciPy, OpenMM Reference platform; one thread; about two minutes on the Pi; no quantum chemistry)

Builds on md/qm/gr_points.py (imported, not edited) through md/qm/n2_common.py.  Needs md/qm/n2/ranges.json
(md/qm/n2_ranges.py) and md/qm/n2/ff/small_as_built.prmtop.

THE CONSTRUCTION is that of gate G-R: zinc and the three ethyl-imidazoles fixed, the ligand body B rigid with
the internal geometry of the crystal pose, Zn-N2 1.931 A, N2-H 1.021 A.  A geometry is given by the six
coordinates of n2_common.Geometry: the place of N2 (s on the great circle X -> Rx), th1 = Zn-N2-S, chi =
NE2(His94)-Zn-N2-S, tau = Zn-N2-S-C4, th2 = Zn-N2-H, eta = azimuth of H about N2->Zn counted from S.

THE DESIGN (fixed here, before any energy of a new point exists)
 (D) designed points: a scrambled Sobol sequence in 7 dimensions, seed 20261006, taken in its order:
       v      -> lam = clip((v - 0.35) / 0.35, 0, 1): 35 % of the points in the X-like box (lam = 0), 30 % in
                 the R-like box (lam = 1), 35 % between the two (0 < lam < 1), where the nitrogen passes
                 through the plane;
       layer  -> below 0.5: the half widths the rehearsal runs visited ("core"); else these widened by half
                 ("margin"), the limit of the plan;
       5 numbers -> th1, chi, tau, th2, eta = centre(lam) + scale * half_width(lam) * (2 u - 1), with centre
                 and half width interpolated linearly in lam between the X-like and the R-like class of
                 md/qm/n2/ranges.json (chi and tau the short way round, eta through 180 = planar);
       s = lam, t = 0: N2 stays on the line between its places in X and in Rx (its place on the sphere about
                 the zinc is set by zinc terms, which the plan does not allow to change).
 (F) frames of the rehearsal mapped onto the rigid model: the 19 frames of the conversion of TR0 and the 3
     frames of TN0 that are not X-like (ranges.json).  From the frame, after the fit of the zinc site on the
     model: the direction Zn -> N2 and th1, chi, tau, th2, eta.  From the model: Zn-N2, N2-H, the body B, the
     zinc and the imidazoles.
 EXCLUDED, each with its reason in md/qm/n2/exclusions.csv: a heavy atom of the ligand within 2.79 A of a heavy
 atom of an imidazole (the shortest such contact in the 2,100 frames of the three rehearsal runs; the design had
 2.2 A, changed after the design check and before any energy of a new point existed: closer contacts test the
 Lennard-Jones wall, which no term at N2 can repair); a force-field energy as built more than 20 kcal/mol above
 the lowest (lowest over the
 36 rigid points of gate G-R, the frames (F) and the first 1024 Sobol candidates that pass the contact rule).
 NUMBER AND ORDER: 144 points in the main list (the frames (F) that pass, one after every five designed
 points, and designed points in Sobol order) and 48 more designed points in a reserve list that follows.  The
 batch computes in this order; any beginning of the list is a balanced sample.
 THE SPLIT: the list is cut into blocks of four in its order; in every block one point is drawn for the test
 set by numpy.random.default_rng(20261006).  Written to md/qm/n2/split.json with its checksum in
 split.json.sha256.  A second run of this script must give the same files (--check compares and writes nothing).

Writes md/qm/n2/points.json, points.csv, exclusions.csv, split.json, split.json.sha256, geom/<id>.xyz.
"""
import argparse
import json

import numpy as np

import n2_common as nc

N_MAIN, N_RESERVE, N_CAND = 144, 48, 1024
CONTACT_HEAVY = 2.79      # was 2.2 in the design; raised 6 Oct 2026, 2:32 PM IST, before any energy of a new point existed (md/N2.md, section 5.1)
E_CAP = 20.0
CORE_SHARE, SHARE_X, SHARE_BETWEEN = 0.5, 0.35, 0.35
DESIGN_KEYS = ('th1', 'chi', 'tau', 'th2', 'eta')
PERIODIC = ('chi', 'tau')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help='compare with the files that exist; write nothing')
    args = ap.parse_args()
    from scipy.stats import qmc
    import scipy
    assert nc.gr_points_unchanged(), 'md/qm/gr_points.py is not the version this was written for'
    g = nc.Geometry()
    m = g.m
    ff = nc.FF()
    rng_doc = json.loads((nc.N2 / 'ranges.json').read_text())
    box = rng_doc['boxes']

    # ---- checks on the 36 rigid points of gate G-R
    old_ids = [p for p in g.old if p not in ('R', 'N')]
    mm_gate = json.loads((nc.GR / 'mm_2e_main.json').read_text())['points']
    built = nc.as_built_parameters(ff)
    worst_back, worst_e, const = 0.0, 0.0, []
    old = {}
    for pid in old_ids:
        c = g.old[pid]
        cd = g.coords(c)
        back = g.construct(cd['th1'], cd['chi'], cd['tau'], cd['th2'], cd['eta'], s=cd['s'], t=cd['t'])
        worst_back = max(worst_back, float(np.abs(back - c).max()))
        e = ff.energy(g.full(c))
        worst_e = max(worst_e, abs(e['total'] - mm_gate[pid]['energy']['total']))
        f = nc.features(cd)
        const.append(e['total'] - e['rest'] - sum(built[k] * f[k] for k in f))
        old[pid] = dict(coords=cd, e=e, f=f)
    assert worst_back < 2e-3, f'construct(coords()) does not give the old points back: {worst_back}'
    assert worst_e < 1e-4, f'force-field energies differ from md/qm/gr/mm_2e_main.json: {worst_e}'
    spread = float(max(const) - min(const))
    assert spread < 1e-4, f'rest + model terms does not reproduce the total: spread {spread}'

    # ---- candidates
    def record(c, kind, extra):
        cd = g.coords(c)
        con, flags = m.contacts(c)
        e = ff.energy(g.full(c))
        return dict(kind=kind, coords=cd, contacts=con, protein_contact='protein_contact' in flags,
                    heavy_heavy=con['model']['heavy_heavy']['d'], e=e, f=nc.features(cd), c=c, **extra)

    cx, cr = box['X'], box['R']

    def centre_hw(k, lam):
        a, b = cx[k]['centre'], cr[k]['centre']
        d = nc.wrap(b - a) if k in PERIODIC else b - a
        return a + lam * d, (1 - lam) * cx[k]['visited_half_width'] + lam * cr[k]['visited_half_width']

    sob = qmc.Sobol(d=7, scramble=True, seed=nc.SEED).random(N_CAND)
    cand = []
    for n, row in enumerate(sob):
        lam = float(np.clip((row[0] - SHARE_X) / SHARE_BETWEEN, 0.0, 1.0))
        scale = 1.0 if row[1] < CORE_SHARE else rng_doc['widen']
        q = {}
        for k, u in zip(DESIGN_KEYS, row[2:]):
            c0, hw = centre_hw(k, lam)
            q[k] = float(c0 + scale * hw * (2 * u - 1))
        c = g.construct(q['th1'], q['chi'], q['tau'], q['th2'], q['eta'], s=lam, t=0.0)
        cand.append(record(c, 'design', dict(sobol_index=n, lam=round(lam, 6), layer='core' if scale == 1.0 else 'margin',
                                              design={k: round(v, 4) for k, v in q.items()})))
    frames = []
    for src, rows in (('TR0', rng_doc['TR0_conversion']['frames']), ('TN0', rng_doc['TN0_frames_not_X_like'])):
        for r in rows:
            c = g.construct(r['th1'], r['chi'], r['tau'], r['th2'], r['eta'], w=r['w'])
            frames.append(record(c, 'frame', dict(run=src, stage=r['stage'], frame_in_stage=r['frame_in_stage'],
                                                   frame_class=r['cls'], frame_signed_volume=r['signed_volume'])))
    ok_contact = [p for p in cand + frames if p['heavy_heavy'] >= CONTACT_HEAVY]
    e_low = min([old[p]['e']['total'] for p in old] + [p['e']['total'] for p in ok_contact])
    low_at = ([p for p in old if old[p]['e']['total'] == e_low] or ['a new candidate'])[0]
    excl = []

    def passes(p, label):
        if p['heavy_heavy'] < CONTACT_HEAVY:
            excl.append((label, p['kind'], 'contact', p['heavy_heavy'], p['contacts']['model']['heavy_heavy']['pair'],
                         round(p['e']['total'] - e_low, 3)))
            return False
        if p['e']['total'] - e_low > E_CAP:
            excl.append((label, p['kind'], 'energy', p['heavy_heavy'], '', round(p['e']['total'] - e_low, 3)))
            return False
        return True
    fr_ok = [p for k, p in enumerate(frames)
             if passes(p, f"frame {p['run']} {p['stage']} {p['frame_in_stage']}")]
    de_ok, used = [], 0
    for p in cand:
        used += 1
        if passes(p, f"sobol {p['sobol_index']}"):
            de_ok.append(p)
        if len(de_ok) == N_MAIN + N_RESERVE - len(fr_ok):
            break
    assert len(de_ok) == N_MAIN + N_RESERVE - len(fr_ok), 'not enough candidates pass: raise N_CAND'
    order, fi, di = [], 0, 0
    while len(order) < N_MAIN:
        if fi < len(fr_ok) and len(order) % 6 == 5:
            order.append(fr_ok[fi])
            fi += 1
        else:
            order.append(de_ok[di])
            di += 1
    assert fi == len(fr_ok), 'not every frame found a place in the main list'
    order += de_ok[di:]
    assert len(order) == N_MAIN + N_RESERVE
    nd = nf = 0
    for p in order:
        if p['kind'] == 'design':
            nd += 1
            p['id'] = f'n{nd:03d}'
        else:
            nf += 1
            p['id'] = f'f{nf:02d}'
    rng = np.random.default_rng(nc.SEED)
    sets = {}
    for b in range(0, len(order), 4):
        t = int(rng.integers(0, 4))
        for j, p in enumerate(order[b:b + 4]):
            sets[p['id']] = 'test' if j == t else 'fit'

    # ---- files
    files = {}
    listing = []
    for n, p in enumerate(order, 1):
        pid = p['id']
        what = (f"designed, Sobol {p['sobol_index']}, lam {p['lam']}, {p['layer']}" if p['kind'] == 'design'
                else f"frame of {p['run']} {p['stage']} {p['frame_in_stage']} mapped on the rigid model")
        files[f'geom/{pid}.xyz'] = g.xyz_text(p['c'], f'G-R2 point {pid}: {what}; frame of 6RL9, Angstrom')
        cd = p['coords']
        sv = cd['signed_volume']
        listing.append(dict(
            order=n, id=pid, list='main' if n <= N_MAIN else 'reserve', kind=p['kind'], definition=what,
            **{k: p[k] for k in ('sobol_index', 'lam', 'layer', 'design', 'run', 'stage', 'frame_in_stage',
                                 'frame_class', 'frame_signed_volume') if k in p},
            region='X' if sv > 0.3 else 'R' if sv < -0.3 else 'between',
            coords={k: round(float(v), 5) for k, v in cd.items()},
            features={k: round(v, 8) for k, v in p['f'].items()},
            mm_as_built=dict(total=round(p['e']['total'], 6), rest=round(p['e']['rest'], 6),
                             above_lowest=round(p['e']['total'] - e_low, 4),
                             rel_X=round(p['e']['total'] - old['X']['e']['total'], 4)),
            contacts=p['contacts'], protein_contact=p['protein_contact'],
            xyz_sha256=nc.hashlib.sha256(files[f'geom/{pid}.xyz'].encode()).hexdigest(),
            ligand_coordinates={nm: [round(float(t), 8) for t in v] for nm, v in zip(m.names, p['c'])}))
    count = lambda key, val, lst: sum(1 for q in listing if q['list'] in lst and q[key] == val)   # noqa: E731
    summary = {}
    for lst in (('main',), ('main', 'reserve')):
        tag = '+'.join(lst)
        summary[tag] = dict(
            n=sum(q['list'] in lst for q in listing),
            fit=sum(sets[q['id']] == 'fit' for q in listing if q['list'] in lst),
            test=sum(sets[q['id']] == 'test' for q in listing if q['list'] in lst),
            by_region={r: count('region', r, lst) for r in ('X', 'between', 'R')},
            by_kind={r: count('kind', r, lst) for r in ('design', 'frame')},
            test_by_region={r: sum(1 for q in listing if q['list'] in lst and q['region'] == r and sets[q['id']] == 'test')
                            for r in ('X', 'between', 'R')},
            protein_contact=sum(q['protein_contact'] for q in listing if q['list'] in lst),
            mm_above_lowest=dict(max=max(q['mm_as_built']['above_lowest'] for q in listing if q['list'] in lst),
                                 median=float(np.median([q['mm_as_built']['above_lowest'] for q in listing if q['list'] in lst]))))
    doc = dict(
        what='new points of gate G-R2 (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 3); written by md/qm/n2_points.py',
        model=dict(natoms=64, charge=1, multiplicity=1, fixed_atoms='1-46', ligand_atoms='47-64', ligand_names=m.names),
        seed=nc.SEED, sobol=dict(d=7, scramble=True, candidates=N_CAND, used=used, scipy=scipy.__version__),
        shares=dict(X=SHARE_X, between=SHARE_BETWEEN, R=round(1 - SHARE_X - SHARE_BETWEEN, 2), core=CORE_SHARE),
        limits=dict(contact_heavy_heavy_A=CONTACT_HEAVY, mm_above_lowest_kcal=E_CAP),
        mm_lowest=dict(total=round(e_low, 6), at=low_at),
        boxes={k: box[k] for k in ('X', 'R', 'all')}, ranges_sha256=nc.sha(nc.N2 / 'ranges.json'),
        force_field=dict(prmtop='md/qm/n2/ff/small_as_built.prmtop', sha256=nc.sha(ff.prmtop),
                         terms_taken_out_for_rest=[[round(v, 5) if isinstance(v, float) else v for v in r] for r in ff.removed['rest']]),
        checks=dict(old_points_back_from_their_coordinates_A=worst_back,
                    mm_total_against_gate_file_kcal=worst_e, rest_plus_model_terms_spread_kcal=spread,
                    gr_points_sha256=nc.GR_POINTS_SHA),
        old_rigid_points={pid: dict(coords={k: round(float(v), 5) for k, v in old[pid]['coords'].items()},
                                    features={k: round(v, 8) for k, v in old[pid]['f'].items()},
                                    mm_as_built=dict(total=round(old[pid]['e']['total'], 6),
                                                     rest=round(old[pid]['e']['rest'], 6))) for pid in old_ids},
        n_main=N_MAIN, n_reserve=N_RESERVE, excluded=len(excl), frames_given=len(frames), frames_kept=len(fr_ok),
        summary=summary, order=[q['id'] for q in listing], points=listing,
        versions=dict(numpy=np.__version__, scipy=scipy.__version__, openmm=ff.mm.__version__))
    files['points.json'] = json.dumps(doc, indent=1) + '\n'
    split = dict(
        what='fitting and test set of gate G-R2, drawn before any energy of a new point was computed '
             '(md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 3); written by md/qm/n2_points.py',
        seed=nc.SEED, rule='the list in its order is cut into blocks of four; in each block the point with index '
                           'numpy.random.default_rng(seed).integers(0, 4), drawn block after block, is the test point',
        prefix_rule='the batch computes in this order and starts no new point after its deadline; the data of the '
                    'gate are the points finished then; nothing is left out afterwards',
        order=[q['id'] for q in listing], n_main=N_MAIN, set=sets,
        xyz_sha256={q['id']: q['xyz_sha256'] for q in listing},
        old_rigid_points=old_ids, old_points_not_used=['R', 'N'],
        counts=dict(fit=sum(v == 'fit' for v in sets.values()), test=sum(v == 'test' for v in sets.values())))
    files['split.json'] = json.dumps(split, indent=1) + '\n'
    files['split.json.sha256'] = nc.hashlib.sha256(files['split.json'].encode()).hexdigest() + '  split.json\n'
    cols = ['order', 'id', 'list', 'set', 'kind', 'region', 'lam', 'layer', 's', 't', 'th1', 'chi', 'tau', 'th2', 'eta',
            'th3', 'angle_sum', 'signed_volume', 'H_N2_S_C4', 'omega_improper', 'mm_above_lowest', 'mm_rel_X',
            'heavy_heavy', 'protein_contact']
    rows = [','.join(cols)]
    for q in listing:
        v = dict(q['coords'], **q, set=sets[q['id']], mm_above_lowest=q['mm_as_built']['above_lowest'],
                 mm_rel_X=q['mm_as_built']['rel_X'], heavy_heavy=q['contacts']['model']['heavy_heavy']['d'])
        rows.append(','.join(str(v.get(k, '')) for k in cols))
    files['points.csv'] = '\n'.join(rows) + '\n'
    files['exclusions.csv'] = '\n'.join(['candidate,kind,reason,closest_heavy_heavy_A,pair,mm_above_lowest_kcal'] +
                                        [','.join(str(t) for t in e) for e in excl]) + '\n'
    if args.check:
        same = {k: (nc.N2 / k).exists() and (nc.N2 / k).read_text() == v for k, v in files.items()}
        print('identical to the files on disk:', all(same.values()), [k for k, v in same.items() if not v][:5])
        raise SystemExit(0 if all(same.values()) else 1)
    if (nc.N2 / 'qm').exists() and any((nc.N2 / 'qm').glob('*.json')):
        same = (nc.N2 / 'split.json').exists() and (nc.N2 / 'split.json').read_text() == files['split.json']
        if not same:
            raise SystemExit('energies of new points exist and the split would change: nothing written')
    (nc.N2 / 'geom').mkdir(parents=True, exist_ok=True)
    for k, v in files.items():
        (nc.N2 / k).write_text(v)
    print(f"{len(listing)} points ({N_MAIN} main, {N_RESERVE} reserve); Sobol candidates used {used}; excluded {len(excl)} "
          f"(contact {sum(e[2] == 'contact' for e in excl)}, energy {sum(e[2] == 'energy' for e in excl)}); "
          f"frames kept {len(fr_ok)} of {len(frames)}")
    print('lowest force-field energy as built:', round(e_low, 4), 'at', low_at)
    print(json.dumps(summary, indent=1))
    print('checks:', doc['checks'])
    print('split.json sha256', files['split.json.sha256'].split()[0])


if __name__ == '__main__':
    main()
