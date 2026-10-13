"""Gate G-R2, step 0: what the three rehearsal runs visited at the sulfonamide nitrogen (md/N2.md, section 3.1).

    /root/cqc-sfa-tools/md/bin/python md/qm/n2_ranges.py      (MDAnalysis, NumPy; one thread; about a minute)

Reads the solute frames (every 10 ps) of equilibration stages 1 to 6 and of the 5 ns of production of the runs
TX0, TR0 and TN0 of md/test/REHEARSAL_2e.md (/root/cqc-sfa-md/runs_test/bonded_2e/).  Every frame is put into
the frame of the small model by a least-squares fit of 19 atoms of the zinc site (zinc; CB, CG, ND1, CE1, NE2,
CD2 of His94, His96, His119) on their places in md/qm/mcpb/CA2_small.pdb; then the six coordinates of
md/qm/n2_common.py are read from the ligand (and Zn-N2, the three angles at N2, the signed volume).
A frame is "X-like" when the signed volume at N2 is above +0.3, "R-like" below -0.3, "between" otherwise
(the limit 0.3 is the one of md/qm/gr_mm.py and gr_queue.py for a planar nitrogen).

Writes md/qm/n2/ranges.json (ranges per class and run, the design boxes of section 3.2 of md/N2.md) and
md/qm/n2/frames.csv (one line per frame).  The rehearsal is a test of the pipeline with terms that failed gate
G-R: the frames say where the force field as built goes, nothing about the molecule.
"""
import json
import pathlib

import numpy as np

import n2_common as nc

RUNS = pathlib.Path('/root/cqc-sfa-md/runs_test/bonded_2e')
STAGES = ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06', 'seg0001']
KEYS = ('s', 't', 'th1', 'chi', 'tau', 'th2', 'eta', 'th3')
PERIODIC = dict(chi=360.0, tau=360.0, eta=360.0)
PLANAR = 0.3
WIDEN = 1.5                                                                # "widened by half" (plan, point 3)


def circ_range(v, period):
    """Smallest arc that holds all values: (centre, half width)."""
    v = np.sort(np.mod(np.asarray(v, float), period))
    if len(v) == 1:
        return float(v[0]), 0.0
    gaps = np.diff(np.concatenate([v, v[:1] + period]))
    k = int(np.argmax(gaps))
    lo, hi = v[(k + 1) % len(v)], v[k]
    if hi < lo:
        hi += period
    return float((lo + hi) / 2.0), float((hi - lo) / 2.0)


def lin_range(v):
    return float((min(v) + max(v)) / 2.0), float((max(v) - min(v)) / 2.0)


def main():
    import MDAnalysis as mda
    g = nc.Geometry()
    m = g.m
    model = {(a['resi'], a['name']): a['xyz'] for a in m.atoms}
    fit_names = [(r, n) for r in (93, 95, 118) for n in ('CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2')] + [(260, 'ZN')]
    ref = np.array([model[k] for k in fit_names])
    rows = []
    for run in ('TX0', 'TR0', 'TN0'):
        d = RUNS / run
        u = mda.Universe(str(d / 'solute.pdb'), [str(d / f'solute_{s}.dcd') for s in STAGES])
        nfr = [mda.Universe(str(d / 'solute.pdb'), str(d / f'solute_{s}.dcd')).trajectory.n_frames for s in STAGES]
        sel = []
        for r, n in fit_names:
            a = u.select_atoms(f'resid {r} and name {n}')
            assert len(a) == 1, (r, n, len(a))
            sel.append(a[0].index)
        lig = []
        for n in m.names:
            a = u.select_atoms(f'resname SN1 and name {n}')
            assert len(a) == 1, n
            lig.append(a[0].index)
        stage_of = [s for s, k in zip(STAGES, nfr) for _ in range(k)]
        idx_in = [j for k in nfr for j in range(k)]
        for ts in u.trajectory:
            x = ts.positions.astype(float)
            r, pc, qc = nc.gp.kabsch(x[sel], ref)
            y = (r @ (x[sel] - pc).T).T + qc
            c = (r @ (x[lig] - pc).T).T + qc
            cd = g.coords(c)
            st = stage_of[ts.frame]
            # time after the release of the ligand (start of stage 4), ps; frames every 10 ps, written at the end of each interval
            t_rel = None
            if st in ('eq04', 'eq05', 'eq06', 'seg0001'):
                before = sum(k for s, k in zip(STAGES, nfr) if STAGES.index(s) < 3)
                t_rel = (ts.frame - before + 1) * 10.0
            sv = cd['signed_volume']
            rows.append(dict(run=run, stage=st, frame_in_stage=idx_in[ts.frame], frame=ts.frame, ps_after_release=t_rel,
                             cls='X' if sv > PLANAR else 'R' if sv < -PLANAR else 'between',
                             fit_rmsd=float(np.sqrt(((y - ref) ** 2).sum(1).mean())),
                             **{k: cd[k] for k in KEYS + ('signed_volume', 'Zn_N2', 'N2_S', 'N2_H', 'angle_sum',
                                                           'H_N2_S_C4', 'omega_improper', 'N94_Zn_N2', 'N96_Zn_N2',
                                                           'N119_Zn_N2')},
                             w=[float(v) for v in nc.unit(c[m.i['N2']] - m.zn)]))
    cols = [k for k in rows[0] if k != 'w'] + ['wx', 'wy', 'wz']
    out = [','.join(cols)]
    for r in rows:
        v = dict(r, wx=r['w'][0], wy=r['w'][1], wz=r['w'][2])
        out.append(','.join('' if v[k] is None else (f'{v[k]:.4f}' if isinstance(v[k], float) else str(v[k])) for k in cols))
    (nc.N2).mkdir(exist_ok=True)
    (nc.N2 / 'frames.csv').write_text('\n'.join(out) + '\n')

    def ranges(sel_rows):
        o = {}
        for k in KEYS + ('signed_volume', 'Zn_N2', 'angle_sum', 'N94_Zn_N2', 'N96_Zn_N2', 'N119_Zn_N2'):
            v = [r[k] for r in sel_rows]
            if not v:
                continue
            c, h = circ_range(v, PERIODIC[k]) if k in PERIODIC else lin_range(v)
            o[k] = dict(centre=round(c, 3), half_width=round(h, 3), low=round(c - h, 3), high=round(c + h, 3),
                        mean=round(float(np.mean(v)), 3) if k not in PERIODIC else None)
        return o
    doc = dict(what='what the rehearsal runs TX0, TR0, TN0 (terms at N2 as built, gate G-R failed) visited at N2; '
                    'written by md/qm/n2_ranges.py', frames_per_run={}, planar_limit=PLANAR, widen=WIDEN,
               fit_atoms=[f'{r}:{n}' for r, n in fit_names], classes={}, per_run={}, boxes={})
    for run in ('TX0', 'TR0', 'TN0'):
        rr = [r for r in rows if r['run'] == run]
        doc['frames_per_run'][run] = dict(n=len(rr), by_stage={s: sum(r['stage'] == s for r in rr) for s in STAGES},
                                          by_class={c: sum(r['cls'] == c for r in rr) for c in ('X', 'R', 'between')},
                                          largest_fit_rmsd=round(max(r['fit_rmsd'] for r in rr), 3))
        doc['per_run'][run] = {c: ranges([r for r in rr if r['cls'] == c]) for c in ('X', 'R', 'between')}
    for c in ('X', 'R', 'between'):
        doc['classes'][c] = dict(n=sum(r['cls'] == c for r in rows), ranges=ranges([r for r in rows if r['cls'] == c]))
    # the frames of the conversion of TR0: from the last frame held (end of stage 3) until 100 ps after the first
    # X-like frame that is followed by X-like frames only
    tr = [r for r in rows if r['run'] == 'TR0']
    last_not_x = max(k for k, r in enumerate(tr) if r['cls'] != 'X')
    first_free = min(k for k, r in enumerate(tr) if r['stage'] == 'eq04')
    conv = [r for k, r in enumerate(tr) if first_free - 3 <= k <= last_not_x + 10]
    doc['TR0_conversion'] = dict(
        rule='frames of TR0 from three frames before the release of the ligand (start of stage 4) to ten frames '
             'after the last frame that is not X-like',
        last_frame_not_X_like=dict(stage=tr[last_not_x]['stage'], frame_in_stage=tr[last_not_x]['frame_in_stage'],
                                   ps_after_release=tr[last_not_x]['ps_after_release']),
        frames=[dict(stage=r['stage'], frame_in_stage=r['frame_in_stage'], ps_after_release=r['ps_after_release'],
                     cls=r['cls'], signed_volume=round(r['signed_volume'], 4), w=r['w'],
                     **{k: round(r[k], 4) for k in KEYS}) for r in conv])
    tn = [r for r in rows if r['run'] == 'TN0' and r['cls'] != 'X']
    doc['TN0_frames_not_X_like'] = [dict(stage=r['stage'], frame_in_stage=r['frame_in_stage'], cls=r['cls'],
                                         signed_volume=round(r['signed_volume'], 4), w=r['w'],
                                         **{k: round(r[k], 4) for k in KEYS}) for r in tn]
    # design boxes: centre and half width of the range each class visited, and the same widened by half
    # (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 3); 'all': the range over every frame of the three runs
    for c, sel in (('X', [r for r in rows if r['cls'] == 'X']), ('R', [r for r in rows if r['cls'] == 'R']),
                   ('all', rows)):
        rg = ranges(sel)
        doc['boxes'][c] = {k: dict(centre=rg[k]['centre'], visited_half_width=rg[k]['half_width'],
                                   widened_half_width=round(WIDEN * rg[k]['half_width'], 3))
                           for k in ('s', 't', 'th1', 'chi', 'tau', 'th2', 'eta', 'th3')}
    nc.write_json(nc.N2 / 'ranges.json', doc)
    print(json.dumps(doc['frames_per_run'], indent=1))
    for c in ('X', 'R', 'between'):
        print(c, doc['classes'][c]['n'], {k: (v['low'], v['high']) for k, v in doc['classes'][c]['ranges'].items()})
    print('conversion frames:', len(conv), 'TN0 not X-like:', len(tn))
    for r in doc['TR0_conversion']['frames']:
        print(r['stage'], r['frame_in_stage'], r['ps_after_release'], r['cls'], r['signed_volume'],
              {k: round(r[k], 1) for k in KEYS})
    print('boxes', json.dumps(doc['boxes']))


if __name__ == '__main__':
    main()
