"""Gate G-R2, a look that is NOT a result: what a candidate set of terms at N2 would change in the energies of the
frames of the three rehearsal runs (md/N2.md, section 6.7).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_frames.py --model MODEL.json [--out FILE]

For each of the 2,100 solute frames that md/qm/n2_ranges.py read (equilibration stages 1 to 6 and production of TX0,
TR0, TN0; the frame is put into the frame of the small model in the same way): the energy of the terms at N2 with
the parameters of MODEL.json minus their energy with the parameters as built, from the functions of
n2_common.features.  Only the terms that the plan allows to change enter, so this is exactly the energy by which
the candidate force field differs from the force field as built at that frame, apart from one constant (and
apart from the charges, if the second step is taken: then this number is the bonded part alone).  Reported: the
mean and spread of that difference for the frames with the configuration of X, for the R-like frames and for those
between, and the difference of the means (how much the candidate shifts the R-like frames against the X-like ones,
on the geometries the runs with the OLD terms visited).  The runs were made with the old terms; a run with the new
terms would visit other geometries.  No quantum-chemical energy is read.
"""
import argparse
import json
import pathlib

import numpy as np

import n2_common as nc

RUNS = pathlib.Path('/root/cqc-sfa-md/runs_test/bonded_2e')
STAGES = ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06', 'seg0001']
PLANAR = 0.3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    import MDAnalysis as mda
    model = json.loads(pathlib.Path(args.model).read_text())
    par = model['parameters']
    built = nc.as_built_parameters(nc.FF(nc.N2 / 'ff' / 'small_as_built.prmtop'))
    dp = {k: par.get(k, 0.0) - built[k] for k, _ in nc.FEATURES}
    g = nc.Geometry()
    m = g.m
    mdl = {(a['resi'], a['name']): a['xyz'] for a in m.atoms}
    fit_names = [(r, n) for r in (93, 95, 118) for n in ('CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2')] + [(260, 'ZN')]
    ref = np.array([mdl[k] for k in fit_names])
    rows = []
    for run in ('TX0', 'TR0', 'TN0'):
        d = RUNS / run
        u = mda.Universe(str(d / 'solute.pdb'), [str(d / f'solute_{s}.dcd') for s in STAGES])
        nfr = [mda.Universe(str(d / 'solute.pdb'), str(d / f'solute_{s}.dcd')).trajectory.n_frames for s in STAGES]
        stage_of = [s for s, k in zip(STAGES, nfr) for _ in range(k)]
        sel = [u.select_atoms(f'resid {r} and name {n}')[0].index for r, n in fit_names]
        lig = [u.select_atoms(f'resname SN1 and name {n}')[0].index for n in m.names]
        for ts in u.trajectory:
            x = ts.positions.astype(float)
            r, pc, qc = nc.gp.kabsch(x[sel], ref)
            c = (r @ (x[lig] - pc).T).T + qc
            cd = g.coords(c)
            f = nc.features(cd)
            sv = cd['signed_volume']
            rows.append(dict(run=run, stage=stage_of[ts.frame], frame=ts.frame, signed_volume=sv,
                             cls='X' if sv > PLANAR else 'R' if sv < -PLANAR else 'between',
                             de=float(sum(dp[k] * f[k] for k in dp)),
                             parts={grp: float(sum(dp[k] * f[k] for k, gg in nc.FEATURES if gg == grp)) for grp in
                                    ('angle', 'tor_H', 'tor_Zn', 'tor_ZnN_S', 'improper', 'tor_ZnN_H')}))
    de = np.array([r['de'] for r in rows])
    cls = np.array([r['cls'] for r in rows])
    base = float(de[cls == 'X'].mean())
    summ = {}
    for c in ('X', 'between', 'R'):
        v = de[cls == c] - base
        if len(v):
            summ[c] = dict(frames=int(len(v)), mean=round(float(v.mean()), 3), sd=round(float(v.std()), 3),
                           lowest=round(float(v.min()), 3), highest=round(float(v.max()), 3),
                           parts_mean={grp: round(float(np.mean([r['parts'][grp] for r in rows if r['cls'] == c])
                                                        - np.mean([r['parts'][grp] for r in rows if r['cls'] == 'X'])), 3)
                                       for grp in rows[0]['parts']})
    tr0 = [r for r in rows if r['run'] == 'TR0' and r['stage'] in ('eq03', 'eq04')]
    k0 = next(i for i, r in enumerate(tr0) if r['stage'] == 'eq04')
    conv = [dict(ps_after_release=(i - k0 + 1) * 10, cls=r['cls'], signed_volume=round(r['signed_volume'], 3),
                 change=round(r['de'] - base, 3)) for i, r in enumerate(tr0) if -3 <= i - k0 <= 12]
    out = dict(what='energy by which a candidate set of terms at N2 differs from the terms as built, at the frames of the three '
                    'rehearsal runs (made with the terms as built); zero = mean over the frames with the configuration of X; '
                    'a look, not a result; written by md/qm/n2_frames.py', written=nc.now(),
               model=dict(file=str(args.model), sha256=nc.sha(args.model), name=model.get('name')),
               frames=len(rows), by_configuration_at_N2=summ,
               R_like_minus_X_like=summ.get('R', {}).get('mean'),
               TR0_around_its_conversion=conv,
               note='positive: the candidate raises these frames against the X-like frames; the terms as built put the '
                    'R-like placement Rx 3.38 kcal/mol too high against X in the small model (gate G-R)')
    path = pathlib.Path(args.out) if args.out else nc.N2 / 'behaviour' / pathlib.Path(args.model).stem / 'frames.json'
    nc.write_json(path, out)
    print(json.dumps({k: out[k] for k in ('frames', 'by_configuration_at_N2', 'R_like_minus_X_like')}, indent=1))
    print('TR0 around its conversion:', [(c['ps_after_release'], c['cls'], c['change']) for c in conv])


if __name__ == '__main__':
    main()
