"""Gate G-R2, diagnosis from the data of the failed gate G-R (md/N2.md, section 1).  No new energy.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_diagnosis.py

Reads md/qm/gr/gr_table.csv (38 points: quantum chemistry with and without D3, force field split by part),
md/qm/n2/points.json (features and `rest` of the 36 rigid points) and the optimised small model
(/root/cqc-sfa-md/qm/prod/opt_final.npy).  Writes
  md/qm/n2/qm_optimum.json  the quantum-chemical optimum pinned for the bounds and for condition 6 of the gate;
  md/qm/n2/diagnosis.json   (a) per family of points: the deviation force field minus quantum chemistry against
                            the base of the family, and the change of every part of the force-field energy;
                            (b) what the allowed terms can absorb on the 36 old points: the ladder of models
                            fitted to them, without and with the bounds, residual and leave-one-out error, in all
                            and per family.  36 points for up to 17 parameters: this is a diagnosis, not a model.
"""
import csv
import json
import pathlib

import numpy as np

import n2_common as nc

OPT = pathlib.Path('/root/cqc-sfa-md/qm/prod/opt_final.npy')
PARTS = ('MM_bond', 'MM_angle_N2', 'MM_angle_other', 'MM_torsion_N2S', 'MM_torsion_other', 'MM_nb_ligand',
         'MM_lj_ligand_frame', 'MM_coulomb_ligand_frame')
FAMILIES = (('a: inversion path a2-a6', lambda r: r['group'] == 'a'),
            ('d: hydrogen moved on X', lambda r: r['group'] == 'd'),
            ('g: hydrogen moved on Rx', lambda r: r['group'] == 'g'),
            ('e: turns about Zn-N2', lambda r: r['group'] == 'e'),
            ('f: bend of Zn-N2-S', lambda r: r['group'] == 'f'),
            ('c1: turns about N2-S from X', lambda r: r['id'].startswith('c1_')),
            ('c2: turns about N2-S from Rx', lambda r: r['id'].startswith('c2_')),
            ('Rx and Nx against X', lambda r: r['id'] in ('Rx', 'Nx')))


def main():
    g = nc.Geometry()
    m = g.m
    # ---- the quantum-chemical optimum, pinned
    x = np.load(OPT)
    assert x.shape == (64, 3)
    zn, lig = x[45], x[46:]
    i = m.i
    a = dict(Zn_N2_S=nc.angle(zn, lig[i['N2']], lig[i['S']]), Zn_N2_H=nc.angle(zn, lig[i['N2']], lig[i['HN21']]),
             S_N2_H=nc.angle(lig[i['S']], lig[i['N2']], lig[i['HN21']]))
    z, h, u = (nc.unit(v - lig[i['N2']]) for v in (zn, lig[i['HN21']], lig[i['S']]))
    nc.write_json(nc.N2 / 'qm_optimum.json', dict(
        what='optimised small model [Zn(ethyl-imidazole)3(sulfonamidate)]+, all atoms free, B3LYP/6-31G* without a '
             'dispersion term (task qm_opt of the rented machine, md/qm/run_qm.py); pinned by md/qm/n2_diagnosis.py '
             'for the bounds on the equilibrium angles (15 deg) and for condition 6 of gate G-R2 (5 deg)',
        source=str(OPT), source_sha256=nc.sha(OPT), written=nc.now(),
        angles_at_N2={k: round(v, 4) for k, v in a.items()}, angle_sum=round(sum(a.values()), 3),
        signed_volume=round(float(np.dot(z, np.cross(h, u))), 4),
        Zn_N2=round(float(np.linalg.norm(zn - lig[i['N2']])), 4), N2_S=round(float(np.linalg.norm(lig[i['S']] - lig[i['N2']])), 4),
        C4_S_N2=round(nc.angle(lig[i['C4']], lig[i['S']], lig[i['N2']]), 3),
        Zn_N2_S_C4=round(nc.dihedral(zn, lig[i['N2']], lig[i['S']], lig[i['C4']]), 3),
        windows_15deg={k: [round(v - 15, 3), round(v + 15, 3)] for k, v in a.items()},
        as_built_theta0=dict(Zn_N2_S=117.55, Zn_N2_H=107.406, S_N2_H=118.059),
        xyz_angstrom=[[round(float(t), 8) for t in v] for v in x]))
    rec, th_opt = nc.qm_optimum()

    rows = list(csv.DictReader(open(nc.GR / 'gr_table.csv')))
    by = {r['id']: r for r in rows}
    fam = {}
    for name, sel in FAMILIES:
        pts = [r for r in rows if sel(r)]
        tab = []
        for r in pts:
            b = by[r['base']] if name != 'Rx and Nx against X' else by['X']
            dev = (float(r['E_MM']) - float(b['E_MM'])) - (float(r['E_QM']) - float(b['E_QM']))
            tab.append(dict(id=r['id'], base=b['id'], QM=round(float(r['E_QM']) - float(b['E_QM']), 3),
                            MM=round(float(r['E_MM']) - float(b['E_MM']), 3), MM_minus_QM=round(dev, 3),
                            QM_without_D3=round(float(r['E_QM_no_D3']) - float(b['E_QM_no_D3']), 3),
                            parts={p[3:]: round(float(r[p]) - float(b[p]), 3) for p in PARTS}))
        d = np.array([t['MM_minus_QM'] for t in tab])
        fam[name] = dict(n=len(tab), rms_MM_minus_QM=round(float(np.sqrt(np.mean(d ** 2))), 3),
                         mean_MM_minus_QM=round(float(d.mean()), 3),
                         rms_of_each_part={p[3:]: round(float(np.sqrt(np.mean([t['parts'][p[3:]] ** 2 for t in tab]))), 3)
                                           for p in PARTS},
                         rms_MM=round(float(np.sqrt(np.mean([t['MM'] ** 2 for t in tab]))), 3),
                         rms_QM=round(float(np.sqrt(np.mean([t['QM'] ** 2 for t in tab]))), 3), points=tab)

    # ---- (b) what the allowed terms can absorb on the old points
    pj = json.loads((nc.N2 / 'points.json').read_text())
    old = pj['old_rigid_points']
    data = []
    for pid, o in old.items():
        e_qm = float(by[pid]['E_QM'])                                      # kcal/mol relative to X, with D3
        data.append(dict(id=pid, features=o['features'], y=e_qm - o['mm_as_built']['rest'],
                         group=by[pid]['group'], dev_as_built=float(by[pid]['E_MM']) - e_qm))
    d0 = np.array([r['dev_as_built'] for r in data])
    ladder = {}
    for name, groups in nc.LADDER:
        out = {}
        for bounded in (False, True):
            f = nc.fit_model(data, groups, 0.0, bounded, th_opt)
            lo = nc.leave_one_out(data, groups, 0.0, bounded, th_opt)
            res = np.array(f['residuals'])
            out['bounded' if bounded else 'free'] = dict(
                n_parameters=f['n_parameters'], points_per_parameter=round(len(data) / f['n_parameters'], 2),
                rms=round(f['rms'], 3), largest=round(f['largest'], 3), leave_one_out_rms=round(lo['rms'], 3),
                leave_one_out_largest=round(lo['largest'], 3), solver=f['solver'],
                amber=f['amber'], outside_bounds=nc.within_bounds(f['amber'], th_opt),
                residual_rms_by_group={gname: round(float(np.sqrt(np.mean(res[[k for k, r in enumerate(data) if r['group'] == gname]] ** 2))), 3)
                                       for gname in sorted({r['group'] for r in data})},
                residual_Rx_minus_X=round(float(res[[r['id'] for r in data].index('Rx')] - res[[r['id'] for r in data].index('X')]), 3))
        ladder[name] = dict(groups=list(groups), **out)
    doc = dict(what='diagnosis of the failed gate G-R for gate G-R2; written by md/qm/n2_diagnosis.py', written=nc.now(),
               source=dict(table='md/qm/gr/gr_table.csv', sha256=nc.sha(nc.GR / 'gr_table.csv')),
               families=fam,
               old_points_single_constant=dict(n=len(data), rms_as_built_about_its_mean=round(float(np.sqrt(np.mean((d0 - d0.mean()) ** 2))), 3),
                                               rms_as_built_against_X=round(float(np.sqrt(np.mean(d0 ** 2))), 3)),
               ladder_on_old_points=ladder, qm_optimum=rec['angles_at_N2'])
    nc.write_json(nc.N2 / 'diagnosis.json', doc)
    print('optimum:', rec['angles_at_N2'], 'sum', rec['angle_sum'], 'sv', rec['signed_volume'])
    for name, v in fam.items():
        print(f"{name:32s} n {v['n']:2d}  rms dev {v['rms_MM_minus_QM']:5.2f}  mean {v['mean_MM_minus_QM']:+5.2f}  parts rms {v['rms_of_each_part']}")
    print('as built, 36 points, one constant:', doc['old_points_single_constant'])
    for name, v in ladder.items():
        for b in ('free', 'bounded'):
            w = v[b]
            print(f"{name} {b:8s} P {w['n_parameters']:2d}  rms {w['rms']:5.2f}  LOO {w['leave_one_out_rms']:5.2f}  Rx-X {w['residual_Rx_minus_X']:+5.2f}  "
                  f"by group {w['residual_rms_by_group']}  out of bounds {len(w['outside_bounds'])}  {w['solver'].get('success')}")
    print(json.dumps(ladder['L3']['bounded']['amber']))


if __name__ == '__main__':
    main()
