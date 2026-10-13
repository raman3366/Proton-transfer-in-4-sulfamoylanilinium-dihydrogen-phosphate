"""Gate G-R2, a look that is NOT the gate: the force field as built and a candidate along the two QUANTUM-CHEMICAL
RELAXATIONS of gate G-R (md/N2.md, section 6.7).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_relaxpath.py --model MODEL.json [--out FILE]

The relaxations of the ligand with zinc and imidazoles fixed, from X and from the R-like placement Rx (tasks
gr_relax_X, gr_relax_Rx and gr_relax_Rx_cont of the rented machine; B3LYP/6-31G* with D3(BJ) in every step), leave
one geometry and one energy per step.  Copies of the four files are read from md/qm/n2/relax_paths/ (X_traj.npz,
X_steps.jsonl, Rx_traj.npz, Rx_steps.jsonl; copied from /root/cqc-sfa-md/qm/gr/relax/ of the machine, not changed).
These geometries are RELAXED, are neither fitting nor test points of gate G-R2, and were not seen by any fit: an
independent set, though not a designed one (two paths, one of them possibly unfinished).

For each step: E_QM (with D3) and the force-field energy of the same 64 atoms with the terms as built and with the
terms of MODEL.json, all against ONE reference, the last step of the relaxation from X.  Reported: the r.m.s. of
force field minus quantum chemistry along each path for both force fields, and the same with the mean of each path
removed (shape only); the values at the start of the run from Rx, at its flattest nitrogen and at its last step.
The two force fields differ only in the terms at N2 (and in the charges, when N2_FF names the second copy), so
which of them follows the paths better is a statement about these terms.
"""
import argparse
import json
import pathlib

import numpy as np

import n2_common as nc

SRC = nc.N2 / 'relax_paths'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    model = json.loads(pathlib.Path(args.model).read_text())
    if model.get('force_field', 'ff') != nc.FF_NAME:
        raise SystemExit(f"the model was fitted with the force-field copy {model.get('force_field')!r}, this run uses {nc.FF_NAME!r} (N2_FF)")
    sfx = '' if nc.FF_NAME == 'ff' else f'_{nc.FF_NAME}'
    out = pathlib.Path(args.out) if args.out else nc.N2 / f'behaviour{sfx}' / pathlib.Path(args.model).stem / 'relax_paths.json'
    work = out.parent / 'work_relaxpath'
    frc = out.parent / 'n2_terms_candidate.frcmod'
    out.parent.mkdir(parents=True, exist_ok=True)
    nc.write_frcmod(model['parameters'], frc, f"terms at N2, candidate {model.get('name', '')}; for a look, NOT ACCEPTED")
    ffs = dict(as_built=nc.FF(nc.N2 / 'ff' / 'small_as_built.prmtop'), candidate=nc.FF(nc.build_topology(frc, work)))
    if nc.FF_NAME != 'ff':                                                 # the new charges alone, terms at N2 as built
        ffs['new_charges_terms_as_built'] = nc.FF(nc.FFDIR / 'small_as_built.prmtop')
    _, x_X = nc.read_xyz(nc.GR / 'geom' / 'X.xyz')
    _, x_Rx = nc.read_xyz(nc.GR / 'geom' / 'Rx.xyz')
    paths = {}
    for start, x0 in (('X', x_X), ('Rx', x_Rx)):
        xyz = np.load(SRC / f'{start}_traj.npz')['xyz']
        steps = [json.loads(l) for l in (SRC / f'{start}_steps.jsonl').read_text().splitlines() if l.strip()]
        n = min(len(xyz), len(steps))
        assert float(np.abs(xyz[0] - x0).max()) < 2e-3, f'the first geometry of the run from {start} is not point {start}'
        paths[start] = dict(xyz=xyz[:n], steps=steps[:n],
                            files={f: nc.sha(SRC / f) for f in (f'{start}_traj.npz', f'{start}_steps.jsonl')})
    ref_x = paths['X']['xyz'][-1]
    e_ref_qm = paths['X']['steps'][-1]['e_total'] * nc.HARTREE
    e_ref = {k: ff.energy(ref_x)['total'] for k, ff in ffs.items()}
    doc = dict(what='force field as built and a candidate along the two quantum-chemical relaxations of gate G-R (zinc and '
                    'imidazoles fixed); a look on relaxed geometries that no fit has seen, not the gate; written by '
                    'md/qm/n2_relaxpath.py', written=nc.now(),
               model=dict(file=str(args.model), sha256=nc.sha(args.model), name=model.get('name')), force_field_copy=nc.FF_NAME,
               reference='the last step of the relaxation from X', paths={})
    for start, p in paths.items():
        rows = []
        for xyz, st in zip(p['xyz'], p['steps']):
            r = dict(step=st['step'], signed_volume=st['signed_volume'], Zn_N2_S_C4=st['Zn_N2_S_C4'],
                     grad_rms_ligand=st['grad_rms_ligand'], QM=round(st['e_total'] * nc.HARTREE - e_ref_qm, 3))
            for k, ff in ffs.items():
                r[k] = round(ff.energy(xyz)['total'] - e_ref[k], 3)
            rows.append(r)
        q = np.array([r['QM'] for r in rows])
        summ = dict(n_steps=len(rows), last_step=p['steps'][-1]['step'], last_time=p['steps'][-1]['time'],
                    gradient_rms_at_the_last_step=p['steps'][-1]['grad_rms_ligand'],
                    signed_volume_first_and_last=[rows[0]['signed_volume'], rows[-1]['signed_volume']],
                    QM_first_and_last=[rows[0]['QM'], rows[-1]['QM']])
        for k in ffs:
            d = np.array([r[k] for r in rows]) - q
            summ[k] = dict(rms_of_force_field_minus_QM=round(float(np.sqrt((d ** 2).mean())), 3), mean=round(float(d.mean()), 3),
                           rms_about_the_mean_of_the_path=round(float(np.sqrt(((d - d.mean()) ** 2).mean())), 3),
                           largest=round(float(d[np.argmax(np.abs(d))]), 3), at_step=rows[int(np.argmax(np.abs(d)))]['step'],
                           first_and_last=[rows[0][k], rows[-1][k]])
        if start == 'Rx':
            # the flattest nitrogen among the steps that are near a stationary point (gradient below 6e-4)
            near = [k for k, r in enumerate(rows) if r['grad_rms_ligand'] < 6e-4] or list(range(len(rows)))
            k_flat = min(near, key=lambda k: abs(rows[k]['signed_volume']))
            summ['at_the_flattest_nitrogen'] = rows[k_flat]
            summ['at_the_first_step'] = rows[0]
            summ['at_the_last_step'] = rows[-1]
        doc['paths'][start] = dict(files_sha256=p['files'], summary=summ, steps=rows)
    nc.write_json(out, doc)
    for start in ('X', 'Rx'):
        s = doc['paths'][start]['summary']
        print(f"run from {start}: {s['n_steps']} steps, signed volume {s['signed_volume_first_and_last']}, QM first and last "
              f"{s['QM_first_and_last']} (gradient at the last step {s['gradient_rms_at_the_last_step']:.1e})")
        for k in ffs:
            print(f"   {k:10s} minus QM: r.m.s. {s[k]['rms_of_force_field_minus_QM']}, mean {s[k]['mean']}, about the mean "
                  f"{s[k]['rms_about_the_mean_of_the_path']}, largest {s[k]['largest']} at step {s[k]['at_step']}; first and last {s[k]['first_and_last']}")
    s = doc['paths']['Rx']['summary']
    print('   run from Rx, first step:', s['at_the_first_step'])
    print('   run from Rx, flattest nitrogen:', s['at_the_flattest_nitrogen'])
    print('   run from Rx, last step:', s['at_the_last_step'])


if __name__ == '__main__':
    main()
