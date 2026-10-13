"""Gate G-R2, step 4: the files of terms at N2 in the form md/build_topology.py --n2-frcmod reads.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_deliver.py

Writes, from the first candidate of each fit of the rule (md/qm/n2_fit.py; no energy is read here):

    md/qm/n2/n2_terms_new.frcmod                      first candidate of the declared second step (md/qm/n2/fit_ff2,
                                                      charges fitted to two placements); to be built ONLY together
                                                      with --mcpb-dir md/qm/n2/charges2/mcpb_2e_two (absolute path)
    md/qm/n2/n2_terms_new_charges_as_built.frcmod     first candidate with the charges as built (md/qm/n2/fit)
    md/qm/n2/n2_terms_new.json                        which model each file is, checksums, the state, the options

THE STATE IS WRITTEN INTO THE FIRST LINE OF EACH FILE.  A candidate that the fit program does not expect to pass
the preview of the gate is marked "NOT ACCEPTED ... not for a run"; this script never writes the word "accepted"
alone: acceptance is the business of md/qm/n2_gate.py final, which writes its own file
(md/qm/n2/gate/n2_terms_accepted.frcmod).
"""
import json
import pathlib

import n2_common as nc

MCPB_TWO = nc.N2 / 'charges2' / 'mcpb_2e_two'


def one(fit_dir, out, charges):
    doc = json.loads((fit_dir / 'candidates.json').read_text())
    order = doc['order_of_preference']
    if not order:
        raise SystemExit(f'{fit_dir}: no eligible candidate')
    c = next(c for c in doc['candidates'] if c['name'] == order[0])
    mfile = fit_dir / f"MODEL_01_{c['name']}.json"
    model = json.loads(mfile.read_text())
    assert model['name'] == c['name'] and model['rank'] == 1 and model['n_points'] == c['n_points']
    exp = c['expected_in_the_preview_of_the_gate']
    ok = c['expected_to_pass_the_preview_without_condition_6']
    state = ('candidate for the final test; NOT ACCEPTED yet' if ok else
             'NOT ACCEPTED and not expected to pass the gate (leave-one-out '
             f"{c['leave_one_out_rms']:.2f} against 1.0): a record, not for a run")
    title = (f"terms at N2, gate G-R2 (md/N2.md): candidate {c['name']}, {charges}; {state}; "
             f"written by md/qm/n2_deliver.py {nc.now()}")
    nc.write_frcmod(model['parameters'], out, title)
    return dict(file=str(out.relative_to(nc.REPO)), sha256=nc.sha(out), candidate=c['name'], charges=charges,
                model_file=str(mfile.relative_to(nc.REPO)), model_sha256=nc.sha(mfile),
                fit=str((fit_dir / 'candidates.json').relative_to(nc.REPO)), fit_sha256=nc.sha(fit_dir / 'candidates.json'),
                n_fitting_rows=c['n_points'], n_parameters=c['n_parameters'], residual_rms=c['residual_rms'],
                leave_one_out_rms=c['leave_one_out_rms'], one_reference=c['one_reference'],
                part_A_in_sample=c['part_A_in_sample']['rms'],
                part_A_from_leave_one_out_errors=c['part_A_from_leave_one_out_errors']['rms'],
                expected_in_the_preview_of_the_gate=exp, expected_to_pass_the_preview_without_condition_6=ok,
                any_candidate_of_this_fit_expected_to_pass=doc['any_candidate_expected_to_pass'], state=state,
                amber_as_written=c['amber_as_written'])


def main():
    two = one(nc.N2 / 'fit_ff2', nc.N2 / 'n2_terms_new.frcmod', 'charges fitted to two placements (md/qm/n2/charges2/mcpb_2e_two)')
    built = one(nc.N2 / 'fit', nc.N2 / 'n2_terms_new_charges_as_built.frcmod', 'charges as built')
    mol2 = {f.name: nc.sha(f) for f in sorted(MCPB_TWO.glob('*.mol2'))}
    rec = dict(
        what='the files of terms at N2 for md/build_topology.py --n2-frcmod; written by md/qm/n2_deliver.py; a file '
             'is usable for a run ONLY after md/qm/n2_gate.py final has accepted its model', written=nc.now(),
        model_ready_for_the_final_test=bool(two['expected_to_pass_the_preview_without_condition_6']
                                            or built['expected_to_pass_the_preview_without_condition_6']),
        second_step=two, charges_as_built=built,
        mol2_set_of_the_second_step=dict(folder=str(MCPB_TWO.relative_to(nc.REPO)), sha256=mol2,
                                         note='the five residues of the zinc site; ALL of them carry new charges, not the ligand alone'),
        builder_options=dict(
            second_step='md/build_topology.py --model bonded --set 2e --start <X|N|R> --n2-frcmod md/qm/n2/n2_terms_new.frcmod '
                        f'--mcpb-dir {MCPB_TWO}',
            charges_as_built='md/build_topology.py --model bonded --set 2e --start <X|N|R> --n2-frcmod '
                             'md/qm/n2/n2_terms_new_charges_as_built.frcmod',
            note='--mcpb-dir must be an absolute path; the options are a record, NOT to be run for a real system '
                 'while no model is accepted'))
    nc.write_json(nc.N2 / 'n2_terms_new.json', rec)
    for r in (two, built):
        print(r['file'], r['sha256'][:12], '|', r['state'])
    print('model ready for the final test:', rec['model_ready_for_the_final_test'])


if __name__ == '__main__':
    main()
