#!/usr/bin/env python
"""Gate G-R2: the verdict file md/qm/n2/n2_gate_result.json, put together from the result files of the examiner
(md/qm/n2/final_test/).  It computes nothing new; every number is read from a file named beside it.

    nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_result.py

If md/qm/n2/final_test/test_set_read.json exists (n2_exam_final.py --final was run on the user's decision), conditions
1 and 2 are taken from it; otherwise they are written as NOT EVALUATED and the test set as unread.
"""
import hashlib
import json
import pathlib
import time

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
N2 = REPO / 'md' / 'qm' / 'n2'
OUT = N2 / 'final_test'


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def load(name):
    f = OUT / name
    return json.loads(f.read_text()) if f.exists() else None


def record(beh):
    """What section 9.4 and 9.5 of md/N2.md report beside the conditions; read from the files named."""
    out = {}
    if beh:
        a, o = beh['along_the_quantum_chemical_relaxations'], beh['own_relaxations_zinc_and_imidazoles_fixed']
        out['relaxed_geometries'] = dict(
            file='md/qm/n2/final_test/behaviour.json', topologies=beh['topologies'],
            flattest_nitrogen_of_the_quantum_chemical_run_from_Rx=a['flattest_nitrogen_of_the_run_from_Rx'],
            quantum_chemical_run_from_Rx_last_step=dict(step=a['runs']['Rx']['last_step'], time=a['runs']['Rx']['last_time'], **a['runs']['Rx']['last']),
            rms_along_the_run_from_X=a['runs']['X']['rms_force_field_minus_QM'], rms_along_the_run_from_Rx=a['runs']['Rx']['rms_force_field_minus_QM'],
            own_relaxation_from_Rx={t: dict(R_like_minimum_of_its_own=v['R_like_placement_has_a_minimum_of_its_own'], signed_volume_at_the_end=v['from_Rx']['signed_volume'],
                                            end_minus_the_minimum_from_X_kcal=v['end_from_Rx_minus_end_from_X_kcal'],
                                            lowest_eigenvalue_of_the_ligand_hessian=v['from_Rx']['lowest_eigenvalues_of_the_ligand_hessian'][0])
                                    for t, v in o.items()})
    for name in ('second_step', 'charges_as_built'):
        b = load(f'built_system_{name}.json')
        if b:
            out[f'built_protein_systems_{name}'] = dict(
                file=f'md/qm/n2/final_test/built_system_{name}.json',
                starts={s: dict(atoms=len(v['comparison']) and None, differences_the_plan_does_not_allow=v['differences_the_plan_does_not_allow'],
                                angles=v['comparison']['angles']['n_different'], torsion_quadruples=v['comparison']['proper']['n_different'],
                                impropers=v['comparison']['improper']['n_different'], bonds=v['comparison']['bonds']['n_different'],
                                pairs_1_4=v['comparison']['pairs14']['n_different'], charges_changed=v['comparison']['charges_changed']['n'],
                                residues_with_changed_charges=v['comparison']['charges_changed']['residues'],
                                builder_checks_failed=v.get('checks_failed'), real_system_checksums=[v['real_system_checksums_before'], v['real_system_checksums_after']])
                        for s, v in b['starts'].items()})
            for v in out[f'built_protein_systems_{name}']['starts'].values():
                v.pop('atoms')
    for name in ('second_step', 'charges_as_built', 'as_built'):
        st = load(f'stability_{name}.json')
        if st:
            out[f'smoke_test_20_ps_{name}'] = dict(
                file=f'md/qm/n2/final_test/stability_{name}.json', what=st['what'],
                starts={s: {k: v.get(k) for k in ('done', 'finite', 'ps_at_4_fs', 'sign_changes_of_the_signed_volume', 'frames_with_signed_volume_below_0_3',
                                                  'frames_with_signed_volume_above_minus_0_3', 'frames_d1_above_6_7_A', 'frames_d2_below_5_2_A', 'frames')}
                        | dict(means={k: m['mean'] for k, m in (v.get('summary_of_the_4_fs_part') or {}).items()},
                               signed_volume_min_max=[(v.get('summary_of_the_4_fs_part') or {}).get('signed_volume', {}).get('min'),
                                                      (v.get('summary_of_the_4_fs_part') or {}).get('signed_volume', {}).get('max')])
                        for s, v in st['starts'].items()})
    return out


def main():
    cond, integ, mm, d3 = load('conditions_fitting_data.json'), load('integrity.json'), load('mm_energies.json'), load('d3_check.json')
    comp, beh, test = load('comparison_with_the_designer.json'), load('behaviour.json'), load('test_set_read.json')
    names = dict(second_step='candidate of the declared second step: md/qm/n2/n2_terms_new.frcmod (L5, ridge 0.01) with the charges fitted to two placements',
                 charges_as_built='candidate with the charges as built: md/qm/n2/n2_terms_new_charges_as_built.frcmod (L5, ridge 0.01)')
    cands, failed_all = {}, {}
    for key, c in cond['candidates'].items():
        c3, c4, c5, c6 = c['condition_3_part_A'], c['condition_4_bounds'], c['condition_5_leave_one_out'], c['condition_6_minimum']
        worst6 = max((abs(v) for m in c6['minima'].values() for v in m['difference'].values()))
        t = test['candidates'][key] if test else None
        conds = {
            '1 test set, r.m.s. against one reference': dict(
                limit='at most 1.0 kcal/mol, at least 20 test points', value=t['condition_1']['rms'] if t else None,
                met=t['condition_1']['holds'] if t else 'NOT EVALUATED: the test set was left unread', file='md/qm/n2/final_test/test_set_read.json' if t else None),
            '2 test set, largest single deviation': dict(
                limit='at most 2.5 kcal/mol', value=t['condition_2']['largest'] if t else None,
                met=t['condition_2']['holds'] if t else 'NOT EVALUATED: the test set was left unread', file='md/qm/n2/final_test/test_set_read.json' if t else None),
            '3 part A of gate G-R, 23 comparisons, two families': dict(
                limit='at most 1.0 kcal/mol', value=c3['rms'], single_reference_beside=c3['single_reference']['rms'], largest=c3['largest'],
                references=[c3['families']['X']['reference'], c3['families']['Rx']['reference']], as_built=c3['as_built'],
                from_leave_one_out_errors=c3['from_my_leave_one_out_errors'], in_sample=True, met=c3['holds'], file='md/qm/n2/final_test/conditions_fitting_data.json'),
            '4 every parameter within the bounds of point 2, nothing else changed': dict(
                limit='K 10 to 300; theta0 within 15 deg of the optimum; torsion and improper coefficients at most 5; periodicity 1 to 3',
                value=dict(angles=c4['angles'], torsion_terms=c4['torsion_terms'], improper=c4['improper'], outside_bounds=c4['outside_bounds'],
                           differences_the_plan_does_not_allow=c4['differences_the_plan_does_not_allow']),
                met=c4['holds'], file='md/qm/n2/final_test/conditions_fitting_data.json, topology_differences.json'),
            '5 leave-one-out error on the fitting set': dict(
                limit='at most 1.0 kcal/mol', value=c5['leave_one_out_rms'], largest=c5['largest'], largest_at=c5['largest_at'], residual_rms=c5['residual_rms'],
                n_fitting_rows=c['n_fitting_rows'], n_parameters=c['n_parameters'], against_one_reference=c5['against_one_reference'],
                the_delivered_file_is_this_fit=c5['the_delivered_file_is_this_fit'], met=c5['holds'], file='md/qm/n2/final_test/conditions_fitting_data.json'),
            '6 force-field minimum of the small model against its quantum-chemical optimum': dict(
                limit='within 5 deg in Zn-N2-S, Zn-N2-H and S-N2-H', value={k: v['difference'] for k, v in c6['minima'].items()},
                largest_abs_difference_deg=worst6, optimum=c6['optimum'], as_built={k: v['difference'] for k, v in c6['as_built'].items()},
                met=c6['holds'], file='md/qm/n2/final_test/conditions_fitting_data.json'),
        }
        failed = [k for k, v in conds.items() if v['met'] is False]
        failed_all[key] = failed
        cands[key] = dict(what=names[key], frcmod_sha256=c['frcmod_sha256'], conditions=conds, failed=failed,
                          not_evaluated=[k for k, v in conds.items() if isinstance(v['met'], str)],
                          accepted=bool(all(v['met'] is True for v in conds.values())), for_the_record=c['for_the_record'],
                          fitting_rows_in_sample=c['fitting_rows_built_topology'])
    doc = dict(
        what='gate G-R2 (md/PLAN.md, entry of 6 Oct 2026, 1:34 PM IST, point 4): verdict of the independent examiner',
        written=now(),
        verdict='NOT ACCEPTED',
        accepted=False,
        verdict_in_words='No model of the sulfonamide nitrogen is accepted; the force field stays as built. Both candidates miss conditions 5 and 6, '
                         'the candidate of the second step also condition 3. These three conditions use no test point, so the verdict does not depend on the test set.',
        test_set=('READ ONCE (md/qm/n2/final_test/test_set_read.json)' if test else
                  'NOT READ. 47 test points are computed and their result files have not been opened (rented machine and the watchdog copy on the Pi). '
                  'Conditions 1 and 2 are therefore not evaluated. Reason: the verdict is fixed without them, the builder put no candidate forward '
                  '(none passed the preview that the design requires before a final test), and a read cannot be taken back. '
                  'To read it: N2_FINAL_TEST=1 python md/qm/n2/final_test/n2_exam_final.py --final (once; a decision for the user).'),
        candidates=cands,
        every_candidate_of_the_ladder={k: dict(smallest_leave_one_out_rms=v['smallest_leave_one_out_rms'], any_at_most_1=v['any_at_most_1'],
                                               n=len(v['candidates']), largest_difference_from_the_builder=v['largest_difference_from_the_builder'])
                                       for k, v in cond['all_candidates_of_the_ladder'].items()},
        integrity=dict(file='md/qm/n2/final_test/integrity.json', split_sha256=integ['split']['sha256'], split_equal_to_the_recorded=integ['split']['equal'],
                       split_equal_on_the_machine=integ['split']['equal_on_the_machine'], split_drawn_again_from_seed=integ['split']['drawn_again_here_from_seed_and_rule_equals_the_file'],
                       split_written=integ['split']['files_written']['split.json'], first_new_single_point_ended=integ['split']['first_single_point_of_a_new_geometry_ended'],
                       points_finished=integ['batch']['finished'], fit=integ['batch']['fit'], test=integ['batch']['test'], all_converged_first_attempt=integ['batch']['all_at_the_first_attempt'],
                       fitting_files_with_other_settings=integ['settings']['fitting_files_that_differ_or_are_not_converged'],
                       test_files_in_the_repository=integ['repository']['flagged_test'], copies_of_test_results_elsewhere=integ['test_result_files']['files_elsewhere_with_the_bytes_of_a_test_result'],
                       fitting_points_left_out_of_the_fits={k: v['fitting_points_left_out'] for k, v in integ['repository']['fits_of_the_rule'].items()}),
        own_route=dict(file='md/qm/n2/final_test/mm_energies.json', route_A_against_route_B=mm['route_A_against_route_B'],
                       sander_against_openmm_largest=max(abs(v['difference']) for v in mm['sander_against_openmm'].values()),
                       as_built_against_the_table_of_gate_G_R=mm['as_built_against_the_table_of_gate_G_R'],
                       d3_gate_routine_against_s_dftd3_largest_hartree=d3['largest_abs_difference_hartree'],
                       d3_stored_against_recomputed_largest_hartree=cond['data_notes']['d3_stored_minus_recomputed_largest_hartree']),
        designer_gate_script_in_preview_mode={k: dict(verdict=v['designer_verdict'], failed=v['designer_failed'], largest_difference_from_the_examiner=v['largest_difference'])
                                              for k, v in comp['candidates'].items()} if comp else None,
        for_the_record_no_weight_in_the_verdict=record(beh),
        files={f.name: sha(f) for f in sorted(OUT.glob('*')) if f.is_file()},
        limits=cond['limits'])
    (N2 / 'n2_gate_result.json').write_text(json.dumps(doc, indent=1) + '\n')
    print('verdict:', doc['verdict'])
    for k, v in cands.items():
        print(k, '| failed:', v['failed'], '| not evaluated:', [x[:1] for x in v['not_evaluated']])
    print('test set:', doc['test_set'][:60])


if __name__ == '__main__':
    main()
