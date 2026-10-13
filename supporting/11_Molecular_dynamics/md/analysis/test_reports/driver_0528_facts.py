# throwaway (version 1.3): facts of the five test runs of the driver's maintainer that the text of section 16.4
# quotes, read from their manifests (read only) and from the summaries that analyse_md.py wrote for them in the
# scratch folder (run_all.sh, section "the test runs of the driver's maintainer").   usage: driver_0528_facts.py SCRATCH OUT.json
import json, sys
S, OUT = sys.argv[1], sys.argv[2]
HYD = '/root/cqc-sfa-md/test_driver/hydration/runs'
out = {'what': "the test runs of the driver's maintainer under analyse_md.py: steps, settings, warnings; and why a DCD header can differ from the stage table", 'runs': {}}
for s in ('dry', 'arrive', 'wet', 'skip', 'late'):
    man = json.load(open(f'{HYD}/{s}/manifest.json'))
    summ = json.load(open(f'{S}/driver_0528/{s}/summary.json'))
    st = {x['tag']: x for x in man['equilibration']['stages']}
    hs = man['hydration_check']['settings']
    w, e = summ['warnings'], summ['equilibration']
    dt4 = st['eq04']['settings']['timestep_fs']
    steps4 = st['eq04']['time_start_ps'] * 1000.0 / dt4
    out['runs'][s] = {
        'driver_script_sha256': man['script']['sha256'], 'analysis_script_sha256': summ['script']['sha256'], 'spec': summ['spec'],
        'created': man['created']['ist'], 'last_written': man['last_written']['ist'], 'equil_scale': man['equil_scale'],
        'settings_in_the_manifest': {'site_atoms': hs['site_atoms'], 'water_cutoff_A': hs['water_cutoff_A']},
        'settings_used_by_the_analysis': {k: summ['water_ring_place'][k] for k in ('site_atoms', 'water_cutoff_A', 'source', 'settings_are_the_defaults')},
        'full_frame_interval_ps': man['protocol']['output']['full_frame_interval_ps'],
        'length_ps_of_eq03': st['eq03']['settings']['length_ps'], 'length_ps_of_an_extension_step': hs['extension_length_ps_in_this_run'],
        'steps_in_the_manifest': list(st), 'steps_read': [a['tag'] for a in e['stages']],
        'n_extension_steps': sum(1 for t in st if 'x' in t),
        'frames_read_as_recorded': all([a['solute']['n_frames'], a['full']['n_frames']] == [st[a['tag']]['frames']['solute'], st[a['tag']]['frames']['full']] for a in e['stages']),
        'counts_of_the_driver': [a['driver_hydration_check']['n_water_O'] for a in e['stages'] if a['driver_hydration_check']],
        'compared': [a['driver_hydration_check']['compared'] for a in e['stages'] if a['driver_hydration_check']],
        'eq04_start_ps': st['eq04']['time_start_ps'], 'eq04_timestep_fs': dt4, 'eq04_start_in_its_time_steps': round(steps4, 6),
        'eq04_starts_on_a_whole_time_step': abs(steps4 - round(steps4)) < 1e-6,
        'header_first_frame_ps_of_solute_eq04': next(a['solute']['header_first_frame_ps'] for a in e['stages'] if a['tag'] == 'eq04'),
        'solute_interval_ps': summ['run']['solute_interval_ps'],
        'n_warnings': len(w), 'warnings_settings_not_those_of_the_plan': sum('hydration check of this run' in x for x in w),
        'warnings_header': sum('time axis of the equilibration is wrong' in x for x in w),
        'warnings_other': sum(('hydration check of this run' not in x) and ('time axis of the equilibration is wrong' not in x) for x in w)}
json.dump(out, open(OUT, 'w'), indent=1)
for s, r in out['runs'].items():
    print(s, r['created'], '| eq04 starts at', r['eq04_start_ps'], 'ps =', r['eq04_start_in_its_time_steps'], 'steps of', r['eq04_timestep_fs'], 'fs | header', r['header_first_frame_ps_of_solute_eq04'],
          '| header warnings', r['warnings_header'], '| eq03', r['length_ps_of_eq03'], 'ps, extension', r['length_ps_of_an_extension_step'], 'ps, full frames every', r['full_frame_interval_ps'], 'ps | compared', r['compared'])
