# throwaway (versions 1.2 and 1.3): the analysis against the monitor values that md/run_md.py wrote for its own stand-in test
# run (a copy of /root/cqc-sfa-md/test_driver/runs/ref and systems/standin taken on 6 Oct 2026): the last frame of
# every production segment and, new, the last frame of every equilibration stage (timeseries_equil_solute.csv).
# Version 1.3, second part: two short runs of the driver as it stands now, with the hydration block of md/protocol.json
# as it stands now (six ring carbons, 2.5 A), on the stand-in systems `dry` and `wet` of the driver's own tests
# (driver_now/, made by run_driver.sh): the steps (with the extension steps eq03x1 to eq03x4 of `dry`), their times and
# frames, the driver's hydration records and its count of water in the place of the ring at the end of EVERY step
# (monitor, key crystal_ring_site) against n_water_ring_place and d_water_ring_place_min_A of the analysis.
import csv, hashlib, json, sys
S = sys.argv[1]
pairs = (('zn_to_lig_N2_A', 'd_Zn_N2_A'), ('lig_N2_to_thr199_OG1_A', 'd_N2_T199OG1_A'), ('lig_N1_to_phe131_CZ_A', 'd1_N1_F131CZ_A'),
         ('lig_N1_to_asn67_OD1_A', 'd2_N1_N67OD1_A'), ('zn_to_his94_NE2_A', 'd_Zn_H94NE2_A'), ('zn_to_his96_NE2_A', 'd_Zn_H96NE2_A'),
         ('zn_to_his119_ND1_A', 'd_Zn_H119ND1_A'))
man = json.load(open(f'{S}/driver_snapshot/runs/standin/ref/manifest.json'))
ts = list(csv.DictReader(open(f'{S}/driver_snapshot/analysis/timeseries_solute.csv')))
te = list(csv.DictReader(open(f'{S}/driver_snapshot/analysis/timeseries_equil_solute.csv')))
tef = list(csv.DictReader(open(f'{S}/driver_snapshot/analysis/timeseries_equil_full.csv')))
summ = json.load(open(f'{S}/driver_snapshot/analysis/summary.json'))
rows, worst, n = [], 0.0, 0
for seg in man['production']['segments']:
    n += seg['frames']['solute']
    fr = ts[n - 1]
    for a, b in pairs:
        d = abs(float(seg['monitor'][a]) - float(fr[b]))
        worst = max(worst, d)
        rows.append({'step': seg['tag'], 'frame': int(fr['frame']), 'driver_key': a, 'driver_value_A': seg['monitor'][a],
                     'analysis_column': b, 'analysis_value_A': float(fr[b]), 'abs_diff_A': round(d, 4)})
n_prod = len(rows)
eq_total = man['equilibration']['stages'][-1]['time_end_ps']
times_ok, n = True, 0
for st in man['equilibration']['stages']:                     # the monitor is taken at the end of the stage = its last frame
    n += st['frames']['solute']
    fr = te[n - 1]
    times_ok &= abs(float(fr['time_ns']) * 1000.0 - (st['time_end_ps'] - eq_total)) < 1e-6 and int(fr['stage']) == int(st['tag'][2:])
    for a, b in pairs:
        d = abs(float(st['monitor'][a]) - float(fr[b]))
        worst = max(worst, d)
        rows.append({'step': st['tag'], 'frame': int(fr['frame']), 'driver_key': a, 'driver_value_A': st['monitor'][a],
                     'analysis_column': b, 'analysis_value_A': float(fr[b]), 'abs_diff_A': round(d, 4)})
E = summ['equilibration']
out = {'what': 'last frame of each production segment and of each equilibration stage: monitor of md/run_md.py (three decimals) against analyse_md.py',
       'driver_script_sha256': man['script']['sha256'], 'analysis_script_sha256': summ['script']['sha256'],
       'n_compared': len(rows), 'n_compared_production': n_prod, 'n_compared_equilibration': len(rows) - n_prod,
       'largest_abs_diff_A': round(worst, 4), 'tolerance_A': 0.001, 'all_agree': worst <= 0.001,
       'n_frames_solute': summ['run']['n_frames_solute'], 'n_frames_full': summ['run']['n_frames_full'],
       'n_atoms_full': man['system']['n_atoms'], 'n_atoms_solute': summ['run']['n_atoms_solute'],
       'equilibration': {'available': E['available'], 'stage_table_source': E['stage_table_source'], 'total_ps': E['total_ps'],
                         'n_frames_solute': E['n_frames_solute'], 'n_frames_full': E['n_frames_full'],
                         'frames_per_stage': [[a['solute']['n_frames'], a['solute']['n_frames_expected'], a['full']['n_frames'], a['full']['n_frames_expected']] for a in E['stages']],
                         'manifest_frames_per_stage': [[s['frames']['solute'], s['frames']['full']] for s in man['equilibration']['stages']],
                         'header_first_frame_ps': [a['solute']['header_first_frame_ps'] for a in E['stages']],
                         'manifest_stage_start_ps': [s['time_start_ps'] for s in man['equilibration']['stages']],
                         'solute_interval_ps': summ['run']['solute_interval_ps'],
                         'time_first_frame_ns': float(te[0]['time_ns']), 'time_last_frame_ns': float(te[-1]['time_ns']),
                         'time_of_the_last_frame_of_every_stage_equals_the_stage_end_of_the_manifest': bool(times_ok),
                         'lines_of_timeseries_equil_full': len(tef)},
       'checks': {k: summ['run'][k] for k in ('header_check_solute', 'header_check_full', 'manifest_segment_check', 'segment_overlap_check', 'consistency_full_vs_solute', 'sources', 'imaging_solute', 'full_trajectory')},
       'warnings': summ['warnings'], 'rows': rows}
# ---- version 1.3: the driver of today with its hydration check
hyd = {'what': 'two short runs of md/run_md.py as it stands, hydration block of md/protocol.json as it stands, on the stand-in systems of '
               'the driver tests: steps, times, frames, hydration records and the count at the end of every step against analyse_md.py',
       'driver_script_sha256_now': hashlib.sha256(open(sys.argv[2] + '/run_md.py', 'rb').read()).hexdigest() if len(sys.argv) > 2 else None,
       'protocol_sha256_now': hashlib.sha256(open(sys.argv[2] + '/protocol.json', 'rb').read()).hexdigest() if len(sys.argv) > 2 else None,
       'runs': {}}
worst_h, n_cnt, n_cnt_same, ok_h = 0.0, 0, 0, True
for s in ('dry', 'wet'):
    man = json.load(open(f'{S}/driver_now/runs/{s}/R1/manifest.json'))
    summ = json.load(open(f'{S}/driver_now/analysis/{s}/summary.json'))
    tef = list(csv.DictReader(open(f'{S}/driver_now/analysis/{s}/timeseries_equil_full.csv')))
    tes = list(csv.DictReader(open(f'{S}/driver_now/analysis/{s}/timeseries_equil_solute.csv')))
    tf = list(csv.DictReader(open(f'{S}/driver_now/analysis/{s}/timeseries_full.csv')))
    E, W = summ['equilibration'], summ['water_ring_place']
    st_m = man['equilibration']['stages']
    eq_total = st_m[-1]['time_end_ps']
    hs = man['hydration_check']['settings']
    key_n = next(k for k in st_m[0]['monitor']['crystal_ring_site'] if k.startswith('n_water_O_within_'))
    steps, n = [], 0
    for st, a in zip(st_m, E['stages']):
        n += st['frames']['full']
        fr = tef[n - 1]                                       # the last full frame of the step
        mon = st['monitor']['crystal_ring_site']
        d = abs(float(mon['nearest_water_O_A']) - float(fr['d_water_ring_place_min_A']))
        worst_h = max(worst_h, d)
        n_cnt += 1
        n_cnt_same += int(mon[key_n]) == int(fr['n_water_ring_place'])
        hc, dh = st.get('hydration_check'), a['driver_hydration_check']
        ext = int(st['tag'].split('x')[1]) if 'x' in st['tag'] else 0
        steps.append({'tag': st['tag'], 'tag_read': a['tag'], 'stage_and_extension_read': [a['stage'], a['extension']],
                      'stage_and_extension_of_the_tag': [int(st['tag'][2:4]), ext],
                      'time_start_end_ps_manifest': [st['time_start_ps'], st['time_end_ps']], 'time_start_end_ps_read': [a['time_start_ps'], a['time_end_ps']],
                      'frames_manifest': [st['frames']['solute'], st['frames']['full']],
                      'frames_read': [a['solute']['n_frames'], a['full']['n_frames']],
                      'frames_expected_by_the_analysis': [a['solute']['n_frames_expected'], a['full']['n_frames_expected']],
                      'last_full_frame_at_the_end_of_the_step': abs(float(fr['time_ns']) * 1000.0 - (st['time_end_ps'] - eq_total)) < 1e-6,
                      'driver_monitor_count_and_nearest_A': [mon[key_n], mon['nearest_water_O_A']],
                      'analysis_count_and_nearest_A': [int(fr['n_water_ring_place']), float(fr['d_water_ring_place_min_A'])],
                      'abs_diff_nearest_A': round(d, 4),
                      'hydration_record_of_the_driver': None if hc is None else {k: hc[k] for k in ('n_water_O', 'applies', 'decision', 'water_cutoff_A', 'site_atoms', 'next_step')},
                      'analysis_driver_hydration_check': dh})
        x = steps[-1]
        ok_h &= (x['tag'] == x['tag_read'] and x['stage_and_extension_read'] == x['stage_and_extension_of_the_tag']
                 and x['time_start_end_ps_manifest'] == x['time_start_end_ps_read'] and x['frames_manifest'] == x['frames_read'] == x['frames_expected_by_the_analysis']
                 and x['last_full_frame_at_the_end_of_the_step'] and int(mon[key_n]) == int(fr['n_water_ring_place']) and d <= 0.001
                 and ((hc is None) == (dh is None))
                 and (hc is None or (dh['compared'] and dh['count_agrees'] and dh['nearest_agrees'] and dh['n_water_O'] == hc['n_water_O']
                                     and dh['decision'] == hc['decision'] and dh['applies'] == hc['applies'])))
    prod, n = [], 0
    for seg in man['production']['segments']:
        n += seg['frames']['full']
        fr, mon = tf[n - 1], seg['monitor']['crystal_ring_site']
        d = abs(float(mon['nearest_water_O_A']) - float(fr['d_water_ring_place_min_A']))
        worst_h = max(worst_h, d)
        n_cnt += 1
        n_cnt_same += int(mon[key_n]) == int(fr['n_water_ring_place'])
        prod.append({'tag': seg['tag'], 'driver_monitor_count_and_nearest_A': [mon[key_n], mon['nearest_water_O_A']],
                     'analysis_count_and_nearest_A': [int(fr['n_water_ring_place']), float(fr['d_water_ring_place_min_A'])], 'abs_diff_nearest_A': round(d, 4)})
        ok_h &= int(mon[key_n]) == int(fr['n_water_ring_place']) and d <= 0.001
    settings_ok = W['site_atoms'] == hs['site_atoms'] and abs(W['water_cutoff_A'] - hs['water_cutoff_A']) < 1e-9 and W['source'] == 'manifest hydration_check.settings'
    ok_h &= settings_ok and not summ['warnings'] and [x['tag'] for x in steps] == [st['tag'] for st in st_m] and len(E['stages']) == len(st_m)
    hyd['runs'][s] = {'driver_script_sha256': man['script']['sha256'], 'analysis_script_sha256': summ['script']['sha256'],
                      'protocol_of_the_run_sha256': man['inputs']['protocol']['sha256'], 'system_dir': man['system_dir'], 'run_id': man['run_id'],
                      'equil_scale': man['equil_scale'], 'test_only': man['test_only'],
                      'n_atoms_full': man['system']['n_atoms'], 'n_atoms_solute': summ['run']['n_atoms_solute'], 'n_CA': summ['run']['n_CA'],
                      'n_water_O_listed': summ['run']['full_trajectory']['n_water_O'],
                      'output_intervals_ps': [man['protocol']['output']['solute_frame_interval_ps'], man['protocol']['output']['full_frame_interval_ps']],
                      'settings_in_the_manifest': {k: hs[k] for k in ('after_step', 'site_atoms', 'water_cutoff_A', 'max_extensions')},
                      'settings_used_by_the_analysis': {k: W[k] for k in ('site_atoms', 'water_cutoff_A', 'source', 'settings_are_the_defaults')},
                      'settings_agree': bool(settings_ok),
                      'n_extensions_manifest': man['hydration_check']['n_extensions'], 'final_decision_manifest': man['hydration_check']['final_decision'],
                      'steps_in_the_manifest': [st['tag'] for st in st_m], 'steps_read': [a['tag'] for a in E['stages']],
                      'stage_table_source': E['stage_table_source'], 'total_ps': E['total_ps'],
                      'n_frames_equilibration_solute_full': [E['n_frames_solute'], E['n_frames_full']],
                      'time_first_and_last_equilibration_frame_ns': [float(tes[0]['time_ns']), float(tes[-1]['time_ns'])],
                      'header_first_frame_ps_solute': [a['solute']['header_first_frame_ps'] for a in E['stages']],
                      'steps': steps, 'production_segments': prod, 'warnings': summ['warnings']}
hyd.update(n_counts_compared=n_cnt, n_counts_identical=n_cnt_same, largest_abs_diff_nearest_A=round(worst_h, 4), tolerance_A=0.001,
           n_hydration_records=sum(1 for r in hyd['runs'].values() for x in r['steps'] if x['hydration_record_of_the_driver'] is not None),
           n_extension_steps_read=sum(1 for r in hyd['runs'].values() for x in r['steps'] if x['stage_and_extension_read'][1] > 0),
           all_agree=bool(ok_h))
out['hydration_runs_of_the_driver_of_today'] = hyd
json.dump(out, open(f'{S}/driver_snapshot/driver_crosscheck.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in hyd.items() if k != 'runs'}))
print(json.dumps({k: v for k, v in out.items() if k not in ('rows', 'checks')}))
