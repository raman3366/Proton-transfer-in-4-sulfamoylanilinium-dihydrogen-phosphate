# throwaway (version 1.3): the synthetic runs of version 1.2 under version 1.3.  The generator must write the same
# coordinates (production and equilibration files) and the same planted values in the columns of version 1.2; the
# analysis must write the same values in every column and under every key that version 1.2 had.  Two columns and two
# entries of the table row have a new name in version 1.3 (the map RENAMED below); nothing else may differ, except
# what the specification changes: the driver's count is compared with the new quantity, and run ext/R1, whose
# manifest now names the settings of its hydration check, gets one warning.
# usage: compare_v12_v13.py OLD_TEST NEW_TEST OUT.json
import csv, json, sys
from pathlib import Path
import numpy as np
from MDAnalysis.lib.formats.libdcd import DCDFile
OLD, NEW, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
RENAMED = {'n_water_ring_site_3p5A': 'n_water_ring_wide_3p5A', 'd_water_ring_site_min_A': 'd_water_ring_wide_min_A',
           'n_water_ring_site_3p5A_mean': 'n_water_ring_wide_3p5A_mean', 'd_water_ring_site_min_min_A': 'd_water_ring_wide_min_min_A'}
EXPECTED_WARNINGS = {'runs_v12/ext/R1': 1}          # one more than in version 1.2: the settings of its manifest are not the defaults
runs = sorted(p.parent.relative_to(OLD) for p in OLD.glob('runs*/*/*/solute.pdb'))
res = {'old_test': str(OLD), 'new_test': str(NEW), 'renamed': RENAMED, 'runs': {}, 'new_columns_solute': None, 'new_columns_full': None,
       'new_keys_of_table_row': None, 'new_keys_of_summary': None, 'new_keys_of_last_full_frame': None}
all_ok = True


def rows(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def new_name(c):
    return RENAMED.get(c, c)


def same_rows(a, b):
    return len(a) == len(b) and all(x[c] == y[new_name(c)] for x, y in zip(a, b) for c in x)


def old_view(e):
    """The block `equilibration` of version 1.3 with the keys that version 1.2 had (last_full_frame cut down to them)."""
    e = json.loads(json.dumps(e))
    for st in e.get('stages', []):
        lf = st.get('last_full_frame')
        if lf is not None:
            st['last_full_frame'] = {'frame': lf['frame'], 'time_ns': lf['time_ns'], 'n_water_ring_site_3p5A': lf['n_water_ring_wide_3p5A'],
                                     'd_water_ring_site_min_A': lf['d_water_ring_wide_min_A']}
    return e


for r in runs:
    rec = {}
    name = '/'.join(r.parts)
    worst, nseg, nfr = 0.0, 0, 0
    for p in sorted(list((OLD / r).glob('*_seg*.dcd')) + list((OLD / r).glob('*_eq*.dcd'))):      # every file, every frame
        q = NEW / r / p.name
        if p.stat().st_size == 0 or q.stat().st_size == 0:
            rec.setdefault('frame_counts_equal', True)
            rec['frame_counts_equal'] &= p.stat().st_size == q.stat().st_size
            continue
        try:
            with DCDFile(str(p)) as a, DCDFile(str(q)) as b:
                same_n = len(a) == len(b)
                if len(a) and same_n:
                    xa, xb = a.readframes().xyz, b.readframes().xyz
                    worst = max(worst, float(np.abs(xa - xb).max()))
                    nfr += len(a)
                rec.setdefault('frame_counts_equal', True)
                rec['frame_counts_equal'] &= same_n
        except OSError:                                       # the file with half a frame at its end: compared byte by byte
            rec.setdefault('frame_counts_equal', True)
            rec['frame_counts_equal'] &= p.stat().st_size == q.stat().st_size
            ba, bb = p.read_bytes()[276:], q.read_bytes()[276:]
            worst = max(worst, 0.0 if ba == bb else float('inf'))
        nseg += 1
    rec.update(files_compared=nseg, frames_compared=nfr, largest_coordinate_difference_A=worst)
    tag = '__'.join(r.parts)
    for kind in ('solute', 'full', 'rmsf', 'equil_solute', 'equil_full'):
        pa, pb = OLD / 'planted' / f'{tag}__{kind}.csv', NEW / 'planted' / f'{tag}__{kind}.csv'
        if pa.is_file():
            rec[f'planted_{kind}_identical_in_the_columns_of_v12'] = same_rows(rows(pa), rows(pb))
    ad, bd = OLD / r / 'analysis', NEW / r / 'analysis'
    for fname in ('timeseries_solute.csv', 'timeseries_full.csv', 'rmsf_ca.csv', 'timeseries_equil_solute.csv', 'timeseries_equil_full.csv'):
        if not (ad / fname).is_file():
            rec[fname + ': written by neither version'] = not (bd / fname).is_file()
            continue
        a, b = rows(ad / fname), rows(bd / fname)
        rec[fname + ': lines'] = [len(a), len(b)]
        rec[fname + ': identical in the columns of v12'] = same_rows(a, b)
        if a and b:
            new = [c for c in b[0] if c not in [new_name(x) for x in a[0]]]
            key = {'timeseries_solute.csv': 'new_columns_solute', 'timeseries_full.csv': 'new_columns_full'}.get(fname)
            if key:
                assert res[key] in (None, new), (r, new)
                res[key] = new
            rec[fname + ': order of the columns of v12 kept'] = [c for c in b[0] if c in [new_name(x) for x in a[0]]] == [new_name(x) for x in a[0]]
    sa, sb = json.loads((ad / 'summary.json').read_text()), json.loads((bd / 'summary.json').read_text())
    same = {}
    for k in ('lig_rmsd_labelling_counts', 'asn67_ND2', 'validation', 'rmsf', 'orientation', 'gate_GV', 'n2_configuration', 'hbond_first_passage',
              'start_of_production', 'fractions'):
        same[k] = sa[k] == sb[k]
    same['stats'] = all(sb['stats'].get(new_name(c)) == v for c, v in sa['stats'].items())
    row_b = dict(sb['table_row'])
    exp_w = len(sa['warnings']) + EXPECTED_WARNINGS.get(name, 0)
    same['warnings (number; ext/R1 one more by the specification)'] = len(sb['warnings']) == exp_w
    row_a = dict(sa['table_row'], n_warnings=exp_w)
    same['table_row'] = all(row_b.get(new_name(c)) == v for c, v in row_a.items())
    same['equilibration'] = old_view(sb['equilibration']) == sa['equilibration']
    same['run'] = all(sa['run'][c] == sb['run'][c] for c in ('role', 'start', 'replica', 'model', 'planned_ns', 'seen_ns', 'complete', 'n_frames_solute',
                                                               'n_frames_full', 'imaging_solute', 'header_check_solute', 'manifest_segment_check',
                                                               'segment_overlap_check', 'consistency_full_vs_solute', 'sources'))
    old_const = {k: v for k, v in sa['constants'].items() if k not in ('ring_site_atoms', 'ring_site_water_distance_max_A')}
    same['constants (those of v12 but the two renamed ones)'] = (all(sb['constants'].get(c) == v for c, v in old_const.items())
                                                                 and sb['constants']['ring_wide_atoms'] == sa['constants']['ring_site_atoms']
                                                                 and sb['constants']['ring_wide_water_distance_max_A'] == sa['constants']['ring_site_water_distance_max_A'])
    rec['summary.json, keys of v12 with the same content'] = same
    lfk = sorted({k for st in sb['equilibration'].get('stages', []) if st.get('last_full_frame') for k in st['last_full_frame']}
                 - {'frame', 'time_ns', 'n_water_ring_wide_3p5A', 'd_water_ring_wide_min_A'})
    for key, new in (('new_keys_of_table_row', [c for c in sb['table_row'] if c not in [new_name(x) for x in sa['table_row']]]),
                     ('new_keys_of_summary', [c for c in sb if c not in sa])):
        assert res[key] in (None, new), (r, key, new)
        res[key] = new
    if lfk:
        assert res['new_keys_of_last_full_frame'] in (None, lfk), (r, lfk)
        res['new_keys_of_last_full_frame'] = lfk
    rec['warnings_old_new'] = [len(sa['warnings']), len(sb['warnings'])]
    rec['all_same'] = bool(rec['frame_counts_equal'] and worst == 0.0
                           and all(v is True for k, v in rec.items() if 'identical' in k or 'kept' in k or 'neither' in k) and all(same.values()))
    all_ok &= rec['all_same']
    res['runs'][name] = rec
ca, cb = (json.loads((d / 'runs' / 'analysis_campaign' / 'campaign.json').read_text()) for d in (OLD, NEW))
res['campaign'] = {'rules_per_model identical': ca['rules_per_model'] == cb['rules_per_model'],
                   'folders left out identical': [Path(x).name for x in ca['skipped']] == [Path(x).name for x in cb['skipped']],
                   'runs in the table identical': [x['table_row']['run'] for x in ca['runs']] == [x['table_row']['run'] for x in cb['runs']]}
all_ok &= all(res['campaign'].values())
res.update(n_runs=len(runs), n_files=sum(v['files_compared'] for v in res['runs'].values()),
           n_frames=sum(v['frames_compared'] for v in res['runs'].values()),
           largest_coordinate_difference_A=max(v['largest_coordinate_difference_A'] for v in res['runs'].values()), all_same=bool(all_ok),
           runs_with_another_number_of_warnings={k: v['warnings_old_new'] for k, v in res['runs'].items() if v['warnings_old_new'][0] != v['warnings_old_new'][1]},
           spec_old=sa['spec'], spec_new=sb['spec'], analysis_script_sha256_old=sa['script']['sha256'], analysis_script_sha256_new=sb['script']['sha256'])
OUT.write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps({k: v for k, v in res.items() if k != 'runs'}, indent=1))
for k, v in res['runs'].items():
    if not v['all_same']:
        print('DIFFERENT', k, json.dumps(v))
print('version 1.2 against version 1.3 on', len(runs), 'runs:', 'ALL THE SAME' if all_ok else 'DIFFERENCES')
