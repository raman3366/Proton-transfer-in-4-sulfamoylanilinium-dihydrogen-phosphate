#!/usr/bin/env python
"""Gate G-R2, the examiner's integrity checks before the final test (md/N2.md, section 9).

    nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_integrity.py [--machine-sha SHA]

NO result file of a test point is parsed.  Of the test points this program uses only: the split, the geometry
files, the lines of the batch log (identifier, set, exit code, converged, seconds), file names, sizes and times
(stat), and SHA-256 checksums of the result files in the watchdog's copy (to look for copies elsewhere; a checksum
shows nothing of the content).  Result files of FITTING points and of the 36 old points are parsed.

Writes md/qm/n2/final_test/integrity.json.
"""
import argparse
import collections
import hashlib
import json
import os
import pathlib
import re
import time

import numpy as np

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
QM = REPO / 'md' / 'qm'
N2 = QM / 'n2'
GR = QM / 'gr'
OUT = N2 / 'final_test'
WATCHDOG = pathlib.Path('/root/cqc-sfa-md/from_gpu/repo_md/qm/n2/qm')
SCRATCH = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                       '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad')
SEED = 20261006


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def ist(t):
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(t + 19800))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--machine-sha', default=None, help='sha256 of split.json on the rented machine (read over ssh)')
    args = ap.parse_args()
    doc = dict(what='gate G-R2: integrity checks of the examiner before the final test; no result file of a test point was parsed',
               written=ist(time.time()))
    split = json.loads((N2 / 'split.json').read_text())
    order, flag = split['order'], split['set']
    test = [p for p in order if flag[p] == 'test']

    # ---- 1. the split
    rng = np.random.default_rng(SEED)
    mine = {}
    for b in range(0, len(order), 4):
        k = int(rng.integers(0, 4))
        for j, pid in enumerate(order[b:b + 4]):
            mine[pid] = 'test' if j == k else 'fit'
    log = {}
    first = None
    for l in (N2 / 'batch' / 'batch.log').read_text().splitlines():
        m = re.match(r'(\S+ \S+ \S+ IST)\s+(\w+) \((fit|test)\): exit code (\d+), (\w+), (\d+) s, attempt (\d+)', l)
        if m:
            log[m.group(2)] = dict(time=m.group(1), set=m.group(3), exit=int(m.group(4)), state=m.group(5),
                                   seconds=int(m.group(6)), attempt=int(m.group(7)))
            first = first or m.group(1)
    doc['split'] = dict(
        sha256=sha(N2 / 'split.json'), recorded=(N2 / 'split.json.sha256').read_text().split()[0],
        equal=sha(N2 / 'split.json') == (N2 / 'split.json.sha256').read_text().split()[0],
        on_the_rented_machine=args.machine_sha, equal_on_the_machine=(args.machine_sha == sha(N2 / 'split.json')) if args.machine_sha else None,
        files_written={f: ist(os.stat(N2 / f).st_mtime) for f in ('split.json', 'split.json.sha256', 'points.json')},
        first_single_point_of_a_new_geometry_ended=first,
        drawn_again_here_from_seed_and_rule_equals_the_file=mine == flag,
        counts=dict(collections.Counter(flag.values())),
        geometry_files_with_another_checksum=[p for p in order if sha(N2 / 'geom' / f'{p}.xyz') != split['xyz_sha256'][p]],
        set_labels_of_the_batch_log_equal_the_split=all(flag[k] == v['set'] for k, v in log.items()),
        note='the batch wrote the set of every point into its log line as the point ended (from 2:42:00 PM IST); these '
             'labels are those of split.json, so the split is the one the batch ran with')

    # ---- 2. the batch: which points, convergence
    done = [p for p in order if p in log]
    doc['batch'] = dict(
        finished=len(done), of=len(order), fit=sum(log[p]['set'] == 'fit' for p in done), test=sum(log[p]['set'] == 'test' for p in done),
        not_started=[(p, flag[p]) for p in order if p not in log], finished_points_are_a_prefix_of_the_list=done == order[:len(done)],
        all_exit_code_0_and_converged=all(v['exit'] == 0 and v['state'] == 'converged' for v in log.values()),
        all_at_the_first_attempt=all(v['attempt'] == 1 for v in log.values()),
        test_points_converged_by_the_log=sum(v['set'] == 'test' and v['state'] == 'converged' for v in log.values()),
        last_line=(N2 / 'batch' / 'batch.log').read_text().splitlines()[-1])

    # ---- 3. settings of the fitting result files against the gate's
    ref = json.loads((GR / 'qm' / 'X.json').read_text())
    gate_settings = [json.loads(f.read_text()).get('settings') for f in sorted((GR / 'qm').glob('*.json'))]
    cfg = json.loads((QM / 'config_prod.json').read_text())
    here = sorted(f.stem for f in (N2 / 'qm').glob('*.json'))
    bad, cyc, second, sizes, hosts, vers, func, guess = [], [], [], set(), set(), set(), set(), set()
    for pid in here:
        r = json.loads((N2 / 'qm' / f'{pid}.json').read_text())
        if r.get('settings') != ref['settings'] or r.get('test_level') or not r.get('converged') or \
                r.get('xyz_sha256') != split['xyz_sha256'][pid] or r.get('set') != 'fit' or r.get('functional') != ref['functional'] or \
                (r.get('natm'), r.get('nelectron'), r.get('nao'), r.get('naux')) != (ref['natm'], ref['nelectron'], ref['nao'], ref['naux']):
            bad.append(pid)
        cyc.append(r.get('scf_cycles'))
        second.append(bool(r.get('second_order_solver')))
        sizes.add(r.get('grid_points'))
        hosts.add(r.get('host'))
        vers.add(json.dumps(r.get('versions'), sort_keys=True))
        func.add(r.get('functional'))
        guess.add(json.dumps(r.get('guess'), sort_keys=True))
    doc['settings'] = dict(
        gate=ref['settings'], gate_functional=ref['functional'], all_38_files_of_gate_G_R_have_these_settings=all(s == ref['settings'] for s in gate_settings),
        config_prod_sha256=sha(QM / 'config_prod.json'), config_method=cfg['method'],
        fitting_files_parsed=len(here), fitting_files_that_differ_or_are_not_converged=bad,
        scf_cycles=[min(cyc), max(cyc)], second_order_solver_used=int(sum(second)), grid_points=sorted(sizes), hosts=sorted(hosts),
        versions=[json.loads(v) for v in vers], guess=[json.loads(g) for g in guess],
        check_X=json.loads((N2 / 'check_X.json').read_text()).get('difference_hartree'),
        test_points='NOT parsed: the settings fields of the 47 test result files were not read (the files stay unopened); the same '
                    'two runs of md/qm/n2_batch.py computed fitting and test points interleaved, with one configuration file per run '
                    '(its checksum is in every fitting file), and the log has every test point converged at the first attempt')

    # ---- 4. the repository holds no test result; no exclusion after energies
    cands = {}
    for d in ('fit', 'fit_ff2'):
        c = json.loads((N2 / d / 'candidates.json').read_text())
        ids = c['points']['ids'] if isinstance(c.get('points'), dict) and 'ids' in c['points'] else None
        want = list(split['old_rigid_points']) + [p for p in order if flag[p] == 'fit' and p in log]
        cands[d] = dict(rows=len(ids) if ids else None, rows_are_all_finished_fitting_points_and_the_36_old=(sorted(ids) == sorted(want)) if ids else None,
                        fitting_points_left_out=c.get('fitting_points_left_out'), sensitivity=c.get('sensitivity'),
                        split_sha256_in_the_fit=c.get('split_sha256'), written=c.get('written'))
    doc['repository'] = dict(
        result_files=len(here), flagged_test=[h for h in here if flag.get(h) == 'test'],
        finished_fitting_points_missing=[p for p in done if flag[p] == 'fit' and p not in here],
        absent={f: not (N2 / f).exists() for f in ('test_set_access.log', 'gate', 'collected_test.json', 'gate_preview', 'gate_preview_ff2')},
        fits_of_the_rule=cands)

    # ---- 5. traces of a read of a test result file
    real = {t: sha(WATCHDOG / f'{t}.json') for t in test if (WATCHDOG / f'{t}.json').exists()}
    realset = set(real.values())
    copies, named = [], collections.Counter()
    for root in (SCRATCH, REPO / 'md'):
        for d, _, fs in os.walk(root):
            if str(OUT) in d:
                continue
            for f in fs:
                p = os.path.join(d, f)
                if f.endswith('.json') and f[:-5] in flag and flag[f[:-5]] == 'test':
                    named[d.replace(str(SCRATCH), '<scratch>')] += 1
                try:
                    if f.endswith('.json') and os.path.getsize(p) < 20000 and sha(p) in realset:
                        copies.append(p)
                except OSError:
                    pass
    st = {t: os.stat(WATCHDOG / f'{t}.json') for t in real}
    doc['test_result_files'] = dict(
        where='rented machine (md/qm/n2/qm) and the watchdog copy on the Pi (' + str(WATCHDOG) + '); none in the repository',
        in_the_watchdog_copy=len(real),
        files_elsewhere_with_the_bytes_of_a_test_result=copies,
        folders_with_files_named_like_test_points=dict(named),
        folders_note='these are the synthetic trees of md/qm/n2_fit_selftest.py (written from a known parameter set; key synthetic) and one '
                     'rehearsal file of a test level in the design scratch; none has the bytes of a real result',
        access_times_now=dict(collections.Counter(ist(s.st_atime) if s.st_atime > s.st_ctime + 1.0 else 'at arrival' for s in st.values())),
        access_times_of_the_fitting_files_in_the_same_folder=dict(collections.Counter(
            ist(os.stat(f).st_atime) if os.stat(f).st_atime > os.stat(f).st_ctime + 1.0 else 'at arrival'
            for f in WATCHDOG.glob('*.json') if flag.get(f.stem) == 'fit')),
        first_observation=json.loads((OUT / 'first_observation.json').read_text()),
        access_time_note='file system mounted relatime: the first read after a copy arrived moves its access time past its change '
                         'time. First observation (stat only, before the examiner hashed the files): none of the 47 test files and none '
                         'of the 141 fitting files had been read since they arrived. The single later time of the test files is the '
                         "examiner's own checksum pass. This cannot show a read on the rented machine, where the watchdog's pull reads every file.",
        mount=os.popen(f'findmnt -no OPTIONS -T {WATCHDOG}').read().strip())

    # ---- 5b. energies of new geometries at another level (the designer's rehearsal of the runner); keys only
    reh = []
    for f in sorted((SCRATCH / 'n2' / 'design').glob('reh*_out/qm/*.json')):
        r = json.loads(f.read_text())
        pid = r.get('id')
        reh.append(dict(file=str(f).replace(str(SCRATCH), '<scratch>'), id=pid, set_in_the_file=r.get('set'), test_level=r.get('test_level'),
                        basis=(r.get('settings') or {}).get('basis'), grid=(r.get('settings') or {}).get('grid'), host=r.get('host'),
                        started=r.get('started'), same_geometry_as_the_point_of_the_final_list=r.get('xyz_sha256') == split['xyz_sha256'].get(pid),
                        set_in_the_final_list=flag.get(pid)))
    doc['rehearsal_energies_at_a_test_level'] = dict(
        files=reh, note='energies of the designer\'s rehearsal of the runner (md/N2.md, section 8.3), marked test_level and refused by every '
                        'program of the gate; one of them is of the geometry of a test point, at another level of theory; only the keys '
                        'listed here were read')

    # ---- 6. the code that judges
    g = (QM / 'n2_gate.py').read_text()
    tol = re.search(r'TOL = dict\((.*?)\)', g, re.S).group(1).replace('\n', ' ')
    doc['code'] = dict(
        sha256={f: sha(QM / f) for f in ('n2_common.py', 'n2_fit.py', 'n2_gate.py', 'n2_points.py', 'n2_batch.py', 'n2_collect.py', 'gr_points.py',
                                         'gr_queue.py', 'gr_mm.py', 'gr_fit.py')},
        tolerances_in_n2_gate=' '.join(tol.split()),
        fit_program_reads_energies_through="n2_common.fit_rows('fit') -> load_points('fit'): a point flagged test is skipped before its file is opened",
        design_copies=dict(n2_gate=ist(os.stat(SCRATCH / 'n2/build/before/n2_gate.py').st_mtime),
                           n2_common=ist(os.stat(SCRATCH / 'n2/build/before/n2_common.py').st_mtime)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'integrity.json').write_text(json.dumps(doc, indent=1) + '\n')
    for k in ('split', 'batch', 'repository', 'test_result_files'):
        print(k, json.dumps({a: b for a, b in doc[k].items() if not isinstance(b, str) or len(b) < 120}, indent=None)[:1800])
    s = doc['settings']
    print('settings:', s['fitting_files_parsed'], 'fitting files; differing:', s['fitting_files_that_differ_or_are_not_converged'], '; cycles', s['scf_cycles'],
          '; second-order', s['second_order_solver_used'], '; all 38 gate files same settings:', s['all_38_files_of_gate_G_R_have_these_settings'],
          '; hosts', s['hosts'], '; versions', s['versions'], '; guess', s['guess'])
    print('tolerances:', doc['code']['tolerances_in_n2_gate'])


if __name__ == '__main__':
    main()
