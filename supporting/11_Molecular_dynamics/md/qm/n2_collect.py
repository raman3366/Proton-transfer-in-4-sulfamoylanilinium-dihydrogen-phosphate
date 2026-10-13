# SANITISED COPY (SAPH hand-over package): a network address was masked.
# The original is in the git repository; the key itself never was.
"""Gate G-R2: bring the result files of the single points into md/qm/n2/qm/ WITHOUT opening a test result.

    /root/cqc-sfa-tools/md/bin/python md/qm/n2_collect.py fit  [--from-machine]     the fitting points only
    N2_FINAL_TEST=1 ... md/qm/n2_collect.py test [--from-machine]                   the test points; only for the
                                                                                    final test, immediately before it
    ... md/qm/n2_collect.py aux                                                     from the machine: the log and the
        status file of the batch (they hold no energy) -> md/qm/n2/batch/, and the files of the potential job of the
        R-like placement -> md/qm/n2/esp_Rx/ (the adapter file in/CA2_large_mk.log, results/, status.json, run_qm.log)

Source: the copy the watchdog of the Pi keeps of the machine's repository
(/root/cqc-sfa-md/from_gpu/repo_md/qm/n2/qm/, every 15 minutes), or with --from-machine the machine itself
(rsync over ssh of exactly the listed files).  Which files are copied is decided from md/qm/n2/split.json alone:
with `fit` no file of a point flagged "test" is touched, so no test result lies in the repository before the final
test (design check of 6 Oct 2026, finding on keeping test energies unread).  The files are copied as bytes; this
program reads no energy.  A file already in md/qm/n2/qm/ with other content is reported and NOT replaced.

Writes md/qm/n2/collected_<set>.json: the identifiers copied with the checksum of every file, and what was missing.
"""
import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import n2_common as nc

WATCHDOG = pathlib.Path('/root/cqc-sfa-md/from_gpu/repo_md/qm/n2/qm')
MACHINE = 'ubuntu@<address of the rented GPU machine, which was deleted on 8 Oct 2026>'
SSH = 'ssh -i /root/.ssh/aitf_pi -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15'
REMOTE = '/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/qm/n2/qm/'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('set', choices=['fit', 'test', 'aux'])
    ap.add_argument('--from-machine', action='store_true')
    args = ap.parse_args()
    if args.set == 'aux':
        got = {}
        root = REMOTE[:-len('md/qm/n2/qm/')]
        for src, dst in (('/root/cqc-sfa-md/qm/n2/batch.log', nc.N2 / 'batch' / 'batch.log'),
                         ('/root/cqc-sfa-md/qm/n2/status.json', nc.N2 / 'batch' / 'status.json'),
                         ('/root/cqc-sfa-md/qm/n2/esp_Rx/status.json', nc.N2 / 'esp_Rx' / 'status.json'),
                         ('/root/cqc-sfa-md/qm/n2/esp_Rx/run_qm.log', nc.N2 / 'esp_Rx' / 'run_qm.log'),
                         (root + 'md/qm/n2/esp_Rx/in/CA2_large_mk.log', nc.N2 / 'esp_Rx' / 'in' / 'CA2_large_mk.log'),
                         (root + 'md/qm/n2/esp_Rx/results/', nc.N2 / 'esp_Rx' / 'results'),
                         (root + 'md/qm/n2/check_X.json', nc.N2 / 'check_X.json'),
                         (root + 'md/qm/n2/check_d3.json', nc.N2 / 'check_d3_machine.json')):
            pathlib.Path(dst).parent.mkdir(parents=True, exist_ok=True)
            r = subprocess.run(['rsync', '-a', '--ignore-missing-args', '--rsync-path=sudo rsync', '-e', SSH,
                                f'{MACHINE}:{src}', str(dst) + ('/' if src.endswith('/') else '')], capture_output=True, text=True)
            got[src] = 'copied' if pathlib.Path(dst).exists() else f'not there ({r.returncode})'
        for k, v in got.items():
            print(v, k)
        return 0
    if args.set == 'test' and os.environ.get('N2_FINAL_TEST') != '1':
        raise SystemExit('the test results are copied only for the final test, with N2_FINAL_TEST=1 (md/N2.md, section 3.3)')
    split = nc.split_record()
    ids = [p for p in split['order'] if split['set'][p] == args.set]
    dst = nc.N2 / 'qm'
    dst.mkdir(exist_ok=True)
    if args.from_machine:
        src = pathlib.Path(tempfile.mkdtemp(prefix='n2_collect_'))
        lst = src / 'files.txt'
        lst.write_text(''.join(f'{p}.json\n' for p in ids))
        r = subprocess.run(['rsync', '-a', '--ignore-missing-args', f'--files-from={lst}', '--rsync-path=sudo rsync',
                            '-e', SSH, f'{MACHINE}:{REMOTE}', str(src) + '/'], capture_output=True, text=True)
        if r.returncode not in (0, 23, 24):
            raise SystemExit(f'rsync from the machine failed ({r.returncode}): {r.stderr[-400:]}')
        source = f'{MACHINE}:{REMOTE} (rsync of the listed files)'
    else:
        src, source = WATCHDOG, str(WATCHDOG)
    copied, same, missing, differ = {}, {}, [], []
    for p in ids:
        f, t = src / f'{p}.json', dst / f'{p}.json'
        if not f.exists():
            missing.append(p)
            continue
        h = nc.sha(f)
        if t.exists():
            if nc.sha(t) == h:
                same[p] = h
            else:
                differ.append(p)
            continue
        tmp = t.with_name(t.name + f'.tmp{os.getpid()}')
        shutil.copyfile(f, tmp)
        os.replace(tmp, t)
        copied[p] = h
    if args.from_machine:
        shutil.rmtree(src, ignore_errors=True)
    present = dict(same, **copied)
    rec = dict(what=f'result files of the points flagged "{args.set}" copied into md/qm/n2/qm/ by md/qm/n2_collect.py; no '
                    'energy was read', written=nc.now(), source=source, set=args.set, n_in_the_split=len(ids),
               n_present=len(present), copied_now=sorted(copied), not_there_yet=missing,
               in_the_repository_with_other_content_NOT_replaced=differ, sha256=present,
               split_sha256=nc.sha(nc.N2 / 'split.json'))
    nc.write_json(nc.N2 / f'collected_{args.set}.json', rec)
    other = 'test' if args.set == 'fit' else 'fit'
    stray = [p for p in split['order'] if split['set'][p] == other and (dst / f'{p}.json').exists()]
    print(f"{args.set}: {len(present)} of {len(ids)} present ({len(copied)} copied now), {len(missing)} not there yet, "
          f"{len(differ)} differ and were left; files of the other set in md/qm/n2/qm/: {len(stray)}")
    return 1 if differ else 0


if __name__ == '__main__':
    sys.exit(main())
