# throwaway (version 1.3): the broken copies of mutate.py run in stages, for a machine that is too loaded to run
# every copy on all 23 synthetic runs (6 Oct 2026, 11 AM IST: one copy on all runs took 4 to 10 minutes of wall time,
# about 100 s of processor time).  The copies themselves are those of mutate.py (its dictionary MUT is used as it is).
#
# The rule, fixed here before any copy was run in this way:
#   stage A  the broken copy analyses only the runs named for it by `chosen()` below; every other run keeps the output
#            of the unbroken script (the test of run_all.sh, which returns every planted value).  The comparison with
#            the planted values (make_test_trajectory.py --compare, unchanged) is then made on a folder that holds the
#            planted files of these runs alone, with what the comparison always reads (the summaries of X1, X2 and N1
#            and the campaign table).  A value that is no longer returned can therefore only come from the copy.
#   stage B  only if stage A returns every planted value: the copy on all runs and the comparison on all planted
#            values, exactly as mutate.py does it.
#   A copy is noticed if a planted value is no longer returned or if the analysis stops with an error.  A copy noticed
#   in stage A would also be noticed on all runs: the comparison on all runs contains the comparisons of stage A.  The
#   number of values given for such a copy is that of the runs analysed, a lower bound of the number on all runs.
# Which runs (`chosen`):
#   - the copies of section 9D and of the comparison with the driver's record that are new or whose code changed in
#     version 1.3 (M58 to M65, M76 to M84): the three runs of version 1.3, default/R1, place/R1 and ext/R1;
#   - a copy that version 1.2 had (the others): the two smallest runs (production frames; then by name) among the
#     runs in which that copy showed a difference in the test of version 1.2 (v12/test_reports/mutation_results.json);
#     if it showed only in the campaign table, the campaign table is written again by the copy from the summaries that
#     are there (--reuse); if it stopped the analysis in version 1.2, the runs ext/R1 and glu/X0.
# Copies that mutate.py has already run on all runs (a record without "mode") are kept as they are.
# usage: mutate_staged.py SCRATCH ANALYSIS_DIR PYTHON [M17,M18,...]
import json, os, shutil, subprocess, sys
from pathlib import Path
mdir = Path(__file__).resolve().parent
head = (mdir / 'mutate.py').read_text()
head = head[:head.index("\nfor name, pairs in MUT.items():") + 1]          # all of mutate.py before its own loop: S, A, PY, src, MUT
ONLY_ARG = sys.argv[4] if len(sys.argv) > 4 else None
sys.argv = sys.argv[:4]                                       # mutate.py reads S, A, PY from the command line; ONLY is handled here
exec(head)                                                    # defines S, A, PY, src, MUT
ONLY = ONLY_ARG.split(',') if ONLY_ARG else None
for name, pairs in MUT.items():
    for old, new in pairs:
        assert src.count(old) == 1, (name, src.count(old), old)
REF, T, SLIM = S / 'test', S / 'test_mut', S / 'test_slim'
env = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
out_file = mdir / 'mutation_results.json'
res = json.loads(out_file.read_text()) if out_file.is_file() else {}
for k, v in res.items():
    v.setdefault('mode', 'all runs')
V12 = json.loads((S / 'v12/test_reports/mutation_results.json').read_text())
V12 = {k.split()[0]: v for k, v in V12.items()}
planted = {}                                                  # label of a run as the reports name it -> tag, folder, frames
for p in sorted((REF / 'planted').glob('runs*__*__*.json')):
    P = json.loads(p.read_text())
    planted[P['run'] + (' (overlap)' if p.stem.startswith('runs_overlap__') else '')] = {'tag': p.stem, 'folder': P['folder'], 'frames': P['n_frames_solute']}
V13_RUNS = ['default/R1', 'place/R1', 'ext/R1']
V13_COPIES = {f'M{n}' for n in list(range(58, 66)) + list(range(76, 85))}


def chosen(number):
    """The runs of stage A, and whether the campaign table is written again; see the head of this file."""
    if number in V13_COPIES:
        return list(V13_RUNS), False
    old = V12[number]
    if old['n_failed'] is None:
        return ['ext/R1', 'glu/X0'], False
    runs = sorted((r for r in old['runs_with_a_failure'] if r != 'campaign'), key=lambda r: (planted[r]['frames'], r))
    if not runs:
        return [], True
    return runs[:2], False


def restore(labels, campaign):
    """Put the output of the unbroken script back (the folder analysis/ of a run, the campaign table)."""
    for r in labels:
        d = T / planted[r]['folder'] / 'analysis'
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(REF / planted[r]['folder'] / 'analysis', d)
    if campaign:
        shutil.rmtree(T / 'runs' / 'analysis_campaign', ignore_errors=True)
        shutil.copytree(REF / 'runs' / 'analysis_campaign', T / 'runs' / 'analysis_campaign')


def slim(labels):
    """A folder with the planted files of these runs only, and what the comparison always reads."""
    shutil.rmtree(SLIM, ignore_errors=True)
    (SLIM / 'planted').mkdir(parents=True)
    os.symlink(T / 'planted' / 'planted.json', SLIM / 'planted' / 'planted.json')
    for r in labels:
        tag, folder = planted[r]['tag'], planted[r]['folder']
        for f in (T / 'planted').iterdir():
            if f.name == tag + '.json' or f.name.startswith(tag + '__'):
                os.symlink(f, SLIM / 'planted' / f.name)
        (SLIM / folder).parent.mkdir(parents=True, exist_ok=True)
        os.symlink(T / folder, SLIM / folder)
    for r in ('X1', 'X2', 'N1'):
        if not (SLIM / 'runs' / 'test' / r).exists():
            (SLIM / 'runs' / 'test').mkdir(parents=True, exist_ok=True)
            os.symlink(T / 'runs' / 'test' / r, SLIM / 'runs' / 'test' / r)
    (SLIM / 'runs').mkdir(parents=True, exist_ok=True)
    os.symlink(T / 'runs' / 'analysis_campaign', SLIM / 'runs' / 'analysis_campaign')


RUN_SOME = '''import sys
sys.path.insert(0, %r)
import %s as m
opts = ['--no-figures', '--no-hash-trajectories']
codes = [m.main(a + opts) for a in %r]
print('exit codes', codes)
sys.exit(max(codes))
'''
RUN_ALL = '''import sys, glob
sys.path.insert(0, %r)
import analyse_md_mutant as m
T = %r
opts = ['--no-figures', '--no-hash-trajectories']
codes = [m.main(['--campaign', T + '/runs'] + opts), m.main([T + '/runs/test/X0'] + opts), m.main([T + '/runs_overlap/test/X9'] + opts)]
codes += [m.main([d] + opts) for d in sorted(glob.glob(T + '/runs_gv/*/*/')) + sorted(glob.glob(T + '/runs_v12/*/*/')) + sorted(glob.glob(T + '/runs_v13/*/*/'))]
print('exit codes', codes)
sys.exit(max(codes))
'''


def compare(root):
    for f in (root / 'test_report.json', root / 'test_report.md'):
        f.unlink(missing_ok=True)
    r3 = subprocess.run([PY, str(A / 'make_test_trajectory.py'), '--compare', str(root)], capture_output=True, text=True, env=env)
    d = json.loads((root / 'test_report.json').read_text()) if (root / 'test_report.json').is_file() else None
    return r3, d


def stage_a(module, labels, campaign):
    args = [[str(T / planted[r]['folder'])] for r in labels] + ([['--campaign', str(T / 'runs'), '--reuse']] if campaign else [])
    (mdir / 'run_some.py').write_text(RUN_SOME % (str(mdir if module == 'analyse_md_mutant' else A), module, args))
    shutil.rmtree(mdir / '__pycache__', ignore_errors=True)
    r1 = subprocess.run([PY, str(mdir / 'run_some.py')], capture_output=True, text=True, env=env)
    slim(labels)
    r3, d = compare(SLIM)
    return r1, r3, d


def record(pairs, r1, r3, d, **extra):
    rec = dict(extra, replacements=len(pairs), analysis_exit=r1.returncode, compare_exit=r3.returncode)
    if d is not None and r1.returncode == 0:
        bad = [v for v in d['values'] if not v['passed']]
        rec.update(n_values=d['n_values'], n_failed=len(bad), runs_with_a_failure=sorted({v['run'] for v in bad}),
                   examples=[f"{v['run']}: {v['quantity']}: planted {v['planted']}, returned {v['returned']}" for v in bad[:3]])
    else:
        rec.update(n_failed=None, note=(r1.stderr.strip().splitlines() or r3.stderr.strip().splitlines() or ['?'])[-1][:300])
    return rec


# ---- the output of the unbroken script in every run, and the check that the staged comparison is silent on it
restore(list(planted), True)
r3, d = compare(T)
print(f"the unbroken script, all runs: {d['n_passed']} of {d['n_values']} planted values returned (comparison on all runs, exit code {r3.returncode})", flush=True)
assert d['n_failed'] == 0
r1, r3, d = stage_a('analyse_md', V13_RUNS, True)             # the unbroken script through the staged machinery
print(f"the unbroken script in stage A on {', '.join(V13_RUNS)} and the campaign table: analysis exit code {r1.returncode}; "
      f"{d['n_passed']} of {d['n_values']} planted values returned", flush=True)
assert r1.returncode == 0 and d['n_failed'] == 0
restore(V13_RUNS, True)

for name, pairs in MUT.items():
    number = name.split()[0]
    if (ONLY and number not in ONLY) or (name in res and res[name].get('mode') == 'all runs'):
        continue
    text = src
    for old, new in pairs:
        text = text.replace(old, new)
    (mdir / 'analyse_md_mutant.py').write_text(text)
    labels, campaign = chosen(number)
    r1, r3, d = stage_a('analyse_md_mutant', labels, campaign)
    rec = record(pairs, r1, r3, d, mode='staged', stage='A', runs_with_the_copy=labels + (['campaign table (from the summaries)'] if campaign else []))
    restore(labels, campaign)
    if rec['n_failed'] == 0:                                  # not noticed on the chosen runs: all runs, as mutate.py does it
        for dd in T.glob('runs*/**/analysis*'):
            shutil.rmtree(dd, ignore_errors=True)
        (mdir / 'run_mutant.py').write_text(RUN_ALL % (str(mdir), str(T)))
        shutil.rmtree(mdir / '__pycache__', ignore_errors=True)
        r1 = subprocess.run([PY, str(mdir / 'run_mutant.py')], capture_output=True, text=True, env=env)
        r3, d = compare(T)
        rec = record(pairs, r1, r3, d, mode='staged', stage='B', runs_with_the_copy=['all runs'], stage_A_runs_without_a_difference=labels)
        restore(list(planted), True)
    res[name] = rec
    print(name, '->', json.dumps(rec)[:900], flush=True)
    out_file.write_text(json.dumps(res, indent=1) + '\n')
for f in ('analyse_md_mutant.py', 'run_some.py'):
    (mdir / f).unlink(missing_ok=True)
shutil.rmtree(mdir / '__pycache__', ignore_errors=True)
shutil.rmtree(SLIM, ignore_errors=True)
missing = [k for k in MUT if k not in res]
detected = sum(1 for r in res.values() if r['n_failed'] is None or r['n_failed'] > 0)
print(f'{detected} of {len(res)} mutations noticed' + (f'; NOT RUN: {len(missing)}' if missing else ''))
print('copies run on all runs:', sum(1 for r in res.values() if r['mode'] == 'all runs'), '| noticed in stage A:', sum(1 for r in res.values() if r.get('stage') == 'A'),
      '| taken to stage B (all runs):', sum(1 for r in res.values() if r.get('stage') == 'B'))
