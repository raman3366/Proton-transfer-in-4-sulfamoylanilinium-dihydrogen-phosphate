"""
Compact progress table across every case that has been started.
Usage: python dft/progress.py [case ...]      (default: every case under molecules/*/results/)
"""
import sys, json, pathlib, datetime, re

ROOT = pathlib.Path(__file__).resolve().parent.parent
STAGES = ['optimize', 'proton_check', 'frequencies', 'fmo', 'mwfn', 'plots', 'render', 'report']


def frames(p):
    if not p.exists(): return 0
    try:
        with open(p) as f:
            n = int(f.readline())
        return sum(1 for _ in open(p)) // (n + 2)
    except Exception:
        return 0


def case_rows(names=None):
    dirs = sorted(ROOT.glob('molecules/*/results/*')) if not names else \
           [d for d in sorted(ROOT.glob('molecules/*/results/*')) if d.name in names]
    rows = []
    for d in dirs:
        if not d.is_dir(): continue
        st = json.load(open(d / 'status.json')) if (d / 'status.json').exists() else {'stages': {}}
        done = [s for s in STAGES if st['stages'].get(s, {}).get('status') == 'done']
        cur = next((s for s in STAGES if st['stages'].get(s, {}).get('status') == 'running'), None)
        failed = [s for s in STAGES if st['stages'].get(s, {}).get('status') == 'failed']
        extra = ''
        if cur == 'optimize':
            log = (d / 'logs' / 'optimize.log')
            txt = log.read_text(errors='replace') if log.exists() else ''
            main_started = 'main optimisation at' in txt
            steps = re.findall(r'^ *Step +(\d+) :', txt, re.M)
            last = int(steps[-1]) if steps else 0
            if main_started:
                pre = re.search(r'pre-optimisation done: (\d+) steps, ([\d.]+) min', txt)
                extra = f'main step {frames(d / "opt_traj.xyz")}' + (f' (pre {pre.group(1)})' if pre else '')
            else:
                extra = f'pre-opt step {last}'
        elif cur == 'frequencies':
            logs = sorted((d / 'logs').glob('fd_worker*.log'))
            tot = sum(len(re.findall(r'column \d+/', l.read_text(errors="replace"))) for l in logs)
            ncol = 0
            m = re.search(r'column \d+/(\d+)', ''.join(l.read_text(errors='replace') for l in logs)) if logs else None
            if m: ncol = int(m.group(1))
            extra = f'{tot}/{ncol} FD columns' if ncol else 'starting'
        elif cur is None and len(done) == len(STAGES):
            extra = 'COMPLETE'
        frq = json.load(open(d / 'freq.json')) if (d / 'freq.json').exists() else None
        pc = json.load(open(d / 'proton_check.json')) if (d / 'proton_check.json').exists() else None
        rows.append(dict(case=d.name, done=len(done), total=len(STAGES), stage=cur or ('done' if len(done) == len(STAGES) else '-'),
                         extra=extra, failed=','.join(failed),
                         imag=(frq or {}).get('n_imaginary'),
                         proton=('transferred' if (pc or {}).get('proton_transferred') else ('retained' if pc else None)),
                         started=st.get('pipeline_start', ''), ended=st.get('pipeline_end', '')))
    return rows


if __name__ == '__main__':
    rows = case_rows(sys.argv[1:] or None)
    w = max([len(r['case']) for r in rows] + [4])
    print(f"{'case'.ljust(w)}  stages  current        detail              proton      imag  ")
    print('-' * (w + 58))
    for r in rows:
        print(f"{r['case'].ljust(w)}  {r['done']}/{r['total']}     {str(r['stage'])[:13].ljust(13)}  "
              f"{str(r['extra'])[:18].ljust(18)}  {str(r['proton'] or '-')[:10].ljust(10)}  "
              f"{str(r['imag'] if r['imag'] is not None else '-').ljust(4)}  {('FAILED: ' + r['failed']) if r['failed'] else ''}")
    # Indian Standard Time (UTC+5:30), 12-hour -- the user reads these reports in India
    ist = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=5, minutes=30)
    print(f"\n{ist:%a %d %b %Y, %I:%M %p} IST")
