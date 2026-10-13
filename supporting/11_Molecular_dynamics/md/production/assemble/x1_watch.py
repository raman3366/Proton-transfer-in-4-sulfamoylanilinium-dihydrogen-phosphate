#!/usr/bin/env python3
"""Makes x1/watch.json for the production document from: the copy of the machine's run_md.log of X1, the pulled
folder of X1 on the Pi, the saved output of the status program and the watchdog log."""
import hashlib, json, os, re, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
S = Path(sys.argv[1]); X1 = Path('/root/cqc-sfa-md/runs/bonded_2e/X1')
IST = timezone(timedelta(hours=5, minutes=30))
now = datetime.now(IST)
fmt = lambda d: d.strftime('%-d %b %Y, %-I:%M %p IST')
log = (S / 'x1/run_md.log.machine').read_text().splitlines()
T = lambda s: datetime.strptime(s.replace(' IST', ''), '%d %b %Y, %I:%M:%S %p').replace(tzinfo=IST)
segs = []
for ln in log:
    m = re.search(r'\[(.*?)\] (seg\d+) \((\d+)-(\d+) ns\): done in (\d+) s, ([\d.]+) ns/day', ln)
    if m:
        segs.append((m.group(2), T(m.group(1)), int(m.group(5)), float(m.group(6))))
seg1_start = [T(re.match(r'\[(.*?)\]', ln).group(1)) for ln in log if 'seg0001 (0-5 ns): NPT' in ln][0]
task_start = datetime(2026, 10, 6, 22, 14, 22, tzinfo=IST)          # queue status: X1 running since 10:14:22 PM IST
before = int((seg1_start - task_start).total_seconds())
eq_s = sum(int(re.search(r'done in (\d+) s', ln).group(1)) for ln in log if re.search(r'\] eq\d+ .*done in', ln))
min_s = before - eq_s
ns_day = sum(s[3] for s in segs) / len(segs)
# exact: wall seconds of the finished segments
ns_day_wall = 5.0 * len(segs) / sum(s[2] for s in segs) * 86400.0
prod_s = round(100.0 / ns_day_wall * 86400.0)
run_s = before + prod_s
eq3 = [float(re.search(r'([\d.]+) ns/day', ln).group(1)) for ln in log if re.search(r'\] eq03 .*done in', ln)][0]
ext_s = round(0.5 / eq3 * 86400.0)
# pulled files
man = json.load(open(X1 / 'manifest.json')) if (X1 / 'manifest.json').exists() else None
pulled = []
if man:
    for seg in man['production'].get('segments') or []:
        for fn, rec in seg['files'].items():
            p = X1 / fn
            if p.exists():
                h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
                pulled.append((seg['tag'], fn, h == rec['sha256']))
segs_ok = sorted({t for t, f, ok in pulled})
all_ok = bool(pulled) and all(ok for _, _, ok in pulled)
wd = Path('/root/cqc-sfa-md/status/watchdog.log').read_text().splitlines()
last_pull = [ln for ln in wd if 'pull now:' in ln][-1]
mp = re.search(r'\| (\d+ \w+ \d+, [\d:]+ [AP]M IST)\].*pull now: (\d+) file\(s\) fetched, ([\d.]+) MB; driver checksums: (\d+) ok, (\d+) bad', last_pull)
n_seg_files = {t: sum(1 for a, f, ok in pulled if a == t) for t in segs_ok}
pull_text = (f"The watchdog's run that ended at {mp.group(1)[:-4].split(', ')[1]} IST (pulls start every 15 minutes) fetched {mp.group(2)} files ({mp.group(3)} MB); its line says "
             f"\"driver checksums: {mp.group(4)} ok, {mp.group(5)} bad\". On the Pi, `/root/cqc-sfa-md/runs/bonded_2e/X1/` then held the "
             f"minimisation, the six equilibration stages and {'segment' if len(segs_ok) == 1 else 'segments'} {', '.join(segs_ok) if segs_ok else 'none yet'} "
             + (f"({sum(n_seg_files.values())} files of {'this segment' if len(segs_ok) == 1 else 'these segments'}; the sha256 of each was computed again here and equals the driver's "
                f"record in the pulled manifest: {str(all_ok).lower()})." if segs_ok else "(no finished segment yet)."))
status = (S / 'x1/status.txt').read_text()
lines = status.splitlines()
i0 = [k for k, ln in enumerate(lines) if ln.startswith('**MD runs**')][0]
exc = [ln for ln in lines[i0:] if ln.strip()]
exc = [ln for ln in exc if not ln.startswith('| runs_test/')]
exc_txt = '\n'.join(exc[:40])
m_end = re.search(r'\*\*Production ends: about (.*?)\*\*', status)
m_cost = re.search(r'Cost at that end: about (US\$\d+\.\d\d) of the cap (US\$\d+\.\d\d)', status)
m_q = re.search(r'\*\*Queue\*\* \((.*?)\): (.*)', status)
w = {'written_ist': fmt(now), 'log_copied_ist': fmt(datetime.fromtimestamp(os.path.getmtime(S / 'x1/run_md.log.machine'), IST)),
     'pull_text': pull_text, 'status_ist': fmt(datetime.fromtimestamp(os.path.getmtime(S / 'x1/status.txt'), IST)),
     'status_excerpt': exc_txt,
     'speed_text': (f"{len(segs)} finished segment(s) of X1 at " + ', '.join(f'{s[3]}' for s in segs) + f" ns per day ({sum(s[2] for s in segs)} s for {5 * len(segs)} ns, "
                    f"that is {ns_day_wall:.1f} ns per day); the two validation runs gave 930.15 and 925.03."),
     'ns_day': ns_day_wall, 'prod_s': prod_s, 'eq_s': eq_s, 'min_s': min_s, 'run_s': run_s, 'ext_s': ext_s,
     'end_nine': fmt(task_start + timedelta(seconds=9 * run_s)), 'end_zn': fmt(task_start + timedelta(seconds=10 * run_s)),
     'end_all': fmt(task_start + timedelta(seconds=11 * run_s)),
     'status_end': ('"Production ends: about ' + m_end.group(1) + '"' if m_end else 'no end time was printed'),
     'cost_text': (f"\nCost: the status program gives about {m_cost.group(1)} at its end time, against the cap of {m_cost.group(2)}."
                   if m_cost else ''),
     'seg_mb': round(sum((X1 / fn).stat().st_size for tg, fn, ok in pulled if tg == 'seg0001') / 1e6, 1),
     'queue_line': (m_q.group(0).replace('**', '') if m_q else None), 'segments_pulled': segs_ok, 'pulled_files_checksums_ok': all_ok,
     'n_segments_done_in_log': len(segs)}
json.dump(w, open(S / 'x1/watch.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in w.items() if k != 'status_excerpt'}, indent=1))
print(exc_txt)
