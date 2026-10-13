# throwaway (version 1.3): the dates of the files under /root/cqc-sfa-md/runs/ (file system only; no file is opened).
# usage: runs_dates.py OUT.json
import json, os, sys
from datetime import datetime, timezone, timedelta
IST = timezone(timedelta(hours=5, minutes=30))
ROOT = '/root/cqc-sfa-md/runs'


def ist(t):
    d = datetime.fromtimestamp(t, IST)
    return f'{d.day} {d.strftime("%b %Y")}, {d.strftime("%I:%M %p").lstrip("0")} IST'


out = {'folder': ROOT, 'exists': os.path.isdir(ROOT), 'looked_at': ist(datetime.now(IST).timestamp()), 'folders': {}}
if out['exists']:
    out['folder_created_or_last_changed'] = ist(os.stat(ROOT).st_mtime)
    for name in sorted(os.listdir(ROOT)):
        p = os.path.join(ROOT, name)
        times = [os.stat(os.path.join(b, f)).st_mtime for b, _, fs in os.walk(p) for f in fs]
        out['folders'][name] = {'n_files': len(times), 'oldest_file': ist(min(times)) if times else None, 'newest_file': ist(max(times)) if times else None,
                                'oldest_file_unix': min(times) if times else None}
# when the watchdog fetched files from the rented machine (its own log; one line per run of the watchdog)
import re
out['watchdog_copies_with_files'] = []
try:
    for line in open('/root/cqc-sfa-md/status/watchdog.log'):
        m = re.match(r'\[\S+ \S+ UTC \| (6 Oct 2026, \d+:\d\d [AP]M IST)\].*pull now: (\d+) file\(s\) fetched', line)
        if m and int(m.group(2)) > 0:
            out['watchdog_copies_with_files'].append({'time': m.group(1), 'files': int(m.group(2))})
except OSError as e:
    out['watchdog_log'] = str(e)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print(json.dumps(out, indent=1))
