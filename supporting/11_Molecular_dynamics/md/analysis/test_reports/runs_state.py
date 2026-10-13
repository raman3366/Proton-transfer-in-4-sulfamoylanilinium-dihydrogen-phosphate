# throwaway (version 1.3): what /root/cqc-sfa-md/runs/ holds on this machine.  Nothing but the listed fields of each
# manifest.json is read (run id, test flags, system folder, status); no trajectory, no log, no number of a run.
# A run of the campaign is a folder named <start><number> (X, N or R with digits: X0, X1, R3, ...), or a folder whose
# manifest has such a run id and is not marked as a test.  Prints one line per folder and the count.
import json, os, re, sys
from datetime import datetime, timezone, timedelta
ROOT = sys.argv[1] if len(sys.argv) > 1 else '/root/cqc-sfa-md/runs'
ist = datetime.now(timezone(timedelta(hours=5, minutes=30)))
print(f'state of {ROOT} on {ist.day} {ist.strftime("%b %Y")}, {ist.strftime("%I:%M %p").lstrip("0")} IST')
if not os.path.isdir(ROOT):
    print('   the folder does not exist')
    print('runs of the campaign found: 0 (folders with a manifest or a trajectory: 0)')
    sys.exit(0)
found, n_campaign = [], 0
for base, dirs, files in os.walk(ROOT):
    dirs.sort()
    if 'manifest.json' in files or any(re.fullmatch(r'solute_(seg|eq)\w*\.dcd', f) for f in files):
        dirs[:] = []
        name = os.path.basename(base)
        man = {}
        if 'manifest.json' in files:
            try:
                man = json.load(open(os.path.join(base, 'manifest.json')))
            except Exception as e:
                man = {'_error': str(e)}
        rid = man.get('run_id')
        by_name = bool(re.fullmatch(r'[XNR]\d+', name))
        by_manifest = isinstance(rid, str) and bool(re.fullmatch(r'[XNR]\d+', rid)) and man.get('test_only') is False
        campaign = by_name or by_manifest
        n_campaign += campaign
        found.append(base)
        print(f'   {os.path.relpath(base, ROOT)}: run id {rid} | test_only {man.get("test_only")} {man.get("test_only_reason")} | equil_scale {man.get("equil_scale")}'
              f' | system {"/".join(str(man.get("system_dir")).split("/")[-2:])} | status {man.get("status")} | host {((man.get("invocations") or [{}])[-1].get("host") or {}).get("hostname") if isinstance((man.get("invocations") or [{}])[-1].get("host"), dict) else None}'
              f' | production files {sum(1 for f in files if re.fullmatch(r"solute_seg[0-9]+[.]dcd", f))} | marker TEST_ONLY {"TEST_ONLY" in files}'
              f' | a run of the campaign: {"YES" if campaign else "no"}')
print(f'runs of the campaign found: {n_campaign} (folders with a manifest or a trajectory: {len(found)})')
