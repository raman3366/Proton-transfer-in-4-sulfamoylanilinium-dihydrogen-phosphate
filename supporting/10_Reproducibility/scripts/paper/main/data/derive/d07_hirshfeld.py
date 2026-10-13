"""d07: Hirshfeld-surface contact percentages re-aggregated from the per-vertex data.

Sums vertex areas by (inside element, outside element) from <ion>_vertices.csv.gz with our own code, so the
percentages in hirshfeld_summary.json are checked against the surface itself; also the smallest d_norm and the
closest contacts.  The surface was built from the X-ray structure with neutron-normalised X-H lengths
(C-H 1.083, N-H 1.009, O-H 0.983 A); the Multiwfn cross-check is read from hirshfeld_compare.json.
"""
import csv
import re
import gzip
from common import *

HD = RES / 'hirshfeld_xray'
out = {}
for ion in ['cation', 'anion']:
    area = {}
    tot = 0.0
    dmin = (9, None)
    with gzip.open(HD / f'{ion}_vertices.csv.gz', 'rt') as f:
        for r in csv.DictReader(f):
            a = float(r['vertex_area_A2'])
            ei = re.sub(r'[^A-Z]', '', r['inside_atom'].split('#')[0])[0]
            eo = re.sub(r'[^A-Z]', '', r['outside_atom'].split('#')[0])[0]
            key = '...'.join(sorted([ei, eo]))
            area[key] = area.get(key, 0) + a
            tot += a
            dn = float(r['d_norm'])
            if dn < dmin[0]:
                dmin = (dn, f"{r['inside_atom']}...{r['outside_atom']}")
    pct = {k: 100 * v / tot for k, v in sorted(area.items(), key=lambda x: -x[1])}
    summ = json.loads((HD / 'hirshfeld_summary.json').read_text())['targets'][ion]
    ref = {c['contact']: c['pct'] for c in summ['contacts_reciprocal']}
    # map the reciprocal-combined names 'H...O/O...H' to our sorted keys
    refk = {'...'.join(sorted(k.split('/')[0].split('...'))): v for k, v in ref.items()}
    out[ion] = dict(area_A2=tot, pct=pct, max_abs_pct_diff_vs_summary=max(abs(pct.get(k, 0) - v) for k, v in refk.items()),
                    d_norm_min=dmin[0], d_norm_min_pair=dmin[1], summary_area_A2=summ['area_A2'],
                    summary_volume_A3=summ['volume_A3'], globularity=summ['globularity'],
                    n_neighbours=summ['n_neighbours_in_contact'])
cmp_ = json.loads((HD / 'hirshfeld_compare.json').read_text())
out['multiwfn_compare'] = cmp_
cl = json.loads((HD / 'cluster.json').read_text())
out['xh_lengths_A'] = cl['xh_lengths_A']
out['cluster_cutoff_A'] = cl['cutoff_A']
dump('d07_hirshfeld.json', out)
