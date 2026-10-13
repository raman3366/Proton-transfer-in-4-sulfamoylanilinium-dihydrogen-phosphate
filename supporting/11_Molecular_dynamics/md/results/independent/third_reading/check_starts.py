#!/usr/bin/env python
"""Check of the third reader's own selections, reference and definitions on the
three starts as built (coordinates of solute.pdb = system.inpcrd), against
md/system/starts.csv.  Also the reference built here against reference.pdb of
the built systems (an input of the runs, not an output of an analysis)."""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import third_lib as T

ref = T.build_reference()
out = {'reference': {'n_ca': int(len(ref['ca'])),
                     'source': 'md/system/protein_md.pdb (CA) and md/system/ligand_X.pdb (11 heavy atoms)'}}

# reference.pdb of the built systems, if there
cmp_ = {}
for s in ('bonded_2e/X', 'bonded_2e/N', 'bonded_2e/R', 'bonded_2s/R', 'bonded_2e_n2fit/R'):
    fn = os.path.join('/root/cqc-sfa-md/systems', s, 'reference.pdb')
    if not os.path.exists(fn):
        cmp_[s] = 'no reference.pdb'
        continue
    atoms, _ = T.read_pdb(fn)
    ca = np.array([a['xyz'] for a in atoms if a['name'] == 'CA' and a['resname'] not in ('SN1', 'SAN', 'ZN1', 'ZN')])
    lig = {a['name']: np.array(a['xyz']) for a in atoms if a['resname'] in ('SN1', 'SAN')}
    d = {'n_ca': int(len(ca))}
    if len(ca) == len(ref['ca']):
        d['max_abs_diff_ca_A'] = float(np.abs(ca - ref['ca']).max())
    d['max_abs_diff_lig_heavy_A'] = float(max(np.abs(lig[n] - ref['lig'][n]).max() for n in T.LIG_HEAVY))
    cmp_[s] = d
out['reference']['against_reference_pdb_of_built_systems'] = cmp_

starts = {}
with open(os.path.join(T.SYS, 'starts.csv')) as f:
    for row in csv.DictReader(f):
        starts[row['start']] = row

rows = {}
for zset, run, start, _ in T.RUN_LIST:
    folder = os.path.join(T.RUNS, zset, run)
    sel = T.select(os.path.join(folder, 'solute.pdb'))
    # CA order: tleap index of solute.pdb against deposited number of the reference
    rmap = T.residue_map()
    order_ok = all(rmap[dep][0] == tl for dep, tl in zip(ref['ca_resseq'], sel['ca_resseq']))
    xyz = sel['xyz'][None, :, :]
    lat = None
    q = T.frame_quantities(xyz, lat, sel, ref)
    st = T.states_of(q['d1'], q['d2'])[0]
    s = starts[start]
    rows[zset + '/' + run] = {
        'n_atoms': sel['n_atoms'], 'n_ca': int(len(sel['ca'])), 'ca_order_matches_residue_map': bool(order_ok),
        'ca_fit_rmsd_A': round(float(q['ca_rmsd'][0]), 4),
        'lig_rmsd_A': [round(float(q['lig_rmsd'][0]), 4), float(s['rmsd_to_crystal_A'])],
        'lig_rmsd_labelling': int(q['lig_rmsd_labelling'][0]),
        'd1_A': [round(float(q['d1'][0]), 3), float(s['N1_Phe131CZ'])],
        'd2_A': [round(float(q['d2'][0]), 3), float(s['N1_Asn67OD1'])],
        'state': str(st),
        'Zn_N2_A': [round(float(q['zn_n2'][0]), 3), float(s['Zn_N2'])],
        'Zn_O1_A': [round(float(q['zn_o1'][0]), 3), float(s['Zn_O1'])],
        'Zn_O2_A': [round(float(q['zn_o2'][0]), 3), float(s['Zn_O2'])],
        'N2_OG1_A': [round(float(q['n2_og1'][0]), 3), float(s['N2_OG1'])],
        'ang_N2_H_OG1_deg': [round(float(q['ang_n2_h_og1'][0]), 2), float(s['N2_H_OG1'])],
        'hb': bool(q['hb'][0]),
        'n2_signed_volume': round(float(q['n2_vol'][0]), 4),
        'Zn_H94NE2_A': round(float(q['zn_h94'][0]), 3), 'Zn_H96NE2_A': round(float(q['zn_h96'][0]), 3),
        'Zn_H119ND1_A': round(float(q['zn_h119'][0]), 3),
        'note': 'pairs are [this reader on solute.pdb, md/system/starts.csv]',
    }
out['starts'] = rows
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'check_starts.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
