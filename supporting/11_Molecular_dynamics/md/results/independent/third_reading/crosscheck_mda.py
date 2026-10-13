#!/usr/bin/env python
"""Cross-check of the third reader's own DCD reader and of d1, d2, the ligand
RMSD and the hydrogen bond against MDAnalysis (its DCD reader, its selection
language, its align/rms functions), on a few files.  Still the third reader's
own code; nothing of the two analysis programs is used."""
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import calc_angles

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import third_lib as T

HERE = os.path.dirname(os.path.abspath(__file__))
ref = T.build_reference()
rmap = T.residue_map()
res = {}
cases = [('bonded_2e', 'X1', 'solute_seg0001.dcd'), ('bonded_2e', 'R1', 'solute_eq04.dcd'),
         ('bonded_2e', 'R3', 'solute_seg0020.dcd'), ('bonded_2s', 'R1', 'solute_seg0009.dcd'),
         ('bonded_2e_n2fit', 'R1', 'solute_seg0015.dcd'), ('bonded_2e_n2fit', 'R1', 'solute_seg0020.dcd')]
for zset, run, fn in cases:
    folder = os.path.join(T.RUNS, zset, run)
    sel = T.select(os.path.join(folder, 'solute.pdb'))
    xyz, lat, head = T.read_dcd(os.path.join(folder, fn), sel['n_atoms'])
    q = T.frame_quantities(xyz, lat, sel, ref)
    u = mda.Universe(os.path.join(folder, 'solute.pdb'), os.path.join(folder, fn))
    n1 = u.select_atoms('resname SN1 and name N1')
    cz = u.select_atoms('resid %d and name CZ' % rmap[131][0])
    od1 = u.select_atoms('resid %d and name OD1' % rmap[67][0])
    n2 = u.select_atoms('resname SN1 and name N2')
    hn = u.select_atoms('resname SN1 and name HN21')
    og1 = u.select_atoms('resid %d and name OG1' % rmap[199][0])
    ca = u.select_atoms('protein and name CA') if False else u.select_atoms('name CA and not resname SN1 ZN1')
    ligh = u.select_atoms('resname SN1 and name ' + ' '.join(T.LIG_HEAVY))
    order = [list(ligh.names).index(n) for n in T.LIG_HEAVY]
    refl = np.array([ref['lig'][n] for n in T.LIG_HEAVY])
    assert len(n1) == len(cz) == len(od1) == len(n2) == len(hn) == len(og1) == 1 and len(ca) == 258
    maxdx = 0.0
    d1 = []; d2 = []; rm = []; hb = []; dims = []
    for i, ts in enumerate(u.trajectory):
        maxdx = max(maxdx, float(np.abs(ts.positions - xyz[i]).max()))
        dims.append(ts.dimensions.copy())
        d1.append(np.linalg.norm(n1.positions[0] - cz.positions[0]))
        d2.append(np.linalg.norm(n1.positions[0] - od1.positions[0]))
        dn = np.linalg.norm(n2.positions[0] - og1.positions[0])
        ang = np.degrees(calc_angles(n2.positions, hn.positions, og1.positions))[0]
        hb.append(bool(dn <= 3.5 and ang >= 135.0))
        # fit by MDAnalysis: rotation matrix from the centred CA sets
        mob = ca.positions.astype(float)
        R, rmsd_ca = align.rotation_matrix(mob - mob.mean(0), ref['ca'] - ref['ca'].mean(0))
        lp = (ligh.positions.astype(float)[order] - mob.mean(0)) @ R.T + ref['ca'].mean(0)
        best = 1e9
        for _, perm in T.LABELLINGS:
            target = np.array([ref['lig'][perm.get(n, n)] for n in T.LIG_HEAVY])
            best = min(best, float(np.sqrt(((lp - target) ** 2).sum(1).mean())))
        rm.append(best)
    dims = np.array(dims)
    mylen = np.linalg.norm(lat, axis=2)
    res['%s/%s/%s' % (zset, run, fn)] = {
        'n_frames_mine': int(xyz.shape[0]), 'n_frames_mdanalysis': int(len(u.trajectory)),
        'max_abs_coordinate_difference_A': maxdx,
        'max_abs_cell_length_difference_A': float(np.abs(dims[:, :3] - mylen).max()),
        'cell_angles_mdanalysis_deg_first_frame': [float(x) for x in dims[0, 3:]],
        'max_abs_d1_difference_A': float(np.abs(np.array(d1) - q['d1']).max()),
        'max_abs_d2_difference_A': float(np.abs(np.array(d2) - q['d2']).max()),
        'max_abs_lig_rmsd_difference_A': float(np.abs(np.array(rm) - q['lig_rmsd']).max()),
        'hbond_flags_differ_in_frames': int(np.sum(np.array(hb) != q['hb'])),
    }
with open(os.path.join(HERE, 'crosscheck_mda.json'), 'w') as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
