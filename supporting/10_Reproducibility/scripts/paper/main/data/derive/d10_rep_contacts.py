"""d10: distances of the docking cluster representatives drawn in Fig. 10, re-derived from the pose files.

For the four sulfanilamide production jobs of Fig. 10 (3TZF and 7S2K with Vina, 6RL9 and 3HS4 with AutoDock4Zn) the
representative is the top pose of the seed that dock/results/production_contacts.json names for cluster 1.  Here the
pose of that seed is read from dock/results/poses/<target>/<job>__top_pose_per_seed.sdf (the <seed> tag, not the
record order) and the receptor from dock/<target>/receptor.pdb, and every distance is recomputed with numpy:
  * each polar pair listed in production_contacts.json (receptor residue + atom name, ligand atom index '#k');
  * the zinc to the nearest ligand heavy atom (hCA II), and which atom that is;
  * as an independent sweep, every receptor N/O within 3.5 A of a ligand N/O (to see that none is missing).
Plain-text parsing only; no dock/ code is imported (derive/ convention).  Output: derived/d10_rep_contacts.json.
"""
from common import *

JOBS = ['T1_3TZF/sfa_neutral__vina__prod', 'T5_7S2K/sfa_neutral__vina__prod',
        'T3_6RL9/sfa_neutral__ad4zn__prod', 'T4_3HS4/sfa_neutral__ad4zn__prod']


def sdf_records(path):
    """[(tags, [(el, xyz)])] of every record of an SDF file (V2000 atom block)."""
    out = []
    for rec in Path(path).read_text().split('$$$$'):
        lines = rec.strip('\n').splitlines()
        if len(lines) < 4:
            continue
        n = int(lines[3][:3])
        atoms = [(ln[31:34].strip(), np.array([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])]))
                 for ln in lines[4:4 + n]]
        tags, i = {}, 4 + n
        while i < len(lines):
            if lines[i].startswith('>'):
                name = lines[i].split('<')[1].split('>')[0]
                tags[name] = lines[i + 1].strip() if i + 1 < len(lines) else ''
                i += 2
            else:
                i += 1
        out.append((tags, atoms))
    return out


def receptor(path):
    rows = []
    for ln in Path(path).read_text().splitlines():
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        el = ln[76:78].strip() or ln[12:16].strip()[0]
        rows.append((f"{ln[17:20].strip()}{int(ln[22:26])}", ln[12:16].strip(), el.capitalize(),
                     np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return rows


pc = json.loads((DOCK / 'results/production_contacts.json').read_text())
res = {}
for job in JOBS:
    target, tag = job.split('/')
    c1 = pc[job]['1']
    seed = int(c1['seed'])
    stem = tag.replace('__prod', '')
    recs = [r for r in sdf_records(DOCK / 'results/poses' / target / f'{stem}__top_pose_per_seed.sdf') if int(r[0].get('seed', -1)) == seed]
    assert len(recs) == 1, (job, seed, len(recs))
    tags, lig = recs[0]
    rec = receptor(DOCK / target / 'receptor.pdb')
    pairs = []
    for p_ in c1['polar_pairs_NO_within_3.5A']:
        k = int(p_['ligand_atom'].split('#')[1])
        ra = [x for x in rec if x[0] == p_['residue'] and x[1] == p_['receptor_atom']]
        assert len(ra) == 1, (job, p_)
        d = dist(ra[0][3], lig[k][1])
        pairs.append(dict(residue=p_['residue'], receptor_atom=p_['receptor_atom'], ligand_atom=p_['ligand_atom'],
                          ligand_element=lig[k][0], distance=d, source_value=p_['distance'], agrees=abs(round(d, 2) - p_['distance']) < 0.006))
    heavy = [(i, el, x) for i, (el, x) in enumerate(lig) if el != 'H']
    sweep = sorted([(x[0], x[1], i, dist(x[3], xyz)) for x in rec if x[2] in ('N', 'O') for i, el, xyz in heavy
                    if el in ('N', 'O') and dist(x[3], xyz) <= 3.5], key=lambda t: t[3])
    zn = None
    zr = [x for x in rec if x[2] == 'Zn']
    if zr:
        i, el, xyz = min(heavy, key=lambda t: dist(zr[0][3], t[2]))
        zn = dict(nearest_ligand_index=i, element=el, distance=dist(zr[0][3], xyz), source_value=(c1.get('zinc') or {}).get('distance'))
        zn['agrees'] = zn['source_value'] is not None and abs(round(zn['distance'], 2) - zn['source_value']) < 0.006
    res[job] = dict(seed=seed, score_tag=float(tags.get('score_kcal_mol', 'nan')), score_source=c1['score'],
                    cluster_tag=tags.get('cluster'), population=c1['population'], pairs=pairs,
                    sweep_NO_within_3p5=[dict(residue=a, receptor_atom=b, ligand_index=i, distance=d) for a, b, i, d in sweep],
                    n_sweep=len(sweep), n_listed=len(pairs), zinc=zn)

out = dict(what=__doc__.strip().split('\n')[0], jobs=res,
           all_pairs_agree=all(p['agrees'] for j in res.values() for p in j['pairs']),
           all_zinc_agree=all(j['zinc']['agrees'] for j in res.values() if j['zinc']))
(OUT / 'd10_rep_contacts.json').write_text(json.dumps(out, indent=1))
print('d10:', {j.split('/')[0]: (v['seed'], v['n_listed'], v['n_sweep'], v['zinc'] and round(v['zinc']['distance'], 3)) for j, v in res.items()},
      'pairs agree', out['all_pairs_agree'], 'zinc agree', out['all_zinc_agree'])
