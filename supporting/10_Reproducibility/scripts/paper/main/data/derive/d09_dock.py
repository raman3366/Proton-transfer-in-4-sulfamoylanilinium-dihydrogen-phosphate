"""d09: docking headline numbers re-derived from the raw run outputs, without analyse.py/production_analysis.py.

Scores: the first 'REMARK VINA RESULT' of every seedNN/out.pdbqt (the top pose); mean and SD (ddof 1) per job.
Ranks: by mean score among the production species of the same target and scoring function (no site filter).
RMSD: heavy atoms of the top pose in the receptor frame (no superposition), symmetry handled by enumerating all
element-preserving graph isomorphisms (RDKit substructure matching on distance-perceived, all-single-bond graphs)
between the pose and the crystal ligand (or the sulfanilamide core inside a larger crystal ligand).
"""
import glob
import re
import statistics as stt
from rdkit import Chem
from rdkit.Chem import rdchem
from common import *

AD2EL = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O', 'S': 'S', 'SA': 'S', 'HD': 'H',
         'H': 'H', 'F': 'F', 'Cl': 'Cl', 'CL': 'Cl', 'P': 'P'}
RC = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05, 'F': 0.57, 'Cl': 1.02}


def top_pose(pdbqt):
    """(score, [(el, xyz)]) of MODEL 1, heavy atoms only."""
    score, atoms = None, []
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith('REMARK VINA RESULT') and score is None:
            score = float(line.split()[3])
        if line.startswith(('ATOM', 'HETATM')):
            t = line[77:79].strip()
            el = AD2EL.get(t, t[0])
            if el != 'H':
                atoms.append((el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
        if line.startswith('ENDMDL'):
            break
    return score, atoms


def graph(atoms):
    m = Chem.RWMol()
    for el, _ in atoms:
        a = Chem.Atom(el); a.SetNoImplicit(True); m.AddAtom(a)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if dist(atoms[i][1], atoms[j][1]) < 1.2 * (RC[atoms[i][0]] + RC[atoms[j][0]]):
                m.AddBond(i, j, rdchem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def sdf_atoms(path):
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
    conf = mol.GetConformer()
    return [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms() if a.GetSymbol() != 'H']


def rmsd_sym(pose, ref, core=False):
    """min RMSD over all mappings of pose onto ref (core=True: pose graph as a substructure of ref)."""
    gp, gr = graph(pose), graph(ref)
    if core:
        matches = gr.GetSubstructMatches(gp, uniquify=False, maxMatches=100000)
        best = min(np.sqrt(np.mean([np.sum((pose[i][1] - ref[j][1]) ** 2) for i, j in enumerate(mt)])) for mt in matches)
    else:
        matches = gp.GetSubstructMatches(gr, uniquify=False, maxMatches=100000)
        best = min(np.sqrt(np.mean([np.sum((pose[i][1] - ref[j][1]) ** 2) for j, i in enumerate(mt)])) for mt in matches)
    return float(best), len(matches)


def job(target, tag):
    seeds = sorted(glob.glob(str(DOCK / target / 'runs' / tag / 'seed*' / 'out.pdbqt')))
    res = [top_pose(s) for s in seeds]
    sc = [r[0] for r in res]
    return dict(n=len(sc), mean=stt.mean(sc), sd=stt.stdev(sc) if len(sc) > 1 else 0.0, best=min(sc),
                per_seed=sc), [r[1] for r in res]


out = {'scores': {}, 'ranks': {}, 'rmsd': {}}
SPECIES = ['sfa_neutral', 'sfa_anion', 'sfa_cation', 'smx_neutral', 'smx_anion', 'azm_neutral', 'azm_anion',
           'paba_anion', 'dhp_anion']
for target in ['T5_7S2K', 'T4_3HS4', 'T3_6RL9', 'T1_3TZF']:
    for scoring in ['vina', 'vinardo', 'ad4zn']:
        rows = {}
        for sp in SPECIES:
            tag = f'{sp}__{scoring}__prod'
            if not (DOCK / target / 'runs' / tag).exists():
                continue
            s, _ = job(target, tag)
            rows[sp] = s
            out['scores'][f'{target}/{tag}'] = s
        if rows:
            order = sorted(rows, key=lambda k: rows[k]['mean'])
            out['ranks'][f'{target}/{scoring}'] = dict(order=[(k, round(rows[k]['mean'], 3)) for k in order],
                                                      n=len(order),
                                                      rank_sfa_neutral=order.index('sfa_neutral') + 1 if 'sfa_neutral' in order else None)

# RMSD re-derivations
checks = [
    ('T3_6RL9', 'sfa_neutral__ad4zn__prod', 'T3_6RL9/ref_SAN.sdf', False, 'production, DFT geometry, own co-crystal SAN'),
    ('T3_6RL9', 'sfa_neutral__ad4zn', 'T3_6RL9/ref_SAN.sdf', False, 'validation, ETKDG geometry'),
    ('T4_3HS4', 'sfa_neutral__ad4zn__prod', 'T4_3HS4/ref_cross_T3_6RL9_SAN.sdf', False, 'production, cross-dock vs 6RL9 SAN carried into 3HS4'),
    ('T4_3HS4', 'sfa_neutral__ad4zn', 'T4_3HS4/ref_cross_T3_6RL9_SAN.sdf', False, 'validation cross-dock'),
    ('T5_7S2K', 'dhp_anion__vina', 'T5_7S2K/ref_78H.sdf', False, 'validation redock 78H, Vina'),
    ('T5_7S2K', 'dhp_anion__vinardo', 'T5_7S2K/ref_78H.sdf', False, 'validation redock 78H, Vinardo'),
    ('T1_3TZF', 'smx_anion__vina', 'T1_3TZF/ref_08D.sdf', False, 'validation redock SMX anion, Vina'),
    ('T1_3TZF', 'sfa_neutral__vina__prod', 'T1_3TZF/ref_08D.sdf', True, 'production; sulfanilamide core inside 08D (SMX)'),
    ('T1_3TZF', 'sfa_neutral__vinardo__prod', 'T1_3TZF/ref_08D.sdf', True, 'production; core inside 08D, Vinardo'),
]
for target, tag, ref, core, what in checks:
    s, poses = job(target, tag)
    refat = sdf_atoms(DOCK / ref)
    r = [rmsd_sym(p, refat, core) for p in poses]
    vals = [x[0] for x in r]
    out['rmsd'][f'{target}/{tag}'] = dict(what=what, ref=f'dock/{ref}', core=core, n=len(vals), mean=stt.mean(vals),
                                         sd=stt.stdev(vals), min=min(vals), max=max(vals),
                                         f_le_1=sum(v <= 1.0 for v in vals) / len(vals),
                                         f_le_2=sum(v <= 2.0 for v in vals) / len(vals), n_mappings=r[0][1],
                                         score=s)
# compare with production_summary.json
ps = json.loads((DOCK / 'results/production_summary.json').read_text())
cmp_ = {}
for k, v in out['scores'].items():
    j = ps['jobs'].get(k)
    if j:
        cmp_[k] = dict(d_mean=v['mean'] - j['score']['mean'], d_sd=v['sd'] - j['score']['sd'])
out['max_abs_score_mean_diff_vs_summary'] = max(abs(c['d_mean']) for c in cmp_.values())
out['max_abs_score_sd_diff_vs_summary'] = max(abs(c['d_sd']) for c in cmp_.values())
out['n_jobs_compared'] = len(cmp_)
dump('d09_dock.json', out)
