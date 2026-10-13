"""
Production docking analysis (called as `python dock/analyse.py --production`).

Reads every production run (dock/T*/runs/<ligand>__<scoring>__prod/results.json) and writes
dock/results/production_summary.json, plus the pose files under dock/results/poses/ and the
geometries for later DFT single points under dock/results/poses_for_dft/.

Per job (target x scoring function x species):

*   **Score** of the top pose of every seed: mean +- SD (ddof 1), median, best, worst.
*   **Ligand efficiency** LE = -score / heavy atoms, and the size-corrected fit quality
    FQ = LE / LE_scale(HA) with LE_scale = 0.0715 + 7.5328/HA + 25.7079/HA^2 - 361.4722/HA^3
    (Reynolds, Tounge & Bembenek 2007).  LE has a known size bias: scores grow more slowly than
    linearly with size, so LE rewards small ligands, and sulfanilamide (11 heavy atoms) is the
    smallest molecule here after pABA.  The bias is measured on this data set too (LE vs HA).
    LE_scale was fitted to experimental affinities, so applied to docking scores FQ is a heuristic.
*   **Clusters** of the seeds' top poses: leader clustering in score order on the
    symmetry-corrected heavy-atom RMSD (spyrmsd, no superposition), cutoff 2.0 A (the AutoDock
    convention).  Population = number of seeds whose top pose falls in the cluster.
*   **RMSD to the co-crystal ligand** where the docked species is the same molecule (same heavy-atom
    element graph) as a crystal ligand of that receptor, or of the other hCA II structure carried
    over by Ca superposition (the cross-docking references of analyse.py).
*   **Core RMSD** where only a moiety is shared: the aniline-SO2N core (11 atoms: sulfanilamide
    in full, the sulfanilyl part of sulfamethoxazole), the aniline-X(O)2 core with X = S or C (10
    atoms: sulfanilamide against the p-aminobenzoyl part of dihydropteroate), and the C-SO2N
    zinc-binding group (5 atoms: sulfanilamide against acetazolamide).  Every chemically equivalent
    atom mapping is tried (ring flip, oxygen swap) and the smallest RMSD is kept.  In the receptor
    frame, no superposition.
*   **Site checks** of the validation (analyse.SITE_CHECKS) plus two generalised ones for species
    the validation SMARTS do not match, each also evaluated on the crystal ligand.
*   **Rank order** of all species per target x scoring function, with the co-crystal ligand as the
    reference score, and for hCA II also over zinc-bound top poses only.
*   **DFT-docking link** for the three sulfanilamide species: every top pose is rebuilt as the full
    DFT molecule (every hydrogen, the DFT bond lengths and angles, the docked torsions) and checked
    to reproduce the docked atoms; the C3-C4-S1-N2 torsion is measured on it (atom labels of the
    crystal, carried through the DFT xyz), and an MMFF94s strain estimate is computed two ways
    (rigid DFT geometry with the polar-H rotors re-optimised; MMFF-relaxed with the torsion fixed).
"""
import hashlib
import json
import math
import pathlib
import shutil
import sys
import time

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolHash
from rdkit.Chem import rdMolTransforms as T

RDLogger.DisableLog('rdApp.warning')
DOCK = pathlib.Path(__file__).resolve().parent
ROOT = DOCK.parent
RES = DOCK / 'results'
sys.path.insert(0, str(DOCK))
import analyse as A  # noqa: E402
from prep_ligands import SPECIES, from_dft_xyz  # noqa: E402

SUFFIX = '__prod'
# diagnostic runs: same analysis, kept out of the rankings and of the DFT pose set.  __csn102 /
# __csn109 are the C4-S1-N2 angle test of the input-geometry sensitivity found in production
# (dock/geometry_variant.py), __etkdg re-runs a production job from the validation conformer.
DIAG_SUFFIXES = ('__csn102', '__csn109', '__etkdg')
CLUSTER_CUTOFF = 2.0
# PoseBusters on the production poses; --no-posebusters (analyse.py) turns it off (review 2026-09-28)
POSEBUSTERS = True
TARGET_ORDER = ['T5_7S2K', 'T4_3HS4', 'T3_6RL9', 'T1_3TZF']
SPECIES_ORDER = ['sfa_neutral', 'sfa_anion', 'sfa_cation', 'smx_anion', 'smx_neutral',
                 'paba_anion', 'azm_neutral', 'azm_anion', 'dhp_anion']
SCORING_ORDER = ['vina', 'vinardo', 'ad4zn']
TARGET_INFO = {
    'T1_3TZF': dict(pdb='3TZF', protein='DHPS, Y. pestis (chain B)', crystal_ligand='08D',
                    crystal_ligand_name='sulfamethoxazole', reference_species=['smx_anion', 'smx_neutral'],
                    validation='passes the 2.0 \u00c5 criterion only (SMX anion, Vina, 1.92 \u00b1 0.04 \u00c5); poses qualitative'),
    'T3_6RL9': dict(pdb='6RL9', protein='hCA II', crystal_ligand='SAN',
                    crystal_ligand_name='sulfanilamide', reference_species=['sfa_neutral'],
                    validation='partially validated with AD4Zn (2.07 \u00b1 1.45 \u00c5, f\u22641 \u00c5 0.5, bimodal: fails the 1.0 \u00c5 gate)'),
    'T4_3HS4': dict(pdb='3HS4', protein='hCA II', crystal_ligand='AZM',
                    crystal_ligand_name='acetazolamide', reference_species=['azm_neutral', 'azm_anion'],
                    validation='validated cross-dock for sulfanilamide with AD4Zn (0.62 \u00b1 0.03 \u00c5, 10/10); the AZM pose itself is not reproduced'),
    'T5_7S2K': dict(pdb='7S2K', protein='Sul2 (sulfonamide-resistant DHPS)', crystal_ligand='78H',
                    crystal_ligand_name='7,8-dihydropteroate', reference_species=['dhp_anion'],
                    validation='validated with Vina and Vinardo (78H redock 0.51 / 0.54 \u00c5, 10/10 seeds within 1 \u00c5)'),
}
OWN_REF = {'T1_3TZF': 'ref_08D.sdf', 'T3_6RL9': 'ref_SAN.sdf', 'T4_3HS4': 'ref_AZM.sdf',
           'T5_7S2K': 'ref_78H.sdf'}
CROSS_REF = {'T3_6RL9': ('T4_3HS4', 'ref_AZM.sdf'), 'T4_3HS4': ('T3_6RL9', 'ref_SAN.sdf')}

CORES = [
    ('aniline-SO2N core (11 atoms)',
     '[#7;!$([#7]-[#16]);!$([#7]~[#8])]-[#6]1~[#6]~[#6]~[#6](~[#6]~[#6]~1)-[#16](~[#8])(~[#8])~[#7]'),
    ('aniline-X(O)2 core, X = S or C (10 atoms)',
     '[#7;!$([#7]-[#16]);!$([#7]~[#8])]-[#6]1~[#6]~[#6]~[#6](~[#6]~[#6]~1)-[#6,#16](~[#8])~[#8]'),
    ('C-SO2N zinc-binding group (5 atoms)', '[#6]-[#16](~[#8])(~[#8])~[#7]'),
]
# Checks for species the validation SMARTS were not written for (the cation's NH3+, the sulfonyl
# oxygens standing in for pABA's carboxylate).  Evaluated on the crystal ligand as well.
EXTRA_CHECKS = {
    'T1_3TZF': [('aniline N (any protonation) to pterin HH2 N/O',
                 '[#7;$([#7]-[#6]1~[#6]~[#6]~[#6](~[#6]~[#6]~1)-[#6,#16]);!$([#7]-[#16]);!$([#7]~[#8])]',
                 ('HH2', None, '*NO'), 3.6)],
    'T5_7S2K': [('sulfonyl or carboxylate O to Ser214 N',
                 '[#8;$([#8]~[#16X4](~[#8])),$([#8]~[#6X3](~[#8])-c)]', ('SER', 214, 'N'), 3.4)],
}
# the DFT atom labels (the xyz files follow the crystal labels C1..C6 N1 N2 O1 O2 S1)
LBL = dict(C1=0, C3=2, C4=3, C5=4, N1=6, N2=7, S1=10)
DFT_XYZ = {'sfa_neutral': 'molecules/sfa/results/sulfanilamide_pcm_d3',
           'sfa_anion': 'molecules/sfa/results/sulfonamidate_pcm_d3',
           'sfa_cation': 'molecules/sfa/results/cation_pcm_d3'}
XRAY_XYZ = ROOT / 'structures' / 'sfa_pair_xray.xyz'


def reynolds_le_scale(ha):
    return 0.0715 + 7.5328 / ha + 25.7079 / ha ** 2 - 361.4722 / ha ** 3


def rnd(x, n=3):
    return None if x is None else round(float(x), n)


def stats(v):
    v = np.asarray(v, float)
    if not len(v):
        return None
    return dict(n=int(len(v)), mean=rnd(v.mean()), sd=rnd(v.std(ddof=1)) if len(v) > 1 else 0.0,
                median=rnd(np.median(v)), min=rnd(v.min()), max=rnd(v.max()))


def heavy_xyz(mol, cid=-1):
    c = mol.GetConformer(cid)
    return np.array([list(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def element_graph(mol):
    return rdMolHash.MolHash(Chem.RemoveHs(mol), rdMolHash.HashFunction.ElementGraph)


# ------------------------------------------------------------------------------------------------
# references
# ------------------------------------------------------------------------------------------------
_REFS = {}


def references(target):
    """[(label, mol, info)] for the crystal ligands of this receptor (own, and cross for hCA II)."""
    if target in _REFS:
        return _REFS[target]
    tdir = DOCK / target
    own = Chem.MolFromMolFile(str(tdir / OWN_REF[target]))
    out = [(f'{TARGET_INFO[target]["crystal_ligand"]} ({TARGET_INFO[target]["pdb"]}, own)', own,
            dict(kind='own', file=str((tdir / OWN_REF[target]).relative_to(DOCK))))]
    if target in CROSS_REF:
        src, fn = CROSS_REF[target]
        sref = Chem.MolFromMolFile(str(DOCK / src / fn))
        (RES / 'refs').mkdir(parents=True, exist_ok=True)
        dst = RES / 'refs' / f'{target}__{src}_{fn[:-4]}_in_{target}_frame.sdf'
        moved, ca_rmsd, n_ca = A.transform_reference(DOCK / src, tdir, sref, dst)
        out.append((f'{TARGET_INFO[src]["crystal_ligand"]} ({TARGET_INFO[src]["pdb"]}, carried into '
                    f'{TARGET_INFO[target]["pdb"]} by Ca superposition)', moved,
                    dict(kind='cross', source=src, file=str(dst.relative_to(DOCK)),
                         superposition_ca_rmsd_A=ca_rmsd, n_ca_atoms=n_ca)))
    _REFS[target] = out
    return out


def core_rmsd_values(ref, lig, lig_confs, smarts):
    """Smallest RMSD over every mapping of the core onto ligand and reference (no superposition)."""
    patt = Chem.MolFromSmarts(smarts)
    rm = ref.GetSubstructMatches(patt, uniquify=False, maxMatches=1000)
    lm = lig.GetSubstructMatches(patt, uniquify=False, maxMatches=1000)
    if not rm or not lm:
        return None
    R = heavy_xyz(ref)
    out = []
    for cid in lig_confs:
        X = heavy_xyz(lig, cid)
        best = min(float(np.sqrt(((X[list(a)] - R[list(b)]) ** 2).sum(1).mean()))
                   for a in lm for b in rm)
        out.append(best)
    return out, patt.GetNumAtoms()


def rmsd_block(vals, n_atoms=None):
    v = np.asarray(vals, float)
    d = dict(per_seed=[rnd(x) for x in v], **stats(v),
             f_le_1A=rnd((v <= 1.0).mean()), f_le_2A=rnd((v <= 2.0).mean()))
    if n_atoms:
        d['n_atoms'] = n_atoms
    return d


# ------------------------------------------------------------------------------------------------
# clustering of the seeds' top poses
# ------------------------------------------------------------------------------------------------
def pairwise_symm_rmsd(mol, cids):
    from spyrmsd import rmsd as srmsd
    m, z, adj = A.heavy_graph(mol)
    # heavy_graph drops the hydrogens; RemoveHs keeps the conformers, so index them the same way
    coords = [heavy_xyz(m, c) for c in cids]
    n = len(cids)
    D = np.zeros((n, n))
    for i in range(n):
        if i + 1 < n:
            vals = srmsd.symmrmsd(coords[i], coords[i + 1:], z, z, adj, adj, center=False,
                                  minimize=False)
            vals = vals if isinstance(vals, list) else [vals]
            D[i, i + 1:] = vals
            D[i + 1:, i] = vals
    return D


def leader_clusters(scores, D, cutoff=CLUSTER_CUTOFF):
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    leaders, members = [], {}
    for i in order:
        for L in leaders:
            if D[i, L] <= cutoff:
                members[L].append(i)
                break
        else:
            leaders.append(i)
            members[i] = [i]
    return [(L, members[L]) for L in leaders]


# ------------------------------------------------------------------------------------------------
# DFT reconstruction of sulfanilamide poses, torsions and MMFF strain
# ------------------------------------------------------------------------------------------------
def parse_pdbqt_atoms(text, model=1):
    """Coordinates (in file order) of one model, plus BRANCH pairs (1-based serials)."""
    cur, m, branches, coords = 0, None, [], []
    has_models = 'MODEL' in text
    for ln in text.splitlines():
        if ln.startswith('MODEL'):
            cur += 1
            continue
        if has_models and cur != model:
            continue
        if ln.startswith(('ATOM', 'HETATM')):
            coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        elif ln.startswith('BRANCH'):
            a, b = ln.split()[1:3]
            branches.append((int(a), int(b)))
        elif ln.startswith('ENDMDL') and has_models:
            break
    return np.array(coords), branches


_DFT = {}


def dft_mol(species):
    if species not in _DFT:
        path = ROOT / DFT_XYZ[species] / 'opt.xyz'
        m = from_dft_xyz(path, SPECIES[species][0], species)
        c = m.GetConformer()
        # the labels must mean what the crystal means: C4 bonded to S1, N2 bonded to S1, C3 to C4
        b = lambda i, j: m.GetBondBetweenAtoms(i, j) is not None  # noqa: E731
        assert m.GetAtomWithIdx(LBL['S1']).GetSymbol() == 'S' and b(LBL['C4'], LBL['S1']) \
            and b(LBL['S1'], LBL['N2']) and b(LBL['C3'], LBL['C4']) and b(LBL['C5'], LBL['C4']) \
            and b(LBL['C1'], LBL['N1']), f'{species}: DFT atom order is not the crystal order'
        summ = json.loads((ROOT / DFT_XYZ[species] / 'opt_summary.json').read_text())
        _DFT[species] = dict(mol=m, xyz=str(path.relative_to(ROOT)),
                             E_Ha=summ.get('E_final_Ha'),
                             tau_C3=rnd(T.GetDihedralDeg(c, LBL['C3'], LBL['C4'], LBL['S1'], LBL['N2']), 2),
                             tau_C5=rnd(T.GetDihedralDeg(c, LBL['C5'], LBL['C4'], LBL['S1'], LBL['N2']), 2))
    return _DFT[species]


def dft_vs_etkdg(species):
    """How far the DFT input geometry is from the ETKDG/MMFF94s conformer the validation runs used:
    bond-length and bond-angle deviations over the matched heavy-atom skeleton, and the torsion.
    This is the size of the input change whose effect on the scores section 7 of the report measures."""
    from prep_ligands import build_species, embed
    d = dft_mol(species)['mol']
    e = embed(build_species(species, SPECIES[species], {}), seed=1)
    dm, em = Chem.RemoveHs(d), Chem.RemoveHs(e)
    match = em.GetSubstructMatch(dm)          # ETKDG index for each DFT index
    if not match:
        return None
    dc, ec = dm.GetConformer(), em.GetConformer()
    bl, ang = [], []
    for b in dm.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bl.append(T.GetBondLength(dc, i, j) - T.GetBondLength(ec, match[i], match[j]))
    for a in dm.GetAtoms():
        nb = [x.GetIdx() for x in a.GetNeighbors()]
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                ang.append(T.GetAngleDeg(dc, nb[i], a.GetIdx(), nb[j])
                           - T.GetAngleDeg(ec, match[nb[i]], match[a.GetIdx()], match[nb[j]]))
    bl, ang = np.array(bl), np.array(ang)
    return dict(n_bonds=len(bl), bond_rmsd_A=rnd(np.sqrt((bl ** 2).mean()), 4),
                bond_max_abs_A=rnd(np.abs(bl).max(), 4), n_angles=len(ang),
                angle_rmsd_deg=rnd(np.sqrt((ang ** 2).mean()), 2),
                angle_max_abs_deg=rnd(np.abs(ang).max(), 2),
                tau_C3_etkdg=rnd(T.GetDihedralDeg(em.GetConformer(), match[LBL['C3']],
                                                  match[LBL['C4']], match[LBL['S1']],
                                                  match[LBL['N2']]), 2))


def csn_angles():
    """The C-S-N angle of the sulfonamide in every geometry that matters, measured in this pass.

    Vina keeps bond angles rigid, and this is the coordinate in which the ETKDG/MMFF94s input and the
    DFT input differ most, so it is the one the AD4Zn ring-orientation diagnostic tests."""
    from prep_ligands import build_species, embed

    def measure(mol, tag):
        c = mol.GetConformer()
        for a in mol.GetAtoms():
            if a.GetSymbol() != 'S':
                continue
            nb = a.GetNeighbors()
            ns = [x.GetIdx() for x in nb if x.GetSymbol() == 'N']
            cs = [x.GetIdx() for x in nb if x.GetSymbol() == 'C']
            os_ = [x.GetIdx() for x in nb if x.GetSymbol() == 'O']
            if ns and cs and len(os_) == 2:
                return dict(source=tag, C_S_N_deg=rnd(T.GetAngleDeg(c, cs[0], a.GetIdx(), ns[0]), 2),
                            O_S_O_deg=rnd(T.GetAngleDeg(c, os_[0], a.GetIdx(), os_[1]), 2),
                            S_N_A=rnd(T.GetBondLength(c, a.GetIdx(), ns[0]), 3))
        return None

    out = dict(dft={}, etkdg={}, crystal={})
    for sp in DFT_XYZ:
        out['dft'][sp] = measure(dft_mol(sp)['mol'], dft_mol(sp)['xyz'])
    for seed in (1, 2, 3):
        out['etkdg'][f'seed{seed}'] = measure(
            embed(build_species('sfa_neutral', SPECIES['sfa_neutral'], {}), seed=seed),
            f'ETKDGv3+MMFF94s sfa_neutral seed {seed}')
    for t in TARGET_ORDER:
        ref = Chem.MolFromMolFile(str(DOCK / t / OWN_REF[t]))
        m = measure(ref, f'{TARGET_INFO[t]["pdb"]} {TARGET_INFO[t]["crystal_ligand"]}')
        if m:
            out['crystal'][TARGET_INFO[t]['pdb']] = m
    # the compound's own crystal structure: the xyz follows the crystal labels
    x = Chem.MolFromXYZFile(str(XRAY_XYZ))
    xc = x.GetConformer()
    out['crystal']['SAPH (this compound, CIF)'] = dict(
        source=str(XRAY_XYZ.relative_to(ROOT)) + ' (cation of the ion pair)',
        C_S_N_deg=rnd(T.GetAngleDeg(xc, LBL['C4'], LBL['S1'], LBL['N2']), 2),
        O_S_O_deg=rnd(T.GetAngleDeg(xc, 8, LBL['S1'], 9), 2),
        S_N_A=rnd(T.GetBondLength(xc, LBL['S1'], LBL['N2']), 3))
    return out


def crystal_torsion():
    x = Chem.MolFromXYZFile(str(XRAY_XYZ))
    c = x.GetConformer()
    return dict(source=str(XRAY_XYZ.relative_to(ROOT)) + ' (CIF coordinates; CIF _geom_torsion '
                'C3-C4-S1-N2 102.68(19), C5-C4-S1-N2 -77.4(2))',
                tau_C3=rnd(T.GetDihedralDeg(c, 2, 3, 10, 7), 2),
                tau_C5=rnd(T.GetDihedralDeg(c, 4, 3, 10, 7), 2))


def side_atoms(mol, a, b):
    """Atoms on b's side of the a-b bond."""
    seen, stack = {b}, [b]
    while stack:
        i = stack.pop()
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j != a and j not in seen:
                seen.add(j)
                stack.append(j)
    return sorted(seen)


def rotate(X, axis_a, axis_b, idx, deg):
    k = X[axis_b] - X[axis_a]
    k = k / np.linalg.norm(k)
    th = math.radians(deg)
    P = X[idx] - X[axis_b]
    P = P * math.cos(th) + np.cross(k, P) * math.sin(th) + np.outer(P @ k, k) * (1 - math.cos(th))
    Y = X.copy()
    Y[idx] = P + X[axis_b]
    return Y


def rebuild_pose(species, lig_pdbqt, out_pdbqt, model=1):
    """The docked pose as the full DFT molecule.  Vina keeps every bond length, bond angle and
    rigid-fragment geometry of its input, so setting the DFT molecule's rotatable torsions to the
    docked values and superposing must reproduce the docked atoms; the residual is returned."""
    d = dft_mol(species)
    m = Chem.Mol(d['mol'])
    D = heavy_xyz(m)
    inp, branches = parse_pdbqt_atoms(pathlib.Path(lig_pdbqt).read_text())
    pose, _ = parse_pdbqt_atoms(pathlib.Path(out_pdbqt).read_text(), model=model)
    if len(inp) != len(pose):
        raise ValueError('pose and input have different atom counts')
    amap = []
    for p in inp:
        dist = np.linalg.norm(D - p, axis=1)
        j = int(np.argmin(dist))
        if dist[j] > 0.002:
            raise ValueError(f'input atom at {p} is not a DFT atom (nearest {dist[j]:.4f} A)')
        amap.append(j)
    present = set(amap)
    rot = {(amap[a - 1], amap[b - 1]) for a, b in branches}
    rotset = rot | {(b, a) for a, b in rot}
    conf = m.GetConformer()
    P = {amap[i]: pose[i] for i in range(len(amap))}

    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1n = b1 / np.linalg.norm(b1)
        v = b0 - (b0 @ b1n) * b1n
        w = b2 - (b2 @ b1n) * b1n
        return math.degrees(math.atan2(np.cross(b1n, v) @ w, v @ w))

    for a, b in rot:
        xs = [n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors()
              if n.GetIdx() != b and n.GetIdx() in present and (a, n.GetIdx()) not in rotset]
        ys = [n.GetIdx() for n in m.GetAtomWithIdx(b).GetNeighbors()
              if n.GetIdx() != a and n.GetIdx() in present and (b, n.GetIdx()) not in rotset]
        if not xs or not ys:
            raise ValueError(f'no reference atoms for torsion {a}-{b}')
        x, y = xs[0], ys[0]
        T.SetDihedralDeg(conf, x, a, b, y, dih(P[x], P[a], P[b], P[y]))
    X = heavy_xyz(m)
    idx = sorted(P)
    Rm, t = A.kabsch(X[idx], np.array([P[i] for i in idx]))
    Y = (Rm @ X.T).T + t
    for i in range(m.GetNumAtoms()):
        conf.SetAtomPosition(i, Y[i].tolist())
    resid = float(np.sqrt(((Y[idx] - np.array([P[i] for i in idx])) ** 2).sum(1).mean()))
    return m, resid, rot


class Strain:
    """MMFF94s strain of the C3-C4-S1-N2 torsion for one species, two ways.

    rigid:   the DFT bond lengths and angles held fixed (as they are in the docked pose); only the
             torsions change.  E(pose, polar-H rotors re-optimised) minus the minimum of the same
             rigid profile over C3-C4-S1-N2.  The polar-H rotors are re-optimised because Vina does
             not score hydrogen positions, so their docked orientation carries no information.
    relaxed: MMFF94s minimisation with only C3-C4-S1-N2 held at the pose value, minus the MMFF94s
             global minimum (best of 36 unconstrained minimisations started around the torsion)."""

    def __init__(self, species):
        d = dft_mol(species)
        self.species = species
        self.m = Chem.Mol(d['mol'])
        self.mp = AllChem.MMFFGetMoleculeProperties(self.m, mmffVariant='MMFF94s')
        self.ff = AllChem.MMFFGetMoleculeForceField(self.m, self.mp)
        m = self.m
        self.tors = (LBL['C3'], LBL['C4'], LBL['S1'], LBL['N2'])
        # polar-H rotors: bonds from a heavy atom to an N carrying H, other than the C4-S1 torsion
        self.rotors = []
        for n_idx in (LBL['N1'], LBL['N2']):
            n = m.GetAtomWithIdx(n_idx)
            hs = [x.GetIdx() for x in n.GetNeighbors() if x.GetSymbol() == 'H']
            heavy = [x.GetIdx() for x in n.GetNeighbors() if x.GetSymbol() != 'H']
            if hs and len(heavy) == 1:
                self.rotors.append((heavy[0], n_idx, hs))
        self.s_side = side_atoms(m, LBL['C4'], LBL['S1'])

    def energy(self, X):
        return self.ff.CalcEnergy(X.ravel().tolist())

    def rotor_opt(self, X, step=10):
        """Grid-optimise the polar-H rotors (all combinations), then refine at 2 degrees."""
        grids = [np.arange(0, 360, step)] * len(self.rotors)
        best = (self.energy(X), X)
        for combo in np.array(np.meshgrid(*grids)).T.reshape(-1, len(self.rotors)):
            Y = X
            for (a, b, hs), ang in zip(self.rotors, combo):
                Y = rotate(Y, a, b, hs, ang)
            e = self.energy(Y)
            if e < best[0]:
                best = (e, Y)
        X0 = best[1]
        for _ in range(2):
            for k, (a, b, hs) in enumerate(self.rotors):
                for ang in np.arange(-10, 10.1, 2):
                    Y = rotate(X0, a, b, hs, ang)
                    e = self.energy(Y)
                    if e < best[0] - 1e-9:
                        best = (e, Y)
                X0 = best[1]
        return best

    def set_tau(self, X, tau):
        c3, c4, s1, n2 = self.tors
        m = Chem.Mol(self.m)
        conf = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, X[i].tolist())
        T.SetDihedralDeg(conf, c3, c4, s1, n2, float(tau))
        return heavy_xyz(m)

    def rigid_profile(self, step=5):
        X = heavy_xyz(self.m)
        prof = []
        for tau in np.arange(-180, 180, step):
            e, _ = self.rotor_opt(self.set_tau(X, tau))
            prof.append((float(tau), e))
        self.rigid_min = min(prof, key=lambda p: p[1])
        return prof

    def _minimise(self, X, tau):
        m = Chem.Mol(self.m)
        conf = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, X[i].tolist())
        if tau is not None:
            T.SetDihedralDeg(conf, *self.tors, float(tau))
        ff = AllChem.MMFFGetMoleculeForceField(m, self.mp)
        if tau is not None:
            # 1e2 holds the torsion to 0.01 deg; a much stiffer constraint (1e5) stalls the
            # minimiser at its starting point while still reporting convergence (checked 2026-09-28)
            ff.MMFFAddTorsionConstraint(*self.tors, False, float(tau), float(tau), 1.0e2)
        ff.Minimize(maxIts=20000)
        free = AllChem.MMFFGetMoleculeForceField(m, self.mp)
        return free.CalcEnergy(), heavy_xyz(m), T.GetDihedralDeg(m.GetConformer(), *self.tors)

    def relaxed_energy(self, X, tau=None):
        """Lowest MMFF94s energy reachable with C3-C4-S1-N2 held at tau (free if None).  A single
        minimisation lands in whichever polar-H rotor basin it starts in (the sulfonamide NH2 has
        two, 3 kcal/mol apart), so it is started from every rotor combination at 120 deg steps."""
        best = None
        for combo in np.array(np.meshgrid(*[[0, 120, 240]] * len(self.rotors))).T.reshape(
                -1, len(self.rotors)):
            Y = X
            for (a, b, hs), ang in zip(self.rotors, combo):
                Y = rotate(Y, a, b, hs, ang)
            r = self._minimise(Y, tau)
            if best is None or r[0] < best[0]:
                best = r
        return best

    def relaxed_profile(self, step=5):
        X = heavy_xyz(self.m)
        prof = [(float(t), self.relaxed_energy(X, t)[0]) for t in np.arange(-180, 180, step)]
        # global minimum: unconstrained minimisations started every 10 degrees
        starts = [self.relaxed_energy(self.set_tau(X, t))[::2] for t in np.arange(-180, 180, 10)]
        e, tau = min(starts, key=lambda s: s[0])
        self.relaxed_min = (float(tau), float(e))
        return prof

    def pose(self, X):
        tau = T.GetDihedralDeg(self._conf(X), *self.tors)
        e_rigid, Xh = self.rotor_opt(X)
        e_rel, _, tau_rel = self.relaxed_energy(Xh, tau)
        return dict(tau=tau, strain_rigid=e_rigid - self.rigid_min[1],
                    strain_relaxed=e_rel - self.relaxed_min[1], X_hrot=Xh)

    def _conf(self, X):
        m = Chem.Mol(self.m)
        c = m.GetConformer()
        for i in range(m.GetNumAtoms()):
            c.SetAtomPosition(i, X[i].tolist())
        return c


def write_xyz(path, mol, X, comment):
    lines = [str(mol.GetNumAtoms()), comment]
    for a, p in zip(mol.GetAtoms(), X):
        lines.append(f'{a.GetSymbol():2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}')
    path.write_text('\n'.join(lines) + '\n')


# ------------------------------------------------------------------------------------------------
# one job
# ------------------------------------------------------------------------------------------------
def analyse_job(rundir, strain_cache, dft_manifest, diagnostic=False):
    res = json.loads((rundir / 'results.json').read_text())
    target, lig, scoring = res['target'], res['ligand'], res['scoring']
    tdir = DOCK / target
    lrec = json.loads((DOCK / 'ligands' / 'ligands.json').read_text())['species'][lig]
    ha = lrec['heavy_atoms']
    seeds = [s for s in res['seeds'] if s.get('ok')]
    out = dict(target=target, ligand=lig, scoring=scoring, tag=res['tag'],
               geometry=res.get('geometry', 'etkdg'), input_geometry=res.get('input_geometry'),
               exhaustiveness=res['exhaustiveness'], cpu=res['cpu'], n_seeds=len(seeds),
               n_failed=len(res['seeds']) - len(seeds), heavy_atoms=ha, charge=lrec['charge'],
               wall_s_total=rnd(sum(s['wall_s'] for s in seeds), 1),
               receptor=res.get('receptor'), maps=res.get('maps'), diagnostic=bool(diagnostic))
    # one molecule holding the top pose of every seed, in seed order
    poses, all_modes = None, []
    for s in seeds:
        pm = A.load_poses(DOCK / s['out'])
        all_modes.append(pm)
        top = Chem.Mol(pm)
        top.RemoveAllConformers()
        top.AddConformer(Chem.Conformer(pm.GetConformer(pm.GetConformers()[0].GetId())), assignId=True)
        if poses is None:
            poses = top
        else:
            poses.AddConformer(top.GetConformer(), assignId=True)
    cids = [c.GetId() for c in poses.GetConformers()]
    sc = np.array([s['best'] for s in seeds])
    out['score'] = dict(**stats(sc), best_seed=int(seeds[int(np.argmin(sc))]['seed']),
                        per_seed=[rnd(x) for x in sc])
    out['ligand_efficiency'] = dict(
        LE_mean=rnd(-sc.mean() / ha), LE_best=rnd(-sc.min() / ha),
        LE_scale_Reynolds=rnd(reynolds_le_scale(ha)),
        FQ_mean=rnd(-sc.mean() / ha / reynolds_le_scale(ha)))

    # RMSD to crystal ligands: full molecule when the element graphs agree, cores otherwise
    rm = {}
    lig_graph = element_graph(poses)
    for label, ref, info in references(target):
        entry = dict(info)
        if element_graph(ref) == lig_graph:
            entry['full'] = rmsd_block(A.symm_rmsd(ref, poses, cids), Chem.RemoveHs(ref).GetNumAtoms())
        for cname, smarts in CORES:
            v = core_rmsd_values(ref, poses, cids, smarts)
            if v:
                entry.setdefault('cores', {})[cname] = rmsd_block(*v)
        if 'full' in entry:
            # Is a near-native pose present but out-ranked?  The same question analyse.py asks of the
            # validation runs: RMSD of every mode of every seed, and the rank at which the first pose
            # within 1 / 2 A appears (production 2026-09-28).
            per_seed_modes = [A.symm_rmsd(ref, pm, [c.GetId() for c in pm.GetConformers()])
                              for pm in all_modes]
            best_of_modes = [min(v) for v in per_seed_modes]
            entry['modes'] = dict(
                n_modes_per_seed=[len(v) for v in per_seed_modes],
                best_of_modes=rmsd_block(best_of_modes),
                rank_of_first_pose_within=dict())
            for lim, tag in ((1.0, '1A'), (2.0, '2A')):
                ranks = [next((i + 1 for i, x in enumerate(v) if x <= lim), None)
                         for v in per_seed_modes]
                found = [r for r in ranks if r]
                # how much score separates the top pose from that first near-native pose: the
                # quantity that decides which orientation is reported (production 2026-09-28)
                gaps = [seeds[i]['scores'][r - 1] - seeds[i]['scores'][0]
                        for i, r in enumerate(ranks) if r]
                entry['modes']['rank_of_first_pose_within'][tag] = dict(
                    fraction_of_seeds_with_any=rnd(len(found) / len(ranks)),
                    median_rank=rnd(np.median(found), 1) if found else None,
                    best_rank=min(found) if found else None,
                    worst_rank=max(found) if found else None,
                    score_gap_from_top_pose=dict(
                        mean=rnd(np.mean(gaps)), min=rnd(min(gaps)), max=rnd(max(gaps)))
                    if gaps else None)
        if 'full' in entry or 'cores' in entry:
            rm[label] = entry
    out['rmsd_to_crystal'] = rm

    # clusters of the top poses
    D = pairwise_symm_rmsd(poses, cids)
    clusters = []
    for k, (L, mem) in enumerate(leader_clusters(list(sc), D)):
        mem_sc = sc[mem]
        cl = dict(id=k + 1, population=len(mem), fraction=rnd(len(mem) / len(sc)),
                  seeds=[int(seeds[i]['seed']) for i in mem],
                  representative_seed=int(seeds[L]['seed']), score_best=rnd(sc[L]),
                  score_mean=rnd(mem_sc.mean()),
                  max_rmsd_to_representative=rnd(max(D[L, i] for i in mem)))
        # where the representative sits relative to the crystal ligand(s)
        for label, e in rm.items():
            if 'full' in e:
                cl.setdefault('rep_rmsd', {})[f'{label} full'] = e['full']['per_seed'][L]
            for cname, blk in e.get('cores', {}).items():
                cl.setdefault('rep_rmsd', {})[f'{label} {cname}'] = blk['per_seed'][L]
        clusters.append(cl)
    out['clusters'] = dict(method=f'leader clustering of the {len(sc)} seeds\' top poses in score '
                           f'order, symmetry-corrected heavy-atom RMSD <= {CLUSTER_CUTOFF} A',
                           n_clusters=len(clusters), clusters=clusters)
    seed_cluster = {}
    for cl in clusters:
        for s in cl['seeds']:
            seed_cluster[s] = cl['id']

    # site checks (validation + generalised), also on the crystal ligand
    base = A.base_target(target)
    checks = A.SITE_CHECKS.get(base, []) + EXTRA_CHECKS.get(base, [])
    per_seed_checks = [A.run_site_checks(target, tdir, poses, cid, checks=checks) for cid in cids]
    ref0 = references(target)[0][1]
    ref_h = Chem.AddHs(ref0, addCoords=True)
    on_ref = {c['check']: c['distance'] for c in A.run_site_checks(target, tdir, ref_h, 0, checks=checks)}
    sites = {}
    for lab in [c[0] for c in checks]:
        rows = [c for row in per_seed_checks for c in row if c['check'] == lab]
        dists = [c['distance'] for c in rows if c['distance'] is not None]
        sites[lab] = dict(cutoff=next(c[3] for c in checks if c[0] == lab),
                          crystal_ligand_value=on_ref.get(lab),
                          pass_fraction=rnd(np.mean([c['passed'] for c in rows if c['passed'] is not None]))
                          if dists else None,
                          distance_mean=rnd(np.mean(dists), 2) if dists else None,
                          applicable=bool(dists))
    out['site_checks'] = sites
    zn_ok = []
    for row in per_seed_checks:
        cs = {c['check']: c for c in row}
        z, t = cs.get('Zn-N(sulfonamide)'), cs.get('Thr199 OG1 H-bond')
        zn_ok.append(None if z is None else bool(z.get('passed') and (t is None or t.get('passed'))))
    if any(v is not None for v in zn_ok):
        keep = [s for s, ok in zip(sc, zn_ok) if ok]
        out['zinc_bound'] = dict(n_seeds=len(keep), of=len(sc),
                                 score_mean=rnd(np.mean(keep)) if keep else None,
                                 score_sd=rnd(np.std(keep, ddof=1)) if len(keep) > 1 else None,
                                 criterion='top pose Zn-N(sulfonamide) <= 2.3 A and Thr199 OG1...N <= 3.3 A')

    per_seed = []
    for i, s in enumerate(seeds):
        row = dict(seed=int(s['seed']), score=rnd(sc[i]), cluster=seed_cluster[int(s['seed'])],
                   n_modes=s['n_modes'], wall_s=s['wall_s'],
                   site_checks={c['check']: c['distance'] for c in per_seed_checks[i]})
        per_seed.append(row)

    # DFT-docking link for the sulfanilamide species.  Not for the diagnostic runs: their input is a
    # deliberately modified geometry, so the torsion would be measured against a molecule that is not
    # the one docked, and their purpose is the pose, not the strain (production 2026-09-28).
    if lig in DFT_XYZ and not diagnostic:
        if lig not in strain_cache:
            st = Strain(lig)
            prof_r = st.rigid_profile()
            prof_x = st.relaxed_profile()
            strain_cache[lig] = st
            strain_cache[lig + ':profile'] = dict(
                rigid=[(t, rnd(e - st.rigid_min[1])) for t, e in prof_r],
                relaxed=[(t, rnd(e - st.relaxed_min[1])) for t, e in prof_x],
                rigid_min_tau=st.rigid_min[0], relaxed_min_tau=rnd(st.relaxed_min[0], 1))
        st = strain_cache[lig]
        tors = []
        for i, s in enumerate(seeds):
            sdir = DOCK / pathlib.Path(s['out']).parent
            m, resid, _ = rebuild_pose(lig, sdir / 'ligand.pdbqt', sdir / 'out.pdbqt', 1)
            c = m.GetConformer()
            t3 = T.GetDihedralDeg(c, LBL['C3'], LBL['C4'], LBL['S1'], LBL['N2'])
            t5 = T.GetDihedralDeg(c, LBL['C5'], LBL['C4'], LBL['S1'], LBL['N2'])
            X = heavy_xyz(m)
            e = st.pose(X)
            rec = dict(seed=int(s['seed']), cluster=seed_cluster[int(s['seed'])], score=rnd(sc[i]),
                       tau_C3=rnd(t3, 1), tau_C5=rnd(t5, 1),
                       tau_fold=rnd(max(abs(t3), abs(t5)), 1),
                       strain_mmff_rigid=rnd(e['strain_rigid'], 2),
                       strain_mmff_relaxed=rnd(e['strain_relaxed'], 2),
                       rebuild_residual_A=rnd(resid, 4))
            tors.append(rec)
            per_seed[i].update(tau_C3=rec['tau_C3'], tau_fold=rec['tau_fold'],
                               strain_mmff_relaxed=rec['strain_mmff_relaxed'])
            # the representative of every cluster goes to the DFT single-point set (production only)
            if not diagnostic and any(cl['representative_seed'] == s['seed'] for cl in clusters):
                cl = next(cl for cl in clusters if cl['representative_seed'] == s['seed'])
                name = f'{target}__{lig}__{scoring}__cl{cl["id"]}_seed{s["seed"]:02d}'
                pdir = RES / 'poses_for_dft'
                pdir.mkdir(parents=True, exist_ok=True)
                cmt = (f'{name}: docked pose rebuilt on the B3LYP-D3(BJ)/PCM geometry '
                       f'({dft_mol(lig)["xyz"]}); charge {lrec["charge"]} mult 1; score {sc[i]:.3f}')
                write_xyz(pdir / f'{name}.xyz', m, X, cmt + '; polar H as docked')
                write_xyz(pdir / f'{name}__Hrot.xyz', m, e['X_hrot'],
                          cmt + '; polar-H rotors re-optimised (MMFF94s) at the docked heavy atoms')
                dft_manifest.append(dict(
                    file=f'{name}.xyz', file_hrot=f'{name}__Hrot.xyz', target=target, species=lig,
                    scoring=scoring, seed=int(s['seed']), cluster=cl['id'],
                    cluster_population=f'{cl["population"]}/{len(sc)}', score=rnd(sc[i]),
                    charge=lrec['charge'], multiplicity=1, tau_C3=rec['tau_C3'], tau_C5=rec['tau_C5'],
                    tau_fold=rec['tau_fold'], strain_mmff_rigid=rec['strain_mmff_rigid'],
                    strain_mmff_relaxed=rec['strain_mmff_relaxed'],
                    rebuild_residual_A=rec['rebuild_residual_A'],
                    reference_opt_xyz=dft_mol(lig)['xyz'], reference_E_Ha=dft_mol(lig)['E_Ha']))
        tf = np.array([r['tau_fold'] for r in tors])
        out['torsion'] = dict(
            definition='C3-C4-S1-N2 with the crystal atom labels (carried through the DFT xyz); '
                       'tau_fold = max(|C3-C4-S1-N2|, |C5-C4-S1-N2|), the ring-flip-invariant value '
                       '(90 = S-N perpendicular to the ring, 180 = in plane)',
            per_seed=tors, tau_fold=stats(tf),
            strain_mmff_rigid=stats([r['strain_mmff_rigid'] for r in tors]),
            strain_mmff_relaxed=stats([r['strain_mmff_relaxed'] for r in tors]),
            max_rebuild_residual_A=rnd(max(r['rebuild_residual_A'] for r in tors), 4))
    out['per_seed'] = per_seed

    # pose files
    pdir = RES / 'poses' / target
    pdir.mkdir(parents=True, exist_ok=True)
    # Review 2026-09-28: the stem used to be '{lig}__{scoring}' for every run, so the diagnostic runs
    # of section 7b (analysed after the production ones) overwrote the production pose files of the
    # same species and scoring: dock/results/poses/T{3_6RL9,4_3HS4}/sfa_neutral__ad4zn__* held the
    # ETKDG diagnostic poses, 2.0-3.1 A from the production ones, while the JSON numbers (computed
    # from the runs themselves) were correct.  A diagnostic run now keeps its own suffix in the file
    # name, so a diagnostic can never overwrite a production pose.
    suffix = res['tag'][len(f'{lig}__{scoring}'):]
    stem = f'{lig}__{scoring}' + (suffix if diagnostic else '')
    w = Chem.SDWriter(str(pdir / f'{stem}__top_pose_per_seed.sdf'))
    for i, cid in enumerate(cids):
        mm = Chem.Mol(poses)
        mm.RemoveAllConformers()
        mm.AddConformer(Chem.Conformer(poses.GetConformer(cid)), assignId=True)
        mm.SetProp('_Name', f'{target} {lig} {scoring} seed {seeds[i]["seed"]}')
        mm.SetProp('seed', str(seeds[i]['seed']))
        mm.SetProp('score_kcal_mol', f'{sc[i]:.3f}')
        mm.SetProp('cluster', str(seed_cluster[int(seeds[i]['seed'])]))
        w.write(mm)
    w.close()

    # PoseBusters on the production poses themselves (review 2026-09-28).  Validation ran it in
    # 'redock' mode on the jobs that have a reference ligand; production poses -- the ones the paper
    # discusses and the ones handed to the DFT single points -- were not checked at all.  'dock' mode
    # needs no reference, so it covers every species: sane chemistry, bond lengths and angles, ring
    # flatness, internal energy ratio, no clash with protein, cofactors or waters.
    if POSEBUSTERS:
        try:
            from posebusters import PoseBusters
            df = PoseBusters(config='dock', max_workers=1).bust(
                [str(pdir / f'{stem}__top_pose_per_seed.sdf')], None, str(tdir / 'receptor.pdb'))
            cols = [c for c in df.columns if df[c].dtype == bool]
            out['posebusters'] = dict(
                config='dock', n_poses=int(len(df)), n_checks=len(cols),
                pass_fraction=rnd(float(df[cols].all(axis=1).mean()), 3),
                failed_checks={c: int((~df[c]).sum()) for c in cols if (~df[c]).any()})
        except Exception as exc:                             # PoseBusters is strict about inputs
            out['posebusters'] = dict(error=f'{type(exc).__name__}: {exc}')

    best = seeds[int(np.argmin(sc))]
    shutil.copyfile(DOCK / best['out'], pdir / f'{stem}__best_seed{best["seed"]:02d}_all_modes.pdbqt')
    out['pose_files'] = [str((pdir / f'{stem}__top_pose_per_seed.sdf').relative_to(DOCK)),
                         str((pdir / f'{stem}__best_seed{best["seed"]:02d}_all_modes.pdbqt').relative_to(DOCK))]
    return out


# ------------------------------------------------------------------------------------------------
def rank_order(jobs):
    """Species ranked by mean top score per target x scoring, with the reference ligand's score."""
    out = {}
    for t in TARGET_ORDER:
        for s in SCORING_ORDER:
            rows = [j for j in jobs.values() if j['target'] == t and j['scoring'] == s]
            if not rows:
                continue
            rows.sort(key=lambda j: j['score']['mean'])
            refs = [j for j in rows if j['ligand'] in TARGET_INFO[t]['reference_species']]
            ref = refs[0] if refs else None
            ranked = []
            for k, j in enumerate(rows):
                r = dict(rank=k + 1, species=j['ligand'], score_mean=j['score']['mean'],
                         score_sd=j['score']['sd'], best=j['score']['min'],
                         LE=j['ligand_efficiency']['LE_mean'], FQ=j['ligand_efficiency']['FQ_mean'],
                         heavy_atoms=j['heavy_atoms'])
                if ref:
                    r['delta_vs_reference'] = rnd(j['score']['mean'] - ref['score']['mean'])
                if 'zinc_bound' in j:
                    r['zinc_bound_seeds'] = f"{j['zinc_bound']['n_seeds']}/{j['zinc_bound']['of']}"
                    r['zinc_bound_score_mean'] = j['zinc_bound']['score_mean']
                ranked.append(r)
            entry = dict(reference_species=ref['ligand'] if ref else None,
                         reference_score=ref['score']['mean'] if ref else None, ranking=ranked)
            if any('zinc_bound_seeds' in r for r in ranked):
                zb = [r for r in ranked if r.get('zinc_bound_score_mean') is not None]
                zb.sort(key=lambda r: r['zinc_bound_score_mean'])
                entry['ranking_zinc_bound_only'] = [dict(rank=k + 1, species=r['species'],
                                                         score_mean=r['zinc_bound_score_mean'],
                                                         seeds=r['zinc_bound_seeds'])
                                                    for k, r in enumerate(zb)]
            sfa = next((r for r in ranked if r['species'] == 'sfa_neutral'), None)
            if sfa:
                entry['sulfanilamide_neutral_rank'] = f"{sfa['rank']}/{len(ranked)}"
            out[f'{t}/{s}'] = entry
    return out


def le_size_bias(jobs):
    """Pearson r between LE and heavy-atom count, and between score and heavy-atom count, per
    scoring function over every (target, species) job: the size bias of LE measured on this set."""
    out = {}
    for s in SCORING_ORDER:
        rows = [j for j in jobs.values() if j['scoring'] == s]
        if len(rows) < 4:
            continue
        ha = np.array([j['heavy_atoms'] for j in rows], float)
        le = np.array([j['ligand_efficiency']['LE_mean'] for j in rows])
        scr = np.array([j['score']['mean'] for j in rows])
        fq = np.array([j['ligand_efficiency']['FQ_mean'] for j in rows])
        if ha.std() == 0:            # only one ligand size ran with this scoring function
            out[s] = dict(n_jobs=len(rows), note='all jobs have the same heavy-atom count; no '
                                                 'correlation with size can be computed')
            continue
        r = lambda v: rnd(np.corrcoef(ha, v)[0, 1]) if np.std(v) > 0 else None  # noqa: E731
        out[s] = dict(n_jobs=len(rows), heavy_atom_range=[int(ha.min()), int(ha.max())],
                      r_LE_vs_HA=r(le), r_score_vs_HA=r(scr), r_FQ_vs_HA=r(fq))
    return out


def reproduction_checks(jobs):
    """Production jobs that repeat a validation job exactly (same ligand input route, target,
    scoring, seeds, exhaustiveness) must reproduce it bit for bit; DFT-geometry jobs are compared
    with their ETKDG validation counterpart as an input-geometry sensitivity."""
    out = {}
    for key, j in jobs.items():
        vdir = DOCK / j['target'] / 'runs' / f"{j['ligand']}__{j['scoring']}"
        if not (vdir / 'results.json').exists():
            continue
        v = json.loads((vdir / 'results.json').read_text())
        p = json.loads((DOCK / j['target'] / 'runs' / j['tag'] / 'results.json').read_text())
        vs = {s['seed']: s.get('scores') for s in v['seeds'] if s.get('ok')}
        ps = {s['seed']: s.get('scores') for s in p['seeds'] if s.get('ok')}
        common = sorted(set(vs) & set(ps))
        vbest = np.array([vs[k][0] for k in common])
        pbest = np.array([ps[k][0] for k in common])
        rec = dict(validation_run=str(vdir.relative_to(DOCK)), n_common_seeds=len(common),
                   production_geometry=j['geometry'], validation_geometry=v.get('geometry', 'etkdg'),
                   validation_score_mean=rnd(vbest.mean()), production_score_mean=rnd(pbest.mean()),
                   delta_mean=rnd(pbest.mean() - vbest.mean()))
        if j['geometry'] == 'etkdg' and v.get('exhaustiveness') == j['exhaustiveness']:
            rec['identical_all_pose_scores'] = all(vs[k] == ps[k] for k in common)
        out[key] = rec
    return out


def receptor_provenance():
    """sha256 of every receptor these runs used, checked against prep_report_recheck.json, so a
    production number can always be traced to the receptor file that produced it."""
    out = {}
    for t in TARGET_ORDER:
        tdir = DOCK / t
        rec = dict(receptor_pdbqt_sha256=hashlib.sha256((tdir / 'receptor.pdbqt').read_bytes()).hexdigest(),
                   receptor_pdb_sha256=hashlib.sha256((tdir / 'receptor.pdb').read_bytes()).hexdigest())
        rc = tdir / 'prep_report_recheck.json'
        if rc.exists():
            h = json.loads(rc.read_text()).get('hashes', {})
            rec['matches_prep_report_recheck'] = (
                h.get('receptor_pdbqt_sha256') == rec['receptor_pdbqt_sha256']
                and h.get('receptor_pdb_sha256') == rec['receptor_pdb_sha256'])
            rec['prep_report_recheck_says_matches_prep_report'] = h.get('matches_prep_report')
        gpf = tdir / 'ad4zn' / 'receptor_tz.gpf'
        if gpf.exists():
            rec['ad4zn'] = dict(gpf_sha256=hashlib.sha256(gpf.read_bytes()).hexdigest(),
                                n_maps=len(list((tdir / 'ad4zn').glob('receptor_tz.*.map'))),
                                maps_report=json.loads((tdir / 'ad4zn' / 'maps_report.json').read_text()).get(
                                    'tz_to_crystal_sulfonamide_N_A'))
        out[t] = rec
    return out


def main(target=None, do_posebusters=True):
    global POSEBUSTERS
    POSEBUSTERS = do_posebusters
    t0 = time.time()
    RES.mkdir(parents=True, exist_ok=True)
    runs = sorted(DOCK.glob(f'T*/runs/*{SUFFIX}/results.json'))
    runs += sorted(r for suf in DIAG_SUFFIXES for r in DOCK.glob(f'T*/runs/*{suf}/results.json'))
    if target:
        runs = [r for r in runs if r.parents[2].name == target]
    jobs, strain_cache, manifest = {}, {}, []
    pdft = RES / 'poses_for_dft'
    if pdft.exists() and not target:
        shutil.rmtree(pdft)          # regenerated in full; no stale poses
    for rj in runs:
        key = f'{rj.parents[2].name}/{rj.parent.name}'
        diag = rj.parent.name.endswith(DIAG_SUFFIXES)
        try:
            jobs[key] = analyse_job(rj.parent, strain_cache, manifest, diagnostic=diag)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            jobs[key] = dict(error=f'{type(exc).__name__}: {exc}')
            print(f'{key}: ERROR {exc}')
            continue
        j = jobs[key]
        fr = next((e.get('full') for e in j['rmsd_to_crystal'].values() if 'full' in e), None)
        print(f"{key:44s} {j['score']['mean']:7.3f} +- {j['score']['sd']:5.3f}  best {j['score']['min']:7.3f}"
              f"  clusters {[c['population'] for c in j['clusters']['clusters']]}"
              + (f"  RMSD {fr['mean']:.2f} f1 {fr['f_le_1A']:.2f}" if fr else '')
              + (f"  tau_fold {j['torsion']['tau_fold']['mean']:.1f}" if 'torsion' in j else ''),
              flush=True)
    ok = {k: v for k, v in jobs.items() if 'error' not in v}
    prod = {k: v for k, v in ok.items() if not v.get('diagnostic')}
    summary = dict(
        generated=time.strftime('%Y-%m-%d %H:%M %Z'),
        generated_by='python dock/analyse.py --production (dock/production_analysis.py)',
        protocol=dict(suffix=SUFFIX, cluster_cutoff_A=CLUSTER_CUTOFF,
                      notes='scores are Vina/Vinardo/AD4Zn scores of the top pose of each seed, '
                            'kcal/mol by convention; never converted to Ki'),
        targets=TARGET_INFO,
        receptors=receptor_provenance(),
        crystal=dict(torsion=crystal_torsion()),
        csn_angle=csn_angles(),
        dft={k: dict(xyz=v['xyz'], E_Ha=v['E_Ha'], tau_C3=v['tau_C3'], tau_C5=v['tau_C5'],
                     tau_fold=rnd(max(abs(v['tau_C3']), abs(v['tau_C5'])), 2),
                     vs_etkdg_input=dft_vs_etkdg(k))
             for k, v in ((s, dft_mol(s)) for s in DFT_XYZ)},
        mmff_torsion_profiles={k.split(':')[0]: v for k, v in strain_cache.items() if k.endswith(':profile')},
        jobs=jobs,
        rank_order=rank_order(prod),
        ligand_efficiency_size_bias=le_size_bias(prod),
        reproduction_of_validation=reproduction_checks(prod),
        diagnostics=sorted(k for k, v in ok.items() if v.get('diagnostic')),
        poses_for_dft=dict(directory='dock/results/poses_for_dft', n=len(manifest)),
        wall_s=rnd(time.time() - t0, 1))
    if manifest:
        pdft.mkdir(parents=True, exist_ok=True)
        (pdft / 'manifest.json').write_text(json.dumps(dict(
            purpose='single-point (and constrained-optimisation) strain of the docked sulfanilamide '
                    'torsion at B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM(water), to be run on the VM',
            how=['E_strain(rigid) = E_SP(pose.xyz or pose__Hrot.xyz) - E(opt.xyz); the __Hrot file '
                 'has the polar-H rotors re-optimised, because Vina does not score H positions',
                 'better: constrained optimisation with C3-C4-S1-N2 (atoms 3,4,11,8 1-based) frozen '
                 'at tau_C3, minus E(opt): the relaxed torsional strain',
                 'atom order = the DFT opt.xyz order (C1..C6 N1 N2 O1 O2 S1 H...)'],
            poses=manifest), indent=2) + '\n')
    (RES / 'production_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'wrote {RES / "production_summary.json"} ({len(ok)} jobs, {len(jobs) - len(ok)} errors, '
          f'{len(manifest)} poses for DFT) in {time.time() - t0:.0f} s')
    return 0 if len(ok) == len(jobs) else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
