"""
Docking validation analysis: gate G2.

    python dock/analyse.py                      # every job under dock/T*/runs/
    python dock/analyse.py --target T3_6RL9
    python dock/analyse.py --no-posebusters     # skip the slow part

What it computes, and why each piece is there:

*   **Symmetry-corrected RMSD (spyrmsd, Meli 2020, 10.1186/s13321-020-00455-2).**  A plain atom-order
    RMSD is wrong for every one of these ligands: sulfanilamide's two sulfonyl oxygens are
    interchangeable and its ring has a two-fold symmetry, so the same pose can score anything from
    0.3 to 2 A depending on how the atoms happen to be numbered.  spyrmsd solves the graph
    isomorphism, so the number reported is the smallest RMSD over all chemically equivalent
    mappings.  Heavy atoms only, no superposition (`minimize=False`): the pose is compared in the
    receptor frame, which is the question docking is being asked.

*   **Success fractions at 1.0 A and 2.0 A.**  2.0 A is the usual literature criterion (Warren 2006,
    Hawkins 2008) but it is far too generous for an 11-heavy-atom ligand that nearly fills its
    pocket, so the gate is set at 1.0 A for the top pose and both fractions over the ten seeds are
    reported.  A protocol that hits 1 A once in ten tries is not a validated protocol.

*   **PoseBusters (Buttenschoen 2024, 10.1039/D3SC04185A).**  Checks that the pose is a physically
    possible molecule in a physically possible place: sane bond lengths and angles, no internal
    clash, inside the box, not overlapping the protein.  A good RMSD with a failed PoseBusters check
    means the reference frame is right and the chemistry is not.

*   **Site checks.**  These are the interactions the literature says must be present, tested
    geometrically rather than taken on trust: for hCA II the Zn-N(sulfonamide) bond (<= 2.3 A) and
    the Thr199 backbone/OG1 hydrogen bond; for DHPS and Sul2 the contact between the ligand's aniline
    nitrogen and the pterin.  The template paper reported contacts 30-40 A from its ligands; a check
    in code is the only defence.

*   **Rank order against measured affinities.**  Scores for 11-13 heavy-atom ligands are compressed
    into about 1 kcal/mol, so the only defensible use is the ranking.  Acetazolamide is a nanomolar
    hCA II inhibitor and sulfanilamide a micromolar one, so acetazolamide must score better; if it
    does not, the protocol cannot rank and the paper must say so.

Writes dock/analysis/{summary.json, gate_G2.json} and prints a table.
"""
import argparse
import json
import pathlib
import sys

import numpy as np
from rdkit import Chem

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
ANALYSIS = DOCK / 'analysis'

RMSD_GATE = 1.0          # A, top pose, symmetry corrected
RMSD_LOOSE = 2.0         # A, the usual literature criterion

# Measured affinities used only for the ranking test.  Kd/Ki for hCA II from the structure papers of
# the two complexes we redock; they are quoted here to fix the *order*, never to be compared with a
# Vina score.  Sources: Glockner 2020 (10.3390/biom10040509) for sulfanilamide, and the long-standing
# acetazolamide value (Krishnamurthy 2008 review, 10.1021/cr050262p).  Marked so the paper cites them.
AFFINITY_ORDER = {
    'hCA_II': [('azm_neutral', 'acetazolamide, low nM'), ('azm_anion', 'acetazolamide, low nM'),
               ('sfa_neutral', 'sulfanilamide, high nM - low uM'),
               ('sfa_anion', 'sulfanilamide, high nM - low uM')],
}

# Site checks per target: (label, ligand SMARTS for the probe atom, receptor atom selector, cutoff A)
# The receptor selector is (resname, resnum_or_None, atomname).
SITE_CHECKS = {
    'T3_6RL9': [('Zn-N(sulfonamide)', '[NX3,NX2-;$(N[SX4](=O)=O)]', ('ZN', None, 'ZN'), 2.3),
                ('Thr199 OG1 H-bond', '[NX3,NX2-;$(N[SX4](=O)=O)]', ('THR', 199, 'OG1'), 3.3),
                ('Thr199 N H-bond to sulfonyl O', '[OX1;$(O=[SX4])]', ('THR', 199, 'N'), 3.3)],
    'T4_3HS4': [('Zn-N(sulfonamide)', '[NX3,NX2-;$(N[SX4](=O)=O)]', ('ZN', None, 'ZN'), 2.3),
                ('Thr199 OG1 H-bond', '[NX3,NX2-;$(N[SX4](=O)=O)]', ('THR', 199, 'OG1'), 3.3),
                ('Thr199 N H-bond to sulfonyl O', '[OX1;$(O=[SX4])]', ('THR', 199, 'N'), 3.3)],
    # DHPS/Sul2: the aniline (or pterin-facing) nitrogen must reach the pterin.  The pterin is a kept
    # cofactor, so the selector names its residue and we take the closest of its N/O atoms.
    'T1_3TZF': [('aniline N to pterin (HH2)', '[NX3H2;$(Nc1ccccc1)]', ('HH2', None, '*NO'), 3.6)],
    'T2_1AJ0': [('aniline N to pterin (PH2)', '[NX3H2;$(Nc1ccccc1)]', ('PH2', None, '*NO'), 3.6)],
    # T5 (Sul2, dihydropteroate product complex).  Verification 2026-09-28: the check that stood here
    # ("ligand carboxylate within 4.0 A of the pyrophosphate") was wrong twice over -- its SMARTS did
    # not match the neutral crystal component at all, and in the crystal that distance is ~7.5 A, so
    # the check could never pass and reported pass_fraction 0.0 for the best-redocking target in the
    # set.  Replaced by the two contacts the deposited pose actually makes, both measured on
    # T5_7S2K/ref_78H.sdf: the pterin 2-amino N to Asp174 carboxylate (2.60 A) and the pABA
    # carboxylate to the Ser214 backbone N (2.94 A).  Every check in this table is now run against the
    # crystal reference itself as well ('site_checks_on_crystal_reference' in summary.json), which is
    # what would have caught the old one.
    'T5_7S2K': [('pterin 2-amino N to Asp174', '[NX3H2;$(N[#6]=,:[#7])]', ('ASP', 174, '*NO'), 3.2),
                ('pABA carboxylate to Ser214 N',
                 '[OX1,OX2H1,OX2-;$(O~[CX3](~[OX1,OX2H1,OX2-])[#6])]', ('SER', 214, 'N'), 3.4)],
}


# Cross-docking: (target docked into, ligand) -> the target whose co-crystal pose is the reference.
# Comparing an acetazolamide pose against 6RL9's sulfanilamide would be comparing two different
# molecules, so the reference is 3HS4's own acetazolamide pose carried into the 6RL9 frame by
# superposing the two hCA II chains.  The Ca RMSD of that superposition is reported with every number
# it produces, because a cross-docking RMSD is only as good as the alignment behind it.
CROSS_REFERENCE = {
    ('T3_6RL9', 'azm_neutral'): 'T4_3HS4', ('T3_6RL9', 'azm_anion'): 'T4_3HS4',
    ('T4_3HS4', 'sfa_neutral'): 'T3_6RL9', ('T4_3HS4', 'sfa_anion'): 'T3_6RL9',
}


def base_target(name):
    """The real target a receptor-variant directory belongs to.

    Receptor variants are named <TARGET><suffix> (T1_3TZF_seed12345, T2_1AJ0_noSO4,
    T3_6RL9_his64HID, ...), so the base target is simply the longest target key the name starts with.
    Deriving it from TARGETS instead of a hard-coded suffix list means a new variant needs no edit
    here, and there is one definition rather than several (verification 2026-09-28: run_site_checks
    had its own copy of the suffix rule and it already lagged behind this one)."""
    from targets import TARGETS
    hits = [k for k in TARGETS if name == k or name.startswith(k + '_')]
    return max(hits, key=len) if hits else name


def ca_atoms(pdb):
    """{(chain, resnum): coordinates} for the Ca atoms, used to superpose two structures."""
    out = {}
    for ln in pdb.read_text().splitlines():
        if ln.startswith('ATOM') and ln[12:16].strip() == 'CA':
            out[(ln[21], int(ln[22:26]))] = np.array([float(ln[30:38]), float(ln[38:46]),
                                                      float(ln[46:54])])
    return out


def kabsch(P, Q):
    """Rotation + translation that best maps P onto Q (Kabsch, with a reflection guard)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def transform_reference(src_dir, dst_dir, ref, out_sdf):
    """Carry a co-crystal ligand from one receptor's frame into another's.

    Both hCA II entries are the same protein with the same author numbering, so the common Ca atoms
    give a clean superposition; the Ca RMSD after fitting is returned so the reader can judge it."""
    src, dst = ca_atoms(src_dir / 'receptor.pdb'), ca_atoms(dst_dir / 'receptor.pdb')
    keys = sorted(set(src) & set(dst))
    if len(keys) < 50:
        raise SystemExit(f'only {len(keys)} common Ca atoms between {src_dir.name} and {dst_dir.name}')
    P = np.array([src[k] for k in keys])
    Q = np.array([dst[k] for k in keys])
    R, t = kabsch(P, Q)
    fitted = (R @ P.T).T + t
    ca_rmsd = float(np.sqrt(((fitted - Q) ** 2).sum(axis=1).mean()))
    moved = Chem.Mol(ref)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        p = np.array(list(conf.GetAtomPosition(i)))
        q = R @ p + t
        conf.SetAtomPosition(i, q.tolist())
    moved.SetProp('_Name', f'{src_dir.name}_reference_in_{dst_dir.name}_frame')
    w = Chem.SDWriter(str(out_sdf))
    w.write(moved)
    w.close()
    return moved, round(ca_rmsd, 3), len(keys)


def pick_reference(target, ligand, pose_mol):
    """The crystal pose this job should be compared against, and how it was obtained."""
    tdir = DOCK / target
    n_pose = Chem.RemoveHs(pose_mol).GetNumAtoms()
    # The target's own co-crystal ligand is ref_<LIG>.sdf; ref_cross_*.sdf are references this script
    # itself wrote when a cross-dock was analysed, and they must never be mistaken for it.  Verification
    # 2026-09-28: this used to take next(glob('ref_*.sdf')), i.e. whichever file the directory listing
    # happened to return first, so the redock/crossdock decision depended on directory order.
    own_sdfs = sorted(p for p in tdir.glob('ref_*.sdf') if not p.name.startswith('ref_cross_'))
    if not own_sdfs:
        raise SystemExit(f'{target}: no co-crystal reference (ref_<LIG>.sdf) in {tdir}')
    for own_sdf in own_sdfs:
        own = Chem.MolFromMolFile(str(own_sdf))
        if own is not None and Chem.RemoveHs(own).GetNumAtoms() == n_pose:
            return own, own_sdf, dict(kind='redock', reference=own_sdf.name)
    src = CROSS_REFERENCE.get((base_target(target), ligand))
    if src:
        sdir = DOCK / src
        src_sdf = next(sdir.glob('ref_*.sdf'))
        src_ref = Chem.MolFromMolFile(str(src_sdf))
        if Chem.RemoveHs(src_ref).GetNumAtoms() == n_pose:
            out_sdf = tdir / f'ref_cross_{src}_{src_sdf.stem.replace("ref_", "")}.sdf'
            moved, ca_rmsd, n_ca = transform_reference(sdir, tdir, src_ref, out_sdf)
            return moved, out_sdf, dict(kind='crossdock', source=src, reference=out_sdf.name,
                                        superposition_ca_rmsd_A=ca_rmsd, n_ca_atoms=n_ca)
    return None, None, dict(kind='no_reference',
                            note=f'no crystal pose of a {n_pose}-heavy-atom ligand for {target}')


def heavy_graph(mol):
    """Atomic numbers and adjacency matrix of the heavy-atom skeleton, for spyrmsd."""
    m = Chem.RemoveHs(mol)
    n = m.GetNumAtoms()
    adj = np.zeros((n, n), dtype=int)
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i, j] = adj[j, i] = 1
    z = np.array([a.GetAtomicNum() for a in m.GetAtoms()])
    return m, z, adj


def load_poses(out_pdbqt):
    """All poses of one Vina run as a single RDKit molecule with one conformer per pose."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(out_pdbqt), skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    if not mols or mols[0] is None:
        raise SystemExit(f'cannot rebuild poses from {out_pdbqt}')
    return mols[0]


def symm_rmsd(ref_mol, pose_mol, conf_ids):
    from spyrmsd import rmsd as srmsd
    rm, rz, radj = heavy_graph(ref_mol)
    pm, pz, padj = heavy_graph(pose_mol)
    rconf = rm.GetConformer()
    cref = np.array([list(rconf.GetAtomPosition(i)) for i in range(rm.GetNumAtoms())])
    coords = []
    for cid in conf_ids:
        c = pm.GetConformer(cid)
        coords.append(np.array([list(c.GetAtomPosition(i)) for i in range(pm.GetNumAtoms())]))
    if rm.GetNumAtoms() != pm.GetNumAtoms():
        raise SystemExit(f'reference has {rm.GetNumAtoms()} heavy atoms, pose has {pm.GetNumAtoms()}')
    vals = srmsd.symmrmsd(cref, coords, rz, pz, radj, padj, center=False, minimize=False)
    return [float(v) for v in (vals if isinstance(vals, list) else [vals])]


def receptor_atoms(tdir):
    """Heavy atoms of the prepared receptor, as (resname, resnum, atomname) -> coordinates."""
    rows = []
    for ln in (tdir / 'receptor.pdb').read_text().splitlines():
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        el = ln[76:78].strip() or ln[12:16].strip()[0]
        if el == 'H':
            continue
        rows.append((ln[17:20].strip(), int(ln[22:26]), ln[12:16].strip(),
                     np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return rows


def run_site_checks(target, tdir, pose_mol, conf_id, checks=None):
    # base_target() is the one place that knows about receptor-variant suffixes (verification
    # 2026-09-28: this function used to strip them itself, so the two rules could drift apart).
    # `checks` (production 2026-09-28) lets production_analysis.py add checks for species the
    # validation SMARTS were not written for; None keeps the validation behaviour.
    checks = SITE_CHECKS.get(base_target(target), []) if checks is None else checks
    if not checks:
        return []
    rows = receptor_atoms(tdir)
    conf = pose_mol.GetConformer(conf_id)
    out = []
    for label, smarts, (resname, resnum, atomname), cutoff in checks:
        patt = Chem.MolFromSmarts(smarts)
        hits = [i for (i,) in pose_mol.GetSubstructMatches(patt)] if patt.GetNumAtoms() == 1 \
            else [m[0] for m in pose_mol.GetSubstructMatches(patt)]
        lig = [np.array(list(conf.GetAtomPosition(i))) for i in hits]
        cand = [c for (rn, ri, an, c) in rows
                if rn == resname and (resnum is None or ri == resnum)
                and (an == atomname if atomname != '*NO' else an[0] in 'NO')]
        if not lig or not cand:
            out.append(dict(check=label, distance=None, cutoff=cutoff, passed=None,
                            note='probe atom or receptor atom not found'))
            continue
        d = min(float(np.linalg.norm(a - b)) for a in lig for b in cand)
        out.append(dict(check=label, distance=round(d, 2), cutoff=cutoff, passed=bool(d <= cutoff)))
    return out


def posebusters(pose_sdf, ref_sdf, receptor_pdb):
    """PoseBusters in 'redock' mode, with its RMSD test separated from its physics tests.

    The redock configuration ends with an "RMSD <= 2 A" test.  Folding that into a single
    pass/fail would make the PoseBusters column say the same thing as the RMSD column and would fail
    a pose whose chemistry and placement are perfectly sound.  The physics tests -- sane chemistry,
    bond lengths and angles, ring flatness, internal energy ratio, and no overlap with protein,
    cofactors or waters -- are what PoseBusters adds beyond the RMSD, so they are reported on their
    own."""
    from posebusters import PoseBusters
    buster = PoseBusters(config='redock', max_workers=1)
    df = buster.bust([str(pose_sdf)], str(ref_sdf), str(receptor_pdb))
    rec = df.iloc[0].to_dict()
    bools = {k: bool(v) for k, v in rec.items() if isinstance(v, (bool, np.bool_))}
    physics = {k: v for k, v in bools.items() if 'rmsd' not in k.lower()}
    return dict(n_checks=len(bools), n_failed=sum(1 for v in bools.values() if not v),
                failed=[k for k, v in bools.items() if not v], all_passed=all(bools.values()),
                n_physics_checks=len(physics), physics_passed=all(physics.values()),
                physics_failed=[k for k, v in physics.items() if not v])


def analyse_job(rundir, do_posebusters=True):
    res = json.loads((rundir / 'results.json').read_text())
    target, ligand = res['target'], res['ligand']
    tdir = DOCK / target
    out = dict(res)
    out.pop('seeds', None)
    per_seed = []
    ref = ref_sdf = None
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    for s in res['seeds']:
        if not s.get('ok'):
            per_seed.append(dict(seed=s['seed'], ok=False))
            continue
        pose_mol = load_poses(DOCK / s['out'])
        cids = [c.GetId() for c in pose_mol.GetConformers()]
        if ref is None:
            # done once per job: redocking compares against this target's own co-crystal ligand,
            # cross-docking against the partner structure's ligand moved into this frame
            ref, ref_sdf, out['reference'] = pick_reference(target, ligand, pose_mol)
            # A site check is only evidence if the crystal pose passes it.  Running every check on the
            # reference itself turns a mis-specified check into a visible failure instead of a run of
            # zeros that reads like a docking failure (which is exactly what happened to T5 before the
            # 2026-09-28 verification).
            if ref is not None:
                ref_h = Chem.AddHs(ref, addCoords=True)
                out['site_checks_on_crystal_reference'] = run_site_checks(
                    target, tdir, ref_h, ref_h.GetConformer().GetId())
        rmsds = symm_rmsd(ref, pose_mol, cids) if ref is not None else []
        row = dict(seed=s['seed'], ok=True, score_top=s['best'], wall_s=s['wall_s'])
        if rmsds:
            row.update(rmsd_top=round(rmsds[0], 3), rmsd_best_of_modes=round(min(rmsds), 3),
                       rank_of_best_rmsd=int(np.argmin(rmsds)) + 1,
                       score_of_best_rmsd=s['scores'][int(np.argmin(rmsds))],
                       _rmsds=[round(v, 3) for v in rmsds])
        row['site_checks'] = run_site_checks(target, tdir, pose_mol, cids[0])
        if do_posebusters:
            top = ANALYSIS / f'_top_{target}_{ligand}_{res["scoring"]}_seed{s["seed"]:02d}.sdf'
            w = Chem.SDWriter(str(top))
            w.write(pose_mol, confId=cids[0])
            w.close()
            try:
                row['posebusters'] = posebusters(top, ref_sdf, tdir / 'receptor.pdb')
            except Exception as exc:                      # PoseBusters is strict about inputs
                row['posebusters'] = dict(error=f'{type(exc).__name__}: {exc}')
            top.unlink(missing_ok=True)
        per_seed.append(row)

    ok = [r for r in per_seed if r.get('ok')]
    if ok and 'rmsd_top' not in ok[0]:
        # no crystal pose to compare against (e.g. pABA docked into hCA II): report scores and the
        # site checks only, and say so, rather than inventing an RMSD
        sc = np.array([r['score_top'] for r in ok])
        out['stats'] = dict(n_seeds=len(ok), rmsd_available=False,
                            score_mean=round(float(sc.mean()), 3),
                            score_sd=round(float(sc.std(ddof=1)), 3) if len(sc) > 1 else 0.0,
                            wall_s_total=round(float(sum(r['wall_s'] for r in ok)), 1))
        out['per_seed'] = per_seed
        return out
    if ok:
        rt = np.array([r['rmsd_top'] for r in ok])
        sc = np.array([r['score_top'] for r in ok])
        out['stats'] = dict(
            n_seeds=len(ok), rmsd_available=True,
            rmsd_top_mean=round(float(rt.mean()), 3), rmsd_top_sd=round(float(rt.std(ddof=1)), 3)
            if len(rt) > 1 else 0.0,
            rmsd_top_min=round(float(rt.min()), 3), rmsd_top_max=round(float(rt.max()), 3),
            rmsd_top_median=round(float(np.median(rt)), 3),
            success_fraction_1A=round(float((rt <= RMSD_GATE).mean()), 3),
            success_fraction_2A=round(float((rt <= RMSD_LOOSE).mean()), 3),
            best_of_modes_fraction_1A=round(
                float((np.array([r['rmsd_best_of_modes'] for r in ok]) <= RMSD_GATE).mean()), 3),
            score_mean=round(float(sc.mean()), 3),
            score_sd=round(float(sc.std(ddof=1)), 3) if len(sc) > 1 else 0.0,
            wall_s_total=round(float(sum(r['wall_s'] for r in ok)), 1))
        # the informative statistic when scoring, not sampling, is the limitation: how far down the
        # ranked list the first near-native pose sits
        for lim, tag in ((RMSD_GATE, '1A'), (RMSD_LOOSE, '2A')):
            ranks = []
            for r in ok:
                rr = r.get('_rmsds', [])
                hit = next((i + 1 for i, v in enumerate(rr) if v <= lim), None)
                ranks.append(hit)
            found = [x for x in ranks if x]
            out['stats'][f'rank_of_first_pose_within_{tag}'] = dict(
                fraction_of_seeds_with_any=round(len(found) / len(ok), 3),
                # not int(): the median of an even number of seeds can be a half-integer, and
                # int(1.5) = 1 would report the near-native pose as top-ranked more often than it is
                # (verification 2026-09-28).
                median_rank=round(float(np.median(found)), 1) if found else None,
                worst_rank=max(found) if found else None,
                best_rank=min(found) if found else None)
        pb = [r.get('posebusters', {}) for r in ok if isinstance(r.get('posebusters'), dict)]
        pb = [p for p in pb if 'all_passed' in p]
        if pb:
            out['stats']['posebusters_pass_fraction'] = round(
                float(np.mean([p['all_passed'] for p in pb])), 3)
            out['stats']['posebusters_physics_pass_fraction'] = round(
                float(np.mean([p['physics_passed'] for p in pb])), 3)
            out['stats']['posebusters_failed_checks'] = sorted(
                {c for p in pb for c in p['failed']})
            out['stats']['posebusters_physics_failed_checks'] = sorted(
                {c for p in pb for c in p['physics_failed']})
        site = [r['site_checks'] for r in ok if r.get('site_checks')]
        if site:
            labels = [c['check'] for c in site[0]]
            out['stats']['site_checks'] = {
                lab: dict(pass_fraction=round(float(np.mean(
                          [1.0 if c['passed'] else 0.0 for row in site for c in row
                           if c['check'] == lab and c['passed'] is not None])), 3)
                          if any(c['passed'] is not None for row in site for c in row
                                 if c['check'] == lab) else None,
                          distance_mean=round(float(np.mean(
                              [c['distance'] for row in site for c in row
                               if c['check'] == lab and c['distance'] is not None])), 2)
                          if any(c['distance'] is not None for row in site for c in row
                                 if c['check'] == lab) else None)
                for lab in labels}
    out['per_seed'] = per_seed
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default=None)
    ap.add_argument('--no-posebusters', action='store_true')
    ap.add_argument('--production', action='store_true',
                    help='analyse the production runs (tag suffix __prod) instead of the validation '
                         'runs: writes dock/results/production_summary.json (production_analysis.py)')
    a = ap.parse_args()
    if a.production:
        import production_analysis
        return production_analysis.main(target=a.target, do_posebusters=not a.no_posebusters)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    # Production runs (<ligand>__<scoring>__prod) have their own analysis and output directory; they
    # are kept out of the validation summary and gate file so that gate G2 stays the record of the
    # validation runs it was decided on (production 2026-09-28).
    jobs = sorted(j for j in DOCK.glob('T*/runs/*/results.json')
                  if not j.parent.name.endswith('__prod'))
    # With --target, keep the entries for the other targets instead of replacing summary.json with a
    # partial file (verification 2026-09-28: a single-target run used to overwrite the whole summary,
    # so the gate file silently lost every other target).  Entries whose run directory has gone are
    # dropped, so the file cannot accumulate stale jobs either.
    summary = {}
    if a.target and (ANALYSIS / 'summary.json').exists():
        live = {f'{j.parents[2].name}/{j.parent.name}' for j in jobs}
        summary = {k: v for k, v in json.loads((ANALYSIS / 'summary.json').read_text()).items()
                   if k in live}
    for j in jobs:
        target = j.parents[2].name
        if a.target and target != a.target:
            continue
        key = f'{target}/{j.parent.name}'
        try:
            summary[key] = analyse_job(j.parent, do_posebusters=not a.no_posebusters)
        except Exception as exc:
            summary[key] = dict(error=f'{type(exc).__name__}: {exc}')
            print(f'{key}: ERROR {exc}')
            continue
        st = summary[key].get('stats')
        if st:
            print(f'{key:52s} RMSD_top {st["rmsd_top_mean"]:5.2f}+-{st["rmsd_top_sd"]:4.2f} '
                  f'(min {st["rmsd_top_min"]:4.2f}) f1A {st["success_fraction_1A"]:.1f} '
                  f'f2A {st["success_fraction_2A"]:.1f} '
                  f'score {st["score_mean"]:7.2f}+-{st["score_sd"]:4.2f} '
                  f'PBphys {st.get("posebusters_physics_pass_fraction", float("nan")):.1f}')
    (ANALYSIS / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    # --- gate G2 ------------------------------------------------------------------------------
    gate = {}
    for key, s in summary.items():
        st = s.get('stats')
        if not st or s.get('from_crystal'):
            continue
        target = key.split('/')[0]
        zn = st.get('site_checks', {}).get('Zn-N(sulfonamide)', {})
        thr = st.get('site_checks', {}).get('Thr199 OG1 H-bond', {})
        on_ref = s.get('site_checks_on_crystal_reference') or []
        crit = dict(
            site_checks_pass_on_crystal_reference=(
                None if not on_ref else all(c['passed'] for c in on_ref)),
            rmsd_top_mean_le_1A=st['rmsd_top_mean'] <= RMSD_GATE,
            success_fraction_1A=st['success_fraction_1A'],
            success_fraction_2A=st['success_fraction_2A'],
            posebusters_physics=st.get('posebusters_physics_pass_fraction'),
            posebusters_all_incl_rmsd=st.get('posebusters_pass_fraction'),
            rank_of_first_pose_within_1A=st.get('rank_of_first_pose_within_1A'),
            rank_of_first_pose_within_2A=st.get('rank_of_first_pose_within_2A'))
        if target.startswith(('T3', 'T4')):
            crit['zn_n_le_2.3A_fraction'] = zn.get('pass_fraction')
            crit['zn_n_mean_A'] = zn.get('distance_mean')
            crit['thr199_hbond_fraction'] = thr.get('pass_fraction')
        site_ok = (not target.startswith(('T3', 'T4'))
                   or ((zn.get('pass_fraction') or 0) >= 0.5
                       and (thr.get('pass_fraction') or 0) >= 0.5))
        crit['passes_G2'] = bool(crit['rmsd_top_mean_le_1A']
                                 and st['success_fraction_1A'] >= 0.5
                                 and (st.get('posebusters_physics_pass_fraction') or 0) >= 0.8
                                 and site_ok)
        # The conventional literature criterion (2.0 A, Hevener 2009 / Warren 2006), reported beside
        # the project's stricter 1.0 A gate rather than replacing it.  Same structure, same site and
        # PoseBusters requirements, only the distance differs.
        crit['rmsd_top_mean_le_2A'] = st['rmsd_top_mean'] <= RMSD_LOOSE
        crit['passes_2A_criterion'] = bool(crit['rmsd_top_mean_le_2A']
                                           and st['success_fraction_2A'] >= 0.5
                                           and (st.get('posebusters_physics_pass_fraction') or 0) >= 0.8
                                           and site_ok)
        gate[key] = crit

    # ranking test for hCA II: acetazolamide must beat sulfanilamide in the same receptor
    # The test is run per scoring function, because the whole point is to find out which of Vina,
    # Vinardo and AD4Zn can put a nanomolar inhibitor above a micromolar one in the same receptor.
    #
    # The ranking is reported twice, because a score only means something if the pose it came from is
    # a zinc-bound pose.  For hCA II with AD4Zn the top pose of several acetazolamide seeds sits with
    # Zn-N above 5 A (see the per-seed site checks), so a mean over all ten seeds mixes bound and
    # unbound poses.  `all_seeds` is the mean over every seed; `site_valid` is the mean over the seeds
    # whose top pose passes the Zn-N and Thr199-OG1 checks, with the number of seeds it rests on.
    # Verification 2026-09-28: before this change only the unfiltered number was reported.
    def site_valid_scores(job):
        rows = [r for r in job['per_seed'] if r.get('ok')]
        keep = []
        for r in rows:
            cs = {c['check']: c for c in r.get('site_checks', [])}
            zn_c, thr_c = cs.get('Zn-N(sulfonamide)'), cs.get('Thr199 OG1 H-bond')
            if zn_c is None:                       # not a zinc target: no filter to apply
                keep.append(r['score_top'])
            elif zn_c.get('passed') and (thr_c is None or thr_c.get('passed')):
                keep.append(r['score_top'])
        return keep, len(rows)

    ranking = {}
    for tgt in ('T3_6RL9', 'T4_3HS4'):
        for scoring in ('vina', 'vinardo', 'ad4zn'):
            # Only the canonical runs enter the ranking: a tag is <ligand>__<scoring> and nothing
            # else, so every diagnostic variant (__fromxtal, __box14, __repro, __exh64, __s20) is
            # excluded by construction rather than by naming each one (verification 2026-09-28: the
            # old 'fromxtal'/'box' test would have let a __repro run silently replace its original).
            jobs = {k.split('/')[1].split('__')[0]: v
                    for k, v in summary.items()
                    if k.startswith(tgt + '/') and v.get('stats')
                    and k.split('/')[1].split('__')[1:] == [scoring]}
            got = {lig: v['stats']['score_mean'] for lig, v in jobs.items()}
            valid = {}
            for lig, v in jobs.items():
                keep, n = site_valid_scores(v)
                valid[lig] = dict(n_site_valid_seeds=len(keep), n_seeds=n,
                                  score_mean=round(float(np.mean(keep)), 3) if keep else None)
            if 'azm_neutral' in got and 'sfa_neutral' in got:
                va, vs = valid['azm_neutral']['score_mean'], valid['sfa_neutral']['score_mean']
                ranking[f'{tgt}/{scoring}'] = dict(
                    all_seeds=dict(
                        scores=got,
                        azm_better_than_sfa=got['azm_neutral'] < got['sfa_neutral'],
                        gap_kcal=round(got['sfa_neutral'] - got['azm_neutral'], 2)),
                    site_valid=dict(
                        per_ligand=valid,
                        azm_better_than_sfa=(None if va is None or vs is None else va < vs),
                        gap_kcal=(None if va is None or vs is None else round(vs - va, 2)),
                        note='mean over the seeds whose top pose passes the zinc-site checks only; '
                             'n_site_valid_seeds says how many seeds that is'),
                    expectation='acetazolamide (low nM) must score below sulfanilamide '
                                '(high nM - low uM); only the order is meaningful')
    (ANALYSIS / 'gate_G2.json').write_text(json.dumps(
        dict(rmsd_gate_A=RMSD_GATE, criteria=gate, ranking=ranking), indent=2) + '\n')
    print(f'\nwrote {ANALYSIS}/summary.json and gate_G2.json')
    for k, v in gate.items():
        print(f'  G2(1.0 A) {"PASS" if v["passes_G2"] else "fail"}   '
              f'2.0 A criterion {"pass" if v["passes_2A_criterion"] else "fail"}   {k}')
    for k, v in ranking.items():
        a, s = v['all_seeds'], v['site_valid']
        ns = {lig: d['n_site_valid_seeds'] for lig, d in s['per_ligand'].items()}
        print(f'  ranking {k}: all seeds azm<sfa={a["azm_better_than_sfa"]} (gap {a["gap_kcal"]}); '
              f'site-valid poses only azm<sfa={s["azm_better_than_sfa"]} (gap {s["gap_kcal"]}, '
              f'seeds {ns})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
