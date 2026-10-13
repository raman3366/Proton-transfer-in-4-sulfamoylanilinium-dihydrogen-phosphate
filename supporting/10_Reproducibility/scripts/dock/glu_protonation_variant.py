"""Deprotonate one protonated glutamate (GLH) in an already prepared receptor, changing nothing else,
and compare the docking runs of the variant with the runs of the receptor it was made from.

    python dock/glu_protonation_variant.py --target T3_6RL9 --resnum 106
    # writes dock/T3_6RL9_glu106GLU/{receptor.pdbqt,receptor.pdb,box.json,ref_*.sdf,variant.json}
    python dock/glu_protonation_variant.py --target T3_6RL9 --resnum 106 --compare --scratch <dir>
    # writes dock/T3_6RL9_glu106GLU/comparison.json (every run in the variant's runs/ against the
    # run of the same name in dock/T3_6RL9/runs/)
    python dock/glu_protonation_variant.py --target T3_6RL9 --resnum 106 --tables T3_6RL9 T4_3HS4
    # prints the tables of dock/GLU106.md from the two comparison.json files

Why this exists (5 Oct 2026).  prep_receptor.py takes the protonation of Asp/Glu from PROPKA, which was
run on the protein without the zinc.  For hCA II it returned pKa 8.91 for Glu106, so both prepared
hCA II receptors carry Glu106 as the neutral acid (GLH: HE2 on OE2).  In the enzyme Glu106 is expected
to be a carboxylate that accepts the hydrogen bond of Thr199 OG1-H, while Thr199 OG1 accepts the N-H
of the zinc-bound sulfonamide.  The proton lies about 4.8 A from the crystal ligand, so its effect on
the docking has to be measured, not assumed.

Why not simply re-run prep_receptor.py with Glu106 deprotonated.  For the reason given in
his_tautomer_variant.py: the rebuild, the hydrogen placement and the restrained minimisation are not
bit-reproducible, so a freshly prepared receptor would differ from the docked one in many places at
once.  This script edits the existing receptor instead and changes one residue:

    receptor.pdbqt   the acid proton (HE2 or HE1) is removed; every remaining atom of the residue gets
                     the AutoDock type and the partial charge that this same file carries on its
                     deprotonated glutamates.  Meeko 0.8.0 types and charges the polymer from residue
                     templates (GLU, GLH; Gasteiger charges computed once on the capped residue), so
                     all deprotonated mid-chain Glu of a receptor carry one pattern; the script reads
                     that pattern from the file, refuses if the reference residues disagree, and copies
                     the charge/type columns as text, so the variant is what Meeko would have written
                     for template GLU on these coordinates.
    receptor.pdb     the acid proton is removed.

Nothing else is touched.  The serial numbers of the other atoms are kept (one number is left unused),
so that every other line of both files is byte-identical to the source; the script checks this and
records it.  The heavy atoms do not move, so the two C-O bonds of the new carboxylate keep the lengths
the minimisation gave the acid (about 1.21 and 1.33 A); that is part of "one change only" and is
recorded in variant.json.

variant.json also records the geometry that decides whether the edit is sensible: the hydrogen bonds
that polar hydrogens of other residues donate to the two carboxylate oxygens (in hCA II: Thr199 HG1
and the backbone N-H of Arg246), where the removed proton pointed, and its distance to the co-crystal
ligand and to the zinc.

--compare uses the validation analysis itself (analyse.analyse_job: symmetry-corrected RMSD, site
checks, PoseBusters) on the baseline and on the variant run, the gate formulae of analyse.main, and
the orientation definitions of paper/main/data/derive/d14_completion.py (near-native: a mode within
1.0 A of the crystal ligand; ring-rotated: beyond 2.0 A; gap = score of the best-ranked near-native
mode minus the score of the top mode; separation = |best near-native score - best ring-rotated score|).
It writes nothing outside the variant directory and the scratch directory.

Scope: GLU carrying exactly one acid proton.  It refuses a residue that is already a carboxylate.
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
ACID_H = {'HE1': 'OE1', 'HE2': 'OE2'}       # acid proton -> the oxygen that carries it


# --------------------------------------------------------------------------------------------------
# the variant
# --------------------------------------------------------------------------------------------------
def coords(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def is_atom(line):
    return line.startswith(('ATOM', 'HETATM'))


def key(line):
    """(chain, resnum, resname)"""
    return line[21], int(line[22:26]), line[17:20].strip()


def residues(lines, resname):
    """{(chain, resnum): {atom name: line index}} for every residue called `resname`."""
    out = {}
    for i, ln in enumerate(lines):
        if is_atom(ln) and ln[17:20].strip() == resname:
            out.setdefault((ln[21], int(ln[22:26])), {})[ln[12:16].strip()] = i
    return out


def angle(a, b, c):
    v, w = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(v @ w / np.linalg.norm(v) / np.linalg.norm(w), -1, 1))))


def dist(a, b):
    return float(np.linalg.norm(a - b))


def deprotonate_pdbqt(lines, chain, resnum):
    """Return (new lines, record) for the PDBQT: proton removed, charge/type columns of the residue
    replaced by the pattern of the deprotonated Glu residues of the same file."""
    glu = residues(lines, 'GLU')
    if (chain, resnum) not in glu:
        raise SystemExit(f'Glu{resnum} chain {chain} not found')
    idx = glu[(chain, resnum)]
    protons = [h for h in ACID_H if h in idx]
    if not protons:
        raise SystemExit(f'Glu{resnum} carries no HE1/HE2: it is already a carboxylate')
    if len(protons) > 1:
        raise SystemExit(f'Glu{resnum} carries both HE1 and HE2: not a state this script handles')
    hname = protons[0]
    names = sorted(n for n in idx if n != hname)
    # reference pattern: every other Glu of the file with exactly these atoms and no acid proton
    refs = {k: v for k, v in glu.items() if k != (chain, resnum) and sorted(v) == names}
    if not refs:
        raise SystemExit('no deprotonated Glu with the same atom set in this receptor: no pattern to copy')
    pattern = {}
    for n in names:
        cols = {lines[v[n]][66:] for v in refs.values()}
        if len(cols) != 1:
            raise SystemExit(f'the deprotonated Glu residues disagree on atom {n}: {sorted(cols)}; '
                             'the charges are not a residue template, refusing to copy one')
        pattern[n] = cols.pop()
    out, changed = [], {}
    for i, ln in enumerate(lines):
        if i == idx[hname]:
            continue                                   # the acid proton goes
        if i in idx.values():
            n = ln[12:16].strip()
            new = ln[:66] + pattern[n]
            if new != ln:
                changed[n] = dict(before=dict(charge=float(ln[70:76]), type=ln[77:79].strip()),
                                  after=dict(charge=float(new[70:76]), type=new[77:79].strip()))
            ln = new
        out.append(ln)
    q_before = round(sum(float(lines[i][70:76]) for i in idx.values()), 3)
    q_after = round(sum(float(pattern[n][4:10]) for n in names), 3)
    hl = lines[idx[hname]]
    rec = dict(residue=f'{chain}:GLU{resnum}', frm='GLH (neutral acid)', to='GLU (carboxylate)',
               proton_removed=dict(name=hname, on=ACID_H[hname], serial=int(hl[6:11]),
                                   xyz=[round(float(v), 3) for v in coords(hl)],
                                   charge=float(hl[70:76]), type=hl[77:79].strip()),
               atoms_retyped_or_recharged=changed,
               atoms_unchanged=[n for n in names if n not in changed],
               residue_charge_before=q_before, residue_charge_after=q_after,
               pattern_source=dict(n_reference_residues=len(refs),
                                   reference_residues=[f'{c}:GLU{n}' for c, n in sorted(refs)],
                                   all_identical=True,
                                   pattern={n: dict(charge=float(pattern[n][4:10]),
                                                    type=pattern[n][11:13].strip()) for n in names}),
               receptor_total_charge_before=round(sum(float(l[70:76]) for l in lines if is_atom(l)), 3),
               receptor_total_charge_after=round(sum(float(l[70:76]) for l in out if is_atom(l)), 3))
    return out, rec, hname


def deprotonate_pdb(lines, chain, resnum, hname):
    hit = [i for i, ln in enumerate(lines) if is_atom(ln) and key(ln) == (chain, resnum, 'GLU')
           and ln[12:16].strip() == hname]
    if len(hit) != 1:
        raise SystemExit(f'receptor.pdb: expected one {hname} on Glu{resnum}, found {len(hit)}')
    return [ln for i, ln in enumerate(lines) if i != hit[0]], dict(
        residue=f'{chain}:GLU{resnum}', proton_removed=dict(name=hname, serial=int(lines[hit[0]][6:11])))


def diff_check(old, new, chain, resnum):
    """Every line that differs must belong to the residue; exactly one line is gone."""
    olds, news = list(old), list(new)
    removed = [l for l in olds if l not in set(news)]
    added = [l for l in news if l not in set(olds)]
    foreign = [l for l in removed + added if not (is_atom(l) and key(l) == (chain, resnum, 'GLU'))]
    # order of the untouched lines must be the same too
    keep_old = [l for l in olds if l in set(news)]
    keep_new = [l for l in news if l in set(olds)]
    return dict(n_lines_source=len(olds), n_lines_variant=len(news),
                n_lines_removed_or_rewritten=len(removed), n_lines_new=len(added),
                all_differences_inside_the_residue=not foreign,
                other_lines_byte_identical_and_in_order=bool(keep_old == keep_new and not foreign))


def geometry_report(src, chain, resnum, hname):
    """Hydrogen bonds to the carboxylate and the surroundings of the removed proton, measured on the
    source receptor (the heavy atoms and the other hydrogens are the same in the variant)."""
    pq = [l for l in (src / 'receptor.pdbqt').read_text().splitlines() if is_atom(l)]
    at = [dict(chain=l[21], resnum=int(l[22:26]), resname=l[17:20].strip(), name=l[12:16].strip(),
               type=l[77:79].strip(), x=coords(l)) for l in pq]
    me = [a for a in at if (a['chain'], a['resnum'], a['resname']) == (chain, resnum, 'GLU')]
    g = {a['name']: a['x'] for a in me}
    heavy = [a for a in at if a['type'] != 'HD']
    lab = lambda a: f"{a['resname']}{a['resnum']} {a['name']}"  # noqa: E731
    donors = []
    for h in at:
        if h['type'] != 'HD' or h in me:
            continue
        d1, d2 = dist(h['x'], g['OE1']), dist(h['x'], g['OE2'])
        if min(d1, d2) > 3.2:
            continue
        d = min((a for a in heavy if a['resnum'] == h['resnum'] and a['chain'] == h['chain']),
                key=lambda a: dist(a['x'], h['x']))
        o = 'OE1' if d1 <= d2 else 'OE2'
        donors.append(dict(hydrogen=lab(h), donor=lab(d), H_to_OE1_A=round(d1, 2), H_to_OE2_A=round(d2, 2),
                           nearer_oxygen=o, donor_to_that_oxygen_A=round(dist(d['x'], g[o]), 2),
                           angle_D_H_O_deg=round(angle(d['x'], h['x'], g[o]), 1),
                           hydrogen_bond=bool(min(d1, d2) <= 2.5 and angle(d['x'], h['x'], g[o]) >= 120)))
    donors.sort(key=lambda r: min(r['H_to_OE1_A'], r['H_to_OE2_A']))
    hx = g[hname]
    near = sorted(((dist(hx, a['x']), lab(a)) for a in heavy if a not in me))[:5]
    rep = dict(measured_on=f'{src.name}/receptor.pdbqt (polar hydrogens; identical in the variant)',
               polar_H_within_3p2A_of_the_carboxylate=donors,
               removed_proton=dict(name=hname, O_H_A=round(dist(hx, g[ACID_H[hname]]), 3),
                                   nearest_heavy_atoms_of_other_residues=[
                                       dict(atom=n, distance_A=round(d, 2)) for d, n in near]),
               carboxylate_bond_lengths_A=dict(CD_OE1=round(dist(g['CD'], g['OE1']), 3),
                                               CD_OE2=round(dist(g['CD'], g['OE2']), 3),
                                               note='heavy atoms not moved: the lengths are those of '
                                                    'the minimised acid'))
    zn = [a for a in at if a['type'].upper() == 'ZN']
    if zn:
        rep['zinc_distances_A'] = {n: round(dist(zn[0]['x'], g[n]), 2) for n in ('OE1', 'OE2', hname)}
    try:
        from rdkit import Chem
        for sdf in sorted(p for p in src.glob('ref_*.sdf') if not p.name.startswith('ref_cross_')):
            m = Chem.MolFromMolFile(str(sdf))
            c = m.GetConformer()
            L = [(a.GetSymbol(), np.array(list(c.GetAtomPosition(a.GetIdx())))) for a in m.GetAtoms()]
            row = {}
            for n in (hname, 'OE1', 'OE2'):
                d, s = min((dist(g[n], x), s) for s, x in L)
                row[n] = dict(distance_A=round(d, 2), ligand_atom=s)
            rep.setdefault('distance_to_co_crystal_ligand', {})[sdf.name] = row
    except Exception as exc:                                   # the report must not block the variant
        rep['ligand_distance_error'] = f'{type(exc).__name__}: {exc}'
    return rep


def make_variant(target, resnum, chain):
    src = DOCK / target
    dst = DOCK / f'{target}_glu{resnum}GLU'
    dst.mkdir(parents=True, exist_ok=True)
    record = dict(source=target, variant=dst.name, chain=chain, resnum=resnum,
                  note='single-residue deprotonation of an existing receptor (GLH -> GLU); nothing else '
                       'changed; serial numbers of the other atoms kept, so one number is unused')
    pq = (src / 'receptor.pdbqt').read_text().splitlines()
    new, rec, hname = deprotonate_pdbqt(pq, chain, resnum)
    (dst / 'receptor.pdbqt').write_text('\n'.join(new) + '\n')
    rec['file_check'] = diff_check(pq, new, chain, resnum)
    record['receptor.pdbqt'] = rec
    pdb = (src / 'receptor.pdb').read_text().splitlines()
    newp, recp = deprotonate_pdb(pdb, chain, resnum, hname)
    (dst / 'receptor.pdb').write_text('\n'.join(newp) + '\n')
    recp['file_check'] = diff_check(pdb, newp, chain, resnum)
    record['receptor.pdb'] = recp
    for name, r in (('receptor.pdbqt', rec), ('receptor.pdb', recp)):
        assert r['file_check']['other_lines_byte_identical_and_in_order'], name
        assert r['file_check']['n_lines_variant'] == r['file_check']['n_lines_source'] - 1, name
        record[f'{name}_sha256'] = dict(
            source=hashlib.sha256((src / name).read_bytes()).hexdigest(),
            variant=hashlib.sha256((dst / name).read_bytes()).hexdigest())
    record['geometry'] = geometry_report(src, chain, resnum, hname)
    # The box, the Vina config and the target's own co-crystal ligand come along; ref_cross_*.sdf are
    # analyse.py's own cross-docking products and are regenerated in the variant directory.
    own_refs = [p.name for p in src.glob('ref_*.sdf') if not p.name.startswith('ref_cross_')]
    for extra in ['box.json', 'vina_config.txt', *own_refs]:
        if (src / extra).exists():
            shutil.copy2(src / extra, dst / extra)
    (dst / 'variant.json').write_text(json.dumps(record, indent=2) + '\n')
    ch = rec['atoms_retyped_or_recharged']
    print(f"{dst.name}: {rec['residue']} {rec['proton_removed']['name']} removed, residue charge "
          f"{rec['residue_charge_before']:+.3f} -> {rec['residue_charge_after']:+.3f}, receptor "
          f"{rec['receptor_total_charge_before']:+.3f} -> {rec['receptor_total_charge_after']:+.3f}")
    print('  recharged: ' + ', '.join(f"{n} {v['before']['charge']:+.3f} {v['before']['type']} -> "
                                      f"{v['after']['charge']:+.3f} {v['after']['type']}"
                                      for n, v in ch.items()))
    for d in record['geometry']['polar_H_within_3p2A_of_the_carboxylate']:
        print(f"  {d['hydrogen']:14s} H...OE1 {d['H_to_OE1_A']:.2f}  H...OE2 {d['H_to_OE2_A']:.2f}  "
              f"{d['donor']}...{d['nearer_oxygen']} {d['donor_to_that_oxygen_A']:.2f} A, angle "
              f"{d['angle_D_H_O_deg']:.1f} deg, hydrogen bond: {d['hydrogen_bond']}")
    return 0


# --------------------------------------------------------------------------------------------------
# the comparison of the runs
# --------------------------------------------------------------------------------------------------
def gate(o):
    """passes_G2 and passes_2A_criterion as analyse.main() computes them for an hCA II job: mean
    top-pose RMSD within the limit, at least half of the seeds within it, PoseBusters physics pass
    fraction at least 0.8, and the Zn-N and Thr199 OG1 checks passed by at least half of the seeds."""
    import analyse as A
    sc = o.get('site_checks') or {}
    zn, thr = sc.get('Zn-N(sulfonamide)', {}), sc.get('Thr199 OG1 H-bond', {})
    site_ok = (zn.get('pass_fraction') or 0) >= 0.5 and (thr.get('pass_fraction') or 0) >= 0.5
    pb = o.get('posebusters_physics_pass_fraction')
    r = o['rmsd']
    g1 = r['mean'] <= A.RMSD_GATE and r['f_le_1A'] >= 0.5 and site_ok
    g2 = r['mean'] <= A.RMSD_LOOSE and r['f_le_2A'] >= 0.5 and site_ok
    return dict(passes_G2=None if pb is None else bool(g1 and pb >= 0.8),
                passes_2A_criterion=None if pb is None else bool(g2 and pb >= 0.8),
                rmsd_and_site_part_G2=bool(g1), rmsd_and_site_part_2A=bool(g2), site_ok=bool(site_ok))


def orientation(per_seed, seeds):
    """Near-native and ring-rotated modes per seed (definitions of d14_completion.py)."""
    rows = []
    for r, s in zip(per_seed, seeds):
        rr, sc = r['_rmsds'], s['scores']
        near = [i for i, v in enumerate(rr) if v <= 1.0]
        far = [i for i, v in enumerate(rr) if v > 2.0]
        rows.append(dict(seed=r['seed'], rmsd_top=rr[0],
                         top=('near-native' if rr[0] <= 1.0 else 'ring-rotated' if rr[0] > 2.0
                              else 'between'),
                         rank_near=(near[0] + 1) if near else None,
                         gap=round(sc[near[0]] - sc[0], 3) if near else None,
                         sep=round(abs(min(sc[i] for i in near) - min(sc[i] for i in far)), 3)
                         if near and far else None,
                         # signed: negative when the near-native orientation scores better
                         near_minus_rotated=round(min(sc[i] for i in near) - min(sc[i] for i in far), 3)
                         if near and far else None))
    gaps = [x['gap'] for x in rows if x['gap'] is not None]
    seps = [x['sep'] for x in rows if x['sep'] is not None]
    sgn = [x['near_minus_rotated'] for x in rows if x['near_minus_rotated'] is not None]
    n = len(rows)
    return dict(n_seeds=n,
                top_near_native=sum(x['top'] == 'near-native' for x in rows),
                top_ring_rotated=sum(x['top'] == 'ring-rotated' for x in rows),
                top_between=sum(x['top'] == 'between' for x in rows),
                seeds_with_near_native_mode=len(gaps), seeds_with_both=len(seps),
                worst_rank_near=max((x['rank_near'] for x in rows if x['rank_near']), default=None),
                gap_min=min(gaps) if gaps else None, gap_max=max(gaps) if gaps else None,
                sep_min=min(seps) if seps else None, sep_max=max(seps) if seps else None,
                near_minus_rotated_mean=round(float(np.mean(sgn)), 3) if sgn else None,
                near_minus_rotated_min=min(sgn) if sgn else None,
                near_minus_rotated_max=max(sgn) if sgn else None,
                per_seed=rows)


def tau_fold(out_pdbqt):
    """tau_fold of the top pose from the raw coordinates (d14_completion.py): the larger of the two
    |ring C - C4 - S - N| dihedrals.  None when the ligand has no C-SO2-N on a six-membered ring."""
    el, xyz, seen = [], [], False
    for ln in pathlib.Path(out_pdbqt).read_text().splitlines():
        if ln.startswith('MODEL'):
            if seen:
                break
            seen = True
        elif ln.startswith(('ATOM', 'HETATM')):
            t = ln[77:79].strip()
            e = {'A': 'C', 'NA': 'N', 'OA': 'O', 'SA': 'S', 'HD': 'H'}.get(t, t)
            if e != 'H':
                el.append(e)
                xyz.append(coords(ln))
    cut = {('C', 'C'): 1.65, ('C', 'S'): 1.9, ('N', 'S'): 1.8}
    at = list(zip(el, xyz))
    bond = lambda a, b: dist(a[1], b[1]) < cut.get(tuple(sorted((a[0], b[0]))), 0)  # noqa: E731
    s = [a for a in at if a[0] == 'S']
    if len(s) != 1:
        return None
    c4 = [a for a in at if a[0] == 'C' and bond(a, s[0])]
    n = [a for a in at if a[0] == 'N' and bond(a, s[0])]
    if len(c4) != 1 or len(n) != 1:
        return None
    ring = [a for a in at if a[0] == 'C' and a is not c4[0] and bond(a, c4[0])]
    if len(ring) != 2:
        return None

    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1)
        v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
        return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))
    return max(abs(dih(r[1], c4[0][1], s[0][1], n[0][1])) for r in ring)


def clusters(A, seeds, site_rows, rmsd_top):
    """Leader clusters of the seeds' top poses, exactly as production_analysis.py forms them (score
    order, symmetry-corrected heavy-atom RMSD <= 2.0 A); each cluster's representative is its
    best-scoring seed."""
    import production_analysis as P
    from rdkit import Chem
    poses = None
    for s in seeds:
        pm = A.load_poses(DOCK / s['out'])
        top = Chem.Mol(pm)
        top.RemoveAllConformers()
        top.AddConformer(Chem.Conformer(pm.GetConformer(pm.GetConformers()[0].GetId())), assignId=True)
        if poses is None:
            poses = top
        else:
            poses.AddConformer(top.GetConformer(), assignId=True)
    cids = [c.GetId() for c in poses.GetConformers()]
    sc = [s['best'] for s in seeds]
    D = P.pairwise_symm_rmsd(poses, cids)
    out = []
    for k, (L, mem) in enumerate(P.leader_clusters(sc, D)):
        out.append(dict(id=k + 1, population=len(mem), seeds=[int(seeds[i]['seed']) for i in mem],
                        representative_seed=int(seeds[L]['seed']), score_best=sc[L],
                        representative_rmsd=rmsd_top[L] if rmsd_top else None,
                        representative_site_checks=site_rows[L]))
    return out


def analyse_run(A, rundir, do_pb):
    res = json.loads((rundir / 'results.json').read_text())
    job = A.analyse_job(rundir, do_posebusters=do_pb)
    st = job.get('stats', {})
    seeds = [s for s in res['seeds'] if s.get('ok')]
    ps = [r for r in job['per_seed'] if r.get('ok')]
    out = dict(run=str(rundir.relative_to(DOCK)), n_seeds=len(seeds),
               settings=dict(scoring=res['scoring'], exhaustiveness=res['exhaustiveness'], cpu=res['cpu'],
                             box=dict(center=res['box']['center'], size=res['box']['size']),
                             num_modes=res.get('num_modes'), energy_range=res.get('energy_range'),
                             geometry=res.get('geometry', 'etkdg'), vina_version=res.get('vina_version'),
                             maps=res.get('maps')),
               reference=job.get('reference'),
               score=dict(mean=st.get('score_mean'), sd=st.get('score_sd'),
                          per_seed=[s['best'] for s in seeds]),
               wall_s_total=round(sum(s['wall_s'] for s in seeds), 1))
    if st.get('rmsd_available'):
        out['rmsd'] = dict(mean=st['rmsd_top_mean'], sd=st['rmsd_top_sd'], best_of_seeds=st['rmsd_top_min'],
                           f_le_1A=st['success_fraction_1A'], f_le_2A=st['success_fraction_2A'],
                           per_seed=[r['rmsd_top'] for r in ps],
                           best_of_modes_per_seed=[r['rmsd_best_of_modes'] for r in ps],
                           rank_of_first_pose_within_1A=st.get('rank_of_first_pose_within_1A'))
        out['orientation'] = orientation(ps, seeds)
    out['site_checks'] = st.get('site_checks')
    out['site_checks_per_seed'] = [{c['check']: c['distance'] for c in r.get('site_checks', [])} for r in ps]
    out['site_checks_on_crystal_reference'] = job.get('site_checks_on_crystal_reference')
    out['posebusters_physics_pass_fraction'] = st.get('posebusters_physics_pass_fraction')
    out['posebusters_source'] = 'computed in this comparison (redock configuration)' if do_pb else None
    out['posebusters_physics_failed_checks'] = st.get('posebusters_physics_failed_checks')
    out['clusters'] = clusters(A, seeds, out['site_checks_per_seed'],
                               out['rmsd']['per_seed'] if 'rmsd' in out else None)
    taus = [tau_fold(DOCK / s['out']) for s in seeds]
    if all(t is not None for t in taus):
        out['tau_fold_deg'] = dict(mean=round(float(np.mean(taus)), 1), min=round(min(taus), 1),
                                   max=round(max(taus), 1), per_seed=[round(t, 1) for t in taus])
    return out, res


def compare(target, resnum, scratch, do_pb, only=None):
    import analyse as A
    var = f'{target}_glu{resnum}GLU'
    scratch = pathlib.Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    # Nothing may be written into the existing directories: the validation analysis writes a
    # cross-docking reference into the target directory and PoseBusters input into dock/analysis/.
    # Both are redirected: the reference of an existing target goes to the scratch directory (the
    # stored file, which the same code wrote, is what the analysis then reads), and so do the
    # PoseBusters input files.
    A.ANALYSIS = scratch
    orig = A.transform_reference

    def transform_reference(src_dir, dst_dir, ref, out_sdf):
        if pathlib.Path(dst_dir).name != var:
            out_sdf = scratch / f'{pathlib.Path(dst_dir).name}__{pathlib.Path(out_sdf).name}'
        return orig(src_dir, dst_dir, ref, out_sdf)
    A.transform_reference = transform_reference
    out_json = DOCK / var / 'comparison.json'
    result = json.loads(out_json.read_text()) if out_json.exists() else {}
    result.update(source=target, variant=var, analysis='analyse.analyse_job on both runs; gate formulae '
                  'of analyse.main; orientation definitions of d14_completion.py',
                  posebusters='variant: computed here for the runs that have a gate (validation names, '
                              'no suffix); baseline: taken from the stored validation analysis '
                              'dock/analysis/summary.json; runs with a suffix (__prod, diagnostics) carry '
                              'no gate in the paper and are not checked here')
    summary_path = DOCK / 'analysis' / 'summary.json'
    stored = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    runs = result.setdefault('runs', {})
    for rj in sorted((DOCK / var / 'runs').glob('*/results.json')):
        tag = rj.parent.name
        if only and tag not in only:
            continue
        base_dir = DOCK / target / 'runs' / tag
        if not (base_dir / 'results.json').exists():
            runs[tag] = dict(error='no baseline run of this name')
            continue
        has_gate = len(tag.split('__')) == 2          # <ligand>__<scoring>: a validation run
        v, vres = analyse_run(A, rj.parent, do_pb and has_gate)
        b, bres = analyse_run(A, base_dir, False)
        if has_gate:
            st = stored.get(f'{target}/{tag}', {}).get('stats', {})
            b['posebusters_physics_pass_fraction'] = st.get('posebusters_physics_pass_fraction')
            b['posebusters_source'] = 'dock/analysis/summary.json (stored validation analysis)'
            # the stored analysis and this one must agree on the baseline before its PoseBusters
            # figure is borrowed
            if st and 'rmsd' in b:
                b['agrees_with_stored_analysis'] = bool(
                    abs(st['rmsd_top_mean'] - b['rmsd']['mean']) < 5e-4
                    and abs(st['score_mean'] - b['score']['mean']) < 5e-4)
        for o in (b, v):
            if 'rmsd' in o:
                o['gate'] = gate(o)
        # are the inputs the same, and are the outputs?
        same_in, same_out = [], []
        for sv, sb in zip(vres['seeds'], bres['seeds']):
            dv, db = DOCK / pathlib.Path(sv['out']).parent, DOCK / pathlib.Path(sb['out']).parent
            same_in.append((dv / 'ligand.pdbqt').read_bytes() == (db / 'ligand.pdbqt').read_bytes())
            same_out.append((dv / 'out.pdbqt').read_bytes() == (db / 'out.pdbqt').read_bytes())
        sdiff = {k: (b['settings'][k], v['settings'][k]) for k in b['settings']
                 if b['settings'][k] != v['settings'][k] and k != 'maps'}
        d = dict(n_seeds=(b['n_seeds'], v['n_seeds']),
                 settings_that_differ=sdiff,
                 ligand_inputs_identical_all_seeds=all(same_in),
                 out_pdbqt_identical_seeds=f'{sum(same_out)}/{len(same_out)}',
                 all_pose_scores_identical=all(sv.get('scores') == sb.get('scores')
                                               for sv, sb in zip(vres['seeds'], bres['seeds'])),
                 score_mean=round(v['score']['mean'] - b['score']['mean'], 3),
                 score_top_max_abs_per_seed=round(max(abs(x - y) for x, y in zip(
                     v['score']['per_seed'], b['score']['per_seed'])), 3))
        if 'rmsd' in v and 'rmsd' in b:
            d['rmsd_mean'] = round(v['rmsd']['mean'] - b['rmsd']['mean'], 3)
            d['rmsd_top_max_abs_per_seed'] = round(max(abs(x - y) for x, y in zip(
                v['rmsd']['per_seed'], b['rmsd']['per_seed'])), 3)
            d['seeds_whose_top_orientation_changed'] = [
                x['seed'] for x, y in zip(v['orientation']['per_seed'], b['orientation']['per_seed'])
                if x['top'] != y['top']]
            d['gate_changed'] = {k: (b['gate'][k], v['gate'][k]) for k in b['gate']
                                 if b['gate'][k] != v['gate'][k]}
        for lab in (v.get('site_checks') or {}):
            bs, vs = (b.get('site_checks') or {}).get(lab), v['site_checks'][lab]
            if bs and bs.get('distance_mean') is not None and vs.get('distance_mean') is not None:
                d.setdefault('site_distance_mean', {})[lab] = round(vs['distance_mean'] - bs['distance_mean'], 2)
                d.setdefault('site_pass_fraction', {})[lab] = (bs['pass_fraction'], vs['pass_fraction'])
        # a production job of an ETKDG control repeats its validation run bit for bit: say so, so that
        # one variant run stands for both baseline rows
        twin = DOCK / target / 'runs' / (tag + '__prod')
        if bres.get('geometry', 'etkdg') == 'etkdg' and not tag.startswith('sfa_') \
                and (twin / 'results.json').exists():
            tr = json.loads((twin / 'results.json').read_text())
            d['baseline_prod_twin_identical'] = [s.get('scores') for s in tr['seeds']] == \
                [s.get('scores') for s in bres['seeds']]
        runs[tag] = dict(baseline=b, variant=v, delta_variant_minus_baseline=d)
        print(f"{var}/{tag:30s} score {b['score']['mean']:7.3f} -> {v['score']['mean']:7.3f}"
              + (f"  RMSD {b['rmsd']['mean']:.3f} -> {v['rmsd']['mean']:.3f}  f1A "
                 f"{b['rmsd']['f_le_1A']:.2f} -> {v['rmsd']['f_le_1A']:.2f}" if 'rmsd' in v else '')
              + f"  identical poses {d['out_pdbqt_identical_seeds']}", flush=True)
    out_json.write_text(json.dumps(result, indent=2) + '\n')
    print(f'wrote {out_json}')
    return 0


LIG = dict(sfa_neutral='sulfanilamide', sfa_anion='sulfonamidate', azm_neutral='acetazolamide',
           azm_anion='acetazolamide anion', smx_neutral='sulfamethoxazole')
SUF = {'': 'ETKDG, validation', '__prod': 'production', '__csn102': 'DFT geometry, C4-S1-N2 102.0',
       '__csn109': 'ETKDG seed 1, C4-S1-N2 109.5', '__etkdg': 'ETKDG seed 1, every seed'}


def tables(targets, resnum):
    """Markdown tables of dock/GLU106.md from the comparison.json of each variant (stdout)."""
    M = '\u2212'
    def f(x, n=2):
        # the analysis stores three decimals; round half up from that text, so that 3.215 prints as
        # 3.22 and not, through its binary value, as 3.21
        if x is None:
            return '\u2013'
        import decimal
        q = decimal.Decimal(str(x)).quantize(decimal.Decimal(1).scaleb(-n), rounding=decimal.ROUND_HALF_UP)
        return f'{q:.{n}f}'.replace('-', M)
    pm = lambda o, n=2: f"{f(o['mean'], n)} \u00b1 {f(o['sd'], n)}"  # noqa: E731
    yn = lambda g: '\u2013' if g is None else ('pass' if g else 'fail')  # noqa: E731
    data = {t: json.loads((DOCK / f'{t}_glu{resnum}GLU' / 'comparison.json').read_text())['runs'] for t in targets}

    def rows(want):
        for t in targets:
            for tag in sorted(data[t], key=lambda k: (k.split('__')[1] if want == 'ad4zn' else '',
                                                      not k.startswith('sfa'), k.split('__')[2:], k)):
                r = data[t][tag]
                if 'error' in r:
                    continue
                parts = tag.split('__')
                if (parts[1] == 'ad4zn') != (want == 'ad4zn'):
                    continue
                yield t, tag, parts, r

    def describe(t, parts, r):
        suf = '__' + parts[2] if len(parts) > 2 else ''
        geo = SUF[suf]
        if suf == '__prod':
            geo = ('DFT geometry' if r['baseline']['settings']['geometry'] == 'dft' else 'ETKDG') + ', production'
        elif suf == '' and r['delta_variant_minus_baseline'].get('baseline_prod_twin_identical'):
            geo = 'ETKDG, validation = production'
        return (f"{t.split('_')[1]} | {LIG.get(parts[0], parts[0])} | "
                f"{dict(vina='Vina', vinardo='Vinardo', ad4zn='AD4Zn')[parts[1]]} | {geo}, "
                f"{r['baseline']['n_seeds']} seeds")

    for want, title in (('vina', 'Vina and Vinardo'), ('ad4zn', 'AutoDock4Zn')):
        print(f'\n#### {title}\n')
        print('| Receptor | Ligand | Scoring | Input, seeds | | RMSD top (\u00c5) | best of seeds | f \u2264 1 \u00c5 '
              '(\u2264 2 \u00c5) | Score (kcal mol\u207b\u00b9) | Zn\u2013N (\u00c5) / pass | '
              'Thr199 OG1\u00b7\u00b7\u00b7N (\u00c5) / pass | Strict gate / 2.0 \u00c5 |')
        print('|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|:--|')
        for t, tag, parts, r in rows(want):
            for which, lab in (('baseline', 'Glu106 acid (reported)'), ('variant', 'Glu106 carboxylate')):
                o = r[which]
                sc = o.get('site_checks') or {}
                zn, th = sc.get('Zn-N(sulfonamide)', {}), sc.get('Thr199 OG1 H-bond', {})
                site = lambda c: f"{f(c.get('distance_mean'))} / {f(c.get('pass_fraction'))}"  # noqa: E731
                has_gate = len(parts) == 2 and 'gate' in o
                g = f"{yn(o['gate']['passes_G2'])} / {yn(o['gate']['passes_2A_criterion'])}" if has_gate else '\u2013'
                head = describe(t, parts, r) if which == 'baseline' else ' | | | '
                if 'rmsd' in o:
                    nd = 1 if o['n_seeds'] == 10 else 2          # 15 of 20 seeds is 0.75, not 0.8
                    rm = (f"{pm(o['rmsd'])} | {f(o['rmsd']['best_of_seeds'])} | {f(o['rmsd']['f_le_1A'], nd)} "
                          f"({f(o['rmsd']['f_le_2A'], nd)})")
                else:
                    rm = '\u2013 | \u2013 | \u2013'
                print(f"| {head} | {lab} | {rm} | {pm(o['score'], 3)} | {site(zn)} | {site(th)} | {g} |")
            d = r['delta_variant_minus_baseline']
            print(f"| | | | | change | {f(d.get('rmsd_mean'), 3) if 'rmsd_mean' in d else ''} | | | "
                  f"{f(d['score_mean'], 3)} | {f((d.get('site_distance_mean') or {}).get('Zn-N(sulfonamide)'))} | "
                  f"{f((d.get('site_distance_mean') or {}).get('Thr199 OG1 H-bond'))} | "
                  f"largest change of one seed: score {f(d['score_top_max_abs_per_seed'], 3)}"
                  + (f", RMSD {f(d['rmsd_top_max_abs_per_seed'], 3)} \u00c5" if 'rmsd_top_max_abs_per_seed' in d else '')
                  + ' |')
    print('\n#### AutoDock4Zn: the two ring orientations of sulfanilamide\n')
    print('| Receptor | Input, seeds | | Top pose near-native / ring-rotated / between | Seeds holding a '
          'near-native mode (worst rank) | Near-native behind the top mode (kcal mol\u207b\u00b9), min\u2013max | '
          'Separation of the two orientations, min\u2013max (seeds with both) | Best near-native minus best '
          'ring-rotated, mean (min to max) | Best mode of each seed (\u00c5), min\u2013max | \u03c4~fold~ of the top '
          'poses (\u00b0), mean (min\u2013max) |')
    print('|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|')
    for t, tag, parts, r in rows('ad4zn'):
        if not parts[0].startswith('sfa_'):
            continue
        for which, lab in (('baseline', 'acid (reported)'), ('variant', 'carboxylate')):
            o = r[which]
            q = o['orientation']
            rng = lambda a, b, n=3: '\u2013' if a is None else f'{f(a, n)}\u2013{f(b, n)}'  # noqa: E731
            head = describe(t, parts, r).split(' | ') if which == 'baseline' else ['', '', '', '']
            name = f"{head[0]} | {head[1] + ', ' + head[3] if head[0] else ''}"
            bm = o['rmsd']['best_of_modes_per_seed']
            tf = o.get('tau_fold_deg')
            print(f"| {name} | {lab} | {q['top_near_native']} / {q['top_ring_rotated']} / {q['top_between']} | "
                  f"{q['seeds_with_near_native_mode']} of {q['n_seeds']} ({q['worst_rank_near'] or M}) | "
                  f"{rng(q['gap_min'], q['gap_max'])} | {rng(q['sep_min'], q['sep_max'])} ({q['seeds_with_both']}) | "
                  + ('\u2013' if q['near_minus_rotated_mean'] is None else
                     f"{f(q['near_minus_rotated_mean'], 3)} ({f(q['near_minus_rotated_min'], 3)} to "
                     f"{f(q['near_minus_rotated_max'], 3)})")
                  + f" | {f(min(bm))}\u2013{f(max(bm))} | "
                  + (f"{f(tf['mean'], 1)} ({f(tf['min'], 1)}\u2013{f(tf['max'], 1)})" if tf else '\u2013') + ' |')
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', required=True, help='an existing prepared receptor directory name')
    ap.add_argument('--resnum', type=int, required=True)
    ap.add_argument('--chain', default=None, help='default: the chain of the target definition')
    ap.add_argument('--compare', action='store_true',
                    help='compare the runs of the variant with the runs of the source receptor')
    ap.add_argument('--scratch', default=None, help='--compare: directory for temporary files (required)')
    ap.add_argument('--no-posebusters', action='store_true', help='--compare: skip PoseBusters')
    ap.add_argument('--only', nargs='*', default=None, help='--compare: only these run names')
    ap.add_argument('--tables', nargs='*', default=None, metavar='TARGET',
                    help='print the markdown tables of dock/GLU106.md from the comparison.json of these '
                         'targets (e.g. --tables T3_6RL9 T4_3HS4)')
    a = ap.parse_args()
    if a.tables is not None:
        return tables(a.tables or [a.target], a.resnum)
    if a.compare:
        if not a.scratch:
            raise SystemExit('--compare needs --scratch <dir> (nothing is written into existing directories)')
        return compare(a.target, a.resnum, a.scratch, not a.no_posebusters, a.only)
    from targets import TARGETS
    chain = a.chain or TARGETS[a.target.split('_seed')[0]]['chain']
    return make_variant(a.target, a.resnum, chain)


if __name__ == '__main__':
    sys.exit(main())
