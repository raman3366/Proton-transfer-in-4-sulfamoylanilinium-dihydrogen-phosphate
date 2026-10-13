"""
Protein-ligand contacts of the top pose of every docking job, with ProLIF and PLIP.

    python dock/contacts.py                    # all jobs with results
    python dock/contacts.py --target T3_6RL9

Two independent tools on purpose.  The template paper (R4) reported contact residues 29-41 A from the
co-crystal ligands and "hydrogen bonds" of 1.06 and 1.20 A, which are physically impossible; a single
tool run once cannot catch that class of mistake.  ProLIF (Bouysset 2021, 10.1186/s13321-021-00548-6)
works on the RDKit/MDAnalysis representation we already have, PLIP (Adasme 2021,
10.1093/nar/gkab294) re-perceives the complex from a PDB file with its own rules, and the script also
prints the plain donor-acceptor distances so that every reported hydrogen bond can be checked by hand
against a number in the JSON.

Distances are reported donor-to-acceptor (heavy atom to heavy atom, 2.6-3.4 A for a real hydrogen
bond), never H...acceptor, because the H positions are ours and not the crystal's.

Writes dock/analysis/contacts.json.
"""
import argparse
import json
import pathlib
import sys
import warnings

import numpy as np
from rdkit import Chem

warnings.filterwarnings('ignore')
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
ANALYSIS = DOCK / 'analysis'
HBOND_MAX = 3.4          # A, donor heavy atom to acceptor heavy atom
CONTACT_MAX = 4.5        # A, any heavy-atom contact worth listing


def top_pose(out_pdbqt):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(out_pdbqt), skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
    cid = mol.GetConformers()[0].GetId()
    keep = Chem.Mol(mol)
    keep.RemoveAllConformers()
    keep.AddConformer(mol.GetConformer(cid), assignId=True)
    return keep


def geometric_contacts(pose, receptor_pdb):
    """Every receptor residue with a heavy atom within CONTACT_MAX of the ligand, and the polar pairs
    within HBOND_MAX.  Deliberately simple: it is the arithmetic a referee would redo."""
    conf = pose.GetConformer()
    lig = [(a.GetIdx(), a.GetSymbol(), np.array(list(conf.GetAtomPosition(a.GetIdx()))))
           for a in pose.GetAtoms() if a.GetSymbol() != 'H']
    rec = []
    for ln in receptor_pdb.read_text().splitlines():
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[0])
        if el == 'H':
            continue
        rec.append((f'{ln[17:20].strip()}{int(ln[22:26])}', ln[12:16].strip(), el,
                    np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    residues, hbonds = {}, []
    for res, aname, el, rp in rec:
        for li, lel, lp in lig:
            d = float(np.linalg.norm(rp - lp))
            if d <= CONTACT_MAX:
                cur = residues.get(res)
                if cur is None or d < cur['min_dist']:
                    residues[res] = dict(min_dist=round(d, 2), receptor_atom=aname,
                                         ligand_atom=f'{lel}{li}')
            if d <= HBOND_MAX and el in 'NOS' and lel in 'NOS':
                hbonds.append(dict(residue=res, receptor_atom=aname, ligand_atom=f'{lel}{li}',
                                   donor_acceptor_distance=round(d, 2)))
    return (dict(sorted(residues.items(), key=lambda kv: kv[1]['min_dist'])),
            sorted(hbonds, key=lambda h: h['donor_acceptor_distance']))


def prolif_contacts(pose, receptor_pdb):
    import MDAnalysis as mda
    import prolif as plf
    u = mda.Universe(str(receptor_pdb))
    prot = plf.Molecule.from_mda(u.select_atoms('protein or not protein'), force=True)
    lig = plf.Molecule.from_rdkit(Chem.AddHs(pose, addCoords=True))
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig], prot, progress=False)
    df = fp.to_dataframe()
    out = {}
    for col in df.columns:
        if bool(df[col].iloc[0]):
            _, res, inter = col
            out.setdefault(str(res), []).append(inter)
    return out


def plip_contacts(pose, receptor_pdb, tmpdir):
    """PLIP needs a single complex PDB; the ligand goes in as a HETATM residue named LIG."""
    from plip.structure.preparation import PDBComplex
    lines = [ln for ln in receptor_pdb.read_text().splitlines()
             if ln.startswith(('ATOM', 'HETATM'))]
    mol = Chem.AddHs(pose, addCoords=True)
    conf = mol.GetConformer()
    n = len(lines)
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        el = atom.GetSymbol()
        lines.append(f'HETATM{n + i + 1:5d} {el + str(i + 1):<4s} LIG Z 900    '
                     f'{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2s}')
    complex_pdb = tmpdir / 'complex.pdb'
    complex_pdb.write_text('\n'.join(lines) + '\nEND\n')
    pc = PDBComplex()
    pc.load_pdb(str(complex_pdb))
    pc.analyze()
    out = {}
    for key, site in pc.interaction_sets.items():
        if 'LIG' not in key:
            continue
        rec = dict(
            hbonds=[dict(residue=f'{h.restype}{h.resnr}', distance=round(h.distance_ad, 2),
                         protein_is_donor=bool(h.protisdon))
                    for h in list(site.hbonds_pdon) + list(site.hbonds_ldon)],
            hydrophobic=[f'{h.restype}{h.resnr}' for h in site.hydrophobic_contacts],
            pi_stacking=[f'{h.restype}{h.resnr}' for h in site.pistacking],
            salt_bridges=[f'{s.restype}{s.resnr}' for s in
                          list(site.saltbridge_lneg) + list(site.saltbridge_pneg)],
            water_bridges=len(site.water_bridges),
            metal_complexes=[dict(residue=f'{m.restype}{m.resnr}', metal=m.metal.type,
                                  target=m.target.atom.type, distance=round(m.distance, 2))
                             for m in site.metal_complexes])
        out[key] = rec
    complex_pdb.unlink(missing_ok=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default=None)
    ap.add_argument('--seed', type=int, default=1, help='which seed\'s top pose to profile')
    ap.add_argument('--production', action='store_true',
                    help='profile the representative pose of every cluster of every production job '
                         '(needs analyse.py --production first); writes '
                         'dock/results/production_contacts.json and merges key contacts into '
                         'production_summary.json')
    a = ap.parse_args()
    if a.production:
        return production_main()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    tmpdir = ANALYSIS / '_tmp'
    tmpdir.mkdir(exist_ok=True)
    # production runs (__prod) are profiled by --production into dock/results/, not here
    jobs = sorted(j for j in DOCK.glob('T*/runs/*/results.json')
                  if not j.parent.name.endswith('__prod'))
    # As in analyse.py: with --target, keep the other targets' entries instead of replacing
    # contacts.json with a single-target file, and drop entries whose run directory has gone
    # (verification 2026-09-28).
    out = {}
    if a.target and (ANALYSIS / 'contacts.json').exists():
        live = {f'{j.parents[2].name}/{j.parent.name}' for j in jobs}
        out = {k: v for k, v in json.loads((ANALYSIS / 'contacts.json').read_text()).items()
               if k in live}
    for res_json in jobs:
        target = res_json.parents[2].name
        if a.target and target != a.target:
            continue
        res = json.loads(res_json.read_text())
        seeds = [s for s in res['seeds'] if s.get('ok')]
        if not seeds:
            continue
        pick = next((s for s in seeds if s['seed'] == a.seed), seeds[0])
        key = f'{target}/{res_json.parent.name}'
        pose = top_pose(DOCK / pick['out'])
        rec = dict(seed=pick['seed'], score=pick['best'], ligand=res['ligand'],
                   scoring=res['scoring'])
        receptor_pdb = DOCK / target / 'receptor.pdb'
        rec['residues_within_4.5A'], rec['polar_pairs_within_3.4A'] = \
            geometric_contacts(pose, receptor_pdb)
        for name, fn in (('prolif', prolif_contacts), ('plip', plip_contacts)):
            try:
                rec[name] = fn(pose, receptor_pdb) if name == 'prolif' \
                    else fn(pose, receptor_pdb, tmpdir)
            except Exception as exc:
                rec[name] = dict(error=f'{type(exc).__name__}: {exc}')
        out[key] = rec
        near = list(rec['residues_within_4.5A'])[:6]
        print(f'{key:52s} score {rec["score"]:7.2f}  nearest: {", ".join(near)}')
        print(f'{"":52s} polar pairs <=3.4 A: '
              f'{[(h["residue"], h["donor_acceptor_distance"]) for h in rec["polar_pairs_within_3.4A"][:5]]}')
    (ANALYSIS / 'contacts.json').write_text(json.dumps(out, indent=2) + '\n')
    print(f'\nwrote {ANALYSIS}/contacts.json')
    return 0


# ------------------------------------------------------------------------------------------------
# Production mode (2026-09-28).  One profile per cluster representative, with every hydrogen bond
# given as donor...acceptor heavy-atom distance and the ligand atom named by its chemical role, the
# zinc coordination measured directly, and pi contacts from both tools.
# ------------------------------------------------------------------------------------------------
RES = DOCK / 'results'
ROLE_SMARTS = [   # first match wins
    ('N(sulfonamide)', '[#7]-[#16](~[#8])~[#8]'),
    ('O(sulfonyl)', '[#8]~[#16X4](~[#8])'),
    ('O(carboxylate)', '[#8]~[#6X3](~[#8])-[#6]'),
    ('N(amide)', '[#7]-[#6]=[#8]'),
    ('O(amide)', '[#8]=[#6]-[#7]'),
    ('N(aniline)', '[#7;!R;$([#7]-c)]'),
    ('N(ring)', '[#7;R]'),
    ('O(ring)', '[#8;R]'),
    ('S(ring)', '[#16;R]'),
    ('C(ring)', '[#6;R]'),
]
_PROT = {}


def atom_roles(mol):
    roles = {}
    for label, sm in ROLE_SMARTS:
        # uniquify=False: a sulfonyl or carboxylate matches once per oxygen, not once per atom set
        for match in mol.GetSubstructMatches(Chem.MolFromSmarts(sm), uniquify=False):
            i = match[0]
            roles.setdefault(i, f'{label}#{i}')
    return {a.GetIdx(): roles.get(a.GetIdx(), f'{a.GetSymbol()}#{a.GetIdx()}') for a in mol.GetAtoms()}


def receptor_rows(receptor_pdb):
    rows = []
    for ln in receptor_pdb.read_text().splitlines():
        if ln.startswith(('ATOM', 'HETATM')):
            rows.append(dict(res=f'{ln[17:20].strip()}{int(ln[22:26])}', atom=ln[12:16].strip(),
                             el=(ln[76:78].strip() or ln[12:16].strip()[0]),
                             xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return rows


def profile_pose(pose, target, receptor_pdb, tmpdir):
    """Everything the report needs about one pose, from three independent routes."""
    mol = Chem.AddHs(pose, addCoords=True)
    roles = atom_roles(mol)
    conf = mol.GetConformer()
    L = {a.GetIdx(): np.array(list(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms()}
    rows = receptor_rows(receptor_pdb)
    heavy = [r for r in rows if r['el'] != 'H']
    out = {}
    # 1. plain geometry: residues within 4.0 A, every N/O...N/O pair within 3.5 A
    near, pairs = {}, []
    for r in heavy:
        for i, p in L.items():
            if mol.GetAtomWithIdx(i).GetSymbol() == 'H':
                continue
            d = float(np.linalg.norm(r['xyz'] - p))
            if d <= 4.0 and (r['res'] not in near or d < near[r['res']]):
                near[r['res']] = round(d, 2)
            if d <= 3.5 and r['el'] in ('N', 'O') and mol.GetAtomWithIdx(i).GetSymbol() in 'NO':
                pairs.append(dict(residue=r['res'], receptor_atom=r['atom'], ligand_atom=roles[i],
                                  distance=round(d, 2)))
    out['residues_within_4A'] = dict(sorted(near.items(), key=lambda kv: kv[1]))
    out['polar_pairs_NO_within_3.5A'] = sorted(pairs, key=lambda x: x['distance'])
    # 2. zinc, measured directly
    zn = [r for r in heavy if r['res'].startswith('ZN')]
    if zn:
        zc = zn[0]['xyz']
        lig = [(float(np.linalg.norm(p - zc)), roles[i]) for i, p in L.items()
               if mol.GetAtomWithIdx(i).GetSymbol() != 'H']
        lig.sort()
        out['zinc'] = dict(nearest_ligand_atom=lig[0][1], distance=round(lig[0][0], 2),
                           ligand_N_O_within_3A=[dict(atom=a, distance=round(d, 2)) for d, a in lig
                                                 if d <= 3.0 and a[0] in 'NO'])
    # 3. PLIP
    try:
        from plip.structure.preparation import PDBComplex
        lines = [ln for ln in receptor_pdb.read_text().splitlines() if ln.startswith(('ATOM', 'HETATM'))]
        n = len(lines)
        for i, atom in enumerate(mol.GetAtoms()):
            p = conf.GetAtomPosition(i)
            el = atom.GetSymbol()
            lines.append(f'HETATM{n + i + 1:5d} {el + str(i + 1):<4s} LIG Z 900    '
                         f'{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2s}')
        cpdb = tmpdir / f'complex_{target}.pdb'
        cpdb.write_text('\n'.join(lines) + '\nEND\n')
        pc = PDBComplex()
        pc.load_pdb(str(cpdb))
        pc.analyze()

        # PLIP reports the PDB atom serial, so resolve it through the serials we wrote rather than
        # through line order (they agree for these receptors, but only by luck)
        by_serial = {int(ln[6:11]): ln for ln in lines}

        def lig_or_rec(idx):
            if idx > n:
                return 'ligand', roles.get(idx - n - 1, f'#{idx - n - 1}')
            ln = by_serial.get(idx)
            if ln is None:
                return 'receptor', f'serial{idx}'
            return 'receptor', f'{ln[17:20].strip()}{int(ln[22:26])}:{ln[12:16].strip()}'

        site = next((v for k, v in pc.interaction_sets.items() if 'LIG' in k), None)
        if site is not None:
            hb = []
            for h in list(site.hbonds_pdon) + list(site.hbonds_ldon):
                dside, dname = lig_or_rec(h.d_orig_idx)
                aside, aname = lig_or_rec(h.a_orig_idx)
                hb.append(dict(residue=f'{h.restype}{h.resnr}', donor=dname, acceptor=aname,
                               ligand_is_donor=dside == 'ligand', D_A=round(h.distance_ad, 2),
                               H_A=round(h.distance_ah, 2), angle_DHA=round(h.angle, 1),
                               sidechain=bool(h.sidechain)))
            out['plip'] = dict(
                hbonds=sorted(hb, key=lambda x: x['D_A']),
                metal=[dict(residue=f'{m.restype}{m.resnr}', metal=m.metal_type,
                            ligand_atom=lig_or_rec(m.target_orig_idx)[1], distance=round(m.distance, 2),
                            geometry=m.geometry, coordination_number=m.coordination_num)
                       for m in site.metal_complexes],
                pi_stacking=[dict(residue=f'{p.restype}{p.resnr}', type=p.type,
                                  distance=round(p.distance, 2), angle=round(p.angle, 1),
                                  offset=round(p.offset, 2)) for p in site.pistacking],
                pi_cation=[dict(residue=f'{p.restype}{p.resnr}', distance=round(p.distance, 2),
                                protein_charged=bool(p.protcharged))
                           for p in list(site.pication_laro) + list(site.pication_paro)],
                salt_bridges=[dict(residue=f'{sb.restype}{sb.resnr}', distance=round(sb.distance, 2))
                              for sb in list(site.saltbridge_lneg) + list(site.saltbridge_pneg)],
                hydrophobic=sorted({f'{h.restype}{h.resnr}' for h in site.hydrophobic_contacts}),
                water_bridges=len(site.water_bridges))
        cpdb.unlink(missing_ok=True)
    except Exception as exc:
        out['plip'] = dict(error=f'{type(exc).__name__}: {exc}')
    # 4. ProLIF, with the distance of every interaction
    try:
        import MDAnalysis as mda
        import prolif as plf
        key = str(receptor_pdb)
        if key not in _PROT:
            u = mda.Universe(key)
            _PROT[key] = plf.Molecule.from_mda(u.select_atoms('protein or not protein'), force=True)
        fp = plf.Fingerprint()
        fp.run_from_iterable([plf.Molecule.from_rdkit(mol)], _PROT[key], progress=False)
        pl = {}
        for (lres, pres), inter in fp.ifp[0].items():
            for name, metas in inter.items():
                if name == 'VdWContact':
                    continue
                md = metas[0]
                li = md['parent_indices']['ligand']
                pl.setdefault(str(pres), []).append(dict(
                    type=name, distance=round(float(md.get('distance', float('nan'))), 2),
                    ligand_atoms=[roles.get(int(i), str(i)) for i in li][:3] if len(li) <= 3
                    else 'ring'))
        out['prolif'] = pl
    except Exception as exc:
        out['prolif'] = dict(error=f'{type(exc).__name__}: {exc}')
    return out


def production_main():
    import time
    t0 = time.time()
    summ_path = RES / 'production_summary.json'
    summ = json.loads(summ_path.read_text())
    tmpdir = RES / '_tmp'
    tmpdir.mkdir(parents=True, exist_ok=True)
    allc = {}
    for key, job in summ['jobs'].items():
        if 'error' in job:
            continue
        target = job['target']
        receptor_pdb = DOCK / target / 'receptor.pdb'
        res = json.loads((DOCK / target / 'runs' / job['tag'] / 'results.json').read_text())
        by_seed = {s['seed']: s for s in res['seeds'] if s.get('ok')}
        allc[key] = {}
        for cl in job['clusters']['clusters']:
            s = by_seed[cl['representative_seed']]
            pose = top_pose(DOCK / s['out'])
            prof = profile_pose(pose, target, receptor_pdb, tmpdir)
            prof.update(cluster=cl['id'], population=cl['population'], seed=s['seed'], score=s['best'])
            allc[key][str(cl['id'])] = prof
            cl['contacts'] = key_contacts(prof)
        pl = allc[key]
        first = pl[min(pl, key=int)]
        hb = [f"{h['residue']} {h['donor'].split(':')[-1]}->{h['acceptor'].split(':')[-1]} {h['D_A']}"
              for h in first.get('plip', {}).get('hbonds', [])][:4]
        print(f'{key:44s} cl1: {hb}' + (f"  Zn {first['zinc']['nearest_ligand_atom']} "
                                        f"{first['zinc']['distance']}" if 'zinc' in first else ''),
              flush=True)
    (RES / 'production_contacts.json').write_text(json.dumps(allc, indent=2) + '\n')
    summ['contacts'] = dict(file='dock/results/production_contacts.json',
                            generated_by='python dock/contacts.py --production',
                            what='one profile per cluster representative (top pose of the best-scoring '
                                 'seed of the cluster): PLIP hydrogen bonds (donor...acceptor D_A), '
                                 'metal, pi-stacking, pi-cation, salt bridges; ProLIF interactions '
                                 'with distances; plain-geometry N/O pairs <= 3.5 A and residues '
                                 '<= 4.0 A; zinc distances measured directly')
    summ_path.write_text(json.dumps(summ, indent=2) + '\n')
    try:
        tmpdir.rmdir()
    except OSError:
        pass
    print(f'wrote {RES / "production_contacts.json"} and merged key contacts into '
          f'{summ_path.name} in {time.time() - t0:.0f} s')
    return 0


def key_contacts(prof):
    """The compact form stored in production_summary.json."""
    p = prof.get('plip', {})
    return dict(
        hbonds=[dict(residue=h['residue'], donor=h['donor'], acceptor=h['acceptor'], D_A=h['D_A'],
                     angle_DHA=h['angle_DHA']) for h in p.get('hbonds', [])],
        # PLIP applies an angle criterion, so a short N/O...N/O pair it rejects (the Thr199 OG1
        # contact of the zinc-bound sulfonamide, for one) would otherwise vanish from the report;
        # the plain distances are listed beside it (production 2026-09-28)
        polar_pairs=prof.get('polar_pairs_NO_within_3.5A', [])[:4],
        zinc=prof.get('zinc'),
        metal_plip=p.get('metal', []),
        pi=[dict(residue=x['residue'], type=f"stacking-{x['type']}", distance=x['distance'])
            for x in p.get('pi_stacking', [])]
        + [dict(residue=x['residue'], type='pi-cation', distance=x['distance'])
           for x in p.get('pi_cation', [])],
        prolif_pi=[dict(residue=r, **i) for r, lst in prof.get('prolif', {}).items()
                   if isinstance(lst, list) for i in lst
                   if i['type'] in ('PiStacking', 'FaceToFace', 'EdgeToFace', 'CationPi', 'PiCation')],
        salt_bridges=p.get('salt_bridges', []),
        hydrophobic=p.get('hydrophobic', []),
        residues_within_4A=list(prof.get('residues_within_4A', {}))[:12])


if __name__ == '__main__':
    sys.exit(main())
