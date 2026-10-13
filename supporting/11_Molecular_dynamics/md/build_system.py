"""Build the MD input structures of the hCA II - sulfanilamide complex (stage S3 of md/PLAN.md, part 1).

    /root/cqc-sfa-tools/dock/bin/python md/build_system.py

From the deposited model (dock/pdb/6RL9.cif.gz, chain A) and the decisions written in this file it
produces, in md/system/:

    protein_md.pdb      heavy atoms of the protein, one conformer per residue, Amber residue names
                        (HID/HIE, GLU, ASP, LYS, CYS ...), ACE cap on the first observed residue;
                        tleap adds the hydrogens and the few side-chain atoms the entry does not model
    zinc.pdb            the zinc ion
    waters_all.pdb      the crystal water kept before any ligand is considered
    waters_X.pdb, waters_N.pdb, waters_R.pdb   the water of each start (water on the ligand removed)
    ligand_X.sdf/.pdb, ligand_N..., ligand_R...   the three start poses of the sulfonamidate anion,
                        all hydrogens, heavy-atom names of SAN in 6RL9
    starts.csv          the numbers of the three starts
    residue_map.csv     deposited residue number -> residue name here -> index after tleap
    decisions.json      every choice with its reason, the sources and their checksums

The evidence behind the decisions is produced by md/system_audit.py (md/system/audit/); the test with
tleap is md/check_tleap.py.  There are no random numbers in this script: running it again gives the
same files, byte for byte.  The parameters of the ligand and of the zinc site are made elsewhere.
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import system_common as C  # noqa: E402

PH = 7.4

# ----------------------------------------------------------------------------------------------------
# DECISIONS.  Each entry: state, the state in the docking receptor, and the reason.  The numbers quoted
# in the reasons are in md/system/audit/ (histidines.json, glu106.json, propka.json, deposited.json).
# References read in full: Fisher, Raum, Weininger 2025, ChemBioChem 26, e202400930
# (10.1002/cbic.202400930; NMR and neutron proton occupancies of the 12 histidines); Raum, Fisher,
# Weininger 2023, Cell. Mol. Life Sci. 80, 286 (10.1007/s00018-023-04936-z; His64); Fisher et al. 2011,
# Biochemistry 50, 9421 (10.1021/bi201487b; neutron structure at pH 7.8, author manuscript); Fisher et
# al. 2012, J. Am. Chem. Soc. 134, 14726 (10.1021/ja3068098; neutron structure with acetazolamide).
# ----------------------------------------------------------------------------------------------------
HIS = {
    3: ('HIE', 'HIE', 'surface, B 23 A2; NE2...HOH486 2.86 A, ND1 without partner: not decided by the structure; '
                      'NMR gives 71-85 % of the neutral form as N(eps)-H and pKa 6.1 (Fisher 2025); the entry carries HE2'),
    4: ('HID', 'HIE', 'ND1 is 2.70 A from the carbonyl O of Trp5, which can only accept, and NE2 has no partner; the '
                      'entry carries HD1. NMR in solution gives 60-72 % N(eps)-H (pKa 6.2), so the tautomer is mixed; '
                      'surface residue, B 29 A2, 15 A from the zinc'),
    10: ('HIE', 'HIE', 'conformer A: no hydrogen-bond partner within 3.4 A; NMR 79-93 % N(eps)-H, pKa 5.6; the entry carries HE2'),
    15: ('HIE', 'HIE', 'conformer A: ND1 accepts from Lys18 NZ (2.94 A), NE2-H donates to the carbonyl O of Lys9 '
                       '(2.74 A); NMR 76-90 % N(eps)-H and the neutron structures agree'),
    17: ('HIE', 'HIE', 'NE2-H donates to the carbonyl O of Glu187 (2.86 A); ND1...HOH525 2.85 A; NMR finds both tautomers (30-70 %)'),
    36: ('HIE', 'HIE', 'conformer A: ND1...HOH412 2.59 A, NE2 only a water of a neighbouring molecule; NMR 57-69 % '
                       'N(eps)-H. Its pKa is 6.94 (NMR), so about a quarter is protonated at pH 7.4: the neutral '
                       'majority state is taken; 25 A from the zinc'),
    64: ('HIE', 'HIE', 'ASSUMPTION. Conformer A (out, occupancy 0.55; in: B 0.39 + C 0.05) as in the docking. NMR: in and '
                       'out about 1:1 at every pH, pKa 6.25 (in) and 7.60 (out), 6.27 and 7.43 with acetazolamide '
                       'bound, the two neutral tautomers about 1:1 (Raum 2023). At pH 7.4 the neutral form is the '
                       'majority (about two thirds if in and out are 1:1) but not the only one; ND1...HOH464 2.73 A '
                       'decides nothing. HIE is '
                       'kept as in the docking receptor and in the entry'),
    94: ('HID', 'HID', 'zinc ligand through NE2 (2.01 A): the hydrogen is on ND1 and donates to Gln92 OE1 (2.79 A); '
                       'NMR and neutron: 100 % N(delta)-H'),
    96: ('HID', 'HID', 'zinc ligand through NE2 (2.02 A): the hydrogen is on ND1 and donates to the carbonyl O of '
                       'Asn244 (2.80 A); NMR and neutron: 100 % N(delta)-H'),
    107: ('HID', 'HIE', 'buried; ND1-H donates to Glu117 OE2 (2.76 A), NE2 accepts from Tyr194 OH (2.63 A); NMR and '
                        'neutron: 100 % N(delta)-H, pKa < 3.5; the entry carries HD1. The docking default (HIE) put '
                        'the bare ND1 against the carboxylate'),
    119: ('HIE', 'HIE', 'zinc ligand through ND1 (2.02 A): the hydrogen is on NE2 and donates to Glu117 OE2 (2.72 A); '
                        'NMR and neutron: 100 % N(eps)-H'),
    122: ('HID', 'HIE', 'buried; ND1-H donates to the carbonyl O of Ala142 (3.02 A) and Leu120 (3.02 A), NE2 accepts '
                        'from Tyr51 OH (2.74 A); NMR and neutron: 100 % N(delta)-H, pKa < 3.5; the entry carries HD1'),
}
OTHER = {
    'GLU106': ('GLU', 'GLH', 'deprotonated. PROPKA gave 8.91 only because the zinc was left out of its input; with the '
                             'zinc it gives 5.33. OE1 accepts from Thr199 OG1 (2.62 A; the riding hydrogen of the '
                             'entry points at it, 1.81 A, 162 deg), OE2 from the N-H of Arg246 (2.79 A) and from '
                             'HOH415 (2.65 A); the zinc is 4.00 A from OE1; C-O 1.246 and 1.249 A. Thr199 OG1 must '
                             'donate to Glu106 because it accepts the sulfonamide N-H (2.86 A), the scheme seen by '
                             'neutron diffraction with acetazolamide (Fisher 2012)'),
    'GLU117': ('GLU', 'GLU', 'deprotonated (PROPKA 0.73 without, 3.60 with the zinc): accepts from His119 NE2-H (2.72 A), '
                             'His107 ND1-H (2.76 A) and the N-H of Glu106 and His107 (2.83 A each)'),
    'TYR7': ('TYR', 'TYR', 'neutral; PROPKA 13.5-14.7; protonated at pH 7.8 and deprotonated at pH 10 in the neutron '
                           'structures (Fisher 2011). Conformer B (0.61); OH...HOH415 2.58 A'),
    'CYS206': ('CYS', 'CYS', 'the only cysteine: free thiol (PROPKA 10.7). In the crystal it carries mercury '
                             '(4-mercuribenzoate, component BE7, 310, occupancy 0.76, with conformer A; Hg 302, 0.24, with '
                             'conformer B), a crystallisation additive that is removed; conformer A is taken; 14 A from the zinc'),
    'other Asp/Glu': ('ASP/GLU', 'ASP/GLU', 'all deprotonated: PROPKA pKa at most 4.83 (Glu187) apart from Glu106'),
    'Lys/Arg': ('LYS/ARG', 'LYS/ARG', 'all protonated: lowest PROPKA pKa 9.49 (Lys154) and 10.38 (Arg246)'),
    'Tyr': ('TYR', 'TYR', 'all neutral: lowest PROPKA pKa 10.43 (Tyr40)'),
    'N-terminus': ('ACE cap on His3', 'NH3+ on His3',
                   'His3 is the first residue with density; seven residues of the construct before it (GSPEF of the '
                   'tag, Met1, Ser2) are not modelled. They are not built. A charged amino group on His3 would be a '
                   'charge the protein does not have there, so His3 is capped with an acetyl group (three atoms '
                   'placed here, no density)'),
    'C-terminus': ('COO- on Lys261', 'COO- on Lys261', 'Lys261 is the true C-terminus and OXT is in the entry'),
}
# NMR values quoted from the two papers (read in full), used for the fractions at pH 7.4 below.
# (pKa, Hill coefficient or None, per cent N(eps)-H of the neutral form)
NMR_HIS = {
    'source': 'Fisher, Raum, Weininger 2025 (10.1002/cbic.202400930), Table 1, NMR at 308 K; His64: Raum, Fisher, '
              'Weininger 2023 (10.1007/s00018-023-04936-z), Table 1',
    'HIS3': (6.1, 0.62, '71-85'), 'HIS4': (6.2, 0.74, '60-72'), 'HIS10': (5.6, 0.79, '79-93'),
    'HIS15': (5.6, 0.92, '76-90'), 'HIS17': (5.6, None, '30-70'), 'HIS36': (6.94, None, '57-69'),
    'HIS64 in': (6.25, None, '44-56'), 'HIS64 out': (7.60, None, '44-56'),
    'HIS64 in, acetazolamide bound': (6.27, None, None), 'HIS64 out, acetazolamide bound': (7.43, None, None),
    'HIS94': ('<3.5', None, '0'), 'HIS96': ('<3.5', None, '0'), 'HIS107': ('<3.5', None, '0'),
    'HIS119': ('<3.5', None, '100'), 'HIS122': ('<3.5', None, '0'),
}


def nmr_histidines():
    """The quoted NMR values with the protonated fraction at pH 7.4 (Henderson-Hasselbalch, with the
    Hill coefficient where the paper gives one)."""
    out = dict(source=NMR_HIS['source'], pH=PH, values={})
    for k, v in NMR_HIS.items():
        if k == 'source':
            continue
        pka, hill, eps = v
        frac = None if isinstance(pka, str) else round(100.0 / (1.0 + 10 ** ((hill or 1.0) * (PH - pka))), 1)
        out['values'][k] = dict(pKa=pka, hill=hill, percent_N_eps_H_of_neutral=eps, percent_protonated_at_pH_7p4=frac)
    a, b = out['values']['HIS64 in']['percent_protonated_at_pH_7p4'], out['values']['HIS64 out']['percent_protonated_at_pH_7p4']
    c, d = (out['values']['HIS64 in, acetazolamide bound']['percent_protonated_at_pH_7p4'],
            out['values']['HIS64 out, acetazolamide bound']['percent_protonated_at_pH_7p4'])
    out['HIS64_percent_protonated_if_in_and_out_are_1_to_1'] = dict(free=round((a + b) / 2, 1), acetazolamide_bound=round((c + d) / 2, 1))
    return out


REMOVED_HETERO = {
    302: 'Hg, occupancy 0.24, on Cys206 conformer B: crystallisation additive',
    303: 'Hg, occupancy 0.05, at His64: crystallisation additive',
    305: 'second sulfanilamide, occupancy 0.68, on the surface at His3/His4/His10/His15, 14 A from the zinc-bound one',
    306: 'third sulfanilamide, occupancy 0.85, at the rim of the active site (Asn67, Glu69, Gln92) and in contact with '
         'a neighbouring molecule of the crystal, 5.7 A from the zinc-bound one; the simulation holds one ligand',
    307: 'glucose (cryoprotectant), occupancy 0.74', 308: 'glucose, occupancy 0.85', 309: 'glucose, occupancy 0.54',
    310: '4-mercuribenzoate on Cys206 (component BE7, p-chloromercuribenzoic acid, modelled without its chlorine), occupancy 0.76: additive',
}
WATER_CLASH = 2.2          # A: a water O this close to a ligand heavy atom is removed for that start
WATER_CLOSE = 2.6          # A: listed as close but kept
NH_ANION = 1.021           # A, N-H of the sulfonamidate in the paper's DFT geometry (M5); checked in main()
TETRA_HALF = 54.75         # deg, half of the tetrahedral angle

# Selection rule of the two docked starts.  It is the rule proposed in the task text; it was written down
# here after all modes had been tabulated (md/system/audit/poses_ad4zn.csv), so it was not fixed blind.
# The clause on Asn67 separates the ring-rotated orientation of the paper from other poses beyond 2 A;
# main() records that the pick is the same without it.
POSE_RULE = dict(
    job='dock/T3_6RL9/runs/sfa_neutral__ad4zn__prod (AutoDock4Zn, neutral sulfanilamide docked from its '
        'B3LYP-D3(BJ)/IEF-PCM geometry, exhaustiveness 32, 20 seeds: the production job of the paper)',
    site_check='Zn-N(sulfonamide) <= 2.3 A and Thr199 OG1...N <= 3.3 A (the site check of the paper)',
    N='the best-scored mode, over all modes of all seeds, that passes the site check and lies within 1.0 A of '
      'the crystal pose (heavy atoms, symmetry-corrected, no superposition); ties: lower seed, then lower mode',
    R='the best-scored mode, over all modes of all seeds, that passes the site check, lies more than 2.0 A from '
      'the crystal pose and has the aniline N within 3.5 A of Asn67 OD1 (the ring-rotated orientation of '
      'Section 3.6); ties: lower seed, then lower mode',
)


# ----------------------------------------------------------------------------------------------------
def build_protein(chain, choice):
    """Atoms of the protein: list of (resname, seqnum, [(name, xyz, occ, B)])."""
    out, missing = [], {}
    for r in C.polymer_residues(chain):
        n = r.seqid.num
        atoms = C.conformer_atoms(r, choice.get(n, '-'))
        name = HIS[n][0] if r.name == 'HIS' else r.name
        have = {a.name for a in atoms}
        lack = [x for x in C.HEAVY[r.name] if x not in have]
        if lack:
            missing[f'{r.name}{n}'] = lack
        extra = sorted(have - set(C.HEAVY[r.name]) - {'OXT'})
        assert not extra, (r.name, n, extra)
        order = {x: i for i, x in enumerate(C.HEAVY[r.name] + ['OXT'])}
        atoms.sort(key=lambda a: order[a.name])
        out.append((name, n, [(a.name, C.xyz(a), a.occ, a.b_iso) for a in atoms]))
    return out, missing


def ace_cap(protein, others):
    """Acetyl cap on the first residue: C on N (1.335 A, C-N-CA 121.9 deg), O and CH3 in the peptide
    plane (trans peptide).  The one free angle, phi of the first residue, is scanned in 10 deg steps and
    the value that keeps the three atoms farthest from every other heavy atom is taken."""
    first = {nm: p for nm, p, _, _ in protein[0][2]}
    N, CA, Cc = first['N'], first['CA'], first['C']
    best = None
    for phi in range(-180, 180, 10):
        c = C.place(Cc, CA, N, 1.335, 121.9, phi)
        o = C.place(CA, N, c, 1.229, 122.9, 0.0)
        m = C.place(CA, N, c, 1.522, 116.6, 180.0)
        d = min(float(np.min(np.linalg.norm(others - p, axis=1))) for p in (c, o, m))
        if best is None or d > best[0] + 1e-9:
            best = (d, phi, c, o, m)
    d, phi, c, o, m = best
    return ('ACE', protein[0][1] - 1, [('CH3', m, 0.0, 0.0), ('C', c, 0.0, 0.0), ('O', o, 0.0, 0.0)]), \
        dict(phi_His3_deg=phi, min_distance_to_other_heavy_atoms_A=round(d, 2),
             geometry='C-N 1.335 A, C-N-CA 121.9 deg, C=O 1.229 A, N-C-O 122.9 deg, C-CH3 1.522 A, N-C-CH3 116.6 deg, omega 180 deg',
             note='three atoms without density (written with occupancy 0.00)')


def choose_waters(chain, protein_xyz):
    """Crystal water kept before the ligand is considered.  A copy with a label is kept when it does not
    come within COUPLING_CUT of an atom of the conformers taken; of several such copies the one with the
    highest occupancy; a site with no free copy is dropped.  Water without a label is kept."""
    kept, record = [], []
    for r in C.waters(chain):
        cops = []
        for a in r:
            d = float(np.min(np.linalg.norm(protein_xyz - C.xyz(a), axis=1)))
            cops.append((a, d))
        if not C.labels_of(r):
            a, d = cops[0]
            assert len(cops) == 1
            kept.append((r.seqid.num, a))
            if d < C.COUPLING_CUT or a.occ < 0.995:
                record.append(dict(water=r.seqid.num, label='-', occupancy=round(a.occ, 2), kept=True,
                                   min_dist_protein_A=round(d, 2), note='unlabelled water closer than 2.4 A to the protein or partly occupied'))
            continue
        free = [(a, d) for a, d in cops if d >= C.COUPLING_CUT]
        pick = max(free, key=lambda t: (round(t[0].occ, 3), -ord(C.alt(t[0])))) if free else None
        for a, d in cops:
            record.append(dict(water=r.seqid.num, label=C.alt(a), occupancy=round(a.occ, 2),
                               min_dist_protein_A=round(d, 2), kept=bool(pick and a is pick[0]),
                               note=('clashes with a conformer taken' if d < C.COUPLING_CUT else
                                     ('kept' if pick and a is pick[0] else 'another copy has the higher occupancy'))))
        if pick:
            kept.append((r.seqid.num, pick[0]))
    return kept, record


def sulfonamide_h(M, site):
    """The one hydrogen of the zinc-bound sulfonamidate nitrogen: N2 is taken as tetrahedral with S and
    Zn as two of its substituents; of the two remaining directions the one nearer to Thr199 OG1 gets the
    hydrogen (N-H 1.021 A, the paper's DFT value for the anion), the other is the lone pair."""
    n2 = M['N2']
    u1 = (M['S'] - n2) / np.linalg.norm(M['S'] - n2)
    u2 = (site['ZN'] - n2) / np.linalg.norm(site['ZN'] - n2)
    b = -(u1 + u2)
    b /= np.linalg.norm(b)
    nrm = np.cross(u1, u2)
    nrm /= np.linalg.norm(nrm)
    t = np.radians(TETRA_HALF)
    cands = [n2 + NH_ANION * (np.cos(t) * b + s * np.sin(t) * nrm) for s in (1, -1)]
    cands.sort(key=lambda h: C.dist(h, site['OG1']))
    return cands[0], cands[1]


def ligand_mol(M, HN21, aniline_h=None, title=''):
    """RDKit molecule of the anion ArSO2NH-, atoms in the order C.LIG_ORDER (the order of the ligand
    parameter files).  Ring hydrogens (and, for the crystal pose, the aniline hydrogens) are placed by
    RDKit on the given heavy atoms; aniline_h, when given, are the two docked positions."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    rw = Chem.RWMol()
    idx = {}
    for n in C.SAN_HEAVY:
        idx[n] = rw.AddAtom(Chem.Atom('S' if n == 'S' else n[0]))
    for a, b, o in C.SAN_BONDS:
        rw.AddBond(idx[a], idx[b], Chem.BondType.DOUBLE if o == 2 else Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(idx['N2']).SetFormalCharge(-1)
    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for n in C.SAN_HEAVY:
        conf.SetAtomPosition(idx[n], Point3D(*[float(v) for v in M[n]]))
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    mh = Chem.AddHs(mol, addCoords=True)
    assert mh.GetNumAtoms() == 18, mh.GetNumAtoms()
    conf = mh.GetConformer()
    hs = {}
    for a in mh.GetAtoms():
        if a.GetAtomicNum() == 1:
            hs.setdefault(C.SAN_HEAVY[a.GetNeighbors()[0].GetIdx()], []).append(a.GetIdx())
    assert {k: len(v) for k, v in hs.items()} == {'N1': 2, 'C2': 1, 'C3': 1, 'C5': 1, 'C6': 1, 'N2': 1}, hs
    conf.SetAtomPosition(hs['N2'][0], Point3D(*[float(v) for v in HN21]))
    if aniline_h is not None:
        for i, p in zip(hs['N1'], aniline_h):
            conf.SetAtomPosition(i, Point3D(*[float(v) for v in p]))
    index = {n: i for i, n in enumerate(C.SAN_HEAVY)}
    index.update(HN11=hs['N1'][0], HN12=hs['N1'][1], H2=hs['C2'][0], H3=hs['C3'][0], H5=hs['C5'][0],
                 H6=hs['C6'][0], HN21=hs['N2'][0])
    names = list(C.LIG_ORDER)
    mh = Chem.RenumberAtoms(mh, [index[n] for n in names])
    for a, nm in zip(mh.GetAtoms(), names):
        a.SetProp('_Name', nm)
    mh.SetProp('_Name', title)
    Chem.SanitizeMol(mh)
    assert Chem.GetFormalCharge(mh) == -1
    return mh, names


def write_ligand(tag, mol, names, props):
    from rdkit import Chem
    mol.SetProp('atom_names', ' '.join(names))
    for k, v in props.items():
        mol.SetProp(k, str(v))
    block = Chem.MolToMolBlock(mol)
    lines = block.split('\n')
    lines[1] = '  md/build_system.py'        # the program line of RDKit carries no date, keep it fixed anyway
    sdf = '\n'.join(lines)
    for k in ['atom_names'] + list(props):
        sdf += f'>  <{k}>\n{mol.GetProp(k)}\n\n'
    sdf += '$$$$\n'
    (C.OUT / f'ligand_{tag}.sdf').write_text(sdf)
    conf = mol.GetConformer()
    out = [f'REMARK   start {tag}: sulfonamidate anion of sulfanilamide (charge -1), frame of 6RL9',
           f'REMARK   {props["source"]}']
    for i, nm in enumerate(names):
        p = conf.GetAtomPosition(i)
        el = mol.GetAtomWithIdx(i).GetSymbol()
        out.append(C.pdb_line('HETATM', i + 1, nm, 'SAN', 'A', C.LIG_SEQ, (p.x, p.y, p.z), 1.0, 0.0, el))
    out.append('END')
    (C.OUT / f'ligand_{tag}.pdb').write_text('\n'.join(out) + '\n')
    return {nm: np.array(list(conf.GetAtomPosition(i))) for i, nm in enumerate(names)}


def pick_poses(crystal, site):
    """Apply POSE_RULE to every mode of the production job."""
    rows = []
    for f in sorted(C.PROD_JOB.glob('seed*/out.pdbqt')):
        seed = int(f.parent.name[4:])
        for i, mode in enumerate(C.read_pdbqt_modes(f), 1):
            rmsd, M, mp = C.name_docked_atoms(mode, crystal)
            mt = C.pose_metrics(M, site)
            rows.append(dict(seed=seed, mode=i, score=mode['score'], rmsd=rmsd, M=M, raw=mode, map=mp, file=f, **mt))
    ok = [r for r in rows if r['Zn_N2'] <= 2.3 and r['N2_OG1'] <= 3.3]
    key = lambda r: (r['score'], r['seed'], r['mode'])  # noqa: E731
    near = sorted([r for r in ok if r['rmsd'] <= 1.0], key=key)
    rot = sorted([r for r in ok if r['rmsd'] > 2.0 and r['N1_Asn67OD1'] <= 3.5], key=key)
    far = sorted([r for r in ok if r['rmsd'] > 2.0], key=key)
    return near[0], rot[0], dict(modes=len(rows), pass_site_check=len(ok), near_native=len(near),
                                 beyond_2A=len(far), rotated_at_Asn67=len(rot), seeds=len({r['seed'] for r in rows}),
                                 best_mode_beyond_2A_without_the_Asn67_clause=dict(seed=far[0]['seed'], mode=far[0]['mode'], score=far[0]['score']),
                                 R_is_the_same_without_the_Asn67_clause=bool(far[0] is rot[0]),
                                 runner_up_N=dict(seed=near[1]['seed'], mode=near[1]['mode'], score=near[1]['score']),
                                 runner_up_R=dict(seed=rot[1]['seed'], mode=rot[1]['mode'], score=rot[1]['score']),
                                 rotated_score_range=[rot[-1]['score'], rot[0]['score']],
                                 rotated_N1_Phe131CZ_range=[round(min(r['N1_Phe131CZ'] for r in rot), 2), round(max(r['N1_Phe131CZ'] for r in rot), 2)],
                                 rotated_N1_Asn67OD1_range=[round(min(r['N1_Asn67OD1'] for r in rot), 2), round(max(r['N1_Asn67OD1'] for r in rot), 2)])


def docking_frame_atoms():
    """Phe131 CZ and Asn67 OD1 in the receptor the docking used (hydrogens added and minimised with the
    deposited heavy atoms restrained), to compare the orientation distances in the two frames."""
    out = {}
    for ln in C.DOCK_RECEPTOR.read_text().splitlines():
        if ln.startswith('ATOM'):
            key = (ln[17:20].strip(), int(ln[22:26]), ln[12:16].strip())
            if key in (('PHE', 131, 'CZ'), ('ASN', 67, 'OD1')):
                out[key[2]] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return out


def dft_anion_nh():
    """N-H distance of the sulfonamide nitrogen in the paper's optimised anion (opt.xyz of M5)."""
    rows = [ln.split() for ln in C.DFT_ANION.read_text().splitlines()[2:] if len(ln.split()) == 4]
    el = [r[0] for r in rows]
    X = np.array([[float(v) for v in r[1:]] for r in rows])
    s = el.index('S')
    n = [i for i, e in enumerate(el) if e == 'N' and C.dist(X[i], X[s]) < 1.8]
    h = [i for i, e in enumerate(el) if e == 'H' and C.dist(X[i], X[n[0]]) < 1.2]
    assert len(n) == 1 and len(h) == 1, (n, h)
    return C.dist(X[n[0]], X[h[0]]), C.angle(X[s], X[n[0]], X[h[0]]), C.dist(X[s], X[n[0]])


def main():
    C.OUT.mkdir(parents=True, exist_ok=True)
    nh, snh, sn_dft = dft_anion_nh()
    assert abs(nh - NH_ANION) < 5e-4, nh
    st, chain = C.read_deposited()
    choice, groups = C.choose_conformers(chain)
    dock_choice = C.docking_conformers()

    # ---- protein ------------------------------------------------------------------------------------
    protein, missing = build_protein(chain, choice)
    prot_xyz = np.array([p for _, _, atoms in protein for _, p, _, _ in atoms])
    kept_w, water_record = choose_waters(chain, prot_xyz)
    wat_xyz = np.array([C.xyz(a) for _, a in kept_w])
    zn = [a for r in chain if r.seqid.num == C.ZN_SEQ for a in r][0]
    first_n = protein[0][1]
    others = np.vstack([np.array([p for _, n, atoms in protein for nm, p, _, _ in atoms
                                  if not (n == first_n and nm in ('N', 'CA'))]), wat_xyz])
    ace, ace_info = ace_cap(protein, others)
    protein_out = [ace] + protein

    # every taken residue pair that is bonded in sequence must have a sane peptide bond; no two
    # non-bonded heavy atoms of the protein may be closer than 2.4 A
    pept = []
    for (n1, s1, a1), (n2, s2, a2) in zip(protein_out, protein_out[1:]):
        c = dict((nm, p) for nm, p, _, _ in a1)['C']
        n = dict((nm, p) for nm, p, _, _ in a2)['N']
        pept.append((C.dist(c, n), s1, s2))
    assert all(1.28 < d < 1.40 for d, _, _ in pept), [p for p in pept if not 1.28 < p[0] < 1.40]

    lines = ['REMARK   hCA II (PDB 6RL9, chain A) for molecular dynamics; written by md/build_system.py',
             'REMARK   heavy atoms only, one conformer per residue, Amber residue names; ACE cap on His3',
             'REMARK   (occupancy 0.00: placed without density).  Zinc, ligand and water are in separate files.']
    serial = 0
    res_map = []
    for k, (name, n, atoms) in enumerate(protein_out, 1):
        for nm, p, occ, b in atoms:
            serial += 1
            lines.append(C.pdb_line('ATOM', serial, nm, name, 'A', n, p, occ, b))
        orig = 'ACE (added)' if name == 'ACE' else [r.name for r in C.polymer_residues(chain) if r.seqid.num == n][0]
        res_map.append(dict(deposited_number=n, deposited_name=orig, name_md=name, conformer=choice.get(n, ''),
                            tleap_index=k))
    lines += ['TER', 'END']
    (C.OUT / 'protein_md.pdb').write_text('\n'.join(lines) + '\n')
    (C.OUT / 'zinc.pdb').write_text('\n'.join([
        'REMARK   zinc ion of 6RL9 chain A (ZN 301, occupancy 1.00)',
        C.pdb_line('HETATM', 1, 'ZN', 'ZN', 'A', C.ZN_SEQ, C.xyz(zn), zn.occ, zn.b_iso, 'ZN'), 'END']) + '\n')
    with open(C.OUT / 'residue_map.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(res_map[0].keys()), lineterminator='\n')
        w.writeheader()
        w.writerows(res_map)

    # ---- ligand starts ------------------------------------------------------------------------------
    site = C.site_atoms(chain, choice)
    crystal = {a.name: C.xyz(a) for r in chain if r.seqid.num == C.LIG_SEQ for a in r}
    assert list(crystal) == C.SAN_HEAVY
    pN, pR, pose_stats = pick_poses(crystal, site)
    dock_frame = docking_frame_atoms()
    res_by_num = {r.seqid.num: r for r in chain}
    his64 = {L: [C.xyz(a) for a in res_by_num[64] if C.alt(a) == L] for L in 'ABC'}
    prot_names = [(name, n, nm) for name, n, atoms in protein_out for nm, _, _, _ in atoms]
    prot_all = np.array([p for _, _, atoms in protein_out for _, p, _, _ in atoms])

    starts, water_sets, lig_xyz = {}, {}, {}
    for tag, pose in (('X', None), ('N', pN), ('R', pR)):
        if pose is None:
            M = crystal
            src = 'crystal pose: SAN 304 of dock/pdb/6RL9.cif.gz (occupancy 1.00, no alternate conformer)'
            aniline = None
            h_report = dict(note='the entry holds no ligand hydrogens; aniline and ring hydrogens placed by RDKit '
                                 '(aniline planar), the sulfonamide hydrogen by the tetrahedral construction')
            score = None
        else:
            M = pose['M']
            rel = pose['file'].relative_to(C.REPO)
            src = f'docked pose: {rel}, MODEL {pose["mode"]} (seed {pose["seed"]}), neutral molecule, AutoDock4Zn score {pose["score"]:.3f} kcal/mol'
            raw = pose['raw']['atoms']
            aniline = [raw[i][1] for i in C.PDBQT_ANILINE_H]
            # the two sulfonamide hydrogens of the docked neutral molecule
            hrep = []
            for i in C.PDBQT_SULFONAMIDE_H:
                h = raw[i][1]
                hrep.append(dict(pdbqt_atom=i, H_Zn=round(C.dist(h, site['ZN']), 2),
                                 Zn_N2_H=round(C.angle(site['ZN'], M['N2'], h), 1),
                                 H_OG1=round(C.dist(h, site['OG1']), 2),
                                 N2_H_OG1=round(C.angle(M['N2'], h, site['OG1']), 1),
                                 S_N2_H=round(C.angle(M['S'], M['N2'], h), 1)))
            to_zn = min(hrep, key=lambda d: d['H_Zn'])
            other = [d for d in hrep if d is not to_zn][0]
            donates = other['H_OG1'] <= 2.5 and other['N2_H_OG1'] >= 120.0
            h_report = dict(docked_sulfonamide_H=hrep, deleted=to_zn['pdbqt_atom'], remaining=other['pdbqt_atom'],
                            remaining_donates_to_Thr199=bool(donates),
                            criterion='H...OG1 <= 2.5 A and N-H...OG1 >= 120 deg')
            score = pose['score']
        HN21, lone = sulfonamide_h(M, site)
        mol, names = ligand_mol(M, HN21, aniline, title=f'6RL9 start {tag}')
        if pose is not None:
            other_xyz = pose['raw']['atoms'][h_report['remaining']][1]
            h_report['constructed_H_from_remaining_docked_H_A'] = round(C.dist(HN21, other_xyz), 2)
            h_report['note'] = ('neither docked hydrogen meets the hydrogen-bond criterion to Thr199 OG1; the one '
                                'pointing at the zinc is deleted and the other replaced by the tetrahedral construction'
                                if not donates else 'the remaining docked hydrogen donates to Thr199 OG1; it is '
                                'nevertheless replaced by the tetrahedral construction so that the three starts are built alike')
        h_report['second_tetrahedral_position_H_OG1_A'] = round(C.dist(lone, site['OG1']), 2)
        props = dict(source=src, charge=-1, start=tag)
        L = write_ligand(tag, mol, names, props)
        lig_xyz[tag] = L
        heavy = np.array([L[n] for n in C.SAN_HEAVY])
        mt = C.pose_metrics(L, site, H=L['HN21'])
        rmsd = 0.0 if pose is None else pose['rmsd']
        # contacts
        d = np.linalg.norm(prot_all[:, None, :] - heavy[None, :, :], axis=2)
        i, j = np.unravel_index(d.argmin(), d.shape)
        dz = [(C.dist(C.xyz(zn), L[n]), n) for n in C.SAN_HEAVY if n != 'N2']
        wd = np.linalg.norm(wat_xyz[:, None, :] - heavy[None, :, :], axis=2)
        removed = [(kept_w[k][0], round(float(wd[k].min()), 2), C.SAN_HEAVY[int(wd[k].argmin())])
                   for k in range(len(kept_w)) if wd[k].min() < WATER_CLASH]
        close = [(kept_w[k][0], round(float(wd[k].min()), 2), C.SAN_HEAVY[int(wd[k].argmin())])
                 for k in range(len(kept_w)) if WATER_CLASH <= wd[k].min() < WATER_CLOSE]
        n1 = L['N1']
        near_n1 = sorted((round(C.dist(n1, wat_xyz[k]), 2), kept_w[k][0]) for k in range(len(kept_w))
                         if C.dist(n1, wat_xyz[k]) <= 4.0)
        lig_shell = sorted((round(float(wd[k].min()), 2), kept_w[k][0], C.SAN_HEAVY[int(wd[k].argmin())])
                           for k in range(len(kept_w)) if wd[k].min() <= 3.6)
        water_sets[tag] = [w for w in kept_w if w[0] not in {r[0] for r in removed}]
        starts[tag] = dict(
            source=src, score_kcal_mol=score, rmsd_to_crystal_A=round(rmsd, 3),
            **{k: round(v, 3) for k, v in mt.items()},
            N1_Phe131CZ_docking_frame=round(C.dist(n1, dock_frame['CZ']), 3),
            N1_Asn67OD1_docking_frame=round(C.dist(n1, dock_frame['OD1']), 3),
            N1_to_His64_ring={Lb: round(min(C.dist(n1, p) for p in his64[Lb]), 2) for Lb in 'ABC'},
            min_heavy_distance_to_protein=dict(A=round(float(d[i, j]), 2), protein_atom='%s%d.%s' % prot_names[i],
                                               ligand_atom=C.SAN_HEAVY[j]),
            min_distance_Zn_to_other_ligand_atom=dict(A=round(min(dz)[0], 2), atom=min(dz)[1]),
            sulfonamide_hydrogen=h_report,
            waters_removed=[dict(water=w, dist_A=dd, ligand_atom=a) for w, dd, a in removed],
            waters_close_kept=[dict(water=w, dist_A=dd, ligand_atom=a) for w, dd, a in close],
            waters_within_4A_of_aniline_N=[dict(water=w, dist_A=dd) for dd, w in near_n1],
            waters_within_3p6A_of_ligand=[dict(water=w, dist_A=dd, ligand_atom=a) for dd, w, a in lig_shell],
            n_waters=len(water_sets[tag]),
        )

    # ---- water files --------------------------------------------------------------------------------
    def write_waters(name, ws, remark):
        out = [f'REMARK   {remark}', 'REMARK   oxygen atoms only, numbers as in 6RL9; tleap adds the hydrogens (TIP3P)']
        for k, (num, a) in enumerate(ws, 1):
            out.append(C.pdb_line('HETATM', k, 'O', 'HOH', 'W', num, C.xyz(a), a.occ, a.b_iso, 'O'))
        out.append('END')
        (C.OUT / name).write_text('\n'.join(out) + '\n')
    write_waters('waters_all.pdb', kept_w, f'crystal water of 6RL9 chain A kept for MD: {len(kept_w)} of {len(C.waters(chain))} sites')
    for tag in 'XNR':
        write_waters(f'waters_{tag}.pdb', water_sets[tag],
                     f'water of start {tag}: {len(water_sets[tag])} molecules (water within {WATER_CLASH} A of a ligand heavy atom removed)')

    # deep water and first shell of the zinc
    t199n = site['T199N']
    deep = sorted((round(C.dist(t199n, wat_xyz[k]), 2), kept_w[k][0]) for k in range(len(kept_w)) if C.dist(t199n, wat_xyz[k]) <= 4.0)
    znw = sorted((round(C.dist(C.xyz(zn), wat_xyz[k]), 2), kept_w[k][0]) for k in range(len(kept_w)) if C.dist(C.xyz(zn), wat_xyz[k]) <= 5.0)
    ww = []
    for i in range(len(kept_w)):
        dd = np.linalg.norm(wat_xyz[i + 1:] - wat_xyz[i], axis=1)
        for j in np.where(dd < 2.4)[0]:
            ww.append(dict(waters=[kept_w[i][0], kept_w[i + 1 + int(j)][0]], dist_A=round(float(dd[j]), 2)))
    dwp = np.linalg.norm(prot_all[:, None, :] - wat_xyz[None, :, :], axis=2)
    wprot = [dict(water=kept_w[k][0], dist_A=round(float(dwp[:, k].min()), 2), protein_atom='%s%d.%s' % prot_names[int(dwp[:, k].argmin())])
             for k in range(len(kept_w)) if dwp[:, k].min() < 2.4]

    # ---- starts.csv ---------------------------------------------------------------------------------
    cols = ['start', 'score_kcal_mol', 'rmsd_to_crystal_A', 'Zn_N2', 'Zn_O_min', 'Zn_O1', 'Zn_O2', 'N2_OG1', 'H_OG1',
            'N2_H_OG1', 'T199N_O_min', 'C3_C4_S_N2', 'C5_C4_S_N2', 'tau_fold', 'C4_S_N2', 'S_N2', 'S_N2_Zn',
            'N1_Phe131CZ', 'N1_Asn67OD1', 'N1_Asn67ND2', 'N1_Phe131CZ_docking_frame', 'N1_Asn67OD1_docking_frame',
            'n_waters']
    with open(C.OUT / 'starts.csv', 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(cols)
        for tag in 'XNR':
            w.writerow([tag] + [starts[tag][c] for c in cols[1:]])

    # ---- comparison with the numbers fixed in section 5 of md/PLAN.md ---------------------------------
    plan = dict(N1_Phe131CZ=dict(crystal=3.9, rotated=[9.3, 9.8]), N1_Asn67OD1=dict(crystal=7.7, rotated=[2.6, 2.7]))
    cmp_plan = dict(
        plan=plan,
        X=dict(N1_Phe131CZ=starts['X']['N1_Phe131CZ'], N1_Asn67OD1=starts['X']['N1_Asn67OD1'],
               matches=bool(round(starts['X']['N1_Phe131CZ'], 1) == 3.9 and round(starts['X']['N1_Asn67OD1'], 1) == 7.7)),
        R=dict(N1_Phe131CZ=starts['R']['N1_Phe131CZ'], N1_Asn67OD1=starts['R']['N1_Asn67OD1'],
               matches=bool(9.3 <= round(starts['R']['N1_Phe131CZ'], 1) <= 9.8 and 2.6 <= round(starts['R']['N1_Asn67OD1'], 1) <= 2.7)),
        N=dict(N1_Phe131CZ=starts['N']['N1_Phe131CZ'], N1_Asn67OD1=starts['N']['N1_Asn67OD1'],
               note='the plan gives no number for the near-native start'),
        midpoints_X_to_R=dict(N1_Phe131CZ=round((starts['X']['N1_Phe131CZ'] + starts['R']['N1_Phe131CZ']) / 2, 2),
                              N1_Asn67OD1=round((starts['X']['N1_Asn67OD1'] + starts['R']['N1_Asn67OD1']) / 2, 2)),
        all_rotated_modes_of_the_production_job=dict(N1_Phe131CZ=pose_stats['rotated_N1_Phe131CZ_range'],
                                                     N1_Asn67OD1=pose_stats['rotated_N1_Asn67OD1_range']),
    )

    # ---- decisions.json -----------------------------------------------------------------------------
    n_res = len(protein)
    counts = {}
    for name, _, _ in protein_out:
        counts[name] = counts.get(name, 0) + 1
    neg = counts.get('ASP', 0) + counts.get('GLU', 0) + 1
    pos = counts.get('LYS', 0) + counts.get('ARG', 0) + counts.get('HIP', 0)
    conf_rows = []
    for g in groups:
        nums = [int(''.join(ch for ch in r if ch.isdigit())) for r in g['residues']]
        conf_rows.append(dict(residues=g['residues'], label=g['label'], occupancy=g['occupancy'],
                              heavy_atoms=g['heavy_atoms'], decided_by=g['decided_by'], coupling=g['coupling'],
                              docking_labels=[dock_choice.get(n) for n in nums],
                              differs_from_docking=[g['residues'][i] for i, n in enumerate(nums) if dock_choice.get(n) != g['label']]))
    diffs = (
        [f'His{n}: {v[1]} in the docking receptor, {v[0]} here' for n, v in HIS.items() if v[0] != v[1]]
        + ['Glu106: GLH (protonated) in the docking receptor, GLU (deprotonated) here',
           'N-terminus: NH3+ on His3 in the docking receptor, acetyl cap here (net charge of the protein -1 against the docking receptor from this, -1 from Glu106)',
           'conformers: ' + ', '.join(r for c in conf_rows for r in c['differs_from_docking']) + ' take another label than in the docking receptor (coupled groups take one label)',
           f'water: none in the docking receptor, {len(kept_w)} crystal waters here',
           'side-chain atoms without density: built by PDBFixer and minimised in the docking receptor, built by tleap here',
           'hydrogens: placed by PDBFixer and minimised in vacuum in the docking receptor (heavy atoms moved by up to 0.175 A); here the deposited heavy-atom coordinates are used unchanged and tleap adds the hydrogens',
           'ligand: neutral sulfanilamide in the AutoDock4Zn docking, the sulfonamidate anion here'])
    decisions = dict(
        written_by='md/build_system.py', stage='S3 (system builder, part 1) of md/PLAN.md', pH=PH,
        sources={str(p.relative_to(C.REPO)): C.sha256(p) for p in
                 (C.CIF, C.DOCK_RECEPTOR, C.DOCK / 'T3_6RL9' / 'prep_report.json', pN['file'], pR['file'], C.DFT_ANION)},
        software=dict(python=sys.version.split()[0], gemmi=__import__('gemmi').__version__,
                      rdkit=__import__('rdkit').__version__, numpy=np.__version__),
        random_numbers='none',
        protein=dict(
            entry='6RL9', chain='A', residues=n_res, first='HIS3', last='LYS261',
            numbering='deposited numbers kept (no residue 126; Thr125 C - Lys127 N is a normal peptide bond)',
            residue_names=dict(sorted(counts.items())),
            formal_charge=dict(negative=neg, positive=pos, net=pos - neg,
                               note='Asp + Glu + C-terminus against Lys + Arg; no charged histidine; ACE cap neutral'),
            histidines={f'HIS{n}': dict(state=v[0], docking=v[1], differs_from_docking=v[0] != v[1], reason=v[2])
                        for n, v in HIS.items()},
            histidines_nmr=nmr_histidines(),
            other={k: dict(state=v[0], docking=v[1], differs_from_docking=v[0] != v[1], reason=v[2]) for k, v in OTHER.items()},
            ace_cap=ace_info,
            conformer_rule=C.choose_conformers.__doc__.split('Rule')[1].split('Returns')[0].strip(),
            conformers=conf_rows,
            heavy_atoms_not_in_the_entry=missing,
            expected_tleap_added_heavy_atoms=sorted(f'{k}.{a}' for k, v in missing.items() for a in v),
            n_heavy_atoms_written=serial,
        ),
        removed_hetero_groups={str(k): v for k, v in REMOVED_HETERO.items()},
        zinc=dict(file='zinc.pdb', note='ZN 301, occupancy 1.00, B 5.5 A2; parameters of the site are made in stage S2'),
        water=dict(
            rule=choose_waters.__doc__.strip(),
            sites_in_entry=len(C.waters(chain)), kept=len(kept_w),
            labelled_or_partial=water_record,
            dropped=[r['water'] for r in water_record if not r['kept'] and
                     not any(q['kept'] for q in water_record if q['water'] == r['water'])],
            ligand_rule=f'for each start a water whose O is within {WATER_CLASH} A of a ligand heavy atom is removed; water between {WATER_CLASH} and {WATER_CLOSE} A is listed and kept',
            kept_pairs_closer_than_2p4A=ww, kept_water_closer_than_2p4A_to_protein=wprot,
            waters_within_4A_of_Thr199_N=[dict(water=w, dist_A=d) for d, w in deep],
            waters_within_5A_of_zinc=[dict(water=w, dist_A=d) for d, w in znw],
            deep_water='no water lies within 4 A of Thr199 N: the sulfonyl O1 of the ligand takes the place of the deep water (Thr199 N...O1 2.94 A)',
            hydrogens='tleap places the two hydrogens of every water in one standard orientation; they are relaxed in the equilibration',
        ),
        ligand=dict(
            species='sulfonamidate anion ArSO2NH- of sulfanilamide, charge -1, 18 atoms',
            residue_name='SAN', residue_number=C.LIG_SEQ,
            atom_names=C.LIG_ORDER,
            atom_order='the order of the paper\'s DFT files and of md/ligand/sfa_anion.mol2 (stage S1); the names are those of SAN in 6RL9 (heavy atoms) and of the chemical component SAN (hydrogens; HN22 is the hydrogen the anion lacks)',
            hydrogen_parents=dict(C.SAN_H),
            paper_labels='C3-C4-S1-N2 of the paper is C3-C4-S-N2 here (S1 = S, N2 = sulfonamide N, N1 = aniline N)',
            docked_atom_names='the ring (C2/C6, C3/C5) and the two sulfonyl oxygens are equivalent: the assignment with the lowest RMSD to the crystal pose is taken',
            sulfonamide_hydrogen=sulfonamide_h.__doc__.strip(),
            dft_anion=dict(file=str(C.DFT_ANION.relative_to(C.REPO)), N_H_A=round(nh, 4), S_N_H_deg=round(snh, 2), S_N_A=round(sn_dft, 4)),
            other_hydrogens='ring hydrogens by RDKit on the given heavy atoms; aniline hydrogens of X by RDKit (planar NH2), of N and R as docked (pyramidal NH2 of the DFT geometry)',
            pose_rule=POSE_RULE, pose_statistics=pose_stats,
            pose_rule_note='the rule proposed in the task; written down after all modes had been tabulated, so not fixed blind; the pick of R is the same without the clause on Asn67',
            note='N and R keep the bond lengths of the docked NEUTRAL molecule (S-N 1.677 A; 1.596 A in the crystal pose, 1.580 A in the DFT anion) and the docked Zn-N distance (2.2 A); the force field relaxes both in the minimisation',
        ),
        starts=starts,
        plan_section_5_check=cmp_plan,
        differences_from_docking_receptor=diffs,
    )
    (C.OUT / 'decisions.json').write_text(json.dumps(decisions, indent=1) + '\n')

    print(f'protein: {n_res} residues + ACE, {serial} heavy atoms; net formal charge {pos - neg}')
    print(f'water: {len(kept_w)} of {len(C.waters(chain))} kept; per start:', {t: len(water_sets[t]) for t in 'XNR'})
    for tag in 'XNR':
        s = starts[tag]
        print(f'{tag}: score {s["score_kcal_mol"]}  RMSD {s["rmsd_to_crystal_A"]}  Zn-N {s["Zn_N2"]}  tau {s["tau_fold"]}  '
              f'C-S-N {s["C4_S_N2"]}  N1-Phe131CZ {s["N1_Phe131CZ"]}  N1-Asn67OD1 {s["N1_Asn67OD1"]}  removed water {[w["water"] for w in s["waters_removed"]]}')
    print('plan check: X', cmp_plan['X']['matches'], ' R', cmp_plan['R']['matches'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
