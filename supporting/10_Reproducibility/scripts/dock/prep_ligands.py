"""
Ligand preparation for the SAPH docking work.

    python dock/prep_ligands.py                       # all species, ETKDG conformers
    python dock/prep_ligands.py --only sfa_neutral     # rebuilds this species, keeps the others
    python dock/prep_ligands.py --dft sfa_neutral=molecules/sfa/results/sulfanilamide_pcm_d3/opt.xyz

Design decisions, and why:

*   **Species, not molecules.**  The crystal is a salt, but at pH 7.4 the anilinium of
    4-sulfamoylanilinium is fully deprotonated (pKa 2.12) and the sulfonamide is almost entirely
    neutral (pKa 10.55; Glockner 2020, 10.3390/biom10040509).  So the ligand that matters is neutral
    sulfanilamide; the sulfonamidate anion is the species that coordinates zinc; the cation is a
    sensitivity case labelled "as crystallised".  **The ion pair is never docked** -- the template
    paper docked one, and a docking box has no counter-ion to balance it.

*   **Geometry now, DFT geometry later.**  Burst 1a is still optimising the fragments, so the default
    input geometry is an ETKDG conformer (Riniker 2015, 10.1021/acs.jcim.5b00654) minimised with
    MMFF94s.  `--dft <species>=<file.xyz>` swaps in the B3LYP-D3(BJ)/6-311++G(d,p) PCM geometry when
    it lands: bonds are perceived from the coordinates and the bond orders and formal charges are then
    taken from the species SMILES, so the DFT file needs no atom-order convention.

*   **Redocking starts from a random conformer, never the crystal pose.**  Redocking that starts from
    the answer measures nothing.  Each seed gets its own ETKDG conformer (seed = the Vina seed), and
    `--randomise` records the RMSD of the starting conformer from the crystal pose so the report can
    show the search had real work to do.

*   **Protonation states are stated, not guessed.**  Each species carries an explicit SMILES with
    explicit formal charges, and the script checks that the neutral parent matches the RCSB chemical
    component for the ligands that have one (SAN, AZM, 08D, PAB, 78H).  A mismatch is a hard error --
    that check is what stops a silently wrong tautomer entering the docking.

Outputs (dock/ligands/):
    <species>.sdf        3D conformer(s), with charge/formula/species metadata
    <species>.pdbqt      Meeko output for Vina
    ligands.json         one record per species: SMILES, charge, heavy atoms, torsions, provenance
"""
import argparse
import json
import pathlib
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
ROOT = DOCK.parent
CCD = DOCK / 'ccd'
OUT = DOCK / 'ligands'

# species -> (SMILES, charge, CCD parent to check against, role)
# The CCD parent is the neutral component as deposited; "None" means there is no PDB component and the
# SMILES is the definition.  Charges are the formal charge of the SMILES and are asserted below.
SPECIES = {
    # --- the compound of the paper -------------------------------------------------------------
    'sfa_neutral': ('Nc1ccc(cc1)S(N)(=O)=O', 0, 'SAN',
                    'neutral sulfanilamide: the species present at pH 7.4 (primary ligand)'),
    'sfa_anion': ('Nc1ccc(cc1)S(=O)(=O)[NH-]', -1, 'SAN',
                  'sulfonamidate anion: the zinc-binding species, and the form used for DHPS/Sul2'),
    'sfa_cation': ('[NH3+]c1ccc(cc1)S(N)(=O)=O', +1, 'SAN',
                   '4-sulfamoylanilinium: the cation of the crystal, sensitivity case only'),
    # --- controls ------------------------------------------------------------------------------
    'smx_neutral': ('Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1', 0, '08D',
                    'sulfamethoxazole, neutral (pKa 5.7, minor at pH 7.4)'),
    'smx_anion': ('Cc1cc([N-]S(=O)(=O)c2ccc(N)cc2)no1', -1, '08D',
                  'sulfamethoxazole anion: the dominant species at pH 7.4'),
    'paba_anion': ('Nc1ccc(cc1)C(=O)[O-]', -1, 'PAB',
                   'p-aminobenzoate: the DHPS substrate, carboxylate at pH 7.4 (pKa 4.85)'),
    'azm_neutral': ('CC(=O)Nc1nnc(s1)S(N)(=O)=O', 0, 'AZM',
                    'acetazolamide, neutral (pKa 7.2, roughly half at pH 7.4)'),
    'azm_anion': ('CC(=O)Nc1nnc(s1)S(=O)(=O)[NH-]', -1, 'AZM',
                  'acetazolamide sulfonamidate: the zinc-binding species'),
    # --- co-crystal ligand of T5, needed for its redocking -------------------------------------
    'dhp_anion': (None, -1, '78H',
                  '7,8-dihydropteroate, carboxylate deprotonated: the Sul2 co-crystal product'),
    # --- control receptor C1 (DNA gyrase B, 4DUH), 2026-10-04 ----------------------------------
    'rli_anion': (None, -1, 'RLI',
                  "the 4,5'-bithiazole inhibitor of 4DUH, benzoate deprotonated: the co-crystal "
                  'ligand of the control receptor (redocking reference)'),
    'cfx_zwit': ('O=C([O-])c1cn(C2CC2)c2cc(N3CC[NH2+]CC3)c(F)cc2c1=O', 0, 'CPF',
                 'ciprofloxacin, zwitterion (carboxylate, protonated piperazine): a fluoroquinolone, '
                 'which acts at the GyrA-DNA cleavage complex and not at the GyrB ATP site; '
                 'comparator of the control receptor'),
}

# Which species is redocked against which target's co-crystal ligand.  This is the G2 gate matrix.
REDOCK = {
    'T1_3TZF': ['smx_anion', 'smx_neutral'],
    'T2_1AJ0': ['sfa_neutral'],
    'T3_6RL9': ['sfa_neutral', 'sfa_anion'],
    'T4_3HS4': ['azm_neutral', 'azm_anion'],
    'T5_7S2K': ['dhp_anion'],
    'C1_4DUH': ['rli_anion'],      # control receptor; not part of gate G2 (dock/control_gyrb.py)
}
CROSSDOCK = {'T3_6RL9': ['azm_neutral', 'azm_anion'], 'T4_3HS4': ['sfa_neutral', 'sfa_anion']}
# The control receptor has no cross-docking ligand: the clorobiocin pose of 1KZN does not fit 4DUH
# (dock/control_gyrb.py, cross_reference), so clorobiocin is not prepared or docked.


def ccd_parent(ccd_id):
    mol = Chem.MolFromMolFile(str(CCD / f'{ccd_id}_ideal.sdf'), removeHs=True)
    if mol is None:
        raise SystemExit(f'cannot read CCD template for {ccd_id}')
    return mol


def neutralise(mol):
    """Canonical SMILES of the neutral parent: every formal charge on N/O/S removed with the matching
    hydrogen count adjusted.  Used only for the consistency check against the CCD component."""
    rw = Chem.RWMol(Chem.RemoveHs(mol))
    for atom in rw.GetAtoms():
        q = atom.GetFormalCharge()
        if q:
            # an anion gains |q| hydrogens on neutralisation, a cation loses q of them
            atom.SetNumExplicitHs(max(0, atom.GetTotalNumHs() - q))
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(False)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def build_species(name, spec, report):
    smiles, charge, parent, role = spec
    if smiles is None:
        # defined by its chemical component (78H); deprotonate the carboxylic acid for pH 7.4
        mol = ccd_parent(parent)
        patt = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        match = mol.GetSubstructMatches(patt)
        rw = Chem.RWMol(mol)
        for _, _, o in match:
            rw.GetAtomWithIdx(o).SetFormalCharge(-1)
            rw.GetAtomWithIdx(o).SetNumExplicitHs(0)
            rw.GetAtomWithIdx(o).SetNoImplicit(True)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f'{name}: bad SMILES {smiles}')
    q = Chem.GetFormalCharge(mol)
    if q != charge:
        raise SystemExit(f'{name}: SMILES charge {q} does not match the declared charge {charge}')
    check = None
    if parent:
        # the neutral parent of the species must be the deposited component, otherwise the SMILES is
        # a different molecule (wrong tautomer, wrong substituent) and everything after it is wrong
        got, want = neutralise(mol), Chem.MolToSmiles(ccd_parent(parent))
        check = dict(ccd=parent, ccd_smiles=want, species_neutral_parent=got, match=got == want)
        if not check['match']:
            raise SystemExit(f'{name}: neutral parent {got} != CCD {parent} {want}')
    report[name] = dict(smiles=smiles, charge=charge, role=role, ccd_check=check,
                        formula=rdMolDescriptors.CalcMolFormula(mol),
                        mw=round(Descriptors.MolWt(mol), 2),
                        heavy_atoms=mol.GetNumHeavyAtoms(),
                        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol))
    return mol


def embed(mol, seed, n_conf=1):
    """ETKDGv3 conformer(s) plus MMFF94s minimisation.  MMFF is only cleaning up bond lengths and
    angles here; Vina keeps them rigid, so the geometry that matters is the one we hand it."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    if not len(cids):
        raise SystemExit('embedding failed')
    AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94s', maxIters=2000)
    return mol


def from_dft_xyz(path, smiles, name):
    """Take a DFT geometry and give it the species' bond orders and formal charges.

    The xyz carries elements and coordinates only.  rdDetermineBonds infers connectivity from the
    geometry; AssignBondOrdersFromTemplate then imposes the bond orders and charges of the species
    SMILES.  No atom-order convention is needed, so the DFT pipeline can write its xyz in whatever
    order it likes."""
    from rdkit.Chem import rdDetermineBonds
    raw = Chem.MolFromXYZFile(str(path))
    if raw is None:
        raise SystemExit(f'{name}: cannot read {path}')
    tpl = Chem.AddHs(Chem.MolFromSmiles(smiles))
    charge = Chem.GetFormalCharge(tpl)
    mol = Chem.Mol(raw)
    rdDetermineBonds.DetermineConnectivity(mol, charge=charge)
    mol = AllChem.AssignBondOrdersFromTemplate(tpl, mol)
    Chem.SanitizeMol(mol)
    if Chem.GetFormalCharge(mol) != charge:
        raise SystemExit(f'{name}: DFT geometry gave charge {Chem.GetFormalCharge(mol)} != {charge}')
    return mol


def input_mol_for_seed(species, seed, geometry='etkdg', xyz_override=None):
    """The 3D molecule (all H) a docking seed starts from, and a provenance record.

    Production 2026-09-28: until this helper existed both runners always rebuilt an ETKDG conformer
    from the species SMILES, so `prep_ligands.py --dft ...` wrote DFT geometries that no docking run
    ever read.  `geometry='dft'` now hands every seed the DFT geometry that prep_ligands wrote to
    dock/ligands/<species>.sdf, after checking in ligands.json that it really is the DFT one.

    Using the same input for every seed is correct, not a shortcut: Vina (and AD4Zn through Vina)
    randomises the position, orientation and every active torsion at the start of each Monte Carlo
    chain from --seed, and keeps bond lengths, bond angles and the geometry inside each rigid fragment
    exactly as given.  So the seed still decides the search, and the DFT file decides the internal
    geometry, which is the point of docking the DFT structure."""
    if geometry == 'etkdg':
        mol = embed(build_species(species, SPECIES[species], {}), seed=seed)
        return mol, dict(source='ETKDGv3+MMFF94s', seed=int(seed))
    if geometry != 'dft':
        raise SystemExit(f'unknown geometry {geometry!r}')
    import hashlib
    if xyz_override is not None:
        # a diagnostic geometry (e.g. the DFT structure with one angle changed): used as given, not
        # cross-checked against dock/ligands/, and recorded by its own sha256 (production 2026-09-28)
        xyz = pathlib.Path(xyz_override)
        mol = from_dft_xyz(xyz, SPECIES[species][0], species)
        return mol, dict(source='xyz_override', xyz=str(xyz),
                         xyz_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                         note='diagnostic geometry, same input for every seed')
    rec = json.loads((OUT / 'ligands.json').read_text())['species'].get(species, {})
    geo = rec.get('geometry', {})
    if geo.get('source') != 'DFT':
        raise SystemExit(f'{species}: ligands.json records geometry {geo}, not DFT; run '
                         f'prep_ligands.py --dft {species}=<opt.xyz> first')
    # rebuilt from the xyz itself (full precision; the SDF rounds to 1e-4 A) and cross-checked
    # against the SDF that prep_ligands wrote, so the two can never silently differ
    xyz = pathlib.Path(geo['file'])
    mol = from_dft_xyz(xyz, rec['smiles'], species)
    sdf = OUT / f'{species}.sdf'
    ref = Chem.MolFromMolFile(str(sdf), removeHs=False)
    want = Chem.MolToSmiles(Chem.MolFromSmiles(rec['smiles']))
    if ref is None or Chem.MolToSmiles(Chem.RemoveHs(mol)) != want \
            or ref.GetNumAtoms() != mol.GetNumAtoms():
        raise SystemExit(f'{species}: {xyz} / {sdf} is not {want}')
    dev = max((mol.GetConformer().GetAtomPosition(i) - ref.GetConformer().GetAtomPosition(i)).Length()
              for i in range(mol.GetNumAtoms()))
    if dev > 1e-3:
        raise SystemExit(f'{species}: {xyz} and {sdf} differ by {dev:.4f} A; rerun prep_ligands --dft')
    return mol, dict(source='DFT', xyz=str(xyz.relative_to(ROOT)) if xyz.is_relative_to(ROOT)
                     else str(xyz), xyz_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                     sdf=str(sdf.relative_to(DOCK)), note='same input for every seed; Vina '
                     'randomises position, orientation and torsions from --seed')


def write_pdbqt(mol, path, name):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    prep = MoleculePreparation()
    setups = prep(mol)
    txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise SystemExit(f'{name}: Meeko refused to write PDBQT: {err}')
    path.write_text(txt)
    return sum(1 for ln in txt.splitlines() if ln.startswith('BRANCH'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', nargs='*', default=None)
    ap.add_argument('--dft', nargs='*', default=[], metavar='SPECIES=PATH',
                    help='use a DFT xyz as the input geometry for this species')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dft = dict(kv.split('=', 1) for kv in a.dft)
    report = {}
    if a.only and (OUT / 'ligands.json').exists():
        # --only rebuilds the named species and keeps every other record (2026-10-04).  Before this,
        # ligands.json was rewritten with the named species alone, so adding one ligand would have
        # dropped the records the production analysis reads for all the others.
        report = json.loads((OUT / 'ligands.json').read_text())['species']
    for name, spec in SPECIES.items():
        if a.only and name not in a.only:
            continue
        mol = build_species(name, spec, report)
        smiles = report[name]['smiles']
        if name in dft:
            path = pathlib.Path(dft[name])
            path = path if path.is_absolute() else ROOT / path
            mol3d = from_dft_xyz(path, smiles, name)
            report[name]['geometry'] = dict(source='DFT', file=str(path))
        else:
            mol3d = embed(mol, seed=a.seed)
            report[name]['geometry'] = dict(source='ETKDGv3+MMFF94s', seed=a.seed)
        mol3d.SetProp('_Name', name)
        mol3d.SetProp('species', name)
        mol3d.SetProp('charge', str(report[name]['charge']))
        w = Chem.SDWriter(str(OUT / f'{name}.sdf'))
        w.write(mol3d)
        w.close()
        nbranch = write_pdbqt(mol3d, OUT / f'{name}.pdbqt', name)
        report[name]['active_torsions'] = nbranch
        print(f'{name:14s} q={report[name]["charge"]:+d} {report[name]["formula"]:14s} '
              f'heavy={report[name]["heavy_atoms"]:2d} torsions={nbranch} '
              f'{report[name]["geometry"]["source"]}')
    (OUT / 'ligands.json').write_text(json.dumps(
        dict(species=report, redock=REDOCK, crossdock=CROSSDOCK), indent=2) + '\n')
    print(f'\nwrote {OUT}/ligands.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
