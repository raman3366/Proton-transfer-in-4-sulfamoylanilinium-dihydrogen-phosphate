"""
Docking target definitions for the SAPH project, plus an inventory/verification pass over the
deposited mmCIFs.

Why this file exists: the template paper (R4) docked into boxes 75-90 A away from the real binding
sites and reported contact residues ~30-40 A from any ligand.  The only defence against repeating
that is to derive every box centre from the co-crystal ligand of the very file we prepare, and to
assert -- in code, not in prose -- which chain, which alternate location, which metals and which
cofactors survive receptor preparation.  `python dock/targets.py` writes dock/targets_inventory.json
and prints a table; every later script imports TARGETS from here instead of hard-coding numbers.

Facts below were checked against RCSB (data.rcsb.org entry + the deposited mmCIF) on 2026-09-27.

Naming: T1 3TZF (DHPS Y. pestis, primary), T2 1AJ0 (DHPS E. coli, secondary), T3 6RL9 (hCA II,
sulfanilamide), T4 3HS4 (hCA II, acetazolamide), T5 7S2K (Sul2, resistance contrast).

Control receptor (2026-10-04): C1 4DUH, the ATP site of E. coli DNA gyrase B, an enzyme that is not a
target of sulfanilamide.  It is named C1, not T6, on purpose: the validation and production analyses
collect their jobs with the glob 'T*/runs/...', so a control directory can never enter gate G2, the
production summary or any statistic of the five targets.  It is run and analysed by
dock/control_gyrb.py only.  X_1KZN supplies the clorobiocin pose, which shows that the control has no
cross-docking test (the pose does not fit 4DUH; clorobiocin is not docked).  Both
entries were checked against RCSB (entry + deposited mmCIF) on 2026-10-04.
"""
import gzip
import json
import pathlib
import sys

import gemmi
import numpy as np

DOCK = pathlib.Path(__file__).resolve().parent
PDBDIR = DOCK / 'pdb'

# Box edge, angstrom.  20 A cubic as decided in PLAN.md section 3 (P3).  Vina's own guidance is that
# the box must enclose the ligand with room for translation; 20 A is generous for 11-17 heavy atoms
# and keeps the search space small enough for exhaustiveness 32 to be meaningful.
BOX_SIZE = 20.0

# Residue names that are, by default, not part of the receptor: bulk solvent and cryo/crystallisation
# additives.  BE7 (p-chloromercuribenzoate) and HG in 6RL9 sit on Cys206, 12 A from the zinc site, and
# are a crystallisation artefact of that particular experiment -- they are not part of hCA II.
#
# Precedence (verification 2026-09-28): a target's `keep` list always overrides JUNK.  Receptor
# selection (prep_receptor.select) keeps exactly the `keep` names and nothing else, so JUNK is only
# consulted when the co-crystal ligand is resolved automatically (T5) and must not pick an additive.
# SO4 is therefore junk in general but deliberately kept in T2_1AJ0, where it occupies the
# pyrophosphate subsite; `junk_for(spec)` is the single place that applies this rule.
JUNK = {
    'HOH', 'DOD', 'GOL', 'EDO', 'PEG', 'PG4', 'MPD', 'DMS', 'SO4', 'PO4', 'ACT', 'ACY', 'CL', 'NA',
    'K', 'TRS', 'EPE', 'MES', 'IMD', 'CIT', 'FMT', 'NO3', 'BME', 'BE7', 'HG', 'AGC', 'GLC', 'BGC',
}

METALS = {'ZN', 'MG', 'MN', 'CA', 'FE', 'NI', 'CO', 'CU', 'NA', 'K'}


def junk_for(spec):
    """JUNK minus whatever this target explicitly keeps (keep always wins; see the note above)."""
    return JUNK - set(spec.get('keep', []))

TARGETS = {
    # --- DHPS, primary -------------------------------------------------------------------------
    'T1_3TZF': dict(
        pdb='3TZF',
        title='DHPS, Yersinia pestis, sulfamethoxazole complex',
        chain='B',                 # SMX (08D) is modelled in chain B only, occupancy 0.82
        ligand=('08D', 'B'),       # sulfamethoxazole = the redocking reference
        keep=['HH2', 'MG'],        # HH2 = 6-hydroxymethylpterin monophosphate/pterin-PP; Mg bridges it
        role='primary',
        # Chain B is missing residues 66-67 (a loop-1 gap next to the pABA pocket).  PDBFixer can
        # rebuild a 2-residue gap reliably, so we do rebuild it, but the run is repeated with the gap
        # left open as a sensitivity check (see prep_receptor.py --no-rebuild).
        rebuild_loops=True,
        notes='Closed, catalytically competent pABA pocket (Yun 2012, 10.1126/science.1214641).',
    ),
    # --- DHPS, secondary ----------------------------------------------------------------------
    'T2_1AJ0': dict(
        pdb='1AJ0',
        title='DHPS, Escherichia coli, sulfanilamide + pterin ternary complex',
        chain='A',
        ligand=('SAN', 'A'),
        # PH2 = 6-hydroxymethylpterin.  This entry has no Mg and no pyrophosphate; a sulfate sits in
        # the pyrophosphate subsite and touches the bound sulfanilamide (2.85 A, measured).  Removing
        # it would leave a hole that a docked ligand could occupy for the wrong reason, so the primary
        # receptor keeps it and `prep_receptor.py T2_1AJ0 --drop SO4 --suffix _noSO4` is the
        # sensitivity case.  Both are redocked and both numbers are reported.
        keep=['PH2', 'SO4'],
        role='secondary',
        rebuild_loops=True,
        notes=('Open-loop state: loops 1 and 2 are ordered but Mg/pyrophosphate are absent, so the '
               'sulfonamide sits ~2.7 A (ring) / ~5.5 A (S) from where it sits in the closed pocket. '
               'Reported with that caveat; not used to argue binding mode.'),
    ),
    # --- hCA II -------------------------------------------------------------------------------
    'T3_6RL9': dict(
        pdb='6RL9',
        title='human carbonic anhydrase II, sulfanilamide complex, 1.0 A',
        chain='A',
        ligand=('SAN', 'A', 304),  # three SAN copies; only A304 is the zinc-bound one
        keep=['ZN'],
        role='zinc',
        rebuild_loops=True,
        # Zinc-site tautomers, verified below from the Zn coordination distances: His94 and His96
        # bind through NE2 (so the proton sits on ND1 -> HID), His119 binds through ND1 (proton on
        # NE2 -> HIE).
        his_templates={94: 'HID', 96: 'HID', 119: 'HIE'},
        # His64 (proton shuttle, ~9-10 A from the crystal ligand, ~11 A from Zn) is neutral by PROPKA
        # (pKa 3.7) but its tautomer is NOT determined by anything in the structure: it is an
        # explicit assumption, HIE.  PDB2PQR's H-bond optimiser picks HIE here and HID in 3HS4; the
        # Vina sensitivity test (HID64 receptors) is reported in dock/VALIDATION.md.
        his_assumed={64: 'HIE'},
        default_altloc='A',        # His64 in/out conformers; take A
        notes='Glockner 2020, 10.3390/biom10040509.  BE7 + Hg on Cys206 are an additive: dropped.',
    ),
    'T4_3HS4': dict(
        pdb='3HS4',
        title='human carbonic anhydrase II, acetazolamide complex, 1.1 A',
        chain='A',
        ligand=('AZM', 'A', 701),
        keep=['ZN'],
        role='zinc',
        rebuild_loops=True,
        his_templates={94: 'HID', 96: 'HID', 119: 'HIE'},
        his_assumed={64: 'HIE'},   # assumption, as in T3 (PDB2PQR would pick HID here)
        default_altloc='A',
        notes='Sippel 2009, 10.1107/S1744309109036665.  Cross-docking partner of T3.',
    ),
    # --- resistance contrast ------------------------------------------------------------------
    'T5_7S2K': dict(
        pdb='7S2K',
        title='Sul2, plasmid-borne sulfonamide-resistant dihydropteroate synthase',
        chain='A',
        ligand=None,               # resolved by the inventory pass: the dihydropteroate product
        keep=['MG', 'POP'],
        role='resistance',
        rebuild_loops=True,
        notes=('Venkatesan 2023, 10.1038/s41467-023-39778-7.  Product complex: 7,8-dihydropteroate '
               '+ pyrophosphate + Mg.'),
    ),
    # --- substrate control donor (not docked into; pABA is lifted out of it) -------------------
    'X_3TYZ': dict(
        pdb='3TYZ',
        title='DHPS, Y. pestis, pABA complex -- source of the pABA control ligand only',
        chain='A',
        ligand=('PAB', 'A'),
        keep=[],
        role='ligand_source',
        rebuild_loops=False,
        notes='Used only to extract the crystallographic pABA geometry.',
    ),
    # --- control receptor: not a target of sulfanilamide (see the module docstring) ---------------
    'C1_4DUH': dict(
        pdb='4DUH',
        title='DNA gyrase B, Escherichia coli, 24 kDa N-terminal (ATPase) fragment, bithiazole inhibitor',
        chain='A',                 # RLI is modelled in both chains; B-factor 11.4 (A) against 23.2 (B)
        ligand=('RLI', 'A', 301),  # the 4,5'-bithiazole inhibitor = the redocking reference
        keep=[],                   # no metal and no cofactor in this fragment
        role='control',
        # Chain A lacks residues 109-116 of the ATP lid (chain B lacks 105-116); the residues that
        # flank the gap lie 9.9-13.6 A from the inhibitor.  PDBFixer rebuilds the eight residues, as it
        # rebuilds the internal gaps of the five targets, but a loop built without density is a guess,
        # so the receptor is prepared a second time from another PDBFixer seed
        # (`prep_receptor.py C1_4DUH --seed 12345 --suffix _seed12345`, the T1_3TZF_seed12345
        # pattern) and both are docked; control_gyrb.py reports how far the built loop lies from
        # every pose and how much the two builds differ.
        rebuild_loops=True,
        notes=('Brvar 2012, 10.1021/jm300395d.  ATP site of GyrB: the site of the aminocoumarins, '
               'not of the sulfonamides (whose target is DHPS) and not of the fluoroquinolones '
               '(which bind the GyrA-DNA cleavage complex).  Negative-control receptor.'),
    ),
    # --- clorobiocin pose for C1: shows that no cross-docking test is possible (not docked into) ---
    'X_1KZN': dict(
        pdb='1KZN',
        title='DNA gyrase B, E. coli, 24 kDa fragment, clorobiocin complex -- source of the '
              'clorobiocin pose only',
        chain='A',
        ligand=('CBN', 'A'),
        keep=[],
        role='ligand_source',
        rebuild_loops=False,
        # The author numbering of 1KZN is not the E. coli numbering between residues 82 and 96
        # (twelve of fifteen residue names differ from 4DUH at the same number), so its pose is
        # carried into the 4DUH frame by a sequence-aligned Ca superposition
        # (control_gyrb.cross_reference), never by residue number.
        notes='Lafitte 2002, 10.1021/bi0159837.  Used only for the clorobiocin crystal pose.',
    ),
}


def read_structure(pdb_id):
    """Read a gzipped mmCIF from dock/pdb/.  Kept gzipped: the repo stores the deposited file as
    downloaded so the sha256 in the data package matches RCSB."""
    path = PDBDIR / f'{pdb_id}.cif.gz'
    if not path.exists():
        raise FileNotFoundError(f'{path} missing; fetch with dock/fetch_pdb.sh')
    with gzip.open(path, 'rt') as fh:
        doc = gemmi.cif.read_string(fh.read())
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    return st


def het_inventory(st):
    """Every non-water, non-polymer residue with its chain, seqid, atom count and mean occupancy."""
    out = []
    for model in st:
        for chain in model:
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                is_poly = bool(info and info.is_standard() and info.is_amino_acid())
                if is_poly or res.name in ('HOH', 'DOD'):
                    continue
                occ = float(np.mean([a.occ for a in res])) if len(res) else 0.0
                bfac = float(np.mean([a.b_iso for a in res])) if len(res) else 0.0
                heavy = [a for a in res if a.element != gemmi.Element('H')]
                out.append(dict(name=res.name, chain=chain.name, seqid=res.seqid.num,
                                natoms=len(res), heavy=len(heavy), occ=round(occ, 3),
                                b_iso=round(bfac, 1),
                                altlocs=sorted({a.altloc for a in res if a.altloc})))
        break
    return out


def centroid(st, resname, chain_name, seqid=None):
    """Heavy-atom centroid and extent of one ligand copy -- the box centre."""
    for model in st:
        for chain in model:
            if chain.name != chain_name:
                continue
            for res in chain:
                if res.name != resname:
                    continue
                if seqid is not None and res.seqid.num != seqid:
                    continue
                pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res
                                if a.element != gemmi.Element('H')])
                return pos.mean(axis=0), pos.max(axis=0) - pos.min(axis=0), len(pos), res.seqid.num
        break
    raise KeyError(f'{resname} not found in chain {chain_name} (seqid={seqid})')


def metal_coordination(st, chain_name, cutoff=2.6):
    """Distances from every metal to nearby N/O/S donors.  Used to fix His tautomers at the zinc
    site: the nitrogen that coordinates carries no proton, so His94/96 binding via NE2 are HID and
    His119 binding via ND1 is HIE."""
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    rows = []
    model = st[0]
    for chain in model:
        for res in chain:
            if res.name not in METALS or len(res) != 1:
                continue
            atom = res[0]
            marks = ns.find_atoms(atom.pos, '\0', radius=cutoff)
            partners = []
            for m in marks:
                cra = m.to_cra(model)
                if cra.residue.name == res.name and cra.residue.seqid.num == res.seqid.num:
                    continue
                if cra.atom.element in (gemmi.Element('C'), gemmi.Element('H')):
                    continue
                d = atom.pos.dist(cra.atom.pos)
                if d <= cutoff:
                    partners.append((f'{cra.chain.name}/{cra.residue.name}{cra.residue.seqid.num}'
                                     f'/{cra.atom.name}', round(d, 2)))
            partners.sort(key=lambda t: t[1])
            rows.append(dict(metal=res.name, chain=chain.name, seqid=res.seqid.num,
                             partners=partners))
    return rows


def deposited_sequence(pdb_id, auth_chain):
    """The crystallised sequence of one author chain, and which of its residues are unobserved, from
    _pdbx_poly_seq_scheme in the deposited mmCIF.

    Not from gemmi's entities: read_structure() calls setup_entities(), which rebuilds the entity list
    from the coordinates and therefore reports only the residues that were *modelled* -- which is
    precisely the information that must not be used to decide what is missing.  This lives here rather
    than in prep_receptor.py so that the inventory and the receptor preparation answer the question
    "which residues are missing" from one function (verification 2026-09-28).
    """
    with gzip.open(PDBDIR / f'{pdb_id}.cif.gz', 'rt') as fh:
        block = gemmi.cif.read_string(fh.read()).sole_block()
    table = block.find('_pdbx_poly_seq_scheme.',
                       ['pdb_strand_id', 'mon_id', 'pdb_seq_num', 'auth_seq_num'])
    seq, unobserved = [], []
    for row in table:
        if row[0] != auth_chain:
            continue
        seq.append(row[1])
        if row[3] in ('?', '.'):
            unobserved.append(row[2])
    return seq, unobserved


def polymer_gaps(st, chain_name):
    """Residue numbers absent from a chain's modelled numbering range.

    Careful: this is *numbering*, not chemistry.  hCA II is numbered without a residue 126 at all
    (_pdbx_poly_seq_scheme goes ... 125, 127 ...), so 6RL9 and 3HS4 show 126 here while their chains
    are chemically continuous.  The list of residues that really are missing is the `unobserved`
    field, which comes from deposited_sequence(); prep_receptor.classify_gaps() uses that one to decide
    what PDBFixer rebuilds."""
    nums = []
    for chain in st[0]:
        if chain.name != chain_name:
            continue
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                nums.append(res.seqid.num)
    if not nums:
        return dict(first=None, last=None, missing=[])
    lo, hi = min(nums), max(nums)
    have = set(nums)
    return dict(first=lo, last=hi, missing=[n for n in range(lo, hi + 1) if n not in have])


def inventory():
    report = {}
    for key, spec in TARGETS.items():
        st = read_structure(spec['pdb'])
        het = het_inventory(st)
        lig = spec['ligand']
        overridden = sorted(JUNK & set(spec.get('keep', [])))
        seq, unobserved = deposited_sequence(spec['pdb'], spec['chain'])
        gaps = polymer_gaps(st, spec['chain'])
        gaps['unobserved'] = unobserved
        gaps['note'] = ('"missing" is numbering-derived and can list a number the entry simply does '
                        'not use (hCA II has no residue 126); "unobserved" is from '
                        '_pdbx_poly_seq_scheme and is the real list of residues not modelled.')
        entry = dict(pdb=spec['pdb'], title=spec['title'], chain=spec['chain'], role=spec['role'],
                     junk_names_kept_by_override=overridden,
                     spacegroup=st.spacegroup_hm, resolution=st.resolution,
                     deposited_sequence_length=len(seq),
                     het=het, gaps=gaps,
                     metals=metal_coordination(st, spec['chain']), notes=spec['notes'])
        if lig is None:
            # T5: the co-crystal ligand is whatever large organic HET sits in chain A and is not in
            # `keep`.  Resolved here rather than hard-coded so the choice is visible in the JSON.
            cands = [h for h in het if h['chain'] == spec['chain'] and h['name'] not in spec['keep']
                     and h['name'] not in junk_for(spec) and h['heavy'] >= 8]
            cands.sort(key=lambda h: -h['heavy'])
            if not cands:
                raise SystemExit(f'{key}: no co-crystal ligand candidate found')
            lig = (cands[0]['name'], spec['chain'], cands[0]['seqid'])
            entry['ligand_resolved_automatically'] = True
        c, ext, nheavy, seqid = centroid(st, lig[0], lig[1], lig[2] if len(lig) > 2 else None)
        entry['ligand'] = dict(name=lig[0], chain=lig[1], seqid=seqid, heavy=nheavy,
                               centre=[round(float(v), 3) for v in c],
                               extent=[round(float(v), 2) for v in ext],
                               fits_box=bool((ext + 6.0 < BOX_SIZE).all()))
        entry['box'] = dict(center=[round(float(v), 3) for v in c],
                            size=[BOX_SIZE, BOX_SIZE, BOX_SIZE])
        report[key] = entry
    return report


def main():
    rep = inventory()
    out = DOCK / 'targets_inventory.json'
    out.write_text(json.dumps(rep, indent=2) + '\n')
    for key, e in rep.items():
        lg = e['ligand']
        print(f"{key:10s} {e['pdb']}  {e['spacegroup']:12s} {e['resolution']:.2f} A  "
              f"chain {e['chain']}  lig {lg['name']}{lg['seqid']} ({lg['heavy']} heavy, "
              f"extent {lg['extent']}, fits={lg['fits_box']})")
        print(f"{'':12s}centre {lg['centre']}  unobserved {e['gaps']['unobserved'] or 'none'} "
              f"(numbering gaps {e['gaps']['missing'] or 'none'}, range "
              f"{e['gaps']['first']}-{e['gaps']['last']})")
        for m in e['metals']:
            print(f"{'':12s}{m['metal']}{m['seqid']}: {m['partners']}")
        hets = sorted({h['name'] for h in e['het']} - {'HOH'})
        print(f"{'':12s}HET present: {hets}")
    print(f'\nwrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
