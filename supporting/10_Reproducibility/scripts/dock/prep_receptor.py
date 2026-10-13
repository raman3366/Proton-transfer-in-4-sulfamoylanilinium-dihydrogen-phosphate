"""
Receptor preparation for the SAPH docking work: one target per run, everything written to
dock/<TARGET>/.

    python dock/prep_receptor.py T3_6RL9
    python dock/prep_receptor.py T1_3TZF --no-rebuild        # sensitivity: leave the 66-67 gap open
    python dock/prep_receptor.py T2_1AJ0 --drop SO4          # sensitivity: remove the sulfate

Why the steps are in this order:

1.  *Selection before anything else.*  We take one chain, one altloc and exactly the HET groups the
    target definition lists.  Everything else -- extra copies of the ligand, sugars, the
    p-chloromercuribenzoate/Hg pair in 6RL9, glycerol, chloride -- is crystallisation chemistry, not
    the protein, and leaving it in would both block the site and inflate the contact lists.

2.  *Geometry repair with PDBFixer, restricted to internal gaps.*  PDBFixer will happily extend a
    chain past its last modelled residue; those atoms are pure guesswork and they sit on the surface,
    so we filter the missing-residue list down to internal gaps (and offer --no-rebuild to test the
    effect).  Missing side-chain atoms are always rebuilt: an incomplete Lys or Arg next to the
    pocket changes the score for no physical reason.

3.  *Protonation decided by PROPKA, applied by OpenMM.*  PDBFixer's own addMissingHydrogens(pH) uses
    fixed rules, so it cannot know that a buried Asp is protonated.  We therefore run PROPKA 3.5 on
    the repaired heavy-atom structure, turn the predicted pKa values into an explicit OpenMM variant
    list (ASH/GLH/HIP/HID/HIE/LYN/CYX), and let Modeller.addHydrogens place the hydrogens.  The whole
    list is written to the report, so the choice for every titratable residue is auditable.

4.  *Metal-site tautomers forced by hand.*  A histidine that coordinates a metal through NE2 carries
    its proton on ND1 (= HID in Amber naming) and vice versa.  PROPKA does not know about the zinc,
    so for hCA II we override His94/His96 to HID and His119 to HIE.  targets.py verified from the
    deposited coordinates that this is indeed the coordination pattern (Zn-NE2 2.0 A for 94 and 96,
    Zn-ND1 2.0 A for 119), so the override is derived from the structure, not from memory.

5.  *Waters: only the first shell of a kept metal survives.*  A docking receptor with crystal waters
    in the pocket scores the wrong thing, because the incoming ligand would displace them.  Water
    bound directly to Mg2+ is different: it never leaves the metal on the timescale of binding and it
    completes the octahedron that orients the pyrophosphate.  The rule implemented here is therefore
    "keep a water only if it is within 2.6 A of a metal we keep", and the kept waters are listed in
    the report.  In hCA II the zinc water/hydroxide is displaced by the sulfonamide nitrogen, so
    nothing is kept there.

6.  *PDBQT in two pieces.*  Meeko 0.8 parameterises the polymer from its residue templates, but it
    has no templates for HH2, PH2, POP or 78H.  Those cofactors are therefore typed through Meeko's
    *ligand* path (RDKit molecule with bond orders taken from the RCSB chemical-component template,
    hydrogens added at the crystal coordinates) and their ATOM records are appended to the rigid
    receptor PDBQT.  Vina's rigid receptor is just a list of typed atoms, so concatenation is exact;
    metals are written directly with their AutoDock type and formal charge.

Outputs in dock/<TARGET>/:
    raw_selected.pdb      chain/altloc selection, heavy atoms, before repair
    fixed_protein.pdb     after PDBFixer + explicit hydrogens
    propka.pka            PROPKA output
    receptor.pdb          protein + kept cofactors/metals/waters, with hydrogens
    receptor.pdbqt        rigid receptor for Vina / AutoDock4
    cofactors.pdbqt       the appended cofactor/metal/water block on its own (for inspection)
    ref_<LIG>.sdf         co-crystal ligand, crystal coordinates, bond orders from the CCD
    box.json              box centre (co-crystal ligand centroid) and size
    prep_report.json      every decision, with counts and the checks that were run
"""
import argparse
import json
import pathlib
import shutil
import subprocess
import sys

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from targets import (TARGETS, PDBDIR, METALS, BOX_SIZE, read_structure,  # noqa: E402
                     deposited_sequence)

DOCK = pathlib.Path(__file__).resolve().parent
CCD = DOCK / 'ccd'
PH = 7.4
METAL_WATER_CUTOFF = 2.6      # A; first coordination shell of Mg2+/Zn2+

# AutoDock atom types for the metals we keep.  Vina 1.2 knows these types; its scoring function has
# no electrostatics, so the charge matters only for AutoDock4/AD4Zn maps.
METAL_AD_TYPE = {'ZN': ('Zn', 2.0), 'MG': ('Mg', 2.0), 'MN': ('Mn', 2.0), 'CA': ('Ca', 2.0),
                 'FE': ('Fe', 3.0), 'NA': ('Na', 1.0), 'K': ('K', 1.0)}


# --------------------------------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------------------------------
def pick_altloc(res, prefer='A'):
    """Choose one alternate conformation per residue and return the atoms to keep.

    gemmi's remove_alternative_conformations() keeps the first altloc of every *atom*, which is wrong
    when a residue models only part of its side chain in one conformer: in 6RL9 (1.0 A, ~90 residues
    with altlocs) that produced a Lys76 with a 2.42 A CG-CD "bond" and Meeko then failed to match any
    lysine template.  Selecting per residue -- the altloc that carries the most atoms, ties broken by
    occupancy and then by `prefer` -- keeps each side chain internally consistent.  Atoms that exist
    only in the rejected conformer are dropped and rebuilt by PDBFixer.
    """
    alts = sorted({a.altloc for a in res if a.altloc not in ('', '\0')})
    if not alts:
        return list(res), None
    best, score = None, None
    for L in alts:
        keep = [a for a in res if a.altloc in ('', '\0', L)]
        s = (len(keep), sum(a.occ for a in keep), L == prefer)
        if score is None or s > score:
            best, score = L, s
    return [a for a in res if a.altloc in ('', '\0', best)], best


def select(spec, drop=()):
    """Split the deposited structure into protein / reference ligand / kept HET / kept waters."""
    st = read_structure(spec['pdb'])
    st.setup_entities()
    st.remove_hydrogens()
    prefer = spec.get('default_altloc', 'A')
    model = st[0]
    chain_name = spec['chain']
    lig = spec['ligand']
    if lig is None:
        # targets.py resolved the co-crystal ligand from the deposited file (T5: the dihydropteroate
        # product).  Read its answer rather than duplicating the rule here.
        inv = json.loads((DOCK / 'targets_inventory.json').read_text())
        key = next(k for k, v in TARGETS.items() if v is spec)
        L = inv[key]['ligand']
        lig = (L['name'], L['chain'], L['seqid'])
    keep_names = [n for n in spec['keep'] if n not in drop]

    protein, cofactors, metals, waters, refres = [], [], [], [], None
    altloc_choices = {}
    for chain in model:
        for res in chain:
            atoms, chosen = pick_altloc(res, prefer)
            if chosen:
                altloc_choices[f'{chain.name}:{res.name}{res.seqid.num}'] = chosen
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                if chain.name == chain_name:
                    protein.append((chain.name, res, atoms))
                continue
            if res.name in ('HOH', 'DOD'):
                if chain.name == chain_name:
                    waters.append((chain.name, res, atoms))
                continue
            # the reference ligand: exact chain (+ seqid when the target pins one)
            if (res.name == lig[0] and chain.name == lig[1]
                    and (len(lig) < 3 or res.seqid.num == lig[2]) and refres is None):
                refres = (chain.name, res, atoms)
                continue
            if res.name in keep_names and chain.name == chain_name:
                (metals if res.name in METALS else cofactors).append((chain.name, res, atoms))
    if refres is None:
        raise SystemExit(f'reference ligand {lig} not found')

    # water triage: first shell of a kept metal only
    mpos = [a.pos for _, _, atoms in metals for a in atoms]
    kept_waters = []
    for cn, res, atoms in waters:
        o = atoms[0].pos
        if any(o.dist(mp) <= METAL_WATER_CUTOFF for mp in mpos):
            kept_waters.append((cn, res, atoms, round(min(o.dist(mp) for mp in mpos), 2)))

    seq, unobs = deposited_sequence(spec['pdb'], chain_name)
    return st, protein, refres, cofactors, metals, kept_waters, altloc_choices, seq, unobs


def seqres_records(chain_id, sequence):
    """SEQRES records for the selected chain.

    PDBFixer finds missing *residues* by aligning the deposited sequence against what is modelled, so
    a PDB file without SEQRES has, as far as PDBFixer is concerned, no gaps at all.  Leaving them out
    silently skipped the 66-67 gap in 3TZF chain B, which is exactly the loop next to the pABA pocket.
    The sequence comes from the mmCIF entity (entity_poly_seq), i.e. the construct that was
    crystallised, tag included; terminal extensions are filtered later."""
    out = []
    for i in range(0, len(sequence), 13):
        chunk = sequence[i:i + 13]
        names = ''.join(f' {n:>3s}' for n in chunk)
        out.append(f'SEQRES {i // 13 + 1:3d} {chain_id:1s} {len(sequence):4d} {names}')
    return out


def write_pdb(path, groups, hetatm=False, seqres=None):
    """Minimal PDB writer.  Used instead of gemmi's so that the cofactor/metal block keeps the exact
    residue names and numbering of the deposited file (PDBFixer never sees this block)."""
    lines, serial = list(seqres or []), 1
    for cn, res, atoms in groups:
        for a in atoms:
            rec = 'HETATM' if (hetatm or not (gemmi.find_tabulated_residue(res.name)
                                              and gemmi.find_tabulated_residue(res.name).is_amino_acid())) else 'ATOM  '
            name = a.name if len(a.name) >= 4 or len(a.element.name) == 2 else f' {a.name}'
            lines.append(f'{rec}{serial:5d} {name:<4s}{"":1s}{res.name:>3s} {cn:1s}'
                         f'{res.seqid.num:4d}{"":4s}{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}'
                         f'{a.occ:6.2f}{a.b_iso:6.2f}{"":10s}{a.element.name.upper():>2s}')
            serial += 1
    lines.append('END')
    path.write_text('\n'.join(lines) + '\n')


# --------------------------------------------------------------------------------------------------
# repair + protonation
# --------------------------------------------------------------------------------------------------
def classify_gaps(unobserved, modelled):
    """Split the unobserved residues of the selected chain into internal gaps and terminal extensions.

    This is the *single* definition of "internal gap" in this script.  Before the 2026-09-28
    verification there were two: PDBFixer's insertion-index test inside run_pdbfixer (which decides
    what is rebuilt) and a separate author-numbering test further down (which decides what the report
    says and whether --blunt_ends is needed).  They could disagree, and for T1_3TZF the report field a
    reader would check said nothing while the rebuild list held the 66-67 loop.  Both callers now use
    this function and the report states the two views side by side.

    `unobserved` is the author numbering from _pdbx_poly_seq_scheme; `modelled` the author numbers
    actually present in the selection.
    """
    nums = sorted(int(u) for u in unobserved if str(u).lstrip('-').isdigit())
    lo, hi = min(modelled), max(modelled)
    return dict(internal=[n for n in nums if lo < n < hi],
                terminal=[n for n in nums if n < lo or n > hi],
                modelled_range=[lo, hi])


def run_pdbfixer(in_pdb, out_pdb, rebuild_loops, report, seed=20260928):
    """PDBFixer places rebuilt atoms from randomly perturbed template positions, so two runs give
    different loops (for the 3TZF Ala66-Ala67 loop the raw energy before minimisation was 9.85e6
    kJ/mol with seed 20260928 and 6.79e4 kJ/mol with seed 12345 -- both values are in the two
    prep_report.json files).  Seeding both RNGs makes the receptor reproducible, which matters
    because every docking number downstream depends on it."""
    import random
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))
    report['pdbfixer_seed'] = seed
    fixer = PDBFixer(filename=str(in_pdb))
    fixer.findMissingResidues()
    # Restrict to internal gaps.  PDBFixer keys missingResidues by (chain index, insertion index);
    # an insertion index of 0 or len(residues) is a terminal extension, which we never build.
    chains = list(fixer.topology.chains())
    terminal, internal, rebuilt_auth = {}, {}, []
    for (ci, ri), names in fixer.missingResidues.items():
        residues = list(chains[ci].residues())
        n = len(residues)
        is_terminal = ri in (0, n)
        (terminal if is_terminal else internal)[f'{ci}:{ri}'] = names
        if not is_terminal:
            # PDBFixer inserts the stretch *before* modelled residue index ri, so in author numbering
            # the gap runs from (id of residue ri-1)+1 upwards.  Recorded so that the rebuild list and
            # the author-numbering gap list in the report can be compared directly.
            first = int(residues[ri - 1].id) + 1
            rebuilt_auth.append(dict(key=f'{ci}:{ri}', names=names,
                                     auth_numbers=list(range(first, first + len(names))),
                                     between=[int(residues[ri - 1].id), int(residues[ri].id)]))
    fixer.missingResidues = {} if not rebuild_loops else {
        k: v for k, v in fixer.missingResidues.items()
        if f'{k[0]}:{k[1]}' in internal}
    report['missing_residues_terminal_skipped'] = terminal
    report['missing_residues_rebuilt'] = internal if rebuild_loops else {}
    report['missing_residues_rebuilt_auth_numbering'] = rebuilt_auth if rebuild_loops else []
    fixer.findNonstandardResidues()
    report['nonstandard_residues'] = [[str(a), b] for a, b in fixer.nonstandardResidues]
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    report['missing_atoms'] = {str(k): [str(a) for a in v] for k, v in fixer.missingAtoms.items()}
    fixer.addMissingAtoms()
    with open(out_pdb, 'w') as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return fixer


def run_propka(pdb, workdir, report):
    """PROPKA 3.5 on the repaired heavy-atom protein.  Returns {(chain, resnum, resname): pKa}."""
    out = workdir / 'propka.pka'
    cmd = ['propka3', '-o', str(PH), '--log-level', 'ERROR', str(pdb)]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    produced = workdir / (pdb.stem + '.pka')
    if produced.exists() and produced != out:
        shutil.move(produced, out)
    if not out.exists():
        report['propka'] = dict(ok=False, stderr=r.stderr[-2000:])
        return {}
    pkas = {}
    in_block = False
    for line in out.read_text().splitlines():
        if line.startswith('SUMMARY OF THIS PREDICTION'):
            in_block = True
            continue
        if in_block:
            if line.startswith('-') or not line.strip():
                if pkas:
                    break
                continue
            f = line.split()
            if len(f) >= 4 and f[0].isalpha() and f[1].lstrip('-').isdigit():
                try:
                    pkas[(f[2], int(f[1]), f[0])] = float(f[3])
                except ValueError:
                    pass
    report['propka'] = dict(ok=True, n_residues=len(pkas), file='propka.pka')
    return pkas


def disulfides(topology, positions):
    """CYS pairs with SG-SG below 2.5 A.  Needed because we pass an explicit variant list to
    Modeller, which then no longer detects disulfides itself."""
    import openmm.unit as unit
    pos = np.array(positions.value_in_unit(unit.angstrom))
    sg = [(res, next((a.index for a in res.atoms() if a.name == 'SG'), None))
          for res in topology.residues() if res.name in ('CYS', 'CYX')]
    sg = [(r, i) for r, i in sg if i is not None]
    bonded = set()
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if np.linalg.norm(pos[sg[i][1]] - pos[sg[j][1]]) < 2.5:
                bonded.add(sg[i][0].index)
                bonded.add(sg[j][0].index)
    return bonded


def variant_list(topology, positions, pkas, his_override, report, his_assumed=None):
    """Turn PROPKA pKa values into the OpenMM Modeller variant list.

    Rules (pH 7.4): an acid whose predicted pKa is above the pH stays protonated; a Lys whose pKa is
    below the pH is neutral; a His above the pH is doubly protonated, otherwise it is neutral and we
    take the HIE tautomer unless the target definition overrides it (metal sites).  Tyr and Arg have
    no OpenMM variant and are left in their standard states; any anomalous prediction for them is
    recorded so it can be discussed rather than silently ignored.

    `his_assumed` (targets.py) names neutral histidines whose tautomer the coordinates do not
    determine but which matter enough to be declared: they get the stated tautomer *and* an entry in
    report['protonation_assumed'], so that nothing site-relevant is decided by a silent default."""
    ss = disulfides(topology, positions)
    variants, decisions, flagged, flagged_assumed = [], [], [], []
    for res in topology.residues():
        num, ch = int(res.id), res.chain.id
        pka = pkas.get((ch, num, res.name))
        v = None
        if res.name == 'HIS':
            if num in his_override:
                v = his_override[num]
                why = 'metal-site override (coordination geometry)'
            elif pka is not None and pka > PH:
                v, why = 'HIP', f'PROPKA pKa {pka:.2f} > {PH}'
            elif his_assumed and num in his_assumed:
                # A neutral histidine whose tautomer nothing in the structure determines, but which is
                # close enough to the site to be worth naming: the choice is an assumption and is
                # labelled as one, and its effect is measured (dock/VALIDATION.md).
                v = his_assumed[num]
                why = (f'PROPKA pKa {pka if pka is None else round(pka, 2)}; neutral; tautomer {v} is '
                       'an ASSUMPTION (targets.py his_assumed), not derived from the coordinates')
                flagged_assumed.append(dict(residue=f'{ch}:{res.name}{num}', variant=v,
                                            pka=None if pka is None else round(pka, 2),
                                            note='tautomer assumed, not derived; sensitivity tested '
                                                 'in dock/VALIDATION.md'))
            else:
                v, why = 'HIE', (f'PROPKA pKa {pka if pka is None else round(pka, 2)}; neutral, HIE '
                                 '(default tautomer; not derived from the coordinates)')
        elif res.name in ('ASP', 'GLU') and pka is not None and pka > PH:
            v = 'ASH' if res.name == 'ASP' else 'GLH'
            why = f'PROPKA pKa {pka:.2f} > {PH}: protonated'
        elif res.name == 'LYS' and pka is not None and pka < PH:
            v, why = 'LYN', f'PROPKA pKa {pka:.2f} < {PH}: neutral'
        elif res.name == 'CYS':
            v = 'CYX' if res.index in ss else 'CYS'
            why = 'disulfide' if res.index in ss else 'free thiol'
        else:
            why = None
        if res.name == 'TYR' and pka is not None and pka < PH:
            flagged.append(dict(residue=f'{ch}:{res.name}{num}', pka=round(pka, 2),
                                note='PROPKA puts this Tyr below pH 7.4; left neutral '
                                     '(no OpenMM variant). Not in any docking box.'))
        variants.append(v)
        if why and v not in (None, 'CYS'):
            decisions.append(dict(residue=f'{ch}:{res.name}{num}', variant=v, reason=why))
    report['protonation_variants'] = decisions
    report['protonation_flagged'] = flagged
    report['protonation_assumed'] = flagged_assumed
    return variants


def add_hydrogens(fixed_pdb, out_pdb, his_override, pkas, report, his_assumed=None):
    from openmm.app import PDBFile, ForceField, Modeller
    pdb = PDBFile(str(fixed_pdb))
    ff = ForceField('amber14-all.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    variants = variant_list(modeller.topology, modeller.positions, pkas, his_override, report,
                            his_assumed)
    applied = modeller.addHydrogens(ff, pH=PH, variants=variants)
    report['n_hydrogens_added'] = sum(1 for v in applied if v is not None) if applied else None
    with open(out_pdb, 'w') as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    return modeller, ff


def minimise(pdb_in, pdb_out, deposited_atoms, report, k=1.0e5, max_iter=2000):
    """Restrained minimisation: deposited heavy atoms held in place, everything else free.

    Two sets of atoms in the file are guesses rather than measurements -- the atoms PDBFixer built
    (the 66-67 loop of 3TZF chain B and any incomplete side chain) and every hydrogen, whose torsion
    about an O-H or N-H bond Modeller has to choose.  Left as built they clash: the raw potential
    energy of the 3TZF receptor was 9.85e6 kJ/mol (T1_3TZF/prep_report.json), and Meeko's per-residue
    RDKit perception then saw a five-bonded carbon in the rebuilt Ala66 and refused the whole file.
    Minimising with the deposited atoms restrained at 1e5 kJ mol-1 nm-2 (about 240 kcal mol-1 A-2)
    relaxes the guessed atoms without moving the crystal structure: the report records the largest
    displacement of a restrained atom, which must stay well below the coordinate error of the
    diffraction data.

    Measured over the seven prepared receptors (their prep_report.json files): restrained atoms move
    at most 0.069-0.197 A (rms 0.015-0.023 A), while the *free* atoms -- hydrogens and rebuilt atoms,
    which is what this step exists to relax -- move up to 1.5-8.0 A.  Both are reported, because
    quoting only the restrained numbers would describe the receptor as barely changed when parts of
    it were rebuilt.
    """
    import numpy as _np
    from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
    from openmm.app import ForceField, NoCutoff, PDBFile, Simulation

    pdb = PDBFile(str(pdb_in))
    ff = ForceField('amber14-all.xml')
    system = ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None,
                             rigidWater=False)
    restraint = CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
    restraint.addGlobalParameter('k', k * unit.kilojoule_per_mole / unit.nanometer ** 2)
    for nm in ('x0', 'y0', 'z0'):
        restraint.addPerParticleParameter(nm)
    pos = pdb.positions
    restrained = []
    for atom in pdb.topology.atoms():
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        if key in deposited_atoms:
            restraint.addParticle(atom.index, pos[atom.index].value_in_unit(unit.nanometer))
            restrained.append(atom.index)
    system.addForce(restraint)
    sim = Simulation(pdb.topology, system,
                     LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                              0.001 * unit.picoseconds),
                     Platform.getPlatformByName('CPU'))
    sim.context.setPositions(pos)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=max_iter)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    new = _np.array(state.getPositions().value_in_unit(unit.angstrom))
    old = _np.array(pos.value_in_unit(unit.angstrom))
    disp = _np.linalg.norm(new - old, axis=1)
    free_idx = _np.setdiff1d(_np.arange(len(old)), _np.array(restrained, dtype=int))
    report['minimisation'] = dict(
        energy_before_kJ_mol=round(float(e0.value_in_unit(unit.kilojoule_per_mole)), 1),
        energy_after_kJ_mol=round(float(state.getPotentialEnergy()
                                        .value_in_unit(unit.kilojoule_per_mole)), 1),
        n_restrained=len(restrained), n_free=int(len(free_idx)),
        max_displacement_restrained_A=round(float(disp[restrained].max()), 3),
        rms_displacement_restrained_A=round(float(_np.sqrt((disp[restrained] ** 2).mean())), 3),
        max_displacement_free_A=round(float(disp[free_idx].max()), 3) if len(free_idx) else None,
        restraint_k_kJ_mol_nm2=k)
    with open(pdb_out, 'w') as fh:
        PDBFile.writeFile(pdb.topology, state.getPositions(), fh, keepIds=True)


# --------------------------------------------------------------------------------------------------
# HET groups -> RDKit -> PDBQT
# --------------------------------------------------------------------------------------------------
def het_to_rdkit(pdb_block, resname, add_h=True, deprotonate_phosphate=True):
    """Crystal coordinates + CCD bond orders.  MolFromPDBBlock gets connectivity from distances but
    no bond orders; AssignBondOrdersFromTemplate transfers them from the RCSB ideal-geometry
    component, which is the reference definition of the component."""
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
    if mol is None:
        raise SystemExit(f'{resname}: RDKit could not read the PDB block')
    tpl_path = CCD / f'{resname}_ideal.sdf'
    if not tpl_path.exists():
        raise SystemExit(f'{resname}: no CCD template at {tpl_path}')
    tpl = Chem.MolFromMolFile(str(tpl_path), removeHs=True)
    if tpl is None:
        raise SystemExit(f'{resname}: could not read {tpl_path}')
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES, catchErrors=True)
    mol = AssignBondOrdersFromTemplate(tpl, mol)
    Chem.SanitizeMol(mol)
    if add_h:
        mol = Chem.AddHs(mol, addCoords=True)
    if deprotonate_phosphate:
        # At pH 7.4 a phosphate/pyrophosphate ester is ionised; the CCD ideal component is drawn
        # neutral.  Removing the P-O-H protons turns each -OH into an acceptor, which is what the
        # cofactor actually presents to the pocket.
        rw = Chem.RWMol(mol)
        drop = []
        for atom in rw.GetAtoms():
            if atom.GetSymbol() != 'O':
                continue
            nbrs = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'P' in nbrs:
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'H':
                        drop.append(n.GetIdx())
                        atom.SetFormalCharge(-1)
                        atom.SetNoImplicit(True)
        for idx in sorted(set(drop), reverse=True):
            rw.RemoveAtom(idx)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
    return mol


def res_pdb_block(cn, res, atoms):
    lines, serial = [], 1
    for a in atoms:
        name = a.name if len(a.name) >= 4 or len(a.element.name) == 2 else f' {a.name}'
        lines.append(f'HETATM{serial:5d} {name:<4s} {res.name:>3s} {cn:1s}{res.seqid.num:4d}    '
                     f'{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}  1.00  0.00'
                     f'{"":10s}{a.element.name.upper():>2s}')
        serial += 1
    return '\n'.join(lines) + '\nEND\n'


def rdkit_to_rigid_pdbqt_atoms(mol, resname, chain, seqid):
    """Type a cofactor with Meeko's ligand machinery and return plain ATOM records for the rigid
    receptor.  All torsions are irrelevant here because we only take the coordinates, charges and
    AutoDock types; the torsion tree is discarded."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    prep = MoleculePreparation(rigid_macrocycles=True)
    setups = prep(mol)
    txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise SystemExit(f'{resname}: Meeko refused to write PDBQT: {err}')
    out = []
    for line in txt.splitlines():
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        # rewrite the residue name / chain / number so PLIP and ProLIF see the real cofactor
        out.append('ATOM  ' + line[6:17] + f'{resname:>3s} {chain:1s}{seqid:4d}' + line[26:])
    return out


def metal_pdbqt_atom(res, chain, atoms):
    t, q = METAL_AD_TYPE[res.name]
    a = atoms[0]
    return (f'ATOM  {1:5d} {res.name:<4s} {res.name:>3s} {chain:1s}{res.seqid.num:4d}    '
            f'{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}  1.00  0.00    {q:+6.3f} {t:<2s}')


def water_pdbqt_atoms(res, chain, atoms):
    """Placing two hydrogens on a metal-bound water would fix a donor direction the crystal does not
    determine, so the oxygen alone is written (type OA, i.e. acceptor).  The report says so."""
    a = atoms[0]
    return [f'ATOM  {1:5d}  O   HOH {chain:1s}{res.seqid.num:4d}    '
            f'{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}  1.00  0.00    -0.411 OA']


def write_vina_config(path, receptor, centre, size=BOX_SIZE, exhaustiveness=32, num_modes=20,
                      energy_range=3.0, cpu=3):
    """A complete, standalone Vina config for the target.

    run_vina.py passes everything on the command line and records it in results.json, so this file is
    not what the reported runs use -- it is here so that a reader (or a referee) can reproduce one
    docking by hand with `vina --config vina_config.txt --ligand <lig>.pdbqt`, without reading any of
    our Python.  Meeko's own --write_vina_box output has the box only."""
    path.write_text(
        f'receptor = {receptor}\n'
        f'center_x = {centre[0]:.3f}\ncenter_y = {centre[1]:.3f}\ncenter_z = {centre[2]:.3f}\n'
        f'size_x = {size:.1f}\nsize_y = {size:.1f}\nsize_z = {size:.1f}\n'
        f'exhaustiveness = {exhaustiveness}\nnum_modes = {num_modes}\n'
        f'energy_range = {energy_range}\ncpu = {cpu}\n')


def renumber(lines):
    out = []
    for i, ln in enumerate(lines, start=1):
        out.append(ln[:6] + f'{i:5d}' + ln[11:])
    return out


def recheck(target, out, spec, drop=()):
    """Re-report an existing receptor without rebuilding it (--recheck).

    Why this exists: prep_report.json is written once, by whichever version of this script prepared
    that receptor, and PDBFixer's rebuild is not bit-reproducible, so re-running the preparation to
    refresh the report would produce a *different* receptor and invalidate every docking run made with
    the old one.  Everything recomputed here is a pure function of files already on disk, and the
    receptor hashes are compared with the ones recorded at preparation time, so the output states
    which receptor the numbers belong to.
    """
    import hashlib
    from scipy.spatial import cKDTree
    old = json.loads((out / 'prep_report.json').read_text())
    modelled = sorted({int(ln[22:26]) for ln in (out / 'raw_selected.pdb').read_text().splitlines()
                       if ln.startswith('ATOM')})
    seq, unobs = deposited_sequence(spec['pdb'], spec['chain'])
    gaps = classify_gaps(unobs, modelled)
    rebuilt_auth = sorted(n for g in old.get('missing_residues_rebuilt_auth_numbering', [])
                          for n in g['auth_numbers'])
    lines = (out / 'receptor.pdbqt').read_text().splitlines()
    heavy = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in lines
                      if ln.startswith('ATOM') and ln[77:79].strip() != 'HD'])
    tree = cKDTree(heavy)
    nn, _ = tree.query(heavy, k=2)
    pairs = tree.query_pairs(1.1, output_type='ndarray')
    rep = dict(
        target=target, directory=out.name, rechecked_at_utc=None,
        note='recomputed from the files on disk; no receptor was rebuilt',
        gap_report=dict(internal_gaps_auth=gaps['internal'],
                        terminal_unobserved_auth=gaps['terminal'],
                        modelled_range=gaps['modelled_range'],
                        rebuilt_auth_from_prep_report=rebuilt_auth,
                        rebuilt_keys_from_prep_report=sorted(
                            old.get('missing_residues_rebuilt', {})),
                        internal_gaps_in_prep_report=old.get('internal_gaps', '<absent>'),
                        # three outcomes, not two: the field can agree, disagree, or predate the code
                        # that writes it -- six of the seven receptors on disk are the third case, and
                        # calling that "inconsistent" would be wrong.
                        prep_report_states_gaps=(
                            'absent (receptor prepared before the field existed)'
                            if 'internal_gaps' not in old
                            else 'agrees' if old['internal_gaps'] == gaps['internal']
                            else 'DISAGREES')),
        geometry=dict(n_heavy_atoms_checked=int(len(heavy)),
                      min_heavy_heavy_distance_A=round(float(nn[:, 1].min()), 3),
                      heavy_pairs_below_1p1A=int(len(pairs)), geometry_ok=bool(len(pairs) == 0)),
        hashes=dict(
            receptor_pdbqt_sha256=hashlib.sha256((out / 'receptor.pdbqt').read_bytes()).hexdigest(),
            receptor_pdb_sha256=hashlib.sha256((out / 'receptor.pdb').read_bytes()).hexdigest()),
        protonation_assumed=old.get('protonation_assumed', '<absent: receptor predates the field>'),
        dropped_het=list(drop) or old.get('dropped_het', []))
    rep['hashes']['matches_prep_report'] = bool(
        rep['hashes']['receptor_pdbqt_sha256'] == old.get('checks', {}).get('receptor_pdbqt_sha256'))
    (out / 'prep_report_recheck.json').write_text(json.dumps(rep, indent=2) + '\n')
    g, ge = rep['gap_report'], rep['geometry']
    print(f"{out.name}: internal gaps {g['internal_gaps_auth'] or 'none'} "
          f"(prep_report said {g['internal_gaps_in_prep_report']}, "
          f"rebuilt {g['rebuilt_auth_from_prep_report'] or g['rebuilt_keys_from_prep_report']}), "
          f"min heavy-heavy {ge['min_heavy_heavy_distance_A']} A over {ge['n_heavy_atoms_checked']} "
          f"atoms, pairs<1.1 A {ge['heavy_pairs_below_1p1A']}, "
          f"hash matches prep_report: {rep['hashes']['matches_prep_report']}")
    return 0


# --------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('target', choices=sorted(TARGETS))
    ap.add_argument('--no-rebuild', action='store_true', help='do not rebuild internal loop gaps')
    ap.add_argument('--drop', nargs='*', default=[], help='HET names to drop from "keep"')
    ap.add_argument('--suffix', default='', help='write to dock/<TARGET><suffix>/')
    ap.add_argument('--seed', type=int, default=20260928,
                    help='PDBFixer rebuild seed; vary it to test how much a rebuilt loop matters')
    ap.add_argument('--recheck', action='store_true',
                    help='do not build anything: recompute the report fields that are pure functions '
                         'of the files already in dock/<TARGET><suffix>/ (gaps, geometry, hashes) and '
                         'write prep_report_recheck.json. Use this to re-report an existing receptor '
                         'without rebuilding it, because a rebuild would change the receptor that the '
                         'docking runs were made with.')
    a = ap.parse_args()

    spec = TARGETS[a.target]
    out = DOCK / (a.target + a.suffix)
    if a.recheck:
        return recheck(a.target, out, spec, a.drop)
    out.mkdir(parents=True, exist_ok=True)
    report = dict(target=a.target, pdb=spec['pdb'], chain=spec['chain'], ph=PH,
                  rebuild_loops=spec['rebuild_loops'] and not a.no_rebuild,
                  dropped_het=a.drop, keep_requested=spec['keep'])

    st, protein, refres, cofactors, metals, waters, altlocs, seq, unobs = select(spec, drop=a.drop)
    report['n_protein_residues'] = len(protein)
    report['cofactors_kept'] = [f'{c}/{r.name}{r.seqid.num}' for c, r, _ in cofactors]
    report['metals_kept'] = [f'{c}/{r.name}{r.seqid.num}' for c, r, _ in metals]
    report['waters_kept'] = [dict(id=f'{c}/HOH{r.seqid.num}', metal_dist=d)
                             for c, r, _, d in waters]
    report['altloc_chosen'] = altlocs
    report['n_residues_with_altlocs'] = len(altlocs)
    report['water_rule'] = (f'kept only waters within {METAL_WATER_CUTOFF} A of a kept metal; '
                            'written as a bare OA oxygen (no hydrogens placed)')

    report['deposited_sequence_length'] = len(seq)
    report['unobserved_residues_auth_numbering'] = unobs
    write_pdb(out / 'raw_selected.pdb', protein,
              seqres=seqres_records(spec['chain'], seq))
    run_pdbfixer(out / 'raw_selected.pdb', out / 'fixed_heavy.pdb',
                 rebuild_loops=report['rebuild_loops'], report=report, seed=a.seed)
    pkas = run_propka(out / 'fixed_heavy.pdb', out, report)
    add_hydrogens(out / 'fixed_heavy.pdb', out / '_with_h.pdb',
                  spec.get('his_templates', {}), pkas, report, spec.get('his_assumed'))
    # Atoms that really were measured: everything in the selection written before PDBFixer ran.
    deposited = {(ln[21], ln[22:26].strip(), ln[12:16].strip())
                 for ln in (out / 'raw_selected.pdb').read_text().splitlines()
                 if ln.startswith('ATOM')}
    minimise(out / '_with_h.pdb', out / 'fixed_protein.pdb', deposited, report)
    (out / '_with_h.pdb').unlink()

    # --- reference ligand, crystal coordinates, CCD bond orders -----------------------------------
    cn, res, refatoms = refres
    ref = het_to_rdkit(res_pdb_block(cn, res, refatoms), res.name, add_h=False,
                       deprotonate_phosphate=False)
    ref.SetProp('_Name', f'{spec["pdb"]}_{res.name}{res.seqid.num}_crystal')
    refpath = out / f'ref_{res.name}.sdf'
    Chem.SDWriter(str(refpath)).write(ref)
    conf = ref.GetConformer()
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(ref.GetNumAtoms())])
    centre = pos.mean(axis=0)
    report['reference_ligand'] = dict(name=res.name, chain=cn, seqid=res.seqid.num,
                                      heavy_atoms=ref.GetNumAtoms(),
                                      smiles=Chem.MolToSmiles(ref),
                                      centre=[round(float(v), 3) for v in centre])

    # --- box ---------------------------------------------------------------------------------------
    box = dict(target=a.target, pdb=spec['pdb'],
               center=[round(float(v), 3) for v in centre],
               size=[BOX_SIZE, BOX_SIZE, BOX_SIZE],
               source=f'heavy-atom centroid of {res.name}{res.seqid.num} chain {cn} in {spec["pdb"]}',
               ligand_extent=[round(float(v), 2) for v in pos.max(0) - pos.min(0)])
    (out / 'box.json').write_text(json.dumps(box, indent=2) + '\n')
    write_vina_config(out / 'vina_config.txt', 'receptor.pdbqt', centre)

    # --- receptor PDB (protein with H + HET) -------------------------------------------------------
    prot_lines = [ln for ln in (out / 'fixed_protein.pdb').read_text().splitlines()
                  if ln.startswith(('ATOM', 'HETATM'))]
    het_groups = list(cofactors) + list(metals) + [(c, r, at) for c, r, at, _ in waters]
    write_pdb(out / '_het.pdb', het_groups, hetatm=True)
    het_lines = [ln for ln in (out / '_het.pdb').read_text().splitlines() if ln.startswith('HETATM')]
    (out / 'receptor.pdb').write_text('\n'.join(renumber(prot_lines + het_lines)) + '\nEND\n')
    (out / '_het.pdb').unlink()

    # --- receptor PDBQT ---------------------------------------------------------------------------
    cmd = ['mk_prepare_receptor.py', '--read_pdb', str(out / 'fixed_protein.pdb'),
           # -p writes the PDBQT we use; Meeko's -j parameterised-receptor JSON is 6.4 MB per target
           # and nothing here reads it, so it is not requested.
           '-o', str(out / 'protein'), '-p',
           '--box_size', *[str(BOX_SIZE)] * 3,
           '--box_center', *[f'{v:.3f}' for v in centre],
           '--default_altloc', spec.get('default_altloc', 'A'), '-v', str(out / 'protein_box.txt')]
    if spec.get('his_templates'):
        cmd += ['--set_template', ','.join(f'{spec["chain"]}:{n}={t}'
                                           for n, t in spec['his_templates'].items())]
    # When a gap is left open (--no-rebuild), the residues flanking it are chain ends with a dangling
    # bond.  Meeko will not match a template for them unless they are declared blunt: index 2 is the
    # backbone C of the residue before the gap, index 0 the backbone N of the residue after it.
    modelled = sorted(int(r.seqid.num) for _, r, _ in protein)
    gaps = classify_gaps(unobs, modelled)
    internal_gaps = gaps['internal']
    report['internal_gaps'] = internal_gaps
    # The two views of the same fact, side by side: what the author numbering says is missing inside
    # the chain, and what PDBFixer actually rebuilt.  They must name the same residues; if they ever
    # disagree, the report says so instead of one of them quietly being empty.
    rebuilt_auth = sorted(n for g in report.get('missing_residues_rebuilt_auth_numbering', [])
                          for n in g['auth_numbers'])
    report['gap_report'] = dict(
        internal_gaps_auth=internal_gaps, terminal_unobserved_auth=gaps['terminal'],
        modelled_range=gaps['modelled_range'], rebuilt_auth=rebuilt_auth,
        rebuilt_matches_internal_gaps=bool(rebuilt_auth == internal_gaps) if report['rebuild_loops']
        else None,
        note=('internal_gaps_auth comes from _pdbx_poly_seq_scheme (classify_gaps); rebuilt_auth is '
              "PDBFixer's own list translated to author numbering. With --no-rebuild the gap is left "
              'open and the flanking residues are declared blunt to Meeko.'))
    if internal_gaps and not report['rebuild_loops']:
        blunt, runs = [], []
        for n in internal_gaps:
            if runs and n == runs[-1][-1] + 1:
                runs[-1].append(n)
            else:
                runs.append([n])
        for run in runs:
            blunt += [f'{spec["chain"]}:{run[0] - 1}=2', f'{spec["chain"]}:{run[-1] + 1}=0']
        cmd += ['--blunt_ends', ','.join(blunt)]
        report['blunt_ends'] = blunt
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=out)
    report['mk_prepare_receptor'] = dict(returncode=r.returncode,
                                         stdout=r.stdout[-4000:], stderr=r.stderr[-4000:],
                                         cmd=' '.join(cmd))
    prot_pdbqt = out / 'protein.pdbqt'
    if r.returncode != 0 or not prot_pdbqt.exists():
        (out / 'prep_report.json').write_text(json.dumps(report, indent=2) + '\n')
        raise SystemExit(f'mk_prepare_receptor failed; see {out}/prep_report.json')

    cof_lines = []
    for c, res_, at in cofactors:
        mol = het_to_rdkit(res_pdb_block(c, res_, at), res_.name, add_h=True)
        cof_lines += rdkit_to_rigid_pdbqt_atoms(mol, res_.name, c, res_.seqid.num)
    for c, res_, at in metals:
        cof_lines.append(metal_pdbqt_atom(res_, c, at))
    for c, res_, at, _ in waters:
        cof_lines += water_pdbqt_atoms(res_, c, at)
    (out / 'cofactors.pdbqt').write_text('\n'.join(cof_lines) + '\n' if cof_lines else '')

    rec_lines = [ln for ln in prot_pdbqt.read_text().splitlines()
                 if ln.startswith(('ATOM', 'HETATM'))]
    (out / 'receptor.pdbqt').write_text('\n'.join(renumber(rec_lines + cof_lines)) + '\n')
    report['n_receptor_atoms'] = dict(protein=len(rec_lines), cofactor_block=len(cof_lines),
                                      total=len(rec_lines) + len(cof_lines))

    # --- checks -----------------------------------------------------------------------------------
    checks = {}
    types = sorted({ln[77:79].strip() for ln in (out / 'receptor.pdbqt').read_text().splitlines()
                    if ln.startswith('ATOM')})
    checks['atom_types'] = types
    checks['all_types_known'] = all(t in {'C', 'A', 'N', 'NA', 'OA', 'SA', 'S', 'P', 'HD', 'F', 'Cl',
                                          'Br', 'I', 'Zn', 'Mg', 'Mn', 'Ca', 'Fe', 'Na', 'K'}
                                    for t in types)
    rpos = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                     for ln in (out / 'receptor.pdbqt').read_text().splitlines()
                     if ln.startswith('ATOM')])
    inbox = np.abs(rpos - centre) <= BOX_SIZE / 2
    checks['receptor_atoms_in_box'] = int(inbox.all(axis=1).sum())
    # geometry sanity on the final receptor: no two heavy atoms closer than a bond, and nothing
    # overlapping.  A rebuilt loop that came out clashing shows up here, not in the docking scores.
    # The minimum heavy-heavy distance is measured over the *whole* receptor, not only over the pairs
    # below the 1.1 A threshold: reporting it only when the check fails (as this did before the
    # 2026-09-28 verification) means a passing receptor states no number at all, and "no clash found"
    # cannot be told apart from "the probe found nothing".
    from scipy.spatial import cKDTree
    heavy = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                      for ln in (out / 'receptor.pdbqt').read_text().splitlines()
                      if ln.startswith('ATOM') and ln[77:79].strip() != 'HD'])
    tree = cKDTree(heavy)
    pairs = tree.query_pairs(1.1, output_type='ndarray')
    nn, _ = tree.query(heavy, k=2)          # k=2: the atom itself, then its nearest neighbour
    checks['n_heavy_atoms_checked'] = int(len(heavy))
    checks['heavy_pairs_below_1.1A'] = int(len(pairs))
    checks['min_heavy_heavy_distance_A'] = round(float(nn[:, 1].min()), 3)
    checks['min_heavy_heavy_distance_threshold_A'] = 1.1
    checks['geometry_ok'] = bool(len(pairs) == 0)
    # Provenance.  Two runs of this script with the same seed give deposited-atom positions that agree
    # to 0.007 A rms (max 0.18 A -- builder's measurement, not re-measured in the 2026-09-28
    # verification), but the atoms PDBFixer *builds* can land in different rotamers --
    # measured up to 7.1 A apart for a surface Arg, and up to 3.4 A for the rebuilt 3TZF Ala66-Ala67
    # loop, 11 of whose atoms are inside the docking box.  PDBFixer's rebuild runs its own OpenMM
    # minimisation, which is not bit-reproducible on the multithreaded CPU platform, so seeding the
    # RNGs is not enough.  The docking is insensitive to this at the level we report (an independent
    # rebuild from seed 12345 changed the 3TZF redocking RMSD from 1.921 to 1.886 A and the score from
    # -7.64 to -7.62 kcal/mol), but the hash pins exactly which receptor produced which numbers.
    import hashlib
    checks['receptor_pdbqt_sha256'] = hashlib.sha256(
        (out / 'receptor.pdbqt').read_bytes()).hexdigest()
    checks['receptor_pdb_sha256'] = hashlib.sha256((out / 'receptor.pdb').read_bytes()).hexdigest()
    # the co-crystal ligand must sit wholly inside the box, with room to spare
    checks['ref_ligand_margin_A'] = round(float((BOX_SIZE / 2 - np.abs(pos - centre).max())), 2)
    report['checks'] = checks
    (out / 'prep_report.json').write_text(json.dumps(report, indent=2) + '\n')

    print(f'{a.target}: receptor {report["n_receptor_atoms"]}, box centre '
          f'{[round(float(v), 2) for v in centre]}, types {types}')
    print(f'  protonation overrides: {len(report["protonation_variants"])}, '
          f'waters kept: {len(report["waters_kept"])}, '
          f'ref margin {checks["ref_ligand_margin_A"]} A, '
          f'atoms in box {checks["receptor_atoms_in_box"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
