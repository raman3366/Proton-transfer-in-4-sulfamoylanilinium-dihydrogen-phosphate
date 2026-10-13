"""
Control receptor: the ATP site of E. coli DNA gyrase B (4DUH), an enzyme that is not a target of
sulfanilamide.

    python dock/control_gyrb.py run        # every control job, one after another (resumes)
    python dock/control_gyrb.py analyse    # dock/control/control_summary.json and the pose files
    python dock/control_gyrb.py report     # dock/CONTROL.md

Why this exists.  A docking score says something about a ligand only if it changes when the receptor
is not a target of the ligand.  The five receptors of the study (dock/targets.py, T1-T5) all belong to
sulfonamide pharmacology (DHPS, hCA II and Sul2, the DHPS variant that rejects sulfonamides), so nothing
in gate G2 or in the production run measures that.  The
gyrase B ATP site is a pocket with its own ligands (the bithiazole of 4DUH, the aminocoumarins) and no
known role for sulfanilamide, whose antibacterial target is DHPS; the fluoroquinolones act on gyrase
too, but at the GyrA-DNA cleavage complex, not here.  Docking the same species into it with the same
protocol turns "the scores cannot argue potency" from a sentence into a measurement.

What is run (JOBS below), all with the protocol of the five receptors (T1-T5) -- AutoDock Vina 1.2.7, Vina and
Vinardo scoring, 20 A cube on the centroid of the crystal ligand, exhaustiveness 32, ten seeds:

*   validation: the 4DUH inhibitor redocked from a fresh ETKDG conformer per seed.  Clorobiocin, the
    aminocoumarin of 1KZN, would be the cross-docking ligand, but its crystal pose, carried into the
    4DUH frame by a sequence-aligned Ca superposition, overlaps residues that are ordered in 4DUH and
    not modelled in 1KZN (cross_reference() measures the overlap).  The receptor cannot hold that
    pose, so clorobiocin is not docked and the control has no cross-docking test;
*   the eight species of the production matrix (sulfanilamide, its anion and cation from their DFT
    geometries; sulfamethoxazole and acetazolamide, neutral and anionic, and p-aminobenzoate from
    ETKDG conformers), so that every one of them has a score in a pocket that is the target of none;
*   ciprofloxacin, which acts on gyrase but not at this site;
    (Vinardo is run for the inhibitor, sulfanilamide and ciprofloxacin only: it refits Vina's terms,
    so a second scoring function is worth its time for the three ligands the argument rests on;)
*   the same key jobs in a receptor whose unobserved lid loop (109-116) was rebuilt from another
    PDBFixer seed, and in the primary receptor with those eight built residues deleted
    (loop_free_variant), because the top pose of sulfanilamide turned out to touch them;
*   the same key jobs in the primary receptor with the ordered water of the crystal put back
    (water_variant): a diagnostic, added when the redocking showed that the top-scored pose of the
    inhibitor lies in a channel that holds water in the crystal.  It is not the protocol of the
    five receptors, which were all docked without such water: its scores are given beside those of
    the receptor as prepared and are always marked as coming from the diagnostic;
*   a search over the whole fragment for sulfanilamide, with the exhaustiveness scaled by the box
    volume, to see whether the ATP pocket is where an unrestricted search ends.

Nothing here touches the five receptors T1-T5: the directories are named C1_*, which the 'T*' globs of
analyse.py, production_analysis.py and contacts.py never match, and the output goes to dock/control/.
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.warning')
DOCK = pathlib.Path(__file__).resolve().parent
ROOT = DOCK.parent
OUT = DOCK / 'control'
PY = sys.executable
C1 = 'C1_4DUH'
C1_ALT = 'C1_4DUH_seed12345'
C1_WAT = 'C1_4DUH_wat'
C1_NOLOOP = 'C1_4DUH_noloop'
WATER_CLASH = 2.4        # A: a crystal water closer than this to a heavy atom overlaps it
BOX = 20.0
EXH = 32
SEEDS = 10
LONG_SEEDS = 5           # the whole-fragment search
PAIR_CUT = 2.5           # A: a heavy-atom pair of the carried-over clorobiocin pose and the receptor closer than this overlaps
WHOLE_MARGIN = 5.0       # A beyond the receptor's extent on every side, whole-fragment search

SFA = ('sfa_neutral', 'sfa_cation', 'sfa_anion')          # docked from their DFT geometries
CONTROLS = ('smx_anion', 'smx_neutral', 'paba_anion', 'azm_neutral', 'azm_anion')


def geom_of(lig):
    return 'dft' if lig in SFA else 'etkdg'


def job(receptor, lig, scoring, kind, seeds=SEEDS, box=None, whole=False):
    """One docking job and the run directory run_vina.py gives it."""
    suffix = {'validation': '', 'production': '__prod', 'whole': '__whole'}[kind]
    tag = f'{lig}__{scoring}' + (f'__box{box:g}' if box else '') + suffix
    return dict(receptor=receptor, ligand=lig, scoring=scoring, kind=kind, seeds=seeds, box=box,
                whole=whole, geometry=geom_of(lig), suffix=suffix, tag=tag,
                rundir=f'{receptor}/runs/{tag}')


def jobs():
    """The control matrix, in the order it is run: what decides the outcome first."""
    J = [job(C1, 'rli_anion', 'vina', 'validation'),
         job(C1, 'sfa_neutral', 'vina', 'production'),
         job(C1, 'cfx_zwit', 'vina', 'production'),
         job(C1, 'rli_anion', 'vinardo', 'validation'),
         job(C1, 'sfa_neutral', 'vinardo', 'production'),
         job(C1, 'cfx_zwit', 'vinardo', 'production')]
    # the ordered water of the crystal put back (diagnostic): the redock, the primary ligand, the comparator
    J += [job(C1_WAT, 'rli_anion', 'vina', 'validation'),
          job(C1_WAT, 'sfa_neutral', 'vina', 'production'),
          job(C1_WAT, 'cfx_zwit', 'vina', 'production')]
    # the rest of the production matrix: Vina only (Vinardo is run for the three ligands above)
    for lig in SFA[1:] + CONTROLS:
        J.append(job(C1, lig, 'vina', 'production'))
    # the lid loop rebuilt from another seed: the redock, the primary ligand, the comparator
    J += [job(C1_ALT, 'rli_anion', 'vina', 'validation'),
          job(C1_ALT, 'sfa_neutral', 'vina', 'production'),
          job(C1_ALT, 'cfx_zwit', 'vina', 'production')]
    # Clorobiocin is not docked: its crystal pose (1KZN) cannot be held by this receptor
    # (cross_reference), so no cross-docking test exists, and a score of a pose that the receptor
    # cannot hold would say nothing.
    # whole-fragment search: five seeds, each of them 21 times the Monte Carlo runs of a focused seed
    J += [job(C1, 'sfa_neutral', 'vina', 'whole', seeds=LONG_SEEDS, whole=True)]
    # ... and the primary receptor with the eight built residues deleted: observed atoms only.  Listed
    # last because these three were run on the spare core beside the queue (`run --only C1_4DUH_noloop
    # --cpu 1`; the result does not depend on --cpu, dock/VALIDATION.md section 7).
    J += [job(C1_NOLOOP, 'rli_anion', 'vina', 'validation'),
          job(C1_NOLOOP, 'sfa_neutral', 'vina', 'production'),
          job(C1_NOLOOP, 'cfx_zwit', 'vina', 'production')]
    # The receptor with the crystal water is the one that redocks the inhibitor (CONTROL.md section 3),
    # so the rest of the production matrix is docked into it as well (spare core, `--only C1_4DUH_wat`).
    for lig in SFA[1:] + CONTROLS:
        J.append(job(C1_WAT, lig, 'vina', 'production'))
    return J


# --------------------------------------------------------------------------------------------------
# boxes
# --------------------------------------------------------------------------------------------------
def receptor_xyz(receptor):
    return np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                     for ln in (DOCK / receptor / 'receptor.pdbqt').read_text().splitlines()
                     if ln.startswith('ATOM')])


def chain_waters():
    """(residue number, B factor, xyz) of every water the deposited 4DUH file assigns to chain A."""
    from targets import read_structure
    st = read_structure('4DUH')
    return [(res.seqid.num, res[0].b_iso, np.array([res[0].pos.x, res[0].pos.y, res[0].pos.z]))
            for ch in st[0] if ch.name == 'A' for res in ch if res.name == 'HOH']


def water_variant():
    """C1_4DUH_wat: the primary receptor, atom for atom, plus the ordered water of the crystal.

    Built from the files of C1_4DUH and not by a second run of prep_receptor.py, because PDBFixer's
    loop rebuild is not bit-reproducible and the comparison should change one thing only.  Every
    water of chain A is added as a bare OA oxygen (no hydrogens, as for the metal-bound waters of the
    targets), except a water closer than WATER_CLASH to a heavy atom of the residues PDBFixer built
    (109-116): it sits where a guessed atom now is.  Short contacts between a water and a deposited
    atom are the crystal's own and are kept.  No water is within 2.8 A of the crystal inhibitor, so
    the crystal pose stays possible.  Diagnostic only (see the module docstring)."""
    import hashlib
    import shutil
    src, dst = DOCK / C1, DOCK / C1_WAT
    dst.mkdir(exist_ok=True)
    built = [r for r in heavy_rows(C1) if r[0] in LOOP]
    P = np.array([r[3] for r in built])
    ref = Chem.MolFromMolFile(str(src / 'ref_RLI.sdf'))
    L = np.array([list(ref.GetConformer().GetAtomPosition(i)) for i in range(ref.GetNumAtoms())])
    kept, dropped = [], []
    for num, b, o in chain_waters():
        d = float(np.linalg.norm(P - o, axis=1).min())
        (kept if d >= WATER_CLASH else dropped).append((num, b, o, d))
    qt = [ln for ln in (src / 'receptor.pdbqt').read_text().splitlines() if ln.startswith('ATOM')]
    pdb = [ln for ln in (src / 'receptor.pdb').read_text().splitlines() if ln.startswith(('ATOM', 'HETATM'))]
    for num, b, o, _ in kept:
        qt.append(f'ATOM  {1:5d}  O   HOH A{num:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00  0.00    -0.411 OA')
        pdb.append(f'HETATM{1:5d}  O   HOH A{num:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00{b:6.2f}           O')
    renum = lambda lines: [ln[:6] + f'{i:5d}' + ln[11:] for i, ln in enumerate(lines, start=1)]  # noqa: E731
    (dst / 'receptor.pdbqt').write_text('\n'.join(renum(qt)) + '\n')
    (dst / 'receptor.pdb').write_text('\n'.join(renum(pdb)) + '\nEND\n')
    shutil.copyfile(src / 'ref_RLI.sdf', dst / 'ref_RLI.sdf')
    box = json.loads((src / 'box.json').read_text())
    box['target'] = C1_WAT
    (dst / 'box.json').write_text(json.dumps(box, indent=2) + '\n')
    W = np.array([k[2] for k in kept])
    dl = np.linalg.norm(W[:, None] - L[None], axis=2).min(1)
    centre = np.array(box['center'])
    rep = dict(target=C1_WAT, derived_from=C1,
               what='the primary control receptor with the ordered water of the crystal added; '
                    'diagnostic only, not prepared by prep_receptor.py',
               protein_identical_to=f'{C1}/receptor.pdbqt',
               protein_receptor_pdbqt_sha256=hashlib.sha256((src / 'receptor.pdbqt').read_bytes()).hexdigest(),
               water_source='4DUH, waters assigned to chain A in the deposited file',
               n_waters_in_file=len(kept) + len(dropped), n_waters_kept=len(kept),
               n_waters_dropped=len(dropped), water_clash_cutoff_A=WATER_CLASH,
               waters_dropped=[dict(id=f'A/HOH{n}', b_iso=round(b, 1), nearest_built_heavy_atom_A=round(d, 2))
                               for n, b, o, d in dropped],
               drop_rule=f'within {WATER_CLASH} A of a heavy atom of the rebuilt residues 109-116',
               water_rule='bare OA oxygen, no hydrogens placed',
               n_waters_in_box=int((np.abs(W - centre) <= BOX / 2).all(axis=1).sum()),
               n_waters_within_4A_of_crystal_inhibitor=int((dl <= 4.0).sum()),
               nearest_water_to_crystal_inhibitor_A=round(float(dl.min()), 2),
               receptor_pdbqt_sha256=hashlib.sha256((dst / 'receptor.pdbqt').read_bytes()).hexdigest(),
               receptor_pdb_sha256=hashlib.sha256((dst / 'receptor.pdb').read_bytes()).hexdigest())
    (dst / 'prep_report.json').write_text(json.dumps(rep, indent=2) + '\n')
    return rep


def loop_free_variant():
    """C1_4DUH_noloop: the primary receptor without the eight residues PDBFixer built (109-116).

    Every other atom is that of C1_4DUH, so the receptor holds observed residues only; Ser108 and
    Gly117 keep the atom types they have in the full chain (Vina reads types and positions, not
    bonds, so the open chain ends need no capping).  Sensitivity case, like the second loop build."""
    import hashlib
    import shutil
    src, dst = DOCK / C1, DOCK / C1_NOLOOP
    dst.mkdir(exist_ok=True)
    keep = lambda ln: not (ln.startswith(('ATOM', 'HETATM')) and int(ln[22:26]) in LOOP)  # noqa: E731
    qt = [ln for ln in (src / 'receptor.pdbqt').read_text().splitlines() if ln.startswith('ATOM') and keep(ln)]
    pdb = [ln for ln in (src / 'receptor.pdb').read_text().splitlines()
           if ln.startswith(('ATOM', 'HETATM')) and keep(ln)]
    n_qt = sum(1 for ln in (src / 'receptor.pdbqt').read_text().splitlines() if ln.startswith('ATOM'))
    renum = lambda lines: [ln[:6] + f'{i:5d}' + ln[11:] for i, ln in enumerate(lines, start=1)]  # noqa: E731
    (dst / 'receptor.pdbqt').write_text('\n'.join(renum(qt)) + '\n')
    (dst / 'receptor.pdb').write_text('\n'.join(renum(pdb)) + '\nEND\n')
    shutil.copyfile(src / 'ref_RLI.sdf', dst / 'ref_RLI.sdf')
    box = json.loads((src / 'box.json').read_text())
    box['target'] = C1_NOLOOP
    (dst / 'box.json').write_text(json.dumps(box, indent=2) + '\n')
    rep = dict(target=C1_NOLOOP, derived_from=C1,
               what='the primary control receptor without the residues PDBFixer built (109-116): observed '
                    'residues only; sensitivity case, not prepared by prep_receptor.py',
               removed_residues=sorted(LOOP), n_atoms_removed=n_qt - len(qt), n_receptor_atoms=len(qt),
               receptor_pdbqt_sha256=hashlib.sha256((dst / 'receptor.pdbqt').read_bytes()).hexdigest(),
               receptor_pdb_sha256=hashlib.sha256((dst / 'receptor.pdb').read_bytes()).hexdigest())
    (dst / 'prep_report.json').write_text(json.dumps(rep, indent=2) + '\n')
    return rep


def whole_box(receptor=C1):
    """A box around the whole fragment, and the exhaustiveness that keeps the search effort per unit
    volume of the 20 A runs.

    Vina searches only inside its box, so a ligand can reach every part of the surface only if the box
    extends beyond the outermost receptor atoms; WHOLE_MARGIN (5 A) leaves room for an 11-heavy-atom
    ligand lying flat on the surface.  The edges are rounded up to 0.5 A.  Exhaustiveness is the
    number of independent Monte Carlo runs, so scaling it with the volume keeps the number of runs per
    cubic angstrom the same as in the focused search."""
    path = DOCK / receptor / 'box_whole.json'
    xyz = receptor_xyz(receptor)
    lo, hi = xyz.min(0), xyz.max(0)
    size = [math.ceil((float(v) + 2 * WHOLE_MARGIN) * 2) / 2 for v in hi - lo]
    vol = size[0] * size[1] * size[2]
    box = dict(target=receptor, center=[round(float(v), 3) for v in (lo + hi) / 2], size=size,
               source=f'bounding box of every atom of {receptor}/receptor.pdbqt plus {WHOLE_MARGIN} A '
                      'on each side, edges rounded up to 0.5 A',
               receptor_extent=[round(float(v), 2) for v in hi - lo],
               volume_A3=round(vol, 1), focused_volume_A3=BOX ** 3,
               volume_ratio=round(vol / BOX ** 3, 2),
               exhaustiveness=int(round(EXH * vol / BOX ** 3)),
               exhaustiveness_rule=f'{EXH} x volume ratio, rounded: the same number of Monte Carlo '
                                   'runs per unit volume as the focused search')
    path.write_text(json.dumps(box, indent=2) + '\n')
    return box


# --------------------------------------------------------------------------------------------------
# run
# --------------------------------------------------------------------------------------------------
def ist():
    return subprocess.run(['date', '+%d %b %I:%M %p IST'], capture_output=True, text=True,
                          env={'TZ': 'Asia/Kolkata', 'PATH': '/usr/bin:/bin'}).stdout.strip()


def say(msg):
    line = f'[{ist()}] {msg}'
    print(line, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / 'control.log', 'a') as fh:
        fh.write(line + '\n')


def complete(j, exh):
    p = DOCK / j['rundir'] / 'results.json'
    if not p.exists():
        return False
    r = json.loads(p.read_text())
    ok = [s for s in r.get('seeds', []) if s.get('ok')]
    return (len(ok) == j['seeds'] and r.get('geometry', 'etkdg') == j['geometry']
            and int(r.get('exhaustiveness', -1)) == exh)


def run(only=None, cpu=3):
    (OUT / 'logs').mkdir(parents=True, exist_ok=True)
    wb = whole_box()
    if not (DOCK / C1_WAT / 'receptor.pdbqt').exists():
        wr = water_variant()
        say(f'=== {C1_WAT}: {wr["n_waters_kept"]} of {wr["n_waters_in_file"]} crystal waters added to the '
            f'primary receptor ({wr["n_waters_dropped"]} overlap built atoms)')
    if not (DOCK / C1_NOLOOP / 'receptor.pdbqt').exists():
        lr = loop_free_variant()
        say(f'=== {C1_NOLOOP}: {lr["n_atoms_removed"]} atoms of the built residues removed from the '
            'primary receptor')
    todo = [j for j in jobs() if not only or any(o in j['rundir'] for o in only)]
    say(f'=== control run: {len(todo)} jobs, cpu {cpu}, whole-fragment box {wb["size"]} A '
        f'(volume ratio {wb["volume_ratio"]}, exhaustiveness {wb["exhaustiveness"]})')
    t_start = time.time()
    for i, j in enumerate(todo, start=1):
        exh = wb['exhaustiveness'] if j['whole'] else EXH
        if complete(j, exh):
            say(f'SKIP  {i}/{len(todo)} {j["rundir"]} (complete)')
            continue
        cmd = ['nice', '-n', '10', PY, str(DOCK / 'run_vina.py'), '--target', j['receptor'],
               '--ligand', j['ligand'], '--scoring', j['scoring'], '--seeds', str(j['seeds']),
               '--exhaustiveness', str(exh), '--cpu', str(cpu), '--geometry', j['geometry']]
        if j['box']:
            cmd += ['--box-size', str(j['box'])]
        if j['whole']:
            cmd += ['--box-json', str(DOCK / C1 / 'box_whole.json')]
        if j['suffix']:
            cmd += ['--tag-suffix', j['suffix']]
        say(f'START {i}/{len(todo)} {j["rundir"]} seeds={j["seeds"]} geometry={j["geometry"]} exh={exh}')
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        log = OUT / 'logs' / (j['rundir'].replace('/runs/', '__') + '.log')
        log.write_text(r.stdout + '\n' + r.stderr)
        tail = [ln for ln in r.stdout.splitlines() if 'kcal/mol (n=' in ln]
        say(f'{"DONE " if r.returncode == 0 and complete(j, exh) else "FAIL "} {i}/{len(todo)} '
            f'{j["rundir"]} {time.time() - t0:.0f} s | {tail[-1].split(": ", 1)[-1] if tail else ""}')
        say(f'      elapsed {(time.time() - t_start) / 60:.0f} min')
    say('CONTROL RUN DONE')
    return 0


# --------------------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------------------
INFO = dict(pdb='4DUH', protein='DNA gyrase B, E. coli, 24 kDa ATPase fragment (chain A)',
            crystal_ligand='RLI', crystal_ligand_name="4,5'-bithiazole inhibitor",
            reference_species=['rli_anion'],
            validation='control receptor; validation in dock/control/control_summary.json')
CROSS_SDF = 'ref_cross_X_1KZN_CBN.sdf'
# The three polar contacts of the inhibitor in the deposited structure (measured on ref_RLI.sdf: 2.82,
# 3.06 and 2.82 A); like every site check of analyse.py they are evaluated on the crystal ligand too.
SITE = [('amide N to Asp73 N/O', '[NX3H1;$(N(C(=O)[CH2])c)]', ('ASP', 73, '*NO'), 3.2),
        ('carboxylate O to Arg136 N/O', '[OX1,OX2H1,OX2-;$(O~[CX3](~[OX1,OX2H1,OX2-])c)]',
         ('ARG', 136, '*NO'), 3.4),
        ('aryl-NH-aryl N to Gly101 O', '[NX3H1;$(N(c)c)]', ('GLY', 101, 'O'), 3.3)]
SITE_CUT = 4.0           # A: a pose is "in the ATP site" when its heavy-atom centroid lies within this
#                          distance of a heavy atom of the crystal inhibitor
LOOP = range(109, 117)   # the residues PDBFixer built
NAMES = dict(sfa_neutral='sulfanilamide', sfa_cation='4-sulfamoylanilinium', sfa_anion='sulfonamidate',
             smx_neutral='sulfamethoxazole', smx_anion='sulfamethoxazole anion',
             paba_anion='p-aminobenzoate', azm_neutral='acetazolamide',
             azm_anion='acetazolamide anion', rli_anion="4DUH inhibitor (crystal ligand, RLI)",
             cfx_zwit='ciprofloxacin')
TARGET_PDB = dict(T1_3TZF='3TZF', T5_7S2K='7S2K', T3_6RL9='6RL9', T4_3HS4='3HS4')


def rnd(x, n=3):
    return None if x is None else round(float(x), n)


def heavy_rows(receptor):
    rows = []
    for ln in (DOCK / receptor / 'receptor.pdb').read_text().splitlines():
        if ln.startswith(('ATOM', 'HETATM')) and (ln[76:78].strip() or ln[12:16].strip()[0]) != 'H':
            rows.append((int(ln[22:26]), ln[17:20].strip(), ln[12:16].strip(),
                         np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return rows


def mol_xyz(mol, cid=-1):
    """Heavy-atom coordinates and symbols of one conformer."""
    c = mol.GetConformer(cid)
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != 'H']
    return np.array([list(c.GetAtomPosition(i)) for i in idx]), [mol.GetAtomWithIdx(i).GetSymbol()
                                                                 for i in idx]


def cross_reference():
    """Clorobiocin of 1KZN in the frame of 4DUH, and whether 4DUH has room for it (it has not, so
    clorobiocin is not docked and the control receptor has no cross-docking test).

    The two entries are the same 24 kDa fragment, but 1KZN is not numbered like E. coli GyrB between
    residues 82 and 96, so residue numbers cannot pair the Ca atoms (analyse.transform_reference would
    give 2.27 A with sixteen outliers).  gemmi aligns the two sequences and superposes the aligned Ca
    atoms; the transformed ligand gets its bond orders from the chemical component, like every other
    reference here.  The result is written beside ref_RLI.sdf with the prefix ref_cross_, which
    analyse.pick_reference never mistakes for the receptor's own ligand."""
    import difflib
    import gemmi
    import prep_receptor as PR
    from targets import read_structure
    src, dst = read_structure('1KZN'), read_structure('4DUH')
    pk, pd_ = src[0]['A'].get_polymer(), dst[0]['A'].get_polymer()
    # Residues of 4DUH that have no counterpart in the deposited 1KZN model: the observed residues of the
    # two chains paired by their common blocks (the pairing the superposition uses: same count, same RMSD).
    obs_d = [r for r in pd_ if r.find_atom('CA', '*')]
    obs_k = [r for r in pk if r.find_atom('CA', '*')]
    blocks = difflib.SequenceMatcher(None, [r.name for r in obs_d], [r.name for r in obs_k],
                                     autojunk=False).get_matching_blocks()
    paired = {obs_d[i + k].seqid.num for i, _, n in blocks for k in range(n)}
    unpaired = sorted(r.seqid.num for r in obs_d if r.seqid.num not in paired)
    sup = gemmi.calculate_superposition(pd_, pk, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
    trim = gemmi.calculate_superposition(pd_, pk, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP,
                                         trim_cycles=3, trim_cutoff=2.0)
    res = next(r for ch in src[0] if ch.name == 'A' for r in ch if r.name == 'CBN').clone()
    for a in res:
        a.pos = gemmi.Position(*sup.transform.apply(a.pos))
    atoms = [a for a in res if a.element != gemmi.Element('H')]
    mol = PR.het_to_rdkit(PR.res_pdb_block('A', res, atoms), 'CBN', add_h=False,
                          deprotonate_phosphate=False)
    mol.SetProp('_Name', 'X_1KZN_CBN1_crystal_in_4DUH_frame')
    for t in (C1, C1_ALT):
        w = Chem.SDWriter(str(DOCK / t / CROSS_SDF))
        w.write(mol)
        w.close()
    X, _ = mol_xyz(mol)
    centre = np.array(json.loads((DOCK / C1 / 'box.json').read_text())['center'])
    reach = float(np.abs(X - centre).max())
    rec = dict(source='X_1KZN (1KZN chain A, CBN 1)', file=f'{C1}/{CROSS_SDF}',
               method='gemmi.calculate_superposition on the sequence-aligned Ca atoms of chain A '
                      '(4DUH fixed, 1KZN moved), no trimming',
               n_ca_atoms=int(sup.count), superposition_ca_rmsd_A=rnd(sup.rmsd),
               trimmed_fit=dict(n_ca_atoms=int(trim.count), ca_rmsd_A=rnd(trim.rmsd),
                                note='3 cycles, pairs beyond 2.0 sigma removed; reported only to show '
                                     'how much of the 0.78 A comes from a few mobile residues'),
               heavy_atoms=int(len(X)),
               max_abs_offset_from_box_centre_A=rnd(reach, 2), margin_in_20A_box_A=rnd(BOX / 2 - reach, 2),
               n_residues_paired_by_sequence=len(paired),
               residues_of_4DUH_without_counterpart_in_1KZN=unpaired,
               pair_cutoff_A=PAIR_CUT, overlap={})
    assert len(paired) == sup.count, (len(paired), sup.count)
    for t in (C1, C1_ALT):
        rows = heavy_rows(t)
        P = np.array([r[3] for r in rows])
        D = np.linalg.norm(P[:, None, :] - X[None, :, :], axis=2)
        i = int(D.min(1).argmin())
        close = sorted({(r[0], r[1]) for r, d in zip(rows, D.min(1)) if d < 3.0})
        hit = [rows[a] for a in np.where(D < PAIR_CUT)[0]]             # one entry per overlapping pair
        by_res = {}
        for r in hit:
            by_res[f'{r[1]}{r[0]}'] = by_res.get(f'{r[1]}{r[0]}', 0) + 1
        rec['overlap'][t] = dict(
            min_heavy_heavy_distance_A=rnd(D.min(), 2),
            closest_receptor_atom=f'{rows[i][1]}{rows[i][0]}:{rows[i][2]}',
            n_pairs_below_2p5A=int((D < PAIR_CUT).sum()), n_pairs_below_3p0A=int((D < 3.0).sum()),
            pairs_below_2p5A_by_residue=dict(sorted(by_res.items(), key=lambda kv: (-kv[1], kv[0]))),
            n_pairs_below_2p5A_without_counterpart_in_1KZN=sum(1 for r in hit if r[0] in unpaired),
            residues_below_2p5A_without_counterpart_in_1KZN=[f'{n}{i_}' for i_, n in sorted(
                {(r[0], r[1]) for r in hit if r[0] in unpaired})],
            residues_within_3A=[f'{n}{i_}' for i_, n in close],
            rebuilt_loop_min_distance_A=rnd(min(d for r, d in zip(rows, D.min(1)) if r[0] in LOOP), 2))
    rec['conclusion'] = ('the transferred pose overlaps residues that are ordered in 4DUH and have no '
                         'counterpart in the 1KZN model, so 4DUH cannot hold clorobiocin in its crystal '
                         'pose; clorobiocin is not docked, and the control receptor has no cross-docking '
                         'test')
    return mol, rec


def patch_pipeline():
    """Teach the two analysis modules about the control receptor, at run time and only in this
    process: nothing is written into their tables, so `analyse.py` and `analyse.py --production`
    behave exactly as before."""
    import analyse as A
    import production_analysis as PA
    A.SITE_CHECKS[C1] = SITE
    PA.DFT_XYZ = {}          # no DFT pose rebuilding or MMFF strain here: the control asks about scores
    PA.RES = OUT
    for t in (C1, C1_ALT, C1_WAT, C1_NOLOOP):
        if not (DOCK / t / 'ref_RLI.sdf').exists():
            continue
        PA.TARGET_INFO[t] = INFO
        PA.OWN_REF[t] = 'ref_RLI.sdf'
        own = Chem.MolFromMolFile(str(DOCK / t / 'ref_RLI.sdf'))
        PA._REFS[t] = [('RLI (4DUH, own)', own, dict(kind='own', file=f'{t}/ref_RLI.sdf'))]
    return A, PA


def gate(st, on_ref):
    """The criteria of gate G2 (analyse.py), applied to a control job; nothing is added to G2."""
    pb = st.get('posebusters_physics_pass_fraction') or 0
    return dict(site_checks_pass_on_crystal_reference=None if not on_ref else all(c['passed'] for c in on_ref),
                rmsd_top_mean_A=st['rmsd_top_mean'], rmsd_top_sd_A=st['rmsd_top_sd'],
                success_fraction_1A=st['success_fraction_1A'],
                success_fraction_2A=st['success_fraction_2A'], posebusters_physics=pb,
                passes_strict_1A=bool(st['rmsd_top_mean'] <= 1.0 and st['success_fraction_1A'] >= 0.5
                                      and pb >= 0.8),
                passes_2A_criterion=bool(st['rmsd_top_mean'] <= 2.0 and st['success_fraction_2A'] >= 0.5
                                         and pb >= 0.8))


def site_geometry(A, receptor, rundir, ref_xyz):
    """Per seed, for the top pose: is it in the ATP site, how close does a ligand N or O come to the
    Asp73 carboxylate (the anchor of adenine and of every ATP-site inhibitor) and to the Arg136
    guanidinium, and how far is the rebuilt loop."""
    res = json.loads((DOCK / rundir / 'results.json').read_text())
    rows = heavy_rows(receptor)
    asp = np.array([r[3] for r in rows if r[0] == 73 and r[2] in ('OD1', 'OD2')])
    arg = np.array([r[3] for r in rows if r[0] == 136 and r[2] in ('NE', 'NH1', 'NH2')])
    loop = np.array([r[3] for r in rows if r[0] in LOOP])
    wat = np.array([w[2] for w in chain_waters()])
    out = []
    for s_ in res['seeds']:
        if not s_.get('ok'):
            continue
        pm = A.load_poses(DOCK / s_['out'])
        X, el = mol_xyz(pm, pm.GetConformers()[0].GetId())
        polar = X[[i for i, e in enumerate(el) if e in 'NO']]
        d_site = float(np.linalg.norm(ref_xyz - X.mean(0), axis=1).min())
        out.append(dict(seed=s_['seed'], centroid_to_inhibitor_A=rnd(d_site, 2),
                        in_site=bool(d_site <= SITE_CUT),
                        crystal_waters_overlapped=int(
                            (np.linalg.norm(wat[:, None] - X[None], axis=2).min(1) < WATER_CLASH).sum()),
                        asp73_O_to_ligand_NO_A=rnd(np.linalg.norm(asp[:, None] - polar[None], axis=2).min(), 2),
                        arg136_N_to_ligand_NO_A=rnd(np.linalg.norm(arg[:, None] - polar[None], axis=2).min(), 2),
                        loop_to_ligand_A=rnd(np.linalg.norm(loop[:, None] - X[None], axis=2).min(), 2)
                        if len(loop) else None))
    a = lambda k: np.array([r[k] for r in out])  # noqa: E731
    return dict(per_seed=out, in_site_fraction=rnd(a('in_site').mean()),
                centroid_to_inhibitor_mean_A=rnd(a('centroid_to_inhibitor_A').mean(), 2),
                centroid_to_inhibitor_range_A=[rnd(a('centroid_to_inhibitor_A').min(), 2),
                                               rnd(a('centroid_to_inhibitor_A').max(), 2)],
                crystal_waters_overlapped_mean=rnd(a('crystal_waters_overlapped').mean(), 1),
                crystal_waters_overlapped_range=[int(a('crystal_waters_overlapped').min()),
                                                 int(a('crystal_waters_overlapped').max())],
                asp73_contact_fraction_3p5A=rnd((a('asp73_O_to_ligand_NO_A') <= 3.5).mean()),
                asp73_distance_mean_A=rnd(a('asp73_O_to_ligand_NO_A').mean(), 2),
                arg136_contact_fraction_3p5A=rnd((a('arg136_N_to_ligand_NO_A') <= 3.5).mean()),
                arg136_distance_mean_A=rnd(a('arg136_N_to_ligand_NO_A').mean(), 2),
                loop_min_distance_A=rnd(a('loop_to_ligand_A').min(), 2) if len(loop) else None,
                loop_distance_range_A=[rnd(a('loop_to_ligand_A').min(), 2), rnd(a('loop_to_ligand_A').max(), 2)]
                if len(loop) else None,
                loop_contact_fraction_4A=rnd((a('loop_to_ligand_A') <= 4.0).mean()) if len(loop) else None)


def whole_fragment(A, j, ref_xyz, same_mol_ref, ref_label):
    """Where an unrestricted search puts the ligand: every mode of every seed is classed as inside or
    outside the ATP site, and the top pose is compared with `same_mol_ref` (the crystal pose for the
    inhibitor, the representative pose of the focused search for sulfanilamide)."""
    import contacts as CT
    res = json.loads((DOCK / j['rundir'] / 'results.json').read_text())
    per, outside = [], {}
    for s_ in res['seeds']:
        if not s_.get('ok'):
            continue
        pm = A.load_poses(DOCK / s_['out'])
        cids = [c.GetId() for c in pm.GetConformers()]
        d = [float(np.linalg.norm(ref_xyz - mol_xyz(pm, c)[0].mean(0), axis=1).min()) for c in cids]
        inside = [x <= SITE_CUT for x in d]
        first = next((i for i, ok in enumerate(inside) if ok), None)
        row = dict(seed=s_['seed'], score_top=s_['best'], n_modes=len(cids), wall_s=s_['wall_s'],
                   top_centroid_to_inhibitor_A=rnd(d[0], 2), top_in_site=bool(inside[0]),
                   rank_of_first_in_site_mode=None if first is None else first + 1,
                   score_of_first_in_site_mode=None if first is None else s_['scores'][first],
                   n_modes_in_site=int(sum(inside)))
        if same_mol_ref is not None:
            row['rmsd_top_to_reference_A'] = rnd(A.symm_rmsd(same_mol_ref, pm, cids[:1])[0])
        if not inside[0]:
            top = Chem.Mol(pm)
            top.RemoveAllConformers()
            top.AddConformer(Chem.Conformer(pm.GetConformer(cids[0])), assignId=True)
            near, _ = CT.geometric_contacts(top, DOCK / j['receptor'] / 'receptor.pdb')
            row['top_pose_nearest_residues'] = list(near)[:6]
            outside[s_['seed']] = row['top_pose_nearest_residues']
        per.append(row)
    sc = np.array([r['score_top'] for r in per])
    ins = np.array([r['top_in_site'] for r in per])
    gap = [r['score_of_first_in_site_mode'] - r['score_top'] for r in per
           if r['score_of_first_in_site_mode'] is not None]
    out = dict(job=j['rundir'], box=res['box'], exhaustiveness=res['exhaustiveness'], n_seeds=len(per),
               site_definition=f'heavy-atom centroid of the pose within {SITE_CUT} A of a heavy atom of '
                               'the crystal inhibitor',
               score_mean=rnd(sc.mean()), score_sd=rnd(sc.std(ddof=1)) if len(sc) > 1 else 0.0,
               top_pose_in_site_fraction=rnd(ins.mean()), n_top_in_site=int(ins.sum()),
               seeds_with_any_in_site_mode=int(sum(r['rank_of_first_in_site_mode'] is not None for r in per)),
               in_site_mode_minus_top_score=dict(mean=rnd(np.mean(gap)), max=rnd(max(gap))) if gap else None,
               score_mean_top_in_site=rnd(sc[ins].mean()) if ins.any() else None,
               score_mean_top_outside=rnd(sc[~ins].mean()) if (~ins).any() else None,
               wall_s_total=rnd(sum(r['wall_s'] for r in per), 1), per_seed=per,
               top_poses_outside_site=outside)
    if same_mol_ref is not None:
        rm = np.array([r['rmsd_top_to_reference_A'] for r in per])
        out['rmsd_top_to_reference'] = dict(reference=ref_label, mean=rnd(rm.mean()),
                                            sd=rnd(rm.std(ddof=1)) if len(rm) > 1 else 0.0,
                                            min=rnd(rm.min()), max=rnd(rm.max()),
                                            f_le_1A=rnd((rm <= 1.0).mean()), f_le_2A=rnd((rm <= 2.0).mean()))
    return out


def analyse(do_posebusters=True):
    import hashlib
    import shutil
    from scipy import stats as sps
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    _, cross = cross_reference()
    A, PA = patch_pipeline()
    PA.POSEBUSTERS = do_posebusters
    import contacts as CT
    ref = Chem.MolFromMolFile(str(DOCK / C1 / 'ref_RLI.sdf'))
    ref_xyz, _ = mol_xyz(ref)
    done = [j for j in jobs() if (DOCK / j['rundir'] / 'results.json').exists()]
    pending = [j['rundir'] for j in jobs() if j not in done]
    tmp = OUT / '_tmp'
    tmp.mkdir(exist_ok=True)

    # --- 1. validation: the gate-G2 analysis of the redocking jobs --------------------------------
    validation = {}
    for j in done:
        if j['kind'] != 'validation':
            continue
        v = A.analyse_job(DOCK / j['rundir'], do_posebusters=do_posebusters)
        st = v['stats']
        entry = dict(job=j['rundir'], ligand=j['ligand'], scoring=j['scoring'], box=v['box'],
                     reference=v['reference'], stats=st,
                     site_checks_on_crystal_reference=v.get('site_checks_on_crystal_reference'),
                     per_seed=[{k: r[k] for k in ('seed', 'score_top', 'rmsd_top', 'rmsd_best_of_modes',
                                                  'rank_of_best_rmsd') if k in r}
                               for r in v['per_seed'] if r.get('ok')])
        if st.get('rmsd_available'):
            entry['gate'] = gate(st, v.get('site_checks_on_crystal_reference') if j['ligand'] == 'rli_anion'
                                 else None)
        validation[j['rundir']] = entry
        print(f"validation {j['rundir']:46s} RMSD {st.get('rmsd_top_mean', float('nan')):6.2f} +- "
              f"{st.get('rmsd_top_sd', float('nan')):5.2f}  f1 {st.get('success_fraction_1A')}  "
              f"score {st['score_mean']:7.3f}", flush=True)

    # --- 2. every focused job: scores, clusters, RMSD to the crystal ligands, pose files, contacts -
    focused = {}
    for j in done:
        if j['kind'] == 'whole':
            continue
        q = PA.analyse_job(DOCK / j['rundir'], {}, [], diagnostic=False)
        q['kind'] = j['kind']
        q['box_edge_A'] = j['box'] or BOX
        q['site'] = site_geometry(A, j['receptor'], j['rundir'], ref_xyz)
        res = json.loads((DOCK / j['rundir'] / 'results.json').read_text())
        by_seed = {s_['seed']: s_ for s_ in res['seeds'] if s_.get('ok')}
        cl = q['clusters']['clusters'][0]
        pose = CT.top_pose(DOCK / by_seed[cl['representative_seed']]['out'])
        prof = CT.profile_pose(pose, j['receptor'], DOCK / j['receptor'] / 'receptor.pdb', tmp)
        q['representative'] = dict(cluster=1, population=cl['population'], seed=cl['representative_seed'],
                                   score=by_seed[cl['representative_seed']]['best'],
                                   contacts=CT.key_contacts(prof),
                                   residues_within_4A=prof['residues_within_4A'],
                                   polar_pairs_NO_within_3p5A=prof['polar_pairs_NO_within_3.5A'])
        focused[j['rundir']] = q
        print(f"focused    {j['rundir']:46s} {q['score']['mean']:7.3f} +- {q['score']['sd']:5.3f}  "
              f"clusters {[c['population'] for c in q['clusters']['clusters']]}  in site "
              f"{q['site']['in_site_fraction']}  Asp73 {q['site']['asp73_distance_mean_A']}", flush=True)

    # --- 3. the score table of the control receptor ------------------------------------------------
    table = {}
    TABLES = (('vina', C1, 'vina'), ('vinardo', C1, 'vinardo'), ('vina_water', C1_WAT, 'vina'))
    for tkey, trec, sc in TABLES:
        rows = [q for k, q in focused.items() if q['target'] == trec and q['scoring'] == sc]
        if not rows:
            continue
        rows.sort(key=lambda q: q['score']['mean'])
        ref_q = next((q for q in rows if q['ligand'] == 'rli_anion'), None)
        ranked = []
        for k, q in enumerate(rows):
            ranked.append(dict(rank=k + 1, species=q['ligand'], name=NAMES[q['ligand']],
                               heavy_atoms=q['heavy_atoms'], charge=q['charge'],
                               score_mean=q['score']['mean'], score_sd=q['score']['sd'],
                               LE=q['ligand_efficiency']['LE_mean'], FQ=q['ligand_efficiency']['FQ_mean'],
                               geometry=q['geometry'], box_edge_A=q['box_edge_A'],
                               n_clusters=q['clusters']['n_clusters'],
                               cluster1_population=q['clusters']['clusters'][0]['population'],
                               in_site_fraction=q['site']['in_site_fraction'],
                               delta_vs_crystal_ligand=None if ref_q is None else
                               rnd(q['score']['mean'] - ref_q['score']['mean']),
                               nearest_residues=list(q['representative']['residues_within_4A'])[:4]))
        # size relation over the species docked into this receptor (all in the 20 A box)
        std = [r for r in ranked if r['box_edge_A'] == BOX]
        assert len(std) == len(ranked)
        ha = np.array([r['heavy_atoms'] for r in std], float)
        sm = np.array([r['score_mean'] for r in std])
        le = np.array([r['LE'] for r in std])
        fit = np.polyfit(ha, sm, 1) if len(std) > 2 else (None, None)
        table[tkey] = dict(receptor=trec, scoring=sc, ranking=ranked, n_species=len(ranked),
                           size_relation=dict(
                               n=len(std), heavy_atom_range=[int(ha.min()), int(ha.max())],
                               r_score_vs_HA=rnd(np.corrcoef(ha, sm)[0, 1]) if len(std) > 2 else None,
                               r_LE_vs_HA=rnd(np.corrcoef(ha, le)[0, 1]) if len(std) > 2 else None,
                               slope_kcal_per_heavy_atom=rnd(fit[0], 4) if len(std) > 2 else None,
                               intercept_kcal=rnd(fit[1], 4) if len(std) > 2 else None,
                               note='Pearson r and least-squares line of the mean score against the number of '
                                    f'heavy atoms, over the species docked into this receptor ({BOX:g} A box)'))

    # --- 4. the same species in the receptors they do bind (production_summary.json, read only) ----
    prod = json.loads((DOCK / 'results' / 'production_summary.json').read_text())['jobs']
    compare = {}
    for tkey, trec, sc in TABLES:
        c1 = {r['species']: r['score_mean'] for r in table.get(tkey, {}).get('ranking', [])}
        rows = {}
        for sp in SFA + CONTROLS:
            if sp not in c1:
                continue
            row = dict(name=NAMES[sp], control_4DUH=c1[sp])
            for t, pdb in TARGET_PDB.items():
                jb = prod.get(f'{t}/{sp}__{sc}__prod')
                if jb and 'score' in jb:
                    row[pdb] = jb['score']['mean']
                    row[f'control_minus_{pdb}'] = rnd(c1[sp] - jb['score']['mean'])
            rows[sp] = row
        entry = dict(scores=rows)
        for pdb in TARGET_PDB.values():
            both = [(r['control_4DUH'], r[pdb]) for r in rows.values() if pdb in r]
            if len(both) >= 4:
                a_, b_ = np.array(both).T
                d_ = a_ - b_
                rho = sps.spearmanr(a_, b_)
                entry[f'vs_{pdb}'] = dict(
                    n_species=len(both), mean_difference=rnd(d_.mean()), min_difference=rnd(d_.min()),
                    max_difference=rnd(d_.max()), pearson_r=rnd(np.corrcoef(a_, b_)[0, 1]),
                    spearman_rho=rnd(rho.statistic), spearman_p=rnd(rho.pvalue, 4),
                    note='difference = control receptor minus target; positive = weaker score in the '
                         'control receptor')
        # the size relation like for like: the same species in every receptor
        common = [sp for sp, r in rows.items() if all(pdb in r for pdb in TARGET_PDB.values())]
        if len(common) >= 4:
            hv = {r['species']: r['heavy_atoms'] for r in table[tkey]['ranking']}
            ha_ = np.array([hv[sp] for sp in common], float)
            entry['size_relation_common_species'] = dict(
                n_species=len(common), species=common,
                r_score_vs_HA=dict({'4DUH': rnd(np.corrcoef(ha_, [rows[sp]['control_4DUH'] for sp in common])[0, 1])},
                                   **{pdb: rnd(np.corrcoef(ha_, [rows[sp][pdb] for sp in common])[0, 1])
                                      for pdb in TARGET_PDB.values()}),
                note='Pearson r between the mean score and the number of heavy atoms over the species '
                     'docked into the control receptor and into all four production structures')
        compare[tkey] = entry
    # neutral sulfanilamide over every variant of the control receptor against the production structures (Vina)
    sv = {t: focused[f'{t}/runs/sfa_neutral__vina__prod']['score']['mean'] for t in (C1, C1_ALT, C1_NOLOOP, C1_WAT)
          if f'{t}/runs/sfa_neutral__vina__prod' in focused}
    st_ = {pdb: prod[f'{t}/sfa_neutral__vina__prod']['score']['mean'] for t, pdb in TARGET_PDB.items()
           if f'{t}/sfa_neutral__vina__prod' in prod}
    compare['sulfanilamide_vina_all_variants'] = dict(
        control=sv, production=st_, control_range=[min(sv.values()), max(sv.values())],
        max_abs_difference=rnd(max(abs(a_ - b_) for a_ in sv.values() for b_ in st_.values())),
        note='mean top-pose Vina scores of neutral sulfanilamide; the largest absolute difference between any '
             'variant of the control receptor and any production structure')

    # --- 5. the second loop build ------------------------------------------------------------------
    loop = dict(receptors={}, jobs={})
    for t in (C1, C1_ALT):
        rep = json.loads((DOCK / t / 'prep_report.json').read_text())
        rows = heavy_rows(t)
        L = np.array([r[3] for r in rows if r[0] in LOOP])
        centre = np.array(json.loads((DOCK / t / 'box.json').read_text())['center'])
        loop['receptors'][t] = dict(
            pdbfixer_seed=rep['pdbfixer_seed'], rebuilt_residues=sorted(LOOP),
            n_loop_heavy_atoms=int(len(L)),
            loop_to_crystal_inhibitor_min_A=rnd(np.linalg.norm(L[:, None] - ref_xyz[None], axis=2).min(), 2),
            loop_heavy_atoms_in_box=int((np.abs(L - centre) <= BOX / 2).all(axis=1).sum()),
            restrained_atoms_max_displacement_A=rep['minimisation']['max_displacement_restrained_A'],
            restrained_atoms_rms_displacement_A=rep['minimisation']['rms_displacement_restrained_A'],
            min_heavy_heavy_distance_A=rep['checks']['min_heavy_heavy_distance_A'],
            receptor_pdbqt_sha256=hashlib.sha256((DOCK / t / 'receptor.pdbqt').read_bytes()).hexdigest(),
            receptor_pdb_sha256=hashlib.sha256((DOCK / t / 'receptor.pdb').read_bytes()).hexdigest(),
            protonation_non_default=[d['residue'] + '=' + d['variant'] for d in rep['protonation_variants']
                                     if d['variant'] != 'HIE'],
            n_histidines_HIE_by_default=sum(1 for d in rep['protonation_variants'] if d['variant'] == 'HIE'))
    a_rows, b_rows = heavy_rows(C1), heavy_rows(C1_ALT)
    bmap = {(r[0], r[2]): r[3] for r in b_rows}
    dl = np.array([np.linalg.norm(r[3] - bmap[(r[0], r[2])]) for r in a_rows if r[0] in LOOP])
    do = np.array([np.linalg.norm(r[3] - bmap[(r[0], r[2])]) for r in a_rows if r[0] not in LOOP])
    loop['build_difference'] = dict(loop_heavy_atoms_rms_A=rnd(np.sqrt((dl ** 2).mean()), 2),
                                    loop_heavy_atoms_max_A=rnd(dl.max(), 2),
                                    other_heavy_atoms_rms_A=rnd(np.sqrt((do ** 2).mean()), 3),
                                    other_heavy_atoms_max_A=rnd(do.max(), 2))
    if (DOCK / C1_NOLOOP / 'prep_report.json').exists():
        loop['receptors'][C1_NOLOOP] = json.loads((DOCK / C1_NOLOOP / 'prep_report.json').read_text())
    for tag, other in (('second_build', C1_ALT), ('no_loop', C1_NOLOOP)):
        for k, q in focused.items():
            if q['target'] != other:
                continue
            twin = focused.get(k.replace(other, C1))
            if twin is None:
                continue
            row = loop['jobs'].setdefault(q['ligand'], dict(
                ligand=q['ligand'], score_primary=twin['score']['mean'], score_sd_primary=twin['score']['sd'],
                loop_contact_fraction_primary=twin['site']['loop_contact_fraction_4A'],
                loop_min_distance_primary_A=twin['site']['loop_min_distance_A'],
                in_site_fraction_primary=twin['site']['in_site_fraction']))
            row.update({f'score_{tag}': q['score']['mean'], f'score_sd_{tag}': q['score']['sd'],
                        f'score_difference_{tag}': rnd(q['score']['mean'] - twin['score']['mean']),
                        f'in_site_fraction_{tag}': q['site']['in_site_fraction'],
                        f'n_clusters_{tag}': q['clusters']['n_clusters']})
            # is the top pose the same pose as in the primary receptor?  (same molecule, same frame)
            ra = json.loads((DOCK / twin['target'] / 'runs' / twin['tag'] / 'results.json').read_text())
            rb = json.loads((DOCK / q['target'] / 'runs' / q['tag'] / 'results.json').read_text())
            pa = CT.top_pose(DOCK / next(s_['out'] for s_ in ra['seeds'] if s_['seed'] ==
                                         twin['clusters']['clusters'][0]['representative_seed']))
            pb_ = A.load_poses(DOCK / next(s_['out'] for s_ in rb['seeds'] if s_['seed'] ==
                                           q['clusters']['clusters'][0]['representative_seed']))
            row[f'representative_rmsd_to_primary_{tag}_A'] = rnd(
                A.symm_rmsd(pa, pb_, [pb_.GetConformers()[0].GetId()])[0])
            v1, v2 = validation.get(k.replace(other, C1)), validation.get(k)
            if v1 and v2 and v1['stats'].get('rmsd_available'):
                row.update({'rmsd_primary': v1['stats']['rmsd_top_mean'],
                            f'rmsd_{tag}': v2['stats']['rmsd_top_mean'],
                            f'rmsd_sd_{tag}': v2['stats']['rmsd_top_sd'],
                            'f1A_primary': v1['stats']['success_fraction_1A'],
                            f'f1A_{tag}': v2['stats']['success_fraction_1A'],
                            f'f2A_{tag}': v2['stats']['success_fraction_2A']})
    dl_ = [q['site']['loop_min_distance_A'] for q in focused.values()
           if q['target'] == C1 and q['site']['loop_min_distance_A'] is not None]
    loop['closest_pose_to_loop_A'] = rnd(min(dl_), 2) if dl_ else None
    loop['poses_touching_loop'] = sorted(k for k, q in focused.items() if q['target'] == C1
                                         and (q['site']['loop_contact_fraction_4A'] or 0) > 0)

    # --- 5b. the near-native mode of the inhibitor, and the ordered water of the crystal ------------
    wat = np.array([w[2] for w in chain_waters()])
    dref = np.linalg.norm(wat[:, None] - ref_xyz[None], axis=2).min(1)
    water = dict(clash_cutoff_A=WATER_CLASH, n_chain_waters=int(len(wat)),
                 crystal_inhibitor=dict(waters_overlapped=int((dref < WATER_CLASH).sum()),
                                        nearest_water_A=rnd(dref.min(), 2),
                                        waters_within_3p3A=int((dref <= 3.3).sum())),
                 near_native_mode={}, jobs={})
    for k, v in validation.items():
        if v['ligand'] != 'rli_anion':
            continue
        res = json.loads((DOCK / k / 'results.json').read_text())
        rows = []
        for s_ in res['seeds']:
            pm = A.load_poses(DOCK / s_['out'])
            cids = [c.GetId() for c in pm.GetConformers()]
            rm = A.symm_rmsd(ref, pm, cids)
            first = next((i for i, x in enumerate(rm) if x <= 1.0), None)
            row = dict(seed=s_['seed'], rmsd_top_A=rnd(rm[0]), score_top=s_['scores'][0],
                       waters_overlapped_top=int((np.linalg.norm(
                           wat[:, None] - mol_xyz(pm, cids[0])[0][None], axis=2).min(1) < WATER_CLASH).sum()))
            if first is not None:
                row.update(rank_near_native=first + 1, rmsd_near_native_A=rnd(rm[first]),
                           score_near_native=s_['scores'][first],
                           gap_near_native_minus_top=rnd(s_['scores'][first] - s_['scores'][0]),
                           waters_overlapped_near_native=int((np.linalg.norm(
                               wat[:, None] - mol_xyz(pm, cids[first])[0][None], axis=2).min(1)
                               < WATER_CLASH).sum()))
            rows.append(row)
        nn = [r for r in rows if 'rank_near_native' in r]
        water['near_native_mode'][k] = dict(
            n_seeds=len(rows), seeds_with_a_mode_within_1A=len(nn),
            rank_min=min(r['rank_near_native'] for r in nn) if nn else None,
            rank_max=max(r['rank_near_native'] for r in nn) if nn else None,
            rmsd_min_A=min(r['rmsd_near_native_A'] for r in nn) if nn else None,
            rmsd_max_A=max(r['rmsd_near_native_A'] for r in nn) if nn else None,
            score_mean=rnd(np.mean([r['score_near_native'] for r in nn])) if nn else None,
            gap_min=min(r['gap_near_native_minus_top'] for r in nn) if nn else None,
            gap_max=max(r['gap_near_native_minus_top'] for r in nn) if nn else None,
            gap_mean=rnd(np.mean([r['gap_near_native_minus_top'] for r in nn])) if nn else None,
            waters_overlapped_top=[min(r['waters_overlapped_top'] for r in rows),
                                   max(r['waters_overlapped_top'] for r in rows)],
            waters_overlapped_near_native=[min(r['waters_overlapped_near_native'] for r in nn),
                                           max(r['waters_overlapped_near_native'] for r in nn)] if nn else None,
            per_seed=rows)
    if (DOCK / C1_WAT / 'prep_report.json').exists():
        water['receptor'] = json.loads((DOCK / C1_WAT / 'prep_report.json').read_text())
        for k, q in focused.items():
            if q['target'] != C1_WAT:
                continue
            twin = focused.get(k.replace(C1_WAT, C1))
            if twin is None:
                continue
            row = dict(ligand=q['ligand'], score_without_water=twin['score']['mean'],
                       score_with_water=q['score']['mean'], score_sd_with_water=q['score']['sd'],
                       score_difference=rnd(q['score']['mean'] - twin['score']['mean']),
                       in_site_fraction_without_water=twin['site']['in_site_fraction'],
                       in_site_fraction_with_water=q['site']['in_site_fraction'],
                       waters_overlapped_without_water=twin['site']['crystal_waters_overlapped_mean'])
            v1, v2 = validation.get(k.replace(C1_WAT, C1)), validation.get(k)
            if v1 and v2 and v2['stats'].get('rmsd_available'):
                row.update(rmsd_without_water=v1['stats']['rmsd_top_mean'], rmsd_with_water=v2['stats']['rmsd_top_mean'],
                           rmsd_sd_with_water=v2['stats']['rmsd_top_sd'],
                           f1A_without_water=v1['stats']['success_fraction_1A'],
                           f1A_with_water=v2['stats']['success_fraction_1A'],
                           f2A_with_water=v2['stats']['success_fraction_2A'],
                           gate_with_water=v2.get('gate'))
            water['jobs'][q['ligand']] = row

    # --- 6. whole-fragment search ------------------------------------------------------------------
    whole = {}
    for j in done:
        if j['kind'] != 'whole':
            continue
        if j['ligand'] == 'rli_anion':
            same, label = ref, 'crystal ligand of 4DUH'
        else:
            fq = focused.get(f'{C1}/runs/{j["ligand"]}__{j["scoring"]}__prod')
            same, label = None, None
            if fq:
                res = json.loads((DOCK / fq['target'] / 'runs' / fq['tag'] / 'results.json').read_text())
                seed = fq['clusters']['clusters'][0]['representative_seed']
                same = CT.top_pose(DOCK / next(s_['out'] for s_ in res['seeds'] if s_['seed'] == seed))
                label = f'representative pose of the focused search ({fq["tag"]}, seed {seed})'
        whole[j['rundir']] = whole_fragment(A, j, ref_xyz, same, label)
        w = whole[j['rundir']]
        print(f"whole      {j['rundir']:46s} {w['score_mean']:7.3f} +- {w['score_sd']:5.3f}  top pose in "
              f"site {w['n_top_in_site']}/{w['n_seeds']}", flush=True)

    summary = dict(
        generated=time.strftime('%Y-%m-%d %H:%M %Z'),
        generated_by='python dock/control_gyrb.py analyse',
        purpose='negative-control receptor: the ATP site of E. coli DNA gyrase B (4DUH), which is not a '
                'target of sulfanilamide; same protocol as the five receptors T1-T5',
        protocol=dict(program='AutoDock Vina 1.2.7', scoring=['vina', 'vinardo'], box_edge_A=BOX,
                      exhaustiveness=EXH, seeds=SEEDS,
                      whole_fragment_box=json.loads((DOCK / C1 / 'box_whole.json').read_text()),
                      geometry='sulfanilamide species from their B3LYP-D3(BJ)/IEF-PCM geometries; every '
                               'other ligand from one ETKDG/MMFF94s conformer per seed',
                      site_definition=f'pose centroid within {SITE_CUT} A of a heavy atom of the crystal '
                                      'inhibitor',
                      scores='empirical docking scores, kcal/mol by convention; never converted to Ki'),
        target=INFO, site_checks=[dict(label=s_[0], smarts=s_[1], receptor=s_[2], cutoff_A=s_[3]) for s_ in SITE],
        clorobiocin_reference=cross, validation=validation, jobs=focused, score_table=table,
        comparison_with_targets=compare, loop_sensitivity=loop, water_diagnostic=water,
        whole_fragment=whole,
        pending_jobs=pending, wall_s=rnd(time.time() - t0, 1))
    (OUT / 'control_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    shutil.rmtree(tmp, ignore_errors=True)
    print(f'wrote {OUT / "control_summary.json"} ({len(focused)} focused jobs, {len(whole)} whole-fragment '
          f'jobs, {len(pending)} pending) in {time.time() - t0:.0f} s')
    return 0


# --------------------------------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------------------------------
def report():
    """dock/CONTROL.md, written from control_summary.json alone: every number in it is read from that
    file, and every sentence that states an outcome is chosen by the numbers, not typed in advance."""
    S = json.loads((OUT / 'control_summary.json').read_text())
    V, Jb, W = S['validation'], S['jobs'], S['water_diagnostic']
    nn, loop, cb = W['near_native_mode'], S['loop_sensitivity'], S['clorobiocin_reference']
    ov = cb['overlap'][C1]
    one = sorted(k for k, q in Jb.items() if q.get('cpu') == 1)
    f1 = lambda x: f'{x:.1f}'  # noqa: E731
    f2 = lambda x: f'{x:.2f}'  # noqa: E731
    f3 = lambda x: f'{x:.3f}'  # noqa: E731
    pm = lambda q: f"{q['score']['mean']:.3f} ± {q['score']['sd']:.3f}"  # noqa: E731
    rname = {C1: 'as prepared (residues 109–116 built, PDBFixer seed 20260928)',
             C1_ALT: 'residues 109–116 built from a second seed (12345)',
             C1_NOLOOP: 'residues 109–116 removed (observed residues only)',
             C1_WAT: 'as prepared, with the ordered water of chain A put back (diagnostic)'}
    L = ['# Control receptor: the ATP site of *E. coli* DNA gyrase B (4DUH)', '',
         f"Generated {S['generated']} by `python dock/control_gyrb.py analyse` and `report`. Every number in this file is "
         'read from `dock/control/control_summary.json`.', '',
         '## 1. Why a control receptor, and why this one', '',
         'The five receptors of the study (T1–T5, `dock/VALIDATION.md` and `dock/RESULTS.md`) all belong to sulfonamide '
         'pharmacology: DHPS, hCA II and Sul2, the DHPS variant that rejects sulfonamides. A docking score says something '
         'about a ligand only if it changes when the receptor is not a target of the ligand, and nothing in gate G2 or in the '
         'production run measures that. The control is the ATP site of DNA gyrase B of '
         '*Escherichia coli*, in the 24 kDa N-terminal fragment crystallised with a 4,5′-bithiazole inhibitor (PDB 4DUH, '
         f"chain A, {json.loads((DOCK / 'targets_inventory.json').read_text())[C1]['resolution']:.2f} Å; Brvar et al. 2012, "
         'doi:10.1021/jm300395d). This site binds ATP and the aminocoumarins. The antibacterial target of sulfanilamide is '
         'dihydropteroate synthase, and the fluoroquinolones act on gyrase at the GyrA–DNA cleavage complex, not here; '
         'neither sulfanilamide nor ciprofloxacin is a known ligand of this pocket. No measurement says that sulfanilamide '
         'does not bind it, and a fragment of 11 heavy atoms may well bind weakly, so the control is described as a '
         'receptor that is not a target, never as one that the ligand does not bind.', '',
         'The control is kept apart from the five receptors by construction: its directories are named `C1_*`, which the '
         '`T*` globs of `analyse.py`, `production_analysis.py` and `contacts.py` never match, and its output is in '
         '`dock/control/`. Gate G2, the production summary and every statistic of the five receptors are unchanged.', '',
         '## 2. Receptors, ligands and protocol', '',
         'Protocol of the five receptors: AutoDock Vina 1.2.7, Vina and Vinardo scoring, a 20 Å cube on the heavy-atom '
         'centroid of the crystal ligand, exhaustiveness 32, 10 seeds, `--cpu 3`; receptor from one chain without water, '
         'PROPKA protonation at pH 7.4, restrained minimisation of the built atoms and hydrogens. The three sulfanilamide '
         'species are docked from their B3LYP-D3(BJ)/IEF-PCM geometries, every other ligand from one ETKDG/MMFF94s conformer '
         'per seed. New ligands (`dock/prep_ligands.py`, SMILES checked against the RCSB chemical components): the 4DUH '
         'inhibitor as its benzoate (`rli_anion`, component RLI) and ciprofloxacin as its zwitterion (`cfx_zwit`, CPF). '
         'Clorobiocin is not docked (section 3.3).', '',
         '| Receptor directory | What it is | Rebuilt residues to the crystal inhibitor (Å) | Rebuilt heavy atoms in the box |',
         '|---|---|---:|---:|']
    for t in (C1, C1_ALT, C1_NOLOOP, C1_WAT):
        r_ = loop['receptors'].get(t, {})
        L.append(f"| `{t}` | {rname[t]} | "
                 f"{f2(r_['loop_to_crystal_inhibitor_min_A']) if 'loop_to_crystal_inhibitor_min_A' in r_ else ('as `C1_4DUH`' if t == C1_WAT else '–')} | "
                 f"{r_.get('loop_heavy_atoms_in_box', 'as `C1_4DUH`' if t == C1_WAT else '–')} |")
    bd = loop['build_difference']
    wr = W.get('receptor', {})
    L += ['',
          f"Chain A of 4DUH lacks residues 109–116 of the ATP lid. PDBFixer builds them, as it builds the internal gaps of the "
          f"five receptors; a loop built without density is a guess, so it was built twice. The two builds differ by "
          f"{f2(bd['loop_heavy_atoms_rms_A'])} Å RMS in the built residues (largest {f2(bd['loop_heavy_atoms_max_A'])} Å) and by "
          f"{f3(bd['other_heavy_atoms_rms_A'])} Å in every other heavy atom. The deposited atoms moved by at most "
          f"{loop['receptors'][C1]['restrained_atoms_max_displacement_A']} Å in the restrained minimisation. PROPKA leaves every "
          f"Asp and Glu charged; non-default protonation: {', '.join(loop['receptors'][C1]['protonation_non_default']) or 'none'}. "
          f"`{C1_NOLOOP}` is `{C1}` with the {loop['receptors'].get(C1_NOLOOP, {}).get('n_atoms_removed', '–')} atoms of the built "
          f"residues deleted. `{C1_WAT}` is `{C1}` with {wr.get('n_waters_kept', '–')} of the {wr.get('n_waters_in_file', '–')} "
          f"waters of chain A added as oxygen atoms ({wr.get('n_waters_in_box', '–')} inside the box; one water that overlaps a "
          'built atom is left out). The water is fixed and has no hydrogens: each oxygen is an acceptor-type atom that '
          'fills space. It is a diagnostic, not the protocol of the five receptors, which were docked without such water; '
          'its scores are given beside those of the receptor as prepared and are marked as such.', '',
          f"{len(one)} of the {len(Jb)} focused jobs ran with `--cpu 1` on the spare core beside the queue "
          f"({', '.join('`' + k.replace('/runs/', ': ') + '`' for k in one)}); the others ran with `--cpu 3`. The result does not "
          f"depend on `--cpu`: the output file of seed 1 of `{C1}/runs/sfa_neutral__vina__prod` was reproduced byte for byte at "
          '`--cpu 1` (sha256 identical), as `dock/VALIDATION.md` found for `--cpu 2`.', '',
          '## 3. Redocking of the crystal inhibitor', '',
          'Same analysis as gate G2 (`analyse.analyse_job`): symmetry-corrected heavy-atom RMSD of the top pose in the receptor '
          'frame, PoseBusters, and the three polar contacts of the deposited pose as site checks. The near-native mode is the '
          'first mode of a seed within 1.0 Å of the crystal pose.', '',
          '| Receptor | Scoring | RMSD (Å) | f ≤ 1 Å | f ≤ 2 Å | Score | Near-native mode: seeds | rank | its RMSD (Å) | its score − top | Waters under the top pose | 1.0 Å gate | 2.0 Å |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|']
    for k, v in V.items():
        if v['ligand'] != 'rli_anion':
            continue
        st, g, n_ = v['stats'], v['gate'], nn[k]
        t = k.split('/runs/')[0]
        has = n_['seeds_with_a_mode_within_1A'] > 0
        L.append(f"| `{t}` | {v['scoring']} | {st['rmsd_top_mean']:.3f} ± {st['rmsd_top_sd']:.3f} | {st['success_fraction_1A']} | "
                 f"{st['success_fraction_2A']} | {st['score_mean']:.3f} ± {st['score_sd']:.3f} | {n_['seeds_with_a_mode_within_1A']}/{n_['n_seeds']} | "
                 + (f"{n_['rank_min']}–{n_['rank_max']} | {n_['rmsd_min_A']:.2f}–{n_['rmsd_max_A']:.2f} | {n_['gap_min']:.3f} to {n_['gap_max']:.3f}"
                    if has else '– | – | –')
                 + f" | {n_['waters_overlapped_top'][0]}–{n_['waters_overlapped_top'][1]} | {'PASS' if g['passes_strict_1A'] else 'fail'} | "
                 f"{'pass' if g['passes_2A_criterion'] else 'fail'} |")
    kv = f'{C1}/runs/rli_anion__vina'
    if kv in V:
        st, n_ = V[kv]['stats'], nn[kv]
        sc_ = st['site_checks']
        xr = {c['check']: c['distance'] for c in V[kv]['site_checks_on_crystal_reference']}
        qv = Jb[kv]
        L += ['',
              f"**The receptor {'passes' if V[kv]['gate']['passes_strict_1A'] else 'fails'} the gate.** With Vina the top pose of every "
              f"seed lies {st['rmsd_top_min']:.2f}–{st['rmsd_top_max']:.2f} Å from the crystal pose, one cluster "
              f"({qv['clusters']['clusters'][0]['population']}/{qv['n_seeds']} seeds), and every top pose passes the PoseBusters physics "
              f"checks ({st['posebusters_physics_pass_fraction']}). The search is not at fault: {n_['seeds_with_a_mode_within_1A']} of "
              f"{n_['n_seeds']} seeds hold a mode within {n_['rmsd_min_A']:.2f}–{n_['rmsd_max_A']:.2f} Å of the crystal pose, at rank "
              f"{n_['rank_min']}–{n_['rank_max']}, {n_['gap_min']:.3f}–{n_['gap_max']:.3f} kcal/mol behind the top pose "
              f"(mean {n_['gap_mean']:.3f}).", '',
              '| Contact of the deposited pose | Crystal (Å) | Top poses, mean (Å) | Pass fraction |', '|---|---:|---:|---:|']
        for lab, q in sc_.items():
            L.append(f"| {lab} | {xr[lab]:.2f} | {q['distance_mean']:.2f} | {q['pass_fraction']} |")
        rep = qv['representative']
        L += ['',
              'The top pose keeps the first contact and gives up the other two: the amide N–H stays on the Asp73 carboxylate, '
              'the bithiazole stays in place, and the benzoate leaves Arg136 and turns under the ATP lid. Residues within 4 Å of '
              f"the representative top pose: {', '.join(f'{k_} ({v_:.2f})' for k_, v_ in list(rep['residues_within_4A'].items())[:10])}. "
              f"The top pose overlaps {n_['waters_overlapped_top'][0]}–{n_['waters_overlapped_top'][1]} crystal waters of chain A "
              f"(oxygen within {W['clash_cutoff_A']} Å of a heavy atom of the pose); the near-native mode overlaps "
              f"{n_['waters_overlapped_near_native'][0]}–{n_['waters_overlapped_near_native'][1]}, and the crystal ligand itself "
              f"{W['crystal_inhibitor']['waters_overlapped']} (its nearest water is {W['crystal_inhibitor']['nearest_water_A']} Å away). "
              'Two things can favour the turned pose: neither Vina nor Vinardo has an electrostatic term, so the salt bridge '
              'between the benzoate and the Arg136 guanidinium group counts only as a hydrogen bond, and the protocol removes the water that '
              'fills the space under the lid in the crystal. Only the second was tested (section 3.1), with Vina, and with all '
              'the ordered water of chain A at once; the three waters under the lid were not tested alone.']
    L += ['', '### 3.1 The ordered water put back (diagnostic)', '']
    wj = W.get('jobs', {})
    if 'rli_anion' in wj:
        r_ = wj['rli_anion']
        L += [f"With the ordered water of chain A in the receptor the top pose of the inhibitor lies "
              f"{r_['rmsd_with_water']:.2f} ± {r_['rmsd_sd_with_water']:.2f} Å from the crystal pose "
              f"({r_['f1A_with_water']} of the seeds within 1 Å, {r_['f2A_with_water']} within 2 Å), against "
              f"{r_['rmsd_without_water']:.2f} Å without it; score {r_['score_with_water']:.3f} against {r_['score_without_water']:.3f} "
              'kcal/mol.']
    L += ['', '| Ligand | Score without water | Score with water | Difference | In the inhibitor site, without / with water |', '|---|---:|---:|---:|---|']
    for sp, r_ in wj.items():
        L.append(f"| {NAMES[sp]} | {r_['score_without_water']:.3f} | {r_['score_with_water']:.3f} ± {r_['score_sd_with_water']:.3f} | "
                 f"{r_['score_difference']:+.3f} | {r_['in_site_fraction_without_water']} / {r_['in_site_fraction_with_water']} |")
    L += ['', '### 3.2 The built residues: a second build, and none', '',
          '| Ligand | Primary | Second build | Difference | Residues removed | Difference | Top pose within 4 Å of a built residue (primary) | Representative pose moved, second build / no loop (Å) |',
          '|---|---:|---:|---:|---:|---:|---:|---|']
    for sp, r_ in loop['jobs'].items():
        g_ = lambda k_, fmt='{:.3f}': fmt.format(r_[k_]) if k_ in r_ else '–'  # noqa: E731
        L.append(f"| {NAMES[sp]} | {r_['score_primary']:.3f} ± {r_['score_sd_primary']:.3f} | {g_('score_second_build')} | "
                 f"{g_('score_difference_second_build', '{:+.3f}')} | {g_('score_no_loop')} | {g_('score_difference_no_loop', '{:+.3f}')} | "
                 f"{r_['loop_contact_fraction_primary']} | {g_('representative_rmsd_to_primary_second_build_A', '{:.2f}')} / "
                 f"{g_('representative_rmsd_to_primary_no_loop_A', '{:.2f}')} |")
    L += ['', f"Jobs of the primary receptor whose top pose lies within 4 Å of a built residue: "
              f"{', '.join('`' + k.split('/runs/')[1] + '`' for k in loop['poses_touching_loop']) or 'none'}; the closest approach of "
              f"any top pose to a built heavy atom is {loop['closest_pose_to_loop_A']} Å.",
          '', '### 3.3 Clorobiocin: no cross-docking test', '',
          f"Clorobiocin is crystallised with the same fragment (1KZN, 2.3 Å; Lafitte et al. 2002, doi:10.1021/bi0159837). 1KZN is "
          'not numbered like *E. coli* GyrB between residues 82 and 96, so its pose is carried into the 4DUH frame by a '
          f"sequence-aligned superposition of {cb['n_ca_atoms']} Cα atoms (RMSD {cb['superposition_ca_rmsd_A']} Å; "
          f"{cb['trimmed_fit']['ca_rmsd_A']} Å over {cb['trimmed_fit']['n_ca_atoms']} atoms after trimming). In 4DUH that pose does not "
          f"fit: its closest heavy-atom contact with the prepared receptor is {ov['min_heavy_heavy_distance_A']} Å "
          f"({ov['closest_receptor_atom']}), with {ov['n_pairs_below_2p5A']} pairs below {cb['pair_cutoff_A']} Å "
          f"({', '.join(f'{k_} {v_}' for k_, v_ in ov['pairs_below_2p5A_by_residue'].items())}) and "
          f"{ov['n_pairs_below_3p0A']} below 3.0 Å, in residues {', '.join(ov['residues_within_3A'])}. "
          f"{ov['n_pairs_below_2p5A_without_counterpart_in_1KZN']} of the {ov['n_pairs_below_2p5A']} pairs are with "
          f"{', '.join(ov['residues_below_2p5A_without_counterpart_in_1KZN'])}, residues that are ordered in 4DUH and have no "
          'counterpart in the 1KZN model (residues of 4DUH without one: '
          f"{', '.join(str(x) for x in cb['residues_of_4DUH_without_counterpart_in_1KZN'])}); the other pairs are with residues "
          'ordered in both crystals and show that the two structures also differ locally. A cross-docking test of 4DUH with '
          'clorobiocin could therefore not pass, and clorobiocin is not docked: a score of a pose that the receptor cannot hold '
          'would say nothing.']
    L += ['', '## 4. Scores of every species in the control receptor', '',
          f"A pose is *in the cavity of the inhibitor* when its heavy-atom centroid lies within {SITE_CUT} Å of a heavy atom of "
          'the crystal inhibitor. This says where a pose is, not that it is right: the top pose of the inhibitor itself in the '
          'receptor as prepared counts in every seed although it lies 6.7 Å from the crystal pose. Scores are empirical '
          'docking scores in kcal/mol, mean ± SD over the seeds, never affinities.', '']
    TT = {'vina': f'Vina, `{C1}` (protocol of the targets)', 'vinardo': f'Vinardo, `{C1}`',
          'vina_water': f'Vina, `{C1_WAT}` (ordered crystal water kept; the receptor that redocks the inhibitor)'}
    for tkey, title in TT.items():
        tb = S['score_table'].get(tkey)
        if not tb:
            continue
        sc, trec = tb['scoring'], tb['receptor']
        L += [f"### {title}", '',
              '| Rank | Ligand | Heavy atoms | Score | Best seed | LE | Clusters (seeds in the first) | Seeds in the cavity of the inhibitor | Centroid to the inhibitor (Å) | Seeds with N/O within 3.5 Å of Asp73 | Waters under the pose | To a built residue (Å) | Nearest residues |',
              '|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|']
        for r_ in tb['ranking']:
            q = next(q_ for q_ in Jb.values() if q_['target'] == trec and q_['ligand'] == r_['species'] and q_['scoring'] == sc)
            L.append(f"| {r_['rank']} | {r_['name']} | {r_['heavy_atoms']} | {pm(q)} | "
                     f"{q['score']['min']:.3f} | {r_['LE']:.3f} | {r_['n_clusters']} ({r_['cluster1_population']}) | "
                     f"{round(q['site']['in_site_fraction'] * q['n_seeds'])}/{q['n_seeds']} | {q['site']['centroid_to_inhibitor_mean_A']} | "
                     f"{round(q['site']['asp73_contact_fraction_3p5A'] * q['n_seeds'])}/{q['n_seeds']} | {q['site']['crystal_waters_overlapped_mean']} | "
                     f"{q['site']['loop_min_distance_A']} | {', '.join(r_['nearest_residues'])} |")
        sr = tb['size_relation']
        if sr.get('r_score_vs_HA') is not None:
            L += ['', f"Pearson r between the mean score and the number of heavy atoms over these {sr['n']} species "
                      f"({sr['heavy_atom_range'][0]}–{sr['heavy_atom_range'][1]} heavy atoms): {sr['r_score_vs_HA']}; least-squares "
                      f"line {sr['slope_kcal_per_heavy_atom']} kcal/mol per heavy atom."]
        L.append('')
    L += ['## 5. The same species in the four production structures', '',
          'Production means of `dock/results/production_summary.json` (read only): 3TZF (DHPS), 7S2K (Sul2, which rejects '
          'sulfonamides), 6RL9 and 3HS4 (hCA II; Vina fails the redocking gate there, having no metal term). Difference = '
          'control receptor minus production structure; positive = weaker score in the control receptor.', '']
    for tkey, title in TT.items():
        cmp_ = S['comparison_with_targets'].get(tkey, {})
        if not cmp_.get('scores'):
            continue
        L += [f"### {title}", '', '| Species | 4DUH | 3TZF | 7S2K | 6RL9 | 3HS4 | 4DUH − 3TZF | 4DUH − 7S2K |', '|---|---:|---:|---:|---:|---:|---:|---:|']
        for sp, r_ in cmp_['scores'].items():
            g_ = lambda k_, fmt='{:.3f}': fmt.format(r_[k_]) if k_ in r_ else '–'  # noqa: E731
            L.append(f"| {r_['name']} | {r_['control_4DUH']:.3f} | {g_('3TZF')} | {g_('7S2K')} | {g_('6RL9')} | {g_('3HS4')} | "
                     f"{g_('control_minus_3TZF', '{:+.3f}')} | {g_('control_minus_7S2K', '{:+.3f}')} |")
        L.append('')
        for pdb in TARGET_PDB.values():
            q = cmp_.get(f'vs_{pdb}')
            if q:
                L.append(f"- against {pdb} ({q['n_species']} species): mean difference {q['mean_difference']:+.3f} kcal/mol "
                         f"({q['min_difference']:+.3f} to {q['max_difference']:+.3f}); Pearson r {q['pearson_r']}, Spearman ρ "
                         f"{q['spearman_rho']} (p = {q['spearman_p']}).")
        sz = cmp_.get('size_relation_common_species')
        if sz:
            L += ['', f"Size relation over the {sz['n_species']} species docked into every receptor (Pearson r between mean score and "
                      'heavy atoms): ' + ', '.join(f'{k_} {v_}' for k_, v_ in sz['r_score_vs_HA'].items()) + '.']
        L.append('')
    av = S['comparison_with_targets'].get('sulfanilamide_vina_all_variants')
    if av:
        L += ['Neutral sulfanilamide with Vina, every variant of the control receptor: '
              + ', '.join(f'`{k_}` {v_:.3f}' for k_, v_ in av['control'].items()) + '; production structures: '
              + ', '.join(f'{k_} {v_:.3f}' for k_, v_ in av['production'].items())
              + f". Largest difference between any variant and any production structure: {av['max_abs_difference']:.3f} kcal/mol.", '']
    L += ['## 6. Whole-fragment search', '']
    wb = S['protocol']['whole_fragment_box']
    L += [f"Box {wb['size'][0]} × {wb['size'][1]} × {wb['size'][2]} Å around the whole fragment ({wb['volume_ratio']} times the volume of "
          f"the 20 Å cube), exhaustiveness {wb['exhaustiveness']} ({wb['exhaustiveness_rule']}).", '']
    if S['whole_fragment']:
        L += ['| Job | Seeds | Score | Top pose in the inhibitor site | Seeds with any mode in the site | Top pose to the reference (Å) |', '|---|---:|---:|---:|---:|---|']
        for k, w in S['whole_fragment'].items():
            rr = w.get('rmsd_top_to_reference')
            L.append(f"| `{k}` | {w['n_seeds']} | {w['score_mean']:.3f} ± {w['score_sd']:.3f} | {w['n_top_in_site']}/{w['n_seeds']} | "
                     f"{w['seeds_with_any_in_site_mode']}/{w['n_seeds']} | "
                     + (f"{rr['mean']:.2f} ± {rr['sd']:.2f} (≤ 2 Å in {round(rr['f_le_2A'] * w['n_seeds'])}/{w['n_seeds']}; {rr['reference']})" if rr else '–') + ' |')
            if w['top_poses_outside_site']:
                L += ['', f"Top poses outside the inhibitor site (`{k}`), by seed, with their nearest residues: "
                          + '; '.join(f"seed {s_}: {', '.join(v_)}" for s_, v_ in w['top_poses_outside_site'].items()) + '.']
    else:
        L.append('Not run.')
    L += ['', '## 7. What the control supports', '',
          '- The level of the scores, not the poses. The receptor as prepared does not pass the redocking gate, so no pose in it '
          'is read as a binding mode; the receptor that passes it holds the crystal water and is not the protocol of the five '
          'receptors.',
          '- A pocket that is no target of sulfanilamide can be compared with the production structures only through the same '
          'scoring function, Vina (section 5); the control was not run with AutoDock4Zn, which is specific to the zinc site.',
          '- It does not show that sulfanilamide does not bind gyrase B. It shows that the Vina score does not tell this pocket '
          'from DHPS, Sul2 or hCA II.',
          '- The diagnostic with the crystal water, the two loop builds and the loop-free receptor say how much of each number '
          'depends on choices of the receptor model (section 3).', '',
          '## 8. Reproduce', '',
          '```',
          'python dock/targets.py',
          'python dock/prep_receptor.py C1_4DUH',
          'python dock/prep_receptor.py C1_4DUH --seed 12345 --suffix _seed12345',
          'python dock/prep_ligands.py --only rli_anion cfx_zwit',
          'python dock/control_gyrb.py run        # builds C1_4DUH_wat and C1_4DUH_noloop, then every job (resumes)',
          'python dock/control_gyrb.py analyse    # dock/control/control_summary.json, poses',
          'python dock/control_gyrb.py report     # this file',
          '```', '']
    (DOCK / 'CONTROL.md').write_text('\n'.join(L))
    print(f'wrote {DOCK / "CONTROL.md"} ({len(L)} lines)')
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['run', 'whole-box', 'list', 'analyse', 'report'])
    ap.add_argument('--no-posebusters', action='store_true')
    ap.add_argument('--only', nargs='*', default=None,
                    help='run only the jobs whose run directory contains one of these strings')
    ap.add_argument('--cpu', type=int, default=3)
    a = ap.parse_args()
    if a.stage == 'whole-box':
        print(json.dumps(whole_box(), indent=2))
        return 0
    if a.stage == 'list':
        for j in jobs():
            print(f'{j["rundir"]:58s} seeds {j["seeds"]:2d} {j["geometry"]:5s} {j["kind"]}')
        return 0
    if a.stage == 'analyse':
        return analyse(do_posebusters=not a.no_posebusters)
    if a.stage == 'report':
        return report()
    return run(a.only, a.cpu)


if __name__ == '__main__':
    sys.exit(main())
