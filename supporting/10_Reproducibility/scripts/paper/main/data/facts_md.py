"""Facts of the molecular-dynamics study and of the two docking tests (8 Oct 2026; section 40 of facts.md), executed
by build_facts.py after facts_jics.py (exec in its namespace: add(), section(), FACTS, J, r, REPO, DOCK, DER,
_stt, _DecB and DISAGREE are the ones defined there).

  40. Molecular dynamics of the carbonic anhydrase II-sulfanilamide complex (md/RESULTS.md) and the two docking
      tests in the hCA II receptors (dock/GLU106.md, dock/WATER4.md).  Key prefixes:
        md.sys.*, md.prot.*, md.def.*, md.ff.*, md.zn.*, md.qm.*, md.sw.*   system, protocol, definitions, force field
        md.run.<run>.*                                                      the eleven runs of 100 ns (X1 ... R3,
                                                                            R1_2s, R1_n2fit)
        md.grp.*                                                            ranges and counts over groups of runs
        md.gate.*                                                           the tests of the force field (G-L, G-Z,
                                                                            G-V, G-R, G-R2) and the relaxation of the
                                                                            bare zinc-site model
        md.agree.*, md.pka.*                                                agreement of the readings; pKa of Glu106
        md.xtal.*, md.rule1.*, md.rule3.*, md.rehearsal.*                   values of 6RL9 that the runs are set against;
                                                                            the rules of the plan; the rehearsal run
        dock.glu106.*                                                       Glu106 as the acid (the receptors of the
                                                                            paper) against Glu106 as a carboxylate
        dock.wat4.*                                                         the four crystal waters of 6RL9

Every value is read from a result file by this code; nothing is typed except constants of a definition, which
name the line of code or of the plan that fixes them.  `check=` is the value that derive/d17_md.py forms with code
of its own from data in git: the time series of md/results/series/, the tracked all-atom series, the final
frames, md/qm/gr/gr_table.csv, the raw out.pdbqt files and the deposited entry 6RL9.  For the Zn-N distances, the
nearer sulfonyl oxygen and the C-alpha RMSD of a run the check is the blind third reading of the trajectories
(md/results/independent/third_reading/, code of its own, read through d17).  Facts without a check are read from
their source only; the section text says which these are.

11 Oct 2026, after the independent check of the numbers (paper/audit/research/md_paper_2026-10-08/FACTS_GUIDE.md):
the mean distances of the docking runs take the full-precision mean of the raw poses where the run record, which
averages per-seed values of two decimals, prints differently (_a, full=True; 8 facts, one of them printed: Zn-N of
sulfanilamide from the DFT geometry in 3HS4 with the carboxylate, 2.31 and not 2.30 A); the value of 6RL9 for the
zinc out of the plane of the His119 ring is 0.15 A, formed by d17 from the coordinates (the plan printed 0.17 A,
which is the value of His94); 41 facts were added (pKa of Glu106 in the preparation of both docking receptors,
numbers of the four waters, the dropped water site, charges, the fitting rows of the refit, the first count of the
hydration check and the water in the place of the ring before the change at N2, values of 6RL9, line totals); the
precision of a few facts was set to the way the number is printed.
"""
MDR_ = REPO / 'md'
D17 = J(DER / 'd17_md.json')
D17S = 'derived/d17_md.json'
D17M, D17D = D17['md'], D17['dock']
RUNS17_ = [('bonded_2e', x_, x_) for x_ in ('X1', 'X2', 'X3', 'N1', 'N2', 'N3', 'R1', 'R2', 'R3')] + \
          [('bonded_2s', 'R1', 'R1_2s'), ('bonded_2e_n2fit', 'R1', 'R1_n2fit')]
RLAB_ = {'R1_2s': 'sensitivity run R1 with zinc set 2s', 'R1_n2fit': 'sensitivity run R1 with the candidate terms at N2'}
G_X_, G_N_, G_R_ = ['X1', 'X2', 'X3'], ['N1', 'N2', 'N3'], ['R1', 'R2', 'R3']
NINE_ = G_X_ + G_N_ + G_R_
ELEVEN_ = NINE_ + ['R1_2s', 'R1_n2fit']


def _tie(x_, p_):
    """True when the decimal x_ lies exactly half way between two numbers of p_ decimals."""
    return (abs(x_) * (10 ** p_)) % 1 == _DecB('0.5')


def _dmean(vals_):
    """Exact mean of values that a file stores with a fixed number of decimals."""
    return sum(_DecB(repr(float(v_))) for v_ in vals_) / len(vals_)


def _rq(x_, p_):
    """x_ rounded half to even to p_ decimals, as a Decimal."""
    return (x_ if isinstance(x_, _DecB) else _DecB(repr(float(x_)))).quantize(_DecB(1).scaleb(-p_), rounding=__import__('decimal').ROUND_HALF_EVEN)


def _a(key_, value_, unit_, prec_, source_, meaning_, check=None, tol=None, note=None, exact=False, full=False):
    """add() with the rounding rule of facts.md applied: a stored value that is itself a tie at the precision of the
    fact is replaced by the full-precision value of the re-derivation, with the stored value kept in the note.
    exact=True: the stored decimals are the whole value (a mean of scores that the program prints with three
    decimals), so a tie is a true tie; it is kept as it is and prints rounded half to even.
    full=True (11 Oct 2026): the source is a mean of per-seed values that the run record stores with fewer decimals
    than the raw poses give (distances with two decimals, RMSDs with three).  Where that mean and the mean formed by
    derive/d17_md.py from the raw coordinates print differently at the precision of the fact, the fact takes the
    full-precision mean and keeps the value of the record in the note, so that the printed number is the one the
    deposited poses give."""
    stored_, why_ = None, 'tie'
    if exact and isinstance(value_, _DecB):
        if isinstance(prec_, int) and _tie(value_, prec_):
            note = ((note + '; ') if note else '') + f'an exact tie at {prec_} decimals (mean of scores printed with three decimals); it prints rounded half to even'
        value_ = float(value_)
    elif isinstance(value_, _DecB):
        if check is not None and isinstance(prec_, int) and _tie(value_, prec_):
            stored_ = float(value_)
            value_ = check
        elif full and check is not None and isinstance(prec_, int) and _rq(value_, prec_) != _rq(check, prec_):
            stored_, why_ = float(value_), 'full'
            value_ = check
        else:
            value_ = float(value_)
    elif isinstance(value_, float) and check is not None and isinstance(prec_, int):
        if _tie(_DecB(repr(value_)), prec_):
            stored_, value_ = value_, check
        elif full and _rq(value_, prec_) != _rq(check, prec_):
            stored_, why_, value_ = value_, 'full', check
    if stored_ is not None and why_ == 'tie':
        note = ((note + '; ') if note else '') + (f'the source gives {stored_!r}, a tie at {prec_} decimals; value = the '
                                                 f'full-precision re-derivation ({D17S}), rounded half to even')
    elif stored_ is not None:
        note = ((note + '; ') if note else '') + (f'the run record gives {stored_!r} (a mean of per-seed values stored with fewer decimals), which prints '
                                                 f'{_rq(stored_, prec_)}; value = the mean formed from the raw poses ({D17S}), which prints {_rq(check, prec_)}')
    add(key_, value_, unit_, prec_, source_, meaning_, note=note, check=check, tol=tol)
    if stored_ is not None:
        FACTS[key_]['stored_value'] = stored_


def _rng(base_, ids_, fn_, unit_, prec_, source_, meaning_, checks=None, tol=None):
    """Smallest and largest of fn_ over the runs ids_, as two facts base_.min and base_.max."""
    v_ = [fn_(i_) for i_ in ids_]
    _a(base_ + '.min', min(v_), unit_, prec_, source_, 'smallest ' + meaning_, check=checks[0] if checks else None, tol=tol)
    _a(base_ + '.max', max(v_), unit_, prec_, source_, 'largest ' + meaning_, check=checks[1] if checks else None, tol=tol)


S1_ = {rid_: J(MDR_ / 'results/runs' / m_ / run_ / 'first_program_summary.json') for m_, run_, rid_ in RUNS17_}
SX0_ = {t_: J(MDR_ / 'results/runs' / f'bonded_{t_}' / 'X0' / 'first_program_summary.json') for t_ in ('2e', '2s')}
LAB_ = {rid_: f'{m_}/{run_}' for m_, run_, rid_ in RUNS17_}
HEALTH_ = J(MDR_ / 'results/run_health.json')
HR_ = {x_['run']: x_ for x_ in HEALTH_['runs']}
EV_ = {x_['run']: x_ for x_ in J(MDR_ / 'results/run_events.json')['runs']}
CT_ = {x_['label']: x_ for x_ in J(MDR_ / 'results/conversion_times.json')['runs']}
SS_ = {x_['run']: x_ for x_ in J(MDR_ / 'results/site_checks_solute.json')['runs']}
SF_ = {x_['run']: x_ for x_ in J(MDR_ / 'results/site_checks_full.json')['runs']}
PROT_ = J(MDR_ / 'protocol.json')
BUILD_ = {(m_, s_): J(MDR_ / 'topology' / m_ / s_ / 'build.json')
          for m_, s_ in (('bonded_2e', 'X'), ('bonded_2e', 'N'), ('bonded_2e', 'R'), ('bonded_2s', 'X'))}
DEC_ = J(MDR_ / 'system/decisions.json')
GR_ = J(MDR_ / 'qm/gr/gr_gate.json')
N2G_ = J(MDR_ / 'qm/n2/n2_gate_result.json')
BEH_ = J(MDR_ / 'qm/n2/final_test/behaviour.json')
RQM_ = J(MDR_ / 'qm/gr/relax_qm.json')
REP_ = J(MDR_ / 'results/campaign/second_program/recheck_report.json')
CAMP_ = J(MDR_ / 'results/campaign/first_program/campaign.json')
STEP1_ = J(MDR_ / 'qm/step1_report.json')
QMC_ = J(MDR_ / 'qm/config_prod.json')
PKA_ = J(MDR_ / 'system/audit/propka.json')
LIGR_ = J(MDR_ / 'ligand/repair.json')
GZ_ = {t_: J(MDR_ / 'production' / f'gate_GZ_bonded_{t_}.json') for t_ in ('2e', '2s')}
_XT = D17M['xtal_6RL9']             # the same quantities in the coordinates of 6RL9 that the systems are built from
_c17 = D17M['comparison']
_ar = D17D['against_the_records']
assert not _c17['mismatches'], _c17['mismatches']
assert not _ar['runs_with_a_difference'], _ar['runs_with_a_difference']

section('40. Molecular dynamics of the carbonic anhydrase II–sulfanilamide complex and the two docking tests (Glu106, four crystal waters)',
        'Molecular dynamics (md/RESULTS.md, record of 8 Oct 2026): hCA II of PDB 6RL9 with the sulfonamidate of sulfanilamide '
        'bonded to the zinc, nine runs of 100 ns with zinc set 2e (three from the crystal pose, start X; three from the docked '
        'near-native pose, start N; three from the docked ring-rotated pose, start R), two sensitivity runs of start R '
        '(R1_2s: the other zinc set; R1_n2fit: candidate terms at the sulfonamide nitrogen N2 that gate G-R2 did not accept) and '
        'two validation runs X0 of 5 ns. A frame is crystal-like if d1 < 6.7 Å and d2 > 5.2 Å, rotated-like if d1 > 6.7 Å and '
        'd2 < 5.2 Å, intermediate otherwise (d1: aniline N1 to Phe131 CZ; d2: N1 to Asn67 OD1); the verdict of a run reads its '
        'last 20 ns. Times "from the release" count from the start of equilibration stage 4, when the ligand is released; '
        'the production starts 1,300 ps later. '
        f"derive/d17_md.py forms the series-based numbers again from md/results/series/ (state of every frame from d1 and d2; "
        f"{_c17['this_script_against_the_first_program']['frames_whose_state_formed_here_differs_from_the_file']} frames classed "
        f"differently from the analysis; means within {_c17['this_script_against_the_first_program']['largest_difference_of_a_mean']:.6f} Å), "
        f"and the all-atom counts, gate G-Z and part A of gate G-R from the tracked files. It also compares {_c17['second_program']['n_values']} "
        f"summary values with the second analysis program (largest difference {_c17['second_program']['largest_difference']:.6f}) and "
        f"{_c17['third_reading']['n_values']} values with the blind third reading (counts, verdicts and times equal; largest difference of a "
        f"mean {_c17['third_reading']['largest_difference_of_a_mean']:.5f} Å). "
        'Read from their source only, because git holds no per-frame data for them: the N–Zn–N angles, C4–S–N2, τ_fold, the '
        'temperature and density of the state logs, His64 and the water of the site, gates G-L, G-V (solute frames) and G-R2, and '
        'the quantum-chemical relaxation. The Zn–N distances of set 2e are put in by the parameter set (their equilibrium values '
        'are those of 6RL9) and are no result. Fractions of time are not populations, and no rate, barrier or free energy '
        'follows from these runs (md/PLAN.md, rule 5). '
        'Docking tests: every AutoDock4Zn, Vina and Vinardo run of dock/GLU106.md (Glu106 as the neutral acid, the state of the '
        'receptors of the paper, against Glu106 as a carboxylate) and of dock/WATER4.md (crystal waters 423, 516, 553 and 568 '
        'of 6RL9 kept as fixed receptor atoms). The runs with the waters are a diagnostic outside the protocol of the paper. '
        f"Our parse of every out.pdbqt reproduces the {_ar['n_runs_compared']} run records of the two tests: mean scores to "
        f"{_ar['max_abs_score_mean_diff']:.4f} kcal/mol, mean RMSDs to {_ar['max_abs_rmsd_mean_diff']:.4f} Å, and every count of seeds "
        f'and modes ({D17S}). Near-native: a pose within 1.0 Å of the crystal ligand; rotated pose: beyond 2.0 Å from the '
        'crystal ligand and within 2.0 Å of the reference rotated pose (seed 14, mode 3 of the production run of the paper).')

# =========================================================================================== system and protocol
_bx = BUILD_[('bonded_2e', 'X')]
_SB = 'md/topology/bonded_2e/{}/build.json'
add('md.sys.pdb', '6RL9', None, None, 'md/system/decisions.json:sources; md/SYSTEM.md section 1', 'crystal structure of the simulated complex (hCA II with sulfanilamide, 1.0 Å)')
for s_ in 'XNR':
    b_ = BUILD_[('bonded_2e', s_)]
    c_ = b_['composition']
    add(f'md.sys.{s_}.atoms', c_['atoms'], None, 0, _SB.format(s_) + ':composition.atoms', f'start {s_}: atoms of the simulated system')
    add(f'md.sys.{s_}.waters', c_['water_molecules'], None, 0, _SB.format(s_) + ':composition.water_molecules', f'start {s_}: water molecules')
    add(f'md.sys.{s_}.xtal_waters', c_['crystal_water_molecules'], None, 0, _SB.format(s_) + ':composition.crystal_water_molecules',
        f'start {s_}: crystal water molecules of 6RL9 kept')
    g_ = b_['geometry_as_built']['values']
    for k_, kk_, u_, p_, m_ in (('d1', 'd1_N1_Phe131CZ_A', 'A', 2, 'd1 (N1 to Phe131 CZ)'), ('d2', 'd2_N1_Asn67OD1_A', 'A', 2, 'd2 (N1 to Asn67 OD1)'),
                                ('zn_n2', 'Zn_N2_A', 'A', 3, 'Zn–N2'), ('rmsd', 'lig_rmsd_symmetry_corrected_A', 'A', 2, 'ligand heavy-atom RMSD to the crystal pose'),
                                ('tau_fold', 'tau_fold_deg', 'deg', 1, 'τ_fold'), ('c4sn2', 'C4_S_N2_deg', 'deg', 1, 'C4–S–N2')):
        add(f'md.sys.{s_}.{k_}', g_[kk_], u_, p_, _SB.format(s_) + f':geometry_as_built.values.{kk_}', f'start {s_} as built: {m_}')
for s_ in 'NR':
    q_ = DEC_['starts'][s_]
    add(f'md.sys.{s_}.dock_score', q_['score_kcal_mol'], 'kcal mol-1', 3, f'md/system/decisions.json:starts.{s_}.score_kcal_mol',
        f'start {s_}: AutoDock4Zn score of the docked pose it was built from ({q_["source"].split(":")[1].strip().split(", neutral")[0]})')
    add(f'md.sys.{s_}.dock_rmsd', q_['rmsd_to_crystal_A'], 'A', 2, f'md/system/decisions.json:starts.{s_}.rmsd_to_crystal_A',
        f'start {s_}: RMSD of that docked pose to the crystal ligand')
add('md.sys.xtal_waters_in_entry', D17D['waters']['n_waters_in_the_entry'], None, 0, f'{D17S}:dock.waters.n_waters_in_the_entry (dock/pdb/6RL9.cif.gz)',
    'water sites of chain A in the deposited entry 6RL9 (215 are kept; site 510, occupancy 0.32, is dropped; start R loses four more)')
add('md.sys.R.waters_removed', len(DEC_['starts']['R']['waters_removed']), None, 0, 'md/system/decisions.json:starts.R.waters_removed',
    'crystal waters removed to build start R because they lie within 2.2 Å of a heavy atom of the rotated ligand (423, 516, 553, 568)')
_w510 = [x_ for x_ in DEC_['water']['labelled_or_partial'] if not x_['kept'] and not any(y_['kept'] for y_ in DEC_['water']['labelled_or_partial'] if y_['water'] == x_['water'])]
assert DEC_['water']['sites_in_entry'] - DEC_['water']['kept'] == 1 and len({x_['water'] for x_ in _w510}) == 1
add('md.sys.water_dropped', _w510[0]['water'], None, 0, 'md/system/decisions.json:water.labelled_or_partial (the one site of which no copy is kept)', 'number of the water site of 6RL9 that is not kept (its only copy clashes with a conformer taken)')
add('md.sys.water_dropped_occ', _w510[0]['occupancy'], None, 2, 'md/system/decisions.json:water.labelled_or_partial[water 510].occupancy', 'occupancy of that water site')
add('md.sys.charge_protein', DEC_['protein']['formal_charge']['net'], 'e', 0, 'md/system/decisions.json:protein.formal_charge.net', 'formal charge of the protein with Glu106 as a carboxylate and the acetyl cap (the docking receptor is neutral)')
add('md.sys.water_removal_cut', 2.2, 'A', 1, 'md/SYSTEM.md section 4; md/build_system.py', 'a crystal water within this distance of a ligand heavy atom is removed when a start is built')
add('md.sys.na', _bx['composition']['Na+'], None, 0, _SB.format('X') + ':composition.Na+', 'sodium ions (the same in every system)')
add('md.sys.cl', _bx['composition']['Cl-'], None, 0, _SB.format('X') + ':composition.Cl-', 'chloride ions (the same in every system)')
assert all(BUILD_[('bonded_2e', s_)]['composition'][i_] == _bx['composition'][i_] for s_ in 'NR' for i_ in ('Na+', 'Cl-', 'solute_atoms'))
add('md.sys.solute_atoms', _bx['composition']['solute_atoms'], None, 0, _SB.format('X') + ':composition.solute_atoms', 'atoms of the solute (protein with its acetyl cap, zinc, ligand)',
    check=D17M['runs']['X1']['final_frame']['n_atoms'], tol=0.01)
add('md.sys.protein_atoms', _bx['composition']['protein_atoms'], None, 0, _SB.format('X') + ':composition.protein_atoms', 'atoms of the protein with its acetyl cap',
    check=D17M['runs']['X1']['final_frame']['n_protein_atoms'], tol=0.01)
add('md.sys.ligand_atoms', _bx['composition']['ligand_atoms'], None, 0, _SB.format('X') + ':composition.ligand_atoms', 'atoms of the ligand (the sulfonamidate anion ArSO₂NH⁻)',
    check=D17M['runs']['X1']['final_frame']['n_ligand_atoms'], tol=0.01)
add('md.sys.ligand_heavy_atoms', len(S1_['X1']['constants']['ring_wide_atoms']) + 4, None, 0, 'md/analysis/SPEC.md section 3 (11 heavy atoms: C1–C6, N1, N2, O1, O2, S); selections.json lig_heavy',
    'heavy atoms of the ligand, over which its RMSD is taken', check=D17M['runs']['X1']['final_frame']['n_heavy_ligand'], tol=0.01)
add('md.sys.residues', S1_['X1']['run']['n_CA'], None, 0, 'md/results/runs/bonded_2e/X1/first_program_summary.json:run.n_CA', 'protein residues (His3 to Lys261 of 6RL9, which has no residue 126); their Cα atoms are fitted',
    check=D17M['runs']['X1']['final_frame']['n_residues_with_CA'], tol=0.01)
with open(MDR_ / 'system/residue_map.csv', newline='') as _fh:
    _rmap = [x_ for x_ in __import__('csv').DictReader(_fh) if '(added)' not in x_['deposited_name']]
_rnum = [int(x_['deposited_number']) for x_ in _rmap]
assert len(_rnum) == S1_['X1']['run']['n_CA']
add('md.sys.res_first', _rnum[0], None, 0, 'md/system/residue_map.csv:deposited_number (first protein residue)', f'first residue of the simulated protein in the numbering of 6RL9 ({_rmap[0]["deposited_name"].capitalize()}{_rnum[0]}; an acetyl cap precedes it)')
add('md.sys.res_last', _rnum[-1], None, 0, 'md/system/residue_map.csv:deposited_number (last protein residue)', f'last residue of the simulated protein ({_rmap[-1]["deposited_name"].capitalize()}{_rnum[-1]})')
add('md.sys.res_missing', [n_ for n_ in range(_rnum[0], _rnum[-1] + 1) if n_ not in _rnum][0], None, 0, 'md/system/residue_map.csv:deposited_number (the one number not used)', 'residue number that 6RL9 does not have')
assert len([n_ for n_ in range(_rnum[0], _rnum[-1] + 1) if n_ not in _rnum]) == 1
add('md.sys.charge_solute', _bx['charge']['solute_with_crystal_water_e'], 'e', 0, _SB.format('X') + ':charge.solute_with_crystal_water_e', 'net charge of the complex (protein −2, zinc +2, ligand −1)')
add('md.sys.box', 'truncated octahedron', None, None, _SB.format('X') + ':box.shape', 'shape of the periodic cell')
add('md.sys.box_buffer', _bx['box']['solvateOct_buffer_A'], 'A', 0, _SB.format('X') + ':box.solvateOct_buffer_A', 'smallest distance from the solute to the cell wall when the cell was built')
add('md.sys.box_edge', _bx['box']['edge_A'], 'A', 1, _SB.format('X') + ':box.edge_A', 'cell edge as built, start X (truncated octahedron, angles 109.47°)')
add('md.sys.salt_target', _bx['ions']['target_mol_per_l'], None, 2, _SB.format('X') + ':ions.target_mol_per_l', 'NaCl concentration aimed at (mol L⁻¹), counted in the volume of the solvent')
add('md.sys.salt_pairs', _bx['ions']['pairs'], None, 0, _SB.format('X') + ':ions.pairs', 'NaCl pairs added (one more Na⁺ neutralises the complex)')
add('md.sys.salt_in_box', _bx['ions']['concentration_in_box_volume_as_built_mol_per_l'], None, 2, _SB.format('X') + ':ions.concentration_in_box_volume_as_built_mol_per_l',
    'the same salt counted in the volume of the cell as built (mol L⁻¹)')
add('md.sys.zn_small_atoms', STEP1_['small']['n_atoms'], None, 0, 'md/qm/step1_report.json:small.n_atoms', 'atoms of the small zinc-site model (zinc, three ethylimidazoles in place of His94, His96 and His119, the ligand); charge +1')
add('md.sys.zn_large_atoms', STEP1_['large']['n_atoms'], None, 0, 'md/qm/step1_report.json:large.n_atoms', 'atoms of the large capped zinc-site model used for the electrostatic potential; charge +1')

_sy, _it, _ba, _eq = PROT_['system'], PROT_['integrator'], PROT_['barostat'], PROT_['equilibration']
_PS = 'md/protocol.json'
add('md.prot.cutoff', _sy['nonbonded_cutoff_nm'], 'nm', 1, f'{_PS}:system.nonbonded_cutoff_nm', 'real-space cut-off of the non-bonded interactions (particle-mesh Ewald)')
add('md.prot.ewald_tol', _sy['ewald_error_tolerance'], None, 4, f'{_PS}:system.ewald_error_tolerance', 'Ewald error tolerance')
add('md.prot.h_mass', _sy['hydrogen_mass_amu'], 'amu', 3, f'{_PS}:system.hydrogen_mass_amu', 'mass of the hydrogens bonded to a heavy atom (hydrogen mass repartitioning; the hydrogens of rigid water keep 1.008 amu)')
add('md.prot.h_mass_water', float(re.search(r'rigid water keep ([\d.]+) amu', _sy['hydrogen_mass_repartitioning']).group(1)), 'amu', 3, f'{_PS}:system.hydrogen_mass_repartitioning (text)',
    'mass of the hydrogens of rigid water, which are not repartitioned')
add('md.prot.temperature', _it['temperature_K'], 'K', 0, f'{_PS}:integrator.temperature_K', 'temperature of the Langevin thermostat (LangevinMiddle integrator)')
add('md.prot.friction', _it['friction_per_ps'], None, 0, f'{_PS}:integrator.friction_per_ps', 'friction coefficient of the thermostat (ps⁻¹)')
add('md.prot.pressure', _ba['pressure_bar'], 'bar', 0, f'{_PS}:barostat.pressure_bar', 'pressure of the Monte Carlo barostat')
add('md.prot.barostat_steps', _ba['frequency_steps'], None, 0, f'{_PS}:barostat.frequency_steps', 'steps between two attempts of the barostat')
add('md.prot.min_k', PROT_['minimisation']['restraint_k_kj_per_mol_nm2'], 'kJ mol-1', 0, f'{_PS}:minimisation.restraint_k_kj_per_mol_nm2',
    'force constant of the position restraints of the minimisation, kJ mol⁻¹ nm⁻² (heavy atoms of the solute and crystal-water oxygens); the unit field leaves out nm⁻²')
for i_, st_ in enumerate(_eq['stages'], 1):
    add(f'md.prot.eq{i_}.length', st_['length_ps'], 'ps', 0, f'{_PS}:equilibration.stages[{i_ - 1}].length_ps', f'equilibration stage {i_} ({st_["name"]}, {st_["ensemble"]}): length')
    add(f'md.prot.eq{i_}.timestep', st_['timestep_fs'], 'fs', 0, f'{_PS}:equilibration.stages[{i_ - 1}].timestep_fs', f'equilibration stage {i_}: time step')
    add(f'md.prot.eq{i_}.k', st_['restraint_k_kj_per_mol_nm2'], 'kJ mol-1', 0, f'{_PS}:equilibration.stages[{i_ - 1}].restraint_k_kj_per_mol_nm2',
        f'equilibration stage {i_}: force constant of the position restraints, kJ mol⁻¹ nm⁻² (atoms held: {st_["restraint_atoms"].replace("_", " ")}); the unit field leaves out nm⁻²')
    add(f'md.prot.eq{i_}.held', st_['restraint_atoms'].replace('_', ' '), None, None, f'{_PS}:equilibration.stages[{i_ - 1}].restraint_atoms', f'equilibration stage {i_}: atoms held')
add('md.prot.eq1.T_start', _eq['stages'][0]['temperature_K'][0], 'K', 0, f'{_PS}:equilibration.stages[0].temperature_K', 'equilibration stage 1: temperature at its start (heated to 300 K)')
add('md.prot.eq_total', _eq['total_ps'], 'ps', 0, f'{_PS}:equilibration.total_ps', 'length of the equilibration without an extension step')
add('md.prot.eq_released', sum(st_['length_ps'] for st_ in _eq['stages'][3:]), 'ps', 0, f'{_PS}:equilibration.stages[3:6].length_ps (sum)',
    'length of the equilibration after the release of the ligand (stages 4, 5 and 6)')
add('md.prot.eq_released_ns', sum(st_['length_ps'] for st_ in _eq['stages'][3:]) / 1000.0, 'ns', 1, f'{_PS}:equilibration.stages[3:6].length_ps (sum, in ns)',
    'the same length in ns: the released part of the equilibration that the two figures draw at negative time')
add('md.prot.eq_stage4_end', _eq['stages'][3]['length_ps'], 'ps', 0, f'{_PS}:equilibration.stages[3].length_ps', 'end of stage 4, counted from the release')
add('md.prot.eq_stage5_end', _eq['stages'][3]['length_ps'] + _eq['stages'][4]['length_ps'], 'ps', 0, f'{_PS}:equilibration.stages[3:5].length_ps (sum)', 'end of stage 5, counted from the release')
_hy = _eq['hydration_check']
add('md.prot.hyd.cut', 10 * _hy['water_cutoff_nm'], 'A', 1, f'{_PS}:equilibration.hydration_check.water_cutoff_nm', 'hydration check at the end of stage 3: a water oxygen within this distance of a crystal position of a ring carbon counts')
add('md.prot.hyd.applies', 10 * _hy['applies_if_ligand_mean_distance_above_nm'], 'A', 1, f'{_PS}:equilibration.hydration_check.applies_if_ligand_mean_distance_above_nm',
    'the hydration check applies only if the ring carbons of the ligand are on average farther than this from their crystal positions')
add('md.prot.hyd.ext_length', _hy['extension_length_ps'], 'ps', 0, f'{_PS}:equilibration.hydration_check.extension_length_ps', 'length of an extension step of stage 3 when the count is zero')
add('md.prot.hyd.max_ext', _hy['max_extensions'], None, 0, f'{_PS}:equilibration.hydration_check.max_extensions', 'largest number of extension steps')
_pd, _ou = PROT_['production'], PROT_['output']
add('md.prot.prod.length', _pd['length_ns'], 'ns', 0, f'{_PS}:production.length_ns', 'length of a production run')
add('md.prot.prod.timestep', _pd['timestep_fs'], 'fs', 0, f'{_PS}:production.timestep_fs', 'time step of the production (bonds to hydrogen constrained, heavy hydrogens)')
add('md.prot.out.solute', _ou['solute_frame_interval_ps'], 'ps', 0, f'{_PS}:output.solute_frame_interval_ps', 'interval of the solute frames')
add('md.prot.out.full', _ou['full_frame_interval_ps'], 'ps', 0, f'{_PS}:output.full_frame_interval_ps', 'interval of the all-atom frames')
_hl = [HR_[l_] for l_ in HR_]
_long = [HR_[LAB_[i_]] for i_ in ELEVEN_]
add('md.prot.n_runs_primary', CAMP_['rules_per_model']['bonded_2e']['n_runs'], None, 0, 'md/results/campaign/first_program/campaign.json:rules_per_model.bonded_2e.n_runs', 'runs of 100 ns of the primary set (three per start)')
add('md.prot.n_runs_per_start', CAMP_['rules_per_model']['bonded_2e']['runs_per_start']['X'], None, 0, 'md/results/campaign/first_program/campaign.json:rules_per_model.bonded_2e.runs_per_start', 'runs per start')
add('md.prot.n_runs_sensitivity', len(ELEVEN_) - len(NINE_), None, 0, 'md/results/run_health.json:runs (bonded_2s/R1, bonded_2e_n2fit/R1)', 'sensitivity runs of 100 ns, both from start R')
add('md.prot.n_runs_long', len(_long), None, 0, 'md/results/run_health.json:runs', 'runs of 100 ns in all')
add('md.prot.n_runs_validation', sum(1 for l_ in HR_ if l_.endswith('/X0')), None, 0, 'md/results/run_health.json:runs (bonded_2e/X0, bonded_2s/X0)', 'validation runs X0 from the crystal pose, one for each zinc set')
add('md.prot.total_production_ns', sum(h_['driver']['production']['completed_ns'] for h_ in _hl), 'ns', 0, 'md/results/run_health.json:runs[].driver.production.completed_ns (sum over 13 runs)',
    'simulated production time of the campaign (900 ns primary, 200 ns sensitivity, 10 ns validation)')
add('md.prot.primary_production_ns', sum(HR_[LAB_[i_]]['driver']['production']['completed_ns'] for i_ in NINE_), 'ns', 0, 'md/results/run_health.json:runs[].driver.production.completed_ns (sum over the nine)', 'simulated production time of the nine primary runs')
add('md.prot.total_equilibration_ns', sum(h_['driver']['equilibration']['completed_ps'] for h_ in _hl) / 1000.0, 'ns', 1, 'md/results/run_health.json:runs[].driver.equilibration.completed_ps (sum over 13 runs)',
    'simulated equilibration time of the campaign')
add('md.prot.validation_ns', HR_['bonded_2e/X0']['driver']['production']['completed_ns'], 'ns', 0, 'md/results/run_health.json:runs[bonded_2e/X0].driver.production.completed_ns', 'length of a validation run X0')
_rng('md.prot.ns_per_day', ELEVEN_, lambda i_: HR_[LAB_[i_]]['driver']['production']['ns_per_day'], None, 0, 'md/results/run_health.json:runs[].driver.production.ns_per_day',
     'speed of a production run over the eleven runs (ns per day, one NVIDIA L40S); the driver writes two decimals, so the fact prints whole numbers')
add('md.prot.frames_solute', S1_['X1']['run']['n_frames_solute'], None, 0, 'md/results/runs/bonded_2e/X1/first_program_summary.json:run.n_frames_solute', 'solute frames of a production run',
    check=D17M['runs']['X1']['n_production'], tol=0.01)
add('md.prot.frames_full', S1_['X1']['run']['n_frames_full'], None, 0, 'md/results/runs/bonded_2e/X1/first_program_summary.json:run.n_frames_full', 'all-atom frames of a production run',
    check=D17M['runs']['X1']['all_atom']['n'], tol=0.01)
add('md.prot.frames_released', EV_['bonded_2e/X1']['released_equilibration']['n'], None, 0, 'md/results/run_events.json:runs[].released_equilibration.n', 'solute frames of the equilibration after the release of the ligand',
    check=D17M['runs']['X1']['n_equilibration'], tol=0.01)
assert all(S1_[i_]['run']['n_frames_solute'] == 10000 and S1_[i_]['run']['n_frames_full'] == 1000 and S1_[i_]['run']['complete'] for i_ in ELEVEN_)

# definitions of the analysis (constants of the specification as the first program writes them)
_k = S1_['X1']['constants']
_KS = 'md/results/runs/bonded_2e/X1/first_program_summary.json:constants'
add('md.def.mid_d1', _k['midpoint_d1_N1_F131CZ_A'], 'A', 1, f'{_KS}.midpoint_d1_N1_F131CZ_A', 'midpoint of d1 (N1 to Phe131 CZ) that separates the two states')
add('md.def.mid_d2', _k['midpoint_d2_N1_N67OD1_A'], 'A', 1, f'{_KS}.midpoint_d2_N1_N67OD1_A', 'midpoint of d2 (N1 to Asn67 OD1)')
add('md.def.window', _k['verdict_window_ns'], 'ns', 0, f'{_KS}.verdict_window_ns', 'the verdict of a run reads its last 20 ns')
add('md.def.verdict_pct', 100.0 * _k['verdict_fraction'][0] / _k['verdict_fraction'][1], '%', 0, f'{_KS}.verdict_fraction', 'share of the frames of the last 20 ns in one state for CONVERTED or STAYED; MIXED otherwise')
add('md.def.dwell', _k['dwell_ns'], 'ns', 0, f'{_KS}.dwell_ns', 'a visit of at least this length in a state is a long visit (the time of a conversion is the first frame of the first long visit in the opposite state)')
add('md.def.hb_dist', _k['hbond_distance_max_A'], 'A', 1, f'{_KS}.hbond_distance_max_A', 'hydrogen bond N2–H···OG1(Thr199): largest N2···OG1 distance')
add('md.def.hb_angle', _k['hbond_angle_min_deg'], 'deg', 0, f'{_KS}.hbond_angle_min_deg', 'hydrogen bond N2–H···OG1(Thr199): smallest angle at H')
add('md.def.fifth', _k['fifth_ligand_distance_max_A'], 'A', 1, f'{_KS}.fifth_ligand_distance_max_A', 'an oxygen within this distance of the zinc is counted as a fifth ligand')
add('md.def.rule1', _k['rule1_rmsd_A'], 'A', 1, f'{_KS}.rule1_rmsd_A', 'rule 1: a crystal-start run has left the crystal pose if its mean ligand RMSD of the last half exceeds this')
# definitions of the checks on the raw frames (md/production/site_checks.py), added on 11 Oct 2026 for Table S38; each
# value is read from the result files that the program wrote, so that a change of the program would stop the build
_sfl = [x_ for x_ in SF_.values() if 'crystal_water_places' in x_]         # the eleven runs of 100 ns (the two validation runs have no such block)
_ssl = [x_ for x_ in SS_.values() if 'his64' in x_]
assert len(_sfl) == 11 and len(_ssl) == 11
_ws = {re.search(r'within ([0-9.]+) A of the place of the crystal water', x_['crystal_water_places']['what']).group(1) for x_ in _sfl}
_h64 = {re.search(r'chi1 between 0 and ([0-9.]+) deg with the distance at most the limit', x_['his64']['what']).group(1) for x_ in _ssl}
assert _ws == {'1.4'} and _h64 == {'120'}, (_ws, _h64)
assert all(set(x_['his64']['deep_in_share_20_to_100_ns']) == {'6.5', '7.0', '7.5'} for x_ in _ssl)
assert all('n_abs_V_below_0.2' in x_['nitrogen']['production'] for x_ in _ssl)
add('md.def.water_site', float(min(_ws)), 'A', 1, 'md/results/site_checks_full.json:runs[].crystal_water_places.what (constant OCC of md/production/site_checks.py)',
    'the place of a crystal water counts as occupied in an all-atom frame when a water oxygen lies within this distance of the crystal position after the fit of the Cα atoms')
add('md.def.his64_deep', 7.0, 'A', 1, 'md/results/site_checks_solute.json:runs[].his64.deep_in_share_20_to_100_ns["7.0"] (constant DEEP of md/production/site_checks.py)',
    'His64 counts as deep in the site when Zn–NE2(His64) is at most this and χ1 lies between 0° and 120°; the limit was chosen in the physical review of the campaign after looking at the data '
    '(6.5 and 7.5 Å are stored beside it)')
add('md.def.his64_chi1_max', float(min(_h64)), 'deg', 0, 'md/results/site_checks_solute.json:runs[].his64.what (md/production/site_checks.py)',
    'upper end of the range of χ1 (N–CA–CB–CG, from 0°) of His64 in the definition of deep in the site')
add('md.def.n2_flat', 0.2, None, 1, 'md/results/site_checks_solute.json:runs[].nitrogen.production.n_abs_V_below_0.2 and md/results/run_events.json:runs[].n2_volume.production.n_abs_below_0.2 (the limit is in the key names)',
    'the nitrogen N2 counts as nearly planar in a frame when the signed volume of the unit vectors from N2 to Zn, HN21 and S is below this in size')

# ----------------------------------------------------------------------------------------------- force field
_ff = _bx['force_field']
add('md.ff.protein', _ff['protein'], None, None, _SB.format('X') + ':force_field.protein', 'force field of the protein')
add('md.ff.water', 'TIP3P', None, None, _SB.format('X') + ':force_field.water', 'water model (rigid)')
add('md.ff.ligand', 'GAFF2 ' + re.search(r'Version ([\d.]+)', _ff['gaff2_version_line']).group(1), None, None, _SB.format('X') + ':force_field.gaff2_version_line',
    'general force field of the ligand, with the terms at sulfur repaired against the DFT of the paper and RESP charges of the zinc site')
add('md.ff.ions', 'Joung–Cheatham', None, None, _SB.format('X') + ':force_field.sodium_chloride', 'parameters of Na⁺ and Cl⁻ (for TIP3P)')
for t_ in _bx['ligand_parameters']['repaired_terms']:
    k_ = re.sub(r'[^A-Za-z0-9]+', '_', t_['term']).strip('_')
    if isinstance(t_.get('in_prmtop'), list) and len(t_['in_prmtop']) == 2 and all(isinstance(x_, (int, float)) for x_ in t_['in_prmtop']):
        add(f'md.ff.lig.{k_}.k', t_['in_prmtop'][0], 'kcal mol-1', 2, _SB.format('X') + f':ligand_parameters.repaired_terms[{t_["term"]}].in_prmtop[0]',
            f'repaired ligand term, {t_["term"]}: force constant in the topology (Amber form; kcal mol⁻¹ Å⁻² for a bond, kcal mol⁻¹ rad⁻² for an angle)')
        add(f'md.ff.lig.{k_}.eq', t_['in_prmtop'][1], 'A' if t_['term'].startswith('bond') else 'deg', 4 if t_['term'].startswith('bond') else 3,
            _SB.format('X') + f':ligand_parameters.repaired_terms[{t_["term"]}].in_prmtop[1]', f'repaired ligand term, {t_["term"]}: equilibrium value in the topology')

# the bonded zinc site: both parameter sets
_NM = {'SN1304:N2': 'N2', 'HD194:NE2': 'His94', 'HD296:NE2': 'His96', 'HE1119:ND1': 'His119'}
for t_ in ('2e', '2s'):
    zb_ = BUILD_[(f'bonded_{t_}', 'X')]['zinc_site_bonded']
    for b_ in zb_['bonds_to_zinc_against_the_parameter_files']:
        n_ = _NM[b_['bond'].split('Zn-')[1]]
        add(f'md.zn.{t_}.bond.{n_}.k', b_['k_kcal_per_mol_A2'], 'kcal mol-1', 1, f'md/topology/bonded_{t_}/X/build.json:zinc_site_bonded.bonds_to_zinc_against_the_parameter_files[{b_["bond"]}].k_kcal_per_mol_A2',
            f'zinc set {t_}: force constant of the bond Zn–{"N2" if n_ == "N2" else "N(" + n_ + ")"}, kcal mol⁻¹ Å⁻² (Amber form, U = k (r − r0)²); the unit field leaves out Å⁻²')
        add(f'md.zn.{t_}.bond.{n_}.r0', b_['r0_A'], 'A', 3 if t_ == '2e' else 4, f'md/topology/bonded_{t_}/X/build.json:zinc_site_bonded.bonds_to_zinc_against_the_parameter_files[{b_["bond"]}].r0_A',
            f'zinc set {t_}: equilibrium length of the bond Zn–{"N2" if n_ == "N2" else "N(" + n_ + ")"}' + (' (the value of 6RL9, put in)' if t_ == '2e' else ' (optimised small model)'))
        if t_ == '2e':
            add(f'md.zn.xtal.{n_}', b_['distance_in_6RL9_A'], 'A', 3, f'md/topology/bonded_2e/X/build.json:zinc_site_bonded.bonds_to_zinc_against_the_parameter_files[{b_["bond"]}].distance_in_6RL9_A',
                f'6RL9: distance Zn–{"N2" if n_ == "N2" else "N(" + n_ + ")"}')
    ka_ = sorted({a_['k_kcal_per_mol_rad2'] for a_ in zb_['angles_with_zinc_against_the_parameter_files'] if a_['kind'] == 'N-Zn-N'})
    kc_ = sorted({a_['k_kcal_per_mol_rad2'] for a_ in zb_['angles_with_zinc_against_the_parameter_files'] if a_['kind'] == 'Zn-N-C'})
    if t_ == '2e':
        assert len(ka_) == 1 and len(kc_) == 1
        add('md.zn.2e.k_N_Zn_N', ka_[0], 'kcal mol-1', 0, 'md/topology/bonded_2e/X/build.json:zinc_site_bonded.angles_with_zinc_against_the_parameter_files (kind N-Zn-N)',
            'zinc set 2e: force constant of the six N–Zn–N angles, kcal mol⁻¹ rad⁻²; their equilibrium values are the angles of 6RL9')
        add('md.zn.2e.k_Zn_N_C', kc_[0], 'kcal mol-1', 0, 'md/topology/bonded_2e/X/build.json:zinc_site_bonded.angles_with_zinc_against_the_parameter_files (kind Zn-N-C)',
            'zinc set 2e: force constant of the six Zn–N–C angles, kcal mol⁻¹ rad⁻²')
    add(f'md.zn.{t_}.n_bonds', len(zb_['bonds_to_zinc_against_the_parameter_files']), None, 0, f'md/topology/bonded_{t_}/X/build.json:zinc_site_bonded.bonds_to_zinc_against_the_parameter_files', f'zinc set {t_}: bonds to the zinc')
    add(f'md.zn.{t_}.n_angles', len(zb_['angles_with_zinc_against_the_parameter_files']), None, 0, f'md/topology/bonded_{t_}/X/build.json:zinc_site_bonded.angles_with_zinc_against_the_parameter_files',
        f'zinc set {t_}: angles with the zinc (six N–Zn–N, six Zn–N–C, Zn–N2–S, Zn–N2–H)')
_zb = _bx['zinc_site_bonded']
add('md.zn.n_set_differences_angles', 12, None, 0, 'md/PRODUCTION.md section 2.1; build.json zinc_site_bonded.n_angles_by_kind (N-Zn-N 6, Zn-N-C 6)', 'angles in which the two zinc sets differ (they also differ in the four Zn–N bonds and in nothing else)')
add('md.zn.n_torsions', _zb['torsions_through_zinc']['n'], None, 0, _SB.format('X') + ':zinc_site_bonded.torsions_through_zinc.n', 'torsions through the zinc, every one with zero force constant')
add('md.zn.charge_zn', _zb['charges']['zinc_e'], 'e', 4, _SB.format('X') + ':zinc_site_bonded.charges.zinc_e', 'RESP charge of the zinc')
add('md.zn.charge_site', _zb['charges']['site_sum_e'], 'e', 0, _SB.format('X') + ':zinc_site_bonded.charges.site_sum_e', 'sum of the RESP charges of the zinc site (zinc, three histidines, ligand)')
add('md.zn.charge_ligand', _zb['charge_of_zinc_site_residues_e']['SN1'], 'e', 4, _SB.format('X') + ':zinc_site_bonded.charge_of_zinc_site_residues_e.SN1', 'sum of the RESP charges of the ligand in the zinc site')
_n2 = _bx['terms_at_N2']
_znh = _n2['in_the_topology']['angles']['Zn-N2-H']
assert abs(_znh['theta0_deg'] - _n2['Zn_N2_H_in_start_X_deg']) < 0.001
add('md.zn.n2.Zn_N2_H.theta0', _znh['theta0_deg'], 'deg', 3, _SB.format('X') + ':terms_at_N2.in_the_topology.angles["Zn-N2-H"].theta0_deg',
    'terms at N2 as built: equilibrium angle Zn–N2–H in the topology (the angle of start X, 107.4065°, written with three decimals)')
add('md.zn.n2.Zn_N2_H.k', _znh['k_kcal_per_mol_rad2'], 'kcal mol-1', 2, _SB.format('X') + ':terms_at_N2.in_the_topology.angles["Zn-N2-H"].k_kcal_per_mol_rad2', 'terms at N2 as built: force constant of Zn–N2–H, kcal mol⁻¹ rad⁻²')
add('md.zn.charge_N2', _bx['ligand_parameters']['charges_e']['N2'], 'e', 4, _SB.format('X') + ':ligand_parameters.charges_e.N2', 'RESP charge of the sulfonamide nitrogen N2 (the potential fitted is that of pose X alone)')
for a_ in _zb['angles_with_zinc']:
    if a_['atoms'] == 'ZN1:ZN-N2-SN1:S':
        add('md.zn.n2.Zn_N2_S.k', a_['k_kcal_per_mol_rad2'], 'kcal mol-1', 2, _SB.format('X') + ':zinc_site_bonded.angles_with_zinc[ZN1:ZN-N2-SN1:S].k_kcal_per_mol_rad2', 'terms at N2 as built: force constant of Zn–N2–S, kcal mol⁻¹ rad⁻²')
        add('md.zn.n2.Zn_N2_S.theta0', a_['theta0_deg'], 'deg', 2, _SB.format('X') + ':zinc_site_bonded.angles_with_zinc[ZN1:ZN-N2-SN1:S].theta0_deg', 'terms at N2 as built: equilibrium angle Zn–N2–S')
add('md.qm.method', 'B3LYP/6-31G* (pure d and f functions, density fitting)', None, None, 'md/qm/config_prod.json:method (xc b3lyp, basis 6-31g*, cart false, density_fit true)',
    'level of the quantum-chemical calculations of the zinc-site model (PySCF; not the Cartesian 6-31G* of Gaussian)')
assert QMC_['method'] == dict(QMC_['method'], xc='b3lyp', basis='6-31g*', cart=False, density_fit=True)
add('md.qm.gate_method', 'B3LYP-D3(BJ)/6-31G*', None, None, 'md/qm/gr/gr_gate.json:d3; md/GR.md section 2', 'level of the reference energies of gates G-R and G-R2 (the D3(BJ) term added to the B3LYP energies)')

# software (text facts)
_sw = S1_['X1']['software']
for k_, v_, src_, m_ in (
        ('openmm', HR_['bonded_2e/X1']['driver']['openmm'], 'md/results/run_health.json:runs[].driver.openmm', 'OpenMM version of the runs (CUDA platform, mixed precision)'),
        ('gpu', HR_['bonded_2e/X1']['driver']['gpu'], 'md/results/run_health.json:runs[].driver.gpu', 'graphics processor of the runs'),
        ('precision', HR_['bonded_2e/X1']['driver']['precision'], 'md/results/run_health.json:runs[].driver.precision', 'precision of the CUDA platform'),
        ('ambertools', re.search(r'ambertools-([\d.]+)', _bx['software']['ambertools']).group(1), _SB.format('X') + ':software.ambertools', 'AmberTools version (tleap, antechamber, MCPB.py)'),
        ('parmed', _bx['software']['parmed'], _SB.format('X') + ':software.parmed', 'ParmEd version'),
        ('mdanalysis', _sw['MDAnalysis'], 'md/results/runs/bonded_2e/X1/first_program_summary.json:software.MDAnalysis', 'MDAnalysis version of the first analysis program'),
        ('numpy', _sw['numpy'], 'md/results/runs/bonded_2e/X1/first_program_summary.json:software.numpy', 'NumPy version of the analysis (the second program uses NumPy alone)'),
        ('python', _sw['python'], 'md/results/runs/bonded_2e/X1/first_program_summary.json:software.python', 'Python version of the analysis'),
        ('pyscf', J(MDR_ / 'qm/gr/relax/qm_Rx.json')['versions']['pyscf'], 'md/qm/gr/relax/qm_Rx.json:versions.pyscf', 'PySCF version of the quantum-chemical calculations'),
        ('geometric', J(MDR_ / 'qm/gr/relax/qm_Rx.json')['versions']['geometric'], 'md/qm/gr/relax/qm_Rx.json:versions.geometric', 'geomeTRIC version of the relaxations'),
        ('propka', PKA_['propka_version'], 'md/system/audit/propka.json:propka_version', 'PROPKA version')):
    add(f'md.sw.{k_}', v_, None, None, src_, m_)
assert all(HR_[LAB_[i_]]['driver']['openmm'] == HR_['bonded_2e/X1']['driver']['openmm'] and HR_[LAB_[i_]]['driver']['gpu'] == HR_['bonded_2e/X1']['driver']['gpu'] for i_ in ELEVEN_)

# pKa of Glu106 by PROPKA, without and with the zinc in the input
for k_, m_ in (('without_zn', 'without the zinc in the input (the preparation of the docking receptors: Glu106 becomes the neutral acid)'),
               ('with_zn', 'with the zinc in the input (below 7.4: carboxylate, the state of the simulations)')):
    add(f'md.pka.glu106.{k_}', PKA_['pKa']['GLU106'][k_], None, 2, f'md/system/audit/propka.json:pKa.GLU106.{k_}', f'pKa of Glu106 of 6RL9 by PROPKA 3.5.1, {m_}',
        check=D17D['propka'][f'md/system/audit/propka_{k_}.pka'], tol=0.001)
add('md.pka.pH', PKA_['pH'], None, 1, 'md/system/audit/propka.json:pH', 'pH at which the protonation states were set')
for pdb_, tg_ in (('6RL9', 'T3_6RL9'), ('3HS4', 'T4_3HS4')):          # the PROPKA run of the docking preparation itself (no zinc in its input)
    v_ = [x_ for x_ in J(DOCK / tg_ / 'prep_report.json')['protonation_variants'] if x_['residue'] == 'A:GLU106']
    assert len(v_) == 1 and v_[0]['variant'] == 'GLH'
    add(f'dock.glu106.{pdb_}.pka_prep', float(re.search(r'pKa ([\d.]+)', v_[0]['reason']).group(1)), None, 2, f'dock/{tg_}/prep_report.json:protonation_variants[A:GLU106].reason',
        f'pKa of Glu106 in the preparation of the docking receptor {pdb_} (PROPKA 3.5.1 without the zinc in its input): above 7.4, so Glu106 was built as the neutral acid'
        + ('' if pdb_ == '6RL9' else '; with the zinc PROPKA was rerun for 6RL9 only'), check=D17D['propka'][f'dock/{tg_}/propka.pka'], tol=0.001)

# ====================================================================================================== the runs
for m_, run_, rid_ in RUNS17_:
    s_, lab_ = S1_[rid_], f'{m_}/{run_}'
    SRC_ = f'md/results/runs/{m_}/{run_}/first_program_summary.json'
    q_, o_, ev_, ct_, ss_, sf_, h_ = D17M['runs'][rid_], s_['orientation'], EV_[lab_], CT_[lab_], SS_[lab_], SF_[lab_], HR_[lab_]
    L_ = RLAB_.get(rid_, f'run {rid_}')
    b_ = f'md.run.{rid_}'
    assert o_['verdict'] == q_['verdict'], rid_
    add(b_ + '.start', s_['run']['start'], None, None, f'{SRC_}:run.start', f'{L_}: start (X crystal pose, N docked near-native pose, R docked ring-rotated pose)')
    add(b_ + '.zinc_set', h_['build']['parameter_set'], None, None, f'md/results/run_health.json:runs[{lab_}].build.parameter_set', f'{L_}: zinc parameter set (2e: equilibrium values of 6RL9, empirical force constants; 2s: Seminario)')
    add(b_ + '.n2_terms', 'candidate, not accepted' if h_['build']['terms_at_N2']['entries_beyond_the_terms_as_built'] else 'as built', None, None,
        f'md/results/run_health.json:runs[{lab_}].build.terms_at_N2.entries_beyond_the_terms_as_built', f'{L_}: force-field terms at the sulfonamide nitrogen N2 (as built: failed gate G-R; candidate: not accepted by gate G-R2)')
    add(b_ + '.verdict', o_['verdict'], None, None, f'{SRC_}:orientation.verdict', f'{L_}: verdict over the last 20 ns', note='the same verdict from the series of md/results/series/ (derive/d17_md.py) and in the blind third reading')
    for w_, kc_, lab2_ in (('last20', '_window', 'last 20 ns'), ('whole', '', 'whole production')):
        n_ = o_['window_frames'] if w_ == 'last20' else s_['run']['n_frames_solute']
        for st_, nm_ in (('C', 'crystal-like'), ('R', 'rotated-like'), ('I', 'intermediate')):
            c_ = o_[f'count_{st_}{kc_}']
            add(f'{b_}.{w_}.n_{st_}', c_, None, 0, f'{SRC_}:orientation.count_{st_}{kc_}', f'{L_}: {nm_} frames, {lab2_} (of {n_:,})',
                check=q_['last_20ns' if w_ == 'last20' else 'whole'][st_], tol=0.01)
            add(f'{b_}.{w_}.pct_{st_}', 100.0 * c_ / n_, '%', 2, f'{SRC_}:orientation.count_{st_}{kc_} / {n_}', f'{L_}: share of {nm_} frames, {lab2_}',
                check=q_['last_20ns' if w_ == 'last20' else 'whole'][f'pct_{st_}'])
    st_ = s_['stats']['lig_rmsd_A']
    _a(b_ + '.rmsd.mean', st_['mean'], 'A', 2, f'{SRC_}:stats.lig_rmsd_A.mean', f'{L_}: ligand heavy-atom RMSD to the crystal pose, mean over the production', check=q_['rmsd']['mean'])
    _a(b_ + '.rmsd.sd', st_['sd'], 'A', 2, f'{SRC_}:stats.lig_rmsd_A.sd', f'{L_}: ligand RMSD, standard deviation', check=q_['rmsd']['sd'])
    _a(b_ + '.rmsd.last_half', st_['mean_last_half'], 'A', 2, f'{SRC_}:stats.lig_rmsd_A.mean_last_half', f'{L_}: ligand RMSD, mean of the last half (50 to 100 ns)', check=q_['rmsd']['mean_last_half'])
    _a(b_ + '.rmsd.last20', ev_['production_blocks_of_20_ns'][4]['lig_rmsd_mean_A'], 'A', 2, f'md/results/run_events.json:runs[{lab_}].production_blocks_of_20_ns[4].lig_rmsd_mean_A',
       f'{L_}: ligand RMSD, mean of the last 20 ns', check=q_['rmsd']['mean_last_20ns'])
    _a(b_ + '.rmsd.max', st_['max'], 'A', 2, f'{SRC_}:stats.lig_rmsd_A.max', f'{L_}: ligand RMSD, largest single frame', check=q_['rmsd']['max'])
    hb_ = s_['fractions']['hb_N2H_T199OG1']
    add(b_ + '.hb.n', hb_['count'], None, 0, f'{SRC_}:fractions.hb_N2H_T199OG1.count', f'{L_}: production frames with the hydrogen bond N2–H···OG1(Thr199) (of {hb_["n"]:,})', check=q_['hb']['count'], tol=0.01)
    add(b_ + '.hb.pct', 100.0 * hb_['count'] / hb_['n'], '%', 2, f'{SRC_}:fractions.hb_N2H_T199OG1.count / n', f'{L_}: share of production frames with the hydrogen bond N2–H···OG1(Thr199)', check=q_['hb']['pct'])
    _a(b_ + '.d1.last_half', s_['stats']['d1_N1_F131CZ_A']['mean_last_half'], 'A', 2, f'{SRC_}:stats.d1_N1_F131CZ_A.mean_last_half', f'{L_}: d1, mean of the last half', check=q_['d1']['mean_last_half'])
    _a(b_ + '.d2.last_half', s_['stats']['d2_N1_N67OD1_A']['mean_last_half'], 'A', 2, f'{SRC_}:stats.d2_N1_N67OD1_A.mean_last_half', f'{L_}: d2, mean of the last half', check=q_['d2']['mean_last_half'])
    # the conversion by the plan's rule and the events after the release
    cv_ = ct_['conversion']
    if cv_:
        add(b_ + '.rule.ps', int(cv_['time_from_release_ps']), 'ps', 0, f'md/results/conversion_times.json:runs[{lab_}].conversion.time_from_release_ps',
            f'{L_}: time of the plan\'s rule from the release of the ligand (first frame of the first visit of at least 1 ns in the crystal-like state)'
            + ('; with the verdict MIXED this is no conversion' if o_['verdict'] == 'MIXED' else '; it lies in the equilibration, before the production'),
            check=q_['rule']['time_from_release_ps'], tol=0.01)
        add(b_ + '.rule.stage', cv_['stage'], None, 0, f'md/results/conversion_times.json:runs[{lab_}].conversion.stage', f'{L_}: equilibration stage of that frame',
            check=int(q_['rule']['part'][2:]), tol=0.01)
        add(b_ + '.rule.text', ct_['reported_as'] if o_['verdict'] != 'MIXED' else 'no conversion (MIXED); first crystal-like visit of at least 1 ns from 60 ps after release', None, None,
            f'md/results/conversion_times.json:runs[{lab_}].reported_as', f'{L_}: the conversion in the wording the plan fixed')
    else:
        add(b_ + '.rule.text', 'none', None, None, f'md/results/conversion_times.json:runs[{lab_}].reported_as', f'{L_}: conversion by the plan\'s rule ({ct_["reported_as"]})',
            note='no visit of at least 1 ns in the opposite state in the series of md/results/series/ either (derive/d17_md.py)')
        assert q_['rule']['time_from_release_ps'] is None
    if s_['run']['start'] == 'R':
        e1_, e2_, e3_ = ev_['first_released_frame_n2_sign_not_of_start'], ev_['first_released_frame_not_in_state_of_start'], ev_['first_released_frame_in_opposite_state']
        add(b_ + '.sign.ps', int(e1_['time_from_release_ps']), 'ps', 0, f'md/results/run_events.json:runs[{lab_}].first_released_frame_n2_sign_not_of_start.time_from_release_ps',
            f'{L_}: first frame after the release with the sign at N2 positive (the configuration of the crystal pose)', check=q_['first_frame_sign_not_of_start']['ps_from_release'], tol=0.01)
        add(b_ + '.sign.stage', e1_['stage'], None, 0, f'md/results/run_events.json:runs[{lab_}].first_released_frame_n2_sign_not_of_start.stage', f'{L_}: equilibration stage of that frame',
            check=int(q_['first_frame_sign_not_of_start']['part'][2:]), tol=0.01)
        add(b_ + '.left.ps', int(e2_['time_from_release_ps']), 'ps', 0, f'md/results/run_events.json:runs[{lab_}].first_released_frame_not_in_state_of_start.time_from_release_ps',
            f'{L_}: first frame after the release that is not rotated-like', check=q_['first_frame_not_in_state_of_start']['ps_from_release'], tol=0.01)
        add(b_ + '.left.stage', e2_['stage'], None, 0, f'md/results/run_events.json:runs[{lab_}].first_released_frame_not_in_state_of_start.stage', f'{L_}: equilibration stage of that frame',
            check=int(q_['first_frame_not_in_state_of_start']['part'][2:]), tol=0.01)
        add(b_ + '.first_C.ps', int(e3_['time_from_release_ps']), 'ps', 0, f'md/results/run_events.json:runs[{lab_}].first_released_frame_in_opposite_state.time_from_release_ps',
            f'{L_}: first crystal-like frame after the release', check=q_['first_frame_in_opposite_state']['ps_from_release'], tol=0.01)
        add(b_ + '.eq.n_R', ev_['released_equilibration']['R'], None, 0, f'md/results/run_events.json:runs[{lab_}].released_equilibration.R', f'{L_}: rotated-like frames among the 130 frames after the release',
            check=q_['released_equilibration']['R'], tol=0.01)
        hc_ = h_['hydration_check']['counts'][-1]
        add(b_ + '.hyd.n', hc_['n_water_O'], None, 0, f'md/results/run_health.json:runs[{lab_}].hydration_check.counts[-1].n_water_O', f'{L_}: water oxygens in the place of the crystal ring when the ligand was released (driver)')
        add(b_ + '.hyd.nearest', hc_['nearest_water_O_A'], 'A', 3, f'md/results/run_health.json:runs[{lab_}].hydration_check.counts[-1].nearest_water_O_A', f'{L_}: nearest water oxygen to a crystal position of a ring carbon at that moment')
        add(b_ + '.hyd.extensions', h_['hydration_check']['n_extensions'], None, 0, f'md/results/run_health.json:runs[{lab_}].hydration_check.n_extensions', f'{L_}: extension steps of stage 3 (500 ps each)')
        hc0_ = h_['hydration_check']['counts'][0]
        assert hc0_['step'] == 'eq03' and len(h_['hydration_check']['counts']) == h_['hydration_check']['n_extensions'] + 1
        add(b_ + '.hyd.n_first', hc0_['n_water_O'], None, 0, f'md/results/run_health.json:runs[{lab_}].hydration_check.counts[0].n_water_O',
            f'{L_}: water oxygens in the place of the crystal ring at the end of stage 3, before any extension step (driver; zero extends the stage)')
        add(b_ + '.hyd.nearest_first', hc0_['nearest_water_O_A'], 'A', 3, f'md/results/run_health.json:runs[{lab_}].hydration_check.counts[0].nearest_water_O_A',
            f'{L_}: nearest water oxygen to a crystal position of a ring carbon at the end of stage 3, before any extension step')
        # the last all-atom frame (one per 100 ps) before the sign at N2 changes: water in the place of the crystal ring
        wb_ = [x_ for x_ in ev_['water_in_the_place_of_the_ring_equilibration'] if x_['stage'] >= 4 and round(1000 * (x_['time_ns'] + 1.3)) < e1_['time_from_release_ps']]
        if wb_:
            add(b_ + '.ring_water_before.ps', int(round(1000 * (wb_[-1]['time_ns'] + 1.3))), 'ps', 0, f'md/results/run_events.json:runs[{lab_}].water_in_the_place_of_the_ring_equilibration[].time_ns',
                f'{L_}: time from the release of the last all-atom frame before the sign at N2 changes (all-atom frames are 100 ps apart)')
            add(b_ + '.ring_water_before.n', wb_[-1]['n_water_ring_place'], None, 0, f'md/results/run_events.json:runs[{lab_}].water_in_the_place_of_the_ring_equilibration[].n_water_ring_place',
                f'{L_}: water oxygens within 2.5 Å of the crystal positions of the ring carbons in that frame')
            add(b_ + '.ring_water_before.nearest', wb_[-1]['d_water_ring_place_min_A'], 'A', 3, f'md/results/run_events.json:runs[{lab_}].water_in_the_place_of_the_ring_equilibration[].d_water_ring_place_min_A',
                f'{L_}: nearest water oxygen to a crystal position of a ring carbon in that frame')
        ca_ = ss_['change_at_N2']
        for k_, kk_, m2_ in (('tors_S.before', 'torsion_S_before_deg', 'NE2(His94)–Zn–N2–S, mean of the released frames before the change of sign'), ('tors_S.prod', 'torsion_S_production_deg', 'NE2(His94)–Zn–N2–S, mean of the production'),
                             ('tors_H.before', 'torsion_H_before_deg', 'NE2(His94)–Zn–N2–H, mean of the released frames before the change of sign'), ('tors_H.prod', 'torsion_H_production_deg', 'NE2(His94)–Zn–N2–H, mean of the production')):
            add(f'{b_}.{k_}', ca_[kk_], 'deg', 0, f'md/results/site_checks_solute.json:runs[{lab_}].change_at_N2.{kk_}', f'{L_}: torsion {m2_}')
    n2_ = s_['n2_configuration']
    add(b_ + '.n2.n_negative', n2_['count_negative'], None, 0, f'{SRC_}:n2_configuration.count_negative', f'{L_}: production frames with the sign at N2 negative (the configuration of start R)', check=q_['n2']['n_negative_production'], tol=0.01)
    add(b_ + '.n2.pct_negative', 100.0 * n2_['count_negative'] / n2_['n'], '%', 2, f'{SRC_}:n2_configuration.count_negative / n', f'{L_}: share of production frames with the sign at N2 negative', check=q_['n2']['pct_negative_production'])
    add(b_ + '.n2.n_changes', ev_['n2_volume']['from_the_release']['n_changes_of_sign'], None, 0, f'md/results/run_events.json:runs[{lab_}].n2_volume.from_the_release.n_changes_of_sign',
        f'{L_}: changes of the sign at N2 from the release to the end of the production', check=q_['n2']['n_changes_from_the_release'], tol=0.01)
    add(b_ + '.n2.n_flat', ev_['n2_volume']['production']['n_abs_below_0.2'], None, 0, f'md/results/run_events.json:runs[{lab_}].n2_volume.production.n_abs_below_0.2',
        f'{L_}: production frames with the signed volume at N2 below 0.2 in size (a nearly planar nitrogen)', check=q_['n2']['n_abs_volume_below_0p2_production'], tol=0.01)
    add(b_ + '.zn_n2_h.mean', ss_['nitrogen']['production']['ang_Zn_N2_H_mean_deg'], 'deg', 1, f'md/results/site_checks_solute.json:runs[{lab_}].nitrogen.production.ang_Zn_N2_H_mean_deg', f'{L_}: angle Zn–N2–H, mean over the production')
    # the ligand and the zinc site (means over the 10,000 solute frames; the zinc geometry of set 2e is put in)
    for k_, kk_, m2_, ref_ in (('c4sn2', 'ang_C4_S_N2_deg', 'angle C4–S–N2', _k['gate_GV']['ang_C4_S_N2_reference_deg']), ('tau_fold', 'tau_fold_deg', 'torsion τ_fold', _k['gate_GV']['tau_fold_reference_deg'])):
        add(f'{b_}.{k_}.mean', s_['stats'][kk_]['mean'], 'deg', 1, f'{SRC_}:stats.{kk_}.mean', f'{L_}: {m2_}, mean over the production')
        add(f'{b_}.{k_}.sd', s_['stats'][kk_]['sd'], 'deg', 1, f'{SRC_}:stats.{kk_}.sd', f'{L_}: {m2_}, standard deviation')
        add(f'{b_}.{k_}.minus_xtal', s_['stats'][kk_]['mean'] - ref_, 'deg', 1, f'{SRC_}:stats.{kk_}.mean − constants.gate_GV reference ({ref_}°)', f'{L_}: {m2_}, mean minus the value of 6RL9')
    tr_ = q_['third_reading']
    for k_, kk_, m2_, t3_ in (('zn_n2', 'd_Zn_N2_A', 'Zn–N2', 'Zn_N2'), ('zn_h94', 'd_Zn_H94NE2_A', 'Zn–NE2(His94)', 'Zn_H94NE2'), ('zn_h96', 'd_Zn_H96NE2_A', 'Zn–NE2(His96)', 'Zn_H96NE2'),
                              ('zn_h119', 'd_Zn_H119ND1_A', 'Zn–ND1(His119)', 'Zn_H119ND1')):
        add(f'{b_}.{k_}.mean', s_['stats'][kk_]['mean'], 'A', 3, f'{SRC_}:stats.{kk_}.mean', f'{L_}: distance {m2_}, mean over the production', check=tr_[t3_], tol=0.0006)
    add(b_ + '.zn_n2.sd', s_['stats']['d_Zn_N2_A']['sd'], 'A', 3, f'{SRC_}:stats.d_Zn_N2_A.sd', f'{L_}: distance Zn–N2, standard deviation')
    gv_ = s_['gate_GV']['GV3_zinc']['criteria']
    dev_ = max((gv_[a_]['deviation_deg'] for a_ in gv_ if a_.startswith('ang_')), key=abs)
    for a_ in gv_:
        if a_.startswith('ang_'):
            add(f'{b_}.ang.{a_[4:-4]}', gv_[a_]['mean_deg'], 'deg', 1, f'{SRC_}:gate_GV.GV3_zinc.criteria.{a_}.mean_deg', f'{L_}: mean angle {a_[4:-4].replace("_", "–")} at the zinc (6RL9: {gv_[a_]["reference_6RL9_deg"]}°)')
    add(b_ + '.ang.max_dev', dev_, 'deg', 1, f'{SRC_}:gate_GV.GV3_zinc.criteria.*.deviation_deg (largest in size)', f'{L_}: largest departure of a mean N–Zn–N angle from its value in 6RL9')
    add(b_ + '.s_n2_zn.mean', s_['stats']['ang_S_N2_Zn_deg']['mean'], 'deg', 1, f'{SRC_}:stats.ang_S_N2_Zn_deg.mean', f'{L_}: angle S–N2–Zn, mean over the production (117.6° in 6RL9)')
    add(b_ + '.so_zn.mean', s_['stats']['d_Zn_Osulfonyl_min_A']['mean'], 'A', 3, f'{SRC_}:stats.d_Zn_Osulfonyl_min_A.mean', f'{L_}: distance of the nearer sulfonyl oxygen from the zinc, mean over the production', check=tr_['Zn_Osulfonyl_min_mean'], tol=0.0006)
    fs_ = s_['fractions']['fifth_sulfonylO']
    add(b_ + '.so_zn.n_within', fs_['count'], None, 0, f'{SRC_}:fractions.fifth_sulfonylO.count', f'{L_}: solute frames with a sulfonyl oxygen within 2.6 Å of the zinc (of {fs_["n"]:,})', check=tr_['sulfonyl_O_within_2p6'], tol=0.01)
    add(b_ + '.so_zn.pct_within', 100.0 * fs_['count'] / fs_['n'], '%', 2, f'{SRC_}:fractions.fifth_sulfonylO.count / n', f'{L_}: share of solute frames with a sulfonyl oxygen within 2.6 Å of the zinc', check=100.0 * tr_['sulfonyl_O_within_2p6'] / fs_['n'])
    fa_ = s_['fractions']['fifth_any']
    add(b_ + '.any_o.n_within', fa_['count'], None, 0, f'{SRC_}:fractions.fifth_any.count', f'{L_}: all-atom frames with any oxygen within 2.6 Å of the zinc (of {fa_["n"]:,}); each is a sulfonyl oxygen', check=q_['all_atom']['any']['count_within_2p6'], tol=0.01)
    add(b_ + '.any_o.pct_within', 100.0 * fa_['count'] / fa_['n'], '%', 1, f'{SRC_}:fractions.fifth_any.count / n', f'{L_}: share of all-atom frames with any oxygen within 2.6 Å of the zinc', check=q_['all_atom']['any']['pct_within_2p6'])
    add(b_ + '.water_o.n_within', s_['fractions']['fifth_water']['count'], None, 0, f'{SRC_}:fractions.fifth_water.count', f'{L_}: all-atom frames with a water oxygen within 2.6 Å of the zinc', check=q_['all_atom']['water']['count_within_2p6'], tol=0.01)
    add(b_ + '.protein_o.n_within', s_['fractions']['fifth_proteinO']['count'], None, 0, f'{SRC_}:fractions.fifth_proteinO.count', f'{L_}: all-atom frames with a protein oxygen within 2.6 Å of the zinc', check=q_['all_atom']['protein']['count_within_2p6'], tol=0.01)
    add(b_ + '.zn_oe1_min', s_['stats']['d_Zn_E106OE1_A']['min'], 'A', 2, f'{SRC_}:stats.d_Zn_E106OE1_A.min', f'{L_}: smallest distance Zn–OE1(Glu106) in a solute frame')
    # health
    pl_ = h_['production_state_log']
    add(b_ + '.temp.mean', pl_['temperature_K']['mean'], 'K', 2, f'md/results/run_health.json:runs[{lab_}].production_state_log.temperature_K.mean', f'{L_}: kinetic temperature of the state log, mean over the production (thermostat at 300 K, time step 4 fs)')
    add(b_ + '.temp.sd', pl_['temperature_K']['sd'], 'K', 2, f'md/results/run_health.json:runs[{lab_}].production_state_log.temperature_K.sd', f'{L_}: kinetic temperature, standard deviation')
    add(b_ + '.density.mean', pl_['density_g_per_cm3']['mean'], 'g cm-3', 4, f'md/results/run_health.json:runs[{lab_}].production_state_log.density_g_per_cm3.mean', f'{L_}: density, mean over the production')
    add(b_ + '.ca_rmsd.mean', s_['stats']['ca_rmsd_A']['mean'], 'A', 2, f'{SRC_}:stats.ca_rmsd_A.mean', f'{L_}: RMSD of the 258 Cα atoms from 6RL9, mean over the production', check=tr_['ca_rmsd_mean'], tol=0.0006)
    add(b_ + '.ca_rmsd.last_half', s_['stats']['ca_rmsd_A']['mean_last_half'], 'A', 2, f'{SRC_}:stats.ca_rmsd_A.mean_last_half', f'{L_}: Cα RMSD, mean of the last half')
    add(b_ + '.ca_rmsd.max', s_['stats']['ca_rmsd_A']['max'], 'A', 2, f'{SRC_}:stats.ca_rmsd_A.max', f'{L_}: Cα RMSD, largest single frame', check=tr_['ca_rmsd_max'], tol=0.0006)
    add(b_ + '.rmsf_max', s_['rmsf']['max_A'], 'A', 2, f'{SRC_}:rmsf.max_A', f'{L_}: largest Cα RMSF (residue {s_["rmsf"]["max_residue_number_deposited"]})')
    # checks on the raw frames (md/production/site_checks.py): the docked rotated pose, His64, water of the site
    lr_ = ss_['ligand_against_the_docked_rotated_pose']
    add(b_ + '.rot_pose.min_rmsd', lr_['production_smallest_A'], 'A', 2, f'md/results/site_checks_solute.json:runs[{lab_}].ligand_against_the_docked_rotated_pose.production_smallest_A', f'{L_}: smallest RMSD of the ligand to the docked rotated pose in the production')
    add(b_ + '.rot_pose.n_within_1', lr_['production_frames_within_1.0_A'], None, 0, f'md/results/site_checks_solute.json:runs[{lab_}].ligand_against_the_docked_rotated_pose.production_frames_within_1.0_A', f'{L_}: production frames with the ligand within 1.0 Å of the docked rotated pose')
    dp_ = ss_['departures']
    add(b_ + '.zn_oop_his119', dp_['zinc_out_of_the_plane_of_the_His119_ring_A'], 'A', 2, f'md/results/site_checks_solute.json:runs[{lab_}].departures.zinc_out_of_the_plane_of_the_His119_ring_A', f'{L_}: zinc out of the plane of the His119 ring, mean over the production ({_XT["zn_oop_his119"]:.2f} Å in 6RL9)')
    add(b_ + '.s_c4.mean', dp_['d_S_C4_A'], 'A', 3, f'md/results/site_checks_solute.json:runs[{lab_}].departures.d_S_C4_A', f'{L_}: bond S–C4, mean over the production (1.771 Å in 6RL9)')
    add(b_ + '.o1_s_c4.mean', dp_['ang_O1_S_C4_deg'], 'deg', 1, f'md/results/site_checks_solute.json:runs[{lab_}].departures.ang_O1_S_C4_deg', f'{L_}: angle O1–S–C4, mean over the production (106.6° in 6RL9)')
    add(b_ + '.his64.deep_pct', ss_['his64']['deep_in_share_20_to_100_ns']['7.0'], '%', 1, f'md/results/site_checks_solute.json:runs[{lab_}].his64.deep_in_share_20_to_100_ns["7.0"]',
        f'{L_}: share of the production from 20 to 100 ns with His64 deep in (Zn–NE2 at most 7.0 Å and χ1 between 0 and 120°; the limit was chosen after looking at the data)')
    cw_ = sf_['crystal_water_places']
    for w_ in ('423', '516', '553', '568'):
        add(f'{b_}.w{w_}.pct', cw_['production'][w_], '%', 1, f'md/results/site_checks_full.json:runs[{lab_}].crystal_water_places.production["{w_}"]',
            f'{L_}: all-atom frames of the production with a water oxygen within 1.4 Å of the place of crystal water {w_}')
    add(b_ + '.w516.pct_20_100', cw_['production_20_to_100_ns']['516'], '%', 1, f'md/results/site_checks_full.json:runs[{lab_}].crystal_water_places.production_20_to_100_ns["516"]', f'{L_}: the place of crystal water 516 occupied, 20 to 100 ns')
    add(b_ + '.w_n2.pct', sf_['water_donating_to_N2']['share_production'], '%', 1, f'md/results/site_checks_full.json:runs[{lab_}].water_donating_to_N2.share_production', f'{L_}: all-atom frames of the production with a water donating a hydrogen bond to N2')

# the candidate run: what its rotated-like frames are
_q, _e = D17M['runs']['R1_n2fit'], EV_['bonded_2e_n2fit/R1']
_ES = 'md/results/run_events.json:runs[bonded_2e_n2fit/R1]'
_vr = _e['production_visits_by_assigned_state']['R']
add('md.run.R1_n2fit.visits_R.n', _vr['n_visits'], None, 0, f'{_ES}.production_visits_by_assigned_state.R.n_visits', 'candidate terms at N2: visits to the rotated-like state in the production', check=_q['production_visits_R']['n'], tol=0.01)
add('md.run.R1_n2fit.visits_R.longest_ns', _vr['longest_visit_frames'] / 100.0, 'ns', 2, f'{_ES}.production_visits_by_assigned_state.R.longest_visit_frames (10 ps each)', 'candidate terms at N2: longest visit to the rotated-like state (none reaches 1 ns)',
    check=_q['production_visits_R']['longest_ns'])
add('md.run.R1_n2fit.visits_R.n_long', _vr['n_long_visits'], None, 0, f'{_ES}.production_visits_by_assigned_state.R.n_long_visits', 'candidate terms at N2: visits of at least 1 ns to the rotated-like state', check=_q['production_visits_R']['n_long'], tol=0.01)
_sn = _e['state_and_n2_sign']['production']
add('md.run.R1_n2fit.n_R_sign_pos', _sn['R_with_sign_positive'], None, 0, f'{_ES}.state_and_n2_sign.production.R_with_sign_positive', 'candidate terms at N2: rotated-like production frames with the sign at N2 of the crystal pose', check=_q['n2']['n_R_with_sign_positive_production'], tol=0.01)
add('md.run.R1_n2fit.n_C_sign_neg', _sn['C_with_sign_negative'], None, 0, f'{_ES}.state_and_n2_sign.production.C_with_sign_negative', 'candidate terms at N2: crystal-like production frames with the sign at N2 of start R', check=_q['n2']['n_C_with_sign_negative_production'], tol=0.01)
add('md.run.R1_n2fit.rot_pose.pct_within_1', 100.0 * SS_['bonded_2e_n2fit/R1']['ligand_against_the_docked_rotated_pose']['production_frames_within_1.0_A'] / 10000, '%', 2,
    'md/results/site_checks_solute.json:runs[bonded_2e_n2fit/R1].ligand_against_the_docked_rotated_pose.production_frames_within_1.0_A / 10000', 'candidate terms at N2: share of production frames within 1.0 Å of the docked rotated pose')
for i_, bl_ in enumerate(_e['production_blocks_of_20_ns']):
    add(f'md.run.R1_n2fit.block{i_ + 1}.n_R', bl_['R'], None, 0, f'{_ES}.production_blocks_of_20_ns[{i_}].R', f'candidate terms at N2: rotated-like frames of the block {int(bl_["from_ns"])} to {int(bl_["to_ns"])} ns (of 2,000)',
        check=_q['blocks_20ns'][i_]['R'], tol=0.01)
add('md.run.R1_n2fit.rmsd.mean_C', _e['production_by_state']['C']['lig_rmsd_mean_A'], 'A', 2, f'{_ES}.production_by_state.C.lig_rmsd_mean_A', 'candidate terms at N2: mean ligand RMSD of the crystal-like frames', check=_q['rmsd']['mean_crystal_like'])
_e2 = EV_['bonded_2s/R1']
add('md.run.R1_2s.n_R_released', _e2['state_and_n2_sign']['from_the_release']['R'], None, 0, 'md/results/run_events.json:runs[bonded_2s/R1].state_and_n2_sign.from_the_release.R', 'zinc set 2s: rotated-like frames from the release to the end')
add('md.run.R1_2s.n_R_released_sign_pos', _e2['state_and_n2_sign']['from_the_release']['R_with_sign_positive'], None, 0, 'md/results/run_events.json:runs[bonded_2s/R1].state_and_n2_sign.from_the_release.R_with_sign_positive',
    'zinc set 2s: of these, frames with the sign at N2 of the crystal pose')

# ----------------------------------------------------------------------------------- groups of runs
_RS = 'first_program_summary.json of the runs named'
_g = D17M['ranges']
for g_, ids_, lab_ in (('X', G_X_, 'X1 to X3'), ('N', G_N_, 'N1 to N3'), ('R', G_R_, 'R1 to R3'), ('XN', G_X_ + G_N_, 'the six runs from X and N')):
    _rng(f'md.grp.{g_}.rmsd_last_half', ids_, lambda i_: S1_[i_]['stats']['lig_rmsd_A']['mean_last_half'], 'A', 2, f'{_RS}: stats.lig_rmsd_A.mean_last_half',
         f'mean ligand RMSD of the last half over {lab_}', checks=_g['rmsd_last_half'][g_])
for g_, ids_, lab_ in (('XN', G_X_ + G_N_, 'the six runs from X and N'), ('R', G_R_, 'R1 to R3'), ('nine', NINE_, 'the nine runs')):
    _rng(f'md.grp.{g_}.hb_pct', ids_, lambda i_: 100.0 * S1_[i_]['fractions']['hb_N2H_T199OG1']['count'] / S1_[i_]['fractions']['hb_N2H_T199OG1']['n'], '%', 2,
         f'{_RS}: fractions.hb_N2H_T199OG1', f'share of production frames with the hydrogen bond to Thr199 over {lab_}', checks=_g['hb_pct'][g_])
    _rng(f'md.grp.{g_}.c4sn2_minus_xtal', ids_, lambda i_: S1_[i_]['stats']['ang_C4_S_N2_deg']['mean'] - _k['gate_GV']['ang_C4_S_N2_reference_deg'], 'deg', 1,
         f'{_RS}: stats.ang_C4_S_N2_deg.mean − 108.8°', f'mean C4–S–N2 minus the value of 6RL9 over {lab_}')
    _rng(f'md.grp.{g_}.c4sn2', ids_, lambda i_: S1_[i_]['stats']['ang_C4_S_N2_deg']['mean'], 'deg', 1, f'{_RS}: stats.ang_C4_S_N2_deg.mean', f'mean C4–S–N2 over {lab_}')
    _rng(f'md.grp.{g_}.tau_fold_minus_xtal', ids_, lambda i_: S1_[i_]['stats']['tau_fold_deg']['mean'] - _k['gate_GV']['tau_fold_reference_deg'], 'deg', 1,
         f'{_RS}: stats.tau_fold_deg.mean − 132.1°', f'mean τ_fold minus the value of 6RL9 over {lab_}')
    _rng(f'md.grp.{g_}.tau_fold', ids_, lambda i_: S1_[i_]['stats']['tau_fold_deg']['mean'], 'deg', 1, f'{_RS}: stats.tau_fold_deg.mean', f'mean τ_fold over {lab_}')
    _rng(f'md.grp.{g_}.his64_deep_pct', ids_, lambda i_: SS_[LAB_[i_]]['his64']['deep_in_share_20_to_100_ns']['7.0'], '%', 1, 'md/results/site_checks_solute.json:runs[].his64.deep_in_share_20_to_100_ns["7.0"]',
         f'share of 20 to 100 ns with His64 deep in over {lab_}')
    _rng(f'md.grp.{g_}.w516_pct_20_100', ids_, lambda i_: SF_[LAB_[i_]]['crystal_water_places']['production_20_to_100_ns']['516'], '%', 1, 'md/results/site_checks_full.json:runs[].crystal_water_places.production_20_to_100_ns["516"]',
         f'the place of crystal water 516 occupied, 20 to 100 ns, over {lab_}')
add('md.grp.XN.rmsd_max', max(S1_[i_]['stats']['lig_rmsd_A']['max'] for i_ in G_X_ + G_N_), 'A', 2, f'{_RS}: stats.lig_rmsd_A.max (largest over the six)', 'largest ligand RMSD of a single frame in the six runs from X and N (in X2)',
    check=_g['rmsd_max']['XN'][1])
for k_, fn_, c_ in (('rule_ps', lambda i_: int(CT_[LAB_[i_]]['conversion']['time_from_release_ps']), _g['rule_ps']['R']),
                    ('sign_ps', lambda i_: int(EV_[LAB_[i_]]['first_released_frame_n2_sign_not_of_start']['time_from_release_ps']), _g['sign_ps']['R']),
                    ('left_ps', lambda i_: int(EV_[LAB_[i_]]['first_released_frame_not_in_state_of_start']['time_from_release_ps']), _g['not_rotated_ps']['R'])):
    m_ = {'rule_ps': 'time of the plan\'s rule from the release over R1 to R3', 'sign_ps': 'time from the release of the first frame with the sign at N2 positive over R1 to R3',
          'left_ps': 'time from the release of the first frame that is not rotated-like over R1 to R3'}[k_]
    _rng(f'md.grp.R.{k_}', G_R_, fn_, 'ps', 0, 'md/results/conversion_times.json; md/results/run_events.json', m_, checks=c_, tol=0.01)
assert all(EV_[LAB_[i_]]['first_released_frame_n2_sign_not_of_start']['stage'] == 4 and EV_[LAB_[i_]]['first_released_frame_not_in_state_of_start']['stage'] == 4 for i_ in G_R_ + ['R1_2s', 'R1_n2fit'])
add('md.grp.R.event_stage', 4, None, 0, 'md/results/run_events.json:runs[].first_released_frame_n2_sign_not_of_start.stage, first_released_frame_not_in_state_of_start.stage',
    'equilibration stage of the first frame with the sign at N2 of the crystal pose and of the first frame that is not rotated-like, in all five runs from start R (protein backbone held with 250 kJ mol⁻¹ nm⁻²)')
for g_, ids_, lab_ in (('nine', NINE_, 'the nine runs'), ('R', G_R_, 'R1 to R3'), ('XN', G_X_ + G_N_, 'the six runs from X and N')):
    add(f'md.grp.{g_}.n_R_frames', sum(S1_[i_]['orientation']['count_R'] for i_ in ids_), None, 0, f'{_RS}: orientation.count_R (sum)', f'rotated-like frames in the production of {lab_}',
        check=_g['n_rotated_like_production'][g_], tol=0.01)
    add(f'md.grp.{g_}.n_frames', sum(S1_[i_]['run']['n_frames_solute'] for i_ in ids_), None, 0, f'{_RS}: run.n_frames_solute (sum)', f'production frames of {lab_}', check=_g['n_frames_production'][g_], tol=0.01)
_rng('md.grp.nine.n2_flat', NINE_, lambda i_: EV_[LAB_[i_]]['n2_volume']['production']['n_abs_below_0.2'], None, 0, 'md/results/run_events.json:runs[].n2_volume.production.n_abs_below_0.2',
     'number of production frames with the signed volume at N2 below 0.2 in size over the nine runs (of 10,000 each)', checks=_g['n_abs_volume_below_0p2']['nine'], tol=0.01)
_ZE = NINE_ + ['X0']
_zs = lambda i_: (SX0_['2e'] if i_ == 'X0' else S1_[i_])   # noqa: E731
_rng('md.grp.2e.zn_n2', _ZE, lambda i_: _zs(i_)['stats']['d_Zn_N2_A']['mean'], 'A', 3, f'{_RS} and of bonded_2e/X0: stats.d_Zn_N2_A.mean', 'mean Zn–N2 distance over the nine runs and X0 of set 2e (1.931 Å in 6RL9; put in by the parameter set)')
_rng('md.grp.2e.ang_max_dev', _ZE, lambda i_: max(abs(v_['deviation_deg']) for a_, v_ in _zs(i_)['gate_GV']['GV3_zinc']['criteria'].items() if a_.startswith('ang_')), 'deg', 1,
     f'{_RS} and of bonded_2e/X0: gate_GV.GV3_zinc.criteria.*.deviation_deg', 'departure in size of a mean N–Zn–N angle from 6RL9 over the nine runs and X0 of set 2e')
_rng('md.grp.nine.so_zn_pct', NINE_, lambda i_: 100.0 * S1_[i_]['fractions']['fifth_sulfonylO']['count'] / S1_[i_]['fractions']['fifth_sulfonylO']['n'], '%', 2, f'{_RS}: fractions.fifth_sulfonylO',
     'share of solute frames with a sulfonyl oxygen within 2.6 Å of the zinc over the nine runs')
_rng('md.grp.nine.any_o_pct', NINE_, lambda i_: 100.0 * S1_[i_]['fractions']['fifth_any']['count'] / S1_[i_]['fractions']['fifth_any']['n'], '%', 1, f'{_RS}: fractions.fifth_any',
     'share of all-atom frames with any oxygen within 2.6 Å of the zinc over the nine runs', checks=_g['pct_any_O_within_2p6']['nine'])
add('md.grp.all.n_water_or_protein_o', sum(S1_[i_]['fractions']['fifth_water']['count'] + S1_[i_]['fractions']['fifth_proteinO']['count'] for i_ in ELEVEN_), None, 0, f'{_RS}: fractions.fifth_water, fifth_proteinO (sum over eleven runs)',
    'all-atom frames of the eleven runs with a water or protein oxygen within 2.6 Å of the zinc', check=_g['n_water_or_protein_O_within_2p6'], tol=0.01)
add('md.grp.all.zn_oe1_min', min([S1_[i_]['stats']['d_Zn_E106OE1_A']['min'] for i_ in ELEVEN_] + [SX0_[t_]['stats']['d_Zn_E106OE1_A']['min'] for t_ in SX0_]), 'A', 2,
    f'{_RS} and of the two X0 runs: stats.d_Zn_E106OE1_A.min', 'smallest distance Zn–OE1(Glu106) in any solute frame of the 13 runs')
_rng('md.grp.nine.s_n2_zn', NINE_, lambda i_: S1_[i_]['stats']['ang_S_N2_Zn_deg']['mean'], 'deg', 1, f'{_RS}: stats.ang_S_N2_Zn_deg.mean', 'mean S–N2–Zn over the nine runs (117.6° in 6RL9)')
_rng('md.grp.nine.zn_oop_his119', NINE_, lambda i_: SS_[LAB_[i_]]['departures']['zinc_out_of_the_plane_of_the_His119_ring_A'], 'A', 2, 'md/results/site_checks_solute.json:runs[].departures.zinc_out_of_the_plane_of_the_His119_ring_A',
     f'zinc out of the plane of the His119 ring over the nine runs ({_XT["zn_oop_his119"]:.2f} Å in 6RL9)')
_rng('md.grp.nine.s_c4', NINE_, lambda i_: SS_[LAB_[i_]]['departures']['d_S_C4_A'], 'A', 3, 'md/results/site_checks_solute.json:runs[].departures.d_S_C4_A', 'mean S–C4 over the nine runs (1.771 Å in 6RL9)')
_rng('md.grp.nine.o1_s_c4', NINE_, lambda i_: SS_[LAB_[i_]]['departures']['ang_O1_S_C4_deg'], 'deg', 1, 'md/results/site_checks_solute.json:runs[].departures.ang_O1_S_C4_deg', 'mean O1–S–C4 over the nine runs (106.6° in 6RL9)')
_rng('md.grp.nine.zn_n2_h', NINE_, lambda i_: SS_[LAB_[i_]]['nitrogen']['production']['ang_Zn_N2_H_mean_deg'], 'deg', 1, 'md/results/site_checks_solute.json:runs[].nitrogen.production.ang_Zn_N2_H_mean_deg', 'mean Zn–N2–H over the nine runs')
add('md.grp.nine.rot_pose_min_rmsd', min(SS_[LAB_[i_]]['ligand_against_the_docked_rotated_pose']['production_smallest_A'] for i_ in NINE_), 'A', 2,
    'md/results/site_checks_solute.json:runs[].ligand_against_the_docked_rotated_pose.production_smallest_A (smallest over the nine)', 'nearest approach of the ligand to the docked rotated pose, as RMSD, in the production of the nine runs')
_rng('md.grp.nine.w553_pct', NINE_, lambda i_: SF_[LAB_[i_]]['crystal_water_places']['production']['553'], '%', 1, 'md/results/site_checks_full.json:runs[].crystal_water_places.production["553"]', 'the place of crystal water 553 occupied over the nine runs')
_rng('md.grp.nine.w_n2_pct', NINE_, lambda i_: SF_[LAB_[i_]]['water_donating_to_N2']['share_production'], '%', 1, 'md/results/site_checks_full.json:runs[].water_donating_to_N2.share_production', 'share of all-atom frames with a water donating a hydrogen bond to N2 over the nine runs')
for w_ in ('423', '516'):
    _rng(f'md.grp.nine.w{w_}_pct', NINE_, lambda i_: SF_[LAB_[i_]]['crystal_water_places']['production'][w_], '%', 1, f'md/results/site_checks_full.json:runs[].crystal_water_places.production["{w_}"]',
         f'the place of crystal water {w_} occupied over the whole production of the nine runs')
_HS = 'md/results/run_health.json:runs[].production_state_log'
_rng('md.grp.eleven.temp', ELEVEN_, lambda i_: HR_[LAB_[i_]]['production_state_log']['temperature_K']['mean'], 'K', 2, f'{_HS}.temperature_K.mean', 'mean kinetic temperature of the production over the eleven runs (thermostat at 300 K; time step 4 fs)')
_rng('md.grp.eleven.density', ELEVEN_, lambda i_: HR_[LAB_[i_]]['production_state_log']['density_g_per_cm3']['mean'], 'g cm-3', 4, f'{_HS}.density_g_per_cm3.mean', 'mean density of the production over the eleven runs')
_rng('md.grp.eleven.ca_rmsd', ELEVEN_, lambda i_: S1_[i_]['stats']['ca_rmsd_A']['mean'], 'A', 2, f'{_RS}: stats.ca_rmsd_A.mean', 'mean Cα RMSD from 6RL9 over the eleven runs',
     checks=[min(D17M['runs'][i_]['third_reading']['ca_rmsd_mean'] for i_ in ELEVEN_), max(D17M['runs'][i_]['third_reading']['ca_rmsd_mean'] for i_ in ELEVEN_)], tol=0.0006)
add('md.grp.eleven.ca_rmsd_max', max(S1_[i_]['stats']['ca_rmsd_A']['max'] for i_ in ELEVEN_), 'A', 2, f'{_RS}: stats.ca_rmsd_A.max (largest over eleven runs)', 'largest Cα RMSD of a single frame in the eleven runs',
    check=max(D17M['runs'][i_]['third_reading']['ca_rmsd_max'] for i_ in ELEVEN_), tol=0.0006)
_tq = [st_['log_mean_temperature_K'] for h_ in _hl for st_ in h_['equilibration_stages'] if st_['timestep_fs'] == 2.0 and st_['ensemble'] == 'NPT']
add('md.grp.temp_2fs.min', min(_tq), 'K', 1, 'md/results/run_health.json:runs[].equilibration_stages[].log_mean_temperature_K (NPT stages at 2 fs)', 'smallest mean kinetic temperature of an equilibration stage at the 2 fs step (13 runs)')
add('md.grp.temp_2fs.max', max(_tq), 'K', 1, 'md/results/run_health.json:runs[].equilibration_stages[].log_mean_temperature_K (NPT stages at 2 fs)', 'largest mean kinetic temperature of an equilibration stage at the 2 fs step (13 runs)')
_r1 = CAMP_['rules_per_model']['bonded_2e']
add('md.rule1.n_left', _r1['rule1']['n_X_left_crystal'], None, 0, 'md/results/campaign/first_program/campaign.json:rules_per_model.bonded_2e.rule1.n_X_left_crystal', 'rule 1: crystal-start runs that left the crystal pose (of 3)',
    check=sum(1 for i_ in G_X_ if D17M['runs'][i_]['rmsd']['mean_last_half'] > 2.0), tol=0.01)
add('md.rule3.n_R_converted', _r1['n_converted']['R'], None, 0, 'md/results/campaign/first_program/campaign.json:rules_per_model.bonded_2e.n_converted.R', 'rule 3: runs from start R with the verdict CONVERTED (of 3)',
    check=sum(1 for i_ in G_R_ if D17M['runs'][i_]['verdict'] == 'CONVERTED'), tol=0.01)
add('md.rule3.n_XN_stayed', _r1['n_stayed']['X'] + _r1['n_stayed']['N'], None, 0, 'md/results/campaign/first_program/campaign.json:rules_per_model.bonded_2e.n_stayed', 'rule 3: runs from starts X and N with the verdict STAYED (of 6)',
    check=sum(1 for i_ in G_X_ + G_N_ if D17M['runs'][i_]['verdict'] == 'STAYED'), tol=0.01)

# the rehearsal run from start N, cited as a limit of "stayed" (a test of the pipeline; no result and in no count)
_tn = J(MDR_ / 'test/rehearsal_2e/runs/TN0/analysis/summary.json')
add('md.rehearsal.TN0.n_R', _tn['orientation']['count_R'], None, 0, 'md/test/rehearsal_2e/runs/TN0/analysis/summary.json:orientation.count_R',
    'rehearsal run of 5 ns from start N with the topology of run N1 (a test of the pipeline, no result): rotated-like frames', check=D17M['rehearsal_TN0']['n_rotated_like'], tol=0.01)
add('md.rehearsal.TN0.n_frames', _tn['run']['n_frames_solute'], None, 0, 'md/test/rehearsal_2e/runs/TN0/analysis/summary.json:run.n_frames_solute', 'the same rehearsal run: solute frames',
    check=D17M['rehearsal_TN0']['n_frames'], tol=0.01)
add('md.rehearsal.TN0.length_ns', _tn['run']['seen_ns'], 'ns', 0, 'md/test/rehearsal_2e/runs/TN0/analysis/summary.json:run.seen_ns', 'the same rehearsal run: length')

# ================================================================================ the tests of the force field
# G-L: the repaired ligand terms against the DFT of the paper (implicit water), anion
for st_, k_, lab_ in (('before', LIGR_['species']['anion']['gate_before']['OBC2'], 'GAFF2 as assigned'), ('after', LIGR_['species']['anion']['results']['OBC2']['gate'], 'repaired terms')):
    src_ = 'md/ligand/repair.json:species.anion.' + ('gate_before.OBC2' if st_ == 'before' else 'results.OBC2.gate')
    add(f'md.gate.GL.{st_}.csn', k_['d_csn'], 'deg', 2, f'{src_}.d_csn', f'gate G-L, sulfonamidate anion, {lab_}, implicit water: C–S–N at the force-field minimum minus DFT')
    add(f'md.gate.GL.{st_}.sn', k_['d_sn'], 'A', 3, f'{src_}.d_sn', f'gate G-L, sulfonamidate anion, {lab_}, implicit water: S–N at the force-field minimum minus DFT')
    add(f'md.gate.GL.{st_}.scan', k_['scan_max_abs_dev'], 'kcal mol-1', 2, f'{src_}.scan_max_abs_dev', f'gate G-L, sulfonamidate anion, {lab_}, implicit water: largest deviation of the torsion profile from the DFT scan')
    add(f'md.gate.GL.{st_}.outcome', 'passed' if k_['pass'] else 'failed', None, None, f'{src_}.pass', f'gate G-L, sulfonamidate anion, {lab_}: outcome')
for k_, v_, u_, p_, m_ in (('csn', 2.0, 'deg', 0, 'C–S–N within this of DFT at the force-field minimum'), ('sn', 0.02, 'A', 2, 'S–N within this of DFT'), ('scan', 0.5, 'kcal mol-1', 1, 'torsion profile within this of the DFT scan')):
    add(f'md.gate.GL.limit.{k_}', v_, u_, p_, 'md/ligand_eval.py:gate (constants); md/PLAN.md section 4, stage S1', f'gate G-L, criterion: {m_}')
# G-Z: the zinc site in the first 100 ps of X0 (heavy atoms of the solute still held)
_cz0 = GZ_['2e']['criteria']
for t_ in ('2e', '2s'):
    g_, d_ = GZ_[t_]['gate_GZ'], D17M['X0'][t_]['gate_GZ']
    src_ = f'md/production/gate_GZ_bonded_{t_}.json:gate_GZ'
    add(f'md.gate.GZ.{t_}.zn_n2', g_['a_Zn_N2']['mean_A'], 'A', 3, f'{src_}.a_Zn_N2.mean_A', f'gate G-Z, zinc set {t_}: mean Zn–N2 of the 10 frames of the first 100 ps of X0', check=d_['means']['d_Zn_N2_A'], tol=0.0006)
    for kk_, nm_ in (('d_Zn_H94NE2_A', 'h94'), ('d_Zn_H96NE2_A', 'h96'), ('d_Zn_H119ND1_A', 'h119')):
        add(f'md.gate.GZ.{t_}.zn_{nm_}', g_['b_Zn_NHis'][kk_]['mean_A'], 'A', 3, f'{src_}.b_Zn_NHis.{kk_}.mean_A', f'gate G-Z, zinc set {t_}: mean {kk_[2:-2].replace("_", "–")}', check=d_['means'][kk_], tol=0.0006)
    for kk_ in g_['c_N_Zn_N']:
        add(f'md.gate.GZ.{t_}.{kk_[:-4]}', g_['c_N_Zn_N'][kk_]['mean_deg'], 'deg', 1, f'{src_}.c_N_Zn_N.{kk_}.mean_deg', f'gate G-Z, zinc set {t_}: mean angle {kk_[4:-4].replace("_", "–")} (within 10° of 6RL9 asked)', check=d_['means'][kk_], tol=0.006)
    add(f'md.gate.GZ.{t_}.ang_max_dev', max((g_['c_N_Zn_N'][kk_]['mean_deg'] - _cz0['c_references_6RL9_deg'][kk_] for kk_ in g_['c_N_Zn_N']), key=abs), 'deg', 1,
        f'{src_}.c_N_Zn_N.*.mean_deg − criteria.c_references_6RL9_deg (largest in size)', f'gate G-Z, zinc set {t_}: largest departure of a mean N–Zn–N angle from its value in 6RL9',
        check=max((d_['means'][kk_] - _cz0['c_references_6RL9_deg'][kk_] for kk_ in g_['c_N_Zn_N']), key=abs), tol=0.006)
    add(f'md.gate.GZ.{t_}.so_n', g_['d_fifth_ligand']['solute_frames']['count'], None, 0, f'{src_}.d_fifth_ligand.solute_frames.count', f'gate G-Z, zinc set {t_}: frames with a sulfonyl oxygen within 2.6 Å of the zinc (of 10)',
        check=d_['sulfonyl_O_within_2p6'], tol=0.01)
    add(f'md.gate.GZ.{t_}.outcome', GZ_[t_]['gate_GZ']['verdict'].lower(), None, None, f'{src_}.verdict', f'gate G-Z, zinc set {t_}: outcome (read on frames in which the zinc site is still held)')
add('md.gate.GZ.n_frames', GZ_['2e']['gate_GZ']['n_solute_frames'], None, 0, 'md/production/gate_GZ_bonded_2e.json:gate_GZ.n_solute_frames', 'gate G-Z: solute frames read (first 100 ps after the minimisation of run X0)', check=D17M['X0']['2e']['gate_GZ']['n_frames'], tol=0.01)
_cz = GZ_['2e']['criteria']
add('md.gate.GZ.limit.zn_n2', _cz['a_window_A'], 'A', 2, 'md/production/gate_GZ_bonded_2e.json:criteria.a_window_A', 'gates G-Z and G-V, criterion: window of the mean Zn–N2 distance')
add('md.gate.GZ.limit.zn_his', _cz['b_window_A'], 'A', 2, 'md/production/gate_GZ_bonded_2e.json:criteria.b_window_A', 'gates G-Z and G-V, criterion: window of each mean Zn–N(His) distance')
add('md.gate.GZ.limit.angle', _cz['c_tolerance_deg'], 'deg', 0, 'md/production/gate_GZ_bonded_2e.json:criteria.c_tolerance_deg', 'gates G-Z and G-V, criterion: each mean N–Zn–N angle within this of its value in 6RL9')
add('md.gate.GZ.limit.fifth_pct', 100.0 * GZ_['2e']['gate_GZ']['d_fifth_ligand']['limit_max_fraction'], '%', 0, 'md/production/gate_GZ_bonded_2e.json:gate_GZ.d_fifth_ligand.limit_max_fraction',
    'gates G-Z and G-V, criterion: an oxygen within 2.6 Å of the zinc in at most this share of the frames')
# G-V: validation run X0 (5 ns from the crystal pose), both zinc sets
for t_ in ('2e', '2s'):
    g_ = SX0_[t_]['gate_GV']
    src_ = f'md/results/runs/bonded_{t_}/X0/first_program_summary.json:gate_GV'
    add(f'md.gate.GV.{t_}.rmsd', g_['GV1_lig_rmsd']['mean_A'], 'A', 2, f'{src_}.GV1_lig_rmsd.mean_A', f'gate G-V, zinc set {t_}: mean ligand RMSD to the crystal pose over the 5 ns of X0 (at most 1.0 Å asked)')
    add(f'md.gate.GV.{t_}.hb_n', g_['GV2_hbond']['count'], None, 0, f'{src_}.GV2_hbond.count', f'gate G-V, zinc set {t_}: frames with the hydrogen bond N2–H···OG1(Thr199) (of {g_["GV2_hbond"]["n"]})')
    add(f'md.gate.GV.{t_}.hb_pct', 100.0 * g_['GV2_hbond']['count'] / g_['GV2_hbond']['n'], '%', 1, f'{src_}.GV2_hbond.count / n', f'gate G-V, zinc set {t_}: share of frames with that hydrogen bond (at least 80 % asked)')
    c_ = g_['GV3_zinc']['criteria']
    add(f'md.gate.GV.{t_}.zn_n2', c_['d_Zn_N2_A']['mean_A'], 'A', 3, f'{src_}.GV3_zinc.criteria.d_Zn_N2_A.mean_A', f'gate G-V, zinc set {t_}: mean Zn–N2')
    for kk_, nm_ in (('d_Zn_H94NE2_A', 'h94'), ('d_Zn_H96NE2_A', 'h96'), ('d_Zn_H119ND1_A', 'h119')):
        add(f'md.gate.GV.{t_}.zn_{nm_}', c_[kk_]['mean_A'], 'A', 3, f'{src_}.GV3_zinc.criteria.{kk_}.mean_A', f'gate G-V, zinc set {t_}: mean {kk_[2:-2].replace("_", "–")}')
    add(f'md.gate.GV.{t_}.ang_max_dev', max((v_['deviation_deg'] for a_, v_ in c_.items() if a_.startswith('ang_')), key=abs), 'deg', 1, f'{src_}.GV3_zinc.criteria.*.deviation_deg (largest in size)',
        f'gate G-V, zinc set {t_}: largest departure of a mean N–Zn–N angle from 6RL9 (within 10° asked)')
    fl_ = c_['fifth_ligand']
    add(f'md.gate.GV.{t_}.so_n', fl_['solute_trajectory']['count'], None, 0, f'{src_}.GV3_zinc.criteria.fifth_ligand.solute_trajectory.count', f'gate G-V, zinc set {t_}: solute frames with a sulfonyl oxygen within 2.6 Å of the zinc (of {fl_["solute_trajectory"]["n"]})')
    add(f'md.gate.GV.{t_}.so_pct', 100.0 * fl_['solute_trajectory']['count'] / fl_['solute_trajectory']['n'], '%', 1, f'{src_}.GV3_zinc.criteria.fifth_ligand.solute_trajectory', f'gate G-V, zinc set {t_}: share of solute frames with a sulfonyl oxygen within 2.6 Å')
    add(f'md.gate.GV.{t_}.any_n', fl_['full_trajectory']['count'], None, 0, f'{src_}.GV3_zinc.criteria.fifth_ligand.full_trajectory.count', f'gate G-V, zinc set {t_}: all-atom frames with any oxygen within 2.6 Å of the zinc (of {fl_["full_trajectory"]["n"]})',
        check=D17M['X0'][t_]['all_atom']['any']['count_within_2p6'], tol=0.01)
    add(f'md.gate.GV.{t_}.any_pct', 100.0 * fl_['full_trajectory']['count'] / fl_['full_trajectory']['n'], '%', 0, f'{src_}.GV3_zinc.criteria.fifth_ligand.full_trajectory', f'gate G-V, zinc set {t_}: share of all-atom frames with any oxygen within 2.6 Å (at most 5 % asked)',
        check=D17M['X0'][t_]['all_atom']['any']['pct_within_2p6'])
    add(f'md.gate.GV.{t_}.so_zn_mean', SX0_[t_]['stats']['d_Zn_Osulfonyl_min_A']['mean'], 'A', 3, f'md/results/runs/bonded_{t_}/X0/first_program_summary.json:stats.d_Zn_Osulfonyl_min_A.mean', f'validation run X0, zinc set {t_}: mean distance of the nearer sulfonyl oxygen from the zinc (2.986 Å in 6RL9)')
    c4_ = g_['GV4_ligand_internal']['criteria']
    add(f'md.gate.GV.{t_}.c4sn2', c4_['ang_C4_S_N2_deg']['mean_deg'], 'deg', 2, f'{src_}.GV4_ligand_internal.criteria.ang_C4_S_N2_deg.mean_deg', f'gate G-V, zinc set {t_}: mean C4–S–N2 (108.8° in 6RL9; within 3° asked)')
    add(f'md.gate.GV.{t_}.tau_fold', c4_['tau_fold_deg']['mean_deg'], 'deg', 2, f'{src_}.GV4_ligand_internal.criteria.tau_fold_deg.mean_deg', f'gate G-V, zinc set {t_}: mean τ_fold (132.1° in 6RL9; within 20° asked)')
    add(f'md.gate.GV.{t_}.outcome', g_['verdict'].lower(), None, None, f'{src_}.verdict', f'gate G-V, zinc set {t_}: outcome' + (' (one part failed: an oxygen within 2.6 Å in 3 of 50 all-atom frames, 2 allowed)' if t_ == '2s' else ''))
_kg = _k['gate_GV']
_fl0 = SX0_['2e']['gate_GV']['GV3_zinc']['criteria']['fifth_ligand']
add('md.gate.GV.n_solute', _fl0['solute_trajectory']['n'], None, 0, 'md/results/runs/bonded_2e/X0/first_program_summary.json:gate_GV.GV3_zinc.criteria.fifth_ligand.solute_trajectory.n', 'gate G-V: solute frames of a validation run')
add('md.gate.GV.n_full', _fl0['full_trajectory']['n'], None, 0, 'md/results/runs/bonded_2e/X0/first_program_summary.json:gate_GV.GV3_zinc.criteria.fifth_ligand.full_trajectory.n', 'gate G-V: all-atom frames of a validation run',
    check=D17M['X0']['2e']['all_atom']['n'], tol=0.01)
add('md.gate.GV.n_full_allowed', int(_fl0['limit_max_fraction'] * _fl0['full_trajectory']['n'] + 1e-9), None, 0, 'md/results/runs/bonded_2e/X0/first_program_summary.json:gate_GV.GV3_zinc.criteria.fifth_ligand (limit_max_fraction × n)',
    'gate G-V: all-atom frames with an oxygen within 2.6 Å of the zinc that the criterion allows')
add('md.gate.GV.limit.rmsd', _kg['lig_rmsd_mean_max_A'], 'A', 1, f'{_KS}.gate_GV.lig_rmsd_mean_max_A', 'gate G-V, criterion: largest mean ligand RMSD')
add('md.gate.GV.limit.hb_pct', 100.0 * _kg['hb_N2H_T199OG1_fraction_min'][0] / _kg['hb_N2H_T199OG1_fraction_min'][1], '%', 0, f'{_KS}.gate_GV.hb_N2H_T199OG1_fraction_min', 'gate G-V, criterion: smallest share of frames with the hydrogen bond to Thr199')
add('md.gate.GV.limit.c4sn2', _kg['ang_C4_S_N2_tolerance_deg'], 'deg', 0, f'{_KS}.gate_GV.ang_C4_S_N2_tolerance_deg', 'gate G-V, criterion: mean C4–S–N2 within this of 108.8°')
add('md.gate.GV.limit.tau_fold', _kg['tau_fold_tolerance_deg'], 'deg', 0, f'{_KS}.gate_GV.tau_fold_tolerance_deg', 'gate G-V, criterion: mean τ_fold within this of 132.1°')
add('md.xtal.c4sn2', _kg['ang_C4_S_N2_reference_deg'], 'deg', 1, f'{_KS}.gate_GV.ang_C4_S_N2_reference_deg', '6RL9: angle C4–S–N2 of the crystal ligand')
add('md.xtal.tau_fold', _kg['tau_fold_reference_deg'], 'deg', 1, f'{_KS}.gate_GV.tau_fold_reference_deg', '6RL9: τ_fold of the crystal ligand')
add('md.xtal.so_zn', DEC_['starts']['X']['Zn_O_min'], 'A', 3, 'md/system/decisions.json:starts.X.Zn_O_min', '6RL9: distance of the nearer sulfonyl oxygen from the zinc', check=_XT['d_Zn_Osulfonyl_min_A'], tol=0.0006)
add('md.xtal.s_c4', DEC_['starts']['X']['C4_S'], 'A', 3, 'md/system/decisions.json:starts.X.C4_S', '6RL9: bond S–C4 of the crystal ligand', check=_XT['d_S_C4_A'], tol=0.0006)
add('md.xtal.o1_s_c4', _XT['ang_O1_S_C4_deg'], 'deg', 1, f'{D17S}:md.xtal_6RL9.ang_O1_S_C4_deg (md/system/ligand_X.pdb)', '6RL9: angle O1–S–C4 of the crystal ligand')
add('md.xtal.zn_oop_his119', _XT['zn_oop_his119'], 'A', 2, f'{D17S}:md.xtal_6RL9.zn_oop_his119 (md/system/protein_md.pdb, md/system/zinc.pdb)',
    '6RL9: zinc out of the least-squares plane of the His119 ring (CG, ND1, CD2, CE1, NE2)',
    note='the plan (7 Oct 2026, 12:10 AM IST, point 8) prints 0.17 Å here, and md/RESULTS.md repeats it; that is the value of His94 (−0.17 Å; His96 0.06 Å, His119 0.15 Å), as the physical review of the '
         'campaign has it (md/results/independent/physical_review/REVIEW.md)')
assert abs(_XT['zinc_out_of_the_ring_plane_A']['His94'] + 0.17) < 0.005 and abs(_XT['zinc_out_of_the_ring_plane_A']['His96'] - 0.06) < 0.005 and abs(_XT['zn_oop_his119'] - 0.15) < 0.005
add('md.xtal.n2_og1', DEC_['starts']['X']['N2_OG1'], 'A', 2, 'md/system/decisions.json:starts.X.N2_OG1', '6RL9: distance N2···OG1(Thr199)')
add('md.xtal.n2_h_og1', DEC_['starts']['X']['N2_H_OG1'], 'deg', 1, 'md/system/decisions.json:starts.X.N2_H_OG1', '6RL9 with the constructed hydrogen: angle N2–H···OG1(Thr199)')
add('md.xtal.s_n2_zn', DEC_['starts']['X']['S_N2_Zn'], 'deg', 1, 'md/system/decisions.json:starts.X.S_N2_Zn', '6RL9: angle S–N2–Zn', check=_XT['ang_S_N2_Zn_deg'], tol=0.006)
for a_, v_ in _cz['c_references_6RL9_deg'].items():
    add(f'md.xtal.{a_[:-4]}', v_, 'deg', 1, f'md/production/gate_GZ_bonded_2e.json:criteria.c_references_6RL9_deg.{a_}', f'6RL9: angle {a_[4:-4].replace("_", "–")} at the zinc')

# G-R: the terms at the sulfonamide nitrogen against quantum chemistry on the small zinc-site model (failed)
_va, _gd = GR_['verdict_part_A'], D17M['gate_GR']
_GS = 'md/qm/gr/gr_gate.json'
add('md.gate.GR.rms', _va['rms'], 'kcal mol-1', 2, f'{_GS}:verdict_part_A.rms', 'gate G-R, part A: root mean square of force field minus quantum chemistry (B3LYP-D3(BJ)/6-31G*) over 23 comparisons on the zinc-site model, each family against its lowest point',
    check=_gd['with_d3']['rms'])
add('md.gate.GR.rms_no_d3', _va['rms_without_d3'], 'kcal mol-1', 2, f'{_GS}:verdict_part_A.rms_without_d3', 'gate G-R, part A: the same without the D3 term of the reference', check=_gd['without_d3']['rms'])
add('md.gate.GR.limit', _va['criterion'], 'kcal mol-1', 1, f'{_GS}:verdict_part_A.criterion', 'gate G-R, criterion: largest root mean square allowed')
add('md.gate.GR.n', _va['n_compared'], None, 0, f'{_GS}:verdict_part_A.n_compared', 'gate G-R, part A: comparisons', check=_gd['with_d3']['n_compared'], tol=0.01)
add('md.gate.GR.largest', _va['largest'], 'kcal mol-1', 2, f'{_GS}:verdict_part_A.largest', f'gate G-R, part A: largest single deviation (point {_va["largest_at"]}, the planar nitrogen of the constructed path)', check=_gd['with_d3']['largest'])
add('md.gate.GR.outcome', 'passed' if _va['passed'] else 'failed', None, None, f'{_GS}:verdict_part_A.passed', 'gate G-R: outcome')
_pa = GR_['part_A_composition']['points']
for k_, pt_, m_ in (('planar_N', 'a4', 'planar nitrogen on the constructed inversion path'), ('H_in_plane', 'd1', 'only the hydrogen moved into the plane of its neighbours'),
                    ('turn_30', 'e_p030', 'ligand turned rigidly by 30° about Zn–N2')):
    add(f'md.gate.GR.{k_}.qm', _pa[pt_]['qm'], 'kcal mol-1', 2, f'{_GS}:part_A_composition.points.{pt_}.qm', f'gate G-R, {m_}: quantum-chemical energy relative to the crystal-like placement X', check=_gd['points'][pt_]['qm'])
    add(f'md.gate.GR.{k_}.mm', _pa[pt_]['mm'], 'kcal mol-1', 2, f'{_GS}:part_A_composition.points.{pt_}.mm', f'gate G-R, {m_}: force-field energy relative to X', check=_gd['points'][pt_]['mm'], tol=0.006)
_xb, _xn = GR_['with_d3']['part_B']['X_to_Rx'], GR_['without_d3']['part_B']['X_to_Rx']
add('md.gate.GR.Rx.qm', _xb['qm'], 'kcal mol-1', 2, f'{_GS}:with_d3.part_B.X_to_Rx.qm', 'gate G-R, part B (reported, no criterion): the unrelaxed R-like placement of the ligand against X in quantum chemistry', check=_gd['points']['Rx']['qm'])
add('md.gate.GR.Rx.mm', _xb['mm'], 'kcal mol-1', 2, f'{_GS}:with_d3.part_B.X_to_Rx.mm', 'gate G-R, part B: the same in the force field as built', check=_gd['points']['Rx']['mm'], tol=0.006)
add('md.gate.GR.Rx.mm_minus_qm', _xb['mm_minus_qm'], 'kcal mol-1', 2, f'{_GS}:with_d3.part_B.X_to_Rx.mm_minus_qm', 'gate G-R, part B: force field minus quantum chemistry for the R-like placement (too high in the force field; the error with the sign of the outcome)',
    check=_gd['points']['Rx']['mm_minus_qm'])
add('md.gate.GR.Rx.mm_minus_qm_no_d3', _xn['mm_minus_qm'], 'kcal mol-1', 2, f'{_GS}:without_d3.part_B.X_to_Rx.mm_minus_qm', 'gate G-R, part B: the same without the D3 term of the reference (the sign rests on the dispersion term)',
    check=_gd['points']['Rx']['mm_minus_qm_without_d3'])
add('md.gate.GR.Rx.qm_no_d3', _xn['qm'], 'kcal mol-1', 2, f'{_GS}:without_d3.part_B.X_to_Rx.qm', 'gate G-R, part B: the R-like placement against X in quantum chemistry without the D3 term', check=_gd['points']['Rx']['qm_without_d3'])
_c2 = N2G_['candidates']['second_step']['for_the_record']['Rx_against_X']
add('md.gate.GR.Rx.mm_two_placements', _c2['as_built'], 'kcal mol-1', 3, 'md/qm/n2/n2_gate_result.json:candidates.second_step.for_the_record.Rx_against_X.as_built',
    'the R-like placement against X with the terms as built and charges fitted to both placements (the charges of the runs are fitted to pose X alone); the record gives three decimals, a tie at two, so the mismatch 0.85 is the number to print')
add('md.gate.GR.Rx.mismatch_two_placements', _c2['mismatch_before'], 'kcal mol-1', 2, 'md/qm/n2/n2_gate_result.json:candidates.second_step.for_the_record.Rx_against_X.mismatch_before',
    'force field minus quantum chemistry for the R-like placement with charges fitted to both placements')
# G-R2: no better model of the nitrogen was accepted
_cb, _cs = N2G_['candidates']['charges_as_built'], N2G_['candidates']['second_step']
_NS = 'md/qm/n2/n2_gate_result.json'
_K5, _K6 = '5 leave-one-out error on the fitting set', '6 force-field minimum of the small model against its quantum-chemical optimum'
add('md.gate.GR2.n_candidates', sum(v_['n'] for v_ in N2G_['every_candidate_of_the_ladder'].values()), None, 0, f'{_NS}:every_candidate_of_the_ladder.*.n (sum)', 'gate G-R2: candidate models of the terms at N2 (five models, three ridge values, two sets of charges)')
add('md.gate.GR2.loo.as_built_charges', N2G_['every_candidate_of_the_ladder']['charges_as_built']['smallest_leave_one_out_rms'], 'kcal mol-1', 2, f'{_NS}:every_candidate_of_the_ladder.charges_as_built.smallest_leave_one_out_rms',
    'gate G-R2: smallest leave-one-out error of a candidate with the charges as built')
add('md.gate.GR2.loo.two_placements', N2G_['every_candidate_of_the_ladder']['second_step']['smallest_leave_one_out_rms'], 'kcal mol-1', 2, f'{_NS}:every_candidate_of_the_ladder.second_step.smallest_leave_one_out_rms',
    'gate G-R2: smallest leave-one-out error of a candidate with charges fitted to two placements')
add('md.gate.GR2.limit.loo', N2G_['limits']['loo_rms'], 'kcal mol-1', 1, f'{_NS}:limits.loo_rms', 'gate G-R2, criterion: largest leave-one-out error allowed')
add('md.gate.GR2.angle_miss', abs(_cb['conditions'][_K6]['largest_abs_difference_deg']), 'deg', 1, f'{_NS}:candidates.charges_as_built.conditions["{_K6}"].largest_abs_difference_deg',
    'gate G-R2: Zn–N2–S of the force-field minimum of the small model against the quantum-chemical optimum, candidate with the charges as built (too small by this)')
add('md.gate.GR2.angle_miss_two_placements', abs(_cs['conditions'][_K6]['largest_abs_difference_deg']), 'deg', 1, f'{_NS}:candidates.second_step.conditions["{_K6}"].largest_abs_difference_deg', 'gate G-R2: the same for the candidate with charges fitted to two placements')
add('md.gate.GR2.limit.angle', N2G_['limits']['minimum_deg'], 'deg', 0, f'{_NS}:limits.minimum_deg', 'gate G-R2, criterion: the angles at N2 of the force-field minimum within this of the quantum-chemical optimum')
add('md.gate.GR2.n_fit', N2G_['integrity']['fit'], None, 0, f'{_NS}:integrity.fit', 'gate G-R2: new quantum-chemical single points of the fitting set (the fits also use 36 points of gate G-R)')
assert _cb['conditions'][_K5]['n_fitting_rows'] == _cs['conditions'][_K5]['n_fitting_rows']
add('md.gate.GR2.n_fit_rows', _cb['conditions'][_K5]['n_fitting_rows'], None, 0, f'{_NS}:candidates.charges_as_built.conditions["{_K5}"].n_fitting_rows', 'gate G-R2: single-point energies that every candidate is fitted to (36 points of gate G-R and 141 new ones)')
add('md.gate.GR2.n_fit_old', _cb['conditions'][_K5]['n_fitting_rows'] - N2G_['integrity']['fit'], None, 0, f'{_NS}:conditions n_fitting_rows − integrity.fit', 'gate G-R2: points of gate G-R among the fitting rows')
add('md.gate.GR2.n_parameters', _cb['conditions'][_K5]['n_parameters'], None, 0, f'{_NS}:candidates.charges_as_built.conditions["{_K5}"].n_parameters', 'gate G-R2: free parameters of the largest model of the ladder (the candidate put to the test)')
add('md.gate.GR2.n_test', N2G_['integrity']['test'], None, 0, f'{_NS}:integrity.test', 'gate G-R2: single points of the test set; their results were never opened, because the verdict was fixed without them')
add('md.gate.GR2.outcome', N2G_['verdict'].lower(), None, None, f'{_NS}:verdict', 'gate G-R2: outcome (the force field stayed as built)')
_fr = _cb['for_the_record']
add('md.gate.GR2.cand.Rx', _fr['Rx_against_X']['new'], 'kcal mol-1', 2, f'{_NS}:candidates.charges_as_built.for_the_record.Rx_against_X.new', 'candidate terms at N2 (run R1_n2fit): the R-like placement against X (14.22 as built, 10.84 in quantum chemistry)')
add('md.gate.GR2.cand.H_in_plane', _fr['hydrogen_moved_into_the_plane']['new'], 'kcal mol-1', 2, f'{_NS}:candidates.charges_as_built.for_the_record.hydrogen_moved_into_the_plane.new', 'candidate terms at N2: the hydrogen moved into the plane (0.46 in quantum chemistry, 5.20 as built)')
_fl = N2G_['for_the_record_no_weight_in_the_verdict']['relaxed_geometries']['flattest_nitrogen_of_the_quantum_chemical_run_from_Rx']['force_field_minus_QM']
_FS = f'{_NS}:for_the_record_no_weight_in_the_verdict.relaxed_geometries.flattest_nitrogen_of_the_quantum_chemical_run_from_Rx.force_field_minus_QM'
add('md.gate.GR2.flat_N.as_built', _fl['T0'], 'kcal mol-1', 1, f'{_FS}.T0', 'at the relaxed near-planar nitrogen of the quantum-chemical relaxation: terms as built minus quantum chemistry')
add('md.gate.GR2.flat_N.candidate', _fl['T1'], 'kcal mol-1', 1, f'{_FS}.T1', 'at the relaxed near-planar nitrogen: candidate terms (charges as built) minus quantum chemistry')
add('md.gate.GR2.flat_N.candidate_two_placements', _fl['T2'], 'kcal mol-1', 1, f'{_FS}.T2', 'at the relaxed near-planar nitrogen: candidate with charges of two placements (not run) minus quantum chemistry')
# the quantum-chemical relaxation of the bare zinc-site model from the R-like placement (converged at step 66)
_rx = BEH_['along_the_quantum_chemical_relaxations']['runs']['Rx']
_BS = 'md/qm/n2/final_test/behaviour.json:along_the_quantum_chemical_relaxations.runs.Rx'
_gl = N2G_['for_the_record_no_weight_in_the_verdict']['relaxed_geometries']['quantum_chemical_run_from_Rx_last_step']
assert _gl['step'] == _rx['last_step'] and _gl['rmsd_ligand_heavy_from_X_run'] == _rx['last']['rmsd_ligand_heavy_from_X_run'] and _gl['e_total_minus_X_run_kcal_mol'] == _rx['last']['e_total_minus_X_run_kcal_mol']
add('md.qm.relax.rmsd', _rx['last']['rmsd_ligand_heavy_from_X_run'], 'A', 2, f'{_BS}.last.rmsd_ligand_heavy_from_X_run',
    'relaxation of the ligand in the bare zinc-site model (64 atoms, vacuum, zinc and imidazoles fixed, B3LYP-D3(BJ)/6-31G*) started from the R-like placement: heavy-atom RMSD of its end point from the end point of the relaxation started from X',
    note=f'the same value in {_NS}:for_the_record_no_weight_in_the_verdict.relaxed_geometries.quantum_chemical_run_from_Rx_last_step; md/qm/gr/relax/qm_Rx.json in git is the record at step 23, before the run converged; the converged record of the machine is md/qm/gr/relax/qm_Rx_converged_step66.json (state done, 66 steps), with its geometry qm_Rx_converged_step66.xyz and the later summaries relax_qm_both_converged.json, gr_gate_after_Rx_converged.json and binding_check_final.json beside it (byte copies added to git on 11 Oct 2026; md/qm/gr/relax/README.md says which record is which)')
add('md.qm.relax.de', _rx['last']['e_total_minus_X_run_kcal_mol'], 'kcal mol-1', 2, f'{_BS}.last.e_total_minus_X_run_kcal_mol', 'the same relaxation: energy of its end point minus that of the relaxation started from X')
add('md.qm.relax.n_steps', _rx['n_steps'], None, 0, f'{_BS}.n_steps', 'the same relaxation: steps to convergence')
_fz = J(MDR_ / 'qm/gr/relax/qm_Rx.json')['frozen_atoms'].split('-')
add('md.qm.relax.n_fixed', int(_fz[1]) - int(_fz[0]) + 1, None, 0, 'md/qm/gr/relax/qm_Rx.json:frozen_atoms', 'the same relaxation: atoms held fixed (the zinc and the three ethylimidazoles; the 18 atoms of the ligand move)')
assert int(_fz[1]) - int(_fz[0]) + 1 + _bx['composition']['ligand_atoms'] == STEP1_['small']['n_atoms']
add('md.qm.relax.volume_start', RQM_['starts']['Rx']['signed_volume_N2']['start'], None, 2, 'md/qm/gr/relax_qm.json:starts.Rx.signed_volume_N2.start', 'the same relaxation: signed volume at N2 at the start (negative: the configuration of the R-like placement)')
add('md.qm.relax.volume_end', _rx['last']['signed_volume'], None, 2, f'{_BS}.last.signed_volume', 'the same relaxation: signed volume at N2 at the end (positive: the nitrogen has inverted to the configuration of the crystal placement)')
_cvf = MDR_ / 'qm/gr/relax/qm_Rx_converged_step66.json'            # the converged record of the machine, brought into git on 11 Oct 2026
if _cvf.exists():
    _cv = J(_cvf)
    assert _cv['state'] == 'done' and _cv['crossed'] and _cv['n_steps'] == _rx['n_steps'] and _cv['final']['rmsd_ligand_heavy_from_X_run'] == _rx['last']['rmsd_ligand_heavy_from_X_run'] \
        and _cv['final']['e_total_minus_X_run_kcal_mol'] == _rx['last']['e_total_minus_X_run_kcal_mol'] and _cv['final']['signed_volume'] == _rx['last']['signed_volume'] \
        and _cv['signed_volume_start'] == RQM_['starts']['Rx']['signed_volume_N2']['start'] and _cv['final']['X_run_at'] == dict(step=RQM_['starts']['X']['n_steps'], converged=True)
add('md.qm.relax.X.n_steps', RQM_['starts']['X']['n_steps'], None, 0, 'md/qm/gr/relax_qm.json:starts.X.n_steps', 'relaxation started from the crystal-like placement X: steps to convergence')
add('md.qm.relax.X.rmsd_from_start', RQM_['starts']['X']['final']['rmsd_ligand_heavy_from_start'], 'A', 2, 'md/qm/gr/relax_qm.json:starts.X.final.rmsd_ligand_heavy_from_start',
    'relaxation started from X: the bare model moves the crystal-like ligand this far (heavy-atom RMSD) from its crystal placement')

# ================================================================================== agreement of the readings
_RPS = 'md/results/campaign/second_program/recheck_report.json'
add('md.agree.n_values', sum(x_['comparison']['n_values'] for x_ in REP_['runs']), None, 0, f'{_RPS}:runs[].comparison.n_values (sum)', 'values of the run folders compared between the two analysis programs')
add('md.agree.n_folders', REP_['n_runs'], None, 0, f'{_RPS}:n_runs', 'run folders compared (the 13 runs and two pre-flight test folders)')
add('md.agree.n_campaign', REP_['campaign']['n_compared'], None, 0, f'{_RPS}:campaign.n_compared', 'values of the campaign rules compared')
add('md.agree.n_disagree', sum(x_['comparison']['n_disagree'] for x_ in REP_['runs']) + REP_['campaign']['n_disagree'], None, 0, f'{_RPS}:runs[].comparison.n_disagree, campaign.n_disagree', 'disagreements beyond tolerance between the two programs')
add('md.agree.n_threshold', REP_['n_threshold_items'], None, 0, f'{_RPS}:n_threshold_items', 'frames that the two programs class differently at a threshold')
add('md.agree.max_diff_series', REP_['largest_differences']['series dist'], 'A', 4, f'{_RPS}:largest_differences["series dist"]', 'largest difference of a distance in a time series between the two programs (one unit of the last decimal written)')
add('md.agree.third.n_values', _c17['third_reading']['n_values'], None, 0, f'{D17S}:md.comparison.third_reading.n_values',
    'values of the blind third reading (own trajectory reader) compared by derive/d17_md.py with the first program and with its own re-derivation; the assembler of md/RESULTS.md compared 176 by a list of its own')
add('md.agree.third.n_different', len(_c17['mismatches']), None, 0, f'{D17S}:md.comparison.mismatches', 'of these, counts, verdicts or times that differ')
add('md.agree.third.max_diff', _c17['third_reading']['largest_difference_of_a_mean'], 'A', 4, f'{D17S}:md.comparison.third_reading.largest_difference_of_a_mean', 'largest difference of a mean between the first program and the third reading (which writes four decimals)')
add('md.agree.raw.n_solute', sum(x_['n_solute_frames'] for x_ in SS_.values()), None, 0, 'md/results/site_checks_solute.json:runs[].n_solute_frames (sum)', 'solute frames of the 13 runs, equilibration included, read once more from the raw trajectories (no frame classed differently)')
add('md.agree.raw.n_full', sum(x_['n_all_atom_frames'] for x_ in SF_.values()), None, 0, 'md/results/site_checks_full.json:runs[].n_all_atom_frames (sum)', 'all-atom frames of the 13 runs read once more from the raw trajectories')
_fv = HEALTH_['final_verification']
add('md.agree.files_checked', _fv['driver_checksums']['files_ok'], None, 0, 'md/results/run_health.json:final_verification.driver_checksums.files_ok', 'trajectory and record files whose checksum equals the record of the driver (none bad, none missing)')
add('md.agree.files_two_machines', _fv['pi_against_machine']['same'], None, 0, 'md/results/run_health.json:final_verification.pi_against_machine.same', f'files with the same checksum on both machines (of {_fv["pi_against_machine"]["files"]:,})')

# ===================================================================================================== Glu106
LIGN_ = {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate', 'azm_neutral': 'acetazolamide', 'azm_anion': 'acetazolamide anion'}
SCN_ = {'ad4zn': 'AutoDock4Zn', 'vina': 'Vina', 'vinardo': 'Vinardo'}
STN_ = {'acid': 'Glu106 acid (receptor of the paper)', 'carb': 'Glu106 carboxylate'}


def _runkey(run_):
    p_ = run_.split('__')
    return p_[0], p_[1], ('prod' if len(p_) > 2 and p_[2] == 'prod' else 'val')


def _runlab(pdb_, run_, n_):
    lig_, sc_, kind_ = _runkey(run_)
    geo_ = ('DFT geometry' if lig_.startswith('sfa_') else 'ETKDG conformers') if kind_ == 'prod' else 'ETKDG conformers'
    return f'{pdb_}, {LIGN_[lig_]}, {SCN_[sc_]}, {geo_} ({"production" if kind_ == "prod" else "validation"} run, {n_} seeds)'


def _dockfacts(b_, rec_, own_, src_, lab_, sfa6_):
    """The facts of one docking run in one receptor: rec_ is its record in comparison(_water).json, own_ the parse
    of d17."""
    n_ = rec_['n_seeds']
    add(b_ + '.n_seeds', n_, None, 0, f'{src_}.n_seeds', f'{lab_}: seeds', check=own_['n_seeds'], tol=0.01)
    _a(b_ + '.rmsd', _dmean(rec_['rmsd']['per_seed']), 'A', 2, f'{src_}.rmsd.per_seed (mean)', f'{lab_}: RMSD of the top pose from the crystal ligand, mean over the seeds', check=own_['rmsd']['mean'], tol=0.0006, full=True)
    _a(b_ + '.rmsd_sd', _stt.stdev(rec_['rmsd']['per_seed']), 'A', 2, f'{src_}.rmsd.per_seed (SD)', f'{lab_}: RMSD of the top pose, SD over the seeds', check=own_['rmsd']['sd'], tol=0.0006, full=True)
    add(b_ + '.n_le_1', round(rec_['rmsd']['f_le_1A'] * n_), None, 0, f'{src_}.rmsd.f_le_1A × n_seeds', f'{lab_}: seeds whose top pose lies within 1.0 Å of the crystal ligand', check=own_['rmsd']['n_le_1'], tol=0.01)
    add(b_ + '.f1', rec_['rmsd']['f_le_1A'], None, 2, f'{src_}.rmsd.f_le_1A', f'{lab_}: fraction of seeds with the top pose within 1.0 Å', check=own_['rmsd']['f_le_1'], tol=0.001)
    add(b_ + '.f2', rec_['rmsd']['f_le_2A'], None, 2, f'{src_}.rmsd.f_le_2A', f'{lab_}: fraction of seeds with the top pose within 2.0 Å', check=own_['rmsd']['f_le_2'], tol=0.001)
    _a(b_ + '.score', _dmean(rec_['score']['per_seed']), 'kcal mol-1', 3, f'{src_}.score.per_seed (mean)', f'{lab_}: score of the top pose, mean over the seeds (an empirical docking score)', check=own_['score']['mean'], tol=1e-6, exact=True)
    _a(b_ + '.score_sd', _stt.stdev(rec_['score']['per_seed']), 'kcal mol-1', 3, f'{src_}.score.per_seed (SD)', f'{lab_}: score of the top pose, SD over the seeds', check=own_['score']['sd'], tol=1e-6)
    zs_ = [x_['Zn-N(sulfonamide)'] for x_ in rec_['site_checks_per_seed']]
    os_ = [x_['Thr199 OG1 H-bond'] for x_ in rec_['site_checks_per_seed']]
    _a(b_ + '.zn_n', _dmean(zs_), 'A', 2, f'{src_}.site_checks_per_seed[]["Zn-N(sulfonamide)"] (mean)', f'{lab_}: distance from the zinc to the sulfonamide nitrogen of the top pose, mean over the seeds', check=own_['zn_n']['mean'], tol=0.006, full=True)
    add(b_ + '.zn_n_pass', round(rec_['site_checks']['Zn-N(sulfonamide)']['pass_fraction'] * n_), None, 0, f'{src_}.site_checks["Zn-N(sulfonamide)"].pass_fraction × n_seeds', f'{lab_}: seeds whose top pose has Zn–N within 2.3 Å',
        check=own_['zn_n']['n_within_2p3'], tol=0.01)
    _a(b_ + '.og1_n', _dmean(os_), 'A', 2, f'{src_}.site_checks_per_seed[]["Thr199 OG1 H-bond"] (mean)', f'{lab_}: distance from Thr199 OG1 to the sulfonamide nitrogen of the top pose, mean over the seeds', check=own_['og1_n']['mean'], tol=0.006, full=True)
    add(b_ + '.og1_n_pass', round(rec_['site_checks']['Thr199 OG1 H-bond']['pass_fraction'] * n_), None, 0, f'{src_}.site_checks["Thr199 OG1 H-bond"].pass_fraction × n_seeds', f'{lab_}: seeds whose top pose has Thr199 OG1···N within 3.3 Å',
        check=own_['og1_n']['n_within_3p3'], tol=0.01)
    o_ = rec_['orientation']
    add(b_ + '.top_near', o_['top_near_native'], None, 0, f'{src_}.orientation.top_near_native', f'{lab_}: seeds whose top pose is near-native (within 1.0 Å)', check=own_['top']['near_native'], tol=0.01)
    add(b_ + '.top_beyond2', o_['top_ring_rotated'], None, 0, f'{src_}.orientation.top_ring_rotated', f'{lab_}: seeds whose top pose lies beyond 2.0 Å from the crystal ligand (the ring-rotated orientation of the paper)', check=own_['top']['beyond_2'], tol=0.01)
    add(b_ + '.top_between', o_['top_between'], None, 0, f'{src_}.orientation.top_between', f'{lab_}: seeds whose top pose lies between 1.0 and 2.0 Å', check=own_['top']['between'], tol=0.01)
    if '.ad4zn.' in b_:
        add(b_ + '.seeds_near', o_['seeds_with_near_native_mode'], None, 0, f'{src_}.orientation.seeds_with_near_native_mode', f'{lab_}: seeds that hold a near-native mode', check=own_['modes']['seeds_holding_near_native'], tol=0.01)
        if o_.get('worst_rank_near') is not None:
            add(b_ + '.worst_rank_near', o_['worst_rank_near'], None, 0, f'{src_}.orientation.worst_rank_near', f'{lab_}: worst rank of the first near-native mode over the seeds', check=own_['orientation']['worst_rank_near'], tol=0.01)
        if o_.get('sep_min') is not None and o_['seeds_with_both'] == n_:
            add(b_ + '.sep_min', o_['sep_min'], 'kcal mol-1', 3, f'{src_}.orientation.sep_min', f'{lab_}: score separation of the best near-native mode and the best mode beyond 2.0 Å, smallest over the seeds', check=own_['orientation']['sep'][0], tol=1e-6)
            add(b_ + '.sep_max', o_['sep_max'], 'kcal mol-1', 3, f'{src_}.orientation.sep_max', f'{lab_}: the same separation, largest over the seeds', check=own_['orientation']['sep'][1], tol=1e-6)
            add(b_ + '.near_minus_far_mean', o_['near_minus_rotated_mean'], 'kcal mol-1', 3, f'{src_}.orientation.near_minus_rotated_mean', f'{lab_}: best near-native minus best mode beyond 2.0 Å, mean over the seeds (negative: near-native scores better)',
                check=own_['orientation']['near_minus_far_mean'], tol=0.0006)
        if o_.get('gap_max') is not None:
            add(b_ + '.gap_max', o_['gap_max'], 'kcal mol-1', 3, f'{src_}.orientation.gap_max', f'{lab_}: first near-native mode behind the top mode, largest over the seeds', check=own_['orientation']['gap'][1], tol=1e-6)
    g_ = rec_['gate']
    if g_.get('passes_G2') is not None:
        add(b_ + '.gate', 'pass' if g_['passes_G2'] else 'fail', None, None, f'{src_}.gate.passes_G2', f'{lab_}: the strict gate (mean RMSD at most 1.0 Å, half of the seeds within 1.0 Å, PoseBusters physics, site checks)')
        add(b_ + '.crit2A', 'pass' if g_['passes_2A_criterion'] else 'fail', None, None, f'{src_}.gate.passes_2A_criterion', f'{lab_}: the conventional 2.0 Å criterion')
        add(b_ + '.pb', rec_['posebusters_physics_pass_fraction'], None, 1, f'{src_}.posebusters_physics_pass_fraction', f'{lab_}: fraction of top poses that pass the PoseBusters physics checks')
        assert (not g_['passes_G2']) or own_['gate_part']['strict'], b_
        assert (not g_['passes_2A_criterion']) or own_['gate_part']['two_A'], b_
    else:
        add(b_ + '.gate_part', 'met' if g_['rmsd_and_site_part_G2'] else 'not met', None, None, f'{src_}.gate.rmsd_and_site_part_G2',
            f'{lab_}: the RMSD and site parts of the strict gate (no PoseBusters run on production runs, so no verdict)', note=f'the same from the raw poses ({D17S})')
        assert g_['rmsd_and_site_part_G2'] == own_['gate_part']['strict'], b_
    if sfa6_ and rec_.get('tau_fold_deg'):
        add(b_ + '.tau_fold', rec_['tau_fold_deg']['mean'], 'deg', 1, f'{src_}.tau_fold_deg.mean', f'{lab_}: τ_fold of the top poses, mean over the seeds')


GLU_ = {}
for base_, var_ in (('T3_6RL9', 'T3_6RL9_glu106GLU'), ('T4_3HS4', 'T4_3HS4_glu106GLU')):
    pdb_ = base_.split('_')[1]
    comp_ = J(DOCK / var_ / 'comparison.json')
    for run_, rec_ in comp_['runs'].items():
        lig_, sc_, kind_ = _runkey(run_)
        GLU_[(pdb_, run_)] = rec_
        for side_, st_ in (('baseline', 'acid'), ('variant', 'carb')):
            _dockfacts(f'dock.glu106.{pdb_}.{lig_}.{sc_}.{kind_}.{st_}', rec_[side_], D17D['glu106'][f'{pdb_}/{run_}']['acid' if st_ == 'acid' else 'carboxylate'],
                       f'dock/{var_}/comparison.json:runs.{run_}.{side_}', f'{_runlab(pdb_, run_, rec_[side_]["n_seeds"])}, {STN_[st_]}', lig_.startswith('sfa_') and pdb_ == '6RL9' or lig_.startswith('sfa_'))
        k_ = f'dock.glu106.{pdb_}.{lig_}.{sc_}.{kind_}'
        d17_ = D17D['glu106'][f'{pdb_}/{run_}']
        add(k_ + '.delta_score', FACTS[k_ + '.carb.score']['value'] - FACTS[k_ + '.acid.score']['value'], 'kcal mol-1', 3, f'dock/{var_}/comparison.json:runs.{run_} (carboxylate minus acid, mean scores)',
            f'{_runlab(pdb_, run_, rec_["baseline"]["n_seeds"])}: mean score with the carboxylate minus that with the acid', check=d17_['carboxylate']['score']['mean'] - d17_['acid']['score']['mean'], tol=1e-6)
_gv = D17D['glu106_vina']
_vr = [(k_, v_) for k_, v_ in GLU_.items() if '__vin' in k_[1]]
add('dock.glu106.n_runs_vina', len(_vr), None, 0, 'dock/T3_6RL9_glu106GLU/comparison.json, dock/T4_3HS4_glu106GLU/comparison.json:runs (Vina and Vinardo)', 'Vina and Vinardo runs repeated with Glu106 as a carboxylate (10 seeds each)', check=_gv['n_runs'], tol=0.01)
add('dock.glu106.n_runs_ad4zn', len(GLU_) - len(_vr), None, 0, 'dock/T3_6RL9_glu106GLU/comparison.json, dock/T4_3HS4_glu106GLU/comparison.json:runs (AutoDock4Zn)', 'AutoDock4Zn runs repeated with Glu106 as a carboxylate (the six input-geometry diagnostics of the paper were not repeated)',
    check=_gv['n_runs_ad4zn'], tol=0.01)
add('dock.glu106.vina.n_verdicts_changed', sum(1 for k_, v_ in _vr if v_['baseline']['gate'].get('passes_G2') != v_['variant']['gate'].get('passes_G2') or v_['baseline']['gate'].get('passes_2A_criterion') != v_['variant']['gate'].get('passes_2A_criterion')),
    None, 0, 'comparison.json of both variants: runs.*.baseline.gate against runs.*.variant.gate (Vina and Vinardo)', 'Vina and Vinardo runs whose verdict on the strict gate or on the 2.0 Å criterion changes with the carboxylate')
add('dock.glu106.vina.max_score_change_sfa', max(abs(FACTS[f'dock.glu106.{p_}.{_runkey(r_)[0]}.{_runkey(r_)[1]}.{_runkey(r_)[2]}.delta_score']['value']) for (p_, r_), v_ in _vr if r_.startswith('sfa_')), 'kcal mol-1', 3,
    'comparison.json of both variants (largest |carboxylate − acid| of a mean score, sulfanilamide and sulfonamidate, Vina and Vinardo)', 'largest change of a mean Vina or Vinardo score of sulfanilamide or its anion with the carboxylate',
    check=_gv['max_abs_score_change_sfa_species'], tol=1e-6)
add('dock.glu106.vina.max_score_change_azm', max(abs(FACTS[f'dock.glu106.{p_}.{_runkey(r_)[0]}.{_runkey(r_)[1]}.{_runkey(r_)[2]}.delta_score']['value']) for (p_, r_), v_ in _vr if r_.startswith('azm_')), 'kcal mol-1', 3,
    'comparison.json of both variants (largest |carboxylate − acid| of a mean score, acetazolamide, Vina and Vinardo)', 'largest change of a mean Vina or Vinardo score of acetazolamide with the carboxylate', check=_gv['max_abs_score_change_azm'], tol=1e-6)
for pdb_ in ('6RL9', '3HS4'):
    for st_, s17_ in (('acid', 'acid'), ('carb', 'carboxylate')):
        add(f'dock.glu106.{pdb_}.azm_minus_sfa.{st_}', FACTS[f'dock.glu106.{pdb_}.sfa_neutral.ad4zn.val.{st_}.score']['value'] - FACTS[f'dock.glu106.{pdb_}.azm_neutral.ad4zn.val.{st_}.score']['value'], 'kcal mol-1', 2,
            f'dock/{"T3_6RL9" if pdb_ == "6RL9" else "T4_3HS4"}_glu106GLU/comparison.json (validation runs, mean scores)', f'{pdb_}, AutoDock4Zn validation runs, {STN_[st_]}: acetazolamide scores better than sulfanilamide by this',
            check=D17D['glu106_azm_minus_sfa'][pdb_][s17_], tol=1e-6)
_sfa_ad = [(p_, r_) for (p_, r_) in GLU_ if r_.startswith('sfa_neutral__ad4zn')]
for st_ in ('acid', 'carb'):
    ks_ = [f'dock.glu106.{p_}.sfa_neutral.ad4zn.{_runkey(r_)[2]}.{st_}' for p_, r_ in _sfa_ad]
    add(f'dock.glu106.sfa_ad4zn.{st_}.sep_max', max(FACTS[k_ + '.sep_max']['value'] for k_ in ks_), 'kcal mol-1', 3, 'comparison.json of both variants: orientation.sep_max of the four AutoDock4Zn runs of sulfanilamide (largest)',
        f'largest score separation of the two orientations in a seed of the four AutoDock4Zn runs of sulfanilamide, {STN_[st_]}')
    add(f'dock.glu106.sfa_ad4zn.{st_}.zn_n_min', min(FACTS[k_ + '.zn_n']['value'] for k_ in ks_), 'A', 2, 'comparison.json of both variants: Zn–N run means of the four AutoDock4Zn runs of sulfanilamide (smallest)', f'smallest run mean of Zn–N over the four AutoDock4Zn runs of sulfanilamide, {STN_[st_]}')
    add(f'dock.glu106.sfa_ad4zn.{st_}.zn_n_max', max(FACTS[k_ + '.zn_n']['value'] for k_ in ks_), 'A', 2, 'comparison.json of both variants: Zn–N run means of the four AutoDock4Zn runs of sulfanilamide (largest)', f'largest run mean of Zn–N over the four AutoDock4Zn runs of sulfanilamide, {STN_[st_]}')
    add(f'dock.glu106.sfa_ad4zn.{st_}.worst_rank_near', max(FACTS[k_ + '.worst_rank_near']['value'] for k_ in ks_), None, 0, 'comparison.json of both variants: orientation.worst_rank_near (largest over the four runs)',
        f'worst rank of the first near-native mode in any seed of the four AutoDock4Zn runs of sulfanilamide, {STN_[st_]} (a near-native mode is in every seed)')
# fixed poses scored in both receptors (no search)
_fx = D17D['fixed_poses']
_gains = []
for var_ in ('T3_6RL9_glu106GLU', 'T4_3HS4_glu106GLU'):
    pdb_ = var_.split('_')[1]
    rs_ = J(DOCK / var_ / 'rescore_fixed_poses.json')
    for run_, rec_ in rs_['runs'].items():
        kind_ = _runkey(run_)[2]
        s_, d_ = rec_['summary'], _fx[f'glu106/{pdb_}/{run_}']
        b_ = f'dock.glu106.fixed.{pdb_}.{kind_}'
        src_ = f'dock/{var_}/rescore_fixed_poses.json:runs.{run_}.summary'
        lab_ = f'stored modes of the {pdb_} {"production" if kind_ == "prod" else "validation"} run of sulfanilamide scored without a search in the maps of both receptors (AutoDock4Zn)'
        _a(b_ + '.near', s_['near-native']['change_mean'], 'kcal mol-1', 3, f'{src_}["near-native"].change_mean', f'{lab_}: change of the near-native modes with the carboxylate, mean', check=d_['near_native']['mean'], tol=0.0006)
        _a(b_ + '.rot', s_['ring-rotated']['change_mean'], 'kcal mol-1', 3, f'{src_}["ring-rotated"].change_mean', f'{lab_}: change of the modes beyond 2.0 Å (ring-rotated), mean', check=d_['ring_rotated']['mean'], tol=0.0006)
        _a(b_ + '.rot_minus_near', s_['ring_rotated_minus_near_native_change'], 'kcal mol-1', 3, f'{src_}.ring_rotated_minus_near_native_change', f'{lab_}: gain of the ring-rotated orientation over the near-native one with the carboxylate',
           check=d_['rotated_minus_near'], tol=0.0006)
        _gains.append(abs(s_['ring_rotated_minus_near_native_change']))
add('dock.glu106.fixed.rot_gain_min', min(_gains), 'kcal mol-1', 3, 'rescore_fixed_poses.json of both variants: ring_rotated_minus_near_native_change (smallest in size over the four runs)',
    'fixed poses, AutoDock4Zn: smallest gain of the ring-rotated orientation over the near-native one with the carboxylate')
add('dock.glu106.fixed.rot_gain_max', max(_gains), 'kcal mol-1', 3, 'rescore_fixed_poses.json of both variants: ring_rotated_minus_near_native_change (largest in size over the four runs)',
    'fixed poses, AutoDock4Zn: largest gain of the ring-rotated orientation over the near-native one with the carboxylate')
_rv = [J(DOCK / v_ / 'rescore_fixed_poses_vina.json') for v_ in ('T3_6RL9_glu106GLU', 'T4_3HS4_glu106GLU')]
add('dock.glu106.fixed.vina.n_poses', sum(x_['n'] for r_ in _rv for x_ in r_['runs'].values()), None, 0, 'dock/T3_6RL9_glu106GLU/rescore_fixed_poses_vina.json, dock/T4_3HS4_glu106GLU/rescore_fixed_poses_vina.json:runs.*.n (sum)',
    'top poses of stored Vina and Vinardo runs of sulfanilamide scored without a search in both receptors', check=D17D['fixed_poses_vina_total']['n_poses'], tol=0.01)
add('dock.glu106.fixed.vina.max_change', max(max(abs(x_['change_min']), abs(x_['change_max'])) for r_ in _rv for x_ in r_['runs'].values()), 'kcal mol-1', 3,
    'rescore_fixed_poses_vina.json of both variants: runs.*.change_min, change_max (largest in size)', 'largest change of the Vina or Vinardo score of a fixed pose with the carboxylate (these functions do not see the proton)',
    check=D17D['fixed_poses_vina_total']['max_abs_change'], tol=1e-6)
for var_ in ('T3_6RL9_glu106GLU', 'T4_3HS4_glu106GLU'):
    pdb_ = var_.split('_')[1]
    vj_ = J(DOCK / var_ / 'variant.json')
    ge_ = vj_['geometry']
    lg_ = list(ge_['distance_to_co_crystal_ligand'].values())[0]
    add(f'dock.glu106.{pdb_}.he2_to_ligand', lg_['HE2']['distance_A'], 'A', 2, f'dock/{var_}/variant.json:geometry.distance_to_co_crystal_ligand[*].HE2.distance_A',
        f'{pdb_} receptor of the paper: distance from the proton HE2 on OE2 of Glu106 to the nearest atom of the crystal ligand (the sulfonamide nitrogen)')
    add(f'dock.glu106.{pdb_}.oe1_to_ligand', lg_['OE1']['distance_A'], 'A', 2, f'dock/{var_}/variant.json:geometry.distance_to_co_crystal_ligand[*].OE1.distance_A',
        f'{pdb_} receptor: distance from OE1 of Glu106 to the nearest atom of the crystal ligand')
    add(f'dock.glu106.{pdb_}.zn_he2', ge_['zinc_distances_A']['HE2'], 'A', 2, f'dock/{var_}/variant.json:geometry.zinc_distances_A.HE2', f'{pdb_} receptor of the paper: distance from the zinc to the proton HE2 on OE2 of Glu106')
    add(f'dock.glu106.{pdb_}.thr199_hg1_oe1', ge_['polar_H_within_3p2A_of_the_carboxylate'][0]['H_to_OE1_A'], 'A', 2, f'dock/{var_}/variant.json:geometry.polar_H_within_3p2A_of_the_carboxylate[0].H_to_OE1_A',
        f'{pdb_} receptor: Thr199 HG1···OE1(Glu106), the hydrogen bond Thr199 → Glu106 present in both protonation states')
    add(f'dock.glu106.{pdb_}.lines_identical', vj_['receptor.pdbqt']['file_check']['n_lines_source'] - vj_['receptor.pdbqt']['file_check']['n_lines_removed_or_rewritten'], None, 0,
        f'dock/{var_}/variant.json:receptor.pdbqt.file_check (n_lines_source − n_lines_removed_or_rewritten)', f'{pdb_}: lines of receptor.pdbqt that are byte-identical in the two receptors (of {vj_["receptor.pdbqt"]["file_check"]["n_lines_source"]:,})')
    add(f'dock.glu106.{pdb_}.lines_total', vj_['receptor.pdbqt']['file_check']['n_lines_source'], None, 0, f'dock/{var_}/variant.json:receptor.pdbqt.file_check.n_lines_source',
        f'{pdb_}: lines of receptor.pdbqt of the receptor of the paper (one line, the proton HE2, is removed and five are rewritten for the carboxylate)')
    assert vj_['receptor.pdbqt']['file_check']['n_lines_removed_or_rewritten'] == 6

# ============================================================================================= four crystal waters
WREC_ = (('acid', 'T3_6RL9', 'T3_6RL9_wat4', 'hbond', 'acid'), ('carb', 'T3_6RL9_glu106GLU', 'T3_6RL9_glu106GLU_wat4', 'hbond', 'carboxylate'),
         ('rev', 'T3_6RL9', 'T3_6RL9_wat4rev', 'rev', 'acid'), ('bare', 'T3_6RL9', 'T3_6RL9_wat4bare', 'bare', 'acid'))
WLAB_ = {'acid': 'four waters, Glu106 acid', 'carb': 'four waters, Glu106 carboxylate', 'rev': 'four waters with the hydrogen-bond ring reversed, Glu106 acid', 'bare': 'four waters as bare oxygens, Glu106 acid'}


def _modefacts(b_, ms_, own_, src_, lab_):
    m_ = ms_['modes']
    add(b_ + '.n_modes', ms_['n_modes'], None, 0, f'{src_}.mode_summary.n_modes', f'{lab_}: modes written over all seeds', check=own_['modes']['n'], tol=0.01)
    add(b_ + '.modes_near', m_['near_native'], None, 0, f'{src_}.mode_summary.modes.near_native', f'{lab_}: near-native modes', check=own_['modes']['near_native'], tol=0.01)
    add(b_ + '.modes_rot', m_['rotated_pose'], None, 0, f'{src_}.mode_summary.modes.rotated_pose', f'{lab_}: modes of the rotated pose (beyond 2.0 Å from the crystal ligand and within 2.0 Å of the reference rotated pose)',
        check=own_['modes']['rotated_pose'], tol=0.01)
    add(b_ + '.modes_beyond2', m_['ring_rotated_paper'], None, 0, f'{src_}.mode_summary.modes.ring_rotated_paper', f'{lab_}: modes beyond 2.0 Å from the crystal ligand', check=own_['modes']['beyond_2'], tol=0.01)
    add(b_ + '.modes_between', m_['between'], None, 0, f'{src_}.mode_summary.modes.between', f'{lab_}: modes between 1.0 and 2.0 Å', check=own_['modes']['between'], tol=0.01)
    add(b_ + '.seeds_rot', ms_['seeds_holding']['rotated_pose'], None, 0, f'{src_}.mode_summary.seeds_holding.rotated_pose', f'{lab_}: seeds that hold a mode of the rotated pose', check=own_['modes']['seeds_holding_rotated_pose'], tol=0.01)
    add(b_ + '.top_rot', ms_['top_pose']['rotated_pose'], None, 0, f'{src_}.mode_summary.top_pose.rotated_pose', f'{lab_}: seeds whose top pose is the rotated pose', check=own_['top']['rotated_pose'], tol=0.01)
    if ms_['best_score']['near_native'] is not None:
        add(b_ + '.best_near', ms_['best_score']['near_native'], 'kcal mol-1', 3, f'{src_}.mode_summary.best_score.near_native', f'{lab_}: best score of a near-native mode', check=own_['modes']['best_score_near_native'], tol=1e-6)
    if ms_['best_score']['rotated_pose'] is not None:
        add(b_ + '.best_rot', ms_['best_score']['rotated_pose'], 'kcal mol-1', 3, f'{src_}.mode_summary.best_score.rotated_pose', f'{lab_}: best score of a mode of the rotated pose', check=own_['modes']['best_score_rotated_pose'], tol=1e-6)
    w_ = ms_['min_dist_water_site']
    add(b_ + '.wat_min', w_['all_modes'], 'A', 2, f'{src_}.mode_summary.min_dist_water_site.all_modes', f'{lab_}: closest approach of a heavy atom of any mode to one of the four water sites', check=own_['water_sites']['min_dist_any_mode'], tol=0.006)
    add(b_ + '.wat_n24', w_['modes_within_2p4A_of_a_site'], None, 0, f'{src_}.mode_summary.min_dist_water_site.modes_within_2p4A_of_a_site', f'{lab_}: modes with a heavy atom within 2.4 Å of a water site', check=own_['water_sites']['modes_within_2p4'], tol=0.01)


WAT_ = {}
for rc_, dry_, var_, rule_, s17_ in WREC_:
    comp_ = J(DOCK / var_ / 'comparison_water.json')
    for run_, rec_ in comp_['runs'].items():
        lig_, sc_, kind_ = _runkey(run_)
        WAT_[(rc_, run_)] = rec_
        own_ = D17D['wat4'][f'{s17_}/{rule_}/{run_}']
        src_ = f'dock/{var_}/comparison_water.json:runs.{run_}'
        lab_ = f'{_runlab("6RL9", run_, rec_["water"]["n_seeds"])}, {WLAB_[rc_]}'
        b_ = f'dock.wat4.{rc_}.{lig_}.{sc_}.{kind_}'
        _dockfacts(b_, rec_['water'], own_['water'], src_ + '.water', lab_, True)
        _modefacts(b_, rec_['water']['mode_summary'], own_['water'], src_ + '.water', lab_)
        if rc_ in ('acid', 'carb'):          # the dry run of the same receptor state: the mode counts that need the water sites
            _modefacts(f'dock.wat4.dry.{rc_}.{lig_}.{sc_}.{kind_}', rec_['dry']['mode_summary'], own_['dry'], src_ + '.dry',
                       f'{_runlab("6RL9", run_, rec_["dry"]["n_seeds"])}, without water, {STN_[rc_]}')
_wc = D17D['wat4_counts']
_rr = J(DOCK / 'T3_6RL9_wat4/comparison_water.json')['rotated_reference']
add('dock.wat4.rot_ref', 'seed %d, mode %d' % (int(_rr['file'].split('/seed')[1].split('/')[0]), _rr['mode']), None, None, 'dock/T3_6RL9_wat4/comparison_water.json:rotated_reference',
    'reference rotated pose: this mode of the production run of sulfanilamide in the dry receptor of the paper (AutoDock4Zn, DFT geometry); it is also the pose from which start R of the simulations was built')
add('dock.wat4.n_runs', len(WAT_), None, 0, 'comparison_water.json of the four water receptors: runs', 'docking runs in a receptor with the four crystal waters (byte-identical ligand inputs to the dry runs)', check=_wc['n_runs'], tol=0.01)
add('dock.wat4.n_runs_ad4zn', sum(1 for k_ in WAT_ if '__ad4zn' in k_[1]), None, 0, 'comparison_water.json of the four water receptors: runs (AutoDock4Zn)',
    'AutoDock4Zn runs with the four waters (both states of Glu106; three placements of the water hydrogens for the production run)', check=_wc['n_runs_ad4zn'], tol=0.01)
add('dock.wat4.ad4zn.modes_rot_total', sum(v_['water']['mode_summary']['modes']['rotated_pose'] for k_, v_ in WAT_.items() if '__ad4zn' in k_[1]), None, 0,
    'comparison_water.json of the four water receptors: runs.*.water.mode_summary.modes.rotated_pose (sum, AutoDock4Zn)', 'modes of the rotated pose in the eight AutoDock4Zn runs with the four waters',
    check=_wc['modes_of_the_rotated_pose_in_the_ad4zn_water_runs'], tol=0.01)
add('dock.wat4.all.modes_rot_total', sum(v_['water']['mode_summary']['modes']['rotated_pose'] for v_ in WAT_.values()), None, 0,
    'comparison_water.json of the four water receptors: runs.*.water.mode_summary.modes.rotated_pose (sum)', 'modes of the rotated pose in all 16 runs with the four waters', check=_wc['modes_of_the_rotated_pose_in_any_water_run'], tol=0.01)
_vz = [FACTS[f'dock.wat4.{rc_}.{_runkey(r_)[0]}.{_runkey(r_)[1]}.{_runkey(r_)[2]}.zn_n']['value'] for (rc_, r_) in WAT_ if '__vin' in r_] + \
      [FACTS[f'dock.glu106.6RL9.{_runkey(r_)[0]}.{_runkey(r_)[1]}.{_runkey(r_)[2]}.{rc_}.zn_n']['value'] for (rc_, r_) in WAT_ if '__vin' in r_]
add('dock.wat4.vina.zn_n_min', min(_vz), 'A', 2, 'comparison_water.json and comparison.json: Zn–N run means of the Vina and Vinardo runs of sulfanilamide, without and with the waters (smallest)',
    'smallest run mean of Zn–N in the 16 Vina and Vinardo runs of sulfanilamide in 6RL9 (dry and with the four waters): beyond the 2.3 Å cut-off in every seed', check=_wc['vina_zn_n_mean_range'][0], tol=0.006)
add('dock.wat4.vina.zn_n_max', max(_vz), 'A', 2, 'comparison_water.json and comparison.json: Zn–N run means of the Vina and Vinardo runs of sulfanilamide, without and with the waters (largest)',
    'largest run mean of Zn–N in the same 16 runs', check=_wc['vina_zn_n_mean_range'][1], tol=0.006)
assert len(_vz) == 16 and sum(1 for (rc_, r_) in WAT_ if '__vin' in r_ and r_.endswith('__prod')) == 4
add('dock.wat4.vina.n_runs', len(_vz), None, 0, 'comparison_water.json and comparison.json: the Vina and Vinardo runs of sulfanilamide in 6RL9, without and with the waters',
    'Vina and Vinardo runs of sulfanilamide in 6RL9 over which the run means of Zn–N are taken: eight with the four waters and their eight dry counterparts '
    '(in each eight, four from the ETKDG conformer, which Table S36 lists, and four from the DFT geometry); added on 11 Oct 2026')
add('dock.wat4.vina.n_seeds_within_1A', _wc['vina_seeds_within_1A'], None, 0, f'{D17S}:dock.wat4_counts.vina_seeds_within_1A (raw out.pdbqt of the 16 runs)', 'seeds of the same 16 Vina and Vinardo runs whose top pose lies within 1.0 Å of the crystal ligand')
# the four waters in the deposited entry
import csv as _csv  # noqa: E402
with open(MDR_ / 'system/audit/waters.csv', newline='') as _fh:
    _wall = list(_csv.DictReader(_fh))              # one row per oxygen position (alternate positions have rows of their own)
_wrows = {int(x_['water']): x_ for x_ in _wall}
assert all(_wrows[w_]['label'] == '-' for w_ in (423, 516, 553, 568))
_wd = D17D['waters']
_wrem = {x_['water']: x_ for x_ in DEC_['starts']['R']['waters_removed']}
for w_ in (423, 516, 553, 568):
    x_, d_ = _wrows[w_], _wd['sites'][str(w_)]
    add(f'dock.wat4.water.{w_}.resid', int(x_['water']), None, 0, f'md/system/audit/waters.csv:water {w_}, column water', f'residue number of crystal water {w_} in 6RL9 (chain A)',
        check=int(d_['number_in_the_entry']), tol=0.01)
    add(f'dock.wat4.water.{w_}.B', float(x_['B']), 'A2', 1, f'md/system/audit/waters.csv:water {w_}, column B', f'crystal water {w_} of 6RL9: B factor', check=d_['B'], tol=0.06)
    add(f'dock.wat4.water.{w_}.occ', float(x_['occupancy']), None, 2, f'md/system/audit/waters.csv:water {w_}, column occupancy', f'crystal water {w_} of 6RL9: occupancy', check=d_['occupancy'], tol=0.001)
    add(f'dock.wat4.water.{w_}.to_ligand', float(x_['min_dist_SAN304_A']), 'A', 2, f'md/system/audit/waters.csv:water {w_}, column min_dist_SAN304_A',
        f'crystal water {w_}: distance of its oxygen from the nearest heavy atom of the crystal ligand ({x_["nearest_SAN304_atom"]})', check=d_['to_crystal_ligand'], tol=0.006)
    add(f'dock.wat4.water.{w_}.to_rotated', _wrem[w_]['dist_A'], 'A', 2, f'md/system/decisions.json:starts.R.waters_removed[water {w_}].dist_A',
        f'crystal water {w_}: distance of its oxygen from the nearest heavy atom of the docked rotated pose ({_wrem[w_]["ligand_atom"]}); the rotated pose lies on the site', check=d_['to_rotated_pose'], tol=0.006)
_Bs = sorted(float(x_['B']) for x_ in _wall)
add('dock.wat4.water.B_min', min(float(_wrows[w_]['B']) for w_ in (423, 516, 553, 568)), 'A2', 1, 'md/system/audit/waters.csv: column B of waters 423, 516, 553, 568 (smallest)', 'smallest B factor of the four waters', check=_wd['B_range_of_the_four'][0], tol=0.06)
add('dock.wat4.water.B_max', max(float(_wrows[w_]['B']) for w_ in (423, 516, 553, 568)), 'A2', 1, 'md/system/audit/waters.csv: column B of waters 423, 516, 553, 568 (largest)', 'largest B factor of the four waters', check=_wd['B_range_of_the_four'][1], tol=0.06)
add('dock.wat4.water.B_median_all', _stt.median(_Bs), 'A2', 1, 'md/system/audit/waters.csv: column B (median over all water sites of chain A)', f'median B factor of the {len(_Bs)} water oxygen positions of chain A of 6RL9 ({len(_wrows)} sites; alternate positions counted)', check=_wd['median_B_all_waters'], tol=0.06)
add('dock.wat4.water.n', 4, None, 0, 'dock/T3_6RL9_wat4/variant.json:waters', 'crystal waters kept in the receptor of the test (423, 516, 553 and 568, all at full occupancy; always kept together)',
    check=len(J(DOCK / 'T3_6RL9_wat4/variant.json')['waters']), tol=0.01)
# stored modes of the dry runs scored in the maps with the waters (no search)
for var_, tag_ in (('T3_6RL9_wat4', 'acid'), ('T3_6RL9_glu106GLU_wat4', 'carb'), ('T3_6RL9_wat4rev', 'rev'), ('T3_6RL9_wat4bare', 'bare'), ('T3_6RL9_wat4fixed', 'fixed')):
    d_ = _fx[f'wat4/{var_}']
    lab_ = f'stored modes of the dry production run scored without a search in the maps with the four waters ({WLAB_.get(tag_, "four waters in one fixed orientation, Glu106 acid")})'
    src_ = f'{D17S}:dock.fixed_poses["wat4/{var_}"] (per-mode values of dock/{var_}/rescore_fixed_poses.json:runs.sfa_neutral__ad4zn__prod.modes)'
    add(f'dock.wat4.fixed.{tag_}.near_change', d_['near_native']['change_mean'], 'kcal mol-1', 3, src_, f'{lab_}: change of the best near-native mode of each seed, mean (the waters do not push the crystal pose out)')
    add(f'dock.wat4.fixed.{tag_}.rot_mean', d_['rotated_pose']['mean'], 'kcal mol-1', 0, src_, f'{lab_}: score of the best mode of the rotated pose of each seed, mean (it overlaps the waters)')
    add(f'dock.wat4.fixed.{tag_}.rot_min', d_['rotated_pose']['min'], 'kcal mol-1', 0, src_, f'{lab_}: the same, smallest')
    add(f'dock.wat4.fixed.{tag_}.rot_max', d_['rotated_pose']['max'], 'kcal mol-1', 0, src_, f'{lab_}: the same, largest')
add('dock.wat4.fixed.rot_min_all', D17D['fixed_poses_rotated_range_all_rules'][0], 'kcal mol-1', 0, f'{D17S}:dock.fixed_poses_rotated_range_all_rules[0]', 'smallest score of a stored mode of the rotated pose in the maps with the four waters, over all four placements of the water hydrogens and both states of Glu106')
add('dock.wat4.fixed.rot_max_all', D17D['fixed_poses_rotated_range_all_rules'][1], 'kcal mol-1', 0, f'{D17S}:dock.fixed_poses_rotated_range_all_rules[1]', 'largest score of a stored mode of the rotated pose in the maps with the four waters')

# ------------------------------------------------------------------------------------ literature value (11 Oct 2026)
# One value of the literature that Section S7 of the supplementary material prints since the last reading of the paper.
# It stands at the end of this section because this plug-in is the last one before the closing check.
add('lit.husak2022.hbond_long.A', 2.613, 'A', 3, 'husak2022potential (Acta Crystallogr., Sect. B 2022, 78, 781), abstract: "When both the cocrystal and salt models converge to salt for a long '
    'hydrogen-bond (longer than 2.613 Å) the structure is most likely salt" (Crossref record, read 11 Oct 2026)',
    'length of the hydrogen bond above which a relaxation that ends in the salt from both starts was found to indicate a salt, in a test of dispersion-corrected PBE on 173 salts and 96 cocrystals '
    '(literature; printed in Section S7 beside N1···O5 of SAPH, 2.734 Å in the relaxed cell and 2.760 Å in the X-ray structure)')
