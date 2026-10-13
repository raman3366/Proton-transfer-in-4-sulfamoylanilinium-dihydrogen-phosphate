"""Facts of the calculations run for the referees (VM burst 3, 29 Sep 2026; sections 32-36 of facts.md), executed by
build_facts.py after facts_revision.py (exec in its namespace: add(), section(), FACTS, J, r, RES, DOCK, DER, CAVEATS,
MISSING, DISAGREE and HA_KCAL are the ones defined there).

  32. M2/M3  Constrained periodic N1-H1A profile (PBE-D3(BJ) and PBE0-D3(BJ)) and the single-transfer control
             (periodic_qe/profile.json; re-derived from the raw pw.x outputs by derive/d13_referee_calcs.py).
  33. M4     Hessian and thermochemistry of the neutral IEF-PCM minimum (sfa_pair_pcm_d3_oform/freq.json; own
             diagonalisation and RRHO of both minima in d13).
  34. M6     QTAIM of the nine crystal hydrogen bonds in a five-pair cluster cut from the periodic cells
             (crystal_cluster_qtaim/qtaim_nine_bonds.json; bond critical points re-assigned by geometry from the raw
             Multiwfn CPprop.txt in d13).
  35. M9     Relaxed C3-C4-S1-N2 torsion scans of sulfanilamide and the sulfonamidate (<case>/torsion_scan/,
             torsion_scan/torsion_strain.json; energies and torsions re-read point by point and the spline redone in d13).
  36. Superseded facts: the constrained single configuration (relax_iirfix) is relabelled, not withdrawn; the rigid DFT
      strain of the docked poses (sections 29 and 31) and the order-of-magnitude checks on relax_iirfix (section 31) are
      retired: value null, the number kept in superseded_value, so numbers_audit.py cannot match them.
The user decided (29 Sep 2026) that the cell-relaxed end points, the molecular PBE calibration and the M5 solvation
checks are not run; they are recorded under "Missing" as stated limitations.
"""
D13 = J(DER / 'd13_referee_calcs.json')
D13S = 'derived/d13_referee_calcs.json'
PRJ = J(RES / 'periodic_qe/profile.json')
PRS = 'molecules/sfa/results/periodic_qe/profile.json'
RY_KJ_FU_ = 13.605693 * 96.485332 / 2
PROF_ = ['prof110', 'prof125', 'prof140', 'prof155', 'prof170']
PKEY = {p: 'p' + p[4:] for p in PROF_}                 # prof110 -> p110
NEW_KEYS_START = len(FACTS)

# =====================================================================================================================
section('32. Periodic N1–H1A profile (referee M2/M3) and the single-transfer control (VM burst 3)',
        'd(N1–H1A) of both formula units held by a pw.x distance constraint at 1.10, 1.25, 1.40, 1.55 and 1.70 Å; H1A across the '
        'bond, O5 and every other atom free; fixed X-ray cell; PBE-D3(BJ)/SSSP, 70/560 Ry, 2×1×1; damped dynamics (pw.x accepts '
        'constraints only with damped or FIRE dynamics), converged on the constrained forces below 10⁻³ Ry Bohr⁻¹. Each point '
        'seeded by interpolating relaxed cell (i) towards the earlier constrained configuration. PBE0-D3(BJ) = the PBE-D3(BJ)/SSSP '
        'difference plus the hybrid correction from SG15 single points (100/400 Ry, q 1×1×1) at the same geometries, as for the '
        'single configuration of section 28. Slope dE/dd = (F_N1 − F_H1A)·u/2 per formula unit from the raw forces at the '
        'constrained minimum (u: N1→H1A): positive means the constraint holds H1A against a force back towards N1. Single-transfer '
        'control relax_iir1: relaxed cell (i) with H1A put 0.98 Å from O5 in formula unit 1 only, unconstrained BFGS without '
        'symmetry. Every number re-derived from the raw pw.x outputs (derive/d13_referee_calcs.py): final energies, final '
        'coordinates and cell, raw and constrained force blocks, pressures, the CONSTRAINTS card of each input.')
PD = D13['periodic']
PR = PD['runs']
add('periodic.profile.method', PRJ['method']['what'], None, None, f'{PRS}:method.what', 'periodic N1–H1A profile protocol', case='periodic_qe')
add('periodic.profile.optimiser', 'pw.x damped dynamics (ion_dynamics = damp) with distance constraints; QE 7.5 refuses CONSTRAINTS with BFGS',
    None, None, f'{PRS}:method.optimiser; runs/relax_prof*/relax_prof*.in (ion_dynamics)', 'optimiser of the constrained profile points', case='periodic_qe')
add('periodic.profile.targets', [PR[p]['constraints'][0]['target_A'] for p in PROF_], 'A', 2,
    'runs/relax_prof*/relax_prof*.in: CONSTRAINTS distance targets (Bohr → Å), atoms N1–H1A and N1_2–H1A_2',
    'held N1–H1A distances of the profile', case='periodic_qe')
add('periodic.profile.forc_conv_thr', PR['prof110']['forc_conv_thr'], 'Ry/bohr', 3, 'runs/relax_prof*/relax_prof*.in: forc_conv_thr',
    'force threshold of the constrained points (on the constrained forces; the unconstrained relaxations used 5 × 10⁻⁴)', case='periodic_qe')
add('periodic.profile.n_points', len(PROF_), None, 0, f'{PRS}:neutral_minimum.complete', 'constrained profile points (PBE-D3(BJ) and PBE0-D3(BJ) each)')
add('periodic.profile.all_converged', PD['verdict']['all_converged'], None, None, f'{D13S}:periodic.verdict.all_converged ("Damped Dynamics: convergence achieved")',
    'every profile point converged', case='periodic_qe')
add('periodic.profile.max_con_force', PD['verdict']['max_constrained_force_Ry_bohr'], 'Ry/bohr', 4, f'{D13S}:periodic.verdict (last "Constrained forces" blocks)',
    'largest constrained force component at the end of any profile point', case='periodic_qe',
    check=max(PRJ['points'][p]['max_constrained_force_Ry_bohr'] for p in PROF_), tol=1e-6)
add('periodic.profile.held_dev', PD['verdict']['held_distance_max_dev_A'], 'A', 5, f'{D13S}:periodic.verdict.held_distance_max_dev_A',
    'largest deviation of a held N1–H1A from its target in the final coordinates', case='periodic_qe')
add('periodic.profile.fu_identical', PD['verdict']['both_fu_identical'], 'A', 6, f'{D13S}:periodic.verdict.both_fu_identical',
    'largest difference between the two formula units in N1–H1A, H1A···O5, N1···O5 and O4···O5 at any point (the profile keeps the Pc symmetry)', case='periodic_qe')
add('periodic.profile.ref_E', PRJ['relax_i']['energy_Ry'], 'Ry', 7, f'{PRS}:relax_i.energy_Ry', 'energy of the salt reference, relaxed cell (i)', case='periodic_qe',
    check=PR['relax_i']['E_Ry'], tol=1e-9)
for p in PROF_:
    k, pj, rr, hy = PKEY[p], PRJ['points'][p], PR[p], PD['hybrid'][p]
    t = f'{rr["constraints"][0]["target_A"]:.2f}'
    src = f'{PRS}:points.{p}'
    g1 = rr['fu']['_1']
    pf = pj['formula_units'][0]
    # geometry (Table 2 keys: proton.periodic.<k>.*)
    add(f'proton.periodic.{k}.d_N1_H1A', g1['d_N1_H1A'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1 (relax_{p}.out, final coordinates)',
        f'periodic profile, N1–H1A held at {t} Å: N1–H1A', case='periodic_qe', check=pf['d_N1_H1A_A'], tol=6e-5)
    add(f'proton.periodic.{k}.d_H1A_O5', g1['d_H1A_O5'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1', f'periodic profile, N1–H1A held at {t} Å: H1A···O5',
        case='periodic_qe', check=pf['d_H1A_O5_A'], tol=6e-5)
    add(f'proton.periodic.{k}.delta', g1['delta'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1', f'periodic profile, N1–H1A held at {t} Å: δ',
        case='periodic_qe', check=pf['delta_A'], tol=6e-5)
    add(f'proton.periodic.{k}.d_N1_O5', g1['d_N1_O5'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1', f'periodic profile, N1–H1A held at {t} Å: N1···O5',
        case='periodic_qe', check=pf['d_N1_O5_A'], tol=6e-5)
    add(f'proton.periodic.{k}.P1O5', g1['d_P1_O5'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1', f'periodic profile, N1–H1A held at {t} Å: P1–O5',
        case='periodic_qe', check=pf['d_P1_O5_A'], tol=6e-5)
    add(f'proton.periodic.{k}.side', g1['side'], None, None, f'{D13S}:periodic.runs.{p}.fu._1.side (shorter of N1–H1A and H1A···O5)',
        f'periodic profile, N1–H1A held at {t} Å: proton nearer N1 or O5', case='periodic_qe')
    add(f'periodic.profile.{k}.d_O4O5', g1['d_O4_O5'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1.d_O4_O5', f'profile point {t} Å: O4···O5 of the anion–anion bond',
        case='periodic_qe', check=pf['O4_H4A_O5']['d_O4_O5_A'], tol=6e-5)
    add(f'periodic.profile.{k}.d_O4H4A', g1['d_O4_H4A'], 'A', 3, f'{D13S}:periodic.runs.{p}.fu._1.d_O4_H4A', f'profile point {t} Å: O4–H4A',
        case='periodic_qe', check=pf['O4_H4A_O5']['d_O4_H4A_A'], tol=6e-5)
    add(f'periodic.profile.{k}.ang', g1['ang_N1_H1A_O5'], 'deg', 1, f'{D13S}:periodic.runs.{p}.fu._1.ang_N1_H1A_O5', f'profile point {t} Å: N1–H1A–O5 angle',
        case='periodic_qe', check=pf['angle_N1_H1A_O5_deg'], tol=0.006)
    for hb, x in rr['hbonds']['_1'].items():
        add(f'periodic.profile.{k}.hb.{hb}.DA', x['DA'], 'A', 3, f'{D13S}:periodic.runs.{p}.hbonds._1.{hb}',
            f'profile point {t} Å: D···A of {hb.replace("_", "···")} (formula unit 1; formula unit 2 identical)', case='periodic_qe')
    # energies
    add(f'periodic.profile.{k}.E', pj['energy_Ry'], 'Ry', 7, f'{src}.energy_Ry', f'energy of profile point {t} Å (PBE-D3(BJ)/SSSP)', case='periodic_qe',
        check=rr['E_Ry'], tol=1e-9)
    dE = (pj['energy_Ry'] - PRJ['relax_i']['energy_Ry']) * RY_KJ_FU_
    add(f'periodic.profile.{k}.PBE_meV', dE * 2 / 96.485332 * 1000, 'meV', 1, f'{src}.energy_Ry − relax_i.energy_Ry (JSON: {pj["E_minus_i_meV_cell"]} meV/cell)',
        f'profile point {t} Å − salt, PBE-D3(BJ), per cell', case='periodic_qe', check=rr['dE_meV_cell'], tol=0.001)
    add(f'periodic.profile.{k}.PBE_kJ', dE, 'kJ mol-1', 1, f'{src}.energy_Ry − relax_i.energy_Ry (JSON: {pj["E_minus_i_kJ_mol_per_formula_unit"]})',
        f'profile point {t} Å − salt, PBE-D3(BJ), per formula unit', case='periodic_qe', check=rr['dE_kJ_fu'], tol=1e-4)
    add(f'periodic.profile.{k}.PBE_kcal', dE / 4.184, 'kcal mol-1', 1, f'periodic.profile.{k}.PBE_kJ / 4.184', 'the same in kcal mol⁻¹ per formula unit',
        case='periodic_qe', check=rr['dE_kcal_fu'], tol=1e-4)
    hj = PRJ['hybrid']['sg15_e100'][p]
    add(f'periodic.profile.{k}.corr_meV', hy['corr_meV_cell'], 'meV', 1, f'{D13S}:periodic.hybrid.{p}.corr_meV_cell (hybrid_{p}_sg15_e100.out, pbe_{p}_sg15_e100.out and the relax_i pair)',
        f'hybrid correction at {t} Å, per cell', case='periodic_qe', check=hj['hybrid_correction_meV_cell'], tol=0.006)
    add(f'periodic.profile.{k}.corr_kJ', hy['corr_kJ_fu'], 'kJ mol-1', 1, f'{D13S}:periodic.hybrid.{p}.corr_kJ_fu', f'hybrid correction at {t} Å, per formula unit',
        case='periodic_qe', check=hj['hybrid_correction_kJ_mol_per_formula_unit'], tol=0.001)
    add(f'periodic.profile.{k}.corr_D3_meV', hy['corr_D3_meV_cell'], 'meV', 1, f'{D13S}:periodic.hybrid.{p}.corr_D3_meV_cell (first "!" energy of the PBE0 runs)',
        f'part of the correction at {t} Å from the pbe0 D3(BJ) parameters, per cell', case='periodic_qe',
        check=hj['split_of_the_correction_meV_cell']['D3_parameters_pbe0_vs_pbe'], tol=0.006)
    add(f'periodic.profile.{k}.PBE0_kJ', hy['best_kJ_fu'], 'kJ mol-1', 1, f'{D13S}:periodic.hybrid.{p}.best_kJ_fu (SSSP PBE-D3 difference + correction)',
        f'profile point {t} Å − salt, PBE0-D3(BJ), per formula unit', case='periodic_qe',
        check=hj['best_estimate_PBE0_D3BJ_kJ_mol_per_formula_unit'], tol=0.001)
    add(f'periodic.profile.{k}.PBE0_kcal', hy['best_kcal_fu'], 'kcal mol-1', 1, f'periodic.profile.{k}.PBE0_kJ / 4.184', 'the same in kcal mol⁻¹ per formula unit',
        case='periodic_qe')
    add(f'periodic.profile.{k}.PBE0_meV', hy['best_meV_cell'], 'meV', 1, f'{D13S}:periodic.hybrid.{p}.best_meV_cell', f'profile point {t} Å − salt, PBE0-D3(BJ), per cell',
        case='periodic_qe', check=hj['best_estimate_PBE0_D3BJ_meV_cell'], tol=0.006)
    add(f'periodic.profile.{k}.exx_converged', hy['exx_converged'], None, None, f'{D13S}:periodic.hybrid.{p}.exx_converged', f'EXX converged at {t} Å (and in the reference)',
        case='periodic_qe')
    s1 = rr['slope']['_1']
    add(f'periodic.profile.{k}.slope_kJ', pj['slope']['dE_dd_kJ_mol_A_per_formula_unit'], 'kJ mol-1 A-1', 0, f'{src}.slope.dE_dd_kJ_mol_A_per_formula_unit',
        f'dE/dd(N1–H1A) at {t} Å, per formula unit (constraint force; positive = pushed back to N1)', case='periodic_qe',
        check=s1['dEdd_kJ_A_fu'], tol=0.01)
    add(f'periodic.profile.{k}.slope_eVA', s1['dEdd_eV_A'], 'eV/A', 2, f'{D13S}:periodic.runs.{p}.slope._1.dEdd_eV_A', f'the same in eV Å⁻¹ per formula unit',
        case='periodic_qe', check=pj['slope']['dE_dd_per_fu_eV_A'][0], tol=6e-5)
    add(f'periodic.profile.{k}.pressure', pj['pressure_kbar_final'], 'kbar', 1, f'{src}.pressure_kbar_final', f'residual pressure at {t} Å (fixed X-ray cell)', case='periodic_qe',
        check=rr['pressure_kbar'], tol=1e-6)
    add(f'periodic.profile.{k}.steps', pj['steps'], None, 0, f'{src}.steps', f'damped-dynamics steps at {t} Å', case='periodic_qe', check=rr['steps'], tol=0)
    add(f'periodic.profile.{k}.rmsd_heavy', rr['rmsd_heavy_vs_xray'], 'A', 3, f'{D13S}:periodic.runs.{p}.rmsd_heavy_vs_xray',
        f'heavy-atom RMSD of the {t} Å cell from the X-ray structure, mean shift removed', case='periodic_qe')
    add(f'periodic.profile.{k}.rmsd_vs_i', pj['rmsd_all_atoms_vs_relax_i_A'], 'A', 2, f'{src}.rmsd_all_atoms_vs_relax_i_A', f'all-atom RMSD of the {t} Å cell from relaxed cell (i)',
        case='periodic_qe')
V_ = PD['verdict']
pb = [FACTS[f'periodic.profile.{PKEY[p]}.PBE_kJ']['value'] for p in PROF_]
p0 = [FACTS[f'periodic.profile.{PKEY[p]}.PBE0_kJ']['value'] for p in PROF_]
add('periodic.profile.PBE_range', [min(pb), max(pb)], 'kJ mol-1', 1, 'periodic.profile.p110/p170.PBE_kJ', 'PBE-D3(BJ) profile, first and last point, per formula unit', case='periodic_qe')
add('periodic.profile.PBE0_range', [min(p0), max(p0)], 'kJ mol-1', 1, 'periodic.profile.p110/p170.PBE0_kJ', 'PBE0-D3(BJ) profile, first and last point, per formula unit', case='periodic_qe')
add('periodic.profile.monotonic_PBE', V_['monotonic_PBE'], None, None, f'{D13S}:periodic.verdict.monotonic_PBE', 'PBE-D3(BJ) energy rises at every point from the salt to 1.70 Å', case='periodic_qe')
add('periodic.profile.monotonic_PBE0', V_['monotonic_PBE0'], None, None, f'{D13S}:periodic.verdict.monotonic_PBE0', 'PBE0-D3(BJ) energy rises at every point', case='periodic_qe')
add('periodic.profile.slopes_positive', V_['all_slopes_positive'], None, None, f'{D13S}:periodic.verdict.all_slopes_positive',
    'the constraint force points back to N1 at every point, in both formula units', case='periodic_qe')
add('periodic.profile.min_slope', V_['min_slope_kJ_A_fu'], 'kJ mol-1 A-1', 0, f'{D13S}:periodic.verdict.min_slope_kJ_A_fu', 'smallest dE/dd over the profile, per formula unit (at 1.55 Å)',
    case='periodic_qe')
add('periodic.profile.max_slope', max(PR[p]['slope']['_1']['dEdd_kJ_A_fu'] for p in PROF_), 'kJ mol-1 A-1', 0, f'{D13S}:periodic.runs.*.slope', 'largest dE/dd over the profile (at 1.25 Å)',
    case='periodic_qe')
add('periodic.profile.verdict_PBE', PRJ['neutral_minimum']['PBE_D3BJ']['text'], None, None, f'{PRS}:neutral_minimum.PBE_D3BJ.text', 'outcome at PBE-D3(BJ): no neutral minimum', case='periodic_qe')
add('periodic.profile.verdict_PBE0', PRJ['neutral_minimum']['PBE0_D3BJ']['text'], None, None, f'{PRS}:neutral_minimum.PBE0_D3BJ.text', 'outcome at PBE0-D3(BJ) (single points): no neutral minimum', case='periodic_qe')
add('periodic.profile.no_neutral_min', PRJ['neutral_minimum']['PBE_D3BJ']['answer'] == 'no' and PRJ['neutral_minimum']['PBE0_D3BJ']['answer'] == 'no'
    and V_['monotonic_PBE'] and V_['monotonic_PBE0'] and V_['all_slopes_positive'], None, None, f'{PRS}:neutral_minimum + {D13S}:periodic.verdict',
    'no neutral minimum between the salt and d(N1–H1A) = 1.70 Å at PBE-D3(BJ) or PBE0-D3(BJ)', case='periodic_qe')
add('periodic.profile.first_O_side', PR[V_['first_O_side_point']]['constraints'][0]['target_A'], 'A', 2, f'{D13S}:periodic.verdict.first_O_side_point',
    'first profile point with H1A nearer O5 than N1', case='periodic_qe')
add('periodic.profile.min_NO', V_['min_N1_O5_A'], 'A', 3, f'{D13S}:periodic.verdict.min_N1_O5_A', 'shortest N1···O5 along the profile (at 1.40 Å)', case='periodic_qe')
cr_ = [PD['hybrid'][p]['corr_kJ_fu'] for p in PROF_]
add('periodic.profile.corr_range', [min(cr_), max(cr_)], 'kJ mol-1', 1, f'{D13S}:periodic.hybrid.*.corr_kJ_fu', 'range of the hybrid correction over the profile, per formula unit '
    '(largest, +4.1, at 1.25 Å, near the shared-proton region; negative only at 1.70 Å)', case='periodic_qe')
add('periodic.profile.corr_D3_max', max(abs(PD['hybrid'][p]['corr_D3_meV_cell']) for p in PROF_), 'meV', 1, f'{D13S}:periodic.hybrid.*.corr_D3_meV_cell',
    'largest |D3(BJ)-parameter part| of the hybrid correction over the profile, per cell (the rest is exact exchange)', case='periodic_qe')
add('periodic.profile.SG15_minus_SSSP_max', max(abs(PD['hybrid'][p]['SG15_minus_SSSP_meV_cell']) for p in PROF_) / 1000 * 96.485332 / 2, 'kJ mol-1', 1,
    f'{D13S}:periodic.hybrid.*.SG15_minus_SSSP_meV_cell', 'largest pseudopotential difference (SG15 − SSSP) of the PBE-D3(BJ) profile energy, per formula unit',
    case='periodic_qe')
add('periodic.profile.iirfix_excess', V_['iirfix_excess_kJ_fu'], 'kJ mol-1', 1, f'{D13S}:periodic.verdict.iirfix_excess_kJ_fu',
    'the constrained single configuration (O5–H1A held at 0.98 Å, d(N1–H1A) 1.616 Å) lies this much above the profile at the same N1–H1A '
    '(PBE-D3(BJ), cubic Hermite interpolation of the 1.55 and 1.70 Å points with their slopes), per formula unit: the cost of the held O–H and the frozen O5',
    case='periodic_qe', status='own derivation (interpolation between two computed points)')
add('periodic.profile.at_iirfix_d', V_['profile_at_iirfix_d_kJ_fu'], 'kJ mol-1', 1, f'{D13S}:periodic.verdict.profile_at_iirfix_d_kJ_fu',
    'PBE-D3(BJ) profile energy interpolated at d(N1–H1A) = 1.616 Å, per formula unit', case='periodic_qe', status='own derivation (interpolation)')
add('periodic.profile.stress_B10', PD['stress_estimate_kJ_fu']['B10'], 'kJ mol-1', 1, f'{D13S}:periodic.stress_estimate_kJ_fu.B10',
    'relaxation energy p²V/(2B) of the 1.70 Å cell minus that of the salt cell at their residual pressures (+2.4 and −4.4 kbar), B = 10 GPa, per formula unit: '
    'negative, so relaxing the cells would lower the salt more and raise the cost', case='periodic_qe', status='order-of-magnitude estimate, not a calculation')
add('periodic.profile.stress_B20', PD['stress_estimate_kJ_fu']['B20'], 'kJ mol-1', 1, f'{D13S}:periodic.stress_estimate_kJ_fu.B20', 'the same with B = 20 GPa',
    case='periodic_qe', status='order-of-magnitude estimate, not a calculation')
add('periodic.salt.pressure', PR['relax_i']['pressure_kbar'], 'kbar', 1, 'runs/relax_i/relax_i.out (final P=)', 'residual pressure of the salt cell (i) at the fixed X-ray cell', case='periodic_qe',
    check=J(RES / 'periodic_qe/periodic_qe.json')['relax']['i']['pressure_kbar_final'], tol=0.006)
# single-transfer control
ic = PRJ['single_transfer_control']
ri1 = PR['iir1']
add('periodic.iir1.start_OH', ri1['start']['_1']['d_H1A_O5'], 'A', 2, f'{D13S}:periodic.runs.iir1.start._1 (relax_iir1.in)', 'start of the single-transfer control: O5–H1A in formula unit 1',
    case='periodic_qe')
add('periodic.iir1.start_sides', [ri1['start']['_1']['side'], ri1['start']['_2']['side']], None, None, f'{D13S}:periodic.runs.iir1.start',
    'start of the control: proton on O5 in formula unit 1, on N1 in formula unit 2', case='periodic_qe')
add('periodic.iir1.final_sides', V_['iir1_final_sides'], None, None, f'{D13S}:periodic.runs.iir1.fu (relax_iir1.out, final coordinates)',
    'end of the control: proton on N1 in both formula units', case='periodic_qe')
add('periodic.iir1.verdict', ic['verdict'], None, None, f'{PRS}:single_transfer_control.verdict', 'outcome of the single-transfer control', case='periodic_qe')
add('periodic.iir1.dE_kJ', ri1['dE_kJ_fu'], 'kJ mol-1', 2, f'{D13S}:periodic.runs.iir1.dE_kJ_fu', 'E(control, relaxed) − E(i), per formula unit: the same salt',
    case='periodic_qe', check=ic['E_minus_i_kJ_mol_per_formula_unit'], tol=0.006)
add('periodic.iir1.E', ri1['E_Ry'], 'Ry', 7, 'runs/relax_iir1/relax_iir1.out (Final energy)', 'relaxed energy of the single-transfer control',
    case='periodic_qe', check=ic['energy_Ry'], tol=1e-9)
add('periodic.iir1.steps', ri1['steps'], None, 0, 'runs/relax_iir1/relax_iir1.out ("bfgs converged in 51 scf cycles and 49 bfgs steps")', 'BFGS steps of the control',
    case='periodic_qe', check=ic['steps'][1], tol=0)
add('periodic.iir1.n_symops', ri1['n_symops'], None, 0, 'runs/relax_iir1/relax_iir1.out (symmetry operations found)', 'symmetry operations used (none beyond identity: the two formula units were free to differ)',
    case='periodic_qe', check=ic['n_symops'], tol=0)
add('periodic.iir1.d_N1_H1A', ri1['fu']['_1']['d_N1_H1A'], 'A', 3, f'{D13S}:periodic.runs.iir1.fu._1', 'control, formula unit 1: final N1–H1A', case='periodic_qe',
    check=ic['formula_units'][0]['d_N1_H1A_A'], tol=6e-5)
add('periodic.iir1.d_N1_O5', ri1['fu']['_1']['d_N1_O5'], 'A', 3, f'{D13S}:periodic.runs.iir1.fu._1', 'control, formula unit 1: final N1···O5', case='periodic_qe',
    check=ic['formula_units'][0]['d_N1_O5_A'], tol=6e-5)
add('periodic.iir1.delta', ri1['fu']['_1']['delta'], 'A', 3, f'{D13S}:periodic.runs.iir1.fu._1', 'control, formula unit 1: final δ', case='periodic_qe',
    check=ic['formula_units'][0]['delta_A'], tol=6e-5)
add('periodic.iir1.max_force', ri1['max_raw_force_Ry_bohr'], 'Ry/bohr', 5, 'runs/relax_iir1/relax_iir1.out (last force block)', 'largest force component at the end of the control',
    case='periodic_qe', check=ic['max_force_Ry_bohr'], tol=1e-6)
add('periodic.iir1.rmsd_heavy', ri1['rmsd_heavy_vs_xray'], 'A', 3, f'{D13S}:periodic.runs.iir1.rmsd_heavy_vs_xray',
    'heavy-atom RMSD of the relaxed control from the X-ray structure, mean shift removed', case='periodic_qe')
add('periodic.iir1.rmsd_vs_i', ic['rmsd_all_atoms_vs_relax_i_A'], 'A', 3, f'{PRS}:single_transfer_control.rmsd_all_atoms_vs_relax_i_A', 'all-atom RMSD of the relaxed control from cell (i)',
    case='periodic_qe')
CAVEATS.append('Periodic profile (section 32): a constrained path in the fixed X-ray cell, with d(N1–H1A) of both formula units held together '
               '(the per-cell energy moves two protons). Converged on the constrained forces below 10⁻³ Ry/Bohr, twice the threshold of the '
               'unconstrained relaxations. The PBE0-D3(BJ) values are single points at the PBE-D3(BJ) geometries (no hybrid relaxation), '
               'so the PBE0 curve shows the functional dependence of the energies, not of the path. The 1.70 Å point is not a minimum either: '
               'its constraint force still points to N1. No cell-relaxed end point was computed (user decision, 29 Sep 2026): the residual '
               'pressures (−4.4 kbar salt, +2.4 kbar at 1.70 Å) give only an order-of-magnitude stress term, −1.1 kJ/mol per formula unit for '
               'B = 10 GPa, which would raise the cost. The single-transfer control is one unconstrained start at PBE-D3(BJ).')

# SCF thresholds as run (the S1 table and Methods quote them)
add('protocol.qe_conv_thr_profile', 1e-9, 'Ry', None, 'runs/relax_prof*/relax_prof*.in: conv_thr = 1.0d-09', 'SCF threshold of the constrained profile points',
    case='periodic_qe')
add('protocol.conv_tol_pcm_opt', 1e-8, 'Ha', None, 'dft/optimize.py (conv_tol 1e-8 with a solvent), dft/proton_tools.py (conv_tol = 1e-8 if solvent)',
    'SCF threshold of the IEF-PCM optimisations, relaxed scans, dielectric ladder and torsion scans (gas phase 10⁻⁹ Ha, protocol.conv_tol)')
add('protocol.conv_tol_final', 1e-10, 'Ha', None, 'dft/fmo.py, dft/frequencies.py, dft/nlo.py, dft/interaction.py (conv_tol=1e-10)',
    'SCF threshold of the final single points (orbitals, Hessian gradients, NLO fields, interaction energies)')
add('protocol.conv_tol_cluster', J(RES / 'sfa_xclu_salt/scf.json')['level']['conv_tol'], 'Ha', None, 'molecules/sfa/results/sfa_xclu_*/scf.json:level.conv_tol',
    'SCF threshold of the five-pair cluster single points')
DISAGREE.append('Methods 2.2 says "each SCF was converged to 10⁻⁹ Ha" (protocol.conv_tol) and Table S1 printed 10⁻⁹ Ha for every molecular '
                'run until 29 Sep 2026. As run: IEF-PCM optimisations, scans, the dielectric ladder and the torsion scans used 10⁻⁸ Ha '
                '(dft/optimize.py, dft/proton_tools.py), gas-phase optimisations, the functional check, ΔSCF and TD-DFT 10⁻⁹ Ha, and the '
                'final single points (orbitals, Hessians, NLO, interaction energies) 10⁻¹⁰ Ha. Table S1 and the Methods (Section 2.2) now '
                'give the thresholds as run.')

# =====================================================================================================================
section('33. Neutral IEF-PCM minimum: Hessian, thermochemistry and IR (referee M4, VM burst 3)',
        'sfa_pair_pcm_d3_oform: the O5-protonated minimum of the relaxed scan (chain B anchor, identical coordinates), finite-difference '
        'Hessian from 162 analytic gradients (step 0.005 Bohr), B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM (water), as for the salt '
        '(section 9). Differences neutral − salt are like for like: the same RRHO ideal-gas treatment on top of the PCM electronic '
        'energy, unscaled frequencies, no low-mode correction. Re-derived (derive/d13_referee_calcs.py): both Hessians re-diagonalised '
        'with translations and rotations projected out, IR intensities from the dipole derivatives, the proton stretch found by '
        'projection on N1–H1A and O5–H1A, and the thermochemistry of both minima recomputed with one code.')
M4 = D13['m4']
FJo = J(RES / 'sfa_pair_pcm_d3_oform/freq.json')
FJs = J(RES / 'sfa_pair_pcm_d3/freq.json')
OS = 'molecules/sfa/results/sfa_pair_pcm_d3_oform/freq.json'
add('pcm.oform.E', FJo['E_Ha'], 'Ha', 6, f'{OS}:E_Ha', 'electronic energy of the neutral IEF-PCM minimum', case='sfa_pair_pcm_d3_oform',
    check=FACTS['scan.B00.E']['value'], tol=1e-8)
add('pcm.oform.same_as_scan_B', M4['neutral_is_scan_B_anchor'] < 1e-8, None, None, f'{D13S}:m4.neutral_is_scan_B_anchor (max coordinate difference {M4["neutral_is_scan_B_anchor"]:.1e} Å)',
    'the Hessian was computed at the chain-B anchor of the relaxed scan', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.n_imag', FJo['n_imaginary'], None, 0, f'{OS}:n_imaginary', 'imaginary frequencies of the neutral minimum', case='sfa_pair_pcm_d3_oform',
    check=M4['neutral']['n_imag'], tol=0)
add('pcm.oform.lowest', FJo['frequencies_cm-1'][0], 'cm-1', 1, f'{OS}:frequencies_cm-1[0]', 'lowest harmonic frequency of the neutral minimum (inter-ion libration)',
    case='sfa_pair_pcm_d3_oform', check=M4['neutral']['lowest3'][0], tol=0.12)
add('pcm.oform.lowest3', FJo['frequencies_cm-1'][:3], 'cm-1', 1, f'{OS}:frequencies_cm-1[0:3]', 'three lowest frequencies of the neutral minimum', case='sfa_pair_pcm_d3_oform')
add('freq.lowest3', FJs['frequencies_cm-1'][:3], 'cm-1', 1, 'molecules/sfa/results/sfa_pair_pcm_d3/freq.json:frequencies_cm-1[0:3]',
    'three lowest frequencies of the salt minimum (the optimisation was accepted on forces after a soft inter-ion crawl; all real)', case='sfa_pair_pcm_d3')
add('pcm.oform.fd_asymmetry', FJo['fd_hessian_asymmetry'], None, 4, f'{OS}:fd_hessian_asymmetry', 'relative asymmetry of the finite-difference Hessian', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.nu_OH', FJo['frequencies_cm-1'][M4['neutral']['proton_mode']['n'] - 1], 'cm-1', 0, f'{OS}:frequencies_cm-1[{M4["neutral"]["proton_mode"]["n"] - 1}]',
    f"O5–H1A stretch of O5–H1A···N1 (mode {M4['neutral']['proton_mode']['n']} of 75; the strongest IR band of the neutral pair)", case='sfa_pair_pcm_d3_oform',
    check=M4['neutral']['proton_mode']['nu'], tol=0.2)
add('pcm.oform.nu_OH.scaled', FJo['frequencies_cm-1'][M4['neutral']['proton_mode']['n'] - 1] * FACTS['freq.scale_factor']['value'], 'cm-1', 0,
    'pcm.oform.nu_OH × freq.scale_factor', 'the same scaled by 0.9688', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.I_OH', FJo['ir'][M4['neutral']['proton_mode']['n'] - 1]['intensity_km_mol'], 'km mol-1', 0, f'{OS}:ir[{M4["neutral"]["proton_mode"]["n"] - 1}]',
    'IR intensity of the O5–H1A stretch', case='sfa_pair_pcm_d3_oform', check=M4['neutral']['proton_mode']['I_km'], tol=1.0)
add('pcm.oform.strongest_is_OH', M4['neutral']['strongest_IR']['n'] == M4['neutral']['proton_mode']['n'], None, None, f'{D13S}:m4.neutral.strongest_IR',
    'the O5–H1A stretch is the strongest IR band of the neutral minimum (the N1–H1A stretch, 2190 cm⁻¹, is that of the salt)', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.zpe_OH', 0.5 * FJo['frequencies_cm-1'][M4['neutral']['proton_mode']['n'] - 1] * 2.859144e-3, 'kcal mol-1', 2, 'pcm.oform.nu_OH / 2 (h c N_A)',
    '½hν of the O5–H1A stretch of the neutral minimum', case='sfa_pair_pcm_d3_oform', check=M4['half_hnu_OH_kcal'], tol=0.001)
add('pcm.oform.anion_OH', sorted(x['nu'] for x in M4['neutral']['anion_OH']), 'cm-1', 0, f'{D13S}:m4.neutral.anion_OH',
    'O3–H3A and O4–H4A stretches of the neutral minimum (H3PO4 half)', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.ZPE_Ha', FJo['ZPE_Ha'], 'Ha', 6, f'{OS}:ZPE_Ha', 'zero-point energy of the neutral minimum (unscaled)', case='sfa_pair_pcm_d3_oform',
    check=M4['neutral']['thermo']['ZPE_Ha'], tol=1e-5)
add('pcm.oform.ZPE', FJo['ZPE_Ha'] * HA_KCAL, 'kcal mol-1', 2, f'{OS}:ZPE_Ha × 627.5095', 'the same in kcal mol⁻¹', case='sfa_pair_pcm_d3_oform')
add('pcm.oform.H298', FJo['H298_Ha'], 'Ha', 6, f'{OS}:H298_Ha', 'enthalpy of the neutral minimum at 298.15 K', case='sfa_pair_pcm_d3_oform',
    check=M4['neutral']['thermo']['H_Ha'], tol=1e-5)
add('pcm.oform.G298', FJo['G298_Ha'], 'Ha', 6, f'{OS}:G298_Ha', 'Gibbs energy of the neutral minimum at 298.15 K, 1 atm', case='sfa_pair_pcm_d3_oform',
    check=M4['neutral']['thermo']['G_Ha'], tol=1e-5)
add('pcm.oform.S', M4['neutral']['thermo']['S_J_molK'], 'J mol-1 K-1', 1, f'{D13S}:m4.neutral.thermo.S_J_molK', 'entropy of the neutral minimum (own RRHO)',
    case='sfa_pair_pcm_d3_oform', status='own derivation')
# Numbers audit burst 3 round 2 (29 Sep 2026): the thermal corrections are read from the raw freq.json (H298_Ha − E_Ha,
# G298_Ha − E_Ha), the same file as the H and G rows of Table S6, so that the row equals H − E as printed. The d13
# value (other atomic masses; 132.8151 against 132.8128 kcal/mol for H) is kept as the check; it printed 132.82.
add('pcm.oform.Hcorr', (FJo['H298_Ha'] - FJo['E_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{OS}:H298_Ha − E_Ha', 'thermal correction to the enthalpy (incl. ZPE)',
    case='sfa_pair_pcm_d3_oform', check=M4['neutral']['thermo']['Hcorr_Ha'] * HA_KCAL, tol=0.005)
add('pcm.oform.Gcorr', (FJo['G298_Ha'] - FJo['E_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{OS}:G298_Ha − E_Ha', 'thermal correction to the Gibbs energy',
    case='sfa_pair_pcm_d3_oform', check=M4['neutral']['thermo']['Gcorr_Ha'] * HA_KCAL, tol=0.005)
Dk = M4['diff_kcal']
dsrc = f'{OS} minus molecules/sfa/results/sfa_pair_pcm_d3/freq.json'
add('pt.pcm.dE', (FJo['E_Ha'] - FJs['E_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{dsrc} (E_Ha)', 'E(neutral minimum) − E(salt minimum), IEF-PCM (water), electronic',
    case='sfa_pair_pcm_d3', check=Dk['dE'], tol=1e-6)
add('pt.pcm.dZPE', (FJo['ZPE_Ha'] - FJs['ZPE_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{dsrc} (ZPE_Ha)', 'ZPE(neutral) − ZPE(salt): the salt has the larger zero-point energy',
    case='sfa_pair_pcm_d3', check=Dk['dZPE'], tol=0.005)
add('pt.pcm.dEZPE', ((FJo['E_Ha'] + FJo['ZPE_Ha']) - (FJs['E_Ha'] + FJs['ZPE_Ha'])) * HA_KCAL, 'kcal mol-1', 2, f'{dsrc} (E_Ha + ZPE_Ha)',
    'Δ(E + ZPE), neutral − salt', case='sfa_pair_pcm_d3', check=Dk['dEZPE'], tol=0.005)
add('pt.pcm.dH298', (FJo['H298_Ha'] - FJs['H298_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{dsrc} (H298_Ha)', 'ΔH(298.15 K), neutral − salt', case='sfa_pair_pcm_d3',
    check=Dk['dH'], tol=0.005)
add('pt.pcm.dG298', (FJo['G298_Ha'] - FJs['G298_Ha']) * HA_KCAL, 'kcal mol-1', 2, f'{dsrc} (G298_Ha)', 'ΔG(298.15 K, 1 atm), neutral − salt', case='sfa_pair_pcm_d3',
    check=Dk['dG'], tol=0.005)
add('pt.pcm.dS', Dk['dS_J_molK'], 'J mol-1 K-1', 1, f'{D13S}:m4.diff_kcal.dS_J_molK', 'S(neutral) − S(salt)', case='sfa_pair_pcm_d3', status='own derivation')
add('pt.pcm.dG298_kJ', (FJo['G298_Ha'] - FJs['G298_Ha']) * HA_KCAL * 4.184, 'kJ mol-1', 2, 'pt.pcm.dG298 × 4.184', 'the same in kJ mol⁻¹', case='sfa_pair_pcm_d3')
add('pt.pcm.nu_shift', FJo['frequencies_cm-1'][M4['neutral']['proton_mode']['n'] - 1] - FACTS['freq.NH_Hbond']['value'], 'cm-1', 0,
    'pcm.oform.nu_OH − freq.NH_Hbond', 'O5–H1A stretch of the neutral minimum minus the N1–H1A stretch of the salt (both unscaled)', case='sfa_pair_pcm_d3')
# Prose audit burst 3 round 2 (29 Sep 2026; re-review O8): the share of the proton stretch in ΔZPE, so that the Table S6
# note can show that the sign of ΔZPE does not rest on that strongly anharmonic mode. Checked against the proton-mode
# wavenumbers of the independent d13 diagonalisation.
_dzpe_p = FACTS['pcm.oform.zpe_OH']['value'] - FACTS['zpe.NH']['value']
_dzpe_p_chk = 0.5 * (M4['neutral']['proton_mode']['nu'] - M4['salt']['proton_mode']['nu']) * 2.859144e-3
add('pt.pcm.dZPE_proton', _dzpe_p, 'kcal mol-1', 2, 'pcm.oform.zpe_OH − zpe.NH',
    '½hν(O5–H1A, neutral minimum) − ½hν(N1–H1A, salt minimum): the proton stretch alone raises ΔZPE, i.e. favours the salt',
    case='sfa_pair_pcm_d3', check=_dzpe_p_chk, tol=0.005)
add('pt.pcm.dZPE_wo_proton', FACTS['pt.pcm.dZPE']['value'] - _dzpe_p, 'kcal mol-1', 2, 'pt.pcm.dZPE − pt.pcm.dZPE_proton',
    'ΔZPE (neutral − salt) without the proton stretch of either minimum', case='sfa_pair_pcm_d3',
    check=M4['diff_kcal']['dZPE'] - _dzpe_p_chk, tol=0.005)
CAVEATS.append('Neutral IEF-PCM minimum (section 33): Δ(E + ZPE), ΔH and ΔG are harmonic RRHO numbers on a continuum solute, with every '
               'mode down to 17.7 cm⁻¹ treated as a harmonic vibration; they sharpen the electronic ordering (the salt has the larger ZPE) '
               'but carry the continuum error of section 31 (the separated-species exchange is 5.6 kcal/mol off the measured ΔpKa). The '
               'proton stretches (2190 and 2567 cm⁻¹) are harmonic values of strongly anharmonic modes.')

# =====================================================================================================================
section('34. QTAIM of the nine crystal hydrogen bonds in a cluster cut from the periodic cells (referee M6, VM burst 3)',
        'Five ion pairs (135 atoms), the smallest connected set found by the cluster builder ("minimal" cover, build.json) that '
        'contains all nine CIF hydrogen bonds, eight of them at the central pair: N1–H1B···O6 occurs only between outer pairs (C1>A4 '
        'and C4>A2; the six-pair "central" cover with all nine at the central pair was not computed). Cut once from relaxed cell (i) (salt) and once from the constrained single configuration relax_iirfix '
        '(H1A held on O5 at 0.98 Å in every pair), every other atom at its periodic PBE-D3(BJ) position. B3LYP-D3(BJ)/6-31+G(d,p) '
        'single points in the gas phase (PySCF, density fitting, grid 3), Multiwfn QTAIM of the whole cluster. The primary value '
        'of each bond is its instance at the central pair (for N1–H1B···O6 the primary instance C1>A4); where a bond occurs several times the spread over the instances shows '
        'the edge sensitivity. In the O5-protonated cluster the N1–H1A···O5 critical point is that of N1···H1A. Re-derived '
        '(derive/d13_referee_calcs.py): the raw CPprop.txt of both clusters parsed, each contact assigned to the (3,−1) point '
        'nearest its H···A segment by our own geometry, D···A of every contact compared with the cluster build.')
M6 = D13['m6']
QN_ = J(RES / 'crystal_cluster_qtaim/qtaim_nine_bonds.json')
QS = 'molecules/sfa/results/crystal_cluster_qtaim/qtaim_nine_bonds.json'
cs_, cn_ = M6['clusters']['salt'], M6['clusters']['neutral']
add('qtaim.xclu.n_atoms', QN_['n_atoms'], None, 0, f'{QS}:n_atoms', 'atoms of the cluster (five ion pairs)', check=cs_['n_atoms'], tol=0)
add('qtaim.xclu.n_pairs', M6['n_pairs'], None, 0, 'crystal_cluster_qtaim/build.json:structures.salt.n_pairs', 'ion pairs in the cluster')
add('qtaim.xclu.nao', QN_['nao'], None, 0, f'{QS}:nao', 'basis functions of the cluster, 6-31+G(d,p) spherical (351 per ion pair)',
    check=cs_['nao_counted'], tol=0)
add('qtaim.xclu.level', 'B3LYP-D3(BJ)/6-31+G(d,p), gas phase, single points at the periodic PBE-D3(BJ) geometries; PySCF 2.14.0 with density fitting, grid 3; '
    'Multiwfn 2026.9.20 QTAIM', None, None, f'{QS}:level; sfa_xclu_*/scf.json', 'level of the cluster QTAIM')
add('qtaim.xclu.basis', cs_['basis'], None, None, 'sfa_xclu_salt/scf.json:level.basis_print', 'basis of the cluster single points')
add('qtaim.xclu.max_DA_dev', max(max(b['salt']['DA_vs_build'], b['neutral']['DA_vs_build']) for b in M6['bonds']), 'A', 5, f'{D13S}:m6.bonds[*].DA_vs_build',
    'largest difference between D···A from the Multiwfn nuclei and the cluster build (rounding of the printed coordinates)')
add('qtaim.xclu.salt_source', M6['salt_source'], None, None, 'crystal_cluster_qtaim/build.json:structures.salt.source', 'geometry of the salt cluster: relaxed cell (i)')
add('qtaim.xclu.neutral_source', M6['neutral_source'], None, None, 'crystal_cluster_qtaim/build.json:structures.neutral.source',
    'geometry of the O5-protonated cluster: the constrained single configuration (O5–H1A 0.98 Å)')
add('qtaim.xclu.E_salt', cs_['E_Ha'], 'Ha', 6, 'molecules/sfa/results/sfa_xclu_salt/scf.json:e_tot_Ha', 'energy of the salt cluster')
add('qtaim.xclu.E_neutral', cn_['E_Ha'], 'Ha', 6, 'molecules/sfa/results/sfa_xclu_neutral/scf.json:e_tot_Ha', 'energy of the O5-protonated cluster')
add('qtaim.xclu.dE_kJ', QN_['cluster_energy']['neutral']['dE_per_fu_kJ'], 'kJ mol-1', 1, f'{QS}:cluster_energy.neutral.dE_per_fu_kJ',
    'E(O5-protonated cluster) − E(salt cluster) per ion pair, gas phase: negative, the finite cluster prefers the neutral arrangement',
    check=M6['dE_kJ_fu'], tol=1e-6)
add('qtaim.xclu.dE_kcal', QN_['cluster_energy']['neutral']['dE_per_fu_kcal'], 'kcal mol-1', 1, f'{QS}:cluster_energy.neutral.dE_per_fu_kcal', 'the same in kcal mol⁻¹',
    check=M6['dE_kcal_fu'], tol=1e-6)
add('qtaim.xclu.swing_kJ', FACTS['periodic.iirfix.dE_kJ']['value'] - QN_['cluster_energy']['neutral']['dE_per_fu_kJ'], 'kJ mol-1', 0,
    'periodic.iirfix.dE_kJ − qtaim.xclu.dE_kJ', 'difference between the periodic (+57.1, PBE-D3(BJ)) and the cluster energy (−54.0, B3LYP-D3(BJ)) of the same '
    'constrained configuration, per formula unit: what the infinite lattice adds (different functionals and basis; qualitative)', status='qualitative comparison')
add('qtaim.xclu.gap_salt', cs_['gap_eV'], 'eV', 2, 'molecules/sfa/results/sfa_xclu_salt/scf.json:gap_eV',
    'HOMO–LUMO gap of the salt cluster: small, the surface anions and cations of a finite ionic cluster in vacuum sit close in energy')
add('qtaim.xclu.gap_neutral', cn_['gap_eV'], 'eV', 2, 'molecules/sfa/results/sfa_xclu_neutral/scf.json:gap_eV', 'HOMO–LUMO gap of the O5-protonated cluster')
for st, c_ in (('salt', cs_), ('neutral', cn_)):
    add(f'qtaim.xclu.{st}.n_cp', sum(c_['counts'].values()), None, 0, f'molecules/sfa/results/sfa_xclu_{st}/mwfn_extra/CPprop.txt (own count)',
        f'critical points of the {st} cluster', check=QN_['qtaim'][st]['n_cp'], tol=0)
    add(f'qtaim.xclu.{st}.poincare_hopf', c_['poincare_hopf'] == 1, None, None, f'{D13S}:m6.clusters.{st}.poincare_hopf', f'Poincaré–Hopf relation satisfied ({st})')
for b in M6['bonds']:
    nm = b['bond'].replace('-', '').replace('...', '_')
    qb = [x for x in QN_['bonds'] if x['bond'] == b['bond']][0]
    lab = b['bond'].replace('-', '–').replace('...', '···')
    for st in ('salt', 'neutral'):
        qs = qb[st]['bcp']
        what = f'{lab}, {"salt" if st == "salt" else "O5-protonated"} cluster ({"primary instance, between outer pairs" if b["bond"] == "N1-H1B...O6" else "central instance"} {b["instance"]})' + (
            ': critical point of N1···H1A, the proton sits on O5' if (st == 'neutral' and b['neutral']['reverse']) else '')
        add(f'qtaim.xclu.{nm}.{st}.rho', qs['rho'], 'a.u.', 4, f'{QS}:bonds[{b["bond"]}].{st}.bcp.rho', f'ρ at the bond critical point: {what}',
            check=b[st]['rho'], tol=1e-7)
        add(f'qtaim.xclu.{nm}.{st}.lap', qs['lap'], 'a.u.', 4, f'{QS}:bonds[{b["bond"]}].{st}.bcp.lap', f'∇²ρ: {what}', check=b[st]['lap'], tol=1e-7)
        add(f'qtaim.xclu.{nm}.{st}.H', qs['H'], 'a.u.', 5, f'{QS}:bonds[{b["bond"]}].{st}.bcp.H', f'total energy density H: {what}', check=b[st]['H'], tol=1e-7)
        add(f'qtaim.xclu.{nm}.{st}.HA', qb[st]['geometry']['H..A'], 'A', 3, f'{QS}:bonds[{b["bond"]}].{st}.geometry["H..A"]', f'H···A of {lab} ({st} cluster)')
    add(f'qtaim.xclu.{nm}.d_rho_pct', qb['delta_neutral_minus_salt']['rho_pct'], '%', 0, f'{QS}:bonds[{b["bond"]}].delta_neutral_minus_salt.rho_pct',
        f'change of ρ at {lab} when H1A sits on O5, % of the salt value', check=b['d_rho_pct'], tol=1e-4)
    add(f'qtaim.xclu.{nm}.d_rho', qb['delta_neutral_minus_salt']['rho'], 'a.u.', 4, f'{QS}:bonds[{b["bond"]}].delta_neutral_minus_salt.rho',
        f'change of ρ at {lab}, O5-protonated − salt', check=b['d_rho'], tol=1e-7)
    add(f'qtaim.xclu.{nm}.d_H', qb['delta_neutral_minus_salt']['H'], 'a.u.', 5, f'{QS}:bonds[{b["bond"]}].delta_neutral_minus_salt.H',
        f'change of H at {lab}, O5-protonated − salt', check=b['d_H'], tol=1e-7)
    if b['n_instances'] > 1:
        add(f'qtaim.xclu.{nm}.rho_range_salt', b['rho_range_salt'], 'a.u.', 4, f'{D13S}:m6.bonds[{b["bond"]}].rho_range_salt',
            f'ρ over the {b["n_instances"]} instances of {lab} in the salt cluster (edge sensitivity)')
        add(f'qtaim.xclu.{nm}.rho_range_neutral', b['rho_range_neutral'], 'a.u.', 4, f'{D13S}:m6.bonds[{b["bond"]}].rho_range_neutral',
            f'the same in the O5-protonated cluster')
    add(f'qtaim.xclu.{nm}.n_instances', b['n_instances'], None, 0, f'{D13S}:m6.bonds[{b["bond"]}].n_instances', f'instances of {lab} in the cluster')
CAVEATS.append('Cluster QTAIM (section 34): a finite, gas-phase cluster of five ion pairs at the periodic geometries, a double-zeta basis and '
               'no embedding. Its descriptors show how each of the nine bonds responds when H1A sits on O5 (ρ at N1–H1B···O6, N1–H1C···O6 '
               'and O4–H4A···O5 falls by 34, 29 and 24 %); its energy does not describe the crystal. The salt cluster has a HOMO–LUMO gap of '
               'only 1.15 eV and prefers the neutral arrangement by 54.0 kJ/mol per ion pair, the opposite sign to the periodic cell (+57.1 '
               'kJ/mol for the same constrained configuration at PBE-D3(BJ)): the under-coordinated surface ions of a finite cluster in vacuum '
               'favour neutral partners, so the local hydrogen-bond network alone, without the rest of the ionic lattice, does not hold the '
               'proton. The O5-protonated cluster uses the constrained single configuration (O5–H1A held at 0.98 Å), not the profile.')

# =====================================================================================================================
section('35. Relaxed C3–C4–S1–N2 torsion scans (referee M9, VM burst 3)',
        'B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM (water), grid 4: C3–C4–S1–N2 held at each value, everything else optimised (geomeTRIC), '
        'from 90° to 180° in 10° steps plus the torsions of the crystal ligands and of docked poses; sulfanilamide 22 points (one on '
        'the C5 branch as a mirror check), the sulfonamidate 18. Strain = E(constrained) − E(unconstrained minimum), the minimum '
        're-optimised from the IEF-PCM optimum with the same settings (ref.json), electronic. Values off the computed points come '
        'from a cubic spline through the converged points and the minimum. This replaces the rigid single points of sections 29 '
        'and 31 (section 36). Re-derived (derive/d13_referee_calcs.py): every point energy re-read, τ_fold re-measured from each '
        'point geometry, the spline redone and evaluated at every docked seed of every production job.')
M9 = D13['m9']
TS_ = J(RES / 'torsion_scan/torsion_strain.json')
TSS = 'molecules/sfa/results/torsion_scan/torsion_strain.json'
SPN = {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate'}
add('tors.scan.level', 'B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM (water), grid 4, relaxed (C3–C4–S1–N2 constrained, all else optimised)', None, None,
    f'{TSS}:per_case.*.level', 'level of the relaxed torsion scans')
for sp, case in (('sfa_neutral', 'sulfanilamide_pcm_d3'), ('sfa_anion', 'sulfonamidate_pcm_d3')):
    m = M9[sp]
    pc = TS_['per_case'][case]
    lab = SPN[sp]
    add(f'tors.ref.{sp}.E', pc['reference']['E_Ha'], 'Ha', 6, f'{TSS}:per_case.{case}.reference.E_Ha', f'{lab}: energy of the unconstrained IEF-PCM minimum (scan reference)',
        check=m['E_ref'], tol=1e-9)
    add(f'tors.ref.{sp}.tau', pc['reference']['tau_fold_deg'], 'deg', 1, f'{TSS}:per_case.{case}.reference.tau_fold_deg', f'{lab}: τ_fold of the scan reference minimum',
        check=m['tau_ref'], tol=1e-4)
    add(f'tors.ref.{sp}.ang_CSN', m['ang_CSN_ref'], 'deg', 1, f'{D13S}:m9.{sp}.ang_CSN_ref', f'{lab}: C4–S1–N2 angle of the scan reference')
    add(f'tors.scan.{sp}.n_points', pc['n_points'], None, 0, f'{TSS}:per_case.{case}.n_points', f'{lab}: constrained points', check=m['n_points'], tol=0)
    add(f'tors.scan.{sp}.n_converged', pc['n_converged'], None, 0, f'{TSS}:per_case.{case}.n_converged', f'{lab}: converged constrained points',
        check=m['n_converged'], tol=0)
    add(f'tors.scan.{sp}.interp_check', pc['interp_check_max_abs_kcal'], 'kcal mol-1', 3, f'{TSS}:per_case.{case}.interp_check_max_abs_kcal',
        f'{lab}: largest difference between the 10° grid spline and the exact off-grid points')
    for pt in sorted([q for q in m['points'] if q['branch'] == 'C3'], key=lambda q: q['target']):
        kk = f"t{int(round(pt['target'] * 10)):04d}"
        prof = [q for q in TS_['profiles'][case] if abs(q['tau_C3_deg'] - pt['target']) < 1e-6]
        add(f'tors.scan.{sp}.{kk}', prof[0]['strain_kcal'] if prof else pt['strain'], 'kcal mol-1', 2,
            f'{TSS}:profiles.{case}[τ = {pt["target"]:g}]', f'{lab}: relaxed strain with C3–C4–S1–N2 held at {pt["target"]:g}°',
            check=pt['strain'], tol=1e-6)
        add(f'tors.scan.{sp}.{kk}.tau', pt['tau_fold'], 'deg', 1, f'{D13S}:m9.{sp}.points ({pt["file"]}, re-measured)', f'{lab}: τ_fold of that point')
    add(f'tors.scan.{sp}.max', m['strain_max'], 'kcal mol-1', 2, f'{D13S}:m9.{sp}.strain_max', f'{lab}: largest relaxed strain between 90° and 180° '
        f'(at {m["tau_at_max"]:.0f}°)')
    add(f'tors.scan.{sp}.min', m['strain_min'], 'kcal mol-1', 2, f'{D13S}:m9.{sp}.strain_min', f'{lab}: smallest relaxed strain between 90° and 180° '
        f'(at {m["tau_at_min"]:.0f}°; negative: below the re-optimised minimum)')
    add(f'tors.scan.{sp}.tau_at_max', m['tau_at_max'], 'deg', 0, f'{D13S}:m9.{sp}.tau_at_max', f'{lab}: τ_fold of the largest strain')
    add(f'tors.scan.{sp}.tau_at_min', m['tau_at_min'], 'deg', 0, f'{D13S}:m9.{sp}.tau_at_min', f'{lab}: τ_fold of the smallest strain')
    # labelled rows: SAPH crystal, crystal ligands, representative poses and job means
    for r_ in [x for x in TS_['table'] if x['case'] == case and x['kind'] != 'grid']:
        L_ = r_['label']
        if r_['kind'] == 'SAPH crystal':
            kk = 'saph_crystal'
        elif r_['kind'] == 'co-crystal':
            kk = 'cocrystal.' + L_[:4]
        else:
            mm = re.match(r'^(\w{4}) \([^)]*\), (Vina|Vinardo|AD4Zn), (representative|mean)', L_)
            kk = ('rep.' if mm.group(3) == 'representative' else 'mean.') + mm.group(1) + '.' + mm.group(2).lower()
        der = m['labelled'][L_]
        add(f'tors.strain.{sp}.{kk}', r_['strain_kcal'], 'kcal mol-1', 2, f'{TSS}:table["{L_}"].strain_kcal ({r_["how"]})',
            f'{lab}: relaxed DFT strain at τ_fold = {r_["tau_fold_deg"]:.1f}° ({L_})', check=der['strain_spline'], tol=0.03)
    # every docked seed: strain range per job (own spline)
    for key_, jb in m['jobs'].items():
        tg_, rest_ = key_.split('/')
        fn_ = rest_.split('__')[1]
        add(f'tors.seeds.{sp}.{tg_.split("_")[1]}.{fn_}', [jb['strain_min'], jb['strain_max']], 'kcal mol-1', 2, f'{D13S}:m9.{sp}.jobs["{key_}"]',
            f'{lab}, {tg_.split("_")[1]} {fn_}: relaxed DFT strain at the τ_fold of the {jb["n"]} docked top poses (min, max)',
            status='own derivation (spline of the relaxed scan at each seed)')
        add(f'tors.seeds.{sp}.{tg_.split("_")[1]}.{fn_}.tau_mean', jb['tau_mean'], 'deg', 1, f'{D13S}:m9.{sp}.jobs["{key_}"].tau_mean',
            f'{lab}, {tg_.split("_")[1]} {fn_}: mean τ_fold of the docked top poses')
m9n, m9a = M9['sfa_neutral'], M9['sfa_anion']
if m9n['mirror']:
    add('tors.scan.sfa_neutral.mirror_diff', abs(m9n['mirror'][0]['C5'] - m9n['mirror'][0]['C3']), 'kcal mol-1', 3, f'{D13S}:m9.sfa_neutral.mirror',
        'sulfanilamide at 150°: C5-branch point minus the C3-branch point (ring-flip check)')
add('tors.scan.sfa_anion.span', m9a['strain_max'] - m9a['strain_min'], 'kcal mol-1', 2, f'{D13S}:m9.sfa_anion',
    'sulfonamidate: spread of the relaxed profile between 90° and 180° (flat within 0.1 kcal/mol)')
reps_n = [v_ for v_ in m9n['labelled'].values() if v_['kind'] == 'docked pose']
reps_a = [v_ for v_ in m9a['labelled'].values() if v_['kind'] == 'docked pose']
add('tors.mmff_minus_dft.sfa_neutral', [min(v_['mmff'] - v_['strain_src'] for v_ in reps_n), max(v_['mmff'] - v_['strain_src'] for v_ in reps_n)],
    'kcal mol-1', 2, f'{TSS}:table (mmff_relaxed_kcal − strain_kcal, representative poses)',
    'sulfanilamide: MMFF94s relaxed strain minus relaxed DFT strain at the ten representative poses (range)')
add('tors.mmff_minus_dft.sfa_anion', [min(v_['mmff'] - v_['strain_src'] for v_ in reps_a), max(v_['mmff'] - v_['strain_src'] for v_ in reps_a)],
    'kcal mol-1', 2, f'{TSS}:table (mmff_relaxed_kcal − strain_kcal, representative poses)',
    'sulfonamidate: MMFF94s relaxed strain minus relaxed DFT strain at the eight representative poses (range): MMFF94s overstates the anion')
g6 = FACTS['tors.strain.sfa_neutral.cocrystal.6RL9']['value']
add('tors.overshoot_cost.6RL9', FACTS['tors.strain.sfa_neutral.mean.6RL9.ad4zn']['value'] - g6, 'kcal mol-1', 2,
    'tors.strain.sfa_neutral.mean.6RL9.ad4zn − tors.strain.sfa_neutral.cocrystal.6RL9',
    'sulfanilamide: relaxed strain at the mean 6RL9 AutoDock4Zn torsion (147.0°) minus that at the 6RL9 crystal-ligand torsion (132.1°)')
add('tors.cocrystal_strain_range', [FACTS['tors.strain.sfa_neutral.cocrystal.3TZF']['value'], g6], 'kcal mol-1', 2,
    'tors.strain.sfa_neutral.cocrystal.3TZF / .6RL9', 'sulfanilamide: relaxed strain at the torsions of the three crystal ligands (105.6–132.1°)')
CAVEATS.append('Relaxed torsion scans (section 35): electronic energies in IEF-PCM of the free ligand, one-dimensional in C3–C4–S1–N2 (the '
               'N–H rotors and the C4–S1–N2 angle relax at each point). They replace the rigid single points of sections 29 and 31. The '
               'sulfonamidate minimum re-optimised at 101.6° is not the lowest point of its own scan: the profile is flat to within 0.1 '
               'kcal/mol, its lowest point at 170° (−0.09). Strain values between computed points are spline interpolations (checked '
               'against exact points to 0.023 kcal/mol).')
MISSING[:] = [m_ for m_ in MISSING if not m_.startswith('No periodic proton-transfer profile') and not m_.startswith('M8 (pair)2 cluster')
              and not m_.startswith('dock.strainref.sfa_cation')]
MISSING.append('Periodic profile (section 32): no cell-relaxed end point (vc-relax) and no hybrid relaxation; the PBE0-D3(BJ) curve is single '
               'points at PBE-D3(BJ) geometries. Stated limitation (user decision, 29 Sep 2026; VM budget closed at 30.0 of 40 VM-h).')
MISSING.append('No molecular calibration of the PBE-to-PBE0 shift (PBE-D3(BJ) single points on the three IEF-PCM scan geometries, referee M3.2): '
               'stated limitation (user decision, 29 Sep 2026).')
MISSING.append('Referee M5 not run: no PCM Hessians of sulfanilamide and H3PO4, no SMD or radius-scaling single points, no larger-basis single '
               'points of the scan geometries; the continuum error (pt.sep.continuum_error) is stated, not corrected (user decision).')
MISSING.append('Referee M1 controls not run: no point-charge embedding of the ion pair and no cluster-continuum pair; the five-pair gas-phase '
               'cluster of section 34 is analysed for its bond descriptors, and only the sign of its energy is read (a finite cluster in '
               'vacuum, not an energy control of the crystal).')
MISSING[:] = [('Frequencies: the IEF-PCM ion pair (salt and neutral minima), the gas-phase sulfanilamide, cation (gas/PCM), H2PO4− (gas/PCM), '
               'H3PO4 (gas) and urea have Hessians. No frequencies for the gas-phase pair, the plain-B3LYP and def2-TZVP pairs, H3PO4 (PCM), and '
               'sulfanilamide (PCM) and the sulfonamidate (PCM): they are stationary points on the force criterion only (SI Section S2).') if m_.startswith('Frequencies:') else m_ for m_ in MISSING]

# =====================================================================================================================
section('36. Superseded facts (29 Sep 2026, after VM burst 3)',
        '(a) Relabelled, still citable as what they are: the energies and geometry of relax_iirfix (sections 17 and 28, and the hybrid '
        'correction of section 28) describe one constrained configuration, O5–H1A held at 0.98 Å with H1A and O5 fixed; they are no '
        'longer the measure of the periodic cost, which is the constrained N1–H1A profile of section 32. (b) Retired, never cite: the '
        'rigid B3LYP-D3(BJ) single points of the docked poses (section 29) and their same-H references (section 31), replaced by the '
        'relaxed torsion scans of section 35 (the rigid values were dominated by the placement of the N–H hydrogens, referee M9); and '
        'the order-of-magnitude estimates on relax_iirfix of section 31 (O–H relief, stress, their differences), replaced by the profile '
        'itself and by the stress estimate of section 32. Retired facts keep their number in superseded_value; value is null, so '
        'numbers_audit.py cannot match a printed number to them.')
RELABEL = [k_ for k_ in FACTS if (k_.startswith('periodic.iirfix.') or k_.startswith('proton.periodic.iirfix.') or k_.startswith('periodic.hybrid.e100.')
                                   or k_.startswith('periodic.hybrid.e80.') or k_.startswith('periodic.hybrid.nq211.')
                                   or k_ in ('periodic.hybrid.corr_spread_kJ', 'periodic.hybrid.verdict', 'periodic.hybrid.PBE_cutoff_change'))]
RETIRE = ([k_ for k_ in FACTS if k_.startswith('dock.strain.') and not (k_.endswith('.tau') or k_.endswith('.seed'))]
          + [k_ for k_ in FACTS if k_.startswith('dock.strainref.')]
          + [f'periodic.iirfix.{s_}' for s_ in ('dx_lo', 'dx_hi', 'OH_relief_lo', 'OH_relief_hi', 'best_minus_OH_hi', 'best_minus_OH_hi_kcal',
                                                 'stress_B10', 'stress_B20', 'best_minus_all', 'best_minus_all_kcal')])
RELABEL = [k_ for k_ in RELABEL if k_ not in RETIRE]
for k_ in RELABEL:
    f_ = FACTS[k_]
    f_['status'] = ((f_.get('status') + '; ') if f_.get('status') else '') + (
        'constrained single configuration (O5–H1A held at 0.98 Å, H1A and O5 fixed): superseded as the periodic cost by the constrained '
        'N1–H1A profile (section 32); cite only as this one configuration')
    if not f_['meaning'].startswith('constrained single configuration'):
        f_['meaning'] = 'constrained single configuration: ' + f_['meaning']
for k_ in RETIRE:
    f_ = FACTS[k_]
    f_['superseded_value'] = f_['value']
    f_['value'] = None
    f_['superseded'] = True
    f_['status'] = ('SUPERSEDED (never cite): ' + ('rigid DFT single point of a docked pose, dominated by the MMFF94s placement of the N–H hydrogens; '
                                                    'replaced by the relaxed torsion scans (section 35)' if k_.startswith('dock.strain')
                                                    else 'order-of-magnitude check on the constrained single configuration; replaced by the constrained '
                                                         'profile and its stress estimate (section 32)'))
FACTS['dock.tors.strain_max']['meaning'] = ('largest MMFF94s relaxed strain of those poses (the relaxed DFT strain of the same torsions is in '
                                            'section 35)')
add('superseded.relabelled', sorted(RELABEL), None, None, 'facts_referee.py section 36 (a)', 'facts relabelled as the constrained single configuration (still citable as such)')
add('superseded.retired', sorted(RETIRE), None, None, 'facts_referee.py section 36 (b)', 'facts retired in favour of sections 32 and 35 (never cite)')
CAVEATS[:] = [c_ for c_ in CAVEATS if not c_.startswith('DFT strain of the docked poses (section 29)')]
CAVEATS.append('SUPERSEDED (section 36): the rigid DFT strain of the docked poses (sections 29, 31) is retired; the relax_iirfix numbers '
               f"(+{FACTS['periodic.iirfix.dE_kJ']['value']:.1f} PBE-D3(BJ), +{FACTS['periodic.hybrid.e100.best_kJ']['value']:.1f} kJ/mol "
               'PBE0-D3(BJ)) are the energy of one constrained single configuration with O5–H1A held at 0.98 Å, '
               f"{FACTS['periodic.profile.iirfix_excess']['value']:.1f} kJ/mol per formula unit above the relaxed profile at the same N1–H1A "
               f"distance; the periodic cost is now the profile of section 32 (+{FACTS['periodic.profile.p170.PBE_kJ']['value']:.1f} PBE-D3(BJ), "
               f"+{FACTS['periodic.profile.p170.PBE0_kJ']['value']:.1f} kJ/mol PBE0-D3(BJ) at 1.70 Å, rising at every point).")
DISAGREE.append('SETTLED (VM burst 3, 29 Sep 2026): referee M2 was right that the single constrained configuration overstated the cost. With '
                'd(N1–H1A) held and O5–H1A free, the profile at the same N1–H1A distance (1.616 Å) lies 10.4 kJ/mol per formula unit below it '
                '(Hermite interpolation of the 1.55 and 1.70 Å points). The sign does not change: the energy rises at every profile point with '
                'PBE-D3(BJ) and with PBE0-D3(BJ), every constraint force points back to N1, and a transfer in one formula unit alone relaxes back '
                'to the salt.')

# rounding ties among the new facts: a stored tie is replaced by its full-precision re-derivation, as in section 30
for _k in list(FACTS)[NEW_KEYS_START:]:
    _f = FACTS[_k]
    _v, _p, _rd = _f.get('value'), _f.get('precision'), _f.get('rederived')
    if not (isinstance(_v, float) and isinstance(_p, int) and _p >= 0 and isinstance(_rd, (int, float)) and not isinstance(_rd, bool)):
        continue
    _t = _Dec(repr(_v)).as_tuple()
    if -_t.exponent != _p + 1 or _t.digits[-1] != 5 or abs(_rd - _v) > 0.5 * 10 ** (-(_p + 1)) + 1e-12:
        continue
    _f['stored_value'] = _v
    _f['value'] = float(_rd)
    _f['note'] = ((_f.get('note') + '; ') if _f.get('note') else '') + f'rounding tie: stored {_v!r}, value = full-precision re-derivation {_rd!r}'
    print('tie replaced:', _k, _v, '->', _rd)
