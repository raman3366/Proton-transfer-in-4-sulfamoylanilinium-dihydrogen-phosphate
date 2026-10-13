"""Facts of burst 2 (folded in 29 Sep 2026; sections 28-30 of facts.md), executed by build_facts.py after
facts_si_extra.py (exec in its namespace, so add(), section(), FACTS, J, r, RES, DOCK, DER, D3, D8, CAVEATS, MISSING and
DISAGREE are the ones defined there).  Kept in its own file so the additions can be reviewed as one block:

  28. Periodic constrained neutral (relax_iirfix: H1A placed 0.98 A from O5 in relaxed cell (i), then H1A and O5 held in
      space) and the PBE0-D3(BJ) hybrid correction: periodic_qe.json ['comparison']['constrained_neutral_iirfix'] and
      ['hybrid'].  Every energy re-derived from the raw Quantum ESPRESSO outputs in the archive (derive/d11_burst2.py);
      the iirfix cell re-analysed with our own geometry code (derive/d08_periodic.py, d03_geometry.py).
  29. DFT strain of the docked poses: B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM single points of the cluster representatives
      (<case>/dock_strain.json of sulfanilamide, the sulfonamidate and the cation), re-derived in d11.
  30. Final passes over all facts: (a) the relax_iifix facts are withdrawn (kept for provenance, never cited): their
      value moves to withdrawn_value and value becomes null, so numbers_audit.py cannot match a printed number to them;
      (b) the rounding rule, round half to even on the full-precision value, applied to the stored rounding ties.
"""
from decimal import Decimal as _Dec, ROUND_HALF_EVEN as _RHE  # noqa: E402

D11 = J(DER / 'd11_burst2.json')
PQJ = J(RES / 'periodic_qe/periodic_qe.json')
RY_KJ_FU = 13.605693 * 96.485332 / 2          # Ry per cell -> kJ/mol per formula unit (Z = 2)
RY_MEV = 13.605693e3
RYB_EVA = 13.605693 / 0.529177210903          # Ry/Bohr -> eV/Å

# =====================================================================================================================
section('28. Periodic constrained neutral (iirfix) and the PBE0-D3(BJ) hybrid correction (burst 2)',
        'relax_iirfix: relaxed cell (i) with H1A of both formula units moved to 0.98 Å from O5 on the N1→O5 line, then '
        'relaxed (PBE-D3(BJ), SSSP, 70/560 Ry, 2×1×1) with H1A and O5 held in space, which keeps the O5–H1A bond and pins the '
        'floating origin of Pc. A constrained structure, not a minimum. Hybrid correction: PBE0-D3(BJ) minus PBE-D3(BJ) single '
        'points at the PBE-D3(BJ) geometries with SG15 ONCV norm-conserving potentials (100/400 Ry, q mesh 1×1×1, primary; '
        'variants 80/320 Ry and q mesh 2×1×1), added to the SSSP PBE-D3(BJ) difference. Source: periodic_qe.json; every '
        'energy re-derived from the archived pw.x outputs (derive/d11_burst2.py: relax "Final energy", PBE "!" energy, PBE0 '
        'first "!" and final "!!" energies), the held atoms and final forces read from relax_iirfix.in/.out, and the single-point '
        'inputs checked to carry the relaxed coordinates (16 inputs, max deviation 1e-8 Å).')
PS = 'molecules/sfa/results/periodic_qe/periodic_qe.json'
cn = PQJ['comparison']['constrained_neutral_iirfix']
dp = D11['periodic']
Ei, Er = PQJ['relax']['i']['energy_Ry'], PQJ['relax']['iirfix']['energy_Ry']
add('periodic.iirfix.E', Er, 'Ry', 7, f'{PS}:relax.iirfix.energy_Ry', 'relaxed energy of the constrained neutral cell (iirfix)', case='periodic_qe',
    check=dp['E_relax_Ry']['iirfix'], tol=1e-8)
add('periodic.iirfix.dE_meV', (Er - Ei) * RY_MEV, 'meV', 1, f'{PS}:relax.iirfix.energy_Ry − relax.i.energy_Ry (JSON comparison.constrained_neutral_iirfix: {cn["E_iirfix_minus_i_meV_cell"]} meV/cell)',
    'E(iirfix) − E(i) per cell, PBE-D3(BJ)/SSSP', case='periodic_qe', check=dp['dE_sssp_meV_cell']['iirfix'], tol=0.01)
add('periodic.iirfix.dE_kJ', (Er - Ei) * RY_KJ_FU, 'kJ mol-1', 1, f'{PS}:relax.iirfix.energy_Ry − relax.i.energy_Ry (JSON: {cn["E_iirfix_minus_i_kJ_mol_per_formula_unit"]} kJ/mol per f.u.)',
    'constrained neutral − salt, PBE-D3(BJ)/SSSP, per formula unit', case='periodic_qe', check=dp['dE_sssp_kJ_fu']['iirfix'], tol=1e-3)
add('periodic.iirfix.dE_kcal', (Er - Ei) * RY_KJ_FU / 4.184, 'kcal mol-1', 1, 'periodic.iirfix.dE_kJ / 4.184',
    'the same in kcal mol⁻¹ per formula unit', case='periodic_qe', check=dp['dE_sssp_kcal_fu']['iirfix'], tol=1e-3)
add('periodic.iirfix.constraint_held', bool(cn['constraint_held']) and D11['periodic']['iirfix_constraint']['held_atoms_max_move_A'] < 1e-6, None, None,
    f'{PS}:comparison.constrained_neutral_iirfix.constraint_held; derived/d11_burst2.json:periodic.iirfix_constraint (held atoms '
    f'{", ".join(dp["iirfix_constraint"]["held_atoms"])} at their input positions)', 'H1A and O5 stayed where they were held; the proton stayed on O5', case='periodic_qe')
rr_ = D8['structures']['relax_iirfix']
add('periodic.iirfix.shift', rr_['mean_heavy_shift_vs_relax_i_norm'], 'A', 2, 'derived/d08_periodic.json:structures.relax_iirfix.mean_heavy_shift_vs_relax_i',
    'mean heavy-atom translation of iirfix relative to relaxed cell (i): no slide of the lattice (withdrawn iifix: 0.66 Å)', case='periodic_qe')
add('periodic.iirfix.N1_moved', rr_['N1_moved_vs_relax_i_A'], 'A', 2, 'derived/d08_periodic.json:structures.relax_iirfix.N1_moved_vs_relax_i_A',
    'displacement of N1 from its position in relaxed cell (i)', case='periodic_qe')
add('periodic.iirfix.rmsd_vs_i', rr_['rmsd_heavy_vs_relax_i_shift_removed'], 'A', 2, 'derived/d08_periodic.json:structures.relax_iirfix.rmsd_heavy_vs_relax_i_shift_removed',
    'heavy-atom RMSD of iirfix against relaxed cell (i), mean shift removed', case='periodic_qe')
add('periodic.iirfix.rmsd_heavy', rr_['rmsd_heavy_shift_removed'], 'A', 3, 'derived/d08_periodic.json:structures.relax_iirfix.rmsd_heavy_shift_removed',
    'heavy-atom RMSD of iirfix against the X-ray, mean shift removed', case='periodic_qe', check=PQJ['relax']['iirfix']['rmsd_heavy_vs_xray']['rmsd_A'], tol=0.0006)
add('periodic.iirfix.ang_N1H1AO5', rr_['fu1']['ang'], 'deg', 1, 'derived/d08_periodic.json:structures.relax_iirfix.fu1.ang',
    'O5–H1A···N1 angle in iirfix', case='periodic_qe', check=cn['proton']['angle_N1_H1A_O5_deg'], tol=0.01)
for h in rr_['hbonds']:
    nm = h['name'].replace('-', '').replace('...', '_')
    add(f'periodic.iirfix.hb.{nm}.DA', h['DA'], 'A', 3, 'derived/d08_periodic.json:structures.relax_iirfix.hbonds',
        f"iirfix: D···A of {h['name']} ({h['sym']}); X-ray {h['cif_DA']:.3f} Å", case='periodic_qe',
        check=[g['D_A_A'] for g in PQJ['relax']['iirfix']['geometry']['hbonds'] if g['bond'].startswith(h['name'])][0], tol=0.0006)
add('periodic.iirfix.bfgs_steps', PQJ['relax']['iirfix']['bfgs_steps'][1], None, 0, f'{PS}:relax.iirfix.bfgs_steps', 'BFGS steps of the iirfix relaxation (64 SCF cycles)', case='periodic_qe')
for k_ in ('i', 'ii'):
    add(f'periodic.bfgs_steps_{k_}', PQJ['relax'][k_]['bfgs_steps'][1], None, 0, f'{PS}:relax.{k_}.bfgs_steps[1] ("bfgs converged in n scf cycles and m bfgs steps")',
        f'BFGS steps of relaxation ({k_}) ({PQJ["relax"][k_]["bfgs_steps"][0]} SCF cycles)', case='periodic_qe')
add('periodic.iirfix.pressure', PQJ['relax']['iirfix']['pressure_kbar_final'], 'kbar', 1, f'{PS}:relax.iirfix.pressure_kbar_final', 'residual pressure of iirfix at the fixed X-ray cell', case='periodic_qe')
ic = dp['iirfix_constraint']
add('periodic.iirfix.max_force_free', ic['max_force_component_free_Ry_bohr'], 'Ry/bohr', 5, 'derived/d11_burst2.json:periodic.iirfix_constraint (last force block of relax_iirfix.out)',
    'largest force component on the 50 free atoms at the end of iirfix (threshold 5e-4)', case='periodic_qe')
add('periodic.iirfix.F_H1A', ic['F_H1A']['fu1']['F_H1A_norm_Ry_bohr'], 'Ry/bohr', 3, 'derived/d11_burst2.json:periodic.iirfix_constraint.F_H1A (relax_iirfix.out, final forces)',
    'force on the held H1A at the end of iirfix (norm; the same in both formula units)', case='periodic_qe')
add('periodic.iirfix.max_force', ic['max_force_component_held_Ry_bohr'], 'Ry/bohr', 5, 'derived/d11_burst2.json:periodic.iirfix_constraint (largest force component, on a held atom)',
    'largest force component in iirfix, on the held H1A/O5 (periodic_qe.json reports the same)', case='periodic_qe', check=PQJ['relax']['iirfix']['max_force_Ry_bohr'], tol=1e-6)
add('periodic.iirfix.F_H1A_toN1', ic['F_H1A']['fu1']['F_H1A_along_H1A_to_N1_Ry_bohr'], 'Ry/bohr', 3, 'derived/d11_burst2.json:periodic.iirfix_constraint.F_H1A.fu1',
    'component of that force along H1A→N1: positive, the held proton is pushed back towards N1', case='periodic_qe')
add('periodic.iirfix.F_H1A_toN1_eVA', ic['F_H1A']['fu1']['F_H1A_along_H1A_to_N1_Ry_bohr'] * RYB_EVA, 'eV/A', 1, 'periodic.iirfix.F_H1A_toN1 × 25.711 (Ry/Bohr → eV/Å)',
    'the same in eV Å⁻¹', case='periodic_qe')
# ---- hybrid correction
HY = PQJ['hybrid']
add('periodic.hybrid.method', 'PBE0 (25 % exact exchange, ACE, Gygi–Baldereschi divergence treatment with x_gamma extrapolation) + two-body D3(BJ) with the pbe0 parameters; SG15 ONCV v1.2 norm-conserving PBE potentials (no nonlinear core correction), ecutrho = ecutfock = 4 ecutwfc; SCF single points at the PBE-D3(BJ) geometries (fixed X-ray cell), k grid 2×1×1; EXX converged when dexx < conv_thr (1e-6 Ry)',
    None, None, f'{PS}:hybrid.method', 'periodic hybrid protocol')
add('periodic.hybrid.ecut', [100, 400], 'Ry', 0, f'{PS}:hybrid.points.sg15_e100.*.ecut_Ry', 'wavefunction/density (= Fock) cutoffs of the primary hybrid variant (SG15)')
add('periodic.hybrid.ecut_alt', [80, 320], 'Ry', 0, f'{PS}:hybrid.points.sg15_e80.*.ecut_Ry', 'cutoffs of the lower-cutoff variant')
add('periodic.hybrid.nq_alt', [2, 1, 1], None, 0, f'{PS}:hybrid.points.sg15_e100_nq211.*.nq', 'EXX q mesh of the denser-q variant (primary 1×1×1)')
add('periodic.hybrid.conv_thr', 1e-6, 'Ry', 6, f'{PS}:hybrid.points.*.*.conv_thr_Ry', 'SCF and EXX convergence threshold of the hybrid single points')
add('periodic.hybrid.exx_converged', all(x['exx_converged'] for v in dp['hybrid'].values() for x in v.values()), None, None,
    'derived/d11_burst2.json:periodic.hybrid (every PBE0 output: "EXX self-consistency reached")', 'all hybrid single points converged in EXX')
VARS = [('e100', 'sg15_e100', '100/400 Ry, q 1×1×1 (primary)'), ('e80', 'sg15_e80', '80/320 Ry, q 1×1×1'),
        ('nq211', 'sg15_e100_nq211', '100/400 Ry, q 2×1×1')]
for vk, vj, vlab in VARS:
    c = HY['comparisons'][vj]['E_iirfix_minus_i']
    d = dp['hybrid'][vj]['iirfix_minus_i']
    b = f'periodic.hybrid.{vk}'
    src = f'{PS}:hybrid.comparisons.{vj}.E_iirfix_minus_i'
    add(b + '.PBE0_meV', c['PBE0_D3BJ_meV_cell'], 'meV', 1, src + '.PBE0_D3BJ_meV_cell', f'E(iirfix) − E(i), PBE0-D3(BJ)/SG15 {vlab}, per cell',
        case='periodic_qe', check=d['PBE0_D3_meV_cell'], tol=0.01)
    add(b + '.PBE_meV', c['PBE_D3BJ_same_setting_meV_cell'], 'meV', 1, src + '.PBE_D3BJ_same_setting_meV_cell', f'E(iirfix) − E(i), PBE-D3(BJ)/SG15 {vlab}, per cell',
        case='periodic_qe', check=d['PBE_D3_same_setting_meV_cell'], tol=0.01)
    add(b + '.corr_meV', c['hybrid_correction_meV_cell'], 'meV', 1, src + '.hybrid_correction_meV_cell', f'hybrid correction (PBE0-D3 − PBE-D3 difference), {vlab}, per cell',
        case='periodic_qe', check=d['correction_meV_cell'], tol=0.01)
    add(b + '.corr_kJ', d['correction_kJ_fu'], 'kJ mol-1', 1, 'derived/d11_burst2.json:periodic.hybrid (raw outputs)',
        f'hybrid correction per formula unit, {vlab}', case='periodic_qe', check=c['hybrid_correction_kJ_mol_per_formula_unit'], tol=1e-3)
    add(b + '.best_kJ', d['best_estimate_kJ_fu'], 'kJ mol-1', 1, 'derived/d11_burst2.json:periodic.hybrid (SSSP PBE-D3 difference + correction)',
        f'constrained neutral − salt at PBE0-D3(BJ), per formula unit ({vlab} correction)', case='periodic_qe',
        check=c['best_estimate_PBE0_D3BJ_kJ_mol_per_formula_unit'], tol=1e-3)
    if vk == 'e100':
        add(b + '.corr_exx_meV', c['split_of_the_correction_meV_cell']['exact_exchange_25pct'], 'meV', 1, src + '.split_of_the_correction_meV_cell.exact_exchange_25pct',
            'exact-exchange part of the correction, per cell', case='periodic_qe', check=d['split_exact_exchange_meV_cell'], tol=0.01)
        add(b + '.corr_d3_meV', c['split_of_the_correction_meV_cell']['D3_parameters_pbe0_vs_pbe'], 'meV', 1, src + '.split_of_the_correction_meV_cell.D3_parameters_pbe0_vs_pbe',
            'part of the correction from the pbe0 instead of pbe D3(BJ) parameters, per cell', case='periodic_qe', check=d['split_D3_params_meV_cell'], tol=0.01)
        add(b + '.corr_exx_kJ', d['split_exact_exchange_meV_cell'] / 1000 * 96.485332 / 2, 'kJ mol-1', 1, 'periodic.hybrid.e100.corr_exx_meV per formula unit',
            'exact-exchange part of the correction per formula unit', case='periodic_qe')
        add(b + '.corr_d3_kJ', d['split_D3_params_meV_cell'] / 1000 * 96.485332 / 2, 'kJ mol-1', 1, 'periodic.hybrid.e100.corr_d3_meV per formula unit',
            'D3-parameter part of the correction per formula unit', case='periodic_qe')
        add(b + '.best_meV', d['best_estimate_meV_cell'], 'meV', 1, 'derived/d11_burst2.json:periodic.hybrid.sg15_e100.iirfix_minus_i.best_estimate_meV_cell',
            'constrained neutral − salt at PBE0-D3(BJ), per cell', case='periodic_qe', check=c['best_estimate_PBE0_D3BJ_meV_cell'], tol=0.01)
        add(b + '.best_kcal', d['best_estimate_kJ_fu'] / 4.184, 'kcal mol-1', 1, 'periodic.hybrid.e100.best_kJ / 4.184',
            'the same in kcal mol⁻¹ per formula unit', case='periodic_qe')
        add(b + '.PBE0_kJ', d['PBE0_D3_kJ_fu'], 'kJ mol-1', 1, 'derived/d11_burst2.json:periodic.hybrid.sg15_e100.iirfix_minus_i.PBE0_D3_kJ_fu',
            'E(iirfix) − E(i) straight from the PBE0-D3(BJ)/SG15 single points, per formula unit', case='periodic_qe',
            check=c['PBE0_D3BJ_kJ_mol_per_formula_unit'], tol=1e-3)
        add(b + '.SG15_minus_SSSP_kJ', d['PBE_SG15_minus_SSSP_meV_cell'] / 1000 * 96.485332 / 2, 'kJ mol-1', 1,
            'derived/d11_burst2.json:periodic.hybrid.sg15_e100.iirfix_minus_i.PBE_SG15_minus_SSSP_meV_cell per formula unit',
            'PBE-D3(BJ) difference with SG15 (100/400 Ry) minus the SSSP relax difference: pseudopotential sensitivity of the PBE value', case='periodic_qe',
            check=c['PBE_D3BJ_SG15_minus_SSSP_meV_cell'] / 1000 * 96.485332 / 2, tol=0.001)
add('periodic.hybrid.corr_spread_kJ', dp['iirfix_correction_spread_kJ_fu'], 'kJ mol-1', 1, 'derived/d11_burst2.json:periodic.iirfix_correction_spread_kJ_fu',
    'spread of the hybrid correction over the three variants (cutoff 80/100 Ry, q mesh 1×1×1/2×1×1)', case='periodic_qe')
cii = HY['comparisons']['sg15_e100']['E_ii_minus_i']
add('periodic.hybrid.ii_best_kJ', dp['hybrid']['sg15_e100']['ii_minus_i']['best_estimate_kJ_fu'], 'kJ mol-1', 2, 'derived/d11_burst2.json:periodic.hybrid.sg15_e100.ii_minus_i',
    'E(ii) − E(i) at PBE0-D3(BJ) per formula unit: the two relaxed salt cells stay degenerate with the hybrid', case='periodic_qe',
    check=cii['best_estimate_PBE0_D3BJ_kJ_mol_per_formula_unit'], tol=1e-3)
add('periodic.hybrid.verdict', HY['verdict'], None, None, f'{PS}:hybrid.verdict', 'outcome of the hybrid check (salt lower at PBE0-D3(BJ))', case='periodic_qe')
CAVEATS.append('Periodic constrained neutral (relax_iirfix): a constrained structure, not a minimum. With O5–H1A held at 0.98 Å the rest of '
               'the cell relaxes (free-atom forces below 5e-4 Ry/Bohr) while the held H1A carries 0.103 Ry/Bohr (2.6 eV/Å) pointing to N1, '
               'so the proton would return to N1. +57.1 (PBE-D3(BJ)) and +52.8 kJ/mol per formula unit (with the PBE0-D3(BJ) correction, '
               '−4.3) are the cost of holding the proton on O5 at a normal O–H length in this lattice; the value depends on the chosen '
               'O5–H1A (no periodic scan of it was run). The hybrid correction is evaluated at PBE-D3(BJ) geometries, with SG15 potentials '
               '(not the SSSP set of the relaxations; the PBE difference itself changes by −0.8 kJ/mol per formula unit between the two '
               'sets) and a 2×1×1 k grid; cutoff (80/100 Ry) and q mesh (1×1×1/2×1×1) move the correction by 0.46 kJ/mol. The correction '
               'is 97 % exact exchange (−87.4 of −90.0 meV/cell) and 3 % the change of D3(BJ) parameters.')
MISSING.append('No periodic proton-transfer profile: the neutral form in the crystal is represented by one constrained point (O5–H1A '
               '0.98 Å, relax_iirfix); no hybrid relaxation of the cell.')

# =====================================================================================================================
section('29. DFT strain of the docked poses (burst 2)',
        'Single points B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM water, grid 4, of the representative pose of the largest cluster '
        'of each production job of the three sulfanilamide species, relative to the IEF-PCM optimum of the same species '
        '(<case>/dock_strain.json). Each pose is the DFT molecule rebuilt on the docked heavy atoms (bond lengths and angles '
        'of the optimum, torsions of the pose; dock/production_analysis.py); "Hrot" = N–H rotors of N1 and N2 re-optimised '
        'with MMFF94s at the docked heavy atoms, because the scoring functions do not score hydrogen positions; "as docked" = '
        'hydrogens where the docking left them (a sensitivity case only). Rigid single points, no constrained relaxation. '
        'Re-derived (derive/d11_burst2.py): strain recomputed from each E_Ha, every geometry re-hashed (sha256 as recorded), '
        'τ_fold re-measured, the heavy atoms of the two variants identical, E(opt) equal to the case optimisation.')
SPL = {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate', 'sfa_cation': 'cation'}
add('dock.strain.method', 'B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM (water), rigid single points of the cluster-representative docked conformer '
    '(DFT bond lengths and angles, docked torsions; N–H rotors re-optimised with MMFF94s), minus the IEF-PCM optimum of the same species', None, None,
    'molecules/sfa/results/<case>/dock_strain.json:provenance; dock/results/poses_for_dft/manifest.json:how', 'DFT strain protocol')
for sp, s in D11['strain'].items():
    src = s['source']
    add(f'dock.strain.{sp}.n', s['n_poses'], None, 0, f'{src}:geometries', f'{SPL[sp]}: docked poses with a DFT single point (one per receptor and scoring function)')
    add(f'dock.strain.{sp}.hrot_min', s['hrot_range'][0], 'kcal mol-1', 2, f'{src}:energies (own min over the __Hrot poses)', f'{SPL[sp]}: smallest DFT strain, N–H rotors re-optimised')
    add(f'dock.strain.{sp}.hrot_max', s['hrot_range'][1], 'kcal mol-1', 2, f'{src}:energies (own max over the __Hrot poses)', f'{SPL[sp]}: largest DFT strain, N–H rotors re-optimised')
    add(f'dock.strain.{sp}.asdocked_min', s['asdocked_range'][0], 'kcal mol-1', 2, f'{src}:energies (own min, as-docked H)', f'{SPL[sp]}: smallest DFT strain, hydrogens as docked (sensitivity only)')
    add(f'dock.strain.{sp}.asdocked_max', s['asdocked_range'][1], 'kcal mol-1', 2, f'{src}:energies (own max, as-docked H)', f'{SPL[sp]}: largest DFT strain, hydrogens as docked (sensitivity only)')
    add(f'dock.strain.{sp}.checks_ok', bool(s['all_sha256_ok'] and abs(s['E_opt_diff_Ha']) < 1e-7 and s['max_tau_dev_vs_manifest'] < 0.05
                                             and all(p_['heavy_max_dev_hrot_vs_asdocked_A'] < 1e-6 for p_ in s['poses'].values())), None, None,
        'derived/d11_burst2.json:strain', f'{SPL[sp]}: geometry hashes, E(opt) identity, τ_fold and heavy atoms of both variants all verified')
    js = J(RES / s['case'] / 'dock_strain.json')
    e0 = js['energies']['b3lyp-d3bj']['opt']['E_Ha']
    for base, p_ in s['poses'].items():
        pdb = p_['target'].split('_')[1]
        k = f"dock.strain.{pdb}.{sp}.{p_['scoring']}"
        for kind, suf, lab in (('hrot', '__Hrot', 'N–H rotors re-optimised'), ('asdocked', '', 'hydrogens as docked')):
            ej = js['energies']['b3lyp-d3bj'][base + suf]['E_Ha']
            add(f'{k}.{kind}', (ej - e0) * HA_KCAL, 'kcal mol-1', 2, f'{src}:energies["b3lyp-d3bj"]["{base}{suf}"].E_Ha − ["opt"].E_Ha',
                f"{pdb} {p_['scoring']} {SPL[sp]}, seed {p_['seed']} (largest cluster, {p_['cluster_population']}): DFT strain, {lab}",
                check=p_[f'strain_{kind}_kcal'], tol=1e-6)
        add(f'{k}.seed', p_['seed'], None, 0, 'dock/results/poses_for_dft/manifest.json', f"{pdb} {p_['scoring']} {SPL[sp]}: seed of the representative pose")
        add(f'{k}.tau', p_['tau_fold_hrot'], 'deg', 1, f'derived/d11_burst2.json:strain.{sp}.poses["{base}"].tau_fold_hrot',
            f"{pdb} {p_['scoring']} {SPL[sp]}: τ_fold of the representative pose", check=p_['tau_fold_manifest'], tol=0.05)
_rg = lambda sp, kind: '–'.join(str(_Dec(repr(x)).quantize(_Dec('0.1'), rounding=_RHE)) for x in D11['strain'][sp][f'{kind}_range'])
CAVEATS.append('DFT strain of the docked poses (section 29) is a rigid single point per job (the cluster representative): the DFT bond '
               'lengths and angles are kept and only the torsions of the pose applied, and the N–H hydrogens are placed by MMFF94s, not '
               'by DFT. It is therefore an upper bound to the strain of that conformer, not a relaxed torsional strain (no constrained '
               'optimisation was run), and it is larger than the MMFF94s estimates (≤ 1.8 kcal/mol) for that reason. Sulfanilamide: '
               f"{_rg('sfa_neutral', 'hrot')} kcal/mol (the largest value is {D11['strain']['sfa_neutral']['hrot_range'][1]:.4f}, so 5.1, not 5.2 from "
               f"re-rounding the 2-decimal 5.15 of dock_strain.md); sulfonamidate {_rg('sfa_anion', 'hrot')}; cation {_rg('sfa_cation', 'hrot')}. "
               f"With the hydrogens as docked the values reach {_rg('sfa_neutral', 'asdocked').split('–')[1]} kcal/mol, which measures the "
               'arbitrary H orientations of the docking output, not the pose: sensitivity only, never quoted as strain.')

# =====================================================================================================================
section('30. Withdrawn facts and the rounding rule',
        'Two final passes over all facts. (a) Withdrawn: the relax_iifix energies, structure metrics and proton metrics (H1A '
        'alone held in space; the lattice slid by 0.66 Å and N1 re-bonded to H1A, so the "constrained neutral" of burst 1 is a '
        'salt). Kept for provenance with value null, the number in withdrawn_value and status WITHDRAWN; numbers_audit.py '
        'therefore cannot match a printed number to them. Never cite. Only the three slide diagnostics (periodic.iifix_shift_mag, '
        '_shift_c, _N1_moved) stay citable, for a methodological note on why O5 must be held as well. (b) Rounding: every printed value '
        'is the full-precision value rounded half to even at its printed precision. Where a source stores a value already rounded '
        'to one digit beyond the printed precision and that digit is a 5 (a tie), the full-precision value re-derived from the same '
        'raw data replaces it (for the X-ray bonds: the full-precision .res coordinates, not the _atom_site coordinates rounded to '
        'their s.u.); the stored value stays in the note.')
# the slide diagnostics stay citable, and only as the methodological note (why O5 must be held as well)
IIFIX_NOTE = ['periodic.iifix_shift_mag', 'periodic.iifix_shift_c', 'periodic.iifix_N1_moved']
for k in IIFIX_NOTE:
    FACTS[k]['status'] = ('methodological note only: the failed constraint of relax_iifix (H1A alone held; the lattice slid and N1 '
                          're-bonded to H1A); its energies and structure metrics are withdrawn')
WITHDRAWN = [k for k in FACTS if ('iifix' in k and 'iirfix' not in k and k not in IIFIX_NOTE)]
for k in WITHDRAWN:
    f_ = FACTS[k]
    f_['withdrawn_value'] = f_['value']
    f_['value'] = None
    f_['withdrawn'] = True
    f_['status'] = ('WITHDRAWN (never cite): relax_iifix held H1A alone, the lattice slid and N1 re-bonded to H1A, so this is a salt '
                    'structure; superseded by relax_iirfix (section 28)')
add('withdrawn.keys', sorted(WITHDRAWN), None, None, 'facts_burst2.py section 30 (keys containing "iifix" except the slide diagnostics)',
    'facts withdrawn in burst 2 (kept for provenance, never cited)')


def rhe(x, p):
    """x rounded half to even at p decimals, as a string (the digits of the float's shortest repr are the value)."""
    return str(_Dec(repr(float(x))).quantize(_Dec(1).scaleb(-p), rounding=_RHE))


TIES = []
for _k, _f in FACTS.items():
    _v, _p, _rd = _f.get('value'), _f.get('precision'), _f.get('rederived')
    if not (isinstance(_v, float) and isinstance(_p, int) and _p >= 0 and isinstance(_rd, (int, float)) and not isinstance(_rd, bool)):
        continue
    _t = _Dec(repr(_v)).as_tuple()
    if -_t.exponent != _p + 1 or _t.digits[-1] != 5:
        continue                                                   # the stored value is not a tie at the printed precision
    if abs(_rd - _v) > 0.5 * 10 ** (-(_p + 1)) + 1e-12:
        continue                                                   # the re-derivation is not the same number at full precision
    _was = f'{_v:.{_p}f}'.replace('-', '−')
    _f['stored_value'] = _v
    _f['value'] = float(_rd)
    _f['note'] = ((_f.get('note') + '; ') if _f.get('note') else '') + (
        f'rounding tie decided 29 Sep 2026: the source stores {_v!r} (a tie at {_p} decimals, which float formatting printed as {_was}); '
        f'value = full-precision re-derivation {_rd:.{_p + 4}f}, printed {rhe(_rd, _p).replace("-", "−")} (round half to even on the full precision)')
    TIES.append((_k, _v, _was, _rd, rhe(_rd, _p)))
DISAGREE[:] = [d_ for d_ in DISAGREE if not d_.startswith('Rounding ties still open')]
if TIES:
    DISAGREE.append('SETTLED 29 Sep 2026, rounding ties (rule: round half to even on the full-precision value; a stored value that is itself '
                    'a tie is replaced by the full-precision re-derivation from the same raw data, the stored value kept in the note): '
                    + '; '.join(f'{k_} stored {v_!r} (printed {w_}) → full precision {rd_:.6f}, prints {n_.replace("-", "−")}' for k_, v_, w_, rd_, n_ in TIES)
                    + '. xray.P1-O6 is re-derived from the full-precision .res coordinates (1.497510 Å): the earlier 1.497493 came from the '
                    '_atom_site coordinates rounded to their s.u., so 1.498 (as in the CIF, 1.4975(14)) stands.')
ROUNDING_RULE = ('Printed value = full-precision value rounded half to even at the fact precision. Stored values that are themselves '
                 'rounded ties are replaced by the full-precision re-derivation from the same raw data (stored_value and note keep the '
                 'original); tables (build_tables.py f) and figure texts (_saph.num) format with the same rule.')
