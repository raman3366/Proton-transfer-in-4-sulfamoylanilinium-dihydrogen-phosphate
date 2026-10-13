"""Facts added for the version of 5 Oct 2026 (Journal of the Indian Chemical Society; section 39 of facts.md), executed
by build_facts.py after facts_control.py (exec in its namespace: add(), section(), FACTS, J, DER are defined there).

  39. Points of the pre-submission review that needed a number the earlier sections did not hold: the construction of
      the one-dimensional proton model and its sensitivity to the continuation of the scanned profile and to the
      effective mass; the zero-point energy of the separated species; the shift of the proton with def2-TZVP; the
      wavelength at which aqueous sulfanilamide is determined by its absorption, against the computed band.  Source:
      derived/d16_jics.json (derive/d16_jics.py, code of our own), which also reproduces the stored 1-D model.
"""
D16 = J(DER / 'd16_jics.json')
D16S = 'derived/d16_jics.json'
O1_ = D16['oned']

section('39. Review before submission (5 Oct 2026): 1-D proton model and its sensitivity, zero-point energy of the '
        'separated species, basis-set shift of the proton, measured absorption of sulfanilamide',
        'The one-dimensional proton model of section 8 (dft/zpe_barrier.py) is rebuilt in derive/d16_jics.py: lower '
        'envelope of the two relaxed scan chains, PCHIP slopes set to zero at the two unconstrained minima, cubic '
        'Hermite interpolation in d(N1–H1A), harmonic walls beyond the scanned range, sinc-DVR on 801 points over the '
        'range plus 0.6 Å on each side, N–H reduced mass. With the stored wall constant on both sides it reproduces '
        f"the stored level, its distance from the barrier top and the O5-side probability to {O1_['max_abs_diff_vs_stored']:.0e} "
        '(kcal mol⁻¹ and probability). The wall constant is the N1–H1A local-mode force constant of the salt; on the '
        'O5 side, where the relaxed coordinate beyond the neutral minimum is the soft N···O stretch, it has no '
        f"physical basis, so that wall, the N1-side wall and the mass were varied ({O1_['n_variants']} variants). The "
        'sign of the distance from the barrier top is not determined by the model; that the ground state lies mainly '
        'on the O5 side is. The other entries are differences of values of earlier sections, recomputed from the raw '
        'files, and one literature value.')

# ------------------------------------------------------------------------------------------------- 1-D proton model
add('oned.re.E0', O1_['reference']['E0'], 'kcal mol-1', 2, f'{D16S}:oned.reference.E0',
    '1-D proton model rebuilt in derive/d16_jics.py: ground level above the salt minimum',
    check=FACTS['zpe.oned_E0']['value'], tol=0.005)
add('oned.re.E0_minus_top', O1_['reference']['E0_minus_top'], 'kcal mol-1', 2, f'{D16S}:oned.reference.E0_minus_top',
    '1-D proton model rebuilt: ground level minus the barrier top', check=FACTS['zpe.oned_E0_minus_top']['value'],
    tol=0.005)
add('oned.re.P_O', O1_['reference']['P_O'], None, 2, f'{D16S}:oned.reference.P_O',
    '1-D proton model rebuilt: ground-state probability on the O5 side', check=FACTS['zpe.oned_P_O']['value'], tol=0.005)
add('oned.wall_k', O1_['wall_k_mdyn_per_A'], 'mdyn A-1', 2,
    'molecules/sfa/results/sfa_pair_pcm_d3/scan/zpe_barrier.json:one_d_model.wall_k_mdyn_per_A',
    'harmonic walls of the 1-D model beyond the scanned range: force constant (the N1–H1A local-mode constant of the salt)')
add('oned.d_first', O1_['d_first_A'], 'A', 3, f'{D16S}:oned.d_first_A', 'first point of the scanned profile, d(N1–H1A) (salt minimum)')
add('oned.d_last', O1_['d_last_A'], 'A', 3, f'{D16S}:oned.d_last_A', 'last point of the scanned profile, d(N1–H1A) (neutral minimum)')
add('oned.grid', O1_['grid_points'], None, 0, 'dft/zpe_barrier.py; derive/d16_jics.py', 'sinc-DVR grid points of the 1-D model')
add('oned.pad', O1_['pad_A'], 'A', 1, 'dft/zpe_barrier.py; derive/d16_jics.py', 'extension of the DVR grid beyond the scanned range on each side')
for k_, v_ in O1_['right_wall'].items():
    add(f'oned.sens.right.k{k_}.k', float(k_), 'mdyn A-1', 1, f'{D16S}:oned.right_wall', '1-D model, sensitivity: force constant of the O5-side wall')
    add(f'oned.sens.right.k{k_}.E0_minus_top', v_['E0_minus_top'], 'kcal mol-1', 2, f'{D16S}:oned.right_wall[{k_}].E0_minus_top',
        f'1-D model with an O5-side wall of {k_} mdyn/Å: ground level minus the barrier top')
    add(f'oned.sens.right.k{k_}.P_O', v_['P_O'], None, 2, f'{D16S}:oned.right_wall[{k_}].P_O',
        f'1-D model with an O5-side wall of {k_} mdyn/Å: probability on the O5 side')
for k_, v_ in O1_['left_wall'].items():
    add(f'oned.sens.left.k{k_}.k', float(k_), 'mdyn A-1', 1, f'{D16S}:oned.left_wall', '1-D model, sensitivity: force constant of the N1-side wall')
    add(f'oned.sens.left.k{k_}.E0_minus_top', v_['E0_minus_top'], 'kcal mol-1', 2, f'{D16S}:oned.left_wall[{k_}].E0_minus_top',
        f'1-D model with an N1-side wall of {k_} mdyn/Å: ground level minus the barrier top')
    add(f'oned.sens.left.k{k_}.P_O', v_['P_O'], None, 2, f'{D16S}:oned.left_wall[{k_}].P_O',
        f'1-D model with an N1-side wall of {k_} mdyn/Å: probability on the O5 side')
for k_, v_ in O1_['mass'].items():
    add(f'oned.sens.mass.m{k_}.m', float(k_), 'u', 3 if k_ == '1.008' else 1, f'{D16S}:oned.mass', '1-D model, sensitivity: effective mass')
    add(f'oned.sens.mass.m{k_}.E0_minus_top', v_['E0_minus_top'], 'kcal mol-1', 2, f'{D16S}:oned.mass[{k_}].E0_minus_top',
        f'1-D model with an effective mass of {k_} u: ground level minus the barrier top')
    add(f'oned.sens.mass.m{k_}.P_O', v_['P_O'], None, 2, f'{D16S}:oned.mass[{k_}].P_O',
        f'1-D model with an effective mass of {k_} u: probability on the O5 side')
add('oned.sens.n_variants', O1_['n_variants'], None, 0, f'{D16S}:oned.n_variants', '1-D model: variants computed beside the reference (three O5-side walls, two N1-side walls, three masses)')
add('oned.sens.E0_minus_top_min', O1_['E0_minus_top_min'], 'kcal mol-1', 2, f'{D16S}:oned.E0_minus_top_min',
    '1-D model, all variants and the reference: lowest ground level relative to the barrier top')
add('oned.sens.E0_minus_top_max', O1_['E0_minus_top_max'], 'kcal mol-1', 2, f'{D16S}:oned.E0_minus_top_max',
    '1-D model, all variants and the reference: highest ground level relative to the barrier top')
add('oned.sens.P_O_min', O1_['P_O_min'], None, 2, f'{D16S}:oned.P_O_min', '1-D model, all variants and the reference: smallest probability on the O5 side')
add('oned.sens.P_O_max', O1_['P_O_max'], None, 2, f'{D16S}:oned.P_O_max', '1-D model, all variants and the reference: largest probability on the O5 side')
add('oned.within', 0.1, 'kcal mol-1', 1, f'{D16S}:oned.reference.E0_minus_top (|−0.06| < 0.1)',
    'bound quoted in the text for the distance of the reference ground level from the barrier top')

# -------------------------------------------------------------------------- zero-point energy of the separated species
add('pt.sep.dZPE', D16['zpe_separated']['dZPE_kcal'], 'kcal mol-1', 2, f'{D16S}:zpe_separated.dZPE_kcal (freq.json:ZPE_Ha of the four gas-phase cases)',
    'ZPE(sulfanilamide) + ZPE(H3PO4) − ZPE(cation) − ZPE(H2PO4−), gas-phase harmonic frequencies (negative = neutral side lower)')
add('pt.sep.dZPE_abs', abs(D16['zpe_separated']['dZPE_kcal']), 'kcal mol-1', 1, f'{D16S}:zpe_separated.dZPE_kcal',
    'the same as a magnitude: zero-point energy widens the gap between the continuum exchange energy and the measured ΔpKa by this much')

# ---------------------------------------------------------------------------------------------- def2-TZVP geometry
add('proton.sens.tz.delta_shift', D16['tz']['delta_shift'], 'A', 3, f'{D16S}:tz.delta_shift (opt.xyz of sfa_pair_pcm_d3 and sfa_pair_pcm_d3_tz)',
    'IEF-PCM ion pair: δ(def2-TZVP) − δ(6-311++G(d,p)); positive = H1A moves towards O5',
    check=FACTS['proton.sens.tz.delta']['value'] - FACTS['proton.pcm.delta']['value'], tol=0.0015)

# ------------------------------------------------------------------------------------------------- TD-DFT anchor
add('lit.uv.sfa_water.nm', 258, 'nm', 0, 'wei2025boron (Molecules 2025, 30, 1735): "The concentration of residual sulfanilamide was quantified using '
    'UV–visible spectroscopy at a wavelength of λ = 258 nm" (aqueous sulfanilamide, no reagent); the same wavelength in bai2023eco '
    '(Polymers 2023, 15, 2342) and, as absorption maximum of sulfanilamide in a reagent assay, in faris2025spectrophotometric',
    'wavelength at which aqueous sulfanilamide is determined by UV–Vis absorption (literature)')
add('td.sfa.cam.bright.minus_exp', FACTS['td.sfa.cam.bright.E']['value'] - D16['uv']['lit_eV'], 'eV', 2,
    f'td.sfa.cam.bright.E − hc/(258 nm) ({D16S}:uv.lit_eV)',
    'sulfanilamide, CAM-B3LYP with non-equilibrium IEF-PCM: computed S2 energy minus hc/λ of the wavelength at which aqueous '
    'sulfanilamide is determined (258 nm)')

# --------------------------------------------------------------------------------------------- docking: review items
# the sulfonamide nitrogen on the zinc in the eight AutoDock4Zn runs of sulfanilamide (two production runs and the six
# input-geometry runs listed under "diagnostics"); receptor minimisation and search settings stated in Section S13
PS_ = J(DOCK / 'results' / 'production_summary.json')
PSS_ = 'dock/results/production_summary.json'
RUNS8_ = ['T3_6RL9/sfa_neutral__ad4zn__prod', 'T4_3HS4/sfa_neutral__ad4zn__prod'] + list(PS_['diagnostics'])
assert len(RUNS8_) == 8, RUNS8_
ZN_ = {k_: PS_['jobs'][k_]['site_checks']['Zn-N(sulfonamide)'] for k_ in RUNS8_}
OG_ = {k_: PS_['jobs'][k_]['site_checks']['Thr199 OG1 H-bond'] for k_ in RUNS8_}
add('dock.ad4zn8.n_runs', len(RUNS8_), None, 0, f'{PSS_}:jobs (two production runs) and diagnostics (six runs)',
    'AutoDock4Zn runs of neutral sulfanilamide in 6RL9 and 3HS4 (production and input-geometry runs)')
add('dock.ad4zn8.ZnN_min', min(v_['distance_mean'] for v_ in ZN_.values()), 'A', 2, f'{PSS_}:jobs[*].site_checks["Zn-N(sulfonamide)"].distance_mean',
    'smallest run mean of Zn–N(sulfonamide) over the eight AutoDock4Zn runs of sulfanilamide')
add('dock.ad4zn8.ZnN_max', max(v_['distance_mean'] for v_ in ZN_.values()), 'A', 2, f'{PSS_}:jobs[*].site_checks["Zn-N(sulfonamide)"].distance_mean',
    'largest run mean of Zn–N(sulfonamide) over the eight AutoDock4Zn runs of sulfanilamide')
add('dock.ad4zn8.OG1_pass_min', min(v_['pass_fraction'] for v_ in OG_.values()), None, 2, f'{PSS_}:jobs[*].site_checks["Thr199 OG1 H-bond"].pass_fraction',
    'smallest fraction of seeds with the Thr199 OG1 contact over the eight runs')
add('dock.ad4zn8.3HS4_dft.ZnN_pass', ZN_['T4_3HS4/sfa_neutral__ad4zn__prod']['pass_fraction'], None, 2,
    f'{PSS_}:jobs["T4_3HS4/sfa_neutral__ad4zn__prod"].site_checks', '3HS4, DFT geometry: fraction of seeds with Zn–N within 2.3 Å')
add('dock.ad4zn8.3HS4_csn102.ZnN_pass', ZN_['T4_3HS4/sfa_neutral__ad4zn__csn102']['pass_fraction'], None, 2,
    f'{PSS_}:jobs["T4_3HS4/sfa_neutral__ad4zn__csn102"].site_checks',
    '3HS4, DFT geometry with the C4–S1–N2 angle of the ETKDG conformer: fraction of seeds with Zn–N within 2.3 Å')


def _heavy_sdf(path_):
    """Heavy atoms of the first molecule of an SDF file (V2000 counts line, atom block)."""
    ls_ = path_.read_text().splitlines()
    n_ = int(ls_[3][:3])
    return sum(1 for l_ in ls_[4:4 + n_] if l_[31:34].strip() != 'H')


add('dock.lig.azm.heavy', _heavy_sdf(DOCK / 'ligands' / 'azm_neutral.sdf'), None, 0, 'dock/ligands/azm_neutral.sdf (atom block)', 'heavy atoms of acetazolamide')
add('dock.lig.sfa.heavy', _heavy_sdf(DOCK / 'ligands' / 'sfa_neutral.sdf'), None, 0, 'dock/ligands/sfa_neutral.sdf (atom block)', 'heavy atoms of sulfanilamide',
    check=FACTS['admet.sfa.heavy_atoms']['value'], tol=0)
PREP_ = {t_: J(DOCK / t_ / 'prep_report.json')['minimisation'] for t_ in ('T1_3TZF', 'T2_1AJ0', 'T3_6RL9', 'T4_3HS4', 'T5_7S2K')}
assert len({m_['restraint_k_kJ_mol_nm2'] for m_ in PREP_.values()}) == 1
add('dock.prep.restraint_k', PREP_['T1_3TZF']['restraint_k_kJ_mol_nm2'], 'kJ mol-1 nm-2', 0, 'dock/T*/prep_report.json:minimisation.restraint_k_kJ_mol_nm2',
    'restraint on the deposited heavy atoms in the minimisation of every receptor (10⁵ kJ mol⁻¹ nm⁻²; Amber14 force field, OpenMM)')
add('dock.prep.restraint_k_exp', 5, None, 0, 'dock/T*/prep_report.json:minimisation.restraint_k_kJ_mol_nm2 = 1e5', 'the same as a power of ten (printed 10⁵)')
add('dock.prep.max_disp', max(m_['max_displacement_restrained_A'] for m_ in PREP_.values()), 'A', 2, 'dock/T*/prep_report.json:minimisation.max_displacement_restrained_A',
    'largest displacement of a deposited atom in the restrained minimisation, over the five receptors of Section 2.5')
add('dock.ad4zn.spacing', 0.375, 'A', 3, 'dock/run_ad4zn.py:SPACING; dock/T3_6RL9/ad4zn/*.gpf', 'AutoGrid spacing of the AutoDock4Zn maps')
add('dock.ad4zn.npts', 54, None, 0, 'dock/run_ad4zn.py:write_gpf; dock/T3_6RL9/ad4zn/*.gpf', 'AutoGrid npts per axis (the number of intervals)')
add('dock.ad4zn.points', 55, None, 0, 'npts + 1 (AutoGrid counts intervals)', 'grid points per axis of the AutoDock4Zn maps')
add('dock.ad4zn.box', 54 * 0.375, 'A', 2, 'dock/run_ad4zn.py:write_gpf (npts × spacing)', 'edge of the AutoDock4Zn grid')
add('dock.vina.num_modes', 20, None, 0, 'dock/T3_6RL9/vina_config.txt:num_modes; dock/run_ad4zn.py', 'output modes per Vina run')
add('dock.vina.energy_range', 3, 'kcal mol-1', 0, 'dock/T3_6RL9/vina_config.txt:energy_range; dock/run_ad4zn.py', 'energy range of the output modes')
add('software.ai.model', 5.5, None, 1, 'git log of the repository (Co-Authored-By trailers of every commit that carries one)', 'version of the Claude Opus model that wrote and ran the scripts (in Claude Code)')
add('protocol.docking.plip', 'PLIP 3.0.1', None, None, 'importlib.metadata.version("plip") in the docking environment, read 5 Oct 2026',
    'version of PLIP used for the contact analysis (the entry protocol.docking names the program without it)')
add('software.figures', 'Matplotlib 3.11.2, PyVista 0.49.0, VTK 9.7.0', None, None,
    'versions of the modules in the paper environment (matplotlib.__version__, pyvista.__version__, vtk.vtkVersion), read 5 Oct 2026',
    'libraries that draw the figures and the graphical abstract from the stored results (figure scripts and figkit)')

# ------------------------------------------------------------------------- the 2017 article on SAPH (cited again, 5 Oct 2026)
add('lit.saph2017.uv_cutoff.nm', 223, 'nm', 0, 'muthuselvi2017structural (Adv. Sci. Eng. Med. 2017, 9, 931), abstract: "lower cut off wavelength at 223 nm" '
    '(archived publisher page, read 5 Oct 2026; the abstract speaks of "the transparency of the crystal" and does not state the form of '
    'the sample: crystal, powder or solution)',
    'reported lower UV–Vis cut-off wavelength of SAPH (literature)')
add('lit.saph2017.optical_gap.eV', 5.57, 'eV', 2, 'muthuselvi2017structural, abstract: "The optical band of the title compound was found to as 5.57 eV"',
    'reported optical band gap of SAPH (literature; not printed in the paper)')
