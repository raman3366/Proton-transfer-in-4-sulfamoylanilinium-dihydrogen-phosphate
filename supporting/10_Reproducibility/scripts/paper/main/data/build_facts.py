#!/usr/bin/env python3
"""Build paper/main/data/facts.json and facts.md -- every number the SAPH paper and SI may print.

Each fact: key -> {value, unit, precision (decimals to print), case, source (file + JSON key or log line),
meaning, note?, status?, rederived?, agrees?}.  `rederived` is the value obtained by the independent scripts in
derive/ (outputs in derived/); `agrees` says whether it matches the source within the tolerance stated in the
script.  Run the derive scripts first (see derive/README in facts.md), then:

    /root/workspace/cqc-dft/tools/paper-env/bin/python paper/main/data/build_facts.py

Reads only; writes only facts.json and facts.md next to this file.
"""
import json
import math
import re
import statistics as _stt
from decimal import Decimal as _DecB
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = REPO / 'molecules/sfa/results'
DOCK = REPO / 'dock'
DER = HERE / 'derived'
ARCH = '/root/cqc-sfa-archive/vm/cqc-sfa'
HA_KCAL, HA_KJ, HA_EV = 627.5095, 2625.4996, 27.211386


def J(p):
    return json.loads(Path(p).read_text())


def r(p):
    return str(Path(p).relative_to(REPO)) if str(p).startswith(str(REPO)) else str(p)


D1, D2, D3, D4 = J(DER / 'd01_xray.json'), J(DER / 'd02_models.json'), J(DER / 'd03_geometry.json'), J(DER / 'd04_scan.json')
D5, D6, D7, D8, D9 = (J(DER / f) for f in ['d05_freq.json', 'd06_electronic.json', 'd07_hirshfeld.json',
                                           'd08_periodic.json', 'd09_dock.json'])

FACTS = {}
SECTIONS = []          # (section title, intro, [keys])
DISAGREE, CAVEATS, MISSING = [], [], []
cur = None


def section(title, intro=''):
    global cur
    cur = (title, intro, [])
    SECTIONS.append(cur)


def su_num(s):
    m = re.match(r'^(-?[\d.]+)', str(s))
    return float(m.group(1)) if m else None


def add(key, value, unit, prec, source, meaning, case=None, note=None, status=None, check=None, tol=None):
    assert key not in FACTS, key
    rec = dict(value=value, unit=unit, precision=prec, case=case, source=source, meaning=meaning)
    if note:
        rec['note'] = note
    if status:
        rec['status'] = status
    if check is not None and isinstance(value, (int, float)) and not isinstance(value, bool):
        t = tol if tol is not None else (0.5 * 10 ** (-prec) if isinstance(prec, int) else 1e-6)
        rec['rederived'] = check
        rec['agrees'] = bool(abs(check - value) <= t)
        if not rec['agrees']:
            DISAGREE.append(f'{key}: source {value} vs re-derived {check} (tolerance {t})')
    FACTS[key] = rec
    cur[2].append(key)


# =====================================================================================================================
section('1. Crystal structure (deposited CIF)',
        'From e-70-0o997-sup1.cif (Acta Cryst. E 2014, 70, o997). Hydrogen bonds recomputed from _atom_site with our '
        'own code (derive/d01_xray.py): every D···A reproduces the deposited value to within '
        f"{D1['max_abs_dDA']:.4f} Å; H···A and angles are checked against the full-precision coordinates of the "
        'embedded .res (the _atom_site H coordinates carry only 3 decimals). N–H hydrogens were refined freely, C–H '
        'riding (AFIX 43, 0.93 Å), O–H restrained with DFIX 0.82 Å: the O–H positions are not experimental.')
C = D1['cif_items']
CIFS = 'e-70-0o997-sup1.cif'
add('xray.ccdc', 1016950, None, 0, f'{CIFS}:_database_code_depnum_ccdc_archive', 'CCDC deposition number')
add('xray.doi', C['_journal_paper_doi'], None, None, f'{CIFS}:_journal_paper_doi', 'DOI of the crystal-structure paper')
add('xray.journal', 'Acta Cryst. E 70 (2014) o997-o998', None, None, f'{CIFS}:_journal_*', 'crystal-structure paper')
add('xray.formula_moiety', C['_chemical_formula_moiety'], None, None, f'{CIFS}:_chemical_formula_moiety', 'moiety formula (1 cation + 1 anion)')
add('xray.formula_sum', C['_chemical_formula_sum'], None, None, f'{CIFS}:_chemical_formula_sum', 'sum formula')
add('xray.mr', float(C['_chemical_formula_weight']), 'g mol-1', 2, f'{CIFS}:_chemical_formula_weight', 'formula weight')
add('xray.crystal_system', C['_space_group_crystal_system'], None, None, f'{CIFS}:_space_group_crystal_system', 'crystal system')
add('xray.space_group', 'Pc', None, None, f"{CIFS}:_space_group_name_H-M_alt '{C['_space_group_name_H-M_alt']}'", 'space group (No. 7, polar, non-centrosymmetric; Hall P -2yc)')
for ax in 'abc':
    s = C[f'_cell_length_{ax}']
    add(f'xray.cell_{ax}', su_num(s), 'A', 4, f'{CIFS}:_cell_length_{ax}', f'cell length {ax}', note=f'raw CIF value {s}')
add('xray.cell_beta', su_num(C['_cell_angle_beta']), 'deg', 3, f'{CIFS}:_cell_angle_beta', 'monoclinic angle', note=f"raw {C['_cell_angle_beta']}")
add('xray.cell_volume', su_num(C['_cell_volume']), 'A3', 2, f'{CIFS}:_cell_volume', 'cell volume', note=f"raw {C['_cell_volume']}",
    check=D1['cell_volume_recomputed_A3'], tol=0.01)
add('xray.z', 2, None, 0, f'{CIFS}:_cell_formula_units_Z', 'formula units per cell (Z′ = 1)')
add('xray.temperature', 294, 'K', 0, f'{CIFS}:_diffrn_ambient_temperature', 'data-collection temperature', note='294(2) K')
add('xray.radiation', 'Mo Kα, λ = 0.71073 Å; Bruker SMART APEX CCD', None, None, f'{CIFS}:_diffrn_radiation_*', 'radiation and diffractometer (the R4 template printed 0.73 Å: wrong)')
add('xray.wavelength', 0.71073, 'A', 5, f'{CIFS}:_diffrn_radiation_wavelength', 'Mo Kα wavelength')
add('xray.crystal', 'colourless needle, 0.28 × 0.18 × 0.10 mm', None, None, f'{CIFS}:_exptl_crystal_*', 'crystal description and size')
add('xray.density', float(C['_exptl_crystal_density_diffrn']), 'g cm-3', 3, f'{CIFS}:_exptl_crystal_density_diffrn', 'calculated density')
add('xray.mu', float(C['_exptl_absorpt_coefficient_mu']), 'mm-1', 3, f'{CIFS}:_exptl_absorpt_coefficient_mu', 'absorption coefficient (multi-scan correction)')
add('xray.f000', 280, None, 0, f'{CIFS}:_exptl_crystal_F_000', 'F(000)')
add('xray.n_refl_measured', int(C['_diffrn_reflns_number']), None, 0, f'{CIFS}:_diffrn_reflns_number', 'reflections measured')
add('xray.n_refl_unique', int(C['_reflns_number_total']), None, 0, f'{CIFS}:_reflns_number_total', 'independent reflections')
add('xray.n_refl_obs', int(C['_reflns_number_gt']), None, 0, f'{CIFS}:_reflns_number_gt', 'reflections with I > 2σ(I)')
add('xray.rint', float(C['_diffrn_reflns_av_R_equivalents']), None, 4, f'{CIFS}:_diffrn_reflns_av_R_equivalents', 'R_int')
add('xray.theta_max', float(C['_diffrn_reflns_theta_max']), 'deg', 3, f'{CIFS}:_diffrn_reflns_theta_max', 'θ_max (θ_min 1.876°)')
add('xray.r1_obs', float(C['_refine_ls_R_factor_gt']), None, 4, f'{CIFS}:_refine_ls_R_factor_gt', 'R1 [I > 2σ(I)]')
add('xray.r1_all', float(C['_refine_ls_R_factor_all']), None, 4, f'{CIFS}:_refine_ls_R_factor_all', 'R1 (all data)')
add('xray.wr2_all', float(C['_refine_ls_wR_factor_ref']), None, 4, f'{CIFS}:_refine_ls_wR_factor_ref', 'wR2 (all data)')
add('xray.gof', float(C['_refine_ls_goodness_of_fit_ref']), None, 3, f'{CIFS}:_refine_ls_goodness_of_fit_ref', 'goodness of fit S')
add('xray.n_params', 174, None, 0, f'{CIFS}:_refine_ls_number_parameters', 'refined parameters')
add('xray.n_restraints', 4, None, 0, f'{CIFS}:_refine_ls_number_restraints', 'restraints (DFIX 0.82 on O3–H3A and O4–H4A)')
add('xray.drho_max', 0.202, 'e A-3', 3, f'{CIFS}:_refine_diff_density_max', 'largest difference peak')
add('xray.drho_min', -0.220, 'e A-3', 3, f'{CIFS}:_refine_diff_density_min', 'deepest difference hole')
add('xray.flack', 0.069, None, 3, f'{CIFS}:_refine_ls_abs_structure_Flack', 'Flack parameter 0.069(16)')
add('xray.refinement', 'SHELXL2013; H treatment mixed: N–H free, C–H riding (AFIX 43), O–H DFIX 0.82 Å', None, None,
    f"{CIFS}:_computing_structure_refinement, _refine_ls_hydrogen_treatment, embedded .res ({'; '.join(D1['dfix'])})", 'refinement and hydrogen treatment')
add('xray.n_hbonds', D1['n_hbonds'], None, 0, f'{CIFS}:_geom_hbond loop', 'rows of the deposited H-bond table (7 N/O–H···O + 2 C–H···O)')
for h in D1['hbonds']:
    k = f"xray.hb.{h['D']}_{h['H']}_{h['A']}"
    c = h['cif']
    add(k + '.DA', su_num(c['DA']), 'A', 3, f"{CIFS}:_geom_hbond {h['D']} {h['H']} {h['A']} {h['sym']}",
        f"X-ray D···A of {h['D']}–{h['H']}···{h['A']} (acceptor at {h['sym']})",
        note=f"CIF: D–H {c['DH']}, H···A {c['HA']}, D···A {c['DA']}, angle {c['DHA']}°",
        check=h['calc_res']['DA'], tol=0.0006)
    add(k + '.HA', su_num(c['HA']), 'A', 2, f"{CIFS}:_geom_hbond", f"X-ray H···A of {h['D']}–{h['H']}···{h['A']}",
        check=h['calc_res']['HA'], tol=0.005)
    add(k + '.DHA', su_num(c['DHA']), 'deg', 0, f"{CIFS}:_geom_hbond", f"X-ray D–H···A angle of {h['D']}–{h['H']}···{h['A']}",
        check=h['calc_res']['DHA'], tol=0.5)
    add(k + '.DH', su_num(c['DH']), 'A', 2, f"{CIFS}:_geom_hbond", f"X-ray D–H of {h['D']}–{h['H']} (refined free / riding / DFIX)")
add('xray.symop_1_554', 'x, y, z − 1', None, None, 'SHELX symmetry code', 'anion partner of the ion pair (anion translated by −c)')
P = D1['pair']
for b, lab in [('C1-N1', 'C1–N1'), ('P1-O3', 'P1–O3 (P–OH)'), ('P1-O4', 'P1–O4 (P–OH)'), ('P1-O5', 'P1–O5 (P–O, accepts H1A)'),
               ('P1-O6', 'P1–O6 (P–O)')]:
    cif = D1['cif_bonds'].get(b) or D1['cif_bonds'].get('-'.join(b.split('-')[::-1]))
    add(f'xray.{b}', su_num(cif), 'A', 3, f'{CIFS}:_geom_bond {b} = {cif}', f'X-ray bond {lab}', check=D1['pair_res'][b], tol=0.0006)
for b, lab in [('C4-S1', 'C4–S1'), ('N2-S1', 'S1–N2'), ('O1-S1', 'S1–O1'), ('O2-S1', 'S1–O2')]:
    cif = D1['cif_bonds'][b]
    add(f'xray.{b}', su_num(cif), 'A', 3, f'{CIFS}:_geom_bond {b} = {cif}', f'X-ray bond {lab}')
add('xray.PO_long_mean', P['PO_long_mean'], 'A', 3, 'derived/d01_xray.json:pair.PO_long_mean', 'mean P–O(H) (O3, O4): protonated oxygens')
add('xray.PO_short_mean', P['PO_short_mean'], 'A', 3, 'derived/d01_xray.json:pair.PO_short_mean', 'mean P–O (O5, O6): deprotonated oxygens')
add('xray.PO_split', P['PO_split'], 'A', 3, 'derived/d01_xray.json:pair.PO_split', 'P–O(H) minus P–O mean split: H-position-independent evidence of H2PO4− (two long, two short)')
add('xray.d_N1_H1A', P['d_N1_H1A'], 'A', 3, 'derived/d01_xray.json:pair.d_N1_H1A (CIF 0.87(3))', 'X-ray N1–H1A (freely refined; X-ray N–H lengths are ~0.15 Å short)')
add('xray.d_H1A_O5', P['d_H1A_O5'], 'A', 3, 'derived/d01_xray.json:pair.d_H1A_O5 (CIF 1.89(3))', 'X-ray H1A···O5')
add('xray.d_N1_O5', P['d_N1_O5'], 'A', 3, 'derived/d01_xray.json:pair.d_N1_O5 (CIF 2.760(2))', 'X-ray N1···O5')
add('xray.delta', P['delta'], 'A', 3, 'derived/d01_xray.json:pair.delta', 'proton metric δ = d(N1–H1A) − d(H1A···O5) at the deposited H1A', note='X-ray H position; with a neutron-normalised N–H (1.009 Å) δ would be ≈ −0.75 Å')
add('xray.tau_C3C4S1N2', su_num(D1['cif_torsions']['C3-C4-S1-N2']), 'deg', 1, f"{CIFS}:_geom_torsion C3-C4-S1-N2 = {D1['cif_torsions']['C3-C4-S1-N2']}",
    'sulfonamide torsion C3–C4–S1–N2', check=P['tau_C3_C4_S1_N2'], tol=0.05)
add('xray.tau_C5C4S1N2', su_num(D1['cif_torsions']['C5-C4-S1-N2']), 'deg', 1, f"{CIFS}:_geom_torsion C5-C4-S1-N2 = {D1['cif_torsions']['C5-C4-S1-N2']}",
    'sulfonamide torsion C5–C4–S1–N2', check=P['tau_C5_C4_S1_N2'], tol=0.05)
add('xray.ang_C4S1N2', P['ang_C4_S1_N2'], 'deg', 1, 'derived/d01_xray.json:pair.ang_C4_S1_N2', 'C4–S1–N2 angle in the crystal (docking input-geometry diagnostic)')
add('xray.ang_O1S1O2', P['ang_O1_S1_O2'], 'deg', 1, 'derived/d01_xray.json:pair.ang_O1_S1_O2', 'O1–S1–O2 angle in the crystal')
o = D1['O4_H4A_O5']
add('xray.O4_O5', 2.500, 'A', 3, f'{CIFS}:_geom_hbond O4 H4A O5 1_455 = 2.500(2)', 'anion–anion O4–H4A···O5 donor–acceptor distance (short, strong)', check=o['OO'], tol=0.0005)
add('xray.O4H4A_frac', o['frac_along'], None, 2, 'derived/d01_xray.json:O4_H4A_O5.frac_along', 'fractional position of H4A along O4→O5 (restrained O–H 0.82 Å; not experimental)')
add('xray.ang_N1H1AO5', P['ang_N1_H1A_O5'], 'deg', 1, 'derived/d01_xray.json:pair.ang_N1_H1A_O5 (CIF 174(3))', 'X-ray N1–H1A···O5 angle from the _atom_site coordinates')
add('xray.ang_C6H6O4', P['ang_C6_H6_O4'], 'deg', 1, 'derived/d01_xray.json:pair.ang_C6_H6_O4 (CIF 156)', 'X-ray C6–H6···O4 angle from the _atom_site coordinates')
add('xray.C6_O4', P['C6_O4'], 'A', 3, 'derived/d01_xray.json:pair.C6_O4 (CIF 3.288(3))', 'X-ray C6···O4 (second leg of the R²₂(8) ring)')
add('xray.r22_8', 'N1–H1A···O5 + C6–H6···O4 close an R²₂(8) ring (H1A, N1, C1, C6, H6, O4, P1, O5)', None, None,
    'derived from the CIF H-bond table (both acceptors on the same anion, 1_554)', 'graph-set motif of the ion pair (ours; not described in the 2014 paper)')

# =====================================================================================================================
section('2. Models, charges and basis-set sizes',
        'Atoms, electrons and basis functions recomputed with PySCF 2.14 (spherical d) from the case inputs '
        '(derive/d02_models.py); every count equals the one the production runs logged. The ion pair has 445 basis '
        'functions, of which PySCF removed 2 near-linear dependencies (443 MOs; molden padded to 445 for Multiwfn).')
MOD = [('sfa_pair_pcm_d3', 'M1', 'ion pair C6H9N2O2S+ · H2PO4−, anion at x, y, z − 1'),
       ('sfa_pair_pcm_d3_tz', 'M1-TZ', 'ion pair, def2-TZVP'),
       ('cation_pcm_d3', 'M2', '4-sulfamoylanilinium cation'), ('h2po4_pcm_d3', 'M3', 'dihydrogen phosphate'),
       ('sulfanilamide_pcm_d3', 'M4', 'neutral sulfanilamide'), ('sulfanilamide_pcm_d3_tz', 'M4-TZ', 'neutral sulfanilamide, def2-TZVP'),
       ('sulfonamidate_pcm_d3', 'M5', 'sulfonamidate anion (N2-deprotonated)'), ('h3po4_pcm_d3', 'M6', 'phosphoric acid'),
       ('urea_gas_d3', 'M7', 'urea (NLO reference)')]
for c, m, what in MOD:
    d = D2[c]
    add(f'model.{m}.n_atoms', d['n_atoms'], None, 0, f"derived/d02_models.json:{c}.n_atoms ({d['xyz']})", f'{m} {what}: atoms', case=c)
    add(f'model.{m}.charge', d['charge'], 'e', 0, f'molecules/sfa/cases/{c}.json:charge', f'{m} {what}: charge (singlet)', case=c)
    add(f'model.{m}.nao', d['nao_spherical'], None, 0, f'derived/d02_models.json:{c}.nao_spherical (PySCF gto.M)', f'{m} {what}: basis functions ({d["basis"]})', case=c)
    add(f'model.{m}.n_electrons', d['n_electrons'], None, 0, f'derived/d02_models.json:{c}.n_electrons', f'{m} {what}: electrons', case=c)
add('model.M1.nmo_used', 443, None, 0, f'{ARCH}/molecules/sfa/results/sfa_pair_pcm_d3/logs/fmo.log ("molden padded from 443 to 445")',
    '6-311++G(d,p) ion pair: MOs after removal of 2 near-linear-dependent combinations')
add('model.M9.n_atoms', 54, None, 0, 'molecules/sfa/results/periodic_qe/build.json:nat', 'periodic Pc cell (Z = 2), plane waves')
add('model.M8', 'not built', None, None, 'PLAN.md P1 (M8 (pair)2 cluster, 890 functions)', 'the (pair)2 cluster of the plan for QTAIM/NCI of the extended network was not built; the label M8 is used in the paper for the five-pair cluster of section 34 (Table 1)', status='missing')
MISSING.append('M8 (pair)2 cluster (PLAN.md P1): not built, so there is no QTAIM/NCI for the anion–anion O4–H4A···O5 bond (2.500 Å) or the other inter-pair H-bonds; the only QTAIM is on the 27-atom ion pair.')

# =====================================================================================================================
section('3. Level of theory and protocol',
        'Settings read from the result JSON provenance blocks, the case files and dft/common.py; not from the draft '
        'text of RESULTS.md, which is stale on the grid (see section 18).')
add('protocol.xc', 'B3LYP (libxc HYB_GGA_XC_B3LYP, VWN-RPA local correlation, Gaussian definition)', None, None, 'dft/common.py:XC and docstring', 'exchange-correlation functional')
add('protocol.dispersion', 'Grimme D3 with Becke–Johnson damping, D3(BJ) (case disp = d3bj); plain B3LYP for the sensitivity case', None, None, 'molecules/sfa/cases/*.json:disp', 'dispersion correction')
add('protocol.basis', '6-311++G(d,p), spherical 5d; def2-TZVP for the *_tz twins', None, None, 'dft/common.py:BASIS; molecules/sfa/cases/*_tz.json:basis', 'basis set')
add('protocol.density_fitting', 'RI-JK density fitting, auxiliary basis chosen by PySCF', None, None, 'scan/scan_summary.json:level.density_fitting; dft/common.py:make_mf', 'density fitting')
add('protocol.grid_level', 4, None, 0, 'scan/scan_summary.json:level.grid_level; dscf/tddft/nlo.json:provenance.grid_level; opt.xyz comment "grid 4"', 'PySCF DFT integration grid level (≈ UltraFine)')
add('protocol.conv_tol', 1e-9, 'Ha', None, 'dft/common.py:make_mf(conv_tol=1e-9)', 'SCF energy convergence threshold')
add('protocol.pcm', 'IEF-PCM, water ε = 78.3553, Lebedev order 29 (302 points per sphere), PySCF default (modified Bondi) radii', None, None,
    'scan_summary.json:level.environment, pcm_lebedev_order; dft/common.py:SOLVENT_EPS', 'continuum model')
add('protocol.eps_water', 78.3553, None, 4, 'dft/common.py:SOLVENT_EPS["water"]', 'dielectric constant of water used')
add('protocol.pcm_lebedev', 29, None, 0, 'scan/scan_summary.json:level.pcm_lebedev_order', 'PCM surface Lebedev order')
add('protocol.pyscf', '2.14.0', None, None, 'dscf.json/nlo.json/tddft json:provenance.pyscf_version', 'PySCF version')
add('protocol.geometric', '1.1.1', None, None, f'{ARCH}/molecules/sfa/results/sfa_pair_pcm_d3/logs/optimize.log', 'geomeTRIC optimiser version')
add('protocol.opt_criteria', 'Gaussian default: max force 4.5e-4, RMS force 3.0e-4 Ha/Bohr (plus displacement criteria); PCM pair and plain-B3LYP pair accepted on negligible forces after the job was killed', None, None,
    'sfa_pair_pcm_d3/opt_summary.json:converged_on; STAGES.md', 'optimisation convergence')
add('protocol.pcm_opt_forces', 'rms 7.35e-5, max 1.85e-4 Ha/Bohr at step 64', None, None, 'sfa_pair_pcm_d3/opt_summary.json:converged_on',
    'final forces of the primary optimisation (exit code −15; accepted on forces, confirmed a minimum by the frequencies)', case='sfa_pair_pcm_d3')
add('protocol.fd_step', 0.005, 'bohr', 3, 'sfa_pair_pcm_d3/freq.json:hessian_method', 'finite-difference step for the Hessian (central differences of analytic gradients)')
add('protocol.fd_gradients', 162, None, 0, 'sfa_pair_pcm_d3/STAGES.md ("162 gradients")', 'gradients in the FD Hessian (2 × 81 coordinates)')
add('protocol.multiwfn', '2026.9.20', None, None, f'{ARCH}/molecules/sfa/results/sfa_pair_pcm_d3/logs/mwfn_mep.out header', 'Multiwfn version')
add('protocol.td', 'TD-DFT (full RPA, not TDA), 20 singlets, CAM-B3LYP/6-311++G(d,p), IEF-PCM water non-equilibrium (ε∞ = 1.78); B3LYP comparison; D3 omitted in the TD step', None, None,
    'sfa_pair_pcm_d3/tddft_camb3lyp_noneq.json:td_method, nstates, environment, notes', 'excited-state protocol')
add('protocol.td_eps_inf', 1.78, None, 2, 'tddft_*_noneq.json:environment', 'optical dielectric constant in the non-equilibrium PCM (PySCF default)')
add('protocol.td_fwhm', 0.3, 'eV', 1, 'tddft_*.json:fwhm_eV', 'Gaussian broadening of the simulated spectra (FWHM)')
add('protocol.nlo', 'finite field on the SCF dipole: ±F, ±2F per axis and (±F, ±F), (±2F, ±2F) per plane, F = 0.001 a.u., one Richardson step; Taylor (T) convention; static, electronic only; PCM field on the solute, static solvent response', None, None,
    'sfa_pair_pcm_d3/nlo.json:method, field_au, beta_convention', 'NLO protocol')
add('protocol.nlo_field', 0.001, 'a.u.', 3, 'nlo.json:field_au[0]', 'base finite-field strength')
add('protocol.nlo_field2', 0.002, 'a.u.', 3, 'nlo.json:field_au[1]', 'doubled field strength used for the Richardson step')
add('protocol.dscf', 'vertical ΔSCF: RKS for N, UKS doublets for N ± 1 at the N-electron geometry; equilibrium PCM for every charge state; stability checked', None, None,
    'sfa_pair_pcm_d3/dscf.json:method, N_minus_1.stability_checked', 'ΔSCF protocol (equilibrium, not non-equilibrium, solvation of the ions: a caveat)')
add('protocol.qe', 'Quantum ESPRESSO 7.5 pw.x, PBE + D3(BJ) (two-body), SSSP PBE efficiency 1.3.0, cell fixed at the X-ray (294 K) values, BFGS with forc_conv_thr 5e-4 Ry/Bohr', None, None,
    'periodic_qe/periodic_qe.json:method', 'periodic protocol')
add('protocol.qe_ecut', '70/560 Ry (wavefunction/density)', None, None, 'periodic_qe/convergence.json:chosen.ecut_Ry', 'plane-wave cutoffs (SSSP floor 60/480 Ry)')
add('protocol.qe_kgrid', '2×1×1 (2 irreducible k-points)', None, None, 'periodic_qe/convergence.json:chosen.k', 'Monkhorst–Pack k-point grid')
add('protocol.docking', 'AutoDock Vina 1.2.7 (Vina and Vinardo scoring), AutoDock4Zn (AutoGrid 4.2.7 maps); exhaustiveness 32; 10 seeds (20 for the sulfanilamide AD4Zn pose jobs); 20 Å box; Meeko 0.8.0, RDKit 2026.03.1, PDBFixer/OpenMM 8.6.1, PROPKA 3.5.1, PoseBusters 0.6.5, ProLIF 2.2.2, PLIP, spyrmsd 0.9.0, gemmi 0.7.5', None, None,
    'dock/VALIDATION.md (Environment), dock/RESULTS.md (What ran)', 'docking software and settings')
add('protocol.dock_exhaustiveness', 32, None, 0, 'dock/*/runs/*/results.json:exhaustiveness', 'Vina/AD4Zn exhaustiveness')
add('protocol.dock_box', 20.0, 'A', 1, 'dock/*/runs/*/results.json:box.size', 'cubic docking box edge (AD4Zn maps 20.25 Å; shown to make no difference)')
add('protocol.dock_seeds', 10, None, 0, 'dock/results/production_summary.json:jobs.*.n_seeds', 'seeds per job (20 for sulfanilamide AD4Zn into 6RL9 and 3HS4)')
add('protocol.dock_ph', 7.4, None, 1, 'dock/VALIDATION.md; dock/prep_ligands.py', 'pH for ligand species and receptor protonation (PROPKA)')

# =====================================================================================================================
section('4. Ion-pair geometries and the comparison with the X-ray structure',
        'All bond statistics recomputed (derive/d03_geometry.py) over the 26 bonds perceived in the X-ray pair; '
        'classes as below. X–H bonds are compared for completeness only: the X-ray C–H are riding at 0.93 Å and the '
        'O–H restrained to 0.82 Å, so their deviations measure the X-ray model, not DFT. In the gas phase H1A has '
        'moved to O5, so N1–H1A is excluded there. "periodic" = relax_i of the PBE-D3(BJ) cell, formula unit 1.')
MODELS = [('pcm', 'sfa_pair_pcm_d3', 'IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p)'),
          ('gas', 'sfa_pair_gas_d3', 'gas phase, B3LYP-D3(BJ)/6-311++G(d,p)'),
          ('sens.b3lyp', 'sfa_pair_pcm_b3lyp', 'IEF-PCM, B3LYP without dispersion'),
          ('sens.tz', 'sfa_pair_pcm_d3_tz', 'IEF-PCM, B3LYP-D3(BJ)/def2-TZVP'),
          ('periodic', 'periodic_qe (relax_i)', 'periodic PBE-D3(BJ)')]
MKEY = {'pcm': 'pcm', 'gas': 'gas', 'sens.b3lyp': 'pcm_b3lyp', 'sens.tz': 'pcm_tz', 'periodic': 'periodic_i'}
CLS = [('C-C (ring)', 'CC'), ('C1-N1', 'C1N1'), ('C4-S1', 'C4S1'), ('S1-O', 'SO'), ('S1-N2', 'S1N2'),
       ('P1-O(H) (O3, O4)', 'POH'), ('P1-O (O5, O6)', 'PO'), ('_heavy', 'heavy'), ('N-H', 'NH'), ('C-H', 'CH'), ('O-H', 'OH'), ('_all', 'all')]
opt_summ = {c: J(RES / c / 'opt_summary.json') for c in ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'sfa_pair_pcm_b3lyp', 'sfa_pair_pcm_d3_tz']}
for pre, case, lab in MODELS:
    m = D3['models'][MKEY[pre]]
    src = f"derived/d03_geometry.json:models.{MKEY[pre]}"
    for cname, ck in CLS:
        s = m['bond_stats'].get(cname)
        if not s:
            continue
        add(f'{pre}.geom.{ck}.mae', s['mae'], 'A', 3, f'{src}.bond_stats[{cname}]', f'{lab}: MAE of {cname.strip("_")} bond lengths vs X-ray (n = {s["n"]})', case=case)
        add(f'{pre}.geom.{ck}.rmsd', s['rmsd'], 'A', 3, f'{src}.bond_stats[{cname}]', f'{lab}: RMSD of {cname.strip("_")} bond lengths vs X-ray (n = {s["n"]})', case=case)
        add(f'{pre}.geom.{ck}.maxdev', s['max_abs'], 'A', 3, f'{src}.bond_stats[{cname}]', f'{lab}: largest |deviation| in {cname.strip("_")} ({s["worst"]})', case=case, note=f'worst bond {s["worst"]}')
    a = m['angle_stats_heavy']
    add(f'{pre}.geom.ang.mae', a['mae'], 'deg', 1, f'{src}.angle_stats_heavy', f'{lab}: MAE of the 22 heavy-atom bond angles vs X-ray', case=case)
    add(f'{pre}.geom.ang.rmsd', a['rmsd'], 'deg', 1, f'{src}.angle_stats_heavy', f'{lab}: RMSD of the heavy-atom angles', case=case)
    add(f'{pre}.geom.ang.maxdev', a['max_abs'], 'deg', 1, f'{src}.angle_stats_heavy', f'{lab}: largest angle deviation ({a["worst"]})', case=case)
    if 'kabsch_rmsd_heavy_pair' in m:
        add(f'{pre}.geom.kabsch_pair', m['kabsch_rmsd_heavy_pair'], 'A', 2, f'{src}.kabsch_rmsd_heavy_pair', f'{lab}: heavy-atom RMSD of the whole pair after superposition (dominated by the anion reorientation)', case=case)
        add(f'{pre}.geom.kabsch_cation', m['kabsch_rmsd_heavy_cation'], 'A', 3, f'{src}.kabsch_rmsd_heavy_cation', f'{lab}: heavy-atom RMSD of the cation alone after superposition', case=case)
        add(f'{pre}.geom.kabsch_anion', m['kabsch_rmsd_heavy_anion'], 'A', 3, f'{src}.kabsch_rmsd_heavy_anion', f'{lab}: heavy-atom RMSD of the anion alone after superposition', case=case)
    if case in opt_summ:
        o = opt_summ[case]
        add(f'{pre}.opt.E', o['E_final_Ha'], 'Ha', 6, f'molecules/sfa/results/{case}/opt_summary.json:E_final_Ha', f'{lab}: final electronic energy (incl. D3 where used)', case=case)
        add(f'{pre}.opt.steps', o['steps'], None, 0, f'molecules/sfa/results/{case}/opt_summary.json:steps', f'{lab}: optimisation steps (after a 6-31G(d) pre-optimisation, except TZ)', case=case)
        add(f'{pre}.opt.bond_rmsd_all', o['bond_rmsd_A'], 'A', 3, f'molecules/sfa/results/{case}/opt_summary.json:bond_rmsd_A',
            f'{lab}: RMSD over all 26 bonds incl. X–H, as reported by the pipeline', case=case,
            check=m['bond_stats']['_all']['rmsd'] if pre != 'gas' else None, tol=0.0005)
for pre, case, lab in MODELS[:4]:
    rows_ = {}
    for l in (RES / case / 'geometry_vs_start.csv').read_text().splitlines()[1:]:
        a_, b_ = l.split(',')[0].replace('#1_554', '').split('-')
        rows_[frozenset((a_, b_))] = (float(l.split(',')[1]), float(l.split(',')[2]))
    mine = D3['models'][MKEY[pre]]['bond_lengths']
    xr = D3['models']['xray']['bond_lengths']
    diffs = [abs(mine[b] - rows_[frozenset(b.split('-'))][1]) for b in mine if frozenset(b.split('-')) in rows_]
    dx = [abs(xr[b] - rows_[frozenset(b.split('-'))][0]) for b in xr if frozenset(b.split('-')) in rows_]
    add(f'{pre}.geom.crosscheck', max(diffs), 'A', 5, f'molecules/sfa/results/{case}/geometry_vs_start.csv vs derived/d03_geometry.json',
        f'{lab}: largest |difference| between our bond lengths and the pipeline CSV over {len(diffs)} bonds (X-ray side: {max(dx):.1e} Å)', case=case)
import math as _m
for pre, case, lab in MODELS:
    bl_ = D3['models'][MKEY[pre]]['bond_lengths']
    xr_ = D3['models']['xray']['bond_lengths']
    dev_ = [bl_[b] - xr_[b] for b in ['C1-N1', 'C4-S1']]
    add(f'{pre}.geom.CN_CS.mae', sum(abs(x) for x in dev_) / 2, 'A', 3, f"derived/d03_geometry.json:models.{MKEY[pre]}.bond_lengths (C1-N1, C4-S1)",
        f'{lab}: MAE of C1–N1 and C4–S1 together (Table 3 class)', case=case)
    add(f'{pre}.geom.CN_CS.rmsd', _m.sqrt(sum(x * x for x in dev_) / 2), 'A', 3, f"derived/d03_geometry.json:models.{MKEY[pre]}.bond_lengths (C1-N1, C4-S1)",
        f'{lab}: RMSD of C1–N1 and C4–S1 together (Table 3 class)', case=case)
add('pcm.opt.accepted', 'negligible forces (rms 7.35e-05, max 1.85e-04 Ha/Bohr) after the job was killed at step 64', None, None,
    'sfa_pair_pcm_d3/opt_summary.json:converged_on', 'how the primary optimisation ended', case='sfa_pair_pcm_d3')
add('sens.b3lyp.opt.accepted', 'negligible forces (rms 5.92e-05, max 1.21e-04 Ha/Bohr) after the job was killed at step 68; no frequencies', None, None,
    'sfa_pair_pcm_b3lyp/opt_summary.json:converged_on', 'how the plain-B3LYP optimisation ended', case='sfa_pair_pcm_b3lyp')
# the two intra-pair contacts in every model
for pre, case, lab in MODELS:
    pm = D3['models'][MKEY[pre]]['proton']
    add(f'{pre}.hb.N1_O5', pm['d_N1_O5'], 'A', 3, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton.d_N1_O5", f'{lab}: N1···O5', case=case)
    add(f'{pre}.hb.C6_O4', pm['d_C6_O4'], 'A', 3, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton.d_C6_O4", f'{lab}: C6···O4 (R²₂(8) second leg)', case=case)
    add(f'{pre}.hb.H6_O4', pm['d_H6_O4'], 'A', 3, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton.d_H6_O4", f'{lab}: H6···O4', case=case)
    if 'ang_N1_H1A_O5' in pm:
        add(f'{pre}.hb.ang_N1H1AO5', pm['ang_N1_H1A_O5'], 'deg', 1, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton", f'{lab}: N1–H1A···O5 angle (O5–H1A···N1 in the gas phase)', case=case)
        add(f'{pre}.hb.ang_C6H6O4', pm['ang_C6_H6_O4'], 'deg', 1, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton", f'{lab}: C6–H6···O4 angle', case=case)
        add(f'{pre}.tau_C3C4S1N2', pm['tau_C3_C4_S1_N2'], 'deg', 1, f"derived/d03_geometry.json:models.{MKEY[pre]}.proton", f'{lab}: C3–C4–S1–N2 torsion', case=case)
CAVEATS.append('Plain B3LYP (no D3) breaks the second leg of the R²₂(8) ring: C6···O4 opens from 3.383 Å (D3) to 4.798 Å, while N1···O5 and δ are unchanged (2.584 Å, −0.387 Å). Dispersion holds the ring closed; the proton position does not depend on it.')
CAVEATS.append('The whole-pair heavy-atom RMSD after superposition (0.54 Å in IEF-PCM) is dominated by the anion rotating about the N1–H1A···O5 axis; each ion alone fits to 0.12 Å (cation) and 0.04 Å (anion).')

# =====================================================================================================================
section('5. Proton metric δ in every environment',
        'δ = d(N1–H1A) − d(H1A···O5), recomputed from each .xyz (derive/d03_geometry.py, d08_periodic.py). Negative δ: '
        'proton on N1 (salt); positive: on O5 (neutral pair).')
ENV = [('xray', 'xray', 'X-ray (deposited H1A)', None),
       ('gas', 'gas', 'gas phase, ε = 1 (started from the X-ray geometry)', 'sfa_pair_gas_d3'),
       ('eps4', 'eps_4', 'IEF-PCM ε = 4 (from the water optimum)', 'sfa_pair_pcm_d3'),
       ('eps8', 'eps_8', 'IEF-PCM ε = 8', 'sfa_pair_pcm_d3'),
       ('eps20', 'eps_20', 'IEF-PCM ε = 20', 'sfa_pair_pcm_d3'),
       ('pcm', 'pcm', 'IEF-PCM water ε = 78.36 (primary optimum)', 'sfa_pair_pcm_d3'),
       ('scanN', 'scan_A_anchor', 'IEF-PCM water, scan anchor A (re-optimised N-side minimum)', 'sfa_pair_pcm_d3'),
       ('scanO', 'scan_B_anchor', 'IEF-PCM water, scan anchor B (O-side minimum)', 'sfa_pair_pcm_d3'),
       ('sens.b3lyp', 'pcm_b3lyp', 'IEF-PCM water, B3LYP without D3', 'sfa_pair_pcm_b3lyp'),
       ('sens.tz', 'pcm_tz', 'IEF-PCM water, def2-TZVP', 'sfa_pair_pcm_d3_tz'),
       ('periodic', 'periodic_i', 'periodic PBE-D3(BJ), relax from the X-ray (i)', 'periodic_qe'),
       ('periodic.ii', 'periodic_ii', 'periodic PBE-D3(BJ), relax from H1A on O5 (ii)', 'periodic_qe'),
       ('periodic.iifix', 'periodic_iifix', 'periodic, H1A fixed at 0.98 Å from the X-ray O5 (iifix; WITHDRAWN, see section 17)', 'periodic_qe'),
       ('periodic.iirfix', 'periodic_iirfix', 'periodic, constrained neutral: H1A placed 0.98 Å from O5 in relaxed cell (i), H1A and O5 held in space (iirfix)', 'periodic_qe')]
ladder_json = J(RES / 'sfa_pair_pcm_d3/eps_ladder/eps_ladder.json')
LADM = {r_['eps']: r_['metrics'] for r_ in ladder_json['rows']}
SCSV = {(l.split(',')[0], l.split(',')[1]): l.split(',') for l in (RES / 'sfa_pair_pcm_d3/scan/scan_profile.csv').read_text().splitlines()[1:]}
PQ = J(RES / 'periodic_qe/periodic_qe.json')


def gvs(case, bond):
    for l in (RES / case / 'geometry_vs_start.csv').read_text().splitlines()[1:]:
        if l.split(',')[0] == bond:
            return float(l.split(',')[2])


def ladder_src(e, label):
    m = LADM[e]
    return ({'d_N1_H1A': m['d_N1_H1A_A'], 'd_H1A_O5': m['d_H1A_O5_A'], 'delta': m['delta_A'], 'd_N1_O5': m['d_N1_O5_A'],
             'P1-O5': m['P_O_A']['P1-O5']}, f'molecules/sfa/results/sfa_pair_pcm_d3/eps_ladder/eps_ladder.json:rows[eps={label}].metrics')


def scan_src(chain):
    r_ = SCSV[(chain, '0')]
    return ({'d_N1_H1A': float(r_[4]), 'd_H1A_O5': float(r_[5]), 'delta': float(r_[6]), 'd_N1_O5': float(r_[7]), 'P1-O5': float(r_[9])},
            f'molecules/sfa/results/sfa_pair_pcm_d3/scan/scan_profile.csv:{chain},0')


def pq_src(k):
    # periodic_qe.json stores 4 decimals; print from our full-precision re-analysis of relax_<k>.extxyz and check
    # against the JSON (avoids double rounding)
    fu = D8['structures'][f'relax_{k}']['fu1']
    g = PQ['relax'][k]['geometry']
    mk = {'i': 'periodic_i', 'ii': 'periodic_ii', 'iifix': 'periodic_iifix', 'iirfix': 'periodic_iirfix'}[k]
    pm = D3['models'][mk]['proton']
    return ({'d_N1_H1A': fu['d_N1_H1A'], 'd_H1A_O5': fu['d_H1A_O5'], 'delta': fu['delta'], 'd_N1_O5': fu['d_N1_O5'],
             'P1-O5': pm['P1-O5'], 'C1-N1': pm['C1-N1']},
            f'molecules/sfa/results/periodic_qe/relax_{k}.extxyz via derived/d08_periodic.json:structures.relax_{k}.fu1 (JSON copy: periodic_qe.json:relax.{k}.geometry.proton, 4 decimals)')


def pc_src(case):
    pc = J(RES / case / 'proton_check.json')
    o = J(RES / case / 'opt_summary.json')
    return ({'d_N1_H1A': pc['d_H_N_A'], 'd_H1A_O5': pc['d_H_O_A'], 'delta': pc['d_H_N_A'] - pc['d_H_O_A'],
             'd_N1_O5': o['distances_A']['N1...O5 (X-ray 2.760 A)']['dft'], 'P1-O5': gvs(case, 'P1#1_554-O5#1_554'),
             'C1-N1': gvs(case, 'C1-N1')}, f'molecules/sfa/results/{case}/proton_check.json, opt_summary.json:distances_A, geometry_vs_start.csv')


XP = D1['pair']
SRC = {'xray': ({'d_N1_H1A': XP['d_N1_H1A'], 'd_H1A_O5': XP['d_H1A_O5'], 'delta': XP['delta'], 'd_N1_O5': XP['d_N1_O5'],
                 'P1-O5': XP['P1-O5'], 'C1-N1': XP['C1-N1']}, 'derived/d01_xray.json:pair (CIF _atom_site coordinates; also periodic_qe.json:xray.proton, 4 decimals)'),
       'gas': (dict(ladder_src(1.0, 1)[0], **{'C1-N1': gvs('sfa_pair_gas_d3', 'C1-N1')}), ladder_src(1.0, 1)[1] + '; sfa_pair_gas_d3/geometry_vs_start.csv'),
       'eps4': ladder_src(4.0, 4), 'eps8': ladder_src(8.0, 8), 'eps20': ladder_src(20.0, 20),
       'pcm': pc_src('sfa_pair_pcm_d3'), 'scanN': scan_src('A'), 'scanO': scan_src('B'),
       'sens.b3lyp': pc_src('sfa_pair_pcm_b3lyp'), 'sens.tz': pc_src('sfa_pair_pcm_d3_tz'),
       'periodic': pq_src('i'), 'periodic.ii': pq_src('ii'), 'periodic.iifix': pq_src('iifix'), 'periodic.iirfix': pq_src('iirfix')}
for pre, mk, lab, case in ENV:
    pm = D3['models'][mk]['proton'] if mk in D3['models'] else D3['proton_extra'][mk]
    dsrc = f"derived/d03_geometry.json:{'models.' + mk + '.proton' if mk in D3['models'] else 'proton_extra.' + mk}"
    sv, ssrc = SRC[pre]
    for q, u, p, what in [('d_N1_H1A', 'A', 3, 'N1–H1A'), ('d_H1A_O5', 'A', 3, 'H1A···O5'), ('delta', 'A', 3, 'δ'),
                          ('d_N1_O5', 'A', 3, 'N1···O5'), ('P1-O5', 'A', 3, 'P1–O5'), ('C1-N1', 'A', 3, 'C1–N1')]:
        k = f'proton.{pre}.{q.replace("-", "")}'
        if q in sv and sv[q] is not None and pre.startswith('periodic'):
            kk = {'periodic': 'i', 'periodic.ii': 'ii', 'periodic.iifix': 'iifix', 'periodic.iirfix': 'iirfix'}[pre]
            g_ = PQ['relax'][kk]['geometry']
            jv = {'d_N1_H1A': g_['proton']['d_N1_H1A_A'], 'd_H1A_O5': g_['proton']['d_H1A_O5_A'], 'delta': g_['proton']['delta_A'],
                  'd_N1_O5': g_['proton']['d_N1_O5_A'], 'P1-O5': g_['proton']['P_O_A']['P1-O5'], 'C1-N1': g_['bonds']['C1-N1']}[q]
            add(k, sv[q], u, p, ssrc, f'{lab}: {what}', case=case, check=jv, tol=0.0001)
        elif q in sv and sv[q] is not None and pre == 'xray':
            jv = {'d_N1_H1A': PQ['xray']['proton']['d_N1_H1A_A'], 'd_H1A_O5': PQ['xray']['proton']['d_H1A_O5_A'], 'delta': PQ['xray']['proton']['delta_A'],
                  'd_N1_O5': PQ['xray']['proton']['d_N1_O5_A'], 'P1-O5': PQ['xray']['proton']['P_O_A']['P1-O5'], 'C1-N1': PQ['xray']['bonds']['C1-N1']}[q]
            add(k, sv[q], u, p, ssrc, f'{lab}: {what}', case=case, check=jv, tol=0.0001)
        elif q in sv and sv[q] is not None:
            add(k, sv[q], u, p, ssrc, f'{lab}: {what}', case=case, check=pm[q], tol=0.0006)
        else:
            add(k, pm[q], u, p, dsrc, f'{lab}: {what}', case=case, status='own derivation')
    add(f'proton.{pre}.side', pm['proton_on'], None, None, dsrc, f'{lab}: proton on N1 or O5')
lad = {r_['eps']: r_ for r_ in ladder_json['rows']}
for e, pre in [(4.0, 'eps4'), (8.0, 'eps8'), (20.0, 'eps20')]:
    add(f'proton.{pre}.steps', lad[e]['steps'], None, 0, f'sfa_pair_pcm_d3/eps_ladder/eps_{int(e)}_N/result.json:steps', f'IEF-PCM ε = {int(e)}: optimisation steps (converged on forces + energy)', case='sfa_pair_pcm_d3',
          check=None)
add('eps.bracket', '4 < ε_crossover < 8', None, None, 'sfa_pair_pcm_d3/eps_ladder/eps_ladder.json:crossover', 'the N-protonated minimum disappears between ε = 4 (proton on O5) and ε = 8 (kept on N1)')
add('eps.values', [1, 4, 8, 20, 78.3553], None, None, 'eps_ladder.csv', 'dielectric constants of the ladder')
CAVEATS.append('Dielectric ladder: every ε > 1 point started from the water optimum with the proton on N1; ε = 1 is the gas optimisation from the X-ray structure. At ε = 4 the proton moved to O5 (80 steps). Only the side of the proton is compared across ε; total energies at different ε are not comparable. The task statement "ε 1 O-side, 8/20/78 N-side" is right but omits ε = 4, which is also O-side.')
# separated-species proton-transfer energetics (own derivation from the fragment energies)
E = {c: J(RES / c / 'opt_summary.json')['E_final_Ha'] for c in ['cation_pcm_d3', 'h2po4_pcm_d3', 'sulfanilamide_pcm_d3', 'h3po4_pcm_d3',
                                                                'cation_gas_d3', 'h2po4_gas_d3', 'sulfanilamide_gas_d3', 'h3po4_gas_d3']}
dpcm = (E['sulfanilamide_pcm_d3'] + E['h3po4_pcm_d3'] - E['cation_pcm_d3'] - E['h2po4_pcm_d3']) * HA_KCAL
dgas = (E['sulfanilamide_gas_d3'] + E['h3po4_gas_d3'] - E['cation_gas_d3'] - E['h2po4_gas_d3']) * HA_KCAL
add('pt.sep.pcm', dpcm, 'kcal mol-1', 2, 'opt_summary.json:E_final_Ha of cation/h2po4/sulfanilamide/h3po4 _pcm_d3',
    'IEF-PCM: E(sulfanilamide + H3PO4) − E(cation + H2PO4−), separated, relaxed species, electronic only (negative = neutral side lower)', note='own derivation; no ZPE, no thermal terms')
add('pt.sep.gas', dgas, 'kcal mol-1', 1, 'opt_summary.json:E_final_Ha of the _gas_d3 fragments',
    'gas phase: the same separated-species proton-transfer energy (dominated by the ion-pair charge separation)', note='own derivation')

# =====================================================================================================================
section('6. Proton scan along N1–H1A in IEF-PCM (water)',
        'Relaxed scan, B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM water, two chains from opposite ends, 0.10 Å grid. '
        'Energies relative to anchor A (re-optimised N-side minimum), recomputed from each point JSON; distances '
        'recomputed from each point .xyz (derive/d04_scan.py). Electronic energies, no ZPE.')
ss = J(RES / 'sfa_pair_pcm_d3/scan/scan_summary.json')
for p in D4['points']:
    k = f"scan.{p['chain']}{p['index']:02d}"
    lab = f"chain {p['chain']} point {p['index']} ({p['kind']}{', target ' + format(p['target'], '.2f') + ' Å' if p['target'] else ''})"
    src = f"molecules/sfa/results/sfa_pair_pcm_d3/scan/{Path(p['file']).name}"
    row = SCSV[(p['chain'], str(p['index']))]
    csvs = f"molecules/sfa/results/sfa_pair_pcm_d3/scan/scan_profile.csv:{p['chain']},{p['index']}"
    add(k + '.d_NH', float(row[4]), 'A', 3, csvs + ' d_N1_H1A_A', f'scan {lab}: d(N1–H1A)', case='sfa_pair_pcm_d3', check=p['d_N1_H1A'], tol=0.0006)
    add(k + '.delta', float(row[6]), 'A', 3, csvs + ' delta_A', f'scan {lab}: δ', case='sfa_pair_pcm_d3', check=p['delta'], tol=0.0006)
    add(k + '.d_NO', float(row[7]), 'A', 3, csvs + ' d_N1_O5_A', f'scan {lab}: N1···O5', case='sfa_pair_pcm_d3', check=p['d_N1_O5'], tol=0.0006)
    add(k + '.Erel', float(row[11]), 'kcal mol-1', 2, csvs + ' E_rel_kcal', f'scan {lab}: E − E(anchor A)', case='sfa_pair_pcm_d3', check=p['E_rel_kcal'], tol=0.001)
    add(k + '.E', p['E_Ha'], 'Ha', 8, src + ':E_Ha', f'scan {lab}: total energy', case='sfa_pair_pcm_d3')
    add(k + '.converged', p['converged'], None, None, src + ':converged', f'scan {lab}: converged flag ({p["steps"]} steps)', case='sfa_pair_pcm_d3')
csv_rows = {(l.split(',')[0], int(l.split(',')[1])): l.split(',') for l in (RES / 'sfa_pair_pcm_d3/scan/scan_profile.csv').read_text().splitlines()[1:]}
env = {(e['chain'], e['index']): e for e in ss['lower_envelope']}
add('scan.barrier', ss['envelope_maxima'][0]['E_rel_kcal'], 'kcal mol-1', 2, 'scan/scan_summary.json:envelope_maxima[0].E_rel_kcal',
    'electronic barrier N→O on the lower envelope (point A_02, d = 1.25 Å)', case='sfa_pair_pcm_d3', check=D4['barrier_from_N_kcal'], tol=1e-4)
add('scan.barrier_kJ', D4['barrier_from_N_kJ'], 'kJ mol-1', 2, 'derived/d04_scan.json:barrier_from_N_kJ', 'same barrier in kJ mol⁻¹', case='sfa_pair_pcm_d3')
add('scan.barrier_spline', J(RES / 'sfa_pair_pcm_d3/scan/zpe_barrier.json')['profile']['barrier_cubic_spline_kcal'], 'kcal mol-1', 2,
    'scan/zpe_barrier.json:profile.barrier_cubic_spline_kcal', 'barrier from a cubic spline through the envelope (at d = 1.277 Å)', case='sfa_pair_pcm_d3')
add('scan.Omin', ss['envelope_minima'][1]['E_rel_kcal'], 'kcal mol-1', 2, 'scan/scan_summary.json:envelope_minima[1]',
    'O-protonated (neutral) pair relative to the salt minimum (anchor B, d = 1.631 Å): negative = neutral pair lower', case='sfa_pair_pcm_d3', check=D4['O_min_rel_kcal'], tol=1e-4)
add('scan.Omin_kJ', D4['O_min_rel_kJ'], 'kJ mol-1', 2, 'derived/d04_scan.json:O_min_rel_kJ', 'same in kJ mol⁻¹', case='sfa_pair_pcm_d3')
add('scan.reverse_barrier', D4['reverse_barrier_kcal'], 'kcal mol-1', 2, 'scan/zpe_barrier.json:profile.reverse_barrier_kcal', 'barrier O→N (top minus O-side minimum)', case='sfa_pair_pcm_d3',
    check=J(RES / 'sfa_pair_pcm_d3/scan/zpe_barrier.json')['profile']['reverse_barrier_kcal'], tol=1e-4)
add('scan.hysteresis_1.35', ss['hysteresis_A_minus_B_kcal']['1.350'], 'kcal mol-1', 3, 'scan/scan_summary.json:hysteresis_A_minus_B_kcal["1.350"]',
    'E(chain A) − E(chain B) at d = 1.35 Å, the one grid point both chains reached', case='sfa_pair_pcm_d3', check=D4['hysteresis']['1.350']['A_minus_B_kcal'], tol=1e-4)
add('scan.hysteresis_1.35_dNO', D4['hysteresis']['1.350']['dNO_A'] - D4['hysteresis']['1.350']['dNO_B'], 'A', 4, 'derived/d04_scan.json:hysteresis',
    'N1···O5 difference between the two chains at d = 1.35 Å (same structure to 2e-4 Å)', case='sfa_pair_pcm_d3')
add('scan.n_points', D4['n_points'], None, 0, 'scan/scan_summary.json:n_points', 'scan points (4 per chain incl. anchors)')
add('scan.n_unconverged', D4['n_unconverged'], None, 0, 'scan/scan_summary.json:n_unconverged', 'points stopped at the step limit without the convergence flag (A_02 at the barrier top; its forces are below the thresholds)')
a02 = [p for p in D4['points'] if p['chain'] == 'A' and p['index'] == 2][0]
add('scan.A02_gmax_proj', a02['last_gmax_conv'], 'Ha/bohr', 5, 'scan/A_02_d1.250.json:history[-1].gmax_conv',
    'barrier-top point: max constrained-projected gradient at the last of 61 steps (below the 4.5e-4 threshold) — the flag is false but the point is essentially relaxed', case='sfa_pair_pcm_d3')
add('scan.A02_last_dE', a02['last_dE_Ha'] * HA_KCAL, 'kcal mol-1', 4, 'scan/A_02_d1.250.json:history[-2:]', 'energy change over the last step of A_02', case='sfa_pair_pcm_d3')
CAVEATS.append('The barrier top (A_02, d = 1.25 Å) is flagged unconverged after 61 steps, but its constrained-projected gradient (max 1.6e-4, rms 7.6e-5 Ha/Bohr) is below the thresholds and the last step lowered E by 0.0008 kcal/mol; the 0.72 kcal/mol barrier is therefore reliable to ~0.01 kcal/mol at this grid, but it is a discrete maximum on a 0.10 Å grid (spline 0.75 kcal/mol).')
CAVEATS.append('Hysteresis: the chains overlap only at d = 1.35 Å; there chain A lies 0.014 kcal/mol above chain B with the same N1···O5 to 2e-4 Å, i.e. no path dependence at this resolution. Chain A at 1.35 Å is therefore off the lower envelope.')

# =====================================================================================================================
section('7. Functional check at the scan end points (fsp_scan.json, five methods)',
        'Single points on the B3LYP-D3(BJ) geometries of A_00 (salt, d = 1.098 Å), A_02 (d = 1.250 Å) and B_00 (neutral, '
        'd = 1.631 Å), IEF-PCM water, 6-311++G(d,p): B3LYP-D3(BJ), ωB97X-D, PBE0-D3(BJ), M06-2X and DF-MP2 (frozen core, '
        'aug-cc-pVTZ-RI fitting basis; PCM as a perturbation, reaction field frozen at the RHF density). Read from '
        'sfa_pair_pcm_d3/fsp_scan.json (burst 2, 28 Sep 2026, all 15 points recomputed). Re-derived (derive/d11_burst2.py): '
        'differences recomputed from each E_Ha; the three geometry files re-hashed (sha256 as recorded) and re-measured; the '
        'four DFT energies reproduce the burst-1 run, whose JSON was lost and whose log had been parsed by derive/d04_scan.py, '
        'to < 1e-10 Ha.')
FSPJ = J(RES / 'sfa_pair_pcm_d3/fsp_scan.json')
FSPS = 'molecules/sfa/results/sfa_pair_pcm_d3/fsp_scan.json'
D11 = J(DER / 'd11_burst2.json')
FC = D11['fcheck']
MN = {'b3lyp-d3bj': ('b3lyp', 'B3LYP-D3(BJ)'), 'wb97x-d': ('wb97xd', 'ωB97X-D'), 'pbe0-d3bj': ('pbe0', 'PBE0-D3(BJ)'),
      'm06-2x': ('m062x', 'M06-2X'), 'df-mp2': ('mp2', 'DF-MP2')}
GN = {'end_N': 'end_N_d1.098', 'top1': 'top1_d1.250', 'end_O': 'end_O_d1.631'}
for jm, (mk, meth) in MN.items():
    rel_ = FSPJ['relative'][jm]['vs_first_kcal']
    der = FC['relative'][jm]
    for g, gl in [('end_N', 'salt end point'), ('top1', 'd = 1.25 Å point'), ('end_O', 'neutral end point')]:
        add(f'fcheck.{mk}.E.{g}', FSPJ['energies'][jm][GN[g]]['E_Ha'], 'Ha', 8, f'{FSPS}:energies["{jm}"]["{GN[g]}"].E_Ha',
            f'{meth}: total energy at the {gl}', case='sfa_pair_pcm_d3')
    add(f'fcheck.{mk}.dE_O_N', rel_['end_O_d1.631'], 'kcal mol-1', 2, f'{FSPS}:relative["{jm}"].vs_first_kcal["end_O_d1.631"]',
        f'{meth}: E(neutral O5–H1A end point) − E(salt end point); negative = neutral lower', case='sfa_pair_pcm_d3', check=der['dE_O_N'], tol=1e-6)
    add(f'fcheck.{mk}.top_N', rel_['top1_d1.250'], 'kcal mol-1', 2, f'{FSPS}:relative["{jm}"].vs_first_kcal["top1_d1.250"]',
        f'{meth}: E(d = 1.25 Å point) − E(salt end point) (barrier estimate on this path; negative = no barrier)', case='sfa_pair_pcm_d3',
        check=der['top_N'], tol=1e-6)
    add(f'fcheck.{mk}.top_O', rel_['top1_d1.250'] - rel_['end_O_d1.631'], 'kcal mol-1', 2,
        f'{FSPS}:relative["{jm}"].vs_first_kcal (top1 − end_O)', f'{meth}: E(d = 1.25 Å point) − E(neutral end point)', case='sfa_pair_pcm_d3',
        check=der['top_O'], tol=1e-6)
    add(f'fcheck.{mk}.dE_O_N_kJ', FSPJ['relative'][jm]['vs_first_kJ']['end_O_d1.631'], 'kJ mol-1', 2, f'{FSPS}:relative["{jm}"].vs_first_kJ["end_O_d1.631"]',
        f'{meth}: the same neutral − salt difference in kJ mol⁻¹', case='sfa_pair_pcm_d3', check=der['dE_O_N_kJ'], tol=1e-5)
    if 'max_abs_E_vs_burst1_log_Ha' in der:
        add(f'fcheck.{mk}.vs_burst1_Ha', der['max_abs_E_vs_burst1_log_Ha'], 'Ha', 11, 'derived/d11_burst2.json:fcheck.relative (fsp_scan.json vs the burst-1 log parse in d04_scan.json)',
            f'{meth}: largest |ΔE| between the burst-2 run and the burst-1 run on the same three geometries (reproducibility)', case='sfa_pair_pcm_d3')
add('fcheck.wb97xd.disp_note', 'ωB97X-D run as libxc HYB_GGA_XC_WB97X_D plus a separately added Chai–Head-Gordon damped C6 term (−0.0163 to −0.0171 Ha at the three geometries); E_Ha includes it', None, None,
    f'{FSPS}:energies["wb97x-d"][*].E_disp_CHG_Ha, notes; dft/fsp.py docstring', 'how ωB97X-D was evaluated')
add('fcheck.mp2.method', 'DF-MP2, frozen core (24 orbitals), aug-cc-pVTZ-RI fitting basis, on the RHF/6-311++G(d,p) reference in IEF-PCM water; PCM enters as a perturbation (reaction field frozen at the RHF density, PTE)', None, None,
    f'{FSPS}:energies["df-mp2"][*].n_frozen, aux_mp2; notes', 'DF-MP2 settings of the functional check')
add('fcheck.mp2.n_frozen', FSPJ['energies']['df-mp2']['end_N_d1.098']['n_frozen'], None, 0, f'{FSPS}:energies["df-mp2"][*].n_frozen',
    'frozen core orbitals in DF-MP2 (1s of C, N, O; 1s2s2p of P and S)', check=FC['frozen_core_expected'], tol=0)
for g in ('end_N', 'top1', 'end_O'):
    e_ = FSPJ['energies']['df-mp2'][GN[g]]
    add(f'fcheck.mp2.Ecorr.{g}', e_['E_corr_Ha'], 'Ha', 6, f'{FSPS}:energies["df-mp2"]["{GN[g]}"].E_corr_Ha', f'DF-MP2 correlation energy at the {g} geometry', case='sfa_pair_pcm_d3')
add('fcheck.vs_burst1_max_Ha', max(r_['max_abs_E_vs_burst1_log_Ha'] for r_ in FC['relative'].values() if 'max_abs_E_vs_burst1_log_Ha' in r_), 'Ha', 12,
    'derived/d11_burst2.json:fcheck.relative[*].max_abs_E_vs_burst1_log_Ha (max over the four DFT methods)',
    'largest |ΔE| between the burst-2 and burst-1 runs of the four DFT methods on the same geometries (the burst-1 log prints 10 decimals)')
add('fcheck.n_methods', FC['n_methods'], None, 0, 'derived/d11_burst2.json:fcheck.n_methods', 'methods in the functional check (four functionals and DF-MP2)')
add('fcheck.n_neutral_lower', FC['n_neutral_lower'], None, 0, 'derived/d11_burst2.json:fcheck.n_neutral_lower', 'methods that put the neutral end point below the salt end point')
add('fcheck.dE_O_N_range', FC['dE_O_N_range'], 'kcal mol-1', 2, 'derived/d11_burst2.json:fcheck.dE_O_N_range', 'range of E(neutral) − E(salt) over the five methods (M06-2X lowest, B3LYP-D3(BJ) highest)')
add('fcheck.top_N_range', FC['top_N_range'], 'kcal mol-1', 2, 'derived/d11_burst2.json:fcheck.top_N_range', 'range of E(1.25 Å) − E(salt) over the five methods (M06-2X lowest, ωB97X-D highest)')
add('fcheck.b3lyp_consistency', max(abs(FSPJ['energies']['b3lyp-d3bj'][GN[g]]['E_Ha'] - [p for p in D4['points'] if (p['chain'], p['index']) == c][0]['E_Ha'])
                                    for g, c in [('end_N', ('A', 0)), ('top1', ('A', 2)), ('end_O', ('B', 0))]) * HA_KCAL, 'kcal mol-1', 6,
    f'{FSPS} vs scan/A_00_anchor.json, A_02_d1.250.json, B_00_anchor.json:E_Ha', 'largest difference between the B3LYP-D3(BJ) single points and the scan energies (identity check)')
for g, gl in GN.items():
    add(f'fcheck.geom.{g}.sha256_ok', FC['geometries'][gl]['sha256_ok'], None, None, f'derived/d11_burst2.json:fcheck.geometries["{gl}"]',
        f'geometry file of the {gl} point hashes to the sha256 recorded in fsp_scan.json')
CAVEATS.append('Functional check (five methods, fsp_scan.json): every method puts the neutral O5-protonated end point below the salt end point '
               '(B3LYP-D3 −0.54, ωB97X-D −0.72, PBE0-D3 −1.27, DF-MP2 −1.42, M06-2X −2.47 kcal/mol). The 1.25 Å point lies +0.72 (B3LYP-D3), '
               '+0.98 (ωB97X-D) and +0.20 (DF-MP2) kcal/mol above the salt but −0.05 (PBE0-D3) and −0.51 (M06-2X) below it, so the barrier '
               'is small or absent depending on the method, and always below the N–H zero-point energy (3.13 kcal/mol). These are single points '
               'at B3LYP-D3 geometries, not re-optimised, so the "barrier" of a method is its energy at one fixed geometry, not its saddle point. '
               'DF-MP2 includes the continuum only through the RHF reaction field (no correlation response of the solvent).')

# =====================================================================================================================
section('8. Zero-point energy of the N–H stretch against the barrier',
        'From zpe_barrier.json (1-D model not re-derived); the N–H stretch frequency re-derived by our own Hessian '
        'diagonalisation (derive/d05_freq.py) and ½hν recomputed.')
zb = J(RES / 'sfa_pair_pcm_d3/scan/zpe_barrier.json')
m65 = [m for m in D5['modes'] if abs(m['nu'] - 2190.2) < 1][0]
add('zpe.nu_NH', zb['nh_stretch']['normal_mode']['freq_cm-1'], 'cm-1', 0, 'scan/zpe_barrier.json:nh_stretch.normal_mode.freq_cm-1',
    'harmonic N1–H1A stretch (normal mode 65 of 75; 97 % N1–H1A stretch character)', case='sfa_pair_pcm_d3', check=m65['nu'], tol=0.5)
add('zpe.NH', zb['nh_stretch']['normal_mode']['zpe_half_h_nu_kcal'], 'kcal mol-1', 2, 'scan/zpe_barrier.json:nh_stretch.normal_mode.zpe_half_h_nu_kcal',
    'harmonic zero-point energy ½hν of the N1–H1A stretch', case='sfa_pair_pcm_d3', check=D4['zpe']['zpe_half_hnu_kcal'], tol=1e-4)
add('zpe.NH_kJ', D4['zpe']['zpe_half_hnu_kJ'], 'kJ mol-1', 1, 'derived/d04_scan.json:zpe.zpe_half_hnu_kJ', 'the same in kJ mol⁻¹', case='sfa_pair_pcm_d3')
add('zpe.ratio', D4['zpe']['ratio_zpe_over_barrier'], None, 1, 'derived/d04_scan.json:zpe.ratio_zpe_over_barrier', 'ZPE(N–H) / electronic barrier (0.72)', case='sfa_pair_pcm_d3')
add('zpe.nh_character', zb['nh_stretch']['normal_mode']['nh_stretch_character'], None, 2, 'scan/zpe_barrier.json:nh_stretch.normal_mode.nh_stretch_character', 'N1–H1A stretch character of the 2190 cm⁻¹ mode')
add('zpe.local_mode', zb['nh_stretch']['local_mode']['freq_cm-1'], 'cm-1', 0, 'scan/zpe_barrier.json:nh_stretch.local_mode', 'Konkoli–Cremer adiabatic local-mode N1–H1A frequency', case='sfa_pair_pcm_d3')
add('zpe.local_mode_zpe', zb['nh_stretch']['local_mode']['zpe_half_h_nu_kcal'], 'kcal mol-1', 2, 'scan/zpe_barrier.json:nh_stretch.local_mode.zpe_half_h_nu_kcal', 'ZPE of the local mode (still 3.7× the barrier)', case='sfa_pair_pcm_d3')
add('zpe.scan_curvature', zb['scan_curvature']['freq_cm-1'], 'cm-1', 0, 'scan/zpe_barrier.json:scan_curvature.freq_cm-1', 'frequency from the curvature of the relaxed scan at the salt minimum', case='sfa_pair_pcm_d3')
add('zpe.oned_E0', zb['one_d_model']['E0_above_reference_kcal'], 'kcal mol-1', 2, 'scan/zpe_barrier.json:one_d_model.E0_above_reference_kcal',
    '1-D DVR ground level on the relaxed profile (m = 0.940 u), above the salt minimum', case='sfa_pair_pcm_d3', status='not independently re-derived')
add('zpe.oned_E0_minus_top', zb['one_d_model']['E0_minus_barrier_top_kcal'], 'kcal mol-1', 2, 'scan/zpe_barrier.json:one_d_model.E0_minus_barrier_top_kcal',
    '1-D ground level minus the barrier top (≈ at the top)', case='sfa_pair_pcm_d3', status='not independently re-derived')
add('zpe.oned_P_O', zb['one_d_model']['levels'][0]['P_O_side'], None, 2, 'scan/zpe_barrier.json:one_d_model.levels[0].P_O_side',
    '1-D ground-state probability on the O side (d > barrier)', case='sfa_pair_pcm_d3', status='not independently re-derived')
add('zpe.oned_mass', zb['effective_mass_amu'], 'u', 3, 'scan/zpe_barrier.json:effective_mass_amu', 'effective mass of the 1-D model (N–H reduced mass)')
CAVEATS.append('ZPE vs barrier: harmonic ½hν of one normal mode (3.13 kcal/mol) against an electronic, B3LYP-D3, 1-D relaxed-scan barrier (0.72); the O-side minimum has had its own Hessian since VM burst 3 (section 33 gives the harmonic Δ(E + ZPE) of the two minima); the 1-D comparison of ½hν with the barrier is unchanged.')

# =====================================================================================================================
section('9. Vibrational frequencies, IR intensities and thermochemistry (IEF-PCM ion pair)',
        'Finite-difference Hessian, 75 real modes, 0 imaginary. Frequencies re-derived by our own diagonalisation of '
        f"the archived hessian.npy (max |Δν| {D5['max_abs_nu_diff_vs_freq_json']:.2f} cm⁻¹ vs freq.json, average atomic masses), "
        f"IR intensities from dipole_derivatives.npy (max |ΔI| {D5['max_abs_I_diff_vs_freq_json']:.2f} km mol⁻¹). Assignments from "
        'bond-stretch/angle-bend projections of each mode (heuristic character, not a PED). Scaled = × 0.9688, the '
        'fundamental-frequency factor of Merrick, Moran & Radom 2007 for B3-LYP/6-311+G(d,p); 6-311++G(d,p) is not in their '
        'table, and 6-311+G(d,p) is the closest basis they list (it lacks only the diffuse s shell on H). Verified 29 Sep 2026: '
        'Crossref record and abstract through lit.py (not retracted; the abstract names 6-311+G(d,p) and a B3-LYP basis-set '
        'series up to 6-311++G(3df,3pd)); the full text is closed (OpenAlex oa_status closed, publisher 403), so the table row '
        'was read from two independent tabulations that agree (see freq.scale_factor). Alternatives for this level are listed '
        'as freq.scale_factor.alt_*; they span 0.959–0.969, i.e. about 1 % (≈ 35 cm⁻¹ at 3500 cm⁻¹), which is the '
        'uncertainty of any scaled wavenumber.')
fr = J(RES / 'sfa_pair_pcm_d3/freq.json')
SF = 0.9688
add('freq.scale_factor', SF, None, 4, 'merrick2007evaluation (doi:10.1021/jp073974n), fundamentals, row B3-LYP/6-311+G(d,p); lit.py check + abstract 29 Sep 2026; '
    'table row read from github.com/liyuanhe211/Collection_of_Frequency_Scale_Factors (row B3LYP/6-311+G(d,p): ZPE 0.9887, fundamentals 0.9688, source doi 10.1021/jp073974n) '
    'and from Table 1 of doi:10.1016/j.newar.2017.01.001 (arXiv:1702.03438: 6-311+G** 0.9688; its Merrick-footnoted 6-311+G(3df,3pd) 0.9672 and 6-311++G(3df,3pd) 0.9673 match the compilation)',
    'fundamental-frequency scale factor (Merrick et al., B3-LYP/6-311+G(d,p), the closest basis they tabulate to 6-311++G(d,p))',
    note='confirmed in the full text (29 Sep 2026; author decision: keep 0.9688): Merrick, Moran & Radom 2007, Table 2 (p. 11684) and Table 6 (p. 11688), '
         'row B3-LYP/6-311+G(d,p), fundamentals 0.9688, rms_ov 35 cm⁻¹; Table 3 (p. 11686) 0.9689/0.9688/0.9688 with default/tight/very tight '
         'geometry criteria. They do not tabulate 6-311++G(d,p); in their own Table 6 adding the diffuse H shell to 6-311+G(3df,3pd) changes '
         'the factor by 0.0001 (0.9672 → 0.9673). Notes: paper/audit/research/irscale.md')
add('freq.scale_factor.merrick_3df3pd_plus', 0.9672, None, 4, 'Merrick 2007 row B3-LYP/6-311+G(3df,3pd), via the two tabulations named in freq.scale_factor', 'Merrick fundamental factor, B3-LYP/6-311+G(3df,3pd) (provenance only)')
add('freq.scale_factor.merrick_3df3pd_plusplus', 0.9673, None, 4, 'Merrick 2007 row B3-LYP/6-311++G(3df,3pd), via the two tabulations named in freq.scale_factor', 'Merrick fundamental factor, B3-LYP/6-311++G(3df,3pd): the diffuse H shell changes the factor by 0.0001 (provenance only)')
add('freq.scale_factor.alt_cccbdb', 0.967, None, 3, 'NIST CCCBDB Release 22 (May 2022), cccbdb.nist.gov/vibscalejustx.asp, B3LYP/6-311G** (read 29 Sep 2026; no 6-311+G** or 6-311++G** entry)', 'alternative fundamental factor, B3LYP/6-311G(d,p) (provenance only)')
add('freq.scale_factor.alt_truhlar', 0.959, None, 3, 'Truhlar group Database of Frequency Scale Factors v5 (22 Jul 2021), comp.chem.umn.edu/freqscale, row B3LYP/6-311++G(d,p): ZPE 0.985 (direct), F 0.959 (reduced model, F/ZPE = 0.974; ref. 27 unpublished)', 'alternative fundamental factor, B3LYP/6-311++G(d,p) (provenance only)')
add('freq.scale_factor.alt_kashinski', 0.963, None, 3, 'Kashinski et al. 2017 (doi:10.1021/acs.jpca.6b12147), B3LYP/6-311++G(d,p), as listed by the liyuanhe211 compilation (26 molecules, 99 frequencies; the value as listed there)', 'alternative fundamental factor, B3LYP/6-311++G(d,p) (provenance only)')
add('freq.n_modes', 75, None, 0, 'sfa_pair_pcm_d3/freq.json:frequencies_cm-1', '3N − 6 modes of the ion pair')
add('freq.n_imag', fr['n_imaginary'], None, 0, 'sfa_pair_pcm_d3/freq.json:n_imaginary', 'imaginary frequencies (true minimum)', check=D5['n_imag'], tol=0)
add('freq.lowest', fr['frequencies_cm-1'][0], 'cm-1', 0, 'sfa_pair_pcm_d3/freq.json:frequencies_cm-1[0]', 'lowest harmonic frequency (inter-ion libration)', check=D5['modes'][0]['nu'], tol=0.5)
add('freq.fd_asymmetry', fr['fd_hessian_asymmetry'], None, 4, 'sfa_pair_pcm_d3/freq.json:fd_hessian_asymmetry', 'relative asymmetry of the FD Hessian')
BANDS = [(65, 'NH_Hbond', 'N1–H1A stretch of the charge-assisted N1–H1A···O5 bond (coupled to H1A···O5); strongest IR band'),
         (70, 'NH3_1', 'NH3⁺ stretch, mainly N1–H1B'), (72, 'NH3_2', 'NH3⁺ stretch, mainly N1–H1C'),
         (71, 'NH2_s', 'NH2 symmetric stretch (N2)'), (73, 'NH2_as', 'NH2 antisymmetric stretch (N2)'),
         (74, 'OH_1', 'O–H stretches of H2PO4− (O3–H3A/O4–H4A, in phase)'), (75, 'OH_2', 'O–H stretches of H2PO4− (out of phase)'),
         (66, 'CH_1', 'aromatic C–H stretch (lowest)'), (69, 'CH_4', 'aromatic C–H stretch (highest)'),
         (62, 'ring', 'ring C=C stretch + C–H bend'), (60, 'NH2_sc', 'NH2 scissoring'),
         (59, 'NH3_def1', 'NH3⁺ deformation (involving H1A)'), (61, 'NH3_def2', 'NH3⁺ deformation (H1B–N1–H1C)'),
         (54, 'SO2_as', 'SO2 antisymmetric stretch coupled with NH2 wag'), (45, 'SO2_s', 'SO2 symmetric stretch + ring'),
         (53, 'CN', 'C1–N1 stretch + ring C–H bend'), (47, 'CS', 'C4–S1 stretch + ring'),
         (51, 'PO2_as', 'terminal P–O (O5, O6) stretch + P–O–H bends; second-strongest IR band'),
         (42, 'POH_b', 'P–O3–H3A in-plane bend'), (41, 'POH_b2', 'P–O4–H4A in-plane bend'),
         (34, 'P_OH_s1', 'P–(OH)2 stretch (P–O3 dominant)'), (33, 'P_OH_s2', 'P–(OH)2 stretch (P–O4 dominant)'),
         (31, 'SN', 'S1–N2 stretch + NH2 wag'), (26, 'OPO_b', 'O–P–O bend + H1A···O5 stretch'),
         (14, 'HB_str', 'inter-ion N1–H1A···O5 stretch/POH torsion region')]
for n, key, what in BANDS:
    m = D5['modes'][n - 1]
    fj = fr['ir'][n - 1]
    add(f'freq.{key}', fj['freq_cm-1'], 'cm-1', 0, f'sfa_pair_pcm_d3/freq.json:ir[{n - 1}].freq_cm-1',
        f'mode {n}: {what} (unscaled harmonic)', case='sfa_pair_pcm_d3', check=m['nu'], tol=0.5,
        note='character: ' + ', '.join(f'{a} {b:.2f}' for a, b in m['top_internal'][:3]))
    add(f'freq.{key}.scaled', fj['freq_cm-1'] * SF, 'cm-1', 0, f'freq.json:ir[{n - 1}] × freq.scale_factor', f'mode {n}: {what}, scaled × {SF}', case='sfa_pair_pcm_d3')
    add(f'freq.{key}.I', fj['intensity_km_mol'], 'km mol-1', 0, f'sfa_pair_pcm_d3/freq.json:ir[{n - 1}].intensity_km_mol',
        f'mode {n}: IR intensity', case='sfa_pair_pcm_d3', check=m['I_km'], tol=0.5)
CAVEATS.append('The 2190 cm⁻¹ N1–H1A stretch is a harmonic value for a very short, anharmonic N⁺–H···O⁻ bond (N1···O5 2.581 Å in IEF-PCM, 2.760 Å in the crystal); a scaled harmonic number is not a prediction of the observed band, which in the crystal would be broad and at higher wavenumber (longer N···O). The O–H stretches at ~3800 cm⁻¹ are of free O–H groups in the ion pair: in the crystal both O–H donate strong H-bonds (O4–H4A···O5 2.500 Å), so these bands do not describe the solid.')
th = D5['thermo']
add('thermo.ZPE', fr['ZPE_Ha'] * HA_KCAL, 'kcal mol-1', 2, 'sfa_pair_pcm_d3/freq.json:ZPE_Ha', 'total harmonic ZPE of the ion pair (unscaled)', case='sfa_pair_pcm_d3', check=th['ZPE_kcal'], tol=0.01)
add('thermo.ZPE_Ha', fr['ZPE_Ha'], 'Ha', 6, 'sfa_pair_pcm_d3/freq.json:ZPE_Ha', 'total ZPE in hartree', case='sfa_pair_pcm_d3', check=th['ZPE_Ha'], tol=1e-5)
add('thermo.H298', fr['H298_Ha'], 'Ha', 6, 'sfa_pair_pcm_d3/freq.json:H298_Ha', 'enthalpy at 298.15 K (RRHO ideal gas + PCM electronic energy)', case='sfa_pair_pcm_d3', check=th['H298_Ha'], tol=1e-5)
add('thermo.G298', fr['G298_Ha'], 'Ha', 6, 'sfa_pair_pcm_d3/freq.json:G298_Ha', 'Gibbs energy at 298.15 K, 1 atm', case='sfa_pair_pcm_d3', check=th['G298_Ha'], tol=1e-5)
add('thermo.E', fr['E_Ha'], 'Ha', 6, 'sfa_pair_pcm_d3/freq.json:E_Ha', 'electronic energy at the Hessian geometry', case='sfa_pair_pcm_d3')
add('thermo.S', th['S_J_molK'], 'J mol-1 K-1', 1, 'derived/d05_freq.json:thermo.S_J_molK', 'total entropy (own RRHO; trans 178.6, rot 141.5, vib 267.8)', case='sfa_pair_pcm_d3', status='own derivation')
add('thermo.Hcorr', th['H_corr_Ha'] * HA_KCAL, 'kcal mol-1', 2, 'derived/d05_freq.json:thermo.H_corr_Ha', 'thermal correction to the enthalpy (incl. ZPE)', case='sfa_pair_pcm_d3', status='own derivation')
add('thermo.Gcorr', th['G_corr_Ha'] * HA_KCAL, 'kcal mol-1', 2, 'derived/d05_freq.json:thermo.G_corr_Ha', 'thermal correction to the Gibbs energy', case='sfa_pair_pcm_d3', status='own derivation')
MISSING.append('Frequencies: only the IEF-PCM ion pair, the gas-phase sulfanilamide, cation (gas/PCM), H2PO4− (gas/PCM), H3PO4 (gas) and urea have Hessians. No frequencies for the gas-phase pair, the plain-B3LYP and def2-TZVP pairs, the O-side scan minimum, H3PO4 (PCM), and — relevant for docking — sulfanilamide (PCM) and the sulfonamidate (PCM): their minima are unconfirmed.')

# =====================================================================================================================
section('10. Frontier orbitals, DOS and reactivity descriptors',
        'Orbital energies re-read from the archived final.molden files (derive/d06_electronic.py), identical to fmo.json. '
        'ΔSCF: I = E(N−1) − E(N), A = E(N) − E(N+1), vertical, same functional/basis/grid/PCM (equilibrium solvation of '
        'the ions); η = (I − A)/2, χ = (I + A)/2, ω = χ²/2η. Koopmans-type values use −ε(HOMO), −ε(LUMO).')
FS = [('pcm', 'sfa_pair_pcm_d3', 'ion pair, IEF-PCM'), ('gas', 'sfa_pair_gas_d3', 'gas-phase pair (proton on O5: sulfanilamide·H3PO4)'),
      ('cation.pcm', 'cation_pcm_d3', 'cation, IEF-PCM'), ('cation.gas', 'cation_gas_d3', 'cation, gas'),
      ('sfa.pcm', 'sulfanilamide_pcm_d3', 'sulfanilamide, IEF-PCM'), ('sfa.gas', 'sulfanilamide_gas_d3', 'sulfanilamide, gas'),
      ('sfaan.pcm', 'sulfonamidate_pcm_d3', 'sulfonamidate, IEF-PCM'), ('h2po4.pcm', 'h2po4_pcm_d3', 'H2PO4−, IEF-PCM'),
      ('urea.gas', 'urea_gas_d3', 'urea, gas')]
for pre, c, lab in FS:
    f = J(RES / c / 'fmo.json')
    dd = D6['fmo'][c]
    add(f'fmo.{pre}.homo', f['HOMO_eV'], 'eV', 3, f'molecules/sfa/results/{c}/fmo.json:HOMO_eV', f'{lab}: E(HOMO) (B3LYP-D3)', case=c, check=dd['HOMO_eV'], tol=1e-4)
    add(f'fmo.{pre}.lumo', f['LUMO_eV'], 'eV', 3, f'molecules/sfa/results/{c}/fmo.json:LUMO_eV', f'{lab}: E(LUMO)', case=c, check=dd['LUMO_eV'], tol=1e-4)
    add(f'fmo.{pre}.gap', f['gap_eV'], 'eV', 3, f'molecules/sfa/results/{c}/fmo.json:gap_eV', f'{lab}: Kohn–Sham gap', case=c, check=dd['gap_eV'], tol=1e-4)
    if c in ('sfa_pair_pcm_d3', 'sulfanilamide_pcm_d3', 'sulfanilamide_gas_d3', 'urea_gas_d3', 'sfa_pair_gas_d3'):
        add(f'fmo.{pre}.dipole', f['dipole_Debye'], 'D', 2, f'molecules/sfa/results/{c}/fmo.json:dipole_Debye', f'{lab}: dipole moment', case=c)
fc = J(RES / 'sfa_pair_pcm_d3/fmo.json')['fragment_composition']
for orb in ['HOMO-1', 'HOMO', 'LUMO', 'LUMO+1']:
    add(f'fmo.pcm.comp.{orb}.cation', fc[orb]['4-sulfamoylanilinium cation'], '%', 1, f'sfa_pair_pcm_d3/fmo.json:fragment_composition.{orb}', f'{orb} on the cation (Mulliken fragment partition)', case='sfa_pair_pcm_d3')
dos = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/dos.json')
for orb in ['HOMO-1', 'HOMO', 'LUMO']:
    add(f'dos.pcm.{orb}.cation', dos['frontier'][orb]['4-sulfamoylanilinium cation'], '%', 1, f'sfa_pair_pcm_d3/mwfn_extra/dos.json:frontier.{orb}',
        f'{orb} on the cation (Hirshfeld orbital composition, Multiwfn)', case='sfa_pair_pcm_d3')
add('fmo.pcm.homo_m1', fc and J(RES / 'sfa_pair_pcm_d3/fmo.json')['HOMO-1_eV'], 'eV', 3, 'sfa_pair_pcm_d3/fmo.json:HOMO-1_eV', 'E(HOMO−1), 80 % on H2PO4−', case='sfa_pair_pcm_d3')
add('fmo.pcm.lumo_p1', J(RES / 'sfa_pair_pcm_d3/fmo.json')['LUMO+1_eV'], 'eV', 3, 'sfa_pair_pcm_d3/fmo.json:LUMO+1_eV', 'E(LUMO+1), cation', case='sfa_pair_pcm_d3')
add('dos.fwhm', 0.3, 'eV', 1, 'sfa_pair_pcm_d3/mwfn_extra/dos.json:method', 'Gaussian broadening of the TDOS/PDOS')
gfc = J(RES / 'sfa_pair_gas_d3/fmo.json')['fragment_composition']
add('fmo.gas.comp.HOMO.cation', gfc['HOMO']['4-sulfamoylanilinium cation'], '%', 1, 'sfa_pair_gas_d3/fmo.json:fragment_composition.HOMO',
    'gas-phase pair: HOMO on the sulfanilamide part (proton transferred; labels are the case fragments)', case='sfa_pair_gas_d3')
for pre, c, lab in [('pcm', 'sfa_pair_pcm_d3', 'ion pair, IEF-PCM'), ('cation.pcm', 'cation_pcm_d3', 'cation, IEF-PCM'),
                    ('cation.gas', 'cation_gas_d3', 'cation, gas'), ('sfa.pcm', 'sulfanilamide_pcm_d3', 'sulfanilamide, IEF-PCM'),
                    ('sfa.gas', 'sulfanilamide_gas_d3', 'sulfanilamide, gas')]:
    dj = J(RES / c / 'dscf.json')
    d6 = D6['dscf'][c]
    for q, jq, p, what in [('I', 'I_eV', 3, 'vertical ionisation energy'), ('A', 'A_eV', 3, 'vertical electron affinity'),
                           ('eta', 'hardness_eta_eV', 3, 'hardness η'), ('chi', 'electronegativity_chi_eV', 3, 'electronegativity χ'),
                           ('omega', 'electrophilicity_omega_eV', 3, 'electrophilicity ω'), ('S', 'softness_S_eV-1', 3, 'softness S = 1/2η')]:
        add(f'dscf.{pre}.{q}', dj['delta_scf'][jq], 'eV' if q != 'S' else 'eV-1', p, f'molecules/sfa/results/{c}/dscf.json:delta_scf.{jq}',
            f'{lab}: ΔSCF {what}', case=c, check=d6[q], tol=1e-4)
        add(f'dscf.{pre}.koop.{q}', dj['koopmans'][jq], 'eV' if q != 'S' else 'eV-1', p, f'molecules/sfa/results/{c}/dscf.json:koopmans.{jq}',
            f'{lab}: Koopmans-type {what} (from ε_HOMO/ε_LUMO)', case=c, check=d6['koopmans'][q], tol=1e-4)
    add(f'dscf.{pre}.fund_gap', dj['fundamental_gap_I_minus_A_eV'], 'eV', 3, f'molecules/sfa/results/{c}/dscf.json:fundamental_gap_I_minus_A_eV', f'{lab}: I − A', case=c)
add('dscf.sfa.gas.unbound', True, None, None, 'sulfanilamide_gas_d3/dscf.json:flags.electron_unbound_in_gas_phase', 'gas-phase sulfanilamide: A < 0 (−0.099 eV), the extra electron is unbound; the ΔSCF anion is a basis-set artefact')
MISSING.append('No ΔSCF for the sulfonamidate, H2PO4− or the gas-phase pair; no condensed Fukui functions for any species (SI section S9 "Fukui functions" has no data).')

# =====================================================================================================================
section('11. TD-DFT (IEF-PCM water, non-equilibrium)',
        'Wavelengths recomputed as 1239.842/E. "CT number" = fraction of the excitation that moves charge between '
        'cation and anion (Plasser–Lischka, Mulliken partition of X+Y); D_CT = hole–electron centroid distance.')
for c, fn, pre, lab in [('sfa_pair_pcm_d3', 'tddft_camb3lyp_noneq', 'td.pair.cam', 'ion pair, CAM-B3LYP'),
                        ('sfa_pair_pcm_d3', 'tddft_b3lyp_noneq', 'td.pair.b3lyp', 'ion pair, B3LYP'),
                        ('cation_pcm_d3', 'tddft_camb3lyp_noneq', 'td.cation.cam', 'cation, CAM-B3LYP'),
                        ('cation_pcm_d3', 'tddft_camb3lyp_eq', 'td.cation.cam_eq', 'cation, CAM-B3LYP equilibrium PCM'),
                        ('cation_pcm_d3', 'tddft_b3lyp_noneq', 'td.cation.b3lyp', 'cation, B3LYP'),
                        ('sulfanilamide_pcm_d3', 'tddft_camb3lyp_noneq', 'td.sfa.cam', 'sulfanilamide, CAM-B3LYP'),
                        ('sulfanilamide_pcm_d3', 'tddft_camb3lyp_eq', 'td.sfa.cam_eq', 'sulfanilamide, CAM-B3LYP equilibrium PCM'),
                        ('sulfanilamide_pcm_d3', 'tddft_b3lyp_noneq', 'td.sfa.b3lyp', 'sulfanilamide, B3LYP')]:
    t = D6['tddft'][f'{c}/{fn}']
    src = f'molecules/sfa/results/{c}/{fn}.json:states'
    for tag, st in [('S1', t['S1']), ('bright', t['first_f_ge_0p1']), ('max', t['brightest'])]:
        what = {'S1': 'S1', 'bright': f"first state with f ≥ 0.1 (S{st['state']})", 'max': f"strongest state (S{st['state']})"}[tag]
        add(f'{pre}.{tag}.E', st['E_eV'], 'eV', 3, f"{src}[{st['state'] - 1}].E_eV", f'{lab}: {what} excitation energy', case=c)
        add(f'{pre}.{tag}.nm', st['nm'], 'nm', 0, f"{src}[{st['state'] - 1}].lambda_nm", f'{lab}: {what} wavelength', case=c)
        add(f'{pre}.{tag}.f', st['f'], None, 3, f"{src}[{st['state'] - 1}].f", f'{lab}: {what} oscillator strength', case=c)
        add(f'{pre}.{tag}.assign', ', '.join(f'{a} ({100 * w:.0f} %)' for a, w in st['dominant']), None, None, f"{src}[{st['state'] - 1}].dominant", f'{lab}: {what} main configurations', case=c)
        if st['CT'] is not None:
            add(f'{pre}.{tag}.CT', st['CT'], None, 2, f"{src}[{st['state'] - 1}].ct.CT_number", f'{lab}: {what} inter-ion CT number', case=c)
            add(f'{pre}.{tag}.DCT', st['D_CT'], 'A', 2, f"{src}[{st['state'] - 1}].ct.D_CT_A", f'{lab}: {what} hole–electron distance', case=c)
    add(f'{pre}.peak', t['spectrum_peak_nm'], 'nm', 0, f'molecules/sfa/results/{c}/{fn}.json:spectrum_peak_nm', f'{lab}: maximum of the broadened spectrum (0.3 eV FWHM)', case=c)
    if t['max_CT_first5'] is not None:
        add(f'{pre}.maxCT5', t['max_CT_first5'], None, 2, f'derived/d06_electronic.json:tddft', f'{lab}: largest CT number among S1–S5', case=c)
add('td.all_converged', True, None, None, 'tddft_*.json:all_converged', 'all 20 roots converged in every TD-DFT run')
bm = [json.loads(l) for l in (REPO / 'data/benchmark/results.jsonl').read_text().splitlines() if l.strip()]
bt = [b for b in bm if b['test'] == 'tddft10_gas_24t'][0]
add('td.benchmark_gas', bt['per_copy'][0]['exc_ev_first3'][0], 'eV', 2, 'data/benchmark/results.jsonl: test tddft10_gas_24t, per_copy[0].exc_ev_first3[0]',
    'first excitation of the gas-phase ion pair at the X-ray geometry in the VM benchmark (B3LYP/6-311++G(d,p), no D3, 10 states; TDA per PLAN.md): a CT artefact for a colourless crystal, the reason CAM-B3LYP in water was chosen',
    status='benchmark, not a production result', note=f"first three: {bt['per_copy'][0]['exc_ev_first3']} eV")
CAVEATS.append('TD-DFT: with B3LYP four of the first five ion-pair states are anion→cation CT (CT number 0.61–0.91, D_CT 3.3–5.0 Å); with CAM-B3LYP none is (CT ≤ 0.04). The first strong CAM-B3LYP band of the pair (S2, 211 nm, f = 0.18) sits at the cation value (208 nm) — the pair spectrum is the cation spectrum. The released neutral sulfanilamide absorbs at 246 nm (S2, 246.45 nm, f = 0.50, non-equilibrium; 256 nm equilibrium). No experimental UV-Vis for SAPH is in hand.')

# =====================================================================================================================
section('12. Static NLO response, against urea at the same level',
        'Recomputed from the raw finite-field dipoles in nlo.json (derive/d06_electronic.py): α and β by central '
        'differences at F = 0.001 and 0.002 a.u. with one Richardson step; β_i = ⅓Σ_j(β_ijj + β_jij + β_jji), '
        'β_tot = |β|. Convention: Taylor (T) series μ = μ0 + αF + ½βF²; β(B convention) = β(T)/2. 1 a.u. of β = '
        '8.6392 × 10⁻³³ esu. Urea is gas phase; the like-for-like ratio uses the gas-phase pair at the IEF-PCM geometry.')
for pre, c, lab in [('pcm', 'sfa_pair_pcm_d3', 'ion pair, IEF-PCM (water)'), ('gasAtPcm', 'sfa_pair_gas_d3', 'ion pair at the IEF-PCM geometry, gas phase'),
                    ('sfa', 'sulfanilamide_pcm_d3', 'sulfanilamide, IEF-PCM'), ('urea', 'urea_gas_d3', 'urea, gas phase')]:
    n = D6['nlo'][c]
    src = n['file']
    add(f'nlo.{pre}.mu', n['json']['mu_D'], 'D', 2, f'{src}:mu_total_debye', f'{lab}: dipole moment', case=c, check=n['mu_D'], tol=1e-3)
    add(f'nlo.{pre}.alpha', n['json']['alpha_iso_au'], 'a.u.', 1, f'{src}:alpha_iso_au', f'{lab}: isotropic polarisability ᾱ', case=c, check=n['alpha_iso_au'], tol=0.01)
    add(f'nlo.{pre}.alpha_A3', n['alpha_iso_A3'], 'A3', 2, f'{src}:alpha_iso_A3', f'{lab}: ᾱ in Å³', case=c)
    add(f'nlo.{pre}.dalpha', n['alpha_aniso_au'], 'a.u.', 1, f'{src}:alpha_aniso_au', f'{lab}: polarisability anisotropy Δα', case=c)
    add(f'nlo.{pre}.beta', n['json']['beta_tot_au'], 'a.u.', 1, f'{src}:beta_tot_au', f'{lab}: β_tot (T convention)', case=c, check=n['beta_tot_au'], tol=0.05)
    add(f'nlo.{pre}.beta_esu', n['json']['beta_tot_esu'] * 1e30, '1e-30 esu', 3, f'{src}:beta_tot_esu', f'{lab}: β_tot in 10⁻³⁰ esu (T convention)', case=c)
    add(f'nlo.{pre}.beta_B', n['beta_tot_B_conv_au'], 'a.u.', 1, 'derived/d06_electronic.json:nlo', f'{lab}: β_tot in the B convention (= T/2)', case=c)
    add(f'nlo.{pre}.beta_par', n['beta_par_mu_au'], 'a.u.', 1, 'derived/d06_electronic.json:nlo.beta_par_mu_au', f'{lab}: β component along μ', case=c)
    add(f'nlo.{pre}.ratio_urea', n['ratio_to_urea'], None, 1, 'derived/d06_electronic.json:nlo.ratio_to_urea', f'{lab}: β_tot / β_tot(urea)', case=c)
add('nlo.richardson_effect', (D6['nlo']['sfa_pair_pcm_d3']['beta_F_only_au'] - D6['nlo']['sfa_pair_pcm_d3']['beta_tot_au']), 'a.u.', 1,
    'derived/d06_electronic.json:nlo.sfa_pair_pcm_d3.beta_F_only_au', 'change of β_tot of the PCM pair from the Richardson step (numerical stability)')
nj = J(RES / 'sfa_pair_pcm_d3/nlo.json')
add('nlo.pcm.alpha_cpks_check', nj['analytic_alpha_check']['rel_diff_iso'], None, 8, 'sfa_pair_pcm_d3/nlo.json:analytic_alpha_check.rel_diff_iso', 'relative difference of ᾱ between finite field and analytic CPKS')
CAVEATS.append('NLO: the IEF-PCM pair against gas-phase urea mixes environments (ratio 4.5); the like-for-like gas-phase ratio at the same geometry is 6.4. Static, electronic, single ion pair: no local-field factors, no crystal packing; Pc has no inversion centre, so SHG is symmetry-allowed, but a molecular β is not a crystal χ(2).')

# =====================================================================================================================
section('13. Interaction energies of the ion pair',
        'Counterpoise-corrected, recomputed from the stored fragment energies (derive/d06_electronic.py). E_int(CP) = '
        'E_AB − E_A(AB basis) − E_B(AB basis) + ΔD3; E_def = fragment in the complex minus the relaxed fragment; '
        'E_bind = E_int(CP) + ΣE_def. In PCM the counterpoise is approximate (each species in its own cavity).')
for key, fn, lab, case in [('int.pcm', 'sfa_pair_pcm_d3/interaction.json', 'IEF-PCM, cation + H2PO4− (ionic partition)', 'sfa_pair_pcm_d3'),
                           ('int.gas_neutral', 'sfa_pair_gas_d3/interaction.json', 'gas phase, sulfanilamide + H3PO4 (neutral partition, own geometry)', 'sfa_pair_gas_d3'),
                           ('int.gas_ionic', 'sfa_pair_gas_d3/interaction_geom-sfa_pair_pcm_d3.json', 'gas phase at the IEF-PCM geometry, cation + H2PO4− (ionic)', 'sfa_pair_gas_d3')]:
    i = D6['interaction'][fn]
    js = i['json']
    add(f'{key}.raw', js['E_int_raw_kcal'], 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.E_int_raw_kcal', f'{lab}: uncorrected E_int', case=case, check=i['E_int_raw_kcal'], tol=0.005)
    add(f'{key}.cp', js['E_int_CP_kcal'], 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.E_int_CP_kcal', f'{lab}: counterpoise-corrected E_int', case=case, check=i['E_int_CP_kcal'], tol=0.005)
    add(f'{key}.bsse', js['BSSE_kcal'], 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.BSSE_kcal', f'{lab}: BSSE', case=case, check=i['BSSE_kcal'], tol=0.005)
    add(f'{key}.d3', js['E_int_D3_kcal'], 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.E_int_D3_kcal', f'{lab}: D3(BJ) part of E_int', case=case, check=i['E_int_D3_kcal'], tol=0.005)
    for fr_, v in js['E_def_kcal'].items():
        short = {'4-sulfamoylanilinium cation': 'cation', 'dihydrogen phosphate anion': 'anion', 'sulfanilamide (neutral)': 'sulfanilamide',
                 'phosphoric acid (H2PO4- + proton on O5)': 'h3po4'}[fr_]
        add(f'{key}.def.{short}', v, 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.E_def_kcal["{fr_}"]', f'{lab}: deformation energy of the {fr_}', case=case, check=i['E_def_kcal'][fr_], tol=0.005)
    add(f'{key}.bind', js['E_bind_kcal'], 'kcal mol-1', 2, f'molecules/sfa/results/{fn}:results.E_bind_kcal', f'{lab}: E_bind = E_int(CP) + ΣE_def', case=case, check=i['E_bind_kcal'], tol=0.005)
add('int.pcm.bsse_gas_same_geom', J(RES / 'sfa_pair_pcm_d3/interaction.json')['gas_crosscheck_for_pcm']['BSSE_gas_same_geometry_kcal'], 'kcal mol-1', 2,
    'sfa_pair_pcm_d3/interaction.json:gas_crosscheck_for_pcm.BSSE_gas_same_geometry_kcal', 'BSSE of the same geometry in the gas phase (cross-check of the PCM counterpoise)')

# =====================================================================================================================
section('14. QTAIM bond critical points and IGMH (ion pair only)',
        'Multiwfn QTAIM on the B3LYP-D3 wavefunctions; E_HB = V(r_BCP)/2 (Espinosa 1998) recomputed here; binding '
        'energies from the Emamian et al. 2019 fits (neutral: −223.08ρ + 0.7423; charged: −332.34ρ − 1.0661 kcal/mol). '
        'An ion pair is outside both training sets: all E_HB values are indicative only.')
for c, pre in [('sfa_pair_pcm_d3', 'pcm'), ('sfa_pair_gas_d3', 'gas')]:
    q = D6['qtaim'][c]
    qj = J(RES / c / 'mwfn_extra/qtaim.json')
    add(f'qtaim.{pre}.n_cp', q['n_cp'], None, 0, f'molecules/sfa/results/{c}/mwfn_extra/qtaim.json:n_cp', f'{pre}: critical points found ({q["counts"]}); Poincaré–Hopf satisfied', case=c)
    for h, hj in zip(q['hbonds'], qj['hbonds']):
        nm = {'N1-H1A...O5': 'N1H1A_O5', 'C6-H6...O4': 'C6H6_O4', 'N1-H1B...O4#1_554': 'N1H1B_O4'}.get(h['name'], 'O5H1A_N1')
        base = f'qtaim.{pre}.{nm}'
        lab = h['name'].split(' (')[0].replace('#1_554', '')
        src = f"molecules/sfa/results/{c}/mwfn_extra/qtaim.json:hbonds[{h['name'].split(' (')[0]}]"
        add(base + '.rho', h['rho'], 'a.u.', 4, src + '.rho', f'{pre}: ρ at the {lab} BCP', case=c)
        add(base + '.lap', h['lap'], 'a.u.', 4, src + '.lap', f'{pre}: ∇²ρ at the {lab} BCP', case=c)
        add(base + '.H', h['Hd'], 'a.u.', 4, src + '.H', f'{pre}: total energy density H at the {lab} BCP (H < 0: partly covalent)', case=c)
        add(base + '.GV', h['mGV'], None, 2, src + '.minusG_over_V', f'{pre}: −G/V at the {lab} BCP', case=c)
        add(base + '.EHB', hj['E_HB_Espinosa_kcal'], 'kcal mol-1', 1, src + '.E_HB_Espinosa_kcal', f'{pre}: E_HB = V/2 of {lab} (indicative)', case=c,
            check=h['E_HB_kcal'], tol=0.01, status='indicative')
        add(base + '.BEn', h['BE_em_neutral'], 'kcal mol-1', 1, src + '.BE_Emamian_neutral_kcal', f'{pre}: Emamian neutral-fit BE of {lab} (indicative)', case=c,
            check=hj['BE_Emamian_neutral_kcal'], tol=0.01, status='indicative')
        add(base + '.BEc', h['BE_em_charged'], 'kcal mol-1', 1, src + '.BE_Emamian_charged_kcal', f'{pre}: Emamian charged-fit BE of {lab} (indicative)', case=c,
            check=hj['BE_Emamian_charged_kcal'], tol=0.01, status='indicative')
        add(base + '.DA', h['DA'], 'A', 3, src + '.D..A_A', f'{pre}: D···A of {lab}', case=c)
        add(base + '.nature', h['nature'], None, None, src + '.nature_Rozas2000', f'{pre}: Rozas–Alkorta–Elguero class of {lab}', case=c)
for c, pre in [('sfa_pair_pcm_d3', 'pcm'), ('sfa_pair_gas_d3', 'gas')]:
    ig = J(RES / c / 'mwfn_extra/igmh.json')
    add(f'igmh.{pre}.int_dg_inter', ig['int_dg_inter_au'], 'a.u.', 3, f'molecules/sfa/results/{c}/mwfn_extra/igmh.json:int_dg_inter_au', f'{pre}: integrated δg_inter between the two fragments', case=c)
    for t in ig['top_atom_pairs'][:4]:
        nm = t['pair'].replace('#1_554', '').replace('...', '_')
        add(f'igmh.{pre}.{nm}', t['pct'], '%', 1, f'molecules/sfa/results/{c}/mwfn_extra/igmh.json:top_atom_pairs[{t["pair"]}].pct',
            f'{pre}: share of δg_inter from the atom pair {t["pair"].replace("#1_554", "")} ({t["dist_A"]:.3f} Å)', case=c)
add('igmh.grid', 'Multiwfn grid option 3 (89 × 125 × 159 points)', None, None, 'sfa_pair_pcm_d3/mwfn_extra/igmh.json:grid_option, grid_points', 'IGMH grid')
nci = (RES / 'sfa_pair_pcm_d3/nci_spikes.csv').read_text().splitlines()
add('nci.pcm.hb_spike', float(nci[1].split(',')[0]), 'a.u.', 3, 'sfa_pair_pcm_d3/nci_spikes.csv row 1', 'sign(λ2)ρ of the strong H-bond spike in the RDG scatter (RDG cutoff raised to ρ < 0.1 a.u.)', case='sfa_pair_pcm_d3')
add('nci.rdg_iso', 0.5, 'a.u.', 1, 'sfa_pair_pcm_d3/RESULTS.md NCI caption', 'RDG isovalue of the NCI isosurface (colour scale −0.035 to +0.02 a.u.)')
CAVEATS.append('QTAIM/IGMH exist only for the 27-atom ion pair (two contacts: N1–H1A···O5 and C6–H6···O4; in the gas phase the transferred O5–H1A···N1, N1–H1B···O4 and C6–H6···O4). The strongest crystal contact between anions, O4–H4A···O5 (2.500 Å), and the other six crystal H-bonds have BCP data only from the five-pair gas-phase cluster of section 34 (M8, VM burst 3).')

# =====================================================================================================================
section('15. Atomic charges and the molecular electrostatic potential',
        'NPA from PySCF natural atomic orbitals (not NBO 7); Hirshfeld, ADCH and CM5 from Multiwfn. Fragment sums '
        'recomputed from the atomic charges. MEP on the 0.001 a.u. density isosurface (Multiwfn), IEF-PCM-polarised density.')
ch = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.json')
for sch in ['NPA', 'Hirshfeld', 'ADCH', 'CM5']:
    add(f'chg.pcm.cation_sum.{sch}', ch['fragment_sums']['4-sulfamoylanilinium cation'][sch], 'e', 3, f'sfa_pair_pcm_d3/mwfn_extra/charges.json:fragment_sums.cation.{sch}',
        f'IEF-PCM pair: {sch} charge on the cation (formal +1); 1 − q is the charge transferred to it', case='sfa_pair_pcm_d3',
        check=D6['charges']['sfa_pair_pcm_d3']['sums'][sch]['cation_sum'], tol=0.0005)
chg = J(RES / 'sfa_pair_gas_d3/mwfn_extra/charges.json')
for sch in ['NPA', 'Hirshfeld', 'ADCH', 'CM5']:
    k = [x for x in chg['fragment_sums'] if 'cation' in x][0]
    add(f'chg.gas.sfa_sum.{sch}', chg['fragment_sums'][k][sch], 'e', 3, f'sfa_pair_gas_d3/mwfn_extra/charges.json:fragment_sums["{k}"].{sch}',
        f'gas-phase pair: {sch} charge on the sulfanilamide part (cation minus H1A): a neutral pair', case='sfa_pair_gas_d3',
        check=D6['charges']['sfa_pair_gas_d3']['sums'][sch]['sulfanilamide_sum'], tol=0.0005)
labs = [l.split(',')[1] for l in (RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.csv').read_text().splitlines()[1:]]
for atom in ['N1', 'H1A', 'O5#1_554', 'O6#1_554', 'P1#1_554', 'S1', 'N2', 'O1']:
    i = labs.index(atom)
    for sch in ['NPA', 'CM5']:
        add(f'chg.pcm.{atom.replace("#1_554", "")}.{sch}', ch['atomic'][sch][i], 'e', 3, f'sfa_pair_pcm_d3/mwfn_extra/charges.json:atomic.{sch}[{i}]',
            f'IEF-PCM pair: {sch} charge on {atom.replace("#1_554", "")}', case='sfa_pair_pcm_d3')
CAVEATS.append('ADCH atomic charges in the IEF-PCM pair are unphysical on the ring C–H pairs (C6 +0.78, C3 +0.52 with H6 −0.22, H3 −0.20 e); use NPA/CM5 atom values and ADCH only as fragment sums.')
ms = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json')
w = ms['whole']
add('mep.pcm.min', w['global_min_au'] * HA_KCAL, 'kcal mol-1', 1, 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json:whole.global_min_au', 'IEF-PCM pair: global MEP minimum (near O6, 1.86 Å)', case='sfa_pair_pcm_d3')
add('mep.pcm.max', w['global_max_au'] * HA_KCAL, 'kcal mol-1', 1, 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json:whole.global_max_au', 'IEF-PCM pair: global MEP maximum (near H1C, 1.20 Å)', case='sfa_pair_pcm_d3')
add('mep.pcm.min_au', w['global_min_au'], 'a.u.', 4, 'mep_summary.json:whole.global_min_au', 'IEF-PCM pair: global minimum in a.u.', case='sfa_pair_pcm_d3')
add('mep.pcm.max_au', w['global_max_au'], 'a.u.', 4, 'mep_summary.json:whole.global_max_au', 'IEF-PCM pair: global maximum in a.u.', case='sfa_pair_pcm_d3')
for e in ms['extrema']:
    if e['type'] == 'minimum' and e['global'] != 'yes':
        near = e['nearest_atoms'].split(' (')[0].replace('#1_554', '')
        add(f'mep.pcm.min_{near}', e['ESP_kcal_mol'], 'kcal mol-1', 1, f'mep_summary.json:extrema (minimum, {e["nearest_atoms"]})', f'IEF-PCM pair: local MEP minimum near {near}', case='sfa_pair_pcm_d3')
    if e['type'] == 'maximum' and e['global'] != 'yes' and e['ESP_kcal_mol'] > 50:
        near = e['nearest_atoms'].split(' (')[0].replace('#1_554', '')
        add(f'mep.pcm.max_{near}', e['ESP_kcal_mol'], 'kcal mol-1', 1, f'mep_summary.json:extrema (maximum, {e["nearest_atoms"]})', f'IEF-PCM pair: local MEP maximum near {near}', case='sfa_pair_pcm_d3')
for q, u, p, what in [('area_A2', 'A2', 1, 'surface area'), ('positive_area_A2', 'A2', 1, 'positive area'), ('negative_area_A2', 'A2', 1, 'negative area'),
                      ('polar_area_pct', '%', 1, 'polar surface fraction (|V| > 10 kcal/mol)'), ('MPI_eV', 'eV', 2, 'molecular polarity index'),
                      ('sigma2_tot_kcal2', 'kcal2 mol-2', 0, 'σ²_tot'), ('nu', None, 3, 'balance ν'), ('Pi_kcal', 'kcal mol-1', 1, 'average deviation Π')]:
    add(f'mep.pcm.{q}', w[q], u, p, f'mep_summary.json:whole.{q}', f'IEF-PCM pair MEP surface: {what}', case='sfa_pair_pcm_d3')
sm = J(RES / 'sulfanilamide_pcm_d3/mwfn_extra/mep_summary.json')['whole']
add('mep.sfa.min', sm['global_min_au'] * HA_KCAL, 'kcal mol-1', 1, 'sulfanilamide_pcm_d3/mwfn_extra/mep_summary.json:whole.global_min_au', 'sulfanilamide (IEF-PCM): global MEP minimum', case='sulfanilamide_pcm_d3')
add('mep.sfa.max', sm['global_max_au'] * HA_KCAL, 'kcal mol-1', 1, 'sulfanilamide_pcm_d3/mwfn_extra/mep_summary.json:whole.global_max_au', 'sulfanilamide (IEF-PCM): global MEP maximum', case='sulfanilamide_pcm_d3')

# =====================================================================================================================
section('16. Hirshfeld surfaces of the crystal (X-ray structure, neutron-normalised X–H)',
        'Own Python implementation (spherical ROHF/aug-cc-pVQZ free atoms), cross-checked against Multiwfn; contact '
        'percentages re-aggregated from the per-vertex data (derive/d07_hirshfeld.py). X–H set to 1.083 (C), 1.009 (N), '
        '0.983 (O) Å; cluster cut-off 6.0 Å.')
hs = J(RES / 'hirshfeld_xray/hirshfeld_summary.json')
for ion in ['cation', 'anion']:
    t = hs['targets'][ion]
    add(f'hs.{ion}.area', t['area_A2'], 'A2', 1, f'hirshfeld_xray/hirshfeld_summary.json:targets.{ion}.area_A2', f'{ion}: Hirshfeld surface area (Multiwfn {D7["multiwfn_compare"][ion]["area_A2"][1]:.1f})', check=D7[ion]['area_A2'], tol=0.05)
    add(f'hs.{ion}.volume', t['volume_A3'], 'A3', 1, f'hirshfeld_summary.json:targets.{ion}.volume_A3', f'{ion}: enclosed volume')
    add(f'hs.{ion}.globularity', t['globularity'], None, 3, f'hirshfeld_summary.json:targets.{ion}.globularity', f'{ion}: globularity')
    add(f'hs.{ion}.dnorm_min', t['d_norm_min'], None, 3, f'hirshfeld_summary.json:targets.{ion}.d_norm_min', f'{ion}: smallest d_norm ({D7[ion]["d_norm_min_pair"]})', check=D7[ion]['d_norm_min'], tol=0.001)
    add(f'hs.{ion}.n_neighbours', t['n_neighbours_in_contact'], None, 0, f'hirshfeld_summary.json:targets.{ion}.n_neighbours_in_contact', f'{ion}: molecules in contact')
    for cc in t['contacts_reciprocal']:
        key = '...'.join(sorted(cc['contact'].split('/')[0].split('...')))
        if cc['pct'] < 0.1:
            continue
        add(f'hs.{ion}.{key.replace("...", "")}', cc['pct'], '%', 1, f'hirshfeld_summary.json:targets.{ion}.contacts_reciprocal[{cc["contact"]}]',
            f'{ion}: share of the surface in {cc["contact"]} contacts', check=D7[ion]['pct'][key], tol=0.01)
    for rs in t['red_spots'][:4]:
        nm = f"{rs['inside'].split('#')[0]}_{rs['outside'].split('#')[0]}"
        add(f'hs.{ion}.spot.{nm}', rs['dnorm'], None, 3, f'hirshfeld_xray/{ion}_red_spots.csv', f"{ion}: d_norm red spot {rs['inside']}···{rs['outside']} ({rs['dist_A']:.3f} Å)")
CAVEATS.append('Hirshfeld: our surface agrees with Multiwfn on area (0.01 %/1.2 %) and contact percentages (≤ 0.1 point), but the d_norm minima differ (cation −0.721 vs −0.383; anion −0.842 vs −0.446) because the Multiwfn d_norm does not follow the CrystalExplorer formula with Bondi radii (same cluster, same normalised X–H; d_i and d_e agree to ~0.01 Å; dft/hirshfeld.py docstring). Quote our d_norm (CrystalExplorer definition) and say so. CrystalExplorer energy frameworks were not computed.')

# =====================================================================================================================
section('17. Periodic PBE-D3(BJ) model of the crystal',
        'Quantum ESPRESSO, cell fixed at the X-ray values, 54 atoms. Energies from periodic_qe.json (Ry) converted '
        'here; structures re-analysed (derive/d08_periodic.py). Heavy-atom RMSD with the mean shift removed (Pc has a '
        'floating origin in the a–c plane).')
e8 = D8['energies']
pj = J(RES / 'periodic_qe/periodic_qe.json')
add('periodic.E_i', e8['E_Ry']['i'], 'Ry', 7, 'periodic_qe/periodic_qe.json:relax.i.energy_Ry', 'relaxed energy from the X-ray start (i)', case='periodic_qe')
add('periodic.E_ii', e8['E_Ry']['ii'], 'Ry', 7, 'periodic_qe/periodic_qe.json:relax.ii.energy_Ry', 'relaxed energy from the O5-protonated start (ii)', case='periodic_qe')
add('periodic.E_iifix', e8['E_Ry']['iifix'], 'Ry', 7, 'periodic_qe/periodic_qe.json:relax.iifix.energy_Ry', 'energy with H1A fixed at 0.98 Å from the X-ray O5 (iifix)', case='periodic_qe')
add('periodic.dE_ii_i_meV', e8['dE_ii_minus_i_meV_cell'], 'meV', 2, 'periodic_qe.json:comparison.E_ii_minus_i_meV_cell', 'E(ii) − E(i) per cell: the same minimum', case='periodic_qe',
    check=pj['comparison']['E_ii_minus_i_meV_cell'], tol=0.01)
add('periodic.dE_ii_i_kJ', e8['dE_ii_minus_i_kJ_fu'], 'kJ mol-1', 2, 'derived/d08_periodic.json:energies.dE_ii_minus_i_kJ_fu', 'E(ii) − E(i) per formula unit', case='periodic_qe')
add('periodic.rmsd_ii_i', pj['comparison']['rmsd_all_atoms_final_ii_vs_final_i_A'], 'A', 2, 'periodic_qe.json:comparison.rmsd_all_atoms_final_ii_vs_final_i_A', 'RMSD between the relaxed (ii) and (i) cells', case='periodic_qe')
add('periodic.verdict', 'the O5-protonated start (ii) relaxes back to the salt: single minimum, proton on N1', None, None, 'periodic_qe.json:comparison.verdict', 'outcome of the two starts', case='periodic_qe')
add('periodic.dE_iifix_i_meV', e8['dE_iifix_minus_i_meV_cell'], 'meV', 1, 'periodic_qe.json:comparison.constrained_neutral.E_iifix_minus_i_meV_cell',
    'E(iifix) − E(i) per cell', case='periodic_qe', check=pj['comparison']['constrained_neutral']['E_iifix_minus_i_meV_cell'], tol=0.05, status='NOT a neutral-form energy (see note)',
    note='in iifix all non-fixed atoms slid by 0.66 Å (mean heavy-atom shift 0.10, 0, −0.65 Å; floating origin of Pc) so that N1 moved to the fixed H1A: the relaxed structure has the proton on N1 (N1–H1A 1.068, H1A···O5 1.660 Å, 165.7°); max force 0.013 Ry/Bohr (fixed atoms)')
add('periodic.dE_iifix_i_kJ', e8['dE_iifix_minus_i_kJ_fu'], 'kJ mol-1', 2, 'periodic_qe.json:comparison.constrained_neutral.E_iifix_minus_i_kJ_mol_per_formula_unit',
    'E(iifix) − E(i) per formula unit', case='periodic_qe', check=pj['comparison']['constrained_neutral']['E_iifix_minus_i_kJ_mol_per_formula_unit'], tol=0.005,
    status='NOT a neutral-form energy (see periodic.dE_iifix_i_meV note)')
add('periodic.iifix_shift_c', D8['structures']['relax_iifix']['mean_heavy_shift_vs_start_ii'][2], 'A', 2, 'derived/d08_periodic.json:structures.relax_iifix.mean_heavy_shift_vs_start_ii',
    'mean z-shift of all heavy atoms in iifix relative to its start (whole-crystal translation towards the fixed H1A)', case='periodic_qe')
add('periodic.iifix_N1_moved', D8['structures']['relax_iifix']['N1_moved_A'], 'A', 2, 'derived/d08_periodic.json:structures.relax_iifix.N1_moved_A', 'distance N1 moved in iifix (H1A fixed)', case='periodic_qe')
DISAGREE.append('SETTLED (burst 2, 28 Sep 2026): the old KEY STORY "constrained neutral +3.9 kJ/mol per formula unit" (relax_iifix, H1A alone held in space) is WITHDRAWN. Its optimiser translated every other atom by ~0.66 Å (floating origin of Pc) until N1 re-bonded to H1A (final N1–H1A 1.068 Å, H1A···O5 1.660 Å, δ = −0.592 Å, proton on N1), so the +3.87 kJ/mol (80.2 meV/cell) compares two salt structures. The iifix facts are kept for provenance with value null and status WITHDRAWN (their numbers sit in withdrawn_value, where numbers_audit cannot match them): never cite them. The constrained neutral that holds is relax_iirfix (H1A and O5 both held, O5–H1A 0.98 Å, no slide): +57.1 kJ/mol per formula unit at PBE-D3(BJ), +52.8 with the PBE0-D3(BJ) correction (section 28). iifix may be mentioned only as a methodological note (why O5 must be held too).')
ri = D8['structures']['relax_i']
add('periodic.rmsd_heavy', pj['relax']['i']['rmsd_heavy_vs_xray']['rmsd_A'], 'A', 3, 'periodic_qe.json:relax.i.rmsd_heavy_vs_xray.rmsd_A',
    'heavy-atom RMSD of the relaxed cell (i) against the X-ray, mean shift removed (32 atoms)', case='periodic_qe', check=ri['rmsd_heavy_shift_removed'], tol=0.0006)
add('periodic.rmsd_heavy_raw', ri['rmsd_heavy_raw'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.rmsd_heavy_raw', 'same without removing the mean shift', case='periodic_qe')
add('periodic.maxdev_heavy', ri['max_dev_heavy'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.max_dev_heavy', f"largest heavy-atom displacement ({ri['max_dev_atom']})", case='periodic_qe')
add('periodic.rmsd_all', ri['rmsd_all_shift_removed'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.rmsd_all_shift_removed', 'all-atom RMSD incl. H (X-ray H short)', case='periodic_qe')
add('periodic.pressure', pj['relax']['i']['pressure_kbar_final'], 'kbar', 1, 'periodic_qe.json:relax.i.pressure_kbar_final', 'residual pressure of the relaxed cell at the fixed X-ray lattice', case='periodic_qe')
add('periodic.max_force_i', pj['relax']['i']['max_force_Ry_bohr'], 'Ry/bohr', 5, 'periodic_qe.json:relax.i.max_force_Ry_bohr', 'largest residual force in (i)', case='periodic_qe')
for h in ri['hbonds']:
    nm = h['name'].replace('-', '').replace('...', '_')
    add(f'periodic.hb.{nm}.DA', h['DA'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.hbonds', f"periodic: D···A of {h['name']} ({h['sym']})", case='periodic_qe')
    add(f'periodic.hb.{nm}.HA', h['HA'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.hbonds', f"periodic: H···A of {h['name']}", case='periodic_qe')
    add(f'periodic.hb.{nm}.DHA', h['DHA'], 'deg', 1, 'derived/d08_periodic.json:structures.relax_i.hbonds', f"periodic: D–H···A of {h['name']}", case='periodic_qe')
    add(f'periodic.hb.{nm}.dDA', h['dDA'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.hbonds', f"periodic − X-ray D···A of {h['name']}", case='periodic_qe')
add('periodic.hb.maxabs_dDA', max(abs(h['dDA']) for h in ri['hbonds']), 'A', 3, 'derived/d08_periodic.json', 'largest |ΔD···A| over the nine crystal H-bonds', case='periodic_qe')
o = ri['fu1_O4H4AO5']
add('periodic.tau_C3C4S1N2', D3['models']['periodic_i']['proton']['tau_C3_C4_S1_N2'], 'deg', 1, 'derived/d03_geometry.json:models.periodic_i.proton.tau_C3_C4_S1_N2 (chained minimum image)',
    'periodic: C3–C4–S1–N2 torsion (X-ray 102.7°)', case='periodic_qe')
add('periodic.O4H4A', o['OH'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.fu1_O4H4AO5.OH', 'periodic: O4–H4A (DFT position of the restrained hydrogen)', case='periodic_qe')
add('periodic.H4AO5', o['HO'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.fu1_O4H4AO5.HO', 'periodic: H4A···O5', case='periodic_qe')
add('periodic.O4O5', o['OO'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.fu1_O4H4AO5.OO', 'periodic: O4···O5 (X-ray 2.500 Å)', case='periodic_qe')
add('periodic.O4H4A_frac', o['frac'], None, 2, 'derived/d08_periodic.json:structures.relax_i.fu1_O4H4AO5.frac', 'periodic: H4A fractional position along O4→O5 (X-ray restrained 0.30)', case='periodic_qe')
add('periodic.O3H3A', ri['hbonds'][5]['DH'], 'A', 3, 'derived/d08_periodic.json:structures.relax_i.hbonds[O3-H3A...O1].DH', 'periodic: O3–H3A', case='periodic_qe')
conv = D8['convergence']
add('periodic.conv_ecut_dE', [r_['dE_meV_atom'] for r_ in conv['ecut']], 'meV', 2, 'periodic_qe/convergence.json:ecut[*].dE_meV_atom', 'energy change per atom vs the 70/560 Ry reference at 50/400, 60/480, 70/560 Ry (k 3×2×2)')
add('periodic.conv_k_dE', [r_['dE_meV_atom'] for r_ in conv['kgrid']], 'meV', 3, 'periodic_qe/convergence.json:kgrid[*].dE_meV_atom', 'energy change per atom vs 5×2×3 at 2×1×1, 3×2×2, 4×2×2, 5×2×3 (50/400 Ry)')
add('periodic.conv_rule', conv['chosen']['rule'], None, None, 'periodic_qe/convergence.json:chosen.rule', 'convergence rule (1 meV/atom, 1e-3 Ry/Bohr)')

# =====================================================================================================================
section('18. Docking: validation (gate G2)',
        'From dock/analysis/summary.json and gate_G2.json (validation runs, ETKDG ligand geometries). Symmetry-'
        'corrected heavy-atom RMSD of the top pose in the receptor frame, mean ± SD over 10 seeds. Headline numbers '
        're-derived from the raw out.pdbqt files (derive/d09_dock.py).')
vs = J(DOCK / 'analysis/summary.json')
gate = J(DOCK / 'analysis/gate_G2.json')
VAL = [('T5_7S2K/dhp_anion__vina', '7S2K', 'dihydropteroate (78H)', 'Vina', 'redock'),
       ('T5_7S2K/dhp_anion__vinardo', '7S2K', 'dihydropteroate (78H)', 'Vinardo', 'redock'),
       ('T4_3HS4/sfa_neutral__ad4zn', '3HS4', 'sulfanilamide (6RL9 SAN carried over)', 'AD4Zn', 'cross-dock'),
       ('T3_6RL9/sfa_neutral__ad4zn', '6RL9', 'sulfanilamide (SAN)', 'AD4Zn', 'redock'),
       ('T3_6RL9/sfa_anion__ad4zn', '6RL9', 'sulfonamidate (SAN)', 'AD4Zn', 'redock'),
       ('T3_6RL9/sfa_neutral__vina', '6RL9', 'sulfanilamide (SAN)', 'Vina', 'redock'),
       ('T1_3TZF/smx_anion__vina', '3TZF', 'sulfamethoxazole anion (08D)', 'Vina', 'redock'),
       ('T4_3HS4/azm_neutral__ad4zn', '3HS4', 'acetazolamide (AZM)', 'AD4Zn', 'redock'),
       ('T4_3HS4/azm_neutral__vinardo', '3HS4', 'acetazolamide (AZM)', 'Vinardo', 'redock'),
       ('T2_1AJ0/sfa_neutral__vina', '1AJ0', 'sulfanilamide (SAN)', 'Vina', 'redock')]
for job, pdb, lig, prog, kind in VAL:
    s = vs[job]['stats']
    g = gate['criteria'].get(job, {})
    base = f"dock.val.{job.replace('/', '.')}"
    der = D9['rmsd'].get(job)
    add(base + '.rmsd', s['rmsd_top_mean'], 'A', 2, f'dock/analysis/summary.json:{job}.stats.rmsd_top_mean', f'validation {kind} {pdb}, {lig}, {prog}: top-pose RMSD mean', check=der['mean'] if der else None, tol=0.006)
    add(base + '.rmsd_sd', s['rmsd_top_sd'], 'A', 2, f'dock/analysis/summary.json:{job}.stats.rmsd_top_sd', f'validation {kind} {pdb}, {lig}, {prog}: RMSD SD over seeds', check=der['sd'] if der else None, tol=0.006)
    add(base + '.f1', s['success_fraction_1A'], None, 1, f'dock/analysis/summary.json:{job}.stats.success_fraction_1A', f'validation {kind} {pdb}, {lig}, {prog}: fraction of seeds with top pose ≤ 1 Å', check=der['f_le_1'] if der else None, tol=0.001)
    add(base + '.f2', s['success_fraction_2A'], None, 1, f'dock/analysis/summary.json:{job}.stats.success_fraction_2A', f'validation {kind} {pdb}, {lig}, {prog}: fraction ≤ 2 Å')
    add(base + '.score', s['score_mean'], 'kcal mol-1', 2, f'dock/analysis/summary.json:{job}.stats.score_mean', f'validation {kind} {pdb}, {lig}, {prog}: mean top score')
    add(base + '.pb', s['posebusters_physics_pass_fraction'], None, 1, f'dock/analysis/summary.json:{job}.stats.posebusters_physics_pass_fraction', f'validation {kind}: PoseBusters physics pass fraction')
    if g:
        add(base + '.G2', 'PASS' if g['passes_G2'] else 'fail', None, None, f'dock/analysis/gate_G2.json:{job}.passes_G2', f'validation {kind} {pdb}, {prog}: gate G2 (≤ 1.0 Å)')
        add(base + '.crit2A', 'pass' if g['passes_2A_criterion'] else 'fail', None, None, f'dock/analysis/gate_G2.json:{job}.passes_2A_criterion', f'validation {kind} {pdb}, {prog}: conventional 2.0 Å criterion')
    sc = s.get('site_checks') or {}
    for chk, v in sc.items():
        if chk.startswith('Zn-N') or chk.startswith('Thr199 OG1'):
            nm = 'ZnN' if chk.startswith('Zn') else 'Thr199'
            add(base + f'.{nm}', v['distance_mean'], 'A', 2, f'dock/analysis/summary.json:{job}.stats.site_checks["{chk}"]', f'validation {kind} {pdb}, {prog}: mean {chk} distance (pass fraction {v["pass_fraction"]})')
add('dock.val.n_pass_G2', 3, None, 0, 'dock/VALIDATION.md section 4; gate_G2.json', 'distinct protocols passing the 1.0 Å gate (7S2K Vina, 7S2K Vinardo, 3HS4 sulfanilamide AD4Zn cross-dock)')
add('dock.val.crystal_ZnN_6RL9', 1.93, 'A', 2, 'dock/VALIDATION.md section 2 (6RL9 SAN304)', 'crystal Zn–N(sulfonamide) in 6RL9')
add('dock.val.crystal_Thr199_6RL9', 2.86, 'A', 2, 'dock/VALIDATION.md section 2', 'crystal Thr199 OG1···N in 6RL9')
add('dock.val.crystal_ZnN_3HS4', 1.94, 'A', 2, 'dock/VALIDATION.md section 2 (3HS4 AZM701)', 'crystal Zn–N(sulfonamide) in 3HS4')
add('dock.val.tz_offset_6RL9', 0.17, 'A', 2, 'dock/VALIDATION.md section 2', 'AD4Zn TZ pseudo-atom distance from the crystal sulfonamide N (6RL9; 0.07 Å in 3HS4)')
add('dock.val.ca_rmsd_6RL9_3HS4', 0.328, 'A', 3, 'dock/analysis/summary.json:T4_3HS4/sfa_neutral__ad4zn.reference.superposition_ca_rmsd_A', 'Cα superposition RMSD 6RL9↔3HS4 (257 Cα) for the cross-dock references')
add('dock.val.ad4zn_gap_6RL9', '0.000–0.046 (mean 0.016)', 'kcal mol-1', None, 'dock/VALIDATION.md section 3 (6RL9 AD4Zn bimodality)', 'AD4Zn score gap between the near-native and the ring-rotated pose in validation')
rk = gate['ranking']
for key_, lab in [('T3_6RL9/ad4zn', '6RL9 AD4Zn'), ('T4_3HS4/ad4zn', '3HS4 AD4Zn'), ('T3_6RL9/vina', '6RL9 Vina'), ('T4_3HS4/vina', '3HS4 Vina')]:
    add(f"dock.val.rank.{key_.replace('/', '.')}.gap", rk[key_]['all_seeds']['gap_kcal'], 'kcal mol-1', 2, f'dock/analysis/gate_G2.json:ranking["{key_}"].all_seeds.gap_kcal',
        f'{lab}: acetazolamide scores better than sulfanilamide by this margin (all seeds; literature order reproduced)')
add('dock.val.rank_note', 'acetazolamide ranks above sulfanilamide in 4/4 receptor × scoring combinations by score, but in both Vina combinations no seed gives a zinc-bound pose for either ligand', None, None,
    'dock/VALIDATION.md section 5; gate_G2.json:ranking', 'ranking test caveat')

# =====================================================================================================================
section('19. Docking: production (DFT geometries of the sulfanilamide species)',
        'dock/results/production_summary.json; 72 production jobs. Scores are docking scores in kcal mol⁻¹ (never '
        'affinities), mean ± SD over seeds; ranks by mean score among all species docked to that receptor with that '
        f"scoring function, no site filter. Our parse of every out.pdbqt reproduces all {D9['n_jobs_compared']} job means to "
        f"{D9['max_abs_score_mean_diff_vs_summary']:.4f} kcal/mol.")
ps = J(DOCK / 'results/production_summary.json')
SPN = {'sfa_neutral': 'sulfanilamide', 'sfa_anion': 'sulfonamidate', 'sfa_cation': 'cation', 'smx_neutral': 'SMX⁰', 'smx_anion': 'SMX⁻',
       'azm_neutral': 'AZM⁰', 'azm_anion': 'AZM⁻', 'paba_anion': 'pABA⁻', 'dhp_anion': 'DHP (78H)'}
for job, j in ps['jobs'].items():
    t, tag = job.split('/')
    if j.get('diagnostic') or not tag.endswith('__prod'):
        continue
    sp, scoring = tag.split('__')[0], tag.split('__')[1]
    base = f"dock.prod.{t.split('_')[1]}.{sp}.{scoring}"
    der = D9['scores'].get(job)
    add(base + '.score', j['score']['mean'], 'kcal mol-1', 1, f'dock/results/production_summary.json:jobs["{job}"].score.mean',
        f"{t.split('_')[1]}, {SPN.get(sp, sp)}, {scoring}: mean top-pose score over {j['n_seeds']} seeds", check=der['mean'] if der else None, tol=0.0006)
    add(base + '.sd', j['score']['sd'], 'kcal mol-1', 2, f'dock/results/production_summary.json:jobs["{job}"].score.sd',
        f"{t.split('_')[1]}, {SPN.get(sp, sp)}, {scoring}: SD over seeds", check=der['sd'] if der else None, tol=0.0006)
for key_, v in ps['rank_order'].items():
    t, sc = key_.split('/')
    rkd = D9['ranks'][key_]
    rank = [x for x in v['ranking'] if x['species'] == 'sfa_neutral'][0]['rank']
    add(f"dock.prod.{t.split('_')[1]}.rank_sfa.{sc}", rank, None, 0, f'dock/results/production_summary.json:rank_order["{key_}"]',
        f"{t.split('_')[1]} {sc}: rank of neutral sulfanilamide among {len(v['ranking'])} species (no site filter)", check=rkd['rank_sfa_neutral'], tol=0)
    add(f"dock.prod.{t.split('_')[1]}.n_species.{sc}", len(v['ranking']), None, 0, f'dock/results/production_summary.json:rank_order["{key_}"]', f"{t.split('_')[1]} {sc}: species ranked")
# headline RMSDs
hl = [('T3_6RL9/sfa_neutral__ad4zn__prod', 'SAN (6RL9, own)', 'full', 'dock.prod.6RL9.sfa.ad4zn.rmsd', '6RL9 AD4Zn, sulfanilamide (DFT geometry): RMSD to the co-crystal SAN'),
      ('T4_3HS4/sfa_neutral__ad4zn__prod', 'SAN (6RL9, carried into 3HS4 by Ca superposition)', 'full', 'dock.prod.3HS4.sfa.ad4zn.rmsd', '3HS4 AD4Zn, sulfanilamide (DFT geometry): RMSD to 6RL9 SAN carried into 3HS4 (ring rotated)'),
      ('T1_3TZF/sfa_neutral__vina__prod', '08D (3TZF, own)', 'aniline-SO2N core (11 atoms)', 'dock.prod.3TZF.sfa.vina.core_rmsd', '3TZF Vina, sulfanilamide: RMSD of its aniline-SO2N core to the same core of co-crystal sulfamethoxazole (08D)'),
      ('T5_7S2K/sfa_neutral__vina__prod', '78H (7S2K, own)', 'aniline-X(O)2 core, X = S or C (10 atoms)', 'dock.prod.7S2K.sfa.vina.core_rmsd', '7S2K Vina, sulfanilamide: RMSD of the aniline-S(O)2 core to the aniline-C(O)2 core of 78H')]
for job, refname, meas, key_, what in hl:
    j = ps['jobs'][job]
    refk = [k for k in j['rmsd_to_crystal'] if k.startswith(refname.split(' (')[0])][0]
    rr = j['rmsd_to_crystal'][refk]
    blk = rr['full'] if meas == 'full' else rr['cores'][meas]
    der = D9['rmsd'].get(job)
    val_, note_ = blk['mean'], None
    # production_summary.json stores means rounded to 3 decimals; when that value sits exactly on the rounding tie at the
    # printed 2 decimals (3.215 for 3HS4), the full-precision re-derivation (derive/d09_dock.py, same poses) decides.
    if der and abs(val_ * 1000 - round(val_ * 1000)) < 1e-6 and round(val_ * 1000) % 10 == 5:
        note_ = (f'production_summary.json stores {val_:.3f} (mean of per-seed values rounded to 3 decimals, a tie at 2 decimals); '
                 f'value is the full-precision mean of derive/d09_dock.py over the same 20 poses, {der["mean"]:.5f}, so it prints {der["mean"]:.2f}')
        val_ = der['mean']
    add(key_, val_, 'A', 2, f'dock/results/production_summary.json:jobs["{job}"].rmsd_to_crystal["{refk}"].{"full" if meas == "full" else "cores[" + meas + "]"}.mean',
        what, check=der['mean'] if der else None, tol=0.006, note=note_)
    add(key_ + '_sd', blk['sd'], 'A', 2, f'dock/results/production_summary.json:jobs["{job}"]...sd', what + ' (SD over seeds)', check=der['sd'] if der else None, tol=0.006)
    add(key_ + '_f1', blk['f_le_1A'], None, 2, f'dock/results/production_summary.json:jobs["{job}"]...f_le_1A', what + ' (fraction of seeds ≤ 1 Å)', check=der['f_le_1'] if der else None, tol=0.001)
    add(key_ + '_n', blk['n'], None, 0, f'dock/results/production_summary.json:jobs["{job}"]...n', what + ' (seeds)')
j = ps['jobs']['T3_6RL9/sfa_neutral__ad4zn__prod']
add('dock.prod.6RL9.sfa.ad4zn.ZnN', j['site_checks']['Zn-N(sulfonamide)']['distance_mean'], 'A', 2, 'production_summary.json:jobs[T3_6RL9/sfa_neutral__ad4zn__prod].site_checks', '6RL9 AD4Zn sulfanilamide: mean Zn–N (crystal 1.93 Å; pass 20/20)')
add('dock.prod.6RL9.sfa.ad4zn.Thr199', j['site_checks']['Thr199 OG1 H-bond']['distance_mean'], 'A', 2, 'production_summary.json:jobs[...].site_checks', '6RL9 AD4Zn sulfanilamide: mean Thr199 OG1···N (crystal 2.86 Å; pass 20/20)')
add('dock.prod.6RL9.sfa.ad4zn.zinc_bound', f"{j['zinc_bound']['n_seeds']}/{j['zinc_bound']['of']}", None, None, 'production_summary.json:jobs[...].zinc_bound', '6RL9 AD4Zn sulfanilamide: seeds with a zinc-bound top pose')
j = ps['jobs']['T4_3HS4/sfa_neutral__ad4zn__prod']
add('dock.prod.3HS4.sfa.ad4zn.ZnN', j['site_checks']['Zn-N(sulfonamide)']['distance_mean'], 'A', 2, 'production_summary.json:jobs[T4_3HS4/sfa_neutral__ad4zn__prod].site_checks',
    f"3HS4 AD4Zn sulfanilamide: mean Zn–N (pass fraction {j['site_checks']['Zn-N(sulfonamide)']['pass_fraction']})")
add('dock.prod.3HS4.sfa.ad4zn.zinc_bound', f"{j['zinc_bound']['n_seeds']}/{j['zinc_bound']['of']}", None, None, 'production_summary.json:jobs[...].zinc_bound', '3HS4 AD4Zn sulfanilamide: seeds with a zinc-bound top pose')
for job, j in ps['jobs'].items():
    if j['scoring'] != 'ad4zn' or j.get('diagnostic') or not job.endswith('__prod'):
        continue
    t_ = job.split('/')[0].split('_')[1]
    for chk, nm in [('Zn-N(sulfonamide)', 'ZnN'), ('Thr199 OG1 H-bond', 'Thr199')]:
        sc_ = j['site_checks'][chk]
        k_ = f"dock.site.{t_}.{j['ligand']}.{nm}"
        if k_ in FACTS:
            continue
        add(k_, sc_['distance_mean'], 'A', 2, f'dock/results/production_summary.json:jobs["{job}"].site_checks["{chk}"].distance_mean',
            f"{t_} AD4Zn {SPN.get(j['ligand'], j['ligand'])}: mean {chk} distance of the top poses (pass fraction {sc_['pass_fraction']}; crystal {sc_['crystal_ligand_value']} Å)")
        add(k_ + '.pass', sc_['pass_fraction'], None, 2, f'dock/results/production_summary.json:jobs["{job}"].site_checks["{chk}"].pass_fraction',
            f"{t_} AD4Zn {SPN.get(j['ligand'], j['ligand'])}: fraction of seeds passing {chk}")
add('dock.prod.ad4zn_orientation_gap', '0.000–0.219', 'kcal mol-1', None, 'dock/RESULTS.md section 7b', 'AD4Zn score gap between near-native and ring-rotated sulfanilamide poses across all eight AD4Zn runs')
add('dock.prod.csn_dft', ps['csn_angle']['dft']['sfa_neutral']['C_S_N_deg'], 'deg', 2,
    'dock/results/production_summary.json:csn_angle.dft.sfa_neutral; dock/RESULTS.md section 7b', 'C4–S1–N2 angle of the DFT sulfanilamide (ETKDG/MMFF94s 101.4–102.0°; crystal SAPH 109.2°)')
add('dock.prod.posebusters', 0.993, None, 3, 'dock/VALIDATION.md section 14 / RESULTS.md section 8', 'PoseBusters pass fraction over 740 production top poses (failures: 5 Zn-coordination "clashes" of AZM⁻ Vinardo in 6RL9)')
add('dock.prod.n_jobs', 72, None, 0, 'dock/RESULTS.md (What ran)', 'production docking jobs (plus 6 diagnostic AD4Zn runs)')
# torsions and strain
cr = ps['crystal']['torsion']
add('dock.tors.crystal', cr['tau_C3'], 'deg', 1, 'production_summary.json:crystal.torsion.tau_C3', 'C3–C4–S1–N2 in the crystal (τ_fold 102.7°)')
for sp in ['sfa_neutral', 'sfa_anion', 'sfa_cation']:
    add(f'dock.tors.dft.{sp}', ps['dft'][sp]['tau_fold'], 'deg', 1, f'production_summary.json:dft.{sp}.tau_fold', f'τ_fold of the DFT (IEF-PCM) {SPN[sp]}', check=D3['fragments'][{'sfa_neutral': 'sulfanilamide_pcm_d3', 'sfa_anion': 'sulfonamidate_pcm_d3', 'sfa_cation': 'cation_pcm_d3'}[sp]]['tau_fold'], tol=0.05)
tf_all, st_all = [], []
for job, j in ps['jobs'].items():
    if j['ligand'].startswith('sfa_') and j.get('torsion') and not j.get('diagnostic'):
        tf_all += [x['tau_fold'] for x in j['torsion']['per_seed']]
        st_all += [x['strain_mmff_relaxed'] for x in j['torsion']['per_seed']]
add('dock.tors.range_min', min(tf_all), 'deg', 1, 'production_summary.json:jobs[sfa_*].torsion.per_seed.tau_fold (own min over all)', f'smallest τ_fold over all {len(tf_all)} docked sulfanilamide-species top poses')
add('dock.tors.range_max', max(tf_all), 'deg', 1, 'production_summary.json:jobs[sfa_*].torsion.per_seed.tau_fold (own max)', 'largest τ_fold over the same poses')
add('dock.tors.n_poses', len(tf_all), None, 0, 'production_summary.json:jobs[sfa_*].torsion.per_seed', 'docked top poses of the three sulfanilamide species with a torsion')
add('dock.tors.strain_max', max(st_all), 'kcal mol-1', 2, 'production_summary.json:jobs[sfa_*].torsion.per_seed.strain_mmff_relaxed (own max)', 'largest MMFF94s relaxed strain of those poses (B3LYP-D3(BJ) single points of the cluster representatives: section 29)')
j = ps['jobs']['T3_6RL9/sfa_neutral__ad4zn__prod']['torsion']
add('dock.tors.6RL9_ad4zn_mean', j['tau_fold']['mean'], 'deg', 1, 'production_summary.json:jobs[T3_6RL9/sfa_neutral__ad4zn__prod].torsion.tau_fold.mean', '6RL9 AD4Zn sulfanilamide: mean τ_fold (range 140.6–151.8°)')
add('dock.tors.6RL9_ad4zn_strain', j['strain_mmff_relaxed']['mean'], 'kcal mol-1', 2, 'production_summary.json:...strain_mmff_relaxed.mean', '6RL9 AD4Zn sulfanilamide: mean MMFF94s relaxed strain')
# contacts of the headline poses
pc = J(DOCK / 'results/production_contacts.json')
for job, lab in [('T3_6RL9/sfa_neutral__ad4zn__prod', '6RL9 AD4Zn sulfanilamide'), ('T5_7S2K/sfa_neutral__vina__prod', '7S2K Vina sulfanilamide'),
                 ('T1_3TZF/sfa_neutral__vina__prod', '3TZF Vina sulfanilamide'), ('T4_3HS4/sfa_neutral__ad4zn__prod', '3HS4 AD4Zn sulfanilamide')]:
    c1 = pc[job]['1']
    res4 = ', '.join(f'{k} {v:.2f}' for k, v in list(c1['residues_within_4A'].items())[:8])
    polar = '; '.join(f"{p['residue']} {p['receptor_atom']}···{p['ligand_atom'].split('#')[0]} {p['distance']:.2f}" for p in c1['polar_pairs_NO_within_3.5A'][:5])
    add(f"dock.contacts.{job.split('/')[0].split('_')[1]}.{job.split('__')[1]}.residues", res4, 'A', None, f'dock/results/production_contacts.json:["{job}"]["1"].residues_within_4A',
        f'{lab}: residues within 4 Å of the cluster-1 representative (closest heavy-atom distance)')
    add(f"dock.contacts.{job.split('/')[0].split('_')[1]}.{job.split('__')[1]}.polar", polar, 'A', None, f'dock/results/production_contacts.json:["{job}"]["1"].polar_pairs_NO_within_3.5A',
        f'{lab}: polar N/O···N/O pairs ≤ 3.5 Å')
# caveats (text facts)
for i, cv in enumerate([
        'Docking scores are empirical scores (Vina, Vinardo, AutoDock4Zn), not binding free energies; never converted to Ki/Kd; only the order against the controls is read.',
        'Rigid receptors (one crystal conformation), ligands flexible in torsions only (bond lengths/angles fixed at the input geometry); no explicit water except the Mg coordination spheres (3 in 7S2K, 3 in 3TZF) and, in one diagnostic outside the protocol, four crystal waters of 6RL9 (section 40). No free energy is computed. The hCA II–sulfanilamide complex was simulated by molecular dynamics (section 40: nine runs of 100 ns and two sensitivity runs, with force-field terms at the sulfonamide nitrogen that failed two quantum-chemical tests); the simulations give no rate, barrier or free energy and replace no docking number.',
        'No experimental affinity or bioassay for SAPH; validation is structural only; acetazolamide above sulfanilamide at hCA II is the only external anchor.',
        '1AJ0 failed validation (sulfanilamide 11.6 Å away in every seed) and is excluded from production.',
        '3TZF passes only the 2.0 Å criterion (SMX anion 1.92 ± 0.04 Å): 3TZF poses are qualitative, the 0.20 Å core RMSD of sulfanilamide notwithstanding.',
        '6RL9 AD4Zn failed the strict gate in validation (2.07 ± 1.45 Å, bimodal); with the DFT input geometry production gives 0.61 ± 0.04 Å (20/20), but the ring orientation flips with the input geometry (3HS4 cross-dock moves from 0.62 to 3.22 Å): report with the sensitivity and the validation figure, not as a passed gate.',
        'Plain Vina/Vinardo have no metal term: no zinc-bound top pose in either hCA II receptor; their hCA II scores are a scoring-function comparison, not a binding mode.',
        'His64 tautomer (HIE) is an assumption (≤ 0.03 kcal/mol effect with Vina; untested with AD4Zn). Species fixed at pH 7.4; the ion pair is never docked.',
        'LE/FQ are size-biased (r(score, heavy atoms) = −0.87 Vina) and are not used to claim potency.']):
    add(f'dock.caveat.{i + 1}', cv, None, None, 'dock/RESULTS.md section 8; dock/VALIDATION.md sections 3–10', 'docking caveat to carry into the paper')

# =====================================================================================================================
section('20. Sensitivity of the proton position to dispersion and basis (summary)',
        'Only the facts not already in sections 4–5.')
add('sens.b3lyp.C6_O4_open', D3['models']['pcm_b3lyp']['proton']['d_C6_O4'] - D3['models']['pcm']['proton']['d_C6_O4'], 'A', 2, 'derived/d03_geometry.json',
    'opening of C6···O4 without D3 (IEF-PCM)', case='sfa_pair_pcm_b3lyp')
add('sens.tz.delta_shift', D3['models']['pcm_tz']['proton']['delta'] - D3['models']['pcm']['proton']['delta'], 'A', 3, 'derived/d03_geometry.json',
    'δ(def2-TZVP) − δ(6-311++G(d,p)) in IEF-PCM: the larger basis moves H1A 0.09 Å towards O5, still on N1', case='sfa_pair_pcm_d3_tz')
add('sens.tz.nao', 574, None, 0, 'derived/d02_models.json:sfa_pair_pcm_d3_tz.nao_spherical', 'basis functions of the def2-TZVP pair', case='sfa_pair_pcm_d3_tz')
add('sens.tz.steps', 47, None, 0, 'sfa_pair_pcm_d3_tz/opt_summary.json:steps', 'def2-TZVP optimisation steps (converged normally)', case='sfa_pair_pcm_d3_tz')
FR = D3['fragments']
for c, lab in [('sulfanilamide_pcm_d3', 'sulfanilamide (IEF-PCM)'), ('sulfanilamide_pcm_d3_tz', 'sulfanilamide (IEF-PCM, def2-TZVP)'),
               ('sulfonamidate_pcm_d3', 'sulfonamidate (IEF-PCM)'), ('cation_pcm_d3', 'cation (IEF-PCM)'), ('cation_gas_d3', 'cation (gas)')]:
    f = FR[c]
    add(f'frag.{c}.heavy_mae', f['heavy_bond_stats']['mae'], 'A', 3, f'derived/d03_geometry.json:fragments.{c}.heavy_bond_stats', f'{lab}: heavy-atom bond MAE vs the same atoms in the X-ray cation', case=c)
    add(f'frag.{c}.tau_fold', f['tau_fold'], 'deg', 1, f'derived/d03_geometry.json:fragments.{c}.tau_fold', f'{lab}: τ_fold = max(|C3–C4–S1–N2|, |C5–C4–S1–N2|)', case=c)
    add(f'frag.{c}.ang_CSN', f['ang_C4_S1_N2'], 'deg', 1, f'derived/d03_geometry.json:fragments.{c}.ang_C4_S1_N2', f'{lab}: C4–S1–N2 angle', case=c)
    add(f'frag.{c}.C1N1', f['d_C1_N1'], 'A', 3, f'derived/d03_geometry.json:fragments.{c}.d_C1_N1', f'{lab}: C1–N1', case=c)

# =====================================================================================================================
section('21. Bulk data behind the SI tables',
        'Every number printed in the SI tables that is not already a fact above, so that a sentence quoting one can be '
        'traced. Listed in facts.json with "bulk": true; summarised here by key prefix only.')
BULK_PREFIX = {}


def addb(key, value, unit, prec, source, meaning, case=None):
    add(key, value, unit, prec, source, meaning, case=case)
    FACTS[key]['bulk'] = True
    parts = key.split('.')
    n = {'bond': 3, 'angdev': 3, 'td': 4, 'chg': 3}.get(parts[1], 2)
    pre = '.'.join(parts[:n])
    what = {'freq': 'harmonic wavenumbers, scaled wavenumbers and IR intensities of the 75 modes (Table S5)',
            'mep': 'MEP extrema of the IEF-PCM pair (Table S17)', 'dockval': 'all docking validation runs (Table S25)',
            'dockprod': 'all production scores, ranks, LE, FQ (Table S26)', 'docktors': 'torsion and strain ranges of the docked poses (Table S27)',
            'dockcont': 'polar contacts and Zn distances of the cluster representatives (Table S28)', 'hs': 'Hirshfeld red spots and contact shares (Tables S20, S21)',
            'opt': 'energies, steps and ZPE of all optimised structures (Table S2)', 'qeconv': 'periodic convergence tests (Table S12)'}.get(parts[1])
    if what is None:
        what = meaning.split(':')[0] + {'bond': ': bond lengths (Table S3)', 'angdev': ': heavy-atom angle deviations (Table S4)',
                                       'td': ': TD-DFT states (Tables S22, S23)', 'chg': ' charges per atom (Table S16)'}[parts[1]]
    BULK_PREFIX.setdefault(pre, [0, source.split(':')[0].split('[')[0], what])
    BULK_PREFIX[pre][0] += 1


XRB = D3['models']['xray']['bond_lengths']
for mk, lab, case in [('xray', 'X-ray', None), ('pcm', 'IEF-PCM', 'sfa_pair_pcm_d3'), ('gas', 'gas phase', 'sfa_pair_gas_d3'),
                      ('pcm_b3lyp', 'IEF-PCM no D3', 'sfa_pair_pcm_b3lyp'), ('pcm_tz', 'IEF-PCM def2-TZVP', 'sfa_pair_pcm_d3_tz'),
                      ('periodic_i', 'periodic', 'periodic_qe')]:
    for b, x in D3['models'][mk]['bond_lengths'].items():
        addb(f'bulk.bond.{mk}.{b}', x, 'A', 3, f'derived/d03_geometry.json:models.{mk}.bond_lengths', f'{lab}: bond {b}', case=case)
    if mk != 'xray':
        for a_, x in D3['models'][mk]['angle_devs'].items():
            addb(f'bulk.angdev.{mk}.{a_}', x, 'deg', 1, f'derived/d03_geometry.json:models.{mk}.angle_devs', f'{lab}: angle {a_} minus X-ray', case=case)
frj = J(RES / 'sfa_pair_pcm_d3/freq.json')
for m in D5['modes']:
    fj = frj['ir'][m['n'] - 1]
    addb(f'bulk.freq.{m["n"]:02d}.nu', fj['freq_cm-1'], 'cm-1', 0, f'sfa_pair_pcm_d3/freq.json:ir[{m["n"] - 1}]', f'mode {m["n"]}: harmonic wavenumber', case='sfa_pair_pcm_d3')
    addb(f'bulk.freq.{m["n"]:02d}.scaled', fj['freq_cm-1'] * 0.9688, 'cm-1', 0, f'sfa_pair_pcm_d3/freq.json:ir[{m["n"] - 1}] x 0.9688', f'mode {m["n"]}: scaled wavenumber', case='sfa_pair_pcm_d3')
    addb(f'bulk.freq.{m["n"]:02d}.I', fj['intensity_km_mol'], 'km mol-1', 1, f'sfa_pair_pcm_d3/freq.json:ir[{m["n"] - 1}]', f'mode {m["n"]}: IR intensity', case='sfa_pair_pcm_d3')
chj = J(RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.json')
labs_ = [l.split(',')[1].replace('#1_554', '') for l in (RES / 'sfa_pair_pcm_d3/mwfn_extra/charges.csv').read_text().splitlines()[1:]]
for sch in ['NPA', 'Hirshfeld', 'ADCH', 'CM5']:
    for i, a_ in enumerate(labs_):
        addb(f'bulk.chg.{sch}.{a_}', chj['atomic'][sch][i], 'e', 3, f'sfa_pair_pcm_d3/mwfn_extra/charges.json:atomic.{sch}[{i}]', f'IEF-PCM pair: {sch} charge of {a_}', case='sfa_pair_pcm_d3')
for i, e in enumerate(J(RES / 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json')['extrema']):
    addb(f'bulk.mep.{i:02d}', e['ESP_kcal_mol'], 'kcal mol-1', 1, f'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json:extrema[{i}]', f"IEF-PCM pair: MEP {e['type']} near {e['nearest_atoms']}", case='sfa_pair_pcm_d3')
    addb(f'bulk.mep.{i:02d}.au', e['ESP_au'], 'a.u.', 4, f'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json:extrema[{i}]', f"IEF-PCM pair: MEP {e['type']} (a.u.)", case='sfa_pair_pcm_d3')
for key_, t in D6['tddft'].items():
    c, fn = key_.split('/')
    tag = fn.replace('tddft_', '')
    for st in t['states']:
        b = f"bulk.td.{c.split('_')[0]}.{tag}.S{st['state']:02d}"
        addb(b + '.E', st['E_eV'], 'eV', 3, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}]', f'{c} {tag}: S{st["state"]} energy', case=c)
        addb(b + '.nm', st['nm'], 'nm', 0, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}]', f'{c} {tag}: S{st["state"]} wavelength', case=c)
        addb(b + '.f', st['f'], None, 4, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}]', f'{c} {tag}: S{st["state"]} oscillator strength', case=c)
        if st['CT'] is not None:
            addb(b + '.CT', st['CT'], None, 2, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}].ct', f'{c} {tag}: S{st["state"]} CT number', case=c)
            addb(b + '.DCT', st['D_CT'], 'A', 2, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}].ct', f'{c} {tag}: S{st["state"]} D_CT', case=c)
def _full_stats(s_):
    """Mean (exact decimal arithmetic) and SD (ddof 1) of the per-seed top-pose score and RMSD: summary.json stores them
    rounded to 3 decimals, which double-rounds the ties at 2 decimals (rounding rule of 29 Sep 2026, section 30)."""
    ok = [x for x in s_['per_seed'] if x.get('ok')]
    out_ = {}
    for q, fld in (('score', 'score_top'), ('rmsd_top', 'rmsd_top')):
        vals = [x[fld] for x in ok if x.get(fld) is not None]
        if vals and s_['stats'].get(q + '_mean') is not None:
            out_[q + '_mean'] = float(sum(_DecB(repr(x)) for x in vals) / len(vals))
            out_[q + '_sd'] = _stt.stdev(vals) if len(vals) > 1 else 0.0
    return out_


for job, s_ in vs.items():
    st = s_['stats']
    full_ = _full_stats(s_)
    b = f"bulk.dockval.{job.replace('/', '.')}"
    for q, u, pr in [('rmsd_top_mean', 'A', 2), ('rmsd_top_sd', 'A', 2), ('success_fraction_1A', None, 1), ('success_fraction_2A', None, 1),
                     ('score_mean', 'kcal mol-1', 2), ('score_sd', 'kcal mol-1', 2), ('posebusters_physics_pass_fraction', None, 1)]:
        if st.get(q) is not None:
            src_ = f'dock/analysis/summary.json:{job}.stats.{q}' + (f' (recomputed at full precision from per_seed; stored {st[q]})' if q in full_ else '')
            addb(f'{b}.{q}', full_.get(q, st[q]), u, pr, src_, f'validation run {job}: {q}')
            if q in full_ and abs(full_[q] - st[q]) > 1e-12:
                FACTS[f'{b}.{q}']['stored_value'] = st[q]
for key_, ro in ps['rank_order'].items():
    for rr in ro['ranking']:
        b = f"bulk.dockprod.{key_.replace('/', '.')}.{rr['species']}"
        # mean and SD at full precision from the per-seed top scores (exact 3-decimal numbers, so the mean is computed in
        # decimal arithmetic); production_summary.json stores both rounded to 3 decimals, which double-rounds the exact
        # ties (e.g. -9.2575). Rounding rule of 29 Sep 2026 (facts_burst2.py section 30).
        per_ = ps['jobs'][f"{key_.split('/')[0]}/{rr['species']}__{key_.split('/')[1]}__prod"]['score']['per_seed']
        full_ = {'score_mean': float(sum(_DecB(repr(x)) for x in per_) / len(per_)), 'score_sd': _stt.stdev(per_)}
        _ha = rr['heavy_atoms']
        full_['LE'] = -full_['score_mean'] / _ha                                    # dock/production_analysis.py: LE = -score/HA
        full_['FQ'] = full_['LE'] / (0.0715 + 7.5328 / _ha + 25.7079 / _ha ** 2 - 361.4722 / _ha ** 3)   # FQ = LE/LE_scale(HA)
        for q, u, pr in [('score_mean', 'kcal mol-1', 3), ('score_sd', 'kcal mol-1', 3), ('best', 'kcal mol-1', 3), ('rank', None, 0),
                         ('LE', None, 3), ('FQ', None, 2)]:
            src_ = f'dock/results/production_summary.json:rank_order["{key_}"]'
            if q in full_:
                src_ += f' (recomputed at full precision from jobs[...].score.per_seed; stored {rr[q]})'
            addb(f'{b}.{q}', full_.get(q, rr[q]), u, pr, src_, f'production {key_} {rr["species"]}: {q}')
            if q in full_ and abs(full_[q] - rr[q]) > 1e-12:
                FACTS[f'{b}.{q}']['stored_value'] = rr[q]
for job, j in ps['jobs'].items():
    if j['ligand'].startswith('sfa_') and j.get('torsion') and not j.get('diagnostic'):
        tt = j['torsion']
        for q in ['tau_fold', 'strain_mmff_rigid', 'strain_mmff_relaxed']:
            for mm in ['min', 'max', 'mean']:
                addb(f"bulk.docktors.{job.replace('/', '.')}.{q}.{mm}", tt[q][mm], 'deg' if q == 'tau_fold' else 'kcal mol-1', 1 if q == 'tau_fold' else 2,
                     f'dock/results/production_summary.json:jobs["{job}"].torsion.{q}.{mm}', f'{job}: {q} {mm}')
pcj_ = J(DOCK / 'results/production_contacts.json')
for job, cl in pcj_.items():
    if '__prod' not in job:
        continue
    c1 = cl.get('1', {})
    for i, p_ in enumerate(c1.get('polar_pairs_NO_within_3.5A', [])[:6]):
        addb(f"bulk.dockcont.{job.replace('/', '.')}.{i}", p_['distance'], 'A', 2, f'dock/results/production_contacts.json:["{job}"]["1"].polar_pairs_NO_within_3.5A[{i}]',
             f"{job}: {p_['residue']} {p_['receptor_atom']}···{p_['ligand_atom']}")
    zn = c1.get('zinc') or {}
    if zn.get('distance') is not None:
        addb(f"bulk.dockcont.{job.replace('/', '.')}.zn", zn['distance'], 'A', 2, f'dock/results/production_contacts.json:["{job}"]["1"].zinc', f"{job}: nearest ligand atom to Zn ({zn.get('nearest_ligand_atom')})")
hsj = J(RES / 'hirshfeld_xray/hirshfeld_summary.json')
for ion in ['cation', 'anion']:
    for i, rs in enumerate(hsj['targets'][ion]['red_spots']):
        addb(f'bulk.hs.{ion}.spot{i:02d}.dnorm', rs['dnorm'], None, 3, f'hirshfeld_xray/hirshfeld_summary.json:targets.{ion}.red_spots[{i}]', f"{ion}: red spot {rs['inside']}···{rs['outside']}")
        addb(f'bulk.hs.{ion}.spot{i:02d}.dist', rs['dist_A'], 'A', 3, f'hirshfeld_xray/hirshfeld_summary.json:targets.{ion}.red_spots[{i}]', f"{ion}: red spot distance {rs['inside']}···{rs['outside']}")
    for c_ in hsj['targets'][ion]['contacts_reciprocal']:
        key_ = f"bulk.hs.{ion}.pct.{c_['contact'].replace('/', '_')}"
        addb(key_, c_['pct'], '%', 1, f'hirshfeld_xray/hirshfeld_summary.json:targets.{ion}.contacts_reciprocal', f"{ion}: {c_['contact']} share")
for c in ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'sfa_pair_pcm_b3lyp', 'sfa_pair_pcm_d3_tz', 'cation_pcm_d3', 'cation_gas_d3', 'h2po4_pcm_d3', 'h2po4_gas_d3',
          'h3po4_pcm_d3', 'h3po4_gas_d3', 'sulfanilamide_pcm_d3', 'sulfanilamide_gas_d3', 'sulfanilamide_pcm_d3_tz', 'sulfonamidate_pcm_d3', 'urea_gas_d3']:
    o = J(RES / c / 'opt_summary.json')
    addb(f'bulk.opt.{c}.E', o['E_final_Ha'], 'Ha', 6, f'molecules/sfa/results/{c}/opt_summary.json:E_final_Ha', f'{c}: final energy', case=c)
    addb(f'bulk.opt.{c}.steps', o['steps'], None, 0, f'molecules/sfa/results/{c}/opt_summary.json:steps', f'{c}: optimisation steps', case=c)
    if (RES / c / 'freq.json').exists():
        addb(f'bulk.opt.{c}.zpe', J(RES / c / 'freq.json')['ZPE_Ha'], 'Ha', 6, f'molecules/sfa/results/{c}/freq.json:ZPE_Ha', f'{c}: ZPE', case=c)
_cv = D8['convergence']
for i, c_ in enumerate(_cv['ecut'] + _cv['kgrid']):
    # dE per atom at full precision from the stored total energies (convergence.json rounds it to 4 decimals); reference:
    # the 70/560 Ry point of the cutoff series and the 5x2x3 point of the k series (rounding rule, section 30)
    _ref = [x for x in (_cv['ecut'] if i < len(_cv['ecut']) else _cv['kgrid']) if x['dE_meV_atom'] == 0.0][0]
    _dE = (c_['E_Ry'] - _ref['E_Ry']) * 13605.693 / 54
    assert abs(_dE - c_['dE_meV_atom']) < 6e-5, (_dE, c_['dE_meV_atom'])
    addb(f'bulk.qeconv.{i}.dE', _dE, 'meV', 3, f'periodic_qe/convergence.json (E_Ry − reference E_Ry, per atom; stored {c_["dE_meV_atom"]})',
         f"QE convergence {c_['ecut_Ry']} {c_['k']}: dE per atom")
    addb(f'bulk.qeconv.{i}.dF', c_['max_dF_Ry_bohr'], 'Ry/bohr', 6, 'periodic_qe/convergence.json', f"QE convergence {c_['ecut_Ry']} {c_['k']}: max dF")
    addb(f'bulk.qeconv.{i}.dP', c_['dP_kbar'], 'kbar', 2, 'periodic_qe/convergence.json', f"QE convergence {c_['ecut_Ry']} {c_['k']}: dP")

# ---- facts for numbers that the SI tables print in captions/notes
section('22. Additional values printed in the SI tables', 'Values that the SI tables print in captions, notes or label columns.')
add('protocol.opt_fmax', 4.5e-4, 'Ha/bohr', 5, 'sfa_pair_pcm_d3/opt_summary.json:converged_on (Gaussian thresholds)', 'max-force convergence threshold (Gaussian default)')
add('protocol.opt_frms', 3.0e-4, 'Ha/bohr', 5, 'sfa_pair_pcm_d3/opt_summary.json:converged_on', 'RMS-force convergence threshold')
add('pcm.opt.frms_final', 7.35e-5, 'Ha/bohr', 7, 'sfa_pair_pcm_d3/opt_summary.json:converged_on', 'final RMS force of the IEF-PCM pair', case='sfa_pair_pcm_d3')
add('pcm.opt.fmax_final', 1.85e-4, 'Ha/bohr', 6, 'sfa_pair_pcm_d3/opt_summary.json:converged_on', 'final max force of the IEF-PCM pair', case='sfa_pair_pcm_d3')
add('sens.b3lyp.opt.frms_final', 5.92e-5, 'Ha/bohr', 7, 'sfa_pair_pcm_b3lyp/opt_summary.json:converged_on', 'final RMS force of the dispersion-free pair', case='sfa_pair_pcm_b3lyp')
add('sens.b3lyp.opt.fmax_final', 1.21e-4, 'Ha/bohr', 6, 'sfa_pair_pcm_b3lyp/opt_summary.json:converged_on', 'final max force of the dispersion-free pair', case='sfa_pair_pcm_b3lyp')
add('protocol.qe_conv_thr', 1e-8, 'Ry', 8, 'periodic_qe/periodic_qe.json:method.relax', 'QE SCF convergence threshold')
add('protocol.qe_ecut_tested', [50, 400, 60, 480, 70, 560], 'Ry', 0, 'periodic_qe/convergence.json:ecut[*].ecut_Ry', 'wavefunction/density cutoffs tested')
add('thermo.T', 298.15, 'K', 2, 'freq.json thermochemistry (PySCF thermo default)', 'temperature of the thermochemistry')
add('thermo.p_atm', 1, 'atm', 0, 'freq.json thermochemistry (101325 Pa)', 'pressure of the thermochemistry')
add('thermo.S_tr', D5['thermo']['S_tr'], 'J mol-1 K-1', 1, 'derived/d05_freq.json:thermo.S_tr', 'translational entropy (own RRHO)', status='own derivation')
add('thermo.S_rot', D5['thermo']['S_rot'], 'J mol-1 K-1', 1, 'derived/d05_freq.json:thermo.S_rot', 'rotational entropy (own RRHO)', status='own derivation')
add('thermo.S_vib', D5['thermo']['S_vib'], 'J mol-1 K-1', 1, 'derived/d05_freq.json:thermo.S_vib', 'vibrational entropy (own RRHO)', status='own derivation')
add('zpe.nh_character_pct', 100 * zb['nh_stretch']['normal_mode']['nh_stretch_character'], '%', 0, 'scan/zpe_barrier.json:nh_stretch.normal_mode.nh_stretch_character', 'N1–H1A stretch character of mode 65, in %')
for m in D5['modes']:
    addb(f'bulk.freq.{m["n"]:02d}.index', m['n'], None, 0, 'sfa_pair_pcm_d3/freq.json:ir (order)', f'mode number {m["n"]} (ascending wavenumber)')
for p in D4['points']:
    row = SCSV[(p['chain'], str(p['index']))]
    add(f"scan.{p['chain']}{p['index']:02d}.d_HO", float(row[5]), 'A', 3, f"molecules/sfa/results/sfa_pair_pcm_d3/scan/scan_profile.csv:{p['chain']},{p['index']} d_H1A_O5_A",
        f"scan chain {p['chain']} point {p['index']}: d(H1A···O5)", case='sfa_pair_pcm_d3', check=p['d_H1A_O5'], tol=0.0006)
for k in ['ii', 'iifix']:
    add(f'periodic.max_force_{k}', PQ['relax'][k]['max_force_Ry_bohr'], 'Ry/bohr', 5, f'periodic_qe.json:relax.{k}.max_force_Ry_bohr', f'largest residual force in ({k})', case='periodic_qe')
    add(f'periodic.rmsd_heavy_{k}', D8['structures'][f'relax_{k}']['rmsd_heavy_shift_removed'], 'A', 3, f'derived/d08_periodic.json:structures.relax_{k}.rmsd_heavy_shift_removed',
        f'heavy-atom RMSD of ({k}) vs X-ray, mean shift removed', case='periodic_qe', check=PQ['relax'][k]['rmsd_heavy_vs_xray']['rmsd_A'], tol=0.0006)
    add(f'periodic.E_{k}_rel_kJ', (PQ['relax'][k]['energy_Ry'] - PQ['relax']['i']['energy_Ry']) * 13.605693 * 96.485332 / 2, 'kJ mol-1', 2,
        f'periodic_qe.json:relax.{k}.energy_Ry − relax.i.energy_Ry', f'E({k}) − E(i) per formula unit', case='periodic_qe')
sh_ = D8['structures']['relax_iifix']['mean_heavy_shift_vs_start_ii']
add('periodic.iifix_shift_mag', (sum(x * x for x in sh_)) ** 0.5, 'A', 2, 'derived/d08_periodic.json:structures.relax_iifix.mean_heavy_shift_vs_start_ii', 'magnitude of the mean heavy-atom translation in iifix', case='periodic_qe')
for x_, v_ in [(k_, w_) for k_, w_ in D7['xh_lengths_A'].items() if k_ in ('C', 'N', 'O')]:
    add(f'hs.xh_{x_}', v_, 'A', 3, 'hirshfeld_xray/cluster.json:xh_lengths_A', f'neutron-normalised {x_}–H length used for the Hirshfeld surfaces')
for i, e in enumerate(J(RES / 'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json')['extrema']):
    addb(f'bulk.mep.{i:02d}.near', e['nearest_atoms'], 'A', None, f'sfa_pair_pcm_d3/mwfn_extra/mep_summary.json:extrema[{i}].nearest_atoms', 'nearest atoms (distances in Å) of a MEP extremum')
for c, pre in [('sfa_pair_pcm_d3', 'pcm'), ('sfa_pair_gas_d3', 'gas')]:
    for h in D6['qtaim'][c]['hbonds']:
        nm = {'N1-H1A...O5': 'N1H1A_O5', 'C6-H6...O4': 'C6H6_O4', 'N1-H1B...O4#1_554': 'N1H1B_O4'}.get(h['name'], 'O5H1A_N1')
        add(f'qtaim.{pre}.{nm}.DHA', h['DHA'], 'deg', 1, f'molecules/sfa/results/{c}/mwfn_extra/qtaim.json:hbonds.DHA_deg', f'{pre}: D–H···A angle at the {nm} BCP contact', case=c)
for key_, t in D6['tddft'].items():
    c, fn = key_.split('/')
    tag = fn.replace('tddft_', '')
    for st in t['states']:
        b = f"bulk.td.{c.split('_')[0]}.{tag}.S{st['state']:02d}"
        if st['CT'] is None and st['D_CT'] is not None:
            addb(b + '.DCT', st['D_CT'], 'A', 2, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}].ct', f'{c} {tag}: S{st["state"]} D_CT', case=c)
        for j_, (lab_, w_) in enumerate(st['dominant'][:3]):
            addb(f'{b}.w{j_}', 100 * w_, '%', 0, f'molecules/sfa/results/{key_}.json:states[{st["state"] - 1}].dominant[{j_}].weight', f'{c} {tag}: S{st["state"]} weight of {lab_}', case=c)

# ---- settings and literature values printed in the Introduction and Methods (added 28 Sep 2026 with the draft of
#      sections 1 and 2; every value is read from the file named in its source, or quoted from an abstract read
#      with `lit.py abstract`)
section('23. Settings and literature values quoted in the Introduction and Methods',
        'Protocol settings read from the scripts and result files that ran them, and the few literature values the '
        'Introduction quotes (source: the abstract of the cited paper, read with lit.py abstract on 28 Sep 2026).')
import ast as _ast  # noqa: E402
_src = lambda p: (REPO / p).read_text(encoding='utf-8')          # noqa: E731
_py_pcm = Path('/root/cqc-sfa-tools/env/lib/python3.1/site-packages/pyscf/solvent/pcm.py')
_pcm_txt = _py_pcm.read_text() if _py_pcm.exists() else ''
_m_h = re.search(r'modified_Bondi\[1\]\s*=\s*([\d.]+)/radii\.BOHR', _pcm_txt)
_m_s = re.search(r'self\.vdw_scale\s*=\s*([\d.]+)', _pcm_txt)
add('protocol.pcm_rH', float(_m_h.group(1)) if _m_h else 1.1, 'A', 2, 'pyscf 2.14.0 solvent/pcm.py: modified_Bondi[1] = 1.1/radii.BOHR (dft/common.py sets no radii_table)',
    'hydrogen radius of the default PySCF PCM cavity ("modified Bondi")')
add('protocol.pcm_vdw_scale', float(_m_s.group(1)) if _m_s else 1.2, None, 1, 'pyscf 2.14.0 solvent/pcm.py: self.vdw_scale = 1.2 (dft/common.py sets no vdw_scale)',
    'scale factor applied to the cavity radii')
_opt = J(RES / 'sfa_pair_pcm_d3/opt_summary.json')
add('protocol.accept_dE', 1e-6, 'Ha', 6, 'dft/proton_tools.py: ACCEPTANCE comment ("forces + energy" rule, convergence_energy 1e-6 Ha)',
    'energy-change limit of the "forces + energy" acceptance rule (scan and dielectric ladder)')
add('protocol.accept_n', 3, None, 0, 'dft/proton_tools.py: NEGLIGIBLE_STEPS = 3', 'consecutive gradients required by either extra acceptance rule')
add('protocol.accept_factor', 10, None, 0, 'dft/proton_tools.py: "negligible forces" rule, 10x inside the rms/max force thresholds',
    'factor inside the force thresholds for the "negligible forces" rule')
add('pcm.opt.crawl_from', 51, None, 0, 'dft/proton_tools.py comment on burst 1a sfa_pair_pcm_d3 (logs/optimize.log): steps 51-64 crawl',
    'first step of the slow inter-ion crawl of the IEF-PCM ion-pair optimisation', case='sfa_pair_pcm_d3')
add('pcm.opt.crawl_dE', 1.6e-5, 'Ha', 6, 'dft/proton_tools.py comment: "The whole 51->64 crawl gained 1.6e-5 Ha (0.01 kcal/mol)"',
    'energy lowered from step 51 to step 64 of the IEF-PCM ion-pair optimisation', case='sfa_pair_pcm_d3')
add('pcm.opt.crawl_dE_kcal', 1.6e-5 * 627.5095, 'kcal mol-1', 2, 'pcm.opt.crawl_dE converted', 'the same in kcal mol⁻¹', case='sfa_pair_pcm_d3')
_scan_doc = _src('dft/scan.py')
_m_g = re.search(r'--grid LO:HI:STEP \(default ([\d.]+):([\d.]+):([\d.]+)\)', _scan_doc)
add('protocol.scan_step', float(_m_g.group(3)) if _m_g else 0.10, 'A', 2, 'dft/scan.py: --grid default 1.05:1.75:0.10 (one absolute grid for both chains)',
    'grid step of the relaxed N1–H1A scan')
add('protocol.scan_B_start_OH', 1.00, 'A', 2, 'dft/scan.py docstring: chain B anchor starts with H1A on the O5->N1 line 1.00 A from O5',
    'O5–H1A distance of the O-side start of scan chain B')
_qe_doc = _src('dft/periodic_qe.py')
_m_oh = re.search(r'^D_OH_START\s*=\s*([\d.]+)', _qe_doc, re.M)
add('protocol.qe_start_OH', float(_m_oh.group(1)) if _m_oh else 0.98, 'A', 2, 'dft/periodic_qe.py: D_OH_START (start ii; also the iirfix constrained neutral form)',
    'O5–H1A distance at which H1A was placed in the periodic neutral starts')
add('protocol.qe_forc_conv', 5e-4, 'Ry/bohr', 4, 'periodic_qe.json:method.relax (BFGS, forc_conv_thr 5e-4 Ry/bohr)', 'BFGS force threshold of the periodic relaxations')
_ct = PQ['convergence']['targets']
add('periodic.conv_dE_meV', _ct['dE_meV_atom'], 'meV', 1, 'periodic_qe.json:convergence.targets.dE_meV_atom', 'energy tolerance per atom of the cutoff and k-point choice')
add('periodic.conv_dF', _ct['dF_Ry_bohr'], 'Ry/bohr', 3, 'periodic_qe.json:convergence.targets.dF_Ry_bohr', 'force tolerance of the cutoff and k-point choice')
add('protocol.mep_iso', 0.001, 'a.u.', 3, 'dft/mwfn.py: MEP on the 0.001 a.u. electron-density surface (Multiwfn main function 12)',
    'electron-density isovalue of the MEP surface')
add('hs.cluster_cutoff', D7['cluster_cutoff_A'], 'A', 1, 'derived/d07_hirshfeld.json:cluster_cutoff_A (dft/hirshfeld.py --cutoff)',
    'Hirshfeld cluster: every ion with an atom within this distance of the central ion')
_m_b = re.search(r'^BONDI\s*=\s*(\{[^}]*\})', _src('dft/hirshfeld.py'), re.M)
_bondi = {k_: v_ for k_, v_ in _ast.literal_eval(_m_b.group(1)).items() if k_ in ('H', 'C', 'N', 'O', 'P', 'S')}
add('hs.bondi', _bondi, 'A', 2, 'dft/hirshfeld.py: BONDI (radii of the CrystalExplorer d_norm)', 'Bondi van der Waals radii used in d_norm')
_inv = J(REPO / 'dock/targets_inventory.json')
for t_ in ('T1_3TZF', 'T2_1AJ0', 'T3_6RL9', 'T4_3HS4', 'T5_7S2K'):
    add(f"dock.res.{_inv[t_]['pdb']}", _inv[t_]['resolution'], 'A', 2, f'dock/targets_inventory.json:{t_}.resolution (deposited mmCIF)',
        f"crystallographic resolution of {_inv[t_]['pdb']} ({_inv[t_]['title']})")
_an = _src('dock/analyse.py')
add('dock.gate.rmsd_strict', float(re.search(r'^RMSD_GATE\s*=\s*([\d.]+)', _an, re.M).group(1)), 'A', 1, 'dock/analyse.py: RMSD_GATE',
    'gate G2: mean top-pose RMSD limit (and success-fraction distance)')
add('dock.gate.rmsd_conv', float(re.search(r'^RMSD_LOOSE\s*=\s*([\d.]+)', _an, re.M).group(1)), 'A', 1, 'dock/analyse.py: RMSD_LOOSE',
    'conventional 2.0 Å criterion, evaluated beside gate G2 (also the clustering radius of the production poses, dock/RESULTS.md section 2)')
add('dock.gate.min_fraction', 0.5, None, 1, 'dock/analyse.py: passes_G2 (success fraction and site-check pass fractions >= 0.5)', 'minimum fraction of seeds for the success and site criteria')
add('dock.gate.pb_min', 0.8, None, 1, 'dock/analyse.py: passes_G2 (posebusters_physics_pass_fraction >= 0.8)', 'minimum PoseBusters physics pass fraction')
add('dock.gate.ZnN_max', 2.3, 'A', 1, "dock/analyse.py: SITE_CHECKS 'Zn-N(sulfonamide)' limit 2.3", 'hCA II site check: Zn–N(sulfonamide) limit')
add('dock.gate.Thr199_max', 3.3, 'A', 1, "dock/analyse.py: SITE_CHECKS 'Thr199 OG1 H-bond' limit 3.3", 'hCA II site check: Thr199 OG1···N(sulfonamide) limit')
add('dock.water_cutoff', float(re.search(r'^METAL_WATER_CUTOFF\s*=\s*([\d.]+)', _src('dock/prep_receptor.py'), re.M).group(1)), 'A', 1,
    'dock/prep_receptor.py: METAL_WATER_CUTOFF (prep_report.json water_rule)', 'crystal waters kept only within this distance of a kept metal')
for t_ in ('T1_3TZF', 'T5_7S2K'):
    add(f"dock.waters_kept.{t_.split('_')[1]}", len(J(REPO / f'dock/{t_}/prep_report.json')['waters_kept']), None, 0, f'dock/{t_}/prep_report.json:waters_kept',
        f"Mg²⁺-bound waters kept in {t_.split('_')[1]}")
add('dock.val.tz_offset_3HS4', 0.07, 'A', 2, 'dock/VALIDATION.md section 2 (receptor_tz.pdbqt vs the co-crystal ligand)', 'AD4Zn TZ pseudo-atom distance from the crystal sulfonamide N (3HS4)')
add('dock.val.n_ca', 257, None, 0, 'dock/VALIDATION.md section 3 (cross-docking references: 257 common Cα atoms)', 'Cα atoms used to superpose 6RL9 and 3HS4')
_lit = 'lit.py abstract (read 28 Sep 2026) of '
add('lit.dpka_cocrystal_below', -1, None, 0, _lit + 'das2024crystals (doi:10.1021/acs.jpcc.4c05301) and cruzcabeza2012acidbase',
    'ΔpKa rule: below this value a cocrystal is expected')
add('lit.dpka_salt_above', 4, None, 0, _lit + 'das2024crystals and cruzcabeza2012acidbase (linear relation between −1 and 4)', 'ΔpKa rule: above this value a salt is expected')
add('lit.dpka_n_complexes', 6465, None, 0, _lit + 'cruzcabeza2012acidbase (doi:10.1039/c2ce26055g)', 'CSD acid–base complexes behind the quantified ΔpKa rule')
add('lit.dpka_unreliable', [0, 3], None, 0, _lit + 'peluzo2025thiabendazole (doi:10.1016/j.cplett.2024.141772): "0 ≤ ΔpKa ≤ 3"', 'ΔpKa range in which the rule is said to fail')
add('lit.abramov_NO_short', 2.57, 'A', 2, _lit + 'abramov2024it (doi:10.1021/acs.cgd.4c00195): short N···O "below ∼2.57 Å at 100–120 K"',
    'N···O length below which efficient periodic DFT failed to describe continuum configurations')
add('saph.dpka_approx', 0, None, 0, 'saph.dpka = −0.13 (anilinium 2.02, greenbaum1955comparison via zheng2026digitized; H3PO4 2.148, bates1951first; author decision 29 Sep 2026)',
    'ΔpKa of sulfanilamide (anilinium) and H3PO4 (first dissociation), ≈ 0')

# =====================================================================================================================
# report-text consistency checks (numbers printed in the pipeline's draft RESULTS.md vs the files)
DISAGREE.append('molecules/sfa/results/sfa_pair_pcm_d3/RESULTS.md (draft by dft/report.py) says "integration grid level 3" and "Multiwfn ." with no version; the runs used grid level 4 (scan_summary.json level.grid_level, every provenance block, the opt.xyz comment "grid 4") and Multiwfn 2026.9.20. Use grid level 4.')
DISAGREE.append('RESULTS.md (draft) reports "ZPE 5.299 eV" — correct (0.19472 Ha = 122.19 kcal/mol) but not a quantity for the paper; its "IP/EA/χ/η/ω" table is Koopmans-type, while dscf.json gives the ΔSCF values (I 7.488 vs 7.642 eV; A 1.857 vs 1.666 eV; η 2.816 vs 2.988 eV). Use ΔSCF in the paper, Koopmans only as the comparison.')
DISAGREE.append('KEY STORY "barrier (~0.6–0.7 kcal/mol)": the B3LYP-D3 electronic barrier on the scan grid is 0.72 kcal/mol (0.75 by spline; 1-D DVR ground level 0.66 kcal/mol above the salt). Quote 0.72 (or "about 0.7"), not 0.6.')
DISAGREE.append('KEY STORY "lower with every other functional": true for the end points with all five methods (ωB97X-D −0.72, PBE0-D3 −1.27, DF-MP2 −1.42, M06-2X −2.47 kcal/mol; B3LYP-D3 −0.54), but PBE0-D3 and M06-2X also put the d = 1.25 Å point below the salt (no barrier on this path), and DF-MP2 puts it only +0.20 kcal/mol above; the salt minimum is a shallow feature of B3LYP-D3, ωB97X-D and DF-MP2 at these geometries, below the N–H zero-point energy in every case.')
DISAGREE.append('KEY STORY "dielectric ladder ε 1 O-side, 8/20/78 N-side": correct, and ε = 4 is also O-side (δ +0.675 Å); a salt minimum first appears between ε = 4 and 8 (one start per ε: this brackets where the minimum appears, not a thermodynamic crossover).')
DISAGREE.append('RESULTS.md (draft) calls ADCH/CM5/Hirshfeld "fragment charge transfer" only; ADCH atomic charges on C3/C6/H3/H6 are sign-inverted (unphysical) in the PCM pair.')
DISAGREE.append('3HS4 cross-dock from the DFT geometry: the RMSD to 6RL9 sulfanilamide carried into 3HS4 is 3.2150 Å at full precision (derive/d09_dock.py, 20 poses; per-seed values '
                'rounded to 3 decimals also average exactly 3.2150) and prints as 3.22 ± 0.02 Å; production_summary.json stores 3.215, which a binary float '
                'rounds down to 3.21. Use 3.22 everywhere (facts dock.prod.3HS4.sfa.ad4zn.rmsd, Fig. 9b, captions, Table 7). RESOLVED: dock/RESULTS.md and dock/VALIDATION.md '
                'print 3.22 ± 0.02 Å for this RMSD (checked in the numbers audit, round 2, 29 Sep 2026).')

# =====================================================================================================================
# SI completion (29 Sep 2026): Fukui functions, ADMET, figure parameters, Fig. 10 representative poses (sections 24-27).
exec(compile((HERE / 'facts_si_extra.py').read_text(encoding='utf-8'), str(HERE / 'facts_si_extra.py'), 'exec'))

# =====================================================================================================================
# Burst 2 (folded in 29 Sep 2026): periodic constrained neutral + PBE0 correction, DFT strain of the docked poses, the
# withdrawn iifix facts and the rounding rule (sections 28-30; the functional check with DF-MP2 is section 7).
exec(compile((HERE / 'facts_burst2.py').read_text(encoding='utf-8'), str(HERE / 'facts_burst2.py'), 'exec'))

# =====================================================================================================================
# Revision after the referee reports (29 Sep 2026): section 31 (derive/d12_revision.py; Pi reference single points).
exec(compile((HERE / 'facts_revision.py').read_text(encoding='utf-8'), str(HERE / 'facts_revision.py'), 'exec'))

# =====================================================================================================================
# Referee calculations (VM burst 3, 29 Sep 2026): sections 32-36 (derive/d13_referee_calcs.py): periodic N1-H1A profile,
# neutral-minimum Hessian, nine-bond cluster QTAIM, relaxed torsion scans; superseded facts relabelled or retired.
exec(compile((HERE / 'facts_referee.py').read_text(encoding='utf-8'), str(HERE / 'facts_referee.py'), 'exec'))

# =====================================================================================================================
# Closing items (29 Sep 2026, night): section 37 (derive/d14_completion.py): the PLBD constants of the hCA II controls, one
# definition of "near-native" for the AutoDock4Zn orientation gaps, and the raw-pose torsions of Table S31.
exec(compile((HERE / 'facts_completion.py').read_text(encoding='utf-8'), str(HERE / 'facts_completion.py'), 'exec'))

# =====================================================================================================================
# Control receptor (4 Oct 2026): section 38 (derive/d15_control.py): the ATP site of E. coli DNA gyrase B (4DUH), which
# is not a target of sulfanilamide, docked with the protocol of the five targets (dock/control_gyrb.py).
exec(compile((HERE / 'facts_control.py').read_text(encoding='utf-8'), str(HERE / 'facts_control.py'), 'exec'))
# Review before submission to the Journal of the Indian Chemical Society (5 Oct 2026): section 39 (derive/d16_jics.py):
# the 1-D proton model rebuilt with its sensitivity, zero-point energy of the separated species, def2-TZVP shift of
# the proton, wavelength of the UV determination of aqueous sulfanilamide.
exec(compile((HERE / 'facts_jics.py').read_text(encoding='utf-8'), str(HERE / 'facts_jics.py'), 'exec'))
# Molecular dynamics of the hCA II-sulfanilamide complex and the two docking tests, Glu106 and four crystal waters
# (8 Oct 2026): section 40 (derive/d17_md.py; md/RESULTS.md, dock/GLU106.md, dock/WATER4.md).
exec(compile((HERE / 'facts_md.py').read_text(encoding='utf-8'), str(HERE / 'facts_md.py'), 'exec'))
# Final wording of the notes, titles and statuses (ACS integration review, 30 Sep 2026): no value changes.
exec(compile((HERE / 'facts_final_notes.py').read_text(encoding='utf-8'), str(HERE / 'facts_final_notes.py'), 'exec'))

# =====================================================================================================================
meta = dict(generated_by='paper/main/data/build_facts.py', n_bulk=sum(1 for v in FACTS.values() if v.get('bulk')), inputs_derived=sorted(p.name for p in DER.glob('*.json')),
            n_facts=len(FACTS), n_rederived=sum(1 for v in FACTS.values() if 'rederived' in v),
            n_agree=sum(1 for v in FACTS.values() if v.get('agrees') is True),
            n_disagree=sum(1 for v in FACTS.values() if v.get('agrees') is False),
            disagreements=DISAGREE, caveats=CAVEATS, missing=MISSING,
            sections=[dict(title=t, keys=k) for t, _, k in SECTIONS],
            conventions='units: A = ångström, deg = degree; energies as labelled; precision = decimals to print; '
                        'withdrawn = true: value is null and the number is kept in withdrawn_value for provenance only (not cited); '
                        'superseded = true: value is null and the number is kept in superseded_value (replaced by a later calculation; not cited)',
            rounding_rule=ROUNDING_RULE)
out = {'_meta': meta}
out.update(FACTS)
(HERE / 'facts.json').write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))


# ---------------------------------------------------------------------------------------------------------------- md
def fmt(v, p):
    if isinstance(v, bool) or v is None:
        return str(v)
    if isinstance(v, (int, float)) and isinstance(p, int):
        if isinstance(v, float) and (abs(v) < 1e-3 and v != 0) and p > 6:
            return f'{v:.3e}'
        return f'{v:.{p}f}'.replace('-', '−') if isinstance(v, float) or p > 0 else str(v)
    if isinstance(v, float):
        return f'{v:.3g}'
    if isinstance(v, list):
        return ', '.join(fmt(x, p if isinstance(p, int) else 3) for x in v)
    return str(v)


L = ['# Facts dossier: 4-sulfamoylanilinium dihydrogen phosphate (SAPH)', '',
     'Every number the main text and the SI may print, with its unit, the precision to print, its source (file + JSON key or log line) '
     'and what it means. Same content as `facts.json` (flat key → {value, unit, precision, case, source, meaning, note, status, '
     'rederived, agrees}). Generated by `paper/main/data/build_facts.py` from the result files and from the independent '
     're-derivations in `derive/d01–d15` (outputs in `derived/`). Regenerate:', '',
     '    P=/root/workspace/cqc-dft/tools/paper-env/bin/python; cd paper/main/data/derive',
     '    for s in d01_xray d03_geometry d04_scan d05_freq d06_electronic d07_hirshfeld d08_periodic d09_dock d10_rep_contacts d11_burst2 d12_revision d13_referee_calcs d14_completion d15_control; do nice $P $s.py; done',
     '    nice /root/cqc-sfa-tools/env/bin/python d02_models.py   # PySCF basis counts',
     '    cd ../../../.. && nice $P paper/main/data/build_facts.py && nice $P paper/main/data/build_tables.py', '',
     f"**{meta['n_facts']} facts ({meta['n_facts'] - meta['n_bulk']} itemised below, {meta['n_bulk']} bulk SI-table values in section 21); {meta['n_rederived']} independently re-derived, {meta['n_agree']} agree, {meta['n_disagree']} disagree.** "
     'Column "check": ✓ = re-derived by our own script and equal within the printed precision; ✗ = disagreement (listed below); '
     'blank = read from the source only. Model names: ion pair (27 atoms, anion at x, y, z − 1), gas phase, IEF-PCM (water), '
     'IEF-PCM (ε = n), periodic model. Level: B3LYP-D3(BJ)/6-311++G(d,p) unless stated.', '',
     f"**Rounding rule**: {ROUNDING_RULE}", '',
     '**Withdrawn facts**: the relax_iifix energies and structure metrics (keys containing "iifix", except the three slide '
     'diagnostics kept for a methodological note) are kept for provenance with a null value and the number struck through '
     'below; they are not cited (section 30).', '',
     '**Superseded facts** (section 36): the rigid DFT strain of the docked poses and the order-of-magnitude checks on the '
     'constrained configuration are retired (value null, the number struck through below); the relax_iirfix energies stay citable '
     'only as the constrained single configuration (O5–H1A held at 0.98 Å), not as the periodic cost.', '',
     '## Notes on the values', '']
L += [f'- {d}' for d in DISAGREE] + ['', '## Caveats of the values', ''] + [f'- {c}' for c in CAVEATS]
L += ['', '## Not computed (stated limitations)', ''] + [f'- {m}' for m in MISSING] + ['']
for title, intro, keys in SECTIONS:
    L += ['---', '', f'## {title}', '']
    if intro:
        L += [intro, '']
    if title.startswith('21.'):
        L += ['| key prefix | facts | source | what |', '|---|---|---|---|']
        for pre, (n, src, what) in BULK_PREFIX.items():
            L.append(f'| `{pre}.*` | {n} | `{src}` | {what} |')
        L.append('')
        continue
    L += ['| key | value | unit | meaning | source | check |', '|---|---|---|---|---|---|']
    for k in keys:
        f = FACTS[k]
        if f.get('bulk'):
            continue
        chk = '' if 'agrees' not in f else ('✓' if f['agrees'] else f"✗ ({fmt(f['rederived'], f['precision'] if isinstance(f['precision'], int) else 3)})")
        mean = f['meaning'] + (f" — *{f['status']}*" if f.get('status') else '') + (f" ({f['note']})" if f.get('note') else '')
        val = (f"~~{fmt(f['withdrawn_value'], f['precision'])}~~ (withdrawn)" if f.get('withdrawn') else
               f"~~{fmt(f['superseded_value'], f['precision'])}~~ (superseded)" if f.get('superseded') else fmt(f['value'], f['precision']))
        L.append(f"| `{k}` | {val} | {f['unit'] or ''} | {mean.replace('|', '/')} | `{str(f['source']).replace('|', '/')}` | {chk} |")
    L.append('')
(HERE / 'facts.md').write_text('\n'.join(L))
print(f"facts: {meta['n_facts']}  rederived: {meta['n_rederived']}  agree: {meta['n_agree']}  disagree: {meta['n_disagree']}")
for d in DISAGREE:
    if 'source' in d and 're-derived' in d:
        print('DISAGREE', d)
