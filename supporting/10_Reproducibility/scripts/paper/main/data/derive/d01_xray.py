"""d01: crystal facts from the deposited CIF, recomputed with our own fractional->Cartesian code.

Reads e-70-0o997-sup1.cif (gemmi.cif only as a tokenizer), rebuilds every hydrogen bond of the CIF table from
the _atom_site coordinates and the SHELX symmetry codes, compares with the deposited values, and measures the
ion-pair (anion at x, y, z-1) geometry that the DFT models start from.
"""
import re
import gemmi
import numpy as np
from common import *

blk = gemmi.cif.read(str(CIF)).sole_block()


def val(tag):
    v = blk.find_value(tag)
    return None if v is None else gemmi.cif.as_string(v)


def num_su(s):
    """'4.8041(7)' -> (4.8041, 0.0007)"""
    m = re.match(r'^(-?[\d.]+)(?:\((\d+)\))?$', s.strip())
    if not m:
        return None, None
    x = float(m.group(1))
    if m.group(2) is None:
        return x, None
    dec = len(m.group(1).split('.')[1]) if '.' in m.group(1) else 0
    return x, int(m.group(2)) * 10 ** (-dec)


a, _ = num_su(val('_cell_length_a')); b, _ = num_su(val('_cell_length_b')); c, _ = num_su(val('_cell_length_c'))
al, _ = num_su(val('_cell_angle_alpha')); be, _ = num_su(val('_cell_angle_beta')); ga, _ = num_su(val('_cell_angle_gamma'))
alr, ber, gar = np.radians([al, be, ga])
# PDB/Cambridge convention: a along x, b in the xy plane
cx = c * np.cos(ber)
cy = c * (np.cos(alr) - np.cos(ber) * np.cos(gar)) / np.sin(gar)
cz = np.sqrt(c * c - cx * cx - cy * cy)
M = np.array([[a, b * np.cos(gar), cx], [0, b * np.sin(gar), cy], [0, 0, cz]])   # columns = a, b, c vectors
vol = float(np.linalg.det(M))

sites = {}
tab = blk.find('_atom_site_', ['label', 'type_symbol', 'fract_x', 'fract_y', 'fract_z', 'U_iso_or_equiv', 'calc_flag'])
for r in tab:
    lab = r[0]
    fx = [num_su(r[i])[0] for i in (2, 3, 4)]
    sites[lab] = dict(el=r[1], f=np.array(fx), uiso=num_su(r[5])[0], calc=r[6],
                      su=[num_su(r[i])[1] for i in (2, 3, 4)])

SYMOPS = {1: (np.eye(3), np.zeros(3)), 2: (np.diag([1, -1, 1]), np.array([0, 0, 0.5]))}   # Pc: x,y,z ; x,-y,z+1/2


def cart_sym(lab, code='1_555'):
    n, klm = code.split('_')
    R, t = SYMOPS[int(n)]
    tr = np.array([int(klm[0]) - 5, int(klm[1]) - 5, int(klm[2]) - 5])
    f = R @ sites[lab]['f'] + t + tr
    return M @ f


out = {'source': rel(CIF)}
cd = {}
for tag in ['_database_code_depnum_ccdc_archive', '_journal_paper_doi', '_journal_year', '_journal_volume',
            '_journal_page_first', '_chemical_formula_moiety', '_chemical_formula_sum', '_chemical_formula_weight',
            '_space_group_crystal_system', '_space_group_name_H-M_alt', '_space_group_name_Hall',
            '_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_beta', '_cell_volume',
            '_cell_formula_units_Z', '_cell_measurement_temperature', '_diffrn_ambient_temperature',
            '_diffrn_radiation_type', '_diffrn_radiation_wavelength', '_diffrn_measurement_device_type',
            '_exptl_crystal_description', '_exptl_crystal_colour', '_exptl_crystal_size_max',
            '_exptl_crystal_size_mid', '_exptl_crystal_size_min', '_exptl_crystal_density_diffrn',
            '_exptl_crystal_F_000', '_exptl_absorpt_coefficient_mu', '_exptl_absorpt_correction_type',
            '_diffrn_reflns_number', '_diffrn_reflns_av_R_equivalents', '_diffrn_reflns_theta_min',
            '_diffrn_reflns_theta_max', '_reflns_number_total', '_reflns_number_gt', '_refine_ls_R_factor_gt',
            '_refine_ls_R_factor_all', '_refine_ls_wR_factor_ref', '_refine_ls_wR_factor_gt',
            '_refine_ls_goodness_of_fit_ref', '_refine_ls_number_parameters', '_refine_ls_number_restraints',
            '_refine_ls_hydrogen_treatment', '_refine_diff_density_max', '_refine_diff_density_min',
            '_refine_ls_abs_structure_Flack', '_computing_structure_refinement']:
    cd[tag] = val(tag)
out['cif_items'] = cd
out['cell_volume_recomputed_A3'] = vol

# hydrogen treatment, from the embedded .res and the atom_site flags
res_txt = CIF.read_text()
out['dfix'] = sorted(set(re.findall(r'^DFIX .*$', res_txt, flags=re.M)))
out['h_sites'] = {l: {'refined_xyz': s['su'][0] is not None, 'calc_flag': s['calc']}
                  for l, s in sites.items() if s['el'] == 'H'}

# CIF hydrogen-bond table vs recomputation
hb = []
t = blk.find('_geom_hbond_', ['atom_site_label_D', 'atom_site_label_H', 'atom_site_label_A', 'site_symmetry_A',
                              'distance_DH', 'distance_HA', 'distance_DA', 'angle_DHA'])
for r in t:
    D, H, A, sym = r[0], r[1], r[2], r[3]
    dD, dH, dA = cart_sym(D), cart_sym(H), cart_sym(A, sym)
    rec = dict(D=D, H=H, A=A, sym=sym,
               cif=dict(DH=r[4], HA=r[5], DA=r[6], DHA=r[7]),
               calc=dict(DH=dist(dD, dH), HA=dist(dH, dA), DA=dist(dD, dA), DHA=angle(dD, dH, dA)))
    rec['dDA_calc_minus_cif'] = rec['calc']['DA'] - num_su(r[6])[0]
    hb.append(rec)
# the same with the full-precision coordinates of the embedded SHELX .res (the _atom_site H coordinates carry only
# 3 decimals, which moves H...A by up to ~0.007 A)
res_part = res_txt.split('_shelx_res_file')[1]
resxyz = {}
for lab in sites:
    m_ = re.search(r'^' + re.escape(lab) + r'\s+\d+\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)', res_part, re.M)
    if m_:
        resxyz[lab] = np.array([float(x) for x in m_.groups()])


def cart_res(lab, code='1_555'):
    n, klm = code.split('_')
    R, t = SYMOPS[int(n)]
    tr = np.array([int(klm[0]) - 5, int(klm[1]) - 5, int(klm[2]) - 5])
    return M @ (R @ resxyz[lab] + t + tr)


for rec in hb:
    D, H, A, sym = rec['D'], rec['H'], rec['A'], rec['sym']
    dD, dH, dA = cart_res(D), cart_res(H), cart_res(A, sym)
    rec['calc_res'] = dict(DH=dist(dD, dH), HA=dist(dH, dA), DA=dist(dD, dA), DHA=angle(dD, dH, dA))
out['hbonds'] = hb
# bonds of the ion pair from the full-precision .res coordinates (the _atom_site coordinates are rounded to their s.u.,
# which moves a bond by up to ~2e-5 A; enough to decide a rounding tie such as P1-O6 = 1.4975 in the CIF; 29 Sep 2026)
out['pair_res'] = {b: dist(cart_res(b.split('-')[0]), cart_res(b.split('-')[1]))
                   for b in ['C1-N1', 'P1-O3', 'P1-O4', 'P1-O5', 'P1-O6']}
out['n_hbonds'] = len(hb)
out['max_abs_dDA'] = max(abs(h['dDA_calc_minus_cif']) for h in hb)

# selected deposited bonds / angles / torsions (with s.u.)
geo = {}
t = blk.find('_geom_bond_', ['atom_site_label_1', 'atom_site_label_2', 'distance', 'site_symmetry_2'])
for r in t:
    geo[f'{r[0]}-{r[1]}'] = r[2]
out['cif_bonds'] = geo
ang = {}
t = blk.find('_geom_angle', ['_atom_site_label_1', '_atom_site_label_2', '_atom_site_label_3', ''])
for r in t:
    ang[f'{r[0]}-{r[1]}-{r[2]}'] = r[3]
out['cif_angles'] = ang
tor = {}
t = blk.find('_geom_torsion', ['_atom_site_label_1', '_atom_site_label_2', '_atom_site_label_3', '_atom_site_label_4', ''])
for r in t:
    tor[f'{r[0]}-{r[1]}-{r[2]}-{r[3]}'] = r[4]
out['cif_torsions'] = tor

# the ion pair (anion at x, y, z-1) rebuilt here and compared with structures/sfa_pair_xray.xyz
pair = np.array([cart_sym(l) for l in CATION] + [cart_sym(l, '1_554') for l in ANION])
syms, ref, _ = read_xyz(REPO / 'structures/sfa_pair_xray.xyz')
# same frame? compare internal distances (frame-independent)
dm1 = np.linalg.norm(pair[:, None] - pair[None], axis=2)
dm2 = np.linalg.norm(ref[:, None] - ref[None], axis=2)
out['pair_model_max_abs_distance_diff_A'] = float(np.abs(dm1 - dm2).max())
P = {l: pair[i] for i, l in enumerate(PAIR_LABELS)}
g = {}
g['d_N1_H1A'] = dist(P['N1'], P['H1A']); g['d_H1A_O5'] = dist(P['H1A'], P['O5'])
g['delta'] = g['d_N1_H1A'] - g['d_H1A_O5']; g['d_N1_O5'] = dist(P['N1'], P['O5'])
g['ang_N1_H1A_O5'] = angle(P['N1'], P['H1A'], P['O5'])
for x in ['O3', 'O4', 'O5', 'O6']:
    g[f'P1-{x}'] = dist(P['P1'], P[x])
g['C1-N1'] = dist(P['C1'], P['N1'])
g['PO_long_mean'] = (g['P1-O3'] + g['P1-O4']) / 2
g['PO_short_mean'] = (g['P1-O5'] + g['P1-O6']) / 2
g['PO_split'] = g['PO_long_mean'] - g['PO_short_mean']
g['C6_O4'] = dist(P['C6'], P['O4']); g['H6_O4'] = dist(P['H6'], P['O4'])
g['ang_C6_H6_O4'] = angle(P['C6'], P['H6'], P['O4'])
g['tau_C3_C4_S1_N2'] = dihedral(P['C3'], P['C4'], P['S1'], P['N2'])
g['tau_C5_C4_S1_N2'] = dihedral(P['C5'], P['C4'], P['S1'], P['N2'])
g['ang_C4_S1_N2'] = angle(P['C4'], P['S1'], P['N2'])
g['ang_O1_S1_O2'] = angle(P['O1'], P['S1'], P['O2'])
out['pair'] = g
# the O4-H4A...O5 anion-anion bond (O5 at 1_455 of the anion's own frame)
o4, h4, o5 = cart_sym('O4'), cart_sym('H4A'), cart_sym('O5', '1_455')
out['O4_H4A_O5'] = dict(OH=dist(o4, h4), HO=dist(h4, o5), OO=dist(o4, o5), ang=angle(o4, h4, o5),
                        frac_along=float((h4 - o4) @ (o5 - o4) / dist(o4, o5) ** 2))
# neutron-normalised O-H (0.983 A) position of H4A: where would the proton sit on the O4...O5 line?
out['xh_note'] = 'O-H restrained to 0.82 A (DFIX); C-H riding (AFIX 43, 0.93 A); N-H refined freely'
dump('d01_xray.json', out)
