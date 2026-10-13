"""d11: burst-2 results re-derived from the raw outputs (29 Sep 2026).

(1) Periodic constrained neutral (relax_iirfix) and the PBE0-D3(BJ) hybrid correction.  Energies are parsed from the
    archived Quantum ESPRESSO outputs, not from periodic_qe.json: the 'Final energy' of each relaxation (SSSP, 70/560
    Ry), the '!' energy of each PBE-D3(BJ) single point with SG15 potentials, and the first '!' (PBE electrons, pbe0
    D3(BJ) parameters) and the '!!' (EXX self-consistent) energies of each PBE0-D3(BJ) run.  The two atoms held in
    space (if_pos 0 0 0 in the relax input) are identified from the input and located in the final coordinates of the
    output; their final positions are compared with the start, and the final forces on held and free atoms are read
    from the last force block.  The single-point inputs are checked to use the relaxed coordinates.
(2) Functional check (fsp_scan.json, five methods incl. DF-MP2): relative energies recomputed from each E_Ha with our
    own conversion; the three geometry files re-hashed (sha256) and re-measured; the four DFT energies compared with
    the burst-1 log parse of d04 (an independent earlier run on the same geometries); DF-MP2 E = E_RHF + E_corr and the
    frozen-core count checked.
(3) DFT strain of the docked poses (dock_strain.json of the sulfanilamide, sulfonamidate and cation cases): E(pose) -
    E(opt) recomputed; every geometry file re-hashed; tau_fold re-measured from each pose .xyz; the heavy atoms of the
    as-docked and the polar-H re-optimised (Hrot) geometries compared; E(opt) compared with the case optimisation.
Reads only; writes derived/d11_burst2.json.
"""
import hashlib
import re

from common import *

PQA = ARES / 'periodic_qe' / 'runs'
PQ = RES / 'periodic_qe'
Z = 2
RY_MEV = RY_EV * 1000.0
RY_KJ_FU = RY_EV * EV_KJ / Z


def energies_qe(path):
    """All '!' and '!!' total energies (Ry) and the BFGS 'Final energy' of a pw.x output."""
    txt = Path(path).read_text(errors='replace')
    bang = [float(m.group(1)) for m in re.finditer(r'^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry', txt, re.M)]
    bang2 = [float(m.group(1)) for m in re.finditer(r'^!!\s+total energy\s+=\s+(-?\d+\.\d+) Ry', txt, re.M)]
    fin = re.search(r'Final energy\s+=\s+(-?\d+\.\d+) Ry', txt)
    disp = [float(m.group(1)) for m in re.finditer(r'DFT-D3 Dispersion\s+=\s+(-?\d+\.\d+) Ry', txt)]
    exx_ok = 'EXX self-consistency reached' in txt
    return dict(bang=bang, bang2=bang2, final=float(fin.group(1)) if fin else None, disp=disp, exx_converged=exx_ok,
                job_done='JOB DONE' in txt)


def run(name):
    return PQA / name / f'{name}.out'


out = {}
# ---------------------------------------------------------------------------------------------------- (1) periodic
relax = {k: energies_qe(run(f'relax_{k}')) for k in ['i', 'ii', 'iifix', 'iirfix']}
E_relax = {k: r['final'] for k, r in relax.items()}
per = dict(E_relax_Ry=E_relax, relax_job_done={k: r['job_done'] for k, r in relax.items()})
per['dE_sssp_meV_cell'] = {k: (E_relax[k] - E_relax['i']) * RY_MEV for k in E_relax}
per['dE_sssp_kJ_fu'] = {k: (E_relax[k] - E_relax['i']) * RY_KJ_FU for k in E_relax}
per['dE_sssp_kcal_fu'] = {k: v / 4.184 for k, v in per['dE_sssp_kJ_fu'].items()}
hyb = {}
for var in ['sg15_e100', 'sg15_e80', 'sg15_e100_nq211']:
    for s in ['i', 'ii', 'iifix', 'iirfix']:
        hp, pp = run(f'hybrid_{s}_{var}'), run(f'pbe_{s}_{var}')
        if not hp.exists():
            continue
        h, p = energies_qe(hp), energies_qe(pp)
        hyb.setdefault(var, {})[s] = dict(E_pbe0_Ry=h['bang2'][-1], E_first_scf_Ry=h['bang'][0], exx_converged=h['exx_converged'],
                                          n_outer=len(h['bang']), E_pbe_Ry=p['bang'][-1], disp_pbe_Ry=p['disp'][-1] if p['disp'] else None)
comp = {}
for var, d in hyb.items():
    for s in ['ii', 'iifix', 'iirfix']:
        if s not in d or 'i' not in d:
            continue
        a, b = d[s], d['i']
        dh = (a['E_pbe0_Ry'] - b['E_pbe0_Ry']) * RY_MEV
        dp = (a['E_pbe_Ry'] - b['E_pbe_Ry']) * RY_MEV
        dfirst = (a['E_first_scf_Ry'] - b['E_first_scf_Ry']) * RY_MEV
        corr = dh - dp
        d3 = dfirst - dp
        comp.setdefault(var, {})[f'{s}_minus_i'] = dict(
            PBE0_D3_meV_cell=dh, PBE_D3_same_setting_meV_cell=dp, correction_meV_cell=corr,
            correction_kJ_fu=corr / 1000 * EV_KJ / Z, split_D3_params_meV_cell=d3, split_exact_exchange_meV_cell=corr - d3,
            PBE0_D3_kJ_fu=dh / 1000 * EV_KJ / Z, PBE_SG15_minus_SSSP_meV_cell=dp - per['dE_sssp_meV_cell'][s],
            best_estimate_meV_cell=per['dE_sssp_meV_cell'][s] + corr,
            best_estimate_kJ_fu=per['dE_sssp_kJ_fu'][s] + corr / 1000 * EV_KJ / Z,
            exx_converged=bool(a['exx_converged'] and b['exx_converged']))
per['hybrid_raw'] = hyb
per['hybrid'] = comp
cr = [comp[v]['iirfix_minus_i']['correction_kJ_fu'] for v in comp if 'iirfix_minus_i' in comp[v]]
per['iirfix_correction_spread_kJ_fu'] = max(cr) - min(cr)
per['iirfix_best_estimate_spread_kJ_fu'] = (max(comp[v]['iirfix_minus_i']['best_estimate_kJ_fu'] for v in comp)
                                           - min(comp[v]['iirfix_minus_i']['best_estimate_kJ_fu'] for v in comp))

# held atoms, final coordinates and final forces of relax_iirfix
rin = (PQA / 'relax_iirfix' / 'relax_iirfix.in').read_text()
rout = run('relax_iirfix').read_text(errors='replace')
blk = rin.split('ATOMIC_POSITIONS')[1].splitlines()[1:]
atoms_in = []
for l in blk:
    p = l.split()
    if len(p) < 4 or not re.match(r'^[A-Z][a-z]?$', p[0]):
        break
    atoms_in.append(dict(el=p[0], frac=np.array([float(x) for x in p[1:4]]), fixed=(len(p) >= 7 and p[4:7] == ['0', '0', '0'])))
fb = rout.split('Begin final coordinates')[1].split('End final coordinates')[0].splitlines()
fin = []
for l in fb:
    p = l.split()
    if len(p) >= 4 and re.match(r'^[A-Z][a-z]?$', p[0]):
        fin.append(np.array([float(x) for x in p[1:4]]))
assert len(fin) == len(atoms_in) == 54
fblocks = rout.split('Forces acting on atoms')
last = fblocks[-1].split('Total force')[0]
forces = np.array([[float(x) for x in m.groups()] for m in
                   re.finditer(r'atom\s+\d+\s+type\s+\d+\s+force =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', last)])
assert forces.shape == (54, 3)
syms, Xr, lab_r, cell, _ = read_extxyz(PQ / 'relax_iirfix.extxyz')
A = cell
fr_ext = np.linalg.solve(A.T, Xr.T).T                       # fractional coordinates of the extxyz atoms
lab_of = []
for q in fin:                                               # map each output atom onto the extxyz label (mod 1)
    dd = [np.linalg.norm(((q - f + 0.5) % 1.0 - 0.5) @ A) for f in fr_ext]
    j = int(np.argmin(dd))
    assert dd[j] < 1e-4, dd[j]
    lab_of.append(lab_r[j])
held = [lab_of[i] for i, a in enumerate(atoms_in) if a['fixed']]
held_move = max(float(np.linalg.norm(((fin[i] - a['frac'] + 0.5) % 1.0 - 0.5) @ A)) for i, a in enumerate(atoms_in) if a['fixed'])
fn = np.linalg.norm(forces, axis=1)
free = [i for i, a in enumerate(atoms_in) if not a['fixed']]
fixd = [i for i, a in enumerate(atoms_in) if a['fixed']]
P = dict(zip(lab_r, Xr))
imgs = [np.array(v) @ A for v in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)]
proj = {}
for suf in ['', '_2']:
    iH = lab_of.index('H1A' + suf)
    h, n = P['H1A' + suf], min((P['N1' + suf] + t for t in imgs), key=lambda q: np.linalg.norm(q - P['H1A' + suf]))
    u = (n - h) / np.linalg.norm(n - h)
    proj['fu1' if not suf else 'fu2'] = dict(F_H1A_along_H1A_to_N1_Ry_bohr=float(forces[iH] @ u),
                                             F_H1A_norm_Ry_bohr=float(fn[iH]))
per['iirfix_constraint'] = dict(held_atoms=held, held_atoms_max_move_A=held_move,
                                max_force_free_Ry_bohr=float(fn[free].max()), max_force_held_Ry_bohr=float(fn[fixd].max()),
                                max_force_component_free_Ry_bohr=float(np.abs(forces[free]).max()),
                                max_force_component_held_Ry_bohr=float(np.abs(forces[fixd]).max()),
                                bfgs_converged='bfgs converged' in rout, forc_conv_thr_Ry_bohr=float(re.search(r'forc_conv_thr\s*=\s*([\d.]+)d?(-?\d+)?', rin).group(1)) * 10 ** int(re.search(r'forc_conv_thr\s*=\s*[\d.]+d(-?\d+)', rin).group(1)),
                                F_H1A=proj, n_free=len(free), n_held=len(fixd))
# the hybrid single points were run at the relaxed coordinates
sp_dev = {}
for var in hyb:
    for s in hyb[var]:
        for kind in ('hybrid', 'pbe'):
            spin = (PQA / f'{kind}_{s}_{var}' / f'{kind}_{s}_{var}.in')
            if not spin.exists():
                continue
            t_ = spin.read_text().split('ATOMIC_POSITIONS')[1].splitlines()[1:]
            fr = [np.array([float(x) for x in l.split()[1:4]]) for l in t_ if len(l.split()) >= 4 and re.match(r'^[A-Z][a-z]?$', l.split()[0])][:54]
            _, Xs, ls, cs, _ = read_extxyz(PQ / f'relax_{s}.extxyz')
            fe = np.linalg.solve(cs.T, Xs.T).T
            dev = max(min(float(np.linalg.norm(((q - f + 0.5) % 1.0 - 0.5) @ cs)) for f in fe) for q in fr)
            sp_dev[f'{kind}_{s}_{var}'] = dev
per['single_point_geometry_max_dev_A'] = max(sp_dev.values())
per['n_single_points_checked'] = len(sp_dev)
out['periodic'] = per

# ---------------------------------------------------------------------------------------------------- (2) functional check
FSP = RES / 'sfa_pair_pcm_d3' / 'fsp_scan.json'
fs = json.loads(FSP.read_text())
d4 = json.loads((OUT / 'd04_scan.json').read_text())
LOGNAME = {'b3lyp-d3bj': 'B3LYP-D3(BJ)', 'wb97x-d': 'wB97X-D', 'm06-2x': 'M06-2X', 'pbe0-d3bj': 'PBE0-D3(BJ)'}
geoms = {}
for g in fs['geometries']:
    p = REPO / g['file']
    s, X, _ = read_xyz(p)
    Pp = {l: X[i] for i, l in enumerate(PAIR_LABELS)}
    geoms[g['name']] = dict(file=g['file'], sha256_ok=hashlib.sha256(p.read_bytes()).hexdigest() == g['sha256'],
                            d_N1_H1A=dist(Pp['N1'], Pp['H1A']), d_H1A_O5=dist(Pp['H1A'], Pp['O5']))
rel_ = {}
for m, e in fs['energies'].items():
    eN, eT, eO = e['end_N_d1.098']['E_Ha'], e['top1_d1.250']['E_Ha'], e['end_O_d1.631']['E_Ha']
    r_ = dict(top_N=(eT - eN) * HA_KCAL, dE_O_N=(eO - eN) * HA_KCAL, top_O=(eT - eO) * HA_KCAL,
              top_N_kJ=(eT - eN) * HA_KJ, dE_O_N_kJ=(eO - eN) * HA_KJ, top_is_max=bool(eT > max(eN, eO)),
              E_Ha={'end_N': eN, 'top1': eT, 'end_O': eO})
    if m in LOGNAME:
        lg = d4['functional_check']['raw'][LOGNAME[m]]
        r_['max_abs_E_vs_burst1_log_Ha'] = max(abs(e[g]['E_Ha'] - lg[g]['E_Ha']) for g in lg)
    if m == 'df-mp2':
        r_['max_abs_RHF_plus_corr_minus_E_Ha'] = max(abs(e[g]['E_RHF_Ha'] + e[g]['E_corr_Ha'] - e[g]['E_Ha']) for g in e)
        r_['n_frozen'] = sorted({e[g]['n_frozen'] for g in e})
        r_['aux'] = sorted({e[g]['aux_mp2'] for g in e})
    r_['json_relative_max_abs_diff_kcal'] = max(abs(fs['relative'][m]['vs_first_kcal'][g] - v) for g, v in
                                                [('top1_d1.250', r_['top_N']), ('end_O_d1.631', r_['dE_O_N'])])
    rel_[m] = r_
# frozen-core count from the element list of the pair: 1 core orbital per C/N/O, 5 per P/S
s_, _, _ = read_xyz(REPO / fs['geometries'][0]['file'])
ncore = sum(1 if el in ('C', 'N', 'O') else 5 if el in ('P', 'S') else 0 for el in s_)
out['fcheck'] = dict(source=rel(FSP), geometries=geoms, relative=rel_, frozen_core_expected=ncore,
                     n_methods=len(rel_), n_neutral_lower=sum(1 for r_ in rel_.values() if r_['dE_O_N'] < 0),
                     dE_O_N_range=[min(r_['dE_O_N'] for r_ in rel_.values()), max(r_['dE_O_N'] for r_ in rel_.values())],
                     top_N_range=[min(r_['top_N'] for r_ in rel_.values()), max(r_['top_N'] for r_ in rel_.values())],
                     notes=fs.get('notes'), provenance=fs['provenance'])

# ---------------------------------------------------------------------------------------------------- (3) docking strain
man = json.loads((DOCK / 'results/poses_for_dft/manifest.json').read_text())
MB = {p_['file'].replace('.xyz', ''): p_ for p_ in man['poses']}


def tau_fold(X):
    # atom order of the DFT opt.xyz: C1..C6 N1 N2 O1 O2 S1 H...
    c3, c4, c5, n2, s1 = X[2], X[3], X[4], X[7], X[10]
    return max(abs(dihedral(c3, c4, s1, n2)), abs(dihedral(c5, c4, s1, n2)))


strain = {}
for case, sp in [('sulfanilamide_pcm_d3', 'sfa_neutral'), ('sulfonamidate_pcm_d3', 'sfa_anion'), ('cation_pcm_d3', 'sfa_cation')]:
    js = json.loads((RES / case / 'dock_strain.json').read_text())
    e = js['energies']['b3lyp-d3bj']
    E0 = e['opt']['E_Ha']
    Eopt_case = json.loads((RES / case / 'opt_summary.json').read_text())['E_final_Ha']
    sha = {g['name']: g for g in js['geometries']}
    s0, X0, _ = read_xyz(REPO / sha['opt']['file'])
    poses = {}
    for name, g in sha.items():
        if name == 'opt':
            continue
        p = REPO / g['file']
        s, X, _ = read_xyz(p)
        base = name.replace('__Hrot', '')
        rec = poses.setdefault(base, dict(target=MB[base]['target'], scoring=MB[base]['scoring'], seed=MB[base]['seed'],
                                          score=MB[base]['score'], cluster_population=MB[base]['cluster_population'],
                                          tau_fold_manifest=MB[base]['tau_fold'], mmff_rigid=MB[base]['strain_mmff_rigid'],
                                          mmff_relaxed=MB[base]['strain_mmff_relaxed']))
        kind = 'hrot' if name.endswith('__Hrot') else 'asdocked'
        rec[f'strain_{kind}_kcal'] = (e[name]['E_Ha'] - E0) * HA_KCAL
        rec[f'sha256_ok_{kind}'] = hashlib.sha256(p.read_bytes()).hexdigest() == g['sha256']
        rec[f'tau_fold_{kind}'] = tau_fold(X)
        rec[f'X_{kind}'] = X
        rec['symbols'] = s
    for base, rec in poses.items():
        s = rec.pop('symbols')
        Xa, Xh = rec.pop('X_asdocked'), rec.pop('X_hrot')
        heavy = [i for i, el in enumerate(s) if el != 'H']
        rec['heavy_max_dev_hrot_vs_asdocked_A'] = float(np.linalg.norm(Xa[heavy] - Xh[heavy], axis=1).max())
        moved = [i for i, el in enumerate(s) if el == 'H' and np.linalg.norm(Xa[i] - Xh[i]) > 1e-3]
        # hydrogens that moved: which heavy atom carries them
        rec['moved_H_on'] = sorted({s[int(np.argmin([np.linalg.norm(Xh[i] - Xh[j]) if j != i else 9 for j in range(len(s))]))] for i in moved})
        rec['n_H_moved'] = len(moved)
    hr = [r_['strain_hrot_kcal'] for r_ in poses.values()]
    ad = [r_['strain_asdocked_kcal'] for r_ in poses.values()]
    strain[sp] = dict(case=case, source=rel(RES / case / 'dock_strain.json'), E_opt_Ha=E0, E_opt_case_Ha=Eopt_case,
                      E_opt_diff_Ha=E0 - Eopt_case, tau_fold_opt=tau_fold(X0), n_poses=len(poses), poses=poses,
                      hrot_range=[min(hr), max(hr)], asdocked_range=[min(ad), max(ad)],
                      all_sha256_ok=all(r_['sha256_ok_hrot'] and r_['sha256_ok_asdocked'] for r_ in poses.values()),
                      max_tau_dev_vs_manifest=max(abs(r_['tau_fold_hrot'] - r_['tau_fold_manifest']) for r_ in poses.values()),
                      provenance=js['provenance'])
out['strain'] = strain
dump('d11_burst2.json', out)
