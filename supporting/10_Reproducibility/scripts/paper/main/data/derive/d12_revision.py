"""d12: numbers added in the revision after the two referee reports (29 Sep 2026), each from stored results; no new
VM calculation.  Reads only; writes derived/d12_revision.json.

(1) Sulfonamide torsion of the ligands in the deposited protein complexes (referees: dock M6, chem M9c).  Read from the
    mmCIF files in dock/pdb/ with gemmi (every copy of the ligand, every altloc) and, as a cross-check, from the
    reference SDF files the docking used (dock/T*/ref_*.sdf).  tau_fold = max(|C3-C4-S-N|, |C5-C4-S-N|) with C4 the
    ring carbon bonded to S and C3/C5 its ring neighbours, found by distance (own code, no RDKit), so atom names do not
    matter.  The copy used as the docking reference is the one closest to the zinc (hCA II) or the one in the chain
    that was docked (3TZF chain B).
(2) Docking scores relative to p-aminobenzoate (referees: dock M1, chem M7): Delta = S(ligand) - S(pABA) within one
    receptor and scoring function, and the double difference DD = Delta(7S2K, Sul2) - Delta(3TZF, DHPS), from the
    independent score re-derivation of d09 (derived/d09_dock.json, mean of the top-pose score over seeds).
(3) hCA II ranks over zinc-bound top poses only (dock M3, chem M8d): a seed counts when its top pose has
    Zn-N(sulfonamide) <= 2.3 A and Thr199 OG1...N <= 3.3 A (the gate definition), from the per-seed site checks in
    dock/results/production_summary.json.
(4) Core RMSDs that the referees asked for (chem m23, dock M1.6): read from production_summary.json.
(5) Proton transfer between the separated species in IEF-PCM against the measured Delta-pKa (chem M5): electronic
    energies of the four optimised IEF-PCM species (opt_summary.json) and Delta-G = -RT ln10 Delta-pKa, with the
    measured aqueous pKa values chosen by the authors on 29 Sep 2026: anilinium of sulfanilamide 2.02 at 25 C (Greenbaum
    1955, UV absorbance, as evaluated in Perrin's IUPAC compilation, digitised by Zheng et al. 2026) and H3PO4 pK1
    2.148 at 25 C (Bates 1951, emf, full text read). The earlier 2.37 (14 - Merck pKb 11.63) was an indirect value.
(6) Delivered dose: mass fraction of sulfanilamide in SAPH from IUPAC 2021 standard atomic weights (dock M8).
(7) Rigid DFT strain of the docked poses referred to the least-rotated pose of the same species (removes the offset
    common to all poses; dock M5, chem M9a), from dock_strain.json via derived/d11_burst2.json; and, when the Pi single
    points of the rotor-reoptimised optimum exist (derived/rev_strainref/), the strain against that reference.
(8) Order-of-magnitude checks on the constrained periodic neutral (chem M2): the energy held in a compressed O5-H1A
    bond, 1/2 F dx, for the residual force on H1A and a relaxation dx between the O5-H1A of the IEF-PCM neutral
    minimum minus 0.98 A and 0.09 A; and the relaxation energy P^2 V / (2 B) of each cell at its residual pressure for
    bulk moduli of 10 and 20 GPa (difference constrained - salt).
"""
import glob
import statistics as stt

import gemmi

from common import *

R_GAS = 1.98720425864083e-3   # kcal/(mol K)
T298 = 298.15

# ------------------------------------------------------------------ (1) co-crystal torsions
CUT = {('C', 'C'): 1.65, ('C', 'S'): 1.9, ('N', 'S'): 1.8, ('O', 'S'): 1.6}


def bonded(e1, x1, e2, x2):
    c = CUT.get(tuple(sorted((e1, e2))))
    return c is not None and dist(x1, x2) < c


def tau_of(atoms):
    """atoms: [(element, xyz)] of one ligand copy (heavy atoms). Returns tau_fold, C4-S-N angle, both torsions."""
    S = [a for a in atoms if a[0] == 'S']
    out = []
    for s in S:
        cs = [a for a in atoms if a[0] == 'C' and bonded('C', a[1], 'S', s[1])]
        ns = [a for a in atoms if a[0] == 'N' and bonded('N', a[1], 'S', s[1])]
        if len(cs) != 1 or len(ns) != 1:
            continue
        c4, n = cs[0], ns[0]
        ring = [a for a in atoms if a[0] == 'C' and a is not c4 and bonded('C', a[1], 'C', c4[1])]
        if len(ring) != 2:
            continue
        t1 = dihedral(ring[0][1], c4[1], s[1], n[1])
        t2 = dihedral(ring[1][1], c4[1], s[1], n[1])
        out.append(dict(tau_fold=max(abs(t1), abs(t2)), tau_a=t1, tau_b=t2, CSN=angle(c4[1], s[1], n[1]),
                        N_xyz=n[1].tolist()))
    assert len(out) == 1, out
    return out[0]


def sdf_heavy(path):
    lines = Path(path).read_text().splitlines()
    na = int(lines[3][:3])
    at = []
    for l in lines[4:4 + na]:
        p = l.split()
        if p[3] != 'H':
            at.append((p[3], np.array([float(p[0]), float(p[1]), float(p[2])])))
    return at


cocrystal = {}
for code, res, ref, note in [('6RL9', 'SAN', 'T3_6RL9/ref_SAN.sdf', 'hCA II + sulfanilamide'),
                             ('1AJ0', 'SAN', 'T2_1AJ0/ref_SAN.sdf', 'E. coli DHPS + sulfanilamide'),
                             ('3TZF', '08D', 'T1_3TZF/ref_08D.sdf', 'Y. pestis DHPS + sulfamethoxazole')]:
    st = gemmi.read_structure(str(DOCK / 'pdb' / f'{code}.cif.gz'))
    st.setup_entities()
    model = st[0]
    zn = [a.pos for ch in model for r in ch if r.name == 'ZN' for a in r]
    copies = []
    for ch in model:
        for r in ch:
            if r.name != res:
                continue
            alts = sorted({a.altloc for a in r if a.altloc != '\x00'}) or ['']
            for alt in alts:
                at = [(a.element.name, np.array([a.pos.x, a.pos.y, a.pos.z])) for a in r
                      if a.element.name != 'H' and (a.altloc in ('\x00', alt) or alt == '')]
                t = tau_of(at)
                d_zn = min((dist(t['N_xyz'], [p.x, p.y, p.z]) for p in zn), default=None)
                copies.append(dict(chain=ch.name, seqid=str(r.seqid), altloc=alt or None, occupancy=float(
                    min(a.occ for a in r)), Zn_N=d_zn, **{k: v for k, v in t.items() if k != 'N_xyz'}))
    tr = tau_of(sdf_heavy(DOCK / ref))
    # the copy that the docking used = the one whose torsion and angle match the reference SDF
    match = min(copies, key=lambda c: abs(c['tau_fold'] - tr['tau_fold']) + abs(c['CSN'] - tr['CSN']))
    cocrystal[code] = dict(ligand=res, what=note, copies=copies, ref_sdf=f'dock/{ref}', tau_fold_ref_sdf=tr['tau_fold'],
                           CSN_ref_sdf=tr['CSN'], used=match, tau_fold=match['tau_fold'], CSN=match['CSN'],
                           agree_sdf_cif=abs(match['tau_fold'] - tr['tau_fold']) < 0.05)

# ------------------------------------------------------------------ (2) scores relative to pABA
D9 = json.loads((OUT / 'd09_dock.json').read_text())
SC = D9['scores']


def score(target, lig, fn):
    s = SC[f'{target}/{lig}__{fn}__prod']
    return s['mean'], s['sd']


rel_paba = {}
for fn in ('vina', 'vinardo'):
    for lig in ('sfa_neutral', 'sfa_anion', 'sfa_cation', 'smx_neutral', 'smx_anion', 'azm_neutral', 'azm_anion'):
        row = {}
        for tg in ('T1_3TZF', 'T5_7S2K'):
            (a, sa), (p, sp) = score(tg, lig, fn), score(tg, 'paba_anion', fn)
            row[tg] = dict(delta=a - p, sd=float(np.hypot(sa, sp)), ligand=a, paba=p)
        row['DD'] = row['T5_7S2K']['delta'] - row['T1_3TZF']['delta']
        row['DD_sd'] = float(np.hypot(row['T5_7S2K']['sd'], row['T1_3TZF']['sd']))
        rel_paba[f'{lig}.{fn}'] = row
# the score change of each ligand between the two receptors
shift = {f'{lig}.{fn}': score('T5_7S2K', lig, fn)[0] - score('T1_3TZF', lig, fn)[0]
         for fn in ('vina', 'vinardo')
         for lig in ('sfa_neutral', 'smx_neutral', 'smx_anion', 'paba_anion', 'azm_neutral')}
# ties: species within 0.1 kcal/mol of sulfanilamide
ties = {}
for key, r in D9['ranks'].items():
    order = r['order']
    s0 = dict(order)['sfa_neutral'] if 'sfa_neutral' in dict(order) else None
    if s0 is None:
        continue
    ties[key] = dict(within_0p1=[k for k, v in order if abs(v - s0) <= 0.1 and k != 'sfa_neutral'],
                     rank=r.get('rank_sfa_neutral'), n=r['n'],
                     rank_best=1 + sum(1 for k, v in order if v < s0 - 0.1),
                     rank_worst=sum(1 for k, v in order if v <= s0 + 0.1))

# ------------------------------------------------------------------ (3) zinc-filtered hCA II ranks
PS = json.loads((DOCK / 'results' / 'production_summary.json').read_text())
zinc = {}
for tg in ('T3_6RL9', 'T4_3HS4'):
    for lig in ('sfa_neutral', 'smx_neutral', 'azm_neutral'):
        J = PS['jobs'][f'{tg}/{lig}__ad4zn__prod']
        ok = [p['score'] for p in J['per_seed']
              if p['site_checks']['Zn-N(sulfonamide)'] <= 2.3 and p['site_checks']['Thr199 OG1 H-bond'] <= 3.3]
        allsc = [p['score'] for p in J['per_seed']]
        zinc[f'{tg}.{lig}'] = dict(n_zinc=len(ok), n=len(allsc), mean_zinc=stt.mean(ok) if ok else None,
                                   mean_all=stt.mean(allsc),
                                   ZnN_mean=stt.mean(p['site_checks']['Zn-N(sulfonamide)'] for p in J['per_seed']))
    ranked = sorted([(k.split('.')[1], v['mean_zinc']) for k, v in zinc.items()
                     if k.startswith(tg) and v['mean_zinc'] is not None], key=lambda x: x[1])
    zinc[f'{tg}.order'] = ranked
    zinc[f'{tg}.gap_azm_sfa'] = dict(ranked)['azm_neutral'] - dict(ranked)['sfa_neutral']

# AutoDock4Zn validation of the anions: mean Zn-N of the top poses (dock/analysis/summary.json)
VS = json.loads((DOCK / 'analysis' / 'summary.json').read_text())
zinc['val_anions'] = {run: dict(ZnN=VS[run]['stats']['site_checks']['Zn-N(sulfonamide)']['distance_mean'],
                                pass_fraction=VS[run]['stats']['site_checks']['Zn-N(sulfonamide)']['pass_fraction'])
                      for run in ('T4_3HS4/azm_anion__ad4zn', 'T4_3HS4/azm_neutral__ad4zn', 'T3_6RL9/sfa_anion__ad4zn')
                      if run in VS}

# ------------------------------------------------------------------ (4) core RMSDs from the production summary
def core(job, ref, name):
    c = PS['jobs'][job]['rmsd_to_crystal'][ref]['cores'][name]
    return dict(mean=c['mean'], sd=c['sd'], f1=c['f_le_1A'], n=c['n'])


cores = {
    'T1_3TZF.smx_anion.vina.aniline_SO2N': core('T1_3TZF/smx_anion__vina__prod', '08D (3TZF, own)',
                                                'aniline-SO2N core (11 atoms)'),
    'T1_3TZF.sfa_neutral.vinardo.aniline_SO2N': core('T1_3TZF/sfa_neutral__vinardo__prod', '08D (3TZF, own)',
                                                     'aniline-SO2N core (11 atoms)'),
    'T1_3TZF.sfa_neutral.vina.aniline_SO2N': core('T1_3TZF/sfa_neutral__vina__prod', '08D (3TZF, own)',
                                                  'aniline-SO2N core (11 atoms)'),
    'T1_3TZF.paba_anion.vina.aniline_XO2': core('T1_3TZF/paba_anion__vina__prod', '08D (3TZF, own)',
                                                'aniline-X(O)2 core, X = S or C (10 atoms)'),
    'T5_7S2K.paba_anion.vina.aniline_XO2': core('T5_7S2K/paba_anion__vina__prod', '78H (7S2K, own)',
                                                'aniline-X(O)2 core, X = S or C (10 atoms)'),
    'T5_7S2K.sfa_neutral.vina.aniline_XO2': core('T5_7S2K/sfa_neutral__vina__prod', '78H (7S2K, own)',
                                                 'aniline-X(O)2 core, X = S or C (10 atoms)'),
}

# ------------------------------------------------------------------ (5) separated species vs Delta-pKa
E = {c: json.loads((RES / c / 'opt_summary.json').read_text())['E_final_Ha']
     for c in ('sulfanilamide_pcm_d3', 'h3po4_pcm_d3', 'cation_pcm_d3', 'h2po4_pcm_d3')}
dE_sep = ((E['sulfanilamide_pcm_d3'] + E['h3po4_pcm_d3']) - (E['cation_pcm_d3'] + E['h2po4_pcm_d3'])) * HA_KCAL
pKa_anilinium = 2.02   # greenbaum1955comparison via the IUPAC compilation (zheng2026digitized, entry perrin613)
pKa_h3po4 = 2.148      # bates1951first: "2.148 +/- 0.005 at 25 deg"
dpKa = pKa_anilinium - pKa_h3po4
RTln10 = R_GAS * T298 * math.log(10)
# proton transfer cation + H2PO4- -> sulfanilamide + H3PO4 in water: K = Ka(anilinium)/Ka1(H3PO4) = 10^(pKa1 - pKa)
dG_exp = RTln10 * dpKa     # negative: the neutral pair lies (slightly) below the ions
sep = dict(E_Ha=E, dE_sep_kcal=dE_sep, pKa_anilinium=pKa_anilinium, pKa_h3po4=pKa_h3po4, dpKa=dpKa,
           RTln10_kcal=RTln10, dG_exp_kcal=dG_exp, continuum_error_kcal=dE_sep - dG_exp,
           note='dG_exp is Delta-G(neutral pair - ions) in water from the two pKa values; dE_sep is the IEF-PCM '
                'electronic energy of the same exchange between separated, relaxed species (no ZPE or thermal terms)')

# ------------------------------------------------------------------ (6) delivered dose
AW = dict(C=12.011, H=1.008, N=14.007, O=15.999, P=30.974, S=32.06)   # IUPAC 2021 abridged standard atomic weights


def mw(f):
    return sum(AW[e] * n for e, n in f.items())


M_sfa = mw(dict(C=6, H=8, N=2, O=2, S=1))
M_h3po4 = mw(dict(H=3, P=1, O=4))
M_saph = mw(dict(C=6, H=11, N=2, O=6, P=1, S=1))
dose = dict(M_sulfanilamide=M_sfa, M_H3PO4=M_h3po4, M_SAPH=M_saph, fraction=M_sfa / M_saph,
            check_sum=M_saph - (M_sfa + M_h3po4))

# ------------------------------------------------------------------ (7) strain relative to the least-rotated pose
D11 = json.loads((OUT / 'd11_burst2.json').read_text())
strain = {}
for sp, blk in D11['strain'].items():
    poses = blk['poses']
    lo_key = min(poses, key=lambda k: poses[k]['tau_fold_hrot'])
    e0 = poses[lo_key]['strain_hrot_kcal']
    inc = {k: v['strain_hrot_kcal'] - e0 for k, v in poses.items()}
    strain[sp] = dict(reference_pose=lo_key, reference_tau=poses[lo_key]['tau_fold_hrot'], reference_E=e0,
                      reference_mmff_relaxed=poses[lo_key].get('mmff_relaxed'),
                      increment_max=max(inc.values()), increment_min=min(inc.values()),
                      per_pose={k: dict(tau=poses[k]['tau_fold_hrot'], strain=poses[k]['strain_hrot_kcal'],
                                        increment=inc[k], mmff_relaxed=poses[k].get('mmff_relaxed'),
                                        mmff_rigid=poses[k].get('mmff_rigid')) for k in poses})
# Pi single points of the IEF-PCM optimum with its polar-H rotors re-optimised by MMFF94s (same procedure as the poses)
SR = OUT / 'rev_strainref'
strainref = {}
for sp, case in (('sfa_neutral', 'sulfanilamide_pcm_d3'), ('sfa_anion', 'sulfonamidate_pcm_d3'),
                 ('sfa_cation', 'cation_pcm_d3')):
    f = SR / f'strainref_{case}.json'
    if not f.exists():
        continue
    J = json.loads(f.read_text())
    en = J['energies']['b3lyp-d3bj']
    e_opt_pi, e_ref = en['opt']['E_Ha'], en['opt__Hrot']['E_Ha']
    e_opt_vm = D11['strain'][sp]['E_opt_Ha']
    off = (e_ref - e_opt_pi) * HA_KCAL
    per = {k: v['strain_hrot_kcal'] - off for k, v in D11['strain'][sp]['poses'].items()}
    strainref[sp] = dict(source=rel(f), E_opt_pi_Ha=e_opt_pi, E_opt_vm_Ha=e_opt_vm,
                         E_opt_diff_Ha=e_opt_pi - e_opt_vm, E_ref_Ha=e_ref, offset_kcal=off,
                         strain_vs_ref=per, strain_vs_ref_min=min(per.values()), strain_vs_ref_max=max(per.values()),
                         host=J['provenance'].get('host'), threads=J['provenance'].get('threads'))

# ------------------------------------------------------------------ (8) periodic order-of-magnitude checks
P = D11['periodic']
PQJ = json.loads((RES / 'periodic_qe' / 'periodic_qe.json').read_text())
F_eVA = P['iirfix_constraint']['F_H1A']['fu1']['F_H1A_along_H1A_to_N1_Ry_bohr'] * RY_EV / BOHR_A
import re as _re
V_cell = None
cifc = CIF.read_text()


def cifnum(tag):
    m = _re.search(rf'^{tag}\s+([0-9.]+)', cifc, _re.M)
    return float(m.group(1))


a, b, c, beta = (cifnum('_cell_length_a'), cifnum('_cell_length_b'), cifnum('_cell_length_c'),
                 cifnum('_cell_angle_beta'))
V_cell = a * b * c * math.sin(math.radians(beta))
p_salt = PQJ['relax']['i'].get('pressure_kbar_final')
p_cn = PQJ['relax']['iirfix'].get('pressure_kbar_final')
OH_pcm = json.loads((OUT / 'd04_scan.json').read_text())
OH_pcm = [p for p in OH_pcm['points'] if p.get('chain') == 'B' and p.get('index') == 0]
OH_pcm = OH_pcm[0]['d_H1A_O5'] if OH_pcm else 1.0264
dx_lo = OH_pcm - 0.98      # O5-H1A of the IEF-PCM neutral minimum (scan end point B00) minus the held 0.98 A
dx_hi = 0.09
peri = dict(F_H1A_eVA=F_eVA, OH_pcm_neutral=OH_pcm, dx_lo=dx_lo, dx_hi=dx_hi,
            OH_relief_lo_kJ=0.5 * F_eVA * dx_lo * EV_KJ if F_eVA else None,
            OH_relief_hi_kJ=0.5 * F_eVA * dx_hi * EV_KJ if F_eVA else None,
            V_cell_A3=V_cell, p_salt_kbar=p_salt, p_constrained_kbar=p_cn)
for B in (10.0, 20.0):
    # relaxation energy of a cell at residual pressure p: p^2 V / (2 B); per formula unit (Z = 2), kJ/mol
    def erel(p_kbar):
        p = p_kbar * 1e8     # Pa
        return p * p * V_cell * 1e-30 / (2 * B * 1e9) * 6.02214076e23 / 1000 / 2
    peri[f'stress_B{int(B)}_kJ'] = erel(p_cn) - erel(p_salt)

dump('d12_revision.json', dict(cocrystal=cocrystal, rel_paba=rel_paba, shift_7S2K_minus_3TZF=shift, ties=ties,
                               zinc=zinc, cores=cores, separated=sep, dose=dose, strain=strain, strainref=strainref,
                               periodic=peri))
