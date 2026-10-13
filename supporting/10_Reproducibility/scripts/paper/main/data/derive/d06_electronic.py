"""d06: frontier orbitals, Delta-SCF descriptors, TD-DFT, NLO, interaction energies, QTAIM energies, charges.

Own recomputation from the raw quantities: orbital energies parsed from the archived final.molden files; Delta-SCF
from the three total energies; TD-DFT wavelengths from the excitation energies; NLO alpha/beta rebuilt by finite
differences from the raw field/dipole list in nlo.json (Richardson-extrapolated, then beta_vec and beta_tot);
interaction energies from the counterpoise energies; E_HB = V/2 and the Emamian fits from the BCP data;
fragment charge sums from the atomic charges.
"""
from common import *

out = {}
# ---------------- FMO from molden
fmo = {}
for c in ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'cation_pcm_d3', 'cation_gas_d3', 'sulfanilamide_pcm_d3',
          'sulfanilamide_gas_d3', 'sulfonamidate_pcm_d3', 'h2po4_pcm_d3', 'h2po4_gas_d3', 'h3po4_pcm_d3',
          'h3po4_gas_d3', 'urea_gas_d3']:
    p = ARES / c / 'final.molden'
    if not p.exists():
        continue
    ene, occ = [], []
    for line in p.read_text().splitlines():
        s = line.strip()
        if s.startswith('Ene='):
            ene.append(float(s.split('=')[1]))
        elif s.startswith('Occup='):
            occ.append(float(s.split('=')[1]))
    ene, occ = np.array(ene), np.array(occ)
    homo = ene[occ > 1.0].max(); lumo = ene[occ < 0.5].min()
    ref = json.loads((RES / c / 'fmo.json').read_text())
    fmo[c] = dict(n_mo=len(ene), n_occ=int((occ > 1).sum()), HOMO_eV=homo * HA_EV, LUMO_eV=lumo * HA_EV,
                  gap_eV=(lumo - homo) * HA_EV, dHOMO_vs_fmo_json=homo * HA_EV - ref['HOMO_eV'],
                  dLUMO_vs_fmo_json=lumo * HA_EV - ref['LUMO_eV'])
out['fmo'] = fmo

# ---------------- Delta-SCF
dscf = {}
for c in ['sfa_pair_pcm_d3', 'cation_pcm_d3', 'cation_gas_d3', 'sulfanilamide_pcm_d3', 'sulfanilamide_gas_d3']:
    d = json.loads((RES / c / 'dscf.json').read_text())
    E0, Em, Ep = d['E_N_Ha'], d['N_minus_1']['E_Ha'], d['N_plus_1']['E_Ha']
    I, A = (Em - E0) * HA_EV, (E0 - Ep) * HA_EV
    chi = (I + A) / 2; eta = (I - A) / 2
    rec = dict(I=I, A=A, chi=chi, mu=-chi, eta=eta, S=1 / (2 * eta), omega=chi ** 2 / (2 * eta), dNmax=chi / eta,
               fundamental_gap=I - A, S2_cation=d['N_minus_1']['S2'], S2_anion=d['N_plus_1']['S2'],
               environment=d['provenance']['environment'])
    kI, kA = -d['HOMO_eV'], -d['LUMO_eV']
    kchi, keta = (kI + kA) / 2, (kI - kA) / 2
    rec['koopmans'] = dict(I=kI, A=kA, chi=kchi, eta=keta, S=1 / (2 * keta), omega=kchi ** 2 / (2 * keta),
                           gap=kI - kA)
    rec['max_abs_diff_vs_json'] = max(abs(I - d['delta_scf']['I_eV']), abs(A - d['delta_scf']['A_eV']),
                                      abs(eta - d['delta_scf']['hardness_eta_eV']),
                                      abs(rec['omega'] - d['delta_scf']['electrophilicity_omega_eV']))
    dscf[c] = rec
out['dscf'] = dscf

# ---------------- TD-DFT
td = {}
for c, fname in [('sfa_pair_pcm_d3', 'tddft_camb3lyp_noneq'), ('sfa_pair_pcm_d3', 'tddft_b3lyp_noneq'),
                 ('cation_pcm_d3', 'tddft_camb3lyp_noneq'), ('cation_pcm_d3', 'tddft_camb3lyp_eq'),
                 ('cation_pcm_d3', 'tddft_b3lyp_noneq'), ('sulfanilamide_pcm_d3', 'tddft_camb3lyp_noneq'),
                 ('sulfanilamide_pcm_d3', 'tddft_camb3lyp_eq'), ('sulfanilamide_pcm_d3', 'tddft_b3lyp_noneq')]:
    d = json.loads((RES / c / f'{fname}.json').read_text())
    st = []
    for s in d['states']:
        ct = s.get('ct') or {}
        st.append(dict(state=s['state'], E_eV=s['E_eV'], nm=1239.841984 / s['E_eV'], f=s['f'],
                       dominant=[(x['label'], x['weight']) for x in s['dominant'][:3]],
                       nto1=s.get('nto_first_pair_weight'), CT=ct.get('CT_number'), D_CT=ct.get('D_CT_A'),
                       hole_on_anion=(ct.get('hole_on') or {}).get('dihydrogen phosphate anion'),
                       electron_on_anion=(ct.get('electron_on') or {}).get('dihydrogen phosphate anion')))
    strong = [s for s in st if s['f'] >= 0.1]
    brightest = max(st, key=lambda s: s['f'])
    td[f'{c}/{fname}'] = dict(functional=d['functional'], solvation=d['solvation'], nstates=d['nstates'],
                              all_converged=d['all_converged'], states=st, S1=st[0],
                              first_f_ge_0p1=strong[0] if strong else None, brightest=brightest,
                              spectrum_peak_nm=d['spectrum_peak_nm'], fwhm_eV=d.get('fwhm_eV'),
                              environment=d['environment'], ground_HOMO_eV=d.get('HOMO_eV'), ground_LUMO_eV=d.get('LUMO_eV'),
                              max_CT_first5=max((s['CT'] or 0) for s in st[:5]) if st[0]['CT'] is not None else None)
out['tddft'] = td

# ---------------- NLO from the raw finite-field dipoles
BETA_ESU = 8.639220664380406e-33


def nlo_from_fields(fields):
    mu0 = None
    tab = {}
    for f in fields:
        F = tuple(np.round(np.array(f['F_au']), 6))
        tab[F] = np.array(f['mu_au'])
    mu0 = tab[(0.0, 0.0, 0.0)]

    def mu(v):
        return tab[tuple(np.round(np.array(v, float), 6) + 0.0)]

    def tensors(F):
        a = np.zeros((3, 3)); b = np.zeros((3, 3, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = F
            mp, mm = mu(e), mu(-e)
            a[:, j] = (mp - mm) / (2 * F)
            b[:, j, j] = (mp + mm - 2 * mu0) / F ** 2
        for j in range(3):
            for k in range(j + 1, 3):
                ej = np.zeros(3); ej[j] = F
                ek = np.zeros(3); ek[k] = F
                v = (mu(ej + ek) - mu(ej - ek) - mu(-ej + ek) + mu(-ej - ek)) / (4 * F ** 2)
                b[:, j, k] = v; b[:, k, j] = v
        return a, b
    a1, b1 = tensors(0.001)
    a2, b2 = tensors(0.002)
    a = (4 * a1 - a2) / 3
    b = (4 * b1 - b2) / 3
    a = 0.5 * (a + a.T)
    # full permutational symmetrisation (static, Kleinman)
    bs = sum(b.transpose(p) for p in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]) / 6
    bvec = np.array([sum(bs[i, j, j] for j in range(3)) for i in range(3)])
    # the unsymmetrised definition beta_i = 1/3 sum_j (b_ijj + b_jij + b_jji)
    bvec_u = np.array([sum(b[i, j, j] + b[j, i, j] + b[j, j, i] for j in range(3)) / 3 for i in range(3)])
    ev = np.linalg.eigvalsh(a)
    aiso = np.trace(a) / 3
    daniso = np.sqrt(0.5 * ((ev[0] - ev[1]) ** 2 + (ev[1] - ev[2]) ** 2 + (ev[2] - ev[0]) ** 2))
    mun = np.linalg.norm(mu0)
    return dict(mu_au=float(mun), mu_D=float(mun * AU_DEBYE), alpha_iso_au=float(aiso), alpha_aniso_au=float(daniso),
                alpha_iso_A3=float(aiso * 0.1481847114721629), beta_tot_au=float(np.linalg.norm(bvec_u)),
                beta_tot_sym_au=float(np.linalg.norm(bvec)), beta_tot_esu=float(np.linalg.norm(bvec_u) * BETA_ESU),
                beta_par_mu_au=float(bvec_u @ mu0 / mun), beta_tot_B_conv_au=float(np.linalg.norm(bvec_u) / 2),
                beta_F_only_au=float(np.linalg.norm([sum(b1[i, j, j] + b1[j, i, j] + b1[j, j, i] for j in range(3)) / 3
                                                     for i in range(3)])))


nlo = {}
for c, fname in [('sfa_pair_pcm_d3', 'nlo.json'), ('sfa_pair_gas_d3', 'nlo_geom-sfa_pair_pcm_d3.json'),
                 ('sulfanilamide_pcm_d3', 'nlo.json'), ('urea_gas_d3', 'nlo.json')]:
    d = json.loads((RES / c / fname).read_text())
    r = nlo_from_fields(d['fields'])
    r['json'] = dict(mu_D=d['mu_total_debye'], alpha_iso_au=d['alpha_iso_au'], beta_tot_au=d['beta_tot_au'],
                     beta_tot_esu=d['beta_tot_esu'], environment=d['provenance']['environment'],
                     geometry=d['provenance']['geometry'])
    r['file'] = f'molecules/sfa/results/{c}/{fname}'
    nlo[c] = r
u = nlo['urea_gas_d3']['beta_tot_au']
for c in nlo:
    nlo[c]['ratio_to_urea'] = nlo[c]['beta_tot_au'] / u
    nlo[c]['alpha_ratio_to_urea'] = nlo[c]['alpha_iso_au'] / nlo['urea_gas_d3']['alpha_iso_au']
out['nlo'] = nlo

# ---------------- interaction energies
inter = {}
for c, fname in [('sfa_pair_pcm_d3', 'interaction.json'), ('sfa_pair_gas_d3', 'interaction.json'),
                 ('sfa_pair_gas_d3', 'interaction_geom-sfa_pair_pcm_d3.json')]:
    d = json.loads((RES / c / fname).read_text())
    env = 'pcm' if 'pcm' in d['energies_Ha'] else 'gas'
    E = d['energies_Ha'][env]
    fr = list(d['fragments'].keys())
    D3 = d['D3_Ha']
    d3int = D3['AB'] - sum(D3[f] for f in fr)
    raw = (E['AB'] - sum(E[f'{f}|own'] for f in fr)) * HA_KCAL + d3int * HA_KCAL
    cp = (E['AB'] - sum(E[f'{f}|AB'] for f in fr)) * HA_KCAL + d3int * HA_KCAL
    edef = {}
    for f in fr:
        # deformation: fragment in the complex geometry (own basis, with its D3) minus the relaxed fragment (with its D3)
        edef[f] = ((E[f'{f}|own'] + D3[f]) - (E[f'{f}|relaxed'] + E[f'{f}|relaxed_D3'])) * HA_KCAL
    inter[f'{c}/{fname}'] = dict(partition=d['partition'], environment=env, E_int_raw_kcal=raw, E_int_CP_kcal=cp,
                                 BSSE_kcal=cp - raw, E_int_D3_kcal=d3int * HA_KCAL, E_def_kcal=edef,
                                 E_bind_kcal=cp + sum(edef.values()),
                                 json=d['results'], gas_bsse_crosscheck=d.get('gas_crosscheck_for_pcm'))
out['interaction'] = inter

# ---------------- QTAIM H-bond energies
qt = {}
for c in ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3']:
    d = json.loads((RES / c / 'mwfn_extra/qtaim.json').read_text())
    rows = []
    for h in d['hbonds']:
        rows.append(dict(name=h['name'], D=h['donor'], H=h['hydrogen'], A=h['acceptor'], DA=h['D..A_A'],
                         HA=h['H..A_A'], DHA=h['DHA_deg'], rho=h['rho'], lap=h['lap'], G=h['G'], V=h['V'], Hd=h['H'],
                         mGV=-h['G'] / h['V'], E_HB_kcal=h['V'] / 2 * HA_KCAL, E_HB_kJ=h['V'] / 2 * HA_KJ,
                         BE_em_neutral=-223.08 * h['rho'] + 0.7423, BE_em_charged=-332.34 * h['rho'] - 1.0661,
                         fit_used=h['Emamian_fit_used'], json_E_HB=h['E_HB_Espinosa_kcal'], nature=h['nature_Rozas2000']))
    qt[c] = dict(n_cp=d['n_cp'], counts=d['counts'], poincare_hopf=d['poincare_hopf'], hbonds=rows)
out['qtaim'] = qt

# ---------------- charges
ch = {}
for c in ['sfa_pair_pcm_d3', 'sfa_pair_gas_d3', 'cation_pcm_d3', 'sulfanilamide_pcm_d3', 'sulfonamidate_pcm_d3',
          'h2po4_pcm_d3', 'urea_gas_d3']:
    p = RES / c / 'mwfn_extra/charges.json'
    if not p.exists():
        continue
    d = json.loads(p.read_text())
    rec = {}
    for sch, q in d['atomic'].items():
        q = np.array(q)
        if c.startswith('sfa_pair'):
            rec[sch] = dict(cation_sum=float(q[:20].sum()), anion_sum=float(q[20:].sum()), total=float(q.sum()))
            if c == 'sfa_pair_gas_d3':           # after transfer: cation minus H1A, anion plus H1A
                hq = q[15]
                rec[sch]['sulfanilamide_sum'] = float(q[:20].sum() - hq)
                rec[sch]['h3po4_sum'] = float(q[20:].sum() + hq)
        else:
            rec[sch] = dict(total=float(q.sum()))
    ch[c] = dict(sums=rec, fragment_sums_json=d.get('fragment_sums'))
out['charges'] = ch
dump('d06_electronic.json', out)
