"""Facts added in the revision (29 Sep 2026, after the referee reports ref_chem.md and ref_dock.md; section 31 of
facts.md), executed by build_facts.py after facts_burst2.py (exec in its namespace: add(), section(), FACTS, J, r, RES,
DOCK, DER, CAVEATS, MISSING are the ones defined there).  Every number comes from derive/d12_revision.py
(derived/d12_revision.json), which reads stored results only.  The one new calculation is a pair of B3LYP-D3(BJ)/
6-311++G(d,p) IEF-PCM single points per sulfanilamide species run on the Pi with the project's own dft/fsp.py
(derived/rev_strainref/): the IEF-PCM optimum as stored (identity check against the VM energy) and the same geometry with
its polar-H rotors re-optimised by the MMFF94s procedure that placed the hydrogens of the docked poses.
"""
D12 = J(DER / 'd12_revision.json')
D12S = 'derived/d12_revision.json'

section('31. Revision after the referee reports (29 Sep 2026)',
        'Numbers asked for by the referees and derived from stored results (derive/d12_revision.py): the sulfonamide torsion of '
        'the ligands in the deposited protein complexes (gemmi on dock/pdb/*.cif.gz, cross-checked against dock/T*/ref_*.sdf); '
        'docking scores relative to p-aminobenzoate and their DHPS-to-Sul2 double difference (d09 re-derived scores); hCA II ranks '
        'over zinc-bound top poses; core RMSDs; the separated-species proton-transfer energy in IEF-PCM against the measured '
        'ΔpKa; the delivered dose; the DFT strain referred to the least-rotated pose and to a reference with the same hydrogen '
        'placement (Pi single points); and two order-of-magnitude estimates for the constrained periodic neutral.')

# ---------------------------------------------------------------- co-crystal torsions
CC = D12['cocrystal']
for code, lab in (('6RL9', 'hCA II + sulfanilamide (SAN 304, zinc-bound)'), ('1AJ0', 'E. coli DHPS + sulfanilamide'),
                  ('3TZF', 'Y. pestis DHPS chain B + sulfamethoxazole (08D)')):
    c = CC[code]
    add(f'dock.cocrystal.{code}.tau', c['tau_fold'], 'deg', 1, f'{D12S}:cocrystal.{code}.tau_fold (dock/pdb/{code}.cif.gz; '
        f'{c["ref_sdf"]} gives {c["tau_fold_ref_sdf"]:.2f}°)', f'τ_fold of the ligand in the deposited complex {code}: {lab}',
        check=c['tau_fold_ref_sdf'], tol=0.05)
    add(f'dock.cocrystal.{code}.CSN', c['CSN'], 'deg', 1, f'{D12S}:cocrystal.{code}.CSN', f'C4–S–N angle of the ligand in {code}',
        check=c['CSN_ref_sdf'], tol=0.05)
u = CC['6RL9']['used']
add('dock.cocrystal.6RL9.ZnN', u['Zn_N'], 'A', 2, f'{D12S}:cocrystal.6RL9.used.Zn_N', 'Zn–N of the docking reference copy SAN 304 in 6RL9')
others = [c for c in CC['6RL9']['copies'] if c['seqid'] != u['seqid']]
add('dock.cocrystal.6RL9.other_tau', [c['tau_fold'] for c in others], 'deg', 1, f'{D12S}:cocrystal.6RL9.copies',
    'τ_fold of the two further sulfanilamide copies in 6RL9 (SAN 305, altloc A, and SAN 306), which lie away from the zinc')
add('dock.cocrystal.6RL9.other_ZnN', [c['Zn_N'] for c in others], 'A', 1, f'{D12S}:cocrystal.6RL9.copies',
    'N···Zn of those two copies (secondary sites)')
add('dock.cocrystal.range', [min(CC[k]['tau_fold'] for k in CC), max(CC[k]['tau_fold'] for k in CC)], 'deg', 1,
    f'{D12S}:cocrystal', 'range of τ_fold over the three deposited complexes')
add('dock.cocrystal.overshoot_6RL9', FACTS['dock.tors.6RL9_ad4zn_mean']['value'] - CC['6RL9']['tau_fold'], 'deg', 0,
    'dock.tors.6RL9_ad4zn_mean − dock.cocrystal.6RL9.tau', 'mean τ_fold of the 6RL9 AutoDock4Zn poses minus the co-crystal value')
add('dock.cocrystal.diff_3TZF', FACTS['dock.strain.3TZF.sfa_neutral.vina.tau']['value'] - CC['3TZF']['tau_fold'], 'deg', 1,
    'dock.strain.3TZF.sfa_neutral.vina.tau − dock.cocrystal.3TZF.tau',
    'τ_fold of the representative sulfanilamide pose in 3TZF (Vina) minus that of crystal sulfamethoxazole')

# ---------------------------------------------------------------- scores relative to pABA
NAMES = dict(sfa_neutral='sulfanilamide', sfa_anion='sulfonamidate', sfa_cation='cation', smx_neutral='SMX⁰',
             smx_anion='SMX⁻', azm_neutral='AZM⁰', azm_anion='AZM⁻')
for key, row in D12['rel_paba'].items():
    lig, fn = key.split('.')
    for tg, code in (('T1_3TZF', '3TZF'), ('T5_7S2K', '7S2K')):
        add(f'dock.relpaba.{code}.{lig}.{fn}', row[tg]['delta'], 'kcal mol-1', 1, f'{D12S}:rel_paba.{key}.{tg}.delta',
            f'{code}, {fn}: score({NAMES[lig]}) − score(pABA⁻), means over seeds (seed SD of the difference {row[tg]["sd"]:.2f})')
    add(f'dock.dd.{lig}.{fn}', row['DD'], 'kcal mol-1', 1, f'{D12S}:rel_paba.{key}.DD',
        f'{fn}: double difference [score({NAMES[lig]}) − score(pABA⁻)]_7S2K − [same]_3TZF (positive: loses ground against pABA in Sul2)')
for key, v in D12['shift_7S2K_minus_3TZF'].items():
    lig, fn = key.split('.')
    add(f'dock.shift.{lig}.{fn}', v, 'kcal mol-1', 1, f'{D12S}:shift_7S2K_minus_3TZF.{key}',
        f'{fn}: score of {NAMES.get(lig, "pABA⁻")} in 7S2K minus that in 3TZF')
for key, t in D12['ties'].items():
    tg, fn = key.split('/')
    code = tg.split('_')[1]
    add(f'dock.ties.{code}.{fn}', [t['rank_best'], t['rank_worst']], None, 0, f'{D12S}:ties.{key}',
        f'{code} {fn}: range of ranks sulfanilamide can take when species within 0.1 kcal mol⁻¹ of it count as ties '
        f'(tied with {", ".join(t["within_0p1"]) or "none"}; printed rank {t["rank"]} of {t["n"]})')

# ---------------------------------------------------------------- zinc-filtered hCA II ranks
Z = D12['zinc']
for tg, code in (('T3_6RL9', '6RL9'), ('T4_3HS4', '3HS4')):
    for lig in ('sfa_neutral', 'azm_neutral', 'smx_neutral'):
        z = Z[f'{tg}.{lig}']
        add(f'dock.zinc.{code}.{lig}.n', z['n_zinc'], None, 0, f'{D12S}:zinc.{tg}.{lig}.n_zinc',
            f'{code} AutoDock4Zn {NAMES[lig]}: seeds whose top pose is zinc-bound (Zn–N ≤ 2.3 Å and Thr199 OG1···N ≤ 3.3 Å), of {z["n"]}')
        if z['mean_zinc'] is not None:
            add(f'dock.zinc.{code}.{lig}.score', z['mean_zinc'], 'kcal mol-1', 1, f'{D12S}:zinc.{tg}.{lig}.mean_zinc',
                f'{code} AutoDock4Zn {NAMES[lig]}: mean score over the zinc-bound top poses only')
    add(f'dock.zinc.{code}.gap', -Z[f'{tg}.gap_azm_sfa'], 'kcal mol-1', 2, f'{D12S}:zinc.{tg}.gap_azm_sfa',
        f'{code} AutoDock4Zn: acetazolamide scores better than sulfanilamide by this much over zinc-bound top poses')

for run, v in Z['val_anions'].items():
    tg, job = run.split('/')
    add(f'dock.val.{tg}.{job}.ZnN_mean', v['ZnN'], 'A', 2, f'{D12S}:zinc.val_anions.{run} (dock/analysis/summary.json)',
        f'validation {run}: mean Zn–N(sulfonamide) of the top poses (pass fraction {v["pass_fraction"]})')

# ---------------------------------------------------------------- core RMSDs
for key, c in D12['cores'].items():
    add(f'dock.core.{key}', c['mean'], 'A', 2, f'{D12S}:cores.{key} (dock/results/production_summary.json)',
        f'core RMSD {key} (mean over {c["n"]} seeds; SD {c["sd"]:.2f}; fraction ≤ 1 Å {c["f1"]:.1f})')
    add(f'dock.core.{key}.sd', c['sd'], 'A', 2, f'{D12S}:cores.{key}.sd', f'SD of the core RMSD {key}')

# ---------------------------------------------------------------- separated species vs Delta-pKa
S = D12['separated']
add('lit.pka_h3po4', S['pKa_h3po4'], None, 3, 'bates1951first (doi:10.6028/jres.047.017), full text read 29 Sep 2026: "2.148 ± 0.005 at 25°" '
    '(mean of 20 values, Table 3; emf without liquid junction, extrapolated to zero ionic strength)', 'first pKa of phosphoric acid in water at 25 °C')
add('saph.dpka', S['dpKa'], None, 2, f'{D12S}:separated.dpKa = admet.exp.pKa_anilinium − lit.pka_h3po4',
    'ΔpKa = pKa(anilinium of sulfanilamide, 2.02, Greenbaum 1955 via IUPAC) − pKa1(H3PO4, 2.148, Bates 1951); ±0.04 from the IUPAC "Approximate" rating')
add('pt.sep.dG_exp', S['dG_exp_kcal'], 'kcal mol-1', 2, f'{D12S}:separated.dG_exp_kcal = RT ln10 ΔpKa (298.15 K)',
    'ΔG°(sulfanilamide + H3PO4 − cation − H2PO4⁻) in water implied by the two pKa values (negative: the neutral pair is very slightly favoured)')
add('pt.sep.RTln10', S['RTln10_kcal'], 'kcal mol-1', 3, f'{D12S}:separated.RTln10_kcal', 'RT ln 10 at 298.15 K')
add('pt.sep.continuum_error', S['continuum_error_kcal'], 'kcal mol-1', 1, f'{D12S}:separated.continuum_error_kcal',
    'IEF-PCM electronic energy of the separated-species exchange (pt.sep.pcm) minus the ΔG° implied by the pKa values: the '
    'continuum favours the neutral side by this much (no ZPE or thermal terms in the computed value)', check=S['dE_sep_kcal'] - S['dG_exp_kcal'], tol=1e-6)

# ---------------------------------------------------------------- dose
Dd = D12['dose']
add('saph.dose.M_sfa', Dd['M_sulfanilamide'], 'g mol-1', 1, f'{D12S}:dose.M_sulfanilamide (IUPAC standard atomic weights)', 'molar mass of sulfanilamide')
add('saph.dose.M_saph', Dd['M_SAPH'], 'g mol-1', 1, f'{D12S}:dose.M_SAPH', 'molar mass of SAPH (C6H9N2O2S⁺·H2PO4⁻)',
    check=FACTS['xray.mr']['value'], tol=0.01)
add('saph.dose.fraction', Dd['fraction'], None, 3, f'{D12S}:dose.fraction', 'mass of sulfanilamide delivered per unit mass of SAPH')

# ---------------------------------------------------------------- strain relative to the least-rotated pose / same-H reference
ST = D12['strain']
for sp, lab in (('sfa_neutral', 'sulfanilamide'), ('sfa_anion', 'sulfonamidate'), ('sfa_cation', 'cation')):
    s = ST[sp]
    add(f'dock.strain.{sp}.inc_max', s['increment_max'], 'kcal mol-1', 1, f'{D12S}:strain.{sp}.increment_max',
        f'{lab}: largest rise of the rigid DFT strain over the representative poses, relative to the least-rotated pose '
        f'({s["reference_pose"]}, τ_fold {s["reference_tau"]:.1f}°)')
SR = D12.get('strainref', {})
for sp, lab in (('sfa_neutral', 'sulfanilamide'), ('sfa_anion', 'sulfonamidate'), ('sfa_cation', 'cation')):
    if sp not in SR:
        MISSING.append(f'dock.strainref.{sp}: Pi reference single points not found (derived/rev_strainref/)')
        continue
    s = SR[sp]
    add(f'dock.strainref.{sp}.Ediff', s['E_opt_diff_Ha'], 'Ha', 12, f'{D12S}:strainref.{sp}.E_opt_diff_Ha ({s["source"]})',
        f'{lab}: E(IEF-PCM optimum) on the Pi minus the VM value of dock_strain.json (identity check of the Pi run)')
    add(f'dock.strainref.{sp}.offset', s['offset_kcal'], 'kcal mol-1', 1, f'{D12S}:strainref.{sp}.offset_kcal',
        f'{lab}: E(optimum with MMFF94s-placed polar-H rotors) − E(optimum), B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM: the offset '
        'common to the rigid pose single points')
    add(f'dock.strainref.{sp}.min', s['strain_vs_ref_min'], 'kcal mol-1', 1, f'{D12S}:strainref.{sp}.strain_vs_ref_min',
        f'{lab}: smallest rigid DFT strain of the representative poses against the same-H reference')
    add(f'dock.strainref.{sp}.max', s['strain_vs_ref_max'], 'kcal mol-1', 1, f'{D12S}:strainref.{sp}.strain_vs_ref_max',
        f'{lab}: largest rigid DFT strain of the representative poses against the same-H reference')
    for k, v in s['strain_vs_ref'].items():
        tg, _, fn = k.split('__')[:3]
        add(f'dock.strainref.{tg.split("_")[1]}.{sp}.{fn}', v, 'kcal mol-1', 1, f'{D12S}:strainref.{sp}.strain_vs_ref.{k}',
            f'{lab}, {tg.split("_")[1]} {fn}: rigid DFT strain of the representative pose against the same-H reference')
    add(f'dock.strainref.{sp}.host', s['host'], None, None, f'{D12S}:strainref.{sp}.host',
        f'{lab}: machine of the reference single points (Raspberry Pi, {s["threads"]} threads; PySCF 2.14.0, dft/fsp.py)')

# ---------------------------------------------------------------- periodic order-of-magnitude checks
P = D12['periodic']
add('periodic.iirfix.dx_lo', P['dx_lo'], 'A', 3, f'{D12S}:periodic.dx_lo',
    'O5–H1A of the IEF-PCM neutral minimum (1.026 Å) minus the held 0.98 Å: lower bound of the O–H relaxation assumed below')
add('periodic.iirfix.dx_hi', P['dx_hi'], 'A', 2, f'{D12S}:periodic.dx_hi', 'upper O–H relaxation assumed (referee estimate)')
add('periodic.iirfix.OH_relief_lo', P['OH_relief_lo_kJ'], 'kJ mol-1', 0, f'{D12S}:periodic.OH_relief_lo_kJ = ½ F dx_lo',
    'energy held in a compressed O5–H1A if the residual force on H1A (2.6 eV Å⁻¹) came wholly from it, dx = dx_lo, per formula unit')
add('periodic.iirfix.OH_relief_hi', P['OH_relief_hi_kJ'], 'kJ mol-1', 0, f'{D12S}:periodic.OH_relief_hi_kJ = ½ F dx_hi',
    'the same for dx = 0.09 Å, per formula unit')
add('periodic.iirfix.best_minus_OH_hi', FACTS['periodic.hybrid.e100.best_kJ']['value'] - P['OH_relief_hi_kJ'], 'kJ mol-1', 0,
    'periodic.hybrid.e100.best_kJ − periodic.iirfix.OH_relief_hi', 'PBE0-D3(BJ) cost of the constrained cell less the largest O–H relief estimate')
add('periodic.iirfix.best_minus_OH_hi_kcal', (FACTS['periodic.hybrid.e100.best_kJ']['value'] - P['OH_relief_hi_kJ']) / 4.184,
    'kcal mol-1', 0, 'periodic.iirfix.best_minus_OH_hi / 4.184', 'the same in kcal mol⁻¹')
add('periodic.iirfix.stress_B10', P['stress_B10_kJ'], 'kJ mol-1', 1, f'{D12S}:periodic.stress_B10_kJ',
    'difference (constrained − salt) of the relaxation energies p²V/(2B) at the residual pressures (+7.1, −4.4 kbar), B = 10 GPa, per formula unit')
add('periodic.iirfix.stress_B20', P['stress_B20_kJ'], 'kJ mol-1', 1, f'{D12S}:periodic.stress_B20_kJ', 'the same with B = 20 GPa')
add('periodic.iirfix.best_minus_all', FACTS['periodic.hybrid.e100.best_kJ']['value'] - P['OH_relief_hi_kJ'] - P['stress_B10_kJ'],
    'kJ mol-1', 0, 'periodic.hybrid.e100.best_kJ − periodic.iirfix.OH_relief_hi − periodic.iirfix.stress_B10',
    'PBE0-D3(BJ) cost of the constrained cell less the largest O–H relief and stress estimates, per formula unit')
add('periodic.iirfix.best_minus_all_kcal', (FACTS['periodic.hybrid.e100.best_kJ']['value'] - P['OH_relief_hi_kJ'] - P['stress_B10_kJ']) / 4.184,
    'kcal mol-1', 1, 'periodic.iirfix.best_minus_all / 4.184', 'the same in kcal mol⁻¹')
add('periodic.assumed_B', [10, 20], 'GPa', 0, 'assumption', 'bulk moduli assumed for the stress estimate (organic salt crystals)')
e80, e100 = FACTS['periodic.hybrid.e80.PBE_meV']['value'], FACTS['periodic.hybrid.e100.PBE_meV']['value']
add('periodic.hybrid.PBE_cutoff_change', (e100 - e80) * 96.485332 / 1000 / 2, 'kJ mol-1', 1,
    'periodic.hybrid.e100.PBE_meV − periodic.hybrid.e80.PBE_meV (SG15, 100/400 vs 80/320 Ry)',
    'change of the PBE-D3(BJ) difference E(constrained) − E(salt) between 80 and 100 Ry wavefunction cutoff, per formula unit')

CAVEATS.append('Revision (section 31): the co-crystal torsions are those of the ligand copy used as the docking reference (6RL9 SAN 304, '
               'zinc-bound; the two further SAN copies of 6RL9 lie in secondary sites); the ΔΔ scores and zinc-filtered ranks are orders of '
               'empirical scores, well inside the scoring error of the functions; the O–H relief and stress estimates are order-of-magnitude '
               'checks with stated assumptions, not calculations; the strain against the same-H reference is still rigid (no constrained '
               'relaxation), so it bounds the torsional strain from above.')
