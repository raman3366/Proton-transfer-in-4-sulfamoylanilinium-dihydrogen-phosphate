"""Facts added when the last open audit items were closed (29 Sep 2026, night; section 37 of facts.md), executed by
build_facts.py after facts_referee.py (exec in its namespace: add(), section(), FACTS, J, r, RES, DOCK, DER, CAVEATS,
MISSING and DISAGREE are the ones defined there).

  37. (a) Measured binding of acetazolamide and sulfanilamide to hCA II from PLBD (snapshot
          paper/audit/research/plbd_hca2_azm_sfm.json; notes paper/audit/research/plbd.md), re-derived in
          derive/d14_completion.py from the experiment records and the intrinsic binding constants.
      (b) The AutoDock4Zn orientation gaps with one definition of "near-native" (RMSD <= 1.0 A to the crystal ligand):
          the 6RL9 validation run and the eight AutoDock4Zn runs of sulfanilamide, from every mode of every seed.
      (c) The sulfonamide torsion of every docked top pose measured on the raw pose coordinates (Table S31), with the
          relaxed DFT strain at those torsions; checked against the rebuilt-pose values of production_summary.json and
          the strain ranges of section 35.
"""
D14 = J(DER / 'd14_completion.json')
D14S = 'derived/d14_completion.json'
PL_ = D14['plbd']
PLS_ = 'paper/audit/research/plbd_hca2_azm_sfm.json'

section('37. Closing items: measured hCA II binding of the controls, one "near-native" definition, raw-pose torsions',
        'Added on 29 Sep 2026 (night) to close the last audit items; no new calculation. (a) PLBD (Lingė et al. 2023), the curated '
        'database of the Matulis laboratory, retrieved through its REST interface on 29 Sep 2026 (rows kept in '
        f'{PLS_}): observed dissociation constants of hCA II from fluorescent thermal shift assays (table kd_obs_tsa_all_cas, '
        'column ca02, 37 °C; the experiments are those of Linkuvienė et al. 2018 at pH 7.0), and the intrinsic constants, from '
        'which PLBD has subtracted the protonation reactions of ligand and protein linked to binding (table intrinsic_parameters). '
        '(b) A docked mode is near-native when its symmetry-corrected heavy-atom RMSD to the crystal ligand is at most 1.0 Å '
        '(receptor frame); gap = score of the best-ranked near-native mode − score of the top mode (0 when the top mode is '
        'near-native); separation = |best near-native score − best score beyond 2.0 Å|. (c) τ_fold measured on the top pose of '
        'each seed as written by the docking program (out.pdbqt), heavy atoms only; strain from the spline of section 35.')

for nm_, lab_ in (('AZM', 'acetazolamide'), ('SFM', 'sulfanilamide')):
    q_ = PL_[nm_]
    k_ = nm_.lower()
    add(f'lit.plbd.{k_}.kd_obs', q_['kd_obs_37_M'] * 1e9 if nm_ == 'AZM' else q_['kd_obs_37_M'] * 1e6,
        'nM' if nm_ == 'AZM' else 'µM', 0,
        f'{PLS_}: kd_obs_tsa_all_cas["{"1783" if nm_ == "AZM" else "2255"}-37"].ca02 (PLBD, retrieved 29 Sep 2026); primary data '
        'linkuviene2018thermodynamic (PLBD publication Linkuviene2018b)',
        f'{lab_}: observed dissociation constant of hCA II, FTSA, pH 7.0, 37 °C (mean over {q_["n_exp_37"]} experiments)',
        check=(q_['kd_obs_37_rederived_M'] * (1e9 if nm_ == 'AZM' else 1e6)), tol=0.005,
        note=f'Re-derived as the mean 1/K_b,obs of the 37 °C, pH 7.0 experiments {", ".join(e["label"] for e in q_["exp_37"])} '
             f'({D14S}:plbd.{nm_}). With every 37 °C record, including those without a primary publication: '
             f'{q_["kd_obs_37_all_exp_M"] * (1e9 if nm_ == "AZM" else 1e6):.2f} {"nM" if nm_ == "AZM" else "µM"}. '
             f'At 25 °C PLBD gives {q_["kd_obs_25_M"] * (1e9 if nm_ == "AZM" else 1e6):.2f} {"nM" if nm_ == "AZM" else "µM"}.')
    add(f'lit.plbd.{k_}.kd_int', q_['kd_int_37_M'] * 1e9, 'nM', 1,
        f'{PLS_}: intrinsic_parameters (compound {"1783" if nm_ == "AZM" else "2255"}, protein 2, 37 °C).kb_int, K_d,int = 1/K_b,int',
        f'{lab_}: intrinsic dissociation constant of hCA II at 37 °C (binding-linked protonation of ligand and protein subtracted)',
        check=1e9 / q_['kb_int_37'], tol=1e-6,
        note=f'n = {q_["n_int_37"]} experiments (FTSA {q_["int_37_tsa_ids"]}; ITC {q_["int_37_itc_ids"] or "none"}); FTSA-only '
             f'intrinsic value (kd_int_tsa_all_cas) {q_["kd_int_tsa_37_M"] * 1e9:.2f} nM.')
    add(f'lit.plbd.{k_}.dG_int', q_['dG_int_37_kJ'], 'kJ mol-1', 2, f'{PLS_}: intrinsic_parameters (37 °C).dbg_int',
        f'{lab_}: intrinsic Gibbs energy of binding to hCA II at 37 °C', check=q_['dG_int_37_rederived_kJ'], tol=0.01)
add('lit.plbd.T', 37, '°C', 0, f'{PLS_}: experiment_tsa temperature', 'temperature of the PLBD hCA II constants quoted in the text')
add('lit.plbd.pH', 7.0, None, 1, f'{PLS_}: experiment_tsa ph', 'pH of the FTSA experiments behind the observed constants')
add('lit.plbd.ratio_obs', PL_['ratio_obs_37'], None, 0, f'{D14S}:plbd.ratio_obs_37',
    'observed K_d(sulfanilamide) / K_d(acetazolamide), hCA II, 37 °C')
assert abs(PL_['ratio_int_37'] - 1) < 0.01

VAL_ = D14['orientation_validation']
add('dock.val.6RL9.ad4zn.orientation_sep', [VAL_['sep_min'], VAL_['sep_max']], 'kcal mol-1', 3,
    f'{D14S}:orientation_validation (T3_6RL9/runs/sfa_neutral__ad4zn, every mode of the 10 seeds)',
    '6RL9 AutoDock4Zn validation run (ETKDG input): score separation of the near-native (≤ 1.0 Å) and ring-rotated (> 2.0 Å) '
    'orientations, min and max over the seeds',
    note=f'Every seed holds both orientations ({VAL_["n_with_both"]} of {VAL_["n_seeds"]}). RMSDs agree with '
         f'dock/analysis/summary.json within {VAL_["check_rmsd_vs_summary"]:.4f} Å. Replaces the 0.000–0.046 of dock/VALIDATION.md, '
         'which used the lowest-RMSD mode (numbers audit round 2, m5).')
add('dock.val.6RL9.ad4zn.orientation_gap', [VAL_['gap_min'], VAL_['gap_max']], 'kcal mol-1', 3, f'{D14S}:orientation_validation',
    '6RL9 AutoDock4Zn validation run: score of the best-ranked near-native mode minus that of the top mode, min and max over the seeds')
E8_ = D14['orientation_eight_all']
add('dock.ad4zn8.orientation_gap', [E8_['gap_min'], E8_['gap_max']], 'kcal mol-1', 3, f'{D14S}:orientation_eight_all',
    'eight AutoDock4Zn runs of sulfanilamide (6RL9 and 3HS4; production and input-geometry diagnostics): best-ranked near-native '
    'mode behind the top mode, min and max over all seeds',
    check=float(FACTS['dock.prod.ad4zn_orientation_gap']['value'].split('–')[1]), tol=0.0005,
    note=f'Near-native mode present in every seed: {E8_["all_seeds_near"]}; its worst rank {E8_["worst_rank"]}.')

TR_ = D14['tau_raw']
for job_, q_ in TR_.items():
    tg_, rest_ = job_.split('/')
    sp_, fn_ = rest_.split('__')[0], rest_.split('__')[1]
    k_ = f'dock.tors.raw.{tg_.split("_")[1]}.{sp_}.{fn_}'
    addbulk(k_, [q_['tau_min'], q_['tau_max']], 'deg', 1, f'{D14S}:tau_raw["{job_}"]',
            f'{tg_.split("_")[1]} {sp_} {fn_}: τ_fold of the {q_["n"]} docked top poses, raw pose coordinates (min, max)',
            'dock.tors.raw', 'τ_fold range of the docked top poses on the raw pose coordinates (Table S31)')
    FACTS[k_]['note'] = f'rebuilt-pose value of production_summary.json within {q_["max_diff_vs_rebuilt"]:.2f}°'
    if 'strain_min' in q_:
        ks_ = f'tors.seeds.{sp_}.{tg_.split("_")[1]}.{fn_}'
        for i_, x_ in enumerate((q_['strain_min'], q_['strain_max'])):
            if round(x_ + 0.0, 2) != round(FACTS[ks_]['value'][i_] + 0.0, 2) and abs(x_ - FACTS[ks_]['value'][i_]) > 0.005:
                DISAGREE.append(f'{ks_}[{i_}]: raw-pose strain {x_:.4f} vs section 35 {FACTS[ks_]["value"][i_]:.4f}')
        addbulk(k_ + '.strain', [q_['strain_min'], q_['strain_max']], 'kcal mol-1', 2, f'{D14S}:tau_raw["{job_}"].strain_min/max',
                f'{tg_.split("_")[1]} {sp_} {fn_}: relaxed DFT strain at the raw-pose τ_fold (min, max)',
                'dock.tors.raw', 'τ_fold range of the docked top poses on the raw pose coordinates (Table S31)')
RA_ = D14['tau_raw_all']
add('dock.tors.raw.all', [RA_['tau_min'], RA_['tau_max']], 'deg', 1, f'{D14S}:tau_raw_all',
    f'τ_fold over all {RA_["n_poses"]} docked top poses of the sulfanilamide species (raw pose coordinates)',
    note=f'largest difference from the rebuilt-pose values {RA_["max_diff_vs_rebuilt"]:.2f}°')
for i_, kk_ in enumerate(('dock.tors.range_min', 'dock.tors.range_max')):
    if f'{RA_["tau_min" if i_ == 0 else "tau_max"]:.1f}' != f'{FACTS[kk_]["value"]:.1f}':
        DISAGREE.append(f'{kk_}: raw-pose extreme prints differently from the rebuilt-pose value')
