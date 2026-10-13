"""Facts of the control receptor (4 Oct 2026; section 38 of facts.md), executed by build_facts.py after
facts_completion.py (exec in its namespace: add(), addbulk(), section(), FACTS, J, r, DOCK, DER, CAVEATS, MISSING
and DISAGREE are the ones defined there).

  38. The ATP site of Escherichia coli DNA gyrase B (PDB 4DUH), an enzyme that is not a target of sulfanilamide,
      docked with the protocol of the five receptors: redocking of the crystal inhibitor, the scores of the species
      of the production matrix and of ciprofloxacin, the receptor variants (second loop build, built residues removed,
      ordered crystal water put back), the clorobiocin pose of 1KZN (which the receptor cannot hold, so that
      clorobiocin is not docked) and the whole-fragment search.  Source:
      dock/control/control_summary.json (dock/control_gyrb.py); re-derived from the raw out.pdbqt files in
      derive/d15_control.py.
"""
CS_ = J(DOCK / 'control' / 'control_summary.json')
D15 = J(DER / 'd15_control.json')
CSS_ = 'dock/control/control_summary.json'
D15S = 'derived/d15_control.json'
C1_, C1B_, C1W_, C1N_ = 'C1_4DUH', 'C1_4DUH_seed12345', 'C1_4DUH_wat', 'C1_4DUH_noloop'
INV_ = J(DOCK / 'targets_inventory.json')['C1_4DUH']
NAME_ = {'sfa_neutral': 'sulfanilamide', 'sfa_cation': '4-sulfamoylanilinium', 'sfa_anion': 'sulfonamidate',
         'smx_neutral': 'sulfamethoxazole', 'smx_anion': 'sulfamethoxazole anion', 'paba_anion': 'pABA⁻',
         'azm_neutral': 'acetazolamide', 'azm_anion': 'acetazolamide anion', 'rli_anion': '4DUH inhibitor',
         'cfx_zwit': 'ciprofloxacin'}
RV_ = {C1_: 'primary receptor', C1B_: 'second loop build', C1W_: 'ordered crystal water kept',
       C1N_: 'built residues removed'}
TGT_ = {'3TZF': 'T1_3TZF', '7S2K': 'T5_7S2K', '6RL9': 'T3_6RL9', '3HS4': 'T4_3HS4'}


def _ms(vals_):
    """Mean and SD (ddof 1) at full precision: the summary stores them rounded to three decimals, and a value
    rounded to a tie would break the rounding rule (half to even on the full-precision value)."""
    return _stt.mean(vals_), (_stt.stdev(vals_) if len(vals_) > 1 else 0.0)


def _score(job_):
    return _ms(CS_['jobs'][job_]['score']['per_seed'])

section('38. Control receptor: the ATP site of E. coli DNA gyrase B (4DUH), not a target of sulfanilamide',
        'Docked with the protocol of the five '
        'receptors of Section 2.5 (AutoDock Vina 1.2.7, Vina and Vinardo scoring, 20 Å cube on the centroid of the crystal ligand, '
        'exhaustiveness 32, 10 seeds; sulfanilamide species from their DFT geometries, every other ligand from one '
        'ETKDG/MMFF94s conformer per seed). The runs are kept out of gate G2 and of the production summary: the '
        f'directories are named C1_*, and every number below comes from {CSS_} (dock/control_gyrb.py analyse). '
        'A pose is "in the cavity of the inhibitor" when its heavy-atom centroid lies within 4.0 Å of a heavy atom of '
        'the crystal inhibitor (this says where a pose is, not that it is right); a crystal water is "overlapped" when '
        'its oxygen lies within 2.4 Å of a heavy atom of the pose. Scores are empirical docking scores, never '
        'affinities. No measurement says that sulfanilamide does not bind this pocket: the control is a receptor that '
        'is not a target, and the texts say so. '
        f"Our parse of every out.pdbqt reproduces all {D15['n_jobs_compared']} job means to "
        f"{D15['max_abs_score_mean_diff_vs_summary']:.4f} kcal/mol and the redocking RMSDs to "
        f"{D15['max_abs_rmsd_mean_diff_vs_summary']:.4f} Å ({D15S}).")

# ---------------------------------------------------------------------------------------------- receptor and boxes
add('ctrl.pdb', '4DUH', None, None, 'dock/targets.py:C1_4DUH; dock/targets_inventory.json', 'control receptor: PDB entry (E. coli DNA gyrase B, 24 kDa N-terminal fragment, chain A)')
add('ctrl.resolution', INV_['resolution'], 'A', 2, 'dock/targets_inventory.json:C1_4DUH.resolution', 'resolution of 4DUH')
add('ctrl.fragment_kda', 24, 'kDa', 0, 'RCSB entry 4DUH, title ("24 kDa domain of E. coli DNA gyrase B"); brvar2012structure', 'size of the N-terminal gyrase B fragment crystallised in 4DUH and 1KZN')
add('ctrl.ligand_heavy', INV_['ligand']['heavy'], None, 0, 'dock/targets_inventory.json:C1_4DUH.ligand.heavy', 'heavy atoms of the crystal inhibitor of 4DUH (component RLI)')
add('ctrl.loop.first', 109, None, 0, 'dock/targets_inventory.json:C1_4DUH.gaps.unobserved', 'first unobserved residue of the ATP lid in chain A of 4DUH')
add('ctrl.loop.last', 116, None, 0, 'dock/targets_inventory.json:C1_4DUH.gaps.unobserved', 'last unobserved residue of the ATP lid in chain A of 4DUH')
add('ctrl.loop.n_residues', 8, None, 0, 'dock/C1_4DUH/prep_report.json:missing_residues_rebuilt_auth_numbering', 'residues of the ATP lid rebuilt by PDBFixer (109–116)')
LS_ = CS_['loop_sensitivity']
for t_, tag_ in ((C1_, 'primary'), (C1B_, 'second')):
    q_ = LS_['receptors'][t_]
    add(f'ctrl.loop.{tag_}.to_inhibitor', q_['loop_to_crystal_inhibitor_min_A'], 'A', 1, f'{CSS_}:loop_sensitivity.receptors["{t_}"].loop_to_crystal_inhibitor_min_A',
        f'{RV_[t_]}: shortest distance from a heavy atom of the rebuilt residues to the crystal inhibitor', check=D15['loop'][t_]['to_inhibitor'], tol=0.006)
    add(f'ctrl.loop.{tag_}.atoms_in_box', q_['loop_heavy_atoms_in_box'], None, 0, f'{CSS_}:loop_sensitivity.receptors["{t_}"].loop_heavy_atoms_in_box',
        f'{RV_[t_]}: heavy atoms of the rebuilt residues inside the 20 Å box', check=D15['loop'][t_]['in_box'], tol=0.01)
    add(f'ctrl.prep.{tag_}.max_disp', q_['restrained_atoms_max_displacement_A'], 'A', 2, f'dock/{t_}/prep_report.json:minimisation.max_displacement_restrained_A',
        f'{RV_[t_]}: largest displacement of a deposited atom in the restrained minimisation')
add('ctrl.loop.build_rms', LS_['build_difference']['loop_heavy_atoms_rms_A'], 'A', 1, f'{CSS_}:loop_sensitivity.build_difference.loop_heavy_atoms_rms_A',
    'RMS difference of the rebuilt residues between the two loop builds (heavy atoms)', check=D15['loop']['build_difference']['loop_rms'], tol=0.006)
add('ctrl.loop.n_built_atoms', LS_['receptors'][C1N_]['n_atoms_removed'], None, 0, f'dock/{C1N_}/prep_report.json:n_atoms_removed', 'atoms of the built residues 109–116 removed in the loop-free receptor (hydrogens included)')
WB_ = CS_['protocol']['whole_fragment_box']
add('ctrl.whole.box', WB_['size'], 'A', 1, f'dock/{C1_}/box_whole.json:size', 'edges of the whole-fragment search box (receptor extent plus 5 Å on every side)')
add('ctrl.whole.volume_ratio', WB_['volume_ratio'], None, 0, f'dock/{C1_}/box_whole.json:volume_ratio', 'volume of the whole-fragment box over that of the 20 Å cube')
add('ctrl.whole.exhaustiveness', WB_['exhaustiveness'], None, 0, f'dock/{C1_}/box_whole.json:exhaustiveness', 'exhaustiveness of the whole-fragment search (32 × volume ratio)')
add('ctrl.site_cut', 4.0, 'A', 1, 'dock/control_gyrb.py:SITE_CUT', 'a pose is in the inhibitor site when its centroid lies within this distance of a heavy atom of the crystal inhibitor')
add('ctrl.water_clash', 2.4, 'A', 1, 'dock/control_gyrb.py:WATER_CLASH', 'a crystal water is overlapped when its oxygen lies within this distance of a heavy atom of a pose')
WD_ = CS_['water_diagnostic']
add('ctrl.water.n_chain', WD_['n_chain_waters'], None, 0, f'{CSS_}:water_diagnostic.n_chain_waters', 'waters assigned to chain A in the deposited 4DUH file')
add('ctrl.water.n_kept', WD_['receptor']['n_waters_kept'], None, 0, f'dock/{C1W_}/prep_report.json:n_waters_kept', 'crystal waters added to the primary receptor in the water diagnostic')
add('ctrl.water.n_in_box', WD_['receptor']['n_waters_in_box'], None, 0, f'dock/{C1W_}/prep_report.json:n_waters_in_box', 'of these, waters inside the 20 Å box')
add('ctrl.water.nearest_to_inhibitor', WD_['crystal_inhibitor']['nearest_water_A'], 'A', 2, f'{CSS_}:water_diagnostic.crystal_inhibitor.nearest_water_A', 'shortest distance from a crystal water to the crystal inhibitor')

# ---------------------------------------------------------------------------------------------- crystal contacts
for v_ in CS_['validation'].values():
    if v_['ligand'] == 'rli_anion' and v_['job'].startswith(C1_ + '/'):
        for c_, k_ in zip(v_['site_checks_on_crystal_reference'], ('asp73', 'arg136', 'gly101')):
            add(f'ctrl.xtal.{k_}', c_['distance'], 'A', 2, f'{CSS_}:validation["{v_["job"]}"].site_checks_on_crystal_reference',
                f'crystal inhibitor in 4DUH: {c_["check"]}')
        break

# ---------------------------------------------------------------------------------------------- redocking
NN_ = WD_['near_native_mode']
for job_, v_ in CS_['validation'].items():
    if v_['ligand'] != 'rli_anion':
        continue
    t_, tag_ = job_.split('/runs/')
    sc_ = v_['scoring']
    b_ = f"ctrl.val.{ {C1_: 'primary', C1B_: 'second', C1W_: 'water', C1N_: 'noloop'}[t_] }.{sc_}"
    st_, g_, d_ = v_['stats'], v_['gate'], D15['rmsd'].get(job_)
    lab_ = f'redock of the 4DUH inhibitor, {RV_[t_]}, {sc_}'
    rm_, rs_ = _ms([x_['rmsd_top'] for x_ in v_['per_seed']])
    sm_, ss_ = _score(job_)
    add(b_ + '.rmsd', rm_, 'A', 2, f'{CSS_}:validation["{job_}"].per_seed[].rmsd_top (mean)', f'{lab_}: top-pose RMSD, mean over seeds', check=d_['mean'] if d_ else None, tol=0.0006)
    add(b_ + '.rmsd_sd', rs_, 'A', 2, f'{CSS_}:validation["{job_}"].per_seed[].rmsd_top (SD)', f'{lab_}: RMSD SD over seeds', check=d_['sd'] if d_ else None, tol=0.0006)
    add(b_ + '.f1', st_['success_fraction_1A'], None, 1, f'{CSS_}:validation["{job_}"].stats.success_fraction_1A', f'{lab_}: fraction of seeds with the top pose within 1 Å', check=d_['f_le_1'] if d_ else None, tol=0.001)
    add(b_ + '.f2', st_['success_fraction_2A'], None, 1, f'{CSS_}:validation["{job_}"].stats.success_fraction_2A', f'{lab_}: fraction within 2 Å', check=d_['f_le_2'] if d_ else None, tol=0.001)
    add(b_ + '.score', sm_, 'kcal mol-1', 1, f'{CSS_}:jobs["{job_}"].score.per_seed (mean)', f'{lab_}: mean top-pose score', check=D15['scores'][job_]['mean'], tol=1e-6)
    add(b_ + '.score_sd', ss_, 'kcal mol-1', 2, f'{CSS_}:jobs["{job_}"].score.per_seed (SD)', f'{lab_}: SD of the top-pose score', check=D15['scores'][job_]['sd'], tol=1e-6)
    add(b_ + '.pb', st_['posebusters_physics_pass_fraction'], None, 1, f'{CSS_}:validation["{job_}"].stats.posebusters_physics_pass_fraction', f'{lab_}: PoseBusters physics pass fraction')
    add(b_ + '.gate', 'pass' if g_['passes_strict_1A'] else 'fail', None, None, f'{CSS_}:validation["{job_}"].gate.passes_strict_1A', f'{lab_}: the 1.0 Å gate')
    add(b_ + '.crit2A', 'pass' if g_['passes_2A_criterion'] else 'fail', None, None, f'{CSS_}:validation["{job_}"].gate.passes_2A_criterion', f'{lab_}: the conventional 2.0 Å criterion')
    for chk_, q_ in st_['site_checks'].items():
        k_ = 'asp73' if 'Asp73' in chk_ else 'arg136' if 'Arg136' in chk_ else 'gly101'
        add(f'{b_}.{k_}', q_['distance_mean'], 'A', 2, f'{CSS_}:validation["{job_}"].stats.site_checks["{chk_}"]', f'{lab_}: mean {chk_} of the top poses (pass fraction {q_["pass_fraction"]})')
    n_ = NN_[job_]
    add(b_ + '.nn.n_seeds', n_['seeds_with_a_mode_within_1A'], None, 0, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].seeds_with_a_mode_within_1A', f'{lab_}: seeds that hold a mode within 1.0 Å of the crystal pose')
    if n_['seeds_with_a_mode_within_1A']:
        add(b_ + '.nn.rank', [n_['rank_min'], n_['rank_max']], None, 0, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].rank_min/max', f'{lab_}: rank of the first mode within 1.0 Å (min, max over seeds)')
        add(b_ + '.nn.rmsd', [n_['rmsd_min_A'], n_['rmsd_max_A']], 'A', 2, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].rmsd_min_A/max_A', f'{lab_}: RMSD of that mode (min, max over seeds)')
        add(b_ + '.nn.gap', [n_['gap_min'], n_['gap_max']], 'kcal mol-1', 1, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].gap_min/max', f'{lab_}: score of that mode minus the top score (min, max over seeds)')
        add(b_ + '.nn.gap_mean', _ms([x_['gap_near_native_minus_top'] for x_ in n_['per_seed'] if 'rank_near_native' in x_])[0], 'kcal mol-1', 2,
            f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].per_seed[].gap_near_native_minus_top (mean)', f'{lab_}: the same, mean over seeds')
        add(b_ + '.nn.score', _ms([x_['score_near_native'] for x_ in n_['per_seed'] if 'rank_near_native' in x_])[0], 'kcal mol-1', 1,
            f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].per_seed[].score_near_native (mean)', f'{lab_}: mean score of the first mode within 1.0 Å')
        add(b_ + '.nn.waters', n_['waters_overlapped_near_native'][1], None, 0, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].waters_overlapped_near_native', f'{lab_}: crystal waters overlapped by that mode (maximum over seeds)')
    add(b_ + '.top_waters', n_['waters_overlapped_top'], None, 0, f'{CSS_}:water_diagnostic.near_native_mode["{job_}"].waters_overlapped_top', f'{lab_}: crystal waters overlapped by the top pose (min, max over seeds)')

# ---------------------------------------------------------------------------------------------- clorobiocin
CB_ = CS_['clorobiocin_reference']
OV_ = CB_['overlap'][C1_]
add('ctrl.cbn.heavy', CB_['heavy_atoms'], None, 0, f'{CSS_}:clorobiocin_reference.heavy_atoms', 'heavy atoms of clorobiocin', check=D15['clorobiocin']['n_heavy'], tol=0.01)
add('ctrl.cbn.n_ca', CB_['n_ca_atoms'], None, 0, f'{CSS_}:clorobiocin_reference.n_ca_atoms', 'Cα atoms paired by sequence alignment in the superposition of 1KZN on 4DUH (chain A)',
    note=f"Our own alignment (difflib) pairs {D15['clorobiocin']['n_ca_pairs']} Cα atoms and gives {D15['clorobiocin']['ca_rmsd']:.3f} Å; the ligand it carries over lies {D15['clorobiocin']['own_vs_file_rmsd']:.3f} Å (RMSD) from the reference file.")
add('ctrl.cbn.ca_rmsd', CB_['superposition_ca_rmsd_A'], 'A', 2, f'{CSS_}:clorobiocin_reference.superposition_ca_rmsd_A', 'Cα RMSD of that superposition')
add('ctrl.cbn.min_dist', OV_['min_heavy_heavy_distance_A'], 'A', 2, f'{CSS_}:clorobiocin_reference.overlap["{C1_}"].min_heavy_heavy_distance_A',
    f"shortest heavy-atom distance between the clorobiocin pose of 1KZN carried into 4DUH and the prepared receptor ({OV_['closest_receptor_atom']})", check=D15['clorobiocin']['min_heavy_heavy'], tol=0.006)
add('ctrl.cbn.n_pairs_2p5', OV_['n_pairs_below_2p5A'], None, 0, f'{CSS_}:clorobiocin_reference.overlap["{C1_}"].n_pairs_below_2p5A', 'heavy-atom pairs of that pose and the receptor closer than 2.5 Å',
    check=D15['clorobiocin']['n_pairs_below_2p5'], tol=0.01)
add('ctrl.cbn.pair_cut', 2.5, 'A', 1, 'dock/control_gyrb.py:cross_reference (n_pairs_below_2p5A)', 'distance below which a heavy-atom pair of the carried-over clorobiocin pose and the receptor is counted as an overlap')
add('ctrl.cbn.n_pairs_unpaired', OV_['n_pairs_below_2p5A_without_counterpart_in_1KZN'], None, 0, f'{CSS_}:clorobiocin_reference.overlap["{C1_}"].n_pairs_below_2p5A_without_counterpart_in_1KZN',
    'of these pairs, those with a residue of 4DUH that has no counterpart in the 1KZN model (' + ', '.join(OV_['residues_below_2p5A_without_counterpart_in_1KZN']) + ')',
    check=D15['clorobiocin']['n_pairs_below_2p5_without_counterpart'], tol=0.01)
add('ctrl.cbn.unpaired_residues', ' and '.join(x_[:1] + x_[1:3].lower() + x_[3:] for x_ in OV_['residues_below_2p5A_without_counterpart_in_1KZN']), None, None,
    f'{CSS_}:clorobiocin_reference.overlap["{C1_}"].residues_below_2p5A_without_counterpart_in_1KZN', 'the residues of those pairs: ordered in 4DUH, not modelled in 1KZN')
assert D15['clorobiocin']['residues_below_2p5_without_counterpart'] == OV_['residues_below_2p5A_without_counterpart_in_1KZN']
assert D15['clorobiocin']['residues_of_4DUH_without_counterpart_in_1KZN'] == CB_['residues_of_4DUH_without_counterpart_in_1KZN']
add('ctrl.cbn.n_res_unpaired', len(CB_['residues_of_4DUH_without_counterpart_in_1KZN']), None, 0, f'{CSS_}:clorobiocin_reference.residues_of_4DUH_without_counterpart_in_1KZN',
    'observed residues of 4DUH chain A without a counterpart in the 1KZN model (' + ', '.join(str(x_) for x_ in CB_['residues_of_4DUH_without_counterpart_in_1KZN']) + ')')

# ---------------------------------------------------------------------------------------------- scores of every species
# Three sets: the receptor as prepared with Vina and with Vinardo (keys ctrl.score.*), and the receptor with the ordered
# crystal water kept, the one that redocks the inhibitor, with Vina (keys ctrl.wat.*).
WHAT_ = 'control receptor 4DUH: scores, ranks, sites and contacts of every docked species (Tables S33 and S34)'
TABS_ = (('vina', C1_, 'vina', 'ctrl.score', '4DUH'), ('vinardo', C1_, 'vinardo', 'ctrl.score', '4DUH'),
         ('vina_water', C1W_, 'vina', 'ctrl.wat', '4DUH with the crystal water'))
for tkey_, rec_, sc_, pre_, rlab_ in TABS_:
    tb_ = CS_['score_table'].get(tkey_)
    if not tb_:
        continue
    for row_ in tb_['ranking']:
        sp_ = row_['species']
        job_ = next(k for k, q in CS_['jobs'].items() if q['target'] == rec_ and q['ligand'] == sp_ and q['scoring'] == sc_)
        q_ = CS_['jobs'][job_]
        b_ = f'{pre_}.{sp_}.{sc_}'
        lab_ = f'{rlab_}, {NAME_[sp_]}, {sc_}'
        sm_, ss_ = _score(job_)
        add(b_, sm_, 'kcal mol-1', 1, f'{CSS_}:jobs["{job_}"].score.per_seed (mean)', f'{lab_}: mean top-pose score over {q_["n_seeds"]} seeds', check=D15['scores'][job_]['mean'], tol=1e-6)
        add(b_ + '.sd', ss_, 'kcal mol-1', 2, f'{CSS_}:jobs["{job_}"].score.per_seed (SD)', f'{lab_}: SD of the top-pose score', check=D15['scores'][job_]['sd'], tol=1e-6)
        addbulk(b_ + '.n_seeds', q_['n_seeds'], None, 0, f'{CSS_}:jobs["{job_}"].n_seeds', f'{lab_}: seeds', pre_, WHAT_)
        addbulk(b_ + '.rank', row_['rank'], None, 0, f'{CSS_}:score_table.{tkey_}.ranking', f'{lab_}: rank by mean score among the {tb_["n_species"]} species docked into this receptor', pre_, WHAT_)
        addbulk(b_ + '.best', q_['score']['min'], 'kcal mol-1', 1, f'{CSS_}:jobs["{job_}"].score.min', f'{lab_}: best top-pose score of a seed', pre_, WHAT_)
        addbulk(b_ + '.worst', q_['score']['max'], 'kcal mol-1', 1, f'{CSS_}:jobs["{job_}"].score.max', f'{lab_}: weakest top-pose score of a seed', pre_, WHAT_)
        addbulk(b_ + '.heavy', q_['heavy_atoms'], None, 0, f'{CSS_}:jobs["{job_}"].heavy_atoms', f'{NAME_[sp_]}: heavy atoms', pre_, WHAT_, check=D15['scores'][job_]['heavy_atoms'], tol=0.01)
        addbulk(b_ + '.LE', -sm_ / q_['heavy_atoms'], 'kcal mol-1', 2, f'{CSS_}:jobs["{job_}"].score.per_seed, heavy_atoms', f'{lab_}: ligand efficiency (−score per heavy atom)', pre_, WHAT_)
        addbulk(b_ + '.n_clusters', q_['clusters']['n_clusters'], None, 0, f'{CSS_}:jobs["{job_}"].clusters.n_clusters', f'{lab_}: clusters of the top poses (2.0 Å)', pre_, WHAT_)
        addbulk(b_ + '.cl1', q_['clusters']['clusters'][0]['population'], None, 0, f'{CSS_}:jobs["{job_}"].clusters.clusters[0].population', f'{lab_}: seeds in the best-scoring cluster', pre_, WHAT_)
        addbulk(b_ + '.in_site', round(q_['site']['in_site_fraction'] * q_['n_seeds']), None, 0, f'{CSS_}:jobs["{job_}"].site.in_site_fraction', f'{lab_}: seeds whose top pose lies in the inhibitor site', pre_, WHAT_,
                check=D15['site'][job_]['in_site_fraction'] * q_['n_seeds'], tol=0.01)
        addbulk(b_ + '.site_dist', q_['site']['centroid_to_inhibitor_mean_A'], 'A', 1, f'{CSS_}:jobs["{job_}"].site.centroid_to_inhibitor_mean_A', f'{lab_}: mean distance of the pose centroid from the nearest heavy atom of the crystal inhibitor', pre_, WHAT_)
        addbulk(b_ + '.site_range', q_['site']['centroid_to_inhibitor_range_A'], 'A', 1, f'{CSS_}:jobs["{job_}"].site.centroid_to_inhibitor_range_A', f'{lab_}: the same distance, smallest and largest over the seeds', pre_, WHAT_,
                check=[D15['site'][job_]['centroid_min'], D15['site'][job_]['centroid_max']], tol=0.006)
        addbulk(b_ + '.asp73', q_['site']['asp73_distance_mean_A'], 'A', 1, f'{CSS_}:jobs["{job_}"].site.asp73_distance_mean_A', f'{lab_}: mean shortest distance from a ligand N or O to the Asp73 carboxylate oxygens', pre_, WHAT_,
                check=D15['site'][job_]['asp73_mean'], tol=0.006)
        addbulk(b_ + '.asp73_seeds', round(q_['site']['asp73_contact_fraction_3p5A'] * q_['n_seeds']), None, 0, f'{CSS_}:jobs["{job_}"].site.asp73_contact_fraction_3p5A', f'{lab_}: seeds whose top pose has a ligand N or O within 3.5 Å of the Asp73 carboxylate', pre_, WHAT_,
                check=D15['site'][job_]['asp73_contact_fraction'] * q_['n_seeds'], tol=0.01)
        addbulk(b_ + '.waters', q_['site']['crystal_waters_overlapped_mean'], None, 1, f'{CSS_}:jobs["{job_}"].site.crystal_waters_overlapped_mean', f'{lab_}: crystal waters overlapped by the top pose, mean over seeds', pre_, WHAT_)
        addbulk(b_ + '.loop', q_['site']['loop_min_distance_A'], 'A', 1, f'{CSS_}:jobs["{job_}"].site.loop_min_distance_A', f'{lab_}: shortest distance from a top pose to a heavy atom of the rebuilt residues', pre_, WHAT_,
                check=D15['site'][job_]['loop_min'], tol=0.006)
        if q_['site'].get('loop_distance_range_A'):
            addbulk(b_ + '.loop_range', q_['site']['loop_distance_range_A'], 'A', 1, f'{CSS_}:jobs["{job_}"].site.loop_distance_range_A', f'{lab_}: closest distance from the top pose to a heavy atom of the rebuilt residues, smallest and largest over the seeds', pre_, WHAT_,
                    check=[D15['site'][job_]['loop_min'], D15['site'][job_]['loop_max']], tol=0.006)
        addbulk(b_ + '.waters_range', q_['site']['crystal_waters_overlapped_range'], None, 0, f'{CSS_}:jobs["{job_}"].site.crystal_waters_overlapped_range', f'{lab_}: crystal waters overlapped by the top pose, smallest and largest over the seeds', pre_, WHAT_)
        addbulk(b_ + '.loop_seeds', round(q_['site']['loop_contact_fraction_4A'] * q_['n_seeds']), None, 0, f'{CSS_}:jobs["{job_}"].site.loop_contact_fraction_4A', f'{lab_}: seeds whose top pose lies within 4.0 Å of a rebuilt residue', pre_, WHAT_)
        addbulk(b_ + '.pb', q_['posebusters']['pass_fraction'], None, 1, f'{CSS_}:jobs["{job_}"].posebusters.pass_fraction', f'{lab_}: PoseBusters pass fraction of the top poses (dock configuration)', pre_, WHAT_)
    sr_ = tb_['size_relation']
    if sr_.get('r_score_vs_HA') is not None:
        add(f'ctrl.size.{tkey_}.r', sr_['r_score_vs_HA'], None, 2, f'{CSS_}:score_table.{tkey_}.size_relation.r_score_vs_HA', f'{rlab_}, {sc_}: Pearson r between the mean score and the number of heavy atoms over the {sr_["n"]} species',
            check=D15['size'][tkey_]['r_score_vs_HA'], tol=0.0006)
    if sr_.get('slope_kcal_per_heavy_atom') is not None:
        add(f'ctrl.size.{tkey_}.slope', sr_['slope_kcal_per_heavy_atom'], 'kcal mol-1', 2, f'{CSS_}:score_table.{tkey_}.size_relation.slope_kcal_per_heavy_atom',
            f'{rlab_}, {sc_}: slope of the least-squares line of the mean score against the number of heavy atoms (per heavy atom)', check=D15['size'][tkey_]['slope'], tol=0.0006)
    add(f'ctrl.size.{tkey_}.n', sr_['n'], None, 0, f'{CSS_}:score_table.{tkey_}.size_relation.n', f'{rlab_}, {sc_}: species docked into this receptor')
    zc_ = CS_['comparison_with_targets'].get(tkey_, {}).get('size_relation_common_species')
    if zc_:
        add(f'ctrl.size8.{tkey_}.n', zc_['n_species'], None, 0, f'{CSS_}:comparison_with_targets.{tkey_}.size_relation_common_species.n_species', f'{sc_}: species docked into {rlab_} and into all four production structures')
        for rk_, rv_ in zc_['r_score_vs_HA'].items():
            add(f'ctrl.size8.{tkey_}.{rk_}.r', rv_, None, 2, f'{CSS_}:comparison_with_targets.{tkey_}.size_relation_common_species.r_score_vs_HA.{rk_}',
                f"{sc_}: Pearson r between the mean score and the number of heavy atoms over these {zc_['n_species']} species, in {rlab_ if rk_ == '4DUH' else rk_}",
                check=D15['size'][tkey_]['common']['r'][rk_], tol=0.0006)
    add(f'ctrl.size.{tkey_}.ha_range', sr_['heavy_atom_range'], None, 0, f'{CSS_}:score_table.{tkey_}.size_relation.heavy_atom_range', f'{rlab_}, {sc_}: smallest and largest ligand (heavy atoms)')

# ---------------------------------------------------------------------------------------------- representative poses
RTAG_ = {C1_: 'primary', C1B_: 'second', C1W_: 'water', C1N_: 'noloop'}
for job_, q_ in CS_['jobs'].items():
    t_, run_ = job_.split('/runs/')
    kk_ = f"ctrl.cont.{RTAG_[t_]}.{run_.replace('__prod', '').replace('__', '.')}"
    for i_, pp_ in enumerate(q_['representative']['polar_pairs_NO_within_3p5A'][:8]):
        addbulk(f'{kk_}.{i_}', pp_['distance'], 'A', 2, f'{CSS_}:jobs["{job_}"].representative.polar_pairs_NO_within_3p5A[{i_}]',
                f"{RV_[t_]}, {NAME_[q_['ligand']]}, {q_['scoring']}, representative pose: {pp_['residue']} {pp_['receptor_atom']}···{pp_['ligand_atom'].split('#')[0]}",
                'ctrl.cont', 'control receptor 4DUH: polar N/O···N/O contacts (≤ 3.5 Å) of the representative poses')
jr_ = f'{C1_}/runs/rli_anion__vina'
cc_ = [pp_ for pp_ in CS_['jobs'][jr_]['representative']['polar_pairs_NO_within_3p5A']
       if pp_['receptor_atom'] == 'O' and pp_['ligand_atom'].startswith('O(carboxylate)')]
if cc_:
    add('ctrl.val.primary.vina.carboxylate_carbonyl', [min(x_['distance'] for x_ in cc_), max(x_['distance'] for x_ in cc_)], 'A', 1,
        f'{CSS_}:jobs["{jr_}"].representative.polar_pairs_NO_within_3p5A (receptor atom O, ligand carboxylate O)',
        'redock of the 4DUH inhibitor, primary receptor, vina, representative top pose: shortest and longest distance from a carboxylate oxygen of the turned benzoate '
        'to a backbone carbonyl oxygen within 3.5 Å (' + ', '.join(sorted({x_['residue'] for x_ in cc_}, key=lambda r_: int(r_[3:]))) + ')')
for tag_, kr_, kc_ in (('primary', 'ctrl.val.primary.vina.score', 'ctrl.score.cfx_zwit.vina'), ('water', 'ctrl.val.water.vina.score', 'ctrl.wat.cfx_zwit.vina')):
    if kr_ in FACTS and kc_ in FACTS:
        add(f'ctrl.gap.cfx.{tag_}', FACTS[kc_]['value'] - FACTS[kr_]['value'], 'kcal mol-1', 1, f'facts {kc_} − {kr_}',
            f"Vina, {'receptor as prepared' if tag_ == 'primary' else 'crystal water kept'}: score of ciprofloxacin minus score of the 4DUH inhibitor (positive: the inhibitor scores better)")

if 'ctrl.size.vina.slope' in FACTS and 'ctrl.score.cfx_zwit.vina.heavy' in FACTS:
    dha_ = FACTS['ctrl.ligand_heavy']['value'] - FACTS['ctrl.score.cfx_zwit.vina.heavy']['value']
    add('ctrl.gap.cfx.dheavy', dha_, None, 0, 'facts ctrl.ligand_heavy − ctrl.score.cfx_zwit.vina.heavy', 'heavy atoms of the 4DUH inhibitor minus those of ciprofloxacin')
    add('ctrl.gap.cfx.size_part', -FACTS['ctrl.size.vina.slope']['value'] * dha_, 'kcal mol-1', 1, 'facts −ctrl.size.vina.slope × ctrl.gap.cfx.dheavy',
        'Vina, receptor as prepared: the part of the score gap between the 4DUH inhibitor and ciprofloxacin that the size line gives for their difference in heavy atoms')
if 'ctrl.val.primary.vina.nn.score' in FACTS and 'ctrl.score.cfx_zwit.vina' in FACTS:
    add('ctrl.gap.cfx.nn', FACTS['ctrl.score.cfx_zwit.vina']['value'] - FACTS['ctrl.val.primary.vina.nn.score']['value'], 'kcal mol-1', 1, 'facts ctrl.score.cfx_zwit.vina − ctrl.val.primary.vina.nn.score',
        'Vina, receptor as prepared: score of ciprofloxacin minus the mean score of the near-native mode of the 4DUH inhibitor (positive: that mode scores better)')

# ---------------------------------------------------------------------------------------------- against the targets
# Differences at full precision: the control mean from the per-seed scores, the target mean from the bulk facts of the
# production table (section 21), so that a printed difference never comes from two rounded numbers.
for tkey_, rec_, sc_, pre_, rlab_ in TABS_:
    for pdb_, tg_ in TGT_.items():
        pairs_ = []
        for sp_ in ('sfa_neutral', 'sfa_cation', 'sfa_anion', 'smx_anion', 'smx_neutral', 'paba_anion', 'azm_neutral', 'azm_anion'):
            kc_, kt_ = f'{pre_}.{sp_}.{sc_}', f'bulk.dockprod.{tg_}.{sc_}.{sp_}.score_mean'
            if kc_ not in FACTS or kt_ not in FACTS:
                continue
            a_, b2_ = FACTS[kc_]['value'], FACTS[kt_]['value']
            pairs_.append((a_, b2_))
            addbulk(f'ctrl.cmp.{tkey_}.{pdb_}.{sp_}', a_ - b2_, 'kcal mol-1', 1, f'facts {kc_} − {kt_}',
                    f'{NAME_[sp_]}, {sc_}: score in {rlab_} minus score in {pdb_} (positive: weaker in the control receptor)', 'ctrl.cmp',
                    'control receptor 4DUH against the targets: score differences of the species of the production matrix (Table S34)',
                    check=(D15['targets'].get(f'{tkey_}/{pdb_}', {}).get('per_species', {}).get(sp_)), tol=1e-6)
        if len(pairs_) < 4:
            continue
        df_ = [a_ - b2_ for a_, b2_ in pairs_]
        d_ = D15['targets'].get(f'{tkey_}/{pdb_}', {})
        q_ = CS_['comparison_with_targets'][tkey_][f'vs_{pdb_}']
        assert q_['n_species'] == len(pairs_)
        add(f'ctrl.cmp.{tkey_}.{pdb_}.n', len(pairs_), None, 0, f'{CSS_}:comparison_with_targets.{tkey_}.vs_{pdb_}.n_species', f'{sc_}: species docked into both {rlab_} and {pdb_}')
        add(f'ctrl.cmp.{tkey_}.{pdb_}.mean', _stt.mean(df_), 'kcal mol-1', 1, f'facts ctrl.cmp.{tkey_}.{pdb_}.* (mean)',
            f'{sc_}: mean of (score in {rlab_} − score in {pdb_}) over the {len(pairs_)} species of the production matrix', check=d_.get('mean_difference'), tol=1e-6)
        add(f'ctrl.cmp.{tkey_}.{pdb_}.range', [min(df_), max(df_)], 'kcal mol-1', 1, f'facts ctrl.cmp.{tkey_}.{pdb_}.* (min, max)',
            f'{sc_}: smallest and largest of these differences')
        add(f'ctrl.cmp.{tkey_}.{pdb_}.r', q_['pearson_r'], None, 2, f'{CSS_}:comparison_with_targets.{tkey_}.vs_{pdb_}.pearson_r',
            f'{sc_}: Pearson r between the scores of the {q_["n_species"]} species in {rlab_} and in {pdb_}', check=d_.get('pearson_r'), tol=0.0006)
        add(f'ctrl.cmp.{tkey_}.{pdb_}.rho', q_['spearman_rho'], None, 2, f'{CSS_}:comparison_with_targets.{tkey_}.vs_{pdb_}.spearman_rho',
            f'{sc_}: Spearman rank correlation of the same scores (p = {q_["spearman_p"]})')
V4_ = [FACTS[f'dock.prod.{p_}.sfa_neutral.vina.score']['value'] for p_ in ('3TZF', '7S2K', '6RL9', '3HS4')]
add('ctrl.sfa.targets_vina_range', [min(V4_), max(V4_)], 'kcal mol-1', 1, 'facts dock.prod.{3TZF,7S2K,6RL9,3HS4}.sfa_neutral.vina.score',
    'neutral sulfanilamide, Vina: lowest and highest mean score over the four production structures (3TZF, 7S2K, 6RL9, 3HS4)')
AV_ = CS_['comparison_with_targets'].get('sulfanilamide_vina_all_variants')
if AV_:
    vv_ = [_score(f'{t_}/runs/sfa_neutral__vina__prod')[0] for t_ in AV_['control']]
    add('ctrl.sfa.variants_vina_range', [min(vv_), max(vv_)], 'kcal mol-1', 1, f'{CSS_}:jobs[<variant>/runs/sfa_neutral__vina__prod].score.per_seed (means)',
        f"neutral sulfanilamide, Vina: lowest and highest mean score over the {len(vv_)} variants of the control receptor ({', '.join(RV_[t_] for t_ in AV_['control'])})",
        check=[D15['sfa_variants']['control_min'], D15['sfa_variants']['control_max']], tol=1e-6)
    add('ctrl.sfa.n_variants', len(vv_), None, 0, f'{CSS_}:comparison_with_targets.sulfanilamide_vina_all_variants.control', 'variants of the control receptor into which neutral sulfanilamide was docked with Vina')
    V4F_ = [FACTS[f'bulk.dockprod.{t_}.vina.sfa_neutral.score_mean']['value'] for t_ in TGT_.values()]      # full precision, as for ctrl.cmp
    add('ctrl.sfa.max_gap', max(abs(a_ - b2_) for a_ in vv_ for b2_ in V4F_), 'kcal mol-1', 1, 'facts: variants of the control receptor against bulk.dockprod.{T1_3TZF,T5_7S2K,T3_6RL9,T4_3HS4}.vina.sfa_neutral.score_mean',
        'neutral sulfanilamide, Vina: largest absolute difference between its mean score in any variant of the control receptor and in any production structure',
        check=D15['sfa_variants']['max_abs_difference'], tol=1e-6)
if 'ctrl.score.sfa_neutral.vinardo' in FACTS:
    N4_ = [FACTS[f'dock.prod.{p_}.sfa_neutral.vinardo.score']['value'] for p_ in ('3TZF', '7S2K', '6RL9', '3HS4')]
    add('ctrl.sfa.targets_vinardo_range', [min(N4_), max(N4_)], 'kcal mol-1', 1, 'facts dock.prod.{3TZF,7S2K,6RL9,3HS4}.sfa_neutral.vinardo.score',
        'neutral sulfanilamide, Vinardo: lowest and highest mean score over the four production structures')

# ---------------------------------------------------------------------------------------------- receptor variants
for sp_, row_ in LS_['jobs'].items():
    for tag_, lab_ in (('second_build', 'second loop build'), ('no_loop', 'built residues removed')):
        if f'score_{tag_}' not in row_:
            continue
        t_ = C1B_ if tag_ == 'second_build' else C1N_
        job_ = next(k for k, q in CS_['jobs'].items() if q['target'] == t_ and q['ligand'] == sp_)
        k_ = f"ctrl.var.{'second' if tag_ == 'second_build' else 'noloop'}.{sp_}"
        sm_, ss_ = _score(job_)
        p_ = FACTS[f'ctrl.val.primary.vina.score' if sp_ == 'rli_anion' else f'ctrl.score.{sp_}.vina']['value']
        add(k_ + '.score', sm_, 'kcal mol-1', 1, f'{CSS_}:jobs["{job_}"].score.per_seed (mean)', f'{NAME_[sp_]}, Vina, {lab_}: mean top-pose score', check=D15['scores'][job_]['mean'], tol=1e-6)
        add(k_ + '.sd', ss_, 'kcal mol-1', 2, f'{CSS_}:jobs["{job_}"].score.per_seed (SD)', f'{NAME_[sp_]}, Vina, {lab_}: SD of the top-pose score', check=D15['scores'][job_]['sd'], tol=1e-6)
        add(k_ + '.diff', sm_ - p_, 'kcal mol-1', 2, f'facts {k_}.score − primary receptor', f'{NAME_[sp_]}, Vina: score with {lab_} minus score in the primary receptor')
        add(k_ + '.rep_rmsd', row_[f'representative_rmsd_to_primary_{tag_}_A'], 'A', 1, f'{CSS_}:loop_sensitivity.jobs.{sp_}.representative_rmsd_to_primary_{tag_}_A',
            f'{NAME_[sp_]}, Vina: RMSD between the representative poses with {lab_} and in the primary receptor (same frame)')
        add(k_ + '.in_site', round(row_[f'in_site_fraction_{tag_}'] * 10), None, 0, f'{CSS_}:loop_sensitivity.jobs.{sp_}.in_site_fraction_{tag_}', f'{NAME_[sp_]}, Vina, {lab_}: seeds whose top pose lies in the inhibitor site')
for sp_, row_ in WD_['jobs'].items():
    kw_ = f'ctrl.wat.{sp_}.vina'
    p_ = FACTS['ctrl.val.primary.vina.score' if sp_ == 'rli_anion' else f'ctrl.score.{sp_}.vina']['value']
    add(f'ctrl.var.water.{sp_}.diff', FACTS[kw_]['value'] - p_, 'kcal mol-1', 2, f'facts {kw_} − primary receptor',
        f'{NAME_[sp_]}, Vina: score with the crystal water kept minus score in the receptor as prepared')

# ---------------------------------------------------------------------------------------------- whole-fragment search
for job_, w_ in CS_['whole_fragment'].items():
    sp_ = job_.split('/runs/')[1].split('__')[0]
    k_, d_ = f'ctrl.whole.{sp_}', D15['whole'][job_]
    add(k_ + '.n_seeds', w_['n_seeds'], None, 0, f'{CSS_}:whole_fragment["{job_}"].n_seeds', f'whole-fragment search, {NAME_[sp_]}: seeds')
    sm_, ss_ = _ms([x_['score_top'] for x_ in w_['per_seed']])
    add(k_ + '.score', sm_, 'kcal mol-1', 1, f'{CSS_}:whole_fragment["{job_}"].per_seed[].score_top (mean)', f'whole-fragment search, {NAME_[sp_]}, Vina: mean top-pose score', check=d_['mean'], tol=1e-6)
    add(k_ + '.sd', ss_, 'kcal mol-1', 2, f'{CSS_}:whole_fragment["{job_}"].per_seed[].score_top (SD)', f'whole-fragment search, {NAME_[sp_]}, Vina: SD of the top-pose score', check=d_['sd'], tol=1e-6)
    add(k_ + '.n_in_site', w_['n_top_in_site'], None, 0, f'{CSS_}:whole_fragment["{job_}"].n_top_in_site', f'whole-fragment search, {NAME_[sp_]}: seeds whose top pose lies in the inhibitor site', check=d_['n_top_in_site'], tol=0.01)
    add(k_ + '.n_any_in_site', w_['seeds_with_any_in_site_mode'], None, 0, f'{CSS_}:whole_fragment["{job_}"].seeds_with_any_in_site_mode', f'whole-fragment search, {NAME_[sp_]}: seeds with any mode in the inhibitor site', check=d_['n_any_in_site'], tol=0.01)
    if w_.get('rmsd_top_to_reference'):
        q_ = w_['rmsd_top_to_reference']
        add(k_ + '.rmsd_to_focused', _ms([x_['rmsd_top_to_reference_A'] for x_ in w_['per_seed']])[0], 'A', 2, f'{CSS_}:whole_fragment["{job_}"].per_seed[].rmsd_top_to_reference_A (mean)',
            f'whole-fragment search, {NAME_[sp_]}: RMSD of the top poses to the {q_["reference"]}, mean over seeds')
        add(k_ + '.rmsd_to_focused_max', q_['max'], 'A', 2, f'{CSS_}:whole_fragment["{job_}"].rmsd_top_to_reference.max', f'whole-fragment search, {NAME_[sp_]}: the same, largest over seeds')
        add(k_ + '.n_same_pose', round(q_['f_le_2A'] * w_['n_seeds']), None, 0, f'{CSS_}:whole_fragment["{job_}"].rmsd_top_to_reference.f_le_2A', f'whole-fragment search, {NAME_[sp_]}: seeds whose top pose lies within 2.0 Å of the pose of the focused search')
add('ctrl.n_jobs', len(CS_['jobs']) + len(CS_['whole_fragment']), None, 0, f'{CSS_}:jobs, whole_fragment', 'docking jobs of the control receptor (focused and whole-fragment; all variants)')
if CS_['pending_jobs']:
    MISSING.append('control receptor: jobs not finished when the facts were built: ' + ', '.join(CS_['pending_jobs']))
