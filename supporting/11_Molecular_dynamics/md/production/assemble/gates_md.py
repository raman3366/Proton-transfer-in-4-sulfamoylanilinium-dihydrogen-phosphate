import json, csv, statistics as st
from pathlib import Path
R = Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md')
RUNS = Path('/root/cqc-sfa-md/runs')
GZ = {s: json.load(open(R/'production'/f'gate_GZ_bonded_{s}.json')) for s in ('2s','2e')}
SM = {s: json.load(open(RUNS/f'bonded_{s}'/'X0'/'analysis'/'summary.json')) for s in ('2s','2e')}
SM2 = {s: json.load(open(RUNS/f'bonded_{s}'/'X0'/'analysis_recheck'/'summary.json')) for s in ('2s','2e')}
RR = {s: json.load(open(RUNS/f'bonded_{s}'/'X0'/'analysis_recheck'/'recheck_report.json')) for s in ('2s','2e')}
MAN = {s: json.load(open(RUNS/f'bonded_{s}'/'X0'/'manifest.json')) for s in ('2s','2e')}
out = []; P = out.append
ok = lambda b: 'met' if b else '**NOT MET**'
names = {'d_Zn_N2_A': 'Zn–N2 (Å)', 'd_Zn_H94NE2_A': 'Zn–NE2 His94 (Å)', 'd_Zn_H96NE2_A': 'Zn–NE2 His96 (Å)', 'd_Zn_H119ND1_A': 'Zn–ND1 His119 (Å)',
         'ang_N2_Zn_H94NE2_deg': 'N2–Zn–NE2 His94 (°)', 'ang_N2_Zn_H96NE2_deg': 'N2–Zn–NE2 His96 (°)', 'ang_N2_Zn_H119ND1_deg': 'N2–Zn–ND1 His119 (°)',
         'ang_H94NE2_Zn_H96NE2_deg': 'NE2 His94–Zn–NE2 His96 (°)', 'ang_H94NE2_Zn_H119ND1_deg': 'NE2 His94–Zn–ND1 His119 (°)', 'ang_H96NE2_Zn_H119ND1_deg': 'NE2 His96–Zn–ND1 His119 (°)'}
# ---------- G-Z
P('### G-Z table')
P('| Part | Quantity | Must lie in | 2s: minimised | 2s: mean of the 10 frames | 2s: frame at 100 ps | 2s | 2e: minimised | 2e: mean of the 10 frames | 2e: frame at 100 ps | 2e |')
P('|---|---|---|---|---|---|---|---|---|---|---|')
def gzrow(part, grp, key, unit):
    cells = []
    for s in ('2s','2e'):
        g = GZ[s]['gate_GZ']; f = GZ[s]['for_information']
        g1 = g[grp] if key is None else g[grp][key]
        h1 = f['frame_at_100_ps_alone'][grp] if key is None else f['frame_at_100_ps_alone'][grp][key]
        col = key or 'd_Zn_N2_A'
        mv = f['minimisation']['structure'][col]
        w = g1['window_'+unit]
        cells += [f"{mv:.3f}", f"{g1['mean_'+unit]:.3f}", f"{h1['mean_'+unit]:.3f}", ok(g1['passed'])]
    P(f"| {part} | {names[key or 'd_Zn_N2_A']} | {w[0]:g} to {w[1]:g} | " + ' | '.join(cells) + ' |')
gzrow('a', 'a_Zn_N2', None, 'A')
for k in ('d_Zn_H94NE2_A','d_Zn_H96NE2_A','d_Zn_H119ND1_A'): gzrow('b', 'b_Zn_NHis', k, 'A')
for k in list(names)[4:]: gzrow('c', 'c_N_Zn_N', k, 'deg')
for s in ('2s','2e'):
    g = GZ[s]['gate_GZ']; d = g['d_fifth_ligand']; f = GZ[s]['for_information']
    P(f"{s}: part d: sulfonyl O within 2.6 in {d['solute_frames']['count']} of {d['solute_frames']['n']} solute frames (nearest {d['solute_frames']['nearest_sulfonyl_O_min_A']}); any O in {d['full_frames']['count']} of {d['full_frames']['n']} full frames (nearest any {d['full_frames']['nearest_O_of_any_kind_min_A']}, water {d['full_frames']['nearest_water_O_min_A']}, protein {d['full_frames']['nearest_protein_O_min_A']}); passed {d['passed']}")
    gl = g['reported_zn_glu106']
    P(f"{s}: Glu106 OE1 mean {gl['d_Zn_E106OE1_A']['mean_A']} min {gl['d_Zn_E106OE1_A']['min_A']} within {gl['d_Zn_E106OE1_A']['count_within']}/{gl['d_Zn_E106OE1_A']['n']}; OE2 mean {gl['d_Zn_E106OE2_A']['mean_A']} min {gl['d_Zn_E106OE2_A']['min_A']}")
    P(f"{s}: verdict {g['verdict']}; frames {GZ[s]['frames_of_the_gate']['solute_frames']['frame_numbers_of_the_equilibration_series']} times {GZ[s]['frames_of_the_gate']['solute_frames']['time_ns_of_the_series']}; full {GZ[s]['frames_of_the_gate']['full_frames']['frame_numbers_of_the_equilibration_series']} {GZ[s]['frames_of_the_gate']['full_frames']['time_ns_of_the_series']}; T end {GZ[s]['frames_of_the_gate']['temperature_at_the_end_K']}; second program {GZ[s]['second_program']['agree']} n={GZ[s]['second_program']['values_compared_on_the_gate_frames']} worst={GZ[s]['second_program']['largest_difference']}")
    ms = f['minimisation']
    P(f"{s}: minimised: O1 {ms['structure']['d_Zn_O1_A']} O2 {ms['structure']['d_Zn_O2_A']} E106OE1 {ms['structure']['d_Zn_E106OE1_A']}; each in window {ms['structure']['each_value_in_its_window']}; lbfgs {ms['lbfgs_iterations_per_pass']}; restrained displacement {ms['restrained_atoms_displacement_A']}")
    P(f"{s}: script sha {GZ[s]['script']['sha256'][:16]} written {GZ[s]['written']['ist']}")
P('')
P('### by stage (means)')
P('| Stage | Held atoms, k | 2s: Zn–N2 | 2s: Zn–N His94/96/119 | 2s: six angles | 2s: O within 2.6 Å, solute / all-atom | 2s a–d | 2e: Zn–N2 | 2e: Zn–N His94/96/119 | 2e: six angles | 2e: O within 2.6 Å | 2e a–d |')
P('|---|---|---|---|---|---|---|---|---|---|---|---|')
stn = sorted(GZ['2s']['for_information']['by_equilibration_stage'], key=lambda k: int(k.split('_')[1]))
for sk in stn:
    cells = []
    for s in ('2s','2e'):
        p = GZ[s]['for_information']['by_equilibration_stage'][sk]
        b = p['b_Zn_NHis']; c = p['c_N_Zn_N']; d = p['d_fifth_ligand']
        cells += [f"{p['a_Zn_N2']['mean_A']:.3f}", ' / '.join(f"{b[k]['mean_A']:.3f}" for k in b), ' / '.join(f"{c[k]['mean_deg']:.1f}" for k in c),
                  f"{d['solute_frames']['count']} of {d['solute_frames']['n']} / {d['full_frames'].get('count')} of {d['full_frames']['n']}", ok(p['passed']) + ('' if p['passed'] else ' (' + '; '.join(p['failed_parts']) + ')')]
    stg = GZ['2s']['for_information']['by_equilibration_stage'][sk]['settings']
    P(f"| {sk.split('_')[1]} ({stg['name']}, {stg['length_ps']:g} ps) | {stg['restraint_atoms']}, {stg['restraint_k_kj_per_mol_nm2']:g} | " + ' | '.join(cells) + ' |')
P('')
# ---------- G-V
P('### G-V table')
P('| Criterion | Quantity | Must be | 2s (X0) | 2s | 2e (X0) | 2e |')
P('|---|---|---|---|---|---|---|')
g = {s: SM[s]['gate_GV'] for s in SM}
P(f"| GV1 | mean heavy-atom RMSD of the ligand from the crystal pose (Å) | ≤ 1.0 | {g['2s']['GV1_lig_rmsd']['mean_A']:.6f} | {ok(g['2s']['GV1_lig_rmsd']['passed'])} | {g['2e']['GV1_lig_rmsd']['mean_A']:.6f} | {ok(g['2e']['GV1_lig_rmsd']['passed'])} |")
P(f"| GV2 | frames with the hydrogen bond N2–H···OG1 Thr199 | ≥ 80 % | {g['2s']['GV2_hbond']['count']} of {g['2s']['GV2_hbond']['n']} ({100*g['2s']['GV2_hbond']['fraction']:.1f} %) | {ok(g['2s']['GV2_hbond']['passed'])} | {g['2e']['GV2_hbond']['count']} of {g['2e']['GV2_hbond']['n']} ({100*g['2e']['GV2_hbond']['fraction']:.1f} %) | {ok(g['2e']['GV2_hbond']['passed'])} |")
for k in names:
    a, b = g['2s']['GV3_zinc']['criteria'][k], g['2e']['GV3_zinc']['criteria'][k]
    part = 'GV3 a' if k == 'd_Zn_N2_A' else ('GV3 b' if k.startswith('d_') else 'GV3 c')
    if k.startswith('d_'):
        P(f"| {part} | mean {names[k]} | {a['window_A'][0]:g} to {a['window_A'][1]:g} | {a['mean_A']:.6f} | {ok(a['passed'])} | {b['mean_A']:.6f} | {ok(b['passed'])} |")
    else:
        P(f"| {part} | mean {names[k]}; 6RL9 {a['reference_6RL9_deg']:g} | {a['window_deg'][0]:g} to {a['window_deg'][1]:g} | {a['mean_deg']:.6f} ({a['deviation_deg']:+.6f}) | {ok(a['passed'])} | {b['mean_deg']:.6f} ({b['deviation_deg']:+.6f}) | {ok(b['passed'])} |")
fa, fb = g['2s']['GV3_zinc']['criteria']['fifth_ligand'], g['2e']['GV3_zinc']['criteria']['fifth_ligand']
P(f"| GV3 d | solute frames with a sulfonyl O within 2.6 Å of the zinc | ≤ 5 % | {fa['solute_trajectory']['count']} of {fa['solute_trajectory']['n']} ({100*fa['solute_trajectory']['fraction']:.1f} %) | {ok(fa['solute_trajectory']['passed'])} | {fb['solute_trajectory']['count']} of {fb['solute_trajectory']['n']} ({100*fb['solute_trajectory']['fraction']:.1f} %) | {ok(fb['solute_trajectory']['passed'])} |")
P(f"| GV3 d | all-atom frames with any O (sulfonyl, water, protein) within 2.6 Å of the zinc | ≤ 5 %, at least one frame | {fa['full_trajectory']['count']} of {fa['full_trajectory']['n']} ({100*fa['full_trajectory']['fraction']:.1f} %) | {ok(fa['full_trajectory']['passed'])} | {fb['full_trajectory']['count']} of {fb['full_trajectory']['n']} ({100*fb['full_trajectory']['fraction']:.1f} %) | {ok(fb['full_trajectory']['passed'])} |")
P(f"| GV3 | all parts | | | {ok(g['2s']['GV3_zinc']['passed'])} | | {ok(g['2e']['GV3_zinc']['passed'])} |")
for k, nmk in (('ang_C4_S_N2_deg', 'GV4 a | mean C4–S–N2 (°), s.d.'), ('tau_fold_deg', 'GV4 b | mean τ_fold (°), s.d.')):
    a, b = g['2s']['GV4_ligand_internal']['criteria'][k], g['2e']['GV4_ligand_internal']['criteria'][k]
    P(f"| {nmk} | {a['window_deg'][0]:g} to {a['window_deg'][1]:g} | {a['mean_deg']:.6f}, {a['sd_deg']:.6f} | {ok(a['passed'])} | {b['mean_deg']:.6f}, {b['sd_deg']:.6f} | {ok(b['passed'])} |")
P(f"| **Gate G-V** | `gate_GV.verdict` | | | **{g['2s']['verdict']}** | | **{g['2e']['verdict']}** |")
for s in ('2s','2e'):
    gl = g[s]['reported_zn_glu106']
    P(f"{s}: Zn-Glu106 OE1 mean {gl['d_Zn_E106OE1_A']['mean_A']} min {gl['d_Zn_E106OE1_A']['min_A']} within {gl['d_Zn_E106OE1_A']['count_within']} of {gl['d_Zn_E106OE1_A']['n']}; OE2 mean {gl['d_Zn_E106OE2_A']['mean_A']} min {gl['d_Zn_E106OE2_A']['min_A']} within {gl['d_Zn_E106OE2_A']['count_within']}")
    P(f"{s}: provisional {g[s]['provisional']}, applies {g[s]['applies_to_this_run']}, role {SM[s]['run']['role']}; second program verdict {SM2[s]['gate_GV']['verdict']} passed {SM2[s]['gate_GV']['passed']}; warnings {SM[s].get('warnings')} / {SM2[s].get('warnings')}")
    # recheck report
    rr = RR[s]; print_keys = [k for k in rr.keys()]
    P(f"{s}: recheck report keys {print_keys}")
    m = MAN[s]; pr = m['production']; eq = m['equilibration']
    P(f"{s}: manifest: run_seed {m['seed']['run_seed']}; created {m['created']['ist']}; last written {m['last_written']['ist']}; min wall {m['minimisation']['wall_s']} s; equil wall {eq['wall_s']} s {eq['ns_per_day']} ns/day; prod wall {pr['wall_s']} s, {pr['ns_per_day']} ns/day; frames {pr['frames_solute']}/{pr['frames_full']}; hydration {json.dumps({k: v for k, v in m['hydration_check'].items() if k not in ('settings','rule','checks')})[:600]}")
    P(f"{s}: inputs prmtop {m['inputs']['system.prmtop']['sha256']}; script {m['script']['sha256'][:16]}; model {m['build_json']['model']}")
    sol = list(csv.DictReader(open(RUNS/f'bonded_{s}'/'X0'/'analysis'/'timeseries_solute.csv')))
    full = list(csv.DictReader(open(RUNS/f'bonded_{s}'/'X0'/'analysis'/'timeseries_full.csv')))
    v = [float(r['d_Zn_Osulfonyl_min_A']) for r in sol]
    P(f"{s}: nearest sulfonyl O in the solute frames: mean {st.mean(v):.4f} sd {st.stdev(v):.4f} min {min(v):.4f} max {max(v):.4f}; <=2.6: {sum(x<=2.6 for x in v)}; <=2.7: {sum(x<=2.7 for x in v)}; <=2.8: {sum(x<=2.8 for x in v)}; O1<=2.6 {sum(float(r['d_Zn_O1_A'])<=2.6 for r in sol)} O2<=2.6 {sum(float(r['d_Zn_O2_A'])<=2.6 for r in sol)}; mean O1 {st.mean(float(r['d_Zn_O1_A']) for r in sol):.4f} O2 {st.mean(float(r['d_Zn_O2_A']) for r in sol):.4f}; mean S-N2-Zn {st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol):.3f}")
    hits = [(r['time_ns'], r['d_Zn_Osulfonyl_min_full_A'], r['d_Zn_Owater_min_A'], r['d_Zn_Oprotein_min_A']) for r in full if r['fifth_any'] == '1']
    P(f"{s}: all-atom frames with fifth_any: {hits}; water min {min(float(r['d_Zn_Owater_min_A']) for r in full):.4f}; protein O min {min(float(r['d_Zn_Oprotein_min_A']) for r in full):.4f}; fifth_water {sum(int(r['fifth_water']) for r in full)} fifth_proteinO {sum(int(r['fifth_proteinO']) for r in full)}")
    P(f"{s}: solute hits at ns {[r['time_ns'] for r in sol if r['fifth_sulfonylO']=='1']}")
    tr = SM[s].get('table_row') or {}
    P(f"{s}: table_row {json.dumps(tr)[:900]}")
    P(f"{s}: orientation {json.dumps(SM[s].get('orientation'))[:500]}")
print('\n'.join(out))
