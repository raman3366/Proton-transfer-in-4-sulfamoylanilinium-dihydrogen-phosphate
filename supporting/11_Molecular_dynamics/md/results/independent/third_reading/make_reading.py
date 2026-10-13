#!/usr/bin/env python
"""Puts third_reading.json together from runs_raw.json (third_reader.py), the
per-run series (npz), the manifests and the two check files.  No trajectory is
read here."""
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import third_lib as T

HERE = os.path.dirname(os.path.abspath(__file__))
raw = json.load(open(os.path.join(HERE, 'runs_raw.json')))
now = subprocess.run(['date', '+%-d %b %Y, %-I:%M %p IST'], env=dict(os.environ, TZ='Asia/Kolkata'),
                     capture_output=True, text=True).stdout.strip()

runs = {}
for zset, run, start, role in T.RUN_LIST:
    key = zset + '/' + run
    v = raw[key]
    z = np.load(os.path.join(HERE, 'series_%s_%s.npz' % (zset, run)))
    st = z['state']
    N = len(st)
    s_state = v['state_of_start']
    opp = v['opposite_state']
    # extras from the series
    eq_tags = ['eq04', 'eq05', 'eq06']
    eq_states = {}
    hb_eq = []
    for t in eq_tags:
        e = T.states_of(z['eq_%s_d1' % t], z['eq_%s_d2' % t])
        eq_states[t] = ''.join(e)
        hb_eq.append(z['eq_%s_hb' % t])
    hb_eq = np.concatenate(hb_eq)
    a = list(z['assigned'])
    vis = T.visits(a)
    start_vis = [x for x in vis if x[0] == s_state]
    w = min(N, T.WINDOW)
    vis_w = T.visits(a[N - w:])
    near = int(np.sum((np.abs(z['d1'] - T.M1) < 1e-4) | (np.abs(z['d2'] - T.M2) < 1e-4)))
    man = T.load_manifest(v['folder'])
    hy = man.get('hydration_check', {})
    v['extras'] = {
        'states_frame_by_frame_stage4_5_6': eq_states,
        'hbond_in_stage4_5_6': {'count': int(hb_eq.sum()), 'n': int(len(hb_eq))},
        'visits_in_start_state_production': {
            'n': len(start_vis), 'longest_frames': max([x[2] for x in start_vis] or [0]),
            'total_assigned_frames': int(sum(x[2] for x in start_vis))},
        'visits_in_last_20ns_by_assigned_state': {s: sum(1 for x in vis_w if x[0] == s) for s in 'CR'},
        'frames_within_0.0001A_of_a_midpoint': near,
        'hydration_check_of_the_driver_from_manifest': {
            'applies': hy.get('applies'), 'counts': hy.get('counts'), 'n_extensions': hy.get('n_extensions'),
            'final_decision': hy.get('final_decision'),
            'site_without_water_when_the_solute_was_released': hy.get('site_without_water_when_the_solute_was_released')},
    }
    c = v['conversion']
    vd = v['verdict']
    rr = v['lig_rmsd_A']
    zn = v['zinc']
    n2 = v['n2_configuration']
    j = n2.get('joined_release_series', {})
    flc = j.get('first_lasting_change')
    v['key_numbers_one_line'] = (
        '%s start %s (%s): N=%d; d1 mean %.3f, d2 mean %.3f A; states C/R/I whole %d/%d/%d, last 20 ns %d/%d/%d; verdict %s; '
        'visits %d (long %d), transitions %d, first passage in production %s; conversion time: %s; '
        'lig RMSD mean %.3f, last half %.3f A; H-bond %.4f; Zn-N2 %.3f, Zn-N(His94/96/119) %.3f/%.3f/%.3f A; '
        'sulfonyl O within 2.6 A of Zn %.4f (%d frames); N2 sign +/- %.4f/%.4f, first lasting change %s'
        % (key, v['start'], v['role'], N, v['d1_A']['whole']['mean'], v['d2_A']['whole']['mean'],
           v['states']['whole']['C'], v['states']['whole']['R'], v['states']['whole']['I'],
           v['states']['last_20ns']['C'], v['states']['last_20ns']['R'], v['states']['last_20ns']['I'],
           vd['verdict'], v['visits']['n_visits'], v['visits']['n_long_visits'], v['visits']['n_transitions'],
           ('none' if v['visits']['first_passage_ns'] is None else '%.2f ns' % v['visits']['first_passage_ns']),
           c['report'], rr['mean_whole'], rr['mean_last_half'], v['hbond_N2H_T199OG1']['fraction'],
           zn['mean_Zn_N2_A'], zn['mean_Zn_H94NE2_A'], zn['mean_Zn_H96NE2_A'], zn['mean_Zn_H119ND1_A'],
           zn['sulfonyl_O_within_2.6A']['fraction'], zn['sulfonyl_O_within_2.6A']['count'],
           n2['fraction_plus'], n2['fraction_minus'],
           ('none' if flc is None else '%d ps after release (%s)' % (round(flc['t_from_release_ps']), flc['where']))))
    runs[key] = v

# ---------------------------------------------------------------- rule 1
r1 = {}
left = 0
for run in ('X1', 'X2', 'X3'):
    v = runs['bonded_2e/' + run]
    m = v['lig_rmsd_A']['mean_last_half']
    yes = bool(m > 2.0)
    left += yes
    r1[run] = {'mean_lig_rmsd_last_half_A': m, 'n_frames_last_half': v['lig_rmsd_A']['n_frames_last_half'],
               'above_2A': 'yes' if yes else 'no',
               'frames_above_2A_in_last_half': v['lig_rmsd_A']['frames_above_2A_last_half'],
               'max_lig_rmsd_whole_run_A': v['lig_rmsd_A']['max']}
rule1 = {
    'rule': 'PLAN.md section 5 rule 1, read by the entry of 5 Oct 2026, 8:40 PM IST, point 2: the mean ligand RMSD over the last half of the run above 2.0 A',
    'per_run': r1, 'runs_that_left_the_crystal_orientation': int(left), 'of': 3,
    'reading': ('The crystal-start runs leave the crystal orientation in %d of three runs. The condition of rule 1 '
                '(two or three of three) is %s: the validation is %s.' %
                (left, 'met' if left >= 2 else 'not met', 'FAILED' if left >= 2 else 'not failed')),
}

# ---------------------------------------------------------------- rule 3
ver = {k.split('/')[1]: runs[k]['verdict']['verdict'] for k in runs if k.startswith('bonded_2e/')}
nR = sum(ver[r] == 'CONVERTED' for r in ('R1', 'R2', 'R3'))
xn = ['X1', 'X2', 'X3', 'N1', 'N2', 'N3']
xn_stay = sum(ver[r] == 'STAYED' for r in xn)
n_conv_all = sum(v == 'CONVERTED' for v in ver.values())
if nR >= 2 and xn_stay == 6:
    branch = 'first'
    words = ('Rotated starts convert in three of three runs' if nR == 3 else 'Rotated starts convert in two of three runs') + \
            ', crystal and near-native starts stay: the rotated orientation is not stable outside the scoring function.'
elif n_conv_all == 0 and all(v == 'STAYED' for v in ver.values()):
    branch = 'second'
    words = 'No run converts: both orientations persisted for 100 ns; the simulation does not decide.'
else:
    branch = 'third'
    words = 'Mixed results, reported run by run; no average is taken over runs that disagree.'
rule3 = {
    'rule': 'PLAN.md section 5 rule 3; a run converts by the 90 % rule on its last 20 ns (entry of 5 Oct 2026, 8:40 PM IST, point 1); '
            'a rotated run whose production holds no rotated frame counts as converted (entry of 6 Oct 2026, 11:40 AM IST, point 2)',
    'verdicts_of_the_nine_runs': ver,
    'rotated_starts_converted': int(nR), 'crystal_and_near_native_starts_stayed': int(xn_stay),
    'branch': branch, 'reading_in_the_plans_words': words,
    'conversions_of_the_rotated_starts': {
        r: {'report': runs['bonded_2e/' + r]['conversion']['report'],
            't_from_release_ps': runs['bonded_2e/' + r]['conversion']['t_from_release_ps'],
            'stage': runs['bonded_2e/' + r]['conversion'].get('stage_of_that_frame'),
            'rotated_frames_in_production': runs['bonded_2e/' + r]['states']['whole']['R'],
            'last_20ns_C_R_I': [runs['bonded_2e/' + r]['states']['last_20ns'][s] for s in 'CRI']}
        for r in ('R1', 'R2', 'R3')},
    'sensitivity_runs_in_no_count_of_rule_3': {
        k: {'verdict': runs[k]['verdict']['verdict'],
            'last_20ns_C_R_I': [runs[k]['states']['last_20ns'][s] for s in 'CRI'],
            'whole_run_C_R_I': [runs[k]['states']['whole'][s] for s in 'CRI'],
            'time_by_the_rule_of_7_Oct': runs[k]['conversion']['report'],
            'visits_in_R_in_production': runs[k]['extras']['visits_in_start_state_production']}
        for k in ('bonded_2s/R1', 'bonded_2e_n2fit/R1')},
    'rule_4': 'The outcome is a conversion, so rule 4 applies; the two sensitivity runs of start R are those fixed on 6 Oct 2026, 9:15 PM IST, point 3.',
    'said_beside_the_reading': [
        'All three conversions of the primary set happened in the equilibration (stage 4 or 5, backbone still restrained), 160 to 540 ps after the ligand was released; no production frame of R1, R2 or R3 is in the rotated state. The production therefore shows the crystal-like state being kept after the conversion, not the conversion itself.',
        'The sensitivity run with the candidate terms at N2 (bonded_2e_n2fit/R1; terms not accepted by gate G-R2) is MIXED by the 90 % rule: it leaves the rotated state 60 ps after release by the time rule, but returns to it again and again (no visit there of 1 ns or more). Under these terms the rotated state is revisited, so the reading of rule 3 depends on the terms at N2, which failed gates G-R and G-R2.',
        'The sensitivity run with zinc set 2s (bonded_2s/R1; the set that failed gate G-V narrowly) is CONVERTED.',
    ],
}

definitions = {
    'blind': 'No file of any analysis/, analysis_recheck/ or analysis_campaign*/ folder, of md/results/ or md/RESULTS.md was opened; md/analysis/analyse_md.py and recheck_md.py were neither opened nor imported. Read: md/PLAN.md, md/analysis/SPEC.md sections 2 to 6 and 9E, md/protocol.json, md/system/{residue_map.csv,starts.csv,protein_md.pdb,ligand_X.pdb,zinc.pdb}, reference.pdb of the built systems (for one comparison), and of each run folder solute.pdb, manifest.json and solute_*.dcd.',
    'code': 'third_lib.py (PDB by fixed columns, own DCD reader with numpy, imaging, Kabsch fit, states, visits), third_reader.py (per run), make_reading.py (this file), check_starts.py and crosscheck_mda.py (checks). One process, nice -n 10, one thread.',
    'atoms': 'By name in solute.pdb, residues by the tleap index of md/system/residue_map.csv: Phe131 CZ (index 129), Asn67 OD1 (66), Thr199 OG1 (197), His94 NE2 (93, residue HD1), His96 NE2 (95, HD2), His119 ND1 (118, HE1), zinc (atom ZN), ligand residue SN1 (18 atoms; HN21 is the hydrogen on N2), 258 CA atoms.',
    'reference': 'Built by this reader from md/system/protein_md.pdb (258 CA, frame of 6RL9) and md/system/ligand_X.pdb (crystal pose, 11 heavy atoms). It is identical (0.0 A) to the CA and ligand atoms of reference.pdb of the built systems (check_starts.json).',
    'frames_and_time': 'All complete frames of solute_seg0001.dcd to solute_seg0020.dcd in order; production frame n (from 0) has the time (n + 1) * 0.01 ns. The DCD headers agree (first step, interval, step length).',
    'images': 'SPEC section 2: the zinc moved by the lattice vector that makes |Zn - NE2(His94)| smallest, then the whole ligand by the one that makes |N2 - Zn| smallest; search over 125 cells; lattice from the cell record of each frame (lengths and cosines). No frame of any of the eleven runs needed a move, and d1 and d2 equal their minimum-image values in every frame.',
    'fit': 'Least-squares superposition of the 258 CA atoms on the reference, equal weights, rotation without reflection (own Kabsch by SVD); the ligand is not fitted on itself.',
    'lig_rmsd': 'sqrt(mean over the 11 heavy atoms of |x_a - r_sigma(a)|^2) after the fit; smallest of the four labellings identity, ring_flip (C2<->C6, C3<->C5), oxygen_swap (O1<->O2), both. Mean over all frames and over the last half = frames N//2 to N-1 (5000 frames, 50.01 to 100.00 ns).',
    'd1_d2_state': 'd1 = |N1 - CZ(Phe131)|, d2 = |N1 - OD1(Asn67)|; C if d1 < 6.7 and d2 > 5.2; R if d1 > 6.7 and d2 < 5.2; I otherwise (a frame exactly on a midpoint would be I; there is none, and none within 0.0001 A unless the run entry says so).',
    'verdict': 'W = min(N, 2000) last frames; CONVERTED if 10*c_opposite >= 9*W, else STAYED if 10*c_start >= 9*W, else MIXED; integers.',
    'visits': 'Assigned state: the state of the frame if C or R, else that of the frame before; before the first frame the state of the start (C for X and N, R for R). Visit: maximal stretch of one assigned state; long: at least 100 frames. First passage: time of the first frame of the first long visit in the opposite state, in the production. Transitions: changes in (state of the start, states of the long visits in order).',
    'conversion_time': 'Plan, 7 Oct 2026, 12:10 AM IST, point 1: the solute frames of solute_eq04.dcd, solute_eq05.dcd, solute_eq06.dcd (40 + 40 + 50) joined before the production; same rules on the joined series; joined frame k (from 0) lies (k + 1) * 10 ps after the start of stage 4. Extension step eq03x1 of bonded_2e/R1 lies before the release and is not joined.',
    'hbond': 'Present if |N2 - OG1(Thr199)| <= 3.5 A and the angle N2-HN21-OG1 at the hydrogen >= 135 deg; fraction over all production frames.',
    'zinc': 'Means over all production frames of |Zn - N2|, |Zn - NE2(His94)|, |Zn - NE2(His96)|, |Zn - ND1(His119)|; fraction of solute frames with min(|Zn - O1|, |Zn - O2|) <= 2.6 A.',
    'n2_sign': 'Signed volume = u(N2->Zn) . [u(N2->HN21) x u(N2->S)] with unit vectors; it gives +0.7238, +0.7854 and -0.7102 on the three starts as built (plan: +0.724, +0.786, -0.710). Sign = +1 if the volume is above zero, -1 if below.',
    'choices_where_the_wording_left_one': [
        'Configuration at N2: section 9A of the specification was not read (only sections 2 to 6 and 9E were open to this reader). The sign is the plain sign of the volume, with no dead band. "First lasting change" is this reader\'s own definition, made like a first passage: assigned sign (a zero keeps the sign before; before the first frame the sign of the start: +1 for X and N, -1 for R), and the first frame of the first stretch of at least 100 frames with the other sign. It is given in the production alone and on the joined series of stages 4 to 6 plus production (time from release). The first single frame with the other sign is given beside it.',
        'Conversion time of a run that is not CONVERTED: the rule of 7 Oct gives a time to every run with a long visit in the opposite state and words it "converted during the equilibration, t ps after release" if that frame lies in the equilibration. For bonded_2e_n2fit/R1 the rule as worded gives that phrase (60 ps) although the verdict of the run is MIXED and its production holds 393 rotated frames; the entry of 6 Oct, 11:40 AM, point 2 ties the phrase to a production without a rotated frame. Reported by the rule as worded, with the verdict beside it; the two do not fit together for this run.',
        'The time of a conversion is the start of the first LONG visit, which can lie well after the first frame in the opposite state when the run goes back and forth first (bonded_2e/R1: first C frame 270 ps, time by the rule 540 ps; bonded_2s/R1: 70 ps against 890 ps). Both are in the run entries.',
        'Rule 3 counts a run as converted or stayed by its verdict (90 % of the last 20 ns). For R1 to R3 the entry of 6 Oct, 11:40 AM (no rotated frame in the production) gives the same count.',
        'Last half with N = 10000 is frames 5000 to 9999; with an odd N this reader would take N - N//2 frames.',
        '"Mean Zn-N distances" and the hydrogen-bond fraction are taken over the whole production (10000 frames), not over a window.',
        'The reference was built from md/system and not taken from reference.pdb; the two are identical.',
    ],
    'checks_of_this_reader': {
        'starts': 'check_starts.json: on the coordinates of solute.pdb the three starts give ligand RMSD 0.0005, 0.6106, 3.5732 A (starts.csv 0.0, 0.61, 3.573), d1 3.897, 3.356, 9.585 A and d2 7.744, 7.726, 2.597 A (starts.csv 3.897, 3.356, 9.585; 7.744, 7.725, 2.597), Zn-N2, Zn-O, N2-OG1 and the N2-H-OG1 angle within 0.001 A and 0.1 deg.',
        'mdanalysis': 'crosscheck_mda.json: on six files (X1 seg 1, R1 eq04, R3 seg 20, 2s/R1 seg 9, n2fit/R1 seg 15 and 20) the own DCD reader and MDAnalysis give the same coordinates (difference 0.0), d1 and d2 within 1e-6 A, ligand RMSD within 1e-13 A (fit by MDAnalysis rotation_matrix), and the same hydrogen-bond flag in every frame.',
    },
}

out = {
    'what': 'Blind third reading of the MD campaign (hCA II with the sulfanilamide anion, 6RL9): eleven runs of 100 ns',
    'written': now,
    'by': 'third reader (own code; scripts in this folder)',
    'runs': runs,
    'rule1': rule1,
    'rule3': rule3,
    'definitions_used': definitions,
}
with open(os.path.join(HERE, 'third_reading.json'), 'w') as f:
    json.dump(out, f, indent=1)
for k in runs:
    print(runs[k]['key_numbers_one_line'])
print(json.dumps(rule1, indent=1))
print(json.dumps(rule3, indent=1))
for k in runs:
    e = runs[k]['extras']
    print(k, e['states_frame_by_frame_stage4_5_6'] if runs[k]['start'] == 'R' else '', e['hbond_in_stage4_5_6'],
          e['visits_in_start_state_production'], e['visits_in_last_20ns_by_assigned_state'], e['frames_within_0.0001A_of_a_midpoint'],
          json.dumps(e['hydration_check_of_the_driver_from_manifest'])[:600])
print([(s['tag'], s['restraint_atoms']) for s in runs['bonded_2e/R1']['equilibration_steps']])
