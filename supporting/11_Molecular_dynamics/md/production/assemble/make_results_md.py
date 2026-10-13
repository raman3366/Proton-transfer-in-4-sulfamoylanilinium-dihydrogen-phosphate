#!/usr/bin/env python3
"""Assembles md/RESULTS.md and md/results/tables.md from the records of the campaign.

Every result is read from a file; the sentences of the plan and of the other records are cut out of those files;
the script stops if one of its own statements does not hold in the records. It reads only files that lie in the
repository: md/results/ (made by md/production/collect_results.py, conversion_time.py, make_result_figures.py and
site_checks.py; md/results/independent/ holds the physical review and the blind third reading), md/PLAN.md,
md/analysis/SPEC.md, md/GR.md, md/N2.md, md/PRODUCTION.md, md/protocol.json, md/system/audit/, the series of the
rehearsal run TN0 (md/test/rehearsal_2e/runs/TN0) and the records of the two gates on the terms at N2
(md/qm/gr/gr_gate.json, md/qm/n2/n2_gate_result.json).

Revised on 8 Oct 2026 after two audits of the first assembly and a physical review of the trajectories; what was
changed is listed in section 11 of the file it writes. Revised on 11 Oct 2026: the two places of the record of gate
G-R2 quoted in section 6 are named by their full keys, and section 10 says that the series and the final frames are
in git; listed in section 12 of the file it writes. Later on 11 Oct 2026: three dated notes on numbers of the
record that differ from the facts of the paper (the crystal value of the zinc out of the His119 plane in the plan's
sentence, section 5; the kinetic temperature at 2 fs, section 7; the count of third-reading values, section 9);
also listed in section 12.

    /root/cqc-sfa-tools/md/bin/python md/production/assemble/make_results_md.py
"""
import csv
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path

REPO = Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
R = REPO / 'md'
RES = R / 'results'
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))
J = lambda p: json.load(open(p))
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()

PROD = ['bonded_2e/' + r for r in ('X1', 'X2', 'X3', 'N1', 'N2', 'N3', 'R1', 'R2', 'R3')]
SENS = ['bonded_2s/R1', 'bonded_2e_n2fit/R1']
VAL = ['bonded_2e/X0', 'bonded_2s/X0']
ELEVEN = PROD + SENS
ALL = ELEVEN + VAL
XR, NR, RR = PROD[0:3], PROD[3:6], PROD[6:9]
RID = lambda lab: lab.split('/')[1]                      # run id inside its set
NAME = {lab: RID(lab) for lab in PROD}
NAME.update({'bonded_2s/R1': '2s/R1', 'bonded_2e_n2fit/R1': '2e_n2fit/R1', 'bonded_2e/X0': '2e/X0', 'bonded_2s/X0': '2s/X0'})

S1 = {lab: J(RES / 'runs' / lab / 'first_program_summary.json') for lab in ALL}
S2 = {lab: J(RES / 'runs' / lab / 'second_program_summary.json') for lab in ALL}
C1 = J(RES / 'campaign/first_program/campaign.json')
C2 = J(RES / 'campaign/second_program/campaign.json')
REP = J(RES / 'campaign/second_program/recheck_report.json')
CT = {r['label']: r for r in J(RES / 'conversion_times.json')['runs']}
HF = J(RES / 'run_health.json')
H = {r['run']: r for r in HF['runs']}
E = {r['run']: r for r in J(RES / 'run_events.json')['runs']}
SRC = J(RES / 'SOURCES.json')
FIG = J(RES / 'fig_sources.json')
GR = J(R / 'qm/gr/gr_gate.json')
N2G = J(R / 'qm/n2/n2_gate_result.json')
PLAN = re.sub(r'\s+', ' ', (R / 'PLAN.md').read_text(encoding='utf-8'))
SPEC = re.sub(r'\s+', ' ', (R / 'analysis/SPEC.md').read_text(encoding='utf-8'))
GATE_STATUS = int((RES / 'gate/gate_S6.exit_status').read_text().strip())
# added after the audits of 8 Oct 2026
SC = {r['run']: r for r in J(RES / 'site_checks_solute.json')['runs']}       # checks on the raw solute frames
SF = {r['run']: r for r in J(RES / 'site_checks_full.json')['runs']}         # checks on the raw all-atom frames
IND = J(RES / 'independent/SOURCES.json')
THIRD = J(RES / 'independent/third_reading/third_reading.json')
REVIEW_LINES = (RES / 'independent/physical_review/REVIEW.md').read_text(encoding='utf-8').splitlines()
REVIEW = re.sub(r'\s+', ' ', ' '.join(REVIEW_LINES))
RCMP = J(RES / 'independent/physical_review/compare.json')
GRMD = re.sub(r'\s+', ' ', (R / 'GR.md').read_text(encoding='utf-8'))
N2MD = re.sub(r'\s+', ' ', (R / 'N2.md').read_text(encoding='utf-8'))
PRODMD = re.sub(r'\s+', ' ', (R / 'PRODUCTION.md').read_text(encoding='utf-8'))
PROT = J(R / 'protocol.json')
HIS = J(R / 'system/audit/histidines.json')
WATCSV = {r['water']: r for r in csv.DictReader(open(R / 'system/audit/waters.csv', newline=''))}
TN0_DIR = R / 'test/rehearsal_2e/runs/TN0'
TN0 = list(csv.DictReader(open(TN0_DIR / 'analysis/timeseries_solute.csv', newline='')))
TN0_MAN = J(TN0_DIR / 'manifest.json')
FIRST_ASSEMBLY = '8 Oct 2026, 5:20 AM IST'          # the file that the two audits read (sha256 8a21f609f69f5ee7...)
# the records must be those of the programs as they stand now (a record left over from an older version stops the script)
for rec_, prog_ in ((RES / 'site_checks_solute.json', 'site_checks.py'), (RES / 'site_checks_full.json', 'site_checks.py'),
                    (RES / 'conversion_times.json', 'conversion_time.py'), (RES / 'fig_sources.json', 'make_result_figures.py')):
    assert J(rec_)['script']['sha256'] == sha(R / 'production' / prog_), 'not written by the present %s: %s' % (prog_, rec_)
GATE_LOG = (RES / 'gate/gate_S6.log').read_text()
GATE_TIMES = (RES / 'gate/gate_S6.times').read_text().splitlines()


def cut(text, a, b=None):
    """The passage of a text from the words a to the words b, as it stands there."""
    b = b or a
    i = text.find(a)
    assert i >= 0, 'not in the text: ' + a
    j = text.find(b, i)
    assert j >= 0, 'not in the text: ' + b
    return text[i:j + len(b)].replace('**', '')


q = lambda a, b=None: cut(PLAN, a, b)
RULE1_READ_A = '"ligand RMSD above 2 Å for the last half of the run" is read as'
PARTA_A = "The failure concerns the force field's whole description of the nitrogen"
for entry in ('5 Oct 2026, 8:40 PM IST', '6 Oct 2026, 1:03 AM IST', '6 Oct 2026, 3:02 AM IST', '6 Oct 2026, 11:40 AM IST',
              '6 Oct 2026, 12:19 PM IST', '6 Oct 2026, 1:09 PM IST', '6 Oct 2026, 1:34 PM IST', '6 Oct 2026, 8:41 PM IST',
              '6 Oct 2026, 9:15 PM IST', '6 Oct 2026, 10:50 PM IST', '7 Oct 2026, 12:10 AM IST'):
    assert entry in PLAN, entry


def f(x, n=2):
    return 'none' if x is None else ('%.' + str(n) + 'f') % x


def pc(count, n, d=1):
    """count / n in per cent with d decimals, rounded half up in whole numbers (so that 1487 of 2000 is 74.4)."""
    if isinstance(count, float) or isinstance(n, float):
        return ('%.' + str(d) + 'f') % (100.0 * count / n)
    k = (2 * count * 100 * 10 ** d + n) // (2 * n)
    return ('%d.%0' + str(d) + 'd') % (k // 10 ** d, k % 10 ** d) if d else str(k)


def th(n):
    return f'{n:,}'


def sg(x, n=1):
    """A signed number; a value that rounds to zero is written +0.0, never -0.0."""
    v = round(x, n)
    return ('%+.' + str(n) + 'f') % (0.0 if v == 0 else v)


# ------------------------------------------------------------------ the two programs give the same numbers
def same(a, b, tol=0.0011):
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(same(a[k], b[k], tol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        return a == b
    return abs(a - b) <= tol


for lab in ALL:
    for key in ('orientation', 'validation', 'hbond_first_passage', 'start_of_production', 'lig_rmsd_labelling_counts',
                'fractions', 'gate_GV'):
        a, b = S1[lab][key], S2[lab][key]
        if key == 'gate_GV':
            a = {k: v for k, v in a.items() if k != 'note'}
            b = {k: v for k, v in b.items() if k != 'note'}
            a['reported_zn_glu106'] = {k: v for k, v in a['reported_zn_glu106'].items() if k != 'note'}
            b['reported_zn_glu106'] = {k: v for k, v in b['reported_zn_glu106'].items() if k != 'note'}
        elif key in ('hbond_first_passage', 'start_of_production'):
            a = {k: v for k, v in a.items() if k != 'note'}
            b = {k: v for k, v in b.items() if k != 'note'}
        assert same(a, b, 0.0101), (lab, key)
    for col, st in S1[lab]['stats'].items():
        for k in ('mean', 'sd', 'min', 'max', 'mean_last_half', 'n'):
            assert same(st[k], S2[lab]['stats'][col][k], 0.0101), (lab, col, k)
    assert same(S1[lab]['table_row'], S2[lab]['table_row'], 0.0101), lab
    assert S1[lab]['script']['sha256'] == REP['primary_script_sha256'] and S2[lab]['script']['sha256'] == REP['script_sha256']
    assert '1.3' in S1[lab]['spec'] and '1.3' in S2[lab]['spec']
    assert S1[lab]['warnings'] == [] and S2[lab]['warnings'] == [], lab
assert C1['rules_per_model'] == C2['rules_per_model']

o = {lab: S1[lab]['orientation'] for lab in ALL}
st = {lab: S1[lab]['stats'] for lab in ALL}
fr = {lab: S1[lab]['fractions'] for lab in ALL}
tr = {lab: S1[lab]['table_row'] for lab in ALL}
run = {lab: S1[lab]['run'] for lab in ALL}
K = S1[PROD[0]]['constants']
M1, M2 = K['midpoint_d1_N1_F131CZ_A'], K['midpoint_d2_N1_N67OD1_A']
RULE1 = K['rule1_rmsd_A']
assert all(S1[lab]['constants'] == K for lab in ALL)
REF = H['bonded_2e/X1']['build']['geometry_as_built']          # start X as built = 6RL9
REFR = H['bonded_2e/R1']['build']['geometry_as_built']
REFN = H['bonded_2e/N1']['build']['geometry_as_built']
assert REF['lig_rmsd_symmetry_corrected_A'] == 0.0

# ------------------------------------------------------------------ statements that must hold
rules = C1['rules_per_model']
b2e = rules['bonded_2e']
assert GATE_STATUS == 0 and REP['agree'] is True and REP['n_runs_disagree'] == 0 and REP['n_runs_not_compared'] == 0
assert REP['n_threshold_items'] == 0 and REP['verdicts_that_differ_through_threshold_frames'] == []
assert REP['campaign']['agree'] is True and REP['campaign']['n_disagree'] == 0
assert 'Exit status 0 ' in GATE_LOG.splitlines()[-1]
assert [r['label'] for r in REP['runs']] == sorted(ALL + ['preflight_a', 'preflight_b'])
assert all(r['comparison']['agree'] and r['comparison']['n_disagree'] == 0 and r['comparison']['n_threshold'] == 0
           and r['comparison']['n_keys_one_side'] == 0 and r['comparison']['warning_count_equal'] for r in REP['runs'])
assert b2e['outcome'] == 'rotated_converts' and b2e['final'] is True and b2e['all_nine_complete'] is True
assert b2e['rule1'] == {'n_X_runs': 3, 'n_X_left_crystal': 0, 'failed': False}
assert b2e['n_converted'] == {'X': 0, 'N': 0, 'R': 3} and b2e['n_stayed'] == {'X': 3, 'N': 3, 'R': 0}
assert b2e['rule4_sensitivity_run_required'] is True
assert [o[l]['verdict'] for l in PROD] == ['STAYED'] * 6 + ['CONVERTED'] * 3
assert o['bonded_2s/R1']['verdict'] == 'CONVERTED' and o['bonded_2e_n2fit/R1']['verdict'] == 'MIXED'
assert all(run[l]['complete'] and run[l]['n_frames_solute'] == 10000 and run[l]['n_frames_full'] == 1000 for l in ELEVEN)
assert all(not o[l]['provisional'] and o[l]['window_frames'] == 2000 for l in ELEVEN)
assert all(CT[l]['same_frame'] for l in ALL)
assert all(CT[l]['conversion'] is None for l in XR + NR + VAL)
assert all(CT[l]['conversion'] and CT[l]['conversion']['in_equilibration'] for l in RR + SENS)
assert all(o[l]['count_R'] == 0 for l in RR)                    # "whose production holds no rotated frame"
assert all(o[l]['first_passage_frame'] == 0 for l in RR + SENS)
assert all(CT[l]['programs'][p]['n_released_equilibration_frames'] == 130 and
           CT[l]['programs'][p]['released_part_of_equilibration_ps'] == 1300.0
           for l in ALL for p in ('analysis', 'analysis_recheck'))
assert all(fr[l]['fifth_water']['count'] == 0 and fr[l]['fifth_proteinO']['count'] == 0 and fr[l]['fifth_E106OE']['count'] == 0
           for l in ALL)
assert all(S1[l]['gate_GV']['verdict'] is None for l in ELEVEN)       # a production run has the numbers of gate G-V without a verdict
assert S1['bonded_2e/X0']['gate_GV']['verdict'] == 'PASSED' and S1['bonded_2s/X0']['gate_GV']['verdict'] == 'FAILED'
assert all(H[l]['driver']['status'] == 'complete' and H[l]['driver']['all_steps_complete'] and H[l]['driver']['n_invocations'] == 1
           and H[l]['driver']['invocations'][0]['exit_code'] == 0 and not H[l]['driver']['test_only']
           and H[l]['driver']['equil_scale'] == 1.0 for l in ALL)
assert all(H[l]['trajectory_checksums']['n_dcd_files_read'] == H[l]['trajectory_checksums']['n_equal_to_the_driver_record']
           for l in ALL)
assert all(H[l]['driver_log']['lines_with_warning_error_or_failure'] == [] for l in ALL)
assert all(H[l]['production_state_log']['non_finite_values'] == 0 and H[l]['production_state_log']['thermostat_K'] == [300.0]
           and H[l]['production_state_log']['timestep_fs'] == [4.0] for l in ALL)
assert HF['final_verification']['verified'] is True and HF['final_verification']['driver_checksums']['bad'] == []
assert HF['final_verification']['pi_against_machine']['differ'] == []
asb = H['bonded_2e/R1']['build']['terms_at_N2']['sha256']
assert all(H[l]['build']['terms_at_N2']['sha256'] == asb for l in PROD + ['bonded_2s/R1'] + VAL)
assert H['bonded_2e_n2fit/R1']['build']['terms_at_N2']['sha256'] == N2G['candidates']['charges_as_built']['frcmod_sha256']
assert [H[l]['build']['parameter_set'] for l in SENS] == ['2s', '2e']
assert H['bonded_2s/R1']['driver']['run_seed'] == H['bonded_2e/R1']['driver']['run_seed'] == H['bonded_2e_n2fit/R1']['driver']['run_seed']
assert len({H[l]['driver']['run_seed'] for l in PROD}) == 9
assert N2G['accepted'] is False and GR['verdict_part_A']['passed'] is False
for l in ALL:
    for S in (S1, S2):
        ru = S[l]['run']
        want = None if l in VAL else True      # a validation run has one segment: nothing to compare
        assert ru['header_check_solute']['contiguous'] is want and ru['header_check_full']['contiguous'] is want
        assert len(ru['segments_solute']) == (1 if l in VAL else 20)
        assert ru['header_check_solute']['interval_matches'] and ru['header_check_full']['interval_matches']
        assert ru['segment_overlap_check']['suspected'] is False and ru['manifest_segment_check']['consistent'] is True
        assert ru['system_files_match_manifest'] == {'selections.json': True, 'reference.pdb': True}
        assert ru['imaging_solute']['frames_without_box'] == 0
x0fl = S1['bonded_2s/X0']['gate_GV']['GV3_zinc']['criteria']['fifth_ligand']['full_trajectory']
assert (x0fl['count'], x0fl['n']) == (3, 50) and x0fl['n'] // 20 == 2            # one frame more than allowed
s2s = E['bonded_2s/R1']['released_frames_n2_sign']
assert '-' not in s2s[s2s.index('+'):]

n_values_runs = sum(r['comparison']['n_values'] for r in REP['runs'])
n_values_13 = sum(r['comparison']['n_values'] for r in REP['runs'] if r['label'] in ALL)
n_values_camp = REP['campaign']['n_compared']
conv = {l: CT[l]['conversion'] for l in ALL}
cps = {l: conv[l]['time_from_release_ps'] for l in RR + SENS}


def conv_words(l):
    c = conv[l]
    if c is None:
        return 'none'
    assert c['in_equilibration']
    return 'converted during the equilibration, %.0f ps after release' % c['time_from_release_ps']


for l in RR + SENS:
    assert CT[l]['reported_as'] == conv_words(l)

# ------------------------------------------------------------------ added after the audits of 8 Oct 2026
FR = RR + SENS                                              # the five runs from start R
stg = {s['tag']: s for s in H[RR[0]]['equilibration_stages']}
t_sign = {l: E[l]['first_released_frame_n2_sign_not_of_start']['time_from_release_ps'] for l in FR}
t_notr = {l: E[l]['first_released_frame_not_in_state_of_start']['time_from_release_ps'] for l in FR}
t_firstc = {l: E[l]['first_released_frame_in_opposite_state']['time_from_release_ps'] for l in FR}
for l in FR:
    # the change of sign at N2 and the first frame outside the rotated-like state lie in stage 4 in every run from R
    assert E[l]['first_released_frame_n2_sign_not_of_start']['stage'] == 4 and E[l]['first_released_frame_not_in_state_of_start']['stage'] == 4, l
    assert t_sign[l] <= t_notr[l] <= stg['eq04']['length_ps'], l
    assert {s['tag']: s for s in H[l]['equilibration_stages']}['eq04']['restraint_k_kj_per_mol_nm2'] == stg['eq04']['restraint_k_kj_per_mol_nm2']
    assert H[l]['hydration_check']['reference_guard']['ca_fit_rmsd_of_system_inpcrd_on_reference_pdb_A'] == 0.0, l
assert stg['eq04']['restraint_atoms'] == 'protein_backbone' == stg['eq05']['restraint_atoms'] and stg['eq06']['restraint_atoms'] == 'none'
K4, K5 = f(stg['eq04']['restraint_k_kj_per_mol_nm2'], 0), f(stg['eq05']['restraint_k_kj_per_mol_nm2'], 0)
REFPOS = PROT['position_restraints']['reference_positions']
assert REFPOS.startswith('the coordinates of system.inpcrd')
RXM = N2G['candidates']['charges_as_built']['for_the_record']['Rx_against_X']['mismatch_before']
assert sum(o[l]['count_R'] for l in PROD) == 0
assert all(S1[l]['n2_configuration']['sign_of_start'] == 1 for l in XR + NR) and all(S1[l]['n2_configuration']['sign_of_start'] == -1 for l in FR)


def rng(vals, n=0, sep=' to '):
    """smallest and largest of a list, written as a range"""
    return f(min(vals), n) + sep + f(max(vals), n)


def lst(vals, n=0):
    """a list written as 'a, b and c'"""
    w = [f(v, n) if not isinstance(v, str) else v for v in vals]
    return ', '.join(w[:-1]) + ' and ' + w[-1]


GFR = N2G['candidates']['charges_as_built']['for_the_record']                 # old points of gate G-R, for the record
G_PL, G_HY, G_TN, G_RX = GFR['inversion_path_against_X']['cost_of_the_planar_nitrogen'], GFR['hydrogen_moved_into_the_plane'], GFR['turn_about_Zn_N2_against_X']['e_p030'], GFR['Rx_against_X']
# the checks on the raw frames (site_checks.py) against the two programs and against run_events.json
assert set(SC) == set(ALL) == set(SF)
for l in ALL:
    for S_ in (SC, SF):
        for sub in ('analysis', 'analysis_recheck'):
            assert not any(S_[l]['beside_the_two_programs'][sub]['frames_classed_differently'].values()), l
    assert SC[l]['n_production_frames'] == run[l]['n_frames_solute'] and SF[l]['n_production_frames'] == run[l]['n_frames_full']
for l in ELEVEN:
    e_, s_ = E[l], SC[l]
    a_, b_ = s_['rotated_like_frames'], e_['state_and_n2_sign']
    assert a_['production']['n'] == o[l]['count_R'] == b_['production']['R']
    assert a_['production']['n_with_V_positive'] == b_['production']['R_with_sign_positive']
    assert a_['production']['n_crystal_like_with_V_negative'] == b_['production']['C_with_sign_negative']
    assert a_['from_the_release']['n'] == b_['from_the_release']['R'] and a_['from_the_release']['n_with_V_positive'] == b_['from_the_release']['R_with_sign_positive']
    assert s_['nitrogen']['from_the_release']['n_changes_of_sign'] == e_['n2_volume']['from_the_release']['n_changes_of_sign']
    assert s_['nitrogen']['production']['n_abs_V_below_0.2'] == e_['n2_volume']['production']['n_abs_below_0.2']
    assert abs(s_['departures']['ang_S_N2_Zn_deg'] - st[l]['ang_S_N2_Zn_deg']['mean']) < 0.011 and abs(s_['departures']['d_S_N2_A'] - st[l]['d_S_N2_A']['mean']) < 0.00011
    assert abs(s_['ligand_against_the_docked_rotated_pose']['rmsd_to_the_crystal_pose_production_mean_A'] - st[l]['lig_rmsd_A']['mean']) < 0.00011
    assert SF[l]['n_water_O'] == H[l]['driver']['n_water_molecules']
for l in RR:
    assert SC[l]['change_at_N2']['first_frame_with_V_positive_ps_from_release'] == t_sign[l]
    assert E[l]['n2_volume']['from_the_release']['n_changes_of_sign'] == 1
# the table of the physical review (its section 2.2) against these checks: the same shares
RV_H64 = {}
for ln in REVIEW_LINES:
    m_ = re.fullmatch(r'\| (2e|2s|n2fit)/([XNR]\d) \| ([\d.]+) \| ([\d.]+) \| (\d+) \|', ln.strip())
    if m_:
        RV_H64[{'2e': 'bonded_2e/', '2s': 'bonded_2s/', 'n2fit': 'bonded_2e_n2fit/'}[m_.group(1)] + m_.group(2)] = (float(m_.group(3)), float(m_.group(4)), int(m_.group(5)))
assert set(RV_H64) == set(ELEVEN)
RVN = lambda l: {'bonded_2e': '2e', 'bonded_2s': '2s', 'bonded_2e_n2fit': 'n2fit'}[l.split('/')[0]] + '/' + RID(l)
rv = lambda key, l: RCMP[key]['X_run_means'][XR.index(l)] if l in XR else RCMP[key]['runs'][RVN(l)]['mean']
for l in ELEVEN:
    assert abs(RV_H64[l][0] - SC[l]['his64']['chi1_gauche_plus_share_production']) < 0.051, l
    assert abs(RV_H64[l][1] - SC[l]['his64']['deep_in_share_20_to_100_ns']['7.0']) < 0.051, l
    for w_ in ('423', '516', '553', '568'):
        assert abs(rv('W%s occupied (%%)' % w_, l) - SF[l]['crystal_water_places']['production_20_to_100_ns'][w_]) < 0.051, (l, w_)
    assert abs(rv('water O within 8 A of Zn (n)', l) - SF[l]['water_within_8_A_of_the_zinc']['mean_production_20_to_100_ns']) < 0.0051, l
# the blind third reading against the record of the first program, value by value
n_third = n_third_diff = 0


def third(a, b, tol=None):
    global n_third, n_third_diff
    n_third += 1
    if not ((a == b) if tol is None else (a is not None and b is not None and abs(a - b) <= tol)):
        n_third_diff += 1


assert set(THIRD['runs']) == set(ELEVEN)
for l in ELEVEN:
    T_ = THIRD['runs'][l]
    third(T_['verdict']['verdict'], o[l]['verdict'])
    for k_ in 'CRI':
        third(T_['states']['whole'][k_], o[l]['count_' + k_])
        third(T_['states']['last_20ns'][k_], o[l]['count_%s_window' % k_])
    third(T_['visits']['first_passage_frame'], o[l]['first_passage_frame'])
    third(T_['visits']['n_transitions'], o[l]['n_transitions'])
    third(T_['conversion']['report'] if T_['conversion']['has_conversion_time'] else None, CT[l]['reported_as'] if conv[l] else None)
    third(T_['lig_rmsd_A']['mean_whole'], st[l]['lig_rmsd_A']['mean'], 0.00011)
    third(T_['lig_rmsd_A']['mean_last_half'], st[l]['lig_rmsd_A']['mean_last_half'], 0.00011)
    third(T_['hbond_N2H_T199OG1']['count'], fr[l]['hb_N2H_T199OG1']['count'])
    third(T_['zinc']['sulfonyl_O_within_2.6A']['count'], fr[l]['fifth_sulfonylO']['count'])
    third(T_['n2_configuration']['frames_minus'], S1[l]['n2_configuration']['count_negative'])
    third(T_['zinc']['mean_Zn_N2_A'], st[l]['d_Zn_N2_A']['mean'], 0.00011)
assert n_third_diff == 0, n_third_diff
assert THIRD['rule3']['verdicts_of_the_nine_runs'] == {RID(l): o[l]['verdict'] for l in PROD}
# the rehearsal run TN0 (a test of the pipeline, no result): the same topology and start as N1
assert TN0_MAN['inputs']['system.prmtop']['sha256'] == H['bonded_2e/N1']['driver']['inputs_sha256']['system.prmtop']
assert TN0_MAN['inputs']['system.inpcrd']['sha256'] == H['bonded_2e/N1']['driver']['inputs_sha256']['system.inpcrd']
assert TN0_MAN['script']['sha256'] == H['bonded_2e/N1']['driver']['script_sha256'] and TN0_MAN['inputs']['protocol']['sha256'] == H['bonded_2e/N1']['driver']['protocol_sha256']
TN0_R = [r for r in TN0 if r['state'] == 'R']
k_tn0 = [k for k, r in enumerate(TN0) if r['state'] == 'R']
assert k_tn0 == list(range(k_tn0[0], k_tn0[-1] + 1)) and TN0[k_tn0[-1] + 1]['state'] == 'C'

L = []
A = L.append

# =============================================================================================== header
now = datetime.datetime.now(IST).strftime('%-d %b %Y, %-I:%M %p IST')
A('# Molecular dynamics of hCA II with the sulfanilamide anion on the zinc (PDB 6RL9): results of the campaign')
A('')
A(f'Stage S6 of `md/PLAN.md`. Assembled on {now} by `md/production/assemble/make_results_md.py` from the records '
  'named beside each number. No result is typed: every count, mean, time and verdict is read from a file. Constants of the protocol and of the definitions are typed in the script '
  '(the 10 ps between solute frames, the 2.6 Å about the zinc, the 90 % of the verdict and the like), and the script stops if a statement it makes does not hold in the records. '
  'The sentences of the plan and of the other records are cut out of those files by the script. '
  'Times are IST (UTC+5:30). This file is a record for the authors: no text of the paper was written or changed '
  '(plan, section 7: "No text is changed before the user has seen the results").')
A('')
A('The campaign: nine production runs of 100 ns with zinc set 2e (three from the crystal pose, start X; three from the '
  'docked near-native pose, start N; three from the docked ring-rotated pose, start R), two sensitivity runs of '
  'start R, and the two validation runs X0 of 5 ns. '
  f'The last run ended on {H["bonded_2e_n2fit/R1"]["driver"]["invocations"][0]["ended_ist"]}; '
  f'the guard\'s final check of {HF["final_verification"]["ist"]} reports `verified` '
  f'{str(HF["final_verification"]["verified"]).lower()} (section 8).')
A('')
A(f'**State of this record.** The two analysis programs agree on every run and on the campaign rules (gate of stage S6, exit status {GATE_STATUS}; section 9). '
  f'The first assembly of this file ({FIRST_ASSEMBLY}) was read by two audits, one of the numbers and one of the rules and the wording, and the trajectories were read by an independent '
  'physical review and by a blind third reading (`results/independent/`). Both audits found that the numbers, the verdicts and the application of rules 1 and 3 stand, and both asked for corrections: '
  'last digits in section 8, and sentences that said more than the data show or that put a limit more softly than the data give it. No finding of the review invalidates a run. '
  'Each point was checked again from the files before it was entered. The corrections are listed in section 11, and this file replaces the first assembly. '
  'Three limits that the first assembly put too softly are stated here before anything else:')
A('')
A(f'1. In R1, R2 and R3 the configuration at N2 changed {rng([t_sign[l] for l in RR])} ps after the release of the ligand, and the first frame outside the rotated-like state came at {rng([t_notr[l] for l in RR])} ps: '
  f'in all three runs in equilibration stage 4, while the protein backbone was held with {K4} kJ mol⁻¹ nm⁻² to the coordinates of the crystal structure. '
  f'No run of the nine holds the rotated configuration at N2 with that restraint lowered or removed, and the {th(sum(run[l]["n_frames_solute"] for l in PROD))} production frames of the nine runs hold '
  f'{sum(o[l]["count_R"] for l in PROD)} rotated-like frames (section 7, point 4).')
A(f'2. Of the known errors of the force-field terms at N2, one has the sign of the outcome: the force field puts the R-like placement {RXM:+.2f} kcal mol⁻¹ too high against the crystal-like one, '
  'and a placement drawn too high is left more readily and returned to less readily. How much of this acts in the protein is not known, its sign rests on the dispersion term of the reference, '
  'and most of it comes from charges fitted to the crystal pose alone, which all eleven long runs share (section 7, point 3).')
A(f'3. That no run returned to the rotated-like state is what the terms as built give, and they hold the nitrogen too stiffly by the numbers of gate G-R. With the candidate terms, which err the other way, '
  f'the one run from start R is {o[SENS[1]]["verdict"]} by the same rule (section 7, points 3 and 12).')
A('')

# =============================================================================================== 1
A('## 1. Outcome in short')
A('')
x_half = ', '.join(f'{NAME[l]} {f(S1[l]["validation"]["lig_rmsd_mean_last_half_A"])}' for l in XR)
A('| | Result | Read from |')
A('|---|---|---|')
A(f'| Gate of stage S6 (the two programs agree) | **Exit status {GATE_STATUS}.** The two programs agree on all {REP["n_runs"]} folders with a '
  f'trajectory below `/root/cqc-sfa-md/runs` (the {len(ALL)} runs of this file and the two pre-flight test folders) and on the campaign '
  f'rules: {th(n_values_runs)} values of the runs and {n_values_camp} of the campaign compared, {sum(r["comparison"]["n_disagree"] for r in REP["runs"])} '
  f'disagreements beyond tolerance, {REP["n_threshold_items"]} threshold items, '
  f'{sum(r["comparison"]["n_keys_one_side"] for r in REP["runs"])} keys in one program only (section 9) | '
  '`results/gate/gate_S6.exit_status`, `results/gate/gate_S6.log`, `results/campaign/second_program/recheck_report.json` |')
A(f'| Rule 1 (validation by the crystal starts) | Plan: "{q("If the crystal-start runs leave the crystal orientation", "reported as a failed validation.")}" '
  f'Read as fixed on 5 Oct 2026, 8:40 PM IST: {q(RULE1_READ_A, "over the last half of the run.")} '
  f'**Found: {b2e["rule1"]["n_X_left_crystal"]} of {b2e["rule1"]["n_X_runs"]} crystal-start runs left the crystal orientation** '
  f'(mean ligand RMSD of the last half: {x_half} Å; limit {f(RULE1, 1)} Å). Rule 1 is not failed (section 3) | '
  '`results/campaign/first_program/campaign.json` (`rules_per_model.bonded_2e.rule1`), the `summary.json` of X1 to X3 (`validation`) |')
A(f'| Rule 3 (the nine runs) | Plan, first case of rule 3: "{q("Rotated starts convert in two or three of three runs", "not stable outside the scoring function.")}" '
  f'**Found: R1, R2 and R3 CONVERTED; X1, X2, X3, N1, N2 and N3 STAYED** (verdicts over the last 20 ns). By the rule as worded '
  f'this is its first case (`{b2e["outcome"]}` in `campaign.json` of both programs, final: nine complete runs). '
  f'All three conversions lie before time zero, in the equilibration after the release of the ligand; in the plan\'s wording '
  f'(entry of 6 Oct 2026, 11:40 AM IST, point 2): R1 "{conv_words(RR[0])}", R2 "{conv_words(RR[1])}", R3 "{conv_words(RR[2])}". '
  f'These are the times of the plan\'s rule, the start of the first crystal-like visit of at least 1 ns. The sign at N2 turned positive at {lst([t_sign[l] for l in RR])} ps and the first frame that is not '
  f'rotated-like came at {lst([t_notr[l] for l in RR])} ps: in stage 4 in all three runs, with the protein backbone held at {K4} kJ mol⁻¹ nm⁻² (section 7, point 4). '
  f'The same entry: "{q("A rotated run whose production holds no rotated frame", "it counts as converted in rule 3.")}" '
  f'The production of R1, R2 and R3 holds {sum(o[l]["count_R"] for l in RR)} rotated-like frames of {th(sum(run[l]["n_frames_solute"] for l in RR))} | '
  '`results/campaign/first_program/campaign.json`, `results/campaign/second_program/campaign.json`, `results/conversion_times.json` |')
e2s, en2 = E['bonded_2s/R1'], E['bonded_2e_n2fit/R1']
A(f'| Sensitivity run with zinc set 2s (`bonded_2s/R1`; in no count of rule 3) | **CONVERTED**; "{conv_words("bonded_2s/R1")}" '
  f'(the time of the plan\'s rule, in stage {conv["bonded_2s/R1"]["stage"]}; the sign at N2 turned positive at {f(t_sign["bonded_2s/R1"], 0)} ps and the first frame that is not rotated-like came at {f(t_notr["bonded_2s/R1"], 0)} ps, in stage 4); '
  f'{o["bonded_2s/R1"]["count_R"]} rotated-like frame of {th(run["bonded_2s/R1"]["n_frames_solute"])} in the production. Set 2s failed gate G-V in its '
  'validation run by one all-atom frame (section 6) | `results/runs/bonded_2s/R1/first_program_summary.json`, `results/conversion_times.json` |')
w = o['bonded_2e_n2fit/R1']
A(f'| Sensitivity run with the candidate terms at N2 (`bonded_2e_n2fit/R1`; in no count of rule 3) | **MIXED**: in the last 20 ns '
  f'{pc(w["count_C_window"], w["window_frames"], 2)} % of the frames are crystal-like, {pc(w["count_R_window"], w["window_frames"], 2)} % rotated-like and '
  f'{pc(w["count_I_window"], w["window_frames"], 2)} % intermediate. The first crystal-like frame came {f(en2["first_released_frame_in_opposite_state"]["time_from_release_ps"], 0)} ps after the release and the first visit of at least 1 ns in the '
  f'crystal-like state began at {f(cps["bonded_2e_n2fit/R1"], 0)} ps; with the verdict MIXED this is no conversion. In the production the ring went back to the rotated-like state '
  f'{en2["production_visits_by_assigned_state"]["R"]["n_visits"]} times, never for 1 ns (longest visit '
  f'{f(en2["production_visits_by_assigned_state"]["R"]["longest_visit_frames"] * run["bonded_2e_n2fit/R1"]["solute_interval_ps"] / 1000.0)} ns). '
  f'A rotated-like frame is read from d1 and d2 only: of the {en2["state_and_n2_sign"]["production"]["R"]} rotated-like frames of the production {en2["state_and_n2_sign"]["production"]["R_with_sign_positive"]} have the sign at N2 of the crystal pose, '
  f'and the ligand is within 1.0 Å of the docked rotated pose in {SC["bonded_2e_n2fit/R1"]["ligand_against_the_docked_rotated_pose"]["production_frames_within_1.0_A"]} frames, so these counts are not counts of returns to the docked pose. '
  'These terms were not accepted by gate G-R2 and are used only to bracket (section 6) | '
  '`results/runs/bonded_2e_n2fit/R1/first_program_summary.json`, `results/run_events.json` |')
A(f'| What the outcome rests on | The force-field terms at the sulfonamide nitrogen N2 (plan, 6 Oct 2026, 3:02 AM IST, point 1: "{q("these terms and the H–N2–S–C4 series decide whether a run converts")}") failed two '
  'quantum-chemical gates (G-R and G-R2) and are used as built by the user\'s decision of 6 Oct 2026, 9:15 PM IST. By the plan\'s rule the three conversions of the nine runs lie '
  f'{rng([cps[l] for l in RR])} ps after the release of the ligand, before the production began; the change of sign at N2 lies at {rng([t_sign[l] for l in RR])} ps, in stage 4, with the protein backbone held. '
  f'(The sensitivity run with set 2s, in no count: {f(cps["bonded_2s/R1"], 0)} ps by the rule, change of sign at {f(t_sign["bonded_2s/R1"], 0)} ps.) '
  f'Of the known errors of these terms the R-like placement drawn {RXM:+.2f} kcal mol⁻¹ too high has the sign of the outcome, and with the candidate terms the one run from start R is {o[SENS[1]]["verdict"]}. '
  'Section 7 lists these and the other limits with the direction of each '
  f'error as the records give it. The plan allows "{q("no statement on rates or on the height of a barrier")}" | `md/PLAN.md`, `md/GR.md`, `md/N2.md`, `results/run_events.json` |')
A('')
A(f'The campaign table of the first program also prints an outcome line for the two sets with one run each (`bonded_2s`, `bonded_2e_n2fit`). '
  f'By the plan (entry of 7 Oct 2026, 12:10 AM IST, point 9) "{q("these lines are not outcomes of rule 3 and are not quoted")}"; they are not given here. '
  f'Of the sensitivity runs the plan says (entry of 6 Oct 2026, 9:15 PM IST, point 3): "{q("They are reported as sensitivity runs", "enter no count of rule 3.")}"')
A('')
assert FIG['conversion_marker_drawn_open_because_the_verdict_is_not_CONVERTED'] == ['bonded_2e_n2fit/R1']
assert FIG['stage_ends_ps_after_the_release'] == [round(stg['eq04']['length_ps']), round(stg['eq04']['length_ps'] + stg['eq05']['length_ps'])]
A('![d1 and d2 against time for the eleven runs](results/fig_orientation.png)')
A('')
A(f'*Figure 1.* The two distances of the orientation coordinate (rule 2 of the plan) against time: d1, aniline N1 to Phe131 CZ, and d2, '
  f'aniline N1 to Asn67 OD1, with their midpoints {f(M1, 1)} and {f(M2, 1)} Å. Grey strip: the {f(CT[PROD[0]]["programs"]["analysis"]["released_part_of_equilibration_ps"] / 1000.0, 1)} ns of the equilibration after the release of the '
  f'ligand (stages 4 to 6), time axis stretched. In the production the pale line is every frame (10 ps) and the strong line the mean over '
  f'{FIG["running_mean_frames"]} frames. A filled triangle marks the conversion by the plan\'s rule. In the panel of the candidate terms the triangle is open: the verdict of that run is MIXED, and the triangle marks the start of its first '
  f'crystal-like visit of at least 1 ns, not a conversion. The two thin broken lines in the grey strip are the ends of stages 4 and 5, {FIG["stage_ends_ps_after_the_release"][0]} and {FIG["stage_ends_ps_after_the_release"][1]} ps after the release: '
  f'the protein backbone is held in stage 4 ({K4} kJ mol⁻¹ nm⁻²) and in stage 5 ({K5}) and free in stage 6. File `results/fig_orientation.png` (and `.pdf`), '
  'drawn by `md/production/make_result_figures.py`; sources and checks in `results/fig_sources.json`. The ligand RMSD in the same layout: '
  '`results/fig_ligand_rmsd.png`.')
A('')

# =============================================================================================== 2
A('## 2. The eleven runs')
A('')
A('| Run | Start | Zinc set; terms at N2 | Length (ns); frames, solute / all-atom | Verdict (last 20 ns) | Last 20 ns: crystal-like / rotated-like / intermediate (%) | '
  'Whole run: crystal-like / rotated-like / intermediate (%) | First passage in the production (ns) | Conversion, time from the release (plan\'s rule) | '
  'Transitions | Mean d1; d2 of the last half (Å) |')
A('|---|---|---|---|---|---|---|---|---|---|---|')


def set_words(l):
    b = H[l]['build']
    n2 = 'as built' if b['terms_at_N2']['sha256'] == asb else 'candidate not accepted by G-R2'
    return f'{b["parameter_set"]}; {n2}'


def row2(l):
    x = o[l]
    W, N = x['window_frames'], S1[l]['run']['n_frames_solute']
    c = conv[l]
    if c is None:
        ctext = 'none'
    elif x['verdict'] == 'CONVERTED':
        ctext = f'{f(c["time_from_release_ps"], 0)} ps, in the equilibration (stage {c["stage"]})'
    else:       # the rule gives a time although the verdict of the run is not CONVERTED
        ctext = (f'{f(c["time_from_release_ps"], 0)} ps (stage {c["stage"]}): start of the first crystal-like visit of at least 1 ns; '
                 f'no conversion, the verdict is {x["verdict"]}')
    return (f'| {NAME[l]} | {run[l]["start"]} | {set_words(l)} | {f(run[l]["seen_ns"])}; {th(N)} / {th(run[l]["n_frames_full"])} | {x["verdict"]} | '
            f'{pc(x["count_C_window"], W, 2)} / {pc(x["count_R_window"], W, 2)} / {pc(x["count_I_window"], W, 2)} | '
            f'{pc(x["count_C"], N, 2)} / {pc(x["count_R"], N, 2)} / {pc(x["count_I"], N, 2)} | {f(x["first_passage_ns"])} | {ctext} | {x["n_transitions"]} | '
            f'{f(st[l]["d1_N1_F131CZ_A"]["mean_last_half"])}; {f(st[l]["d2_N1_N67OD1_A"]["mean_last_half"])} |')


for l in PROD:
    A(row2(l))
A('| *sensitivity runs, in no count of rule 3* | | | | | | | | | | |')
for l in SENS:
    A(row2(l))
A('')
A(f'Starts as built (`md/system/starts.csv`; here from the build record in each manifest, `results/run_health.json`): X d1 {f(REF["d1_N1_Phe131CZ_A"])}, d2 {f(REF["d2_N1_Asn67OD1_A"])} Å; '
  f'N {f(REFN["d1_N1_Phe131CZ_A"])}, {f(REFN["d2_N1_Asn67OD1_A"])} Å; R {f(REFR["d1_N1_Phe131CZ_A"])}, {f(REFR["d2_N1_Asn67OD1_A"])} Å. '
  f'A frame is crystal-like if d1 < {f(M1, 1)} Å and d2 > {f(M2, 1)} Å, rotated-like if d1 > {f(M1, 1)} Å and d2 < {f(M2, 1)} Å, intermediate otherwise '
  f'(`md/analysis/SPEC.md`, section 4). Verdict, in the plan\'s words (5 Oct 2026, 8:40 PM IST): "{q("A run has converted if at least 90 %", "and is mixed otherwise.")}"')
A('')
A('How to read three columns, and what the time of the rule is not:')
A('')
A(f'- **First passage in the production** is counted by the two programs from the start of the production only. The value {f(o[RR[0]]["first_passage_ns"])} ns of the five runs from '
  'start R is the time of the first production frame: the run was already in a visit of at least 1 ns in the crystal-like state when the production began.')
A(f'- **Conversion, time from the release** is the number the plan asks for (entry of 7 Oct 2026, 12:10 AM IST, point 1): "{q("The time of the conversion is the time of the first frame of the first long visit", "counted from the start of stage 4.")}" '
  f'It is formed on the frames of equilibration stages 4, 5 and 6 joined with the production, from the `state` written by both programs; the two give the same frame in all {len(ALL)} runs '
  '(`md/production/conversion_time.py`, tested by `md/production/test_conversion_time.py`; `results/conversion_times.json`). The plan: the time describes one run; '
  f'"{q("it is not a rate")}".')
A(f'- For `2e_n2fit/R1` the rule gives a time although the verdict of the run is MIXED: {f(cps["bonded_2e_n2fit/R1"], 0)} ps is when its first visit of at least 1 ns in the crystal-like state began, '
  'not the time of a conversion. The rule was applied as worded. Its one transition counts visits of at least 1 ns only; the shorter visits to the rotated-like state are in section 6.')
A('- **The time of the rule is not the time at which the ligand left the rotated-like state.** First frame with the sign at N2 positive and first frame that is not rotated-like, from the release: '
  + '; '.join(f'{NAME[l]} {f(t_sign[l], 0)} and {f(t_notr[l], 0)} ps' for l in FR)
  + f'. In every one of the five runs both lie in stage 4 (0 to {f(stg["eq04"]["length_ps"], 0)} ps after the release; protein backbone held at {K4} kJ mol⁻¹ nm⁻²). '
  f'The time of the rule is later where the ring went back for single frames: R1 {f(cps[RR[0]], 0)} ps (stage {conv[RR[0]]["stage"]}), 2s/R1 {f(cps[SENS[0]], 0)} ps (stage {conv[SENS[0]]["stage"]}). Sections 4 and 6 give the frames.')
A('')
A('The sign at N2, used here and below, is the sign of the signed volume of Zn, H and S about N2 (specification, section 9A): positive in the crystal pose and in start N, negative in start R '
  f'(`n2_configuration.sign_of_start` of each `summary.json`). Its limit, in the words of the specification (section 15): "{cut(SPEC, "The sign at N2 (section 9A) says on which side", "No rule uses it.")}"')
A('')

# =============================================================================================== 3
A('## 3. Validation by the crystal starts (rule 1)')
A('')
A('| Run | Ligand RMSD to the crystal pose, mean ± SD (Å) | Mean of the last half (Å) | Largest (Å) | Hydrogen bond N2–H···OG1(Thr199), % of frames | '
  'C4–S–N2, mean ± SD (°); mean minus 6RL9 | τ_fold, mean ± SD (°); mean minus 6RL9 | Last half above ' + f(RULE1, 1) + ' Å (rule 1) |')
A('|---|---|---|---|---|---|---|---|')
A(f'| 6RL9 (start X as built) | {f(REF["lig_rmsd_symmetry_corrected_A"])} | | | N2–H···OG1 {f(REF["N2_HN21_Thr199OG1_deg"], 1)}°, N2···OG1 {f(REF["N2_Thr199OG1_A"])} Å | '
  f'{f(REF["C4_S_N2_deg"], 1)} | {f(REF["tau_fold_deg"], 1)} | |')


def row3(l):
    v = S1[l]['validation']
    g4 = S1[l]['gate_GV']['GV4_ligand_internal']['criteria']
    a, t = g4['ang_C4_S_N2_deg'], g4['tau_fold_deg']
    applies = 'no' if not v['rule1_left_crystal'] else '**yes**'
    if not v['rule1_applies_to_this_run']:
        applies += ' (not counted by rule 1)'
    return (f'| {NAME[l]} | {f(v["lig_rmsd_mean_A"])} ± {f(v["lig_rmsd_sd_A"])} | {f(v["lig_rmsd_mean_last_half_A"])} | {f(st[l]["lig_rmsd_A"]["max"])} | '
            f'{pc(fr[l]["hb_N2H_T199OG1"]["count"], fr[l]["hb_N2H_T199OG1"]["n"], 2)} | {f(a["mean_deg"], 1)} ± {f(a["sd_deg"], 1)}; {sg(a["deviation_deg"])} | '
            f'{f(t["mean_deg"], 1)} ± {f(t["sd_deg"], 1)}; {sg(t["deviation_deg"])} | {applies} |')


for l in XR + NR:
    A(row3(l))
A(row3('bonded_2e/X0').replace('| 2e/X0 |', '| 2e/X0 (validation run, 5 ns; gate G-V ' + S1['bonded_2e/X0']['gate_GV']['verdict'] + ') |'))
A('')
assert all(S1[l]['validation']['rule1_applies_to_this_run'] for l in XR) and not any(S1[l]['validation']['rule1_applies_to_this_run'] for l in NR + VAL)
assert all(S1[l]['lig_rmsd_labelling_counts']['identity'] == run[l]['n_frames_solute'] for l in PROD + VAL)
assert all(S1[l]['n2_configuration']['count_negative'] == 0 and S1[l]['n2_configuration']['count_zero'] == 0 for l in XR + NR + VAL)
A(f'**Rule 1 as worded:** "{q("If the crystal-start runs leave the crystal orientation", "reported as a failed validation.")}" '
  f'**Verdict:** the mean ligand RMSD of the last half is {x_half} Å in the three crystal-start runs, against the limit of {f(RULE1, 1)} Å; '
  f'{b2e["rule1"]["n_X_left_crystal"]} of {b2e["rule1"]["n_X_runs"]} left the crystal orientation, so rule 1 is not failed. The three near-native starts are shown beside; rule 1 does not count them.')
A('')
xn = XR + NR
A(f'In the six runs the largest RMSD of a single frame is {f(max(st[l]["lig_rmsd_A"]["max"] for l in xn))} Å ({NAME[max(xn, key=lambda l: st[l]["lig_rmsd_A"]["max"])]}), '
  f'the hydrogen bond to Thr199 is present in {pc(min(fr[l]["hb_N2H_T199OG1"]["count"] for l in xn), 10000, 2)} to {pc(max(fr[l]["hb_N2H_T199OG1"]["count"] for l in xn), 10000, 2)} % of the frames, '
  f'no frame is rotated-like ({sum(o[l]["count_R"] for l in xn)} of {th(sum(run[l]["n_frames_solute"] for l in xn))}), the sign at N2 is that of the start in every frame, and '
  f'the smallest RMSD is always found with the atom labels as they are (labelling `identity` in all frames). '
  f'The mean C4–S–N2 angle lies {sg(min(S1[l]["gate_GV"]["GV4_ligand_internal"]["criteria"]["ang_C4_S_N2_deg"]["deviation_deg"] for l in xn))} to '
  f'{sg(max(S1[l]["gate_GV"]["GV4_ligand_internal"]["criteria"]["ang_C4_S_N2_deg"]["deviation_deg"] for l in xn))}° from the crystal value and the mean τ_fold '
  f'{sg(min(S1[l]["gate_GV"]["GV4_ligand_internal"]["criteria"]["tau_fold_deg"]["deviation_deg"] for l in xn))} to '
  f'{sg(max(S1[l]["gate_GV"]["GV4_ligand_internal"]["criteria"]["tau_fold_deg"]["deviation_deg"] for l in xn))}°. '
  'The numbers of gate G-V are written for every run, without a verdict: for a production run they are '
  'not a second validation (`md/analysis/SPEC.md`, sections 11 and 15), and the gate was decided on run X0 alone (`md/PRODUCTION.md`, section 5).')
A('')
vfr = [E[l]['n2_volume']['from_the_release'] for l in xn]
vpr = [E[l]['n2_volume']['production'] for l in xn]
assert all(v['n_changes_of_sign'] == 0 and v['min'] > 0 for v in vfr)
tn_d1, tn_d2, tn_v = [float(r['d1_N1_F131CZ_A']) for r in TN0_R], [float(r['d2_N1_N67OD1_A']) for r in TN0_R], [float(r['n2_volume']) for r in TN0_R]
A(f'**What "stayed" shows and what it does not.** In the six runs the signed volume at N2 lies between {sg(min(v["min"] for v in vfr), 2)} and {sg(max(v["max"] for v in vfr), 2)} from the release of the ligand to the end, '
  f'and its size is below 0.2 in {min(v["n_abs_below_0.2"] for v in vpr)} to {max(v["n_abs_below_0.2"] for v in vpr)} of {th(vpr[0]["n"])} production frames per run: the nitrogen stays pyramidal throughout. '
  f'By the numbers of gate G-R the terms as built overstate the cost of flattening it (the hydrogen moved into the plane {f(G_HY["as_built"])} against {f(G_HY["QM"])} kcal mol⁻¹ in quantum chemistry, '
  f'the planar nitrogen on the constructed path {f(G_PL["as_built"])} against {f(G_PL["QM"])}; section 7, point 1). So "{sum(o[l]["count_R"] for l in xn)} of {th(sum(run[l]["n_frames_solute"] for l in xn))}" is what these terms give in these six runs. '
  'It is not evidence about the molecule, and rule 1 tests whether the crystal pose is kept, not whether it is kept too firmly (section 7, point 12). '
  'It is also a statement about these six runs only. In the rehearsal of 6 Oct 2026, which was a test of the pipeline, is no result and enters no count, the run of '
  f'{f(TN0_MAN["production"]["target_ns"], 0)} ns from start N (TN0; the same topology, start coordinates, protocol and driver as N1: `system.prmtop` sha256 `{TN0_MAN["inputs"]["system.prmtop"]["sha256"][:12]}…` in both manifests) '
  f'had {len(TN0_R)} rotated-like frames of {len(TN0)}, at {f(float(TN0_R[0]["time_ns"]))} to {f(float(TN0_R[-1]["time_ns"]))} ns (d1 {", ".join(f(v) for v in tn_d1)} Å; d2 {", ".join(f(v) for v in tn_d2)} Å; '
  f'signed volume at N2 {", ".join(sg(v, 2) for v in tn_v)}, so negative in {sum(1 for v in tn_v if v < 0)} of them), and was crystal-like again in the next frame '
  '(`md/test/rehearsal_2e/runs/TN0/analysis/timeseries_solute.csv`; `md/test/REHEARSAL_2e.md`).')
A('')
A(f'N3 has the largest mean of the last half ({f(S1["bonded_2e/N3"]["validation"]["lig_rmsd_mean_last_half_A"])} Å): between 60 and 80 ns its mean RMSD is '
  f'{f(E["bonded_2e/N3"]["production_blocks_of_20_ns"][3]["lig_rmsd_mean_A"])} Å, against {f(min(b["lig_rmsd_mean_A"] for i, b in enumerate(E["bonded_2e/N3"]["production_blocks_of_20_ns"]) if i != 3))} to '
  f'{f(max(b["lig_rmsd_mean_A"] for i, b in enumerate(E["bonded_2e/N3"]["production_blocks_of_20_ns"]) if i != 3))} Å in its other four blocks of 20 ns (`results/run_events.json`; figure `results/fig_ligand_rmsd.png`).')
A('')

# =============================================================================================== 4
A('## 4. The rotated starts, run by run')
A('')
stg = {s['tag']: s for s in H[RR[0]]['equilibration_stages']}
for l in RR + SENS:
    for t in ('eq04', 'eq05', 'eq06'):
        s_ = {s['tag']: s for s in H[l]['equilibration_stages']}[t]
        assert (s_['restraint_atoms'], s_['restraint_k_kj_per_mol_nm2'], s_['length_ps']) == (stg[t]['restraint_atoms'], stg[t]['restraint_k_kj_per_mol_nm2'], stg[t]['length_ps'])
assert stg['eq03']['restraint_atoms'] == 'solute_heavy' and stg['eq06']['restraint_atoms'] == 'none'
A(f'The ligand is held in stages 1 to 3 of the equilibration (heavy atoms of the solute, {f(stg["eq03"]["restraint_k_kj_per_mol_nm2"], 0)} kJ mol⁻¹ nm⁻² in stage 3) and is released at the start of stage 4. '
  f'After the release: stage 4, {f(stg["eq04"]["length_ps"], 0)} ps with the protein backbone held at {f(stg["eq04"]["restraint_k_kj_per_mol_nm2"], 0)} kJ mol⁻¹ nm⁻²; '
  f'stage 5, {f(stg["eq05"]["length_ps"], 0)} ps at {f(stg["eq05"]["restraint_k_kj_per_mol_nm2"], 0)} kJ mol⁻¹ nm⁻²; stage 6, {f(stg["eq06"]["length_ps"], 0)} ps with nothing held; then the production '
  '(each manifest, `equilibration.stages[].settings`; `results/run_health.json`). A solute frame is written every 10 ps, an all-atom frame every 100 ps. '
  'No restraint and no collective variable acts on the orientation of the ring (plan, section 4).')
A('')


def hyd_words(l):
    hc = H[l]['hydration_check']
    parts = []
    for c in hc['counts']:
        parts.append(f'end of `{c["step"]}`: {c["n_water_O"]} (nearest {f(c["nearest_water_O_A"], 3)} Å), "{c["decision"]}"')
    return '; '.join(parts)


def strip_pos(l, cond):
    s_, n_, h_ = E[l]['released_frames_states'], E[l]['released_frames_n2_sign'], E[l]['released_frames_hb_N2H_T199OG1']
    for k in range(len(s_)):
        if cond(s_[k], n_[k], h_[k]):
            return k
    return None


def ev(l, key):
    v = E[l][key]
    return 'none' if v is None else f'{f(v["time_from_release_ps"], 0)} ps'


def last_full_before(l):
    """The last all-atom frame of the released equilibration that lies before the first frame with the sign at N2
    changed (all-atom frames are 100 ps apart)."""
    rel_ps = SC[l]['n_released_equilibration_frames'] * run[l]['solute_interval_ps']
    xs = [dict(x, ps=round(x['time_ns'] * 1000.0 + rel_ps, 3)) for x in E[l]['water_in_the_place_of_the_ring_equilibration'] if x['stage'] >= 4]
    assert [x['ps'] for x in xs] == [100.0 * (k + 1) for k in range(len(xs))], l
    before = [x for x in xs if x['ps'] < t_sign[l]]
    return before[-1] if before else None


rows4 = [
    ('Held part of the equilibration: solute frames; rotated-like; sign at N2 negative; hydrogen bond to Thr199 present',
     lambda l: f'{E[l]["held_part"]["n"]}; {E[l]["held_part"]["R"]}; {E[l]["held_part"]["n2_negative"]}; {E[l]["held_part"]["hb_N2H_T199OG1"]}'),
    ('Hydration check of the driver (water oxygens within ' + f(H[RR[0]]['hydration_check']['settings']['water_cutoff_A'], 1) + ' Å of the crystal positions of the six ring carbons)', hyd_words),
    ('Extension steps of stage 3 (500 ps each); length of the equilibration (ps)',
     lambda l: f'{H[l]["hydration_check"]["n_extensions"]}; {f(H[l]["driver"]["equilibration"]["completed_ps"], 0)}'),
    ('Place of the ring without water when the ligand was released (driver)', lambda l: 'no' if H[l]['hydration_check']['site_without_water_when_the_solute_was_released'] is False else 'yes'),
    ('First released frame: state; signed volume at N2; ligand RMSD (Å)',
     lambda l: f'{E[l]["first_released_frame"]["state"]}; {E[l]["first_released_frame"]["n2_volume"]:+.2f}; {f(E[l]["first_released_frame"]["lig_rmsd_A"])}'),
    ('First frame with the sign at N2 positive (time from the release; state of that frame)',
     lambda l: f'{ev(l, "first_released_frame_n2_sign_not_of_start")}; {E[l]["first_released_frame_n2_sign_not_of_start"]["state"]}'),
    ('First frame that is not rotated-like', lambda l: f'{ev(l, "first_released_frame_not_in_state_of_start")}; {E[l]["first_released_frame_not_in_state_of_start"]["state"]}'),
    ('Stage of these two frames; protein backbone held in that stage at (kJ mol⁻¹ nm⁻²)',
     lambda l: f'{E[l]["first_released_frame_n2_sign_not_of_start"]["stage"]} and {E[l]["first_released_frame_not_in_state_of_start"]["stage"]}; {K4}'),
    ('Torsion NE2(His94)–Zn–N2–S, mean: released frames before the change of sign; production (°)',
     lambda l: f'{sg(SC[l]["change_at_N2"]["torsion_S_before_deg"], 0)}; {sg(SC[l]["change_at_N2"]["torsion_S_production_deg"], 0)}'),
    ('Torsion NE2(His94)–Zn–N2–H, mean: the same (°)',
     lambda l: f'{sg(SC[l]["change_at_N2"]["torsion_H_before_deg"], 0)}; {sg(SC[l]["change_at_N2"]["torsion_H_production_deg"], 0)}'),
    ('Last all-atom frame before the change of sign: time from the release; water oxygens in the place of the ring; nearest (Å)',
     lambda l: f'{f(last_full_before(l)["ps"], 0)} ps; {last_full_before(l)["n_water_ring_place"]}; {f(last_full_before(l)["d_water_ring_place_min_A"], 3)}'),
    ('First crystal-like frame', lambda l: ev(l, 'first_released_frame_in_opposite_state')),
    ('Last rotated-like frame of the equilibration', lambda l: ev(l, 'last_released_frame_in_state_of_start')),
    ('**Conversion by the plan\'s rule** (first frame of the first visit of at least 1 ns in the crystal-like state); stage',
     lambda l: f'**{f(cps[l], 0)} ps**; stage {conv[l]["stage"]}'),
    ('Released frames before the conversion: number; rotated-like / crystal-like / intermediate; hydrogen bond present',
     lambda l: f'{E[l]["released_before_the_conversion"]["n"]}; {E[l]["released_before_the_conversion"]["R"]} / {E[l]["released_before_the_conversion"]["C"]} / {E[l]["released_before_the_conversion"]["I"]}; {E[l]["released_before_the_conversion"]["hb_N2H_T199OG1"]}'),
    ('Released frames from the conversion on: number; rotated-like; hydrogen bond present',
     lambda l: f'{E[l]["released_from_the_conversion_on"]["n"]}; {E[l]["released_from_the_conversion_on"]["R"]}; {E[l]["released_from_the_conversion_on"]["hb_N2H_T199OG1"]}'),
    ('Water oxygens in the place of the ring in the 13 all-atom frames after the release (every 100 ps)',
     lambda l: ' '.join(str(x['n_water_ring_place']) for x in E[l]['water_in_the_place_of_the_ring_equilibration'] if x['stage'] >= 4)),
    ('Production: crystal-like / rotated-like / intermediate frames',
     lambda l: f'{th(o[l]["count_C"])} / {o[l]["count_R"]} / {o[l]["count_I"]}'),
    ('Production: frames with the sign at N2 negative; hydrogen bond to Thr199, % of frames',
     lambda l: f'{S1[l]["n2_configuration"]["count_negative"]}; {pc(fr[l]["hb_N2H_T199OG1"]["count"], fr[l]["hb_N2H_T199OG1"]["n"], 2)}'),
    ('Production: all-atom frames with a water oxygen in the place of the ring (of ' + th(run[RR[0]]['n_frames_full']) + ')',
     lambda l: str(E[l]['water_in_the_place_of_the_ring_production']['frames_with_one_or_more'])),
    ('Production: ligand RMSD to the crystal pose, mean; mean of the last half; mean of the last 20 ns (Å)',
     lambda l: f'{f(st[l]["lig_rmsd_A"]["mean"])}; {f(st[l]["lig_rmsd_A"]["mean_last_half"])}; {f(E[l]["production_blocks_of_20_ns"][-1]["lig_rmsd_mean_A"])}'),
    ('Production: d1; d2, mean of the last half (Å)',
     lambda l: f'{f(st[l]["d1_N1_F131CZ_A"]["mean_last_half"])}; {f(st[l]["d2_N1_N67OD1_A"]["mean_last_half"])}'),
    ('Production: C4–S–N2, mean ± SD (°); τ_fold, mean ± SD (°)',
     lambda l: f'{f(gv4(l, "ang_C4_S_N2_deg")["mean_deg"], 1)} ± {f(gv4(l, "ang_C4_S_N2_deg")["sd_deg"], 1)}; {f(gv4(l, "tau_fold_deg")["mean_deg"], 1)} ± {f(gv4(l, "tau_fold_deg")["sd_deg"], 1)}'),
]
gv4 = lambda l, k: S1[l]['gate_GV']['GV4_ligand_internal']['criteria'][k]
A('| | R1 | R2 | R3 |')
A('|---|---|---|---|')
for name, fn in rows4:
    A(f'| {name} | ' + ' | '.join(fn(l) for l in RR) + ' |')
A('')
A('Sources: `results/run_events.json` (counted from the equilibration and production series of both programs, which give the same counts), `results/run_health.json` '
  '(the driver\'s manifest), `results/conversion_times.json`, the `summary.json` of each run, and for the two torsions `results/site_checks_solute.json` (own reading of the raw frames). The equilibration series themselves are in '
  '`results/runs/bonded_2e/R1/first_program_timeseries_equil_solute.csv` (and R2, R3).')
A('')
A('The 130 solute frames after the release, one character per frame (10 ps): `R` rotated-like, `C` crystal-like, `I` intermediate; small letters where the sign at N2 is negative '
  '(the configuration of start R), capitals where it is positive (the configuration of the crystal pose). `|` marks the ends of stages 4 and 5.')
A('')
A('```')


def strip_text(l):
    s_, n_ = E[l]['released_frames_states'], E[l]['released_frames_n2_sign']
    t = ''.join(c if g == '+' else c.lower() for c, g in zip(s_, n_))
    return t[:40] + '|' + t[40:80] + '|' + t[80:]


for l in RR:
    A(f'{NAME[l]:<4} {strip_text(l)}')
A('```')
A('')
# statements about the three runs, each checked
for l in RR:
    k_sign = strip_pos(l, lambda s_, n_, h_: n_ == '+')
    k_opp = strip_pos(l, lambda s_, n_, h_: s_ == 'C')
    k_not = strip_pos(l, lambda s_, n_, h_: s_ != 'R')
    assert k_sign is not None and k_sign <= k_not <= k_opp and k_opp - k_sign <= 2, l
    assert set(E[l]['released_frames_hb_N2H_T199OG1'][max(0, k_sign - 1):k_opp + 1]) == {'1'}, l
    assert '-' not in E[l]['released_frames_n2_sign'][k_sign:], l
    assert E[l]['held_part']['R'] == E[l]['held_part']['n'] == E[l]['held_part']['n2_negative'], l
    assert S1[l]['n2_configuration']['count_negative'] == 0
A('What the series show, run by run:')
A('')
for l in RR:
    e = E[l]
    hc = H[l]['hydration_check']
    k_sign = e['first_released_frame_n2_sign_not_of_start']
    k_opp = e['first_released_frame_in_opposite_state']
    last_r = e['last_released_frame_in_state_of_start']
    stage_k = stg['eq0%d' % conv[l]['stage']]['restraint_k_kj_per_mol_nm2']
    txt = f'- **{NAME[l]}.** '
    if hc['n_extensions']:
        c0, c1 = hc['counts'][0], hc['counts'][-1]
        txt += (f'At the end of stage 3 the driver found {c0["n_water_O"]} water oxygens in the place of the crystal ring (nearest {f(c0["nearest_water_O_A"], 3)} Å) and extended the held part by '
                f'{hc["n_extensions"]} step of {f(H[l]["hydration_check"]["counts"][1]["time_end_ps"] - H[l]["hydration_check"]["counts"][0]["time_end_ps"], 0)} ps, after which it found {c1["n_water_O"]} (nearest {f(c1["nearest_water_O_A"], 3)} Å) and went on. ')
    else:
        c1 = hc['counts'][-1]
        txt += f'At the end of stage 3 the driver found {c1["n_water_O"]} water oxygen{"s" if c1["n_water_O"] != 1 else ""} in the place of the crystal ring (nearest {f(c1["nearest_water_O_A"], 3)} Å); no extension step. '
    txt += (f'The ligand stayed rotated-like with the sign at N2 negative for {k_sign["released_frame"]} frames after the release. The sign at N2 turned positive at {f(k_sign["time_from_release_ps"], 0)} ps '
            f'and the first crystal-like frame followed at {f(k_opp["time_from_release_ps"], 0)} ps, with the hydrogen bond N2–H···OG1(Thr199) present in every frame from the one before the change of sign to the first crystal-like frame; '
            f'the sign stayed positive in every later frame of the equilibration and in all {th(run[l]["n_frames_solute"])} frames of the production. ')
    seg = e['released_frames_states'][k_opp['released_frame'] + 1: conv[l]['joined_index']]
    nR = seg.count('R')
    if nR:
        assert last_r['time_from_release_ps'] > k_opp['time_from_release_ps'] and last_r['n2_sign'] == 1
        f530 = SC[l]['change_at_N2']['frame_at_530_ps']
        assert last_r['time_from_release_ps'] == 530.0 and f530['state'] == 'R' and f530['V'] > 0 and nR == 1
        txt += (f'In the {len(seg)} frames after it the ring was crystal-like in {seg.count("C")}, intermediate in {seg.count("I")} and rotated-like in {nR} '
                f'(at {f(last_r["time_from_release_ps"], 0)} ps, with the sign at N2 positive), so the first visit of at least 1 ns in the crystal-like state, and with it the conversion by the plan\'s rule, begins at {f(cps[l], 0)} ps. '
                f'In that one frame d1 is {f(f530["d1_A"])} Å and d2 {f(f530["d2_A"])} Å with the signed volume at N2 {sg(f530["V"], 2)} and the sum of the three angles at N2 {f(f530["sum_of_angles_deg"], 1)}°: '
                'the ring swung back for one frame without a change of the configuration at N2. The rule was applied as worded; the time of the rule is '
                f'{f(cps[l] - k_sign["time_from_release_ps"], 0)} ps later than the change of sign. ')
    else:
        assert last_r['time_from_release_ps'] < k_opp['time_from_release_ps'] and len(seg) == 0
        txt += f'No rotated-like frame came after it, so the conversion by the plan\'s rule is at {f(cps[l], 0)} ps. '
    wrel = [x for x in e['water_in_the_place_of_the_ring_equilibration'] if x['stage'] >= 4]
    t_conv_ns = conv[l]['time_ns_in_file']
    wb = [x for x in wrel if x['time_ns'] < t_conv_ns - 1e-9]
    wa = [x for x in wrel if x['time_ns'] >= t_conv_ns - 1e-9]
    assert len(wrel) == 13
    if conv[l]['stage'] == 4:
        txt += f'The change of sign, the first frame outside the rotated-like state and the time of the rule all lie in stage 4 (protein backbone held at {K4} kJ mol⁻¹ nm⁻²). '
    else:
        txt += (f'The time of the rule lies in stage {conv[l]["stage"]} (protein backbone held at {f(stage_k, 0)} kJ mol⁻¹ nm⁻²); the change of sign and the first frame outside the rotated-like state '
                f'lie in stage 4 (held at {K4} kJ mol⁻¹ nm⁻²). ')
    txt += (f'All-atom frames after the release with at least one water oxygen in the place of the ring: {sum(1 for x in wb if x["n_water_ring_place"] >= 1)} of {len(wb)} before the conversion, '
            f'{sum(1 for x in wa if x["n_water_ring_place"] >= 1)} of {len(wa)} after it; the last {next(i for i, x in enumerate(reversed(wrel)) if x["n_water_ring_place"] >= 1)} before the production have none. ')
    lf = last_full_before(l)
    wseq = [x for x in wrel if round(x['time_ns'] * 1000.0 + 1300.0, 3) <= lf['ps']]
    if lf['n_water_ring_place'] == 0:
        assert l == RR[1] and hc['counts'][-1]['n_water_O'] >= 1
        txt += (f'In the last all-atom frame before the change of sign ({f(lf["ps"], 0)} ps after the release) the place of the crystal ring held no water oxygen (nearest {f(lf["d_water_ring_place_min_A"], 3)} Å), after '
                f'{lst([str(x["n_water_ring_place"]) for x in wseq[:-1]])} in the frames before it. The hydration check of the plan was met at the end of the held part ({hc["counts"][-1]["n_water_O"]} water oxygens), so the rule was kept; '
                f'but the concern behind the check (plan, 6 Oct 2026, 1:03 AM IST, point 3: "{q("An empty pocket beside the ligand could draw the ring back", "look easier than it is.")}") applies to this moment of R2. '
                'All-atom frames are 100 ps apart, so nothing finer can be said.')
    else:
        txt += (f'In the last all-atom frame before the change of sign ({f(lf["ps"], 0)} ps after the release) the place of the crystal ring held {lf["n_water_ring_place"]} water oxygen{"s" if lf["n_water_ring_place"] != 1 else ""} '
                f'(nearest {f(lf["d_water_ring_place_min_A"], 3)} Å).')
    A(txt)
A('')
x_lo, x_hi = min(st[l]['lig_rmsd_A']['mean_last_half'] for l in XR), max(st[l]['lig_rmsd_A']['mean_last_half'] for l in XR)
n_lo, n_hi = min(st[l]['lig_rmsd_A']['mean_last_half'] for l in NR), max(st[l]['lig_rmsd_A']['mean_last_half'] for l in NR)
r_lo, r_hi = min(st[l]['lig_rmsd_A']['mean_last_half'] for l in RR), max(st[l]['lig_rmsd_A']['mean_last_half'] for l in RR)
chg = {l: SC[l]['change_at_N2'] for l in RR}
A('**The event is a change of the configuration at N2 with the hydrogen bond kept.** In each of the three runs the signed volume at N2 changes its sign once, between two frames 10 ps apart, and not again '
  f'({max(E[l]["n2_volume"]["from_the_release"]["n_changes_of_sign"] for l in RR)} change of sign from the release to the end of the production in each of R1, R2 and R3). '
  'The sulfur changes its place about the Zn–N2 bond and the hydrogen keeps its side: the torsion NE2(His94)–Zn–N2–S goes from '
  f'{lst([sg(chg[l]["torsion_S_before_deg"], 0) for l in RR])}° (mean of the released frames before the change) to {lst([sg(chg[l]["torsion_S_production_deg"], 0) for l in RR])}° (production), '
  f'the torsion NE2(His94)–Zn–N2–H from {lst([sg(chg[l]["torsion_H_before_deg"], 0) for l in RR])}° to {lst([sg(chg[l]["torsion_H_production_deg"], 0) for l in RR])}°, '
  'and N2–H···OG1(Thr199) is present in every frame from the one before the change of sign to the first crystal-like frame (`results/site_checks_solute.json`; `results/tables.md`, table D). '
  f'The plan named this path before any run (6 Oct 2026, 3:02 AM IST, point 1): "{q("Going from R to X needs an inversion at N2 or the loss of that hydrogen bond.")}" '
  'The event therefore runs through the force-field terms at N2 that failed gates G-R and G-R2 (section 7, points 1 to 3).')
A('')
c4 = lambda ls: rng([gv4(l, 'ang_C4_S_N2_deg')['mean_deg'] for l in ls], 1)
tf = lambda ls: rng([gv4(l, 'tau_fold_deg')['mean_deg'] for l in ls], 1)
A(f'**Does a converted run end like the crystal-start runs?** In the pose of the ligand, its labelling, its internal angles and the hydrogen bond to Thr199, yes: mean ligand RMSD of the last half {f(r_lo)} to {f(r_hi)} Å in R1 to R3, against {f(x_lo)} to {f(x_hi)} Å in X1 to X3 and '
  f'{f(n_lo)} to {f(n_hi)} Å in N1 to N3; the smallest RMSD is found with the atom labels as they are in all frames (no turned ring, no exchanged oxygens); '
  f'the hydrogen bond to Thr199 is present in {pc(min(fr[l]["hb_N2H_T199OG1"]["count"] for l in RR), 10000, 2)} to {pc(max(fr[l]["hb_N2H_T199OG1"]["count"] for l in RR), 10000, 2)} % of the production frames '
  f'({pc(min(fr[l]["hb_N2H_T199OG1"]["count"] for l in xn), 10000, 2)} to {pc(max(fr[l]["hb_N2H_T199OG1"]["count"] for l in xn), 10000, 2)} % in the six X and N runs); '
  f'the mean C4–S–N2 angle is {c4(RR)}° in R1 to R3 against {c4(xn)}° in the six X and N runs and the mean τ_fold {tf(RR)}° against {tf(xn)}° (table above and section 3); '
  'the four Zn–N distances and the N–Zn–N angles are those of the other runs of set 2e (section 5; `results/tables.md`, table A). '
  f'The site as a whole was not compared by the two programs; the specification (section 15): "{cut(SPEC, "The state of His64 and the charge of Glu106 are fixed by the system, not examined here.")}" It does differ, in His64 and in the water of the site.')
A('')
# ---- His64 and the water of the site (check on the raw frames; found by the physical review)
h64 = {x['label']: x for x in HIS if x['residue'] == 'HIS64'}
assert h64['A']['taken'] and not h64['B']['taken'] and not h64['C']['taken']
w516 = re.search(r'HIS64\.ND1\[B\]\(occ [\d.]+\):([\d.]+)', WATCSV['516']['polar_partners']).group(1)
dp = lambda l, c='7.0': SC[l]['his64']['deep_in_share_20_to_100_ns'][c]
wl = lambda l, w: SF[l]['crystal_water_places']['production_20_to_100_ns'][w]
deep_mean = [SC[l]['his64']['Zn_NE2_mean_while_deep_in_7.0_A'] for l in RR]
A('**What differs: His64 and the water of the site.** The physical review found this (`results/independent/physical_review/REVIEW.md`, sections 2.2 and 4.3); the numbers here are those of a check on the raw frames made for this record '
  '(`md/production/site_checks.py`), and the table of the review has the same shares. His64 was built in its outer position (conformer A of 6RL9, occupancy '
  f'{f(h64["A"]["occupancy"])}; Zn to NE2 {f(h64["A"]["Zn_NE2"], 1)} Å). In the runs it moves between outside and inside over tens of ns. Deep inside, its NE2 lies {rng(deep_mean, 1)} Å from the zinc (means of R1 to R3 while deep in), '
  f'nearer than the nearer ring nitrogen of either inner conformer of 6RL9 ({f(min(h64["B"]["Zn_ND1"], h64["B"]["Zn_NE2"]), 1)} Å in conformer B, occupancy {f(h64["B"]["occupancy"])}; '
  f'{f(min(h64["C"]["Zn_ND1"], h64["C"]["Zn_NE2"]))} Å in C, occupancy {f(h64["C"]["occupancy"])}; `md/system/audit/histidines.json`), and the place of crystal water 516, which in 6RL9 is fully occupied and '
  f'{w516} Å from ND1 of conformer B (`md/system/audit/waters.csv`), is then mostly empty. '
  f'In R1, R2 and R3 His64 is deep in for {lst([dp(l) for l in RR], 1)} % of the production from 20 to 100 ns, and the place of water 516 is occupied in {lst([wl(l, "516") for l in RR], 1)} % of the all-atom frames of that time. '
  f'In the six runs from X and N the shares are {rng([dp(l) for l in xn], 1)} % and {rng([wl(l, "516") for l in xn], 1)} %. The place of water 423 is occupied in {rng([wl(l, "423") for l in RR], 1)} % against {rng([wl(l, "423") for l in xn], 1)} %, '
  f'those of waters 553 and 568 alike in the nine runs ({rng([wl(l, "553") for l in PROD], 1)} and {rng([wl(l, "568") for l in PROD], 1)} %). With His64 goes the number of water oxygens within 8 Å of the zinc: '
  f'{rng([SF[l]["water_within_8_A_of_the_zinc"]["mean_production_20_to_100_ns"] for l in RR], 1)} in R1 to R3 against {rng([SF[l]["water_within_8_A_of_the_zinc"]["mean_production_20_to_100_ns"] for l in xn], 1)} in the six '
  f'(in the first 10 ns of the production: {rng([SF[l]["water_within_8_A_of_the_zinc"]["mean_first_10_ns"] for l in RR], 1)} in R1 to R3).')
A('')
A('| Run | His64 deep in, % of the production from 20 to 100 ns: Zn–NE2 at most 6.5 / 7.0 / 7.5 Å, each with χ1 between 0 and 120° | χ1 of His64 between 0 and 120°, % of the production | '
  'Places of crystal waters 516 / 423 / 553 / 568 occupied, % of the all-atom frames from 20 to 100 ns | Water oxygens within 8 Å of the zinc, mean, 20 to 100 ns | A water donating a hydrogen bond to N2, % of the all-atom frames of the production |')
A('|---|---|---|---|---|---|')
for l in ELEVEN:
    if l == SENS[0]:
        A('| *sensitivity runs* | | | | | |')
    A(f'| {NAME[l]} | {f(dp(l, "6.5"), 1)} / {f(dp(l), 1)} / {f(dp(l, "7.5"), 1)} | {f(SC[l]["his64"]["chi1_gauche_plus_share_production"], 1)} | '
      f'{f(wl(l, "516"), 1)} / {f(wl(l, "423"), 1)} / {f(wl(l, "553"), 1)} / {f(wl(l, "568"), 1)} | {f(SF[l]["water_within_8_A_of_the_zinc"]["mean_production_20_to_100_ns"], 2)} | {f(SF[l]["water_donating_to_N2"]["share_production"], 1)} |')
A('')
A('A place of a crystal water counts as occupied when a water oxygen lies within 1.4 Å of it after the fit of the Cα atoms; a water donates to N2 when its oxygen is within 3.5 Å of N2 with O–H···N2 at least 135° '
  '(`results/site_checks_solute.json`, `results/site_checks_full.json`). Four things limit what this table says. '
  'The limit of 7.0 Å for "deep in" was chosen by the review after looking at the data; the shares with 6.5 and 7.5 Å are beside it. '
  f'His64 is a slow coordinate that 100 ns do not sample: X3 and N3 are deep in for {f(dp(XR[2]), 1)} and {f(dp(NR[2]), 1)} % of the same time. '
  'Whether three of three runs from start R against none to about half in the six others is chance in a slow coordinate or a memory of the rotated start cannot be told from these runs, and the review found no cause. '
  f'And the two sensitivity runs from start R do not show it ({f(dp(SENS[0]), 1)} and {f(dp(SENS[1]), 1)} %). '
  'The water structure of the site is therefore not called reproduced or converged in this record, and "ends like the crystal-start runs" holds for the ligand only.')
A('')
A(f'The two programs give the hydrogen bond "before and after the first passage" for the production only (`md/analysis/SPEC.md`, section 9B); there the passage is frame 0 in all three runs, so "before" holds no frame and "after" is the whole run '
  f'({", ".join(pc(S1[l]["hbond_first_passage"]["count_after"], S1[l]["hbond_first_passage"]["n_after"], 2) for l in RR)} %). The counts of the table before and after the conversion in the equilibration were made by `md/production/collect_results.py` from the column `hb_N2H_T199OG1` of the equilibration series.')
A('')
A(f'No statement on rates or barriers is made (plan, entry of 6 Oct 2026, 9:15 PM IST, point 4: "{q("no statement on rates or on the height of a barrier")}"). '
  'Frames 10 ps apart show the order of events only to that resolution.')
A('')

# =============================================================================================== 5
A('## 5. The zinc site in all runs')
A('')
ANG = ['ang_N2_Zn_H94NE2_deg', 'ang_N2_Zn_H96NE2_deg', 'ang_N2_Zn_H119ND1_deg', 'ang_H94NE2_Zn_H96NE2_deg', 'ang_H94NE2_Zn_H119ND1_deg', 'ang_H96NE2_Zn_H119ND1_deg']
ANGN = {'ang_N2_Zn_H94NE2_deg': 'N2–Zn–NE2(His94)', 'ang_N2_Zn_H96NE2_deg': 'N2–Zn–NE2(His96)', 'ang_N2_Zn_H119ND1_deg': 'N2–Zn–ND1(His119)',
        'ang_H94NE2_Zn_H96NE2_deg': 'NE2(His94)–Zn–NE2(His96)', 'ang_H94NE2_Zn_H119ND1_deg': 'NE2(His94)–Zn–ND1(His119)', 'ang_H96NE2_Zn_H119ND1_deg': 'NE2(His96)–Zn–ND1(His119)'}
A('| Run | Zn–N2 | Zn–NE2(His94) | Zn–NE2(His96) | Zn–ND1(His119) | Largest departure of a mean N–Zn–N angle from 6RL9 (°) | Nearer sulfonyl oxygen to Zn, mean (Å) | '
  'Sulfonyl oxygen within 2.6 Å, solute frames (10 ps) | Any oxygen within 2.6 Å, all-atom frames (100 ps) | of these water / protein oxygen |')
A('|---|---|---|---|---|---|---|---|---|---|')
A(f'| 6RL9 (start X as built) | {f(REF["Zn_N2_A"], 3)} | {f(REF["Zn_His94NE2_A"], 3)} | {f(REF["Zn_His96NE2_A"], 3)} | {f(REF["Zn_His119ND1_A"], 3)} | | {f(min(REF["Zn_O1_A"], REF["Zn_O2_A"]), 3)} | | | |')


def row5(l):
    cr = S1[l]['gate_GV']['GV3_zinc']['criteria']
    worst = max(ANG, key=lambda a: abs(cr[a]['deviation_deg']))
    fl = cr['fifth_ligand']
    s_, f_ = fl['solute_trajectory'], fl['full_trajectory']
    return (f'| {NAME[l]} | {f(st[l]["d_Zn_N2_A"]["mean"], 3)} ± {f(st[l]["d_Zn_N2_A"]["sd"], 3)} | {f(st[l]["d_Zn_H94NE2_A"]["mean"], 3)} ± {f(st[l]["d_Zn_H94NE2_A"]["sd"], 3)} | '
            f'{f(st[l]["d_Zn_H96NE2_A"]["mean"], 3)} ± {f(st[l]["d_Zn_H96NE2_A"]["sd"], 3)} | {f(st[l]["d_Zn_H119ND1_A"]["mean"], 3)} ± {f(st[l]["d_Zn_H119ND1_A"]["sd"], 3)} | '
            f'{sg(cr[worst]["deviation_deg"])} ({ANGN[worst]}) | {f(st[l]["d_Zn_Osulfonyl_min_A"]["mean"], 3)} | {s_["count"]} of {th(s_["n"])} ({pc(s_["count"], s_["n"], 2)} %) | '
            f'{f_["count"]} of {th(f_["n"])} ({pc(f_["count"], f_["n"], 1)} %) | {fr[l]["fifth_water"]["count"]} / {fr[l]["fifth_proteinO"]["count"]} |')


for l in ELEVEN + VAL:
    A(row5(l))
A('')
A('Distances are means ± SD over the solute frames (Å). The six angles of every run are in `results/tables.md`, table A. Source: `stats`, `fractions` and `gate_GV` of each `summary.json`.')
A('')
e_runs = PROD + ['bonded_2e/X0']
dev_all = [abs(S1[l]['gate_GV']['GV3_zinc']['criteria'][a]['deviation_deg']) for l in e_runs for a in ANG]
A(f'**For set 2e the Zn–N distances and the N–Zn–N angles are put in, not reproduced.** The plan (entry of 7 Oct 2026, 12:10 AM IST, point 3): "{q("In set 2e the equilibrium values of the four Zn–N bonds", "independent reproduction of the zinc geometry.")}" '
  f'In the nine runs and X0 of set 2e the mean Zn–N2 distance is {f(min(st[l]["d_Zn_N2_A"]["mean"] for l in e_runs), 3)} to {f(max(st[l]["d_Zn_N2_A"]["mean"] for l in e_runs), 3)} Å and no mean angle is more than {f(max(dev_all), 1)}° from its value in 6RL9. '
  'What is not put in: the place of the sulfonyl oxygens, of water and of protein oxygens about the zinc.')
A('')
A(f'**An oxygen within 2.6 Å of the zinc.** In no all-atom frame of any run is a water oxygen or a protein oxygen within 2.6 Å of the zinc, and neither carboxylate oxygen of Glu106 comes within 2.6 Å in any solute frame '
  f'(smallest Zn–OE1 distance over all runs {f(min(S1[l]["gate_GV"]["reported_zn_glu106"]["d_Zn_E106OE1_A"]["min_A"] for l in ALL))} Å). Every count of the table is a sulfonyl oxygen of the ligand. '
  f'In the nine runs of set 2e the nearer sulfonyl oxygen is within 2.6 Å in {pc(min(fr[l]["fifth_sulfonylO"]["count"] for l in PROD), 10000, 2)} to {pc(max(fr[l]["fifth_sulfonylO"]["count"] for l in PROD), 10000, 2)} % of the solute frames; '
  f'in the two sensitivity runs in {pc(fr[SENS[0]]["fifth_sulfonylO"]["count"], 10000, 2)} % (set 2s) and {pc(fr[SENS[1]]["fifth_sulfonylO"]["count"], 10000, 2)} % (candidate terms at N2). '
  f'The plan on this count: "{q("In the bonded model the sulfonyl oxygens and the zinc are three bonds apart", "the scaled 1–4 one")}".')
A('')
# added on 11 Oct 2026: the crystal value in the plan's sentence (0.17 Å) is that of His94; the value of His119 is
# formed here from the coordinates the systems are built from (appended to the paragraph, so that every line of the
# file keeps its number)
def _zn_out_of_ring_planes():
    import numpy as np
    def atoms(path):
        return [(ln[12:16].strip(), int(ln[22:26]), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
                for ln in open(path) if ln.startswith(('ATOM', 'HETATM'))]
    zn = [a for a in atoms(R / 'system/zinc.pdb') if a[0].upper() == 'ZN']
    assert len(zn) == 1
    prot, out = atoms(R / 'system/protein_md.pdb'), {}
    for res in (94, 96, 119):
        ring = np.array([a[2] for a in prot if a[1] == res and a[0] in ('CG', 'ND1', 'CD2', 'CE1', 'NE2')])
        assert ring.shape == (5, 3), res
        c = ring.mean(0)
        out[res] = abs(float(np.dot(zn[0][2] - c, np.linalg.svd(ring - c)[2][2])))
    return out


OOP = _zn_out_of_ring_planes()
assert f(OOP[119]) == '0.15' and f(OOP[94]) == '0.17' and f(OOP[96]) == '0.06', OOP
assert '(0.17 Å in 6RL9)' in q("the zinc lies 0.5 Å out of the plane of the His119 ring", "O1–S–C4 101.6° (106.6°).")
NOTE_OOP = (' **Added on 11 Oct 2026.** The crystal value in the sentence of the plan, 0.17 Å, is that of the His94 ring and not of His119. In the coordinates of 6RL9 that the systems are built from '
            '(`md/system/protein_md.pdb`, `md/system/zinc.pdb`; least-squares plane of the five ring atoms, formed by this script) the zinc lies '
            f'{f(OOP[119])} Å out of the plane of the His119 ring ({f(OOP[119], 3)} Å), {f(OOP[94])} Å out of that of His94 and {f(OOP[96])} Å out of that of His96; the physical review has 0.15 for His119, as quoted above. '
            f'The right crystal value for this departure is therefore {f(OOP[119])} Å. The sentence of the plan is quoted as it stands; the means of the runs and the reading of the departure do not change.')
A(f'Not tested by any gate and named by the plan (entry of 7 Oct 2026, 12:10 AM IST, point 8) for set 2e: "{q("the zinc lies 0.5 Å out of the plane of the His119 ring", "O1–S–C4 101.6° (106.6°).")}" '
  f'Of these the analysis follows S–N2–Zn: its mean is {f(min(st[l]["ang_S_N2_Zn_deg"]["mean"] for l in PROD), 1)} to {f(max(st[l]["ang_S_N2_Zn_deg"]["mean"] for l in PROD), 1)}° in the nine runs '
  f'({f(REF["S_N2_Zn_deg"], 1)}° in 6RL9), and the mean S–N2 distance {f(min(st[l]["d_S_N2_A"]["mean"] for l in PROD), 4)} to {f(max(st[l]["d_S_N2_A"]["mean"] for l in PROD), 4)} Å ({f(REF["S_N2_A"], 4)} Å in 6RL9). '
  'The other three are not followed by the two programs and were not measured again by them over the nine runs. The check on the raw frames made for this record gives, as means over the production of each of the nine runs: '
  f'the zinc {rng([SC[l]["departures"]["zinc_out_of_the_plane_of_the_His119_ring_A"] for l in PROD], 2)} Å out of the plane of the His119 ring, S–C4 {rng([SC[l]["departures"]["d_S_C4_A"] for l in PROD], 3)} Å, '
  f'O1–S–C4 {rng([SC[l]["departures"]["ang_O1_S_C4_deg"] for l in PROD], 1)}° (`results/site_checks_solute.json`). The physical review, with its own code, has the first two: "{cut(REVIEW, "zinc 0.46 to 0.53 Å out of the His119 plane (0.15)")}" and '
  f'"{cut(REVIEW, "S-C4 1.821 to 1.822 Å (1.770)")}" (its section 6, point 7; its reference values are its own measurement of start X). So the four departures that the plan named from X0 and the first 20 ns of X1 hold over the nine runs of 100 ns.'
  + NOTE_OOP)
A('')
A(f'**Gate G-Z** (stage S2) is no result of this stage and is not evaluated by the two programs as a gate of its own (`md/analysis/SPEC.md`, section 15). It passed for both zinc sets on the first 100 ps of X0, frames in which the zinc site is still held; '
  f'the plan (7 Oct 2026, 12:10 AM IST, point 5): "{q("Both sets pass under five other readings tried by the checker", "it is reported as a check on held frames.")}" Its numbers, with those of the free stages beside them, are in `md/PRODUCTION.md`, section 4.')
A('')

# =============================================================================================== 6
A('## 6. The two sensitivity runs')
A('')
x0s, x0e = S1['bonded_2s/X0'], S1['bonded_2e/X0']
fl_s = x0s['gate_GV']['GV3_zinc']['criteria']['fifth_ligand']
fl_e = x0e['gate_GV']['GV3_zinc']['criteria']['fifth_ligand']
assert not fl_s['passed'] and fl_e['passed'] and x0s['gate_GV']['GV1_lig_rmsd']['passed'] and x0s['gate_GV']['GV2_hbond']['passed'] and x0s['gate_GV']['GV4_ligand_internal']['passed']
assert all(v['passed'] for k, v in x0s['gate_GV']['GV3_zinc']['criteria'].items() if k != 'fifth_ligand')
A('Both are single runs of start R of 100 ns with the run id R1, so with the seeds of R1 of the primary set '
  f'(run seed {H["bonded_2e/R1"]["driver"]["run_seed"]} in all three manifests). They follow rule 4 of the plan as changed on 6 Oct 2026 '
  f'(3:02 AM IST, point 4: "{q("The model with a free Zn²⁺ and four distance restraints does not keep the site", "Rule 4 uses the second bonded set of (3).")}"; '
  f'9:15 PM IST, point 3: "{q("with the other zinc set; and with the best candidate of the refit", "which was not accepted by gate G-R2 and is used only to bracket")}"). '
  'They are in no count of rule 3.')
A('')
bn2s, bn2f = H['bonded_2s/R1']['build'], H['bonded_2e_n2fit/R1']['build']
note2s = bn2s['status_note'].split(' Whether 2s or 2e')[0]
assert note2s.startswith('PRODUCTION system') and 'Zinc set 2s: MCPB.py step 2s' in note2s and 'the same file as in set 2e' in note2s
assert 'did NOT ACCEPT' in bn2f['status_note'] and bn2f['terms_at_N2']['given_file'] == 'md/qm/n2/n2_terms_new_charges_as_built.frcmod'
A(f'- **`bonded_2s/R1`: the other zinc set.** The plan (6 Oct 2026, 3:02 AM IST, point 3): "{q("For the Zn–N bonds and the N–Zn–N and Zn–N–C angles there are two sets", "the empirical constants (step 2e).")}" '
  f'The build record of the system (in its manifest): "{note2s}" (terms at N2: sha256 `{asb[:12]}…` in this run and in the nine). **What failed for it:** gate G-V in its validation run X0, in one part: an oxygen within 2.6 Å of the zinc in {fl_s["full_trajectory"]["count"]} of {fl_s["full_trajectory"]["n"]} all-atom frames '
  f'({pc(fl_s["full_trajectory"]["count"], fl_s["full_trajectory"]["n"], 0)} %; at most {pc(fl_s["limit_max_fraction"], 1, 0)} % allowed), a sulfonyl oxygen each time; {fl_s["solute_trajectory"]["count"]} of {fl_s["solute_trajectory"]["n"]} solute frames '
  f'({pc(fl_s["solute_trajectory"]["count"], fl_s["solute_trajectory"]["n"], 1)} %). Every other part of the gate is met. The plan (entry of 7 Oct 2026, 12:10 AM IST, point 4): "{q("set 2s missed the limit by one frame of 50", "two give 3 or more.")}" and '
  f'"{q("Set 2s is not to be described as five-coordinate or as an invalid model.")}" Mean distance of the nearer sulfonyl oxygen from the zinc in X0: {f(x0s["stats"]["d_Zn_Osulfonyl_min_A"]["mean"], 3)} Å with 2s, {f(x0e["stats"]["d_Zn_Osulfonyl_min_A"]["mean"], 3)} Å with 2e, '
  f'{f(min(REF["Zn_O1_A"], REF["Zn_O2_A"]), 3)} Å in 6RL9 (`results/runs/bonded_2s/X0/first_program_summary.json`, `results/runs/bonded_2e/X0/first_program_summary.json`). '
  '**What this run varies and what it does not:** the Zn–N bonds and the N–Zn–N and Zn–N–C angles. The terms at N2 and the charges are those of the nine runs, so this run tests neither (section 7, point 3).')
cb = N2G['candidates']['charges_as_built']
c5 = cb['conditions']['5 leave-one-out error on the fitting set']
c6 = cb['conditions']['6 force-field minimum of the small model against its quantum-chemical optimum']
flat = N2G['for_the_record_no_weight_in_the_verdict']['relaxed_geometries']['flattest_nitrogen_of_the_quantum_chemical_run_from_Rx']
topo = N2G['for_the_record_no_weight_in_the_verdict']['relaxed_geometries']['topologies']
assert topo['T0'] == 'terms and charges as built' and topo['T1'] == 'candidate, charges as built'
assert cb['failed'] == ['5 leave-one-out error on the fitting set', '6 force-field minimum of the small model against its quantum-chemical optimum']
A(f'- **`bonded_2e_n2fit/R1`: candidate terms at N2.** The build record of the system: "{bn2f["status_note"]}" The terms file (sha256 `{cb["frcmod_sha256"][:12]}…`) is, in the builder\'s words, "{bn2f["terms_at_N2"]["class_of_this_set"]}"; '
  f'beyond the terms as built it holds {len(bn2f["terms_at_N2"]["entries_beyond_the_terms_as_built"])} entries (torsions through N2 and one improper). **What failed for it:** gate G-R2 did not accept it (`md/N2.md`; `md/qm/n2/n2_gate_result.json`, verdict "{N2G["verdict"]}"): leave-one-out error {f(c5["value"])} kcal mol⁻¹ against the limit of {f(N2G["limits"]["loo_rms"], 1)}, '
  f'and the force-field minimum of the small model {f(c6["largest_abs_difference_deg"], 1)}° from the quantum-chemical optimum in Zn–N2–S against the limit of {f(N2G["limits"]["minimum_deg"], 0)}°; the test set was not read. '
  f'The plan uses it "only to bracket": "{q("at the relaxed near-planar nitrogen the terms as built are", "too low.")}" In the record (`for_the_record_no_weight_in_the_verdict.relaxed_geometries`), at the flattest nitrogen of the quantum-chemical relaxation (signed volume {flat["signed_volume"]:+.2f}): '
  f'terms as built {flat["force_field_minus_QM"]["T0"]:+.1f} kcal mol⁻¹ against quantum chemistry, this candidate {flat["force_field_minus_QM"]["T1"]:+.1f} (the candidate with the second set of charges, which was not run, {flat["force_field_minus_QM"]["T2"]:+.1f}). '
  '**The candidate differs from the terms as built in more than this one number.** In the record of gate G-R2 (`md/qm/n2/n2_gate_result.json`, `candidates.charges_as_built.for_the_record`, kcal mol⁻¹ relative to X; `md/N2.md`, sections 6.6, 6.7 and 9.4): '
  f'it puts the R-like placement {f(G_RX["new"])} above X, where the terms as built put it {f(G_RX["as_built"])} and quantum chemistry {f(G_RX["QM"])}; '
  f'it holds the hydrogen too near the plane (the hydrogen moved into the plane costs {sg(G_HY["new"], 2)} against {sg(G_HY["QM"], 2)} in quantum chemistry and {sg(G_HY["as_built"], 2)} as built); '
  f'its force-field minimum of the small model has Zn–N2–S {f(abs(c6["value"]["from_X"]["Zn_N2_S"]), 1)}° below the quantum-chemical optimum of {f(c6["optimum"]["ZnN2S"], 1)}°; '
  f'and, in the words of `md/N2.md`, "{cut(N2MD, "at the geometries the rehearsal visited, the candidate with the charges as built lowers the R-like frames", "the frame of the conversion by 6.59")}". '
  'So the MIXED run cannot be put down to the planar nitrogen alone: the candidate also changes the balance between the two placements and the geometry at N2.')
A('')
A('| | R1 | R2 | R3 | 2s/R1 | 2e_n2fit/R1 |')
A('|---|---|---|---|---|---|')
rows6 = [
    ('Zinc set; terms at N2', set_words),
    ('Water oxygens in the place of the ring when the ligand was released (driver); extension steps', lambda l: f'{H[l]["hydration_check"]["counts"][-1]["n_water_O"]}; {H[l]["hydration_check"]["n_extensions"]}'),
    ('First frame with the sign at N2 positive (from the release)', lambda l: ev(l, 'first_released_frame_n2_sign_not_of_start')),
    ('First frame that is not rotated-like; stage of these two frames', lambda l: f'{ev(l, "first_released_frame_not_in_state_of_start")}; {E[l]["first_released_frame_n2_sign_not_of_start"]["stage"]} and {E[l]["first_released_frame_not_in_state_of_start"]["stage"]}'),
    ('First crystal-like frame', lambda l: ev(l, 'first_released_frame_in_opposite_state')),
    ('Time of the plan\'s rule (start of the first crystal-like visit of at least 1 ns); stage', lambda l: f'{f(cps[l], 0)} ps; {conv[l]["stage"]}'),
    ('Verdict (last 20 ns)', lambda l: o[l]['verdict']),
    ('Changes of the sign at N2 from the release to the end of the production', lambda l: th(E[l]['n2_volume']['from_the_release']['n_changes_of_sign'])),
    ('Rotated-like frames from the release; of these with the sign at N2 positive (the configuration of the crystal pose); within 1.0 Å of the docked rotated pose',
     lambda l: f'{E[l]["state_and_n2_sign"]["from_the_release"]["R"]}; {E[l]["state_and_n2_sign"]["from_the_release"]["R_with_sign_positive"]}; {SC[l]["rotated_like_frames"]["from_the_release"]["n_within_1.0_A_of_the_docked_rotated_pose"]}'),
    ('Production: rotated-like frames with the sign at N2 positive; crystal-like frames with the sign negative',
     lambda l: f'{E[l]["state_and_n2_sign"]["production"]["R_with_sign_positive"]} of {E[l]["state_and_n2_sign"]["production"]["R"]}; {E[l]["state_and_n2_sign"]["production"]["C_with_sign_negative"]} of {th(E[l]["state_and_n2_sign"]["production"]["C"])}'),
    ('Production: frames with the ligand within 1.0 Å of the docked rotated pose; smallest RMSD to that pose (Å)',
     lambda l: f'{SC[l]["ligand_against_the_docked_rotated_pose"]["production_frames_within_1.0_A"]}; {f(SC[l]["ligand_against_the_docked_rotated_pose"]["production_smallest_A"])}'),
    ('Production: frames with the signed volume at N2 below 0.2 in size; with the sum of the three angles at N2 above 355°; Zn–N2–H, mean (°)',
     lambda l: f'{SC[l]["nitrogen"]["production"]["n_abs_V_below_0.2"]}; {th(SC[l]["nitrogen"]["production"]["n_sum_of_angles_above_355_deg"])}; {f(SC[l]["nitrogen"]["production"]["ang_Zn_N2_H_mean_deg"], 1)}'),
    ('Last 20 ns: crystal-like / rotated-like / intermediate (%)', lambda l: f'{pc(o[l]["count_C_window"], 2000, 2)} / {pc(o[l]["count_R_window"], 2000, 2)} / {pc(o[l]["count_I_window"], 2000, 2)}'),
    ('Whole production: crystal-like / rotated-like / intermediate frames', lambda l: f'{th(o[l]["count_C"])} / {o[l]["count_R"]} / {o[l]["count_I"]}'),
    ('Visits to the rotated-like state in the production (assigned state); longest (ns)', lambda l: f'{E[l]["production_visits_by_assigned_state"]["R"]["n_visits"]}; {f(E[l]["production_visits_by_assigned_state"]["R"]["longest_visit_frames"] * 0.01)}'),
    ('Frames with the sign at N2 negative, production', lambda l: f'{S1[l]["n2_configuration"]["count_negative"]} ({pc(S1[l]["n2_configuration"]["count_negative"], 10000, 2)} %)'),
    ('Hydrogen bond N2–H···OG1(Thr199), % of production frames', lambda l: pc(fr[l]['hb_N2H_T199OG1']['count'], 10000, 2)),
    ('Ligand RMSD to the crystal pose: mean ± SD; mean of the last half; mean of the crystal-like frames (Å)',
     lambda l: f'{f(st[l]["lig_rmsd_A"]["mean"])} ± {f(st[l]["lig_rmsd_A"]["sd"])}; {f(st[l]["lig_rmsd_A"]["mean_last_half"])}; {f(E[l]["production_by_state"]["C"]["lig_rmsd_mean_A"])}'),
    ('Frames whose smallest RMSD needs the ring turned by 180° (labelling `ring_flip`)', lambda l: th(S1[l]['lig_rmsd_labelling_counts']['ring_flip'])),
    ('τ_fold, mean ± SD (°)', lambda l: f'{f(st[l]["tau_fold_deg"]["mean"], 1)} ± {f(st[l]["tau_fold_deg"]["sd"], 1)}'),
    ('Zn–N2, mean (Å); S–N2–Zn, mean (°)', lambda l: f'{f(st[l]["d_Zn_N2_A"]["mean"], 3)}; {f(st[l]["ang_S_N2_Zn_deg"]["mean"], 1)}'),
    ('Sulfonyl oxygen within 2.6 Å of the zinc: % of solute frames; % of all-atom frames', lambda l: f'{pc(fr[l]["fifth_sulfonylO"]["count"], 10000, 2)}; {pc(fr[l]["fifth_any"]["count"], 1000, 1)}'),
    ('All-atom frames with a water oxygen in the place of the ring, production (of 1,000)', lambda l: str(E[l]['water_in_the_place_of_the_ring_production']['frames_with_one_or_more'])),
    ('Zinc out of the plane of the His119 ring, mean (Å)', lambda l: f(SC[l]['departures']['zinc_out_of_the_plane_of_the_His119_ring_A'])),
    ('Place of crystal water 553 occupied; a water donating a hydrogen bond to N2 (% of the all-atom frames of the production)',
     lambda l: f'{f(SF[l]["crystal_water_places"]["production"]["553"], 1)}; {f(SF[l]["water_donating_to_N2"]["share_production"], 1)}'),
    ('His64 deep in (Zn–NE2 at most 7.0 Å, χ1 between 0 and 120°), % of the production from 20 to 100 ns', lambda l: f(SC[l]['his64']['deep_in_share_20_to_100_ns']['7.0'], 1)),
]
for name, fn in rows6:
    A(f'| {name} | ' + ' | '.join(fn(l) for l in RR + SENS) + ' |')
A('')
A('The released equilibration of the two runs, as in section 4:')
A('')
A('```')
for l in SENS:
    A(f'{NAME[l]:<12} {strip_text(l)}')
A('```')
A('')
vs, vn = e2s['production_visits_by_assigned_state'], en2['production_visits_by_assigned_state']
ks2 = e2s['first_released_frame_n2_sign_not_of_start']
b2 = e2s['released_before_the_conversion']
assert b2['n2_negative'] == 2 and s2s[:2] == '--' and e2s['released_frames_states'][:2] == 'RR' and s2s.index('+') == 2
ssn2 = {l: E[l]['state_and_n2_sign'] for l in FR}
a2s = S1['bonded_2s/R1']['gate_GV']['GV3_zinc']['criteria']
a94, a96 = a2s['ang_N2_Zn_H94NE2_deg'], a2s['ang_N2_Zn_H96NE2_deg']
assert a94['mean_deg'] > a96['mean_deg'] and a94['reference_6RL9_deg'] < a96['reference_6RL9_deg']       # exchanged in size
oop9 = [SC[l]['departures']['zinc_out_of_the_plane_of_the_His119_ring_A'] for l in PROD]
A(f'**Zinc set 2s.** The verdict of the run is CONVERTED. The sign at N2 turned positive {f(ks2["time_from_release_ps"], 0)} ps after the release '
  f'({e2s["held_part"]["n2_positive"]} of the {e2s["held_part"]["n"]} held frames already had it positive; the hydrogen is not held) and the first frame that is not rotated-like came at {f(t_notr["bonded_2s/R1"], 0)} ps, both in stage 4; '
  f'the first crystal-like frame came at {ev("bonded_2s/R1", "first_released_frame_in_opposite_state")}. '
  f'Until {f(cps["bonded_2s/R1"], 0)} ps the ring moved between the two sides ({b2["C"]} crystal-like, {b2["R"]} rotated-like and {b2["I"]} intermediate frames before {f(cps["bonded_2s/R1"], 0)} ps), '
  f'with the sign at N2 positive in every frame from {f(ks2["time_from_release_ps"], 0)} ps on (the {b2["n2_negative"]} frames before it are rotated-like with the sign negative). '
  f'The time of the plan\'s rule is {f(cps["bonded_2s/R1"], 0)} ps, in stage {conv["bonded_2s/R1"]["stage"]}, with nothing held: later than in R1 to R3 by the rule ({rng([cps[l] for l in RR])} ps), '
  f'while the change of sign came earlier ({f(ks2["time_from_release_ps"], 0)} ps against {rng([t_sign[l] for l in RR])} ps). '
  f'In the production {o["bonded_2s/R1"]["count_R"]} frame is rotated-like (at {f(vs["R"]["longest_visit_starts_ns"])} ns), with the sign at N2 positive. '
  f'Of the {ssn2["bonded_2s/R1"]["from_the_release"]["R"]} rotated-like frames from the release {ssn2["bonded_2s/R1"]["from_the_release"]["R_with_sign_positive"]} have the sign at N2 of the crystal pose: '
  f'a ring turned towards the other side on a nitrogen that has already changed. {SC["bonded_2s/R1"]["rotated_like_frames"]["from_the_release"]["n_within_1.0_A_of_the_docked_rotated_pose"]} of the {ssn2["bonded_2s/R1"]["from_the_release"]["R"]} lie within 1.0 Å of the docked rotated pose. '
  f'The ligand lies further from the crystal pose than with set 2e (mean RMSD {f(st["bonded_2s/R1"]["lig_rmsd_A"]["mean"])} Å against {f(min(st[l]["lig_rmsd_A"]["mean"] for l in RR))} to {f(max(st[l]["lig_rmsd_A"]["mean"] for l in RR))} Å), as in the validation run of the set '
  f'({f(x0s["validation"]["lig_rmsd_mean_A"])} Å in X0 of 2s against {f(x0e["validation"]["lig_rmsd_mean_A"])} Å in X0 of 2e); a sulfonyl oxygen is within 2.6 Å of the zinc more often (table); '
  f'the zinc lies {f(SC["bonded_2s/R1"]["departures"]["zinc_out_of_the_plane_of_the_His119_ring_A"])} Å out of the plane of the His119 ring ({rng(oop9, 2)} Å in the nine runs); '
  f'and the mean angles N2–Zn–NE2(His94) and N2–Zn–NE2(His96) are {f(a94["mean_deg"], 1)}° and {f(a96["mean_deg"], 1)}°, exchanged in size against 6RL9 ({f(a94["reference_6RL9_deg"], 1)}° and {f(a96["reference_6RL9_deg"], 1)}°; `results/tables.md`, table A). '
  'Over 100 ns set 2s holds the site and the ligand less well than set 2e, as it did in its validation run.')
A('')
bl = en2['production_blocks_of_20_ns']
nn2 = SC['bonded_2e_n2fit/R1']
assert vs['R']['frames_assigned'] == o['bonded_2s/R1']['count_R'] == 1
assert all(S1[l]['lig_rmsd_labelling_counts']['ring_flip'] == 0 for l in ALL if l != SENS[1])       # "0 frames in every other run"
A(f'**Candidate terms at N2.** The verdict of the run is MIXED, so by the rule it has no conversion. The sign at N2 turned positive and the first crystal-like frame came {f(en2["first_released_frame_in_opposite_state"]["time_from_release_ps"], 0)} ps after the release, '
  f'and the first visit of at least 1 ns in the crystal-like state began at {f(cps["bonded_2e_n2fit/R1"], 0)} ps (stage {conv["bonded_2e_n2fit/R1"]["stage"]}). '
  f'In the production the ring is crystal-like in {pc(o[SENS[1]]["count_C"], 10000, 2)} % of the frames, rotated-like in {pc(o[SENS[1]]["count_R"], 10000, 2)} % and intermediate in {pc(o[SENS[1]]["count_I"], 10000, 2)} %. It visited the rotated-like state {vn["R"]["n_visits"]} times '
  f'({vn["R"]["frames_assigned"]} frames by assigned state), the longest visit lasting {f(vn["R"]["longest_visit_frames"] * 0.01)} ns (from {f(vn["R"]["longest_visit_starts_ns"])} ns), so no visit reaches the 1 ns of a long visit and the run has one transition. '
  f'The sign at N2 is negative in {pc(S1[SENS[1]]["n2_configuration"]["count_negative"], 10000, 2)} % of the frames ({S1[SENS[1]]["n2_configuration"]["n_stretches"]} stretches of one sign in the production, and '
  f'{th(en2["n2_volume"]["from_the_release"]["n_changes_of_sign"])} changes of sign from the release to the end; one change in each of R1 to R3), '
  f'the signed volume is below 0.2 in size in {nn2["nitrogen"]["production"]["n_abs_V_below_0.2"]} production frames ({min(SC[l]["nitrogen"]["production"]["n_abs_V_below_0.2"] for l in RR)} to {max(SC[l]["nitrogen"]["production"]["n_abs_V_below_0.2"] for l in RR)} in R1 to R3), and '
  f'the ring is turned by 180° about its axis against the crystal labels in {pc(S1[SENS[1]]["lig_rmsd_labelling_counts"]["ring_flip"], 10000, 2)} % of the frames (an equivalent labelling of the same molecule; 0 frames in every other run).')
A('')
A(f'**A rotated-like frame in the two sensitivity runs is mostly not the docked rotated pose.** The state is read from d1 and d2 only. In the nine runs of set 2e it goes with the sign at N2 except in single frames at a change: '
  f'rotated-like with the sign positive in {lst([str(ssn2[l]["from_the_release"]["R_with_sign_positive"]) for l in RR])} frames of R1, R2 and R3 from the release, crystal-like with the sign negative in {sum(E[l]["state_and_n2_sign"]["from_the_release"]["C_with_sign_negative"] for l in PROD)} frames of the nine. '
  f'In the two sensitivity runs it does not. With set 2s, {ssn2[SENS[0]]["from_the_release"]["R_with_sign_positive"]} of {ssn2[SENS[0]]["from_the_release"]["R"]} rotated-like frames from the release have the sign of the crystal pose. '
  f'With the candidate terms, {ssn2[SENS[1]]["production"]["R_with_sign_positive"]} of the {ssn2[SENS[1]]["production"]["R"]} rotated-like frames of the production have the sign of the crystal pose and {ssn2[SENS[1]]["production"]["C_with_sign_negative"]} crystal-like frames have that of start R '
  f'({ssn2[SENS[1]]["from_the_release"]["R_with_sign_positive"]} of {ssn2[SENS[1]]["from_the_release"]["R"]} and {ssn2[SENS[1]]["from_the_release"]["C_with_sign_negative"]} from the release); its rotated-like frames have d1 {f(nn2["rotated_like_frames"]["from_the_release"]["d1_mean_A"])} Å and d2 {f(nn2["rotated_like_frames"]["from_the_release"]["d2_mean_A"])} Å on average, just beyond the midpoints '
  f'({f(M1, 1)} and {f(M2, 1)} Å; start R as built {f(REFR["d1_N1_Phe131CZ_A"])} and {f(REFR["d2_N1_Asn67OD1_A"])} Å), and the ligand is within 1.0 Å of the docked rotated pose in {nn2["ligand_against_the_docked_rotated_pose"]["production_frames_within_1.0_A"]} production frames '
  f'({pc(nn2["ligand_against_the_docked_rotated_pose"]["production_frames_within_1.0_A"], 10000, 2)} %). A count of rotated-like frames in these two runs is therefore not a count of the docked rotated pose (`results/run_events.json`, `results/site_checks_solute.json`). '
  f'The candidate also changes the crystal-like pose itself: the mean Zn–N2–H angle is {f(nn2["nitrogen"]["production"]["ang_Zn_N2_H_mean_deg"], 1)}° against {rng([SC[l]["nitrogen"]["production"]["ang_Zn_N2_H_mean_deg"] for l in PROD], 1)}° in the nine runs, '
  f'the mean τ_fold {f(st[SENS[1]]["tau_fold_deg"]["mean"], 1)}° against {f(REF["tau_fold_deg"], 1)}° in 6RL9, and the mean ligand RMSD of its crystal-like frames {f(en2["production_by_state"]["C"]["lig_rmsd_mean_A"])} Å '
  f'against {rng([E[l]["production_by_state"]["C"]["lig_rmsd_mean_A"] for l in RR], 2)} Å in R1 to R3; a water donates a hydrogen bond to N2 in {f(SF[SENS[1]]["water_donating_to_N2"]["share_production"], 1)} % of its all-atom frames '
  f'({rng([SF[l]["water_donating_to_N2"]["share_production"] for l in PROD], 1)} % in the nine runs). The share of rotated-like frames is not steady over the run:')
A('')
A('| Block of the production (ns) | ' + ' | '.join(f'{f(b["from_ns"], 0)} to {f(b["to_ns"], 0)}' for b in bl) + ' |')
A('|---|' + '---|' * len(bl))
A('| crystal-like / rotated-like / intermediate frames (of 2,000) | ' + ' | '.join(f'{b["C"]} / {b["R"]} / {b["I"]}' for b in bl) + ' |')
A('| frames with the sign at N2 negative | ' + ' | '.join(str(b['n2_negative']) for b in bl) + ' |')
A('| frames with the hydrogen bond to Thr199 | ' + ' | '.join(str(b['hb_N2H_T199OG1']) for b in bl) + ' |')
A('| mean ligand RMSD to the crystal pose (Å) | ' + ' | '.join(f(b['lig_rmsd_mean_A']) for b in bl) + ' |')
A('')
assert bl[-1]['C'] == o[SENS[1]]['count_C_window'] and bl[-1]['R'] == o[SENS[1]]['count_R_window']
A(f'By the rule of the plan the last block gives the verdict: {pc(bl[-1]["C"], 2000, 2)} % crystal-like is below the 90 % of a converted run and {pc(bl[-1]["R"], 2000, 2)} % rotated-like below the 90 % of a run that stayed, so the run is MIXED. '
  f'The plan on mixed results: "{q("Mixed results are reported run by run", "runs that disagree.")}" (the same blocks for all eleven runs: `results/tables.md`, table C).')
A('')

# =============================================================================================== 7
A('## 7. What the result rests on, and its limits')
A('')
A('Each point names its record. The direction of an error is given as the record gives it. Where an error has the sign of the outcome, this is said (point 3); beyond that, no point argues that an error helps or harms the outcome.')
A('')
va = GR['verdict_part_A']
ftr = cb['for_the_record']
pl = ftr['inversion_path_against_X']['cost_of_the_planar_nitrogen']
hy = ftr['hydrogen_moved_into_the_plane']
tn = ftr['turn_about_Zn_N2_against_X']['e_p030']
rx = ftr['Rx_against_X']
assert (round(va['rms'], 2), round(pl['QM'], 1), round(pl['as_built'], 1), round(hy['QM'], 1), round(hy['as_built'], 1), round(tn['QM'], 1), round(tn['as_built'], 1), round(rx['mismatch_before'], 1)) == (3.08, 6.6, 13.8, 0.5, 5.2, 4.2, 1.8, 3.4)
assert round(GR['with_d3']['part_B']['X_to_Rx']['mm_minus_qm'], 3) == round(rx['mismatch_before'], 3)
A(f'1. **The terms at N2 of the nine runs failed gate G-R** (`md/GR.md`; `md/qm/gr/gr_gate.json`). Part A: the force field misses the quantum-chemical energies of the zinc-site model by {f(va["rms"])} kcal mol⁻¹ root mean square over {va["n_compared"]} comparisons '
  f'against the criterion of {f(va["criterion"], 1)} ({f(va["rms_without_d3"])} without the dispersion term). Direction, in the plan\'s words (6 Oct 2026, 12:19 PM IST): "{q("The direction of the failure: along the constructed inversion path", "so that rotation is too free in the force field.")}" '
  f'In the record (kcal mol⁻¹, relative to X, quantum chemistry against force field): planar nitrogen on the constructed path {f(pl["QM"])} against {f(pl["as_built"])}; only the hydrogen moved into the plane {f(hy["QM"])} against {f(hy["as_built"])}; '
  f'the ligand turned by 30° about Zn–N2 {f(tn["QM"])} against {f(tn["as_built"])}; the R-like placement of the crystal ligand against X {f(rx["QM"])} against {f(rx["as_built"])}, that is {rx["mismatch_before"]:+.2f} in the force field '
  f'(part B of the gate, reported and not a criterion; plan, 1:34 PM IST: "{q("including the energy of the R-like configuration against the crystal-like one", "overstate the cost of a planar nitrogen.")}"). '
  f'What part A measures (plan, 1:09 PM IST, point 2): "{q(PARTA_A, "(0.5 against 5.2).")}"')
A(f'2. **No better model of the nitrogen was accepted (gate G-R2)** (`md/N2.md`; `md/qm/n2/n2_gate_result.json`). The plan (6 Oct 2026, 8:41 PM IST): "{q("No model of the allowed class meets the conditions on the fitting data alone", "(limit 5°).")}" '
  f'In the record: smallest leave-one-out error of the candidates {f(N2G["every_candidate_of_the_ladder"]["charges_as_built"]["smallest_leave_one_out_rms"])} (charges as built) and {f(N2G["every_candidate_of_the_ladder"]["second_step"]["smallest_leave_one_out_rms"])} (charges of two placements) against {f(N2G["limits"]["loo_rms"], 1)}. '
  f'The force field therefore stayed as built, by the user\'s decision (9:15 PM IST): "{q("production with the terms at N2 as built, both failed gates", "direction of the error.")}"')
gr_no = GR['without_d3']['part_B']['X_to_Rx']
assert round(gr_no['mm'], 2) == round(rx['as_built'], 2) and gr_no['mm_minus_qm'] < 0 < rx['mismatch_before']
w553_n2 = re.search(r'SAN304\.N2:([\d.]+)', WATCSV['553']['polar_partners']).group(1)
w553_t200 = re.search(r'THR200\.OG1:([\d.]+)', WATCSV['553']['polar_partners']).group(1)
assert 'the same in every system of this file' in PRODMD and 'Zinc terms and charges as in the set named in the folder.' in bn2f['status_note']
assert 'the same file as in set 2e' in note2s
don9 = [SF[l]['water_donating_to_N2']['share_production'] for l in PROD]
occ9 = [SF[l]['crystal_water_places']['production']['553'] for l in PROD]
hb9 = [100.0 * fr[l]['hb_N2H_T199OG1']['count'] / fr[l]['hb_N2H_T199OG1']['n'] for l in PROD]
star = lambda t: t.replace('6-31G*', '6-31G\\*')
A(f'3. **What these two records say about a conversion, and what they do not.** The conversion of a rotated start goes with a change of the configuration at N2 (section 4), so it runs through the terms of points 1 and 2. Their known errors, each with its direction: '
  f'(a) on the constructed path the force field as built puts the planar nitrogen about twice as high as quantum chemistry ({f(pl["as_built"])} against {f(pl["QM"])} kcal mol⁻¹), and at the relaxed near-planar nitrogen {flat["force_field_minus_QM"]["T0"]:+.1f} kcal mol⁻¹ too high, so it resists the change; '
  f'(b) it makes the turn about Zn–N2 too easy ({f(tn["as_built"])} against {f(tn["QM"])} at 30°), a point of the gate that is a rigid turn of the whole ligand and not the motion seen at the change (`md/GR.md`: "{cut(GRMD, "In the four turns about Zn–N2 the fitted terms do not change at all (0.00)", "electrostatics with the fixed atoms.")}"); '
  f'(c) it puts the unrelaxed R-like placement {rx["mismatch_before"]:+.2f} kcal mol⁻¹ too high against X. '
  '**Of these, (c) is the one with the sign of the outcome:** a placement that the force field draws too high is left more readily and returned to less readily than quantum chemistry would have it. Three things qualify the number. '
  'Its size in the protein is not known: it comes from rigid, unrelaxed points of the bare zinc-site model. '
  f'Its sign rests on the dispersion term of the reference: without D3(BJ) the same comparison gives {gr_no["mm_minus_qm"]:+.2f} (`md/qm/gr/gr_gate.json`, `without_d3`); `md/GR.md`: "{cut(GRMD, "Quantum chemistry puts it 10.84 above X with D3 and 19.93 without", "the sign of the mismatch depends on it.")}" '
  f'and "{star(cut(GRMD, "Two errors of the reference are not quantified", "(section 9)."))}" '
  f'And most of it comes from the charges, not from the terms at N2. The charges were fitted to the potential of pose X alone (plan, 6 Oct 2026, 3:02 AM IST, point 5: "{q("The potential is that of pose X.")}"), and `md/N2.md` records: '
  f'"{cut(N2MD, "The charges as built describe the potential of the placement they were fitted in", "was worried about, measured.")}" and "{cut(N2MD, "the energy of Rx against X goes from 14.22 to 11.70", "from +3.38 to +0.85.")}" '
  f'and "{cut(N2MD, "So a good part of what gate G-R found as", "was the charges.")}" '
  'All eleven long runs have these charges (`md/PRODUCTION.md`, section 2.6: "the same in every system of this file"; build record of the candidate system: "Zinc terms and charges as in the set named in the folder."). '
  'So no run of the campaign varies the charges, and the run with zinc set 2s shares both the terms at N2 and the charges with the nine. '
  'These are errors of the bare zinc-site model in vacuum; no record shows what they do, taken together, to a run in the protein, and none is claimed here. In particular the too high cost of a planar nitrogen does not allow the statement that a better model would convert as readily or more readily: '
  f'`md/GR.md` records "{cut(GRMD, "Nothing is relaxed in part A.", "the barrier of the model, is not known.")}", and the plan allows no statement on a barrier. '
  f'**What the gates could not see.** The model of gates G-R and G-R2 has, in the words of `md/GR.md`, "{cut(GRMD, "no Thr199, no Thr200, no Leu198, no water")}". In the runs the nitrogen of the crystal-like state has two partners that this model lacks. '
  f'One is OG1 of Thr199, which accepts N2–H in {rng(hb9, 2)} % of the production frames of the nine runs. The other is a water at the place of crystal water 553, which in 6RL9 lies {w553_n2} Å from N2 and {w553_t200} Å from OG1 of Thr200 (`md/system/audit/waters.csv`): '
  f'in the nine runs that place is occupied in {rng(occ9, 1)} % of the all-atom frames of the production, and a water donates a hydrogen bond to N2 in {rng(don9, 1)} % (`results/site_checks_full.json`; the review, section 3, has the same occupancies). '
  'Water 553 is one of the four waters removed to build start R (point 5), and the rotated ligand lies on its place. The effect of both partners is in the simulation through fixed charges and TIP3P water and was tested by no gate; '
  'the gates therefore say nothing about this part of the cost of leaving the crystal configuration. '
  f'Two records bear on the whole without settling it. (a) Of the rehearsal with the same terms the plan wrote (6 Oct 2026, 11:40 AM IST, point 1): "{q("The rotated start converted to the crystal-like pose 70 ps after the ligand was released", "these terms decide the outcome.")}" '
  f'(b) Of the relaxed zinc-site model (8:41 PM IST): "{q("In the zinc-site model in vacuum the R-like placement therefore has no minimum", "without Thr199, the pocket wall or water.")}" '
  f'The sensitivity run with the candidate terms (at the relaxed near-planar nitrogen {flat["force_field_minus_QM"]["T1"]:+.1f} kcal mol⁻¹, on the other side of quantum chemistry) is the bracket the plan asked for: there the nitrogen changes its configuration back and forth and the run is MIXED (section 6). '
  f'The verdict for start R is thus {o[RR[0]]["verdict"]} with the terms as built and {o[SENS[1]]["verdict"]} with the candidate: it is not the same on the two sides of the bracket.')
assert [o[l]['verdict'] for l in FR] == ['CONVERTED'] * 4 + ['MIXED']
rmin9 = min(SC[l]['ligand_against_the_docked_rotated_pose']['production_smallest_A'] for l in PROD)
A(f'4. **The three conversions happened in the equilibration, with the protein backbone held.** By the plan\'s rule the three runs of the nine that converted did so {rng([cps[l] for l in RR])} ps after the release '
  f'({", ".join(NAME[l] + " " + f(cps[l], 0) + " ps in stage " + str(conv[l]["stage"]) for l in RR)}). The rule gives the start of the first crystal-like visit of at least 1 ns. '
  f'The configuration at N2 changed earlier, at {lst([t_sign[l] for l in RR])} ps, and the first frame outside the rotated-like state came at {lst([t_notr[l] for l in RR])} ps: in stage 4 in all three runs. '
  f'In stage 4 the protein backbone is held with {K4} kJ mol⁻¹ nm⁻² (in stage 5 with {K5}; in stage 6 nothing is held) to its positions in `system.inpcrd` (`md/protocol.json`, `position_restraints.reference_positions`: "{REFPOS}"). '
  'For the protein these are the coordinates of the crystal structure (the driver\'s check in each manifest: Cα atoms of `system.inpcrd` on `reference.pdb` '
  f'{f(H[RR[0]]["hydration_check"]["reference_guard"]["ca_fit_rmsd_of_system_inpcrd_on_reference_pdb_A"], 1)} Å), that is, of the protein as it is with the ligand in the crystal pose. '
  '**So no run of the nine holds the rotated configuration at N2 with the backbone restraint lowered or removed**, and whether the rotated pose would last longer in a protein that is free to adapt was not tested. '
  f'The two sensitivity runs, which enter no count: `2s/R1` (CONVERTED) has the time of the rule at {f(cps[SENS[0]], 0)} ps in stage {conv[SENS[0]]["stage"]}, the change of sign at {f(t_sign[SENS[0]], 0)} ps and the first frame outside the rotated-like state at {f(t_notr[SENS[0]], 0)} ps; '
  f'`2e_n2fit/R1` (MIXED, so without a conversion) began its first crystal-like visit of at least 1 ns at {f(cps[SENS[1]], 0)} ps, had its first change of sign at {f(t_sign[SENS[1]], 0)} ps, and went back to the rotated-like state {en2["production_visits_by_assigned_state"]["R"]["n_visits"]} times in the production. '
  'The plan foresaw a conversion before time zero after the rehearsal and fixed how it is counted (6 Oct 2026, 11:40 AM IST, point 2), so rule 3 stands as worded. '
  f'It follows that the 100 ns of production of R1 to R3 sample the crystal-like state only ({sum(o[l]["count_R"] for l in RR)} rotated-like frames of {th(sum(run[l]["n_frames_solute"] for l in RR))}; in the production of the nine runs the ligand never comes nearer than {f(rmin9)} Å, '
  'as RMSD, to the docked rotated pose; `results/site_checks_solute.json`): for the ligand, the production is nine runs of the crystal-like pose. They do not show how long a rotated pose can last, and no rate is derived.')
r2w = [dict(x, ps=round(x['time_ns'] * 1000.0 + 1300.0, 3)) for x in E[RR[1]]['water_in_the_place_of_the_ring_equilibration'] if x['stage'] >= 4]
r2w = [x for x in r2w if x['ps'] < t_sign[RR[1]]]
assert r2w[-1]['n_water_ring_place'] == 0 and r2w[-1] == last_full_before(RR[1]) and all(x['n_water_ring_place'] > 0 for x in r2w[:-1])
assert all(H[l]['hydration_check']['counts'][-1]['n_water_O'] == 1 for l in FR[2:]) and all(H[l]['hydration_check']['counts'][-1]['n_water_O'] == 2 for l in FR[:2])
A(f'5. **Start R was built without four ordered waters, and with the place of the crystal ring empty.** The plan (7 Oct 2026, 12:10 AM IST, point 2): "{q("Start R holds 211: waters 423, 516, 553 and 568", "(the hydration check).")}" '
  f'(crystal water molecules in the built systems: {H[XR[0]]["build"]["composition"]["crystal_water_molecules"]} in X, {H[NR[0]]["build"]["composition"]["crystal_water_molecules"]} in N, {H[RR[0]]["build"]["composition"]["crystal_water_molecules"]} in R; build record in each manifest). '
  f'Direction named by the plan (6 Oct 2026, 1:03 AM IST, point 3): "{q("An empty pocket beside the ligand could draw the ring back", "look easier than it is.")}" What was found: when the ligand was released, the driver counted water oxygens within '
  f'{f(H[RR[0]]["hydration_check"]["settings"]["water_cutoff_A"], 1)} Å of a ring-carbon position of the crystal pose in every run from start R: '
  f'{", ".join(str(H[l]["hydration_check"]["counts"][-1]["n_water_O"]) + " in " + NAME[l] for l in RR + SENS)} '
  f'(R1 after one extension step of the held part; nearest water oxygen {lst([H[l]["hydration_check"]["counts"][-1]["nearest_water_O_A"] for l in FR], 3)} Å, cut-off {f(H[RR[0]]["hydration_check"]["settings"]["water_cutoff_A"], 1)} Å), '
  'and both programs reproduce each count of the driver (`results/tables.md`, table F). By this count no run was released with the place of the ring empty, and the rule was kept. '
  f'In three of the five runs (R3 and the two sensitivity runs) the count of 1 rests on a single water oxygen {rng([H[l]["hydration_check"]["counts"][-1]["nearest_water_O_A"] for l in FR[2:]], 3)} Å from a ring-carbon position, close to the cut-off. '
  f'The count is also one instant. In R2 the place held {lst([str(x["n_water_ring_place"]) for x in r2w[:-1]])} water oxygens at {lst([x["ps"] for x in r2w[:-1]])} ps after the release and none at {f(r2w[-1]["ps"], 0)} ps (nearest {f(r2w[-1]["d_water_ring_place_min_A"], 3)} Å), '
  f'the last all-atom frame before its change of sign at {f(t_sign[RR[1]], 0)} ps: the situation that the plan named applies to that moment of R2. In R1 and R3 the last all-atom frame before the change held {last_full_before(RR[0])["n_water_ring_place"]} and {last_full_before(RR[2])["n_water_ring_place"]}. '
  'All-atom frames are 100 ps apart, and nothing finer is claimed. '
  'After the conversions the count is 0 in most all-atom frames (section 4). The count sees water oxygens within a cut-off of six points every 100 ps and says nothing of how that water is held (`md/analysis/SPEC.md`, section 15); '
  'it is not a statement that the place of the ring was fully solvated. '
  f'No run from X or N, which start with the four ordered waters in the place of the rotated ring, has a rotated-like frame ({sum(o[l]["count_R"] for l in XR + NR)} of {th(sum(run[l]["n_frames_solute"] for l in XR + NR))}).')
A(f'6. **Zinc set 2s failed its validation narrowly**; its run is a sensitivity run of a set that did not pass (section 6). The plan (7 Oct 2026, 12:10 AM IST, point 4): "{q("What is solid is the mean distance of the nearer sulfonyl oxygen from the zinc", "2.986 Å in 6RL9")}".')
A('7. **In set 2e the zinc geometry is put in** (section 5): the Zn–N distances and N–Zn–N angles of the runs are not an independent reproduction of the crystal structure.')
assert 'Runs | 3 starts × 3 independent runs (new velocities) × 100 ns' in PLAN and 'protein ff14SB, water TIP3P' in PLAN
A(f'8. **Design.** Three runs per start, 100 ns each, one force field (protein ff14SB, water TIP3P, ligand GAFF2 with RESP charges and repaired terms at sulfur, bonded zinc site; plan, section 3). The specification of the analysis: "{cut(SPEC, "Fractions of time in a state are not populations", "(plan, rule 5).")}" '
  f'And the plan, section 2: "{q("What the simulation cannot show, and will not be said to show:", "keeps by construction.")}"')
tm = [H[l]['production_state_log']['temperature_K']['mean'] for l in ELEVEN]
sm_lo = min(H[l]['production_state_log']['temperature_K']['segment_means_min'] for l in ELEVEN)
sm_hi = max(H[l]['production_state_log']['temperature_K']['segment_means_max'] for l in ELEVEN)
# added on 11 Oct 2026: the range of the plan at the 2 fs step against all 13 runs (appended to the item)
t2 = {l: [s_['log_mean_temperature_K'] for s_ in H[l]['equilibration_stages'] if s_['timestep_fs'] == 2.0 and s_['ensemble'] == 'NPT'] for l in ALL}
t2all = [v_ for l in ALL for v_ in t2[l]]
assert len(ALL) == 13 and all(len(t2[l]) >= 2 for l in ALL)
assert [f(v_, 1) for v_ in sorted(t2['bonded_2e/X0'])] == ['300.2', '300.6'], t2['bonded_2e/X0']
assert (f(min(t2all), 1), f(max(t2all), 1)) == ('299.8', '300.9')
NOTE_T2 = (f' **Added on 11 Oct 2026.** "300.2 to 300.6 K in the stages at 2 fs" in the sentence of the plan equals the two stage means of the validation run X0 of set 2e ({f(min(t2["bonded_2e/X0"]))} and {f(max(t2["bonded_2e/X0"]))} K). '
           f'Over all {len(ALL)} runs of the campaign the means of the {len(t2all)} equilibration stages at the 2 fs step with the barostat on are {f(min(t2all), 1)} to {f(max(t2all), 1)} K '
           f'({f(min(t2all))} to {f(max(t2all))} K; `results/run_health.json`, `equilibration_stages[].log_mean_temperature_K`). '
           'The plan is right for the run it was written from; for the campaign the range over the 13 runs is the right one, and it is the one in the facts of the paper (`md.grp.temp_2fs.min`, `md.grp.temp_2fs.max`). '
           'The numbers at the 4 fs step are not touched by this.')
A(f'9. **Kinetic temperature.** The plan (7 Oct 2026, 12:10 AM IST, point 6): "{q("The kinetic temperature of the log averages 301.0 K", "was not established by a test.")}" Over the whole campaign: mean of the state log {f(min(tm))} to {f(max(tm))} K in the eleven runs '
  f'(segment means {f(sm_lo)} to {f(sm_hi)} K), thermostat at {f(H[PROD[0]]["production_state_log"]["thermostat_K"][0], 0)} K, time step {f(H[PROD[0]]["production_state_log"]["timestep_fs"][0], 0)} fs (section 8).'
  + NOTE_T2)
A(f'10. **The docked starts.** Starts N and R are docked poses of the neutral molecule rebuilt as the anion; as built, Zn–N2 is {f(REFR["Zn_N2_A"], 3)} Å in R and {f(REFN["Zn_N2_A"], 3)} Å in N against {f(REF["Zn_N2_A"], 3)} Å in the crystal, and the hydrogen on N2 was placed by construction '
  '(plan, 6 Oct 2026, 3:02 AM IST, point 7; build record in each manifest). The plan on the minimisation (11:40 AM IST, point 4): the displacement limit of the dropped restrained model '
  f'"{q("is exceeded by the bonded model for the docked starts", "it is reported and not used as a gate.")}"')
A('11. **The coordinate.** The two distances separate two states only; a third orientation would show as intermediate frames. Water is seen every 100 ps; the hydrogen bonds are geometric criteria (`md/analysis/SPEC.md`, section 15). '
  'A rotated-like frame is a frame with d1 above and d2 below their midpoints; it is not by itself the docked rotated pose and not by itself the configuration of start R at N2. In the nine runs the three go together except in single frames at a change; '
  'in the two sensitivity runs they do not (section 6).')
nit = {l: SC[l]['nitrogen'] for l in ELEVEN}
n355 = [nit[l]['production']['n_sum_of_angles_above_355_deg'] for l in PROD]
nv02 = [nit[l]['production']['n_abs_V_below_0.2'] for l in PROD]
assert all(E[l]['n2_volume']['from_the_release']['n_changes_of_sign'] == 0 for l in XR + NR) and all(E[l]['n2_volume']['production']['n_changes_of_sign'] == 0 for l in PROD)
A(f'12. **The stiffness of the nitrogen: "stayed" and "no return" belong to the terms as built.** In the nine runs the signed volume at N2 is below 0.2 in size in {min(nv02)} to {max(nv02)} of {th(run[PROD[0]]["n_frames_solute"])} production frames per run, '
  f'and the sum of the three angles at N2 exceeds 355° in at most {max(n355)} ({pc(max(n355), run[PROD[0]]["n_frames_solute"], 2)} %). The six runs from X and N have no change of sign from the release to the end '
  f'({th(sum(E[l]["n2_volume"]["from_the_release"]["n"] for l in XR + NR))} frames), and R1 to R3 none after their one change. The gates found the terms as built too stiff in just this coordinate: '
  f'the hydrogen moved into the plane costs {f(hy["as_built"])} against {f(hy["QM"])} kcal mol⁻¹, the planar nitrogen on the constructed path {f(pl["as_built"])} against {f(pl["QM"])}, and at the relaxed near-planar nitrogen the terms are {abs(flat["force_field_minus_QM"]["T0"]):.1f} kcal mol⁻¹ too high (points 1 and 3). '
  '"Stayed" for the starts X and N and "no return" for the start R are therefore what these terms give in these nine runs. They are not statements about the molecule, and rule 1 cannot detect a force field that holds the crystal configuration too firmly: it tests whether the crystal pose is kept. '
  f'Two records show the other side. With the candidate terms the signed volume is below 0.2 in size in {nit[SENS[1]]["production"]["n_abs_V_below_0.2"]} production frames, the sum of the angles exceeds 355° in {th(nit[SENS[1]]["production"]["n_sum_of_angles_above_355_deg"])}, '
  f'the sign changes {th(E[SENS[1]]["n2_volume"]["from_the_release"]["n_changes_of_sign"])} times, and the one run from start R is MIXED by the same rule. And in the rehearsal, with the terms and the topology of the nine runs, the run from start N had {len(TN0_R)} rotated-like frames (section 3). '
  f'The crystal-like orientation holds the larger share on both sides of the bracket ({pc(o[SENS[1]]["count_C"], 10000, 2)} % of the production frames with the candidate); that the rotated-like one is absent holds only with the terms as built. '
  'The share with the candidate is that of one run of a model that failed its gate and that drifts over the run (section 6); it is not a population.')
A(f'13. **His64 and the water of the site** (section 4). Not examined by the two programs. His64 reaches an inner position that 6RL9 does not have and that empties the place of crystal water 516; it does so for {rng([dp(l) for l in RR], 1)} % of the production from 20 to 100 ns in R1 to R3, '
  f'for {rng([dp(l) for l in xn], 1)} % in the six runs from X and N, and for {f(dp(SENS[0]), 1)} and {f(dp(SENS[1]), 1)} % in the two sensitivity runs from start R. One hundred ns do not sample this coordinate, chance and a memory of the start cannot be told apart, '
  'and the limit that defines the inner position was chosen after looking at the data. The water structure of the site is not called reproduced or converged, and His64 was run in one tautomer and one protonation state only (HIE, built in conformer A).')
rv_q = lambda a, b=None: cut(REVIEW, a, b)
A('14. **The physical review: what is taken from it, and what is not.** The review (`results/independent/physical_review/REVIEW.md`) is in git as its author wrote it; it used its own code and did not read the output of the two programs. '
  'Its observations were checked from the files and are entered above: the stage and the place of the change at N2 (point 4), the water at N2 (point 3), R2 without water before its change (point 5), the stiffness of the nitrogen (point 12), '
  'His64 and water 516 (point 13), the time of the rule against the change of sign in R1, set 2s, and what a rotated-like frame is in the sensitivity runs (sections 4 and 6). Six of its statements go beyond a record or against the plan and are not taken over. '
  f'(1) "{rv_q("A lower cost would make it faster, not rarer.")}" and "{rv_q("The observation", "does not rest on error (a).")}": statements on a rate and on a barrier, which the plan excludes and which the gate record does not carry (point 3). '
  f'(2) "{rv_q("Error (c) helps the event a little.")}": the 30° point of the gate is a rigid turn of the whole ligand, not the swing of the hydrogen seen at the change (point 3). '
  f'(3) "{rv_q("did not persist in any of the five runs started from it")}", "{rv_q("within 0.4 ns")}" and "{rv_q("the two sensitivity runs converted even earlier")}": the sensitivity runs enter no count, one of them is MIXED, '
  f'0.4 ns is the time of the change of sign, and the times of the plan\'s rule are {f(cps[RR[0]], 0)} ps for R1 and {f(cps[SENS[0]], 0)} ps for 2s/R1. '
  f'(4) "{rv_q("Errors (a) and (b) are why set 2e shows no fluctuation at all.")}": true of the nine runs; the rehearsal run with the same topology did cross (section 3). '
  f'(5) "{rv_q("the runs from R spend 87 to 97 % of their last 80 ns there")}": this holds for R1 to R3 and not for the two sensitivity runs from start R. '
  f'(6) The {rx["mismatch_before"]:.1f} kcal mol⁻¹ by which the rotated placement is too high is given there without its dependence on the dispersion term (point 3).')
A('')

# =============================================================================================== 8
A('## 8. Health of the runs')
A('')
A('| Run | Driver: started; ended (IST) | Starts of the driver; exit code | Frames read / expected: production solute, all-atom; equilibration solute | DCD files with the driver\'s checksum | '
  'Temperature of the state log, mean ± SD (K); range of the segment means | Density, mean ± SD (g cm⁻³) | Cα RMSD to 6RL9: mean; last half; largest (Å) | Largest Cα RMSF (Å; residue) | Warnings, first / second program |')
A('|---|---|---|---|---|---|---|---|---|---|')


def row8(l):
    d = H[l]['driver']
    inv = d['invocations'][0]
    p = H[l]['production_state_log']
    eq = S1[l]['equilibration']
    eq_read = sum(s['solute']['n_frames'] for s in eq['stages'])
    eq_exp = sum(s['solute']['n_frames_expected'] for s in eq['stages'])
    exp_s = round(run[l]['planned_ns'] * 1000.0 / run[l]['solute_interval_ps'])
    exp_f = round(run[l]['planned_ns'] * 1000.0 / run[l]['full_interval_ps'])
    assert eq_read == eq_exp == eq['n_frames_solute'] and run[l]['n_frames_solute'] == exp_s == d['production']['frames_solute'] and run[l]['n_frames_full'] == exp_f == d['production']['frames_full']
    assert p['rows'] == exp_s
    tc = H[l]['trajectory_checksums']
    return (f'| {NAME[l]} | {inv["started_ist"].replace(" IST", "")}; {inv["ended_ist"].replace(" IST", "")} | {d["n_invocations"]}; {inv["exit_code"]} | '
            f'{th(run[l]["n_frames_solute"])} / {th(exp_s)}, {th(run[l]["n_frames_full"])} / {th(exp_f)}; {eq_read} / {eq_exp} | {tc["n_equal_to_the_driver_record"]} of {tc["n_dcd_files_read"]} | '
            f'{f(p["temperature_K"]["mean"])} ± {f(p["temperature_K"]["sd"])}; {f(p["temperature_K"]["segment_means_min"])} to {f(p["temperature_K"]["segment_means_max"])} | '
            f'{f(p["density_g_per_cm3"]["mean"], 4)} ± {f(p["density_g_per_cm3"]["sd"], 4)} | {f(st[l]["ca_rmsd_A"]["mean"])}; {f(st[l]["ca_rmsd_A"]["mean_last_half"])}; {f(st[l]["ca_rmsd_A"]["max"])} | '
            f'{f(S1[l]["rmsf"]["max_A"])}; {S1[l]["rmsf"]["max_residue_number_deposited"]} | {len(S1[l]["warnings"])} / {len(S2[l]["warnings"])} |')


order8 = sorted(ALL, key=lambda l: datetime.datetime.strptime(H[l]['driver']['invocations'][0]['started_ist'].replace(' IST', ''), '%d %b %Y, %I:%M:%S %p'))
for l in order8:
    A(row8(l))
A('')
fv = HF['final_verification']
A(f'- **Complete.** Every run has all its frames: {th(run[PROD[0]]["n_frames_solute"])} solute and {th(run[PROD[0]]["n_frames_full"])} all-atom frames in each of the eleven long runs, {run[VAL[0]]["n_frames_solute"]} and {run[VAL[0]]["n_frames_full"]} in the two validation runs, and every equilibration step the frames expected for it '
  f'(R1 has {S1[RR[0]]["equilibration"]["n_frames_solute"]} equilibration frames, with its extension step). Both programs find the manifest and the files consistent, the DCD headers of the {len(run[PROD[0]]["segments_solute"])} segments of each long run contiguous (a validation run has one segment) and no frame repeated at a segment boundary, and give no warning for any of the {len(ALL)} runs.')
A(f'- **Nothing failed and nothing was repeated.** Each run folder records one start of the driver, ended with exit code 0 ("{H[PROD[0]]["driver"]["invocations"][0]["outcome"]}"); the driver\'s logs hold no line with a warning, an error or a failure; '
  f'the state logs hold no value that is not finite. All runs were made on {H[PROD[0]]["driver"]["gpu"]} ({H[PROD[0]]["driver"]["platform"]}, precision {H[PROD[0]]["driver"]["precision"]}, OpenMM {H[PROD[0]]["driver"]["openmm"]}) with the driver '
  f'`md/run_md.py` sha256 `{H[PROD[0]]["driver"]["script_sha256"][:12]}…` in every manifest.')
assert len({H[l]['driver']['script_sha256'] for l in ALL}) == 1 and len({H[l]['driver']['gpu'] for l in ALL}) == 1
A(f'- **Pulled and verified.** The guard\'s record of the end of the work (`/root/cqc-sfa-md/status/final_verification.json`, {fv["ist"]}; copied into `results/run_health.json`): verified {str(fv["verified"]).lower()}; the driver\'s checksums of {th(fv["driver_checksums"]["files_ok"])} files in {fv["driver_checksums"]["runs"]} run folders right, '
  f'{len(fv["driver_checksums"]["bad"])} bad, {len(fv["driver_checksums"]["missing"])} missing; {th(fv["pi_against_machine"]["same"])} of {th(fv["pi_against_machine"]["files"])} files the same on the Pi and on the machine. '
  'Here, in addition, the checksum of every DCD file as the first program read it equals the driver\'s record (table).')
dm = [H[l]['production_state_log']['density_g_per_cm3']['mean'] for l in ELEVEN]
A(f'- **Temperature and density.** Mean kinetic temperature {f(min(tm))} to {f(max(tm))} K and mean density {f(min(dm), 4)} to {f(max(dm), 4)} g cm⁻³ over the eleven long runs; first and second half of each run agree within '
  f'{f(max(abs(H[l]["production_state_log"]["temperature_K"]["mean_first_half"] - H[l]["production_state_log"]["temperature_K"]["mean_last_half"]) for l in ELEVEN))} K and '
  f'{f(max(abs(H[l]["production_state_log"]["density_g_per_cm3"]["mean_first_half"] - H[l]["production_state_log"]["density_g_per_cm3"]["mean_last_half"]) for l in ELEVEN), 4)} g cm⁻³. The temperature of the equilibration stages is in `results/tables.md`, table B. '
  'The 1 K above the thermostat is point 9 of section 7.')
A(f'- **Protein.** Mean Cα RMSD to the crystal structure {f(min(st[l]["ca_rmsd_A"]["mean"] for l in ELEVEN))} to {f(max(st[l]["ca_rmsd_A"]["mean"] for l in ELEVEN))} Å, largest single frame {f(max(st[l]["ca_rmsd_A"]["max"] for l in ELEVEN))} Å; '
  f'largest Cα RMSF {f(min(S1[l]["rmsf"]["max_A"] for l in ELEVEN))} to {f(max(S1[l]["rmsf"]["max_A"] for l in ELEVEN))} Å. Neighbouring Cα atoms are at most {f(max(run[l]["largest_neighbour_CA_distance_A"] for l in ALL))} Å apart in every frame (the protein is written whole).')
A('')

# =============================================================================================== 9
A('## 9. Agreement of the two programs')
A('')
t0 = [x for x in GATE_TIMES if x.startswith('start ')][0]
t1 = [x for x in GATE_TIMES if x.startswith('end ')][0]
dur = int(t1.split('epoch ')[1]) - int(t0.split('epoch ')[1])
tfmt = lambda x: datetime.datetime.fromtimestamp(int(x.split('epoch ')[1]), IST).strftime('%-d %b %Y, %-I:%M:%S %p IST')
A('The gate of stage S6 is one command (`md/analysis/RECHECK.md`, section 1; `md/PRODUCTION.md`, "After the runs"):')
A('')
A('```')
A('nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare-campaign /root/cqc-sfa-md/runs --run-primary')
A('```')
A('')
calls = REP['primary_calls']
assert all(c['returncode'] == 0 for c in calls)
called_for = [c['command'].split()[-1].replace('/root/cqc-sfa-md/runs/', '') for c in calls if '--campaign' not in c['command']]
assert sorted(called_for) == sorted(ELEVEN + ['preflight_a', 'preflight_b']) and sum(1 for c in calls if '--campaign' in c['command']) == 1
w1 = {c['copy']: c['source_written'] for c in SRC['copied']}
A(f'It was run once, from {tfmt(t0)} to {tfmt(t1)} ({dur // 60} min {dur % 60} s; one thread per process), and ended by itself with **exit status {GATE_STATUS}** '
  f'("{GATE_LOG.splitlines()[-3].replace("Result: ", "").rstrip(".")}"). Its full output is `results/gate/gate_S6.log`. '
  f'The second program (`recheck_md.py`, sha256 `{REP["script_sha256"][:12]}…`) analysed all {REP["n_runs"]} folders in this call. The first program (`analyse_md.py`, sha256 `{REP["primary_script_sha256"][:12]}…`) was called by the command {len(calls)} times, '
  f'all with return code 0: for the eleven long runs, for the two pre-flight test folders and once for the campaign table. For the two validation runs the output that the first program wrote when gate G-V was read was compared '
  f'(`summary.json` of X0 of 2e written {w1["md/results/runs/bonded_2e/X0/first_program_summary.json"]}, of X0 of 2s {w1["md/results/runs/bonded_2s/X0/first_program_summary.json"]}; not older than the trajectories, which is the rule of `RECHECK.md` for calling the first program again); it carries the same checksum of the program. '
  'Both checksums are those on record in `md/analysis/RECHECK.md`.')
RECHECK = (R / 'analysis/RECHECK.md').read_text(encoding='utf-8')
assert REP['script_sha256'][:16] in RECHECK and REP['primary_script_sha256'][:16] in RECHECK
A('')
A('| Run | Solute / all-atom frames | Values compared | Disagreements beyond tolerance | Threshold items | Keys in one program only | Warnings, first / second | Agree |')
A('|---|---|---|---|---|---|---|---|')
for r in REP['runs']:
    c = r['comparison']
    A(f'| {r["label"]} | {th(r["n_frames_solute"])} / {th(r["n_frames_full"])} | {th(c["n_values"])} | {c["n_disagree"]} | {c["n_threshold"]} | {c["n_keys_one_side"]} | '
      f'{c["summary"]["n_warnings_primary"]} / {c["summary"]["n_warnings_recheck"]} | {"yes" if c["agree"] else "**no**"} |')
A(f'| **all {REP["n_runs"]}** | | **{th(n_values_runs)}** | **{sum(r["comparison"]["n_disagree"] for r in REP["runs"])}** | {REP["n_threshold_items"]} | {sum(r["comparison"]["n_keys_one_side"] for r in REP["runs"])} | | {REP["n_runs_agree"]} of {REP["n_runs"]} |')
A(f'| the {len(ALL)} runs of this file | | {th(n_values_13)} | | | | | |')
A(f'| campaign rules (`rules_per_model`, folders left out, order of the runs) | | {n_values_camp} | {REP["campaign"]["n_disagree"]} | | | | yes |')
A('')
ld = REP['largest_differences']
pre = [r for r in REP['runs'] if r['label'].startswith('preflight')]
PRE_WARN = pre[0]['comparison']['summary']['warnings_primary'][0]
assert all(len(r['comparison']['summary']['warnings_primary']) == 1 == len(r['comparison']['summary']['warnings_recheck']) and
           r['comparison']['summary']['warnings_primary'][0] == PRE_WARN for r in pre)
figdiff = FIG['largest_difference_between_programs_in_the_columns_drawn_A']
figwords = ('are the same in the files of the two programs, digit for digit' if figdiff == 0 else
            'agree between the programs within %.4f Å' % figdiff)
nlab = sum(r['comparison']['summary']['n_labels'] for r in REP['runs'])
nlabdiff = sum(len(r['comparison']['summary']['labels_differ']) for r in REP['runs'])
assert all(v == 0 for k, v in ld.items() if k not in ('series dist', 'summary dist')), ld          # "0 in every angle, torsion, ..."
tsol = {k: sum(SC[l]['threshold'][k] for l in ALL) for k in ('n_d1', 'n_d2', 'n_hb_distance', 'n_hb_angle', 'n_sulfonyl_O', 'n_n2_sign', 'n_frames_any', 'n_frames_any_production')}
tful = {k: sum(SF[l]['threshold'][k] for l in ALL) for k in ('n_zinc_oxygen', 'n_ring_place', 'n_water_N1', 'n_ring_wide', 'n_frames_any', 'n_frames_any_production')}
marg = min(min(v for k, v in SC[l]['threshold']['nearest_margin'].items() if k.endswith('_A')) for l in ALL)
vsm = [l for l in ALL if SC[l]['threshold']['n_n2_sign']]
wd = max(max(v for k, v in SC[l]['beside_the_two_programs'][sub]['largest_difference'].items() if 'deg' not in k) for l in ALL for sub in ('analysis', 'analysis_recheck'))
wa = max(SC[l]['beside_the_two_programs'][sub]['largest_difference']['ang_N2_HN2_T199OG1_deg'] for l in ALL for sub in ('analysis', 'analysis_recheck'))
wf = max(max(SF[l]['beside_the_two_programs'][sub]['largest_difference'].values()) for l in ALL for sub in ('analysis', 'analysis_recheck'))
A(f'Largest differences over all values: {f(ld["series dist"], 4)} Å in a series (one unit of the last decimal written) and {ld["summary dist"]:.6f} Å in a `summary.json`; 0 in every angle, torsion, signed volume, time, count and fraction '
  f'(`largest_differences` of `recheck_report.json`). The two programs class no frame differently: that is what {REP["n_threshold_items"]} threshold items means (`RECHECK.md` uses the term for a frame that the two programs class differently at a threshold). '
  f'Frames that lie within the tolerance of a threshold (`md/analysis/SPEC.md`, section 14) do occur. By the check on the raw frames, {tsol["n_frames_any"]} solute frames of the {len(ALL)} runs do '
  f'({tsol["n_frames_any_production"]} of the production, {tsol["n_frames_any"] - tsol["n_frames_any_production"]} of the equilibration): {tsol["n_d1"]} within 0.001 Å of {f(M1, 1)} Å in d1, {tsol["n_d2"]} of {f(M2, 1)} Å in d2, '
  f'{tsol["n_hb_distance"]} of 3.5 Å in N2···OG1(Thr199), {tsol["n_hb_angle"]} within 0.01° of 135° in N2–H···OG1, {tsol["n_sulfonyl_O"]} within 0.001 Å of 2.6 Å in the nearer sulfonyl oxygen, '
  f'and {tsol["n_n2_sign"]} with the signed volume at N2 below 0.001 in size (all in {", ".join(NAME[l] for l in vsm)}); the nearest margin is {f(marg, 4)} Å. '
  f'In the all-atom frames, {tful["n_zinc_oxygen"]} have the nearest oxygen within 0.001 Å of 2.6 Å from the zinc and {tful["n_ring_place"]} have a water oxygen within 0.001 Å of the 2.5 Å of the place of the ring, all in the production; '
  f'the two wider water counts (3.5 Å about N1 and about the ring with N1) have such a water in {tful["n_water_N1"]} and {tful["n_ring_wide"]} frames ({tful["n_frames_any"]} all-atom frames in all). '
  'This reading and the two programs class every one of these frames alike, so no state, flag, count, verdict or gate part differs between the readings (`results/site_checks_solute.json`, `results/site_checks_full.json`; per run in `results/tables.md`, table G). '
  f'The two pre-flight folders are tests of the rented machine and no results; each has one warning in both programs (first program: "{PRE_WARN}"). '
  f'Texts that the specification does not fix differ in {nlabdiff} places and are not counted (the word for the kind of a column, and the description of the order of the atoms at N2; `RECHECK.md`, section 6, point 18).')
A('')
A(f'Beside the gate: the conversion times of section 2 are formed twice, from the `state` of each program, and give the same frame in all {len(ALL)} runs; the columns drawn in the figures {figwords}; the counts of `results/run_events.json` are the same from both. What the comparison cannot show: both programs follow the same specification, so a definition that is wrong in the text is wrong in both (`RECHECK.md`, section 7).')
A('')
# added on 11 Oct 2026: the count of this script against the count of the re-derivation of the paper's facts
assert n_third == 16 * len(ELEVEN) == 176
NOTE_THIRD = ('(Added on 11 Oct 2026: 176 is the count of the list of this script, 16 values for each of the 11 runs: the verdict, six state counts, the first passage, the number of transitions, the wording of the conversion, three counts and three means. '
              'The re-derivation of the facts of the paper, `paper/main/data/derive/d17_md.py`, compares another list with the same file, 20 values for each run, 220 in all: seven counts, eleven means, the verdict and the time of the rule; it finds no count, verdict or time different, and its largest difference of a mean is in `paper/main/data/derived/d17_md.json`. '
              'Both counts are right, each for its list, and neither list holds the other; the paper prints no count.) ')
A('**Two further readings, made with other code.** '
  f'(1) A blind third reading of the eleven long runs (`results/independent/third_reading/third_reading.json`; own DCD reader, no file of the two programs or of `md/results/` opened) gives the same verdicts, the same state counts of the whole run and of the last 20 ns, '
  'the same first passage, transitions and wording of the conversion, the same counts of the hydrogen bond, of the sulfonyl oxygen within 2.6 Å and of the frames with the sign at N2 negative, and the same mean ligand RMSD and Zn–N2 distance within 0.0001 Å: '
  f'{n_third} values compared by this script, {n_third_diff} different. '
  + NOTE_THIRD +
  f'(2) The check on the raw frames made for this record (`md/production/site_checks.py`; MDAnalysis reads the files, atoms by residue number and name) forms d1, d2, the hydrogen bond to Thr199, the nearer sulfonyl oxygen and the signed volume at N2 for every solute frame of the {len(ALL)} runs, '
  f'equilibration included ({th(sum(SC[l]["n_solute_frames"] for l in ALL))} frames), and the oxygens about the zinc and the three water counts for every all-atom frame ({th(sum(SF[l]["n_all_atom_frames"] for l in ALL))} frames). '
  f'Its largest difference from the series of either program is {f(wd, 5)} Å, {f(wa, 4)}° and {f(wf, 5)} Å in the all-atom series, about half a unit of the last digit written, and it classes no frame differently. '
  'Neither reading replaces the gate; both follow the same specification for the quantities named.')
A('')

# =============================================================================================== 10
A('## 10. Where every file is')
A('')
A('| What | Where |')
A('|---|---|')
A('| This file; its long tables | `md/RESULTS.md`; `md/results/tables.md` (both written by `md/production/assemble/make_results_md.py`) |')
A('| Contract and definitions | `md/PLAN.md` (rules of section 5; dated entries under "Changes to this plan"), `md/analysis/SPEC.md` (version 1.3), `md/analysis/RECHECK.md`, `md/PRODUCTION.md` |')
A('| The two analysis programs | `md/analysis/analyse_md.py` (first), `md/analysis/recheck_md.py` (second); not changed by this stage |')
A('| Gate of stage S6: wrapper, output, times, exit status | `md/results/gate/run_gate_S6.sh`, `gate_S6.log`, `gate_S6.times`, `gate_S6.exit_status` (copies of `/root/cqc-sfa-md/analysis_logs/`) |')
A('| Campaign tables and records of both programs; comparison report | `md/results/campaign/first_program/` (`campaign_table.md`, `campaign_table.csv`, `campaign.json`), `md/results/campaign/second_program/` (`campaign_table.csv`, `campaign.json`, `recheck_report.md`, `recheck_report.json`) |')
A(f'| Per run ({len(ALL)} runs): `summary.json` of both programs; the small series of the first program (equilibration solute and all-atom, production all-atom, RMSF) | `md/results/runs/<set>/<run>/first_program_summary.json`, `second_program_summary.json`, `first_program_timeseries_equil_solute.csv`, `first_program_timeseries_equil_full.csv`, `first_program_timeseries_full.csv`, `first_program_rmsf_ca.csv` |')
tlog = (RES / 'conversion_time_test.log').read_text().strip().splitlines()
assert re.fullmatch(r'(\d+) of \1 checks passed', tlog[-1]) and not any(ln.startswith('FAIL') for ln in tlog)
A(f'| Conversion times: program, test with made-up series ({tlog[-1]}), result | `md/production/conversion_time.py`, `md/production/test_conversion_time.py` (output: `md/results/conversion_time_test.log`), `md/results/conversion_times.json` |')
A('| Figures: program, files, sources | `md/production/make_result_figures.py`; `md/results/fig_orientation.png`, `.pdf`, `md/results/fig_ligand_rmsd.png`, `.pdf`; `md/results/fig_sources.json` |')
A('| Copies, health and event counts: program and records | `md/production/collect_results.py`; `md/results/SOURCES.json` (every copy with the checksum of its source, and every file left out), `md/results/run_health.json`, `md/results/run_events.json` |')
A('| Checks on the raw frames, added after the audits (His64, water at the crystal places and at N2, the nitrogen, the docked rotated pose, frames near a threshold): program and records | `md/production/site_checks.py`; `md/results/site_checks_solute.json`, `md/results/site_checks_full.json` |')
ind_p, ind_t = IND['parts']['physical_review'], IND['parts']['third_reading']
assert any(c['file'].endswith('physical_review/REVIEW.md') for c in ind_p['copied']) and any(c['file'].endswith('third_reading/third_reading.json') for c in ind_t['copied'])
A(f'| Independent readings, copied on {IND["copied_on"]} from the scratch folder of the session (checksums in `md/results/independent/SOURCES.json`) | physical review: `md/results/independent/physical_review/REVIEW.md` with its scripts, tables and records ({ind_p["n_copied"]} files); '
  f'blind third reading: `md/results/independent/third_reading/third_reading.json` with its scripts and checks ({ind_t["n_copied"]} files). Their per-frame arrays are not in git ({ind_p["n_left_out"] + ind_t["n_left_out"]} files, {f((ind_p["bytes_left_out"] + ind_t["bytes_left_out"]) / 1e6, 1)} MB, listed with checksums in the same file). '
  'The scripts of the two audits were in the same scratch folder and are not in git; what the audits found is in section 11, each point checked again from the files |')
A('| The rehearsal run cited in sections 3 and 7 (a test of the pipeline, no result) | `md/test/rehearsal_2e/runs/TN0/`, `md/test/REHEARSAL_2e.md` |')
big = [x for x in SRC['left_out'] if x['file'].endswith('timeseries_solute.csv')]
A(f'| **Not in git** ({SRC["n_left_out"]} files, {f(SRC["bytes_left_out"] / 1e6, 1)} MB, listed with checksums in `SOURCES.json`): the solute series of the production of both programs ({len(big)} files, up to {f(max(x["bytes"] for x in big) / 1e6, 1)} MB each; one line per 10 ps with every column of the specification), '
  'the second program\'s copies of the small series, the per-run figures of the first program. Also not in git: the final frame of each run (`final_solute.pdb` in each run folder) | `/root/cqc-sfa-md/runs/<set>/<run>/analysis/` and `analysis_recheck/`; campaign figures in `/root/cqc-sfa-md/runs/analysis_campaign/` |')
A(f'| Trajectories, manifests, state logs, driver logs ({f(HF["final_verification"]["driver_checksums"]["bytes_hashed"] / 1e9, 1)} GB in the files that the driver checksums) | `/root/cqc-sfa-md/runs/bonded_2e/` (X0, X1 to X3, N1 to N3, R1 to R3), `/root/cqc-sfa-md/runs/bonded_2s/` (X0, R1), `/root/cqc-sfa-md/runs/bonded_2e_n2fit/R1`; not changed by this stage beyond the `analysis/` and `analysis_recheck/` folders that the gate writes |')
A('| Records of the two gates on the terms at N2 | `md/GR.md`, `md/qm/gr/gr_gate.json`; `md/N2.md`, `md/qm/n2/n2_gate_result.json` |')
A('| Logs of this stage; all steps after the gate in their order | `/root/cqc-sfa-md/analysis_logs/`; `md/production/make_results.sh` (conversion times with their test, figures, copies and counts, the checks on the raw frames, this file) |')
A('')
A(f'**What the plan asks for and git does not yet hold.** Plan, section 6: "{q("inputs, parameters, final frames and analysis", "supplementary material.")}" Rule 6: "{q("All nine runs enter the supplementary material as time series", "whichever way they fall.")}" '
  'In git the production series of d1 and d2 exist as Figure 1 only (the equilibration series and the all-atom series are there in full), and the final frames are not there. '
  'Whether a compressed copy of the nine production series (time, d1, d2, state, signed volume at N2, ligand RMSD) and the final frames go into git is put to the user; nothing was added.')
# added on 11 Oct 2026: the series and the final frames are in git since stage S7 of the plan (appended to the paragraph
# above, so that every line of the file keeps its number)
ser = sorted((RES / 'series').glob('*.csv.gz'))
ffr = sorted((RES / 'final_frames').glob('*_final_solute.pdb.gz'))
assert [p_.name[:-len('.csv.gz')] for p_ in ser] == [p_.name[:-len('_final_solute.pdb.gz')] for p_ in ffr] == sorted(l.replace('/', '_') for l in ELEVEN)
for d_ in ('series', 'final_frames'):
    sums = dict(reversed(ln.split(None, 1)) for ln in (RES / d_ / 'SHA256SUMS').read_text().splitlines() if ln.strip())
    assert all(sums[p_.name.strip()] == sha(p_) for p_ in (ser if d_ == 'series' else ffr)), 'checksum of a file of results/%s' % d_
L[-1] += (f' **Added on 11 Oct 2026.** Both are in git since stage S7 of the plan: `md/results/series/` holds one gzip-compressed series for each of the {len(ser)} runs of 100 ns (the nine runs and the two sensitivity runs; '
  f'the released part of the equilibration and the production, one line per 10 ps), written by `md/production/make_series.py`, and `md/results/final_frames/` the last solute frame of each of the {len(ffr)} runs, '
  'written by `md/production/make_final_frames.py`; each folder has a `README.md`, `SOURCES.json` and `SHA256SUMS`, and this script checked the sums. The rows of the table above describe what git held on 8 Oct 2026.')
A('')
A(f'Checksums of the records this file was assembled from: `conversion_times.json` `{sha(RES / "conversion_times.json")[:16]}…`, `run_health.json` `{sha(RES / "run_health.json")[:16]}…`, `run_events.json` `{sha(RES / "run_events.json")[:16]}…`, '
  f'`site_checks_solute.json` `{sha(RES / "site_checks_solute.json")[:16]}…`, `site_checks_full.json` `{sha(RES / "site_checks_full.json")[:16]}…`, `independent/SOURCES.json` `{sha(RES / "independent/SOURCES.json")[:16]}…`, '
  f'`recheck_report.json` `{sha(RES / "campaign/second_program/recheck_report.json")[:16]}…`, first program\'s `campaign.json` `{sha(RES / "campaign/first_program/campaign.json")[:16]}…`, `md/PLAN.md` `{sha(R / "PLAN.md")[:16]}…`.')
A('')

# =============================================================================================== 11
A('## 11. Corrections after the two audits and the physical review of 8 Oct 2026')
A('')
A(f'The first assembly ({FIRST_ASSEMBLY}; sha256 `8a21f609f69f5ee7…`) was read by an audit of the numbers and an audit of the rules and the wording; the physical review ended after that assembly and none of its findings was in it. '
  'Each point below was checked from the files before anything was changed. No verdict, no rule, no threshold, no definition and no time of the plan\'s rule was changed; every correction is a last digit, a wording, or a limit that was missing. '
  'The corrections were made in the scripts and the file was assembled again; nothing was typed into it.')
A('')
A('| Point | What was found | What was done |')
A('|---|---|---|')
CORR = [
    ('N1', 'The audit found eleven numbers in eight cells of the table of section 8 (temperature and density of the state logs) one unit of the last digit off, because they were rounded when collected and again when printed.',
     '`collect_results.py` now writes these means and standard deviations as computed; they are rounded once, here. The new checks on the raw frames are written the same way.'),
    ('N2', f'Section 4, R1: the nearest water oxygen of the driver\'s check was printed as {"%.2f" % H[RR[0]]["hydration_check"]["counts"][-1]["nearest_water_O_A"]} Å for the {f(H[RR[0]]["hydration_check"]["counts"][-1]["nearest_water_O_A"], 3)} Å of the record.', 'The driver\'s three decimals are printed wherever this distance is given.'),
    ('N3', 'Section 9 said that no frame lies at a threshold. The report of the comparison shows only that the two programs class no frame differently.',
     'Worded as found; the frames within the tolerance of a threshold were counted on the raw frames and are given (section 9; `results/tables.md`, table G).'),
    ('N4, D9', f'Sections 1 and 7 counted five "conversions", among them the MIXED run, took the {f(max(cps.values()) / 1000.0)} ns of section 1 from a sensitivity run, and section 3 said that all eleven long runs meet every part of gate G-V, for which a production run has numbers without a verdict.',
     'The three conversions of the nine runs are named with their times; the two sensitivity runs are given apart, the MIXED one as a first crystal-like visit and not as a conversion; the sentence on gate G-V gives no verdict.'),
    ('N5', f'Section 6, set 2s: "the sign at N2 positive throughout" covered {b2["n2_negative"]} frames in which it is negative.', f'Worded as found: positive from {f(ks2["time_from_release_ps"], 0)} ps on; the {b2["n2_negative"]} frames before it are named.'),
    ('N6, D7', f'Section 7, point 5: the count of the driver was given without its margin (in three runs one water oxygen at {rng([H[l]["hydration_check"]["counts"][-1]["nearest_water_O_A"] for l in FR[2:]], 3)} Å, cut-off {f(H[RR[0]]["hydration_check"]["settings"]["water_cutoff_A"], 1)} Å), and it was not said that R2 had no water in the place of the ring in the last all-atom frame before its change.',
     'Both are in point 5 and in the paragraph on R2 of section 4.'),
    ('N7', 'Section 4 pointed for the ligand angles of R1 to R3 to places that do not hold them.', 'C4–S–N2 and τ_fold of R1 to R3 are in the table of section 4.'),
    ('N8', 'The header said that no number is typed; constants of the protocol and of the definitions are. Two statements were not guarded, and one count was read from a neighbouring quantity.',
     'The header says what is typed; the two statements are checked by the script; the count is read from `orientation.count_R`.'),
    ('N9', '`conversion_time.py` ended with an error instead of a message for a folder without an equilibration series (no run of the campaign).', 'It reports the production time alone with a plain message; a case was added to its test.'),
    ('D1', 'The time of the plan\'s rule was worded as the time at which the runs left the rotated state, and its stage as the stage of the conversion. In all five runs from start R the sign at N2 changed and the ligand left the rotated-like state in stage 4, with the backbone held.',
     'The times of the rule are kept as the plan fixed them; the time of the change of sign and of the first frame outside the rotated-like state, and their stage, are given beside them everywhere (header, sections 1, 2, 4, 6 and 7, point 4).'),
    ('D2', 'The direction of the errors of the terms at N2 was given without saying which error has the sign of the outcome, that its sign rests on the dispersion term, and that most of it comes from charges fitted to pose X alone, which no run varies.',
     'Section 7, point 3; section 6 says what the run with set 2s does not test.'),
    ('D3', 'The absence of any excursion of the nitrogen was reported without the limit that the terms as built hold it too stiffly, and the consequence of the bracket run was not stated.', 'Section 3 and section 7, point 12.'),
    ('D4', 'The rehearsal run from start N with the same topology, which had three rotated-like frames, was left out, and the plan\'s sentence on it was quoted without its first half.', 'Section 3 and section 7, points 3 and 12; the sentence is quoted in full.'),
    ('D5', '"Ends like the crystal-start runs" was wider than what was compared: His64 and the water of the site differ.', 'The answer is bounded to the ligand; His64, the crystal-water places and their limits are in section 4 and in section 7, point 13.'),
    ('D6', 'The water that donates a hydrogen bond to N2 (crystal water 553) is absent from the model of the gates; this was not reported.', 'Section 7, point 3; numbers in the tables of sections 4 and 6.'),
    ('D8', 'A rotated-like frame in the two sensitivity runs is mostly not the configuration of start R at N2 and not the docked pose; the record did not say so.', 'Section 1 (row of the candidate run), section 6 and section 7, point 11.'),
    ('D10', 'The candidate terms were described as a bracket in one number; the record of gate G-R2 shows that they also change the balance of the two placements and the geometry at N2.', 'Section 6, in the description of the candidate and in the paragraph on its run.'),
    ('D11', 'Gate G-Z was not mentioned; the sign at N2 was used without the limit that the specification gives it; three of the four departures of the plan were not measured again; the production series and final frames are not in git.',
     'Section 5 (gate G-Z; the three departures, measured on the raw frames), section 2 (the limit of the sign), section 10 (what git does not hold, put to the user).'),
    ('D12', 'Six statements of the physical review go beyond a record or against the plan.', 'They are named in section 7, point 14, and not taken over; the review itself is in git unchanged.'),
    ('D13', 'Figure 1 drew the triangle of a conversion on the MIXED run without a word, and did not show where stages 4 and 5 end.', 'The triangle of that run is open and explained in legend and caption; the ends of stages 4 and 5 are drawn.'),
    ('Pipeline', 'Found while these corrections were made: `make_results.sh` pipes every step to `tee` and did not stop when a step failed, so the record could have been assembled from the output of an older version of a program.',
     '`make_results.sh` now stops at a failed step (`pipefail`), and this script stops if a record was not written by the present version of its program (checksums of `conversion_time.py`, `make_result_figures.py` and `site_checks.py`).'),
    ('Review', 'Findings that limit a claim or must be reported (none invalidates a run): the conversions under backbone restraints; the event as an inversion at N2 through the failed terms; MIXED with the candidate terms; the stiffness of the nitrogen; water 553; His64 and water 516; the ligand after the conversion; the time of the rule in R1; R2 without water; set 2s; the state R in the sensitivity runs.',
     'Each was confirmed on the raw frames (`md/production/site_checks.py`) or from the series of the two programs and is entered: section 4, section 6 and section 7, points 3, 4, 5, 12 and 13.'),
]
for a_, b_, c_ in CORR:
    A(f'| {a_} | {b_} | {c_} |')
A('')
A('One number of the audit of the numbers is not repeated here as it was given: the audit counted five all-atom frames within the tolerance of a threshold; the check made for this record counts '
  f'{tful["n_zinc_oxygen"]} at the 2.6 Å about the zinc and {tful["n_ring_place"]} at the 2.5 Å of the place of the ring ({tful["n_zinc_oxygen"] + tful["n_ring_place"]} together), and more at the 3.5 Å of the two wider water counts (section 9). '
  'The classes of the frames do not depend on which count is taken: every reading classes them alike.')
A('')

# =============================================================================================== 12 (11 Oct 2026)
ASSEMBLY_8_OCT = ('8 Oct 2026, 6:45 AM IST', 'dfdb8452cbf69fa8')    # the file of commit a2bcbec; time from its header, first 16 digits of its sha256
A('## 12. Changes of 11 Oct 2026')
A('')
A(f'This file was assembled again on 11 Oct 2026 by the same script from the same records. Against the assembly of {ASSEMBLY_8_OCT[0]} (sha256 `{ASSEMBLY_8_OCT[1]}…`; compared line by line) it differs in the eight places listed here and in nothing else. '
  '(1) The time of the assembly in the header. '
  '(2) Section 6, the description of the candidate terms: the two places of `md/qm/n2/n2_gate_result.json` that the paragraph quotes are named by their keys. The numbers at the flattest nitrogen of the relaxation are under `for_the_record_no_weight_in_the_verdict.relaxed_geometries`; '
  'the old points of gate G-R (the R-like placement against X, the hydrogen moved into the plane) are under `candidates.charges_as_built.for_the_record`, which was printed as `for_the_record` alone, a name the record has no top-level key for. '
  '(3) Section 10: the sentences added to its last paragraph, on the series and the final frames that are in git now. '
  '(4) Section 10: the checksum of `md/PLAN.md`, because the plan received entries after that assembly (stage S6 closed and the rented machine deleted, both of 8 Oct 2026, and any dated entry that follows them in the plan). '
  '(5) This section. '
  '(6) Section 5, the paragraph on the departures that no gate tests: a dated note that the crystal value of the zinc out of the plane of the His119 ring is 0.15 Å, and that the 0.17 Å in the quoted sentence of the plan is the value of His94. '
  '(7) Section 7, point 9: a dated note that the range of the plan at the 2 fs step is that of one validation run, with the range over the 13 runs beside it. '
  '(8) Section 9: a dated note on the two counts of values of the third reading, 176 by this script and 220 by the re-derivation of the facts of the paper. '
  'No number of a run, verdict, rule or limit changed; the three notes set a reference value or a count right or beside another, and `md/results/tables.md` changed in its date line only.')
A('')

text = '\n'.join(L)
assert chr(0x2014) not in text, 'an em dash'
assert '{q(' not in text and '{f(' not in text and '{pc(' not in text, 'a passage was not filled in'
(R / 'RESULTS.md').write_text(text, encoding='utf-8')

# =============================================================================================== tables.md
T = []
B = T.append
B('# Results of the MD campaign: the long tables')
B('')
B(f'Written with `md/RESULTS.md` by `md/production/assemble/make_results_md.py` on {now}; every number is read from the files named.')
B('')
B('## A. The six N–Zn–N angles: mean of each run and its departure from 6RL9 (°)')
B('')
B('Source: `gate_GV.GV3_zinc.criteria` of each `summary.json` (first program; the second gives the same). Reference: the angle in 6RL9 as the specification prints it (section 11).')
B('')
B('| Run | ' + ' | '.join(ANGN[a] for a in ANG) + ' | S–N2–Zn |')
B('|---|' + '---|' * (len(ANG) + 1))
cr0 = S1[PROD[0]]['gate_GV']['GV3_zinc']['criteria']
B('| 6RL9 | ' + ' | '.join(f(cr0[a]['reference_6RL9_deg'], 1) for a in ANG) + f' | {f(REF["S_N2_Zn_deg"], 1)} |')
for l in ELEVEN + VAL:
    cr = S1[l]['gate_GV']['GV3_zinc']['criteria']
    B(f'| {NAME[l]} | ' + ' | '.join(f'{f(cr[a]["mean_deg"], 1)} ({sg(cr[a]["deviation_deg"])})' for a in ANG) + f' | {f(st[l]["ang_S_N2_Zn_deg"]["mean"], 1)} |')
B('')
B('## B. Equilibration: steps, restraints and mean temperature of the state log (K)')
B('')
B('Source: `results/run_health.json` (manifest of each run and its `log_<step>.csv`).')
B('')
B('| Step | Length (ps) | Time step (fs) | Ensemble | Atoms held; force constant (kJ mol⁻¹ nm⁻²) | ' + ' | '.join(NAME[l] for l in order8) + ' |')
B('|---|---|---|---|---|' + '---|' * len(order8))
tags = [s['tag'] for s in H[RR[0]]['equilibration_stages']]
for t in tags:
    ref_s = stg[t]
    cells = []
    for l in order8:
        m = {s['tag']: s for s in H[l]['equilibration_stages']}
        cells.append(f(m[t]['log_mean_temperature_K'], 1) if t in m else 'not run')
    B(f'| {t} ({ref_s["name"]}) | {f(ref_s["length_ps"], 0)} | {f(ref_s["timestep_fs"], 0)} | {ref_s["ensemble"]} | {ref_s["restraint_atoms"]}; {f(ref_s["restraint_k_kj_per_mol_nm2"], 0)} | ' + ' | '.join(cells) + ' |')
B('')
B('## C. The production in blocks of 20 ns')
B('')
B('Frames of each state (crystal-like / rotated-like / intermediate, of 2,000), frames with the sign at N2 negative, frames with the hydrogen bond N2–H···OG1(Thr199), mean ligand RMSD to the crystal pose (Å). Source: `results/run_events.json`, counted from the production series of both programs.')
B('')
B('| Run | 0 to 20 ns | 20 to 40 ns | 40 to 60 ns | 60 to 80 ns | 80 to 100 ns |')
B('|---|---|---|---|---|---|')
for l in ELEVEN:
    B(f'| {NAME[l]} | ' + ' | '.join(f'{b["C"]} / {b["R"]} / {b["I"]}; N2− {b["n2_negative"]}; HB {b["hb_N2H_T199OG1"]}; {f(b["lig_rmsd_mean_A"])}' for b in E[l]['production_blocks_of_20_ns']) + ' |')
B('')
B('## D. The frames around the change of sign at N2 in R1, R2 and R3')
B('')
B('From two frames before the first frame with the sign at N2 positive to two frames after the first crystal-like frame. Source: `results/run_events.json` (`frames_around_the_conversion`), from `timeseries_equil_solute.csv`.')
B('')
for l in RR:
    e = E[l]
    k0 = e['first_released_frame_n2_sign_not_of_start']['time_from_release_ps']
    k1 = e['first_released_frame_in_opposite_state']['time_from_release_ps']
    B(f'**{NAME[l]}**')
    B('')
    B('| Time from the release (ps) | State | Signed volume at N2 | N2–H···OG1 present | N2···OG1 (Å) | N2–H···OG1 (°) | Ligand RMSD (Å) | d1 (Å) | d2 (Å) |')
    B('|---|---|---|---|---|---|---|---|---|')
    for x in e['frames_around_the_conversion']:
        if k0 - 20 <= x['time_from_release_ps'] <= k1 + 20:
            B(f'| {f(x["time_from_release_ps"], 0)} | {x["state"]} | {x["n2_volume"]:+.3f} | {"yes" if x["hb_N2H_T199OG1"] else "no"} | {f(x["d_N2_T199OG1_A"])} | {f(x["ang_N2_HN2_T199OG1_deg"], 0)} | {f(x["lig_rmsd_A"])} | {f(x["d1_A"])} | {f(x["d2_A"])} |')
    B('')
B('## E. Water in the place of the crystal ring through the equilibration of the runs from start R')
B('')
B('All-atom frames every 100 ps; time from the start of the production (ns); count of water oxygens within 2.5 Å of the crystal positions of the six ring carbons (`n_water_ring_place`), distance of the nearest one (Å), and the state of the solute frame of the same time. Source: `results/run_events.json`.')
B('')
for l in RR + SENS:
    B(f'**{NAME[l]}**: ' + '; '.join(f'{f(x["time_ns"], 1)}: {x["n_water_ring_place"]} ({f(x["d_water_ring_place_min_A"])}, {x["state_of_the_solute_frame"]})' for x in E[l]['water_in_the_place_of_the_ring_equilibration']))
    B('')
B('## F. Hydration check of the driver and its reproduction by the two programs')
B('')
B('Source: `equilibration.stages[].driver_hydration_check` of each `summary.json`; `count_agrees` and `nearest_agrees` are true in both programs for every record.')
B('')
B('| Run | Step | Driver: count; nearest water oxygen (Å); applies; decision | Compared; count agrees; nearest agrees |')
B('|---|---|---|---|')
for l in order8:
    for s in S1[l]['equilibration']['stages']:
        d = s.get('driver_hydration_check')
        if d:
            d2 = [x for x in S2[l]['equilibration']['stages'] if x['tag'] == s['tag']][0]['driver_hydration_check']
            assert d['compared'] and d['count_agrees'] and d['nearest_agrees'] and d2['count_agrees'] and d2['nearest_agrees']
            B(f'| {NAME[l]} | {s["tag"]} | {d["n_water_O"]}; {f(d["nearest_water_O_A"], 3)}; {str(d["applies"]).lower()}; {d["decision"]} | {str(d["compared"]).lower()}; {str(d["count_agrees"]).lower()}; {str(d["nearest_agrees"]).lower()} |')
B('')
B('## G. Frames within the tolerance of a threshold, per run')
B('')
B('Own reading of the raw frames (`md/production/site_checks.py`; `results/site_checks_solute.json`, `results/site_checks_full.json`), equilibration and production together. Tolerances of `md/analysis/SPEC.md`, section 14: 0.001 Å, 0.01°, 0.001 for the signed volume. '
  'This reading, the first program and the second class every one of these frames alike.')
B('')
B('| Run | Solute frames read | d1 at 6.7 Å | d2 at 5.2 Å | N2···OG1 at 3.5 Å | N2–H···OG1 at 135° | nearer sulfonyl O at 2.6 Å | signed volume at N2 below 0.001 in size | Solute frames, any | All-atom frames read | nearest oxygen at 2.6 Å from the zinc | water at the 2.5 Å of the place of the ring | water at 3.5 Å of N1 | water at 3.5 Å of the ring with N1 | All-atom frames, any |')
B('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for l in ALL:
    a_, b_ = SC[l]['threshold'], SF[l]['threshold']
    B(f'| {NAME[l]} | {th(SC[l]["n_solute_frames"])} | {a_["n_d1"]} | {a_["n_d2"]} | {a_["n_hb_distance"]} | {a_["n_hb_angle"]} | {a_["n_sulfonyl_O"]} | {a_["n_n2_sign"]} | {a_["n_frames_any"]} | '
      f'{th(SF[l]["n_all_atom_frames"])} | {b_["n_zinc_oxygen"]} | {b_["n_ring_place"]} | {b_["n_water_N1"]} | {b_["n_ring_wide"]} | {b_["n_frames_any"]} |')
B(f'| all {len(ALL)} | {th(sum(SC[l]["n_solute_frames"] for l in ALL))} | {tsol["n_d1"]} | {tsol["n_d2"]} | {tsol["n_hb_distance"]} | {tsol["n_hb_angle"]} | {tsol["n_sulfonyl_O"]} | {tsol["n_n2_sign"]} | {tsol["n_frames_any"]} | '
  f'{th(sum(SF[l]["n_all_atom_frames"] for l in ALL))} | {tful["n_zinc_oxygen"]} | {tful["n_ring_place"]} | {tful["n_water_N1"]} | {tful["n_ring_wide"]} | {tful["n_frames_any"]} |')
B('')
B('## H. The places of four crystal waters and the water at N2, whole production')
B('')
B('Own reading of the raw all-atom frames (`results/site_checks_full.json`): a water oxygen within 1.4 Å of the place of the crystal water after the fit of the Cα atoms, per cent of the 1,000 all-atom frames of the production; a water donating a hydrogen bond to N2 (O within 3.5 Å, O–H···N2 at least 135°); '
  'the same two for the frames in which the solute frame of the same time is crystal-like, rotated-like or intermediate. Waters 423, 516, 553 and 568 are the four removed to build start R.')
B('')
B('| Run | 423 | 516 | 553 | 568 | water donating to N2 | all-atom frames C / R / I | place 553 occupied in C / R / I frames (%) | water in the place of the crystal ring in C / R / I frames (%) |')
B('|---|---|---|---|---|---|---|---|---|')
nz = lambda v: 'none' if v is None else f(v, 1)
for l in ELEVEN:
    c_ = SF[l]['crystal_water_places']['production']
    b_ = SF[l]['by_state_of_the_solute_frame']
    B(f'| {NAME[l]} | {f(c_["423"], 1)} | {f(c_["516"], 1)} | {f(c_["553"], 1)} | {f(c_["568"], 1)} | {f(SF[l]["water_donating_to_N2"]["share_production"], 1)} | '
      f'{th(b_["C"]["n_all_atom_frames"])} / {b_["R"]["n_all_atom_frames"]} / {b_["I"]["n_all_atom_frames"]} | {nz(b_["C"]["place_553_occupied_share"])} / {nz(b_["R"]["place_553_occupied_share"])} / {nz(b_["I"]["place_553_occupied_share"])} | '
      f'{nz(b_["C"]["water_in_the_place_of_the_crystal_ring_share"])} / {nz(b_["R"]["water_in_the_place_of_the_crystal_ring_share"])} / {nz(b_["I"]["water_in_the_place_of_the_crystal_ring_share"])} |')
B('')
ttext = '\n'.join(T)
assert chr(0x2014) not in ttext and '{f(' not in ttext
(RES / 'tables.md').write_text(ttext, encoding='utf-8')
print('written: md/RESULTS.md (%d lines, %d words), md/results/tables.md (%d lines)' % (len(L), len(text.split()), len(T)))
