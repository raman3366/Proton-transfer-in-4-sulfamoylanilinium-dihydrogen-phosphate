"""d17: the molecular-dynamics campaign of md/ and the two docking tests (Glu106, four crystal waters), re-derived
with code of our own from data in git.  Reads only; writes derived/d17_md.json.  Nothing of md/analysis/,
md/production/ or dock/*.py is imported.

Molecular dynamics (md/RESULTS.md)
  series   md/results/series/<set>_<run>.csv.gz (md/production/make_series.py): one line per 10 ps with d1, d2, the
           ligand RMSD, the state, the signed volume at N2 and the hydrogen-bond flag, the released part of the
           equilibration at negative time.  From these the script forms, for each of the eleven runs: the state
           counts of the whole production, of its last half and of its last 20 ns, the verdict by the 90 % rule,
           the means of the ligand RMSD and of d1 and d2, the share of frames with the hydrogen bond to Thr199,
           the counts on the sign at N2, the first frames after the release (sign positive, not rotated-like,
           crystal-like), the time of the conversion by the plan's rule (first frame of the first visit of at
           least 1 ns in the opposite state on the joined series, assigned states as in section 4 of
           md/analysis/SPEC.md), the visits to the rotated-like state and the blocks of 20 ns.
  frames   md/results/final_frames/*.pdb.gz: d1 and d2 measured on the last frame against the last line of the
           series; atoms of the solute counted.
  all-atom md/results/runs/<set>/<run>/first_program_timeseries_full.csv (tracked): frames with any oxygen, a water
           oxygen or a protein oxygen within 2.6 A of the zinc.
  gate G-Z md/results/runs/<set>/X0/first_program_timeseries_equil_solute.csv (tracked): means over the ten frames
           of the first 100 ps.
  rehearsal md/test/rehearsal_2e/runs/TN0/analysis/timeseries_solute.csv (tracked): the rotated-like frames of the
           5 ns test run from start N, counted from d1 and d2.
  crystal  md/system/protein_md.pdb, zinc.pdb, ligand_X.pdb (the coordinates of 6RL9 that the systems are built
           from): the zinc out of the ring planes of His94, His96 and His119, S-C4, O1-S-C4, S-N2-Zn, C4-S-N2, Zn-N2
           and the nearer sulfonyl oxygen, the values of 6RL9 that the departures of the runs are set against.
  gate G-R md/qm/gr/gr_table.csv: the energies of every point relative to X from the Hartree columns, and the root
           mean square of part A (force field minus quantum chemistry, each family against its lowest
           quantum-chemical point; 23 comparisons), with and without the D3 term.
  beside   every summary value that the facts take is compared with the second program
           (second_program_summary.json) and, where the blind third reading has it
           (md/results/independent/third_reading/third_reading.json), with that reading.
  Not re-derived here (no per-frame data in git): the Zn-N distances and angles, C4-S-N2, tau_fold, the nearer
  sulfonyl oxygen and the C-alpha RMSD over the 10,000 solute frames, temperature and density of the state logs,
  gates G-L, G-V (solute part) and G-R2.  For these the script records the agreement of the readings that exist.

Docking tests (dock/GLU106.md, dock/WATER4.md)
  Every out.pdbqt of the runs in the reported receptors (T3_6RL9, T4_3HS4), in the receptors with Glu106 as a
  carboxylate (T3_6RL9_glu106GLU, T4_3HS4_glu106GLU) and in the receptors with the four crystal waters
  (T3_6RL9_wat4, T3_6RL9_glu106GLU_wat4, T3_6RL9_wat4rev, T3_6RL9_wat4bare) is parsed here: score of the top pose
  per seed; symmetry-corrected heavy-atom RMSD of every mode from the crystal ligand in the receptor frame (own
  graph isomorphism, as in d14) and, for the sulfanilamide species in 6RL9, from the rotated reference (seed 14,
  mode 3 of the reported production run); near-native (within 1.0 A), beyond 2.0 A, rotated pose (beyond 2.0 A and
  within 2.0 A of the rotated reference); Zn-N and Thr199 OG1...N of the top poses; score separation of the two
  orientations per seed; closest approach of a mode to a water site.
  The four waters: B factors and occupancies from dock/pdb/6RL9.cif.gz, their distances to the crystal ligand
  and to the docked rotated pose (md/system/ligand_R.pdb).
  Glu106: its pKa read from the raw PROPKA output (*.pka) of the protonation audit and of the docking preparation.
  Fixed-pose rescoring: the means are formed again from the per-pose values of rescore_fixed_poses*.json (the
  rescoring itself needs the docking program and is not repeated).
"""
import csv
import glob
import gzip
import itertools
import statistics as stt

import networkx as nx
from networkx.algorithms import isomorphism

from common import *

MD = REPO / 'md'
SER = MD / 'results' / 'series'
FRA = MD / 'results' / 'final_frames'
RUNS = [('bonded_2e', r, r) for r in ('X1', 'X2', 'X3', 'N1', 'N2', 'N3', 'R1', 'R2', 'R3')] + \
       [('bonded_2s', 'R1', 'R1_2s'), ('bonded_2e_n2fit', 'R1', 'R1_n2fit')]
START_STATE = {'X': 'C', 'N': 'C', 'R': 'R'}
OPP = {'C': 'R', 'R': 'C'}
M1, M2 = 6.7, 5.2                 # midpoints of d1 and d2 (plan, rule 2)
DWELL = 100                       # frames of a long visit (1 ns)
WINDOW = 2000                     # frames of the last 20 ns
out = {'md': {'runs': {}}, 'dock': {}}


# ====================================================================================================== series
def read_series(model, run):
    with gzip.open(SER / f'{model}_{run}.csv.gz', 'rt') as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r['t'] = float(r['time_ns'])
        for k in ('d1_A', 'd2_A', 'lig_rmsd_A', 'n2_volume'):
            r[k] = float(r[k])
        r['sign'] = (r['n2_volume'] > 0) - (r['n2_volume'] < 0)
        r['hb'] = int(r['hb_T199'])
        # the state formed again from d1 and d2 (strict inequalities, as the specification words them)
        r['state_own'] = 'C' if (r['d1_A'] < M1 and r['d2_A'] > M2) else 'R' if (r['d1_A'] > M1 and r['d2_A'] < M2) else 'I'
    return rows


def counts(rows):
    c = {s: sum(1 for r in rows if r['state_own'] == s) for s in 'CRI'}
    c['n'] = len(rows)
    return c


def pct(c, n):
    return 100.0 * c / n


def visits(states, start_state):
    """Assigned states (the state of a frame if C or R, else that of the frame before; before the first frame the
    state of the start) and the visits: [(state, first index, n frames)]."""
    assigned, prev = [], start_state
    for s in states:
        prev = s if s in 'CR' else prev
        assigned.append(prev)
    vis = []
    for k, s in enumerate(assigned):
        if vis and vis[-1][0] == s:
            vis[-1][2] += 1
        else:
            vis.append([s, k, 1])
    return assigned, vis


def run_from_series(model, run):
    rows = read_series(model, run)
    eq, pr = [r for r in rows if r['part'] != 'prod'], [r for r in rows if r['part'] == 'prod']
    start = START_STATE[run[0]]
    opp = OPP[start]
    start_sign = 1 if start == 'C' else -1
    q = dict(n_equilibration=len(eq), n_production=len(pr), start_state=start,
             state_differs_from_the_file=sum(1 for r in rows if r['state_own'] != r['state']),
             sign_differs_from_the_file=sum(1 for r in rows if r['sign'] != int(r['n2_sign'])))
    for name, part in (('whole', pr), ('last_half', pr[len(pr) // 2:]), ('last_20ns', pr[-WINDOW:])):
        c = counts(part)
        q[name] = dict(c, pct_C=pct(c['C'], c['n']), pct_R=pct(c['R'], c['n']), pct_I=pct(c['I'], c['n']))
    w = q['last_20ns']
    c_start, c_opp = w[start], w[opp]
    q['verdict'] = 'CONVERTED' if 10 * c_opp >= 9 * w['n'] else 'STAYED' if 10 * c_start >= 9 * w['n'] else 'MIXED'
    for col, key in (('lig_rmsd_A', 'rmsd'), ('d1_A', 'd1'), ('d2_A', 'd2')):
        x = np.array([r[col] for r in pr])
        q[key] = dict(mean=float(x.mean()), sd=float(x.std()), min=float(x.min()), max=float(x.max()),
                      mean_last_half=float(x[len(x) // 2:].mean()), mean_last_20ns=float(x[-WINDOW:].mean()))
    q['rmsd']['mean_crystal_like'] = float(np.mean([r['lig_rmsd_A'] for r in pr if r['state_own'] == 'C']))
    hb = sum(r['hb'] for r in pr)
    q['hb'] = dict(count=hb, n=len(pr), pct=pct(hb, len(pr)))
    # the sign at N2
    sg = [r['sign'] for r in rows]
    q['n2'] = dict(n_negative_production=sum(1 for r in pr if r['sign'] < 0),
                   pct_negative_production=pct(sum(1 for r in pr if r['sign'] < 0), len(pr)),
                   n_changes_from_the_release=sum(1 for a, b in zip(sg, sg[1:]) if a != b),
                   n_abs_volume_below_0p2_production=sum(1 for r in pr if abs(r['n2_volume']) < 0.2),
                   volume_min_from_the_release=min(r['n2_volume'] for r in rows),
                   volume_max_from_the_release=max(r['n2_volume'] for r in rows),
                   n_R_with_sign_positive_production=sum(1 for r in pr if r['state_own'] == 'R' and r['sign'] > 0),
                   n_C_with_sign_negative_production=sum(1 for r in pr if r['state_own'] == 'C' and r['sign'] < 0))

    def first(cond):
        for k, r in enumerate(rows):
            if cond(r):
                return dict(index=k, ps_from_release=10 * (k + 1), part=r['part'], state=r['state_own'])
        return None
    q['first_frame_sign_not_of_start'] = first(lambda r: r['sign'] == -start_sign)
    q['first_frame_not_in_state_of_start'] = first(lambda r: r['state_own'] != start)
    q['first_frame_in_opposite_state'] = first(lambda r: r['state_own'] == opp)
    last_start_eq = [k for k, r in enumerate(eq) if r['state_own'] == start]
    q['last_equilibration_frame_in_state_of_start_ps'] = 10 * (last_start_eq[-1] + 1) if last_start_eq else None
    ce = counts(eq)
    q['released_equilibration'] = ce
    # the conversion by the plan's rule, on the joined series
    assigned, vis = visits([r['state_own'] for r in rows], start)
    longv = [v for v in vis if v[2] >= DWELL]
    conv = next((v for v in longv if v[0] == opp), None)
    q['rule'] = dict(n_visits=len(vis), n_long_visits=len(longv),
                     time_from_release_ps=10 * (conv[1] + 1) if conv else None,
                     part=rows[conv[1]]['part'] if conv else None,
                     in_equilibration=(conv[1] < len(eq)) if conv else None,
                     n_frames_of_the_visit=conv[2] if conv else None,
                     released_frames_before=conv[1] if conv else None,
                     rotated_like_frames_before=sum(1 for r in rows[:conv[1]] if r['state_own'] == 'R') if conv else None)
    # visits of the production alone (assigned states), as md/RESULTS.md counts the returns of the candidate run
    a_pr, v_pr = visits([r['state_own'] for r in pr], start)
    vr = [v for v in v_pr if v[0] == 'R']
    q['production_visits_R'] = dict(n=len(vr), longest_frames=max((v[2] for v in vr), default=0),
                                    longest_ns=0.01 * max((v[2] for v in vr), default=0),
                                    longest_starts_ns=(0.01 * (max(vr, key=lambda v: v[2])[1] + 1)) if vr else None,
                                    frames_assigned=sum(v[2] for v in vr), n_long=sum(1 for v in vr if v[2] >= DWELL))
    q['blocks_20ns'] = []
    for b in range(5):
        blk = pr[b * WINDOW:(b + 1) * WINDOW]
        c = counts(blk)
        q['blocks_20ns'].append(dict(c, n2_negative=sum(1 for r in blk if r['sign'] < 0), hb=sum(r['hb'] for r in blk),
                                     rmsd_mean=float(np.mean([r['lig_rmsd_A'] for r in blk]))))
    q['last_line'] = dict(time_ns=rows[-1]['t'], d1=rows[-1]['d1_A'], d2=rows[-1]['d2_A'], state=rows[-1]['state_own'])
    return q


def final_frame(model, run):
    with gzip.open(FRA / f'{model}_{run}_final_solute.pdb.gz', 'rt') as fh:
        at = [(ln[12:16].strip(), ln[17:20].strip(), int(ln[22:26]),
               np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]), ln[76:78].strip())
              for ln in fh if ln.startswith('ATOM')]
    # residues are numbered in sequence: Phe131 is residue 129, Asn67 residue 66 (6RL9 has no residue 126)
    n1 = next(a[3] for a in at if a[1] == 'SN1' and a[0] == 'N1')
    cz = next(a[3] for a in at if a[1] == 'PHE' and a[2] == 129 and a[0] == 'CZ')
    od1 = next(a[3] for a in at if a[1] == 'ASN' and a[2] == 66 and a[0] == 'OD1')
    zn = next(a[3] for a in at if a[1] == 'ZN1')
    n2 = next(a[3] for a in at if a[1] == 'SN1' and a[0] == 'N2')
    return dict(n_atoms=len(at), n_ligand_atoms=sum(1 for a in at if a[1] == 'SN1'),
                n_zinc=sum(1 for a in at if a[1] == 'ZN1'), n_protein_atoms=sum(1 for a in at if a[1] not in ('SN1', 'ZN1')),
                n_residues_with_CA=sum(1 for a in at if a[0] == 'CA'), n_heavy_ligand=sum(1 for a in at if a[1] == 'SN1' and a[4] != 'H'),
                d1=dist(n1, cz), d2=dist(n1, od1), Zn_N2=dist(zn, n2))


def full_series(model, run):
    p = MD / 'results' / 'runs' / model / run / 'first_program_timeseries_full.csv'
    with open(p, newline='') as fh:
        rows = list(csv.DictReader(fh))
    n = len(rows)
    d = dict(n=n, file=rel(p))
    for col, key in (('d_Zn_O_any_min_A', 'any'), ('d_Zn_Owater_min_A', 'water'), ('d_Zn_Oprotein_min_A', 'protein'),
                     ('d_Zn_Osulfonyl_min_full_A', 'sulfonyl')):
        x = np.array([float(r[col]) for r in rows])
        d[key] = dict(count_within_2p6=int((x <= 2.6).sum()), pct_within_2p6=pct(int((x <= 2.6).sum()), n),
                      min=float(x.min()), mean=float(x.mean()))
    d['flag_any'] = sum(int(r['fifth_any']) for r in rows)
    ring = np.array([int(r['n_water_ring_place']) for r in rows])
    d['frames_with_water_in_the_place_of_the_ring'] = int((ring > 0).sum())
    return d


def gate_gz(model):
    p = MD / 'results' / 'runs' / model / 'X0' / 'first_program_timeseries_equil_solute.csv'
    with open(p, newline='') as fh:
        rows = [r for r in csv.DictReader(fh) if int(r['stage']) == 1 and int(r['extension']) == 0]
    cols = ['d_Zn_N2_A', 'd_Zn_H94NE2_A', 'd_Zn_H96NE2_A', 'd_Zn_H119ND1_A', 'ang_N2_Zn_H94NE2_deg', 'ang_N2_Zn_H96NE2_deg',
            'ang_N2_Zn_H119ND1_deg', 'ang_H94NE2_Zn_H96NE2_deg', 'ang_H94NE2_Zn_H119ND1_deg', 'ang_H96NE2_Zn_H119ND1_deg']
    d = dict(n_frames=len(rows), file=rel(p), means={c: float(np.mean([float(r[c]) for r in rows])) for c in cols})
    so = np.array([float(r['d_Zn_Osulfonyl_min_A']) for r in rows])
    d['sulfonyl_O_within_2p6'] = int((so <= 2.6).sum())
    d['sulfonyl_O_nearest'] = float(so.min())
    return d


def dig(d, path):
    for k in path.split('/'):
        d = d[k] if not isinstance(d, list) else d[int(k)]
    return d


SUMMARY_KEYS = [   # values of first_program_summary.json that the facts take: compared with the second program
    'stats/lig_rmsd_A/mean', 'stats/lig_rmsd_A/sd', 'stats/lig_rmsd_A/mean_last_half', 'stats/lig_rmsd_A/max',
    'stats/d1_N1_F131CZ_A/mean_last_half', 'stats/d2_N1_N67OD1_A/mean_last_half',
    'stats/d_Zn_N2_A/mean', 'stats/d_Zn_N2_A/sd', 'stats/d_Zn_H94NE2_A/mean', 'stats/d_Zn_H96NE2_A/mean',
    'stats/d_Zn_H119ND1_A/mean', 'stats/ang_N2_Zn_H94NE2_deg/mean', 'stats/ang_N2_Zn_H96NE2_deg/mean',
    'stats/ang_N2_Zn_H119ND1_deg/mean', 'stats/ang_H94NE2_Zn_H96NE2_deg/mean', 'stats/ang_H94NE2_Zn_H119ND1_deg/mean',
    'stats/ang_H96NE2_Zn_H119ND1_deg/mean', 'stats/ang_S_N2_Zn_deg/mean', 'stats/d_Zn_Osulfonyl_min_A/mean',
    'stats/tau_fold_deg/mean', 'stats/tau_fold_deg/sd', 'stats/ang_C4_S_N2_deg/mean', 'stats/ang_C4_S_N2_deg/sd',
    'stats/d_S_N2_A/mean', 'stats/ca_rmsd_A/mean', 'stats/ca_rmsd_A/mean_last_half', 'stats/ca_rmsd_A/max',
    'stats/d_Zn_E106OE1_A/min', 'fractions/hb_N2H_T199OG1/count', 'fractions/fifth_sulfonylO/count',
    'fractions/fifth_any/count', 'fractions/fifth_water/count', 'fractions/fifth_proteinO/count',
    'orientation/count_C', 'orientation/count_R', 'orientation/count_I', 'orientation/count_C_window',
    'orientation/count_R_window', 'orientation/count_I_window', 'n2_configuration/count_negative', 'rmsf/max_A']

third = json.loads((MD / 'results/independent/third_reading/third_reading.json').read_text())['runs']
conv_rec = {r['label']: r for r in json.loads((MD / 'results/conversion_times.json').read_text())['runs']}
worst_second, worst_third, n_second, n_third, mism = 0.0, 0.0, 0, 0, []
for model, run, rid in RUNS:
    q = run_from_series(model, run)
    q['final_frame'] = final_frame(model, run)
    q['final_frame']['d1_minus_series'] = q['final_frame']['d1'] - q['last_line']['d1']
    q['final_frame']['d2_minus_series'] = q['final_frame']['d2'] - q['last_line']['d2']
    q['all_atom'] = full_series(model, run)
    s1 = json.loads((MD / 'results/runs' / model / run / 'first_program_summary.json').read_text())
    s2 = json.loads((MD / 'results/runs' / model / run / 'second_program_summary.json').read_text())
    for k in SUMMARY_KEYS:
        a, b = dig(s1, k), dig(s2, k)
        n_second += 1
        worst_second = max(worst_second, abs(a - b))
    # this script against the first program
    own = [(q['whole']['C'], dig(s1, 'orientation/count_C')), (q['whole']['R'], dig(s1, 'orientation/count_R')),
           (q['last_20ns']['C'], dig(s1, 'orientation/count_C_window')), (q['last_20ns']['R'], dig(s1, 'orientation/count_R_window')),
           (q['hb']['count'], dig(s1, 'fractions/hb_N2H_T199OG1/count')),
           (q['n2']['n_negative_production'], dig(s1, 'n2_configuration/count_negative')),
           (q['all_atom']['any']['count_within_2p6'], dig(s1, 'fractions/fifth_any/count')),
           (q['all_atom']['water']['count_within_2p6'], dig(s1, 'fractions/fifth_water/count')),
           (q['all_atom']['protein']['count_within_2p6'], dig(s1, 'fractions/fifth_proteinO/count'))]
    ownf = [(q['rmsd']['mean'], dig(s1, 'stats/lig_rmsd_A/mean')), (q['rmsd']['sd'], dig(s1, 'stats/lig_rmsd_A/sd')),
            (q['rmsd']['mean_last_half'], dig(s1, 'stats/lig_rmsd_A/mean_last_half')),
            (q['d1']['mean_last_half'], dig(s1, 'stats/d1_N1_F131CZ_A/mean_last_half')),
            (q['d2']['mean_last_half'], dig(s1, 'stats/d2_N1_N67OD1_A/mean_last_half'))]
    if q['verdict'] != dig(s1, 'orientation/verdict') or any(a != b for a, b in own):
        mism.append(f'{model}/{run}: counts or verdict differ from the first program')
    q['largest_difference_from_the_first_program'] = max(abs(a - b) for a, b in ownf)
    c1 = conv_rec[f'{model}/{run}']['conversion']
    if (c1['time_from_release_ps'] if c1 else None) != q['rule']['time_from_release_ps']:
        mism.append(f'{model}/{run}: time of the rule differs from conversion_times.json')
    # the blind third reading
    t = third[f'{model}/{run}']
    tv = [(q['whole']['C'], t['states']['whole']['C']), (q['whole']['R'], t['states']['whole']['R']),
          (q['last_20ns']['C'], t['states']['last_20ns']['C']), (q['last_20ns']['R'], t['states']['last_20ns']['R']),
          (q['hb']['count'], t['hbond_N2H_T199OG1']['count']), (q['n2']['n_negative_production'], t['n2_configuration']['frames_minus']),
          (dig(s1, 'fractions/fifth_sulfonylO/count'), t['zinc']['sulfonyl_O_within_2.6A']['count'])]
    tf = [(q['rmsd']['mean'], t['lig_rmsd_A']['mean_whole']), (q['rmsd']['mean_last_half'], t['lig_rmsd_A']['mean_last_half']),
          (dig(s1, 'stats/d_Zn_N2_A/mean'), t['zinc']['mean_Zn_N2_A']), (dig(s1, 'stats/d_Zn_H94NE2_A/mean'), t['zinc']['mean_Zn_H94NE2_A']),
          (dig(s1, 'stats/d_Zn_H96NE2_A/mean'), t['zinc']['mean_Zn_H96NE2_A']), (dig(s1, 'stats/d_Zn_H119ND1_A/mean'), t['zinc']['mean_Zn_H119ND1_A']),
          (dig(s1, 'stats/d_Zn_Osulfonyl_min_A/mean'), t['zinc']['mean_Zn_Osulfonyl_min_A']),
          (dig(s1, 'stats/ca_rmsd_A/mean'), t['imaging_and_fit']['ca_rmsd_to_reference_A']['mean']),
          (dig(s1, 'stats/ca_rmsd_A/max'), t['imaging_and_fit']['ca_rmsd_to_reference_A']['max']),
          (dig(s1, 'stats/d1_N1_F131CZ_A/mean_last_half'), t['d1_A']['last_half']['mean']),
          (dig(s1, 'stats/d2_N1_N67OD1_A/mean_last_half'), t['d2_A']['last_half']['mean'])]
    n_third += len(tv) + len(tf) + 2
    worst_third = max([worst_third] + [abs(a - b) for a, b in tf])
    if any(a != b for a, b in tv) or t['verdict']['verdict'] != q['verdict'] or \
            (t['conversion'].get('t_from_release_ps') if t['conversion'].get('has_conversion_time') else None) != \
            (float(q['rule']['time_from_release_ps']) if q['rule']['time_from_release_ps'] else None):
        mism.append(f'{model}/{run}: a count, the verdict or the time of the rule differs from the third reading')
    q['third_reading'] = dict(Zn_N2=t['zinc']['mean_Zn_N2_A'], Zn_H94NE2=t['zinc']['mean_Zn_H94NE2_A'],
                              Zn_H96NE2=t['zinc']['mean_Zn_H96NE2_A'], Zn_H119ND1=t['zinc']['mean_Zn_H119ND1_A'],
                              Zn_Osulfonyl_min_mean=t['zinc']['mean_Zn_Osulfonyl_min_A'],
                              sulfonyl_O_within_2p6=t['zinc']['sulfonyl_O_within_2.6A']['count'],
                              ca_rmsd_mean=t['imaging_and_fit']['ca_rmsd_to_reference_A']['mean'],
                              ca_rmsd_max=t['imaging_and_fit']['ca_rmsd_to_reference_A']['max'],
                              rmsd_mean=t['lig_rmsd_A']['mean_whole'], rmsd_mean_last_half=t['lig_rmsd_A']['mean_last_half'],
                              rmsd_mean_last_20ns=t['lig_rmsd_A']['mean_last_20ns'], rmsd_max=t['lig_rmsd_A']['max'])
    out['md']['runs'][rid] = q

out['md']['comparison'] = dict(
    second_program=dict(n_values=n_second, largest_difference=worst_second, keys=SUMMARY_KEYS),
    third_reading=dict(n_values=n_third, largest_difference_of_a_mean=worst_third,
                       note='counts, verdicts and times of the rule equal; means are written with four decimals there'),
    this_script_against_the_first_program=dict(
        largest_difference_of_a_mean=max(q['largest_difference_from_the_first_program'] for q in out['md']['runs'].values()),
        frames_whose_state_formed_here_differs_from_the_file=sum(q['state_differs_from_the_file'] for q in out['md']['runs'].values()),
        frames_whose_sign_formed_here_differs_from_the_file=sum(q['sign_differs_from_the_file'] for q in out['md']['runs'].values())),
    mismatches=mism)

# the two validation runs: all-atom frames and gate G-Z
out['md']['X0'] = {}
for model, tag in (('bonded_2e', '2e'), ('bonded_2s', '2s')):
    out['md']['X0'][tag] = dict(all_atom=full_series(model, 'X0'), gate_GZ=gate_gz(model))

# ranges over groups of runs that the text names
R_ = out['md']['runs']


def rng(ids, fn):
    v = [fn(R_[i]) for i in ids]
    return [min(v), max(v)]


G_X, G_N, G_R = ['X1', 'X2', 'X3'], ['N1', 'N2', 'N3'], ['R1', 'R2', 'R3']
NINE = G_X + G_N + G_R
out['md']['ranges'] = dict(
    rmsd_last_half={'X': rng(G_X, lambda q: q['rmsd']['mean_last_half']), 'N': rng(G_N, lambda q: q['rmsd']['mean_last_half']),
                    'R': rng(G_R, lambda q: q['rmsd']['mean_last_half']), 'XN': rng(G_X + G_N, lambda q: q['rmsd']['mean_last_half']),
                    'nine': rng(NINE, lambda q: q['rmsd']['mean_last_half'])},
    rmsd_max={'XN': rng(G_X + G_N, lambda q: q['rmsd']['max']), 'nine': rng(NINE, lambda q: q['rmsd']['max'])},
    hb_pct={'XN': rng(G_X + G_N, lambda q: q['hb']['pct']), 'R': rng(G_R, lambda q: q['hb']['pct']),
            'nine': rng(NINE, lambda q: q['hb']['pct'])},
    rule_ps={'R': rng(G_R, lambda q: q['rule']['time_from_release_ps'])},
    sign_ps={'R': rng(G_R, lambda q: q['first_frame_sign_not_of_start']['ps_from_release'])},
    not_rotated_ps={'R': rng(G_R, lambda q: q['first_frame_not_in_state_of_start']['ps_from_release'])},
    n_rotated_like_production={'nine': sum(R_[i]['whole']['R'] for i in NINE), 'R': sum(R_[i]['whole']['R'] for i in G_R),
                               'XN': sum(R_[i]['whole']['R'] for i in G_X + G_N)},
    n_frames_production={'nine': sum(R_[i]['whole']['n'] for i in NINE), 'R': sum(R_[i]['whole']['n'] for i in G_R),
                         'XN': sum(R_[i]['whole']['n'] for i in G_X + G_N)},
    n_abs_volume_below_0p2={'nine': rng(NINE, lambda q: q['n2']['n_abs_volume_below_0p2_production'])},
    pct_any_O_within_2p6={'nine': rng(NINE, lambda q: q['all_atom']['any']['pct_within_2p6'])},
    n_water_or_protein_O_within_2p6=sum(q['all_atom']['water']['count_within_2p6'] + q['all_atom']['protein']['count_within_2p6']
                                        for q in R_.values()))

# the rehearsal run from start N (a test of the pipeline with the topology of N1; no result): its rotated-like frames
with open(MD / 'test/rehearsal_2e/runs/TN0/analysis/timeseries_solute.csv', newline='') as fh:
    tn0 = list(csv.DictReader(fh))
tn0_R = [float(r['time_ns']) for r in tn0 if float(r['d1_N1_F131CZ_A']) > M1 and float(r['d2_N1_N67OD1_A']) < M2]
out['md']['rehearsal_TN0'] = dict(n_frames=len(tn0), n_rotated_like=len(tn0_R), times_ns=tn0_R,
                                  file=rel(MD / 'test/rehearsal_2e/runs/TN0/analysis/timeseries_solute.csv'))

# the same quantities in the crystal structure, from the coordinates of 6RL9 that the systems are built from
def _pdb_atoms(path):
    return [(ln[12:16].strip(), int(ln[22:26]), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
            for ln in path.read_text().splitlines() if ln.startswith(('ATOM', 'HETATM'))]


_prot, _zn = _pdb_atoms(MD / 'system/protein_md.pdb'), _pdb_atoms(MD / 'system/zinc.pdb')[0][2]
_lx = {a[0]: a[2] for a in _pdb_atoms(MD / 'system/ligand_X.pdb')}


def _out_of_plane(resnum):
    """Signed distance of the zinc from the least-squares plane of the imidazole ring (CG, ND1, CD2, CE1, NE2); the
    sense of the normal is that of (ND1 - CG) x (CD2 - CG), as in md/production/site_checks.py."""
    ring = np.array([next(a[2] for a in _prot if a[1] == resnum and a[0] == nme) for nme in ('CG', 'ND1', 'CD2', 'CE1', 'NE2')])
    c = ring.mean(axis=0)
    nrm = np.linalg.svd(ring - c)[2][2]
    return float((_zn - c) @ nrm * np.sign(nrm @ np.cross(ring[1] - ring[0], ring[2] - ring[0])))


out['md']['xtal_6RL9'] = dict(
    files=[rel(MD / 'system/protein_md.pdb'), rel(MD / 'system/zinc.pdb'), rel(MD / 'system/ligand_X.pdb')],
    zinc_out_of_the_ring_plane_A={f'His{k}': _out_of_plane(k) for k in (94, 96, 119)},
    zn_oop_his119=abs(_out_of_plane(119)), d_S_C4_A=dist(_lx['S'], _lx['C4']), ang_O1_S_C4_deg=angle(_lx['O1'], _lx['S'], _lx['C4']),
    ang_S_N2_Zn_deg=angle(_lx['S'], _lx['N2'], _zn), d_Zn_N2_A=dist(_zn, _lx['N2']), ang_C4_S_N2_deg=angle(_lx['C4'], _lx['S'], _lx['N2']),
    d_Zn_Osulfonyl_min_A=min(dist(_zn, _lx['O1']), dist(_zn, _lx['O2'])))

# ------------------------------------------------------------------------------------------------- gate G-R
with open(MD / 'qm/gr/gr_table.csv', newline='') as fh:
    GT = {r['id']: r for r in csv.DictReader(fh)}
e_d3 = {k: (float(r['QM_hartree']) + float(r['D3_hartree'])) * HA_KCAL for k, r in GT.items()}
e_no = {k: float(r['QM_hartree']) * HA_KCAL for k, r in GT.items()}
e_mm = {k: float(r['E_MM']) for k, r in GT.items()}


def part_a(eq):
    """r.m.s. of force field minus quantum chemistry over the points of part A, each family (base X, base Rx)
    against its lowest quantum-chemical point."""
    dev = {}
    for base in ('X', 'Rx'):
        fam = [k for k, r in GT.items() if 'A' in r['part'].split('+') and r['base'] == base]
        ref = min(fam, key=lambda k: eq[k])
        for k in fam:
            if k != ref:
                dev[k] = (e_mm[k] - e_mm[ref]) - (eq[k] - eq[ref])
    v = np.array(list(dev.values()))
    return dict(n_compared=len(v), rms=float(np.sqrt((v ** 2).mean())), largest=float(v[np.abs(v).argmax()]),
                largest_at=list(dev)[int(np.abs(v).argmax())], n_points=len(v) + 2)


def rel_x(eq, k):
    return eq[k] - eq['X']


out['md']['gate_GR'] = dict(
    with_d3=part_a(e_d3), without_d3=part_a(e_no),
    points={k: dict(qm=rel_x(e_d3, k), mm=e_mm[k] - e_mm['X'], mm_minus_qm=(e_mm[k] - e_mm['X']) - rel_x(e_d3, k),
                    qm_without_d3=rel_x(e_no, k), mm_minus_qm_without_d3=(e_mm[k] - e_mm['X']) - rel_x(e_no, k))
            for k in ('a4', 'd1', 'e_p030', 'Rx', 'Nx')},
    n_points_in_the_table=len(GT), source=rel(MD / 'qm/gr/gr_table.csv'))


# ====================================================================================================== docking
AD2EL = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O', 'O': 'O', 'S': 'S', 'SA': 'S',
         'H': 'H', 'HD': 'H', 'HS': 'H', 'P': 'P'}
RCOV = {'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07, 'H': 0.31}


def pdbqt_models(path):
    """[(score, [(element, xyz), ...]), ...] in the order of the file (mode 1 first), heavy atoms only."""
    models, cur, score = [], None, None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = []
        elif ln.startswith('REMARK VINA RESULT'):
            score = float(ln.split()[3])
        elif ln.startswith(('ATOM', 'HETATM')) and cur is not None:
            el = AD2EL[ln[77:79].strip()]
            if el != 'H':
                cur.append((el, np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
        elif ln.startswith('ENDMDL'):
            models.append((score, cur))
            cur = None
    return models


def graph(atoms):
    G = nx.Graph()
    for i, (e, x) in enumerate(atoms):
        G.add_node(i, el=e)
    for i, j in itertools.combinations(range(len(atoms)), 2):
        if dist(atoms[i][1], atoms[j][1]) < RCOV[atoms[i][0]] + RCOV[atoms[j][0]] + 0.45:
            G.add_edge(i, j)
    return G


def sdf_atoms(path):
    lines = Path(path).read_text().splitlines()
    na = int(lines[3][:3])
    at = [(l[31:34].strip(), np.array([float(l[0:10]), float(l[10:20]), float(l[20:30])])) for l in lines[4:4 + na]]
    return [a for a in at if a[0] != 'H']


_ISO = {}


def mappings(ref, pose, key):
    """All element-preserving isomorphisms reference atom -> pose atom (cached per reference and ligand)."""
    if key not in _ISO:
        gm = isomorphism.GraphMatcher(graph(ref), graph(pose), node_match=lambda a, b: a['el'] == b['el'])
        _ISO[key] = [np.array([m[i] for i in range(len(ref))]) for m in gm.isomorphisms_iter()]
        assert _ISO[key], f'no isomorphism for {key}'
    return _ISO[key]


def rmsd_sym(ref, pose, key):
    R = np.array([a[1] for a in ref])
    P = np.array([a[1] for a in pose])
    return float(min(np.sqrt(((P[m] - R) ** 2).sum(1).mean()) for m in mappings(ref, pose, key)))


def receptor_atoms(target):
    zn = og1 = None
    waters = {}
    for ln in (DOCK / target / 'receptor.pdb').read_text().splitlines():
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        name, resn, resi = ln[12:16].strip(), ln[17:20].strip(), int(ln[22:26])
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        if name == 'ZN':
            zn = xyz
        elif resn == 'THR' and resi == 199 and name == 'OG1':
            og1 = xyz
        elif resn == 'HOH' and name == 'O':
            waters[resi] = xyz
    return zn, og1, waters


def sulfonamide_n(atoms):
    """Index of the sulfonamide nitrogen: the nitrogen whose only heavy neighbour is a sulfur."""
    G = graph(atoms)
    c = [i for i, (e, _) in enumerate(atoms) if e == 'N' and [atoms[j][0] for j in G[i]] == ['S']]
    assert len(c) == 1, c
    return c[0]


_REC, _RUN = {}, {}
ROT_REF = pdbqt_models(DOCK / 'T3_6RL9/runs/sfa_neutral__ad4zn__prod/seed14/out.pdbqt')[2][1]     # mode 3
WAT_SITES = None


def analyse(target, run, ref_file, frame_6rl9):
    """One run of one receptor from its out.pdbqt files."""
    if (target, run) in _RUN:
        return _RUN[(target, run)]
    if target not in _REC:
        _REC[target] = receptor_atoms(target)
    zn, og1, _ = _REC[target]
    ref = sdf_atoms(DOCK / ref_file)
    lig = run.split('__')[0]
    seeds = sorted(glob.glob(str(DOCK / target / 'runs' / run / 'seed*' / 'out.pdbqt')))
    per, modes_all = [], []
    for s in seeds:
        ms = pdbqt_models(s)
        rows = []
        for k, (sc, at) in enumerate(ms):
            m = dict(mode=k + 1, score=sc, rmsd=rmsd_sym(ref, at, (ref_file, lig)))
            if frame_6rl9 and lig.startswith('sfa_'):
                m['rmsd_rot'] = rmsd_sym(ROT_REF, at, ('rot', lig))
                m['cls'] = ('near-native' if m['rmsd'] <= 1.0 else
                            'rotated pose' if (m['rmsd'] > 2.0 and m['rmsd_rot'] <= 2.0) else
                            'other beyond 2.0' if m['rmsd'] > 2.0 else 'between')
                X = np.array([a[1] for a in at])
                dw = np.linalg.norm(X[:, None, :] - WAT_SITES[None, :, :], axis=2)
                m['min_dist_water_site'] = float(dw.min())
                m['n_sites_within_2p4'] = int((dw.min(0) <= 2.4).sum())
            else:
                m['cls'] = 'near-native' if m['rmsd'] <= 1.0 else 'beyond 2.0' if m['rmsd'] > 2.0 else 'between'
            rows.append(m)
        top = ms[0][1]
        i_n = sulfonamide_n(top)
        near = [m for m in rows if m['rmsd'] <= 1.0]
        far = [m for m in rows if m['rmsd'] > 2.0]
        rot = [m for m in rows if m['cls'] == 'rotated pose']
        d = dict(seed=int(Path(s).parent.name[4:]), score_top=ms[0][0], rmsd_top=rows[0]['rmsd'], cls_top=rows[0]['cls'],
                 n_modes=len(rows), zn_n=dist(zn, top[i_n][1]), og1_n=dist(og1, top[i_n][1]),
                 rmsd_best_mode=min(m['rmsd'] for m in rows),
                 rank_near=near[0]['mode'] if near else None,
                 gap=(near[0]['score'] - rows[0]['score']) if near else None,
                 sep=abs(min(m['score'] for m in near) - min(m['score'] for m in far)) if near and far else None,
                 near_minus_far=(min(m['score'] for m in near) - min(m['score'] for m in far)) if near and far else None,
                 rot_behind_near=(min(m['score'] for m in rot) - min(m['score'] for m in near)) if near and rot else None)
        per.append(d)
        modes_all += [dict(m, seed=d['seed']) for m in rows]
    n = len(per)
    sc = [p['score_top'] for p in per]
    rm = [p['rmsd_top'] for p in per]

    def nn(key):
        return [p[key] for p in per if p[key] is not None]
    q = dict(target=target, run=run, reference=ref_file, n_seeds=n,
             score=dict(mean=stt.mean(sc), sd=stt.stdev(sc), per_seed=sc),
             rmsd=dict(mean=stt.mean(rm), sd=stt.stdev(rm), min=min(rm), max=max(rm), n_le_1=sum(r <= 1.0 for r in rm),
                       n_le_2=sum(r <= 2.0 for r in rm), f_le_1=sum(r <= 1.0 for r in rm) / n,
                       f_le_2=sum(r <= 2.0 for r in rm) / n, per_seed=rm,
                       best_mode_min=min(p['rmsd_best_mode'] for p in per), best_mode_max=max(p['rmsd_best_mode'] for p in per)),
             zn_n=dict(mean=stt.mean(p['zn_n'] for p in per), n_within_2p3=sum(p['zn_n'] <= 2.3 for p in per),
                       min=min(p['zn_n'] for p in per), max=max(p['zn_n'] for p in per)),
             og1_n=dict(mean=stt.mean(p['og1_n'] for p in per), n_within_3p3=sum(p['og1_n'] <= 3.3 for p in per)),
             top=dict(near_native=sum(p['cls_top'] == 'near-native' for p in per),
                      beyond_2=sum(p['rmsd_top'] > 2.0 for p in per),
                      between=sum(1.0 < p['rmsd_top'] <= 2.0 for p in per),
                      rotated_pose=sum(p['cls_top'] == 'rotated pose' for p in per)),
             modes=dict(n=len(modes_all), near_native=sum(m['rmsd'] <= 1.0 for m in modes_all),
                        beyond_2=sum(m['rmsd'] > 2.0 for m in modes_all),
                        between=sum(1.0 < m['rmsd'] <= 2.0 for m in modes_all),
                        rotated_pose=sum(m['cls'] == 'rotated pose' for m in modes_all),
                        seeds_holding_rotated_pose=len({m['seed'] for m in modes_all if m['cls'] == 'rotated pose'}),
                        seeds_holding_near_native=len({m['seed'] for m in modes_all if m['rmsd'] <= 1.0}),
                        best_score_near_native=min((m['score'] for m in modes_all if m['rmsd'] <= 1.0), default=None),
                        best_score_rotated_pose=min((m['score'] for m in modes_all if m['cls'] == 'rotated pose'), default=None)),
             orientation=dict(worst_rank_near=max(nn('rank_near'), default=None),
                              gap=[min(nn('gap')), max(nn('gap'))] if nn('gap') else None,
                              sep=[min(nn('sep')), max(nn('sep'))] if nn('sep') else None, n_seeds_with_both=len(nn('sep')),
                              near_minus_far=[min(nn('near_minus_far')), max(nn('near_minus_far'))] if nn('near_minus_far') else None,
                              near_minus_far_mean=stt.mean(nn('near_minus_far')) if nn('near_minus_far') else None,
                              rot_behind_near=[min(nn('rot_behind_near')), max(nn('rot_behind_near'))] if nn('rot_behind_near') else None),
             # the parts of the two criteria that need no PoseBusters run (analyse.py: mean RMSD, half of the seeds,
             # both site checks passed by half of the seeds)
             gate_part=dict(strict=bool(stt.mean(rm) <= 1.0 and sum(r <= 1.0 for r in rm) / n >= 0.5
                                        and sum(p['zn_n'] <= 2.3 for p in per) / n >= 0.5
                                        and sum(p['og1_n'] <= 3.3 for p in per) / n >= 0.5),
                            two_A=bool(stt.mean(rm) <= 2.0 and sum(r <= 2.0 for r in rm) / n >= 0.5
                                       and sum(p['zn_n'] <= 2.3 for p in per) / n >= 0.5
                                       and sum(p['og1_n'] <= 3.3 for p in per) / n >= 0.5)))
    if frame_6rl9 and lig.startswith('sfa_'):
        q['water_sites'] = dict(min_dist_any_mode=min(m['min_dist_water_site'] for m in modes_all),
                                modes_within_2p4=sum(m['n_sites_within_2p4'] > 0 for m in modes_all))
    _RUN[(target, run)] = q
    return q


# the four waters: coordinates as deposited (the receptors are in the frame of the deposited file)
_, _, wsites = receptor_atoms('T3_6RL9_wat4')
WATERS = sorted(wsites)
assert WATERS == [423, 516, 553, 568], WATERS
WAT_SITES = np.array([wsites[w] for w in WATERS])

# ---------------------------------------------------------------------------------------------------- Glu106
glu = {}
for base, var in (('T3_6RL9', 'T3_6RL9_glu106GLU'), ('T4_3HS4', 'T4_3HS4_glu106GLU')):
    comp = json.loads((DOCK / var / 'comparison.json').read_text())
    for run, rec in comp['runs'].items():
        ref = rec['baseline']['reference']['reference']
        a = analyse(base, run, f'{base}/{ref}', base == 'T3_6RL9')
        b = analyse(var, run, f'{base}/{ref}', base == 'T3_6RL9')
        glu[f"{base.split('_')[1]}/{run}"] = dict(acid=a, carboxylate=b)
out['dock']['glu106'] = glu
vv = [(k, abs(v['carboxylate']['score']['mean'] - v['acid']['score']['mean'])) for k, v in glu.items() if '__vin' in k]
out['dock']['glu106_vina'] = dict(
    n_runs=len(vv), n_runs_ad4zn=sum(1 for k in glu if '__ad4zn' in k),
    max_abs_score_change_sfa=max(d for k, d in vv if '/sfa_neutral' in k),
    max_abs_score_change_sfa_species=max(d for k, d in vv if '/sfa_' in k),
    max_abs_score_change_azm=max(d for k, d in vv if '/azm_' in k),
    zn_n_within_2p3_any_vina_run=sum(v[s]['zn_n']['n_within_2p3'] for k, v in glu.items() if k.endswith('__vina') or '__vina__' in k
                                     for s in ('acid', 'carboxylate')))
# acetazolamide ahead of sulfanilamide, validation runs (score of sulfanilamide minus score of acetazolamide)
out['dock']['glu106_azm_minus_sfa'] = {
    pdb: {s: glu[f'{pdb}/sfa_neutral__ad4zn'][s]['score']['mean'] - glu[f'{pdb}/azm_neutral__ad4zn'][s]['score']['mean']
          for s in ('acid', 'carboxylate')} for pdb in ('6RL9', '3HS4')}

# ---------------------------------------------------------------------------------------------- four waters
wat = {}
for dry, var, rule, state in (('T3_6RL9', 'T3_6RL9_wat4', 'hbond', 'acid'),
                              ('T3_6RL9_glu106GLU', 'T3_6RL9_glu106GLU_wat4', 'hbond', 'carboxylate'),
                              ('T3_6RL9', 'T3_6RL9_wat4rev', 'rev', 'acid'), ('T3_6RL9', 'T3_6RL9_wat4bare', 'bare', 'acid')):
    comp = json.loads((DOCK / var / 'comparison_water.json').read_text())
    for run in comp['runs']:
        wat[f'{state}/{rule}/{run}'] = dict(dry=analyse(dry, run, 'T3_6RL9/ref_SAN.sdf', True),
                                           water=analyse(var, run, 'T3_6RL9/ref_SAN.sdf', True))
out['dock']['wat4'] = wat
out['dock']['wat4_counts'] = dict(
    n_runs=len(wat), n_runs_ad4zn=sum(1 for k in wat if '__ad4zn' in k),
    modes_of_the_rotated_pose_in_any_water_run=sum(v['water']['modes']['rotated_pose'] for v in wat.values()),
    modes_of_the_rotated_pose_in_the_ad4zn_water_runs=sum(v['water']['modes']['rotated_pose'] for k, v in wat.items() if '__ad4zn' in k),
    vina_zn_n_mean_range=[min(v[s]['zn_n']['mean'] for k, v in wat.items() if '__vin' in k for s in ('dry', 'water')),
                          max(v[s]['zn_n']['mean'] for k, v in wat.items() if '__vin' in k for s in ('dry', 'water'))],
    vina_seeds_within_1A=sum(v[s]['rmsd']['n_le_1'] for k, v in wat.items() if '__vin' in k for s in ('dry', 'water')),
    vina_seeds_zn_within_2p3=sum(v[s]['zn_n']['n_within_2p3'] for k, v in wat.items() if '__vin' in k for s in ('dry', 'water')))

# the four waters in the deposited entry, and their distances to the two placements of the ligand
import gemmi  # noqa: E402

with gzip.open(DOCK / 'pdb' / '6RL9.cif.gz', 'rt') as fh:
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(fh.read()).sole_block())
hoh, san, allB = {}, [], []
for ch in st[0]:
    for res in ch:
        if res.name == 'HOH':
            a = res[0]
            hoh[res.seqid.num] = dict(B=a.b_iso, occ=a.occ, xyz=np.array([a.pos.x, a.pos.y, a.pos.z]))
            allB += [x.b_iso for x in res]            # every oxygen position, alternate positions included
        elif res.name == 'SAN' and res.seqid.num == 304:
            san = [(a.name, np.array([a.pos.x, a.pos.y, a.pos.z])) for a in res if a.element.name != 'H']
ligR = [(ln[12:16].strip(), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
        for ln in (MD / 'system/ligand_R.pdb').read_text().splitlines()
        if ln.startswith(('ATOM', 'HETATM')) and ln[76:78].strip() != 'H' and not ln[12:16].strip().startswith('H')]
wrec = {}
for w in WATERS:
    dx = [(dist(hoh[w]['xyz'], x), nme) for nme, x in san]
    dr = [(dist(hoh[w]['xyz'], x), nme) for nme, x in ligR]
    wrec[str(w)] = dict(number_in_the_entry=next(k for k in hoh if k == w), B=hoh[w]['B'], occupancy=hoh[w]['occ'], to_crystal_ligand=min(dx)[0], to_crystal_ligand_atom=min(dx)[1],
                        to_rotated_pose=min(dr)[0], to_rotated_pose_atom=min(dr)[1],
                        deposited_minus_receptor_file=dist(hoh[w]['xyz'], wsites[w]))
out['dock']['waters'] = dict(sites=wrec, n_waters_in_the_entry=len(hoh), n_water_oxygen_positions=len(allB),
                             median_B_all_waters=float(np.median(allB)),
                             B_range_of_the_four=[min(v['B'] for v in wrec.values()), max(v['B'] for v in wrec.values())],
                             n_heavy_atoms_crystal_ligand=len(san), n_heavy_atoms_rotated_pose=len(ligR))

# pKa of Glu106 in the raw PROPKA output: the two runs of the protonation audit of the simulation system (without and
# with the zinc in the input, 6RL9) and the run of the docking preparation of each receptor (no zinc in the input)
pka = {}
for f in ('md/system/audit/propka_without_zn.pka', 'md/system/audit/propka_with_zn.pka', 'dock/T3_6RL9/propka.pka', 'dock/T4_3HS4/propka.pka'):
    rows_ = [ln.split() for ln in (REPO / f).read_text().splitlines() if ln.startswith('GLU 106 A') and len(ln.split()) > 5 and ln.split()[5] == '%']
    summ_ = [ln.split() for ln in (REPO / f).read_text().splitlines() if ln.split()[:3] == ['GLU', '106', 'A'] and len(ln.split()) == 5]
    assert len(summ_) == 1 and len(rows_) == 1 and float(summ_[0][3]) == float(rows_[0][3]), (f, summ_, rows_[:1])
    pka[f] = float(summ_[0][3])
out['dock']['propka'] = pka

# fixed poses: the means formed again from the per-pose values of the rescoring records
fix = {}
for var in ('T3_6RL9_glu106GLU', 'T4_3HS4_glu106GLU'):
    r = json.loads((DOCK / var / 'rescore_fixed_poses.json').read_text())
    for run, rec in r['runs'].items():
        ch = {k: [m['change'] for m in rec['modes'] if m['kind'] == k] for k in ('near-native', 'ring-rotated')}
        fix[f"glu106/{var.split('_')[1]}/{run}"] = dict(
            near_native=dict(n=len(ch['near-native']), mean=stt.mean(ch['near-native']), min=min(ch['near-native']), max=max(ch['near-native'])),
            ring_rotated=dict(n=len(ch['ring-rotated']), mean=stt.mean(ch['ring-rotated']), min=min(ch['ring-rotated']), max=max(ch['ring-rotated'])),
            rotated_minus_near=stt.mean(ch['ring-rotated']) - stt.mean(ch['near-native']))
    r = json.loads((DOCK / var / 'rescore_fixed_poses_vina.json').read_text())
    allc = [p['change'] for rec in r['runs'].values() for p in rec['poses']]
    fix[f"glu106_vina/{var.split('_')[1]}"] = dict(n_poses=len(allc), n_runs=len(r['runs']), max_abs_change=max(abs(c) for c in allc))
for var in ('T3_6RL9_wat4', 'T3_6RL9_glu106GLU_wat4', 'T3_6RL9_wat4rev', 'T3_6RL9_wat4bare', 'T3_6RL9_wat4fixed'):
    r = json.loads((DOCK / var / 'rescore_fixed_poses.json').read_text())
    rec = r['runs']['sfa_neutral__ad4zn__prod']
    nn_ = [m['water_maps'] - m['dry_maps'] for m in rec['modes'] if m['kind'] == 'near-native']
    rp = [m['water_maps'] for m in rec['modes'] if m['kind'] == 'rotated pose']
    fix[f'wat4/{var}'] = dict(near_native=dict(n=len(nn_), change_mean=stt.mean(nn_), change_min=min(nn_), change_max=max(nn_)),
                              rotated_pose=dict(n=len(rp), mean=stt.mean(rp), min=min(rp), max=max(rp)))
out['dock']['fixed_poses'] = fix
out['dock']['fixed_poses_vina_total'] = dict(
    n_poses=sum(v['n_poses'] for k, v in fix.items() if k.startswith('glu106_vina/')),
    max_abs_change=max(v['max_abs_change'] for k, v in fix.items() if k.startswith('glu106_vina/')))
out['dock']['fixed_poses_rotated_range_all_rules'] = [min(v['rotated_pose']['min'] for k, v in fix.items() if k.startswith('wat4/')),
                                                     max(v['rotated_pose']['max'] for k, v in fix.items() if k.startswith('wat4/'))]

# ---------------------------------------------------------------------------- against the records of the tests
cmp_ = []
for base, var in (('T3_6RL9', 'T3_6RL9_glu106GLU'), ('T4_3HS4', 'T4_3HS4_glu106GLU')):
    comp = json.loads((DOCK / var / 'comparison.json').read_text())
    for run, rec in comp['runs'].items():
        for side, s in (('baseline', 'acid'), ('variant', 'carboxylate')):
            own = glu[f"{base.split('_')[1]}/{run}"][s]
            cmp_.append(dict(run=f'{var}:{run}:{s}',
                             d_score=abs(own['score']['mean'] - stt.mean(rec[side]['score']['per_seed'])),
                             d_rmsd=abs(own['rmsd']['mean'] - stt.mean(rec[side]['rmsd']['per_seed'])),
                             d_rmsd_seed=max(abs(a - b) for a, b in zip(own['rmsd']['per_seed'], rec[side]['rmsd']['per_seed'])),
                             d_zn=abs(round(own['zn_n']['mean'], 2) - rec[side]['site_checks']['Zn-N(sulfonamide)']['distance_mean']),
                             pass_equal=bool(abs(own['zn_n']['n_within_2p3'] / own['n_seeds'] - rec[side]['site_checks']['Zn-N(sulfonamide)']['pass_fraction']) < 1e-9
                                             and abs(own['og1_n']['n_within_3p3'] / own['n_seeds'] - rec[side]['site_checks']['Thr199 OG1 H-bond']['pass_fraction']) < 1e-9
                                             and own['rmsd']['n_le_1'] == round(rec[side]['rmsd']['f_le_1A'] * own['n_seeds'])
                                             and own['top']['near_native'] == rec[side]['orientation']['top_near_native']
                                             and own['top']['beyond_2'] == rec[side]['orientation']['top_ring_rotated'])))
for dry, var, rule, state in (('T3_6RL9', 'T3_6RL9_wat4', 'hbond', 'acid'),
                              ('T3_6RL9_glu106GLU', 'T3_6RL9_glu106GLU_wat4', 'hbond', 'carboxylate'),
                              ('T3_6RL9', 'T3_6RL9_wat4rev', 'rev', 'acid'), ('T3_6RL9', 'T3_6RL9_wat4bare', 'bare', 'acid')):
    comp = json.loads((DOCK / var / 'comparison_water.json').read_text())
    for run, rec in comp['runs'].items():
        for side in ('dry', 'water'):
            own = wat[f'{state}/{rule}/{run}'][side]
            ms = rec[side]['mode_summary']
            cmp_.append(dict(run=f'{var}:{run}:{side}',
                             d_score=abs(own['score']['mean'] - stt.mean(rec[side]['score']['per_seed'])),
                             d_rmsd=abs(own['rmsd']['mean'] - stt.mean(rec[side]['rmsd']['per_seed'])),
                             d_rmsd_seed=max(abs(a - b) for a, b in zip(own['rmsd']['per_seed'], rec[side]['rmsd']['per_seed'])),
                             d_zn=abs(round(own['zn_n']['mean'], 2) - rec[side]['site_checks']['Zn-N(sulfonamide)']['distance_mean']),
                             pass_equal=bool(abs(own['zn_n']['n_within_2p3'] / own['n_seeds'] - rec[side]['site_checks']['Zn-N(sulfonamide)']['pass_fraction']) < 1e-9
                                             and abs(own['og1_n']['n_within_3p3'] / own['n_seeds'] - rec[side]['site_checks']['Thr199 OG1 H-bond']['pass_fraction']) < 1e-9),
                             modes_equal=bool(own['modes']['n'] == ms['n_modes'] and own['modes']['near_native'] == ms['modes']['near_native']
                                              and own['modes']['rotated_pose'] == ms['modes']['rotated_pose']
                                              and own['modes']['beyond_2'] == ms['modes']['ring_rotated_paper']
                                              and own['top']['near_native'] == ms['top_pose']['near_native']),
                             d_water=abs(round(own['water_sites']['min_dist_any_mode'], 2) - ms['min_dist_water_site']['all_modes'])))
out['dock']['against_the_records'] = dict(
    n_runs_compared=len(cmp_), max_abs_score_mean_diff=max(c['d_score'] for c in cmp_),
    max_abs_rmsd_mean_diff=max(c['d_rmsd'] for c in cmp_), max_abs_rmsd_seed_diff=max(c['d_rmsd_seed'] for c in cmp_),
    max_abs_zn_n_mean_diff_at_two_decimals=max(c['d_zn'] for c in cmp_),
    mode_counts_equal_in_every_water_comparison=all(c.get('modes_equal', True) for c in cmp_),
    seed_counts_equal_in_every_comparison=all(c['pass_equal'] for c in cmp_),
    max_abs_water_site_distance_diff_at_two_decimals=max(c.get('d_water', 0.0) for c in cmp_),
    runs_with_a_difference=[c['run'] for c in cmp_ if c['d_score'] > 1e-9 or c['d_rmsd'] > 6e-4 or c['d_zn'] > 0.0101
                            or not c.get('modes_equal', True) or not c['pass_equal']])

dump('d17_md.json', out)
m = out['md']['comparison']
a = out['dock']['against_the_records']
print(f"MD: {len(R_)} runs; mismatches {len(m['mismatches'])}; largest difference from the first program "
      f"{m['this_script_against_the_first_program']['largest_difference_of_a_mean']:.1e}; second program "
      f"{m['second_program']['largest_difference']:.1e} over {m['second_program']['n_values']} values; third reading "
      f"{m['third_reading']['largest_difference_of_a_mean']:.1e}; G-R rms {out['md']['gate_GR']['with_d3']['rms']:.4f} "
      f"({out['md']['gate_GR']['without_d3']['rms']:.4f} without D3)")
print(f"docking: {a['n_runs_compared']} run records compared; max |d score| {a['max_abs_score_mean_diff']:.1e}; "
      f"max |d RMSD mean| {a['max_abs_rmsd_mean_diff']:.1e}; runs with a difference {len(a['runs_with_a_difference'])}")
