#!/usr/bin/env python3
"""Assembles md/PRODUCTION.md from the records (every number is read from a file)."""
import csv, hashlib, json, re, statistics as st, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
R = REPO / 'md'
BULK = Path('/root/cqc-sfa-md')
RUNS = BULK / 'runs'
S = Path(sys.argv[1])                       # scratch folder with frags.md, x1/...
IST = timezone(timedelta(hours=5, minutes=30))


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def ist(dtobj):
    return dtobj.astimezone(IST).strftime('%-d %b %Y, %-I:%M %p IST')


J = lambda p: json.load(open(p))
B = {m: {s: J(R / 'topology' / m / s / 'build.json') for s in ('X', 'N', 'R') if (R / 'topology' / m / s / 'build.json').exists()}
     for m in ('bonded_2e', 'bonded_2s', 'bonded_2e_n2fit', 'bonded_2s_n2fit')}
GZ = {s: J(R / 'production' / f'gate_GZ_bonded_{s}.json') for s in ('2s', '2e')}
SM = {s: J(RUNS / f'bonded_{s}' / 'X0' / 'analysis' / 'summary.json') for s in ('2s', '2e')}
SM2 = {s: J(RUNS / f'bonded_{s}' / 'X0' / 'analysis_recheck' / 'summary.json') for s in ('2s', '2e')}
RR = {s: J(RUNS / f'bonded_{s}' / 'X0' / 'analysis_recheck' / 'recheck_report.json') for s in ('2s', '2e')}
MAN = {s: J(RUNS / f'bonded_{s}' / 'X0' / 'manifest.json') for s in ('2s', '2e')}
SB = J(R / 'production' / 'seminario_blocks.json')
AQ = J(BULK / 'qm/gr/after_qm_report_2s/after_qm_check_2s.json')
frags = (S / 'frags.md').read_text()
gates = (S / 'gates.md').read_text()


def table_after(text, first_cell, n=1):
    """The n-th Markdown table of text whose header starts with first_cell."""
    lines, out, k, on = text.splitlines(), [], 0, False
    for ln in lines:
        if ln.startswith('| ' + first_cell + ' |') and not on:
            k += 1
            on = k == n
        if on:
            if not ln.startswith('|'):
                break
            out.append(ln)
    assert out, first_cell
    return '\n'.join(out)


ok = lambda b: 'met' if b else '**NOT MET**'
g = {s: SM[s]['gate_GV'] for s in SM}
fl = {s: g[s]['GV3_zinc']['criteria']['fifth_ligand'] for s in g}
sol = {s: list(csv.DictReader(open(RUNS / f'bonded_{s}' / 'X0' / 'analysis' / 'timeseries_solute.csv'))) for s in SM}
full = {s: list(csv.DictReader(open(RUNS / f'bonded_{s}' / 'X0' / 'analysis' / 'timeseries_full.csv'))) for s in SM}
osul = {s: [float(r['d_Zn_Osulfonyl_min_A']) for r in sol[s]] for s in sol}
hits_full = {s: [r for r in full[s] if r['fifth_any'] == '1'] for s in full}
z2s = B['bonded_2s']['X']['zinc_site_bonded']
ch = z2s['charges']
prim, other = '2e', '2s'
assert g['2e']['verdict'] == 'PASSED' and g['2s']['verdict'] == 'FAILED'
assert all(SM2[s]['gate_GV']['verdict'] == g[s]['verdict'] for s in g)
assert all(RR[s]['agree'] is True and RR[s]['n_runs_disagree'] == 0 for s in RR)
nval = {s: RR[s]['runs'][0]['comparison'] for s in RR}


def find_num(d, *names):
    for n in names:
        if n in d:
            return d[n]
    raise KeyError(names)


ncomp = {s: nval[s]['n_values'] for s in nval}
assert all(nval[s]['n_disagree'] == 0 and nval[s]['n_threshold'] == 0 and nval[s]['n_keys_one_side'] == 0 and nval[s]['agree'] for s in nval)
import os
_fe = {os.path.basename(k): v for k, v in B['bonded_2e']['X']['mcpb']['files_sha256'].items()}
_fs = {os.path.basename(k): v for k, v in B['bonded_2s']['X']['mcpb']['files_sha256'].items()}
assert [k for k in _fe if _fe[k] != _fs[k]] == ['CA2_mcpbpy.frcmod'], 'the work folders differ in more than the frcmod'
for _t in SB['step_2s']['terms']:
    assert abs(_t['own_code']['k_block_as_mcpb'] - _t['k_mcpb']) < 0.006, _t
    assert all(b['computed'] and not b['eigenvalues_complex'] for b in _t['blocks_read']), _t
_hb = [t for t in SB['step_2s']['terms'] if t['kind'] == 'bond' and 'SAN261@N2' not in t['atoms']]
assert len(_hb) == 3 and all(abs(t['own_code']['k_block_as_mcpb'] - t['own_code']['k_other_block']) < 0.3 for t in _hb)
_n2b = [t for t in SB['step_2s']['terms'] if t['kind'] == 'bond' and 'SAN261@N2' in t['atoms']][0]['own_code']
assert (_n2b['k_block_as_mcpb'], _n2b['k_other_block'], _n2b['k_mean_of_the_two']) == (102.44, 99.53, 100.99)
for _s in MAN:
    _pl = MAN[_s]['minimisation']['platform']
    assert _pl['name'] == 'CUDA' and _pl['properties_actual']['DeviceName'] == 'NVIDIA L40S' and _pl['properties_actual']['Precision'] == 'mixed'
    assert MAN[_s]['equilibration']['n_extension_steps'] == 0 and MAN[_s]['hydration_check']['skipped'] is True
_sol_hits = {float(r['time_ns']) for r in sol['2s'] if r['fifth_sulfonylO'] == '1'}
assert all(float(r['time_ns']) in _sol_hits for r in hits_full['2s'])
assert all(float(r['d_Zn_O2_A']) <= 2.6 for r in sol['2s'] if r['fifth_sulfonylO'] == '1')
assert sum(int(r['fifth_water']) + int(r['fifth_proteinO']) for r in full['2s']) == 0
assert all(not B[m][s]['checks_failed'] for m in B for s in B[m])
assert len({B[m][s]['script_sha256'] for m in B for s in B[m]}) == 1 and B['bonded_2e']['X']['script_sha256'] == sha(R / 'build_topology.py')
keep_ist = datetime.fromtimestamp(os.path.getmtime(BULK / 'status/KEEP_RUNNING'), IST).strftime('%-I:%M %p IST on %-d %b %Y')
tasks_x0 = J(R / 'gpu/tasks_x0.json')
tasks_p = J(R / 'gpu/tasks_production.json')
tasks_s = J(R / 'gpu/tasks_sensitivity.json')

# ---------------------------------------------------------------- additions of 7 Oct 2026, after the independent check of the
# set-up (md/PLAN.md, entry of 7 Oct 2026, 12:10 AM IST, points 2 and 3; paper/audit/research/md_production_check_2026-10-07.json)
# (1) crystal water in the systems: md/system/decisions.json, md/SYSTEM.md section 4, the eight build.json
DEC = J(R / 'system/decisions.json')
_w, _wst = DEC['water'], DEC['starts']
_w510 = [x for x in _w['labelled_or_partial'] if x['water'] == 510]
_wr = _wst['R']['waters_removed']
assert (_w['sites_in_entry'], _w['kept'], _w['dropped']) == (216, 215, [510]) and len(_w510) == 1 and _w510[0]['kept'] is False
assert [x['water'] for x in _wr] == [423, 516, 553, 568] and not _wst['X']['waters_removed'] and not _wst['N']['waters_removed']
assert (_wst['X']['n_waters'], _wst['N']['n_waters'], _wst['R']['n_waters']) == (_w['kept'], _w['kept'], _w['kept'] - len(_wr))
assert 'within 2.2 A of a ligand heavy atom is removed' in _w['ligand_rule']
assert all(B[m][s]['composition']['crystal_water_molecules'] == _wst[s]['n_waters'] for m in B for s in B[m])
assert 'All four are fully occupied with B 11.5–17.6 Å² (median of all water 22.5 Å²)' in (R / 'SYSTEM.md').read_text()
_wd = [f"{x['dist_A']:.2f} Å ({x['ligand_atom']})" for x in _wr]
_wr_txt = ', '.join(_wd[:-1]) + ' and ' + _wd[-1]
# (2) set 2e: the equilibrium values of the zinc bonds and angles are the values of 6RL9 (start X) to the digits of the
# parameter file, so the means of gate G-V are compared with the terms' own minima; the largest shifts of the means
_z2e, _c2e = B['bonded_2e']['X']['zinc_site_bonded'], g['2e']['GV3_zinc']['criteria']
_zb = {b['bond']: b for b in _z2e['bonds_to_zinc_against_the_parameter_files']}
_za = {a['angle']: a for a in _z2e['angles_with_zinc_against_the_parameter_files'] if a['kind'] == 'N-Zn-N'}
_gvb = {'d_Zn_N2_A': ('Zn-SN1304:N2', 'Zn–N2'), 'd_Zn_H94NE2_A': ('Zn-HD194:NE2', 'Zn–NE2 His94'),
        'd_Zn_H96NE2_A': ('Zn-HD296:NE2', 'Zn–NE2 His96'), 'd_Zn_H119ND1_A': ('Zn-HE1119:ND1', 'Zn–ND1 His119')}
_gva = {'ang_N2_Zn_H94NE2_deg': ('HD194:NE2-Zn-SN1304:N2', 'N2–Zn–NE2 His94'),
        'ang_N2_Zn_H96NE2_deg': ('HD296:NE2-Zn-SN1304:N2', 'N2–Zn–NE2 His96'),
        'ang_N2_Zn_H119ND1_deg': ('HE1119:ND1-Zn-SN1304:N2', 'N2–Zn–ND1 His119'),
        'ang_H94NE2_Zn_H96NE2_deg': ('HD194:NE2-Zn-HD296:NE2', 'NE2 His94–Zn–NE2 His96'),
        'ang_H94NE2_Zn_H119ND1_deg': ('HD194:NE2-Zn-HE1119:ND1', 'NE2 His94–Zn–ND1 His119'),
        'ang_H96NE2_Zn_H119ND1_deg': ('HD296:NE2-Zn-HE1119:ND1', 'NE2 His96–Zn–ND1 His119')}
assert len(_zb) == 4 and len(_za) == 6 and {v[0] for v in _gvb.values()} == set(_zb) and {v[0] for v in _gva.values()} == set(_za)
assert all(abs(b['r0_A'] - b['distance_in_6RL9_A']) < 0.0006 for b in _zb.values())
assert all(abs(a['theta0_deg'] - a['angle_in_this_start_deg']) < 0.006 for a in _za.values())
_zns = B['bonded_2e']['X']['terms_at_N2']['in_the_topology']['angles']['Zn-N2-S']
assert abs(_zns['theta0_deg'] - _zns['angle_in_this_start_deg']) < 0.006
_sh_bond = max((abs(_c2e[k]['mean_A'] - _zb[b]['r0_A']), n) for k, (b, n) in _gvb.items())
_sh_ang = max((abs(_c2e[k]['mean_deg'] - _za[a]['theta0_deg']), n) for k, (a, n) in _gva.items())
_sh_zns = _zns['theta0_deg'] - st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol['2e'])
# the same three numbers worked out by hand from the tables of sections 2.2, 2.5 and 5.3 of the file as it stood on 6 Oct 2026:
# 2.047101 - 2.0210 = 0.026 A; 108.92 - 105.617480 = 3.30 deg; 117.55 - 113.99 = 3.56 deg
assert (f'{_sh_bond[0]:.3f}', _sh_bond[1]) == ('0.026', 'Zn–NE2 His96') and (f'{_sh_ang[0]:.2f}', _sh_ang[1]) == ('3.30', 'N2–Zn–NE2 His94')
assert f'{_sh_zns:.2f}' == '3.56'
# the references of gate G-V for the six angles are these equilibrium values rounded to 0.1 degree, and the deviation of the
# gate (the bracket in the table of section 5.3) is taken against the rounded reference: 108.9 - 105.617480 = 3.28, not 3.30
assert all(abs(_c2e[k]['reference_6RL9_deg'] - round(_za[a]['theta0_deg'], 1)) < 1e-9
           and abs(_c2e[k]['deviation_deg'] - (_c2e[k]['mean_deg'] - _c2e[k]['reference_6RL9_deg'])) < 1e-5 for k, (a, n) in _gva.items())
# (3) the parameter files of both zinc sets in the repository are the files of the work folders
_frc = {s: B[f'bonded_{s}']['X']['mcpb']['files_sha256'][B[f'bonded_{s}']['X']['mcpb']['dir'] + '/CA2_mcpbpy.frcmod'] for s in ('2e', '2s')}
assert all(sha(R / f'qm/gr/charges/CA2_mcpbpy_{s}.frcmod') == _frc[s] for s in _frc), 'a parameter file in md/qm/gr/charges is not that of its work folder'

# ---------------------------------------------------------------- the start of production (inputs made while watching)
x1log = (S / 'x1/run_md.log.machine').read_text().splitlines() if (S / 'x1/run_md.log.machine').exists() else []
segs = []
for ln in x1log:
    m = re.search(r'\[(.*?)\] (seg\d+) \((\d+)-(\d+) ns\): done in (\d+) s, ([\d.]+) ns/day', ln)
    if m:
        segs.append({'when': m.group(1), 'seg': m.group(2), 'to_ns': int(m.group(4)), 'wall_s': int(m.group(5)), 'ns_day': float(m.group(6))})
eq_lines = [ln for ln in x1log if re.search(r'\] (min|eq\d+)', ln)]
watch = J(S / 'x1/watch.json') if (S / 'x1/watch.json').exists() else {}

P = []
A = P.append
A(f"""# Production: parameter sets, built systems, validation runs and the queue

Written on 6 Oct 2026 (IST) by the integration task of the evening; state of {watch.get('written_ist', 'the time of writing')}.
The contract is the plan's entry of 6 Oct 2026, 9:15 PM IST (`md/PLAN.md`): production with the terms at N2 as
built, by the user's decision. Every number below is read from a file named beside it; this file is assembled
from those files by a script, not typed. Times are IST (UTC+5:30). What was added on 7 Oct 2026 is marked where
it stands and listed at the end of the file.

## 1. Outcome in short

| | |
|---|---|
| Primary zinc set | **2e** (empirical constants, MCPB.py step 2e), systems `bonded_2e` |
| Why not 2s | Set 2s (Seminario) passed gate G-Z, but its validation run X0 **failed gate G-V**: part d of GV3, an oxygen within 2.6 Å of the zinc in {fl['2s']['full_trajectory']['count']} of {fl['2s']['full_trajectory']['n']} all-atom frames ({100 * fl['2s']['full_trajectory']['fraction']:.0f} %; at most 5 % allowed). Every other number of the gate is met (section 5) |
| Set 2e | passed gate G-Z and gate G-V with every number (sections 4 and 5; what parts a to c of GV3 can show for this set: section 5.3, paragraph added on 7 Oct 2026) |
| The two analysis programs | agree on every value of both runs: {ncomp['2s']:,} values compared for X0 of 2s and {ncomp['2e']:,} for X0 of 2e, 0 disagreements, 0 threshold frames, exit status 0 both times |
| Queued | nine runs X1 N1 R1, X2 N2 R2, X3 N3 R3 with `bonded_2e` (100 ns each), then two sensitivity runs of start R: `R_zn` (zinc set 2s) and `R_n2` (set 2e with the candidate terms at N2 that gate G-R2 did not accept) |
| Terms at N2 | as built, the same file in both zinc sets (sha256 `{B['bonded_2e']['X']['terms_at_N2']['sha256']}`); gates G-R and G-R2 on these terms failed (`md/GR.md`, `md/N2.md`) |

The plan's rule (entry of 9:15 PM, point 1): "the Seminario set (`bonded_2s`) is the primary one, provided it
passes gate G-Z and its validation run X0 passes gate G-V; if it fails either, that is reported and the empirical
set (`bonded_2e`) … is primary, again only after its own X0." That rule was applied as written. The failing part
of set 2s is narrow (three frames where two are allowed; the same event is seen in
{fl['2s']['solute_trajectory']['count']} of {fl['2s']['solute_trajectory']['n']} of the ten times denser solute frames,
{100 * fl['2s']['solute_trajectory']['fraction']:.1f} %); it is reported with all its numbers in section 5.3 and the rule
was not read again after the result.

One thing was done ahead of the rule: **X0 of both sets was put into the queue at once** (X0 of 2s first, X0 of
2e behind it), before the result of 2s was known. `md/RUN.md` gives every simulated zinc model its own
validation run, set 2e is simulated in either case (as the primary set or as the sensitivity set), and the two
runs take 12 minutes each. Nothing of production was queued before both gates were read.

## 2. Parameter sets

### 2.1 Where they come from

| Set | MCPB.py step | Work folder (Pi) | `CA2_mcpbpy.frcmod` sha256 | Made |
|---|---|---|---|---|
| 2e | 2e, empirical zinc constants | `{B['bonded_2e']['X']['mcpb']['dir']}` | `{B['bonded_2e']['X']['mcpb']['files_sha256'][B['bonded_2e']['X']['mcpb']['dir'] + '/CA2_mcpbpy.frcmod']}` | 6 Oct 2026, 7:10 AM IST (date of the folder), unchanged since |
| 2s | 2s, Seminario's method on the Hessian of the optimised small model | `{B['bonded_2s']['X']['mcpb']['dir']}` | `{B['bonded_2s']['X']['mcpb']['files_sha256'][B['bonded_2s']['X']['mcpb']['dir'] + '/CA2_mcpbpy.frcmod']}` | 6 Oct 2026, 9:18 PM IST, this task |

Set 2s was made with `bash md/qm/after_qm.sh 2s /root/cqc-sfa-md/qm/gr/mcpb_2s <report folder>` in a copy of the
step-1 files of `md/qm/mcpb/` with `CA2_small_opt.fchk` (sha256 `{sha(R / 'qm/mcpb/CA2_small_opt.fchk')}`, the
file the machine wrote at 8:47 PM IST). `md/qm/after_qm.sh` (sha256 `{sha(R / 'qm/after_qm.sh')[:16]}…`) and the
guard `md/qm/mcpb_step2s.py` (sha256 `{sha(R / 'qm/mcpb_step2s.py')[:16]}…`) were used as they are; **neither
needed a change**. The guard's line: "36 values passed to MCPB.py's averaging; 18 were of complex type, largest
|imaginary|/|value| 0.0e+00 (real part taken)". The four hard checks of the script passed
(`/root/cqc-sfa-md/qm/gr/after_qm_report_2s/after_qm_check_2s.json`):

""")
for h in AQ['hard']:
    A(f"- {h['check']}: {'passed' if h['passed'] else 'FAILED'}\n")
A(f"""
The two sets differ in the four zinc bonds and the fourteen angles with the zinc in `CA2_mcpbpy.frcmod`, and in
nothing else: the residue templates (`HD1.mol2`, `HD2.mol2`, `HE1.mol2`, `ZN1.mol2`, `SN1.mol2`, with the
equal-charge RESP set), `CA2_tleap.in` and `CA2_mcpbpy.pdb` of the two work folders are identical byte for byte
(sha256 in `md/topology/bonded_2e/X/build.json` and `md/topology/bonded_2s/X/build.json`, key `mcpb.files_sha256`).
Step 2ms (modified Seminario) is tabulated beside them and not used (plan, 6 Oct 2026, 3:02 AM, point 3).

### 2.2 The zinc terms of both sets (and 2ms beside them)

Amber form: bonds U = k (r − r0)², k in kcal mol⁻¹ Å⁻²; angles U = k (θ − θ0)², k in kcal mol⁻¹ rad⁻².
The columns 2e and 2s are read back from the built topologies of start X and equal the parameter files
(`build.json`, key `zinc_site_bonded`, check `zinc_bonds_and_angles_as_the_parameter_files`); the column 2ms is
from the frcmod that step 2ms writes (`md/production/seminario_blocks.json`). "In 6RL9" is the value in the
crystal start X.

**Bonds to the zinc**

{table_after(frags, 'Bond')}

**Angles with the zinc** (2ms has the equilibrium values of 2s)

{table_after(frags, 'Angle')}

**The two angle lines at N2 that MCPB.py writes are not in force** in any system: the file of the terms at N2 is
loaded last and replaces them (section 2.5). For the record:

{table_after(frags, 'Types')}

Torsions through the zinc: {z2s['torsions_through_zinc']['n']} in each set, every one with zero force constant;
no improper torsion holds the zinc (`torsions_through_zinc` in both `build.json`).

### 2.3 What step 2s read from the Hessian

The Hessian is a partial one: the six atoms that were displaced are
{', '.join(a['atom'] + ' (' + str(a['number']) + ')' for a in SB['displaced_atoms_found_in_the_matrix'])} of the
{SB['n_atoms']}-atom small model (MCPB.py's numbering: residues 93, 95 and 118 are His94, His96 and His119 of 6RL9, 260 is the
zinc and 261 the ligand), so {SB['n_blocks_non_zero']} of the {SB['n_blocks']} 3 × 3 blocks are computed
and the others are zero. `md/production/seminario_blocks.py` runs MCPB.py's own step 2s (and 2ms) in a
temporary copy with the reading of the Hessian instrumented and records every block that is read
(`md/production/seminario_blocks.json`):

- step 2s reads {SB['step_2s']['n_bonds']} bonds and {SB['step_2s']['n_angles']} angles; for a bond a–b the block
  [a, b], for an angle a–b–c the blocks [a, b] and [c, b] (the averaging over both directions is off);
- **every block read for every bond and every angle is a computed one: {str(SB['step_2s']['every_block_read_was_computed']).lower()}**
  (each has a displaced atom on at least one side; none is a zero block); the same holds for step 2ms
  ({str(SB['step_2ms']['every_block_read_was_computed']).lower()});
- the constants that the instrumented run gives equal the frcmod of the work folder
  ({str(SB['step_2s_here_equals_the_frcmod_of_the_work_folder']).lower()}), and a second calculation with code of
  its own (same blocks) gives the same numbers to the printed digits.

So nothing blocks set 2s on the side of the Hessian.

### 2.4 Are the constants of 2s physical?

In size, yes: bonds {min(b['k_kcal_per_mol_A2'] for b in z2s['bonds_to_zinc'])} to {max(b['k_kcal_per_mol_A2'] for b in z2s['bonds_to_zinc'])} kcal mol⁻¹ Å⁻², angles at the zinc
{min(a['k_kcal_per_mol_rad2'] for a in z2s['angles_with_zinc_against_the_parameter_files'] if a['kind'] == 'N-Zn-N'):.2f} to {max(a['k_kcal_per_mol_rad2'] for a in z2s['angles_with_zinc_against_the_parameter_files'] if a['kind'] == 'N-Zn-N'):.2f} and Zn–N–C
{min(a['k_kcal_per_mol_rad2'] for a in z2s['angles_with_zinc_against_the_parameter_files'] if a['kind'] == 'Zn-N-C'):.2f} to {max(a['k_kcal_per_mol_rad2'] for a in z2s['angles_with_zinc_against_the_parameter_files'] if a['kind'] == 'Zn-N-C'):.2f} kcal mol⁻¹ rad⁻²; no negative
or near-zero constant, no complex eigenvalue with an imaginary part, and MCPB.py's rules of thumb hold (bond
constants below 200 and lengths below 2.8 Å; angle constants below 100 and angles above 100°;
`after_qm_check_2s.json`, key `mcpb_rules_of_thumb`). Three things have to be said plainly:

1. **The equilibrium values of 2s are those of the small model optimised in vacuum, not of the crystal.** Two of
   them are far from 6RL9: N2–Zn–NE2 His96 has θ0 = 100.64° against 114.28° in the crystal (−13.64°), and
   N2–Zn–NE2 His94 118.19° against 108.92° (+9.27°); the two Zn–N–C angles of His96 are 132.81° and 120.86°
   against 122.91° and 129.31°. The gate asks for means within 10° of 6RL9, so set 2s was at risk at that
   angle before any run; in X0 the mean came out at {g['2s']['GV3_zinc']['criteria']['ang_N2_Zn_H96NE2_deg']['mean_deg']:.2f}° ({g['2s']['GV3_zinc']['criteria']['ang_N2_Zn_H96NE2_deg']['deviation_deg']:+.2f}°, inside the window).
2. **The bond Zn–NE2 His96 is soft and long in 2s**: k = 51.2 with r0 = 2.0840 Å (2e: 73.3, 2.0210 Å; crystal
   2.0206 Å). In X0 its mean is {g['2s']['GV3_zinc']['criteria']['d_Zn_H96NE2_A']['mean_A']:.3f} Å, inside the window that ends at 2.12 Å.
3. **Seminario's method is not unique for Zn–N2**: with the block [Zn, N2], which MCPB.py reads, k = 102.44;
   with the block [N2, Zn] 99.53 (mean 100.99). For the three histidine bonds the two blocks agree within 0.3.

What followed in the run is in section 5.3: with set 2s the sulfonyl oxygen O2 is nearer the zinc than with
set 2e.

### 2.5 The terms at N2

Both zinc sets are built with **the same file of the terms at N2, as built**
(`md/topology/bonded_2e/<start>/n2_terms.frcmod`, sha256 `{B['bonded_2s']['X']['terms_at_N2']['sha256']}`), loaded after every
file of MCPB.py. Read back from the topology of each of the six systems (`build.json`, key
`terms_at_N2.in_the_topology.entries_of_the_file`; check `terms_at_N2_are_those_of_the_n2_file`):

| Term | Types | k | Equilibrium or phase |
|---|---|---|---|
| angle Zn–N2–S | M1-Y4-sy | 50.00 kcal mol⁻¹ rad⁻² | 117.55° |
| angle Zn–N2–H | M1-Y4-hn | 50.00 kcal mol⁻¹ rad⁻² | 107.407° (its value in start X, not the 120° of step 2e) |
| angle S–N2–H | hn-Y4-sy | 61.550 kcal mol⁻¹ rad⁻² | 118.059° |
| torsion H–N2–S–C4 | ca-sy-Y4-hn | 0.2624 kcal mol⁻¹, n = 1 | 180° |
| | | 2.703 kcal mol⁻¹, n = 2 | 0° |

These terms did not pass gate G-R (`md/GR.md`) and the refit did not pass gate G-R2 (`md/N2.md`); they are used
by the user's decision of 6 Oct 2026 (plan, entry of 9:15 PM).

**The candidate of the sensitivity run** (`md/qm/n2/n2_terms_new_charges_as_built.frcmod`, sha256
`{B['bonded_2e_n2fit']['R']['terms_at_N2']['sha256']}`; its first line says "NOT ACCEPTED"): read back from the topology of
`bonded_2e_n2fit/R`, {sum(1 for e in B['bonded_2e_n2fit']['R']['terms_at_N2']['in_the_topology']['entries_of_the_file'] if e['as_the_file'])} of {len(B['bonded_2e_n2fit']['R']['terms_at_N2']['in_the_topology']['entries_of_the_file'])} entries of the file are in force as written:

| Entry of the file | In force (k and equilibrium; for torsions k, n, phase) |
|---|---|
""")
for e in B['bonded_2e_n2fit']['R']['terms_at_N2']['in_the_topology']['entries_of_the_file']:
    v = e.get('in_the_file') or e.get('in_the_file_k_n_phase') or e.get('in_the_file_k_phase_n')
    if 'in_the_file_k_n_phase' in e:
        txt = '; '.join(f"{k:g} kcal mol⁻¹, n = {n}, {ph:g}°" for k, n, ph in sorted(v, key=lambda t: t[1]))
    elif 'in_the_file_k_phase_n' in e:
        txt = f"{v[0]:g} kcal mol⁻¹, n = {v[2]}, {v[1]:g}° (N2 the central atom)"
    else:
        txt = f"{v[0]:g} kcal mol⁻¹ rad⁻², {v[1]:g}°"
    A(f"| {e['entry']} | {txt} |\n")
A(f"""
### 2.6 The charges of the site

The equal-charge RESP set (`md/qm/resp_equiv.py`; the six pairs of chemically equal ligand atoms made equal), the
same in every system of this file: the residue templates of the two work folders are identical, and the charges
of the five site residues in each topology equal the templates (`build.json`, key `zinc_site_bonded.charges`;
checks `site_charges_are_those_of_the_templates`, `site_charge_plus_one`, `ligand_equal_atoms_have_equal_charges`).
Record of the fit: `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/resp_equiv_2s.json`, sha256 `{ch['resp_equiv_json_sha256']}`.

{table_after(frags, 'Residue')}

{table_after(frags, 'Atom')}

## 3. The built systems

Builder: `md/build_topology.py`, sha256 `{sha(R / 'build_topology.py')}` (the same file built all eight systems:
`script_sha256` in each `build.json`). Built on the Pi on 6 Oct 2026 between
9:27 PM and 9:31 PM IST with

```bash
cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
export AMBERHOME=/root/cqc-sfa-tools/md OMP_NUM_THREADS=2 OPENMM_CPU_THREADS=2
PY="nice -n 10 /root/cqc-sfa-tools/md/bin/python"
for s in X N R; do
  $PY md/build_topology.py --model bonded --set 2e --start $s --note "…"
  $PY md/build_topology.py --model bonded --set 2s --start $s --n2-frcmod md/topology/bonded_2e/$s/n2_terms.frcmod \\
      --compare-with bonded_2e --note "…"
done
for b in 2e 2s; do
  $PY md/build_topology.py --model bonded --set $b --variant n2fit --start R \\
      --n2-frcmod md/qm/n2/n2_terms_new_charges_as_built.frcmod --compare-with bonded_$b --note "…"
done
```

(the notes are in `STATUS_NOTE.txt` of each system). Systems: `/root/cqc-sfa-md/systems/<system>/<start>/`
(`system.prmtop`, `system.inpcrd`, `selections.json`, `reference.pdb`, `extra_forces.json` with an empty list,
`build.json`, `STATUS_NOTE.txt`); records: `md/topology/<system>/<start>/` (`build.json`, `leap.in`,
`n2_terms.frcmod`, `SHA256SUMS`, `openmm_check.json`, `compare_with_<other>.json`).

{table_after(frags, 'System')}

No check failed in any system (`checks_failed` is empty in all eight `build.json`). The checks that do not apply
to a bonded system are counted as n.a. The three systems of `bonded_2e` have the checksums of the five system
files that the rehearsal ran with (`md/test/REHEARSAL_2e.md`): only `build.json` and `STATUS_NOTE.txt` are new
(the note and the new checks).

**Crystal water** (added on 7 Oct 2026 after the independent check of the set-up; plan, entry of 7 Oct 2026,
12:10 AM, point 2). The systems hold {_w['kept']} of the {_w['sites_in_entry']} water sites of chain A of 6RL9: site 510 (occupancy {_w510[0]['occupancy']}) is
dropped in all starts. Start R holds {_wst['R']['n_waters']}, because four ordered waters (423, 516, 553 and 568; all fully occupied,
B 11.5 to 17.6 Å² against a median of 22.5 Å² for all water) lie on the rotated ligand,
{_wr_txt} from its heavy atoms, and were removed when the start was built, by
the rule that a water whose oxygen is within 2.2 Å of a heavy atom of the ligand is removed (`md/SYSTEM.md`,
section 4; `md/system/decisions.json`, keys `water` and `starts.R.waters_removed`; `composition.crystal_water_molecules`
in the eight `build.json`: {B['bonded_2e']['X']['composition']['crystal_water_molecules']} in X and N, {B['bonded_2e']['R']['composition']['crystal_water_molecules']} in R). So the runs from X and N start with these four waters in the place
that the rotated ring would take, and the runs from R start with the ligand there and without them. The docking
test `dock/WATER4.md` is built on the same four waters.

**Each start has the same coordinates, box, ions and selections in every set**: `system.inpcrd`,
`selections.json` and `reference.pdb` have the same sha256 as in `bonded_2e` (table, and
`compare_with_*.json`, keys `system_inpcrd_same_sha256`, `selections_json_same_sha256`,
`reference_pdb_same_sha256`).

**What differs between the sets, by read-back of the two topologies** (`compare_with_<other>.json` in each
record folder; every bond, angle, torsion and improper term compared atom by atom, with the atom names, types,
charges, masses, Lennard-Jones parameters, the excluded-atom lists and the 1–4 pairs with their scaling):

| System | Compared with | Terms that differ | Any other difference | 1–4 pairs | Sections of the prmtop that differ as text |
|---|---|---|---|---|---|
""")
for m in ('bonded_2s', 'bonded_2e_n2fit', 'bonded_2s_n2fit'):
    for s in B[m]:
        c = B[m][s]['compared_with']
        dif = ', '.join(f"{v} × {k}" for k, v in sorted(c['differences_by_class'].items()))
        nother = sum(c['atoms_that_differ_in'].values()) + len(c['classes_that_differ_but_may_not']) + sum(v['in_one_only'] for v in c['terms'].values() if m == 'bonded_2s')
        A(f"| {m}/{s} | {c['other_model']} | {dif} | {'none' if c['only_the_allowed_terms_differ'] else '**YES**'} (atoms {sum(c['atoms_that_differ_in'].values())}, classes not allowed {len(c['classes_that_differ_but_may_not'])}) | {c['pairs_1_4']['this']:,} = {c['pairs_1_4']['other']:,}, identical with their scaling: {str(c['pairs_1_4']['identical_with_their_scaling']).lower()} | {', '.join(c['prmtop_sections']['that_differ'])} |\n")
A(f"""
So `bonded_2s` differs from `bonded_2e` in the 4 zinc bonds and the 12 zinc angles of section 2.2 only, and each
`n2fit` system differs from its base set only in the terms at N2: the 3 angles, the 8 torsion quadruples (2 about
N2–S: H–N2–S–C4 and Zn–N2–S–C4; 6 about Zn–N2: each histidine nitrogen with S and with H) and the 1 improper at
N2; the 1–4 pairs are identical. (In the `n2fit` systems the pointers and the 1–4 scale-factor tables differ as
text because the new torsion types add entries to them; the pairs themselves and their scaling are compared one
by one and are identical.)

The builder's checks were made **set-specific, not weaker**: `torsions_through_zinc_all_zero` holds for the sets
as built and is replaced, only in a variant, by `torsions_through_zinc_zero_or_as_the_n2_file` (a torsion through
the zinc is zero unless the file of the terms at N2 carries it, and then it is the file's);
`n2_file_is_of_the_class_of_the_set` refuses a file that holds anything but the terms of its class;
`terms_at_N2_are_those_of_the_n2_file` reads every entry of the file back from the topology;
`differs_from_the_compared_set_only_where_the_sets_differ` is the comparison above.

## 4. Gate G-Z

The plan (stage S2; entries of 5 Oct 2026, 9:41 PM, and 6 Oct 2026, 3:02 AM): "after minimisation and 100 ps",
(a) Zn–N2 within 0.10 Å of 1.93 Å, (b) each Zn–N(His) within 0.10 Å of 2.01–2.02 Å, (c) each mean N–Zn–N angle
within 10° of its value in 6RL9, (d) an oxygen (sulfonyl, water or protein) within 2.6 Å of the zinc in at most
5 % of the frames; Zn···Glu106 OE1 reported by name. These are parts a to d of GV3 of `md/analysis/SPEC.md`,
section 11, and their numbers are taken as printed there.

**Which frames.** The gate is read from the validation run X0 of each set, on **the first 100 ps after the
minimisation**: step `eq01` (heating from 50 K to 300 K, NVT), frames 0 to 9 of the equilibration series of the
analysis (10 solute frames at 10, 20, … 100 ps after the minimisation; −1.99 to −1.90 ns on the time axis of the
series) and frame 0 of the all-atom series (at 100 ps). Program: `md/production/gate_GZ.py` (sha256
`{GZ['2e']['script']['sha256'][:16]}…`), which reads `timeseries_equil_solute.csv` and `timeseries_equil_full.csv`
of the analysis; records `md/production/gate_GZ_bonded_2s.json` and `md/production/gate_GZ_bonded_2e.json`. The two
analysis programs give the same numbers on these frames ({GZ['2e']['second_program']['values_compared_on_the_gate_frames']} values compared in each run, largest
difference {GZ['2e']['second_program']['largest_difference']:g}; the two series files are identical byte for byte).

{table_after(gates, 'Part')}

| | 2s | 2e |
|---|---|---|
| d: sulfonyl O within 2.6 Å of the zinc, solute frames | {GZ['2s']['gate_GZ']['d_fifth_ligand']['solute_frames']['count']} of {GZ['2s']['gate_GZ']['d_fifth_ligand']['solute_frames']['n']} (nearest {GZ['2s']['gate_GZ']['d_fifth_ligand']['solute_frames']['nearest_sulfonyl_O_min_A']:.3f} Å) | {GZ['2e']['gate_GZ']['d_fifth_ligand']['solute_frames']['count']} of {GZ['2e']['gate_GZ']['d_fifth_ligand']['solute_frames']['n']} (nearest {GZ['2e']['gate_GZ']['d_fifth_ligand']['solute_frames']['nearest_sulfonyl_O_min_A']:.3f} Å) |
| d: any O within 2.6 Å, all-atom frame | {GZ['2s']['gate_GZ']['d_fifth_ligand']['full_frames']['count']} of {GZ['2s']['gate_GZ']['d_fifth_ligand']['full_frames']['n']} (nearest: sulfonyl {GZ['2s']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_O_of_any_kind_min_A']:.3f}, water {GZ['2s']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_water_O_min_A']:.3f}, protein {GZ['2s']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_protein_O_min_A']:.3f} Å) | {GZ['2e']['gate_GZ']['d_fifth_ligand']['full_frames']['count']} of {GZ['2e']['gate_GZ']['d_fifth_ligand']['full_frames']['n']} (nearest: sulfonyl {GZ['2e']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_O_of_any_kind_min_A']:.3f}, water {GZ['2e']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_water_O_min_A']:.3f}, protein {GZ['2e']['gate_GZ']['d_fifth_ligand']['full_frames']['nearest_protein_O_min_A']:.3f} Å) |
| d | {ok(GZ['2s']['gate_GZ']['d_fifth_ligand']['passed'])} | {ok(GZ['2e']['gate_GZ']['d_fifth_ligand']['passed'])} |
| Zn···Glu106 OE1 (reported, no criterion): mean, smallest | {GZ['2s']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['mean_A']:.3f}, {GZ['2s']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['min_A']:.3f} Å; within 2.6 Å in {GZ['2s']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['count_within']} of {GZ['2s']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['n']} | {GZ['2e']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['mean_A']:.3f}, {GZ['2e']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['min_A']:.3f} Å; within 2.6 Å in {GZ['2e']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['count_within']} of {GZ['2e']['gate_GZ']['reported_zn_glu106']['d_Zn_E106OE1_A']['n']} |
| **Gate G-Z** | **{GZ['2s']['gate_GZ']['verdict']}** | **{GZ['2e']['gate_GZ']['verdict']}** |

The columns "minimised" (the structure after the minimisation, from `state_min.xml`) and "frame at 100 ps" are
for information; the verdict is taken from the mean of the 10 frames and from part d.

**What this gate can and cannot show.** In `eq01` the heavy atoms of the solute, the zinc and the ligand
included, are held to the crystal positions (U = ½ k d², k = 1000 kJ mol⁻¹ nm⁻² = 2.39 kcal mol⁻¹ Å⁻²). That is
far weaker than a zinc bond term, but not negligible beside an angle term at the zinc, so the angles of these
frames lean towards the crystal. The same four parts for every stage of the equilibration, for information and
without a verdict (stages 4 to 6 hold no atom of the zinc site):

{table_after(gates, 'Stage')}

The six angles are in the order of the table above. For set 2s one stage (the second, 30 frames) has the
sulfonyl oxygen within 2.6 Å in 2 of 30 solute frames, more than 5 %; it is not one of the frames of the gate,
and it is the same behaviour that fails part d of gate G-V in the 5 ns (section 5.3).

## 5. Validation run X0 and gate G-V

### 5.1 The runs

| | X0 of 2s | X0 of 2e |
|---|---|---|
| Folder (machine and, pulled, Pi) | `/root/cqc-sfa-md/runs/bonded_2s/X0` | `/root/cqc-sfa-md/runs/bonded_2e/X0` |
| System, `system.prmtop` sha256 | `{MAN['2s']['inputs']['system.prmtop']['sha256']}` | `{MAN['2e']['inputs']['system.prmtop']['sha256']}` |
| Driver `md/run_md.py` sha256 | `{MAN['2s']['script']['sha256'][:16]}…` | `{MAN['2e']['script']['sha256'][:16]}…` |
| Run seed (from the run id X0) | {MAN['2s']['seed']['run_seed']} | {MAN['2e']['seed']['run_seed']} |
| Started, ended | {MAN['2s']['created']['ist']}; {MAN['2s']['last_written']['ist']} | {MAN['2e']['created']['ist']}; {MAN['2e']['last_written']['ist']} |
| Minimisation | {MAN['2s']['minimisation']['wall_s']} s, L-BFGS iterations per pass {MAN['2s']['minimisation']['lbfgs_iterations_per_pass']} | {MAN['2e']['minimisation']['wall_s']} s, {MAN['2e']['minimisation']['lbfgs_iterations_per_pass']} |
| Equilibration, 6 stages, 2 ns | {MAN['2s']['equilibration']['wall_s']} s; no extension step | {MAN['2e']['equilibration']['wall_s']} s; no extension step |
| Hydration check at the end of stage 3 | skipped, the ligand itself is at the site ({MAN['2s']['hydration_check']['counts'][0]['ligand_mean_distance_from_site_A']} Å from it on average) | skipped ({MAN['2e']['hydration_check']['counts'][0]['ligand_mean_distance_from_site_A']} Å) |
| Production, 5 ns in one segment | {MAN['2s']['production']['wall_s']} s, **{MAN['2s']['production']['ns_per_day']} ns per day**; {MAN['2s']['production']['frames_solute']} solute and {MAN['2s']['production']['frames_full']} all-atom frames | {MAN['2e']['production']['wall_s']} s, **{MAN['2e']['production']['ns_per_day']} ns per day**; {MAN['2e']['production']['frames_solute']} and {MAN['2e']['production']['frames_full']} |
| Checksums of the pulled files against the driver's records | ok (watchdog log, "driver checksums … 0 bad") | ok |

Both runs have the same seeds (the seed follows from the run id), so they start with the same velocities and
differ in the zinc terms alone. Full protocol of `md/protocol.json` (sha256 `{sha(R / 'protocol.json')[:16]}…`), on the L40S, CUDA, mixed precision.

### 5.2 The analysis and the comparison of the two programs

```bash
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare /root/cqc-sfa-md/runs/bonded_2s/X0 --run-primary
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare /root/cqc-sfa-md/runs/bonded_2e/X0 --run-primary
```

run on the Pi after the watchdog had pulled each run; `analyse_md.py` sha256 `{sha(R / 'analysis/analyse_md.py')[:16]}…`,
`recheck_md.py` sha256 `{sha(R / 'analysis/recheck_md.py')[:16]}…`, specification 1.3. Both programs analysed each run in this call.

| | X0 of 2s | X0 of 2e |
|---|---|---|
| Values compared | {ncomp['2s']:,} | {ncomp['2e']:,} |
| Disagreements beyond tolerance | {RR['2s']['n_runs_disagree']} runs; 0 values | {RR['2e']['n_runs_disagree']} runs; 0 values |
| Threshold frames; verdicts that differ through them | {RR['2s']['n_threshold_items']}; {len(RR['2s']['verdicts_that_differ_through_threshold_frames'])} | {RR['2e']['n_threshold_items']}; {len(RR['2e']['verdicts_that_differ_through_threshold_frames'])} |
| Warnings (first / second program) | {len(SM['2s'].get('warnings') or [])} / {len(SM2['2s'].get('warnings') or [])} | {len(SM['2e'].get('warnings') or [])} / {len(SM2['2e'].get('warnings') or [])} |
| Exit status | 0 | 0 |
| Gate verdict, first / second program | {g['2s']['verdict']} / {SM2['2s']['gate_GV']['verdict']} | {g['2e']['verdict']} / {SM2['2e']['gate_GV']['verdict']} |
| Report | `runs/bonded_2s/X0/analysis_recheck/recheck_report.md` | `runs/bonded_2e/X0/analysis_recheck/recheck_report.md` |

**The two programs agree on every value of both runs.**

### 5.3 Gate G-V, every number

Gate G-V as `md/analysis/SPEC.md`, section 11, version 1.3, words it (the plan's entries of 5 Oct 2026, 9:41 PM);
means over the 500 solute frames of the 5 ns as written to `summary.json` (key `gate_GV`); the verdict is the
program's own (`gate_GV.verdict`; role `validation`, not provisional). In brackets: the deviation from the
reference.

{table_after(gates, 'Criterion')}

**What parts a to c of GV3 show for set 2e** (added on 7 Oct 2026 after the independent check of the set-up;
plan, entry of 7 Oct 2026, 12:10 AM, point 3). For set 2e the equilibrium values of the four Zn–N bonds and the
six N–Zn–N angles are those of 6RL9 itself (MCPB.py step 2e; section 2.2; `md/ZINC.md`, section 8). Parts a to c of GV3 therefore show only
that the rest of the force field does not move the site away from these values (largest shifts of the means of
X0 from the equilibrium values: {_sh_bond[0]:.3f} Å, {_sh_bond[1]}, and {_sh_ang[0]:.2f}°, {_sh_ang[1]}; the brackets of the table are taken
against the references of the gate, which are rounded to 0.1°); they are not an independent reproduction of the zinc geometry.
The equilibrium value of the angle Zn–N2–S in the file of the terms at N2 is the value of 6RL9 as well ({_zns['theta0_deg']:.2f}°, section 2.5);
its mean in X0, which is no number of the gate, lies {_sh_zns:.2f}° below it.
The parts of the gate that are not put in are GV1, GV2, GV3 d and GV4.

Reported beside the gate, no criterion: Zn···Glu106 OE1 mean {g['2s']['reported_zn_glu106']['d_Zn_E106OE1_A']['mean_A']:.3f} Å, smallest {g['2s']['reported_zn_glu106']['d_Zn_E106OE1_A']['min_A']:.3f} Å, within
2.6 Å in {g['2s']['reported_zn_glu106']['d_Zn_E106OE1_A']['count_within']} of {g['2s']['reported_zn_glu106']['d_Zn_E106OE1_A']['n']} frames (2s); mean {g['2e']['reported_zn_glu106']['d_Zn_E106OE1_A']['mean_A']:.3f} Å, smallest {g['2e']['reported_zn_glu106']['d_Zn_E106OE1_A']['min_A']:.3f} Å, within 2.6 Å in {g['2e']['reported_zn_glu106']['d_Zn_E106OE1_A']['count_within']} of {g['2e']['reported_zn_glu106']['d_Zn_E106OE1_A']['n']} (2e).
Zn···Glu106 OE2: mean {g['2s']['reported_zn_glu106']['d_Zn_E106OE2_A']['mean_A']:.3f} and {g['2e']['reported_zn_glu106']['d_Zn_E106OE2_A']['mean_A']:.3f} Å, smallest {g['2s']['reported_zn_glu106']['d_Zn_E106OE2_A']['min_A']:.3f} and {g['2e']['reported_zn_glu106']['d_Zn_E106OE2_A']['min_A']:.3f} Å.
Orientation of both runs: {SM['2s']['orientation']['verdict']} (2s: {SM['2s']['orientation']['count_C_window']} of {SM['2s']['orientation']['window_frames']} frames in the crystal state, {SM['2s']['orientation']['count_R_window']} rotated) and {SM['2e']['orientation']['verdict']} (2e: {SM['2e']['orientation']['count_C_window']} of {SM['2e']['orientation']['window_frames']}, {SM['2e']['orientation']['count_R_window']} rotated).

**The failure of set 2s, in full.** Part d of GV3 has two counts, each of which may be at most 5 %: the solute
frames (every 10 ps; the sulfonyl oxygens) and the all-atom frames (every 100 ps; any oxygen).

| | 2s | 2e |
|---|---|---|
| Solute frames with a sulfonyl O within 2.6 Å of the zinc | {fl['2s']['solute_trajectory']['count']} of {fl['2s']['solute_trajectory']['n']} ({100 * fl['2s']['solute_trajectory']['fraction']:.1f} %): met | {fl['2e']['solute_trajectory']['count']} of {fl['2e']['solute_trajectory']['n']} ({100 * fl['2e']['solute_trajectory']['fraction']:.1f} %): met |
| All-atom frames with any O within 2.6 Å | **{fl['2s']['full_trajectory']['count']} of {fl['2s']['full_trajectory']['n']} ({100 * fl['2s']['full_trajectory']['fraction']:.0f} %): not met** (2 of 50 would be the most allowed) | {fl['2e']['full_trajectory']['count']} of {fl['2e']['full_trajectory']['n']}: met |
| The all-atom frames concerned: time, distance of the sulfonyl O | {'; '.join(f"{float(r['time_ns']):.1f} ns, {float(r['d_Zn_Osulfonyl_min_full_A']):.3f} Å" for r in hits_full['2s'])} | none |
| Which oxygen | the sulfonyl oxygen O2 every time (O1 never: {sum(float(r['d_Zn_O1_A']) <= 2.6 for r in sol['2s'])} frames); no water and no protein oxygen within 2.6 Å in any all-atom frame (nearest water O {min(float(r['d_Zn_Owater_min_A']) for r in full['2s']):.3f} Å, nearest protein O {min(float(r['d_Zn_Oprotein_min_A']) for r in full['2s']):.3f} Å) | nearest water O {min(float(r['d_Zn_Owater_min_A']) for r in full['2e']):.3f} Å, nearest protein O {min(float(r['d_Zn_Oprotein_min_A']) for r in full['2e']):.3f} Å |
| Nearest sulfonyl O in the 500 solute frames: mean, s.d., smallest | {st.mean(osul['2s']):.3f}, {st.stdev(osul['2s']):.3f}, {min(osul['2s']):.3f} Å | {st.mean(osul['2e']):.3f}, {st.stdev(osul['2e']):.3f}, {min(osul['2e']):.3f} Å |
| Solute frames with it within 2.7 Å; within 2.8 Å | {sum(x <= 2.7 for x in osul['2s'])}; {sum(x <= 2.8 for x in osul['2s'])} | {sum(x <= 2.7 for x in osul['2e'])}; {sum(x <= 2.8 for x in osul['2e'])} |
| Mean angle S–N2–Zn | {st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol['2s']):.2f}° | {st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol['2e']):.2f}° (start X: 117.55°) |

Read plainly: the three all-atom frames are three of the thirteen solute frames, so the two counts measure the
same event, a sulfonyl oxygen that comes within 2.6 Å of the zinc for single frames; no water and no protein
oxygen is involved, and the zinc does not gain a fifth ligand that stays. With 50 frames the limit of 5 % is
passed by the third frame. By the rule as it was fixed before the run this is a failed gate, and it is counted
as one. It is also a real difference between the sets and not only a matter of counting: with set 2s the
sulfonyl oxygen is on average {st.mean(osul['2e']) - st.mean(osul['2s']):.2f} Å nearer the zinc than with set 2e
({st.mean(osul['2s']):.3f} against {st.mean(osul['2e']):.3f} Å) and within 2.7 Å in {sum(x <= 2.7 for x in osul['2s'])} frames against {sum(x <= 2.7 for x in osul['2e'])},
and it goes with smaller mean angles N2–Zn–NE2 His96 ({g['2s']['GV3_zinc']['criteria']['ang_N2_Zn_H96NE2_deg']['mean_deg']:.2f}° against {g['2e']['GV3_zinc']['criteria']['ang_N2_Zn_H96NE2_deg']['mean_deg']:.2f}°) and S–N2–Zn
({st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol['2s']):.2f}° against {st.mean(float(r['ang_S_N2_Zn_deg']) for r in sol['2e']):.2f}°), the direction of the equilibrium angles of section 2.4.

## 6. The queue as added

Machine `cqc-md-gpu-01`, task folder `/root/cqc-sfa-md/queue/tasks.d/`; one GPU task at a time; each task waits
for the one named in `depends_on`. The files are in `md/gpu/` and have the same sha256 in the task folder.

| File | sha256 | Added (IST) | Tasks |
|---|---|---|---|
| `md/gpu/tasks_x0.json` | `{sha(R / 'gpu/tasks_x0.json')}` | 6 Oct 2026, 9:36:57 PM | `X0_2s`, then `X0_2e` (both done) |
| `md/gpu/tasks_production.json` | `{sha(R / 'gpu/tasks_production.json')}` | 6 Oct 2026, 10:14:10 PM | {', '.join('`' + t['id'] + '`' for t in tasks_p['tasks'])} |
| `md/gpu/tasks_sensitivity.json` | `{sha(R / 'gpu/tasks_sensitivity.json')}` | 6 Oct 2026, 10:14:10 PM | {', '.join('`' + t['id'] + '`' for t in tasks_s['tasks'])} |

| Order | Task | Waits for | System | Run folder on the machine | Run id | Length |
|---|---|---|---|---|---|---|
""")
k = 0
for t in tasks_p['tasks'] + tasks_s['tasks']:
    k += 1
    m = re.search(r'--system "\$SYS/([^"]+)" --out "\$RUNS/([^"]+)" --run-id (\S+) .*--prod-ns (\d+) --segment-ns (\d+)', t['command'])
    A(f"| {k} | `{t['id']}` | {', '.join(t['depends_on'])} | `systems/{m.group(1)}` | `/root/cqc-sfa-md/runs/{m.group(2)}` | {m.group(3)} | {m.group(4)} ns in segments of {m.group(5)} ns |\n")
A(f"""
Every command is that of `md/RUN.md` and `md/gpu/GPU_RUN.md`, section 5:
`run_md.py --system … --out … --run-id … --platform CUDA --stage all --prod-ns 100 --segment-ns 5 --resume`; the
test of a task is the run's `manifest.json` with `status=complete`, `production.complete=true`,
`production.completed_ns=100` and `test_only=false`. Each run has a time limit of {tasks_p['tasks'][0]['timeout_s'] // 3600} hours for one attempt
(`timeout_s`; a run takes about {100 / MAN['2e']['production']['ns_per_day'] * 24 + 0.08:.1f} hours): a run that hangs is ended and continued by the next attempt, three
attempts in all.

**The two sensitivity runs** (plan, entry of 9:15 PM, point 3; each once, start R, 100 ns, after the nine; in no
count of rule 3):

- `R_zn`: the other zinc set, here **2s**, which did not pass gate G-V (section 5.3); terms at N2 as built.
  Folder `runs/bonded_2s/R1`.
- `R_n2`: zinc set 2e as in the nine runs, with the terms at N2 of
  `md/qm/n2/n2_terms_new_charges_as_built.frcmod`, the candidate that gate G-R2 did **not accept**; it brackets
  the terms as built (plan: at the relaxed near-planar nitrogen the terms as built are 2.4 kcal mol⁻¹ too high
  and this candidate 1.9 to 2.7 too low). Folder `runs/bonded_2e_n2fit/R1`.

Both have run id R1 in a folder of their own model, so the analysis reads them as start R, replica 1, of the
models `bonded_2s` and `bonded_2e_n2fit` (the model is taken from the manifest), never as runs of `bonded_2e`;
and both have the seeds of run R1 of the primary set, so the three runs R1 start with the same velocities.
Each task writes a file `SENSITIVITY_RUN.txt` into its run folder before the run starts; its text is in
`md/gpu/tasks_sensitivity.json`. The system `bonded_2s_n2fit/R` is built and on the machine but in no task (it
would have been the system of `R_n2` had 2s been the primary set).

**Not queued and not decided here:** what follows if a sensitivity run and the nine runs disagree; that is the
plan's rule 4 and the user's to read.
""")

# ---------------------------------------------------------------- section 7
A("""
## 7. The start of production, the speed and the end
""")
if watch:
    A(f"""
Run X1 (start X, replica 1, `bonded_2e`) started on the machine at 10:14:22 PM IST on 6 Oct 2026 and was watched
from the Pi (the driver's log, copied to the Pi at {watch['log_copied_ist']}):

| Step | Driver's line |
|---|---|
""")
    for ln in eq_lines:
        m = re.match(r'\[(?:\d+ \w+ \d+, )?(.*?)\] (.*)', ln)
        if 'done in' in ln or 'potential energy after' in ln or 'hydration check' in ln:
            A(f"| {m.group(1)} | {m.group(2)} |\n")
    A(f"""
| Segment | Ended (IST) | Wall time | ns per day |
|---|---|---|---|
""")
    for s_ in segs:
        A(f"| {s_['seg']} (to {s_['to_ns']} ns) | {s_['when']} | {s_['wall_s']} s | {s_['ns_day']} |\n")
    A(f"""
The hydration check does not apply to a run from start X (the ligand itself is at the site); it applies to the
runs from start R, where the held part of the equilibration is extended by up to 4 × 500 ps if no water is at
the place of the crystal ring (`md/RUN.md`, section 1A).

**Pulled to the Pi.** {watch['pull_text']}

**The status program.** `python3 ops/pi/cqc_gpu_status.py` at {watch['status_ist']} shows the progress of every run
and an end time ({watch['queue_line']}; the three lines of the rehearsal folders are left out
here, see section 9, point 6):

```
{watch['status_excerpt']}
```

**Speed and end.** Measured: {watch['speed_text']} One run of 100 ns then takes
{watch['run_s']:,} s = {watch['run_s'] / 3600:.2f} h ({watch['prod_s']:,} s of production at {watch['ns_day']:.1f} ns per day, {watch['eq_s']} s of equilibration,
{watch['min_s']} s of minimisation and start). From the start of X1 at 10:14:22 PM IST on 6 Oct 2026:

| | Estimated end (IST) |
|---|---|
| the nine runs (R3 done) | **{watch['end_nine']}** |
| `R_zn` done | {watch['end_zn']} |
| `R_n2` done, the queue empty | **{watch['end_all']}** |

The status program's own estimate at that time: {watch['status_end']}. Not in these times: extension steps of
the equilibration of the R runs (at most 4 × 500 ps, about {watch['ext_s']} s each), a pre-emption of the spot machine
and the wait for its restart, and the time the Pi needs to pull the last files (the files of one segment
are {watch['seg_mb']} MB on the Pi, so a run of 20 segments is about {20 * watch['seg_mb'] / 1000:.1f} GB).
{watch.get('cost_text', '')}
""")
else:
    A("\n(not written yet)\n")

A(f"""
## 8. How to follow and how to resume

On the Pi:

```bash
python3 ops/pi/cqc_gpu_status.py            # machine, cost, queue, every run with ns done, speed and the end time
tail -n 3 /root/cqc-sfa-md/status/watchdog.log
ls /root/cqc-sfa-md/runs/bonded_2e/         # the pulled runs; a finished segment arrives within about 15 minutes
```

On the machine (read only; `$SSH` as in `md/gpu/GPU_RUN.md`, section 2):

```bash
$SSH "sudo python3 /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md/gpu/queue.py status"
$SSH "sudo tail -n 3 /root/cqc-sfa-md/runs/bonded_2e/X1/run_md.log"     # last finished step and its ns per day
$SSH "sudo tail -n 5 /root/cqc-sfa-md/logs/X1.log"                      # the task's own output
```

**Nothing has to be done to resume.** A machine that was stopped or pre-empted is started again by the Pi's
watchdog, the runner starts by itself and starts the interrupted task again, and the run goes on after its last
complete segment (`--resume`; at most 5 ns are computed twice). This path was tested on the real machine at
9:21 PM IST on 6 Oct 2026 with the queue empty (plan, entry of 9:22 PM). A run whose command fails is tried
three times; after that the task is `failed` and the chain behind it waits: read `/root/cqc-sfa-md/logs/<task>.log`
and the run's `run_md.log`, then `sudo python3 …/md/gpu/queue.py retry <task>` on the machine. To take the
sensitivity runs out before they start: `queue.py cancel R_zn R_n2`.

**After the runs** (on the Pi, when every folder is pulled; `md/analysis/RECHECK.md`):

```bash
nice -n 10 /root/cqc-sfa-tools/md/bin/python md/analysis/recheck_md.py --compare-campaign /root/cqc-sfa-md/runs --run-primary
```

The campaign is evaluated per zinc model (`md/analysis/SPEC.md`, section 12), and the model of a run is read
from its manifest: `bonded_2s` and `bonded_2e_n2fit` are therefore models of their own there, with one run each
(R1) and an outcome that is marked as not final. They are the sensitivity runs and belong in no count of rule 3
of `bonded_2e`. The two validation runs X0 are left out of the campaign by the programs.

The file `/root/cqc-sfa-md/status/KEEP_RUNNING` (last touched at {keep_ist}, not by this task, which
touched it last at 10:13 PM) holds back the watchdog's stops for "all done" and "idle" for 6 hours after the last
touch (`ops/pi/GPU_WATCH.md`). It was not set by this task and is not removed by it; it runs out long before the
queue is empty, so the watchdog stops the machine by itself when every task is done and pulled, unless someone
touches the file again.

## 9. Defects found and fixed, and things to know

1. **Builder: a refused build could destroy a system that existed** (the accident of `md/N2.md`, section 5.6):
   the folder was emptied before the inputs were judged. Now every refusal comes before anything is touched,
   the system is written to `<start>/build/new/` and moved into place with atomic renames only when every step
   has ended; tested with a refused build in a scratch folder, which left the system as it was.
2. **Builder: the checks did not know a second zinc set or a variant of the terms at N2.** Added `--variant`,
   `--compare-with`, the reading of `IMPROPER` entries, the set-specific checks of section 3 and the
   term-by-term comparison of two systems. `--mcpb-dir` is made absolute; `--set 2s` refuses to run without
   `--n2-frcmod`, because MCPB.py's own lines at N2 differ between 2e and 2s (section 2.2).
3. **`md/qm/after_qm.sh`: a relative report folder is taken relative to the work folder**, because the script
   changes into the work folder first. The reports of set 2s first landed in
   `/root/cqc-sfa-md/qm/gr/mcpb_2s/md/qm/results/`; they were moved to
   `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/` beside those of 2e. The script was not changed (give it an
   absolute path).
4. **Set 2s fails gate G-V** by one all-atom frame in part d of GV3 (section 5.3). Not a defect of the pipeline;
   listed here because it decides the primary set.
5. **Gate G-Z as worded can hardly fail on angles**: in the first 100 ps the zinc site is held to the crystal
   positions (section 4). The numbers of the free stages and of the 5 ns are given beside it for that reason.
6. **The status table on the Pi shows the queue as of the last pull** (up to 15 minutes old), and it lists the
   three rehearsal runs of `runs_test/` as "done, pull pending" with 0 of 5 ns, because the task file of the
   rehearsal is still in the queue and `runs_test/` is not pulled. These three lines are not production and the
   totals of the table include their 15 ns as planned. Taking `tasks_rehearsal_2e.json` out of the queue (and
   only then deleting the `runs_test` folders on the machine) removes them; a done task is tested again at every
   start of the runner, so the folders must not be deleted while the file is in the queue.

## 10. Files

| What | Where |
|---|---|
| This file | `md/PRODUCTION.md` |
| Builder | `md/build_topology.py` |
| Records of the eight systems | `md/topology/bonded_2e/`, `bonded_2s/` (X, N, R), `bonded_2e_n2fit/R/`, `bonded_2s_n2fit/R/` |
| Systems | `/root/cqc-sfa-md/systems/` on the Pi and on the machine (copied with a checksum list: `/root/cqc-sfa-md/status/push_production.sha256` on the machine, 57 files, all verified) |
| Step 2s | `/root/cqc-sfa-md/qm/gr/mcpb_2s/`, reports `/root/cqc-sfa-md/qm/gr/after_qm_report_2s/` |
| Parameter files of the two zinc sets in the repository | `md/qm/gr/charges/CA2_mcpbpy_2e.frcmod` and, added on 7 Oct 2026, `md/qm/gr/charges/CA2_mcpbpy_2s.frcmod`: copies of `CA2_mcpbpy.frcmod` of the two work folders, with the sha256 of section 2.1 (compared each time this file is assembled) |
| Blocks of the Hessian, table 2e / 2s / 2ms | `md/production/seminario_blocks.py`, `md/production/seminario_blocks.json` |
| Gate G-Z | `md/production/gate_GZ.py`, `md/production/gate_GZ_bonded_2s.json`, `md/production/gate_GZ_bonded_2e.json` |
| Validation runs and their analyses | `/root/cqc-sfa-md/runs/bonded_2s/X0/`, `/root/cqc-sfa-md/runs/bonded_2e/X0/` (`analysis/`, `analysis_recheck/`) |
| Task files | `md/gpu/tasks_x0.json`, `md/gpu/tasks_production.json`, `md/gpu/tasks_sensitivity.json` |
| How this file was assembled | `md/production/assemble/`: `frags.py` and `gates_md.py` (tables from the records, written to `frags.md` and `gates.md`), `x1_watch.py` (section 7, from `x1/run_md.log.machine`, the copy of the machine's log of X1, and `x1/status.txt`, the saved status table), `make_production_md.py` (this file; it stops if one of its own statements does not hold in the records), `make_tasks.py` (the two task files), `push/` (the log and the checksum lists of the copies to the machine) |
| Documents brought up to date | `md/RUN.md` (section 2; the path in section 4; section 7, points 4 and 11), `md/gpu/GPU_RUN.md` (sections 3 and 5); on 7 Oct 2026 `md/RUN.md` again, to the hydration check as it is run (six ring carbons, 2.5 Å), to the runs made on the GPU with their drivers and protocols (the outcome of the pre-flight test among them) and to the analysis specification as fixed (its header and state paragraph; section 1A; section 5, rows "Hydration check" and "GPU precision"; the first lines of section 6.2; section 7, points 1, 2, 6 and 10) |

This file as first written (state of {watch.get('written_ist', 'the time of writing')}) is in commit `d517798` (6 Oct 2026, 10:47 PM IST),
which also holds the files of section 10 that lie under `md/` in their state of 6 Oct 2026, all but
`md/qm/gr/charges/CA2_mcpbpy_2s.frcmod`. Added on 7 Oct 2026, after the independent check of the set-up
(`paper/audit/research/md_production_check_2026-10-07.json`; plan, entry of 7 Oct 2026, 12:10 AM IST): the file
`md/qm/gr/charges/CA2_mcpbpy_2s.frcmod` and its row in section 10, the paragraph on crystal water in section 3,
the paragraph after the table of section 5.3, the pointer to it in the row "Set 2e" of section 1, the second
half of the row "Documents brought up to date" (the changes of `md/RUN.md` are named there), the sentence on
these additions at the top, and this paragraph, which replaces the line "Nothing of this is committed".
`md/production/assemble/make_production_md.py` was changed to write and to check these additions.
""")
out = ''.join(P)
(R / 'PRODUCTION.md').write_text(out)
print('written', R / 'PRODUCTION.md', len(out.splitlines()), 'lines')
