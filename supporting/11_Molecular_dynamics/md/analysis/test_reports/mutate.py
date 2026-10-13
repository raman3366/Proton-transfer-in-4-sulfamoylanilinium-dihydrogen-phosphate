# throwaway: does the synthetic test notice a wrong analysis?  Each entry breaks one rule of the SPEC in a copy of
# analyse_md.py (one or several text replacements); the copy is run on a copy of the synthetic test and compared
# with the planted values.  M01-M14: the rules of version 1.0.  M16-M41: the rules changed or added in version 1.1.
# M42-M75: the reported quantities added in version 1.2 (sections 9A to 9E of the SPEC).
# Version 1.3 (section 9D, the water count that the driver acts on): M58-M64 follow the new code; M65 is reversed (in
# version 1.2 "the six ring carbons only" was the broken copy; now the count on seven positions is); M76-M84 are new.
# M15 (a standard deviation with N - 1) is run with the others and reported as found; its difference is small.
import json, os, shutil, subprocess, sys
from pathlib import Path
S = Path(sys.argv[1]); A = Path(sys.argv[2]); PY = sys.argv[3]
ONLY = sys.argv[4].split(',') if len(sys.argv) > 4 else None
src = (A / 'analyse_md.py').read_text()
MUT = {
 # ---- version 1.0 -------------------------------------------------------------------------------------------
 'M01 no symmetry correction of the ligand RMSD': [("    lab = np.argmin(rms, axis=0)   ", "    lab = np.zeros(rms.shape[1], int)   ")],
 'M02 no fit on the CA atoms': [("    return U @ D @ Vt, cx, cr", "    return np.broadcast_to(np.eye(3), H.shape).copy(), cx, cr")],
 'M03 no imaging of zinc and ligand': [("\n    imaging = image_solute(A, lig_all_xyz, D)", "\n    imaging = {'frames_with_box': n, 'frames_without_box': 0, 'frames_zinc_moved': 0, 'frames_ligand_moved': 0}")],
 'M04 water distances without minimum image': [("            box = ts.dimensions\n", "            box = None\n")],
 'M05 90 % rule strict (more than 90 %)': [("    if den * cw[opp] >= num * w:", "    if den * cw[opp] > num * w:")],
 'M06 hydrogen-bond angle at the donor, not at the hydrogen': [("angle(A['N2'], A['HN2'], A['thr199_OG1'])", "angle(A['HN2'], A['N2'], A['thr199_OG1'])")],
 'M07 time axis from zero': [("'t_first_ns': r6((i + 1) * dt_ns)", "'t_first_ns': r6(i * dt_ns)")],
 'M08 dwell of 0.4 ns instead of 1 ns': [("    dwell = max(1, int(round(CONST['dwell_ns'] / dt_ns)))", "    dwell = max(1, int(round(0.4 / dt_ns)))")],
 'M09 tau_fold from one torsion only': [("np.maximum(np.abs(q['tors_C3_C4_S_N2_deg']), np.abs(q['tors_C5_C4_S_N2_deg']))", "np.abs(q['tors_C3_C4_S_N2_deg'])")],
 'M10 intermediate frames not carried': [("        if s in ('C', 'R'):\n            cur = s\n", "        cur = s\n")],
 'M11 last half taken as the first half': [("    last = x[n - n // 2:]\n    blocks", "    last = x[:n // 2]\n    blocks")],
 'M12 midpoint of d1 at 7.0 instead of 6.7 A': [("    'midpoint_d1_N1_F131CZ_A': 6.7,", "    'midpoint_d1_N1_F131CZ_A': 7.0,")],
 'M13 fifth ligand at 2.3 instead of 2.6 A': [("    'fifth_ligand_distance_max_A': 2.6,", "    'fifth_ligand_distance_max_A': 2.3,")],
 'M14 solute index taken as the index of the full system': [("    pos = {p: k for k, p in enumerate(sol)}  ", "    pos = {p: p for k, p in enumerate(sol)}  ")],
 'M15 standard deviation with N - 1 in the denominator': [("        out['mean'], out['sd'] = r6(x.mean()), r6(x.std())           # population standard deviation",
                                                           "        out['mean'], out['sd'] = r6(x.mean()), r6(x.std(ddof=1))")],
 # ---- version 1.1: N-Zn-N angles ------------------------------------------------------------------------------
 'M16 N-Zn-N angles in the window 100-120 degrees (the rule of version 1.0)': [("        lo, hi = _window(ref, tol)\n        z[c] =", "        lo, hi = 100.0, 120.0\n        z[c] =")],
 'M17 N-Zn-N angles within 12 instead of 10 degrees': [("        'ang_N_Zn_N_tolerance_deg': 10.0,", "        'ang_N_Zn_N_tolerance_deg': 12.0,")],
 'M18 N-Zn-N angles within 8 instead of 10 degrees': [("        'ang_N_Zn_N_tolerance_deg': 10.0,", "        'ang_N_Zn_N_tolerance_deg': 8.0,")],
 'M19 references of two N-Zn-N angles exchanged (108.9 and 98.4)': [("            'ang_N2_Zn_H94NE2_deg': 108.9,", "            'ang_N2_Zn_H94NE2_deg': 98.4,"),
                                                                     ("            'ang_H96NE2_Zn_H119ND1_deg': 98.4},", "            'ang_H96NE2_Zn_H119ND1_deg': 108.9},")],
 # ---- version 1.1: fifth ligand -------------------------------------------------------------------------------
 'M20 fifth ligand: no frame allowed (the rule of version 1.0)': [("        'fifth_ligand_fraction_max': [1, 20],", "        'fifth_ligand_fraction_max': [0, 1],")],
 'M21 fifth ligand: up to 10 % instead of 5 % of the frames': [("        'fifth_ligand_fraction_max': [1, 20],", "        'fifth_ligand_fraction_max': [1, 10],")],
 'M22 fifth ligand, solute trajectory: fewer than 5 % (limit not included)': [("    ok_s = den * fs['count'] <= num * fs['n']", "    ok_s = den * fs['count'] < num * fs['n']")],
 'M23 fifth ligand, full trajectory: fewer than 5 % (limit not included)': [("    ok_f = bool(have_full and den * ff['count'] <= num * ff['n'])", "    ok_f = bool(have_full and den * ff['count'] < num * ff['n'])")],
 'M24 fifth ligand: the solute trajectory is not counted': [("        'passed': bool(ok_s and ok_f)}", "        'passed': bool(ok_f)}")],
 'M25 fifth ligand: the full trajectory is not counted': [("        'passed': bool(ok_s and ok_f)}", "        'passed': bool(ok_s)}")],
 'M26 fifth ligand: one fraction over the frames of both trajectories together': [
     ("    ok_s = den * fs['count'] <= num * fs['n']", "    ok_s = den * (fs['count'] + ff['count']) <= num * (fs['n'] + ff['n'])"),
     ("    ok_f = bool(have_full and den * ff['count'] <= num * ff['n'])", "    ok_f = bool(have_full and ok_s)")],
 'M27 fifth ligand, full trajectory: water and protein oxygen only, no sulfonyl oxygen': [("            cols['fifth_any'].append(int(np.nanmin([dw, dp, ds]) <= cut5))", "            cols['fifth_any'].append(int(np.nanmin([dw, dp]) <= cut5))")],
 'M28 fifth ligand: fractions over the last half of the run': [("    fs, ff = fr['fifth_sulfonylO'], fr['fifth_any']",
     "    fs, ff = ({'count': x['count_last_half'], 'n': x['n_last_half'], 'fraction': x['fraction_last_half']} for x in (fr['fifth_sulfonylO'], fr['fifth_any']))")],
 # ---- version 1.1: ligand internal coordinates -------------------------------------------------------------------
 'M29 C4-S-N2 within 4 instead of 3 degrees': [("        'ang_C4_S_N2_tolerance_deg': 3.0,", "        'ang_C4_S_N2_tolerance_deg': 4.0,")],
 'M30 C4-S-N2 within 2 instead of 3 degrees': [("        'ang_C4_S_N2_tolerance_deg': 3.0,", "        'ang_C4_S_N2_tolerance_deg': 2.0,")],
 'M31 tau_fold within 25 instead of 20 degrees': [("        'tau_fold_tolerance_deg': 20.0,", "        'tau_fold_tolerance_deg': 25.0,")],
 'M32 tau_fold within 15 instead of 20 degrees': [("        'tau_fold_tolerance_deg': 20.0,", "        'tau_fold_tolerance_deg': 15.0,")],
 'M33 C4-S-N2 and tau_fold reported only, not part of the verdict (version 1.0)': [("out['GV3_zinc']['passed']\n                         and out['GV4_ligand_internal']['passed'])", "out['GV3_zinc']['passed'])")],
 # ---- version 1.1: the validation run X0 --------------------------------------------------------------------------
 'M34 run X0 not treated as the validation run': [("    role = 'validation' if (start == vr['start'] and rep_int == vr['replica']) else 'production'", "    role = 'production'")],
 'M35 validation run: planned length from the manifest, then 100 ns': [("    if role == 'validation':                                  # fixed by the specification, not by the manifest", "    if False:")],
 'M36 run X0 entered in the campaign table': [("        if not re.fullmatch(r'[XNR][1-9]\\d*', d.name): ", "        if not re.fullmatch(r'[XNR]\\d+', d.name): "),
                                               ("        if s['run']['role'] != 'production':  ", "        if False:  ")],
 'M37 gate G-V with a verdict for every run of start X (version 1.0)': [("    gv = gate_gv(stats, fracs, R['role'], complete)", "    gv = gate_gv(stats, fracs, 'validation' if R['start'] == 'X' else 'production', complete)")],
 'M38 rule 1 applied to the validation run as well': [("R['start'] == 'X' and R['role'] == 'production',", "R['start'] == 'X',")],
 # ---- version 1.1: Asn67 ND2 --------------------------------------------------------------------------------------
 'M39 Asn67 flag reversed (set where OD1 is nearer)': [("(q['d_N1_N67ND2_A'] < d2).astype(int)", "(q['d_N1_N67ND2_A'] > d2).astype(int)")],
 'M40 distance to Asn67 ND2 taken to OD1': [("        q['d_N1_N67ND2_A'] = dist(A['N1'], A['asn67_ND2'])", "        q['d_N1_N67ND2_A'] = dist(A['N1'], A['asn67_OD1'])")],
 'M41 orientation distance d2 to the nearer of OD1 and ND2': [("    d2 = dist(A['N1'], A['asn67_OD1'])\n",
     "    d2 = np.minimum(dist(A['N1'], A['asn67_OD1']), dist(A['N1'], A['asn67_ND2'])) if 'asn67_ND2' in A else dist(A['N1'], A['asn67_OD1'])\n")],
 # ---- version 1.2, section 9A: configuration at N2 ----------------------------------------------------------------
 'M42 signed volume at N2 with HN21 and S exchanged (the sign reversed)': [("    q['n2_volume'] = (e_zn * np.cross(e_h, e_s)).sum(-1)", "    q['n2_volume'] = (e_zn * np.cross(e_s, e_h)).sum(-1)")],
 'M43 signed volume at N2 from vectors that are not normalised': [("    e_zn, e_h, e_s = (unit_vector(A[k] - A['N2']) for k in ('zn', 'HN2', 'S'))", "    e_zn, e_h, e_s = ((A[k] - A['N2']) for k in ('zn', 'HN2', 'S'))")],
 'M44 sign at N2 of the rotated start taken as +1': [("    'n2_sign_of_start': {'X': 1, 'N': 1, 'R': -1},", "    'n2_sign_of_start': {'X': 1, 'N': 1, 'R': 1},")],
 'M45 a sign change at N2 counts after 0.5 ns instead of 1 ns': [("    'n2_dwell_ns': 1.0,", "    'n2_dwell_ns': 0.5,")],
 'M46 a stretch of exactly 1 ns is not long (more than 100 frames asked)': [("                          'long': (b - a + 1) >= dwell})", "                          'long': (b - a + 1) > dwell})")],
 'M47 every change of sign at N2 counted (short stretches not left out)': [("    longs = [x for x in stretches if x['long']]\n    changes, prev", "    longs = list(stretches)\n    changes, prev")],
 'M48 time of a sign change at N2 counted from zero': [("'t_first_ns': r6((a + 1) * dt_ns)", "'t_first_ns': r6(a * dt_ns)")],
 'M49 sign at N2 counted per assigned state, not per state of the frame': [("n2_configuration(q['n2_volume'], q['state'], start_sign, dt_ns)", "n2_configuration(q['n2_volume'], q['state_assigned'], start_sign, dt_ns)")],
 'M50 sign at N2 of the last 20 ns counted over the whole run': [("    ca, cw = counts(sign), counts(sign[n - w:])", "    ca, cw = counts(sign), counts(sign)")],
 # ---- version 1.2, section 9B: the hydrogen bond around the first passage ---------------------------------------------
 'M51 hydrogen bond split at the first frame in the opposite state (not at the first long visit)': [
     ("        'first_passage_frame': None if first is None else first['first_frame'],", "        'first_passage_frame': next((v['first_frame'] for v in visits if v['state'] == opp), None),")],
 'M52 the frame of the first passage counted with the frames before it': [("    before, after = hb[:p], hb[p:]", "    before, after = hb[:p + 1], hb[p + 1:]")],
 'M53 hydrogen bond of Thr199 N-H to the sulfonyl oxygens split instead of N2-H...OG1': [("hb_split = hbond_first_passage(q['hb_N2H_T199OG1'], orient)", "hb_split = hbond_first_passage(q['hb_T199NH_Osulfonyl'], orient)")],
 # ---- version 1.2, section 9C: zinc to Glu106 -----------------------------------------------------------------------
 'M54 OE1 and OE2 of Glu106 exchanged': [("    q['d_Zn_E106OE1_A'], q['d_Zn_E106OE2_A'] = dist(zn, A['glu106_OE1']), dist(zn, A['glu106_OE2'])",
                                          "    q['d_Zn_E106OE1_A'], q['d_Zn_E106OE2_A'] = dist(zn, A['glu106_OE2']), dist(zn, A['glu106_OE1'])")],
 'M55 Glu106 counted within 2.8 instead of 2.6 A': [("    'glu106_distance_max_A': 2.6,", "    'glu106_distance_max_A': 2.8,")],
 'M56 Glu106 counted within 2.52 instead of 2.6 A': [("    'glu106_distance_max_A': 2.6,", "    'glu106_distance_max_A': 2.52,")],
 'M57 Glu106 made a criterion: counted in the solute fraction of the fifth ligand of gate G-V': [
     ("    gv = gate_gv(stats, fracs, R['role'], complete)", "    gv = gate_gv(stats, dict(fracs, fifth_sulfonylO=fraction(np.maximum(q['fifth_sulfonylO'], q['fifth_E106OE']), blk)), R['role'], complete)")],
 # ---- version 1.2, section 9D: water in the place of the crystal ring -------------------------------------------------
 'M58 crystal ring positions not carried into the frame (no fit of the full frame)': [("                site = (ref_site - crf) @ Qf[0].T + cxf[0, 0]", "                site = ref_site")],
 'M59 the fit applied the wrong way round (rotation not inverted)': [("                site = (ref_site - crf) @ Qf[0].T + cxf[0, 0]", "                site = (ref_site - crf) @ Qf[0] + cxf[0, 0]")],
 'M60 ring positions of the ligand in the frame instead of the crystal positions': [("                site = (ref_site - crf) @ Qf[0].T + cxf[0, 0]", "                site = x[[sel['lig_heavy'][a] for a in site_atoms]].astype(float)")],
 'M61 the wider count within 3.3 instead of 3.5 A': [("    'ring_wide_water_distance_max_A': 3.5,", "    'ring_wide_water_distance_max_A': 3.3,")],
 'M62 the wider count within 3.7 instead of 3.5 A': [("    'ring_wide_water_distance_max_A': 3.5,", "    'ring_wide_water_distance_max_A': 3.7,")],
 'M63 pairs of water and position counted instead of water molecules': [
     ("                    dsite = distance_array(site.astype(np.float32), x[wat], box=box).min(axis=0)\n                    cols[col_n].append(int(np.sum(dsite <= cut_site)))",
      "                    dall = distance_array(site.astype(np.float32), x[wat], box=box)\n                    dsite = dall.min(axis=0)\n                    cols[col_n].append(int(np.sum(dall <= cut_site)))")],
 'M64 water at the ring positions without the minimum image': [("x[wat], box=box).min(axis=0)", "x[wat], box=None).min(axis=0)")],
 # version 1.3: reversed.  The count that the driver acts on is taken on the six ring carbons; seven positions is the broken copy
 'M65 the count of the driver on seven positions (the aniline N added, as version 1.2 counted)': [
     ("    'ring_place_atoms_default': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],", "    'ring_place_atoms_default': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'],")],
 # ---- version 1.2, section 9E: equilibration frames and the first production frame -------------------------------------
 'M66 equilibration time counted from the start of the equilibration (positive times)': [("(s0 + (k + 1) * interval_ps - t_eq) / 1000.0", "(s0 + (k + 1) * interval_ps) / 1000.0")],
 'M67 first frame of a stage put at the start of the stage': [("(s0 + (k + 1) * interval_ps - t_eq) / 1000.0", "(s0 + k * interval_ps - t_eq) / 1000.0")],
 'M68 stage table of the protocol used although the manifest lists the stages': [("    recs = _dig(man, 'equilibration', 'stages')\n    if isinstance(recs, list):\n        tab = {}", "    recs = None\n    if isinstance(recs, list):\n        tab = {}")],
 'M69 equilibration frames read as frames of the production': [
     ("    for p in sorted(Path(run_dir).glob(f'{kind}_seg*.dcd')):\n        m = re.fullmatch(rf'{kind}_seg(\\d+)\\.dcd', p.name)",
      "    for p in sorted(Path(run_dir).glob(f'{kind}_*.dcd')):\n        m = re.fullmatch(rf'{kind}_(?:eq|seg)(\\d+)\\.dcd', p.name)")],
 'M70 first production frame said to differ only when it is in the opposite state': [("        'state_differs': bool(q['state'][0] != start_state),", "        'state_differs': bool(q['state'][0] == orient['opposite_state']),")],
 'M71 sign at N2 of the first production frame compared with +1 for every start': [("        'n2_sign_differs': bool(q['n2_sign'][0] != start_sign),", "        'n2_sign_differs': bool(q['n2_sign'][0] != 1),")],
 'M72 no imaging of zinc and ligand in the equilibration frames': [("        imaging = image_solute(A, lig_all_xyz, D)\n        q, _ = solute_quantities", "        imaging = {'frames_with_box': ns, 'frames_without_box': 0, 'frames_zinc_moved': 0, 'frames_ligand_moved': 0}\n        q, _ = solute_quantities")],
 'M73 extension steps of the held part left out (only the files eqNN read)': [("EQ_STEP = r'eq(\\d+)(?:x([1-9]\\d*))?'", "EQ_STEP = r'eq(\\d+)(x[1-9]\\d*){0}'")],
 'M74 extension steps read after all stages, not after their own stage': [("    return 1000 * int(nn) + int(mm or 0)", "    return 1000 * int(mm or 0) + int(nn)")],
 'M75 the count of the driver compared with the first full frame of the step, not the last': [
     ("        k = np.where(fstage == nn)[0]\n        if len(k):\n            i = int(k[-1])", "        k = np.where(fstage == nn)[0]\n        if len(k):\n            i = int(k[0])")],
 # ---- version 1.3, section 9D: the count that the driver acts on, its settings, the comparison with the driver ----------
 'M76 the wider count without the aniline N (six positions instead of seven)': [
     ("    'ring_wide_atoms': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'],", "    'ring_wide_atoms': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],")],
 'M77 the count of the driver within 2.7 instead of 2.5 A (default cut-off)': [("    'ring_place_water_cutoff_default_A': 2.5,", "    'ring_place_water_cutoff_default_A': 2.7,")],
 'M78 the count of the driver within 2.3 instead of 2.5 A (default cut-off)': [("    'ring_place_water_cutoff_default_A': 2.5,", "    'ring_place_water_cutoff_default_A': 2.3,")],
 'M79 site atoms and cut-off of the manifest not read (the defaults for every run)': [
     ("    st = _dig(man, 'hydration_check', 'settings')\n    if st is None:\n        return out", "    st = _dig(man, 'hydration_check', 'settings')\n    if True:\n        return out")],
 'M80 the cut-off of the manifest read, its site atoms not': [
     ("    out.update(site_atoms=[str(a) for a in atoms], water_cutoff_A=float(cut),", "    out.update(site_atoms=atoms0, water_cutoff_A=float(cut),")],
 'M81 the site atoms of the manifest read, its cut-off not': [
     ("    out.update(site_atoms=[str(a) for a in atoms], water_cutoff_A=float(cut),", "    out.update(site_atoms=[str(a) for a in atoms], water_cutoff_A=cut0,")],
 'M82 the count of the driver compared with the wider count (the comparison of version 1.2)': [
     ("count_agrees=bool(lf['n_water_ring_place'] == d['n_water_O'])", "count_agrees=bool(lf['n_water_ring_wide_3p5A'] == d['n_water_O'])")],
 'M83 the nearest distance of the driver compared with that of the wider count': [
     ("abs(lf['d_water_ring_place_min_A'] - float(d['nearest_water_O_A'])) <= 0.002", "abs(lf['d_water_ring_wide_min_A'] - float(d['nearest_water_O_A'])) <= 0.002")],
 'M84 no warning when the settings of a run are not those of the plan': [
     ("    if not same:\n        warn(f'the hydration check of this run", "    if False:\n        warn(f'the hydration check of this run")],
}
for name, pairs in MUT.items():                               # every replacement must hit exactly one place
    for old, new in pairs:
        assert src.count(old) == 1, (name, src.count(old), old)
if ONLY and ONLY[0] == 'check':
    print(len(MUT), 'mutations, every replacement found once'); sys.exit(0)

T = S / 'test_mut'
env = dict(os.environ, OMP_NUM_THREADS='1')
mdir = S / 'mutation'
RUNNER = '''import sys
sys.path.insert(0, %r)
import analyse_md_mutant as m
T = %r
opts = ['--no-figures', '--no-hash-trajectories']
import glob
codes = [m.main(['--campaign', T + '/runs'] + opts), m.main([T + '/runs/test/X0'] + opts), m.main([T + '/runs_overlap/test/X9'] + opts)]
codes += [m.main([d] + opts) for d in sorted(glob.glob(T + '/runs_gv/*/*/')) + sorted(glob.glob(T + '/runs_v12/*/*/')) + sorted(glob.glob(T + '/runs_v13/*/*/'))]
print('exit codes', codes)
sys.exit(max(codes))
''' % (str(mdir), str(T))
(mdir / 'run_mutant.py').write_text(RUNNER)
out_file = mdir / 'mutation_results.json'
res = json.loads(out_file.read_text()) if (ONLY and out_file.is_file()) else {}
for name, pairs in MUT.items():
    if ONLY and name.split()[0] not in ONLY:
        continue
    text = src
    for old, new in pairs:
        text = text.replace(old, new)
    (mdir / 'analyse_md_mutant.py').write_text(text)
    shutil.rmtree(mdir / '__pycache__', ignore_errors=True)
    for d in T.glob('runs*/**/analysis*'):
        shutil.rmtree(d, ignore_errors=True)
    for f in (T / 'test_report.json', T / 'test_report.md'):
        f.unlink(missing_ok=True)
    # no second "nice" here: this script is started with nice -n 10 and its children inherit that (with a second one
    # they ran at nice 19 and, on the loaded machine, five times slower: the first start of 6 Oct 2026, stopped after M01)
    r1 = subprocess.run([PY, str(mdir / 'run_mutant.py')], capture_output=True, text=True, env=env)
    r3 = subprocess.run([PY, str(A / 'make_test_trajectory.py'), '--compare', str(T)], capture_output=True, text=True, env=env)
    rec = {'replacements': len(pairs), 'analysis_exit': r1.returncode, 'compare_exit': r3.returncode}
    if (T / 'test_report.json').is_file():
        d = json.loads((T / 'test_report.json').read_text())
        bad = [v for v in d['values'] if not v['passed']]
        rec.update(n_values=d['n_values'], n_failed=len(bad), runs_with_a_failure=sorted({v['run'] for v in bad}),
                   examples=[f"{v['run']}: {v['quantity']}: planted {v['planted']}, returned {v['returned']}" for v in bad[:3]])
        if name.startswith('M15'):
            sd = [v for v in d['values'] if v['quantity'].endswith(': sd') and v['abs_diff'] is not None]
            w = max(sd, key=lambda v: v['abs_diff'])
            rec.update(largest_difference_of_a_standard_deviation=w['abs_diff'], where=w['run'] + ' ' + w['quantity'],
                       n_frames_there=next(x['planted'] for x in d['values'] if x['run'] == w['run'] and x['quantity'] == 'solute frames seen'),
                       all_failures_are_standard_deviations=all(v['quantity'].endswith((': sd', ': sd_last_half')) for v in bad),
                       failures=[{'run': v['run'], 'quantity': v['quantity'], 'unit': v['unit'], 'planted': v['planted'], 'returned': v['returned'],
                                  'frames': next(x['planted'] for x in d['values'] if x['run'] == v['run'] and x['quantity'] == 'solute frames seen')} for v in bad])
    else:
        rec.update(n_failed=None, note=(r1.stderr.strip().splitlines() or r3.stderr.strip().splitlines() or ['?'])[-1][:300])
    res[name] = rec
    print(name, '->', json.dumps(rec)[:900], flush=True)
    out_file.write_text(json.dumps(res, indent=1) + '\n')
(mdir / 'analyse_md_mutant.py').unlink(missing_ok=True)
shutil.rmtree(mdir / '__pycache__', ignore_errors=True)
detected = sum(1 for r in res.values() if r['n_failed'] is None or r['n_failed'] > 0)
print(f'{detected} of {len(res)} mutations noticed')
