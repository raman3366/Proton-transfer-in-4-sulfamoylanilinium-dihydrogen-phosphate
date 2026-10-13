### G-Z table
| Part | Quantity | Must lie in | 2s: minimised | 2s: mean of the 10 frames | 2s: frame at 100 ps | 2s | 2e: minimised | 2e: mean of the 10 frames | 2e: frame at 100 ps | 2e |
|---|---|---|---|---|---|---|---|---|---|---|
| a | Zn–N2 (Å) | 1.83 to 2.03 | 1.949 | 1.911 | 1.898 | met | 1.927 | 1.904 | 1.951 | met |
| b | Zn–NE2 His94 (Å) | 1.91 to 2.12 | 1.993 | 1.998 | 2.053 | met | 2.008 | 2.011 | 2.039 | met |
| b | Zn–NE2 His96 (Å) | 1.91 to 2.12 | 2.105 | 2.109 | 2.110 | met | 2.034 | 2.018 | 1.970 | met |
| b | Zn–ND1 His119 (Å) | 1.91 to 2.12 | 2.056 | 2.048 | 1.999 | met | 2.042 | 2.024 | 2.005 | met |
| c | N2–Zn–NE2 His94 (°) | 98.9 to 118.9 | 111.845 | 113.080 | 114.053 | met | 104.331 | 104.119 | 106.049 | met |
| c | N2–Zn–NE2 His96 (°) | 104.3 to 124.3 | 109.855 | 110.451 | 110.191 | met | 115.722 | 117.204 | 114.036 | met |
| c | N2–Zn–ND1 His119 (°) | 106.8 to 126.8 | 117.737 | 116.390 | 110.970 | met | 118.769 | 118.422 | 115.926 | met |
| c | NE2 His94–Zn–NE2 His96 (°) | 93.8 to 113.8 | 102.797 | 103.436 | 108.361 | met | 102.677 | 101.979 | 103.616 | met |
| c | NE2 His94–Zn–ND1 His119 (°) | 103.6 to 123.6 | 115.094 | 114.933 | 116.148 | met | 114.822 | 115.073 | 117.658 | met |
| c | NE2 His96–Zn–ND1 His119 (°) | 88.4 to 108.4 | 97.181 | 96.046 | 95.433 | met | 99.753 | 99.236 | 99.002 | met |
2s: part d: sulfonyl O within 2.6 in 0 of 10 solute frames (nearest 2.6418); any O in 0 of 1 full frames (nearest any 2.8948, water 3.5169, protein 3.9063); passed True
2s: Glu106 OE1 mean 4.04856 min 3.8744 within 0/10; OE2 mean 5.60183 min 5.4055
2s: verdict PASSED; frames [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] times [-1.99, -1.98, -1.97, -1.96, -1.95, -1.94, -1.93, -1.92, -1.91, -1.9]; full [0] [-1.9]; T end 296.29; second program True n=144 worst=0.0
2s: minimised: O1 3.9718 O2 2.8386 E106OE1 4.0435; each in window True; lbfgs [906, 3312, 6, 7]; restrained displacement {'rms': 0.2842, 'max': 1.9312}
2s: script sha 680c974ea76ac5ee written 6 Oct 2026, 10:01:51 PM IST
2e: part d: sulfonyl O within 2.6 in 0 of 10 solute frames (nearest 2.8104); any O in 0 of 1 full frames (nearest any 2.9727, water 3.6476, protein 3.9373); passed True
2e: Glu106 OE1 mean 3.92777 min 3.78 within 0/10; OE2 mean 5.5028 min 5.2949
2e: verdict PASSED; frames [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] times [-1.99, -1.98, -1.97, -1.96, -1.95, -1.94, -1.93, -1.92, -1.91, -1.9]; full [0] [-1.9]; T end 299.05; second program True n=144 worst=0.0
2e: minimised: O1 3.98 O2 2.976 E106OE1 3.9195; each in window True; lbfgs [879, 4362, 6, 6]; restrained displacement {'rms': 0.2839, 'max': 1.7878}
2e: script sha 680c974ea76ac5ee written 6 Oct 2026, 10:12:56 PM IST

### by stage (means)
| Stage | Held atoms, k | 2s: Zn–N2 | 2s: Zn–N His94/96/119 | 2s: six angles | 2s: O within 2.6 Å, solute / all-atom | 2s a–d | 2e: Zn–N2 | 2e: Zn–N His94/96/119 | 2e: six angles | 2e: O within 2.6 Å | 2e a–d |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 (heat, 100 ps) | solute_heavy, 1000 | 1.911 | 1.998 / 2.109 / 2.048 | 113.1 / 110.5 / 116.4 / 103.4 / 114.9 / 96.0 | 0 of 10 / 0 of 1 | met | 1.904 | 2.011 / 2.018 / 2.024 | 104.1 / 117.2 / 118.4 / 102.0 / 115.1 / 99.2 | 0 of 10 / 0 of 1 | met |
| 2 (density, 300 ps) | solute_heavy, 1000 | 1.926 | 1.991 / 2.103 / 2.023 | 112.6 / 110.1 / 117.2 / 102.3 / 115.8 / 95.7 | 2 of 30 / 0 of 3 | **NOT MET** (d: fifth ligand) | 1.912 | 1.998 / 2.052 / 2.002 | 105.0 / 115.2 / 118.7 / 102.8 / 115.2 / 98.8 | 0 of 30 / 0 of 3 | met |
| 3 (loosen, 300 ps) | solute_heavy, 250 | 1.962 | 2.006 / 2.078 / 2.041 | 111.0 / 110.2 / 115.8 / 103.0 / 117.2 / 97.3 | 1 of 30 / 0 of 3 | met | 1.935 | 2.023 / 2.020 / 2.040 | 104.0 / 116.6 / 116.9 / 104.1 / 115.3 / 99.3 | 0 of 30 / 0 of 3 | met |
| 4 (release_side_chains_zinc_ligand, 400 ps) | protein_backbone, 250 | 1.949 | 2.007 / 2.106 / 2.078 | 112.7 / 107.6 / 119.0 / 100.0 / 115.6 / 98.0 | 0 of 40 / 0 of 4 | met | 1.933 | 2.020 / 2.039 / 2.039 | 104.1 / 116.4 / 118.5 / 101.9 / 115.6 / 99.1 | 0 of 40 / 0 of 4 | met |
| 5 (loosen_backbone, 400 ps) | protein_backbone, 50 | 1.951 | 2.008 / 2.086 / 2.063 | 112.1 / 108.5 / 119.7 / 102.1 / 114.8 / 96.0 | 1 of 40 / 0 of 4 | met | 1.917 | 2.025 / 2.041 / 2.037 | 105.3 / 117.1 / 118.5 / 102.7 / 114.3 / 97.7 | 0 of 40 / 0 of 4 | met |
| 6 (free, 500 ps) | none, 0 | 1.969 | 2.029 / 2.106 / 2.057 | 112.2 / 108.4 / 118.5 / 102.6 / 114.9 / 96.9 | 1 of 50 / 0 of 5 | met | 1.927 | 2.016 / 2.039 / 2.040 | 104.8 / 115.9 / 118.4 / 103.0 / 114.9 / 98.7 | 0 of 50 / 0 of 5 | met |

### G-V table
| Criterion | Quantity | Must be | 2s (X0) | 2s | 2e (X0) | 2e |
|---|---|---|---|---|---|---|
| GV1 | mean heavy-atom RMSD of the ligand from the crystal pose (Å) | ≤ 1.0 | 0.712205 | met | 0.550160 | met |
| GV2 | frames with the hydrogen bond N2–H···OG1 Thr199 | ≥ 80 % | 500 of 500 (100.0 %) | met | 499 of 500 (99.8 %) | met |
| GV3 a | mean Zn–N2 (Å) | 1.83 to 2.03 | 1.952221 | met | 1.924440 | met |
| GV3 b | mean Zn–NE2 His94 (Å) | 1.91 to 2.12 | 2.017209 | met | 2.016962 | met |
| GV3 b | mean Zn–NE2 His96 (Å) | 1.91 to 2.12 | 2.096659 | met | 2.047101 | met |
| GV3 b | mean Zn–ND1 His119 (Å) | 1.91 to 2.12 | 2.054899 | met | 2.034849 | met |
| GV3 c | mean N2–Zn–NE2 His94 (°); 6RL9 108.9 | 98.9 to 118.9 | 113.239917 (+4.339917) | met | 105.617480 (-3.282520) | met |
| GV3 c | mean N2–Zn–NE2 His96 (°); 6RL9 114.3 | 104.3 to 124.3 | 107.957375 (-6.342625) | met | 115.783846 (+1.483846) | met |
| GV3 c | mean N2–Zn–ND1 His119 (°); 6RL9 116.8 | 106.8 to 126.8 | 118.493403 (+1.693403) | met | 118.314403 (+1.514403) | met |
| GV3 c | mean NE2 His94–Zn–NE2 His96 (°); 6RL9 103.8 | 93.8 to 113.8 | 101.700238 (-2.099762) | met | 102.413849 (-1.386151) | met |
| GV3 c | mean NE2 His94–Zn–ND1 His119 (°); 6RL9 113.6 | 103.6 to 123.6 | 115.046929 (+1.446929) | met | 114.710536 (+1.110536) | met |
| GV3 c | mean NE2 His96–Zn–ND1 His119 (°); 6RL9 98.4 | 88.4 to 108.4 | 96.633566 (-1.766434) | met | 98.718929 (+0.318929) | met |
| GV3 d | solute frames with a sulfonyl O within 2.6 Å of the zinc | ≤ 5 % | 13 of 500 (2.6 %) | met | 1 of 500 (0.2 %) | met |
| GV3 d | all-atom frames with any O (sulfonyl, water, protein) within 2.6 Å of the zinc | ≤ 5 %, at least one frame | 3 of 50 (6.0 %) | **NOT MET** | 0 of 50 (0.0 %) | met |
| GV3 | all parts | | | **NOT MET** | | met |
| GV4 a | mean C4–S–N2 (°), s.d. | 105.8 to 111.8 | 109.968094, 2.532340 | met | 109.562265, 2.506503 | met |
| GV4 b | mean τ_fold (°), s.d. | 112.1 to 152.1 | 126.777612, 14.613415 | met | 131.375949, 13.634488 | met |
| **Gate G-V** | `gate_GV.verdict` | | | **FAILED** | | **PASSED** |
2s: Zn-Glu106 OE1 mean 4.412004 min 3.556918 within 0 of 500; OE2 mean 6.225113 min 5.167509 within 0
2s: provisional False, applies True, role validation; second program verdict FAILED passed False; warnings [] / []
2s: recheck report keys ['what', 'spec', 'script_sha256', 'primary_script_sha256', 'runs', 'primary_calls', 'largest_differences', 'n_runs', 'n_runs_agree', 'n_runs_disagree', 'n_runs_not_compared', 'agree', 'verdicts_that_differ_through_threshold_frames', 'n_threshold_items', 'result_text']
2s: manifest: run_seed 1793718340; created 6 Oct 2026, 9:37:14 PM IST; last written 6 Oct 2026, 9:49:27 PM IST; min wall 21.64 s; equil wall 238.6 s 724.35 ns/day; prod wall 464.4 s, 930.15 ns/day; frames 500/50; hydration {"requested": true, "reference_guard": {"n_ca": 258, "ca_fit_rmsd_of_system_inpcrd_on_reference_pdb_A": 0.0, "n_ca_residue_names_that_differ": 0, "ca_residue_names_that_differ_first_5": []}, "status": "decided", "applies": false, "skipped": true, "counts": [{"step": "eq03", "time_end_ps": 700.0, "n_water_O": 0, "nearest_water_O_A": 3.103, "ligand_mean_distance_from_site_A": 0.282, "applies": false, "decision": "go on"}], "n_extensions": 0, "extension_ps": 0, "final_decision": "go on: check skipped: the ligand itself is at the site", "site_without_water_when_the_solute_was_released": false}
2s: inputs prmtop 5c1e5543950f05d4db0a707f6e739657ff506c7b1166047e50af3fdc9674b319; script 6d797f080ff75679; model bonded_2s
2s: nearest sulfonyl O in the solute frames: mean 2.8348 sd 0.1239 min 2.4644 max 3.2261; <=2.6: 13; <=2.7: 74; <=2.8: 200; O1<=2.6 0 O2<=2.6 13; mean O1 4.0021 O2 2.8348; mean S-N2-Zn 109.341
2s: all-atom frames with fifth_any: [('1.600', '2.5820', '3.7580', '3.8535'), ('4.600', '2.5954', '3.8770', '4.0599'), ('4.900', '2.5676', '3.4667', '4.0426')]; water min 3.1801; protein O min 3.7027; fifth_water 0 fifth_proteinO 0
2s: solute hits at ns ['1.250', '1.450', '1.600', '1.890', '2.050', '2.450', '2.660', '3.130', '3.260', '4.120', '4.470', '4.600', '4.900']
2s: table_row {"run": "bonded_2s/X0", "model": "bonded_2s", "start": "X", "replica": 0, "seen_ns": 5.0, "planned_ns": 5.0, "complete": true, "verdict": "STAYED", "provisional": false, "frac_start_state_last20ns": 0.996, "frac_opposite_state_last20ns": 0.0, "frac_intermediate_last20ns": 0.004, "first_passage_ns": null, "n_transitions": 0, "d1_N1_F131CZ_mean_last_half_A": 4.731975, "d2_N1_N67OD1_mean_last_half_A": 7.376453, "d_N1_N67ND2_mean_last_half_A": 6.667422, "asn67_ND2_nearer_N1_fraction": 0.964, "lig_rmsd_mean_A": 0.712205, "lig_rmsd_sd_A": 0.254105, "lig_rmsd_mean_last_half_A": 0.704676, "lig_rmsd_sd_last_half_A": 0.235031, "rule1_left_crystal": false, "hb_N2H_T199OG1_fraction": 1.0, "hb_T199NH_Osulfonyl_fraction": 0.994, "hb_T199OG1H_E106OE_fraction": 1.0, "d_Zn_N2_mean_A": 1.952221, "d_Zn_N2_sd_A": 0.051669, "d_Zn_H94NE2_mean_A": 2.017209, "d_Zn_H96NE2_mean_A": 2.096659, "d_Zn_H119ND1_mean_A"
2s: orientation {"start_state": "C", "opposite_state": "R", "midpoint_d1_A": 6.7, "midpoint_d2_A": 5.2, "window_frames": 500, "window_ns": 5.0, "window_complete": false, "count_C_window": 498, "count_R_window": 0, "count_I_window": 2, "fraction_start_state_window": 0.996, "fraction_opposite_state_window": 0.0, "fraction_intermediate_window": 0.004, "verdict": "STAYED", "provisional": false, "dwell_frames": 100, "first_passage_ns": null, "first_passage_frame": null, "n_transitions": 0, "n_visits": 1, "n_long_vis
2e: Zn-Glu106 OE1 mean 4.278796 min 3.405378 within 0 of 500; OE2 mean 6.201289 min 5.092557 within 0
2e: provisional False, applies True, role validation; second program verdict PASSED passed True; warnings [] / []
2e: recheck report keys ['what', 'spec', 'script_sha256', 'primary_script_sha256', 'runs', 'primary_calls', 'largest_differences', 'n_runs', 'n_runs_agree', 'n_runs_disagree', 'n_runs_not_compared', 'agree', 'verdicts_that_differ_through_threshold_frames', 'n_threshold_items', 'result_text']
2e: manifest: run_seed 1793718340; created 6 Oct 2026, 9:49:30 PM IST; last written 6 Oct 2026, 10:01:50 PM IST; min wall 26.2 s; equil wall 239.0 s 722.95 ns/day; prod wall 467.0 s, 925.03 ns/day; frames 500/50; hydration {"requested": true, "reference_guard": {"n_ca": 258, "ca_fit_rmsd_of_system_inpcrd_on_reference_pdb_A": 0.0, "n_ca_residue_names_that_differ": 0, "ca_residue_names_that_differ_first_5": []}, "status": "decided", "applies": false, "skipped": true, "counts": [{"step": "eq03", "time_end_ps": 700.0, "n_water_O": 0, "nearest_water_O_A": 3.186, "ligand_mean_distance_from_site_A": 0.228, "applies": false, "decision": "go on"}], "n_extensions": 0, "extension_ps": 0, "final_decision": "go on: check skipped: the ligand itself is at the site", "site_without_water_when_the_solute_was_released": false}
2e: inputs prmtop a6c57a7c20aaea633317173b691837de56301feeed6f5cf38d687e91c9085aad; script 6d797f080ff75679; model bonded_2e
2e: nearest sulfonyl O in the solute frames: mean 2.9730 sd 0.1230 min 2.5469 max 3.3297; <=2.6: 1; <=2.7: 6; <=2.8: 45; O1<=2.6 0 O2<=2.6 1; mean O1 4.0085 O2 2.9730; mean S-N2-Zn 113.991
2e: all-atom frames with fifth_any: []; water min 3.1536; protein O min 3.4054; fifth_water 0 fifth_proteinO 0
2e: solute hits at ns ['4.870']
2e: table_row {"run": "bonded_2e/X0", "model": "bonded_2e", "start": "X", "replica": 0, "seen_ns": 5.0, "planned_ns": 5.0, "complete": true, "verdict": "STAYED", "provisional": false, "frac_start_state_last20ns": 0.998, "frac_opposite_state_last20ns": 0.0, "frac_intermediate_last20ns": 0.002, "first_passage_ns": null, "n_transitions": 0, "d1_N1_F131CZ_mean_last_half_A": 4.303084, "d2_N1_N67OD1_mean_last_half_A": 7.239734, "d_N1_N67ND2_mean_last_half_A": 6.870635, "asn67_ND2_nearer_N1_fraction": 0.782, "lig_rmsd_mean_A": 0.55016, "lig_rmsd_sd_A": 0.216361, "lig_rmsd_mean_last_half_A": 0.574989, "lig_rmsd_sd_last_half_A": 0.25256, "rule1_left_crystal": false, "hb_N2H_T199OG1_fraction": 0.998, "hb_T199NH_Osulfonyl_fraction": 0.998, "hb_T199OG1H_E106OE_fraction": 1.0, "d_Zn_N2_mean_A": 1.92444, "d_Zn_N2_sd_A": 0.044397, "d_Zn_H94NE2_mean_A": 2.016962, "d_Zn_H96NE2_mean_A": 2.047101, "d_Zn_H119ND1_mean_A":
2e: orientation {"start_state": "C", "opposite_state": "R", "midpoint_d1_A": 6.7, "midpoint_d2_A": 5.2, "window_frames": 500, "window_ns": 5.0, "window_complete": false, "count_C_window": 499, "count_R_window": 0, "count_I_window": 1, "fraction_start_state_window": 0.998, "fraction_opposite_state_window": 0.0, "fraction_intermediate_window": 0.002, "verdict": "STAYED", "provisional": false, "dwell_frames": 100, "first_passage_ns": null, "first_passage_frame": null, "n_transitions": 0, "n_visits": 1, "n_long_vis
