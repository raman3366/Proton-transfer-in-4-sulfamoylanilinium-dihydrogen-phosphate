# Production docking: results

Generated 2026-09-28 14:14 IST by `python dock/analyse.py --production (dock/production_analysis.py)` and `python dock/production_report.py`. Every number in this file is read from `dock/results/production_summary.json` and `dock/results/production_contacts.json`; the validation numbers quoted for context come from `dock/analysis/{summary,gate_G2}.json` and are discussed in `dock/VALIDATION.md`.

**What ran.** 72 production jobs (tag suffix `__prod`): 4 targets × 8–9 species × {Vina, Vinardo}, plus AutoDock4Zn for the two hCA II targets with the species whose sulfonamide is neutral. Exhaustiveness 32, 20 Å box, 10/20 seeds per job, `--cpu 3`, `nice -n 10`, one docking process tree at a time on the Raspberry Pi 5. Total docking time 4.6 h.

Six further AD4Zn runs (0.6 h) are the input-geometry diagnostic of section 7b; they are kept out of the tables, the rankings and the DFT pose set, and are listed in `production_summary.json` → `diagnostics`.

The three sulfanilamide species are docked **from their DFT geometries** (B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM water): `molecules/sfa/results/sulfanilamide_pcm_d3/opt.xyz`, `molecules/sfa/results/sulfonamidate_pcm_d3/opt.xyz`, `molecules/sfa/results/cation_pcm_d3/opt.xyz`. Vina and AD4Zn hold every bond length, bond angle and rigid-fragment geometry of their input fixed and randomise position, orientation and the active torsions from `--seed`, so the DFT structure is what is docked and the seed still drives the search. The controls (sulfamethoxazole, acetazolamide, p-aminobenzoate, dihydropteroate) keep their ETKDG/MMFF94s geometries: they were never optimised at the DFT level, and using the same route for them as in validation keeps those runs comparable. The ion pair is never docked.

**Scores are docking scores** (kcal/mol by the programs' convention), reported as mean ± SD over seeds with the best value; they are never converted to a Ki or a Kd, and only the ordering is interpreted.

**What was changed in the code to run this, and why.** Each change is commented at the point of change with the date.

1. `prep_ligands.py`, `run_vina.py`, `run_ad4zn.py`: **the DFT geometries were not reaching the docking.** `prep_ligands.py --dft` wrote them to `dock/ligands/`, but both runners rebuilt an ETKDG conformer from the species SMILES for every seed and never read that file, so the production run as scripted would have docked ETKDG geometries and reported them as DFT ones. A new `prep_ligands.input_mol_for_seed()` and a `--geometry {etkdg,dft}` flag on both runners fix it; the geometry actually used is recorded in each run's `results.json` (`geometry`, `input_geometry` with the xyz path and its sha256).
2. `run_production.sh`: production jobs are tagged `__prod`. As scripted they carried no suffix, so they would have overwritten the validation runs that gate G2 rests on.
3. `run_production.sh`: AD4Zn is run only for the species whose sulfonamide is neutral (sulfanilamide, sulfamethoxazole, acetazolamide). Its TZ pseudo-atom models the Zn–N(H) bond of the neutral form; in validation the anion gave Zn–N 4.86 Å for acetazolamide against 2.08 Å for the neutral form. 20 seeds for the two sulfanilamide AD4Zn jobs whose poses are discussed, 10 elsewhere.
4. `run_production.sh`: each target's own co-crystal ligand is docked as the reference score (7S2K's dihydropteroate was added for this; the others were already species of the matrix), jobs run in priority order one process tree at a time, a completed job is skipped on a rerun, and the progress log with an ETA is written to `dock/results/production.log`.
5. `analyse.py --production` (new `production_analysis.py`) and `contacts.py --production` write to `dock/results/`; the validation files `dock/analysis/{summary,gate_G2}.json` are left alone, and production runs are excluded from them. Re-running `analyse.py` after these edits reproduced both files byte for byte.

## 1. Headline numbers

| target | protocol | validated? (dock/VALIDATION.md) | sulfanilamide⁰ score | rank of sulfanilamide⁰ | pose vs the crystal ligand |
|---|---|---|---|---|---|
| 7S2K (Sul2 (sulfonamide-resistant DHPS)) | Vina | validated with Vina and Vinardo (78H redock 0.51 / 0.54 Å, 10/10 seeds within 1 Å) | -6.409 ± 0.018 (best -6.436) | 8/9 by score | 1.04 Å over the aniline-X(O)2 core, X = S or C vs 78H |
| 7S2K (Sul2 (sulfonamide-resistant DHPS)) | Vinardo | validated with Vina and Vinardo (78H redock 0.51 / 0.54 Å, 10/10 seeds within 1 Å) | -4.991 ± 0.003 (best -4.996) | 9/9 by score | 0.91 Å over the aniline-X(O)2 core, X = S or C vs 78H |
| 3HS4 (hCA II) | Vina | validated cross-dock for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10); the AZM pose itself is not reproduced | -5.767 ± 0.008 (best -5.778) | 5/8 by score | 1.42 ± 0.26 Å to SAN (f≤1 Å 0.0) |
| 3HS4 (hCA II) | Vinardo | validated cross-dock for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10); the AZM pose itself is not reproduced | -4.987 ± 0.004 (best -4.993) | 6/8 by score | 1.63 ± 0.00 Å to SAN (f≤1 Å 0.0) |
| 3HS4 (hCA II) | AD4Zn | validated cross-dock for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10); the AZM pose itself is not reproduced | -8.817 ± 0.007 (best -8.830) | 3/3 by score | 3.22 ± 0.02 Å to SAN (f≤1 Å 0.0) |
| 6RL9 (hCA II) | Vina | partially validated with AD4Zn (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal: fails the 1.0 Å gate) | -5.880 ± 0.012 (best -5.897) | 5/8 by score | 1.08 ± 0.00 Å to SAN (f≤1 Å 0.0) |
| 6RL9 (hCA II) | Vinardo | partially validated with AD4Zn (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal: fails the 1.0 Å gate) | -5.077 ± 0.004 (best -5.083) | 5/8 by score | 1.11 ± 0.00 Å to SAN (f≤1 Å 0.0) |
| 6RL9 (hCA II) | AD4Zn | partially validated with AD4Zn (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal: fails the 1.0 Å gate) | -8.836 ± 0.004 (best -8.843) | 3/3 by score | 0.61 ± 0.04 Å to SAN (f≤1 Å 1.0) |
| 3TZF (DHPS, Y. pestis (chain B)) | Vina | passes the 2.0 Å criterion only (SMX anion, Vina, 1.92 ± 0.04 Å); poses qualitative | -6.423 ± 0.009 (best -6.438) | 3/8 by score | 0.20 Å over the aniline-SO2N core vs 08D |
| 3TZF (DHPS, Y. pestis (chain B)) | Vinardo | passes the 2.0 Å criterion only (SMX anion, Vina, 1.92 ± 0.04 Å); poses qualitative | -5.310 ± 0.019 (best -5.352) | 4/8 by score | 0.96 Å over the aniline-SO2N core vs 08D |

The rank column is the rank by mean score among all species docked to that target with that scoring function, **with no site filter**: for hCA II it can rank above sulfanilamide a species whose top pose does not reach the zinc (sulfamethoxazole with AD4Zn is exactly that case — Zn–N 3.28–3.56 Å, the zinc-site check passing in 0.00 of the seeds in both hCA II targets). The ranking restricted to zinc-bound top poses is in section 2, and section 4 has the distances.

## 2. Scores, clusters and ligand efficiency, per target × scoring × species

`score` = top pose of each seed, mean ± SD (SD over seeds, ddof 1). `clusters` = populations of the leader clusters of those top poses (symmetry-corrected heavy-atom RMSD ≤ 2.0 Å, seeds ordered by score), largest cluster first as `population/seeds`. `LE` = −score/heavy atoms; `FQ` = LE / LE_scale(HA) with Reynolds' size correction — read section 5 before using either. `Δ ref` = score minus the score of the receptor's own co-crystal ligand in the same scoring function.

### 7S2K — Sul2 (sulfonamide-resistant DHPS)

Co-crystal ligand 78H (7,8-dihydropteroate). Validation: validated with Vina and Vinardo (78H redock 0.51 / 0.54 Å, 10/10 seeds within 1 Å).

| species | scoring | score (kcal/mol) | best | seeds | clusters | LE | FQ | Δ ref |
|---|---|---|---|---|---|---|---|---|
| sulfanilamide (neutral) | Vina | -6.409 ± 0.018 | -6.436 | 10 | 10/10 | 0.583 | 0.84 | 3.70 |
| sulfanilamide (neutral) | Vinardo | -4.991 ± 0.003 | -4.996 | 10 | 10/10 | 0.454 | 0.65 | 3.10 |
| sulfonamidate anion | Vina | -6.678 ± 0.015 | -6.700 | 10 | 10/10 | 0.607 | 0.87 | 3.43 |
| sulfonamidate anion | Vinardo | -5.268 ± 0.001 | -5.269 | 10 | 10/10 | 0.479 | 0.69 | 2.82 |
| 4-sulfamoylanilinium cation | Vina | -6.376 ± 0.008 | -6.390 | 10 | 10/10 | 0.580 | 0.83 | 3.73 |
| 4-sulfamoylanilinium cation | Vinardo | -5.026 ± 0.004 | -5.034 | 10 | 10/10 | 0.457 | 0.66 | 3.06 |
| 7,8-dihydropteroate (co-crystal ligand) | Vina | -10.106 ± 0.036 | -10.202 | 10 | 10/10 | 0.439 | 1.05 | ref |
| 7,8-dihydropteroate (co-crystal ligand) | Vinardo | -8.090 ± 0.077 | -8.186 | 10 | 10/10 | 0.352 | 0.84 | ref |
| sulfamethoxazole anion | Vina | -7.853 ± 0.099 | -8.049 | 10 | 10/10 | 0.462 | 0.87 | 2.25 |
| sulfamethoxazole anion | Vinardo | -6.311 ± 0.094 | -6.489 | 10 | 10/10 | 0.371 | 0.70 | 1.78 |
| sulfamethoxazole (neutral) | Vina | -8.170 ± 0.035 | -8.224 | 10 | 10/10 | 0.481 | 0.91 | 1.94 |
| sulfamethoxazole (neutral) | Vinardo | -6.659 ± 0.043 | -6.698 | 10 | 10/10 | 0.392 | 0.74 | 1.43 |
| acetazolamide (neutral) | Vina | -6.971 ± 0.089 | -7.046 | 10 | 10/10 | 0.536 | 0.84 | 3.13 |
| acetazolamide (neutral) | Vinardo | -5.251 ± 0.035 | -5.294 | 10 | 7/10, 3/10 | 0.404 | 0.63 | 2.84 |
| acetazolamide anion | Vina | -6.826 ± 0.083 | -6.931 | 10 | 10/10 | 0.525 | 0.82 | 3.28 |
| acetazolamide anion | Vinardo | -5.154 ± 0.091 | -5.331 | 10 | 3/10, 7/10 | 0.396 | 0.62 | 2.94 |
| p-aminobenzoate | Vina | -7.501 ± 0.009 | -7.516 | 10 | 10/10 | 0.750 | 1.04 | 2.60 |
| p-aminobenzoate | Vinardo | -6.422 ± 0.003 | -6.427 | 10 | 10/10 | 0.642 | 0.89 | 1.67 |

**Rank order, Vina** (best first): dhp_anion (-10.11) < smx_neutral (-8.17) < smx_anion (-7.85) < paba_anion (-7.50) < azm_neutral (-6.97) < azm_anion (-6.83) < sfa_anion (-6.68) < sfa_neutral (-6.41) < sfa_cation (-6.38). Reference ligand: `dhp_anion` at -10.11.

**Rank order, Vinardo** (best first): dhp_anion (-8.09) < smx_neutral (-6.66) < paba_anion (-6.42) < smx_anion (-6.31) < sfa_anion (-5.27) < azm_neutral (-5.25) < azm_anion (-5.15) < sfa_cation (-5.03) < sfa_neutral (-4.99). Reference ligand: `dhp_anion` at -8.09.

### 3HS4 — hCA II

Co-crystal ligand AZM (acetazolamide). Validation: validated cross-dock for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10); the AZM pose itself is not reproduced.

| species | scoring | score (kcal/mol) | best | seeds | clusters | LE | FQ | Δ ref |
|---|---|---|---|---|---|---|---|---|
| sulfanilamide (neutral) | Vina | -5.767 ± 0.008 | -5.778 | 10 | 10/10 | 0.524 | 0.75 | 0.71 |
| sulfanilamide (neutral) | Vinardo | -4.987 ± 0.004 | -4.993 | 10 | 10/10 | 0.453 | 0.65 | 0.28 |
| sulfanilamide (neutral) | AD4Zn | -8.817 ± 0.007 | -8.830 | 20 | 20/20 | 0.802 | 1.15 | 0.44 |
| sulfonamidate anion | Vina | -5.758 ± 0.006 | -5.766 | 10 | 10/10 | 0.523 | 0.75 | 0.72 |
| sulfonamidate anion | Vinardo | -5.024 ± 0.004 | -5.029 | 10 | 10/10 | 0.457 | 0.66 | 0.25 |
| 4-sulfamoylanilinium cation | Vina | -5.763 ± 0.009 | -5.774 | 10 | 10/10 | 0.524 | 0.75 | 0.71 |
| 4-sulfamoylanilinium cation | Vinardo | -4.959 ± 0.007 | -4.968 | 10 | 10/10 | 0.451 | 0.65 | 0.31 |
| sulfamethoxazole anion | Vina | -7.721 ± 0.029 | -7.741 | 10 | 10/10 | 0.454 | 0.86 | -1.24 |
| sulfamethoxazole anion | Vinardo | -5.543 ± 0.086 | -5.639 | 10 | 4/10, 2/10, 1/10 (+1 more) | 0.326 | 0.61 | -0.27 |
| sulfamethoxazole (neutral) | Vina | -7.432 ± 0.019 | -7.457 | 10 | 10/10 | 0.437 | 0.82 | -0.95 |
| sulfamethoxazole (neutral) | Vinardo | -5.589 ± 0.060 | -5.747 | 10 | 1/10, 6/10, 3/10 | 0.329 | 0.62 | -0.32 |
| sulfamethoxazole (neutral) | AD4Zn | -9.668 ± 0.038 | -9.712 | 10 | 10/10 | 0.569 | 1.07 | -0.41 |
| acetazolamide (neutral) | Vina | -6.409 ± 0.141 | -6.574 | 10 | 10/10 | 0.493 | 0.77 | 0.07 |
| acetazolamide (neutral) | Vinardo | -5.272 ± 0.032 | -5.322 | 10 | 10/10 | 0.406 | 0.64 | ref |
| acetazolamide (neutral) | AD4Zn | -9.257 ± 0.018 | -9.296 | 10 | 8/10, 2/10 | 0.712 | 1.11 | ref |
| acetazolamide anion | Vina | -6.477 ± 0.017 | -6.502 | 10 | 10/10 | 0.498 | 0.78 | ref |
| acetazolamide anion | Vinardo | -5.112 ± 0.020 | -5.147 | 10 | 10/10 | 0.393 | 0.62 | 0.16 |
| p-aminobenzoate | Vina | -5.418 ± 0.013 | -5.442 | 10 | 10/10 | 0.542 | 0.75 | 1.06 |
| p-aminobenzoate | Vinardo | -4.789 ± 0.003 | -4.794 | 10 | 10/10 | 0.479 | 0.67 | 0.48 |

**Rank order, Vina** (best first): smx_anion (-7.72) < smx_neutral (-7.43) < azm_anion (-6.48) < azm_neutral (-6.41) < sfa_neutral (-5.77) < sfa_cation (-5.76) < sfa_anion (-5.76) < paba_anion (-5.42). Reference ligand: `azm_anion` at -6.48.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): **no species has a zinc-bound top pose in any seed**, so this scoring function ranks these species on poses that do not reproduce the zinc anchor.

**Rank order, Vinardo** (best first): smx_neutral (-5.59) < smx_anion (-5.54) < azm_neutral (-5.27) < azm_anion (-5.11) < sfa_anion (-5.02) < sfa_neutral (-4.99) < sfa_cation (-4.96) < paba_anion (-4.79). Reference ligand: `azm_neutral` at -5.27.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): sfa_anion -5.02 (10/10 seeds), sfa_neutral -4.99 (10/10 seeds), sfa_cation -4.96 (9/10 seeds). Species with no zinc-bound seed are absent from that list.

**Rank order, AD4Zn** (best first): smx_neutral (-9.67) < azm_neutral (-9.26) < sfa_neutral (-8.82). Reference ligand: `azm_neutral` at -9.26.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): azm_neutral -9.26 (10/10 seeds), sfa_neutral -8.81 (6/20 seeds). Species with no zinc-bound seed are absent from that list.

### 6RL9 — hCA II

Co-crystal ligand SAN (sulfanilamide). Validation: partially validated with AD4Zn (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal: fails the 1.0 Å gate).

| species | scoring | score (kcal/mol) | best | seeds | clusters | LE | FQ | Δ ref |
|---|---|---|---|---|---|---|---|---|
| sulfanilamide (neutral) | Vina | -5.880 ± 0.012 | -5.897 | 10 | 10/10 | 0.535 | 0.77 | ref |
| sulfanilamide (neutral) | Vinardo | -5.077 ± 0.004 | -5.083 | 10 | 10/10 | 0.462 | 0.66 | ref |
| sulfanilamide (neutral) | AD4Zn | -8.836 ± 0.004 | -8.843 | 20 | 20/20 | 0.803 | 1.15 | ref |
| sulfonamidate anion | Vina | -5.838 ± 0.013 | -5.855 | 10 | 10/10 | 0.531 | 0.76 | 0.04 |
| sulfonamidate anion | Vinardo | -5.030 ± 0.003 | -5.033 | 10 | 10/10 | 0.457 | 0.66 | 0.05 |
| 4-sulfamoylanilinium cation | Vina | -5.873 ± 0.007 | -5.878 | 10 | 10/10 | 0.534 | 0.77 | 0.01 |
| 4-sulfamoylanilinium cation | Vinardo | -5.059 ± 0.003 | -5.066 | 10 | 10/10 | 0.460 | 0.66 | 0.02 |
| sulfamethoxazole anion | Vina | -7.728 ± 0.012 | -7.753 | 10 | 10/10 | 0.455 | 0.86 | -1.85 |
| sulfamethoxazole anion | Vinardo | -6.015 ± 0.022 | -6.044 | 10 | 10/10 | 0.354 | 0.67 | -0.94 |
| sulfamethoxazole (neutral) | Vina | -7.529 ± 0.149 | -7.764 | 10 | 3/10, 7/10 | 0.443 | 0.84 | -1.65 |
| sulfamethoxazole (neutral) | Vinardo | -6.007 ± 0.061 | -6.097 | 10 | 4/10, 6/10 | 0.353 | 0.67 | -0.93 |
| sulfamethoxazole (neutral) | AD4Zn | -9.847 ± 0.062 | -9.925 | 10 | 10/10 | 0.579 | 1.09 | -1.01 |
| acetazolamide (neutral) | Vina | -6.653 ± 0.025 | -6.686 | 10 | 10/10 | 0.512 | 0.80 | -0.77 |
| acetazolamide (neutral) | Vinardo | -5.235 ± 0.048 | -5.302 | 10 | 3/10, 7/10 | 0.403 | 0.63 | -0.16 |
| acetazolamide (neutral) | AD4Zn | -9.484 ± 0.026 | -9.519 | 10 | 2/10, 3/10, 5/10 | 0.730 | 1.14 | -0.65 |
| acetazolamide anion | Vina | -6.637 ± 0.016 | -6.664 | 10 | 10/10 | 0.511 | 0.80 | -0.76 |
| acetazolamide anion | Vinardo | -5.133 ± 0.062 | -5.259 | 10 | 2/10, 7/10, 1/10 | 0.395 | 0.62 | -0.06 |
| p-aminobenzoate | Vina | -5.449 ± 0.005 | -5.456 | 10 | 8/10, 2/10 | 0.545 | 0.76 | 0.43 |
| p-aminobenzoate | Vinardo | -4.978 ± 0.002 | -4.980 | 10 | 10/10 | 0.498 | 0.69 | 0.10 |

**Rank order, Vina** (best first): smx_anion (-7.73) < smx_neutral (-7.53) < azm_neutral (-6.65) < azm_anion (-6.64) < sfa_neutral (-5.88) < sfa_cation (-5.87) < sfa_anion (-5.84) < paba_anion (-5.45). Reference ligand: `sfa_neutral` at -5.88.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): **no species has a zinc-bound top pose in any seed**, so this scoring function ranks these species on poses that do not reproduce the zinc anchor.

**Rank order, Vinardo** (best first): smx_anion (-6.01) < smx_neutral (-6.01) < azm_neutral (-5.24) < azm_anion (-5.13) < sfa_neutral (-5.08) < sfa_cation (-5.06) < sfa_anion (-5.03) < paba_anion (-4.98). Reference ligand: `sfa_neutral` at -5.08.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): **no species has a zinc-bound top pose in any seed**, so this scoring function ranks these species on poses that do not reproduce the zinc anchor.

**Rank order, AD4Zn** (best first): smx_neutral (-9.85) < azm_neutral (-9.48) < sfa_neutral (-8.84). Reference ligand: `sfa_neutral` at -8.84.

Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å): azm_neutral -9.50 (5/10 seeds), sfa_neutral -8.84 (20/20 seeds). Species with no zinc-bound seed are absent from that list.

### 3TZF — DHPS, Y. pestis (chain B)

Co-crystal ligand 08D (sulfamethoxazole). Validation: passes the 2.0 Å criterion only (SMX anion, Vina, 1.92 ± 0.04 Å); poses qualitative.

| species | scoring | score (kcal/mol) | best | seeds | clusters | LE | FQ | Δ ref |
|---|---|---|---|---|---|---|---|---|
| sulfanilamide (neutral) | Vina | -6.423 ± 0.009 | -6.438 | 10 | 10/10 | 0.584 | 0.84 | 1.38 |
| sulfanilamide (neutral) | Vinardo | -5.310 ± 0.019 | -5.352 | 10 | 10/10 | 0.483 | 0.69 | 0.69 |
| sulfonamidate anion | Vina | -6.400 ± 0.006 | -6.408 | 10 | 10/10 | 0.582 | 0.83 | 1.40 |
| sulfonamidate anion | Vinardo | -5.267 ± 0.007 | -5.273 | 10 | 10/10 | 0.479 | 0.69 | 0.73 |
| 4-sulfamoylanilinium cation | Vina | -6.403 ± 0.009 | -6.422 | 10 | 10/10 | 0.582 | 0.83 | 1.40 |
| 4-sulfamoylanilinium cation | Vinardo | -5.236 ± 0.004 | -5.242 | 10 | 10/10 | 0.476 | 0.68 | 0.76 |
| sulfamethoxazole anion | Vina | -7.636 ± 0.091 | -7.695 | 10 | 10/10 | 0.449 | 0.85 | 0.16 |
| sulfamethoxazole anion | Vinardo | -5.926 ± 0.200 | -6.078 | 10 | 9/10, 1/10 | 0.349 | 0.66 | 0.07 |
| sulfamethoxazole (neutral) | Vina | -7.799 ± 0.022 | -7.834 | 10 | 10/10 | 0.459 | 0.87 | ref |
| sulfamethoxazole (neutral) | Vinardo | -5.997 ± 0.053 | -6.054 | 10 | 10/10 | 0.353 | 0.67 | ref |
| acetazolamide (neutral) | Vina | -5.582 ± 0.057 | -5.640 | 10 | 10/10 | 0.429 | 0.67 | 2.22 |
| acetazolamide (neutral) | Vinardo | -4.495 ± 0.151 | -4.623 | 10 | 7/10, 3/10 | 0.346 | 0.54 | 1.50 |
| acetazolamide anion | Vina | -5.552 ± 0.028 | -5.589 | 10 | 10/10 | 0.427 | 0.67 | 2.25 |
| acetazolamide anion | Vinardo | -4.450 ± 0.078 | -4.506 | 10 | 8/10, 2/10 | 0.342 | 0.54 | 1.55 |
| p-aminobenzoate | Vina | -6.348 ± 0.015 | -6.369 | 10 | 10/10 | 0.635 | 0.88 | 1.45 |
| p-aminobenzoate | Vinardo | -5.384 ± 0.002 | -5.389 | 10 | 10/10 | 0.538 | 0.75 | 0.61 |

**Rank order, Vina** (best first): smx_neutral (-7.80) < smx_anion (-7.64) < sfa_neutral (-6.42) < sfa_cation (-6.40) < sfa_anion (-6.40) < paba_anion (-6.35) < azm_neutral (-5.58) < azm_anion (-5.55). Reference ligand: `smx_neutral` at -7.80.

**Rank order, Vinardo** (best first): smx_neutral (-6.00) < smx_anion (-5.93) < paba_anion (-5.38) < sfa_neutral (-5.31) < sfa_anion (-5.27) < sfa_cation (-5.24) < azm_neutral (-4.50) < azm_anion (-4.45). Reference ligand: `smx_neutral` at -6.00.

## 3. Where the poses sit relative to the crystal ligands

Symmetry-corrected heavy-atom RMSD in the receptor frame, no superposition (spyrmsd). `full` is used only where the docked species is the same molecule as the crystal ligand; otherwise the shared moiety is compared, over every chemically equivalent mapping of that moiety onto both molecules. These are **not** redocking RMSDs for the sulfanilamide species (only 6RL9 co-crystallises sulfanilamide): they answer "does it bind where the reference ligand binds".

| target | species | scoring | reference | measure | RMSD (Å) | best | f≤1 Å | f≤2 Å |
|---|---|---|---|---|---|---|---|---|
| 7S2K | sfa_neutral | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.04 ± 0.07 | 1.00 | 0.1 | 1.0 |
| 7S2K | sfa_neutral | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.91 ± 0.00 | 0.91 | 1.0 | 1.0 |
| 7S2K | sfa_anion | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.21 ± 0.01 | 1.20 | 0.0 | 1.0 |
| 7S2K | sfa_anion | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.23 ± 0.00 | 1.23 | 0.0 | 1.0 |
| 7S2K | sfa_cation | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.87 ± 0.01 | 0.86 | 1.0 | 1.0 |
| 7S2K | sfa_cation | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.86 ± 0.00 | 0.86 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vina | 78H (7S2K, own) | full molecule (23 atoms) | 0.52 ± 0.03 | 0.47 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.71 ± 0.04 | 0.65 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vina | 78H (7S2K, own) | best pose of the ranked list (any rank) | 0.52 ± 0.03 | 0.47 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vinardo | 78H (7S2K, own) | full molecule (23 atoms) | 0.54 ± 0.01 | 0.52 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.76 ± 0.02 | 0.73 | 1.0 | 1.0 |
| 7S2K | dhp_anion | Vinardo | 78H (7S2K, own) | best pose of the ranked list (any rank) | 0.54 ± 0.01 | 0.52 | 1.0 | 1.0 |
| 7S2K | smx_anion | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.93 ± 0.03 | 0.86 | 1.0 | 1.0 |
| 7S2K | smx_anion | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.86 ± 0.01 | 0.85 | 1.0 | 1.0 |
| 7S2K | smx_neutral | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.90 ± 0.01 | 0.88 | 1.0 | 1.0 |
| 7S2K | smx_neutral | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.85 ± 0.02 | 0.83 | 1.0 | 1.0 |
| 7S2K | paba_anion | Vina | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.43 ± 0.01 | 0.42 | 1.0 | 1.0 |
| 7S2K | paba_anion | Vinardo | 78H (7S2K, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.38 ± 0.00 | 0.37 | 1.0 | 1.0 |
| 3HS4 | sfa_neutral | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.60 ± 0.04 | 1.56 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.42 ± 0.26 | 1.12 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.42 ± 0.26 | 1.12 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.35 ± 0.37 | 0.91 | 0.4 | 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.63 ± 0.04 | 1.58 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.07 ± 0.11 | 0.77 | 0.2 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.64 ± 0.01 | 1.63 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.63 ± 0.00 | 1.63 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.63 ± 0.00 | 1.63 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.68 ± 0.00 | 1.68 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.69 ± 0.01 | 1.69 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.32 ± 0.27 | 1.05 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.69 ± 0.03 | 1.63 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 3.22 ± 0.02 | 3.18 | 0.0 | 0.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 3.22 ± 0.02 | 3.18 | 0.0 | 0.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 3.37 ± 0.02 | 3.33 | 0.0 | 0.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.82 ± 0.03 | 1.77 | 0.0 | 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 0.56 ± 0.05 | 0.50 | 1.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.56 ± 0.00 | 1.56 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.12 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.12 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.90 ± 0.01 | 0.90 | 1.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.57 ± 0.00 | 1.57 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.11 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.66 ± 0.00 | 1.65 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.65 ± 0.00 | 1.65 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.65 ± 0.00 | 1.65 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.71 ± 0.01 | 1.70 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.72 ± 0.00 | 1.71 | 0.0 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.52 ± 0.23 | 1.07 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.64 ± 0.02 | 1.60 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.64 ± 0.02 | 1.61 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.64 ± 0.02 | 1.61 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.65 ± 0.02 | 1.61 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.68 ± 0.02 | 1.64 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.05 ± 0.23 | 0.38 | 0.1 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.64 ± 0.00 | 1.63 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | full molecule (11 atoms) | 1.63 ± 0.01 | 1.62 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 1.63 ± 0.01 | 1.62 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 1.68 ± 0.01 | 1.67 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.69 ± 0.00 | 1.69 | 0.0 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | best pose of the ranked list (any rank) | 1.39 ± 0.24 | 1.14 | 0.0 | 1.0 |
| 3HS4 | smx_anion | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 4.77 ± 0.05 | 4.68 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 4.14 ± 0.01 | 4.12 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.07 ± 0.01 | 4.05 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.81 ± 0.05 | 4.72 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 4.71 ± 0.20 | 4.50 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 4.66 ± 0.45 | 4.33 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.59 ± 0.59 | 4.20 | 0.0 | 0.0 |
| 3HS4 | smx_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.76 ± 0.22 | 4.54 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 4.81 ± 0.03 | 4.77 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 4.15 ± 0.01 | 4.14 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.09 ± 0.01 | 4.08 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.86 ± 0.03 | 4.82 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 4.97 ± 0.18 | 4.52 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 5.24 ± 0.33 | 4.30 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 5.34 ± 0.42 | 4.17 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 5.04 ± 0.20 | 4.56 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 4.08 ± 0.04 | 4.03 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-SO2N core (11 atoms) | 4.19 ± 0.30 | 3.73 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.32 ± 0.32 | 3.84 | 0.0 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.15 ± 0.04 | 4.09 | 0.0 | 0.0 |
| 3HS4 | azm_neutral | Vina | AZM (3HS4, own) | full molecule (13 atoms) | 5.39 ± 0.09 | 5.29 | 0.0 | 0.0 |
| 3HS4 | azm_neutral | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 6.37 ± 0.27 | 6.14 | 0.0 | 0.0 |
| 3HS4 | azm_neutral | Vina | AZM (3HS4, own) | best pose of the ranked list (any rank) | 1.40 ± 0.24 | 0.95 | 0.1 | 1.0 |
| 3HS4 | azm_neutral | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 6.42 ± 0.28 | 6.20 | 0.0 | 0.0 |
| 3HS4 | azm_neutral | Vinardo | AZM (3HS4, own) | full molecule (13 atoms) | 1.35 ± 0.01 | 1.33 | 0.0 | 1.0 |
| 3HS4 | azm_neutral | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.62 ± 0.00 | 1.61 | 0.0 | 1.0 |
| 3HS4 | azm_neutral | Vinardo | AZM (3HS4, own) | best pose of the ranked list (any rank) | 1.35 ± 0.01 | 1.33 | 0.0 | 1.0 |
| 3HS4 | azm_neutral | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.65 ± 0.00 | 1.64 | 0.0 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | AZM (3HS4, own) | full molecule (13 atoms) | 2.70 ± 0.20 | 2.31 | 0.0 | 0.0 |
| 3HS4 | azm_neutral | AD4Zn | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 0.57 ± 0.02 | 0.54 | 1.0 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | AZM (3HS4, own) | best pose of the ranked list (any rank) | 1.04 ± 0.60 | 0.42 | 0.5 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 0.49 ± 0.02 | 0.46 | 1.0 | 1.0 |
| 3HS4 | azm_anion | Vina | AZM (3HS4, own) | full molecule (13 atoms) | 5.37 ± 0.02 | 5.35 | 0.0 | 0.0 |
| 3HS4 | azm_anion | Vina | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 6.26 ± 0.04 | 6.18 | 0.0 | 0.0 |
| 3HS4 | azm_anion | Vina | AZM (3HS4, own) | best pose of the ranked list (any rank) | 1.34 ± 0.17 | 1.02 | 0.0 | 1.0 |
| 3HS4 | azm_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 6.31 ± 0.04 | 6.24 | 0.0 | 0.0 |
| 3HS4 | azm_anion | Vinardo | AZM (3HS4, own) | full molecule (13 atoms) | 1.38 ± 0.01 | 1.36 | 0.0 | 1.0 |
| 3HS4 | azm_anion | Vinardo | AZM (3HS4, own) | C-SO2N zinc-binding group (5 atoms) | 1.64 ± 0.01 | 1.63 | 0.0 | 1.0 |
| 3HS4 | azm_anion | Vinardo | AZM (3HS4, own) | best pose of the ranked list (any rank) | 1.38 ± 0.01 | 1.36 | 0.0 | 1.0 |
| 3HS4 | azm_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.67 ± 0.01 | 1.66 | 0.0 | 1.0 |
| 3HS4 | paba_anion | Vina | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.13 ± 0.00 | 4.13 | 0.0 | 0.0 |
| 3HS4 | paba_anion | Vinardo | SAN (6RL9, carried into 3HS4 by Ca superposition) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.14 ± 0.00 | 4.14 | 0.0 | 0.0 |
| 6RL9 | sfa_neutral | Vina | SAN (6RL9, own) | full molecule (11 atoms) | 1.08 ± 0.00 | 1.08 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vina | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.08 ± 0.00 | 1.08 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.86 ± 0.01 | 0.85 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.58 ± 0.00 | 1.57 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vina | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.00 ± 0.24 | 0.33 | 0.1 | 1.0 |
| 6RL9 | sfa_neutral | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.56 ± 0.00 | 1.56 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN (6RL9, own) | full molecule (11 atoms) | 1.11 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.11 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.98 ± 0.00 | 0.98 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.60 ± 0.00 | 1.60 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.11 ± 0.01 | 1.10 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.58 ± 0.00 | 1.58 | 0.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN (6RL9, own) | full molecule (11 atoms) | 0.61 ± 0.04 | 0.56 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 0.61 ± 0.04 | 0.56 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.63 ± 0.04 | 0.58 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 0.45 ± 0.02 | 0.42 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN (6RL9, own) | best pose of the ranked list (any rank) | 0.57 ± 0.06 | 0.47 | 1.0 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 0.58 ± 0.02 | 0.56 | 1.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | SAN (6RL9, own) | full molecule (11 atoms) | 1.07 ± 0.00 | 1.07 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.07 ± 0.00 | 1.07 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.84 ± 0.01 | 0.83 | 1.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.57 ± 0.00 | 1.56 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.07 ± 0.00 | 1.07 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.56 ± 0.00 | 1.56 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN (6RL9, own) | full molecule (11 atoms) | 1.09 ± 0.00 | 1.09 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.09 ± 0.00 | 1.09 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.95 ± 0.00 | 0.94 | 1.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.59 ± 0.00 | 1.58 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.09 ± 0.00 | 1.09 | 0.0 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.57 ± 0.00 | 1.57 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vina | SAN (6RL9, own) | full molecule (11 atoms) | 1.08 ± 0.00 | 1.08 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vina | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.08 ± 0.00 | 1.08 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.85 ± 0.01 | 0.85 | 1.0 | 1.0 |
| 6RL9 | sfa_cation | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.58 ± 0.00 | 1.57 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vina | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.00 ± 0.25 | 0.28 | 0.1 | 1.0 |
| 6RL9 | sfa_cation | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.57 ± 0.00 | 1.56 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN (6RL9, own) | full molecule (11 atoms) | 1.11 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 1.11 ± 0.00 | 1.11 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.98 ± 0.00 | 0.98 | 1.0 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.60 ± 0.00 | 1.60 | 0.0 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN (6RL9, own) | best pose of the ranked list (any rank) | 1.08 ± 0.10 | 0.79 | 0.1 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.58 ± 0.00 | 1.58 | 0.0 | 1.0 |
| 6RL9 | smx_anion | Vina | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 4.14 ± 0.01 | 4.13 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.06 ± 0.01 | 4.05 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 4.86 ± 0.01 | 4.84 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.81 ± 0.01 | 4.79 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vinardo | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 5.36 ± 0.01 | 5.34 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 5.44 ± 0.01 | 5.43 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 5.13 ± 0.02 | 5.09 | 0.0 | 0.0 |
| 6RL9 | smx_anion | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 5.04 ± 0.02 | 5.00 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vina | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 4.61 ± 0.71 | 4.16 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.60 ± 0.79 | 4.09 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 5.05 ± 0.18 | 4.90 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.98 ± 0.16 | 4.85 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 5.42 ± 0.06 | 5.34 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 5.56 ± 0.11 | 5.43 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 5.12 ± 0.03 | 5.09 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 5.05 ± 0.01 | 5.03 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | SAN (6RL9, own) | aniline-SO2N core (11 atoms) | 4.25 ± 0.38 | 3.98 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.40 ± 0.40 | 4.12 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 4.23 ± 0.14 | 4.08 | 0.0 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 4.16 ± 0.14 | 4.02 | 0.0 | 0.0 |
| 6RL9 | azm_neutral | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 5.73 ± 0.02 | 5.68 | 0.0 | 0.0 |
| 6RL9 | azm_neutral | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | full molecule (13 atoms) | 5.23 ± 0.01 | 5.22 | 0.0 | 0.0 |
| 6RL9 | azm_neutral | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 5.69 ± 0.02 | 5.64 | 0.0 | 0.0 |
| 6RL9 | azm_neutral | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | best pose of the ranked list (any rank) | 1.51 ± 0.15 | 1.35 | 0.0 | 1.0 |
| 6RL9 | azm_neutral | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 3.10 ± 2.26 | 1.66 | 0.0 | 0.7 |
| 6RL9 | azm_neutral | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | full molecule (13 atoms) | 2.64 ± 1.97 | 1.37 | 0.0 | 0.7 |
| 6RL9 | azm_neutral | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 3.06 ± 2.25 | 1.62 | 0.0 | 0.7 |
| 6RL9 | azm_neutral | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | best pose of the ranked list (any rank) | 1.41 ± 0.03 | 1.36 | 0.0 | 1.0 |
| 6RL9 | azm_neutral | AD4Zn | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 1.37 ± 0.78 | 0.60 | 0.5 | 0.5 |
| 6RL9 | azm_neutral | AD4Zn | AZM (3HS4, carried into 6RL9 by Ca superposition) | full molecule (13 atoms) | 2.68 ± 0.81 | 1.55 | 0.0 | 0.3 |
| 6RL9 | azm_neutral | AD4Zn | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 1.40 ± 0.77 | 0.51 | 0.5 | 0.5 |
| 6RL9 | azm_neutral | AD4Zn | AZM (3HS4, carried into 6RL9 by Ca superposition) | best pose of the ranked list (any rank) | 0.99 ± 0.42 | 0.45 | 0.5 | 1.0 |
| 6RL9 | azm_anion | Vina | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 5.74 ± 0.02 | 5.71 | 0.0 | 0.0 |
| 6RL9 | azm_anion | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | full molecule (13 atoms) | 5.22 ± 0.01 | 5.21 | 0.0 | 0.0 |
| 6RL9 | azm_anion | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 5.70 ± 0.02 | 5.67 | 0.0 | 0.0 |
| 6RL9 | azm_anion | Vina | AZM (3HS4, carried into 6RL9 by Ca superposition) | best pose of the ranked list (any rank) | 1.46 ± 0.32 | 0.90 | 0.1 | 1.0 |
| 6RL9 | azm_anion | Vinardo | SAN (6RL9, own) | C-SO2N zinc-binding group (5 atoms) | 2.66 ± 1.96 | 1.70 | 0.0 | 0.8 |
| 6RL9 | azm_anion | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | full molecule (13 atoms) | 2.33 ± 1.68 | 1.45 | 0.0 | 0.8 |
| 6RL9 | azm_anion | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | C-SO2N zinc-binding group (5 atoms) | 2.62 ± 1.95 | 1.68 | 0.0 | 0.8 |
| 6RL9 | azm_anion | Vinardo | AZM (3HS4, carried into 6RL9 by Ca superposition) | best pose of the ranked list (any rank) | 1.40 ± 0.17 | 1.05 | 0.0 | 1.0 |
| 6RL9 | paba_anion | Vina | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 3.39 ± 1.53 | 0.48 | 0.2 | 0.2 |
| 6RL9 | paba_anion | Vinardo | SAN (6RL9, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 4.13 ± 0.00 | 4.13 | 0.0 | 0.0 |
| 3TZF | sfa_neutral | Vina | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.20 ± 0.03 | 0.15 | 1.0 | 1.0 |
| 3TZF | sfa_neutral | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.18 ± 0.02 | 0.13 | 1.0 | 1.0 |
| 3TZF | sfa_neutral | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 0.28 ± 0.04 | 0.20 | 1.0 | 1.0 |
| 3TZF | sfa_neutral | Vinardo | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.96 ± 0.32 | 0.34 | 0.2 | 1.0 |
| 3TZF | sfa_neutral | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.70 ± 0.19 | 0.32 | 1.0 | 1.0 |
| 3TZF | sfa_neutral | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 1.38 ± 0.52 | 0.38 | 0.2 | 1.0 |
| 3TZF | sfa_anion | Vina | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.22 ± 0.01 | 0.20 | 1.0 | 1.0 |
| 3TZF | sfa_anion | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.19 ± 0.01 | 0.18 | 1.0 | 1.0 |
| 3TZF | sfa_anion | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 0.30 ± 0.02 | 0.28 | 1.0 | 1.0 |
| 3TZF | sfa_anion | Vinardo | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.32 ± 0.01 | 0.32 | 1.0 | 1.0 |
| 3TZF | sfa_anion | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.29 ± 0.01 | 0.29 | 1.0 | 1.0 |
| 3TZF | sfa_anion | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 0.42 ± 0.02 | 0.41 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vina | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.21 ± 0.03 | 0.17 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.18 ± 0.02 | 0.15 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 0.29 ± 0.04 | 0.23 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vinardo | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.30 ± 0.00 | 0.30 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.28 ± 0.00 | 0.27 | 1.0 | 1.0 |
| 3TZF | sfa_cation | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 0.38 ± 0.00 | 0.37 | 1.0 | 1.0 |
| 3TZF | smx_anion | Vina | 08D (3TZF, own) | full molecule (17 atoms) | 1.92 ± 0.04 | 1.80 | 0.0 | 1.0 |
| 3TZF | smx_anion | Vina | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 1.09 ± 0.01 | 1.09 | 0.0 | 1.0 |
| 3TZF | smx_anion | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.83 ± 0.02 | 0.81 | 1.0 | 1.0 |
| 3TZF | smx_anion | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 1.61 ± 0.01 | 1.60 | 0.0 | 1.0 |
| 3TZF | smx_anion | Vina | 08D (3TZF, own) | best pose of the ranked list (any rank) | 1.43 ± 0.44 | 0.56 | 0.1 | 1.0 |
| 3TZF | smx_anion | Vinardo | 08D (3TZF, own) | full molecule (17 atoms) | 2.18 ± 0.40 | 1.42 | 0.0 | 0.1 |
| 3TZF | smx_anion | Vinardo | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 1.05 ± 0.22 | 0.47 | 0.2 | 1.0 |
| 3TZF | smx_anion | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.80 ± 0.15 | 0.39 | 1.0 | 1.0 |
| 3TZF | smx_anion | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 1.50 ± 0.35 | 0.58 | 0.1 | 1.0 |
| 3TZF | smx_anion | Vinardo | 08D (3TZF, own) | best pose of the ranked list (any rank) | 1.44 ± 0.31 | 0.98 | 0.1 | 0.9 |
| 3TZF | smx_neutral | Vina | 08D (3TZF, own) | full molecule (17 atoms) | 3.00 ± 0.04 | 2.94 | 0.0 | 0.0 |
| 3TZF | smx_neutral | Vina | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.78 ± 0.01 | 0.77 | 1.0 | 1.0 |
| 3TZF | smx_neutral | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.64 ± 0.01 | 0.62 | 1.0 | 1.0 |
| 3TZF | smx_neutral | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 1.11 ± 0.01 | 1.09 | 0.0 | 1.0 |
| 3TZF | smx_neutral | Vina | 08D (3TZF, own) | best pose of the ranked list (any rank) | 0.98 ± 0.37 | 0.54 | 0.4 | 1.0 |
| 3TZF | smx_neutral | Vinardo | 08D (3TZF, own) | full molecule (17 atoms) | 3.00 ± 0.03 | 2.94 | 0.0 | 0.0 |
| 3TZF | smx_neutral | Vinardo | 08D (3TZF, own) | aniline-SO2N core (11 atoms) | 0.81 ± 0.08 | 0.58 | 1.0 | 1.0 |
| 3TZF | smx_neutral | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.68 ± 0.05 | 0.55 | 1.0 | 1.0 |
| 3TZF | smx_neutral | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 1.12 ± 0.12 | 0.78 | 0.1 | 1.0 |
| 3TZF | smx_neutral | Vinardo | 08D (3TZF, own) | best pose of the ranked list (any rank) | 1.18 ± 0.37 | 0.52 | 0.2 | 1.0 |
| 3TZF | azm_neutral | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 6.85 ± 0.06 | 6.74 | 0.0 | 0.0 |
| 3TZF | azm_neutral | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 6.37 ± 0.23 | 6.20 | 0.0 | 0.0 |
| 3TZF | azm_anion | Vina | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 6.78 ± 0.12 | 6.49 | 0.0 | 0.0 |
| 3TZF | azm_anion | Vinardo | 08D (3TZF, own) | C-SO2N zinc-binding group (5 atoms) | 6.29 ± 0.17 | 6.16 | 0.0 | 0.0 |
| 3TZF | paba_anion | Vina | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.31 ± 0.03 | 0.28 | 1.0 | 1.0 |
| 3TZF | paba_anion | Vinardo | 08D (3TZF, own) | aniline-X(O)2 core, X = S or C (10 atoms) | 0.44 ± 0.00 | 0.44 | 1.0 | 1.0 |


Where a whole-molecule reference exists, the rank at which the first near-native pose appears in each seed's ranked list says whether a poor top-pose RMSD is a sampling failure or a ranking one:

| target | species | scoring | reference | seeds with a pose ≤1 Å anywhere | median rank of it | same for ≤2 Å |
|---|---|---|---|---|---|---|
| 7S2K | dhp_anion | Vina | 78H | 1.00 | 1.0 | 1.00 at rank 1.0 |
| 7S2K | dhp_anion | Vinardo | 78H | 1.00 | 1.0 | 1.00 at rank 1.0 |
| 3HS4 | sfa_neutral | Vina | SAN | 0.20 | 7.5 | 1.00 at rank 1.0 |
| 3HS4 | sfa_neutral | Vinardo | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | SAN | 1.00 | 2.0 | 1.00 at rank 2.0 |
| 3HS4 | sfa_anion | Vina | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 3HS4 | sfa_anion | Vinardo | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 3HS4 | sfa_cation | Vina | SAN | 0.10 | 3.0 | 1.00 at rank 1.0 |
| 3HS4 | sfa_cation | Vinardo | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 3HS4 | azm_neutral | Vina | AZM | 0.10 | 6.0 | 1.00 at rank 5.5 |
| 3HS4 | azm_neutral | Vinardo | AZM | 0.00 | – | 1.00 at rank 1.0 |
| 3HS4 | azm_neutral | AD4Zn | AZM | 0.50 | 7.0 | 1.00 at rank 2.0 |
| 3HS4 | azm_anion | Vina | AZM | 0.00 | – | 1.00 at rank 6.0 |
| 3HS4 | azm_anion | Vinardo | AZM | 0.00 | – | 1.00 at rank 1.0 |
| 6RL9 | sfa_neutral | Vina | SAN | 0.10 | 3.0 | 1.00 at rank 1.0 |
| 6RL9 | sfa_neutral | Vinardo | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | SAN | 1.00 | 1.0 | 1.00 at rank 1.0 |
| 6RL9 | sfa_anion | Vina | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 6RL9 | sfa_anion | Vinardo | SAN | 0.00 | – | 1.00 at rank 1.0 |
| 6RL9 | sfa_cation | Vina | SAN | 0.10 | 3.0 | 1.00 at rank 1.0 |
| 6RL9 | sfa_cation | Vinardo | SAN | 0.10 | 7.0 | 1.00 at rank 1.0 |
| 6RL9 | azm_neutral | Vina | AZM | 0.00 | – | 1.00 at rank 7.0 |
| 6RL9 | azm_neutral | Vinardo | AZM | 0.00 | – | 1.00 at rank 1.0 |
| 6RL9 | azm_neutral | AD4Zn | AZM | 0.50 | 7.0 | 1.00 at rank 2.0 |
| 6RL9 | azm_anion | Vina | AZM | 0.10 | 13.0 | 1.00 at rank 9.0 |
| 6RL9 | azm_anion | Vinardo | AZM | 0.00 | – | 1.00 at rank 1.0 |
| 3TZF | smx_anion | Vina | 08D | 0.10 | 3.0 | 1.00 at rank 1.0 |
| 3TZF | smx_anion | Vinardo | 08D | 0.10 | 9.0 | 0.90 at rank 2.0 |
| 3TZF | smx_neutral | Vina | 08D | 0.40 | 4.0 | 1.00 at rank 2.0 |
| 3TZF | smx_neutral | Vinardo | 08D | 0.20 | 3.0 | 1.00 at rank 2.0 |

### The validation gate, for the same targets

Gate G2 as defined before the runs: top-pose symmetry-corrected RMSD ≤ **1.0 Å** (mean over seeds) with the success fraction reported, PoseBusters physics ≥ 0.8, and for hCA II also Zn–N(sulfonamide) ≤ 2.3 Å and the Thr199 hydrogen bond in ≥ half the seeds. The conventional **2.0 Å** criterion (Hevener 2009; Warren 2006) is evaluated with the same structure. Both are reported; the gate was not changed.

| validation run (redock / cross-dock) | RMSD top (Å) | f≤1 Å | f≤2 Å | G2 (1.0 Å) | 2.0 Å criterion |
|---|---|---|---|---|---|
| `T1_3TZF/smx_anion__vina` | 1.92 ± 0.04 | 0.0 | 1.0 | fail | pass |
| `T1_3TZF/smx_anion__vinardo` | 2.18 ± 0.40 | 0.0 | 0.1 | fail | fail |
| `T1_3TZF/smx_neutral__vina` | 3.00 ± 0.04 | 0.0 | 0.0 | fail | fail |
| `T3_6RL9/azm_anion__vina` | 5.22 ± 0.01 | 0.0 | 0.0 | fail | fail |
| `T3_6RL9/azm_neutral__ad4zn` | 2.68 ± 0.81 | 0.0 | 0.3 | fail | fail |
| `T3_6RL9/azm_neutral__vina` | 5.23 ± 0.01 | 0.0 | 0.0 | fail | fail |
| `T3_6RL9/sfa_anion__ad4zn` | 1.12 ± 0.32 | 0.3 | 1.0 | fail | pass |
| `T3_6RL9/sfa_anion__vina` | 1.09 ± 0.00 | 0.0 | 1.0 | fail | fail |
| `T3_6RL9/sfa_neutral__ad4zn` | 2.07 ± 1.45 | 0.5 | 0.5 | fail | fail |
| `T3_6RL9/sfa_neutral__vina` | 1.10 ± 0.01 | 0.0 | 1.0 | fail | fail |
| `T3_6RL9/sfa_neutral__vinardo` | 1.13 ± 0.00 | 0.0 | 1.0 | fail | fail |
| `T4_3HS4/azm_anion__ad4zn` | 3.65 ± 0.44 | 0.0 | 0.0 | fail | fail |
| `T4_3HS4/azm_anion__vina` | 5.37 ± 0.02 | 0.0 | 0.0 | fail | fail |
| `T4_3HS4/azm_neutral__ad4zn` | 2.70 ± 0.20 | 0.0 | 0.0 | fail | fail |
| `T4_3HS4/azm_neutral__vina` | 5.39 ± 0.09 | 0.0 | 0.0 | fail | fail |
| `T4_3HS4/azm_neutral__vinardo` | 1.35 ± 0.01 | 0.0 | 1.0 | fail | fail |
| `T4_3HS4/sfa_anion__vina` | 1.17 ± 0.16 | 0.0 | 1.0 | fail | fail |
| `T4_3HS4/sfa_neutral__ad4zn` | 0.62 ± 0.03 | 1.0 | 1.0 | **PASS** | pass |
| `T4_3HS4/sfa_neutral__vina` | 1.62 ± 0.01 | 0.0 | 1.0 | fail | fail |
| `T5_7S2K/dhp_anion__vina` | 0.51 ± 0.03 | 1.0 | 1.0 | **PASS** | pass |
| `T5_7S2K/dhp_anion__vinardo` | 0.54 ± 0.01 | 1.0 | 1.0 | **PASS** | pass |

## 4. Site checks and key contacts

Site checks are the geometric tests of `analyse.py` (plus two generalised ones for species whose SMARTS the validation set did not cover), each also evaluated on the receptor's own crystal ligand so that a mis-specified check shows up as a failing crystal pose. Distances are heavy-atom to heavy-atom.

| target | species | scoring | check | crystal ligand | docked mean | pass fraction |
|---|---|---|---|---|---|---|
| 7S2K | sfa_neutral | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.61 | 0.1 |
| 7S2K | sfa_neutral | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.13 | 1.0 |
| 7S2K | sfa_anion | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 2.98 | 1.0 |
| 7S2K | sfa_anion | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 2.94 | 1.0 |
| 7S2K | sfa_cation | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.06 | 1.0 |
| 7S2K | sfa_cation | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.04 | 1.0 |
| 7S2K | dhp_anion | Vina | pterin 2-amino N to Asp174 (≤ 3.2 Å) | 2.60 | 2.65 | 1.0 |
| 7S2K | dhp_anion | Vina | pABA carboxylate to Ser214 N (≤ 3.4 Å) | 2.94 | 2.90 | 1.0 |
| 7S2K | dhp_anion | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 2.90 | 1.0 |
| 7S2K | dhp_anion | Vinardo | pterin 2-amino N to Asp174 (≤ 3.2 Å) | 2.60 | 2.68 | 1.0 |
| 7S2K | dhp_anion | Vinardo | pABA carboxylate to Ser214 N (≤ 3.4 Å) | 2.94 | 2.96 | 1.0 |
| 7S2K | dhp_anion | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 2.96 | 1.0 |
| 7S2K | smx_anion | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.26 | 1.0 |
| 7S2K | smx_anion | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.15 | 1.0 |
| 7S2K | smx_neutral | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.16 | 1.0 |
| 7S2K | smx_neutral | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.06 | 1.0 |
| 7S2K | azm_neutral | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 5.06 | 0.0 |
| 7S2K | azm_neutral | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 5.39 | 0.0 |
| 7S2K | azm_anion | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 5.15 | 0.0 |
| 7S2K | azm_anion | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 5.41 | 0.0 |
| 7S2K | paba_anion | Vina | pABA carboxylate to Ser214 N (≤ 3.4 Å) | 2.94 | 3.04 | 1.0 |
| 7S2K | paba_anion | Vina | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 3.04 | 1.0 |
| 7S2K | paba_anion | Vinardo | pABA carboxylate to Ser214 N (≤ 3.4 Å) | 2.94 | 2.96 | 1.0 |
| 7S2K | paba_anion | Vinardo | sulfonyl or carboxylate O to Ser214 N (≤ 3.4 Å) | 2.94 | 2.96 | 1.0 |
| 3HS4 | sfa_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.59 | 0.0 |
| 3HS4 | sfa_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 4.01 | 0.0 |
| 3HS4 | sfa_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.88 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.30 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 3.13 | 1.0 |
| 3HS4 | sfa_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.73 | 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.32 | 0.3 |
| 3HS4 | sfa_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 2.46 | 1.0 |
| 3HS4 | sfa_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.84 | 1.0 |
| 3HS4 | sfa_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.74 | 0.0 |
| 3HS4 | sfa_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 4.71 | 0.0 |
| 3HS4 | sfa_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.96 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.30 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 3.10 | 1.0 |
| 3HS4 | sfa_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.74 | 1.0 |
| 3HS4 | sfa_cation | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.50 | 0.0 |
| 3HS4 | sfa_cation | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 3.58 | 0.0 |
| 3HS4 | sfa_cation | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.83 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.30 | 0.9 |
| 3HS4 | sfa_cation | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 3.17 | 1.0 |
| 3HS4 | sfa_cation | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.73 | 1.0 |
| 3HS4 | smx_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 5.16 | 0.0 |
| 3HS4 | smx_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 7.28 | 0.0 |
| 3HS4 | smx_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 6.76 | 0.0 |
| 3HS4 | smx_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 6.21 | 0.0 |
| 3HS4 | smx_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 7.20 | 0.0 |
| 3HS4 | smx_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 7.11 | 0.0 |
| 3HS4 | smx_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 5.12 | 0.0 |
| 3HS4 | smx_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 7.14 | 0.0 |
| 3HS4 | smx_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 6.82 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 5.10 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 5.65 | 0.0 |
| 3HS4 | smx_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 7.40 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 3.28 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 4.84 | 0.0 |
| 3HS4 | smx_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 6.86 | 0.0 |
| 3HS4 | azm_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 6.09 | 0.0 |
| 3HS4 | azm_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 8.22 | 0.0 |
| 3HS4 | azm_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 8.76 | 0.0 |
| 3HS4 | azm_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.54 | 0.0 |
| 3HS4 | azm_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 4.54 | 0.0 |
| 3HS4 | azm_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.81 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.08 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 2.58 | 1.0 |
| 3HS4 | azm_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.88 | 1.0 |
| 3HS4 | azm_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 6.00 | 0.0 |
| 3HS4 | azm_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 8.15 | 0.0 |
| 3HS4 | azm_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 8.63 | 0.0 |
| 3HS4 | azm_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.94 | 2.60 | 0.0 |
| 3HS4 | azm_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.82 | 4.63 | 0.0 |
| 3HS4 | azm_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.95 | 2.76 | 1.0 |
| 6RL9 | sfa_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.75 | 0.0 |
| 6RL9 | sfa_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.68 | 0.0 |
| 6RL9 | sfa_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.96 | 1.0 |
| 6RL9 | sfa_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.40 | 0.0 |
| 6RL9 | sfa_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.35 | 0.0 |
| 6RL9 | sfa_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.76 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.20 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 2.55 | 1.0 |
| 6RL9 | sfa_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.97 | 1.0 |
| 6RL9 | sfa_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.75 | 0.0 |
| 6RL9 | sfa_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.72 | 0.0 |
| 6RL9 | sfa_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.94 | 1.0 |
| 6RL9 | sfa_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.43 | 0.0 |
| 6RL9 | sfa_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.41 | 0.0 |
| 6RL9 | sfa_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.73 | 1.0 |
| 6RL9 | sfa_cation | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.76 | 0.0 |
| 6RL9 | sfa_cation | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.70 | 0.0 |
| 6RL9 | sfa_cation | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.99 | 1.0 |
| 6RL9 | sfa_cation | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 2.41 | 0.0 |
| 6RL9 | sfa_cation | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.36 | 0.0 |
| 6RL9 | sfa_cation | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.77 | 1.0 |
| 6RL9 | smx_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 5.18 | 0.0 |
| 6RL9 | smx_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 7.21 | 0.0 |
| 6RL9 | smx_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 6.66 | 0.0 |
| 6RL9 | smx_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 6.01 | 0.0 |
| 6RL9 | smx_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 6.00 | 0.0 |
| 6RL9 | smx_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 7.16 | 0.0 |
| 6RL9 | smx_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 5.42 | 0.0 |
| 6RL9 | smx_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 6.83 | 0.0 |
| 6RL9 | smx_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 7.00 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 4.97 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 5.48 | 0.0 |
| 6RL9 | smx_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 7.50 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 3.56 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.21 | 0.0 |
| 6RL9 | smx_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 6.71 | 0.0 |
| 6RL9 | azm_neutral | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 5.68 | 0.0 |
| 6RL9 | azm_neutral | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 7.99 | 0.0 |
| 6RL9 | azm_neutral | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 7.93 | 0.0 |
| 6RL9 | azm_neutral | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 3.62 | 0.0 |
| 6RL9 | azm_neutral | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 5.64 | 0.0 |
| 6RL9 | azm_neutral | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 4.54 | 0.7 |
| 6RL9 | azm_neutral | AD4Zn | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 3.68 | 0.5 |
| 6RL9 | azm_neutral | AD4Zn | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 4.22 | 0.5 |
| 6RL9 | azm_neutral | AD4Zn | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 2.93 | 1.0 |
| 6RL9 | azm_anion | Vina | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 5.61 | 0.0 |
| 6RL9 | azm_anion | Vina | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 7.89 | 0.0 |
| 6RL9 | azm_anion | Vina | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 7.88 | 0.0 |
| 6RL9 | azm_anion | Vinardo | Zn-N(sulfonamide) (≤ 2.3 Å) | 1.93 | 3.24 | 0.0 |
| 6RL9 | azm_anion | Vinardo | Thr199 OG1 H-bond (≤ 3.3 Å) | 2.86 | 5.15 | 0.0 |
| 6RL9 | azm_anion | Vinardo | Thr199 N H-bond to sulfonyl O (≤ 3.3 Å) | 2.94 | 3.93 | 0.8 |
| 3TZF | sfa_neutral | Vina | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.16 | 1.0 |
| 3TZF | sfa_neutral | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.16 | 1.0 |
| 3TZF | sfa_neutral | Vinardo | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.30 | 1.0 |
| 3TZF | sfa_neutral | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.30 | 1.0 |
| 3TZF | sfa_anion | Vina | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.15 | 1.0 |
| 3TZF | sfa_anion | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.15 | 1.0 |
| 3TZF | sfa_anion | Vinardo | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.26 | 1.0 |
| 3TZF | sfa_anion | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.26 | 1.0 |
| 3TZF | sfa_cation | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.13 | 1.0 |
| 3TZF | sfa_cation | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.26 | 1.0 |
| 3TZF | smx_anion | Vina | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.20 | 1.0 |
| 3TZF | smx_anion | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.20 | 1.0 |
| 3TZF | smx_anion | Vinardo | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.48 | 0.8 |
| 3TZF | smx_anion | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.48 | 0.8 |
| 3TZF | smx_neutral | Vina | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 2.97 | 1.0 |
| 3TZF | smx_neutral | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 2.97 | 1.0 |
| 3TZF | smx_neutral | Vinardo | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 2.92 | 1.0 |
| 3TZF | smx_neutral | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 2.92 | 1.0 |
| 3TZF | paba_anion | Vina | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.30 | 1.0 |
| 3TZF | paba_anion | Vina | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.30 | 1.0 |
| 3TZF | paba_anion | Vinardo | aniline N to pterin (HH2) (≤ 3.6 Å) | 3.10 | 3.19 | 1.0 |
| 3TZF | paba_anion | Vinardo | aniline N (any protonation) to pterin HH2 N/O (≤ 3.6 Å) | 3.10 | 3.19 | 1.0 |

### Contacts of the cluster representatives

One profile per cluster (the top pose of the cluster's best-scoring seed), from PLIP, ProLIF and plain geometry independently; full records in `dock/results/production_contacts.json`. Hydrogen bonds are given as donor···acceptor heavy-atom distances (the hydrogens are ours, not the crystal's).

**7S2K**

| species | scoring | cluster | population | best score | PLIP H-bonds (D···A, Å) | closest N/O···N/O pairs (Å) | nearest ligand atom to Zn | π contacts |
|---|---|---|---|---|---|---|---|---|
| sfa_neutral | Vina | 1 | 10/10 | -6.44 | SER214 N···N(sulfonamide)#6 3.34; PHE179 N···N(sulfonamide)#6 3.49; SER52 N(aniline)#0···O 3.57; SER214 N(sulfonamide)#6···OG 3.81 | SER52:OG···N(aniline)#0 2.83; POP302:O1···N(aniline)#0 3.01; SER214:N···O(sulfonyl)#8 3.02; SER214:OG···O(sulfonyl)#7 3.15 | – | PHE18 stacking-T 4.79 |
| sfa_neutral | Vinardo | 1 | 10/10 | -5.00 | SER214 N···O(sulfonyl)#7 3.12; SER52 N(aniline)#0···O 3.68; PHE179 N···O(sulfonyl)#7 3.71; ARG212 N(sulfonamide)#6···O 3.73 | SER214:N···N(sulfonamide)#6 2.87; SER214:OG···O(sulfonyl)#8 2.99; SER52:OG···N(aniline)#0 3.06; POP302:O1···N(aniline)#0 3.11 | – | PHE18 stacking-T 4.76 |
| sfa_anion | Vina | 1 | 10/10 | -6.70 | SER214 N···N(sulfonamide)#0 3.32; PHE179 N···N(sulfonamide)#0 3.50; SER52 N(aniline)#8···O 3.68; SER214 N(sulfonamide)#0···OG 3.82 | SER52:OG···N(aniline)#8 2.86; POP302:O1···N(aniline)#8 2.96; SER214:N···O(sulfonyl)#3 3.00; SER214:OG···O(sulfonyl)#2 3.23 | – | – |
| sfa_anion | Vinardo | 1 | 10/10 | -5.27 | SER214 N···N(sulfonamide)#0 3.12; SER52 N(aniline)#8···O 3.54; PHE179 N···N(sulfonamide)#0 3.63; SER214 N(sulfonamide)#0···OG 3.68 | SER52:OG···N(aniline)#8 2.92; SER214:N···O(sulfonyl)#3 2.94; SER214:OG···O(sulfonyl)#2 3.04; POP302:O1···N(aniline)#8 3.05 | – | – |
| sfa_cation | Vina | 1 | 10/10 | -6.39 | SER52 N(aniline)#8···OG 2.80; PHE179 N···O(sulfonyl)#3 3.33; SER214 N···O(sulfonyl)#3 3.71 | SER52:OG···N(aniline)#8 2.80; POP302:O1···N(aniline)#8 2.87; SER214:N···O(sulfonyl)#2 3.07; SER214:OG···N(sulfonamide)#0 3.23 | – | – |
| sfa_cation | Vinardo | 1 | 10/10 | -5.03 | PHE179 N···O(sulfonyl)#3 3.48; SER214 N···O(sulfonyl)#3 3.50; SER52 N(aniline)#8···O 3.58 | SER52:OG···N(aniline)#8 2.77; POP302:O1···N(aniline)#8 2.95; SER214:OG···N(sulfonamide)#0 3.01; SER214:N···O(sulfonyl)#2 3.04 | – | – |
| dhp_anion | Vina | 1 | 10/10 | -10.20 | ASP85 OD2···N(ring)#22 2.71; SER214 N···O(carboxylate)#18 2.89; ASN104 N(aniline)#0···OD1 2.93; LYS213 NZ···O#6 2.98 | SER214:OG···O(carboxylate)#18 2.55; ASP174:OD2···N(aniline)#0 2.64; ASP174:OD2···N(ring)#7 2.70; ASP85:OD2···N(ring)#22 2.71 | – | PHE18 stacking-T 5.01; LYS213 pi-cation 4.05; ARG247 pi-cation 3.91; ARG247 pi-cation 3.71 |
| dhp_anion | Vinardo | 1 | 10/10 | -8.19 | ASP85 OD2···N(ring)#22 2.70; ASN104 N(aniline)#0···OD1 2.88; LYS213 NZ···O#6 3.04; SER214 N···O(carboxylate)#18 3.06 | SER214:OG···O(carboxylate)#18 2.52; ASP174:OD2···N(ring)#7 2.68; ASP174:OD2···N(aniline)#0 2.69; ASP85:OD2···N(ring)#22 2.70 | – | PHE18 stacking-T 4.90; LYS213 pi-cation 4.11; ARG247 pi-cation 3.91; ARG247 pi-cation 3.71 |
| smx_anion | Vina | 1 | 10/10 | -8.05 | SER52 N(aniline)#12···OG 2.99; ARG217 NH2···N(ring)#15 3.07; SER214 N···O(sulfonyl)#7 3.55; PHE179 N···O(sulfonyl)#7 3.55 | SER214:OG···N(ring)#15 2.77; POP302:O1···N(aniline)#12 2.97; SER52:OG···N(aniline)#12 2.99; ARG217:NH2···N(ring)#15 3.07 | – | – |
| smx_anion | Vinardo | 1 | 10/10 | -6.49 | ARG217 NH2···N(ring)#15 2.84; SER214 N···O(sulfonyl)#7 3.40; PHE179 N···O(sulfonyl)#7 3.58; SER52 N(aniline)#12···O 3.72 | SER214:OG···N(ring)#15 2.61; ARG217:NH2···N(ring)#15 2.84; SER52:OG···N(aniline)#12 2.92; ARG217:NH2···O(ring)#16 2.93 | – | – |
| smx_neutral | Vina | 1 | 10/10 | -8.22 | ARG217 NH2···N(ring)#15 2.99; SER214 N···O(sulfonyl)#6 3.38; PHE179 N···O(sulfonyl)#6 3.62; SER52 N(aniline)#12···O 3.76 | SER214:OG···N(ring)#15 2.72; SER52:OG···N(aniline)#12 2.88; ARG217:NH2···N(ring)#15 2.99; POP302:O1···N(aniline)#12 3.01 | – | – |
| smx_neutral | Vinardo | 1 | 10/10 | -6.70 | ARG217 NH2···N(ring)#15 2.86; SER214 N···O(sulfonyl)#6 3.39; PHE179 N···O(sulfonyl)#6 3.54; SER52 N(aniline)#12···O 3.67 | SER214:OG···N(ring)#15 2.59; ARG217:NH2···N(ring)#15 2.86; SER52:OG···N(aniline)#12 2.89; ARG217:NH2···O(ring)#16 2.94 | – | – |
| azm_neutral | Vina | 1 | 10/10 | -7.05 | LYS213 NZ···N(ring)#5 2.92; PHE179 N···O(sulfonyl)#11 3.02; SER214 N···N(sulfonamide)#9 3.27; SER52 N(amide)#3···OG 3.45 | LYS213:NZ···N(ring)#5 2.92; POP302:O6···N(ring)#5 2.93; LYS213:NZ···O(amide)#2 2.96; PHE179:N···O(sulfonyl)#11 3.02 | – | LYS213 pi-cation 3.40 |
| azm_neutral | Vinardo | 1 | 7/10 | -5.29 | LYS213 NZ···N(ring)#5 2.83; PHE179 N···O(sulfonyl)#10 3.07; SER214 N···N(sulfonamide)#9 3.15; SER52 N(amide)#3···OG 3.49 | LYS213:NZ···N(ring)#5 2.83; POP302:O6···N(ring)#5 2.93; LYS213:NZ···O(amide)#2 2.95; ARG247:NH1···O(amide)#2 2.95 | – | – |
| azm_neutral | Vinardo | 2 | 3/10 | -5.22 | MET177 N(sulfonamide)#9···O 2.85; LYS213 NZ···O(sulfonyl)#11 2.96; SER214 N···N(ring)#5 3.01; PHE180 N···N(sulfonamide)#9 3.08 | MET177:O···N(sulfonamide)#9 2.85; SER214:OG···N(amide)#3 2.87; PHE179:N···N(sulfonamide)#9 2.88; SER214:OG···N(ring)#5 2.90 | – | – |
| azm_anion | Vina | 1 | 10/10 | -6.93 | LYS213 NZ···N(ring)#5 2.90; PHE179 N···O(sulfonyl)#11 3.05; SER214 N···N(sulfonamide)#9 3.36; SER52 N(amide)#3···OG 3.46 | POP302:O6···N(ring)#5 2.88; LYS213:NZ···N(ring)#5 2.90; LYS213:NZ···O(amide)#2 3.02; PHE179:N···O(sulfonyl)#11 3.05 | – | LYS213 pi-cation 3.38 |
| azm_anion | Vinardo | 1 | 3/10 | -5.33 | MET177 N(sulfonamide)#9···O 2.95; LYS213 NZ···O(sulfonyl)#10 2.97; SER214 N···N(ring)#5 3.03; PHE180 N···N(sulfonamide)#9 3.16 | PHE179:N···N(sulfonamide)#9 2.86; SER214:OG···N(amide)#3 2.86; SER214:OG···N(ring)#5 2.92; MET177:O···N(sulfonamide)#9 2.95 | – | – |
| azm_anion | Vinardo | 2 | 7/10 | -5.14 | LYS213 NZ···N(ring)#5 2.82; PHE179 N···O(sulfonyl)#10 3.06; SER214 N···N(sulfonamide)#9 3.25; SER52 N(amide)#3···OG 3.50 | LYS213:NZ···N(ring)#5 2.82; POP302:O6···N(ring)#5 2.94; LYS213:NZ···O(amide)#2 2.97; ARG247:NH1···O(amide)#2 2.97 | – | – |
| paba_anion | Vina | 1 | 10/10 | -7.52 | SER214 N···O(carboxylate)#7 3.03; PHE179 N···O(carboxylate)#6 3.16; SER52 N(aniline)#0···O 3.58 | SER214:OG···O(carboxylate)#7 2.89; POP302:O1···N(aniline)#0 2.95; SER214:N···O(carboxylate)#7 3.03; SER52:OG···N(aniline)#0 3.09 | – | PHE18 stacking-T 4.84; LYS213 pi-cation 4.22 |
| paba_anion | Vinardo | 1 | 10/10 | -6.43 | PHE179 N···O(carboxylate)#6 2.97; SER214 N···O(carboxylate)#7 2.97; SER52 N(aniline)#0···O 3.47 | SER214:OG···O(carboxylate)#7 2.80; PHE179:N···O(carboxylate)#6 2.97; SER214:N···O(carboxylate)#7 2.97; POP302:O1···N(aniline)#0 3.04 | – | PHE18 stacking-T 4.95; LYS213 pi-cation 4.31 |

**3HS4**

| species | scoring | cluster | population | best score | PLIP H-bonds (D···A, Å) | closest N/O···N/O pairs (Å) | nearest ligand atom to Zn | π contacts |
|---|---|---|---|---|---|---|---|---|
| sfa_neutral | Vina | 1 | 10/10 | -5.78 | THR199 N···O(sulfonyl)#8 2.81; THR200 N···O(sulfonyl)#8 2.94; THR200 N(sulfonamide)#6···OG1 3.94 | THR199:N···O(sulfonyl)#8 2.81; THR200:N···O(sulfonyl)#8 2.94; HIS94:NE2···N(sulfonamide)#6 2.98; THR200:OG1···O(sulfonyl)#8 3.06 | N(sulfonamide)#6 2.50 Å | – |
| sfa_neutral | Vinardo | 1 | 10/10 | -4.99 | THR200 N···O(sulfonyl)#8 2.71; THR199 N···O(sulfonyl)#7 3.61 | THR200:N···O(sulfonyl)#8 2.71; THR199:N···O(sulfonyl)#8 2.73; THR200:OG1···O(sulfonyl)#8 2.84; HIS94:NE2···N(sulfonamide)#6 3.12 | N(sulfonamide)#6 2.30 Å | HIS94 pi-cation 4.70 |
| sfa_neutral | AD4Zn | 1 | 20/20 | -8.83 | THR199 N(sulfonamide)#6···OG1 2.46; THR200 N···O(sulfonyl)#8 2.87; THR199 N···N(sulfonamide)#6 3.63 | THR199:OG1···N(sulfonamide)#6 2.46; ASN67:OD1···N(aniline)#0 2.57; THR199:N···O(sulfonyl)#8 2.84; THR200:N···O(sulfonyl)#8 2.87 | N(sulfonamide)#6 2.32 Å | HIS94 stacking-P 3.47 |
| sfa_anion | Vina | 1 | 10/10 | -5.77 | THR199 N···O(sulfonyl)#2 2.91 | THR199:OG1···O(sulfonyl)#3 2.74; THR199:N···O(sulfonyl)#2 2.91; HIS119:ND1···N(sulfonamide)#0 2.93; HIS94:NE2···N(sulfonamide)#0 3.40 | O(sulfonyl)#3 2.34 Å | – |
| sfa_anion | Vinardo | 1 | 10/10 | -5.03 | THR200 N···O(sulfonyl)#3 2.71; THR199 N···O(sulfonyl)#2 3.54 | THR200:N···O(sulfonyl)#3 2.71; THR199:N···O(sulfonyl)#3 2.74; THR200:OG1···O(sulfonyl)#3 2.83; THR199:OG1···N(sulfonamide)#0 3.09 | N(sulfonamide)#0 2.30 Å | HIS94 pi-cation 4.66 |
| sfa_cation | Vina | 1 | 10/10 | -5.77 | THR199 N···O(sulfonyl)#3 2.83; THR200 N···O(sulfonyl)#3 2.94; THR200 N(sulfonamide)#0···OG1 3.96 | THR199:N···O(sulfonyl)#3 2.83; THR200:N···O(sulfonyl)#3 2.94; HIS94:NE2···N(sulfonamide)#0 2.96; THR200:OG1···O(sulfonyl)#3 3.05 | N(sulfonamide)#0 2.49 Å | – |
| sfa_cation | Vinardo | 1 | 10/10 | -4.97 | THR200 N···O(sulfonyl)#3 2.70; THR199 N···O(sulfonyl)#2 3.63 | THR200:N···O(sulfonyl)#3 2.70; THR199:N···O(sulfonyl)#3 2.72; THR200:OG1···O(sulfonyl)#3 2.85; HIS94:NE2···N(sulfonamide)#0 3.10 | N(sulfonamide)#0 2.30 Å | HIS94 pi-cation 4.70 |
| smx_anion | Vina | 1 | 10/10 | -7.74 | ASN62 ND2···N(ring)#15 2.92; ASN67 N(sulfonamide)#4···OD1 3.80; ASN67 ND2···N(ring)#15 4.07 | ASN62:ND2···N(ring)#15 2.92; GLN92:NE2···O(sulfonyl)#6 3.09; THR200:OG1···O(sulfonyl)#7 3.15; GLN92:NE2···N(sulfonamide)#4 3.23 | C(ring)#10 3.61 Å | HIS94 pi-cation 4.45 |
| smx_anion | Vinardo | 1 | 4/10 | -5.64 | THR199 N(aniline)#12···OG1 2.79; GLN92 NE2···N(ring)#15 2.89 | THR199:OG1···N(aniline)#12 2.79; GLN92:NE2···N(ring)#15 2.89; GLN92:NE2···N(sulfonamide)#4 2.95; HIS94:ND1···N(ring)#15 3.09 | N(aniline)#12 2.74 Å | HIS94 stacking-T 4.24; HIS94 stacking-T 4.73 |
| smx_anion | Vinardo | 2 | 2/10 | -5.56 | THR200 N(sulfonamide)#4···OG1 2.74; GLN92 NE2···O(sulfonyl)#7 3.06 | THR200:OG1···N(sulfonamide)#4 2.74; THR199:OG1···O(ring)#16 2.83; THR199:OG1···N(ring)#15 3.05; GLN92:NE2···O(sulfonyl)#7 3.06 | O(ring)#16 2.37 Å | HIS94 pi-cation 4.23 |
| smx_anion | Vinardo | 4 | 3/10 | -5.44 | THR199 N(aniline)#12···OG1 2.75; THR200 OG1···N(ring)#15 2.77 | THR200:OG1···O(ring)#16 2.72; THR199:OG1···N(aniline)#12 2.75; THR200:OG1···N(ring)#15 2.77; PRO201:O···O(ring)#16 2.80 | N(aniline)#12 2.75 Å | HIS94 stacking-T 4.59 |
| smx_neutral | Vina | 1 | 10/10 | -7.46 | ASN62 ND2···N(ring)#15 2.92; ASN67 N(sulfonamide)#4···OD1 3.95 | ASN62:ND2···N(ring)#15 2.92; GLN92:NE2···O(sulfonyl)#6 3.04; THR200:OG1···O(sulfonyl)#7 3.04; ASN62:ND2···O(ring)#16 3.39 | C(ring)#10 3.62 Å | HIS94 pi-cation 4.43 |
| smx_neutral | Vinardo | 2 | 6/10 | -5.60 | THR200 N(sulfonamide)#4···OG1 2.80; GLN92 NE2···O(sulfonyl)#6 3.01 | THR199:OG1···O(ring)#16 2.80; THR200:OG1···N(sulfonamide)#4 2.80; GLN92:NE2···O(sulfonyl)#6 3.01; THR199:OG1···N(ring)#15 3.08 | O(ring)#16 2.40 Å | HIS94 pi-cation 4.24 |
| smx_neutral | Vinardo | 3 | 3/10 | -5.57 | THR200 N···N(ring)#15 2.92; THR199 N···N(ring)#15 3.01; ASN62 N(aniline)#12···OD1 3.38; THR200 OG1···O(sulfonyl)#6 3.94 | THR200:OG1···N(sulfonamide)#4 2.70; THR200:OG1···N(ring)#15 2.76; THR199:N···O(ring)#16 2.79; THR199:OG1···O(ring)#16 2.89 | C#0 3.04 Å | HIS94 stacking-T 4.93 |
| smx_neutral | AD4Zn | 1 | 10/10 | -9.71 | GLN92 NE2···O(sulfonyl)#7 3.10 | THR199:OG1···N(ring)#15 2.78; PRO201:O···N(aniline)#12 2.78; THR199:OG1···O(ring)#16 2.87; HIS94:NE2···N(sulfonamide)#4 2.95 | N(ring)#15 2.33 Å | – |
| azm_neutral | Vina | 1 | 10/10 | -6.57 | ASN67 N(sulfonamide)#9···OD1 2.83; GLN92 NE2···O(sulfonyl)#10 2.92; ASN67 ND2···O(sulfonyl)#10 3.05; ASN62 ND2···N(sulfonamide)#9 3.17 | ASN67:OD1···N(sulfonamide)#9 2.83; GLN92:NE2···O(sulfonyl)#10 2.92; ASN67:ND2···O(sulfonyl)#10 3.05; THR199:OG1···O(amide)#2 3.10 | O(amide)#2 2.56 Å | HIS94 pi-cation 3.85 |
| azm_neutral | Vinardo | 1 | 10/10 | -5.32 | THR199 N···O(sulfonyl)#10 2.80; THR200 N···N(ring)#6 3.72; SER197 O(sulfonyl)#10···O 4.02 | THR199:OG1···O(sulfonyl)#11 2.49; THR200:OG1···N(ring)#6 2.70; THR199:N···O(sulfonyl)#10 2.80; GLN92:NE2···O(amide)#2 2.82 | O(sulfonyl)#11 2.20 Å | – |
| azm_neutral | AD4Zn | 1 | 8/10 | -9.30 | THR200 OG1···O(amide)#2 2.60; THR199 N···N(sulfonamide)#9 3.85 | THR199:OG1···N(sulfonamide)#9 2.59; THR200:OG1···O(amide)#2 2.60; THR199:N···O(sulfonyl)#10 2.85; THR200:O···O(amide)#2 3.01 | N(sulfonamide)#9 2.08 Å | – |
| azm_neutral | AD4Zn | 2 | 2/10 | -9.26 | GLN92 NE2···O(amide)#2 2.90; THR199 N···O(sulfonyl)#11 3.02; ASN67 ND2···O(amide)#2 3.86 | THR199:OG1···N(sulfonamide)#9 2.58; GLN92:NE2···O(amide)#2 2.90; HIS119:ND1···O(sulfonyl)#10 2.92; THR199:N···O(sulfonyl)#11 3.02 | N(sulfonamide)#9 2.07 Å | HIS94 pi-cation 4.17 |
| azm_anion | Vina | 1 | 10/10 | -6.50 | ASN67 N(sulfonamide)#9···OD1 2.80; ASN62 ND2···N(sulfonamide)#9 3.06; ASN67 ND2···O(sulfonyl)#11 3.14; THR199 N···O(amide)#2 3.77 | ASN67:OD1···N(sulfonamide)#9 2.80; GLN92:NE2···O(sulfonyl)#11 2.80; ASN62:ND2···N(sulfonamide)#9 3.06; THR199:OG1···O(amide)#2 3.12 | O(amide)#2 2.55 Å | HIS94 pi-cation 3.95 |
| azm_anion | Vinardo | 1 | 10/10 | -5.15 | THR199 N···O(sulfonyl)#10 2.78; THR199 O(sulfonyl)#10···OG1 3.38; THR200 N···N(ring)#6 3.64 | THR199:OG1···O(sulfonyl)#11 2.39; THR200:OG1···N(ring)#6 2.70; GLN92:NE2···O(amide)#2 2.73; THR200:OG1···N(ring)#5 2.74 | O(sulfonyl)#11 2.21 Å | – |
| paba_anion | Vina | 1 | 10/10 | -5.44 | THR199 N(aniline)#0···OG1 2.89 | THR199:OG1···N(aniline)#0 2.89; GLN92:NE2···O(carboxylate)#6 3.13; HIS119:ND1···N(aniline)#0 3.27 | N(aniline)#0 2.81 Å | HIS94 stacking-T 4.68 |
| paba_anion | Vinardo | 1 | 10/10 | -4.79 | THR199 N(aniline)#0···OG1 2.70 | THR199:OG1···N(aniline)#0 2.70; GLN92:NE2···O(carboxylate)#7 2.96; HIS119:ND1···N(aniline)#0 3.20 | N(aniline)#0 2.75 Å | HIS94 stacking-T 4.64 |

**6RL9**

| species | scoring | cluster | population | best score | PLIP H-bonds (D···A, Å) | closest N/O···N/O pairs (Å) | nearest ligand atom to Zn | π contacts |
|---|---|---|---|---|---|---|---|---|
| sfa_neutral | Vina | 1 | 10/10 | -5.90 | THR199 N···O(sulfonyl)#7 2.94; SER197 O(sulfonyl)#7···O 3.92 | THR199:OG1···O(sulfonyl)#8 2.75; THR199:N···O(sulfonyl)#7 2.94; HIS119:ND1···N(sulfonamide)#6 3.02; HIS94:NE2···N(sulfonamide)#6 3.45 | O(sulfonyl)#8 2.28 Å | – |
| sfa_neutral | Vinardo | 1 | 10/10 | -5.08 | THR199 N···O(sulfonyl)#7 2.76 | THR199:OG1···O(sulfonyl)#8 2.61; THR199:N···O(sulfonyl)#7 2.76; HIS119:ND1···N(sulfonamide)#6 2.80; HIS94:NE2···N(sulfonamide)#6 3.25 | O(sulfonyl)#8 2.35 Å | – |
| sfa_neutral | AD4Zn | 1 | 20/20 | -8.84 | THR199 N···N(sulfonamide)#6 3.69 | THR199:OG1···N(sulfonamide)#6 2.55; HIS119:ND1···O(sulfonyl)#7 2.91; THR199:N···O(sulfonyl)#8 2.98; HIS119:ND1···N(sulfonamide)#6 3.45 | N(sulfonamide)#6 2.21 Å | – |
| sfa_anion | Vina | 1 | 10/10 | -5.86 | THR199 N···O(sulfonyl)#2 2.89 | THR199:OG1···O(sulfonyl)#3 2.74; THR199:N···O(sulfonyl)#2 2.89; HIS119:ND1···N(sulfonamide)#0 3.00; HIS94:NE2···N(sulfonamide)#0 3.41 | O(sulfonyl)#3 2.31 Å | – |
| sfa_anion | Vinardo | 1 | 10/10 | -5.03 | THR199 N···O(sulfonyl)#2 2.72 | THR199:OG1···O(sulfonyl)#3 2.62; THR199:N···O(sulfonyl)#2 2.72; HIS119:ND1···N(sulfonamide)#0 2.80; HIS94:NE2···N(sulfonamide)#0 3.23 | O(sulfonyl)#3 2.32 Å | – |
| sfa_cation | Vina | 1 | 10/10 | -5.88 | THR199 N···O(sulfonyl)#2 3.00; SER197 O(sulfonyl)#2···O 3.92 | THR199:OG1···O(sulfonyl)#3 2.75; THR199:N···O(sulfonyl)#2 3.00; HIS119:ND1···N(sulfonamide)#0 3.04; HIS94:NE2···N(sulfonamide)#0 3.45 | O(sulfonyl)#3 2.27 Å | – |
| sfa_cation | Vinardo | 1 | 10/10 | -5.07 | THR199 N···O(sulfonyl)#2 2.77; THR200 N···O(sulfonyl)#3 4.10 | THR199:OG1···O(sulfonyl)#3 2.61; THR199:N···O(sulfonyl)#2 2.77; HIS119:ND1···N(sulfonamide)#0 2.80; HIS94:NE2···N(sulfonamide)#0 3.23 | O(sulfonyl)#3 2.36 Å | – |
| smx_anion | Vina | 1 | 10/10 | -7.75 | GLN92 NE2···O(sulfonyl)#7 2.94; ASN62 ND2···N(ring)#15 2.96; ASN67 N(sulfonamide)#4···OD1 3.91; ASN67 ND2···N(ring)#15 4.02 | THR200:OG1···O(sulfonyl)#6 2.92; GLN92:NE2···O(sulfonyl)#7 2.94; ASN62:ND2···N(ring)#15 2.96; GLN92:NE2···N(sulfonamide)#4 3.30 | C(ring)#13 3.65 Å | HIS94 pi-cation 4.47 |
| smx_anion | Vinardo | 1 | 10/10 | -6.04 | THR200 N···N(ring)#15 2.94; GLN92 NE2···O(sulfonyl)#6 3.04; THR199 N···N(ring)#15 3.04; HIS64 N(aniline)#12···O 4.06 | THR200:OG1···N(ring)#15 2.70; THR199:N···O(ring)#16 2.77; THR199:OG1···O(ring)#16 2.87; THR200:OG1···N(sulfonamide)#4 2.91 | C#0 3.07 Å | HIS94 stacking-T 4.93 |
| smx_neutral | Vina | 1 | 3/10 | -7.76 | GLN92 NE2···O(sulfonyl)#7 2.83; THR200 N···N(ring)#15 3.27; THR199 N···N(ring)#15 3.29; HIS64 N(aniline)#12···O 3.34 | THR200:OG1···N(ring)#15 2.82; GLN92:NE2···O(sulfonyl)#7 2.83; THR199:N···O(ring)#16 2.91; THR200:OG1···N(sulfonamide)#4 3.01 | C#0 3.38 Å | HIS94 stacking-T 4.82 |
| smx_neutral | Vina | 2 | 7/10 | -7.46 | ASN62 ND2···N(ring)#15 2.98; ASN67 N(sulfonamide)#4···OD1 4.04 | THR200:OG1···O(sulfonyl)#7 2.88; GLN92:NE2···O(sulfonyl)#6 2.92; ASN62:ND2···N(ring)#15 2.98; GLN92:NE2···N(sulfonamide)#4 3.49 | C(ring)#10 3.67 Å | HIS94 pi-cation 4.50 |
| smx_neutral | Vinardo | 1 | 4/10 | -6.10 | THR200 N···N(ring)#15 2.94; THR199 N···N(ring)#15 3.05; HIS64 N(aniline)#12···O 4.07 | THR200:OG1···N(ring)#15 2.70; THR199:N···O(ring)#16 2.74; THR199:OG1···O(ring)#16 2.87; THR200:OG1···N(sulfonamide)#4 2.88 | C#0 3.08 Å | HIS94 stacking-T 4.97 |
| smx_neutral | Vinardo | 2 | 6/10 | -5.98 | THR200 N(sulfonamide)#4···OG1 2.69; GLN92 NE2···O(sulfonyl)#6 2.84 | THR200:OG1···N(sulfonamide)#4 2.69; THR199:OG1···O(ring)#16 2.70; GLN92:NE2···O(sulfonyl)#6 2.84; THR199:OG1···N(ring)#15 3.02 | O(ring)#16 2.33 Å | HIS94 pi-cation 4.17 |
| smx_neutral | AD4Zn | 1 | 10/10 | -9.93 | THR200 N(sulfonamide)#4···OG1 2.91; GLN92 NE2···O(sulfonyl)#7 3.11; THR199 N···O(ring)#16 3.54 | THR199:OG1···N(ring)#15 2.72; PRO201:O···N(aniline)#12 2.86; THR200:OG1···N(sulfonamide)#4 2.91; THR199:OG1···O(ring)#16 3.07 | N(ring)#15 2.13 Å | – |
| azm_neutral | Vina | 1 | 10/10 | -6.69 | GLN92 NE2···N(sulfonamide)#9 2.97; ASN67 N(sulfonamide)#9···OD1 3.09; ASN67 ND2···N(sulfonamide)#9 3.12; THR200 N···N(ring)#5 3.88 | GLN92:NE2···N(sulfonamide)#9 2.97; ASN67:OD1···N(sulfonamide)#9 3.09; ASN67:ND2···N(sulfonamide)#9 3.12; THR199:OG1···N(ring)#5 3.12 | N(ring)#5 2.82 Å | HIS94 stacking-P 4.06 |
| azm_neutral | Vinardo | 1 | 3/10 | -5.30 | ASN67 N(sulfonamide)#9···OD1 2.71; ASN67 ND2···O(sulfonyl)#11 2.74; GLN92 NE2···O(sulfonyl)#11 2.80; ASN62 ND2···N(sulfonamide)#9 2.81 | ASN67:OD1···N(sulfonamide)#9 2.71; ASN67:ND2···O(sulfonyl)#11 2.74; GLN92:NE2···O(sulfonyl)#11 2.80; ASN62:ND2···N(sulfonamide)#9 2.81 | O(amide)#2 2.27 Å | HIS94 pi-cation 3.95 |
| azm_neutral | Vinardo | 2 | 7/10 | -5.26 | THR199 N···O(sulfonyl)#10 2.79; THR200 N···N(ring)#6 3.58; SER197 O(sulfonyl)#10···O 4.04; ASN67 ND2···O(amide)#2 4.08 | THR199:OG1···O(sulfonyl)#11 2.44; THR200:OG1···N(ring)#6 2.70; GLN92:NE2···O(amide)#2 2.75; THR199:N···O(sulfonyl)#10 2.79 | O(sulfonyl)#11 2.20 Å | – |
| azm_neutral | AD4Zn | 1 | 2/10 | -9.52 | THR199 N···O(sulfonyl)#10 3.03; THR200 OG1···O(amide)#2 3.07; THR200 N···N(ring)#6 3.09 | PRO201:O···O(amide)#2 2.44; THR199:OG1···N(sulfonamide)#9 2.46; THR200:OG1···N(ring)#5 2.63; THR199:N···N(ring)#6 2.95 | N(sulfonamide)#9 1.98 Å | – |
| azm_neutral | AD4Zn | 2 | 3/10 | -9.51 | GLN92 NE2···O(amide)#2 2.86; THR199 N···N(sulfonamide)#9 3.68 | THR199:OG1···N(sulfonamide)#9 2.58; THR199:N···O(sulfonyl)#10 2.79; HIS119:ND1···O(sulfonyl)#11 2.83; GLN92:NE2···O(amide)#2 2.86 | N(sulfonamide)#9 2.24 Å | HIS94 stacking-T 5.16 |
| azm_neutral | AD4Zn | 3 | 5/10 | -9.49 | THR200 N(amide)#3···OG1 2.93; THR199 N···O(sulfonyl)#11 2.98; GLN92 NE2···O(amide)#2 3.05; ASN67 ND2···O(amide)#2 3.88 | THR199:OG1···N(ring)#6 2.77; THR200:OG1···N(amide)#3 2.93; THR199:N···O(sulfonyl)#11 2.98; GLN92:NE2···O(amide)#2 3.05 | N(ring)#6 2.08 Å | HIS94 pi-cation 3.87 |
| azm_anion | Vina | 1 | 10/10 | -6.66 | GLN92 NE2···N(sulfonamide)#9 3.03; ASN67 N(sulfonamide)#9···OD1 3.15; ASN67 ND2···N(sulfonamide)#9 3.24; THR200 N···N(ring)#5 3.93 | GLN92:NE2···N(sulfonamide)#9 3.03; ASN67:OD1···N(sulfonamide)#9 3.15; THR199:OG1···N(ring)#5 3.19; THR200:OG1···N(ring)#6 3.22 | N(ring)#5 2.81 Å | HIS94 stacking-P 4.04 |
| azm_anion | Vinardo | 1 | 2/10 | -5.26 | ASN67 N(sulfonamide)#9···OD1 2.79; ASN62 ND2···N(sulfonamide)#9 2.82; ASN67 ND2···O(sulfonyl)#11 2.85; THR199 N···O(amide)#2 3.75 | THR199:OG1···O(amide)#2 2.77; GLN92:NE2···O(sulfonyl)#11 2.78; ASN67:OD1···N(sulfonamide)#9 2.79; THR200:OG1···N(ring)#5 2.81 | O(amide)#2 2.23 Å | HIS94 pi-cation 3.95 |
| azm_anion | Vinardo | 2 | 7/10 | -5.13 | THR199 N···O(sulfonyl)#10 2.70; THR199 O(sulfonyl)#10···OG1 3.16; THR200 N···N(ring)#6 3.36 | THR199:OG1···O(sulfonyl)#11 2.31; THR199:N···O(sulfonyl)#10 2.70; THR200:OG1···N(ring)#6 2.70; THR200:OG1···N(ring)#5 2.75 | O(sulfonyl)#11 2.15 Å | – |
| paba_anion | Vina | 1 | 8/10 | -5.46 | THR199 N(aniline)#0···OG1 2.93; GLN92 NE2···O(carboxylate)#6 2.96 | THR199:OG1···N(aniline)#0 2.93; GLN92:NE2···O(carboxylate)#6 2.96; HIS119:ND1···N(aniline)#0 3.24 | N(aniline)#0 2.79 Å | HIS94 stacking-T 4.98 |
| paba_anion | Vina | 2 | 2/10 | -5.45 | THR199 N···O(carboxylate)#6 3.11 | THR199:OG1···O(carboxylate)#6 3.02; THR199:N···O(carboxylate)#6 3.11; HIS119:ND1···O(carboxylate)#7 3.36 | O(carboxylate)#7 3.24 Å | – |
| paba_anion | Vinardo | 1 | 10/10 | -4.98 | GLN92 NE2···O(carboxylate)#7 2.71; THR199 N(aniline)#0···OG1 2.71 | GLN92:NE2···O(carboxylate)#7 2.71; THR199:OG1···N(aniline)#0 2.71; HIS119:ND1···N(aniline)#0 3.10 | N(aniline)#0 2.68 Å | HIS94 stacking-T 4.83 |

**3TZF**

| species | scoring | cluster | population | best score | PLIP H-bonds (D···A, Å) | closest N/O···N/O pairs (Å) | nearest ligand atom to Zn | π contacts |
|---|---|---|---|---|---|---|---|---|
| sfa_neutral | Vina | 1 | 10/10 | -6.44 | SER222 OG···N(sulfonamide)#6 3.16; SER222 N···O(sulfonyl)#8 3.28; THR62 N(aniline)#0···OG1 3.29 | SER222:N···O(sulfonyl)#7 2.91; SER222:OG···O(sulfonyl)#8 2.98; SER222:OG···N(sulfonamide)#6 3.16; HH2278:O6P···N(aniline)#0 3.16 | – | PHE190 stacking-T 5.10; PHE28 stacking-T 5.06; LYS221 pi-cation 4.99 |
| sfa_neutral | Vinardo | 1 | 10/10 | -5.35 | SER222 OG···N(sulfonamide)#6 2.91; SER222 N···O(sulfonyl)#8 3.05; THR62 N(aniline)#0···OG1 3.49 | SER222:OG···O(sulfonyl)#8 2.77; SER222:N···O(sulfonyl)#7 2.79; SER222:OG···N(sulfonamide)#6 2.91; SER222:N···O(sulfonyl)#8 3.05 | – | PHE190 stacking-T 5.16; PHE28 stacking-T 5.12; LYS221 pi-cation 5.00 |
| sfa_anion | Vina | 1 | 10/10 | -6.41 | SER222 OG···N(sulfonamide)#0 3.17; SER222 N···O(sulfonyl)#3 3.25; THR62 N(aniline)#8···OG1 3.29 | SER222:OG···O(sulfonyl)#3 2.95; SER222:N···O(sulfonyl)#2 2.99; HH2278:O6P···N(aniline)#8 3.14; SER222:OG···N(sulfonamide)#0 3.17 | – | PHE190 stacking-T 5.10; PHE28 stacking-T 5.10; LYS221 pi-cation 5.02 |
| sfa_anion | Vinardo | 1 | 10/10 | -5.27 | SER222 OG···N(sulfonamide)#0 2.93; SER222 N···O(sulfonyl)#3 3.13; THR62 N(aniline)#8···OG1 3.44 | SER222:N···O(sulfonyl)#2 2.79; SER222:OG···O(sulfonyl)#3 2.87; SER222:OG···N(sulfonamide)#0 2.93; SER222:N···O(sulfonyl)#3 3.13 | – | PHE190 stacking-T 5.17; PHE28 stacking-T 5.08; LYS221 pi-cation 5.01 |
| sfa_cation | Vina | 1 | 10/10 | -6.42 | SER222 OG···N(sulfonamide)#0 3.15; THR62 N(aniline)#8···OG1 3.26; SER222 N···O(sulfonyl)#3 3.29 | SER222:N···O(sulfonyl)#2 2.91; SER222:OG···O(sulfonyl)#3 2.98; HH2278:O6P···N(aniline)#8 3.14; SER222:OG···N(sulfonamide)#0 3.15 | – | PHE190 stacking-T 5.12; PHE28 stacking-T 5.06; LYS221 pi-cation 5.01 |
| sfa_cation | Vinardo | 1 | 10/10 | -5.24 | SER222 OG···N(sulfonamide)#0 2.93; SER222 N···O(sulfonyl)#3 3.16; THR62 N(aniline)#8···OG1 3.41 | SER222:N···O(sulfonyl)#2 2.79; SER222:OG···O(sulfonyl)#3 2.89; SER222:OG···N(sulfonamide)#0 2.93; SER222:N···O(sulfonyl)#3 3.16 | – | PHE190 stacking-T 5.18; PHE28 stacking-T 5.07; LYS221 pi-cation 5.01 |
| smx_anion | Vina | 1 | 10/10 | -7.70 | SER222 N···O(sulfonyl)#6 3.27; THR62 N(aniline)#12···OG1 3.27; SER222 OG···O(sulfonyl)#7 3.31; ARG235 NH2···N(ring)#15 3.38 | SER222:N···N(sulfonamide)#4 2.90; SER222:OG···O(sulfonyl)#6 3.08; HH2278:O6P···N(aniline)#12 3.19; THR62:O···N(aniline)#12 3.22 | – | PHE190 stacking-T 5.09; PHE28 stacking-T 5.11; LYS221 pi-cation 5.05 |
| smx_anion | Vinardo | 1 | 9/10 | -6.08 | ARG235 NH2···N(ring)#15 2.80; SER222 N···O(sulfonyl)#7 3.07; THR62 N(aniline)#12···OG1 3.43 | ARG235:NH2···N(ring)#15 2.80; SER222:OG···O(sulfonyl)#7 2.85; SER222:N···N(sulfonamide)#4 2.87; SER222:N···O(sulfonyl)#7 3.07 | – | PHE190 stacking-T 5.17; PHE28 stacking-T 5.13; LYS221 pi-cation 5.09 |
| smx_neutral | Vina | 1 | 10/10 | -7.83 | SER222 OG···N(ring)#15 2.80; GLY189 N(sulfonamide)#4···O 3.23; SER222 N···O(sulfonyl)#6 3.32; THR62 N(aniline)#12···OG1 3.37 | SER222:OG···N(ring)#15 2.80; GLN149:NE2···O(ring)#16 2.88; HH2278:O6P···N(aniline)#12 2.99; SER222:OG···O(sulfonyl)#6 3.03 | – | PHE190 stacking-T 4.99; PHE28 stacking-T 5.28; LYS221 pi-cation 5.05 |
| smx_neutral | Vinardo | 1 | 10/10 | -6.05 | SER222 OG···N(ring)#15 2.77; GLY189 N(sulfonamide)#4···O 3.07; SER222 N···O(sulfonyl)#7 3.08; THR62 N(aniline)#12···OG1 3.58 | SER222:OG···N(ring)#15 2.77; HH2278:O6P···N(aniline)#12 2.84; SER222:OG···O(sulfonyl)#7 2.89; GLN149:NE2···O(ring)#16 2.95 | – | PHE190 stacking-T 5.06; PHE28 stacking-T 5.15; LYS221 pi-cation 4.94 |
| azm_neutral | Vina | 1 | 10/10 | -5.64 | GLY189 N(amide)#3···O 2.93; SER222 OG···N(ring)#6 2.95; SER222 N···O(amide)#2 3.19; GLN226 NE2···N(sulfonamide)#9 3.25 | SER222:OG···N(ring)#5 2.85; GLY189:O···N(amide)#3 2.93; GLN149:NE2···N(ring)#6 2.94; SER222:OG···N(ring)#6 2.95 | – | – |
| azm_neutral | Vinardo | 1 | 7/10 | -4.62 | ARG235 NH2···N(ring)#6 2.80; SER222 N···O(amide)#2 3.07; ARG235 NH1···N(sulfonamide)#9 3.10 | SER222:N···N(ring)#5 2.79; ARG235:NH2···N(ring)#6 2.80; SER222:OG···O(amide)#2 2.86; ARG235:NH1···O(sulfonyl)#10 2.89 | – | – |
| azm_neutral | Vinardo | 2 | 3/10 | -4.29 | SER222 OG···N(ring)#6 2.97; GLY189 N(amide)#3···O 2.97; SER222 N···O(amide)#2 3.05; GLN149 NE2···O(sulfonyl)#11 3.83 | SER222:OG···O(amide)#2 2.79; SER222:OG···N(ring)#5 2.87; GLN149:NE2···N(ring)#6 2.91; GLY189:O···N(amide)#3 2.97 | – | – |
| azm_anion | Vina | 1 | 10/10 | -5.59 | GLY189 N(amide)#3···O 2.92; SER222 OG···N(ring)#6 2.95; SER222 N···O(amide)#2 3.16; GLN149 NE2···N(sulfonamide)#9 3.32 | SER222:OG···N(ring)#5 2.81; GLN149:NE2···N(ring)#6 2.90; GLY189:O···N(amide)#3 2.92; SER222:OG···N(ring)#6 2.95 | – | – |
| azm_anion | Vinardo | 1 | 8/10 | -4.51 | ARG235 NH2···N(ring)#6 2.80; SER222 N···O(amide)#2 3.08; ARG235 NH1···N(sulfonamide)#9 3.18 | SER222:N···N(ring)#5 2.79; ARG235:NH2···N(ring)#6 2.80; SER222:OG···O(amide)#2 2.86; ARG235:NH1···O(sulfonyl)#11 2.93 | – | – |
| azm_anion | Vinardo | 2 | 2/10 | -4.34 | SER222 N···O(amide)#2 2.98; SER222 OG···N(ring)#6 2.98; GLY189 N(amide)#3···O 3.03; GLN149 NE2···N(sulfonamide)#9 3.17 | SER222:OG···O(amide)#2 2.82; SER222:OG···N(ring)#5 2.85; GLN149:NE2···N(ring)#6 2.87; GLN226:NE2···O(sulfonyl)#11 2.90 | – | – |
| paba_anion | Vina | 1 | 10/10 | -6.37 | SER222 N···O(carboxylate)#6 3.23; THR62 N(aniline)#0···OG1 3.42 | SER222:OG···O(carboxylate)#6 2.95; SER222:N···O(carboxylate)#7 3.01; THR62:O···N(aniline)#0 3.22; SER222:N···O(carboxylate)#6 3.23 | – | PHE190 stacking-T 5.16; PHE28 stacking-T 4.99; LYS221 pi-cation 4.89 |
| paba_anion | Vinardo | 1 | 10/10 | -5.39 | SER222 N···O(carboxylate)#6 2.71; THR62 N(aniline)#0···OG1 3.52 | SER222:N···O(carboxylate)#6 2.71; SER222:OG···O(carboxylate)#7 2.82; HH2278:O6P···N(aniline)#0 3.19; SER222:OG···O(carboxylate)#6 3.31 | – | PHE190 stacking-T 5.16; PHE28 stacking-T 5.07; LYS221 pi-cation 4.94 |

## 5. Ligand efficiency and its size bias

LE = −score / heavy atoms. It has a well known size bias: docking scores grow sub-linearly with molecular size, so LE systematically favours small ligands. In this set sulfanilamide has 11 heavy atoms, p-aminobenzoate 10, acetazolamide 13, sulfamethoxazole 17 and dihydropteroate 23, so an LE comparison across them is largely a comparison of size. Measured here, over every production job of each scoring function:

| scoring | jobs | heavy atoms | r(score, heavy atoms) | r(LE, heavy atoms) | r(FQ, heavy atoms) |
|---|---|---|---|---|---|
| Vina | 33 | 10–23 | -0.87 | -0.72 | 0.45 |
| Vinardo | 33 | 10–23 | -0.76 | -0.77 | 0.12 |
| AD4Zn | 6 | 11–17 | -0.94 | -1.00 | -0.94 |

Read it like this: r(score, HA) far below zero means bigger ligands simply score better in this set, and r(LE, HA) far below zero means the per-atom measure then favours the small ones — the two halves of the size bias, measured rather than asserted. r(FQ, HA) near zero means the size correction has removed most of it. FQ = LE / LE_scale(HA), LE_scale = 0.0715 + 7.5328/HA + 25.7079/HA² − 361.4722/HA³ (Reynolds, Tounge & Bembenek 2007, doi:10.1016/j.bmcl.2007.07.071), removes most of it; but that scale was fitted to *experimental* affinities, so applied to docking scores it is a heuristic, not a calibrated quantity. **Neither LE nor FQ is used here to claim potency**; both are reported as the common size-corrected measures, and the score ranking against the controls in section 2 is the defensible comparison.

## 6. The DFT–docking link: the C3–C4–S1–N2 torsion and the strain of the docked pose

The crystal value is **102.66°** for C3–C4–S1–N2 and -77.38° for C5–C4–S1–N2 (structures/sfa_pair_xray.xyz (CIF coordinates; CIF _geom_torsion C3-C4-S1-N2 102.68(19), C5-C4-S1-N2 -77.4(2))). Which of the two a given structure reports depends on the arbitrary choice of ring carbon, so the table also gives the ring-flip-invariant `τ_fold = max(|C3–C4–S1–N2|, |C5–C4–S1–N2|)`: 90° means the S–N bond is perpendicular to the ring plane, 180° means it lies in it.

| geometry | C3–C4–S1–N2 | C5–C4–S1–N2 | τ_fold |
|---|---|---|---|
| crystal (CIF) | 102.66° | -77.38° | 102.66° |
| DFT optimum, sulfanilamide (neutral) (PCM water) | 87.16° | -90.44° | 90.44° |
| DFT optimum, sulfonamidate anion (PCM water) | 101.70° | -74.50° | 101.70° |
| DFT optimum, 4-sulfamoylanilinium cation (PCM water) | 88.95° | -89.60° | 89.60° |

The torsions below are measured on the rebuilt poses (see "How the torsion of a pose is obtained" under the table). The paper (Section 3.8, SI Table S31) measures them on the raw pose coordinates as the docking program wrote them, which differ from the values here by at most 0.18°, and takes the strain from the relaxed DFT torsion scans instead of MMFF94s.

| target | species | scoring | cluster | population | τ_fold of the top poses | MMFF94s strain, rigid (kcal/mol) | MMFF94s strain, relaxed (kcal/mol) |
|---|---|---|---|---|---|---|---|
| 7S2K | sfa_neutral | Vina | 1 | 10/10 | 119.4–168.4° | 0.36–1.06 | 0.28–1.76 |
| 7S2K | sfa_neutral | Vinardo | 1 | 10/10 | 170.1–170.5° | 1.09–1.10 | 1.79–1.80 |
| 7S2K | sfa_anion | Vina | 1 | 10/10 | 113.6–118.1° | 0.13–0.19 | 0.20–0.27 |
| 7S2K | sfa_anion | Vinardo | 1 | 10/10 | 115.3–115.3° | 0.15–0.16 | 0.22–0.22 |
| 7S2K | sfa_cation | Vina | 1 | 10/10 | 122.3–124.0° | 0.39–0.41 | 0.34–0.38 |
| 7S2K | sfa_cation | Vinardo | 1 | 10/10 | 121.9–122.1° | 0.38–0.39 | 0.33–0.33 |
| 3HS4 | sfa_neutral | Vina | 1 | 10/10 | 137.0–155.5° | 0.39–0.76 | 0.75–1.40 |
| 3HS4 | sfa_neutral | Vinardo | 1 | 10/10 | 126.8–131.9° | 0.44–0.51 | 0.39–0.51 |
| 3HS4 | sfa_neutral | AD4Zn | 1 | 20/20 | 139.9–146.6° | 0.43–0.55 | 0.84–1.07 |
| 3HS4 | sfa_anion | Vina | 1 | 10/10 | 148.0–150.4° | 0.55–0.60 | 0.83–0.87 |
| 3HS4 | sfa_anion | Vinardo | 1 | 10/10 | 127.9–132.2° | 0.85–0.94 | 0.37–0.45 |
| 3HS4 | sfa_cation | Vina | 1 | 10/10 | 138.8–145.5° | 0.59–0.72 | 0.78–0.94 |
| 3HS4 | sfa_cation | Vinardo | 1 | 10/10 | 124.7–129.7° | 0.37–0.44 | 0.39–0.51 |
| 6RL9 | sfa_neutral | Vina | 1 | 10/10 | 159.5–161.5° | 0.86–0.91 | 1.53–1.59 |
| 6RL9 | sfa_neutral | Vinardo | 1 | 10/10 | 155.6–156.7° | 0.77–0.79 | 1.41–1.44 |
| 6RL9 | sfa_neutral | AD4Zn | 1 | 20/20 | 140.6–151.8° | 0.44–0.79 | 0.86–1.25 |
| 6RL9 | sfa_anion | Vina | 1 | 10/10 | 153.5–158.1° | 0.68–0.94 | 0.91–0.96 |
| 6RL9 | sfa_anion | Vinardo | 1 | 10/10 | 152.6–153.7° | 0.65–0.68 | 0.90–0.91 |
| 6RL9 | sfa_cation | Vina | 1 | 10/10 | 160.0–161.0° | 1.05–1.10 | 1.48–1.50 |
| 6RL9 | sfa_cation | Vinardo | 1 | 10/10 | 155.2–156.7° | 0.93–0.97 | 1.34–1.38 |
| 3TZF | sfa_neutral | Vina | 1 | 10/10 | 107.7–111.8° | 0.10–0.14 | 0.10–0.15 |
| 3TZF | sfa_neutral | Vinardo | 1 | 10/10 | 104.4–129.3° | 0.08–0.49 | 0.06–0.52 |
| 3TZF | sfa_anion | Vina | 1 | 10/10 | 114.1–115.9° | 0.46–0.49 | 0.20–0.23 |
| 3TZF | sfa_anion | Vinardo | 1 | 10/10 | 113.6–114.4° | 0.45–0.46 | 0.20–0.21 |
| 3TZF | sfa_cation | Vina | 1 | 10/10 | 109.1–112.9° | 0.26–0.30 | 0.11–0.16 |
| 3TZF | sfa_cation | Vinardo | 1 | 10/10 | 110.0–110.9° | 0.27–0.28 | 0.12–0.13 |

**The link in one sentence.** Across all 280 docked top poses of the three sulfanilamide species, τ_fold is 104.4–170.5°, i.e. the S–N bond is rotated away from the perpendicular arrangement of the isolated DFT molecule (90.4° for the neutral species) and past the crystal value (102.7°) towards the ring plane, at an MMFF94s cost of at most 1.80 kcal/mol. The receptors therefore ask for a modest, affordable change in this torsion rather than an unusual conformer; the relaxed DFT torsion scans (paper Section 3.8) give the cost at the B3LYP-D3(BJ) level.

**How the torsion of a pose is obtained.** Vina writes only the atoms it moved, so each top pose was rebuilt as the full DFT molecule: the DFT geometry's rotatable torsions are set to the docked values and the result is superposed on the docked atoms. The residual of that fit is recorded for every pose (`rebuild_residual_A`, ≤ 0.0014 Å over all poses), which is what makes the torsion a measurement rather than a re-perception.

**The two MMFF94s strain estimates.** `rigid` holds the DFT bond lengths and angles fixed, as the docked pose does, and varies only the torsions: E(pose, polar-H rotors re-optimised) minus the minimum of the same rigid profile. `relaxed` minimises MMFF94s with C3–C4–S1–N2 restrained at the pose value and subtracts the MMFF94s global minimum. The polar-H rotors are re-optimised because Vina does not score hydrogen positions, so their docked orientation carries no information. Both profiles (every 5°) are in `production_summary.json` → `mmff_torsion_profiles`. MMFF94s is a cheap estimate: the profile is shallow — over sfa_anion: 1.08 kcal/mol, sfa_cation: 1.77 kcal/mol, sfa_neutral: 1.85 kcal/mol between the perpendicular minimum and the in-plane maximum (relaxed profiles) — so the ranking of poses by strain is safe, and the paper takes the absolute numbers from the relaxed DFT torsion scans (06_DFT_results/torsion_scan of the data package). `dock/results/poses_for_dft/` holds the 26 pose geometries (one per cluster representative, with and without the polar-H re-optimisation) plus a `manifest.json` giving the charge, multiplicity, torsion, MMFF strain and the reference `opt.xyz` energy for each, from which the rigid B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM single points of burst 2 were computed; the relaxed torsion scans replaced those single points.

## 7. Production against validation, for the jobs that repeat one

A production job with an ETKDG species repeats its validation run exactly and reproduces it bit for bit; a DFT-geometry job differs from its validation counterpart only in the input geometry, so the difference is the input-geometry sensitivity.

How big that input change is, DFT optimum against the ETKDG/MMFF94s conformer the validation used (matched heavy-atom skeleton):

| species | bonds: RMSD / max | angles: RMSD / max | C3–C4–S1–N2, DFT vs ETKDG |
|---|---|---|---|
| sfa_neutral | 0.0134 / 0.0247 Å | 2.05 / 7.48° | 87.2° vs 91.0° |
| sfa_anion | 0.0170 / 0.0408 Å | 2.05 / 6.16° | 101.7° vs -90.5° |
| sfa_cation | 0.0157 / 0.0321 Å | 1.96 / 6.15° | 89.0° vs -98.0° |

Both score columns are means over the seeds the two runs **share** (`n_common_seeds` in the JSON), so for the two 20-seed AD4Zn jobs the production column is the mean over seeds 1–10 and differs in the third decimal from the 20-seed mean in sections 1 and 2.

| production job | validation run | production geometry | validation score | production score | Δ | all pose scores identical |
|---|---|---|---|---|---|---|
| `T1_3TZF/smx_anion__vina__prod` | `T1_3TZF/runs/smx_anion__vina` | etkdg | -7.636 | -7.636 | 0.000 | True |
| `T1_3TZF/smx_anion__vinardo__prod` | `T1_3TZF/runs/smx_anion__vinardo` | etkdg | -5.926 | -5.926 | 0.000 | True |
| `T1_3TZF/smx_neutral__vina__prod` | `T1_3TZF/runs/smx_neutral__vina` | etkdg | -7.799 | -7.799 | 0.000 | True |
| `T3_6RL9/azm_anion__vina__prod` | `T3_6RL9/runs/azm_anion__vina` | etkdg | -6.637 | -6.637 | 0.000 | True |
| `T3_6RL9/azm_neutral__ad4zn__prod` | `T3_6RL9/runs/azm_neutral__ad4zn` | etkdg | -9.484 | -9.484 | 0.000 | True |
| `T3_6RL9/azm_neutral__vina__prod` | `T3_6RL9/runs/azm_neutral__vina` | etkdg | -6.653 | -6.653 | 0.000 | True |
| `T3_6RL9/sfa_anion__vina__prod` | `T3_6RL9/runs/sfa_anion__vina` | dft | -5.881 | -5.838 | 0.043 | – (different input geometry) |
| `T3_6RL9/sfa_neutral__ad4zn__prod` | `T3_6RL9/runs/sfa_neutral__ad4zn` | dft | -8.862 | -8.835 | 0.027 | – (different input geometry) |
| `T3_6RL9/sfa_neutral__vina__prod` | `T3_6RL9/runs/sfa_neutral__vina` | dft | -5.870 | -5.880 | -0.010 | – (different input geometry) |
| `T3_6RL9/sfa_neutral__vinardo__prod` | `T3_6RL9/runs/sfa_neutral__vinardo` | dft | -5.075 | -5.077 | -0.002 | – (different input geometry) |
| `T4_3HS4/azm_anion__vina__prod` | `T4_3HS4/runs/azm_anion__vina` | etkdg | -6.477 | -6.477 | 0.000 | True |
| `T4_3HS4/azm_neutral__ad4zn__prod` | `T4_3HS4/runs/azm_neutral__ad4zn` | etkdg | -9.257 | -9.257 | 0.000 | True |
| `T4_3HS4/azm_neutral__vina__prod` | `T4_3HS4/runs/azm_neutral__vina` | etkdg | -6.409 | -6.409 | 0.000 | True |
| `T4_3HS4/azm_neutral__vinardo__prod` | `T4_3HS4/runs/azm_neutral__vinardo` | etkdg | -5.272 | -5.272 | 0.000 | True |
| `T4_3HS4/sfa_anion__vina__prod` | `T4_3HS4/runs/sfa_anion__vina` | dft | -5.765 | -5.758 | 0.007 | – (different input geometry) |
| `T4_3HS4/sfa_neutral__ad4zn__prod` | `T4_3HS4/runs/sfa_neutral__ad4zn` | dft | -8.789 | -8.816 | -0.028 | – (different input geometry) |
| `T4_3HS4/sfa_neutral__vina__prod` | `T4_3HS4/runs/sfa_neutral__vina` | dft | -5.818 | -5.767 | 0.051 | – (different input geometry) |
| `T5_7S2K/dhp_anion__vina__prod` | `T5_7S2K/runs/dhp_anion__vina` | etkdg | -10.106 | -10.106 | 0.000 | True |
| `T5_7S2K/dhp_anion__vinardo__prod` | `T5_7S2K/runs/dhp_anion__vinardo` | etkdg | -8.090 | -8.090 | 0.000 | True |

## 7b. Input-geometry sensitivity of the AD4Zn ring orientation (diagnostic)

The two hCA II AD4Zn results moved when the input geometry changed from the ETKDG/MMFF94s conformer used in validation to the DFT structure, in opposite directions and with the zinc anchor untouched in both: 6RL9 went from bimodal (0.5 of seeds within 1 Å) to the near-native pose in every seed, and 3HS4 went from the near-native pose in every seed to the pose rotated about the Zn–N/S axis in every seed. AD4Zn separates those two orientations by hundredths of a kcal/mol (`dock/VALIDATION.md` section 3), so a small change to a rigid internal coordinate can decide which of them ranks first. The two inputs differ most in the C4–S1–N2 angle, and Vina holds bond angles rigid, so that was the hypothesis, and it was tested as a 2 × 2: each of the two geometries docked as it is and with that one angle set to the other's value, nothing else changed (`dock/geometry_variant.py`, `dock/run_geometry_diagnostic.sh`, geometries in `dock/results/diagnostics/`).

C4–S1–N2 as measured in this pass (`production_summary.json` → `csn_angle`):

| geometry | C–S–N | O–S–O | S–N |
|---|---|---|---|
| DFT B3LYP-D3(BJ)/PCM, sfa_neutral | **109.45°** | 119.46° | 1.677 Å |
| DFT B3LYP-D3(BJ)/PCM, sfa_anion | **109.83°** | 113.75° | 1.580 Å |
| DFT B3LYP-D3(BJ)/PCM, sfa_cation | **108.29°** | 120.13° | 1.664 Å |
| ETKDGv3+MMFF94s, seed1 | **101.97°** | 122.40° | 1.679 Å |
| ETKDGv3+MMFF94s, seed2 | **101.44°** | 122.63° | 1.685 Å |
| ETKDGv3+MMFF94s, seed3 | **101.44°** | 122.63° | 1.685 Å |
| crystal, 3HS4 | **106.90°** | 118.78° | 1.563 Å |
| crystal, 6RL9 | **108.79°** | 116.83° | 1.596 Å |
| crystal, 3TZF | **108.11°** | 109.25° | 1.605 Å |
| crystal, SAPH (this compound, CIF) | **109.21°** | 118.67° | 1.605 Å |

The MMFF94s value is the outlier: 101.44–101.97° against 106.90–109.21° in four crystal structures (including this compound's own) and 109.45° at B3LYP-D3(BJ)/PCM — about 7° tighter than both.

| target | input geometry | C4–S1–N2 | RMSD to the crystal SAN pose (Å) | f≤1 Å | near-native pose in the ranked list | Zn–N (Å, pass fraction) | Thr199 OG1 | clusters |
|---|---|---|---|---|---|---|---|---|
| 6RL9 | DFT (production) | 109.5° | 0.61 ± 0.04 | 1.00 | 1.00 of seeds, median rank 1.0 (best of modes 0.57 Å) | 2.20 (1.00) | 2.55 (1.00) | 20/20 |
| 6RL9 | DFT, C4–S1–N2 → 102.0° | 102.0° | 3.48 ± 0.03 | 0.00 | 1.00 of seeds, median rank 2.5 (best of modes 0.66 Å) | 2.23 (1.00) | 2.58 (1.00) | 10/10 |
| 6RL9 | ETKDG/MMFF94s, C4–S1–N2 → 109.5° | 109.5° | 0.66 ± 0.03 | 1.00 | 1.00 of seeds, median rank 1.0 (best of modes 0.63 Å) | 2.17 (1.00) | 2.50 (1.00) | 10/10 |
| 6RL9 | ETKDG/MMFF94s (validation conformer) | 102.0° | 0.97 ± 0.85 | 0.90 | 1.00 of seeds, median rank 1.0 (best of modes 0.70 Å) | 2.22 (1.00) | 2.57 (1.00) | 1/10, 9/10 |
| 3HS4 | DFT (production) | 109.5° | 3.22 ± 0.02 | 0.00 | 1.00 of seeds, median rank 2.0 (best of modes 0.56 Å) | 2.32 (0.30) | 2.46 (1.00) | 20/20 |
| 3HS4 | DFT, C4–S1–N2 → 102.0° | 102.0° | 3.31 ± 0.01 | 0.00 | 1.00 of seeds, median rank 2.0 (best of modes 0.60 Å) | 2.37 (0.00) | 2.48 (1.00) | 10/10 |
| 3HS4 | ETKDG/MMFF94s, C4–S1–N2 → 109.5° | 109.5° | 0.57 ± 0.04 | 1.00 | 1.00 of seeds, median rank 1.0 (best of modes 0.54 Å) | 2.09 (1.00) | 2.59 (1.00) | 10/10 |
| 3HS4 | ETKDG/MMFF94s (validation conformer) | 102.0° | 0.59 ± 0.03 | 1.00 | 1.00 of seeds, median rank 1.0 (best of modes 0.57 Å) | 2.20 (1.00) | 2.59 (1.00) | 10/10 |

For comparison, the validation run `T3_6RL9/sfa_neutral__ad4zn` (ETKDG conformer, a fresh one per seed, 10 seeds): RMSD 2.07 ± 1.45 Å, f≤1 Å 0.50.
For comparison, the validation run `T4_3HS4/sfa_neutral__ad4zn` (ETKDG conformer, a fresh one per seed, 10 seeds): RMSD 0.62 ± 0.03 Å, f≤1 Å 1.00.

**The hypothesis is only half right, and the table says so.** In 6RL9 the angle does decide: the DFT geometry gives the near-native pose at 109.5° and the rotated pose at 102.0°, and the MMFF94s conformer gives the near-native pose either way (9 of 10 seeds at its own 102.0°). In 3HS4 the angle changes nothing: the DFT geometry gives the rotated pose at both values and the MMFF94s conformer the near-native pose at both. So the orientation is not controlled by this one coordinate; it responds to the difference between the two rigid geometries as a whole, of which the angle is the largest part.

**What is robust, in all eight AD4Zn runs of sulfanilamide across the two receptors** (the six above and the two production runs): the near-native pose is *found* in every single seed, at rank 1 or rank 2, and the zinc chemistry is right in every seed — Zn–N(sulfonamide) 2.09–2.37 Å and the Thr199 OG1 bond 2.46–2.59 Å, against 1.93–1.94 Å and 2.82–2.86 Å in the two crystal structures. What is not robust is which of the two orientations the scoring function puts first, and the score gap between them says why:

| run | median rank of the near-native pose | score gap to the top pose (kcal/mol) |
|---|---|---|
| `T3_6RL9/sfa_neutral__ad4zn__csn102` | 2.5 | 0.187 (range 0.166–0.219) |
| `T3_6RL9/sfa_neutral__ad4zn__csn109` | 1.0 | 0.000 (range 0.000–0.000) |
| `T3_6RL9/sfa_neutral__ad4zn__etkdg` | 1.0 | 0.001 (range 0.000–0.012) |
| `T3_6RL9/sfa_neutral__ad4zn__prod` | 1.0 | 0.000 (range 0.000–0.000) |
| `T4_3HS4/sfa_neutral__ad4zn__csn102` | 2.0 | 0.156 (range 0.134–0.166) |
| `T4_3HS4/sfa_neutral__ad4zn__csn109` | 1.0 | 0.000 (range 0.000–0.000) |
| `T4_3HS4/sfa_neutral__ad4zn__etkdg` | 1.0 | 0.000 (range 0.000–0.000) |
| `T4_3HS4/sfa_neutral__ad4zn__prod` | 2.0 | 0.027 (range 0.008–0.048) |

**How the paper reports it (Section 3.6).** Statements about the zinc anchor and the Thr199 hydrogen bonds are supported by every run. The orientation of the aniline ring is stated with this sensitivity: the near-native pose trails the top pose by 0.000–0.219 kcal/mol in these runs, and which orientation comes first changes with the input geometry. The paper reports the DFT geometry, which the rest of the paper uses and which agrees with the crystal structures on the angle; with it the compound's own co-crystal structure, 6RL9, is reproduced at 0.61 ± 0.04 Å in 20 of 20 seeds, while the 3HS4 cross-dock lands on the rotated orientation with the same zinc anchor.

## 8. What these numbers do and do not support

- **7S2K (Sul2)** is the only target validated with both scoring functions (0.51 / 0.54 Å, 10/10 seeds within 1 Å), so its poses and contacts are discussed quantitatively.
- **3HS4 (hCA II)** was validated as a cross-docking target for sulfanilamide with AD4Zn (0.62 ± 0.03 Å, 10/10 within 1 Å, Zn–N 2.18 Å and the Thr199 bond in 10/10) — but that validation used the ETKDG/MMFF94s conformer. With the DFT geometry of production the same protocol gives 3.22 ± 0.02 Å: the sulfonamide still points at the zinc (Zn–N 2.32 Å on average, within the 2.3 Å cutoff in only 0.30 of the seeds, against 1.94 Å in the crystal) and the Thr199 OG1 bond holds in 1.00 of them, while the aniline ring takes the rotated orientation (section 7b). Acetazolamide's own pose is reproduced by no scoring function. So 3HS4 supports the zinc-binding chemistry and the score ranking; it does not, on its own, fix the ring orientation.
- **6RL9 (hCA II)**, the compound's own co-crystal structure, failed the strict gate in validation (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal). With the DFT geometry of production it gives 0.61 ± 0.04 Å with f≤1 Å 1.00 over 20 seeds and the zinc chemistry right in every seed. That is a better result than the gate was given, and section 7b shows it is not robust to the input geometry, so the paper reports 6RL9 as partially validated, with the sensitivity and the validation figure (Figure 8) beside it, not as a passed gate.
- **3TZF (DHPS)** passes the conventional 2.0 Å criterion and not the 1.0 Å gate (sulfamethoxazole anion 1.92 ± 0.04 Å), so its poses are qualitative.
- **1AJ0 was excluded from production docking**: in validation sulfanilamide landed 11.6 Å away in a surface cleft in every seed. Nothing in this file comes from 1AJ0, which served only for the structural comparison with 3TZF.
- **Plain Vina and Vinardo have no metal term**, and it shows in the zinc site. Over the 12 hCA II jobs of the three sulfanilamide species scored with them: 6RL9 Zn–N 2.40–2.76 Å, Thr199 OG1 4.35–4.72 Å, both checks passing in 0.00–0.00 of seeds; 3HS4 Zn–N 2.30–2.74 Å, Thr199 OG1 3.10–4.71 Å, both checks passing in 0.00–1.00 of seeds — against Zn–N 1.93–1.94 Å and Thr199 OG1 2.82–2.86 Å in the two crystal structures. What those poses put nearest the zinc is O(sulfonyl) in 7 of 12 cluster representatives, N(sulfonamide) in 5 of 12 cluster representatives (section 4), and the split is by receptor: in 6RL9 all six give the zinc a sulfonyl oxygen at 2.27–2.36 Å, in 3HS4 five of six give it the sulfonamide nitrogen but at 2.30–2.50 Å, against 1.93–1.94 Å in the crystals. So these poses find the right pocket with a loose or wrong anchor, and only the Vinardo poses in 3HS4 (Zn–N exactly at the 2.30 Å cutoff) pass the zinc-site checks at all. Their hCA II scores are a scoring-function comparison, not a binding mode; AD4Zn is the scoring function that reproduces the zinc chemistry.
- Scores of 10–23 heavy-atom ligands lie within about 1 kcal/mol of each other. Only the ranking against the controls is interpreted, and section 5 says why LE does not rescue the comparison.

**Limitations of the method (paper Section 4.3).**

- **Docking scores are not binding free energies.** Vina, Vinardo and AutoDock4Zn are empirical scoring functions whose reported accuracy is of the order of 1–2 kcal/mol on benchmark sets; no score here is converted to a Ki, a Kd or an affinity, and no compound is called potent or a drug candidate. Only the ordering against the co-crystal ligand and the controls is read, and section 5 shows that the ordering itself carries a size bias.
- **The receptor is rigid.** Every run uses one crystal conformation with fixed side chains: no induced fit, no loop or side-chain rearrangement, no receptor ensemble. The ligand is flexible in its torsions only — bond lengths, bond angles and rigid fragments are held at the input geometry, which is why the input geometry changes the AD4Zn ring orientation (section 7b).
- **Water and entropy are implicit at best.** No explicit solvation model, no conformational entropy, no protonation change on binding, and all crystallographic waters were removed except the metal coordination sphere (3 in 7S2K, 3 in 3TZF, each within 2.4 Å of the Mg). The zinc chemistry is modelled by the AutoDock4Zn TZ pseudo-atom, a parameterised stand-in for the Zn–N bond, not by quantum mechanics.
- **No molecular dynamics and no free-energy calculation** (MM/PBSA, MM/GBSA, alchemical) was run: every pose is a single static structure.
- **No experimental data on this compound.** There is no measured affinity, inhibition or bioassay for 4-sulfamoylanilinium dihydrogen phosphate here; the validation is structural (redocking and cross-docking against deposited complexes) and the literature ordering acetazolamide before sulfanilamide at hCA II is the only external anchor.
- **Protonation and tautomers are assumptions.** Species are fixed at pH 7.4 by pKa argument (`dock/prep_ligands.py`), the His64 tautomer of hCA II is an explicit assumption (`dock/targets.py: his_assumed`; PDB2PQR picks HIE for 6RL9 but HID for 3HS4, and the Vina baselines moved by ≤0.03 kcal/mol when it was flipped; the His64 test was run with Vina only), and the ion pair of the crystal is never docked because it dissociates at pH 7.4.
- **Pose physics.** PoseBusters (`dock` configuration, 72 production jobs, 740 top poses, 22 checks each) passes on 99.3% of the poses, lowest per job 0.50. The only check that ever fails is `minimum_distance_to_protein`, in 5 of the 10 poses of `T3_6RL9/azm_anion__vinardo__prod` — poses whose closest receptor atom is the catalytic zinc, so what the generic van der Waals test counts as a clash is the metal-coordination distance itself. That tests the physical sanity of a pose, not its correctness.

## 9. Files

| file | what |
|---|---|
| `dock/results/production_summary.json` | every number in this report, per job: scores, clusters, RMSDs, site checks, ligand efficiency, torsions, strain, per-seed rows |
| `dock/results/production_contacts.json` | PLIP / ProLIF / geometric contacts of every cluster representative |
| `dock/results/poses/<target>/<species>__<scoring>__top_pose_per_seed.sdf` | the top pose of every seed, with score, seed and cluster as SDF properties |
| `dock/results/poses/<target>/<species>__<scoring>__best_seedNN_all_modes.pdbqt` | all poses of the best-scoring seed, as Vina wrote them |
| `dock/results/poses/<target>/<species>__<scoring>__{csn102,csn109,etkdg}__*` | the same two files for the diagnostic runs of section 7b, kept apart from the production poses by their suffix |
| `dock/results/poses_for_dft/*.xyz` + `manifest.json` | cluster-representative geometries of the rigid DFT single points of burst 2 (replaced by the relaxed torsion scans) |
| `dock/<target>/runs/<species>__<scoring>__prod/seedNN/` | the raw runs: input conformer, Vina log, all poses, `results.json` with the exact settings |
| `dock/results/production.log` | the run log with per-job timings |
| `dock/results/diagnostics/` + `dock/run_geometry_diagnostic.sh`, `dock/geometry_variant.py` | the input-geometry diagnostic of section 7b: the modified geometries, its log and its runs (`dock/T{3_6RL9,4_3HS4}/runs/sfa_neutral__ad4zn__{csn102,csn109,etkdg}/`) |
| `dock/VALIDATION.md`, `dock/analysis/{summary,gate_G2}.json` | the validation this production run rests on |

Receptors used (sha256 of the file each run read, checked against `dock/T*/prep_report_recheck.json`):

| target | receptor.pdbqt sha256 | matches the prepared receptor | AD4Zn maps |
|---|---|---|---|
| 7S2K | `81531f8ef481089c…` | True | – |
| 3HS4 | `5e9ade1962a707f5…` | True | 9 maps, TZ pseudo-atom 0.07 Å from the crystal sulfonamide N |
| 6RL9 | `f335c7365cedd676…` | True | 9 maps, TZ pseudo-atom 0.17 Å from the crystal sulfonamide N |
| 3TZF | `03a431de350ed509…` | True | – |

