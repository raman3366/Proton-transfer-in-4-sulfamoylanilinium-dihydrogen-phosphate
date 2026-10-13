# MD campaign: one line per run

Written by `md/analysis/analyse_md.py` after `md/analysis/SPEC.md`, version 1.3 (6 Oct 2026). Nothing is averaged over runs. All columns with full precision: `campaign_table.csv`; every number: `campaign.json` and the `summary.json` of each run.

**Zinc model `bonded_2e`**: 9 runs found, 9 complete (X 3, N 3, R 3). Rule 1: 0 of 3 crystal-start runs left the crystal pose (not failed). Converted: X 0, N 0, R 3. Outcome: `rotated_converts` – rule 3: the rotated starts converted in two or three runs and the crystal and near-native starts stayed; the rotated orientation is not stable outside the scoring function.

**Zinc model `bonded_2e_n2fit`**: 1 runs found, 1 complete (X 0, N 0, R 1). Rule 1: 0 of 0 crystal-start runs left the crystal pose (not failed). Converted: X 0, N 0, R 0. Outcome (NOT FINAL: the nine runs are not all present and complete): `mixed` – rule 3: the runs disagree; they are reported run by run and nothing is averaged over them.

**Zinc model `bonded_2s`**: 1 runs found, 1 complete (X 0, N 0, R 1). Rule 1: 0 of 0 crystal-start runs left the crystal pose (not failed). Converted: X 0, N 0, R 1. Outcome (NOT FINAL: the nine runs are not all present and complete): `mixed` – rule 3: the runs disagree; they are reported run by run and nothing is averaged over them.

| Run | ns seen / planned | Verdict | Start state, last 20 ns (%) | Opposite state, last 20 ns (%) | First passage (ns) | Transitions | Ligand RMSD, mean ± SD (Å) | Ligand RMSD, last half (Å) | Left crystal (rule 1) | N2–H···OG1 Thr199 (% of frames) | Zn–N2 (Å) | τ_fold, mean ± SD (°) | C4–S–N2, mean ± SD (°) | Fifth ligand (% of frames, 100 ps) | CA RMSD (Å) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bonded_2e/X1 | 100.00 / 100 | STAYED | 99.8 | 0.0 | – | 0 | 0.56 ± 0.23 | 0.58 | no | 99.9 | 1.926 | 132.4 ± 13.6 | 109.5 ± 2.7 | 0.3 | 0.88 |
| bonded_2e/X2 | 100.00 / 100 | STAYED | 99.7 | 0.0 | – | 0 | 0.54 ± 0.22 | 0.53 | no | 99.9 | 1.928 | 130.9 ± 12.9 | 109.4 ± 2.7 | 0.1 | 0.91 |
| bonded_2e/X3 | 100.00 / 100 | STAYED | 98.8 | 0.0 | – | 0 | 0.54 ± 0.21 | 0.56 | no | 99.8 | 1.928 | 130.3 ± 11.8 | 109.3 ± 2.6 | 0.1 | 0.96 |
| bonded_2e/N1 | 100.00 / 100 | STAYED | 99.6 | 0.0 | – | 0 | 0.53 ± 0.22 | 0.53 | no | 99.8 | 1.926 | 131.5 ± 13.0 | 109.3 ± 2.6 | 0.0 | 0.87 |
| bonded_2e/N2 | 100.00 / 100 | STAYED | 99.9 | 0.0 | – | 0 | 0.54 ± 0.24 | 0.53 | no | 99.9 | 1.926 | 132.9 ± 14.4 | 109.5 ± 2.7 | 0.2 | 0.87 |
| bonded_2e/N3 | 100.00 / 100 | STAYED | 99.4 | 0.0 | – | 0 | 0.61 ± 0.27 | 0.67 | no | 99.8 | 1.927 | 135.0 ± 17.3 | 109.4 ± 2.7 | 0.1 | 0.90 |
| bonded_2e/R1 | 100.00 / 100 | CONVERTED | 0.0 | 99.7 | 0.01 | 1 | 0.54 ± 0.21 | 0.50 | no | 99.9 | 1.926 | 129.6 ± 11.2 | 109.4 ± 2.6 | 0.1 | 1.01 |
| bonded_2e/R2 | 100.00 / 100 | CONVERTED | 0.0 | 99.7 | 0.01 | 1 | 0.55 ± 0.24 | 0.51 | no | 99.8 | 1.925 | 131.7 ± 14.1 | 109.4 ± 2.6 | 0.1 | 0.95 |
| bonded_2e/R3 | 100.00 / 100 | CONVERTED | 0.0 | 100.0 | 0.01 | 1 | 0.53 ± 0.21 | 0.53 | no | 99.9 | 1.927 | 130.4 ± 12.0 | 109.3 ± 2.6 | 0.0 | 0.91 |
| bonded_2e_n2fit/R1 | 100.00 / 100 | MIXED | 9.2 | 74.4 | 0.01 | 1 | 0.99 ± 0.54 | 1.08 | no | 98.3 | 1.924 | 149.1 ± 20.4 | 110.6 ± 2.7 | 1.8 | 0.98 |
| bonded_2s/R1 | 100.00 / 100 | CONVERTED | 0.0 | 99.5 | 0.01 | 1 | 0.80 ± 0.38 | 0.85 | no | 99.9 | 1.951 | 130.7 ± 18.4 | 110.2 ± 2.6 | 2.6 | 0.90 |

Folders with a trajectory that are not production runs (by their name, or a validation run by its manifest); left out, the validation run X0 is analysed on its own:

- `/root/cqc-sfa-md/runs/bonded_2e/X0`
- `/root/cqc-sfa-md/runs/bonded_2s/X0`
- `/root/cqc-sfa-md/runs/preflight_a`
- `/root/cqc-sfa-md/runs/preflight_b`
