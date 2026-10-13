# MD campaign: one line per run

Written by `md/analysis/analyse_md.py` after `md/analysis/SPEC.md`, version 1.2 (6 Oct 2026). Nothing is averaged over runs. All columns with full precision: `campaign_table.csv`; every number: `campaign.json` and the `summary.json` of each run.

| Run | ns seen / planned | Verdict | Start state, last 20 ns (%) | Opposite state, last 20 ns (%) | First passage (ns) | Transitions | Ligand RMSD, mean ± SD (Å) | Ligand RMSD, last half (Å) | Left crystal (rule 1) | N2–H···OG1 Thr199 (% of frames) | Zn–N2 (Å) | τ_fold, mean ± SD (°) | C4–S–N2, mean ± SD (°) | Fifth ligand (% of frames, 100 ps) | CA RMSD (Å) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

Folders with a trajectory that are not production runs (by their name, or a validation run by its manifest); left out, the validation run X0 is analysed on its own:

- `/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/rehearsal/runs_test/bonded_2e/TN0`
- `/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/rehearsal/runs_test/bonded_2e/TR0`
- `/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/rehearsal/runs_test/bonded_2e/TX0`
