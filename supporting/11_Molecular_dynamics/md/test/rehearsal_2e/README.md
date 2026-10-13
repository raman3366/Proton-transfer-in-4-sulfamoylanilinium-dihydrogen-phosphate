# md/test/rehearsal_2e/: files of the rehearsal of 6 Oct 2026

**A REHEARSAL, NOT A RESULT FOR THE PAPER.** Three runs of 5 ns of the bonded zinc model with parameter set 2e and
terms at the ligand nitrogen N2 that gate G-R has not validated. They test the pipeline. The account is
`md/test/REHEARSAL_2e.md`; every number there is in a file of this folder or of `md/topology/bonded_2e/`.
Nothing here was committed.

| Path | What |
|---|---|
| `look/TX0.json`, `TR0.json`, `TN0.json` | per run: the driver's records (minimisation, every stage, temperature, density, box, hydration check, speed) and the numbers of the primary analysis (gate G-V block, orientation, configuration at N2, statistics, the series stage by stage and nanosecond by nanosecond); made by `scripts/look_run.py` |
| `look/min_site_<run>.json` | zinc site and ligand as built, after the minimisation on the GPU and at the end of every step; the amide hydrogens next to the zinc histidines; made by `scripts/min_site.py` |
| `look/nearest_oxygens.json` | which water oxygen and which protein oxygen is nearest to the zinc, by name |
| `look/qm_speed.json` | seconds per optimisation step of the quantum-chemistry task before, during and after the GPU tasks |
| `look/n2_frcmod_test.json`, `n2_test.frcmod` | test of the builder with a SYNTHETIC file of terms at N2 (built in scratch; not a model) |
| `look/set_2s_test.json` | test of `--set 2s` on a relabelled copy of the 2e folder (built in scratch; not a Seminario set, not a model) |
| `look/analysis_v13_against_v12.json` | the first analysis program of specification 1.3 (another agent's work, as it stood at 11:32 AM IST) on the three runs, against the output of version 1.2 that the account uses |
| `runs/<run>/` | `manifest.json`, the records `rec_*.json`, the state logs `log_*.csv`, `run_md.log`, `final_solute.pdb`, the marker `REHEARSAL_NOT_A_RESULT.txt`; `analysis/` (primary program: `summary.json`, the four series, the RMSF table, five figures); `analysis_recheck/` (second program) |
| `comparison/` | the report of `recheck_md.py --compare-campaign` on the three runs, and the campaign files of both programs (the three runs are left out of the campaign table by both, as their names ask) |
| `push/` | log and checksum list of the copy to the machine |
| `scripts/` | `pull_run.sh`, `check_sums.py`, `look_run.py`, `min_site.py` |
| `logs/` | console output of the analysis programs and of the local checks; what was read on the machine (queue, quantum-chemistry log) |
| `build_topology_6oct.diff` | `md/build_topology.py` before and after this work |

Not here: the trajectories and saved states (124.5 MB per run). They are on the rented machine in
`/root/cqc-sfa-md/runs_test/bonded_2e/` (377 MB with the checkpoints) and, until the session's scratch folder is
cleared, on the Pi in `…/scratchpad/rehearsal/runs_test/bonded_2e/`. The watchdog does not pull `runs_test/`.
