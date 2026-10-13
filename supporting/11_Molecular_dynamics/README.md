# 11 Molecular dynamics

Molecular dynamics of the complex of human carbonic anhydrase II with sulfanilamide (PDB 6RL9; the
sulfonamidate anion on the zinc), reported in Section S16 of the supplementary material and in Sections 2.5 and 3.6 of
the main text. `md/` here is the folder `md/` of the repository, file for file (what is left out, and which copies
have an address masked, is said below): a source path `md/...` in 05_Tables/facts.md or in one of the records is
`11_Molecular_dynamics/md/...` in this package.

## What was run

- 9 runs of 100 ns with one force field (ff14SB for the protein, TIP3P water, GAFF2 2.2.30 for
  the ligand with its terms at sulfur repaired against the DFT of the paper, and a bonded model of the zinc site, set
  2e), 3 from each of three starts: the crystal pose (X), the docked near-native pose (N) and the docked
  ring-rotated pose (R).
- 2 sensitivity runs of 100 ns from start R: one with the second set of zinc terms (2s), one with candidate terms
  at the sulfonamide nitrogen. They enter no count of the rules.
- 2 validation runs of 5 ns from the crystal pose, one for each zinc set.
- OpenMM 8.6.1 on one rented GPU (NVIDIA L40S), 6 to 8 Oct 2026; 1,110 ns of production in
  all. The machine cost US$37.97 for 43.60 h of running by the project's ledger and was deleted after the last run
  (`md/PLAN.md`, entry of 8 Oct 2026, 8:30 AM IST).

## What the runs show, by the rules fixed before they were made

The rules are in `md/PLAN.md`, section 5; every number is in `md/RESULTS.md`, section 1, and in 05_Tables/facts.md.

- None of the three crystal-start runs left the crystal pose (0 of 3 by rule 1), and the
  6 runs from starts X and N stayed.
- The 3 runs from start R left the rotated pose during the equilibration, 140 to
  360 ps after the release of the ligand, while the protein backbone was still held: the configuration at the
  sulfonamide nitrogen changed, and the hydrogen bond to Thr199 was kept. By the rule fixed beforehand the times are
  160 to 540 ps. The production of these three runs holds 0 rotated-like frames of
  30,000.
- Sensitivity runs: with zinc set 2s the run converted; with the candidate terms at the nitrogen the run is MIXED,
  crystal-like in 74.35 % of the frames of its last 20 ns.

## Limits, which belong to every statement above

- The force-field terms at the sulfonamide nitrogen failed two quantum-chemical tests (gates G-R and G-R2; `md/GR.md`,
  `md/N2.md`) and were used as built. One of their errors has the sign of the outcome: the force field puts the
  rotated placement 3.38 kcal mol⁻¹ too high against the crystal-like one.
- The changes happened with the protein backbone restrained. There are three runs per start, of 100 ns, with one
  force field.
- Zinc set 2s failed its validation narrowly, and the candidate terms at the nitrogen were not accepted by their
  test; "stayed" and "no return" therefore belong to the terms as built.
- No rate, no barrier and no free energy is given, and no end-point free energy (MM-GBSA, MM-PBSA) was computed.
  The simulations replace no docking number.
- The simulations, the quantum-chemical relaxation of the bare zinc-site model (`md/qm/gr/relax/`: started from the
  rotated placement it ends 0.12 Å from the end point reached from the crystal-like one) and the docking with
  four crystal waters (07_Docking, `WATER4.md`) point the same way, in line with the crystal structure. None of the
  three decides.

## What is where

| Path | What |
|---|---|
| `md/RESULTS.md` | The results record of the campaign: the outcome by the rules, every run, the validation, the conversions run by run, the zinc site, the two sensitivity runs, what the result rests on and its limits, the health of the runs, the agreement of the two analysis programs, and where every file is. It is written by `md/production/assemble/make_results_md.py` from the records; no result is typed. |
| `md/results/series/` | Rule 6 of the plan: the time series of the 11 runs of 100 ns, one gzip-compressed CSV file per run and one line per 10 ps (time, d1, d2, ligand RMSD, state, the configuration at the nitrogen, the hydrogen bond to Thr199), with `README.md`, `SOURCES.json` and `SHA256SUMS`. |
| `md/results/final_frames/` | The last solute frame of each of the 11 runs (gzip-compressed PDB), with `README.md`, `SOURCES.json` and `SHA256SUMS`. |
| `md/results/runs/<set>/<run>/` | Per run: `summary.json` of both analysis programs and the small series of the first (equilibration, all-atom frames, Cα fluctuations). |
| `md/results/campaign/`, `gate/`, `conversion_times.json`, `run_health.json`, `run_events.json`, `site_checks_*.json`, `tables.md` | The campaign tables of both programs and the report of their comparison; the gate of the analysis stage; the conversion times; the driver's records and checksums; event counts; checks on the raw frames; the long tables of the results record. |
| `md/results/independent/` | A blind third reading of the trajectories with its own reader (`third_reading/`) and a physical review (`physical_review/REVIEW.md` with its scripts and tables), kept as its author wrote it; `md/RESULTS.md`, section 7, says what was taken from it and what was not. |
| `md/results/fig_orientation.*`, `fig_ligand_rmsd.*` | The two figures of the results record. The figures of the supplementary material were redrawn from the same series (04_Figures/SI). |
| `md/PLAN.md` | The plan: the question, the design, the stages with their gates, the rules for reading the outcome, and every later change as a dated entry. |
| `md/SYSTEM.md`, `md/LIGAND.md`, `md/ZINC.md`, `md/TOPOLOGY.md`, `md/RUN.md`, `md/PRODUCTION.md` | The method records: the system and its protonation audit, the ligand terms, the zinc model, the topologies, the driver and its protocol, the production. |
| `md/GR.md`, `md/N2.md` | The two quantum-chemical tests of the terms at the sulfonamide nitrogen, with what failed and in which direction. |
| `md/system/`, `md/ligand/`, `md/topology/`, `md/protocol.json` | Inputs and parameters: the prepared protein, the three starts, the crystal water and the audits (`system/`); the ligand parameters with their checks against DFT (`ligand/`); the tleap inputs, added terms, build records and checksums of every built system (`topology/`); the protocol that the driver read. |
| `md/qm/` | The quantum-chemical records: the zinc-site models and their charges (`mcpb/`, `results/`), gate G-R (`gr/`: geometries, single points, the gate, the two relaxations in `gr/relax/`) and gate G-R2 (`n2/`), with the scripts. |
| `md/*.py`, `md/analysis/`, `md/production/`, `md/gpu/`, `md/test/`, `md/env/` | The scripts: system and topology builders, the driver `run_md.py`, the two analysis programs with their specification and test reports, the programs that collect the results, the queue and task files of the rented machine, the tests and the rehearsal, and the software environment. |
| `trajectories_not_included/` | 13 run manifests (`<set>_<run>_manifest.json`, written by the driver: protocol, seeds, every step, and the sha256 and size of every file it wrote) and `trajectories.csv`, the list of the 602 trajectory files (DCD) of these runs with size and sha256. |

## Not included

The trajectories are not in the package and not in git: 13.4 GB in the 1284 files that the driver
checksums (16 run folders with the test folders), of which the DCD files of the 13 runs listed in
`trajectories_not_included/trajectories.csv` are 11.2 GB. They are kept on the Pi under
`/root/cqc-sfa-md/runs/`. The guard's final check of 8 Oct 2026, 4:12 AM IST compared 3128 files on the machine with the copies on the Pi by sha256 (3128 the same) before the machine was stopped. Each manifest here has the sha256 that
`md/results/run_health.json` records for it. Also not included: the built topology and coordinate files of the
systems (`system.prmtop`, `system.inpcrd`; rebuilt by `md/build_topology.py`, their sha256 in each manifest under
`inputs`), and the files of this table.

| Left out of `md/` | Files | Size |
|---|---:|---:|
| Python bytecode | 26 | 1.6 MB |
| formatted checkpoint of a quantum-chemical run | 1 | 310.1 kB |
| log or scratch folder of the simulation work that is not in git | 18 | 32.6 kB |

## How the records are to be read

The Markdown and text files under `md/` are the dated records of the work and are kept as written. A record states
what was known at its date; where a later result changed a statement, a dated note says so in the record (for
example in `md/GR.md` and `md/qm/gr/relax/README.md` on the relaxation that finished after the record was written).
The check of the package for unfinished text (VERIFICATION.md) reads these records with the same patterns as every
other document. 43 occurrences in 16 of them were read one by one and are statements of a record, not open
items of the package: the name of the queue state of a task that has not started, the provisional charges of a
self-test, a test set whose results were never opened because the verdict was fixed without it, and the like. One more
occurrence than was read stops the package.
Sanitised copies, with a network address masked: `md/qm/n2/collected_fit.json`, `md/qm/n2_collect.py`, `md/test/rehearsal_2e/scripts/pull_run.sh`.
