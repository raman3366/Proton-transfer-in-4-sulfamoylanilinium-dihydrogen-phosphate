# md/: molecular dynamics of the hCA II–sulfanilamide complex

Work decided by the user on 5 Oct 2026 (root `PLAN.md`, row "Molecular dynamics (5 Oct 2026, night)"). Why it is
done, what the scoping found and what it left unverified: `paper/audit/research/md_scoping_2026-10-05.md`.

## PLAN.md

The plan, and the file that binds all work in this folder. It fixes, before any production run: the question (does
the ring-rotated orientation at the zinc persist or relax to the crystal one), the design (hCA II 6RL9 only, three
starts, three runs of 100 ns each, bonded zinc model, no MM-GBSA), the stages S0 to S8 with their gates, and the rules
for reading the outcome. A gate that fails is reported as failed. A change to the plan is a dated line at its end.

## Sub-folders

Each sub-folder belongs to one stage of the plan and is filled as that stage is done.

| Folder | Stage | What it will hold |
|---|---|---|
| `env/` | S0 | the software environment of the set-up: packages with versions (`md_env.yml`), the exact file list (`md_env_explicit.txt`) and how to rebuild it |
| `ligand/` | S1 | sulfanilamide and its sulfonamidate anion: the GAFF2 terms checked against the paper's DFT geometries and torsion scans (`molecules/sfa/results/`), the refitted terms at sulfur, the charges, the result of gate G-L |
| `qm/` | S2 | the zinc site: the small model (zinc, three imidazoles, ligand) with its optimised geometry and Hessian, the large capped model for the electrostatic potential, the PySCF inputs and small outputs, the bonded terms and RESP charges derived from them, the result of gate G-Z |
| `system/` | S3 | the system builder: protonation audit, crystal water, the three starts, solvation, topology and the equilibration protocol |

Folders for the later stages (validation run, production, analysis) are added when their stage begins. The docking
sensitivity test of stage S3b is docking work and has its files under `dock/`.

## What goes where

- In git, in this folder: scripts, small inputs, parameters, final frames and result tables, each file under about 5 MB.
- Not in git: checkpoints, trajectories and quantum-chemical scratch files. They live on the Pi in `/root/cqc-sfa-md/`,
  one sub-folder per folder here (for example `/root/cqc-sfa-md/system/`), and go to the R2 archive of the project at
  the next sync.
- One script per step, written so that a second run gives the same files; seeds and software versions are recorded
  in the outputs.
- Nothing here changes the paper. Results are shown to the user first (`PLAN.md`, section 7).

## Added on 11 Oct 2026: the campaign is finished

The text above is the description of this folder as it stood when the work began (5 and 6 Oct 2026) and is kept as
written. What has happened since:

- **All stages of the plan up to S6 are done, and S7 (the paper) is being closed.** The systems were built and
  tested on the Pi; the nine runs of 100 ns, the two sensitivity runs of 100 ns and the two validation runs of
  5 ns ran from 6 to 8 Oct 2026 on one rented GPU machine; the two analysis programs agree on every value
  (`md/PLAN.md`, entries of 8 Oct 2026).
- **The rented machine was deleted on 8 Oct 2026** at 8:20 AM IST, at the user's word, after every file had been
  pulled to the Pi and verified. Cost US$37.97 against the cap of US$60. Nothing of this project runs in the
  account (`md/PLAN.md`, entry of 8 Oct 2026, 8:30 AM IST; `ops/README.md`).
- **The results are in `md/RESULTS.md`** (outcome by the rules fixed before the runs, every number with its file,
  the limits, the corrections after two audits and a physical review) **and in `md/results/`**: per-run summaries
  of both programs, `tables.md`, the series of every run of 100 ns (`series/`, one line per 10 ps), the last solute
  frame of every run (`final_frames/`), the health and event records, and the two readings made with other code
  (`independent/`). `md/RESULTS.md` is written by `md/production/assemble/make_results_md.py` and is not edited
  by hand.
- **In the paper** the simulations are Section S16 of the supplementary material (Tables S38 to S40, Figs. S9 and
  S10), with a paragraph each in Sections 2.5 and 3.6 of the main text; the two docking tests on hCA II
  (`dock/GLU106.md`, `dock/WATER4.md`) are Section S15 (Tables S36 and S37). The user approved their entry on
  8 Oct 2026, so the last line of "What goes where" above ("Nothing here changes the paper") describes the state
  before that approval. What the text may say of the simulations is fixed in
  `paper/audit/research/md_paper_2026-10-08/DESIGN.md`, D4. The numbers of the paper are facts of section 40 of
  `paper/main/data/facts.md`, re-derived by `paper/main/data/derive/d17_md.py`.

Folders and records added after the table above was written:

| Path | Stage | What it holds |
|---|---|---|
| `SYSTEM.md`, `LIGAND.md`, `ZINC.md`, `TOPOLOGY.md`, `RUN.md`, `GR.md`, `N2.md`, `PRODUCTION.md` | S1 to S5 | the method records, one per step, each dated |
| `topology/` | S3 | the built systems of the zinc sets (`bonded_2e`, `bonded_2s`, and the two with the candidate terms at the sulfonamide nitrogen), each with its build record |
| `test/` | S4 | the tests on the Pi and the rehearsal run |
| `gpu/` | S5 | what ran on the rented machine: the queue, the task lists, the switch that shuts the machine down when the Pi is not seen (`GPU_RUN.md`) |
| `analysis/` | S6 | the specification of the analysis (`SPEC.md`), the two programs and the record of their comparison (`RECHECK.md`) |
| `production/` | S5, S6 | the gate on the zinc geometry, the collection of the results, the conversion times, the site checks, the series and final frames, and the assembler of `RESULTS.md` |
| `results/` | S6 | the results, as listed above |
| `qm/gr/`, `qm/n2/` | S2 | the two quantum-chemical tests of the terms at the sulfonamide nitrogen, both failed; `qm/gr/relax/README.md` says which record of the relaxations is which |

Where the large files are: the trajectories (13.4 GB) are on the Pi in `/root/cqc-sfa-md/runs/` and are not in
git; their checksums and sizes are in the manifests named in `md/results/run_health.json`.
