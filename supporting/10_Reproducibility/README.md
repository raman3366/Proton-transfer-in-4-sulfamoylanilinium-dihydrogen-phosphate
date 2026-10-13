# 10 Reproducibility

The complete source tree is the git repository; this package was assembled from its working tree, whose last commit
was **`672a72e29e2567d47f3c522c130805b7198937be`** (`refs/heads/main`). This folder holds a snapshot of the code, the inputs, the software versions and the compute records, so that every step can be read and
repeated.

| Path | What |
|---|---|
| `SOFTWARE_VERSIONS.md` | Programs and versions, hardware, levels of theory. |
| `inputs/` | The deposited CIF (CCDC 1016950, with the embedded reflections), the start geometries cut from it (`start_geometries/`), the case definitions (`case_definitions/`: charge, level, solvent, atom labels, fragments), and `plbd_hca2_azm_sfm.json`, the rows of the Protein-Ligand Binding Database (PLBD) for acetazolamide and sulfanilamide with carbonic anhydrase II that `derive/d14_completion.py` reads (retrieved 29 Sep 2026; sources in 09_References/SOURCE_NOTES.md). |
| `scripts/dft/` | The DFT code: `pipeline.py` (optimise → proton check → frequencies → FMO → Multiwfn → plots → report), `scan.py`, `eps_ladder.py`, `fsp.py`, `zpe_barrier.py`, `tddft.py`, `dscf.py`, `nlo.py`, `interaction.py`, `mwfn*.py`, `hirshfeld.py`, `periodic_qe.py`; shared settings in `common.py`. |
| `scripts/structures/` | `build_models.py`: builds the models from the CIF. |
| `scripts/dock/` | The docking scripts (`targets.py`, `prep_receptor.py`, `prep_ligands.py`, `run_vina.py`, `run_ad4zn.py`, `analyse.py`, `contacts.py`, `production_*.py`, `control_gyrb.py` for the control receptor, the `run_*.sh` drivers) and `vendor/` (AutoDock4Zn parameters and helpers). |
| `scripts/ops/` | The compute set-up: VM set-up scripts (`vm/setup_*.sh`), `run_case.sh`, the burst runner `vm/burst.py` with the burst definitions `vm/bursts/*.json`, and the Pi watchdog and status scripts; for the rented GPU machine of the molecular dynamics, the scripts that created, guarded and deleted it (`ibm/cqc_gpu.py`, `ibm/cqc_gpu_up.py`, `pi/cqc_gpu_watchdog.py`, `pi/cqc_gpu_capguard.py`, `pi/cqc_gpu_status.py`, `vm/setup_gpu.sh`) with their documentation (`ibm/GPU.md`, `pi/GPU_WATCH.md`); and `make_submission.py`, which builds the Word files of the cover letter, the highlights and the declaration of interests and assembles the files to upload in `submission/JMGM/` of the repository (the packaging script is in the repository, not here). **Sanitised:** the IBM Cloud login argument is redacted and network addresses are masked in `ibm/cqc_vm.py`, `ibm/cqc_vm_up.py`, `pi/GPU_WATCH.md`, `pi/cqc_credit.py`, `pi/cqc_gpu_watchdog.py` and `pi/cqc_watchdog.py`; no key, machine address or instance id is included. The scripts of the simulations themselves are in 11_Molecular_dynamics/md. |
| `scripts/paper/` | The paper sources: manuscripts (`main/manuscript.md`, `si/manuscript.md`), `highlights.md`, `cover_letter.md` and `declaration_of_interest.md`, the facts and table builders (`main/data/*.py`, `derive/`), the figure scripts, the ADMET and Fukui scripts (`si/compute/`), and the shared tools (`_shared/tools`: build, lint, numbers audit, overlap check, Word check, `lit.py`), figure toolkit, citation styles and Word template. `paper/README.md` has the exact build commands. |
| `benchmark/` | The VM timing benchmark on the real ion pair. |
| `vm_burst_status/` | Status snapshots of the VM bursts. |
| `vm_ledger_sanitised.csv` | VM start/stop events (instance id and price removed). |
| `gpu_ledger_sanitised.csv` | The ledger of the rented GPU machine of the molecular dynamics: its events with their times (5 lines; instance id, price and note removed). |
| `vm_archive_MANIFEST.sha256` | sha256 of the 1657 files pulled from the VM into the Pi archive. |
| `excluded_files.csv` | Every file left out of the directories this package copies (the result trees, the docking tree, the folder of the molecular dynamics, the code and the paper sources), the internal review and audit records of `paper/` and the submission folder of the repository (`submission/`, which repeats the documents of 01 to 03), with size, source location and reason. |
| `source_map.csv` | Every file of this package with the file it was copied from, or how it was generated. |
| `git_commit.txt` | The last commit id at packaging time. |

## Compute

VM `cqc-dft-01` (IBM Cloud, bx3d-48x240, us-east): **30.04 VM-hours** in 8 sessions (from the
ledger). The VM was stopped after the last burst (29 Sep 2026 13:30 IST) and deleted on 29 Sep 2026 16:21 IST, with its boot volume and network; every result is in the sha256-verified Pi archive and in this package.

| Started (IST) | Stopped (IST) | Hours |
|---|---|---:|
| 27 Sep 2026 17:45 | 27 Sep 2026 18:25 | 0.67 |
| 27 Sep 2026 19:50 | 28 Sep 2026 03:21 | 7.52 |
| 28 Sep 2026 05:35 | 28 Sep 2026 19:25 | 13.83 |
| 28 Sep 2026 22:43 | 29 Sep 2026 00:10 | 1.45 |
| 29 Sep 2026 05:55 | 29 Sep 2026 08:11 | 2.26 |
| 29 Sep 2026 09:04 | 29 Sep 2026 09:05 | 0.01 |
| 29 Sep 2026 09:11 | 29 Sep 2026 11:15 | 2.07 |
| 29 Sep 2026 11:17 | 29 Sep 2026 13:30 | 2.22 |

Molecular dynamics: one rented GPU machine `cqc-md-gpu-01` (IBM Cloud, NVIDIA L40S), 43.60 h of running and
US$37.97 by the project's ledger, as the plan of the simulations records it at the deletion of the machine
(11_Molecular_dynamics/md/PLAN.md, entry of 8 Oct 2026, 8:30 AM IST). The machine was deleted after the last run, with its boot volume and network; its
results were compared with the copies on the Pi by sha256 before that (11_Molecular_dynamics/README.md).

## How to reproduce (from a checkout of the repository)

1. **Models (any machine).** `python3 structures/build_models.py` writes `molecules/sfa/input/*.xyz` and
   `molecules/sfa/cases/*.json` from the CIF.
2. **Molecular DFT (24-core Linux machine).** Install PySCF 2.14 + geomeTRIC, D3(BJ), Multiwfn and Quantum ESPRESSO
   7.5 with `ops/vm/setup_linux.sh`, `setup_extra.sh` and `setup_qe.sh`. One case: `ops/vm/run_case.sh <case>
   <threads>` (it runs `python dft/pipeline.py <case>`). The production runs were the bursts in `ops/vm/bursts/*.json`,
   run with `python ops/vm/burst.py ops/vm/bursts/<burst>.json` (burst1a, burst1a_cont, burst1b, burst1b_cont, burst1b_cont2, burst1b_props, burst1c_mwfn, burst2, burst3a, burst3b, burst3c, burst_periodic). The scan, dielectric
   ladder, functional check, TD-DFT, ΔSCF, NLO and interaction energies are separate scripts in `dft/` called by the
   bursts; the periodic model, with the constrained N1–H1A profile, is `dft/periodic_qe.py`; the five-pair cluster
   and its QTAIM `dft/crystal_cluster_qtaim.py`; the relaxed torsion scans `dft/torsion_scan.py`. Burst 3
   consisted of 3a (Hessian of the neutral minimum), 3b (profile, single-transfer control, cluster QTAIM, first
   torsion points) and 3c (the remaining torsion points); `burst3b` notes how it was assembled from
   `ops/vm/bursts/fragments/`.
3. **Docking (Raspberry Pi 5 or any Linux machine).** Vina 1.2.7, Meeko 0.8.0, RDKit, PDBFixer/OpenMM, PROPKA,
   spyrmsd, PoseBusters, ProLIF, PLIP. In order: `python dock/targets.py`; `python dock/prep_receptor.py <target>`;
   `python dock/prep_ligands.py --dft ...`; `bash dock/run_redock_all.sh` then `python dock/analyse.py` (gate G2);
   `bash dock/run_production.sh` then `python dock/analyse.py --production`; `python dock/contacts.py`;
   `python dock/production_report.py`. AutoDock4Zn maps need x86 AutoGrid 4.2.7:
   `python dock/run_ad4zn.py --target <T> --stage maps`. The control receptor: `python dock/control_gyrb.py run`, then
   `analyse` and `report` (the commands are at the end of `07_Docking/control_gyrase_B/CONTROL.md`).
4. **Molecular dynamics (set-up and analysis on any Linux machine; the runs on a GPU).** The environment is
   `md/env/md_env.yml` (OpenMM, AmberTools, MDAnalysis). In order: the system builder `md/build_system.py` and the
   topology builder `md/build_topology.py`; the driver `md/run_md.py`, one call per run (the task files
   `md/gpu/tasks_x0.json`, `tasks_production.json` and `tasks_sensitivity.json` hold the command of every run, and
   `md/gpu/queue.py` ran them); the two analysis programs `md/analysis/analyse_md.py` and `recheck_md.py`, whose
   agreement is the gate of the analysis; then `md/production/make_results.sh`, which runs every later step in its
   order and ends with `md/RESULTS.md`. The quantum-chemical tests of the force field are the scripts of `md/qm/`.
   The method records in 11_Molecular_dynamics/md (`SYSTEM.md`, `LIGAND.md`, `ZINC.md`, `TOPOLOGY.md`, `RUN.md`, `GR.md`, `N2.md`,
   `PRODUCTION.md`) give the commands with their arguments and what each step must show before the next.
5. **Numbers, tables, figures and documents.** `paper/main/data/derive/d*.py`, then `build_facts.py` and
   `build_tables.py`; the figure scripts; the gates (`lint_prose.py`, `numbers_audit.py`, `overlap.py`); then
   `paper/_shared/tools/build.py main` and `si` (add `--spacing 2` for a double-spaced version) and
   `docx_check.py`; then `ops/make_submission.py --letters` for the Word files of the cover letter, the highlights
   and the declaration of interests. The exact commands are in `scripts/paper/README.md`.

Rebuilding the paper needs a checkout of the repository with its toolchain (pandoc, the Python environment of
`paper/README.md` and LuaLaTeX for the PDF): the tables (05_Tables) and figures (04_Figures) sit in
`paper/main/tables`, `paper/si/tables`, `paper/main/figures` and `paper/si/figures`, and the results in
`molecules/sfa/results`, `dock/` and `md/`.
