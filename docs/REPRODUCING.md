# Reproducing and inspecting the records

Commands below assume a checkout and Python 3.10+. Repository-maintenance tools use the standard library unless explicitly marked optional. Installing the computational stack or running a simulation is a separate task. Do all writable work under `work/` or another new directory; preserve the archived files.

## 1. Verify the package and regenerate numerical exports

```sh
python3 tools/validate_repository.py
python3 tools/check_results.py
python3 tools/export_results.py
git diff -- data/derived
```

The export should be deterministic. `validate_repository.py` checks the inventory and the byte identity of all original ZIP members and imported handover files. `check_results.py` recomputes selected differences from stored energies, checks their printed tables, and compares 602 historical trajectory hashes with the newly hashed external inventory. These are checks on stored records, not reruns of scientific calculations.

Optional format checks:

```sh
python3 -m venv tmp/validation-env
tmp/validation-env/bin/python -m pip install -r tools/requirements-validation.txt
tmp/validation-env/bin/python tools/validate_repository.py --extended
```

The [inventory](../metadata/inventory.csv) contains relative paths, byte sizes and SHA-256 hashes of every non-ignored publication file except itself. Regenerate it with `python3 tools/make_inventory.py` only after reviewing intentional changes. It does not recursively include `.git/`, `work/`, `tmp/` or the external MD payload. Existing historical manifests retain their own original scopes.

## 2. Reconstruct the original project paths

```sh
python3 tools/prepare_workspace.py molecular --dest work/molecular
python3 tools/prepare_workspace.py project --dest work/project
```

Both refuse to overwrite an existing destination and execute no scientific program. `molecular` copies the 18 cases and their starting geometries into the layout expected by `common.py`:

```text
work/molecular/
  dft/                         archived PySCF scripts
  molecules/sfa/cases/          case JSONs
  molecules/sfa/input/          starting XYZs and torsion plan
```

Real copies are used because the archived scripts resolve `__file__`; symlinking `common.py` would resolve back to the archive and break its root calculation. The helpers verify every case's `xyz` path. `DFT_THREADS` and `DFT_MAX_MEMORY_MB` control resources. The pipeline can launch additional stages and automatic PCM calculations; inspect its flags before running it. A single-stage invocation after installing the stated dependencies is `python dft/pipeline.py sfa_pair_pcm_d3 --stage optimize`, **from the reconstructed directory**. This command was not executed during repository preparation.

`project` uses the original handover `source_map.csv` to restore 10,142 verified relative paths, including `dft/`, `molecules/sfa/`, `dock/`, `md/`, structure builders and available paper scripts/tables. VM-archive paths are mapped back to the same project root. It checks every copied file's source checksum, decompresses the two large Molden copies without changing their contents, rejects conflicting/path-traversing mappings, and records omissions in `RESTORATION.json`. Generated-only records and external MD locations (121 entries) are not guessed into paths. The public original tables remain in `supporting/05_Tables/`.

This resolves relative file placement. It does **not** resolve hard-coded `/root/`, `/home/ubuntu/` or `/mnt/scratch/` paths, absent binaries, map files, trajectories or all manuscript-build assets.

## 3. Quantum ESPRESSO and pseudopotentials

The 45 root QE inputs retain original `outdir`/`pseudo_dir` values. The recovered output tree has 49 runs, including historical alternatives; use [SCIENTIFIC_RECORDS.md](SCIENTIFIC_RECORDS.md) to select the reported runs.

The required potentials are SSSP **PBE efficiency 1.3.0** and **SG15 ONCV PBE v1.2**, as specified by the archived driver. See [pseudopotentials.json](../metadata/pseudopotentials.json) for all twelve filenames and recorded MD5 values. SSSP is distributed through [Materials Cloud](https://www.materialscloud.org/discover/sssp); the historical source record is `rcyfm-68h65`, with DOI [10.24435/materialscloud:f3-ym](https://doi.org/10.24435/materialscloud:f3-ym). SG15's source is [quantum-simulation.org](http://www.quantum-simulation.org/potentials/sg15_oncv/); the driver names `sg15_oncv_upf_2020-02-06.tar.gz`, archive MD5 `0d72299c3962fb9de17a22d81cb79734`. Do not substitute a newer library merely because it is now the default download.

After obtaining the exact potential files, prepare a disposable input with:

```sh
python3 tools/prepare_workspace.py qe inputs/quantum_espresso/relax_i.in \
  --dest work/qe-relax-i --pseudo-dir /absolute/path/to/sssp-1.3.0 \
  --scratch-dir /absolute/path/to/scratch/relax_i
```

This verifies the six named files against recorded MD5 hashes and changes only `pseudo_dir` and `outdir` in a new input. It refuses missing/mismatched potentials. Create the scratch directory on the compute machine and use the recorded QE 7.5 settings. No QE run or pseudopotential download was performed here; positive execution with real potentials remains untested. The MD5 values are historical identity records, not a claim of independently recomputed pseudopotential SHA-256 hashes.

## 4. Docking and MD

For browsing, use root `docking/` and `molecular_dynamics/`. For the source layout use `work/project/dock/` and `work/project/md/` after reconstruction. Read the preserved method records before executing any script:

- [Docking results and source map](../supporting/07_Docking/README.md), with deposited receptors, raw runs and preparation diagnostics.
- [MD methods overview](../supporting/11_Molecular_dynamics/README.md), [topology construction](../supporting/11_Molecular_dynamics/md/TOPOLOGY.md) and [analysis specification](../supporting/11_Molecular_dynamics/md/analysis/SPEC.md).
- [MD environment](../molecular_dynamics/scripts/env/README.md): the explicit lock is Linux **aarch64**, not a portable macOS/x86 installation recipe. Remove the original `prefix` only in a working copy of the YAML when creating another environment. Production used a separate CUDA machine.

The root Vina configurations resolve `receptor.pdbqt` from their receptor directory. Ligand JSON paths and source docking scripts refer to the old project root. AutoDock4Zn map generation also depends on vendor files and an appropriate AutoGrid executable; its source helper contains remote-execution machinery. Configure local binaries/paths in a copied workspace, and do not run the historical GPU/cloud queue/deadman scripts against current infrastructure.

MD additionally needs the external matching `systems/` and `runs/` data. `AMBERHOME` and `CQC_SFA_MD_SYSTEMS` are supported by the topology builder, but other helpers still hard-code bulk/tool paths and need local configuration. Rebuilding ligand parameters needs the restored DFT reference files. The complete source `md/` adds documents/QM scripts omitted from the smaller ZIP, but does not provide the absent full trajectories or all solvated systems in git. The few ligand-only `.prmtop`/`.inpcrd` files included in the archive are not the complete solvated run topologies.

## 5. Source recovery

The original ZIP is preserved at the root. Additional files can be re-imported into a fresh checkout without `supporting/` using `python3 tools/import_handover.py /path/to/SAPH_deliverables_2026-10-12.zip`. This verifies the pinned package SHA-256, ZIP CRCs, safe paths and each selected member against the source manifest before copying. The large handover ZIP remains outside git; [source_packages.json](../metadata/source_packages.json) records its identity. Do not use this command to overwrite an existing import.
