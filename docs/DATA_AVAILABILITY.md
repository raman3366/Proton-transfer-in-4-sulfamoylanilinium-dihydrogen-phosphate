# Data availability

The public repository is [raman3366/Proton-transfer-in-4-sulfamoylanilinium-dihydrogen-phosphate](https://github.com/raman3366/Proton-transfer-in-4-sulfamoylanilinium-dihydrogen-phosphate). Its root contains the unchanged supplementary PDF and ZIP and all 749 extracted scientific files. `supporting/` adds selected research records from the complete handover package, including the main and SI tables, raw/processed DFT outputs, docking runs and MD method/QM records.

This is a curated research package, not a claim that every file from every calculation is included. See [DATA_GAPS.md](DATA_GAPS.md) for exclusions and [SCIENTIFIC_RECORDS.md](SCIENTIFIC_RECORDS.md) for provisional and validation records.

## Full MD data

The full trajectories exist on the author's research machine but **have not yet been deposited publicly**. No trajectory DOI or download link exists in this repository. A new [external inventory](../metadata/external_md_inventory.csv) records relative paths, sizes, SHA-256 hashes and roles for the `runs/` and `systems/` trees: 2,490 files, 14,091,291,987 bytes, of which 674 trajectory files account for 11,192,273,928 bytes. This is broader than the thirteen scientifically reported/validation run manifests because it also includes benchmark/preflight records and system files. The original archive's “about 13 GB” is an approximate description of excluded MD data, not a measured size of this repository.

The [deposit plan](TRAJECTORY_DEPOSIT.md) specifies what to upload and how to connect the deposit to this repository. The author retains the full data; their location is not a public download service.

## Suggested manuscript wording after this commit is public

> The supplementary structures, calculation inputs, available molecular and periodic DFT outputs, numerical tables, prepared docking systems and docking outputs, and molecular-dynamics preparation, analysis and quantum-chemical records are available in the public research-data repository “Proton-transfer-in-4-sulfamoylanilinium-dihydrogen-phosphate” at https://github.com/raman3366/Proton-transfer-in-4-sulfamoylanilinium-dihydrogen-phosphate. The repository includes a SHA-256 inventory and documents retained provisional records, software requirements and remaining data gaps. Compact MD time series and final solute frames are included. Full MD trajectories and associated bulk system/restart files are inventoried separately and have not yet been deposited publicly; they are retained by the corresponding author. File-level and selected stored-result consistency checks have been performed; the scientific calculations have not been rerun as part of repository preparation.

Add the actual release/commit identifier and deposit DOI when available. The original manuscript and handover statements that no package is public are historical and need this update; their scientific files have not been edited to change that history.
