# Proton transfer in 4-sulfamoylanilinium dihydrogen phosphate: a combined molecular and periodic DFT study with molecular docking

Research-data repository accompanying the manuscript by **C. Muthuselvi, Anantha Raman Lakshmipathi, D. Saravanakkumar and Venkatramanan Varadharajan**.

**Status:** manuscript in preparation for submission to *Journal of Molecular Graphics and Modelling* (Elsevier). The supplied records establish a prepared submission package; submission, acceptance and a publication DOI have not been independently confirmed. Repository preparation date: 13 October 2026.

SAPH is ionic in the crystal although an isolated gas-phase pair favours a neutral sulfanilamide–phosphoric-acid adduct. Molecular DFT, continuum-solvent calculations, periodic DFT and hydrogen-bond analyses examine this difference. Docking and molecular dynamics examine sulfanilamide binding to carbonic anhydrase II. The MD interpretation is conditional on a force field whose sulfonamide-nitrogen terms failed two quantum-chemical tests; docking scores are not binding affinities.

## Start here

- [Supplementary PDF: Sections S1–S16, Tables S1–S40](16_SAPH_Supplementary_Material.pdf)
- [All 40 supplementary tables, with Markdown captions/footnotes and CSV copies](supporting/05_Tables/SI/)
- [All seven main-text tables](supporting/05_Tables/main/)
- [Numerical exports with source files, source keys and units](data/derived/README.md)
- [Reproduction guide](docs/REPRODUCING.md), [data availability](docs/DATA_AVAILABILITY.md), [remaining gaps](docs/DATA_GAPS.md)

The original [749-file ZIP](17_SAPH_Supplementary_structures.zip) and its [README.txt](README.txt) are preserved. Its scientific directories were extracted directly at the root without changing their contents. Additional records were recovered from the author's verified 12 October handover package on the research machine. The smaller ZIP is **not** the complete research dataset.

## Contents and manuscript mapping

| Directory | Included records | Verified SI mapping |
|---|---|---|
| [molecular/](molecular/) | 16 optimised molecular XYZ structures, energies in comments | S1–S2; Table S2 |
| [periodic/](periodic/) | Nine fixed-cell structures as CIF and extended XYZ | S7; Tables S13 and S15 |
| [clusters/](clusters/) | Two five-ion-pair clusters, 135 atoms each | S9; Table S22 |
| [inputs/](inputs/) | 18 PySCF case definitions, starts and scripts; 45 QE inputs | S1, S5–S12 |
| [docking/](docking/) | Prepared receptors/ligands, boxes, poses, control and hCA II sensitivity tests | S13, S15; Tables S28–S34, S36–S37 |
| [molecular_dynamics/](molecular_dynamics/) | Inputs, scripts, summaries, compact series, final frames and force-field tests | S16; Tables S38–S40 |
| [tddft/](tddft/) | Excitation tables in CSV | Tables S25–S26 |
| [supporting/06_DFT_results/](supporting/06_DFT_results/) | Recovered molecular outputs, scan/ladder/method results, 49 raw QE outputs, wavefunctions, cluster and torsion results | S2–S12, including missing-from-ZIP Tables S7–S9 and S15 sources |
| [supporting/07_Docking/](supporting/07_Docking/) | Raw docking runs, logs, deposited PDB sources and preparation records, including failed/diagnostic records | S13 and S15 |
| [supporting/08_ADMET_and_Fukui/](supporting/08_ADMET_and_Fukui/) | Stored descriptor and screening outputs | S9 and S14 |
| [supporting/10_Reproducibility/](supporting/10_Reproducibility/) | Code, source map, original CIF, software records and historical exclusions | Cross-workflow provenance |
| [supporting/11_Molecular_dynamics/](supporting/11_Molecular_dynamics/) | Detailed MD methods, analysis specification, QM tests and historical trajectory manifests | S16 |
| [metadata/](metadata/) | File inventory, import provenance, pseudopotential records and external MD inventory | Repository provenance |

The nine primary MD runs are `bonded_2e/{X1,X2,X3,N1,N2,N3,R1,R2,R3}`. `bonded_2s/R1` and `bonded_2e_n2fit/R1` are the two sensitivity runs and enter no primary count. Two `X0` runs validate the parameter sets; benchmark/preflight records have separate roles. See the [run index](data/derived/md_run_index.csv) and [scientific-record status](docs/SCIENTIFIC_RECORDS.md).

## Reuse and verification

From a checkout, using Python 3.10 or newer:

```sh
python3 tools/validate_repository.py
python3 tools/check_results.py
python3 tools/export_results.py
python3 tools/prepare_workspace.py molecular --dest work/molecular
```

The last command reconstructs the molecular input layout without executing scientific software. A fuller `project` reconstruction and QE path configuration are described in [REPRODUCING.md](docs/REPRODUCING.md). Never run the archived scientific scripts directly in the preserved source directories.

**Tested during repository preparation:** archive/source-byte preservation, inventory SHA-256 checks, selected structured-file parsers, input references, stored-results arithmetic, compact MD series/final-frame consistency, numerical export regeneration, and directory reconstruction. The optional extended validation parses CIF with Gemmi and YAML with PyYAML. **DFT, docking and MD calculations were not rerun.** Detailed limits and commands are in [VALIDATION.md](docs/VALIDATION.md).

Software versions below are those reported in the supplied records, not an environment installed or independently reproduced here:

| Workflow | Recorded software |
|---|---|
| Molecular DFT | PySCF 2.14.0, geomeTRIC 1.1.1, D3(BJ) |
| Periodic DFT | Quantum ESPRESSO 7.5; SSSP efficiency 1.3.0; SG15 ONCV v1.2 |
| Wavefunction analyses | Multiwfn 2026.9.20 |
| Docking | AutoDock Vina 1.2.7, AutoGrid 4.2.7/AutoDock4Zn, Meeko 0.8.0 |
| MD | OpenMM 8.6.1, AmberTools 26.0, ParmEd 4.3.1, MDAnalysis 2.10.0 |

See the [software record](supporting/10_Reproducibility/SOFTWARE_VERSIONS.md) and [MD environment](molecular_dynamics/scripts/env/README.md) for further packages and platform constraints.

## Availability, citation and licenses

Full MD trajectories remain outside git and **have not been deposited in a public research-data archive**. A freshly verified inventory covers 14,091,291,987 bytes of external MD run/system files, including 11,192,273,928 bytes in 674 trajectories. This includes primary, sensitivity, validation and other records. The [deposit plan](docs/TRAJECTORY_DEPOSIT.md) distinguishes them. Pseudopotentials, software installations and selected bulky/intermediate records are also excluded; see [DATA_GAPS.md](docs/DATA_GAPS.md).

Use [CITATION.cff](CITATION.cff) for author names and ORCIDs verified against the supplied SI. Until a DOI exists, cite this repository URL and the commit used. No manuscript DOI, acceptance status or publication year is asserted.

Original code is [MIT licensed](LICENSES/MIT.txt); original data and documentation are [CC BY 4.0](LICENSES/CC-BY-4.0.txt). Third-party code, structural data and parameter material retain their own terms. Read the [license scope](LICENSE) and [provenance/third-party notices](docs/PROVENANCE_AND_LICENSES.md).
