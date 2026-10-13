# Validation record

Preparation date: 13 October 2026. The checks below concern file preservation, stored data and packaging. **Scientific calculations reproduced during this work: none.** No DFT, docking, MD, AutoGrid or topology-generation job was run, and no full scientific stack was installed.

## Commands

```sh
python3 tools/export_results.py
python3 tools/check_results.py
python3 tools/test_repository_tools.py
python3 tools/prepare_workspace.py molecular --dest work/molecular
python3 tools/prepare_workspace.py project --dest work/project
python3 tools/scan_publication.py --report tmp/privacy_review.json
python3 tools/make_inventory.py
tmp/validation-env/bin/python tools/validate_repository.py --extended --report tmp/validation.json
```

The optional validation environment contains only Gemmi, PyYAML and JSON Schema validation dependencies. The core repository checks also run with standard-library Python. Citation metadata was validated against the CFF 1.2.0 JSON schema and all four ORCID check digits were verified against the values printed in the SI.

## Scope

| Check | Meaning |
|---|---|
| Original ZIP | CRC verified; all 749 extracted scientific files equal original members byte for byte |
| Fuller handover | Source ZIP SHA-256 verified; 10,288 imported records match source-manifest hashes, including decompressed hashes of two large wavefunctions |
| Inventory | Every published file except the inventory itself has a size and SHA-256; no full trajectory or file ≥50 MiB is permitted |
| Structured files | Python JSON parsing (six documented nonstandard-NaN files), CSV readability, XYZ frame counts/finite coordinates, PDB/PDBQT coordinates, V2000 SDF counts/bond indices, Python syntax, gzip and ZIP/NPZ CRCs; optional Gemmi CIF and PyYAML parsing |
| Format distinctions | Seven AutoGrid `.maps.xyz` files are grid extents, not molecular XYZ; two empty `cofactors.pdbqt` files are valid “no cofactors” placeholders; variable-width CSVs are counted without rewriting them |
| References | Case geometry/atom labels, receptor paths, QE potential filenames, new documentation links and MD index targets |
| Stored DFT results | S7 energies/geometries/convergence flags, S8 geometries, S9 five-method differences and S15 raw QE energies plus the hybrid correction/per-formula-unit conversion checked independently |
| MD compact data | Eleven series × 10,130 rows; 10 ps spacing, production/equilibration counts, state classification from stored distances; 4,094-atom final frames checked against final distances within 0.002 Å |
| MD roles and provenance | Nine primary, two sensitivity, two validation runs; provisional/completed Rx distinction; all 602 trajectories in the historical thirteen-run manifest match fresh remote hashes and sizes |
| Workspace reconstruction | 18 molecular cases restored; full source-map reconstruction restores 10,142 original relative paths without running the scripts |
| Negative/synthetic checks | Corrupted extracted bytes and truncated XYZ rejected; a temporary non-scientific QE fixture verifies that only two path values change and that mismatched potential hashes are rejected |
| Privacy | Heuristic scan plus review of 33 address-like matches; no credential findings; published scientific attribution retained |

The source package's historical verification report is not presented as a new reproduction. Checksums do not establish chemical validity, force-field quality or optimiser convergence. The MD distance/state checks operate on supplied reduced records; no omitted trajectory was integrated or reanalysed. The historical nonstandard JSON, old absolute paths and missing dependencies remain documented limitations.

The extended preparation run passed with zero errors: 11,079 inventory records; 30 CIF parses, 290 SDF files, 197 PDB and 4,191 PDBQT files, 601 XYZ and 21 extended XYZ files, 1,322 JSON files (six containing documented nonstandard constants), four JSONL files, 453 CSV files (eight variable-width), 319 Python syntax checks, 56 gzip streams and 33 ZIP/NPZ containers. These are file counts, not independent scientific calculations. The inventory excludes its own CSV, so the published tree contains one additional file.

## Publication completion

The final staged-tree review, Git commit and remote verification are completed after the inventory and validation pass. The commit identity is available from Git history; it is not embedded into a self-referential file checksum.
