# Remaining data gaps and limitations

Initial inspection covered the working directory and relevant SAPH submission folders in the local Downloads directory. With the author's authorization, the source repository and MD data on the connected research machine were then inspected. Its clean source revision was `672a72e29e2567d47f3c522c130805b7198937be`; the 12 October handover package reports the same revision. The older JICS package and unrelated projects were not mixed into the published scientific data.

| Initial concern | Finding after recovery |
|---|---|
| Supplementary PDF absent from ZIP | Present beside the ZIP; the supplied 96-page PDF is included unchanged |
| Scan, dielectric ladder and method results missing | Recovered JSON/CSV/geometry records, original SI tables S7–S9 and source-linked numerical exports are now included |
| Periodic inputs without raw output | Recovered 49 QE `.out` files plus processed convergence, relaxation and profile records |
| Provisional versus final QM records | Both retained and explicitly mapped in [SCIENTIFIC_RECORDS.md](SCIENTIFIC_RECORDS.md) |
| Original project paths | Restorable relative paths via the source map; absolute machine paths and external programs still require configuration |
| Missing license | Author authorized MIT for original code and CC BY 4.0 for original data/documentation, with third-party exclusions |

Still excluded or unresolved:

- Full MD trajectories, full system builds and restart payloads remain on the research machine. They are inventoried but not publicly deposited. Compact series/final frames cannot replace full trajectories for arbitrary reanalysis.
- Pseudopotential files and complete software installations are absent. The 12 required SSSP/SG15 file names and recorded MD5 values are available in [pseudopotentials.json](../metadata/pseudopotentials.json); no SHA-256 of an unavailable potential is invented. Downloaded potentials must match these records.
- AutoGrid maps, many scratch grids/checkpoints and some original stdout streams are excluded by the source package. Authentic processed JSON is available for workflows whose complete engine stdout is not included. See the historical [excluded-file ledger](../supporting/10_Reproducibility/excluded_files.csv). That ledger describes the source handover's exclusions, not the new Git inventory.
- This repository intentionally omits submission letters, reviewer suggestions, journal forms, cloud provisioning ledgers, literature copies and manuscript Word files. Figure/table builders are preserved as available; rebuilding every final manuscript document also requires omitted document assets and the original publishing toolchain.
- Since 13 October 2026 the current tree also omits four text sources that the handover package carried under `supporting/10_Reproducibility/scripts/paper/`: the manuscript (`main/manuscript.md`), its highlights (`main/highlights.md`), the cover letter (`main/cover_letter.md`) and the supplementary text (`si/manuscript.md`). Release v1.0.0 and its archive copy contain them and are unchanged. The package's own lists ([README](../supporting/10_Reproducibility/README.md), [source map](../supporting/10_Reproducibility/source_map.csv)) are kept byte for byte and still name them; [handover_import.csv](../metadata/handover_import.csv) and the [inventory](../metadata/inventory.csv) list only what is present. The supplementary material itself remains as the PDF at the repository root.
- Six historical JSON records contain nonstandard `NaN` tokens (timing or unavailable-analysis fields): the two pair `opt_summary.json` files, `mep_summary.json` for cation gas and phosphate gas/PCM, and the recovered MD `results/independent/physical_review/compare.json`. Their source bytes are retained. Python can read them; strict JSON consumers must handle this explicitly. No replacement values have been inserted.
- Some CSV files are variable-width records rather than rectangular tables; the validator reports this count. CIF/SDF/PDB readability and checksum validity do not establish chemical correctness or force-field validity.
- End-to-end DFT, docking, MD, map generation and topology rebuilding have not been tested here. In particular, historical claims of byte-identical rebuilt topologies and agreement between two MD analysis programs are source records, not fresh reproductions by this repository-preparation session.

The manuscript's central numerical records are now present, but the package should not be described as a fully reproduced or exhaustive computational archive. Scientific qualifications (including incomplete optimisation flags and failed force-field tests) remain part of the data.
