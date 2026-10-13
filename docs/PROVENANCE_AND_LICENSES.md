# Provenance and licensing

## Original contributions

The corresponding researcher authorized MIT for original code and CC BY 4.0 for original data/documentation on 13 October 2026, confirming authority to grant those licenses. The grant does not invent a sole copyright holder or transfer rights between authors. The [MIT text](../LICENSES/MIT.txt) refers to the respective original authors/contributors; the verified manuscript author list is in [CITATION.cff](../CITATION.cff). The [root license](../LICENSE) states the scope.

The scientific scripts' own comments identify adaptation from the authors' related DFT project; the supplied SI states that scripts were written under the authors' direction with AI assistance. These provenance statements are retained. Third-party implementations, parameters and embedded copied portions are excluded from the original-code MIT grant.

## Third-party notices

| Material | Provenance and applicable terms |
|---|---|
| Original SAPH crystallography / CIF, CCDC 1016950 | C. Muthuselvi, N. Mala, N. Srinivasan, S. Pandiarajan and R. V. Krishnakumar, *Acta Cryst.* E70 (2014), o997–o998, [10.1107/S1600536814017462](https://doi.org/10.1107/S1600536814017462). The [IUCr original article and supporting data](https://journals.iucr.org/e/issues/2014/09/00/hb7251/) identify CC BY availability and the deposited CIF. Retain its original attribution/terms; the new license does not relicense the original experiment. This is a single published supporting structure, not a redistribution of the CSD database. |
| Protein structures and extracted reference ligands | PDB entries include 3TZF, 6RL9, 3HS4, 7S2K, 4DUH and further validation/control references. `supporting/07_Docking/pdb_sources/` and the target inventory preserve entry identities. wwPDB archive data are [CC0](https://www.wwpdb.org/about/usage-policies); cite the original entry authors when using their structures. Computational modifications are described in preparation/variant records. |
| `scripts/dock/vendor/AD4Zn.dat` | Its preserved header names Garrett M. Morris, David S. Goodsell, Ruth Huey and Arthur J. Olson and explicitly grants **GPL-2.0-or-later**. [License text](../LICENSES/GPL-2.0-or-later.txt) is included. This file is not MIT/CC BY relicensed. |
| `scripts/dock/vendor/prepare_gpf4zn.py`, `zinc_pseudo.py` | Both are byte-identical to the files at AutoDock-Vina upstream commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`. Their exact URLs and hashes are in [vendor_provenance.json](../metadata/vendor_provenance.json). The upstream Vina project supplies [Apache 2.0 terms](https://vina.scripps.edu/license/); a [license copy](../LICENSES/Apache-2.0.txt) is included. Preserve the original AutoDockTools/CVS notices. Any more specific embedded notice takes precedence; these files are outside the original-code MIT grant. |
| Amber/GAFF2/ff14SB, water/ion and zinc parameter material | The package contains assigned/derived parameters, source checksums and MCPB.py output, not an original claim to the underlying parameter libraries. Existing `.frcmod`, `.mol2`, topology and template files retain third-party portions and citations. [AmberTools](https://ambermd.org/AmberTools.php) comprises components under multiple licenses; do not infer one blanket license for every parameter file or redistribute the full installation under this repository's license. Original fitted/repair contributions are CC BY only to the extent of the authors' rights. |
| Pseudopotentials, computational engines and dependencies | Not distributed here as installed software. Obtain the specified versions from their upstream providers and observe their individual licenses; the environment records are specifications, not software license grants. |

Upstream web policies and the vendor matches were checked during preparation. They do not replace component-specific notices or determine ownership of every underlying third-party parameter. The license scope deliberately retains such rights instead of claiming them for the manuscript authors.

## Source integrity and publication review

- The original ZIP contains 749 files and 21,881,155 uncompressed bytes. All extracted members remain byte-identical. The original PDF and ZIP are also identified by SHA-256 in [source_packages.json](../metadata/source_packages.json).
- The fuller 12 October handover ZIP was copied from the authorized research machine and verified against its source SHA-256. Every imported file matched its original handover manifest. [handover_import.csv](../metadata/handover_import.csv) maps 10,288 repository paths to archive members and source hashes. Only the two approximately 86 MB Molden files were losslessly gzip-compressed; decompressed bytes match their source hashes.
- Imported handover records may already contain the source package's documented sanitization. They were preserved as supplied, with that history in its `source_map.csv`. No scientific values were changed during this preparation.
- The publication screen searched text and supported compressed containers for credential signatures and address-like text. The 33 address-like matches were reviewed: package versions, enzyme classification numbers, crystallographic table numbers, a CVS revision and an OS package version. No credential match or unrelated private address was found. Published scientific author/contact information remains in the SI/CIF. No local key store, credentials, `.git` history from the source machine, private correspondence, reviewer suggestions or submission forms was imported.

See [VALIDATION.md](VALIDATION.md) for the limits of these checks. The inventory is evidence of byte identity, not proof of scientific validity or blanket ownership.
