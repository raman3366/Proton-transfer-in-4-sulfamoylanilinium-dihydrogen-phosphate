# External MD deposit plan

**Status on 13 October 2026: files located and hashed; no public deposit, DOI or upload has been created.** The source is the author's research-machine directory `/root/cqc-sfa-md`, with relative paths rooted there. No remote files were changed.

The [fresh inventory](../metadata/external_md_inventory.csv) covers `runs/` and `systems/`, including trajectory, topology, coordinate, restart, log and manifest files. It records SHA-256 hashes computed from the source bytes, with a size/mtime stability check during hashing. It is separate from [the Git inventory](../metadata/inventory.csv). The historical [thirteen-run manifests](../supporting/11_Molecular_dynamics/trajectories_not_included/) are preserved for comparison.

## Upload units

| Upload group | Contents | Role |
|---|---|---|
| `primary-bonded_2e` | X1–X3, N1–N3, R1–R3 run directories; matching X/N/R system files | Nine primary 100 ns runs |
| `sensitivity-bonded_2s` | R1 and matching system files | Zinc-set sensitivity |
| `sensitivity-bonded_2e_n2fit` | R1 and matching system files | Rejected-candidate N2 sensitivity |
| `validation` | Both X0 directories and their systems | Two 5 ns validation runs |
| `other-records` | Bench/preflight/support records and associated systems | Preserve separately; never count as primary |

Prefer one archive per run or clearly labelled group, plus a shared systems archive. Include all equilibration and production segments, atom-order/topology files, selections, seeds/protocols, driver manifests and analysis metadata. Include the scientific run index, force-field limitations and repository commit in the deposit description. Keep the internal inventory alongside the bundles, and generate an additional SHA-256 inventory of the uploaded archive files themselves.

1. Choose Zenodo or an equivalent general research-data archive and inspect its current upload limits. Reserve a DOI in a draft record using the verified authors/ORCIDs from `CITATION.cff`.
2. On the source machine, regenerate the inventory with `python3 inventory_external.py /root/cqc-sfa-md > external_md_inventory.csv` using [the supplied tool](../tools/inventory_external.py); compare with this snapshot before packaging.
3. Create bundles outside the calculation directories. Preserve relative paths and file bytes; do not rewrite trajectories or mix sensitivity/validation runs into a primary bundle.
4. Check every bundle against its per-file manifest after extraction, then upload. Set licensing for original data consistently with this repository while retaining third-party parameter/structure terms.
5. Verify the remote file list, sizes and download checksums before publishing the deposit. Record the version DOI, repository commit and related-identifier links in both locations.
6. Update `README.md`, `docs/DATA_AVAILABILITY.md`, `CITATION.cff` if appropriate, and the manuscript's data-availability statement with the real DOI and access conditions. Do not turn a reserved DOI into a claim that files are publicly available before publication.

The MD inventory totals 14.09 GB (decimal) of run/system data; trajectories alone total 11.19 GB. These are uncompressed source sizes. Compressed bundle sizes have not been measured. No user account or archive was contacted to create a deposit during preparation.
