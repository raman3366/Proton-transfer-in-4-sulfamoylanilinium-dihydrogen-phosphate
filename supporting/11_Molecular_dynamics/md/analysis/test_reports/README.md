# Test reports of the analysis script, version 1.3 (6 Oct 2026, 11:26 AM IST)

Copies of the reports that `md/analysis/analyse_md.py` (sha256 `f75e6155de700700…`) produced on the
synthetic trajectories of `md/analysis/make_test_trajectory.py` (sha256 `b98c4a8d18a2b832…`),
made before any validation or production run existed. They replace the reports of version 1.2 in this folder
(6232 planted values, 74 broken copies); the numbers are also in section 16 of `../SPEC.md`.

| File | What it holds |
|---|---|
| `test_report.md`, `test_report.json` | planted against returned: 7339 of 7339 values within tolerance, 0 failed, over 23 synthetic runs and the campaign |
| `planted.json` | what the generator planted: the systems, the bent poses, the gate G-V of the validation runs by the generator's own evaluation |
| `check_starts.json` | the three starts against `md/system/starts.csv` (51 of 51 values agree); the references of gate G-V against the crystal start; the signed volume at N2 and the zinc to Glu106 of the starts; new: the crystal water of each start against the place of the ring (within 2.5 Å of the six ring-carbon positions: X 0, N 0, R 0; the wider count: X 1, N 1, R 0) |
| `crystal_site.json`, `crystal_site.py` | the six N–Zn–N angles of 6RL9 from `md/system/` by a separate short program with plain formulas |
| `unit_checks.json`, `unit_checks.py` | rules of sections 9A, 9B, 9D and 9E on values typed by hand (32 of 32 pass); new: the settings of the water count of section 9D from a manifest or from the defaults |
| `built_starts.json`, `built_starts.py` | the three built systems of the restrained-zinc model as runs of one frame (no simulation), analysed and compared with plain formulas, now with both water counts (largest difference 5.1e-05 Å; all agree: yes) |
| `v12_against_v13.json`, `compare_v12_v13.py` | new: the 21 synthetic runs of version 1.2 under version 1.3: the same coordinates (24625 frames in 113 files), the same values in every column and key of version 1.2, four of them under a new name (all the same: yes) |
| `mutation_results.json`, `mutate.py`, `mutate_staged.py` | the deliberately broken copies of the analysis (15 of them run on all runs, 68 noticed on the runs named for them in stage A, 0 taken on to all runs in stage B; the machine was too loaded for all copies on all runs): 83 of 83 noticed (14 of version 1.0, 26 of version 1.1, 34 of version 1.2, 9 new); and M15, a standard deviation with N − 1, which changes 7 of 7339 values by more than the tolerance |
| `driver_crosscheck.json`, `driver_crosscheck.py` | the analysis against trajectories written by `md/run_md.py`: its first test run on a stand-in system (56 monitored distances; largest difference 0.0005 Å), and, new, two runs of the driver as it stands with its hydration check (4 extension steps read, 6 hydration records compared, 18 of 18 counts at the ends of steps identical; all agree: yes) |
| `driver_now_run_driver.sh`, `driver_now_make_protocol.py`, `driver_now_protocol_test.json`, `driver_now_run_driver.log` | new: how these two driver runs were made (stand-in systems of the driver's own tests, one thread, a test protocol derived from `md/protocol.json` as it stands), and the log |
| `appendix_check.py`, `appendix_check.json`, `appendix_check.log` | new: a reader written from appendix A and sections 2, 9D and 9E of the specification alone (no MDAnalysis) makes both water counts in 123 full frames of two driver runs and one synthetic run; all agree with `analyse_md.py`: yes |
| `driver_diff.txt`, `driver_0528_facts.py`, `driver_0528.json` | new: the one line by which `md/run_md.py` differs from the file that version 1.2 named; the five test runs of the driver's maintainer (made with that earlier file and the hydration block of that hour) as the analysis reads them |
| `vm_check.txt`, `scripts_since.txt` | new: a read-only look at the rented machine (checksums of its driver and protocol, its queue with the times of the rehearsal runs, the names in its runs folders), and the lines of an earlier log that show since when the two scripts have their final checksums |
| `runs_state.py`, `runs_dates.py`, `runs_dates.json` | new: what `/root/cqc-sfa-md/runs/` held while this version was tested (3 folders of tests of the rented machine: bench, preflight_a, preflight_b; no run of the campaign), from the manifests' run ids and test flags and from file dates only |
| `run_all.sh`, `run_test.sh`, `run_tail.sh`, `run_all.log`, `timing_one_copy.txt` | how the tests were run, and the log of that run (edge cases included). `run_all.sh` was stopped in its loop over the broken copies after M16; `run_tail.sh` ran the other copies in stages and the last lines of `run_all.sh`, into the same log; `timing_one_copy.txt` has the time of one pass over all runs on the loaded machine |

New against version 1.2: the water count of section 9D is the count that the driver acts on (water oxygens within
2.5 Å of the crystal positions of the six ring carbons; settings from the manifest of the run where it gives them),
beside the wider count of version 1.2 under its new name; two synthetic runs for it, the run ext/R1 with the
settings of an older driver in its manifest, the broken copy "the six ring carbons only" reversed, nine new broken
copies, and two runs of the driver as it stands. Files of version 1.2 that are gone from this folder
(2): `compare_v11_v12.py`, `v11_against_v12.json`; `v12_against_v13.json` takes the place of
`v11_against_v12.json`.

The tests were run on the two scripts as they stood in the scratch folder (`dev/md/analysis/`); the files in
`md/analysis/` are the same files, byte for byte (the checksums above). The scripts and the log name scratch
folders that are not kept; `make_test_trajectory.py` writes the same runs again (seed 20261005).
