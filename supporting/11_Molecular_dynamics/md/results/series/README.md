# Time series of the eleven runs of 100 ns (rule 6 of `md/PLAN.md`)

Written by `md/production/make_series.py` from the output of the first analysis program
(`md/analysis/analyse_md.py`; specification `md/analysis/SPEC.md`, version 1.3). One file per run, gzip-compressed
CSV, one line per solute frame (10 ps): 130 frames of the equilibration after the release of the ligand
(stages 4, 5 and 6; 1.3 ns) followed by the 10,000 frames of the production. 11 files, 10,130 lines of
data each, 1,630 kB together.

No number is formed in these files. Every cell is the text of a cell of `timeseries_equil_solute.csv` (rows of
stages 4 to 6) or `timeseries_solute.csv` in the folder `analysis/` of the run, except the first two columns as
described below. The same columns of the second program (`md/analysis/recheck_md.py`, folder
`analysis_recheck/`) were read and compared before each file was written: the largest difference of a number is
0.0000 (tolerance 0.001), and every state, sign and flag is the same. `SOURCES.json` names every file
read with its sha256; `SHA256SUMS` holds the sums of the files of this folder.

## Columns

| Column | What it is | Source column |
|---|---|---|
| `time_ns` | Time in ns; 0 is the start of the production. The released part of the equilibration runs from −1.29 to 0.00, the production from 0.01 to 100.00. Written with two decimals | `time_ns` |
| `part` | `eq4`, `eq5`, `eq6`: equilibration stage 4 (protein backbone held with 250 kJ mol⁻¹ nm⁻²), stage 5 (50 kJ mol⁻¹ nm⁻²) and stage 6 (nothing held); `prod`: production. The ligand is released at the start of stage 4 | `stage` |
| `d1_A` | d1, distance from the aniline nitrogen N1 of the ligand to CZ of Phe131 (Å) | `d1_N1_F131CZ_A` |
| `d2_A` | d2, distance from N1 to OD1 of Asn67 (Å) | `d2_N1_N67OD1_A` |
| `lig_rmsd_A` | RMSD of the 11 heavy atoms of the ligand from the crystal pose of 6RL9 after the fit of the 258 Cα atoms; smallest of the four equivalent labellings (Å) | `lig_rmsd_A` |
| `state` | `C` crystal-like (d1 < 6.7 Å and d2 > 5.2 Å), `R` rotated-like (d1 > 6.7 Å and d2 < 5.2 Å), `I` intermediate (neither) | `state` |
| `n2_volume` | Signed volume of the unit vectors from the sulfonamide nitrogen N2 to Zn, H and S; positive in the crystal pose | `n2_volume` |
| `n2_sign` | Sign of `n2_volume`: 1, −1 or 0 | `n2_sign` |
| `hb_T199` | 1 if the hydrogen bond N2–H···OG1(Thr199) is present (N2···OG1 at most 3.5 Å and the angle at H at least 135°), else 0 | `hb_N2H_T199OG1` |

The verdict of a run reads the last 20 ns (the last 2,000 lines): CONVERTED or STAYED if at least 90 % of these
frames are in the state opposite to the start or in the state of the start, MIXED otherwise. The starts X and N
are crystal-like, the start R is rotated-like. The definitions are those of `md/analysis/SPEC.md` (sections 4,
9A and 9B) and of the plan (rules 1 to 3).

## Files

| File | Run | Lines of data (equilibration + production) | Size (bytes) | End of stage 4; 5 (ns) |
|---|---|---|---|---|
| `bonded_2e_X1.csv.gz` | bonded_2e/X1 | 10130 (130 + 10000) | 147272 | -0.90; -0.50 |
| `bonded_2e_X2.csv.gz` | bonded_2e/X2 | 10130 (130 + 10000) | 147480 | -0.90; -0.50 |
| `bonded_2e_X3.csv.gz` | bonded_2e/X3 | 10130 (130 + 10000) | 147300 | -0.90; -0.50 |
| `bonded_2e_N1.csv.gz` | bonded_2e/N1 | 10130 (130 + 10000) | 146778 | -0.90; -0.50 |
| `bonded_2e_N2.csv.gz` | bonded_2e/N2 | 10130 (130 + 10000) | 147437 | -0.90; -0.50 |
| `bonded_2e_N3.csv.gz` | bonded_2e/N3 | 10130 (130 + 10000) | 147666 | -0.90; -0.50 |
| `bonded_2e_R1.csv.gz` | bonded_2e/R1 | 10130 (130 + 10000) | 146895 | -0.90; -0.50 |
| `bonded_2e_R2.csv.gz` | bonded_2e/R2 | 10130 (130 + 10000) | 146941 | -0.90; -0.50 |
| `bonded_2e_R3.csv.gz` | bonded_2e/R3 | 10130 (130 + 10000) | 146724 | -0.90; -0.50 |
| `bonded_2s_R1.csv.gz` | bonded_2s/R1 | 10130 (130 + 10000) | 149981 | -0.90; -0.50 |
| `bonded_2e_n2fit_R1.csv.gz` | bonded_2e_n2fit/R1 | 10130 (130 + 10000) | 155258 | -0.90; -0.50 |

`bonded_2s_R1` and `bonded_2e_n2fit_R1` are the two sensitivity runs of start R (the other zinc set; the
candidate terms at N2, which gate G-R2 did not accept). They are in no count of rule 3.

## To read a file

```python
import csv, gzip
with gzip.open('md/results/series/bonded_2e_R1.csv.gz', 'rt') as fh:
    rows = list(csv.DictReader(fh))
```
