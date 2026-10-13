# Last solute frame of the eleven runs of 100 ns

Written by `md/production/make_final_frames.py`. Each file is the `final_solute.pdb` that the driver
(`md/run_md.py`) wrote at the end of the last production segment of a run, compressed with gzip and otherwise
unchanged: the 4,094 atoms of the solute (protein with its acetyl cap, zinc, ligand) at 100 ns, coordinates in
Å as integrated, with the cell in the `CRYST1` line. Water and ions are not in these files. 11 files,
758 kB together. `SOURCES.json` names every source with its sha256, which is also the sha256 of the
uncompressed content of the copy; `SHA256SUMS` holds the sums of the files of this folder.

Check made for every file: d1 (N1 of the ligand to CZ of Phe131) and d2 (N1 to OD1 of Asn67) measured on the
frame equal the last line of the time series of the run (`md/results/series/`, 100.00 ns) within 0.0008 Å
(limit 0.002 Å; the coordinates of a PDB file carry three decimals).

## Numbering in the files

Atoms carry the serial number "index in the topology + 1". Residues are numbered in sequence: the acetyl cap is
residue 1 and His3 of the deposited structure 6RL9 is residue 2, so a residue with the deposited number *n*
is residue *n* − 1 up to residue 125 and *n* − 2 from residue 127 on (6RL9 has no residue 126; Phe131 is
residue 129, Asn67 residue 66). The zinc is residue 260 (`ZN1`) and the ligand residue 261 (`SN1`; N1 is the
aniline nitrogen, N2 the sulfonamide nitrogen, `HN21` its hydrogen). Histidines are named by their tautomer
(`HID`, `HIE`); the three zinc ligands His94, His96 and His119 carry the names `HD1`, `HD2` and `HE1` of the
bonded zinc model.

## Files

| File | Run | Size (bytes) | d1; d2 on the frame (Å) | State at 100 ns |
|---|---|---|---|---|
| `bonded_2e_X1_final_solute.pdb.gz` | bonded_2e/X1 | 68657 | 4.192; 8.175 | crystal-like |
| `bonded_2e_X2_final_solute.pdb.gz` | bonded_2e/X2 | 68764 | 4.389; 7.555 | crystal-like |
| `bonded_2e_X3_final_solute.pdb.gz` | bonded_2e/X3 | 69060 | 3.810; 7.417 | crystal-like |
| `bonded_2e_N1_final_solute.pdb.gz` | bonded_2e/N1 | 68701 | 5.074; 7.380 | crystal-like |
| `bonded_2e_N2_final_solute.pdb.gz` | bonded_2e/N2 | 68718 | 4.526; 7.679 | crystal-like |
| `bonded_2e_N3_final_solute.pdb.gz` | bonded_2e/N3 | 69812 | 5.101; 7.588 | crystal-like |
| `bonded_2e_R1_final_solute.pdb.gz` | bonded_2e/R1 | 68722 | 4.797; 7.643 | crystal-like |
| `bonded_2e_R2_final_solute.pdb.gz` | bonded_2e/R2 | 68802 | 4.931; 7.473 | crystal-like |
| `bonded_2e_R3_final_solute.pdb.gz` | bonded_2e/R3 | 69423 | 5.272; 7.979 | crystal-like |
| `bonded_2s_R1_final_solute.pdb.gz` | bonded_2s/R1 | 68581 | 4.928; 9.506 | crystal-like |
| `bonded_2e_n2fit_R1_final_solute.pdb.gz` | bonded_2e_n2fit/R1 | 68766 | 3.650; 7.455 | crystal-like |

`bonded_2s_R1` and `bonded_2e_n2fit_R1` are the two sensitivity runs of start R (the other zinc set; the
candidate terms at N2, which gate G-R2 did not accept). They are in no count of rule 3 of the plan.
