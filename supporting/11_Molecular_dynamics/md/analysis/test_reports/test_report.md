# Synthetic test of analyse_md.py: planted against returned

7339 of 7339 values within tolerance (0 failed). Tolerances: 0.001 Å, 0.01°, identical counts; against `md/system/starts.csv` (three decimals) 0.0015 Å and 0.06°.

| Run | Quantity | Planted | Returned | Difference | Tolerance | Source of the planted value | Result |
|---|---|---|---|---|---|---|---|
| test/N1 | solute frames seen | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| test/N1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| test/N1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| test/N1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/N1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/N1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/N1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/N1 | full frames seen | 123 | 123 | 0 | 0 | construction | passed |
| test/N1 | time seen (ns) | 12.34 | 12.34 | 0 | 1e-09 | construction | passed |
| test/N1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/N1 | frames in which the ligand was in another periodic image | 176 | 176 | 0 | 0 | construction | passed |
| test/N1 | frames in which the zinc was in another periodic image | 112 | 112 | 0 | 0 | construction | passed |
| test/N1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | lines of timeseries_solute.csv | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/N1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 8e-06 | 8e-06 | 0.001 | second formulas | passed |
| test/N1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000821 | 0.000821 | 0.01 | second formulas | passed |
| test/N1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 2.1e-05 | 2.1e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000941 | 0.000941 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000617 | 0.000617 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000682 | 0.000682 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000526 | 0.000526 | 0.01 | second formulas | passed |
| test/N1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5.6e-05 | 5.6e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.2e-05 | 5.2e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.2e-05 | 5.2e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 1.7e-05 | 1.7e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.00066 | 0.00066 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| test/N1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/N1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000698 | 0.000698 | 0.01 | second formulas | passed |
| test/N1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.001079 | 0.001079 | 0.01 | second formulas | passed |
| test/N1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000698 | 0.000698 | 0.01 | second formulas | passed |
| test/N1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000512 | 0.000512 | 0.01 | second formulas | passed |
| test/N1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | every frame: n2_volume, largest difference | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/N1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | lines of timeseries_full.csv | 123 | 123 | 0 | 0 | construction | passed |
| test/N1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/N1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/N1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| test/N1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| test/N1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| test/N1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.4e-05 | 2.4e-05 | 0.001 | second formulas | passed |
| test/N1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| test/N1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/N1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | verdict | STAYED | STAYED | – | equal | construction | passed |
| test/N1 | verdict provisional (run not finished) | yes | yes | – | equal | construction | passed |
| test/N1 | run complete | no | no | – | equal | construction | passed |
| test/N1 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| test/N1 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | frames of the verdict window | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | count_C_window | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | count_C | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| test/N1 | count_C (design against the planted frames) | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/N1 | frames with the hydrogen bond N2-H...OG1 | 988 | 988 | 0 | 0 | second formulas | passed |
| test/N1 | full frames with a water within 2.6 A of the zinc | 12 | 12 | 0 | 0 | construction | passed |
| test/N1 | frames with labelling identity | 1234 | 1234 | 0 | 0 | second formulas | passed |
| test/N1 | frames with labelling ring_flip | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | frames with labelling oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | frames with labelling ring_flip+oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | lig_rmsd_A: mean (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| test/N1 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | lig_rmsd_A: mean_last_half (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| test/N1 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | lig_rmsd_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | lig_rmsd_A: block means, largest difference (A) | 0 | 2.91617e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_N2_A: mean (A) | 2.23727 | 2.23727 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_N2_A: sd (A) | 0.090194 | 0.090194 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_N2_A: mean_last_half (A) | 2.23727 | 2.23727 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_N2_A: sd_last_half (A) | 0.090194 | 0.090194 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_N2_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 2.40782e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | tau_fold_deg: mean (deg) | 148.632 | 148.632 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | tau_fold_deg: sd (deg) | 0 | 7.1e-05 | 7.1e-05 | 0.01 | second formulas | passed |
| test/N1 | tau_fold_deg: mean_last_half (deg) | 148.632 | 148.632 | 4e-06 | 0.01 | second formulas | passed |
| test/N1 | tau_fold_deg: sd_last_half (deg) | 0 | 6.8e-05 | 6.8e-05 | 0.01 | second formulas | passed |
| test/N1 | tau_fold_deg: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | tau_fold_deg: block means, largest difference (deg) | 0 | 4.81458e-07 | 0 | 0.01 | second formulas | passed |
| test/N1 | ang_C4_S_N2_deg: mean (deg) | 109.483 | 109.483 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.2e-05 | 5.2e-05 | 0.01 | second formulas | passed |
| test/N1 | ang_C4_S_N2_deg: mean_last_half (deg) | 109.483 | 109.483 | 0 | 0.01 | second formulas | passed |
| test/N1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.3e-05 | 5.3e-05 | 0.01 | second formulas | passed |
| test/N1 | ang_C4_S_N2_deg: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 5.96364e-07 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | d_S_N2_A: mean (A) | 1.67745 | 1.67745 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | d_S_N2_A: mean_last_half (A) | 1.67745 | 1.67745 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | d_S_N2_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d_S_N2_A: block means, largest difference (A) | 0 | 2.55567e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/N1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/N1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/N1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/N1 | ca_rmsd_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d1_N1_F131CZ_A: mean (A) | 3.35615 | 3.35615 | 0 | 0.001 | second formulas | passed |
| test/N1 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | d1_N1_F131CZ_A: mean_last_half (A) | 3.35615 | 3.35615 | 0 | 0.001 | second formulas | passed |
| test/N1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | d1_N1_F131CZ_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 1.10093e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d2_N1_N67OD1_A: mean (A) | 7.72547 | 7.72547 | 0 | 0.001 | second formulas | passed |
| test/N1 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | d2_N1_N67OD1_A: mean_last_half (A) | 7.72547 | 7.72547 | 0 | 0.001 | second formulas | passed |
| test/N1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/N1 | d2_N1_N67OD1_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 7.58427e-08 | 0 | 0.001 | second formulas | passed |
| test/N1 | n2_volume: mean | 0.509469 | 0.509469 | 0 | 0.001 | second formulas | passed |
| test/N1 | n2_volume: sd | 0.553213 | 0.553213 | 0 | 0.001 | second formulas | passed |
| test/N1 | n2_volume: mean_last_half | 0.509469 | 0.509469 | 0 | 0.001 | second formulas | passed |
| test/N1 | n2_volume: sd_last_half | 0.553213 | 0.553213 | 0 | 0.001 | second formulas | passed |
| test/N1 | n2_volume: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | n2_volume: block means, largest difference | 0 | 2.46347e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE1_A: mean (A) | 4.01154 | 4.01154 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE1_A: sd (A) | 0.023387 | 0.023387 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01154 | 4.01154 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.023387 | 0.023387 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE1_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 4.72222e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE2_A: mean (A) | 5.477 | 5.477 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE2_A: sd (A) | 0.020665 | 0.020665 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.477 | 5.477 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.020665 | 0.020665 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE2_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 3.714e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE_min_A: mean (A) | 4.01154 | 4.01154 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE_min_A: sd (A) | 0.023387 | 0.023387 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01154 | 4.01154 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.023387 | 0.023387 | 0 | 0.001 | second formulas | passed |
| test/N1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 4.72222e-07 | 0 | 0.001 | second formulas | passed |
| test/N1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/N1 | nearest water to the zinc, mean (A) | 3.56415 | 3.56415 | 0 | 0.001 | construction | passed |
| test/N1 | water molecules within 3.5 A of the aniline N, mean | 1.48781 | 1.48781 | 0 | 1e-06 | construction | passed |
| test/N1 | frame 1, lig_rmsd_A against starts.csv (N, rmsd_to_crystal_A) (A) | 0.61 | 0.6104 | 0.0004 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, d_Zn_N2_A against starts.csv (N, Zn_N2) (A) | 2.207 | 2.2071 | 0.0001 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, d_Zn_Osulfonyl_min_A against starts.csv (N, Zn_O_min) (A) | 2.775 | 2.7751 | 0.0001 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, d_N2_T199OG1_A against starts.csv (N, N2_OG1) (A) | 2.552 | 2.5522 | 0.0002 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, ang_N2_HN2_T199OG1_deg against starts.csv (N, N2_H_OG1) (deg) | 149.15 | 149.113 | 0.037 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, d_T199N_Osulfonyl_min_A against starts.csv (N, T199N_O_min) (A) | 2.983 | 2.9835 | 0.0005 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, tors_C3_C4_S_N2_deg against starts.csv (N, C3_C4_S_N2) (deg) | 148.632 | 148.632 | 0 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, tors_C5_C4_S_N2_deg against starts.csv (N, C5_C4_S_N2) (deg) | -28.91 | -28.91 | 0 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, tau_fold_deg against starts.csv (N, tau_fold) (deg) | 148.632 | 148.632 | 0 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, ang_C4_S_N2_deg against starts.csv (N, C4_S_N2) (deg) | 109.483 | 109.482 | 0.001 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, d_S_N2_A against starts.csv (N, S_N2) (A) | 1.677 | 1.6775 | 0.0005 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, ang_S_N2_Zn_deg against starts.csv (N, S_N2_Zn) (deg) | 105.908 | 105.908 | 0 | 0.06 | starts.csv | passed |
| test/N1 | frame 1, d1_N1_F131CZ_A against starts.csv (N, N1_Phe131CZ) (A) | 3.356 | 3.3561 | 0.0001 | 0.0015 | starts.csv | passed |
| test/N1 | frame 1, d2_N1_N67OD1_A against starts.csv (N, N1_Asn67OD1) (A) | 7.725 | 7.7255 | 0.0005 | 0.0015 | starts.csv | passed |
| test/N1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| test/N1 | planned length (X0: 5 ns by the specification) (ns) | 100 | 100 | 0 | 1e-09 | construction | passed |
| test/N1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| test/N1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| test/N1 | gate G-V: the criteria that are not met | GV3 Zn-N2, GV3 fifth ligand, full trajectory | GV3 Zn-N2, GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| test/N1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| test/N1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| test/N1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/N1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| test/N1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| test/N1 | gate G-V: provisional (run not finished) | yes | yes | – | equal | construction | passed |
| test/N1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 7.56654 | 7.56653 | 3e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -3.52282 | -3.52282 | 0 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.3486 | -7.34861 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.861245 | 0.861247 | 2e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.986639 | 0.98664 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.627846 | 0.627847 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.682512 | 0.682511 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 16.5317 | 16.5317 | 1e-06 | 0.01 | second formulas | passed |
| test/N1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | gate G-V: solute frames counted for that fraction | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 12 | 12 | 0 | 0 | second formulas | passed |
| test/N1 | gate G-V: full frames counted for that fraction | 123 | 123 | 0 | 0 | construction | passed |
| test/N1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 7.56654 | 7.56653 | 3e-06 | 0.01 | second formulas | passed |
| test/N1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/N1 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/N1 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| test/N1 | N2: frames with positive, negative and zero signed volume | [988, 246, 0] | [988, 246, 0] | – | equal | second formulas | passed |
| test/N1 | N2: fraction of frames with positive signed volume | 0.800648 | 0.800648 | 0 | 1e-06 | second formulas | passed |
| test/N1 | N2: fraction of frames with negative signed volume | 0.199352 | 0.199352 | 0 | 1e-06 | second formulas | passed |
| test/N1 | N2: frames of the window (last 20 ns) | 1234 | 1234 | 0 | 0 | construction | passed |
| test/N1 | N2: positive, negative and zero in the last 20 ns | [988, 246, 0] | [988, 246, 0] | – | equal | second formulas | passed |
| test/N1 | N2: fraction positive in the last 20 ns | 0.800648 | 0.800648 | 0 | 1e-06 | second formulas | passed |
| test/N1 | N2, frames of state C: number, positive, negative, zero | [1234, 988, 246, 0] | [1234, 988, 246, 0] | – | equal | second formulas | passed |
| test/N1 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.800648 | 0.800648 | 0 | 1e-06 | second formulas | passed |
| test/N1 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/N1 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/N1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/N1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/N1 | N2: stretches of one sign, and those of at least 1 ns | [493, 0] | [493, 0] | – | equal | second formulas | passed |
| test/N1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/N1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/N1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/N1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/N1 | N2: sign of the first and of the last frame | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| test/N1 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| test/N1 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| test/N1 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| test/N1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| test/N1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/N1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| test/N1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| test/N1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01154 | 4.01154 | 0 | 0.001 | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 1234] | [0, 1234] | – | equal | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.477 | 5.477 | 0 | 0.001 | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| test/N1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 1234] | [0, 1234] | – | equal | second formulas | passed |
| test/N1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.011543, 4.003722, 0.0] | [4.011543, 4.003722, 0.0] | – | equal | second formulas | passed |
| test/N1 | water at the crystal ring positions: mean number within 3.5 A | 1.48781 | 1.48781 | 0 | 1e-06 | second formulas | passed |
| test/N1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 2.32712 | 2.32711 | 9e-06 | 0.001 | second formulas | passed |
| test/N1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.24566 | 3.24566 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.48781 | 1.48781 | 0 | 1e-06 | second formulas | passed |
| test/N1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/N1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/N1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/N1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/N1 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.48668 | 3.48667 | 1e-05 | 0.001 | second formulas | passed |
| test/N1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.39285 | 4.39285 | 1e-06 | 0.001 | second formulas | passed |
| test/N1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.486674] | [0.0, 3.486674] | – | equal | second formulas | passed |
| test/N1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.32711 | 2.32711 | – | equal | second formulas | passed |
| test/N1 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| test/N1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/N1 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| test/N1 | first production frame: signed volume at N2 | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| test/N1 | first production frame: ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| test/N1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| test/N1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.800648, 0.800648, 0] | [0.800648, 0.800648, 0] | – | equal | second formulas | passed |
| test/N1 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/N1 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/N1 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | solute frames seen | 3000 | 3000 | 0 | 0 | construction | passed |
| test/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| test/R1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| test/R1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| test/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/R1 | where the system folder was found | convention (the system_dir of the manifest is not on this machine) | convention (the system_dir of the manifest is not on this machine) | – | equal | construction | passed |
| test/R1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/R1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/R1 | full frames seen | 300 | 300 | 0 | 0 | construction | passed |
| test/R1 | time seen (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/R1 | frames in which the ligand was in another periodic image | 429 | 429 | 0 | 0 | construction | passed |
| test/R1 | frames in which the zinc was in another periodic image | 273 | 273 | 0 | 0 | construction | passed |
| test/R1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | lines of timeseries_solute.csv | 3000 | 3000 | 0 | 0 | construction | passed |
| test/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| test/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000963 | 0.000963 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000808 | 0.000808 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000847 | 0.000847 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| test/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000763 | 0.000763 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| test/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/R1 | every frame: n2_volume, largest difference | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| test/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | lines of timeseries_full.csv | 300 | 300 | 0 | 0 | construction | passed |
| test/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| test/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 4.3e-05 | 4.3e-05 | 0.001 | second formulas | passed |
| test/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | verdict | CONVERTED | CONVERTED | – | equal | construction | passed |
| test/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| test/R1 | run complete | yes | yes | – | equal | construction | passed |
| test/R1 | first passage into the opposite state (ns) | 8.01 | 8.01 | 0 | 1e-09 | construction | passed |
| test/R1 | transitions | 1 | 1 | 0 | 0 | construction | passed |
| test/R1 | long visits (at least 1 ns) | 2 | 2 | 0 | 0 | construction | passed |
| test/R1 | frames of the verdict window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/R1 | count_C_window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/R1 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | count_C | 2200 | 2200 | 0 | 0 | construction | passed |
| test/R1 | count_R | 800 | 800 | 0 | 0 | construction | passed |
| test/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| test/R1 | count_C (design against the planted frames) | 2200 | 2200 | 0 | 0 | construction | passed |
| test/R1 | count_R (design against the planted frames) | 800 | 800 | 0 | 0 | construction | passed |
| test/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/R1 | frames with the hydrogen bond N2-H...OG1 | 2400 | 2400 | 0 | 0 | second formulas | passed |
| test/R1 | full frames with a water within 2.6 A of the zinc | 30 | 30 | 0 | 0 | construction | passed |
| test/R1 | frames with labelling identity | 1350 | 1350 | 0 | 0 | second formulas | passed |
| test/R1 | frames with labelling ring_flip | 550 | 550 | 0 | 0 | second formulas | passed |
| test/R1 | frames with labelling oxygen_swap | 550 | 550 | 0 | 0 | second formulas | passed |
| test/R1 | frames with labelling ring_flip+oxygen_swap | 550 | 550 | 0 | 0 | second formulas | passed |
| test/R1 | lig_rmsd_A: mean (A) | 0.95286 | 0.952861 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | lig_rmsd_A: sd (A) | 1.58014 | 1.58014 | 0 | 0.001 | second formulas | passed |
| test/R1 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/R1 | lig_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_N2_A: mean (A) | 2.04243 | 2.04243 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_N2_A: sd (A) | 0.161671 | 0.161671 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 3.79398e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | tau_fold_deg: mean (deg) | 134.906 | 134.906 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | tau_fold_deg: sd (deg) | 4.72659 | 4.72659 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | tau_fold_deg: sd_last_half (deg) | 0 | 7.2e-05 | 7.2e-05 | 0.01 | second formulas | passed |
| test/R1 | tau_fold_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 3.30303e-06 | 3e-06 | 0.01 | second formulas | passed |
| test/R1 | ang_C4_S_N2_deg: mean (deg) | 108.965 | 108.965 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | ang_C4_S_N2_deg: sd (deg) | 0.281836 | 0.281835 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 6e-05 | 6e-05 | 0.01 | second formulas | passed |
| test/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 2.63153e-06 | 3e-06 | 0.01 | second formulas | passed |
| test/R1 | d_S_N2_A: mean (A) | 1.61791 | 1.61791 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_S_N2_A: sd (A) | 0.035783 | 0.035783 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | d_S_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 3.27904e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/R1 | ca_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | d1_N1_F131CZ_A: mean (A) | 5.41349 | 5.41349 | 0 | 0.001 | second formulas | passed |
| test/R1 | d1_N1_F131CZ_A: sd (A) | 2.51548 | 2.51548 | 0 | 0.001 | second formulas | passed |
| test/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| test/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 3.06539e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | d2_N1_N67OD1_A: mean (A) | 6.37164 | 6.37164 | 0 | 0.001 | second formulas | passed |
| test/R1 | d2_N1_N67OD1_A: sd (A) | 2.27632 | 2.27632 | 0 | 0.001 | second formulas | passed |
| test/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| test/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 2.48237e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | n2_volume: mean | 0.203448 | 0.203448 | 0 | 0.001 | second formulas | passed |
| test/R1 | n2_volume: sd | 0.666698 | 0.666698 | 0 | 0.001 | second formulas | passed |
| test/R1 | n2_volume: mean_last_half | 0.451372 | 0.451373 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| test/R1 | n2_volume: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | n2_volume: block means, largest difference | 0 | 5.20361e-07 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE1_A: mean (A) | 4.00856 | 4.00856 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE1_A: sd (A) | 0.01521 | 0.01521 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 3.25403e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE2_A: mean (A) | 5.47469 | 5.47469 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE2_A: sd (A) | 0.014189 | 0.014189 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 3.27302e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE_min_A: mean (A) | 4.00856 | 4.00856 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE_min_A: sd (A) | 0.01521 | 0.01521 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 3.25403e-07 | 0 | 0.001 | second formulas | passed |
| test/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/R1 | nearest water to the zinc, mean (A) | 3.561 | 3.561 | 0 | 0.001 | construction | passed |
| test/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| test/R1 | frame 1, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| test/R1 | frame 1, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| test/R1 | frame 1, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| test/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| test/R1 | planned length (X0: 5 ns by the specification) (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| test/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| test/R1 | gate G-V: the criteria that are not met | GV3 Zn-N2, GV3 fifth ligand, full trajectory | GV3 Zn-N2, GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| test/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| test/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| test/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| test/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| test/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| test/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 1.16175 | 1.16175 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -1.29869 | -1.29869 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -2.47232 | -2.47232 | 0 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.764976 | 0.764975 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.938973 | 0.938974 | 1e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.76465 | 0.76465 | 0 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.164552 | 0.16455 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 2.80631 | 2.80631 | 2e-06 | 0.01 | second formulas | passed |
| test/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | gate G-V: solute frames counted for that fraction | 3000 | 3000 | 0 | 0 | construction | passed |
| test/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 30 | 30 | 0 | 0 | second formulas | passed |
| test/R1 | gate G-V: full frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| test/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 2.47232 | 2.47232 | 0 | 0.01 | second formulas | passed |
| test/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/R1 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| test/R1 | N2: frames with positive, negative and zero signed volume | [1920, 1080, 0] | [1920, 1080, 0] | – | equal | second formulas | passed |
| test/R1 | N2: fraction of frames with positive signed volume | 0.64 | 0.64 | 0 | 1e-06 | second formulas | passed |
| test/R1 | N2: fraction of frames with negative signed volume | 0.36 | 0.36 | 0 | 1e-06 | second formulas | passed |
| test/R1 | N2: frames of the window (last 20 ns) | 2000 | 2000 | 0 | 0 | construction | passed |
| test/R1 | N2: positive, negative and zero in the last 20 ns | [1600, 400, 0] | [1600, 400, 0] | – | equal | second formulas | passed |
| test/R1 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/R1 | N2, frames of state C: number, positive, negative, zero | [2200, 1760, 440, 0] | [2200, 1760, 440, 0] | – | equal | second formulas | passed |
| test/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/R1 | N2, frames of state R: number, positive, negative, zero | [800, 160, 640, 0] | [800, 160, 640, 0] | – | equal | second formulas | passed |
| test/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| test/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/R1 | N2: stretches of one sign, and those of at least 1 ns | [1199, 0] | [1199, 0] | – | equal | second formulas | passed |
| test/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/R1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/R1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/R1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/R1 | N2: sign of the first and of the last frame | [-1, -1] | [-1, -1] | – | equal | second formulas | passed |
| test/R1 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| test/R1 | hydrogen bond around the first passage: frame of the first passage | 800 | 800 | – | equal | second formulas | passed |
| test/R1 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [800, 640] | [800, 640] | – | equal | second formulas | passed |
| test/R1 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [2200, 1760] | [2200, 1760] | – | equal | second formulas | passed |
| test/R1 | hydrogen bond N2-H...OG1, fraction before the first passage | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/R1 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.8, 0.8] | [0.8, 0.8] | – | equal | construction | passed |
| test/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/R1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/R1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| test/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00856 | 4.00856 | 0 | 0.001 | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47469 | 5.47469 | 0 | 0.001 | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| test/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.00856, 4.00372, 0.0] | [4.00856, 4.00372, 0.0] | – | equal | second formulas | passed |
| test/R1 | water at the crystal ring positions: mean number within 3.5 A | 1.1 | 1.1 | 0 | 1e-06 | second formulas | passed |
| test/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.3e-05 | 0.001 | second formulas | passed |
| test/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 4.72609 | 4.72609 | 1e-06 | 0.001 | second formulas | passed |
| test/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.1 | 1.1 | 0 | 1e-06 | second formulas | passed |
| test/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/R1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/R1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1e-05 | 0.001 | second formulas | passed |
| test/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 5.53234 | 5.53234 | 0 | 0.001 | second formulas | passed |
| test/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020123] | [0.0, 4.020123] | – | equal | second formulas | passed |
| test/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| test/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| test/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| test/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| test/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| test/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| test/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.64, 0.8, 0] | [0.64, 0.8, 0] | – | equal | second formulas | passed |
| test/R1 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/R1 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/R1 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| test/X0 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| test/X0 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| test/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/X0 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/X0 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| test/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| test/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| test/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| test/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| test/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| test/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000104 | 0.000104 | 0.01 | second formulas | passed |
| test/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| test/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| test/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| test/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| test/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X0 | run complete | yes | yes | – | equal | construction | passed |
| test/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| test/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| test/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| test/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| test/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| test/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| test/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| test/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| test/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| test/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | tau_fold_deg: sd (deg) | 0 | 6.6e-05 | 6.6e-05 | 0.01 | second formulas | passed |
| test/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 6.1e-05 | 6.1e-05 | 0.01 | second formulas | passed |
| test/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| test/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.7e-05 | 5.7e-05 | 0.01 | second formulas | passed |
| test/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.8e-05 | 5.8e-05 | 0.01 | second formulas | passed |
| test/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| test/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| test/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| test/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| test/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| test/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| test/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| test/X0 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| test/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| test/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 0 | 0.001 | construction | passed |
| test/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| test/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| test/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| test/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| test/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| test/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| test/X0 | gate G-V: the criteria that are not met | none | none | – | equal | second formulas | passed |
| test/X0 | gate G-V: all criteria met | yes | yes | – | equal | second formulas | passed |
| test/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| test/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| test/X0 | gate G-V: verdict (none for a run that is not X0) | PASSED | PASSED | – | equal | second formulas | passed |
| test/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823962 | 2e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840542 | 3e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776332 | 2e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731084 | 3e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807202 | 3e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005403 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043932 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| test/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| test/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| test/X0 | gate G-V by design: verdict | PASSED | PASSED | – | equal | construction | passed |
| test/X0 | gate G-V by design: the criteria that are not met | none | none | – | equal | construction | passed |
| test/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| test/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| test/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| test/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/X0 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| test/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| test/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| test/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| test/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| test/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| test/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| test/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| test/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| test/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| test/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| test/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| test/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| test/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| test/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| test/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| test/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003722, 0.0] | [4.007685, 4.003722, 0.0] | – | equal | second formulas | passed |
| test/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| test/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.1e-05 | 0.001 | second formulas | passed |
| test/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.78 | 3.78 | 2e-06 | 0.001 | second formulas | passed |
| test/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| test/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.1e-05 | 0.001 | second formulas | passed |
| test/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.83703 | 4.83703 | 0 | 0.001 | second formulas | passed |
| test/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020122] | [0.0, 4.020122] | – | equal | second formulas | passed |
| test/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| test/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| test/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| test/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| test/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| test/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| test/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | solute frames seen | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X1 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| test/X1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| test/X1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| test/X1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/X1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/X1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/X1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/X1 | full frames seen | 300 | 300 | 0 | 0 | construction | passed |
| test/X1 | time seen (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/X1 | frames in which the ligand was in another periodic image | 429 | 429 | 0 | 0 | construction | passed |
| test/X1 | frames in which the zinc was in another periodic image | 273 | 273 | 0 | 0 | construction | passed |
| test/X1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | lines of timeseries_solute.csv | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/X1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 3.9e-05 | 3.9e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001128 | 0.001128 | 0.01 | second formulas | passed |
| test/X1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.4e-05 | 4.4e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| test/X1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000588 | 0.000588 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000641 | 0.000641 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000791 | 0.000791 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| test/X1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/X1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| test/X1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| test/X1 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X1 | every frame: n2_volume, largest difference | 0 | 1.7e-05 | 1.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | lines of timeseries_full.csv | 300 | 300 | 0 | 0 | construction | passed |
| test/X1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/X1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/X1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/X1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | verdict | CONVERTED | CONVERTED | – | equal | construction | passed |
| test/X1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X1 | run complete | yes | yes | – | equal | construction | passed |
| test/X1 | first passage into the opposite state (ns) | 9.01 | 9.01 | 0 | 1e-09 | construction | passed |
| test/X1 | transitions | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | long visits (at least 1 ns) | 5 | 5 | 0 | 0 | construction | passed |
| test/X1 | frames of the verdict window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X1 | count_C_window | 200 | 200 | 0 | 0 | construction | passed |
| test/X1 | count_R_window | 1800 | 1800 | 0 | 0 | construction | passed |
| test/X1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X1 | count_C | 1050 | 1050 | 0 | 0 | construction | passed |
| test/X1 | count_R | 1930 | 1930 | 0 | 0 | construction | passed |
| test/X1 | count_I | 20 | 20 | 0 | 0 | construction | passed |
| test/X1 | rule 1: left the crystal pose | yes | yes | – | equal | construction | passed |
| test/X1 | count_C (design against the planted frames) | 1050 | 1050 | 0 | 0 | construction | passed |
| test/X1 | count_R (design against the planted frames) | 1930 | 1930 | 0 | 0 | construction | passed |
| test/X1 | count_I (design against the planted frames) | 20 | 20 | 0 | 0 | construction | passed |
| test/X1 | frames with the hydrogen bond N2-H...OG1 | 2384 | 2384 | 0 | 0 | second formulas | passed |
| test/X1 | full frames with a water within 2.6 A of the zinc | 30 | 30 | 0 | 0 | construction | passed |
| test/X1 | frames with labelling identity | 1247 | 1247 | 0 | 0 | second formulas | passed |
| test/X1 | frames with labelling ring_flip | 1227 | 1227 | 0 | 0 | second formulas | passed |
| test/X1 | frames with labelling oxygen_swap | 263 | 263 | 0 | 0 | second formulas | passed |
| test/X1 | frames with labelling ring_flip+oxygen_swap | 263 | 263 | 0 | 0 | second formulas | passed |
| test/X1 | lig_rmsd_A: mean (A) | 2.31128 | 2.31128 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | lig_rmsd_A: sd (A) | 1.69663 | 1.69663 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | lig_rmsd_A: mean_last_half (A) | 3.32413 | 3.32413 | 0 | 0.001 | second formulas | passed |
| test/X1 | lig_rmsd_A: sd_last_half (A) | 0.888411 | 0.888411 | 0 | 0.001 | second formulas | passed |
| test/X1 | lig_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | lig_rmsd_A: block means, largest difference (A) | 0 | 1.15581e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_N2_A: mean (A) | 2.14124 | 2.14124 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_N2_A: sd (A) | 0.1731 | 0.1731 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_N2_A: mean_last_half (A) | 2.21175 | 2.21175 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_N2_A: sd_last_half (A) | 0.11214 | 0.11214 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 5.17667e-07 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X1 | tau_fold_deg: sd (deg) | 0 | 7.2e-05 | 7.2e-05 | 0.01 | second formulas | passed |
| test/X1 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X1 | tau_fold_deg: sd_last_half (deg) | 0 | 7.2e-05 | 7.2e-05 | 0.01 | second formulas | passed |
| test/X1 | tau_fold_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | tau_fold_deg: block means, largest difference (deg) | 0 | 5.04222e-06 | 5e-06 | 0.01 | second formulas | passed |
| test/X1 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 0 | 0.01 | second formulas | passed |
| test/X1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.8e-05 | 5.8e-05 | 0.01 | second formulas | passed |
| test/X1 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| test/X1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.9e-05 | 5.9e-05 | 0.01 | second formulas | passed |
| test/X1 | ang_C4_S_N2_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 2.49546e-06 | 2e-06 | 0.01 | second formulas | passed |
| test/X1 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | d_S_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d_S_N2_A: block means, largest difference (A) | 0 | 3.27904e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X1 | ca_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | d1_N1_F131CZ_A: mean (A) | 7.55731 | 7.55731 | 0 | 0.001 | second formulas | passed |
| test/X1 | d1_N1_F131CZ_A: sd (A) | 2.6949 | 2.6949 | 0 | 0.001 | second formulas | passed |
| test/X1 | d1_N1_F131CZ_A: mean_last_half (A) | 9.17846 | 9.17846 | 0 | 0.001 | second formulas | passed |
| test/X1 | d1_N1_F131CZ_A: sd_last_half (A) | 1.41164 | 1.41164 | 0 | 0.001 | second formulas | passed |
| test/X1 | d1_N1_F131CZ_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 4.37034e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | d2_N1_N67OD1_A: mean (A) | 4.45022 | 4.45022 | 0 | 0.001 | second formulas | passed |
| test/X1 | d2_N1_N67OD1_A: sd (A) | 2.44001 | 2.44001 | 0 | 0.001 | second formulas | passed |
| test/X1 | d2_N1_N67OD1_A: mean_last_half (A) | 2.97511 | 2.97511 | 0 | 0.001 | second formulas | passed |
| test/X1 | d2_N1_N67OD1_A: sd_last_half (A) | 1.27462 | 1.27462 | 0 | 0.001 | second formulas | passed |
| test/X1 | d2_N1_N67OD1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 3.96598e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | n2_volume: mean | -0.013422 | -0.013422 | 0 | 0.001 | second formulas | passed |
| test/X1 | n2_volume: sd | 0.576491 | 0.576491 | 0 | 0.001 | second formulas | passed |
| test/X1 | n2_volume: mean_last_half | -0.221972 | -0.221972 | 0 | 0.001 | second formulas | passed |
| test/X1 | n2_volume: sd_last_half | 0.457509 | 0.457509 | 0 | 0.001 | second formulas | passed |
| test/X1 | n2_volume: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | n2_volume: block means, largest difference | 0 | 1.72083e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE1_A: mean (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE1_A: sd (A) | 0.016208 | 0.016208 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00932 | 4.00932 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.016848 | 0.016848 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 4.36552e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE2_A: mean (A) | 5.4745 | 5.4745 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE2_A: sd (A) | 0.014321 | 0.014321 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47449 | 5.47449 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.013193 | 0.013193 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 1.24432e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE_min_A: mean (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE_min_A: sd (A) | 0.016208 | 0.016208 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00932 | 4.00932 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.016848 | 0.016848 | 0 | 0.001 | second formulas | passed |
| test/X1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 4.36552e-07 | 0 | 0.001 | second formulas | passed |
| test/X1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/X1 | nearest water to the zinc, mean (A) | 3.561 | 3.561 | 0 | 0.001 | construction | passed |
| test/X1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| test/X1 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| test/X1 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X1 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| test/X1 | planned length (X0: 5 ns by the specification) (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| test/X1 | rule 1 applies to this run (production runs of start X) | yes | yes | – | equal | construction | passed |
| test/X1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | GV1 ligand RMSD, GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| test/X1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| test/X1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| test/X1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/X1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| test/X1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| test/X1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 2.64267 | 2.64267 | 2e-06 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.118427 | 0.118427 | 0 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -5.58019 | -5.58018 | 2e-06 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.831615 | 0.831615 | 0 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.894684 | 0.894684 | 0 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.731171 | 0.731172 | 1e-06 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005402 | 0 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.04393 | 1e-06 | 0.01 | second formulas | passed |
| test/X1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 20 | 20 | 0 | 0 | second formulas | passed |
| test/X1 | gate G-V: solute frames counted for that fraction | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 32 | 32 | 0 | 0 | second formulas | passed |
| test/X1 | gate G-V: full frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| test/X1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 5.58019 | 5.58018 | 2e-06 | 0.01 | second formulas | passed |
| test/X1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/X1 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/X1 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| test/X1 | N2: frames with positive, negative and zero signed volume | [1242, 1758, 0] | [1242, 1758, 0] | – | equal | second formulas | passed |
| test/X1 | N2: fraction of frames with positive signed volume | 0.414 | 0.414 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2: fraction of frames with negative signed volume | 0.586 | 0.586 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2: frames of the window (last 20 ns) | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X1 | N2: positive, negative and zero in the last 20 ns | [520, 1480, 0] | [520, 1480, 0] | – | equal | second formulas | passed |
| test/X1 | N2: fraction positive in the last 20 ns | 0.26 | 0.26 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2, frames of state C: number, positive, negative, zero | [1050, 840, 210, 0] | [1050, 840, 210, 0] | – | equal | second formulas | passed |
| test/X1 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2, frames of state R: number, positive, negative, zero | [1930, 386, 1544, 0] | [1930, 386, 1544, 0] | – | equal | second formulas | passed |
| test/X1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2, frames of state I: number, positive, negative, zero | [20, 16, 4, 0] | [20, 16, 4, 0] | – | equal | second formulas | passed |
| test/X1 | N2, frames of state I: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X1 | N2: stretches of one sign, and those of at least 1 ns | [1193, 0] | [1193, 0] | – | equal | second formulas | passed |
| test/X1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/X1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/X1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/X1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/X1 | N2: sign of the first and of the last frame | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| test/X1 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| test/X1 | hydrogen bond around the first passage: frame of the first passage | 900 | 900 | – | equal | second formulas | passed |
| test/X1 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [900, 720] | [900, 720] | – | equal | second formulas | passed |
| test/X1 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [2100, 1664] | [2100, 1664] | – | equal | second formulas | passed |
| test/X1 | hydrogen bond N2-H...OG1, fraction before the first passage | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X1 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.792381 | 0.792381 | 0 | 1e-06 | second formulas | passed |
| test/X1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.8, 0.792381] | [0.8, 0.792381] | – | equal | construction | passed |
| test/X1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/X1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/X1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 7e-06 | 0.001 | second formulas | passed |
| test/X1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.4745 | 5.4745 | 0 | 0.001 | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 7e-06 | 0.001 | second formulas | passed |
| test/X1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.008917, 4.00372, 0.0] | [4.008917, 4.00372, 0.0] | – | equal | second formulas | passed |
| test/X1 | water at the crystal ring positions: mean number within 3.5 A | 0.516667 | 0.516667 | 0 | 1e-06 | second formulas | passed |
| test/X1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.3e-05 | 0.001 | second formulas | passed |
| test/X1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 6.10123 | 6.10123 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 0.516667 | 0.516667 | 0 | 1e-06 | second formulas | passed |
| test/X1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/X1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/X1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X1 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.7e-05 | 0.001 | second formulas | passed |
| test/X1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 6.55991 | 6.55991 | 1e-06 | 0.001 | second formulas | passed |
| test/X1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020116] | [0.0, 4.020116] | – | equal | second formulas | passed |
| test/X1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| test/X1 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| test/X1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/X1 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| test/X1 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| test/X1 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| test/X1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.414, 0.26, 0] | [0.414, 0.26, 0] | – | equal | second formulas | passed |
| test/X1 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/X1 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/X1 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | solute frames seen | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X2 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| test/X2 | manifest lists the segments | no | no | – | equal | construction | passed |
| test/X2 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| test/X2 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/X2 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/X2 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| test/X2 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| test/X2 | full frames seen | 300 | 300 | 0 | 0 | construction | passed |
| test/X2 | time seen (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X2 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/X2 | frames in which the ligand was in another periodic image | 429 | 429 | 0 | 0 | construction | passed |
| test/X2 | frames in which the zinc was in another periodic image | 273 | 273 | 0 | 0 | construction | passed |
| test/X2 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | lines of timeseries_solute.csv | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X2 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/X2 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001128 | 0.001128 | 0.01 | second formulas | passed |
| test/X2 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| test/X2 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000847 | 0.000847 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000755 | 0.000755 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| test/X2 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X2 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/X2 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X2 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X2 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X2 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| test/X2 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| test/X2 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | lines of timeseries_full.csv | 300 | 300 | 0 | 0 | construction | passed |
| test/X2 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/X2 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/X2 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X2 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X2 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/X2 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | verdict | STAYED | STAYED | – | equal | construction | passed |
| test/X2 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X2 | run complete | yes | yes | – | equal | construction | passed |
| test/X2 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| test/X2 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| test/X2 | frames of the verdict window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X2 | count_C_window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X2 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | count_C | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X2 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| test/X2 | count_C (design against the planted frames) | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X2 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| test/X2 | frames with the hydrogen bond N2-H...OG1 | 2400 | 2400 | 0 | 0 | second formulas | passed |
| test/X2 | full frames with a water within 2.6 A of the zinc | 30 | 30 | 0 | 0 | construction | passed |
| test/X2 | frames with labelling identity | 750 | 750 | 0 | 0 | second formulas | passed |
| test/X2 | frames with labelling ring_flip | 750 | 750 | 0 | 0 | second formulas | passed |
| test/X2 | frames with labelling oxygen_swap | 750 | 750 | 0 | 0 | second formulas | passed |
| test/X2 | frames with labelling ring_flip+oxygen_swap | 750 | 750 | 0 | 0 | second formulas | passed |
| test/X2 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | lig_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | lig_rmsd_A: block means, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d_Zn_N2_A: block means, largest difference (A) | 0 | 2.07264e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 2e-06 | 0.01 | second formulas | passed |
| test/X2 | tau_fold_deg: sd (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| test/X2 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 3e-06 | 0.01 | second formulas | passed |
| test/X2 | tau_fold_deg: sd_last_half (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| test/X2 | tau_fold_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | tau_fold_deg: block means, largest difference (deg) | 0 | 6.04222e-06 | 6e-06 | 0.01 | second formulas | passed |
| test/X2 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 0 | 0.01 | second formulas | passed |
| test/X2 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.9e-05 | 5.9e-05 | 0.01 | second formulas | passed |
| test/X2 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| test/X2 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.6e-05 | 5.6e-05 | 0.01 | second formulas | passed |
| test/X2 | ang_C4_S_N2_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 3.49546e-06 | 3e-06 | 0.01 | second formulas | passed |
| test/X2 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | d_S_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d_S_N2_A: block means, largest difference (A) | 0 | 3.27904e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X2 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X2 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | ca_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| test/X2 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| test/X2 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | d1_N1_F131CZ_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 1.84776e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| test/X2 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| test/X2 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | d2_N1_N67OD1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 1.18924e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| test/X2 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| test/X2 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| test/X2 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| test/X2 | n2_volume: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | n2_volume: block means, largest difference | 0 | 4.79639e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 1.47441e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 3.27302e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| test/X2 | d_Zn_E106OE_min_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X2 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 1.47441e-07 | 0 | 0.001 | second formulas | passed |
| test/X2 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/X2 | nearest water to the zinc, mean (A) | 3.561 | 3.561 | 0 | 0.001 | construction | passed |
| test/X2 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| test/X2 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| test/X2 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X2 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X2 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| test/X2 | planned length (X0: 5 ns by the specification) (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X2 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| test/X2 | rule 1 applies to this run (production runs of start X) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V: the criteria that are not met | GV3 fifth ligand, full trajectory | GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| test/X2 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| test/X2 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| test/X2 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/X2 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| test/X2 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| test/X2 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X2 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823962 | 2e-06 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840546 | 1e-06 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776335 | 1e-06 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731087 | 0 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930668 | 1e-06 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807205 | 0 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005402 | 0 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043929 | 2e-06 | 0.01 | second formulas | passed |
| test/X2 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | gate G-V: solute frames counted for that fraction | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X2 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 30 | 30 | 0 | 0 | second formulas | passed |
| test/X2 | gate G-V: full frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| test/X2 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930668 | 1e-06 | 0.01 | second formulas | passed |
| test/X2 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/X2 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/X2 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| test/X2 | N2: frames with positive, negative and zero signed volume | [2400, 600, 0] | [2400, 600, 0] | – | equal | second formulas | passed |
| test/X2 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X2 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| test/X2 | N2: frames of the window (last 20 ns) | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X2 | N2: positive, negative and zero in the last 20 ns | [1600, 400, 0] | [1600, 400, 0] | – | equal | second formulas | passed |
| test/X2 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X2 | N2, frames of state C: number, positive, negative, zero | [3000, 2400, 600, 0] | [3000, 2400, 600, 0] | – | equal | second formulas | passed |
| test/X2 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X2 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/X2 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/X2 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| test/X2 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| test/X2 | N2: stretches of one sign, and those of at least 1 ns | [1200, 0] | [1200, 0] | – | equal | second formulas | passed |
| test/X2 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/X2 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/X2 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/X2 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/X2 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| test/X2 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| test/X2 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| test/X2 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| test/X2 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| test/X2 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/X2 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/X2 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| test/X2 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 7e-06 | 0.001 | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| test/X2 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X2 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.00372, 0.0] | [4.007685, 4.00372, 0.0] | – | equal | second formulas | passed |
| test/X2 | water at the crystal ring positions: mean number within 3.5 A | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| test/X2 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.2e-05 | 0.001 | second formulas | passed |
| test/X2 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.75 | 3.75 | 2e-06 | 0.001 | second formulas | passed |
| test/X2 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| test/X2 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/X2 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/X2 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X2 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X2 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.2e-05 | 0.001 | second formulas | passed |
| test/X2 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.80561 | 4.80561 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020121] | [0.0, 4.020121] | – | equal | second formulas | passed |
| test/X2 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| test/X2 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| test/X2 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/X2 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| test/X2 | first production frame: signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| test/X2 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X2 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| test/X2 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| test/X2 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/X2 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/X2 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | solute frames seen | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X3 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| test/X3 | manifest lists the segments | no | no | – | equal | construction | passed |
| test/X3 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| test/X3 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/X3 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/X3 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| test/X3 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| test/X3 | full frames seen | 300 | 300 | 0 | 0 | construction | passed |
| test/X3 | time seen (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X3 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| test/X3 | frames in which the ligand was in another periodic image | 429 | 429 | 0 | 0 | construction | passed |
| test/X3 | frames in which the zinc was in another periodic image | 273 | 273 | 0 | 0 | construction | passed |
| test/X3 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | lines of timeseries_solute.csv | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X3 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | no | no | – | equal | construction | passed |
| test/X3 | every frame: lig_rmsd_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 6.1e-05 | 6.1e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001057 | 0.001057 | 0.01 | second formulas | passed |
| test/X3 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.4e-05 | 4.4e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| test/X3 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000588 | 0.000588 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000641 | 0.000641 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000606 | 0.000606 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000689 | 0.000689 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| test/X3 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| test/X3 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X3 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| test/X3 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X3 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X3 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| test/X3 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| test/X3 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| test/X3 | every frame: n2_volume, largest difference | 0 | 1.7e-05 | 1.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| test/X3 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | lines of timeseries_full.csv | 300 | 300 | 0 | 0 | construction | passed |
| test/X3 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/X3 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| test/X3 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| test/X3 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X3 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 5.7e-05 | 5.7e-05 | 0.001 | second formulas | passed |
| test/X3 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| test/X3 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | verdict | MIXED | MIXED | – | equal | construction | passed |
| test/X3 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X3 | run complete | yes | yes | – | equal | construction | passed |
| test/X3 | first passage into the opposite state (ns) | 9.01 | 9.01 | 0 | 1e-09 | construction | passed |
| test/X3 | transitions | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | long visits (at least 1 ns) | 5 | 5 | 0 | 0 | construction | passed |
| test/X3 | frames of the verdict window | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X3 | count_C_window | 201 | 201 | 0 | 0 | construction | passed |
| test/X3 | count_R_window | 1799 | 1799 | 0 | 0 | construction | passed |
| test/X3 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| test/X3 | count_C | 1051 | 1051 | 0 | 0 | construction | passed |
| test/X3 | count_R | 1929 | 1929 | 0 | 0 | construction | passed |
| test/X3 | count_I | 20 | 20 | 0 | 0 | construction | passed |
| test/X3 | rule 1: left the crystal pose | yes | yes | – | equal | construction | passed |
| test/X3 | count_C (design against the planted frames) | 1051 | 1051 | 0 | 0 | construction | passed |
| test/X3 | count_R (design against the planted frames) | 1929 | 1929 | 0 | 0 | construction | passed |
| test/X3 | count_I (design against the planted frames) | 20 | 20 | 0 | 0 | construction | passed |
| test/X3 | frames with the hydrogen bond N2-H...OG1 | 2384 | 2384 | 0 | 0 | second formulas | passed |
| test/X3 | full frames with a water within 2.6 A of the zinc | 30 | 30 | 0 | 0 | construction | passed |
| test/X3 | frames with labelling identity | 1247 | 1247 | 0 | 0 | second formulas | passed |
| test/X3 | frames with labelling ring_flip | 1226 | 1226 | 0 | 0 | second formulas | passed |
| test/X3 | frames with labelling oxygen_swap | 263 | 263 | 0 | 0 | second formulas | passed |
| test/X3 | frames with labelling ring_flip+oxygen_swap | 264 | 264 | 0 | 0 | second formulas | passed |
| test/X3 | lig_rmsd_A: mean (A) | 2.31009 | 2.31009 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | lig_rmsd_A: sd (A) | 1.697 | 1.697 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | lig_rmsd_A: mean_last_half (A) | 3.32176 | 3.32176 | 0 | 0.001 | second formulas | passed |
| test/X3 | lig_rmsd_A: sd_last_half (A) | 0.892523 | 0.892523 | 0 | 0.001 | second formulas | passed |
| test/X3 | lig_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | lig_rmsd_A: block means, largest difference (A) | 0 | 1.15581e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_N2_A: mean (A) | 2.14115 | 2.14115 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_N2_A: sd (A) | 0.17314 | 0.17314 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_N2_A: mean_last_half (A) | 2.21157 | 2.21157 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_N2_A: sd_last_half (A) | 0.112372 | 0.112372 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d_Zn_N2_A: block means, largest difference (A) | 0 | 4.82333e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | tau_fold_deg: sd (deg) | 0 | 7.2e-05 | 7.2e-05 | 0.01 | second formulas | passed |
| test/X3 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | tau_fold_deg: sd_last_half (deg) | 0 | 7.1e-05 | 7.1e-05 | 0.01 | second formulas | passed |
| test/X3 | tau_fold_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | tau_fold_deg: block means, largest difference (deg) | 0 | 1.95778e-06 | 2e-06 | 0.01 | second formulas | passed |
| test/X3 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.9e-05 | 5.9e-05 | 0.01 | second formulas | passed |
| test/X3 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.8e-05 | 5.8e-05 | 0.01 | second formulas | passed |
| test/X3 | ang_C4_S_N2_deg: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 1.50454e-06 | 2e-06 | 0.01 | second formulas | passed |
| test/X3 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | d_S_N2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d_S_N2_A: block means, largest difference (A) | 0 | 3.27904e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X3 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X3 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| test/X3 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X3 | ca_rmsd_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | d1_N1_F131CZ_A: mean (A) | 7.55543 | 7.55543 | 0 | 0.001 | second formulas | passed |
| test/X3 | d1_N1_F131CZ_A: sd (A) | 2.69548 | 2.69548 | 0 | 0.001 | second formulas | passed |
| test/X3 | d1_N1_F131CZ_A: mean_last_half (A) | 9.17469 | 9.17469 | 0 | 0.001 | second formulas | passed |
| test/X3 | d1_N1_F131CZ_A: sd_last_half (A) | 1.41817 | 1.41817 | 0 | 0.001 | second formulas | passed |
| test/X3 | d1_N1_F131CZ_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 4.37034e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | d2_N1_N67OD1_A: mean (A) | 4.45192 | 4.45192 | 0 | 0.001 | second formulas | passed |
| test/X3 | d2_N1_N67OD1_A: sd (A) | 2.44053 | 2.44053 | 0 | 0.001 | second formulas | passed |
| test/X3 | d2_N1_N67OD1_A: mean_last_half (A) | 2.97851 | 2.97851 | 0 | 0.001 | second formulas | passed |
| test/X3 | d2_N1_N67OD1_A: sd_last_half (A) | 1.28052 | 1.28052 | 0 | 0.001 | second formulas | passed |
| test/X3 | d2_N1_N67OD1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 3.94481e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | n2_volume: mean | -0.013818 | -0.013818 | 0 | 0.001 | second formulas | passed |
| test/X3 | n2_volume: sd | 0.576512 | 0.576512 | 0 | 0.001 | second formulas | passed |
| test/X3 | n2_volume: mean_last_half | -0.222765 | -0.222765 | 0 | 0.001 | second formulas | passed |
| test/X3 | n2_volume: sd_last_half | 0.4572 | 0.4572 | 0 | 0.001 | second formulas | passed |
| test/X3 | n2_volume: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | n2_volume: block means, largest difference | 0 | 2.33511e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE1_A: mean (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE1_A: sd (A) | 0.016208 | 0.016208 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00932 | 4.00932 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE1_A: sd_last_half (A) | 0.016848 | 0.016848 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE1_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 4.36552e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE2_A: mean (A) | 5.4745 | 5.4745 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE2_A: sd (A) | 0.014321 | 0.014321 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47449 | 5.47449 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE2_A: sd_last_half (A) | 0.013193 | 0.013193 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE2_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 1.24432e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE_min_A: mean (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE_min_A: sd (A) | 0.016208 | 0.016208 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00932 | 4.00932 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.016848 | 0.016848 | 0 | 0.001 | second formulas | passed |
| test/X3 | d_Zn_E106OE_min_A: number of 10 ns blocks | 3 | 3 | 0 | 0 | construction | passed |
| test/X3 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 4.36552e-07 | 0 | 0.001 | second formulas | passed |
| test/X3 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| test/X3 | nearest water to the zinc, mean (A) | 3.561 | 3.561 | 0 | 0.001 | construction | passed |
| test/X3 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| test/X3 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| test/X3 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X3 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| test/X3 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| test/X3 | planned length (X0: 5 ns by the specification) (ns) | 30 | 30 | 0 | 1e-09 | construction | passed |
| test/X3 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| test/X3 | rule 1 applies to this run (production runs of start X) | yes | yes | – | equal | construction | passed |
| test/X3 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | GV1 ligand RMSD, GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| test/X3 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| test/X3 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| test/X3 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| test/X3 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| test/X3 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| test/X3 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| test/X3 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 2.6409 | 2.6409 | 0 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.117692 | 0.117691 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -5.57768 | -5.57768 | 0 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.831615 | 0.831616 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.894684 | 0.894684 | 0 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.731171 | 0.731171 | 0 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005401 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043932 | 1e-06 | 0.01 | second formulas | passed |
| test/X3 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 20 | 20 | 0 | 0 | second formulas | passed |
| test/X3 | gate G-V: solute frames counted for that fraction | 3000 | 3000 | 0 | 0 | construction | passed |
| test/X3 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 32 | 32 | 0 | 0 | second formulas | passed |
| test/X3 | gate G-V: full frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| test/X3 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 5.57768 | 5.57768 | 0 | 0.01 | second formulas | passed |
| test/X3 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | no | no | – | equal | construction | passed |
| test/X3 | campaign-table value: fraction of frames with ND2 nearer (none without the key) | none | none | – | equal | construction | passed |
| test/X3 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| test/X3 | N2: frames with positive, negative and zero signed volume | [1241, 1759, 0] | [1241, 1759, 0] | – | equal | second formulas | passed |
| test/X3 | N2: fraction of frames with positive signed volume | 0.413667 | 0.413667 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2: fraction of frames with negative signed volume | 0.586333 | 0.586333 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2: frames of the window (last 20 ns) | 2000 | 2000 | 0 | 0 | construction | passed |
| test/X3 | N2: positive, negative and zero in the last 20 ns | [519, 1481, 0] | [519, 1481, 0] | – | equal | second formulas | passed |
| test/X3 | N2: fraction positive in the last 20 ns | 0.2595 | 0.2595 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2, frames of state C: number, positive, negative, zero | [1051, 840, 211, 0] | [1051, 840, 211, 0] | – | equal | second formulas | passed |
| test/X3 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.799239 | 0.799239 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2, frames of state R: number, positive, negative, zero | [1929, 385, 1544, 0] | [1929, 385, 1544, 0] | – | equal | second formulas | passed |
| test/X3 | N2, frames of state R: fraction positive (none without a frame of that state) | 0.199585 | 0.199585 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2, frames of state I: number, positive, negative, zero | [20, 16, 4, 0] | [20, 16, 4, 0] | – | equal | second formulas | passed |
| test/X3 | N2, frames of state I: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X3 | N2: stretches of one sign, and those of at least 1 ns | [1192, 0] | [1192, 0] | – | equal | second formulas | passed |
| test/X3 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| test/X3 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| test/X3 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| test/X3 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| test/X3 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| test/X3 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| test/X3 | hydrogen bond around the first passage: frame of the first passage | 900 | 900 | – | equal | second formulas | passed |
| test/X3 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [900, 720] | [900, 720] | – | equal | second formulas | passed |
| test/X3 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [2100, 1664] | [2100, 1664] | – | equal | second formulas | passed |
| test/X3 | hydrogen bond N2-H...OG1, fraction before the first passage | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| test/X3 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.792381 | 0.792381 | 0 | 1e-06 | second formulas | passed |
| test/X3 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.8, 0.792381] | [0.8, 0.792381] | – | equal | construction | passed |
| test/X3 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| test/X3 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| test/X3 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 7e-06 | 0.001 | second formulas | passed |
| test/X3 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00892 | 4.00892 | 0 | 0.001 | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.4745 | 5.4745 | 0 | 0.001 | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 7e-06 | 0.001 | second formulas | passed |
| test/X3 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 3000] | [0, 3000] | – | equal | second formulas | passed |
| test/X3 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.008917, 4.003721, 0.0] | [4.008917, 4.003721, 0.0] | – | equal | second formulas | passed |
| test/X3 | water at the crystal ring positions: mean number within 3.5 A | 0.526667 | 0.526667 | 0 | 1e-06 | second formulas | passed |
| test/X3 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.2e-05 | 0.001 | second formulas | passed |
| test/X3 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 6.08635 | 6.08635 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | campaign-table value: mean number of water oxygens at the crystal ring positions | 0.526667 | 0.526667 | 0 | 1e-06 | second formulas | passed |
| test/X3 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| test/X3 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| test/X3 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X3 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| test/X3 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.2e-05 | 0.001 | second formulas | passed |
| test/X3 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 6.54828 | 6.54828 | 1e-06 | 0.001 | second formulas | passed |
| test/X3 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020121] | [0.0, 4.020121] | – | equal | second formulas | passed |
| test/X3 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| test/X3 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| test/X3 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| test/X3 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| test/X3 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| test/X3 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| test/X3 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| test/X3 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.413667, 0.2595, 0] | [0.413667, 0.2595, 0] | – | equal | second formulas | passed |
| test/X3 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| test/X3 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| test/X3 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| angle/X0 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| angle/X0 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| angle/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| angle/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| angle/X0 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| angle/X0 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| angle/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| angle/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| angle/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| angle/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| angle/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| angle/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| angle/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| angle/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.00072 | 0.00072 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.00068 | 0.00068 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000104 | 0.000104 | 0.01 | second formulas | passed |
| angle/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| angle/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| angle/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| angle/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| angle/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| angle/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| angle/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| angle/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| angle/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| angle/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| angle/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| angle/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| angle/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| angle/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| angle/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| angle/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| angle/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| angle/X0 | run complete | yes | yes | – | equal | construction | passed |
| angle/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| angle/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| angle/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| angle/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| angle/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| angle/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| angle/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| angle/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| angle/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| angle/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 8e-06 | 0.01 | second formulas | passed |
| angle/X0 | tau_fold_deg: sd (deg) | 0 | 7e-05 | 7e-05 | 0.01 | second formulas | passed |
| angle/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1.1e-05 | 0.01 | second formulas | passed |
| angle/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 7.4e-05 | 7.4e-05 | 0.01 | second formulas | passed |
| angle/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| angle/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 0 | 0.01 | second formulas | passed |
| angle/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.1e-05 | 5.1e-05 | 0.01 | second formulas | passed |
| angle/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| angle/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 4.5e-05 | 4.5e-05 | 0.01 | second formulas | passed |
| angle/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| angle/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| angle/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| angle/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | n2_volume: mean | 0.451372 | 0.451373 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| angle/X0 | n2_volume: mean_last_half | 0.451372 | 0.451373 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| angle/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| angle/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| angle/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 0 | 0.001 | construction | passed |
| angle/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| angle/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| angle/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| angle/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| angle/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| angle/X0 | gate G-V: the criteria that are not met | GV3 ang_H96NE2_Zn_H119ND1_deg | GV3 ang_H96NE2_Zn_H119ND1_deg | – | equal | second formulas | passed |
| angle/X0 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| angle/X0 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| angle/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| angle/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| angle/X0 | gate G-V: verdict (none for a run that is not X0) | FAILED | FAILED | – | equal | second formulas | passed |
| angle/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| angle/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823966 | 2e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840548 | 3e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.08366 | -7.08365 | 5e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731086 | 1e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -3.28926 | -3.28926 | 1e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 11.8451 | 11.8451 | 0 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005402 | 0 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043939 | 8e-06 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| angle/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| angle/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 11.8451 | 11.8451 | 0 | 0.01 | second formulas | passed |
| angle/X0 | gate G-V by design: verdict | FAILED | FAILED | – | equal | construction | passed |
| angle/X0 | gate G-V by design: the criteria that are not met | GV3 ang_H96NE2_Zn_H119ND1_deg | GV3 ang_H96NE2_Zn_H119ND1_deg | – | equal | construction | passed |
| angle/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| angle/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| angle/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| angle/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| angle/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| angle/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| angle/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| angle/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| angle/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| angle/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| angle/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| angle/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| angle/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| angle/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| angle/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| angle/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| angle/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| angle/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| angle/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| angle/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| angle/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| angle/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| angle/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| angle/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| angle/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| angle/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| angle/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| angle/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| angle/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| angle/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| angle/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| angle/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| angle/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| angle/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003723, 0.0] | [4.007685, 4.003723, 0.0] | – | equal | second formulas | passed |
| angle/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.2e-05 | 0.001 | second formulas | passed |
| angle/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.78 | 3.78 | 2e-06 | 0.001 | second formulas | passed |
| angle/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| angle/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| angle/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| angle/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.1e-05 | 0.001 | second formulas | passed |
| angle/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.83703 | 4.83703 | 1e-06 | 0.001 | second formulas | passed |
| angle/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020122] | [0.0, 4.020122] | – | equal | second formulas | passed |
| angle/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| angle/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| angle/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| angle/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| angle/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| angle/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| angle/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| angle/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| angle/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| angle/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| angle/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| csn/X0 | manifest lists the segments | no | no | – | equal | construction | passed |
| csn/X0 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| csn/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| csn/X0 | where the system folder was found | convention | convention | – | equal | construction | passed |
| csn/X0 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| csn/X0 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| csn/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| csn/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| csn/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| csn/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| csn/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| csn/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| csn/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| csn/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.9e-05 | 3.9e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 6e-06 | 6e-06 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| csn/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| csn/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000322 | 0.000322 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000104 | 0.000104 | 0.01 | second formulas | passed |
| csn/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| csn/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| csn/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| csn/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| csn/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| csn/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| csn/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| csn/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| csn/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| csn/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| csn/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| csn/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| csn/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| csn/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| csn/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| csn/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| csn/X0 | run complete | yes | yes | – | equal | construction | passed |
| csn/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| csn/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| csn/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| csn/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| csn/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| csn/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| csn/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| csn/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| csn/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| csn/X0 | lig_rmsd_A: mean (A) | 0.187409 | 0.187409 | 0 | 0.001 | second formulas | passed |
| csn/X0 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | lig_rmsd_A: mean_last_half (A) | 0.187409 | 0.187409 | 0 | 0.001 | second formulas | passed |
| csn/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | tau_fold_deg: sd (deg) | 0 | 7e-05 | 7e-05 | 0.01 | second formulas | passed |
| csn/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| csn/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| csn/X0 | ang_C4_S_N2_deg: mean (deg) | 112.295 | 112.295 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.4e-05 | 5.4e-05 | 0.01 | second formulas | passed |
| csn/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 112.295 | 112.295 | 1e-06 | 0.01 | second formulas | passed |
| csn/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.6e-05 | 5.6e-05 | 0.01 | second formulas | passed |
| csn/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| csn/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| csn/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| csn/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d1_N1_F131CZ_A: mean (A) | 3.56521 | 3.56521 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.56521 | 3.56521 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d2_N1_N67OD1_A: mean (A) | 7.96554 | 7.96554 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.96554 | 7.96554 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| csn/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| csn/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| csn/X0 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| csn/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| csn/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| csn/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| csn/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 1e-06 | 0.001 | construction | passed |
| csn/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| csn/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| csn/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| csn/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| csn/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| csn/X0 | gate G-V: the criteria that are not met | GV4 ang_C4_S_N2_deg | GV4 ang_C4_S_N2_deg | – | equal | second formulas | passed |
| csn/X0 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| csn/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| csn/X0 | gate G-V: GV4 as a whole | no | no | – | equal | second formulas | passed |
| csn/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| csn/X0 | gate G-V: verdict (none for a run that is not X0) | FAILED | FAILED | – | equal | second formulas | passed |
| csn/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| csn/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823966 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840542 | 3e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776336 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731089 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807206 | 1e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 3.4946 | 3.4946 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043933 | 2e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| csn/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| csn/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| csn/X0 | gate G-V by design: verdict | FAILED | FAILED | – | equal | construction | passed |
| csn/X0 | gate G-V by design: the criteria that are not met | GV4 ang_C4_S_N2_deg | GV4 ang_C4_S_N2_deg | – | equal | construction | passed |
| csn/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| csn/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| csn/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| csn/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| csn/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| csn/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| csn/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| csn/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| csn/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| csn/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| csn/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| csn/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| csn/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| csn/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| csn/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| csn/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| csn/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| csn/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| csn/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| csn/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| csn/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| csn/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| csn/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| csn/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| csn/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| csn/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| csn/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| csn/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| csn/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| csn/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| csn/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| csn/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| csn/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| csn/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003722, 0.0] | [4.007685, 4.003722, 0.0] | – | equal | second formulas | passed |
| csn/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 2.88812 | 2.88811 | 9e-06 | 0.001 | second formulas | passed |
| csn/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.72966 | 3.72966 | 0 | 0.001 | second formulas | passed |
| csn/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| csn/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| csn/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| csn/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.95591 | 3.9559 | 9e-06 | 0.001 | second formulas | passed |
| csn/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.82053 | 4.82053 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.955898] | [0.0, 3.955898] | – | equal | second formulas | passed |
| csn/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.88811 | 2.88811 | – | equal | second formulas | passed |
| csn/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| csn/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| csn/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| csn/X0 | first production frame: signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| csn/X0 | first production frame: ligand RMSD (A) | 0.187409 | 0.187409 | 0 | 0.001 | second formulas | passed |
| csn/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| csn/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| csn/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| csn/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| csn/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| edge/X0 | manifest lists the segments | no | no | – | equal | construction | passed |
| edge/X0 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| edge/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| edge/X0 | where the system folder was found | convention | convention | – | equal | construction | passed |
| edge/X0 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| edge/X0 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| edge/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| edge/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| edge/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| edge/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| edge/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| edge/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| edge/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 6e-06 | 6e-06 | 0.001 | second formulas | passed |
| edge/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.7e-05 | 3.7e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.4e-05 | 5.4e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| edge/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000488 | 0.000488 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000586 | 0.000586 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000222 | 0.000222 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| edge/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| edge/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| edge/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| edge/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| edge/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| edge/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| edge/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| edge/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| edge/X0 | every frame: n2_volume, largest difference | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 4.6e-05 | 4.6e-05 | 0.001 | second formulas | passed |
| edge/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| edge/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| edge/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| edge/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| edge/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.9e-05 | 3.9e-05 | 0.001 | second formulas | passed |
| edge/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| edge/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| edge/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| edge/X0 | run complete | yes | yes | – | equal | construction | passed |
| edge/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| edge/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| edge/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| edge/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| edge/X0 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| edge/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| edge/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| edge/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| edge/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| edge/X0 | lig_rmsd_A: mean (A) | 0.272894 | 0.272894 | 0 | 0.001 | second formulas | passed |
| edge/X0 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| edge/X0 | lig_rmsd_A: mean_last_half (A) | 0.272894 | 0.272894 | 0 | 0.001 | second formulas | passed |
| edge/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| edge/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_N2_A: mean (A) | 1.95008 | 1.95008 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_N2_A: sd (A) | 0.105993 | 0.105993 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_N2_A: mean_last_half (A) | 1.95053 | 1.95053 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_N2_A: sd_last_half (A) | 0.105424 | 0.105424 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | tau_fold_deg: mean (deg) | 151.056 | 151.056 | 1e-06 | 0.01 | second formulas | passed |
| edge/X0 | tau_fold_deg: sd (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| edge/X0 | tau_fold_deg: mean_last_half (deg) | 151.056 | 151.056 | 0 | 0.01 | second formulas | passed |
| edge/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 7.7e-05 | 7.7e-05 | 0.01 | second formulas | passed |
| edge/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| edge/X0 | ang_C4_S_N2_deg: mean (deg) | 111.295 | 111.295 | 5e-06 | 0.01 | second formulas | passed |
| edge/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 6.6e-05 | 6.6e-05 | 0.01 | second formulas | passed |
| edge/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 111.295 | 111.295 | 6e-06 | 0.01 | second formulas | passed |
| edge/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 6.3e-05 | 6.3e-05 | 0.01 | second formulas | passed |
| edge/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| edge/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| edge/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| edge/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| edge/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| edge/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d1_N1_F131CZ_A: mean (A) | 3.66633 | 3.66633 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| edge/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.66633 | 3.66633 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| edge/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d2_N1_N67OD1_A: mean (A) | 7.86093 | 7.86093 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d2_N1_N67OD1_A: sd (A) | 0.235252 | 0.235252 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.86093 | 7.86093 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0.235252 | 0.235252 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | n2_volume: mean | 0.454384 | 0.454384 | 0 | 0.001 | second formulas | passed |
| edge/X0 | n2_volume: sd | 0.547205 | 0.547205 | 0 | 0.001 | second formulas | passed |
| edge/X0 | n2_volume: mean_last_half | 0.45401 | 0.45401 | 0 | 0.001 | second formulas | passed |
| edge/X0 | n2_volume: sd_last_half | 0.547294 | 0.547294 | 0 | 0.001 | second formulas | passed |
| edge/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE1_A: mean (A) | 4.01917 | 4.01917 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE1_A: sd (A) | 0.050562 | 0.050562 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01871 | 4.01871 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.049649 | 0.049649 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE2_A: mean (A) | 5.48976 | 5.48976 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE2_A: sd (A) | 0.068532 | 0.068532 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.48913 | 5.48913 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.06726 | 0.06726 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE_min_A: mean (A) | 4.01917 | 4.01917 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE_min_A: sd (A) | 0.050562 | 0.050562 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01871 | 4.01871 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.049649 | 0.049649 | 0 | 0.001 | second formulas | passed |
| edge/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| edge/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| edge/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| edge/X0 | nearest water to the zinc, mean (A) | 3.6642 | 3.6642 | 1e-06 | 0.001 | construction | passed |
| edge/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| edge/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| edge/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| edge/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| edge/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| edge/X0 | gate G-V: the criteria that are not met | none | none | – | equal | second formulas | passed |
| edge/X0 | gate G-V: all criteria met | yes | yes | – | equal | second formulas | passed |
| edge/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| edge/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| edge/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| edge/X0 | gate G-V: verdict (none for a run that is not X0) | PASSED | PASSED | – | equal | second formulas | passed |
| edge/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| edge/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.400692 | -0.400691 | 1e-06 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -1.41523 | -1.41523 | 1e-06 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -5.09239 | -5.09239 | 0 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.058368 | -0.05837 | 2e-06 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -2.3708 | -2.3708 | 0 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 9.15181 | 9.15181 | 0 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 2.4946 | 2.4946 | 5e-06 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 18.9561 | 18.9561 | 1e-06 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 25 | 25 | 0 | 0 | second formulas | passed |
| edge/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| edge/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| edge/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 9.15181 | 9.15181 | 0 | 0.01 | second formulas | passed |
| edge/X0 | gate G-V by design: verdict | PASSED | PASSED | – | equal | construction | passed |
| edge/X0 | gate G-V by design: the criteria that are not met | none | none | – | equal | construction | passed |
| edge/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| edge/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 25 | 25 | 0 | 0 | construction | passed |
| edge/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| edge/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| edge/X0 | frames with ND2 nearer to the aniline N than OD1 | 480 | 480 | 0 | 0 | second formulas | passed |
| edge/X0 | state C: frames, and frames with ND2 nearer | [500, 480] | [500, 480] | – | equal | second formulas | passed |
| edge/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| edge/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| edge/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [480, 500, 0, 0] | [480, 500, 0, 0] | – | equal | construction | passed |
| edge/X0 | campaign-table value: fraction of frames with ND2 nearer | 0.96 | 0.96 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| edge/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| edge/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| edge/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| edge/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| edge/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| edge/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| edge/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| edge/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| edge/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| edge/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| edge/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| edge/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| edge/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| edge/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| edge/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| edge/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| edge/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| edge/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| edge/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| edge/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| edge/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| edge/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01917 | 4.01917 | 0 | 0.001 | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.48976 | 5.48976 | 0 | 0.001 | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| edge/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| edge/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.019171, 4.003722, 0.0] | [4.019171, 4.003722, 0.0] | – | equal | second formulas | passed |
| edge/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 2.91197 | 2.91196 | 9e-06 | 0.001 | second formulas | passed |
| edge/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.74072 | 3.74072 | 0 | 0.001 | second formulas | passed |
| edge/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| edge/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| edge/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| edge/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.9782 | 3.97819 | 1.1e-05 | 0.001 | second formulas | passed |
| edge/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.82709 | 4.82709 | 0 | 0.001 | second formulas | passed |
| edge/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.978193] | [0.0, 3.978193] | – | equal | second formulas | passed |
| edge/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.91196 | 2.91196 | – | equal | second formulas | passed |
| edge/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| edge/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| edge/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| edge/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| edge/X0 | first production frame: ligand RMSD (A) | 0.272894 | 0.272893 | 1e-06 | 0.001 | second formulas | passed |
| edge/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| edge/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| edge/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| edge/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| edge/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | solute frames seen | 400 | 400 | 0 | 0 | construction | passed |
| gate/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| gate/R1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| gate/R1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| gate/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| gate/R1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| gate/R1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| gate/R1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| gate/R1 | full frames seen | 40 | 40 | 0 | 0 | construction | passed |
| gate/R1 | time seen (ns) | 4 | 4 | 0 | 1e-09 | construction | passed |
| gate/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| gate/R1 | frames in which the ligand was in another periodic image | 57 | 57 | 0 | 0 | construction | passed |
| gate/R1 | frames in which the zinc was in another periodic image | 36 | 36 | 0 | 0 | construction | passed |
| gate/R1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | lines of timeseries_solute.csv | 400 | 400 | 0 | 0 | construction | passed |
| gate/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| gate/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| gate/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000269 | 0.000269 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000322 | 0.000322 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| gate/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| gate/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| gate/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| gate/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| gate/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| gate/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| gate/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| gate/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| gate/R1 | every frame: n2_volume, largest difference | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| gate/R1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | lines of timeseries_full.csv | 40 | 40 | 0 | 0 | construction | passed |
| gate/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| gate/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| gate/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| gate/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| gate/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| gate/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| gate/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 4.3e-05 | 4.3e-05 | 0.001 | second formulas | passed |
| gate/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| gate/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | verdict | MIXED | MIXED | – | equal | construction | passed |
| gate/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| gate/R1 | run complete | yes | yes | – | equal | construction | passed |
| gate/R1 | first passage into the opposite state (ns) | 2.51 | 2.51 | 0 | 1e-09 | construction | passed |
| gate/R1 | transitions | 1 | 1 | 0 | 0 | construction | passed |
| gate/R1 | long visits (at least 1 ns) | 2 | 2 | 0 | 0 | construction | passed |
| gate/R1 | frames of the verdict window | 400 | 400 | 0 | 0 | construction | passed |
| gate/R1 | count_C_window | 150 | 150 | 0 | 0 | construction | passed |
| gate/R1 | count_R_window | 250 | 250 | 0 | 0 | construction | passed |
| gate/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | count_C | 150 | 150 | 0 | 0 | construction | passed |
| gate/R1 | count_R | 250 | 250 | 0 | 0 | construction | passed |
| gate/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| gate/R1 | count_C (design against the planted frames) | 150 | 150 | 0 | 0 | construction | passed |
| gate/R1 | count_R (design against the planted frames) | 250 | 250 | 0 | 0 | construction | passed |
| gate/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | frames with the hydrogen bond N2-H...OG1 | 320 | 320 | 0 | 0 | second formulas | passed |
| gate/R1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| gate/R1 | frames with labelling identity | 287 | 287 | 0 | 0 | second formulas | passed |
| gate/R1 | frames with labelling ring_flip | 37 | 37 | 0 | 0 | second formulas | passed |
| gate/R1 | frames with labelling oxygen_swap | 38 | 38 | 0 | 0 | second formulas | passed |
| gate/R1 | frames with labelling ring_flip+oxygen_swap | 38 | 38 | 0 | 0 | second formulas | passed |
| gate/R1 | lig_rmsd_A: mean (A) | 2.23326 | 2.23327 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | lig_rmsd_A: sd (A) | 1.72988 | 1.72988 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | lig_rmsd_A: mean_last_half (A) | 0.893306 | 0.893307 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | lig_rmsd_A: sd_last_half (A) | 1.54725 | 1.54725 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_N2_A: mean (A) | 2.15126 | 2.15126 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_N2_A: sd (A) | 0.172389 | 0.172389 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_N2_A: mean_last_half (A) | 2.03737 | 2.03737 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_N2_A: sd_last_half (A) | 0.159356 | 0.159356 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | tau_fold_deg: mean (deg) | 138.736 | 138.736 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | tau_fold_deg: sd (deg) | 5.1745 | 5.1745 | 4e-06 | 0.01 | second formulas | passed |
| gate/R1 | tau_fold_deg: mean_last_half (deg) | 134.728 | 134.728 | 9e-06 | 0.01 | second formulas | passed |
| gate/R1 | tau_fold_deg: sd_last_half (deg) | 4.62822 | 4.62822 | 3e-06 | 0.01 | second formulas | passed |
| gate/R1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| gate/R1 | ang_C4_S_N2_deg: mean (deg) | 109.193 | 109.193 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | ang_C4_S_N2_deg: sd (deg) | 0.308544 | 0.308545 | 1e-06 | 0.01 | second formulas | passed |
| gate/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.954 | 108.954 | 1e-06 | 0.01 | second formulas | passed |
| gate/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0.27597 | 0.275967 | 3e-06 | 0.01 | second formulas | passed |
| gate/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| gate/R1 | d_S_N2_A: mean (A) | 1.64691 | 1.64691 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_S_N2_A: sd (A) | 0.039174 | 0.039174 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_S_N2_A: mean_last_half (A) | 1.61656 | 1.61656 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_S_N2_A: sd_last_half (A) | 0.035038 | 0.035039 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| gate/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| gate/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d1_N1_F131CZ_A: mean (A) | 7.45182 | 7.45182 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d1_N1_F131CZ_A: sd (A) | 2.75386 | 2.75386 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 5.31869 | 5.31869 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 2.46313 | 2.46313 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d2_N1_N67OD1_A: mean (A) | 4.5291 | 4.5291 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | d2_N1_N67OD1_A: sd (A) | 2.29742 | 2.29742 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 6.28306 | 6.28306 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 2.07559 | 2.07559 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | n2_volume: mean | -0.129702 | -0.129702 | 0 | 0.001 | second formulas | passed |
| gate/R1 | n2_volume: sd | 0.669859 | 0.669859 | 0 | 0.001 | second formulas | passed |
| gate/R1 | n2_volume: mean_last_half | 0.218943 | 0.218943 | 0 | 0.001 | second formulas | passed |
| gate/R1 | n2_volume: sd_last_half | 0.662486 | 0.662486 | 0 | 0.001 | second formulas | passed |
| gate/R1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE1_A: mean (A) | 4.00974 | 4.00974 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE1_A: sd (A) | 0.018716 | 0.018716 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00851 | 4.00851 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.015025 | 0.015025 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE2_A: mean (A) | 5.47553 | 5.47553 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE2_A: sd (A) | 0.016716 | 0.016716 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47465 | 5.47465 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.014059 | 0.014059 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE_min_A: mean (A) | 4.00974 | 4.00974 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE_min_A: sd (A) | 0.018716 | 0.018716 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00851 | 4.00851 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.015025 | 0.015025 | 0 | 0.001 | second formulas | passed |
| gate/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| gate/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| gate/R1 | nearest water to the zinc, mean (A) | 3.65775 | 3.65775 | 0 | 0.001 | construction | passed |
| gate/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| gate/R1 | frame 2, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| gate/R1 | frame 2, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| gate/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (A) | 3.31 | 3.3102 | 0.0002 | 0.0015 | starts.csv | passed |
| gate/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| gate/R1 | planned length (X0: 5 ns by the specification) (ns) | 4 | 4 | 0 | 1e-09 | construction | passed |
| gate/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| gate/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| gate/R1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | – | equal | second formulas | passed |
| gate/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| gate/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| gate/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| gate/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| gate/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| gate/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| gate/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 3.83005 | 3.83006 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -1.91433 | -1.91433 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -4.7513 | -4.75131 | 5e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.810514 | 0.810516 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.950134 | 0.950136 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.707467 | 0.707466 | 1e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.392927 | 0.392929 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 6.63633 | 6.63633 | 2e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | gate G-V: solute frames counted for that fraction | 400 | 400 | 0 | 0 | construction | passed |
| gate/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| gate/R1 | gate G-V: full frames counted for that fraction | 40 | 40 | 0 | 0 | construction | passed |
| gate/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 4.7513 | 4.75131 | 5e-06 | 0.01 | second formulas | passed |
| gate/R1 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| gate/R1 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| gate/R1 | by design: full frames with any oxygen within 2.6 A | 1 | 1 | 0 | 0 | construction | passed |
| gate/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| gate/R1 | frames with ND2 nearer to the aniline N than OD1 | 176 | 176 | 0 | 0 | second formulas | passed |
| gate/R1 | state C: frames, and frames with ND2 nearer | [150, 113] | [150, 113] | – | equal | second formulas | passed |
| gate/R1 | state R: frames, and frames with ND2 nearer | [250, 63] | [250, 63] | – | equal | second formulas | passed |
| gate/R1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| gate/R1 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [113, 150, 63, 250] | [113, 150, 63, 250] | – | equal | construction | passed |
| gate/R1 | campaign-table value: fraction of frames with ND2 nearer | 0.44 | 0.44 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| gate/R1 | N2: frames with positive, negative and zero signed volume | [170, 230, 0] | [170, 230, 0] | – | equal | second formulas | passed |
| gate/R1 | N2: fraction of frames with positive signed volume | 0.425 | 0.425 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2: fraction of frames with negative signed volume | 0.575 | 0.575 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2: frames of the window (last 20 ns) | 400 | 400 | 0 | 0 | construction | passed |
| gate/R1 | N2: positive, negative and zero in the last 20 ns | [170, 230, 0] | [170, 230, 0] | – | equal | second formulas | passed |
| gate/R1 | N2: fraction positive in the last 20 ns | 0.425 | 0.425 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2, frames of state C: number, positive, negative, zero | [150, 120, 30, 0] | [150, 120, 30, 0] | – | equal | second formulas | passed |
| gate/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2, frames of state R: number, positive, negative, zero | [250, 50, 200, 0] | [250, 50, 200, 0] | – | equal | second formulas | passed |
| gate/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| gate/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| gate/R1 | N2: stretches of one sign, and those of at least 1 ns | [159, 0] | [159, 0] | – | equal | second formulas | passed |
| gate/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| gate/R1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| gate/R1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| gate/R1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| gate/R1 | N2: sign of the first and of the last frame | [-1, -1] | [-1, -1] | – | equal | second formulas | passed |
| gate/R1 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| gate/R1 | hydrogen bond around the first passage: frame of the first passage | 250 | 250 | – | equal | second formulas | passed |
| gate/R1 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [250, 200] | [250, 200] | – | equal | second formulas | passed |
| gate/R1 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [150, 120] | [150, 120] | – | equal | second formulas | passed |
| gate/R1 | hydrogen bond N2-H...OG1, fraction before the first passage | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.8, 0.8] | [0.8, 0.8] | – | equal | construction | passed |
| gate/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| gate/R1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| gate/R1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| gate/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00974 | 4.00974 | 0 | 0.001 | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 400] | [0, 400] | – | equal | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47553 | 5.47553 | 0 | 0.001 | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| gate/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 400] | [0, 400] | – | equal | second formulas | passed |
| gate/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.009736, 4.003722, 0.0] | [4.009736, 4.003722, 0.0] | – | equal | second formulas | passed |
| gate/R1 | water at the crystal ring positions: mean number within 3.5 A | 0.6 | 0.6 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99998 | 1.6e-05 | 0.001 | second formulas | passed |
| gate/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 6.01826 | 6.01826 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 0.6 | 0.6 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| gate/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| gate/R1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| gate/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 9e-06 | 0.001 | second formulas | passed |
| gate/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 6.50447 | 6.50447 | 0 | 0.001 | second formulas | passed |
| gate/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020124] | [0.0, 4.020124] | – | equal | second formulas | passed |
| gate/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99998 | 2.99998 | – | equal | second formulas | passed |
| gate/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| gate/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| gate/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| gate/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| gate/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| gate/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| gate/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.425, 0.425, 0] | [0.425, 0.425, 0] | – | equal | second formulas | passed |
| gate/R1 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| gate/R1 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| gate/R1 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | solute frames seen | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| long/X0 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| long/X0 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| long/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| long/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| long/X0 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| long/X0 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| long/X0 | full frames seen | 60 | 60 | 0 | 0 | construction | passed |
| long/X0 | time seen (ns) | 6 | 6 | 0 | 1e-09 | construction | passed |
| long/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| long/X0 | frames in which the ligand was in another periodic image | 86 | 86 | 0 | 0 | construction | passed |
| long/X0 | frames in which the zinc was in another periodic image | 55 | 55 | 0 | 0 | construction | passed |
| long/X0 | warnings | 2 | 2 | 0 | 0 | construction | passed |
| long/X0 | lines of timeseries_solute.csv | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| long/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| long/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| long/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| long/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| long/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| long/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| long/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| long/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| long/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| long/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| long/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | lines of timeseries_full.csv | 60 | 60 | 0 | 0 | construction | passed |
| long/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| long/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| long/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| long/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| long/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| long/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| long/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| long/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| long/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| long/X0 | run complete | yes | yes | – | equal | construction | passed |
| long/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| long/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| long/X0 | frames of the verdict window | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | count_C_window | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | count_C | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| long/X0 | count_C (design against the planted frames) | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | frames with the hydrogen bond N2-H...OG1 | 480 | 480 | 0 | 0 | second formulas | passed |
| long/X0 | full frames with a water within 2.6 A of the zinc | 3 | 3 | 0 | 0 | construction | passed |
| long/X0 | frames with labelling identity | 150 | 150 | 0 | 0 | second formulas | passed |
| long/X0 | frames with labelling ring_flip | 150 | 150 | 0 | 0 | second formulas | passed |
| long/X0 | frames with labelling oxygen_swap | 150 | 150 | 0 | 0 | second formulas | passed |
| long/X0 | frames with labelling ring_flip+oxygen_swap | 150 | 150 | 0 | 0 | second formulas | passed |
| long/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 0 | 0.01 | second formulas | passed |
| long/X0 | tau_fold_deg: sd (deg) | 0 | 7.5e-05 | 7.5e-05 | 0.01 | second formulas | passed |
| long/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 5e-06 | 0.01 | second formulas | passed |
| long/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| long/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| long/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| long/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| long/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 4e-06 | 0.01 | second formulas | passed |
| long/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.9e-05 | 5.9e-05 | 0.01 | second formulas | passed |
| long/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| long/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| long/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| long/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| long/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| long/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| long/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| long/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| long/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| long/X0 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| long/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| long/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| long/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| long/X0 | nearest water to the zinc, mean (A) | 3.6255 | 3.6255 | 0 | 0.001 | construction | passed |
| long/X0 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| long/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| long/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| long/X0 | frame 8, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| long/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| long/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| long/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| long/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| long/X0 | gate G-V: the criteria that are not met | none | none | – | equal | second formulas | passed |
| long/X0 | gate G-V: all criteria met | yes | yes | – | equal | second formulas | passed |
| long/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| long/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| long/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| long/X0 | gate G-V: verdict (none for a run that is not X0) | PASSED | PASSED | – | equal | second formulas | passed |
| long/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| long/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823966 | 2e-06 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840543 | 2e-06 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776334 | 0 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731087 | 0 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930668 | 1e-06 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807205 | 0 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.0054 | 2e-06 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043931 | 0 | 0.01 | second formulas | passed |
| long/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | gate G-V: solute frames counted for that fraction | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 3 | 3 | 0 | 0 | second formulas | passed |
| long/X0 | gate G-V: full frames counted for that fraction | 60 | 60 | 0 | 0 | construction | passed |
| long/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930668 | 1e-06 | 0.01 | second formulas | passed |
| long/X0 | gate G-V by design: verdict | PASSED | PASSED | – | equal | construction | passed |
| long/X0 | gate G-V by design: the criteria that are not met | none | none | – | equal | construction | passed |
| long/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| long/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| long/X0 | by design: full frames with any oxygen within 2.6 A | 3 | 3 | 0 | 0 | construction | passed |
| long/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| long/X0 | frames with ND2 nearer to the aniline N than OD1 | 600 | 600 | 0 | 0 | second formulas | passed |
| long/X0 | state C: frames, and frames with ND2 nearer | [600, 600] | [600, 600] | – | equal | second formulas | passed |
| long/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| long/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| long/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [600, 600, 0, 0] | [600, 600, 0, 0] | – | equal | construction | passed |
| long/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| long/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| long/X0 | N2: frames with positive, negative and zero signed volume | [480, 120, 0] | [480, 120, 0] | – | equal | second formulas | passed |
| long/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| long/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| long/X0 | N2: frames of the window (last 20 ns) | 600 | 600 | 0 | 0 | construction | passed |
| long/X0 | N2: positive, negative and zero in the last 20 ns | [480, 120, 0] | [480, 120, 0] | – | equal | second formulas | passed |
| long/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| long/X0 | N2, frames of state C: number, positive, negative, zero | [600, 480, 120, 0] | [600, 480, 120, 0] | – | equal | second formulas | passed |
| long/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| long/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| long/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| long/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| long/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| long/X0 | N2: stretches of one sign, and those of at least 1 ns | [240, 0] | [240, 0] | – | equal | second formulas | passed |
| long/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| long/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| long/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| long/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| long/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| long/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| long/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| long/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| long/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| long/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| long/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| long/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| long/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 600] | [0, 600] | – | equal | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| long/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 600] | [0, 600] | – | equal | second formulas | passed |
| long/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003721, 0.0] | [4.007685, 4.003721, 0.0] | – | equal | second formulas | passed |
| long/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| long/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.2e-05 | 0.001 | second formulas | passed |
| long/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.75 | 3.75 | 2e-06 | 0.001 | second formulas | passed |
| long/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| long/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| long/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| long/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| long/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| long/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1e-05 | 0.001 | second formulas | passed |
| long/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.80561 | 4.80561 | 1e-06 | 0.001 | second formulas | passed |
| long/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020123] | [0.0, 4.020123] | – | equal | second formulas | passed |
| long/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| long/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| long/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| long/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| long/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| long/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| long/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| long/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| long/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| long/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| long/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | solute frames seen | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| short/X0 | manifest lists the segments | no | no | – | equal | construction | passed |
| short/X0 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| short/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| short/X0 | where the system folder was found | convention | convention | – | equal | construction | passed |
| short/X0 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| short/X0 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| short/X0 | full frames seen | 30 | 30 | 0 | 0 | construction | passed |
| short/X0 | time seen (ns) | 3 | 3 | 0 | 1e-09 | construction | passed |
| short/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| short/X0 | frames in which the ligand was in another periodic image | 43 | 43 | 0 | 0 | construction | passed |
| short/X0 | frames in which the zinc was in another periodic image | 27 | 27 | 0 | 0 | construction | passed |
| short/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | lines of timeseries_solute.csv | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| short/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| short/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| short/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000322 | 0.000322 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| short/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| short/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| short/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| short/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| short/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.001069 | 0.001069 | 0.01 | second formulas | passed |
| short/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| short/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| short/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | lines of timeseries_full.csv | 30 | 30 | 0 | 0 | construction | passed |
| short/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| short/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| short/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| short/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| short/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| short/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| short/X0 | verdict provisional (run not finished) | yes | yes | – | equal | construction | passed |
| short/X0 | run complete | no | no | – | equal | construction | passed |
| short/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| short/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| short/X0 | frames of the verdict window | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | count_C_window | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | count_C | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| short/X0 | count_C (design against the planted frames) | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | frames with the hydrogen bond N2-H...OG1 | 240 | 240 | 0 | 0 | second formulas | passed |
| short/X0 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| short/X0 | frames with labelling identity | 75 | 75 | 0 | 0 | second formulas | passed |
| short/X0 | frames with labelling ring_flip | 75 | 75 | 0 | 0 | second formulas | passed |
| short/X0 | frames with labelling oxygen_swap | 75 | 75 | 0 | 0 | second formulas | passed |
| short/X0 | frames with labelling ring_flip+oxygen_swap | 75 | 75 | 0 | 0 | second formulas | passed |
| short/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 3e-06 | 0.01 | second formulas | passed |
| short/X0 | tau_fold_deg: sd (deg) | 0 | 8.1e-05 | 8.1e-05 | 0.01 | second formulas | passed |
| short/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| short/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 6.6e-05 | 6.6e-05 | 0.01 | second formulas | passed |
| short/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| short/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 0 | 0.01 | second formulas | passed |
| short/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 6e-05 | 6e-05 | 0.01 | second formulas | passed |
| short/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 3e-06 | 0.01 | second formulas | passed |
| short/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.3e-05 | 5.3e-05 | 0.01 | second formulas | passed |
| short/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| short/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| short/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| short/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| short/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| short/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| short/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| short/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| short/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| short/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| short/X0 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| short/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| short/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011875 | 1e-06 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| short/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| short/X0 | nearest water to the zinc, mean (A) | 3.647 | 3.647 | 0 | 0.001 | construction | passed |
| short/X0 | water molecules within 3.5 A of the aniline N, mean | 1.43333 | 1.43333 | 0 | 1e-06 | construction | passed |
| short/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| short/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| short/X0 | frame 8, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| short/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| short/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| short/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| short/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| short/X0 | gate G-V: the criteria that are not met | none | none | – | equal | second formulas | passed |
| short/X0 | gate G-V: all criteria met | yes | yes | – | equal | second formulas | passed |
| short/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| short/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| short/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| short/X0 | gate G-V: verdict (none for a run that is not X0) | PASSED | PASSED | – | equal | second formulas | passed |
| short/X0 | gate G-V: provisional (run not finished) | yes | yes | – | equal | construction | passed |
| short/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823964 | 0 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840544 | 1e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776335 | 1e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731088 | 1e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930669 | 2e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807203 | 2e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005402 | 0 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043934 | 3e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | gate G-V: solute frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| short/X0 | gate G-V: full frames counted for that fraction | 30 | 30 | 0 | 0 | construction | passed |
| short/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930669 | 2e-06 | 0.01 | second formulas | passed |
| short/X0 | gate G-V by design: verdict | PASSED | PASSED | – | equal | construction | passed |
| short/X0 | gate G-V by design: the criteria that are not met | none | none | – | equal | construction | passed |
| short/X0 | gate G-V by design: provisional | yes | yes | – | equal | construction | passed |
| short/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| short/X0 | by design: full frames with any oxygen within 2.6 A | 1 | 1 | 0 | 0 | construction | passed |
| short/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| short/X0 | frames with ND2 nearer to the aniline N than OD1 | 300 | 300 | 0 | 0 | second formulas | passed |
| short/X0 | state C: frames, and frames with ND2 nearer | [300, 300] | [300, 300] | – | equal | second formulas | passed |
| short/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| short/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| short/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [300, 300, 0, 0] | [300, 300, 0, 0] | – | equal | construction | passed |
| short/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| short/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| short/X0 | N2: frames with positive, negative and zero signed volume | [240, 60, 0] | [240, 60, 0] | – | equal | second formulas | passed |
| short/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| short/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| short/X0 | N2: frames of the window (last 20 ns) | 300 | 300 | 0 | 0 | construction | passed |
| short/X0 | N2: positive, negative and zero in the last 20 ns | [240, 60, 0] | [240, 60, 0] | – | equal | second formulas | passed |
| short/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| short/X0 | N2, frames of state C: number, positive, negative, zero | [300, 240, 60, 0] | [300, 240, 60, 0] | – | equal | second formulas | passed |
| short/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| short/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| short/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| short/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| short/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| short/X0 | N2: stretches of one sign, and those of at least 1 ns | [120, 0] | [120, 0] | – | equal | second formulas | passed |
| short/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| short/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| short/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| short/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| short/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| short/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| short/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| short/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| short/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| short/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| short/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| short/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 3e-06 | 0.001 | second formulas | passed |
| short/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 300] | [0, 300] | – | equal | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 3e-06 | 0.001 | second formulas | passed |
| short/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 300] | [0, 300] | – | equal | second formulas | passed |
| short/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003723, 0.0] | [4.007685, 4.003723, 0.0] | – | equal | second formulas | passed |
| short/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.43333 | 1.43333 | 0 | 1e-06 | second formulas | passed |
| short/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 9e-06 | 0.001 | second formulas | passed |
| short/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.8 | 3.8 | 3e-06 | 0.001 | second formulas | passed |
| short/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.43333 | 1.43333 | 0 | 1e-06 | second formulas | passed |
| short/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| short/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| short/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| short/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| short/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02013 | 7e-06 | 0.001 | second formulas | passed |
| short/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.85798 | 4.85798 | 0 | 0.001 | second formulas | passed |
| short/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020126] | [0.0, 4.020126] | – | equal | second formulas | passed |
| short/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| short/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| short/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| short/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| short/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| short/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| short/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| short/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| short/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| short/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| short/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| sulfonyl/X0 | manifest lists the segments | no | no | – | equal | construction | passed |
| sulfonyl/X0 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| sulfonyl/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| sulfonyl/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| sulfonyl/X0 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| sulfonyl/X0 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| sulfonyl/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| sulfonyl/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| sulfonyl/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| sulfonyl/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| sulfonyl/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| sulfonyl/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| sulfonyl/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001128 | 0.001128 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000512 | 0.000512 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000354 | 0.000354 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| sulfonyl/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| sulfonyl/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| sulfonyl/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| sulfonyl/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| sulfonyl/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| sulfonyl/X0 | run complete | yes | yes | – | equal | construction | passed |
| sulfonyl/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| sulfonyl/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| sulfonyl/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| sulfonyl/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| sulfonyl/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_N2_A: mean (A) | 1.94963 | 1.94963 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_N2_A: sd (A) | 0.106558 | 0.106558 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_N2_A: mean_last_half (A) | 1.95053 | 1.95053 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_N2_A: sd_last_half (A) | 0.105424 | 0.105424 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | tau_fold_deg: sd (deg) | 0 | 8e-05 | 8e-05 | 0.01 | second formulas | passed |
| sulfonyl/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 6e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 8.4e-05 | 8.4e-05 | 0.01 | second formulas | passed |
| sulfonyl/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 4e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.4e-05 | 5.4e-05 | 0.01 | second formulas | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 4e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.3e-05 | 5.3e-05 | 0.01 | second formulas | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| sulfonyl/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | n2_volume: mean | 0.454758 | 0.454758 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | n2_volume: sd | 0.547115 | 0.547115 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | n2_volume: mean_last_half | 0.454498 | 0.454498 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | n2_volume: sd_last_half | 0.546971 | 0.546971 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: mean (A) | 4.01963 | 4.01963 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: sd (A) | 0.051455 | 0.051455 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01871 | 4.01871 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.049649 | 0.049649 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: mean (A) | 5.49039 | 5.49039 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: sd (A) | 0.069775 | 0.069775 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.48913 | 5.48913 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.06726 | 0.06726 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: mean (A) | 4.01963 | 4.01963 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: sd (A) | 0.051455 | 0.051455 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01871 | 4.01871 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.049649 | 0.049649 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| sulfonyl/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| sulfonyl/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 0 | 0.001 | construction | passed |
| sulfonyl/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| sulfonyl/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | frame 8, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| sulfonyl/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| sulfonyl/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| sulfonyl/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| sulfonyl/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| sulfonyl/X0 | gate G-V: the criteria that are not met | GV3 fifth ligand, solute trajectory | GV3 fifth ligand, solute trajectory | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: verdict (none for a run that is not X0) | FAILED | FAILED | – | equal | second formulas | passed |
| sulfonyl/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| sulfonyl/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.383761 | -0.383761 | 0 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -1.43822 | -1.43822 | 3e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.106532 | -0.10653 | 2e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.089947 | -0.089946 | 1e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.872159 | 0.872159 | 0 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.094466 | 0.094465 | 1e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005406 | 4e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.04393 | 1e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 26 | 26 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| sulfonyl/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 1.43822 | 1.43822 | 3e-06 | 0.01 | second formulas | passed |
| sulfonyl/X0 | gate G-V by design: verdict | FAILED | FAILED | – | equal | construction | passed |
| sulfonyl/X0 | gate G-V by design: the criteria that are not met | GV3 fifth ligand, solute trajectory | GV3 fifth ligand, solute trajectory | – | equal | construction | passed |
| sulfonyl/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| sulfonyl/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 26 | 26 | 0 | 0 | construction | passed |
| sulfonyl/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| sulfonyl/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| sulfonyl/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| sulfonyl/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| sulfonyl/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| sulfonyl/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| sulfonyl/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| sulfonyl/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| sulfonyl/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| sulfonyl/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| sulfonyl/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| sulfonyl/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| sulfonyl/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| sulfonyl/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01963 | 4.01963 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.49039 | 5.49039 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| sulfonyl/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.01963, 4.003721, 0.0] | [4.01963, 4.003721, 0.0] | – | equal | second formulas | passed |
| sulfonyl/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.78 | 3.78 | 2e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| sulfonyl/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| sulfonyl/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| sulfonyl/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1e-05 | 0.001 | second formulas | passed |
| sulfonyl/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.83703 | 4.83703 | 1e-06 | 0.001 | second formulas | passed |
| sulfonyl/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020123] | [0.0, 4.020123] | – | equal | second formulas | passed |
| sulfonyl/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| sulfonyl/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| sulfonyl/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| sulfonyl/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| sulfonyl/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| sulfonyl/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| sulfonyl/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| sulfonyl/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| sulfonyl/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| sulfonyl/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| tau/X0 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| tau/X0 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| tau/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| tau/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| tau/X0 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| tau/X0 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| tau/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| tau/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| tau/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| tau/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| tau/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| tau/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| tau/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.1e-05 | 2.1e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 3.5e-05 | 3.5e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| tau/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| tau/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000355 | 0.000355 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000104 | 0.000104 | 0.01 | second formulas | passed |
| tau/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| tau/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| tau/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| tau/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| tau/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| tau/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| tau/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| tau/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| tau/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| tau/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| tau/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| tau/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| tau/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| tau/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| tau/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| tau/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 2.2e-05 | 2.2e-05 | 0.001 | second formulas | passed |
| tau/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| tau/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| tau/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| tau/X0 | run complete | yes | yes | – | equal | construction | passed |
| tau/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| tau/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| tau/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| tau/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| tau/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| tau/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| tau/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| tau/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| tau/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| tau/X0 | lig_rmsd_A: mean (A) | 0.265179 | 0.265179 | 0 | 0.001 | second formulas | passed |
| tau/X0 | lig_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | lig_rmsd_A: mean_last_half (A) | 0.265179 | 0.265179 | 0 | 0.001 | second formulas | passed |
| tau/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | tau_fold_deg: mean (deg) | 153.056 | 153.056 | 7e-06 | 0.01 | second formulas | passed |
| tau/X0 | tau_fold_deg: sd (deg) | 0 | 7e-05 | 7e-05 | 0.01 | second formulas | passed |
| tau/X0 | tau_fold_deg: mean_last_half (deg) | 153.056 | 153.056 | 4e-06 | 0.01 | second formulas | passed |
| tau/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 6.8e-05 | 6.8e-05 | 0.01 | second formulas | passed |
| tau/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| tau/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| tau/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.3e-05 | 5.3e-05 | 0.01 | second formulas | passed |
| tau/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| tau/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 4.9e-05 | 4.9e-05 | 0.01 | second formulas | passed |
| tau/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| tau/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| tau/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| tau/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d1_N1_F131CZ_A: mean (A) | 3.90553 | 3.90553 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| tau/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.90553 | 3.90553 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d2_N1_N67OD1_A: mean (A) | 7.75215 | 7.75215 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| tau/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.75215 | 7.75215 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | n2_volume: mean | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| tau/X0 | n2_volume: sd | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| tau/X0 | n2_volume: mean_last_half | 0.451372 | 0.451372 | 0 | 0.001 | second formulas | passed |
| tau/X0 | n2_volume: sd_last_half | 0.545234 | 0.545234 | 0 | 0.001 | second formulas | passed |
| tau/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE1_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE1_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE_min_A: sd (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| tau/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| tau/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| tau/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 0 | 0.001 | construction | passed |
| tau/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| tau/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| tau/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| tau/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| tau/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| tau/X0 | gate G-V: the criteria that are not met | GV4 tau_fold_deg | GV4 tau_fold_deg | – | equal | second formulas | passed |
| tau/X0 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| tau/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| tau/X0 | gate G-V: GV4 as a whole | no | no | – | equal | second formulas | passed |
| tau/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| tau/X0 | gate G-V: verdict (none for a run that is not X0) | FAILED | FAILED | – | equal | second formulas | passed |
| tau/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| tau/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823967 | 3e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840544 | 1e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776332 | 2e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731089 | 2e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807205 | 0 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005404 | 2e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 20.9561 | 20.9561 | 7e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| tau/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| tau/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| tau/X0 | gate G-V by design: verdict | FAILED | FAILED | – | equal | construction | passed |
| tau/X0 | gate G-V by design: the criteria that are not met | GV4 tau_fold_deg | GV4 tau_fold_deg | – | equal | construction | passed |
| tau/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| tau/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| tau/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| tau/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| tau/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| tau/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| tau/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| tau/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| tau/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| tau/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| tau/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| tau/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| tau/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| tau/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| tau/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| tau/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| tau/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| tau/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| tau/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| tau/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| tau/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| tau/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| tau/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| tau/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| tau/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| tau/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| tau/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| tau/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| tau/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| tau/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| tau/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| tau/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47008 | 6e-06 | 0.001 | second formulas | passed |
| tau/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| tau/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.007685, 4.003721, 0.0] | [4.007685, 4.003721, 0.0] | – | equal | second formulas | passed |
| tau/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 2.99678 | 2.99677 | 8e-06 | 0.001 | second formulas | passed |
| tau/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.77905 | 3.77905 | 2e-06 | 0.001 | second formulas | passed |
| tau/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| tau/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| tau/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| tau/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02842 | 4.02841 | 1e-05 | 0.001 | second formulas | passed |
| tau/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.84239 | 4.84239 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.028406] | [0.0, 4.028406] | – | equal | second formulas | passed |
| tau/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99677 | 2.99677 | – | equal | second formulas | passed |
| tau/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| tau/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| tau/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| tau/X0 | first production frame: signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| tau/X0 | first production frame: ligand RMSD (A) | 0.265179 | 0.265179 | 0 | 0.001 | second formulas | passed |
| tau/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| tau/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| tau/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| tau/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| tau/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| water/X0 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| water/X0 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| water/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| water/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| water/X0 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| water/X0 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| water/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| water/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| water/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| water/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| water/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| water/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| water/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000999 | 0.000999 | 0.01 | second formulas | passed |
| water/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000512 | 0.000512 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000354 | 0.000354 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000437 | 0.000437 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000678 | 0.000678 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000206 | 0.000206 | 0.01 | second formulas | passed |
| water/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| water/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| water/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| water/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| water/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| water/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| water/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | every frame: n2_volume, largest difference | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| water/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| water/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| water/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| water/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| water/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| water/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| water/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| water/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| water/X0 | run complete | yes | yes | – | equal | construction | passed |
| water/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| water/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| water/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| water/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| water/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| water/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| water/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| water/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| water/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| water/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_N2_A: mean (A) | 1.95917 | 1.95917 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_N2_A: sd (A) | 0.093528 | 0.093528 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_N2_A: mean_last_half (A) | 1.95962 | 1.95962 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_N2_A: sd_last_half (A) | 0.092838 | 0.092838 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | tau_fold_deg: sd (deg) | 0 | 7.3e-05 | 7.3e-05 | 0.01 | second formulas | passed |
| water/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 0 | 0.01 | second formulas | passed |
| water/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 7.3e-05 | 7.3e-05 | 0.01 | second formulas | passed |
| water/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| water/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.1e-05 | 5.1e-05 | 0.01 | second formulas | passed |
| water/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.1e-05 | 5.1e-05 | 0.01 | second formulas | passed |
| water/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| water/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| water/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| water/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| water/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| water/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| water/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| water/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | n2_volume: mean | 0.45178 | 0.45178 | 0 | 0.001 | second formulas | passed |
| water/X0 | n2_volume: sd | 0.545759 | 0.545759 | 0 | 0.001 | second formulas | passed |
| water/X0 | n2_volume: mean_last_half | 0.451405 | 0.451405 | 0 | 0.001 | second formulas | passed |
| water/X0 | n2_volume: sd_last_half | 0.545846 | 0.545846 | 0 | 0.001 | second formulas | passed |
| water/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE1_A: mean (A) | 4.00998 | 4.00998 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE1_A: sd (A) | 0.025401 | 0.025402 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00952 | 4.00952 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.02335 | 0.02335 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE2_A: mean (A) | 5.4772 | 5.4772 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE2_A: sd (A) | 0.033054 | 0.033054 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47658 | 5.47658 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.030068 | 0.030068 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE_min_A: mean (A) | 4.00998 | 4.00998 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE_min_A: sd (A) | 0.025401 | 0.025402 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00952 | 4.00952 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.02335 | 0.02335 | 0 | 0.001 | second formulas | passed |
| water/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| water/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| water/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 1e-06 | 0.001 | construction | passed |
| water/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| water/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| water/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| water/X0 | frame 8, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| water/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| water/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| water/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| water/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| water/X0 | gate G-V: the criteria that are not met | GV3 fifth ligand, full trajectory | GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| water/X0 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| water/X0 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| water/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| water/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| water/X0 | gate G-V: verdict (none for a run that is not X0) | FAILED | FAILED | – | equal | second formulas | passed |
| water/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| water/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.73931 | -0.73931 | 0 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.955482 | -0.955488 | 6e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.647526 | -0.647524 | 2e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.573196 | 0.573198 | 2e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.919415 | 0.919418 | 3e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.67014 | 0.670139 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005403 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.04393 | 1e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 5 | 5 | 0 | 0 | second formulas | passed |
| water/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 3 | 3 | 0 | 0 | second formulas | passed |
| water/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| water/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.955482 | 0.955488 | 6e-06 | 0.01 | second formulas | passed |
| water/X0 | gate G-V by design: verdict | FAILED | FAILED | – | equal | construction | passed |
| water/X0 | gate G-V by design: the criteria that are not met | GV3 fifth ligand, full trajectory | GV3 fifth ligand, full trajectory | – | equal | construction | passed |
| water/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| water/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 5 | 5 | 0 | 0 | construction | passed |
| water/X0 | by design: full frames with any oxygen within 2.6 A | 3 | 3 | 0 | 0 | construction | passed |
| water/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| water/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| water/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| water/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| water/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| water/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| water/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| water/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| water/X0 | N2: frames with positive, negative and zero signed volume | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| water/X0 | N2: fraction of frames with positive signed volume | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| water/X0 | N2: fraction of frames with negative signed volume | 0.2 | 0.2 | 0 | 1e-06 | second formulas | passed |
| water/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| water/X0 | N2: positive, negative and zero in the last 20 ns | [400, 100, 0] | [400, 100, 0] | – | equal | second formulas | passed |
| water/X0 | N2: fraction positive in the last 20 ns | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| water/X0 | N2, frames of state C: number, positive, negative, zero | [500, 400, 100, 0] | [500, 400, 100, 0] | – | equal | second formulas | passed |
| water/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 0.8 | 0.8 | 0 | 1e-06 | second formulas | passed |
| water/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| water/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| water/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| water/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| water/X0 | N2: stretches of one sign, and those of at least 1 ns | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| water/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [] | [] | – | equal | second formulas | passed |
| water/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| water/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| water/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| water/X0 | N2: sign of the first and of the last frame | [1, -1] | [1, -1] | – | equal | second formulas | passed |
| water/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| water/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| water/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| water/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| water/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| water/X0 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| water/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| water/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00998 | 4.00998 | 0 | 0.001 | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 5e-06 | 0.001 | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.4772 | 5.4772 | 0 | 0.001 | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| water/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| water/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.009982, 4.003722, 0.0] | [4.009982, 4.003722, 0.0] | – | equal | second formulas | passed |
| water/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| water/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.4e-05 | 0.001 | second formulas | passed |
| water/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.78 | 3.78 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| water/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| water/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| water/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| water/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| water/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 4.02013 | 4.02012 | 1.1e-05 | 0.001 | second formulas | passed |
| water/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 4.83727 | 4.83727 | 2e-06 | 0.001 | second formulas | passed |
| water/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 4.020122] | [0.0, 4.020122] | – | equal | second formulas | passed |
| water/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| water/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| water/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| water/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| water/X0 | first production frame: signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| water/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| water/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| water/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.8, 0.8, 0] | [0.8, 0.8, 0] | – | equal | second formulas | passed |
| water/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| water/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| water/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| test/X9 (overlap) | solute frames seen | 250 | 250 | 0 | 0 | construction | passed |
| test/X9 (overlap) | DCD headers: segments contiguous (none = the headers tell nothing) | no | no | – | equal | construction | passed |
| test/X9 (overlap) | manifest lists the segments | yes | yes | – | equal | construction | passed |
| test/X9 (overlap) | manifest and files agree (none = no list) | no | no | – | equal | construction | passed |
| test/X9 (overlap) | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| test/X9 (overlap) | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| test/X9 (overlap) | selections.json is the file of the manifest (none = no checksum there) | no | no | – | equal | construction | passed |
| test/X9 (overlap) | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| test/X9 (overlap) | repeated frames found at the segment boundary | yes | yes | – | equal | construction | passed |
| test/X9 (overlap) | frame of the first segment that is repeated | 150 | 150 | 0 | 0 | construction | passed |
| test/X9 (overlap) | warnings written (headers, repeated frame, manifest, system file) | 4 | 4 | 0 | 0 | construction | passed |
| conv/R1 | solute frames seen | 2200 | 2200 | 0 | 0 | construction | passed |
| conv/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | yes | yes | – | equal | construction | passed |
| conv/R1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| conv/R1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| conv/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| conv/R1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| conv/R1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| conv/R1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| conv/R1 | full frames seen | 220 | 220 | 0 | 0 | construction | passed |
| conv/R1 | time seen (ns) | 22 | 22 | 0 | 1e-09 | construction | passed |
| conv/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| conv/R1 | frames in which the ligand was in another periodic image | 314 | 314 | 0 | 0 | construction | passed |
| conv/R1 | frames in which the zinc was in another periodic image | 200 | 200 | 0 | 0 | construction | passed |
| conv/R1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | lines of timeseries_solute.csv | 2200 | 2200 | 0 | 0 | construction | passed |
| conv/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| conv/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000872 | 0.000872 | 0.01 | second formulas | passed |
| conv/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000847 | 0.000847 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000755 | 0.000755 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| conv/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| conv/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.9e-05 | 5.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000393 | 0.000393 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 0.000474 | 0.000474 | 0.01 | second formulas | passed |
| conv/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| conv/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| conv/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| conv/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| conv/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| conv/R1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | lines of timeseries_full.csv | 220 | 220 | 0 | 0 | construction | passed |
| conv/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| conv/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| conv/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| conv/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | verdict | MIXED | MIXED | – | equal | construction | passed |
| conv/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| conv/R1 | run complete | yes | yes | – | equal | construction | passed |
| conv/R1 | first passage into the opposite state (ns) | 6.01 | 6.01 | 0 | 1e-09 | construction | passed |
| conv/R1 | transitions | 1 | 1 | 0 | 0 | construction | passed |
| conv/R1 | long visits (at least 1 ns) | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | frames of the verdict window | 2000 | 2000 | 0 | 0 | construction | passed |
| conv/R1 | count_C_window | 1600 | 1600 | 0 | 0 | construction | passed |
| conv/R1 | count_R_window | 400 | 400 | 0 | 0 | construction | passed |
| conv/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | count_C | 1600 | 1600 | 0 | 0 | construction | passed |
| conv/R1 | count_R | 600 | 600 | 0 | 0 | construction | passed |
| conv/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| conv/R1 | count_C (design against the planted frames) | 1600 | 1600 | 0 | 0 | construction | passed |
| conv/R1 | count_R (design against the planted frames) | 600 | 600 | 0 | 0 | construction | passed |
| conv/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | frames with the hydrogen bond N2-H...OG1 | 1101 | 1101 | 0 | 0 | second formulas | passed |
| conv/R1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| conv/R1 | frames with labelling identity | 1000 | 1000 | 0 | 0 | second formulas | passed |
| conv/R1 | frames with labelling ring_flip | 400 | 400 | 0 | 0 | second formulas | passed |
| conv/R1 | frames with labelling oxygen_swap | 400 | 400 | 0 | 0 | second formulas | passed |
| conv/R1 | frames with labelling ring_flip+oxygen_swap | 400 | 400 | 0 | 0 | second formulas | passed |
| conv/R1 | lig_rmsd_A: mean (A) | 0.974516 | 0.974517 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | lig_rmsd_A: sd (A) | 1.59138 | 1.59138 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | lig_rmsd_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_N2_A: mean (A) | 2.04427 | 2.04427 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_N2_A: sd (A) | 0.162465 | 0.162465 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_N2_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 3.36364e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | tau_fold_deg: mean (deg) | 134.971 | 134.971 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | tau_fold_deg: sd (deg) | 4.76021 | 4.76021 | 0 | 0.01 | second formulas | passed |
| conv/R1 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 2e-06 | 0.01 | second formulas | passed |
| conv/R1 | tau_fold_deg: sd_last_half (deg) | 0 | 6.6e-05 | 6.6e-05 | 0.01 | second formulas | passed |
| conv/R1 | tau_fold_deg: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 1.95778e-06 | 2e-06 | 0.01 | second formulas | passed |
| conv/R1 | ang_C4_S_N2_deg: mean (deg) | 108.968 | 108.968 | 0 | 0.01 | second formulas | passed |
| conv/R1 | ang_C4_S_N2_deg: sd (deg) | 0.28384 | 0.283842 | 2e-06 | 0.01 | second formulas | passed |
| conv/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 2e-06 | 0.01 | second formulas | passed |
| conv/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| conv/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 2.49546e-06 | 2e-06 | 0.01 | second formulas | passed |
| conv/R1 | d_S_N2_A: mean (A) | 1.6184 | 1.6184 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_S_N2_A: sd (A) | 0.036038 | 0.036038 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | d_S_N2_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 3.27904e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| conv/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| conv/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| conv/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| conv/R1 | ca_rmsd_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 1.98745e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d1_N1_F131CZ_A: mean (A) | 5.44797 | 5.44797 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d1_N1_F131CZ_A: sd (A) | 2.53337 | 2.53337 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 4.73901e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d2_N1_N67OD1_A: mean (A) | 6.34044 | 6.34044 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d2_N1_N67OD1_A: sd (A) | 2.29251 | 2.29251 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 1.18924e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | n2_volume: mean | 0.461309 | 0.461309 | 0 | 0.001 | second formulas | passed |
| conv/R1 | n2_volume: sd | 0.553214 | 0.553214 | 0 | 0.001 | second formulas | passed |
| conv/R1 | n2_volume: mean_last_half | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| conv/R1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | n2_volume: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | n2_volume: block means, largest difference | 0 | 4.64274e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE1_A: mean (A) | 3.96547 | 3.96547 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE1_A: sd (A) | 0.24701 | 0.24701 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 3.96457 | 3.96457 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.246667 | 0.246667 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 4.93638e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE2_A: mean (A) | 5.4455 | 5.4455 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE2_A: sd (A) | 0.291357 | 0.291357 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.44487 | 5.44487 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.291189 | 0.291189 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 4.28772e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE_min_A: mean (A) | 3.95093 | 3.95093 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE_min_A: sd (A) | 0.284295 | 0.284295 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 3.95004 | 3.95004 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.283951 | 0.283951 | 0 | 0.001 | second formulas | passed |
| conv/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 2 | 2 | 0 | 0 | construction | passed |
| conv/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 2.27247e-07 | 0 | 0.001 | second formulas | passed |
| conv/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| conv/R1 | nearest water to the zinc, mean (A) | 3.68414 | 3.68414 | 0 | 0.001 | construction | passed |
| conv/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| conv/R1 | frame 2, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| conv/R1 | frame 2, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| conv/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (A) | 3.31 | 3.3102 | 0.0002 | 0.0015 | starts.csv | passed |
| conv/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| conv/R1 | planned length (X0: 5 ns by the specification) (ns) | 22 | 22 | 0 | 1e-09 | construction | passed |
| conv/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| conv/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| conv/R1 | gate G-V: the criteria that are not met | GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | GV2 hydrogen bond, GV3 Zn-N2, GV3 fifth ligand, full trajectory | – | equal | second formulas | passed |
| conv/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| conv/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| conv/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| conv/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | no | no | – | equal | second formulas | passed |
| conv/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| conv/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| conv/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 1.20688 | 1.20688 | 0 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -1.3091 | -1.3091 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -2.51086 | -2.51086 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.765746 | 0.765746 | 0 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.939162 | 0.939163 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.763683 | 0.763683 | 0 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.168414 | 0.168414 | 0 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 2.87109 | 2.87109 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | gate G-V: solute frames counted for that fraction | 2200 | 2200 | 0 | 0 | construction | passed |
| conv/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 21 | 21 | 0 | 0 | second formulas | passed |
| conv/R1 | gate G-V: full frames counted for that fraction | 220 | 220 | 0 | 0 | construction | passed |
| conv/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 2.51086 | 2.51086 | 1e-06 | 0.01 | second formulas | passed |
| conv/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| conv/R1 | frames with ND2 nearer to the aniline N than OD1 | 1600 | 1600 | 0 | 0 | second formulas | passed |
| conv/R1 | state C: frames, and frames with ND2 nearer | [1600, 1600] | [1600, 1600] | – | equal | second formulas | passed |
| conv/R1 | state R: frames, and frames with ND2 nearer | [600, 0] | [600, 0] | – | equal | second formulas | passed |
| conv/R1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| conv/R1 | campaign-table value: fraction of frames with ND2 nearer | 0.727273 | 0.727273 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| conv/R1 | N2: frames with positive, negative and zero signed volume | [1799, 401, 0] | [1799, 401, 0] | – | equal | second formulas | passed |
| conv/R1 | N2: fraction of frames with positive signed volume | 0.817727 | 0.817727 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2: fraction of frames with negative signed volume | 0.182273 | 0.182273 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2: frames of the window (last 20 ns) | 2000 | 2000 | 0 | 0 | construction | passed |
| conv/R1 | N2: positive, negative and zero in the last 20 ns | [1799, 201, 0] | [1799, 201, 0] | – | equal | second formulas | passed |
| conv/R1 | N2: fraction positive in the last 20 ns | 0.8995 | 0.8995 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2, frames of state C: number, positive, negative, zero | [1600, 1600, 0, 0] | [1600, 1600, 0, 0] | – | equal | second formulas | passed |
| conv/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2, frames of state R: number, positive, negative, zero | [600, 199, 401, 0] | [600, 199, 401, 0] | – | equal | second formulas | passed |
| conv/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0.331667 | 0.331667 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| conv/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| conv/R1 | N2: stretches of one sign, and those of at least 1 ns | [6, 5] | [6, 5] | – | equal | second formulas | passed |
| conv/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[-1, 0, 200], [-1, 299, 101], [1, 400, 100], [-1, 500, 100], [1, 600, 1600]] | [[-1, 0, 200], [-1, 299, 101], [1, 400, 100], [-1, 500, 100], [1, 600, 1600]] | – | equal | second formulas | passed |
| conv/R1 | N2: number of sign changes that last at least 1 ns | 3 | 3 | 0 | 0 | second formulas | passed |
| conv/R1 | N2: these sign changes (frame, from, to) | [[400, -1, 1], [500, 1, -1], [600, -1, 1]] | [[400, -1, 1], [500, 1, -1], [600, -1, 1]] | – | equal | second formulas | passed |
| conv/R1 | N2: times of these sign changes (ns) | [4.01, 5.01, 6.01] | [4.01, 5.01, 6.01] | – | equal | second formulas | passed |
| conv/R1 | N2: time of the first lasting sign change (ns) | 4.01 | 4.01 | 0 | 1e-09 | second formulas | passed |
| conv/R1 | N2: sign of the first and of the last frame | [-1, 1] | [-1, 1] | – | equal | second formulas | passed |
| conv/R1 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| conv/R1 | hydrogen bond around the first passage: frame of the first passage | 600 | 600 | – | equal | second formulas | passed |
| conv/R1 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [600, 301] | [600, 301] | – | equal | second formulas | passed |
| conv/R1 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [1600, 800] | [1600, 800] | – | equal | second formulas | passed |
| conv/R1 | hydrogen bond N2-H...OG1, fraction before the first passage | 0.501667 | 0.501667 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.5 | 0.5 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.501667, 0.5] | [0.501667, 0.5] | – | equal | construction | passed |
| conv/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [44, 22, 66] | [44, 22, 66] | – | equal | second formulas | passed |
| conv/R1 | Glu106: smallest distance zinc to OE1 (A) | 2.5 | 2.5 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | Glu106: smallest distance zinc to OE2 (A) | 2.55 | 2.55 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 20 | 20 | 0 | 0 | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 3.96547 | 3.96547 | 0 | 0.001 | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 2.5 | 2.5 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [44, 2200] | [44, 2200] | – | equal | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.4455 | 5.4455 | 0 | 0.001 | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 2.55 | 2.55 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [22, 2200] | [22, 2200] | – | equal | second formulas | passed |
| conv/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [3.965468, 2.499999, 0.02] | [3.965468, 2.499999, 0.02] | – | equal | second formulas | passed |
| conv/R1 | water at the crystal ring positions: mean number within 3.5 A | 1.36364 | 1.36364 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 1.38237 | 1.38236 | 8e-06 | 0.001 | second formulas | passed |
| conv/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.05294 | 3.05294 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.36364 | 1.36364 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| conv/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| conv/R1 | water in the place of the ring: mean number | 0.090909 | 0.090909 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | water in the place of the ring: largest number in a full frame | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| conv/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 1.38237 | 1.38236 | 8e-06 | 0.001 | second formulas | passed |
| conv/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 3.3984 | 3.3984 | 0 | 0.001 | second formulas | passed |
| conv/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.090909, 1.382364] | [0.090909, 1.382364] | – | equal | second formulas | passed |
| conv/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 1.38236 | 1.38236 | – | equal | second formulas | passed |
| conv/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| conv/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| conv/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| conv/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| conv/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.817727, 0.8995, 3] | [0.817727, 0.8995, 3] | – | equal | second formulas | passed |
| conv/R1 | equilibration: frames read (files solute_eqNN.dcd present) | yes | yes | – | equal | construction | passed |
| conv/R1 | equilibration: time series written (only with equilibration files) | [True, True] | [True, True] | – | equal | construction | passed |
| conv/R1 | campaign-table value: solute frames of the equilibration | 200 | 200 | 0 | 0 | construction | passed |
| conv/R1 | equilibration: solute frames and full frames read | [200, 20] | [200, 20] | – | equal | construction | passed |
| conv/R1 | equilibration: length of the equilibration (ps) | 2000 | 2000 | 0 | 1e-09 | construction | passed |
| conv/R1 | equilibration: where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| conv/R1 | equilibration: lines of timeseries_equil_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| conv/R1 | equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned) | yes | yes | – | equal | construction | passed |
| conv/R1 | equilibration, every frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | equilibration, every frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| conv/R1 | equilibration: no time above zero | yes | yes | – | equal | construction | passed |
| conv/R1 | equilibration, every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000269 | 0.000269 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_O2_A, largest difference (A) | 0 | 4.9e-05 | 4.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 4.9e-05 | 4.9e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration: lines of timeseries_equil_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| conv/R1 | equilibration, every full frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | equilibration, every full frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| conv/R1 | equilibration, every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| conv/R1 | equilibration, every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| conv/R1 | equilibration: the steps with a file, in the order of the run | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 1: start and end (ps) | [0.0, 100.0] | [0.0, 100.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 1: solute and full frames read, and expected | [10, 10, 1, 1] | [10, 10, 1, 1] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 10, 0, 0, 10] | [0, 10, 0, 0, 10] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 1: last solute frame (frame, time, state, sign at N2) | [9, -1.9, 'R', -1] | [9, -1.9, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 1: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 1: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [0, -1.9, 0] | [0, -1.9, 0] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 3e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 3e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 2: start and end (ps) | [100.0, 400.0] | [100.0, 400.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 2: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 2: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 2: last solute frame (frame, time, state, sign at N2) | [39, -1.6, 'R', -1] | [39, -1.6, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 2: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 2: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [3, -1.6, 0] | [3, -1.6, 0] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 2: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 3: start and end (ps) | [400.0, 700.0] | [400.0, 700.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 3: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 3: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 3: last solute frame (frame, time, state, sign at N2) | [69, -1.3, 'R', -1] | [69, -1.3, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 3: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 3: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 3: last full frame (frame, time, water oxygens at the crystal ring positions) | [6, -1.3, 1] | [6, -1.3, 1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 3: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.59999 | 5e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.4 | 3.4 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 4: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 4: start and end (ps) | [700.0, 1100.0] | [700.0, 1100.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 4: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 4: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 4: last solute frame (frame, time, state, sign at N2) | [109, -0.9, 'R', -1] | [109, -0.9, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 4: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 4: last solute frame, signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 4: last full frame (frame, time, water oxygens at the crystal ring positions) | [10, -0.9, 2] | [10, -0.9, 2] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 4: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the place of the ring (A) | 1.38237 | 1.38237 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 1.38237 | 1.38237 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 5: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 5: start and end (ps) | [1100.0, 1500.0] | [1100.0, 1500.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 5: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 5: last solute frame (frame, time, state, sign at N2) | [149, -0.5, 'R', -1] | [149, -0.5, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 5: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 5: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 5: last full frame (frame, time, water oxygens at the crystal ring positions) | [14, -0.5, 2] | [14, -0.5, 2] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 5: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the place of the ring (A) | 1.38237 | 1.38237 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 1.38237 | 1.38237 | 2e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 6: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 6: start and end (ps) | [1500.0, 2000.0] | [1500.0, 2000.0] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 6: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| conv/R1 | equilibration, stage 6: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 6: last solute frame (frame, time, state, sign at N2) | [199, 0.0, 'R', -1] | [199, 0.0, 'R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 6: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 3e-06 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 6: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 6: last full frame (frame, time, water oxygens at the crystal ring positions) | [19, 0.0, 2] | [19, 0.0, 2] | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 6: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| conv/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the place of the ring (A) | 1.38237 | 1.38237 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 1.38237 | 1.38237 | 0 | 0.001 | second formulas | passed |
| conv/R1 | equilibration: frames of state C, R, I; positive and negative sign at N2 | [0, 200, 0, 0, 200] | [0, 200, 0, 0, 200] | – | equal | second formulas | passed |
| conv/R1 | equilibration: first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| conv/R1 | equilibration: first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| conv/R1 | equilibration: first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| conv/R1 | equilibration: state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | second formulas | passed |
| conv/R1 | equilibration: time of the last solute frame (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| conv/R1 | equilibration: frames in which the ligand and the zinc were in another periodic image | [29, 18] | [29, 18] | – | equal | construction | passed |
| conv/R1 | by design: frames with positive and with negative signed volume at N2 | [1799, 401] | [1799, 401] | – | equal | construction | passed |
| conv/R1 | by design: positive and negative in the last 20 ns | [1799, 201] | [1799, 201] | – | equal | construction | passed |
| conv/R1 | by design, frames of state C: number, positive, negative | [1600, 1600, 0] | [1600, 1600, 0] | – | equal | construction | passed |
| conv/R1 | by design, frames of state R: number, positive, negative | [600, 199, 401] | [600, 199, 401] | – | equal | construction | passed |
| conv/R1 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| conv/R1 | by design: stretches of one sign, and those of at least 1 ns | [6, 5] | [6, 5] | – | equal | construction | passed |
| conv/R1 | by design: times of the sign changes that last at least 1 ns (ns) | [4.01, 5.01, 6.01] | [4.01, 5.01, 6.01] | – | equal | construction | passed |
| conv/R1 | by design: sign after each of these changes | [1, -1, 1] | [1, -1, 1] | – | equal | construction | passed |
| conv/R1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [600, 301, 1600, 800] | [600, 301, 1600, 800] | – | equal | construction | passed |
| conv/R1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [44, 22, 66] | [44, 22, 66] | – | equal | construction | passed |
| conv/R1 | by design: smallest distance zinc to OE1 of Glu106 (A) | 2.5 | 2.5 | 1e-06 | 0.001 | construction | passed |
| conv/R1 | by design: smallest distance zinc to OE2 of Glu106 (A) | 2.55 | 2.55 | 1e-06 | 0.001 | construction | passed |
| conv/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.36364 | 1.36364 | 0 | 1e-06 | construction | passed |
| conv/R1 | by design: smallest distance of a water oxygen to the crystal ring positions (A) | 1.382 | 1.38236 | 0.000364 | 0.001 | construction | passed |
| conv/R1 | by design: water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] | – | equal | construction | passed |
| conv/R1 | by design: mean number of water oxygens in the place of the ring | 0.090909 | 0.090909 | 0 | 1e-06 | construction | passed |
| conv/R1 | by design: smallest distance of a water oxygen to the place of the ring (A) | 1.382 | 1.38236 | 0.000364 | 0.001 | construction | passed |
| conv/R1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['R', False, -1, False] | ['R', False, -1, False] | – | equal | construction | passed |
| conv/R1 | by design: equilibration frames read | yes | yes | – | equal | construction | passed |
| conv/R1 | by design: equilibration, solute frames, full frames, length (ps) | [200, 20, 2000.0] | [200, 20, 2000.0] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| conv/R1 | by design: equilibration, solute and full frames of each stage | [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]] | [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | [-1.99, 0.0] | [-1.99, 0.0] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, time of solute frame 10 (the first of stage 2) (ns) | [-1.89, '2'] | [-1.89, '2'] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, water oxygens at the crystal ring positions in every full frame | [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] | [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, water oxygens at the crystal ring positions at the end of stage 3 (the held part) | 1 | 1 | 0 | 0 | construction | passed |
| conv/R1 | by design: equilibration, water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | – | equal | construction | passed |
| conv/R1 | by design: equilibration, water oxygens in the place of the ring at the end of stage 3 | 0 | 0 | 0 | 0 | construction | passed |
| conv/R1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| conv/R1 | by design: equilibration, first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| conv/R1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| conv/R1 | by design: equilibration, state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | construction | passed |
| ext/R1 | solute frames seen | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| ext/R1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| ext/R1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| ext/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| ext/R1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| ext/R1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| ext/R1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| ext/R1 | full frames seen | 20 | 20 | 0 | 0 | construction | passed |
| ext/R1 | time seen (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| ext/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| ext/R1 | frames in which the ligand was in another periodic image | 29 | 29 | 0 | 0 | construction | passed |
| ext/R1 | frames in which the zinc was in another periodic image | 18 | 18 | 0 | 0 | construction | passed |
| ext/R1 | warnings | 1 | 1 | 0 | 0 | construction | passed |
| ext/R1 | lines of timeseries_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| ext/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| ext/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000731 | 0.000731 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000434 | 0.000434 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| ext/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| ext/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| ext/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| ext/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| ext/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| ext/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| ext/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | lines of timeseries_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| ext/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| ext/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| ext/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| ext/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | verdict | STAYED | STAYED | – | equal | construction | passed |
| ext/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| ext/R1 | run complete | yes | yes | – | equal | construction | passed |
| ext/R1 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| ext/R1 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| ext/R1 | frames of the verdict window | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | count_C_window | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | count_R_window | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | count_C | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | count_R | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | rule 1: left the crystal pose | yes | yes | – | equal | construction | passed |
| ext/R1 | count_C (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | count_R (design against the planted frames) | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | frames with the hydrogen bond N2-H...OG1 | 200 | 200 | 0 | 0 | second formulas | passed |
| ext/R1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| ext/R1 | frames with labelling identity | 200 | 200 | 0 | 0 | second formulas | passed |
| ext/R1 | frames with labelling ring_flip | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | frames with labelling oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | frames with labelling ring_flip+oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | lig_rmsd_A: mean (A) | 3.57322 | 3.57322 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | lig_rmsd_A: mean_last_half (A) | 3.57322 | 3.57322 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_N2_A: mean (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_N2_A: mean_last_half (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | tau_fold_deg: mean (deg) | 142.744 | 142.744 | 2e-06 | 0.01 | second formulas | passed |
| ext/R1 | tau_fold_deg: sd (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| ext/R1 | tau_fold_deg: mean_last_half (deg) | 142.744 | 142.744 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | tau_fold_deg: sd_last_half (deg) | 0 | 5.8e-05 | 5.8e-05 | 0.01 | second formulas | passed |
| ext/R1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| ext/R1 | ang_C4_S_N2_deg: mean (deg) | 109.432 | 109.432 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.3e-05 | 5.3e-05 | 0.01 | second formulas | passed |
| ext/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 109.432 | 109.432 | 4e-06 | 0.01 | second formulas | passed |
| ext/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| ext/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| ext/R1 | d_S_N2_A: mean (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d_S_N2_A: mean_last_half (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| ext/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| ext/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d1_N1_F131CZ_A: mean (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d2_N1_N67OD1_A: mean (A) | 2.59678 | 2.59678 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d2_N1_N67OD1_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 2.59678 | 2.59678 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | n2_volume: mean | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | n2_volume: sd | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | n2_volume: mean_last_half | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE1_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE1_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE2_A: mean (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE2_A: sd (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE_min_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE_min_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| ext/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| ext/R1 | nearest water to the zinc, mean (A) | 3.6255 | 3.6255 | 2e-06 | 0.001 | construction | passed |
| ext/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| ext/R1 | frame 2, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| ext/R1 | frame 2, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| ext/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (A) | 3.31 | 3.3102 | 0.0002 | 0.0015 | starts.csv | passed |
| ext/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| ext/R1 | planned length (X0: 5 ns by the specification) (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| ext/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| ext/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| ext/R1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | – | equal | second formulas | passed |
| ext/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| ext/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| ext/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| ext/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| ext/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| ext/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| ext/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 6.62247 | 6.62246 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -2.55859 | -2.55859 | 2e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.13628 | -7.13628 | 4e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.85817 | 0.858173 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.961814 | 0.961815 | 1e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.647624 | 0.647627 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.631924 | 0.631921 | 3e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 10.6445 | 10.6445 | 2e-06 | 0.01 | second formulas | passed |
| ext/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | gate G-V: solute frames counted for that fraction | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| ext/R1 | gate G-V: full frames counted for that fraction | 20 | 20 | 0 | 0 | construction | passed |
| ext/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 7.13628 | 7.13628 | 4e-06 | 0.01 | second formulas | passed |
| ext/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| ext/R1 | frames with ND2 nearer to the aniline N than OD1 | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | state C: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| ext/R1 | state R: frames, and frames with ND2 nearer | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| ext/R1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| ext/R1 | campaign-table value: fraction of frames with ND2 nearer | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| ext/R1 | N2: frames with positive, negative and zero signed volume | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| ext/R1 | N2: fraction of frames with positive signed volume | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | N2: fraction of frames with negative signed volume | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | N2: frames of the window (last 20 ns) | 200 | 200 | 0 | 0 | construction | passed |
| ext/R1 | N2: positive, negative and zero in the last 20 ns | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| ext/R1 | N2: fraction positive in the last 20 ns | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | N2, frames of state C: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| ext/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| ext/R1 | N2, frames of state R: number, positive, negative, zero | [200, 0, 200, 0] | [200, 0, 200, 0] | – | equal | second formulas | passed |
| ext/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| ext/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| ext/R1 | N2: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| ext/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[-1, 0, 200]] | [[-1, 0, 200]] | – | equal | second formulas | passed |
| ext/R1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| ext/R1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| ext/R1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| ext/R1 | N2: sign of the first and of the last frame | [-1, -1] | [-1, -1] | – | equal | second formulas | passed |
| ext/R1 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| ext/R1 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| ext/R1 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| ext/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| ext/R1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| ext/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.010967, 4.003724, 0.0] | [4.010967, 4.003724, 0.0] | – | equal | second formulas | passed |
| ext/R1 | water at the crystal ring positions: mean number within 3.5 A | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3.4 | 3.39999 | 7e-06 | 0.001 | second formulas | passed |
| ext/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.4 | 3.4 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'], 3.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'], 3.5] | – | equal | construction | passed |
| ext/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['manifest hydration_check.settings', False] | ['manifest hydration_check.settings', False] | – | equal | construction | passed |
| ext/R1 | water in the place of the ring: mean number | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | water in the place of the ring: largest number in a full frame | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| ext/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.4 | 3.39999 | 7e-06 | 0.001 | second formulas | passed |
| ext/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 3.4 | 3.4 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [1.0, 3.399993] | [1.0, 3.399993] | – | equal | second formulas | passed |
| ext/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 3.39999 | 3.39999 | – | equal | second formulas | passed |
| ext/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| ext/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| ext/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| ext/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| ext/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.0, 0.0, 0] | [0.0, 0.0, 0] | – | equal | second formulas | passed |
| ext/R1 | equilibration: frames read (files solute_eqNN.dcd present) | yes | yes | – | equal | construction | passed |
| ext/R1 | equilibration: time series written (only with equilibration files) | [True, True] | [True, True] | – | equal | construction | passed |
| ext/R1 | campaign-table value: solute frames of the equilibration | 300 | 300 | 0 | 0 | construction | passed |
| ext/R1 | equilibration: solute frames and full frames read | [300, 30] | [300, 30] | – | equal | construction | passed |
| ext/R1 | equilibration: length of the equilibration (ps) | 3000 | 3000 | 0 | 1e-09 | construction | passed |
| ext/R1 | equilibration: where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| ext/R1 | equilibration: lines of timeseries_equil_solute.csv | 300 | 300 | 0 | 0 | construction | passed |
| ext/R1 | equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned) | yes | yes | – | equal | construction | passed |
| ext/R1 | equilibration, every frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | equilibration, every frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| ext/R1 | equilibration: no time above zero | yes | yes | – | equal | construction | passed |
| ext/R1 | equilibration, every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000284 | 0.000284 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000269 | 0.000269 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000566 | 0.000566 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration: lines of timeseries_equil_full.csv | 30 | 30 | 0 | 0 | construction | passed |
| ext/R1 | equilibration, every full frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | equilibration, every full frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| ext/R1 | equilibration, every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| ext/R1 | equilibration, every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration: the steps with a file, in the order of the run | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 1: start and end (ps) | [0.0, 100.0] | [0.0, 100.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 1: solute and full frames read, and expected | [10, 10, 1, 1] | [10, 10, 1, 1] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 10, 0, 0, 10] | [0, 10, 0, 0, 10] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 1: last solute frame (frame, time, state, sign at N2) | [9, -2.9, 'R', -1] | [9, -2.9, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 1: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 1: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [0, -2.9, 0] | [0, -2.9, 0] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 2: start and end (ps) | [100.0, 400.0] | [100.0, 400.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 2: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 2: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 2: last solute frame (frame, time, state, sign at N2) | [39, -2.6, 'R', -1] | [39, -2.6, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 2: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 2: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [3, -2.6, 0] | [3, -2.6, 0] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 2: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [0, True] | [0, True] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3: start and end (ps) | [400.0, 700.0] | [400.0, 700.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3: last solute frame (frame, time, state, sign at N2) | [69, -2.3, 'R', -1] | [69, -2.3, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3: last full frame (frame, time, water oxygens at the crystal ring positions) | [6, -2.3, 0] | [6, -2.3, 0] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [0, True] | [0, True] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 1: start and end (ps) | [700.0, 1200.0] | [700.0, 1200.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 1: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last solute frame (frame, time, state, sign at N2) | [119, -1.8, 'R', -1] | [119, -1.8, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last solute frame, signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [11, -1.8, 0] | [11, -1.8, 0] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 3e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [1, True] | [1, True] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 2: start and end (ps) | [1200.0, 1700.0] | [1200.0, 1700.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 2: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 3, extension step 2: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last solute frame (frame, time, state, sign at N2) | [169, -1.3, 'R', -1] | [169, -1.3, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last solute frame, signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [16, -1.3, 1] | [16, -1.3, 1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last full frame, nearest water oxygen to the place of the ring (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 3, extension step 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 4: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 4: start and end (ps) | [1700.0, 2100.0] | [1700.0, 2100.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 4: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 4: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 4: last solute frame (frame, time, state, sign at N2) | [209, -0.9, 'R', -1] | [209, -0.9, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 4: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 4: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 4: last full frame (frame, time, water oxygens at the crystal ring positions) | [20, -0.9, 1] | [20, -0.9, 1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 4: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the place of the ring (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 5: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 5: start and end (ps) | [2100.0, 2500.0] | [2100.0, 2500.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 5: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 5: last solute frame (frame, time, state, sign at N2) | [249, -0.5, 'R', -1] | [249, -0.5, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 5: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 5: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 5: last full frame (frame, time, water oxygens at the crystal ring positions) | [24, -0.5, 1] | [24, -0.5, 1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 5: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the place of the ring (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.4 | 3.4 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 6: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 6: start and end (ps) | [2500.0, 3000.0] | [2500.0, 3000.0] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 6: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| ext/R1 | equilibration, stage 6: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 6: last solute frame (frame, time, state, sign at N2) | [299, 0.0, 'R', -1] | [299, 0.0, 'R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 6: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 6: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 6: last full frame (frame, time, water oxygens at the crystal ring positions) | [29, 0.0, 1] | [29, 0.0, 1] | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 6: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| ext/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the place of the ring (A) | 3.4 | 3.4 | 5e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.4 | 3.4 | 5e-06 | 0.001 | second formulas | passed |
| ext/R1 | equilibration: frames of state C, R, I; positive and negative sign at N2 | [0, 300, 0, 0, 300] | [0, 300, 0, 0, 300] | – | equal | second formulas | passed |
| ext/R1 | equilibration: first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| ext/R1 | equilibration: first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| ext/R1 | equilibration: first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| ext/R1 | equilibration: state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | second formulas | passed |
| ext/R1 | equilibration: time of the last solute frame (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| ext/R1 | equilibration: frames in which the ligand and the zinc were in another periodic image | [43, 27] | [43, 27] | – | equal | construction | passed |
| ext/R1 | by design: frames with positive and with negative signed volume at N2 | [0, 200] | [0, 200] | – | equal | construction | passed |
| ext/R1 | by design: positive and negative in the last 20 ns | [0, 200] | [0, 200] | – | equal | construction | passed |
| ext/R1 | by design, frames of state C: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| ext/R1 | by design, frames of state R: number, positive, negative | [200, 0, 200] | [200, 0, 200] | – | equal | construction | passed |
| ext/R1 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| ext/R1 | by design: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | construction | passed |
| ext/R1 | by design: times of the sign changes that last at least 1 ns (ns) | [] | [] | – | equal | construction | passed |
| ext/R1 | by design: sign after each of these changes | [] | [] | – | equal | construction | passed |
| ext/R1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [None, None, None, None] | [None, None, None, None] | – | equal | construction | passed |
| ext/R1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| ext/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | by design: mean number of water oxygens at the crystal ring positions | 1 | 1 | 0 | 1e-06 | construction | passed |
| ext/R1 | by design: smallest distance of a water oxygen to the crystal ring positions (A) | 3.4 | 3.39999 | 7e-06 | 0.001 | construction | passed |
| ext/R1 | by design: water oxygens in the place of the ring in every full frame | [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | – | equal | construction | passed |
| ext/R1 | by design: mean number of water oxygens in the place of the ring | 1 | 1 | 0 | 1e-06 | construction | passed |
| ext/R1 | by design: smallest distance of a water oxygen to the place of the ring (A) | 3.4 | 3.39999 | 7e-06 | 0.001 | construction | passed |
| ext/R1 | by design: site atoms and cut-off of the count that the driver acts on (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'], 3.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'], 3.5] | – | equal | construction | passed |
| ext/R1 | by design: where these settings came from, and whether they are the defaults of the specification | ['manifest hydration_check.settings', False] | ['manifest hydration_check.settings', False] | – | equal | construction | passed |
| ext/R1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['R', False, -1, False] | ['R', False, -1, False] | – | equal | construction | passed |
| ext/R1 | by design: equilibration frames read | yes | yes | – | equal | construction | passed |
| ext/R1 | by design: equilibration, solute frames, full frames, length (ps) | [300, 30, 3000.0] | [300, 30, 3000.0] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| ext/R1 | by design: equilibration, solute and full frames of each stage | [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]] | [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | [-2.99, 0.0] | [-2.99, 0.0] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, time of solute frame 10 (the first of stage 2) (ns) | [-2.89, '2'] | [-2.89, '2'] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, water oxygens at the crystal ring positions in every full frame | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, water oxygens at the crystal ring positions at the end of stage 3 (the held part) | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | by design: equilibration, water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, water oxygens in the place of the ring at the end of stage 3 | 0 | 0 | 0 | 0 | construction | passed |
| ext/R1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| ext/R1 | by design: equilibration, first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| ext/R1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| ext/R1 | by design: equilibration, state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, the steps in the order of the run | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, stage, extension step and time of solute frame 70 (the first of eq03x1) | [3, 1, -2.29] | [3, 1, -2.29] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, stage, extension step and time of solute frame 170 (the first of eq04) | [4, 0, -1.29] | [4, 0, -1.29] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, water oxygens in the place of the ring (the count the driver acts on) at the end of eq03, eq03x1, eq03x2 | [0, 0, 1] | [0, 0, 1] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, these steps were compared with the driver (their last full frame is at the end of the step) | [True, True, True] | [True, True, True] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, the count of the driver at the end of these steps agrees with the analysis | [True, True, True] | [True, True, True] | – | equal | construction | passed |
| ext/R1 | by design: equilibration, the nearest water oxygen of the driver at the end of these steps agrees with the analysis | [True, True, True] | [True, True, True] | – | equal | construction | passed |
| flip/X1 | solute frames seen | 300 | 300 | 0 | 0 | construction | passed |
| flip/X1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| flip/X1 | manifest lists the segments | no | no | – | equal | construction | passed |
| flip/X1 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| flip/X1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| flip/X1 | where the system folder was found | convention | convention | – | equal | construction | passed |
| flip/X1 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| flip/X1 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| flip/X1 | full frames seen | 30 | 30 | 0 | 0 | construction | passed |
| flip/X1 | time seen (ns) | 3 | 3 | 0 | 1e-09 | construction | passed |
| flip/X1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| flip/X1 | frames in which the ligand was in another periodic image | 43 | 43 | 0 | 0 | construction | passed |
| flip/X1 | frames in which the zinc was in another periodic image | 27 | 27 | 0 | 0 | construction | passed |
| flip/X1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | lines of timeseries_solute.csv | 300 | 300 | 0 | 0 | construction | passed |
| flip/X1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| flip/X1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 3.9e-05 | 3.9e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000448 | 0.000448 | 0.01 | second formulas | passed |
| flip/X1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.4e-05 | 4.4e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| flip/X1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000588 | 0.000588 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000394 | 0.000394 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000355 | 0.000355 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000758 | 0.000758 | 0.01 | second formulas | passed |
| flip/X1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000567 | 0.000567 | 0.01 | second formulas | passed |
| flip/X1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 6e-05 | 6e-05 | 0.01 | second formulas | passed |
| flip/X1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| flip/X1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| flip/X1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| flip/X1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| flip/X1 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | lines of timeseries_full.csv | 30 | 30 | 0 | 0 | construction | passed |
| flip/X1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| flip/X1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| flip/X1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| flip/X1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | verdict | MIXED | MIXED | – | equal | construction | passed |
| flip/X1 | verdict provisional (run not finished) | yes | yes | – | equal | construction | passed |
| flip/X1 | run complete | no | no | – | equal | construction | passed |
| flip/X1 | first passage into the opposite state (ns) | 0.06 | 0.06 | 0 | 1e-09 | construction | passed |
| flip/X1 | transitions | 2 | 2 | 0 | 0 | construction | passed |
| flip/X1 | long visits (at least 1 ns) | 2 | 2 | 0 | 0 | construction | passed |
| flip/X1 | frames of the verdict window | 300 | 300 | 0 | 0 | construction | passed |
| flip/X1 | count_C_window | 150 | 150 | 0 | 0 | construction | passed |
| flip/X1 | count_R_window | 145 | 145 | 0 | 0 | construction | passed |
| flip/X1 | count_I_window | 5 | 5 | 0 | 0 | construction | passed |
| flip/X1 | count_C | 150 | 150 | 0 | 0 | construction | passed |
| flip/X1 | count_R | 145 | 145 | 0 | 0 | construction | passed |
| flip/X1 | count_I | 5 | 5 | 0 | 0 | construction | passed |
| flip/X1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| flip/X1 | count_C (design against the planted frames) | 150 | 150 | 0 | 0 | construction | passed |
| flip/X1 | count_R (design against the planted frames) | 145 | 145 | 0 | 0 | construction | passed |
| flip/X1 | count_I (design against the planted frames) | 5 | 5 | 0 | 0 | construction | passed |
| flip/X1 | frames with the hydrogen bond N2-H...OG1 | 245 | 245 | 0 | 0 | second formulas | passed |
| flip/X1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| flip/X1 | frames with labelling identity | 114 | 114 | 0 | 0 | second formulas | passed |
| flip/X1 | frames with labelling ring_flip | 110 | 110 | 0 | 0 | second formulas | passed |
| flip/X1 | frames with labelling oxygen_swap | 38 | 38 | 0 | 0 | second formulas | passed |
| flip/X1 | frames with labelling ring_flip+oxygen_swap | 38 | 38 | 0 | 0 | second formulas | passed |
| flip/X1 | lig_rmsd_A: mean (A) | 1.77143 | 1.77143 | 0 | 0.001 | second formulas | passed |
| flip/X1 | lig_rmsd_A: sd (A) | 1.77286 | 1.77286 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_N2_A: mean (A) | 2.10923 | 2.10923 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_N2_A: sd (A) | 0.204071 | 0.204071 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 0 | 0.01 | second formulas | passed |
| flip/X1 | tau_fold_deg: sd (deg) | 0 | 6.1e-05 | 6.1e-05 | 0.01 | second formulas | passed |
| flip/X1 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 6e-06 | 0.01 | second formulas | passed |
| flip/X1 | tau_fold_deg: sd_last_half (deg) | 0 | 6e-05 | 6e-05 | 0.01 | second formulas | passed |
| flip/X1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| flip/X1 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| flip/X1 | ang_C4_S_N2_deg: sd (deg) | 0 | 4.9e-05 | 4.9e-05 | 0.01 | second formulas | passed |
| flip/X1 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 5e-06 | 0.01 | second formulas | passed |
| flip/X1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.2e-05 | 5.2e-05 | 0.01 | second formulas | passed |
| flip/X1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| flip/X1 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| flip/X1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| flip/X1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d1_N1_F131CZ_A: mean (A) | 6.68185 | 6.68185 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d1_N1_F131CZ_A: sd (A) | 2.80563 | 2.80563 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d2_N1_N67OD1_A: mean (A) | 5.25767 | 5.25767 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d2_N1_N67OD1_A: sd (A) | 2.54045 | 2.54045 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | n2_volume: mean | 0.12303 | 0.12303 | 0 | 0.001 | second formulas | passed |
| flip/X1 | n2_volume: sd | 0.601049 | 0.601049 | 0 | 0.001 | second formulas | passed |
| flip/X1 | n2_volume: mean_last_half | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| flip/X1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE1_A: mean (A) | 4.00901 | 4.00901 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE1_A: sd (A) | 0.01807 | 0.01807 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE2_A: mean (A) | 5.47499 | 5.47499 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE2_A: sd (A) | 0.019032 | 0.019032 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE_min_A: mean (A) | 4.00901 | 4.00901 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE_min_A: sd (A) | 0.01807 | 0.01807 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.00769 | 4.00769 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.011876 | 0.011876 | 0 | 0.001 | second formulas | passed |
| flip/X1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| flip/X1 | nearest water to the zinc, mean (A) | 3.647 | 3.647 | 1e-06 | 0.001 | construction | passed |
| flip/X1 | water molecules within 3.5 A of the aniline N, mean | 1.43333 | 1.43333 | 0 | 1e-06 | construction | passed |
| flip/X1 | frame 152, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| flip/X1 | frame 152, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| flip/X1 | frame 152, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| flip/X1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| flip/X1 | planned length (X0: 5 ns by the specification) (ns) | 100 | 100 | 0 | 1e-09 | construction | passed |
| flip/X1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| flip/X1 | rule 1 applies to this run (production runs of start X) | yes | yes | – | equal | construction | passed |
| flip/X1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | – | equal | second formulas | passed |
| flip/X1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| flip/X1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| flip/X1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| flip/X1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| flip/X1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| flip/X1 | gate G-V: provisional (run not finished) | yes | yes | – | equal | construction | passed |
| flip/X1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 1.8489 | 1.8489 | 0 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.932817 | -0.93282 | 3e-06 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -4.29923 | -4.29923 | 1e-06 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.771002 | 0.771004 | 2e-06 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.915209 | 0.915209 | 0 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.721267 | 0.721268 | 1e-06 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005401 | 1e-06 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043931 | 0 | 0.01 | second formulas | passed |
| flip/X1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 5 | 5 | 0 | 0 | second formulas | passed |
| flip/X1 | gate G-V: solute frames counted for that fraction | 300 | 300 | 0 | 0 | construction | passed |
| flip/X1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| flip/X1 | gate G-V: full frames counted for that fraction | 30 | 30 | 0 | 0 | construction | passed |
| flip/X1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 4.29923 | 4.29923 | 1e-06 | 0.01 | second formulas | passed |
| flip/X1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| flip/X1 | frames with ND2 nearer to the aniline N than OD1 | 155 | 155 | 0 | 0 | second formulas | passed |
| flip/X1 | state C: frames, and frames with ND2 nearer | [150, 150] | [150, 150] | – | equal | second formulas | passed |
| flip/X1 | state R: frames, and frames with ND2 nearer | [145, 0] | [145, 0] | – | equal | second formulas | passed |
| flip/X1 | state I: frames, and frames with ND2 nearer | [5, 5] | [5, 5] | – | equal | second formulas | passed |
| flip/X1 | campaign-table value: fraction of frames with ND2 nearer | 0.516667 | 0.516667 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| flip/X1 | N2: frames with positive, negative and zero signed volume | [150, 150, 0] | [150, 150, 0] | – | equal | second formulas | passed |
| flip/X1 | N2: fraction of frames with positive signed volume | 0.5 | 0.5 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2: fraction of frames with negative signed volume | 0.5 | 0.5 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2: frames of the window (last 20 ns) | 300 | 300 | 0 | 0 | construction | passed |
| flip/X1 | N2: positive, negative and zero in the last 20 ns | [150, 150, 0] | [150, 150, 0] | – | equal | second formulas | passed |
| flip/X1 | N2: fraction positive in the last 20 ns | 0.5 | 0.5 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2, frames of state C: number, positive, negative, zero | [150, 150, 0, 0] | [150, 150, 0, 0] | – | equal | second formulas | passed |
| flip/X1 | N2, frames of state C: fraction positive (none without a frame of that state) | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2, frames of state R: number, positive, negative, zero | [145, 0, 145, 0] | [145, 0, 145, 0] | – | equal | second formulas | passed |
| flip/X1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2, frames of state I: number, positive, negative, zero | [5, 0, 5, 0] | [5, 0, 5, 0] | – | equal | second formulas | passed |
| flip/X1 | N2, frames of state I: fraction positive (none without a frame of that state) | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | N2: stretches of one sign, and those of at least 1 ns | [2, 2] | [2, 2] | – | equal | second formulas | passed |
| flip/X1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[-1, 0, 150], [1, 150, 150]] | [[-1, 0, 150], [1, 150, 150]] | – | equal | second formulas | passed |
| flip/X1 | N2: number of sign changes that last at least 1 ns | 2 | 2 | 0 | 0 | second formulas | passed |
| flip/X1 | N2: these sign changes (frame, from, to) | [[0, 1, -1], [150, -1, 1]] | [[0, 1, -1], [150, -1, 1]] | – | equal | second formulas | passed |
| flip/X1 | N2: times of these sign changes (ns) | [0.01, 1.51] | [0.01, 1.51] | – | equal | second formulas | passed |
| flip/X1 | N2: time of the first lasting sign change (ns) | 0.01 | 0.01 | 0 | 1e-09 | second formulas | passed |
| flip/X1 | N2: sign of the first and of the last frame | [-1, 1] | [-1, 1] | – | equal | second formulas | passed |
| flip/X1 | hydrogen bond around the first passage: reported (a first passage exists) | yes | yes | – | equal | second formulas | passed |
| flip/X1 | hydrogen bond around the first passage: frame of the first passage | 5 | 5 | – | equal | second formulas | passed |
| flip/X1 | hydrogen bond N2-H...OG1 before the first passage: frames, frames with the bond | [5, 0] | [5, 0] | – | equal | second formulas | passed |
| flip/X1 | hydrogen bond N2-H...OG1 from the first passage on: frames, frames with the bond | [295, 245] | [295, 245] | – | equal | second formulas | passed |
| flip/X1 | hydrogen bond N2-H...OG1, fraction before the first passage | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | hydrogen bond N2-H...OG1, fraction from the first passage on | 0.830508 | 0.830508 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [0.0, 0.830508] | [0.0, 0.830508] | – | equal | construction | passed |
| flip/X1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| flip/X1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| flip/X1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| flip/X1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.00901 | 4.00901 | 0 | 0.001 | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 4e-06 | 0.001 | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 300] | [0, 300] | – | equal | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47499 | 5.47499 | 0 | 0.001 | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| flip/X1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 300] | [0, 300] | – | equal | second formulas | passed |
| flip/X1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.009015, 4.003723, 0.0] | [4.009015, 4.003723, 0.0] | – | equal | second formulas | passed |
| flip/X1 | water at the crystal ring positions: mean number within 3.5 A | 0.733333 | 0.733333 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 6e-06 | 0.001 | second formulas | passed |
| flip/X1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.38 | 3.38 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 0.733333 | 0.733333 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| flip/X1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| flip/X1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| flip/X1 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.6 | 3.59999 | 1.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.599987] | [0.0, 3.599987] | – | equal | second formulas | passed |
| flip/X1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| flip/X1 | first production frame: state, and whether it differs from the state of the start | ['I', True] | ['I', True] | – | equal | second formulas | passed |
| flip/X1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| flip/X1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, True] | [-1, True] | – | equal | second formulas | passed |
| flip/X1 | first production frame: signed volume at N2 | -0.556792 | -0.556792 | 0 | 0.001 | second formulas | passed |
| flip/X1 | first production frame: ligand RMSD (A) | 3 | 3 | 0 | 0.001 | second formulas | passed |
| flip/X1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['I', True, -1] | ['I', True, -1] | – | equal | second formulas | passed |
| flip/X1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.5, 0.5, 2] | [0.5, 0.5, 2] | – | equal | second formulas | passed |
| flip/X1 | equilibration: frames read (files solute_eqNN.dcd present) | yes | yes | – | equal | construction | passed |
| flip/X1 | equilibration: time series written (only with equilibration files) | [True, True] | [True, True] | – | equal | construction | passed |
| flip/X1 | campaign-table value: solute frames of the equilibration | 200 | 200 | 0 | 0 | construction | passed |
| flip/X1 | equilibration: solute frames and full frames read | [200, 20] | [200, 20] | – | equal | construction | passed |
| flip/X1 | equilibration: length of the equilibration (ps) | 2000 | 2000 | 0 | 1e-09 | construction | passed |
| flip/X1 | equilibration: where the stage table came from | protocol (constants of the specification) | protocol (constants of the specification) | – | equal | construction | passed |
| flip/X1 | equilibration: lines of timeseries_equil_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| flip/X1 | equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned) | yes | yes | – | equal | construction | passed |
| flip/X1 | equilibration, every frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | equilibration, every frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| flip/X1 | equilibration: no time above zero | yes | yes | – | equal | construction | passed |
| flip/X1 | equilibration, every frame: lig_rmsd_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 3.9e-05 | 3.9e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000128 | 0.000128 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.4e-05 | 4.4e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000588 | 0.000588 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000209 | 0.000209 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000517 | 0.000517 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000758 | 0.000758 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 4.3e-05 | 4.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000739 | 0.000739 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000931 | 0.000931 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 3.8e-05 | 3.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration: lines of timeseries_equil_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| flip/X1 | equilibration, every full frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| flip/X1 | equilibration, every full frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| flip/X1 | equilibration, every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| flip/X1 | equilibration, every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration: the steps with a file, in the order of the run | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 1: start and end (ps) | [0.0, 100.0] | [0.0, 100.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 1: solute and full frames read, and expected | [10, 10, 1, 1] | [10, 10, 1, 1] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 1: frames of state C, R, I; frames with positive and negative sign at N2 | [10, 0, 0, 10, 0] | [10, 0, 0, 10, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 1: last solute frame (frame, time, state, sign at N2) | [9, -1.9, 'C', 1] | [9, -1.9, 'C', 1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 1: last solute frame, ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 1: last solute frame, signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [0, -1.9, 0] | [0, -1.9, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 1: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.59999 | 8e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.59999 | 8e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 2: start and end (ps) | [100.0, 400.0] | [100.0, 400.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 2: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 2: frames of state C, R, I; frames with positive and negative sign at N2 | [30, 0, 0, 30, 0] | [30, 0, 0, 30, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 2: last solute frame (frame, time, state, sign at N2) | [39, -1.6, 'C', 1] | [39, -1.6, 'C', 1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 2: last solute frame, ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 2: last solute frame, signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [3, -1.6, 3] | [3, -1.6, 3] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 2: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 2: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.59999 | 7e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3 | 2.99999 | 6e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 3: start and end (ps) | [400.0, 700.0] | [400.0, 700.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 3: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 3: frames of state C, R, I; frames with positive and negative sign at N2 | [30, 0, 0, 30, 0] | [30, 0, 0, 30, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 3: last solute frame (frame, time, state, sign at N2) | [69, -1.3, 'C', 1] | [69, -1.3, 'C', 1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 3: last solute frame, ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 3: last solute frame, signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 3: last full frame (frame, time, water oxygens at the crystal ring positions) | [6, -1.3, 2] | [6, -1.3, 2] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 3: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 3: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 2e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 3: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3 | 3 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 4: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 4: start and end (ps) | [700.0, 1100.0] | [700.0, 1100.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 4: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 4: frames of state C, R, I; frames with positive and negative sign at N2 | [40, 0, 0, 40, 0] | [40, 0, 0, 40, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 4: last solute frame (frame, time, state, sign at N2) | [109, -0.9, 'C', 1] | [109, -0.9, 'C', 1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 4: last solute frame, ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 4: last solute frame, signed volume at N2 | 0.723989 | 0.72399 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 4: last full frame (frame, time, water oxygens at the crystal ring positions) | [10, -0.9, 2] | [10, -0.9, 2] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 4: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 4: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 4: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3 | 3 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 5: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 5: start and end (ps) | [1100.0, 1500.0] | [1100.0, 1500.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 5: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | [10, 30, 0, 10, 30] | [10, 30, 0, 10, 30] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 5: last solute frame (frame, time, state, sign at N2) | [149, -0.5, 'R', -1] | [149, -0.5, 'R', -1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 5: last solute frame, ligand RMSD (A) | 3.56157 | 3.56157 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 5: last solute frame, signed volume at N2 | -0.47521 | -0.47521 | 0 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 5: last full frame (frame, time, water oxygens at the crystal ring positions) | [14, -0.5, 0] | [14, -0.5, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 5: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 5: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 5: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.6 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 6: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 6: start and end (ps) | [1500.0, 2000.0] | [1500.0, 2000.0] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 6: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| flip/X1 | equilibration, stage 6: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 6: last solute frame (frame, time, state, sign at N2) | [199, 0.0, 'R', -1] | [199, 0.0, 'R', -1] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 6: last solute frame, ligand RMSD (A) | 3.56157 | 3.56156 | 3e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 6: last solute frame, signed volume at N2 | -0.47521 | -0.475209 | 1e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 6: last full frame (frame, time, water oxygens at the crystal ring positions) | [19, 0.0, 0] | [19, 0.0, 0] | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 6: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| flip/X1 | equilibration, stage 6: last full frame, nearest water oxygen to the place of the ring (A) | 3.6 | 3.59999 | 6e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration, stage 6: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 3.6 | 3.59999 | 6e-06 | 0.001 | second formulas | passed |
| flip/X1 | equilibration: frames of state C, R, I; positive and negative sign at N2 | [120, 80, 0, 120, 80] | [120, 80, 0, 120, 80] | – | equal | second formulas | passed |
| flip/X1 | equilibration: first_frame_not_in_state_of_start (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | second formulas | passed |
| flip/X1 | equilibration: first_frame_in_opposite_state (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | second formulas | passed |
| flip/X1 | equilibration: first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | second formulas | passed |
| flip/X1 | equilibration: state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | second formulas | passed |
| flip/X1 | equilibration: time of the last solute frame (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| flip/X1 | equilibration: frames in which the ligand and the zinc were in another periodic image | [29, 18] | [29, 18] | – | equal | construction | passed |
| flip/X1 | by design: frames with positive and with negative signed volume at N2 | [150, 150] | [150, 150] | – | equal | construction | passed |
| flip/X1 | by design: positive and negative in the last 20 ns | [150, 150] | [150, 150] | – | equal | construction | passed |
| flip/X1 | by design, frames of state C: number, positive, negative | [150, 150, 0] | [150, 150, 0] | – | equal | construction | passed |
| flip/X1 | by design, frames of state R: number, positive, negative | [145, 0, 145] | [145, 0, 145] | – | equal | construction | passed |
| flip/X1 | by design, frames of state I: number, positive, negative | [5, 0, 5] | [5, 0, 5] | – | equal | construction | passed |
| flip/X1 | by design: stretches of one sign, and those of at least 1 ns | [2, 2] | [2, 2] | – | equal | construction | passed |
| flip/X1 | by design: times of the sign changes that last at least 1 ns (ns) | [0.01, 1.51] | [0.01, 1.51] | – | equal | construction | passed |
| flip/X1 | by design: sign after each of these changes | [-1, 1] | [-1, 1] | – | equal | construction | passed |
| flip/X1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [5, 0, 295, 245] | [5, 0, 295, 245] | – | equal | construction | passed |
| flip/X1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| flip/X1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['I', True, -1, True] | ['I', True, -1, True] | – | equal | construction | passed |
| flip/X1 | by design: equilibration frames read | yes | yes | – | equal | construction | passed |
| flip/X1 | by design: equilibration, solute frames, full frames, length (ps) | [200, 20, 2000.0] | [200, 20, 2000.0] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, where the stage table came from | protocol (constants of the specification) | protocol (constants of the specification) | – | equal | construction | passed |
| flip/X1 | by design: equilibration, solute and full frames of each stage | [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]] | [[10, 1], [30, 3], [30, 3], [40, 4], [40, 4], [50, 5]] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, time of the first and of the last solute frame (ns) | [-1.99, 0.0] | [-1.99, 0.0] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, time of solute frame 10 (the first of stage 2) (ns) | [-1.89, '2'] | [-1.89, '2'] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, first_frame_in_opposite_state (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | [120, 5, 0, -0.79] | [120, 5, 0, -0.79] | – | equal | construction | passed |
| flip/X1 | by design: equilibration, state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | construction | passed |
| glu/X0 | solute frames seen | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| glu/X0 | manifest lists the segments | no | no | – | equal | construction | passed |
| glu/X0 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| glu/X0 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| glu/X0 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| glu/X0 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| glu/X0 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| glu/X0 | full frames seen | 50 | 50 | 0 | 0 | construction | passed |
| glu/X0 | time seen (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| glu/X0 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| glu/X0 | frames in which the ligand was in another periodic image | 71 | 71 | 0 | 0 | construction | passed |
| glu/X0 | frames in which the zinc was in another periodic image | 45 | 45 | 0 | 0 | construction | passed |
| glu/X0 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | lines of timeseries_solute.csv | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| glu/X0 | every frame: lig_rmsd_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.3e-05 | 1.3e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001074 | 0.001074 | 0.01 | second formulas | passed |
| glu/X0 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.8e-05 | 5.8e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 9e-06 | 9e-06 | 0.001 | second formulas | passed |
| glu/X0 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 3.7e-05 | 3.7e-05 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000638 | 0.000638 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000153 | 0.000153 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000896 | 0.000896 | 0.01 | second formulas | passed |
| glu/X0 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000261 | 0.000261 | 0.01 | second formulas | passed |
| glu/X0 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 6e-05 | 6e-05 | 0.01 | second formulas | passed |
| glu/X0 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| glu/X0 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000556 | 0.000556 | 0.01 | second formulas | passed |
| glu/X0 | every frame: tau_fold_deg, largest difference (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| glu/X0 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000598 | 0.000598 | 0.01 | second formulas | passed |
| glu/X0 | every frame: d_S_N2_A, largest difference (A) | 0 | 3.3e-05 | 3.3e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | every frame: n2_volume, largest difference | 0 | 1.1e-05 | 1.1e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| glu/X0 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | lines of timeseries_full.csv | 50 | 50 | 0 | 0 | construction | passed |
| glu/X0 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| glu/X0 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| glu/X0 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| glu/X0 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| glu/X0 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| glu/X0 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| glu/X0 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | verdict | STAYED | STAYED | – | equal | construction | passed |
| glu/X0 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| glu/X0 | run complete | yes | yes | – | equal | construction | passed |
| glu/X0 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| glu/X0 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| glu/X0 | frames of the verdict window | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | count_C_window | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | count_C | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| glu/X0 | count_C (design against the planted frames) | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | frames with the hydrogen bond N2-H...OG1 | 400 | 400 | 0 | 0 | second formulas | passed |
| glu/X0 | full frames with a water within 2.6 A of the zinc | 2 | 2 | 0 | 0 | construction | passed |
| glu/X0 | frames with labelling identity | 125 | 125 | 0 | 0 | second formulas | passed |
| glu/X0 | frames with labelling ring_flip | 125 | 125 | 0 | 0 | second formulas | passed |
| glu/X0 | frames with labelling oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| glu/X0 | frames with labelling ring_flip+oxygen_swap | 125 | 125 | 0 | 0 | second formulas | passed |
| glu/X0 | lig_rmsd_A: mean (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | lig_rmsd_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | lig_rmsd_A: mean_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | lig_rmsd_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_N2_A: mean (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_N2_A: mean_last_half (A) | 1.96144 | 1.96144 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | tau_fold_deg: mean (deg) | 132.056 | 132.056 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | tau_fold_deg: sd (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | second formulas | passed |
| glu/X0 | tau_fold_deg: mean_last_half (deg) | 132.056 | 132.056 | 9e-06 | 0.01 | second formulas | passed |
| glu/X0 | tau_fold_deg: sd_last_half (deg) | 0 | 7.3e-05 | 7.3e-05 | 0.01 | second formulas | passed |
| glu/X0 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| glu/X0 | ang_C4_S_N2_deg: mean (deg) | 108.795 | 108.795 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | ang_C4_S_N2_deg: sd (deg) | 0 | 6.2e-05 | 6.2e-05 | 0.01 | second formulas | passed |
| glu/X0 | ang_C4_S_N2_deg: mean_last_half (deg) | 108.795 | 108.795 | 6e-06 | 0.01 | second formulas | passed |
| glu/X0 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 6.1e-05 | 6.1e-05 | 0.01 | second formulas | passed |
| glu/X0 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| glu/X0 | d_S_N2_A: mean (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | d_S_N2_A: mean_last_half (A) | 1.59633 | 1.59633 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| glu/X0 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| glu/X0 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d1_N1_F131CZ_A: mean (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | d1_N1_F131CZ_A: mean_last_half (A) | 3.8966 | 3.8966 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d2_N1_N67OD1_A: mean (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | d2_N1_N67OD1_A: mean_last_half (A) | 7.74431 | 7.74431 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| glu/X0 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | n2_volume: mean | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| glu/X0 | n2_volume: sd | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | n2_volume: mean_last_half | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| glu/X0 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE1_A: mean (A) | 3.91746 | 3.91746 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE1_A: sd (A) | 0.358311 | 0.358311 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE1_A: mean_last_half (A) | 3.97761 | 3.97761 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE1_A: sd_last_half (A) | 0.21142 | 0.21142 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE2_A: mean (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE2_A: sd (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE2_A: sd_last_half (A) | 0.011929 | 0.011929 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE_min_A: mean (A) | 3.91746 | 3.91746 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE_min_A: sd (A) | 0.358311 | 0.358311 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE_min_A: mean_last_half (A) | 3.97761 | 3.97761 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.21142 | 0.21142 | 0 | 0.001 | second formulas | passed |
| glu/X0 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| glu/X0 | nearest water to the zinc, mean (A) | 3.6384 | 3.6384 | 1e-06 | 0.001 | construction | passed |
| glu/X0 | water molecules within 3.5 A of the aniline N, mean | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| glu/X0 | frame 8, lig_rmsd_A against starts.csv (X, rmsd_to_crystal_A) (A) | 0 | 0 | 0 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, d_Zn_N2_A against starts.csv (X, Zn_N2) (A) | 1.931 | 1.9314 | 0.0004 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, d_Zn_Osulfonyl_min_A against starts.csv (X, Zn_O_min) (A) | 2.986 | 2.9856 | 0.0004 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, d_N2_T199OG1_A against starts.csv (X, N2_OG1) (A) | 2.861 | 2.861 | 0 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, ang_N2_HN2_T199OG1_deg against starts.csv (X, N2_H_OG1) (deg) | 163.678 | 163.69 | 0.012 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, d_T199N_Osulfonyl_min_A against starts.csv (X, T199N_O_min) (A) | 2.941 | 2.9406 | 0.0004 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, tors_C3_C4_S_N2_deg against starts.csv (X, C3_C4_S_N2) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, tors_C5_C4_S_N2_deg against starts.csv (X, C5_C4_S_N2) (deg) | -50.092 | -50.092 | 0 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, tau_fold_deg against starts.csv (X, tau_fold) (deg) | 132.056 | 132.056 | 0 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, ang_C4_S_N2_deg against starts.csv (X, C4_S_N2) (deg) | 108.795 | 108.795 | 0 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, d_S_N2_A against starts.csv (X, S_N2) (A) | 1.596 | 1.5963 | 0.0003 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, ang_S_N2_Zn_deg against starts.csv (X, S_N2_Zn) (deg) | 117.552 | 117.552 | 0 | 0.06 | starts.csv | passed |
| glu/X0 | frame 8, d1_N1_F131CZ_A against starts.csv (X, N1_Phe131CZ) (A) | 3.897 | 3.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, d2_N1_N67OD1_A against starts.csv (X, N1_Asn67OD1) (A) | 7.744 | 7.7443 | 0.0003 | 0.0015 | starts.csv | passed |
| glu/X0 | frame 8, d_N1_N67ND2_A against starts.csv (X, N1_Asn67ND2) (A) | 6.551 | 6.5511 | 0.0001 | 0.0015 | starts.csv | passed |
| glu/X0 | role of the run (X0 = validation) | validation | validation | – | equal | construction | passed |
| glu/X0 | planned length (X0: 5 ns by the specification) (ns) | 5 | 5 | 0 | 1e-09 | construction | passed |
| glu/X0 | gate G-V applies to this run (only to X0) | yes | yes | – | equal | construction | passed |
| glu/X0 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| glu/X0 | gate G-V: the criteria that are not met | none | none | – | equal | second formulas | passed |
| glu/X0 | gate G-V: all criteria met | yes | yes | – | equal | second formulas | passed |
| glu/X0 | gate G-V: GV3 as a whole | yes | yes | – | equal | second formulas | passed |
| glu/X0 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| glu/X0 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| glu/X0 | gate G-V: verdict (none for a run that is not X0) | PASSED | PASSED | – | equal | second formulas | passed |
| glu/X0 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| glu/X0 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | -0.823964 | -0.823962 | 2e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -0.840545 | -0.840543 | 2e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -0.776334 | -0.776335 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.731087 | 0.731086 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.807205 | 0.807203 | 2e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | -0.005402 | -0.005401 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: tau_fold_deg, mean minus the reference (deg) | -0.043931 | -0.043932 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | gate G-V: solute frames counted for that fraction | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 2 | 2 | 0 | 0 | second formulas | passed |
| glu/X0 | gate G-V: full frames counted for that fraction | 50 | 50 | 0 | 0 | construction | passed |
| glu/X0 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 0.930667 | 0.930666 | 1e-06 | 0.01 | second formulas | passed |
| glu/X0 | gate G-V by design: verdict | PASSED | PASSED | – | equal | construction | passed |
| glu/X0 | gate G-V by design: the criteria that are not met | none | none | – | equal | construction | passed |
| glu/X0 | gate G-V by design: provisional | no | no | – | equal | construction | passed |
| glu/X0 | by design: solute frames with a sulfonyl oxygen within 2.6 A | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | by design: full frames with any oxygen within 2.6 A | 2 | 2 | 0 | 0 | construction | passed |
| glu/X0 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| glu/X0 | frames with ND2 nearer to the aniline N than OD1 | 500 | 500 | 0 | 0 | second formulas | passed |
| glu/X0 | state C: frames, and frames with ND2 nearer | [500, 500] | [500, 500] | – | equal | second formulas | passed |
| glu/X0 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| glu/X0 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| glu/X0 | by design: ND2 nearer in state C, frames of state C, ND2 nearer in state R, frames of state R | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | construction | passed |
| glu/X0 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| glu/X0 | N2: frames with positive, negative and zero signed volume | [500, 0, 0] | [500, 0, 0] | – | equal | second formulas | passed |
| glu/X0 | N2: fraction of frames with positive signed volume | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | N2: fraction of frames with negative signed volume | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | N2: frames of the window (last 20 ns) | 500 | 500 | 0 | 0 | construction | passed |
| glu/X0 | N2: positive, negative and zero in the last 20 ns | [500, 0, 0] | [500, 0, 0] | – | equal | second formulas | passed |
| glu/X0 | N2: fraction positive in the last 20 ns | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | N2, frames of state C: number, positive, negative, zero | [500, 500, 0, 0] | [500, 500, 0, 0] | – | equal | second formulas | passed |
| glu/X0 | N2, frames of state C: fraction positive (none without a frame of that state) | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| glu/X0 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| glu/X0 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| glu/X0 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| glu/X0 | N2: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| glu/X0 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[1, 0, 500]] | [[1, 0, 500]] | – | equal | second formulas | passed |
| glu/X0 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| glu/X0 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| glu/X0 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| glu/X0 | N2: sign of the first and of the last frame | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| glu/X0 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| glu/X0 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| glu/X0 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| glu/X0 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| glu/X0 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [30, 0, 30] | [30, 0, 30] | – | equal | second formulas | passed |
| glu/X0 | Glu106: smallest distance zinc to OE1 (A) | 2.5 | 2.5 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| glu/X0 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE1: mean distance (A) | 3.91746 | 3.91746 | 0 | 0.001 | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE1: smallest distance (A) | 2.5 | 2.5 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [30, 500] | [30, 500] | – | equal | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47407 | 5.47407 | 0 | 0.001 | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| glu/X0 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 500] | [0, 500] | – | equal | second formulas | passed |
| glu/X0 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [3.917462, 2.499999, 0.06] | [3.917462, 2.499999, 0.06] | – | equal | second formulas | passed |
| glu/X0 | water at the crystal ring positions: mean number within 3.5 A | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 3 | 2.99999 | 1.3e-05 | 0.001 | second formulas | passed |
| glu/X0 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 3.156 | 3.156 | 1e-06 | 0.001 | second formulas | passed |
| glu/X0 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| glu/X0 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| glu/X0 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| glu/X0 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.6 | 3.59999 | 9e-06 | 0.001 | second formulas | passed |
| glu/X0 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 3.6 | 3.6 | 0 | 0.001 | second formulas | passed |
| glu/X0 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.599991] | [0.0, 3.599991] | – | equal | second formulas | passed |
| glu/X0 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.99999 | 2.99999 | – | equal | second formulas | passed |
| glu/X0 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| glu/X0 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| glu/X0 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| glu/X0 | first production frame: signed volume at N2 | 0.723989 | 0.723989 | 0 | 0.001 | second formulas | passed |
| glu/X0 | first production frame: ligand RMSD (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| glu/X0 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| glu/X0 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [1.0, 1.0, 0] | [1.0, 1.0, 0] | – | equal | second formulas | passed |
| glu/X0 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| glu/X0 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| glu/X0 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | by design: frames with positive and with negative signed volume at N2 | [500, 0] | [500, 0] | – | equal | construction | passed |
| glu/X0 | by design: positive and negative in the last 20 ns | [500, 0] | [500, 0] | – | equal | construction | passed |
| glu/X0 | by design, frames of state C: number, positive, negative | [500, 500, 0] | [500, 500, 0] | – | equal | construction | passed |
| glu/X0 | by design, frames of state R: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| glu/X0 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| glu/X0 | by design: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | construction | passed |
| glu/X0 | by design: times of the sign changes that last at least 1 ns (ns) | [] | [] | – | equal | construction | passed |
| glu/X0 | by design: sign after each of these changes | [] | [] | – | equal | construction | passed |
| glu/X0 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [None, None, None, None] | [None, None, None, None] | – | equal | construction | passed |
| glu/X0 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [30, 0, 30] | [30, 0, 30] | – | equal | construction | passed |
| glu/X0 | by design: smallest distance zinc to OE1 of Glu106 (A) | 2.5 | 2.5 | 1e-06 | 0.001 | construction | passed |
| glu/X0 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| glu/X0 | by design: mean number of water oxygens at the crystal ring positions | 1.46 | 1.46 | 0 | 1e-06 | construction | passed |
| glu/X0 | by design: smallest distance of a water oxygen to the crystal ring positions (A) | 3 | 2.99999 | 1.3e-05 | 0.001 | construction | passed |
| glu/X0 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['C', False, 1, False] | ['C', False, 1, False] | – | equal | construction | passed |
| glu/X0 | by design: equilibration frames read | no | no | – | equal | construction | passed |
| scaled/N1 | solute frames seen | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| scaled/N1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| scaled/N1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| scaled/N1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| scaled/N1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| scaled/N1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| scaled/N1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| scaled/N1 | full frames seen | 20 | 20 | 0 | 0 | construction | passed |
| scaled/N1 | time seen (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| scaled/N1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| scaled/N1 | frames in which the ligand was in another periodic image | 29 | 29 | 0 | 0 | construction | passed |
| scaled/N1 | frames in which the zinc was in another periodic image | 18 | 18 | 0 | 0 | construction | passed |
| scaled/N1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | lines of timeseries_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| scaled/N1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 8e-06 | 8e-06 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000821 | 0.000821 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 2.1e-05 | 2.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000293 | 0.000293 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000383 | 0.000383 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000682 | 0.000682 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000691 | 0.000691 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000526 | 0.000526 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 5.6e-05 | 5.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 4.8e-05 | 4.8e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 4.8e-05 | 4.8e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 1.7e-05 | 1.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.00066 | 0.00066 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000698 | 0.000698 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000921 | 0.000921 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000698 | 0.000698 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000512 | 0.000512 | 0.01 | second formulas | passed |
| scaled/N1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: n2_volume, largest difference | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | lines of timeseries_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| scaled/N1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| scaled/N1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| scaled/N1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.4e-05 | 2.4e-05 | 0.001 | second formulas | passed |
| scaled/N1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 1.6e-05 | 1.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| scaled/N1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | verdict | STAYED | STAYED | – | equal | construction | passed |
| scaled/N1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| scaled/N1 | run complete | yes | yes | – | equal | construction | passed |
| scaled/N1 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| scaled/N1 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| scaled/N1 | frames of the verdict window | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | count_C_window | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | count_R_window | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | count_C | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | count_R | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | rule 1: left the crystal pose | no | no | – | equal | construction | passed |
| scaled/N1 | count_C (design against the planted frames) | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | count_R (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | frames with the hydrogen bond N2-H...OG1 | 200 | 200 | 0 | 0 | second formulas | passed |
| scaled/N1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| scaled/N1 | frames with labelling identity | 200 | 200 | 0 | 0 | second formulas | passed |
| scaled/N1 | frames with labelling ring_flip | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | frames with labelling oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | frames with labelling ring_flip+oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | lig_rmsd_A: mean (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | lig_rmsd_A: mean_last_half (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_N2_A: mean (A) | 2.23712 | 2.23712 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_N2_A: mean_last_half (A) | 2.23712 | 2.23712 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | tau_fold_deg: mean (deg) | 148.632 | 148.632 | 3e-06 | 0.01 | second formulas | passed |
| scaled/N1 | tau_fold_deg: sd (deg) | 0 | 9.5e-05 | 9.5e-05 | 0.01 | second formulas | passed |
| scaled/N1 | tau_fold_deg: mean_last_half (deg) | 148.632 | 148.632 | 5e-06 | 0.01 | second formulas | passed |
| scaled/N1 | tau_fold_deg: sd_last_half (deg) | 0 | 8.4e-05 | 8.4e-05 | 0.01 | second formulas | passed |
| scaled/N1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| scaled/N1 | ang_C4_S_N2_deg: mean (deg) | 109.483 | 109.483 | 1e-06 | 0.01 | second formulas | passed |
| scaled/N1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5e-05 | 5e-05 | 0.01 | second formulas | passed |
| scaled/N1 | ang_C4_S_N2_deg: mean_last_half (deg) | 109.483 | 109.483 | 5e-06 | 0.01 | second formulas | passed |
| scaled/N1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 3.1e-05 | 3.1e-05 | 0.01 | second formulas | passed |
| scaled/N1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| scaled/N1 | d_S_N2_A: mean (A) | 1.67745 | 1.67745 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d_S_N2_A: mean_last_half (A) | 1.67745 | 1.67745 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d1_N1_F131CZ_A: mean (A) | 3.35615 | 3.35615 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d1_N1_F131CZ_A: mean_last_half (A) | 3.35615 | 3.35615 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d2_N1_N67OD1_A: mean (A) | 7.72547 | 7.72547 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d2_N1_N67OD1_A: mean_last_half (A) | 7.72547 | 7.72547 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | n2_volume: mean | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | n2_volume: sd | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | n2_volume: mean_last_half | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE1_A: mean (A) | 4.01151 | 4.01151 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE1_A: sd (A) | 0.023337 | 0.023337 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01151 | 4.0115 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.023337 | 0.023337 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE2_A: mean (A) | 5.47696 | 5.47696 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE2_A: sd (A) | 0.020621 | 0.020621 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47696 | 5.47696 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.020621 | 0.020621 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE_min_A: mean (A) | 4.01151 | 4.01151 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE_min_A: sd (A) | 0.023337 | 0.023337 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01151 | 4.0115 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.023337 | 0.023337 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| scaled/N1 | nearest water to the zinc, mean (A) | 3.6255 | 3.6255 | 0 | 0.001 | construction | passed |
| scaled/N1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| scaled/N1 | frame 1, lig_rmsd_A against starts.csv (N, rmsd_to_crystal_A) (A) | 0.61 | 0.6104 | 0.0004 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, d_Zn_N2_A against starts.csv (N, Zn_N2) (A) | 2.207 | 2.2071 | 0.0001 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, d_Zn_Osulfonyl_min_A against starts.csv (N, Zn_O_min) (A) | 2.775 | 2.7751 | 0.0001 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, d_N2_T199OG1_A against starts.csv (N, N2_OG1) (A) | 2.552 | 2.5522 | 0.0002 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, ang_N2_HN2_T199OG1_deg against starts.csv (N, N2_H_OG1) (deg) | 149.15 | 149.113 | 0.037 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, d_T199N_Osulfonyl_min_A against starts.csv (N, T199N_O_min) (A) | 2.983 | 2.9835 | 0.0005 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, tors_C3_C4_S_N2_deg against starts.csv (N, C3_C4_S_N2) (deg) | 148.632 | 148.632 | 0 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, tors_C5_C4_S_N2_deg against starts.csv (N, C5_C4_S_N2) (deg) | -28.91 | -28.91 | 0 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, tau_fold_deg against starts.csv (N, tau_fold) (deg) | 148.632 | 148.632 | 0 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, ang_C4_S_N2_deg against starts.csv (N, C4_S_N2) (deg) | 109.483 | 109.483 | 0 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, d_S_N2_A against starts.csv (N, S_N2) (A) | 1.677 | 1.6775 | 0.0005 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, ang_S_N2_Zn_deg against starts.csv (N, S_N2_Zn) (deg) | 105.908 | 105.908 | 0 | 0.06 | starts.csv | passed |
| scaled/N1 | frame 1, d1_N1_F131CZ_A against starts.csv (N, N1_Phe131CZ) (A) | 3.356 | 3.3561 | 0.0001 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, d2_N1_N67OD1_A against starts.csv (N, N1_Asn67OD1) (A) | 7.725 | 7.7255 | 0.0005 | 0.0015 | starts.csv | passed |
| scaled/N1 | frame 1, d_N1_N67ND2_A against starts.csv (N, N1_Asn67ND2) (A) | 6.557 | 6.5574 | 0.0004 | 0.0015 | starts.csv | passed |
| scaled/N1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| scaled/N1 | planned length (X0: 5 ns by the specification) (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| scaled/N1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| scaled/N1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| scaled/N1 | gate G-V: the criteria that are not met | GV3 Zn-N2 | GV3 Zn-N2 | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| scaled/N1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| scaled/N1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 7.57048 | 7.57048 | 3e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -3.51878 | -3.51878 | 1e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.34454 | -7.34454 | 2e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.857168 | 0.857168 | 0 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.98182 | 0.981818 | 2e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.62485 | 0.624848 | 2e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.682512 | 0.682513 | 1e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 16.5317 | 16.5317 | 3e-06 | 0.01 | second formulas | passed |
| scaled/N1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | gate G-V: solute frames counted for that fraction | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| scaled/N1 | gate G-V: full frames counted for that fraction | 20 | 20 | 0 | 0 | construction | passed |
| scaled/N1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 7.57048 | 7.57048 | 3e-06 | 0.01 | second formulas | passed |
| scaled/N1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| scaled/N1 | frames with ND2 nearer to the aniline N than OD1 | 200 | 200 | 0 | 0 | second formulas | passed |
| scaled/N1 | state C: frames, and frames with ND2 nearer | [200, 200] | [200, 200] | – | equal | second formulas | passed |
| scaled/N1 | state R: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| scaled/N1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| scaled/N1 | campaign-table value: fraction of frames with ND2 nearer | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | N2: sign of the start as built (X and N +1, R -1) | 1 | 1 | 0 | 0 | construction | passed |
| scaled/N1 | N2: frames with positive, negative and zero signed volume | [200, 0, 0] | [200, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | N2: fraction of frames with positive signed volume | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | N2: fraction of frames with negative signed volume | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | N2: frames of the window (last 20 ns) | 200 | 200 | 0 | 0 | construction | passed |
| scaled/N1 | N2: positive, negative and zero in the last 20 ns | [200, 0, 0] | [200, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | N2: fraction positive in the last 20 ns | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | N2, frames of state C: number, positive, negative, zero | [200, 200, 0, 0] | [200, 200, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | N2, frames of state C: fraction positive (none without a frame of that state) | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | N2, frames of state R: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | N2, frames of state R: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| scaled/N1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| scaled/N1 | N2: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| scaled/N1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[1, 0, 200]] | [[1, 0, 200]] | – | equal | second formulas | passed |
| scaled/N1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| scaled/N1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| scaled/N1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| scaled/N1 | N2: sign of the first and of the last frame | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| scaled/N1 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| scaled/N1 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| scaled/N1 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| scaled/N1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01151 | 4.01151 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47696 | 5.47696 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 4e-06 | 0.001 | second formulas | passed |
| scaled/N1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| scaled/N1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.011506, 4.003724, 0.0] | [4.011506, 4.003724, 0.0] | – | equal | second formulas | passed |
| scaled/N1 | water at the crystal ring positions: mean number within 3.5 A | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 2.32712 | 2.32711 | 9e-06 | 0.001 | second formulas | passed |
| scaled/N1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 2.81475 | 2.81475 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.5 | 1.5 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| scaled/N1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| scaled/N1 | water in the place of the ring: mean number | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | water in the place of the ring: largest number in a full frame | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| scaled/N1 | water in the place of the ring: smallest distance of a water oxygen (A) | 3.48668 | 3.48667 | 1.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 3.54334 | 3.54334 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [0.0, 3.486673] | [0.0, 3.486673] | – | equal | second formulas | passed |
| scaled/N1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 2.32711 | 2.32711 | – | equal | second formulas | passed |
| scaled/N1 | first production frame: state, and whether it differs from the state of the start | ['C', False] | ['C', False] | – | equal | second formulas | passed |
| scaled/N1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| scaled/N1 | first production frame: sign at N2, and whether it differs from the sign of the start | [1, False] | [1, False] | – | equal | second formulas | passed |
| scaled/N1 | first production frame: signed volume at N2 | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | first production frame: ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['C', False, 1] | ['C', False, 1] | – | equal | second formulas | passed |
| scaled/N1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [1.0, 1.0, 0] | [1.0, 1.0, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration: frames read (files solute_eqNN.dcd present) | yes | yes | – | equal | construction | passed |
| scaled/N1 | equilibration: time series written (only with equilibration files) | [True, True] | [True, True] | – | equal | construction | passed |
| scaled/N1 | campaign-table value: solute frames of the equilibration | 8 | 8 | 0 | 0 | construction | passed |
| scaled/N1 | equilibration: solute frames and full frames read | [8, 0] | [8, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration: length of the equilibration (ps) | 100 | 100 | 0 | 1e-09 | construction | passed |
| scaled/N1 | equilibration: where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| scaled/N1 | equilibration: lines of timeseries_equil_solute.csv | 8 | 8 | 0 | 0 | construction | passed |
| scaled/N1 | equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned) | yes | yes | – | equal | construction | passed |
| scaled/N1 | equilibration, every frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | equilibration, every frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| scaled/N1 | equilibration: no time above zero | yes | yes | – | equal | construction | passed |
| scaled/N1 | equilibration, every frame: lig_rmsd_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 4.5e-05 | 4.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 2.9e-05 | 2.9e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 8e-06 | 8e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000179 | 0.000179 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_N2_A, largest difference (A) | 0 | 2.3e-05 | 2.3e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 4.7e-05 | 4.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 2.1e-05 | 2.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000293 | 0.000293 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000383 | 0.000383 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000481 | 0.000481 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000355 | 0.000355 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000526 | 0.000526 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_O1_A, largest difference (A) | 0 | 4.4e-05 | 4.4e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_O2_A, largest difference (A) | 0 | 4.8e-05 | 4.8e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 4.8e-05 | 4.8e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199N_O1_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199N_O2_A, largest difference (A) | 0 | 1.7e-05 | 1.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.00034 | 0.00034 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000302 | 0.000302 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 7.9e-05 | 7.9e-05 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000302 | 0.000302 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 0.000512 | 0.000512 | 0.01 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_S_N2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: n2_volume, largest difference | 0 | 1.5e-05 | 1.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 2.7e-05 | 2.7e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration, every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| scaled/N1 | equilibration: lines of timeseries_equil_full.csv | 0 | 0 | 0 | 0 | construction | passed |
| scaled/N1 | equilibration: the steps with a file, in the order of the run | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 1: start and end (ps) | [0.0, 5.0] | [0.0, 5.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 1: solute and full frames read, and expected | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 0, 0, 0, 0] | [0, 0, 0, 0, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 1: last solute frame (frame, time, state, sign at N2) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 1: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 1: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 2: start and end (ps) | [5.0, 20.0] | [5.0, 20.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 2: solute and full frames read, and expected | [1, 1, 0, 0] | [1, 1, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 2: frames of state C, R, I; frames with positive and negative sign at N2 | [1, 0, 0, 1, 0] | [1, 0, 0, 1, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 2: last solute frame (frame, time, state, sign at N2) | [0, -0.085, 'C', 1] | [0, -0.085, 'C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 2: last solute frame, ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 2: last solute frame, signed volume at N2 | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 2: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 2: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 3: start and end (ps) | [20.0, 35.0] | [20.0, 35.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 3: solute and full frames read, and expected | [1, 1, 0, 0] | [1, 1, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 3: frames of state C, R, I; frames with positive and negative sign at N2 | [1, 0, 0, 1, 0] | [1, 0, 0, 1, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 3: last solute frame (frame, time, state, sign at N2) | [1, -0.07, 'C', 1] | [1, -0.07, 'C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 3: last solute frame, ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 3: last solute frame, signed volume at N2 | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 3: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 3: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 4: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 4: start and end (ps) | [35.0, 55.0] | [35.0, 55.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 4: solute and full frames read, and expected | [2, 2, 0, 0] | [2, 2, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 4: frames of state C, R, I; frames with positive and negative sign at N2 | [2, 0, 0, 2, 0] | [2, 0, 0, 2, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 4: last solute frame (frame, time, state, sign at N2) | [3, -0.045, 'C', 1] | [3, -0.045, 'C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 4: last solute frame, ligand RMSD (A) | 0.610426 | 0.610429 | 3e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 4: last solute frame, signed volume at N2 | 0.785515 | 0.785513 | 2e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 4: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 4: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 5: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 5: start and end (ps) | [55.0, 75.0] | [55.0, 75.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 5: solute and full frames read, and expected | [2, 2, 0, 0] | [2, 2, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | [2, 0, 0, 2, 0] | [2, 0, 0, 2, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 5: last solute frame (frame, time, state, sign at N2) | [5, -0.025, 'C', 1] | [5, -0.025, 'C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 5: last solute frame, ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 5: last solute frame, signed volume at N2 | 0.785515 | 0.785516 | 1e-06 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 5: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 5: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 6: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 6: start and end (ps) | [75.0, 100.0] | [75.0, 100.0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 6: solute and full frames read, and expected | [2, 2, 0, 0] | [2, 2, 0, 0] | – | equal | construction | passed |
| scaled/N1 | equilibration, stage 6: frames of state C, R, I; frames with positive and negative sign at N2 | [2, 0, 0, 2, 0] | [2, 0, 0, 2, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 6: last solute frame (frame, time, state, sign at N2) | [7, -0.005, 'C', 1] | [7, -0.005, 'C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 6: last solute frame, ligand RMSD (A) | 0.610426 | 0.610426 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 6: last solute frame, signed volume at N2 | 0.785515 | 0.785515 | 0 | 0.001 | second formulas | passed |
| scaled/N1 | equilibration, stage 6: last full frame (frame, time, water oxygens at the crystal ring positions) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration, stage 6: last full frame, water oxygens in the place of the ring | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration: frames of state C, R, I; positive and negative sign at N2 | [8, 0, 0, 8, 0] | [8, 0, 0, 8, 0] | – | equal | second formulas | passed |
| scaled/N1 | equilibration: first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration: first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration: first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| scaled/N1 | equilibration: state and sign at N2 of the last frame | ['C', 1] | ['C', 1] | – | equal | second formulas | passed |
| scaled/N1 | equilibration: time of the last solute frame (ns) | -0.005 | -0.005 | 0 | 1e-09 | construction | passed |
| scaled/N1 | equilibration: frames in which the ligand and the zinc were in another periodic image | [1, 1] | [1, 1] | – | equal | construction | passed |
| scaled/N1 | by design: frames with positive and with negative signed volume at N2 | [200, 0] | [200, 0] | – | equal | construction | passed |
| scaled/N1 | by design: positive and negative in the last 20 ns | [200, 0] | [200, 0] | – | equal | construction | passed |
| scaled/N1 | by design, frames of state C: number, positive, negative | [200, 200, 0] | [200, 200, 0] | – | equal | construction | passed |
| scaled/N1 | by design, frames of state R: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| scaled/N1 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| scaled/N1 | by design: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | construction | passed |
| scaled/N1 | by design: times of the sign changes that last at least 1 ns (ns) | [] | [] | – | equal | construction | passed |
| scaled/N1 | by design: sign after each of these changes | [] | [] | – | equal | construction | passed |
| scaled/N1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [None, None, None, None] | [None, None, None, None] | – | equal | construction | passed |
| scaled/N1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| scaled/N1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['C', False, 1, False] | ['C', False, 1, False] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration frames read | yes | yes | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, solute frames, full frames, length (ps) | [8, 0, 100.0] | [8, 0, 100.0] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, solute and full frames of each stage | [[0, 0], [1, 0], [1, 0], [2, 0], [2, 0], [2, 0]] | [[0, 0], [1, 0], [1, 0], [2, 0], [2, 0], [2, 0]] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, time of the first and of the last solute frame (ns) | [-0.085, -0.005] | [-0.085, -0.005] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, the times of all solute frames (ns) | [-0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005] | [-0.085, -0.07, -0.055, -0.045, -0.035, -0.025, -0.015, -0.005] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, water oxygens at the crystal ring positions in every full frame | [] | [] | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| scaled/N1 | by design: equilibration, state and sign at N2 of the last frame | ['C', 1] | ['C', 1] | – | equal | construction | passed |
| default/R1 | solute frames seen | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| default/R1 | manifest lists the segments | no | no | – | equal | construction | passed |
| default/R1 | manifest and files agree (none = no list) | none | none | – | equal | construction | passed |
| default/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| default/R1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| default/R1 | selections.json is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| default/R1 | reference.pdb is the file of the manifest (none = no checksum there) | none | none | – | equal | construction | passed |
| default/R1 | full frames seen | 20 | 20 | 0 | 0 | construction | passed |
| default/R1 | time seen (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| default/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| default/R1 | frames in which the ligand was in another periodic image | 29 | 29 | 0 | 0 | construction | passed |
| default/R1 | frames in which the zinc was in another periodic image | 18 | 18 | 0 | 0 | construction | passed |
| default/R1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | lines of timeseries_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| default/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| default/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000284 | 0.000284 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000269 | 0.000269 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| default/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| default/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| default/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| default/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| default/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| default/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| default/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| default/R1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | lines of timeseries_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| default/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| default/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| default/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| default/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| default/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| default/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| default/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | verdict | STAYED | STAYED | – | equal | construction | passed |
| default/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| default/R1 | run complete | yes | yes | – | equal | construction | passed |
| default/R1 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| default/R1 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| default/R1 | frames of the verdict window | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | count_C_window | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | count_R_window | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | count_C | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | count_R | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | rule 1: left the crystal pose | yes | yes | – | equal | construction | passed |
| default/R1 | count_C (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | count_R (design against the planted frames) | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | frames with the hydrogen bond N2-H...OG1 | 200 | 200 | 0 | 0 | second formulas | passed |
| default/R1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| default/R1 | frames with labelling identity | 200 | 200 | 0 | 0 | second formulas | passed |
| default/R1 | frames with labelling ring_flip | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | frames with labelling oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | frames with labelling ring_flip+oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | lig_rmsd_A: mean (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| default/R1 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | lig_rmsd_A: mean_last_half (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| default/R1 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_N2_A: mean (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_N2_A: mean_last_half (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | tau_fold_deg: mean (deg) | 142.744 | 142.744 | 0 | 0.01 | second formulas | passed |
| default/R1 | tau_fold_deg: sd (deg) | 0 | 8.5e-05 | 8.5e-05 | 0.01 | second formulas | passed |
| default/R1 | tau_fold_deg: mean_last_half (deg) | 142.744 | 142.744 | 1e-06 | 0.01 | second formulas | passed |
| default/R1 | tau_fold_deg: sd_last_half (deg) | 0 | 9.2e-05 | 9.2e-05 | 0.01 | second formulas | passed |
| default/R1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| default/R1 | ang_C4_S_N2_deg: mean (deg) | 109.432 | 109.432 | 3e-06 | 0.01 | second formulas | passed |
| default/R1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.6e-05 | 5.6e-05 | 0.01 | second formulas | passed |
| default/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 109.432 | 109.432 | 1e-06 | 0.01 | second formulas | passed |
| default/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 5.7e-05 | 5.7e-05 | 0.01 | second formulas | passed |
| default/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| default/R1 | d_S_N2_A: mean (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | d_S_N2_A: mean_last_half (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| default/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| default/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d1_N1_F131CZ_A: mean (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| default/R1 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| default/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d2_N1_N67OD1_A: mean (A) | 2.59678 | 2.59678 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 2.59678 | 2.59678 | 0 | 0.001 | second formulas | passed |
| default/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | n2_volume: mean | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| default/R1 | n2_volume: sd | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | n2_volume: mean_last_half | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| default/R1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE1_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE1_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE2_A: mean (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE2_A: sd (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE_min_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE_min_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| default/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| default/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| default/R1 | nearest water to the zinc, mean (A) | 3.6255 | 3.6255 | 1e-06 | 0.001 | construction | passed |
| default/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| default/R1 | frame 2, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| default/R1 | frame 2, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| default/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (A) | 3.31 | 3.3102 | 0.0002 | 0.0015 | starts.csv | passed |
| default/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| default/R1 | planned length (X0: 5 ns by the specification) (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| default/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| default/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| default/R1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | – | equal | second formulas | passed |
| default/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| default/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| default/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| default/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| default/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| default/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| default/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 6.62247 | 6.62246 | 1e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -2.55859 | -2.55859 | 2e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.13628 | -7.13628 | 5e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.85817 | 0.858167 | 3e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.961814 | 0.961813 | 1e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.647624 | 0.647624 | 0 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.631924 | 0.631921 | 3e-06 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 10.6445 | 10.6445 | 0 | 0.01 | second formulas | passed |
| default/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | gate G-V: solute frames counted for that fraction | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| default/R1 | gate G-V: full frames counted for that fraction | 20 | 20 | 0 | 0 | construction | passed |
| default/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 7.13628 | 7.13628 | 5e-06 | 0.01 | second formulas | passed |
| default/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| default/R1 | frames with ND2 nearer to the aniline N than OD1 | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | state C: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| default/R1 | state R: frames, and frames with ND2 nearer | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| default/R1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| default/R1 | campaign-table value: fraction of frames with ND2 nearer | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| default/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| default/R1 | N2: frames with positive, negative and zero signed volume | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| default/R1 | N2: fraction of frames with positive signed volume | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| default/R1 | N2: fraction of frames with negative signed volume | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| default/R1 | N2: frames of the window (last 20 ns) | 200 | 200 | 0 | 0 | construction | passed |
| default/R1 | N2: positive, negative and zero in the last 20 ns | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| default/R1 | N2: fraction positive in the last 20 ns | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| default/R1 | N2, frames of state C: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| default/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| default/R1 | N2, frames of state R: number, positive, negative, zero | [200, 0, 200, 0] | [200, 0, 200, 0] | – | equal | second formulas | passed |
| default/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| default/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| default/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| default/R1 | N2: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| default/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[-1, 0, 200]] | [[-1, 0, 200]] | – | equal | second formulas | passed |
| default/R1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| default/R1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| default/R1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| default/R1 | N2: sign of the first and of the last frame | [-1, -1] | [-1, -1] | – | equal | second formulas | passed |
| default/R1 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| default/R1 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| default/R1 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| default/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| default/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| default/R1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| default/R1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 3e-06 | 0.001 | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 2e-06 | 0.001 | second formulas | passed |
| default/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| default/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.010967, 4.003724, 0.0] | [4.010967, 4.003724, 0.0] | – | equal | second formulas | passed |
| default/R1 | water at the crystal ring positions: mean number within 3.5 A | 1.75 | 1.75 | 0 | 1e-06 | second formulas | passed |
| default/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 1.38237 | 1.38237 | 3e-06 | 0.001 | second formulas | passed |
| default/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 2.17059 | 2.17059 | 1e-06 | 0.001 | second formulas | passed |
| default/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | 0 | 1e-06 | second formulas | passed |
| default/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| default/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| default/R1 | water in the place of the ring: mean number | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| default/R1 | water in the place of the ring: largest number in a full frame | 2 | 2 | 0 | 1e-06 | second formulas | passed |
| default/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 1.38237 | 1.38237 | 3e-06 | 0.001 | second formulas | passed |
| default/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 2.19559 | 2.19559 | 0 | 0.001 | second formulas | passed |
| default/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [1.0, 1.382369] | [1.0, 1.382369] | – | equal | second formulas | passed |
| default/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 1.38237 | 1.38237 | – | equal | second formulas | passed |
| default/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| default/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| default/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| default/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| default/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| default/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| default/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.0, 0.0, 0] | [0.0, 0.0, 0] | – | equal | second formulas | passed |
| default/R1 | equilibration: frames read (files solute_eqNN.dcd present) | no | no | – | equal | construction | passed |
| default/R1 | equilibration: time series written (only with equilibration files) | [False, False] | [False, False] | – | equal | construction | passed |
| default/R1 | campaign-table value: solute frames of the equilibration | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | by design: frames with positive and with negative signed volume at N2 | [0, 200] | [0, 200] | – | equal | construction | passed |
| default/R1 | by design: positive and negative in the last 20 ns | [0, 200] | [0, 200] | – | equal | construction | passed |
| default/R1 | by design, frames of state C: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| default/R1 | by design, frames of state R: number, positive, negative | [200, 0, 200] | [200, 0, 200] | – | equal | construction | passed |
| default/R1 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| default/R1 | by design: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | construction | passed |
| default/R1 | by design: times of the sign changes that last at least 1 ns (ns) | [] | [] | – | equal | construction | passed |
| default/R1 | by design: sign after each of these changes | [] | [] | – | equal | construction | passed |
| default/R1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [None, None, None, None] | [None, None, None, None] | – | equal | construction | passed |
| default/R1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| default/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| default/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | 0 | 1e-06 | construction | passed |
| default/R1 | by design: smallest distance of a water oxygen to the crystal ring positions (A) | 1.382 | 1.38237 | 0.000369 | 0.001 | construction | passed |
| default/R1 | by design: water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2] | [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2] | – | equal | construction | passed |
| default/R1 | by design: mean number of water oxygens in the place of the ring | 1 | 1 | 0 | 1e-06 | construction | passed |
| default/R1 | by design: smallest distance of a water oxygen to the place of the ring (A) | 1.382 | 1.38237 | 0.000369 | 0.001 | construction | passed |
| default/R1 | by design: nearest water oxygen to the place of the ring in every full frame, largest difference (A) | 0 | 0.0004 | 0.0004 | 0.001 | construction | passed |
| default/R1 | by design: site atoms and cut-off of the count that the driver acts on (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| default/R1 | by design: where these settings came from, and whether they are the defaults of the specification | ['specification (defaults)', True] | ['specification (defaults)', True] | – | equal | construction | passed |
| default/R1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['R', False, -1, False] | ['R', False, -1, False] | – | equal | construction | passed |
| default/R1 | by design: equilibration frames read | no | no | – | equal | construction | passed |
| place/R1 | solute frames seen | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | DCD headers: segments contiguous (none = the headers tell nothing) | none | none | – | equal | construction | passed |
| place/R1 | manifest lists the segments | yes | yes | – | equal | construction | passed |
| place/R1 | manifest and files agree (none = no list) | yes | yes | – | equal | construction | passed |
| place/R1 | interval in the DCD headers as used | yes | yes | – | equal | construction | passed |
| place/R1 | where the system folder was found | manifest system_dir | manifest system_dir | – | equal | construction | passed |
| place/R1 | selections.json is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| place/R1 | reference.pdb is the file of the manifest (none = no checksum there) | yes | yes | – | equal | construction | passed |
| place/R1 | full frames seen | 20 | 20 | 0 | 0 | construction | passed |
| place/R1 | time seen (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| place/R1 | repeated frames found at a segment boundary | no | no | – | equal | construction | passed |
| place/R1 | frames in which the ligand was in another periodic image | 29 | 29 | 0 | 0 | construction | passed |
| place/R1 | frames in which the zinc was in another periodic image | 18 | 18 | 0 | 0 | construction | passed |
| place/R1 | warnings | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | lines of timeseries_solute.csv | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | columns of timeseries_solute.csv for Asn67 ND2 (only with the key asn67_ND2) | yes | yes | – | equal | construction | passed |
| place/R1 | every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.000745 | 0.000745 | 0.01 | second formulas | passed |
| place/R1 | every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000731 | 0.000731 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00033 | 0.00033 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000295 | 0.000295 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| place/R1 | every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| place/R1 | every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| place/R1 | every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| place/R1 | every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| place/R1 | every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| place/R1 | every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| place/R1 | every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: n2_sign_assigned, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| place/R1 | every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | lines of timeseries_full.csv | 20 | 20 | 0 | 0 | construction | passed |
| place/R1 | every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| place/R1 | every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| place/R1 | every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| place/R1 | every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| place/R1 | every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 2.8e-05 | 2.8e-05 | 0.001 | second formulas | passed |
| place/R1 | CA RMSF of 258 residues, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | construction | passed |
| place/R1 | residue numbers of rmsf_ca.csv that differ (3 to 125, 127 to 261) | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | verdict | STAYED | STAYED | – | equal | construction | passed |
| place/R1 | verdict provisional (run not finished) | no | no | – | equal | construction | passed |
| place/R1 | run complete | yes | yes | – | equal | construction | passed |
| place/R1 | first passage into the opposite state (ns) | none | none | – | equal | construction | passed |
| place/R1 | transitions | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | long visits (at least 1 ns) | 1 | 1 | 0 | 0 | construction | passed |
| place/R1 | frames of the verdict window | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | count_C_window | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | count_R_window | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | count_I_window | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | count_C | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | count_R | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | count_I | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | rule 1: left the crystal pose | yes | yes | – | equal | construction | passed |
| place/R1 | count_C (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | count_R (design against the planted frames) | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | count_I (design against the planted frames) | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | frames with the hydrogen bond N2-H...OG1 | 200 | 200 | 0 | 0 | second formulas | passed |
| place/R1 | full frames with a water within 2.6 A of the zinc | 1 | 1 | 0 | 0 | construction | passed |
| place/R1 | frames with labelling identity | 200 | 200 | 0 | 0 | second formulas | passed |
| place/R1 | frames with labelling ring_flip | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | frames with labelling oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | frames with labelling ring_flip+oxygen_swap | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | lig_rmsd_A: mean (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | lig_rmsd_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | lig_rmsd_A: mean_last_half (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | lig_rmsd_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | lig_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | lig_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_N2_A: mean (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_N2_A: sd (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_N2_A: mean_last_half (A) | 2.26515 | 2.26515 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_N2_A: sd_last_half (A) | 0.09 | 0.09 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d_Zn_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | tau_fold_deg: mean (deg) | 142.744 | 142.744 | 3e-06 | 0.01 | second formulas | passed |
| place/R1 | tau_fold_deg: sd (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| place/R1 | tau_fold_deg: mean_last_half (deg) | 142.744 | 142.744 | 2e-06 | 0.01 | second formulas | passed |
| place/R1 | tau_fold_deg: sd_last_half (deg) | 0 | 8.1e-05 | 8.1e-05 | 0.01 | second formulas | passed |
| place/R1 | tau_fold_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | tau_fold_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| place/R1 | ang_C4_S_N2_deg: mean (deg) | 109.432 | 109.432 | 2e-06 | 0.01 | second formulas | passed |
| place/R1 | ang_C4_S_N2_deg: sd (deg) | 0 | 5.5e-05 | 5.5e-05 | 0.01 | second formulas | passed |
| place/R1 | ang_C4_S_N2_deg: mean_last_half (deg) | 109.432 | 109.432 | 6e-06 | 0.01 | second formulas | passed |
| place/R1 | ang_C4_S_N2_deg: sd_last_half (deg) | 0 | 4.4e-05 | 4.4e-05 | 0.01 | second formulas | passed |
| place/R1 | ang_C4_S_N2_deg: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | ang_C4_S_N2_deg: block means, largest difference (deg) | 0 | 0 | 0 | 0.01 | second formulas | passed |
| place/R1 | d_S_N2_A: mean (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_S_N2_A: sd (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | d_S_N2_A: mean_last_half (A) | 1.67725 | 1.67725 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_S_N2_A: sd_last_half (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | d_S_N2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d_S_N2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | ca_rmsd_A: mean (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| place/R1 | ca_rmsd_A: sd (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | ca_rmsd_A: mean_last_half (A) | 0.171798 | 0.171798 | 0 | 0.001 | second formulas | passed |
| place/R1 | ca_rmsd_A: sd_last_half (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | ca_rmsd_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | ca_rmsd_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d1_N1_F131CZ_A: mean (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| place/R1 | d1_N1_F131CZ_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | d1_N1_F131CZ_A: mean_last_half (A) | 9.58495 | 9.58495 | 0 | 0.001 | second formulas | passed |
| place/R1 | d1_N1_F131CZ_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | d1_N1_F131CZ_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d1_N1_F131CZ_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d2_N1_N67OD1_A: mean (A) | 2.59678 | 2.59678 | 0 | 0.001 | second formulas | passed |
| place/R1 | d2_N1_N67OD1_A: sd (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | d2_N1_N67OD1_A: mean_last_half (A) | 2.59678 | 2.59678 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | d2_N1_N67OD1_A: sd_last_half (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | d2_N1_N67OD1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d2_N1_N67OD1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | n2_volume: mean | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | n2_volume: sd | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | n2_volume: mean_last_half | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | n2_volume: sd_last_half | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | n2_volume: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | n2_volume: block means, largest difference | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE1_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE1_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE1_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE1_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE1_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d_Zn_E106OE1_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE2_A: mean (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE2_A: sd (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE2_A: mean_last_half (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE2_A: sd_last_half (A) | 0.018963 | 0.018963 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE2_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d_Zn_E106OE2_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE_min_A: mean (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE_min_A: sd (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE_min_A: mean_last_half (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE_min_A: sd_last_half (A) | 0.021721 | 0.021721 | 0 | 0.001 | second formulas | passed |
| place/R1 | d_Zn_E106OE_min_A: number of 10 ns blocks | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | d_Zn_E106OE_min_A: block means, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | CA RMSD, mean (0.01 x radius of gyration of the CA atoms) (A) | 0.171798 | 0.171798 | 0 | 0.001 | construction | passed |
| place/R1 | nearest water to the zinc, mean (A) | 3.6255 | 3.6255 | 1e-06 | 0.001 | construction | passed |
| place/R1 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-06 | construction | passed |
| place/R1 | frame 2, lig_rmsd_A against starts.csv (R, rmsd_to_crystal_A) (A) | 3.573 | 3.5732 | 0.0002 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, d_Zn_N2_A against starts.csv (R, Zn_N2) (A) | 2.235 | 2.2351 | 0.0001 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, d_Zn_Osulfonyl_min_A against starts.csv (R, Zn_O_min) (A) | 3.907 | 3.9068 | 0.0002 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, d_N2_T199OG1_A against starts.csv (R, N2_OG1) (A) | 2.585 | 2.5852 | 0.0002 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, ang_N2_HN2_T199OG1_deg against starts.csv (R, N2_H_OG1) (deg) | 140.178 | 140.205 | 0.027 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, d_T199N_Osulfonyl_min_A against starts.csv (R, T199N_O_min) (A) | 3.02 | 3.0203 | 0.0003 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, tors_C3_C4_S_N2_deg against starts.csv (R, C3_C4_S_N2) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, tors_C5_C4_S_N2_deg against starts.csv (R, C5_C4_S_N2) (deg) | -34.865 | -34.865 | 0 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, tau_fold_deg against starts.csv (R, tau_fold) (deg) | 142.744 | 142.744 | 0 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, ang_C4_S_N2_deg against starts.csv (R, C4_S_N2) (deg) | 109.432 | 109.432 | 0 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, d_S_N2_A against starts.csv (R, S_N2) (A) | 1.677 | 1.6773 | 0.0003 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, ang_S_N2_Zn_deg against starts.csv (R, S_N2_Zn) (deg) | 119.529 | 119.529 | 0 | 0.06 | starts.csv | passed |
| place/R1 | frame 2, d1_N1_F131CZ_A against starts.csv (R, N1_Phe131CZ) (A) | 9.585 | 9.5849 | 0.0001 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, d2_N1_N67OD1_A against starts.csv (R, N1_Asn67OD1) (A) | 2.597 | 2.5968 | 0.0002 | 0.0015 | starts.csv | passed |
| place/R1 | frame 2, d_N1_N67ND2_A against starts.csv (R, N1_Asn67ND2) (A) | 3.31 | 3.3102 | 0.0002 | 0.0015 | starts.csv | passed |
| place/R1 | role of the run (X0 = validation) | production | production | – | equal | construction | passed |
| place/R1 | planned length (X0: 5 ns by the specification) (ns) | 2 | 2 | 0 | 1e-09 | construction | passed |
| place/R1 | gate G-V applies to this run (only to X0) | no | no | – | equal | construction | passed |
| place/R1 | rule 1 applies to this run (production runs of start X) | no | no | – | equal | construction | passed |
| place/R1 | gate G-V: the criteria that are not met | GV1 ligand RMSD, GV3 Zn-N2 | GV1 ligand RMSD, GV3 Zn-N2 | – | equal | second formulas | passed |
| place/R1 | gate G-V: all criteria met | no | no | – | equal | second formulas | passed |
| place/R1 | gate G-V: GV3 as a whole | no | no | – | equal | second formulas | passed |
| place/R1 | gate G-V: GV4 as a whole | yes | yes | – | equal | second formulas | passed |
| place/R1 | gate G-V: fifth-ligand criterion as a whole (both fractions) | yes | yes | – | equal | second formulas | passed |
| place/R1 | gate G-V: verdict (none for a run that is not X0) | none | none | – | equal | second formulas | passed |
| place/R1 | gate G-V: provisional (run not finished) | no | no | – | equal | construction | passed |
| place/R1 | gate G-V: ang_N2_Zn_H94NE2_deg, mean minus its value in 6RL9 (deg) | 6.62247 | 6.62246 | 3e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_N2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | -2.55859 | -2.55859 | 0 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_N2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | -7.13628 | -7.13628 | 3e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_H94NE2_Zn_H96NE2_deg, mean minus its value in 6RL9 (deg) | 0.85817 | 0.858171 | 1e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_H94NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.961814 | 0.961813 | 1e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_H96NE2_Zn_H119ND1_deg, mean minus its value in 6RL9 (deg) | 0.647624 | 0.647626 | 2e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: ang_C4_S_N2_deg, mean minus the reference (deg) | 0.631924 | 0.631926 | 2e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: tau_fold_deg, mean minus the reference (deg) | 10.6445 | 10.6445 | 3e-06 | 0.01 | second formulas | passed |
| place/R1 | gate G-V: solute frames with a sulfonyl oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | gate G-V: solute frames counted for that fraction | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | gate G-V: full frames with any oxygen within 2.6 A of the zinc | 1 | 1 | 0 | 0 | second formulas | passed |
| place/R1 | gate G-V: full frames counted for that fraction | 20 | 20 | 0 | 0 | construction | passed |
| place/R1 | campaign-table value: largest deviation of a mean N-Zn-N angle from 6RL9 (deg) | 7.13628 | 7.13628 | 3e-06 | 0.01 | second formulas | passed |
| place/R1 | Asn67 ND2 reported (the key asn67_ND2 is in selections.json) | yes | yes | – | equal | construction | passed |
| place/R1 | frames with ND2 nearer to the aniline N than OD1 | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | state C: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| place/R1 | state R: frames, and frames with ND2 nearer | [200, 0] | [200, 0] | – | equal | second formulas | passed |
| place/R1 | state I: frames, and frames with ND2 nearer | [0, 0] | [0, 0] | – | equal | second formulas | passed |
| place/R1 | campaign-table value: fraction of frames with ND2 nearer | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| place/R1 | N2: sign of the start as built (X and N +1, R -1) | -1 | -1 | 0 | 0 | construction | passed |
| place/R1 | N2: frames with positive, negative and zero signed volume | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| place/R1 | N2: fraction of frames with positive signed volume | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| place/R1 | N2: fraction of frames with negative signed volume | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| place/R1 | N2: frames of the window (last 20 ns) | 200 | 200 | 0 | 0 | construction | passed |
| place/R1 | N2: positive, negative and zero in the last 20 ns | [0, 200, 0] | [0, 200, 0] | – | equal | second formulas | passed |
| place/R1 | N2: fraction positive in the last 20 ns | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| place/R1 | N2, frames of state C: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| place/R1 | N2, frames of state C: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| place/R1 | N2, frames of state R: number, positive, negative, zero | [200, 0, 200, 0] | [200, 0, 200, 0] | – | equal | second formulas | passed |
| place/R1 | N2, frames of state R: fraction positive (none without a frame of that state) | 0 | 0 | 0 | 1e-06 | second formulas | passed |
| place/R1 | N2, frames of state I: number, positive, negative, zero | [0, 0, 0, 0] | [0, 0, 0, 0] | – | equal | second formulas | passed |
| place/R1 | N2, frames of state I: fraction positive (none without a frame of that state) | none | none | – | equal | second formulas | passed |
| place/R1 | N2: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | second formulas | passed |
| place/R1 | N2: the stretches of at least 1 ns (sign, first frame, frames) | [[-1, 0, 200]] | [[-1, 0, 200]] | – | equal | second formulas | passed |
| place/R1 | N2: number of sign changes that last at least 1 ns | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | N2: these sign changes (frame, from, to) | [] | [] | – | equal | second formulas | passed |
| place/R1 | N2: times of these sign changes (ns) | [] | [] | – | equal | second formulas | passed |
| place/R1 | N2: time of the first lasting sign change (ns) | none | none | – | equal | second formulas | passed |
| place/R1 | N2: sign of the first and of the last frame | [-1, -1] | [-1, -1] | – | equal | second formulas | passed |
| place/R1 | hydrogen bond around the first passage: reported (a first passage exists) | no | no | – | equal | second formulas | passed |
| place/R1 | hydrogen bond around the first passage: frame of the first passage | none | none | – | equal | second formulas | passed |
| place/R1 | hydrogen bond around the first passage: the six numbers (none without a first passage) | [None, None, None, None, None, None] | [None, None, None, None, None, None] | – | equal | construction | passed |
| place/R1 | campaign-table values: hydrogen-bond fraction before and after the first passage | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | Glu106: solute frames with OE1, with OE2, with either within 2.6 A of the zinc | [0, 0, 0] | [0, 0, 0] | – | equal | second formulas | passed |
| place/R1 | Glu106: smallest distance zinc to OE1 (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| place/R1 | Glu106: smallest distance zinc to OE2 (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| place/R1 | Glu106: full frames with a protein oxygen within 2.6 A of the zinc | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE1: mean distance (A) | 4.01097 | 4.01097 | 0 | 0.001 | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE1: smallest distance (A) | 4.00373 | 4.00372 | 6e-06 | 0.001 | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE1: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE2: mean distance (A) | 5.47641 | 5.47641 | 0 | 0.001 | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE2: smallest distance (A) | 5.47009 | 5.47009 | 5e-06 | 0.001 | second formulas | passed |
| place/R1 | Glu106 beside gate G-V, OE2: frames within 2.6 A, frames counted | [0, 200] | [0, 200] | – | equal | second formulas | passed |
| place/R1 | campaign-table values: zinc to OE1 of Glu106, mean, smallest, fraction within 2.6 A | [4.010967, 4.003721, 0.0] | [4.010967, 4.003721, 0.0] | – | equal | second formulas | passed |
| place/R1 | water at the crystal ring positions: mean number within 3.5 A | 1.75 | 1.75 | 0 | 1e-06 | second formulas | passed |
| place/R1 | water at the crystal ring positions: smallest distance of a water oxygen (A) | 1.38237 | 1.38237 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | water at the crystal ring positions: mean distance of the nearest water oxygen (A) | 2.17059 | 2.17059 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | campaign-table value: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | 0 | 1e-06 | second formulas | passed |
| place/R1 | water in the place of the ring: site atoms and cut-off used (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| place/R1 | water in the place of the ring: where the settings came from, and whether they are the defaults | ['manifest hydration_check.settings', True] | ['manifest hydration_check.settings', True] | – | equal | construction | passed |
| place/R1 | water in the place of the ring: mean number | 1 | 1 | 0 | 1e-06 | second formulas | passed |
| place/R1 | water in the place of the ring: largest number in a full frame | 2 | 2 | 0 | 1e-06 | second formulas | passed |
| place/R1 | water in the place of the ring: smallest distance of a water oxygen (A) | 1.38237 | 1.38237 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | water in the place of the ring: mean distance of the nearest water oxygen (A) | 2.19559 | 2.19559 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | campaign-table values: water in the place of the ring, mean number and smallest distance | [1.0, 1.38237] | [1.0, 1.38237] | – | equal | second formulas | passed |
| place/R1 | campaign-table value: smallest distance of a water oxygen to the seven positions of the wider count | 1.38237 | 1.38237 | – | equal | second formulas | passed |
| place/R1 | first production frame: state, and whether it differs from the state of the start | ['R', False] | ['R', False] | – | equal | second formulas | passed |
| place/R1 | first production frame: in the opposite state | no | no | – | equal | second formulas | passed |
| place/R1 | first production frame: sign at N2, and whether it differs from the sign of the start | [-1, False] | [-1, False] | – | equal | second formulas | passed |
| place/R1 | first production frame: signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | first production frame: ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | campaign-table values: state of the first frame, differs from the start, sign at N2 | ['R', False, -1] | ['R', False, -1] | – | equal | second formulas | passed |
| place/R1 | campaign-table values: N2 fraction positive, the same in the last 20 ns, lasting sign changes | [0.0, 0.0, 0] | [0.0, 0.0, 0] | – | equal | second formulas | passed |
| place/R1 | equilibration: frames read (files solute_eqNN.dcd present) | yes | yes | – | equal | construction | passed |
| place/R1 | equilibration: time series written (only with equilibration files) | [True, True] | [True, True] | – | equal | construction | passed |
| place/R1 | campaign-table value: solute frames of the equilibration | 350 | 350 | 0 | 0 | construction | passed |
| place/R1 | equilibration: solute frames and full frames read | [350, 35] | [350, 35] | – | equal | construction | passed |
| place/R1 | equilibration: length of the equilibration (ps) | 3500 | 3500 | 0 | 1e-09 | construction | passed |
| place/R1 | equilibration: where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| place/R1 | equilibration: lines of timeseries_equil_solute.csv | 350 | 350 | 0 | 0 | construction | passed |
| place/R1 | equilibration: no column that depends on the frames before (state_assigned, n2_sign_assigned) | yes | yes | – | equal | construction | passed |
| place/R1 | equilibration, every frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | equilibration, every frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| place/R1 | equilibration: no time above zero | yes | yes | – | equal | construction | passed |
| place/R1 | equilibration, every frame: lig_rmsd_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: lig_rmsd_labelling, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: lig_rmsd_identity_A, largest difference (A) | 0 | 2.5e-05 | 2.5e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d1_N1_F131CZ_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d2_N1_N67OD1_A, largest difference (A) | 0 | 1.9e-05 | 1.9e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: state, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_N2_T199OG1_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_N2_HN2_T199OG1_deg, largest difference (deg) | 0 | 0.001255 | 0.001255 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: hb_N2H_T199OG1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_N2_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_H94NE2_A, largest difference (A) | 0 | 3.6e-05 | 3.6e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_H96NE2_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_H119ND1_A, largest difference (A) | 0 | 4.2e-05 | 4.2e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_N2_Zn_H94NE2_deg, largest difference (deg) | 0 | 0.000716 | 0.000716 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_N2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.000589 | 0.000589 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_N2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000269 | 0.000269 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_H94NE2_Zn_H96NE2_deg, largest difference (deg) | 0 | 0.00067 | 0.00067 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_H94NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000705 | 0.000705 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_H96NE2_Zn_H119ND1_deg, largest difference (deg) | 0 | 0.000645 | 0.000645 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_S_N2_Zn_deg, largest difference (deg) | 0 | 0.000562 | 0.000562 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_O1_A, largest difference (A) | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_O2_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_Osulfonyl_min_A, largest difference (A) | 0 | 5.1e-05 | 5.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: fifth_sulfonylO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199N_O1_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199N_O2_A, largest difference (A) | 0 | 3.2e-05 | 3.2e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_T199N_H_O1_deg, largest difference (deg) | 0 | 0.000237 | 0.000237 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_T199N_H_O2_deg, largest difference (deg) | 0 | 0.000132 | 0.000132 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199N_Osulfonyl_min_A, largest difference (A) | 0 | 5.5e-05 | 5.5e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: hb_T199NH_Osulfonyl, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199OG1_E106OE1_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199OG1_E106OE2_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE1_deg, largest difference (deg) | 0 | 0.000382 | 0.000382 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_T199OG1_HG1_E106OE2_deg, largest difference (deg) | 0 | 1e-05 | 1e-05 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: d_T199OG1_E106OE_min_A, largest difference (A) | 0 | 4.1e-05 | 4.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: hb_T199OG1H_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: tors_C3_C4_S_N2_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: tors_C5_C4_S_N2_deg, largest difference (deg) | 0 | 0.000522 | 0.000522 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: tau_fold_deg, largest difference (deg) | 0 | 0.000519 | 0.000519 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: ang_C4_S_N2_deg, largest difference (deg) | 0 | 7.6e-05 | 7.6e-05 | 0.01 | second formulas | passed |
| place/R1 | equilibration, every frame: d_S_N2_A, largest difference (A) | 0 | 5e-05 | 5e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: ca_rmsd_A, largest difference (A) | 0 | 2e-06 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: n2_volume, largest difference | 0 | 1.8e-05 | 1.8e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: n2_sign, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_E106OE1_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_E106OE2_A, largest difference (A) | 0 | 1e-05 | 1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: d_Zn_E106OE_min_A, largest difference (A) | 0 | 3.1e-05 | 3.1e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: fifth_E106OE1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: fifth_E106OE2, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: fifth_E106OE, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every frame: d_N1_N67ND2_A, largest difference (A) | 0 | 2.6e-05 | 2.6e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every frame: asn67_ND2_nearer_N1, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration: lines of timeseries_equil_full.csv | 35 | 35 | 0 | 0 | construction | passed |
| place/R1 | equilibration, every full frame: stage and extension step, frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | equilibration, every full frame: time, largest difference (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| place/R1 | equilibration, every full frame: d_Zn_Owater_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: water_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_Zn_Oprotein_min_A, largest difference (A) | 0 | 3e-06 | 3e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: protein_O_index_nearest_Zn, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_Zn_Osulfonyl_min_full_A, largest difference (A) | 0 | 1e-06 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_Zn_O_any_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: fifth_water, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: fifth_proteinO, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: fifth_any, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: n_water_N1_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_Zn_N2_full_A, largest difference (A) | 0 | 5.3e-05 | 5.3e-05 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: n_water_ring_wide_3p5A, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_water_ring_wide_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, every full frame: n_water_ring_place, frames that differ | 0 | 0 | 0 | 0 | second formulas | passed |
| place/R1 | equilibration, every full frame: d_water_ring_place_min_A, largest difference (A) | 0 | 0 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration: the steps with a file, in the order of the run | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq03x3', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq03x3', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| place/R1 | equilibration, stage 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | equilibration, stage 1: start and end (ps) | [0.0, 100.0] | [0.0, 100.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 1: solute and full frames read, and expected | [10, 10, 1, 1] | [10, 10, 1, 1] | – | equal | construction | passed |
| place/R1 | equilibration, stage 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 10, 0, 0, 10] | [0, 10, 0, 0, 10] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 1: last solute frame (frame, time, state, sign at N2) | [9, -3.4, 'R', -1] | [9, -3.4, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 1: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 1: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [0, -3.4, 1] | [0, -3.4, 1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the place of the ring (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | equilibration, stage 2: start and end (ps) | [100.0, 400.0] | [100.0, 400.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 2: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| place/R1 | equilibration, stage 2: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 2: last solute frame (frame, time, state, sign at N2) | [39, -3.1, 'R', -1] | [39, -3.1, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 2: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 2: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [3, -3.1, 1] | [3, -3.1, 1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 2: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the place of the ring (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [0, True] | [0, True] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3: start and end (ps) | [400.0, 700.0] | [400.0, 700.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3: solute and full frames read, and expected | [30, 30, 3, 3] | [30, 30, 3, 3] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 30, 0, 0, 30] | [0, 30, 0, 0, 30] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3: last solute frame (frame, time, state, sign at N2) | [69, -2.8, 'R', -1] | [69, -2.8, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3: last full frame (frame, time, water oxygens at the crystal ring positions) | [6, -2.8, 1] | [6, -2.8, 1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the place of the ring (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.6 | 2.6 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [0, True] | [0, True] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 1: start and end (ps) | [700.0, 1200.0] | [700.0, 1200.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 1: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 1: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last solute frame (frame, time, state, sign at N2) | [119, -2.3, 'R', -1] | [119, -2.3, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last full frame (frame, time, water oxygens at the crystal ring positions) | [11, -2.3, 1] | [11, -2.3, 1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last full frame, nearest water oxygen to the place of the ring (A) | 2.6 | 2.59999 | 8e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 1: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.6 | 2.59999 | 8e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [0, True] | [0, True] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 2: start and end (ps) | [1200.0, 1700.0] | [1200.0, 1700.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 2: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 2: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last solute frame (frame, time, state, sign at N2) | [169, -1.8, 'R', -1] | [169, -1.8, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last solute frame, signed volume at N2 | -0.710418 | -0.710417 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last full frame (frame, time, water oxygens at the crystal ring positions) | [16, -1.8, 1] | [16, -1.8, 1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last full frame, water oxygens in the place of the ring | 0 | 0 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last full frame, nearest water oxygen to the place of the ring (A) | 2.6 | 2.6 | 3e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 2: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.6 | 2.6 | 3e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [1, True] | [1, True] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 3: start and end (ps) | [1700.0, 2200.0] | [1700.0, 2200.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 3: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| place/R1 | equilibration, stage 3, extension step 3: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last solute frame (frame, time, state, sign at N2) | [219, -1.3, 'R', -1] | [219, -1.3, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last full frame (frame, time, water oxygens at the crystal ring positions) | [21, -1.3, 2] | [21, -1.3, 2] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last full frame, nearest water oxygen to the place of the ring (A) | 2.4 | 2.4 | 4e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 3, extension step 3: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.3 | 2.3 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 4: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | equilibration, stage 4: start and end (ps) | [2200.0, 2600.0] | [2200.0, 2600.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 4: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| place/R1 | equilibration, stage 4: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 4: last solute frame (frame, time, state, sign at N2) | [259, -0.9, 'R', -1] | [259, -0.9, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 4: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 4: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 4: last full frame (frame, time, water oxygens at the crystal ring positions) | [25, -0.9, 2] | [25, -0.9, 2] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 4: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the place of the ring (A) | 2.4 | 2.4 | 2e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 4: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.3 | 2.3 | 1e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 5: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | equilibration, stage 5: start and end (ps) | [2600.0, 3000.0] | [2600.0, 3000.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 5: solute and full frames read, and expected | [40, 40, 4, 4] | [40, 40, 4, 4] | – | equal | construction | passed |
| place/R1 | equilibration, stage 5: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 40, 0, 0, 40] | [0, 40, 0, 0, 40] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 5: last solute frame (frame, time, state, sign at N2) | [299, -0.5, 'R', -1] | [299, -0.5, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 5: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 5: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 5: last full frame (frame, time, water oxygens at the crystal ring positions) | [29, -0.5, 2] | [29, -0.5, 2] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 5: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the place of the ring (A) | 2.4 | 2.4 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 5: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.3 | 2.3 | 3e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 6: the count of the driver in the manifest (none: not recorded), and whether the analysis agrees with it | [None, None] | [None, None] | – | equal | construction | passed |
| place/R1 | equilibration, stage 6: start and end (ps) | [3000.0, 3500.0] | [3000.0, 3500.0] | – | equal | construction | passed |
| place/R1 | equilibration, stage 6: solute and full frames read, and expected | [50, 50, 5, 5] | [50, 50, 5, 5] | – | equal | construction | passed |
| place/R1 | equilibration, stage 6: frames of state C, R, I; frames with positive and negative sign at N2 | [0, 50, 0, 0, 50] | [0, 50, 0, 0, 50] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 6: last solute frame (frame, time, state, sign at N2) | [349, 0.0, 'R', -1] | [349, 0.0, 'R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 6: last solute frame, ligand RMSD (A) | 3.57322 | 3.57322 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 6: last solute frame, signed volume at N2 | -0.710418 | -0.710418 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 6: last full frame (frame, time, water oxygens at the crystal ring positions) | [34, 0.0, 2] | [34, 0.0, 2] | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 6: last full frame, water oxygens in the place of the ring | 1 | 1 | – | equal | second formulas | passed |
| place/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the place of the ring (A) | 2.4 | 2.4 | 0 | 0.001 | second formulas | passed |
| place/R1 | equilibration, stage 6: last full frame, nearest water oxygen to the seven positions of the wider count (A) | 2.3 | 2.3 | 4e-06 | 0.001 | second formulas | passed |
| place/R1 | equilibration: frames of state C, R, I; positive and negative sign at N2 | [0, 350, 0, 0, 350] | [0, 350, 0, 0, 350] | – | equal | second formulas | passed |
| place/R1 | equilibration: first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| place/R1 | equilibration: first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| place/R1 | equilibration: first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | second formulas | passed |
| place/R1 | equilibration: state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | second formulas | passed |
| place/R1 | equilibration: time of the last solute frame (ns) | 0 | 0 | 0 | 1e-09 | construction | passed |
| place/R1 | equilibration: frames in which the ligand and the zinc were in another periodic image | [50, 32] | [50, 32] | – | equal | construction | passed |
| place/R1 | by design: frames with positive and with negative signed volume at N2 | [0, 200] | [0, 200] | – | equal | construction | passed |
| place/R1 | by design: positive and negative in the last 20 ns | [0, 200] | [0, 200] | – | equal | construction | passed |
| place/R1 | by design, frames of state C: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| place/R1 | by design, frames of state R: number, positive, negative | [200, 0, 200] | [200, 0, 200] | – | equal | construction | passed |
| place/R1 | by design, frames of state I: number, positive, negative | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| place/R1 | by design: stretches of one sign, and those of at least 1 ns | [1, 1] | [1, 1] | – | equal | construction | passed |
| place/R1 | by design: times of the sign changes that last at least 1 ns (ns) | [] | [] | – | equal | construction | passed |
| place/R1 | by design: sign after each of these changes | [] | [] | – | equal | construction | passed |
| place/R1 | by design: hydrogen bond before the first passage (frames, with the bond) and after it (frames, with the bond) | [None, None, None, None] | [None, None, None, None] | – | equal | construction | passed |
| place/R1 | by design: solute frames with OE1, with OE2, with either oxygen of Glu106 within 2.6 A | [0, 0, 0] | [0, 0, 0] | – | equal | construction | passed |
| place/R1 | by design: water oxygens at the crystal ring positions, full frames that differ | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | by design: mean number of water oxygens at the crystal ring positions | 1.75 | 1.75 | 0 | 1e-06 | construction | passed |
| place/R1 | by design: smallest distance of a water oxygen to the crystal ring positions (A) | 1.382 | 1.38237 | 0.00037 | 0.001 | construction | passed |
| place/R1 | by design: water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2] | [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2] | – | equal | construction | passed |
| place/R1 | by design: mean number of water oxygens in the place of the ring | 1 | 1 | 0 | 1e-06 | construction | passed |
| place/R1 | by design: smallest distance of a water oxygen to the place of the ring (A) | 1.382 | 1.38237 | 0.00037 | 0.001 | construction | passed |
| place/R1 | by design: nearest water oxygen to the place of the ring in every full frame, largest difference (A) | 0 | 0.0004 | 0.0004 | 0.001 | construction | passed |
| place/R1 | by design: site atoms and cut-off of the count that the driver acts on (A) | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | [['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], 2.5] | – | equal | construction | passed |
| place/R1 | by design: where these settings came from, and whether they are the defaults of the specification | ['manifest hydration_check.settings', True] | ['manifest hydration_check.settings', True] | – | equal | construction | passed |
| place/R1 | by design: state of the first production frame, differs from the start; sign at N2, differs from the start | ['R', False, -1, False] | ['R', False, -1, False] | – | equal | construction | passed |
| place/R1 | by design: equilibration frames read | yes | yes | – | equal | construction | passed |
| place/R1 | by design: equilibration, solute frames, full frames, length (ps) | [350, 35, 3500.0] | [350, 35, 3500.0] | – | equal | construction | passed |
| place/R1 | by design: equilibration, where the stage table came from | manifest equilibration.stages | manifest equilibration.stages | – | equal | construction | passed |
| place/R1 | by design: equilibration, solute and full frames of each stage | [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]] | [[10, 1], [30, 3], [30, 3], [50, 5], [50, 5], [50, 5], [40, 4], [40, 4], [50, 5]] | – | equal | construction | passed |
| place/R1 | by design: equilibration, time of the first and of the last solute frame (ns) | [-3.49, 0.0] | [-3.49, 0.0] | – | equal | construction | passed |
| place/R1 | by design: equilibration, time of solute frame 10 (the first of stage 2) (ns) | [-3.39, '2'] | [-3.39, '2'] | – | equal | construction | passed |
| place/R1 | by design: equilibration, water oxygens at the crystal ring positions in every full frame | [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] | [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] | – | equal | construction | passed |
| place/R1 | by design: equilibration, water oxygens at the crystal ring positions at the end of stage 3 (the held part) | 1 | 1 | 0 | 0 | construction | passed |
| place/R1 | by design: equilibration, water oxygens in the place of the ring in every full frame | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] | – | equal | construction | passed |
| place/R1 | by design: equilibration, nearest water oxygen to the place of the ring in every full frame, largest difference (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| place/R1 | by design: equilibration, water oxygens in the place of the ring at the end of stage 3 | 0 | 0 | 0 | 0 | construction | passed |
| place/R1 | by design: equilibration, first_frame_not_in_state_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| place/R1 | by design: equilibration, first_frame_in_opposite_state (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| place/R1 | by design: equilibration, first_frame_n2_sign_not_of_start (frame, stage, extension step, time) | none | none | – | equal | construction | passed |
| place/R1 | by design: equilibration, state and sign at N2 of the last frame | ['R', -1] | ['R', -1] | – | equal | construction | passed |
| place/R1 | by design: equilibration, the steps in the order of the run | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq03x3', 'eq04', 'eq05', 'eq06'] | ['eq01', 'eq02', 'eq03', 'eq03x1', 'eq03x2', 'eq03x3', 'eq04', 'eq05', 'eq06'] | – | equal | construction | passed |
| place/R1 | by design: equilibration, stage, extension step and time of solute frame 70 (the first of eq03x1) | [3, 1, -2.79] | [3, 1, -2.79] | – | equal | construction | passed |
| place/R1 | by design: equilibration, stage, extension step and time of solute frame 170 (the first of eq03x3) | [3, 3, -1.79] | [3, 3, -1.79] | – | equal | construction | passed |
| place/R1 | by design: equilibration, water oxygens in the place of the ring (the count the driver acts on) at the end of eq03, eq03x1, eq03x2, eq03x3 | [0, 0, 0, 1] | [0, 0, 0, 1] | – | equal | construction | passed |
| place/R1 | by design: equilibration, these steps were compared with the driver (their last full frame is at the end of the step) | [True, True, True, True] | [True, True, True, True] | – | equal | construction | passed |
| place/R1 | by design: equilibration, the count of the driver at the end of these steps agrees with the analysis | [True, True, True, True] | [True, True, True, True] | – | equal | construction | passed |
| place/R1 | by design: equilibration, the nearest water oxygen of the driver at the end of these steps agrees with the analysis | [True, True, True, True] | [True, True, True, True] | – | equal | construction | passed |
| test/X2 | ligand RMSD, mean (labels swapped in 3 of 4 frames) (A) | 0 | 1e-06 | 1e-06 | 0.001 | construction | passed |
| test/X2 | ligand RMSD, largest value (A) | 0 | 1.1e-05 | 1.1e-05 | 0.001 | construction | passed |
| test/X2 | ligand RMSD without the symmetry correction, largest value (must be large) | yes | yes | – | equal | construction | passed |
| test/X2 | hydrogen bond N2-H...OG1, fraction of frames | 0.8 | 0.8 | 0 | 1e-09 | construction | passed |
| test/X2 | hydrogen bond N2-H...OG1, fraction in the last half | 0.8 | 0.8 | 0 | 1e-09 | construction | passed |
| test/X2 | hydrogen bond N2-H...OG1, the three block fractions | [0.8, 0.8, 0.8] | [0.8, 0.8, 0.8] | – | equal | construction | passed |
| test/X2 | Zn-N2, mean (starts.csv 1.931 + 0.300 in every tenth frame) (A) | 1.961 | 1.96144 | 0.000442 | 0.0015 | starts.csv | passed |
| test/X2 | Zn-N2, largest value (starts.csv 1.931 + 0.300) (A) | 2.231 | 2.23145 | 0.000451 | 0.0015 | starts.csv | passed |
| test/X2 | nearest water to the zinc, mean (0.1 x 2.400 + 0.9 x 3.690) (A) | 3.561 | 3.561 | 0 | 0.001 | construction | passed |
| test/X2 | nearest water to the zinc, smallest (A) | 2.4 | 2.39999 | 1.2e-05 | 0.001 | construction | passed |
| test/X2 | fraction of full frames with a fifth ligand (water) | 0.1 | 0.1 | 0 | 1e-09 | construction | passed |
| test/X2 | water molecules within 3.5 A of the aniline N, mean | 1.5 | 1.5 | 0 | 1e-09 | construction | passed |
| test/X2 | tau_fold, mean (starts.csv, unchanged by the label swaps) (deg) | 132.056 | 132.056 | 7.1e-05 | 0.06 | starts.csv | passed |
| test/X2 | tau_fold, standard deviation (deg) | 0 | 6.9e-05 | 6.9e-05 | 0.01 | construction | passed |
| test/X2 | gate G-V, criterion 1 (mean ligand RMSD at most 1.0 A) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 2 (hydrogen bond in at least 80 % of frames; here exactly 80 %) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 3, no fifth ligand (planted: water at 2.400 A in 10 % of the full frames) | no | no | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 3, fifth ligand in the solute trajectory alone (no sulfonyl oxygen near) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 3, full frames with an oxygen within 2.6 A, of 300 | 30 | 30 | 0 | 0 | construction | passed |
| test/X2 | gate G-V, criterion 3, Zn-N2 window | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 3, angle His96 NE2-Zn-His119 ND1 within 10 deg of 98.4 (version 1.0: in 100-120 deg, not met) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V, criterion 4 (C4-S-N2 within 3 deg of 108.8, tau_fold within 20 deg of 132.1) | yes | yes | – | equal | construction | passed |
| test/X2 | gate G-V passed | no | no | – | equal | construction | passed |
| test/X2 | gate G-V verdict (a production run has none) | none | none | – | equal | construction | passed |
| test/X1 | fraction of the last 20 ns in the opposite state (boundary of the 90 % rule) | 0.9 | 0.9 | 0 | 1e-09 | construction | passed |
| test/X1 | tau_fold, mean (rigid moves and label swaps keep the crystal value of starts.csv) (deg) | 132.056 | 132.056 | 7e-05 | 0.06 | starts.csv | passed |
| test/X1 | S-N2, standard deviation (rigid moves) (A) | 0 | 1e-06 | 1e-06 | 0.001 | construction | passed |
| test/X1 | frame 960 (ligand translated by 3.000 A): ligand RMSD (A) | 3 | 3 | 0 | 0.001 | construction | passed |
| test/X1 | frame 960: d1 (starts.csv 3.897 + 3.000) (A) | 6.897 | 6.8966 | 0.0004 | 0.0015 | starts.csv | passed |
| test/X1 | frame 960: state | I | I | – | equal | construction | passed |
| test/X1 | frame 960: assigned state (carried from the visit it lies in) | R | R | – | equal | construction | passed |
| test/X1 | frames 1602 and 1603 (labels swapped in one): same ligand RMSD (A) | 0 | 0 | 0 | 0.001 | construction | passed |
| test/X1 | frames 1602 and 1603: labelling differs | yes | yes | – | equal | construction | passed |
| test/X1 | frame 320 (short excursion, 0.5 ns): assigned state | R | R | – | equal | construction | passed |
| test/N1 | verdict window complete (20 ns) | no | no | – | equal | construction | passed |
| test/N1 | ligand RMSD, mean (starts.csv, near-native pose) (A) | 0.61 | 0.610426 | 0.000426 | 0.0015 | starts.csv | passed |
| campaign | n_runs | 5 | 5 | 0 | 0 | construction | passed |
| campaign | n_complete | 4 | 4 | 0 | 0 | construction | passed |
| campaign | all_nine_complete | no | no | – | equal | construction | passed |
| campaign | final | no | no | – | equal | construction | passed |
| campaign | rule1.n_X_runs | 3 | 3 | 0 | 0 | construction | passed |
| campaign | rule1.n_X_left_crystal | 2 | 2 | 0 | 0 | construction | passed |
| campaign | rule1.failed | yes | yes | – | equal | construction | passed |
| campaign | n_converted.X | 1 | 1 | 0 | 0 | construction | passed |
| campaign | n_converted.N | 0 | 0 | 0 | 0 | construction | passed |
| campaign | n_converted.R | 1 | 1 | 0 | 0 | construction | passed |
| campaign | n_stayed.X | 1 | 1 | 0 | 0 | construction | passed |
| campaign | n_stayed.N | 1 | 1 | 0 | 0 | construction | passed |
| campaign | n_stayed.R | 0 | 0 | 0 | 0 | construction | passed |
| campaign | n_mixed.X | 1 | 1 | 0 | 0 | construction | passed |
| campaign | outcome | validation_failed | validation_failed | – | equal | construction | passed |
| campaign | rule4_sensitivity_run_required | no | no | – | equal | construction | passed |
| campaign | folders left out of the campaign (the validation run X0) | test/X0 | test/X0 | – | equal | construction | passed |
| campaign | runs that could not be analysed | 0 | 0 | 0 | 0 | construction | passed |
| campaign | lines of campaign_table.csv | 5 | 5 | 0 | 0 | construction | passed |
| campaign | order of the runs | test/X1 test/X2 test/X3 test/N1 test/R1 | test/X1 test/X2 test/X3 test/N1 test/R1 | – | equal | construction | passed |
