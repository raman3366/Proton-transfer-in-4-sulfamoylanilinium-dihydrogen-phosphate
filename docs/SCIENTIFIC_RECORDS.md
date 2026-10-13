# Which records support the reported interpretation

## Molecular and periodic DFT

The reported molecular scan is the eight completed records under `supporting/06_DFT_results/proton_scan/`, summarised in `scan_summary.json` and Table S7. The A-chain point at 1.25 Å has `converged: false`, 61 steps, and is retained with the SI's force/energy qualification. The primary PCM pair and the no-dispersion pair were accepted on negligible forces despite missing full displacement convergence; their `opt_summary.json` and SI S1 describe that decision. They must not be described as unconditionally optimiser-converged.

`functional_check/fsp_scan.json` is the finished five-method comparison. `fsp_scan_partial.json` is a checkpointed intermediate. `eps_ladder/eps_ladder.json` records the geometry change between ε = 4 and 8; it does not establish a thermodynamic crossover by comparing energies at different ε.

For the periodic conclusion use `periodic_QE/profile.json`, `periodic_qe.json`, `convergence.json` and their linked `runs/*/*.out`. Table S15 is the constrained N–H profile. `relax_iirfix` is a separate O–H-constrained configuration, not the profile's neutral minimum. `_pre` records are coarse pre-relaxations; `hybrid_estimate.json` is a historical cost estimate. Neither supports final energy conclusions. The fuller package contains 49 QE outputs, four more than the original ZIP's 45 inputs: `relax_iifix`, `relax_iifix_pre`, `pbe_iifix_sg15_e100`, `hybrid_iifix_sg15_e100`. These additional records do not create four additional reported production results.

The two 135-atom clusters, their results in `sfa_xclu_{salt,neutral}/`, and `crystal_cluster_qtaim/qtaim_nine_bonds.*` support the cluster/interaction analysis (Table S22). A finite vacuum cluster is not the periodic ionic lattice. Two large Molden files are losslessly gzip-compressed in this repository; their decompressed bytes are checked against the original source manifest.

## MD run roles and clocks

| Role | Runs | Duration | Interpretation |
|---|---|---|---|
| Primary | `bonded_2e/X1–X3`, `N1–N3`, `R1–R3` | 9 × 100 ns | Three starts × three replicas; only these nine enter the primary count |
| Sensitivity | `bonded_2s/R1` | 100 ns | Alternative zinc set, which failed part of its validation |
| Sensitivity | `bonded_2e_n2fit/R1` | 100 ns | Candidate nitrogen terms that were not accepted; MIXED verdict |
| Validation | `bonded_2e/X0`, `bonded_2s/X0` | 2 × 5 ns | Parameter-set validation, not production |
| Other/support | `bench`, `preflight_a`, `preflight_b`, system builds and historical test records | Various | Do not pool with reported runs |

See [md_run_index.csv](../data/derived/md_run_index.csv). The eleven compact series each contain 130 released-equilibration frames followed by 10,000 production frames. Production time starts at 0; ligand release occurred 1.3 ns earlier. Table S39 uses time since release. Read `conversion_times.json` for that clock, not just the campaign CSV's `first_passage_ns`. Fractions of frames are not equilibrium populations or rate constants.

## Provisional and completed zinc-model quantum chemistry

All paths below are under `molecular_dynamics/quantum_chemistry/`:

| Record | Status / use |
|---|---|
| `gate_G-R/relax/qm_Rx.json`, `qm_Rx_last.xyz` | Provisional, step 23; `state: running`; retained for history |
| `gate_G-R/relax/qm_Rx_before_crossing.xyz` | Earlier step 15, not an endpoint |
| `gate_G-R/relax/qm_Rx_step63_before_convergence.xyz` | Later but still unconverged geometry |
| `gate_G-R/relax/qm_X.json`, `qm_X.xyz` | Completed X relaxation, 55 steps |
| `gate_G-R/relax/qm_Rx_converged_step66.json`, corresponding XYZ | Completed Rx relaxation, 66 steps; inverted nitrogen |
| `gate_G-R/relax/relax_qm_both_converged.json` | Completed endpoint comparison |
| `gate_G-R/relax/gr_gate_after_Rx_converged.json`, `binding_check_final.json` | Later gate and non-provisional failed binding check |
| `gate_G-R2/behaviour.json`, `n2_gate_result.json` | Finished G-R2 interpretation used for manuscript facts |

Earlier `gate_G-R/relax_qm.json`, `gr_gate.json` and `binding_check.json` remain historical records. The detailed [relaxation README](../molecular_dynamics/quantum_chemistry/gate_G-R/relax/README.md) explains their relationship. G-R and G-R2 failed; no refitted nitrogen terms were accepted for the nine primary runs. The small vacuum zinc-site model contains neither the pocket nor Thr199/water, so its relaxation is not a simulation of protein binding.

The recovered handover preserves older labels in scripts: internal “Table S27” can mean printed Table S31, and “Fig. 11” can mean manuscript Fig. 10. Use the printed-number mapping in [the table catalogue](../supporting/05_Tables/README.md).
