# md/qm/gr/relax/: the relaxations of part B of gate G-R

Written on 11 Oct 2026. The method and the reading of these files are in `md/GR.md`, sections 1 and 6.6; this file
only says which record is which.

The ligand of the 64-atom zinc-site model (zinc, three ethyl-imidazoles, the sulfonamidate anion; vacuum; zinc and
imidazoles fixed) was relaxed at B3LYP-D3(BJ)/6-31G* from the crystal-like placement X and from the R-like placement
Rx, on the rented machine, by `md/qm/gr_queue.py relax`.

| File | What |
|---|---|
| `qm_X.json`, `qm_X.xyz` | The relaxation from X: every step, and the converged geometry (55 steps, 6 Oct 2026, 3:11 PM IST). |
| `qm_Rx.json` | The relaxation from Rx **as it stood at step 23** (6 Oct 2026, 3:45 PM IST; `state` "running", `crossed` false). This is the record that `md/GR.md` describes; it is kept as it was. |
| `qm_Rx_converged_step66.json` | The same relaxation **finished**: `state` "done", 66 steps, converged on 6 Oct 2026, 6:49:52 PM IST, `crossed` true. Added to git on 11 Oct 2026 as a byte copy of the record the watchdog brought from the machine (`/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/relax/qm_Rx.json` on the Pi). Steps 1 to 23 are those of `qm_Rx.json`. |
| `qm_Rx_converged_step66.xyz` | The **converged geometry** of that relaxation (step 66; its comment line gives E+D3 = −3584.6714605524 Hartree). Added to git on 11 Oct 2026 as a byte copy of `qm_Rx.xyz` of the machine. |
| `qm_Rx_step63_before_convergence.xyz` | The geometry of step 63 of the same relaxation, the last one the machine wrote as "not converged" (6:41 PM IST). A byte copy of `qm_Rx_last.xyz` of the machine, added on 11 Oct 2026. |
| `qm_Rx_last.xyz` | The geometry of step 23, the last step of `qm_Rx.json`. |
| `qm_Rx_before_crossing.xyz` | The geometry of step 15, written by the rule in force until 1:09 PM IST of 6 Oct 2026 and kept. |
| `mm_X.xyz`, `mm_Rx.xyz` | The force-field relaxations from the same two starts (`md/qm/gr/relax_mm.json`). |
| `relax_qm_both_converged.json` | The summary of the two relaxations written when the one from Rx had finished (6 Oct 2026, 6:49:52 PM IST; `both_converged` true). The file of the same kind one folder up, `md/qm/gr/relax_qm.json`, is the state of 3:48 PM IST and is kept. A byte copy of `relax_qm.json` of the machine, added on 11 Oct 2026. |
| `gr_gate_after_Rx_converged.json` | The record of gate G-R written by the machine at 8:47:30 PM IST of 6 Oct 2026, with part B read on both converged relaxations. `md/qm/gr/gr_gate.json` is the state of 3:48 PM IST, the one `md/GR.md` describes, and is kept. A byte copy of `gr_gate.json` of the machine, added on 11 Oct 2026. |
| `binding_check_final.json` | The binding check written at the same time (8:47:30 PM IST) from the formatted checkpoint of the finished optimisation: verdict "failed", no longer marked provisional. `md/qm/gr/binding_check.json` is the state of 3:48 PM IST ("failed (provisional)") and is kept. A byte copy of `binding_check.json` of the machine, added on 11 Oct 2026. |
| `fixed_atom_shift.json`, `qm_turns.json` | Measured on the Pi on the steps fetched from the machine up to 3:48 PM IST of 6 Oct 2026. |

## The finished relaxation from Rx

From `qm_Rx_converged_step66.json`, key `final`:

- signed volume at N2 −0.6783 at the start, +0.6192 at the end: the nitrogen has inverted and is pyramidal in the
  configuration of the crystal placement (`nitrogen_at_the_end`);
- heavy atoms of the ligand 0.1237 Å (r.m.s., no fitting) from the end point of the relaxation from X
  (`rmsd_ligand_heavy_from_X_run`), and 5.2597 Å from their start;
- energy −0.011 kcal mol⁻¹ against the end point of the relaxation from X (`e_total_minus_X_run_kcal_mol`);
- angles at N2 104.604, 125.547 and 109.857° (Zn–N2–S, Zn–N2–H, S–N2–H; sum 340.009°), Zn–N2–S–C4 −147.761°;
- last step: r.m.s. gradient on the ligand 1.3e-04 Hartree bohr⁻¹; fixed atoms: `largest_move_of_a_fixed_atom_A`
  5e-08 Å, and `largest_shift_of_a_fixed_atom_over_the_run_A` 0.0047 Å, the value the record of step 23 already had.

In the bare zinc-site model in vacuum the R-like placement therefore ends in the minimum reached from X. The model has
no Thr199, no pocket wall and no water, so this says nothing on what the ligand does in the protein. The same numbers
are in `md/qm/n2/final_test/behaviour.json` (`along_the_quantum_chemical_relaxations.runs.Rx`) and
`md/qm/n2/n2_gate_result.json`, which were computed from the finished run and are the sources of the paper's facts
`md.qm.relax.*`.

## The later records of the machine (added on 11 Oct 2026, 10:25 AM IST)

`md/qm/gr/relax_qm.json`, `md/qm/gr/gr_gate.json` and `md/qm/gr/binding_check.json` in git are those of 3:48 PM IST of
6 Oct 2026 (`both_converged` false; binding check "failed (provisional)"). The versions the machine wrote after the
relaxation from Rx had finished existed only on the Pi, in `/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/`, where the
watchdog had put them. On 11 Oct 2026 they were copied into this folder as byte copies under names of their own, so
that no tracked file is overwritten: the five files marked "added on 11 Oct 2026" in the table above, beside
`qm_Rx_converged_step66.json`.

What the later records hold beyond the tracked ones:

- `relax_qm_both_converged.json`: `both_converged` true; `Rx_minus_X_with_d3` −0.011 and `Rx_minus_X_without_d3`
  −0.05 kcal mol⁻¹; `end_points.rmsd_ligand_heavy_A` 0.1237 and `end_points.largest_distance_of_a_ligand_heavy_atom_A`
  0.2204 Å; `starts.Rx.n_steps` 66, `starts.Rx.state` "done".
- `gr_gate_after_Rx_converged.json`: the same end points under `with_d3.part_B.relaxed.quantum_chemistry` and
  `without_d3.part_B.relaxed.quantum_chemistry`. The force-field descents were run again on the machine and differ
  from the tracked record in the last digits (for example `careful_descent.Rx.energy_drop` 24.2023 against 24.2027).
- `binding_check_final.json`: `verdict` "failed" and `provisional` false; the source of the quantum-chemical optimum
  is the formatted checkpoint and no longer the array of the finished optimisation.

Not copied, because the tracked file is the later one or the same: `qm_X.json` and `qm_X_last.xyz` of the machine
(3:03 PM IST; git has `qm_X.json` and `qm_X.xyz` of 3:11 PM IST), `mm_X.xyz`, `mm_Rx.xyz`, `relax_mm.json`, the two
`gr_*_gate.json` variants, the three `mm_2e_*.json` files and `fit/` (all of the morning of 6 Oct 2026 on the machine,
of the afternoon in git), and `gr_table.csv`, `points.csv`, `points.json` and the two scan tables (byte-identical).

The facts of the paper `md.qm.relax.*` are read from `md/qm/n2/final_test/behaviour.json` and
`md/qm/n2/n2_gate_result.json`, not from the files of this folder; the copies are the records of the same run and
give the same numbers (0.1237 Å, −0.011 kcal mol⁻¹, signed volume +0.6192).

## Checksums (sha256)

```
35bcb734564f7ca67fe18950434fa0d61645af77d9396d44e735a298fbcc4a2c  qm_Rx_converged_step66.json
1aa946bbc6d8f7843c6a1b4007bc693365f642a1c8276f2dd21b165fd1d2a38c  qm_Rx_converged_step66.xyz
f8e8ca8d0e8883904addeecccf179fb92cb4f1ecae0c050250a86a0d69b4425a  qm_Rx_step63_before_convergence.xyz
403a39577181991af6deb1d54918e1b338b005d6a8d4affd2b5a63688d377b12  relax_qm_both_converged.json
ed6a8b450e714178c3d1957bdcf7f176993a33ccac76c2b17c931ca483ea59f7  gr_gate_after_Rx_converged.json
f8f763b5ac2daf875b8bff5c6cd5f127dd21f7214e9dbb84c2258b4f523015b4  binding_check_final.json
aed93649ee970e9c16930806e2f853ea65e3276f1233ec39aeabc42c35c8aa8b  qm_Rx.json
```

The first six are byte copies of files of `/root/cqc-sfa-md/from_gpu/repo_md/qm/gr/` on the Pi (`relax/qm_Rx.json`,
`relax/qm_Rx.xyz`, `relax/qm_Rx_last.xyz`, `relax_qm.json`, `gr_gate.json`, `binding_check.json`), compared with
`cmp` after the copy; the last is the tracked record of step 23. Check: `sha256sum -c` on the block above, run in
this folder.
