# Second analysis against the first: /root/cqc-sfa-md/runs

Written by `md/analysis/recheck_md.py` (sha256 `5b4c0595c067ff57…`); first program `analyse_md.py` (sha256 `f75e6155de700700…`). Specification: md/analysis/SPEC.md version 1.3. Tolerances: section 14 (distances 0.001 Å, angles 0.01°, signed volume 0.001, times 0.000001 ns, counts, flags, states and verdicts identical).

**Result: the two programs agree on all 15 run(s) and on the campaign rules.**

| Run | Solute / full frames | Values compared | Disagreements | Threshold frames and what follows | Keys in one program only | Warnings (first / second) | Agree |
|---|---|---|---|---|---|---|---|
| bonded_2e/N1 | 10000 / 1000 | 602070 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/N2 | 10000 / 1000 | 602065 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/N3 | 10000 / 1000 | 602065 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/R1 | 10000 / 1000 | 605010 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/R2 | 10000 / 1000 | 602084 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/R3 | 10000 / 1000 | 602084 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/X0 | 500 / 50 | 42845 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/X1 | 10000 / 1000 | 602070 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/X2 | 10000 / 1000 | 602065 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e/X3 | 10000 / 1000 | 602065 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2e_n2fit/R1 | 10000 / 1000 | 602383 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2s/R1 | 10000 / 1000 | 602103 | 0 | 0 | 0 | 0 / 0 | yes |
| bonded_2s/X0 | 500 / 50 | 42850 | 0 | 0 | 0 | 0 / 0 | yes |
| preflight_a | 40 / 4 | 4622 | 0 | 0 | 0 | 1 / 1 | yes |
| preflight_b | 40 / 4 | 4627 | 0 | 0 | 0 | 1 / 1 | yes |

Campaign rules (`rules_per_model`, folders left out, order of the runs): 65 values compared, 0 disagreements; agree.

Largest differences found, over all values compared: series ang 0; series dist 0.0001; series time_ns 0; series tors 0; series vol 0; summary ang 0; summary dist 9e-06; summary fraction 0; summary number 0; summary time_ns 0; summary time_ps 0; summary tors 0; summary vol 0.

No disagreement to list.

Texts that the specification does not fix and that differ (not counted): `kind` 'A' against 'dist' (465 times); `kind` 'deg' against 'ang' (210 times); `n2_configuration/order` 'e(Zn) . (e(HN21) x e(S)), unit vectors from N2' against 'Zn, HN21, S' (15 times).

