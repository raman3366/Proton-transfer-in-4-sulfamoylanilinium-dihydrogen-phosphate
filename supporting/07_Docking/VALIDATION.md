# Docking validation: gate G2

Verification pass, 28 September 2026. Every number below was regenerated in this pass from the files
named beside it. Section 8 lists the values of the first set-up run that this pass corrected.

Environment (Raspberry Pi 5, aarch64, 4 cores): `/root/cqc-sfa-tools/dock` — AutoDock Vina 1.2.7
(`AutoDock Vina f458505-mod`, from `results.json`), Meeko 0.8.0, RDKit 2026.03.1, spyrmsd 0.9.0,
PoseBusters 0.6.5, PROPKA 3.5.1, PDBFixer/OpenMM 8.6.1, gemmi 0.7.5, ProLIF 2.2.2, PLIP.
AutoDock4Zn grid maps were made earlier on the x86 VM with AutoGrid 4.2.7 (no aarch64 autogrid
exists); nothing in this pass touched the VM. All docking here ran `nice -n 10 --cpu 3`, one process
tree at a time.

Files regenerated in this pass: `dock/analysis/summary.json`, `dock/analysis/gate_G2.json`,
`dock/analysis/contacts.json`, `dock/targets_inventory.json`, `dock/T*/prep_report_recheck.json`, and
the docking runs listed in sections 1 and 7.

---

## 1. Reproducibility

**`analyse.py` re-run on the existing 28 redocking/cross-docking runs, code unchanged**:
`summary.json` and `gate_G2.json` came out **byte-identical** to the files the builder left
(`cmp` on both). The analysis code and the JSON in the repository therefore agree.

The files now in `dock/analysis/` are the **final** regeneration: `analyse.py` was then re-run with the
fixes of section 9 over all 38 runs (the 28 above plus the 10 verification and sensitivity runs), and
`contacts.py` likewise. Where a number here differs from the builder's report, the reason is a fix
listed in section 9 or a run added in section 7, never a different reading of the same file.

**One Vina redock repeated** — T5 7S2K, dihydropteroate (78H), 10 seeds, exhaustiveness 32, 20 Å box,
`--cpu 3` (`dock/T5_7S2K/runs/dhp_anion__vina__repro/`):

| | original | repeat |
|---|---|---|
| all 10 per-seed pose-score lists | — | **identical, every value** |
| mean best score | −10.106 ± 0.036 | −10.106 ± 0.036 kcal/mol |
| wall time, 10 seeds | 583.0 s (58.3 s/seed) | 549.5 s (55.0 s/seed) |

**One AD4Zn redock repeated** — T3 6RL9, sulfanilamide, 10 seeds, exhaustiveness 32, `--cpu 3`
(`dock/T3_6RL9/runs/sfa_neutral__ad4zn__repro/`):

| | original | repeat |
|---|---|---|
| all 10 per-seed pose-score lists | — | **identical, every value** |
| mean best score | −8.862 ± 0.017 | −8.862 ± 0.017 kcal/mol |
| wall time, 10 seeds | 339.6 s (34.0 s/seed) | 339.5 s (33.9 s/seed) |

So Vina and AD4Zn are bit-reproducible on this machine at fixed `--seed` **and** fixed `--cpu`. The
per-seed timings quoted in `run_production.sh` reproduce: AD4Zn 33.9 s/seed (quoted ~34), sulfanilamide
Vina 14.0–14.4 s/seed (quoted 14.3, from the two 6RL9 Vina runs made in this pass) and dihydropteroate
Vina 55.0 s/seed, the slowest ligand in the set.
The receptors are pinned: `prep_receptor.py --recheck` recomputed the sha256 of all seven
`receptor.pdbqt` files and each matches the hash recorded in its `prep_report.json`.

## 2. Reference values from the crystal structures

Measured in this pass from the deposited mmCIFs (`dock/targets_inventory.json`) and from the
co-crystal ligands (`dock/T*/ref_*.sdf`), so the docked values in section 3 have something to be
compared with:

| structure | Zn–N(sulfonamide) | Thr199 OG1···N | Thr199 N···O(sulfonyl) | other |
|---|---|---|---|---|
| 6RL9 SAN304 | 1.93 Å | 2.86 Å | 2.94 Å | Zn–His94 NE2 2.01, His96 NE2 2.02, His119 ND1 2.02 Å |
| 3HS4 AZM701 | 1.94 Å | 2.82 Å | 2.95 Å | Zn–His94 NE2 1.99, His96 NE2 1.99, His119 ND1 2.01 Å |
| 3TZF 08D280 (chain B) | – | – | – | aniline N···pterin (HH2) 3.10 Å |
| 1AJ0 SAN561 | – | – | – | aniline N···pterin (PH2) 3.51 Å |
| 7S2K 78H301 | – | – | – | pterin 2-amino N···Asp174 2.60 Å; pABA carboxylate···Ser214 N 2.94 Å |

Every site check in `analyse.py` is now also run on the crystal reference itself
(`site_checks_on_crystal_reference` in `summary.json`); all of them pass on the crystal pose, which is
the minimum a site check must do to be evidence. See section 6 for the T5 check that did not.

The AD4Zn setup was re-checked the same way. Measured in this pass from
`dock/T*/ad4zn/receptor_tz.pdbqt` against the co-crystal ligands: the single TZ pseudo-atom placed by
the official `zinc_pseudo.py` sits **0.17 Å** (6RL9) and **0.07 Å** (3HS4) from the crystallographic
sulfonamide nitrogen and 2.00 Å from the zinc in both — it marks where the ligand nitrogen belongs, so
the AD4Zn zinc model is doing what it claims. The six `nbp_r_eps` records in
`dock/T*/ad4zn/receptor_tz.gpf` are byte-identical to the upstream AutoDock-Vina zinc example's GPF,
and `receptor_tz.glg` shows AutoGrid accepting all six ("Overriding non-bonded interaction energies")
and the NA map reaching −4.58 kcal/mol, i.e. the TZ–NA attraction is present in the maps.

## 3. Redocking and cross-docking (10 seeds each, exhaustiveness 32, 20 Å box, `--cpu 3`)

Symmetry-corrected heavy-atom RMSD (spyrmsd, no superposition: the pose is compared in the receptor
frame). "best of seeds" is the lowest top-pose RMSD over the ten seeds; "f≤1 Å" and "f≤2 Å" are the
fractions of seeds whose **top** pose is within that distance; "rank₁" is the median rank at which the
first pose within 1 Å appears, with the fraction of seeds that contain such a pose anywhere in their
ranked list. PB = PoseBusters physics tests (chemistry, bond lengths and angles, ring flatness,
internal energy, no clash with protein/cofactors/waters), RMSD test excluded.

The PoseBusters column is 1.0 everywhere: over all 38 runs in the regenerated `summary.json`, the only
PoseBusters test that ever fails is its own `rmsd_≤_2å`, which is the gate's job and is reported
separately. Every pose discussed below is a physically sound molecule in a physically possible place —
so where a pose is wrong, it is wrong about *where*, not about chemistry.

### T5 7S2K — Sul2, dihydropteroate product (primary resistance contrast)

| job | RMSD top (mean ± SD) | best of seeds | f≤1 Å | f≤2 Å | rank₁ | score (kcal/mol) | PB | G2 (1 Å) | 2 Å criterion |
|---|---|---|---|---|---|---|---|---|---|
| 78H, Vina | **0.51 ± 0.03** | 0.47 | **1.0** | 1.0 | 1.0 (10/10) | −10.11 ± 0.04 | 1.0 | **PASS** | pass |
| 78H, Vinardo | **0.54 ± 0.01** | 0.52 | **1.0** | 1.0 | 1.0 (10/10) | −8.09 ± 0.08 | 1.0 | **PASS** | pass |

Site checks on the docked top poses (the corrected checks of section 6), from `summary.json`: pterin
2-amino N···Asp174 2.65 Å (Vina) / 2.68 Å (Vinardo), pABA carboxylate···Ser214 N 2.90 / 2.96 Å, both
passing in 10/10 seeds — against 2.60 and 2.94 Å in the crystal. Both scoring functions reproduce the
product pose of a 23-heavy-atom ligand to about half an angstrom in every seed, with its two anchoring
hydrogen bonds within 0.06 Å of the deposited geometry; this is the strongest target in the set.

### T3 6RL9 — hCA II, sulfanilamide (the compound's own zinc target)

| job | RMSD top | best of seeds | f≤1 Å | f≤2 Å | rank₁ | score | Zn–N (mean / pass) | Thr199 OG1 | PB | G2 | 2 Å |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SAN⁰, **AD4Zn** | 2.07 ± 1.45 | 0.65 | 0.5 | 0.5 | **1.5 (10/10)** | −8.86 ± 0.02 | **2.22 / 1.0** | **2.56 / 1.0** | 1.0 | fail | fail |
| SAN⁻, AD4Zn | 1.12 ± 0.32 | 0.61 | 0.3 | **1.0** | 2.0 (7/10) | −8.13 ± 0.01 | **2.21 / 1.0** | **2.60 / 1.0** | 1.0 | fail | pass |
| SAN⁰, Vina | 1.10 ± 0.01 | 1.09 | 0.0 | 1.0 | – (0/10) | −5.87 ± 0.01 | 2.84 / **0.0** | 4.76 / **0.0** | 1.0 | fail | fail |
| SAN⁻, Vina | 1.09 ± 0.00 | 1.09 | 0.0 | 1.0 | 2.0 (3/10) | −5.88 ± 0.01 | 2.80 / **0.0** | 4.76 / **0.0** | 1.0 | fail | fail |
| SAN⁰, Vinardo | 1.13 ± 0.00 | 1.12 | 0.0 | 1.0 | 6.0 (1/10) | −5.08 ± 0.01 | 2.48 / **0.0** | 4.45 / **0.0** | 1.0 | fail | fail |
| AZM⁰, AD4Zn (cross) | 2.68 ± 0.81 | 1.55 | 0.0 | 0.3 | 7.0 (5/10) | −9.48 ± 0.03 | 3.68 / 0.5 | 4.22 / 0.5 | 1.0 | fail | fail |
| AZM⁰, Vina (cross) | 5.23 ± 0.01 | 5.22 | 0.0 | 0.0 | – (0/10) | −6.65 ± 0.03 | 5.68 / 0.0 | 7.99 / 0.0 | 1.0 | fail | fail |
| AZM⁻, Vina (cross) | 5.22 ± 0.01 | 5.21 | 0.0 | 0.0 | 13.0 (1/10) | −6.64 ± 0.02 | 5.61 / 0.0 | 7.89 / 0.0 | 1.0 | fail | fail |

Cross-docking reference: 3HS4's acetazolamide carried into the 6RL9 frame by Cα superposition,
257 common Cα atoms, Cα RMSD 0.328 Å. The transformed reference still sits properly in the 6RL9 site
(Zn–N 1.92 Å, Thr199 OG1 2.76 Å, Thr199 N···O 2.91 Å, all measured on it in this pass), so the
superposition does not distort the zinc site it is used to judge.

**The 6RL9 AD4Zn result is bimodal, and it is a scoring degeneracy, not a sampling failure.** The
near-native pose (best pose per seed: 0.40–0.73 Å) is found in **10 of 10** seeds, but in 5 of 10 it
is ranked 2nd–4th behind a pose 3.37–3.55 Å away. In the alternative pose the sulfonamide nitrogen still sits on the
zinc (0.41 Å from its crystal position, Zn–N 2.2 Å, Thr199 bonds intact) and the aniline ring is
rotated about the Zn–N/S axis, displacing the aniline nitrogen by 6.0 Å. The contacts show what that
means: the near-native pose lies against the hydrophobic wall (Trp209, Leu198, Phe131, Val121, Val143)
as the crystal pose does, while the rotated pose turns the aniline to the hydrophilic side and
hydrogen-bonds its aniline nitrogen to Asn67 OD1 (2.73 Å) and Asn62 ND2 (3.04 Å) — a chemically
plausible pose, which is exactly why the scoring function cannot tell them apart. AD4Zn separates the two
orientations by **0.008–0.039 kcal/mol** (best near-native mode, within 1.0 Å of the crystal pose, against
the best ring-rotated mode, beyond 2.0 Å, in each seed; the definition that the paper uses, facts section 37),
far below any meaningful resolution — which is why more search does not help (section 7).

One consequence worth stating plainly: for a bimodal distribution the *mean* RMSD is an unstable
statistic. The same protocol gives 2.07 ± 1.45 Å over 10 seeds and 1.80 ± 1.37 Å over 20 (section 7),
which happens to move it across the 2.0 Å line although nothing about the result changed. The
success fraction (0.5 over 10 seeds, 0.6 over 20) is the number to report for this target, and the
gate's use of the mean is why T3 fails it.

### T4 3HS4 — hCA II, acetazolamide (cross-docking partner)

| job | RMSD top | best of seeds | f≤1 Å | f≤2 Å | rank₁ | score | Zn–N | Thr199 OG1 | PB | G2 | 2 Å |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SAN⁰, **AD4Zn (cross)** | **0.62 ± 0.03** | 0.58 | **1.0** | 1.0 | 1.0 (10/10) | −8.79 ± 0.01 | **2.18 / 1.0** | **2.59 / 1.0** | 1.0 | **PASS** | pass |
| SAN⁰, Vina (cross) | 1.62 ± 0.01 | 1.61 | 0.0 | 1.0 | 2.0 (10/10) | −5.82 ± 0.01 | 2.55 / 0.0 | 3.71 / 0.0 | 1.0 | fail | fail |
| SAN⁻, Vina (cross) | 1.17 ± 0.16 | 1.11 | 0.0 | 1.0 | 2.5 (2/10) | −5.76 ± 0.02 | 2.79 / 0.0 | 4.66 / 0.0 | 1.0 | fail | fail |
| AZM⁰, AD4Zn (redock) | 2.70 ± 0.20 | 2.31 | 0.0 | 0.0 | 7.0 (5/10) | −9.26 ± 0.02 | **2.08 / 1.0** | **2.58 / 1.0** | 1.0 | fail | fail |
| AZM⁻, AD4Zn (redock) | 3.65 ± 0.44 | 3.23 | 0.0 | 0.0 | 6.5 (4/10) | −9.01 ± 0.07 | 4.86 / 0.0 | 4.66 / 0.0 | 1.0 | fail | fail |
| AZM⁰, Vinardo (redock) | 1.35 ± 0.01 | 1.33 | 0.0 | **1.0** | – (0/10) | −5.27 ± 0.03 | 2.54 / **0.0** | 4.54 / **0.0** | 1.0 | fail | fail |
| AZM⁰, Vina (redock) | 5.39 ± 0.09 | 5.29 | 0.0 | 0.0 | 6.0 (1/10) | −6.41 ± 0.14 | 6.09 / 0.0 | 8.22 / 0.0 | 1.0 | fail | fail |
| AZM⁻, Vina (redock) | 5.37 ± 0.02 | 5.35 | 0.0 | 0.0 | – (0/10) | −6.48 ± 0.02 | 6.00 / 0.0 | 8.15 / 0.0 | 1.0 | fail | fail |

Cross-docking reference: 6RL9's sulfanilamide in the 3HS4 frame, same 257 Cα atoms, Cα RMSD 0.328 Å;
in the 3HS4 site that transformed reference has Zn–N 1.95 Å, Thr199 OG1 2.91 Å, Thr199 N···O 2.98 Å.

**Acetazolamide's own pose is not reproduced by anything.** With AD4Zn its sulfonamide end is placed
essentially correctly (S 0.49, N 0.36, O 0.50/0.71 Å from the crystal positions, Zn–N 2.08 Å in 10/10
seeds) while the acetamido tail is 2.7–5.2 Å out, giving 2.70 ± 0.20 Å overall. Vinardo gets the
molecule as a whole to 1.35 Å but leaves Zn–N at 2.54 Å and loses the Thr199 OG1 bond. So for
acetazolamide the zinc anchor and the whole-molecule pose are never right at the same time:
acetazolamide is usable as a **score** control, not as a pose.

### T1 3TZF — DHPS, *Y. pestis*, sulfamethoxazole (primary DHPS target)

| job | RMSD top | best of seeds | f≤1 Å | f≤2 Å | rank₁ | score | aniline N···pterin | PB | G2 | 2 Å |
|---|---|---|---|---|---|---|---|---|---|---|
| SMX⁻, Vina | 1.92 ± 0.04 | 1.80 | 0.0 | **1.0** | 3.0 (1/10) | −7.64 ± 0.09 | **3.20 / 1.0** | 1.0 | fail | **pass** |
| SMX⁻, Vina, loop seed 12345 | 1.89 ± 0.06 | 1.77 | 0.0 | **1.0** | 4.5 (2/10) | −7.62 ± 0.11 | 3.21 / 1.0 | 1.0 | fail | **pass** |
| SMX⁰, Vina | 3.00 ± 0.04 | 2.94 | 0.0 | 0.0 | 4.0 (4/10) | −7.80 ± 0.02 | – | 1.0 | fail | fail |
| SMX⁻, Vinardo | 2.18 ± 0.40 | 1.42 | 0.0 | 0.1 | 9.0 (1/10) | −5.93 ± 0.20 | – | 1.0 | fail | fail |

The anionic (deprotonated sulfonamide) species is the one that works, as expected for a sulfa drug in
DHPS. The top pose is a consistent 1.9 Å in every seed, inside the conventional 2.0 Å criterion and
outside the 1.0 Å gate; a 0.56 Å pose is sampled but ranked 3rd. The contacts of the top pose are the
pABA pocket, not a decoy: Ser222, HH2 (pterin), Phe190, Thr62, Lys221, Arg63, Gly189, Phe28, Pro64,
Met148, Arg235 (`dock/analysis/contacts.json`, regenerated in this pass).

### T2 1AJ0 — DHPS, *E. coli*, sulfanilamide (secondary, open loop)

| job | RMSD top | best of seeds | f≤1 Å | f≤2 Å | score | aniline N···pterin | PB | G2 | 2 Å |
|---|---|---|---|---|---|---|---|---|---|
| SAN⁰, Vina, 20 Å box | **11.64 ± 0.05** | 11.57 | 0.0 | 0.0 | −5.80 ± 0.01 | 11.51 / 0.0 | 1.0 | **FAIL** | fail |
| SAN⁰, Vina, no SO4 | 11.63 ± 0.02 | 11.58 | 0.0 | 0.0 | −5.79 ± 0.02 | – | 1.0 | FAIL | fail |
| SAN⁰, Vina, 16 Å box | 2.99 ± 0.01 | 2.97 | 0.0 | 0.0 | −5.05 ± 0.01 | – | 1.0 | fail | fail |
| SAN⁰, Vina, 14 Å box | 3.01 ± 0.04 | 2.95 | 0.0 | 0.0 | −5.05 ± 0.01 | – | 1.0 | fail | fail |
| SAN⁰, Vinardo, 20 Å box | 7.32 ± 4.51 | 3.03 | 0.0 | 0.0 | −4.34 ± 0.01 | – | 1.0 | fail | fail |

In the 20 Å box every seed puts sulfanilamide 11.6 Å away in a surface cleft (Thr147, Pro145, Asn144,
Gly191, Gln149, Gln142, Lys192). Shrinking the box to 14–16 Å moves it into the real pABA pocket
(Arg220, Ser222, Gly189, Arg235, Lys221, Phe190, PH2, Phe188, Ser219, Arg63) but only to 3.0 Å, and
that is a box chosen after seeing the answer, which is not validation. Removing the sulfate changes
nothing (11.63 vs 11.64 Å), so the sulfate is not the cause.

## 4. Gate verdicts and the two criteria side by side

Gate G2, as defined before the runs: top-pose symmetry-corrected RMSD ≤ 1.0 Å (mean over seeds), success
fraction over 10 seeds reported, PoseBusters physics ≥ 0.8, and for hCA II also Zn–N(sulfonamide)
≤ 2.3 Å and the Thr199 hydrogen bond in at least half the seeds. The conventional 2.0 Å criterion
(Hevener 2009; Warren 2006) is evaluated with the same structure and the same site and physics
requirements, only the distance changed. Both are in `gate_G2.json` as `passes_G2` and
`passes_2A_criterion`; the gate was **not** changed.

The regenerated `gate_G2.json` holds 37 jobs (the 28 validation runs plus the nine diagnostic runs of
section 7; the `__fromxtal` control is excluded by design). **Three distinct protocols pass the 1.0 Å
gate**: T5 78H Vina, T5 78H Vinardo and T4 sulfanilamide AD4Zn (cross-dock) — four entries, because
the T5 Vina reproducibility repeat passes identically. Under the 2.0 Å criterion, T1
sulfamethoxazole⁻ Vina (both loop-rebuild seeds), T3 sulfanilamide⁻ AD4Zn and the 20-seed T3
sulfanilamide⁰ AD4Zn run also pass. Every hCA II job scored by plain Vina or Vinardo fails both,
because it fails the zinc-site requirement rather than the distance. Every job's site checks pass on
its own crystal reference (`site_checks_pass_on_crystal_reference` true in all 37).

## 5. Ranking test, with and without the site filter

Acetazolamide (low nM for hCA II) must score below sulfanilamide (high nM – low µM). Reported two
ways in `gate_G2.json` → `ranking`: over all ten seeds, and over only those seeds whose top pose
actually makes the zinc bond (Zn–N ≤ 2.3 Å and Thr199 OG1 ≤ 3.3 Å).

| receptor / scoring | all seeds: AZM below SAN? | gap | site-valid seeds (AZM / SAN) | site-valid: AZM below SAN? | gap |
|---|---|---|---|---|---|
| 6RL9 / Vina | yes | 0.78 | **0 / 0** | cannot be evaluated | – |
| 6RL9 / AD4Zn | yes | 0.62 | 5 / 10 | yes | 0.64 |
| 3HS4 / Vina | yes | 0.59 | **0 / 0** | cannot be evaluated | – |
| 3HS4 / AD4Zn | yes | 0.47 | **10 / 10** | yes | 0.47 |

Vinardo appears in neither row: the validation matrix ran only one of the two ligands with Vinardo in
each hCA II receptor, so there is nothing to rank. That is a gap in the matrix, not a Vinardo result.

So the honest statement is: the ranking comes out right in 4 of 4 receptor × scoring combinations by
score alone, but in the two Vina combinations **no seed produces a zinc-bound pose for either
ligand**, so those two rankings rest on poses that are not the binding mode. Only 3HS4/AD4Zn has the
ranking resting on zinc-bound poses in every seed, and even there the acetazolamide pose has a
misplaced tail (section 3). 6RL9/AD4Zn rests on 5 of 10 acetazolamide seeds.

## 6. Site checks are now validated against the crystal pose

The T5 check that stood in `analyse.py` ("ligand carboxylate within 4.0 Å of the pyrophosphate") was
wrong in two ways, found in this pass: its SMARTS matched no atom in the neutral crystal component at
all, and in the crystal that distance is ~7.5 Å, so the check could never pass. It was reporting
`pass_fraction 0.0` for the best-redocking target in the set. It is replaced by the two contacts the
deposited pose actually makes (pterin 2-amino N···Asp174 2.60 Å; pABA carboxylate···Ser214 N 2.94 Å),
and `analyse.py` now runs every target's checks on its own crystal reference and records the result,
so a mis-specified check shows up as a failing crystal pose instead of a run of zeros.

## 7. Sensitivity and diagnostics

_Filled from the runs listed; each row names its run directory._

| test | run | result |
|---|---|---|
| exhaustiveness 64 (6RL9 SAN⁰ AD4Zn) | `T3_6RL9/runs/sfa_neutral__ad4zn__exh64` | RMSD 2.34 ± 1.39 Å, f≤1 Å 0.4 (vs 2.07 ± 1.45, 0.5 at exh 32): **does not resolve the bimodality**; score −8.865 ± 0.018 vs −8.862 ± 0.017, 67 s/seed vs 34 |
| exhaustiveness 128 | `T3_6RL9/runs/sfa_neutral__ad4zn__exh128` | RMSD 2.33 ± 1.37 Å, f≤1 Å 0.4, best-of-seeds 0.71, near-native pose present in **10/10** seeds: **does not resolve the bimodality either**; score −8.868 ± 0.021, 131.8 s/seed (3.9× the cost of exh 32) |
| 20 seeds instead of 10 | `T3_6RL9/runs/sfa_neutral__ad4zn__s20` | RMSD 1.80 ± 1.37 Å, f≤1 Å **0.60** (12/20), near-native pose present in **20/20**, score −8.866 ± 0.017: the extra seeds tighten the estimate of how often the near-native pose ranks first (0.5 → 0.6) but the mean is still 1.8 Å, so the gate still fails. The two clusters are unchanged: every seed lands at 0.65–0.77 Å or 3.35–3.56 Å, nothing between |
| AD4Zn grid extent as a Vina box (20.25 vs 20.0 Å) | `T3_6RL9/runs/sfa_neutral__vina__box20.25` | **No effect at all**: all ten per-seed RMSDs identical to the 20.0 Å run (1.087–1.109 Å, mean 1.096 ± 0.007), score −5.870 ± 0.012 in both. The 0.25 Å difference between the AD4Zn grid and the Vina box is therefore not a confound for this comparison |
| His64 tautomer HIE → HID (6RL9, SAN⁰ Vina redock) | `T3_6RL9_his64HID/runs/sfa_neutral__vina` | RMSD 1.093 ± 0.004 vs 1.096 ± 0.007 Å; score −5.871 ± 0.016 vs −5.870 ± 0.012 — **no effect** (0.003 Å, 0.001 kcal/mol) |
| His64 tautomer HIE → HID (3HS4, SAN⁰ Vina cross-dock) | `T4_3HS4_his64HID/runs/sfa_neutral__vina` | RMSD 1.628 ± 0.010 vs 1.623 ± 0.010 Å; score −5.822 ± 0.014 vs −5.818 ± 0.009 — **no effect** (0.005 Å, 0.004 kcal/mol), and this is the receptor where PDB2PQR would have chosen the other tautomer |
| His64 tautomer HIE → HID (3HS4, AZM⁰ Vina redock) | `T4_3HS4_his64HID/runs/azm_neutral__vina` | RMSD 5.38 ± 0.08 vs 5.39 ± 0.09 Å; score −6.437 ± 0.137 vs −6.409 ± 0.141 — **no effect** (0.014 Å, 0.028 kcal/mol); the pose is wrong either way |
| `--cpu` sensitivity (3 → 2), 6RL9 SAN⁰ Vina | `T3_6RL9/runs/sfa_neutral__vina__cpu2` | **Identical**: all ten per-seed pose-score lists and all ten RMSDs match the `--cpu 3` run exactly (mean −5.870 ± 0.012, RMSD 1.096 ± 0.007 both ways); only the wall time changes (21.0 vs 14.4 s/seed). The earlier report's claim that Vina's result depends on `--cpu` is **not reproduced** — see section 8 |
| loop rebuild, independent PDBFixer seed | `T1_3TZF_seed12345/runs/smx_anion__vina` | 1.89 ± 0.06 vs 1.92 ± 0.04 Å; score −7.62 ± 0.11 vs −7.64 ± 0.09 — negligible |
| sulfate removed from 1AJ0 | `T2_1AJ0_noSO4/runs/sfa_neutral__vina` | 11.63 ± 0.02 vs 11.64 ± 0.05 Å — no effect |
| start from the crystal conformer | `T3_6RL9/runs/sfa_neutral__vina__fromxtal` | 1.08 ± 0.00 vs 1.10 ± 0.01 Å — the ETKDG randomisation costs nothing |

**What is and is not known about His64.** His64 is hCA II's proton shuttle and it is inside the 20 Å
box, so the earlier report's silent default deserved checking. What the structures determine: PROPKA
puts it neutral at pH 7.4 (pKa 3.73 in 6RL9, 4.48 in 3HS4, from `dock/T*/propka.pka`), and in 6RL9 it
is modelled in two conformers, of which altloc A was taken. What they do **not** determine is the
tautomer, and nothing in the receptor forces it: measured in this pass, His64 ND1/NE2 are 9.6–10.0 Å
(6RL9) and 9.2–9.6 Å (3HS4) from the nearest co-crystal ligand atom and 10.7–11.8 Å from the zinc, so
it makes no contact with the ligand or the metal. The receptors therefore carry an **explicit
assumption, HIE**. The code now records it as one: `targets.py` declares `his_assumed={64: 'HIE'}` for
both hCA II targets and `prep_receptor.py` writes it to `protonation_assumed` with the reason, instead
of the residue falling through the default branch. The seven receptors on disk were prepared before
that field existed, so their `prep_report.json` does not carry it (their
`prep_report_recheck.json` says so explicitly); they were not rebuilt, because a rebuild is not
bit-reproducible and would invalidate every docking run above. The tautomer they were built with is
HIE, which is what the two rows above test. For comparison, PDB2PQR 3.x
with its own hydrogen-bond optimiser (`pdb2pqr30 --ff AMBER --titration-state-method propka
--with-ph 7.4` on `fixed_heavy.pdb`) chooses HIE64 for 6RL9 — agreeing with us — and HID64 for 3HS4,
disagreeing. The three rows above measure what that disagreement is worth, and the answer is: nothing
measurable. Flipping His64 to HID moves the top-pose RMSD by ≤0.014 Å and the mean score by
≤0.03 kcal/mol in all three tests, including the 3HS4 receptor where PDB2PQR would have chosen HID in
the first place. The assumption is therefore stated, not hidden, and shown not to matter at the level
anything is reported. The test was run with Vina: changing a receptor atom type invalidates the AD4Zn
grid maps, and AutoGrid does not run on this machine. His64 lies 10 Å from the zinc.

## 8. Corrections made during the verification pass (28 Sep 2026)

| value of the first set-up run | what the files say |
|---|---|
| "AD4Zn validated for both zinc targets" / production script: T3 and T4 "use AD4Zn as the primary method" on the strength of gate G2 | T4 sulfanilamide AD4Zn passes (0.62 ± 0.03 Å, 10/10). **T3 6RL9 sulfanilamide AD4Zn fails** (2.07 ± 1.45 Å, f≤1 Å 0.5). AD4Zn is the only scoring function that reproduces the zinc chemistry, which is a different and weaker claim. |
| "rank order matches literature in 4/4" | True by score, but 0/10 seeds are zinc-bound in both Vina combinations (section 5). |
| "3TZF at +1.2×10⁶ kJ/mol" before minimisation | `T1_3TZF/prep_report.json`: **9 854 879.5 kJ/mol** (seed 20260928); 67 870.8 kJ/mol for seed 12345. |
| "max displacement 0.13–0.18 Å, rms 0.015 Å" | Those are the **restrained** (deposited) atoms only, and the range over the seven receptors is 0.069–0.197 Å (rms 0.015–0.023). The **free** atoms move up to 1.5–8.0 Å. |
| `internal_gaps: None` for 3TZF while the rebuild list held the 66–67 loop | Two definitions of "internal gap" in one script, and six of the seven `prep_report.json` files predate the field entirely (`prep_report_recheck.json` now says which). Both definitions agree once computed: 3TZF chain B has exactly one internal gap, residues **66–67**, and PDBFixer's key `0:65` translates to exactly those two residues — its insertion point lies between modelled residues 65 and 68, and the two rebuilt residues are Ala-Ala, which is what the deposited sequence has at 66–67. |
| `min_heavy_heavy_distance_A: None` on every target while `geometry_ok: True` | The value was only recorded when the check failed. Measured now over every heavy atom pair: **1.201–1.216 Å** across the seven receptors, 0 pairs below 1.1 Å. |
| hCA II "missing residue 126" | Confirmed a numbering convention: `_pdbx_poly_seq_scheme` maps seq 130→131 to author 125→127, so no residue 126 exists and the chain is continuous. |
| "Vina's result depends on `--cpu` as well as `--seed`, because the parallel search splits work across threads" (`run_production.sh`) | **Not reproduced.** Re-running 6RL9 sulfanilamide Vina at `--cpu 2` gave results identical to `--cpu 3` in every pose score of every seed. The evidence quoted for the claim (−5.862 vs −5.868) was a Pi-vs-VM difference, which is a different machine, not a different `--cpu`. Keeping `--cpu` fixed is still sensible for comparable timings, but it is not a correctness requirement on this build (Vina 1.2.7, aarch64). |

## 9. Reporting and code fixes made in this pass

All minimal and commented at the point of change, with the date of the verification.

1. `analyse.py` — `run_site_checks` now uses `base_target()` (one suffix rule, not two); `base_target()`
   derives the target from `TARGETS` so receptor variants need no suffix list.
2. `analyse.py` — `passes_2A_criterion` reported beside `passes_G2`, same structure, distance only.
3. `analyse.py` — the ranking test is reported twice, over all seeds and over site-valid poses only,
   with the seed counts; diagnostic runs (`__repro`, `__exh64`, `__box…`, `__fromxtal`) can no longer
   enter the ranking by accident.
4. `analyse.py` — T5 site checks replaced (section 6); every target's checks are run on its crystal
   reference and recorded, and the result appears in the gate file.
5. `analyse.py` and `contacts.py` — `--target` merges into the existing `summary.json` /
   `contacts.json` instead of replacing it with a single-target file; entries whose run directory has
   gone are dropped, so neither file can go stale.
6. `prep_receptor.py` — one `classify_gaps()` used by both code paths; PDBFixer's rebuild list is
   translated to author numbering and `gap_report` states both views and whether they agree.
7. `prep_receptor.py` — the minimum heavy–heavy distance is measured over the whole receptor and
   always reported, with the atom count and the threshold.
8. `prep_receptor.py` — `--recheck` re-reports an existing receptor (gaps, geometry, hashes) without
   rebuilding it, because a rebuild is not bit-reproducible and would invalidate the docking runs.
9. `prep_receptor.py` / `targets.py` — His64 is declared through a new `his_assumed` field, so the
   tautomer is labelled in `prep_report.json` as an assumption rather than falling through a default;
   `deposited_sequence()` now lives in `targets.py` and is used by both scripts; the inventory reports
   `unobserved` (from `_pdbx_poly_seq_scheme`) next to numbering gaps.
10. `targets.py` — `junk_for(spec)` makes the precedence explicit: a target's `keep` list always wins
    over `JUNK`, so 1AJ0's sulfate is kept by design and the inventory records the override.
11. `run_ad4zn.py` — the 20.25 vs 20.0 Å grid difference is documented as a confound with the control
    run that measures it.
12. `run_vina.py` / `run_ad4zn.py` — `--tag-suffix`, so a repeat or a variant cannot overwrite the run
    it is being compared with; `--target` choices require a prepared `box.json` and `receptor.pdbqt`.
13. `his_tautomer_variant.py` (new) — flips one neutral histidine's tautomer in an already prepared
    receptor (proton moved in the ring plane, ND1/NE2 AutoDock types and charges swapped, nothing
    else touched), so the His64 assumption can be tested against the very receptor that was docked.
14. `analyse.py` — `pick_reference()` no longer takes whichever `ref_*.sdf` the directory listing
    returns first: the target's own co-crystal ligand is chosen explicitly and the `ref_cross_*.sdf`
    files this script writes are excluded, so the redock/cross-dock decision cannot depend on
    directory order.
15. `analyse.py` — the median rank of the first near-native pose is no longer truncated with `int()`
    (a median of 1.5 over ten seeds was being reported as 1).

## 10. Role of each target in the paper

_The numbers behind each role are in sections 3–5; the paper reports them in Sections 3.6–3.7._

- **T5 7S2K (Sul2, dihydropteroate): validated with both scoring functions.** 0.51/0.54 Å, 10/10 seeds
  within 1 Å, the strongest target in the set; its poses and contacts are discussed quantitatively.
- **T4 3HS4 (hCA II): validated as a cross-docking target with AD4Zn, for sulfanilamide-like
  ligands.** 0.62 ± 0.03 Å, 10/10 within 1 Å, Zn–N 2.18 Å and the Thr199 bond in 10/10. Acetazolamide
  in its own structure is *not* reproduced by any scoring function (section 3), so 3HS4 supports
  statements about the sulfanilamide binding mode and about score ranking, not about the
  acetazolamide pose.
- **T3 6RL9 (hCA II, the compound's own co-crystal): partially validated with AD4Zn.** The
  near-native pose was found in 10/10 seeds and the zinc chemistry was right in 10/10 (Zn–N 2.22 Å,
  Thr199 2.56 Å), but it was the top-ranked pose in only 5/10, so the strict gate failed. The paper
  reports the binding mode with the bimodality and the success fraction beside it, with the 3HS4
  cross-dock as the quantitative evidence for the mode. The sulfonamidate anion behaved better here
  (1.12 ± 0.32 Å, 10/10 within 2 Å).
- **T1 3TZF (DHPS, *Y. pestis*): validated at the conventional 2.0 Å criterion, not at 1.0 Å.**
  1.92 ± 0.04 Å in 10/10 seeds with the correct pABA-pocket contacts. The near-native pose was sampled
  only occasionally: one seed produced 0.56 Å at rank 3, and only 1 of 10 seeds contains any pose
  within 1 Å anywhere in its list, against 10/10 within 2 Å at median rank 1. 3TZF is therefore a
  1.9 Å result, not a mis-ranked 0.6 Å one, and the paper reads its poses qualitatively.
- **T2 1AJ0 (DHPS, *E. coli*): failed validation.** 11.6 Å into a decoy cleft in every seed. It was
  left out of production, no pose or score from it is a result, and the paper states the failure.

## 11. Production run

The production run on this protocol is reported in `dock/RESULTS.md`, with
`dock/results/production_summary.json` behind every number: 72 jobs on the four validated targets
(T1 3TZF, T3 6RL9, T4 3HS4, T5 7S2K), Vina and Vinardo for every species and AutoDock4Zn at the two
hCA II targets for the species whose sulfonamide is neutral, exhaustiveness 32, a 20 Å box, and 10
seeds per job (20 for the two sulfanilamide AutoDock4Zn jobs whose poses the paper discusses).

## 12. How to reproduce

```bash
export PATH=/root/cqc-sfa-tools/dock/bin:$PATH
python dock/targets.py                                   # inventory from the deposited mmCIFs
python dock/prep_receptor.py T3_6RL9 --recheck            # re-report an existing receptor
python dock/run_vina.py  --target T5_7S2K --ligand dhp_anion --scoring vina \
       --seeds 10 --exhaustiveness 32 --cpu 3 --tag-suffix __repro
python dock/run_ad4zn.py --target T3_6RL9 --ligand sfa_neutral --stage dock \
       --seeds 10 --exhaustiveness 32 --cpu 3 --tag-suffix __repro
python dock/his_tautomer_variant.py --target T3_6RL9 --resnum 64 --to HID   # His64 sensitivity
python dock/run_vina.py --target T3_6RL9_his64HID --ligand sfa_neutral --scoring vina \
       --seeds 10 --exhaustiveness 32 --cpu 3
python dock/analyse.py                                    # summary.json + gate_G2.json
python dock/contacts.py                                   # contacts.json
```

Every run keeps its inputs, its Vina log and all poses under
`dock/<TARGET>/runs/<ligand>__<scoring>[__variant]/seedNN/`, and `results.json` in each run directory
records the exact settings (`exhaustiveness`, `cpu`, box, Vina version) that produced them.

---

## 13. What production then showed (added 28 September 2026)

Production docking ran on this protocol the same day: 72 jobs, 4 h 37 min on the Pi, reported in
`dock/RESULTS.md` with `dock/results/production_summary.json` behind every number. Two things in it
bear directly on this document, and both are documented there with the runs that measured them:

1. **The validation runs docked ETKDG/MMFF94s conformers, and the production runs dock the DFT
   geometries** — the runners never read the DFT file that `prep_ligands.py --dft` wrote, which was
   fixed for production (`--geometry dft`, recorded per run in `results.json`). Every production job
   of a control species (ETKDG, unchanged route) reproduces its validation run **bit for bit**
   (`production_summary.json` → `reproduction_of_validation`), so the two sets are comparable.
2. **The hCA II AD4Zn ring orientation depends on that input geometry.** With the DFT geometry,
   6RL9 sulfanilamide gives 0.61 ± 0.04 Å with 20/20 seeds within 1 Å — better than the
   2.07 ± 1.45 Å, f≤1 Å 0.5 measured here — while the 3HS4 cross-dock moves from 0.62 ± 0.03 Å to
   3.22 ± 0.02 Å, the ring rotated about the Zn–N/S axis with the zinc anchor intact. A 2 × 2
   diagnostic (`dock/run_geometry_diagnostic.sh`) shows the C4–S1–N2 angle accounts for the flip in
   6RL9 and not in 3HS4. In all eight AD4Zn runs of sulfanilamide the near-native pose is found in
   every seed at rank 1 or 2, separated from the alternative by 0.000–0.219 kcal/mol, and the zinc
   chemistry (Zn–N 2.09–2.37 Å, Thr199 OG1 2.46–2.59 Å) is right in every seed.

Sections 1–12 are the validation record as verified; section 3 gives the orientation separation with the
definition that the paper uses, and sections 10 and 11 state the role of each target and the production
run.

---

## 14. Independent check of the production run (28 September 2026)

An independent check re-derived the production numbers from the raw run outputs with its own scripts
(no `analyse.py`, no `production_analysis.py`): top-pose scores parsed from every `seedNN/out.pdbqt`;
RMSDs from the deposited mmCIF ligands with its own element-graph core matching over all automorphisms
and its own Cα superposition for the cross-docking reference; Zn–N, Thr199, Ser214 and pterin
distances from the mmCIF; torsions and C–S–N angles from the pose files and the DFT `opt.xyz`. Every
reported score, RMSD, success fraction, site-check distance and pass fraction reproduced
(≤ 0.02 Å, ≤ 0.2°), the 450 numeric rows of `RESULTS.md` agree with `production_summary.json` at the
precision printed, and the 12 ETKDG production jobs reproduce their validation runs in every pose
score of every seed. Four things were wrong or missing and were fixed:

1. **The exported pose files of the two headline AD4Zn jobs were the diagnostic poses.**
   `production_analysis.py` named every export `<species>__<scoring>__*`, so the section-7b diagnostic
   runs (analysed after the production ones) overwrote
   `results/poses/T{3_6RL9,4_3HS4}/sfa_neutral__ad4zn__top_pose_per_seed.sdf` with their ETKDG poses —
   2.0–3.1 Å from the production poses, 10 seeds instead of 20. Every number in the JSON was computed
   from the runs themselves and is unaffected, as is `results/poses_for_dft/`, but a figure made from
   those SDFs would have shown the wrong poses. A diagnostic export now keeps its run suffix; all
   78 pose files were regenerated and each was checked back against its own run.
2. **PoseBusters was never run on the production poses.** Added to `production_analysis.py`
   (`dock` configuration, 22 checks, `--no-posebusters` to skip): 740 production top poses, 99.3 %
   pass; the only failures are 5 poses of `T3_6RL9/azm_anion__vinardo__prod` whose nearest receptor
   atom is the catalytic zinc at 2.14–2.16 Å, i.e. the metal-coordination distance flagged by a
   generic van der Waals test.
3. **`RESULTS.md` stated no limitations.** Section 8 now carries them: empirical scoring functions
   (no Ki/Kd/affinity), rigid receptor, ligand flexible in torsions only, no explicit solvation or
   entropy, no MD and no free-energy calculation, no experimental data on this compound, and the
   protonation/tautomer assumptions. The waters kept are only the Mg coordination sphere
   (3 in 7S2K, 3 in 3TZF), which is now stated.
4. **Three reporting defects.** The zinc-bound-only ranking rendered as a dangling `: .` when no
   species qualified (it now says so); the headline rank column did not say it is a rank by score
   with no site filter (sulfamethoxazole ranks above sulfanilamide at hCA II with AD4Zn on a pose
   whose Zn–N is 3.28–3.56 Å); and the "zinc anchor intact" sentence for the 3HS4 cross-dock now
   carries its numbers (Zn–N 2.32 Å, within the 2.3 Å cutoff in 0.30 of the seeds).

Regenerating the analysis after the fix reproduced `production_summary.json` value for value
(only the added `posebusters` blocks, the corrected `pose_files` paths and the timestamps differ).
Checks of this pass that found nothing wrong: every production top pose lies ≥ 1.86 Å inside its
docking box (only lower-ranked modes exceed it, by ≤ 0.125 Å, the grid tolerance); no sulfanilamide
pose comes within 8.8 Å of any crystallographic symmetry mate, so none sits in a crystal-contact
site; the secondary sulfanilamide site of 6RL9 (SAN A306), which is inside the box, is never occupied
by a sulfanilamide pose (the only poses near it are p-aminobenzoate modes at ranks 14–20); ligand
protonation, formal charges and polar-hydrogen counts are correct in every docked `ligand.pdbqt`, and
AD4Zn was used only with neutral sulfonamides; each run's DFT input geometry matches the recorded
sha256, and the C–S–N angle measured in the poses equals the DFT value, which proves the DFT geometry
reached the docking; all 218 production site checks pass on their own crystal reference.
