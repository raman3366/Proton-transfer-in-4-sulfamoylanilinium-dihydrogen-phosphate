# Four ordered crystal waters in the 6RL9 receptor: does the ring-rotated pose of sulfanilamide survive them?

6 October 2026. Every number below comes from a file written in this test or from a stored file of the
repository; the file is named beside it.
Times are IST. Nothing that existed before was edited (`GLU106.md`, the receptors, the stored runs, the
manuscript and `md/` are as they were).

New files: `dock/water_variant.py`, `dock/run_water4_sensitivity.sh`, this file, and the variant
directories `dock/T3_6RL9_wat4/`, `dock/T3_6RL9_glu106GLU_wat4/`, `dock/T3_6RL9_wat4rev/`,
`dock/T3_6RL9_wat4bare/` and `dock/T3_6RL9_wat4fixed/`.

**Status (6 October 2026, 11:05 PM IST).** Everything that was queued is complete: the variants, their maps, the fixed-pose rescoring under all four hydrogen rules, and 16 docking runs in four water receptors. The runs are the two key runs (sulfanilamide from the DFT geometry, AutoDock4Zn, 20 seeds, both states of Glu106); the validation runs (ETKDG conformers, 10 seeds, both states) with the gates and PoseBusters; the key run under the hydrogen rules `rev` and `bare` (Glu106 as reported); the sulfonamidate with AutoDock4Zn (10 seeds, both states); and Vina and Vinardo for sulfanilamide (validation and production, 10 seeds each, both states). The queue finished at 11:38 AM IST on 6 October 2026, every job with return code 0 (`dock/T3_6RL9_wat4/run_water4.log`); the comparison files and the final tables (`dock/T3_6RL9_wat4/tables_after_queue.md`) were written between 11:39 and 11:41 AM (file times). A first version of this file was written at 4:46 AM IST, when only the two key runs were done; this version replaces it and holds every finished run. What was not tested is in section 9.

## 1. The question

AutoDock4Zn cannot separate two orientations of sulfanilamide in hCA II: the crystal one and one with
the ring turned about the Zn–N/S axis to the hydrophilic side (paper, Section 3.6; `dock/GLU106.md`).
The preparation of the molecular-dynamics system (`md/SYSTEM.md`, section 4; `md/system/audit/waters.csv`)
found that the ring-rotated pose lies on four water sites of the 6RL9 crystal: waters 423, 516, 553 and
568 of chain A, all at full occupancy with B 11.5–17.6 Å², present in the crystal together with the
ligand in its crystal pose. The reported docking used receptors without any water, so it could not
see them. The question here: does keeping these four waters remove the ring-rotated orientation from
the AutoDock4Zn results?

## 2. What was changed

`dock/water_variant.py --target <receptor> --waters 423 516 553 568 [--rule ...]` writes a variant that
is the prepared receptor, line for line, plus the named waters of the deposited file
(`dock/pdb/6RL9.cif.gz`) appended after the zinc. Nothing else changes: the 2510 lines of
`receptor.pdbqt` (2509 for the carboxylate variant) and the lines of `receptor.pdb` are byte-identical
and in order (`variant.json` → `file_check`), the box and the reference ligand are copied.

- **Frame.** The prepared receptor is in the frame of the deposited file: without any superposition
  the zinc is 0.000 Å from the deposited zinc, and 258 CA atoms differ by 0.013 Å rms, 0.032 Å at most,
  which is what the restrained minimisation of the receptor preparation did (`variant.json` →
  `frame_check`). The water oxygens are therefore written at their deposited coordinates.
- **Atoms.** Oxygen: AutoDock type OA, charge −0.411, the way `prep_receptor.py` and the gyrase B
  control (`control_gyrb.py`, `C1_4DUH_wat`) write a kept water. Hydrogens: type HD, charge +0.206, the
  hydrogen of Meeko's HOH template (0.2058, beside −0.4115 on its oxygen); O–H 0.9572 Å, H–O–H
  104.52°. A water is neutral to 0.001 e. The hydrogens follow their oxygen in the file, which is how
  AutoGrid finds the atoms bonded to an acceptor oxygen.
- **The precedent has no hydrogens** (a bare OA oxygen; it was used with Vina, which has no charge
  term). With AutoDock4Zn a bare oxygen leaves −0.411 e per water uncompensated (AutoGrid reports a
  total charge of −1.64 e for the receptor), and AutoGrid cannot give a hydrogen-bond direction to an
  acceptor oxygen without bonded atoms: it warns "Oxygen atom found with no bonded atoms" and then, 153 861
  times over the grid, that a division by zero "was just prevented" at these four oxygens
  (`T3_6RL9_wat4bare/ad4zn/maps_report.json`; that log was cut from 68 MB to its head, the counts and
  its tail). The waters therefore carry hydrogens in the main variants, and the bare oxygen is kept as
  one of the alternative rules. The maps of the variants with hydrogens were made without any warning.

Four hydrogen rules, each a stated construction (`water_variant.py`, module text):

| Rule | Variant | What it is |
|:--|:--|:--|
| `hbond` | `T3_6RL9_wat4`, `T3_6RL9_glu106GLU_wat4` | the orientation of the four waters together that makes the most hydrogen bonds with the atoms of the receptor and among themselves (H···A ≤ 2.5 Å, D–H···A ≥ 130°); exhaustive over 97–111 candidate orientations per water, 1.16 × 10⁸ combinations; two waters may not donate to each other along one line. The crystal ligand and the waters that are not kept take no part |
| `rev` | `T3_6RL9_wat4rev` | the best network of the same kind in which no water–water hydrogen bond of `hbond` is kept |
| `bare` | `T3_6RL9_wat4bare` | no hydrogens: the precedent |
| `fixed` | `T3_6RL9_wat4fixed` | every water in one and the same orientation, which looks at nothing in particular (H1 along +x, H2 in the xy plane); the kind of start a set-up program gives a crystal water before minimisation |

The four waters form a ring, 423 → 568 → 553 → 516 → 423, each donating one hydrogen to the next.
Water 423 also donates to Asn67 OD1 and accepts from Asn62 ND2–H; 553 donates to Thr200 OG1; 568
accepts from Gln92 NE2–H; one hydrogen of 516 and one of 568 point at nothing in the receptor (in the
crystal waters 430 and 582, which are not kept, lie 2.54 and 2.73 Å from these two;
`md/system/audit/waters.csv`). The receptor does not decide the direction
of the ring: turned the other way (`rev`) it makes the same number of hydrogen bonds, six from the
waters and two to them, and scores 6.572 against 6.599 on the tie-break of linearity
(`variant.json` → `hydrogen_search`). The two networks are therefore the two orientations that a
hydrogen-bond rule can give, and `rev` is the natural second rule. The orientation is the same for
both states of Glu106 (the water atoms of the two `hbond` variants are byte-identical); Glu106 is not
a neighbour of any of the four.

The orientations and the partners of every water under each rule (`variant.json` → `hydrogen_bonds`):

#### The waters and their hydrogen bonds, by hydrogen rule

| Variant (rule) | Water | B (Å²) | Donates to: acceptor, H···A (Å), O–H···A (°) | Accepts from: donor, H···O (Å), D–H···O (°) | Hydrogens without a partner | Hydrogen pairs closer than 2.0 Å |
|:--|:--|--:|:--|:--|--:|:--|
| `T3_6RL9_wat4` (hydrogen-bond network) | 423 | 12.7 | HOH568 O, 1.79, 176.8; ASN67 OD1, 1.66, 176.7 | ASN62 ND2-H, 2.05, 155.2; HOH516 H1, 1.71, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH423 O, 1.71, 180.0 | HOH553 H1, 1.70, 173.9 | 1 | – |
|  | 553 | 11.5 | HOH516 O, 1.70, 173.9; THR200 OG1, 2.06, 174.3 | HOH568 H1, 1.93, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH553 O, 1.93, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH423 H1, 1.79, 176.8 | 1 | – |
| `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 423 | 12.7 | HOH568 O, 1.79, 176.8; ASN67 OD1, 1.66, 176.7 | ASN62 ND2-H, 2.05, 155.2; HOH516 H1, 1.71, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH423 O, 1.71, 180.0 | HOH553 H1, 1.70, 173.9 | 1 | – |
|  | 553 | 11.5 | HOH516 O, 1.70, 173.9; THR200 OG1, 2.06, 174.3 | HOH568 H1, 1.93, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH553 O, 1.93, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH423 H1, 1.79, 176.8 | 1 | – |
| `T3_6RL9_wat4rev` (network, ring reversed) | 423 | 12.7 | HOH516 O, 1.79, 149.8; ASN67 OD1, 1.75, 149.5 | ASN62 ND2-H, 2.05, 155.2; HOH568 H1, 1.79, 180.0 | 0 | – |
|  | 516 | 17.6 | HOH553 O, 1.69, 180.0 | HOH423 H1, 1.79, 149.8 | 1 | – |
|  | 553 | 11.5 | HOH568 O, 1.93, 174.7; THR200 OG1, 2.06, 174.9 | HOH516 H1, 1.69, 180.0 | 0 | – |
|  | 568 | 12.0 | HOH423 O, 1.79, 180.0 | GLN92 NE2-H, 2.15, 144.7; HOH553 H1, 1.93, 174.7 | 1 | H2···GLN92 HE22 1.90 |
| `T3_6RL9_wat4bare` (bare oxygen (precedent)) | 423 | 12.7 | no hydrogens | ASN62 ND2-H, 2.05, 155.2 | – | – |
|  | 516 | 17.6 | no hydrogens | – | – | – |
|  | 553 | 11.5 | no hydrogens | – | – | – |
|  | 568 | 12.0 | no hydrogens | GLN92 NE2-H, 2.15, 144.7 | – | – |
| `T3_6RL9_wat4fixed` (one fixed orientation) | 423 | 12.7 | – | ASN62 ND2-H, 2.05, 155.2 | 2 | H1···ASN62 HD22 1.84; H2···ASN67 HD21 1.96 |
|  | 516 | 17.6 | – | – | 2 | – |
|  | 553 | 11.5 | THR200 OG1, 2.21, 141.9 | – | 1 | – |
|  | 568 | 12.0 | – | GLN92 NE2-H, 2.15, 144.7 | 2 | H2···GLN92 HE22 1.66 |

## 3. What was run

**Maps.** For every variant the zinc pseudo-atom was placed by the vendor script and AutoGrid was run
with the stored grid parameter file of the reported maps, unchanged, and the binary that
`dock/GLU106.md` validated against the stored maps (18 of 18 byte-identical; its sha256 is checked
before use). In every variant the pseudo-atom does not move (0.000 Å), the only atoms added to the
receptor are the water atoms, and no atom is lost (`<variant>/ad4zn/maps_report.json`).

**Settings.** Those of the stored runs, read from their `results.json` and used through the unchanged
drivers `run_ad4zn.py` and `run_vina.py`: exhaustiveness 32, the 20 Å box, 20 modes within 3 kcal mol⁻¹,
`--cpu 2` (shown in `GLU106.md` to give the same bytes as `--cpu 3`). Seeds 1–20 and the DFT geometry
for the AutoDock4Zn production run (the key run); seeds 1–10 for every other run, with ETKDG conformers
for the validation runs and the DFT geometry for the Vina and Vinardo production runs. AutoDock4Zn
docks into the maps of the variant, Vina and Vinardo into its `receptor.pdbqt`. In all 16 runs the
ligand input of every seed is byte-identical to that of the dry run it is compared with, and no setting
differs except `--cpu` against the stored runs of the paper (`comparison_water.json` → `check`).

**Dry runs.** The stored runs: `T3_6RL9` (Glu106 as reported, the runs of the paper) and
`T3_6RL9_glu106GLU` (carboxylate, the runs of `GLU106.md`). Nothing was docked again in a dry receptor.

**Runs analysed in this file** (16; `<variant>/runs/<run>/results.json`). In `T3_6RL9_wat4` and in
`T3_6RL9_glu106GLU_wat4`, seven each: `sfa_neutral__ad4zn__prod` (the key run: DFT geometry, 20 seeds),
`sfa_neutral__ad4zn` (validation: ETKDG conformers, 10 seeds), `sfa_anion__ad4zn` (the sulfonamidate:
ETKDG conformers, 10 seeds), `sfa_neutral__vina` and `sfa_neutral__vinardo` (validation, 10 seeds),
`sfa_neutral__vina__prod` and `sfa_neutral__vinardo__prod` (DFT geometry, 10 seeds). In `T3_6RL9_wat4rev`
and in `T3_6RL9_wat4bare`: `sfa_neutral__ad4zn__prod` (20 seeds). `T3_6RL9_wat4fixed` has maps and
fixed-pose scores and no docking run.

## 4. Does the crystal pose still dock?

A water variant that spoils the crystal pose would prove nothing, so this comes first.

- **Fixed poses.** The best near-native mode of each of the 20 seeds of the stored dry production run,
  scored without a search in the maps with the waters, changes by −0.089 (−0.100 to −0.070) kcal mol⁻¹ (Glu106 acid)
  and −0.084 (−0.098 to −0.067) kcal mol⁻¹ (carboxylate): the waters do not push the crystal pose out, they add to
  its score (`<variant>/rescore_fixed_poses.json`). The nearest of the four waters, 553, is 2.93 Å
  from the sulfonamide nitrogen of the crystal ligand (`md/system/audit/waters.csv`).
- **Docking, key run, Glu106 as reported.** RMSD of the top pose to the crystal ligand 0.61 ± 0.04 → 0.51 ± 0.04 Å (dry →
  water); top pose near-native in 20 → 20 of 20 seeds; top poses at 0.45–0.59 Å; Zn–N 2.20 → 2.16 Å,
  within the 2.3 Å cutoff in 20 of 20 seeds; Thr199 OG1···N 2.55 → 2.66 Å, within 3.3 Å in
  20 of 20. The sulfonamide torsion τ~fold~ of the top poses goes from 147.0 to 140.3° on average; the
  crystal ligand has 132.1° (`md/system/starts.csv`).
- **Docking, key run, Glu106 carboxylate.** RMSD 1.34 ± 1.30 → 0.51 ± 0.03 Å; top pose near-native in 15 → 20 of 20 seeds; top poses at 0.47–0.57 Å; Zn–N 2.21 → 2.16 Å, within the cutoff in 20 of 20 seeds; Thr199 OG1···N 2.55 → 2.65 Å, within 3.3 Å in 20 of 20.
- **Gates.** The paper applies its gates to the validation runs (ETKDG conformers, 10 seeds). The
  strict gate asks for a mean RMSD of the top pose of at most 1.0 Å, at least half of the seeds within
  1.0 Å, PoseBusters physics passed by at least 0.8 of the top poses, and the Zn–N check (2.3 Å) and
  the Thr199 OG1 check (3.3 Å) each passed by at least half of the seeds; the 2.0 Å criterion is the
  same with 2.0 Å in place of 1.0 Å (`analyse.py`; `glu_protonation_variant.py`, `gate`). The validation
  runs in the two water receptors are done; their verdicts stand beside those of the dry runs in the
  table below.

#### Gates of the validation runs (ETKDG conformers, 10 seeds; every cell is dry → water)

| Validation run | Glu106 | Strict gate | 2.0 Å criterion | PoseBusters physics | RMSD of the top pose (Å) | Seeds within 1.0 Å | Zn–N (Å) / fraction of seeds within 2.3 Å | Thr199 OG1···N (Å) / fraction of seeds within 3.3 Å |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| sulfanilamide, AutoDock4Zn | as reported | fail → pass | fail → pass | 1.0 → 1.0 | 2.07 ± 1.45 → 0.58 ± 0.03 | 5 of 10 → 10 of 10 | 2.22 / 1.00 → 2.14 / 1.00 | 2.56 / 1.00 → 2.68 / 1.00 |
|  | carboxylate | fail → pass | fail → pass | 1.0 → 1.0 | 2.70 ± 1.38 → 0.60 ± 0.03 | 3 of 10 → 10 of 10 | 2.23 / 1.00 → 2.15 / 1.00 | 2.57 / 1.00 → 2.67 / 1.00 |
| sulfonamidate, AutoDock4Zn | as reported | fail → pass | pass → pass | 1.0 → 1.0 | 1.12 ± 0.32 → 0.61 ± 0.04 | 3 of 10 → 10 of 10 | 2.21 / 1.00 → 2.23 / 1.00 | 2.60 / 1.00 → 2.59 / 1.00 |
|  | carboxylate | pass → pass | pass → pass | 1.0 → 1.0 | 0.86 ± 0.33 → 0.58 ± 0.04 | 7 of 10 → 10 of 10 | 2.22 / 1.00 → 2.23 / 1.00 | 2.60 / 1.00 → 2.59 / 1.00 |
| sulfanilamide, Vina | as reported | fail → fail | fail → fail | 1.0 → 1.0 | 1.10 ± 0.01 → 1.10 ± 0.01 | 0 of 10 → 0 of 10 | 2.84 / 0.00 → 2.93 / 0.00 | 4.76 / 0.00 → 4.81 / 0.00 |
|  | carboxylate | fail → fail | fail → fail | 1.0 → 1.0 | 1.09 ± 0.01 → 1.10 ± 0.01 | 0 of 10 → 0 of 10 | 2.83 / 0.00 → 2.90 / 0.00 | 4.77 / 0.00 → 4.80 / 0.00 |
| sulfanilamide, Vinardo | as reported | fail → fail | fail → fail | 1.0 → 1.0 | 1.13 ± 0.00 → 1.12 ± 0.00 | 0 of 10 → 0 of 10 | 2.48 / 0.00 → 2.50 / 0.00 | 4.45 / 0.00 → 4.47 / 0.00 |
|  | carboxylate | fail → fail | fail → fail | 1.0 → 1.0 | 1.12 ± 0.00 → 1.12 ± 0.01 | 0 of 10 → 0 of 10 | 2.48 / 0.00 → 2.50 / 0.00 | 4.45 / 0.00 → 4.47 / 0.00 |

Source: `<variant>/comparison_water.json` → `runs` → `<run>` → `dry` and `water` → `gate`,
`posebusters_physics_pass_fraction`, `rmsd`, `site_checks`.

PoseBusters was run on the top pose of each seed. For the water receptors it was computed in this
comparison against the `receptor.pdb` of the variant, which holds the four waters, and no physics check
failed in any of the eight runs (`posebusters_physics_failed_checks` is empty). The figures of the dry
runs are the stored ones (`dock/analysis/summary.json` for Glu106 as reported,
`dock/T3_6RL9_glu106GLU/comparison.json` for the carboxylate).

With AutoDock4Zn the validation run of sulfanilamide fails both criteria in both dry receptors and
passes both in both water receptors. The sulfonamidate goes from fail to pass on the strict gate with
Glu106 as reported and passes both criteria in the dry and in the water receptor with the carboxylate.
Vina and Vinardo fail both criteria in every receptor, dry and water: no seed is within 1.0 Å, and no
seed passes the Zn–N check or the Thr199 OG1 check, so the 2.0 Å criterion fails too although every top
pose is within 2.0 Å (section 5). A pass in a receptor that holds waters of the complex's own crystal
is not a validation of the protocol (section 8).

The production runs (DFT geometry) have no gate verdict: PoseBusters is run only on the validation
runs. The AutoDock4Zn production run meets the other criteria of the strict gate (mean RMSD at most
1.0 Å, at least half of the seeds within 1.0 Å, both site checks passed by at least half of the seeds)
in the four water receptors in which it was docked: `hbond` for both states of Glu106, `rev` and `bare`
for Glu106 as reported (`gate` → `rmsd_and_site_part_G2`). The Vina and Vinardo production runs do not
meet them in any receptor.

## 5. Dry against water, for both states of Glu106

Definitions. Near-native: a mode within 1.0 Å of the crystal ligand. Ring-rotated in the paper's
sense: a mode beyond 2.0 Å (`paper/main/data/derive/d14_completion.py`). In the dry AutoDock4Zn runs
of sulfanilamide that is mostly, but not only, the rotated pose (40 of 62 and 44 of 66 modes in the key
runs); in the dry sulfonamidate run with the Glu106 carboxylate (20 of 50 modes) and in the dry Vina
and Vinardo runs (14 to 23 of 167 to 175 modes) it is mostly something else. The tables therefore
split the modes beyond 2.0 Å into "rotated
pose" (within 2.0 Å of the reference rotated pose, the mode that `md/SYSTEM.md` takes as its start R:
seed 14, mode 3 of the reported production run, 3.57 Å from the crystal ligand) and "other". All RMSDs
are symmetry-corrected over heavy atoms in the receptor frame, by the validation code
(`analyse.symm_rmsd`); the gates and site checks are those of `analyse.py`. Each seed reports at most
20 modes within 3 kcal mol⁻¹ of its best one, so "all modes" means all modes that the program wrote.
Source: `<variant>/comparison_water.json`; the tables are printed by `water_variant.py --tables`.

The four tables below are the tables printed after the queue finished
(`dock/T3_6RL9_wat4/tables_after_queue.md`, from "Orientations" on) and hold every finished run: seven
runs for each state of Glu106 in the `hbond` receptors, and the key run in `rev` and in `bare`. Every
water row stands under the dry row it is compared with, so the dry row of the key run with Glu106 as
reported appears three times. The AutoDock4Zn production run has 20 seeds, every other run 10.
PoseBusters and the gate verdicts exist for the validation runs only; the family table lists the
AutoDock4Zn runs only.

#### Orientations: every mode of every seed

| Run | Glu106 | Receptor (hydrogen rule) | Seeds | Top pose: near-native / rotated pose / other beyond 2.0 Å / between | All modes: near-native / rotated pose / other beyond 2.0 Å / between | Seeds holding a mode of the rotated pose | Best score: near-native | Best score: rotated pose | Rotated pose behind the best near-native mode of its seed, min–max | Modes beyond 2.0 Å within 2 kcal mol⁻¹ of the seed best (of them the rotated pose) | Closest heavy atom of any mode to a water site (Å) |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|:--|:--|--:|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 20 | 20 / 0 / 0 / 0 | 39 / 0 / 1 / 11 | 0 | −9.238 | – | – | 0 (0) | 2.79 |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 5 / 5 / 0 / 0 | 20 / 20 / 5 / 15 | 10 | −8.874 | −8.888 | −0.039–0.020 | 23 (20) | 0.20 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 20 / 0 / 0 / 6 | 0 | −9.295 | – | – | 0 (0) | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 3 / 0 / 0 / 7 | 9 / 20 / 8 / 21 | 10 | −8.138 | −8.068 | 0.061–0.097 | 28 (20) | 0.34 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 18 / 0 / 0 / 4 | 0 | −8.142 | – | – | 0 (0) | 2.79 |
| sulfanilamide, Vina, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 160 / 26 | 9 | – | −5.517 | – | 174 (14) | 0.23 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 185 / 15 | 0 | – | – | – | 155 (0) | 2.68 |
| sulfanilamide, Vinardo, ETKDG (validation) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 2 / 23 / 145 / 30 | 10 | −4.679 | −4.442 | 0.237–0.237 | 168 (23) | 0.18 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 5 / 0 / 179 / 15 | 0 | −5.160 | – | – | 139 (0) | 2.46 |
| sulfanilamide, Vina, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 1 / 16 / 154 / 29 | 10 | −5.714 | −5.567 | 0.200–0.200 | 170 (16) | 0.22 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 16 | 0 | – | – | – | 137 (0) | 2.69 |
| sulfanilamide, Vinardo, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 156 / 30 | 9 | – | −4.473 | – | 170 (14) | 0.29 |
|  |  | water: `T3_6RL9_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 1 / 0 / 184 / 15 | 0 | −5.078 | – | – | 128 (0) | 2.53 |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 20 | 15 / 5 / 0 / 0 | 42 / 44 / 22 / 29 | 20 | −8.886 | −8.897 | −0.017–0.050 | 61 (41) | 0.09 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 20 | 20 / 0 / 0 / 0 | 38 / 0 / 1 / 10 | 0 | −9.247 | – | – | 1 (0) | 2.78 |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 3 / 7 / 0 / 0 | 21 / 20 / 7 / 16 | 10 | −8.919 | −8.945 | −0.046–0.005 | 27 (20) | 0.19 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 17 / 0 / 0 / 7 | 0 | −9.300 | – | – | 0 (0) | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 7 / 0 / 0 / 3 | 16 / 20 / 30 / 14 | 10 | −7.770 | −7.708 | 0.030–0.095 | 30 (20) | 0.32 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 10 / 0 / 0 / 0 | 18 / 0 / 0 / 5 | 0 | −7.772 | – | – | 0 (0) | 2.79 |
| sulfanilamide, Vina, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 0 / 19 / 156 / 25 | 10 | – | −5.517 | – | 175 (19) | 0.19 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 2 / 0 / 183 / 15 | 0 | −6.184 | – | – | 150 (0) | 2.68 |
| sulfanilamide, Vinardo, ETKDG (validation) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 1 / 19 / 148 / 31 | 10 | −4.582 | −4.457 | 0.151–0.151 | 167 (19) | 0.18 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 2 / 0 / 181 / 17 | 0 | −5.107 | – | – | 138 (0) | 2.52 |
| sulfanilamide, Vina, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 0 / 14 / 157 / 29 | 10 | – | −5.547 | – | 171 (14) | 0.24 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 17 | 0 | – | – | – | 131 (0) | 2.70 |
| sulfanilamide, Vinardo, DFT geometry (production) | carboxylate | dry: `T3_6RL9_glu106GLU` | 10 | 0 / 0 / 0 / 10 | 1 / 20 / 153 / 26 | 10 | −4.668 | −4.443 | 0.245–0.245 | 173 (20) | 0.25 |
|  |  | water: `T3_6RL9_glu106GLU_wat4` (hydrogen-bond network) | 10 | 0 / 0 / 0 / 10 | 0 / 0 / 182 / 18 | 0 | – | – | – | 130 (0) | 2.49 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4rev` (network, ring reversed) | 20 | 20 / 0 / 0 / 0 | 39 / 0 / 1 / 7 | 0 | −9.249 | – | – | 0 (0) | 2.79 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry: `T3_6RL9` | 20 | 20 / 0 / 0 / 0 | 41 / 40 / 22 / 35 | 20 | −8.843 | −8.823 | 0.014–0.070 | 57 (39) | 0.07 |
|  |  | water: `T3_6RL9_wat4bare` (bare oxygen (precedent)) | 20 | 20 / 0 / 0 / 0 | 41 / 0 / 0 / 6 | 0 | −9.542 | – | – | 0 (0) | 2.51 |

#### Does the crystal pose still dock: top poses, contacts, gates

| Run | Glu106 | Receptor | RMSD top (Å) | best of seeds | f ≤ 1 Å (≤ 2 Å) | Score (kcal mol⁻¹) | Zn–N (Å) / pass | Thr199 OG1···N (Å) / pass | Best mode of each seed (Å), min–max | PoseBusters physics | Strict gate / 2.0 Å |
|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, hydrogen-bond network | 0.51 ± 0.04 | 0.45 | 1.00 (1.00) | −9.232 ± 0.004 | 2.16 / 1.00 | 2.66 / 1.00 | 0.32–0.51 | – | – |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry | 2.07 ± 1.45 | 0.65 | 0.5 (0.5) | −8.862 ± 0.017 | 2.22 / 1.00 | 2.56 / 1.00 | 0.40–0.73 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 0.58 ± 0.03 | 0.53 | 1.0 (1.0) | −9.286 ± 0.005 | 2.14 / 1.00 | 2.68 / 1.00 | 0.48–0.58 | 1.0 | pass / pass |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1.12 ± 0.32 | 0.61 | 0.3 (1.0) | −8.133 ± 0.006 | 2.21 / 1.00 | 2.60 / 1.00 | 0.60–1.29 | 1.0 | fail / pass |
|  |  | water, hydrogen-bond network | 0.61 ± 0.04 | 0.52 | 1.0 (1.0) | −8.131 ± 0.009 | 2.23 / 1.00 | 2.59 / 1.00 | 0.46–0.62 | 1.0 | pass / pass |
| sulfanilamide, Vina, ETKDG (validation) | acid (reported) | dry | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −5.870 ± 0.012 | 2.84 / 0.00 | 4.76 / 0.00 | 1.09–1.11 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −6.279 ± 0.061 | 2.93 / 0.00 | 4.81 / 0.00 | 1.09–1.11 | 1.0 | fail / fail |
| sulfanilamide, Vinardo, ETKDG (validation) | acid (reported) | dry | 1.13 ± 0.00 | 1.12 | 0.0 (1.0) | −5.075 ± 0.010 | 2.48 / 0.00 | 4.45 / 0.00 | 0.35–1.13 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.12 ± 0.00 | 1.11 | 0.0 (1.0) | −5.498 ± 0.013 | 2.50 / 0.00 | 4.47 / 0.00 | 0.26–1.12 | 1.0 | fail / fail |
| sulfanilamide, Vina, DFT geometry (production) | acid (reported) | dry | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −5.880 ± 0.012 | 2.75 / 0.00 | 4.68 / 0.00 | 0.33–1.08 | – | – |
|  |  | water, hydrogen-bond network | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −6.348 ± 0.004 | 2.75 / 0.00 | 4.67 / 0.00 | 1.08–1.08 | – | – |
| sulfanilamide, Vinardo, DFT geometry (production) | acid (reported) | dry | 1.11 ± 0.00 | 1.11 | 0.0 (1.0) | −5.077 ± 0.004 | 2.40 / 0.00 | 4.35 / 0.00 | 1.10–1.12 | – | – |
|  |  | water, hydrogen-bond network | 1.10 ± 0.00 | 1.10 | 0.0 (1.0) | −5.537 ± 0.004 | 2.43 / 0.00 | 4.37 / 0.00 | 0.25–1.11 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry | 1.34 ± 1.30 | 0.41 | 0.75 (0.75) | −8.882 ± 0.007 | 2.21 / 1.00 | 2.55 / 1.00 | 0.41–0.68 | – | – |
|  |  | water, hydrogen-bond network | 0.51 ± 0.03 | 0.47 | 1.00 (1.00) | −9.241 ± 0.004 | 2.16 / 1.00 | 2.65 / 1.00 | 0.30–0.56 | – | – |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry | 2.70 ± 1.38 | 0.67 | 0.3 (0.3) | −8.911 ± 0.021 | 2.23 / 1.00 | 2.57 / 1.00 | 0.66–0.73 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 0.60 ± 0.03 | 0.55 | 1.0 (1.0) | −9.291 ± 0.004 | 2.15 / 1.00 | 2.67 / 1.00 | 0.43–0.62 | 1.0 | pass / pass |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry | 0.86 ± 0.33 | 0.58 | 0.7 (1.0) | −7.764 ± 0.003 | 2.22 / 1.00 | 2.60 / 1.00 | 0.45–0.71 | 1.0 | pass / pass |
|  |  | water, hydrogen-bond network | 0.58 ± 0.04 | 0.53 | 1.0 (1.0) | −7.764 ± 0.006 | 2.23 / 1.00 | 2.59 / 1.00 | 0.52–0.64 | 1.0 | pass / pass |
| sulfanilamide, Vina, ETKDG (validation) | carboxylate | dry | 1.09 ± 0.01 | 1.09 | 0.0 (1.0) | −5.871 ± 0.011 | 2.83 / 0.00 | 4.77 / 0.00 | 1.09–1.10 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.10 ± 0.01 | 1.09 | 0.0 (1.0) | −6.302 ± 0.054 | 2.90 / 0.00 | 4.80 / 0.00 | 0.33–1.11 | 1.0 | fail / fail |
| sulfanilamide, Vinardo, ETKDG (validation) | carboxylate | dry | 1.12 ± 0.00 | 1.12 | 0.0 (1.0) | −5.078 ± 0.010 | 2.48 / 0.00 | 4.45 / 0.00 | 0.46–1.13 | 1.0 | fail / fail |
|  |  | water, hydrogen-bond network | 1.12 ± 0.01 | 1.11 | 0.0 (1.0) | −5.498 ± 0.012 | 2.50 / 0.00 | 4.47 / 0.00 | 0.26–1.12 | 1.0 | fail / fail |
| sulfanilamide, Vina, DFT geometry (production) | carboxylate | dry | 1.08 ± 0.00 | 1.08 | 0.0 (1.0) | −5.878 ± 0.012 | 2.75 / 0.00 | 4.68 / 0.00 | 1.08–1.08 | – | – |
|  |  | water, hydrogen-bond network | 1.09 ± 0.00 | 1.08 | 0.0 (1.0) | −6.348 ± 0.004 | 2.74 / 0.00 | 4.66 / 0.00 | 1.08–1.09 | – | – |
| sulfanilamide, Vinardo, DFT geometry (production) | carboxylate | dry | 1.11 ± 0.00 | 1.11 | 0.0 (1.0) | −5.076 ± 0.003 | 2.40 / 0.00 | 4.35 / 0.00 | 0.37–1.11 | – | – |
|  |  | water, hydrogen-bond network | 1.11 ± 0.00 | 1.10 | 0.0 (1.0) | −5.541 ± 0.003 | 2.42 / 0.00 | 4.37 / 0.00 | 1.10–1.11 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, network, ring reversed | 0.51 ± 0.03 | 0.45 | 1.00 (1.00) | −9.242 ± 0.003 | 2.16 / 1.00 | 2.65 / 1.00 | 0.37–0.53 | – | – |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 0.61 ± 0.04 | 0.56 | 1.00 (1.00) | −8.836 ± 0.004 | 2.20 / 1.00 | 2.55 / 1.00 | 0.47–0.66 | – | – |
|  |  | water, bare oxygen (precedent) | 0.52 ± 0.04 | 0.42 | 1.00 (1.00) | −9.538 ± 0.003 | 2.12 / 1.00 | 2.70 / 1.00 | 0.36–0.53 | – | – |

#### What the search finds besides the crystal pose: families of modes beyond 1.0 Å (2.0 Å clusters, all seeds pooled)

| Run | Glu106 | Receptor | Family | Modes (seeds; top poses) | Best score | Behind the best mode of the run | RMSD to the crystal ligand (Å) | RMSD to the rotated reference (Å) | Zn–N (Å) | Aniline N···Asn67 OD1 (Å) | Nearest water site (Å) |
|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, hydrogen-bond network | 1: between | 11 (11; 0) | −8.857 | 0.381 | 1.57 | 3.61 | 2.28 | 7.12 | 2.79 |
|  |  |  | 2: between | 1 (1; 0) | −7.721 | 1.517 | 1.65 | 3.74 | 2.23 | 8.88 | 2.80 |
| sulfanilamide, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1: rotated pose | 21 (10; 5) | −8.888 | 0.000 | 3.37 | 0.26 | 2.23 | 2.70 | 0.23 |
|  |  |  | 2: between | 14 (8; 0) | −8.400 | 0.488 | 1.15 | 2.77 | 2.13 | 5.86 | 1.61 |
|  |  |  | 3: other, beyond 2.0 A | 2 (2; 0) | −7.728 | 1.160 | 3.16 | 5.09 | 2.21 | 10.62 | 1.32 |
|  |  |  | 4: other, beyond 2.0 A | 3 (3; 0) | −7.443 | 1.445 | 2.17 | 3.72 | 2.00 | 7.74 | 0.91 |
|  |  | water, hydrogen-bond network | 1: between | 6 (5; 0) | −8.509 | 0.786 | 1.58 | 3.58 | 2.32 | 6.95 | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | acid (reported) | dry | 1: between | 19 (9; 7) | −8.140 | 0.000 | 1.37 | 3.86 | 2.18 | 7.21 | 2.48 |
|  |  |  | 2: rotated pose | 20 (10; 0) | −8.068 | 0.072 | 3.32 | 0.32 | 2.24 | 2.81 | 0.41 |
|  |  |  | 3: between | 3 (3; 0) | −7.310 | 0.830 | 1.67 | 3.06 | 2.13 | 7.41 | 1.89 |
|  |  |  | 4: other, beyond 2.0 A | 4 (4; 0) | −6.937 | 1.203 | 3.12 | 5.00 | 1.99 | 10.45 | 1.31 |
|  |  |  | 5: other, beyond 2.0 A | 3 (3; 0) | −6.875 | 1.265 | 2.20 | 3.54 | 2.36 | 7.23 | 0.84 |
|  |  | water, hydrogen-bond network | 1: between | 3 (3; 0) | −7.656 | 0.486 | 1.59 | 3.55 | 2.28 | 7.00 | 2.79 |
|  |  |  | 2: between | 1 (1; 0) | −6.910 | 1.232 | 1.52 | 3.70 | 2.25 | 8.73 | 2.80 |
| sulfanilamide, AD4Zn, DFT geometry (production) | carboxylate | dry | 1: rotated pose | 45 (20; 5) | −8.897 | 0.000 | 3.51 | 0.08 | 2.24 | 2.58 | 0.16 |
|  |  |  | 2: between | 26 (18; 0) | −8.445 | 0.452 | 1.02 | 2.81 | 2.17 | 5.97 | 1.92 |
|  |  |  | 3: between | 1 (1; 0) | −8.021 | 0.876 | 1.85 | 3.28 | 2.11 | 8.00 | 1.98 |
|  |  |  | 4: other, beyond 2.0 A | 12 (12; 0) | −7.722 | 1.175 | 3.07 | 4.98 | 2.04 | 10.39 | 1.23 |
|  |  |  | 5: other, beyond 2.0 A | 11 (11; 0) | −7.644 | 1.253 | 2.17 | 3.61 | 2.12 | 7.44 | 0.84 |
|  |  | water, hydrogen-bond network | 1: between | 9 (9; 0) | −8.858 | 0.389 | 1.59 | 3.53 | 2.31 | 6.99 | 2.79 |
|  |  |  | 2: other, beyond 2.0 A | 1 (1; 0) | −7.770 | 1.477 | 2.16 | 4.88 | 2.21 | 9.11 | 3.01 |
|  |  |  | 3: between | 1 (1; 0) | −7.621 | 1.626 | 1.57 | 3.65 | 2.16 | 8.62 | 2.78 |
| sulfanilamide, AD4Zn, ETKDG (validation) | carboxylate | dry | 1: rotated pose | 21 (10; 7) | −8.945 | 0.000 | 3.41 | 0.22 | 2.24 | 2.71 | 0.20 |
|  |  |  | 2: between | 14 (10; 0) | −8.408 | 0.537 | 1.09 | 2.73 | 2.15 | 5.95 | 1.75 |
|  |  |  | 3: other, beyond 2.0 A | 5 (4; 0) | −7.669 | 1.276 | 3.16 | 5.08 | 2.20 | 10.61 | 1.33 |
|  |  |  | 4: other, beyond 2.0 A | 3 (3; 0) | −7.358 | 1.587 | 2.14 | 3.65 | 1.94 | 7.54 | 0.87 |
|  |  | water, hydrogen-bond network | 1: between | 7 (7; 0) | −8.483 | 0.817 | 1.59 | 3.51 | 2.31 | 6.84 | 2.79 |
| sulfonamidate, AD4Zn, ETKDG (validation) | carboxylate | dry | 1: between | 6 (4; 3) | −7.767 | 0.003 | 1.39 | 3.93 | 2.15 | 7.34 | 2.57 |
|  |  |  | 2: rotated pose | 21 (10; 0) | −7.708 | 0.062 | 3.28 | 0.38 | 2.24 | 2.95 | 0.39 |
|  |  |  | 3: between | 8 (8; 0) | −7.041 | 0.729 | 1.34 | 3.22 | 2.15 | 7.60 | 2.16 |
|  |  |  | 4: other, beyond 2.0 A | 5 (5; 0) | −6.584 | 1.186 | 3.11 | 4.99 | 2.03 | 10.40 | 1.33 |
|  |  |  | 5: other, beyond 2.0 A | 6 (5; 0) | −6.558 | 1.212 | 2.15 | 3.62 | 2.01 | 7.45 | 0.84 |
|  |  |  | 6: other, beyond 2.0 A | 18 (10; 0) | −4.869 | 2.901 | 5.51 | 4.48 | 6.12 | 8.38 | 0.57 |
|  |  | water, hydrogen-bond network | 1: between | 3 (3; 0) | −7.353 | 0.419 | 1.67 | 3.67 | 2.30 | 7.15 | 2.79 |
|  |  |  | 2: between | 2 (2; 0) | −7.034 | 0.738 | 1.10 | 3.65 | 2.01 | 8.40 | 3.17 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, network, ring reversed | 1: between | 8 (8; 0) | −8.822 | 0.427 | 1.59 | 3.58 | 2.29 | 7.10 | 2.79 |
| sulfanilamide, AD4Zn, DFT geometry (production) | acid (reported) | dry | 1: rotated pose | 40 (20; 0) | −8.823 | 0.020 | 3.57 | 0.00 | 2.24 | 2.57 | 0.19 |
|  |  |  | 2: between | 39 (20; 0) | −8.404 | 0.439 | 1.17 | 2.75 | 2.19 | 5.79 | 1.49 |
|  |  |  | 3: other, beyond 2.0 A | 8 (7; 0) | −7.657 | 1.186 | 2.14 | 3.64 | 2.00 | 7.52 | 0.90 |
|  |  |  | 4: other, beyond 2.0 A | 9 (9; 0) | −7.635 | 1.208 | 3.10 | 4.99 | 1.95 | 10.42 | 1.18 |
|  |  |  | 5: other, beyond 2.0 A | 1 (1; 0) | −6.115 | 2.728 | 2.91 | 2.15 | 1.86 | 4.59 | 0.55 |
|  |  | water, bare oxygen (precedent) | 1: between | 4 (4; 0) | −9.144 | 0.398 | 1.60 | 3.59 | 2.25 | 7.08 | 2.79 |
|  |  |  | 2: between | 2 (2; 0) | −8.219 | 1.323 | 1.62 | 3.72 | 2.25 | 8.76 | 2.80 |

#### Stored modes of the dry runs scored in the maps with the waters (no search)

| Maps | Stored dry run | Near-native modes: n, score in the dry maps → in the water maps (change) | Modes of the rotated pose: n, score in the water maps, mean (min to max) |
|:--|:--|:--|:--|
| `T3_6RL9_wat4` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.925 (−0.100 to −0.070) | 20, +6840 (+6000 to +7432) |
| `T3_6RL9_glu106GLU_wat4` | `T3_6RL9_glu106GLU`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.878 → −8.962 (−0.098 to −0.067) | 20, +6818 (+6472 to +7360) |
| `T3_6RL9_wat4rev` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.926 (−0.102 to −0.073) | 20, +6867 (+5846 to +7209) |
| `T3_6RL9_wat4bare` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −9.175 (−0.358 to −0.313) | 20, +5877 (+5346 to +6858) |
| `T3_6RL9_wat4fixed` | `T3_6RL9`, sulfanilamide, AD4Zn, DFT geometry (production) | 20, −8.836 → −8.923 (−0.103 to −0.068) | 20, +7151 (+6410 to +7644) |

Reading of the key run (sulfanilamide from the DFT geometry, AutoDock4Zn, 20 seeds, hydrogen rule `hbond`):

| | Glu106 as reported | Glu106 carboxylate |
|:--|:--|:--|
| modes of the rotated pose, dry → water | 40 → 0 | 44 → 0 |
| modes beyond 2.0 Å (the paper's ring-rotated), dry → water | 62 → 1 | 66 → 1 |
| of them within 2 kcal mol⁻¹ of the best mode of their seed | 57 → 0 | 61 → 1 |
| the modes beyond 2.0 Å that remain with the waters | seed 13, mode 3: 2.23 Å from the crystal ligand, 4.60 Å from the rotated reference, 2.311 kcal mol⁻¹ behind the top mode of its seed, Zn–N 1.84 Å, aniline N 8.96 Å from Asn67 OD1 | seed 2, mode 3: 2.16 Å from the crystal ligand, 4.88 Å from the rotated reference, 1.477 kcal mol⁻¹ behind the top mode of its seed, Zn–N 2.21 Å, aniline N 9.11 Å from Asn67 OD1 |
| seeds whose top pose is near-native | 20 → 20 of 20 | 15 → 20 of 20 |
| seeds whose top pose is the rotated pose | 0 → 0 | 5 → 0 |
| RMSD of the top pose (Å) | 0.61 ± 0.04 → 0.51 ± 0.04 | 1.34 ± 1.30 → 0.51 ± 0.03 |
| score of the top pose (kcal mol⁻¹) | −8.836 ± 0.004 → −9.232 ± 0.004 | −8.882 ± 0.007 → −9.241 ± 0.004 |
| best near-native mode (kcal mol⁻¹) | −8.843 → −9.238 | −8.886 → −9.247 |
| best mode of the rotated pose in the dry receptor (kcal mol⁻¹) | −8.823 | −8.897 |
| all modes written, and per seed | 138 → 51; 6.90 → 2.55 | 137 → 49; 6.85 → 2.45 |
| largest score span of a seed's list (kcal mol⁻¹), dry / water | 2.72 / 2.31 | 2.38 / 2.72 |
| modes with a heavy atom within 2.4 Å of a water site, dry → water | 84 → 0 | 86 → 0 |
| closest heavy atom of any mode to a water site (Å), dry / water | 0.07 / 2.79 | 0.09 / 2.78 |

Reading of the validation runs with AutoDock4Zn (ETKDG conformers, 10 seeds, hydrogen rule `hbond`;
every cell is dry → water unless the row says otherwise):

| | Sulfanilamide, Glu106 as reported | Sulfanilamide, Glu106 carboxylate | Sulfonamidate, Glu106 as reported | Sulfonamidate, Glu106 carboxylate |
|:--|:--|:--|:--|:--|
| seeds whose top pose is near-native | 5 → 10 of 10 | 3 → 10 of 10 | 3 → 10 of 10 | 7 → 10 of 10 |
| seeds whose top pose is the rotated pose | 5 → 0 | 7 → 0 | 0 → 0 | 0 → 0 |
| seeds whose top pose lies between 1.0 and 2.0 Å | 0 → 0 | 0 → 0 | 7 → 0 | 3 → 0 |
| modes of the rotated pose | 20 → 0 | 20 → 0 | 20 → 0 | 20 → 0 |
| seeds holding a mode of the rotated pose | 10 → 0 | 10 → 0 | 10 → 0 | 10 → 0 |
| modes beyond 2.0 Å (the paper's ring-rotated) | 25 → 0 | 27 → 0 | 28 → 0 | 50 → 0 |
| modes between 1.0 and 2.0 Å | 15 → 6 | 16 → 7 | 21 → 4 | 14 → 5 |
| near-native modes | 20 → 20 | 21 → 17 | 9 → 18 | 16 → 18 |
| all modes written | 60 → 26 | 64 → 24 | 58 → 22 | 80 → 23 |
| RMSD of the top pose (Å) | 2.07 ± 1.45 → 0.58 ± 0.03 | 2.70 ± 1.38 → 0.60 ± 0.03 | 1.12 ± 0.32 → 0.61 ± 0.04 | 0.86 ± 0.33 → 0.58 ± 0.04 |
| score of the top pose (kcal mol⁻¹) | −8.862 ± 0.017 → −9.286 ± 0.005 | −8.911 ± 0.021 → −9.291 ± 0.004 | −8.133 ± 0.006 → −8.131 ± 0.009 | −7.764 ± 0.003 → −7.764 ± 0.006 |
| best near-native mode (kcal mol⁻¹) | −8.874 → −9.295 | −8.919 → −9.300 | −8.138 → −8.142 | −7.770 → −7.772 |
| best mode of the rotated pose in the dry receptor (kcal mol⁻¹) | −8.888 | −8.945 | −8.068 | −7.708 |
| rotated pose behind the best near-native mode of its seed in the dry receptor, min to max (seeds that hold both) | −0.039 to 0.020 (10) | −0.046 to 0.005 (10) | 0.061 to 0.097 (7) | 0.030 to 0.095 (10) |
| modes with a heavy atom within 2.4 Å of a water site | 36 → 0 | 37 → 0 | 38 → 0 | 59 → 0 |
| closest heavy atom of any mode to a water site (Å), dry / water | 0.20 / 2.79 | 0.19 / 2.79 | 0.34 / 2.79 | 0.32 / 2.79 |

**What the search does instead (AutoDock4Zn).** No pose other than the crystal pose is the top pose of
any seed, and in the key runs the largest family that remains is, by its stored distances, of a kind
the dry runs also hold. With the waters each seed
writes fewer modes (6.90 → 2.55 per seed in the key run for Glu106 as reported), because the rotated
pose and most of the other alternatives are no longer there. Beyond 1.0 Å from the crystal ligand what
remains in the key run, pooled over the seeds in 2.0 Å clusters, is, for Glu106 as reported: 11 modes in 11 seeds (representative 1.57 Å from the crystal ligand, 3.61 Å from the rotated reference; best −8.857, 0.381 behind the best mode; Zn–N 2.28 Å; aniline N to Asn67 OD1 7.12 Å; nearest water site 2.79 Å); 1 mode in 1 seed (representative 1.65 Å from the crystal ligand, 3.74 Å from the rotated reference; best −7.721, 1.517 behind the best mode; Zn–N 2.23 Å; aniline N to Asn67 OD1 8.88 Å; nearest water site 2.80 Å). For the carboxylate:
9 modes in 9 seeds (representative 1.59 Å from the crystal ligand, 3.53 Å from the rotated reference; best −8.858, 0.389 behind the best mode; Zn–N 2.31 Å; aniline N to Asn67 OD1 6.99 Å; nearest water site 2.79 Å); 1 mode in 1 seed (representative 2.16 Å from the crystal ligand, 4.88 Å from the rotated reference; best −7.770, 1.477 behind the best mode; Zn–N 2.21 Å; aniline N to Asn67 OD1 9.11 Å; nearest water site 3.01 Å); 1 mode in 1 seed (representative 1.57 Å from the crystal ligand, 3.65 Å from the rotated reference; best −7.621, 1.626 behind the best mode; Zn–N 2.16 Å; aniline N to Asn67 OD1 8.62 Å; nearest water site 2.78 Å). In these modes the sulfonamide nitrogen is still on the zinc and the ring is
still on the side of the crystal pose: the aniline N is 7–9 Å from Asn67 OD1, as in the crystal pose
(7.74 Å, `md/system/starts.csv`), not 2.6 Å as in the rotated pose. The dry key runs hold modes with
distances close to those of the largest of these families: for Glu106 as reported, for example, seed 8,
mode 6 (1.59 Å from the crystal ligand, 3.50 Å from the rotated reference, Zn–N 2.27 Å, aniline N
6.96 Å from Asn67 OD1, nearest water site 2.79 Å), and for the carboxylate seed 7, mode 6 (1.70 Å,
3.65 Å, 2.30 Å, 7.06 Å, 2.79 Å); 14 of the 35 dry modes between 1.0 and 2.0 Å, and 10 of 29 for the
carboxylate, have no heavy atom within 2.4 Å of a water site (`comparison_water.json` → `dry` →
`modes`). No RMSD between a mode of a water run and a mode of its dry run was computed, and the
families of one mode were not compared with the dry lists. In the validation runs and in the sulfonamidate
runs nothing else is left beyond 1.0 Å with the waters: one family for sulfanilamide (6 modes in 5
seeds, 0.786 kcal mol⁻¹ behind the best mode of the run, for Glu106 as reported; 7 modes in 7 seeds,
0.817 behind, for the carboxylate) and two for the sulfonamidate (3 modes and 1 mode; 3 modes and 2
modes), all between 1.0 and 2.0 Å from the crystal ligand, none of them a top pose, with Zn–N 2.01 to
2.32 Å and the aniline N 6.84 to 8.73 Å from Asn67 OD1 at their representatives (family table).

**The validation run of sulfanilamide (AutoDock4Zn, ETKDG conformers, 10 seeds).** This is the run to
which the paper applies its gates, and in the dry receptors it is the run in which the two orientations
compete most: the rotated pose is the top pose in 5 of 10 seeds (Glu106 as reported) and in 7 of 10
(carboxylate), and in a seed the best mode of the rotated pose scores from 0.039 kcal mol⁻¹ better to
0.020 worse than the best near-native mode (0.046 better to 0.005 worse for the carboxylate). With the waters
neither receptor gives a mode of the rotated pose or any mode beyond 2.0 Å. The top pose is near-native
in 10 of 10 seeds in both, at 0.58 ± 0.03 and 0.60 ± 0.03 Å (2.07 ± 1.45 and 2.70 ± 1.38 Å dry), with
Zn–N 2.14 and 2.15 Å, and the run passes both criteria in both water receptors (section 4). The answer
of the key run is thus also obtained from the ETKDG conformers, with 10 seeds.

**The sulfonamidate (AutoDock4Zn, ETKDG conformers, 10 seeds, both states of Glu106).** In the dry
receptors the rotated pose is never the top pose of the sulfonamidate, but it is in the list of every
seed (20 modes in 10 seeds), 0.061 to 0.097 kcal mol⁻¹ behind the best near-native mode of its seed
with Glu106 as reported (the 7 seeds that hold a near-native mode) and 0.030 to 0.095 behind with the
carboxylate (10 seeds). What takes the top pose from the crystal pose there is a pose between 1.0 and
2.0 Å, in 7 of 10 seeds (as reported) and 3 of 10 (carboxylate); the representative of its family is
1.37 and 1.39 Å from the crystal ligand and has a heavy atom 2.48 and 2.57 Å from a water site (family
table). What changes with the waters: modes of the rotated pose 20 → 0 in both receptors; modes beyond
2.0 Å 28 → 0 and 50 → 0; top pose near-native in 3 → 10 of 10 and 7 → 10 of 10 seeds; RMSD of the top
pose 1.12 ± 0.32 → 0.61 ± 0.04 Å and 0.86 ± 0.33 → 0.58 ± 0.04 Å; the strict gate goes from fail to
pass with Glu106 as reported. What does not change: the score of the top pose (−8.133 ± 0.006 →
−8.131 ± 0.009 and −7.764 ± 0.003 → −7.764 ± 0.006 kcal mol⁻¹, where that of neutral sulfanilamide in
the validation run goes from −8.862 to −9.286 and from −8.911 to −9.291), the Zn–N distance (2.21 →
2.23 and 2.22 → 2.23 Å) and the distance to Thr199 OG1 (2.60 → 2.59 Å in both). The difference between
the two states of Glu106 in the dry receptors (3 against 7 of 10 seeds near-native) is not there with
the waters (10 and 10). The sulfonamidate was run with the waters only from ETKDG conformers, only with
AutoDock4Zn and only under the rule `hbond`.

**Vina and Vinardo (sulfanilamide; validation runs from ETKDG conformers and production runs from the
DFT geometry, 10 seeds each, both states of Glu106, rule `hbond`).** Every cell of the table is dry →
water (`comparison_water.json` → `mode_summary`, `rmsd`, `site_checks`):

| Run (10 seeds) | Glu106 | Top pose: near-native / rotated pose / other beyond 2.0 Å / between | Modes of the rotated pose (seeds holding one) | Seeds holding a near-native mode | All modes written | RMSD of the top pose (Å) | Zn–N (Å); seeds within 2.3 Å | Thr199 OG1···N (Å); seeds within 3.3 Å | Modes with a heavy atom within 2.4 Å of a water site |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| Vina, ETKDG (validation) | as reported | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 14 (9) → 0 (0) | 0 → 0 | 200 → 200 | 1.10 ± 0.01 → 1.10 ± 0.01 | 2.84 → 2.93; 0 → 0 of 10 | 4.76 → 4.81; 0 → 0 of 10 | 167 → 0 |
|  | carboxylate | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 19 (10) → 0 (0) | 0 → 2 | 200 → 200 | 1.09 ± 0.01 → 1.10 ± 0.01 | 2.83 → 2.90; 0 → 0 of 10 | 4.77 → 4.80; 0 → 0 of 10 | 163 → 0 |
| Vinardo, ETKDG (validation) | as reported | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 23 (10) → 0 (0) | 1 → 5 | 200 → 199 | 1.13 ± 0.00 → 1.12 ± 0.00 | 2.48 → 2.50; 0 → 0 of 10 | 4.45 → 4.47; 0 → 0 of 10 | 152 → 0 |
|  | carboxylate | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 19 (10) → 0 (0) | 1 → 2 | 199 → 200 | 1.12 ± 0.00 → 1.12 ± 0.01 | 2.48 → 2.50; 0 → 0 of 10 | 4.45 → 4.47; 0 → 0 of 10 | 153 → 0 |
| Vina, DFT geometry (production) | as reported | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 16 (10) → 0 (0) | 1 → 0 | 200 → 198 | 1.08 ± 0.00 → 1.08 ± 0.00 | 2.75 → 2.75; 0 → 0 of 10 | 4.68 → 4.67; 0 → 0 of 10 | 149 → 0 |
|  | carboxylate | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 14 (10) → 0 (0) | 0 → 0 | 200 → 199 | 1.08 ± 0.00 → 1.09 ± 0.00 | 2.75 → 2.74; 0 → 0 of 10 | 4.68 → 4.66; 0 → 0 of 10 | 154 → 0 |
| Vinardo, DFT geometry (production) | as reported | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 14 (9) → 0 (0) | 0 → 1 | 200 → 200 | 1.11 ± 0.00 → 1.10 ± 0.00 | 2.40 → 2.43; 0 → 0 of 10 | 4.35 → 4.37; 0 → 0 of 10 | 154 → 0 |
|  | carboxylate | 0 / 0 / 0 / 10 → 0 / 0 / 0 / 10 | 20 (10) → 0 (0) | 1 → 0 | 200 → 200 | 1.11 ± 0.00 → 1.11 ± 0.00 | 2.40 → 2.42; 0 → 0 of 10 | 4.35 → 4.37; 0 → 0 of 10 | 154 → 0 |

What changes with the waters: the rotated pose leaves the mode lists (14 to 23 modes in 9 or 10 seeds
of a dry run, none in any water run), and no mode keeps a heavy atom within 2.4 Å of a water site (149
to 167 modes per dry run did). What does not change: the top pose. In all eight dry runs and in all
eight water runs the top pose of every seed lies between 1.0 and 2.0 Å from the crystal ligand, at a
mean RMSD of 1.08 to 1.13 Å, with the sulfonamide nitrogen at a mean of 2.40 to 2.93 Å from the zinc,
beyond the 2.3 Å cutoff in every seed, and 4.35 to 4.81 Å from Thr199 OG1, beyond 3.3 Å in every seed.
The validation runs therefore fail the strict gate and the 2.0 Å criterion with the waters as they did
without them (section 4).

The waters did not help these two scores in these runs (sulfanilamide, 10 seeds per run, rule
`hbond`), and the numbers give two reasons. First, the rotated
pose was never their top pose (0 of 10 seeds in each of the eight dry runs), so taking it away removes
nothing that stood in the way of the crystal pose. Second, what fails is the distance of the nitrogen
from the zinc, and the four waters do not bring it within the 2.3 Å cutoff: they lie where the rotated pose puts the ring
(section 1), the nearest of
them, 553, 2.93 Å from the sulfonamide nitrogen of the crystal ligand (section 4), and the mean Zn–N
distance of the top pose changes by −0.01 to +0.09 Å between the dry and the water receptor (table;
the largest change, 2.84 → 2.93 Å, is away from the zinc). A near-native mode is written in 0 to 5 of
the 10 seeds of a water run (0 or 1 of a dry run), never as the top pose.

With Vina and Vinardo the number of modes beyond 2.0 Å rises with the waters (167 to 175 per dry run,
179 to 185 per water run). These modes are not the rotated orientation. The list of 20 modes of a seed
is full, or nearly, in both receptors (198 to 200 modes per run), so a mode that leaves is replaced by
another; and in the water runs every mode beyond 2.0 Å has its sulfonamide nitrogen more than 5 Å from
the zinc (124 to 148 of them in a dry run; counted here from `comparison_water.json` → `modes` →
`zn_n`). For these two scores the count of modes beyond 2.0 Å is therefore not a measure of the rotated
orientation; the column "rotated pose" is. The score of the top pose is more negative with the waters
in each of the eight Vina and Vinardo runs (table "Does the crystal pose still dock").

**The stored rotated poses in the water maps.** The best mode of the rotated pose of each seed of the
two dry production runs, scored without a search in the maps with the waters, gets +5346 to +7644 kcal mol⁻¹
over the 5 variants (in the dry runs the best of these modes scored −8.823 kcal mol⁻¹ with Glu106 as
reported and −8.897 with the carboxylate): a hard
overlap with the water oxygens, under every hydrogen rule (table "Stored modes of the dry runs scored
in the maps with the waters"). This rescoring was made for the key run only.

## 6. Dependence on the hydrogen rule

At fixed poses the answer does not depend on the rule. The 20 stored modes of the rotated pose of the dry production run (Glu106 as reported) score +6000 to +7432 kcal mol⁻¹ in the maps of `hbond`, +5846 to +7209 in `rev`, +5346 to +6858 in `bare` and +6410 to +7644 in `fixed`; the 20 near-native modes change by −0.089 (−0.100 to −0.070), −0.090 (−0.102 to −0.073), −0.339 (−0.358 to −0.313) and −0.087 (−0.103 to −0.068) kcal mol⁻¹ (`<variant>/rescore_fixed_poses.json`). The overlap is with the oxygens, which are at the same place under every rule; the hydrogens move the near-native score by a few thousandths between `hbond`, `rev` and `fixed`, and only the bare oxygens, with their uncompensated −1.64 e, shift it by more.

The search gives the same answer under the three rules that were docked. The key run (sulfanilamide
from the DFT geometry, AutoDock4Zn, 20 seeds, Glu106 as reported) was docked in `hbond`, `rev` and
`bare` (`<variant>/comparison_water.json` → `runs` → `sfa_neutral__ad4zn__prod` → `water`; the dry
column is the stored run):

| | dry: `T3_6RL9` | `hbond`: `T3_6RL9_wat4` | `rev`: `T3_6RL9_wat4rev` | `bare`: `T3_6RL9_wat4bare` |
|:--|:--|:--|:--|:--|
| seeds whose top pose is near-native | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 |
| modes of the rotated pose | 40 | 0 | 0 | 0 |
| modes beyond 2.0 Å (the paper's ring-rotated) | 62 | 1 | 1 | 0 |
| of them within 2 kcal mol⁻¹ of the best mode of their seed | 57 | 0 | 0 | 0 |
| modes between 1.0 and 2.0 Å | 35 | 11 | 7 | 6 |
| near-native modes | 41 | 39 | 39 | 41 |
| all modes written | 138 | 51 | 47 | 47 |
| RMSD of the top pose (Å) | 0.61 ± 0.04 | 0.51 ± 0.04 | 0.51 ± 0.03 | 0.52 ± 0.04 |
| score of the top pose (kcal mol⁻¹) | −8.836 ± 0.004 | −9.232 ± 0.004 | −9.242 ± 0.003 | −9.538 ± 0.003 |
| best near-native mode (kcal mol⁻¹) | −8.843 | −9.238 | −9.249 | −9.542 |
| Zn–N of the top pose (Å); seeds within 2.3 Å | 2.20; 20 of 20 | 2.16; 20 of 20 | 2.16; 20 of 20 | 2.12; 20 of 20 |
| Thr199 OG1···N of the top pose (Å); seeds within 3.3 Å | 2.55; 20 of 20 | 2.66; 20 of 20 | 2.65; 20 of 20 | 2.70; 20 of 20 |
| modes with a heavy atom within 2.4 Å of a water site | 84 | 0 | 0 | 0 |
| closest heavy atom of any mode to a water site (Å) | 0.07 | 2.79 | 2.79 | 2.51 |
| largest score span of a seed's list (kcal mol⁻¹) | 2.72 | 2.31 | 2.74 | 2.44 |

Under all three rules no mode of the rotated pose is written, the top pose is near-native in 20 of 20
seeds with the nitrogen on the zinc and on Thr199 in every seed, and no mode has a heavy atom closer
than 2.51 Å to a water site. `hbond` and `rev`, the two directions of the ring of hydrogen bonds, agree
to 0.010 kcal mol⁻¹ in the mean score of the top pose and to 0.01 Å in its RMSD and in its two
contacts. The bare oxygens give a more negative score (−9.538 ± 0.003 kcal mol⁻¹) and a Zn–N distance
of 2.12 Å, as they gave the more negative score at fixed poses. Under `rev` one mode beyond 2.0 Å
remains, as under `hbond`: seed 2, mode 3, 2.36 Å from the crystal ligand, 4.81 Å from the rotated
reference, 2.743 kcal mol⁻¹ behind the top mode of its seed, Zn–N 1.87 Å, aniline N 9.39 Å from Asn67
OD1. It is not the rotated pose. Under `bare` none remains.

Limits of this section. `fixed` was rescored at fixed poses only and was not docked. `rev` and `bare`
were docked for the key run with Glu106 as reported only: the carboxylate receptor, the validation
runs, the sulfonamidate, Vina and Vinardo were run under `hbond` alone.

## 7. The same sites in 3HS4

The four sites were carried into 3HS4 by the CA superposition that the cross-docking analysis uses
(257 atoms, 0.328 Å) and compared with the deposited 3HS4 (`T3_6RL9_wat4/sites_in_T4_3HS4.json`):

| Site (6RL9 water) | Nearest water of 3HS4 (Å) | Nearest other deposited atom of 3HS4 (Å) | Nearest atom of acetazolamide (Å) |
|:--|:--|:--|:--|
| HOH423 | HOH1013, 4.36 | GOL400 O2, 0.21 | AZM C2, 5.62 |
| HOH516 | HOH1216, 2.28 | GOL400 C3, 1.09 | AZM N5, 5.15 |
| HOH553 | HOH1013, 3.83 | GOL400 C3, 1.92 | AZM N7, 2.85 |
| HOH568 | HOH1342, 4.11 | GOL400 C1, 1.47 | AZM C2, 2.98 |

None of the four sites holds a water in 3HS4. They are taken by a glycerol molecule of the
crystallisation (GOL 400, full occupancy), which the receptor preparation removes with the other
additives. A variant `T4_3HS4_wat4` built from 3HS4's own waters therefore does not exist, and none was
made; putting the 6RL9 waters into 3HS4 would be a different experiment. What 3HS4 adds is an
observation: in that crystal the water of these sites has been displaced by another small molecule.

## 8. Conclusion

With waters 423, 516, 553 and 568 kept in the 6RL9 receptor the ring-rotated orientation is gone from the AutoDock4Zn results. This holds in every AutoDock4Zn run that was made, eight in all: the key run of sulfanilamide (DFT geometry, 20 seeds) for both protonation states of Glu106 and, for Glu106 as reported, under three hydrogen rules; the validation run of sulfanilamide (ETKDG conformers, 10 seeds) for both states; and the sulfonamidate (ETKDG conformers, 10 seeds) for both states. In none of the eight is a mode of the rotated pose written (40 or 44 modes in the dry key runs, 20 in each dry run of 10 seeds), and in all eight the top pose is near-native in every seed. Dry → water, for Glu106 as reported and for the carboxylate: key run, top pose near-native in 20 → 20 of 20 and 15 → 20 of 20 seeds, RMSD of the top pose 0.61 ± 0.04 → 0.51 ± 0.04 Å and 1.34 ± 1.30 → 0.51 ± 0.03 Å; validation run, 5 → 10 of 10 and 3 → 10 of 10 seeds, 2.07 ± 1.45 → 0.58 ± 0.03 Å and 2.70 ± 1.38 → 0.60 ± 0.03 Å; sulfonamidate, 3 → 10 of 10 and 7 → 10 of 10 seeds, 1.12 ± 0.32 → 0.61 ± 0.04 Å and 0.86 ± 0.33 → 0.58 ± 0.04 Å. The crystal pose still docks, more closely than without the waters, with the sulfonamide nitrogen on the zinc and on Thr199 in every seed. The stored rotated poses of the key run score +5346 to +7644 kcal mol⁻¹ in the water maps under every hydrogen rule. The differences between the two states of Glu106 that the dry receptors show (key run, 20 against 15 of 20 seeds near-native, `GLU106.md`; validation run, 5 against 3 of 10; sulfonamidate, 3 against 7 of 10) are not there with the waters. By the paper's literal definition (any mode beyond 2.0 Å) one such mode is still written in three of the eight runs, the key run in `hbond` for both states of Glu106 and in `rev`: it is 2.31, 1.48 and 2.74 kcal mol⁻¹ behind the top mode of its seed and is not the rotated pose (4.6, 4.9 and 4.8 Å from it; sections 5 and 6).

The validation runs with AutoDock4Zn pass the paper's strict gate in the water receptors, where the dry run of sulfanilamide fails it for both states of Glu106 and the dry run of the sulfonamidate fails it for Glu106 as reported; PoseBusters physics is passed by every top pose, in the dry and in the water receptors. Vina and Vinardo, which were run with the waters for sulfanilamide only (10 seeds per run, rule `hbond`), fail the strict gate and the 2.0 Å criterion with the waters as without them, in the validation runs of both states. The rotated pose leaves their mode lists too (14 to 23 modes per dry run, none with the waters), but it was never their top pose, and their top pose keeps the sulfonamide nitrogen at a mean of 2.40 to 2.93 Å from the zinc, beyond the 2.3 Å cutoff in every seed, with or without the waters. The waters remove a competitor of the crystal pose. They do not correct a score that places the nitrogen too far from the zinc.

What the result shows and what it does not. A rigid receptor with four waters fixed in place forbids
their displacement by construction: a pose that needs their space scores thousands of kcal mol⁻¹ and
cannot appear. The test therefore shows this: with the four ordered waters that coexist with the bound
ligand in the 6RL9 crystal held in place together, the ring-rotated orientation cannot be docked and
the crystal-like orientation can. The four were always kept together, so the test does not say which
of them, or how many, the rotated orientation has to displace: in the dry key runs 38 of 40 (Glu106 as
reported) and 39 of 44 (carboxylate) modes of the rotated pose have a heavy atom within 2.4 Å of all
four sites, and the other 2 and 5 modes within 2.4 Å of two or three of them (counted here from
`comparison_water.json` → `dry` → `modes` → `water_sites_within_2p4A`). The test does not give the
free-energy cost of displacing them, and so
it does not show that the rotated orientation is impossible in solution; ligands do displace water,
and in 3HS4 a glycerol molecule sits on these very sites. It is also not the protocol of the paper,
which keeps water only on Mg²⁺: a receptor with the waters of the complex's own crystal has been told
part of the answer, as the gyrase B control was when its ordered water was put back
(`dock/CONTROL.md`, section 3.1). This applies to the gate verdicts of section 4 as well: a pass in
this receptor shows which poses made the dry run fail, not that the protocol is validated. For
sulfanilamide that pose is the rotated pose. For the sulfonamidate with Glu106 as reported it is a pose
1.27 to 1.37 Å from the crystal ligand, the top pose of 7 of 10 dry seeds, which is the top pose of no
seed with the waters. That pose has no heavy atom within 2.4 Å of a water site (the nearest is 2.48 to
2.66 Å away; `comparison_water.json` → `sfa_anion__ad4zn` → `dry` → `modes`), so by the criterion used
here it does not overlap a water oxygen as the rotated pose does, and why it loses the top rank was not
measured (the fixed-pose rescoring covers the key run only). The sulfonamidate was docked under the
rule `hbond` only, so whether this second effect depends on the placement of the water hydrogens was
not tested. The numbers here are a diagnostic of the degeneracy, not validation results.

## 9. Not tested

- The hydrogen rule `fixed` in a search: it was rescored at fixed poses only. The rules `rev` and `bare` were docked for the key run with Glu106 as reported only, not for the carboxylate receptor and not for the validation runs, the sulfonamidate, Vina or Vinardo (section 6).
- The sulfonamidate from a DFT geometry, and the sulfonamidate with Vina or Vinardo. With the waters it was run only with AutoDock4Zn from ETKDG conformers, 10 seeds.
- More seeds than those of the stored runs: 20 for the key run, 10 for every other run.
- PoseBusters on the production runs. It was run on the validation runs only, so the production runs have no gate verdict (section 4).
- The fixed-pose rescoring for any run other than the key run.
- 3HS4: no water variant exists (section 7), so the 3HS4 cross-dock, where the DFT geometry gives the rotated pose in 20 of 20 seeds (`dock/GLU106.md`, Glu106 as reported), was not tested with water.
- Subsets of the four waters. Every variant holds all four; no run and no rescoring kept one, two or three of them.
- Other waters. Only the four waters under the rotated pose were added. The 212 other water sites of chain A (216 in the entry, 215 of them kept in the molecular-dynamics system, `md/SYSTEM.md` section 4), some of them near the crystal pose (582 is 3.80 Å from the aniline N, `md/system/audit/waters.csv`), were not; a receptor with all ordered water is a different test.
- Any treatment in which the waters can leave: displaceable or scored waters, flexible receptor, or simulation. That, not this test, is what gives the cost of the rotated orientation; the molecular-dynamics study starts from both orientations for this reason.
- Acetazolamide and the other species in the water receptors; PLIP/ProLIF contacts; torsional strain.

The runs after the two key runs were made by the same driver in the stages `val`, `rule2`, `anion` and `vina`, one job at a time at nice 15 with two threads, the driver waiting up to 20 minutes before a job while the 5-minute load was above 8; the comparison and the tables followed the last job (times in section 10).

## 10. How to reproduce, and how long it took

```bash
export SCRATCH=<a scratch directory>
bash dock/run_water4_sensitivity.sh variants maps      # needs the autogrid4 of dock/GLU106.md
bash dock/run_water4_sensitivity.sh key                 # the two key runs
bash dock/run_water4_sensitivity.sh val rule2 anion vina
bash dock/run_water4_sensitivity.sh rescore compare     # rescore_fixed_poses.json, comparison_water.json
/root/cqc-sfa-tools/dock/bin/python dock/water_variant.py --tables T3_6RL9_wat4 T3_6RL9_glu106GLU_wat4 \
    T3_6RL9_wat4rev T3_6RL9_wat4bare T3_6RL9_wat4fixed  # the tables of section 5
```

The stored `rescore_fixed_poses.json` files hold the key run only. The `rescore` stage as written also rescores the stored validation run `sfa_neutral__ad4zn`, which is the default of `water_variant.py --rescore`; `water_variant.py --rescore --scratch $SCRATCH/rescore --only sfa_neutral__ad4zn__prod` restricts it to the key run, as in the stored files. The log holds no output of the `rescore` stage, only the note that the rescoring was done in two processes.

Wall time (`dock/T3_6RL9_wat4/run_water4.log`, `comparison_water.json`; all times on 6 October 2026, IST): the five variants were written 1:05 to 1:07 AM (file times of `variant.json`), their maps made 1:08 to 1:10 AM (log) and the last map check written at 1:11 AM (file time of `T3_6RL9_wat4rev/ad4zn/maps_report.json`); the fixed-pose rescoring 8 min on two threads (2:27–2:35 AM); the 16 docking jobs 1:12 AM to 11:38 AM, as in the table; the final comparison of the 16 runs, with PoseBusters on the eight validation runs, and the tables 11:38 to 11:41 AM. Everything at `--cpu 2`, nice 15, one job at a time.

| Run | Seeds | Started | Finished | Wall time (s) | Wait before the job (min; 5-minute load at the start of the wait → at the start of the job) |
|:--|--:|:--|:--|--:|:--|
| `T3_6RL9_wat4` `sfa_neutral__ad4zn__prod` | 20 | 1:12 AM | 2:26 AM | 4466 | 1; 8.13 → 7.97 |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__ad4zn__prod` | 20 | 2:46 AM | 4:43 AM | 7015 | 10; 10.35 → 8.00 |
| `T3_6RL9_wat4` `sfa_neutral__ad4zn` | 10 | 5:03 AM | 6:06 AM | 3797 | 20; 11.37 → 9.40 |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__ad4zn` | 10 | 6:26 AM | 7:19 AM | 3215 | 20; 12.14 → 10.37 |
| `T3_6RL9_wat4rev` `sfa_neutral__ad4zn__prod` | 20 | 7:28 AM | 8:39 AM | 4233 | 9; 10.87 → 7.89 |
| `T3_6RL9_wat4bare` `sfa_neutral__ad4zn__prod` | 20 | 8:41 AM | 9:19 AM | 2285 | 2; 8.57 → 7.68 |
| `T3_6RL9_wat4` `sfa_anion__ad4zn` | 10 | 9:19 AM | 9:29 AM | 576 | none |
| `T3_6RL9_glu106GLU_wat4` `sfa_anion__ad4zn` | 10 | 9:29 AM | 9:46 AM | 1033 | none |
| `T3_6RL9_wat4` `sfa_neutral__vina__prod` | 10 | 9:46 AM | 9:52 AM | 349 | none |
| `T3_6RL9_wat4` `sfa_neutral__vinardo__prod` | 10 | 9:52 AM | 9:59 AM | 446 | none |
| `T3_6RL9_wat4` `sfa_neutral__vina` | 10 | 9:59 AM | 10:06 AM | 385 | none |
| `T3_6RL9_wat4` `sfa_neutral__vinardo` | 10 | 10:06 AM | 10:25 AM | 1158 | none |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__vina__prod` | 10 | 10:45 AM | 11:20 AM | 2112 | 20; 9.86 → 9.75 |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__vinardo__prod` | 10 | 11:24 AM | 11:30 AM | 327 | 4; 8.66 → 7.66 |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__vina` | 10 | 11:30 AM | 11:33 AM | 230 | none |
| `T3_6RL9_glu106GLU_wat4` `sfa_neutral__vinardo` | 10 | 11:33 AM | 11:38 AM | 264 | none |

The 16 jobs add up to 31891 s (8.86 h) of docking, and the driver waited 86 min in all before eight of them. These times depend on the load of the machine and are not a measure of the cost of the waters. A seed took 3.7 and 5.8 minutes in the two key runs, against 34 s and 92 s in the stored dry runs they are compared with (`comparison_water.json` → `dry` → `wall_s_total`), and the same Vina production run took 349 s in one water receptor and 2112 s in the other. Other jobs held the machine (during the key runs a quantum-chemistry job and a molecular-dynamics run): the 5-minute load read before a job was above 8 eight times (8.13 to 12.14), and three jobs started after the full wait of 20 minutes with the load still above 8 (9.40, 10.37 and 9.75).
