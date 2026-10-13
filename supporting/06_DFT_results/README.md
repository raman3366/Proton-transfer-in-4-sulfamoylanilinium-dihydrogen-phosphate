# 06 DFT results

All molecular calculations used PySCF 2.14.0 (density fitting) with geomeTRIC 1.1.1. The primary level is
**B3LYP-D3(BJ)/6-311++G(d,p)**; IEF-PCM means water (*ε* = 78.3553), Lebedev order 29. The periodic model used
Quantum ESPRESSO 7.5. The wavefunction analyses used Multiwfn 2026.9.20. Table S1 of the supplementary
material (05_Tables/SI/Table_S01_settings) lists every calculation. Tables and figures
numbered with an S are those of the supplementary material.

The calculations of VM burst 3 (facts sections 32–35) are: the Hessian of the neutral
IEF-PCM minimum (`sfa_pair_pcm_d3_oform/`), the constrained periodic N1–H1A profile and the single-transfer control
(`periodic_QE/profile.json`, `relax_prof*`, `relax_iir1`), the five-pair cluster and its nine-bond QTAIM
(`crystal_cluster_qtaim/`, `sfa_xclu_salt/`, `sfa_xclu_neutral/`) and the relaxed torsion scans
(`sulfanilamide_pcm_d3/torsion_scan/`, `sulfonamidate_pcm_d3/torsion_scan/`, `torsion_scan/`).

## One folder per case

| Folder | Model (Table 1) | Level | Environment | Used for | Contents |
|---|---|---|---|---|---|
| `sfa_pair_pcm_d3` | M1 ion pair | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | primary model: salt minimum; frequencies, FMO, TD-DFT, NLO, ΔSCF, interaction energy, QTAIM/IGMH/NCI/ELF/LOL, MEP, charges | opt, freq+IR, Hessian, FMO, TD-DFT, ΔSCF, NLO, interaction, Multiwfn, molden |
| `sfa_pair_gas_d3` | M1 ion pair | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | optimises to the neutral sulfanilamide···H3PO4 adduct (proton on O5) | opt, FMO, interaction, Multiwfn, molden |
| `sfa_pair_pcm_b3lyp` | M1 ion pair | B3LYP/6-311++G(d,p), no dispersion | IEF-PCM (water) | dispersion sensitivity check (geometry) | opt, molden (converted from the checkpoint, see "Wavefunctions") |
| `sfa_pair_pcm_d3_tz` | M1 ion pair | B3LYP-D3(BJ)/def2-TZVP | IEF-PCM (water) | basis-set check (geometry) | opt, molden (converted from the checkpoint, see "Wavefunctions") |
| `sfa_pair_pcm_d3_oform` | M1 ion pair, neutral minimum (H1A on O5) | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | the O5-protonated (neutral) minimum of the scan (chain B end point): Hessian, thermochemistry and IR (VM burst 3) | opt, freq+IR, Hessian, molden (converted from the checkpoint, see "Wavefunctions") |
| `cation_pcm_d3` | M2 cation | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | fragment: interaction energy, TD-DFT, descriptors, Fukui, ADMET reference; docking sensitivity species | opt, freq+IR, Hessian, FMO, TD-DFT, ΔSCF, Multiwfn, rigid dock strain (retired), molden |
| `cation_gas_d3` | M2 cation | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | fragment | opt, freq+IR, Hessian, FMO, ΔSCF, Multiwfn, molden |
| `h2po4_pcm_d3` | M3 dihydrogen phosphate | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | fragment: interaction energy | opt, freq+IR, Hessian, FMO, Multiwfn, molden |
| `h2po4_gas_d3` | M3 dihydrogen phosphate | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | fragment | opt, freq+IR, Hessian, FMO, Multiwfn, molden |
| `sulfanilamide_pcm_d3` | M4 sulfanilamide | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | neutral species at pH 7.4: primary docking ligand, proton-transfer energetics, Fukui, ADMET, relaxed C3–C4–S1–N2 torsion scan (`torsion_scan/`, VM burst 3) | opt, FMO, ΔSCF, NLO, TD-DFT, Multiwfn, relaxed torsion scan, rigid dock strain (retired), molden |
| `sulfanilamide_gas_d3` | M4 sulfanilamide | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | proton-transfer energetics | opt, freq+IR, Hessian, FMO, ΔSCF, NLO, Multiwfn, molden |
| `sulfanilamide_pcm_d3_tz` | M4 sulfanilamide | B3LYP-D3(BJ)/def2-TZVP | IEF-PCM (water) | basis-set check | opt, molden (converted from the checkpoint, see "Wavefunctions") |
| `sulfonamidate_pcm_d3` | M5 sulfonamidate | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | docking species for the zinc site; relaxed C3–C4–S1–N2 torsion scan (`torsion_scan/`, VM burst 3) | opt, FMO, Multiwfn, relaxed torsion scan, rigid dock strain (retired), molden |
| `h3po4_pcm_d3` | M6 phosphoric acid | B3LYP-D3(BJ)/6-311++G(d,p) | IEF-PCM (water) | proton-transfer energetics | opt, FMO, Multiwfn, molden |
| `h3po4_gas_d3` | M6 phosphoric acid | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | proton-transfer energetics | opt, freq+IR, Hessian, FMO, Multiwfn, molden |
| `urea_gas_d3` | M7 urea | B3LYP-D3(BJ)/6-311++G(d,p) | gas phase | NLO reference | opt, freq+IR, Hessian, FMO, NLO, Multiwfn, molden |

Special folders:

| Folder | What |
|---|---|
| `proton_scan/` | Relaxed scan of *d*(N1–H1A) in IEF-PCM (water) on the ion pair: chain A from the salt minimum, chain B from the neutral minimum (`A_*`, `B_*`: `.xyz` geometries, `.json` energies and convergence, `.partial.xyz` last steps, `*_traj.xyz` trajectories, `.constraints` geomeTRIC inputs); `scan_profile.csv`/`.png`, `scan_summary.json`; `zpe_barrier.json`/`.png` (zero-point energy of the N–H stretch against the barrier). |
| `eps_ladder/` | The ion pair optimised in IEF-PCM with *ε* = 4, 8 and 20 (`eps_*_N/`: start, optimised and partial geometries, trajectory, `result.json`); summary `eps_ladder.csv`/`.json`/`.png`. |
| `functional_check/` | Single points with ωB97X-D, M06-2X, PBE0-D3(BJ), DF-MP2 and B3LYP-D3(BJ) (IEF-PCM water, 6-311++G(d,p)) on three B3LYP-D3(BJ) scan geometries: the salt end point, the 1.25 Å point and the neutral end point (`fsp_scan.json`, `fsp_scan.md`; `fsp_scan_partial.json` is the run's checkpointed partial output). |
| `periodic_QE/` | The periodic *Pc* cell (54 atoms, cell fixed at the X-ray values). `periodic_qe.json` holds the relaxations, the comparison and the hybrid correction of the single configuration (`relax`, `comparison`, `comparison.constrained_neutral_iirfix`, `hybrid`); `profile.json` holds the constrained N1–H1A profile with its PBE0-D3(BJ) single points and the single-transfer control (`points`, `hybrid`, `neutral_minimum`, `single_transfer_control`; VM burst 3); `hybrid_estimate.json` the cost estimate made before the hybrid runs (a record of the planning step; its own status string calls it an order-of-magnitude estimate, and the hybrid single points in `periodic_qe.json` and `profile.json` replaced it). `runs/<run>/<run>.in` and `.out` are the pw.x inputs and outputs of every run: the relaxations (`relax_*`, with `_pre` coarse pre-relaxations), the constrained profile points (`relax_prof*`), the convergence tests (`scf_*`) and the PBE and PBE0 single points (`pbe_*`, `hybrid_*`, `hybrid_prof*`); `runs.jsonl` is the run log, `convergence.json` the cutoff and k-grid tests, `build.json` the cell construction. Structures as `.extxyz` and, converted here, as `cif/*.cif` (P1, all 54 atoms; the comment lines of each CIF say what the structure is). |
| `crystal_cluster_qtaim/` | QTAIM of the nine crystal hydrogen bonds: the salt cluster against the cluster with H1A on O5, bond by bond (`qtaim_nine_bonds.json/.csv/.md`), and how the clusters were cut (`build.json`) (VM burst 3). |
| `sfa_xclu_salt/` | Five-pair cluster (M8) cut from the relaxed periodic salt, relax_i: B3LYP-D3(BJ)/6-31+G(d,p) gas-phase single point (`scf.json`, `final.molden`), Multiwfn QTAIM (`mwfn_extra/`: `CPs.txt`, `CPprop.txt`, `paths.txt`, `qtaim*.json/csv`, keystrokes in `inputs/`, output in `logs/`) and the nine bond critical points (`nine_bcps.json`). |
| `sfa_xclu_neutral/` | The same cluster cut from the periodic cell with H1A held on O5 (relax_iirfix, O5–H1A 0.98 Å); same files. |
| `torsion_scan/` | Relaxed torsional strain of sulfanilamide and the sulfonamidate at the crystal-ligand torsions and at the docked torsions (`torsion_strain.json/.csv/.md`, `torsion_profile.png`), from the two scans in `sulfanilamide_pcm_d3/torsion_scan/` and `sulfonamidate_pcm_d3/torsion_scan/` (VM burst 3). |
| `hirshfeld_surface/` | Hirshfeld surfaces of the cation and the anion in the X-ray structure (neutron-normalised X–H): contacts, red spots, neighbours (`*.csv`), fingerprint and surface arrays (`*.npz`), surface vertices (`*_vertices.csv.gz`), the cluster used (`cluster.xyz`, `.json`), `hirshfeld_summary.json`, and the Multiwfn cross-check (`mwfn_*/`, `mwfn_*.json`, `hirshfeld_compare.json`). |

Periodic model, what to cite: `relax_i` is the salt (from the X-ray start); `relax_ii` started with H1A on O5 in both
formula units and `relax_iir1` in one only, and both relaxed back to the salt. The periodic cost of moving the proton is
the constrained N1–H1A profile `relax_prof110` … `relax_prof170` with its PBE0-D3(BJ) single points (`profile.json`;
facts section 32). `relax_iirfix` is one constrained configuration (H1A and O5 held, O5–H1A 0.98 Å), cited only as such
(facts section 36). **`relax_iifix` is withdrawn** (its
lattice slid back to a salt); its files are kept only for provenance. The QE potentials are not included (SSSP
efficiency 1.3.0 and SG15 ONCV v1.2; file names and md5 sums in `periodic_qe.json` → `method`); the runs kept no
`save/` wavefunction directories. Note that `comparison.verdict` in `periodic_qe.json` says "single minimum"; the paper's
claim is only that no neutral minimum was found.

## File types

| File | Content |
|---|---|
| `opt.xyz` | Optimised geometry (Å); the comment line gives the level and the energy. `preopt.xyz`: B3LYP/6-31G(d) pre-optimisation. `opt_traj.xyz`, `geomeTRIC_optim_backup.xyz`: optimisation trajectory. |
| `opt_summary.json`, `status.json`, `STAGES.md`, `proton_check.json`, `geometry_vs_start.csv` | Convergence, run log, where the proton ended, geometry against the start. `status.json` and `STAGES.md` are the pipeline's own log, written during the run and kept as written: an optimisation whose job was stopped (exit code −15) is logged there as failed and was accepted on negligible forces, and the line after every optimisation names the frequency stage, which decides whether the point is a minimum. Table S1 gives the final forces, and Table S2 which structures have a Hessian (the others are stationary points on the force criterion only). |
| `freq.json` | Harmonic frequencies, IR intensities, zero-point energy, H and G at 298 K. `hessian.npy`, `dipole_derivatives.npy`: the raw arrays (NumPy). |
| `ir_sticks.csv`, `ir_spectrum.csv` | Made here from `freq.json`: every mode (harmonic and scaled by 0.9688, the factor of Figure S1) and a Lorentzian spectrum (FWHM 10 cm⁻¹, 400–4000 cm⁻¹) as in Figure S1. Harmonic values for a single molecule or ion pair; see the Figure S1 caption before comparing with a measured solid-state spectrum. |
| `fmo.json`, `fmo_table.md` | Orbital energies, HOMO/LUMO, Koopmans-type descriptors, dipole, fragment composition of the frontier orbitals. |
| `dscf.json`, `.md` | Vertical ΔSCF ionisation energy, electron affinity and descriptors (used in the paper instead of the Koopmans values). |
| `tddft_<functional>_<eq/noneq>.json`, `.md`; `spectrum_*.csv` | TD-DFT excitations and the broadened UV-Vis spectrum; `nto_*.npz`: natural transition orbitals (NumPy). |
| `nlo.json`, `.md` (`nlo_geom-*`: at another case's geometry) | Finite-field dipole moment, polarisability and first hyperpolarisability. |
| `interaction.json`, `.md` | Counterpoise-corrected interaction and deformation energies. |
| `mwfn_extra/` | Multiwfn analyses: `charges.*`, `qtaim*.json/csv` (bond critical points; `CPs.txt`, `CPprop.txt` raw), `igmh*`, `dos*`, `mep_*`, `mo_composition.csv`, `summary.json`; keystroke inputs in `inputs/`, outputs in `logs/`. |
| `mep_extrema.csv`, `mep_surface.json`, `mep_range.json`, `surfanalysis.txt`, `vtx.txt` | Electrostatic potential on the 0.001 a.u. surface: extrema, statistics, raw Multiwfn output and surface vertices. |
| `nci_spikes.csv`, `rdg_scatter.npz`, `elf_*.npz`, `lol_*.npz` | NCI spikes and the plotted RDG scatter, ELF and LOL planes (NumPy archives behind Figures S2 and S3). |
| `final.molden` | The SCF wavefunction at the optimised geometry (from `dft/fmo.py`). It opens in Multiwfn, Avogadro, Jmol, IboView or Molden. |
| `mwfn_*.in`, `mwfn_home/settings.ini` | Multiwfn keystroke inputs and the settings used. |
| `quicklook_renders/` | Automatic renders made by the pipeline (orbitals, MEP, NCI, ELF/LOL); these are not the paper figures. |
| `torsion_scan/` | Relaxed scan of the C3–C4–S1–N2 torsion (sulfanilamide and sulfonamidate cases): one constrained optimisation per angle (`tau<angle×100>.json/.xyz`, `.partial.xyz`, `*_traj.xyz`), the unconstrained reference (`ref.*`), `torsion_scan.csv` and `torsion_scan_summary.json`. |
| `superseded/dock_strain.json`, `.md`, `dock_strain_partial.json` | **Superseded:** rigid DFT single points of the docked poses (burst 2), replaced by the relaxed torsion scans (`<case>/torsion_scan/`; the strain values are in the top-level `torsion_scan/` folder; facts section 36). Kept in a `superseded/` folder for provenance, with a first line in the `.md` file saying so; not used in the paper. |
| `case_definition.json` | The case as built from the CIF: charge, CIF atom labels, fragments, hydrogen bonds tracked. |
| `logs/` | Run logs. |

## Wavefunctions

`final.molden` (written by `dft/fmo.py` at the optimised geometry) is included for 12 of the 16 cases, and a molden file converted from the stored checkpoint for the other 4 (`sfa_pair_pcm_b3lyp`, `sfa_pair_pcm_d3_tz`, `sfa_pair_pcm_d3_oform`, `sulfanilamide_pcm_d3_tz`). The special folders `sfa_xclu_neutral/`, `sfa_xclu_salt/` carry the `final.molden` of their own single point.
- `sfa_pair_pcm_b3lyp/laststep_from_opt_chk.molden`: converted by the packaging script from the stored checkpoint `opt.chk` (PySCF `molden.from_mo`, no new SCF; 443 of 445 MOs, the rest padded with empty vectors as in `dft/fmo.py` so that Multiwfn reads it); the SCF of the last evaluated optimisation step, E = -1536.05465988 Ha, atoms at most 0.0040 Å from `opt.xyz` (the optimisation was stopped on negligible forces, so the last evaluated step differs slightly from the accepted geometry). Basis 6-311++g(d,p).
- `sfa_pair_pcm_d3_tz/laststep_from_opt_chk.molden`: converted by the packaging script from the stored checkpoint `opt.chk` (PySCF `molden.from_mo`, no new SCF; all 574 MOs); the SCF of the last evaluated optimisation step, E = -1536.28934544 Ha, at the geometry of `opt.xyz`. Basis def2-tzvp.
- `sfa_pair_pcm_d3_oform/final_from_chk.molden`: converted by the packaging script from the stored checkpoint `final.chk` (PySCF `molden.from_mo`, no new SCF; 444 of 445 MOs, the rest padded with empty vectors as in `dft/fmo.py` so that Multiwfn reads it); the SCF of the final single point at the optimised geometry, E = -1536.11151121 Ha, at the geometry of `opt.xyz`. Basis 6-311++g(d,p).
- `sulfanilamide_pcm_d3_tz/laststep_from_opt_chk.molden`: converted by the packaging script from the stored checkpoint `opt.chk` (PySCF `molden.from_mo`, no new SCF; all 395 MOs); the SCF of the last evaluated optimisation step, E = -891.88365832 Ha, at the geometry of `opt.xyz`. Basis def2-tzvp.
- The converted files were not used for any analysis in the paper.

## Not included, and how to regenerate it

| Excluded | Files | Size |
|---|---:|---:|
| cube files | 56 | 485.2 MB |
| rdg_scatter.txt | 2 | 244.6 MB |
| VTK render meshes | 26 | 142.3 MB |
| checkpoints (.chk) | 64 | 140.8 MB |
| ELF/LOL *_plane.txt | 12 | 76.7 MB |
| periodic_qe_sandbox (coarse test) | 9 | 32.9 kB |
| other (lock files, bytecode) | 5 | 22.0 kB |
| **Total** | | **1.1 GB** |

Every excluded file is listed with its size, its location in the Pi archive (`/root/cqc-sfa-archive/vm/cqc-sfa/`) and
the reason in 10_Reproducibility/excluded_files.csv. The archive is sha256-verified (10_Reproducibility/vm_archive_MANIFEST.sha256).

**Cube files** (HOMO, LUMO, density, ESP; NCI and IGMH grids) can be regenerated from `final.molden`:

- In Multiwfn: `Multiwfn final.molden`, main function 5 (grid data), then property 1 (electron density), 4 (an
  orbital; give its index), 12 (total ESP), 9 (ELF) or 10 (LOL), choose the grid, and export the grid as a
  Gaussian-type cube file. The RDG and sign(*λ*~2~)*ρ* cubes come from main function 20 → 1 (NCI) and the IGMH cubes
  from 20 → 11; the exact keystrokes used are in `mwfn_nci.in` and `mwfn_extra/inputs/igmh.in`.
- Or with PySCF as in the pipeline (10_Reproducibility/scripts/dft/fmo.py): `pyscf.tools.cubegen` on an 80³ grid
  with a 4.5 Å margin (ESP on 64³).

**Checkpoints** are PySCF restart files. They are needed only to restart a calculation; every result is in the JSON,
CSV and molden files.
