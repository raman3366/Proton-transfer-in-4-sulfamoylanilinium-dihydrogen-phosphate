---
# Supplementary material (assembled 29 Sep 2026; named and formatted for JICS, Elsevier, 5 Oct 2026; Sections S15 and
# S16 added 11 Oct 2026). Build: tools/paper-env/bin/python paper/_shared/tools/build.py si
# Since 8 Oct 2026 this file is the source and is edited directly (the conversion script apply_jics.py is retired).
# First page: title, the author block of the main text
# and a contents list that crossref.lua writes into the empty {#si-contents} div from the section headings and the
# tables and figures of each section.
# Numbers carry the prefix S (Fig. S1, Table S1, Section S1); numbers without it refer to the main text.
# Tables and the ADMET text are pulled in by <!-- include: ... --> lines (tables/, sections/); every number in them
# comes from paper/main/data/facts.json. Tables are numbered by first appearance, so file names and printed numbers
# differ after Table S14 (VM burst 3, 29 Sep 2026): tableS13b_periodic_profile.md prints as Table S15, tableS15 as S16,
# tableS_fukui as S17, tableS16-S19 as S18-S21, tableS18b_cluster_qtaim as S22, tableS20-S26 as S23-S29,
# tableS26b_torsion_scan as S30, tableS27 and S28 as S31 and S32, the two control-receptor tables
# (tableS29_control_receptors, tableS30_control_scores; 4 Oct 2026) as S33 and S34, tableS_admet as S35, and the five
# tables of the two docking tests and of the simulations (11 Oct 2026), included after it so that no earlier number
# moved: tableS31_dock_glu106_water as S36, tableS32_dock_glu106_other as S37, tableS33_md_settings as S38,
# tableS34_md_runs as S39 and tableS35_md_gates as S40. The main text cites supplementary tables by number (text and
# table notes: S1-S3, S6-S11, S13, S15, S16, S22, S30, S31, S33 and S34); check them after moving any table. It cites
# the two new sections by number (Section S15, Section S16); any table or figure number it quotes from them
# (Tables S36-S40, Figs. S9 and S10) is typed by hand there and must be checked against the order of the includes here.
# Figures: figS_mep prints as Fig. S4, the cluster QTAIM figure as S5, the Hirshfeld decomposition as S6, the
# docking seeds as S7, the control receptor (figS_control_gyrb) as S8, the orientation coordinate of the simulations
# (figS_md_orientation) as S9 and their ligand RMSD (figS_md_rmsd) as S10. Main-text figures quoted here (Figs. 5, 6
# and 8) follow the printed order of the main text (fig09_redock prints as Fig. 8).
# In the text of Sections S1 and S13 the new tables are pointed to by section ("Section S15") and not by label: a
# label used there would be the first mention of a late table and the linter reads first mentions in order.
title: "Supplementary material for: Proton transfer in 4-sulfamoylanilinium dihydrogen phosphate: a combined molecular and periodic DFT study with molecular docking"
author-block: |
  C. Muthuselvi^a^, Anantha Raman Lakshmipathi^b,\*^, D. Saravanakkumar^c^, Venkatramanan Varadharajan^d^

  ^a^ Department of Physics, Sri Ramanas College of Arts and Science for Women, Chidambarapuram, Aruppukottai 626134, Tamil Nadu, India

  ^b^ Department of Mechanical Engineering, Rajalakshmi Institute of Technology, Chennai 600124, Tamil Nadu, India

  ^c^ PG and Research Department of Physics, Thiagarajar College (Autonomous), affiliated to Madurai Kamaraj University, Theppakulam, Madurai 625009, Tamil Nadu, India

  ^d^ Department of Biotechnology, PSG College of Technology, Peelamedu, Coimbatore 641004, Tamil Nadu, India

  \*Corresponding author. E-mail address: raman3366@gmail.com (A.R. Lakshmipathi).

  E-mail addresses and ORCID iDs:\
  C. Muthuselvi: muthuselviganesan100@gmail.com, ORCID 0000-0003-2935-9714\
  A.R. Lakshmipathi: raman3366@gmail.com, ORCID 0000-0003-3465-286X\
  D. Saravanakkumar: saravanakkumar_physf@tcarts.in, ORCID 0009-0006-1437-0266\
  V. Varadharajan: mail4venkat1992@gmail.com, ORCID 0000-0003-1532-7762
# PDF only (crossref.lua): cap figure height at 18.9 cm; the Hirshfeld decomposition figure and its caption overflow
# an A4 page at full size. The Word file keeps full size.
pdf-max-figure-height: 18.9
number-prefix: S
---

<!-- Citation keys used only inside the included files, listed here so that lit.py bib puts them in
     references.json (the bibliography scan does not follow includes): [@merrick2007evaluation; @parr1984density;
     @yang1986use; @reed1985natural; @hirshfeld1977bonded; @wildman1999prediction; @ertl2000fast;
     @lipinski2001experimental; @veber2002molecular; @ghose1999knowledge; @egan2000prediction; @swanson2024admet;
     @kim2025pubchem; @knox2024drugbank; @sorkun2019aqsoldb; @zhu2009quantitative; @hansen2009benchmark;
     @thakkar2020drug; @davis1943binding; @taylor1940comparative; @feyereisen1993use; @weigend2002efficient;
     @greenbaum1955comparison; @klotz1945isoelectric; @zheng2026digitized; @koizumi1964absorption;
     @nogami1968studies] -->

This document holds the data behind the main text on 4-sulfamoylanilinium dihydrogen phosphate (SAPH): settings and
provenance, coordinates, complete tables, the analyses that the main text only summarises, and the docking and ADMET
details. The model names are those of the main text. The ion pair is the 27-atom cation–anion pair of the crystal,
computed in the gas phase, in IEF-PCM (water) and along a ladder of IEF-PCM (*ε* = 4, 8, 20); the periodic model is the
full *Pc* cell with PBE-D3(BJ). Unless a table says otherwise, the molecular level is B3LYP-D3(BJ)/6-311++G(d,p). The
docked species are neutral sulfanilamide, its sulfonamidate and, as a sensitivity case, the cation; no docking run used
the ion pair. The same species are also docked into a control receptor that is not a sulfonamide target (Section S13).
Two further sections report work done after the docking. Section S15 repeats part of the human carbonic anhydrase II
(hCA II) runs with Glu106 as a carboxylate and with four crystal waters kept in the receptor, and Section S16 reports molecular
dynamics of the carbonic anhydrase II–sulfanilamide complex, in which the ligand is the sulfonamidate bonded to the
zinc, with the tests of its force field, among them the two quantum-chemical tests that its terms at the sulfonamide
nitrogen failed, and its limits. Every table was written by script from the stored results.

# Contents

::: {#si-contents}
:::

# S1. Computational settings, software and provenance {#sec:s-settings}

@tbl:s-settings lists every calculation with its level, environment, integration grid and SCF threshold, and the
directory of the data package that holds its inputs and results. The molecular calculations used PySCF 2.14.0
[@sun2020recent; @sun2026python] with the Libxc functional library [@lehtola2018recent] and geomeTRIC 1.1.1
[@wang2016geometry] as optimiser, and Multiwfn 2026.9.20 [@lu2012multiwfn] for the wavefunction analyses. The periodic
model was run with Quantum ESPRESSO 7.5 [@giannozzi2020q], with SSSP efficiency 1.3.0 potentials
[@prandini2018precision] for PBE-D3(BJ) and SG15 ONCV v1.2 potentials [@schlipf2015optimization] for the PBE0-D3(BJ)
single points. The ±1 ions of the ΔSCF and Fukui runs are unrestricted doublets. All other runs are closed-shell.
Docking used AutoDock Vina 1.2.7 [@eberhardt2021autodock] with the Vina and Vinardo functions and AutoDock4Zn
[@santosmartins2014autodock4zn] on AutoGrid 4.2.7 maps. Receptors and ligands were prepared with PDBFixer (OpenMM 8.6.1
[@eastman2024openmm]), PROPKA 3.5.1 [@olsson2011propka3], RDKit 2026.03.1 and Meeko 0.8.0, and the poses were analysed
with spyrmsd 0.9.0 [@meli2020spyrmsd], PoseBusters 0.6.5 [@buttenschoen2024posebusters], PLIP 3.0.1 [@adasme2021plip],
ProLIF 2.2.2 [@bouysset2021prolif] and gemmi 0.7.5 [@wojdyr2022gemmi]. The analyses used NumPy [@harris2020array] and
SciPy [@virtanen2020scipy]; the figures and the graphical abstract were drawn from the stored results with Matplotlib
3.11.2 [@hunter2007matplotlib] and, for the three-dimensional views, PyVista 0.49.0 [@sullivan2019pyvista] (VTK 9.7.0).
The scripts of the calculations, the analyses and the figures, those of this document among them, were written and run
with the generative artificial-intelligence tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors'
direction, as Section 2.6 of the main text states; no image-generation tool was used.
For the ADMET check of Section S14 we used ADMET-AI 1.4.0 [@swanson2024admet] (Chemprop 1.6.1, PyTorch 2.5.0) and RDKit
2026.03.6. The docking tests of Section S15 used the docking programs named above; their AutoDock4Zn maps were made with
an AutoGrid binary built from the official source, which regenerates the stored maps of the unmodified receptors byte
for byte. The molecular dynamics of Section S16 was run with OpenMM 8.6.1 [@eastman2024openmm] on the CUDA platform in
mixed precision, on one NVIDIA L40S graphics processor. Its topologies were built with AmberTools 26.0
[@case2023ambertools] (tleap, antechamber and MCPB.py [@li2016mcpb]) and read back with ParmEd 4.3.1
[@shirts2017lessons]; the quantum chemistry of the zinc site, whose settings are given in Section S16 and not in
@tbl:s-settings, used PySCF 2.14.0 and geomeTRIC 1.1.1 as above; and the trajectories were analysed by two programs
written to one specification, the first on MDAnalysis 2.10.0 [@michaudagrawal2011mdanalysis; @gowers2016mdanalysis] and
the second on NumPy 2.5.3 alone (Python 3.12.14).

The following details complete the Computational methods of the main text. B3LYP [@becke1993density;
@lee1988development] was used in the definition that PySCF shares with Gaussian, with VWN-RPA local correlation; a
VWN5-based B3LYP gives different absolute energies [@hertwig1997parameterization]. For S and P the 6-311G part of
6-311++G(d,p) is the set of McLean and Chandler [@mclean1980contracted]. Of the 445 basis functions of the ion pair,
PySCF removed two near-linear combinations (443 orbitals), and Coulomb and exchange terms were density-fitted (RI-JK)
with the auxiliary basis that PySCF assigns, aug-cc-pVTZ-JKFIT (def2-TZVP-JKFIT with def2-TZVP, heavy-aug-cc-pVDZ-JKFIT
for the five-pair cluster). The IEF-PCM cavity used a SWIG-discretised surface [@lange2010smooth] on Lebedev spheres of
order 29 (302 points each) and PySCF's default cavity of Bondi radii [@bondi1964van] with H at 1.10 Å, all scaled by
1.2. ωB97X-D was run as the Libxc functional HYB_GGA_XC_WB97X_D with the Chai–Head-Gordon damped dispersion term
[@chai2008long] added separately (−0.0163 to −0.0171 Ha at the three geometries of the functional check, Section S5).

Optimisations used geomeTRIC 1.1.1 [@wang2016geometry] with Gaussian's default criteria (maximum and RMS forces 4.5 ×
10⁻⁴ and 3.0 × 10⁻⁴ Ha Bohr⁻¹, plus displacement limits), after a B3LYP/6-31G(d) pre-optimisation (none for def2-TZVP).
Two optimisations ended without the convergence flag. The IEF-PCM (water) ion pair met the force criteria and then crept
along a nearly free inter-ion rotation, in steps too long for the displacement limit, its energy falling by 1.6 × 10⁻⁵
Ha (0.01 kcal mol⁻¹) between steps 51 and 64; the job was stopped at step 64 with RMS and maximum forces of 7.35 × 10⁻⁵
and 1.85 × 10⁻⁴ Ha Bohr⁻¹. The ion pair without dispersion ended alike (5.92 × 10⁻⁵ and 1.21 × 10⁻⁴ Ha Bohr⁻¹). Both
were accepted on their forces, and the Hessian of the primary model has no imaginary frequency. Scan and ladder points
were accepted by a rule fixed in advance: geomeTRIC convergence, or three consecutive gradients either ten times inside
the force thresholds or inside them with an energy change below 10⁻⁶ Ha. Finite differences for the Hessian used a step
of 0.005 Bohr (Hessian asymmetry 0.0006). The top of the relaxed scan, at 1.25 Å, stopped at the step limit with a
constrained maximum gradient of 1.6 × 10⁻⁴ Ha Bohr⁻¹ and a last-step energy change of 0.0008 kcal mol⁻¹, and is kept
with a flag. The B3LYP-D3(BJ) single points reproduce the scan energies, and the DFT single points reproduce an earlier
run on the same geometries within 3.9 × 10⁻¹¹ Ha. The zero-point energy of the barrier comparison is that of the normal
mode with the largest N1–H1A stretch character, checked against the curvature of the scan; the one-dimensional proton
model follows treatments of strong N–H···O and O–H···O bonds [@kaplanskiy2025spectralstructural] and solves the
Schrödinger equation on the relaxed profile in a sinc discrete variable representation [@colbert1992novel]. The profile
is the lower envelope of the two scan chains in a cubic Hermite interpolation with zero slope at the two minima. Beyond
the scanned range (1.098–1.631 Å) it is continued by harmonic walls with the N1–H1A local-mode force constant of the
salt (1.99 mdyn Å⁻¹), and the grid has 801 points over that range extended by 0.6 Å on each side (sensitivity in Section
S5).

<!-- include: tables/tableS01_settings.md -->

# S2. Cartesian coordinates of the optimised structures {#sec:s-coords}

@tbl:s-energies gives, for each optimised structure, the electronic energy, the number of optimisation steps and,
where a Hessian was computed, the number of imaginary frequencies and the zero-point energy. The coordinates follow in
the order of the table, in ångström, one block per case. Structures without a Hessian are stationary on the force
criterion only; this includes the IEF-PCM sulfanilamide and sulfonamidate used as docking inputs. The same
coordinates, the relaxed periodic cells and the docked poses are supplied in machine-readable form (XYZ, CIF and SDF
files) in the ZIP file of the Supplementary material, together with the input files of the calculations; the periodic
cells and the docked poses are not printed here.

<!-- include: tables/tableS02_energies.md -->

<!-- include: tables/coordinates_S2.md -->

# S3. Complete comparison with the X-ray geometry {#sec:s-geometry}

@tbl:s-bonds lists every bond of the ion pair in the crystal and in each model, and @tbl:s-angles the deviations of
the 22 heavy-atom bond angles. The bonds to hydrogen cannot be compared on equal terms, because the X-ray C–H
hydrogens ride at 0.93 Å and the O–H hydrogens were restrained to 0.82 Å. For the heavy-atom bonds the periodic model
is closest to the crystal, with a mean absolute error of 0.017 Å against 0.024 Å for the IEF-PCM ion pair
(Section S6).

<!-- include: tables/tableS03_bonds.md -->

<!-- include: tables/tableS04_angles.md -->

# S4. Harmonic frequencies and IR intensities {#sec:s-freq}

The Hessian of the IEF-PCM ion pair comes from central differences of 162 analytic gradients and has no imaginary
frequency. @tbl:s-freq lists the 75 modes with scaled wavenumbers, IR intensities and the dominant stretch and bend
character, @tbl:s-thermo the harmonic thermochemistry of both IEF-PCM minima, and @fig:s-ir the simulated spectrum of
the salt. FT-IR and FT-Raman spectra of SAPH have been measured and compared with Hartree–Fock and B3LYP/6-311++G(d,p)
wavenumbers [@muthuselvi2017structural]; the modes listed here are those of the ion pair in the continuum and are not a
new assignment of those spectra. The neutral, O5-protonated minimum at the end of the relaxed scan has a Hessian of its
own, again without an imaginary frequency (lowest mode 17.7 cm⁻¹). The salt carries the larger zero-point energy, so the
neutral pair falls further below it once vibrations are counted: by 0.54 kcal mol⁻¹ in the electronic energy, 1.66 kcal
mol⁻¹ with the zero-point energy and 1.64 kcal mol⁻¹ in the Gibbs energy at 298.15 K. The strongest band of the neutral
pair is the O5–H1A stretch of O5–H1A···N1 at 2567 cm⁻¹ (unscaled), 377 cm⁻¹ above the N1–H1A stretch of the salt. Two
groups of bands describe the model and not the solid: the N1–H1A stretch at 2190 cm⁻¹ (unscaled), which belongs to an
N1···O5 contact much shorter than in the crystal, and the O–H stretches of the anion, whose hydrogens are free in the
ion pair but donate hydrogen bonds in the crystal.

<!-- include: tables/tableS05_freq.md -->

<!-- include: tables/tableS06_thermo.md -->

![Harmonic IR spectrum of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p) (finite-difference Hessian, 75
modes, no imaginary frequency), wavenumbers scaled by 0.9688: sticks, IR intensities (right axis); curve, Lorentzian
bands with 10 cm⁻¹ FWHM. The labelled bands carry their dominant assignment. The scale factor is the
fundamental-frequency factor of Merrick et al. for B3LYP/6-311+G(d,p) [@merrick2007evaluation], the closest basis in
their table (6-311++G(d,p) adds only diffuse s functions on H). The 2122 cm⁻¹ N1–H1A stretch (2190 cm⁻¹ unscaled)
belongs to the very short N1–H1A···O5 bond of the ion pair in the continuum (N1···O5 2.581 Å, against 2.760 Å in the
crystal), and the O–H stretches at 3679 and 3682 cm⁻¹ belong to O–H groups that are free in the ion pair but donate
hydrogen bonds in the crystal; neither band describes the solid. The computed spectrum is not compared here with the
reported FT-IR and FT-Raman spectra of the salt [@muthuselvi2017structural]. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_ir.png){#fig:s-ir}

# S5. Proton scan, dielectric ladder, functional check and zero-point energy {#sec:s-scan}

The relaxed N1–H1A scan in IEF-PCM (water) is given point by point in @tbl:s-scan, with the convergence status and the
number of steps of each point. The two chains share only the 1.35 Å point, where they differ by 0.014 kcal mol⁻¹, so
the relaxed coordinates show no hysteresis at this grid. @tbl:s-eps follows the proton along the dielectric ladder
together with the two heavy-atom markers that do not depend on its position: C1–N1 is 1.413 Å at *ε* = 4, with the
proton on O5, and 1.455 Å at *ε* = 8, with the proton on N1. @tbl:s-fcheck repeats three scan geometries with four
further methods, and @tbl:s-zpe compares the barrier with the N1–H1A zero-point energy and with the ground level of
a one-dimensional proton model.

How close the ground level of the one-dimensional model lies to the barrier top depends on how the profile is continued
and on the mass. With the walls of Section S1 the level lies 0.06 kcal mol⁻¹ below the top, with a probability of 0.84
on the O5 side. The wall on the O5 side has no physical basis, because beyond the neutral minimum the relaxed coordinate
is the soft N···O stretch: walls of 0.5, 1.0 and 5.0 mdyn Å⁻¹ there put the level 0.29 and 0.17 kcal mol⁻¹ below and
0.06 kcal mol⁻¹ above the top, with 0.91, 0.88 and 0.79 on the O5 side. Walls of 1.0 and 5.0 mdyn Å⁻¹ on the N1 side
give 0.09 and 0.03 kcal mol⁻¹ below the top (0.80 and 0.87), and effective masses of 1.008, 1.5 and 2.0 u in place of
the N–H reduced mass (0.940 u) give 0.09, 0.26 and 0.36 kcal mol⁻¹ below it (0.85, 0.89 and 0.92). The sign of the
distance from the top is therefore not determined by the model, whereas a ground state that lies mainly on the O5 side
and is not confined to the salt well is found in all eight variants.

<!-- include: tables/tableS07_scan.md -->

<!-- include: tables/tableS08_eps.md -->

<!-- include: tables/tableS09_fcheck.md -->

<!-- include: tables/tableS10_zpe.md -->

# S6. Sensitivity to dispersion and basis set {#sec:s-sensitivity}

@tbl:s-sensitivity compares the two contacts of the R²₂(8) ring and the agreement with the crystal in every model.
Without the D3(BJ) term the proton stays where it was (*δ* = −0.386 Å against −0.387 Å), but C6···O4 opens from 3.383 to
4.798 Å: in the ion pair, dispersion, not the proton, closes the second leg of the ring. The def2-TZVP basis
[@weigend2005balanced] raises *δ* by 0.095 Å, to −0.292 Å, without transferring the proton and lowers the heavy-atom
bond MAE to 0.015 Å. For the docking input the basis matters little: def2-TZVP gives sulfanilamide a *τ*~fold~ of 90.5°
(90.4° with 6-311++G(d,p)) and a C4–S1–N2 angle of 108.9°.

<!-- include: tables/tableS11_sensitivity.md -->

# S7. Periodic calculations: convergence, relaxations and the constrained proton-transfer profile {#sec:s-periodic}

The cutoff and *k*-grid tests at the X-ray geometry (@tbl:s-periodic-conv) fixed the production settings at 70/560 Ry
and a 2×1×1 grid. @tbl:s-periodic-relax lists the relaxations, and @tbl:s-periodic-hbonds the nine hydrogen bonds of
relaxed cell (i). Relaxations from the deposited structure, from a cell with H1A moved to O5 in both formula units and
from one in which only formula unit 1 had its proton moved converge to one salt, 0.04 and 0.05 kJ mol⁻¹ per formula
unit apart at PBE-D3(BJ); the last run imposed no symmetry, so each proton was free to stay on O5. A relaxation from
such a start has been tested as a way to tell salts from cocrystals: with dispersion-corrected PBE every salt of the
test set returned to the salt, but so did some cocrystals, and where both starts gave the salt the result was judged
reliable only for hydrogen bonds longer than 2.613 Å [@husak2022potential]; N1···O5 is 2.734 Å in the relaxed cell and
2.760 Å in the X-ray structure. The periodic model
also locates the proton of the short anion–anion bond, which the X-ray refinement restrained: H4A sits 1.031 Å from O4
and 1.488 Å from O5.

@tbl:s-periodic-profile follows the proton across the N1···O5 bond in the crystal. Both N1–H1A bonds of the cell were
constrained to five lengths between 1.10 and 1.70 Å, and every other atom, O5 and the H1A–O5 separation included,
relaxed. Every step costs energy: at 1.70 Å the cell lies +51.6 kJ mol⁻¹ per formula unit above the salt with
PBE-D3(BJ) and +50.0 kJ mol⁻¹ after the PBE0-D3(BJ) correction, and the constraint force points back to N1 at every
point (56–133 kJ mol⁻¹ Å⁻¹), so no neutral minimum lies on this path with either functional. Exact exchange raises the
middle of the path most (+4.1 kJ mol⁻¹ at 1.25 Å, where the proton is shared) and lowers only the last point (−1.6 kJ
mol⁻¹). Along the path N1···O5 first contracts to 2.568 Å at 1.40 Å and then opens again, while N1–H1B···O6 and
N1–H1C···O6 lengthen to 2.926 and 2.968 Å and O4–H4A···O5 to 2.566 Å. An earlier calculation held H1A and O5 in space
with the O5–H1A bond at 0.98 Å; that single configuration lies +57.1 kJ mol⁻¹ (PBE-D3(BJ)) and +52.8 kJ mol⁻¹ (with
the PBE0-D3(BJ) correction, whose cutoff and *q*-mesh variants span 0.5 kJ mol⁻¹) above the salt, 10.4 kJ mol⁻¹ higher
than the PBE-D3(BJ) profile interpolated to its N1–H1A distance: the price of the compressed O–H bond and the frozen
O5.

All relaxations used an SCF threshold of 10⁻⁸ Ry and followed a pre-relaxation at 50/400 Ry. Quantum ESPRESSO applies
distance constraints only with damped dynamics, so the profile points were relaxed that way (SCF threshold 10⁻⁹ Ry),
from cells interpolated between relaxed cell (i) and the constrained single configuration, to a residual constrained
force under 10⁻³ Ry Bohr⁻¹; the cell stayed at the X-ray values. The single configuration held both H1A and O5 because
an earlier run that fixed H1A alone let the lattice slide by 0.66 Å until N1 bonded to the proton again; that run is
not used. The PBE0 single points used 25 % exact exchange, cutoffs of 100 Ry for the wavefunctions and 400 Ry for the
density and the exchange operator, a 1×1×1 *q* mesh for the exchange and the 2×1×1 *k* grid; an 80/320 Ry set and a
2×1×1 *q* mesh were run as checks. For the single configuration the correction is almost entirely exact exchange
(−87.4 of −90.0 meV per cell).

<!-- include: tables/tableS12_periodic_conv.md -->

<!-- include: tables/tableS13_periodic_relax.md -->

<!-- include: tables/tableS14_periodic_hbonds.md -->

<!-- include: tables/tableS13b_periodic_profile.md -->

# S8. ELF, LOL and reduced density gradient {#sec:s-elf}

@fig:s-elf-lol maps the ELF and LOL of the IEF-PCM ion pair in the planes of the benzene ring, the sulfonamide group
and the R²₂(8) ring. @fig:s-rdg plots the reduced density gradient [@johnson2010revealing; @lu2025visualization] of the
ion pair in both environments; in IEF-PCM
(water) the N1–H1A···O5 spike lies at sign(*λ*~2~)*ρ* = −0.0775 a.u.

![ELF (left) and LOL (right) of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p), in three planes, each
through three atoms: (a, b) the benzene ring (C1, C3, C5), (c, d) the sulfonamide group (S1, N2, O1) and (e, f) the
R²₂(8) ring (N1, O5, O4). Common colour scale with contour lines; every atom less than 0.6 Å from a plane is shown
with its label, and N1–H1A···O5 is dashed. In (e, f) the basin of H1A is joined to N1 and separated from O5 by a
region of lower ELF and LOL, as expected for a proton that stays on N1. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_elf_lol.png){#fig:s-elf-lol}

![Reduced density gradient against sign(*λ*~2~)*ρ* for the ion pair (a) in IEF-PCM (water) and (b) in the gas phase,
B3LYP-D3(BJ)/6-311++G(d,p), from the Multiwfn grid with the density cut-off raised so that the strong hydrogen bond is
kept. Colours as in Fig. 5 of the main text (−0.08 a.u. blue, 0 green, +0.02 a.u. red); the bar above marks the
hydrogen-bond, van der Waals and steric regions. Numbered spikes are named by the atom pair nearest to their
lowest-RDG point. The N1–H1A···O5 spike of the salt lies at −0.0775 a.u.; in the gas phase the proton has moved to
O5, and the corresponding spike, O5–H1A···N1, lies at weaker sign(*λ*~2~)*ρ*. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_rdg.png){#fig:s-rdg}

# S9. Reactivity descriptors, Fukui functions, charges, electrostatic potential, QTAIM and interaction energies {#sec:s-fukui}

@tbl:s-descriptors sets the vertical ΔSCF descriptors, with hardness and electrophilicity as defined by Parr and
co-workers [@parr1983absolute; @parr1999electrophilicity], beside the orbital-energy estimates. In IEF-PCM (water) the
two stay close, but in the gas phase they part: the ionisation energy of the cation is 13.139 eV by ΔSCF and 11.185 eV
from −*E*~HOMO~. The condensed Fukui functions [@parr1984density; @yang1986use] of @tbl:s-fukui show what protonation
does to the nucleophilic site. In sulfanilamide *f*⁻ peaks on the aniline nitrogen N1 (0.293 with NPA
charges, 0.187 with Hirshfeld charges; 0.36 on the whole amino group with NPA), whereas in the cation the NH~3~⁺ group
carries almost none of it (group sum 0.04) and the weight moves to the sulfonamide group (0.30), with its largest value
on N2 (0.169).

@tbl:s-charges gives four sets of atomic charges for the IEF-PCM ion pair, NPA [@reed1985natural], Hirshfeld
[@hirshfeld1977bonded], ADCH [@lu2012atomic] and CM5 [@marenich2012charge], and they show how much of the formal ±1
survives in the pair. The cation fragment carries +0.891 e by NPA, +0.758 e by ADCH, +0.729 e by CM5 and +0.671 e by
Hirshfeld partitioning, so some charge flows back from the anion through N1–H1A···O5, by an amount that depends on the
scheme; H1A carries +0.481 e and O5 −1.189 e (NPA). ADCH gives unphysical charges on the ring C–H pairs of the IEF-PCM
pair and is used only for the fragment sum. Mulliken charges of SAPH at the Hartree–Fock and B3LYP levels were reported
earlier [@muthuselvi2017structural].

The molecular electrostatic potential was mapped on the *ρ* = 0.001 a.u. isosurface [@lu2012quantitative;
@zhang2021efficient; @murray2011electrostatic] (@fig:s-mep, @tbl:s-mep). On the IEF-PCM ion pair its extrema lie on
the atoms that make the crystal contacts the pair lacks. The deepest minimum, −74.8 kcal mol⁻¹, is next to O6, the
anion oxygen that is free in the pair but accepts both N1–H1B···O6 and N1–H1C···O6 in the crystal, and the two highest
maxima sit on exactly those donors, H1C (+97.2 kcal mol⁻¹) and H1B (+90.4 kcal mol⁻¹). The pair also has minima near
O2 (−43.5) and O1 (−42.7 kcal mol⁻¹), the acceptors of N2–H2B···O2 and O3–H3A···O1, and its next maxima, at H2B
(+64.2), H2A (+61.8), H4A (+56.1) and H3A (+51.9 kcal mol⁻¹), are the other donors of the crystal: most hydrogen bonds
that the lattice adds join a maximum of the potential to a minimum [@haritha2023fundamental]. Neutral sulfanilamide
has a much flatter potential, from −57.0 kcal mol⁻¹ at the sulfonyl O1 to +58.9 kcal mol⁻¹ at an amino hydrogen.

At the N1–H1A···O5 bond critical point (@tbl:s-qtaim) the total energy density is negative (−0.0238 a.u.), the mark of a
partly covalent contact, while that of C6–H6···O4 is positive; the integrated IGMH interaction between the ions is 0.280
a.u., and the H1A···O5 pair carries 24.4 % of the summed inter-fragment atom-pair *δG* index. The *V*/2 and
Emamian bond energies [@espinosa1998hydrogen; @emamian2019exploring] of @tbl:s-qtaim are listed only for comparison with
other work; for N1–H1A···O5 the *V*/2 value is larger in magnitude than the counterpoise-corrected interaction energy of
the entire ion pair given below, so none of them is read as the energy of a bond. @tbl:s-interaction gives
the interaction energy of the ion pair, −16.72 kcal mol⁻¹ in IEF-PCM (water), of which −4.58 kcal mol⁻¹ is the D3(BJ)
term. Interaction energies were counterpoise-corrected [@boys1970calculation], the D3(BJ) term being added separately,
and deformation energies refer to the fragments optimised in the same environment. In IEF-PCM each species has its own
cavity, which makes the correction approximate; where H1A had moved to O5 (gas phase), the pair was divided into
sulfanilamide and H~3~PO~4~. Separated, the neutral molecules lie 123.0 kcal mol⁻¹ below the gas-phase ions; each in its
own cavity, sulfanilamide and H~3~PO~4~ lie 5.82 kcal mol⁻¹ below the separated ions (electronic energies).

<!-- include: tables/tableS15_descriptors.md -->

<!-- include: tables/tableS_fukui.md -->

<!-- include: tables/tableS16_charges.md -->

<!-- include: tables/tableS17_mep.md -->

![Molecular electrostatic potential mapped on the *ρ* = 0.001 a.u. electron-density surface, B3LYP-D3(BJ)/6-311++G(d,p)
in IEF-PCM (water). (a) Ion pair in the view of Fig. 5 of the main text (front view) and (b) seen from the opposite
side (back view); (c) neutral sulfanilamide, the species used for docking. One colour scale for all panels, from −80
(red) through 0 (green) to +100 kcal mol⁻¹ (blue). Beads mark the surface minima (cyan) and maxima (amber) found by
Multiwfn; those of the ion pair with |*V*| ≥ 40 kcal mol⁻¹ and the global extrema of sulfanilamide are labelled with
their value (kcal mol⁻¹) and the nearest atom, each in the view where it faces the reader. The most negative site of
the pair is the terminal O6 of the anion (−74.8 kcal mol⁻¹); the most positive are the two N–H of the NH~3~⁺ group
that do not bind the anion, H1C (+97.2) and H1B (+90.4). In neutral sulfanilamide the range narrows to −57.0 (sulfonyl
O1) and +58.9 kcal mol⁻¹ (aniline H1C). (d) Surface area per 2.5 kcal mol⁻¹ bin of the electrostatic potential; the
surface of the ion pair is 153.4 Å² positive and 117.8 Å² negative. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_mep.png){#fig:s-mep}

<!-- include: tables/tableS18_qtaim.md -->

<!-- include: tables/tableS19_interaction.md -->

Only two of the crystal's nine hydrogen bonds are present in the ion pair. To describe all of them, five ion pairs
(135 atoms) that contain all nine bonds, eight of them at the central pair (N1–H1B···O6 joins two outer pairs), were
cut once from relaxed periodic cell (i) and once from the periodic cell with H1A held on O5 at 0.98 Å, and analysed
with B3LYP-D3(BJ)/6-31+G(d,p) [@hariharan1973influence; @francl1982self; @clark1983efficient] in the gas phase at
those geometries (@tbl:s-cluster-qtaim, @fig:s-cluster-qtaim). With the proton on O5, *ρ* at the critical points of
N1–H1B···O6 and N1–H1C···O6 falls by 34 and 29 %, and that of the anion–anion O4–H4A···O5 bond by 24 %, its total
energy density rising from −0.01293 to −0.00226 a.u.; O3–H3A···O1, whose acceptor is the sulfonyl O1 of the cation,
gains 13 %. The same cluster, taken as a finite molecule in vacuum, puts the O5-protonated arrangement 54.0 kJ mol⁻¹
per ion pair below the salt, whereas the periodic cell puts it 57.1 kJ mol⁻¹ above. The HOMO–LUMO gap of the salt
cluster is just 1.15 eV: its outer ions lack the neighbours, and part of the hydrogen bonds, that the lattice gives
them, and all five pairs were switched together. The sign of the cluster energy therefore shows only that a finite
piece of the network in vacuum does not keep the proton on N1; we use the cluster for its bond descriptors and do not
interpret the size of its energy.

<!-- include: tables/tableS18b_cluster_qtaim.md -->

![QTAIM of all nine crystal hydrogen bonds in a five-pair cluster (135 atoms) cut from the relaxed periodic salt
(circles) and from the periodic cell with H1A held on O5 at 0.98 Å in every pair (diamonds), B3LYP-D3(BJ)/6-31+G(d,p)
single points in the gas phase at the periodic PBE-D3(BJ) geometries; values at the contact of the central ion pair,
except N1–H1B···O6, which joins two outer pairs. (a) Electron density *ρ* at the bond critical point, with its change
in per cent (a: with H1A on O5 the critical point is that of N1···H1A). (b) Total energy density *H* at the same
points; *H* < 0 marks a partly covalent contact. Moving the proton to O5 weakens the other two N–H···O bonds donated
by NH~3~⁺ (−34 and −29 %) and the anion–anion O4–H4A···O5 bond (−24 %), whose *H* rises from −0.01293 to −0.00226
a.u.; O3–H3A···O1, whose acceptor is the sulfonyl O1 of the cation, gains 13
%. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_cluster_qtaim.png){#fig:s-cluster-qtaim}

# S10. Hirshfeld surface decomposition {#sec:s-hirshfeld}

@tbl:s-hirshfeld splits the Hirshfeld surfaces of the two ions into reciprocal element-pair contacts, and
@fig:s-hirshfeld-decomp shows the matching fingerprint plots. H···O/O···H contacts cover 51.1 % of the cation surface
and 80.3 % of the anion surface. @tbl:s-hirshfeld-contacts lists every red spot with its contact distance; the deepest,
*d*~norm~ = −0.842, is the H4A···O5 contact of 1.595 Å between neighbouring anions.

Hirshfeld surfaces [@spackman2009hirshfeld] of the two ions were computed with our own code (checked against Multiwfn)
from the X-ray structure (the atom-site coordinates of the deposited CIF), with X–H bonds set to neutron-normalised
lengths (C–H 1.083, N–H 1.009, O–H 0.983 Å [@allen2010bond]), every ion with an atom within 6.0 Å of the central one,
and spherical free atoms (ROHF/aug-cc-pVQZ) as the promolecule. *d*~norm~ follows the CrystalExplorer definition
[@spackman2021crystalexplorer] with Bondi radii (H 1.20 Å); the *d*~norm~ of Multiwfn uses another reference radius,
so only ours is reported. Fingerprint plots [@spackman2002fingerprinting] were split by element pair with reciprocal
contacts combined.

<!-- include: tables/tableS20_hirshfeld.md -->

<!-- include: tables/tableS21_hirshfeld_contacts.md -->

![Complete decomposition of the Hirshfeld surfaces of Fig. 6 of the main text (X-ray structure, neutron-normalised
X–H). (a) Share of every reciprocal contact type in the surface of the cation (blue) and of H~2~PO~4~⁻ (orange); shares
below 0.1 % are not labelled. (b) Cation and (c) anion: fingerprint plots of every contact type with at least 0.5 % of
the surface, each over the full plot of the same ion (grey); bins of 0.02 Å, colours as in Fig. 6 of the main text.
Contacts of the S and P atoms cover at most 0.2 % of either
surface. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_hirshfeld_decomp.png){#fig:s-hirshfeld-decomp}

# S11. TD-DFT excitations {#sec:s-tddft}

@tbl:s-tddft-pair lists the first ten singlet states of the ion pair with CAM-B3LYP and with B3LYP, and
@tbl:s-tddft-fragments those of the cation and of neutral sulfanilamide. The CT number tells the two functionals apart.
Under B3LYP, states S2 to S5 of the ion pair shift electron density from anion to cation (CT number 0.61–0.91); under
CAM-B3LYP none of the ten does so to a large extent, and the most intense low-lying band of the pair (S2, 211 nm) sits
close to S2 of the isolated cation (208 nm). The lower cut-off wavelength reported for SAPH, 223 nm
[@muthuselvi2017structural], lies between S1 of the ion pair (232 nm, *f* = 0.015) and its first strong band (211 nm); a
cut-off is an absorption edge that depends on the sample, so the comparison is not taken further. Both tables are also
supplied as CSV files in the ZIP file of the Supplementary material.

The states are linear-response TD-DFT excitations [@casida1995time]; non-equilibrium solvation used an optical
dielectric constant of 1.78. In a Tamm–Dancoff test on the gas-phase pair at
the X-ray geometry, B3LYP placed the first excitation at 2.44 eV, in the visible, for a crystal that is colourless, and
the first B3LYP state of appreciable intensity of the ion pair in IEF-PCM (S2, 227 nm, *f* = 0.14) is a charge-transfer
state whose hole and electron lie 3.39 Å apart.

<!-- include: tables/tableS22_tddft_pair.md -->

<!-- include: tables/tableS23_tddft_fragments.md -->

# S12. Finite-field NLO details {#sec:s-nlo}

*Pc* has no inversion centre, so a SAPH crystal may generate a second harmonic, which a centrosymmetric crystal cannot
[@kurtz1968powder], and the salt has been described as a nonlinear optical (NLO) material [@muthuselvi2017structural].
@tbl:s-nlo gives the dipole moment, the isotropic and anisotropic polarisability and the first hyperpolarisability in
the Taylor and Buckingham conventions (1 a.u. = 8.639 × 10⁻³³ esu; the B convention halves every *β*), with *β*~tot~ the
norm of $\beta_i = \tfrac{1}{3}\sum_j(\beta_{ijj} + \beta_{jij} + \beta_{jji})$, for the ion pair in IEF-PCM (water) and
in the gas phase at the same geometry, for sulfanilamide and for urea. The ion pair in water has *μ* = 7.25 D, *ᾱ* =
202.7 a.u. and *β*~tot~ = 406.6 a.u. (3.512 × 10⁻³⁰ esu). Referenced to gas-phase urea, the gas-phase ion pair gives the
like-for-like ratio, 6.4 (573.2 a.u.); the value of 4.5 for the pair in IEF-PCM (water) mixes two environments. Neutral sulfanilamide
in water reaches 1807.7 a.u., several times the value of the pair, the drop expected once the amino donor of this
push–pull system is protonated. None of these numbers is a property of the crystal: a molecular *β* at one rigid
geometry leaves out vibrational and frequency-dependent terms, local fields and the orientation of the ions relative to
the polar axis, and structural fluctuations can spread such responses widely around their rigid-geometry value
[@castet2022predicting].

The finite fields [@kurtz1990calculation] were applied as SCF dipoles at ±0.001 and ±0.002 a.u. along each axis and each
pair of axes, and the derivatives taken by central differences with one Richardson step. The values are static,
electronic and molecular, without local-field factors or packing. Urea, in the gas phase, is the reference. Whether the
hydrogen-bond network of SAPH aligns the cation dipoles as efficiently as in crystals engineered for frequency doubling
[@zhao2022engineering] could be settled only by a powder SHG measurement or a calculation of the crystal susceptibility,
and neither is reported here. For two other polar salts of the same cation, one of them the sulfamate, calculated
susceptibilities gave a second-harmonic efficiency comparable with that of potassium dihydrogen phosphate
[@wojnarska2018crystal].

<!-- include: tables/tableS24_nlo.md -->

# S13. Docking: preparation, search settings, validation, sensitivity, torsional strain and a control receptor {#sec:s-docking}

Each search box is a 20 Å cube centred on the heavy-atom centroid of the crystal ligand in the prepared receptor. The
AutoDock4Zn maps of 6RL9 and 3HS4 put the zinc pseudo-atom 0.17 and 0.07 Å from the sulfonamide nitrogen of the crystal
ligand. Every job used exhaustiveness 32 and 10 seeds, 20 for the two production AutoDock4Zn jobs of sulfanilamide.
Scores repeat exactly for a repeated seed. The cross-docking references are the crystal ligands of 6RL9 and 3HS4
carried into the other structure by superposing 257 Cα atoms (root-mean-square deviation, RMSD, 0.328 Å), and every site
check passes on them.

@tbl:s-dock-validation collects all validation and sensitivity runs behind Table 6 of the main text, each changing one
setting. A larger exhaustiveness did not remove the bimodality of the 6RL9 AutoDock4Zn redock (2.34 Å at 64 and 2.33 Å
at 128, against 2.07 Å at 32), and 20 seeds instead of 10 raised its success fraction only from 0.5 to 0.6. The His64
tautomer, a second loop rebuild of 3TZF chain B, the sulfate of 1AJ0 and a start from the crystal conformer left the
RMSDs of their runs within a few hundredths of an ångström. In 1AJ0 only a box made smaller after the failure brought
sulfanilamide close to the pocket, and that box was not used. The input-geometry test of Section 3.6 of the main text
was extended by six runs, counted among the eight of the note of its Table 6: into both hCA II structures were docked
the DFT geometry of sulfanilamide with its C4–S1–N2 angle set to the ETKDG value, the ETKDG conformer with the DFT
angle, and the ETKDG conformer unchanged. The angle alone switched the top-ranked ring
orientation in 6RL9 but not in 3HS4, while the near-native pose appeared in the ranked list of every seed and the zinc
and Thr199 contacts held throughout: the sulfonamide nitrogen is the atom on the zinc in all eight AutoDock4Zn runs of
sulfanilamide, these six and the two production runs (run means 2.09–2.37 Å), although the 2.3 Å criterion of the site
check is met by only 0.30 and 0.00 of the seeds in the two 3HS4 runs from the DFT geometry. These six runs are kept out
of the production tables and the strain analysis. An earlier comparison of AutoDock4Zn with MM/GBSA for carbonic
anhydrase inhibitors also found docked poses that deviate from the crystal structures [@taylor2023mm]. Docked poses of
benzenesulfonamides in hCA II, bound to the zinc through the deprotonated nitrogen, have been checked by cross-docking
before the complexes were treated with a quantum-mechanical/molecular-mechanical method [@samanta2016prediction].

Three settings complete the description of the protocol. After the hydrogen atoms were added, the protein of every
receptor was energy-minimised in OpenMM with the Amber14 force field, in vacuum and before the metal ions, cofactors and
kept waters were put back, the deposited heavy atoms restrained with 10⁵ kJ mol⁻¹ nm⁻² and the hydrogen and rebuilt
atoms free; no deposited atom moved by more than 0.20 Å in the five receptors of Section 2.5 of the main text and 0.22 Å
in the control receptor. The AutoDock4Zn maps have a spacing of 0.375 Å and 54 intervals (55 points) per axis, a 20.25 Å
grid on the centre of the 20 Å Vina box. Every search wrote up to 20 modes within 3 kcal mol⁻¹ of the best one.

Production (@tbl:s-dock-scores, @fig:s-dock-seeds) docked every species into each validated receptor with Vina and
Vinardo, and the neutral sulfonamides into hCA II with AutoDock4Zn as well, 72 jobs in all. The spread between seeds is
small: the largest SD of a mean score is 0.20 kcal mol⁻¹, so the orders discussed in the main text do not come from
sampling noise. LE, the score per heavy atom with its sign reversed, and the fit quality FQ, which scales LE by the size
dependence found for measured affinities [@reynolds2007role], are listed for completeness only; applied to docking
scores, the scale is a heuristic. PoseBusters passes a fraction of 0.993 of the production top poses. The only check
that fails is the minimum distance to the protein, for half of the acetazolamide-anion poses scored by Vinardo in 6RL9,
whose nearest receptor atom is the zinc.

The cost of the docked torsions comes from relaxed scans of C3–C4–S1–N2 in IEF-PCM (water), with the torsion held and
everything else optimised (@tbl:s-torsion-scan). For neutral sulfanilamide the strain grows steadily from the
perpendicular minimum (90.4°), reaching 1.89 kcal mol⁻¹ when S1–N2 lies in the plane of the ring, and the torsions of
the crystal ligands, 105.6° (3TZF), 111.3° (1AJ0) and 132.1° (6RL9), cost 0.18, 0.28 and 0.84 kcal mol⁻¹. The
sulfonamidate is flat within 0.10 kcal mol⁻¹ over the whole range. @tbl:s-dock-torsion sets the docked torsions of
every production job beside the crystal ligand of its receptor and evaluates the scan at the torsion of each seed:
0.23–0.28 kcal mol⁻¹ for sulfanilamide in 3TZF with Vina, 1.09–1.43 kcal mol⁻¹ in 6RL9 with AutoDock4Zn, where the
poses turn 15° further than the crystal ligand, and at most 1.81 kcal mol⁻¹ in any job. The MMFF94s estimate, the
molecule minimised with the torsion restrained and referred to the global MMFF94s minimum, follows the DFT values for
sulfanilamide within 0.31 kcal mol⁻¹ but overstates the sulfonamidate by 0.20–1.02 kcal mol⁻¹; the cation was not
scanned. @tbl:s-dock-contacts gives the polar contacts of the cluster representatives. It includes the Vina poses in
hCA II, which show the missing metal term: the atom nearest to the zinc is a sulfonyl oxygen (2.28 Å) in 6RL9 and the
sulfonamide nitrogen at 2.50 Å in 3HS4.

Receptor preparation kept one chain and altloc, with the metals and cofactors of the site; PDBFixer rebuilt missing side
chains and one internal gap (Ala66–Ala67 in 3TZF chain B). In hCA II, His94 and His96 were set to HID and His119 to HIE
from their coordination; His64, the proton shuttle, touches neither ligand nor metal and was assumed HIE; HID, tested
with Vina only, moved the scores by at most 0.03 kcal mol⁻¹. PROPKA was run on the protein without the zinc. This gave
Glu106 of both hCA II structures a pK~a~ of 8.91, above the pH of 7.4, so the residue was built as the neutral acid;
with the zinc in the input PROPKA gives 5.33 for 6RL9, a carboxylate. Every hCA II docking number of the main text and
of this section comes from the receptors with the acid unless the carboxylate is named, and Section S15 repeats eight
AutoDock4Zn and 21 Vina and Vinardo runs with the carboxylate. Waters were kept only within 2.6 Å of Mg²⁺: three in 3TZF and
three in 7S2K. hCA II was therefore docked without any water. In 6RL9 and 3HS4 the sulfonamide nitrogen of the crystal
ligand is the fourth ligand of the zinc, so no zinc-bound water or hydroxide had to make way for the docked ligand; what
the rule removed from 6RL9 includes four fully occupied waters of the active site (423, 516, 553 and 568 of the entry),
on which the ring-rotated docked pose lies. Section S15 tests what keeping them does. All ligands were defined by SMILES
with explicit charges checked against the RCSB chemical components. Beyond the RMSD criteria, the validation gate
required a PoseBusters physics pass fraction of at least 0.8 (RMSD test excluded) and, for hCA II, Zn–N(sulfonamide) ≤
2.3 Å and Thr199 OG1···N ≤ 3.3 Å in at least half of the seeds; every site check had to pass on the crystal ligand
itself. Top poses were clustered by a leader algorithm, best score first.

Vina, without a metal term, finds the pocket but not the anchor: its 6RL9 poses lie 1.10 ± 0.01 Å from the crystal
ligand, yet the sulfonamide nitrogen stays 2.84 Å from the zinc. The AutoDock4Zn ranks do not check the site either:
with AutoDock4Zn sulfamethoxazole ranks first in both hCA II structures, although its top poses stay away from the zinc
(mean Zn–N 3.56 and 3.28 Å). In the cluster representatives the ring of sulfanilamide lies against Trp209, Leu198,
Phe131 and Val121 in 6RL9, and in 3HS4 its aniline nitrogen meets Asn67 OD1 (2.57 Å). Because the programs keep the
input bond angles, the DFT geometry enters the docking directly: with its C4–S1–N2 angle of 109.5°, close to the 109.2°
of the SAPH crystal, the 6RL9 crystal pose is the top pose of all 20 seeds in the receptor as prepared, with Glu106 as
the neutral acid. The input geometry is not the only thing that selects the orientation: with Glu106 as a carboxylate
the same geometry gives the crystal pose as the top pose of 15 of the 20 seeds, although a near-native mode stays in
the ranked list of every seed (Section S15).

A docked mode counts as near-native when its heavy-atom RMSD to the crystal ligand is at most 1.0 Å. In the bimodal
6RL9 validation run every seed holds both orientations, and the best near-native and the best ring-rotated mode
(beyond 2.0 Å) differ by 0.008–0.039 kcal mol⁻¹; over the eight AutoDock4Zn runs of sulfanilamide the best-ranked
near-native mode trails the top mode by 0.000–0.219 kcal mol⁻¹. The 2.0 Å limit takes in mostly, but not only, the
ring-rotated pose: counted as the modes within 2.0 Å of one reference mode of that pose (Section S15), it accounts for
20 of the 25 modes beyond 2.0 Å in the validation run and for 40 of the 62 in the 20-seed run from the DFT geometry.
The acetazolamide check of the validation rests on the
binding constants compiled in the Protein–Ligand Binding Database (PLBD) [@linge2023plbd] from thermal shift assays of
human carbonic anhydrase II [@linkuviene2018thermodynamic]. At 37 °C and pH 7.0 the observed dissociation constants
are 46 nM for acetazolamide (three experiments) and 11 µM for sulfanilamide (two experiments). PLBD also gives
intrinsic constants, from which the protonation reactions of ligand and protein linked to binding are subtracted, and
these are equal: 9.8 nM for both, an intrinsic Gibbs energy of binding of −47.55 kJ mol⁻¹ for acetazolamide and −47.56
kJ mol⁻¹ for sulfanilamide. The measured order that the docking scores follow is therefore that of the observed
constants, set by how much of each sulfonamide is ionised at pH 7.0, which a score of the neutral ligand does not
contain.

A control receptor shows what the scoring returns for a pocket that is not a target of sulfanilamide
(@tbl:s-control-receptors, @tbl:s-control-scores, @fig:s-control). It is the ATP site of DNA gyrase B of *Escherichia
coli*, in the 24 kDa N-terminal fragment crystallised with a 4,5′-bithiazole inhibitor of low micromolar activity (4DUH,
chain A, 1.50 Å) [@brvar2012structure]. This is the site of ATP and of the ATPase inhibitors of gyrase, of which only
novobiocin reached clinical use [@bisacchi2015class]. The antibacterial target of the sulfonamides is DHPS (Section 1 of
the main text), and the fluoroquinolones, the class of ciprofloxacin, bind gyrase in its DNA cleavage core together with
the DNA [@blower2016crystal; @bush2020quinolones]; neither sulfanilamide nor ciprofloxacin is a known ligand of this
site. No measurement shows that sulfanilamide does not bind there, so the control tests the scores and not the binding.
The 4-bromobenzylidene imine of sulfanilamide has been docked into the ATPase domain of the *Staphylococcus aureus*
enzyme (PDB 3U2D) [@muthukumar2022synthesis]. The 4DUH receptor was prepared like the five receptors of Section 2.5 of
the main text: one chain, no water (the protocol keeps water only on Mg²⁺, and this site holds no metal), PROPKA
protonation at pH 7.4, a 20 Å box on the centroid of the crystal inhibitor, exhaustiveness 32 and 10 seeds. Vina docked
the eight species of the production runs, ciprofloxacin as its zwitterion and the crystal inhibitor as its benzoate;
Vinardo docked sulfanilamide, ciprofloxacin and the inhibitor. Chain A lacks residues Tyr109–His116 of the ATP lid.
PDBFixer built them, no closer than 8.3 Å to the crystal inhibitor, and because a loop built without density is a guess,
a second receptor was prepared from another seed (the two builds differ by 2.2 Å RMS in the built residues) and a third
is the first without the eight built residues. A fourth is the first with 322 ordered waters of chain A added as fixed
oxygen atoms, 37 of them inside the box; it tests the water rule of the protocol and is not part of it.

Under the protocol, redocking does not return the crystal inhibitor to its crystal pose. With Vina the top pose of every
seed lies 6.68 ± 0.01 Å from it, with Vinardo 6.77 ± 0.02 Å, so the receptor fails the 1.0 Å gate and the 2.0 Å
criterion, as 1AJ0 did. The search is not at fault: with Vina every seed holds a mode 0.45–0.53 Å from the crystal pose,
ranked second or third, 0.5–0.8 kcal mol⁻¹ behind the top pose (Vinardo: ranked second, 0.1–0.4 kcal mol⁻¹ behind). Both
poses keep the amide N–H of the inhibitor on the carboxylate of Asp73 (2.90 Å in the top poses, 2.82 Å in the crystal).
They differ in the benzoate, which reaches the guanidinium group of Arg136 in the crystal (3.06 Å) and turns under the
ATP lid in the top pose, 15.6 Å from Arg136 (panel a of @fig:s-control). Two things can favour the turned pose. The
functions have no electrostatic term [@trott2010autodock; @quiroga2016vinardo], so the salt bridge to Arg136 counts only
as a hydrogen bond and the carboxylate beside the carbonyl oxygens of Lys103 and Gly117 (3.2–3.4 Å) costs nothing; and
the protocol empties the space under the lid, where the turned benzoate lies on three ordered waters of the crystal and
the crystal pose on none. The second was tested, with Vina. With the ordered water kept as fixed oxygen atoms, the top
pose of all ten seeds is the crystal pose (0.52 ± 0.02 Å), and this receptor passes the gate. The test is a lenient one:
all the ordered water went in at once, and it is the water of the inhibitor complex itself. The built residues are not
the cause: the second build gives 6.67 ± 0.02 Å and the receptor without them 6.68 ± 0.02 Å, with the crystal pose again
among the modes of every seed.

Clorobiocin, a coumarin antibiotic crystallised with the same fragment (1KZN) [@lafitte2002dna], gives no cross-docking
test. Carried into the frame of 4DUH by superposing 183 sequence-aligned Cα atoms (RMSD 0.78 Å), its crystal pose comes
within 1.65 Å of Lys103 NZ and has 11 heavy-atom contacts below 2.5 Å with the receptor, eight of them with Gly101 and
Lys103, lid residues that are ordered in 4DUH and not modelled in 1KZN. 4DUH cannot hold that pose, so clorobiocin was
not docked.

The receptor as prepared fails the gate, so none of its poses is read as a binding mode, and the receptor that passes it
holds water that the protocol removes. The control is therefore used only for the level of the scores, which are given
for both (@tbl:s-control-scores, panel c of @fig:s-control). With Vina, neutral sulfanilamide scores −5.96 ± 0.01 kcal
mol⁻¹ in the receptor as prepared and −5.93 ± 0.01 kcal mol⁻¹ with the crystal water, against −6.42 in 3TZF, −6.41 in
7S2K, −5.88 in 6RL9 and −5.77 kcal mol⁻¹ in 3HS4. The control receptor thus scores it 0.4–0.5 kcal mol⁻¹ short of DHPS
and of Sul2, the enzyme that rejects sulfonamides, and level with the Vina scores in hCA II, inside the range of the
four production structures. The second build gives −5.85 and the receptor without the built residues −5.83 kcal mol⁻¹,
so every variant of the control receptor lies within 0.6 kcal mol⁻¹ of every production structure. This is less than the
0.7 kcal mol⁻¹ by which acetazolamide scores better here than in DHPS. Vinardo gives −4.97 kcal mol⁻¹ here, against
−4.99 to −5.31 kcal mol⁻¹ in the four structures, and sulfanilamide in the surface cleft of 1AJ0, 11.64 Å from its
crystal pose, scored −5.80 kcal mol⁻¹ with Vina (@tbl:s-dock-validation).

The poses of sulfanilamide differ between the two receptors (panel b). As prepared, all ten seeds place it under the
lid, 5.0 Å from the nearest atom of the crystal inhibitor (pose centroid), 3.5–3.6 Å from the built residues and on four
crystal waters. This site depends on the built residues: the second build returns it in seven seeds and the cavity of
the inhibitor in three, and without the built residues all ten seeds lie in the cavity. With the crystal water, the four
waters occupy the first site and all ten seeds place sulfanilamide in the cavity, its aniline nitrogen 3.05 Å from Asp73
OD1, where the amide of the inhibitor sits. That pose has the look of a binding mode, in a receptor that is not known to
bind sulfanilamide. A search over the whole fragment in the receptor as prepared, outside the protocol (box 53.5 × 56.0
× 57.0 Å, exhaustiveness 683, five seeds), scores sulfanilamide −5.96 ± 0.01 kcal mol⁻¹ and returns the pose of the 20 Å
box in all of them.

Over the eight species of the production runs the Vina score in 4DUH as prepared differs from that in 3TZF by +0.16 kcal
mol⁻¹ on average (−0.74 to +0.81), from that in 7S2K by +0.74 (+0.43 to +1.96), from that in 6RL9 by +0.09 and from that
in 3HS4 by −0.02 kcal mol⁻¹, a positive difference being a weaker score in the control receptor; with the crystal water
the four averages are +0.25, +0.83, +0.18 and +0.08 kcal mol⁻¹. *p*-Aminobenzoate, the substrate of DHPS and Sul2, loses
most (0.81 and 1.96 kcal mol⁻¹). Acetazolamide, which hCA II binds more tightly than sulfanilamide, scores −6.27 kcal
mol⁻¹ in the control receptor against −6.65 and −6.41 kcal mol⁻¹ in 6RL9 and 3HS4. Ligand size accounts for the control
scores: over these eight species the score follows the number of heavy atoms with *r* = −0.99 in 4DUH (−0.99 with the
crystal water), as in 6RL9 and 3HS4 (both −0.99), where Vina has no metal term and fails the gate, and less closely in
3TZF (−0.68) and 7S2K (−0.76). In a recent benchmark several docking programs predicted affinity no better than a
regression on molecular weight [@tornesakis2026evaluating].

Ciprofloxacin (24 heavy atoms) scores −8.29 ± 0.64 kcal mol⁻¹ as prepared, in three clusters, and −7.77 ± 0.26 kcal
mol⁻¹ with the crystal water, against −9.27 and −9.04 kcal mol⁻¹ for the inhibitor (26 heavy atoms). The inhibitor of
the site thus scores better than a drug that binds gyrase elsewhere, by 1.0 kcal mol⁻¹ as prepared and 1.3 kcal mol⁻¹
with the crystal water, but this is weak evidence. About 0.4 kcal mol⁻¹ of the gap is what two more heavy atoms give on
the size line of the control receptor (−0.21 kcal mol⁻¹ per heavy atom over its ten species, *r* = −0.99); as prepared,
the gap lies within the spread of ciprofloxacin over its seeds (−7.6 to −9.2 kcal mol⁻¹) and belongs to a top pose of
the inhibitor that is not the crystal pose; and it is one ordered pair. The control does not show that sulfanilamide
fails to bind gyrase B. It shows that a Vina score near −6 kcal mol⁻¹ and a pose with polar contacts are what the
protocol returns for sulfanilamide in a pocket that is not its target, so that neither identifies a target.

<!-- include: tables/tableS25_dock_validation.md -->

<!-- include: tables/tableS26_dock_scores.md -->

![Seed-to-seed spread of the docking results. (a–j) Top-pose score of every seed (dots; 10 seeds, 20 for sulfanilamide
with AutoDock4Zn) and the mean (bar) for every species docked to each receptor with each scoring function. The three
sulfanilamide species (green) were docked from their B3LYP-D3(BJ)/IEF-PCM geometries, the controls (grey) from
ETKDG/MMFF94s conformers. SMX, sulfamethoxazole; AZM, acetazolamide; pABA, *p*-aminobenzoate; DHP,
7,8-dihydropteroate, the crystal ligand of 7S2K; superscripts give the charge. Scores are empirical docking scores
in kcal mol⁻¹, not free energies of binding. (k) Top-pose RMSD of every seed from the crystal reference (heavy atoms,
symmetry-corrected, no superposition) for the redocking and cross-docking jobs of Fig. 8 of the main text, including
the bimodal validation run in 6RL9 (ETKDG input) and the two AutoDock4Zn runs from the DFT geometry. Solid line, the
1.0 Å gate; dashed line, the conventional 2.0 Å criterion. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_dock_seeds.png){#fig:s-dock-seeds}

<!-- include: tables/tableS26b_torsion_scan.md -->

<!-- include: tables/tableS27_dock_torsion.md -->

<!-- include: tables/tableS28_dock_contacts.md -->

<!-- Caption copied from figures/captions.md (group B, figS_control_gyrb, 4 Oct 2026); keep the two in step. -->

![Control receptor: the ATP site of *Escherichia coli* DNA gyrase B (4DUH, chain A), which is not a target of
sulfanilamide; AutoDock Vina, 20 Å box. (a) Redocking of the crystal inhibitor (grey carbons). In the receptor as
prepared the top pose (orange; seed of the median RMSD, 6.68 Å) turns its benzoate under the ATP lid, onto three ordered
waters of the crystal (red spheres); in the receptor with the crystal water kept the top pose (purple, 0.52 Å) is the
crystal pose. Both keep the amide on Asp73. (b) Representative poses of neutral sulfanilamide (B3LYP-D3(BJ)/IEF-PCM
geometry), seen across the line that joins them: as prepared (orange) it lies under the lid, beside the built residues
and on four crystal waters (red spheres); with the crystal water kept (purple) it lies in the cavity of the inhibitor,
its aniline nitrogen 3.05 Å from Asp73 OD1 (dashed line). The panel title gives the mean scores over the seeds (−5.96
and −5.93 kcal mol⁻¹), and the crystal inhibitor is drawn thin. Thin sticks: pocket residues within 3.8 Å of a docked
pose, with the carbons of the residues that PDBFixer built (Tyr109–His116) in pink. (c) Mean top-pose Vina score against
the number of heavy atoms of the ligand: the species docked into 4DUH as prepared (circles) and with the crystal water
kept (diamonds), bars the SD over the seeds, and the same species in the four production structures (open symbols, set
side by side at each ligand size). Dashed line: linear fit to the scores in 4DUH as prepared (*r* = −0.99). Scores are
empirical docking scores, not free energies of binding. SMX, sulfamethoxazole; AZM, acetazolamide; pABA,
*p*-aminobenzoate; superscripts give the charge. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_control_gyrb.png){#fig:s-control}

<!-- include: tables/tableS29_control_receptors.md -->

<!-- include: tables/tableS30_control_scores.md -->

# S14. ADMET and toxicity sanity check {#sec:s-admet}

The descriptors and ADMET-AI predictions for sulfanilamide, with the cation for reference, are compared with
experimental data in @tbl:s-admet.

<!-- include: sections/admet.md -->

<!-- include: tables/tableS_admet.md -->

# S15. Docking into carbonic anhydrase II: the protonation of Glu106 and four ordered waters {#sec:s-dock-tests}

Two properties of the hCA II receptors of Section S13 were tested after the docking was complete: the protonation of
Glu106 and the absence of water. @tbl:s-dock-glu106-water gives the runs of sulfanilamide in 6RL9 for both states of
Glu106, without and with four crystal waters, and @tbl:s-dock-glu106-other the other AutoDock4Zn runs that were
repeated. The docking numbers of the main text and of Section S13 are left as computed, for the receptors with the
acid: we know of no measurement of the charge of Glu106 in this complex, and
the test changes this one residue in an otherwise identical receptor, so that every difference between the two states
in the tables can be assigned to it. The whole docking set was not run again with the carboxylate.

PROPKA [@olsson2011propka3] had been run without the zinc in its input. It then gives Glu106 a pK~a~ of 8.91 in both
structures, above the pH of 7.4, so the residue was built as the neutral acid; with the zinc it gives 5.33 for 6RL9
(3HS4 was not rerun), a carboxylate. The acid is an effect of the preparation and not a choice. Structural work points
to the carboxylate: in the unliganded enzyme Thr199 donates its hydrogen to Glu106, taken there as deprotonated
[@eriksson1988refined], and neutron diffraction of the acetazolamide complex shows the sulfonamide N–H donating to
Thr199, which donates to Glu106 [@fisher2012neutron]; neither measures the charge of the residue. In the docked
receptors the extra proton sat 4.84 and 4.68 Å from the nearest atom of the crystal ligand.

For the test the proton was deleted and five atoms were recharged, with no heavy atom moved. Eight AutoDock4Zn runs
(those without water in the two tables) and 21 Vina and Vinardo runs were repeated with identical ligand inputs. Vina
and Vinardo do not see the proton: 80 stored top poses scored without a search change by 0.000 kcal mol⁻¹ and no
verdict changes, so the His64 test of Section S13, made with Vina, could not have seen a receptor proton either.
AutoDock4Zn does. Scored without a search in the maps of both receptors, the stored modes of the four sulfanilamide
runs shift with the carboxylate so that the ring-rotated orientation gains 0.013 to 0.018 kcal mol⁻¹ on the
near-native one, the size of the separation of the two orientations, and the ranks follow. From the DFT geometry the
top pose in 6RL9 lies within 1.0 Å of the crystal ligand in 15 of 20 seeds, where the acid gave 20 of 20 (1.34 ± 1.30
Å against 0.61 ± 0.04 Å), and the sulfonamidate rises from 3 to 7 of 10 seeds and now meets the strict gate.

Other results do not move: in the four repeated sulfanilamide runs the sulfonamide nitrogen stays on the zinc (run
means of Zn–N 2.13 to 2.31 Å), a near-native mode is present in every seed at rank 1 to 3, and the two orientations
stay within 0.058 kcal mol⁻¹. What was not repeated is listed under @tbl:s-dock-glu106-other. The test does not
establish the protonation state of Glu106. It shows that the orientation which AutoDock4Zn ranks first depends on that
state as well as on the input geometry.

The second test came from the preparation of the simulations (Section S16). The ring-rotated pose lies on four water
sites of the 6RL9 crystal: the oxygens of waters 423, 516, 553 and 568 are 0.19 to 1.43 Å from one of its heavy atoms.
All four are fully occupied, with B factors of 11.5 to 17.6 Å² against a median of 22.5 Å² for the water of the chain,
and coexist in the crystal with the ligand in its crystal pose. The protocol keeps no water in hCA II, so the docking
could not see them, although docking of other complexes has been found more accurate with crystal water kept
[@roberts2008ligandprotein]. With the four waters put back as fixed atoms 16 runs were repeated, eight of them with
AutoDock4Zn (its rows with waters in @tbl:s-dock-glu106-water), the 20-seed run under three placements of the water
hydrogens.

None of these eight runs writes a mode of the rotated pose, against 40 and 44 such modes in the 20-seed run without
the waters (acid and carboxylate), and the top pose lies within 1.0 Å of the crystal ligand in every seed: 0.51 ± 0.04
Å with the acid and 0.51 ± 0.03 Å with the carboxylate in the 20-seed run, so the difference between the two states
disappears. The validation runs of sulfanilamide, which failed both criteria without the waters, pass them. Stored
modes of the rotated pose scored in the maps with the waters return +5,346 to +7,644 kcal mol⁻¹, a clash and not an
energy to interpret. Vina and Vinardo fail without and with the waters, with run means of Zn–N of 2.40 to 2.93 Å over
their 16 runs, of which @tbl:s-dock-glu106-water lists the eight from the ETKDG conformer.

The rows with the waters do not validate the protocol. A rigid receptor cannot let a ligand displace water, so the
rotated pose is excluded by construction, and a receptor that holds the water of the complex's own crystal has been
told part of the answer, as the control receptor of Section S13 was. The test is a diagnostic of the degeneracy: the
pose that competes with the crystal pose in the water-free receptor needs the space of these waters. It does not show
that the rotated orientation cannot occur. The four were always kept together, no cost of displacing them is obtained,
and ligands do displace water: acetazolamide displaces four waters of the unbound active site [@fisher2012neutron], and
a glycerol molecule occupies these very sites in 3HS4, for which no receptor with water was made. Nor are the four
places always filled: in the production of the nine simulations of Section S16, in which no frame is rotated-like,
those of waters 423 and 516 are occupied in only 12.4 to 32.8 % and 10.6 to 79.5 % of the all-atom frames, while those
of waters 553 and 568 are occupied in 86.8 to 96.9 % and 71.5 to 88.2 %; the water of the site is not converged in
these runs.

<!-- include: tables/tableS31_dock_glu106_water.md -->

<!-- include: tables/tableS32_dock_glu106_other.md -->

# S16. Molecular dynamics of the carbonic anhydrase II–sulfanilamide complex {#sec:s-md}

The simulations ask whether the complex keeps the orientation of the ring from which it starts: the crystal pose, or
the ring-rotated pose that AutoDock4Zn ranks level with it (Sections S13 and S15). @tbl:s-md-settings lists the
settings, @tbl:s-md-runs the eleven runs and @tbl:s-md-gates the tests of the force field, with those that failed;
[@fig:s-md-orientation; @fig:s-md-rmsd] show every run as a time series. The complex has been simulated before, from
three docked conformations and with no bonded zinc terms described [@nguyen2025understanding].

*Rules fixed beforehand.* The rules for reading the runs were written down, with criteria and dates, before the
production (file `md/PLAN.md` of the data package). Two points were settled after a rehearsal of the pipeline, three
runs of 5 ns that enter no count, made with the zinc set and the terms at N2 of the nine runs. In it the run from start R
changed to the crystal-like orientation during the equilibration, and the plan then laid down, before the validation
runs and the production, that such a change counts as a conversion in rule 3. The definition of its time, read on the
equilibration and the production joined, was completed while the first production run from the crystal pose was
running and before any production run from a docked pose began. By rule 1 the simulation supports nothing if two or three of the three runs
from the crystal pose leave it (mean ligand RMSD of the last half above 2.0 Å). Rule 2 fixes the coordinate and the
verdict of a run, given under Analysis below. Rule 3 fixes the reading of the nine verdicts: if the rotated starts
convert in two or three of three runs while the others stay, the rotated orientation is read as one that does not
persist outside the scoring function; if no run converts, the simulation does not decide; mixed results are reported
run by run. The first case is the one found; the Limits below say how far its reading can be taken. No restraint or
bias acts on the orientation.

*System.* The protein is chain A of 6RL9 (1.0 Å) [@glockner2020influence], 258 residues from His3 to Lys261 with an
acetyl cap. The ligand is the sulfonamidate bonded to the zinc through its nitrogen N2, the form in which sulfonamides
bind this zinc [@linkuviene2018thermodynamic; @kovalevsky2018be]. Against the docking receptor of Section S13, Glu106
is a carboxylate (Section S15), His4, His107 and His122 carry their proton on Nδ, as found for the two buried ones
[@zoefisher2025proton], the amino terminus is capped (protein −2, complex −1), ten residues far from the ligand take
another conformer, and 215 of the 216 water sites are kept (site 510, of occupancy 0.32, clashes with a conformer
taken). The other hetero groups of the entry were removed, as in the
docking receptor: two further sulfanilamide molecules, neither at the zinc (one on the protein surface, one at the rim
of the active site), three glucose molecules and the mercury compounds at Cys206 and His64. His64 is neutral in its
outward conformer, one fixed state
of a residue that populates two positions and two neutral tautomers in solution [@raum2023energetics].

Three starts differ in the ligand: X is the crystal pose, and N and R are the best-scored near-native and ring-rotated
modes of the 20-seed AutoDock4Zn run from the DFT geometry in the receptor of Section S13, with Glu106 as the acid
(0.61 and 3.57 Å from the crystal ligand; scores −8.843 and
−8.823 kcal mol⁻¹), picked by a rule written after the modes had been tabulated. Both are poses of the neutral
molecule, so the two sulfonamide hydrogens were removed and one was rebuilt towards OG1 of Thr199. As docked, Zn–N2 is
2.21 Å in N and 2.24 Å in R (1.93 Å in the crystal); the bond of the model draws it in during the minimisation. For
start R the waters within 2.2 Å of a heavy atom of the rotated ligand were removed: waters 423, 516, 553 and 568, the
four that the docking test of Section S15 adds. The place of the ring in the crystal pose is then empty. The runs are
named after their starts, X1 to X3, N1 to N3 and R1 to R3; the run names N1 and N2 are not the atom labels N1 (aniline
nitrogen) and N2 (sulfonamide nitrogen) of the ligand, which are meant wherever a bond, a term or a distance is named.

*Force field and zinc model.* The protein is described by ff14SB [@maier2015ff14sb], water by TIP3P
[@jorgensen1983comparison], the ions by the Joung–Cheatham parameters [@joung2008determination] and the ligand by the
general Amber force field [@wang2004development] in its GAFF2 form (version 2.2.30). GAFF2 as assigned failed the test
of the free sulfonamidate against the DFT of this work, so the equilibrium values of S–N2, C4–S–N2 and S–N2–H and the
torsion about S–N2 were refitted in implicit water (a generalised Born model, OBC2).

The zinc is bonded to N2 and to the three histidine nitrogens in a model built by the MCPB.py procedure [@li2016mcpb].
A bonded zinc model with restrained electrostatic potential (RESP) charges has been tested in simulations of a carbonic
anhydrase II complex [@lin2010systematic], and a model with bonds from the zinc to the three histidines only has been
used for sulfonamide complexes of carbonic anhydrases [@taylor2023mm]. A bond between the metal and a bound ligand can
hinder the motion of the ligand, which led to a modified bonded model for a zinc β-lactamase [@eshtiwi2023modified];
here the bond to N2 is kept, so the runs test the orientation of the ring and not whether the ligand stays on the zinc.
The quantum chemistry of the model was done with
PySCF
at B3LYP/6-31G\*, with the set of Rassolov and co-workers for zinc [@rassolov19986], pure d and f functions and
Coulomb and exchange terms density-fitted in the def2-universal-JKFIT auxiliary basis [@weigend2008hartreefock]; this
is not the Cartesian 6-31G\* of Gaussian. The
RESP charges [@bayly1993well] of the zinc site (the zinc, the three histidines and the ligand) were fitted on
Merz–Singh–Kollman points [@singh1984approach] to the potential of a 101-atom capped model with the ligand in the
crystal pose alone, with equal charges on the chemically equal atoms of the ligand, and all eleven runs share them. Two
sets of zinc terms differ in the four Zn–N bonds and in 12 angles only. Set 2e has the empirical force constants of the
extended zinc Amber force field [@yu2018extended] with the Zn–N distances and N–Zn–N angles of 6RL9 as equilibrium
values; set 2s has force constants from Seminario's method [@seminario1996calculation] on a partial Hessian of the
optimised 64-atom model (zinc, three ethylimidazoles in place of the histidines, ligand). By a rule fixed beforehand
set 2s was primary unless it failed its validation; it failed one part, so the nine runs use set 2e. The terms at N2
are the same in both sets (@tbl:s-md-settings), and the zinc keeps the Lennard-Jones terms that MCPB.py assigns, the
set fitted to ion–oxygen distances by Li and co-workers [@li2013rational]. A zinc ion without bonds, held by
restraints on its four Zn–N distances, was the sensitivity model of the first plan. In tests of a few picoseconds it
did not keep the site (a carboxylate oxygen of Glu106 or a sulfonyl oxygen of the ligand came within 2.6 Å of the
zinc), so it was dropped before any run and the second bonded set took its place. Non-polarisable zinc models have
been compared at the ligand site of this enzyme and differ in the ligating atoms they prefer and in the water they
attract [@melse2023benchmarking].

*Protocol.* The runs were made with OpenMM 8.6.1 [@eastman2024openmm] with particle-mesh Ewald [@essmann1995smooth],
hydrogen masses repartitioned to 3.024 u for a step of 4 fs [@hopkins2015long], a Langevin integrator of the middle
scheme [@zhang2019unified] and a Monte Carlo barostat [@chow1995isothermal; @aqvist2004molecular]
(@tbl:s-md-settings). The heavy atoms of the solute are held in stages 1 to 3 of the equilibration. At the start of
stage 4 the ligand, the zinc site and the side chains are released, while the protein backbone stays held through
stages 4 and 5; nothing is held in stage 6. The production starts 1,300 ps after the release, so every time given
below as counted from the release lies in the equilibration. For start R stage 3 was extended by 500 ps if no water
oxygen lay within 2.5 Å of the crystal positions of the ring carbons at its end (once, in run R1). The production is
100 ns without restraint, with a solute frame every 10 ps: nine runs with set 2e, three per start with different
seeds, and two sensitivity runs from start R, which take the seeds of run R1. The mean kinetic temperature of the
production is 300.96 to 301.00 K with the thermostat at 300 K (299.8 to 300.9 K in stages 2 and 3, with a step of 2
fs); the cause of the excess was not tested.

*Analysis.* The orientation is read from d1, the distance from the aniline nitrogen N1 to CZ of Phe131, and d2, from
N1 to OD1 of Asn67: a frame is crystal-like if d1 is below 6.7 Å and d2 above 5.2 Å, and rotated-like in the opposite
case. A run has CONVERTED or STAYED if at least 90 % of the frames of its last 20 ns are in the state opposite to its
start or in that of its start, and is MIXED otherwise. The time of a conversion is the start of the first visit of at
least 1 ns in the opposite state, counted from the release; it is not a rate. Two programs written to one
specification agree on all 6,721,008 values compared, and a third reading, made without sight of their results, gives
every count, verdict and time alike.

*Tests of the force field.* @tbl:s-md-gates gives each test with its criterion. The repaired ligand terms meet the DFT
of the free sulfonamidate (G-L); that test was made with the AM1-BCC charges of the free anion and not with the RESP
charges of the runs, and it does not reach the coordinate which decides the simulations. Both zinc
sets passed the check of the site (G-Z), read on frames in which the site is still held. In its validation run of 5 ns
(G-V) set 2e kept the ligand at a mean RMSD of 0.55 Å with no oxygen within 2.6 Å of the zinc in its 50 all-atom
frames; set 2s (0.71 Å) had a sulfonyl oxygen there in 3 of the 50 where 2 were allowed. Set 2s is not a
five-coordinate model, but the rule made set 2e the primary set.

The force field of the nine runs failed two tests, both on the terms at the sulfonamide nitrogen, through which a
change between the two orientations must pass if the hydrogen bond to Thr199 is kept. The first (G-R) compared it with
B3LYP-D3(BJ)/6-31G\* energies [@grimme2010consistent; @grimme2011effect] of the 64-atom zinc-site model in vacuum: over
23 comparisons the root mean square difference is 3.08 kcal mol⁻¹ against a criterion of 1.0 kcal mol⁻¹. On a
constructed path between the two configurations, with no point relaxed, the planar nitrogen costs 13.76 kcal mol⁻¹ in
the force field and 6.62 kcal mol⁻¹ in quantum chemistry, so on this path the force field resists an inversion about
twice too strongly; the lowest path is not known, and the comparison reads the whole description of the nitrogen,
charges and 1–4 terms included, not the angle terms alone. The unrelaxed rotated placement is drawn 3.38 kcal mol⁻¹ too
high, and this error has the sign of the outcome reported below: a placement drawn too high is left more readily and
returned to less readily. Its sign rests on the dispersion term of the reference (−5.71 kcal mol⁻¹ without D3), against
which the basis-set superposition error of the 6-31G\* basis acts by an amount that was not determined; most of it
comes from the charges fitted to the crystal pose alone (+0.85 kcal mol⁻¹ with charges fitted to both placements), and
its size in the protein is not known. The second test was a refit: the best candidates have leave-one-out errors of 1.29
and 1.18 kcal mol⁻¹ against a limit of 1.0 kcal mol⁻¹ and miss the Zn–N2–S angle of the quantum-chemical optimum by
8.4° and 8.5° against a limit of 5°, so none was accepted and the nine runs use the terms as built. By the plan a
failed test was to stop the work. After the second failure we went on with the terms as built, a departure from the
plan that was recorded in it before the validation runs and the production, and we report both failures with the
direction of their errors.

*Results.* Rule 1 is met: none of the three runs from the crystal pose left it. The mean ligand
RMSD of the last half is 0.53 to 0.58 Å in runs X1 to X3, against the 2.0 Å of the rule, and 0.53 to 0.67 Å in runs N1
to N3 (@tbl:s-md-runs, @fig:s-md-rmsd). In these six runs the hydrogen bond N2–H···OG1(Thr199) is present in 99.77 to
99.86 % of the production frames, none of the 60,000 frames is rotated-like, and all six have the verdict STAYED
(panels a to f of @fig:s-md-orientation).

The three runs from the rotated pose left it during the equilibration. The configuration at the sulfonamide nitrogen
(the side of the plane of its three neighbours on which it lies) changed to that of the crystal pose 260, 360 and 140 ps
after the release of the ligand in R1, R2 and R3, and the first frame that is not rotated-like came at 270, 360 and 150
ps, in all three runs in stage 4, with the protein backbone still held. By the rule, fixed before any production run
from the rotated pose began, the conversions are at 540, 370 and 160 ps. At the change the sulfur moves about the Zn–N2 bond
while the hydrogen keeps its side (in R1 the torsion NE2(His94)–Zn–N2–S goes from −8° to +78° and NE2(His94)–Zn–N2–H
from −144° to −155°), and the hydrogen bond to Thr199 is kept: an inversion at the nitrogen, once in each run. The
production of R1 to R3 holds no rotated-like frame among 30,000, and its ligand is as in the runs from the crystal pose
(mean RMSD of the last half 0.50 to 0.53 Å, hydrogen bond in 99.82 to 99.91 % of the frames).

The two sensitivity runs enter no count. With zinc set 2s the run from start R is CONVERTED (nitrogen changed at 30
ps, time of the rule 890 ps); it varies the bond and angle terms of the zinc only and tests neither the terms at N2
nor the charges. With the candidate terms at N2 the run is MIXED: 74.35 % of the frames of its last 20 ns are
crystal-like, 9.15 % rotated-like and 16.50 % intermediate, and the rotated-like frames per 20 ns grew over the run
(23, 11, 50, 126 and 183), so 100 ns give no share that could be read as a population. A rotated-like frame of this run
is mostly not the docked pose: of its 393 rotated-like production frames 179 have the nitrogen in the configuration of
the crystal pose, and the ligand lies within 1.0 Å of the docked rotated pose in 80 frames. The candidate is no clean
bracket either: at the relaxed near-planar nitrogen it errs by −1.9 kcal mol⁻¹ where the terms as built err by +2.4
kcal mol⁻¹, but it also puts the rotated placement 11.07 kcal mol⁻¹ above the crystal-like one (14.22 as built, 10.84
in quantum chemistry).

In none of the eleven runs does a water or protein oxygen come within 2.6 Å of the zinc in an all-atom frame, and
neither carboxylate oxygen of Glu106 comes within 2.6 Å of it in a solute frame (smallest Zn–OE1 distance 3.25 Å). A
sulfonyl oxygen of the ligand, which interacts with the zinc through the scaled 1–4 terms of the bonded model, is
within 2.6 Å in 0.07 to 0.17 % of the solute frames of the nine runs and in 2.54 and 2.55 % in the two sensitivity runs.

One quantum-chemical calculation without a force field bears on the same question. In the 64-atom zinc-site model with
the zinc and the imidazoles fixed, a B3LYP-D3(BJ)/6-31G\* relaxation of the ligand from the rotated placement inverts
at the nitrogen and ends 0.12 Å and −0.01 kcal mol⁻¹ from the end point reached from the crystal placement. It is one
relaxation in vacuum, in a model without Thr199, pocket wall or water that also moves the crystal-like ligand 2.29 Å
from its crystal placement.

*Limits.* These are the results of one force field, whose terms at the sulfonamide nitrogen failed both tests made of
them; one of their errors, the rotated placement drawn 3.38 kcal mol⁻¹ too high, has the sign of the outcome. All
three changes happened while the protein backbone was held to its crystal coordinates; whether the rotated pose would
last longer in a protein free to adapt was not tested. That the crystal-like state was kept and not left again belongs
to the terms as built, which hold the nitrogen too stiffly (nearly planar in 0 to 1 of the 10,000 production frames of
a run); rule 1 cannot detect a force field that holds the crystal pose too firmly. The rehearsal run of 5 ns from start
N, with the same topology, had 3 rotated-like frames among 500, and with the candidate terms the run from start R was
MIXED.

The zinc geometry of set 2e is put in and is no result: the mean Zn–N2 distance of 1.924 to 1.928 Å (1.931 Å in 6RL9)
follows from the equilibrium values. Other parts of the site depart from the crystal. The zinc lies 0.46 to 0.53 Å out
of the plane of the His119 ring (0.15 Å in 6RL9) and S–C4 is 1.821 to 1.822 Å (1.771 Å), neither of which is put in;
S–N2–Zn is 113.7 to 114.5°, although its equilibrium value in the force field is that of 6RL9 (117.6°).

Start R was built without the four waters under the rotated ligand and with the place of the crystal ring empty, which
could make a conversion look easier than it is: the hydration check was met in every run from R (in R3 by one water
oxygen 2.312 Å from a ring-carbon position, with the cut-off at 2.5 Å), but in R2 the last all-atom frame before the
change, 300 ps after the release, shows no water oxygen in that place (nearest 2.901 Å). The site around the ligand is
not converged. From 20 to 100 ns His64 lies deep in the site, its ring nearer to the zinc than in either inward
conformer of 6RL9, for 96.8, 87.0 and 89.3 % of the time in R1, R2 and R3, against 0.0 to 52.9 % in the six other runs,
and the place of crystal water 516 is occupied in 1.1 to 8.2 % of the all-atom frames, against 30.4 to 81.8 %. The limit
that defines the inner position was chosen after looking at the data, and chance cannot be told from a memory of the
start. That the runs from R end like those from the crystal pose holds for the ligand only.

The design is small: three runs per start, which share their start structure, 100 ns each, one force field, docked
starts with a constructed hydrogen, and a coordinate of two distances that can separate two states only. The times are
single events seen in frames 10 ps apart and the shares are not populations, so no rate, no barrier and no free energy
follows. No end-point free energy (MM/GBSA or MM/PBSA) was computed: the question is the orientation of the ring and
not an affinity, and runs with terms that failed their tests could not carry such a number.

Three things point the same way: the simulations, the relaxation of the bare zinc-site model, and the docking with the
four waters of Section S15. None reproduces the near-degeneracy of the two orientations that AutoDock4Zn shows in the
water-free receptor, and each is in line with the crystal structure. None decides. The simulations rest on terms that
failed their tests, the relaxation lacks protein and water, the docking excludes the rotated pose by construction, and
the three are not independent of one another: the force field was judged against the same zinc-site model, and start R
was built by removing the waters that the docking test puts back. The 1.0 Å crystal structure, with the zinc-bound
ligand in one conformer at full occupancy, remains the primary evidence for the orientation of the ring.

The ZIP file of the Supplementary material holds the prepared protein with its zinc and crystal water, the three
starts, the ligand parameter files, the terms at N2, the zinc parameter file of set 2e, the tleap inputs and, for each
system, a build record that lists every zinc term and every charge of the site, the protocol, the scripts, the
analysis output of every run, the records of the quantum-chemical tests, a time series of each of the eleven runs
every 10 ps and the last solute frame of each, with the receptors and comparison files of the docking tests of Section
S15. The built topology files and the trajectories are not in it.

![Orientation of the ring in the eleven molecular dynamics runs of 100 ns of the hCA II–sulfanilamide complex (6RL9,
sulfonamidate bonded to the zinc; ff14SB, GAFF2 and a bonded zinc model, OpenMM): the distances d1 from the aniline
nitrogen N1 to Phe131 CZ (blue) and d2 from N1 to Asn67 OD1 (orange) against time. (a–i) The nine runs with zinc set
2e: three from the crystal pose (X1 to X3), three from the docked near-native pose (N1 to N3) and three from the
docked ring-rotated pose (R1 to R3). (j, k) The two sensitivity runs from the ring-rotated pose, with zinc set 2s and
with the candidate terms at N2 that were not accepted; they enter no count. Time zero is the start of the production.
The grey strip, on a stretched scale, is the 1.3 ns of the equilibration after the release of the ligand; its vertical
broken lines mark the ends of stages 4 and 5, 400 and 800 ps after the release, and the protein backbone is held in
both stages. In the production the pale line is every frame (10 ps) and the strong line the mean over 1 ns. The
horizontal lines are the midpoints of d1 and d2, 6.7 and 5.2 Å: a frame is crystal-like with d1 below and d2 above its
midpoint and rotated-like in the opposite case. The word at the top right of a panel is the verdict of the run over
its last 20 ns. A filled triangle marks the conversion by the rule fixed beforehand, the start of the first visit of
at least 1 ns in the crystal-like state; the open triangle in (k) marks the same frame in the run whose verdict is
MIXED and is not a conversion. Every conversion lies in the equilibration; the times are single events and give no
rate. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_md_orientation.png){#fig:s-md-orientation pdf-max-height=16.5}

![Ligand RMSD in the eleven molecular dynamics runs of 100 ns of the hCA II–sulfanilamide complex: the 11 heavy atoms
of the ligand against the crystal pose of 6RL9 after the fit of the 258 Cα atoms, smallest of the four equivalent
labellings. Panels, time axis, grey strip, pale and strong lines, verdicts and triangles as in @fig:s-md-orientation:
(a–i) the nine runs with zinc set 2e, (j, k) the two sensitivity runs from the ring-rotated pose, which enter no
count. The horizontal broken line is the 2.0 Å of the rule fixed beforehand for the runs from the crystal pose, which
reads the mean of the last half of a run. The open triangle in (k) is not a conversion: the verdict of that run is
MIXED. Drawn by a script written and run with the generative AI tool Claude (Section S1).](figures/figS_md_rmsd.png){#fig:s-md-rmsd}

<!-- 11 Oct 2026: Figs. S9 and S10 stand before Tables S38 to S40, with a page break after them, so that the two
     full-page figures do not cut the note of Table S40 (PDF only; the Word build drops the raw line). -->

\clearpage

<!-- include: tables/tableS33_md_settings.md -->

<!-- include: tables/tableS34_md_runs.md -->

<!-- include: tables/tableS35_md_gates.md -->

# References

::: {#refs}
:::
