---
# Main text of the SAPH paper. Target journal since 12 Oct 2026: Journal of Molecular Graphics and Modelling (JMGM,
# Elsevier); before, from 5 Oct 2026, the Journal of the Indian Chemical Society (JICS), which declined the manuscript.
# The molecular dynamics of the hCA II complex and the two docking tests (Sections S15 and S16) entered on 11 Oct 2026.
# JMGM guide of 12 Oct 2026: abstract at most 250 words; 1 to 7 keywords; highlights and graphical abstract required;
# every figure also as a separate file; AI tools disclosed in each figure caption and in the declaration.
# Since 8 Oct 2026 this file is the source and is edited directly (the conversion script apply_jics.py is retired).
# Build: tools/paper-env/bin/python paper/_shared/tools/build.py main   (journal profile jics, the default)
# Writing rules: paper/_shared/STYLE.md. References: only keys from lit/candidates.jsonl (lit.py list / lit/REFS.md).
# Numbers: every number must come from paper/main/data/facts.json (numbers_audit.py main).
# Figure captions are copied from figures/captions.md; keep the two in step.
# Length: at most 8150 words of running prose (abstract and body paragraphs; lint_prose.py prints the count); the cap
# was 8000 until 8 Oct 2026. Sections S15 and S16, Tables S36-S40 and Figs. S9 and S10 of the Supplementary material
# are cited here by number, typed by hand: check them against paper/si/manuscript.md after any change there.
# Journal format (JMGM as JICS: single anonymized review, subscription route). Elsevier numeric citations "[1]" and the
# reference list of elsevier-with-titles.csl. Front matter: title, authors with affiliation letters, corresponding
# author, e-mail addresses, ORCID iDs, abstract, at most six keywords. Back matter: CRediT authorship contribution
# statement, Declaration of competing interest, Funding, Data availability, Appendix A. Supplementary material,
# Declaration of generative AI (the last section before the reference list, as Elsevier asks), References. Highlights
# (paper/main/highlights.md) and the graphical abstract (figures/graphical_abstract.tif) are separate files.
title: "Proton transfer in 4-sulfamoylanilinium dihydrogen phosphate: a combined molecular and periodic DFT study with molecular docking"
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
# PDF only (crossref.lua): cap figure height at 15.5 cm; Fig. 10 and its caption overflow an A4 page at full size. The
# Word file keeps full size.
pdf-max-figure-height: 15.5
abstract: |
  Protonated sulfanilamide and phosphoric acid have almost equal pK~a~ values (ΔpK~a~ ≈ 0), yet they crystallise as the
  salt 4-sulfamoylanilinium dihydrogen phosphate (SAPH). With B3LYP-D3(BJ)/6-311++G(d,p) the gas-phase ion pair
  optimises to a neutral sulfanilamide···H~3~PO~4~ adduct. In a water continuum (IEF-PCM) the salt is the higher of two
  minima, 0.54 kcal mol⁻¹ above the neutral pair, with the proton's ground level near the top of the 0.72 kcal mol⁻¹
  barrier (one-dimensional model). Four more methods favour the neutral pair, and along a dielectric ladder a salt
  minimum first appears between *ε* = 4 and 8. In the periodic PBE-D3(BJ) crystal, starts with the proton on the
  phosphate relax to the salt, and along a constrained N–H path the energy rises at every point, to 51.6 kJ mol⁻¹ per
  formula unit at 1.70 Å (50.0 with PBE0-D3(BJ)). The transfer weakens the two other anilinium N–H···O bonds and the
  2.500 Å anion–anion bond, yet five ion pairs cut from the lattice prefer the neutral form in vacuum. The crystal,
  hydrogen bonds and ionic lattice together, holds the proton. A gas-phase ion pair is no model of this salt. At pH 7.4
  SAPH releases neutral sulfanilamide. Docking tested by redocking keeps its zinc and Thr199 contacts in carbonic
  anhydrase II, but the top-ranked ring orientation depends on input geometry and Glu106 protonation; in nine molecular
  dynamics runs only the crystal pose persists, though force-field terms failed two quantum-chemical tests.
keywords:
  - "Salt–cocrystal continuum"
  - "Charge-assisted hydrogen bond"
  - "Sulfanilamide"
  - "Periodic DFT"
  - "Molecular docking"
  - "Molecular dynamics"
---

# 1. Introduction {#sec:intro}

Sulfanilamide is the parent compound of the sulfonamides, the first class of synthetic antibacterial drugs
[@ovung2021sulfonamide; @venkatesan2023molecular]. Sulfa drugs compete with *p*-aminobenzoic acid (pABA) for its pocket
in dihydropteroate synthase (DHPS), a pocket that two conserved loops of the enzyme build during catalysis
[@yun2012catalysis], and a bound sulfonamide is turned into a false product that blocks the formation of
7,8-dihydropteroate [@snoke2025development]. Resistance comes from mutations in *folP*, the gene of DHPS, or from
plasmid-borne *sul* genes encoding Sul1, Sul2 and Sul3, divergent enzymes that keep binding pABA but reject sulfonamides
through a reorganised pABA region with a Phe–Gly sequence [@venkatesan2023molecular]. A Sul enzyme is therefore the
natural counterpart to DHPS.

Sulfanilamide also inhibits carbonic anhydrases: human carbonic anhydrase II (hCA II) has been crystallised with it at
1.0 Å resolution [@glockner2020influence], and the complex has been simulated by molecular dynamics
[@nguyen2025understanding]. The sulfonamide nitrogen binds the zinc after losing its proton, so measured binding data
combine the intrinsic binding step with the protonation equilibria of ligand and protein [@linkuviene2018thermodynamic;
@jiang2022quantum]. With a 1.1 Å structure of the acetazolamide complex [@sippel2009high] and a curated database of
measured sulfonamide binding [@linge2023plbd], hCA II is where a docking protocol for sulfanilamide can be checked
against both a crystal pose and an experimental order of binding strength.

An acid and a base form a salt or a cocrystal according to where one proton sits in one hydrogen bond, and between the
two limits lies a continuum [@childs2007saltcocrystal; @radekstocek2024hydrogenbond]. The usual guide is the ΔpK~a~
rule, with ΔpK~a~ = pK~a~(protonated base) − pK~a~(acid): below −1 a cocrystal is expected, above 4 a salt, and in
between both are found [@cruzcabeza2012acidbase]. Cruz-Cabeza and co-workers traced the rule to solvation: the lattice
plays the part of the solvent, and in the intermediate range a salt forms when ancillary hydrogen bonds can stabilise
the ions, so that the crystal, and not only the pair of partners, decides [@cruzcabeza2022role].

X-ray diffraction locates hydrogen atoms poorly [@malaspina2021advanced], so proton transfer is often inferred from
heavy-atom bond lengths. Hirshfeld atom refinement of X-ray data can place the proton directly, and in very short
hydrogen bonds nuclear quantum effects can smear it between donor and acceptor [@radekstocek2024hydrogenbond].
Calculations err in opposite directions. An isolated acid–base adduct can stay neutral in the gas phase
[@omalley2021crystallization] and transfer its proton only once a polarisable continuum is switched on
[@fedorov2025influence]. Periodic density functional theory (DFT) calculations with generalised-gradient-approximation
(GGA) functionals err the other way, because their delocalisation error stabilises the salt and can move a proton that
should stay put [@leblanc2018pervasive; @bridgersanguiano2026salt; @streek2026best]. Hybrid functionals reduce the
problem [@bridgersanguiano2026salt; @bal2024charting], and a benchmark against neutron and solid-state NMR proton
positions found efficient periodic methods unreliable for N···O contacts shorter than about 2.57 Å [@abramov2024it].

Salts of the 4-sulfamoylanilinium cation are known with hydrochloric acid [@zaoualizgolli20104], where proton transfer
is not in doubt; a citrate and a sulfamate have been reported [@muthuselvi2022molecular; @wojnarska2018crystal].
Phosphoric acid is far weaker than hydrochloric acid. Its first pK~a~ lies close to the pK~a~ of the anilinium group of sulfanilamide, so the pair sits at
ΔpK~a~ ≈ 0 [@greenbaum1955comparison; @bates1951first].

4-Sulfamoylanilinium dihydrogen phosphate (SAPH) is the salt of these two partners. Its crystal structure (294 K, *Pc*,
*Z* = 2) consists of cation chains held by N–H···O bonds of the sulfonamide group and anion chains held by O–H···O
bonds, both along *a* [@muthuselvi2014crystal]. The deposited data establish the salt without relying on hydrogen
positions: C1–N1 is 1.463 Å, as in an anilinium ion, and the phosphate has two short (1.498 and 1.503 Å) and two long
(1.558 and 1.577 Å) P–O bonds, the pattern of H~2~PO~4~⁻. The same file holds the shortest hydrogen bond of the crystal,
O4–H4A···O5 (2.500 Å) between neighbouring anions, and a weak R²₂(8) ring, closed by N1–H1A···O5 and a long C6–H6···O4
contact, which the 2014 report does not discuss. A later study presented SAPH as a nonlinear optical (NLO) material,
with FT-IR, FT-Raman and UV–Vis spectra and Hartree–Fock and B3LYP/6-311++G(d,p) calculations
[@muthuselvi2017structural]. To our knowledge no periodic calculation or docking study of the salt has been reported.

At ΔpK~a~ ≈ 0 the proton on N1 cannot be taken for granted in any model of SAPH. An optimisation of the isolated ion
pair may return a neutral sulfanilamide···H~3~PO~4~ adduct, and then every property computed for it belongs to a
species that is not in the crystal. The crystal does hold the proton on N1, and it is not obvious whether bulk
screening or the lattice itself does the holding. The proton also sets what a biological test of SAPH measures: once
the salt dissolves at pH 7.4 the anilinium group gives up its proton, and what reaches DHPS, hCA II or a Sul enzyme is
neutral sulfanilamide, with its sulfonamidate anion as the form that binds zinc [@linkuviene2018thermodynamic].

Here we follow H1A from the isolated ion pair to the solid and find that the crystal, not screening alone, keeps it on
N1. We optimise the ion pair in the gas phase, in a continuum model of water and along a ladder of dielectric constants,
scan the N1–H1A distance from both ends, and test the result with further functionals, density-fitted second-order
Møller–Plesset perturbation theory (DF-MP2), the thermochemistry of both minima and a one-dimensional proton model. In
the *Pc* cell we relax salt and neutral starts and follow the transfer along a constrained N1–H1A path with PBE and
PBE0, and a cluster cut from the periodic cells gives density descriptors for all nine hydrogen bonds. Finally, we dock
neutral sulfanilamide and its sulfonamidate into DHPS, hCA II and Sul2 with a protocol tested by redocking and
cross-docking, and probe its limits at hCA II by molecular dynamics.

# 2. Computational methods {#sec:methods}

## 2.1. Crystal data and models {#sec:models}

Every model was cut from the deposited structure of SAPH [@muthuselvi2014crystal], Cambridge Crystallographic Data
Centre (CCDC) deposition 1016950: *Pc*, *a* = 4.8041, *b* = 10.8564, *c* = 10.3862 Å, *β* = 101.067°, *Z* = 2, Mo Kα
radiation (*λ* = 0.71073 Å), 294 K, R1 = 0.021 for 2502 reflections with *I* > 2σ(*I*). The N–H hydrogens were refined
freely, the C–H hydrogens rode at 0.93 Å and the two O–H hydrogens were restrained to 0.82 Å, so the O–H positions are a
modelling choice and the N–H bonds come out short [@malaspina2021advanced].

The primary model is the 27-atom ion pair (M1): the cation and the anion at *x*, *y*, *z* − 1, joined by N1–H1A···O5 and
closed into the R²₂(8) ring by C6–H6···O4 (panel a of @fig:scheme). Its two ions (M2, M3), the neutral products of a
proton transfer (sulfanilamide, M4, and H~3~PO~4~, M6), the sulfonamidate anion deprotonated at N2 (M5, the zinc-binding
form) and urea (M7, the NLO reference) complete the molecular set (@tbl:models). The periodic model (M9) is the full
*Pc* cell, 54 atoms.

<!-- include: tables/table01_models.md -->

We measure the position of the proton by *δ* of @eq:delta, which is negative with H1A on N1 (the salt) and positive
with H1A on O5 (the neutral pair):

$$\delta = d(\mathrm{N1{-}H1A}) - d(\mathrm{H1A{\cdots}O5})$$ {#eq:delta}

The X-ray value, *δ* = −1.02 Å with the refined H1A, is a lower bound because X-ray N–H bonds are short. C1–N1 and
P1–O5 carry the same information without using H1A and are reported beside *δ*.

<!-- Caption copied from figures/captions.md (group A, fig01_scheme, 28 Sep 2026); keep the two in step. -->

![4-Sulfamoylanilinium dihydrogen phosphate and the hydrogen bonds of its crystal (CCDC 1016950, space group *Pc*,
*Z* = 2, X-ray data at 294 K). (a) The ion pair used as the primary model: the cation and the anion at *x*, *y*,
*z* − 1 are joined by N1–H1A···O5 (bond 1) and C6–H6···O4 (bond 9), which close an R²₂(8) ring (shaded; the ring is
drawn regular, not to scale). CIF atom labels in grey. (b) The nine hydrogen bonds of the deposited structure,
numbered in the order of the CIF, with their D···A distances and D–H···A angles and the panel in
which each is drawn. (c) One hydrogen-bonded layer, made of ions related by lattice translations only, viewed along
*b*: bonds 1 and 9 join each cation to its anion along *c*, bonds 5 (N2–H2B···O2) and 7 (O4–H4A···O5, 2.500 Å, the
shortest) run along *a*, and bonds 6 and 8 link neighbouring pairs diagonally. (d) Packing viewed along *a*: layers
related by the *c*-glide are joined by bonds 2 and 3 (N1–H1B···O6 and N1–H1C···O6, whose acceptors lie one
translation apart along *a* and so coincide in this projection) and bond 4 (N2–H2A···O3); bonds 5–8 run along *a* and
are seen end-on here. Teal dashes: the two bonds kept in the 27-atom ion-pair model; brown dashes: bonds present only
in the crystal; yellow halo: the ion pair modelled by DFT; shading: R²₂(8) rings; thin grey lines in (d): unit-cell
edges. Each bond type is numbered once per panel. Atoms: C dark grey, N blue, O red, S yellow, P orange, H light
grey. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig01_scheme.png){#fig:scheme}

## 2.2. Level of theory and environments {#sec:level}

Molecular calculations used PySCF 2.14.0 [@sun2020recent]; its libraries and the other programs are listed, with their
references, in Section S1 of the Supplementary material. B3LYP [@becke1993density; @lee1988development] was combined
with Grimme's D3 correction with Becke–Johnson damping, two-body terms only [@grimme2010consistent; @grimme2011effect].
The basis is 6-311++G(d,p) [@krishnan1980self; @frisch1984self] in spherical form, 445 functions for the ion pair.
Coulomb and exchange terms were density-fitted, the integration grid was PySCF level 4, and the self-consistent field
(SCF) was converged to 10⁻⁸ Ha in the continuum optimisations and scans and to 10⁻⁹ or 10⁻¹⁰ Ha elsewhere (Table S1).
This level is written B3LYP-D3(BJ)/6-311++G(d,p) below.

Gas phase means the isolated model. IEF-PCM (water) is the integral-equation formalism of the polarisable continuum
model [@cances1997integral] with *ε* = 78.3553 and a switching/Gaussian (SWIG) discretised surface built from Bondi
radii (Section S1). A
continuum makes no hydrogen bonds [@herbert2021dielectric], and that is why it is used here: it isolates bulk
dielectric screening, whereas explicit water would add specific solvation. The ion pair was also optimised at *ε* = 4,
8 and 20 from the water optimum with H1A on N1; *ε* = 1 is the gas-phase optimisation, which started from the X-ray
structure, and along this dielectric ladder only the side on which H1A ends is compared. B3LYP without dispersion in
IEF-PCM (water), and B3LYP-D3(BJ)/def2-TZVP for the ion pair and sulfanilamide (Section S6), test the level itself.

## 2.3. Optimisations, frequencies, proton scan and periodic model {#sec:protocol}

Optimisations used geomeTRIC 1.1.1 with Gaussian's default criteria; two optimisations that crept along a nearly free
inter-ion rotation were accepted on their forces (Section S1). Harmonic frequencies of both IEF-PCM (water) minima of
the ion pair, salt and neutral, come from central differences of 162 analytic gradients each, with IR intensities from
the dipole derivatives; neither Hessian has an imaginary frequency (lowest modes 20 and 18 cm⁻¹). Frequencies set
against fundamentals are scaled by 0.9688 (Section S4), and thermochemistry is harmonic (298.15 K).

In IEF-PCM (water) the proton-transfer profile is a relaxed scan of *d*(N1–H1A) on a 0.10 Å grid, run outwards from
the salt minimum (chain A) and inwards from the neutral minimum (chain B, reached from H1A placed 1.00 Å from O5) so
that hysteresis would show; the lower envelope is the profile (convergence in Section S1).

Density functionals can disagree on whether an acid–base system is a salt or a cocrystal [@radekstocek2024hydrogenbond;
@seyedraoufi2026engineering]. The salt end point (*d*(N1–H1A) = 1.098 Å), the 1.25 Å point and the neutral end point
(1.631 Å) were therefore recomputed in IEF-PCM (water) with ωB97X-D [@chai2008long], M06-2X [@zhao2008m06], PBE0-D3(BJ)
[@adamo1999reliable] and frozen-core DF-MP2 [@mller1934note] with the aug-cc-pVTZ-RI fitting basis (Section S5), the
last with the continuum acting only through the reaction field of the Hartree–Fock density. A one-dimensional
Schrödinger equation for the proton on the relaxed profile (sinc discrete variable representation, N–H reduced mass
0.940 u; Section S1) places the ground level relative to the barrier; it ignores the other modes and multidimensional
tunnelling.

For the crystal we used Quantum ESPRESSO 7.5 [@giannozzi2017advanced], PBE [@perdew1996generalized] with two-body
D3(BJ), SSSP efficiency 1.3.0 pseudopotentials and fixed occupations, relaxing the atoms in the X-ray cell by BFGS to
forces below 5 × 10⁻⁴ Ry Bohr⁻¹ (70/560 Ry, 2×1×1 Monkhorst–Pack grid; pseudopotentials in Section S1, convergence tests
in Section S7). Relaxations started (i) from the deposited structure, (ii) with H1A moved to 0.98 Å
from O5 on the N1···O5 line in both formula units, and (iii) with H1A moved in formula unit 1 only, without symmetry.
Heavy-atom root-mean-square deviations (RMSDs) are given after removing the mean shift (the origin of *Pc* floats in the
*a*–*c* plane).

Because the delocalisation error of GGA functionals favours the salt [@leblanc2018pervasive;
@bridgersanguiano2026salt; @husak2022potential], the transfer was also followed along a constrained path. The N1–H1A distance of both
formula units was held at 1.10, 1.25, 1.40, 1.55 and 1.70 Å while every other coordinate relaxed, O5 and the H1A–O5
separation included; Quantum ESPRESSO applies such constraints only with damped dynamics, which ran until the
constrained forces fell below 10⁻³ Ry Bohr⁻¹, and the constraint force gives the slope of the profile. Each point and
cell (i) received a hybrid correction, *E*(PBE0-D3(BJ)) − *E*(PBE-D3(BJ)) from single points with 25 % exact exchange
and norm-conserving SG15 pseudopotentials (Section S1), as in hybrid-functional treatments of acid–base crystals
[@bal2024charting]; cutoff and *q*-mesh checks are in Section S7.

## 2.4. Wavefunction and crystal-surface analyses {#sec:analyses}

Frontier orbitals and densities of states come from the final SCF of each model and were projected on the two ions by
Hirshfeld orbital composition [@lu2012multiwfn]. The highest occupied and lowest unoccupied molecular orbitals (HOMO and
LUMO) of the ion pair are split differently between the ions, which gives its orbital gap charge-transfer character, so
vertical ΔSCF energies replace it: $I = E_{N-1} - E_N$ and $A = E_N - E_{N+1}$, with stability-checked unrestricted
doublets and equilibrium IEF-PCM for every charge state, so that the solvent relaxes around each ion. Hardness $\eta =
(I - A)/2$, electronegativity $\chi = (I + A)/2$, softness $S = 1/(2\eta)$ and electrophilicity $\omega =
\chi^2/(2\eta)$ follow [@geerlings2020conceptual]; Koopmans-type values serve only for comparison
(Section S9).

Twenty singlet excitations of the ion pair, the cation and neutral sulfanilamide were computed at their B3LYP-D3(BJ)
geometries by linear-response time-dependent DFT (TD-DFT) with the range-separated hybrid CAM-B3LYP [@yanai2004hybrid]
and the 6-311++G(d,p) basis, in IEF-PCM (water) with non-equilibrium solvation (Section S11). B3LYP was run alongside to
expose the spurious charge-transfer states of a global hybrid [@maitra2022double], and the inter-ion charge transfer of
each state is its charge-transfer (CT) number [@plasser2012analysis]. Dipole moments, polarisabilities and first
hyperpolarisabilities were obtained by finite fields (Section S12).

Multiwfn 2026.9.20 [@lu2012multiwfn; @lu2024comprehensive] located the critical points of *ρ* (quantum theory of atoms
in molecules, QTAIM). At each intermolecular bond critical point we report *ρ*, ∇²*ρ* and the total energy density *H*,
negative for partly covalent bonds [@cremer1984chemical]. The ion pair holds two of the nine hydrogen bonds of the
crystal, so five ion pairs
that hold all nine, eight at the central pair (135 atoms, M8 in @tbl:models), were cut from cell (i) and from the
periodic configuration with O5–H1A held at 0.98 Å (@sec:proton), and their B3LYP-D3(BJ)/6-31+G(d,p) densities
(gas-phase single points) were analysed in the same way. Non-covalent regions were mapped with the reduced density
gradient (Section S8) and with the independent gradient model based on Hirshfeld partition (IGMH), whose
inter-fragment *δg*^inter^ was integrated and split into atom-pair contributions [@lu2022independent]. Atomic charges
are natural population analysis (NPA) charges from PySCF's natural atomic orbitals [@reed1985natural] and Hirshfeld,
atomic dipole moment corrected Hirshfeld (ADCH) and charge model 5 (CM5) charges from Multiwfn; ADCH is used only for
fragment sums. The cluster basis, the three Hirshfeld-based charge schemes and the electrostatic potential on the *ρ* =
0.001 a.u. density surface are described, with their references, in Section S9.


## 2.5. Molecular docking and molecular dynamics {#sec:docking}

The docked ligands are the species present at pH 7.4: neutral sulfanilamide (M4) is the primary ligand, its
sulfonamidate (M5) is the zinc-binding form [@linkuviene2018thermodynamic; @samanta2016prediction], and the cation (M2) is only a sensitivity
case. The controls were sulfamethoxazole and acetazolamide (neutral and anionic), *p*-aminobenzoate and each receptor's
crystal ligand. The sulfanilamide species were docked from their B3LYP-D3(BJ)/IEF-PCM (water) geometries (stationary
on forces, without a frequency check; only bond lengths and angles enter the docking), rigid apart from torsions;
controls and validation ligands were ETKDG conformers [@riniker2015better] minimised with MMFF94s [@halgren1999mmff].

The receptors [@burley2023rcsb] were DHPS of *Yersinia pestis* with sulfamethoxazole, pterin pyrophosphate and Mg²⁺
(3TZF, chain B, 2.1 Å [@yun2012catalysis]); DHPS of *Escherichia coli* with sulfanilamide (1AJ0, 2.0 Å
[@achari1997crystal]); hCA II with sulfanilamide (6RL9, 1.0 Å [@glockner2020influence]) and with acetazolamide (3HS4,
1.1 Å [@sippel2009high]); and Sul2 with 7,8-dihydropteroate, pyrophosphate and Mg²⁺ (7S2K, 1.74 Å
[@venkatesan2023molecular]). PDBFixer (OpenMM 8.6.1 [@eastman2024openmm]) rebuilt missing side chains and internal gaps,
and PROPKA 3.5.1 [@olsson2011propka3] set the protonation states at pH 7.4; run without the zinc, it made Glu106 of hCA
II the neutral acid (pK~a~ 8.91; 5.33 for 6RL9 with the zinc; Section S15). The three zinc histidines of hCA II, which
keep one neutral tautomer at every pH studied [@zoefisher2025proton], were assigned from their coordination, and His64
was set to HIE (Nε–H). Waters were kept only within 2.6 Å of Mg²⁺, and Meeko 0.8.0 wrote the PDBQT files with Gasteiger
charges [@gasteiger1980iterative]; Section S13 gives the remaining details. The ligands above, ciprofloxacin and the
crystal inhibitor were docked with Vina into a control receptor that is not a sulfonamide target, the ATP site of *E.
coli* gyrase B (4DUH, 1.50 Å [@brvar2012structure]), prepared as above and, separately, with its crystal water kept.

AutoDock Vina 1.2.7 [@trott2010autodock; @eberhardt2021autodock] was run with the Vina and Vinardo [@quiroga2016vinardo]
scoring functions, and for hCA II also with AutoDock4Zn [@santosmartins2014autodock4zn] on AutoGrid 4.2.7 maps
[@morris2009autodock4]. Its zinc potential was built to let an uncharged sulfonamide –NH~2~ group coordinate the metal
[@santosmartins2014autodock4zn], and in validation it kept the sulfonamidate on the zinc of 6RL9 (2.21 Å) but left the
acetazolamide anion 4.86 Å away in 3HS4, so it was used in production only with neutral sulfonamides. Each job used a 20
Å cube centred on the crystal ligand, exhaustiveness 32 and 10 seeds (20 for the two sulfanilamide AutoDock4Zn jobs
whose poses are discussed).

Validation required the protocol to reproduce the crystal ligands, judged by symmetry-corrected heavy-atom RMSDs in
the receptor frame (spyrmsd 0.9.0 [@meli2020spyrmsd]). The strict gate was a mean top-pose RMSD of at most 1.0 Å with at
least half of the seeds within 1.0 Å; the conventional 2.0 Å criterion [@buttenschoen2024posebusters] was evaluated
beside it. Both also required the PoseBusters 0.6.5 physics checks and, for hCA II, zinc and Thr199 contacts (Section
S13). Redocking alone is not a sufficient test [@kittelson2026reproducibility], so the ligands of 6RL9 and 3HS4 were
also cross-docked after superposing 257 Cα atoms, and acetazolamide, the stronger hCA II inhibitor [@linge2023plbd], had
to score better than sulfanilamide.

Top poses of all seeds were clustered at 2.0 Å and the cluster representatives profiled with PLIP [@adasme2021plip] and
ProLIF 2.2.2 [@bouysset2021prolif]. Scores are given in kcal mol⁻¹ as mean ± standard deviation (SD) over the seeds.
They are empirical docking scores, not free energies of binding, and only their order against the controls is read
[@paggi2024art; @yang2022proteinligand]. The C3–C4–S1–N2 torsion of the poses is compared with that of the ligands in
the deposited complexes, and its cost is read from relaxed scans of the free ligand in IEF-PCM (water), with the torsion
held between 90° and 180° and every other coordinate optimised at B3LYP-D3(BJ)/6-311++G(d,p).

Molecular dynamics of the 6RL9 complex (sulfonamidate on the zinc, Glu106 a carboxylate, crystal water kept) used
OpenMM [@eastman2024openmm] with ff14SB [@maier2015ff14sb], GAFF2 (general Amber force field [@wang2004development]) and
a bonded zinc model [@li2016mcpb], in explicit water at 300 K and 1 bar. Nine 100 ns runs, three each from the crystal
pose and the docked near-native and ring-rotated poses (@sec:validation), were classified by ring orientation with
criteria fixed beforehand (Section S16).

## 2.6. ADMET check and software {#sec:software}

Sulfanilamide is an old drug with measured pharmacokinetics, so predicted ADMET (absorption, distribution, metabolism,
excretion and toxicity) properties add little; a check with ADMET-AI 1.4.0 [@swanson2024admet] is in Section S14. Ligand
preparation and the docking analyses used RDKit 2026.03.1 and gemmi 0.7.5, with NumPy and SciPy; the figures were drawn
from the stored results with Matplotlib and PyVista (references in Section S1). All scripts (calculations, analysis and
figures) were written and run with the generative artificial-intelligence tool Claude (Claude Opus 5.5 in Claude Code;
Anthropic) under the authors' direction.

# 3. Results {#sec:results}

## 3.1. Where the proton sits: ion pair, continuum and crystal {#sec:proton}

Optimised in the gas phase from the X-ray structure, the ion pair does not remain an ion pair. H1A moves to O5 (O5–H1A
0.998 Å, *δ* = +0.760 Å; @tbl:proton, [@fig:proton b]), P1–O5 stretches from 1.503 to 1.577 Å, the length of a P–OH
bond, and C1–N1 contracts from 1.463 to 1.416 Å. What remains is a neutral sulfanilamide···H~3~PO~4~ adduct held by
O5–H1A···N1, whose sulfanilamide half carries an NPA charge of only +0.062 e. Complexes of a sulfonic acid with
bipyridines also stay neutral in the gas phase [@fedorov2025influence].

In IEF-PCM (water) the proton stays on N1 ([@fig:proton a]): N1–H1A is 1.097 Å (*δ* = −0.387 Å), and P1–O5 (1.520 Å)
and C1–N1 (1.459 Å) stay close to their crystal values. N1···O5, however, contracts to 2.581 Å (2.760 Å in the
crystal), and the relaxed scan from the other end finds a second minimum, confirmed by its Hessian, with the proton on
O5 (*d*(N1–H1A) = 1.631 Å, *δ* = +0.604 Å; [@fig:proton d], Table S7). The neutral pair lies 0.54 kcal mol⁻¹ below the
salt in electronic energy and further once vibrations are counted, because the salt has the larger zero-point energy
(ZPE): Δ(*E* + ZPE) = −1.66 and Δ*G*~298~ = −1.64 kcal mol⁻¹ (@tbl:proton, Table S6).

The highest point of the scan, 0.72 kcal mol⁻¹ above the salt at *d*(N1–H1A) = 1.25 Å, lies near the closest
donor–acceptor approach (N1···O5 2.495 Å); it is the maximum on a 0.10 Å grid (0.75 kcal mol⁻¹ by spline), and no saddle
point was located. A one-dimensional Schrödinger equation on the relaxed profile puts the ground level within 0.1 kcal
mol⁻¹ of the barrier top, with probability 0.84 on the O5 side (Table S10 and Section S5); ½*hν* of the N1–H1A stretch
alone is 3.13 kcal mol⁻¹. The model moves the proton with the N–H reduced mass (0.940 u), although N1···O5 shortens from
2.581 to 2.495 Å on the way to the barrier and opens to 2.657 Å beyond it; the heavy atoms that move with the proton
raise its effective mass, so the model overstates the delocalisation, and no two-dimensional treatment was attempted.
The conclusion does not need its precision: in the continuum the proton is not confined to the salt well.

<!-- Caption copied from figures/captions.md (group A, fig02_proton, 29 Sep 2026, VM burst 3); keep the two in step. -->

![Position of the proton H1A. (a–c) The X-ray ion pair (thin black) overlaid on (a) the IEF-PCM (water) and (b) the
gas-phase optimum, both B3LYP-D3(BJ)/6-311++G(d,p) and superposed on the cation heavy atoms, and on (c) formula unit 1
of the relaxed periodic PBE-D3(BJ) cell, drawn in the cell frame after removal of the mean heavy-atom shift; heavy atoms
with H1A (disc) and H6. *δ* = *d*(N1–H1A) − *d*(H1A···O5), −1.02 Å in the X-ray structure with H1A as refined. (d)
Relaxed scan of *d*(N1–H1A) in IEF-PCM (water): electronic energy relative to the salt minimum. The open circle at 1.25
Å stopped at the step limit with its forces below the convergence thresholds (Table S7); the two chains meet at 1.35 Å
within 0.014 kcal mol⁻¹. (e) Single points with four functionals and DF-MP2 (6-311++G(d,p), IEF-PCM (water)) at the
B3LYP-D3(BJ) scan geometries, relative to the salt end point. (f) Constrained proton transfer in the periodic cell:
*d*(N1–H1A) of both formula units held, every other atom relaxed in the fixed X-ray cell; energy per formula unit (f.u.)
relative to the salt. The black strokes have the slope given by the constraint force, which is positive at every point,
so no point is a minimum. Grey open markers: an earlier configuration with O5–H1A held at 0.98 Å. (g) Donor–acceptor
distances along the profile: N1···O5 (teal) and the bonds that the crystal adds (brown). Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig02_proton.png){#fig:proton pdf-max-height=12}

<!-- include: tables/table02_proton.md -->

Energies this small are within the error of a single functional, so the three scan points were recomputed with four
more methods at the B3LYP-D3(BJ) geometries ([@fig:proton e]; Table S9). All of them put the neutral end point below
the salt, and by more than B3LYP-D3(BJ) does: ωB97X-D by 0.72, PBE0-D3(BJ) by 1.27, DF-MP2 by 1.42 and M06-2X by 2.47
kcal mol⁻¹. The barrier is less certain. ωB97X-D and DF-MP2 put the 1.25 Å point 0.98 and 0.20 kcal mol⁻¹ above the
salt, while PBE0-D3(BJ) (−0.05 kcal mol⁻¹) and M06-2X (−0.51 kcal mol⁻¹) put it at or below the salt, so that for them
the salt is not a minimum along this path.

The dielectric ladder measures how much screening the salt needs (@fig:dielectric; Table S8). Started from the water
optimum, the pair gives the proton to O5 at *ε* = 4 (*δ* = +0.675 Å, P1–O5 1.578 Å, C1–N1 1.413 Å) and keeps it on N1 at
*ε* = 8 (*δ* = −0.257 Å): a salt minimum first exists between *ε* = 4 and 8. With one start per *ε* this brackets where
the minimum appears, not a thermodynamic crossover. Stronger screening lengthens N1···O5 again (2.526, 2.563 and 2.581 Å
at *ε* = 8, 20 and 78.4), since the anion then needs to pull less on H1A, yet even in the water continuum the salt is
the higher of two minima. Removing D3 leaves *δ* at −0.386 Å, and def2-TZVP raises it to −0.292 Å without transferring
the proton (Table S11).

<!-- Caption copied from figures/captions.md (group A, fig03_dielectric, 29 Sep 2026); keep the two in step. -->

![Dielectric ladder: the ion pair optimised in IEF-PCM with *ε* = 1 (gas phase), 4, 8, 20 and 78.4 (water),
B3LYP-D3(BJ)/6-311++G(d,p). (a) Proton metric *δ*; (b) *d*(P1–O5) and (c) *d*(C1–N1), two heavy-atom bonds that report
the protonation state without relying on the position of H1A. Circles: proton on N1 (salt); diamonds: proton on O5
(neutral pair); lines join optimisations with the proton on the same atom and only guide the eye. The gas-phase
optimisation started from the X-ray structure and the *ε* = 4–20 optimisations from the water optimum with the proton
on N1; at *ε* = 4 the proton moved to O5 (*δ* = +0.675 Å), at *ε* = 8 it stayed on N1 (*δ* = −0.257 Å), so a salt
minimum first exists between *ε* = 4 and 8 (grey band; one start per *ε*). In water both minima exist: the salt (*δ* =
−0.387 Å) and the neutral pair at the end of the relaxed scan (*δ* = +0.604 Å), −0.54 kcal mol⁻¹ relative to the salt.
P1–O5 is 1.503 Å in the crystal, 1.520 Å in the salt minimum in water and 1.578 Å when O5 carries the proton (*ε* =
4); C1–N1 is 1.463, 1.459 and 1.413 Å. Dashed: X-ray structure (its *δ* uses the refined H1A; X-ray N–H bonds are
short, so this *δ* is a lower bound); dotted: periodic PBE-D3(BJ) model (*δ* = −0.610
Å). Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig03_dielectric.png){#fig:dielectric}

Relaxed from the deposited structure, the periodic *Pc* cell keeps H1A on N1 in both formula units (*δ* = −0.610 Å,
N1···O5 2.734 Å; [@fig:proton c]). Starts with H1A on O5 relax back to the same salt, whether the proton was moved in
both formula units or, with no symmetry imposed, in one only (Table S13). The constrained path shows why ([@fig:proton
f]; @tbl:proton, Table S15). With PBE-D3(BJ) the energy rises at every point, to +51.6 kJ mol⁻¹ (12.3 kcal mol⁻¹) per
formula unit at *d*(N1–H1A) = 1.70 Å, and everywhere the constraint force pushes the proton back towards N1 (56–133 kJ
mol⁻¹ Å⁻¹). Exact exchange changes the curve most where the proton is shared: the PBE0-D3(BJ) correction peaks at +4.1
kJ mol⁻¹ at 1.25 Å, raises the four inner points and lowers only the last, to +50.0 kJ mol⁻¹ (12.0 kcal mol⁻¹).
Neither functional gives a neutral minimum. The inner points are the least certain, because N1···O5 contracts to 2.568
Å at 1.40 Å, at the limit of about 2.57 Å below which efficient periodic methods proved unreliable in a benchmark
[@abramov2024it]. A structure with O5 and H1A both frozen and O5–H1A at 0.98 Å lies 10.4 kJ mol⁻¹ above the PBE-D3(BJ)
profile at the same N1–H1A distance, so freezing O5 and compressing the O–H bond overstate the cost.

In the continuum the whole attraction between the ions runs through one bond, and the two proton wells come close
together. The crystal spreads the load: the anilinium group donates to three different anions (N1–H1A···O5,
N1–H1B···O6 and N1–H1C···O6; @tbl:hbonds), O5 also accepts the short O4–H4A···O5 bond of the neighbouring anion, and
N1···O5 stays at 2.734 Å. Moving H1A strains all of these at once ([@fig:proton g]). By 1.70 Å, N1–H1B···O6 and
N1–H1C···O6 have lengthened from 2.836 and 2.829 Å to 2.926 and 2.968 Å, as expected once the amino group has lost its
charge, and O4–H4A···O5 from 2.503 to 2.566 Å, as O5 now carries a proton of its own.

<!-- include: tables/table04_hbonds.md -->

## 3.2. X-ray versus optimised geometry {#sec:geometry}

Of the three models, the periodic one reproduces the 15 heavy-atom bonds of the ion pair best (@fig:parity, @tbl:geom),
with a mean absolute error (MAE) of 0.017 Å, no bond off by more than 0.033 Å and heavy-atom angles within 1.0°. The
IEF-PCM ion pair is close on average (MAE 0.024 Å), but its largest errors, up to 0.083 Å on P1–O3 and P1–O4 and 0.062 Å
on S1–N2, fall on bonds whose groups donate hydrogen bonds in the crystal but find no acceptor in the pair. Part of
these errors is basis-set incompleteness. With def2-TZVP they fall to 0.034 Å on S1–N2 and 0.062 Å on P1–O4, and the
heavy-atom MAE to 0.015 Å (Tables S3 and S11). The heavy atoms alone would have exposed the transfer in the gas phase.

<!-- Caption copied from figures/captions.md (group A, fig04_parity, 29 Sep 2026); keep the two in step. -->

![Optimised against X-ray lengths of the 15 heavy-atom bonds of the ion pair. IEF-PCM (water) and gas phase:
B3LYP-D3(BJ)/6-311++G(d,p); periodic model: PBE-D3(BJ) with the cell fixed at the X-ray values (formula unit 1 of the
relaxed cell). (a) Parity plot with *y* = *x* (line) and ±0.02 Å (grey band); the bond classes of @tbl:geom are labelled
at their X-ray length and the legend gives the mean absolute error (MAE) over the 15 bonds. (b) Signed deviation (model
− X-ray) of every bond, grouped by class. The largest ion-pair errors are on P1–O3 and P1–O4 (up to +0.083 Å in
IEF-PCM), whose O–H groups donate the O3–H3A···O1 and O4–H4A···O5 bonds in the crystal but are free in the ion pair, and
on S1–N2 (+0.062 Å); the periodic model, which keeps all nine hydrogen bonds, has no deviation above 0.033 Å. In the gas
phase H1A sits on O5, which lengthens P1–O5 (+0.075 Å) and shortens C1–N1 (−0.046 Å). Bonds to hydrogen are omitted: the
X-ray C–H are riding and the O–H restrained. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig04_parity.png){#fig:parity}

<!-- include: tables/table03_geom.md -->

Superposed as a whole on the crystal pair, the IEF-PCM pair has a heavy-atom RMSD of 0.54 Å, but the cation and the
anion fitted separately give 0.115 and 0.044 Å ([@fig:proton a]): each ion keeps its crystal shape, and the anion
turns about the N1···O5 axis once no other anion holds it.

## 3.3. The hydrogen-bond network {#sec:hbond}

The periodic model keeps all nine hydrogen bonds of the deposited structure, with donor–acceptor distances within
0.036 Å of the X-ray values (@tbl:hbonds). In graph-set terms the network combines the R²₂(8) ring of the ion pair
([@fig:scheme a]) with a pair of C(4) chains along *a*, one of anions linked by O4–H4A···O5 and one of cations linked
by N2–H2B···O2, which the N1–H···O bonds cross-link. The second leg of the R²₂(8) ring is a long C–H···O contact
(C6···O4 3.288 Å), so the ring rests on its N–H leg.

At the N1–H1A···O5 critical point of the IEF-PCM pair (@fig:qtaim) *ρ* = 0.0783 a.u. and ∇²*ρ* = +0.1558 a.u., but the
total energy density is negative (*H* = −0.0238 a.u.), the mark of partial covalent character [@cremer1984chemical];
C6–H6···O4 is an ordinary closed-shell contact (*ρ* = 0.0100 a.u., *H* = +0.0010 a.u.). In the isolated pair this second
leg depends on dispersion: without D3 the IEF-PCM optimum opens C6···O4 from 3.383 to 4.798 Å; in the crystal the
packing fixes it (3.289 Å in the periodic model). These numbers describe the contracted continuum pair (N1···O5 2.581 Å;
crystal 2.760 Å). The bond-energy estimate from the potential energy density would, for N1–H1A···O5, exceed the
interaction energy of the whole pair (Section S9). We therefore use critical-point descriptors only to compare contacts.

In the five-pair cluster cut from the periodic salt all nine hydrogen bonds have a bond critical point (@tbl:hbonds;
Table S22). The densest is the anion–anion O4–H4A···O5 bond (*ρ* = 0.0713 a.u., *H* =
−0.0129 a.u.), ahead of N1–H1A···O5 (0.0503 a.u.) and the two other N–H···O6 bonds of the anilinium group (0.0373 and
0.0368 a.u.). Placing the proton on O5 lowers *ρ* by 34 and 29 % at N1–H1B···O6 and N1–H1C···O6 and by 24 % at
O4–H4A···O5, whose *H* rises to −0.0023 a.u.; the sulfonamide N–H···O bonds change by 11 % or less, and O3–H3A···O1
gains 13 %. The transfer thus loosens exactly the contacts that tie the NH~3~⁺ group and O5 into the network; the
extrema of the electrostatic potential of the IEF-PCM pair mark the same atoms (Section S9).

<!-- Caption copied from figures/captions.md (group B, fig05_qtaim_nci, 29 Sep 2026); keep the two in step. -->

![Interactions between the two ions of the ion pair in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p). (a) IGMH
isosurface of the inter-fragment *δg* (*δg*^inter^ = 0.005 a.u.) between the cation and H~2~PO~4~⁻, coloured by
sign(*λ*~2~)*ρ* from −0.08 a.u. (blue, attraction) through 0 (green) to +0.02 a.u. (red, steric repulsion). The disc
between H1A and O5 is the N1–H1A···O5 bond: the atom pairs H1A···O5, N1···O5 and H6···O4 contribute 24.4, 11.9 and 5.6
% of the summed inter-fragment atom-pair *δG* index (integrated *δg*^inter^ 0.280 a.u.). The small green disc between
H6 and O4 is C6–H6···O4, the second leg of the R²₂(8) ring. (b) QTAIM molecular graph in the same view: bond critical
points (red dots) and ring critical points (open squares, one in the benzene ring and one in the R²₂(8) ring), with
the two inter-ion bond paths dotted and *ρ* and ∇²*ρ* (a.u.) at their critical points. N1–H1A···O5 (N1···O5 2.581 Å)
has a negative total energy density (*H* = −0.024 a.u.), the signature of a partly covalent, charge-assisted bond;
C6–H6···O4 (C6···O4 3.383 Å) is a weak closed-shell contact (*H* = +0.001 a.u.). Only the 27-atom ion pair is drawn;
@tbl:hbonds also gives *ρ* at all nine hydrogen bonds of the crystal, among them the anion–anion O4–H4A···O5 (O4···O5
2.500 Å), in a cluster cut from the periodic cells. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig05_qtaim_nci.png){#fig:qtaim}

The shortest hydrogen bond of the crystal, O4–H4A···O5 (2.500 Å), joins two anions. The refinement restrained O4–H4A
to 0.82 Å, which puts H4A 0.30 of the way from O4 to O5; the periodic model, free of that restraint, moves it to 0.41
(O4–H4A 1.031 Å, H4A···O5 1.488 Å) and keeps O4···O5 at 2.503 Å. The proton stays on O4, but its DFT position is an
estimate: a GGA can shift protons spuriously along hydrogen bonds [@leblanc2018pervasive], and a classical relaxation
lacks the nuclear quantum effects that spread protons in such bonds [@radekstocek2024hydrogenbond].

## 3.4. Hirshfeld surfaces {#sec:hirshfeld}

Hirshfeld surfaces (@fig:hirshfeld; method in Section S10) rank the contacts of the crystal from the X-ray
structure alone. The deepest red spot belongs to the anion, *d*~norm~ = −0.842 at O4–H4A···O5 (H···O 1.595 Å), and the
deepest on the cation, −0.721, to H1A···O5 (1.756 Å). Next come O1···H3A (−0.712) and the two N1–H···O6 contacts
(−0.647 and −0.644), so the four deepest spots on the cation are the three N–H···O bonds of the anilinium group and
the O3–H3A···O1 bond that its sulfonyl group accepts. H···O contacts make up 51.1 % of the cation surface (189.6 Å²)
and 80.3 % of the anion surface (90.4 Å²), while the cation also has a hydrophobic share, 15.4 % of its surface in
H···C contacts of the ring. A related sulfonamide salt shows the same dominance of H···O contacts
[@moroni2024synthesis].

<!-- Caption copied from figures/captions.md (group B, fig06_hirshfeld, 29 Sep 2026); keep the two in step. -->

![Hirshfeld surfaces of the two ions in the crystal, computed from the X-ray structure with the X–H bonds set to the
neutron lengths (C–H 1.083 Å, N–H 1.009 Å, O–H 0.983 Å) and a promolecule of spherical free atoms. (a) Cation and (b)
H~2~PO~4~⁻ mapped with *d*~norm~ (CrystalExplorer definition with Bondi radii), from −0.85 (red) through 0 (white) to
+2.1 (blue); the red spots that face the viewer are named by their contact. The deepest spot of the anion, *d*~norm~ =
−0.842, belongs to the anion–anion O4–H4A···O5 bond (H···O 1.595 Å with the normalised O–H); that of the cation, −0.721,
to H1A···O5 (1.756 Å). The patchy blue shading is part of the *d*~norm~ definition, which changes its reference radius
where the nearest atom inside or outside the surface changes. (c) Cation and (d) anion: fingerprint plots (*d*~i~ and *d*~e~, the distances from the surface to the nearest nucleus inside and outside it; bins of 0.02 Å, colour from blue to red with increasing surface area per bin), full and for the three largest
contact types, each over the full plot in grey, with its share of the surface. H···O/O···H contacts make up 51.1 % of
the cation surface (189.6 Å²) and 80.3 % of the anion surface (90.4 Å²); the two sharp spikes are the donor and acceptor
halves of the strong O–H···O and N–H···O bonds. Every contact type is shown in the Supplementary
material. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig06_hirshfeld.png){#fig:hirshfeld}

## 3.5. Frontier orbitals and electronic excitations {#sec:fmo}

In the IEF-PCM ion pair the HOMO (−7.64 eV) is shared between the cation (55 %) and the anion, while the LUMO (−1.67 eV)
lies almost entirely on the cation ([@fig:fmo a; @fig:fmo b; @fig:fmo c]). The Kohn–Sham gap of 5.98 eV thus has partial
inter-ion character, and it is neither the fundamental nor the optical gap. The vertical ΔSCF energies (@tbl:fmo) give
the pair in water *I* = 7.49 eV and *A* = 1.86 eV, *η* = 2.82 eV and *ω* = 3.88 eV, close to the orbital-based values
(Table S16) because the continuum stabilises the ions that ΔSCF creates.

<!-- Caption copied from figures/captions.md (group B, fig07_fmo_tddft, 29 Sep 2026); keep the two in step. -->

![Frontier orbitals, density of states and electronic excitations of the ion pair in IEF-PCM (water). (a) HOMO (−7.64
eV) and (b) LUMO (−1.67 eV) at the B3LYP-D3(BJ)/6-311++G(d,p) level, isosurfaces at ±0.03 a.u. (orange positive, violet
negative), in the view of @fig:qtaim. By Hirshfeld orbital composition the HOMO is shared between the cation (55 %) and
H~2~PO~4~⁻, while the LUMO lies on the cation (99 %). (c) Total and fragment-projected density of states (Gaussian
broadening, FWHM 0.3 eV), normalised to the highest peak in the window; Kohn–Sham gap 5.98 eV. (d)
CAM-B3LYP/6-311++G(d,p) absorption spectra in IEF-PCM (water) with non-equilibrium solvation (20 singlet states,
Gaussian FWHM 0.3 eV) for the ion pair, the cation and neutral sulfanilamide, with the ion-pair states as sticks
(oscillator strength, right axis). The first strong band of the pair (S2, 211 nm, *f* = 0.18) and its strongest state
(S3, 185 nm, *f* = 0.87) sit where the cation absorbs (S2, 208 nm). The
neutral sulfanilamide released at pH 7.4 absorbs at 246 nm (S2, *f* = 0.50). (e) Inter-ion charge-transfer number of all
20 ion-pair states with CAM-B3LYP (filled) and B3LYP (open; same basis and solvation), marker area proportional to *f*.
With B3LYP four of the first five states are anion-to-cation charge-transfer states (CT number up to 0.91); with
CAM-B3LYP the largest CT number among the first five is 0.04, and the states with appreciable charge transfer all lie
above S5 and are weak. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig07_fmo_tddft.png){#fig:fmo}

<!-- include: tables/table05_fmo.md -->

Deprotonation makes the drug a better electron donor. Sulfanilamide in water has *I* = 5.98 eV against 7.88 eV for the
cation and *η* = 2.38 against 2.91 eV: the amino lone pair, tied up in an N–H bond in the anilinium ion, is returned to
the ring (Fukui functions in Section S9). Salt formation lowers the static first hyperpolarisability, *β*~tot~, from
1807.7 (sulfanilamide) to 406.6 a.u. (ion pair) (Section S12); neither value is a property of the crystal.

Excited states show why a range-separated hybrid was needed ([@fig:fmo d; @fig:fmo e]). With B3LYP four of the first
five states of the ion pair carry charge from the anion to the cation (CT numbers up to 0.91), the spurious low-lying
charge-transfer states of a global hybrid [@maitra2022double]. With CAM-B3LYP no state among the first five has a CT
number above 0.04. The first strong band of the pair (S2, 211 nm, *f* = 0.18) and its strongest state (S3, 185 nm, *f* =
0.87) sit where the cation alone absorbs (208 and 184 nm): the spectrum of the salt is that of its cation. The reported
lower cut-off wavelength of SAPH, an absorption edge and not a band maximum, is 223 nm [@muthuselvi2017structural].
Neutral sulfanilamide absorbs further to the red, at 246 nm (S2, *f* = 0.50). Aqueous sulfanilamide is determined by its
absorption at 258 nm [@wei2025boron], 0.23 eV below this band.

## 3.6. Docking validation and molecular dynamics {#sec:validation}

In the receptors as prepared, three protocols pass the strict gate (@tbl:redock, @fig:redock): 7,8-dihydropteroate
redocked into Sul2 with Vina and with Vinardo, and sulfanilamide cross-docked into the hCA II structure 3HS4 with
AutoDock4Zn. The conventional 2.0 Å criterion adds the sulfamethoxazole anion in DHPS and the sulfonamidate in 6RL9.

<!-- include: tables/table07_redock.md -->

In Sul2 the product returns to the crystal pose with Vina and Vinardo, every seed inside 1 Å (panel d of @fig:redock),
with its hydrogen bonds to Asp174 and Ser214; docked pABA lies on the pABA half of the product at 0.43 Å. Among the
receptors as prepared, only 7S2K passes the strict gate with its own crystal ligand, although that ligand, which fills
both the pterin and the pABA sites, is twice the size of sulfanilamide.

In 6RL9, hCA II crystallised with sulfanilamide itself, the validation run found the near-native pose in every seed but
ranked it first in only half of them and failed both criteria (2.07 ± 1.45 Å, success fraction 0.5). The competing pose
keeps the sulfonamide nitrogen on the zinc and on Thr199 and turns the aniline ring about the Zn–N/S axis, from the
hydrophobic wall of the site to its hydrophilic side, the two halves into which the site is classically divided
[@bonardi2020sulfonamide]; AutoDock4Zn separates the two orientations by only 0.008–0.039 kcal mol⁻¹.

<!-- Caption copied from figures/captions.md (group B, fig09_redock, 29 Sep 2026; Glu106 sentence 11 Oct 2026); keep the
     two in step. -->

![Validation of the docking protocol. The crystal ligand (grey carbons) and the top pose of the seed whose RMSD is the
median of its job (orange carbons) are drawn in the receptor frame, without superposition, with the pocket residues
within 3.8 Å of either ligand (thin sticks; backbone drawn only where its N or O is within 3.5 Å of a ligand), metal
ions (spheres) and cofactors (tan). RMSD: symmetry-corrected heavy-atom RMSD of the top pose, mean ± SD over the seeds
(exhaustiveness 32, 20 Å box). (a) Human carbonic anhydrase II (hCA II), 6RL9: neutral sulfanilamide docked from its
B3LYP-D3(BJ)/IEF-PCM geometry with AutoDock4Zn (AD4Zn in the panel titles), 0.61 ± 0.04 Å over 20 seeds, all within 1.0 Å; the sulfonamide N sits on
the zinc (Zn–N 2.20 Å, 1.93 Å in the crystal). From an ETKDG/MMFF94s conformer the validation run was bimodal (2.07 ±
1.45 Å) and failed the 1.0 Å gate. (b) hCA II, 3HS4, cross-docking: reference is the 6RL9 sulfanilamide carried into
3HS4 by superposition of 257 Cα atoms (RMSD 0.328 Å). AutoDock4Zn gives 0.62 ± 0.03 Å from the ETKDG conformer (orange;
10 of 10 seeds within 1.0 Å) and 3.22 ± 0.02 Å from the DFT geometry (purple), where the ring takes the rotated
orientation with the sulfonamide still on the zinc. (c) Dihydropteroate synthase (DHPS) of *Yersinia pestis*, 3TZF chain
B: sulfamethoxazole anion with AutoDock Vina, 1.92 ± 0.04 Å, inside the conventional 2.0 Å criterion in every seed but
outside the 1.0 Å gate; pterin pyrophosphate (pterin-PP) and Mg²⁺ are kept. (d) Sul2, 7S2K: 7,8-dihydropteroate with AutoDock Vina,
0.51 ± 0.03 Å, 10 of 10 seeds within 1.0 Å; PPi, pyrophosphate. Both hCA II receptors have Glu106 as the neutral acid; with the carboxylate
the run of panel (a) gives 1.34 ± 1.30 Å, 15 of 20 seeds within 1.0 Å (Section S15). The seed-by-seed RMSDs are in the
Supplementary material. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig09_redock.png){#fig:redock}

The input geometry changes which orientation ranks first. Docked from its DFT geometry instead of an ETKDG/MMFF94s
conformer, sulfanilamide takes the crystal pose in all 20 seeds of 6RL9 (0.61 ± 0.04 Å; panel a of @fig:redock). The
same change moves the 3HS4 cross-dock the other way, from 0.62 ± 0.03 Å to 3.22 ± 0.02 Å, with the ring rotated and the
sulfonamide still on the zinc (panel b). We therefore report 6RL9 as partially validated: the sulfonamide nitrogen
stays on the zinc (eight runs, Section S13: means 2.09–2.37 Å; 1.93 Å in the crystal) and on Thr199; the ring
orientation is reproduced only for a given input and receptor.

With Glu106 as a carboxylate the DFT geometry gives the crystal pose in 15 of 20 seeds. With four crystal waters kept,
no AutoDock4Zn run returns the rotated pose (20-seed runs without them: 40 and 44 poses): the waters occupy its place,
so a rigid receptor excludes it by construction (Section S15). In the simulations (Section S16) the six runs from
crystal and near-native poses kept their orientation, and the three from the rotated pose left it for the crystal one
by an inversion at the sulfonamide nitrogen 140–360 ps after the ligand restraints were lifted, in the equilibration
with the backbone restrained, without return in production. The force-field terms at this nitrogen failed two
quantum-chemical tests, and one error (the rotated placement 3.38 kcal mol⁻¹ too high in a bare zinc-site model)
favours this outcome. Neither docking test nor simulation decides the orientation.

Acetazolamide is reproduced by no scoring function (AutoDock4Zn places its sulfonamide on the zinc but not its tail,
2.70 ± 0.20 Å), so it checks only the scores. AutoDock4Zn scores it better than sulfanilamide by 0.62 kcal mol⁻¹ in 6RL9
and 0.47 kcal mol⁻¹ in 3HS4 (validation runs), the order of the observed dissociation constants (46 nM against 11 µM at
37 °C and pH 7.0) [@linge2023plbd; @linkuviene2018thermodynamic], but the scores also follow the heavy-atom counts (13
against 11), so the order is weak evidence. The same function ranks sulfamethoxazole first in both structures although
none of its top poses reaches the zinc; over zinc-bound top poses only, acetazolamide leads.

In DHPS the sulfamethoxazole anion lands 1.92 ± 0.04 Å from its crystal pose in every seed (panel c of @fig:redock),
inside the conventional criterion and outside the strict one, so the 3TZF poses are read qualitatively. 1AJ0, the *E.
coli* enzyme crystallised with sulfanilamide, failed outright: every seed put the ligand 11.64 ± 0.05 Å away in a
surface cleft. We attribute this to the two loops that build the pABA pocket during catalysis [@yun2012catalysis],
closed in 3TZF and open in 1AJ0, which was therefore left out of production (Section S13).

## 3.7. Poses of the sulfanilamide species and the controls {#sec:poses}

Vina hardly distinguishes the protonation states, having no explicit electrostatic term [@trott2010autodock]: in 3TZF
the neutral molecule, the cation and the sulfonamidate all score −6.4 kcal mol⁻¹ (@tbl:scores), so the species had to be
fixed from the pH (@sec:docking).

<!-- include: tables/table08_scores.md -->

In DHPS (panel a of @fig:poses) sulfanilamide lies where the aniline and sulfonamide groups of the crystal
sulfamethoxazole lie (RMSD over these atoms 0.20 ± 0.03 Å with Vina and 0.96 Å with Vinardo; the redocked
sulfamethoxazole anion overlays its own core at 1.09 Å). The sulfonyl oxygens accept from the backbone NH and the
hydroxyl of Ser222 (2.91 and 2.98 Å), and the aniline nitrogen points at the pterin pyrophosphate (3.16 Å). The score,
−6.4 kcal mol⁻¹, ties for third of eight, behind sulfamethoxazole (−7.8) and its anion (−7.6) and level with pABA
(−6.3 kcal mol⁻¹).

<!-- Caption copied from figures/captions.md (group B, fig10_poses, 29 Sep 2026); keep the two in step. The panel
     reference is written [@fig:redock b] here so that crossref.lua prints "Fig. 8b". -->

![Production poses of neutral sulfanilamide (B3LYP-D3(BJ)/IEF-PCM geometry; orange, thick sticks) and of the control of
each receptor (blue, thin sticks): cluster-1 representatives, i.e. the top pose of the best-scoring seed of the largest
cluster, with the pocket residues within 3.8 Å of either ligand; right, contacts of the sulfanilamide pose. The panel
titles give the mean top-pose scores over the seeds (AD4Zn, AutoDock4Zn). (a) DHPS, 3TZF chain B: the aniline-SO₂N core
of sulfanilamide lies on that of crystal sulfamethoxazole (core RMSD 0.20 Å); pterin-PP, pterin pyrophosphate. (b) Sul2,
7S2K: sulfamethoxazole, not drawn, scores −8.2 kcal mol⁻¹; PPi, pyrophosphate. (c) hCA II, 6RL9, and (d) hCA II, 3HS4:
the sulfonamide N of the poses shown binds the zinc and meets Thr199 OG1; in 3HS4 the ring takes the rotated orientation
of [@fig:redock b]. Contact maps: receptor N or O within 3.2 Å of a ligand N or O (distances in Å; two values when two
atoms of one residue reach the same ligand atom); π contacts found by PLIP (T, T-shaped; P, parallel; π–cation); grey
boxes, the other residues within 4.0 Å, with an arc where PLIP finds a hydrophobic
contact. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig10_poses.png){#fig:poses}

In Sul2 (panel b) sulfanilamide docks into the pABA region of the product, between Ser52, the pyrophosphate and Ser214.
Its score, −6.4 kcal mol⁻¹, ties with that of the cation at the bottom of the nine species and trails pABA by 1.1 kcal
mol⁻¹, whereas in DHPS the two are level. The shift comes from pABA, which scores 1.2 kcal mol⁻¹ better in 7S2K than in
3TZF while sulfanilamide scores the same in both, and sulfamethoxazole, the drug that Sul2 makes bacteria resist, still
outscores pABA in Sul2 (−8.2 kcal mol⁻¹). Against pABA, sulfanilamide loses 1.2, the sulfamethoxazole anion 0.9 and
neutral sulfamethoxazole 0.8 kcal mol⁻¹ from DHPS to Sul2 with Vina (acetazolamide gains 0.2), shifts within the error
of the scoring function. The pockets also differ. Sul2 was docked without the pterin half of its product, whereas the
pterin pyrophosphate caps the pABA pocket of 3TZF.

At the zinc site (panel c) the AutoDock4Zn pose in 6RL9 makes the contacts of the crystal complex: the sulfonamide
nitrogen binds the zinc (2.21 Å) and donates to Thr199 OG1 (2.55 Å), and a sulfonyl oxygen accepts from the backbone NH
of Thr199 (2.98 Å). In 3HS4 (panel d) the representative pose keeps the nitrogen on the zinc (2.32 Å) and on Thr199
(2.46 Å) while the ring turns to the hydrophilic side, but only 6 of the 20 top poses meet the 2.3 Å site criterion.

In the control receptor, gyrase B, Vina scores sulfanilamide −6.0 kcal mol⁻¹ (−5.9 with the crystal water; Tables S33
and S34), within 0.5 kcal mol⁻¹ of its Vina scores in DHPS, Sul2 and hCA II. Only the level of these scores is
interpreted, because the protocol redocks the gyrase inhibitor only with the crystal water kept (0.52 Å; 6.68 Å
without).

## 3.8. From the crystal to the binding site: the sulfonamide torsion {#sec:torsion}

The C3–C4–S1–N2 torsion links the crystal to the binding sites (@fig:torsion). We follow it as *τ*~fold~, the larger
of |C3–C4–S1–N2| and |C5–C4–S1–N2|, so that small values put S1–N2 perpendicular to the ring plane and large values in
it. The crystal has 102.7° and the periodic model keeps it (103.5°), but the isolated models relax towards the
perpendicular: the ion pair to 90.9° in IEF-PCM (water), neutral sulfanilamide to 90.4° and the cation to 89.6°. The
sulfonamidate stops at 101.7°, but its torsion is free (below). The crystal value is therefore set by the packing.

<!-- Caption copied from figures/captions.md (group B, fig11_torsion, 29 Sep 2026); keep the two in step (rewrapped so
     that no line starts with a vertical bar). -->

![The C3–C4–S1–N2 torsion from the crystal to the docked poses, given as *τ*~fold~ =
max(|C3–C4–S1–N2|, |C5–C4–S1–N2|), which does not depend on which ring carbon is called C3: small values put S1–N2
perpendicular to the ring plane, large values in it. (a) Relaxed torsional strain of neutral sulfanilamide (green) and
the sulfonamidate (pink): C3–C4–S1–N2 held, everything else optimised, B3LYP-D3(BJ)/6-311++G(d,p) in IEF-PCM (water),
relative to the unconstrained minimum (open markers); filled markers are the computed points, the line a
shape-preserving interpolation. Sulfanilamide rises steadily to 1.89 kcal mol⁻¹ at 180°, while the sulfonamidate stays
within 0.10 kcal mol⁻¹. Thin dashed lines: MMFF94s profiles of the same species (torsion restrained, relative to the
MMFF94s minimum). Black triangles: the torsions of the ligands in the deposited complexes, sulfamethoxazole in 3TZF
(105.6°) and sulfanilamide in 1AJ0 (111.3°) and 6RL9 (132.1°), which cost sulfanilamide 0.18, 0.28 and 0.84 kcal
mol⁻¹. Dashed: the SAPH crystal (102.7°). (b) *τ*~fold~ of the top pose of every seed of every production job of the
three species, all docked from their DFT geometries, one row per receptor and scoring function; dotted lines repeat
the DFT minima; triangles: the torsion of the crystal ligand of that receptor (open in the 3HS4 rows: the 6RL9
sulfanilamide carried into 3HS4, the cross-docking reference; the 7S2K ligand has no sulfonamide). All 280 poses lie
between 104.4° and 170.5°. With AutoDock4Zn the 6RL9 poses average 147.0°, 15° beyond the crystal ligand, which raises
the relaxed strain of sulfanilamide by 0.45 kcal mol⁻¹. Drawn by a script written and run with the generative AI tool Claude (Claude Opus 5.5 in Claude Code; Anthropic) under the authors' direction; no image-generation tool was used.](figures/fig11_torsion.png){#fig:torsion}

No docked pose returns to the perpendicular minimum, but the proteins do not hold the sulfonamide perpendicular
either: in the deposited complexes *τ*~fold~ is 105.6° for sulfamethoxazole in 3TZF, 111.3° for sulfanilamide in 1AJ0
and 132.1° for sulfanilamide in 6RL9 ([@fig:torsion b]). The Vina poses in DHPS match the crystal value (107.7° in the
representative pose; Vinardo scatters to 129.3°), whereas the AutoDock4Zn poses in 6RL9 overshoot it by 15° (mean
147.0°) and the Sul2 poses turn furthest. Neither
Vina nor AutoDock4Zn has a torsion-angle energy term [@trott2010autodock; @santosmartins2014autodock4zn], so contacts,
mostly with the pocket, set these values.

Relaxed scans of the free ligand ([@fig:torsion a]; Table S30) put a price on the rotation. For neutral sulfanilamide
the energy rises steadily from the minimum at 90.4° to 1.89 kcal mol⁻¹ with S1–N2 in the ring plane. The torsions of
the crystal ligands cost it 0.18 (3TZF), 0.28 (1AJ0) and 0.84 kcal mol⁻¹ (6RL9), and those of the docked poses
0.16–1.81 kcal mol⁻¹ (Table S31), of which the AutoDock4Zn overshoot in 6RL9 adds 0.45 kcal mol⁻¹. The sulfonamidate,
the form that binds the zinc, is flat within 0.10 kcal mol⁻¹ (MMFF94s comparison in Section S13).

Internal strain is rarely part of a docking score [@gu2021ligand]. Here, at up to 1.81 kcal mol⁻¹, it is of the same
order as the score differences in @tbl:scores, which therefore cannot be read as differences in the stability of the
bound states.

# 4. Discussion {#sec:discussion}

## 4.1. What keeps the proton on N1 {#sec:disc-proton}

The three environments separate bulk screening from the crystal as a whole. A dielectric continuum supplies screening
and no hydrogen bonds [@herbert2021dielectric]; the periodic cell supplies both, together with the Madelung field of
an ionic lattice. Screening alone turns the neutral adduct of the gas phase into a metastable salt between *ε* = 4 and
8, yet in the water continuum the salt still lies 0.54–2.47 kcal mol⁻¹ above the neutral pair in electronic energy,
depending on the method, and 1.64 kcal mol⁻¹ in Gibbs energy. In the lattice no neutral minimum was found, and the
constrained transfer costs 50.0 kJ mol⁻¹ (12.0 kcal mol⁻¹) per formula unit at 1.70 Å with PBE0-D3(BJ). The models
differ in more than their environment, yet with PBE0-D3(BJ) single points in both, the signs stay opposite.

Which part of the crystal holds the proton is harder to isolate. The transfer strains the hydrogen-bond network: the two
further N–H···O bonds of the NH~3~⁺ group and the anion–anion bond to O5 lengthen along the periodic path, and in the
cluster the density at their critical points falls by 24–34 %. A finite piece of the network does not keep H1A on N1,
however. The same five ion pairs, cut from the lattice and computed in vacuum with every proton moved, prefer the
O5-protonated arrangement by 54.0 kJ mol⁻¹ per pair, whereas the periodic cell puts that arrangement 57.1 kJ mol⁻¹ per
formula unit above the salt. A finite cluster of ions is a poor model of an ionic solid (its salt form has a HOMO–LUMO
gap of only 1.15 eV, and its outer ions miss some of their hydrogen bonds), so the reversal shows only that what sets
the sign lies outside the cluster; clusters of other acid–base pairs, by contrast, can form the salt spontaneously
[@omalley2021crystallization]. We therefore attribute the salt to the crystal environment as a whole, hydrogen bonds and
ionic lattice together, and cannot say which contributes more.

This is what the solvation reading of the ΔpK~a~ rule predicts for the intermediate range, where the lattice acts as the
solvent and a salt forms only if ancillary hydrogen bonds stabilise the ions [@cruzcabeza2022role]. Even the protonation
site of sulfanilamide changes with its surroundings: the isolated cation carries the proton on the sulfonamide nitrogen,
whereas in solution the aniline nitrogen is protonated [@uhlemann2021preferred]. With a strong acid, as in the chloride
of the same cation, the hydrogen bonds arrange ions that would exist anyway; in SAPH the environment determines whether
the proton is transferred at all. A gas-phase ion pair is accordingly no model of this salt: what is computed for it
describes the neutral adduct.

## 4.2. The released drug and its targets {#sec:disc-drug}

At pH 7.4, far above the pK~a~ of the anilinium group, the dissolved salt loses its proton, so what reaches a target is
sulfanilamide, with phosphate as a buffer anion; a gram of SAPH delivers 0.637 g of the drug. SAPH is a solid form of
sulfanilamide. Its biological activity should be that of sulfanilamide at the delivered dose, and models of its
complexes with a protein must start from the neutral drug or its anion.

Against the deposited complexes the docking recovers the pABA-site and zinc contacts and leaves the ring orientation at
hCA II to the crystal structure. It does not resolve the resistance of Sul2: the shifts of sulfanilamide and
sulfamethoxazole against pABA from DHPS to Sul2 (0.8–1.2 kcal mol⁻¹) have the direction of that resistance
[@venkatesan2023molecular] but lie within the error of the scores.

## 4.3. What the calculations cannot settle {#sec:limitations}

We computed the periodic profile in the fixed X-ray cell; by a rough isotropic estimate from the residual pressures, a
variable-cell relaxation, which was not run, would raise the cost at 1.70 Å by 1.1 kJ mol⁻¹ per formula unit (Table
S15). The hybrid values are single points at PBE-D3(BJ) geometries, the choice of pseudopotentials moves the PBE
energy difference by at most 0.8 kJ mol⁻¹ per formula unit, and the PBE-to-PBE0 shift was not calibrated on the
molecular ion pair. No point-charge embedding separates the hydrogen bonds from the Madelung field. The protons were
treated classically in the crystal and with a one-dimensional model in the continuum; a path-integral treatment
[@radekstocek2024hydrogenbond] would show how far H1A and H4A spread at 294 K.

The continuum is the weakest part of the molecular model. With electronic energies it places separated sulfanilamide and
H~3~PO~4~ 5.82 kcal mol⁻¹ lower than the separated ions, whereas the measured ΔpK~a~ of −0.13 implies −0.17 kcal mol⁻¹,
an error of 5.6 kcal mol⁻¹ that exceeds every energy difference between the two minima. Gas-phase zero-point energies
(Table S2) widen it by 1.3 kcal mol⁻¹; other solvation models and larger-basis energies were not tested, and no explicit
water was added, so IEF-PCM (water) means a dielectric continuum without specific solvation, whose message is
qualitative: screening alone leaves the neutral pair competitive. Hirshfeld atom refinement against the deposited
reflections [@radekstocek2024hydrogenbond; @chalupna2026disorder] or neutron diffraction would place H1A and H4A
experimentally.

We read orders within one receptor and one scoring function; the Sul2–DHPS comparison, the one exception, uses pABA as
the internal reference. Vinardo refits Vina's terms and drops one [@quiroga2016vinardo], so their agreement is weak
corroboration [@masters2020evaluation]. The scores also grow with size (correlation coefficient −0.87 between Vina score and heavy-atom
count), and per-atom measures reverse the bias [@reynolds2007role], so neither is used to argue potency. Nor does a Vina
score identify a target (Section 3.7), as a docking benchmark over the essential proteome of *E. coli* also found
[@wong2022benchmarking]. Each receptor is one rigid crystal conformation and the ligands move only in their torsions.
Molecular dynamics, which commonly follows docking [@guterres2020improving; @shimizu2026molecular;
@gopinath2022molecular], covers hCA II alone, with one force field and three 100 ns runs per start; a further run with
refitted sulfonamide-nitrogen terms, themselves not accepted, visited both orientations, so the persistent crystal
orientation is specific to the force field. Protonation states are assumptions: in the docking His64 was HIE and Glu106 the neutral acid, and
AutoDock4Zn was run with the neutral sulfonamide although the zinc binds the anion [@linkuviene2018thermodynamic], so
the poses reproduce the geometry of the bound state, not its protonation. We know of no affinity or bioassay data for
SAPH (an antibacterial test exists for the reported citrate of its cation [@muthuselvi2022molecular]), and the observed
order of acetazolamide before sulfanilamide at hCA II is the only external anchor. That order comes from the protonation
reactions linked to binding: once they are subtracted, both share an intrinsic dissociation constant of 9.8 nM
[@linge2023plbd].

# 5. Conclusions {#sec:conclusions}

An isolated ion pair of SAPH does not keep its proton on N1. In the gas phase it becomes a neutral
sulfanilamide···H~3~PO~4~ adduct, and in a water continuum the salt survives only as the higher of two minima, 1.64 kcal
mol⁻¹ above the neutral pair in Gibbs energy, with the ground level of the proton near the top of the barrier between
them; four more methods favour the neutral pair, and screening creates the salt minimum between *ε* = 4 and 8 without
making it the lower one. In this continuum, dielectric screening alone does not make SAPH a salt, and properties
computed for the gas-phase pair belong to a species that the crystal does not contain.

The crystal does hold the proton: every start with H1A on O5 returned to the salt, and a constrained transfer raised
the energy at every step, to 50.0 kJ mol⁻¹ per formula unit at 1.70 Å with PBE0-D3(BJ). The transfer loosens the other
two N–H···O bonds of the anilinium group and the anion–anion O4–H4A···O5 bond, but five ion pairs cut from the lattice
do not hold the proton in vacuum. At ΔpK~a~ ≈ 0 it is the crystal environment, its hydrogen bonds and ionic lattice
together, that sets the protonation state, as the solvation picture of the ΔpK~a~ rule proposes for the intermediate
range.

Dissolved at pH 7.4, SAPH releases neutral sulfanilamide, whose sulfonamide torsion the proteins turn away from the
perpendicular at a cost of at most 0.84 kcal mol⁻¹. Docking tested by redocking places it on the zinc and
Thr199 of hCA II and on the pABA site of DHPS, but its ring orientation at hCA II changes with the input geometry and
the protonation of Glu106, and the scores do not separate DHPS from the resistant Sul2 beyond their error. Molecular
dynamics, whose force-field terms failed two quantum-chemical tests, is in line with the crystal pose.

# CRediT authorship contribution statement

**C. Muthuselvi:** Conceptualization, Resources, Investigation, Writing – review & editing. **Anantha Raman
Lakshmipathi:** Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization,
Writing – original draft, Project administration. **D. Saravanakkumar:** Validation, Writing – review & editing.
**Venkatramanan Varadharajan:** Validation, Writing – review & editing.

# Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have
appeared to influence the work reported in this paper.

# Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit
sectors.

# Data availability

The data that support the findings of this study are available within the article and its Supplementary material
(Appendix A), which includes the coordinates of the optimised structures, the periodic cells and the docked poses, the
input files of the calculations and, for the molecular dynamics, the inputs, parameters, time series and final frames of
every run in machine-readable form. The complete inputs and outputs of all calculations, the prepared receptors and
ligands, and the scripts that produce every number, table and figure form a data package with a SHA-256 manifest; the
trajectories of the simulations, which the package lists with their checksums, are too large to be submitted with the
article. The package and the trajectories have not been deposited in a public repository: the files that reproduce the
reported results are in the Supplementary material, and the package and the trajectories are available from the
corresponding author on reasonable request. The starting structures are public: the crystal structure in the
Cambridge Structural Database (CCDC 1016950) and the proteins in the Protein Data Bank under the codes given in Section
2.5 and Section S13.

# Appendix A. Supplementary material

Supplementary material to this article is submitted with it in two files. The document (PDF) holds the computational
settings, software versions and provenance of every calculation, energies and Cartesian coordinates, the complete
comparison with the X-ray geometry, frequencies and thermochemistry, the proton scan, the dielectric ladder and the
functional check, the periodic convergence tests and the constrained proton-transfer profile, the wavefunction and
Hirshfeld analyses, TD-DFT and NLO details, the docking preparation, validation, scores, torsion scans, contacts and
control receptor, the ADMET check, the docking tests of the protonation of Glu106 and of four crystal waters, and the
molecular dynamics with its settings, runs and force-field tests (Sections S1–S16, Tables S1–S40, Figs. S1–S10). The
archive (ZIP) holds the coordinates of the optimised structures, periodic cells and docked poses (XYZ, CIF, SDF), the
PySCF, Quantum ESPRESSO and docking input files, the TD-DFT excitations (CSV), the receptors and comparison files of the
two docking tests, and for the molecular dynamics the inputs, parameters, scripts, time series and last frames of every
run, the analysis tables and the quantum-chemical records; the trajectories are not in it.

# Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the authors used Claude (Anthropic), under their direction, in order to write and
run the scripts of the calculations, the analysis and the figures (Section 2.6), to draft and edit the text, and to
search for, select and check the references, each of which was taken from its Crossref record. Every number in the text
is read by script from the stored results, and those that carry the argument were re-derived by independent code; the
figures and the graphical abstract were drawn by these scripts, the plots and structure views from the computed data and the
deposited crystal structures, and no image-generation tool was used. After using this tool, the authors reviewed and
edited the content as needed and take full responsibility for the content of the published article.

# References

::: {#refs}
:::
