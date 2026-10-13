# Software versions

From Section S1 of the supplementary material and the run records.

| Task | Software |
|---|---|
| Molecular DFT (optimisation, frequencies, FMO, TD-DFT, ΔSCF, NLO, interaction energies, functional check, Fukui, relaxed torsion scans, five-pair cluster single points) | PySCF 2.14.0 (density fitting, IEF-PCM with Lebedev order 29), geomeTRIC 1.1.1, D3(BJ) dispersion |
| Wavefunction analysis (QTAIM, IGMH, NCI, ELF/LOL, MEP, charges, DOS) | Multiwfn 2026.9.20 |
| Periodic DFT | Quantum ESPRESSO 7.5 (pw.x): SSSP efficiency 1.3.0 potentials for PBE-D3(BJ) (BFGS relaxations; the constrained N1–H1A profile with distance constraints and damped dynamics); SG15 ONCV v1.2 for the PBE0-D3(BJ) single points |
| Docking | AutoDock Vina 1.2.7 (Vina and Vinardo), AutoDock4Zn on AutoGrid 4.2.7 maps (maps made on the x86 VM) |
| Receptor and ligand preparation | PDBFixer (OpenMM 8.6.1), PROPKA 3.5.1, RDKit 2026.03.1, Meeko 0.8.0 |
| Pose analysis | spyrmsd 0.9.0, PoseBusters 0.6.5, PLIP, ProLIF 2.2.2, gemmi 0.7.5 |
| ADMET | ADMET-AI 1.4.0 (Chemprop 1.6.1, PyTorch 2.5.0), RDKit 2026.03.6 |
| Paper build | pandoc 3.11, Python 3.12 (matplotlib 3.11.2, pyvista 0.49.0, RDKit 2026.03.6, python-docx, lxml), LuaLaTeX for the PDF previews |

Hardware: DFT and the periodic runs on an IBM Cloud VM (bx3d-48x240, 24 cores used, x86); docking preparation and
runs, analysis, figures and the paper build on a Raspberry Pi 5 (4 cores, aarch64).

Level of theory: B3LYP-D3(BJ)/6-311++G(d,p) is the primary level (Table S1 of the supplementary material lists
every calculation). B3LYP is the definition PySCF shares with Gaussian (VWN-RPA local correlation). TD-DFT: CAM-B3LYP (and B3LYP) in IEF-PCM (water),
non-equilibrium. Periodic energies: PBE-D3(BJ)/SSSP energy difference plus a PBE0-D3(BJ) − PBE-D3(BJ) correction
with SG15 potentials, at every point of the constrained N1–H1A profile. Five-pair cluster: B3LYP-D3(BJ)/6-31+G(d,p),
gas phase, at the periodic geometries.
