| Case | Model | Functional | Basis | Environment | Grid | SCF convergence | Source |
|:----------------|:--------|:-----------|:-------------|:-----------|:-----|:-------|:-------------|
| sfa_pair_pcm_d3 | Ion pair | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water) | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_d3/ |
| sfa_pair_gas_d3 | Ion pair | B3LYP-D3(BJ) | 6-311++G(d,p) | Gas phase | level 4 | 10⁻⁹ Ha | sfa_pair_gas_d3/ |
| sfa_pair_pcm_b3lyp | Ion pair | B3LYP | 6-311++G(d,p) | IEF-PCM (water) | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_b3lyp/ |
| sfa_pair_pcm_d3_tz | Ion pair | B3LYP-D3(BJ) | def2-TZVP | IEF-PCM (water) | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_d3_tz/ |
| sfa_pair_pcm_d3 (eps_ladder) | Ion pair | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM, ε = 4, 8, 20 | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_d3/ |
| sfa_pair_pcm_d3 (scan) | Ion pair | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water), relaxed scan | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_d3/ |
| sfa_pair_pcm_d3 (fsp) | Ion pair | ωB97X-D, M06-2X, PBE0-D3(BJ), DF-MP2 | 6-311++G(d,p) | IEF-PCM (water), single points | level 4 | 10⁻⁹ Ha | sfa_pair_pcm_d3/ |
| sfa_pair_pcm_d3 (TD) | Ion pair | CAM-B3LYP, B3LYP | 6-311++G(d,p) | IEF-PCM (water), non-equilibrium | level 4 | 10⁻⁹ Ha | sfa_pair_pcm_d3/ |
| cation / h2po4 / h3po4 / sulfanilamide _pcm_d3, _gas_d3 | Fragments | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water) and gas | level 4 | 10⁻⁸ Ha (IEF-PCM), 10⁻⁹ Ha (gas) | (fragment)_pcm_d3/, (fragment)_gas_d3/ |
| sulfonamidate_pcm_d3 | Sulfonamidate | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water) | level 4 | 10⁻⁸ Ha | sulfonamidate_pcm_d3/ |
| sulfanilamide_pcm_d3_tz | Sulfanilamide | B3LYP-D3(BJ) | def2-TZVP | IEF-PCM (water) | level 4 | 10⁻⁸ Ha | sulfanilamide_pcm_d3_tz/ |
| urea_gas_d3 | Urea | B3LYP-D3(BJ) | 6-311++G(d,p) | Gas phase | level 4 | 10⁻⁹ Ha | urea_gas_d3/ |
| periodic_qe | Periodic cell | PBE-D3(BJ) | SSSP 1.3.0, 70/560 Ry | Crystal, 2×1×1 k-points | – | 10⁻⁸ Ry | periodic_qe/ |
| periodic_qe (hybrid) | Periodic cell | PBE0-D3(BJ), PBE-D3(BJ) | SG15 ONCV, 100/400 Ry (also 80/320) | Crystal, single points at the PBE-D3(BJ) geometries | – | 10⁻⁶ Ry | periodic_qe/ |
| periodic_qe (profile) | Periodic cell | PBE-D3(BJ); PBE0-D3(BJ) single points | SSSP 1.3.0, 70/560 Ry; SG15 ONCV, 100/400 Ry | Crystal, N1–H1A held (damped dynamics) | – | 10⁻⁹ Ry | periodic_QE/ profile.json |
| sfa_pair_pcm_d3_oform | Ion pair, neutral minimum | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water), Hessian | level 4 | 10⁻⁸ Ha | sfa_pair_pcm_d3_oform/ |
| sfa_xclu_salt, sfa_xclu_neutral | Five-pair cluster (M8) | B3LYP-D3(BJ) | 6-31+G(d,p) | Gas phase, single points at periodic geometries | level 3 | 10⁻⁸ Ha | sfa_xclu_salt/, sfa_xclu_neutral/, crystal_cluster_qtaim/ |
| sulfanilamide / sulfonamidate _pcm_d3 (torsion_scan) | Sulfanilamide, sulfonamidate | B3LYP-D3(BJ) | 6-311++G(d,p) | IEF-PCM (water), relaxed torsion scan | level 4 | 10⁻⁸ Ha | (species)_pcm_d3/ torsion_scan/ |

Table: Settings of every calculation, with the case name and the folder of the data package (06_DFT_results/) that holds its results. {#tbl:s-settings}

PySCF 2.14.0 with density fitting and geomeTRIC 1.1.1 for the molecular models; IEF-PCM with *ε* = 78.3553 and Lebedev order 29; Quantum ESPRESSO 7.5 for the periodic model; Multiwfn 2026.9.20 for the wavefunction analyses. The IEF-PCM ion pair and the dispersion-free pair were accepted on negligible forces (rms 7.35 × 10⁻⁵ and 5.92 × 10⁻⁵ Ha Bohr⁻¹) after the job was stopped; the IEF-PCM pair is a confirmed minimum (0 imaginary frequencies). DF-MP2 is frozen-core with the aug-cc-pVTZ-RI fitting basis; its IEF-PCM term is the reaction field of the RHF density. The periodic PBE0-D3(BJ) runs use exact exchange with the adaptively compressed operator and a 1×1×1 (also 2×1×1) *q* mesh. Final single points for orbitals, Hessians, NLO and interaction energies used 10⁻¹⁰ Ha. The constrained periodic profile converged on the constrained forces below 10⁻³ Ry Bohr⁻¹ (damped dynamics; the unconstrained relaxations used BFGS and 5 × 10⁻⁴ Ry Bohr⁻¹).
