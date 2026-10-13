---
# Cover letter for the Journal of Molecular Graphics and Modelling (Elsevier; target since 12 Oct 2026). Addressed to
# the journal's Editor, Dr. Michelle A. Sahai (one Editor and four Associate Editors in the journal's list of 12 Oct
# 2026, supplied by the corresponding author); under one page. Order: the question, the new insight in the journal's
# terms ("draw conclusions and explore implications"), what the two earlier reports on the salt covered (both by the
# first author) and what is new against them, said positively (the corresponding author's wish of 12 Oct 2026; "first
# DFT study" is NOT claimed, because the 2017 article has B3LYP calculations; every general statement is one that
# the manuscript itself makes), the docking and simulation part with its limits, the reproducibility material, and
# the statement on prior and concurrent submission. Funding, declarations and suggested reviewers are entered in the
# submission form and are not repeated here. Every number comes from paper/main/data/facts.json
# (numbers_audit.py main --manuscript this file). The date is printed from the field `date` below.
# Build: python3 ops/make_submission.py --letters   (writes paper/main/SAPH_cover_letter.docx, no line numbers)
title: "Cover letter"
date: "13 October 2026"
---

Dr. Michelle A. Sahai\
Editor, Journal of Molecular Graphics and Modelling

Dear Dr. Sahai,

We submit the manuscript "Proton transfer in 4-sulfamoylanilinium dihydrogen phosphate: a combined molecular and
periodic DFT study with molecular docking" for consideration as a research article in the Journal of Molecular
Graphics and Modelling.

Protonated sulfanilamide and phosphoric acid have almost equal pK~a~ values (ΔpK~a~ ≈ 0), where the ΔpK~a~ rule cannot
tell a salt from a neutral cocrystal, yet they crystallise as the salt 4-sulfamoylanilinium dihydrogen phosphate. The
new insight is that, for this salt, molecular and periodic models give opposite answers on where the proton sits, and
the manuscript shows where the answer changes. The gas-phase ion pair optimises to a neutral adduct; in a water continuum five
electronic-structure methods place the neutral pair below the salt, and a salt minimum first appears between *ε* = 4
and 8. In the periodic crystal the energy rises at every point of a constrained transfer, to 50.0 kJ mol⁻¹ per
formula unit with PBE0-D3(BJ), while five ion pairs cut from the lattice prefer the neutral form in vacuum: the crystal
as a whole holds the proton.

Two earlier reports on this salt, both by the first author of this manuscript and both cited, gave its crystal
structure and, later, its measured spectra with molecular Hartree–Fock and B3LYP calculations. New here are the
periodic calculations, the comparison of the molecular with the periodic model, and the docking: to our knowledge no
periodic calculation and no docking study of the salt has been reported. The implication concerns how such salts are
modelled: at ΔpK~a~ ≈ 0 the proton cannot be taken for granted in an isolated ion pair, and here the properties
computed for the gas-phase pair belong to a species that the crystal does not contain.

At pH 7.4 the salt releases neutral sulfanilamide. In carbonic anhydrase II, docking tested by redocking keeps the zinc
and Thr199 contacts, but the top-ranked ring orientation depends on the input geometry and on the protonation of
Glu106. In nine molecular dynamics simulations of 100 ns only the crystal pose persists. The force-field terms at the
sulfonamide nitrogen failed two quantum-chemical tests, so the result is reported as in line with the crystal pose
and not as independent evidence for it.

For reproducibility the supplementary material gives, also in machine-readable form, the settings of every
calculation, the coordinates of the optimised structures, the input files and, for the simulations, the force-field
parameter files with the inputs and scripts that rebuild each topology. The manuscript has not been published before
and is not under consideration by another journal, and all authors have approved its submission.

Yours sincerely,

Anantha Raman Lakshmipathi (corresponding author), on behalf of all authors\
Associate Professor, Department of Mechanical Engineering\
Rajalakshmi Institute of Technology, Chennai, Tamil Nadu, India\
E-mail: raman3366@gmail.com
