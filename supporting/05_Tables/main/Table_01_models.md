| Model | Composition | Charge | Atoms | Basis functions | Use |
|:--------------|:----------------------|------:|-----:|---------------:|:-------------------------|
| Ion pair (M1) | C₆H₉N₂O₂S⁺ · H₂PO₄⁻ (anion at *x*, *y*, *z* − 1) | 0 | 27 | 445 | Primary model: proton position, geometry, spectra, wavefunction analyses |
| Ion pair, def2-TZVP | as M1 | 0 | 27 | 574 | Basis-set check (IEF-PCM) |
| Cation (M2) | C₆H₉N₂O₂S⁺ | +1 | 20 | 313 | Fragment: interaction energy, TD-DFT, descriptors; docking sensitivity case |
| Dihydrogen phosphate (M3) | H₂PO₄⁻ | −1 | 7 | 132 | Fragment: interaction energy |
| Sulfanilamide (M4) | C₆H₈N₂O₂S | 0 | 19 | 306 | Neutral species at pH 7.4; primary docking ligand; proton-transfer energetics |
| Sulfanilamide, def2-TZVP | as M4 | 0 | 19 | 395 | Basis-set check (IEF-PCM) |
| Sulfonamidate (M5) | C₆H₇N₂O₂S⁻ | −1 | 18 | 299 | Docking species for the zinc site; ligand of the simulations |
| Phosphoric acid (M6) | H₃PO₄ | 0 | 8 | 139 | Proton-transfer energetics |
| Urea (M7) | CH₄N₂O | 0 | 8 | 116 | NLO reference (gas phase) |
| Five-pair cluster (M8) | 5 × (C₆H₉N₂O₂S⁺ · H₂PO₄⁻), cut from the periodic cells | 0 | 135 | 1755 (6-31+G(d,p)) | QTAIM of the nine crystal hydrogen bonds, salt and proton on O5 (gas phase) |
| Periodic model (M9) | 2 × (C₆H₉N₂O₂S⁺ · H₂PO₄⁻), *Pc*, *Z* = 2 | 0 | 54 | plane waves (70/560 Ry) | Proton position in the crystal: relaxations and the constrained N1–H1A profile (PBE-D3(BJ), PBE0-D3(BJ)) |

Table: Models, their charges and sizes, and what each is used for. {#tbl:models}

Basis functions are spherical 6-311++G(d,p) unless another basis is stated; all molecular models are closed-shell singlets cut from the crystal coordinates (urea from a force-field conformer). IEF-PCM, integral-equation-formalism polarisable continuum model; TD-DFT, time-dependent density functional theory; QTAIM, quantum theory of atoms in molecules; NLO, nonlinear optical. The ion pair uses 443 of its 445 functions after removal of two near-linear dependencies. The cluster (M8) holds 5 ion pairs that contain all nine hydrogen bonds of the crystal, eight of them at the central pair (N1–H1B···O6 joins two outer pairs), cut from relaxed periodic cell (i) and from the cell with H1A held on O5 at 0.98 Å; it is used for bond descriptors, and its energy only for its sign.
