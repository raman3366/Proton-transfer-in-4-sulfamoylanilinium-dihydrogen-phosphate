| Property | Sulfanilamide | Cation (reference) | Experimental, sulfanilamide | Source |
|:------------------------------|-------------:|------------------:|:------------------------------|:------------------------------|
| *Drug-likeness (RDKit)* |  |  |  |  |
| Molecular weight (g mol⁻¹) | 172.21 | 173.22 | – |  |
| cLogP (Wildman–Crippen) | −0.08 | −0.79 | log *P* −0.62 | [@knox2024drugbank] |
| Molar refractivity (cm³ mol⁻¹) | 42.2 | 40.4 | – |  |
| TPSA, N and O (Å²) | 86.2 | 87.8 | – |  |
| H-bond donors / acceptors (NH + OH / N + O) | 4 / 4 | 5 / 4 | – |  |
| Rotatable bonds | 1 | 1 | – |  |
| Atoms, including H | 19 | 20 | – |  |
| Violations: Lipinski / Veber / Ghose / Egan | 0 / 0 / 1 / 0 | 0 / 0 / 1 / 0 | – |  |
| *Absorption and distribution (ADMET-AI)* |  |  |  |  |
| Aqueous solubility, log *S* (mol L⁻¹) | −1.67† | −1.21 | −1.36 (7.5 g L⁻¹ at 25 °C) | [@knox2024drugbank; @sorkun2019aqsoldb] |
| Lipophilicity, log *D*~7.4~ | 0.51 | −0.04 | log *P* −0.62 | [@knox2024drugbank] |
| Human intestinal absorption (probability) | 0.97 | 0.73 | 70–100 % of an oral dose absorbed (sulfonamides) | [@kim2025pubchem] |
| Oral bioavailability (probability) | 0.73 | 0.50 | blood level 4–8 mg per 100 mL 6 h after an oral dose (no bioavailability measured) | [@taylor1940comparative] |
| Plasma protein binding (%) | 54 | 41 | about 22 (human plasma, dialysis) | [@davis1943binding] |
| Blood–brain barrier penetration (probability) | 0.74 | 0.42 | – |  |
| *Excretion and metabolism (ADMET-AI)* |  |  |  |  |
| Half-life (h) | 16.3 | 10.6 | ≈10 | [@koizumi1964absorption; @nogami1968studies] |
| CYP2C9 / CYP3A4 inhibition (probability) | 0.01† / 0.00† | 0.00 / 0.00 | – |  |
| *Toxicity (ADMET-AI)* |  |  |  |  |
| hERG blocking (probability) | 0.02 | 0.01 | – |  |
| Ames mutagenicity (probability) | 0.01† | 0.01 | negative | [@hansen2009benchmark] |
| Drug-induced liver injury (probability) | 0.95† | 0.86 | positive (DILIst) | [@thakkar2020drug] |
| Rat oral LD~50~, −log(LD~50~/mol kg⁻¹) | 1.66† | 1.94 | 1.645 (3.9 g kg⁻¹) | [@zhu2009quantitative] |
| *Ionisation (experimental only)* |  |  |  |  |
| pK~a~, sulfonamide N–H | – | – | 10.43; 10.6 (20 °C) | [@kim2025pubchem; @knox2024drugbank] |
| pK~a~, anilinium | – | – | 2.02 (25 °C); 2.06 (24 °C) | [@greenbaum1955comparison; @klotz1945isoelectric; @zheng2026digitized] |

Table: Drug-likeness descriptors and ADMET-AI predictions for neutral sulfanilamide, the species present at pH 7.4, with the 4-sulfamoylanilinium cation for reference, against experimental values for sulfanilamide. {#tbl:s-admet}

Descriptors from RDKit (cLogP and molar refractivity [@wildman1999prediction], TPSA over N and O [@ertl2000fast]); filters of Lipinski [@lipinski2001experimental], Veber [@veber2002molecular], Ghose [@ghose1999knowledge] and Egan [@egan2000prediction]; ADMET-AI [@swanson2024admet] probabilities are classifier outputs, not measured frequencies. † Sulfanilamide is in the training data of this endpoint, so the prediction is not independent of the experimental value. The cation lies outside the training domain of ADMET-AI. Sulfanilamide fails only the Ghose atom count (19 atoms, range 20–70); the cation fails only the Ghose logP range. The oral-absorption range is a class statement for sulfonamides; the blood level was measured 6 h after an oral dose. Plasma protein binding: bound/free ratio 4.6 % per g/100 mL of protein at 10 mg/100 mL drug, evaluated at 6 g/100 mL protein. Half-life: ln 2/*k* from the overall elimination rate constant of sulfanilamide in man (about 0.067 h⁻¹, urinary excretion after a 1 g oral dose); the conversion is ours, neither source prints a half-life. Anilinium p*K*~a~: UV absorbance in dilute acid, as evaluated in the IUPAC compilation.
