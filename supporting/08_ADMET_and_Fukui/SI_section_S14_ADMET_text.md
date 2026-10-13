<!-- SI text block for the ADMET section of paper/si/manuscript.md (heading "ADMET and toxicity sanity check");
     manuscript.md includes tables/tableS_admet.md right after it. Every number is a fact of paper/main/data/facts.json (ADMET
     section, keys admet.*), checked with numbers_audit.py; every citation key is in lit/candidates.jsonl. -->

Sulfanilamide has a long clinical record, so its absorption, binding and toxicity are known from measurement, and a
prediction for it tests the predictor more than it informs about the molecule. We computed descriptors and ADMET
endpoints locally for neutral sulfanilamide, the species present at pH 7.4 (pK~a~ 2.02 at 25 °C for the
anilinium group [@greenbaum1955comparison; @zheng2026digitized] and 10.43 for the sulfonamide [@kim2025pubchem]), and,
for reference only, for the cation. Descriptors come from RDKit 2026.03.6, with the Wildman–Crippen log *P* [@wildman1999prediction] and the
topological polar surface area [@ertl2000fast]; the endpoints come from ADMET-AI 1.4.0 [@swanson2024admet], whose
Chemprop models were trained on the TDC data sets. @tbl:s-admet collects the values.

Sulfanilamide passes the Lipinski [@lipinski2001experimental], Veber [@veber2002molecular] and Egan
[@egan2000prediction] filters and misses one Ghose criterion [@ghose1999knowledge], the atom count: 19 atoms against a
lower limit of 20. Its Crippen log *P*, −0.08, is higher than the measured octanol–water value of −0.62
[@knox2024drugbank], an overestimate that also appears in the ADMET-AI lipophilicity (0.51).

A comparison with experiment says little where the drug itself was used to train the model, and sulfanilamide is in the
training data of 26 of the 41 ADMET-AI endpoints. For solubility, rat oral toxicity, mutagenicity and liver injury the
model returns labels it has seen: log *S* −1.67 against −1.34 in AqSolDB [@sorkun2019aqsoldb] (7.5 g L⁻¹ measured,
log *S* −1.36), −log(LD~50~/mol kg⁻¹) 1.66 against 1.645 [@zhu2009quantitative], an Ames probability of 0.01 for a
compound tested negative [@hansen2009benchmark], and a liver-injury probability of 0.95 for a drug classified DILI
positive [@thakkar2020drug]. Where the model had not seen the drug it can be tested, and it overestimates plasma protein
binding: it predicts 54 %, whereas equilibrium dialysis of sulfanilamide against human plasma at pH 7.4, run at
refrigerator temperature, gives about 22 % bound [@davis1943binding]. The high absorption probabilities, 0.97 for intestinal absorption and 0.73 for oral
bioavailability, agree with the rapid oral uptake of sulfanilamide, which reached 4–8 mg per 100 mL of blood 6 h after
a dose [@taylor1940comparative]; that is a blood level, not a measured bioavailability. The predicted half-life, 16.3 h, is longer than the ≈10 h
derived from urinary excretion after a 1 g oral dose in healthy men, whose overall elimination rate constant of about
0.067 h⁻¹ gives *t*~1/2~ = ln 2/*k* (our conversion; neither source states a half-life)
[@koizumi1964absorption; @nogami1968studies].

None of these numbers bears on potency or efficacy. They show that the predictor agrees with the known profile of the
released drug where that drug was part of its training data, and overestimates its plasma binding where it was not. The
cation rows only show how the charge moves the descriptors (Crippen log *P* −0.79); ADMET-AI was not trained on such
species, and the salt dissociates at pH 7.4.
