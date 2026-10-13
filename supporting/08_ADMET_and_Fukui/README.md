# 08 ADMET and Fukui

| Path | What |
|---|---|
| `sulfanilamide/admet/` | ADMET-AI 1.4.0 predictions and RDKit 2026.03.6 descriptors for neutral sulfanilamide, the species present at pH 7.4 (`admet.csv`, `admet.json` with provenance, `run.log`, `tdc_dataverse_ids.json`: the Harvard Dataverse file ids of the Therapeutics Data Commons training sets, downloaded to check which endpoints were trained on sulfanilamide). |
| `cation/admet/` | The same for the 4-sulfamoylanilinium cation (reference only). |
| `sulfanilamide/fukui/`, `cation/fukui/` | Condensed Fukui functions (f⁺, f⁻, f⁰ and the dual descriptor Δf, from NPA and Hirshfeld charges) from single points on the N, N+1 and N−1 electron systems in IEF-PCM (water), B3LYP-D3(BJ)/6-311++G(d,p): `fukui.csv`, `fukui.json`, `fukui.md`, `LOG.md`, `run.log`. |
| `SI_section_S14_ADMET_text.md` | The source text of Section S14 of the supplementary material (pandoc Markdown). |

- The printed tables are 05_Tables/SI/Table_S17_fukui and Table_S35_admet (Markdown and CSV).
- The scripts are `10_Reproducibility/scripts/paper/si/compute/admet.py` and `fukui.py`.
- Not included: the checkpoints of the Fukui single points (`N.chk`, `Nm1.chk`, `Np1.chk`; 10_Reproducibility/excluded_files.csv).
- ADMET here is a sanity check, not a result: sulfanilamide is in the training data of 26 of the 41 ADMET-AI
  endpoints. The human half-life that the supplementary material compares with the prediction (about 10 h) is our
  ln 2/*k* conversion of the urinary-excretion rate constant measured in healthy men (Koizumi 1964; Nogami 1968), as
  Section S14 states.
