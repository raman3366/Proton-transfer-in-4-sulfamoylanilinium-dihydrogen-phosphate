# Writing guide for the SAPH manuscript and its Supplementary material

Applies to every sentence in `paper/main/manuscript.md` and `paper/si/manuscript.md`. Adapted from the guide
of the four cqc-dft papers (`/root/workspace/cqc-dft/papers/_shared/STYLE.md`); the rules on voice, banned
phrases, dashes, notation and citations are unchanged, the compound-specific parts are rewritten for
4-sulfamoylanilinium dihydrogen phosphate (SAPH) and extended to docking. The paper will be run through
similarity checkers against the colleague's R4 template and
the four sister papers, so it must read as its own piece of work, written by computational chemists who looked
hard at their own numbers.

## Voice

- Write like an experienced physical chemist writing for specialists in the same field. Plain, specific,
  confident where the data are clear, openly uncertain where they are not.
- First person plural is fine ("we optimised", "we attribute"). Past tense for what was done; present tense
  for what the results show and for established knowledge.
- Every paragraph should carry at least one concrete number, comparison or named atom/interaction. A
  paragraph that could be pasted into a paper about a different compound is a paragraph to rewrite.
- Vary sentence length. Mix short declarative sentences with longer ones that carry a qualification.
  Vary paragraph length too. Do not open consecutive sentences or paragraphs the same way.
- Interpret, don't narrate. "N1···O5 shortens by *x* Å in the continuum, so the field of the solvent, not the
  pair itself, holds H1A on N1" beats "Fig. 5 shows the optimised structure." (The example is a pattern, not
  a result.)
- Say what is surprising, what disagrees with experiment, and what the model cannot capture. Specialist readers trust
  papers that state their limitations (an ion pair rather than the lattice, a continuum rather than explicit
  water, PBE-D3 in the periodic model, harmonic frequencies, rigid receptors and empirical scoring in docking).

## Words and constructions to avoid

These are the most recognisable markers of generated text. Do not use them at all:

delve, delving, intricate, intricacies, pivotal, crucial, crucially, underscore(s), showcase(s),
leverage (as a verb), robust (unless statistical), comprehensive, comprehensively, meticulous(ly),
notably, noteworthy, furthermore, moreover, additionally, in addition, in conclusion, overall,
it is worth noting, it is important to note, it should be noted, plays a vital/crucial/key/pivotal role,
sheds light, shed light, paves the way, holistic, tapestry, realm, landscape, testament, seamless(ly),
cutting-edge, state-of-the-art, game-changer, unveil(s), unravel(s), harness(es), elucidate (max once per
paper), facilitate(s), myriad, plethora, a wide range of, a variety of, various, significant(ly) (unless a
statistical test was done), remarkable, remarkably, intriguing, fascinating, novel (let the reader judge),
"In recent years", "In this study, we" (at most once), "This study aims to", "This work provides valuable
insights", "valuable insights", "deeper understanding", "rich", "vibrant", "nuanced", "multifaceted",
"serves as", "stands as", "not only ... but also", "a key", "key role", "key insights", "in the context of",
"in terms of" (max twice), "with respect to" (max twice), "Taken together", "Collectively", "Notably".

Also avoid:
- Em dashes (—) as a sentence device. Use commas, parentheses, semicolons or a new sentence. En dashes in
  ranges (5–10 Å) and bond notation (N–H) are required and are fine.
- Triplets by reflex ("stability, reactivity and selectivity"). Lists of three are fine when there really
  are three things.
- Rhetorical questions, exclamation marks, bold or italic for emphasis in running text.
- Bullet lists in the body. Everything is prose except tables.
- Closing a paragraph with a generic summary sentence that restates it.
- Symmetric "on one hand / on the other hand" framing.
- Stacked hedges ("may potentially suggest"). One hedge per claim, and only when needed.
- Docking clichés: "binding energy" or "binding affinity" for a docking score, "strong/excellent binding",
  "potent inhibitor", "drug candidate", "promising lead". A Vina or AutoDock4Zn number is a score in
  kcal mol⁻¹; call it that.

## Structure

Target journal: the Journal of Molecular Graphics and Modelling (JMGM, Elsevier; user's decision of 12 Oct 2026 after
the Journal of the Indian Chemical Society, the target since 5 Oct 2026, declined the manuscript; subscription route;
the structure below was set for the earlier journal and holds for this one, see
`paper/audit/research/journal_jmgm_2026-10-12.md`). Title (specific, names the compound and the main finding; fixed by the user); under it the author
block (`author-block` in the YAML: names with lower-case affiliation letters, full postal addresses with country, the
corresponding author with e-mail, then every author's e-mail and ORCID iD on a line of its own); Abstract (single
paragraph, 180–250 words, with the main numbers; the journal's limit is 250 words); at most six keywords (`keywords` in
the YAML; the journal takes 1 to 7; no "and" or "of"); 1. Introduction; 2. Computational methods (numbered subsections);
3. Results; 4. Discussion; 5. Conclusions (prose, no bullets, 250–400 words, no new material). The guide's wording is
that of Elsevier's unnumbered-sections template, but the journal prints numbered headings, so the numbers stay and a
section is cited by number ("Section 2.5"). Back matter: CRediT authorship contribution statement; Declaration of
competing interest (Elsevier's sentence, the same as in the separate declaration file); Funding; Data availability;
Appendix A. Supplementary material (what each supplementary file holds); Declaration of generative AI and AI-assisted
technologies in the manuscript preparation process (Elsevier's two sentences, the last section before the list);
References. An Acknowledgements section, if one is needed, goes between Data availability and Appendix A. The use of
the AI tool for the scripts is also stated in the Methods (Section 2.6) with the model and the developer, as Elsevier's
policy asks. Highlights (3 to 5, at most 85 characters each, no abbreviations where a word will do; they state
results), the graphical abstract (500:200, separate TIFF), the cover letter and the declaration of interests (the
document of Elsevier's declarations tool) are separate files, and so is every figure of the main text (`Figure_N`;
line drawings at 1000 dpi or more) (`ops/make_submission.py`). Every figure caption of the main text ends with the
sentence that names the AI tool with which its script was written (guide of 12 Oct 2026); captions explain every symbol
and abbreviation that the figure shows. Citations print as numbers in square brackets (`build.py --journal jics`, the default) and
the reference list in the Elsevier numbered style with titles and ISO 4 journal abbreviations; cite after the word the
reference supports. The supplementary file is called "Supplementary material", never "Supporting Information".
Main text at most 8,150 words of running prose (the count of `lint_prose.py main`; 8,000 until 8 Oct 2026, raised
by about 150 words when the simulations and the two docking tests entered the paper). Everything else (coordinates, full tables, ELF/LOL, Fukui functions,
descriptor tables, ADMET) goes to the SI, whose sections are numbered S1, S2, ... and whose figures and tables
are Fig. S1 and Table S1. Labels follow Elsevier: "Fig. 1" in captions and text, "Eq. (1)" for equations
(crossref.lua prints them). A table title is one phrase or sentence; details go in the notes under the table (build_tables.py).
Abbreviations are defined at their first use in the text (Table 1 notes cover those it uses first).

The Introduction moves from why this salt matters (a sulfa drug protonated by phosphoric acid at
ΔpK~a~ ≈ 0, where the pK~a~ rule makes no prediction) to what is known about SAPH, to the gap, to what this
paper does. It must cite the crystal structure [@muthuselvi2014crystal] (flagged must-cite in `lit/candidates.jsonl`) and
say plainly what the present work adds: proton position across environments, the solid-state H-bond network,
periodic DFT and validated docking of the species that exists at pH 7.4.

## The compound, its models and its species

- Name: 4-sulfamoylanilinium dihydrogen phosphate, abbreviated SAPH after first use. The cation is
  4-sulfamoylanilinium (protonated sulfanilamide), the anion dihydrogen phosphate. It is a salt; do not call the
  crystal a co-crystal, and never call any of these compounds Schiff bases.
- Model names, defined once in Computational methods and then used unchanged: "ion pair" (the 27-atom
  cation–anion pair joined by N1–H1A···O5), "gas phase" and "IEF-PCM (water)" for the two environments,
  "IEF-PCM (ε = 4)" etc. for the dielectric ladder, "periodic model" for the plane-wave PBE-D3 cell.
- Level of theory, written exactly: B3LYP-D3(BJ)/6-311++G(d,p); the sensitivity check is "B3LYP without
  dispersion"; TD-DFT is CAM-B3LYP/6-311++G(d,p) in IEF-PCM (water).
- Docking species: neutral sulfanilamide (primary), its sulfonamidate anion (zinc site), the cation
  (sensitivity case only). Never write that "SAPH" or "the salt" was docked, and never dock or describe the ion
  pair as a ligand. Controls: sulfamethoxazole, pABA, acetazolamide.
- Receptors by PDB code and chain, with organism on first mention: DHPS from *Yersinia pestis* (3TZF, chain B),
  DHPS from *Escherichia coli* (1AJ0), human carbonic anhydrase II (hCA II; 6RL9, 3HS4), Sul2 (7S2K).
- The control receptor (4 Oct 2026) is the ATP site of DNA gyrase B of *Escherichia coli* (4DUH, chain A). It is "the
  control receptor", never a target. Write that it is "no target of sulfanilamide", not that sulfanilamide "does not
  bind" it, which no measurement shows. Nothing docked into it is called an inhibitor
  of gyrase except its crystal ligand ("the crystal inhibitor"). Its two forms are "the receptor as prepared" (the
  protocol of the five receptors) and "the receptor with the crystal water kept". It is read for the level of the
  scores only, and its scores are compared with "the four production structures" (3TZF, 7S2K, 6RL9, 3HS4) through the
  same scoring function, Vina. Never write "its targets" for these four: Sul2 (7S2K) rejects sulfonamides.
  Ciprofloxacin is a comparator that binds gyrase elsewhere (SI only, stated as weak evidence), the dihydrogen
  phosphate anion is not docked, and clorobiocin is not docked (its 1KZN pose does not fit 4DUH).

## Notation and units

- Units with a thin space: 1.084 Å, 5.70 eV, 10.6 kcal mol⁻¹, 3450 cm⁻¹, 16.4 D, 298.15 K.
  Write kcal mol⁻¹ (not kcal/mol) and cm⁻¹ with superscript minus.
- Hydrogen bonds: N1–H1A···O5, written with three middle dots (···, U+00B7 ×3), used consistently.
  Charge-assisted: N⁺–H···O⁻ only where the charges are real in that model. The proton metric is
  δ = *d*(N1–H1A) − *d*(H1A···O5), defined once as an equation.
- Graph sets as R²₂(8) (Bernstein notation; not cited since the reference cut of 30 Sep 2026).
- Symbols in italics via Markdown: *E*~HOMO~, *E*~LUMO~, Δ*E*~gap~, *η*, *S*, *μ*, *χ*, *ω*, *ρ*, *λ*~2~,
  sign(*λ*~2~)*ρ*, *β*~tot~. Equations in LaTeX math ($...$ or display $$...$$ with an {#eq:label}).
- Atom labels exactly as in the CIF (C1–C6, N1, N2, O1, O2, S1 in the cation; P1, O3–O6, H3A, H4A in the anion).
  Never show the `#1_554` symmetry suffixes used internally; mention symmetry-generated fragments in words
  and give the symmetry operation (the anion at *x*, *y*, *z* − 1) once in Computational methods.
- Residues as Arg63, His94 (three-letter code + number, chain on first mention if needed); zinc as Zn²⁺ only
  when its charge matters, otherwise "the zinc ion".
- Numbers: quote to the precision the method supports. Bond lengths 3 decimals (Å), angles 1 decimal,
  energies in eV to 2 decimals in text (3 in tables), frequencies to the nearest cm⁻¹, docking scores to
  1 decimal as mean ± SD over seeds, RMSD to 2 decimals. Never convert a docking score into a K~i~ or K~d~.
- Every number in the text must be traceable to `paper/data/facts.json` (checked by `numbers_audit.py`) or to
  a cited source in the literature notes.

## Citations

- Pandoc keys: `[@lu2012multiwfn]`, `[@becke1993density; @lee1988development]`. Keys come only from
  `lit/candidates.jsonl` (see `python3 paper/_shared/tools/lit.py list`). Never type a reference by hand;
  cite only keys that `lit.py` created, and `build.py` stops on an unknown key.
- Cite for a specific statement the source actually supports (check its abstract; claims marked "(title only)"
  in the pool need the full text first). No citation dumps: at most three references in one bracket, and only
  when they all support the same point.
- Distribute citations through the paper. Introduction and discussion both need them: compare our numbers
  with published values for related salts wherever a fair comparison exists.
- At least 50 distinct references in the final paper, at least half from 2020 onwards. Foundational method
  papers (B3LYP, D3(BJ), basis sets, PCM, QTAIM, NCI, ELF, LOL, conceptual DFT, software) are cited in their
  original form.
- Reference budget (set on 30 Sep 2026 at the co-authors' request): the main text cites about 40 DFT and
  structural-chemistry references beside its docking and biology references (74 in all at the cut; 75 since the
  control receptor of 5 Oct 2026 added the reference of its crystal structure). For the Journal of the Indian Chemical
  Society, which sets no limit, the review before submission restored three classic sources (the salt–cocrystal
  continuum, the ΔpK~a~ rule, proton transfer in gas-phase adducts and clusters) and added recent work and two
  articles of the journal where a sentence uses them; a reference of the journal is cited for what it shows, never
  for the journal's sake. A setting,
  numerical detail or utility program is cited in the SI section that describes it, and the main text points
  there ("(Section S1)"). Record: `paper/audit/research/refs_trim_2026-09-30.md`. Software that has no Crossref DOI (RDKit, PyMOL, Meeko, PDBFixer) is named with its version in
  the text (Section 2.6 and Section S1), not given an invented reference.
- No reference that is retracted or carries an expression of concern (`lit.py` refuses them; `lit.py verify`
  repeats the check before submission).

## Figures and tables

- Refer to every figure and table in the text, in order, as Fig. 1, Figs. 2a and 2b, Table 3.
- Captions are self-contained: what is shown, model (gas phase / IEF-PCM (water) / periodic), level of
  theory, isovalue or colour scale, and what the colours mean. Docking captions name the receptor (PDB code,
  chain), the ligand species, the program and force field, and whether the pose is the top-scored or the
  cluster representative.
- Tables in pipe-table Markdown with a caption line `Table: ... {#tbl:label}`.

## Originality

- Write every sentence fresh for this compound. Do not reuse sentences from the R4 template, from the four
  cqc-dft manuscripts, from planning documents or from any source abstract
  (`overlap.py` checks all of them).
- The Computational methods section is the highest-risk place for overlap with the sister papers: describe the
  protocol in different words and order, and emphasise what was specific to SAPH (the anion translated by one
  cell along *c*, the restrained O–H hydrogens, the proton scan and dielectric ladder, the periodic cell, the
  docking species chosen for pH 7.4).
