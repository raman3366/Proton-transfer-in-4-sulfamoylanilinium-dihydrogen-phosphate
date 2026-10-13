"""Final wording of the notes of facts.md / facts.json (ACS integration review, 30 Sep 2026), executed by build_facts.py
after every section has been added (exec in its namespace: FACTS, SECTIONS, DISAGREE, CAVEATS, MISSING).

The dossier ships in the deliverables package, so its notes read as final notes on the values: past-tense statements
of what was computed, found or corrected, with no instructions to a writer, no review-round history and no open items.
No value, unit, precision or source file changes here; only the wording of the notes, the section titles and the
statuses. The last block checks that none of the old process wording is left.
"""
import re as _re


def _sub_all(lst, pairs):
    """Apply (old, new) substring replacements to every string of lst; each old must occur somewhere."""
    for old, new in pairs:
        hit = False
        for i_, s_ in enumerate(lst):
            if old in s_:
                lst[i_] = s_.replace(old, new)
                hit = True
        if not hit:
            raise SystemExit(f'facts_final_notes.py: text not found: {old[:80]!r}')


# ---------------------------------------------------------------------------------------- notes on the values
_sub_all(DISAGREE, [
    ('SETTLED (burst 2, 28 Sep 2026): the old KEY STORY "constrained neutral +3.9 kJ/mol per formula unit" (relax_iifix, '
     'H1A alone held in space) is WITHDRAWN. Its optimiser',
     'The first constrained-neutral periodic structure (relax_iifix, burst 2, H1A alone held in space, +3.9 kJ/mol per '
     'formula unit) is withdrawn. Its optimiser'),
    ('(their numbers sit in withdrawn_value, where numbers_audit cannot match them): never cite them. The constrained '
     'neutral that holds is relax_iirfix',
     '(their numbers sit in withdrawn_value, where numbers_audit cannot match them) and are not cited. The constrained '
     'configuration with the proton on O5 is relax_iirfix'),
    ('iifix may be mentioned only as a methodological note (why O5 must be held too).',
     'iifix appears only in a methodological note: O5 has to be held as well.'),
    ('molecules/sfa/results/sfa_pair_pcm_d3/RESULTS.md (draft by dft/report.py) says "integration grid level 3" and '
     '"Multiwfn ." with no version; the runs used',
     'The automatic draft report of the pipeline (molecules/sfa/results/sfa_pair_pcm_d3/RESULTS.md, written by '
     'dft/report.py; not part of the paper) gives "integration grid level 3" and no Multiwfn version; the runs used'),
    ('and Multiwfn 2026.9.20. Use grid level 4.', 'and Multiwfn 2026.9.20, as the paper and SI Table S1 state.'),
    ('RESULTS.md (draft) reports "ZPE 5.299 eV" — correct (0.19472 Ha = 122.19 kcal/mol) but not a quantity for the '
     'paper; its "IP/EA/χ/η/ω" table is Koopmans-type,',
     'The same draft gives the zero-point energy of the pair as 5.299 eV (0.19472 Ha = 122.19 kcal/mol), which the paper '
     'does not use, and a Koopmans-type "IP/EA/χ/η/ω" table,'),
    ('η 2.816 vs 2.988 eV). Use ΔSCF in the paper, Koopmans only as the comparison.',
     'η 2.816 vs 2.988 eV). The paper gives the ΔSCF values, with the Koopmans-type values as the comparison.'),
    ('KEY STORY "barrier (~0.6–0.7 kcal/mol)": the B3LYP-D3 electronic barrier',
     'Barrier: the B3LYP-D3 electronic barrier'),
    ('Quote 0.72 (or "about 0.7"), not 0.6.', 'The paper gives 0.72 kcal/mol.'),
    ('KEY STORY "lower with every other functional": true for the end points with all five methods',
     'Functional check: every method puts the neutral end point below the salt end point'),
    ('KEY STORY "dielectric ladder ε 1 O-side, 8/20/78 N-side": correct, and ε = 4 is also O-side (δ +0.675 Å);',
     'Dielectric ladder: the proton ends on O5 at ε = 1 and 4 (δ +0.675 Å at ε = 4) and on N1 at ε = 8, 20 and 78.4;'),
    ('RESULTS.md (draft) calls ADCH/CM5/Hirshfeld "fragment charge transfer" only; ADCH atomic charges',
     'ADCH, CM5 and Hirshfeld charges are used for the fragment charge transfer; the ADCH atomic charges'),
    ('Use 3.22 everywhere (facts dock.prod.3HS4.sfa.ad4zn.rmsd, Fig. 9b, captions, Table 7). RESOLVED: dock/RESULTS.md '
     'and dock/VALIDATION.md print 3.22 ± 0.02 Å for this RMSD (checked in the numbers audit, round 2, 29 Sep 2026).',
     'The paper (Figure 8b, Table 6), the SI, dock/RESULTS.md and dock/VALIDATION.md print 3.22 ± 0.02 Å (fact '
     'dock.prod.3HS4.sfa.ad4zn.rmsd).'),
    ('SETTLED 29 Sep 2026, rounding ties (rule:', 'Rounding ties (rule:'),
    ('Methods 2.2 says "each SCF was converged to 10⁻⁹ Ha" (protocol.conv_tol) and Table S1 printed 10⁻⁹ Ha for every '
     'molecular run until 29 Sep 2026. As run:',
     'SCF convergence as run:'),
    ('Table S1 and the Methods (Section 2.2) now give the thresholds as run.',
     'Table S1 and Section 2.2 give these thresholds.'),
    ('SETTLED (VM burst 3, 29 Sep 2026): referee M2 was right that the single constrained configuration overstated the '
     'cost.',
     'The single constrained configuration (relax_iirfix) overstates the periodic cost (VM burst 3).'),
    ('The sign does not change:', 'The sign is the same:'),
])

_sub_all(CAVEATS, [
    (' The task statement "ε 1 O-side, 8/20/78 N-side" is right but omits ε = 4, which is also O-side.', ''),
    ('The barrier top (A_02, d = 1.25 Å) is flagged unconverged after 61 steps, but its',
     'The barrier top (A_02, d = 1.25 Å) stopped at the step limit (61 steps) without the convergence flag; its'),
    ('No experimental UV-Vis for SAPH is in hand.', 'Section 3.5 of the main text quotes two measured values beside the computed bands: the wavelength at which aqueous sulfanilamide is determined (258 nm) and the reported UV–Vis cut-off of SAPH (223 nm); no spectrum is compared in detail.'),
    ('Quote our d_norm (CrystalExplorer definition) and say so. CrystalExplorer energy frameworks were not computed.',
     'The paper gives our d_norm (CrystalExplorer definition) and says so.'),
    ('use NPA/CM5 atom values and ADCH only as fragment sums.',
     'the paper uses NPA/CM5 atom values and ADCH only as fragment sums.'),
    ('Revision (section 31): the co-crystal torsions', 'Section 31: the co-crystal torsions'),
    ('No cell-relaxed end point was computed (user decision, 29 Sep 2026): the residual pressures',
     'No cell-relaxed end point was computed: the residual pressures'),
    ('SUPERSEDED (section 36):', 'Superseded (section 36):'),
])

_sub_all(MISSING, [
    ('Stated limitation (user decision, 29 Sep 2026; VM budget closed at 30.0 of 40 VM-h).',
     'Stated as a limitation in Section 4.3 of the main text.'),
    ('No molecular calibration of the PBE-to-PBE0 shift (PBE-D3(BJ) single points on the three IEF-PCM scan geometries, '
     'referee M3.2): stated limitation (user decision, 29 Sep 2026).',
     'No molecular calibration of the PBE-to-PBE0 shift (PBE-D3(BJ) single points on the three IEF-PCM scan geometries): '
     'stated as a limitation in Section 4.3.'),
    ('Referee M5 not run: no PCM Hessians of sulfanilamide and H3PO4, no SMD or radius-scaling single points, no '
     'larger-basis single points of the scan geometries; the continuum error (pt.sep.continuum_error) is stated, not '
     'corrected (user decision).',
     'Continuum: no PCM Hessians of sulfanilamide and H3PO4, no SMD or radius-scaling single points and no larger-basis '
     'single points of the scan geometries were computed; the continuum error (pt.sep.continuum_error) is stated in '
     'Section 4.3, not corrected.'),
    ('Referee M1 controls not run: no point-charge embedding of the ion pair and no cluster-continuum pair;',
     'No point-charge embedding of the ion pair and no cluster-continuum pair were computed;'),
])

# ------------------------------------------------------------------------------------------- section titles
_TITLES = {
    '27. Docking cluster representatives of Fig. 10 (sulfanilamide)':
        '27. Docking cluster representatives of the pose figure (Figure 9; sulfanilamide)',
    '31. Revision after the referee reports (29 Sep 2026)':
        '31. Ligand torsions in the deposited complexes, scores relative to p-aminobenzoate and reference single points',
    '32. Periodic N1–H1A profile (referee M2/M3) and the single-transfer control (VM burst 3)':
        '32. Periodic N1–H1A profile and the single-transfer control (VM burst 3)',
    '33. Neutral IEF-PCM minimum: Hessian, thermochemistry and IR (referee M4, VM burst 3)':
        '33. Neutral IEF-PCM minimum: Hessian, thermochemistry and IR (VM burst 3)',
    '34. QTAIM of the nine crystal hydrogen bonds in a cluster cut from the periodic cells (referee M6, VM burst 3)':
        '34. QTAIM of the nine crystal hydrogen bonds in a cluster cut from the periodic cells (VM burst 3)',
    '35. Relaxed C3–C4–S1–N2 torsion scans (referee M9, VM burst 3)': '35. Relaxed C3–C4–S1–N2 torsion scans (VM burst 3)',
    '37. Closing items: measured hCA II binding of the controls, one "near-native" definition, raw-pose torsions':
        '37. Measured hCA II binding of the controls, one "near-native" definition, raw-pose torsions',
}
_INTRO = [
    ('Numbers asked for by the referees and derived from stored results', 'Numbers derived from stored results'),
    ('placement of the N–H hydrogens, referee M9)', 'placement of the N–H hydrogens)'),
    ('Added on 29 Sep 2026 (night) to close the last audit items; no new calculation.', 'No new calculation.'),
    ('Verified 29 Sep 2026: Crossref record and abstract through lit.py (not retracted; the abstract names 6-311+G(d,p) '
     'and a B3-LYP basis-set series up to 6-311++G(3df,3pd)); the full text is closed (OpenAlex oa_status closed, '
     'publisher 403), so the table row was read from two independent tabulations that agree (see freq.scale_factor).',
     'The Crossref record was checked through lit.py (not retracted), and the table row was read in the full text '
     '(Tables 2 and 6) and in two independent tabulations that agree (see freq.scale_factor).'),
    ('(rows kept in paper/audit/research/plbd_hca2_azm_sfm.json)',
     '(rows kept in paper/audit/research/plbd_hca2_azm_sfm.json, 10_Reproducibility/inputs/ of the data package; source '
     'notes in lit/SOURCE_NOTES.md)'),
]
_seen = set()
for i_, (t_, intro_, keys_) in enumerate(SECTIONS):
    if t_ in _TITLES:
        _seen.add(t_)
        t_ = _TITLES[t_]
    for old, new in _INTRO:
        if intro_ and old in intro_:
            intro_ = intro_.replace(old, new)
            _seen.add(old)
    if intro_:
        intro_ = _re.sub(r'\bNever cited?\b', 'Not cited', intro_)
        intro_ = _re.sub(r'\bnever cited?\b', 'not cited', intro_)
    SECTIONS[i_] = (t_, intro_, keys_)
_missing = [k for k in list(_TITLES) + [o for o, _ in _INTRO] if k not in _seen]
if _missing:
    raise SystemExit(f'facts_final_notes.py: section titles or intros not found: {_missing}')

# ------------------------------------------------------------------------------ single facts and fact wording
_m8 = FACTS['model.M8']
_m8['value'] = 'label reused'
_m8['status'] = None
_m8['meaning'] = ('M8 names the five-pair cluster of section 34 (Table 1); the (pair)2 cluster of the first plan was not '
                  'built')

# The 6RL9 orientation gap of dock/VALIDATION.md used the lowest-RMSD mode; the paper uses the d14 definition
# (near-native <= 1.0 A, best mode of each orientation): dock.val.6RL9.ad4zn.orientation_sep, section 37.
_g = FACTS['dock.val.ad4zn_gap_6RL9']
_g['superseded_value'] = _g['value']
_g['value'] = None
_g['superseded'] = True
_g['status'] = ('SUPERSEDED (not cited): measured on the lowest-RMSD mode (dock/VALIDATION.md section 3, first '
                'definition); replaced by dock.val.6RL9.ad4zn.orientation_sep (best near-native mode against best '
                'ring-rotated mode, section 37)')

_FACT_TEXT = [
    ('; author decision 29 Sep 2026', ''),
    ('confirmed in the full text (29 Sep 2026; author decision: keep 0.9688): ', 'Read in the full text: '),
    ('Notes: paper/audit/research/irscale.md', 'Source notes: lit/SOURCE_NOTES.md'),
    ('(referee estimate)', '(order-of-magnitude assumption)'),
    ('Replaces the 0.000–0.046 of dock/VALIDATION.md, which used the lowest-RMSD mode (numbers audit round 2, m5).',
     'Replaces the 0.000–0.046 first measured on the lowest-RMSD mode (fact dock.val.ad4zn_gap_6RL9, superseded); '
     'dock/VALIDATION.md section 3 gives this value.'),
    (' (referee M2/M3)', ''),
]
for _k, _f in FACTS.items():
    for _fld in ('source', 'meaning', 'note', 'status', 'value'):   # 'value': descriptive strings such as a method
        _v = _f.get(_fld)
        if isinstance(_v, str):
            for old, new in _FACT_TEXT:
                _v = _v.replace(old, new)
            _v = _re.sub(r'\bNever cited?\b', 'Not cited', _v)
            _f[_fld] = _re.sub(r'\bnever cited?\b', 'not cited', _v)

# ------------------------------------------------------------------------------------------------- the check
_PROCESS = _re.compile(r'KEY STORY|\breferees?\b|user decision|author decision|\bQuote\b|\bUse (?:grid|ΔSCF|3\.22)|'
                       r'in hand|rerun pending|UNVERIFIED|never cite|SETTLED|RESOLVED|read first|must travel|full text is closed|'
                       r'not flagged converged', _re.I)
_bad = [s_ for s_ in DISAGREE + CAVEATS + MISSING if _PROCESS.search(s_)]
_bad += [f'section {t_}' for t_, i_, _ in SECTIONS if _PROCESS.search(t_ + ' ' + (i_ or ''))]
_bad += [f'{k_}.{fl_}' for k_, f_ in FACTS.items() for fl_ in ('source', 'meaning', 'note', 'status', 'value')
         if isinstance(f_.get(fl_), str) and _PROCESS.search(f_[fl_])]
if _bad:
    raise SystemExit('facts_final_notes.py: process wording left:\n  ' + '\n  '.join(b_[:160] for b_ in _bad))
