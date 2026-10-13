"""Facts added with the SI completion of 29 Sep 2026 (sections 24-27 of facts.md), executed by build_facts.py
(exec in its namespace, after section 23, so add(), section(), FACTS, BULK_PREFIX, J, r, RES, DOCK, DER, CAVEATS,
MISSING and DISAGREE are the ones defined there).  Kept in its own file so the additions can be read and reviewed
as one block:

  24. Condensed Fukui functions of neutral sulfanilamide and the cation in IEF-PCM (water), from
      molecules/sfa/results/<case>/fukui/fukui.json (paper/si/compute/fukui.py, run on the Pi).
  25. ADMET / drug-likeness sanity check (SI only): RDKit descriptors, rule filters and ADMET-AI predictions from
      molecules/sfa/results/<case>/admet/admet.json (paper/si/compute/admet.py), and the experimental values of
      sulfanilamide they are compared with, each with the source it was read from.
  26. Figure rendering parameters (isovalues, colour ranges, cut-offs): figure_params*.json merged as fig.* keys.
  27. Docking cluster representatives drawn in Fig. 10: every distance, re-derived by derive/d10_rep_contacts.py.
"""


def addbulk(key, value, unit, prec, source, meaning, prefix, what, case=None, check=None, tol=None):
    """A bulk fact (listed in facts.json, summarised in facts.md section 21 by its key prefix)."""
    add(key, value, unit, prec, source, meaning, case=case, check=check, tol=tol)
    FACTS[key]['bulk'] = True
    BULK_PREFIX.setdefault(prefix, [0, source.split(':')[0].split('[')[0], what])
    BULK_PREFIX[prefix][0] += 1


# =====================================================================================================================
section('24. Condensed Fukui functions (IEF-PCM water; SI section S9)',
        'Finite differences of atomic charges of the N, N−1 and N+1 electron systems at the IEF-PCM optimum of each '
        'species (vertical; N: RKS, N±1: UKS doublets; equilibrium IEF-PCM for every state), B3LYP-D3(BJ)/6-311++G(d,p), '
        'grid level 4, run on the Pi with 2 threads (paper/si/compute/fukui.py). f+ = q(N) − q(N+1), f− = q(N−1) − q(N), '
        'f0 = [q(N−1) − q(N+1)]/2 (Yang & Mortier 1986), with NPA charges (PySCF natural atomic orbitals of each state, as '
        'dft/mwfn_extra.py) and Hirshfeld charges (promolecule of spherically averaged neutral B3LYP/6-311++G(d,p) atoms). '
        'Checks: the three total energies equal those of dscf.json (VM run, same code path; its UKS stability analysis found '
        'no instability), and the N-state NPA charges equal mwfn_extra/charges.json.')
FK_SPECIES = [('sfa', 'sulfanilamide_pcm_d3', 'sulfanilamide'), ('cation', 'cation_pcm_d3', 'cation')]
FK_GROUPS = {
    'sfa': {'amino': ['N1', 'H1B', 'H1C'], 'ring': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'H2', 'H3', 'H5', 'H6'],
            'sulfonamide': ['S1', 'O1', 'O2', 'N2', 'H2A', 'H2B']},
    'cation': {'amino': ['N1', 'H1A', 'H1B', 'H1C'], 'ring': ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'H2', 'H3', 'H5', 'H6'],
               'sulfonamide': ['S1', 'O1', 'O2', 'N2', 'H2A', 'H2B']}}
FK_NAMES = {'f_plus': ('fp', 'f+ (nucleophilic attack: where an added electron goes)'),
            'f_minus': ('fm', 'f− (electrophilic attack: where an electron is removed from)'),
            'f_zero': ('f0', 'f0 (radical attack)'), 'dual_descriptor': ('df', 'dual descriptor f+ − f−')}
FK_SCHEME = {'NPA': 'npa', 'Hirshfeld': 'hirsh'}
for sp, case, lab in FK_SPECIES:
    p_ = RES / case / 'fukui' / 'fukui.json'
    if not p_.exists():
        MISSING.append(f'Condensed Fukui functions of the {lab} (IEF-PCM): {r(p_)} not written yet (paper/si/compute/fukui.py).')
        continue
    fk = J(p_)
    src = r(p_)
    atoms = fk['atoms']
    st = fk['states']
    dj = J(RES / case / 'dscf.json')
    add(f'fukui.{sp}.pcm.E_N', st['N']['E_Ha'], 'Ha', 8, src + ':states.N.E_Ha', f'{lab}, IEF-PCM: N-electron energy of the Fukui run (Pi)',
        case=case, check=dj['E_N_Ha'], tol=1e-6)
    add(f'fukui.{sp}.pcm.E_Nm1', st['N-1']['E_Ha'], 'Ha', 8, src + ':states["N-1"].E_Ha', f'{lab}, IEF-PCM: N−1 (UKS doublet, vertical) energy',
        case=case, check=dj['N_minus_1']['E_Ha'], tol=1e-6)
    add(f'fukui.{sp}.pcm.E_Np1', st['N+1']['E_Ha'], 'Ha', 8, src + ':states["N+1"].E_Ha', f'{lab}, IEF-PCM: N+1 (UKS doublet, vertical) energy',
        case=case, check=dj['N_plus_1']['E_Ha'], tol=1e-6)
    add(f'fukui.{sp}.pcm.S2_Nm1', st['N-1']['S2'], None, 3, src + ':states["N-1"].S2', f'{lab}: <S²> of the N−1 doublet (ideal 0.75)', case=case)
    add(f'fukui.{sp}.pcm.S2_Np1', st['N+1']['S2'], None, 3, src + ':states["N+1"].S2', f'{lab}: <S²> of the N+1 doublet (ideal 0.75)', case=case)
    add(f'fukui.{sp}.pcm.I', fk['I_eV'], 'eV', 3, src + ':I_eV', f'{lab}, IEF-PCM: vertical I from the Fukui run (= dscf.{sp}.pcm.I)', case=case,
        check=FACTS[f'dscf.{sp}.pcm.I']['value'], tol=1e-4)
    add(f'fukui.{sp}.pcm.A', fk['A_eV'], 'eV', 3, src + ':A_eV', f'{lab}, IEF-PCM: vertical A from the Fukui run (= dscf.{sp}.pcm.A)', case=case,
        check=FACTS[f'dscf.{sp}.pcm.A']['value'], tol=1e-4)
    ck = fk.get('check_vs_mwfn_extra_N_charges') or {}
    add(f'fukui.{sp}.check.npa_vs_mwfn', ck.get('max_abs_diff_NPA'), 'e', 5, src + ':check_vs_mwfn_extra_N_charges.max_abs_diff_NPA',
        f'{lab}: max |ΔNPA| of the N state against mwfn_extra/charges.json (same recipe, VM run)', case=case)
    add(f'fukui.{sp}.check.hirsh_vs_mwfn', ck.get('max_abs_diff_Hirshfeld'), 'e', 3, src + ':check_vs_mwfn_extra_N_charges.max_abs_diff_Hirshfeld',
        f'{lab}: max |ΔHirshfeld| of the N state against Multiwfn (Multiwfn uses its own built-in free-atom densities, '
        'here B3LYP/6-311++G(d,p) spherical atoms: the promolecules differ, the Fukui differences much less)', case=case)
    add(f'fukui.{sp}.hirsh.grid_points', fk['hirshfeld_grid']['points'], None, 0, src + ':hirshfeld_grid.points', f'{lab}: integration points of the Hirshfeld partition (SCF grid, level 4)', case=case)
    for scheme, s in FK_SCHEME.items():
        for kind, (kk, what) in FK_NAMES.items():
            vals = fk['fukui'][scheme][kind]
            for a_, v_ in zip(atoms, vals):
                addbulk(f'fukui.{sp}.{s}.{kk}.{a_}', v_, None, 3, f'{src}:fukui.{scheme}.{kind}[{a_}]', f'{lab}, {scheme}: condensed {what} on {a_}',
                        f'fukui.{sp}.{s}', f'condensed Fukui functions of the {lab} per atom, {scheme} (Table S9b)', case=case)
            if kind == 'dual_descriptor':
                continue
            add(f'fukui.{sp}.{s}.{kk}.sum', float(sum(vals)), None, 3, f'{src}:sums.{scheme}.{kind}', f'{lab}, {scheme}: Σ {kk} over atoms (must be 1)', case=case)
            i_ = max(range(len(vals)), key=lambda i: vals[i])
            add(f'fukui.{sp}.{s}.{kk}.max', vals[i_], None, 3, f'{src}:fukui.{scheme}.{kind} (largest)', f'{lab}, {scheme}: largest condensed {kk}, on {atoms[i_]}', case=case)
            add(f'fukui.{sp}.{s}.{kk}.max_atom', atoms[i_], None, None, f'{src}:fukui.{scheme}.{kind} (argmax)', f'{lab}, {scheme}: atom with the largest {kk}', case=case)
            mn = min(vals)
            add(f'fukui.{sp}.{s}.{kk}.min', mn, None, 3, f'{src}:fukui.{scheme}.{kind} (smallest)', f'{lab}, {scheme}: smallest condensed {kk}, on {atoms[vals.index(mn)]}'
                + (' (negative condensed Fukui values are a known feature of finite-difference populations)' if mn < 0 else ''), case=case)
            for g, members in FK_GROUPS[sp].items():
                gv = float(sum(vals[atoms.index(m)] for m in members))
                add(f'fukui.{sp}.{s}.{kk}.grp.{g}', gv, None, 2, f'{src}:fukui.{scheme}.{kind} summed over {"+".join(members)}',
                    f'{lab}, {scheme}: {kk} summed over the {g} group ({", ".join(members)})', case=case)
CAVEATS.append('Condensed Fukui functions (section 24) are vertical finite differences in IEF-PCM (water) with each charge state '
               'solvated at equilibrium, at the neutral (or cation) optimum; they describe the isolated, solvated species, not the '
               'crystal. NPA and Hirshfeld give the same ordering of the leading sites but different magnitudes (Hirshfeld '
               'flattens them); a few condensed values are slightly negative (NPA f− on C3, C5 and S1 of sulfanilamide), a known '
               'artefact of condensing finite differences, not a physical negative reactivity.')

# =====================================================================================================================
_ad0 = RES / 'sulfanilamide_pcm_d3' / 'admet' / 'admet.json'
_ad0 = J(_ad0)['admet_ai'] if _ad0.exists() else {}
N_EP = sum(1 for e in _ad0.values() if e['category'] != 'Physicochemical')
N_TR = sum(1 for e in _ad0.values() if (e.get('training_set') or {}).get('sulfanilamide_in_set'))
section('25. ADMET and drug-likeness sanity check (SI only; section S14)',
        'Computed locally (no web server): RDKit 2026.03 descriptors and rule filters, and ADMET-AI 1.4.0 (Chemprop-RDKit '
        'ensembles on the TDC ADMET datasets) in /root/cqc-sfa-tools/admet (aarch64; torch 2.5.0 CPU wheel), by '
        'paper/si/compute/admet.py. The SMILES of each species was checked against the bond graph of its DFT optimum. '
        'ADMET-AI output for sulfanilamide reproduces the prediction shipped in its DrugBank reference set (DB00259) to 2e-6. '
        f'Sulfanilamide itself is in the training data of {N_TR} of the {N_EP} ADMET-AI endpoints (AqSolDB, LD50, AMES, DILI, ClinTox, '
        'CYP inhibition and substrate sets, Tox21, skin reaction); for those a prediction is a recall of its label, not an independent test. The '
        'experimental values are read from the source named in each row; the cation rows are for reference only (it does '
        'not exist at pH 7.4 and is outside the training domain of ADMET-AI). No potency or efficacy claim is attached.')
ADM_SP = [('sfa', 'sulfanilamide_pcm_d3', 'sulfanilamide'), ('cation', 'cation_pcm_d3', 'cation')]
RD_KEYS = [('MW', 'MW', 'g mol-1', 2, 'molecular weight (average masses)'), ('exact_mass', 'exact_mass', 'Da', 3, 'monoisotopic mass'),
           ('cLogP', 'cLogP_Crippen', None, 2, 'Wildman–Crippen logP (RDKit)'), ('MR', 'MR_Crippen', None, 1, 'Wildman–Crippen molar refractivity (cm³ mol⁻¹)'),
           ('TPSA', 'TPSA_NO', 'A2', 2, 'Ertl TPSA over N and O (original definition; used for Veber and Egan)'),
           ('TPSA_SP', 'TPSA_NOSP', 'A2', 2, 'Ertl TPSA including S and P contributions'),
           ('HBD', 'HBD_Lipinski_NHOH', None, 0, 'H-bond donors, Lipinski count (NH + OH)'), ('HBA', 'HBA_Lipinski_NO', None, 0, 'H-bond acceptors, Lipinski count (N + O)'),
           ('HBD_rdkit', 'HBD_RDKit', None, 0, 'H-bond donors, RDKit pharmacophore definition (as ADMET-AI prints)'),
           ('HBA_rdkit', 'HBA_RDKit', None, 0, 'H-bond acceptors, RDKit pharmacophore definition (as ADMET-AI prints)'),
           ('rotb', 'rotatable_bonds', None, 0, 'rotatable bonds (RDKit strict definition)'), ('heavy_atoms', 'heavy_atoms', None, 0, 'heavy atoms'),
           ('atoms', 'atoms_incl_H', None, 0, 'atoms including hydrogen (Ghose count)'), ('QED', 'QED', None, 2, 'quantitative estimate of drug-likeness (RDKit)'),
           ('charge', 'formal_charge', None, 0, 'formal charge')]
AD_PRINT = ['HIA_Hou', 'Bioavailability_Ma', 'Solubility_AqSolDB', 'Lipophilicity_AstraZeneca', 'Caco2_Wang', 'BBB_Martins', 'PPBR_AZ',
            'VDss_Lombardo', 'Half_Life_Obach', 'Clearance_Hepatocyte_AZ', 'CYP2C9_Veith', 'CYP3A4_Veith', 'hERG', 'AMES', 'DILI',
            'ClinTox', 'Carcinogens_Lagunin', 'LD50_Zhu', 'Skin_Reaction']
AD_UNITS = {'%': '%', 'hr': 'h', 'L/kg': 'L kg-1', 'log(mol/L)': 'log(mol L-1)', 'log-ratio': None, 'log(10^-6 cm/s)': 'log(1e-6 cm s-1)',
            'uL/min/10^6 cells': 'uL min-1 (1e6 cells)-1', 'log(1/(mol/kg))': 'log(kg mol-1)', '-': None}
for sp, case, lab in ADM_SP:
    p_ = RES / case / 'admet' / 'admet.json'
    if not p_.exists():
        MISSING.append(f'ADMET for the {lab}: {r(p_)} missing (paper/si/compute/admet.py).')
        continue
    ad = J(p_)
    src = r(p_)
    add(f'admet.{sp}.identity_ok', ad['identity']['identity_ok'], None, None, src + ':identity.identity_ok',
        f'{lab}: SMILES {ad["identity"]["smiles_canonical"]} has the bond graph of the DFT optimum (InChIKey {ad["identity"]["inchikey"]})', case=case)
    for k_, jk, u_, pr, what in RD_KEYS:
        add(f'admet.{sp}.{k_}', ad['rdkit_descriptors'][jk], u_, pr, f'{src}:rdkit_descriptors.{jk}', f'{lab}: {what}', case=case)
    for rule, rs in ad['rule_summary'].items():
        add(f'admet.{sp}.{rule.lower()}.n_viol', rs['n_violations'], None, 0, f'{src}:rule_summary.{rule}.n_violations',
            f'{lab}: violations of the {rule} filter' + (f' ({"; ".join(rs["violations"])})' if rs['violations'] else ''), case=case)
    for ep, e in ad['admet_ai'].items():
        if e['category'] == 'Physicochemical':
            continue
        u_ = AD_UNITS.get(e['units'], e['units'])
        kind = 'probability (classifier output)' if e['task'] == 'classification' else f'predicted value ({e["units"]})'
        pairs_ = [(f'admet.{sp}.ai.{ep}', e['value'], u_, 2, 'value', f'{lab}, ADMET-AI {e["name"]}: {kind}'),
                  (f'admet.{sp}.ai.{ep}.pct', e['drugbank_percentile'], '%', 0, 'drugbank_percentile',
                   f'{lab}, ADMET-AI {e["name"]}: percentile among the approved DrugBank drugs')]
        for k2, v2, uu, pp, fld, what2 in pairs_:
            if ep in AD_PRINT:
                add(k2, v2, uu, pp, f'{src}:admet_ai.{ep}.{fld}', what2, case=case)
            else:
                addbulk(k2, v2, uu, pp, f'{src}:admet_ai.{ep}.{fld}', what2, f'admet.{sp}.ai',
                        f'ADMET-AI endpoints of the {lab} not printed in Table S14 (full list in admet.json)', case=case)
        ts = e.get('training_set')
        if ts is not None:
            add(f'admet.{sp}.train.{ep}', ts['sulfanilamide_in_set'], None, None, f'{src}:admet_ai.{ep}.training_set (TDC {ts["dataset"]})',
                f'whether sulfanilamide is in the TDC training set of {ep}', case=case)
            if ts['sulfanilamide_in_set'] and ts['label'] is not None:
                add(f'admet.{sp}.train.{ep}.label', ts['label'], u_ if e['task'] == 'regression' else None, 3 if e['task'] == 'regression' else 0,
                    f'{src}:tdc_training_sets.{ts["dataset"]}.rows (Harvard Dataverse file {J(RES / case / "admet/tdc_dataverse_ids.json").get(ts["dataset"], [None])[0]})',
                    f'training label of sulfanilamide in TDC {ts["dataset"]}', case=case)
    if sp == 'sfa':
        n_in = sum(1 for e in ad['admet_ai'].values() if (e.get('training_set') or {}).get('sulfanilamide_in_set'))
        n_ep = sum(1 for e in ad['admet_ai'].values() if e['category'] != 'Physicochemical')
        add('admet.sfa.n_endpoints', n_ep, None, 0, src + ':admet_ai (non-physicochemical endpoints)', 'ADMET-AI endpoints (TDC datasets)')
        add('admet.sfa.n_in_training', n_in, None, 0, src + ':admet_ai[*].training_set.sulfanilamide_in_set', 'endpoints whose training set contains sulfanilamide')
        st_ = ad['admet_ai_vs_shipped_drugbank_row']
        add('admet.sfa.vs_shipped', st_['max_abs_diff_vs_shipped_prediction'], None, 6, src + ':admet_ai_vs_shipped_drugbank_row',
            f'max |difference| from the ADMET-AI prediction shipped for {st_["drugbank_id"]} ({st_["n_compared"]} values)')
        add('admet.n_drugbank', int(re.search(r'(\d+) molecules', ad['provenance']['admet_ai_reference_set']).group(1)), None, 0,
            src + ':provenance.admet_ai_reference_set', 'approved DrugBank drugs in the ADMET-AI reference set (percentiles)')
_pv = J(RES / 'sulfanilamide_pcm_d3' / 'admet' / 'admet.json')['provenance']['versions']
for _k in ('admet_ai', 'rdkit', 'chemprop', 'torch'):
    add(f'admet.version.{_k}', _pv[_k], None, None, f'molecules/sfa/results/sulfanilamide_pcm_d3/admet/admet.json:provenance.versions.{_k}', f'{_k} version of the ADMET run')
# Rule thresholds (printed in the table caption).
add('admet.rule.lipinski', {'MW': 500, 'logP': 5, 'HBD': 5, 'HBA': 10}, None, 0, 'lipinski2001experimental (the rule-of-five thresholds as implemented in RDKit / ADMET-AI physchem.py)', 'Lipinski rule-of-five limits (at most one violation allowed)')
add('admet.rule.veber', {'rotb': 10, 'TPSA': 140}, None, 0, 'veber2002molecular (abstract: "10 or fewer rotatable bonds", "polar surface area equal to or less than 140 A2")', 'Veber limits')
add('admet.rule.ghose', {'logP_lo': -0.4, 'logP_hi': 5.6, 'MR_lo': 40, 'MR_hi': 130, 'MW_lo': 160, 'MW_hi': 480, 'atoms_lo': 20, 'atoms_hi': 70}, None, 1,
    'ghose1999knowledge (abstract: qualifying ranges of the CMC drugs)', 'Ghose qualifying ranges')
add('admet.rule.egan', {'logP': 5.88, 'TPSA': 131.6}, None, 2, 'egan2000prediction model; the numeric limits are those of its SwissADME implementation (daina2017swissadme)',
    'Egan limits', status='threshold values of the SwissADME implementation of the Egan model')
# Experimental values of sulfanilamide (literature). kim2025pubchem = PubChem 2025 update (source of the compiled records);
# knox2024drugbank = DrugBank 6.0 (DB00259 as deposited in PubChem CID 5333). Records fetched with PUG-View on 28 Sep 2026.
PC = 'PubChem CID 5333 (PUG-View, 28 Sep 2026)'
MW_sfa = FACTS['admet.sfa.MW']['value'] if 'admet.sfa.MW' in FACTS else 172.21
add('admet.exp.logP', -0.62, None, 2, f'{PC}: DrugBank DB00259 record [knox2024drugbank] "LogP -0.62 (HANSCH,C ET AL. (1995))"; HSDB 223 "log Kow = -0.62" (Hansch, Leo & Hoekman, Exploring QSAR, ACS 1995, p. 22)',
    'sulfanilamide: measured octanol/water log P')
add('admet.exp.solubility', 7.5, 'g L-1', 1, f'{PC}: DrugBank DB00259 [knox2024drugbank] water solubility 7500 mg/L (Merck Index 1976); HSDB 223 "7.5 g/L at 25 °C" (Merck Index 1989, p. 1409)',
    'sulfanilamide: water solubility at 25 °C')
add('admet.exp.T_sol', 25, '°C', 0, f'{PC}: HSDB 223 and DrugBank DB00259 (water solubility at 25 °C)', 'temperature of the water-solubility value')
add('admet.exp.logS', math.log10(7.5 / MW_sfa), 'log(mol L-1)', 2, 'log10(admet.exp.solubility / admet.sfa.MW)', 'sulfanilamide: log S from the 7.5 g/L solubility')
add('admet.exp.solubility_aquasol', 8.36, 'g L-1', 2, f'{PC}: HSDB 223 "Water solubility = 8.36X10+3 mg/l at 25 °C" (Yalkowsky & Dannenfelser, AQUASOL database, 1992)',
    'sulfanilamide: water solubility at 25 °C (second compilation)')
add('admet.exp.pKa_SO2NH2', 10.43, None, 2, f'{PC}: HSDB 223 "pKa = 10.43, pKb = 11.63" (Merck Index 1989, p. 1409)', 'sulfanilamide: acid dissociation of the sulfonamide N–H')
add('admet.exp.pKa_SO2NH2_20C', 10.6, None, 1, f'{PC}: DrugBank DB00259 [knox2024drugbank] "10.6 (at 20 °C)" (Serjeant & Dempsey 1979); HSDB 223 10.58 (same source)',
    'sulfanilamide: sulfonamide pKa at 20 °C')
add('admet.exp.T_pKa', 20, '°C', 0, f'{PC}: DrugBank DB00259 "10.6 (at 20 °C)"', 'temperature of the Serjeant & Dempsey sulfonamide pKa')
# Anilinium pKa (author decision 29 Sep 2026, paper/audit/research/pka.md): the direct aqueous UV values of Perrin's IUPAC
# compilation replace 14 - Merck pKb 11.63 = 2.37 (Bell & Roblin 1942, indirect, rated Uncertain).
add('admet.exp.pKa_anilinium', 2.02, None, 2, 'greenbaum1955comparison (UV absorbance in dilute mineral acid, 5×10⁻⁵ M, 25 °C) as evaluated in '
    'Perrin 1965 (IUPAC), digitised by zheng2026digitized (dataset v2.4, entry perrin613, rated Approximate, ±0.04); rows read 29 Sep 2026; '
    'the value is cited as evaluated in that compilation', 'sulfanilamide: pKa of the conjugate acid (anilinium, the SAPH cation) in water at 25 °C')
add('admet.exp.T_pKa_anilinium', 25, '°C', 0, 'greenbaum1955comparison via zheng2026digitized (perrin613)', 'temperature of the anilinium pKa')
add('admet.exp.pKa_anilinium_klotz', 2.058, None, 2, 'klotz1945isoelectric (UV absorbance, 6×10⁻⁵ M, 24 °C) via zheng2026digitized (perrin613, rated Approximate); the value as evaluated in that compilation',
    'sulfanilamide: anilinium pKa, second direct aqueous measurement')
add('admet.exp.T_pKa_anilinium_klotz', 24, '°C', 0, 'klotz1945isoelectric via zheng2026digitized (perrin613)', 'temperature of the second anilinium pKa')
add('admet.exp.davis.pK', 10.5, None, 1, 'davis1943binding p. 759 (full text read): "from pK 10.5 for sulfanilamide", quoting Fox & Rose and Schmelkes et al.', 'sulfanilamide pK quoted by Davis 1943')
add('admet.exp.ppb.ratio', 4.6, '%', 1, 'davis1943binding Tables I and VI (full text read): ratio of bound to free drug per g/100 mL of plasma protein, human plasma, pH 7.4, 1.72e-4 M (10 mg/100 mL)',
    'sulfanilamide: bound/free ratio per gram per cent of protein (equilibrium dialysis at refrigerator temperature)')
add('admet.exp.ppb.drug', 10, 'mg/100 mL', 0, 'davis1943binding Table VI (plasma level 10 mgm. per cent = 1.72e-4 M)', 'sulfanilamide plasma level of the Davis 1943 binding estimate')
add('admet.exp.ppb.protein', 6, 'g/100 mL', 0, 'davis1943binding Table VI (typical clinical condition: 6 per cent protein)', 'plasma protein concentration of the Davis 1943 calculation')
_rb = 4.6 / 100 * 6
add('admet.exp.ppb.bound', 100 * _rb / (1 + _rb), '%', 0, 'derived: bound/free = admet.exp.ppb.ratio × admet.exp.ppb.protein = 0.276; bound fraction = 0.276/1.276',
    'sulfanilamide: fraction bound to plasma protein at 10 mg/100 mL and 6 g/100 mL protein (Davis 1943)')
add('admet.exp.ppb.dialysate_plasma', 82, '%', 0, 'davis1943binding Table VI: dialysate/plasma 82 per cent', 'sulfanilamide: dialysate to plasma concentration ratio (human plasma, pH 7.4)')
add('admet.exp.oral_absorbed', [70, 100], '%', 0, f'{PC}: HSDB 223 "Approximately 70 to 100% of an oral dose is absorbed" (sulfonamides as a class; Goodman & Gilman, 8th ed., 1990, p. 1049)',
    'sulfonamides (class statement): fraction of an oral dose absorbed')
add('admet.exp.blood_6h', [4, 8], 'mg/100 mL', 0, 'taylor1940comparative p. 208 (full text read): "Values of 4 to 8 mgm. of free sulfapyridine or sulfanilamide per 100 ml. of blood were obtained 6 hours after the beginning of ingestion"',
    'sulfanilamide: free drug in blood 6 h after an oral dose (determinable at 24 h)')
add('admet.exp.blood_6h_time', 6, 'h', 0, 'taylor1940comparative p. 208: "6 hours after the beginning of ingestion"', 'time after the oral dose of the Taylor 1940 blood level')
# Human half-life (author decision 29 Sep 2026, paper/audit/research/halflife.md): derived by us as ln 2 / k from urinary-
# excretion rate constants; neither source prints a half-life.
add('admet.exp.k_el', 0.067, 'h-1', 3, 'nogami1968studies Fig. 3 (full text read; human data of Kakemi et al. 1965): overall elimination rate constant '
    'K = k1 (acetylation) + k3 (renal) of sulfanilamide, read off the figure (reading error about ±0.3 h in t½)', 'sulfanilamide: overall elimination rate constant in man')
add('admet.exp.k3', 0.0317, 'h-1', 4, 'koizumi1964absorption Table I (full text read): K3 (Nelson), sulfanilamide, healthy adult men, 1 g oral',
    'sulfanilamide: urinary excretion rate constant of the unchanged drug in man')
add('admet.exp.f_unchanged', 0.48, None, 2, 'nogami1968studies Fig. 2 (read off the figure)', 'sulfanilamide: fraction excreted unchanged in man')
add('admet.exp.dose_oral', 1, 'g', 0, 'koizumi1964absorption p. 430: "Adult men in apparent good health ingested 1 g. of finely powdered drug"',
    'oral dose of the human excretion study')
add('admet.exp.half_life', math.log(2) / 0.067, 'h', 0, 'derived: ln 2 / admet.exp.k_el = 10.3 h (our conversion); cross-check ln 2 / (admet.exp.k3 / '
    'admet.exp.f_unchanged) = 10.5 h', 'sulfanilamide: elimination half-life in man (about 10 h)', check=math.log(2) / (0.0317 / 0.48), tol=0.3)
add('admet.exp.LD50_rat_oral_log', 1.645, 'log(kg mol-1)', 3, 'TDC LD50_Zhu label of sulfanilamide (zhu2009quantitative; admet.sfa.train.LD50_Zhu.label)', 'sulfanilamide: rat oral LD50 as −log(LD50 / mol kg⁻¹)')
add('admet.exp.LD50_rat_oral', 10 ** (-1.645) * MW_sfa, 'g kg-1', 1, 'derived: 10^−1.645 mol/kg × admet.sfa.MW (g/mol) = g/kg', 'sulfanilamide: rat oral LD50 (Zhu et al. data set)')
add('admet.exp.DILI', 'DILI positive', None, None, f'{PC}: FDA Liver Toxicity Knowledge Base, DILIst classification "DILI Positive" (thakkar2020drug); TDC DILI label 1',
    'sulfanilamide: human drug-induced liver injury classification')
add('admet.exp.AMES', 'negative', None, None, 'TDC AMES label 0 for sulfanilamide (hansen2009benchmark)', 'sulfanilamide: Ames mutagenicity')
CAVEATS.append(f'ADMET (section 25): sulfanilamide is in the training data of {N_TR} of the {N_EP} ADMET-AI endpoints, so agreement there '
               '(solubility, LD50, AMES, DILI, CYP inhibition) is recall, not validation. Where it is not (plasma protein binding, '
               'half-life, bioavailability, lipophilicity), the model can be tested: PPBR 54 % predicted against about 22 % bound '
               '(Davis 1943, human plasma, dialysis at refrigerator temperature). Crippen cLogP (−0.08) is 0.5 above the measured '
               'log P (−0.62). The predicted half-life (16.3 h) is longer than the about 10 h derived from urinary excretion in man '
               '(Koizumi 1964, Nogami 1968; t½ = ln 2/k is our conversion).')
MISSING[:] = [m for m in MISSING if not m.startswith('No ΔSCF for the sulfonamidate')]
MISSING.append('No ΔSCF for the sulfonamidate, H2PO4− or the gas-phase pair. Condensed Fukui functions exist for sulfanilamide and the '
               'cation in IEF-PCM (section 24), not for the ion pair, the sulfonamidate or the gas phase.')

# =====================================================================================================================
section('26. Figure rendering parameters (fig.*)',
        'Drawing choices printed in the figure captions (isovalues, colour ranges, cut-offs, bins): figure_params.json '
        '(group B) and figure_params_groupA.json (group A) merged as fig.<name>, plus the output resolution of figkit. '
        'These are settings, not results.')
for pf in ['paper/main/figures/figure_params_groupA.json', 'paper/main/figures/figure_params.json']:
    for k_, e in J(REPO / pf).items():
        if k_ == '_meta':
            continue
        add(f'fig.{k_}', e['value'], e.get('unit') or None, e.get('precision'), f"{pf}:{k_} ({e.get('source', '')})", e['meaning'])
_dpi = re.search(r'^DPI\s*=\s*(\d+)', (REPO / 'paper/_shared/figkit/style.py').read_text(), re.M)
add('fig.dpi', int(_dpi.group(1)), 'dpi', 0, 'paper/_shared/figkit/style.py: DPI', 'resolution of every PNG written by figkit (vector PDF alongside)')

# =====================================================================================================================
section('27. Docking cluster representatives of Fig. 10 (sulfanilamide)',
        'The pose drawn for each receptor in Fig. 10: the top pose of the best-scoring seed of the largest cluster '
        '(dock/results/production_contacts.json, cluster "1"). Every distance re-derived from the pose SDF and '
        'receptor.pdb by derive/d10_rep_contacts.py (✓ = equal to 0.01 Å). Heavy-atom distances; N/O pairs within 3.5 Å.')
D10 = J(DER / 'd10_rep_contacts.json')
RN = {'HH2': 'pterin (HH2)', 'POP': 'PPi'}
for job, d10 in D10['jobs'].items():
    t_, tag_ = job.split('/')
    pdb_ = t_.split('_')[1]
    sc_ = tag_.split('__')[1]
    base = f'dock.rep.{pdb_}.sfa.{sc_}'
    pcsrc = f'dock/results/production_contacts.json:["{job}"]["1"]'
    add(base + '.seed', d10['seed'], None, 0, pcsrc + '.seed', f'{pdb_} {sc_}: seed of the representative pose drawn in Fig. 10')
    add(base + '.population', d10['population'], None, 0, pcsrc + '.population', f'{pdb_} {sc_}: seeds whose top pose falls in the largest cluster')
    add(base + '.score', d10['score_source'], 'kcal mol-1', 1, pcsrc + '.score', f'{pdb_} {sc_}: docking score of the representative pose',
        check=d10['score_tag'], tol=0.0005)
    for p in d10['pairs']:
        res_ = RN.get(p['residue'][:3], p['residue'][:1] + p['residue'][1:3].lower() + p['residue'][3:])
        la = p['ligand_atom'].split('#')[0]
        key_ = f"{base}.{p['residue']}_{p['receptor_atom']}_{la.split('(')[1].rstrip(')')}"
        if key_ in FACTS:
            key_ += '_' + p['ligand_atom'].split('#')[1]
        add(key_, p['source_value'], 'A', 2, pcsrc + '.polar_pairs_NO_within_3.5A', f'{pdb_} {sc_} representative: {res_} {p["receptor_atom"]}···{la}',
            check=p['distance'], tol=0.006)
    if d10['zinc']:
        add(base + '.ZnN', d10['zinc']['source_value'], 'A', 2, pcsrc + '.zinc.distance',
            f'{pdb_} {sc_} representative: Zn···N(sulfonamide), the ligand atom nearest to the zinc', check=d10['zinc']['distance'], tol=0.006)

# Rounding ties: decided in facts_burst2.py (section 30; round half to even on the full-precision value).
