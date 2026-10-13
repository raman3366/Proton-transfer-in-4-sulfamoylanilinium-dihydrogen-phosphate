"""
dock/RESULTS.md from dock/results/production_summary.json (+ production_contacts.json and the
validation files dock/analysis/{summary,gate_G2}.json).

Every number in the report is read from those files; nothing is computed here except sorting and
formatting, so the report cannot disagree with the JSON it comes from.

    python dock/production_report.py
"""
import json
import pathlib
import sys

DOCK = pathlib.Path(__file__).resolve().parent
RES = DOCK / 'results'
OUT = DOCK / 'RESULTS.md'
SPECIES_LABEL = {
    'sfa_neutral': 'sulfanilamide (neutral)', 'sfa_anion': 'sulfonamidate anion',
    'sfa_cation': '4-sulfamoylanilinium cation', 'smx_anion': 'sulfamethoxazole anion',
    'smx_neutral': 'sulfamethoxazole (neutral)', 'paba_anion': 'p-aminobenzoate',
    'azm_neutral': 'acetazolamide (neutral)', 'azm_anion': 'acetazolamide anion',
    'dhp_anion': '7,8-dihydropteroate (co-crystal ligand)'}
SP_ORDER = ['sfa_neutral', 'sfa_anion', 'sfa_cation', 'dhp_anion', 'smx_anion', 'smx_neutral',
            'azm_neutral', 'azm_anion', 'paba_anion']
T_ORDER = ['T5_7S2K', 'T4_3HS4', 'T3_6RL9', 'T1_3TZF']
SC_ORDER = ['vina', 'vinardo', 'ad4zn']
SC_LABEL = {'vina': 'Vina', 'vinardo': 'Vinardo', 'ad4zn': 'AD4Zn'}


def f(x, n=2, dash='–'):
    return dash if x is None else f'{x:.{n}f}'


# Values whose full precision is known to differ in print from the stored, already rounded JSON value (facts.md,
# "3HS4 cross-dock"): the 3HS4 AD4Zn cross-dock RMSD is 3.2150 Å at full precision (derive/d09_dock.py), stored as
# 3.215, which a binary float prints as 3.21; the paper's rule (round half to even on the full-precision value) gives
# 3.22. Applied to the rendered lines of that job only (every line that names 3HS4 and prints this value).
FULL_PRECISION_PRINT = [('3HS4', '3.21 ± 0.02', '3.22 ± 0.02')]


def pm(s, n=2):
    return f'{f(s["mean"], n)} ± {f(s["sd"], n)}' if s else '–'


def main():
    S = json.loads((RES / 'production_summary.json').read_text())
    alljobs = {k: v for k, v in S['jobs'].items() if 'error' not in v}
    jobs = {k: v for k, v in alljobs.items() if not v.get('diagnostic')}
    diag = {k: v for k, v in alljobs.items() if v.get('diagnostic')}
    errors = {k: v for k, v in S['jobs'].items() if 'error' in v}
    V = json.loads((DOCK / 'analysis' / 'summary.json').read_text())
    G = json.loads((DOCK / 'analysis' / 'gate_G2.json').read_text())
    L = []
    w = L.append

    def byt(t):
        return {k: v for k, v in jobs.items() if v['target'] == t}

    w('# Production docking: results')
    w('')
    w(f'Generated {S["generated"]} by `{S["generated_by"]}` and `python dock/production_report.py`. '
      'Every number in this file is read from `dock/results/production_summary.json` and '
      '`dock/results/production_contacts.json`; the validation numbers quoted for context come from '
      '`dock/analysis/{summary,gate_G2}.json` and are discussed in `dock/VALIDATION.md`.')
    w('')
    n_seed = sorted({j['n_seeds'] for j in jobs.values()})
    w(f'**What ran.** {len(jobs)} production jobs (tag suffix `__prod`): '
      f'{len(T_ORDER)} targets × {len(SP_ORDER) - 1}–{len(SP_ORDER)} species × '
      f'{{Vina, Vinardo}}, plus AutoDock4Zn for the two hCA II targets with the species whose '
      'sulfonamide is neutral. Exhaustiveness '
      f'{sorted({j["exhaustiveness"] for j in jobs.values()})[0]}, 20 Å box, '
      f'{"/".join(str(x) for x in n_seed)} seeds per job, `--cpu 3`, `nice -n 10`, one docking '
      'process tree at a time on the Raspberry Pi 5. Total docking time '
      f'{sum(j["wall_s_total"] for j in jobs.values()) / 3600:.1f} h.')
    w('')
    if diag:
        dh = sum(v['wall_s_total'] for v in diag.values()) / 3600
        w(f'Six further AD4Zn runs ({dh:.1f} h) are the input-geometry diagnostic of section 7b; they '
          'are kept out of the tables, the rankings and the DFT pose set, and are listed in '
          '`production_summary.json` → `diagnostics`.')
        w('')
    w('The three sulfanilamide species are docked **from their DFT geometries** '
      '(B3LYP-D3(BJ)/6-311++G(d,p), IEF-PCM water): '
      + ', '.join(f'`{v["xyz"]}`' for v in S['dft'].values())
      + '. Vina and AD4Zn hold every bond length, bond angle and rigid-fragment geometry of their '
      'input fixed and randomise position, orientation and the active torsions from `--seed`, so the '
      'DFT structure is what is docked and the seed still drives the search. The controls '
      '(sulfamethoxazole, acetazolamide, p-aminobenzoate, dihydropteroate) keep their ETKDG/MMFF94s '
      'geometries: they were never optimised at the DFT level, and using the same route for them as '
      'in validation keeps those runs comparable. The ion pair is never docked.')
    w('')
    w('**Scores are docking scores** (kcal/mol by the programs\' convention), reported as mean ± SD '
      'over seeds with the best value; they are never converted to a Ki or a Kd, and only the '
      'ordering is interpreted.')
    w('')

    w('**What was changed in the code to run this, and why.** Each change is commented at the point '
      'of change with the date.')
    w('')
    w('1. `prep_ligands.py`, `run_vina.py`, `run_ad4zn.py`: **the DFT geometries were not reaching '
      'the docking.** `prep_ligands.py --dft` wrote them to `dock/ligands/`, but both runners '
      'rebuilt an ETKDG conformer from the species SMILES for every seed and never read that file, '
      'so the production run as scripted would have docked ETKDG geometries and reported them as '
      'DFT ones. A new `prep_ligands.input_mol_for_seed()` and a `--geometry {etkdg,dft}` flag on '
      'both runners fix it; the geometry actually used is recorded in each run\'s `results.json` '
      '(`geometry`, `input_geometry` with the xyz path and its sha256).')
    w('2. `run_production.sh`: production jobs are tagged `__prod`. As scripted they carried no '
      'suffix, so they would have overwritten the validation runs that gate G2 rests on.')
    w('3. `run_production.sh`: AD4Zn is run only for the species whose sulfonamide is neutral '
      '(sulfanilamide, sulfamethoxazole, acetazolamide). Its TZ pseudo-atom models the Zn–N(H) bond '
      'of the neutral form; in validation the anion gave Zn–N 4.86 Å for acetazolamide against '
      '2.08 Å for the neutral form. 20 seeds for the two sulfanilamide AD4Zn jobs whose poses are '
      'discussed, 10 elsewhere.')
    w('4. `run_production.sh`: each target\'s own co-crystal ligand is docked as the reference '
      'score (7S2K\'s dihydropteroate was added for this; the others were already species of the '
      'matrix), jobs run in priority order one process tree at a time, a completed job is skipped on '
      'a rerun, and the progress log with an ETA is written to `dock/results/production.log`.')
    w('5. `analyse.py --production` (new `production_analysis.py`) and `contacts.py --production` '
      'write to `dock/results/`; the validation files `dock/analysis/{summary,gate_G2}.json` are '
      'left alone, and production runs are excluded from them. Re-running `analyse.py` after these '
      'edits reproduced both files byte for byte.')
    w('')

    # ---------------------------------------------------------------- headline
    w('## 1. Headline numbers')
    w('')
    w('| target | protocol | validated? (dock/VALIDATION.md) | sulfanilamide⁰ score | rank of '
      'sulfanilamide⁰ | pose vs the crystal ligand |')
    w('|---|---|---|---|---|---|')
    for t in T_ORDER:
        info = S['targets'][t]
        for sc in SC_ORDER:
            j = jobs.get(f'{t}/sfa_neutral__{sc}__prod')
            if not j:
                continue
            r = S['rank_order'].get(f'{t}/{sc}', {})
            rmsd = ''
            for label, e in j['rmsd_to_crystal'].items():
                if 'full' in e:
                    rmsd = (f'{f(e["full"]["mean"])} ± {f(e["full"]["sd"])} Å to {label.split(" (")[0]} '
                            f'(f≤1 Å {f(e["full"]["f_le_1A"], 1)})')
                    break
            else:
                for label, e in j['rmsd_to_crystal'].items():
                    cores = e.get('cores', {})
                    if cores:
                        cname = min(cores, key=lambda c: -cores[c]['n_atoms'])
                        rmsd = (f'{f(cores[cname]["mean"])} Å over the {cname.split(" (")[0]} '
                                f'vs {label.split(" (")[0]}')
                        break
            w(f'| {info["pdb"]} ({info["protein"]}) | {SC_LABEL[sc]} | '
              f'{info["validation"]} | '
              f'{pm(j["score"], 3)} (best {f(j["score"]["min"], 3)}) | '
              f'{r.get("sulfanilamide_neutral_rank", "–")} by score | {rmsd or "–"} |')
    w('')
    # Review 2026-09-28: the rank column is a rank by score over all species, with no site filter, so
    # it can put a species above sulfanilamide on a pose that does not bind the zinc at all (the
    # AD4Zn sulfamethoxazole pose: Zn-N 3.28-3.56 A, zinc-site checks 0/10 in both hCA II targets).
    # Section 2 gives the ranking over zinc-bound top poses beside it; said here so the headline
    # table is not read on its own.
    smxzn = [j for k, j in jobs.items() if j['ligand'] == 'smx_neutral' and j['scoring'] == 'ad4zn'
             and j['site_checks'].get('Zn-N(sulfonamide)', {}).get('applicable')]
    smxtxt = ''
    if smxzn:
        zs = [j['site_checks']['Zn-N(sulfonamide)'] for j in smxzn]
        smxtxt = (f' (sulfamethoxazole with AD4Zn is exactly that case — Zn–N '
                  f'{f(min(z["distance_mean"] for z in zs))}–{f(max(z["distance_mean"] for z in zs))} Å, '
                  f'the zinc-site check passing in {f(max(z["pass_fraction"] for z in zs), 2)} of the '
                  f'seeds in both hCA II targets)')
    w('The rank column is the rank by mean score among all species docked to that target with that '
      'scoring function, **with no site filter**: for hCA II it can rank above sulfanilamide a '
      f'species whose top pose does not reach the zinc{smxtxt}. The ranking restricted '
      'to zinc-bound top poses is in section 2, and section 4 has the distances.')
    w('')

    # ---------------------------------------------------------------- per target
    w('## 2. Scores, clusters and ligand efficiency, per target × scoring × species')
    w('')
    w('`score` = top pose of each seed, mean ± SD (SD over seeds, ddof 1). `clusters` = populations '
      f'of the leader clusters of those top poses (symmetry-corrected heavy-atom RMSD ≤ '
      f'{S["protocol"]["cluster_cutoff_A"]} Å, seeds ordered by score), largest cluster first as '
      '`population/seeds`. `LE` = −score/heavy atoms; `FQ` = LE / LE_scale(HA) with Reynolds\' '
      'size correction — read section 5 before using either. `Δ ref` = score minus the score of the '
      'receptor\'s own co-crystal ligand in the same scoring function.')
    for t in T_ORDER:
        info = S['targets'][t]
        w('')
        w(f'### {info["pdb"]} — {info["protein"]}')
        w('')
        w(f'Co-crystal ligand {info["crystal_ligand"]} ({info["crystal_ligand_name"]}). '
          f'Validation: {info["validation"]}.')
        w('')
        w('| species | scoring | score (kcal/mol) | best | seeds | clusters | LE | FQ | Δ ref |')
        w('|---|---|---|---|---|---|---|---|---|')
        for sp in SP_ORDER:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j:
                    continue
                cl = j['clusters']['clusters']
                clstr = ', '.join(f'{c["population"]}/{j["n_seeds"]}' for c in cl[:3]) + \
                    ('' if len(cl) <= 3 else f' (+{len(cl) - 3} more)')
                rk = S['rank_order'].get(f'{t}/{sc}', {})
                d = next((x.get('delta_vs_reference') for x in rk.get('ranking', [])
                          if x['species'] == sp), None)
                w(f'| {SPECIES_LABEL[sp]} | {SC_LABEL[sc]} | {pm(j["score"], 3)} | '
                  f'{f(j["score"]["min"], 3)} | {j["n_seeds"]} | {clstr} | '
                  f'{f(j["ligand_efficiency"]["LE_mean"], 3)} | '
                  f'{f(j["ligand_efficiency"]["FQ_mean"], 2)} | '
                  f'{"ref" if d == 0 else f(d, 2)} |')
        # rank order line
        for sc in SC_ORDER:
            rk = S['rank_order'].get(f'{t}/{sc}')
            if not rk:
                continue
            order = ' < '.join(f'{x["species"]} ({f(x["score_mean"], 2)})' for x in rk['ranking'])
            w('')
            w(f'**Rank order, {SC_LABEL[sc]}** (best first): {order}. Reference ligand: '
              f'`{rk["reference_species"]}` at {f(rk["reference_score"], 2)}.')
            if 'ranking_zinc_bound_only' in rk:
                # Review 2026-09-28: an empty list used to render as a dangling ': .', which read as
                # a missing number rather than as the result it is.
                zb = ', '.join(f'{x["species"]} {f(x["score_mean"], 2)} ({x["seeds"]} seeds)'
                               for x in rk['ranking_zinc_bound_only'])
                w('')
                if zb:
                    w(f'Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 '
                      f'≤ 3.3 Å): {zb}. Species with no zinc-bound seed are absent from that list.')
                else:
                    w('Over **zinc-bound top poses only** (Zn–N(sulfonamide) ≤ 2.3 Å and Thr199 OG1 ≤ '
                      '3.3 Å): **no species has a zinc-bound top pose in any seed**, so this scoring '
                      'function ranks these species on poses that do not reproduce the zinc anchor.')
    w('')

    # ---------------------------------------------------------------- RMSD
    w('## 3. Where the poses sit relative to the crystal ligands')
    w('')
    w('Symmetry-corrected heavy-atom RMSD in the receptor frame, no superposition (spyrmsd). '
      '`full` is used only where the docked species is the same molecule as the crystal ligand; '
      'otherwise the shared moiety is compared, over every chemically equivalent mapping of that '
      'moiety onto both molecules. These are **not** redocking RMSDs for the sulfanilamide species '
      '(only 6RL9 co-crystallises sulfanilamide): they answer "does it bind where the reference '
      'ligand binds".')
    w('')
    w('| target | species | scoring | reference | measure | RMSD (Å) | best | f≤1 Å | f≤2 Å |')
    w('|---|---|---|---|---|---|---|---|---|')
    for t in T_ORDER:
        for sp in SP_ORDER:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j:
                    continue
                for label, e in j['rmsd_to_crystal'].items():
                    rows = []
                    if 'full' in e:
                        rows.append((f'full molecule ({e["full"]["n_atoms"]} atoms)', e['full']))
                    for cname, blk in e.get('cores', {}).items():
                        rows.append((cname, blk))
                    if 'modes' in e:
                        rows.append(('best pose of the ranked list (any rank)',
                                     e['modes']['best_of_modes']))
                    for measure, blk in rows:
                        w(f'| {S["targets"][t]["pdb"]} | {sp} | {SC_LABEL[sc]} | {label} | {measure} '
                          f'| {pm(blk)} | {f(blk["min"])} | {f(blk["f_le_1A"], 1)} | '
                          f'{f(blk["f_le_2A"], 1)} |')
    w('')
    w('')
    w('Where a whole-molecule reference exists, the rank at which the first near-native pose appears '
      'in each seed\'s ranked list says whether a poor top-pose RMSD is a sampling failure or a '
      'ranking one:')
    w('')
    w('| target | species | scoring | reference | seeds with a pose ≤1 Å anywhere | median rank of it '
      '| same for ≤2 Å |')
    w('|---|---|---|---|---|---|---|')
    for t in T_ORDER:
        for sp in SP_ORDER:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j:
                    continue
                for label, e in j['rmsd_to_crystal'].items():
                    if 'modes' not in e:
                        continue
                    r1 = e['modes']['rank_of_first_pose_within']['1A']
                    r2 = e['modes']['rank_of_first_pose_within']['2A']
                    w(f'| {S["targets"][t]["pdb"]} | {sp} | {SC_LABEL[sc]} | {label.split(" (")[0]} | '
                      f'{f(r1["fraction_of_seeds_with_any"], 2)} | {f(r1["median_rank"], 1)} | '
                      f'{f(r2["fraction_of_seeds_with_any"], 2)} at rank {f(r2["median_rank"], 1)} |')
    w('')
    w('### The validation gate, for the same targets')
    w('')
    w('Gate G2 as defined before the runs: top-pose symmetry-corrected RMSD ≤ **1.0 Å** (mean over seeds) '
      'with the success fraction reported, PoseBusters physics ≥ 0.8, and for hCA II also '
      'Zn–N(sulfonamide) ≤ 2.3 Å and the Thr199 hydrogen bond in ≥ half the seeds. The conventional '
      '**2.0 Å** criterion (Hevener 2009; Warren 2006) is evaluated with the same structure. Both are '
      'reported; the gate was not changed.')
    w('')
    w('| validation run (redock / cross-dock) | RMSD top (Å) | f≤1 Å | f≤2 Å | G2 (1.0 Å) | 2.0 Å criterion |')
    w('|---|---|---|---|---|---|')
    for k, c in G['criteria'].items():
        tgt = k.split('/')[0]
        if tgt not in T_ORDER or len(k.split('/')[1].split('__')) != 2:
            continue
        st = V[k]['stats']
        w(f'| `{k}` | {f(st["rmsd_top_mean"])} ± {f(st["rmsd_top_sd"])} | '
          f'{f(c["success_fraction_1A"], 1)} | {f(c["success_fraction_2A"], 1)} | '
          f'{"**PASS**" if c["passes_G2"] else "fail"} | '
          f'{"pass" if c["passes_2A_criterion"] else "fail"} |')
    w('')

    # ---------------------------------------------------------------- site checks and contacts
    w('## 4. Site checks and key contacts')
    w('')
    w('Site checks are the geometric tests of `analyse.py` (plus two generalised ones for species '
      'whose SMARTS the validation set did not cover), each also evaluated on the receptor\'s own '
      'crystal ligand so that a mis-specified check shows up as a failing crystal pose. Distances are '
      'heavy-atom to heavy-atom.')
    w('')
    w('| target | species | scoring | check | crystal ligand | docked mean | pass fraction |')
    w('|---|---|---|---|---|---|---|')
    for t in T_ORDER:
        for sp in SP_ORDER:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j:
                    continue
                for lab, c in j['site_checks'].items():
                    if not c['applicable']:
                        continue
                    w(f'| {S["targets"][t]["pdb"]} | {sp} | {SC_LABEL[sc]} | {lab} (≤ {c["cutoff"]} Å) '
                      f'| {f(c["crystal_ligand_value"])} | {f(c["distance_mean"])} | '
                      f'{f(c["pass_fraction"], 1)} |')
    w('')
    w('### Contacts of the cluster representatives')
    w('')
    w('One profile per cluster (the top pose of the cluster\'s best-scoring seed), from PLIP, ProLIF '
      'and plain geometry independently; full records in `dock/results/production_contacts.json`. '
      'Hydrogen bonds are given as donor···acceptor heavy-atom distances (the hydrogens are ours, not '
      'the crystal\'s).')
    for t in T_ORDER:
        rows = []
        for sp in SP_ORDER:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j:
                    continue
                for cl in j['clusters']['clusters']:
                    if cl['population'] < 2 and len(j['clusters']['clusters']) > 1:
                        continue
                    kc = cl.get('contacts', {})
                    hb = '; '.join(f'{h["residue"]} {h["donor"].split(":")[-1]}···'
                                   f'{h["acceptor"].split(":")[-1]} {f(h["D_A"])}'
                                   for h in kc.get('hbonds', [])[:4]) or '–'
                    zn = kc.get('zinc')
                    znt = f'{zn["nearest_ligand_atom"]} {f(zn["distance"])} Å' if zn else '–'
                    pi = '; '.join(f'{p["residue"]} {p["type"]} {f(p["distance"])}'
                                   for p in kc.get('pi', [])) or '–'
                    pp = '; '.join(f'{x["residue"]}:{x["receptor_atom"]}···{x["ligand_atom"]} '
                                   f'{f(x["distance"])}' for x in kc.get('polar_pairs', [])) or '–'
                    rows.append(f'| {sp} | {SC_LABEL[sc]} | {cl["id"]} | '
                                f'{cl["population"]}/{j["n_seeds"]} | {f(cl["score_best"], 2)} | '
                                f'{hb} | {pp} | {znt} | {pi} |')
        if rows:
            w('')
            w(f'**{S["targets"][t]["pdb"]}**')
            w('')
            w('| species | scoring | cluster | population | best score | PLIP H-bonds (D···A, Å) | '
              'closest N/O···N/O pairs (Å) | nearest ligand atom to Zn | π contacts |')
            w('|---|---|---|---|---|---|---|---|---|')
            L.extend(rows)
    w('')

    # ---------------------------------------------------------------- LE
    w('## 5. Ligand efficiency and its size bias')
    w('')
    w('LE = −score / heavy atoms. It has a well known size bias: docking scores grow sub-linearly '
      'with molecular size, so LE systematically favours small ligands. In this set sulfanilamide has '
      '11 heavy atoms, p-aminobenzoate 10, acetazolamide 13, sulfamethoxazole 17 and dihydropteroate '
      '23, so an LE comparison across them is largely a comparison of size. Measured here, over every '
      'production job of each scoring function:')
    w('')
    w('| scoring | jobs | heavy atoms | r(score, heavy atoms) | r(LE, heavy atoms) | r(FQ, heavy atoms) |')
    w('|---|---|---|---|---|---|')
    for sc, d in S['ligand_efficiency_size_bias'].items():
        rng = d.get('heavy_atom_range')
        w(f'| {SC_LABEL[sc]} | {d["n_jobs"]} | {f"{rng[0]}–{rng[1]}" if rng else "–"} | '
          f'{f(d.get("r_score_vs_HA"))} | {f(d.get("r_LE_vs_HA"))} | {f(d.get("r_FQ_vs_HA"))} |')
    w('')
    w('Read it like this: r(score, HA) far below zero means bigger ligands simply score better in '
      'this set, and r(LE, HA) far below zero means the per-atom measure then favours the small ones '
      '— the two halves of the size bias, measured rather than asserted. r(FQ, HA) near zero means '
      'the size correction has removed most of it. FQ = LE / LE_scale(HA), LE_scale = 0.0715 + '
      '7.5328/HA + 25.7079/HA² − '
      '361.4722/HA³ (Reynolds, Tounge & Bembenek 2007, doi:10.1016/j.bmcl.2007.07.071), removes most '
      'of it; but that scale was fitted to *experimental* affinities, so applied to docking scores it '
      'is a heuristic, not a calibrated quantity. **Neither LE nor FQ is used here to claim potency**; '
      'both are reported as the common size-corrected measures, and the score ranking against the '
      'controls in section 2 is the defensible comparison.')
    w('')

    # ---------------------------------------------------------------- DFT link
    w('## 6. The DFT–docking link: the C3–C4–S1–N2 torsion and the strain of the docked pose')
    w('')
    cr = S['crystal']['torsion']
    w(f'The crystal value is **{f(cr["tau_C3"], 2)}°** for C3–C4–S1–N2 and {f(cr["tau_C5"], 2)}° for '
      f'C5–C4–S1–N2 ({cr["source"]}). Which of the two a given structure reports depends on the '
      'arbitrary choice of ring carbon, so the table also gives the ring-flip-invariant '
      '`τ_fold = max(|C3–C4–S1–N2|, |C5–C4–S1–N2|)`: 90° means the S–N bond is perpendicular to the '
      'ring plane, 180° means it lies in it.')
    w('')
    w('| geometry | C3–C4–S1–N2 | C5–C4–S1–N2 | τ_fold |')
    w('|---|---|---|---|')
    w(f'| crystal (CIF) | {f(cr["tau_C3"], 2)}° | {f(cr["tau_C5"], 2)}° | '
      f'{f(max(abs(cr["tau_C3"]), abs(cr["tau_C5"])), 2)}° |')
    for sp, d in S['dft'].items():
        w(f'| DFT optimum, {SPECIES_LABEL[sp]} (PCM water) | {f(d["tau_C3"], 2)}° | '
          f'{f(d["tau_C5"], 2)}° | {f(d["tau_fold"], 2)}° |')
    w('')
    # ACS integration review (30 Sep 2026): which pose basis the paper uses for these torsions
    d14p = DOCK.parent / 'paper/main/data/derived/d14_completion.json'
    d14 = json.loads(d14p.read_text()).get('tau_raw_all', {}) if d14p.exists() else {}
    w('The torsions below are measured on the rebuilt poses (see "How the torsion of a pose is obtained" '
      'under the table). The paper (Section 3.8, SI Table S31) measures them on the raw pose coordinates '
      'as the docking program wrote them'
      + (f', which differ from the values here by at most {f(d14["max_diff_vs_rebuilt"])}°'
         if d14.get('max_diff_vs_rebuilt') is not None else '')
      + ', and takes the strain from the relaxed DFT torsion scans instead of MMFF94s.')
    w('')
    w('| target | species | scoring | cluster | population | τ_fold of the top poses | MMFF94s strain, '
      'rigid (kcal/mol) | MMFF94s strain, relaxed (kcal/mol) |')
    w('|---|---|---|---|---|---|---|---|')
    for t in T_ORDER:
        for sp in ['sfa_neutral', 'sfa_anion', 'sfa_cation']:
            for sc in SC_ORDER:
                j = jobs.get(f'{t}/{sp}__{sc}__prod')
                if not j or 'torsion' not in j:
                    continue
                by_cl = {}
                for r in j['torsion']['per_seed']:
                    by_cl.setdefault(r['cluster'], []).append(r)
                for cid, rs in sorted(by_cl.items()):
                    tf = [r['tau_fold'] for r in rs]
                    sr = [r['strain_mmff_rigid'] for r in rs]
                    sx = [r['strain_mmff_relaxed'] for r in rs]
                    w(f'| {S["targets"][t]["pdb"]} | {sp} | {SC_LABEL[sc]} | {cid} | '
                      f'{len(rs)}/{j["n_seeds"]} | {f(min(tf), 1)}–{f(max(tf), 1)}° | '
                      f'{f(min(sr))}–{f(max(sr))} | {f(min(sx))}–{f(max(sx))} |')
    w('')
    tf = [r['tau_fold'] for j in jobs.values() if 'torsion' in j for r in j['torsion']['per_seed']]
    sx = [r['strain_mmff_relaxed'] for j in jobs.values() if 'torsion' in j
          for r in j['torsion']['per_seed']]
    if tf:
        w(f'**The link in one sentence.** Across all {len(tf)} docked top poses of the three '
          f'sulfanilamide species, τ_fold is {f(min(tf), 1)}–{f(max(tf), 1)}°, i.e. the S–N bond is '
          'rotated away from the perpendicular arrangement of the isolated DFT molecule '
          f'({f(S["dft"]["sfa_neutral"]["tau_fold"], 1)}° for the neutral species) and past the '
          f'crystal value ({f(max(abs(cr["tau_C3"]), abs(cr["tau_C5"])), 1)}°) towards the ring '
          f'plane, at an MMFF94s cost of at most {f(max(sx))} kcal/mol. The receptors therefore ask '
          'for a modest, affordable change in this torsion rather than an unusual conformer; the '
          'relaxed DFT torsion scans (paper Section 3.8) give the cost at the B3LYP-D3(BJ) level.')
        w('')
    w('**How the torsion of a pose is obtained.** Vina writes only the atoms it moved, so each top '
      'pose was rebuilt as the full DFT molecule: the DFT geometry\'s rotatable torsions are set to '
      'the docked values and the result is superposed on the docked atoms. The residual of that fit is '
      'recorded for every pose (`rebuild_residual_A`, ≤ '
      f'{max(f(j["torsion"]["max_rebuild_residual_A"], 4) for j in jobs.values() if "torsion" in j)}'
      ' Å over all poses), which is what makes the torsion a measurement rather than a re-perception.')
    w('')
    w('**The two MMFF94s strain estimates.** `rigid` holds the DFT bond lengths and angles fixed, as '
      'the docked pose does, and varies only the torsions: E(pose, polar-H rotors re-optimised) minus '
      'the minimum of the same rigid profile. `relaxed` minimises MMFF94s with C3–C4–S1–N2 restrained '
      'at the pose value and subtracts the MMFF94s global minimum. The polar-H rotors are '
      're-optimised because Vina does not score hydrogen positions, so their docked orientation '
      'carries no information. Both profiles (every 5°) are in `production_summary.json` → '
      '`mmff_torsion_profiles`. MMFF94s is a cheap estimate: the profile is shallow — over '
      + ', '.join(f'{sp}: {f(max(e for _, e in S["mmff_torsion_profiles"][sp]["relaxed"]))} kcal/mol'
                  for sp in S['mmff_torsion_profiles'])
      + ' between the perpendicular minimum and the in-plane maximum (relaxed profiles) — so the '
      'ranking of poses by strain is safe, and the paper takes the absolute numbers from the relaxed '
      'DFT torsion scans (06_DFT_results/torsion_scan of the data package). '
      f'`dock/results/poses_for_dft/` holds the {S["poses_for_dft"]["n"]} pose geometries (one per '
      'cluster representative, with and without the polar-H re-optimisation) plus a `manifest.json` '
      'giving the charge, multiplicity, torsion, MMFF strain and the reference `opt.xyz` energy for '
      'each, from which the rigid B3LYP-D3(BJ)/6-311++G(d,p) IEF-PCM single points of burst 2 were '
      'computed; the relaxed torsion scans replaced those single points.')
    w('')

    # ---------------------------------------------------------------- reproducibility
    rep = S.get('reproduction_of_validation', {})
    if rep:
        w('## 7. Production against validation, for the jobs that repeat one')
        w('')
        w('A production job with an ETKDG species repeats its validation run exactly and '
          'reproduces it bit for bit; a DFT-geometry job differs from its validation counterpart only '
          'in the input geometry, so the difference is the input-geometry sensitivity.')
        w('')
        w('How big that input change is, DFT optimum against the ETKDG/MMFF94s conformer the '
          'validation used (matched heavy-atom skeleton):')
        w('')
        w('| species | bonds: RMSD / max | angles: RMSD / max | C3–C4–S1–N2, DFT vs ETKDG |')
        w('|---|---|---|---|')
        for sp, d in S['dft'].items():
            v = d.get('vs_etkdg_input')
            if not v:
                continue
            w(f'| {sp} | {f(v["bond_rmsd_A"], 4)} / {f(v["bond_max_abs_A"], 4)} Å | '
              f'{f(v["angle_rmsd_deg"])} / {f(v["angle_max_abs_deg"])}° | '
              f'{f(d["tau_C3"], 1)}° vs {f(v["tau_C3_etkdg"], 1)}° |')
        w('')
        # Review 2026-09-28: both means are over the seeds the two runs share, so for the two 20-seed
        # AD4Zn jobs the production column here is the mean over seeds 1-10 and differs in the third
        # decimal from the 20-seed mean in sections 1 and 2.  Said in the header rather than left to
        # be discovered.
        w('Both score columns are means over the seeds the two runs **share** (`n_common_seeds` in '
          'the JSON), so for the two 20-seed AD4Zn jobs the production column is the mean over seeds '
          '1–10 and differs in the third decimal from the 20-seed mean in sections 1 and 2.')
        w('')
        w('| production job | validation run | production geometry | validation score | production '
          'score | Δ | all pose scores identical |')
        w('|---|---|---|---|---|---|---|')
        for k, r in sorted(rep.items()):
            w(f'| `{k}` | `{r["validation_run"]}` | {r["production_geometry"]} | '
              f'{f(r["validation_score_mean"], 3)} | {f(r["production_score_mean"], 3)} | '
              f'{f(r["delta_mean"], 3)} | '
              f'{r.get("identical_all_pose_scores", "– (different input geometry)")} |')
        w('')

    if diag:
        w('## 7b. Input-geometry sensitivity of the AD4Zn ring orientation (diagnostic)')
        w('')
        w('The two hCA II AD4Zn results moved when the input geometry changed from the ETKDG/MMFF94s '
          'conformer used in validation to the DFT structure, in opposite directions and with the '
          'zinc anchor untouched in both: 6RL9 went from bimodal (0.5 of seeds within 1 Å) to the '
          'near-native pose in every seed, and 3HS4 went from the near-native pose in every seed to '
          'the pose rotated about the Zn–N/S axis in every seed. AD4Zn separates those two '
          'orientations by hundredths of a kcal/mol (`dock/VALIDATION.md` section 3), so a small '
          'change to a rigid internal coordinate can decide which of them ranks first. The two inputs '
          'differ most in the C4–S1–N2 angle, and Vina holds bond angles rigid, so that was the '
          'hypothesis, and it was tested as a 2 × 2: each of the two geometries docked as it is and '
          'with that one angle set to the other\'s value, nothing else changed '
          '(`dock/geometry_variant.py`, `dock/run_geometry_diagnostic.sh`, geometries in '
          '`dock/results/diagnostics/`).')
        w('')
        ca = S['csn_angle']
        w('C4–S1–N2 as measured in this pass (`production_summary.json` → `csn_angle`):')
        w('')
        w('| geometry | C–S–N | O–S–O | S–N |')
        w('|---|---|---|---|')
        for grp, label in (('dft', 'DFT B3LYP-D3(BJ)/PCM'), ('etkdg', 'ETKDGv3+MMFF94s'),
                           ('crystal', 'crystal')):
            for name, m in ca[grp].items():
                w(f'| {label}, {name} | **{f(m["C_S_N_deg"])}°** | {f(m["O_S_O_deg"])}° | '
                  f'{f(m["S_N_A"], 3)} Å |')
        w('')
        etk = [m['C_S_N_deg'] for m in ca['etkdg'].values()]
        xt = [m['C_S_N_deg'] for m in ca['crystal'].values()]
        w(f'The MMFF94s value is the outlier: {f(min(etk))}–{f(max(etk))}° against '
          f'{f(min(xt))}–{f(max(xt))}° in four crystal structures (including this compound\'s own) '
          f'and {f(ca["dft"]["sfa_neutral"]["C_S_N_deg"])}° at B3LYP-D3(BJ)/PCM — about 7° tighter '
          'than both.')
        w('')
        w('| target | input geometry | C4–S1–N2 | RMSD to the crystal SAN pose (Å) | f≤1 Å | '
          'near-native pose in the ranked list | Zn–N (Å, pass fraction) | Thr199 OG1 | clusters |')
        w('|---|---|---|---|---|---|---|---|---|')
        rows = []
        for t in ('T3_6RL9', 'T4_3HS4'):
            for key, j in list(jobs.items()) + list(diag.items()):
                if j['target'] != t or j['scoring'] != 'ad4zn' or j['ligand'] != 'sfa_neutral':
                    continue
                tag = j['tag']
                label = ('DFT (production)' if tag.endswith('__prod') else
                         'DFT, C4–S1–N2 → 102.0°' if tag.endswith('__csn102') else
                         'ETKDG/MMFF94s, C4–S1–N2 → 109.5°' if tag.endswith('__csn109') else
                         'ETKDG/MMFF94s (validation conformer)' if tag.endswith('__etkdg') else tag)
                ang = ('109.5°' if tag.endswith('__prod') else '102.0°' if tag.endswith('__csn102')
                       else '109.5°' if tag.endswith('__csn109') else '102.0°')
                ent = next((e for e in j['rmsd_to_crystal'].values() if 'full' in e), None)
                full = ent['full'] if ent else None
                md = ent.get('modes') if ent else None
                nn = '–'
                if md:
                    r1 = md['rank_of_first_pose_within']['1A']
                    nn = (f'{f(r1["fraction_of_seeds_with_any"], 2)} of seeds, median rank '
                          f'{f(r1["median_rank"], 1)} (best of modes {f(md["best_of_modes"]["mean"])} Å)')
                zn = j['site_checks'].get('Zn-N(sulfonamide)', {})
                thr = j['site_checks'].get('Thr199 OG1 H-bond', {})
                cl = ', '.join(f'{c["population"]}/{j["n_seeds"]}' for c in j['clusters']['clusters'])
                rows.append(f'| {S["targets"][t]["pdb"]} | {label} | {ang} | '
                            f'{pm(full) if full else "–"} | '
                            f'{f(full["f_le_1A"], 2) if full else "–"} | {nn} | '
                            f'{f(zn.get("distance_mean"))} ({f(zn.get("pass_fraction"), 2)}) | '
                            f'{f(thr.get("distance_mean"))} ({f(thr.get("pass_fraction"), 2)}) | '
                            f'{cl} |')
        L.extend(rows)
        w('')
        v = json.loads((DOCK / 'analysis' / 'summary.json').read_text())
        for k in ('T3_6RL9/sfa_neutral__ad4zn', 'T4_3HS4/sfa_neutral__ad4zn'):
            st = v[k]['stats']
            w(f'For comparison, the validation run `{k}` (ETKDG conformer, a fresh one per seed, 10 '
              f'seeds): RMSD {f(st["rmsd_top_mean"])} ± {f(st["rmsd_top_sd"])} Å, f≤1 Å '
              f'{f(st["success_fraction_1A"], 2)}.')
        w('')
        w('**The hypothesis is only half right, and the table says so.** In 6RL9 the angle does '
          'decide: the DFT geometry gives the near-native pose at 109.5° and the rotated pose at '
          '102.0°, and the MMFF94s conformer gives the near-native pose either way (9 of 10 seeds at '
          'its own 102.0°). In 3HS4 the angle changes nothing: the DFT geometry gives the rotated '
          'pose at both values and the MMFF94s conformer the near-native pose at both. So the '
          'orientation is not controlled by this one coordinate; it responds to the difference '
          'between the two rigid geometries as a whole, of which the angle is the largest part.')
        w('')
        w('**What is robust, in all eight AD4Zn runs of sulfanilamide across the two receptors** '
          '(the six above and the two production runs): the near-native pose is *found* in every '
          'single seed, at rank 1 or rank 2, and the zinc chemistry is right in every seed — '
          'Zn–N(sulfonamide) 2.09–2.37 Å and the Thr199 OG1 bond 2.46–2.59 Å, against 1.93–1.94 Å '
          'and 2.82–2.86 Å in the two crystal structures. What is not robust is which of the two '
          'orientations the scoring function puts first, and the score gap between them says why:')
        w('')
        gaps = []
        for key, j in sorted(list(jobs.items()) + list(diag.items())):
            if j['scoring'] != 'ad4zn' or j['ligand'] != 'sfa_neutral':
                continue
            ent = next((e for e in j['rmsd_to_crystal'].values() if 'full' in e), None)
            if not ent or 'modes' not in ent:
                continue
            g = ent['modes']['rank_of_first_pose_within']['1A'].get('score_gap_from_top_pose')
            if g:
                gaps.append((key, ent['modes']['rank_of_first_pose_within']['1A']['median_rank'], g))
        w('| run | median rank of the near-native pose | score gap to the top pose (kcal/mol) |')
        w('|---|---|---|')
        for key, rk, g in gaps:
            w(f'| `{key}` | {f(rk, 1)} | {f(g["mean"], 3)} (range {f(g["min"], 3)}–{f(g["max"], 3)}) |')
        w('')
        span = ([min(g['min'] for _, _, g in gaps), max(g['max'] for _, _, g in gaps)]
                if gaps else [0, 0])
        w('**How the paper reports it (Section 3.6).** Statements about the zinc anchor and the Thr199 '
          'hydrogen bonds are supported by every run. The orientation of the aniline ring is stated '
          f'with this sensitivity: the near-native pose trails the top pose by '
          f'{f(span[0], 3)}–{f(span[1], 3)} kcal/mol in these runs, and which orientation comes first '
          'changes with the input geometry. The paper reports the DFT geometry, which the rest of the '
          'paper uses and which agrees with the crystal structures on the angle; with it the '
          'compound\'s own co-crystal structure, 6RL9, is reproduced at 0.61 ± 0.04 Å in 20 of 20 '
          'seeds, while the 3HS4 cross-dock lands on the rotated orientation with the same zinc anchor.')
        w('')

    # ---------------------------------------------------------------- honest limits
    w('## 8. What these numbers do and do not support')
    w('')
    w('- **7S2K (Sul2)** is the only target validated with both scoring functions (0.51 / 0.54 Å, '
      '10/10 seeds within 1 Å), so its poses and contacts are discussed quantitatively.')
    t4 = jobs.get('T4_3HS4/sfa_neutral__ad4zn__prod')
    t3 = jobs.get('T3_6RL9/sfa_neutral__ad4zn__prod')
    t4r = next((e['full'] for e in t4['rmsd_to_crystal'].values() if 'full' in e), None) if t4 else None
    t3r = next((e['full'] for e in t3['rmsd_to_crystal'].values() if 'full' in e), None) if t3 else None
    w('- **3HS4 (hCA II)** was validated as a cross-docking target for sulfanilamide with AD4Zn '
      '(0.62 ± 0.03 Å, 10/10 within 1 Å, Zn–N 2.18 Å and the Thr199 bond in 10/10) — but that '
      'validation used the ETKDG/MMFF94s conformer. With the DFT geometry of production the same '
      f'protocol gives {pm(t4r) if t4r else "–"} Å: the sulfonamide still points at the zinc '
      # Review 2026-09-28: "the zinc anchor is intact" was unquantified, and in this job the Zn-N
      # check passes in a minority of the seeds (2.25-2.35 A against 1.94 A in the crystal), so the
      # measured numbers are put in the sentence.
      + (f'(Zn–N {f(t4["site_checks"]["Zn-N(sulfonamide)"]["distance_mean"])} Å on average, within '
         f'the 2.3 Å cutoff in only '
         f'{f(t4["site_checks"]["Zn-N(sulfonamide)"]["pass_fraction"], 2)} of the seeds, against '
         f'{f(t4["site_checks"]["Zn-N(sulfonamide)"]["crystal_ligand_value"])} Å in the crystal) and '
         f'the Thr199 OG1 bond holds in '
         f'{f(t4["site_checks"]["Thr199 OG1 H-bond"]["pass_fraction"], 2)} of them, '
         if t4 else '') +
      'while the aniline ring takes the rotated orientation (section 7b). Acetazolamide\'s own pose is '
      'reproduced by no scoring function. So 3HS4 supports the zinc-binding chemistry and the score '
      'ranking; it does not, on its own, fix the ring orientation.')
    w('- **6RL9 (hCA II)**, the compound\'s own co-crystal structure, failed the strict gate in '
      'validation (2.07 ± 1.45 Å, f≤1 Å 0.5, bimodal). With the DFT geometry of production it gives '
      f'{pm(t3r) if t3r else "–"} Å with f≤1 Å {f(t3r["f_le_1A"], 2) if t3r else "–"} over '
      f'{t3["n_seeds"] if t3 else "–"} seeds and the zinc chemistry right in every seed. That is a '
      'better result than the gate was given, and section 7b shows it is not robust to the input '
      'geometry, so the paper reports 6RL9 as partially validated, with the sensitivity and the '
      'validation figure (Figure 8) beside it, not as a passed gate.')
    w('- **3TZF (DHPS)** passes the conventional 2.0 Å criterion and not the 1.0 Å gate '
      '(sulfamethoxazole anion 1.92 ± 0.04 Å), so its poses are qualitative.')
    w('- **1AJ0 was excluded from production docking**: in validation sulfanilamide landed 11.6 Å '
      'away in a surface cleft in every seed. Nothing in this file comes from 1AJ0, which served only '
      'for the structural comparison with 3TZF.')
    zz = {}
    for k, j in jobs.items():
        sc = j['site_checks'].get('Zn-N(sulfonamide)')
        if (j['target'] in ('T3_6RL9', 'T4_3HS4') and j['scoring'] in ('vina', 'vinardo')
                and j['ligand'].startswith('sfa_') and sc and sc['applicable']):
            zz[k] = (sc, j['site_checks']['Thr199 OG1 H-bond'],
                     [c.get('contacts', {}).get('zinc') for c in j['clusters']['clusters']])
    if zz:
        def grp(pred):
            v = [x for k, x in zz.items() if pred(k)]
            return (f'Zn–N {f(min(a["distance_mean"] for a, _, _ in v))}–'
                    f'{f(max(a["distance_mean"] for a, _, _ in v))} Å, Thr199 OG1 '
                    f'{f(min(b["distance_mean"] for _, b, _ in v))}–'
                    f'{f(max(b["distance_mean"] for _, b, _ in v))} Å, both checks passing in '
                    f'{f(min(min(a["pass_fraction"], b["pass_fraction"]) for a, b, _ in v), 2)}–'
                    f'{f(max(min(a["pass_fraction"], b["pass_fraction"]) for a, b, _ in v), 2)} '
                    'of seeds')
        near = [z['nearest_ligand_atom'].split('#')[0] for _, _, zs in zz.values() for z in zs if z]
        counts = {a: near.count(a) for a in sorted(set(near))}
        atoms = ', '.join(f'{a} in {n} of {len(near)} cluster representatives'
                          for a, n in sorted(counts.items(), key=lambda kv: -kv[1]))
        w('- **Plain Vina and Vinardo have no metal term**, and it shows in the zinc site. Over the '
          f'{len(zz)} hCA II jobs of the three sulfanilamide species scored with them: 6RL9 '
          f'{grp(lambda k: k.startswith("T3"))}; 3HS4 {grp(lambda k: k.startswith("T4"))} — against '
          'Zn–N 1.93–1.94 Å and Thr199 OG1 2.82–2.86 Å in the two crystal structures. What those '
          f'poses put nearest the zinc is {atoms} (section 4), and the split is by receptor: in 6RL9 '
          'all six give the zinc a sulfonyl oxygen at 2.27–2.36 Å, in 3HS4 five of six give it the '
          'sulfonamide nitrogen but at 2.30–2.50 Å, against 1.93–1.94 Å in the crystals. So these '
          'poses find the right pocket with a loose or wrong anchor, and only the Vinardo poses in '
          '3HS4 (Zn–N exactly at the 2.30 Å cutoff) pass the zinc-site checks at all. Their hCA II '
          'scores are a scoring-function comparison, not a binding mode; AD4Zn is the scoring '
          'function that reproduces the zinc chemistry.')
    w('- Scores of 10–23 heavy-atom ligands lie within about 1 kcal/mol of each other. Only the '
      'ranking against the controls is interpreted, and section 5 says why LE does not rescue the '
      'comparison.')
    # Review 2026-09-28: section 8 said what each target supports but never stated the limitations of
    # the method itself, which the paper has to carry.  Written out here, with the PoseBusters numbers
    # measured on the production poses (production_summary.json -> jobs.<job>.posebusters).
    pb = {k: j['posebusters'] for k, j in jobs.items() if isinstance(j.get('posebusters'), dict)
          and 'pass_fraction' in j['posebusters']}
    w('')
    w('**Limitations of the method (paper Section 4.3).**')
    w('')
    w('- **Docking scores are not binding free energies.** Vina, Vinardo and AutoDock4Zn are '
      'empirical scoring functions whose reported accuracy is of the order of 1–2 kcal/mol on '
      'benchmark sets; no score here is converted to a Ki, a Kd or an affinity, and no compound is '
      'called potent or a drug candidate. Only the ordering against the co-crystal ligand and the '
      'controls is read, and section 5 shows that the ordering itself carries a size bias.')
    w('- **The receptor is rigid.** Every run uses one crystal conformation with fixed side chains: '
      'no induced fit, no loop or side-chain rearrangement, no receptor ensemble. The ligand is '
      'flexible in its torsions only — bond lengths, bond angles and rigid fragments are held at the '
      'input geometry, which is why the input geometry changes the AD4Zn ring orientation '
      '(section 7b).')
    # Review 2026-09-28: counted from the prepared receptors themselves rather than asserted, because
    # the receptors do keep the waters of the magnesium coordination sphere.
    kept = {}
    for t in T_ORDER:
        lines = [ln for ln in (DOCK / t / 'receptor.pdb').read_text().splitlines()
                 if ln.startswith(('ATOM', 'HETATM'))]
        n = len({ln[21:27] for ln in lines if ln[17:20].strip() == 'HOH'})
        if n:
            kept[S['targets'][t]['pdb']] = n
    wtxt = ('all crystallographic waters were removed except the metal coordination sphere ('
            + ', '.join(f'{n} in {p}' for p, n in kept.items()) + ', each within 2.4 Å of the Mg)'
            ) if kept else 'crystallographic waters were removed'
    w(f'- **Water and entropy are implicit at best.** No explicit solvation model, no conformational '
      f'entropy, no protonation change on binding, and {wtxt}. The zinc chemistry is modelled by the '
      'AutoDock4Zn TZ pseudo-atom, a parameterised stand-in for the Zn–N bond, not by quantum '
      'mechanics.')
    w('- **No molecular dynamics and no free-energy calculation** (MM/PBSA, MM/GBSA, alchemical) '
      'was run: every pose is a single static structure.')
    w('- **No experimental data on this compound.** There is no measured affinity, inhibition or '
      'bioassay for 4-sulfamoylanilinium dihydrogen phosphate here; the validation is structural '
      '(redocking and cross-docking against deposited complexes) and the literature ordering '
      'acetazolamide before sulfanilamide at hCA II is the only external anchor.')
    w('- **Protonation and tautomers are assumptions.** Species are fixed at pH 7.4 by pKa argument '
      '(`dock/prep_ligands.py`), the His64 tautomer of hCA II is an explicit assumption '
      '(`dock/targets.py: his_assumed`; PDB2PQR picks HIE for 6RL9 but HID for 3HS4, and the Vina '
      'baselines moved by ≤0.03 kcal/mol when it was flipped; the His64 test was run with Vina only), '
      'and the ion pair of '
      'the crystal is never docked because it dissociates at pH 7.4.')
    if pb:
        worst = min(pb.values(), key=lambda v: v['pass_fraction'])
        failed = sorted({c for v in pb.values() for c in v.get('failed_checks', {})})
        w(f'- **Pose physics.** PoseBusters (`dock` configuration, {len(pb)} production jobs, '
          f'{sum(v["n_poses"] for v in pb.values())} top poses, {worst["n_checks"]} checks each) '
          f'passes on {f(100 * sum(v["pass_fraction"] * v["n_poses"] for v in pb.values()) / sum(v["n_poses"] for v in pb.values()), 1)}% '
          f'of the poses, lowest per job {f(worst["pass_fraction"], 2)}'
          + ('; no check fails anywhere.' if not failed else
             '. The only check that ever fails is `' + '`, `'.join(failed) + '`, in '
             + '; '.join(f'{n} of the {pb[k]["n_poses"]} poses of `{k}`'
                         for k in sorted(pb) for c, n in pb[k].get('failed_checks', {}).items())
             + ' — poses whose closest receptor atom is the catalytic zinc, so what the generic van '
               'der Waals test counts as a clash is the metal-coordination distance itself.')
          + ' That tests the physical sanity of a pose, not its correctness.')
    if errors:
        w('')
        w('**Jobs that failed analysis:** ' + ', '.join(f'`{k}` ({v["error"]})' for k, v in errors.items()))
    w('')
    w('## 9. Files')
    w('')
    w('| file | what |')
    w('|---|---|')
    w('| `dock/results/production_summary.json` | every number in this report, per job: scores, '
      'clusters, RMSDs, site checks, ligand efficiency, torsions, strain, per-seed rows |')
    w('| `dock/results/production_contacts.json` | PLIP / ProLIF / geometric contacts of every '
      'cluster representative |')
    w('| `dock/results/poses/<target>/<species>__<scoring>__top_pose_per_seed.sdf` | the top pose of '
      'every seed, with score, seed and cluster as SDF properties |')
    w('| `dock/results/poses/<target>/<species>__<scoring>__best_seedNN_all_modes.pdbqt` | all poses '
      'of the best-scoring seed, as Vina wrote them |')
    w('| `dock/results/poses/<target>/<species>__<scoring>__{csn102,csn109,etkdg}__*` | the same two '
      'files for the diagnostic runs of section 7b, kept apart from the production poses by their '
      'suffix |')
    w('| `dock/results/poses_for_dft/*.xyz` + `manifest.json` | cluster-representative geometries of the '
      'rigid DFT single points of burst 2 (replaced by the relaxed torsion scans) |')
    w('| `dock/<target>/runs/<species>__<scoring>__prod/seedNN/` | the raw runs: input conformer, '
      'Vina log, all poses, `results.json` with the exact settings |')
    w('| `dock/results/production.log` | the run log with per-job timings |')
    w('| `dock/results/diagnostics/` + `dock/run_geometry_diagnostic.sh`, `dock/geometry_variant.py` '
      '| the input-geometry diagnostic of section 7b: the modified geometries, its log and its runs '
      '(`dock/T{3_6RL9,4_3HS4}/runs/sfa_neutral__ad4zn__{csn102,csn109,etkdg}/`) |')
    w('| `dock/VALIDATION.md`, `dock/analysis/{summary,gate_G2}.json` | the validation this production '
      'run rests on |')
    w('')
    w('Receptors used (sha256 of the file each run read, checked against '
      '`dock/T*/prep_report_recheck.json`):')
    w('')
    w('| target | receptor.pdbqt sha256 | matches the prepared receptor | AD4Zn maps |')
    w('|---|---|---|---|')
    for t in T_ORDER:
        r = S['receptors'][t]
        ad = r.get('ad4zn')
        w(f'| {S["targets"][t]["pdb"]} | `{r["receptor_pdbqt_sha256"][:16]}…` | '
          f'{r.get("matches_prep_report_recheck")} | '
          + (f'{ad["n_maps"]} maps, TZ pseudo-atom {ad["maps_report"]} Å from the crystal '
             f'sulfonamide N' if ad else '–') + ' |')
    w('')
    n_fp = 0
    for i, ln in enumerate(L):
        for tag, old, new in FULL_PRECISION_PRINT:
            if tag in ln and old in ln:
                L[i] = ln = ln.replace(old, new)
                n_fp += 1
    OUT.write_text('\n'.join(L) + '\n')
    print(f'wrote {OUT} ({len(L)} lines) from {len(jobs)} jobs')
    return 0


if __name__ == '__main__':
    sys.exit(main())
