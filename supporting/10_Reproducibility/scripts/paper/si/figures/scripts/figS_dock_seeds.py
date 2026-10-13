"""SI figure (figS_dock_seeds): seed-to-seed spread of the docking results.

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_dock_seeds.py

(a-j) Top-pose score of every seed (dots) and the mean (bar), per target (rows) and scoring function (columns), for
      every species of the production matrix (dock/results/production_summary.json -> jobs[...].score.per_seed; the
      means are facts dock.prod.<T>.<species>.<scoring>.score).  Sulfanilamide species from their DFT geometries,
      controls from ETKDG/MMFF94s conformers.  Scores are empirical docking scores, not free energies of binding.
(k)   Top-pose RMSD of every seed from the crystal (or cross-docking) reference for the redocking jobs behind
      gate G2 (dock/analysis/summary.json -> per_seed.rmsd_top for the validation runs; production_summary.json ->
      rmsd_to_crystal.<ref>.full.per_seed for the sulfanilamide AD4Zn production runs), with the 1.0 A gate and the
      conventional 2.0 A criterion.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import compose, style  # noqa: E402

STEM = 'figS_dock_seeds'
PROD = json.loads((B.DOCK / 'results' / 'production_summary.json').read_text())
VAL = json.loads((B.DOCK / 'analysis' / 'summary.json').read_text())
TARGETS = [('T1_3TZF', '3TZF (DHPS)'), ('T5_7S2K', '7S2K (Sul2)'), ('T3_6RL9', '6RL9 (hCA II)'),
           ('T4_3HS4', '3HS4 (hCA II)')]
SCORINGS = [('vina', 'Vina'), ('vinardo', 'Vinardo'), ('ad4zn', 'AD4Zn')]
SPECIES = [('sfa_neutral', 'sulfanilamide'), ('sfa_anion', 'sulfonamidate'), ('sfa_cation', 'cation'),
           ('smx_neutral', 'SMX$^0$'), ('smx_anion', 'SMX$^-$'), ('azm_neutral', 'AZM$^0$'),
           ('azm_anion', 'AZM$^-$'), ('paba_anion', 'pABA$^-$'), ('dhp_anion', 'DHP$^-$')]
SFA = {'sfa_neutral', 'sfa_anion', 'sfa_cation'}
RMSD_JOBS = [('val', 'T5_7S2K/dhp_anion__vina', None, '7S2K, DHP$^-$, Vina'),
             ('val', 'T5_7S2K/dhp_anion__vinardo', None, '7S2K, DHP$^-$, Vinardo'),
             ('val', 'T1_3TZF/smx_anion__vina', None, '3TZF, SMX$^-$, Vina'),
             ('val', 'T3_6RL9/sfa_neutral__ad4zn', None, '6RL9, SAN, AD4Zn (ETKDG)'),
             ('prod', 'T3_6RL9/sfa_neutral__ad4zn__prod', 'SAN (6RL9, own)', '6RL9, SAN, AD4Zn (DFT)'),
             ('val', 'T4_3HS4/sfa_neutral__ad4zn', None, '3HS4, SAN cross, AD4Zn (ETKDG)'),
             ('prod', 'T4_3HS4/sfa_neutral__ad4zn__prod', 'SAN (6RL9, carried into 3HS4 by Ca superposition)',
              '3HS4, SAN cross, AD4Zn (DFT)')]


def main():
    style.apply()
    MM = style.MM
    W = style.DOUBLE
    head = compose._header_mm() * MM
    cell_h = 30 * MM
    left, gap_x = 19.5 * MM, 19.5 * MM
    cw = (W - left - 2 * gap_x - 1.5 * MM) / 3
    rmsd_h = 38 * MM
    H = len(TARGETS) * (head + cell_h + 8.5 * MM) + head + rmsd_h + 10 * MM
    fig = plt.figure(figsize=(W, H))
    letters = iter('abcdefghijklmnop')
    y = H
    for tg, tname in TARGETS:
        pdb = tg.split('_')[1]
        for c_, (sc, sname) in enumerate(SCORINGS):
            rows = []
            for sp, slab in SPECIES:
                key = f'{tg}/{sp}__{sc}__prod'
                if key in PROD['jobs']:
                    rows.append((sp, slab, PROD['jobs'][key]['score']['per_seed'],
                                 f'dock.prod.{pdb}.{sp}.{sc}.score'))
            if not rows:
                continue
            x0 = left + c_ * (cw + gap_x)
            ax = fig.add_axes([x0 / W, (y - head - cell_h) / H, cw / W, cell_h / H])
            compose._header(fig, (x0 - left + 0.6 * MM) / W, (y - 0.2 * MM) / H,
                            f'({next(letters)})', f'{tname}, {sname}')
            for k, (sp, slab, per, fkey) in enumerate(rows):
                yy = len(rows) - 1 - k
                col = B.FRAG_COLOURS['sulfanilamide'] if sp in SFA else '0.35'
                jit = (np.random.default_rng(k).random(len(per)) - 0.5) * 0.3
                ax.plot(per, yy + jit, 'o', ms=2.2, mfc=col, mec='none', alpha=0.8, zorder=3)
                mean = B.F(fkey)
                assert abs(mean - np.mean(per)) < 5e-3, (fkey, mean, np.mean(per))
                ax.plot([mean, mean], [yy - 0.32, yy + 0.32], color='k', lw=1.0, zorder=4)
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels([r[1] for r in rows][::-1], fontsize=style.TICK_PT)
            ax.set_ylim(-0.6, len(rows) - 0.4)
            ax.set_xlabel('Top-pose score (kcal mol$^{-1}$)', fontsize=style.TICK_PT, labelpad=1)
            ax.tick_params(axis='x', labelsize=style.TICK_PT)
            ax.grid(axis='x', color='0.9', lw=0.5, zorder=0)
            if c_ > 0:
                ax.tick_params(labelleft=True)
        y -= head + cell_h + 8.5 * MM
    # ---------------- RMSD per seed
    ax = fig.add_axes([48 * MM / W, (y - head - rmsd_h) / H, (W - 50 * MM) / W, rmsd_h / H])
    compose._header(fig, 0.6 * MM / W, (y - 0.2 * MM) / H, f'({next(letters)})',
                    'Top-pose RMSD of every seed from the crystal reference (redocking and cross-docking)')
    labels = []
    for k, (kind, job, ref, lab) in enumerate(RMSD_JOBS):
        if kind == 'val':
            r = [p['rmsd_top'] for p in VAL[job]['per_seed']]
        else:
            r = PROD['jobs'][job]['rmsd_to_crystal'][ref]['full']['per_seed']
        yy = len(RMSD_JOBS) - 1 - k
        jit = (np.random.default_rng(k).random(len(r)) - 0.5) * 0.3
        ax.plot(r, yy + jit, 'o', ms=2.4, mfc='#2B6CB0', mec='none', alpha=0.85, zorder=3)
        labels.append(lab)
    for xv, ls, lab in ((1.0, '-', '1.0 Å gate (G2)'), (2.0, (0, (3, 2)), '2.0 Å criterion')):
        ax.axvline(xv, color='0.25', lw=0.8, ls=ls, zorder=2)
        ax.text(xv + 0.04, len(RMSD_JOBS) - 0.55, lab, fontsize=style.SMALL_PT, ha='left', va='top', color='0.25')
    ax.set_yticks(range(len(RMSD_JOBS)))
    ax.set_yticklabels(labels[::-1], fontsize=style.TICK_PT)
    ax.set_ylim(-0.6, len(RMSD_JOBS) - 0.4)
    ax.set_xlim(0, 4.0)
    ax.set_xlabel('RMSD (Å)', fontsize=style.BASE_PT)
    B.save(fig, STEM, si=True, width_mm=W / MM)


if __name__ == '__main__':
    main()
