"""Shared drawing code of the two molecular-dynamics figures of the SI (figS_md_orientation, figS_md_rmsd).

    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_md_orientation.py
    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/figS_md_rmsd.py

Both figures have one layout: the nine runs with zinc set 2e in a 3 x 3 block (rows: start X, N, R; columns: runs 1
to 3) and, set apart below, the two sensitivity runs of start R, which are in no count of rule 3 of md/PLAN.md.
Every panel has the same axes.  The 1.3 ns of the equilibration after the release of the ligand (stages 4, 5 and
6) are drawn at negative time on a stretched strip at the left of each panel, because every conversion of the
campaign lies there; the production runs from 0 to 100 ns.

Data: md/results/series/<set>_<run>.csv.gz (md/production/make_series.py; one line per 10 ps, written from the
output of the first analysis program and compared with the second).  Nothing is computed from coordinates.  In
the production the pale line is every frame and the strong line the mean over 100 frames (1 ns); in the strip
every frame is drawn.  The ends of stages 4 and 5 are read from the column `part` of the series.  Verdicts,
midpoints, the 2.0 A of rule 1 and the times of the plan's rule are facts of section 40 (facts.json, md.*) and
are written to the sidecar <stem>.facts.json; the script stops if a verdict or a time of the rule formed here
from the series differs from its fact.  The design follows md/results/fig_orientation.png
(md/production/make_result_figures.py), redrawn to the house rules (figkit: 175 mm wide, Liberation Sans, 600 dpi).
Colours: blue and orange of that figure, a pair checked there for colour-vision deficiencies.
"""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'main' / 'figures' / 'scripts'))
import _figb as B  # noqa: E402
from figkit import compose, style  # noqa: E402

SER = B.REPO / 'md' / 'results' / 'series'
BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, STRIP = '#0b0b0b', '#52514e', '#ecebe6'
SMOOTH = 100                      # frames of the running mean in the production (1 ns)
NINE = [['X1', 'X2', 'X3'], ['N1', 'N2', 'N3'], ['R1', 'R2', 'R3']]
ROW_LABEL = {'X': 'Start X: crystal pose', 'N': 'Start N: docked, near-native', 'R': 'Start R: docked, ring-rotated'}
SENS = [('R1_2s', 'bonded_2s_R1', 'R1, zinc set 2s'), ('R1_n2fit', 'bonded_2e_n2fit_R1', 'R1, candidate N2 terms')]
OPP = {'C': 'R', 'R': 'C'}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(stem):
    """One series: time (ns), d1, d2, RMSD, state and part of every frame."""
    with gzip.open(SER / f'{stem}.csv.gz', 'rt') as fh:
        rows = list(csv.DictReader(fh))
    part = np.array([r['part'] for r in rows])
    d = dict(t=np.array([float(r['time_ns']) for r in rows]), d1=np.array([float(r['d1_A']) for r in rows]),
             d2=np.array([float(r['d2_A']) for r in rows]), rmsd=np.array([float(r['lig_rmsd_A']) for r in rows]),
             state=[r['state'] for r in rows], part=part)
    d['eq'] = part != 'prod'
    d['stage_ends'] = [float(d['t'][part == p].max()) for p in ('eq4', 'eq5')]
    return d


def verdict_and_rule(d, start_state):
    """The verdict over the last 20 ns and the time of the plan's rule (ps from the release), formed here from the
    states of the series; compared with the facts before anything is drawn."""
    st = d['state']
    w = st[-2000:]
    opp = OPP[start_state]
    verdict = 'CONVERTED' if 10 * w.count(opp) >= 9 * len(w) else 'STAYED' if 10 * w.count(start_state) >= 9 * len(w) else 'MIXED'
    assigned, prev = [], start_state
    for s in st:
        prev = s if s in 'CR' else prev
        assigned.append(prev)
    k, rule = 0, None
    while k < len(assigned):
        j = k
        while j < len(assigned) and assigned[j] == assigned[k]:
            j += 1
        if assigned[k] == opp and j - k >= 100:
            rule = 10 * (k + 1)
            break
        k = j
    return verdict, rule


def running_mean(t, y, n=SMOOTH):
    k = np.ones(n) / n
    return np.convolve(t, k, mode='valid'), np.convolve(y, k, mode='valid')


def panel(fig, W, H, x0, ytop, pw, ah, d, kind, letter, title, rid, show_x, show_y, ylim, yticks):
    """One panel with its header; x0, ytop, pw, ah in inches (left edge, top of the header, width, axes height)."""
    MM = style.MM
    head = 4.6 * MM
    sw = pw / 4.3 - 0.3 * MM
    y0 = ytop - head - ah
    ax0 = fig.add_axes([x0 / W, y0 / H, sw / W, ah / H])
    ax1 = fig.add_axes([(x0 + sw + 0.6 * MM) / W, y0 / H, (pw - sw - 0.6 * MM) / W, ah / H])
    eq, pr = d['eq'], ~d['eq']
    for key, col in ([('d1', BLUE), ('d2', ORANGE)] if kind == 'orientation' else [('rmsd', BLUE)]):
        ax0.plot(d['t'][eq], d[key][eq], color=col, lw=0.7, solid_joinstyle='round')
        ax1.plot(d['t'][pr], d[key][pr], color=col, lw=0.3, alpha=0.33, rasterized=True)
        tm, ym = running_mean(d['t'][pr], d[key][pr])
        ax1.plot(tm, ym, color=col, lw=0.95, solid_joinstyle='round')
    for ax in (ax0, ax1):
        if kind == 'orientation':
            ax.axhline(B.F('md.def.mid_d1'), color=INK, lw=0.6, ls=(0, (4, 2)), zorder=4)
            ax.axhline(B.F('md.def.mid_d2'), color=INK, lw=0.6, ls=(0, (1, 1.6)), zorder=4)
        else:
            ax.axhline(B.F('md.def.rule1'), color=INK, lw=0.6, ls=(0, (4, 2)), zorder=4)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.tick_params(axis='both', labelsize=style.TICK_PT, length=2.2, width=0.5, colors=INK2, labelcolor=INK2, pad=1.5)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'):
            ax.spines[s].set_linewidth(0.5)
            ax.spines[s].set_color(INK2)
    released = B.F('md.prot.eq_released') / 1000.0                 # ns of the equilibration after the release
    ax0.set_facecolor(STRIP)
    ax0.set_xlim(-released, 0.0)
    ax0.set_xticks([-released, 0.0])
    ax1.set_xlim(0.0, B.F('md.prot.prod.length'))
    ax1.set_xticks([0, 25, 50, 75, 100])
    if show_x:
        ax0.set_xticklabels([B.num(-released, 1), ''])
        ax1.set_xticklabels(['0', '25', '50', '75', '100'])
        ax0.get_xticklabels()[0].set_ha('left')
        ax1.get_xticklabels()[-1].set_ha('right')
    else:
        ax0.set_xticklabels(['', ''])
        ax1.set_xticklabels([])
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(axis='y', left=False, labelleft=False)
    if not show_y:
        ax0.set_yticklabels([])
    for t_end in d['stage_ends']:
        ax0.axvline(t_end, color=INK2, lw=0.45, ls=(0, (2, 2)), alpha=0.8, zorder=3)
    # the verdict and the time of the rule: facts, checked against the series
    start_state = 'C' if rid[0] in 'XN' else 'R'
    verdict, rule = verdict_and_rule(d, start_state)
    assert verdict == B.F(f'md.run.{rid}.verdict'), (rid, verdict)
    fact_rule = B.F(f'md.run.{rid}.rule.ps') if f'md.run.{rid}.rule.ps' in B.facts() else None
    assert rule == fact_rule, (rid, rule, fact_rule)
    if rule is not None:
        tr = -released + rule / 1000.0
        ytop_m = ylim[0] + 0.945 * (ylim[1] - ylim[0])
        ax0.plot([tr, tr], [ylim[0], ytop_m], color=INK, lw=0.5, alpha=0.55, zorder=5)
        if verdict == 'CONVERTED':
            ax0.plot([tr], [ytop_m], marker='v', ms=4.0, color=INK, mec='white', mew=0.4, zorder=6, clip_on=False)
        else:                       # the rule gives a time, the verdict of the run is MIXED: no conversion
            ax0.plot([tr], [ytop_m], marker='v', ms=4.0, mfc='white', mec=INK, mew=0.7, zorder=6, clip_on=False)
    compose._header(fig, x0 / W, (ytop - 0.2 * MM) / H, f'({letter})', title)
    fig.text((x0 + pw) / W, (ytop - 0.9 * MM) / H, verdict, fontsize=style.TICK_PT, color=INK2, ha='right', va='top')
    return ax0, ax1


def make(kind, stem):
    style.apply()
    MM = style.MM
    W = style.DOUBLE - 1.0 * MM                # the saved file adds 0.5 mm on each side: 175 mm at most
    left, right, gap_x = 15.5 * MM, 0.8 * MM, 4.0 * MM
    pw = (W - left - right - 2 * gap_x) / 3
    head, ah, gap_y = 4.6 * MM, 33.0 * MM, 2.6 * MM
    ticks_h, xlab_h, sep_h, top, bottom = 4.2 * MM, 4.4 * MM, 13.0 * MM, 0.6 * MM, 0.8 * MM
    xlabel = 'Time (ns); 0 = start of the production'
    H = top + 3 * (head + ah) + 2 * gap_y + ticks_h + xlab_h + sep_h + head + ah + ticks_h + xlab_h + bottom
    fig = plt.figure(figsize=(W, H))
    ylim, yticks = ((2.0, 12.0), [2, 4, 6, 8, 10, 12]) if kind == 'orientation' else ((0.0, 4.5), [0, 1, 2, 3, 4])
    ylabel = 'Distance (Å)' if kind == 'orientation' else 'Ligand RMSD (Å)'
    letters = iter('abcdefghijk')
    series, outside = {}, {}
    y = H - top
    for i, row in enumerate(NINE):
        for j, rid in enumerate(row):
            d = series[rid] = load(f'bonded_2e_{rid}')
            x0 = left + j * (pw + gap_x)
            ax0, _ = panel(fig, W, H, x0, y, pw, ah, d, kind, next(letters), rid, rid, show_x=(i == 2), show_y=(j == 0), ylim=ylim, yticks=yticks)
            if j == 0:
                ax0.set_ylabel(ylabel, fontsize=style.BASE_PT, color=INK2, labelpad=1.5)
                fig.text(0.3 * MM / W, (y - head - ah / 2) / H, ROW_LABEL[rid[0]], rotation=90, va='center', ha='left', fontsize=style.BASE_PT, color=INK)
        y -= head + ah + gap_y
    y += gap_y
    y -= ticks_h
    fig.text((left + (W - left - right) / 2) / W, (y - 0.4 * MM) / H, xlabel, fontsize=style.BASE_PT, color=INK2, ha='center', va='top')
    y -= xlab_h + 3.2 * MM
    y_line = y
    fig.add_artist(Line2D([0.3 * MM / W, 1 - right / W], [y / H, y / H], color=INK2, lw=0.5))
    fig.text(0.3 * MM / W, (y - 1.2 * MM) / H, 'Sensitivity runs of start R, one run each (not among the nine)', fontsize=style.BASE_PT, fontweight='bold',
             color=INK, va='top', ha='left')
    y -= sep_h - 3.2 * MM
    for j, (rid, stem_s, title) in enumerate(SENS):
        d = series[rid] = load(stem_s)
        x0 = left + j * (pw + gap_x)
        ax0, ax1 = panel(fig, W, H, x0, y, pw, ah, d, kind, next(letters), title, rid, show_x=True, show_y=(j == 0), ylim=ylim, yticks=yticks)
        if j == 0:
            ax0.set_ylabel(ylabel, fontsize=style.BASE_PT, color=INK2, labelpad=1.5)
            fig.text(0.3 * MM / W, (y - head - ah / 2) / H, 'Start R', rotation=90, va='center', ha='left', fontsize=style.BASE_PT, color=INK)
    fig.text((left + pw + gap_x / 2) / W, (y - head - ah - ticks_h - 0.4 * MM) / H, xlabel, fontsize=style.BASE_PT, color=INK2, ha='center', va='top')
    for rid, d in series.items():
        outside[rid] = (int(((d['d1'] < ylim[0]) | (d['d1'] > ylim[1]) | (d['d2'] < ylim[0]) | (d['d2'] > ylim[1])).sum()) if kind == 'orientation'
                        else int((d['rmsd'] > ylim[1]).sum()))
    # legend, in the free cell of the lower row
    x0 = left + 2 * (pw + gap_x)
    y_bot = y - head - ah - ticks_h - xlab_h
    lax = fig.add_axes([x0 / W, y_bot / H, pw / W, (y_line - 2.0 * MM - y_bot) / H])
    lax.axis('off')
    e4, e5 = round(B.F('md.prot.eq_stage4_end')), round(B.F('md.prot.eq_stage5_end'))
    assert all(abs((d['stage_ends'][0] + B.F('md.prot.eq_released') / 1000) * 1000 - e4) < 1e-6 and
               abs((d['stage_ends'][1] + B.F('md.prot.eq_released') / 1000) * 1000 - e5) < 1e-6 for d in series.values())
    if kind == 'orientation':
        handles = [Line2D([], [], color=BLUE, lw=1.5), Line2D([], [], color=ORANGE, lw=1.5),
                   Line2D([], [], color=INK, lw=0.7, ls=(0, (4, 2))), Line2D([], [], color=INK, lw=0.7, ls=(0, (1, 1.6)))]
        labels = ['d1: aniline N1 to Phe131 CZ', 'd2: aniline N1 to Asn67 OD1', f"midpoint of d1, {B.fmt('md.def.mid_d1')} Å", f"midpoint of d2, {B.fmt('md.def.mid_d2')} Å"]
    else:
        handles = [Line2D([], [], color=BLUE, lw=1.5), Line2D([], [], color=INK, lw=0.7, ls=(0, (4, 2)))]
        labels = ['ligand heavy atoms against the\ncrystal pose (6RL9), after the\nfit of the Cα atoms',
                  f"{B.fmt('md.def.rule1')} Å (rule 1 reads the mean of the\nlast half of a crystal-start run)"]
    handles += [Line2D([], [], color=INK2, lw=0.5, alpha=0.5), Line2D([], [], color=INK2, lw=1.2), Patch(facecolor=STRIP, edgecolor='none'),
                Line2D([], [], color=INK2, lw=0.6, ls=(0, (2, 2))),
                Line2D([], [], marker='v', ms=4.0, color=INK, mec='white', mew=0.4, ls='none'),
                Line2D([], [], marker='v', ms=4.0, mfc='white', mec=INK, mew=0.7, ls='none')]
    labels += [f"pale: every frame ({B.fmt('md.prot.out.solute')} ps)", f"strong: mean over {B.fmt('md.def.dwell')} ns (production)",
               f"grey strip: equilibration after the release\nof the ligand ({B.num(-B.F('md.prot.eq_released') / 1000, 1)} to 0 ns, stretched)",
               f"ends of stages 4 and 5 ({e4} and {e5} ps\nafter the release; backbone held in\nstages 4 and 5, free in stage 6)",
               'conversion by the rule fixed beforehand\n(start of the first visit of at least 1 ns\nin the crystal-like state)',
               'the same frame in the run whose\nverdict is MIXED: not a conversion']
    leg = lax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=style.TICK_PT, frameon=False, handlelength=1.5, handletextpad=0.6, labelspacing=0.42,
                     borderaxespad=0.0)
    for txt in leg.get_texts():
        txt.set_color(INK)
    extra = {'series': {rid: {'file': str((SER / f'{s}.csv.gz').relative_to(B.REPO)), 'sha256': sha256_file(SER / f'{s}.csv.gz')}
                        for rid, s in [(r, f'bonded_2e_{r}') for row in NINE for r in row] + [(r, s) for r, s, _ in SENS]},
             'running_mean_frames': SMOOTH, 'y_axis': list(ylim), 'frames_with_a_value_outside_the_y_axis': outside,
             'canvas_mm': [round(W / MM, 1), round(H / MM, 1)], 'colours': {'d1 and RMSD': BLUE, 'd2': ORANGE},
             'conversion_marker_open_because_the_verdict_is_MIXED': [rid for rid in series if B.F(f'md.run.{rid}.verdict') == 'MIXED']}
    B.save(fig, stem, si=True, width_mm=style.DOUBLE / MM - 1.0, extra=extra)       # QA allows 1 mm over: 175 mm
    # the size of the saved PNG, read back from the file: the white border is trimmed on saving, so it is not the canvas
    from PIL import Image
    side = B.FIG_SI / f'{stem}.facts.json'
    with Image.open(B.FIG_SI / f'{stem}.png') as im:
        px, dpi = im.size, float(im.info['dpi'][0])
    rec = json.loads(side.read_text(encoding='utf-8'))
    rec['png_px'], rec['png_dpi'], rec['png_mm'] = list(px), round(dpi), [round(px[0] / dpi * 25.4, 1), round(px[1] / dpi * 25.4, 1)]
    side.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str) + '\n', encoding='utf-8')


if __name__ == '__main__':
    print('shared code of figS_md_orientation.py and figS_md_rmsd.py; run one of those')
