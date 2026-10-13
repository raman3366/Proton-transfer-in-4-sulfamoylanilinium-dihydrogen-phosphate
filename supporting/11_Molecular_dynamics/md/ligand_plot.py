"""Stage S1 of md/PLAN.md: the figure of the ligand check.

  /root/cqc-sfa-tools/md/bin/python md/ligand_plot.py        (after md/ligand_params.py and md/ligand_fit.py)

Reads the tables in md/ligand/ and writes md/ligand/torsion_profiles.png:
 top    relaxed scans about C3-C4-S-N2, DFT against GAFF2 as assigned and repaired (anion, neutral molecule);
 bottom MM against DFT for the single points on rigidly rotated geometries (rotation about S-N).
Colours: two categorical hues checked with the palette validator of the charting guide (blue = repaired, orange =
as assigned); the DFT reference is drawn in grey ink.  Text is never coloured.
"""
import csv, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = pathlib.Path(__file__).resolve().parent / 'ligand'
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE = '#2a78d6', '#eb6834'


def read(name):
    with open(OUT / name) as f:
        return list(csv.DictReader(f))


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc='left', fontsize=9.5, color=INK, pad=8)
    ax.set_xlabel(xlabel, fontsize=9, color=INK2)
    ax.set_ylabel(ylabel, fontsize=9, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=0)


def main():
    plt.rcParams['font.family'] = 'DejaVu Sans'
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    mk = dict(markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.3)
    for col, (sp, label) in enumerate((('anion', 'Sulfonamidate anion'), ('neutral', 'Sulfanilamide, neutral'))):
        ax = axes[0][col]
        b = read(f'scan_before_{sp}_obc2.csv'); a = read(f'scan_after_{sp}.csv'); v = read(f'scan_after_{sp}_vacuum.csv')
        x = [float(r['tau_C3_C4_S_N2']) for r in a]
        ax.plot([float(r['tau_C3_C4_S_N2']) for r in b], [float(r['E_mm_rel']) for r in b], '-o', color=ORANGE,
                linewidth=2, label='GAFF2 as assigned, implicit water', **mk)
        ax.plot(x, [float(r['E_mm_rel']) for r in v], '--', color=BLUE, linewidth=1.4,
                label='repaired, vacuum (not like for like)')
        ax.plot(x, [float(r['E_dft_rel']) for r in a], '-s', color=INK2, linewidth=2, markersize=8.5,
                markeredgecolor=SURFACE, markeredgewidth=1.3, label='DFT, relaxed scan (IEF-PCM water)')
        ax.plot(x, [float(r['E_mm_rel']) for r in a], '-o', color=BLUE, linewidth=2, markersize=5.5,
                markeredgecolor=SURFACE, markeredgewidth=1.0, label='repaired, implicit water')
        style(ax, f'{label}: rotation about C4–S, relaxed', 'C3–C4–S–N2 (degrees)',
              'energy above the lowest scan point (kcal/mol)')
        ax.set_xlim(86, 184); ax.set_xticks(range(90, 181, 15))
        ax.set_ylim(bottom=-0.05 * ax.get_ylim()[1])
        ax.legend(fontsize=7.5, frameon=False, loc='upper left' if sp == 'neutral' else 'upper right', labelcolor=INK2)
        # rigid single points
        ax = axes[1][col]
        rows = [r for r in read('rigid_single_points.csv') if r['species'] == sp and r['geometry'] != 'opt']
        lim = max(max(float(r['E_dft']) for r in rows), max(float(r['mm_repaired_OBC2']) for r in rows)) + 1.5
        lo = min(min(float(r['mm_gaff2_OBC2']) for r in rows), 0.0) - 1.0
        ax.plot([lo, lim], [lo, lim], color=INK2, linewidth=1.0)
        for key, colr, lab in (('mm_gaff2_OBC2', ORANGE, 'GAFF2 as assigned'), ('mm_repaired_OBC2', BLUE, 'repaired')):
            cl = [r for r in rows if abs(float(r['aniline_rotation'])) < 5]
            ro = [r for r in rows if abs(float(r['aniline_rotation'])) >= 5]
            ax.plot([float(r['E_dft']) for r in cl], [float(r[key]) for r in cl], 'o', color=colr,
                    label=f'{lab}; aniline as in the minimum', **mk)
            ax.plot([float(r['E_dft']) for r in ro], [float(r[key]) for r in ro], 'o', markerfacecolor=SURFACE,
                    markeredgecolor=colr, markeredgewidth=1.6, markersize=6,
                    label=f'{lab}; aniline also rotated')
        style(ax, f'{label}: single points on rigidly rotated geometries',
              'DFT energy above the minimum geometry (kcal/mol)', 'MM energy, same geometries (kcal/mol)')
        ax.set_xlim(lo, lim); ax.set_ylim(lo, lim)
        ax.legend(fontsize=7.5, frameon=False, loc='upper left', labelcolor=INK2)
        ax.annotate('MM = DFT', xy=(lim - 0.4, lim - 0.4), xytext=(-10, -22), textcoords='offset points',
                    ha='right', fontsize=7.5, color=INK2)
    fig.tight_layout(h_pad=2.2, w_pad=2.0)
    fig.savefig(OUT / 'torsion_profiles.png', facecolor=SURFACE)
    print('wrote', OUT / 'torsion_profiles.png')


if __name__ == '__main__':
    main()
