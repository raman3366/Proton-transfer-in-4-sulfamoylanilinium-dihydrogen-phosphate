#!/usr/bin/env python3
"""The two figures of md/RESULTS.md, drawn from the time series files of the analysis.

    md/results/fig_orientation.png, .pdf    d1 and d2 against time, eleven panels
    md/results/fig_ligand_rmsd.png, .pdf    ligand heavy-atom RMSD to the crystal pose, same layout
    md/results/fig_sources.json             the files read, with checksums, and what was checked

Layout: the nine production runs in a 3 x 3 block (rows: start X, N, R; columns: runs 1 to 3); below them, set
apart, the two sensitivity runs of start R, which are in no count of rule 3 of the plan. Every panel has the
same axes. The part of the equilibration after the release of the ligand (stages 4, 5 and 6; 1.3 ns) is drawn
at negative time on a stretched strip at the left of each panel, because every conversion of the campaign lies
there; the production runs from 0 to 100 ns.

What is drawn: the columns `d1_N1_F131CZ_A`, `d2_N1_N67OD1_A` and `lig_rmsd_A` of
`<run>/analysis/timeseries_equil_solute.csv` (stages 4 to 6) and `<run>/analysis/timeseries_solute.csv`, the
output of the first program (md/analysis/analyse_md.py). Before drawing, the same columns of the second
program (`<run>/analysis_recheck/`) are read and must agree within 0.001 A in every frame (the tolerance of
section 14 of md/analysis/SPEC.md); the script stops if they do not. Nothing is computed from coordinates.
In the production the pale line is every frame (10 ps) and the strong line the mean over 1 ns (100 frames,
centred); in the strip every frame is drawn. The midpoints 6.7 A and 5.2 A are those of rule 2 of the plan.
The marker of a conversion is read from md/results/conversion_times.json (md/production/conversion_time.py): a
filled triangle in a run whose verdict is CONVERTED; an open triangle in a run whose verdict is MIXED, where the
plan's rule gives the start of the first visit of at least 1 ns in the opposite state and that is no conversion.
Two thin lines in the strip mark the ends of equilibration stages 4 and 5, read from the column `stage` of the
equilibration series (the protein backbone is held in stages 4 and 5 and free in stage 6).

    /root/cqc-sfa-tools/md/bin/python md/production/make_result_figures.py [--runs DIR] [--out DIR]
"""
import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MD = os.path.dirname(HERE)

# colours: slots 1 and 2 of a validated categorical palette (blue, orange); checked for colour-vision
# deficiencies with the palette validator (worst pair: protan, OKLab difference 24.7; at least 8 is asked)
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, STRIP = "#0b0b0b", "#52514e", "#ecebe6"
M1, M2 = 6.7, 5.2                 # midpoints of rule 2 (md/PLAN.md, section 5; md/analysis/SPEC.md, section 4)
RULE1 = 2.0                       # rule 1: mean ligand RMSD of the last half above 2.0 A
SMOOTH = 100                      # frames of the running mean in the production (1 ns)
TOL = 0.001                       # A; agreement of the two programs (SPEC.md, section 14)

NINE = [["X1", "X2", "X3"], ["N1", "N2", "N3"], ["R1", "R2", "R3"]]
ROW_LABEL = {"X": "start X: crystal pose", "N": "start N: docked, near-native", "R": "start R: docked, ring rotated"}
SENS = [("bonded_2s/R1", "R1, zinc set 2s"), ("bonded_2e_n2fit/R1", "R1, candidate N2 terms")]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_cols(path, names, released_only=False):
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        idx = [header.index(n) for n in names]
        st = header.index("stage") if released_only else None
        ex = header.index("extension") if released_only else None
        rows = []
        for row in rd:
            if not row:
                continue
            if released_only and not (int(row[st]) >= 4 and int(row[ex]) == 0):
                continue
            rows.append([row[i] for i in idx])
    return rows


def load_run(run_dir, sources):
    """Released equilibration and production of one run: time (ns), d1, d2, RMSD, state; first program,
    checked against the second."""
    cols = ["time_ns", "d1_N1_F131CZ_A", "d2_N1_N67OD1_A", "lig_rmsd_A", "state"]
    out = {}
    # ends of stages 4 and 5 on the time axis of the strip: the time of the last frame of each
    st_rows = read_cols(os.path.join(run_dir, "analysis", "timeseries_equil_solute.csv"), ["time_ns", "stage"], released_only=True)
    out["stage_ends_ns"] = [max(float(t) for t, st in st_rows if int(st) == k) for k in (4, 5)]
    for part, fname, rel in (("equil", "timeseries_equil_solute.csv", True), ("prod", "timeseries_solute.csv", False)):
        data = {}
        for sub in ("analysis", "analysis_recheck"):
            p = os.path.join(run_dir, sub, fname)
            rows = read_cols(p, cols, released_only=rel)
            data[sub] = (np.array([[float(x) for x in r[:4]] for r in rows]), [r[4] for r in rows])
            sources.append({"file": p, "sha256": sha256_file(p), "rows_used": len(rows)})
        a, b = data["analysis"], data["analysis_recheck"]
        if a[0].shape != b[0].shape:
            raise SystemExit("%s %s: the two programs wrote %s and %s rows" % (run_dir, fname, a[0].shape, b[0].shape))
        worst = float(np.max(np.abs(a[0] - b[0]))) if a[0].size else 0.0
        if worst > TOL or a[1] != b[1]:
            raise SystemExit("%s %s: the two programs differ in the columns drawn (largest difference %.4f)"
                             % (run_dir, fname, worst))
        out[part] = {"t": a[0][:, 0], "d1": a[0][:, 1], "d2": a[0][:, 2], "rmsd": a[0][:, 3], "state": a[1],
                     "max_abs_diff_between_programs": worst}
    with open(os.path.join(run_dir, "analysis", "summary.json")) as fh:
        summ = json.load(fh)
    out["verdict"] = summ["orientation"]["verdict"]
    out["start"] = summ["run"]["start"]
    return out


def running_mean(t, y, n=SMOOTH):
    if len(y) < n:
        return t[:0], y[:0]
    k = np.ones(n) / n
    return np.convolve(t, k, mode="valid"), np.convolve(y, k, mode="valid")


def conversion_marker(conv):
    """(part, time in ns on the axis of that part) of the conversion of a run, or None."""
    if not conv:
        return None
    return ("equil" if conv["in_equilibration"] else "prod", conv["time_ns_in_file"])


def label_time_axis(ax0, ax1, show):
    """Tick labels of the two parts of the time axis; the outer ones are turned inwards so that the labels
    of neighbouring panels do not meet."""
    if not show:
        ax0.set_xticklabels(["", ""])
        ax1.set_xticklabels([])
        return
    ax0.set_xticklabels(["−1.3", ""])
    ax1.set_xticklabels(["0", "25", "50", "75", "100"])
    ax0.get_xticklabels()[0].set_ha("left")
    ax1.get_xticklabels()[-1].set_ha("right")


def draw_panel(fig, cell, run, kind, title, conv, ylim, yticks, show_x, show_y):
    inner = GridSpecFromSubplotSpec(1, 2, subplot_spec=cell, width_ratios=[1.0, 3.3], wspace=0.035)
    ax0 = fig.add_subplot(inner[0])
    ax1 = fig.add_subplot(inner[1])
    eq, pr = run["equil"], run["prod"]
    series = [("d1", BLUE), ("d2", ORANGE)] if kind == "orientation" else [("rmsd", BLUE)]
    for key, col in series:
        ax0.plot(eq["t"], eq[key], color=col, lw=0.75, solid_joinstyle="round")
        ax1.plot(pr["t"], pr[key], color=col, lw=0.35, alpha=0.33, rasterized=True)
        tm, ym = running_mean(pr["t"], pr[key])
        ax1.plot(tm, ym, color=col, lw=1.05, solid_joinstyle="round")
    for ax in (ax0, ax1):
        if kind == "orientation":
            ax.axhline(M1, color=INK, lw=0.7, ls=(0, (4, 2)), zorder=4)
            ax.axhline(M2, color=INK, lw=0.7, ls=(0, (1, 1.6)), zorder=4)
        else:
            ax.axhline(RULE1, color=INK, lw=0.7, ls=(0, (4, 2)), zorder=4)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.tick_params(axis="both", labelsize=6.5, length=2.5, width=0.5, colors=INK2, labelcolor=INK2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_linewidth(0.5)
            ax.spines[s].set_color(INK2)
    ax0.set_facecolor(STRIP)
    ax0.set_xlim(-1.3, 0.0)
    ax0.set_xticks([-1.3, 0.0])
    ax1.set_xlim(0.0, 100.0)
    ax1.set_xticks([0, 25, 50, 75, 100])
    label_time_axis(ax0, ax1, show_x)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", left=False, labelleft=False)
    if not show_y:
        ax0.set_yticklabels([])
    for t_end in run["stage_ends_ns"]:
        ax0.axvline(t_end, color=INK2, lw=0.45, ls=(0, (2, 2)), alpha=0.8, zorder=3)
    mk = conversion_marker(conv)
    if mk:
        ax = ax0 if mk[0] == "equil" else ax1
        ytop = ylim[0] + 0.955 * (ylim[1] - ylim[0])
        ax.plot([mk[1], mk[1]], [ylim[0], ytop], color=INK, lw=0.5, alpha=0.55, zorder=5)
        if run["verdict"] == "CONVERTED":
            ax.plot([mk[1]], [ytop], marker="v", ms=4.5, color=INK, mec="white", mew=0.4, zorder=6)
        else:                       # the rule gives a time, the verdict of the run is not CONVERTED
            ax.plot([mk[1]], [ytop], marker="v", ms=4.5, mfc="white", mec=INK, mew=0.7, zorder=6)
    ax0.set_title(title[0], loc="left", fontsize=7.5, fontweight="bold", color=INK, pad=3)
    ax1.set_title(title[1], loc="right", fontsize=6.5, color=INK2, pad=3)
    return ax0, ax1


def make_figure(kind, runs, conv, out_base):
    ylim, yticks = ((2.0, 12.0), [2, 4, 6, 8, 10, 12]) if kind == "orientation" else ((0.0, 4.5), [0, 1, 2, 3, 4])
    fig = plt.figure(figsize=(7.2, 8.3))
    outer = GridSpec(5, 3, figure=fig, height_ratios=[1, 1, 1, 0.03, 1], hspace=0.34, wspace=0.10,
                     left=0.105, right=0.985, top=0.968, bottom=0.060)
    ylabel = "distance (Å)" if kind == "orientation" else "ligand RMSD (Å)"
    for i, row in enumerate(NINE):
        for j, rid in enumerate(row):
            label = "bonded_2e/" + rid
            ax0, ax1 = draw_panel(fig, outer[i, j], runs[label], kind, (rid, runs[label]["verdict"]),
                                  conv.get(label), ylim, yticks, show_x=(i == 2), show_y=(j == 0))
            if j == 0:
                ax0.set_ylabel(ylabel, fontsize=7, color=INK2, labelpad=2)
                bb = ax0.get_position()
                fig.text(0.012, 0.5 * (bb.y0 + bb.y1), ROW_LABEL[rid[0]], rotation=90, va="center", ha="left",
                         fontsize=7.5, color=INK)
    # the two sensitivity runs, set apart
    y_rule = outer[2, 0].get_position(fig).y0 - 0.042
    fig.add_artist(Line2D([0.012, 0.985], [y_rule, y_rule], color=INK2, lw=0.5))
    fig.text(0.012, y_rule - 0.009, "Sensitivity runs of start R, one run each (not among the nine; in no count of rule 3)",
             fontsize=7.5, fontweight="bold", color=INK, va="top", ha="left")
    for j, (label, name) in enumerate(SENS):
        ax0, ax1 = draw_panel(fig, outer[4, j], runs[label], kind, (name, runs[label]["verdict"]),
                              conv.get(label), ylim, yticks, show_x=True, show_y=(j == 0))
        if j == 0:
            ax0.set_ylabel(ylabel, fontsize=7, color=INK2, labelpad=2)
            bb = ax0.get_position()
            fig.text(0.012, 0.5 * (bb.y0 + bb.y1), "start R", rotation=90, va="center", ha="left", fontsize=7.5, color=INK)
        ax1.set_xlabel("time (ns); 0 = start of the production", fontsize=7, color=INK2, labelpad=2)
    # (the time axis is the same in every panel; its labels stand under the bottom row of each block)
    # legend, in the free cell of the lower row
    lax = fig.add_subplot(outer[4, 2])
    lax.axis("off")
    if kind == "orientation":
        handles = [Line2D([], [], color=BLUE, lw=1.6), Line2D([], [], color=ORANGE, lw=1.6),
                   Line2D([], [], color=INK, lw=0.8, ls=(0, (4, 2))), Line2D([], [], color=INK, lw=0.8, ls=(0, (1, 1.6)))]
        labels = ["d1: aniline N1 to Phe131 CZ", "d2: aniline N1 to Asn67 OD1",
                  "midpoint of d1, 6.7 Å", "midpoint of d2, 5.2 Å"]
    else:
        handles = [Line2D([], [], color=BLUE, lw=1.6), Line2D([], [], color=INK, lw=0.8, ls=(0, (4, 2)))]
        labels = ["ligand heavy atoms against the\ncrystal pose (6RL9), after the\nfit of the Cα atoms",
                  "2.0 Å (rule 1 reads the mean of the\nlast half of a crystal-start run)"]
    ends = sorted({round((t + 1.3) * 1000.0) for lab in runs for t in runs[lab]["stage_ends_ns"]})
    assert len(ends) == 2, ends                      # the same two times in every run
    handles += [Line2D([], [], color=INK2, lw=0.5, alpha=0.5), Line2D([], [], color=INK2, lw=1.2),
                Patch(facecolor=STRIP, edgecolor="none"),
                Line2D([], [], color=INK2, lw=0.6, ls=(0, (2, 2))),
                Line2D([], [], marker="v", ms=4.5, color=INK, mec="white", mew=0.4, ls="none"),
                Line2D([], [], marker="v", ms=4.5, mfc="white", mec=INK, mew=0.7, ls="none")]
    labels += ["pale: every frame (10 ps)", "strong: mean over 1 ns (production)",
               "grey strip: equilibration after the\nrelease of the ligand (−1.3 to 0 ns,\nstretched; every frame)",
               "ends of stages 4 and 5 (%d and %d ps\nafter the release; backbone held in\nstages 4 and 5, free in stage 6)" % tuple(ends),
               "conversion by the plan's rule: first\nframe of the first visit of at least\n1 ns in the opposite state",
               "the same frame in a run whose verdict\nis MIXED: not a conversion"]
    leg = lax.legend(handles, labels, loc="center left", bbox_to_anchor=(0.0, 0.5), fontsize=6.0, frameon=False,
                     handlelength=1.9, labelspacing=0.5, borderaxespad=0.0)
    for txt in leg.get_texts():
        txt.set_color(INK)
    fig.savefig(out_base + ".png", dpi=300)
    fig.savefig(out_base + ".pdf", dpi=300)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", default="/root/cqc-sfa-md/runs")
    ap.add_argument("--out", default=os.path.join(REPO_MD, "results"))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "conversion_times.json")) as fh:
        ct = json.load(fh)
    conv = {r["label"]: r.get("conversion") for r in ct["runs"]}
    labels = ["bonded_2e/" + rid for row in NINE for rid in row] + [lab for lab, _ in SENS]
    sources, runs = [], {}
    for lab in labels:
        runs[lab] = load_run(os.path.join(args.runs, lab), sources)
        assert len(runs[lab]["equil"]["t"]) == 130 and abs(runs[lab]["equil"]["t"][0] + 1.29) < 1e-9, lab
    make_figure("orientation", runs, conv, os.path.join(args.out, "fig_orientation"))
    make_figure("rmsd", runs, conv, os.path.join(args.out, "fig_ligand_rmsd"))
    rec = {"what": "Sources of md/results/fig_orientation.(png|pdf) and md/results/fig_ligand_rmsd.(png|pdf)",
           "script": {"name": "md/production/make_result_figures.py", "sha256": sha256_file(os.path.abspath(__file__))},
           "drawn_from": "analysis/ (first program); analysis_recheck/ (second program) read and compared",
           "tolerance_between_programs_A": TOL,
           "largest_difference_between_programs_in_the_columns_drawn_A":
               max(runs[lab][part]["max_abs_diff_between_programs"] for lab in labels for part in ("equil", "prod")),
           "running_mean_frames": SMOOTH, "midpoints_A": [M1, M2], "rule1_line_A": RULE1,
           "y_axis_A": {"orientation": [2.0, 12.0], "ligand_rmsd": [0.0, 4.5]},
           "values_outside_the_y_axis": {
               lab: {"d_below_2_or_above_12": int(sum(int(np.sum((runs[lab][p][k] < 2.0) | (runs[lab][p][k] > 12.0)))
                                                      for p in ("equil", "prod") for k in ("d1", "d2"))),
                     "rmsd_above_4.5": int(sum(int(np.sum(runs[lab][p]["rmsd"] > 4.5)) for p in ("equil", "prod")))}
               for lab in labels},
           "conversion_markers": {lab: conversion_marker(conv.get(lab)) for lab in labels},
           "conversion_marker_drawn_open_because_the_verdict_is_not_CONVERTED":
               [lab for lab in labels if conversion_marker(conv.get(lab)) and runs[lab]["verdict"] != "CONVERTED"],
           "stage_ends_drawn_ns": {lab: runs[lab]["stage_ends_ns"] for lab in labels},
           "stage_ends_ps_after_the_release": sorted({round((t + 1.3) * 1000.0) for lab in labels for t in runs[lab]["stage_ends_ns"]}),
           "colours": {"d1 and RMSD": BLUE, "d2": ORANGE},
           "files": sources,
           "figures": {f: sha256_file(os.path.join(args.out, f))
                       for f in ("fig_orientation.png", "fig_ligand_rmsd.png")}}
    with open(os.path.join(args.out, "fig_sources.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
        fh.write("\n")
    print("written: fig_orientation.png/.pdf, fig_ligand_rmsd.png/.pdf, fig_sources.json in %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
