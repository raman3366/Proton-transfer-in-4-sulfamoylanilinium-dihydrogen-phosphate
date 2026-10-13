#!/usr/bin/env python3
"""Compact time series of the eleven long runs for git and for the supplementary archive (rule 6 of md/PLAN.md:
"All nine runs enter the supplementary material as time series, whichever way they fall").

For every run one gzip-compressed CSV in md/results/series/, one line per solute frame (10 ps): the 130 frames
of the equilibration after the release of the ligand (stages 4, 5 and 6) at negative time, then the 10,000
frames of the production.  Nothing is computed from coordinates and no number is formed here: every cell is
the text of a cell of the analysis output of the first program (md/analysis/analyse_md.py),

    <runs>/<set>/<run>/analysis/timeseries_equil_solute.csv     rows with stage 4, 5 or 6
    <runs>/<set>/<run>/analysis/timeseries_solute.csv           all rows

except `time_ns`, which is written with two decimals (the files give it with six and with three), and `part`,
which is formed from the column `stage`.  Before a file is written the same columns of the second program
(<run>/analysis_recheck/) are read and must agree: numbers within 0.001 (the tolerance of section 14 of
md/analysis/SPEC.md), states, signs and flags exactly.  The script stops if they do not.

Written beside the series: README.md (every column and its source), SOURCES.json (every file read, with its
sha256, and what was compared) and SHA256SUMS.  The run folders are read only; nothing is written there.

    python3 md/production/make_series.py [--runs DIR] [--out DIR]
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MD = os.path.dirname(HERE)
RUNS = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + \
       [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
RELEASED_STAGES = (4, 5, 6)
TOL = 0.001                     # agreement of the two programs (SPEC.md, section 14)
N_EQUIL, N_PROD = 130, 10000    # frames expected: 1.3 ns and 100 ns at 10 ps
# (column of this file, column of the analysis output, kind)
COLUMNS = [("d1_A", "d1_N1_F131CZ_A", "number"),
           ("d2_A", "d2_N1_N67OD1_A", "number"),
           ("lig_rmsd_A", "lig_rmsd_A", "number"),
           ("state", "state", "text"),
           ("n2_volume", "n2_volume", "number"),
           ("n2_sign", "n2_sign", "text"),
           ("hb_T199", "hb_N2H_T199OG1", "text")]
HEADER = ["time_ns", "part"] + [c[0] for c in COLUMNS]
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ist_now():
    t = datetime.datetime.now(IST)
    return "%d %s %d, %d:%02d %s IST" % (t.day, t.strftime("%b"), t.year, (t.hour % 12) or 12, t.minute,
                                         "AM" if t.hour < 12 else "PM")


def read_rows(path, released_only):
    """[(time text, stage or None, [cells of COLUMNS])] of one analysis file."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        it = head.index("time_ns")
        ist_ = head.index("stage") if released_only else None
        idx = [head.index(src) for _, src, _ in COLUMNS]
        out = []
        for row in rd:
            if not row:
                continue
            stage = int(row[ist_]) if released_only else None
            if released_only and stage not in RELEASED_STAGES:
                continue
            out.append((row[it], stage, [row[i] for i in idx]))
    return out


def compare(a, b, what):
    """The two programs on the columns written: largest difference of a number; texts must be equal."""
    if len(a) != len(b):
        raise SystemExit("%s: the two programs wrote %d and %d rows" % (what, len(a), len(b)))
    worst = 0.0
    for (ta, sa, ca), (tb, sb, cb) in zip(a, b):
        if abs(float(ta) - float(tb)) > 1e-9 or sa != sb:
            raise SystemExit("%s: the two programs differ in time or stage" % what)
        for (name, _, kind), x, y in zip(COLUMNS, ca, cb):
            if kind == "number":
                worst = max(worst, abs(float(x) - float(y)))
            elif x != y:
                raise SystemExit("%s: the two programs differ in column %s (%s against %s)" % (what, name, x, y))
    if worst > TOL:
        raise SystemExit("%s: the two programs differ by %.4f in a number (tolerance %.3f)" % (what, worst, TOL))
    return worst


def one_run(runs_dir, model, run, out_dir):
    run_dir = os.path.join(runs_dir, model, run)
    rec = {"run": "%s/%s" % (model, run), "files_read": [], "largest_difference_between_the_two_programs": {}}
    parts = {}
    for part, fname, rel in (("equilibration", "timeseries_equil_solute.csv", True),
                             ("production", "timeseries_solute.csv", False)):
        data = {}
        for sub in ("analysis", "analysis_recheck"):
            p = os.path.join(run_dir, sub, fname)
            data[sub] = read_rows(p, rel)
            rec["files_read"].append({"file": p, "sha256": sha256_file(p), "rows_used": len(data[sub]),
                                      "use": "written" if sub == "analysis" else "compared"})
        rec["largest_difference_between_the_two_programs"][part] = compare(
            data["analysis"], data["analysis_recheck"], "%s/%s %s" % (model, run, fname))
        parts[part] = data["analysis"]
    eq, pr = parts["equilibration"], parts["production"]
    if len(eq) != N_EQUIL or len(pr) != N_PROD:
        raise SystemExit("%s/%s: %d released equilibration frames and %d production frames, %d and %d expected"
                         % (model, run, len(eq), len(pr), N_EQUIL, N_PROD))
    lines = []
    for t, stage, cells in eq:
        lines.append(["%.2f" % float(t), "eq%d" % stage] + cells)
    for t, _, cells in pr:
        lines.append(["%.2f" % float(t), "prod"] + cells)
    # the time axis: one frame every 10 ps from -1.29 ns to 100.00 ns without a gap
    for k, ln in enumerate(lines):
        if abs(float(ln[0]) - (-1.29 + 0.01 * k)) > 1e-6:
            raise SystemExit("%s/%s: line %d has time %s, %.2f expected" % (model, run, k, ln[0], -1.29 + 0.01 * k))
    buf = io.StringIO()
    wr = csv.writer(buf, lineterminator="\n")
    wr.writerow(HEADER)
    wr.writerows(lines)
    raw = buf.getvalue().encode("ascii")
    name = "%s_%s.csv.gz" % (model, run)
    out = os.path.join(out_dir, name)
    with open(out, "wb") as fh:                                   # no time stamp and no name in the gzip header:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=9) as gz:   # same bytes each time
            gz.write(raw)
    rec.update(file=name, rows=len(lines), rows_equilibration=len(eq), rows_production=len(pr),
               bytes=os.path.getsize(out), bytes_uncompressed=len(raw), sha256=sha256_file(out),
               sha256_uncompressed=hashlib.sha256(raw).hexdigest(),
               stage_ends_ns={"eq%d" % s: max(float(t) for t, st, _ in eq if st == s) for s in RELEASED_STAGES})
    return rec


README = """# Time series of the eleven runs of 100 ns (rule 6 of `md/PLAN.md`)

Written by `md/production/make_series.py` from the output of the first analysis program
(`md/analysis/analyse_md.py`; specification `md/analysis/SPEC.md`, version 1.3). One file per run, gzip-compressed
CSV, one line per solute frame (10 ps): {n_eq} frames of the equilibration after the release of the ligand
(stages 4, 5 and 6; 1.3 ns) followed by the {n_pr} frames of the production. {n_files} files, {n_rows} lines of
data each, {total_kb} kB together.

No number is formed in these files. Every cell is the text of a cell of `timeseries_equil_solute.csv` (rows of
stages 4 to 6) or `timeseries_solute.csv` in the folder `analysis/` of the run, except the first two columns as
described below. The same columns of the second program (`md/analysis/recheck_md.py`, folder
`analysis_recheck/`) were read and compared before each file was written: the largest difference of a number is
{worst:.4f} (tolerance 0.001), and every state, sign and flag is the same. `SOURCES.json` names every file
read with its sha256; `SHA256SUMS` holds the sums of the files of this folder.

## Columns

| Column | What it is | Source column |
|---|---|---|
| `time_ns` | Time in ns; 0 is the start of the production. The released part of the equilibration runs from −1.29 to 0.00, the production from 0.01 to 100.00. Written with two decimals | `time_ns` |
| `part` | `eq4`, `eq5`, `eq6`: equilibration stage 4 (protein backbone held with 250 kJ mol⁻¹ nm⁻²), stage 5 (50 kJ mol⁻¹ nm⁻²) and stage 6 (nothing held); `prod`: production. The ligand is released at the start of stage 4 | `stage` |
| `d1_A` | d1, distance from the aniline nitrogen N1 of the ligand to CZ of Phe131 (Å) | `d1_N1_F131CZ_A` |
| `d2_A` | d2, distance from N1 to OD1 of Asn67 (Å) | `d2_N1_N67OD1_A` |
| `lig_rmsd_A` | RMSD of the 11 heavy atoms of the ligand from the crystal pose of 6RL9 after the fit of the 258 Cα atoms; smallest of the four equivalent labellings (Å) | `lig_rmsd_A` |
| `state` | `C` crystal-like (d1 < 6.7 Å and d2 > 5.2 Å), `R` rotated-like (d1 > 6.7 Å and d2 < 5.2 Å), `I` intermediate (neither) | `state` |
| `n2_volume` | Signed volume of the unit vectors from the sulfonamide nitrogen N2 to Zn, H and S; positive in the crystal pose | `n2_volume` |
| `n2_sign` | Sign of `n2_volume`: 1, −1 or 0 | `n2_sign` |
| `hb_T199` | 1 if the hydrogen bond N2–H···OG1(Thr199) is present (N2···OG1 at most 3.5 Å and the angle at H at least 135°), else 0 | `hb_N2H_T199OG1` |

The verdict of a run reads the last 20 ns (the last 2,000 lines): CONVERTED or STAYED if at least 90 % of these
frames are in the state opposite to the start or in the state of the start, MIXED otherwise. The starts X and N
are crystal-like, the start R is rotated-like. The definitions are those of `md/analysis/SPEC.md` (sections 4,
9A and 9B) and of the plan (rules 1 to 3).

## Files

| File | Run | Lines of data (equilibration + production) | Size (bytes) | End of stage 4; 5 (ns) |
|---|---|---|---|---|
{rows}

`bonded_2s_R1` and `bonded_2e_n2fit_R1` are the two sensitivity runs of start R (the other zinc set; the
candidate terms at N2, which gate G-R2 did not accept). They are in no count of rule 3.

## To read a file

```python
import csv, gzip
with gzip.open('md/results/series/bonded_2e_R1.csv.gz', 'rt') as fh:
    rows = list(csv.DictReader(fh))
```
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", default="/root/cqc-sfa-md/runs")
    ap.add_argument("--out", default=os.path.join(REPO_MD, "results", "series"))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    recs = [one_run(args.runs, m, r, args.out) for m, r in RUNS]
    worst = max(v for rec in recs for v in rec["largest_difference_between_the_two_programs"].values())
    rows = "\n".join("| `%s` | %s | %d (%d + %d) | %d | %.2f; %.2f |"
                     % (rec["file"], rec["run"], rec["rows"], rec["rows_equilibration"], rec["rows_production"],
                        rec["bytes"], rec["stage_ends_ns"]["eq4"], rec["stage_ends_ns"]["eq5"]) for rec in recs)
    with open(os.path.join(args.out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(README.format(n_eq=N_EQUIL, n_pr="{:,}".format(N_PROD), n_files=len(recs),
                               n_rows="{:,}".format(N_EQUIL + N_PROD),
                               total_kb="{:,.0f}".format(sum(r["bytes"] for r in recs) / 1000.0),
                               worst=worst, rows=rows))
    with open(os.path.join(args.out, "SOURCES.json"), "w", encoding="utf-8") as fh:
        json.dump({"what": "Sources of the compact time series of md/results/series/ (rule 6 of md/PLAN.md)",
                   "script": {"name": "md/production/make_series.py", "sha256": sha256_file(os.path.abspath(__file__))},
                   "written": ist_now(), "runs_dir": args.runs,
                   "columns": [{"column": c, "source_column": s} for c, s, _ in COLUMNS],
                   "tolerance_between_the_two_programs": TOL,
                   "largest_difference_between_the_two_programs": worst,
                   "runs": recs}, fh, indent=1)
        fh.write("\n")
    with open(os.path.join(args.out, "SHA256SUMS"), "w", encoding="ascii") as fh:
        for rec in recs:
            fh.write("%s  %s\n" % (rec["sha256"], rec["file"]))
        for extra in ("README.md", "SOURCES.json"):
            fh.write("%s  %s\n" % (sha256_file(os.path.join(args.out, extra)), extra))
    print("written: %d series in %s (%d bytes); largest difference between the two programs %.4f"
          % (len(recs), args.out, sum(r["bytes"] for r in recs), worst))
    return 0


if __name__ == "__main__":
    sys.exit(main())
