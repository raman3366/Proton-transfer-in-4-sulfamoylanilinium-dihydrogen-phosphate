#!/usr/bin/env python3
"""The time of a conversion, counted from the release of the ligand.

Rule: md/PLAN.md, "Changes to this plan", entry of 7 Oct 2026, 12:10 AM IST, point 1 (it completes the entry
of 6 Oct 2026, 11:40 AM IST, point 2). In the words of the plan:

  The solute frames of equilibration stages 4, 5 and 6 (extension steps of stage 3 come before the release and
  do not count) and the solute frames of the production are joined in order of time. On the joined series the
  assigned states and the visits are formed by the rules of section 4 of the specification, unchanged (the
  state before the first frame is the state of the start; a long visit has at least 100 frames, 1 ns). The
  time of the conversion is the time of the first frame of the first long visit in the opposite state, counted
  from the start of stage 4. If that frame lies in the equilibration the run is reported as "converted during
  the equilibration, t ps after release"; if it lies in the production both times are given. A run without
  such a visit has no time. The number is taken from the `state` of each frame as written by both programs
  (timeseries_equil_solute.csv, timeseries_solute.csv of analysis/ and of analysis_recheck/), and the two must
  give the same frame.

Nothing is computed here from coordinates: only the column `state` (C, R or I) of the four files is read, with
`frame`, `stage`, `extension` and `time_ns`. The verdict of a run and rule 3 are untouched (they read the last
20 ns); the time describes one run and is not a rate.

Use:
    python md/production/conversion_time.py                       all runs below /root/cqc-sfa-md/runs
    python md/production/conversion_time.py --runs DIR --out FILE
The test with made-up series of known answer is md/production/test_conversion_time.py.
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys

DWELL_FRAMES = 100          # specification, section 4: a long visit has at least D = 100 frames (1 ns)
RELEASE_STAGE = 4           # the ligand is released at the start of equilibration stage 4
RELEASED_STAGES = (4, 5, 6)
OPPOSITE = {"C": "R", "R": "C"}
PROGRAMS = (("analysis", "first program (analyse_md.py)"),
            ("analysis_recheck", "second program (recheck_md.py)"))
EQUIL_FILE = "timeseries_equil_solute.csv"
PROD_FILE = "timeseries_solute.csv"


class Stop(Exception):
    pass


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_columns(path, names):
    """Rows of a CSV file as dicts of the named columns (text)."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        try:
            idx = [header.index(n) for n in names]
        except ValueError as err:
            raise Stop("%s: column missing (%s)" % (path, err))
        return [{n: row[i] for n, i in zip(names, idx)} for row in rd if row]


def joined_series(folder):
    """The joined series of one program's output folder.

    Returns a list of frames in order of time, each a dict with `part` ("equilibration" or "production"),
    `frame` (the frame number of its own file), `stage` and `extension` (equilibration only), `time_ns` (the
    time written in the file: zero or negative in the equilibration, positive in the production) and `state`.
    Of the equilibration only stages 4, 5 and 6 with extension 0 are taken."""
    out = []
    eq_path = os.path.join(folder, EQUIL_FILE)
    n_equil_all = 0
    if os.path.isfile(eq_path):
        rows = read_columns(eq_path, ["frame", "stage", "extension", "time_ns", "state"])
        n_equil_all = len(rows)
        for r in rows:
            stage, ext = int(r["stage"]), int(r["extension"])
            if stage in RELEASED_STAGES and ext == 0:
                out.append({"part": "equilibration", "frame": int(r["frame"]), "stage": stage, "extension": ext,
                            "time_ns": float(r["time_ns"]), "state": r["state"]})
    n_eq = len(out)
    rows = read_columns(os.path.join(folder, PROD_FILE), ["frame", "time_ns", "state"])
    for r in rows:
        out.append({"part": "production", "frame": int(r["frame"]), "stage": None, "extension": None,
                    "time_ns": float(r["time_ns"]), "state": r["state"]})
    # the series must be in order of time, without a doubled time
    for a, b in zip(out, out[1:]):
        if not b["time_ns"] > a["time_ns"]:
            raise Stop("%s: the joined series is not in order of time at %s frame %d (%.6f ns after %.6f ns)"
                       % (folder, b["part"], b["frame"], b["time_ns"], a["time_ns"]))
    for f in out:
        if f["state"] not in ("C", "R", "I"):
            raise Stop("%s: state %r in %s frame %d" % (folder, f["state"], f["part"], f["frame"]))
    return out, n_eq, n_equil_all


def assigned_states(states, start_state):
    """Section 4, point 1: the state of the frame if it is C or R; otherwise the assigned state of the frame
    before; before the first frame the state of the start."""
    out = []
    last = start_state
    for s in states:
        if s in ("C", "R"):
            last = s
        out.append(last)
    return out


def visits(assigned):
    """Section 4, point 2: maximal stretches of frames with the same assigned state, as (state, first index,
    number of frames)."""
    out = []
    i = 0
    n = len(assigned)
    while i < n:
        j = i
        while j + 1 < n and assigned[j + 1] == assigned[i]:
            j += 1
        out.append((assigned[i], i, j - i + 1))
        i = j + 1
    return out


def first_long_visit_opposite(states, start_state, dwell=DWELL_FRAMES):
    """Index (in the series given) of the first frame of the first long visit in the opposite state, or None.
    Also returns the list of all visits and the number of transitions between long visits (section 4, point 4)."""
    if start_state not in OPPOSITE:
        raise Stop("state of the start %r" % (start_state,))
    opp = OPPOSITE[start_state]
    vs = visits(assigned_states(states, start_state))
    first = None
    seq = [start_state]
    for st, i0, n in vs:
        if n >= dwell:
            seq.append(st)
            if first is None and st == opp:
                first = i0
    n_trans = sum(1 for a, b in zip(seq, seq[1:]) if a != b)
    return first, vs, n_trans


def analyse_folder(folder, start_state, dwell=DWELL_FRAMES, solute_interval_ps=10.0):
    """The conversion of one program's output folder (analysis/ or analysis_recheck/ of a run)."""
    series, n_eq, n_equil_all = joined_series(folder)
    states = [f["state"] for f in series]
    idx, vs, n_trans = first_long_visit_opposite(states, start_state, dwell)
    rec = {"folder": folder,
           "files": {EQUIL_FILE: (sha256_file(os.path.join(folder, EQUIL_FILE))
                                  if os.path.isfile(os.path.join(folder, EQUIL_FILE)) else None),
                     PROD_FILE: sha256_file(os.path.join(folder, PROD_FILE))},
           "n_equilibration_frames_in_file": n_equil_all,
           "n_released_equilibration_frames": n_eq,
           "n_production_frames": len(series) - n_eq,
           "n_joined_frames": len(series)}
    # the release: one frame interval before the first frame of stage 4
    if n_eq:
        f0 = series[0]
        if f0["stage"] != RELEASE_STAGE:
            raise Stop("%s: the first released frame is of stage %s, not %d" % (folder, f0["stage"], RELEASE_STAGE))
        t_release_ns = f0["time_ns"] - solute_interval_ps / 1000.0
    else:
        t_release_ns = None
    rec["release_time_ns_on_production_clock"] = None if t_release_ns is None else round(t_release_ns, 6)
    rec["released_part_of_equilibration_ps"] = None if t_release_ns is None else round(-t_release_ns * 1000.0, 3)
    # the time axis of the joined series must be (k + 1) intervals after the release
    if n_eq:
        worst = max(abs((f["time_ns"] - t_release_ns) * 1000.0 - (k + 1) * solute_interval_ps)
                    for k, f in enumerate(series))
        rec["time_axis_max_deviation_ps"] = round(worst, 6)
        if worst > 1e-3:
            raise Stop("%s: the joined time axis departs from (k + 1) x %.0f ps after the release by %.4f ps"
                       % (folder, solute_interval_ps, worst))
    long_vs = [{"state": st, "first_joined_index": i0, "n_frames": n,
                "part_of_first_frame": series[i0]["part"],
                "time_from_release_ps": (round((series[i0]["time_ns"] - t_release_ns) * 1000.0, 3)
                                         if t_release_ns is not None else None)}
               for st, i0, n in vs if n >= dwell]
    rec["n_visits"] = len(vs)
    rec["n_long_visits"] = len(long_vs)
    rec["long_visits"] = long_vs
    rec["n_transitions_joined"] = n_trans
    rec["state_first_released_frame"] = series[0]["state"] if series else None
    rec["counts_released_equilibration"] = {s: sum(1 for f in series[:n_eq] if f["state"] == s) for s in "CRI"}
    if idx is None:
        rec["conversion"] = None
    else:
        f = series[idx]
        conv = {"joined_index": idx, "part": f["part"], "in_equilibration": f["part"] == "equilibration",
                "frame_in_its_file": f["frame"], "stage": f["stage"], "extension": f["extension"],
                "time_ns_in_file": round(f["time_ns"], 6),
                "n_frames_of_the_visit": next(n for st, i0, n in vs if i0 == idx)}
        if t_release_ns is not None:
            conv["time_from_release_ps"] = round((f["time_ns"] - t_release_ns) * 1000.0, 3)
            conv["time_from_release_ns"] = round(f["time_ns"] - t_release_ns, 6)
        else:
            conv["time_from_release_ps"] = conv["time_from_release_ns"] = None
        conv["time_in_production_ns"] = round(f["time_ns"], 6) if f["part"] == "production" else None
        rec["conversion"] = conv
    return rec, states


def conversion_key(conv):
    return None if conv is None else (conv["part"], conv["frame_in_its_file"], conv["joined_index"])


def wording(conv):
    """The plan's wording for one run."""
    if conv is None:
        return "none (no visit of at least 1 ns in the opposite state after the release)"
    if conv["in_equilibration"]:
        return "converted during the equilibration, %.0f ps after release" % conv["time_from_release_ps"]
    if conv["time_from_release_ns"] is None:
        # a folder without timeseries_equil_solute.csv: the release is not on record (no run of the campaign)
        return ("%.2f ns in the production (no equilibration series in the folder: the time from the release "
                "is not known)" % conv["time_in_production_ns"])
    return "%.2f ns in the production, %.2f ns from the release" % (conv["time_in_production_ns"],
                                                                 conv["time_from_release_ns"])


def analyse_run(run_dir, label=None):
    """Both programs' series of one run folder."""
    rec = {"run_dir": run_dir, "label": label or os.path.basename(run_dir)}
    per = {}
    states = {}
    for sub, name in PROGRAMS:
        folder = os.path.join(run_dir, sub)
        with open(os.path.join(folder, "summary.json")) as fh:
            summ = json.load(fh)
        start_state = summ["orientation"]["start_state"]
        interval = float(summ["run"]["solute_interval_ps"])
        r, st = analyse_folder(folder, start_state, DWELL_FRAMES, interval)
        r["program"] = name
        r["start"] = summ["run"]["start"]
        r["model"] = summ["run"]["model"]
        r["role"] = summ["run"]["role"]
        r["state_of_start"] = start_state
        r["opposite_state"] = summ["orientation"]["opposite_state"]
        r["verdict_of_the_run_last_20_ns"] = summ["orientation"]["verdict"]
        r["first_passage_ns_production_only"] = summ["orientation"]["first_passage_ns"]
        r["first_passage_frame_production_only"] = summ["orientation"]["first_passage_frame"]
        r["equilibration_total_ps"] = (summ.get("equilibration") or {}).get("total_ps")
        r["equilibration_steps"] = [s.get("tag") for s in (summ.get("equilibration") or {}).get("stages", [])]
        per[sub] = r
        states[sub] = st
    a, b = per["analysis"], per["analysis_recheck"]
    rec["start"], rec["model"], rec["role"] = a["start"], a["model"], a["role"]
    rec["state_of_start"], rec["opposite_state"] = a["state_of_start"], a["opposite_state"]
    rec["programs"] = per
    same_len = len(states["analysis"]) == len(states["analysis_recheck"])
    rec["n_joined_frames_equal"] = same_len
    rec["n_frames_with_different_state"] = (sum(1 for x, y in zip(states["analysis"], states["analysis_recheck"])
                                                if x != y) if same_len else None)
    rec["same_frame"] = conversion_key(a["conversion"]) == conversion_key(b["conversion"])
    rec["same_start_state"] = a["state_of_start"] == b["state_of_start"]
    ok = rec["same_frame"] and rec["same_start_state"]
    rec["conversion"] = a["conversion"] if ok else None
    rec["reported_as"] = wording(a["conversion"]) if ok else "NOT REPORTED: the two programs give different frames"
    # the production-only first passage of the programs, set beside (it must be the same frame when the
    # conversion lies in the production and the production begins in the state of the start)
    rec["first_passage_ns_production_only"] = a["first_passage_ns_production_only"]
    return rec


def find_runs(top):
    """Run folders <top>/<model>/<run> that hold the output of both programs. Folders directly below <top>
    (the pre-flight tests of the rented machine) are not runs of the campaign and are passed over."""
    out = []
    for model in sorted(os.listdir(top)):
        mdir = os.path.join(top, model)
        if not os.path.isdir(mdir) or model.startswith("analysis"):
            continue
        for run in sorted(os.listdir(mdir)):
            d = os.path.join(mdir, run)
            if all(os.path.isfile(os.path.join(d, sub, PROD_FILE)) and os.path.isfile(os.path.join(d, sub, "summary.json"))
                   for sub, _ in PROGRAMS):
                out.append((model + "/" + run, d))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", default="/root/cqc-sfa-md/runs")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                                                  "conversion_times.json"))
    args = ap.parse_args(argv)
    top = os.path.abspath(args.runs)
    runs = find_runs(top)
    if not runs:
        print("no run folder with the output of both programs below %s" % top, file=sys.stderr)
        return 2
    ist = datetime.timezone(datetime.timedelta(hours=5, minutes=30))
    out = {"what": "Time of a conversion from the release of the ligand, for every run, by the rule of md/PLAN.md "
                   "(entry of 7 Oct 2026, 12:10 AM IST, point 1). Read from the column `state` of "
                   "timeseries_equil_solute.csv (stages 4, 5, 6) and timeseries_solute.csv of both programs.",
           "rule": {"released_stages": list(RELEASED_STAGES), "dwell_frames": DWELL_FRAMES,
                    "assigned_state": "the state of the frame if C or R, else that of the frame before; before "
                                      "the first joined frame the state of the start (specification, section 4)",
                    "time": "first frame of the first long visit in the opposite state, counted from the start "
                            "of stage 4"},
           "script": {"name": "md/production/conversion_time.py", "sha256": sha256_file(os.path.abspath(__file__))},
           "runs_dir": top,
           "written": datetime.datetime.now(ist).strftime("%-d %b %Y, %-I:%M %p IST"),
           "runs": []}
    status = 0
    for label, d in runs:
        try:
            rec = analyse_run(d, label)
        except Stop as err:
            rec = {"run_dir": d, "label": label, "error": str(err)}
            status = 1
        out["runs"].append(rec)
        if "error" in rec:
            print("%-22s ERROR %s" % (label, rec["error"]))
        else:
            if not rec["same_frame"]:
                status = 1
            print("%-22s start %s (%s)  %-66s same frame in both programs: %s"
                  % (label, rec["start"], rec["state_of_start"], rec["reported_as"], rec["same_frame"]))
    out["all_runs_same_frame_in_both_programs"] = all(r.get("same_frame") for r in out["runs"])
    path = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print("written to %s" % path)
    return status


if __name__ == "__main__":
    sys.exit(main())
