#!/usr/bin/env python
"""pi_test_report.py: speed of the test runs, the processor time they really got, and what that means
for a validation run on the Raspberry Pi.

TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper (md/test/PI_TEST.md).

  python pi_test_report.py --runs /root/cqc-sfa-md/runs_test --out md/test/results/pi_test_speed.json

Reads the step records of every run below RUNS/restrained/ (rec_<step>.json: wall time, number of steps)
and RUNS/cpu_samples.csv (pi_test_cpu_sampler.sh: processor time of the run_md.py process and of the whole
machine, once a minute), and RUNS/bench_*/benchmark.json if pi_test_bench.sh has been run.

For every step it gives: steps per second, ns per day, the cores the process received (its processor time
divided by the wall time), the cores the rest of the machine used, the load average, and steps per second
of processor time. The last number is what carries over to another load on the same machine.

Estimates (labelled as such; nothing here was measured with 4 threads on an idle machine):
  the plan's equilibration (2 ns: 700 ps at 2 fs and 1,300 ps at 4 fs = 675,000 steps) and its 5 ns
  validation run (1,250,000 steps at 4 fs), for a given number of cores received by the process.

Units: s, steps, ns per day, cores (processor seconds per second of wall time).
"""
import argparse
import bisect
import calendar
import glob
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_test_common as C  # noqa: E402


def epoch(utc):
    return calendar.timegm(time.strptime(utc, "%Y-%m-%dT%H:%M:%SZ"))


def load_samples(path):
    by_pid, machine = {}, []
    if not os.path.isfile(path):
        return by_pid, machine
    with open(path) as f:
        f.readline()
        for ln in f:
            p = ln.strip().split(",")
            if len(p) < 10:
                continue
            t, run, pid = float(p[0]), p[1], int(p[2])
            if run == "machine":
                machine.append((t, float(p[4]), float(p[5]), float(p[6])))
            else:
                by_pid.setdefault(pid, []).append((t, float(p[3])))
    return by_pid, machine


def interp(samples, t, col=1):
    """Linear interpolation of a cumulative counter; None outside the sampled time."""
    ts = [s[0] for s in samples]
    if not ts or t < ts[0] - 90 or t > ts[-1] + 90:
        return None
    k = bisect.bisect_left(ts, t)
    if k == 0:
        return samples[0][col]
    if k >= len(ts):
        return samples[-1][col]
    a, b = samples[k - 1], samples[k]
    return a[col] + (b[col] - a[col]) * (t - a[0]) / (b[0] - a[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default="/root/cqc-sfa-md/runs_test")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    by_pid, machine = load_samples(os.path.join(args.runs, "cpu_samples.csv"))
    out = {"note": C.NOTE, "written": C.now_ist(), "runs_dir": os.path.abspath(args.runs),
           "units": "wall_s and cpu_s in seconds; cores = processor seconds per second of wall time; ns_per_day",
           "software": C.software(), "script_sha256": C.sha256_file(os.path.abspath(__file__)), "runs": {}}
    pooled = {"2.0": [], "4.0": []}
    for run_dir in sorted(glob.glob(os.path.join(args.runs, "restrained", "*"))):
        if not os.path.isfile(os.path.join(run_dir, "manifest.json")):
            continue
        man = C.read_json(os.path.join(run_dir, "manifest.json"))
        plan = C.expected_steps(man)
        inv = [(epoch(i["started"]["utc"]), i.get("pid")) for i in man["invocations"]]
        steps = []
        for st in plan:
            rp = os.path.join(run_dir, f"rec_{st['tag']}.json")
            if not os.path.isfile(rp):
                continue
            rec = C.read_json(rp)
            t0, t1 = epoch(rec["started"]["utc"]), epoch(rec["ended"]["utc"])
            pid = [p for t, p in inv if t <= t0 + 1][-1]
            row = {"step": st["tag"], "started_ist": rec["started"]["ist"], "ended_ist": rec["ended"]["ist"],
                   "wall_s_record": rec["wall_s"], "wall_s_start_to_end_of_step": t1 - t0, "pid": pid,
                   "threads": rec["platform"]["properties_set"].get("Threads"), "platform": rec["platform"]["name"]}
            samples = by_pid.get(pid, [])
            c0, c1 = interp(samples, t0), interp(samples, t1)
            if t0 - inv[[p for _t, p in inv].index(pid)][0] < 30 and c0 is None:
                c0 = 0.0                                    # the step began with the process
            if c0 is not None and c1 is not None and t1 > t0:
                row["process_cpu_s"] = round(c1 - c0, 1)
                row["cores_received"] = round((c1 - c0) / (t1 - t0), 3)
            b0, b1 = interp(machine, t0, 1), interp(machine, t1, 1)
            i0, i1 = interp(machine, t0, 2), interp(machine, t1, 2)
            if None not in (b0, b1, i0, i1) and t1 > t0:
                row["machine_busy_cores"] = round((b1 - b0) / (t1 - t0), 3)
                row["machine_idle_cores"] = round((i1 - i0) / (t1 - t0), 3)
                if "cores_received" in row:
                    row["cores_used_by_everything_else"] = round(row["machine_busy_cores"] - row["cores_received"], 3)
                lo = [m[3] for m in machine if t0 <= m[0] <= t1]
                row["load_average_1min_mean"] = round(float(np.mean(lo)), 2) if lo else None
            if st["kind"] == "minimisation":
                n_it = rec["minimiser_calls"]
                row.update({"kind": "minimisation", "lbfgs_iterations": n_it,
                            "lbfgs_iterations_per_pass": rec["lbfgs_iterations_per_pass"],
                            "s_per_iteration": round(rec["wall_s"] / max(n_it, 1), 4)})
                if "process_cpu_s" in row:
                    row["cpu_s_per_iteration"] = round(row["process_cpu_s"] / max(n_it, 1), 4)
            else:
                n = st["n_steps"]
                row.update({"kind": st["kind"], "timestep_fs": st["timestep_fs"], "n_steps": n,
                            "steps_per_s": round(n / rec["wall_s"], 4), "ns_per_day": rec["ns_per_day"],
                            "ensemble": rec["settings"]["ensemble"], "n_restrained_atoms": rec["settings"]["n_restrained_atoms"]})
                if "process_cpu_s" in row:
                    row["steps_per_cpu_s"] = round(n / row["process_cpu_s"], 4)
                    pooled[str(st["timestep_fs"])].append((n, rec["wall_s"], row["process_cpu_s"], st["kind"], man["run_id"]))
            steps.append(row)
        # a step that was cut off: its state log gives the steps done and, through the cumulative speed column
        # (two decimals), the wall time of every row
        cut = next((st for st in plan if not os.path.isfile(os.path.join(run_dir, f"rec_{st['tag']}.json"))), None)
        cut_row = None
        lp = os.path.join(run_dir, f"log_{cut['tag']}.csv") if cut else None
        if cut and cut["kind"] != "minimisation" and os.path.isfile(lp) and os.path.getsize(lp) > 0:
            header, log = C.read_log_csv(lp)
            t_start = None
            with open(os.path.join(run_dir, "run_md.log")) as f:
                for ln in f:
                    if f"] {cut['tag']} (" in ln and " steps of " in ln:
                        stamp = ln[1:ln.index("]")].replace(" IST", "")
                        t_start = calendar.timegm(time.strptime(stamp, "%d %b %Y, %I:%M:%S %p")) - 19800
            if len(log) and t_start is not None:
                dt_fs, every = cut["timestep_fs"], cut["every_steps"]["state_log_interval_ps"]
                ns = (log[:, header.index("time_ps")] - cut["time_start_ps"]) * 1e-3
                elapsed = ns / log[:, header.index("speed_ns_per_day")] * 86400.0
                pid = [p for t, p in inv if t <= t_start + 1][-1]
                samples = by_pid.get(pid, [])
                rows_, prev_e = [], 0.0
                for k in range(len(log)):
                    a, b = t_start + prev_e, t_start + float(elapsed[k])
                    r_ = {"to_time_ps": float(log[k, header.index("time_ps")]), "wall_s": round(float(elapsed[k]) - prev_e, 0),
                          "steps_per_s": round(every / (float(elapsed[k]) - prev_e), 3),
                          "ns_per_day": round(every * dt_fs * 1e-6 / (float(elapsed[k]) - prev_e) * 86400.0, 3)}
                    c0, c1 = interp(samples, a), interp(samples, b)
                    if c0 is not None and c1 is not None:
                        r_["cores_received"] = round((c1 - c0) / (b - a), 3)
                        b0, b1 = interp(machine, a, 1), interp(machine, b, 1)
                        if b0 is not None and b1 is not None:
                            r_["cores_used_by_everything_else"] = round((b1 - b0) / (b - a) - r_["cores_received"], 3)
                        lo = [m[3] for m in machine if a <= m[0] <= b]
                        r_["load_average_1min_mean"] = round(float(np.mean(lo)), 2) if lo else None
                    rows_.append(r_)
                    prev_e = float(elapsed[k])
                n_done = every * len(log)
                c0, c1 = interp(samples, t_start), interp(samples, t_start + prev_e)
                cut_row = {"step": cut["tag"], "kind": cut["kind"] + " (cut off; from its state log)", "timestep_fs": dt_fs,
                           "n_steps": n_done, "n_steps_planned": cut["n_steps"], "wall_s_to_the_last_row": round(prev_e, 0),
                           "steps_per_s": round(n_done / prev_e, 4), "ns_per_day": round(n_done * dt_fs * 1e-6 / prev_e * 86400.0, 3),
                           "threads": man["invocations"][-1]["platform"]["properties_set"].get("Threads"),
                           "wall_time_source": "cumulative speed column of the log (two decimals): each wall time is good to about 1 %",
                           "rows": rows_}
                if c0 is not None and c1 is not None:
                    cut_row["process_cpu_s"] = round(c1 - c0, 1)
                    cut_row["cores_received"] = round((c1 - c0) / prev_e, 3)
                    cut_row["steps_per_cpu_s"] = round(n_done / (c1 - c0), 4)
                    pooled[str(dt_fs)].append((n_done, prev_e, c1 - c0, cut["kind"], man["run_id"]))
        prod = [s for s in steps if s["kind"] == "production"]
        summary = {}
        if prod:
            n = sum(s["n_steps"] for s in prod)
            w = sum(s["wall_s_record"] for s in prod)
            summary["production"] = {"segments": len(prod), "n_steps": n, "wall_s": round(w, 1), "steps_per_s": round(n / w, 4),
                                     "ns_per_day": round(n * 4e-6 / w * 86400.0, 3),
                                     "cores_received": (round(sum(s["process_cpu_s"] for s in prod) / w, 3)
                                                        if all("process_cpu_s" in s for s in prod) else None)}
        out["runs"][man["run_id"]] = {"run_dir": run_dir, "n_atoms": man["system"]["n_atoms"], "steps": steps,
                                      "cut_off_step": cut_row, "summary": summary,
                                      "invocations": [{"pid": i.get("pid"), "started_ist": i["started"]["ist"],
                                                       "exit_code": i.get("exit_code"), "outcome": i.get("outcome")}
                                                      for i in man["invocations"]]}

    # steps per second of processor time, pooled over the MD steps of all runs
    rate = {}
    for dt, rows in pooled.items():
        if rows:
            n, w, c = (sum(r[k] for r in rows) for k in (0, 1, 2))
            rate[dt] = {"n_steps": n, "wall_s": round(w, 1), "process_cpu_s": round(c, 1), "steps_per_s": round(n / w, 4),
                        "steps_per_cpu_s": round(n / c, 4), "cores_received": round(c / w, 3)}
    out["pooled_over_md_steps_by_timestep_fs"] = rate

    # the barostat in a run with a state-log row at every barostat attempt (the dense-output test of md/TOPOLOGY.md,
    # start N, read by pi_test_look.py): how often and by how much the volume changed, per step of the chain
    here = os.path.dirname(os.path.abspath(__file__))
    dense = os.path.join(here, "results", "N_min300_of_TOPOLOGY_log.csv")
    if os.path.isfile(dense):
        import csv
        with open(dense) as f:
            rows = list(csv.DictReader(ln for ln in f if not ln.startswith("#")))
        table = {}
        for ph in dict.fromkeys(r["phase"] for r in rows):
            v = np.array([float(r["volume_nm3"]) for r in rows if r["phase"] == ph])
            t = np.array([float(r["time_ps"]) for r in rows if r["phase"] == ph])
            if len(v) < 2:
                continue
            dv = np.diff(v) / v[:-1]
            ch = np.abs(dv) > 1e-9
            table[ph] = {"log_interval_ps": round(float(t[1] - t[0]), 3), "intervals": int(len(dv)),
                         "intervals_with_a_volume_change": int(ch.sum()),
                         "relative_volume_changes_percent": [round(float(x) * 100.0, 3) for x in dv[ch]]}
        out["barostat_volume_changes_in_the_dense_output_test_of_TOPOLOGY"] = {
            "file": dense, "sha256": C.sha256_file(dense),
            "note": "start N, restrained model, state log every 0.1 ps: one barostat attempt per row in the 4 fs steps "
                    "(eq04 on), two in the 2 fs steps; not a run of this test", "by_step": table}
    bench = []
    for bp in sorted(glob.glob(os.path.join(args.runs, "bench_*", "benchmark.json"))):
        for b in C.read_json(bp):
            entry = {"folder": os.path.basename(os.path.dirname(bp)), "threads": b["platform"]["properties_set"].get("Threads"),
                     "ns_per_day": b["ns_per_day"], "steps": b["steps"], "wall_s": b["wall_s"],
                     "steps_per_s": round(b["steps"] / b["wall_s"], 4), "when_ist": b["when"]["ist"]}
            side = os.path.join(os.path.dirname(bp), "bench_cpu.json")
            if os.path.isfile(side):
                entry.update(C.read_json(side))
            bench.append(entry)
    out["benchmarks"] = bench

    # what the plan's validation run would cost on this machine
    all_rows = [r for rows in pooled.values() for r in rows]
    if all_rows:
        n, c = sum(r[0] for r in all_rows), sum(r[2] for r in all_rows)
        k = n / c                                            # steps per second of processor time, 2 threads
        eq_steps = 0
        proto = None
        for run in out["runs"].values():
            proto = C.read_json(os.path.join(run["run_dir"], "manifest.json"))["protocol"]
            break
        for st in proto["equilibration"]["stages"]:
            eq_steps += int(round(st["length_ps"] * 1000.0 / st["timestep_fs"]))
        val_steps = int(round(5.0 * 1e6 / proto["production"]["timestep_fs"]))
        mins = [s for run in out["runs"].values() for s in run["steps"] if s["kind"] == "minimisation" and "process_cpu_s" in s]
        min_cpu = float(np.mean([s["process_cpu_s"] for s in mins])) if mins else None
        est = {"basis": "steps per second of processor time measured here with 2 threads, pooled over all MD steps of the "
                        "test runs; the time for another load is steps / (that rate x cores received). ESTIMATE, not measured.",
               "steps_per_cpu_s_measured_2_threads": round(k, 4),
               "equilibration_steps_of_the_plan": eq_steps, "validation_steps_5_ns": val_steps,
               "minimisation_cpu_s_measured_mean": round(min_cpu, 1) if min_cpu else None, "scenarios": []}
        for label, cores, eff in (
                ("4 threads, nothing else running: 4 cores, threading as efficient as with 2 threads (upper bound)", 4.0, 1.00),
                ("4 threads, nothing else running: 4 cores, 15 % lost to threading (assumed)", 4.0, 0.85),
                ("2 threads beside a 4-thread quantum-chemical job, two cores for each job (the scheduler's even split between the two jobs)", 2.0, 1.00),
                ("2 threads beside a 4-thread quantum-chemical job, cores shared evenly per thread (2 of 6 threads on 4 cores)", 4.0 * 2 / 6, 1.00),
                ("2 threads with the cores this test really received", None, 1.00)):
            if cores is None:
                w = sum(r[1] for r in all_rows)
                cores = c / w
            sps = k * cores * eff
            est["scenarios"].append({
                "scenario": label, "cores_received": round(cores, 3), "threading_factor": eff,
                "steps_per_s": round(sps, 3), "ns_per_day_at_4_fs": round(sps * 4e-6 * 86400.0, 2),
                "equilibration_2_ns_hours": round(eq_steps / sps / 3600.0, 1),
                "validation_5_ns_hours": round(val_steps / sps / 3600.0, 1),
                "minimisation_hours": round(min_cpu / (cores * eff) / 3600.0, 2) if min_cpu else None,
                "all_three_days": round(((eq_steps + val_steps) / sps + (min_cpu / (cores * eff) if min_cpu else 0.0)) / 86400.0, 2)})
        out["estimate_validation_run_on_the_pi"] = est
    C.write_json(args.out, out)
    for rid, run in out["runs"].items():
        for s in run["steps"]:
            print(rid, s["step"], s.get("steps_per_s", s.get("s_per_iteration")), "ns/day", s.get("ns_per_day"),
                  "cores", s.get("cores_received"), "others", s.get("cores_used_by_everything_else"), "load", s.get("load_average_1min_mean"))
    print("written", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
