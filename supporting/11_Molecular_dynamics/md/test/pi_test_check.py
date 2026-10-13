#!/usr/bin/env python
"""pi_test_check.py: the files of a test run against what the run should hold.

TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper (md/test/PI_TEST.md).

  python pi_test_check.py snapshot RUN_DIR --killed-pid PID --signal SIGKILL --out FILE
      what an interrupted run left behind: every file with size and sha256, and for each step without a
      record the frames and log rows that were on disk.
  python pi_test_check.py check RUN_DIR [--snapshot FILE] [--incomplete] --out FILE
      a finished run: frame counts, log rows, time axis, the chain of saved states, checksums, the last
      frame against the saved state, the two kinds of trajectory against each other, repeated frames, and
      (with --snapshot) what the interruption did. With --incomplete the run may have been cut off: the
      complete steps are checked in the same way and the cut-off step is described.
  python pi_test_check.py recovery DIR --snapshot FILE --out FILE
      copies of a killed run folder (DIR/<label>/X1) after "run_md.py --stage min --resume": did the driver
      take the lock, mark the killed invocation, keep the complete steps and remove the cut-off files?

The expected numbers are computed here from the protocol recorded in the manifest (pi_test_common.py,
expected_steps); the DCD files are read with this folder's own reader and, as a second opinion, with
MDAnalysis. Units: Angstrom, ps.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pi_test_common as C  # noqa: E402


def step_of(name):
    """The step tag of a file name of a run folder (None for manifest, PDB files and the like)."""
    base = os.path.splitext(name)[0]
    if "_" not in base:
        return None
    tag = base.rsplit("_", 1)[1]
    return tag if (tag == "min" or tag[:2] == "eq" or tag[:3] == "seg") else None


def list_files(run_dir):
    out = {}
    for p in sorted(glob.glob(os.path.join(run_dir, "*"))):
        if os.path.isfile(p):
            out[os.path.basename(p)] = {"bytes": os.path.getsize(p), "sha256": C.sha256_file(p)}
    return out


def snapshot(args):
    run_dir = os.path.abspath(args.run_dir)
    files = list_files(run_dir)
    complete = sorted(os.path.basename(p)[4:-5] for p in glob.glob(os.path.join(run_dir, "rec_*.json")))
    open_steps = {}
    for name in files:
        tag = step_of(name)
        if tag is None or tag in complete:
            continue
        entry = open_steps.setdefault(tag, {"files": {}})
        info = {"bytes": files[name]["bytes"]}
        if name.endswith(".dcd"):
            info.update(C.read_dcd_header(os.path.join(run_dir, name)))
        if name.endswith(".csv"):
            with open(os.path.join(run_dir, name)) as f:
                lines = [ln for ln in f.read().split("\n")]
            info["state_log_rows"] = max(0, len([ln for ln in lines if ln.strip()]) - 1)
            info["last_line_complete"] = lines[-1] == ""
            info["last_row"] = [ln for ln in lines if ln.strip()][-1] if len(lines) > 1 else None
        entry["files"][name] = info
    man = C.read_json(os.path.join(run_dir, "manifest.json"))
    lock = os.path.join(run_dir, "run_md.lock")
    result = {
        "note": C.NOTE, "what": "state of the run folder right after run_md.py was killed, before --resume",
        "run_dir": run_dir, "taken": C.now_ist(), "killed_pid": args.killed_pid, "signal": args.signal,
        "load_average_1_5_15_min": open("/proc/loadavg").read().split()[:3],
        "steps_with_a_record": complete, "steps_without_a_record": open_steps,
        "lock_file_left": open(lock).read().strip() if os.path.isfile(lock) else None,
        "manifest_status": man.get("status"),
        "manifest_last_invocation": {k: man["invocations"][-1].get(k) for k in ("pid", "started", "ended", "exit_code", "outcome")},
        "files": files, "software": C.software(),
    }
    C.write_json(args.out, result)
    for tag, entry in open_steps.items():
        for name, info in entry["files"].items():
            extra = (f", {info.get('whole_frames_by_size')} whole frames (header says {info.get('header_n_frames')}), "
                     f"{info.get('bytes_beyond_whole_frames')} bytes beyond" if name.endswith(".dcd") else
                     f", {info.get('state_log_rows')} rows" if name.endswith(".csv") else "")
            print(f"left without a record: {name}, {info['bytes']} bytes{extra}")
    print(f"snapshot written to {args.out}")
    return 0


def check(args):
    run_dir = os.path.abspath(args.run_dir)
    man = C.read_json(os.path.join(run_dir, "manifest.json"))
    plan = C.expected_steps(man)
    sel = C.read_json(os.path.join(man["system_dir"], "selections.json"))
    solute = np.array(sel["solute"])
    tests, per_step = [], []

    def add(name, value, passed):
        tests.append({"test": name, "value": value, "passed": bool(passed)})

    inc = bool(getattr(args, "incomplete", False))
    add("manifest status" + (" (a run that was cut off: 'running' is what a killed process leaves)" if inc else ""),
        f"status {man.get('status')}, all_steps_complete {man.get('all_steps_complete')}, "
        f"test_only {man.get('test_only')} ({'; '.join(man.get('test_only_reason', []))}), TEST_ONLY file "
        f"{'present' if os.path.isfile(os.path.join(run_dir, 'TEST_ONLY')) else 'absent'}",
        (inc or (man.get("status") == "complete" and man.get("all_steps_complete"))) and man.get("test_only")
        and os.path.isfile(os.path.join(run_dir, "TEST_ONLY")))
    cut_off = None

    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.formats.libdcd import DCDFile
    except Exception:                                         # pragma: no cover
        mda = None

    prev = None
    bad_counts, bad_sha, bad_chain, bad_time, bad_mda = [], [], [], [], []
    last_frame_vs_state, full_vs_solute = {}, {}
    seg_data = {}
    for st in plan:
        tag = st["tag"]
        rp = os.path.join(run_dir, f"rec_{tag}.json")
        if not os.path.isfile(rp):
            if not inc:
                bad_counts.append(f"{tag}: no record")
                prev = st
                continue
            cut_off = {"step": tag, "expected": {k: st.get(k) for k in ("n_steps", "state_log_rows", "solute", "full")}, "files": {}}
            for name in sorted(os.listdir(run_dir)):
                if step_of(name) != tag:
                    continue
                info = {"bytes": os.path.getsize(os.path.join(run_dir, name))}
                if name.endswith(".dcd"):
                    info.update(C.read_dcd_header(os.path.join(run_dir, name)))
                if name.endswith(".csv"):
                    hd, lg = C.read_log_csv(os.path.join(run_dir, name))
                    info["state_log_rows"] = len(lg)
                    if len(lg):
                        info["last_time_ps"] = float(lg[-1, hd.index("time_ps")])
                cut_off["files"][name] = info
            break
        rec = C.read_json(rp)
        row = {"step": tag, "kind": st["kind"], "wall_s": rec.get("wall_s"), "ns_per_day": rec.get("ns_per_day"),
               "started_ist": rec["started"]["ist"], "ended_ist": rec["ended"]["ist"]}
        for name, info in rec["files"].items():
            p = os.path.join(run_dir, name)
            if not os.path.isfile(p) or os.path.getsize(p) != info["bytes"] or C.sha256_file(p) != info["sha256"]:
                bad_sha.append(name)
        if st["kind"] != "minimisation":
            # the chain of saved states
            prev_rec = C.read_json(os.path.join(run_dir, f"rec_{prev['tag']}.json"))
            want = prev_rec["files"][f"state_{prev['tag']}.xml"]["sha256"]
            if rec["state_in"]["sha256"] != want or rec["state_in"]["file"] != f"state_{prev['tag']}.xml":
                bad_chain.append(tag)
            # counts
            hs = C.read_dcd_header(os.path.join(run_dir, f"solute_{tag}.dcd"))
            hf = C.read_dcd_header(os.path.join(run_dir, f"full_{tag}.dcd"))
            header, log = C.read_log_csv(os.path.join(run_dir, f"log_{tag}.csv"))
            got = {"state_log_rows": len(log), "solute": hs["whole_frames_by_size"], "full": hf["whole_frames_by_size"]}
            row.update({"n_steps": st["n_steps"], "timestep_fs": st["timestep_fs"], "length_ps": st["length_ps"],
                        "time_start_ps": st["time_start_ps"],
                        "expected": {k: st[k] for k in ("state_log_rows", "solute", "full")}, "found": got,
                        "record_frames": rec["frames"],
                        "dcd_header_n_frames": {"solute": hs["header_n_frames"], "full": hf["header_n_frames"]},
                        "bytes_beyond_whole_frames": {"solute": hs["bytes_beyond_whole_frames"], "full": hf["bytes_beyond_whole_frames"]}})
            for k in ("state_log_rows", "solute", "full"):
                if got[k] != st[k] or rec["frames"][k] != st[k]:
                    bad_counts.append(f"{tag} {k}: expected {st[k]}, file {got[k]}, record {rec['frames'][k]}")
            for kind, h in (("solute", hs), ("full", hf)):
                if h["header_n_frames"] != h["whole_frames_by_size"] or h["bytes_beyond_whole_frames"] != 0:
                    bad_counts.append(f"{tag} {kind}: header {h['header_n_frames']} frames, size {h['whole_frames_by_size']} "
                                      f"and {h['bytes_beyond_whole_frames']} bytes beyond")
            if rec["settings"]["n_steps"] != st["n_steps"] or abs(rec["time_start_ps"] - st["time_start_ps"]) > 1e-6:
                bad_counts.append(f"{tag}: record has {rec['settings']['n_steps']} steps from {rec['time_start_ps']} ps, "
                                  f"expected {st['n_steps']} from {st['time_start_ps']}")
            # time axis
            every = st["every_steps"]
            dt_ps = st["timestep_fs"] / 1000.0
            want_t = st["time_start_ps"] + (np.arange(st["state_log_rows"]) + 1) * every["state_log_interval_ps"] * dt_ps
            if len(log) != len(want_t) or (len(log) and np.abs(log[:, header.index("time_ps")] - want_t).max() > 1e-6):
                bad_time.append(f"log_{tag}.csv")
            head0 = int(round(st["time_start_ps"] / dt_ps))
            for kind, h, key in (("solute", hs, "solute_frame_interval_ps"), ("full", hf, "full_frame_interval_ps")):
                if (h["header_first_step"] != head0 + every[key] or h["header_step_interval"] != every[key]
                        or abs(h["header_timestep_ps"] - dt_ps) > 1e-7):
                    bad_time.append(f"{kind}_{tag}.dcd header: first step {h['header_first_step']}, interval "
                                    f"{h['header_step_interval']}, step {h['header_timestep_ps']} ps")
            # second opinion on the frame counts
            if mda is not None:
                for kind, h in (("solute", hs), ("full", hf)):
                    if h["whole_frames_by_size"] == 0:
                        continue
                    with DCDFile(os.path.join(run_dir, f"{kind}_{tag}.dcd")) as d:
                        if d.n_frames != st[kind] or d.header["natoms"] != h["n_atoms"]:
                            bad_mda.append(f"{kind}_{tag}.dcd: MDAnalysis reads {d.n_frames} frames")
            # coordinates
            _h, xs, cs = C.read_dcd(os.path.join(run_dir, f"solute_{tag}.dcd"))
            _h, xf, cf = C.read_dcd(os.path.join(run_dir, f"full_{tag}.dcd"))
            if len(xs) and st["n_steps"] % every["solute_frame_interval_ps"] == 0:
                pos, _vel, _cell = C.read_state_xml(os.path.join(run_dir, f"state_{tag}.xml"))
                last_frame_vs_state[tag] = float(np.abs(xs[-1] - pos[solute]).max())
            if len(xf):
                ratio = every["full_frame_interval_ps"] // every["solute_frame_interval_ps"]
                d = [float(np.abs(xf[j][solute] - xs[ratio * (j + 1) - 1]).max()) for j in range(len(xf))]
                full_vs_solute[tag] = max(d)
            if st["kind"] == "production":
                seg_data[tag] = (hs, xs)
            row["mean_temperature_K"] = (rec.get("mean_over_state_log") or {}).get("temperature_K")
            row["mean_density_g_per_cm3"] = (rec.get("mean_over_state_log") or {}).get("density_g_per_cm3")
            row["final_temperature_K"] = rec["final"]["temperature_K"]
            row["final_density_g_per_cm3"] = rec["final"]["density_g_per_cm3"]
        else:
            row.update({"lbfgs_iterations_per_pass": rec.get("lbfgs_iterations_per_pass"),
                        "minimiser_calls": rec.get("minimiser_calls"),
                        "passes_stopped_by_the_iteration_limit": rec.get("passes_stopped_by_the_iteration_limit"),
                        "last_pass_stopped_by_the_iteration_limit": rec.get("last_pass_stopped_by_the_iteration_limit"),
                        "energy_before_kj_per_mol": rec["energies_before"]["potential_energy_kj_per_mol"],
                        "energy_after_kj_per_mol": rec["energies_after"]["potential_energy_kj_per_mol"]})
        per_step.append(row)
        prev = st

    done_tags = {r["step"] for r in per_step}
    plan_all = plan
    if inc:
        plan = [s for s in plan if s["tag"] in done_tags]
    n_md = sum(1 for s in plan if s["kind"] != "minimisation")
    tot = {k: sum(s[k] for s in plan if s["kind"] == "production") for k in ("solute", "full", "state_log_rows")}
    tot_eq = {k: sum(s[k] for s in plan if s["kind"] == "equilibration") for k in ("solute", "full", "state_log_rows")}
    add(("every complete step" if inc else "every step has a record")
        + "; frames and log rows as expected (file size, DCD header, record)",
        f"{len(per_step)} of {len(plan_all)} steps complete; production {tot['solute']} solute frames, {tot['full']} full frames, "
        f"{tot['state_log_rows']} log rows; equilibration {tot_eq['solute']}, {tot_eq['full']}, {tot_eq['state_log_rows']}"
        + ("; PROBLEMS: " + "; ".join(bad_counts) if bad_counts else ""), not bad_counts and len(per_step) == len(plan))
    if not inc:
        add("manifest totals of the production",
            f"frames_solute {man['production']['frames_solute']}, frames_full {man['production']['frames_full']}, "
            f"completed_ns {man['production']['completed_ns']}, segments {man['production']['n_segments_complete']} of "
            f"{man['production']['n_segments_target']}",
            man["production"]["frames_solute"] == tot["solute"] and man["production"]["frames_full"] == tot["full"]
            and man["production"]["complete"])
    else:
        ok = True
        if cut_off:
            for name, info in cut_off["files"].items():
                if name.endswith(".dcd") and info["bytes"] > 0:
                    ok = ok and info.get("readable") and info["bytes_beyond_whole_frames"] == 0
        add("the cut-off step: its files hold whole frames only",
            {n: ({"bytes": i["bytes"], "whole_frames": i.get("whole_frames_by_size"), "header_frames": i.get("header_n_frames"),
                  "bytes_beyond": i.get("bytes_beyond_whole_frames")} if n.endswith(".dcd") else
                 {"bytes": i["bytes"], "rows": i.get("state_log_rows"), "last_time_ps": i.get("last_time_ps")})
             for n, i in (cut_off or {"files": {}})["files"].items()} or "no cut-off step", ok)
    add("time axis: log times and DCD headers continue from step to step",
        f"{n_md} steps" + ("; PROBLEMS: " + "; ".join(bad_time) if bad_time else "; none off"), not bad_time)
    add("MDAnalysis reads the same numbers of frames",
        "MDAnalysis not importable" if mda is None else ("; ".join(bad_mda) if bad_mda else "all files with a frame agree"),
        mda is not None and not bad_mda)
    add("sha256 and size of every file of every step equal the record",
        ("differ: " + ", ".join(bad_sha)) if bad_sha else f"{sum(len(C.read_json(os.path.join(run_dir, 'rec_' + s['tag'] + '.json'))['files']) for s in plan)} files",
        not bad_sha)
    add("each step started from the recorded state of the step before",
        ("broken at " + ", ".join(bad_chain)) if bad_chain else f"{n_md} links", not bad_chain)
    worst = max(last_frame_vs_state.values()) if last_frame_vs_state else None
    add("last solute frame of a step = its saved state (largest coordinate difference, A)",
        {k: round(v, 7) for k, v in last_frame_vs_state.items()}, worst is not None and worst < 1e-4)
    not_checked = []
    if full_vs_solute:
        add("solute atoms of each full frame = the solute frame of the same time (largest difference, A)",
            {k: round(v, 7) for k, v in full_vs_solute.items()}, max(full_vs_solute.values()) == 0.0)
    else:
        not_checked.append("solute atoms of each full frame against the solute frame of the same time: the run holds no "
                           "all-atom frame (no step reached the 100 ps interval of the all-atom file)")
    tags = sorted(seg_data)
    boundaries = []
    last_seg = ([s["tag"] for s in plan if s["kind"] == "production"] or [None])[-1]
    for a, b in zip(tags[:-1], tags[1:]):
        ha, xa = seg_data[a]
        hb, xb = seg_data[b]
        contiguous = hb["header_first_step"] == ha["header_first_step"] + len(xa) * ha["header_step_interval"]
        nearest = float(min(np.abs(xb[0] - f).max() for f in xa))
        boundaries.append({"between": f"{a}|{b}", "header_contiguous": contiguous,
                           "first_frame_differs_from_every_frame_before_by_at_least_A": round(nearest, 4)})
    if tot["solute"]:
        add("segment boundaries: headers contiguous, no frame repeated",
            boundaries if boundaries else "one segment only", all(b["header_contiguous"] and
            b["first_frame_differs_from_every_frame_before_by_at_least_A"] > 0.001 for b in boundaries))
    left = [os.path.basename(p) for p in glob.glob(os.path.join(run_dir, "*.tmp"))]
    if os.path.isfile(os.path.join(run_dir, "run_md.lock")):
        left.append("run_md.lock")
    chk = sorted(os.path.basename(p) for p in glob.glob(os.path.join(run_dir, "checkpoint_*.chk")))
    known = {f"rec_{s['tag']}.json" for s in plan}
    for s in plan:
        known |= set(C.read_json(os.path.join(run_dir, f"rec_{s['tag']}.json"))["files"]) if os.path.isfile(os.path.join(run_dir, f"rec_{s['tag']}.json")) else set()
    stray = [os.path.basename(p) for p in glob.glob(os.path.join(run_dir, "*")) if os.path.isfile(p)
             and step_of(os.path.basename(p)) is not None and os.path.basename(p) not in known and not p.endswith(".chk")]
    if not inc:
        add("nothing left behind: no .tmp file, no lock, only the newest checkpoint",
            f"left: {left or 'nothing'}; checkpoints: {chk}", not left and chk == [f"checkpoint_{last_seg}.chk"])
        add("no file of a step that the records do not list", stray or "none", not stray)
    else:
        expected_stray = sorted((cut_off or {"files": {}})["files"])
        add("what the kill left behind: the lock of the dead process and the files of the cut-off step, nothing else",
            f"left: {left or 'nothing'}; files without a record: {sorted(stray) or 'none'}; checkpoints: {chk or 'none'}",
            sorted(stray) == expected_stray and not [x for x in left if x != "run_md.lock"])

    invocations = [{"pid": i.get("pid"), "started_ist": i["started"]["ist"],
                    "ended_ist": (i.get("ended") or {}).get("ist") if isinstance(i.get("ended"), dict) else None,
                    "exit_code": i.get("exit_code"), "outcome": i.get("outcome"), "resume": i.get("resume"),
                    "threads": (i.get("platform") or {}).get("properties_set", {}).get("Threads"),
                    "platform": (i.get("platform") or {}).get("name")} for i in man["invocations"]]
    result = {"note": C.NOTE, "run_dir": run_dir, "run_id": man["run_id"], "checked": C.now_ist(),
              "inputs_sha256": {"manifest.json": C.sha256_file(os.path.join(run_dir, "manifest.json")),
                                "run_md.py (as recorded in the manifest)": man["script"]["sha256"],
                                "protocol (as recorded)": man["inputs"]["protocol"]["sha256"],
                                "system.prmtop (as recorded)": man["inputs"]["system.prmtop"]["sha256"],
                                "system.inpcrd (as recorded)": man["inputs"]["system.inpcrd"]["sha256"]},
              "run_seed": man["seed"]["run_seed"], "equil_scale": man["equil_scale"], "segment_ns": man["segment_ns"],
              "n_atoms": man["system"]["n_atoms"], "invocations": invocations, "steps": per_step,
              "cut_off_step": cut_off, "tests": tests, "not_checked": not_checked,
              "software": C.software(), "script_sha256": C.sha256_file(os.path.abspath(__file__))}

    if args.snapshot:
        snap = C.read_json(args.snapshot)
        itests = []
        killed = [i for i in man["invocations"] if i.get("pid") == snap["killed_pid"]]
        itests.append({"test": "the manifest records the killed invocation as one without an end",
                       "value": [{k: i.get(k) for k in ("pid", "ended", "exit_code", "outcome")} for i in killed],
                       "passed": len(killed) == 1 and killed[0].get("ended") is None})
        unchanged, changed = 0, []
        for name, info in snap["files"].items():
            tag = step_of(name)
            if tag is None or tag not in snap["steps_with_a_record"]:
                continue
            p = os.path.join(run_dir, name)
            if name.endswith(".chk") and not os.path.isfile(p):
                continue                                      # an older checkpoint is removed by design
            if os.path.isfile(p) and C.sha256_file(p) == info["sha256"]:
                unchanged += 1
            else:
                changed.append(name)
        itests.append({"test": "files of the steps complete before the kill are unchanged (sha256)",
                       "value": f"{unchanged} unchanged" + (f"; CHANGED: {changed}" if changed else ""), "passed": not changed})
        for tag, entry in snap["steps_without_a_record"].items():
            st = next(s for s in plan if s["tag"] == tag)
            rec = C.read_json(os.path.join(run_dir, f"rec_{tag}.json"))
            before = {n: (i.get("whole_frames_by_size") if n.endswith(".dcd") else i.get("state_log_rows"))
                      for n, i in entry["files"].items()}
            hs = C.read_dcd_header(os.path.join(run_dir, f"solute_{tag}.dcd"))
            hf = C.read_dcd_header(os.path.join(run_dir, f"full_{tag}.dcd"))
            _hd, log = C.read_log_csv(os.path.join(run_dir, f"log_{tag}.csv"))
            after = {f"solute_{tag}.dcd": hs["whole_frames_by_size"], f"full_{tag}.dcd": hf["whole_frames_by_size"],
                     f"log_{tag}.csv": len(log)}
            itests.append({"test": f"{tag}: frames and rows left by the kill, and after --resume",
                           "value": {"left_by_the_kill": before, "after_resume": after,
                                     "expected": {f"solute_{tag}.dcd": st["solute"], f"full_{tag}.dcd": st["full"],
                                                  f"log_{tag}.csv": st["state_log_rows"]},
                                     "step_started_again_ist": rec["started"]["ist"], "killed_at_ist": snap["taken"]},
                           "passed": after == {f"solute_{tag}.dcd": st["solute"], f"full_{tag}.dcd": st["full"],
                                               f"log_{tag}.csv": st["state_log_rows"]}})
            itests.append({"test": f"{tag}: started from the state recorded before the kill",
                           "value": {"state_in": rec["state_in"]["file"], "sha256_matches_snapshot":
                                     rec["state_in"]["sha256"] == snap["files"][rec["state_in"]["file"]]["sha256"]},
                           "passed": rec["state_in"]["sha256"] == snap["files"][rec["state_in"]["file"]]["sha256"]})
        result["interruption"] = {"snapshot_file": os.path.abspath(args.snapshot), "killed_pid": snap["killed_pid"],
                                  "signal": snap["signal"], "killed_at": snap["taken"],
                                  "lock_file_left_by_the_kill": snap["lock_file_left"], "tests": itests}
        tests.extend(itests)

    result["all_passed"] = all(t["passed"] for t in tests)
    C.write_json(args.out, result)
    for t in tests:
        print(("passed  " if t["passed"] else "FAILED  ") + t["test"] + ": " + str(t["value"])[:400])
    for t in not_checked:
        print("NOT CHECKED  " + t)
    print(("ALL PASSED" if result["all_passed"] else "SOME TESTS FAILED")
          + (f" ({len(not_checked)} not checked)" if not_checked else "") + f"; written to {args.out}")
    return 0 if result["all_passed"] else 1


def recovery(args):
    snap = C.read_json(args.snapshot)
    out = {"note": C.NOTE, "what": "copies of the killed run folder after 'run_md.py --stage min --resume': the driver's "
           "recovery up to the point where it would repeat the cut-off segment (the repeat itself was not run)",
           "snapshot_file": os.path.abspath(args.snapshot), "checked": C.now_ist(), "copies": {}, "software": C.software(),
           "script_sha256": C.sha256_file(os.path.abspath(__file__))}
    all_ok = True
    for copy in sorted(glob.glob(os.path.join(args.dir, "*", "X1"))):
        label = os.path.basename(os.path.dirname(copy))
        man = C.read_json(os.path.join(copy, "manifest.json"))
        inv = man["invocations"]
        killed = [i for i in inv if i.get("pid") == snap["killed_pid"]]
        last = inv[-1]
        cut = sorted(snap["steps_without_a_record"])
        gone = [n for tag in cut for n in snap["steps_without_a_record"][tag]["files"] if not os.path.isfile(os.path.join(copy, n))]
        still = [n for tag in cut for n in snap["steps_without_a_record"][tag]["files"] if os.path.isfile(os.path.join(copy, n))]
        same, changed = 0, []
        for name, info in snap["files"].items():
            p = os.path.join(copy, name)
            if step_of(name) in snap["steps_with_a_record"] and os.path.isfile(p):
                if C.sha256_file(p) == info["sha256"]:
                    same += 1
                else:
                    changed.append(name)
        log = os.path.join(copy, "run_md.log")
        lines = open(log).read().split("\n") if os.path.isfile(log) else []
        tests = [
            {"test": "exit code 0 and nothing run", "value": {"exit_code": last.get("exit_code"), "outcome": last.get("outcome"),
                                                              "status": man.get("status")},
             "passed": last.get("exit_code") == 0 and man.get("status") == "incomplete"},
            {"test": "the killed invocation is marked as one without an end",
             "value": [{k: i.get(k) for k in ("pid", "ended", "exit_code", "outcome")} for i in killed],
             "passed": len(killed) == 1 and killed[0].get("ended") is None and "killed" in str(killed[0].get("outcome"))},
            {"test": "the lock of the dead process was taken and released",
             "value": {"lock_left_by_the_kill": snap["lock_file_left"], "lock_now": os.path.isfile(os.path.join(copy, "run_md.lock"))},
             "passed": snap["lock_file_left"] is not None and not os.path.isfile(os.path.join(copy, "run_md.lock"))},
            {"test": "the files of the cut-off step were removed", "value": {"removed": gone, "still_there": still},
             "passed": bool(gone) and not still},
            {"test": "the files of the complete steps present in the copy are unchanged (sha256)",
             "value": f"{same} unchanged" + (f"; CHANGED {changed}" if changed else ""), "passed": same > 0 and not changed},
            {"test": "the driver names the step it would continue with",
             "value": [ln for ln in lines if "continuing with" in ln or "removed the incomplete" in ln or "WARNING" in ln or "NOTE" in ln],
             "passed": any("continuing with " + cut[0] in ln for ln in lines) if cut else False},
        ]
        out["copies"][label] = {"folder": copy, "driver_sha256": last.get("script_sha256"),
                                "driver_of_the_run_sha256": man["script"]["sha256"], "tests": tests,
                                "all_passed": all(t["passed"] for t in tests)}
        all_ok = all_ok and out["copies"][label]["all_passed"]
        for t in tests:
            print(label, "passed " if t["passed"] else "FAILED ", t["test"] + ":", str(t["value"])[:300])
    out["all_passed"] = all_ok and bool(out["copies"])
    C.write_json(args.out, out)
    print(("ALL PASSED" if out["all_passed"] else "SOME TESTS FAILED") + f"; written to {args.out}")
    return 0 if out["all_passed"] else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot")
    s.add_argument("run_dir")
    s.add_argument("--killed-pid", type=int, required=True)
    s.add_argument("--signal", default="SIGKILL")
    s.add_argument("--out", required=True)
    c = sub.add_parser("check")
    c.add_argument("run_dir")
    c.add_argument("--snapshot")
    c.add_argument("--incomplete", action="store_true")
    c.add_argument("--out", required=True)
    r = sub.add_parser("recovery")
    r.add_argument("dir")
    r.add_argument("--snapshot", required=True)
    r.add_argument("--out", required=True)
    args = ap.parse_args()
    return {"snapshot": snapshot, "check": check, "recovery": recovery}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
