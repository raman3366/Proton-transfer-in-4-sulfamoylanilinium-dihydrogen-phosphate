#!/usr/bin/env python3
"""queue.py: a small task runner for the rented GPU machine (stage S5 of md/PLAN.md; how to use: md/gpu/GPU_RUN.md).

Python 3.8 or newer, standard library only, Linux only (it reads /proc). No project code is imported.

    queue.py run      [--root DIR] [--cpus N] [--gpus N]     the runner (started at boot by cqc-md-queue.service)
    queue.py status   [--root DIR] [--json]                  table of the tasks, age of the heartbeat
    queue.py check    FILE...  [--root DIR] [--cpus N]       validate task files, also against the tasks already there
    queue.py add      FILE...  [--root DIR] [--replace]      validate, then copy into the task folder (atomic)
    queue.py retry    ID...    [--root DIR]                  give a failed or cancelled task new attempts
    queue.py cancel   ID...    [--root DIR]                  end a running task and do not start it again
    queue.py stop | drain | resume  [--root DIR]             write or remove the STOP and DRAIN files
    queue.py json     FILE KEY=VALUE ...                     exit code 0 if every condition holds in a JSON file
                                                             (for success tests; see `json` below)
    queue.py selftest [--dir DIR] [--keep]                   dummy tasks (sleep, echo): the tests of GPU_RUN.md

WHAT IT DOES
  It reads every *.json file of <root>/queue/tasks.d (again every few seconds, so tasks can be ADDED while it
  runs), starts each task when the tasks it depends on are done and the machine has room, at most ONE task
  with resource "gpu" at a time and CPU tasks beside it within the CPUs of the machine, and writes what
  happens to <root>/status/queue.json (every change of state with UTC and IST time), <root>/status/queue.log
  (one line per event), <root>/status/heartbeat.json (every minute) and <root>/logs/<id>.log (output of the
  task, appended).

WHEN A TASK IS DONE
  Only when its success test holds: the file `done_marker` exists and/or the shell command `success_test`
  exits with 0. A task whose test holds already is marked done WITHOUT being run. A task that ran is done
  when its command ended with exit code 0 AND the test holds. At every start of the runner the test of every
  task recorded as done is made again.

INTERRUPTIONS
  The runner assumes that every task can be started again (run_md.py --resume, run_qm.py from its
  checkpoints). A task that was running when the runner or the machine went down is started again. A task
  whose command fails is tried `max_attempts` times (default 3) and is then "failed": it and everything that
  depends on it wait for a person (`queue.py retry ID`).

STOPPING
  <root>/queue/STOP   running tasks get SIGTERM (SIGKILL after `stop_grace_s`), the state is written, the
                      runner ends with exit code 0. While the file exists the runner starts nothing, also
                      not after a reboot. Remove it and start the service to go on.
  <root>/queue/DRAIN  nothing new is started; what runs is left to finish; the runner stays alive.
  SIGTERM, SIGINT     as STOP, without the file (what systemd sends at shutdown).

TASK FILE (JSON): {"env": {...}, "tasks": [ {...}, ... ]}  or just the list
  id            letters, digits, _ . - ; unique over all files (the first file in name order wins)
  command       text, or a list of lines; run by  bash -eo pipefail -c  in `cwd`
  cwd           working directory (absolute; created if absent); $VARS of the task's environment are expanded
  resources     "gpu" (the GPU and one CPU), "cpu:N" (N CPUs), or "gpu+cpu:N"
  depends_on    list of ids that must be done first (default none)
  done_marker   a file that exists when the task is done, or "auto" (the runner writes its own marker when
                the command ends with exit code 0)
  success_test  shell command, run in `cwd`; exit code 0 means done. At least one of the two is required.
  optional      title; env {NAME: value}; nice (default 0 for gpu tasks, 10 for cpu tasks);
                max_attempts (3); max_interruptions (30); retry_delay_s (60); stop_grace_s (90);
                timeout_s (none: wall-time limit of one attempt); test_timeout_s (120)
  Every task gets OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and NUMEXPR_NUM_THREADS set to its
  number of CPUs unless its env names them, and QUEUE_TASK_ID, QUEUE_ATTEMPT, QUEUE_CPUS, QUEUE_ROOT.

JSON TESTS:  queue.py json manifest.json status=complete production.completed_ns=100 test_only=false
  KEY is a path of names and list indices joined by dots (jobs.opt.state, invocations.-1.exit_code); VALUE is
  read as JSON (true, 100, "text") and otherwise as text; numbers are compared as numbers. KEY!=VALUE asks
  for a difference. Exit code 0: all hold; 1: at least one does not (it is printed); 2: file or key missing.

Exit codes of `run`: 0 stopped as asked (STOP, signal, --exit-when-idle); 2 wrong call; 3 another runner holds
the lock of this root; 1 anything else.  This program prints no secrets and reads none.
"""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
import pwd
import re
import shutil
import signal
import socket
import string
import subprocess
import sys
import tempfile
import time
import traceback
import uuid

SCHEMA = "cqc-md queue state, version 1"
ROOT_DEFAULT = "/root/cqc-sfa-md"
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), "IST")
TOKEN_VAR = "CQC_QUEUE_TOKEN"
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
DEFAULTS = {"max_attempts": 3, "max_interruptions": 30, "retry_delay_s": 60, "stop_grace_s": 90,
            "timeout_s": None, "test_timeout_s": 120}
KNOWN_KEYS = {"id", "title", "command", "cwd", "resources", "depends_on", "done_marker", "success_test", "env",
              "nice", *DEFAULTS}
HISTORY_MAX = 200


class TaskError(Exception):
    pass


# ------------------------------------------------------------------------------------------------ small tools
def stamp(t=None):
    now = (datetime.datetime.fromtimestamp(t, datetime.timezone.utc) if t is not None
           else datetime.datetime.now(datetime.timezone.utc))
    ist = now.astimezone(IST)
    hour = ist.hour % 12 or 12
    return {"utc": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "ist": f"{ist.day} {ist.strftime('%b %Y')}, {hour}:{ist.minute:02d}:{ist.second:02d} "
                   f"{'AM' if ist.hour < 12 else 'PM'} IST"}


def write_json_atomic(path, obj):
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def read_json(path, default=None):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return default


def boot_id():
    try:
        with open("/proc/sys/kernel/random/boot_id") as f:
            return f.read().strip()
    except OSError:
        return None


def token_pids(token):
    """PIDs of all living processes that carry the token in their environment (a task and its children)."""
    needle = f"{TOKEN_VAR}={token}".encode()
    out = []
    for name in os.listdir("/proc"):
        if not name.isdigit() or int(name) == os.getpid():
            continue
        try:
            with open(f"/proc/{name}/environ", "rb") as f:
                if needle in f.read().split(b"\0"):
                    out.append(int(name))
        except OSError:
            continue
    return sorted(out)


def send(pids, sig):
    for p in pids:
        try:
            os.kill(p, sig)
        except OSError:
            pass


class Paths:
    def __init__(self, root):
        self.root = os.path.abspath(root)
        self.queue = os.path.join(self.root, "queue")
        self.tasks = os.path.join(self.queue, "tasks.d")
        self.control = os.path.join(self.queue, "control")
        self.done = os.path.join(self.queue, "done")
        self.stop = os.path.join(self.queue, "STOP")
        self.drain = os.path.join(self.queue, "DRAIN")
        self.lock = os.path.join(self.queue, "queue.lock")
        self.status = os.path.join(self.root, "status")
        self.state = os.path.join(self.status, "queue.json")
        self.heartbeat = os.path.join(self.status, "heartbeat.json")
        self.events = os.path.join(self.status, "queue.log")
        self.logs = os.path.join(self.root, "logs")

    def make(self):
        for d in (self.tasks, self.control, self.done, self.status, self.logs):
            os.makedirs(d, exist_ok=True)


# ------------------------------------------------------------------------------------------------ task files
def parse_resources(text):
    """'gpu' -> (1, 1); 'cpu:N' -> (0, N); 'gpu+cpu:N' -> (1, N)."""
    if not isinstance(text, str) or not text.strip():
        raise TaskError('resources must be "gpu", "cpu:N" or "gpu+cpu:N"')
    gpu, cpu = 0, None
    for tok in re.split(r"[+,\s]+", text.strip()):
        if tok == "gpu":
            gpu = 1
        elif re.fullmatch(r"cpu:[0-9]+", tok):
            cpu = int(tok[4:])
        else:
            raise TaskError(f'resources "{text}": expected "gpu", "cpu:N" or "gpu+cpu:N"')
    if cpu is None:
        if not gpu:
            raise TaskError(f'resources "{text}" names neither gpu nor cpu:N')
        cpu = 1
    if cpu < 1:
        raise TaskError(f'resources "{text}": at least one CPU')
    return gpu, cpu


def normalise(raw, file_env, source):
    """One task of a file -> the definition the runner works with. Raises TaskError."""
    if not isinstance(raw, dict):
        raise TaskError("a task must be a JSON object")
    tid = raw.get("id")
    if not isinstance(tid, str) or not ID_RE.match(tid):
        raise TaskError(f"bad or missing id: {tid!r}")
    unknown = sorted(set(raw) - KNOWN_KEYS)
    if unknown:
        raise TaskError(f"task {tid}: unknown key(s) {', '.join(unknown)}")
    cmd = raw.get("command")
    if isinstance(cmd, list) and cmd and all(isinstance(x, str) for x in cmd):
        cmd = "\n".join(cmd)
    if not isinstance(cmd, str) or not cmd.strip():
        raise TaskError(f"task {tid}: command must be text or a list of lines")
    cwd = raw.get("cwd")
    if not isinstance(cwd, str) or not cwd.startswith(("/", "$")):
        raise TaskError(f"task {tid}: cwd must be an absolute path")
    try:
        gpu, cpu = parse_resources(raw.get("resources"))
    except TaskError as e:
        raise TaskError(f"task {tid}: {e}")
    deps = raw.get("depends_on", [])
    if not isinstance(deps, list) or not all(isinstance(x, str) and ID_RE.match(x) for x in deps):
        raise TaskError(f"task {tid}: depends_on must be a list of ids")
    if tid in deps:
        raise TaskError(f"task {tid}: depends on itself")
    marker, test = raw.get("done_marker"), raw.get("success_test")
    if isinstance(test, list) and test and all(isinstance(x, str) for x in test):
        test = "\n".join(test)
    if marker is None and test is None:
        raise TaskError(f'task {tid}: needs done_marker (a path or "auto") or success_test')
    if marker is not None and (not isinstance(marker, str) or not (marker == "auto" or marker.startswith(("/", "$")))):
        raise TaskError(f'task {tid}: done_marker must be an absolute path or "auto"')
    if test is not None and (not isinstance(test, str) or not test.strip()):
        raise TaskError(f"task {tid}: success_test must be a shell command")
    env = dict(file_env)
    tenv = raw.get("env", {})
    if not isinstance(tenv, dict):
        raise TaskError(f"task {tid}: env must be an object")
    env.update(tenv)
    for k, v in env.items():
        if not isinstance(k, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", k) or not isinstance(v, (str, int, float)):
            raise TaskError(f"task {tid}: env entry {k!r} must be NAME: text or number")
    env = {k: str(v) for k, v in env.items()}
    d = {"id": tid, "title": str(raw.get("title", "")), "command": cmd, "cwd": cwd, "resources": raw["resources"],
         "gpu": gpu, "cpu": cpu, "depends_on": list(deps), "done_marker": marker, "success_test": test, "env": env,
         "nice": raw.get("nice", 0 if gpu else 10), "source": source}
    for k, default in DEFAULTS.items():
        v = raw.get(k, default)
        if v is not None and (isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0):
            raise TaskError(f"task {tid}: {k} must be a number")
        d[k] = v
    if isinstance(d["nice"], bool) or not isinstance(d["nice"], int) or not 0 <= d["nice"] <= 19:
        raise TaskError(f"task {tid}: nice must be 0..19")
    if d["max_attempts"] < 1:
        raise TaskError(f"task {tid}: max_attempts must be at least 1")
    core = {k: d[k] for k in ("command", "cwd", "resources", "depends_on", "done_marker", "success_test", "env",
                              "nice", *DEFAULTS)}
    d["def_sha256"] = hashlib.sha256(json.dumps(core, sort_keys=True).encode()).hexdigest()
    return d


def read_task_file(path):
    """-> (list of definitions, list of error texts). A file that is not JSON gives no tasks and one error."""
    name = os.path.basename(path)
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        return None, [f"{name}: cannot be read as JSON ({e})"]
    if isinstance(data, list):
        file_env, raw_tasks = {}, data
    elif isinstance(data, dict) and isinstance(data.get("tasks"), list):
        file_env, raw_tasks = data.get("env", {}), data["tasks"]
        extra = sorted(set(data) - {"env", "tasks", "comment"})
        if extra:
            return None, [f"{name}: unknown top-level key(s) {', '.join(extra)}"]
        if not isinstance(file_env, dict):
            return None, [f"{name}: env must be an object"]
    else:
        return None, [f'{name}: expected a list of tasks or {{"env": ..., "tasks": [...]}}']
    out, errors = [], []
    for raw in raw_tasks:
        try:
            out.append(normalise(raw, file_env, name))
        except TaskError as e:
            errors.append(f"{name}: {e}")
    return out, errors


def list_task_files(folder):
    try:
        names = sorted(n for n in os.listdir(folder) if n.endswith(".json") and not n.startswith("."))
    except OSError:
        return []
    return [os.path.join(folder, n) for n in names]


def merge_files(per_file):
    """per_file: list of (file name, definitions). -> (ordered definitions, errors); the first id wins."""
    seen, out, errors = {}, [], []
    for name, defs in per_file:
        for d in defs:
            if d["id"] in seen:
                errors.append(f"{name}: id {d['id']} is already defined in {seen[d['id']]}; this one is ignored")
                continue
            seen[d["id"]] = name
            out.append(d)
    return out, errors


def structural_problems(defs, cpus, gpus):
    """-> (errors that make a task unrunnable {id: text}, warnings [text])."""
    by_id = {d["id"]: d for d in defs}
    bad, warn = {}, []
    for d in defs:
        if d["cpu"] > cpus:
            bad[d["id"]] = f"needs {d['cpu']} CPUs, the machine offers {cpus}"
        elif d["gpu"] > gpus:
            bad[d["id"]] = f"needs a GPU, the runner was started with --gpus {gpus}"
        for dep in d["depends_on"]:
            if dep not in by_id:
                warn.append(f"task {d['id']} waits for {dep}, which is not in the list (yet)")
    # cycles: depth-first search
    colour = {}

    def visit(tid, stack):
        colour[tid] = 1
        for dep in by_id[tid]["depends_on"]:
            if dep not in by_id:
                continue
            if colour.get(dep) == 1:
                cyc = stack[stack.index(dep):] + [dep] if dep in stack else [tid, dep]
                for t in cyc:
                    bad.setdefault(t, "depends on itself through " + " -> ".join(cyc))
            elif dep not in colour:
                visit(dep, stack + [dep])
        colour[tid] = 2

    for d in defs:
        if d["id"] not in colour:
            visit(d["id"], [d["id"]])
    return bad, warn


# ------------------------------------------------------------------------------------------------ the runner
class Runner:
    def __init__(self, args):
        self.a = args
        self.p = Paths(args.root)
        self.p.make()
        self.cpus = args.cpus or os.cpu_count() or 1
        self.gpus = args.gpus
        self.t = {}                  # id -> record (definition, state, history)
        self.order = []
        self.procs = {}              # id -> Popen
        self.prechecked = set()
        self.file_cache = {}         # file -> (signature, definitions)
        self.file_errors, self.warnings = [], []
        self.stop_signal = None
        self.started = stamp()
        self.idle_since = None
        self.last_heartbeat = 0.0
        self.last_reload = 0.0
        self.dirty = True
        self.lock_fd = None
        self.qstate = "starting"
        self.boot_enabled = None
        self.removed = {}            # id -> last state of tasks that left the task files
        self.gpu_last = None

    # ---- bookkeeping
    def event(self, text):
        st = stamp()
        line = f"{st['utc']} | {st['ist']} | {text}"
        print(line, flush=True)
        try:
            with open(self.p.events, "a") as f:
                f.write(line + "\n")
        except OSError as e:
            print(f"WARNING: cannot append to {self.p.events}: {e}", file=sys.stderr, flush=True)

    def set_state(self, tid, new, note=""):
        rec = self.t[tid]
        old = rec["state"]
        st = stamp()
        rec.update(state=new, state_since=st, note=note)
        rec["history"].append({"utc": st["utc"], "ist": st["ist"], "from": old, "to": new, "note": note})
        del rec["history"][:-HISTORY_MAX]
        self.event(f"task {tid}: {old} -> {new}" + (f" ({note})" if note else ""))
        self.dirty = True

    def counts(self):
        c = {}
        for tid in self.order:
            c[self.t[tid]["state"]] = c.get(self.t[tid]["state"], 0) + 1
        return c

    def in_use(self):
        cpu = sum(self.t[t]["def"]["cpu"] for t in self.procs)
        gpu = sum(self.t[t]["def"]["gpu"] for t in self.procs)
        return cpu, gpu

    def attention(self):
        out = []
        for tid in self.order:
            rec = self.t[tid]
            if rec["state"] in ("failed", "cancelled"):
                out.append(f"{tid}: {rec['state']} ({rec.get('note', '')})")
            elif rec["state"] == "pending" and rec.get("blocked_by"):
                out.append(f"{tid}: waits for {', '.join(rec['blocked_by'])}")
        return out + [f"task file: {e}" for e in self.file_errors]

    def write_state(self, heartbeat=False):
        st = stamp()
        cpu, gpu = self.in_use()
        running = []
        for tid, proc in self.procs.items():
            rec = self.t[tid]
            running.append({"id": tid, "pid": proc.pid, "resources": rec["def"]["resources"],
                            "since_utc": rec["started"]["utc"], "since_ist": rec["started"]["ist"],
                            "elapsed_s": round(time.time() - rec["started_epoch"])})
        counts = self.counts()
        all_done = bool(self.order) and all(self.t[t]["state"] == "done" for t in self.order)
        q = {"state": self.qstate, "pid": os.getpid(), "host": socket.gethostname(), "boot_id": boot_id(),
             "started_utc": self.started["utc"], "started_ist": self.started["ist"],
             "heartbeat_utc": st["utc"], "heartbeat_ist": st["ist"], "heartbeat_epoch": round(time.time()),
             "started_by": "systemd" if os.getppid() == 1 and os.environ.get("INVOCATION_ID") else "shell",
             "starts_at_boot": self.boot_enabled,
             "cpus_total": self.cpus, "gpus_total": self.gpus, "cpus_in_use": cpu, "gpus_in_use": gpu,
             "task_folder": self.p.tasks, "stop_file": self.p.stop, "stop_file_present": os.path.exists(self.p.stop),
             "drain_file_present": os.path.exists(self.p.drain),
             "idle": not self.procs and not counts.get("retry"), "all_done": all_done,
             "idle_since_utc": self.idle_since["utc"] if self.idle_since else None,
             "idle_since_ist": self.idle_since["ist"] if self.idle_since else None,
             "task_file_errors": self.file_errors, "warnings": self.warnings}
        state = {"schema": SCHEMA, "queue": q, "summary": counts, "running": running, "order": self.order,
                 "tasks": {tid: {k: v for k, v in self.t[tid].items() if k not in ("started_epoch", "not_before")}
                           for tid in self.order},
                 "removed": self.removed}
        try:
            write_json_atomic(self.p.state, state)
            if heartbeat:
                hb = {"schema": "cqc-md queue heartbeat, version 1", "utc": st["utc"], "ist": st["ist"],
                      "epoch": round(time.time()), "host": q["host"], "pid": q["pid"], "boot_id": q["boot_id"],
                      "queue_state": self.qstate, "running": running, "counts": counts, "n_tasks": len(self.order),
                      "idle": q["idle"], "all_done": all_done, "idle_since_utc": q["idle_since_utc"],
                      "attention": self.attention(), "cpus_total": self.cpus, "cpus_in_use": cpu,
                      "gpus_in_use": gpu, "load_average": [round(x, 2) for x in os.getloadavg()],
                      "disk_free_gb": round(shutil.disk_usage(self.p.root).free / 1e9, 1),
                      "uptime_s": uptime_s(), "gpu": self.gpu_last, "starts_at_boot": self.boot_enabled}
                write_json_atomic(self.p.heartbeat, hb)
            self.dirty = False
        except OSError as e:
            print(f"WARNING: state not written: {e}", file=sys.stderr, flush=True)

    # ---- the task list
    def reload(self):
        self.last_reload = time.time()
        per_file, errors = [], []
        for path in list_task_files(self.p.tasks):
            name = os.path.basename(path)
            try:
                s = os.stat(path)
                sig = (s.st_mtime_ns, s.st_size, s.st_ino)
            except OSError:
                continue
            cached = self.file_cache.get(name)
            if cached and cached[0] == sig:
                defs, errs = cached[1], cached[2]
            else:
                defs, errs = read_task_file(path)
                if defs is None:                    # not readable: keep what this file gave the last time
                    defs = cached[1] if cached else []
                    errs = errs + ([f"{name}: its last readable version stays in force"] if cached else [])
                    if not cached or cached[2] != errs:
                        self.event(f"task file problem: {'; '.join(errs)}")
                elif not cached or cached[1] != defs:
                    self.event(f"task file {name}: {len(defs)} task(s) read"
                               + (f"; problems: {'; '.join(errs)}" if errs else ""))
                self.file_cache[name] = (sig, defs, errs)
            per_file.append((name, defs))
            errors += errs
        defs, dup = merge_files(per_file)
        errors += dup
        bad, warn = structural_problems(defs, self.cpus, self.gpus)
        if errors != self.file_errors or warn != self.warnings:
            self.file_errors, self.warnings, self.dirty = errors, warn, True
        new_order = [d["id"] for d in defs]
        for d in defs:
            tid = d["id"]
            rec = self.t.get(tid)
            if rec is None:
                self.t[tid] = {"def": d, "state": "new", "attempts": 0, "interruptions": 0, "history": [],
                               "log": os.path.join(self.p.logs, f"{tid}.log"), "note": ""}
                self.set_state(tid, "pending", f"read from {d['source']}")
            elif rec["def"]["def_sha256"] != d["def_sha256"]:
                if tid in self.procs:
                    if rec.get("next_def", {}).get("def_sha256") != d["def_sha256"]:
                        rec["next_def"] = d
                        self.event(f"task {tid}: its definition changed while it runs; the new one is used at its next start")
                        self.dirty = True
                else:
                    rec["def"] = d
                    rec.pop("next_def", None)
                    rec.update(attempts=0, interruptions=0)
                    self.prechecked.discard(tid)
                    self.set_state(tid, "pending", "definition changed; attempts counted from zero")
            else:
                rec["def"] = d           # title or source may have changed
                if tid not in self.procs:
                    rec.pop("next_def", None)
            if tid in bad and self.t[tid]["state"] in ("pending", "retry"):
                self.set_state(tid, "failed", bad[tid])
        gone = [tid for tid in self.order if tid not in new_order]
        for tid in gone:
            if tid in self.procs:
                new_order.append(tid)     # stays until it has ended
                if not self.t[tid].get("removed_from_list"):
                    self.t[tid]["removed_from_list"] = True
                    self.event(f"task {tid}: no longer in the task files; it runs on and is dropped when it ends")
            else:
                self.removed[tid] = {"state": self.t[tid]["state"], "removed": stamp()}
                self.event(f"task {tid}: no longer in the task files, dropped (was {self.t[tid]['state']})")
                del self.t[tid]
                self.dirty = True
        if new_order != self.order:
            self.order, self.dirty = new_order, True

    def restore(self):
        """Take over the states of the last run of this root, deal with tasks that were running, read the task
        files, and make the success test of every task recorded as done again."""
        old = read_json(self.p.state)
        try:
            self.take_over(old)
        except Exception as e:                                  # noqa: BLE001
            self.t, self.order, self.removed = {}, [], {}
            self.event(f"the state of the last run could not be taken over ({type(e).__name__}: {e}); starting from the task files")
        self.reload()
        for tid in self.order:
            rec = self.t[tid]
            if rec["state"] == "done":
                ok, detail = self.check_success(rec["def"])
                rec["last_test"] = {"ok": ok, "detail": detail, **stamp()}
                if ok:
                    self.prechecked.add(tid)
                else:
                    rec["attempts"] = 0
                    self.set_state(tid, "pending", f"was done, but its success test does not hold any more: {detail}")

    def take_over(self, old):
        if old and old.get("schema") == SCHEMA:
            self.removed = old.get("removed", {})
            same_boot = old.get("queue", {}).get("boot_id") == boot_id()
            for tid in old.get("order", []):
                o = old.get("tasks", {}).get(tid)
                if not o or "def" not in o or "state" not in o:
                    continue
                rec = dict(o)
                for k in ("pid", "next_def", "blocked_by", "waiting_for", "removed_from_list"):
                    rec.pop(k, None)
                rec.setdefault("history", [])
                rec.setdefault("attempts", 0)
                rec.setdefault("interruptions", 0)
                rec["log"] = os.path.join(self.p.logs, f"{tid}.log")
                self.t[tid] = rec
                self.order.append(tid)
            for tid in self.order:
                rec = self.t[tid]
                if rec["state"] == "running":
                    left = token_pids(rec["token"]) if same_boot and rec.get("token") else []
                    if left:
                        self.event(f"task {tid}: {len(left)} process(es) of the last runner are still alive "
                                   f"({' '.join(map(str, left))}); ending them before the task is started again")
                        self.event(f"task {tid}: {self.end_processes(rec['token'], rec['def']['stop_grace_s'])}")
                    rec["interruptions"] += 1
                    self.set_state(tid, "interrupted", "was running when the runner or the machine went down"
                                   + ("" if same_boot else "; the machine was restarted"))
                    if rec["interruptions"] > rec["def"]["max_interruptions"]:
                        self.set_state(tid, "failed", f"interrupted {rec['interruptions']} times, more than max_interruptions")
                    else:
                        self.set_state(tid, "pending", "to be started again")
                elif rec["state"] in ("retry", "interrupted"):
                    self.set_state(tid, "pending", "the runner was started again")

    # ---- environment, tests, processes
    def task_env(self, d, attempt=0):
        env = dict(os.environ)
        env.setdefault("HOME", pwd.getpwuid(os.getuid()).pw_dir)
        env.setdefault("PYTHONUNBUFFERED", "1")
        for v in THREAD_VARS:
            env[v] = str(d["cpu"])
        env.update(d["env"])
        env.update(QUEUE_TASK_ID=d["id"], QUEUE_ATTEMPT=str(attempt), QUEUE_CPUS=str(d["cpu"]), QUEUE_ROOT=self.p.root)
        env.pop(TOKEN_VAR, None)
        return env

    def marker_path(self, d, env):
        if d["done_marker"] == "auto":
            return os.path.join(self.p.done, f"{d['id']}.done")
        return string.Template(d["done_marker"]).safe_substitute(env)

    def check_success(self, d):
        env = self.task_env(d)
        if d["done_marker"]:
            path = self.marker_path(d, env)
            if not os.path.exists(path):
                return False, f"no file {path}"
        if d["success_test"]:
            cwd = string.Template(d["cwd"]).safe_substitute(env)
            try:
                r = subprocess.run(["/bin/bash", "-c", d["success_test"]], cwd=cwd if os.path.isdir(cwd) else "/",
                                   env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   timeout=d["test_timeout_s"])
            except subprocess.TimeoutExpired:
                return False, f"success_test did not end in {d['test_timeout_s']} s"
            except OSError as e:
                return False, f"success_test could not be run: {e}"
            if r.returncode != 0:
                tail = r.stdout.decode(errors="replace").strip().splitlines()[-1:] or [""]
                return False, f"success_test exit code {r.returncode}" + (f": {tail[0][:200]}" if tail[0] else "")
        return True, "holds"

    def end_processes(self, token, grace):
        """SIGTERM to every process of a task, SIGKILL to what is left after `grace` seconds."""
        send(token_pids(token), signal.SIGTERM)
        deadline = time.time() + grace
        while time.time() < deadline:
            for proc in self.procs.values():
                proc.poll()
            if not token_pids(token):
                return "ended on SIGTERM"
            time.sleep(0.2)
        left = token_pids(token)
        send(left, signal.SIGKILL)
        for _ in range(50):
            for proc in self.procs.values():
                proc.poll()
            if not token_pids(token):
                break
            time.sleep(0.2)
        return f"SIGKILL after {grace:g} s to {len(left)} process(es)"

    def start(self, tid):
        rec = self.t[tid]
        if rec.get("next_def"):
            rec["def"] = rec.pop("next_def")
        d = rec["def"]
        attempt = rec["attempts"] + 1
        env = self.task_env(d, attempt)
        token = uuid.uuid4().hex
        env[TOKEN_VAR] = token
        cwd = string.Template(d["cwd"]).safe_substitute(env)
        st = stamp()
        try:
            os.makedirs(cwd, exist_ok=True)
            with open(rec["log"], "ab") as logf:
                logf.write((f"\n===== {st['utc']} | {st['ist']} | task {tid}: attempt {attempt} of {d['max_attempts']}"
                            f", interrupted {rec['interruptions']} time(s) before | resources {d['resources']} "
                            f"| cwd {cwd} =====\n").encode())
                logf.flush()
                proc = subprocess.Popen(["nice", "-n", str(d["nice"]), "/bin/bash", "-eo", "pipefail", "-c", d["command"]],
                                        cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=logf,
                                        stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as e:
            rec["attempts"] = attempt
            self.fail(tid, f"could not be started: {e}")
            return
        self.procs[tid] = proc
        rec.update(pid=proc.pid, token=token, started=st, started_epoch=time.time(), attempt_running=attempt)
        rec.pop("ended", None)
        rec.pop("exit_code", None)
        self.set_state(tid, "running", f"attempt {attempt}, pid {proc.pid}"
                       + (f", started again after {rec['interruptions']} interruption(s)" if rec["interruptions"] else ""))

    def fail(self, tid, why):
        rec = self.t[tid]
        d = rec["def"]
        if rec["attempts"] >= d["max_attempts"]:
            self.set_state(tid, "failed", f"{why}; {rec['attempts']} attempt(s) made, none left")
        else:
            rec["not_before"] = time.time() + d["retry_delay_s"]
            self.set_state(tid, "retry", f"{why}; attempt {rec['attempts'] + 1} of {d['max_attempts']} "
                                         f"in {d['retry_delay_s']:g} s")

    def finish(self, tid, rc, stopped_by=None):
        """The command of a task has ended (rc: exit code, negative for a signal)."""
        rec = self.t[tid]
        d = rec["def"]
        self.procs.pop(tid, None)
        st = stamp()
        wall = time.time() - rec.get("started_epoch", time.time())
        rec.update(ended=st, exit_code=rc, wall_s_last=round(wall, 1),
                   wall_s_total=round(rec.get("wall_s_total", 0) + wall, 1))
        rec.pop("pid", None)
        left = token_pids(rec["token"]) if rec.get("token") else []
        if left:                               # children that outlived the command: none may keep the GPU
            self.event(f"task {tid}: {len(left)} process(es) outlived the command; ending them")
            self.end_processes(rec["token"], 10)
        try:
            with open(rec["log"], "ab") as logf:
                logf.write(f"===== {st['utc']} | {st['ist']} | task {tid}: exit code {rc}, {wall:.0f} s"
                           f"{' | ' + stopped_by if stopped_by else ''} =====\n".encode())
        except OSError:
            pass
        if rec.get("removed_from_list"):
            self.set_state(tid, "cancelled", f"ended with exit code {rc} after it was removed from the task files")
            return
        if stopped_by in ("STOP", "signal"):
            rec["interruptions"] += 1
            self.set_state(tid, "interrupted", f"stopped by the runner ({stopped_by}), exit code {rc}")
            return
        if stopped_by == "cancel":
            self.set_state(tid, "cancelled", f"ended on request, exit code {rc}")
            return
        rec["attempts"] = rec.get("attempt_running", rec["attempts"] + 1)
        if stopped_by == "timeout":
            self.fail(tid, f"no end within timeout_s = {d['timeout_s']:g} s")
            return
        if rc != 0:
            self.fail(tid, f"exit code {rc}" if rc > 0 else f"ended by signal {-rc}")
            return
        if d["done_marker"] == "auto":
            try:
                with open(self.marker_path(d, {}), "w") as f:
                    f.write(f"{st['utc']} task {tid} ended with exit code 0\n")
            except OSError as e:
                self.fail(tid, f"exit code 0, but the marker could not be written: {e}")
                return
        ok, detail = self.check_success(d)
        rec["last_test"] = {"ok": ok, "detail": detail, **stamp()}
        if ok:
            self.set_state(tid, "done", f"exit code 0 and the success test holds; {wall:.0f} s")
        else:
            self.fail(tid, f"exit code 0, but the success test does not hold: {detail}")

    def reap(self):
        for tid, proc in list(self.procs.items()):
            rc = proc.poll()
            rec = self.t[tid]
            limit = rec["def"]["timeout_s"]
            if rc is None and limit and time.time() - rec["started_epoch"] > limit:
                how = self.end_processes(rec["token"], rec["def"]["stop_grace_s"])
                self.event(f"task {tid}: over its time limit, {how}")
                self.finish(tid, proc.wait(), "timeout")
            elif rc is not None:
                self.finish(tid, rc)

    def control(self):
        """retry.<id> and cancel.<id> files in <root>/queue/control."""
        try:
            names = sorted(os.listdir(self.p.control))
        except OSError:
            return
        for name in names:
            verb, _, tid = name.partition(".")
            path = os.path.join(self.p.control, name)
            if verb not in ("retry", "cancel") or not tid:
                continue
            try:
                os.remove(path)
            except OSError:
                continue
            rec = self.t.get(tid)
            if rec is None:
                self.event(f"control file {name}: no task {tid}")
            elif verb == "retry":
                if rec["state"] in ("failed", "cancelled", "retry"):
                    rec.update(attempts=0, interruptions=0)
                    rec.pop("not_before", None)
                    self.prechecked.discard(tid)
                    self.set_state(tid, "pending", "retry asked for; attempts counted from zero")
                else:
                    self.event(f"control file {name}: task {tid} is {rec['state']}, nothing to retry")
            elif verb == "cancel":
                if tid in self.procs:
                    how = self.end_processes(rec["token"], rec["def"]["stop_grace_s"])
                    self.event(f"task {tid}: cancel asked for, {how}")
                    self.finish(tid, self.procs[tid].wait(), "cancel")
                elif rec["state"] in ("pending", "retry"):
                    self.set_state(tid, "cancelled", "cancelled before it was started")
                else:
                    self.event(f"control file {name}: task {tid} is {rec['state']}, nothing to cancel")

    def schedule(self):
        """Start what can be started, in the order of the list. Returns True if something changed."""
        changed = False
        for tid in list(self.order):
            rec = self.t[tid]
            d = rec["def"]
            if rec["state"] == "retry" and time.time() >= rec.get("not_before", 0):
                self.set_state(tid, "pending", "next attempt")
            if rec["state"] != "pending":
                continue
            states = {dep: self.t[dep]["state"] if dep in self.t else "not in the list" for dep in d["depends_on"]}
            blocked = [f"{dep} ({s})" for dep, s in states.items() if s in ("failed", "cancelled", "not in the list")]
            waiting = [dep for dep, s in states.items() if s != "done"]
            if rec.get("blocked_by", []) != blocked or rec.get("waiting_for", []) != waiting:
                rec["blocked_by"], rec["waiting_for"], self.dirty = blocked, waiting, True
            if waiting:
                continue
            if tid not in self.prechecked:
                self.prechecked.add(tid)
                ok, detail = self.check_success(d)
                rec["last_test"] = {"ok": ok, "detail": detail, **stamp()}
                if ok:
                    self.set_state(tid, "done", "its success test holds already; not run")
                    changed = True
                    continue
            cpu, gpu = self.in_use()
            if cpu + d["cpu"] > self.cpus or gpu + d["gpu"] > self.gpus:
                continue
            self.start(tid)
            changed = True
        return changed

    def shutdown(self, why):
        """SIGTERM to every running task; SIGKILL to each one that is still there after its own stop_grace_s."""
        self.qstate = "stopping"
        tag = "STOP" if why.startswith("STOP") else "signal"
        self.event(f"stopping ({why}); {len(self.procs)} task(s) running")
        for tid in self.procs:
            send(token_pids(self.t[tid]["token"]), signal.SIGTERM)
        t0, killed = time.time(), {}
        while self.procs:
            for tid, proc in list(self.procs.items()):
                rec = self.t[tid]
                alive = token_pids(rec["token"])
                waited = time.time() - t0
                if proc.poll() is not None and not alive:
                    self.finish(tid, proc.returncode, tag)
                elif tid not in killed and waited > rec["def"]["stop_grace_s"]:
                    send(alive, signal.SIGKILL)
                    killed[tid] = time.time()
                    self.event(f"task {tid}: still there {waited:.0f} s after SIGTERM; SIGKILL to {len(alive)} process(es)")
                elif tid in killed and time.time() - killed[tid] > 20:
                    self.event(f"task {tid}: {len(alive)} process(es) did not go away after SIGKILL")
                    self.procs.pop(tid)
                    rec["interruptions"] += 1
                    self.set_state(tid, "interrupted", f"stopped by the runner ({tag}); processes left behind")
            time.sleep(0.2)
        self.qstate = "stopped"
        self.write_state(heartbeat=True)
        self.event(f"stopped ({why})")

    def run(self):
        self.lock_fd = os.open(self.p.lock, os.O_CREAT | os.O_RDWR, 0o644)
        try:
            fcntl.flock(self.lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print(f"another runner works on {self.p.root} (lock {self.p.lock}); not started", file=sys.stderr)
            return 3
        os.ftruncate(self.lock_fd, 0)
        os.write(self.lock_fd, f"{os.getpid()}\n".encode())
        signal.signal(signal.SIGTERM, self.on_signal)
        signal.signal(signal.SIGINT, self.on_signal)
        self.boot_enabled = service_enabled()
        self.event(f"runner started: pid {os.getpid()}, host {socket.gethostname()}, root {self.p.root}, "
                   f"{self.cpus} CPUs, {self.gpus} GPU(s), tasks from {self.p.tasks}")
        if os.path.exists(self.p.stop):
            self.qstate = "stopped"
            self.restore_quietly()
            self.write_state(heartbeat=True)
            self.event(f"the file {self.p.stop} exists: nothing is started. Remove it and start the runner again.")
            return 0
        self.restore()
        self.qstate = "running"
        errors = 0
        while True:
            try:
                code = self.one_pass()
                errors = 0
            except Exception:                                       # noqa: BLE001
                errors += 1
                self.event(f"ERROR in the runner's loop ({errors} in a row): " + traceback.format_exc().strip().replace("\n", " | ")[-600:])
                if errors >= 20:
                    self.event("20 errors in a row: the runner ends with exit code 1; running tasks are not signalled by it")
                    return 1
                code = None
            if code is not None:
                return code
            time.sleep(self.a.poll)

    def one_pass(self):
        """One look at everything. Returns the exit code when the runner has to end, else None."""
        if self.stop_signal is not None:
            self.shutdown(f"signal {self.stop_signal}")
            return 0
        if os.path.exists(self.p.stop):
            self.shutdown(f"STOP file {self.p.stop}")
            return 0
        self.reap()
        self.control()
        if time.time() - self.last_reload >= self.a.reload:
            self.reload()
        draining = os.path.exists(self.p.drain)
        want = "draining" if draining else "running"
        if want != self.qstate:
            self.qstate, self.dirty = want, True
            self.event("DRAIN file present: nothing new is started" if draining else "DRAIN file gone: starting tasks again")
        if not draining:
            while self.schedule():
                pass
        idle = not self.procs and not any(self.t[t]["state"] == "retry" for t in self.order)
        if idle and self.idle_since is None:
            self.idle_since, self.dirty = stamp(), True
            self.event("idle: nothing runs and nothing can be started ("
                       + ", ".join(f"{n} {s}" for s, n in sorted(self.counts().items())) + ")")
        elif not idle and self.idle_since is not None:
            self.idle_since, self.dirty = None, True
        beat = time.time() - self.last_heartbeat >= self.a.heartbeat
        if beat:
            self.gpu_last = gpu_snapshot()
        if self.dirty or beat:
            self.write_state(heartbeat=True)
            if beat:
                self.last_heartbeat = time.time()
        if idle and self.a.exit_when_idle and not draining:
            self.qstate = "stopped"
            self.write_state(heartbeat=True)
            self.event("idle and --exit-when-idle: runner ends")
            return 0
        return None

    def restore_quietly(self):
        """With a STOP file at the start: show the last known states, change nothing."""
        old = read_json(self.p.state)
        if old and old.get("schema") == SCHEMA:
            self.order = [t for t in old.get("order", []) if t in old.get("tasks", {})]
            self.t = {t: old["tasks"][t] for t in self.order}
            self.removed = old.get("removed", {})

    def on_signal(self, signum, _frame):
        self.stop_signal = signum


def uptime_s():
    try:
        with open("/proc/uptime") as f:
            return round(float(f.read().split()[0]))
    except (OSError, ValueError):
        return None


def gpu_snapshot():
    """One line of nvidia-smi for the heartbeat; None where there is no nvidia-smi or it does not answer."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return None
    try:
        r = subprocess.run([exe, "--query-gpu=name,utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw",
                            "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=15)
        if r.returncode != 0 or not r.stdout.strip():
            return {"error": (r.stderr or r.stdout).strip()[:200]}
        v = [x.strip() for x in r.stdout.strip().splitlines()[0].split(",")]
        return dict(zip(("name", "utilization_percent", "memory_used_mib", "memory_total_mib", "temperature_c",
                         "power_w"), v))
    except (OSError, subprocess.TimeoutExpired) as e:
        return {"error": str(e)[:200]}


def service_enabled():
    """True/False: is cqc-md-queue.service enabled for the next boot; None where systemctl cannot tell."""
    exe = shutil.which("systemctl")
    if not exe:
        return None
    try:
        r = subprocess.run([exe, "is-enabled", "cqc-md-queue.service"], capture_output=True, text=True, timeout=15)
        return r.stdout.strip() == "enabled"
    except (OSError, subprocess.TimeoutExpired):
        return None


# ------------------------------------------------------------------------------------------------ other commands
def cmd_status(a):
    p = Paths(a.root)
    st = read_json(p.state)
    if not st:
        print(f"no state file {p.state}")
        return 1
    if a.json:
        print(json.dumps(st, indent=1))
        return 0
    q = st["queue"]
    age = time.time() - q.get("heartbeat_epoch", 0)
    print(f"runner: {q['state']}, pid {q['pid']} on {q['host']}, started {q['started_ist']} by {q.get('started_by')}; "
          f"last heartbeat {q['heartbeat_ist']} ({age:.0f} s ago); starts at boot: {q.get('starts_at_boot')}")
    print(f"in use: {q['cpus_in_use']} of {q['cpus_total']} CPUs, {q['gpus_in_use']} of {q['gpus_total']} GPU; "
          f"tasks: " + ", ".join(f"{n} {s}" for s, n in sorted(st["summary"].items()))
          + ("; STOP file present" if q.get("stop_file_present") else "")
          + ("; DRAIN file present" if q.get("drain_file_present") else ""))
    print(f"{'task':28s} {'state':12s} {'resources':10s} {'tries':5s} {'intr':4s} since (IST)                  note")
    for tid in st["order"]:
        r = st["tasks"][tid]
        print(f"{tid:28s} {r['state']:12s} {r['def']['resources']:10s} {r['attempts']:<5d} {r['interruptions']:<4d} "
              f"{r.get('state_since', {}).get('ist', ''):28s} {r.get('note', '')}"
              + (f" [waits for {', '.join(r['waiting_for'])}]" if r["state"] == "pending" and r.get("waiting_for") else ""))
    for e in q.get("task_file_errors", []):
        print(f"TASK FILE PROBLEM: {e}")
    for w in q.get("warnings", []):
        print(f"note: {w}")
    return 0


def collect(paths_new, folder, replace=()):
    """Definitions of the folder (without the files named in `replace`) followed by those of the new files."""
    per_file, errors = [], []
    for path in list_task_files(folder):
        if os.path.basename(path) in replace:
            continue
        defs, errs = read_task_file(path)
        per_file.append((os.path.basename(path), defs or []))
        errors += [f"(already there) {e}" for e in errs]
    n_old = len(per_file)
    for path in paths_new:
        defs, errs = read_task_file(path)
        per_file.append((os.path.basename(path), defs or []))
        errors += errs
    defs, dup = merge_files(per_file)
    return defs, errors + dup, n_old


def cmd_check(a, quiet=False):
    p = Paths(a.root)
    cpus = a.cpus or os.cpu_count() or 1
    defs, errors, _ = collect(a.files, p.tasks, replace={os.path.basename(f) for f in a.files})
    bad, warn = structural_problems(defs, cpus, a.gpus)
    new = {os.path.basename(f) for f in a.files}
    hard = [e for e in errors if not e.startswith("(already there)")] + [f"task {t}: {w}" for t, w in bad.items()]
    for e in errors:
        print(("ERROR: " if not e.startswith("(already there)") else "note: ") + e)
    for t, w in bad.items():
        print(f"ERROR: task {t}: {w}")
    for w in warn:
        print(f"note: {w}")
    if not quiet or hard:
        for d in defs:
            if d["source"] in new:
                print(f"  {d['id']:28s} {d['resources']:10s} after {', '.join(d['depends_on']) or '-':30s} "
                      f"done when: {d['done_marker'] or ''}{' and ' if d['done_marker'] and d['success_test'] else ''}"
                      f"{'success_test' if d['success_test'] else ''}")
    print(f"{'NOT accepted' if hard else 'accepted'}: {len([d for d in defs if d['source'] in new])} task(s) in "
          f"{len(a.files)} file(s), checked for {cpus} CPUs and {a.gpus} GPU(s)"
          + (f", together with the tasks in {p.tasks}" if os.path.isdir(p.tasks) else ""))
    return 1 if hard else 0


def cmd_add(a):
    p = Paths(a.root)
    p.make()
    for f in a.files:
        target = os.path.join(p.tasks, os.path.basename(f))
        if os.path.exists(target) and not a.replace:
            print(f"ERROR: {target} exists; give the file another name, or --replace")
            return 1
        if not f.endswith(".json"):
            print(f"ERROR: {f}: the name must end in .json")
            return 1
    if cmd_check(a, quiet=True) != 0:
        return 1
    for f in a.files:
        target = os.path.join(p.tasks, os.path.basename(f))
        tmp = os.path.join(p.tasks, f".{os.path.basename(f)}.tmp{os.getpid()}")
        shutil.copyfile(f, tmp)
        os.replace(tmp, target)
        print(f"added: {target}")
    return 0


def cmd_control(a):
    p = Paths(a.root)
    p.make()
    if a.cmd in ("retry", "cancel"):
        for tid in a.ids:
            if not ID_RE.match(tid):
                print(f"bad id {tid!r}")
                return 2
            open(os.path.join(p.control, f"{a.cmd}.{tid}"), "w").close()
            print(f"{a.cmd} asked for task {tid} (the runner reads the request within seconds; see `queue.py status`)")
    elif a.cmd == "stop":
        with open(p.stop, "w") as f:
            f.write(f"{stamp()['utc']} written by queue.py stop\n")
        print(f"written: {p.stop} (running tasks get SIGTERM; the runner ends and starts nothing while the file exists)")
    elif a.cmd == "drain":
        with open(p.drain, "w") as f:
            f.write(f"{stamp()['utc']} written by queue.py drain\n")
        print(f"written: {p.drain} (nothing new is started; running tasks finish)")
    elif a.cmd == "resume":
        for f in (p.stop, p.drain):
            if os.path.exists(f):
                os.remove(f)
                print(f"removed: {f}")
        print("if the runner had ended, start it: systemctl start cqc-md-queue.service")
    return 0


def cmd_json(a):
    """Conditions on a JSON file, for success tests."""
    try:
        with open(a.file) as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        print(f"{a.file}: cannot be read as JSON ({e})")
        return 2
    code = 0
    for cond in a.conditions:
        m = re.fullmatch(r"([^=!]+)(!=|=)(.*)", cond, re.S)
        if not m:
            print(f"bad condition {cond!r}: expected KEY=VALUE or KEY!=VALUE")
            return 2
        key, op, text = m.groups()
        node = data
        try:
            for part in key.split("."):
                node = node[int(part)] if isinstance(node, list) else node[part]
        except (KeyError, IndexError, ValueError, TypeError):
            print(f"{a.file}: no key {key}")
            return 2
        try:
            want = json.loads(text)
        except ValueError:
            want = text
        if isinstance(node, bool) or isinstance(want, bool):
            same = node is want
        elif isinstance(node, (int, float)) and isinstance(want, (int, float)):
            same = abs(node - want) <= 1e-9 * max(1.0, abs(node), abs(want))
        else:
            same = node == want
        if same != (op == "="):
            shown = json.dumps(node)
            print(f"{a.file}: {key} is {shown if len(shown) < 200 else shown[:200] + '...'}, asked for: {cond}")
            code = 1
    return code


# ------------------------------------------------------------------------------------------------ self-test
def selftest(a):
    """Dummy tasks (sleep, echo) against a runner started as a separate process. Prints one line per check."""
    base = a.dir or tempfile.mkdtemp(prefix="cqc_queue_selftest_")
    os.makedirs(base, exist_ok=True)
    results, me = [], os.path.abspath(__file__)

    def check(name, ok, detail=""):
        if callable(ok):                       # a missing trace line must give FAIL, not a crash of the test
            try:
                ok = ok()
            except Exception as e:             # noqa: BLE001
                ok, detail = False, f"{detail} ({type(e).__name__}: {e})"
        results.append((name, bool(ok), detail))
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""), flush=True)

    def fresh(name):
        root = os.path.join(base, name)
        shutil.rmtree(root, ignore_errors=True)
        p = Paths(root)
        p.make()
        return root, p, os.path.join(root, "trace.txt")

    def dummy(tid, seconds, res, deps=(), extra="", **kw):
        # every dummy writes "start <id> <pid> <time>" and "end <id> <pid> <time>" to the trace
        d = {"id": tid, "cwd": "$ROOTDIR/work", "resources": res, "depends_on": list(deps), "done_marker": f"$ROOTDIR/marks/{tid}",
             "command": [f'echo "start {tid} $$ $(date +%s.%N)" >> "$TRACE"', extra or ":", f"sleep {seconds}",
                         f'echo "end {tid} $$ $(date +%s.%N)" >> "$TRACE"', f'mkdir -p "$ROOTDIR/marks"; touch "$ROOTDIR/marks/{tid}"'],
             "retry_delay_s": 0.5, "stop_grace_s": 10}
        d.update(kw)
        return d

    def write_tasks(p, name, tasks, root, trace):
        tmp = os.path.join(p.queue, f".{name}.tmp")
        with open(tmp, "w") as f:
            json.dump({"env": {"ROOTDIR": root, "TRACE": trace}, "tasks": tasks}, f, indent=1)
        os.replace(tmp, os.path.join(p.tasks, name))

    def launch(root, cpus=4, extra=()):
        out = open(os.path.join(root, "runner.out"), "ab")
        return subprocess.Popen([sys.executable, me, "run", "--root", root, "--cpus", str(cpus), "--poll", "0.1",
                                 "--reload", "0.3", "--heartbeat", "1", *extra], stdout=out, stderr=subprocess.STDOUT)

    def trace_of(trace):
        ev = []
        if os.path.exists(trace):
            for line in open(trace):
                f = line.split()
                if len(f) == 4:
                    ev.append((f[0], f[1], int(f[2]), float(f[3])))
        return ev

    def wait_for(cond, seconds):
        end = time.time() + seconds
        while time.time() < end:
            if cond():
                return True
            time.sleep(0.1)
        return cond()

    def state_of(p, tid=None):
        st = read_json(p.state, {}) or {}
        if tid is None:
            return st
        return st.get("tasks", {}).get(tid, {}).get("state")

    def when(ev, kind, tid, n=0):
        hits = [t for k, i, _pid, t in ev if k == kind and i == tid]
        return hits[n] if len(hits) > n else None

    # ---- 1 and 2: dependencies; one gpu task at a time; CPU tasks beside it within the CPUs
    root, p, trace = fresh("t1_deps_resources")
    write_tasks(p, "10.json", [dummy("A", 1, "cpu:1"), dummy("B", 1, "cpu:1", ["A"]), dummy("C", 0.3, "cpu:1", ["B"]),
                               dummy("G1", 1.5, "gpu"), dummy("G2", 1.5, "gpu"), dummy("G3", 1.5, "gpu", ["A"]),
                               dummy("W1", 2, "cpu:2"), dummy("W2", 2, "cpu:2")], root, trace)
    r = launch(root, cpus=4, extra=["--exit-when-idle"])
    rc = r.wait(timeout=120)
    ev = trace_of(trace)
    st = state_of(p)
    ids = ["A", "B", "C", "G1", "G2", "G3", "W1", "W2"]
    check("1a all eight dummy tasks done, runner ended with 0", rc == 0 and all(st["tasks"][i]["state"] == "done" for i in ids),
          f"exit {rc}; " + " ".join(f"{i}={st['tasks'][i]['state']}" for i in ids))
    check("1b dependencies: B starts after A ended, C after B ended, G3 after A ended",
          lambda: when(ev, "end", "A") <= when(ev, "start", "B") and when(ev, "end", "B") <= when(ev, "start", "C")
          and when(ev, "end", "A") <= when(ev, "start", "G3"))
    spans = {i: (when(ev, "start", i) or 0.0, when(ev, "end", i) or 0.0) for i in ids}
    res = {"A": (0, 1), "B": (0, 1), "C": (0, 1), "G1": (1, 1), "G2": (1, 1), "G3": (1, 1), "W1": (0, 2), "W2": (0, 2)}
    points = sorted(t for s in spans.values() for t in s)
    max_gpu = max(sum(res[i][0] for i in ids if spans[i][0] <= t < spans[i][1]) for t in points)
    max_cpu = max(sum(res[i][1] for i in ids if spans[i][0] <= t < spans[i][1]) for t in points)
    beside = any(spans[g][0] < spans[c][1] and spans[c][0] < spans[g][1] for g in ("G1", "G2", "G3") for c in ("A", "B", "W1", "W2"))
    check("2a never more than one gpu task at a time (the three gpu tasks did not overlap)", max_gpu == 1 and len(ev) == 16,
          f"largest number of gpu tasks at one time: {max_gpu}; {len(ev)} trace lines of 16")
    check("2b CPU tasks ran beside a gpu task, and never more than the 4 CPUs given", beside and 3 <= max_cpu <= 4,
          f"largest CPU sum {max_cpu} of 4; a CPU task overlapped a gpu task: {beside}")
    hist = st["tasks"]["B"]["history"]
    check("2c every change of state is recorded with UTC and IST time",
          [h["to"] for h in hist] == ["pending", "running", "done"] and all(h["utc"].endswith("Z") and h["ist"].endswith("IST") for h in hist),
          " -> ".join(h["to"] for h in hist))
    logtxt = open(os.path.join(p.logs, "A.log")).read()
    check("2d each task has its own log with a start and an end line", "attempt 1 of 3" in logtxt and "exit code 0" in logtxt)

    # ---- 3: the runner is killed (SIGKILL) while a task runs, and started again
    root, p, trace = fresh("t3_kill_restart")
    write_tasks(p, "10.json", [dummy("L", 4, "gpu"), dummy("D", 0.3, "cpu:1", ["L"]), dummy("G", 0.3, "gpu")], root, trace)
    r = launch(root)
    ok_started = wait_for(lambda: when(trace_of(trace), "start", "L") is not None and state_of(p, "L") == "running", 30)
    first_pid = [pid for k, i, pid, _t in trace_of(trace) if k == "start" and i == "L"][0] if ok_started else None
    os.kill(r.pid, signal.SIGKILL)
    r.wait()
    alive_after_kill = first_pid is not None and os.path.exists(f"/proc/{first_pid}")
    check("3a setting: runner killed with SIGKILL in the middle of task L; the state file still says running; L's process is an orphan",
          ok_started and state_of(p, "L") == "running" and alive_after_kill, f"L pid {first_pid} alive: {alive_after_kill}")
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=120)
    ev = trace_of(trace)
    st = state_of(p)
    starts_l = [(pid, t) for k, i, pid, t in ev if k == "start" and i == "L"]
    ends_l = [(pid, t) for k, i, pid, t in ev if k == "end" and i == "L"]
    hist = [h["to"] for h in st["tasks"]["L"]["history"]]
    check("3b the orphan was ended before L was started again (it never wrote its end line; its process is gone)",
          lambda: len(starts_l) == 2 and len(ends_l) == 1 and ends_l[0][0] == starts_l[1][0] and not os.path.exists(f"/proc/{first_pid}"),
          f"starts {len(starts_l)}, ends {len(ends_l)}")
    check("3c L was run again and is done; history pending, running, interrupted, pending, running, done; 1 interruption, 1 attempt",
          rc == 0 and hist == ["pending", "running", "interrupted", "pending", "running", "done"]
          and st["tasks"]["L"]["interruptions"] == 1 and st["tasks"]["L"]["attempts"] == 1, " -> ".join(hist))
    check("3d D (after L) and G (gpu) ran only after the second L had ended; never two gpu tasks at once",
          lambda: st["tasks"]["D"]["state"] == "done" and st["tasks"]["G"]["state"] == "done"
          and when(ev, "start", "D") >= ends_l[0][1] and when(ev, "start", "G") >= ends_l[0][1])

    # ---- 4: a task added while the runner works; a broken file; a duplicate id; queue.py add
    root, p, trace = fresh("t4_add")
    write_tasks(p, "10.json", [dummy("L", 8, "gpu")], root, trace)
    r = launch(root)
    wait_for(lambda: state_of(p, "L") == "running", 30)
    newfile = os.path.join(root, "20_added.json")
    with open(newfile, "w") as f:
        json.dump({"env": {"ROOTDIR": root, "TRACE": trace}, "tasks": [dummy("N", 0.3, "cpu:1"), dummy("M", 0.3, "gpu", ["L"])]}, f)
    rc_add = subprocess.run([sys.executable, me, "add", newfile, "--root", root, "--cpus", "4"], capture_output=True, text=True)
    with open(os.path.join(p.tasks, "30_broken.json"), "w") as f:
        f.write('{"tasks": [ {"id": "X", ')
    dupfile = os.path.join(root, "40_dup.json")
    with open(dupfile, "w") as f:
        json.dump({"env": {"ROOTDIR": root, "TRACE": trace}, "tasks": [dummy("N", 0.3, "cpu:1")]}, f)
    rc_dup = subprocess.run([sys.executable, me, "add", dupfile, "--root", root, "--cpus", "4"], capture_output=True, text=True)
    big = os.path.join(root, "50_big.json")
    with open(big, "w") as f:
        json.dump([dummy("BIG", 0.1, "cpu:64")], f)
    rc_big = subprocess.run([sys.executable, me, "check", big, "--root", root, "--cpus", "4"], capture_output=True, text=True)
    wait_for(lambda: state_of(p, "M") == "done", 60)
    st = state_of(p)
    ev = trace_of(trace)
    check("4a a file added with `queue.py add` while L ran: N was run beside L, M (after L) when L had ended",
          lambda: rc_add.returncode == 0 and st["tasks"]["N"]["state"] == "done" and st["tasks"]["M"]["state"] == "done"
          and when(ev, "end", "N") < when(ev, "end", "L") and when(ev, "start", "M") >= when(ev, "end", "L"),
          rc_add.stdout.strip().splitlines()[-1] if rc_add.stdout else rc_add.stderr.strip()[-200:])
    check("4b a file that is not JSON is reported in queue.json and the runner works on",
          any("30_broken.json" in e for e in st["queue"]["task_file_errors"]) and r.poll() is None)
    check("4c `queue.py add` refuses a second task with an id that exists, and `check` a task that needs more CPUs than there are",
          rc_dup.returncode == 1 and "already defined" in rc_dup.stdout and rc_big.returncode == 1 and "needs 64 CPUs" in rc_big.stdout)
    hb1 = read_json(p.heartbeat)
    time.sleep(2.2)
    hb2 = read_json(p.heartbeat)
    check("4d the heartbeat file is written again at its interval (1 s in this test) while nothing happens",
          hb1 and hb2 and hb2["epoch"] > hb1["epoch"] and time.time() - hb2["epoch"] < 3 and hb2["idle"] is True
          and hb2["counts"].get("done") == 3, f"{hb1 and hb1['epoch']} -> {hb2 and hb2['epoch']}")
    r.send_signal(signal.SIGTERM)
    check("4e SIGTERM ends an idle runner with exit code 0", r.wait(timeout=30) == 0)

    # ---- 5: STOP file
    root, p, trace = fresh("t5_stop")
    trap = "trap 'echo \"term S $$ $(date +%s.%N)\" >> \"$TRACE\"; exit 143' TERM"
    write_tasks(p, "10.json", [dummy("S", "5 & wait $!", "gpu", extra=trap), dummy("K", "30", "cpu:1", extra="trap '' TERM", stop_grace_s=1.5),
                               dummy("T", 0.2, "cpu:1", ["S"])], root, trace)
    r = launch(root)
    wait_for(lambda: state_of(p, "S") == "running" and state_of(p, "K") == "running"
             and when(trace_of(trace), "start", "S") is not None and when(trace_of(trace), "start", "K") is not None, 30)
    time.sleep(0.5)
    t0 = time.time()
    subprocess.run([sys.executable, me, "stop", "--root", root], capture_output=True)
    rc = r.wait(timeout=60)
    took = time.time() - t0
    ev = trace_of(trace)
    st = state_of(p)
    k_pid = ([pid for k, i, pid, _t in ev if k == "start" and i == "K"] or [0])[0]
    check("5a STOP file: the runner ended with exit code 0; S got SIGTERM and ended itself; K, which ignores SIGTERM, got SIGKILL after its grace time",
          lambda: rc == 0 and when(ev, "term", "S") is not None and when(ev, "end", "S") is None and k_pid and not os.path.exists(f"/proc/{k_pid}")
          and st["tasks"]["S"]["state"] == "interrupted" and st["tasks"]["K"]["state"] == "interrupted" and st["queue"]["state"] == "stopped"
          and st["tasks"]["T"]["state"] == "pending" and took < 9,
          f"exit {rc} after {took:.1f} s; S {st['tasks']['S']['state']}, K {st['tasks']['K']['state']}, T {st['tasks']['T']['state']}")
    n_before = len(trace_of(trace))
    r = launch(root)
    rc = r.wait(timeout=30)
    check("5b started again with the STOP file in place: ends at once with 0 and starts nothing",
          rc == 0 and len(trace_of(trace)) == n_before and state_of(p)["queue"]["state"] == "stopped")
    subprocess.run([sys.executable, me, "resume", "--root", root], capture_output=True)
    subprocess.run([sys.executable, me, "cancel", "K", "--root", root], capture_output=True)
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=60)
    st = state_of(p)
    check("5c STOP removed, K cancelled: S is run again and done, T after it; K is cancelled and was not started again",
          rc == 0 and st["tasks"]["S"]["state"] == "done" and st["tasks"]["T"]["state"] == "done" and st["tasks"]["K"]["state"] == "cancelled"
          and len([1 for k, i, _p, _t in trace_of(trace) if k == "start" and i == "K"]) == 1,
          f"S {st['tasks']['S']['state']}, T {st['tasks']['T']['state']}, K {st['tasks']['K']['state']}")

    # ---- 6: the success test decides
    root, p, trace = fresh("t6_success_test")
    os.makedirs(os.path.join(root, "marks"))
    open(os.path.join(root, "marks", "PRE"), "w").close()
    count = os.path.join(root, "count_FLAKY")
    tasks = [dummy("PRE", 0.2, "cpu:1"),                                        # its marker exists already
             {"id": "LIE", "cwd": "$ROOTDIR/work", "resources": "cpu:1", "command": 'echo "start LIE $$ $(date +%s.%N)" >> "$TRACE"',
              "success_test": "test -f $ROOTDIR/never_written", "max_attempts": 2, "retry_delay_s": 0.3},
             dummy("AFTER_LIE", 0.2, "cpu:1", ["LIE"]),
             {"id": "FLAKY", "cwd": "$ROOTDIR/work", "resources": "gpu", "done_marker": "auto", "retry_delay_s": 0.3,
              "command": [f'echo "start FLAKY $$ $(date +%s.%N)" >> "$TRACE"', f'echo x >> {count}', f'test $(wc -l < {count}) -ge 2']},
             {"id": "ENVT", "cwd": "$ROOTDIR/work", "resources": "cpu:3", "done_marker": "$ROOTDIR/marks/ENVT",
              "command": 'echo "$OMP_NUM_THREADS $OPENBLAS_NUM_THREADS $MKL_NUM_THREADS $QUEUE_TASK_ID $QUEUE_CPUS $PWD" > "$ROOTDIR/marks/ENVT"'},
             dummy("SLOW", 30, "cpu:1", timeout_s=1, max_attempts=1, stop_grace_s=2)]
    write_tasks(p, "10.json", tasks, root, trace)
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=120)
    st = state_of(p)
    ev = trace_of(trace)
    T = st["tasks"]
    check("6a a task whose success test holds already is done without being run",
          T["PRE"]["state"] == "done" and when(ev, "start", "PRE") is None, T["PRE"]["note"])
    check("6b exit code 0 without the success test holding is a failure: two attempts, then failed; the task after it is not started",
          T["LIE"]["state"] == "failed" and T["LIE"]["attempts"] == 2 and len([1 for k, i, _p, _t in ev if i == "LIE"]) == 2
          and T["AFTER_LIE"]["state"] == "pending" and when(ev, "start", "AFTER_LIE") is None and "LIE (failed)" in " ".join(T["AFTER_LIE"]["blocked_by"]),
          T["LIE"]["note"])
    check("6c a command that fails once is tried again and then done (marker written by the runner)",
          T["FLAKY"]["state"] == "done" and T["FLAKY"]["attempts"] == 2 and os.path.exists(os.path.join(p.done, "FLAKY.done")))
    envline = open(os.path.join(root, "marks", "ENVT")).read().split()
    check("6d a cpu:3 task gets the thread variables set to 3 and runs in its cwd",
          envline == ["3", "3", "3", "ENVT", "3", os.path.join(root, "work")], " ".join(envline))
    check("6e a task over its timeout_s is ended and counted as failed", T["SLOW"]["state"] == "failed" and "timeout_s" in T["SLOW"]["note"],
          T["SLOW"]["note"])
    hb = read_json(p.heartbeat)
    check("6f the runner reports idle with the failed tasks under `attention`; exit code 0",
          rc == 0 and hb["idle"] and not hb["all_done"] and any(x.startswith("LIE: failed") for x in hb["attention"])
          and any(x.startswith("AFTER_LIE: waits for") for x in hb["attention"]), "; ".join(hb["attention"])[:160])
    subprocess.run([sys.executable, me, "retry", "LIE", "--root", root], capture_output=True)
    open(os.path.join(root, "never_written"), "w").close()
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=60)
    st = state_of(p)
    check("6g `queue.py retry` after the cause is gone: LIE done (its test holds now), AFTER_LIE run and done",
          rc == 0 and st["tasks"]["LIE"]["state"] == "done" and st["tasks"]["AFTER_LIE"]["state"] == "done")

    # ---- 7: DRAIN, a second runner on the same root, a changed definition
    root, p, trace = fresh("t7_drain_lock")
    write_tasks(p, "10.json", [dummy("P", 2, "gpu"), dummy("Q", 0.2, "gpu")], root, trace)
    r = launch(root)
    wait_for(lambda: state_of(p, "P") == "running", 30)
    subprocess.run([sys.executable, me, "drain", "--root", root], capture_output=True)
    second = subprocess.run([sys.executable, me, "run", "--root", root, "--cpus", "4"], capture_output=True, text=True, timeout=30)
    wait_for(lambda: state_of(p, "P") == "done", 30)
    time.sleep(1.0)
    st = state_of(p)
    check("7a DRAIN: the running task finished, the next one was not started, the runner stays alive",
          st["tasks"]["P"]["state"] == "done" and st["tasks"]["Q"]["state"] == "pending" and r.poll() is None and st["queue"]["state"] == "draining")
    check("7b a second runner on the same root refuses to start (exit code 3)", second.returncode == 3, second.stderr.strip()[:100])
    write_tasks(p, "10.json", [dummy("P", 2, "gpu"), dummy("Q", 0.4, "gpu")], root, trace)     # Q changed
    wait_for(lambda: any("definition changed" in h["note"] for h in state_of(p)["tasks"]["Q"]["history"]), 20)
    subprocess.run([sys.executable, me, "resume", "--root", root], capture_output=True)
    wait_for(lambda: state_of(p, "Q") == "done", 30)
    st = state_of(p)
    check("7c DRAIN removed: Q is run, with its changed definition; P, unchanged and done, is not run again",
          st["tasks"]["Q"]["state"] == "done" and any("definition changed" in h["note"] for h in st["tasks"]["Q"]["history"])
          and len([1 for k, i, _p, _t in trace_of(trace) if k == "start" and i == "P"]) == 1)
    r.send_signal(signal.SIGTERM)
    r.wait(timeout=30)
    # the marker of P is removed while the runner is down: at the next start P is run again
    os.remove(os.path.join(root, "marks", "P"))
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=60)
    st = state_of(p)
    check("7d a task recorded as done whose success test no longer holds at the start of the runner is run again",
          rc == 0 and st["tasks"]["P"]["state"] == "done" and len([1 for k, i, _p, _t in trace_of(trace) if k == "start" and i == "P"]) == 2
          and any("does not hold any more" in h["note"] for h in st["tasks"]["P"]["history"]))

    # ---- 8: the json helper used in success tests
    jf = os.path.join(base, "json_test.json")
    with open(jf, "w") as f:
        json.dump({"status": "complete", "test_only": False, "production": {"complete": True, "completed_ns": 100.0},
                   "invocations": [{"exit_code": 143}, {"exit_code": 0}], "n": 0}, f)

    def jrun(*conds):
        return subprocess.run([sys.executable, me, "json", jf, *conds], capture_output=True, text=True).returncode

    check("8a `queue.py json`: 0 when all conditions hold, 1 when one does not, 2 for a missing key or file",
          jrun("status=complete", "production.complete=true", "production.completed_ns=100", "test_only=false",
               "invocations.-1.exit_code=0", "invocations.0.exit_code!=0") == 0
          and jrun("status=complete", "production.completed_ns=95") == 1 and jrun("test_only=true") == 1
          and jrun("test_only=0") == 1 and jrun("n=false") == 1 and jrun("status!=complete") == 1
          and jrun("production.absent=1") == 2
          and subprocess.run([sys.executable, me, "json", jf + ".none", "a=1"], capture_output=True).returncode == 2)

    # ---- 9: a state file of the last run that cannot be used
    root, p, trace = fresh("t9_bad_state")
    write_tasks(p, "10.json", [dummy("A", 0.2, "cpu:1")], root, trace)
    with open(p.state, "w") as f:
        json.dump({"schema": SCHEMA, "queue": {}, "order": ["A", "B"], "tasks": {"A": {"def": "not a definition", "state": "running"}, "B": 5}}, f)
    r = launch(root, extra=["--exit-when-idle"])
    rc = r.wait(timeout=60)
    st = state_of(p)
    check("9a a state file with unusable content does not keep the runner from starting: it says so and works from the task files",
          lambda: rc == 0 and st["tasks"]["A"]["state"] == "done" and "could not be taken over" in open(p.events).read(), f"exit {rc}")

    n_fail = sum(1 for _n, ok, _d in results if not ok)
    summary = {"when": stamp(), "host": socket.gethostname(), "python": sys.version.split()[0], "machine": os.uname().machine,
               "script_sha256": hashlib.sha256(open(me, "rb").read()).hexdigest(), "folder": base,
               "checks": [{"check": n, "passed": ok, "detail": d} for n, ok, d in results],
               "n_checks": len(results), "n_failed": n_fail}
    write_json_atomic(os.path.join(base, "selftest.json"), summary)
    print(f"{len(results) - n_fail} of {len(results)} checks passed; record: {os.path.join(base, 'selftest.json')}")
    if not a.keep and not a.dir and n_fail == 0:
        shutil.rmtree(base, ignore_errors=True)
    return 1 if n_fail else 0


# ------------------------------------------------------------------------------------------------ command line
def main():
    ap = argparse.ArgumentParser(description="Task runner of the rented GPU machine (md/gpu/GPU_RUN.md).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(sp, cpus=True):
        sp.add_argument("--root", default=ROOT_DEFAULT, help=f"folder with queue/, status/ and logs/ (default {ROOT_DEFAULT})")
        if cpus:
            sp.add_argument("--cpus", type=int, default=0, help="CPUs to give out (default: all of the machine)")
            sp.add_argument("--gpus", type=int, default=1, help="GPU tasks at one time (default 1)")

    sp = sub.add_parser("run", help="the runner")
    common(sp)
    sp.add_argument("--poll", type=float, default=2.0, help="seconds between two looks at the tasks (default 2)")
    sp.add_argument("--reload", type=float, default=10.0, help="seconds between two readings of the task files (default 10)")
    sp.add_argument("--heartbeat", type=float, default=60.0, help="seconds between two heartbeats (default 60)")
    sp.add_argument("--exit-when-idle", action="store_true", help="end when nothing runs and nothing can be started (tests)")
    sp = sub.add_parser("status", help="table of the tasks")
    common(sp, cpus=False)
    sp.add_argument("--json", action="store_true")
    for name in ("check", "add"):
        sp = sub.add_parser(name, help="validate task files" + (" and copy them into the task folder" if name == "add" else ""))
        common(sp)
        sp.add_argument("files", nargs="+")
        if name == "add":
            sp.add_argument("--replace", action="store_true", help="replace a file of the same name")
    for name in ("retry", "cancel"):
        sp = sub.add_parser(name)
        common(sp, cpus=False)
        sp.add_argument("ids", nargs="+")
    for name in ("stop", "drain", "resume"):
        common(sub.add_parser(name), cpus=False)
    sp = sub.add_parser("json", help="conditions on a JSON file (for success tests)")
    sp.add_argument("file")
    sp.add_argument("conditions", nargs="+", metavar="KEY=VALUE")
    sp = sub.add_parser("selftest", help="dummy tasks; about a minute")
    sp.add_argument("--dir", help="folder for the test (default: a new temporary folder, removed if all passes)")
    sp.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    if a.cmd == "run":
        return Runner(a).run()
    if a.cmd == "status":
        return cmd_status(a)
    if a.cmd == "check":
        return cmd_check(a)
    if a.cmd == "add":
        return cmd_add(a)
    if a.cmd == "json":
        return cmd_json(a)
    if a.cmd == "selftest":
        return selftest(a)
    return cmd_control(a)


if __name__ == "__main__":
    sys.exit(main())
