#!/usr/bin/env python3
"""deadman.py: the switch ON the rented GPU machine that ends the billing when the Pi has gone silent.

Run by systemd every 10 minutes, the first time 45 minutes after boot (cqc-md-deadman.timer starts
cqc-md-deadman.service; ops/vm/setup_gpu.sh installs and starts both; md/gpu/GPU_RUN.md, section 4):

    python3 md/gpu/deadman.py --yes     what the timer runs
    python3 md/gpu/deadman.py           DRY RUN: says what it sees and what it would do; changes nothing
    python3 md/gpu/deadman.py --state   the dry run in one word: armed | silent | off  (ops/vm/setup_gpu.sh asks this
                                        BEFORE it enables the timer: "silent" would shut the machine down at once)

WHY.  Everything that enforces the cost ceiling runs on the Pi.  If the Pi, its network, its cron or its IBM
session fails, the machine would go on costing US$0.84 per hour until the credit is gone.  The Pi's watchdog
(ops/pi/cqc_gpu_watchdog.py, every 5 minutes) therefore touches /root/cqc-sfa-md/status/PI_SEEN at every SSH
contact, and this script looks at the age of that file.

WHAT.  PI_SEEN older than 240 minutes (never touched: 240 minutes after boot, or after this switch first ran,
whichever is later; if the switch can neither read nor write its own file status/deadman.json, so that it cannot
know when it first ran: 240 minutes after boot) means that nothing on the Pi watches the machine any more.  Then:
  1. the task queue is stopped cleanly:  python3 md/gpu/queue.py stop  (every task gets SIGTERM; run_md.py ends
     within 500 steps and leaves a consistent folder);
  2. the end of the task runner is waited for, 10 minutes at most;
  3. the queue's STOP file is removed again if this run made it, so that the work goes on by itself when the
     machine is started again (a STOP file that a person put there stays);
  4. the machine is shut down:  shutdown -h now.  Its boot volume keeps everything.
At the next boot the timer gives the Pi 45 minutes to be seen; if it is not, the machine shuts down again.

OFF:  touch /root/cqc-sfa-md/status/NO_DEADMAN   (remove the file to arm the switch again).
LOG:  /root/cqc-sfa-md/status/deadman.log, a line whenever the state changes and for every action, UTC and IST;
      status/deadman.json is rewritten at every run (state, age of PI_SEEN).  The Pi pulls both with status/.

It shuts down only with --yes and only on an x86-64 machine (the rented machine is one; the Pi is aarch64, so a
test there can never switch the Pi off).  Exit codes: 0 nothing to do, or shut down; 1 the shutdown command
failed (the timer tries again in 10 minutes); 3 --yes on a machine that is not x86-64.
NOT VERIFIED on 6 Oct 2026 (no machine had this file yet): that IBM reports an instance that was shut down from
inside as "stopped" and stops billing it as running.  Test it once in the first hour (GPU_RUN.md, section 4).
Standard library only; reads and prints no secrets.
"""
import argparse, datetime, json, os, platform, subprocess, sys, time
from pathlib import Path

ROOT = Path(os.environ.get("CQC_DEADMAN_ROOT", "/root/cqc-sfa-md"))      # the variable is for tests
LIMIT_MIN, WAIT_MIN = 240.0, 10.0
QUEUE = [sys.executable, str(Path(__file__).resolve().with_name("queue.py"))]
SHUTDOWN = ["shutdown", "-h", "now"]
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30))
sleep = time.sleep                                           # replaced in the tests, as the four functions below


def now():
    return time.time()


def uptime_s():
    try:
        with open("/proc/uptime") as f:
            return float(f.read().split()[0])
    except (OSError, ValueError, IndexError):
        return 0.0


def machine_ok():
    return platform.machine() == "x86_64"


def run(cmd, timeout=120):
    """(exit code, last line of the output).  Never raises."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, ((r.stdout + r.stderr).strip().splitlines() or [""])[-1][:200]
    except (OSError, subprocess.TimeoutExpired) as e:
        return 999, e.__class__.__name__


def runner_pid():
    """PID of the living task runner (md/gpu/queue.py writes it into its lock file), or None.  The lock itself is
    not tried: that could keep a runner that starts in the same moment from getting it."""
    try:
        pid = int((ROOT / "queue" / "queue.lock").read_text().split()[0])
        with open(f"/proc/{pid}/cmdline", "rb") as f:
            return pid if b"queue.py" in f.read() else None
    except (OSError, ValueError, IndexError):
        return None


def stamp(t):
    u = datetime.datetime.fromtimestamp(t, datetime.timezone.utc)
    return f"{u:%Y-%m-%d %H:%M:%S} UTC | {u.astimezone(IST).strftime('%-d %b %Y, %-I:%M %p IST')}"


def say(msg, live):
    """Printed; with --yes also appended to status/deadman.log."""
    line = f"[{stamp(now())}] {msg}"
    try:
        print(line, flush=True)
    except (OSError, ValueError):                            # a line that cannot be printed must not hold up the shutdown
        pass
    if live:
        try:
            (ROOT / "status").mkdir(parents=True, exist_ok=True)
            with (ROOT / "status" / "deadman.log").open("a") as f:
                f.write(line + "\n")
        except OSError:
            pass


def main(argv=None):
    ap = argparse.ArgumentParser(description="shuts the rented machine down when the Pi has not been seen; a dry run unless --yes")
    ap.add_argument("--yes", action="store_true", help="act (what the systemd timer runs)")
    ap.add_argument("--limit-min", type=float, default=LIMIT_MIN, help="age of PI_SEEN that counts as silence (default 240)")
    ap.add_argument("--wait-min", type=float, default=WAIT_MIN, help="how long the end of the task runner is waited for (default 10)")
    ap.add_argument("--state", action="store_true", help="a dry run that prints one word: armed, silent or off")
    a = ap.parse_args(argv)
    live, t, status = a.yes and not a.state, now(), ROOT / "status"
    old = {}
    try:
        old = json.loads((status / "deadman.json").read_text())
        old = old if isinstance(old, dict) else {}
    except (OSError, ValueError):
        pass
    kept = isinstance(old.get("first_run_epoch"), (int, float)) and not isinstance(old.get("first_run_epoch"), bool)
    try:
        seen = (status / "PI_SEEN").stat().st_mtime
    except OSError:
        seen = None

    def look(first, counted_from):
        """(age of the last sign of the Pi in minutes, what it was counted from, state, sentence)."""
        if seen is not None:
            age, basis = (t - seen) / 60, "age of status/PI_SEEN"
        else:
            age, basis = (t - max(t - uptime_s(), first)) / 60, counted_from
        if (status / "NO_DEADMAN").exists():
            return age, basis, "off", "the file status/NO_DEADMAN exists: the switch is OFF, nothing on this machine ends the billing"
        if age < a.limit_min:
            return age, basis, "armed", f"armed: the Pi was seen {age:.0f} min ago ({basis}); the limit is {a.limit_min:g} min"
        return age, basis, "silent", (f"THE PI HAS NOT BEEN SEEN FOR {age:.0f} MIN ({basis}); the limit is {a.limit_min:g} min: "
                                      "the task queue is stopped and the machine is shut down")
    first = old["first_run_epoch"] if kept else t
    age, basis, state, text = look(first, "no status/PI_SEEN yet: counted from the boot or from the first run of this switch")
    if a.state:
        print(state)
        return 0
    if not live:
        print(f"[{stamp(t)}] DRY RUN: {text}" + (" (with --yes)" if state == "silent" else ""))
        return 0

    def keep():
        try:
            status.mkdir(parents=True, exist_ok=True)
            tmp = status / "deadman.json.tmp"
            tmp.write_text(json.dumps({"state": state, "checked": stamp(t), "epoch": int(t), "pi_seen_age_min": round(age, 1),
                                       "basis": basis, "limit_min": a.limit_min, "first_run_epoch": first}, indent=1) + "\n")
            os.replace(tmp, status / "deadman.json")
            return True
        except OSError:
            return False
    if not keep() and not kept and seen is None:
        # Never seen by the Pi, and the time of this switch's first run can neither be read nor kept (the status folder
        # cannot be written): the count would begin again at every run and the switch would never act.  The boot counts.
        first = 0.0
        age, basis, state, text = look(first, "no status/PI_SEEN, and status/deadman.json can neither be read nor written: "
                                              "counted from the boot")
    if state != old.get("state") or state == "silent":
        say(text, True)
    if state != "silent":
        return 0
    if not machine_ok():
        say(f"this machine is {platform.machine()}, not x86-64 (the rented machine is; the Pi is not): nothing is stopped, "
            "nothing is shut down", True)
        return 3
    stop_file = ROOT / "queue" / "STOP"
    ours = not stop_file.exists()
    rc, text = run(QUEUE + ["stop", "--root", str(ROOT)], 60)
    say(f"queue.py stop: exit code {rc}" + (f" ({text})" if rc else "") + "; waiting for the task runner to end, "
        f"{a.wait_min:g} min at most", True)
    t0 = now()
    while runner_pid() and now() - t0 < a.wait_min * 60:
        sleep(5)
    left = runner_pid()
    say(f"the task runner (pid {left}) is still there after {(now() - t0) / 60:.0f} min: shutting down all the same (systemd ends it)"
        if left else f"the task runner has ended ({now() - t0:.0f} s)", True)
    if ours:
        try:
            stop_file.unlink()
            say("the STOP file of the queue, made by this run, is removed again: the work goes on at the next start of the machine", True)
        except OSError:
            pass
    else:
        say("the STOP file of the queue was there before this run: left in place", True)
    rc, text = run(SHUTDOWN, 60)
    say(f"shutdown -h now: exit code {rc}" + (f" ({text}); the timer tries again in 10 minutes" if rc else ""), True)
    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
