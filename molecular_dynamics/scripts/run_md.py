#!/usr/bin/env python
"""run_md.py: one molecular-dynamics trajectory of a built system, in steps that survive interruptions.

Part of md/ (plan: md/PLAN.md; how to use it: md/RUN.md; every setting: md/protocol.json).
Needs OpenMM 8.1 or newer, NumPy and the Python standard library; nothing else.

  python run_md.py --system DIR --out DIR --run-id X1 [--seed N] [--platform auto|CUDA|OpenCL|CPU]
                   [--threads N] [--stage all|min|equil|prod] [--prod-ns 100] [--segment-ns 5]
                   [--resume] [--benchmark-seconds S] [--equil-scale F] [--protocol FILE]

A run is a chain of steps: min, eq01..eqNN (the equilibration stages of the protocol; after eq03
up to four extension steps eq03x1.., see below), seg0001.. (production segments). Every step
starts from the saved state (XML) of the step before it, in a new OpenMM Context, with seeds that
depend only on the run seed and the name of the step. A step is complete when its record
rec_<step>.json exists; the record is written last. --resume removes the files of a step without
a record and repeats it from the state before it, so an interruption leaves no duplicated and no
missing frame, and on a deterministic platform the repeated step is the step an uninterrupted run
would have made.

Hydration check (protocol key equilibration.hydration_check; md/PLAN.md, entries of 6 Oct 2026,
1:03 AM point 3 and 3:02 AM point 8). At the end of the last stage in which the solute is held
(eq03) the water oxygens near the place of the ring in the crystal pose are counted. If the
ligand is elsewhere and the count is zero, the held part is extended by the steps eq03x1, eq03x2,
.. (at most four), each an ordinary step with its own record, files and seeds, and the count is
made again after each; then the run goes on in any case. The count and the decision are part of
the record of the step, so the chain of steps is fixed by the records alone and --resume finds it
again. The run id is a name only: X0, X1 and any other id are treated alike.

Units in the files: state logs and JSON files name the unit in every column and key (ps, nm,
kJ mol-1, K, g cm-3; keys ending in _A are Angstrom). DCD and PDB files are in Angstrom.

Exit codes: 0 done; 1 unexpected error; 2 wrong command line or input files; 3 the output folder
does not fit the request (not empty without --resume, changed inputs, missing earlier stage);
4 the simulation failed (no usable platform, energy not finite); 130 or 143 stopped by a signal.
The last line printed on success is "DONE <run-id>".
"""
import argparse
import datetime
import glob
import hashlib
import itertools
import json
import math
import os
import platform as py_platform
import signal
import socket
import struct
import subprocess
import sys
import time

import numpy as np

try:
    import openmm
    from openmm import app, unit
except Exception as exc:  # pragma: no cover
    sys.stderr.write(f"ERROR: OpenMM cannot be imported ({exc}).\n")
    sys.exit(2)

PROGRAM = "cqc-sfa-md run_md.py"
MANIFEST_SCHEMA = "cqc-sfa-md run manifest, version 2"      # version 2: key hydration_check, extension steps among the stages
EXIT_OK, EXIT_ERROR, EXIT_INPUT, EXIT_STATE, EXIT_SIM = 0, 1, 2, 3, 4
GROUP_FF, GROUP_EXTRA, GROUP_POSRES = 0, 1, 2
GAS_CONSTANT_KJ_PER_MOL_K = 8.31446261815324e-3
AMU_PER_NM3_TO_G_PER_CM3 = 1.66053906660e-3
CHUNK_STEPS = 500           # longest stretch between two looks at the stop flag
INPUT_FILES = ("system.prmtop", "system.inpcrd", "selections.json", "extra_forces.json")
OPTIONAL_INPUT_FILES = ("reference.pdb", "build.json")
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), "IST")

_STOP = {"signal": None}


class DriverError(Exception):
    def __init__(self, message, code):
        super().__init__(message)
        self.code = code


class Stopped(Exception):
    pass


# ----------------------------------------------------------------------------- small tools

def stamp():
    now = datetime.datetime.now(datetime.timezone.utc)
    ist = now.astimezone(IST)
    hour = ist.hour % 12 or 12
    return {"utc": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "ist": f"{ist.day} {ist.strftime('%b %Y')}, {hour}:{ist.minute:02d}:{ist.second:02d} "
                   f"{'AM' if ist.hour < 12 else 'PM'} IST"}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def file_info(path):
    return {"sha256": sha256_file(path), "bytes": os.path.getsize(path)}


def write_atomic(path, data):
    """Write text or bytes to path so that the file is either complete or absent."""
    tmp = path + ".tmp"
    with open(tmp, "wb" if isinstance(data, bytes) else "w") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    try:                                   # make the new name durable as well
        fd = os.open(os.path.dirname(path) or ".", os.O_RDONLY)
        os.fsync(fd)
        os.close(fd)
    except OSError:
        pass


def write_json(path, obj):
    write_atomic(path, json.dumps(obj, indent=1) + "\n")


def read_json(path, what, code=EXIT_INPUT):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        raise DriverError(f"{what} is missing: {path}", code)
    except json.JSONDecodeError as exc:
        raise DriverError(f"{what} is not valid JSON: {path} ({exc})", code)


def derive_seed(*parts):
    """A seed in 1..2147483646 from the first 4 bytes of the sha256 of the parts joined by '|'."""
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big") % 2147483646 + 1


def run_seed_from_id(run_id):
    return derive_seed("cqc-sfa-md run seed", run_id)


def step_seed(run_seed, tag, purpose):
    return derive_seed("cqc-sfa-md step seed", run_seed, tag, purpose)


def steps_for(interval_ps, dt_fs, what):
    n = interval_ps * 1000.0 / dt_fs
    if n < 1 - 1e-9 or abs(n - round(n)) > 1e-6:
        raise DriverError(f"{what}: {interval_ps} ps is not a whole number of {dt_fs} fs steps", EXIT_INPUT)
    return int(round(n))


def pid_is_run_md(pid):
    """True if a process with this id is alive on this machine and (where /proc tells) is a run_md.py."""
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, ValueError, TypeError):
        return False
    except PermissionError:
        return True
    try:
        with open(f"/proc/{int(pid)}/cmdline", "rb") as f:
            return b"run_md" in f.read()
    except OSError:
        return True


def read_reference_pdb(path):
    """reference.pdb of a built system: the CA atoms in file order (coordinates in nm, residue names) and
    the coordinates (nm) of every other atom by its name."""
    ca, ca_res, others = [], [], {}
    try:
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    name = line[12:16].strip()
                    xyz = [float(line[30:38]) / 10.0, float(line[38:46]) / 10.0, float(line[46:54]) / 10.0]
                    if name == "CA":
                        ca.append(xyz)
                        ca_res.append(line[17:20].strip())
                    else:
                        others.setdefault(name, []).append(xyz)
    except (OSError, ValueError) as exc:
        raise DriverError(f"reference.pdb cannot be read: {path} ({exc})", EXIT_INPUT)
    return np.array(ca, dtype=np.float64).reshape(-1, 3), ca_res, others


def superpose(ref, cur):
    """Least-squares fit (no weights, proper rotation only) of the points ref onto the points cur, both
    (m, 3). Returns q, c_ref, c_cur and the r.m.s. deviation after the fit: a point y of the frame of ref
    is at (y - c_ref) @ q + c_cur in the frame of cur."""
    c_ref, c_cur = ref.mean(axis=0), cur.mean(axis=0)
    u, _s, vt = np.linalg.svd((ref - c_ref).T @ (cur - c_cur))
    q = u @ np.diag([1.0, 1.0, float(np.sign(np.linalg.det(u @ vt)))]) @ vt
    dev = (ref - c_ref) @ q + c_cur - cur
    return q, c_ref, c_cur, float(np.sqrt((dev ** 2).sum(axis=1).mean()))


def hydration_decision(applies, n_water, n_extensions_done, max_extensions):
    """The rule of the hydration check (md/PLAN.md, 6 Oct 2026, 1:03 AM point 3 and 3:02 AM point 8).
    Returns ("extend" or "go on", the reason in words)."""
    if not applies:
        return "go on", "check skipped: the ligand itself is at the site"
    if n_water > 0:
        return "go on", f"{n_water} water oxygen{'s' if n_water != 1 else ''} at the site"
    if n_extensions_done < max_extensions:
        return "extend", (f"no water oxygen at the site after {n_extensions_done} of at most {max_extensions} "
                          "extension steps: the held part is extended")
    return "go on", (f"still no water oxygen at the site after {n_extensions_done} extension steps, the most allowed: "
                     "the run goes on as built, with an empty site")


def on_signal(signum, _frame):
    _STOP["signal"] = signum


def check_stop():
    if _STOP["signal"] is not None:
        raise Stopped()


# ----------------------------------------------------------------------------- prmtop names

def read_prmtop_sections(path, wanted):
    """The raw fields of some %FLAG sections of an Amber prmtop (fixed-width Fortran formats)."""
    out, flag, width, keep = {}, None, None, False
    with open(path) as f:
        for line in f:
            if line.startswith("%FLAG"):
                flag = line.split()[1]
                keep = flag in wanted
                if keep:
                    out[flag] = []
            elif line.startswith("%FORMAT"):
                fmt = line[line.index("(") + 1:line.index(")")]
                letter = next(i for i, ch in enumerate(fmt) if ch.isalpha())
                width = int(fmt[letter + 1:].split(".")[0])
            elif line.startswith("%"):
                continue
            elif keep:
                text = line.rstrip("\n")
                out[flag].extend(text[i:i + width] for i in range(0, len(text), width))
    missing = [w for w in wanted if w not in out]
    if missing:
        raise DriverError(f"{path} lacks the section(s) {', '.join(missing)}", EXIT_INPUT)
    return out


def prmtop_names(path):
    """Atom names, residue names and the zero-based residue index of every atom, as written by tleap."""
    sec = read_prmtop_sections(path, ("ATOM_NAME", "RESIDUE_LABEL", "RESIDUE_POINTER"))
    atom_names = [s.strip() for s in sec["ATOM_NAME"]]
    res_names = [s.strip() for s in sec["RESIDUE_LABEL"]]
    first = [int(s) - 1 for s in sec["RESIDUE_POINTER"]]
    res_of_atom = np.zeros(len(atom_names), dtype=np.int64)
    for r, start in enumerate(first):
        end = first[r + 1] if r + 1 < len(first) else len(atom_names)
        res_of_atom[start:end] = r
    return atom_names, res_names, res_of_atom


# ----------------------------------------------------------------------------- geometry, files

def box_lengths_angles(box):
    """Lengths (nm) and angles (degrees) of a cell given as rows a, b, c (nm)."""
    a, b, c = (np.linalg.norm(box[i]) for i in range(3))
    alpha = math.degrees(math.acos(np.dot(box[1], box[2]) / (b * c)))
    beta = math.degrees(math.acos(np.dot(box[0], box[2]) / (a * c)))
    gamma = math.degrees(math.acos(np.dot(box[0], box[1]) / (a * b)))
    return a, b, c, alpha, beta, gamma


def min_image(d, box):
    """Shortest periodic images of the vectors d (n, 3; nm) in a cell with rows a, b, c (reduced form)."""
    d = d - np.outer(np.round(d[:, 2] / box[2, 2]), box[2])
    d = d - np.outer(np.round(d[:, 1] / box[1, 1]), box[1])
    d = d - np.outer(np.round(d[:, 0] / box[0, 0]), box[0])
    best, best_r2 = d.copy(), np.einsum("ij,ij->i", d, d)
    for i, j, k in itertools.product((-1, 0, 1), repeat=3):
        if i == j == k == 0:
            continue
        cand = d + i * box[0] + j * box[1] + k * box[2]
        r2 = np.einsum("ij,ij->i", cand, cand)
        closer = r2 < best_r2
        best[closer], best_r2[closer] = cand[closer], r2[closer]
    return best


class SolventWrapper:
    """Moves every residue outside the solute by cell vectors to the image nearest to the solute's centre."""

    def __init__(self, n_atoms, solute, res_of_atom):
        mask = np.ones(n_atoms, dtype=bool)
        mask[solute] = False
        self.solute = np.asarray(solute, dtype=np.int64)
        self.solvent = np.nonzero(mask)[0]
        _, first, self.group = np.unique(res_of_atom[self.solvent], return_index=True, return_inverse=True)
        self.anchor = self.solvent[first]          # first atom of each residue (the oxygen of a water)

    def __call__(self, pos, box):
        out = pos.copy()
        if len(self.solvent):
            d = pos[self.anchor] - pos[self.solute].mean(axis=0)
            out[self.solvent] += (min_image(d, box) - d)[self.group]
        return out


class DcdWriter:
    """DCD file in the CHARMM layout that OpenMM writes, but without a time stamp in the header,
    so that the same trajectory gives the same bytes. Coordinates in Angstrom, single precision."""

    def __init__(self, path, n_atoms, first_step, interval_steps, dt_ps, titles):
        self.path, self.n_atoms, self.first, self.interval, self.count = path, n_atoms, first_step, interval_steps, 0
        self.f = open(path, "wb")
        head = struct.pack("<i4c9if", 84, b"C", b"O", b"R", b"D", 0, first_step, interval_steps,
                           0, 0, 0, 0, 0, 0, dt_ps / 0.04888821)
        head += struct.pack("<13i", 1, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
        for line in (list(titles) + ["", ""])[:2]:
            head += struct.pack("<80s", line.encode("ascii", "replace")[:80].ljust(80))
        head += struct.pack("<4i", 164, 4, n_atoms, 4)
        self.f.write(head)
        self.f.flush()          # a new file is a valid header-only DCD from the start (a run in progress can be analysed)

    def write(self, pos_nm, box_nm):
        if pos_nm.shape != (self.n_atoms, 3):
            raise ValueError("wrong number of atoms for " + self.path)
        if not np.isfinite(pos_nm).all():
            raise DriverError("a coordinate is not finite; the simulation has failed", EXIT_SIM)
        self.count += 1
        f = self.f
        f.seek(8)
        f.write(struct.pack("<i", self.count))
        f.seek(20)
        f.write(struct.pack("<i", self.first + (self.count - 1) * self.interval))
        f.seek(0, os.SEEK_END)
        a, b, c, alpha, beta, gamma = box_lengths_angles(box_nm)
        f.write(struct.pack("<i6di", 48, 10 * a, math.cos(math.radians(gamma)), 10 * b,
                            math.cos(math.radians(beta)), math.cos(math.radians(alpha)), 10 * c, 48))
        marker = struct.pack("<i", 4 * self.n_atoms)
        for dim in range(3):
            f.write(marker)
            f.write((10.0 * pos_nm[:, dim]).astype("<f4").tobytes())
            f.write(marker)
        f.flush()

    def close(self):
        if not self.f.closed:
            self.f.flush()
            os.fsync(self.f.fileno())
            self.f.close()


def pdb_text(atoms, pos_nm, box_nm, remarks):
    """PDB text (Angstrom) for the atoms given as (prmtop index, atom name, residue name, residue
    index, element). The serial number is the prmtop index + 1 and the residue number the prmtop
    residue index + 1 (both cut to the width of the PDB columns)."""
    lines = [f"REMARK   {r}" for r in remarks]
    a, b, c, alpha, beta, gamma = box_lengths_angles(box_nm)
    lines.append("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1" % (10 * a, 10 * b, 10 * c, alpha, beta, gamma))
    for row, (index, name, resname, res, element) in enumerate(atoms):
        shown = name if len(name) >= 4 or len(element) == 2 else " " + name
        x, y, z = 10.0 * pos_nm[row]
        lines.append("ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                     % ((index + 1) % 100000, shown[:4], resname[:3], (res + 1) % 10000, x, y, z, element.upper()))
    lines.append("END")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------- platform

def host_info():
    info = {"hostname": socket.gethostname(), "os": py_platform.platform(), "machine": py_platform.machine(),
            "cpu_count": os.cpu_count(), "python": sys.version.split()[0], "numpy": np.__version__,
            "openmm": openmm.Platform.getOpenMMVersion(),
            "openmm_git_revision": getattr(openmm.version, "git_revision", None),
            "openmm_platforms": [openmm.Platform.getPlatform(i).getName()
                                 for i in range(openmm.Platform.getNumPlatforms())],
            "openmm_plugin_load_failures": list(openmm.Platform.getPluginLoadFailures())}
    try:
        res = subprocess.run(["nvidia-smi", "--query-gpu=index,name,driver_version,memory.total",
                              "--format=csv,noheader"], capture_output=True, text=True, timeout=20)
        info["nvidia_smi"] = [ln.strip() for ln in res.stdout.splitlines() if ln.strip()] if res.returncode == 0 else None
    except Exception:
        info["nvidia_smi"] = None
    return info


class PlatformChoice:
    """Picks the OpenMM platform once and makes every Context of the run with the same properties."""

    def __init__(self, requested, threads, proto, say):
        self.requested, self.attempts = requested, []
        wanted = ["CUDA", "OpenCL", "CPU"] if requested == "auto" else [requested]
        for name in wanted:
            try:
                plat = openmm.Platform.getPlatformByName(name)
                props = self._properties(plat, threads, proto)
                probe = openmm.System()
                probe.addParticle(1.0)
                ctx = openmm.Context(probe, openmm.VerletIntegrator(0.001), plat, props)
                del ctx
            except Exception as exc:
                self.attempts.append({"platform": name, "usable": False, "error": str(exc)})
                continue
            self.attempts.append({"platform": name, "usable": True})
            self.platform, self.name, self.properties_set = plat, name, props
            break
        else:
            tried = "; ".join(f"{a['platform']}: {a['error']}" for a in self.attempts)
            raise DriverError(f"no usable OpenMM platform for --platform {requested} ({tried})", EXIT_SIM)
        if requested == "auto" and self.name == "CPU":
            say("WARNING: --platform auto found no GPU platform and uses the CPU platform.")
        self.properties_actual, self.gpu_name = {}, None

    @staticmethod
    def _properties(plat, threads, proto):
        names, props = set(plat.getPropertyNames()), {}
        if plat.getName() in ("CUDA", "OpenCL") and "Precision" in names:
            props["Precision"] = proto["platform"]["gpu_precision"]
        if "DeterministicForces" in names and proto["platform"]["deterministic_forces"]:
            props["DeterministicForces"] = "true"
        if plat.getName() == "CPU" and threads and "Threads" in names:
            props["Threads"] = str(threads)
        return props

    def context(self, system, integrator):
        ctx = openmm.Context(system, integrator, self.platform, self.properties_set)
        self.properties_actual = {n: self.platform.getPropertyValue(ctx, n) for n in self.platform.getPropertyNames()}
        self.gpu_name = self.properties_actual.get("DeviceName")
        return ctx

    def record(self):
        return {"requested": self.requested, "name": self.name, "properties_set": self.properties_set,
                "properties_actual": self.properties_actual, "gpu_name": self.gpu_name,
                "deterministic_forces_set": self.properties_set.get("DeterministicForces") == "true",
                "precision_set": self.properties_set.get("Precision", "double (CPU platform)" if self.name == "CPU" else None),
                "repeatability": {
                    "CPU": "bit for bit only with one thread and DeterministicForces (measured on the Raspberry Pi, md/RUN.md)",
                    "CUDA": "DeterministicForces is set; whether a repeated segment is bit for bit the same has to be "
                            "checked on the machine (pre-flight test of md/RUN.md)",
                    "OpenCL": "no deterministic mode; a repeated segment is a new, equally valid trajectory"}[self.name],
                "attempts": self.attempts}


# ----------------------------------------------------------------------------- the run

class Run:
    def __init__(self, args):
        self.args = args
        self.out = os.path.abspath(args.out)
        self.sysdir = os.path.abspath(args.system)
        self.log_handle = None
        self.invocation = None
        self.manifest = None
        self.lock = None
        self.absent = set()

    def take_lock(self):
        lock = self.path("run_md.lock")
        me = {"pid": os.getpid(), "host": socket.gethostname(), "started": stamp()}
        for _attempt in (1, 2):
            try:
                fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                try:
                    with open(lock) as f:
                        other = json.load(f)
                except Exception:
                    other = {}
                if other.get("host") == me["host"] and pid_is_run_md(other.get("pid")):
                    raise DriverError(f"another run_md.py (pid {other.get('pid')}) is working in {self.out}; "
                                      f"if that is wrong, delete {lock}", EXIT_STATE)
                os.remove(lock)                 # left by a process that was killed, or on another machine
                continue
            os.write(fd, (json.dumps(me) + "\n").encode())
            os.close(fd)
            self.lock = lock
            return
        raise DriverError(f"the lock {lock} could not be taken", EXIT_STATE)

    def release_lock(self):
        if self.lock and os.path.isfile(self.lock):
            os.remove(self.lock)
        self.lock = None

    # ---- messages
    def say(self, text):
        line = f"[{stamp()['ist']}] {text}"
        print(line, flush=True)
        if self.log_handle:
            self.log_handle.write(line + "\n")
            self.log_handle.flush()

    def path(self, name):
        return os.path.join(self.out, name)

    # ---- inputs
    def load_inputs(self):
        args = self.args
        if not os.path.isdir(self.sysdir):
            raise DriverError(f"--system is not a folder: {self.sysdir}", EXIT_INPUT)
        for name in INPUT_FILES:
            if not os.path.isfile(os.path.join(self.sysdir, name)):
                raise DriverError(f"the built system lacks {name}: {self.sysdir}", EXIT_INPUT)
        self.protocol_path = os.path.abspath(args.protocol)
        self.proto = read_json(self.protocol_path, "the protocol file")
        for key in ("purpose", "system", "integrator", "barostat", "minimisation", "equilibration",
                    "production", "output", "platform"):
            if key not in self.proto:
                raise DriverError(f"the protocol file has no '{key}': {self.protocol_path}", EXIT_INPUT)
        if self.proto["system"]["nonbonded_method"] != "PME":
            raise DriverError("this driver supports nonbonded_method PME only", EXIT_INPUT)
        if not 0.0 < args.equil_scale <= 1.0:
            raise DriverError("--equil-scale must be above 0 and at most 1", EXIT_INPUT)

        self.inputs = {name: file_info(os.path.join(self.sysdir, name)) for name in INPUT_FILES}
        for name in OPTIONAL_INPUT_FILES:
            p = os.path.join(self.sysdir, name)
            self.inputs[name] = file_info(p) if os.path.isfile(p) else None
        self.inputs["protocol"] = dict(file_info(self.protocol_path), path=self.protocol_path)
        self.protocol_kept = self.keep_protocol_of_the_run()
        self.script_info = dict(file_info(os.path.abspath(__file__)), path=os.path.abspath(__file__))

        try:
            self.prmtop = app.AmberPrmtopFile(os.path.join(self.sysdir, "system.prmtop"))
            inpcrd = app.AmberInpcrdFile(os.path.join(self.sysdir, "system.inpcrd"))
        except Exception as exc:
            raise DriverError(f"system.prmtop or system.inpcrd cannot be read: {exc}", EXIT_INPUT)
        self.n_atoms = self.prmtop.topology.getNumAtoms()
        self.pos0 = np.array(inpcrd.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=np.float64)
        if self.pos0.shape != (self.n_atoms, 3):
            raise DriverError(f"system.inpcrd has {len(self.pos0)} atoms, system.prmtop {self.n_atoms}", EXIT_INPUT)
        if inpcrd.boxVectors is None:
            raise DriverError("system.inpcrd holds no periodic box", EXIT_INPUT)
        self.box0 = np.array([[float(c) for c in vec.value_in_unit(unit.nanometer)] for vec in inpcrd.boxVectors],
                             dtype=np.float64)
        self.atom_names, self.res_names, self.res_of_atom = prmtop_names(os.path.join(self.sysdir, "system.prmtop"))
        if len(self.atom_names) != self.n_atoms:
            raise DriverError("ATOM_NAME of system.prmtop does not match its number of atoms", EXIT_INPUT)
        self.elements = [a.element.symbol if a.element is not None else "X" for a in self.prmtop.topology.atoms()]

        sel = read_json(os.path.join(self.sysdir, "selections.json"), "selections.json")
        for key in ("zn", "lig_heavy", "protein_heavy", "solute"):
            if key not in sel:
                raise DriverError(f"selections.json has no '{key}'", EXIT_INPUT)
        self.sel = sel
        used = [sel["zn"]] + list(sel["lig_heavy"].values()) + list(sel["protein_heavy"]) + list(sel["solute"])
        used += list(sel.get("crystal_water_O", [])) + list(sel.get("water_O", [])) + list(sel.get("protein_CA", []))
        if min(used) < 0 or max(used) >= self.n_atoms:
            raise DriverError("selections.json holds an atom index outside the prmtop", EXIT_INPUT)
        if self.elements[sel["zn"]] != "Zn":
            raise DriverError(f"selections.json 'zn' points at element {self.elements[sel['zn']]}, not Zn", EXIT_INPUT)
        if any(self.elements[i] == "H" for i in list(sel["lig_heavy"].values()) + list(sel["protein_heavy"])):
            raise DriverError("selections.json lists a hydrogen among lig_heavy or protein_heavy", EXIT_INPUT)
        self.solute = np.array(sel["solute"], dtype=np.int64)
        if len(set(sel["solute"])) != len(self.solute) or list(self.solute) != sorted(sel["solute"]):
            raise DriverError("selections.json 'solute' must be in topology order, without repeats", EXIT_INPUT)

        heavy = sorted(set(sel["protein_heavy"]) | {sel["zn"]} | set(sel["lig_heavy"].values()))
        crystal_o = sorted(set(sel.get("crystal_water_O", [])))
        self.atom_sets = {
            "solute_heavy": heavy,
            "solute_heavy_and_crystal_water_O": sorted(set(heavy) | set(crystal_o)),
            "protein_backbone": [i for i in sorted(set(sel["protein_heavy"]))
                                 if self.atom_names[i] in ("N", "CA", "C", "O")],
            "none": [],
        }
        self.crystal_water_listed = len(crystal_o)

        extra = read_json(os.path.join(self.sysdir, "extra_forces.json"), "extra_forces.json")
        if "harmonic_bonds" not in extra or not isinstance(extra["harmonic_bonds"], list):
            raise DriverError("extra_forces.json needs a list 'harmonic_bonds' (empty if there are none)", EXIT_INPUT)
        self.extra_bonds = []
        for n, b in enumerate(extra["harmonic_bonds"]):
            try:
                i, j, r0, k = int(b["i"]), int(b["j"]), float(b["r0_nm"]), float(b["k_kj_per_mol_nm2"])
            except (KeyError, TypeError, ValueError):
                raise DriverError(f"harmonic bond {n} of extra_forces.json needs i, j, r0_nm and k_kj_per_mol_nm2", EXIT_INPUT)
            if not (0 <= i < self.n_atoms and 0 <= j < self.n_atoms and i != j and r0 > 0 and k >= 0):
                raise DriverError(f"harmonic bond {n} of extra_forces.json has impossible values", EXIT_INPUT)
            self.extra_bonds.append({"i": i, "j": j, "r0_nm": r0, "k_kj_per_mol_nm2": k, "label": str(b.get("label", "")),
                                     "atoms": f"{self.atom_label(i)} - {self.atom_label(j)}"})

        self.wrapper = SolventWrapper(self.n_atoms, self.solute, self.res_of_atom)
        self.hyd = self.load_hydration_check()
        self.build_info = None
        bpath = os.path.join(self.sysdir, "build.json")
        if os.path.isfile(bpath) and os.path.getsize(bpath) < 200000:
            try:
                with open(bpath) as f:
                    self.build_info = json.load(f)
            except Exception:
                self.build_info = "build.json is not valid JSON"

    def atom_label(self, i):
        r = int(self.res_of_atom[i])
        return f"{self.res_names[r]}{r + 1}:{self.atom_names[i]}"

    def keep_protocol_of_the_run(self):
        """A run that began with a protocol without the block equilibration.hydration_check keeps that
        protocol (the copy in its manifest) if the protocol file has since gained this block and nothing
        else: the run continues as it began, without the check. Any other change of the protocol file is
        refused by check_manifest, as before. Returns a note for the manifest, or None."""
        try:
            with open(self.path("manifest.json")) as f:
                man = json.load(f)
            old, old_info = man["protocol"], man["inputs"]["protocol"]
        except Exception:
            return None
        if old_info.get("sha256") == self.inputs["protocol"]["sha256"]:
            return None
        if not isinstance(old.get("equilibration"), dict) or "hydration_check" in old["equilibration"]:
            return None
        now = json.loads(json.dumps(self.proto))
        if not isinstance(now.get("equilibration"), dict) or now["equilibration"].pop("hydration_check", None) is None:
            return None
        if now != old:
            return None
        note = {"protocol_file_now": dict(self.inputs["protocol"]),
                "protocol_file_of_the_run": dict(old_info),
                "difference": "the file now holds the block equilibration.hydration_check and is otherwise the same",
                "consequence": "the run keeps the protocol recorded in its manifest: no hydration check, no extension step"}
        self.proto = old
        self.inputs["protocol"] = dict(old_info)
        return note

    def load_hydration_check(self):
        """The settings of the hydration check (protocol key equilibration.hydration_check) together with
        what it needs of the built system; None if the protocol does not ask for the check."""
        block = self.proto["equilibration"].get("hydration_check")
        if block is None:
            return None
        stages, sel = self.proto["equilibration"]["stages"], self.sel
        try:
            after, names = int(block["after_stage"]), [str(n) for n in block["site_atoms"]]
            cutoff, apart = float(block["water_cutoff_nm"]), float(block["applies_if_ligand_mean_distance_above_nm"])
            ext_ps, max_ext = float(block["extension_length_ps"]), int(block["max_extensions"])
        except (KeyError, TypeError, ValueError):
            raise DriverError("equilibration.hydration_check of the protocol needs after_stage, site_atoms, water_cutoff_nm, "
                              "applies_if_ligand_mean_distance_above_nm, extension_length_ps and max_extensions", EXIT_INPUT)
        if not (1 <= after <= len(stages) and names and len(set(names)) == len(names) and cutoff > 0 and apart >= 0
                and ext_ps > 0 and 0 <= max_ext <= 99):
            raise DriverError("equilibration.hydration_check of the protocol holds an impossible value", EXIT_INPUT)
        ref_path = os.path.join(self.sysdir, "reference.pdb")
        if not os.path.isfile(ref_path):
            raise DriverError(f"the hydration check of the protocol needs reference.pdb in the built system: {self.sysdir}", EXIT_INPUT)
        for key in ("protein_CA", "water_O"):
            if not sel.get(key):
                raise DriverError(f"the hydration check of the protocol needs '{key}' in selections.json", EXIT_INPUT)
        ca = [int(i) for i in sel["protein_CA"]]
        ref_ca, ref_res, others = read_reference_pdb(ref_path)
        if len(ca) < 3 or len(ref_ca) != len(ca):
            raise DriverError(f"reference.pdb holds {len(ref_ca)} CA atoms, selections.json 'protein_CA' {len(ca)} "
                              "(the same number, at least 3, is needed)", EXIT_INPUT)
        wrong = [i for i in ca if self.atom_names[i] != "CA"]
        if wrong:
            raise DriverError(f"selections.json 'protein_CA' points at an atom that is not named CA ({self.atom_label(wrong[0])})", EXIT_INPUT)
        bad = [n for n in names if len(others.get(n, [])) != 1 or n not in sel["lig_heavy"]]
        if bad:
            raise DriverError(f"the hydration check needs the atom(s) {', '.join(bad)} exactly once in reference.pdb and in "
                              "selections.json 'lig_heavy'", EXIT_INPUT)
        differ = [f"{self.res_names[self.res_of_atom[i]]} in system.prmtop, {r} in reference.pdb (CA {k + 1})"
                  for k, (i, r) in enumerate(zip(ca, ref_res)) if self.res_names[self.res_of_atom[i]] != r]
        rms_start = superpose(ref_ca, self.pos0[ca])[3]
        base = stages[after - 1]
        temps = base["temperature_K"]
        ext_stage = dict(base, length_ps=ext_ps, temperature_K=float(temps[-1]) if isinstance(temps, list) else float(temps))
        ext_stage.pop("temperature_ramp", None)
        later = stages[after]["name"] if after < len(stages) else "production"
        return {"block": block, "after_stage": after, "base_tag": f"eq{after:02d}", "site_atoms": names,
                "cutoff_nm": cutoff, "apart_nm": apart, "extension_ps": ext_ps, "max_extensions": max_ext,
                "extension_stage": ext_stage, "next_tag": f"eq{after + 1:02d}" if after < len(stages) else "seg0001",
                "next_name": later,
                "ca": np.array(ca, dtype=np.int64), "ref_ca_nm": ref_ca,
                "ref_site_nm": np.array([others[n][0] for n in names], dtype=np.float64),
                "lig": np.array([int(sel["lig_heavy"][n]) for n in names], dtype=np.int64),
                "water": np.array(sel["water_O"], dtype=np.int64),
                "reference_guard": {
                    "n_ca": len(ca),
                    "ca_fit_rmsd_of_system_inpcrd_on_reference_pdb_A": round(rms_start * 10.0, 4),
                    "n_ca_residue_names_that_differ": len(differ), "ca_residue_names_that_differ_first_5": differ[:5]}}

    # ---- OpenMM objects
    def build_system(self, restraint_set, restraint_k, barostat_seed):
        """The System of one step. Force groups: 0 force field, 1 extra bonds, 2 position restraints."""
        p = self.proto["system"]
        constraints = {"HBonds": app.HBonds, "AllBonds": app.AllBonds, "HAngles": app.HAngles, "None": None}[p["constraints"]]
        system = self.prmtop.createSystem(
            nonbondedMethod=app.PME, nonbondedCutoff=p["nonbonded_cutoff_nm"] * unit.nanometer,
            constraints=constraints, rigidWater=bool(p["rigid_water"]), removeCMMotion=bool(p["remove_cm_motion"]),
            hydrogenMass=p["hydrogen_mass_amu"] * unit.amu, ewaldErrorTolerance=p["ewald_error_tolerance"])
        system.setDefaultPeriodicBoxVectors(*[openmm.Vec3(*row) for row in self.box0])
        for force in system.getForces():
            force.setForceGroup(GROUP_FF)
            if isinstance(force, openmm.NonbondedForce):
                force.setUseDispersionCorrection(bool(p["long_range_dispersion_correction"]))
                force.setUseSwitchingFunction(bool(p["lennard_jones_switching_function"]))
            elif isinstance(force, openmm.CustomNonbondedForce):
                force.setUseLongRangeCorrection(bool(p["long_range_dispersion_correction"]))
        if self.extra_bonds:
            bonds = openmm.HarmonicBondForce()          # U = 0.5 k (r - r0)^2
            bonds.setUsesPeriodicBoundaryConditions(True)
            bonds.setForceGroup(GROUP_EXTRA)
            for b in self.extra_bonds:
                bonds.addBond(b["i"], b["j"], b["r0_nm"], b["k_kj_per_mol_nm2"])
            system.addForce(bonds)
        atoms = self.atom_sets[restraint_set]
        if atoms and restraint_k > 0:
            posres = openmm.CustomExternalForce("0.5*k_posres*periodicdistance(x,y,z,x0,y0,z0)^2")
            for name in ("k_posres", "x0", "y0", "z0"):
                posres.addPerParticleParameter(name)
            for i in atoms:
                posres.addParticle(int(i), [float(restraint_k)] + [float(v) for v in self.pos0[i]])
            posres.setForceGroup(GROUP_POSRES)
            system.addForce(posres)
        if barostat_seed is not None:
            b = self.proto["barostat"]
            baro = openmm.MonteCarloBarostat(b["pressure_bar"] * unit.bar, b["temperature_K"] * unit.kelvin,
                                             int(b["frequency_steps"]))
            baro.setRandomNumberSeed(int(barostat_seed))
            system.addForce(baro)
        return system

    def make_context(self, system, dt_fs, temperature_K, integrator_seed):
        g = self.proto["integrator"]
        integ = openmm.LangevinMiddleIntegrator(temperature_K * unit.kelvin, g["friction_per_ps"] / unit.picosecond,
                                                dt_fs * unit.femtosecond)
        integ.setConstraintTolerance(self.proto["system"]["constraint_tolerance"])
        integ.setRandomNumberSeed(int(integrator_seed))
        try:
            ctx = self.plat.context(system, integ)
        except Exception as exc:
            raise DriverError(f"the {self.plat.name} platform could not make a Context: {exc}", EXIT_SIM)
        return ctx, integ

    def describe_system(self):
        system = self.build_system("none", 0.0, 1)
        masses = np.array([system.getParticleMass(i).value_in_unit(unit.amu) for i in range(self.n_atoms)])
        self.total_mass_amu = float(masses.sum())
        n_massive = int((masses > 0).sum())
        n_cons = system.getNumConstraints()
        has_cmm = any(isinstance(f, openmm.CMMotionRemover) for f in system.getForces())
        self.dof = 3 * n_massive - n_cons - (3 if has_cmm else 0)
        is_h = np.array([e == "H" for e in self.elements])
        water = np.array([a.residue.name == "HOH" for a in self.prmtop.topology.atoms()])
        hm = self.proto["system"]["hydrogen_mass_amu"]
        a, b, c, alpha, beta, gamma = box_lengths_angles(self.box0)
        volume = float(abs(np.linalg.det(self.box0)))
        return {
            "n_atoms": self.n_atoms, "n_residues": len(self.res_names), "n_solute_atoms": int(len(self.solute)),
            "n_water_molecules": int((water & ~is_h).sum()),
            "total_mass_amu": round(self.total_mass_amu, 4),
            "hydrogen_masses": {
                "n_hydrogens_not_in_water": int((is_h & ~water).sum()),
                "n_with_repartitioned_mass": int((is_h & ~water & (np.abs(masses - hm) < 1e-9)).sum()),
                "repartitioned_mass_amu": hm,
                "water_hydrogen_masses_amu": sorted({round(float(m), 4) for m in masses[is_h & water]}),
                "lightest_heavy_atom_amu": round(float(masses[~is_h].min()), 4),
            },
            "n_constraints": n_cons, "degrees_of_freedom": self.dof,
            "forces_in_production_system": [f"{type(f).__name__} (group {f.getForceGroup()})" for f in system.getForces()],
            "start_box_vectors_nm": [[round(float(v), 7) for v in row] for row in self.box0],
            "start_box_lengths_nm_angles_deg": [round(float(v), 5) for v in (a, b, c, alpha, beta, gamma)],
            "start_volume_nm3": round(volume, 4),
            "start_density_g_per_cm3": round(self.total_mass_amu * AMU_PER_NM3_TO_G_PER_CM3 / volume, 5),
            "restraint_atom_sets_n_atoms": {k: len(v) for k, v in self.atom_sets.items()},
            "crystal_water_O_listed_in_selections": self.crystal_water_listed,
        }

    # ---- states
    @staticmethod
    def group_energy(ctx, group):
        return ctx.getState(getEnergy=True, groups={group}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    def save_state(self, ctx, name):
        st = ctx.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
        write_atomic(self.path(name), openmm.XmlSerializer.serialize(st))
        return st

    def read_state(self, name):
        with open(self.path(name)) as f:
            st = openmm.XmlSerializer.deserialize(f.read())
        pos = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=np.float64)
        vel = np.array(st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond), dtype=np.float64)
        box = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer), dtype=np.float64)
        if pos.shape != (self.n_atoms, 3):
            raise DriverError(f"{name} does not belong to this system (number of atoms)", EXIT_STATE)
        return pos, vel, box

    @staticmethod
    def set_state(ctx, pos, vel, box):
        ctx.setPeriodicBoxVectors(*[openmm.Vec3(*row) for row in box])
        ctx.setPositions(pos)
        if vel is not None:
            ctx.setVelocities(vel)

    # ---- monitoring (reported, never used to steer the run)
    def monitor(self, pos, box):
        sel, out = self.sel, {}

        def dist(i, j):
            return float(np.linalg.norm(min_image((pos[i] - pos[j])[None, :], box)[0])) * 10.0

        lig = sel["lig_heavy"]
        pairs = [("zn_to_lig_N2_A", sel["zn"], lig.get("N2")),
                 ("lig_N2_to_thr199_OG1_A", lig.get("N2"), sel.get("thr199_OG1")),
                 ("lig_N1_to_phe131_CZ_A", lig.get("N1"), sel.get("phe131_CZ")),
                 ("lig_N1_to_asn67_OD1_A", lig.get("N1"), sel.get("asn67_OD1")),
                 ("zn_to_his94_NE2_A", sel["zn"], sel.get("his94_NE2")),
                 ("zn_to_his96_NE2_A", sel["zn"], sel.get("his96_NE2")),
                 ("zn_to_his119_ND1_A", sel["zn"], sel.get("his119_ND1"))]
        for key, i, j in pairs:
            if i is not None and j is not None:
                out[key] = round(dist(i, j), 3)
        if sel.get("water_O"):
            cutoff = self.proto["output"]["monitor_water_cutoff_nm"]
            ligpos = pos[list(lig.values())]
            centre = ligpos.mean(axis=0)
            d = min_image(pos[sel["water_O"]] - centre, box)
            near = d[np.einsum("ij,ij->i", d, d) < (cutoff + 1.0) ** 2] + centre
            if len(near):
                dmin = np.sqrt(((near[:, None, :] - ligpos[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
                out[f"n_water_O_within_{10 * cutoff:.1f}_A_of_ligand_heavy_atoms"] = int((dmin < cutoff).sum())
            else:
                out[f"n_water_O_within_{10 * cutoff:.1f}_A_of_ligand_heavy_atoms"] = 0
        ref = self.pos0[self.atom_sets["solute_heavy"]]
        now = pos[self.atom_sets["solute_heavy"]]
        shift = (now - ref) - (now - ref).mean(axis=0)
        out["solute_heavy_rms_displacement_from_start_A"] = round(float(np.sqrt((shift ** 2).sum(axis=1).mean())) * 10.0, 3)
        out["solute_heavy_rms_displacement_note"] = "translation removed, no rotational fit"
        if self.extra_bonds:
            out["extra_bond_distances_A"] = [round(dist(b["i"], b["j"]), 3) for b in self.extra_bonds]
        if self.hyd:
            r = self.ring_site_numbers(pos, box)
            out["crystal_ring_site"] = {
                f"n_water_O_within_{10 * self.hyd['cutoff_nm']:.1f}_A": r["n_water_O"],
                "nearest_water_O_A": round(r["nearest_water_O_nm"] * 10.0, 3),
                "ligand_mean_distance_from_site_A": round(r["ligand_mean_nm"] * 10.0, 3)}
        return out

    # ---- hydration check (md/PLAN.md, 6 Oct 2026, 1:03 AM point 3 and 3:02 AM point 8)
    def ring_site_numbers(self, pos, box):
        """Water at the place the ring has in the crystal pose. The positions of the site atoms in
        reference.pdb are carried into the frame of pos by the least-squares fit of the CA atoms (the
        protein is whole in pos: coordinates as integrated); every distance to them is a minimum-image
        distance in the cell box. pos and box in nm."""
        h = self.hyd
        q, c_ref, c_cur, rms = superpose(h["ref_ca_nm"], pos[h["ca"]])
        site = (h["ref_site_nm"] - c_ref) @ q + c_cur
        lig = np.linalg.norm(min_image(pos[h["lig"]] - site, box), axis=1)
        dist = np.stack([np.linalg.norm(min_image(pos[h["water"]] - point, box), axis=1) for point in site])
        nearest = dist.min(axis=0)                       # of every water oxygen to the nearest site position
        inside = np.nonzero(nearest <= h["cutoff_nm"])[0]
        return {"n_water_O": int(len(inside)), "water_O": [int(h["water"][i]) for i in inside],
                "water_O_nm": [float(nearest[i]) for i in inside],
                "nearest_water_O_nm": float(nearest.min()), "nearest_water_O": int(h["water"][int(nearest.argmin())]),
                "nearest_per_site_nm": [float(v) for v in dist.min(axis=1)],
                "ligand_per_atom_nm": [float(v) for v in lig], "ligand_mean_nm": float(lig.mean()),
                "ca_fit_rmsd_nm": rms}

    def hydration_record(self, pos, box, tag, n_done):
        """The count at the end of step tag (the stage after which the check is made, or an extension
        step; n_done extension steps are complete) and what follows from it. Part of the record of the
        step: the chain of steps is read from these records."""
        h, r = self.hyd, self.ring_site_numbers(pos, box)
        applies = r["ligand_mean_nm"] > h["apart_nm"]
        decision, reason = hydration_decision(applies, r["n_water_O"], n_done, h["max_extensions"])
        return {
            "made_on": f"the coordinates at the end of this step (state_{tag}.xml)",
            "site_atoms": h["site_atoms"],
            "site_positions": "reference.pdb, carried into the frame of the coordinates by the fit of the protein CA atoms",
            "ca_fit": {"n_ca": int(len(h["ca"])), "rmsd_A": round(r["ca_fit_rmsd_nm"] * 10.0, 4)},
            "ligand_distance_from_site_A": {n: round(v * 10.0, 3) for n, v in zip(h["site_atoms"], r["ligand_per_atom_nm"])},
            "ligand_mean_distance_from_site_A": round(r["ligand_mean_nm"] * 10.0, 3),
            "applies_if_ligand_mean_distance_above_A": round(h["apart_nm"] * 10.0, 6),
            "applies": bool(applies), "skipped": not applies,
            "water_cutoff_A": round(h["cutoff_nm"] * 10.0, 6),
            "n_water_O": r["n_water_O"],
            "water_O": [{"atom_index": i, "atom": self.atom_label(i), "distance_A": round(d * 10.0, 3)}
                        for i, d in zip(r["water_O"], r["water_O_nm"])],
            "nearest_water_O": {"atom_index": r["nearest_water_O"], "atom": self.atom_label(r["nearest_water_O"]),
                                "distance_A": round(r["nearest_water_O_nm"] * 10.0, 3)},
            "nearest_water_O_per_site_atom_A": {n: round(v * 10.0, 3) for n, v in zip(h["site_atoms"], r["nearest_per_site_nm"])},
            "n_extensions_done": int(n_done), "max_extensions": h["max_extensions"],
            "decision": decision, "reason": reason,
            "next_step": f"{h['base_tag']}x{n_done + 1}" if decision == "extend" else h["next_tag"],
        }

    # ---- bookkeeping of steps
    def plan(self, n_segments):
        """The chain of steps as far as it is known. The extension steps of the hydration check follow
        from the records: step eqNNx<n> is in the chain if the record of the step before it says
        "extend". While that record does not exist the chain goes on with the next stage."""
        steps = [{"tag": "min", "kind": "minimisation"}]
        for k, stage in enumerate(self.proto["equilibration"]["stages"], start=1):
            steps.append({"tag": f"eq{k:02d}", "kind": "equilibration", "stage": stage, "number": k, "first": k == 1})
            if self.hyd and k == self.hyd["after_stage"]:
                n = 0
                while True:
                    rec = self.record_of(steps[-1]["tag"])
                    if rec is None:
                        break
                    check = rec.get("hydration_check")
                    if check is None:
                        raise DriverError(f"rec_{steps[-1]['tag']}.json holds no hydration check although the protocol asks for "
                                          "one: the step was made by an earlier run_md.py. Start this run again in an empty folder",
                                          EXIT_STATE)
                    if check["decision"] != "extend":
                        break
                    n += 1
                    if n > self.hyd["max_extensions"]:
                        raise DriverError(f"rec_{steps[-1]['tag']}.json asks for more extension steps than the protocol allows; "
                                          f"the folder {self.out} is damaged", EXIT_STATE)
                    name = f"{stage['name']}, extension {n} of at most {self.hyd['max_extensions']}"
                    steps.append({"tag": f"eq{k:02d}x{n}", "kind": "equilibration", "number": k, "first": False, "extension": n,
                                  "stage": dict(self.hyd["extension_stage"], name=name)})
        for n in range(1, n_segments + 1):
            steps.append({"tag": f"seg{n:04d}", "kind": "production", "number": n})
        return steps

    def stray_extension_steps(self, steps):
        """Files of extension steps that are not in the chain: none can exist in a sound folder."""
        known, stray = {s["tag"] for s in steps}, set()
        for fp in glob.glob(self.path("*_eq[0-9]*x[0-9]*.*")):
            tag = os.path.basename(fp).split("_", 1)[1].split(".", 1)[0]
            if tag not in known:
                stray.add(tag)
        for tag in sorted(stray):
            if self.record_of(tag) is not None:
                raise DriverError(f"step {tag} is recorded as complete but the records before it do not ask for it; "
                                  f"the folder {self.out} is damaged", EXIT_STATE)
            self.clean_step(tag)

    def record_of(self, tag):
        p = self.path(f"rec_{tag}.json")
        if not os.path.isfile(p):
            return None
        try:
            with open(p) as f:
                rec = json.load(f)
        except Exception:
            return None
        if not rec.get("complete"):
            return None
        for name, info in rec.get("files", {}).items():
            fp = self.path(name)
            if not os.path.isfile(fp):
                self.absent.add(name)           # allowed: finished files may have been moved to the archive
            elif os.path.getsize(fp) != info["bytes"]:
                raise DriverError(f"step {tag} is recorded as complete but {name} has another size than recorded; "
                                  f"the folder {self.out} is damaged", EXIT_STATE)
        return rec

    def files_of(self, tag):
        return sorted(glob.glob(self.path(f"*_{tag}.*")))

    def clean_step(self, tag):
        for fp in self.files_of(tag):
            os.remove(fp)
            self.say(f"  removed the incomplete file {os.path.basename(fp)}")

    def finish_step(self, tag, rec, files):
        rec["files"] = {name: file_info(self.path(name)) for name in files}
        rec["complete"] = True
        write_json(self.path(f"rec_{tag}.json"), rec)

    # ---- manifest
    def new_manifest(self, run_seed, seed_source, segment_ns):
        a = self.args
        test_only = a.equil_scale != 1.0 or self.proto.get("purpose") != "production"
        return {
            "schema": MANIFEST_SCHEMA, "written_by": PROGRAM,
            "run_id": a.run_id,
            "test_only": test_only,
            "test_only_reason": ([f"--equil-scale {a.equil_scale}"] if a.equil_scale != 1.0 else [])
                                + ([f"protocol purpose is '{self.proto.get('purpose')}'"] if self.proto.get("purpose") != "production" else []),
            "equil_scale": a.equil_scale,
            "seed": {"run_seed": run_seed, "source": seed_source,
                     "rule_run_seed": "first 4 bytes of sha256('cqc-sfa-md run seed|<run-id>') mod 2147483646, plus 1",
                     "rule_step_seed": "first 4 bytes of sha256('cqc-sfa-md step seed|<run seed>|<step tag>|<purpose>') mod 2147483646, plus 1"},
            "system_dir": self.sysdir, "out_dir": self.out,
            "inputs": self.inputs, "script": self.script_info,
            "build_json": self.build_info,
            "protocol": self.proto,
            "system": self.system_summary,
            "extra_forces": {"energy": "U = 0.5 k (r - r0)^2, minimum image", "harmonic_bonds": self.extra_bonds},
            "segment_ns": segment_ns,
            "created": stamp(),
            "invocations": [],
        }

    def check_manifest(self, man, run_seed_arg, segment_ns):
        a = self.args
        problems = []
        if man.get("run_id") != a.run_id:
            problems.append(f"run id is {man.get('run_id')}, not {a.run_id}")
        if run_seed_arg is not None and man["seed"]["run_seed"] != run_seed_arg:
            problems.append(f"run seed is {man['seed']['run_seed']}, not {run_seed_arg}")
        if man.get("equil_scale") != a.equil_scale:
            problems.append(f"--equil-scale was {man.get('equil_scale')}, now {a.equil_scale}")
        if abs(man.get("segment_ns", 0) - segment_ns) > 1e-12:
            problems.append(f"--segment-ns was {man.get('segment_ns')}, now {segment_ns}")
        for name, info in self.inputs.items():
            old = man["inputs"].get(name)
            if (old or {}).get("sha256") != (info or {}).get("sha256"):
                problems.append(f"{name} has changed since the run began")
        if problems:
            raise DriverError("the folder " + self.out + " holds another run: " + "; ".join(problems), EXIT_STATE)
        if man["script"]["sha256"] != self.script_info["sha256"]:
            self.say("WARNING: run_md.py has changed since this run began; both checksums are in the manifest.")

    def write_manifest(self, status):
        man = self.manifest
        steps = self.plan(self.n_segments_target)
        recs = {s["tag"]: self.record_of(s["tag"]) for s in steps}
        man["status"] = status
        man["minimisation"] = recs["min"]
        eq_steps = [s for s in steps if s["kind"] == "equilibration"]
        eq = [recs[s["tag"]] for s in eq_steps]
        eq_done = [r for r in eq if r is not None]
        eq_wall = sum(r["wall_s"] for r in eq_done)
        eq_ps = sum(r["settings"]["length_ps"] for r in eq_done)
        man["equilibration"] = {"complete": len(eq_done) == len(eq),
                                "n_stages": sum(1 for s in eq_steps if not s.get("extension")),
                                "n_stages_complete": sum(1 for s, r in zip(eq_steps, eq) if r is not None and not s.get("extension")),
                                "n_extension_steps": sum(1 for s in eq_steps if s.get("extension")),
                                "n_extension_steps_complete": sum(1 for s, r in zip(eq_steps, eq) if r is not None and s.get("extension")),
                                "completed_ps": round(eq_ps, 6), "wall_s": round(eq_wall, 1),
                                "ns_per_day": round(eq_ps * 1e-3 / eq_wall * 86400.0, 2) if eq_wall > 0 else None,
                                "stages_note": "the records of all complete equilibration steps in the order of the chain; "
                                               "the extension steps of the hydration check (eqNNx<n>) are among them",
                                "stages": eq_done}
        man["hydration_check"] = self.hydration_summary(eq_steps, recs)
        segs = [recs[s["tag"]] for s in steps if s["kind"] == "production"]
        done = [r for r in segs if r is not None]
        wall = sum(r["wall_s"] for r in done)
        man["production"] = {
            "target_ns": self.prod_ns, "segment_ns": man["segment_ns"], "n_segments_target": len(segs),
            "n_segments_complete": len(done),
            "completed_ns": round(sum(r["settings"]["length_ps"] for r in done) / 1000.0, 9),
            "complete": len(done) == len(segs) and len(segs) > 0,
            "frames_solute": sum(r["frames"]["solute"] for r in done),
            "frames_full": sum(r["frames"]["full"] for r in done),
            "wall_s": round(wall, 1),
            "ns_per_day": (round(sum(r["settings"]["length_ps"] for r in done) * 1e-3 / wall * 86400.0, 2) if wall > 0 else None),
            "segments": done,
        }
        man["all_steps_complete"] = all(r is not None for r in recs.values())
        man["last_written"] = stamp()
        write_json(self.path("manifest.json"), man)

    def hydration_summary(self, eq_steps, recs):
        """The hydration check of the run for the manifest, gathered from the records of the steps."""
        if not self.hyd:
            if self.protocol_kept:
                return dict(self.protocol_kept, requested=False,
                            status="not made: the run began with a protocol without the block equilibration.hydration_check")
            return {"requested": False, "status": "not made: the protocol has no block equilibration.hydration_check"}
        h = self.hyd
        dt_fs = float(h["extension_stage"]["timestep_fs"])
        ext_steps = max(1, int(round(h["extension_ps"] * 1000.0 / dt_fs * self.args.equil_scale)))
        done = [(s, recs[s["tag"]]) for s in eq_steps if s["number"] == h["after_stage"] and recs[s["tag"]] is not None]
        counts = [{"step": s["tag"], "time_end_ps": r["time_end_ps"], "n_water_O": r["hydration_check"]["n_water_O"],
                   "nearest_water_O_A": r["hydration_check"]["nearest_water_O"]["distance_A"],
                   "ligand_mean_distance_from_site_A": r["hydration_check"]["ligand_mean_distance_from_site_A"],
                   "applies": r["hydration_check"]["applies"], "decision": r["hydration_check"]["decision"]} for s, r in done]
        last = done[-1][1]["hydration_check"] if done else None
        decided = last is not None and last["decision"] != "extend"
        ext = [r for s, r in done if s.get("extension")]
        return {
            "requested": True,
            "rule": "md/PLAN.md, entries of 6 Oct 2026, 1:03 AM (point 3) and 3:02 AM (point 8); settings in the protocol, "
                    "key equilibration.hydration_check",
            "settings": {"after_step": h["base_tag"], "site_atoms": h["site_atoms"],
                         "water_cutoff_A": round(h["cutoff_nm"] * 10.0, 6),
                         "applies_if_ligand_mean_distance_above_A": round(h["apart_nm"] * 10.0, 6),
                         "max_extensions": h["max_extensions"], "extension_length_ps_of_the_protocol": h["extension_ps"],
                         "extension_length_ps_in_this_run": round(ext_steps * dt_fs / 1000.0, 6),
                         "extension_steps_in_this_run": ext_steps},
            "reference_guard": h["reference_guard"],
            "status": "decided" if decided else "pending: the held part of the equilibration is not complete",
            "applies": last["applies"] if last else None,
            "skipped": last["skipped"] if last else None,
            "counts": counts,
            "n_extensions": len(ext),
            "extension_ps": round(sum(r["settings"]["length_ps"] for r in ext), 6),
            "final_decision": f"{last['decision']}: {last['reason']}" if decided else None,
            "site_without_water_when_the_solute_was_released": (bool(last["applies"] and last["n_water_O"] == 0)
                                                                if decided else None),
        }

    # ---- the steps
    def write_topology_files(self):
        full = self.wrapper(self.pos0, self.box0)
        rows_all = [(i, self.atom_names[i], self.res_names[self.res_of_atom[i]], int(self.res_of_atom[i]), self.elements[i])
                    for i in range(self.n_atoms)]
        note = ["coordinates: those of system.inpcrd (Angstrom); serial = prmtop atom index + 1",
                f"written by {PROGRAM} for run {self.args.run_id}"]
        if not os.path.isfile(self.path("full.pdb")):
            write_atomic(self.path("full.pdb"), pdb_text(rows_all, full, self.box0,
                         ["topology of full_*.dcd: all atoms, in the order of system.prmtop"] + note))
        if not os.path.isfile(self.path("solute.pdb")):
            write_atomic(self.path("solute.pdb"), pdb_text([rows_all[i] for i in self.solute], self.pos0[self.solute], self.box0,
                         ["topology of solute_*.dcd: the atoms of selections.json 'solute', in that order"] + note))

    def run_min(self, step):
        m = self.proto["minimisation"]
        tag = step["tag"]
        rec = {"tag": tag, "kind": "minimisation", "started": stamp(), "host": socket.gethostname(),
               "openmm": openmm.Platform.getOpenMMVersion()}
        system = self.build_system(m["restraint_atoms"], m["restraint_k_kj_per_mol_nm2"], None)
        ctx, integ = self.make_context(system, 2.0, self.proto["integrator"]["temperature_K"], 1)
        self.set_state(ctx, self.pos0, None, self.box0)

        def energies():
            return {"potential_energy_kj_per_mol": ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                    "force_field_kj_per_mol": self.group_energy(ctx, GROUP_FF),
                    "extra_bonds_kj_per_mol": self.group_energy(ctx, GROUP_EXTRA),
                    "position_restraints_kj_per_mol": self.group_energy(ctx, GROUP_POSRES)}

        before = energies()
        self.say(f"min: potential energy before {before['potential_energy_kj_per_mol']:.1f} kJ/mol; minimising "
                 f"({len(self.atom_sets[m['restraint_atoms']])} atoms restrained, k = {m['restraint_k_kj_per_mol_nm2']} kJ mol-1 nm-2)")

        class Counter(openmm.MinimizationReporter):
            calls, runs, last = 0, [], None

            def report(self, iteration, x, grad, args):
                Counter.calls += 1
                if iteration == 0:
                    Counter.runs.append(0)          # the minimiser starts L-BFGS again with stiffer constraint terms
                if Counter.runs:
                    Counter.runs[-1] = iteration + 1
                Counter.last = {str(k): float(v) for k, v in dict(args).items()}
                return _STOP["signal"] is not None

        t0 = time.time()
        openmm.LocalEnergyMinimizer.minimize(ctx, m["tolerance_kj_per_mol_nm"], int(m["max_iterations"]), Counter())
        check_stop()
        wall = time.time() - t0
        after = energies()
        if not all(math.isfinite(v) for v in after.values()):
            raise DriverError("the energy after the minimisation is not finite", EXIT_SIM)
        st = ctx.getState(getPositions=True, enforcePeriodicBox=False)
        pos = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        ctx.setVelocities(np.zeros((self.n_atoms, 3)))
        self.save_state(ctx, f"state_{tag}.xml")
        moved = np.linalg.norm(pos - self.pos0, axis=1)
        restrained = self.atom_sets[m["restraint_atoms"]]
        rec.update({
            "platform": self.plat.record(),
            "settings": dict(m, n_restrained_atoms=len(restrained)),
            "energies_before": {k: round(v, 3) for k, v in before.items()},
            "energies_after": {k: round(v, 3) for k, v in after.items()},
            "energy_change_kj_per_mol": round(after["potential_energy_kj_per_mol"] - before["potential_energy_kj_per_mol"], 3),
            "minimiser_calls": Counter.calls,
            "lbfgs_iterations_per_pass": Counter.runs,
            "passes_stopped_by_the_iteration_limit": sum(1 for r in Counter.runs if r >= int(m["max_iterations"])),
            "last_pass_stopped_by_the_iteration_limit": bool(Counter.runs) and Counter.runs[-1] >= int(m["max_iterations"]),
            "minimiser_last_report": Counter.last,
            "restrained_atoms_displacement_A": {"rms": round(float(np.sqrt((moved[restrained] ** 2).mean())) * 10, 4),
                                                "max": round(float(moved[restrained].max()) * 10, 4)} if restrained else None,
            "all_atoms_max_displacement_A": round(float(moved.max()) * 10, 3),
            "monitor": self.monitor(pos, self.box0),
            "wall_s": round(wall, 2), "ended": stamp(),
        })
        del ctx, integ
        self.finish_step(tag, rec, [f"state_{tag}.xml"])
        self.say(f"min: potential energy after {after['potential_energy_kj_per_mol']:.1f} kJ/mol "
                 f"(change {rec['energy_change_kj_per_mol']:.1f}), {Counter.calls} calls, {wall:.0f} s; L-BFGS iterations per pass "
                 + ", ".join(str(r) for r in Counter.runs))
        if rec["last_pass_stopped_by_the_iteration_limit"]:
            self.say(f"WARNING: the last pass of the minimiser was stopped by the limit of {int(m['max_iterations'])} iterations "
                     "per pass; the tolerance was not reached. The run goes on (rec_min.json has the numbers).")

    def run_md_step(self, step, prev_tag):
        """One equilibration stage or one production segment."""
        tag, kind = step["tag"], step["kind"]
        out = self.proto["output"]
        seeds = {"integrator": step_seed(self.run_seed, tag, "integrator")}
        if kind == "equilibration":
            stage = step["stage"]
            dt_fs = float(stage["timestep_fs"])
            n_steps = max(1, int(round(stage["length_ps"] * 1000.0 / dt_fs * self.args.equil_scale)))
            ensemble = stage["ensemble"]
            temps = stage["temperature_K"]
            ramp = list(temps) if isinstance(temps, list) else None
            t_start = ramp[0] if ramp else float(temps)
            rset, rk = stage["restraint_atoms"], float(stage["restraint_k_kj_per_mol_nm2"])
            first = bool(step.get("first"))
            time0_ps, step0 = self.equil_time_ps, self.equil_steps
            label = f"{tag} ({stage['name']})"
        else:
            prod = self.proto["production"]
            dt_fs = float(prod["timestep_fs"])
            n_steps = self.steps_per_segment
            ensemble, ramp, t_start = prod["ensemble"], None, float(prod["temperature_K"])
            rset, rk, first = "none", 0.0, False
            step0 = (step["number"] - 1) * n_steps
            time0_ps = step0 * dt_fs / 1000.0
            label = f"{tag} ({time0_ps / 1000.0:g}-{(time0_ps + n_steps * dt_fs / 1000.0) / 1000.0:g} ns)"
        if ensemble not in ("NVT", "NPT"):
            raise DriverError(f"{tag}: ensemble must be NVT or NPT", EXIT_INPUT)
        if ensemble == "NPT":
            seeds["barostat"] = step_seed(self.run_seed, tag, "barostat")
        if first:
            seeds["velocities"] = step_seed(self.run_seed, tag, "velocities")
        every_log = steps_for(out["state_log_interval_ps"], dt_fs, "state log interval")
        every_sol = steps_for(out["solute_frame_interval_ps"], dt_fs, "solute frame interval")
        every_full = steps_for(out["full_frame_interval_ps"], dt_fs, "full frame interval")
        ramp_block = steps_for(1.0, dt_fs, "temperature ramp block") if ramp else None
        n_blocks = int(math.ceil(n_steps / ramp_block)) if ramp else 0

        state_in = f"state_{prev_tag}.xml"
        rec = {"tag": tag, "kind": kind, "started": stamp(), "host": socket.gethostname(),
               "openmm": openmm.Platform.getOpenMMVersion(),
               "state_in": dict(file_info(self.path(state_in)), file=state_in), "seeds": seeds}
        if self.state_in_expected.get(prev_tag) not in (None, rec["state_in"]["sha256"]):
            raise DriverError(f"{state_in} is not the file recorded at the end of step {prev_tag}", EXIT_STATE)
        pos, vel, box = self.read_state(state_in)

        system = self.build_system(rset, rk, seeds.get("barostat"))
        ctx, integ = self.make_context(system, dt_fs, t_start, seeds["integrator"])
        self.set_state(ctx, pos, None if first else vel, box)
        if first:
            ctx.applyConstraints(self.proto["system"]["constraint_tolerance"])
            ctx.setVelocitiesToTemperature(t_start * unit.kelvin, seeds["velocities"])
        ctx.setTime(time0_ps * unit.picosecond)
        ctx.setStepCount(step0)
        has_posres = bool(self.atom_sets[rset]) and rk > 0
        self.say(f"{label}: {ensemble}, {n_steps} steps of {dt_fs:g} fs, restraints {rset} k = {rk:g}, on {self.plat.name}")

        log_name, sol_name, full_name = f"log_{tag}.csv", f"solute_{tag}.dcd", f"full_{tag}.dcd"
        titles = [f"cqc-sfa-md run {self.args.run_id} step {tag}; coordinates in Angstrom"[:80],
                  f"time zero: start of the {'equilibration' if kind == 'equilibration' else 'production'}; step {dt_fs:g} fs"]
        # the header counts steps of this step's length since time zero, so that readers show the right time
        head0 = int(round(time0_ps * 1000.0 / dt_fs))
        dcd_sol = DcdWriter(self.path(sol_name), len(self.solute), head0 + every_sol, every_sol, dt_fs / 1000.0, titles)
        dcd_full = DcdWriter(self.path(full_name), self.n_atoms, head0 + every_full, every_full, dt_fs / 1000.0, titles)
        columns = ["step", "time_ps", "potential_energy_kj_per_mol", "kinetic_energy_kj_per_mol", "temperature_K",
                   "volume_nm3", "density_g_per_cm3", "position_restraint_energy_kj_per_mol",
                   "extra_bond_energy_kj_per_mol", "thermostat_temperature_K", "speed_ns_per_day"]
        rows = []
        t_wall = time.time()
        target_T = t_start
        try:
            with open(self.path(log_name), "w") as log:
                log.write(",".join(columns) + "\n")
                done = 0
                while done < n_steps:
                    if ramp:
                        block = done // ramp_block
                        target_T = ramp[0] + (ramp[1] - ramp[0]) * (block + 1) / n_blocks
                        integ.setTemperature(target_T * unit.kelvin)
                    nxt = min([n_steps] + [(done // m + 1) * m for m in (every_log, every_sol, every_full)]
                              + ([(done // ramp_block + 1) * ramp_block] if ramp else []))
                    while done < nxt:
                        n = min(CHUNK_STEPS, nxt - done)
                        integ.step(n)
                        done += n
                        check_stop()
                    want_log, want_sol, want_full = done % every_log == 0, done % every_sol == 0, done % every_full == 0
                    if not (want_log or want_sol or want_full):
                        continue
                    st = ctx.getState(getPositions=want_sol or want_full, getEnergy=want_log, enforcePeriodicBox=False)
                    fbox = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
                    if want_log:
                        pe = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                        ke = st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                        if not (math.isfinite(pe) and math.isfinite(ke)):
                            raise DriverError(f"{tag}: the energy is not finite at step {step0 + done}; the simulation has failed", EXIT_SIM)
                        volume = float(abs(np.linalg.det(fbox)))
                        elapsed = max(time.time() - t_wall, 1e-9)
                        row = [step0 + done, time0_ps + done * dt_fs / 1000.0, pe, ke,
                               2.0 * ke / (self.dof * GAS_CONSTANT_KJ_PER_MOL_K), volume,
                               self.total_mass_amu * AMU_PER_NM3_TO_G_PER_CM3 / volume,
                               self.group_energy(ctx, GROUP_POSRES) if has_posres else 0.0,
                               self.group_energy(ctx, GROUP_EXTRA) if self.extra_bonds else 0.0,
                               target_T, done * dt_fs * 1e-6 / elapsed * 86400.0]
                        rows.append(row)
                        log.write("%d,%.3f,%.4f,%.4f,%.3f,%.5f,%.6f,%.4f,%.4f,%.2f,%.2f\n" % tuple(row))
                        log.flush()
                    if want_sol or want_full:
                        fpos = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
                        if want_sol:
                            dcd_sol.write(fpos[self.solute], fbox)
                        if want_full:
                            dcd_full.write(self.wrapper(fpos, fbox), fbox)
                log.flush()
                os.fsync(log.fileno())
        except openmm.OpenMMException as exc:
            raise DriverError(f"{tag}: OpenMM stopped the simulation ({exc})", EXIT_SIM)
        finally:
            dcd_sol.close()
            dcd_full.close()
        st = ctx.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=False)
        pe = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(pe):
            raise DriverError(f"{tag}: the energy is not finite at the end of the step", EXIT_SIM)
        wall = time.time() - t_wall
        write_atomic(self.path(f"state_{tag}.xml"), openmm.XmlSerializer.serialize(st))
        files = [log_name, sol_name, full_name, f"state_{tag}.xml"]
        fpos = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        fbox = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
        volume = float(abs(np.linalg.det(fbox)))
        arr = np.array(rows) if rows else np.zeros((0, len(columns)))
        length_ps = n_steps * dt_fs / 1000.0
        rec.update({
            "platform": self.plat.record(),
            "settings": {"name": stage["name"] if kind == "equilibration" else "production",
                         "ensemble": ensemble, "timestep_fs": dt_fs, "n_steps": n_steps, "length_ps": round(length_ps, 6),
                         "temperature_K": ramp if ramp else t_start,
                         "restraint_atoms": rset, "restraint_k_kj_per_mol_nm2": rk,
                         "n_restrained_atoms": len(self.atom_sets[rset]) if rk > 0 else 0,
                         "barostat": ensemble == "NPT"},
            "first_step": step0, "last_step": step0 + n_steps,
            "time_start_ps": round(time0_ps, 6), "time_end_ps": round(time0_ps + length_ps, 6),
            "frames": {"solute": dcd_sol.count, "full": dcd_full.count, "state_log_rows": len(rows)},
            "final": {"potential_energy_kj_per_mol": round(pe, 3),
                      "kinetic_energy_kj_per_mol": round(st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole), 3),
                      "temperature_K": round(2.0 * st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                                             / (self.dof * GAS_CONSTANT_KJ_PER_MOL_K), 2),
                      "volume_nm3": round(volume, 4),
                      "density_g_per_cm3": round(self.total_mass_amu * AMU_PER_NM3_TO_G_PER_CM3 / volume, 5)},
            "mean_over_state_log": ({"temperature_K": round(float(arr[:, 4].mean()), 2),
                                     "density_g_per_cm3": round(float(arr[:, 6].mean()), 5),
                                     "potential_energy_kj_per_mol": round(float(arr[:, 2].mean()), 2)} if len(rows) else None),
            "monitor": self.monitor(fpos, fbox),
            "wall_s": round(wall, 2),
            "ns_per_day": round(length_ps * 1e-3 / max(wall, 1e-9) * 86400.0, 2),
            "ended": stamp(),
        })
        check = None
        if kind == "equilibration" and self.hyd and step["number"] == self.hyd["after_stage"]:
            if step.get("extension"):
                rec["settings"].update(extension_of=self.hyd["base_tag"], extension_number=step["extension"])
            check = rec["hydration_check"] = self.hydration_record(fpos, fbox, tag, step.get("extension", 0))
        if kind == "production":
            chk = f"checkpoint_{tag}.chk"
            write_atomic(self.path(chk), ctx.createCheckpoint())
            rec["checkpoint"] = dict(file_info(self.path(chk)), file=chk,
                                     note="binary, tied to this machine and platform; --resume uses the XML state")
        del ctx, integ
        self.finish_step(tag, rec, files)
        if kind == "production" and out.get("keep_checkpoints", "last") == "last" and step["number"] > 1:
            old = self.path(f"checkpoint_seg{step['number'] - 1:04d}.chk")
            if os.path.isfile(old):
                os.remove(old)
        self.state_in_expected[tag] = rec["files"][f"state_{tag}.xml"]["sha256"]
        if kind == "equilibration":
            self.equil_time_ps += length_ps
            self.equil_steps += n_steps
        self.say(f"{label}: done in {wall:.0f} s, {rec['ns_per_day']:.1f} ns/day; T {rec['final']['temperature_K']:.1f} K, "
                 f"density {rec['final']['density_g_per_cm3']:.4f} g/cm3")
        if check:
            self.say(f"{tag}: hydration check: {check['n_water_O']} water oxygen(s) within {check['water_cutoff_A']:g} A of the "
                     f"crystal ring site (nearest {check['nearest_water_O']['distance_A']:.2f} A); ligand "
                     f"{check['ligand_mean_distance_from_site_A']:.2f} A from the site on average; {check['reason']}; "
                     f"next step {check['next_step']}")

    def run_benchmark(self):
        a, prod = self.args, self.proto["production"]
        dt_fs = float(prod["timestep_fs"])
        eq_tags = [s["tag"] for s in self.plan(0) if s["kind"] == "equilibration"]
        if self.record_of(eq_tags[-1]):
            source = f"state_{eq_tags[-1]}.xml"
        else:
            if not self.record_of("min"):
                self.clean_step("min")
                self.run_min({"tag": "min"})
            source = "state_min.xml"
        pos, vel, box = self.read_state(source)
        equilibrated = source != "state_min.xml"
        seeds = {p: step_seed(self.run_seed, "benchmark", p) for p in ("integrator", "barostat", "velocities")}
        system = self.build_system("none", 0.0, seeds["barostat"])
        ctx, integ = self.make_context(system, dt_fs, float(prod["temperature_K"]), seeds["integrator"])
        self.set_state(ctx, pos, vel if equilibrated else None, box)
        warm = 0
        if not equilibrated:
            ctx.applyConstraints(self.proto["system"]["constraint_tolerance"])
            ctx.setVelocitiesToTemperature(float(prod["temperature_K"]) * unit.kelvin, seeds["velocities"])
            integ.setStepSize(1.0 * unit.femtosecond)
            warm = 1000
            integ.step(warm)                      # 1 ps at 1 fs, not timed: the structure is only minimised
            integ.setStepSize(dt_fs * unit.femtosecond)
        self.say(f"benchmark: {self.n_atoms} atoms from {source}, {dt_fs:g} fs steps, NPT, on {self.plat.name}; "
                 f"timing for {a.benchmark_seconds:g} s")
        try:
            integ.step(10)
            ctx.getState(getEnergy=True)
            chunk, steps, t0 = 10, 0, time.time()
            while True:
                integ.step(chunk)
                pe = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                steps += chunk
                elapsed = time.time() - t0
                check_stop()
                if not math.isfinite(pe):
                    raise DriverError("benchmark: the energy is not finite"
                                      + ("" if equilibrated else "; run the equilibration first"), EXIT_SIM)
                if elapsed >= a.benchmark_seconds:
                    break
                rate = steps / max(elapsed, 1e-9)
                chunk = int(max(10, min(rate * 2.0, rate * (a.benchmark_seconds - elapsed) + 1)))
        except openmm.OpenMMException as exc:
            raise DriverError(f"benchmark: OpenMM stopped the simulation ({exc})", EXIT_SIM)
        ns_day = steps * dt_fs * 1e-6 / elapsed * 86400.0
        result = {"run_id": a.run_id, "ns_per_day": round(ns_day, 2), "steps": steps, "wall_s": round(elapsed, 2),
                  "timestep_fs": dt_fs, "ensemble": "NPT", "n_atoms": self.n_atoms, "start_state": source,
                  "warm_up_steps_of_1_fs_not_timed": warm, "trajectory_output": "none (integration only)",
                  "hours_for_100_ns": round(100.0 / ns_day * 24.0, 2),
                  "platform": self.plat.record(), "host": host_info(), "when": stamp()}
        bpath = self.path("benchmark.json")
        allb = read_json(bpath, "benchmark.json") if os.path.isfile(bpath) else []
        allb.append(result)
        write_json(bpath, allb)
        del ctx, integ
        self.say(f"BENCHMARK {a.run_id}: {ns_day:.2f} ns/day ({self.plat.name}"
                 + (f", {self.plat.gpu_name}" if self.plat.gpu_name else "")
                 + f", {self.n_atoms} atoms, {dt_fs:g} fs, {steps} steps in {elapsed:.1f} s); 100 ns would take {result['hours_for_100_ns']:.1f} h")

    # ---- main sequence
    def main(self):
        a = self.args
        self.load_inputs()
        prod = self.proto["production"]
        self.prod_ns = float(a.prod_ns if a.prod_ns is not None else prod["length_ns"])
        segment_ns = float(a.segment_ns)
        dt_fs = float(prod["timestep_fs"])
        self.steps_per_segment = steps_for(segment_ns * 1000.0, dt_fs, "--segment-ns")
        out = self.proto["output"]
        for key in ("state_log_interval_ps", "solute_frame_interval_ps", "full_frame_interval_ps"):
            if self.steps_per_segment % steps_for(out[key], dt_fs, key) != 0:
                raise DriverError(f"--segment-ns {segment_ns:g} is not a whole multiple of {key} ({out[key]} ps)", EXIT_INPUT)
        n_seg = self.prod_ns / segment_ns
        if abs(n_seg - round(n_seg)) > 1e-6 or round(n_seg) < 1:
            raise DriverError(f"--prod-ns {self.prod_ns:g} is not a whole multiple of --segment-ns {segment_ns:g}", EXIT_INPUT)
        self.n_segments_target = int(round(n_seg))

        os.makedirs(self.out, exist_ok=True)
        man_path = self.path("manifest.json")
        have_manifest = os.path.isfile(man_path)
        have_records = bool(glob.glob(self.path("rec_*.json")))
        if (have_manifest or have_records) and not a.resume and a.benchmark_seconds is None:
            raise DriverError(f"{self.out} already holds a run; add --resume to continue it or choose another --out", EXIT_STATE)
        if have_records and not have_manifest:
            raise DriverError(f"{self.out} holds step records but no manifest.json; the folder is damaged", EXIT_STATE)
        self.take_lock()
        self.log_handle = open(self.path("run_md.log"), "a")
        self.say(f"{PROGRAM}: run {a.run_id}, stage {a.stage}, system {self.sysdir}, out {self.out}"
                 + (", --resume" if a.resume else ""))

        self.plat = PlatformChoice(a.platform, a.threads, self.proto, self.say)
        self.system_summary = self.describe_system()
        if have_manifest:
            self.manifest = read_json(man_path, "manifest.json", EXIT_STATE)
            self.check_manifest(self.manifest, a.seed, segment_ns)
            self.run_seed = self.manifest["seed"]["run_seed"]
            for inv in self.manifest["invocations"]:
                if "ended" not in inv:
                    inv["ended"] = None
                    inv["outcome"] = "no end recorded (killed or crashed)"
        else:
            self.run_seed = a.seed if a.seed is not None else run_seed_from_id(a.run_id)
            if not 1 <= self.run_seed <= 2147483646:
                raise DriverError("--seed must be between 1 and 2147483646", EXIT_INPUT)
            self.manifest = self.new_manifest(self.run_seed, "--seed" if a.seed is not None else "derived from the run id", segment_ns)
        if self.manifest["test_only"]:
            write_atomic(self.path("TEST_ONLY"), "This folder holds a test run, not a production run: "
                         + "; ".join(self.manifest["test_only_reason"]) + "\n")
            self.say("TEST RUN (not a production run): " + "; ".join(self.manifest["test_only_reason"]))
        self.invocation = {"started": stamp(), "argv": sys.argv, "pid": os.getpid(), "stage": a.stage,
                           "resume": bool(a.resume), "prod_ns": self.prod_ns,
                           "benchmark_seconds": a.benchmark_seconds, "host": host_info(),
                           "platform": self.plat.record(), "script_sha256": self.script_info["sha256"]}
        if self.protocol_kept:
            self.invocation["protocol_kept"] = self.protocol_kept
            self.say("NOTE: the protocol file has gained the block equilibration.hydration_check since this run began and is "
                     "otherwise the same; the run keeps the protocol of its manifest (no hydration check).")
        self.manifest["invocations"].append(self.invocation)
        self.write_topology_files()
        self.write_manifest("running")

        if a.benchmark_seconds is not None:
            self.run_benchmark()
            self.invocation["platform"] = self.plat.record()
            return "benchmark done"

        steps = self.plan(self.n_segments_target)
        self.stray_extension_steps(steps)
        wanted = {"all": ("minimisation", "equilibration", "production"), "min": ("minimisation",),
                  "equil": ("equilibration",), "prod": ("production",)}[a.stage]
        # where the chain stands
        self.state_in_expected, self.equil_time_ps, self.equil_steps = {}, 0.0, 0
        first_open = None
        for n, s in enumerate(steps):
            rec = self.record_of(s["tag"])
            if rec is None:
                if first_open is None:
                    first_open = n
            else:
                if first_open is not None:
                    raise DriverError(f"step {s['tag']} is complete but the earlier step {steps[first_open]['tag']} is not; "
                                      f"the folder {self.out} is damaged", EXIT_STATE)
                self.state_in_expected[s["tag"]] = rec["files"][f"state_{s['tag']}.xml"]["sha256"]
                if s["kind"] == "equilibration":
                    self.equil_time_ps += rec["settings"]["length_ps"]
                    self.equil_steps += rec["settings"]["n_steps"]
        if self.absent:
            some = sorted(self.absent)
            self.say(f"NOTE: {len(some)} file(s) of complete steps are not in this folder (moved away?), for example {some[0]}; "
                     "they are not needed to continue")
        if first_open is None:
            self.say("every step of the request is already complete")
        else:
            if first_open > 0:
                prev = steps[first_open - 1]["tag"]
                if not os.path.isfile(self.path(f"state_{prev}.xml")):
                    raise DriverError(f"state_{prev}.xml is needed to continue and is not in {self.out}", EXIT_STATE)
                if sha256_file(self.path(f"state_{prev}.xml")) != self.state_in_expected[prev]:
                    raise DriverError(f"state_{prev}.xml does not have the recorded checksum; the folder {self.out} is damaged", EXIT_STATE)
                self.say(f"steps complete up to {prev}; continuing with {steps[first_open]['tag']}")
            for s in steps[first_open:]:
                if self.files_of(s["tag"]):
                    self.clean_step(s["tag"])
            n = first_open
            while n < len(steps):
                s = steps[n]
                if s["kind"] not in wanted:
                    if any(t["kind"] in wanted for t in steps[n + 1:]):
                        raise DriverError(f"--stage {a.stage} needs step {s['tag']} first; run --stage all, or the earlier stage", EXIT_STATE)
                    break
                check_stop()
                if s["kind"] == "minimisation":
                    self.run_min(s)
                else:
                    self.run_md_step(s, steps[n - 1]["tag"])
                self.invocation["platform"] = self.plat.record()
                steps = self.plan(self.n_segments_target)      # the new record may put an extension step next
                if steps[n]["tag"] != s["tag"]:
                    raise DriverError(f"the chain of steps changed under step {s['tag']}; the folder {self.out} is damaged", EXIT_STATE)
                self.write_manifest("running")
                n += 1
        last_seg = self.record_of(steps[-1]["tag"])
        if last_seg is not None and os.path.isfile(self.path(f"state_{steps[-1]['tag']}.xml")):
            pos, _vel, box = self.read_state(f"state_{steps[-1]['tag']}.xml")
            rows = [(int(i), self.atom_names[i], self.res_names[self.res_of_atom[i]], int(self.res_of_atom[i]), self.elements[i])
                    for i in self.solute]
            write_atomic(self.path("final_solute.pdb"), pdb_text(rows, pos[self.solute], box, [
                f"solute at the end of step {steps[-1]['tag']} of run {self.args.run_id} ({self.prod_ns:g} ns of production)",
                "coordinates in Angstrom, as integrated; serial = prmtop atom index + 1"]))
        p = self.manifest_production_line()
        if p:
            self.say(p)
        return "stage complete"

    def manifest_production_line(self):
        done = [r for r in (self.record_of(f"seg{n:04d}") for n in range(1, self.n_segments_target + 1)) if r]
        if not done:
            return None
        wall = sum(r["wall_s"] for r in done)
        ns = sum(r["settings"]["length_ps"] for r in done) / 1000.0
        return (f"production: {len(done)} of {self.n_segments_target} segments complete, {ns:g} of {self.prod_ns:g} ns, "
                f"{ns / max(wall, 1e-9) * 86400.0:.1f} ns/day")


def parse_args(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="One MD trajectory of a built system, in steps that survive interruptions "
                                            "(see md/RUN.md).")
    p.add_argument("--system", required=True, help="folder of a built system (system.prmtop, system.inpcrd, selections.json, extra_forces.json)")
    p.add_argument("--out", required=True, help="folder of the run (created if absent)")
    p.add_argument("--run-id", required=True, help="name of the run, for example X1")
    p.add_argument("--seed", type=int, default=None, help="run seed, 1..2147483646; default: derived from the run id")
    p.add_argument("--platform", default="auto", choices=["auto", "CUDA", "OpenCL", "CPU"])
    p.add_argument("--threads", type=int, default=None, help="threads of the CPU platform (default: OpenMM's choice)")
    p.add_argument("--stage", default="all", choices=["all", "min", "equil", "prod"])
    p.add_argument("--prod-ns", type=float, default=None, help="length of the production (default: the protocol's, 100)")
    p.add_argument("--segment-ns", type=float, default=5.0, help="length of one production segment (default 5)")
    p.add_argument("--resume", action="store_true", help="continue the run in --out (starts it if the folder is empty)")
    p.add_argument("--benchmark-seconds", type=float, default=None, help="time the equilibrated (or minimised) system for this wall time and print ns per day")
    p.add_argument("--equil-scale", type=float, default=1.0, help="tests only: shorten every equilibration stage and extension step by this factor (recorded as a test run)")
    p.add_argument("--protocol", default=os.path.join(here, "protocol.json"), help="protocol file (default: protocol.json beside this script)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    version = tuple(int(x) for x in openmm.Platform.getOpenMMVersion().split(".")[:2])
    if version < (8, 1):
        sys.stderr.write(f"ERROR: OpenMM {openmm.Platform.getOpenMMVersion()} is too old; 8.1 or newer is needed.\n")
        return EXIT_INPUT
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)
    run = Run(args)
    code, outcome = EXIT_OK, "done"
    try:
        outcome = run.main()
    except Stopped:
        signum = _STOP["signal"]
        code, outcome = 128 + signum, f"stopped by signal {signum}; continue with --resume"
    except DriverError as exc:
        code, outcome = exc.code, str(exc)
    except Exception as exc:
        import traceback
        traceback.print_exc()
        code, outcome = EXIT_ERROR, f"unexpected error: {type(exc).__name__}: {exc}"
    if run.invocation is not None:
        run.invocation["ended"] = stamp()
        run.invocation["exit_code"] = code
        run.invocation["outcome"] = outcome
        try:
            if code == EXIT_OK:
                status = "complete" if run.record_of(run.plan(run.n_segments_target)[-1]["tag"]) else "incomplete"
            else:
                status = "stopped" if code > 128 else "failed"
            run.write_manifest(status)
        except Exception as exc:
            sys.stderr.write(f"ERROR: the manifest could not be written: {exc}\n")
            code = code or EXIT_ERROR
    run.release_lock()
    if code != EXIT_OK:
        text = f"ERROR ({code}): {outcome}"
        sys.stderr.write(text + "\n")
        sys.stderr.flush()
        if run.log_handle:
            run.log_handle.write(f"[{stamp()['ist']}] {text}\n")
            run.log_handle.close()
        return code
    if run.log_handle:
        run.log_handle.write(f"[{stamp()['ist']}] DONE {args.run_id}\n")
        run.log_handle.close()
    print(f"DONE {args.run_id}", flush=True)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
