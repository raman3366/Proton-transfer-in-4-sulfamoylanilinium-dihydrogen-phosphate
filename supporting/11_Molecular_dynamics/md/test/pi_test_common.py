"""pi_test_common.py: readers and geometry shared by the scripts of the pipeline test on the Raspberry Pi.

TEST OF THE PIPELINE with the restrained-zinc model. NOT a result for the paper (md/test/PI_TEST.md).

Written independently of md/run_md.py and md/analysis/analyse_md.py (nothing is imported from them), so
that the checks do not repeat a mistake of the code they check. Needs NumPy; OpenMM only to read a saved
state (state_<step>.xml).

Units: coordinates and distances in Angstrom everywhere in this module (DCD and inpcrd files are in
Angstrom; a saved state is in nm and is converted on reading); angles in degrees; times in ps.
"""
import datetime
import hashlib
import itertools
import json
import os
import struct

import numpy as np

IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), "IST")
NOTE = ("TEST OF THE PIPELINE with the restrained-zinc model (nonbonded Zn2+, four Zn-N distance restraints); "
        "NOT a result for the paper")
DCD_HEADER_BYTES = 276


def now_ist():
    t = datetime.datetime.now(IST)
    hour = t.hour % 12 or 12
    return f"{t.day} {t.strftime('%b %Y')}, {hour}:{t.minute:02d}:{t.second:02d} {'AM' if t.hour < 12 else 'PM'} IST"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    with open(path) as f:
        return json.load(f)


def write_json(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")


def software():
    import platform
    out = {"python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine(),
           "host": platform.node()}
    try:
        import openmm
        out["openmm"] = openmm.Platform.getOpenMMVersion()
    except Exception as exc:                                  # pragma: no cover
        out["openmm"] = f"not importable ({exc})"
    try:
        import MDAnalysis
        out["MDAnalysis"] = MDAnalysis.__version__
    except Exception as exc:                                  # pragma: no cover
        out["MDAnalysis"] = f"not importable ({exc})"
    return out


# ----------------------------------------------------------------------------- files

def dcd_frame_bytes(n_atoms):
    """Bytes of one frame of a DCD file with a unit cell: the cell record and three coordinate records."""
    return (4 + 48 + 4) + 3 * (4 + 4 * n_atoms + 4)


def read_dcd_header(path):
    """Header fields of a DCD file in the CHARMM layout, and the number of whole frames by the file size."""
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        head = f.read(DCD_HEADER_BYTES)
    if len(head) < DCD_HEADER_BYTES or head[4:8] != b"CORD":
        return {"readable": False, "bytes": size}
    nset, istart, nsavc = struct.unpack("<3i", head[8:20])
    delta = struct.unpack("<f", head[44:48])[0]
    has_cell = struct.unpack("<i", head[48:52])[0]
    n_atoms = struct.unpack("<i", head[268:272])[0]
    per = dcd_frame_bytes(n_atoms)
    body = size - DCD_HEADER_BYTES
    return {"readable": True, "bytes": size, "header_n_frames": nset, "header_first_step": istart,
            "header_step_interval": nsavc, "header_timestep_ps": round(delta * 0.04888821, 9),
            "header_has_cell": bool(has_cell), "n_atoms": n_atoms,
            "whole_frames_by_size": body // per, "bytes_beyond_whole_frames": body % per,
            "titles": [head[100 + 80 * k:180 + 80 * k].decode("ascii", "replace").strip() for k in (0, 1)]}


def read_dcd(path):
    """All whole frames of a DCD file: header, coordinates (n, atoms, 3; float32; Angstrom) and cells
    (n, 3, 3; rows a, b, c; Angstrom)."""
    h = read_dcd_header(path)
    if not h["readable"]:
        raise ValueError(f"{path} is not a DCD file of the expected layout")
    n, na = h["whole_frames_by_size"], h["n_atoms"]
    xyz = np.zeros((n, na, 3), dtype=np.float32)
    cells = np.zeros((n, 3, 3))
    with open(path, "rb") as f:
        f.seek(DCD_HEADER_BYTES)
        for k in range(n):
            rec = f.read(56)
            a, cosg, b, cosb, cosa, c = struct.unpack("<6d", rec[4:52])
            cells[k] = cell_vectors(a, b, c, cosa, cosb, cosg)
            for dim in range(3):
                f.read(4)
                xyz[k, :, dim] = np.frombuffer(f.read(4 * na), dtype="<f4")
                f.read(4)
    return h, xyz, cells


def cell_vectors(a, b, c, cos_alpha, cos_beta, cos_gamma):
    """Rows a, b, c of a cell from its lengths and the cosines of its angles (a along x, b in the xy plane)."""
    sin_gamma = np.sqrt(max(0.0, 1.0 - cos_gamma ** 2))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * cos_gamma, b * sin_gamma, 0.0])
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    vc = np.array([cx, cy, np.sqrt(max(0.0, c * c - cx * cx - cy * cy))])
    return np.array([va, vb, vc])


def read_inpcrd(path):
    """Coordinates (Angstrom) and the cell (rows a, b, c; Angstrom) of an Amber inpcrd file."""
    with open(path) as f:
        lines = f.read().split("\n")
    n = int(lines[1].split()[0])
    vals = []
    for ln in lines[2:]:
        vals.extend(float(ln[i:i + 12]) for i in range(0, len(ln), 12) if ln[i:i + 12].strip())
    xyz = np.array(vals[:3 * n]).reshape(n, 3)
    box = vals[3 * n:3 * n + 6]
    a, b, c, alpha, beta, gamma = box
    cell = cell_vectors(a, b, c, np.cos(np.radians(alpha)), np.cos(np.radians(beta)), np.cos(np.radians(gamma)))
    return xyz, cell


def read_state_xml(path):
    """Positions (Angstrom), velocities (Angstrom per ps) and cell (rows; Angstrom) of a saved OpenMM state."""
    import openmm
    from openmm import unit
    with open(path) as f:
        st = openmm.XmlSerializer.deserialize(f.read())
    pos = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
    try:
        vel = np.array(st.getVelocities(asNumpy=True).value_in_unit(unit.angstrom / unit.picosecond))
    except Exception:
        vel = None
    cell = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
    return pos, vel, cell


def read_pdb_atoms(path):
    """(serial, atom name, residue name, residue number, xyz) of the ATOM and HETATM lines of a PDB file."""
    atoms = []
    with open(path) as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                atoms.append((int(ln[6:11]), ln[12:16].strip(), ln[17:20].strip(), int(ln[22:26]),
                              (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return atoms


def read_log_csv(path):
    with open(path) as f:
        header = f.readline().strip().split(",")
        rows = [[float(v) for v in ln.strip().split(",")] for ln in f if ln.strip()]
    return header, np.array(rows).reshape(len(rows), len(header))


# ----------------------------------------------------------------------------- geometry

def min_image(d, cell):
    """Shortest periodic images of the vectors d (n, 3) in a cell with rows a, b, c in reduced form."""
    d = np.array(d, dtype=np.float64)
    d = d - np.outer(np.round(d[:, 2] / cell[2, 2]), cell[2])
    d = d - np.outer(np.round(d[:, 1] / cell[1, 1]), cell[1])
    d = d - np.outer(np.round(d[:, 0] / cell[0, 0]), cell[0])
    best, best_r2 = d.copy(), np.einsum("ij,ij->i", d, d)
    for i, j, k in itertools.product((-1, 0, 1), repeat=3):
        if i == j == k == 0:
            continue
        cand = d + i * cell[0] + j * cell[1] + k * cell[2]
        r2 = np.einsum("ij,ij->i", cand, cand)
        closer = r2 < best_r2
        best[closer], best_r2[closer] = cand[closer], r2[closer]
    return best


def dist(a, b):
    return float(np.linalg.norm(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)))


def angle(a, b, c):
    """Angle at b in degrees."""
    u = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    v = np.asarray(c, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))))


def torsion(a, b, c, d):
    """Torsion a-b-c-d in degrees, IUPAC sign, in (-180, 180]."""
    a, b, c, d = (np.asarray(x, dtype=np.float64) for x in (a, b, c, d))
    b1 = (c - b) / np.linalg.norm(c - b)
    v = (a - b) - np.dot(a - b, b1) * b1
    w = (d - c) - np.dot(d - c, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def kabsch(mobile, target):
    """Rotation R and translation t that put mobile (n, 3) on target (n, 3) by least squares: x' = x R^T + t."""
    mobile, target = np.asarray(mobile, dtype=np.float64), np.asarray(target, dtype=np.float64)
    cm, ct = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - cm).T @ (target - ct)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, ct - cm @ rot.T


def rmsd(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


LIG_NAMES = ("C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S")
LABELLINGS = {"identity": {}, "ring_flip": {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3"},
              "oxygen_swap": {"O1": "O2", "O2": "O1"},
              "ring_flip+oxygen_swap": {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3", "O1": "O2", "O2": "O1"}}


def ligand_rmsd_symmetric(lig_now, lig_ref):
    """Smallest heavy-atom RMSD over the four equivalent labellings; both arguments are dicts name -> xyz."""
    best = None
    for name, swap in LABELLINGS.items():
        v = rmsd([lig_now[a] for a in LIG_NAMES], [lig_ref[swap.get(a, a)] for a in LIG_NAMES])
        if best is None or v < best[0] - 1e-12:
            best = (v, name)
    return best


# ----------------------------------------------------------------------------- the plan of a run

def expected_steps(manifest):
    """What the run should hold, computed here from the protocol in the manifest, the equilibration scale
    and the segment length: per step the number of MD steps, the time step, the start time and the numbers of
    state-log rows, solute frames and full frames (counted inside each step, as md/RUN.md describes)."""
    proto, scale = manifest["protocol"], manifest["equil_scale"]
    out_ps = {k: proto["output"][k] for k in ("state_log_interval_ps", "solute_frame_interval_ps", "full_frame_interval_ps")}
    steps, t_ps = [{"tag": "min", "kind": "minimisation"}], 0.0
    for k, st in enumerate(proto["equilibration"]["stages"], start=1):
        dt = float(st["timestep_fs"])
        n = max(1, int(round(st["length_ps"] * 1000.0 / dt * scale)))
        every = {key: int(round(v * 1000.0 / dt)) for key, v in out_ps.items()}
        steps.append({"tag": f"eq{k:02d}", "kind": "equilibration", "name": st["name"], "timestep_fs": dt, "n_steps": n,
                      "length_ps": n * dt / 1000.0, "time_start_ps": t_ps,
                      "state_log_rows": n // every["state_log_interval_ps"],
                      "solute": n // every["solute_frame_interval_ps"], "full": n // every["full_frame_interval_ps"],
                      "every_steps": every})
        t_ps += n * dt / 1000.0
    dt = float(proto["production"]["timestep_fs"])
    seg_ns = manifest["segment_ns"]
    n = int(round(seg_ns * 1e6 / dt))
    every = {key: int(round(v * 1000.0 / dt)) for key, v in out_ps.items()}
    n_seg = int(round(manifest["production"]["target_ns"] / seg_ns))
    for s in range(1, n_seg + 1):
        steps.append({"tag": f"seg{s:04d}", "kind": "production", "timestep_fs": dt, "n_steps": n,
                      "length_ps": n * dt / 1000.0, "time_start_ps": (s - 1) * n * dt / 1000.0,
                      "state_log_rows": n // every["state_log_interval_ps"],
                      "solute": n // every["solute_frame_interval_ps"], "full": n // every["full_frame_interval_ps"],
                      "every_steps": every})
    return steps
