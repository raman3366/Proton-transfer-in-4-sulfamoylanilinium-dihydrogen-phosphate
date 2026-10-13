#!/usr/bin/env python
"""Simulation topologies of hCA II (6RL9) with the sulfonamidate anion on the zinc: stage S3 of md/PLAN.md.

    AMBERHOME=<prefix> python md/build_topology.py --model restrained --start X|N|R [--buffer-A 10.0]
    AMBERHOME=<prefix> python md/build_topology.py --model bonded --set 2e|2s --start X|N|R [--mcpb-dir DIR]
                                                   [--n2-frcmod FILE] [--variant NAME] [--compare-with MODEL] [--note TEXT]
    AMBERHOME=<prefix> python md/build_topology.py --model restrained --start X --min-check [--short-md] [--reuse-run]
    AMBERHOME=<prefix> python md/build_topology.py --model restrained --start X --min-check --min-full
    AMBERHOME=<prefix> python md/build_topology.py --model restrained --start X --k-scan
    AMBERHOME=<prefix> python md/build_topology.py --selftest-bonded [--start X] [--mcpb-dir DIR]
    AMBERHOME=<prefix> python md/build_topology.py --buffer-scan [--start X]

Reads (changes nothing there)
    md/system/protein_md.pdb, zinc.pdb, ligand_<start>.pdb, waters_<start>.pdb, residue_map.csv, starts.csv
    md/ligand/sfa_anion.mol2, sfa_anion_gaff2.frcmod, sfa_anion_repair.frcmod (loaded after the first, so
    that the repaired terms win), sfa_anion_repaired.prmtop (only to compare the ligand terms)
    model "bonded" only: the output of MCPB.py step 4 in --mcpb-dir (default for a named set: the first of
    /root/cqc-sfa-md/qm/gr/mcpb_<set> and md/qm/mcpb that holds a <group>_tleap.in):
    <group>_tleap.in, <group>_mcpbpy.pdb, <group>_mcpbpy.frcmod and the mol2 files named in <group>_tleap.in
    the checks: md/run_md.py, md/protocol.json, md/analysis/analyse_md.py (run, not changed)

The bonded model has two parameter sets (md/PLAN.md, 6 Oct 2026, 3:02 AM, point 3): --set 2e (MCPB.py step 2e,
empirical) and --set 2s (step 2s, Seminario).  The folder is named <model>_<set> (bonded_2e, bonded_2s), and the
builder stops if the frcmod of --mcpb-dir was not made by the step that the set names.  The terms at the ligand
nitrogen N2 (Zn-N2-S, Zn-N2-H, S-N2-H, H-N2-S-C4) are the same in both sets and are loaded last, from one file:
--n2-frcmod FILE (what gate G-R settles, md/GR.md; or n2_terms.frcmod of an earlier build).  Without it, and for
set 2e only, the builder writes that file itself with the terms AS BUILT BEFORE THE GATE: those of
<group>_mcpbpy.frcmod, but Zn-N2-H at its value in start X (107.41 deg) instead of the 120 deg of step 2e.
--variant NAME builds a set with OTHER terms at N2 (the file of --n2-frcmod may then also hold torsions about N2-S
and about Zn-N2 and one improper at N2); its folder is <model>_<set>_<NAME>, and the checks that assume the terms as
built (every torsion through the zinc zero, no improper there) ask instead for exactly what the file holds.
--compare-with MODEL compares the new system term by term with <out-root>/MODEL/<start>: only the terms in which the
two sets differ may differ (zinc bonds and angles between 2e and 2s; the terms at N2 between a set and its variant).
The files of a system are written to <start>/build/new/ and moved into place at the end, each by one rename: a build
that is refused or stops leaves the system that was there untouched.

Writes
    <out-root>/<model>/<start>/   system.prmtop, system.inpcrd (Amber, truncated octahedron), selections.json,
                                  extra_forces.json, reference.pdb, build.json, and build/ (work files of tleap);
                                  <model> is bonded_<set> for a named set; with --note also STATUS_NOTE.txt
    md/topology/<model>/<start>/  leap_stage1.in, leap.in, leap.log, ions.pdb, build.json (copy),
                                  selections_table.txt, openmm_check.json, SHA256SUMS; bonded: n2_terms.frcmod
    --min-check     <out-root>/_checks/<model>/<start>/min_<n>/ (or min_full/): a run folder of md/run_md.py;
                    md/topology/<model>/<start>/min_check.json (or min_check_full.json).  With --short-md the six
                    equilibration stages shortened to 4 ps and two production segments of 1 ps follow, and
                    md/analysis/analyse_md.py is run on the folder.
    --k-scan        md/topology/<model>/<start>/k_scan.json
    --buffer-scan   md/topology/buffer_scan_<start>.json
    --selftest-bonded  <out-root>/_selftest_bonded_mock/ (a MOCK, never to be simulated) and
                    md/topology/selftest_bonded_mock/; with --mcpb-dir: _selftest_bonded_synthetic/ in both places

How the system is made (two passes of tleap, AmberTools; leaprc.protein.ff14SB, leaprc.gaff2,
leaprc.water.tip3p):
  1. One PDB file is assembled from the four files of md/system/ (coordinates copied as text; residues
     numbered 1..N in the order protein, zinc, ligand, crystal water).
  2. Pass 1: tleap loads it, adds the hydrogens, and solvates with "solvateOct <unit> TIP3PBOX <buffer>".
  3. The ions: N0 = round(c * N_water / c_water) pairs of Na+ and Cl- plus the ions that neutralise the
     solute.  Each ion takes the place of one water molecule added in pass 1 (never a crystal water).  The
     molecules are chosen by a fixed rule without chance: see choose_ion_sites().
  4. Pass 2: the same input as pass 1, then the chosen water molecules are removed and the ions added at the
     positions of their oxygen atoms.  The builder verifies that every other atom is where pass 1 put it.
  5. The date in the first line of system.prmtop is replaced by a fixed stamp, so that a second run gives the
     same file (the driver records the sha256 of the topology and refuses to continue a run on another one).

Units: Angstrom in PDB, mol2, inpcrd and in every key that ends in _A; nm and kJ mol-1 nm-2 in
extra_forces.json (stated in the key names); kcal mol-1 where a key says so (Amber force-field units).
Indices in selections.json and extra_forces.json are zero-based atom indices of system.prmtop.
No random numbers are used anywhere.  A second build gives the same six files of the system (checked by
sha256).  The check runs use the OpenMM CPU platform with two threads and repeat only approximately.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import itertools
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
MD = REPO / "md"
SYSTEM = MD / "system"
LIGAND = MD / "ligand"
RECORDS = MD / "topology"
AMBERHOME = os.environ.get("AMBERHOME") or sys.prefix
DEFAULT_OUT_ROOT = os.environ.get("CQC_SFA_MD_SYSTEMS", "/root/cqc-sfa-md/systems")
THREADS = 2

# ---- settings fixed for the project (each is recorded in build.json)
BUFFER_A = 10.0                      # argument of solvateOct
SALT_MOL_PER_L = 0.15
WATER_MOL_PER_L = 55.3               # liquid water at 300 K and 1 bar: 0.9965 g cm-3 / 18.015 g mol-1
ION_SEED = 20261005                  # enters only the fixed order in which water molecules are offered
ION_MIN_DIST_TO_SOLUTE_A = 6.0       # ion to any atom of protein, zinc, ligand or crystal water
ION_MIN_SEPARATION_A = 6.0           # ion to ion, minimum image
KT_300_KJ = 0.0083144626 * 300.0     # kJ mol-1
ZN_RESTRAINT_K_KJ_NM2 = 500000.0     # restrained model: U = 0.5 k (r - r0)^2 on each of the four Zn-N distances
FIXED_PRMTOP_DATE = "DATE = 01/01/00  00:00:00"

LIG_HEAVY = ["C1", "C2", "C3", "C4", "C5", "C6", "N1", "N2", "O1", "O2", "S"]
LIG_H = ["H2", "H3", "H5", "H6", "HN11", "HN12", "HN21"]
LABELLINGS = [("identity", {}), ("ring_flip", {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3"}),
              ("oxygen_swap", {"O1": "O2", "O2": "O1"}),
              ("ring_flip+oxygen_swap", {"C2": "C6", "C6": "C2", "C3": "C5", "C5": "C3", "O1": "O2", "O2": "O1"})]
# key of selections.json -> (deposited residue number, atom name)
SITE = {"his94_NE2": (94, "NE2"), "his96_NE2": (96, "NE2"), "his119_ND1": (119, "ND1"),
        "thr199_OG1": (199, "OG1"), "thr199_HG1": (199, "HG1"), "thr199_N": (199, "N"), "thr199_H": (199, "H"),
        "glu106_OE1": (106, "OE1"), "glu106_OE2": (106, "OE2"), "phe131_CZ": (131, "CZ"),
        "asn67_OD1": (67, "OD1"), "asn67_ND2": (67, "ND2")}
ZN_PARTNERS = [("N2", None, "N2 of the ligand (SAN 304)"), ("NE2", 94, "NE2 of His94"),
               ("NE2", 96, "NE2 of His96"), ("ND1", 119, "ND1 of His119")]
# the bonded parameter sets (md/PLAN.md, 6 Oct 2026, 3:02 AM, point 3): the word MCPB.py writes behind every zinc
# term of its frcmod ("Created by <method> method using MCPB.py"), and where the MCPB.py folder of a set is looked for
SET_METHOD = {"2e": "empirical", "2s": "Seminario"}
MCPB_SET_DIRS = (str(pathlib.Path(DEFAULT_OUT_ROOT).parent / "qm" / "gr" / "mcpb_{set}"), str(MD / "qm" / "mcpb"))
LIG_EQUAL_PAIRS = [("O1", "O2"), ("C2", "C6"), ("C3", "C5"), ("H2", "H6"), ("H3", "H5"), ("HN11", "HN12")]
SCEE, SCNB = 1.2, 2.0                # Amber's scaling of 1-4 electrostatics and 1-4 Lennard-Jones (ff14SB, GAFF2)
# limits of the minimisation check, fixed before the first minimisation was run
MIN_CHECK = {"zn_distance_tolerance_A": 0.10,
             "max_displacement_zinc_and_ligand_heavy_A": {"X": 0.30, "N": 0.50, "R": 0.50},
             "why": "heavy atoms are held at k = 1000 kJ mol-1 nm-2; in N and R the docked Zn-N2 distance "
                    "(2.21, 2.24 A) is drawn to 1.93 A and the S-N bond of the neutral molecule (1.677 A) "
                    "shortens, so zinc and ligand must move by a few tenths of an Angstrom"}


class BuildError(Exception):
    pass


# --------------------------------------------------------------------------------------------- small tools
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path, obj):
    pathlib.Path(path).write_text(json.dumps(obj, indent=1) + "\n")


def rel(path):
    p = pathlib.Path(path).resolve()
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def now_ist():
    t = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=5, minutes=30)))
    return t.strftime("%d %b %Y, %I:%M %p IST").lstrip("0")


def dist(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def angle(a, b, c):
    u, v = np.asarray(a) - np.asarray(b), np.asarray(c) - np.asarray(b)
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))))


def torsion(a, b, c, d):
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    b1 = (c - b) / np.linalg.norm(c - b)
    v = (a - b) - np.dot(a - b, b1) * b1
    w = (d - c) - np.dot(d - c, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def kabsch(P, Q):
    """Rotation R and translation t with Q ~ P @ R.T + t (least squares, no reflection)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def box_vectors(box):
    """Cell vectors (Angstrom) of an Amber box (a, b, c, alpha, beta, gamma): a along x, b in the xy plane."""
    a, b, c, al, be, ga = [float(x) for x in box]
    al, be, ga = (math.radians(x) for x in (al, be, ga))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(ga), b * math.sin(ga), 0.0])
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    vc = np.array([cx, cy, math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))])
    return np.array([va, vb, vc])


def image_shifts(vecs):
    return [n[0] * vecs[0] + n[1] * vecs[1] + n[2] * vecs[2]
            for n in itertools.product((-1, 0, 1), repeat=3) if n != (0, 0, 0)]


def read_pdb(path):
    atoms = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                resi = int(line[22:26])
            except ValueError:                       # tleap writes residue numbers above 9999 in another notation
                resi = None
            atoms.append({"line": line.ljust(80), "name": line[12:16].strip(), "resn": line[17:20].strip(),
                          "resi": resi,
                          "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])})
    return atoms


def versions():
    import parmed
    import scipy
    amber = sorted(glob.glob(os.path.join(AMBERHOME, "conda-meta", "ambertools-*.json")))
    return {"ambertools": os.path.basename(amber[-1])[:-5] if amber else "unknown (no conda-meta entry)",
            "tleap": os.path.join(AMBERHOME, "bin", "tleap"), "python": sys.version.split()[0],
            "parmed": parmed.__version__, "numpy": np.__version__, "scipy": scipy.__version__}


def nice_prefix():
    """Child processes run at niceness 10; nothing is added if this process is already there (nice adds up)."""
    return ["nice", "-n", str(10 - os.nice(0))] if os.nice(0) < 10 else []


def run_tleap(script_path, cwd):
    env = dict(os.environ, AMBERHOME=AMBERHOME, OMP_NUM_THREADS=str(THREADS))
    p = subprocess.run(nice_prefix() + [os.path.join(AMBERHOME, "bin", "tleap"), "-f", str(script_path)],
                       cwd=str(cwd), env=env, capture_output=True, text=True)
    log = p.stdout + p.stderr
    m = re.search(r"Exiting LEaP: Errors = (\d+); Warnings = (\d+); Notes = (\d+)", log)
    if p.returncode != 0 or not m or int(m.group(1)) != 0:
        pathlib.Path(str(script_path) + ".failed.log").write_text(log)
        raise BuildError(f"tleap failed for {script_path} (log: {script_path}.failed.log)")
    return log, {"errors": int(m.group(1)), "warnings": int(m.group(2)), "notes": int(m.group(3))}


def leap_warnings(log):
    """The warnings of a tleap log, sorted into the expected kinds and the rest."""
    lines = log.splitlines()
    out = {"close_contacts_with_a_hydrogen": [], "close_contacts_between_heavy_atoms": [],
           "charge_not_zero": [], "c_terminal_name_for_savepdb": [], "one_sided_connection": [], "other": []}
    for k, line in enumerate(lines):
        if not line.rstrip().endswith("Warning!") or k + 1 >= len(lines):
            continue
        msg = lines[k + 1].strip()
        m = re.match(r"Close contact of ([\d.]+) angstroms between nonbonded atoms (\S+) and (\S+)", msg)
        if m:
            nxt = lines[k + 2] if k + 2 < len(lines) else ""
            where = re.search(r"\.R<(\S+) (\d+)>\.A<(\S+) \d+> and \.R<(\S+) (\d+)>\.A<(\S+) \d+>", nxt)
            item = {"distance_A": float(m.group(1)),
                    "atoms": f"{where.group(1)}{where.group(2)}:{where.group(3)} - "
                             f"{where.group(4)}{where.group(5)}:{where.group(6)}" if where else msg}
            heavy = not (m.group(2).startswith("H") or m.group(3).startswith("H"))
            out["close_contacts_between_heavy_atoms" if heavy else "close_contacts_with_a_hydrogen"].append(item)
        elif "charge of the unit" in msg and "is not zero" in msg:
            out["charge_not_zero"].append(msg)
        elif msg.startswith("Converting C-terminal residue name to PDB format"):
            out["c_terminal_name_for_savepdb"].append(msg)
        elif msg.startswith("One sided connection."):
            out["one_sided_connection"].append(msg)
        else:
            out["other"].append(msg)
    return out


def solvation_lines(log):
    keep = ("Scaling up box", "Solute vdw bounding box", "Total bounding box", "box expansion", "Volume:",
            "Total mass", "Added ")
    return [ln.strip() for ln in log.splitlines() if any(ln.strip().startswith(k) or k in ln[:40] for k in keep)
            and "missing" not in ln]


# ------------------------------------------------------------------------------------------- force-field files
def read_frcmod(path):
    """Sections of an Amber frcmod file: {'BOND': {(t1, t2): (k, r0)}, 'ANGLE': {...}, 'DIHE': {types: [(k, phase, n)]},
    'IMPROPER': {types as written, the third is the central atom: (k, phase, n)}}."""
    sec, out = None, {"BOND": {}, "ANGLE": {}, "DIHE": {}, "IMPROPER": {}}
    for line in pathlib.Path(path).read_text().splitlines()[1:]:
        word = line.strip()
        if word[:4] in ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB"):
            sec = {"ANGL": "ANGLE", "IMPR": "IMPROPER", "NONB": "NONBON"}.get(word[:4], word[:4])
            continue
        if not word or sec not in out:
            continue
        n = {"BOND": 2, "ANGLE": 3, "DIHE": 4, "IMPROPER": 4}[sec]
        types = tuple(t.strip() for t in line[:3 * n - 1].split("-"))
        num = line[3 * n - 1:].split()
        if sec == "BOND":
            out[sec][types] = (float(num[0]), float(num[1]))
        elif sec == "ANGLE":
            out[sec][types] = (float(num[0]), float(num[1]))
        elif sec == "IMPROPER":
            out[sec][types] = (float(num[0]), float(num[1]), abs(int(float(num[2]))))
        else:
            out[sec].setdefault(types, []).append((float(num[1]) / float(num[0]), float(num[2]), abs(int(float(num[3])))))
    return out


def frcmod_nonbon(path, atom_type):
    """Lennard-Jones line of an atom type in the NONBON section of a frcmod file (Rmin/2 in A, epsilon in kcal mol-1)."""
    inside = False
    for line in pathlib.Path(path).read_text().splitlines():
        w = line.split()
        if w and w[0][:4] == "NONB":
            inside = True
            continue
        if inside and len(w) >= 3 and w[0] == atom_type:
            try:
                return {"rmin_half_A": float(w[1]), "epsilon_kcal_per_mol": float(w[2]), "line": line.strip()}
            except ValueError:
                continue
    return None


def ligand_terms(parm, atoms):
    """Bonded terms among the given atoms of a ParmEd topology, keyed by atom names (direction-free)."""
    ids = {a.idx for a in atoms}

    def key(*aa):
        names = tuple(a.name for a in aa)
        return names if names <= names[::-1] else names[::-1]

    out = {"bonds": {}, "angles": {}, "dihedrals": {}, "impropers": {}}
    for b in parm.bonds:
        if b.atom1.idx in ids and b.atom2.idx in ids:
            out["bonds"]["-".join(key(b.atom1, b.atom2))] = [round(b.type.k, 4), round(b.type.req, 5)]
    for a in parm.angles:
        if all(x.idx in ids for x in (a.atom1, a.atom2, a.atom3)):
            out["angles"]["-".join(key(a.atom1, a.atom2, a.atom3))] = [round(a.type.k, 4), round(a.type.theteq, 4)]
    for d in parm.dihedrals:
        aa = (d.atom1, d.atom2, d.atom3, d.atom4)
        if all(x.idx in ids for x in aa):
            term = [round(d.type.phi_k, 5), int(d.type.per), round(d.type.phase, 3) % 360.0,
                    round(d.type.scee, 4), round(d.type.scnb, 4)]
            if d.improper:
                out["impropers"].setdefault("-".join(key(*aa)), []).append(term[:3])
            else:
                out["dihedrals"].setdefault("-".join(key(*aa)), []).append(term)
    for sec in ("dihedrals", "impropers"):
        for k in out[sec]:
            out[sec][k].sort()
    return out


# ------------------------------------------------------------------------------------------------- MCPB.py
def read_mcpb(mcpb_dir, group=None):
    """What the builder needs from the output of MCPB.py step 4 (model "bonded")."""
    d = pathlib.Path(mcpb_dir)
    cands = sorted(d.glob(f"{group}_tleap.in" if group else "*_tleap.in"))
    if len(cands) != 1:
        raise BuildError(f"model 'bonded' needs exactly one <group>_tleap.in in {d} (found {len(cands)}); it is "
                         "written by 'MCPB.py -i <group>.in -s 4'. See md/TOPOLOGY.md, section 5.")
    tl = cands[0]
    group = tl.name[:-len("_tleap.in")]
    text = tl.read_text()
    out = {"dir": str(d), "group": group, "tleap_in": str(tl), "sources": re.findall(r"(?m)^source\s+(\S+)", text),
           "atom_types": re.findall(r'\{\s*"(\S+?)"\s+"(\S+?)"\s+"(\S+?)"\s*\}', text),
           "mol2": re.findall(r"(?m)^(\S+)\s*=\s*loadmol2\s+(\S+)", text),
           "frcmods": re.findall(r"(?mi)^loadamberparams\s+(\S+)", text),
           "pdb": re.findall(r"(?mi)^(\S+)\s*=\s*loadpdb\s+(\S+)", text),
           "bonds": re.findall(r"(?mi)^bond\s+(\S+?)\.(\d+)\.(\S+)\s+(\S+?)\.(\d+)\.(\S+)", text)}
    if len(out["pdb"]) != 1:
        raise BuildError(f"{tl}: one 'loadpdb' line expected")
    need = [d / out["pdb"][0][1], d / f"{group}_mcpbpy.frcmod"] + [d / f for _, f in out["mol2"]]
    missing = [str(p) for p in need if not p.is_file()]
    if missing:
        raise BuildError("files of MCPB.py step 4 are missing: " + ", ".join(missing))
    if not out["mol2"] or not out["bonds"]:
        raise BuildError(f"{tl}: no loadmol2 or no bond command; is this the bonded model of MCPB.py?")
    residues = {}
    for a in read_pdb(d / out["pdb"][0][1]):
        r = residues.setdefault(a["resi"], {"name": a["resn"], "atoms": {}})
        r["atoms"][a["name"]] = a["xyz"]
    out["residues"] = residues
    out["frcmod"] = d / f"{group}_mcpbpy.frcmod"
    out["methods"] = sorted(set(re.findall(r"Created by (\S+) method using MCPB\.py", out["frcmod"].read_text())))
    return out


def frcmod_lines(path, section, types):
    """The lines of one section of a frcmod file that carry these atom types (in either direction), as written."""
    n = {"BOND": 2, "ANGLE": 3, "DIHE": 4}[section]
    want = {tuple(types), tuple(types)[::-1]}
    sec, hits = None, []
    for line in pathlib.Path(path).read_text().splitlines()[1:]:
        word = line.strip()
        if word[:4] in ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB"):
            sec = {"ANGL": "ANGLE"}.get(word[:4], word[:4])
            continue
        if word and sec == section and tuple(t.strip() for t in line[:3 * n - 1].split("-")) in want:
            hits.append(line.rstrip())
    return hits


def n2_terms_file(mcpb, types, theta_x, path):
    """The terms at the ligand nitrogen N2 AS BUILT BEFORE GATE G-R, in a frcmod file of their own that is loaded
    last: Zn-N2-S, S-N2-H and the H-N2-S-C4 series as <group>_mcpbpy.frcmod has them (copied line by line), and
    Zn-N2-H with the force constant of that file and its angle in start X instead of the 120 deg that MCPB.py
    step 2e writes (md/PLAN.md, 6 Oct 2026, 3:02 AM point 2 and 8:01 AM point 3; md/qm/gr_mm.py does the same for
    the small model of the gate).  types: atom types of Zn, N2, S, HN21 and C4."""
    zn, n2, s, hn, c4 = (types[k] for k in ("ZN", "N2", "S", "HN21", "C4"))
    src = mcpb["frcmod"]
    zn_s, s_h = frcmod_lines(src, "ANGLE", (zn, n2, s)), frcmod_lines(src, "ANGLE", (hn, n2, s))
    zn_h, tor = frcmod_lines(src, "ANGLE", (zn, n2, hn)), frcmod_lines(src, "DIHE", (c4, s, n2, hn))
    if not (len(zn_s) == len(s_h) == len(zn_h) == 1 and tor):
        raise BuildError(f"{src}: the terms at N2 ({zn}-{n2}-{s}, {hn}-{n2}-{s}, {zn}-{n2}-{hn}, {c4}-{s}-{n2}-{hn}) "
                         "were not found once each")
    k_h = float(zn_h[0][8:].split()[0])
    text = [f"Terms at N2 AS BUILT BEFORE GATE G-R (md/build_topology.py): Zn-N2-H at its value in start X, not the "
            f"{float(zn_h[0][8:].split()[1]):g} deg of {src.name}",
            "MASS", "", "BOND", "", "ANGLE", zn_s[0],
            f"{zn}-{n2}-{hn}   {k_h:6.2f}    {theta_x:7.3f}    angle of start X (md/system); force constant of {src.name}",
            s_h[0], "", "DIHE"] + tor + ["", "IMPROPER", "", "NONBON", "", ""]
    pathlib.Path(path).write_text("\n".join(text))
    return {"replaced_line_of_mcpb_frcmod": zn_h[0].strip(), "lines_copied_from_mcpb_frcmod": [x.strip() for x in zn_s + s_h + tor]}


def zinc_pair_sets(parm, zn):
    """Pairs of atoms whose shortest path of bonds runs through the zinc, by the number of bonds between them
    (1 to 4).  The zinc has four arms (its nitrogen atoms with what hangs on them); a pair is the zinc with an
    atom of an arm, or two atoms of different arms.  Returns {n_bonds: [(i, j), ...]} with i < j."""
    arms = []
    for n in parm.atoms[zn].bond_partners:
        depth, layer = {n.idx: 1}, [n]
        for d in (2, 3, 4):
            nxt = []
            for a in layer:
                for b in a.bond_partners:
                    if b.idx != zn and b.idx not in depth:
                        depth[b.idx] = d
                        nxt.append(b)
            layer = nxt
        arms.append(depth)
    pairs = {1: [], 2: [], 3: [], 4: []}
    for depth in arms:
        for i, d in depth.items():
            pairs[d].append(tuple(sorted((zn, i))))
    for p, q in itertools.combinations(range(len(arms)), 2):
        if set(arms[p]) & set(arms[q]):
            raise BuildError("two arms of the zinc share an atom within four bonds: the pair rules cannot be checked")
        for i, di in arms[p].items():
            for j, dj in arms[q].items():
                if di + dj <= 4:
                    pairs[di + dj].append(tuple(sorted((i, j))))
    return {d: sorted(set(v)) for d, v in pairs.items()}


# --------------------------------------------------------------------------------------------- the assembly
def load_inputs(start):
    files = {"protein": SYSTEM / "protein_md.pdb", "zinc": SYSTEM / "zinc.pdb",
             "ligand": SYSTEM / f"ligand_{start}.pdb", "waters": SYSTEM / f"waters_{start}.pdb",
             "ligand_crystal": SYSTEM / "ligand_X.pdb", "residue_map": SYSTEM / "residue_map.csv",
             "starts": SYSTEM / "starts.csv", "mol2": LIGAND / "sfa_anion.mol2",
             "frcmod_gaff2": LIGAND / "sfa_anion_gaff2.frcmod", "frcmod_repair": LIGAND / "sfa_anion_repair.frcmod",
             "ligand_repaired_prmtop": LIGAND / "sfa_anion_repaired.prmtop"}
    for k, p in files.items():
        if not p.is_file():
            raise BuildError(f"input missing: {p}")
    resmap = {}
    with open(files["residue_map"]) as f:
        for row in csv.DictReader(f):
            resmap[int(row["tleap_index"])] = {"deposited": int(row["deposited_number"]), "name": row["name_md"]}
    with open(files["starts"]) as f:
        starts = {row["start"]: row for row in csv.DictReader(f)}
    prot = read_pdb(files["protein"])
    residues = []                                    # protein residues in file order
    for a in prot:
        if not residues or residues[-1]["deposited"] != a["resi"]:
            residues.append({"deposited": a["resi"], "name": a["resn"], "atoms": []})
        residues[-1]["atoms"].append(a)
    if len(residues) != len(resmap):
        raise BuildError("protein_md.pdb and residue_map.csv differ in the number of residues")
    for k, r in enumerate(residues, start=1):
        if resmap[k]["deposited"] != r["deposited"] or resmap[k]["name"] != r["name"]:
            raise BuildError(f"residue {k}: protein_md.pdb has {r['name']}{r['deposited']}, residue_map.csv "
                             f"{resmap[k]['name']}{resmap[k]['deposited']}")
    zn = read_pdb(files["zinc"])
    lig = read_pdb(files["ligand"])
    wat = read_pdb(files["waters"])
    ligx = read_pdb(files["ligand_crystal"])
    if len(zn) != 1 or [a["name"] for a in lig] != LIG_HEAVY + LIG_H or [a["name"] for a in ligx] != LIG_HEAVY + LIG_H:
        raise BuildError("zinc.pdb or the ligand files are not as expected (1 atom; 18 atoms in the order of the mol2 file)")
    if any(a["name"] != "O" for a in wat):
        raise BuildError("the water file should hold oxygen atoms only")
    return {"files": files, "resmap": resmap, "starts": starts, "protein": residues, "zinc": zn, "ligand": lig,
            "waters": wat, "ligand_crystal": ligx}


def write_assembly(inp, path, rename):
    """One PDB file: protein, zinc, ligand, crystal water; residues 1..N; coordinates copied as text.
    rename: {residue ordinal (1-based): new residue name} for the zinc-site residues of the bonded model."""
    out, table, serial = [], [], 0

    def put(rec, atom, resn, chain, ordinal):
        nonlocal serial
        serial += 1
        ln = atom["line"]
        out.append("%-6s%5d %4s %3s %1s%4d    %24s%12s          %2s" %
                   (rec, serial, ln[12:16], resn, chain, ordinal, ln[30:54], ln[54:66], ln[76:78]))

    ordinal = 0
    for r in inp["protein"]:
        ordinal += 1
        name = rename.get(ordinal, r["name"])
        for a in r["atoms"]:
            put("ATOM", a, name, "A", ordinal)
        table.append({"ordinal": ordinal, "kind": "protein", "name": name, "name_md": r["name"],
                      "deposited": r["deposited"]})
    out.append("TER")
    for kind, atoms, dep in (("zinc", inp["zinc"], 301), ("ligand", inp["ligand"], 304)):
        ordinal += 1
        name = rename.get(ordinal, atoms[0]["resn"])
        for a in atoms:
            put("HETATM", a, name, " ", ordinal)
        out.append("TER")
        table.append({"ordinal": ordinal, "kind": kind, "name": name, "name_md": atoms[0]["resn"], "deposited": dep})
    for a in inp["waters"]:
        ordinal += 1
        put("HETATM", a, "HOH", " ", ordinal)
        out.append("TER")
        table.append({"ordinal": ordinal, "kind": "crystal_water", "name": "HOH", "name_md": "HOH", "deposited": a["resi"]})
    out.append("END")
    pathlib.Path(path).write_text("\n".join(out) + "\n")
    return table


def leap_head(model, mcpb, assembly, build, buffer_A=BUFFER_A, last_frcmods=(), chain_units=()):
    """The lines of the tleap input up to and including the solvation.  last_frcmods (bonded model): parameter
    files loaded after those of MCPB.py, so that their terms are the ones in force (the terms at N2).
    chain_units (bonded model): the units of MCPB.py that are residues of the protein chain.  They get head, tail
    and connect atoms (N and C), so that loadpdb joins them to their neighbours as it joins library residues.
    Without that the peptide bonds are made by the 'bond' commands only AFTER loadpdb has built the hydrogens, and
    the amide hydrogen of the residue after each of them is built without its C-N bond, out of the peptide plane
    (found on 6 Oct 2026 with the real output: Phe95, Trp97, Leu120, by 0.4 to 1.2 A).  The peptide 'bond'
    commands of MCPB.py for these units are then left out (the builder checks that the chain is unbroken)."""
    L = ["source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p"]
    notes = []
    if model == "restrained":
        L += [f"loadamberparams {LIGAND / 'sfa_anion_gaff2.frcmod'}",
              f"loadamberparams {LIGAND / 'sfa_anion_repair.frcmod'}",
              f"SAN = loadmol2 {LIGAND / 'sfa_anion.mol2'}",
              f"sys = loadpdb {assembly}"]
    else:
        d = pathlib.Path(mcpb["dir"])
        if mcpb["atom_types"]:
            L.append("addAtomTypes {")
            L += ['    { "%s" "%s" "%s" }' % t for t in mcpb["atom_types"]]
            L.append("}")
        for unit, fname in mcpb["mol2"]:
            L.append(f"{unit} = loadmol2 {d / fname}")
        for unit in chain_units:
            L += [f"set {unit} head {unit}.1.N", f"set {unit} tail {unit}.1.C",
                  f"set {unit}.1 connect0 {unit}.1.N", f"set {unit}.1 connect1 {unit}.1.C"]
        for fname in mcpb["frcmods"]:
            if (d / fname).is_file():
                L.append(f"loadamberparams {d / fname}")
            elif fname.startswith("frcmod.ions"):
                notes.append(f"{fname} of {mcpb['group']}_tleap.in is not loaded: Na+ and Cl- keep the parameters of "
                             "leaprc.water.tip3p (Joung-Cheatham), as in the restrained model")
            else:
                L.append(f"loadamberparams {fname}")
        L += [f"loadamberparams {f}" for f in last_frcmods]
        L.append(f"sys = loadpdb {assembly}")
        for u1, r1, a1, u2, r2, a2 in mcpb["bonds"]:
            if chain_units and a1 == "C" and a2 == "N" and int(r2) == int(r1) + 1:
                continue                             # a peptide bond: made by loadpdb through the head and tail atoms
            L.append(f"bond sys.{r1}.{a1} sys.{r2}.{a2}")
    L += ["check sys", "charge sys",
          f"saveamberparm sys {build / 'dry.prmtop'} {build / 'dry.inpcrd'}",
          f"solvateOct sys TIP3PBOX {buffer_A:.1f}"]
    return L, notes


# ------------------------------------------------------------------------------------------------- the ions
def choose_ion_sites(xyz, o_bulk, solute_atoms, n_ions, vecs):
    """Water molecules that give way to ions.  Rule, without chance:
    1. candidates are the water molecules added by solvateOct whose oxygen is at least
       ION_MIN_DIST_TO_SOLUTE_A from every atom of protein, zinc, ligand and crystal water;
    2. they are put in the order of sha256("cqc-sfa-md ion sites|<ION_SEED>|<k>"), k being the place of the
       molecule among the added water molecules in the topology of pass 1 (a fixed pseudo-random order);
    3. going through that order, a molecule is taken if its oxygen is at least ION_MIN_SEPARATION_A (minimum
       image) from every oxygen already taken, until n_ions are taken.
    Returns the positions (in o_bulk) of the molecules taken, in the order taken."""
    from scipy.spatial import cKDTree
    d, _ = cKDTree(xyz[solute_atoms]).query(xyz[o_bulk], k=1)
    cand = [k for k in range(len(o_bulk)) if d[k] >= ION_MIN_DIST_TO_SOLUTE_A]
    cand.sort(key=lambda k: hashlib.sha256(f"cqc-sfa-md ion sites|{ION_SEED}|{k}".encode()).hexdigest())
    shifts = [np.zeros(3)] + image_shifts(vecs)
    taken = []
    for k in cand:
        p = xyz[o_bulk[k]]
        ok = True
        for j in taken:
            q = xyz[o_bulk[j]]
            if min(np.linalg.norm(p - q + s) for s in shifts) < ION_MIN_SEPARATION_A:
                ok = False
                break
        if ok:
            taken.append(k)
            if len(taken) == n_ions:
                break
    if len(taken) < n_ions:
        raise BuildError(f"only {len(taken)} of {n_ions} ion sites found")
    return taken, len(cand)


# --------------------------------------------------------------------------------- two built systems compared
ZINC_SET_CLASSES = ["bond Zn-N", "angle N-Zn-N", "angle Zn-N-C"]
N2_VARIANT_CLASSES = ["angle at N2", "torsion about N2-S", "torsion about Zn-N2", "improper at N2"]


def prmtop_sections(path):
    """The %FLAG sections of an Amber prmtop as text (the first line with its date left out)."""
    out, flag = {}, None
    for line in pathlib.Path(path).read_text().splitlines()[1:]:
        if line.startswith("%FLAG"):
            flag = line.split()[1]
            out[flag] = []
        elif flag is not None:
            out[flag].append(line)
    return {k: "\n".join(v) for k, v in out.items()}


def topology_tables(parm):
    """Every term of a ParmEd topology in a form in which two systems of the same atoms can be compared."""
    t = {"atoms": [(a.name, a.residue.name, a.type, round(a.charge, 6), round(a.mass, 4), round(a.rmin, 5),
                    round(a.epsilon, 8)) for a in parm.atoms],
         "bond": {}, "angle": {}, "torsion": {}, "improper": {}, "pairs_1_4": {}}
    for b in parm.bonds:
        i, j = sorted((b.atom1.idx, b.atom2.idx))
        t["bond"][(i, j)] = (round(b.type.k, 4), round(b.type.req, 5))
    for a in parm.angles:
        i, j, k = a.atom1.idx, a.atom2.idx, a.atom3.idx
        t["angle"][(i, j, k) if i < k else (k, j, i)] = (round(a.type.k, 4), round(a.type.theteq, 4))
    for d in parm.dihedrals:
        q = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
        term = (round(d.type.phi_k, 5), int(d.type.per), round(d.type.phase, 3) % 360.0)
        if d.improper:                               # Amber: the third atom is the central one
            t["improper"].setdefault((tuple(sorted((q[0], q[1], q[3]))), q[2]), []).append(term)
        else:
            key = q if q[0] < q[3] else q[::-1]
            t["torsion"].setdefault(key, [])
            if abs(term[0]) > 1e-9:
                t["torsion"][key].append(term)
            if not d.ignore_end:
                t["pairs_1_4"][tuple(sorted((q[0], q[3])))] = (round(d.type.scee, 4), round(d.type.scnb, 4))
    for sec in ("torsion", "improper"):
        for k in t[sec]:
            t[sec][k].sort()
    t["excluded"] = (list(parm.parm_data["NUMBER_EXCLUDED_ATOMS"]), list(parm.parm_data["EXCLUDED_ATOMS_LIST"]))
    return t


def term_class(kind, key, sel):
    """Where a bonded term sits: at the zinc, at N2, or elsewhere."""
    zn, n2, s_ = sel["zn"], sel["lig_heavy"]["N2"], sel["lig_heavy"]["S"]
    his = {sel["his94_NE2"], sel["his96_NE2"], sel["his119_ND1"]}
    if kind == "bond":
        return "bond Zn-N" if zn in key else "bond elsewhere"
    if kind == "angle":
        if key[1] == zn:
            return "angle N-Zn-N"
        if key[1] == n2:
            return "angle at N2"
        if key[1] in his and zn in (key[0], key[2]):
            return "angle Zn-N-C"
        return "angle elsewhere"
    if kind == "torsion":
        mid = {key[1], key[2]}
        return "torsion about N2-S" if mid == {n2, s_} else "torsion about Zn-N2" if mid == {zn, n2} else "torsion elsewhere"
    return "improper at N2" if key[1] == n2 else "improper elsewhere"


def compare_systems(parm, new_dir, other_dir, sel, this_set, other_set):
    """The system just built (topology parm, files in new_dir) against the built system in other_dir, term by term.
    this_set, other_set: (parameter set, variant).  Terms may differ only where the two sets differ: the bonds and
    angles of the zinc (other parameter set) and the terms at N2 (other variant); everything else must be the same."""
    import parmed
    other = parmed.load_file(str(other_dir / "system.prmtop"), xyz=str(other_dir / "system.inpcrd"))
    xyz_new = parmed.load_file(str(new_dir / "system.inpcrd")).coordinates.reshape(-1, 3)
    a, b = topology_tables(parm), topology_tables(other)
    allowed = (ZINC_SET_CLASSES if this_set[0] != other_set[0] else []) + (N2_VARIANT_CLASSES if this_set[1] != other_set[1] else [])
    res = {"other_system": str(other_dir), "this_set": list(this_set), "other_set": list(other_set),
           "classes_that_may_differ": allowed, "same_number_of_atoms": len(a["atoms"]) == len(b["atoms"])}
    if not res["same_number_of_atoms"]:
        res.update({"only_the_allowed_terms_differ": False, "note": "the two systems do not have the same atoms"})
        return res

    def name(i):
        at = parm.atoms[i]
        return "Zn" if i == sel["zn"] else f"{at.residue.name}{at.residue.idx + 1}:{at.name}"

    fields = ("name", "residue name", "type", "charge", "mass", "Rmin/2", "epsilon")
    atom_diff = {f: sum(1 for x, y in zip(a["atoms"], b["atoms"]) if x[n] != y[n]) for n, f in enumerate(fields)}
    type_diff = sorted({(x[2], y[2]) for x, y in zip(a["atoms"], b["atoms"]) if x[2] != y[2]})
    counts, rows, classes = {}, [], {}
    for kind in ("bond", "angle", "torsion", "improper"):
        ka, kb = a[kind], b[kind]
        counts[kind] = {"this": len(ka), "other": len(kb), "in_one_only": len(set(ka) ^ set(kb)),
                        "different": 0}
        for key in sorted(set(ka) | set(kb), key=str):
            if ka.get(key) == kb.get(key):
                continue
            if kind == "torsion" and not ka.get(key) and not kb.get(key):
                continue                             # no term with a coefficient on either side
            counts[kind]["different"] += 1
            ckey = (key[0][0], key[1]) if kind == "improper" else key
            cls = term_class(kind, ckey, sel)
            classes[cls] = classes.get(cls, 0) + 1
            flat = list(key[0]) + [key[1]] if kind == "improper" else list(key)
            rows.append({"kind": kind, "class": cls, "atoms": "-".join(name(i) for i in flat),
                         "this": ka.get(key), "other": kb.get(key)})
    pa, pb = a["pairs_1_4"], b["pairs_1_4"]
    sa, sb = prmtop_sections(new_dir / "system.prmtop"), prmtop_sections(other_dir / "system.prmtop")
    coord = float(np.abs(xyz_new - other.coordinates).max())
    not_allowed = sorted(c for c in classes if c not in allowed)
    res.update({
        "atoms": len(a["atoms"]),
        "atoms_that_differ_in": atom_diff, "atom_types_that_differ": [list(x) for x in type_diff],
        "terms": counts, "differences_by_class": dict(sorted(classes.items())),
        "classes_that_differ_but_may_not": not_allowed,
        "pairs_1_4": {"this": len(pa), "other": len(pb), "identical_with_their_scaling": pa == pb},
        "excluded_atom_lists_identical": a["excluded"] == b["excluded"],
        "largest_coordinate_difference_A": round(coord, 6), "box_identical": bool(np.allclose(parm.box, other.box, atol=1e-6)),
        "system_inpcrd_same_sha256": sha256(new_dir / "system.inpcrd") == sha256(other_dir / "system.inpcrd"),
        "selections_json_same_sha256": sha256(new_dir / "selections.json") == sha256(other_dir / "selections.json"),
        "reference_pdb_same_sha256": sha256(new_dir / "reference.pdb") == sha256(other_dir / "reference.pdb"),
        "system_prmtop_same_sha256": sha256(new_dir / "system.prmtop") == sha256(other_dir / "system.prmtop"),
        "prmtop_sections": {"n": len(sa), "identical_as_text": sum(1 for k in sa if sa[k] == sb.get(k)),
                            "that_differ": sorted(k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k))},
        "terms_that_differ": rows})
    res["only_the_allowed_terms_differ"] = bool(
        not not_allowed and not any(atom_diff.values()) and pa == pb and res["excluded_atom_lists_identical"]
        and coord < 1e-6 and res["box_identical"]
        and all(c["in_one_only"] == 0 for k, c in counts.items() if k in ("bond", "angle"))
        and (counts["improper"]["in_one_only"] == 0 or "improper at N2" in allowed))
    return res


# ------------------------------------------------------------------------------------------------ the build
def folder_name(model, pset=None, variant=None):
    """Name of the model folder under <out-root>, md/topology and the runs: bonded_2e, bonded_2s for a named set,
    bonded_2e_n2fit for a variant of a set (other terms at N2)."""
    name = f"{model}_{pset}" if pset else model
    return f"{name}_{variant}" if variant else name


def default_mcpb_dir(pset):
    """The MCPB.py folder of a named set: the first of MCPB_SET_DIRS that holds the tleap input of step 4."""
    cands = [pathlib.Path(p.format(set=pset)) for p in MCPB_SET_DIRS] if pset else [MD / "qm" / "mcpb"]
    return next((p for p in cands if sorted(p.glob("*_tleap.in"))), cands[0])


def build(model, start, out_root, mcpb_dir=None, mcpb_group=None, records_root=RECORDS, label=None, openmm_check=True,
          buffer_A=BUFFER_A, pset=None, n2_frcmod=None, note=None, variant=None, compare_with=None):
    import parmed
    from scipy.spatial import cKDTree, ConvexHull
    from scipy.spatial.distance import pdist

    inp = load_inputs(start)
    kind, model = model, folder_name(model, pset, variant)   # kind: restrained or bonded; model: the name of the folder
    final = pathlib.Path(out_root) / model / start
    bdir = final / "build"
    out = bdir / "new"                               # the files of the system are written here and moved into place at the
    rec = pathlib.Path(records_root) / model / start  # very end, so a build that stops leaves the system that was there

    # -- everything that can be refused is refused before a folder is made or a file is touched
    mcpb = None
    if kind == "bonded":
        mcpb = read_mcpb(pathlib.Path(mcpb_dir or default_mcpb_dir(pset)).resolve(), mcpb_group)
        if pset and mcpb["methods"] != [SET_METHOD[pset]]:
            raise BuildError(f"--set {pset} needs zinc terms made by the {SET_METHOD[pset]} method, but "
                             f"{mcpb['frcmod']} names {mcpb['methods'] or 'no method'}; give the right --mcpb-dir")
        if n2_frcmod and not pathlib.Path(n2_frcmod).is_file():
            raise BuildError(f"--n2-frcmod: no such file: {n2_frcmod}")
        if pset and pset != "2e" and not n2_frcmod:
            raise BuildError(f"--set {pset} needs --n2-frcmod: the terms at N2 are the same in both parameter sets "
                             "(md/PLAN.md, 6 Oct 2026, 3:02 AM, point 3) and are not taken from this step of MCPB.py. "
                             "Give md/topology/bonded_2e/<start>/n2_terms.frcmod (the terms as built)")
        if variant and not (pset and n2_frcmod):
            raise BuildError("--variant names a set with other terms at N2: it needs --set and --n2-frcmod")
    elif variant or compare_with:
        raise BuildError("--variant and --compare-with belong to --model bonded")
    if compare_with and not (pathlib.Path(out_root) / compare_with / start / "system.prmtop").is_file():
        raise BuildError(f"--compare-with {compare_with}: no built system {pathlib.Path(out_root) / compare_with / start}")
    for d in (bdir, rec):
        d.mkdir(parents=True, exist_ok=True)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir()
    print(f"== {model}/{start}: building in {final}  ({now_ist()})", flush=True)

    # -- model-specific input
    rename, mcpb_record, n2_record, last_frcmods = {}, None, None, []
    n_prot = len(inp["protein"])
    if kind == "bonded":
        units = {u for u, _ in mcpb["mol2"]}
        other_names = []
        if sorted(mcpb["residues"])[:n_prot + 2] != list(range(1, n_prot + 3)):
            raise BuildError(f"{mcpb['pdb'][0][1]}: residues 1..{n_prot + 2} expected (protein as in "
                             "md/system/residue_map.csv, then the zinc, then the ligand)")
        mine = [(k, r["name"], {a["name"]: a["xyz"] for a in r["atoms"]}) for k, r in enumerate(inp["protein"], 1)]
        mine.append((n_prot + 1, "ZN", {"ZN": inp["zinc"][0]["xyz"]}))
        mine.append((n_prot + 2, "SAN", None))
        for k, name, atoms in mine:
            theirs = mcpb["residues"][k]
            if atoms is not None:
                probe = "CA" if "CA" in atoms else next(iter(atoms))
                if probe not in theirs["atoms"] or dist(theirs["atoms"][probe], atoms[probe]) > 0.5:
                    raise BuildError(f"residue {k} of {mcpb['pdb'][0][1]} ({theirs['name']}) is not residue {k} of "
                                     f"the system ({name}): atom {probe} is missing or elsewhere")
            elif not set(LIG_HEAVY) <= set(theirs["atoms"]):
                raise BuildError(f"residue {k} of {mcpb['pdb'][0][1]} is not the ligand")
            if theirs["name"] != name:
                if theirs["name"] in units:
                    rename[k] = theirs["name"]
                else:
                    other_names.append(f"residue {k}: {theirs['name']} in {mcpb['pdb'][0][1]}, {name} here (kept)")
        if n_prot + 1 not in rename or n_prot + 2 not in rename:
            raise BuildError("MCPB.py should have renamed the zinc and the ligand residue")
        used = sorted({int(r) for b in mcpb["bonds"] for r in (b[1], b[4])})
        if max(used) > n_prot + 2:
            raise BuildError("a bond command of MCPB.py names a residue beyond the ligand")
        mcpb_files = [mcpb["tleap_in"], str(pathlib.Path(mcpb["dir"]) / mcpb["pdb"][0][1]),
                      str(pathlib.Path(mcpb["dir"]) / f"{mcpb['group']}_mcpbpy.frcmod")]
        mcpb_files += [str(pathlib.Path(mcpb["dir"]) / f) for _, f in mcpb["mol2"]]
        mcpb_files += [str(pathlib.Path(mcpb["dir"]) / f) for f in mcpb["frcmods"] if (pathlib.Path(mcpb["dir"]) / f).is_file()]
        mcpb_record = {"dir": rel(mcpb["dir"]), "group": mcpb["group"], "atom_types_added": mcpb["atom_types"],
                       "renamed_residues": {str(k): v for k, v in sorted(rename.items())},
                       "bond_commands": [f"bond sys.{b[1]}.{b[2]} sys.{b[4]}.{b[5]}" for b in mcpb["bonds"]],
                       "bond_commands_used": "the four to the zinc; the peptide bonds of the renamed residues are made by "
                                             "loadpdb through head and tail atoms set on their units (see leap.in)",
                       "frcmod_files_named_by_mcpb": mcpb["frcmods"],
                       "residue_names_that_differ_without_a_mol2_unit": other_names,
                       "files_sha256": {rel(p): sha256(p) for p in sorted(set(mcpb_files))}}
        for name in ("sfa_anion_gaff2.frcmod", "sfa_anion_repair.frcmod"):
            p = pathlib.Path(mcpb["dir"]) / name
            mcpb_record[f"{name}_same_as_md_ligand"] = (sha256(p) == sha256(LIGAND / name)) if p.is_file() else None
        mcpb_record["parameter_set"] = {"named": pset, "variant": variant, "method_named_in_the_frcmod": mcpb["methods"],
                                        "expected_method": SET_METHOD.get(pset)}
        # the templates of the five site residues as MCPB.py (and md/qm/resp_equiv.py) wrote them
        mol2_of = dict(mcpb["mol2"])
        site_mol2 = {k: parmed.load_file(str(pathlib.Path(mcpb["dir"]) / mol2_of[rename[k]]), structure=True)
                     for k in sorted(rename)}
        # -- the terms at N2: one file, loaded last, the same for both parameter sets
        lig_types = {a.name: a.type for a in site_mol2[n_prot + 2].atoms}
        n2_types = dict({k: lig_types[k] for k in ("N2", "S", "HN21", "C4")}, ZN=site_mol2[n_prot + 1].atoms[0].type)
        ligx = {a["name"]: a["xyz"] for a in inp["ligand_crystal"]}
        theta_x = angle(inp["zinc"][0]["xyz"], ligx["N2"], ligx["HN21"])
        if n2_frcmod:
            shutil.copy(n2_frcmod, bdir / "n2_terms.frcmod")
            n2_record = {"source": "given with --n2-frcmod", "given_file": rel(n2_frcmod),
                         "first_line_of_the_file": pathlib.Path(n2_frcmod).read_text().splitlines()[0],
                         "validated_by_gate_GR": "not known to the builder: the file was given"}
        elif pset == "2e":
            n2_record = dict({"source": "written by the builder: the terms AS BUILT BEFORE GATE G-R (Zn-N2-S, S-N2-H and "
                                        "H-N2-S-C4 as in the frcmod of MCPB.py; Zn-N2-H at its angle in start X with the "
                                        "force constant of that frcmod)", "validated_by_gate_GR": False},
                             **n2_terms_file(mcpb, n2_types, theta_x, bdir / "n2_terms.frcmod"))
        if n2_record:
            last_frcmods.append(bdir / "n2_terms.frcmod")
            n2_record.update({"file": "n2_terms.frcmod (in build/ of the system and in md/topology)",
                              "sha256": sha256(bdir / "n2_terms.frcmod"), "atom_types": n2_types,
                              "Zn_N2_H_in_start_X_deg": round(theta_x, 4),
                              "loaded": "after every parameter file of MCPB.py, so its terms are the ones in force"})
            # what the file holds, against what a file of this set may hold.  A set without a variant has the terms
            # AS BUILT: the three angles at N2 and the series on H-N2-S-C4, nothing else.  A variant (other terms at
            # N2; md/PLAN.md, 6 Oct 2026, 1:34 PM, point 2) may also hold torsions about N2-S and about Zn-N2 and
            # one improper at N2.
            f = read_frcmod(bdir / "n2_terms.frcmod")
            T = dict(n2_types, O=lig_types["O1"])
            his = sorted({a.type for k in sorted(rename) if k <= n_prot for a in site_mol2[k].atoms
                          if any(int(b[1]) == k and b[2] == a.name or int(b[4]) == k and b[5] == a.name
                                 for b in mcpb["bonds"] if n_prot + 1 in (int(b[1]), int(b[4])))})

            def is_(types, *want):
                return tuple(types) in (tuple(want), tuple(want)[::-1])

            as_built = {"ANGLE": [(T["ZN"], T["N2"], T["S"]), (T["ZN"], T["N2"], T["HN21"]), (T["HN21"], T["N2"], T["S"])],
                        "DIHE": [(T["C4"], T["S"], T["N2"], T["HN21"])]}
            wider = {"DIHE": [(T["ZN"], T["N2"], T["S"], T["C4"]), (T["ZN"], T["N2"], T["S"], T["O"]),
                              (T["HN21"], T["N2"], T["S"], T["O"])]
                     + [(h, T["ZN"], T["N2"], x) for h in his for x in (T["S"], T["HN21"])]}
            outside, beyond_as_built = [], []
            for sec in ("BOND", "ANGLE", "DIHE"):
                for types in f[sec]:
                    if any(is_(types, *w) for w in as_built.get(sec, [])):
                        continue
                    beyond_as_built.append(f"{sec} {'-'.join(types)}")
                    if not (variant and any(is_(types, *w) for w in wider.get(sec, []))):
                        outside.append(f"{sec} {'-'.join(types)}")
            for types in f["IMPROPER"]:
                beyond_as_built.append(f"IMPROPER {'-'.join(types)}")
                if not (variant and types[2] == T["N2"] and sorted(types) == sorted((T["ZN"], T["HN21"], T["N2"], T["S"]))):
                    outside.append(f"IMPROPER {'-'.join(types)}")
            missing = [f"{sec} {'-'.join(w)}" for sec, ws in as_built.items() for w in ws
                       if not any(is_(types, *w) for types in f[sec])]
            n2_record["holds"] = {
                "angles": ["-".join(t) for t in f["ANGLE"]], "torsions": ["-".join(t) for t in f["DIHE"]],
                "impropers": ["-".join(t) for t in f["IMPROPER"]], "bonds": ["-".join(t) for t in f["BOND"]],
                "zinc_nitrogen_types_of_the_histidines": his,
                "class_of_this_set": ("a variant: the three angles at N2, torsions about N2-S and about Zn-N2, one improper "
                                      "at N2" if variant else
                                      "the terms as built: the three angles at N2 and the series on H-N2-S-C4"),
                "entries_beyond_the_terms_as_built": beyond_as_built,
                "entries_outside_the_class_of_this_set": outside, "entries_of_the_class_that_are_missing": missing,
                "of_the_class_of_this_set": not outside and not missing}

    assembly = bdir / "assembly.pdb"
    table = write_assembly(inp, assembly, rename)
    n_xtal = len(inp["waters"])
    n_dry_res = n_prot + 2 + n_xtal

    # -- pass 1: solvated, no ions
    chain_units = [rename[k] for k in sorted(rename) if k <= n_prot]
    head, leap_notes = leap_head(kind, mcpb, assembly, bdir, buffer_A, last_frcmods, chain_units)
    s1 = head + ["charge sys", f"saveamberparm sys {bdir / 'stage1.prmtop'} {bdir / 'stage1.inpcrd'}",
                 f"savepdb sys {bdir / 'stage1.pdb'}", "quit"]
    (bdir / "leap_stage1.in").write_text("\n".join(s1) + "\n")
    log1, count1 = run_tleap(bdir / "leap_stage1.in", bdir)
    (bdir / "leap_stage1.log").write_text(log1)
    p1 = parmed.load_file(str(bdir / "stage1.prmtop"), xyz=str(bdir / "stage1.inpcrd"))
    xyz1 = p1.coordinates.copy()
    if len(p1.residues) < n_dry_res or [r.name for r in p1.residues[n_prot + 2:n_dry_res]] != ["WAT"] * n_xtal:
        raise BuildError("pass 1: the residues are not in the order protein, zinc, ligand, crystal water")
    if any(r.name != "WAT" for r in p1.residues[n_dry_res:]):
        raise BuildError("pass 1: a residue after the crystal water is not water")
    dry = parmed.load_file(str(bdir / "dry.prmtop"))
    q_dry = float(sum(a.charge for a in dry.atoms))
    q_int = int(round(q_dry))
    if abs(q_dry - q_int) > 0.01:
        raise BuildError(f"the charge of the dry system is {q_dry:.4f} e, not an integer")
    n_water1 = len(p1.residues) - (n_prot + 2)
    n_pairs = int(round(SALT_MOL_PER_L * n_water1 / WATER_MOL_PER_L))
    n_na = n_pairs + max(0, -q_int)
    n_cl = n_pairs + max(0, q_int)
    vecs = box_vectors(p1.box)
    volume1 = float(abs(np.linalg.det(vecs)))
    solute_atoms = np.array([a.idx for a in p1.atoms if a.residue.idx < n_dry_res])
    bulk_res = list(p1.residues[n_dry_res:])
    o_bulk = np.array([r.atoms[0].idx for r in bulk_res])
    if any(p1.atoms[i].name != "O" for i in o_bulk):
        raise BuildError("pass 1: the first atom of a water molecule is not O")
    taken, n_cand = choose_ion_sites(xyz1, o_bulk, solute_atoms, n_na + n_cl, vecs)
    # the order taken alternates Na+, Cl-; the ions that neutralise the solute are the last ones taken
    kinds = []
    na_left, cl_left = n_na, n_cl
    for _ in taken:
        if na_left and (na_left >= cl_left):
            kinds.append("Na+")
            na_left -= 1
        else:
            kinds.append("Cl-")
            cl_left -= 1
    ion_sites = [{"ion": kinds[n], "replaces_added_water_number": int(k) + 1,
                  "pass1_residue_number": n_dry_res + int(k) + 1,
                  "xyz_inpcrd_frame_A": [round(float(v), 4) for v in xyz1[o_bulk[k]]]} for n, k in enumerate(taken)]
    ion_sites.sort(key=lambda s: (s["ion"] != "Na+",))          # Na+ first, each kind in the order taken
    # frame of tleap's own coordinates (savepdb) against the inpcrd: a constant shift
    pdb1 = np.array([a["xyz"] for a in read_pdb(bdir / "stage1.pdb")])
    if pdb1.shape != xyz1.shape:
        raise BuildError("pass 1: savepdb and saveamberparm differ in the number of atoms")
    shift = np.median(xyz1 - pdb1, axis=0)
    if np.abs(xyz1 - pdb1 - shift).max() > 0.002:
        raise BuildError("pass 1: the inpcrd is not the PDB frame plus one shift")
    ion_lines = []
    for n, s in enumerate(ion_sites, start=1):
        x = np.array(s["xyz_inpcrd_frame_A"]) - shift
        el = "NA" if s["ion"] == "Na+" else "CL"
        ion_lines += ["HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" %
                      (n, s["ion"], s["ion"], " ", n, x[0], x[1], x[2], el), "TER"]
    (bdir / "ions.pdb").write_text("REMARK   ions of %s/%s in the frame of tleap after solvateOct (Angstrom)\n" % (model, start)
                                   + "\n".join(ion_lines) + "\nEND\n")

    # -- pass 2: the same, water molecules removed, ions added
    s2 = list(head)
    s2 += [f"remove sys sys.{s['pass1_residue_number']}" for s in ion_sites]
    s2 += [f"ions = loadpdb {bdir / 'ions.pdb'}", "sys = combine { sys ions }", "charge sys",
           f"saveamberparm sys {bdir / 'system_tleap.prmtop'} {out / 'system.inpcrd'}", "quit"]
    (bdir / "leap.in").write_text("\n".join(s2) + "\n")
    log2, count2 = run_tleap(bdir / "leap.in", bdir)
    (bdir / "leap.log").write_text(log2)
    first, rest = (bdir / "system_tleap.prmtop").read_text().split("\n", 1)
    stamped = re.sub(r"DATE = \S+\s+\S+", FIXED_PRMTOP_DATE, first)
    if stamped == first or "DATE" not in first:
        raise BuildError("the first line of the prmtop holds no date to replace")
    (out / "system.prmtop").write_text(stamped.ljust(len(first)) + "\n" + rest)
    (bdir / "system_tleap.prmtop").unlink()

    # -- read the result back
    parm = parmed.load_file(str(out / "system.prmtop"), xyz=str(out / "system.inpcrd"))
    xyz = parm.coordinates.copy()
    atoms = parm.atoms
    resn = [r.name for r in parm.residues]
    n_ion = n_na + n_cl
    exp_order = resn[n_prot + 2:n_prot + 2 + n_ion] == ["Na+"] * n_na + ["Cl-"] * n_cl and \
        all(x == "WAT" for x in resn[n_prot + 2 + n_ion:])
    if not exp_order:
        raise BuildError("pass 2: the residues are not in the order protein, zinc, ligand, Na+, Cl-, water")
    removed = {s["pass1_residue_number"] - 1 for s in ion_sites}
    keep1 = np.array([a.idx for a in p1.atoms if a.residue.idx not in removed])
    keep2 = np.array([a.idx for a in atoms if a.residue.name not in ("Na+", "Cl-")])
    if len(keep1) != len(keep2) or [p1.atoms[i].name for i in keep1] != [atoms[i].name for i in keep2]:
        raise BuildError("pass 2 did not remove the chosen water molecules")
    same_as_pass1 = float(np.abs(xyz1[keep1] - xyz[keep2]).max())
    ion_idx = [a.idx for a in atoms if a.residue.name in ("Na+", "Cl-")]
    ion_off = float(max(np.linalg.norm(xyz[i] - np.array(s["xyz_inpcrd_frame_A"])) for i, s in zip(ion_idx, ion_sites)))
    if same_as_pass1 > 1e-6 or ion_off > 0.002 or not np.allclose(parm.box, p1.box):
        raise BuildError(f"pass 2 differs from pass 1 (atoms {same_as_pass1:.2e} A, ions {ion_off:.4f} A)")

    # -- groups
    res_of = np.array([a.residue.idx for a in atoms])
    protein_all = [a.idx for a in atoms if a.residue.idx < n_prot]
    zn_res, lig_res = parm.residues[n_prot], parm.residues[n_prot + 1]
    if len(zn_res.atoms) != 1 or zn_res.atoms[0].atomic_number != 30:
        raise BuildError("the residue after the protein is not the zinc")
    zn = zn_res.atoms[0].idx
    lig = {a.name: a.idx for a in lig_res.atoms}
    if sorted(lig) != sorted(LIG_HEAVY + LIG_H):
        raise BuildError(f"the ligand residue has the atoms {sorted(lig)}")
    first_water = n_prot + 2 + n_ion
    xtal_res = parm.residues[first_water:first_water + n_xtal]
    crystal_o = [r.atoms[0].idx for r in xtal_res]
    water_o = [r.atoms[0].idx for r in parm.residues[first_water:]]
    dep_of_res = {k - 1: inp["resmap"][k]["deposited"] for k in inp["resmap"]}

    def site_atom(dep, name):
        k = next(i for i, d in dep_of_res.items() if d == dep)
        hits = [a.idx for a in parm.residues[k].atoms if a.name == name]
        if len(hits) != 1:
            raise BuildError(f"residue {dep} has {len(hits)} atoms named {name}")
        return hits[0]

    sel = {"note": "zero-based atom indices of system.prmtop; residue numbers in the key names are those of the "
                   "deposited structure 6RL9. Keys beyond the shared conventions: crystal_water_O, "
                   "crystal_water_numbers_deposited, asn67_ND2, ions_Na, ions_Cl.",
           "zn": zn, "lig_heavy": {n: lig[n] for n in LIG_HEAVY}, "lig_HN2": lig["HN21"],
           "lig_all": [a.idx for a in lig_res.atoms]}
    for key, (dep, name) in SITE.items():
        sel[key] = site_atom(dep, name)
    ca = [i for i in protein_all if atoms[i].name == "CA"]
    sel.update({"protein_CA": ca, "protein_heavy": [i for i in protein_all if atoms[i].atomic_number > 1],
                "protein_all": protein_all, "water_O": water_o, "ions": ion_idx,
                "solute": sorted(protein_all + [zn] + sel["lig_all"]),
                "residue_numbers_deposited_for_CA": [dep_of_res[int(res_of[i])] for i in ca],
                "crystal_water_O": crystal_o,
                "crystal_water_numbers_deposited": [a["resi"] for a in inp["waters"]],
                "ions_Na": [i for i in ion_idx if atoms[i].residue.name == "Na+"],
                "ions_Cl": [i for i in ion_idx if atoms[i].residue.name == "Cl-"]})
    write_json(out / "selections.json", sel)

    # -- reference.pdb: crystal coordinates of the CA atoms and of the ligand heavy atoms of SAN 304
    ref = ["REMARK   reference for fitting and RMSD: crystal coordinates of 6RL9 (deposited frame, Angstrom)",
           "REMARK   CA atoms of the protein in the order of selections.json protein_CA (residue names of the",
           "REMARK   built system, deposited residue numbers) and the 11 heavy atoms of sulfanilamide SAN 304",
           "REMARK   (the crystal pose, for every start). Written by md/build_topology.py."]
    ref_ca = []
    for n, i in enumerate(ca, start=1):
        r = inp["protein"][int(res_of[i])]
        a = next(x for x in r["atoms"] if x["name"] == "CA")
        ref_ca.append(a["xyz"])
        ref.append("ATOM  %5d %4s %-4s%1s%4d    %24s  1.00  0.00           C" %
                   (n, a["line"][12:16], atoms[i].residue.name[:4], "A", r["deposited"], a["line"][30:54]))
    ref_lig = {}
    for n, a in enumerate(inp["ligand_crystal"][:len(LIG_HEAVY)], start=len(ca) + 1):
        ref_lig[a["name"]] = a["xyz"]
        ref.append("HETATM%5d %4s %-4s%1s%4d    %24s  1.00  0.00          %2s" %
                   (n, a["line"][12:16], lig_res.name[:4], "A", 304, a["line"][30:54], a["line"][76:78]))
    ref.append("END")
    (out / "reference.pdb").write_text("\n".join(ref) + "\n")
    ref_ca = np.array(ref_ca)

    # -- the built coordinates in the frame of 6RL9 (fit on the CA atoms, as the analysis does)
    R, t = kabsch(xyz[ca], ref_ca)
    fit = xyz @ R.T + t
    ca_fit_rms = float(np.sqrt(((fit[ca] - ref_ca) ** 2).sum(axis=1).mean()))
    dev = []                                         # every heavy atom of the input files against the built system
    for k, r in enumerate(inp["protein"]):
        byname = {a.name: a.idx for a in parm.residues[k].atoms}
        dev += [dist(fit[byname[a["name"]]], a["xyz"]) for a in r["atoms"]]
    dev.append(dist(fit[zn], inp["zinc"][0]["xyz"]))
    dev += [dist(fit[lig[a["name"]]], a["xyz"]) for a in inp["ligand"]]
    dev += [dist(fit[i], a["xyz"]) for i, a in zip(crystal_o, inp["waters"])]
    input_dev = float(max(dev))

    def geometry(X):
        g = {"lig_rmsd_A": {}}
        for name, perm in LABELLINGS:
            g["lig_rmsd_A"][name] = math.sqrt(sum(float(((X[lig[a]] - ref_lig[perm.get(a, a)]) ** 2).sum())
                                                  for a in LIG_HEAVY) / len(LIG_HEAVY))
        g["lig_rmsd_symmetry_corrected_A"] = min(g["lig_rmsd_A"].values())
        P = {n: X[lig[n]] for n in LIG_HEAVY + ["HN21"]}
        zx = X[zn]
        g.update({"d1_N1_Phe131CZ_A": dist(P["N1"], X[sel["phe131_CZ"]]),
                  "d2_N1_Asn67OD1_A": dist(P["N1"], X[sel["asn67_OD1"]]),
                  "N1_Asn67ND2_A": dist(P["N1"], X[sel["asn67_ND2"]]),
                  "Zn_N2_A": dist(zx, P["N2"]), "Zn_His94NE2_A": dist(zx, X[sel["his94_NE2"]]),
                  "Zn_His96NE2_A": dist(zx, X[sel["his96_NE2"]]), "Zn_His119ND1_A": dist(zx, X[sel["his119_ND1"]]),
                  "Zn_O1_A": dist(zx, P["O1"]), "Zn_O2_A": dist(zx, P["O2"]),
                  "N2_Thr199OG1_A": dist(P["N2"], X[sel["thr199_OG1"]]),
                  "HN21_Thr199OG1_A": dist(P["HN21"], X[sel["thr199_OG1"]]),
                  "N2_HN21_Thr199OG1_deg": angle(P["N2"], P["HN21"], X[sel["thr199_OG1"]]),
                  "tau_fold_deg": max(abs(torsion(P["C3"], P["C4"], P["S"], P["N2"])),
                                      abs(torsion(P["C5"], P["C4"], P["S"], P["N2"]))),
                  "C4_S_N2_deg": angle(P["C4"], P["S"], P["N2"]), "S_N2_A": dist(P["S"], P["N2"]),
                  "S_N2_Zn_deg": angle(P["S"], P["N2"], zx), "Zn_N2_S_C4_deg": torsion(zx, P["N2"], P["S"], P["C4"])})
        four = [P["N2"], X[sel["his94_NE2"]], X[sel["his96_NE2"]], X[sel["his119_ND1"]]]
        names = ["N2", "H94NE2", "H96NE2", "H119ND1"]
        g["N_Zn_N_deg"] = {f"{names[i]}-Zn-{names[j]}": angle(four[i], zx, four[j])
                           for i in range(4) for j in range(i + 1, 4)}
        return g

    geo = geometry(fit)
    row = inp["starts"][start]
    compare = {"lig_rmsd_symmetry_corrected_A": "rmsd_to_crystal_A", "d1_N1_Phe131CZ_A": "N1_Phe131CZ",
               "d2_N1_Asn67OD1_A": "N1_Asn67OD1", "N1_Asn67ND2_A": "N1_Asn67ND2", "Zn_N2_A": "Zn_N2",
               "Zn_O1_A": "Zn_O1", "Zn_O2_A": "Zn_O2", "N2_Thr199OG1_A": "N2_OG1", "HN21_Thr199OG1_A": "H_OG1",
               "N2_HN21_Thr199OG1_deg": "N2_H_OG1", "tau_fold_deg": "tau_fold", "C4_S_N2_deg": "C4_S_N2",
               "S_N2_A": "S_N2", "S_N2_Zn_deg": "S_N2_Zn"}
    against = {}
    for mine, theirs in compare.items():
        tol = 0.1 if mine.endswith("_deg") else (0.006 if theirs == "rmsd_to_crystal_A" else 0.002)
        against[mine] = {"built": round(geo[mine], 4), "starts_csv": float(row[theirs]),
                         "difference": round(geo[mine] - float(row[theirs]), 4), "tolerance": tol,
                         "agrees": abs(geo[mine] - float(row[theirs])) <= tol}

    # -- selections: what each index is, and where it is
    lines = [f"selections of {model}/{start}: key, index (zero-based), residue name, deposited number, atom name, "
             f"distance (A) from the same atom of the md/system files after the fit on CA"]
    sel_ok = True

    def describe(key, i, dep, exp_name, where):
        nonlocal sel_ok
        a = atoms[i]
        d = None if where is None else dist(fit[i], where)
        ok = a.name == exp_name and (d is None or d < 0.01)
        sel_ok = sel_ok and ok
        lines.append(f"{key:14s} {i:6d}  {a.residue.name:4s} {dep:4d}  {a.name:5s} "
                     f"{'' if d is None else '%.4f' % d}  {'ok' if ok else 'WRONG'}")

    describe("zn", zn, 301, "ZN", inp["zinc"][0]["xyz"])
    ligxyz = {a["name"]: a["xyz"] for a in inp["ligand"]}
    for n in LIG_HEAVY:
        describe(f"lig_heavy.{n}", lig[n], 304, n, ligxyz[n])
    describe("lig_HN2", sel["lig_HN2"], 304, "HN21", ligxyz["HN21"])
    for key, (dep, name) in SITE.items():
        k = next(i for i, d in dep_of_res.items() if d == dep)
        src = next((a["xyz"] for a in inp["protein"][k]["atoms"] if a["name"] == name), None)
        describe(key, sel[key], dep, name, src)
        if src is None:                              # a hydrogen: it must sit on the right heavy atom
            partner = {"HG1": "OG1", "H": "N"}[name]
            bonded = [b.name for b in atoms[sel[key]].bond_partners]
            ok = bonded == [partner]
            sel_ok = sel_ok and ok
            lines[-1] += f"  (hydrogen; bonded to {','.join(bonded)})"
    hn_partner = [b.name for b in atoms[sel["lig_HN2"]].bond_partners]
    sel_ok = sel_ok and hn_partner == ["N2"]
    names_ok = (all(atoms[i].name == "CA" for i in ca) and len(ca) == n_prot - 1
                and all(atoms[i].name == "O" and atoms[i].residue.name == "WAT" for i in water_o)
                and sel["solute"] == list(range(len(sel["solute"])))
                and sel["residue_numbers_deposited_for_CA"] == [r["deposited"] for r in inp["protein"][1:]])
    lines.append(f"protein_CA: {len(ca)} atoms, all named CA: {all(atoms[i].name == 'CA' for i in ca)}; deposited "
                 f"numbers {sel['residue_numbers_deposited_for_CA'][0]}..{sel['residue_numbers_deposited_for_CA'][-1]}")
    lines.append(f"protein_heavy {len(sel['protein_heavy'])}, protein_all {len(protein_all)}, lig_all {len(sel['lig_all'])}, "
                 f"solute {len(sel['solute'])} (indices 0..{sel['solute'][-1]}), water_O {len(water_o)}, "
                 f"crystal_water_O {len(crystal_o)}, ions {len(ion_idx)} (Na+ {n_na}, Cl- {n_cl})")
    (rec / "selections_table.txt").write_text("\n".join(lines) + "\n")

    # -- force-field read-back
    lig_atoms = list(lig_res.atoms)
    terms = ligand_terms(parm, lig_atoms)
    refparm = parmed.load_file(str(inp["files"]["ligand_repaired_prmtop"]))
    ref_terms = ligand_terms(refparm, list(refparm.atoms))
    repair = read_frcmod(inp["files"]["frcmod_repair"])
    assigned = read_frcmod(inp["files"]["frcmod_gaff2"])
    mol2_types = {a.name: a.type for a in parmed.load_file(str(inp["files"]["mol2"]), structure=True).atoms}
    mol2_charge = {a.name: a.charge for a in parmed.load_file(str(inp["files"]["mol2"]), structure=True).atoms}

    def by_types(table, names):
        tt = tuple(mol2_types[n] for n in names)
        return table.get(tt, table.get(tt[::-1]))

    repaired_check = []
    for names in (("N2", "S"),):
        exp = by_types(repair["BOND"], names)
        got = terms["bonds"].get("-".join(names if names <= names[::-1] else names[::-1]))
        repaired_check.append({"term": "bond " + "-".join(names), "types_free_anion": "-".join(mol2_types[n] for n in names),
                               "repair_frcmod": list(exp), "in_prmtop": got,
                               "in_force": got is not None and abs(got[0] - exp[0]) < 5e-3 and abs(got[1] - exp[1]) < 5e-5})
    for names in (("C4", "S", "N2"), ("HN21", "N2", "S")):
        exp = by_types(repair["ANGLE"], names)
        got = terms["angles"].get("-".join(names if names <= names[::-1] else names[::-1]))
        repaired_check.append({"term": "angle " + "-".join(names), "types_free_anion": "-".join(mol2_types[n] for n in names),
                               "repair_frcmod": list(exp), "in_prmtop": got,
                               "in_force": got is not None and abs(got[0] - exp[0]) < 5e-3 and abs(got[1] - exp[1]) < 5e-3})
    for names in (("HN21", "N2", "S", "C4"), ("HN21", "N2", "S", "O1"), ("HN21", "N2", "S", "O2")):
        exp = sorted([round(k, 5), n, round(ph, 3) % 360.0] for k, ph, n in by_types(repair["DIHE"], names))
        old = sorted([round(k, 5), n, round(ph, 3) % 360.0] for k, ph, n in by_types(assigned["DIHE"], names))
        got = sorted(tm[:3] for tm in terms["dihedrals"].get("-".join(names if names <= names[::-1] else names[::-1]), []))
        nz = lambda L: [x for x in L if abs(x[0]) > 1e-9]
        repaired_check.append({"term": "torsion " + "-".join(names),
                               "types_free_anion": "-".join(mol2_types[n] for n in names),
                               "repair_frcmod_k_per_phase": exp, "gaff2_as_assigned_k_per_phase": old,
                               "in_prmtop_k_per_phase": got,
                               "in_force": nz(got) == nz(exp) and len(got) == len(exp)})
    same_as_free = {sec: terms[sec] == ref_terms[sec] for sec in terms}
    lig_charge = {a.name: round(a.charge, 6) for a in lig_atoms}
    charges_as_mol2 = all(abs(lig_charge[n] - mol2_charge[n]) < 2e-6 for n in mol2_charge)
    znatom = atoms[zn]
    zn_rec = {"atom_type": znatom.type, "charge_e": round(znatom.charge, 6), "rmin_half_A": round(znatom.rmin, 5),
              "epsilon_kcal_per_mol": round(znatom.epsilon, 8), "mass_amu": round(znatom.mass, 3),
              "bonded_to": sorted(f"{b.residue.name}{dep_of_res.get(b.residue.idx, 304)}:{b.name}" for b in znatom.bond_partners),
              "lennard_jones_12_6_4_terms_in_prmtop": "LENNARD_JONES_CCOEF" in parm.parm_data}
    ion_rec = {}
    for nm, lst in (("Na+", sel["ions_Na"]), ("Cl-", sel["ions_Cl"])):
        if lst:
            a = atoms[lst[0]]
            ion_rec[nm] = {"n": len(lst), "atom_type": a.type, "charge_e": round(a.charge, 6),
                           "rmin_half_A": round(a.rmin, 5), "epsilon_kcal_per_mol": round(a.epsilon, 8)}
    parmdir = pathlib.Path(AMBERHOME) / "dat" / "leap" / "parm"
    ffiles = {"frcmod.ions234lm_126_tip3p": frcmod_nonbon(parmdir / "frcmod.ions234lm_126_tip3p", "Zn2+"),
              "frcmod.ionsjc_tip3p Na+": frcmod_nonbon(parmdir / "frcmod.ionsjc_tip3p", "Na+"),
              "frcmod.ionsjc_tip3p Cl-": frcmod_nonbon(parmdir / "frcmod.ionsjc_tip3p", "Cl-")}
    ions_as_jc = all(abs(ion_rec[k]["rmin_half_A"] - ffiles[f"frcmod.ionsjc_tip3p {k}"]["rmin_half_A"]) < 1e-4 and
                     abs(ion_rec[k]["epsilon_kcal_per_mol"] - ffiles[f"frcmod.ionsjc_tip3p {k}"]["epsilon_kcal_per_mol"]) < 1e-7
                     for k in ion_rec)

    # -- extra forces
    r0 = {}
    for name, dep, text in ZN_PARTNERS:
        if dep is None:
            j, x = lig[name], next(a["xyz"] for a in inp["ligand_crystal"] if a["name"] == name)
        else:
            j = site_atom(dep, name)
            k = next(i for i, d in dep_of_res.items() if d == dep)
            x = next(a["xyz"] for a in inp["protein"][k]["atoms"] if a["name"] == name)
        r0[text] = (j, dist(inp["zinc"][0]["xyz"], x))
    sigma_A = 10.0 * math.sqrt(KT_300_KJ / ZN_RESTRAINT_K_KJ_NM2)
    if kind == "restrained":
        bonds = [{"i": zn, "j": j, "r0_nm": round(d / 10.0, 5), "k_kj_per_mol_nm2": ZN_RESTRAINT_K_KJ_NM2,
                  "label": f"Zn to {text}; r0 = {d:.3f} A, the distance in 6RL9"} for text, (j, d) in r0.items()]
        extra = {"energy": "U = 0.5 k (r - r0)^2 for each entry; r0 in nm, k in kJ mol-1 nm-2",
                 "model": "restrained: nonbonded Zn2+ with four harmonic distance restraints; nothing acts on "
                          "angles, torsions or the orientation of the ring",
                 "expected_thermal_spread_A_at_300K": round(sigma_A, 4), "harmonic_bonds": bonds}
    else:
        extra = {"energy": "U = 0.5 k (r - r0)^2 for each entry; r0 in nm, k in kJ mol-1 nm-2",
                 "model": "bonded: the zinc-nitrogen bonds are in system.prmtop; no extra force", "harmonic_bonds": []}
    write_json(out / "extra_forces.json", extra)

    # -- how fast are the restraints, against the bonds the 4 fs step already carries?
    # masses as the driver sets them: every hydrogen outside water has hydrogen_mass_amu, taken from its heavy atom
    h_mass = json.loads((MD / "protocol.json").read_text())["system"]["hydrogen_mass_amu"]
    hmr = np.array([a.mass for a in atoms], dtype=float)
    for a in atoms:
        if a.atomic_number == 1 and a.residue.name != "WAT" and len(a.bond_partners) == 1:
            hmr[a.bond_partners[0].idx] -= h_mass - a.mass
            hmr[a.idx] = h_mass

    def period_fs(k_kj_nm2, i, j):                   # harmonic period of two masses on a spring, fs
        mu = hmr[i] * hmr[j] / (hmr[i] + hmr[j])
        return 2000.0 * math.pi * math.sqrt(mu / k_kj_nm2)

    fastest = None
    for b in parm.bonds:
        if b.atom1.atomic_number > 1 and b.atom2.atomic_number > 1 and b.atom1.residue.name != "WAT":
            t_fs = period_fs(2.0 * b.type.k * 418.4, b.atom1.idx, b.atom2.idx)      # Amber: U = K (r - r0)^2
            if fastest is None or t_fs < fastest[0]:
                fastest = (t_fs, b)
    stiffness = {"fastest_bond_of_the_topology_without_hydrogen": {
        "atoms": f"{fastest[1].atom1.residue.name}:{fastest[1].atom1.name}-{fastest[1].atom2.name}",
        "types": f"{fastest[1].atom1.type}-{fastest[1].atom2.type}",
        "k_kj_per_mol_nm2_in_the_form_half_k_x2": round(2.0 * fastest[1].type.k * 418.4, 0),
        "masses_amu": [round(float(hmr[fastest[1].atom1.idx]), 3), round(float(hmr[fastest[1].atom2.idx]), 3)],
        "period_fs": round(fastest[0], 1), "steps_of_4_fs_per_period": round(fastest[0] / 4.0, 1)},
        "masses": f"as the driver sets them (hydrogen {h_mass} amu, taken from the heavy atom); bonds to hydrogen are "
                  "constrained and not counted",
        "restraints": {text: {"masses_amu": [round(float(hmr[zn]), 3), round(float(hmr[j]), 3)],
                              "period_fs": round(period_fs(ZN_RESTRAINT_K_KJ_NM2, zn, j), 1),
                              "steps_of_4_fs_per_period": round(period_fs(ZN_RESTRAINT_K_KJ_NM2, zn, j) / 4.0, 1)}
                       for text, (j, d) in r0.items()},
        "note": "period of two free masses on the spring alone; the nitrogen atoms sit in a ring or in the anion, and "
                "the nonbonded attraction to the zinc adds stiffness, so the real motion differs somewhat"}

    # -- zinc site of the bonded model, read back
    bonded_rec = None
    if kind == "bonded":
        want = sorted([f"{atoms[j].residue.name}{dep_of_res.get(atoms[j].residue.idx, 304)}:{atoms[j].name}"
                       for j, _ in r0.values()])
        zb = [{"to": f"{(b.atom2 if b.atom1.idx == zn else b.atom1).residue.name}:{(b.atom2 if b.atom1.idx == zn else b.atom1).name}",
               "k_kcal_per_mol_A2": round(b.type.k, 3), "r0_A": round(b.type.req, 4)}
              for b in parm.bonds if zn in (b.atom1.idx, b.atom2.idx)]
        za = [{"atoms": f"{a.atom1.residue.name}:{a.atom1.name}-{a.atom2.name}-{a.atom3.residue.name}:{a.atom3.name}",
               "k_kcal_per_mol_rad2": round(a.type.k, 3), "theta0_deg": round(a.type.theteq, 3)}
              for a in parm.angles if zn in (a.atom1.idx, a.atom2.idx, a.atom3.idx)]
        zd = [abs(d.type.phi_k) for d in parm.dihedrals if zn in (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)]
        site_res = sorted({int(k) - 1 for k in rename})
        bonded_rec = {"zinc_bonded_to": zn_rec["bonded_to"], "expected": want,
                      "zinc_has_the_four_bonds": zn_rec["bonded_to"] == want,
                      "bonds_to_zinc": zb, "angles_with_zinc": za,
                      "n_angles_centred_on_zinc": sum(1 for a in parm.angles if a.atom2.idx == zn),
                      "torsions_through_zinc": {"n": len(zd), "largest_abs_k_kcal_per_mol": round(max(zd), 6) if zd else None},
                      "charge_of_zinc_site_residues_e": {parm.residues[k].name: round(sum(a.charge for a in parm.residues[k].atoms), 5)
                                                           for k in site_res},
                      "charge_of_zinc_site_e": round(sum(a.charge for k in site_res for a in parm.residues[k].atoms), 5),
                      "ligand_atom_types": {a.name: a.type for a in lig_atoms}}

        # ---- every zinc term against the parameter files: the terms-at-N2 file (loaded last), then MCPB.py's frcmod
        frc = read_frcmod(mcpb["frcmod"])
        n2f = read_frcmod(bdir / "n2_terms.frcmod") if n2_record else {"BOND": {}, "ANGLE": {}, "DIHE": {}, "IMPROPER": {}}

        def look(section, aa):
            tt = tuple(a.type for a in aa)
            for table, src in ((n2f[section], "n2_terms.frcmod"), (frc[section], mcpb["frcmod"].name)):
                hit = table.get(tt, table.get(tt[::-1]))
                if hit is not None:
                    return hit, src
            return None, None

        def tag(a):
            return f"{a.residue.name}{dep_of_res.get(a.residue.idx, 301 if a.idx == zn else 304)}:{a.name}"

        crystal_d = {j: d for j, d in r0.values()}
        zbonds, zangles = [], []
        for b in parm.bonds:
            if zn in (b.atom1.idx, b.atom2.idx):
                o = b.atom2 if b.atom1.idx == zn else b.atom1
                exp, src = look("BOND", (b.atom1, b.atom2))
                zbonds.append({"bond": f"Zn-{tag(o)}", "types": f"{b.atom1.type}-{b.atom2.type}",
                               "k_kcal_per_mol_A2": round(b.type.k, 3), "r0_A": round(b.type.req, 4),
                               "parameter_file": list(exp) if exp else None, "from": src,
                               "as_file": exp is not None and abs(b.type.k - exp[0]) < 5e-3 and abs(b.type.req - exp[1]) < 5e-5,
                               "distance_in_6RL9_A": round(crystal_d[o.idx], 4),
                               "distance_in_this_start_A": round(dist(xyz[zn], xyz[o.idx]), 4)})
        for a in parm.angles:
            trio = (a.atom1, a.atom2, a.atom3)
            if zn not in [x.idx for x in trio]:
                continue
            if a.atom2.idx == zn:
                group = "N-Zn-N"
            else:
                end = a.atom3 if a.atom1.idx == zn else a.atom1
                group = ("Zn-N2-S" if end.name == "S" else "Zn-N2-H") if end.residue.idx == lig_res.idx else "Zn-N-C"
            exp, src = look("ANGLE", trio)
            zangles.append({"angle": "-".join("Zn" if x.idx == zn else tag(x) for x in trio),
                            "kind": group, "types": "-".join(x.type for x in trio),
                            "k_kcal_per_mol_rad2": round(a.type.k, 3), "theta0_deg": round(a.type.theteq, 3),
                            "parameter_file": list(exp) if exp else None, "from": src,
                            "as_file": exp is not None and abs(a.type.k - exp[0]) < 5e-3 and abs(a.type.theteq - exp[1]) < 5e-3,
                            "angle_in_this_start_deg": round(angle(*(xyz[x.idx] for x in trio)), 2)})
        zangles.sort(key=lambda z: (["N-Zn-N", "Zn-N-C", "Zn-N2-S", "Zn-N2-H"].index(z["kind"]), z["angle"]))
        ztor = [d for d in parm.dihedrals if zn in (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)]

        # ---- the terms at N2 as they stand in the topology
        n2a, sa, ha, c4a = (atoms[lig[n]] for n in ("N2", "S", "HN21", "C4"))

        def angle_term(i, j, k):
            hits = [a for a in parm.angles if a.atom2.idx == j.idx and {a.atom1.idx, a.atom3.idx} == {i.idx, k.idx}]
            return hits[0] if len(hits) == 1 else None

        def series(table_entry):
            return sorted([round(k, 5), n, round(ph, 3) % 360.0] for k, ph, n in table_entry if abs(k) > 1e-9)

        n2_angles = {}
        for label_, trio in (("Zn-N2-S", (znatom, n2a, sa)), ("Zn-N2-H", (znatom, n2a, ha)), ("S-N2-H", (sa, n2a, ha))):
            a = angle_term(*trio)
            exp = n2f["ANGLE"].get(tuple(x.type for x in trio), n2f["ANGLE"].get(tuple(x.type for x in trio)[::-1]))
            n2_angles[label_] = {"types": "-".join(x.type for x in trio),
                                 "k_kcal_per_mol_rad2": round(a.type.k, 3) if a else None,
                                 "theta0_deg": round(a.type.theteq, 3) if a else None,
                                 "n2_terms_file": list(exp) if exp else None,
                                 "as_n2_terms_file": (a is not None and exp is not None and abs(a.type.k - exp[0]) < 5e-3
                                                      and abs(a.type.theteq - exp[1]) < 5e-3) if n2_record else None,
                                 "angle_in_this_start_deg": round(angle(*(xyz[x.idx] for x in trio)), 2)}
        tor_got = sorted([round(d.type.phi_k, 5), int(d.type.per), round(d.type.phase, 3) % 360.0] for d in parm.dihedrals
                         if not d.improper and abs(d.type.phi_k) > 1e-9 and
                         [d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx] in ([ha.idx, n2a.idx, sa.idx, c4a.idx],
                                                                                  [c4a.idx, sa.idx, n2a.idx, ha.idx]))
        tt = (c4a.type, sa.type, n2a.type, ha.type)
        tor_exp = n2f["DIHE"].get(tt, n2f["DIHE"].get(tt[::-1]))
        # every entry of the terms-at-N2 file against the topology (angles, torsion series, impropers), and the
        # torsions through the zinc quadruple by quadruple
        def typed(aa):
            return tuple(a.type for a in aa)

        def same_types(tt_, want):
            return tt_ == tuple(want) or tt_ == tuple(want)[::-1]

        quads = {}                                   # proper torsions: atoms of the quadruple -> its terms
        for d in parm.dihedrals:
            if not d.improper:
                aa = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
                quads.setdefault(aa if aa[0] < aa[3] else aa[::-1], []).append(
                    [round(d.type.phi_k, 5), int(d.type.per), round(d.type.phase, 3) % 360.0])
        file_entries = []
        if n2_record:
            for types, (k_, th_) in n2f["ANGLE"].items():
                hits = [a for a in parm.angles if same_types(typed((a.atom1, a.atom2, a.atom3)), types)]
                file_entries.append({"entry": "ANGLE " + "-".join(types), "in_the_file": [k_, th_],
                                     "n_in_the_topology": len(hits),
                                     "as_the_file": bool(hits) and all(abs(a.type.k - k_) < 5e-3 and
                                                                       abs(a.type.theteq - th_) < 5e-3 for a in hits)})
            for types, terms_ in n2f["DIHE"].items():
                hits = [q for q in quads if same_types(tuple(atoms[i].type for i in q), types)]
                file_entries.append({"entry": "DIHE " + "-".join(types), "in_the_file_k_n_phase": series(terms_),
                                     "n_in_the_topology": len(hits),
                                     "as_the_file": bool(hits) and all(sorted(t for t in quads[q] if abs(t[0]) > 1e-9)
                                                                       == series(terms_) for q in hits)})
            for types, (k_, ph_, n_) in n2f["IMPROPER"].items():
                hits = [d for d in parm.dihedrals if d.improper and d.atom3.type == types[2]
                        and sorted(typed((d.atom1, d.atom2, d.atom3, d.atom4))) == sorted(types)]
                file_entries.append({"entry": "IMPROPER " + "-".join(types), "in_the_file_k_phase_n": [k_, ph_, n_],
                                     "n_in_the_topology": len(hits),
                                     "as_the_file": bool(hits) and all(abs(d.type.phi_k - k_) < 5e-4 and int(d.type.per) == n_
                                                                       and abs(d.type.phase % 360.0 - ph_ % 360.0) < 5e-3
                                                                       for d in hits)})
        zquads = {q: t for q, t in quads.items() if zn in q}
        zrows, n_carried, tz_all_zero, tz_zero_or_file = [], 0, bool(zquads), bool(zquads)
        for q, terms_ in sorted(zquads.items()):
            tt_ = tuple(atoms[i].type for i in q)
            exp_ = n2f["DIHE"].get(tt_, n2f["DIHE"].get(tt_[::-1]))
            nonzero = sorted(t for t in terms_ if abs(t[0]) > 1e-9)
            ok_ = (nonzero == series(exp_)) if exp_ is not None else (not nonzero)
            n_carried += exp_ is not None
            tz_all_zero = tz_all_zero and not nonzero
            tz_zero_or_file = tz_zero_or_file and ok_
            if nonzero or exp_ is not None:
                zrows.append({"atoms": "-".join("Zn" if i == zn else tag(atoms[i]) for i in q), "types": "-".join(tt_),
                              "terms_k_n_phase": nonzero, "carried_by_the_n2_file": exp_ is not None,
                              "zero_or_as_the_n2_file": bool(ok_)})
        zimp = [d for d in ztor if d.improper]
        zimp_rows = [{"atoms": "-".join("Zn" if x.idx == zn else tag(x) for x in (d.atom1, d.atom2, d.atom3, d.atom4)),
                      "central_atom": tag(d.atom3), "types": "-".join(typed((d.atom1, d.atom2, d.atom3, d.atom4))),
                      "k_phase_n": [round(d.type.phi_k, 5), round(d.type.phase, 3) % 360.0, int(d.type.per)]} for d in zimp]
        zimp_ok = len(zimp) == len(n2f["IMPROPER"]) and all(
            any(d.atom3.type == types[2] and sorted(typed((d.atom1, d.atom2, d.atom3, d.atom4))) == sorted(types)
                and abs(d.type.phi_k - v[0]) < 5e-4 for d in zimp) for types, v in n2f["IMPROPER"].items())
        if n2_record:
            n2_record["in_the_topology"] = {
                "angles": n2_angles,
                "sum_of_the_three_equilibrium_angles_deg": round(sum(v["theta0_deg"] or 0.0 for v in n2_angles.values()), 2),
                "H-N2-S-C4_k_n_phase": tor_got, "H-N2-S-C4_n2_terms_file": series(tor_exp) if tor_exp else None,
                "entries_of_the_file": file_entries,
                "all_as_n2_terms_file": bool(all(v["as_n2_terms_file"] for v in n2_angles.values())
                                             and tor_exp is not None and tor_got == series(tor_exp)
                                             and all(e["as_the_file"] for e in file_entries)),
                "Zn-N2-H_is_not_120_deg": abs((n2_angles["Zn-N2-H"]["theta0_deg"] or 120.0) - 120.0) > 1.0}
            for c in repaired_check:                 # a repaired ligand term that the terms-at-N2 file carries
                what, names = c["term"].split(" ", 1)
                sec = {"bond": "BOND", "angle": "ANGLE", "torsion": "DIHE"}[what]
                types_now = tuple(atoms[lig[n]].type for n in names.split("-"))
                hit = n2f[sec].get(types_now, n2f[sec].get(types_now[::-1]))
                if hit is None:
                    continue
                c["types_in_this_system"] = "-".join(types_now)
                c["expected_from"] = "n2_terms.frcmod, which is loaded after the repair frcmod and MCPB.py's frcmod"
                if sec == "DIHE":
                    c["n2_terms_file_k_per_phase"] = series(hit)
                    c["in_force"] = [x for x in c["in_prmtop_k_per_phase"] if abs(x[0]) > 1e-9] == series(hit)
                else:
                    c["n2_terms_file"] = list(hit)
                    got = c["in_prmtop"]
                    c["in_force"] = got is not None and abs(got[0] - hit[0]) < 5e-3 and abs(got[1] - hit[1]) < 5e-3

        # ---- charges of the five residues of the site against the templates (mol2 files) they came from
        gate_copy = MD / "qm" / "gr" / "charges"
        site_q, q_ok, pair_rows = {}, True, []
        for k in sorted(rename):
            res = parm.residues[k - 1]
            fname = mol2_of[rename[k]]
            tmpl = {a.name: a.charge for a in site_mol2[k].atoms}
            got = {a.name: a.charge for a in res.atoms}
            worst = max(abs(got[n] - tmpl[n]) for n in tmpl) if set(got) == set(tmpl) else None
            q_ok = q_ok and worst is not None and worst < 2e-6
            row = {"residue": f"{res.name} (number {dep_of_res.get(k - 1, 301 if k == n_prot + 1 else 304)} of 6RL9)",
                   "template": fname,
                   "template_charge_line": (pathlib.Path(mcpb["dir"]) / fname).read_text().splitlines()[4].strip(),
                   "sum_e": round(sum(got.values()), 6), "largest_difference_from_template_e": worst,
                   "atoms_e": {a.name: round(a.charge, 6) for a in res.atoms}}
            free = pathlib.Path(mcpb["dir"]) / "charges_unconstrained" / fname
            if free.is_file():
                fq = {a.name: a.charge for a in parmed.load_file(str(free), structure=True).atoms}
                row["largest_difference_from_the_unconstrained_fit_e"] = round(max(abs(got[n] - fq[n]) for n in fq), 6)
            if (gate_copy / fname).is_file():
                row["template_same_as_md_qm_gr_charges"] = sha256(gate_copy / fname) == sha256(pathlib.Path(mcpb["dir"]) / fname)
            site_q[res.name] = row
        ligq = site_q[lig_res.name]["atoms_e"]
        pair_rows = [{"pair": f"{a}, {b}", "charges_e": [ligq[a], ligq[b]], "equal": abs(ligq[a] - ligq[b]) < 1e-6}
                     for a, b in LIG_EQUAL_PAIRS]

        # ---- exclusions and 1-4 pairs across the zinc, as the topology has them
        ps = zinc_pair_sets(parm, zn)
        near = {i for d in ps for pr in ps[d] for i in pr}
        excl, ptr = set(), 0
        for i, n in enumerate(parm.parm_data["NUMBER_EXCLUDED_ATOMS"]):
            if i in near:
                excl |= {(i, j - 1) for j in parm.parm_data["EXCLUDED_ATOMS_LIST"][ptr:ptr + n] if j > 0}
            ptr += n
        one4 = {}
        for d in parm.dihedrals:
            if not d.improper and not d.ignore_end:
                key = tuple(sorted((d.atom1.idx, d.atom4.idx)))
                if key[0] in near and key[1] in near:
                    one4.setdefault(key, []).append((d.type.scee, d.type.scnb))
        rule12 = all(p in excl and p not in one4 for p in ps[1] + ps[2])
        rule3 = all(p in excl and len(one4.get(p, [])) == 1 and abs(one4[p][0][0] - SCEE) < 1e-6
                    and abs(one4[p][0][1] - SCNB) < 1e-6 for p in ps[3])
        rule4 = not any(p in excl or p in one4 for p in ps[4])

        def pair_name(p):
            return " ... ".join("Zn" if i == zn else tag(atoms[i]) for i in p)

        bonded_rec.update({
            "bonds_to_zinc_against_the_parameter_files": zbonds,
            "angles_with_zinc_against_the_parameter_files": zangles,
            "n_angles_by_kind": {g: sum(1 for z in zangles if z["kind"] == g) for g in ("N-Zn-N", "Zn-N-C", "Zn-N2-S", "Zn-N2-H")},
            "all_zinc_bonds_and_angles_as_the_parameter_files": bool(zbonds and zangles and all(z["as_file"] for z in zbonds + zangles)),
            "torsions_through_zinc": {"n": len(zquads), "n_terms": sum(len(t) for t in zquads.values()),
                                      "n_proper": len(zquads), "n_improper": len(zimp),
                                      "largest_abs_k_kcal_per_mol": round(max(abs(d.type.phi_k) for d in ztor), 6) if ztor else None,
                                      "all_zero": bool(tz_all_zero and not zimp),
                                      "n_carried_by_the_n2_file": int(n_carried),
                                      "zero_or_as_the_n2_file": bool(tz_zero_or_file),
                                      "non_zero_or_carried_by_the_n2_file": zrows,
                                      "impropers_with_zinc": zimp_rows,
                                      "impropers_with_zinc_none_or_those_of_the_n2_file": bool(zimp_ok)},
            "charges": {"residues": site_q, "site_sum_e": bonded_rec["charge_of_zinc_site_e"],
                        "zinc_e": round(znatom.charge, 6), "all_as_the_templates": bool(q_ok),
                        "equal_pairs_of_the_ligand": pair_rows, "pairs_equal": all(r_["equal"] for r_ in pair_rows),
                        "resp_equiv_json_sha256": (sha256(pathlib.Path(mcpb["dir"]) / "resp_equiv.json")
                                                   if (pathlib.Path(mcpb["dir"]) / "resp_equiv.json").is_file() else None)},
            "pairs_across_the_zinc": {
                "rule": "Amber: atoms one or two bonds apart do not interact (excluded); atoms three bonds apart interact "
                        f"with the electrostatics divided by {SCEE} and the Lennard-Jones term by {SCNB} (1-4); atoms four "
                        "or more bonds apart interact in full. Checked here for every pair whose bonds run through the zinc",
                "n_pairs_1_or_2_bonds_apart": len(ps[1]) + len(ps[2]), "n_pairs_3_bonds_apart": len(ps[3]),
                "n_pairs_4_bonds_apart": len(ps[4]),
                "1_or_2_bonds_excluded_and_no_1_4_term": bool(rule12),
                "3_bonds_one_1_4_term_each_with_the_amber_scaling": bool(rule3),
                "4_bonds_neither_excluded_nor_1_4": bool(rule4),
                "pairs_1_or_2_bonds_apart": [pair_name(p) for p in ps[1] + ps[2]],
                "pairs_3_bonds_apart_1_4": [pair_name(p) for p in ps[3]],
                "zinc_to_sulfonyl_oxygens": "three bonds apart (Zn-N2-S-O): a 1-4 pair, scaled as above"},
        })

    # -- water at the ligand, and in the place the crystal ligand takes (empty in start R)
    wat_fit = fit[water_o]
    is_xtal = np.arange(len(water_o)) < n_xtal
    lig_now = fit[[lig[n] for n in LIG_HEAVY]]
    lig_xtal = np.array([ref_lig[n] for n in LIG_HEAVY])
    d_now = cKDTree(lig_now).query(wat_fit, k=1)[0]
    d_xtal = cKDTree(lig_xtal).query(wat_fit, k=1)[0]
    site_water = {
        "crystal_water_O_within_3.5_A_of_a_ligand_heavy_atom": int((d_now[is_xtal] <= 3.5).sum()),
        "added_water_O_within_3.5_A_of_a_ligand_heavy_atom": int((d_now[~is_xtal] <= 3.5).sum()),
        "crystal_water_O_within_2.5_A_of_a_position_of_the_crystal_ligand": int((d_xtal[is_xtal] <= 2.5).sum()),
        "added_water_O_within_2.5_A_of_a_position_of_the_crystal_ligand": int((d_xtal[~is_xtal] <= 2.5).sum()),
        "nearest_added_water_O_to_a_ligand_heavy_atom_A": round(float(d_now[~is_xtal].min()), 2),
        "nearest_added_water_O_to_the_zinc_A": round(float(np.linalg.norm(wat_fit[~is_xtal] - fit[zn], axis=1).min()), 2),
        "nearest_atom_to_each_position_of_the_crystal_ligand_A": {
            n: round(float(d), 2) for n, d in zip(LIG_HEAVY, cKDTree(fit).query(lig_xtal, k=1)[0])},
        "positions_of_the_crystal_ligand_with_no_atom_within_2.0_A": [
            n for n, d in zip(LIG_HEAVY, cKDTree(fit).query(lig_xtal, k=1)[0]) if d > 2.0],
        "note": "positions of the crystal ligand: the 11 heavy atoms of SAN 304; a water molecule there in start R "
                "means that solvateOct filled the place the rotated ligand leaves empty"}

    # -- box, clearance, density
    vecs = box_vectors(parm.box)
    volume = float(abs(np.linalg.det(vecs)))
    mass = float(sum(a.mass for a in atoms))
    sol = np.array(sel["solute"])
    tree = cKDTree(xyz[sol])
    image_min = min(float(tree.query(xyz[sol] + s, k=1)[0].min()) for s in image_shifts(vecs))
    hull = xyz[sol][ConvexHull(xyz[sol]).vertices]
    extent = float(pdist(hull).max())
    near = cKDTree(xyz[sol]).query(xyz[ion_idx], k=1)[0] if ion_idx else np.array([])
    shifts = [np.zeros(3)] + image_shifts(vecs)
    ion_sep = min((min(float(np.linalg.norm(xyz[i] - xyz[j] + s)) for s in shifts)
                   for i, j in itertools.combinations(ion_idx, 2)), default=None)
    q_total = float(sum(a.charge for a in atoms))
    n_bulk = len(water_o) - n_xtal
    avogadro_per_A3 = 6.02214076e-4                      # molecules per cubic Angstrom at 1 mol per litre
    warn1, warn2 = leap_warnings(log1), leap_warnings(log2)
    unexpected = warn2["other"] + [f"close contact between heavy atoms: {c['atoms']} {c['distance_A']} A"
                                   for c in warn2["close_contacts_between_heavy_atoms"]]
    n_one_sided = 0                                  # the renamed residues have head and tail atoms (leap_head)
    if len(warn2["one_sided_connection"]) != n_one_sided:
        unexpected += warn2["one_sided_connection"]
    bonded_pairs = {(b.atom1.idx, b.atom2.idx) for b in parm.bonds} | {(b.atom2.idx, b.atom1.idx) for b in parm.bonds}

    def joined(k1, k2):                              # peptide bond between residue ordinals k1 and k2 = k1 + 1
        c = [a.idx for a in parm.residues[k1 - 1].atoms if a.name == "C"]
        n = [a.idx for a in parm.residues[k2 - 1].atoms if a.name == "N"]
        return len(c) == 1 and len(n) == 1 and (c[0], n[0]) in bonded_pairs

    chain_ok = all(joined(k, k + 1) for k in range(1, n_prot))
    cc = warn2["close_contacts_with_a_hydrogen"]
    cc_lig = [c for c in cc if lig_res.name in c["atoms"] or zn_res.name in c["atoms"].replace("Zn2+", "")]

    # -- the same solvation and the same ions as the restrained system of this start?
    against_restrained = None
    other = pathlib.Path(out_root) / "restrained" / start
    if kind == "bonded" and (other / "build.json").is_file() and (other / "system.inpcrd").is_file():
        ob = json.loads((other / "build.json").read_text())
        oc = ob["composition"]
        d_max = None
        if oc["atoms"] == len(atoms):
            d_max = float(np.abs(parmed.load_file(str(other / "system.inpcrd")).coordinates.reshape(-1, 3) - xyz).max())
        against_restrained = {
            "folder": str(other), "note": "each pair is [this system, the restrained system]",
            "atoms": [len(atoms), oc["atoms"]], "water_molecules": [len(water_o), oc["water_molecules"]],
            "crystal_water_molecules": [n_xtal, oc["crystal_water_molecules"]],
            "Na+": [n_na, oc["Na+"]], "Cl-": [n_cl, oc["Cl-"]],
            "cell_edge_A": [round(float(parm.box[0]), 7), ob["box"]["edge_A"]],
            "ions_replace_the_same_water_molecules": [(x["ion"], x["replaces_added_water_number"]) for x in ion_sites]
            == [(x["ion"], x["replaces_added_water_number"]) for x in ob["ions"]["sites"]],
            "largest_coordinate_difference_A": None if d_max is None else round(d_max, 6),
            "same_atoms_at_the_same_places_within_0.005_A": d_max is not None and d_max < 0.005}
        osel = json.loads((other / "selections.json").read_text())
        if osel["solute"] == sel["solute"] and osel["protein_CA"] == ca:
            ox = parmed.load_file(str(other / "system.inpcrd")).coordinates.reshape(-1, 3)
            R2, t2 = kabsch(ox[ca], ref_ca)
            dsol = np.linalg.norm((ox @ R2.T + t2)[sel["solute"]] - fit[sel["solute"]], axis=1)
            worst = int(np.argmax(dsol))
            against_restrained.update({
                "solute_atoms_largest_deviation_after_the_fit_on_CA_A": round(float(dsol.max()), 4),
                "solute_atom_with_the_largest_deviation": f"{atoms[sel['solute'][worst]].residue.name}"
                                                          f"{atoms[sel['solute'][worst]].residue.idx + 1}:{atoms[sel['solute'][worst]].name}",
                "solute_atoms_farther_than_0.01_A": int((dsol > 0.01).sum())})

    # -- against another built set of the same start: only the terms in which the two sets differ may differ
    compared = None
    if compare_with:
        odir = pathlib.Path(out_root) / compare_with / start
        ob = json.loads((odir / "build.json").read_text())
        compared = compare_systems(parm, out, odir, sel, (pset, variant), (ob.get("parameter_set"),
                                                                           (ob.get("mcpb") or {}).get("parameter_set", {}).get("variant")))
        compared["other_model"] = compare_with
        compared["other_system_prmtop_sha256"] = sha256(odir / "system.prmtop")
        if n2_record and (odir / "build" / "n2_terms.frcmod").is_file():
            compared["n2_terms_file_same_sha256"] = n2_record["sha256"] == sha256(odir / "build" / "n2_terms.frcmod")
        write_json(rec / f"compare_with_{compare_with}.json", compared)

    checks = {
        "tleap_errors_0": count1["errors"] == 0 and count2["errors"] == 0,
        "tleap_no_unexpected_warning": not unexpected,
        "protein_chain_unbroken": chain_ok,
        "net_charge_zero": abs(q_total) < 5e-3,
        "ion_numbers_as_computed": len(sel["ions_Na"]) == n_na and len(sel["ions_Cl"]) == n_cl,
        "pass2_equals_pass1_except_ions": same_as_pass1 <= 1e-6 and ion_off <= 0.002,
        "input_heavy_atoms_at_their_coordinates": input_dev < 0.002,
        "selections_resolve_to_the_right_atoms": bool(sel_ok and names_ok),
        "geometry_agrees_with_starts_csv": all(v["agrees"] for v in against.values()),
        "zinc_charge_and_type": (abs(znatom.charge - 2.0) < 1e-6 and znatom.type == "Zn2+") if kind == "restrained" else None,
        "zinc_lj_is_the_frcmod_value": (abs(znatom.rmin - ffiles["frcmod.ions234lm_126_tip3p"]["rmin_half_A"]) < 1e-4 and
                                        abs(znatom.epsilon - ffiles["frcmod.ions234lm_126_tip3p"]["epsilon_kcal_per_mol"]) < 1e-7)
        if kind == "restrained" else None,
        "zinc_not_bonded": (not znatom.bond_partners) if kind == "restrained" else None,
        "zinc_has_the_four_bonds": bonded_rec["zinc_has_the_four_bonds"] if kind == "bonded" else None,
        "zinc_bonds_and_angles_as_the_parameter_files": bonded_rec["all_zinc_bonds_and_angles_as_the_parameter_files"]
        if kind == "bonded" else None,
        "zinc_angles_6_N_Zn_N_6_Zn_N_C_1_Zn_N2_S_1_Zn_N2_H": (bonded_rec["n_angles_by_kind"] ==
                                                              {"N-Zn-N": 6, "Zn-N-C": 6, "Zn-N2-S": 1, "Zn-N2-H": 1})
        if kind == "bonded" else None,
        # the next five are specific to the set: a set without a variant has the terms at N2 as built (every torsion
        # through the zinc zero, no improper with the zinc, Zn-N2-H not the 120 deg of step 2e); a variant has what
        # its terms-at-N2 file holds, and nothing else
        "torsions_through_zinc_all_zero": bonded_rec["torsions_through_zinc"]["all_zero"]
        if kind == "bonded" and not variant else None,
        "torsions_through_zinc_zero_or_as_the_n2_file": bonded_rec["torsions_through_zinc"]["zero_or_as_the_n2_file"]
        if kind == "bonded" and variant else None,
        "impropers_with_zinc_none_or_those_of_the_n2_file":
        bonded_rec["torsions_through_zinc"]["impropers_with_zinc_none_or_those_of_the_n2_file"] if kind == "bonded" else None,
        "n2_file_is_of_the_class_of_the_set": n2_record["holds"]["of_the_class_of_this_set"] if n2_record else None,
        "terms_at_N2_are_those_of_the_n2_file": n2_record["in_the_topology"]["all_as_n2_terms_file"] if n2_record else None,
        "zn_n2_h_is_not_the_120_deg_of_step_2e": n2_record["in_the_topology"]["Zn-N2-H_is_not_120_deg"]
        if n2_record and not variant else None,
        "site_charges_are_those_of_the_templates": bonded_rec["charges"]["all_as_the_templates"] if kind == "bonded" else None,
        "site_charge_plus_one": (abs(bonded_rec["charge_of_zinc_site_e"] - 1.0) < 5e-4) if kind == "bonded" else None,
        "ligand_equal_atoms_have_equal_charges": bonded_rec["charges"]["pairs_equal"] if pset else None,
        "solute_hydrogens_where_the_restrained_build_has_them": (
            against_restrained["solute_atoms_farther_than_0.01_A"] == 0)
        if against_restrained and "solute_atoms_farther_than_0.01_A" in against_restrained else None,
        "differs_from_the_compared_set_only_where_the_sets_differ": compared["only_the_allowed_terms_differ"] if compared else None,
        "exclusions_and_1_4_pairs_across_the_zinc_follow_amber": (
            bonded_rec["pairs_across_the_zinc"]["1_or_2_bonds_excluded_and_no_1_4_term"]
            and bonded_rec["pairs_across_the_zinc"]["3_bonds_one_1_4_term_each_with_the_amber_scaling"]
            and bonded_rec["pairs_across_the_zinc"]["4_bonds_neither_excluded_nor_1_4"]) if kind == "bonded" else None,
        "na_cl_are_joung_cheatham_tip3p": ions_as_jc,
        "ligand_repaired_terms_in_force": all(c["in_force"] for c in repaired_check),
        "ligand_terms_equal_free_anion_topology": all(same_as_free.values()) if kind == "restrained" else None,
        "ligand_charges_are_those_of_the_mol2": charges_as_mol2 if kind == "restrained" else None,
        "ions_keep_their_distances": bool(len(near) == 0 or (near.min() >= ION_MIN_DIST_TO_SOLUTE_A - 0.002 and
                                                               (ion_sep is None or ion_sep >= ION_MIN_SEPARATION_A - 0.002))),
        "no_crystal_water_replaced": len(crystal_o) == n_xtal and
        max(dist(fit[i], a["xyz"]) for i, a in zip(crystal_o, inp["waters"])) < 0.002,
    }
    failed = [k for k, v in checks.items() if v is False]

    inputs = {rel(p): sha256(p) for p in inp["files"].values()}
    amber_files = {}
    for sub, name in (("cmd", "leaprc.protein.ff14SB"), ("cmd", "leaprc.gaff2"), ("cmd", "leaprc.water.tip3p"),
                      ("parm", "gaff2.dat"), ("parm", "frcmod.ff14SB"), ("parm", "frcmod.tip3p"),
                      ("parm", "frcmod.ions234lm_126_tip3p"), ("parm", "frcmod.ionsjc_tip3p"),
                      ("parm", "frcmod.ions1lm_126_tip3p"), ("lib", "atomic_ions.lib"), ("lib", "solvents.lib")):
        p = pathlib.Path(AMBERHOME) / "dat" / "leap" / sub / name
        amber_files[f"dat/leap/{sub}/{name}"] = sha256(p) if p.is_file() else None
    gaff_line = next((ln.strip() for ln in log2.splitlines() if "General Force Field" in ln), None)

    outputs = {name: sha256(out / name) for name in ("system.prmtop", "system.inpcrd", "selections.json",
                                                      "extra_forces.json", "reference.pdb")}
    build_json = {
        "schema": "cqc-sfa-md built system, version 1",
        "what": label or "hCA II (PDB 6RL9, chain A) with the sulfonamidate anion of sulfanilamide on the zinc, "
                         "crystal water, TIP3P water and NaCl in a truncated octahedron",
        "model": model, "model_kind": kind, "parameter_set": pset, "start": start,
        "status_note": note,
        "written_by": "md/build_topology.py", "script_sha256": sha256(__file__),
        "software": versions(),
        "units": "Angstrom, degrees, e, amu, kcal mol-1 (Amber parameters); nm and kJ mol-1 nm-2 only where the key says so",
        "random_numbers": "none; the ion sites follow a fixed rule (ions.rule) and tleap places nothing by chance",
        "inputs_sha256": inputs,
        "amber_data_files_sha256": amber_files,
        "force_field": {
            "tleap_sources": ["leaprc.protein.ff14SB", "leaprc.gaff2", "leaprc.water.tip3p"],
            "gaff2_version_line": gaff_line,
            "protein": "ff14SB", "water": "TIP3P (rigid; the topology holds the H-H bond)",
            "ligand": ("GAFF2 atom types; AM1-BCC charges of the free anion (md/ligand/sfa_anion.mol2, sum -1); "
                       "sfa_anion_gaff2.frcmod, then sfa_anion_repair.frcmod (repaired terms at sulfur)")
            if kind == "restrained" else "as delivered by MCPB.py (charges of the zinc site refitted; see mcpb)",
            "zinc": ("nonbonded Zn2+, 12-6 Lennard-Jones, Li/Merz 'CM' (compromise) set for TIP3P of "
                     "frcmod.ions234lm_126_tip3p, which leaprc.water.tip3p loads last")
            if kind == "restrained" else "bonded to its four nitrogen atoms with the terms of MCPB.py (see mcpb, zinc_site_bonded)",
            "zinc_in_prmtop": zn_rec, "zinc_line_of_frcmod": ffiles["frcmod.ions234lm_126_tip3p"],
            "sodium_chloride": "Joung-Cheatham for TIP3P (frcmod.ionsjc_tip3p, loaded by leaprc.water.tip3p)",
            "ions_in_prmtop": ion_rec, "ion_lines_of_frcmod": {k: v for k, v in ffiles.items() if "ionsjc" in k},
            "notes": leap_notes},
        "mcpb": mcpb_record,
        "terms_at_N2": n2_record,
        "zinc_site_bonded": bonded_rec,
        "against_the_restrained_system_of_this_start": against_restrained,
        "compared_with": None if not compared else dict(
            {k: v for k, v in compared.items() if k != "terms_that_differ"},
            terms_that_differ=f"listed in md/topology/{model}/{start}/compare_with_{compare_with}.json "
                              f"({len(compared.get('terms_that_differ', []))} terms)"),
        "ligand_parameters": {
            "repaired_terms": repaired_check,
            "all_ligand_terms_equal_to_md_ligand_sfa_anion_repaired_prmtop": same_as_free,
            "terms_that_differ_from_that_topology": {sec: sorted(k for k in set(terms[sec]) | set(ref_terms[sec])
                                                                 if terms[sec].get(k) != ref_terms[sec].get(k)) for sec in terms},
            "terms": terms if kind == "bonded" else "as md/ligand/sfa_anion_repaired.prmtop (see the line above)",
            "charges_e": lig_charge, "charge_sum_e": round(sum(lig_charge.values()), 6),
            "charges_are_those_of_sfa_anion_mol2": charges_as_mol2,
            "n_terms": {k: len(v) for k, v in terms.items()}},
        "composition": {
            "atoms": len(atoms), "residues": len(parm.residues),
            "protein_atoms": len(protein_all), "protein_residues_with_ACE": n_prot, "zinc_atoms": 1,
            "ligand_atoms": len(sel["lig_all"]), "solute_atoms": len(sel["solute"]),
            "water_molecules": len(water_o), "crystal_water_molecules": n_xtal, "added_water_molecules": n_bulk,
            "Na+": n_na, "Cl-": n_cl,
            "plan_atoms": 25600, "atoms_over_plan_percent": round(100.0 * (len(atoms) / 25600.0 - 1.0), 1)},
        "charge": {"solute_with_crystal_water_e": round(q_dry, 6), "net_e": round(q_total, 6)},
        "box": {
            "shape": "truncated octahedron (tleap solvateOct; triclinic cell with three equal edges)",
            "solvateOct_buffer_A": buffer_A,
            "tleap_solvation_lines": solvation_lines(log2),
            "edge_A": round(float(parm.box[0]), 7), "angles_deg": [round(float(x), 7) for x in parm.box[3:]],
            "volume_A3": round(volume, 1), "volume_nm3": round(volume / 1000.0, 3),
            "mass_amu": round(mass, 2), "density_as_built_g_per_cm3": round(mass * 1.66053907 / volume, 4),
            "nearest_image_distance_A": round(float(parm.box[0]), 3),
            "solute_largest_atom_atom_distance_A": round(extent, 2),
            "solute_to_own_image_smallest_distance_as_built_A": round(image_min, 2),
            "solute_to_own_image_smallest_distance_any_orientation_A": round(float(parm.box[0]) - extent, 2),
            "note": "as built the box is larger than it will be at 1 bar (see density_as_built); 'any orientation' "
                    "is the cell edge minus the largest extent of the solute, a lower bound for a rigid solute"},
        "ions": {
            "target_mol_per_l": SALT_MOL_PER_L,
            "rule": f"pairs = round(c * N_water / c_water) with c = {SALT_MOL_PER_L} mol/L, c_water = {WATER_MOL_PER_L} mol/L "
                    "(water at 300 K) and N_water the water molecules after solvateOct (crystal water included): the "
                    "concentration in the volume of the solvent. Then one Na+ per negative charge of the solute.",
            "n_water_before_ions": n_water1, "pairs": n_pairs, "neutralising": {"Na+": max(0, -q_int), "Cl-": max(0, q_int)},
            "for_comparison_pairs_from_box_volume_as_built": round(SALT_MOL_PER_L * avogadro_per_A3 * volume1, 1),
            "concentration_in_box_volume_as_built_mol_per_l": round(n_pairs / (avogadro_per_A3 * volume1), 4),
            "placement": {
                "rule": "each ion replaces one water molecule added by solvateOct, never a crystal water: candidates "
                        f"are the added molecules whose O is at least {ION_MIN_DIST_TO_SOLUTE_A} A from every atom of "
                        "protein, zinc, ligand and crystal water; they are ordered by "
                        "sha256('cqc-sfa-md ion sites|<seed>|<k>') (k = place among the added molecules, from 0) and "
                        f"taken in that order if at least {ION_MIN_SEPARATION_A} A (minimum image) from every site "
                        "already taken; the kinds alternate Na+, Cl- in the order taken",
                "seed": ION_SEED, "candidates": n_cand,
                "smallest_ion_to_solute_distance_A": round(float(near.min()), 2) if len(near) else None,
                "smallest_ion_to_ion_distance_A": round(ion_sep, 2) if ion_sep is not None else None,
                "why_not_tleap": "addIonsRand has no seed (its help text: the result differs between invocations) and "
                                 "may replace a crystal water; addIons puts every ion on the protein surface"},
            "sites": ion_sites},
        "zinc_restraints": None if kind != "restrained" else {
            "k_kj_per_mol_nm2": ZN_RESTRAINT_K_KJ_NM2,
            "k_kcal_per_mol_A2_same_form": round(ZN_RESTRAINT_K_KJ_NM2 / 418.4, 3),
            "form": "U = 0.5 k (r - r0)^2",
            "expected_thermal_spread_A_at_300K": round(sigma_A, 4),
            "r0_A": {text: round(d, 4) for text, (j, d) in r0.items()},
            "r0_source": "distances in the deposited coordinates (md/system/zinc.pdb, protein_md.pdb, ligand_X.pdb); "
                         "the same four values for the three starts",
            "stiffness": stiffness},
        "tleap": {"pass1": count1, "pass2": count2,
                  "warnings_pass2": {"close_contacts_with_a_hydrogen": len(cc),
                                     "shortest_close_contact_A": min((c["distance_A"] for c in cc), default=None),
                                     "close_contacts_with_zinc_or_ligand": cc_lig,
                                     "close_contacts_between_heavy_atoms": warn2["close_contacts_between_heavy_atoms"],
                                     "charge_not_zero": warn2["charge_not_zero"],
                                     "one_sided_connection": warn2["one_sided_connection"],
                                     "one_sided_connection_expected": n_one_sided,
                                     "other": warn2["other"]},
                  "warnings_pass1_other": warn1["other"]},
        "water_at_the_site": site_water,
        "geometry_as_built": {"frame": "fitted on the CA atoms of reference.pdb, as md/analysis/SPEC.md does",
                              "ca_fit_rms_A": round(ca_fit_rms, 5),
                              "largest_deviation_of_an_input_heavy_atom_A": round(input_dev, 5),
                              "values": json.loads(json.dumps(geo, default=float), parse_float=lambda x: round(float(x), 4)),
                              "against_starts_csv": against},
        "checks": checks, "checks_failed": failed,
        "outputs_sha256": outputs,
    }
    write_json(out / "build.json", build_json)
    shutil.copy(out / "build.json", rec / "build.json")
    if n2_record:
        shutil.copy(bdir / "n2_terms.frcmod", rec / "n2_terms.frcmod")
    for d in (out, rec):
        if note:
            (d / "STATUS_NOTE.txt").write_text(f"{model}/{start}, built {now_ist()} by md/build_topology.py\n{note}\n")
        elif (d / "STATUS_NOTE.txt").exists():
            (d / "STATUS_NOTE.txt").unlink()
    for name in ("leap_stage1.in", "leap.in", "leap.log", "ions.pdb"):
        shutil.copy(bdir / name, rec / name)
    sums = [f"{sha256(out / n)}  {n}" for n in ("system.prmtop", "system.inpcrd", "selections.json", "extra_forces.json",
                                                "reference.pdb", "build.json")]
    (rec / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    for f in ("stage1.pdb",):
        (bdir / f).unlink()
    # -- into place: each file by one rename, so the system folder never lacks a file while a build runs
    for name in ("system.prmtop", "system.inpcrd", "selections.json", "extra_forces.json", "reference.pdb", "build.json"):
        os.replace(out / name, final / name)
    if note:
        os.replace(out / "STATUS_NOTE.txt", final / "STATUS_NOTE.txt")
    elif (final / "STATUS_NOTE.txt").exists():
        (final / "STATUS_NOTE.txt").unlink()
    out.rmdir()
    out = final

    print(f"   atoms {len(atoms)} (plan 25,600), water {len(water_o)} ({n_xtal} crystal), Na+ {n_na}, Cl- {n_cl}, "
          f"net charge {q_total:+.6f} e")
    print(f"   cell edge {parm.box[0]:.3f} A, volume {volume / 1000:.1f} nm3, density as built "
          f"{mass * 1.66053907 / volume:.3f} g/cm3, solute to its image {image_min:.1f} A")
    print(f"   tleap warnings: {len(cc)} close contacts with a hydrogen, {len(unexpected)} unexpected")
    print(f"   ligand RMSD {geo['lig_rmsd_symmetry_corrected_A']:.4f} A, d1 {geo['d1_N1_Phe131CZ_A']:.3f} A, "
          f"d2 {geo['d2_N1_Asn67OD1_A']:.3f} A, Zn-N2 {geo['Zn_N2_A']:.3f} A")
    print(f"   checks failed: {failed if failed else 'none'}")

    if openmm_check:
        om = openmm_load_check(out, rec)
        if (om.get("pairs_across_the_zinc") or {}).get("as_amber_rules") is False:
            build_json["checks_failed"].append("openmm_check.json: pairs_across_the_zinc")
    return build_json


def openmm_load_check(out, rec):
    """The topology loads in OpenMM with the settings of the driver; one energy on the CPU platform."""
    import openmm
    from openmm import app, unit
    proto = json.loads((MD / "protocol.json").read_text())["system"]
    top = app.AmberPrmtopFile(str(out / "system.prmtop"))
    crd = app.AmberInpcrdFile(str(out / "system.inpcrd"))
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=proto["nonbonded_cutoff_nm"] * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True, hydrogenMass=proto["hydrogen_mass_amu"] * unit.amu,
                              ewaldErrorTolerance=proto["ewald_error_tolerance"])
    for k, f in enumerate(system.getForces()):
        f.setForceGroup(k)
    integ = openmm.VerletIntegrator(0.001)
    plat = openmm.Platform.getPlatformByName("CPU")
    ctx = openmm.Context(system, integ, plat, {"Threads": str(THREADS)})
    ctx.setPeriodicBoxVectors(*crd.boxVectors)
    ctx.setPositions(crd.positions)
    kj = unit.kilojoule_per_mole
    total = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kj)
    parts = {type(f).__name__: round(ctx.getState(getEnergy=True, groups={k}).getPotentialEnergy().value_in_unit(kj), 1)
             for k, f in enumerate(system.getForces())}
    masses = [system.getParticleMass(i).value_in_unit(unit.amu) for i in range(system.getNumParticles())]
    res = {"what": "system.prmtop and system.inpcrd read by OpenMM with the settings of md/protocol.json; energy of the "
                   "built coordinates (before any minimisation; hydrogens as tleap placed them)",
           "openmm": openmm.Platform.getOpenMMVersion(), "platform": "CPU", "threads": THREADS,
           "n_particles": system.getNumParticles(), "n_constraints": system.getNumConstraints(),
           "total_mass_amu": round(float(sum(masses)), 3),
           "box_vectors_nm": [[round(float(c), 6) for c in v.value_in_unit(unit.nanometer)] for v in crd.boxVectors],
           "potential_energy_kj_per_mol": round(total, 1), "by_force_kj_per_mol": parts,
           "energy_is_finite": bool(math.isfinite(total)),
           "note": "the energy is not repeatable to the last digit with more than one thread; this file is not part "
                   "of the built system"}
    # bonded model: the pairs across the zinc as the engine will treat them (exceptions of the NonbondedForce)
    import parmed
    zn = json.loads((out / "selections.json").read_text())["zn"]
    parm = parmed.load_file(str(out / "system.prmtop"))
    if parm.atoms[zn].bond_partners:
        ps = zinc_pair_sets(parm, zn)
        near = {i for d in ps for pr in ps[d] for i in pr}
        nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        e2 = unit.elementary_charge ** 2
        exc = {}
        for n in range(nb.getNumExceptions()):
            i, j, qq, sig, eps = nb.getExceptionParameters(n)
            if i in near and j in near:
                exc[tuple(sorted((i, j)))] = (qq.value_in_unit(e2), sig.value_in_unit(unit.nanometer), eps.value_in_unit(kj))
        par = {i: nb.getParticleParameters(i) for i in near}

        def scaled(p):                               # what Amber's 1-4 rule gives for the pair
            (q1, s1, e1), (q2, s2, e2_) = par[p[0]], par[p[1]]
            return (q1.value_in_unit(unit.elementary_charge) * q2.value_in_unit(unit.elementary_charge) / SCEE,
                    0.5 * (s1.value_in_unit(unit.nanometer) + s2.value_in_unit(unit.nanometer)),
                    math.sqrt(e1.value_in_unit(kj) * e2_.value_in_unit(kj)) / SCNB)

        def close(x, y):
            return abs(x - y) <= 1e-6 + 1e-5 * abs(y)

        r12 = all(p in exc and exc[p][0] == 0.0 and exc[p][2] == 0.0 for p in ps[1] + ps[2])
        r3 = all(p in exc and close(exc[p][0], scaled(p)[0]) and close(exc[p][2], scaled(p)[2])
                 and (scaled(p)[2] == 0.0 or close(exc[p][1], scaled(p)[1])) for p in ps[3])
        r4 = not any(p in exc for p in ps[4])
        res["pairs_across_the_zinc"] = {
            "what": "exceptions of OpenMM's NonbondedForce for the pairs whose bonds run through the zinc: one or two "
                    f"bonds apart, no interaction; three bonds apart, charge product / {SCEE} and epsilon / {SCNB}; four "
                    "bonds apart, no exception (full interaction)",
            "n_pairs_1_or_2_bonds_apart": len(ps[1]) + len(ps[2]), "n_pairs_3_bonds_apart": len(ps[3]),
            "n_pairs_4_bonds_apart": len(ps[4]), "1_or_2_bonds_no_interaction": bool(r12),
            "3_bonds_scaled_as_amber_1_4": bool(r3), "4_bonds_no_exception": bool(r4),
            "as_amber_rules": bool(r12 and r3 and r4)}
        print(f"   OpenMM: pairs across the zinc ({len(ps[1]) + len(ps[2])} excluded, {len(ps[3])} scaled 1-4, "
              f"{len(ps[4])} in full) follow the Amber rules: {res['pairs_across_the_zinc']['as_amber_rules']}")
    write_json(rec / "openmm_check.json", res)
    print(f"   OpenMM: {res['n_particles']} particles, {res['n_constraints']} constraints, energy as built "
          f"{total:.4g} kJ/mol (finite: {res['energy_is_finite']})")
    return res


# ------------------------------------------------------------------------------------------- the minimisation
def state_positions_A(path):
    import openmm
    from openmm import unit
    st = openmm.XmlSerializer.deserialize(pathlib.Path(path).read_text())
    return np.array(st.getPositions(asNumpy=True).value_in_unit(unit.angstrom), dtype=float)


def site_report(X, sel, ref_ca, ref_lig, r0_A):
    """Numbers of the zinc site and of the ligand for coordinates X (Angstrom), after the fit on CA."""
    ca = sel["protein_CA"]
    R, t = kabsch(X[ca], ref_ca)
    F = X @ R.T + t
    lig = sel["lig_heavy"]
    rm = {name: math.sqrt(sum(float(((F[lig[a]] - ref_lig[perm.get(a, a)]) ** 2).sum()) for a in LIG_HEAVY) / 11.0)
          for name, perm in LABELLINGS}
    zx = X[sel["zn"]]
    four = {"N2": X[lig["N2"]], "H94NE2": X[sel["his94_NE2"]], "H96NE2": X[sel["his96_NE2"]], "H119ND1": X[sel["his119_ND1"]]}
    names = list(four)
    out = {"lig_rmsd_symmetry_corrected_A": round(min(rm.values()), 4),
           "ca_rmsd_A": round(float(np.sqrt(((F[ca] - ref_ca) ** 2).sum(axis=1).mean())), 4),
           "d1_N1_Phe131CZ_A": round(dist(X[lig["N1"]], X[sel["phe131_CZ"]]), 3),
           "d2_N1_Asn67OD1_A": round(dist(X[lig["N1"]], X[sel["asn67_OD1"]]), 3),
           "zn_distances_A": {n: round(dist(zx, four[n]), 4) for n in names},
           "zn_distance_minus_r0_A": {n: round(dist(zx, four[n]) - r0_A[n], 4) for n in names} if r0_A else None,
           "N_Zn_N_deg": {f"{names[i]}-Zn-{names[j]}": round(angle(four[names[i]], zx, four[names[j]]), 2)
                          for i in range(4) for j in range(i + 1, 4)},
           "Zn_O1_A": round(dist(zx, X[lig["O1"]]), 3), "Zn_O2_A": round(dist(zx, X[lig["O2"]]), 3),
           "S_N2_A": round(dist(X[lig["S"]], X[lig["N2"]]), 4),
           "C4_S_N2_deg": round(angle(X[lig["C4"]], X[lig["S"]], X[lig["N2"]]), 2),
           "tau_fold_deg": round(max(abs(torsion(X[lig["C3"]], X[lig["C4"]], X[lig["S"]], X[lig["N2"]])),
                                     abs(torsion(X[lig["C5"]], X[lig["C4"]], X[lig["S"]], X[lig["N2"]]))), 2),
           "Zn_N2_S_C4_deg": round(torsion(zx, X[lig["N2"]], X[lig["S"]], X[lig["C4"]]), 2),
           "N2_Thr199OG1_A": round(dist(X[lig["N2"]], X[sel["thr199_OG1"]]), 3),
           "N2_HN21_Thr199OG1_deg": round(angle(X[lig["N2"]], X[sel["lig_HN2"]], X[sel["thr199_OG1"]]), 1)}
    e = [(X[i] - X[lig["N2"]]) / np.linalg.norm(X[i] - X[lig["N2"]]) for i in (sel["zn"], sel["lig_HN2"], lig["S"])]
    out["n2_signed_volume_Zn_H_S"] = round(float(np.dot(e[0], np.cross(e[1], e[2]))), 4)   # md/analysis/SPEC.md, 9A
    out["Zn_N2_S_deg"] = round(angle(zx, X[lig["N2"]], X[lig["S"]]), 2)
    out["Zn_N2_HN21_deg"] = round(angle(zx, X[lig["N2"]], X[sel["lig_HN2"]]), 2)
    out["S_N2_HN21_deg"] = round(angle(X[lig["S"]], X[lig["N2"]], X[sel["lig_HN2"]]), 2)
    return out


def read_reference(path):
    ref_ca, ref_lig = [], {}
    for a in read_pdb(path):
        if a["name"] == "CA":
            ref_ca.append(a["xyz"])
        elif a["name"] in LIG_HEAVY:
            ref_lig[a["name"]] = a["xyz"]
    return np.array(ref_ca), ref_lig


def system_numbers(out):
    import parmed
    sel = json.loads((out / "selections.json").read_text())
    extra = json.loads((out / "extra_forces.json").read_text())
    ref_ca, ref_lig = read_reference(out / "reference.pdb")
    parm = parmed.load_file(str(out / "system.prmtop"), xyz=str(out / "system.inpcrd"))
    key = {sel["lig_heavy"]["N2"]: "N2", sel["his94_NE2"]: "H94NE2", sel["his96_NE2"]: "H96NE2", sel["his119_ND1"]: "H119ND1"}
    r0 = {key[b["j"]]: 10.0 * b["r0_nm"] for b in extra["harmonic_bonds"]} or None
    if r0 is None:                                   # bonded model: the equilibrium lengths of the four zinc bonds
        zn = sel["zn"]
        r0 = {key[(b.atom2 if b.atom1.idx == zn else b.atom1).idx]: float(b.type.req) for b in parm.bonds
              if zn in (b.atom1.idx, b.atom2.idx) and (b.atom2 if b.atom1.idx == zn else b.atom1).idx in key} or None
    return sel, extra, ref_ca, ref_lig, parm, r0


def min_check(model, start, out_root, records_root=RECORDS, short_md=False, iterations=300, threads=THREADS,
              reuse=False, pset=None, variant=None):
    """The minimisation of md/run_md.py on the built system, and what it does to the zinc site and the ligand.
    iterations = None: md/protocol.json as it is (up to 5000 L-BFGS iterations per pass; an hour or more on the
    shared Pi).  iterations = n: the same protocol with at most n iterations per pass and output every 0.1 ps,
    marked as a test; short_md then continues with the six equilibration stages shortened to 4 ps in all.
    For the bonded model r0 are the equilibrium lengths of the four zinc bonds of the topology."""
    model = folder_name(model, pset, variant)
    out = pathlib.Path(out_root) / model / start
    rec = pathlib.Path(records_root) / model / start
    tag = "min_full" if iterations is None else f"min_{iterations}"
    run = pathlib.Path(out_root) / "_checks" / model / start / tag
    have_min = reuse and (run / "rec_min.json").is_file()
    have_md = reuse and (run / "rec_eq06.json").is_file()
    if run.exists() and not have_min:
        shutil.rmtree(run)
    run.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    base = nice_prefix() + [sys.executable, str(MD / "run_md.py"), "--system", str(out), "--out", str(run),
                            "--run-id", f"{start}check", "--platform", "CPU", "--threads", str(threads)]
    protocol = MD / "protocol.json"
    if iterations is not None:
        proto = json.loads(protocol.read_text())
        proto["purpose"] = "test"
        proto["minimisation"]["max_iterations"] = int(iterations)
        proto["output"].update({"state_log_interval_ps": 0.1, "solute_frame_interval_ps": 0.1, "full_frame_interval_ps": 1.0})
        proto["test_note"] = (f"md/protocol.json with purpose 'test', at most {iterations} minimiser iterations per pass "
                              "and output every 0.1 ps (all atoms every 1 ps); nothing else is changed")
        protocol = run.parent / f"protocol_{tag}.json"
        if not have_min:
            write_json(protocol, proto)
        base += ["--protocol", str(protocol), "--equil-scale", "0.002"]
        if short_md and not have_min:
            base += ["--segment-ns", "0.001", "--prod-ns", "0.002"]      # two production segments of 1 ps at the end
    if have_min:
        print(f"== {model}/{start}: reading the finished minimisation in {run}", flush=True)
    else:
        print(f"== {model}/{start}: minimisation with md/run_md.py in {run}  ({now_ist()})", flush=True)
        p = subprocess.run(base + ["--stage", "min"], env=env, capture_output=True, text=True)
        (run.parent / f"{tag}.console").write_text(p.stdout + p.stderr)
        if p.returncode != 0:
            raise BuildError(f"run_md.py --stage min ended with exit code {p.returncode}; see {run.parent / (tag + '.console')}")
    r = json.loads((run / "rec_min.json").read_text())
    sel, extra, ref_ca, ref_lig, parm, r0 = system_numbers(out)
    X0 = parm.coordinates.copy()
    X1 = state_positions_A(run / "state_min.xml")
    moved = np.linalg.norm(X1 - X0, axis=1)
    lig = sel["lig_heavy"]
    lig_h = [i for i in sel["lig_all"] if i not in lig.values()]
    before, after = site_report(X0, sel, ref_ca, ref_lig, r0), site_report(X1, sel, ref_ca, ref_lig, r0)
    heavy_move = {"ZN": round(float(moved[sel["zn"]]), 4)}
    heavy_move.update({n: round(float(moved[i]), 4) for n, i in lig.items()})
    dep_ca = sel["residue_numbers_deposited_for_CA"]
    xtal = dict(zip(sel["crystal_water_O"], sel["crystal_water_numbers_deposited"]))

    def label(i):
        a = parm.atoms[i]
        if i in xtal:
            return f"HOH{xtal[i]}:O"
        k = a.residue.idx
        dep = 2 if k == 0 else (dep_ca[k - 1] if k <= len(dep_ca) else {len(dep_ca) + 1: 301, len(dep_ca) + 2: 304}.get(k, k + 1))
        return f"{a.residue.name}{dep}:{a.name}"

    held = sorted(set(sel["protein_heavy"]) | {sel["zn"]} | set(lig.values()) | set(sel["crystal_water_O"]))
    top5 = sorted(held, key=lambda i: -moved[i])[:5]
    limit = MIN_CHECK["max_displacement_zinc_and_ligand_heavy_A"][start]
    tol = MIN_CHECK["zn_distance_tolerance_A"]
    e0, e1 = r["energies_before"]["potential_energy_kj_per_mol"], r["energies_after"]["potential_energy_kj_per_mol"]
    res = {"what": "minimisation of the built system by md/run_md.py --stage min, CPU platform; "
                   + ("md/protocol.json unchanged" if iterations is None else
                      f"md/protocol.json with at most {iterations} L-BFGS iterations per pass (a test protocol; the "
                      "production protocol allows 5000)"),
           "model": model, "start": start, "threads": threads, "run_folder": str(run),
           "system_prmtop_sha256": sha256(out / "system.prmtop"), "system_inpcrd_sha256": sha256(out / "system.inpcrd"),
           "extra_forces_sha256": sha256(out / "extra_forces.json"), "selections_sha256": sha256(out / "selections.json"),
           "run_md_sha256": sha256(MD / "run_md.py"), "protocol_sha256": sha256(protocol),
           "limits_fixed_beforehand": MIN_CHECK,
           "r0_A": {k: round(v, 4) for k, v in r0.items()} if r0 else None,
           "r0_are": ("the targets of the four harmonic restraints of extra_forces.json" if extra["harmonic_bonds"] else
                      "the equilibrium lengths of the four zinc bonds of system.prmtop (bonded model)"),
           "energies_before_kj_per_mol": r["energies_before"], "energies_after_kj_per_mol": r["energies_after"],
           "lbfgs_iterations_per_pass": r["lbfgs_iterations_per_pass"],
           "passes_stopped_by_the_iteration_limit": r["passes_stopped_by_the_iteration_limit"],
           "n_restrained_atoms": r["settings"]["n_restrained_atoms"],
           "restrained_atoms_displacement_A": r["restrained_atoms_displacement_A"],
           "restrained_atoms_that_moved_most_A": {label(i): round(float(moved[i]), 3) for i in top5},
           "all_atoms_max_displacement_A": r["all_atoms_max_displacement_A"],
           "wall_s": r["wall_s"],
           "site_before": before, "site_after": after,
           "displacement_A": {"zinc_and_ligand_heavy": heavy_move,
                              "ligand_hydrogens": {parm.atoms[i].name: round(float(moved[i]), 4) for i in lig_h},
                              "his94_NE2": round(float(moved[sel["his94_NE2"]]), 4),
                              "his96_NE2": round(float(moved[sel["his96_NE2"]]), 4),
                              "his119_ND1": round(float(moved[sel["his119_ND1"]]), 4),
                              "crystal_water_O_max": round(float(moved[sel["crystal_water_O"]].max()), 4),
                              "crystal_water_O_rms": round(float(np.sqrt((moved[sel["crystal_water_O"]] ** 2).mean())), 4)},
           "checks": {"energy_lower": e1 < e0,
                      "zn_distances_within_tolerance_of_r0": (all(abs(v) <= tol for v in after["zn_distance_minus_r0_A"].values())
                                                              if r0 else None),
                      "zinc_and_ligand_heavy_moved_at_most_limit": max(heavy_move.values()) <= limit},
           "note": "with one thread the CPU platform repeats itself (md/RUN.md, test 0); with two threads a second run "
                   "gives slightly different numbers"}
    res["checks_failed"] = [k for k, v in res["checks"].items() if v is False]
    if short_md and iterations is not None:
        code = 0
        if not have_md:
            print(f"== {model}/{start}: 4 ps of the equilibration protocol (scale 0.002) in {run}  ({now_ist()})", flush=True)
            p = subprocess.run(base + ["--stage", "equil", "--resume"], env=env, capture_output=True, text=True)
            (run.parent / f"{tag}_shortmd.console").write_text(p.stdout + p.stderr)
            code = p.returncode
        res["short_md"] = short_md_report(run, sel, ref_ca, ref_lig, r0, code)
        if (run / "state_eq06.xml").is_file():           # what sits on the zinc at the end of the 4 ps
            import openmm
            from openmm import unit
            st = openmm.XmlSerializer.deserialize((run / "state_eq06.xml").read_text())
            XE = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.angstrom), dtype=float)
            cell = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom), dtype=float)
            dv = XE - XE[sel["zn"]]
            dv -= np.round(dv @ np.linalg.inv(cell)) @ cell      # nearest image (good enough within a few Angstrom)
            dz = np.linalg.norm(dv, axis=1)
            near = [i for i in np.argsort(dz) if 0 < dz[i] <= 2.6 and parm.atoms[int(i)].atomic_number > 1]
            water_o = set(sel["water_O"])
            res["short_md"]["heavy_atoms_within_2.6_A_of_the_zinc_at_the_end"] = {
                (label(int(i)) if int(i) not in water_o or int(i) in xtal else f"added water {parm.atoms[int(i)].residue.idx + 1}:O"):
                round(float(dz[i]), 3) for i in near}
            res["short_md"]["n_heavy_atoms_within_2.6_A_of_the_zinc_at_the_end"] = len(near)
        res["checks_failed"] += ["short_md: " + k for k in res["short_md"]["checks_failed"]]
    write_json(rec / ("min_check_full.json" if iterations is None else "min_check.json"), res)
    if short_md and iterations is not None and not res["checks_failed"] and (not reuse or (run / "rec_seg0002.json").is_file()):
        try:
            res["end_to_end"] = end_to_end(run, out, base, env, start, reuse)
        except Exception as exc:                     # the report above stays valid whatever happens here
            res["end_to_end"] = {"failed": f"{type(exc).__name__}: {exc}"}
        if res["end_to_end"].get("failed") or res["end_to_end"].get("checks_failed"):
            res["checks_failed"] += ["end_to_end"]
        write_json(rec / "min_check.json", res)
    print(f"   energy {e0:.6g} -> {e1:.6g} kJ/mol in {r['wall_s']:.0f} s, passes {r['lbfgs_iterations_per_pass']}")
    print(f"   Zn distances minus r0 (A): {after['zn_distance_minus_r0_A']}")
    print(f"   largest move of zinc or a ligand heavy atom {max(heavy_move.values()):.3f} A (limit {limit}); "
          f"ligand RMSD {before['lig_rmsd_symmetry_corrected_A']} -> {after['lig_rmsd_symmetry_corrected_A']} A")
    if "short_md" in res and "zinc_and_ligand_free_stages_4_to_6" in res["short_md"]:
        f = res["short_md"]["zinc_and_ligand_free_stages_4_to_6"]
        print(f"   short MD, stages 4 to 6 ({f['n_frames']} frames): mean minus r0 {f['zn_distance_mean_minus_r0_A']}, "
              f"sd {f['zn_distances_sd_A']}")
    print(f"   checks failed: {res['checks_failed'] if res['checks_failed'] else 'none'}")
    return res


def end_to_end(run, out, base, env, start, reuse):
    """Builder -> driver -> analysis on the real system: two production segments of 1 ps after the shortened
    equilibration, then md/analysis/analyse_md.py on the run folder.  A test that the files fit; the numbers of
    2 ps mean nothing."""
    if not reuse:
        p = subprocess.run(base + ["--stage", "prod", "--resume"], env=env, capture_output=True, text=True)
        (run.parent / f"{run.name}_prod.console").write_text(p.stdout + p.stderr)
        if p.returncode != 0:
            return {"failed": f"run_md.py --stage prod ended with exit code {p.returncode}"}
    adir = run.parent / f"{run.name}_analysis"
    if adir.exists():
        shutil.rmtree(adir)
    cmd = nice_prefix() + [sys.executable, str(MD / "analysis" / "analyse_md.py"), str(run), "--system", str(out),
                           "--start", start, "--planned-ns", "0.002", "--out", str(adir), "--no-figures"]
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    (run.parent / f"{run.name}_analysis.console").write_text(p.stdout + p.stderr)
    res = {"what": end_to_end.__doc__.strip(), "analysis_exit_code": p.returncode, "analysis_folder": str(adir),
           "analyse_md_sha256": sha256(MD / "analysis" / "analyse_md.py"),
           "analysis_last_lines": (p.stdout + p.stderr).strip().splitlines()[-6:]}
    sj = adir / "summary.json"
    if sj.is_file():
        sm = json.loads(sj.read_text())
        res["summary_keys"] = sorted(sm)
        res["warnings"] = sm.get("warnings")
        runinfo = sm.get("run") or {}
        res["run"] = {k: runinfo.get(k) for k in ("n_frames_solute", "n_frames_full", "seen_ns", "complete", "name_mismatches",
                                                    "start", "model", "planned_ns") if k in runinfo}
    res["checks"] = {"analysis_ran": p.returncode == 0, "summary_written": sj.is_file()}
    res["checks_failed"] = [k for k, v in res["checks"].items() if v is False]
    return res


def short_md_report(run, sel, ref_ca, ref_lig, r0, code):
    """Six shortened equilibration stages (4 ps in all): does the system run, and where are the four distances?"""
    out = {"what": "md/run_md.py --equil-scale 0.002: stages of 0.2, 0.6, 0.6 ps at 2 fs and 0.8, 0.8, 1.0 ps at 4 fs; "
                   "a smoke test of the built system, not an equilibration", "exit_code": code, "stages": []}
    names = ["N2", "H94NE2", "H96NE2", "H119ND1"]
    pos = {p: k for k, p in enumerate(sel["solute"])}
    idx = [pos[sel["zn"]], pos[sel["lig_heavy"]["N2"]], pos[sel["his94_NE2"]], pos[sel["his96_NE2"]], pos[sel["his119_ND1"]]]
    io = [pos[sel["lig_heavy"]["O1"]], pos[sel["lig_heavy"]["O2"]]]
    free, free_o, free_ang = [], [], []
    for k in range(1, 7):
        rp = run / f"rec_eq{k:02d}.json"
        if not rp.is_file():
            break
        r = json.loads(rp.read_text())
        st = {"stage": r.get("tag"),
              "settings": {k: r["settings"].get(k) for k in ("ensemble", "timestep_fs", "n_steps", "length_ps",
                                                             "restraint_atoms", "restraint_k_kj_per_mol_nm2")},
              "final_temperature_K": r["final"].get("temperature_K"),
              "final_density_g_per_cm3": r["final"].get("density_g_per_cm3"),
              "final_potential_energy_kj_per_mol": r["final"].get("potential_energy_kj_per_mol"),
              "mean_temperature_K": (r.get("mean_over_state_log") or {}).get("temperature_K"),
              "ns_per_day": r.get("ns_per_day"), "wall_s": r.get("wall_s"), "monitor_at_end": r.get("monitor")}
        dcd = run / f"solute_eq{k:02d}.dcd"
        if dcd.is_file():
            from MDAnalysis.lib.formats.libdcd import DCDFile
            with DCDFile(str(dcd)) as f:
                frames = [fr.xyz.astype(float) for fr in f]
            d = np.array([[np.linalg.norm(x[idx[0]] - x[j]) for j in idx[1:]] for x in frames])
            st["n_frames"] = len(frames)
            if len(frames):
                st["zn_distances_mean_A"] = {n: round(float(d[:, c].mean()), 4) for c, n in enumerate(names)}
                st["zn_distances_sd_A"] = {n: round(float(d[:, c].std()), 4) for c, n in enumerate(names)}
                if k >= 4:
                    free.append(d)
                    free_o.append(np.array([[np.linalg.norm(x[idx[0]] - x[j]) for j in io] for x in frames]))
                    free_ang.append(np.array([[angle(x[idx[1 + i]], x[idx[0]], x[idx[1 + j]])
                                               for i in range(4) for j in range(i + 1, 4)] for x in frames]))
        out["stages"].append(st)
    if free:
        d = np.vstack(free)
        out["zinc_and_ligand_free_stages_4_to_6"] = {
            "n_frames": int(len(d)),
            "zn_distances_mean_A": {n: round(float(d[:, c].mean()), 4) for c, n in enumerate(names)},
            "zn_distances_sd_A": {n: round(float(d[:, c].std()), 4) for c, n in enumerate(names)},
            "zn_distance_mean_minus_r0_A": {n: round(float(d[:, c].mean()) - r0[n], 4) for c, n in enumerate(names)} if r0 else None,
            "largest_abs_deviation_from_r0_A": round(float(max(np.abs(d[:, c] - r0[n]).max() for c, n in enumerate(names))), 4) if r0 else None}
        o = np.vstack(free_o)
        ang = np.vstack(free_ang)
        pairs = [f"{names[i]}-Zn-{names[j]}" for i in range(4) for j in range(i + 1, 4)]
        out["zinc_and_ligand_free_stages_4_to_6"].update({
            "Zn_O1_A": {"mean": round(float(o[:, 0].mean()), 3), "min": round(float(o[:, 0].min()), 3)},
            "Zn_O2_A": {"mean": round(float(o[:, 1].mean()), 3), "min": round(float(o[:, 1].min()), 3)},
            "frames_with_a_sulfonyl_O_within_2.6_A_of_the_zinc": int((o.min(axis=1) <= 2.6).sum()),
            "N_Zn_N_mean_deg": {pairs[c]: round(float(ang[:, c].mean()), 1) for c in range(6)},
            "note": "2.6 ps in all, the first picoseconds after the release: a look, not a measurement"})
    final = run / "state_eq06.xml"
    if final.is_file():
        out["site_at_the_end"] = site_report(state_positions_A(final), sel, ref_ca, ref_lig, r0)
    out["checks"] = {"all_six_stages_ran": code == 0 and len(out["stages"]) == 6,
                     "temperatures_finite": all(isinstance(s["final_temperature_K"], (int, float)) and
                                                math.isfinite(s["final_temperature_K"]) for s in out["stages"])}
    out["checks_failed"] = [k for k, v in out["checks"].items() if v is False]
    return out


def k_scan(model, start, out_root, records_root=RECORDS, ks=(100000.0, 200000.0, 300000.0, 500000.0, 700000.0, 1000000.0)):
    """How far the nonbonded force field pulls the four Zn-N distances from r0, for several force constants.
    Start: the minimised structure of --min-check.  Everything farther than 8 A from the zinc is frozen (mass 0);
    the residues with an atom within 8 A are minimised with the position restraints of the driver's minimisation
    and the four harmonic bonds at each k, to 1.0 kJ mol-1 nm-1, without bond constraints (a frozen atom cannot
    carry one).  A measurement at 0 K; it says how large the systematic pull of the force field is."""
    import openmm
    from openmm import app, unit
    out = pathlib.Path(out_root) / model / start
    rec = pathlib.Path(records_root) / model / start
    runs = [d for d in sorted((pathlib.Path(out_root) / "_checks" / model / start).glob("min_*"))
            if (d / "state_min.xml").is_file()]
    if not runs:
        raise BuildError("run --min-check first")
    run = next((d for d in runs if d.name == "min_full"), runs[-1])
    sel, extra, ref_ca, ref_lig, parm, r0 = system_numbers(out)
    if not extra["harmonic_bonds"]:
        raise BuildError("the k-scan is for the restrained model")
    proto = json.loads((MD / "protocol.json").read_text())
    top = app.AmberPrmtopFile(str(out / "system.prmtop"))
    crd = app.AmberInpcrdFile(str(out / "system.inpcrd"))
    X0 = parm.coordinates.copy() / 10.0
    Xm = state_positions_A(run / "state_min.xml") / 10.0
    near = np.linalg.norm(Xm - Xm[sel["zn"]], axis=1) <= 0.8
    res_near = {a.residue.idx for a in parm.atoms if near[a.idx]}
    mobile = np.array([a.residue.idx in res_near for a in parm.atoms])       # whole residues, so that no constraint is cut
    results = []
    for k in ks:
        system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=proto["system"]["nonbonded_cutoff_nm"] * unit.nanometer,
                                  constraints=None, rigidWater=False, ewaldErrorTolerance=proto["system"]["ewald_error_tolerance"])
        bonds = openmm.HarmonicBondForce()
        bonds.setUsesPeriodicBoundaryConditions(True)
        for b in extra["harmonic_bonds"]:
            bonds.addBond(b["i"], b["j"], b["r0_nm"], k)
        system.addForce(bonds)
        posres = openmm.CustomExternalForce("0.5*kp*periodicdistance(x,y,z,x0,y0,z0)^2")
        for nm in ("kp", "x0", "y0", "z0"):
            posres.addPerParticleParameter(nm)
        held = sorted(set(sel["protein_heavy"]) | {sel["zn"]} | set(sel["lig_heavy"].values()) | set(sel["crystal_water_O"]))
        for i in held:
            posres.addParticle(int(i), [proto["minimisation"]["restraint_k_kj_per_mol_nm2"]] + [float(v) for v in X0[i]])
        system.addForce(posres)
        for i in np.where(~mobile)[0]:
            system.setParticleMass(int(i), 0.0)
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"),
                             {"Threads": str(THREADS)})
        ctx.setPeriodicBoxVectors(*crd.boxVectors)
        ctx.setPositions(Xm)
        openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 1500)
        X = np.array(ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom))
        rep = site_report(X, sel, ref_ca, ref_lig, r0)
        sig = 10.0 * math.sqrt(KT_300_KJ / k)
        mu = 14.007 * 65.38 / (14.007 + 65.38)                 # reduced mass of a free N against Zn, amu
        period_fs = 2.0 * math.pi / math.sqrt(k / mu) * 1000.0  # sqrt(kJ mol-1 nm-2 / amu) = 1/ps
        results.append({"k_kj_per_mol_nm2": k, "thermal_spread_A_at_300K": round(sig, 4),
                        "period_fs_free_N_against_Zn": round(period_fs, 1), "steps_of_4_fs_per_period": round(period_fs / 4.0, 1),
                        "zn_distance_minus_r0_A": rep["zn_distance_minus_r0_A"],
                        "largest_abs_offset_A": round(max(abs(v) for v in rep["zn_distance_minus_r0_A"].values()), 4),
                        "N_Zn_N_deg": rep["N_Zn_N_deg"], "Zn_O2_A": rep["Zn_O2_A"]})
        del ctx
        print(f"   k = {k:9.0f}: offsets (A) {rep['zn_distance_minus_r0_A']}", flush=True)
    res = {"what": k_scan.__doc__.strip(), "model": model, "start": start, "started_from": str(run / "state_min.xml"),
           "threads": THREADS, "n_mobile_atoms": int(mobile.sum()),
           "r0_A": r0, "chosen_k_kj_per_mol_nm2": ZN_RESTRAINT_K_KJ_NM2, "scan": results,
           "period_note": "harmonic period for the reduced mass of a free nitrogen (14.007 amu) against zinc (65.38 amu), "
                          "the fastest the restraint alone can make; in the molecule the nitrogen carries its ring"}
    write_json(rec / "k_scan.json", res)
    return res


def buffer_scan(start, out_root, buffers=(4.5, 6.0, 8.0, 10.0, 12.0), records_root=RECORDS):
    """Size of the box for several buffers of solvateOct (pass 1 only, restrained model, no ions): atoms, cell
    edge, and how near the solute comes to its own image.  For the choice of the buffer; builds nothing."""
    import parmed
    from scipy.spatial import cKDTree, ConvexHull
    from scipy.spatial.distance import pdist
    inp = load_inputs(start)
    work = pathlib.Path(out_root) / "_checks" / "buffer_scan" / start
    work.mkdir(parents=True, exist_ok=True)
    table = write_assembly(inp, work / "assembly.pdb", {})
    n_sol_res = len(inp["protein"]) + 2
    rows = []
    for b in buffers:
        head, _ = leap_head("restrained", None, work / "assembly.pdb", work, b)
        tag = f"b{b:.1f}"
        script = head + [f"saveamberparm sys {work / (tag + '.prmtop')} {work / (tag + '.inpcrd')}", "quit"]
        (work / f"{tag}.in").write_text("\n".join(script) + "\n")
        log, _ = run_tleap(work / f"{tag}.in", work)
        p = parmed.load_file(str(work / f"{tag}.prmtop"), xyz=str(work / f"{tag}.inpcrd"))
        xyz = p.coordinates
        sol = np.array([a.idx for a in p.atoms if a.residue.idx < n_sol_res])
        vecs = box_vectors(p.box)
        tree = cKDTree(xyz[sol])
        image_min = min(float(tree.query(xyz[sol] + s, k=1)[0].min()) for s in image_shifts(vecs))
        extent = float(pdist(xyz[sol][ConvexHull(xyz[sol]).vertices]).max())
        n_water = len(p.residues) - n_sol_res
        pairs = int(round(SALT_MOL_PER_L * n_water / WATER_MOL_PER_L))
        volume = float(abs(np.linalg.det(vecs)))
        mass = float(sum(a.mass for a in p.atoms))
        rows.append({"solvateOct_buffer_A": b, "atoms_before_ions": len(p.atoms), "water_molecules": n_water,
                     "ion_pairs_by_the_rule": pairs,
                     "atoms_with_ions": len(p.atoms) - 3 * (2 * pairs + 1) + (2 * pairs + 1),
                     "cell_edge_A": round(float(p.box[0]), 3), "volume_nm3": round(volume / 1000.0, 1),
                     "density_as_built_g_per_cm3": round(mass * 1.66053907 / volume, 4),
                     "solute_to_own_image_as_built_A": round(image_min, 2),
                     "solute_largest_extent_A": round(extent, 2),
                     "solute_to_own_image_any_orientation_A": round(float(p.box[0]) - extent, 2),
                     "tleap_lines": solvation_lines(log)})
        for ext in (".prmtop", ".inpcrd"):
            (work / (tag + ext)).unlink()
        print(f"   buffer {b:4.1f} A: {rows[-1]['atoms_with_ions']} atoms, edge {rows[-1]['cell_edge_A']} A, "
              f"image {rows[-1]['solute_to_own_image_as_built_A']} A as built, "
              f"{rows[-1]['solute_to_own_image_any_orientation_A']} A in any orientation", flush=True)
    res = {"what": buffer_scan.__doc__.strip(), "start": start, "plan_atoms": 25600, "chosen_buffer_A": BUFFER_A,
           "note": "'any orientation' = cell edge minus the largest atom-atom distance of protein, zinc and ligand: the "
                   "nearest the solute can come to its image if it turns as a rigid body. All values are for the box as "
                   "tleap builds it; at 1 bar the box shrinks (see the density).",
           "inputs_sha256": {rel(p): sha256(p) for k, p in inp["files"].items() if k in ("protein", "zinc", "ligand", "waters")},
           "software": versions(), "rows": rows}
    write_json(pathlib.Path(records_root) / f"buffer_scan_{start}.json", res)
    return res


# ----------------------------------------------------------------------------------------- self-test, bonded
def selftest_bonded(start, out_root):
    """Tests the path of the bonded model on files written by MCPB.py itself, WITHOUT quantum chemistry:
    a copy of the step-1 files of md/qm/mcpb, 'MCPB.py -s 2e' (empirical force constants) and '-s 4'
    (tleap input), and mol2 files of the zinc-site residues with the UNFITTED charges of ff14SB, of the free
    anion and +2 for the zinc.  The result is a MOCK: it tests the builder, it is not a zinc model."""
    src = MD / "qm" / "mcpb"
    root = pathlib.Path(out_root) / "_selftest_bonded_mock"
    work = root / "mcpb"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    group = None
    for p in sorted(src.glob("*.in")):
        if "original_pdb" in p.read_text():
            group = p.stem
    if group is None:
        raise BuildError(f"no MCPB.py input in {src}")
    need = [f"{group}.in", f"{group}_small.pdb", f"{group}_small.res", f"{group}_standard.pdb", f"{group}_standard.fingerprint",
            f"{group}_large.pdb", f"{group}_large.fingerprint"]
    conf = dict(line.split(None, 1) for line in (src / f"{group}.in").read_text().splitlines() if line.strip())
    need += conf["original_pdb"].split() + conf["ion_mol2files"].split() + conf["naa_mol2files"].split() + conf["frcmod_files"].split()
    copied = {}
    for name in need:
        if not (src / name).is_file():
            raise BuildError(f"the self-test needs {src / name} (step 1 of MCPB.py)")
        shutil.copy(src / name, work / name)
        copied[rel(src / name)] = sha256(src / name)
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=os.path.join(AMBERHOME, "bin") + os.pathsep + os.environ.get("PATH", ""))
    logs = {}
    for step in ("2e", "4"):
        p = subprocess.run([sys.executable, os.path.join(AMBERHOME, "bin", "MCPB.py"), "-i", f"{group}.in", "-s", step],
                           cwd=str(work), env=env, capture_output=True, text=True)
        (work / f"selftest_step{step}.log").write_text(p.stdout + p.stderr)
        logs[step] = p.returncode
        if p.returncode != 0:
            raise BuildError(f"MCPB.py -s {step} failed in {work}; see selftest_step{step}.log")
    mc = read_mcpb_names(work, group)
    # mock mol2 files: library residues with the new names and atom types
    newtype = {}
    for line in (work / f"{group}_standard.fingerprint").read_text().splitlines():
        m = re.match(r"(\d+)-(\S+?)-(\S+)\s+\d+\s+(\S+)\s*->\s*(\S+)", line)
        if m and m.group(4) != m.group(5):
            newtype[(int(m.group(1)), m.group(3))] = m.group(5)
    L = ["source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         f"SAN = loadmol2 {LIGAND / 'sfa_anion.mol2'}"]
    for resid, old, new in mc:
        L += [f"u{resid} = copy {old}", f"set u{resid} name \"{new}\"", f"set u{resid}.1 name \"{new}\""]
        L += [f"set u{resid}.1.{an} type \"{tp}\"" for (r, an), tp in sorted(newtype.items()) if r == resid]
        L += [f"savemol2 u{resid} {work / (new + '.mol2')} 1"]
    L.append("quit")
    (work / "selftest_mock_mol2.in").write_text("\n".join(L) + "\n")
    run_tleap(work / "selftest_mock_mol2.in", work)
    (root / "MOCK_NOT_A_MODEL").write_text("Self-test of md/build_topology.py: MCPB.py file formats with unfitted charges and "
                                           "empirical force constants. Not a zinc model; never simulate it.\n")
    bj = build("bonded", start, root, mcpb_dir=work, mcpb_group=group, records_root=RECORDS / "selftest_bonded_mock",
               label="SELF-TEST MOCK of the bonded path (MCPB.py formats, unfitted charges, empirical force constants); not a model")
    res = {"what": selftest_bonded.__doc__.strip(), "start": start, "mcpb_work_folder": str(work),
           "copied_from_md_qm_mcpb_sha256": copied, "mcpb_exit_codes": logs,
           "renamed_residues": [{"residue": r, "from": o, "to": n} for r, o, n in mc],
           "new_atom_types": {f"{r}:{a}": t for (r, a), t in sorted(newtype.items())},
           "build": {k: bj[k] for k in ("composition", "charge", "mcpb", "zinc_site_bonded", "checks", "checks_failed")},
           "repaired_terms": bj["ligand_parameters"]["repaired_terms"],
           "tleap_warnings_other": bj["tleap"]["warnings_pass2"]["other"]}
    write_json(RECORDS / "selftest_bonded_mock" / f"selftest_{start}.json", res)
    print(f"   self-test of the bonded path, start {start}: checks failed: {bj['checks_failed'] or 'none'}; "
          f"zinc bonded to {bj['zinc_site_bonded']['zinc_bonded_to']}")
    return res


def selftest_bonded_given(start, out_root, mcpb_dir, group=None):
    """The bonded path on a folder that already holds the output of MCPB.py steps 2 to 4, for example the one
    md/qm/test_after_qm.py makes from SYNTHETIC quantum-chemical files (mol2 files written by MCPB.py itself).
    The result goes to a folder marked as a test; it is not a zinc model."""
    root = pathlib.Path(out_root) / "_selftest_bonded_synthetic"
    root.mkdir(parents=True, exist_ok=True)
    (root / "MOCK_NOT_A_MODEL").write_text("Self-test of md/build_topology.py on MCPB.py output made from synthetic input. "
                                           "Not a zinc model; never simulate it.\n")
    bj = build("bonded", start, root, mcpb_dir=mcpb_dir, mcpb_group=group, records_root=RECORDS / "selftest_bonded_synthetic",
               label="SELF-TEST of the bonded path on MCPB.py output made from SYNTHETIC quantum-chemical input; not a model")
    res = {"what": selftest_bonded_given.__doc__.strip(), "start": start, "mcpb_folder": str(mcpb_dir),
           "build": {k: bj[k] for k in ("composition", "charge", "mcpb", "zinc_site_bonded", "checks", "checks_failed")},
           "repaired_terms": bj["ligand_parameters"]["repaired_terms"],
           "tleap_warnings_other": bj["tleap"]["warnings_pass2"]["other"]}
    write_json(RECORDS / "selftest_bonded_synthetic" / f"selftest_{start}.json", res)
    print(f"   bonded path on the given MCPB.py folder, start {start}: checks failed: {bj['checks_failed'] or 'none'}; "
          f"zinc bonded to {bj['zinc_site_bonded']['zinc_bonded_to']}")
    return res


def read_mcpb_names(work, group):
    """(residue number, old name, new name) of the residues MCPB.py renames, from its tleap input and its two PDB files."""
    text = (work / f"{group}_tleap.in").read_text()
    conf = dict(line.split(None, 1) for line in (work / f"{group}.in").read_text().splitlines() if line.strip())
    old = {a["resi"]: a["resn"] for a in read_pdb(work / conf["original_pdb"].split()[0])}
    new = {a["resi"]: a["resn"] for a in read_pdb(work / re.search(r"loadpdb\s+(\S+)", text).group(1))}
    units = set(re.findall(r"(?m)^(\S+)\s*=\s*loadmol2", text))
    out = [(r, old[r], new[r]) for r in sorted(new) if new[r] in units and new[r] != old.get(r)]
    if not out:
        raise BuildError("MCPB.py renamed no residue")
    return out


# ------------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description="Solvated topologies of the hCA II-sulfonamidate complex (tleap), with checks.")
    ap.add_argument("--model", choices=["restrained", "bonded"])
    ap.add_argument("--start", choices=["X", "N", "R"])
    ap.add_argument("--out-root", default=DEFAULT_OUT_ROOT, help=f"folder of the built systems (default {DEFAULT_OUT_ROOT})")
    ap.add_argument("--set", dest="pset", choices=sorted(SET_METHOD), default=None,
                    help="model bonded: the parameter set, 2e (MCPB.py step 2e, empirical) or 2s (step 2s, Seminario); "
                         "the folder is then bonded_<set>")
    ap.add_argument("--n2-frcmod", default=None,
                    help="model bonded: frcmod file with the terms at N2 (Zn-N2-S, Zn-N2-H, S-N2-H, H-N2-S-C4), loaded last; "
                         "default for --set 2e: written by the builder, the terms as built before gate G-R")
    ap.add_argument("--variant", default=None,
                    help="model bonded: name of a variant of the set with other terms at N2 (needs --n2-frcmod); the folder "
                         "is then bonded_<set>_<variant>, and the checks ask for what the file holds instead of the terms as built")
    ap.add_argument("--compare-with", default=None, metavar="MODEL",
                    help="model bonded: after the build, compare term by term with the built system <out-root>/MODEL/<start>; "
                         "only the terms in which the two sets differ may differ")
    ap.add_argument("--note", default=None, help="a sentence on the status of this build, written to build.json "
                                                 "(status_note) and to STATUS_NOTE.txt beside the system")
    ap.add_argument("--mcpb-dir", default=None,
                    help="model bonded: folder with the output of MCPB.py step 4 (default: the first of "
                         "/root/cqc-sfa-md/qm/gr/mcpb_<set> and md/qm/mcpb that holds it)")
    ap.add_argument("--mcpb-group", default=None, help="model bonded: group name of MCPB.py (default: the one <group>_tleap.in found)")
    ap.add_argument("--no-openmm", action="store_true", help="skip the OpenMM load check after the build")
    ap.add_argument("--min-check", action="store_true", help="minimise the built system with md/run_md.py and report")
    ap.add_argument("--short-md", action="store_true", help="with --min-check: also 4 ps of the shortened equilibration")
    ap.add_argument("--min-iterations", type=int, default=300, help="with --min-check: L-BFGS iterations per pass (default 300)")
    ap.add_argument("--min-full", action="store_true", help="with --min-check: md/protocol.json unchanged (slow on the Pi)")
    ap.add_argument("--reuse-run", action="store_true", help="with --min-check: evaluate a finished check run again without running it")
    ap.add_argument("--threads", type=int, default=THREADS, help=f"threads of the OpenMM CPU platform in the checks (default {THREADS})")
    ap.add_argument("--k-scan", action="store_true", help="restrained model: offsets of the four distances for several k")
    ap.add_argument("--selftest-bonded", action="store_true",
                    help="test the bonded path on a mock made with MCPB.py (no QM), or, with --mcpb-dir, on that folder")
    ap.add_argument("--buffer-scan", action="store_true", help="box size for several buffers of solvateOct (builds nothing)")
    ap.add_argument("--buffer-A", type=float, default=BUFFER_A, help=f"buffer of solvateOct in Angstrom (default {BUFFER_A}, the plan's)")
    args = ap.parse_args(argv)
    os.environ.setdefault("OMP_NUM_THREADS", str(THREADS))
    try:
        if args.selftest_bonded:
            if args.mcpb_dir:
                selftest_bonded_given(args.start or "X", args.out_root, args.mcpb_dir, args.mcpb_group)
            else:
                selftest_bonded(args.start or "X", args.out_root)
            return 0
        if args.buffer_scan:
            buffer_scan(args.start or "X", args.out_root)
            return 0
        if not args.model or not args.start:
            ap.error("--model and --start are required")
        if args.model == "bonded" and not args.pset:
            ap.error("--model bonded needs --set 2e or --set 2s (md/PLAN.md, 6 Oct 2026, 3:02 AM, point 3)")
        if args.model != "bonded" and (args.pset or args.n2_frcmod or args.variant or args.compare_with):
            ap.error("--set, --n2-frcmod, --variant and --compare-with belong to --model bonded")
        if args.variant and not re.fullmatch(r"[A-Za-z0-9]+", args.variant):
            ap.error("--variant: letters and digits only")
        if args.min_check:
            res = min_check(args.model, args.start, args.out_root, short_md=args.short_md,
                            iterations=None if args.min_full else args.min_iterations, threads=args.threads,
                            reuse=args.reuse_run, pset=args.pset, variant=args.variant)
            return 1 if res["checks_failed"] else 0
        if args.k_scan:
            k_scan(args.model, args.start, args.out_root)
            return 0
        bj = build(args.model, args.start, args.out_root, mcpb_dir=args.mcpb_dir, mcpb_group=args.mcpb_group,
                   openmm_check=not args.no_openmm, buffer_A=args.buffer_A, pset=args.pset, n2_frcmod=args.n2_frcmod,
                   note=args.note, variant=args.variant, compare_with=args.compare_with)
        return 1 if bj["checks_failed"] else 0
    except BuildError as exc:
        print(f"build_topology.py: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
