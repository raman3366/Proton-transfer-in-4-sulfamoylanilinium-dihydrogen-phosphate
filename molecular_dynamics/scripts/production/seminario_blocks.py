#!/usr/bin/env python
"""seminario_blocks.py WORK OUT.json: which blocks of the Hessian MCPB.py reads for the zinc terms, and the constants.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/production/seminario_blocks.py \
        /root/cqc-sfa-md/qm/gr/mcpb_2s md/production/seminario_blocks.json

WORK is a folder in which md/qm/after_qm.sh 2s has run (step-1 files, CA2_small_opt.fchk).  Nothing there is changed:
the folder is copied to a temporary folder twice, and MCPB.py step 2s and step 2ms (modified Seminario; tabulated,
not used: md/PLAN.md, 6 Oct 2026, 3:02 AM, point 3) are run there with their force-constant routines wrapped, so
that every call is recorded with the atoms it was made for.  From the source of MCPB.py 2026
(pymsmt/mcpb/gene_final_frcmod_file.py) a bond (a, b) reads the 3 x 3 block [rows of a, columns of b] and, with
bondfc_avg 1, also [b, a]; an angle (a, b, c) reads [a, b] and [c, b] and, with anglefc_avg 1, also [b, a] and
[b, c].  The fchk file holds a PARTIAL Hessian: columns were computed for the displaced atoms only, and a block
between two atoms that were both not displaced is zero (md/ZINC.md, section 7.3).  A block is "computed" when at
least one of its two atoms was displaced.  The script lists every block read, whether it was computed, and its
norm; then the constants of step 2s beside those of step 2e and step 2ms and beside the project's own Seminario
code (md/qm/seminario.py) for both choices of the block.  No random numbers.
"""
import contextlib
import io
import json
import os
import re
import runpy
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "qm"))
import seminario as own                                         # noqa: E402  (md/qm/seminario.py)
import pymsmt.mcpb.gene_final_frcmod_file as g                  # noqa: E402
from pymsmt.mol.gauio import get_crds_from_fchk, get_matrix_from_fchk   # noqa: E402

work, out_path = os.path.abspath(sys.argv[1]), sys.argv[2]
GROUP = "CA2"
B2A = 0.52917721092


def real(v):
    return float(np.real(v))


def run_step(step):
    """MCPB.py -s <step> in a copy of the work folder, the force-constant routines wrapped. -> (calls, frcmod text)."""
    tmp = tempfile.mkdtemp(prefix=f"seminario_{step}_")
    for name in os.listdir(work):
        p = os.path.join(work, name)
        if os.path.isfile(p):
            shutil.copy(p, tmp)
    calls = []
    orig = {n: getattr(g, n) for n in ("get_bond_fc_with_sem", "get_ang_fc_with_sem", "get_ang_fc_with_mod_sem",
                                       "avg_bond_para", "avg_angle_para")}

    def bond(crds0, crds, fcmatrix, nat1, nat2, scalef, bondavg):
        res = orig["get_bond_fc_with_sem"](crds0, crds, fcmatrix, nat1, nat2, scalef, bondavg)
        calls.append({"kind": "bond", "atoms": [nat1, nat2], "avg": int(bondavg), "scale_factor": float(scalef),
                      "equilibrium": real(res[0]), "k": real(res[1]), "equilibrium_from_crystal": crds0 != [],
                      "imaginary_part_of_k": float(np.imag(res[1]))})
        return res

    def ang(crds0, crds, fcmatrix, nat1, nat2, nat3, scalef, angavg):
        res = orig["get_ang_fc_with_sem"](crds0, crds, fcmatrix, nat1, nat2, nat3, scalef, angavg)
        calls.append({"kind": "angle", "atoms": [nat1, nat2, nat3], "avg": int(angavg), "scale_factor": float(scalef),
                      "equilibrium": real(res[0]), "k": real(res[1]), "equilibrium_from_crystal": crds0 != [],
                      "imaginary_part_of_k": float(np.imag(res[1]))})
        return res

    def modang(crds0, crds, fcmatrix, nat1, nat2, nat3, scalef, angavg, s1, s2):
        res = orig["get_ang_fc_with_mod_sem"](crds0, crds, fcmatrix, nat1, nat2, nat3, scalef, angavg, s1, s2)
        calls.append({"kind": "angle", "atoms": [nat1, nat2, nat3], "avg": int(angavg), "scale_factor": float(scalef),
                      "equilibrium": real(res[0]), "k": real(res[1]), "angle_scalings": [real(s1), real(s2)],
                      "imaginary_part_of_k": float(np.imag(res[1]))})
        return res

    def take_real(values):                       # the guard of md/qm/mcpb_step2s.py: NumPy 2 cannot round a complex number
        for v in values:
            if np.iscomplexobj(v) and abs(np.imag(v)) > 1e-6 * max(abs(v), 1e-300):
                raise SystemExit(f"a force constant has a real imaginary part: {v}")
        return [real(v) for v in values]

    g.get_bond_fc_with_sem, g.get_ang_fc_with_sem, g.get_ang_fc_with_mod_sem = bond, ang, modang
    g.avg_bond_para = lambda t, lens, fcs, sigf=4: orig["avg_bond_para"](t, take_real(lens), take_real(fcs), sigf)
    g.avg_angle_para = lambda t, vals, fcs: orig["avg_angle_para"](t, take_real(vals), take_real(fcs))
    cwd, argv = os.getcwd(), sys.argv
    os.chdir(tmp)
    for f in (f"{GROUP}_mcpbpy.frcmod",):
        if os.path.exists(f):
            os.remove(f)
    sys.argv = ["MCPB.py", "-i", f"{GROUP}.in", "-s", step]
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            try:
                runpy.run_path(sys.prefix + "/bin/MCPB.py", run_name="__main__")
            except SystemExit as e:
                if e.code not in (None, 0):
                    raise
    finally:
        os.chdir(cwd)
        sys.argv = argv
        for n, f in orig.items():
            setattr(g, n, f)
    text = open(os.path.join(tmp, f"{GROUP}_mcpbpy.frcmod")).read()
    shutil.rmtree(tmp)
    return calls, text, buf.getvalue()


def zinc_lines(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r"^(\S+(?:-\S+)+)\s+([\d.]+)\s+([\d.]+)\s+Created by (.+) method using MCPB\.py", line)
        if m:
            out[m.group(1)] = {"k": float(m.group(2)), "equilibrium": float(m.group(3)), "method": m.group(4)}
    return out


# ---- the model, the Hessian, the displaced atoms
names = []
for line in open(os.path.join(work, f"{GROUP}_small.pdb")):
    if line.startswith(("ATOM", "HETATM")):
        names.append(f"{line[17:20].strip()}{int(line[22:26])}@{line[12:16].strip()}")
n = len(names)
fchk = os.path.join(work, f"{GROUP}_small_opt.fchk")
crds = np.array(get_crds_from_fchk(fchk, n), dtype=float).reshape(n, 3)            # Bohr
H = np.array(get_matrix_from_fchk(fchk, 3 * n), dtype=float)
norms = np.array([[np.linalg.norm(H[3 * a:3 * a + 3, 3 * b:3 * b + 3]) for b in range(n)] for a in range(n)])
has_column = [bool(norms[:, b].max() > 0) and bool(norms[[a for a in range(n) if a != b], b].max() > 0) for b in range(n)]
# an atom was displaced if its column of blocks is non-zero against atoms that were certainly not displaced:
# count, for every atom, the atoms with which it has a non-zero block; displaced atoms have all n of them
nonzero_partners = (norms > 0).sum(axis=0)
displaced = [a for a in range(n) if nonzero_partners[a] == n]
title = open(fchk).readline().strip()

calls_2s, frcmod_2s, log_2s = run_step("2s")
calls_2ms, frcmod_2ms, log_2ms = run_step("2ms")
zl = {"2s": zinc_lines(frcmod_2s), "2ms": zinc_lines(frcmod_2ms),
      "2s_in_work_folder": zinc_lines(open(os.path.join(work, f"{GROUP}_mcpbpy.frcmod")).read())}
e_dir = os.path.join(os.path.dirname(work), "mcpb_2e")
if os.path.isfile(os.path.join(e_dir, f"{GROUP}_mcpbpy.frcmod")):
    zl["2e"] = zinc_lines(open(os.path.join(e_dir, f"{GROUP}_mcpbpy.frcmod")).read())


def blocks_of(c):
    a = [x - 1 for x in c["atoms"]]
    if c["kind"] == "bond":
        bl = [(a[0], a[1])] + ([(a[1], a[0])] if c["avg"] else [])
    else:
        bl = [(a[0], a[1]), (a[2], a[1])] + ([(a[1], a[0]), (a[1], a[2])] if c["avg"] else [])
    return bl


def describe(calls):
    rows, all_ok = [], True
    for c in calls:
        a = [x - 1 for x in c["atoms"]]
        bl = []
        for r, col in blocks_of(c):
            computed = (r in displaced) or (col in displaced)
            ev = np.linalg.eigvals(-H[3 * r:3 * r + 3, 3 * col:3 * col + 3])
            bl.append({"rows_of": names[r], "columns_of": names[col], "computed": bool(computed),
                       "norm_hartree_per_bohr2": round(float(norms[r, col]), 6),
                       "row_atom_displaced": r in displaced, "column_atom_displaced": col in displaced,
                       "eigenvalues_complex": bool(np.abs(np.imag(ev)).max() > 1e-12)})
            all_ok = all_ok and computed and norms[r, col] > 0
        row = {"kind": c["kind"], "atoms": [names[i] for i in a], "atom_numbers": c["atoms"],
               "k_mcpb": round(c["k"], 3), "equilibrium_mcpb": round(c["equilibrium"], 4),
               "anglefc_or_bondfc_avg": c["avg"], "blocks_read": bl,
               "all_blocks_computed": all(b["computed"] and b["norm_hartree_per_bohr2"] > 0 for b in bl)}
        if "angle_scalings" in c:
            row["angle_scalings_of_the_modified_method"] = [round(x, 4) for x in c["angle_scalings"]]
        if c["kind"] == "bond":
            r1, k1 = own.bond(H, crds, a[0], a[1])
            r2, k2 = own.bond(H, crds, a[1], a[0])
            row["own_code"] = {"length_A": round(r1, 4), "k_block_as_mcpb": round(k1, 2), "k_other_block": round(k2, 2),
                               "k_mean_of_the_two": round(0.5 * (k1 + k2), 2)}
        else:
            t1, k1 = own.angle(H, crds, a[0], a[1], a[2])
            row["own_code"] = {"angle_deg": round(t1, 2), "k_block_as_mcpb": round(k1, 2)}
        rows.append(row)
    return rows, all_ok


rows_2s, ok_2s = describe(calls_2s)
rows_2ms, ok_2ms = describe(calls_2ms)
res = {
    "what": __doc__.split("\n\n")[0] + " Step 2ms is tabulated only.",
    "work_folder": work, "fchk_title": title,
    "n_atoms": n, "displaced_atoms_found_in_the_matrix": [{"number": a + 1, "atom": names[a]} for a in displaced],
    "n_blocks_non_zero": int((norms > 0).sum()), "n_blocks": n * n,
    "rule": "a block between two atoms is computed if at least one of them was displaced; found from the matrix itself: "
            "a displaced atom has a non-zero block with every atom",
    "step_2s": {"n_bonds": sum(1 for r in rows_2s if r["kind"] == "bond"),
                "n_angles": sum(1 for r in rows_2s if r["kind"] == "angle"),
                "every_block_read_was_computed": bool(ok_2s), "terms": rows_2s},
    "step_2ms": {"n_bonds": sum(1 for r in rows_2ms if r["kind"] == "bond"),
                 "n_angles": sum(1 for r in rows_2ms if r["kind"] == "angle"),
                 "every_block_read_was_computed": bool(ok_2ms), "terms": rows_2ms},
    "frcmod_zinc_lines": zl,
    "step_2s_here_equals_the_frcmod_of_the_work_folder": zl["2s"] == zl["2s_in_work_folder"],
}
with open(out_path, "w") as f:
    json.dump(res, f, indent=1)
    f.write("\n")
print("displaced atoms:", ", ".join(names[a] for a in displaced))
print(f"step 2s: {res['step_2s']['n_bonds']} bonds, {res['step_2s']['n_angles']} angles; every block read was computed: {ok_2s}")
print(f"step 2ms: {res['step_2ms']['n_bonds']} bonds, {res['step_2ms']['n_angles']} angles; every block read was computed: {ok_2ms}")
print("written", out_path)
