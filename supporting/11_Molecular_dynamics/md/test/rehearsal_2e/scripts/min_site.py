#!/usr/bin/env python
"""min_site.py RUN_DIR SYSTEM_DIR OUT.json: the zinc site and the ligand before and after the minimisation of a
REHEARSAL run on the GPU (state_min.xml against system.inpcrd), with the functions of md/build_topology.py; also
the amide hydrogens of the residues after the three zinc histidines.  Not a result for the paper."""
import json, sys, pathlib
import numpy as np
sys.path.insert(0, '/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md')
import build_topology as bt
run, sysdir, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), sys.argv[3]
sel, extra, ref_ca, ref_lig, parm, r0 = bt.system_numbers(sysdir)
X0 = parm.coordinates.copy()
res = {"what": __doc__.strip(), "run": str(run), "system": str(sysdir), "r0_A": r0, "states": {}}
lig = sel["lig_heavy"]
def amide(X):
    o = {}
    for k in (94, 96, 119):                    # tleap numbers of Phe95, Trp97, Leu120 (deposited numbers)
        r, p = parm.residues[k - 1], parm.residues[k - 2]
        n = {a.name: a.idx for a in r.atoms}; c = {a.name: a.idx for a in p.atoms}
        o[f"{r.name}{k + 1}"] = {"C-N-H_deg": round(bt.angle(X[c['C']], X[n['N']], X[n['H']]), 1),
                                 "improper_C_CA_N_H_deg": round(bt.torsion(X[c['C']], X[n['CA']], X[n['N']], X[n['H']]), 1)}
    return o
res["states"]["as built (system.inpcrd)"] = dict(bt.site_report(X0, sel, ref_ca, ref_lig, r0), amide_H=amide(X0))
for tag in ("min", "eq01", "eq02", "eq03", "eq04", "eq05", "eq06", "seg0001"):
    f = run / f"state_{tag}.xml"
    if not f.is_file():
        continue
    X = bt.state_positions_A(f)
    rep = bt.site_report(X, sel, ref_ca, ref_lig, r0)
    rep["amide_H"] = amide(X)
    if tag == "min":
        moved = np.linalg.norm(X - X0, axis=1)
        rep["moved_A"] = {"ZN": round(float(moved[sel["zn"]]), 3), **{n: round(float(moved[i]), 3) for n, i in lig.items()},
                          "HN21": round(float(moved[sel["lig_HN2"]]), 3)}
        rep["largest_move_of_zinc_or_a_ligand_heavy_atom_A"] = round(float(max(moved[[sel["zn"]] + list(lig.values())])), 3)
    res["states"][tag] = rep
json.dump(res, open(out, "w"), indent=1)
print("written", out)
