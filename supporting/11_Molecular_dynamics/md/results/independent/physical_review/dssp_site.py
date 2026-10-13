#!/usr/bin/env python
"""Physical review: secondary structure (three states, MDAnalysis' DSSP after pydssp) on one frame per ns of the
production of every 100 ns run, against the assignment of the start as built (crystal heavy atoms).
The zinc histidines are renamed HID/HIE in memory only, so that the routine sees all 258 residues.
Writes dssp_summary.json and dssp_table.txt beside itself; reads the run folders only."""
import os, glob, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
RUNS = "/root/cqc-sfa-md/runs"; OUT = os.path.dirname(os.path.abspath(__file__))
LIST = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
SITE = [62, 64, 65, 67, 92, 94, 96, 106, 117, 119, 121, 131, 141, 143, 145, 198, 199, 200, 201, 202, 207, 209, 244, 246]
def dep_of(resid): return resid + 1 if resid <= 124 else resid + 2
def prep(u):
    for r in u.residues:
        if r.resname in ("HD1", "HD2"): r.resname = "HID"
        elif r.resname == "HE1": r.resname = "HIE"
    return u.select_atoms("resid 2:259")
res = {}; lines = []
for model, run in LIST:
    rd = f"{RUNS}/{model}/{run}"
    u0 = mda.Universe(f"{rd}/solute.pdb"); ag0 = prep(u0)
    ref = DSSP(ag0).run().results.dssp[0]
    deps = np.array([dep_of(r.resid) for r in ag0.residues])
    assert len(ref) == 258
    u = mda.Universe(f"{rd}/solute.pdb", sorted(glob.glob(f"{rd}/solute_seg*.dcd"))); ag = prep(u)
    d = DSSP(ag).run(start=99, step=100).results.dssp      # frames at 1, 2, ..., 100 ns
    same = (d == ref[None, :])
    struct = ref != "-"
    keep_struct = same[:, struct].mean(0)                    # per structured residue, fraction of frames kept
    site_idx = [int(np.where(deps == s)[0][0]) for s in SITE]
    lost = [(int(deps[i]), str(ref[i]), float(same[:, i].mean())) for i in np.where(struct)[0] if same[:, i].mean() < 0.5]
    j = dict(n_frames=int(d.shape[0]), n_res=int(d.shape[1]), ref_H=int((ref == "H").sum()), ref_E=int((ref == "E").sum()),
             H_mean=float((d == "H").sum(1).mean()), H_min=int((d == "H").sum(1).min()), E_mean=float((d == "E").sum(1).mean()), E_min=int((d == "E").sum(1).min()),
             same_all_mean=float(same.mean()), same_all_min_frame=float(same.mean(1).min()), same_struct_mean=float(same[:, struct].mean()),
             same_first10=float(same[:10].mean()), same_last10=float(same[-10:].mean()),
             E_first10=float((d[:10] == "E").sum(1).mean()), E_last10=float((d[-10:] == "E").sum(1).mean()),
             site={str(s): [str(ref[i]), float(same[:, i].mean())] for s, i in zip(SITE, site_idx)}, structured_residues_kept_below_half=lost)
    res[f"{model}/{run}"] = j
    lines.append(f"{model}/{run}: {j['n_frames']} frames; start: H {j['ref_H']}, E {j['ref_E']} of 258; run: H {j['H_mean']:.1f} (min {j['H_min']}), E {j['E_mean']:.1f} (min {j['E_min']}); E first 10 ns {j['E_first10']:.1f}, last 10 ns {j['E_last10']:.1f}; "
                 f"residues with the assignment of the start: all {100*j['same_all_mean']:.1f} % (worst frame {100*j['same_all_min_frame']:.1f} %), H or E residues of the start {100*j['same_struct_mean']:.1f} %; first 10 ns {100*j['same_first10']:.1f} %, last 10 ns {100*j['same_last10']:.1f} %")
    lines.append("   site residues (start assignment, % of frames kept): " + ", ".join(f"{s} {v[0]} {100*v[1]:.0f}" for s, v in j["site"].items()))
    lines.append("   H or E residues of the start kept in under half of the frames: " + (", ".join(f"{a} {b} {100*c:.0f} %" for a, b, c in lost) if lost else "none"))
json.dump(res, open(f"{OUT}/dssp_summary.json", "w"), indent=1)
open(f"{OUT}/dssp_table.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
