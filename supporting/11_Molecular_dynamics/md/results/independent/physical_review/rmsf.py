#!/usr/bin/env python
"""Physical review: CA fluctuation per residue of every 100 ns run (fit of the CA atoms to the 6RL9 reference, then
fluctuation about the run's own mean), set beside the B factors of 6RL9 (CA atoms of md/system/protein_md.pdb).
Also the mean displacement of each CA from its crystal position. Writes rmsf_summary.json, rmsf_table.txt."""
import os, glob, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
RUNS = "/root/cqc-sfa-md/runs"; SYSTEMS = "/root/cqc-sfa-md/systems"; OUT = os.path.dirname(os.path.abspath(__file__))
REPO = "/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate"
LIST = [("bonded_2e", r) for r in ("X1", "X2", "X3", "N1", "N2", "N3", "R1", "R2", "R3")] + [("bonded_2s", "R1"), ("bonded_2e_n2fit", "R1")]
def dep_of(resid): return resid + 1 if resid <= 124 else resid + 2
def kabsch(X, R):
    gx = X.mean(0); gr = R.mean(0)
    U, S, Vt = np.linalg.svd((X - gx).T @ (R - gr)); d = np.sign(np.linalg.det(U @ Vt))
    Q = U @ np.diag([1.0, 1.0, d]) @ Vt
    return (X - gx) @ Q + gr
# B factors of the CA atoms of 6RL9 as kept for the build
bf = {}
for line in open(f"{REPO}/md/system/protein_md.pdb"):
    if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
        bf[int(line[22:26])] = float(line[60:66])
res = {}; lines = []
for model, run in LIST:
    rd = f"{RUNS}/{model}/{run}"
    ref = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(f"{SYSTEMS}/{model}/{run[0]}/reference.pdb") if l.startswith("ATOM") and l[12:16].strip() == "CA"])
    u = mda.Universe(f"{rd}/solute.pdb", sorted(glob.glob(f"{rd}/solute_seg*.dcd")))
    ca = u.select_atoms("name CA and not resname SN1 ZN1"); assert len(ca) == 258
    deps = np.array([dep_of(a.resid) for a in ca])
    n = 0; s1 = np.zeros((258, 3)); s2 = np.zeros(258)
    for ts in u.trajectory:
        x = kabsch(ca.positions.astype(float), ref)
        s1 += x; s2 += (x * x).sum(1); n += 1
    mean = s1 / n
    rmsf = np.sqrt(np.maximum(0.0, s2 / n - (mean * mean).sum(1)))
    disp = np.sqrt(((mean - ref) ** 2).sum(1))
    b = np.array([bf[d] for d in deps]); rb = np.sqrt(3.0 * b / (8.0 * np.pi ** 2))
    core = (deps >= 6) & (deps <= 258)
    top = np.argsort(-rmsf)[:6]; topd = np.argsort(-disp)[:6]
    site = [64, 67, 92, 94, 96, 106, 119, 121, 131, 143, 198, 199, 200, 209]
    j = dict(n_frames=n, rmsf_mean=float(rmsf.mean()), rmsf_median=float(np.median(rmsf)), rmsf_max=float(rmsf.max()), rmsf_max_res=int(deps[rmsf.argmax()]),
             top_rmsf=[(int(deps[i]), float(rmsf[i])) for i in top], top_mean_displacement=[(int(deps[i]), float(disp[i])) for i in topd],
             pearson_rmsf_vs_Bfactor_rmsf=float(np.corrcoef(rmsf, rb)[0, 1]), pearson_without_termini=float(np.corrcoef(rmsf[core], rb[core])[0, 1]),
             spearman=float(np.corrcoef(np.argsort(np.argsort(rmsf)), np.argsort(np.argsort(rb)))[0, 1]),
             rmsd_of_mean_structure=float(np.sqrt((disp ** 2).mean())), site_rmsf={str(d): float(rmsf[deps == d][0]) for d in site}, site_disp={str(d): float(disp[deps == d][0]) for d in site},
             crystal_B_rmsf_mean=float(rb.mean()))
    res[f"{model}/{run}"] = j
    lines.append(f"{model}/{run}: {n} frames; CA RMSF mean {j['rmsf_mean']:.2f} A, median {j['rmsf_median']:.2f}, max {j['rmsf_max']:.2f} at residue {j['rmsf_max_res']}; six largest: " + ", ".join(f"{a} {b_:.2f}" for a, b_ in j["top_rmsf"]) +
                 f"; correlation with sqrt(3B/8pi^2) of 6RL9 (mean {j['crystal_B_rmsf_mean']:.2f} A): Pearson {j['pearson_rmsf_vs_Bfactor_rmsf']:.2f} (residues 6-258: {j['pearson_without_termini']:.2f}), Spearman {j['spearman']:.2f}; RMSD of the mean structure from the crystal {j['rmsd_of_mean_structure']:.2f} A; largest mean displacements: " + ", ".join(f"{a} {b_:.2f}" for a, b_ in j["top_mean_displacement"]))
    lines.append("   site residues RMSF / mean displacement (A): " + ", ".join(f"{d} {j['site_rmsf'][str(d)]:.2f}/{j['site_disp'][str(d)]:.2f}" for d in site))
json.dump(res, open(f"{OUT}/rmsf_summary.json", "w"), indent=1)
open(f"{OUT}/rmsf_table.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
