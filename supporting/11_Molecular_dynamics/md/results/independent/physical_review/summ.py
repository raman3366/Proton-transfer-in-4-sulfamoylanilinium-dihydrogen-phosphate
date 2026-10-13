#!/usr/bin/env python
"""Physical review: tables over the eleven 100 ns runs from the npz files of extract.py.
Writes tables.txt and tables.json beside itself."""
import os, json, sys
import numpy as np
from common import *
OUT = HERE
cry = json.load(open(f"{HERE}/crystal_values.json"))
runs = [Run(m, r) for m, r in RUNLIST if os.path.exists(f"{HERE}/data_{m}_{r}.npz")]
J = {}
L = []
def P(*a):
    L.append(" ".join(str(x) for x in a))
ST = "CRI"
ANG = ["N2-H94", "N2-H96", "N2-H119", "H94-H96", "H94-H119", "H96-H119"]

# ---------------------------------------------------------------- 1 integrity
P("== 1. INTEGRITY (production frames, 10 ps; full frames 100 ps) ==")
P("run | solute frames, full frames | CA RMSD A: 0-10 ns, 45-55 ns, 90-100 ns, max, largest step between frames, slope over the last 50 ns (A per 100 ns) | Rg heavy atoms A: mean +- sd (min to max) [start X %.3f] | largest CA(i)-CA(i+1) A | frames in which the zinc / the ligand had to be moved by a lattice vector | protein heavy atoms to nearest image A: min" % cry["X"]["rg_heavy"])
for r in runs:
    s = r.s; p = r.prod; t = r.t[p] / 1000.0
    ca = s["ca_rmsd"][p]
    seg = lambda a, b: float(ca[(t > a) & (t <= b)].mean())
    sl = float(np.polyfit(t[t > 50], ca[t > 50], 1)[0] * 100)
    rg = s["rg_heavy"][p]
    img = r.f["prot_image_min"][r.prodf]
    j = dict(n_solute=int(p.sum()), n_full=int(r.prodf.sum()), ca_0_10=seg(0, 10), ca_45_55=seg(45, 55), ca_90_100=seg(90, 100), ca_max=float(ca.max()),
             ca_step_max_between_frames=float(np.abs(np.diff(ca)).max()), ca_slope_last50_A_per_100ns=sl, rg_mean=float(rg.mean()), rg_sd=float(rg.std()), rg_min=float(rg.min()), rg_max=float(rg.max()),
             rg_step_max=float(np.abs(np.diff(rg)).max()),
             ca_neighbour_max=float(max(s["ca_step_max"].max(), r.f["ca_step_max"].max())), zn_shifted=int(s["zn_shifted"].sum()), lig_shifted=int(s["lig_shifted"].sum()),
             raw_Zn_N2_max=float(s["raw_Zn_N2"].max()), raw_Zn_H94_max=float(s["raw_Zn_H94NE2"].max()), image_min=float(np.nanmin(img)), image_mean=float(np.nanmean(img)),
             all_frames_solute=int(len(r.t)), all_frames_full=int(len(r.tf)), frames=r.meta["frames"], sel_mismatch=r.meta["sel_mismatch"],
             box_a_mean=float(s["box_a"][p].mean()), box_a_min=float(s["box_a"][p].min()), box_a_max=float(s["box_a"][p].max()))
    J.setdefault(r.name, {})["integrity"] = j
    P(f"{r.name} | {j['n_solute']}, {j['n_full']} | {j['ca_0_10']:.2f}, {j['ca_45_55']:.2f}, {j['ca_90_100']:.2f}, {j['ca_max']:.2f}, {j['ca_step_max_between_frames']:.2f}, {sl:+.2f} | {j['rg_mean']:.3f} +- {j['rg_sd']:.3f} ({j['rg_min']:.3f} to {j['rg_max']:.3f}) | {j['ca_neighbour_max']:.2f} | {j['zn_shifted']} / {j['lig_shifted']} (raw Zn-N2 max {j['raw_Zn_N2_max']:.2f}) | {j['image_min']:.1f}")
P("")
P("ligand bonds and ring, production, mean (min to max), A")
P("run | S-N2 | S-C4 | S-O1 | S-O2 | C1-N1 | ring C-C min, max | N2-H | ring r.m.s. from plane max | S out of ring plane mean, max | C4-S-N2 deg | O1-S-O2 deg | cis peptides |omega| mean, max: Ser29-Pro30; Pro201-Pro202")
for r in runs:
    s = r.s; p = r.prod
    f = lambda k: "%.3f (%.3f to %.3f)" % (s[k][p].mean(), s[k][p].min(), s[k][p].max())
    o1 = np.abs(s["omega_S29_P30"][p]); o2 = np.abs(s["omega_P201_P202"][p])
    J[r.name]["ligand_internal"] = {k: [float(s[k][p].mean()), float(s[k][p].min()), float(s[k][p].max())] for k in ("d_S_N2", "d_S_C4", "d_S_O1", "d_S_O2", "d_C1_N1", "d_N2_H", "ring_cc_min", "ring_cc_max", "ring_rms_plane", "oop_S", "oop_N1", "ang_C4_S_N2", "ang_O1_S_O2", "ang_O1_S_C4", "ang_O2_S_C4")}
    J[r.name]["cis_peptides"] = dict(S29_P30_mean=float(o1.mean()), S29_P30_max=float(o1.max()), P201_P202_mean=float(o2.mean()), P201_P202_max=float(o2.max()))
    P(f"{r.name} | {f('d_S_N2')} | {f('d_S_C4')} | {f('d_S_O1')} | {f('d_S_O2')} | {f('d_C1_N1')} | {s['ring_cc_min'][p].min():.3f}, {s['ring_cc_max'][p].max():.3f} | {f('d_N2_H')} | {s['ring_rms_plane'][p].max():.3f} | {s['oop_S'][p].mean():.2f}, {s['oop_S'][p].max():.2f} | {s['ang_C4_S_N2'][p].mean():.1f} | {s['ang_O1_S_O2'][p].mean():.1f} | {o1.mean():.1f}, {o1.max():.1f}; {o2.mean():.1f}, {o2.max():.1f}")

# ---------------------------------------------------------------- 2 zinc site
P("")
P("== 2. ZINC SITE (production) ==")
P("crystal (start X as built): Zn-N2 %.3f, Zn-His94 %.3f, Zn-His96 %.3f, Zn-His119 %.3f; angles %s; out of plane His94 %.2f, His96 %.2f, His119 %.2f; Zn-O1 %.3f, Zn-O2 %.3f; Zn-Glu106 OE1 %.3f, OE2 %.3f" % (
    cry["X"]["d_Zn_N2"], cry["X"]["d_Zn_H94"], cry["X"]["d_Zn_H96"], cry["X"]["d_Zn_H119"],
    ", ".join("%s %.1f" % (a, cry["X"]["angZn_%d_%s" % (k, a)]) for k, a in enumerate(ANG)), cry["X"]["zn_oop_94"], cry["X"]["zn_oop_96"], cry["X"]["zn_oop_119"],
    cry["X"]["d_Zn_O1"], cry["X"]["d_Zn_O2"], cry["X"]["d_Zn_E106OE1"], cry["X"]["d_Zn_E106OE2"]))
P("run | Zn-N2 | Zn-His94 NE2 | Zn-His96 NE2 | Zn-His119 ND1   [mean +- sd (min to max), A]")
for r in runs:
    s = r.s; p = r.prod
    J[r.name]["zinc"] = {}
    row = []
    for k in ("d_Zn_N2", "d_Zn_H94", "d_Zn_H96", "d_Zn_H119"):
        a = s[k][p]; J[r.name]["zinc"][k] = [float(a.mean()), float(a.std()), float(a.min()), float(a.max())]
        row.append(ms(a))
    P(r.name, "|", " | ".join(row))
P("run | six N-Zn-N angles, mean +- sd (min to max), deg: " + ", ".join(ANG))
for r in runs:
    s = r.s; p = r.prod; row = []
    for k in range(6):
        a = s[f"angZn_{k}"][p]; J[r.name]["zinc"]["ang_" + ANG[k]] = [float(a.mean()), float(a.std()), float(a.min()), float(a.max())]
        row.append("%.1f +- %.1f (%.1f to %.1f)" % (a.mean(), a.std(), a.min(), a.max()))
    P(r.name, "|", " | ".join(row))
P("run | zinc out of the imidazole plane, signed, mean +- sd (min to max) A: His94 | His96 | His119 | frames with |out of plane| above 1.0 A: 94, 96, 119")
for r in runs:
    s = r.s; p = r.prod; row = []; cnt = []
    for h in (94, 96, 119):
        a = s[f"zn_oop_{h}"][p]; J[r.name]["zinc"][f"oop_{h}"] = [float(a.mean()), float(a.std()), float(a.min()), float(a.max()), int((np.abs(a) > 1.0).sum())]
        row.append("%.2f +- %.2f (%.2f to %.2f)" % (a.mean(), a.std(), a.min(), a.max())); cnt.append(str(int((np.abs(a) > 1.0).sum())))
    P(r.name, "|", " | ".join(row), "|", ", ".join(cnt))
P("run | Zn-O1 mean (min) | Zn-O2 mean (min) | nearer sulfonyl O: mean; frames <= 2.6 A; <= 3.0 A (of the solute frames) | water O nearest to Zn: mean, min; full frames <= 2.6, <= 3.0, <= 3.5 A | protein O nearest to Zn: mean, min, which atom at the minimum; frames <= 2.6, <= 3.0 | zinc displacement from its crystal position after the CA fit, mean (max)")
for r in runs:
    s = r.s; p = r.prod; f = r.f; pf = r.prodf
    omin = np.minimum(s["d_Zn_O1"][p], s["d_Zn_O2"][p])
    dw = f["d_Zn_Ow_min"][pf]; dp = f["d_Zn_Oprot_min"][pf]
    kp = int(np.argmin(dp)); lab = r.meta["prot_O_label"][int(f["i_Zn_Oprot_min"][pf][kp])]
    # which protein oxygen is the nearest, over frames
    ids, cnts = np.unique(f["i_Zn_Oprot_min"][pf], return_counts=True)
    who = ", ".join("%s %.0f%%" % (r.meta["prot_O_label"][int(i)], 100.0 * c / pf.sum()) for i, c in sorted(zip(ids, cnts), key=lambda z: -z[1])[:3])
    j = dict(ZnO1_mean=float(s["d_Zn_O1"][p].mean()), ZnO1_min=float(s["d_Zn_O1"][p].min()), ZnO2_mean=float(s["d_Zn_O2"][p].mean()), ZnO2_min=float(s["d_Zn_O2"][p].min()),
             Osulf_min_mean=float(omin.mean()), Osulf_le_2p6=int((omin <= 2.6).sum()), Osulf_le_3p0=int((omin <= 3.0).sum()), n_solute=int(p.sum()),
             Ow_mean=float(dw.mean()), Ow_min=float(dw.min()), Ow_le_2p6=int((dw <= 2.6).sum()), Ow_le_3p0=int((dw <= 3.0).sum()), Ow_le_3p5=int((dw <= 3.5).sum()),
             Oprot_mean=float(dp.mean()), Oprot_min=float(dp.min()), Oprot_min_atom=lab, Oprot_le_2p6=int((dp <= 2.6).sum()), Oprot_le_3p0=int((dp <= 3.0).sum()), Oprot_nearest_who=who,
             n_full=int(pf.sum()), zn_disp_mean=float(s["zn_disp_fit"][p].mean()), zn_disp_max=float(s["zn_disp_fit"][p].max()),
             n_w_5A_mean=float(f["n_Ow_Zn_5p0"][pf].mean()), n_w_8A_mean=float(f["n_Ow_Zn_8p0"][pf].mean()))
    J[r.name]["zinc"].update(j)
    P(f"{r.name} | {j['ZnO1_mean']:.2f} ({j['ZnO1_min']:.2f}) | {j['ZnO2_mean']:.2f} ({j['ZnO2_min']:.2f}) | {j['Osulf_min_mean']:.3f}; {j['Osulf_le_2p6']} ({100*j['Osulf_le_2p6']/j['n_solute']:.2f} %); {j['Osulf_le_3p0']} ({100*j['Osulf_le_3p0']/j['n_solute']:.1f} %) | {j['Ow_mean']:.2f}, {j['Ow_min']:.2f}; {j['Ow_le_2p6']}, {j['Ow_le_3p0']}, {j['Ow_le_3p5']} of {j['n_full']} | {j['Oprot_mean']:.2f}, {j['Oprot_min']:.2f}, {lab}; {j['Oprot_le_2p6']}, {j['Oprot_le_3p0']} [nearest: {who}] | {j['zn_disp_mean']:.2f} ({j['zn_disp_max']:.2f})")
P("")
P("crystal: Thr199 OG1...Glu106 OE1 %.2f A, angle at HG1 %.0f (hydrogen as built); Glu106 OE2...Arg246 N %.2f; His94 ND1...Gln92 OE1 %.2f; His96 ND1...Asn244 O %.2f; His119 NE2...Glu117 OE %.2f; His64 chi1 %.0f chi2 %.0f, Zn...ND1 %.1f, Zn...NE2 %.1f" % (
    cry["X"]["d_T199OG1_E106OE1"], cry["X"]["a_T199OG1H_E106OE1"], cry["X"]["d_E106OE2_R246N"], cry["X"]["d_H94ND1_Q92OE1"], cry["X"]["d_H96ND1_N244O"], cry["X"]["d_H119NE2_E117OE"], cry["X"]["h64_chi1"], cry["X"]["h64_chi2"], cry["X"]["d_Zn_H64ND1"], cry["X"]["d_Zn_H64NE2"]))
P("run | Thr199 OG1-H...Glu106 OE (either O; d<=3.5, angle>=135): %; to OE1 %; OG1...OE1 mean A | Thr199 HG1 nearest approach to N2 / O1 / O2 A | Zn...Glu106 OE1 mean (min); OE2 mean (min) | Arg246 N-H...Glu106 OE % | His94 ND1-H...Gln92 OE1 % (d mean) | His96 ND1-H...Asn244 O % (d mean) | His119 NE2-H...Glu117 OE % (d mean)")
for r in runs:
    s = r.s; p = r.prod
    h1 = hb(s["d_T199OG1_E106OE1"][p], s["a_T199OG1H_E106OE1"][p]); h2 = hb(s["d_T199OG1_E106OE2"][p], s["a_T199OG1H_E106OE2"][p])
    hr = hb(s["d_E106OE2_R246N"][p], s["a_R246NH_E106OE2"][p]) | hb(s["d_E106OE1_R246N"][p], s["a_R246NH_E106OE1"][p])
    h94 = hb(s["d_H94ND1_Q92OE1"][p], s["a_H94ND1H_Q92OE1"][p]); h96 = hb(s["d_H96ND1_N244O"][p], s["a_H96ND1H_N244O"][p]); h119 = hb(s["d_H119NE2_E117OE"][p], s["a_H119NE2H_E117OE"][p])
    j = dict(T199_E106_either=float((h1 | h2).mean()), T199_E106_OE1=float(h1.mean()), T199_E106_OE2=float(h2.mean()), d_OG1_OE1_mean=float(s["d_T199OG1_E106OE1"][p].mean()), d_OG1_OE1_max=float(s["d_T199OG1_E106OE1"][p].max()),
             HG1_N2_min=float(s["d_T199HG1_N2"][p].min()), HG1_O1_min=float(s["d_T199HG1_O1"][p].min()), HG1_O2_min=float(s["d_T199HG1_O2"][p].min()),
             Zn_OE1_mean=float(s["d_Zn_E106OE1"][p].mean()), Zn_OE1_min=float(s["d_Zn_E106OE1"][p].min()), Zn_OE2_mean=float(s["d_Zn_E106OE2"][p].mean()), Zn_OE2_min=float(s["d_Zn_E106OE2"][p].min()),
             R246_E106=float(hr.mean()), H94_Q92=float(h94.mean()), d_H94_Q92=float(s["d_H94ND1_Q92OE1"][p].mean()), H96_N244=float(h96.mean()), d_H96_N244=float(s["d_H96ND1_N244O"][p].mean()),
             H119_E117=float(h119.mean()), d_H119_E117=float(s["d_H119NE2_E117OE"][p].mean()), d_H94_Q92NE2=float(s["d_H94ND1_Q92NE2"][p].mean()))
    J[r.name]["second_shell"] = j
    P(f"{r.name} | {100*j['T199_E106_either']:.1f}; {100*j['T199_E106_OE1']:.1f}; {j['d_OG1_OE1_mean']:.2f} | {j['HG1_N2_min']:.2f} / {j['HG1_O1_min']:.2f} / {j['HG1_O2_min']:.2f} | {j['Zn_OE1_mean']:.2f} ({j['Zn_OE1_min']:.2f}); {j['Zn_OE2_mean']:.2f} ({j['Zn_OE2_min']:.2f}) | {100*j['R246_E106']:.1f} | {100*j['H94_Q92']:.1f} ({j['d_H94_Q92']:.2f}) | {100*j['H96_N244']:.1f} ({j['d_H96_N244']:.2f}) | {100*j['H119_E117']:.1f} ({j['d_H119_E117']:.2f})")
P("")
P("His64 (built 'out': chi1 %.0f). classes by chi1: g+ (0 to 120, 'in'-like), g- (-120 to 0, 'out'-like), t (else)" % cry["X"]["h64_chi1"])
P("run | % g+ / g- / t | Zn...NE2 mean A in g+ ; in g- | Zn...ND1 mean in g+ ; in g- | number of changes of class that last >= 100 ps | first time in g+ (ns) | class in the last 20 ns: % g+ | ring of His64 to aniline N1: min A, % frames <= 4.0")
for r in runs:
    s = r.s; p = r.prod; c1 = s["h64_chi1"][p]; t = r.t[p] / 1000
    cl = np.where((c1 > 0) & (c1 <= 120), 0, np.where((c1 < 0) & (c1 >= -120), 1, 2))
    v = [x for x in visits(cl) if x[2] >= 10]
    nchg = sum(1 for a, b in zip(v[:-1], v[1:]) if a[0] != b[0])
    gp = cl == 0; gm = cl == 1
    first = float(t[np.argmax(gp)]) if gp.any() else None
    j = dict(frac_gplus=float(gp.mean()), frac_gminus=float(gm.mean()), frac_t=float((cl == 2).mean()),
             ZnNE2_gplus=float(s["d_Zn_H64NE2"][p][gp].mean()) if gp.any() else None, ZnNE2_gminus=float(s["d_Zn_H64NE2"][p][gm].mean()) if gm.any() else None,
             ZnND1_gplus=float(s["d_Zn_H64ND1"][p][gp].mean()) if gp.any() else None, ZnND1_gminus=float(s["d_Zn_H64ND1"][p][gm].mean()) if gm.any() else None,
             n_changes_100ps=nchg, first_gplus_ns=first, last20_gplus=float(gp[t > 80].mean()), N1_ring_min=float(s["d_N1_H64ring"][p].min()), N1_ring_le4=float((s["d_N1_H64ring"][p] <= 4.0).mean()),
             ZnNE2_min=float(s["d_Zn_H64NE2"][p].min()), ZnND1_min=float(s["d_Zn_H64ND1"][p].min()))
    J[r.name]["his64"] = j
    fm = lambda x: "none" if x is None else "%.1f" % x
    P(f"{r.name} | {100*j['frac_gplus']:.1f} / {100*j['frac_gminus']:.1f} / {100*j['frac_t']:.1f} | {fm(j['ZnNE2_gplus'])} ; {fm(j['ZnNE2_gminus'])} | {fm(j['ZnND1_gplus'])} ; {fm(j['ZnND1_gminus'])} | {nchg} | {fm(first) if first is None else '%.2f' % first} | {100*j['last20_gplus']:.1f} | {j['N1_ring_min']:.2f}, {100*j['N1_ring_le4']:.1f}")

# ---------------------------------------------------------------- 3 ligand
P("")
P("== 3. LIGAND (production) ==")
P("crystal/start X: N2...OG1 %.2f, angle %.0f; Thr199 N...O1 %.2f (angle %.0f), N...O2 %.2f; V %+.3f; Zn-N2-S %.1f, Zn-N2-H %.1f, S-N2-H %.1f; Zn-N2-S-C4 %.1f; H-N2-S-C4 %.1f; C3-C4-S-N2 %.1f; His94NE2-Zn-N2-S %.1f, -H %.1f; angle (Zn-N2, C4->C1) %.1f" % (
    cry["X"]["d_N2_T199OG1"], cry["X"]["a_N2H_T199OG1"], cry["X"]["d_T199N_O1"], cry["X"]["a_T199NH_O1"], cry["X"]["d_T199N_O2"], cry["X"]["n2_vol"], cry["X"]["a_Zn_N2_S"], cry["X"]["a_Zn_N2_H"], cry["X"]["a_S_N2_H"],
    cry["X"]["t_Zn_N2_S_C4"], cry["X"]["t_H_N2_S_C4"], cry["X"]["t_C3_C4_S_N2"], cry["X"]["t_H94_Zn_N2_S"], cry["X"]["t_H94_Zn_N2_H"], cry["X"]["a_ZnN2_axis"]))
P("start R as built: V %+.3f; Zn-N2-S-C4 %.1f; His94NE2-Zn-N2-S %.1f, -H %.1f; N2...OG1 %.2f angle %.0f; angle (Zn-N2, C4->C1) %.1f" % (cry["R"]["n2_vol"], cry["R"]["t_Zn_N2_S_C4"], cry["R"]["t_H94_Zn_N2_S"], cry["R"]["t_H94_Zn_N2_H"], cry["R"]["d_N2_T199OG1"], cry["R"]["a_N2H_T199OG1"], cry["R"]["a_ZnN2_axis"]))
P("run | Zn-N2 max A | ligand centre to nearest protein heavy atom max A | N2-H...OG1(Thr199) % (N...O mean, max) | Thr199 N-H...O1 %, ...O2 %, either % | Thr199 N...nearer sulfonyl O mean A | ligand RMSD to the crystal pose mean +- sd (max)")
for r in runs:
    s = r.s; p = r.prod
    h = hb(s["d_N2_T199OG1"][p], s["a_N2H_T199OG1"][p]); h1 = hb(s["d_T199N_O1"][p], s["a_T199NH_O1"][p]); h2 = hb(s["d_T199N_O2"][p], s["a_T199NH_O2"][p])
    dn = np.minimum(s["d_T199N_O1"][p], s["d_T199N_O2"][p])
    j = dict(ZnN2_max=float(s["d_Zn_N2"][p].max()), ligcen_prot_max=float(s["d_ligcen_prot_min"][p].max()), hb_N2H_OG1=float(h.mean()), d_N2_OG1_mean=float(s["d_N2_T199OG1"][p].mean()), d_N2_OG1_max=float(s["d_N2_T199OG1"][p].max()),
             hb_T199NH_O1=float(h1.mean()), hb_T199NH_O2=float(h2.mean()), hb_T199NH_either=float((h1 | h2).mean()), d_T199N_O_mean=float(dn.mean()),
             rmsd_mean=float(s["lig_rmsd_X"][p].mean()), rmsd_sd=float(s["lig_rmsd_X"][p].std()), rmsd_max=float(s["lig_rmsd_X"][p].max()),
             rmsd_last_half=float(s["lig_rmsd_X"][p][r.t[p] > 50000].mean()),
             hb_T199NH_N2=float(hb(s["d_T199N_N2"][p], s["a_T199NH_N2"][p]).mean()))
    J[r.name]["ligand"] = j
    P(f"{r.name} | {j['ZnN2_max']:.2f} | {j['ligcen_prot_max']:.2f} | {100*j['hb_N2H_OG1']:.1f} ({j['d_N2_OG1_mean']:.2f}, {j['d_N2_OG1_max']:.2f}) | {100*j['hb_T199NH_O1']:.1f}, {100*j['hb_T199NH_O2']:.1f}, {100*j['hb_T199NH_either']:.1f} | {j['d_T199N_O_mean']:.2f} | {j['rmsd_mean']:.2f} +- {j['rmsd_sd']:.2f} ({j['rmsd_max']:.2f})")

# ---------------------------------------------------------------- orientation and N2
P("")
P("== ORIENTATION AND THE CONFIGURATION AT N2 ==")
P("run | production frames C / R / I % | last 20 ns C / R / I % and verdict by the rule | V > 0 % (production) | |V| < 0.2 %, |V| < 0.1 % | V mean +- sd | sum of the three angles at N2 mean (max) | changes of the sign of V between consecutive frames | first long visit (>= 100 frames) in the opposite state on the joined series: ps from release (and stage) | frames from release")
for r in runs:
    s = r.s; p = r.prod; st = s["state"][p]; V = s["n2_vol"][p]
    last = st[-2000:]
    cC, cR, cI = [(last == k).sum() for k in (0, 1, 2)]
    start_state = 1 if r.start == "R" else 0; opp = 1 - start_state
    c_opp = (last == opp).sum(); c_st = (last == start_state).sum()
    verdict = "CONVERTED" if 10 * c_opp >= 9 * len(last) else ("STAYED" if 10 * c_st >= 9 * len(last) else "MIXED")
    fr = r.free
    a = assigned(s["state"][fr], start_state); tt = r.t[fr]
    fp = None; fp_stage = None
    for val, i0, n in visits(a):
        if val == opp and n >= 100:
            fp = float(tt[i0] - r.trel); fp_stage = str(s["stage"][fr][i0]); break
    vis = [(ST[v], float(tt[i0]), n) for v, i0, n in visits(a) if n >= 100]
    sang = s["a_Zn_N2_S"][p] + s["a_Zn_N2_H"][p] + s["a_S_N2_H"][p]
    Vfree = s["n2_vol"][fr]
    j = dict(frac_C=float((st == 0).mean()), frac_R=float((st == 1).mean()), frac_I=float((st == 2).mean()), last20=[int(cC), int(cR), int(cI)], verdict=verdict,
             V_pos=float((V > 0).mean()), V_abs_lt_0p2=float((np.abs(V) < 0.2).mean()), V_abs_lt_0p1=float((np.abs(V) < 0.1).mean()), V_abs_lt_0p3=float((np.abs(V) < 0.3).mean()), V_mean=float(V.mean()), V_sd=float(V.std()), V_min=float(V.min()), V_max=float(V.max()),
             sumang_mean=float(sang.mean()), sumang_max=float(sang.max()), sign_changes=int((np.sign(V[1:]) != np.sign(V[:-1])).sum()),
             sign_changes_from_release=int((np.sign(Vfree[1:]) != np.sign(Vfree[:-1])).sum()),
             first_passage_ps_from_release=fp, first_passage_stage=fp_stage, long_visits=vis, n_free=int(fr.sum()),
             aZnN2S=float(s["a_Zn_N2_S"][p].mean()), aZnN2H=float(s["a_Zn_N2_H"][p].mean()), aSN2H=float(s["a_S_N2_H"][p].mean()))
    J[r.name]["orientation"] = j
    P(f"{r.name} | {100*j['frac_C']:.2f} / {100*j['frac_R']:.2f} / {100*j['frac_I']:.2f} | {cC/20:.1f} / {cR/20:.1f} / {cI/20:.1f} {verdict} | {100*j['V_pos']:.2f} | {100*j['V_abs_lt_0p2']:.2f}, {100*j['V_abs_lt_0p1']:.2f} | {j['V_mean']:+.3f} +- {j['V_sd']:.3f} ({j['V_min']:+.2f} to {j['V_max']:+.2f}) | {j['sumang_mean']:.1f} ({j['sumang_max']:.1f}) | {j['sign_changes']} (from release {j['sign_changes_from_release']}) | {fp} ({fp_stage}) | {j['n_free']}")
P("")
P("angles at N2, production means (deg): run | Zn-N2-S | Zn-N2-H | S-N2-H | Zn-N2-S-C4 circular mean +- sd | H-N2-S-C4 | His94NE2-Zn-N2-S | His94NE2-Zn-N2-H | tau_fold mean")
for r in runs:
    s = r.s; p = r.prod
    tf = np.maximum(np.abs(s["t_C3_C4_S_N2"][p]), np.abs(s["t_C5_C4_S_N2"][p]))
    j = dict(t_Zn_N2_S_C4=[circ_mean(s["t_Zn_N2_S_C4"][p]), circ_sd(s["t_Zn_N2_S_C4"][p])], t_H_N2_S_C4=[circ_mean(s["t_H_N2_S_C4"][p]), circ_sd(s["t_H_N2_S_C4"][p])],
             t_H94_Zn_N2_S=[circ_mean(s["t_H94_Zn_N2_S"][p]), circ_sd(s["t_H94_Zn_N2_S"][p])], t_H94_Zn_N2_H=[circ_mean(s["t_H94_Zn_N2_H"][p]), circ_sd(s["t_H94_Zn_N2_H"][p])], tau_fold=float(tf.mean()))
    J[r.name]["torsions"] = j
    o = J[r.name]["orientation"]
    P(f"{r.name} | {o['aZnN2S']:.1f} | {o['aZnN2H']:.1f} | {o['aSN2H']:.1f} | {j['t_Zn_N2_S_C4'][0]:.1f} +- {j['t_Zn_N2_S_C4'][1]:.1f} | {j['t_H_N2_S_C4'][0]:.1f} +- {j['t_H_N2_S_C4'][1]:.1f} | {j['t_H94_Zn_N2_S'][0]:.1f} +- {j['t_H94_Zn_N2_S'][1]:.1f} | {j['t_H94_Zn_N2_H'][0]:.1f} +- {j['t_H94_Zn_N2_H'][1]:.1f} | {j['tau_fold']:.1f}")

# ---------------------------------------------------------------- per state: tilt, partners, water
P("")
P("== LIGAND BY STATE (frames from the release of the ligand; C crystal-like, R rotated-like, I between) ==")
P("run state | frames | ligand RMSD to crystal pose | to the rotated start | ring axis C4->C1 against the crystal axis deg | ring normal against the crystal normal deg | angle between Zn-N2 and the ring axis | V > 0 % | N2-H...OG1 % | Thr199 N-H...O % | Zn-N2-S-C4 circ. mean | nearer sulfonyl O to Zn")
for r in runs:
    s = r.s; fr = r.free
    J[r.name]["by_state"] = {}
    for k in (0, 1, 2):
        m = fr & (s["state"] == k)
        if m.sum() == 0: continue
        h = hb(s["d_N2_T199OG1"][m], s["a_N2H_T199OG1"][m]); h12 = hb(s["d_T199N_O1"][m], s["a_T199NH_O1"][m]) | hb(s["d_T199N_O2"][m], s["a_T199NH_O2"][m])
        j = dict(n=int(m.sum()), rmsd_X=float(s["lig_rmsd_X"][m].mean()), rmsd_X_sd=float(s["lig_rmsd_X"][m].std()), rmsd_R=float(s["lig_rmsd_R"][m].mean()), axis_to_X=float(s["axis_to_X"][m].mean()), axis_to_X_sd=float(s["axis_to_X"][m].std()),
                 normal_to_X=float(s["normal_to_X"][m].mean()), a_ZnN2_axis=float(s["a_ZnN2_axis"][m].mean()), a_ZnN2_axis_sd=float(s["a_ZnN2_axis"][m].std()), V_pos=float((s["n2_vol"][m] > 0).mean()), hb_N2H=float(h.mean()), hb_T199NH=float(h12.mean()),
                 t_ZNSC=circ_mean(s["t_Zn_N2_S_C4"][m]), Osulf=float(np.minimum(s["d_Zn_O1"][m], s["d_Zn_O2"][m]).mean()), axis_to_R=float(s["axis_to_R"][m].mean()))
        J[r.name]["by_state"][ST[k]] = j
        P(f"{r.name} {ST[k]} | {j['n']} | {j['rmsd_X']:.2f} +- {j['rmsd_X_sd']:.2f} | {j['rmsd_R']:.2f} | {j['axis_to_X']:.1f} +- {j['axis_to_X_sd']:.1f} | {j['normal_to_X']:.1f} | {j['a_ZnN2_axis']:.1f} +- {j['a_ZnN2_axis_sd']:.1f} | {100*j['V_pos']:.1f} | {100*j['hb_N2H']:.1f} | {100*j['hb_T199NH']:.1f} | {j['t_ZNSC']:.1f} | {j['Osulf']:.2f}")
P("")
P("partners of the aniline N1 among protein N and O within 3.5 A, by state (frames from release); 'd' = N1-H...acceptor angle >= 135 (N1 donates), 'a' = donor-H...N1 angle >= 135 (N1 accepts), plain = contact only")
for r in runs:
    s = r.s; fr = r.free; pl = r.meta["polar_label"]
    J[r.name]["N1_partners"] = {}
    for k in (0, 1, 2):
        m = np.where(fr & (s["state"] == k))[0]
        if len(m) == 0: continue
        pi = s["partner_idx"][m]; pc = s["partner_code"][m]
        tally = {}
        for q in range(pi.shape[1]):
            ok = pi[:, q] >= 0
            for i, c in zip(pi[ok, q], pc[ok, q]):
                e = tally.setdefault(pl[int(i)], [0, 0, 0]); e[0] += 1
                if c & 1: e[1] += 1
                if c & 2: e[2] += 1
        none = float((pi[:, 0] < 0).mean())
        top = sorted(tally.items(), key=lambda z: -z[1][0])[:6]
        J[r.name]["N1_partners"][ST[k]] = dict(n=int(len(m)), no_protein_partner=none, top={a: [b[0] / len(m), b[1] / len(m), b[2] / len(m)] for a, b in top})
        P(f"{r.name} {ST[k]} ({len(m)} frames; no protein N/O within 3.5 A in {100*none:.1f} %): " + "; ".join(f"{a} {100*b[0]/len(m):.1f} % (d {100*b[1]/len(m):.1f}, a {100*b[2]/len(m):.1f})" for a, b in top))
P("")
P("water at the ligand in the all-atom frames (from release), by state of the solute frame of the same time: run state | full frames | water O within 3.5 A of N1: mean number; of them N1 donates to / accepts from (mean number) | water within 3.5 A of N2 | of O1 | of O2")
for r in runs:
    f = r.f; ff = r.freef
    J[r.name]["water_at_ligand"] = {}
    for k in (0, 1, 2):
        m = ff & (r.statef == k)
        if m.sum() == 0: continue
        j = dict(n=int(m.sum()), N1=float(f["n_w_N1_3p5"][m].mean()), N1_donates=float(f["n_w_N1_donated"][m].mean()), N1_accepts=float(f["n_w_N1_accepted"][m].mean()), N1_none=float((f["n_w_N1_3p5"][m] == 0).mean()),
                 N2=float(f["n_w_N2_3p5"][m].mean()), O1=float(f["n_w_O1_3p5"][m].mean()), O2=float(f["n_w_O2_3p5"][m].mean()))
        J[r.name]["water_at_ligand"][ST[k]] = j
        P(f"{r.name} {ST[k]} | {j['n']} | {j['N1']:.2f}; {j['N1_donates']:.2f} / {j['N1_accepts']:.2f} (no water at N1 in {100*j['N1_none']:.1f} %) | {j['N2']:.2f} | {j['O1']:.2f} | {j['O2']:.2f}")
P("")
P("contacts of the ring atoms (C1..C6, N1) with residues: mean of the smallest heavy-atom distance (A), by state, frames from release")
P("residues: " + ", ".join(str(x) for x in runs[0].meta["contact_res"]))
for r in runs:
    s = r.s; fr = r.free
    J[r.name]["ring_contacts"] = {}
    for k in (0, 1):
        m = fr & (s["state"] == k)
        if m.sum() < 20: continue
        c = s["contacts"][m].mean(0)
        J[r.name]["ring_contacts"][ST[k]] = [float(x) for x in c]
        P(f"{r.name} {ST[k]}: " + " ".join("%.1f" % x for x in c))

# ---------------------------------------------------------------- water sites, ions
P("")
P("== WATER AT THE FOUR CRYSTAL-WATER POSITIONS 423, 516, 553, 568 (place of the rotated ring) AND AT THE TWO RING PLACES; all-atom frames ==")
P("occupied = a water oxygen within 1.4 A of the crystal position after the CA fit (1.0 A in brackets)")
P("run window | full frames | occupancy % of W423, W516, W553, W568 | all four occupied % | mean number of water O within 2.5 A of the crystal ring carbons (place of the crystal ring); frames with at least one % | the same for the ring carbons of the rotated start (place of the rotated ring)")
sn = runs[0].meta["site_names"]; wi = [sn.index(f"W{w}") for w in (423, 516, 553, 568)]
def water_row(r, m, tag):
    f = r.f
    d = f["site_dmin"][m][:, wi]
    occ = (d <= 1.4).mean(0); occ1 = (d <= 1.0).mean(0)
    j = dict(n=int(m.sum()), occ_1p4=[float(x) for x in occ], occ_1p0=[float(x) for x in occ1], all4=float((d <= 1.4).all(1).mean()), n_any=float(((d <= 1.4).sum(1)).mean()),
             Xring_n=float(f["n_w_Xring_2p5"][m].mean()), Xring_any=float((f["n_w_Xring_2p5"][m] > 0).mean()), Rring_n=float(f["n_w_Rring_2p5"][m].mean()), Rring_any=float((f["n_w_Rring_2p5"][m] > 0).mean()))
    J[r.name].setdefault("water_sites", {})[tag] = j
    P(f"{r.name} {tag} | {j['n']} | " + ", ".join("%.0f (%.0f)" % (100 * a, 100 * b) for a, b in zip(occ, occ1)) + f" | {100*j['all4']:.0f} | {j['Xring_n']:.2f}; {100*j['Xring_any']:.1f} | {j['Rring_n']:.2f}; {100*j['Rring_any']:.1f}")
for r in runs:
    tf = r.tf
    water_row(r, tf <= r.trel, "held part (before release)")
    water_row(r, (tf > r.trel) & (tf <= 0), "equilibration after release")
    water_row(r, (tf > 0) & (tf <= 20000), "0-20 ns")
    water_row(r, (tf > 20000) & (tf <= 80000), "20-80 ns")
    water_row(r, tf > 80000, "80-100 ns")
    for k in (0, 1):
        m = r.prodf & (r.statef == k)
        if m.sum() >= 5: water_row(r, m, f"production frames in state {ST[k]}")
P("")
P("ions and hydration of the site (production, all-atom frames): run | nearest Na+ to Zn: mean, min A | nearest Cl- to Zn: mean, min | nearest Na+ / Cl- to N1: min | water O within 5 A of Zn: mean | within 8 A: mean")
for r in runs:
    f = r.f; pf = r.prodf
    j = dict(Na_mean=float(f["d_Zn_Na_min"][pf].mean()), Na_min=float(f["d_Zn_Na_min"][pf].min()), Cl_mean=float(f["d_Zn_Cl_min"][pf].mean()), Cl_min=float(f["d_Zn_Cl_min"][pf].min()),
             N1_Na_min=float(f["d_N1_Na_min"][pf].min()), N1_Cl_min=float(f["d_N1_Cl_min"][pf].min()), w5=float(f["n_Ow_Zn_5p0"][pf].mean()), w8=float(f["n_Ow_Zn_8p0"][pf].mean()))
    J[r.name]["ions"] = j
    P(f"{r.name} | {j['Na_mean']:.1f}, {j['Na_min']:.1f} | {j['Cl_mean']:.1f}, {j['Cl_min']:.1f} | {j['N1_Na_min']:.1f} / {j['N1_Cl_min']:.1f} | {j['w5']:.2f} | {j['w8']:.2f}")

json.dump(J, open(f"{OUT}/tables.json", "w"), indent=1)
open(f"{OUT}/tables.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
