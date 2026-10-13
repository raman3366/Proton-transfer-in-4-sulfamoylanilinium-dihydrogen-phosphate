"""Fig. 9 (fig09_redock): redocking / cross-docking overlays behind the docking protocol (gate G2).

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig09_redock.py

Crystal ligand (grey carbons) against the top pose of the seed whose RMSD is the median of its job (coloured carbons),
in the receptor frame (no superposition), with the pocket residues within PARAM dock_pocket_cut (side chains; backbone only where its N or O is within
3.5 A of a ligand), metals and cofactors.
(a) 6RL9 (hCA II): sulfanilamide, DFT input geometry, AutoDock4Zn (production run; facts dock.prod.6RL9.sfa.ad4zn.*);
    the validation run with the ETKDG conformer was bimodal (facts dock.val.T3_6RL9.sfa_neutral__ad4zn.*; caption).
(b) 3HS4 (hCA II) cross-dock: 6RL9 sulfanilamide carried into 3HS4 by C-alpha superposition as the reference;
    validation pose (ETKDG input, AD4Zn; facts dock.val.T4_3HS4.sfa_neutral__ad4zn.*) and production pose (DFT
    input, AD4Zn, ring rotated about the Zn-N/S axis; facts dock.prod.3HS4.sfa.ad4zn.*).
(c) 3TZF chain B (DHPS, Y. pestis): sulfamethoxazole anion, Vina (facts dock.val.T1_3TZF.smx_anion__vina.*;
    passes the 2.0 A criterion only); the production job reproduces the validation run bit for bit
    (dock/RESULTS.md section 7), so its pose file is used.
(d) 7S2K (Sul2): 7,8-dihydropteroate, Vina (facts dock.val.T5_7S2K.dhp_anion__vina.*).
"""
import json

import numpy as np
from matplotlib.lines import Line2D

import _figb as B
import _dock3d as D
from figkit import compose, render3d as r3, style

STEM = 'fig09_redock'
SUMMARY = json.loads((B.DOCK / 'results' / 'production_summary.json').read_text())
VAL = json.loads((B.DOCK / 'analysis' / 'summary.json').read_text())
C_XTAL = '#D2D2D2'
C_POSE = '#F28E2B'
C_POSE2 = '#9C5FB5'


def median_seed_prod(job, refkey):
    per = SUMMARY['jobs'][job]['rmsd_to_crystal'][refkey]['full']['per_seed']
    seeds = [p['seed'] for p in SUMMARY['jobs'][job]['per_seed']]
    med = np.median(per)
    k = int(np.argmin(np.abs(np.array(per) - med)))
    return seeds[k], per[k]


def median_seed_val(run):
    ps = VAL[run]['per_seed']
    r = np.array([p['rmsd_top'] for p in ps])
    k = int(np.argmin(np.abs(r - np.median(r))))
    return ps[k]['seed'], r[k]


def prod_pose(target, tag, seed):
    return D.read_sdf(B.DOCK / 'results' / 'poses' / target / f'{tag}__top_pose_per_seed.sdf', prop=('seed', seed))


def panels():
    out = []
    # (a) 6RL9
    s, r = median_seed_prod('T3_6RL9/sfa_neutral__ad4zn__prod', 'SAN (6RL9, own)')
    ref = D.read_sdf(B.DOCK / 'T3_6RL9' / 'ref_SAN.sdf')
    pose = prod_pose('T3_6RL9', 'sfa_neutral__ad4zn', s)
    B.log(f'6RL9: median seed {s} ({r:.3f} A; own check {D.rmsd_same_order(ref, pose):.3f})')
    t = ('6RL9, sulfanilamide, AD4Zn: ' + B.fmt('dock.prod.6RL9.sfa.ad4zn.rmsd', 2) + ' ± '
         + B.fmt('dock.prod.6RL9.sfa.ad4zn.rmsd_sd', 2) + ' Å, '
         + f"{round(B.F('dock.prod.6RL9.sfa.ad4zn.rmsd_f1') * B.F('dock.prod.6RL9.sfa.ad4zn.rmsd_n'))}/"
         + B.fmt('dock.prod.6RL9.sfa.ad4zn.rmsd_n', 0) + ' ≤ 1 Å')
    out.append(dict(target='T3_6RL9', ref=ref, poses=[(pose, C_POSE)], title=t))
    # (b) 3HS4 cross-dock
    ref = D.read_sdf(B.DOCK / 'results' / 'refs' / 'T4_3HS4__T3_6RL9_ref_SAN_in_T4_3HS4_frame.sdf')
    sv, rv = median_seed_val('T4_3HS4/sfa_neutral__ad4zn')
    pv_ = D.read_pdbqt(B.DOCK / 'T4_3HS4' / 'runs' / 'sfa_neutral__ad4zn' / f'seed{sv:02d}' / 'out.pdbqt')
    sp, rp = median_seed_prod('T4_3HS4/sfa_neutral__ad4zn__prod', 'SAN (6RL9, carried into 3HS4 by Ca superposition)')
    pp = prod_pose('T4_3HS4', 'sfa_neutral__ad4zn', sp)
    B.log(f'3HS4: validation seed {sv} ({rv:.3f}; own {D.rmsd_same_order(ref, pv_):.3f}), production seed {sp} '
          f'({rp:.3f}; own {D.rmsd_same_order(ref, pp):.3f})')
    n_val = B.F('protocol.dock_seeds')
    t = ('3HS4 cross-dock, AD4Zn: ' + B.fmt('dock.val.T4_3HS4.sfa_neutral__ad4zn.rmsd', 2)
         + ' Å (ETKDG), ' + B.fmt('dock.prod.3HS4.sfa.ad4zn.rmsd', 2) + ' Å (DFT)')
    out.append(dict(target='T4_3HS4', ref=ref, poses=[(pv_, C_POSE), (pp, C_POSE2)], title=t))
    # (c) 3TZF
    s, r = median_seed_prod('T1_3TZF/smx_anion__vina__prod', '08D (3TZF, own)')
    ref = D.read_sdf(B.DOCK / 'T1_3TZF' / 'ref_08D.sdf')
    pose = prod_pose('T1_3TZF', 'smx_anion__vina', s)
    B.log(f'3TZF: median seed {s} ({r:.3f}; own {D.rmsd_same_order(ref, pose):.3f})')
    t = ('3TZF, sulfamethoxazole anion, Vina: ' + B.fmt('dock.val.T1_3TZF.smx_anion__vina.rmsd', 2) + ' ± '
         + B.fmt('dock.val.T1_3TZF.smx_anion__vina.rmsd_sd', 2) + ' Å')
    out.append(dict(target='T1_3TZF', ref=ref, poses=[(pose, C_POSE)], title=t))
    # (d) 7S2K
    s, r = median_seed_prod('T5_7S2K/dhp_anion__vina__prod', '78H (7S2K, own)')
    ref = D.read_sdf(B.DOCK / 'T5_7S2K' / 'ref_78H.sdf')
    pose = prod_pose('T5_7S2K', 'dhp_anion__vina', s)
    B.log(f'7S2K: median seed {s} ({r:.3f}; own {D.rmsd_same_order(ref, pose):.3f})')
    t = ('7S2K, dihydropteroate, Vina: ' + B.fmt('dock.val.T5_7S2K.dhp_anion__vina.rmsd', 2) + ' ± '
         + B.fmt('dock.val.T5_7S2K.dhp_anion__vina.rmsd_sd', 2) + ' Å, '
         + f"{round(B.F('dock.val.T5_7S2K.dhp_anion__vina.f1') * n_val)}/{n_val} ≤ 1 Å")
    out.append(dict(target='T5_7S2K', ref=ref, poses=[(pose, C_POSE)], title=t))
    return out


def main():
    style.apply()
    P = panels()
    rr, mols = [], []
    views = []
    for p in P:
        rec = D.read_receptor(p['target'])
        ligs = [p['ref']] + [x for x, _ in p['poses']]
        keys, dmin = D.pocket(rec, ligs, cutoff=B.PARAM('dock_pocket_cut'), extra=('ZN', 'MG'))
        sc, atoms, metals = D.build_scene(rec, keys, ligs, [C_XTAL] + [c for _, c in p['poses']])
        ctx = np.array([a['xyz'] for a in atoms + metals])
        v, info = D.ligand_view(ligs, ctx)
        B.log(f'{p["target"]}: {len(keys)} pocket residues; view {info}')
        fitted = v.fit(sc, pad=0.4)
        p.update(scene=sc, atoms=atoms, metals=metals, ligs=ligs, view=fitted, keys=keys)
        views.append(fitted)
    pxA = r3.px_per_angstrom(views, 84.0)
    for p in P:
        r = B.render_cached(f'fig09a_{p["target"]}', lambda p=p: p['scene'], p['view'], pxA,
                            files=[B.DOCK / p['target'] / 'receptor.pdb'], extra=('v3', p['title']))
        rr.append(r)
    fig, axes = compose.grid(rr, ncols=2, width=style.DOUBLE, titles=[p['title'] for p in P], extra_bottom_mm=6.0)
    for p, ax in zip(P, axes):
        mol, anchors = D.label_molecule(p['atoms'], p['metals'], p['ligs'])
        placed, dropped = compose.label_atoms(ax, anchors, mol=mol, surface=None, hydrogens_droppable=False,
                                              placer_kw=dict(n_rings=11, n_dirs=36, leader_from=1))
        B.log(f'{p["target"]}: labels {len(placed)}, dropped {[d.text for d in dropped]}')
    x0, y0, w, h = fig._figkit_grid['bottom_band']
    W, H = fig.get_size_inches()
    hs = [Line2D([], [], color=C_XTAL, lw=3.0, label='crystal (reference)'),
          Line2D([], [], color=C_POSE, lw=3.0, label='top pose, median-RMSD seed'),
          Line2D([], [], color=C_POSE2, lw=3.0, label='(b) same, DFT input geometry')]
    fig.legend(handles=hs, loc='center', bbox_to_anchor=((x0 + w / 2) / W, (y0 + h / 2) / H), ncol=3,
               fontsize=style.TICK_PT, handlelength=2.0, columnspacing=1.8, frameon=False)
    B.save(fig, STEM, extra=dict(cameras={p['target']: p['view'].as_dict() for p in P}))


if __name__ == '__main__':
    main()
