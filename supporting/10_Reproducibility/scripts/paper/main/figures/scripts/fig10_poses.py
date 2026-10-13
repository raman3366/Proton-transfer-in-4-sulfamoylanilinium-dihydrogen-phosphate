"""Fig. 10 (fig10_poses): production poses of neutral sulfanilamide (DFT input geometry) against the control of each
target, and the contacts of the sulfanilamide pose.

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/fig10_poses.py

One row per target (validated protocol of dock/VALIDATION.md): 3TZF chain B (DHPS, Vina; control: sulfamethoxazole
anion, the co-crystal drug), 7S2K (Sul2, Vina; control: p-aminobenzoate), 6RL9 and 3HS4 (hCA II, AutoDock4Zn; control:
acetazolamide, neutral sulfonamide as AD4Zn requires).  Left: cluster-1 representative pose of each ligand (top pose of
the best-scoring seed of the largest cluster, dock/results/production_contacts.json) in the pocket (residues within
4 A of either ligand, thin sticks).  Right: interaction map of the sulfanilamide pose (_dockmap: polar N/O...N/O
<= 3.2 A with distances, zinc, PLIP pi contacts, other residues within 4 A; every printed distance checked against
facts bulk.dockcont.*).  Titles: mean top-pose scores over seeds (facts dock.prod.<T>.<species>.<scoring>.score).
"""
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _figb as B
import _dock3d as D
import _dockmap as M
from figkit import compose, render3d as r3, style

STEM = 'fig10_poses'
CONTACTS = json.loads((B.DOCK / 'results' / 'production_contacts.json').read_text())
C_SFA = '#F28E2B'
C_CTRL = '#3A8FB7'
ROWS = [
    dict(target='T1_3TZF', pdb='3TZF', name='DHPS, 3TZF chain B', scoring='vina', sname='Vina',
         control='smx_anion', cname='sulfamethoxazole anion'),
    dict(target='T5_7S2K', pdb='7S2K', name='Sul2, 7S2K', scoring='vina', sname='Vina',
         control='paba_anion', cname='$p$-aminobenzoate'),
    dict(target='T3_6RL9', pdb='6RL9', name='hCA II, 6RL9', scoring='ad4zn', sname='AD4Zn',
         control='azm_neutral', cname='acetazolamide'),
    dict(target='T4_3HS4', pdb='3HS4', name='hCA II, 3HS4', scoring='ad4zn', sname='AD4Zn',
         control='azm_neutral', cname='acetazolamide'),
]


def rep_pose(target, species, scoring):
    job = f'{target}/{species}__{scoring}__prod'
    rep = CONTACTS[job]['1']
    lig = D.read_sdf(B.DOCK / 'results' / 'poses' / target / f'{species}__{scoring}__top_pose_per_seed.sdf',
                     prop=('seed', rep['seed']))
    return lig, rep


def check_distances(row, rep):
    """Every distance the map prints must be a facts.json value (bulk.dockcont.<T>.<job>.k / .zn)."""
    job = f"sfa_neutral__{row['scoring']}__prod"
    pref = f"bulk.dockcont.{row['target']}.{job}."
    fk = B.facts()
    vals = {k: v['value'] for k, v in fk.items() if k.startswith(pref)}
    need = [p['distance'] for p in rep['polar_pairs_NO_within_3.5A'] if p['distance'] <= M.POLAR_MAX]
    if rep.get('zinc') and rep['zinc']['distance'] <= 2.6:
        need.append(rep['zinc']['distance'])
    for d in need:
        hit = [k for k, v in vals.items() if abs(v - d) < 1e-6]
        assert hit, (row['target'], d, 'not in facts')
        B.USED[hit[0]] = f'{d:.2f}'
    assert all(p['distance'] <= M.POLAR_MAX or p['distance'] > M.POLAR_MAX for p in rep['polar_pairs_NO_within_3.5A'])


def main():
    style.apply()
    MM = style.MM
    W = style.DOUBLE
    row_h = 47.0 * MM
    head = 4.8 * MM
    legend_h = 6.0 * MM
    H = len(ROWS) * (row_h + head) + legend_h + 1.0 * MM
    fig = plt.figure(figsize=(W, H))
    r_w = 96.0 * MM
    m_x = 0.6 * MM + r_w + 3.0 * MM
    m_w = W - m_x - 0.6 * MM
    renders = []
    for k, row in enumerate(ROWS):
        sfa, rep_s = rep_pose(row['target'], 'sfa_neutral', row['scoring'])
        ctl, rep_c = rep_pose(row['target'], row['control'], row['scoring'])
        rec = D.read_receptor(row['target'])
        keys, _ = D.pocket(rec, [sfa, ctl], cutoff=B.PARAM('dock_pocket_cut'), extra=('ZN', 'MG'))
        sc, atoms, metals = D.build_scene(rec, keys, [ctl, sfa], [C_CTRL, C_SFA], lig_radius=[0.095, 0.17])
        ctx = np.array([a['xyz'] for a in atoms + metals])
        v, info = D.ligand_view([sfa, ctl], ctx)
        fitted = v.fit(sc, pad=0.4)
        B.log(f"{row['target']}: sfa seed {rep_s['seed']} ({rep_s['score']}), control seed {rep_c['seed']} "
              f"({rep_c['score']}); {len(keys)} residues; view {info}")
        row.update(sfa=sfa, ctl=ctl, rep=rep_s, atoms=atoms, metals=metals, scene=sc, view=fitted)
    # one pixel density that fits the tallest / widest pocket into its box
    pxA = min(min(r_w / MM / 25.4 * 600 / row['view'].width, row_h / MM / 25.4 * 600 / row['view'].height)
              for row in ROWS) * 1.0
    for k, row in enumerate(ROWS):
        r = B.render_cached(f"fig09b_{row['target']}", lambda row=row: row['scene'], row['view'], pxA,
                            files=[B.DOCK / row['target'] / 'receptor.pdb'], extra=('v2', row['control']))
        y_top = H - k * (row_h + head)
        s_f = B.fmt(f"dock.prod.{row['pdb']}.sfa_neutral.{row['scoring']}.score", 1)
        s_c = B.fmt(f"dock.prod.{row['pdb']}.{row['control']}.{row['scoring']}.score", 1)
        title = (f"{row['name']}, {row['sname']}: sulfanilamide {s_f}, {row['cname']} {s_c} kcal mol" + r'$^{-1}$')
        compose._header(fig, 0.6 * MM / W, (y_top - 0.2 * MM) / H, f'({"abcd"[k]})', title)
        ih, iw = r.image.shape[:2]
        sc_ = min(r_w / iw, row_h / ih)
        wi, hi = iw * sc_, ih * sc_
        ax = fig.add_axes([(0.6 * MM + (r_w - wi) / 2) / W, (y_top - head - row_h + (row_h - hi) / 2) / H,
                           wi / W, hi / H])
        ax.imshow(r.image, extent=(0, iw, ih, 0), interpolation='antialiased')
        ax.set_xlim(0, iw)
        ax.set_ylim(ih, 0)
        ax.set_axis_off()
        ax._figkit_render = r
        mol, anchors = D.label_molecule(row['atoms'], row['metals'], [row['ctl'], row['sfa']])
        placed, dropped = compose.label_atoms(ax, anchors, mol=mol, surface=None, hydrogens_droppable=False,
                                              placer_kw=dict(n_rings=15, n_dirs=48, leader_from=1), name=f"residue labels {row['target']}")
        B.log(f"{row['target']}: 3D labels {len(placed)}, dropped {[d.text for d in dropped]}")
        axm = fig.add_axes([m_x / W, (y_top - head - row_h) / H, m_w / W, row_h / H])
        check_distances(row, row['rep'])
        sdf = B.DOCK / 'results' / 'poses' / row['target'] / f"sfa_neutral__{row['scoring']}__top_pose_per_seed.sdf"
        M.draw_map(axm, row['target'], 'sfa_neutral', row['rep'], sdf, left_atom=0, right_atom=6)
        B.log(f'map distances without room: {M.MISSING}')
    hs = [Line2D([], [], color=C_SFA, lw=3.0, label='sulfanilamide (DFT geometry)'),
          Line2D([], [], color=C_CTRL, lw=3.0, label='control'),
          Line2D([], [], color=M.GREEN, lw=0.9, ls=(0, (2.0, 1.4)), label='N/O···N/O ≤ 3.2 Å'),
          Line2D([], [], color=M.PI_COL, lw=0.9, ls=(0, (0.8, 1.2)), label='π contact (PLIP)'),
          Line2D([], [], color='0.45', lw=0.9, marker=None, label='hydrophobic (PLIP), arc')]
    fig.legend(handles=hs[:4], loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=style.TICK_PT,
               handlelength=2.0, columnspacing=1.6, frameon=False)
    B.save(fig, STEM, width_mm=W / MM, extra=dict(cameras={r['target']: r['view'].as_dict() for r in ROWS}))


if __name__ == '__main__':
    main()
