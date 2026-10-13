"""2D interaction maps of docked poses (LigPlot-like, drawn in matplotlib) for figure group B.

The ligand is drawn from the docked pose itself (RDKit 2D coordinates of the same molecule, so every atom index is
the index used by dock/contacts.py for its role labels, e.g. 'N(sulfonamide)#6').  Contacts come from
dock/results/production_contacts.json (one profile per cluster representative):
  * polar contacts: every receptor N/O ... ligand N/O pair <= POLAR_MAX A (3.2) (geometric, 'polar_pairs_NO_within_3.5A'),
    one dashed green line per pair, labelled with the distance;
  * zinc: the directly measured Zn ... nearest ligand atom distance ('zinc'), solid grey line;
  * pi contacts from PLIP (pi_stacking T/P, pi_cation): dotted line to the ring centre;
  * every other residue within 4 A ('residues_within_4A'): grey label placed towards its nearest ligand atom (from
    the receptor and pose coordinates); PLIP hydrophobic contacts are marked with an arc.
Labels sit on an ellipse round the ligand, in angular order of their anchors, spread so they never overlap.
"""
from __future__ import annotations

import re

import numpy as np
from matplotlib.patches import Arc, Circle, FancyBboxPatch

import _dock3d as D
from figkit import style

import _figb as B  # noqa: E402

POLAR_MAX = B.PARAM('dock_polar_max')
MISSING = []
GREEN = '#2E7D32'
PI_COL = '#8E44AD'
ZN_COL = '#6C6F9E'


def _res_key(s):
    """'THR199' / 'THR199.B' / 'THR199:OG1' / 'HH2278' / 'ZN301' -> ('THR', 199) ..."""
    s = s.split('.')[0].split(':')[0]
    for n in (3, 2):
        if s[:n] in D.RES3 or s[:n] in D.COFACTOR_NAMES:
            return s[:n], int(s[n:])
    m = re.match(r'([A-Z]+)(\d+)', s)
    return m.group(1), int(m.group(2))


def res_text(s):
    rn, num = _res_key(s)
    return D.res_label(rn, num)


def ligand_2d(pose_sdf, seed):
    """RDKit mol (heavy atoms, same order as the pose) with 2D coordinates, bond length 1, oriented so that the
    longest heavy-atom axis is horizontal."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    sup = Chem.SDMolSupplier(str(pose_sdf), removeHs=False)
    mol = next(m for m in sup if m is not None and m.GetProp('seed') == str(seed))
    nH = [sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() == 1) for a in mol.GetAtoms()]
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    assert heavy == list(range(len(heavy))), 'heavy atoms must come first in the pose SDF'
    m = Chem.RemoveHs(mol)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    conf = m.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(m.GetNumAtoms())])
    bl = np.mean([np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in m.GetBonds()])
    xy = (xy - xy.mean(0)) / bl
    w, V = np.linalg.eigh(np.cov(xy.T))
    xy = xy @ V[:, ::-1]
    X3 = mol.GetConformer().GetPositions()[:len(heavy)]
    return m, xy, nH[:len(heavy)], X3


def orient(xy, first_left=None, last_right=None):
    """Flip so that atom `first_left` is left of atom `last_right`."""
    if first_left is not None and last_right is not None and xy[first_left, 0] > xy[last_right, 0]:
        xy = xy * np.array([-1, 1])
    return xy


def draw_ligand(ax, m, xy, nH, lw=1.1, fs=None):
    """Skeletal drawing: aromatic ring as Kekule structure (inner second line), heteroatoms as symbols with H
    counts, bonds clipped at the symbols."""
    from rdkit import Chem
    fs = fs or style.SMALL_PT + 0.5
    mk = Chem.Mol(m)
    Chem.Kekulize(mk, clearAromaticFlags=True)
    ring = [list(r) for r in mk.GetRingInfo().AtomRings()]
    centre = {i: xy[r].mean(0) for r in ring for i in r}
    labelled = {i for i, a in enumerate(mk.GetAtoms()) if a.GetSymbol() != 'C'}
    texts = []
    clip = 0.28
    for b in mk.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        p, q = xy[i].copy(), xy[j].copy()
        u = (q - p) / np.linalg.norm(q - p)
        if i in labelled:
            p = p + u * clip
        if j in labelled:
            q = q - u * clip
        ax.plot([p[0], q[0]], [p[1], q[1]], color='0.1', lw=lw, solid_capstyle='round', zorder=5)
        if b.GetBondTypeAsDouble() == 2:
            if b.IsInRing():
                c = centre[i]
                n = np.array([-u[1], u[0]])
                if (c - (p + q) / 2) @ n < 0:
                    n = -n
                s = 0.16
                a2, b2 = p + n * s + u * 0.14, q + n * s - u * 0.14
            else:
                n = np.array([-u[1], u[0]])
                off = 0.07
                ax.plot([p[0] + n[0] * off * 2, q[0] + n[0] * off * 2], [p[1] + n[1] * off * 2, q[1] + n[1] * off * 2],
                        color='0.1', lw=lw, solid_capstyle='round', zorder=5)
                continue
            ax.plot([a2[0], b2[0]], [a2[1], b2[1]], color='0.1', lw=lw, solid_capstyle='round', zorder=5)
    for i in labelled:
        a = mk.GetAtomWithIdx(i)
        sym = a.GetSymbol()
        h = nH[i]
        txt = sym + ('' if h == 0 else ('H' if h == 1 else f'H$_{h}$'))
        if sym == 'N' and h and xy[i, 0] < xy[:, 0].mean():
            txt = ('H' if h == 1 else f'H$_{h}$') + 'N'
        col = style.ELEMENT_COLOURS.get(sym, 'k')
        col = {'S': '#B8860B', 'O': '#C62828', 'N': '#1F4FB0'}.get(sym, col)
        texts.append(ax.text(xy[i, 0], xy[i, 1], txt, ha='center', va='center', fontsize=fs, color=col, zorder=6,
                             fontweight='bold'))
    return ring, texts


def _box_overlap(a, b, m=0.0):
    return not (a[2] + m <= b[0] or b[2] + m <= a[0] or a[3] + m <= b[1] or b[3] + m <= a[1])


def _seg_box(p, q, box, n=16):
    t = np.linspace(0, 1, n)
    P = p[None] + (q - p)[None] * t[:, None]
    return bool(((P[:, 0] > box[0]) & (P[:, 0] < box[2]) & (P[:, 1] > box[1]) & (P[:, 1] < box[3])).any())


def _seg_seg(p1, p2, q1, q2):
    def o(a, b, c):
        return np.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    return o(p1, p2, q1) != o(p1, p2, q2) and o(q1, q2, p1) != o(q1, q2, p2)


def _seg_disc(p, q, c, r):
    d = q - p
    t = np.clip(((c - p) @ d) / max(d @ d, 1e-12), 0, 1)
    return np.linalg.norm(p + t * d - c) < r


def _edge_point(centre, box, towards):
    """Point where the segment centre -> towards leaves the box (so leaders stop at the box edge)."""
    v = towards - centre
    hw, hh = (box[2] - box[0]) / 2, (box[3] - box[1]) / 2
    s = min(hw / abs(v[0]) if abs(v[0]) > 1e-9 else 1e9, hh / abs(v[1]) if abs(v[1]) > 1e-9 else 1e9)
    return centre + v * min(s, 1.0)


def draw_map(ax, target, tag, rep, pose_sdf, left_atom=None, right_atom=None, window=(13.6, 8.6)):
    """Interaction map of one cluster representative `rep` (a production_contacts.json entry) on axes `ax`.
    Residue boxes are placed greedily (zinc, polar, pi, hydrophobic, other) on candidate points round their anchor
    atom, never overlapping the ligand, another box or another leader; distance labels slide along their line
    until free."""
    m, xy, nH, X3 = ligand_2d(pose_sdf, rep['seed'])
    xy = orient(xy, left_atom, right_atom)
    c = xy.mean(0)
    xy = xy - c
    c = np.zeros(2)
    ax.set_aspect('equal')
    ax.set_xlim(-window[0] / 2, window[0] / 2)
    ax.set_ylim(-window[1] / 2, window[1] / 2)
    ax.axis('off')
    ring, lig_texts = draw_ligand(ax, m, xy, nH)
    rec = D.read_receptor(target)
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    lab_boxes = []
    for t_ in lig_texts:
        bb = t_.get_window_extent(rend)
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        lab_boxes.append((x0 - 0.03, y0 - 0.03, x1 + 0.03, y1 + 0.03))
    items = []
    polar = {}
    for p in rep['polar_pairs_NO_within_3.5A']:
        if p['distance'] > POLAR_MAX:
            continue
        idx = int(p['ligand_atom'].split('#')[1])
        polar.setdefault(p['residue'], []).append((idx, p['receptor_atom'], p['distance']))
    zn = rep.get('zinc')
    if zn and zn['distance'] <= 2.6:
        idx = int(zn['nearest_ligand_atom'].split('#')[1])
        items.append(dict(text='Zn$^{2+}$', kind='zinc', links=[(idx, 'ZN', f"{zn['distance']:.2f}")], anchor=xy[idx]))
    for res, lst in sorted(polar.items(), key=lambda kv: min(d for *_, d in kv[1])):
        # two receptor atoms of one residue at the same ligand atom share one line: '2.57/3.13'
        merged = {}
        for idx, ratom, d in lst:
            merged.setdefault(idx, []).append(d)
        lst = [(idx, None, '/'.join(f'{d:.2f}' for d in sorted(ds))) for idx, ds in merged.items()]
        items.append(dict(text=res_text(res), kind='polar', links=lst, anchor=np.mean([xy[i] for i, *_ in lst], 0)))
    pl = rep['plip']
    rc = xy[ring[0]].mean(0) if ring else c
    done = set(polar)
    for p in pl.get('pi_stacking', []):
        items.append(dict(text=res_text(p['residue']) + f" (π, {p['type']})", kind='pi', anchor=rc, links=[]))
        done.add(p['residue'])
    for p in pl.get('pi_cation', []):
        items.append(dict(text=res_text(p['residue']) + ' (π–cation)', kind='pi', anchor=rc, links=[]))
        done.add(p['residue'])
    hyd = set(pl.get('hydrophobic', []))
    for res, d in rep['residues_within_4A'].items():
        rn, num = _res_key(res)
        if res in done or rn == 'ZN':
            continue
        at = [a for a in rec if a['resname'] == rn and a['resnum'] == num]
        if not at:
            continue
        P = np.array([a['xyz'] for a in at])
        dd = np.linalg.norm(X3[:, None] - P[None], axis=2)
        idx = int(np.unravel_index(np.argmin(dd), dd.shape)[0])
        items.append(dict(text=res_text(res), kind='hyd' if res in hyd else 'near', anchor=xy[idx], links=[]))
    styles = {'polar': (GREEN, '#EAF4EA', GREEN), 'zinc': (ZN_COL, '#ECECF5', ZN_COL),
              'pi': (PI_COL, '#F4ECF8', PI_COL), 'hyd': ('0.3', 'white', '0.45'), 'near': ('0.4', 'white', '0.8')}
    # obstacles: ligand atoms (labelled heteroatoms bigger) and bonds
    discs = [(xy[i], 0.42 if m.GetAtomWithIdx(i).GetSymbol() != 'C' else 0.22) for i in range(len(xy))]
    bonds = [(xy[b.GetBeginAtomIdx()], xy[b.GetEndAtomIdx()]) for b in m.GetBonds()]
    lig_box = (xy[:, 0].min() - 0.45, xy[:, 1].min() - 0.35, xy[:, 0].max() + 0.45, xy[:, 1].max() + 0.35)
    placed_boxes, leaders, dist_boxes = [], [], []
    W2, H2 = window[0] / 2 - 0.05, window[1] / 2 - 0.05
    _sizes = {}

    def text_size(txt):
        if txt in _sizes:
            return _sizes[txt]
        probe = ax.text(0, 0, txt, fontsize=style.SMALL_PT - 0.5, bbox=dict(boxstyle='square,pad=0.04'))
        bb = probe.get_window_extent(rend)
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        probe.remove()
        _sizes[txt] = (x1 - x0, y1 - y0)
        return _sizes[txt]

    def start_of(idx, p0, e):
        u = (e - p0) / np.linalg.norm(e - p0)
        return p0 + u * (0.3 if m.GetAtomWithIdx(idx).GetSymbol() != 'C' else 0.12)

    def find_spot(q0, e, occupied, lines, txt):
        dw, dh = text_size(txt)
        u = (e - q0) / np.linalg.norm(e - q0)
        n = np.array([-u[1], u[0]])
        for f in (0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8):
            for sgn in (1, -1, 2, -2):
                pc = q0 + (e - q0) * f + n * sgn * (dh / 2 + 0.05)
                box = (pc[0] - dw / 2, pc[1] - dh / 2, pc[0] + dw / 2, pc[1] + dh / 2)
                if box[0] < -W2 or box[2] > W2 or box[1] < -H2 or box[3] > H2:
                    continue
                if any(_box_overlap(box, b2, 0.02) for b2 in occupied):
                    continue
                if any(_seg_box(l0, l1, box) for l0, l1 in lines):
                    continue
                if any(_box_overlap(box, (q[0] - r * 0.8, q[1] - r * 0.8, q[0] + r * 0.8, q[1] + r * 0.8))
                       for q, r in discs):
                    continue
                return pc, box
        return None

    for it in items:
        col, fc, ec = styles[it['kind']]
        t = ax.text(0, 0, it['text'], fontsize=style.SMALL_PT, color=col, ha='center', va='center', zorder=8,
                    bbox=dict(boxstyle='round,pad=0.25', fc=fc, ec=ec, lw=0.5))
        bb = t.get_window_extent(rend)                       # text extent + the round box pad (0.25 em)
        padpx = 0.25 * style.SMALL_PT * fig.dpi / 72
        (x0, y0), (x1, y1) = inv.transform([(bb.x0 - padpx, bb.y0 - padpx), (bb.x1 + padpx, bb.y1 + padpx)])
        w, h = x1 - x0, y1 - y0
        a = it['anchor']
        out_dir = a - c if np.linalg.norm(a - c) > 1e-6 else np.array([0.0, 1.0])
        phi0 = np.arctan2(out_dir[1], out_dir[0])
        best = None
        links = it['links'] or [(None, None, None)]
        R0 = {'pi': 2.1, 'polar': 2.0, 'zinc': 2.0}.get(it['kind'], 1.4)
        for R in (R0, R0 + 0.4, R0 + 0.8, R0 + 1.25, R0 + 1.7, R0 + 2.2, R0 + 2.8, R0 + 3.4):
            for dphi in np.radians(np.arange(-180, 180, 5.0)):
                ph = phi0 + dphi
                pc = a + R * np.array([np.cos(ph), np.sin(ph)])
                box = (pc[0] - w / 2, pc[1] - h / 2, pc[0] + w / 2, pc[1] + h / 2)
                if box[0] < -W2 or box[2] > W2 or box[1] < -H2 or box[3] > H2:
                    continue
                if _box_overlap(box, lig_box, 0.0) and any(_box_overlap(box, (q[0] - r, q[1] - r, q[0] + r, q[1] + r))
                                                           for q, r in discs):
                    continue
                if any(_seg_box(p_, q_, box) for p_, q_ in bonds):
                    continue
                if any(_box_overlap(box, b2, 0.08) for b2 in placed_boxes + lab_boxes + dist_boxes):
                    continue
                if any(_seg_box(l0, l1, box) for l0, l1 in leaders):
                    continue
                ok = True
                segs, spots = [], []
                if it['kind'] in ('polar', 'zinc', 'pi'):
                    for idx, _, _ in links:
                        p0 = xy[idx] if idx is not None else a
                        e = _edge_point(pc, box, p0)
                        segs.append((p0, e))
                    for p0, e in segs:
                        if any(_seg_box(p0 + (e - p0) * 0.08, e - (e - p0) * 0.02, b2)
                               for b2 in placed_boxes + dist_boxes):
                            ok = False
                        if any(_seg_seg(p0, e, l0, l1) for l0, l1 in leaders
                               if np.linalg.norm(l0 - p0) > 1e-6):          # leaders from one atom may fan out
                            ok = False
                        # a pi leader starts at the ring centre: only its part outside the ring is checked
                        s0 = p0 + (e - p0) / np.linalg.norm(e - p0) * (1.25 if it['kind'] == 'pi' else 0.0)
                        for k_, (q, r) in enumerate(discs):
                            if np.linalg.norm(q - p0) > 1e-6 and _seg_disc(s0, e, q, r * 0.8):
                                ok = False
                    if ok and it['kind'] in ('polar', 'zinc'):
                        occ = placed_boxes + dist_boxes + lab_boxes + [box]
                        for k_, ((idx, _, dtxt), (p0, e)) in enumerate(zip(links, segs)):
                            q0 = start_of(idx, p0, e)
                            others = leaders + [sg for j, sg in enumerate(segs) if j != k_]
                            sp = find_spot(q0, e, occ, others, dtxt)
                            if sp is None:
                                ok = False
                                break
                            spots.append(sp)
                            occ = occ + [sp[1]]
                if not ok:
                    continue
                cost = (R - R0) + 1.2 * abs(dphi)
                if best is None or cost < best[0]:
                    best = (cost, pc, box, segs, spots)
        assert best is not None, f'no room for {it["text"]} ({target})'
        _, pc, box, segs, spots = best
        t.set_position(pc)
        placed_boxes.append(box)
        it['pos'], it['box'], it['segs'], it['spots'] = pc, box, segs, spots
        for p0, e in segs:
            leaders.append((p0, e))
        for sp in spots:
            dist_boxes.append(sp[1])
    # lines and distance labels
    for it in items:
        col = styles[it['kind']][0]
        if it['kind'] in ('polar', 'zinc'):
            for (idx, ratom, dist), (p0, e), (spc, _) in zip(it['links'], it['segs'], it['spots']):
                q0 = start_of(idx, p0, e)
                ax.plot([q0[0], e[0]], [q0[1], e[1]], color=col, lw=0.8,
                        ls=(0, (2.0, 1.4)) if it['kind'] == 'polar' else '-', zorder=4)
                ax.text(spc[0], spc[1], dist, fontsize=style.SMALL_PT - 0.5, color=col, ha='center',
                        va='center', zorder=7, bbox=dict(boxstyle='square,pad=0.04', fc='white', ec='none'))
        elif it['kind'] == 'pi':
            p0, e = it['segs'][0]
            u = (e - p0) / np.linalg.norm(e - p0)
            ax.plot([p0[0] + u[0] * 0.2, e[0]], [p0[1] + u[1] * 0.2, e[1]], color=col, lw=0.8, ls=(0, (0.8, 1.2)),
                    zorder=4)
            ax.plot(*p0, marker='o', ms=2.0, color=col, zorder=4)
        elif it['kind'] == 'hyd':
            v = it['anchor'] - it['pos']
            v = v / np.linalg.norm(v)
            e = _edge_point(it['pos'], it['box'], it['anchor'])
            ang = np.degrees(np.arctan2(v[1], v[0]))
            ax.add_patch(Arc(e + v * 0.12, 0.34, 0.34, angle=ang, theta1=-65, theta2=65, color='0.35', lw=0.8,
                             zorder=7))
    # shrink the window to the content (the axes keep their aspect and sit at the right edge of their box)
    xs = [b_[0] for b_ in placed_boxes + dist_boxes + lab_boxes] + [b_[2] for b_ in placed_boxes + dist_boxes + lab_boxes]
    ys = [b_[1] for b_ in placed_boxes + dist_boxes + lab_boxes] + [b_[3] for b_ in placed_boxes + dist_boxes + lab_boxes]
    ax.set_xlim(min(xs + list(xy[:, 0])) - 0.1, max(xs + list(xy[:, 0])) + 0.1)
    ax.set_ylim(min(ys + list(xy[:, 1])) - 0.1, max(ys + list(xy[:, 1])) + 0.1)
    ax.set_anchor('E')
    return items
