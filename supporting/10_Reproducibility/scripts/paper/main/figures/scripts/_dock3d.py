"""Docking pictures for figure group B (Figs. 9a, 9b and the SI): receptor pockets, crystal ligands and docked poses
as stick models rendered headless with PyVista (figkit render3d), and 2D interaction maps drawn in matplotlib.

Inputs are read only: dock/T*/receptor.pdb (prepared receptor incl. metals and cofactors), dock/T*/ref_*.sdf and
dock/results/refs/*.sdf (crystal ligands, cross-docking references), dock/results/poses/<T>/*__top_pose_per_seed.sdf
(production top poses, SDF properties seed / score / cluster) and dock/T*/runs/<tag>/seedNN/out.pdbqt (validation
runs, MODEL 1 = top pose).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from matplotlib import colors as mcolors

import _figb as B
from figkit import geometry, render3d as r3, style

DOCK = B.DOCK
AD_ELEMENT = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O', 'O': 'O', 'SA': 'S',
              'S': 'S', 'HD': 'H', 'HS': 'H', 'H': 'H', 'P': 'P', 'F': 'F', 'Cl': 'Cl', 'CL': 'Cl', 'Zn': 'Zn',
              'Mg': 'Mg', 'TZ': None}
METAL_COLOURS = {'ZN': '#7D80B0', 'MG': '#3CA03C'}
RES3 = {'ALA': 'Ala', 'ARG': 'Arg', 'ASN': 'Asn', 'ASP': 'Asp', 'CYS': 'Cys', 'GLN': 'Gln', 'GLU': 'Glu',
        'GLY': 'Gly', 'HIS': 'His', 'HIE': 'His', 'HID': 'His', 'HIP': 'His', 'ILE': 'Ile', 'LEU': 'Leu',
        'LYS': 'Lys', 'MET': 'Met', 'PHE': 'Phe', 'PRO': 'Pro', 'SER': 'Ser', 'THR': 'Thr', 'TRP': 'Trp',
        'TYR': 'Tyr', 'VAL': 'Val'}
COFACTOR_NAMES = {'HH2': 'pterin-PP', 'POP': 'PPi', 'ZN': 'Zn', 'MG': 'Mg'}


@dataclass
class Lig:
    symbols: list
    coords: np.ndarray
    bonds: list
    name: str = ''
    props: dict = field(default_factory=dict)


def res_label(resname, resnum):
    if resname in RES3:
        return f'{RES3[resname]}{resnum}'
    return COFACTOR_NAMES.get(resname, resname)


# ------------------------------------------------------------------------------------------------------------
# reading
# ------------------------------------------------------------------------------------------------------------
def read_receptor(target):
    """Heavy atoms of dock/<target>/receptor.pdb as a list of dicts (waters dropped)."""
    out = []
    for ln in (DOCK / target / 'receptor.pdb').read_text().splitlines():
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        resn = ln[17:20].strip()
        if resn == 'HOH':
            continue
        el = ln[76:78].strip() or re.sub(r'[^A-Za-z]', '', ln[12:16])[:1]
        el = el.capitalize()
        if el == 'H':
            continue
        out.append(dict(name=ln[12:16].strip(), resname=resn, chain=ln[21], resnum=int(ln[22:26]),
                        el=el, xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
                        het=ln.startswith('HETATM')))
    return out


def read_sdf(path, index=None, prop=None, heavy=True):
    """Molecules of an SDF as Lig (heavy atoms); index = record number, or prop=(name, value) to select."""
    from rdkit import Chem
    sup = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    mols = [m for m in sup if m is not None]
    if prop is not None:
        mols = [m for m in mols if m.HasProp(prop[0]) and str(m.GetProp(prop[0])) == str(prop[1])]
        assert len(mols) == 1, (path, prop, len(mols))
        m = mols[0]
    else:
        m = mols[index or 0]
    conf = m.GetConformer()
    keep = [a.GetIdx() for a in m.GetAtoms() if not heavy or a.GetAtomicNum() > 1]
    remap = {k: i for i, k in enumerate(keep)}
    sym = [m.GetAtomWithIdx(k).GetSymbol() for k in keep]
    X = np.array([list(conf.GetAtomPosition(k)) for k in keep])
    bonds = [(remap[b.GetBeginAtomIdx()], remap[b.GetEndAtomIdx()]) for b in m.GetBonds()
             if b.GetBeginAtomIdx() in remap and b.GetEndAtomIdx() in remap]
    props = {p: m.GetProp(p) for p in m.GetPropNames()}
    return Lig(sym, X, bonds, name=str(path.name), props=props)


def read_pdbqt(path, model=1, heavy=True):
    sym, X = [], []
    cur = 0
    for ln in path.read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = int(ln.split()[1])
        elif ln.startswith('ENDMDL') and cur == model:
            break
        elif ln.startswith(('ATOM', 'HETATM')) and (cur == model or cur == 0):
            t = ln[77:79].strip()
            el = AD_ELEMENT.get(t, t[:1])
            if el is None or (heavy and el == 'H'):
                continue
            sym.append(el)
            X.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    X = np.array(X)
    return Lig(sym, X, geometry.detect_bonds(sym, X), name=f'{path.parent.name}/{path.name}#{model}')


def rmsd_same_order(a, b):
    """Plain RMSD for two heavy-atom sets matched by element-graph symmetry (RDKit), for logging only."""
    from itertools import permutations  # noqa: F401 (kept simple: Hungarian on element classes)
    from scipy.optimize import linear_sum_assignment
    tot = 0.0
    n = 0
    for el in set(a.symbols):
        ia = [i for i, s in enumerate(a.symbols) if s == el]
        ib = [i for i, s in enumerate(b.symbols) if s == el]
        if len(ia) != len(ib):
            return float('nan')
        D = np.linalg.norm(a.coords[ia][:, None] - b.coords[ib][None], axis=2) ** 2
        r, c = linear_sum_assignment(D)
        tot += D[r, c].sum()
        n += len(ia)
    return float(np.sqrt(tot / n))


# ------------------------------------------------------------------------------------------------------------
# pocket and scene
# ------------------------------------------------------------------------------------------------------------
def pocket(rec, ligs, cutoff=4.0, extra=()):
    """Residues (chain, resnum, resname) with any heavy atom within `cutoff` of any ligand atom, plus `extra`
    residue names (e.g. metals) wherever they are within 7 A."""
    L = np.vstack([l.coords for l in ligs])
    res = {}
    for a in rec:
        d = np.min(np.linalg.norm(L - a['xyz'], axis=1))
        key = (a['chain'], a['resnum'], a['resname'])
        res[key] = min(res.get(key, 99.0), d)
    keep = [k for k, d in res.items() if d <= cutoff or (k[2] in extra and d <= 7.0)]
    return sorted(keep, key=lambda k: res[k]), res


def _stick_mesh(X, bonds, cols, radius, res=18):
    import pyvista as pv
    cyl = []
    for i, j in bonds:
        a, b = X[i], X[j]
        m = (a + b) / 2
        for p, q, c in ((a, m, cols[i]), (m, b, cols[j])):
            cy = r3._cylinder(p, q, radius, res)
            cy.point_data['rgb'] = np.tile((np.array(mcolors.to_rgb(c)) * 255).astype(np.uint8), (cy.n_points, 1))
            cyl.append(cy)
    joints = r3._sphere_mesh(X, [radius] * len(X), cols, 16)
    parts = cyl + [joints]
    return pv.merge(parts)


def element_colours(symbols, carbon):
    out = []
    for s in symbols:
        if s == 'C':
            out.append(carbon)
        else:
            out.append(style.ELEMENT_COLOURS.get(s, '#FF69B4'))
    return out


def add_sticks(scene, symbols, X, bonds, carbon, radius):
    mesh = _stick_mesh(np.asarray(X), bonds, element_colours(symbols, carbon), radius)
    scene.add(mesh, scalars='rgb', rgb=True, smooth_shading=True, specular=0.3, specular_power=25, diffuse=0.8,
              ambient=0.22)


def pocket_atoms(rec, keys, ligs=None, bb_cut=None):
    """Atoms of the pocket residues (heavy atoms).  With `ligs`, amino-acid backbone N, C, O are kept only for
    residues whose backbone N or O lies within `bb_cut` A of a ligand (e.g. Thr199 N in hCA II); CA is always kept
    so that every side chain is anchored."""
    ks = set(keys)
    bb_cut = B.PARAM('dock_bb_cut') if bb_cut is None else bb_cut
    L = None if not ligs else np.vstack([l.coords for l in ligs])
    keep_bb = set()
    if L is not None:
        for a in rec:
            key = (a['chain'], a['resnum'], a['resname'])
            if key in ks and a['name'] in ('N', 'O') and np.min(np.linalg.norm(L - a['xyz'], axis=1)) <= bb_cut:
                keep_bb.add(key)
    out = []
    for a in rec:
        key = (a['chain'], a['resnum'], a['resname'])
        if key not in ks:
            continue
        if (L is not None and a['resname'] in RES3 and a['resname'] != 'GLY' and a['name'] in ('N', 'C', 'O')
                and key not in keep_bb):
            continue
        out.append(a)
    return out


def residue_bonds(atoms):
    sym = [a['el'] for a in atoms]
    X = np.array([a['xyz'] for a in atoms])
    bonds = geometry.detect_bonds(sym, X)
    same = [(i, j) for i, j in bonds if (atoms[i]['chain'], atoms[i]['resnum']) == (atoms[j]['chain'], atoms[j]['resnum'])
            or {atoms[i]['name'], atoms[j]['name']} == {'C', 'N'}]
    return sym, X, same


def build_scene(rec, keys, ligs, lig_colours, res_carbon='#ABABAB', res_radius=0.065, lig_radius=0.16,
                metal_radius=0.42):
    """Scene: pocket residues (thin sticks: side chains, plus the backbone of residues whose N or O touches a
    ligand), cofactors (thin, tan), metals (spheres), ligands (thick sticks)."""
    sc = r3.Scene(None)
    atoms = [a for a in pocket_atoms(rec, keys, ligs) if a['resname'] not in METAL_COLOURS]
    prot = [a for a in atoms if a['resname'] in RES3]
    cof = [a for a in atoms if a['resname'] not in RES3]
    for group, carbon in ((prot, res_carbon), (cof, '#C9A66B')):
        if group:
            sym, X, bonds = residue_bonds(group)
            add_sticks(sc, sym, X, bonds, carbon, res_radius)
    metals = [a for a in pocket_atoms(rec, keys) if a['resname'] in METAL_COLOURS]
    if metals:
        sc.add(r3._sphere_mesh([a['xyz'] for a in metals], [metal_radius] * len(metals),
                               [METAL_COLOURS[a['resname']] for a in metals], 36),
               scalars='rgb', rgb=True, smooth_shading=True, specular=0.4)
    radii = lig_radius if isinstance(lig_radius, (list, tuple)) else [lig_radius] * len(ligs)
    for lig, col, rad in zip(ligs, lig_colours, radii):
        add_sticks(sc, lig.symbols, lig.coords, lig.bonds, col, rad)
    return sc, atoms, metals


def label_molecule(atoms, metals, ligs):
    """A figkit Molecule over every drawn atom whose labels are the residue names at ONE atom per residue (the
    atom nearest to the ligands; CA-free side chains preferred) and '' elsewhere; for compose.label_atoms.  Amino
    acids farther than PARAM dock_label_cut from every ligand are drawn but not labelled (metals and cofactors are
    always labelled)."""
    L = np.vstack([l.coords for l in ligs])
    allat = atoms + metals
    sym = [a['el'] for a in allat]
    X = np.array([a['xyz'] for a in allat])
    labels = [''] * len(allat)
    best = {}
    cen = L.mean(0)
    for i, a in enumerate(allat):
        key = (a['chain'], a['resnum'], a['resname'])
        # anchor: the atom of the residue farthest from the ligand centre (its outer tip, where there is room for
        # the label), backbone atoms only for glycine
        d = -float(np.linalg.norm(a['xyz'] - cen))
        if a['name'] in ('N', 'C', 'O', 'CA') and a['resname'] in RES3 and a['resname'] != 'GLY':
            d += 50.0
        if key not in best or d < best[key][0]:
            best[key] = (d, i)
    anchors = []
    cut = B.PARAM('dock_label_cut')
    for key, (d, i) in best.items():
        dmin = float(np.min(np.linalg.norm(L - allat[i]['xyz'], axis=1)))
        near = min(float(np.min(np.linalg.norm(L - a['xyz'], axis=1))) for a in allat
                   if (a['chain'], a['resnum'], a['resname']) == key)
        if key[2] in RES3 and near > cut:              # distant residues are drawn but not named
            continue
        labels[i] = res_label(key[2], key[1])
        anchors.append(i)
    raw = [f'{a["resname"]}{a["resnum"]}:{a["name"]}' for a in allat]
    for li, l in enumerate(ligs):                      # ligand atoms: obstacles for the labels, never labelled
        sym = sym + list(l.symbols)
        X = np.vstack([X, l.coords])
        labels += [''] * len(l.symbols)
        raw += [f'lig{li}:{k}' for k in range(len(l.symbols))]
    mol = r3.Molecule(symbols=sym, coords=X, labels=labels, raw_labels=raw, bonds=geometry.detect_bonds(sym, X))
    return mol, anchors


def ligand_view(ligs, context_xyz, prefer_ring=None, n=600):
    """Camera looking at the ligand(s) so that as few ligand atoms as possible are hidden behind pocket atoms
    (a pocket atom nearer to the camera within 0.9 A on screen), with the reference ligand's own plane seen face-on
    where possible (+0.8 |n . d|).  Landscape along the ligand's long axis."""
    L = np.vstack([l.coords for l in ligs])
    C = np.asarray(context_xyz)
    c = L.mean(0)
    ref = ligs[0].coords
    w, V = np.linalg.eigh(np.cov((ref - ref.mean(0)).T))
    nrm = V[:, 0]
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** 0.5) * k
    D = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    best = None
    for d in D:
        v = r3.View(c, d, V[:, 2] if abs(V[:, 2] @ d) < 0.95 else V[:, 1])
        sL = v.screen(L)
        sC = v.screen(C)
        occ = 0
        for p in sL:
            dd = np.hypot(sC[:, 0] - p[0], sC[:, 1] - p[1])
            occ += int(((dd < 0.9) & (sC[:, 2] > p[2] + 0.3)).any())
        score = -occ + 0.8 * abs(nrm @ d)
        if best is None or score > best[0]:
            best = (score, v, occ)
    v = best[1]
    # landscape along the reference ligand's long axis
    long3 = V[:, 2]
    up = np.cross(v.direction, long3)
    v = r3.View(v.focal, v.direction, up)
    return v, dict(occluded_ligand_atoms=best[2], face=round(abs(nrm @ v.direction), 2))
