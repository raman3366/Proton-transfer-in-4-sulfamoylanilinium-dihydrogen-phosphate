"""Add named crystal waters to an already prepared receptor, changing nothing else, and compare the
docking runs of the water-containing receptor with the runs of the dry one.

    python dock/water_variant.py --target T3_6RL9 --waters 423 516 553 568
    # writes dock/T3_6RL9_wat4/{receptor.pdbqt,receptor.pdb,box.json,ref_*.sdf,variant.json}
    python dock/water_variant.py --target T3_6RL9 --waters 423 516 553 568 --rule bare
    # the same waters by another hydrogen rule: dock/T3_6RL9_wat4bare/
    python dock/water_variant.py --target T4_3HS4 --map-sites-from T3_6RL9 --waters 423 516 553 568
    # where those water sites lie in another structure of the protein (prints, writes nothing)
    python dock/water_variant.py --rescore --scratch <dir>
    # stored modes of the dry runs scored in the maps of every water variant (no search):
    # dock/<variant>/rescore_fixed_poses.json
    python dock/water_variant.py --target T3_6RL9_wat4 --compare --scratch <dir>   (or --compare-all)
    # every run of the variant against the run of the same name in its dry source:
    # dock/T3_6RL9_wat4/comparison_water.json
    python dock/water_variant.py --tables T3_6RL9_wat4 T3_6RL9_glu106GLU_wat4 ...
    # the tables of dock/WATER4.md

Why this exists (6 October 2026).  The reported docking used receptors without water (the protocol
keeps water only on Mg2+).  The ring-rotated pose of sulfanilamide in 6RL9, the orientation that
AutoDock4Zn cannot separate from the crystal one, lies on four water sites of the crystal structure
(waters 423, 516, 553 and 568 of chain A, full occupancy, B 11.5-17.6 A2) that are present in 6RL9
together with the ligand in its crystal pose (md/SYSTEM.md section 4).  The dry receptor cannot see
them.  This script puts them back into the very receptor that was docked.

What a variant is.  The prepared receptor, line for line, plus the named waters of the deposited
structure appended after the last atom:

    oxygen      at the deposited coordinates (the prepared receptor is in the deposited frame; the
                script measures the zinc and the CA atoms against the deposited file and records it),
                AutoDock type OA, charge -0.411: the way prep_receptor.py and control_gyrb.py write a
                kept water (the Gasteiger charge of the oxygen in Meeko's HOH template);
    hydrogens   by a stated rule (below), O-H 0.9572 A, H-O-H 104.52 deg, type HD, charge +0.206
                (the hydrogens of Meeko's HOH template, 0.2058): a water is neutral to 0.001 e.

The hydrogens come right after their oxygen, which AutoGrid needs to find them as its bonded atoms.

Hydrogen rules (--rule):

    hbond   the orientation of the set of waters that makes the most hydrogen bonds with the atoms
            of this receptor and among the waters themselves.  A hydrogen bond: H...A <= 2.5 A and
            D-H...A >= 130 deg; acceptors are the OA/NA atoms of the receptor and the water oxygens,
            donors the polar hydrogens (HD) of the receptor, whose positions are those of the
            prepared receptor.  Candidate orientations of each water: both hydrogens towards a pair
            of acceptors within 3.5 A (bisectors aligned), one hydrogen exactly at an acceptor with
            the other turned in 30 deg steps, and 72 orientations that point at nothing.  Every
            combination of candidates is scored: +1 per hydrogen bond, 0.1 x cos(deviation from
            linear) per hydrogen bond to break ties, and -10 per pair of hydrogens closer than
            1.8 A, which forbids two waters donating to each other along one O...O line and a
            water hydrogen pointing at a donor hydrogen of the protein.
            The crystal ligand and the crystal waters that are not kept take no part: the receptor
            must not be told where the ligand is.
    rev     the best network of the same kind in which no water-water hydrogen bond of the `hbond`
            network is kept.  When the waters form a ring whose direction the receptor does not
            decide, this is the ring turned the other way.
    bare    no hydrogens: the oxygen alone, as in the precedent (dock/control_gyrb.py water_variant
            and the Mg-bound waters of prep_receptor.py).  With AutoDock4 the water then carries its
            oxygen charge of -0.411 e uncompensated.
    fixed   every water in one standard orientation that looks at nothing in particular (H1 along
            +x, H2 in the xy plane), which is what tleap gives a crystal water before minimisation.

--compare analyses the dry and the water run with the validation code itself (through
glu_protonation_variant.analyse_run: RMSD, site checks, PoseBusters, gates) and then classifies every
mode of every seed: near-native (within 1.0 A of the crystal ligand), ring-rotated in the paper's
sense (beyond 2.0 A), and among those the rotated pose proper (within 2.0 A of a reference rotated
pose, ROT_REF below) or something else.  It also pools the modes beyond 1.0 A of all seeds into 2.0 A
clusters, to show what the search finds besides the crystal pose.  Nothing is written outside the
variant directory and the scratch directory.

A receptor with waters fixed in place forbids their displacement by construction.  What such a run can
show is whether a pose needs the space of the waters; it cannot show what displacing them would cost.
"""
import argparse
import hashlib
import itertools
import json
import pathlib
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

DOCK = pathlib.Path(__file__).resolve().parent
R_OH, A_HOH = 0.9572, 104.52
Q_O, Q_H = -0.411, 0.206
HB_DIST, HB_ANGLE, HH_CLASH, NEIGHBOUR, CLASH_COST = 2.5, 130.0, 1.8, 3.5, 10.0


# --------------------------------------------------------------------------------------------------
# small geometry
# --------------------------------------------------------------------------------------------------
def unit(v):
    return v / np.linalg.norm(v)


def angle(a, b, c):
    v, w = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(v @ w / np.linalg.norm(v) / np.linalg.norm(w), -1, 1))))


def dist(a, b):
    return float(np.linalg.norm(a - b))


def coords(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def is_atom(line):
    return line.startswith(('ATOM', 'HETATM'))


def rot(axis, deg, v):
    k, t = unit(axis), np.radians(deg)
    return v * np.cos(t) + np.cross(k, v) * np.sin(t) + k * (k @ v) * (1 - np.cos(t))


def water_from_bisector(o, bis, normal):
    """Hydrogens of a water at o whose H-O-H bisector is `bis` and whose plane has normal `normal`."""
    b = unit(bis)
    n = unit(normal - (normal @ b) * b)
    inplane = np.cross(n, b)
    h = np.radians(A_HOH / 2)
    return (o + R_OH * (np.cos(h) * b + np.sin(h) * inplane),
            o + R_OH * (np.cos(h) * b - np.sin(h) * inplane))


def water_from_h1(o, d1, turn):
    """H1 along d1; H2 on the cone of A_HOH about it, turned by `turn` degrees."""
    d1 = unit(d1)
    ref = np.array([1.0, 0, 0]) if abs(d1[0]) < 0.9 else np.array([0, 1.0, 0])
    perp = unit(np.cross(d1, ref))
    d2 = rot(d1, turn, rot(perp, A_HOH, d1))
    return o + R_OH * d1, o + R_OH * unit(d2)


# --------------------------------------------------------------------------------------------------
# the deposited structure
# --------------------------------------------------------------------------------------------------
def base_and_spec(target):
    from analyse import base_target
    from targets import TARGETS
    base = base_target(target)
    return base, TARGETS[base]


def deposited(pdb, chain):
    """Waters, zinc and CA atoms of one author chain of the deposited file."""
    from targets import read_structure
    st = read_structure(pdb)
    waters, ca, zn = {}, {}, []
    for ch in st[0]:
        if ch.name != chain:
            continue
        for res in ch:
            if res.name == 'HOH':
                waters.setdefault(res.seqid.num, []).extend(
                    dict(altloc=a.altloc if a.altloc != '\0' else '', occ=a.occ, b=a.b_iso,
                         x=np.array([a.pos.x, a.pos.y, a.pos.z])) for a in res)
            elif res.name == 'ZN':
                zn += [np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res]
            else:
                for a in res:
                    if a.name == 'CA':
                        ca.setdefault(res.seqid.num, []).append(np.array([a.pos.x, a.pos.y, a.pos.z]))
    return waters, ca, zn


def frame_check(src, pdb, chain):
    """Is the prepared receptor in the frame of the deposited file?  Zinc and CA atoms, no fit."""
    _, ca, zn = deposited(pdb, chain)
    rows = [l for l in (src / 'receptor.pdb').read_text().splitlines() if is_atom(l)]
    d = []
    for l in rows:
        if l[12:16].strip() == 'CA' and l[17:20].strip() != 'CA':
            n = int(l[22:26])
            if n in ca:
                d.append(min(dist(coords(l), x) for x in ca[n]))
    d = np.array(d)
    z = [coords(l) for l in rows if l[17:20].strip() == 'ZN']
    out = dict(deposited_file=f'pdb/{pdb}.cif.gz', chain=chain, n_ca_compared=int(len(d)),
               ca_rms_A=round(float(np.sqrt((d ** 2).mean())), 3), ca_max_A=round(float(d.max()), 3),
               zinc_shift_A=round(min(dist(a, b) for a in z for b in zn), 3) if z and zn else None,
               note='no superposition; the CA shifts are those of the restrained minimisation of the '
                    'receptor preparation (an alternate conformer is compared with its nearest copy)')
    out['same_frame'] = bool(out['ca_max_A'] < 0.5 and (out['zinc_shift_A'] is None or out['zinc_shift_A'] < 0.01))
    return out


# --------------------------------------------------------------------------------------------------
# hydrogens
# --------------------------------------------------------------------------------------------------
def receptor_atoms(src):
    at = []
    for l in (src / 'receptor.pdbqt').read_text().splitlines():
        if is_atom(l):
            at.append(dict(label=f"{l[17:20].strip()}{int(l[22:26])} {l[12:16].strip()}", type=l[77:79].strip(),
                           x=coords(l), res=(l[21], int(l[22:26]))))
    return at


def candidates(o, acc):
    """Candidate hydrogen positions of one water (see the module docstring)."""
    out = []
    for a, b in itertools.combinations(acc, 2):
        da, db = unit(a - o), unit(b - o)
        if 60.0 <= angle(a, o, b) <= 160.0:
            out.append(water_from_bisector(o, da + db, np.cross(da, db)))
    for a in acc:
        out += [water_from_h1(o, a - o, t) for t in range(0, 360, 30)]
    # orientations that look at nothing: H1 along the 12 vertices of an icosahedron, six turns each
    g = (1 + 5 ** 0.5) / 2
    ico = [np.array(v, float) for v in ((0, 1, g), (0, -1, g), (0, 1, -g), (0, -1, -g), (1, g, 0), (-1, g, 0),
                                        (1, -g, 0), (-1, -g, 0), (g, 0, 1), (g, 0, -1), (-g, 0, 1), (-g, 0, -1))]
    for v in ico:
        out += [water_from_h1(o, v, t) for t in range(0, 360, 60)]
    return out


def hb(h, d, a):
    """(is a hydrogen bond, quality) for donor d with hydrogen h and acceptor a."""
    r, ang = dist(h, a), angle(d, h, a)
    return (r <= HB_DIST and ang >= HB_ANGLE), float(np.cos(np.radians(180.0 - ang)))


def water_edges(oxy, hyd):
    """Directed water-water hydrogen bonds (donor number, acceptor number) of a placement."""
    return sorted((n, m) for n in oxy for m in oxy if n != m
                  and any(hb(h, oxy[n], oxy[m])[0] for h in hyd.get(n, ())))


def place_hbond(oxy, rec, forbid=()):
    """The orientation of all waters that maximises the hydrogen-bond count (exhaustive over the
    candidates).  oxy: {number: O position}; rec: receptor atoms; forbid: directed water-water
    hydrogen bonds (donor, acceptor) that the network must not contain.  Returns {number: (H1, H2)}
    and a record of the search."""
    nums = sorted(oxy)
    acc_rec = [a for a in rec if a['type'] in ('OA', 'NA')]
    hd = [a for a in rec if a['type'] == 'HD']
    cand, single = {}, {}
    for n in nums:
        o = oxy[n]
        acc = [a['x'] for a in acc_rec if dist(a['x'], o) <= NEIGHBOUR] + \
              [oxy[m] for m in nums if m != n and dist(oxy[m], o) <= NEIGHBOUR]
        cand[n] = candidates(o, acc)
        near_acc = [a for a in acc_rec if dist(a['x'], o) <= NEIGHBOUR + 1.0]
        near_h = [a for a in hd if dist(a['x'], o) <= 4.0]
        s = np.zeros(len(cand[n]))
        for i, hs in enumerate(cand[n]):
            for h in hs:
                for a in near_acc:
                    ok, q = hb(h, o, a['x'])
                    if ok:
                        s[i] += 1 + 0.1 * q
                for p in near_h:
                    if dist(h, p['x']) < HH_CLASH:
                        s[i] -= CLASH_COST
        single[n] = s
    pair = {}
    for n, m in itertools.combinations(nums, 2):
        p = np.zeros((len(cand[n]), len(cand[m])))
        if dist(oxy[n], oxy[m]) <= NEIGHBOUR + 2.0:
            for i, hn in enumerate(cand[n]):
                for j, hm in enumerate(cand[m]):
                    v = 0.0
                    for h in hn:
                        ok, q = hb(h, oxy[n], oxy[m])
                        if ok:
                            v += 1 + 0.1 * q
                    for h in hm:
                        ok, q = hb(h, oxy[m], oxy[n])
                        if ok:
                            v += 1 + 0.1 * q
                    v -= CLASH_COST * sum(1 for a in hn for b in hm if dist(a, b) < HH_CLASH)
                    if ((n, m) in forbid and any(hb(h, oxy[n], oxy[m])[0] for h in hn)) or \
                            ((m, n) in forbid and any(hb(h, oxy[m], oxy[n])[0] for h in hm)):
                        v -= 1000.0
                    p[i, j] = v
        pair[(n, m)] = p
    # exhaustive maximum over all combinations (the last two waters vectorised)
    best, arg, n_comb = -1e9, None, 1
    for n in nums:
        n_comb *= len(cand[n])
    head, tail = nums[:-2], nums[-2:]
    for combo in itertools.product(*[range(len(cand[n])) for n in head]):
        s0 = sum(single[n][i] for n, i in zip(head, combo))
        s0 += sum(pair[(a, b)][i, j] for (a, i), (b, j) in itertools.combinations(list(zip(head, combo)), 2))
        if len(tail) == 2:
            t1, t2 = tail
            m = s0 + single[t1][:, None] + single[t2][None, :] + pair[(t1, t2)]
            for n, i in zip(head, combo):
                m = m + pair[(n, t1)][i][:, None] + pair[(n, t2)][i][None, :]
            k = np.unravel_index(int(np.argmax(m)), m.shape)
            if m[k] > best:
                best, arg = float(m[k]), dict(zip(head + tail, list(combo) + [int(k[0]), int(k[1])]))
        else:                                             # one or two waters in all
            t1 = tail[-1]
            m = s0 + single[t1]
            for n, i in zip(head, combo):
                m = m + pair[(n, t1)][i]
            k = int(np.argmax(m))
            if m[k] > best:
                best, arg = float(m[k]), dict(zip(head + [t1], list(combo) + [k]))
    hyd = {n: cand[n][arg[n]] for n in nums}
    return hyd, dict(candidates_per_water={str(n): len(cand[n]) for n in nums}, combinations=int(n_comb),
                     best_score=round(best, 3), water_water_hydrogen_bonds=[f'HOH{a} -> HOH{b}' for a, b in water_edges(oxy, hyd)],
                     forbidden=[f'HOH{a} -> HOH{b}' for a, b in forbid])


def place_fixed(oxy):
    h = np.radians(A_HOH)
    return {n: (o + R_OH * np.array([1.0, 0, 0]), o + R_OH * np.array([np.cos(h), np.sin(h), 0]))
            for n, o in oxy.items()}, dict(orientation='H1 along +x, H2 in the xy plane, the same for every water')


def partners(oxy, hyd, rec, extra=None):
    """Hydrogen bonds of every water as placed: what each hydrogen points at, and what points at the
    oxygen.  `extra`: crystal neighbours that are not in the receptor (label -> position), listed by
    distance only."""
    out = {}
    acc_rec = [a for a in rec if a['type'] in ('OA', 'NA')]
    heavy = [a for a in rec if a['type'] != 'HD']
    for n, o in sorted(oxy.items()):
        row = dict(donates=[], accepts=[], hydrogens_without_partner=0, close_H_H=[])
        for k, h in enumerate(hyd.get(n, ())):
            found = False
            for a in acc_rec:
                ok, _ = hb(h, o, a['x'])
                if ok:
                    row['donates'].append(dict(hydrogen=f'H{k + 1}', to=a['label'], H_A=round(dist(h, a['x']), 2),
                                               O_A=round(dist(o, a['x']), 2), angle=round(angle(o, h, a['x']), 1)))
                    found = True
            for m, om in oxy.items():
                if m != n:
                    ok, _ = hb(h, o, om)
                    if ok:
                        row['donates'].append(dict(hydrogen=f'H{k + 1}', to=f'HOH{m} O', H_A=round(dist(h, om), 2),
                                                   O_A=round(dist(o, om), 2), angle=round(angle(o, h, om), 1)))
                        found = True
            row['hydrogens_without_partner'] += 0 if found else 1
            for p in rec:
                if p['type'] == 'HD' and dist(h, p['x']) < 2.0:
                    row['close_H_H'].append(dict(hydrogen=f'H{k + 1}', to=p['label'], distance=round(dist(h, p['x']), 2)))
        for p in rec:
            if p['type'] == 'HD' and dist(p['x'], o) <= HB_DIST:
                d = min((a for a in heavy if a['res'] == p['res']), key=lambda a: dist(a['x'], p['x']))
                if angle(d['x'], p['x'], o) >= HB_ANGLE:
                    row['accepts'].append(dict(frm=f"{d['label']}-H", H_O=round(dist(p['x'], o), 2),
                                               D_O=round(dist(d['x'], o), 2), angle=round(angle(d['x'], p['x'], o), 1)))
        for m, om in oxy.items():
            if m != n:
                for k, h in enumerate(hyd.get(m, ())):
                    ok, _ = hb(h, om, o)
                    if ok:
                        row['accepts'].append(dict(frm=f'HOH{m} H{k + 1}', H_O=round(dist(h, o), 2),
                                                   D_O=round(dist(om, o), 2), angle=round(angle(om, h, o), 1)))
        row['neighbours_within_3p5A'] = sorted(
            [dict(atom=a['label'], distance=round(dist(a['x'], o), 2)) for a in heavy if dist(a['x'], o) <= NEIGHBOUR]
            + [dict(atom=f'HOH{m} O', distance=round(dist(om, o), 2)) for m, om in oxy.items()
               if m != n and dist(om, o) <= NEIGHBOUR]
            + [dict(atom=f'{lab} (crystal, not in the receptor)', distance=round(dist(x, o), 2))
               for lab, x in (extra or {}).items() if dist(x, o) <= NEIGHBOUR], key=lambda r: r['distance'])
        out[str(n)] = row
    return out


# --------------------------------------------------------------------------------------------------
# the variant
# --------------------------------------------------------------------------------------------------
def make_variant(target, numbers, rule, name=None):
    base, spec = base_and_spec(target)
    src = DOCK / target
    name = name or f"{target}_wat{len(numbers)}{'' if rule == 'hbond' else rule}"
    dst = DOCK / name
    waters, _, _ = deposited(spec['pdb'], spec['chain'])
    frame = frame_check(src, spec['pdb'], spec['chain'])
    if not frame['same_frame']:
        raise SystemExit(f'{target} is not in the frame of the deposited file: {frame}')
    oxy, meta = {}, {}
    for n in numbers:
        if n not in waters:
            raise SystemExit(f'water {n} of chain {spec["chain"]} is not in {spec["pdb"]}')
        if len(waters[n]) != 1:
            raise SystemExit(f'water {n} has alternate positions; name one explicitly (not supported here)')
        w = waters[n][0]
        oxy[n], meta[n] = w['x'], dict(occupancy=w['occ'], b_iso=round(w['b'], 1))
    rec = receptor_atoms(src)
    if rule == 'hbond':
        hyd, search = place_hbond(oxy, rec)
    elif rule == 'rev':
        first, s1 = place_hbond(oxy, rec)
        hyd, search = place_hbond(oxy, rec, forbid=tuple(water_edges(oxy, first)))
        search['first_network'] = s1
    elif rule == 'fixed':
        hyd, search = place_fixed(oxy)
    else:
        hyd, search = {}, dict(orientation='no hydrogens (bare oxygen, the precedent)')
    # crystal neighbours that the receptor does not hold, for the record: other waters and the ligand
    extra = {f'HOH{m} O': w[0]['x'] for m, w in waters.items() if m not in oxy and len(w) == 1}
    try:
        from rdkit import Chem
        for sdf in sorted(p for p in src.glob('ref_*.sdf') if not p.name.startswith('ref_cross_')):
            mol = Chem.MolFromMolFile(str(sdf))
            c = mol.GetConformer()
            for a in mol.GetAtoms():
                extra[f'{sdf.stem[4:]} {a.GetSymbol()}{a.GetIdx() + 1} (crystal ligand)'] = np.array(
                    list(c.GetAtomPosition(a.GetIdx())))
    except Exception:
        pass
    dst.mkdir(parents=True, exist_ok=True)
    pq = (src / 'receptor.pdbqt').read_text().splitlines()
    pdb = (src / 'receptor.pdb').read_text().splitlines()
    serial_q = max(int(l[6:11]) for l in pq if is_atom(l))
    serial_p = max(int(l[6:11]) for l in pdb if is_atom(l))
    add_q, add_p = [], []
    for n in sorted(oxy):
        o = oxy[n]
        serial_q += 1
        serial_p += 1
        add_q.append(f'ATOM  {serial_q:5d}  O   HOH {spec["chain"]}{n:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00  0.00    {Q_O:+6.3f} OA')
        add_p.append(f'HETATM{serial_p:5d}  O   HOH {spec["chain"]}{n:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00{meta[n]["b_iso"]:6.2f}           O  ')
        for k, h in enumerate(hyd.get(n, ())):
            serial_q += 1
            serial_p += 1
            add_q.append(f'ATOM  {serial_q:5d}  H{k + 1}  HOH {spec["chain"]}{n:4d}    {h[0]:8.3f}{h[1]:8.3f}{h[2]:8.3f}  1.00  0.00    {Q_H:+6.3f} HD')
            add_p.append(f'HETATM{serial_p:5d}  H{k + 1}  HOH {spec["chain"]}{n:4d}    {h[0]:8.3f}{h[1]:8.3f}{h[2]:8.3f}  1.00{meta[n]["b_iso"]:6.2f}           H  ')
    new_q = pq + add_q
    end = [i for i, l in enumerate(pdb) if l.startswith('END')]
    cut = end[0] if end else len(pdb)
    new_p = pdb[:cut] + add_p + pdb[cut:]
    (dst / 'receptor.pdbqt').write_text('\n'.join(new_q) + '\n')
    (dst / 'receptor.pdb').write_text('\n'.join(new_p) + '\n')
    assert new_q[:len(pq)] == pq and [l for l in new_p if l not in set(add_p)] == pdb
    own_refs = [p.name for p in src.glob('ref_*.sdf') if not p.name.startswith('ref_cross_')]
    for extra_file in ['box.json', 'vina_config.txt', *own_refs]:
        if (src / extra_file).exists():
            shutil.copy2(src / extra_file, dst / extra_file)
    box = json.loads((src / 'box.json').read_text())
    centre, half = np.array(box['center']), np.array(box['size']) / 2
    lig = np.array([x for lab, x in extra.items() if 'crystal ligand' in lab])
    record = dict(
        source=target, variant=name, rule=rule,
        note='the prepared receptor, line for line, plus the named crystal waters appended; nothing else changed',
        deposited=dict(pdb=spec['pdb'], chain=spec['chain']), frame_check=frame,
        water_atoms=dict(oxygen=dict(type='OA', charge=Q_O, source='prep_receptor.py water_pdbqt_atoms / '
                                     'control_gyrb.py water_variant (Gasteiger charge of O in Meeko\'s HOH template)'),
                         hydrogen=dict(type='HD', charge=Q_H, O_H_A=R_OH, H_O_H_deg=A_HOH,
                                       source='Gasteiger charge of H in Meeko\'s HOH template (0.2058)')
                         if hyd else None,
                         charge_per_water=round(Q_O + (2 * Q_H if hyd else 0), 3)),
        hydrogen_search=search,
        waters={str(n): dict(**meta[n], O=[round(float(v), 3) for v in oxy[n]],
                             H=[[round(float(v), 3) for v in h] for h in hyd.get(n, ())],
                             inside_box=bool((np.abs(oxy[n] - centre) <= half).all()),
                             nearest_crystal_ligand_atom_A=round(float(np.linalg.norm(lig - oxy[n], axis=1).min()), 2)
                             if len(lig) else None)
                for n in sorted(oxy)},
        hydrogen_bonds=partners(oxy, hyd, rec, extra),
        file_check=dict(receptor_pdbqt_lines=(len(pq), len(new_q)), receptor_pdb_lines=(len(pdb), len(new_p)),
                        source_lines_byte_identical_and_in_order=True),
        receptor_total_charge=(round(sum(float(l[70:76]) for l in pq if is_atom(l)), 3),
                               round(sum(float(l[70:76]) for l in new_q if is_atom(l)), 3)))
    for f in ('receptor.pdbqt', 'receptor.pdb'):
        record[f'{f}_sha256'] = dict(source=hashlib.sha256((src / f).read_bytes()).hexdigest(),
                                     variant=hashlib.sha256((dst / f).read_bytes()).hexdigest())
    (dst / 'variant.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f"{name}: {len(oxy)} waters of {spec['pdb']} chain {spec['chain']} added to {target} (rule {rule}); frame: CA rms "
          f"{frame['ca_rms_A']} max {frame['ca_max_A']} A, zinc {frame['zinc_shift_A']} A; {search}")
    for n, r in record['hydrogen_bonds'].items():
        print(f"  HOH{n}: donates to " + (', '.join(f"{d['to']} ({d['H_A']} A, {d['angle']} deg)" for d in r['donates']) or 'nothing')
              + '; accepts from ' + (', '.join(f"{d['frm']} ({d['H_O']} A, {d['angle']} deg)" for d in r['accepts']) or 'nothing')
              + (f"; close H...H {r['close_H_H']}" if r['close_H_H'] else ''))
    return 0


def map_sites(target, source, numbers, out_json=None):
    """Where the named water sites of `source` lie in `target` (another structure of the same
    protein): CA superposition of the two prepared receptors, then the nearest deposited water."""
    import analyse as A
    sb, sspec = base_and_spec(source)
    tb, tspec = base_and_spec(target)
    ws, _, _ = deposited(sspec['pdb'], sspec['chain'])
    wt, _, _ = deposited(tspec['pdb'], tspec['chain'])
    src, dst = A.ca_atoms(DOCK / source / 'receptor.pdb'), A.ca_atoms(DOCK / target / 'receptor.pdb')
    keys = sorted(set(src) & set(dst))
    P, Q = np.array([src[k] for k in keys]), np.array([dst[k] for k in keys])
    R, t = A.kabsch(P, Q)
    rms = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    from rdkit import Chem
    lig = []
    for sdf in sorted(q for q in (DOCK / target).glob('ref_*.sdf') if not q.name.startswith('ref_cross_')):
        mol = Chem.MolFromMolFile(str(sdf))
        c = mol.GetConformer()
        lig += [(f'{sdf.stem[4:]} {a.GetSymbol()}{a.GetIdx() + 1}', np.array(list(c.GetAtomPosition(a.GetIdx()))))
                for a in mol.GetAtoms()]
    from targets import read_structure
    dep = [(f'{res.name}{res.seqid.num} {a.name}', a.occ, round(a.b_iso, 1), np.array([a.pos.x, a.pos.y, a.pos.z]))
           for ch in read_structure(tspec['pdb'])[0] if ch.name == tspec['chain'] for res in ch for a in res
           if res.name != 'HOH']
    rows = []
    for n in numbers:
        x = R @ ws[n][0]['x'] + t
        near = sorted(((dist(x, w['x']), m, w) for m, ww in wt.items() for w in ww), key=lambda r: r[0])[:2]
        dl, lab = min((dist(x, y), lab) for lab, y in lig) if lig else (None, None)
        nd = sorted(((dist(x, y), lab_, occ, b) for lab_, occ, b, y in dep), key=lambda r: r[0])[:2]
        rows.append(dict(site=f"{sspec['pdb']} HOH{n}", position_in_target_frame=[round(float(v), 3) for v in x],
                         nearest_waters=[dict(water=m, altloc=w['altloc'], occupancy=w['occ'], b_iso=round(w['b'], 1),
                                              distance_A=round(d, 2)) for d, m, w in near],
                         nearest_crystal_ligand_atom=dict(atom=lab, distance_A=round(dl, 2)) if lig else None,
                         nearest_deposited_atoms_not_water=[dict(atom=a_, occupancy=o_, b_iso=b_, distance_A=round(d_, 2))
                                                            for d_, a_, o_, b_ in nd],
                         a_water_of_the_target_within_1A=bool(near and near[0][0] <= 1.0)))
    out = dict(source=source, target=target, n_ca=len(keys), ca_rmsd_A=round(rms, 3), sites=rows,
               n_sites_present_as_water=sum(r['a_water_of_the_target_within_1A'] for r in rows))
    if out_json:
        pathlib.Path(out_json).write_text(json.dumps(out, indent=1) + '\n')
    for r in rows:
        print(f"{r['site']}: nearest water of {tspec['pdb']} HOH{r['nearest_waters'][0]['water']} at "
              f"{r['nearest_waters'][0]['distance_A']} A (occ {r['nearest_waters'][0]['occupancy']}, B {r['nearest_waters'][0]['b_iso']}); "
              f"nearest atom of its crystal ligand {r['nearest_crystal_ligand_atom']}; nearest deposited atom that is not water "
              f"{r['nearest_deposited_atoms_not_water'][0]}")
    print(f"{out['n_sites_present_as_water']} of {len(rows)} sites hold a water in {tspec['pdb']} (CA fit {len(keys)} atoms, {out['ca_rmsd_A']} A)")
    return out


# --------------------------------------------------------------------------------------------------
# the comparison of the runs: dry receptor against the receptor with the waters
# --------------------------------------------------------------------------------------------------
ROT_REF = ('T3_6RL9/runs/sfa_neutral__ad4zn__prod/seed14/out.pdbqt', 3)
# The ring-rotated pose itself: the mode that md/SYSTEM.md section 5 takes as start R (the best-scored
# mode of the reported production run beyond 2.0 A with the aniline N on Asn67 OD1).  The paper's own
# definition of "ring-rotated" is any mode beyond 2.0 A from the crystal ligand; once the waters are
# in, such a mode need not be that pose, so every mode is also measured against this reference.
ANILINE_N = '[#7;$([#7]-c);!$([#7]-[#16])]'
SULF_N = '[#7;$([#7]-[#16])]'


def one_conformer(mol, cid):
    from rdkit import Chem
    m = Chem.Mol(mol)
    m.RemoveAllConformers()
    m.AddConformer(Chem.Conformer(mol.GetConformer(cid)), assignId=True)
    return m


def rotated_reference(A, base):
    """The reference ring-rotated pose in the frame of `base` (6RL9: the stored mode; 3HS4: the same
    mode carried over by the CA superposition the cross-docking analysis uses)."""
    pm = A.load_poses(DOCK / ROT_REF[0])
    ref = one_conformer(pm, pm.GetConformers()[ROT_REF[1] - 1].GetId())
    if base != 'T3_6RL9':
        src, dst = A.ca_atoms(DOCK / 'T3_6RL9' / 'receptor.pdb'), A.ca_atoms(DOCK / base / 'receptor.pdb')
        keys = sorted(set(src) & set(dst))
        R, t = A.kabsch(np.array([src[k] for k in keys]), np.array([dst[k] for k in keys]))
        c = ref.GetConformer()
        for i in range(ref.GetNumAtoms()):
            c.SetAtomPosition(i, (R @ np.array(list(c.GetAtomPosition(i))) + t).tolist())
    return ref


def classify(r_xtal, r_rot):
    if r_xtal <= 1.0:
        return 'near-native'
    if r_xtal <= 2.0:
        return 'between'
    return 'rotated pose' if r_rot <= 2.0 else 'other, beyond 2.0 A'


def mode_table(A, rundir, ref, rot, sites, tdir):
    """Every mode of every seed: score, RMSD to the crystal ligand and to the rotated reference,
    zinc and Thr199 distances, aniline N to Asn67 OD1, nearest water site."""
    from rdkit import Chem
    res = json.loads((rundir / 'results.json').read_text())
    rows_rec = A.receptor_atoms(tdir)
    pick = lambda rn, ri, an: [x for (a, b, c, x) in rows_rec if a == rn and (ri is None or b == ri) and c == an]  # noqa: E731
    zn, og1, od1 = pick('ZN', None, 'ZN'), pick('THR', 199, 'OG1'), pick('ASN', 67, 'OD1')
    W = np.array(list(sites.values()))
    out = []
    for s in res['seeds']:
        if not s.get('ok'):
            continue
        pm = A.load_poses(DOCK / s['out'])
        cids = [c.GetId() for c in pm.GetConformers()]
        rx = A.symm_rmsd(ref, pm, cids)
        rr = A.symm_rmsd(rot, pm, cids)
        heavy = [a.GetIdx() for a in pm.GetAtoms() if a.GetAtomicNum() > 1]
        nan = [m[0] for m in pm.GetSubstructMatches(Chem.MolFromSmarts(ANILINE_N))]
        nsu = [m[0] for m in pm.GetSubstructMatches(Chem.MolFromSmarts(SULF_N))]
        for k, cid in enumerate(cids):
            X = np.array(pm.GetConformer(cid).GetPositions())
            dw = np.linalg.norm(X[heavy][:, None] - W[None], axis=2)
            out.append(dict(seed=int(s['seed']), mode=k + 1, score=s['scores'][k],
                            rmsd_crystal=round(rx[k], 3), rmsd_rotated_ref=round(rr[k], 3),
                            cls=classify(rx[k], rr[k]),
                            zn_n=round(min(dist(X[i], z) for i in nsu for z in zn), 2) if nsu and zn else None,
                            thr199_og1_n=round(min(dist(X[i], z) for i in nsu for z in og1), 2) if nsu and og1 else None,
                            aniline_n_asn67_od1=round(min(dist(X[i], z) for i in nan for z in od1), 2) if nan and od1 else None,
                            min_dist_water_site=round(float(dw.min()), 2),
                            water_sites_within_2p4A=int((dw.min(0) < 2.4).sum())))
    return out


def summarise_modes(modes, window=2.0):
    """Counts and scores per orientation, by the paper's definition (near-native <= 1.0 A,
    ring-rotated > 2.0 A) and for the rotated pose proper (within 2.0 A of the reference)."""
    seeds = sorted({m['seed'] for m in modes})
    by = {s: [m for m in modes if m['seed'] == s] for s in seeds}
    paper = lambda m: 'near-native' if m['rmsd_crystal'] <= 1.0 else ('ring-rotated' if m['rmsd_crystal'] > 2.0 else 'between')  # noqa: E731
    top = [by[s][0] for s in seeds]
    out = dict(n_seeds=len(seeds), n_modes=len(modes),
               top_pose=dict(near_native=sum(paper(m) == 'near-native' for m in top),
                             ring_rotated_paper=sum(paper(m) == 'ring-rotated' for m in top),
                             between=sum(paper(m) == 'between' for m in top),
                             rotated_pose=sum(m['cls'] == 'rotated pose' for m in top),
                             other_beyond_2A=sum(m['cls'] == 'other, beyond 2.0 A' for m in top)),
               modes=dict(near_native=sum(paper(m) == 'near-native' for m in modes),
                          ring_rotated_paper=sum(paper(m) == 'ring-rotated' for m in modes),
                          between=sum(paper(m) == 'between' for m in modes),
                          rotated_pose=sum(m['cls'] == 'rotated pose' for m in modes),
                          other_beyond_2A=sum(m['cls'] == 'other, beyond 2.0 A' for m in modes)),
               seeds_holding=dict(near_native=sum(any(paper(m) == 'near-native' for m in by[s]) for s in seeds),
                                  ring_rotated_paper=sum(any(paper(m) == 'ring-rotated' for m in by[s]) for s in seeds),
                                  rotated_pose=sum(any(m['cls'] == 'rotated pose' for m in by[s]) for s in seeds)))

    def best(sel):
        v = [m['score'] for m in modes if sel(m)]
        return min(v) if v else None
    out['best_score'] = dict(near_native=best(lambda m: paper(m) == 'near-native'),
                             ring_rotated_paper=best(lambda m: paper(m) == 'ring-rotated'),
                             rotated_pose=best(lambda m: m['cls'] == 'rotated pose'),
                             other_beyond_2A=best(lambda m: m['cls'] == 'other, beyond 2.0 A'),
                             any=min(m['score'] for m in modes))
    gaps = dict(ring_rotated_paper=[], rotated_pose=[])
    for s in seeds:
        nn = [m['score'] for m in by[s] if paper(m) == 'near-native']
        for key, sel in (('ring_rotated_paper', lambda m: paper(m) == 'ring-rotated'),
                         ('rotated_pose', lambda m: m['cls'] == 'rotated pose')):
            ro = [m['score'] for m in by[s] if sel(m)]
            if nn and ro:
                gaps[key].append(round(min(ro) - min(nn), 3))
    # per seed: best mode of that kind minus best near-native mode (positive: near-native scores better)
    out['gap_to_best_near_native'] = {k: (dict(n_seeds=len(v), mean=round(float(np.mean(v)), 3), min=min(v), max=max(v))
                                          if v else None) for k, v in gaps.items()}
    win = dict(ring_rotated_paper=[], rotated_pose=[])
    for s in seeds:
        b = by[s][0]['score']
        for m in by[s]:
            if m['score'] - b <= window:
                if paper(m) == 'ring-rotated':
                    win['ring_rotated_paper'].append((s, m['mode']))
                if m['cls'] == 'rotated pose':
                    win['rotated_pose'].append((s, m['mode']))
    out[f'modes_within_{window:g}_kcal_of_the_seed_best'] = {
        k: dict(n_modes=len(v), n_seeds=len({a for a, _ in v})) for k, v in win.items()}
    out['largest_score_span_of_a_seed'] = round(max(by[s][-1]['score'] - by[s][0]['score'] for s in seeds), 3)
    out['min_dist_water_site'] = dict(
        top_poses=[min(m['min_dist_water_site'] for m in top), max(m['min_dist_water_site'] for m in top)],
        all_modes=min(m['min_dist_water_site'] for m in modes),
        modes_within_2p4A_of_a_site=sum(m['water_sites_within_2p4A'] > 0 for m in modes))
    return out


def families(A, rundir, modes, limit=1.0):
    """Leader clusters (2.0 A, score order) of all modes of all seeds that lie beyond `limit` from
    the crystal ligand: what the search finds besides the crystal pose."""
    import production_analysis as P
    from rdkit import Chem
    res = json.loads((rundir / 'results.json').read_text())
    pool, meta = None, []
    for s in res['seeds']:
        if not s.get('ok'):
            continue
        pm = A.load_poses(DOCK / s['out'])
        for k, c in enumerate(pm.GetConformers()):
            m = next(x for x in modes if x['seed'] == s['seed'] and x['mode'] == k + 1)
            if m['rmsd_crystal'] <= limit:
                continue
            if pool is None:
                pool = one_conformer(pm, c.GetId())
            else:
                pool.AddConformer(Chem.Conformer(c), assignId=True)
            meta.append(m)
    if pool is None:
        return []
    cids = [c.GetId() for c in pool.GetConformers()]
    D = P.pairwise_symm_rmsd(pool, cids)
    best_all = min(x['score'] for x in modes)
    out = []
    for L, mem in P.leader_clusters([m['score'] for m in meta], D):
        r = meta[L]
        out.append(dict(n_modes=len(mem), n_seeds=len({meta[i]['seed'] for i in mem}),
                        n_top_poses=sum(meta[i]['mode'] == 1 for i in mem),
                        best_score=r['score'], behind_best_of_run=round(r['score'] - best_all, 3),
                        representative=dict(seed=r['seed'], mode=r['mode']), rmsd_crystal=r['rmsd_crystal'],
                        rmsd_rotated_ref=r['rmsd_rotated_ref'], cls=r['cls'], zn_n=r['zn_n'],
                        thr199_og1_n=r['thr199_og1_n'], aniline_n_asn67_od1=r['aniline_n_asn67_od1'],
                        min_dist_water_site=r['min_dist_water_site']))
    return out


def water_sites(variant):
    v = json.loads((DOCK / variant / 'variant.json').read_text())
    return {n: np.array(w['O']) for n, w in v['waters'].items()}, v


def compare(variant, scratch, do_pb=True, only=None):
    """Every run of the water variant against the run of the same name in its dry source."""
    import analyse as A
    import glu_protonation_variant as G
    scratch = pathlib.Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    A.ANALYSIS = scratch                       # PoseBusters input files go to the scratch directory
    sites, vj = water_sites(variant)
    dry = vj['source']
    base, _ = base_and_spec(variant)
    rot = rotated_reference(A, base)
    orig = A.transform_reference               # a cross-docking reference of an existing directory is not rewritten

    def transform_reference(src_dir, dst_dir, ref, out_sdf):
        if pathlib.Path(dst_dir).name != variant:
            out_sdf = scratch / f'{pathlib.Path(dst_dir).name}__{pathlib.Path(out_sdf).name}'
        return orig(src_dir, dst_dir, ref, out_sdf)
    A.transform_reference = transform_reference
    out_json = DOCK / variant / 'comparison_water.json'
    result = json.loads(out_json.read_text()) if out_json.exists() else {}
    result.update(variant=variant, dry=dry, rule=vj['rule'], waters=sorted(sites),
                  rotated_reference=dict(file=ROT_REF[0], mode=ROT_REF[1], frame=base),
                  definitions='near-native: within 1.0 A of the crystal ligand; ring-rotated (paper): beyond 2.0 A; '
                              'rotated pose: beyond 2.0 A and within 2.0 A of the rotated reference; RMSDs symmetry-corrected, '
                              'heavy atoms, receptor frame; water sites: the oxygen positions of the kept waters')
    runs = result.setdefault('runs', {})
    for rj in sorted((DOCK / variant / 'runs').glob('*/results.json')):
        tag = rj.parent.name
        if only and tag not in only:
            continue
        ddir = DOCK / dry / 'runs' / tag
        if not (ddir / 'results.json').exists():
            runs[tag] = dict(error=f'no run of this name in {dry}')
            continue
        has_gate = len(tag.split('__')) == 2
        entry = {}
        for which, rd, tname in (('dry', ddir, dry), ('water', rj.parent, variant)):
            o, res = G.analyse_run(A, rd, do_pb and has_gate and which == 'water')
            if has_gate and which == 'dry':
                # the PoseBusters figure of a dry run: the stored validation analysis or dock/GLU106.md
                for src_json, key in ((DOCK / 'analysis' / 'summary.json', None),
                                      (DOCK / dry / 'comparison.json', 'variant')):
                    if not src_json.exists():
                        continue
                    j = json.loads(src_json.read_text())
                    if key is None and f'{dry}/{tag}' in j:
                        o['posebusters_physics_pass_fraction'] = j[f'{dry}/{tag}']['stats'].get('posebusters_physics_pass_fraction')
                        o['posebusters_source'] = 'dock/analysis/summary.json'
                    elif key and tag in j.get('runs', {}):
                        o['posebusters_physics_pass_fraction'] = j['runs'][tag][key].get('posebusters_physics_pass_fraction')
                        o['posebusters_source'] = f'dock/{dry}/comparison.json'
            if 'rmsd' in o:
                o['gate'] = G.gate(o)
            ref, _, _ = A.pick_reference(tname, res['ligand'], A.load_poses(DOCK / res['seeds'][0]['out']))
            modes = mode_table(A, rd, ref, rot, sites, DOCK / tname)
            o['mode_summary'] = summarise_modes(modes)
            o['families_beyond_1A'] = families(A, rd, modes)
            o['modes'] = modes
            entry[which] = o
        d, w = entry['dry'], entry['water']
        same_in = all((DOCK / pathlib.Path(a['out']).parent / 'ligand.pdbqt').read_bytes()
                      == (DOCK / pathlib.Path(b['out']).parent / 'ligand.pdbqt').read_bytes()
                      for a, b in zip(json.loads(rj.read_text())['seeds'], json.loads((ddir / 'results.json').read_text())['seeds']))
        entry['check'] = dict(ligand_inputs_identical_all_seeds=bool(same_in),
                              settings_that_differ={k: (d['settings'][k], w['settings'][k]) for k in d['settings']
                                                    if d['settings'][k] != w['settings'][k] and k != 'maps'})
        runs[tag] = entry
        ms_d, ms_w = d['mode_summary'], w['mode_summary']
        print(f"{variant}/{tag}: RMSD {d['rmsd']['mean']:.2f} -> {w['rmsd']['mean']:.2f}; top near-native "
              f"{ms_d['top_pose']['near_native']} -> {ms_w['top_pose']['near_native']} of {ms_w['n_seeds']}; modes of the rotated pose "
              f"{ms_d['modes']['rotated_pose']} -> {ms_w['modes']['rotated_pose']}; modes beyond 2.0 A "
              f"{ms_d['modes']['ring_rotated_paper']} -> {ms_w['modes']['ring_rotated_paper']}; score {d['score']['mean']:.3f} -> "
              f"{w['score']['mean']:.3f}", flush=True)
    out_json.write_text(json.dumps(result, indent=1) + '\n')
    print(f'wrote {out_json}')
    return 0


def rescore(scratch, variants=None, tags=('sfa_neutral__ad4zn__prod', 'sfa_neutral__ad4zn')):
    """Stored modes of the dry runs scored without a search (vina --score_only --scoring ad4) in the
    dry maps and in the maps of every water variant: per seed the best near-native mode and the best
    mode of the rotated pose; plus the rotated reference itself."""
    import re
    import subprocess
    import analyse as A
    scratch = pathlib.Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    def models(path):
        out, cur = [], None
        for ln in path.read_text().splitlines():
            if ln.startswith('MODEL'):
                cur = []
            elif ln.startswith('ENDMDL'):
                out.append(cur)
                cur = None
            elif cur is not None:
                cur.append(ln)
        return out

    def score(lines, maps):
        import os
        f = scratch / f'pose_water_{os.getpid()}.pdbqt'
        f.write_text('\n'.join(lines) + '\n')
        r = subprocess.run(['vina', '--ligand', str(f), '--maps', str(maps), '--scoring', 'ad4', '--score_only',
                            '--cpu', '1'], capture_output=True, text=True)
        return float(re.search(r'Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)', r.stdout).group(1))
    for vdir in sorted(DOCK.glob('T*wat*')):
        if not (vdir / 'variant.json').exists() or not list((vdir / 'ad4zn').glob('receptor_tz.*.map')):
            continue
        if variants and vdir.name not in variants:
            continue
        variant = vdir.name
        sites, vj = water_sites(variant)
        dry = vj['source']
        base, _ = base_and_spec(variant)
        rot = rotated_reference(A, base)
        out = dict(variant=variant, dry=dry, what='vina --score_only --scoring ad4 on stored modes of the dry runs; near-native: '
                   'RMSD <= 1.0 A; rotated pose: > 2.0 A from the crystal ligand and <= 2.0 A from the rotated reference; '
                   'best mode of each kind per seed', runs={})
        for tag in tags:
            rd = DOCK / dry / 'runs' / tag
            if not (rd / 'results.json').exists():
                continue
            res = json.loads((rd / 'results.json').read_text())
            ref, _, _ = A.pick_reference(dry, res['ligand'], A.load_poses(DOCK / res['seeds'][0]['out'])) \
                if base == 'T3_6RL9' else (__import__('rdkit.Chem', fromlist=['Chem']).MolFromMolFile(
                    str(DOCK / base / 'ref_cross_T3_6RL9_SAN.sdf')), None, None)
            rows = []
            for s in res['seeds']:
                sd = DOCK / pathlib.Path(s['out']).parent
                pm = A.load_poses(sd / 'out.pdbqt')
                cids = [c.GetId() for c in pm.GetConformers()]
                rx, rr = A.symm_rmsd(ref, pm, cids), A.symm_rmsd(rot, pm, cids)
                ms = models(sd / 'out.pdbqt')
                for kind, idx in (('near-native', [i for i, v in enumerate(rx) if v <= 1.0][:1]),
                                  ('rotated pose', [i for i, v in enumerate(rx) if v > 2.0 and rr[i] <= 2.0][:1])):
                    for i in idx:
                        rows.append(dict(seed=int(s['seed']), kind=kind, mode=i + 1, rmsd_crystal=round(rx[i], 3),
                                         stored_score=s['scores'][i],
                                         dry_maps=score(ms[i], DOCK / dry / 'ad4zn' / 'receptor_tz'),
                                         water_maps=score(ms[i], vdir / 'ad4zn' / 'receptor_tz')))
            summ = {}
            for kind in ('near-native', 'rotated pose'):
                x = [r for r in rows if r['kind'] == kind]
                if x:
                    w = np.array([r['water_maps'] for r in x])
                    ch = np.array([r['water_maps'] - r['dry_maps'] for r in x])
                    summ[kind] = dict(n=len(x), water_maps_mean=round(float(w.mean()), 3), water_maps_min=float(w.min()),
                                      water_maps_max=float(w.max()), change_mean=round(float(ch.mean()), 3),
                                      change_min=round(float(ch.min()), 3), change_max=round(float(ch.max()), 3),
                                      max_abs_dry_minus_stored=round(max(abs(r['dry_maps'] - r['stored_score']) for r in x), 3))
            out['runs'][tag] = dict(summary=summ, modes=rows)
            print(f"rescore {variant} <- {dry}/{tag}: " + '; '.join(
                f"{k}: n={v['n']}, in the water maps {v['water_maps_mean']:+.3f} ({v['water_maps_min']:+.3f} to {v['water_maps_max']:+.3f}), "
                f"change {v['change_mean']:+.3f}" for k, v in summ.items()), flush=True)
        (vdir / 'rescore_fixed_poses.json').write_text(json.dumps(out, indent=1) + '\n')
    return 0


# --------------------------------------------------------------------------------------------------
# the tables of dock/WATER4.md
# --------------------------------------------------------------------------------------------------
RUN_LABEL = {'sfa_neutral__ad4zn__prod': 'sulfanilamide, AD4Zn, DFT geometry (production)',
             'sfa_neutral__ad4zn': 'sulfanilamide, AD4Zn, ETKDG (validation)',
             'sfa_anion__ad4zn': 'sulfonamidate, AD4Zn, ETKDG (validation)',
             'sfa_neutral__vina': 'sulfanilamide, Vina, ETKDG (validation)',
             'sfa_neutral__vinardo': 'sulfanilamide, Vinardo, ETKDG (validation)',
             'sfa_neutral__vina__prod': 'sulfanilamide, Vina, DFT geometry (production)',
             'sfa_neutral__vinardo__prod': 'sulfanilamide, Vinardo, DFT geometry (production)'}
RULE_LABEL = dict(hbond='hydrogen-bond network', rev='network, ring reversed', bare='bare oxygen (precedent)',
                  fixed='one fixed orientation')


def tables(variants):
    """Markdown tables of dock/WATER4.md from comparison_water.json / rescore_fixed_poses.json."""
    import decimal
    M = '−'

    def f(x, n=2):
        if x is None:
            return '–'
        q = decimal.Decimal(str(round(x, 6))).quantize(decimal.Decimal(1).scaleb(-n), rounding=decimal.ROUND_HALF_UP)
        return f'{q:.{n}f}'.replace('-', M)
    pm = lambda o, n=2: f"{f(o['mean'], n)} ± {f(o['sd'], n)}"  # noqa: E731
    yn = lambda g: '–' if g is None else ('pass' if g else 'fail')  # noqa: E731
    data = {v: json.loads((DOCK / v / 'comparison_water.json').read_text()) for v in variants
            if (DOCK / v / 'comparison_water.json').exists()}
    glu = lambda v: 'carboxylate' if 'glu106GLU' in v else 'acid (reported)'  # noqa: E731
    order = list(RUN_LABEL)

    def each():
        for v, d in data.items():
            for tag in sorted(d['runs'], key=lambda t: order.index(t) if t in order else 99):
                r = d['runs'][tag]
                if 'error' not in r:
                    yield v, d, tag, r

    print('\n#### The waters and their hydrogen bonds, by hydrogen rule\n')
    print('| Variant (rule) | Water | B (\u00c5\u00b2) | Donates to: acceptor, H\u00b7\u00b7\u00b7A (\u00c5), O\u2013H\u00b7\u00b7\u00b7A (\u00b0) | '
          'Accepts from: donor, H\u00b7\u00b7\u00b7O (\u00c5), D\u2013H\u00b7\u00b7\u00b7O (\u00b0) | Hydrogens without a partner | '
          'Hydrogen pairs closer than 2.0 \u00c5 |')
    print('|:--|:--|--:|:--|:--|--:|:--|')
    for v in variants:
        vj = json.loads((DOCK / v / 'variant.json').read_text())
        for k, (n, r) in enumerate(vj['hydrogen_bonds'].items()):
            nh, dots = len(vj['waters'][n]['H']), '\u00b7\u00b7\u00b7'
            none = '\u2013'
            print(f"| {'`' + v + '` (' + RULE_LABEL[vj['rule']] + ')' if k == 0 else ''} | {n} | {vj['waters'][n]['b_iso']} | "
                  + ('; '.join(f"{x['to']}, {f(x['H_A'])}, {f(x['angle'], 1)}" for x in r['donates']) or (none if nh else 'no hydrogens'))
                  + ' | ' + ('; '.join(f"{x['frm']}, {f(x['H_O'])}, {f(x['angle'], 1)}" for x in r['accepts']) or none)
                  + f" | {r['hydrogens_without_partner'] if nh else none} | "
                  + ('; '.join(f"{x['hydrogen']}{dots}{x['to']} {f(x['distance'])}" for x in r['close_H_H']) or '\u2013') + ' |')

    print('\n#### Orientations: every mode of every seed\n')
    print('| Run | Glu106 | Receptor (hydrogen rule) | Seeds | Top pose: near-native / rotated pose / other beyond 2.0 Å / '
          'between | All modes: near-native / rotated pose / other beyond 2.0 Å / between | Seeds holding a mode of the '
          'rotated pose | Best score: near-native | Best score: rotated pose | Rotated pose behind the best near-native mode of '
          'its seed, min–max | Modes beyond 2.0 Å within 2 kcal mol⁻¹ of the seed best (of them the rotated pose) | '
          'Closest heavy atom of any mode to a water site (Å) |')
    print('|:--|:--|:--|--:|:--|:--|--:|--:|--:|:--|:--|--:|')
    for v, d, tag, r in each():
        for which in ('dry', 'water'):
            ms = r[which]['mode_summary']
            t, m = ms['top_pose'], ms['modes']
            g = ms['gap_to_best_near_native']['rotated_pose']
            w = ms['modes_within_2_kcal_of_the_seed_best']
            rec = d['dry'] if which == 'dry' else f"{v} ({RULE_LABEL[d['rule']]})"
            print(f"| {RUN_LABEL.get(tag, tag) if which == 'dry' else ''} | {glu(v) if which == 'dry' else ''} | "
                  f"{'dry: ' if which == 'dry' else 'water: '}`{rec.split(' ')[0]}`{rec[len(rec.split(' ')[0]):]} | {ms['n_seeds']} | "
                  f"{t['near_native']} / {t['rotated_pose']} / {t['other_beyond_2A']} / {t['between']} | "
                  f"{m['near_native']} / {m['rotated_pose']} / {m['other_beyond_2A']} / {m['between']} | "
                  f"{ms['seeds_holding']['rotated_pose']} | {f(ms['best_score']['near_native'], 3)} | "
                  f"{f(ms['best_score']['rotated_pose'], 3)} | "
                  + (f"{f(g['min'], 3)}–{f(g['max'], 3)}" if g else '–')
                  + f" | {w['ring_rotated_paper']['n_modes']} ({w['rotated_pose']['n_modes']}) | "
                  f"{f(ms['min_dist_water_site']['all_modes'])} |")

    print('\n#### Does the crystal pose still dock: top poses, contacts, gates\n')
    print('| Run | Glu106 | Receptor | RMSD top (Å) | best of seeds | f ≤ 1 Å (≤ 2 Å) | Score (kcal mol⁻¹) | '
          'Zn–N (Å) / pass | Thr199 OG1···N (Å) / pass | Best mode of each seed (Å), min–max | '
          'PoseBusters physics | Strict gate / 2.0 Å |')
    print('|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|')
    for v, d, tag, r in each():
        for which in ('dry', 'water'):
            o = r[which]
            sc = o.get('site_checks') or {}
            site = lambda c: f"{f(c.get('distance_mean'))} / {f(c.get('pass_fraction'))}"  # noqa: E731
            nd = 1 if o['n_seeds'] == 10 else 2
            bm = o['rmsd']['best_of_modes_per_seed']
            has_gate = len(tag.split('__')) == 2
            gt = f"{yn(o['gate']['passes_G2'])} / {yn(o['gate']['passes_2A_criterion'])}" if has_gate else '–'
            print(f"| {RUN_LABEL.get(tag, tag) if which == 'dry' else ''} | {glu(v) if which == 'dry' else ''} | "
                  f"{'dry' if which == 'dry' else 'water, ' + RULE_LABEL[d['rule']]} | {pm(o['rmsd'])} | "
                  f"{f(o['rmsd']['best_of_seeds'])} | {f(o['rmsd']['f_le_1A'], nd)} ({f(o['rmsd']['f_le_2A'], nd)}) | "
                  f"{pm(o['score'], 3)} | {site(sc.get('Zn-N(sulfonamide)', {}))} | {site(sc.get('Thr199 OG1 H-bond', {}))} | "
                  f"{f(min(bm))}–{f(max(bm))} | {f(o.get('posebusters_physics_pass_fraction'), 1) if has_gate else '–'} | {gt} |")

    print('\n#### What the search finds besides the crystal pose: families of modes beyond 1.0 Å (2.0 Å clusters, all seeds pooled)\n')
    print('| Run | Glu106 | Receptor | Family | Modes (seeds; top poses) | Best score | Behind the best mode of the run | RMSD to the '
          'crystal ligand (Å) | RMSD to the rotated reference (Å) | Zn–N (Å) | Aniline N···Asn67 OD1 (Å) | '
          'Nearest water site (Å) |')
    print('|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|')
    for v, d, tag, r in each():
        if tag.split('__')[1] != 'ad4zn':
            continue
        for which in ('dry', 'water'):
            for k, fam in enumerate(r[which]['families_beyond_1A']):
                print(f"| {RUN_LABEL.get(tag, tag) if which == 'dry' and k == 0 else ''} | {glu(v) if which == 'dry' and k == 0 else ''} | "
                      f"{('dry' if which == 'dry' else 'water, ' + RULE_LABEL[d['rule']]) if k == 0 else ''} | {k + 1}: {fam['cls']} | "
                      f"{fam['n_modes']} ({fam['n_seeds']}; {fam['n_top_poses']}) | {f(fam['best_score'], 3)} | "
                      f"{f(fam['behind_best_of_run'], 3)} | {f(fam['rmsd_crystal'])} | {f(fam['rmsd_rotated_ref'])} | {f(fam['zn_n'])} | "
                      f"{f(fam['aniline_n_asn67_od1'])} | {f(fam['min_dist_water_site'])} |")

    print('\n#### Stored modes of the dry runs scored in the maps with the waters (no search)\n')
    print('| Maps | Stored dry run | Near-native modes: n, score in the dry maps → in the water maps (change) | Modes of the rotated pose: '
          'n, score in the water maps, mean (min to max) |')
    print('|:--|:--|:--|:--|')
    for v in variants:
        pth = DOCK / v / 'rescore_fixed_poses.json'
        if not pth.exists():
            continue
        rj = json.loads(pth.read_text())
        for tag, x in rj['runs'].items():
            nn, ro = x['summary'].get('near-native'), x['summary'].get('rotated pose')
            dry_mean = float(np.mean([m['dry_maps'] for m in x['modes'] if m['kind'] == 'near-native'])) if nn else None
            print(f"| `{v}` | `{rj['dry']}`, {RUN_LABEL.get(tag, tag)} | "
                  + (f"{nn['n']}, {f(dry_mean, 3)} → {f(nn['water_maps_mean'], 3)} ({f(nn['change_min'], 3)} to {f(nn['change_max'], 3)})"
                     if nn else '–') + ' | '
                  + (f"{ro['n']}, +{ro['water_maps_mean']:.0f} (+{ro['water_maps_min']:.0f} to +{ro['water_maps_max']:.0f})"
                     if ro else '–') + ' |')
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default=None, help='an existing prepared receptor directory name')
    ap.add_argument('--waters', type=int, nargs='+', default=None, help='residue numbers of the waters')
    ap.add_argument('--rule', default='hbond', choices=['hbond', 'rev', 'bare', 'fixed'])
    ap.add_argument('--name', default=None, help='name of the variant directory (default <target>_wat<N>[rule])')
    ap.add_argument('--map-sites-from', default=None, help='report where the water sites of this target lie in --target')
    ap.add_argument('--out', default=None, help='--map-sites-from: write the result to this JSON file')
    ap.add_argument('--compare', action='store_true', help='compare the runs of the water variant --target with its dry source')
    ap.add_argument('--compare-all', action='store_true', help='--compare for every water variant under dock/')
    ap.add_argument('--rescore', action='store_true', help='stored dry modes scored in the maps of every water variant')
    ap.add_argument('--variants', nargs='*', default=None, help='--rescore: only these variants (--only: only these run names)')
    ap.add_argument('--scratch', default=None, help='directory for temporary files (--compare, --rescore)')
    ap.add_argument('--tables', nargs='*', default=None, metavar='VARIANT', help='print the tables of dock/WATER4.md for these variants')
    ap.add_argument('--no-posebusters', action='store_true')
    ap.add_argument('--only', nargs='*', default=None, help='--compare: only these run names')
    a = ap.parse_args()
    if a.tables is not None:
        return tables(a.tables)
    if a.compare or a.compare_all or a.rescore:
        if not a.scratch:
            raise SystemExit('--scratch <dir> is required (nothing is written into existing directories)')
        if a.rescore:
            return rescore(a.scratch, a.variants, tuple(a.only) if a.only else ('sfa_neutral__ad4zn__prod', 'sfa_neutral__ad4zn'))
        todo = [a.target] if a.compare else sorted(d.name for d in DOCK.glob('T*wat*') if (d / 'runs').is_dir())
        for v in todo:
            compare(v, a.scratch, not a.no_posebusters, a.only)
        return 0
    if a.map_sites_from:
        map_sites(a.target, a.map_sites_from, a.waters, a.out)
        return 0
    return make_variant(a.target, a.waters, a.rule, a.name)


if __name__ == '__main__':
    sys.exit(main())
