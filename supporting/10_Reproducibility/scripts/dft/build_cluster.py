"""
Build a chemically meaningful hydrogen-bonded cluster from a CIF.

Why this exists: for 3ATB the asymmetric unit happened to BE the hydrogen-bonded ion pair (the CIF lists
N4-H4...O11 with symmetry code '.'), so the atoms could be used directly.  That is not true for the
piperazinium/benzotriazolium/benzimidazolium picrate salts: their N-H...A bonds point at symmetry-related
neighbours, so the asymmetric unit is a set of fragments that are NOT bonded to each other.  Optimising it
as-is would be meaningless.  This script expands the crystal symmetry and picks, for each moiety of the
formula unit, the periodic image that is actually hydrogen-bonded to the growing cluster.

Usage: python dft/build_cluster.py <cif> <id> [--core <label>] [--charge N]
Writes molecules/<id>/input/<id>_cluster.xyz and molecules/<id>/cases/<id>_gas.json
"""
import sys, os, re, json, math, itertools, pathlib
import numpy as np

COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'Zn': 1.22}
DONOR = {'N', 'O', 'S'}                 # X-H donors we look for
ACCEPTOR = {'N', 'O', 'S', 'Cl', 'Br', 'F'}
HB_MAX_HA = 2.65                        # H...A cutoff, Angstrom
HB_MIN_ANG = 130.0                      # D-H...A angle cutoff, degrees
CLASH_MIN = 2.40                        # reject images with heavy atoms closer than this

# ---------------------------------------------------------------- CIF parsing
def read_cif(path):
    txt = open(path).read()
    def num(tag):
        m = re.search(rf'^{re.escape(tag)}\s+(\S+)', txt, re.M)
        return float(m.group(1).split('(')[0])
    cell = [num(f'_cell_{k}') for k in ('length_a', 'length_b', 'length_c',
                                        'angle_alpha', 'angle_beta', 'angle_gamma')]
    lines = txt.splitlines()
    atoms, ops = [], []
    i = 0
    while i < len(lines):
        if lines[i].strip() != 'loop_':
            i += 1; continue
        j = i + 1; tags = []
        while j < len(lines) and lines[j].strip().startswith('_'):
            tags.append(lines[j].strip().split()[0]); j += 1
        body = []
        while j < len(lines):
            s = lines[j].strip()
            if not s or s.startswith(('loop_', '_', '#')): break
            body.append(s); j += 1
        if any(t.endswith('_fract_x') for t in tags):
            ix = {t: k for k, t in enumerate(tags)}
            for row in body:
                p = row.split()
                if len(p) < len(tags): continue
                lab = p[ix['_atom_site_label']]
                sym = (p[ix['_atom_site_type_symbol']] if '_atom_site_type_symbol' in ix
                       else re.match(r'[A-Za-z]+', lab).group(0))
                f = [float(p[ix[f'_atom_site_fract_{c}']].split('(')[0]) for c in 'xyz']
                atoms.append({'label': lab.replace('(', '').replace(')', ''),
                              'el': sym.capitalize(), 'f': np.array(f)})
        if any('symop' in t or 'symmetry_equiv_pos' in t for t in tags):
            ops = [r.strip().strip("'\"") for r in body]
        i = j
    if not ops: ops = ['x,y,z']
    return cell, atoms, ops

def frac_to_cart_matrix(cell):
    a, b, c, al, be, ga = cell
    al, be, ga = map(math.radians, (al, be, ga))
    v = math.sqrt(1 - math.cos(al)**2 - math.cos(be)**2 - math.cos(ga)**2
                  + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    return np.array([[a, b * math.cos(ga), c * math.cos(be)],
                     [0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)],
                     [0, 0, c * v / math.sin(ga)]])

def parse_op(op):
    """'-x,1/2+y,1/2-z' -> (3x3 rotation, 3 translation) acting on fractional coordinates."""
    R = np.zeros((3, 3)); t = np.zeros(3)
    for r, part in enumerate(op.lower().replace(' ', '').split(',')):
        for m in re.finditer(r'([+-]?)(\d+/\d+|\d*\.\d+|\d+)?\*?([xyz])?', part):
            sgn, val, var = m.group(1), m.group(2), m.group(3)
            if not val and not var: continue
            s = -1.0 if sgn == '-' else 1.0
            if var:
                coef = 1.0
                if val: coef = eval(val) if '/' in val else float(val)
                R[r, 'xyz'.index(var)] += s * coef
            else:
                t[r] += s * (eval(val) if '/' in val else float(val))
    return R, t

# ---------------------------------------------------------------- geometry helpers
def bonds(els, X, scale=1.25):
    out = []
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            if np.linalg.norm(X[i] - X[j]) < scale * (COV.get(els[i], 0.8) + COV.get(els[j], 0.8)):
                out.append((i, j))
    return out

def fragments(els, X):
    n = len(els); par = list(range(n))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i, j in bonds(els, X):
        par[find(i)] = find(j)
    g = {}
    for i in range(n): g.setdefault(find(i), []).append(i)
    return sorted(g.values(), key=len, reverse=True)

def angle(a, b, c):
    v1, v2 = a - b, c - b
    cs = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cs))))

def find_hbonds(els, X, groupA, groupB):
    """D-H...A hydrogen bonds with D,H in one group and A in the other (both directions)."""
    bs = bonds(els, X)
    hpar = {}
    for i, j in bs:
        if els[i] == 'H' and els[j] in DONOR: hpar[i] = j
        if els[j] == 'H' and els[i] in DONOR: hpar[j] = i
    out = []
    for src, dst in ((groupA, groupB), (groupB, groupA)):
        for h, d in hpar.items():
            if h not in src: continue
            for a in dst:
                if els[a] not in ACCEPTOR: continue
                dha = np.linalg.norm(X[h] - X[a])
                if dha < HB_MAX_HA and angle(X[d], X[h], X[a]) > HB_MIN_ANG:
                    out.append({'D': d, 'H': h, 'A': a, 'H_A': float(dha),
                                'D_A': float(np.linalg.norm(X[d] - X[a])),
                                'angle': float(angle(X[d], X[h], X[a]))})
    return out

def rings(els, X, size=(5, 6)):
    """Small rings of heavy atoms (for the ELF/LOL plane definitions)."""
    heavy = [i for i, e in enumerate(els) if e != 'H']
    adj = {i: set() for i in heavy}
    for i, j in bonds(els, X):
        if i in adj and j in adj: adj[i].add(j); adj[j].add(i)
    found = []
    for start in heavy:
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            for nxt in adj[cur]:
                if nxt == start and len(path) in size:
                    key = frozenset(path)
                    if key not in [f[0] for f in found]: found.append((key, list(path)))
                elif nxt not in path and len(path) < max(size):
                    stack.append((nxt, path + [nxt]))
    return [p for _, p in found]

# ---------------------------------------------------------------- cluster assembly
def strip(l):
    return l.split('#')[0]


def build(cif, mid, core_label=None, charge=0, out_root=None, require=None):
    cell, atoms, ops = read_cif(cif)
    M = frac_to_cart_matrix(cell)
    els0 = [a['el'] for a in atoms]
    F0 = np.array([a['f'] for a in atoms])
    X0 = F0 @ M.T
    frs = fragments(els0, X0)
    print(f'[cluster] {mid}: {len(atoms)} atoms in the asymmetric unit, {len(frs)} fragments, {len(ops)} symmetry operations')
    for k, f in enumerate(frs):
        comp = {}
        for i in f: comp[els0[i]] = comp.get(els0[i], 0) + 1
        print(f'          frag {k}: {"".join(f"{e}{n}" for e, n in sorted(comp.items())):12s} '
              f'({len(f)} atoms) {",".join(atoms[i]["label"] for i in f[:6])}...')

    # all periodic images of every fragment: (frag index, op index, lattice shift)
    images = []
    for fi, f in enumerate(frs):
        for oi, op in enumerate(ops):
            R, t = parse_op(op)
            for sh in itertools.product((-1, 0, 1), repeat=3):
                Ff = F0[f] @ R.T + t + np.array(sh)
                images.append({'frag': fi, 'op': oi, 'shift': sh, 'F': Ff, 'X': Ff @ M.T,
                               'els': [els0[i] for i in f], 'labels': [atoms[i]['label'] for i in f],
                               'code': f'{oi+1}_{"".join(str(s+5) for s in sh)}'})

    def score_add(cur_els, cur_X, img):
        els = cur_els + img['els']
        X = np.vstack([cur_X, img['X']])
        gA = list(range(len(cur_els))); gB = list(range(len(cur_els), len(els)))
        # clash check on heavy atoms
        for i in gA:
            if els[i] == 'H': continue
            for j in gB:
                if els[j] == 'H': continue
                if np.linalg.norm(X[i] - X[j]) < CLASH_MIN: return None, None
        hb = find_hbonds(els, X, gA, gB)
        if not hb: return 0.0, []
        return sum(HB_MAX_HA - h['H_A'] for h in hb) + 0.5 * len(hb), hb

    # ---- exhaustive assembly search -------------------------------------------------
    # Greedy placement fails: it picks the image that scores best against the PARTIAL cluster and can
    # miss a motif that only pays off once everything is placed (in BTP it dropped the N-H...N proton-
    # sharing dimer).  Instead: reduce each fragment to the handful of images that hydrogen-bond to
    # something, then score every complete assembly.
    def pair_images(fi, base_els, base_X):
        """images of fragment fi that hydrogen-bond to the given placed atoms"""
        out = []
        for img in images:
            if img['frag'] != fi: continue
            s, hb = score_add(base_els, base_X, img)
            if s: out.append((s, hb, img))
        return out

    def lab_of(els, groups, combo, seed, i):
        if i < len(groups[0]):
            return atoms[frs[seed][i]]['label']
        off = len(groups[0])
        for gi, im in enumerate(combo):
            if i < off + len(im['els']):
                return im['labels'][i - off]
            off += len(im['els'])
        return '?'

    best = None
    seeds = range(len(frs)) if core_label is None else [
        next(k for k, f in enumerate(frs) if any(atoms[i]['label'] == core_label for i in f))]
    for seed in seeds:
        seed_els = [els0[i] for i in frs[seed]]
        seed_X = X0[frs[seed]].copy()
        others = [fi for fi in range(len(frs)) if fi != seed]
        # candidate images: bonded to the seed, or bonded to a seed-bonded image of another fragment
        cand = {}
        for fi in others:
            cand[fi] = {im['code']: im for _, _, im in pair_images(fi, seed_els, seed_X)}
        for fi in others:
            for fj in others:
                if fi == fj: continue
                for im2 in list(cand[fj].values()):
                    for _, _, im in pair_images(fi, im2['els'], im2['X']):
                        cand[fi].setdefault(im['code'], im)
        for fi in others:
            if not cand[fi]:                       # nothing hydrogen-bonds: keep the closest image
                cand[fi] = {min(((np.linalg.norm(seed_X[:, None] - im['X'][None], axis=-1).min(), im)
                                 for im in images if im['frag'] == fi), key=lambda z: z[0])[1]['code']:
                            min(((np.linalg.norm(seed_X[:, None] - im['X'][None], axis=-1).min(), im)
                                 for im in images if im['frag'] == fi), key=lambda z: z[0])[1]}
        for combo in itertools.product(*[list(cand[fi].values()) for fi in others]):
            els = list(seed_els); X = seed_X.copy()
            groups = [list(range(len(els)))]
            ok = True
            for im in combo:
                s, hb = score_add(els, X, im)
                if s is None: ok = False; break     # steric clash
                groups.append(list(range(len(els), len(els) + len(im['els']))))
                els = els + im['els']; X = np.vstack([X, im['X']])
            if not ok: continue
            hbs, total = [], 0.0
            for ga, gb in itertools.combinations(groups, 2):
                h = find_hbonds(els, X, ga, gb)
                hbs += h
                total += sum(HB_MAX_HA - x['H_A'] for x in h) + 0.5 * len(h)
            # require the hydrogen-bond graph to connect every fragment, and reward it
            linked = set()
            for ga, gb in itertools.combinations(range(len(groups)), 2):
                if find_hbonds(els, X, groups[ga], groups[gb]): linked |= {ga, gb}
            connected = len(linked) == len(groups)
            if require:                       # motifs the chemist has decided must be present
                base = lambda i: lab_of(els, groups, combo, seed, i)
                have = {tuple(sorted((strip(base(h['D'])), strip(base(h['A']))))) for h in hbs}
                if any(tuple(sorted(r)) not in have for r in require): continue
            score = total + (10.0 if connected else 0.0)
            if best is None or score > best['score']:
                lab = [atoms[i]['label'] for i in frs[seed]]
                parts = [{'frag': seed, 'code': '1_555', 'n': len(frs[seed])}]
                for fi, im in zip(others, combo):
                    lab += [f"{l}#{im['code']}" if im['code'] != '1_555' else l for l in im['labels']]
                    parts.append({'frag': fi, 'code': im['code'], 'n': len(im['els'])})
                best = {'score': score, 'total': total, 'connected': connected, 'els': els, 'X': X,
                        'labels': lab, 'parts': parts, 'hb': hbs, 'seed': seed}

    els, X, labs = best['els'], best['X'], best['labels']
    X = X - X.mean(0)
    print(f"[cluster] assembled {len(els)} atoms, seed fragment {best['seed']}, H-bond score {best['total']:.2f}, "
          f"all fragments hydrogen-bonded: {best['connected']}")
    for p in best['parts']:
        print(f"          fragment {p['frag']} ({p['n']} atoms) via symmetry image {p['code']}")
    frs2 = fragments(els, X)
    hb_all = []
    for ga, gb in itertools.combinations(frs2, 2):
        hb_all += find_hbonds(els, X, ga, gb)
    print(f'[cluster] {len(hb_all)} inter-fragment hydrogen bonds in the cluster:')
    for h in sorted(hb_all, key=lambda h: h['H_A']):
        print(f"          {labs[h['D']]:>10s}-{labs[h['H']]:<8s}...{labs[h['A']]:<10s} "
              f"H...A {h['H_A']:.3f} A  D...A {h['D_A']:.3f} A  angle {h['angle']:.1f} deg")
    return cell, els, X, labs, frs2, hb_all, best

def write_outputs(mid, cell, els, X, labs, frs2, hb_all, cif, charge, solvent=None, title=''):
    root = pathlib.Path(__file__).resolve().parent.parent
    mol = root / 'molecules' / mid
    (mol / 'input').mkdir(parents=True, exist_ok=True)
    (mol / 'cases').mkdir(parents=True, exist_ok=True)
    xyz = mol / 'input' / f'{mid}_cluster.xyz'
    with open(xyz, 'w') as f:
        f.write(f'{len(els)}\n{title or mid} hydrogen-bonded cluster built from {os.path.basename(cif)} by dft/build_cluster.py\n')
        for e, (x, y, z) in zip(els, X):
            f.write(f'{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n')
    rs = rings(els, X)
    planes = {}
    for k, r in enumerate(sorted(rs, key=len, reverse=True)):
        comp = ''.join(sorted({els[i] for i in r}))
        planes[f'ring{k+1}_{comp}'] = [int(r[0]), int(r[len(r)//3]), int(r[2*len(r)//3])]
    frag_named = {}
    for k, fr in enumerate(frs2):
        comp = {}
        for i in fr: comp[els[i]] = comp.get(els[i], 0) + 1
        frag_named[''.join(f'{e}{n}' for e, n in sorted(comp.items()))] = [int(i) for i in fr]
    strongest = min(hb_all, key=lambda h: h['H_A']) if hb_all else None
    case = {
        'title': title or mid,
        'xyz': str(xyz.relative_to(root)),
        'charge': charge, 'spin': 0,
        'source_cif': str((mol / 'input' / os.path.basename(cif)).relative_to(root)),
        'cif_labels': labs,
        'fragments': frag_named,
        'hbonds': {f"{labs[h['D']]}-{labs[h['H']]}...{labs[h['A']]}": [int(h['D']), int(h['H']), int(h['A'])]
                   for h in sorted(hb_all, key=lambda h: h['H_A'])},
        'distances': {f"{labs[h['D']]}...{labs[h['A']]} (X-ray {h['D_A']:.3f} A)": [int(h['D']), int(h['A'])]
                      for h in sorted(hb_all, key=lambda h: h['H_A'])[:4]},
        'planes': planes,
        'cell': cell,
    }
    if strongest:
        case['proton'] = {'H': int(strongest['H']), 'donor_N': int(strongest['D']), 'acceptor_O': int(strongest['A'])}
    if solvent: case['solvent'] = solvent
    out = mol / 'cases' / f'{mid}_gas.json'
    json.dump(case, open(out, 'w'), indent=2)
    print(f'[cluster] wrote {xyz.relative_to(root)} and {out.relative_to(root)}')
    print(f'[cluster] rings for ELF/LOL planes: {list(planes)}')
    return case

if __name__ == '__main__':
    cif = sys.argv[1]; mid = sys.argv[2]
    core = sys.argv[sys.argv.index('--core') + 1] if '--core' in sys.argv else None
    chg = int(sys.argv[sys.argv.index('--charge') + 1]) if '--charge' in sys.argv else 0
    title = sys.argv[sys.argv.index('--title') + 1] if '--title' in sys.argv else ''
    # --require-hbond D:A[,D:A...]  force the assembly to contain these hydrogen bonds
    req = ([tuple(pair.split(':')) for pair in sys.argv[sys.argv.index('--require-hbond') + 1].split(',')]
           if '--require-hbond' in sys.argv else None)
    res = build(cif, mid, core, chg, require=req)
    write_outputs(mid, res[0], res[1], res[2], res[3], res[4], res[5], cif, chg, title=title)
