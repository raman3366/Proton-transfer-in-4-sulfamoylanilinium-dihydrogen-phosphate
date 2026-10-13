"""d03: X-ray versus optimised geometry, by bond class, and the proton metric in every environment.

Bonds are taken from the X-ray ion pair (structures/sfa_pair_xray.xyz, CIF coordinates, anion at x, y, z-1) with
a covalent-radius criterion, classified, and measured in each model with our own code.  The periodic model is
read from periodic_qe/relax_*.extxyz (formula unit 1, minimum-image bonds).  delta = d(N1-H1A) - d(H1A...O5).
"""
import itertools
from common import *

RCOV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05}
sx, X, _ = read_xyz(REPO / 'structures/sfa_pair_xray.xyz')
LAB = PAIR_LABELS
bonds = []
for i, j in itertools.combinations(range(len(LAB)), 2):
    if dist(X[i], X[j]) < 1.2 * (RCOV[sx[i]] + RCOV[sx[j]]):
        bonds.append((i, j))


def bclass(i, j):
    a, b = sorted([LAB[i], LAB[j]])
    ea, eb = sorted([sx[i], sx[j]])
    names = {LAB[i], LAB[j]}
    if 'H' in (sx[i], sx[j]):
        other = sx[j] if sx[i] == 'H' else sx[i]
        return f'{other}-H'
    if ea == eb == 'C':
        return 'C-C (ring)'
    if names == {'C1', 'N1'}:
        return 'C1-N1'
    if names == {'C4', 'S1'}:
        return 'C4-S1'
    if names == {'S1', 'N2'}:
        return 'S1-N2'
    if 'S1' in names:
        return 'S1-O'
    if names & {'O3', 'O4'}:
        return 'P1-O(H) (O3, O4)'
    return 'P1-O (O5, O6)'


CLASSES = ['C-C (ring)', 'C1-N1', 'C4-S1', 'S1-O', 'S1-N2', 'P1-O(H) (O3, O4)', 'P1-O (O5, O6)', 'N-H', 'C-H', 'O-H']
bond_list = [dict(i=i, j=j, name=f'{LAB[i]}-{LAB[j]}', cls=bclass(i, j), xray=dist(X[i], X[j])) for i, j in bonds]

# heavy-atom bond angles (both bonds between heavy atoms)
heavy = [k for k, s in enumerate(sx) if s != 'H']
nb = {k: set() for k in range(len(LAB))}
for i, j in bonds:
    nb[i].add(j); nb[j].add(i)
angles = []
for c in heavy:
    hn = sorted(n for n in nb[c] if sx[n] != 'H')
    for a, b in itertools.combinations(hn, 2):
        angles.append((a, c, b))


def stats(dev):
    dev = np.asarray(dev, float)
    if len(dev) == 0:
        return None
    k = int(np.argmax(np.abs(dev)))
    return dict(n=len(dev), mae=float(np.mean(np.abs(dev))), rmsd=float(np.sqrt(np.mean(dev ** 2))),
                mean_signed=float(np.mean(dev)), max_abs=float(np.abs(dev).max()), argmax=k)


def load_periodic(name, fu=1):
    syms, xyz, labels, cell, head = read_extxyz(REPO / f'molecules/sfa/results/periodic_qe/{name}.extxyz')
    lab2 = {}
    for s, p, l in zip(syms, xyz, labels):
        lab2[l] = p
    inv = np.linalg.inv(cell.T)

    def mic(a, b):             # vector b - a, minimum image
        d = b - a
        f = inv @ d
        f -= np.round(f)
        return cell.T @ f
    return lab2, mic, cell, labels


def model_geom(name):
    """Return coordinates in PAIR_LABELS order as a dict label -> xyz, plus a distance function."""
    if name.startswith('periodic:'):
        lab2, mic, cell, labels = load_periodic(name.split(':')[1])
        P = {l: lab2[l] for l in LAB}                      # formula unit 1 labels are the bare CIF labels
        # put every atom of the ion pair next to N1 (minimum image), i.e. rebuild the pair around N1
        base = P['N1']
        Q = {l: base + mic(base, P[l]) for l in LAB}
        # anion: bring it next to O5 of the pair (the O5 nearest to H1A)
        return Q, lab2, mic
    s, xyz, _ = read_xyz(name)
    return {l: xyz[i] for i, l in enumerate(LAB)}, None, None


def pmetrics(P, mic=None, lab2=None):
    d = lambda a, b: float(np.linalg.norm(mic(P[a], P[b]))) if mic else dist(P[a], P[b])
    m = dict(d_N1_H1A=d('N1', 'H1A'), d_H1A_O5=d('H1A', 'O5'), d_N1_O5=d('N1', 'O5'))
    m['delta'] = m['d_N1_H1A'] - m['d_H1A_O5']
    for o in ['O3', 'O4', 'O5', 'O6']:
        m[f'P1-{o}'] = d('P1', o)
    m['C1-N1'] = d('C1', 'N1')
    m['d_C6_O4'] = d('C6', 'O4'); m['d_H6_O4'] = d('H6', 'O4')
    if mic is None:
        m['ang_N1_H1A_O5'] = angle(P['N1'], P['H1A'], P['O5'])
        m['ang_C6_H6_O4'] = angle(P['C6'], P['H6'], P['O4'])
        m['tau_C3_C4_S1_N2'] = dihedral(P['C3'], P['C4'], P['S1'], P['N2'])
        m['tau_C5_C4_S1_N2'] = dihedral(P['C5'], P['C4'], P['S1'], P['N2'])
    if mic is not None:
        # torsions/angles through bonded chains, each step by minimum image
        c3 = P['C3']; c4 = c3 + mic(P['C3'], P['C4']); s1 = c4 + mic(P['C4'], P['S1']); n2 = s1 + mic(P['S1'], P['N2'])
        c5 = c4 + mic(P['C4'], P['C5'])
        m['tau_C3_C4_S1_N2'] = dihedral(c3, c4, s1, n2)
        m['tau_C5_C4_S1_N2'] = dihedral(c5, c4, s1, n2)
    m['proton_on'] = 'N' if m['d_N1_H1A'] < m['d_H1A_O5'] else 'O'
    return m


MODELS = {
    'xray': REPO / 'structures/sfa_pair_xray.xyz',
    'pcm': RES / 'sfa_pair_pcm_d3/opt.xyz',
    'gas': RES / 'sfa_pair_gas_d3/opt.xyz',
    'pcm_b3lyp': RES / 'sfa_pair_pcm_b3lyp/opt.xyz',
    'pcm_tz': RES / 'sfa_pair_pcm_d3_tz/opt.xyz',
    'periodic_i': 'periodic:relax_i',
    'periodic_ii': 'periodic:relax_ii',
    'periodic_iifix': 'periodic:relax_iifix',
    'periodic_iirfix': 'periodic:relax_iirfix',          # constrained neutral (H1A and O5 held), burst 2
    'periodic_start_i': 'periodic:start_i_xray',
}
out = dict(bonds=[b['name'] for b in bond_list], bond_classes={}, models={})
for b in bond_list:
    out['bond_classes'].setdefault(b['cls'], []).append(b['name'])
for key, path in MODELS.items():
    P, lab2, mic = model_geom(str(path))
    rec = dict(source=rel(path) if not str(path).startswith('periodic:') else
               f"molecules/sfa/results/periodic_qe/{str(path).split(':')[1]}.extxyz")
    bl = {}
    for b in bond_list:
        a, c = LAB[b['i']], LAB[b['j']]
        v = float(np.linalg.norm(mic(P[a], P[c]))) if mic else dist(P[a], P[c])
        bl[b['name']] = v
    rec['bond_lengths'] = bl
    excl = {'N1-H1A'} if key == 'gas' else set()      # H1A has moved to O5 in the gas phase
    dev_by = {}
    for b in bond_list:
        if b['name'] in excl:
            continue
        dev_by.setdefault(b['cls'], []).append((b['name'], bl[b['name']] - b['xray']))
    cls_stats = {}
    for cname, lst in dev_by.items():
        s = stats([x[1] for x in lst])
        s['worst'] = lst[s.pop('argmax')][0]
        cls_stats[cname] = s
    heavy_devs = [(b['name'], bl[b['name']] - b['xray']) for b in bond_list if 'H' not in (sx[b['i']], sx[b['j']])]
    xh_devs = [(b['name'], bl[b['name']] - b['xray']) for b in bond_list
               if 'H' in (sx[b['i']], sx[b['j']]) and b['name'] not in excl]
    all_devs = [(b['name'], bl[b['name']] - b['xray']) for b in bond_list if b['name'] not in excl]
    for nm, lst in [('heavy', heavy_devs), ('XH', xh_devs), ('all', all_devs)]:
        s = stats([x[1] for x in lst])
        s['worst'] = lst[s.pop('argmax')][0]
        cls_stats[f'_{nm}'] = s
    rec['bond_stats'] = cls_stats
    # heavy-atom angles
    if mic is None:
        ad = []
        for a, c, b in angles:
            ad.append((f'{LAB[a]}-{LAB[c]}-{LAB[b]}',
                       angle(P[LAB[a]], P[LAB[c]], P[LAB[b]]) - angle(X[a], X[c], X[b])))
    else:
        ad = []
        for a, c, b in angles:
            va, vb = mic(P[LAB[c]], P[LAB[a]]), mic(P[LAB[c]], P[LAB[b]])
            ang = float(np.degrees(np.arccos(va @ vb / np.linalg.norm(va) / np.linalg.norm(vb))))
            ad.append((f'{LAB[a]}-{LAB[c]}-{LAB[b]}', ang - angle(X[a], X[c], X[b])))
    s = stats([x[1] for x in ad]); s['worst'] = ad[s.pop('argmax')][0]
    rec['angle_stats_heavy'] = s
    rec['angle_devs'] = dict(ad)
    # heavy-atom superposition RMSD of the ion pair and of each ion
    if mic is None:
        hv = [l for l, s_ in zip(LAB, sx) if s_ != 'H']
        rec['kabsch_rmsd_heavy_pair'] = kabsch_rmsd([P[l] for l in hv], [X[LAB.index(l)] for l in hv])
        hc = [l for l in hv if l in CATION]
        rec['kabsch_rmsd_heavy_cation'] = kabsch_rmsd([P[l] for l in hc], [X[LAB.index(l)] for l in hc])
        ha = [l for l in hv if l in ANION]
        rec['kabsch_rmsd_heavy_anion'] = kabsch_rmsd([P[l] for l in ha], [X[LAB.index(l)] for l in ha])
    rec['proton'] = pmetrics(P, mic, lab2)
    out['models'][key] = rec

# proton metric in the other environments (eps ladder, scan anchors)
extra = {
    'eps_4': RES / 'sfa_pair_pcm_d3/eps_ladder/eps_4_N/opt.xyz',
    'eps_8': RES / 'sfa_pair_pcm_d3/eps_ladder/eps_8_N/opt.xyz',
    'eps_20': RES / 'sfa_pair_pcm_d3/eps_ladder/eps_20_N/opt.xyz',
    'scan_A_anchor': RES / 'sfa_pair_pcm_d3/scan/A_00_anchor.xyz',
    'scan_B_anchor': RES / 'sfa_pair_pcm_d3/scan/B_00_anchor.xyz',
}
out['proton_extra'] = {}
for k, p in extra.items():
    s, xyz, com = read_xyz(p)
    P = {l: xyz[i] for i, l in enumerate(LAB)}
    out['proton_extra'][k] = dict(source=rel(p), comment=com, **pmetrics(P))

# fragments vs the same atoms in the X-ray pair (cation, sulfanilamide = cation minus H1A, H2PO4-)
frag = {}
for c, labs in [('cation_pcm_d3', CATION), ('cation_gas_d3', CATION),
                ('sulfanilamide_pcm_d3', [l for l in CATION if l != 'H1A']),
                ('sulfanilamide_gas_d3', [l for l in CATION if l != 'H1A']),
                ('sulfanilamide_pcm_d3_tz', [l for l in CATION if l != 'H1A']),
                ('sulfonamidate_pcm_d3', [l for l in CATION if l not in ('H1A', 'H2A')]),
                ('h2po4_pcm_d3', ANION), ('h2po4_gas_d3', ANION)]:
    s, xyz, _ = read_xyz(RES / c / 'opt.xyz')
    P = {l: xyz[i] for i, l in enumerate(labs)}
    devs = []
    for b in bond_list:
        a, d_ = LAB[b['i']], LAB[b['j']]
        if a in P and d_ in P and 'H' not in (sx[b['i']], sx[b['j']]):
            devs.append((b['name'], dist(P[a], P[d_]) - b['xray']))
    st = stats([x[1] for x in devs]); st['worst'] = devs[st.pop('argmax')][0]
    r = dict(heavy_bond_stats=st, bond_lengths={n: v + [b['xray'] for b in bond_list if b['name'] == n][0]
                                                for n, v in devs})
    if 'C3' in P and 'N2' in P:
        r['tau_C3_C4_S1_N2'] = dihedral(P['C3'], P['C4'], P['S1'], P['N2'])
        r['tau_C5_C4_S1_N2'] = dihedral(P['C5'], P['C4'], P['S1'], P['N2'])
        r['tau_fold'] = max(abs(r['tau_C3_C4_S1_N2']), abs(r['tau_C5_C4_S1_N2']))
        r['ang_C4_S1_N2'] = angle(P['C4'], P['S1'], P['N2'])
        r['d_S1_N2'] = dist(P['S1'], P['N2'])
        r['d_C1_N1'] = dist(P['C1'], P['N1'])
    hv = [l for l in labs if not l.startswith('H')]
    r['kabsch_rmsd_heavy'] = kabsch_rmsd([P[l] for l in hv], [X[LAB.index(l)] for l in hv])
    frag[c] = r
out['fragments'] = frag
# torsion also for the pair models
dump('d03_geometry.json', out)
