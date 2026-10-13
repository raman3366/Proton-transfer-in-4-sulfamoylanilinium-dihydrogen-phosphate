"""
M6 (referee chem, M6 item 1): QTAIM of all nine crystal hydrogen bonds, in the salt and in the constrained
neutral cell, on a charge-neutral cluster cut from the periodic PBE-D3(BJ) geometries.

Why a cluster and not the periodic density: the QE runs use pseudopotentials (SSSP US/PAW), so their valence density
has no core and is augmented; critic2 on it would need a core-augmented PAW reconstruction that was never set up or
tested in this project.  A cluster cut from the relaxed cell, with every atom kept where the periodic relaxation put
it, gives an all-electron density that Multiwfn's QTAIM (dft/mwfn_extra.py, validated on the VM in burst 1c) reads.
The referee named this route ("or a cluster cut from the periodic geometry").

The cluster (stage 'build')
  * Geometries: relax_i (salt, proton on N1) and relax_iirfix (constrained neutral: H1A held 0.98 A from O5, O5 held,
    everything else relaxed), molecules/sfa/results/periodic_qe/*.extxyz, fixed X-ray cell.  Nothing is moved:
    molecules are only made whole by lattice translations.  Checked at build time: every cluster atom is its periodic
    atom plus an integer lattice vector (to 1e-15), no covalent bond of the cell is cut, and every H-bond has the
    D-H / H...A / D...A of the periodic report (periodic_qe.json, formula unit 1) to < 1e-4 A.
  * Building block = the M1 ion pair: a cation plus the anion that receives its N1-H1A...O5 bond.  Whole pairs make
    the cluster neutral in BOTH structures (in relax_iirfix every H1A sits on its own pair's O5, so each pair is
    sulfanilamide + H3PO4), give equal numbers of cations and anions and the same atom list, so the salt and neutral
    clusters differ only in coordinates.  (A cation without its O5 partner would leave H1A dangling 0.98 A from a
    missing oxygen in relax_iirfix, and an extra anion would break neutrality: whole pairs are the only consistent cut.)
    The pair of the M1 model (asymmetric-unit cation + anion 1_554) is the central pair P0 (molecules C0 and A0).
  * Selection, exhaustive over the 57 pairs within 9 A of P0 (types found inside every pair and between every two
    pairs tabulated once, combinations scored by set unions):
      --cover minimal (default): the fewest pairs, P0 included, in which each of the nine CIF bonds occurs at least
        once between complete molecules.  Result: 5 pairs, 135 atoms (no 4-pair cover exists).  Ties (704 covers of
        size 5) are broken by the number of bond types that occur AT the central pair (8/9: only N1-H1B...O6 is
        between two outer molecules, C1 > A4), then by the number of H-bond instances (18), then by compactness
        (Rg 6.64 A); every pair is H-bond connected to the rest.
      --cover central: all nine bonds at the central pair: 6 pairs, 162 atoms, 2106 basis functions; its exchange
        build costs ~2.1x per SCF cycle and its DF tensor (182 GB) twice the disk -- too slow for burst 3b (Cost).
  * The same pair images are used for every structure: atoms of the other structures are put at the lattice image
    nearest to the corresponding salt atom (relax_iirfix moves no atom by more than 0.7 A).
  * Every instance of the nine bond types inside the cluster is written to the case file ('hbonds'), so Multiwfn's
    BCPs are named per instance; the 'primary' instance of each type is the one at the central pair (C0/A0 as
    donor first).  Where a type occurs more than once (N1-H1A...O5 and C6-H6...O4 x5, N1-H1B...O6 x2) the spread
    between instances is a direct measure of the cluster-edge effect (compare stage, rho_instances_*).

Level of theory (stage 'scf')
  B3LYP-D3(BJ)/6-31+G(d,p) (spherical 5d; the basis the referee suggested for clusters, M1 item 1), density fitting
  with PySCF's default auxiliary basis for it (heavy-aug-cc-pVDZ-JKFIT), grid level 3 (GRID_LEVEL; ~Gaussian
  FineGrid; on the central ion pair, grid 4 changes rho at the H1A...O5 and H6...O4 midpoints by < 1e-6 a.u., Pi
  test 29 Sep 2026), gas phase, single point at the periodic geometry, no optimisation, conv_tol 1e-8 Ha.  --eps E adds
  IEF-PCM with dielectric constant E (Lebedev 17 cavity, PCM_LEBEDEV_ORDER).  The D3(BJ) term does not change the
  density; it is included so that the cluster energies are B3LYP-D3(BJ) energies.  Initial guess: the block-diagonal
  sum of the converged densities of the separate ion pairs (same functional, basis and grid; a pair that is a lattice
  translation of an earlier one reuses its density, so 2 pair SCFs per structure), MINAO if that fails;
  convergence ladder DIIS (80 cycles) -> level shift 0.25 Ha (60) -> Newton (20 macro-iterations), else the task fails.
  Outputs results/<case>/final.chk, final.molden (padded to a square MO matrix as dft/fmo.py does), scf.json.

QTAIM (stage 'qtaim'): runs 'dft/mwfn_extra.py <case> qtaim' unchanged (the Multiwfn menu sequence every other case of
  the project used: CP searches from nuclei, atom pairs, triangles, pyramids and nuclear spheres; bond paths; all CP
  properties), then picks for every case H-bond instance the BCP that dft/mwfn_extra.py matched to it (the H...A BCP,
  or N1...H1A when the proton sits on O5), plus the H...A and D-H BCPs by atom pair: rho, Laplacian, G, V, H = G + V,
  -G/V, ellipticity, eta.  -> results/<case>/nine_bcps.json.  MWFN_THREADS sets Multiwfn's threads.

Comparison (stage 'compare'): per bond type, primary instance: salt vs constrained neutral (and any extra structure,
  e.g. an M2 profile point built with --extra): values, neutral - salt changes (rho also in %), the instance spread;
  plus the cluster energy difference per formula unit (B3LYP-D3(BJ) at the PBE geometries, all pairs switched; a
  gas-phase cluster with under-coordinated surface ions, so qualitative only next to the periodic 57.1/52.8 kJ/mol).
  -> molecules/sfa/results/crystal_cluster_qtaim/qtaim_nine_bonds{,_epsE}.{json,csv,md}

Cost (measured 29 Sep 2026 on the Pi, 2 threads, with the Pi otherwise loaded; VM figures from data/benchmark)
  * Size (5 pairs): 135 atoms, 700 electrons, 1755 basis functions, 8565 auxiliary; DF tensor 106 GB, so on disk
    (/mnt/scratch, read back with prefetch every cycle) unless DFT_MAX_MEMORY_MB > ~118000.
  * Exchange: 4 nocc nao^2 naux = 3.7e13 flop per SCF cycle.  PySCF's DF exchange ran at 74% of the dgemm rate on
    the Pi (1 pair: 9.0 of 12.1 GFLOP/s), so it is dgemm-bound.  XC (grid 3) grew 5.9x from 1 to 2 pairs (30.7 ->
    181.6 s on the Pi).  Extrapolated to 5 pairs on the Pi: ~4000 s exchange + ~2000 s XC per cycle.
  * Pi -> VM: the 445-function benchmark SCF (B3LYP/6-311++G(d,p), grid 4) took 1373 s on the Pi (2 threads) vs 44 s
    on the VM (12 threads, 31x); the dgemm-bound exchange is taken at ~45x for 16 VM threads.  Hence ~150 s per
    cycle (range 100-220 s on a VM shared with other burst tasks), 16-24 cycles (the Pi test clusters took 9 and 26
    from the pair guess), + ~5 min (pair guess, DF tensor): ~55 min per structure at 16 threads (range 30-90 min).
    Burst est_min 75 (x1.3 over ~55); hard cap 120 min.  The 6-pair 'central' cluster: ~2x per cycle.
  * Multiwfn QTAIM: 27 atoms took 4 s; 135 atoms with ~7x the starting points and ~5x the occupied orbitals per
    density evaluation: ~5 min at 8 threads (est_min 12).  IEF-PCM (Lebedev 17): estimated, not measured, at ~15%
    more per cycle, but it is a second full SCF per structure (optional tasks).

Usage
  python dft/crystal_cluster_qtaim.py build [--cover minimal|central] [--extra TAG=path.extxyz ...] [--only TAG ...]
  python dft/crystal_cluster_qtaim.py estimate                        # basis functions, memory, exchange time
  python dft/crystal_cluster_qtaim.py scf salt|neutral|TAG [--eps E] [--conv-tol T] [--no-frag-guess]
  python dft/crystal_cluster_qtaim.py qtaim salt|neutral|TAG [--eps E]
  python dft/crystal_cluster_qtaim.py compare [--eps E]
  python dft/crystal_cluster_qtaim.py test [--keep]                   # Pi self-test (see stage_test)
Environment: DFT_THREADS (PySCF/BLAS threads, pinned here before numpy loads), DFT_MAX_MEMORY_MB (PySCF memory),
  PYSCF_TMPDIR (default: a per-process directory on /mnt/scratch, removed on exit), GRID_LEVEL (default 3 here),
  PCM_LEBEDEV_ORDER (default 17 here), MWFN_THREADS, MWFN_BIN.
"""
import os, sys, atexit, shutil, signal, tempfile, pathlib


def _pin_environment():
    """Same pinning as dft/proton_tools.py (must run before numpy/PySCF load): every BLAS/OpenMP pool to
    DFT_THREADS, PySCF scratch in a per-process directory on the NVMe volume, removed on exit."""
    t = os.environ.get('DFT_THREADS')
    if t:
        for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                  'VECLIB_MAXIMUM_THREADS'):
            os.environ[v] = t
    os.environ.setdefault('DFT_MAX_MEMORY_MB', '6000')
    if 'PYSCF_TMPDIR' not in os.environ and os.path.isdir('/mnt/scratch') and os.access('/mnt/scratch', os.W_OK):
        tmp = tempfile.mkdtemp(prefix='xclu_', dir='/mnt/scratch')
        os.environ['PYSCF_TMPDIR'] = os.environ['TMPDIR'] = tmp
        atexit.register(shutil.rmtree, tmp, True)
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))


_pin_environment()

import re, json, time, hashlib, itertools, argparse, subprocess, datetime   # noqa: E402
import numpy as np                                                          # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
MOL_DIR = ROOT / 'molecules' / 'sfa'
PQE = MOL_DIR / 'results' / 'periodic_qe'
OUT = MOL_DIR / 'results' / 'crystal_cluster_qtaim'
STRUCTS = {'salt': PQE / 'relax_i.extxyz', 'neutral': PQE / 'relax_iirfix.extxyz'}
STRUCT_TITLES = {'salt': 'periodic salt relax_i (proton on N1)',
                 'neutral': 'periodic constrained neutral relax_iirfix (H1A held 0.98 A from O5)'}
CASE_PREFIX = 'sfa_xclu'
CLUSTER_BASIS = '6-31+g(d,p)'
BASIS_PRINT = {'6-31+g(d,p)': '6-31+G(d,p)', 'sto-3g': 'STO-3G', '3-21g': '3-21G', '6-31g': '6-31G'}

# The nine hydrogen bonds of the deposited _geom_hbond loop (Table 1 of the 2014 paper): D, H, A, CIF code of A,
# X-ray D...A (A).  Bonds are identified here by atom labels and geometry, not by the codes (the periodic cell
# wraps atoms differently from the CIF); the codes and X-ray distances are carried for the report only.
CIF_HBONDS = [('N1', 'H1A', 'O5', '1_554', 2.760), ('N1', 'H1B', 'O6', '2_554', 2.856), ('N1', 'H1C', 'O6', '2_654', 2.855),
              ('N2', 'H2A', 'O3', '2_564', 3.076), ('N2', 'H2B', 'O2', '1_655', 2.886), ('O3', 'H3A', 'O1', '1_455', 2.750),
              ('O4', 'H4A', 'O5', '1_455', 2.500), ('C3', 'H3', 'O6', '1_655', 3.283), ('C6', 'H6', 'O4', '1_554', 3.288)]
TYPE_NAMES = [f'{d}-{h}...{a}' for d, h, a, _, _ in CIF_HBONDS]
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05}
BOND_SCALE = 1.25          # covalent bond if d < 1.25 (r_a + r_b): N1...H1A 1.616 A in relax_iirfix is NOT bonded
HB_MAX_HA, HB_MAX_DA, HB_MIN_ANG = 2.9, 3.6, 110.0   # geometric search window for the nine bonds (all 9 are far inside)
POOL_RANGE = 2             # translations -2..2 along a, b, c for the pool of candidate pairs
POOL_CUTOFF = 9.0          # a candidate pair must have an atom within this distance (A) of the central pair


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')


def save_json(path, obj):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def base_label(lab):
    """'O5_2' (second formula unit of the periodic cell) -> 'O5'."""
    return re.sub(r'_2$', '', lab)


# ------------------------------------------------------------------------------------------------ structures ----
def read_extxyz(path):
    lines = pathlib.Path(path).read_text().splitlines()
    n = int(lines[0])
    lat = np.array([float(v) for v in re.search(r'Lattice="([^"]+)"', lines[1]).group(1).split()]).reshape(3, 3)
    syms, X, labs = [], [], []
    for l in lines[2:2 + n]:
        p = l.split()
        syms.append(p[0]); X.append([float(v) for v in p[1:4]]); labs.append(p[4])
    return syms, np.array(X), labs, lat


SHIFTS = np.array(list(itertools.product((-1, 0, 1), repeat=3)))


def min_image(d, lat):
    """Minimum-image vector of d (small monoclinic cell: brute force over the 27 neighbouring images)."""
    cand = d[None, :] + SHIFTS @ lat
    return cand[np.argmin(np.linalg.norm(cand, axis=1))]


def molecules_in_cell(syms, X, lat):
    """Covalently bonded molecules of the periodic cell, each unwrapped to a contiguous set of positions."""
    n = len(syms)
    nb = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            v = min_image(X[j] - X[i], lat)
            if np.linalg.norm(v) < BOND_SCALE * (COV[syms[i]] + COV[syms[j]]):
                nb[i].append((j, v)); nb[j].append((i, -v))
    seen, mols = {}, []
    for s in range(n):
        if s in seen:
            continue
        pos = {s: X[s].copy()}; stack = [s]; seen[s] = len(mols)
        while stack:
            i = stack.pop()
            for j, v in nb[i]:
                if j not in pos:
                    pos[j] = pos[i] + v; seen[j] = len(mols); stack.append(j)
        mols.append(pos)
    return mols


def cell_molecules(path):
    """The four molecules of the Z=2 cell (cation/anion of formula units 1 and 2), unwrapped, keyed 'cat1' ..."""
    syms, X, labs, lat = read_extxyz(path)
    out = {}
    for pos in molecules_in_cell(syms, X, lat):
        idx = sorted(pos)
        el = [syms[i] for i in idx]
        fu = 2 if labs[idx[0]].endswith('_2') else 1
        kind = 'cat' if 'N' in el and 'S' in el else 'an' if 'P' in el else None
        if kind is None or any(labs[i].endswith('_2') != (fu == 2) for i in idx):
            raise SystemExit(f'[xclu] unexpected molecule in {path}: {[labs[i] for i in idx]}')
        out[f'{kind}{fu}'] = {'idx': idx, 'labels': [labs[i] for i in idx], 'syms': el,
                              'X': np.array([pos[i] for i in idx])}
    if sorted(out) != ['an1', 'an2', 'cat1', 'cat2']:
        raise SystemExit(f'[xclu] {path}: expected 2 cations + 2 anions, found {sorted(out)} (proton transferred? '
                         'molecules are defined on the salt structure)')
    return {'syms': syms, 'X': X, 'labels': labs, 'lat': lat, 'mols': out}


def map_structure(salt, path):
    """Positions of another structure of the same cell, per salt molecule, each atom at the lattice image nearest
    to the unwrapped salt atom (molecule membership is by label, i.e. H1A stays with 'its' cation index-wise)."""
    syms, X, labs, lat = read_extxyz(path)
    if labs != salt['labels'] or syms != salt['syms']:
        raise SystemExit(f'[xclu] {path}: atom labels/order differ from the salt structure')
    if np.abs(lat - salt['lat']).max() > 1e-6:
        raise SystemExit(f'[xclu] {path}: lattice differs from the salt structure (the cluster assumes the fixed X-ray cell)')
    out, maxdev = {}, 0.0
    for k, m in salt['mols'].items():
        P = []
        for i, xs in zip(m['idx'], m['X']):
            v = min_image(X[i] - xs, lat)
            P.append(xs + v)
            maxdev = max(maxdev, float(np.linalg.norm(v)))
        out[k] = dict(m, X=np.array(P))
    return {'syms': syms, 'X': X, 'labels': labs, 'lat': lat, 'mols': out, 'max_atom_shift_vs_salt_A': maxdev}


def pair_partner(cell):
    """For each formula unit's cation: (anion key, lattice translation) of the anion whose O5 takes N1-H1A."""
    lat, res = cell['lat'], {}
    for fu in (1, 2):
        c = cell['mols'][f'cat{fu}']
        h = c['X'][c['labels'].index('H1A' if fu == 1 else 'H1A_2')]
        best = None
        for ak in ('an1', 'an2'):
            a = cell['mols'][ak]
            o = a['X'][[base_label(l) for l in a['labels']].index('O5')]
            for n in itertools.product(range(-2, 3), repeat=3):
                d = np.linalg.norm(o + np.array(n) @ lat - h)
                if best is None or d < best[0]:
                    best = (d, ak, n)
        if best[0] > 2.2:
            raise SystemExit(f'[xclu] no O5 within 2.2 A of H1A of cation {fu}')
        res[fu] = (best[1], tuple(int(v) for v in best[2]))
    return res


class Pool:
    """Candidate pairs (fu, n): cation of formula unit fu translated by n, plus its N1-H1A...O5 anion."""

    def __init__(self, salt):
        self.cell, self.lat = salt, salt['lat']
        self.partner = pair_partner(salt)
        pairs = [(fu, n) for fu in (1, 2) for n in itertools.product(range(-POOL_RANGE, POOL_RANGE + 1), repeat=3)]
        c0 = self.pair_xyz(salt, (1, (0, 0, 0)))
        keep = []
        for p in pairs:
            X = self.pair_xyz(salt, p)
            dmin = np.min(np.linalg.norm(X[:, None, :] - c0[None, :, :], axis=-1))
            if dmin < POOL_CUTOFF:
                keep.append(p)
        self.pairs = keep
        self.center = c0.mean(0)

    def mol_xyz(self, cell, key, n):
        return cell['mols'][key]['X'] + np.array(n) @ self.lat

    def pair_members(self, p):
        fu, n = p
        ak, t = self.partner[fu]
        return [(f'cat{fu}', tuple(n)), (ak, tuple(np.add(n, t)))]

    def pair_xyz(self, cell, p):
        return np.vstack([self.mol_xyz(cell, k, n) for k, n in self.pair_members(p)])


def find_instances(pool, cell, mol_list):
    """All occurrences of the nine CIF bonds between the molecules in mol_list [(key, n), ...] in structure `cell`,
    judged on the SALT geometry (the bond topology of the CIF): D-H in one molecule, A in another."""
    res = []
    info = []
    for key, n in mol_list:
        m = cell['mols'][key]
        info.append({'bl': [base_label(l) for l in m['labels']], 'X': pool.mol_xyz(cell, key, n), 'key': key, 'n': n})
    for t, (D, H, A, _, _) in enumerate(CIF_HBONDS):
        for i, di in enumerate(info):
            if D not in di['bl'] or H not in di['bl']:
                continue
            xd, xh = di['X'][di['bl'].index(D)], di['X'][di['bl'].index(H)]
            for j, aj in enumerate(info):
                if i == j or A not in aj['bl']:
                    continue
                xa = aj['X'][aj['bl'].index(A)]
                ha, da = np.linalg.norm(xa - xh), np.linalg.norm(xa - xd)
                v1, v2 = xd - xh, xa - xh
                ang = np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1)))
                if ha < HB_MAX_HA and da < HB_MAX_DA and ang > HB_MIN_ANG:
                    res.append({'type': t, 'donor_mol': i, 'acceptor_mol': j, 'H_A': float(ha), 'D_A': float(da),
                                'DHA': float(ang)})
    return res


def select_cluster(pool, cover, kmax=6):
    """Fewest extra pairs covering the nine bond types (see the module docstring), deterministic tie-breaks.
    Exhaustive over the pool: the bond types found inside every pair and between every two pairs are tabulated
    once, so a combination is scored by set unions (the 'any' cover also counts bonds between two non-central
    pairs)."""
    cell = pool.cell
    P0 = (1, (0, 0, 0))
    cands = [p for p in pool.pairs if p != P0]
    members = {p: pool.pair_members(p) for p in pool.pairs}

    def inst(plist):
        mols = [m for p in plist for m in members[p]]
        return mols, find_instances(pool, cell, mols)

    internal = {p: {x['type'] for x in inst([p])[1]} for p in pool.pairs}
    between, at_centre = {}, {}
    for p, q in itertools.combinations(pool.pairs, 2):
        mols, ins = inst([p, q])
        n1 = len(members[p])
        cross = [x for x in ins if (x['donor_mol'] < n1) != (x['acceptor_mol'] < n1)]
        between[frozenset((p, q))] = {x['type'] for x in cross}
        if P0 in (p, q):
            at_centre[q if p == P0 else p] = {x['type'] for x in cross}
    central_types = internal[P0]
    need = set(range(9))
    best, tried, n_opt = None, 0, 0
    for k in range(0, kmax + 1):
        sols = []
        for combo in itertools.combinations(cands, k):
            tried += 1
            if cover == 'central':
                got = central_types.union(*[at_centre[p] for p in combo])
            else:                               # 'minimal'
                got = central_types.union(*[internal[p] | at_centre[p] for p in combo],
                                          *[between[frozenset(pq)] for pq in itertools.combinations(combo, 2)])
            if got == need:
                sols.append(combo)
        if sols:
            scored = []
            for combo in sols:
                plist = [P0, *combo]
                mols, ins = inst(plist)
                X = np.vstack([pool.pair_xyz(cell, p) for p in plist])
                rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
                # H-bond connectivity between pairs
                pair_of = [plist.index(p) for p in plist for _ in members[p]]
                par = list(range(len(plist)))

                def find(i):
                    while par[i] != i:
                        par[i] = par[par[i]]; i = par[i]
                    return i
                for x in ins:
                    par[find(pair_of[x['donor_mol']])] = find(pair_of[x['acceptor_mol']])
                connected = len({find(i) for i in range(len(plist))}) == 1
                # bond types with an instance at the central pair (molecules 0 and 1): the more, the better
                n_central = len({x['type'] for x in ins if {x['donor_mol'], x['acceptor_mol']} & {0, 1}})
                scored.append((not connected, -n_central, -len(ins), rg, combo))
            scored.sort(key=lambda s: (s[0], s[1], s[2], s[3], str(s[4])))
            best = scored[0]
            n_opt = len(sols)
            break
    if best is None:
        raise SystemExit('[xclu] no cover of the nine bonds found with <= 6 extra pairs')
    return {'pairs': [P0, *best[4]], 'n_solutions_same_size': n_opt, 'n_combinations_tried': tried,
            'connected': not best[0], 'n_types_at_central_pair': -best[1], 'n_hbond_instances': -best[2],
            'rg_A': best[3], 'types_inside_central_pair': [TYPE_NAMES[t] for t in sorted(central_types)], 'cover': cover}


def build_cluster(pool, cell, plist):
    """Atoms (symbols, xyz, labels, molecule tags) of the pairs in plist for structure `cell`."""
    syms, X, labs, frags, tags = [], [], [], {}, []
    for k, p in enumerate(plist):
        for (key, n), kind in zip(pool.pair_members(p), ('C', 'A')):
            tag = f'{kind}{k}'
            m = cell['mols'][key]
            i0 = len(syms)
            syms += m['syms']; X.append(pool.mol_xyz(cell, key, n))
            labs += [f'{base_label(l)}@{tag}' for l in m['labels']]
            name = f'{"cation" if kind == "C" else "anion"} {tag}'
            frags[name] = list(range(i0, len(syms)))
            tags.append({'tag': tag, 'molecule': key, 'translation': list(n), 'pair': k, 'fragment': name})
    return syms, np.vstack(X), labs, frags, tags


def periodic_check(cell, pool, plist, X, syms):
    """No geometry change: (1) every cluster atom is its periodic atom plus an integer lattice vector; (2) every
    covalent bond of the cell (minimum image) has the same length in the cluster (molecules are whole).
    Returns (largest non-integer part of a translation in fractional units, largest bond-length deviation in A)."""
    idx = []
    for p in plist:
        for key, n in pool.pair_members(p):
            idx += cell['mols'][key]['idx']
    frac = (X - cell['X'][idx]) @ np.linalg.inv(cell['lat'])
    dfrac = float(np.abs(frac - np.round(frac)).max())
    # covalent bonds of this structure's cell (minimum image; cell edges >= 4.8 A, so a bond is always the minimum
    # image); each must be present in the cluster with the partner exactly at the bonded image
    csyms, cX, lat = cell['syms'], cell['X'], cell['lat']
    nb = {i: [] for i in range(len(csyms))}
    for i in range(len(csyms)):
        for j in range(i + 1, len(csyms)):
            v = min_image(cX[j] - cX[i], lat)
            if np.linalg.norm(v) < BOND_SCALE * (COV[csyms[i]] + COV[csyms[j]]):
                nb[i].append((j, v)); nb[j].append((i, -v))
    where = {}
    for k, i in enumerate(idx):
        where.setdefault(i, []).append(k)
    dangling, worst = 0, 0.0
    for k, i in enumerate(idx):
        for j, v in nb[i]:
            dev = min((np.linalg.norm(X[l] - X[k] - v) for l in where.get(j, [])), default=np.inf)
            if dev > 1e-6:
                dangling += 1
            else:
                worst = max(worst, dev)
    return dfrac, dangling, worst


def electron_count(syms):
    Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16}
    return sum(Z[s] for s in syms)


def write_xyz(path, syms, X, comment):
    with open(path, 'w') as f:
        f.write(f'{len(syms)}\n{comment}\n')
        for s, (x, y, z) in zip(syms, X):
            f.write(f'{s:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n')




def case_name(tag, eps=None, prefix=CASE_PREFIX):
    return f'{prefix}_{tag}' + (f'_eps{eps:g}' if eps else '')


def make_cases(pool, salt, plist, structs, prefix=CASE_PREFIX, cover='central', require_all=True, basis=CLUSTER_BASIS):
    """Cluster xyz + case JSON for every structure; returns the per-structure part of build.json."""
    (MOL_DIR / 'input').mkdir(exist_ok=True)
    ref_geom = json.loads((PQE / 'periodic_qe.json').read_text())['relax']
    out = {}
    for tag, path in structs.items():
        cell = salt if tag == 'salt' else map_structure(salt, path)
        syms, X, labs, frags, tags = build_cluster(pool, cell, plist)
        dfrac, dangling, worst = periodic_check(cell, pool, plist, X, syms)
        if dfrac > 1e-8 or dangling:
            raise SystemExit(f'[xclu] {tag}: geometry check failed: non-lattice translation {dfrac:.1e}, '
                             f'{dangling} covalent bond(s) cut or displaced')
        mols = [m for p in plist for m in pool.pair_members(p)]
        ins = find_instances(pool, salt, mols)      # bond topology from the salt; geometry measured in `cell` below
        hb, rows = {}, []
        for x in ins:
            D, H, A = CIF_HBONDS[x['type']][:3]
            td, ta = tags[x['donor_mol']]['tag'], tags[x['acceptor_mol']]['tag']
            iD, iH, iA = labs.index(f'{D}@{td}'), labs.index(f'{H}@{td}'), labs.index(f'{A}@{ta}')
            name = f'{D}-{H}...{A} [{td}>{ta}]'
            hb[name] = [iD, iH, iA]
            g = hb_geometry(X, iD, iH, iA)
            central = bool({td, ta} & {'C0', 'A0'})
            rows.append({'name': name, 'type': TYPE_NAMES[x['type']], 'cif_code': CIF_HBONDS[x['type']][3],
                         'xray_D_A_A': CIF_HBONDS[x['type']][4], 'donor_mol': td, 'acceptor_mol': ta,
                         'at_central_pair': central, 'atoms': [iD, iH, iA], **{k: round(v, 4) for k, v in g.items()}})
        # primary instance per type: at the central pair (C0/A0), C0/A0 as donor first, then by name
        prim = {}
        for r in rows:
            key = (not r['at_central_pair'], r['donor_mol'] not in ('C0', 'A0'), r['name'])
            if r['type'] not in prim or key < prim[r['type']][0]:
                prim[r['type']] = (key, r['name'])
        missing = [t for t in TYPE_NAMES if t not in prim]
        if missing and require_all:
            raise SystemExit(f'[xclu] {tag}: bond types missing from the cluster: {missing}')
        for r in rows:
            r['primary'] = prim[r['type']][1] == r['name']
        # independent check against the periodic geometry report (formula-unit-1 values in periodic_qe.json)
        pchk = []
        for ref in ref_geom.get({'salt': 'i', 'neutral': 'iirfix'}.get(tag, ''), {}).get('geometry', {}).get('hbonds', []):
            t = ref['bond'].split(' ')[0]
            same = [r for r in rows if r['type'] == t]
            if same:
                dev = min(abs(r['D..A'] - ref['D_A_A']) + abs(r['H..A'] - ref['H_A_A']) + abs(r['D-H'] - ref['D_H_A'])
                          for r in same)
                pchk.append({'type': t, 'periodic_D_A': ref['D_A_A'], 'closest_instance_dev_sum_A': round(dev, 4)})
        n_el = electron_count(syms)
        if n_el % 2:
            raise SystemExit(f'[xclu] {tag}: odd electron count {n_el}')
        name = case_name(tag, None, prefix)
        xyz = MOL_DIR / 'input' / f'{name}.xyz'
        write_xyz(xyz, syms, X, f'SAPH crystal cluster ({len(plist)} ion pairs, cover {cover}) cut from '
                                f'{pathlib.Path(path).name}; no atom moved; dft/crystal_cluster_qtaim.py')
        src = pathlib.Path(path)
        case = {
            'title': f'SAPH crystal cluster, {len(plist)} ion pairs, from the {STRUCT_TITLES.get(tag, tag)}; '
                     f'B3LYP-D3(BJ)/{BASIS_PRINT.get(basis, basis)} single point, gas phase',
            'xyz': str(xyz.relative_to(ROOT)), 'charge': 0, 'spin': 0, 'disp': 'd3bj', 'xc': 'b3lyp', 'basis': basis,
            'structure': tag, 'source_extxyz': str(src.relative_to(ROOT)) if src.is_relative_to(ROOT) else str(src),
            'source_sha256': sha256(src), 'cif_labels': labs, 'fragments': frags,
            'fragment_charges': {k: (1 if k.startswith('cation') else -1) for k in frags},
            'hbonds': hb, 'hbond_types': {r['name']: r['type'] for r in rows},
            'primary_hbonds': {v[1]: t for t, v in prim.items()}, 'molecules': tags,
            'pairs': [[p[0], list(p[1])] for p in plist], 'cover': cover,
            'proton': {'H': labs.index('H1A@C0'), 'donor_N': labs.index('N1@C0'), 'acceptor_O': labs.index('O5@A0')},
            'notes': 'Cut by dft/crystal_cluster_qtaim.py (M6). Whole ion pairs (cation + its N1-H1A...O5 anion), so the '
                     'cluster is neutral with the proton on N1 (salt) or on O5 (constrained neutral). Fragment charges '
                     'are those of the salt; dft/mwfn_extra.py re-derives them from the geometry (H1A on O5 -> 0/0).'}
        save_json(MOL_DIR / 'cases' / f'{name}.json', case)
        out[tag] = {'case': name, 'source': case['source_extxyz'], 'source_sha256': case['source_sha256'],
                    'n_atoms': len(syms), 'n_electrons': n_el, 'n_pairs': len(plist),
                    'geometry_check': {'max_nonlattice_translation_frac': dfrac, 'cut_covalent_bonds': dangling,
                                       'max_bond_length_deviation_A': worst},
                    'max_atom_shift_vs_salt_A': cell.get('max_atom_shift_vs_salt_A', 0.0),
                    'hbond_instances': rows, 'periodic_report_check': pchk, 'types_missing': missing,
                    'n_instances_per_type': {t: sum(r['type'] == t for r in rows) for t in TYPE_NAMES}}
        print(f"[xclu] {tag}: {len(syms)} atoms, {n_el} electrons, {len(rows)} H-bond instances "
              f"({', '.join(f'{t} x{n}' for t, n in out[tag]['n_instances_per_type'].items() if n)}); "
              f"geometry check: lattice translations integer to {dfrac:.0e}, 0 cut bonds; wrote {xyz.relative_to(ROOT)}")
        for r in rows:
            if r['primary']:
                print(f"         primary {r['name']:28s} D-H {r['D-H']:.3f}  H..A {r['H..A']:.3f}  D..A {r['D..A']:.3f}"
                      f"  angle {r['DHA']:.1f}")
        bad = [c for c in pchk if c['closest_instance_dev_sum_A'] > 0.003]
        if bad:
            print(f'[xclu] {tag}: NOTE instances differ from the periodic report (formula unit 1) by > 0.003 A: {bad}')
    return out


def hb_geometry(X, iD, iH, iA):
    v1, v2 = X[iD] - X[iH], X[iA] - X[iH]
    return {'D-H': float(np.linalg.norm(v1)), 'H..A': float(np.linalg.norm(v2)),
            'D..A': float(np.linalg.norm(X[iA] - X[iD])),
            'DHA': float(np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))}


def parse_structs(extra):
    structs = dict(STRUCTS)
    for e in extra or []:
        tag, path = e.split('=', 1)
        if tag in structs or not re.fullmatch(r'[a-z0-9]+', tag):
            raise SystemExit(f'[xclu] bad --extra tag {tag!r} (lower-case letters/digits, not salt/neutral)')
        p = pathlib.Path(path)
        structs[tag] = p if p.is_absolute() else ROOT / p
        if not structs[tag].exists():
            raise SystemExit(f'[xclu] --extra {tag}: {structs[tag]} not found')
    return structs


def stage_build(args):
    salt = cell_molecules(STRUCTS['salt'])
    pool = Pool(salt)
    structs = parse_structs(args.extra)
    if args.only:
        structs = {k: v for k, v in structs.items() if k in args.only}
        if 'salt' not in structs and not (MOL_DIR / 'cases' / f'{case_name("salt")}.json').exists():
            raise SystemExit('[xclu] --only without salt needs an earlier full build')
    print(f"[xclu] pool: {len(pool.pairs)} candidate pairs within {POOL_CUTOFF} A of the central pair; "
          f"pair partner of each cation (anion, translation): {pool.partner}")
    sel = {c: select_cluster(pool, c) for c in ('minimal', 'central')}
    for c, s in sel.items():
        print(f"[xclu] cover '{c}': {len(s['pairs'])} pairs = {27 * len(s['pairs'])} atoms, "
              f"{s['n_types_at_central_pair']}/9 bond types at the central pair, "
              f"{s['n_hbond_instances']} H-bond instances, Rg {s['rg_A']:.2f} A, H-bond connected {s['connected']}, "
              f"{s['n_solutions_same_size']} solution(s) of that size; pairs {s['pairs']}")
    print(f"[xclu] selected cover: {args.cover}")
    report = {'generated': now(), 'script': 'dft/crystal_cluster_qtaim.py build', 'cover': args.cover,
              'selection': {c: dict(s, pairs=[[p[0], list(p[1])] for p in s['pairs']]) for c, s in sel.items()},
              'pair_definition': 'cation of formula unit fu (periodic labels without/with _2) translated by n (lattice '
                                 'vectors), plus the anion that receives its N1-H1A...O5 bond',
              'pair_partner': {str(k): [v[0], list(v[1])] for k, v in pool.partner.items()},
              'structures': make_cases(pool, salt, sel[args.cover]['pairs'], structs, cover=args.cover)}
    report['size'] = size_estimate(next(iter(report['structures'].values()))['case'])
    bj = OUT / ('build.json' if not args.only else f'build_{"_".join(args.only)}.json')
    save_json(bj, report)
    print(f"[xclu] basis functions {report['size']['nao']}, auxiliary {report['size']['naux']}, "
          f"3-centre DF tensor {report['size']['cderi_GB']:.0f} GB; wrote {bj.relative_to(ROOT)}")
    return report


# ------------------------------------------------------------------------------------------------- estimate ----
# Time model for one cluster SCF (docstring 'Cost'): exchange at an effective dgemm-bound rate + XC extrapolated from
# the Pi (1 -> 2 pairs grew 5.9x; 5 pairs ~2000 s on 2 Pi threads, /31 for 12-16 VM threads) scaled with natm^2.6.
VM_GFLOPS_EFFECTIVE = {16: 450.0, 20: 520.0, 40: 650.0}     # sustained DF exchange rate, unshared VM
XC_S_PER_CYCLE_5PAIRS_16T = 65.0


def size_estimate(name, basis=None):
    from pyscf import gto
    from pyscf.df import addons
    c = json.loads((MOL_DIR / 'cases' / f'{name}.json').read_text())
    mol = gto.M(atom=str(ROOT / c['xyz']), basis=basis or c['basis'], charge=c['charge'], verbose=0)
    naux = addons.make_auxmol(mol, addons.make_auxbasis(mol)).nao
    nao, nocc = mol.nao, mol.nelectron // 2
    return {'natm': mol.natm, 'nao': nao, 'naux': naux, 'nocc': nocc,
            'cderi_GB': naux * nao * (nao + 1) / 2 * 8 / 1e9,
            'K_flops_per_cycle': 4.0 * nocc * nao ** 2 * naux}


def stage_estimate(args):
    s = size_estimate(case_name('salt'))
    print(json.dumps(s, indent=1))
    xc16 = XC_S_PER_CYCLE_5PAIRS_16T * (s['natm'] / 135) ** 2.6
    for t, g in VM_GFLOPS_EFFECTIVE.items():
        k = s['K_flops_per_cycle'] / g / 1e9
        xc = xc16 * 16 / t / 0.8 if t > 16 else xc16
        per = k + xc
        print(f'[xclu] {t} threads (unshared VM): exchange ~{k:.0f} s + XC ~{xc:.0f} s = ~{per:.0f} s per cycle; '
              f'16-24 cycles + ~5 min setup = {(16 * per + 300) / 60:.0f}-{(24 * per + 300) / 60:.0f} min per structure')


# ------------------------------------------------------------------------------------------------------ SCF ----
def variant_case(tag, eps, prefix=CASE_PREFIX):
    """Case of structure `tag`, or its IEF-PCM(eps) twin (written from the gas-phase case on first use)."""
    base = case_name(tag, None, prefix)
    if not (MOL_DIR / 'cases' / f'{base}.json').exists():
        raise SystemExit(f'[xclu] case {base} missing: run the build stage first')
    if not eps:
        return base
    name = case_name(tag, eps, prefix)
    c = json.loads((MOL_DIR / 'cases' / f'{base}.json').read_text())
    c['title'] = c['title'].replace('gas phase', f'IEF-PCM eps = {eps:g}')
    c['pcm_eps'] = float(eps)
    save_json(MOL_DIR / 'cases' / f'{name}.json', c)
    return name


def pair_atom_lists(case):
    """Atom indices of each ion pair (cation + anion of the same pair number) in cluster order."""
    fr = case['fragments']
    pairs = {}
    for m in case['molecules']:
        pairs.setdefault(m['pair'], []).extend(fr[m['fragment']])
    return [sorted(v) for _, v in sorted(pairs.items())]


def pair_guess(mol, case, grid, log):
    """Block-diagonal initial density from separate SCFs of the ion pairs (same functional, basis, grid; gas)."""
    from pyscf import gto, dft
    ao = mol.aoslice_by_atom()
    dm = np.zeros((mol.nao, mol.nao))
    t0 = time.time()
    done = {}               # pairs that are lattice translations of each other have the same AO density matrix
    for k, atoms in enumerate(pair_atom_lists(case)):
        xyz = np.array([mol.atom_coord(i, unit='Angstrom') for i in atoms])
        key = (tuple(mol.atom_symbol(i) for i in atoms), tuple(np.round(xyz - xyz[0], 6).ravel()))
        idx = np.concatenate([np.arange(ao[i][2], ao[i][3]) for i in atoms])
        if key in done:
            dm[np.ix_(idx, idx)] = done[key][0]
            log.append(dict(done[key][1], pair=k, reused_from=done[key][1]['pair']))
            print(f'[xclu] pair guess {k}: translation of pair {done[key][1]["pair"]}, density reused', flush=True)
            continue
        sub = gto.M(atom=[(mol.atom_symbol(i), c) for i, c in zip(atoms, xyz)], unit='Angstrom',
                    basis=mol.basis, charge=0, spin=0, max_memory=mol.max_memory, verbose=0)
        mf = dft.RKS(sub).density_fit()
        mf.xc = os.environ.get('DFT_XC', 'b3lyp'); mf.grids.level = grid; mf.conv_tol = 1e-7; mf.max_cycle = 80
        mf.kernel()
        if sub.nao != len(idx):
            raise RuntimeError('pair guess: AO count mismatch')
        dmk = mf.make_rdm1()
        dm[np.ix_(idx, idx)] = dmk
        ent = {'pair': k, 'natm': len(atoms), 'nao': sub.nao, 'converged': bool(mf.converged),
               'cycles': int(getattr(mf, 'cycles', -1)), 'e_tot': float(mf.e_tot)}
        log.append(ent)
        done[key] = (dmk, ent)
        print(f'[xclu] pair guess {k}: {len(atoms)} atoms, E = {mf.e_tot:.8f} Ha, converged {mf.converged} '
              f'({time.time() - t0:.0f} s so far)', flush=True)
    return dm


def stage_scf(args, prefix=CASE_PREFIX):
    sys.path.insert(0, str(ROOT / 'dft'))
    from common import load_case, make_mf, stage_note, MAX_MEMORY, HARTREE_EV   # noqa: E402 (after pinning)
    import pyscf
    from pyscf import gto, lib
    from pyscf.tools import molden
    name = variant_case(args.tag, args.eps, prefix)
    case = load_case(name)                       # also exports DFT_XC / DFT_DISP from the case
    d = case['dir']
    nthreads = int(os.environ.get('DFT_THREADS', '8'))
    lib.num_threads(nthreads)
    grid = int(os.environ.get('GRID_LEVEL', '3'))
    basis = args.basis or case['basis']
    t_all = time.time()
    mol = gto.M(atom=str(ROOT / case['xyz']), basis=basis, charge=case['charge'], spin=case['spin'],
                max_memory=MAX_MEMORY, verbose=args.verbose)
    sz = size_estimate(name, basis)
    incore = sz['cderi_GB'] * 1e3 < 0.9 * MAX_MEMORY
    tmp = lib.param.TMPDIR
    free_GB = shutil.disk_usage(tmp).free / 1e9
    print(f'[xclu] {name}: {mol.natm} atoms, {mol.nao} basis functions ({basis}), {sz["naux"]} auxiliary, '
          f'{mol.nelectron} electrons, {nthreads} threads, max_memory {MAX_MEMORY} MB, grid level {grid}; DF tensor '
          f'{sz["cderi_GB"]:.1f} GB {"in memory" if incore else f"on disk in {tmp} ({free_GB:.0f} GB free)"}', flush=True)
    if not incore and free_GB < 1.1 * sz['cderi_GB']:
        raise SystemExit(f'[xclu] {name}: {free_GB:.0f} GB free in {tmp}, the DF tensor needs {sz["cderi_GB"]:.0f} GB')
    stage_note(case, f'SCF started: {mol.natm} atoms, {mol.nao} basis functions, {nthreads} threads')
    mf = make_mf(mol, conv_tol=args.conv_tol, grid_level=grid, chkfile=d / 'final.chk')
    env = 'gas phase'
    if args.eps:
        mf = mf.PCM()
        mf.with_solvent.method = 'IEF-PCM'
        mf.with_solvent.eps = float(args.eps)
        mf.with_solvent.lebedev_order = int(os.environ.get('PCM_LEBEDEV_ORDER', '17'))
        env = f'IEF-PCM eps = {args.eps:g}, Lebedev {mf.with_solvent.lebedev_order}'
    mf.max_cycle = 80
    guess_log, dm0, t_guess = [], None, 0.0
    if not args.no_frag_guess:
        t = time.time()
        try:
            dm0 = pair_guess(mol, case, grid, guess_log)
        except Exception as e:                    # the guess is an accelerator only
            print(f'[xclu] pair guess failed ({type(e).__name__}: {e}); using MINAO', flush=True)
            dm0 = None
        t_guess = time.time() - t
    t = time.time()
    strategy = ['DIIS']
    mf.kernel(dm0=dm0)
    ncyc = int(getattr(mf, 'cycles', 0))
    if not mf.converged:
        print('[xclu] not converged after 80 cycles: level shift 0.25 Ha, 60 more cycles', flush=True)
        strategy.append('level shift 0.25')
        mf.level_shift = 0.25; mf.max_cycle = 60
        mf.kernel(dm0=mf.make_rdm1())
        ncyc += int(getattr(mf, 'cycles', 0))
        mf.level_shift = 0.0
    if not mf.converged:
        print('[xclu] still not converged: second-order (Newton) SCF, 20 macro cycles', flush=True)
        strategy.append('Newton')
        mf2 = mf.newton(); mf2.max_cycle = 20
        mf2.kernel(mf.mo_coeff, mf.mo_occ)
        ncyc += int(getattr(mf2, 'cycles', 0))
        mf = mf2
    t_scf = time.time() - t
    if not mf.converged:
        stage_note(case, f'SCF NOT converged ({" -> ".join(strategy)})')
        raise SystemExit(f'[xclu] {name}: SCF did not converge')
    dm = mf.make_rdm1()
    nel = float(np.trace(dm @ mol.intor('int1e_ovlp')))
    nocc = mol.nelectron // 2
    e = mf.mo_energy * HARTREE_EV
    cderi = getattr(mf.with_df, '_cderi', None)
    res = {'case': name, 'structure': case.get('structure', args.tag), 'generated': now(),
           'level': {'xc': os.environ.get('DFT_XC'), 'disp': os.environ.get('DFT_DISP'), 'basis': basis,
                     'basis_print': BASIS_PRINT.get(basis, basis), 'spherical_d': not mol.cart, 'grid_level': grid,
                     'density_fitting': True, 'auxbasis': 'PySCF default for the basis (heavy-aug-cc-pVDZ-JKFIT for 6-31+G(d,p))',
                     'environment': env, 'conv_tol': args.conv_tol},
           'natm': mol.natm, 'nao': mol.nao, 'nmo': int(mf.mo_coeff.shape[1]), 'naux': int(mf.with_df.get_naoaux()),
           'nelectron': mol.nelectron, 'nelectron_from_density': nel, 'e_tot_Ha': float(mf.e_tot),
           'e_disp_Ha': float(mf.scf_summary.get('dispersion', 0.0)) if hasattr(mf, 'scf_summary') else None,
           'converged': bool(mf.converged), 'cycles': ncyc, 'strategy': strategy,
           'homo_eV': float(e[nocc - 1]), 'lumo_eV': float(e[nocc]), 'gap_eV': float(e[nocc] - e[nocc - 1]),
           'initial_guess': 'block-diagonal ion-pair densities' if dm0 is not None else 'MINAO', 'pair_guess': guess_log,
           'cderi_in_memory': not isinstance(cderi, str), 'threads': nthreads, 'max_memory_MB': MAX_MEMORY,
           't_guess_s': round(t_guess, 1), 't_scf_s': round(t_scf, 1), 'pyscf': pyscf.__version__,
           'xyz_sha256': sha256(ROOT / case['xyz'])}
    # Multiwfn needs a square MO matrix (see dft/fmo.py): pad linearly dependent vectors with empty orbitals
    nao, nmo = mf.mo_coeff.shape
    if nmo < nao:
        C = np.zeros((nao, nao)); C[:, :nmo] = mf.mo_coeff
        en = np.zeros(nao); en[:nmo] = mf.mo_energy; en[nmo:] = mf.mo_energy.max() + 1e3
        oc = np.zeros(nao); oc[:nmo] = mf.mo_occ
        molden.from_mo(mol, str(d / 'final.molden'), C, ene=en, occ=oc)
        res['molden_padded'] = nao - nmo
    else:
        molden.from_scf(mf, str(d / 'final.molden'))
    res['t_total_s'] = round(time.time() - t_all, 1)
    save_json(d / 'scf.json', res)
    msg = (f"SCF {' -> '.join(strategy)} converged in {res['cycles']} cycles, {t_scf / 60:.1f} min (guess {t_guess / 60:.1f} min); "
           f"E = {mf.e_tot:.8f} Ha; HOMO-LUMO gap {res['gap_eV']:.2f} eV; {mol.nao} basis functions, {res['naux']} auxiliary, "
           f"DF tensor {'in memory' if res['cderi_in_memory'] else 'on disk'}")
    stage_note(case, msg)
    print(f'[xclu] {name}: {msg}', flush=True)
    if abs(nel - mol.nelectron) > 1e-5:
        raise SystemExit(f'[xclu] density integrates to {nel:.6f} electrons, not {mol.nelectron}')
    if res['gap_eV'] < 1.0:
        print(f"[xclu] WARNING: HOMO-LUMO gap {res['gap_eV']:.2f} eV is small for a closed-shell salt cluster", flush=True)
    return res


# ---------------------------------------------------------------------------------------------------- QTAIM ----
BCP_KEYS = ('rho', 'lap', 'G', 'V', 'H', 'minusG_over_V', 'ellipticity', 'eta', 'E_HB_Espinosa_kcal')


def read_cps_csv(path):
    import csv
    out = {}
    for r in csv.DictReader(open(path)):
        if r.get('type') != '(3,-1)' or '-' not in (r.get('atoms') or ''):
            continue
        a, b = r['atoms'].split('-', 1)
        row = {}
        for k in ('rho', 'lap', 'G', 'V', 'H', 'ellipticity', 'eta', 'd_AB_A', 'bond_path_A'):
            try:
                row[k] = float(r[k])
            except (TypeError, ValueError):
                row[k] = None
        if row.get('G') is not None and row.get('V'):
            row['minusG_over_V'] = -row['G'] / row['V']
        row.update({'cp': int(r['cp']), 'class': r.get('class'), 'atoms': r['atoms']})
        out[frozenset((a, b))] = row
    return out


def extract_bcps(name):
    """The BCPs of every case H-bond instance from dft/mwfn_extra.py's qtaim.json / qtaim_cps.csv."""
    case = json.loads((MOL_DIR / 'cases' / f'{name}.json').read_text())
    d = MOL_DIR / 'results' / name
    q = json.loads((d / 'mwfn_extra' / 'qtaim.json').read_text())
    summ = json.loads((d / 'mwfn_extra' / 'summary.json').read_text())
    cps = read_cps_csv(d / 'mwfn_extra' / 'qtaim_cps.csv')
    rows = {r['case_hbond']: r for r in q['hbonds'] + q['other_noncovalent'] if r.get('case_hbond')}
    syms, X = [], []
    for l in (ROOT / case['xyz']).read_text().splitlines()[2:]:
        p = l.split()
        if len(p) >= 4:
            syms.append(p[0]); X.append([float(v) for v in p[1:4]])
    X = np.array(X)
    lab = case['cif_labels']
    inst = {}
    for hname, (iD, iH, iA) in case['hbonds'].items():
        r = rows.get(hname)
        g = hb_geometry(X, iD, iH, iA)
        ent = {'type': case['hbond_types'][hname], 'primary': hname in case['primary_hbonds'],
               'geometry': {k: round(v, 4) for k, v in g.items()},
               'hbond_bcp': None, 'bcp_H_A': cps.get(frozenset((lab[iH], lab[iA]))),
               'bcp_D_H': cps.get(frozenset((lab[iD], lab[iH])))}
        if r is not None:
            ent['hbond_bcp'] = {k: r.get(k) for k in BCP_KEYS} | {
                'multiwfn_name': r.get('name'), 'class': r.get('cls'), 'criteria_met': r.get('hbond_criteria_met'),
                'proton_transferred': 'reverse of case H-bond' in (r.get('name') or ''), 'cp': r.get('cp'),
                'bond_path_A': r.get('bond_path_A'), 'nature_Rozas2000': r.get('nature_Rozas2000')}
        inst[hname] = ent
    missing = [k for k, v in inst.items() if v['hbond_bcp'] is None]
    res = {'case': name, 'structure': case.get('structure'), 'pcm_eps': case.get('pcm_eps'), 'generated': now(),
           'wavefunction_sha256': summ.get('wavefunction_sha256'), 'multiwfn_version': summ.get('multiwfn_version'),
           'qtaim': {k: q.get(k) for k in ('n_cp', 'counts', 'poincare_hopf', 'cp_searches_logged',
                                             'n_noncovalent_bcp', 'n_bcp_without_bond_path')},
           'instances': inst, 'instances_without_hbond_bcp': missing,
           'primary': {t: h for h, t in case['primary_hbonds'].items()}}
    save_json(d / 'nine_bcps.json', res)
    return res


def stage_qtaim(args, prefix=CASE_PREFIX):
    name = case_name(args.tag, args.eps, prefix)
    d = MOL_DIR / 'results' / name
    if not (d / 'final.molden').exists() or not (d / 'scf.json').exists():
        raise SystemExit(f'[xclu] {name}: final.molden / scf.json missing: run the scf stage first')
    t0 = time.time()
    env = dict(os.environ)
    env.setdefault('MWFN_THREADS', os.environ.get('DFT_THREADS', '8'))
    env.setdefault('MWFN_BIN', str(ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI'))
    rc = subprocess.run([sys.executable, str(ROOT / 'dft' / 'mwfn_extra.py'), name, 'qtaim'], cwd=ROOT, env=env).returncode
    summ_p = d / 'mwfn_extra' / 'summary.json'
    st = json.loads(summ_p.read_text()).get('status', {}).get('qtaim', {}) if summ_p.exists() else {}
    if rc != 0 or not st.get('ok') or (d / 'mwfn_extra' / 'qtaim.json').stat().st_mtime < t0:
        raise SystemExit(f'[xclu] {name}: dft/mwfn_extra.py qtaim failed (exit {rc}, status {st})')
    res = extract_bcps(name)
    n = len(res['instances'])
    print(f"[xclu] {name}: QTAIM {res['qtaim']['n_cp']} CPs, Poincare-Hopf "
          f"{'OK' if (res['qtaim']['poincare_hopf'] or {}).get('satisfied') else 'NOT satisfied'}; "
          f"{n - len(res['instances_without_hbond_bcp'])}/{n} H-bond instances have a BCP "
          f"({time.time() - t0:.0f} s)", flush=True)
    for t, h in res['primary'].items():
        b = res['instances'][h]['hbond_bcp']
        print(f"         {h:30s} " + (f"rho {b['rho']:.4f}  lap {b['lap']:+.4f}  H {b['H']:+.5f}  eps {b['ellipticity']:.3f}"
                                        if b else 'NO BCP'), flush=True)
    if res['instances_without_hbond_bcp']:
        print(f"[xclu] NOTE: no H-bond BCP for {res['instances_without_hbond_bcp']}", flush=True)
    return res


# -------------------------------------------------------------------------------------------------- compare ----
def stage_compare(args, prefix=CASE_PREFIX, out_dir=None):
    out_dir = pathlib.Path(out_dir or OUT)
    got = {}
    for p in sorted((MOL_DIR / 'results').glob(f'{prefix}_*/nine_bcps.json')):
        r = json.loads(p.read_text())
        if (r.get('pcm_eps') or None) != (args.eps or None):
            continue
        got[r['structure']] = r
    if 'salt' not in got:
        raise SystemExit(f'[xclu] compare: no salt result (nine_bcps.json) for eps={args.eps}')
    others = [s for s in got if s != 'salt']
    scf = {s: json.loads((MOL_DIR / 'results' / got[s]['case'] / 'scf.json').read_text()) for s in got}
    rows, table = [], []
    salt = got['salt']
    for t in TYPE_NAMES:
        if t not in salt['primary']:
            continue
        h = salt['primary'][t]
        row = {'bond': t, 'instance': h, 'cif_code': dict((n, c) for n, (_, _, _, c, _) in zip(TYPE_NAMES, CIF_HBONDS))[t]}
        for s in ['salt'] + others:
            ent = got[s]['instances'].get(h)
            b = (ent or {}).get('hbond_bcp') or {}
            row[s] = {'geometry': (ent or {}).get('geometry'), 'bcp': b or None,
                      'bcp_is_reverse': bool(b.get('proton_transferred')),
                      'bcp_H_A': (ent or {}).get('bcp_H_A'), 'bcp_D_H': (ent or {}).get('bcp_D_H')}
            if s != 'salt' and b and row['salt']['bcp']:
                row[f'delta_{s}_minus_salt'] = {k: (b[k] - row['salt']['bcp'][k]) if (b.get(k) is not None and row['salt']['bcp'].get(k) is not None) else None
                                                for k in ('rho', 'lap', 'G', 'V', 'H', 'ellipticity', 'minusG_over_V')}
                row[f'delta_{s}_minus_salt']['rho_pct'] = 100 * (b['rho'] / row['salt']['bcp']['rho'] - 1)
                row[f'delta_{s}_minus_salt']['H..A_A'] = row[s]['geometry']['H..A'] - row['salt']['geometry']['H..A']
        # instance spread inside each cluster (edge-effect indicator)
        for s in ['salt'] + others:
            vals = [e['hbond_bcp']['rho'] for n, e in got[s]['instances'].items() if e['type'] == t and e['hbond_bcp']]
            row[f'rho_instances_{s}'] = vals
        rows.append(row)
    # cluster energy difference (all pairs transferred), per formula unit
    energy = {}
    for s in others:
        dE = (scf[s]['e_tot_Ha'] - scf['salt']['e_tot_Ha']) / (scf['salt']['natm'] / 27)
        energy[s] = {'dE_per_fu_kJ': dE * 2625.4996394799, 'dE_per_fu_kcal': dE * 627.509474,
                     'n_pairs': scf['salt']['natm'] // 27,
                     'note': 'B3LYP-D3(BJ) cluster energy at the fixed periodic PBE-D3(BJ) geometries, every pair of the '
                             'cluster in the same state; gas-phase cluster (surface ions under-coordinated); compare with '
                             'the periodic +57.12 (PBE-D3(BJ)) / +52.78 (PBE0-D3(BJ)) kJ/mol per f.u. only qualitatively'}
    res = {'generated': now(), 'eps': args.eps, 'structures': list(got), 'cases': {s: got[s]['case'] for s in got},
           'level': scf['salt']['level'], 'n_atoms': scf['salt']['natm'], 'nao': scf['salt']['nao'],
           'qtaim': {s: got[s]['qtaim'] for s in got}, 'bonds': rows, 'cluster_energy': energy,
           'units': 'atomic units (rho e/bohr^3; lap e/bohr^5; G, V, H hartree/bohr^3); E_HB kcal/mol; distances A',
           'notes': ['primary instance = the bond at the central ion pair (C0/A0); rho_instances_* lists every instance '
                     'of the same bond in the cluster (spread = cluster-edge sensitivity)',
                     'N1-H1A...O5 in the constrained neutral cell: the BCP is N1...H1A (proton on O5), marked bcp_is_reverse',
                     'E_HB = V/2 (Espinosa) is carried for completeness only (referee M6 point 2)']}
    tagstr = f'_eps{args.eps:g}' if args.eps else ''
    save_json(out_dir / f'qtaim_nine_bonds{tagstr}.json', res)
    # CSV + Markdown
    cols = ['bond', 'instance']
    for s in ['salt'] + others:
        cols += [f'{s}_H..A_A', f'{s}_rho', f'{s}_lap', f'{s}_G', f'{s}_V', f'{s}_H', f'{s}_ellipticity']
    for s in others:
        cols += [f'd_{s}_rho', f'd_{s}_rho_pct', f'd_{s}_lap', f'd_{s}_H', f'd_{s}_ellipticity']
    lines = [','.join(cols)]
    md = [f"# M6: QTAIM of the nine crystal H-bonds, cluster of {res['n_atoms']} atoms ({res['nao']} basis functions)",
          '', f"Level: {res['level']['xc']}-{res['level']['disp']}/{res['level']['basis_print']}, "
              f"{res['level']['environment']}, single points at the periodic geometries. Units: atomic units.", '',
          '| Bond | ' + ' | '.join(f'{s} rho | {s} lap | {s} H' for s in ['salt'] + others)
          + ''.join(f' | d rho ({s}) | d rho % | d H ({s})' for s in others) + ' |',
          '|---|' + '---|' * (3 * (1 + len(others)) + 3 * len(others))]
    f = lambda v, n=4: '' if v is None else f'{v:.{n}f}'
    for r in rows:
        vals = [r['bond'], r['instance']]
        cells = [r['bond'] + (' (rev.)' if any(r[s]['bcp_is_reverse'] for s in others) else '')]
        for s in ['salt'] + others:
            b, g = r[s]['bcp'] or {}, r[s]['geometry'] or {}
            vals += [f(g.get('H..A')), f(b.get('rho'), 5), f(b.get('lap'), 5), f(b.get('G'), 5), f(b.get('V'), 5),
                     f(b.get('H'), 6), f(b.get('ellipticity'), 4)]
            cells += [f(b.get('rho')), f(b.get('lap')), f(b.get('H'), 5)]
        for s in others:
            dd = r.get(f'delta_{s}_minus_salt') or {}
            vals += [f(dd.get('rho'), 5), f(dd.get('rho_pct'), 1), f(dd.get('lap'), 5), f(dd.get('H'), 6),
                     f(dd.get('ellipticity'), 4)]
            cells += [f(dd.get('rho'), 4), f(dd.get('rho_pct'), 1), f(dd.get('H'), 5)]
        lines.append(','.join(f'"{v}"' if ',' in str(v) else str(v) for v in vals))
        md.append('| ' + ' | '.join(cells) + ' |')
    for s, e in energy.items():
        md += ['', f"Cluster energy {s} - salt: {e['dE_per_fu_kJ']:+.2f} kJ/mol per formula unit "
                   f"({e['n_pairs']} pairs). {e['note']}."]
    (out_dir / f'qtaim_nine_bonds{tagstr}.csv').write_text('\n'.join(lines) + '\n')
    (out_dir / f'qtaim_nine_bonds{tagstr}.md').write_text('\n'.join(md) + '\n')
    print('\n'.join(md))
    print(f'[xclu] wrote {out_dir}/qtaim_nine_bonds{tagstr}.{{json,csv,md}}')
    return res


# ----------------------------------------------------------------------------------------------------- test ----
MOCK_MULTIWFN = r'''#!{python}
# MOCK of Multiwfn for the Pi self-test of dft/crystal_cluster_qtaim.py ONLY (Multiwfn is x86-only).  It ignores
# the menu input and writes CPprop.txt / CPs.txt / a log in Multiwfn's formats: nuclear CPs, and one (3,-1) CP at the
# midpoint of every covalent bond and every H...N/O contact < 2.7 A.  rho at that point is the REAL density of the
# molden wavefunction (PySCF); every other property is synthetic.  Never use these numbers.
import sys, re, numpy as np
from pyscf.tools import molden
wfn = sys.argv[1]
sys.stdin.read()
mol, ene, C, occ, _, _ = molden.load(wfn)
Z = [mol.atom_pure_symbol(i) for i in range(mol.natm)]
R = mol.atom_coords(unit='Bohr'); B = 0.529177210903
cov = {{'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05}}
dm = (C * occ) @ C.T
def rho(p):
    ao = mol.eval_gto('GTOval_cart' if mol.cart else 'GTOval_sph', np.array([p]))
    return float(np.einsum('pi,ij,pj->p', ao, dm, ao)[0])
cps = [('(3,-3)', i, None, R[i]) for i in range(mol.natm)]
for i in range(mol.natm):
    for j in range(i + 1, mol.natm):
        d = np.linalg.norm(R[i] - R[j]) * B
        hb = ('H' in (Z[i], Z[j])) and ({{Z[i], Z[j]}} & {{'N', 'O'}}) and d < 2.7
        if d < 1.25 * (cov[Z[i]] + cov[Z[j]]) or hb:
            cps.append(('(3,-1)', i, j, (R[i] + R[j]) / 2))
with open('CPprop.txt', 'w') as f:
    for k, (t, i, j, p) in enumerate(cps, 1):
        f.write(f' ----------------   CP {{k:5d}},     Type {{t}}   ----------------\n')
        if j is None:
            f.write(f' Corresponding nucleus: {{i+1:5d}}({{Z[i]:2s}})\n')
        else:
            f.write(f' Connected atoms: {{i+1:5d}}({{Z[i]:2s}})   -- {{j+1:5d}}({{Z[j]:2s}})\n')
        a = p * B
        f.write(f' Position (Bohr): {{p[0]:20.12f}}{{p[1]:20.12f}}{{p[2]:20.12f}}\n')
        f.write(f' Position (Angstrom): {{a[0]:20.12f}}{{a[1]:20.12f}}{{a[2]:20.12f}}\n')
        r = rho(p); G = 0.8 * r ** 1.3; V = -1.1 * G; lap = 4 * (G + V / 2)
        for lab, v in (('Density of all electrons', r), ('Lagrangian kinetic energy G(r)', G),
                       ('Hamiltonian kinetic energy K(r)', -(G + V)), ('Potential energy density V(r)', V),
                       ('Energy density E(r) or H(r)', G + V), ('Laplacian of electron density', lap),
                       ('Electron localization function (ELF)', 0.1), ('Localized orbital locator (LOL)', 0.2),
                       ('Reduced density gradient (RDG)', 0.0), ('Sign(lambda2)*rho', -r),
                       ('Delta-g (under Hirshfeld partition)', 0.01)):
            f.write(f' {{lab}}: {{v:18.10E}}\n')
        f.write(f' Total ESP: {{0.1:18.10E}} a.u. ( 0.0 eV, 0.0 kcal/mol)\n')
        f.write(f' Eigenvalues of Hessian: {{-r:18.10E}}{{-0.9 * r:18.10E}}{{3 * r:18.10E}}\n')
        f.write(f' Ellipticity of electron density: {{0.05:12.6f}}\n eta index: {{0.3:12.6f}}\n\n')
with open('CPs.txt', 'w') as f:
    f.write(f'{{len(cps):7d}}\n')
    for k, (t, i, j, p) in enumerate(cps, 1):
        f.write(f'{{k:6d}} {{p[0]:12.6f}} {{p[1]:12.6f}} {{p[2]:12.6f}} 1\n')
nb = sum(t == '(3,-1)' for t, *_ in cps)
print(' Multiwfn MOCK Version test')
for n in (mol.natm, nb, 0, 0, 0):
    print(f' Totally find {{n:5d}} new critical points')
print(f' Poincare-Hopf verification: {{mol.natm}} - {{nb}} + {{nb - mol.natm + 1}} - 0 = 1')
for k, (t, i, j, p) in enumerate(cps, 1):
    if t == '(3,-1)':
        for n in (i, j):
            print(f' Path {{k:4d}},    CP: {{k:4d}} (3,-1) --> CP: {{n + 1:4d}} (3,-3)   Length: {{np.linalg.norm(p - R[n]):8.5f}} Bohr')
'''


def stage_test(args):
    """Pi self-test (seconds to minutes on 2 threads): the build checks for the production clusters, then the whole
    chain scf -> dft/mwfn_extra.py qtaim (MOCK Multiwfn) -> extraction -> compare on a 2-pair test cluster with
    STO-3G, written under the case prefix sfa_xclutest and removed afterwards (--keep keeps it)."""
    prefix = 'sfa_xclutest'
    salt = cell_molecules(STRUCTS['salt'])
    pool = Pool(salt)
    sel = select_cluster(pool, 'minimal')          # the production default
    checks = {}
    # 1. production build checks (in memory only: nothing written)
    for tag, path in STRUCTS.items():
        cell = salt if tag == 'salt' else map_structure(salt, path)
        syms, X, labs, frags, tags = build_cluster(pool, cell, sel['pairs'])
        dfrac, dangling, worst = periodic_check(cell, pool, sel['pairs'], X, syms)
        ins = find_instances(pool, salt, [m for p in sel['pairs'] for m in pool.pair_members(p)])
        types = {x['type'] for x in ins}
        nq = sum(1 if f.startswith('cation') else -1 for f in frags)
        checks[f'build_{tag}'] = (dfrac < 1e-8 and dangling == 0 and types == set(range(9)) and nq == 0
                                  and electron_count(syms) % 2 == 0)
        print(f'[xclu-test] build {tag}: {len(syms)} atoms, translations integer to {dfrac:.0e}, {dangling} cut bonds, '
              f'{len(types)}/9 bond types, net charge {nq} -> {"OK" if checks[f"build_{tag}"] else "FAIL"}')
    # 2. tiny end-to-end chain on the central pair + its first partner pair
    tiny = sel['pairs'][:2]
    made = make_cases(pool, salt, tiny, STRUCTS, prefix=prefix, cover='test', require_all=False, basis='sto-3g')
    mock_dir = pathlib.Path(tempfile.mkdtemp(prefix='xclu_mock_'))
    try:
        (mock_dir / 'Multiwfn_noGUI').write_text(MOCK_MULTIWFN.format(python=sys.executable))
        (mock_dir / 'Multiwfn_noGUI').chmod(0o755)
        (mock_dir / 'settings.ini').write_text('  nthreads=  4\n  iESPcode= 2\n')
        os.environ['MWFN_BIN'] = str(mock_dir / 'Multiwfn_noGUI')
        os.environ['MWFN_THREADS'] = '1'
        os.environ.setdefault('GRID_LEVEL', '1')
        ns = argparse.Namespace(eps=None, basis=None, conv_tol=1e-8, no_frag_guess=False, verbose=0)
        for tag in ('salt', 'neutral'):
            ns.tag = tag
            r = stage_scf(ns, prefix=prefix)
            checks[f'scf_{tag}'] = r['converged'] and abs(r['nelectron_from_density'] - r['nelectron']) < 1e-5
            q = stage_qtaim(ns, prefix=prefix)
            n_case = len(made[tag]['hbond_instances'])
            ok = len(q['instances']) == n_case and not q['instances_without_hbond_bcp']
            if tag == 'neutral':
                h1a = [v for k, v in q['instances'].items() if k.startswith('N1-H1A...O5 [C0>A0]')][0]
                ok = ok and h1a['hbond_bcp']['proton_transferred'] and h1a['bcp_D_H']['class'] != 'covalent'
            checks[f'qtaim_{tag}'] = ok
        test_out = mock_dir / 'compare'
        test_out.mkdir()
        cmp = stage_compare(argparse.Namespace(eps=None), prefix=prefix, out_dir=test_out)
        rho_s = {r['bond']: r['salt']['bcp']['rho'] for r in cmp['bonds']}
        checks['compare'] = len(cmp['bonds']) == len(set(i['type'] for i in made['salt']['hbond_instances'])) and \
            all(r.get('delta_neutral_minus_salt') for r in cmp['bonds'])
        # physics sanity even with STO-3G: the shortest contact (O4-H4A...O5, 1.49 A) has the largest H...A density
        checks['rho_order'] = max(rho_s, key=rho_s.get) in ('O4-H4A...O5', 'N1-H1A...O5')
        print(f'[xclu-test] test rho at the H...A BCPs (STO-3G, mock topology): '
              + ', '.join(f'{k} {v:.4f}' for k, v in rho_s.items()))
    finally:
        shutil.rmtree(mock_dir, ignore_errors=True)
        if not args.keep:
            for tag in STRUCTS:
                n = case_name(tag, None, prefix)
                (MOL_DIR / 'cases' / f'{n}.json').unlink(missing_ok=True)
                (MOL_DIR / 'input' / f'{n}.xyz').unlink(missing_ok=True)
                shutil.rmtree(MOL_DIR / 'results' / n, ignore_errors=True)
    print('[xclu-test] ' + ', '.join(f'{k} {"OK" if v else "FAIL"}' for k, v in checks.items()))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='M6: QTAIM of the nine crystal H-bonds on a cluster cut from the periodic cell')
    ap.add_argument('stage', choices=['build', 'estimate', 'scf', 'qtaim', 'compare', 'test'])
    ap.add_argument('tag', nargs='?', help='structure: salt, neutral, or an --extra tag')
    ap.add_argument('--cover', choices=['minimal', 'central'], default='minimal')
    ap.add_argument('--only', action='append', help='build: write only these structure tags (selection unchanged)')
    ap.add_argument('--extra', action='append', help='TAG=path.extxyz: another structure of the same cell (build)')
    ap.add_argument('--eps', type=float, help='IEF-PCM dielectric constant (scf/qtaim/compare); default gas phase')
    ap.add_argument('--basis', help='override the case basis (tests only)')
    ap.add_argument('--conv-tol', type=float, default=1e-8)
    ap.add_argument('--no-frag-guess', action='store_true')
    ap.add_argument('--verbose', type=int, default=4)
    ap.add_argument('--keep', action='store_true', help='test: keep the sfa_xclutest cases and results')
    args = ap.parse_args()
    if args.stage in ('scf', 'qtaim') and not args.tag:
        ap.error('scf/qtaim need a structure tag (salt, neutral, ...)')
    {'build': stage_build, 'estimate': stage_estimate, 'scf': stage_scf, 'qtaim': stage_qtaim,
     'compare': stage_compare, 'test': stage_test}[args.stage](args)
