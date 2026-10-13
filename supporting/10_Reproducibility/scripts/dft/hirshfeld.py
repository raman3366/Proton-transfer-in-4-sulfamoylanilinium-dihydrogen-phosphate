"""
Crystal Hirshfeld surfaces and 2-D fingerprint plots for SAPH (CCDC 1016950) straight from the deposited CIF.

Why this script exists and why it looks the way it does:
  * Hirshfeld-surface analysis is normally done in CrystalExplorer, which has no headless Linux/ARM build we can
    script.  Multiwfn (main function 12 -> surface type 5 "Hirshfeld surface") computes the surface, d_i, d_e,
    d_norm and element-pair contact areas, but it has no shape index or curvedness and it runs only on x86.
  * The surface needs NO wavefunction: it is the isosurface w(r) = 0.5 of the Hirshfeld weight
        w(r) = rho_in(r) / rho_all(r),  rho = sum of spherical free-atom densities (the "promolecule"),
    so everything can be computed here from the CIF alone, on the Pi, before any DFT is finished.  Multiwfn is
    then run on the VM on the identical cluster as an independent cross-check (step 'mwfn').
  * Free-atom densities: ground-state ROHF atoms (H 2S, C 3P, N 4S, O 3P, P 4S, S 3P) with aug-cc-pVQZ,
    spherically averaged on a 302-point Lebedev shell: near-HF-limit free atoms (E within 2 mHa of the
    numerical HF limits; O and S slightly below them because the ROHF solution breaks spherical symmetry).
    CrystalExplorer also uses tabulated free-atom densities, but which ones was NOT checked (unverified).
    PySCF's own atom_hf.get_atm_nrhf() is NOT used: its spin-restricted fractional-occupation atoms sit
    0.35-0.52 Ha above ROHF for C/O and their O density is 25-45 % too diffuse at 1.5-2 A, i.e. exactly where
    the surface lies (checked 28 Sep 2026).
  * Conventions, chosen to match CrystalExplorer so the numbers are comparable with the literature:
      - X-H bonds normalised to neutron lengths C-H 1.083, N-H 1.009, O-H 0.983 A (CrystalExplorer FAQ).  The CIF
        O-H were DFIX-restrained to 0.82 A, so un-normalised surfaces would be meaningless around H3A/H4A.
      - d_norm = (d_i - r_i)/r_i + (d_e - r_e)/r_e with Bondi radii H 1.20, C 1.70, N 1.55, O 1.52, P 1.80,
        S 1.80 A (the table in CrystalExplorer's back-end library occ, include/occ/core/element.h).
      - shape index S = (2/pi) atan[(k1+k2)/(k1-k2)], k1 >= k2 principal curvatures (1/A) with the sign chosen
        so that a convex bump (seen from outside) has k > 0: S = +1 convex cap, S = -1 concave cup (where a
        neighbour's atom pokes into the surface), S = 0 saddle (tested on a sphere: k = 1/R, S = +1).  The usual
        CrystalExplorer colouring (red hollows, blue bumps) corresponds to red = -1 -- not checked in CE itself.
        curvedness C = (2/pi) ln sqrt[(k1^2 + k2^2)/2] (Spackman & Jayatilaka, CrystEngComm 2009, 11, 19,
        doi:10.1039/B818330A; the log form -- SAPH values span -3.8 .. +0.7).
      - curvatures are ANALYTIC (gradient and Hessian of w from the radial splines), not estimated from the
        triangle mesh, so they do not depend on the grid spacing.
      - fingerprint: area-weighted 2-D histogram of (d_i, d_e), 0.01 A bins; a contact X...Y means nearest
        interior atom X, nearest exterior atom Y; "O...H/H...O" sums both directions (reciprocal contacts).
  * Model: the central formula unit is the asymmetric-unit cation plus the anion at (x, y, z-1), the same pairing
    as DFT model M1 (N1-H1A...O5, labels "#1_554").  Every symmetry-generated ion with any atom within --cutoff of
    that unit is added whole.  The per-ion surfaces do not depend on which formula unit is called central.

Steps (default: build surface):
  build   : gemmi symmetry expansion -> cluster.xyz, cluster.json (labels with SHELX codes n_klm, molecule ids),
            and a check that all 9 CIF hydrogen bonds are reproduced (D...A to 0.005 A).
  surface : for each target ion (cation, anion; 'pair' = whole formula unit on request) -> <t>_surface.npz
            (vertices, faces, normals, d_i, d_e, d_norm, shape index, curvedness, vertex areas, nearest atoms),
            <t>_vertices.csv.gz, <t>_fingerprint.npz, <t>_contacts.csv, <t>_neighbours.csv, <t>_red_spots.csv,
            hirshfeld_summary.json.  With --png also quick-look fingerprint plots.
  mwfn    : Multiwfn 12 -> 1 -> 5 (Hirshfeld surface of the same atoms, d_norm mapped) + 20 -> 3 (element-pair
            contact areas) on cluster.xyz; parsed into mwfn_<t>.json.  x86 only (Multiwfn binary).
  compare : Python vs Multiwfn area, volume and contact percentages -> hirshfeld_compare.json.
            Measured 28 Sep 2026: areas 189.61/189.62 (cation) and 90.36/91.48 A^2 (anion), volumes within 0.6/1.6 %,
            every contact percentage within 0.2 points; d_i and d_e agree to ~0.01 A.  Multiwfn's d_norm does NOT
            follow the CrystalExplorer formula (its minimum is -0.38 where the formula with Bondi radii gives -0.72
            at the same point; a least-squares fit of radii to its values fails, rms 0.13), so only our d_norm is reported.

Usage (from the repo root):
  python dft/hirshfeld.py [build] [surface] [mwfn] [compare] [--cutoff 6] [--spacing 0.12] [--targets cation,anion]
Output: molecules/sfa/results/hirshfeld_xray/ (override with --out)
"""
import argparse, gzip, json, os, re, subprocess, sys, time, pathlib, itertools
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import ROOT, COV_R, BOHR_A, save_json, write_xyz

CIF = ROOT / 'e-70-0o997-sup1.cif'
OUT_DEFAULT = ROOT / 'molecules' / 'sfa' / 'results' / 'hirshfeld_xray'
XH_NEUTRON = {'C': 1.083, 'N': 1.009, 'O': 0.983, 'B': 1.180}          # CrystalExplorer FAQ values, A
BONDI = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80, 'F': 1.47, 'Cl': 1.75}  # occ element.h
ROHF_UNPAIRED = {'H': 1, 'C': 2, 'N': 3, 'O': 2, 'P': 3, 'S': 2, 'F': 1, 'Cl': 1}     # ground-state 2S+1 - 1
ATOM_BASIS = os.environ.get('HS_ATOM_BASIS', 'aug-cc-pvqz')
R_GRID_CUT = 5.0     # A: atoms farther than this from a grid point are ignored when building the w grid
R_POINT_CUT = 7.0    # A: ... and when refining vertices / evaluating curvature (O at 7 A: rho ~ 1e-10 a.u.)
FP_BINS = np.round(np.arange(0.40, 3.0001, 0.01), 2)                  # fingerprint bin edges, A


# ============================================================ free-atom densities
def atomic_density_table(elements, out):
    """ln(rho) of spherically averaged ROHF free atoms on a uniform radial grid (A), cached in out/."""
    cache = out / f'atom_densities_{ATOM_BASIS}.npz'
    have = {}
    if cache.exists():
        z = np.load(cache)
        have = {k: z[k] for k in z.files if k != 'r'}
        r = z['r']
    need = [e for e in elements if e not in have]
    if need:
        from pyscf import gto, scf, lib
        from pyscf.dft import gen_grid
        lib.num_threads(int(os.environ.get('OMP_NUM_THREADS', '2')))
        r = np.linspace(0.0, 12.0, 6001)                       # 0.002 A steps
        dirs = gen_grid.MakeAngularGrid(302)[:, :3]            # Lebedev 302: exact average up to l = 29
        for el in need:
            m = gto.M(atom=f'{el} 0 0 0', basis=ATOM_BASIS, spin=ROHF_UNPAIRED[el], verbose=0)
            mf = scf.ROHF(m); mf.conv_tol = 1e-10; mf.max_cycle = 200; mf.kernel()
            if not mf.converged:
                raise RuntimeError(f'ROHF for free atom {el} did not converge')
            dm = mf.make_rdm1(); dm = dm[0] + dm[1] if dm.ndim == 3 else dm
            rho = np.empty(len(r))
            for k in range(0, len(r), 500):                     # chunks keep the AO arrays small
                rr = r[k:k + 500] / BOHR_A
                pts = (rr[:, None, None] * dirs[None]).reshape(-1, 3)
                ao = m.eval_gto('GTOval_sph', pts)
                rho[k:k + 500] = np.einsum('pi,ij,pj->p', ao, dm, ao).reshape(len(rr), -1).mean(1)
            have[el] = np.log(np.maximum(rho, 1e-300))
            print(f'[hs] free atom {el}: ROHF/{ATOM_BASIS} E = {mf.e_tot:.6f} Ha, '
                  f'N = {4*np.pi*np.trapezoid(rho*(r/BOHR_A)**2, r/BOHR_A):.4f}', flush=True)
        np.savez(cache, r=r, **have)
    return r, have


class Promolecule:
    """Sum of spherical atomic densities; value, gradient and Hessian from cubic splines of ln(rho) (A units)."""

    def __init__(self, r, lnrho):
        from scipy.interpolate import CubicSpline
        self.r, self.ln = r, lnrho
        self.sp = {e: CubicSpline(r, v) for e, v in lnrho.items()}

    def rho_grid(self, el, d):
        return np.exp(np.interp(d, self.r, self.ln[el]))

    def radial(self, el, d):
        """rho, rho', rho'' at distances d (A)."""
        sp = self.sp[el]
        f, f1, f2 = sp(d), sp(d, 1), sp(d, 2)
        rho = np.exp(f)
        return rho, f1 * rho, (f2 + f1 * f1) * rho

    def at_points(self, P, X, els, inside_mask, want_hess=False, rcut=R_POINT_CUT, chunk=3000):
        """rho_in, rho_all (+ gradients, Hessians) at points P from atoms X (els), inside_mask marks interior atoms.
        Processed in chunks of points: the per-pair Hessian arrays would otherwise take GBs."""
        if len(P) > chunk:
            parts = [self.at_points(P[k:k + chunk], X, els, inside_mask, want_hess, rcut, chunk)
                     for k in range(0, len(P), chunk)]
            return {key: np.concatenate([q[key] for q in parts]) for key in parts[0]}
        from scipy.spatial import cKDTree
        n = len(P)
        out = {k: np.zeros(n) for k in ('A', 'T')}
        out.update({k: np.zeros((n, 3)) for k in ('gA', 'gT')})
        if want_hess:
            out.update({k: np.zeros((n, 3, 3)) for k in ('HA', 'HT')})
        tree = cKDTree(X)
        pairs = tree.query_ball_point(P, rcut)
        ip = np.repeat(np.arange(n), [len(p) for p in pairs])
        ja = np.concatenate([np.asarray(p, dtype=int) for p in pairs]) if n else np.zeros(0, int)
        u = P[ip] - X[ja]; d = np.linalg.norm(u, axis=1); d = np.maximum(d, 1e-8); uh = u / d[:, None]
        rho = np.zeros(len(d)); r1 = np.zeros(len(d)); r2 = np.zeros(len(d))
        E = np.asarray(els)[ja]
        for el in np.unique(E):
            s = E == el
            rho[s], r1[s], r2[s] = self.radial(el, d[s])
        g = r1[:, None] * uh
        ins = np.asarray(inside_mask)[ja]
        np.add.at(out['T'], ip, rho); np.add.at(out['gT'], ip, g)
        np.add.at(out['A'], ip[ins], rho[ins]); np.add.at(out['gA'], ip[ins], g[ins])
        if want_hess:
            uu = uh[:, :, None] * uh[:, None, :]
            H = r2[:, None, None] * uu + (r1 / d)[:, None, None] * (np.eye(3)[None] - uu)
            np.add.at(out['HT'], ip, H); np.add.at(out['HA'], ip[ins], H[ins])
        return out


def weight_and_derivs(v, want_hess=False):
    """w = A/T, grad w, Hessian of w from the accumulated promolecular sums (see module docstring)."""
    A, T = v['A'], np.maximum(v['T'], 1e-300)
    w = A / T
    gw = (v['gA'] - w[:, None] * v['gT']) / T[:, None]
    if not want_hess:
        return w, gw, None
    Hw = (v['HA'] - w[:, None, None] * v['HT']
          - gw[:, :, None] * v['gT'][:, None, :] - v['gT'][:, :, None] * gw[:, None, :]) / T[:, None, None]
    return w, gw, Hw


# ============================================================ cluster from the CIF
def build_cluster(out, cutoff, normalise_h=True):
    import gemmi
    ss = gemmi.make_small_structure_from_block(gemmi.cif.read(str(CIF)).sole_block())
    cell = ss.cell
    M = np.array(cell.orth.mat.tolist()); Minv = np.linalg.inv(M)
    labels = [s.label for s in ss.sites]; els = [s.element.name for s in ss.sites]
    F = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in ss.sites])
    X = F @ M.T
    # molecules of the asymmetric unit by covalent connectivity (both ions are connected inside the ASU here)
    n = len(labels); par = list(range(n))
    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]; i = par[i]
        return i
    for i, j in itertools.combinations(range(n), 2):
        if np.linalg.norm(X[i] - X[j]) < 1.25 * (COV_R[els[i]] + COV_R[els[j]]):
            par[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    mols = sorted(groups.values(), key=lambda g: -len(g))
    kinds = ['anion' if any(els[i] == 'P' for i in g) else 'cation' for g in mols]
    if sorted(kinds) != ['anion', 'cation']:
        raise RuntimeError(f'expected one cation and one anion in the asymmetric unit, got {kinds}')
    # neutron-normalise X-H along the bond (in the ASU, before symmetry: isometries preserve it)
    h_moves = []
    if normalise_h:
        for i in range(n):
            if els[i] != 'H':
                continue
            dist = [(np.linalg.norm(X[i] - X[j]), j) for j in range(n) if j != i and els[j] != 'H']
            dmin, j = min(dist)
            if els[j] in XH_NEUTRON and dmin < 1.3:
                new = X[j] + (X[i] - X[j]) * XH_NEUTRON[els[j]] / dmin
                h_moves.append({'H': labels[i], 'parent': labels[j], 'd_xray_A': round(float(dmin), 4),
                                'd_normalised_A': XH_NEUTRON[els[j]]})
                X[i] = new
        F = X @ Minv.T
    ops = [gemmi.Op(s) for s in ss.symops]
    # every image (op, lattice translation) of every ASU molecule
    images = []
    rng = range(-3, 4)
    for k, op in enumerate(ops):
        for m, g in enumerate(mols):
            fr = np.array([op.apply_to_xyz(list(F[i])) for i in g])
            for t in itertools.product(rng, rng, rng):
                images.append((m, k, t, (fr + np.array(t)) @ M.T))
    code = lambda k, t: f'{k + 1}_{5 + t[0]}{5 + t[1]}{5 + t[2]}'
    ci = kinds.index('cation'); ai = kinds.index('anion')
    central = {(ci, 0, (0, 0, 0)): 'cation', (ai, 0, (0, 0, -1)): 'anion'}
    cen_xyz = np.vstack([xyz for (m, k, t, xyz) in images if (m, k, t) in central])
    atoms = []                                   # (label, el, xyz, mol_id, kind, symcode, central_role)
    mol_meta = []
    order = sorted(images, key=lambda im: (0 if (im[0], im[1], im[2]) in central else 1,
                                           0 if (im[0], im[1], im[2]) == (ci, 0, (0, 0, 0)) else 1))
    for m, k, t, xyz in order:
        key = (m, k, t)
        dmin = np.min(np.linalg.norm(xyz[:, None, :] - cen_xyz[None], axis=2))
        role = central.get(key)
        if role is None and dmin > cutoff:
            continue
        sc = code(k, t)
        mid = len(mol_meta)
        mol_meta.append({'id': mid, 'kind': kinds[m], 'symcode': sc, 'central': role,
                         'min_dist_to_central_A': round(float(dmin), 3)})
        for i, x in zip(mols[m], xyz):
            lab = labels[i] if sc == '1_555' else f'{labels[i]}#{sc}'
            atoms.append({'label': lab, 'el': els[i], 'xyz': x, 'mol': mid, 'kind': kinds[m], 'symcode': sc,
                          'role': role})
    syms = [a['el'] for a in atoms]; xyz = np.array([a['xyz'] for a in atoms])
    # clash check: no two atoms of different molecules closer than 1.5 A (would indicate a wrong op)
    from scipy.spatial import cKDTree
    mol_of = np.array([a['mol'] for a in atoms])
    bad = [(i, j) for i, j in cKDTree(xyz).query_pairs(1.5) if mol_of[i] != mol_of[j]]
    if bad:
        raise RuntimeError(f'intermolecular clash in the cluster: {[(atoms[i]["label"], atoms[j]["label"]) for i, j in bad[:5]]}')
    hb = check_cif_hbonds(atoms, xyz)
    write_xyz(out / 'cluster.xyz', syms, xyz, f'SAPH X-ray cluster, cutoff {cutoff} A, X-H normalised={normalise_h}')
    meta = {'cif': CIF.name, 'cell': [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma],
            'cell_volume_A3': cell.volume, 'Z': 2, 'symops': list(ss.symops), 'cutoff_A': cutoff,
            'xh_normalised': normalise_h, 'xh_lengths_A': XH_NEUTRON if normalise_h else None, 'h_moves': h_moves,
            'n_atoms': len(atoms), 'n_molecules': len(mol_meta), 'molecules': mol_meta,
            'atoms': [{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in a.items()} for a in atoms],
            'cif_hbonds_check': hb}
    save_json(out / 'cluster.json', meta)
    nc = sum(m['kind'] == 'cation' for m in mol_meta); na = len(mol_meta) - nc
    print(f'[hs] cluster: {len(atoms)} atoms, {len(mol_meta)} ions ({nc} cations, {na} anions) within {cutoff} A of '
          f'the central formula unit; CIF H-bonds reproduced {sum(h["found"] for h in hb)}/{len(hb)}', flush=True)
    return meta


def check_cif_hbonds(atoms, xyz):
    """Every D-H...A of the CIF _geom_hbond loop with D in the central unit must appear at the CIF D...A distance."""
    txt = CIF.read_text()
    m = re.search(r'_geom_hbond_angle_DHA\s*\n(.*?)\n_', txt, re.S)
    rows = [l.split() for l in m.group(1).strip().splitlines()] if m else []
    res = []
    for D, H, A, sc, dh, ha, da, ang in rows:
        da = float(da.split('(')[0])
        cand = [i for i, a in enumerate(atoms) if a['role'] and a['label'].split('#')[0] == D]
        best = None
        for i in cand:
            for j, a in enumerate(atoms):
                if a['label'].split('#')[0] == A and a['mol'] != atoms[i]['mol']:
                    d = float(np.linalg.norm(xyz[i] - xyz[j]))
                    if best is None or abs(d - da) < abs(best[0] - da):
                        best = (d, atoms[i]['label'], a['label'])
        res.append({'hbond': f'{D}-{H}...{A} [{sc}]', 'DA_cif_A': da, 'DA_cluster_A': round(best[0], 4) if best else None,
                    'pair': f'{best[1]}...{best[2]}' if best else None,
                    'found': bool(best and abs(best[0] - da) < 0.005)})
    return res


# ============================================================ surface
def marching(w, origin, h):
    from skimage.measure import marching_cubes
    verts, faces, _, _ = marching_cubes(w, level=0.5, spacing=(h, h, h), allow_degenerate=False)
    return verts + origin, faces


def curvatures(gw, Hw):
    """Principal curvatures (k1 >= k2, 1/A) of the level set w = 0.5 with outward normal n = -grad w/|grad w|."""
    gn = np.linalg.norm(gw, axis=1)
    n = -gw / gn[:, None]
    a = np.where(np.abs(n[:, 0:1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    t1 = a - np.sum(a * n, 1)[:, None] * n; t1 /= np.linalg.norm(t1, axis=1)[:, None]
    t2 = np.cross(n, t1)
    K = -Hw / gn[:, None, None]                  # F = -w increases outward: K = P (Hess F) P / |grad F|
    k11 = np.einsum('ni,nij,nj->n', t1, K, t1); k22 = np.einsum('ni,nij,nj->n', t2, K, t2)
    k12 = np.einsum('ni,nij,nj->n', t1, K, t2)
    mean = (k11 + k22) / 2; rad = np.sqrt(((k11 - k22) / 2) ** 2 + k12 ** 2)
    return mean + rad, mean - rad, n


def surface_for_target(target, meta, pm, out, spacing, png=False):
    t0 = time.time()
    atoms = meta['atoms']
    X = np.array([a['xyz'] for a in atoms]); E = [a['el'] for a in atoms]
    lab = [a['label'] for a in atoms]
    roles = {'cation': ('cation',), 'anion': ('anion',), 'pair': ('cation', 'anion')}[target]
    inside = np.array([a['role'] in roles for a in atoms])
    Xi = X[inside]
    pad = 3.2                                    # A beyond the interior nuclei; the surface lies at <= ~2.5 A
    lo = Xi.min(0) - pad; hi = Xi.max(0) + pad
    shape = tuple(int(np.ceil((hi[k] - lo[k]) / spacing)) + 1 for k in range(3))
    axes = [lo[k] + spacing * np.arange(shape[k]) for k in range(3)]
    A = np.zeros(shape); T = np.zeros(shape)
    for i in range(len(X)):
        sl = []
        for k in range(3):
            a0 = int(np.searchsorted(axes[k], X[i, k] - R_GRID_CUT)); a1 = int(np.searchsorted(axes[k], X[i, k] + R_GRID_CUT))
            if a1 <= a0:
                break
            sl.append(slice(a0, a1))
        if len(sl) < 3:
            continue
        dx = axes[0][sl[0]] - X[i, 0]; dy = axes[1][sl[1]] - X[i, 1]; dz = axes[2][sl[2]] - X[i, 2]
        d = np.sqrt(dx[:, None, None] ** 2 + dy[None, :, None] ** 2 + dz[None, None, :] ** 2)
        v = pm.rho_grid(E[i], d); v[d > R_GRID_CUT] = 0.0
        T[sl[0], sl[1], sl[2]] += v
        if inside[i]:
            A[sl[0], sl[1], sl[2]] += v
    W = A / np.maximum(T, 1e-300)
    edge = max(W[0].max(), W[-1].max(), W[:, 0].max(), W[:, -1].max(), W[:, :, 0].max(), W[:, :, -1].max())
    if edge >= 0.5:
        raise RuntimeError(f'{target}: Hirshfeld weight reaches {edge:.2f} on the box edge -- the surface is open, '
                           f'i.e. the cluster is missing neighbours; increase --cutoff')
    V, Fc = marching(W, lo, spacing)
    # Newton projection of the linearly interpolated vertices onto the exact w = 0.5 surface
    for it in range(4):
        v = pm.at_points(V, X, E, inside)
        w, gw, _ = weight_and_derivs(v)
        V = V - ((w - 0.5) / np.sum(gw * gw, 1))[:, None] * gw
    v = pm.at_points(V, X, E, inside, want_hess=True)
    w, gw, Hw = weight_and_derivs(v, want_hess=True)
    resid = float(np.abs(w - 0.5).max())
    k1, k2, nrm = curvatures(gw, Hw)
    S = (2 / np.pi) * np.arctan2(k1 + k2, k1 - k2)
    C = (2 / np.pi) * np.log(np.sqrt((k1 ** 2 + k2 ** 2) / 2))
    # geometry of the mesh: triangle areas, vertex areas (1/3 of adjacent triangles), enclosed volume
    p0, p1, p2 = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
    cr = np.cross(p1 - p0, p2 - p0); tri_area = 0.5 * np.linalg.norm(cr, axis=1)
    # orient faces outward (marching_cubes orientation depends on the gradient direction)
    outward = np.sum(cr * (-gw[Fc].mean(1)), 1)
    if np.mean(outward < 0) > 0.5:
        Fc = Fc[:, ::-1].copy()
        p0, p1, p2 = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
    n_inverted = int(np.sum(np.sum(np.cross(p1 - p0, p2 - p0) * (-gw[Fc].mean(1)), 1) < 0))
    vol = float(np.sum(np.einsum('ij,ij->i', p0, np.cross(p1, p2))) / 6.0)
    varea = np.zeros(len(V)); np.add.at(varea, Fc.ravel(), np.repeat(tri_area / 3, 3))
    area = float(tri_area.sum())
    edges = np.unique(np.sort(np.vstack([Fc[:, [0, 1]], Fc[:, [1, 2]], Fc[:, [0, 2]]]), axis=1), axis=0)
    euler = len(V) - len(edges) + len(Fc)
    # nearest nuclei inside / outside
    from scipy.spatial import cKDTree
    idx_in = np.where(inside)[0]; idx_out = np.where(~inside)[0]
    di, ii = cKDTree(X[idx_in]).query(V); de, jj = cKDTree(X[idx_out]).query(V)
    ai = idx_in[ii]; ae = idx_out[jj]
    ri = np.array([BONDI[E[k]] for k in ai]); re_ = np.array([BONDI[E[k]] for k in ae])
    dnorm = (di - ri) / ri + (de - re_) / re_
    ei = np.array([E[k] for k in ai]); ee = np.array([E[k] for k in ae])
    mol = np.array([a['mol'] for a in atoms])
    # ---- contacts by element pair
    rows = []
    for (x, y) in sorted(set(zip(ei, ee))):
        s = (ei == x) & (ee == y)
        rows.append({'contact': f'{x}...{y}', 'area_A2': float(varea[s].sum()), 'pct': 100 * float(varea[s].sum()) / area,
                     'min_dnorm': float(dnorm[s].min()), 'min_di_plus_de_A': float((di + de)[s].min())})
    recip = {}
    for r in rows:
        x, y = r['contact'].split('...')
        key = f'{x}...{y}' if x == y else '/'.join(sorted([f'{x}...{y}', f'{y}...{x}'], key=lambda s: (s[0] != 'H', s)))
        acc = recip.setdefault(key, {'contact': key, 'area_A2': 0.0, 'pct': 0.0, 'min_dnorm': 9.9})
        acc['area_A2'] += r['area_A2']; acc['pct'] += r['pct']; acc['min_dnorm'] = min(acc['min_dnorm'], r['min_dnorm'])
    recip = sorted(recip.values(), key=lambda r: -r['pct'])
    with open(out / f'{target}_contacts.csv', 'w') as f:
        f.write('contact,kind,area_A2,pct_of_surface,min_dnorm,min_di_plus_de_A\n')
        for r in recip:
            f.write(f"{r['contact']},reciprocal-combined,{r['area_A2']:.3f},{r['pct']:.2f},{r['min_dnorm']:.4f},\n")
        for r in sorted(rows, key=lambda r: -r['pct']):
            f.write(f"{r['contact']},inside...outside,{r['area_A2']:.3f},{r['pct']:.2f},{r['min_dnorm']:.4f},{r['min_di_plus_de_A']:.3f}\n")
    # ---- neighbour patches
    neigh = []
    for m in np.unique(mol[ae]):
        s = mol[ae] == m
        k = np.argmin(dnorm[s]); vi = np.where(s)[0][k]
        mm = meta['molecules'][m]
        neigh.append({'molecule': int(m), 'kind': mm['kind'], 'symcode': mm['symcode'], 'central': mm['central'],
                      'area_A2': float(varea[s].sum()), 'pct': 100 * float(varea[s].sum()) / area,
                      'min_dnorm': float(dnorm[vi]), 'closest_pair': f'{lab[ai[vi]]}...{lab[ae[vi]]}',
                      'closest_pair_dist_A': float(np.linalg.norm(X[ai[vi]] - X[ae[vi]]))})
    neigh.sort(key=lambda r: -r['area_A2'])
    with open(out / f'{target}_neighbours.csv', 'w') as f:
        f.write('molecule,kind,symcode,area_A2,pct_of_surface,min_dnorm,closest_pair,closest_pair_dist_A\n')
        for r in neigh:
            f.write(f"{r['molecule']},{r['kind']},{r['symcode']},{r['area_A2']:.3f},{r['pct']:.2f},{r['min_dnorm']:.4f},"
                    f"{r['closest_pair']},{r['closest_pair_dist_A']:.3f}\n")
    # ---- red spots: local minima of d_norm below zero, one per (inside atom, outside atom) pair
    nb = [[] for _ in range(len(V))]
    for a_, b_ in edges:
        nb[a_].append(b_); nb[b_].append(a_)
    spots = {}
    for vi in np.where(dnorm < 0)[0]:
        if all(dnorm[vi] <= dnorm[j] for j in nb[vi]):
            key = (int(ai[vi]), int(ae[vi]))
            if key not in spots or dnorm[vi] < spots[key]['dnorm']:
                spots[key] = {'dnorm': float(dnorm[vi]), 'inside': lab[ai[vi]], 'outside': lab[ae[vi]],
                              'dist_A': float(np.linalg.norm(X[ai[vi]] - X[ae[vi]])),
                              'di_A': float(di[vi]), 'de_A': float(de[vi]), 'xyz': V[vi].round(4).tolist()}
    spots = sorted(spots.values(), key=lambda s: s['dnorm'])
    with open(out / f'{target}_red_spots.csv', 'w') as f:
        f.write('inside_atom,outside_atom,interatomic_dist_A,min_dnorm,d_i_A,d_e_A,x,y,z\n')
        for s in spots:
            f.write(f"{s['inside']},{s['outside']},{s['dist_A']:.3f},{s['dnorm']:.4f},{s['di_A']:.3f},{s['de_A']:.3f},"
                    + ','.join(f'{c:.4f}' for c in s['xyz']) + '\n')
    # ---- fingerprint (total and per element pair), area-weighted
    H_tot, _, _ = np.histogram2d(di, de, bins=[FP_BINS, FP_BINS], weights=varea)
    fp = {'bins_A': FP_BINS, 'total_A2': H_tot}
    for x, y in sorted(set(zip(ei, ee))):
        s = (ei == x) & (ee == y)
        fp[f'{x}...{y}'] = np.histogram2d(di[s], de[s], bins=[FP_BINS, FP_BINS], weights=varea[s])[0]
    np.savez_compressed(out / f'{target}_fingerprint.npz', **fp)
    # ---- vertex data
    np.savez_compressed(out / f'{target}_surface.npz', vertices=V, faces=Fc, normals=nrm, d_i=di, d_e=de, d_norm=dnorm,
                        shape_index=S, curvedness=C, k1=k1, k2=k2, vertex_area=varea, inside_atom=ai, outside_atom=ae,
                        atom_labels=np.array(lab), atom_xyz=X, atom_el=np.array(E), inside_mask=inside)
    with gzip.open(out / f'{target}_vertices.csv.gz', 'wt') as f:
        f.write('x,y,z,d_i,d_e,d_norm,shape_index,curvedness,vertex_area_A2,inside_atom,outside_atom\n')
        for k in range(len(V)):
            f.write(f'{V[k,0]:.4f},{V[k,1]:.4f},{V[k,2]:.4f},{di[k]:.4f},{de[k]:.4f},{dnorm[k]:.4f},{S[k]:.4f},{C[k]:.4f},'
                    f'{varea[k]:.5f},{lab[ai[k]]},{lab[ae[k]]}\n')
    glob = (36 * np.pi * vol ** 2) ** (1 / 3) / area
    summ = {'target': target, 'n_inside_atoms': int(inside.sum()), 'grid_spacing_A': spacing, 'grid_shape': shape,
            'n_vertices': int(len(V)), 'n_faces': int(len(Fc)), 'n_faces_inverted': n_inverted,
            'euler_characteristic': int(euler),
            'max_abs_w_residual': resid, 'area_A2': area, 'volume_A3': vol, 'globularity': float(glob),
            'd_norm_min': float(dnorm.min()), 'd_norm_max': float(dnorm.max()), 'd_i_range_A': [float(di.min()), float(di.max())],
            'd_e_range_A': [float(de.min()), float(de.max())], 'shape_index_range': [float(S.min()), float(S.max())],
            'curvedness_range': [float(C.min()), float(C.max())],
            'contacts_reciprocal': recip, 'contacts_directional': sorted(rows, key=lambda r: -r['pct']),
            'neighbours': neigh, 'n_neighbours_in_contact': len(neigh), 'red_spots': spots,
            'wall_s': round(time.time() - t0, 1)}
    print(f"[hs] {target}: {len(V)} vertices, area {area:.2f} A^2, volume {vol:.2f} A^3, G {glob:.3f}, "
          f"d_norm {dnorm.min():.3f}..{dnorm.max():.3f}, Euler {euler}, |w-0.5|max {resid:.1e}, "
          f"{len(neigh)} neighbours, {len(spots)} red spots ({time.time()-t0:.0f} s)", flush=True)
    print('      contacts: ' + ', '.join(f"{r['contact']} {r['pct']:.1f}%" for r in recip), flush=True)
    if png:
        quicklook_png(out, target, fp, recip)
    return summ


def quicklook_png(out, target, fp, recip):
    """Quick-look fingerprint plots (total + the four largest reciprocal contacts); publication figures come later."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    b = fp['bins_A']; tot = fp['total_A2']
    keys = [r['contact'] for r in recip[:4]]
    fig, axs = plt.subplots(1, 1 + len(keys), figsize=(3.2 * (1 + len(keys)), 3.2), constrained_layout=True)
    def show(ax, Hm, title):
        Hm = np.ma.masked_where(Hm <= 0, Hm)
        ax.pcolormesh(b, b, Hm.T, cmap='jet', shading='flat'); ax.set_xlim(0.6, 2.6); ax.set_ylim(0.6, 2.6)
        ax.set_aspect('equal'); ax.set_title(title, fontsize=9); ax.set_xlabel('d_i (A)'); ax.set_ylabel('d_e (A)')
    show(axs[0], tot, f'{target}: all contacts')
    for ax, k in zip(axs[1:], keys):
        parts = [p for p in k.split('/')]
        Hm = sum(fp[p] for p in parts if p in fp)
        ax.pcolormesh(b, b, np.ma.masked_where(tot.T <= 0, tot.T), cmap='Greys', vmin=0, vmax=tot.max() * 4, shading='flat')
        show(ax, Hm, f"{k} {[r['pct'] for r in recip if r['contact'] == k][0]:.1f}%")
    fig.savefig(out / f'{target}_fingerprint_quicklook.png', dpi=150); plt.close(fig)


# ============================================================ Multiwfn cross-check (x86 VM)
def mwfn_crosscheck(meta, out, targets, spacing_bohr=0.2):
    mw = pathlib.Path(os.environ.get('MWFN_BIN') or ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI')
    if not mw.exists():
        raise SystemExit(f'Multiwfn not found at {mw} (x86 only; run this step on the VM)')
    nt = os.environ.get('MWFN_THREADS', '2')
    res = {}
    for target in targets:
        roles = {'cation': ('cation',), 'anion': ('anion',), 'pair': ('cation', 'anion')}[target]
        ids = [i + 1 for i, a in enumerate(meta['atoms']) if a['role'] in roles]
        sel = ','.join(str(i) for i in ids)
        # 12 surface analysis; 1 -> 5 Hirshfeld surface of the fragment; 2 -> 22 map d_norm; 3 grid spacing (Bohr);
        # 0 start; post-processing 20 -> 3 element-pair contact areas, 0 fingerprint analysis, 5 export d_i/d_e,
        # -1 -1 back to the post-processing menu, -1 -1 back to the main menu, q quit.
        seq = ['12', '1', '5', sel, '2', '22', '3', str(spacing_bohr), '0', '20', '3', '0', '5', '-1', '-1', '-1', '-1', 'q']
        wd = out / f'mwfn_{target}'; wd.mkdir(exist_ok=True)
        (wd / 'mwfn.in').write_text('\n'.join(seq) + '\n')
        # OMP_THREAD_LIMIT: Multiwfn's nested OpenMP regions otherwise use ~2x the -nt thread count
        env = dict(os.environ, Multiwfnpath=str(mw.parent), KMP_STACKSIZE='200M', OMP_THREAD_LIMIT=nt)
        t = time.time()
        with open(wd / 'mwfn.in') as fi, open(wd / 'mwfn.out', 'w') as fo:
            rc = subprocess.run([str(mw), str(out / 'cluster.xyz'), '-nt', nt, '-silent'], stdin=fi, stdout=fo,
                                stderr=subprocess.STDOUT, cwd=wd, env=env, timeout=3600,
                                preexec_fn=_unlimit_stack).returncode
        txt = (wd / 'mwfn.out').read_text(errors='replace')
        g = lambda pat: (lambda m: float(m.group(1)) if m else None)(re.search(pat, txt))
        r = {'exit': rc, 'wall_s': round(time.time() - t, 1), 'menu_sequence': seq,
             'version': (re.search(r'Version (\S+)', txt) or [None, None])[1],
             'area_A2': g(r'Isosurface area:.*?\(\s*([\d.]+) Angstrom\^2'),
             'volume_A3': g(r'Volume enclosed by the isosurface:.*?\(\s*([\d.]+) Angstrom\^3'),
             'd_norm_min': g(r'Global surface minimum:\s*(-?[\d.]+)'), 'd_norm_max': g(r'Global surface maximum:\s*(-?[\d.]+)')}
        m = re.search(r'do not distinguish inside and outside elements\n(.*?)\n\s*Area of total contact surface', txt, re.S)
        r['contacts_reciprocal'] = [{'contact': p[0], 'area_A2': float(p[1]), 'pct': float(p[2])}
                                    for p in (l.split() for l in (m.group(1).splitlines() if m else [])) if len(p) == 3]
        m = re.search(r'Inside element, outside element.*?\n(.*?)\n\s*The same as above', txt, re.S)
        r['contacts_directional'] = [{'contact': p[0], 'area_A2': float(p[1]), 'pct': float(p[2])}
                                     for p in (l.split() for l in (m.group(1).splitlines() if m else [])) if len(p) == 3]
        if r['area_A2'] is None or not r['contacts_reciprocal']:
            raise RuntimeError(f'Multiwfn Hirshfeld surface for {target} failed (exit {rc}); see {wd}/mwfn.out')
        save_json(out / f'mwfn_{target}.json', r)
        res[target] = r
        print(f"[hs] Multiwfn {target}: area {r['area_A2']} A^2, volume {r['volume_A3']} A^3, d_norm min {r['d_norm_min']} "
              f"({r['wall_s']} s)", flush=True)
    return res


def _unlimit_stack():
    import resource
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, OSError):
        pass


def compare(out, targets):
    summ = json.loads((out / 'hirshfeld_summary.json').read_text())
    comp = {}
    for t in targets:
        p = out / f'mwfn_{t}.json'
        if not p.exists() or t not in summ['targets']:
            continue
        m = json.loads(p.read_text()); s = summ['targets'][t]
        # Multiwfn writes element pairs as "H-O/O-H"; ours are "H...O/O...H" -- compare on the sorted element set
        key = lambda c: tuple(sorted(re.split(r'\.\.\.|-|/', c)[:2]))
        py = {key(r['contact']): r['pct'] for r in s['contacts_reciprocal']}
        mw = {key(r['contact']): r['pct'] for r in m['contacts_reciprocal']}
        comp[t] = {'area_A2': [s['area_A2'], m['area_A2']], 'volume_A3': [s['volume_A3'], m['volume_A3']],
                   'd_norm_min': [s['d_norm_min'], m['d_norm_min']],
                   'contact_pct_python_vs_multiwfn': {'...'.join(k): [round(py.get(k, 0), 2), round(mw.get(k, 0), 2)]
                                                      for k in sorted(set(py) | set(mw))}}
    save_json(out / 'hirshfeld_compare.json', comp)
    print(json.dumps(comp, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('steps', nargs='*', default=['build', 'surface'])
    ap.add_argument('--out', default=str(OUT_DEFAULT))
    ap.add_argument('--cutoff', type=float, default=6.0,
                    help='neighbour ions with any atom within this of the central unit (A); 4, 6 and 8 A give identical surfaces')
    ap.add_argument('--spacing', type=float, default=0.12, help='grid spacing for marching cubes (A)')
    ap.add_argument('--targets', default='cation,anion')
    ap.add_argument('--no-normalise-h', action='store_true')
    ap.add_argument('--png', action='store_true')
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    targets = a.targets.split(',')
    meta = None
    if 'build' in a.steps:
        meta = build_cluster(out, a.cutoff, normalise_h=not a.no_normalise_h)
    if meta is None:
        meta = json.loads((out / 'cluster.json').read_text())
    if 'surface' in a.steps:
        r, ln = atomic_density_table(sorted({x['el'] for x in meta['atoms']}), out)
        pm = Promolecule(r, ln)
        summ = {'method': {'promolecule': f'spherically averaged ROHF/{ATOM_BASIS} free atoms (PySCF)',
                           'vdw_radii_A': BONDI, 'xh_normalised': meta['xh_normalised'], 'cluster_cutoff_A': meta['cutoff_A'],
                           'grid_spacing_A': a.spacing, 'grid_neighbour_cut_A': R_GRID_CUT, 'point_neighbour_cut_A': R_POINT_CUT,
                           'shape_index': 'S = (2/pi) atan[(k1+k2)/(k1-k2)], k>0 convex; +1 cap (blue), -1 cup (red)',
                           'curvedness': 'C = (2/pi) ln sqrt[(k1^2+k2^2)/2], k in 1/A',
                           'fingerprint_bins_A': [float(FP_BINS[0]), float(FP_BINS[-1]), 0.01]},
                'cell_volume_A3': meta['cell_volume_A3'], 'Z': meta['Z'], 'targets': {}}
        prev = out / 'hirshfeld_summary.json'
        if prev.exists() and 'build' not in a.steps:     # a later '--targets pair' run must not drop cation/anion
            old = json.loads(prev.read_text())
            if old.get('method') == json.loads(json.dumps(summ['method'])):
                summ['targets'] = old.get('targets', {})
        for t in targets:
            summ['targets'][t] = surface_for_target(t, meta, pm, out, a.spacing, png=a.png)
        if {'cation', 'anion'} <= set(summ['targets']):
            vsum = summ['targets']['cation']['volume_A3'] + summ['targets']['anion']['volume_A3']
            summ['formula_unit_volume_A3'] = vsum
            summ['packing_fraction_of_V_cell_over_Z'] = vsum / (meta['cell_volume_A3'] / meta['Z'])
            print(f"[hs] V(cation)+V(anion) = {vsum:.2f} A^3 = {100*summ['packing_fraction_of_V_cell_over_Z']:.1f} % of V_cell/Z", flush=True)
        save_json(out / 'hirshfeld_summary.json', summ)
    if 'mwfn' in a.steps:
        mwfn_crosscheck(meta, out, targets)
    if 'compare' in a.steps:
        compare(out, targets)


if __name__ == '__main__':
    main()
