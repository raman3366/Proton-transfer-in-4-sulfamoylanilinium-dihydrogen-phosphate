"""Helpers shared by figure group B of the SAPH paper (main Figs. 5-10 and their SI figures).

Each figure script stays runnable on its own, from the repository root:

    nice -n 10 tools/paper-env/bin/python paper/main/figures/scripts/<stem>.py
    nice -n 10 tools/paper-env/bin/python paper/si/figures/scripts/<stem>.py

What lives here
  * thread and priority limits for the shared Pi (<= 2 threads for BLAS/OpenMP/llvmpipe, niceness >= 10),
    set BEFORE numpy / VTK are imported;
  * facts: F(key) reads paper/main/data/facts.json, fmt(key) prints it at its stated precision with a true minus
    sign; every key a figure prints is recorded and written next to the figure as <stem>.facts.json;
  * rendering parameters that are not results (isovalues, colour ranges, broadening) come from
    paper/main/figures/figure_params.json through PARAM(key), recorded the same way;
  * render cache (pickled figkit Render objects under paper/main/figures/work/r3d_cache, git-ignored, keyed by name, parameters and the
    mtimes/sizes of the input files) so that re-running a script after a layout change does not re-render;
  * save(): figkit QA-checked 600 dpi PNG + vector PDF.
Raw inputs (cubes, Multiwfn grids, Hirshfeld surfaces) are read from the archive (/root/cqc-sfa-archive/vm/cqc-sfa)
and never written; outputs go only to paper/main/figures and paper/si/figures.
"""
from __future__ import annotations

import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'LP_NUM_THREADS',
           'VTK_SMP_MAX_THREADS'):
    os.environ.setdefault(_v, '2')
os.environ.setdefault('VTK_SMP_IMPLEMENTATION_TYPE', 'Sequential')
try:                                   # the Pi is shared: never run a render above niceness 10
    _cur = os.nice(0)
    if _cur < 10:
        os.nice(10 - _cur)
except OSError:
    pass

import hashlib  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / 'paper' / '_shared'))

from figkit import compose, geometry, qa, style  # noqa: E402

FIG_MAIN = REPO / 'paper' / 'main' / 'figures'
FIG_SI = REPO / 'paper' / 'si' / 'figures'
FACTS_PATH = REPO / 'paper' / 'main' / 'data' / 'facts.json'
PARAMS_PATH = FIG_MAIN / 'figure_params.json'
CACHE = FIG_MAIN / 'work' / 'r3d_cache'           # git-ignored (paper/**/figures/work/); safe to delete
ARCHIVE = style.DATA                                   # /root/cqc-sfa-archive/vm/cqc-sfa (read only)
RES_RAW = ARCHIVE / 'molecules' / 'sfa' / 'results'
RES = REPO / 'molecules' / 'sfa' / 'results'
DOCK = REPO / 'dock'
MINUS = '−'
T0 = time.time()

_FACTS = None
_PARAMS = None
USED = {}                                              # key -> printed value (facts and params)


def log(msg):
    print(f'[{time.time() - T0:7.1f} s] {msg}', flush=True)


# ------------------------------------------------------------------------------------------------------------
# facts and parameters
# ------------------------------------------------------------------------------------------------------------
def facts():
    global _FACTS
    if _FACTS is None:
        _FACTS = json.loads(FACTS_PATH.read_text())
    return _FACTS


def params():
    global _PARAMS
    if _PARAMS is None:
        _PARAMS = json.loads(PARAMS_PATH.read_text()) if PARAMS_PATH.exists() else {}
    return _PARAMS


def F(key):
    """Raw value of a fact (KeyError if the fact does not exist: a figure never prints an unsourced number)."""
    v = facts()[key]['value']
    USED.setdefault(key, v)
    return v


def PARAM(key):
    """A rendering parameter (not a result) from figure_params.json."""
    v = params()[key]['value']
    USED.setdefault('param:' + key, v)
    return v


def num(v, nd, signed=False):
    """Format a number with nd decimals and a true minus sign ('-0.00' -> '0.00')."""
    s = f'{float(v):+.{nd}f}' if signed else f'{float(v):.{nd}f}'
    if s.lstrip('+-').strip('0.') == '':
        s = s.lstrip('+-')
        if signed:
            s = s
    return s.replace('-', MINUS)


def fmt(key, nd=None, signed=False):
    """A fact printed at `nd` decimals (default: the fact's own precision)."""
    f = facts()[key]
    if nd is None:
        nd = f.get('precision')
        nd = 2 if nd is None else int(nd)
    s = num(f['value'], nd, signed)
    USED[key] = s
    return s


def write_provenance(path_stem, extra=None):
    """<stem>.facts.json next to the figure: every fact / parameter key the figure printed or used."""
    p = Path(str(path_stem) + '.facts.json')
    out = {'script': str(Path(sys.argv[0]).resolve().relative_to(REPO)) if sys.argv and sys.argv[0] else None,
           'facts_file': str(FACTS_PATH.relative_to(REPO)), 'params_file': str(PARAMS_PATH.relative_to(REPO)),
           'used': USED}
    if extra:
        out.update(extra)
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str) + '\n')
    return p


# ------------------------------------------------------------------------------------------------------------
# cache
# ------------------------------------------------------------------------------------------------------------
def _fingerprint(files):
    parts = []
    for f in files or []:
        f = Path(f)
        st = f.stat()
        parts.append(f'{f}|{st.st_size}|{int(st.st_mtime)}')
    return parts


def cached(name, fn, files=(), extra=''):
    """fn() cached on disk by name + input-file fingerprints + extra (any repr-able parameters)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(json.dumps([name, _fingerprint(files), repr(extra)]).encode()).hexdigest()[:16]
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in name)[:60]
    p = CACHE / f'{safe}__{key}.pkl'
    if p.exists():
        try:
            with open(p, 'rb') as fh:
                return pickle.load(fh)
        except Exception:                                 # a stale or truncated cache file is simply rebuilt
            pass
    t = time.time()
    out = fn()
    with open(p, 'wb') as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    log(f'cached {name} ({time.time() - t:.1f} s)')
    return out


# ------------------------------------------------------------------------------------------------------------
# saving
# ------------------------------------------------------------------------------------------------------------
def save(fig, stem, si=False, strict=True, width_mm=None, extra=None):
    """QA-checked 600 dpi PNG + PDF in paper/main/figures (or paper/si/figures) and the provenance sidecar."""
    out_dir = FIG_SI if si else FIG_MAIN
    if os.environ.get('FIGB_DRAFT'):              # layout iterations: write the QA issues, do not stop
        strict = False
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / stem
    if hasattr(fig, '_figkit_grid'):
        paths = compose.save(fig, path, strict=strict)
    else:
        paths = qa.save(fig, path, strict=strict, max_width_mm=width_mm, glyphs=True)
    prov = write_provenance(path, extra)
    for p in list(paths) + [prov]:
        log(f'wrote {p.relative_to(REPO)}')
    return paths


# ------------------------------------------------------------------------------------------------------------
# molecules
# ------------------------------------------------------------------------------------------------------------
def pair_labels():
    """Clean CIF labels of the 27-atom ion pair (cation atoms, then the anion at x, y, z - 1)."""
    return geometry.read_case('sfa_pair_pcm_d3')['labels']


def index_of(labels, lab):
    return list(labels).index(lab)


# ------------------------------------------------------------------------------------------------------------
# 3D: the ion pair, its camera, cached renders, NCI/IGMH colours
# ------------------------------------------------------------------------------------------------------------
_MOL = {}


def molecule(case='sfa_pair_pcm_d3'):
    """figkit Molecule of a case (opt.xyz, clean CIF labels, fragments, the case-JSON hydrogen bonds)."""
    from figkit import render3d as r3
    if case not in _MOL:
        _MOL[case] = r3.load_molecule(case)
    return _MOL[case]


def pair_view():
    """ONE camera for every render of the IEF-PCM ion pair in figure group B (Figs. 5, 7, 8): ring_view on the
    R2,2(8) ring (N1, H1A, O5, P1, O4, H6, C6, C1) about 25-55 deg off its normal, so that the discs of both of
    its hydrogen bonds (N1-H1A...O5, C6-H6...O4) are seen obliquely and no heavy atom (or H1A, H6) hides behind
    another; landscape, cation left.  Unframed (fit per figure)."""
    mol = molecule('sfa_pair_pcm_d3')
    L = mol.labels
    ring = [L.index(x) for x in ('N1', 'H1A', 'O5', 'P1', 'O4', 'H6', 'C6', 'C1')]
    key = [i for i, s in enumerate(mol.symbols) if s != 'H'] + [L.index('H1A'), L.index('H6')]
    return cached('pair_view_pcm_v3', lambda: ring_view(mol, ring, mol.hbonds, key_atoms=key),
                  files=[RES / 'sfa_pair_pcm_d3' / 'opt.xyz'])


def render_cached(name, scene_fn, view, px_per_A, files=(), extra='', **kw):
    """figkit render3d.render of scene_fn() with `view`, cached on disk (the scene is only built when needed)."""
    from figkit import render3d as r3
    key = (repr(np.round(view.focal, 5).tolist()), repr(np.round(view.direction, 6).tolist()),
           repr(np.round(view.up, 6).tolist()), view.width, view.height, round(px_per_A, 4), extra,
           sorted(kw.items()))
    return cached('render_' + name, lambda: r3.render(scene_fn(), view, px_per_A=px_per_A, **kw), files=files,
                  extra=key)


def nci_cmap():
    from figkit import render3d as r3
    return r3.NCI_CMAP


def two_slope_scale(vmin, vmax, label=r'sign($\lambda_2$)$\rho$ (a.u.)', step=0.02, fmt='{:.2f}'):
    """NCI/IGMH colour scale with green exactly at 0 (two linear slopes), for a render3d/compose colour bar."""
    from matplotlib import colors as mcolors
    from figkit import render3d as r3
    norm = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    ticks = sorted({round(t, 4) for t in np.arange(np.ceil(vmin / step) * step, vmax + 1e-12, step)} | {0.0})
    return r3.ColourScale(r3.NCI_CMAP, norm, label, ticks, fmt)


def recolour_two_slope(scene, vmin, vmax, field='sl2r'):
    """Re-map the sign(lambda2)rho colours of an nci_scene to a two-slope norm (green at 0)."""
    scale = two_slope_scale(vmin, vmax)
    for m, k in scene.items:
        if field in m.point_data:
            m.point_data['t'] = np.clip(scale.norm(np.asarray(m.point_data[field], float)), 0, 1)
            k['cmap'] = scale.cmap
    scene.colour = scale
    return scene


def _mesh_neighbours(surf):
    from scipy import sparse
    f = np.asarray(surf.faces).reshape(-1, 4)
    assert (f[:, 0] == 3).all(), 'triangles expected'
    f = f[:, 1:]
    i = np.concatenate([f[:, 0], f[:, 1], f[:, 2], f[:, 1], f[:, 2], f[:, 0]])
    j = np.concatenate([f[:, 1], f[:, 2], f[:, 0], f[:, 0], f[:, 1], f[:, 2]])
    n = surf.n_points
    A = sparse.csr_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    A = (A + sparse.identity(n, format='csr')).tocsr()
    A.sum_duplicates()
    return A.indptr, A.indices


def _vertex_areas(surf):
    f = np.asarray(surf.faces).reshape(-1, 4)[:, 1:]
    P = np.asarray(surf.points, float)
    a = 0.5 * np.linalg.norm(np.cross(P[f[:, 1]] - P[f[:, 0]], P[f[:, 2]] - P[f[:, 0]]), axis=1)
    A = np.zeros(surf.n_points)
    for k in range(3):
        np.add.at(A, f[:, k], a / 3)
    return A


def smooth_sl2r(scene, vmin, vmax, outlier_passes=3, mean_passes=3, taubin_iter=20, field='sl2r'):
    """Clean sign(lambda2)rho on an NCI/IGMH isosurface without posterising it (adapted from the 3ATB paper's
    _r3d3atb.smooth_nci, cqc-dft): (1) a vertex whose sign differs from ALL its 1-ring neighbours (speckle from
    trilinear interpolation of the signed cube across a lambda2 sign change) takes their median; (2) mean_passes
    of area-weighted 1-ring mean smoothing; (3) Taubin smoothing of the vertex positions (serrated rims where a
    thin lens meets the grid).  Colours are then re-mapped with the two-slope norm.  Returns stats."""
    surf, kw = next((m, k) for m, k in scene.items if field in m.point_data)
    v0 = np.asarray(surf.point_data[field], float)
    move = 0.0
    if taubin_iter:
        P0 = np.asarray(surf.points, float).copy()
        sm = surf.smooth_taubin(n_iter=taubin_iter, pass_band=0.1, boundary_smoothing=False,
                                feature_smoothing=False)
        surf.points = np.asarray(sm.points, float)
        move = float(np.linalg.norm(np.asarray(surf.points) - P0, axis=1).max())
    indptr, indices = _mesh_neighbours(surf)
    n = surf.n_points
    deg = np.diff(indptr)
    rows = np.repeat(np.arange(n), deg)
    nb = indices != rows
    v = v0.copy()
    fixed = 0
    for _ in range(outlier_passes):
        sv = np.sign(v)
        same = np.zeros(n)
        np.add.at(same, rows[nb], (np.sign(v[indices[nb]]) == sv[rows[nb]]).astype(float))
        out = np.nonzero((same == 0) & (np.bincount(rows[nb], minlength=n) > 0))[0]
        if not len(out):
            break
        for i in out:
            js = indices[indptr[i]:indptr[i + 1]]
            v[i] = np.median(v[js[js != i]])
        fixed += len(out)
    A = _vertex_areas(surf)
    for _ in range(mean_passes):
        num = np.zeros(n)
        den = np.zeros(n)
        np.add.at(num, rows, A[indices] * v[indices])
        np.add.at(den, rows, A[indices])
        v = num / den
    surf.point_data[field] = v
    scale = two_slope_scale(vmin, vmax)
    surf.point_data['t'] = np.clip(scale.norm(v), 0, 1)
    kw['cmap'] = scale.cmap
    scene.colour = scale
    return dict(n_vertices=int(n), sign_outliers_fixed=int(fixed), max_vertex_move_A=round(move, 4),
                sl2r_min=float(v.min()), sl2r_max=float(v.max()))


def ring_view(mol, ring, hbonds, tilts=range(25, 56, 5), phis=range(0, 360, 10), key_atoms=None):
    """Camera that looks at a hydrogen-bonded ring (e.g. the R2,2(8) ring) about `tilt` degrees off its normal,
    so the H-bond discs (normal H->A, in the ring plane) are seen obliquely, while no key atom hides behind
    another one.  Score per candidate: sum over hbonds of min(|n.d| / 0.5, 1) - 1.5 x occluded key atoms
    (a nearer atom ball covering > 40 % of the radius) - 3 x the hidden fraction of each 0.6 A H-bond disc
    - 0.5 x (tilt - 40)^2 / 225.  Landscape, first fragment
    on the left.  Returns (View, info)."""
    from figkit import render3d as r3
    X = mol.coords
    ring = np.asarray(ring)
    c = X[ring].mean(0)
    w, V = np.linalg.eigh(np.cov((X[ring] - c).T))
    nrm = V[:, 0]
    rad = np.array([style.BALL_RADII.get(s, 0.35) for s in mol.symbols])
    key = np.flatnonzero(mol.heavy) if key_atoms is None else np.asarray(key_atoms)
    best = None
    e1 = V[:, 2]
    for t in tilts:
        for phi in phis:
            ax_ = np.cos(np.radians(phi)) * e1 + np.sin(np.radians(phi)) * np.cross(nrm, e1)
            R = r3._rot(ax_, t)
            for sgn in (1, -1):
                d = R @ (sgn * nrm)
                v = r3.View(c, d, e1)
                s = v.screen(X)
                occ = 0
                for i in key:
                    dd = np.hypot(s[:, 0] - s[i, 0], s[:, 1] - s[i, 1])
                    front = (s[:, 2] > s[i, 2] + 0.1) & (dd < rad + 0.6 * rad[i])
                    front[i] = False
                    occ += int(front.any())
                vis = sum(min(abs(r3._unit(X[a] - X[h]) @ d) / 0.5, 1.0) for _, h, a in hbonds)
                hid = sum(r3._disc_hidden_fraction(v, mol, h, a, 0.6) for _, h, a in hbonds)
                score = vis - 1.5 * occ - 3.0 * hid - 0.5 * (t - 40) ** 2 / 225
                if best is None or score > best[0] + 1e-9:
                    best = (score, v, t, phi, sgn, occ, vis, hid)
    v = r3._orient(best[1], mol, None, True)
    return v, dict(score=round(best[0], 3), tilt=best[2], phi=best[3], side=best[4], occluded=best[5],
                   disc_visibility=round(best[6], 3), disc_hidden=round(best[7], 3))


# fragment / species colours used in every plot of figure group B (colour-blind-safe blue / orange / green)
# sulfonamidate was '#8E44AD' until the QA pass of 28 Sep 2026: under deuteranopia it was indistinguishable from the
# cation blue (OKLab dE x100 = 1.5, Machado 2009 simulation); Okabe-Ito reddish purple '#CC79A7' gives >= 10.9
# against both the cation and sulfanilamide colours for protan/deutan/tritan.  Used only in Fig. 11.
FRAG_COLOURS = {'cation': '#2B6CB0', 'anion': '#D9730D', 'pair': '#222222', 'sulfanilamide': '#2E8B57',
                'sulfonamidate': '#CC79A7'}


def panel_letter(fig, ax, letter, title=None, dx_mm=0.0, dy_mm=0.6):
    """(a) + optional title above the top-left corner of an ordinary matplotlib axes (house style)."""
    from matplotlib.transforms import ScaledTranslation
    MMi = style.MM
    bb = ax.get_position()
    x = bb.x0 + dx_mm * MMi / fig.get_figwidth()
    y = bb.y1 + dy_mm * MMi / fig.get_figheight()
    return compose._header(fig, x, y + (style.PANEL_PT / 72 * 1.25) / fig.get_figheight(), f'({letter})', title)
