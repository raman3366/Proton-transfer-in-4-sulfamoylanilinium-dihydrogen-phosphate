"""Shared helpers for the SAPH figure scripts of group A (structure and proton story).

    cd /root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate
    nice -n 10 env OMP_NUM_THREADS=2 /root/workspace/cqc-dft/tools/paper-env/bin/python \
        paper/main/figures/scripts/<stem>.py

What lives here
    facts        F(key) -> value; txt(key, nd=None, sign=False) -> printed string (precision from facts.json,
                 half-up rounding, U+2212 minus); every key read is recorded and written next to the figure
                 (figures/build/<stem>_facts_used.json) so the numbers audit can trace each printed number
    geometries   xray_pair(), case_pair(case), scan_point(name), periodic_pair(name): the 27-atom ion pair in the
                 CIF atom order (PAIR_LABELS), in Angstrom; the periodic pair is formula unit 1 of the relaxed cell,
                 made whole across the cell boundaries and assembled as in the crystal (anion at x, y, z - 1)
    colours      MODEL_COLOURS (figkit house colours, periodic in Okabe-Ito bluish green; all pairs checked for
                 colour-vision deficiency), SIDE_COLOURS for "proton on N1" / "proton on O5"
    output       save(fig, stem, width_mm): figkit.qa strict save (600 dpi PNG + vector PDF, text-overlap, font-size,
                 glyph and border checks) plus a 200 dpi print preview in figures/build/

Nothing here writes outside paper/main/figures/.  Inputs are read only.
"""
from __future__ import annotations

import json
import math
import os
import sys
from decimal import Decimal, ROUND_HALF_EVEN
from pathlib import Path

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '2')

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PAPER = REPO / 'paper'
sys.path.insert(0, str(PAPER / '_shared'))

from figkit import style, qa  # noqa: E402

FIG = PAPER / 'main' / 'figures'
BUILD = FIG / 'build'
RES = REPO / 'molecules' / 'sfa' / 'results'
CIF = REPO / 'e-70-0o997-sup1.cif'
FACTS_PATH = PAPER / 'main' / 'data' / 'facts.json'

MM = 1 / 25.4
MINUS = '−'

CATION = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S1',
          'H2', 'H3', 'H5', 'H6', 'H1A', 'H1B', 'H1C', 'H2A', 'H2B']
ANION = ['P1', 'O3', 'O4', 'O5', 'O6', 'H3A', 'H4A']
PAIR_LABELS = CATION + ANION
IDX = {l: i for i, l in enumerate(PAIR_LABELS)}
SYMBOLS = [('H' if l.startswith('H') else l[0]) for l in PAIR_LABELS]

# ---------------------------------------------------------------------------------------------------------------
# colours (validated with the dataviz validate_palette.js: every pair of pcm / gas / periodic passes the simulated
# protan / deutan / tritan separation check, worst OKLab dE x100 = 8.9)
# ---------------------------------------------------------------------------------------------------------------
MODEL_COLOURS = dict(style.MODEL_COLOURS)
assert MODEL_COLOURS['periodic'].upper() == '#009E73', 'figkit.style periodic colour changed; re-validate'
MODEL_MARKERS = {'xray': 'o', 'pcm': 'o', 'gas': 's', 'periodic': '^'}
MODEL_LABELS = {'xray': 'X-ray', 'pcm': 'IEF-PCM (water)', 'gas': 'Gas phase', 'periodic': 'Periodic model'}
SIDE_COLOURS = {'N': '#1F5AA6', 'O': '#E07B00'}      # proton on N1 (salt) / on O5 (neutral pair)
HB_PAIR = '#0E7C7B'          # hydrogen bonds kept in the 27-atom ion pair (figkit draw2d teal)
HB_CRYSTAL = '#B4502D'       # hydrogen bonds of the crystal only (brown-vermillion)


# ---------------------------------------------------------------------------------------------------------------
# facts
# ---------------------------------------------------------------------------------------------------------------
_FACTS = None
USED: dict[str, dict] = {}


def facts():
    global _FACTS
    if _FACTS is None:
        _FACTS = json.loads(FACTS_PATH.read_text())
    return _FACTS


def F(key):
    """Raw value of a fact (recorded as used)."""
    f = facts()[key]
    USED.setdefault(key, {'value': f['value'], 'unit': f.get('unit'), 'status': f.get('status'), 'printed': []})
    return f['value']


def rnd(x, nd):
    """Round half to even at nd decimals on the full-precision value (the rounding rule of facts.json
    _meta.rounding_rule, as build_tables.py), returned as a float."""
    q = Decimal(1).scaleb(-nd)
    return float(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_EVEN))


def num(x, nd, sign=False):
    """Format a number with nd decimals, half to even, U+2212 minus, optional explicit plus."""
    v = rnd(x, nd)
    if v == 0:
        v = 0.0
    s = f'{abs(v):.{nd}f}'
    if v < 0:
        return MINUS + s
    return ('+' + s) if sign else s


def txt(key, nd=None, sign=False):
    """Printed form of a fact at its stored precision (or nd); recorded for the audit."""
    v = F(key)
    p = facts()[key].get('precision') if nd is None else nd
    s = num(v, p, sign=sign) if isinstance(v, (int, float)) and p is not None else str(v)
    USED[key]['printed'].append(s)
    return s


def status(key):
    return facts()[key].get('status') or ''


def write_used(stem, extra=None):
    """Two records of the facts a figure used: figures/<stem>.facts.json (key -> printed strings, the same shape as
    the group-B figures) and figures/build/<stem>_facts_used.json (with values, units and status)."""
    BUILD.mkdir(parents=True, exist_ok=True)
    out = {'figure': stem, 'facts_file': str(FACTS_PATH.relative_to(REPO)), 'facts': USED}
    if extra:
        out.update(extra)
    p = BUILD / f'{stem}_facts_used.json'
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    short = {'script': f'paper/main/figures/scripts/{stem}.py', 'facts_file': out['facts_file'],
             'used': {k: (v['printed'] if v['printed'] else 'checked') for k, v in sorted(USED.items())}}
    if extra:
        short.update(extra)
    (FIG / f'{stem}.facts.json').write_text(json.dumps(short, indent=1, ensure_ascii=False, default=float) + '\n')
    return p


# ---------------------------------------------------------------------------------------------------------------
# geometries
# ---------------------------------------------------------------------------------------------------------------
def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    sym, X = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        sym.append(p[0])
        X.append([float(p[1]), float(p[2]), float(p[3])])
    return sym, np.array(X), lines[1]


def _check_symbols(sym, what):
    assert [s.capitalize() for s in sym] == SYMBOLS, f'{what}: atom order differs from the CIF pair order'


def xray_pair():
    """The deposited ion pair (CIF orthogonal frame, anion at x, y, z - 1)."""
    sym, X, _ = read_xyz(REPO / 'molecules' / 'sfa' / 'input' / 'sfa_pair.xyz')
    _check_symbols(sym, 'sfa_pair.xyz')
    case = json.loads((REPO / 'molecules' / 'sfa' / 'cases' / 'sfa_pair_pcm_d3.json').read_text())
    assert [l.split('#')[0] for l in case['cif_labels']] == PAIR_LABELS
    return X


def case_pair(case):
    sym, X, _ = read_xyz(RES / case / 'opt.xyz')
    _check_symbols(sym, case)
    return X


def scan_point(name):
    sym, X, _ = read_xyz(RES / 'sfa_pair_pcm_d3' / 'scan' / f'{name}.xyz')
    _check_symbols(sym, name)
    return X


def eps_point(eps):
    sym, X, _ = read_xyz(RES / 'sfa_pair_pcm_d3' / 'eps_ladder' / f'eps_{eps}_N' / 'opt.xyz')
    _check_symbols(sym, f'eps {eps}')
    return X


def read_extxyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    head = lines[1]
    lat = np.array([float(x) for x in head.split('Lattice="')[1].split('"')[0].split()]).reshape(3, 3)
    sym, X, lab = [], [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        sym.append(p[0])
        X.append([float(p[1]), float(p[2]), float(p[3])])
        lab.append(p[4])
    return sym, np.array(X), lab, lat


def _mic(v, lat):
    """Minimum-image version of one vector v for lattice rows lat (search of the 27 neighbouring images)."""
    f = np.asarray(v, float) @ np.linalg.inv(lat)
    f -= np.round(f)
    cands = [(f + np.array(s) - 1) @ lat for s in np.ndindex(3, 3, 3)]
    return min(cands, key=np.linalg.norm)


def _whole(X, lat, root, bonded):
    """Make a fragment whole: breadth-first along bonds, each atom at the minimum image of its parent."""
    pos = {root: X[root].copy()}
    todo = [root]
    while todo:
        i = todo.pop(0)
        for j in bonded[i]:
            if j not in pos:
                pos[j] = pos[i] + _mic(X[j] - X[i], lat)
                todo.append(j)
    return pos


def periodic_pair(name='relax_i', fu=1):
    """Ion pair of formula unit `fu` from periodic_qe/<name>.extxyz, as a 27 x 3 array in PAIR_LABELS order:
    cation whole, anion whole and placed at the minimum image of O5 from N1 (the N1-H1A...O5 partner).
    Returns (X, lattice)."""
    sym, X, lab, lat = read_extxyz(RES / 'periodic_qe' / f'{name}.extxyz')
    suf = '' if fu == 1 else f'_{fu}'
    L = {l: i for i, l in enumerate(lab)}
    idx = [L[l + suf] for l in PAIR_LABELS]
    Y = X[idx]
    rad = style.COVALENT_RADII
    bonded = {i: [] for i in range(len(Y))}
    for i in range(len(Y)):
        for j in range(i + 1, len(Y)):
            d = np.linalg.norm(_mic(Y[j] - Y[i], lat))
            if d < 1.2 * (rad[SYMBOLS[i]] + rad[SYMBOLS[j]]):
                bonded[i].append(j)
                bonded[j].append(i)
    cat = [IDX[l] for l in CATION]
    an = [IDX[l] for l in ANION]
    pc = _whole(Y, lat, IDX['C1'], {i: [j for j in bonded[i] if j in cat] for i in cat})
    pa = _whole(Y, lat, IDX['P1'], {i: [j for j in bonded[i] if j in an] for i in an})
    assert len(pc) == len(cat) and len(pa) == len(an), 'fragment not connected'
    out = np.zeros_like(Y)
    for i, p in pc.items():
        out[i] = p
    shift = (pc[IDX['N1']] + _mic(pa[IDX['O5']] - pc[IDX['N1']], lat)) - pa[IDX['O5']]
    for i, p in pa.items():
        out[i] = p + shift
    return out, lat


def cell_shift_removed(name='relax_i'):
    """Periodic pair (formula unit 1) translated by minus the mean heavy-atom displacement of the whole cell from
    the X-ray start (start_i_xray), i.e. the frame of the fact periodic.rmsd_heavy (32 atoms, mean shift removed).
    Returns (pair, X-ray pair from the same start file)."""
    sym, X, lab, lat = read_extxyz(RES / 'periodic_qe' / f'{name}.extxyz')
    sym0, X0, lab0, lat0 = read_extxyz(RES / 'periodic_qe' / 'start_i_xray.extxyz')
    assert lab == lab0 and np.allclose(lat, lat0)
    heavy = [k for k, s in enumerate(sym) if s != 'H']
    d = np.array([_mic(X[k] - X0[k], lat) for k in heavy])
    t = d.mean(0)
    rmsd = float(np.sqrt(np.mean(np.sum((d - t) ** 2, axis=1))))
    P, _ = periodic_pair(name)
    R, _ = periodic_pair('start_i_xray')
    # bring the relaxed pair next to the start pair (same lattice image of C1)
    P = P - P[IDX['C1']] + R[IDX['C1']] + _mic(P[IDX['C1']] - R[IDX['C1']], lat)
    return P - t, R, rmsd


def dist(X, a, b):
    return float(np.linalg.norm(X[IDX[a]] - X[IDX[b]]))


def delta(X):
    return dist(X, 'N1', 'H1A') - dist(X, 'H1A', 'O5')


# ---------------------------------------------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------------------------------------------
def save(fig, stem, width_mm, extra_used=None, strict=True):
    """figkit.qa strict save to paper/main/figures/<stem>.png/.pdf + print preview + facts-used record."""
    if os.environ.get('SAPH_QA_LAX'):           # drafting only: report QA issues instead of stopping
        strict = False
    out = qa.save(fig, FIG / stem, strict=strict, max_width_mm=width_mm, glyphs=True)
    BUILD.mkdir(parents=True, exist_ok=True)
    from PIL import Image
    png = FIG / f'{stem}.png'
    with Image.open(png) as im:
        w, h = im.size
        dpi = im.info.get('dpi', (600, 600))[0]
        s = 200.0 / float(dpi)
        im.convert('RGB').resize((max(1, round(w * s)), max(1, round(h * s))), Image.LANCZOS).save(
            BUILD / f'{stem}_print200.png', dpi=(200, 200))
        size_mm = (w / dpi * 25.4, h / dpi * 25.4)
    used = write_used(stem, extra_used)
    print(f'{stem}: {size_mm[0]:.1f} x {size_mm[1]:.1f} mm; wrote', ', '.join(str(p.relative_to(REPO)) for p in out),
          '+', used.relative_to(REPO))
    return out


def setup():
    style.apply()
    import matplotlib.pyplot as plt
    plt.rcParams['savefig.bbox'] = 'tight'
    return plt


def figure_mm(w_mm, h_mm):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM))
    fig._saph_mm = (w_mm, h_mm)
    return fig


def axes_mm(fig, x, y, w, h, **kw):
    """Axes at (x, y) mm from the lower-left corner of the figure, size w x h mm."""
    W, H = fig._saph_mm
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def panel_letter(fig, x_mm, y_mm, letter, title=None, fontsize=None):
    """Bold (a) at (x, y) mm (baseline), optional title after it in regular weight."""
    W, H = fig._saph_mm
    t = fig.text(x_mm / W, y_mm / H, f'({letter})', fontsize=style.PANEL_PT, fontweight='bold', ha='left',
                 va='baseline')
    out = [t]
    if title:
        fig.canvas.draw()
        bb = t.get_window_extent(fig.canvas.get_renderer())
        x1 = bb.x1 / fig.dpi / MM + 1.6
        out.append(fig.text(x1 / W, y_mm / H, title, fontsize=fontsize or style.BASE_PT, ha='left', va='baseline'))
    return out
