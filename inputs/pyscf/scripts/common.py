"""
Shared setup for the DFT project on 4-sulfamoylanilinium dihydrogen phosphate (adapted from the sister
project cqc-dft, which studied 3ATB and three further proton-transfer salts).

Level of theory requested by the collaborator: B3LYP/6-311++G(d,p) (Gaussian 09 style).
Engine used here: PySCF (open source).  Notes on equivalence:
  * PySCF's 'b3lyp' keyword maps to libxc HYB_GGA_XC_B3LYP, which uses the VWN-RPA (VWN "functional III")
    local correlation -- the same definition Gaussian uses.  ('b3lyp5' would be the VWN5 variant; not used.)
  * 6-311++G(d,p) with spherical (5-component) d functions, as Gaussian uses for 6-311G-type basis sets.
  * Density fitting (RI-JK, auxiliary basis chosen automatically by PySCF) speeds up the calculation ~5x
    with energy errors ~1e-5 Hartree, far below anything reported here.
  * Integration grid: PySCF level 3 (~75 radial x 302 angular points per atom, pruned), comparable to
    Gaussian 09's default "FineGrid".
"""
import os, json, time, pathlib
import numpy as np
from pyscf import gto, dft, lib

ROOT = pathlib.Path(__file__).resolve().parent.parent
HARTREE_EV = 27.211386245988
BOHR_A = 0.529177210903
XC = 'b3lyp'
BASIS = '6-311++g(d,p)'
NTHREADS = int(os.environ.get('DFT_THREADS', '8'))
MAX_MEMORY = int(os.environ.get('DFT_MAX_MEMORY_MB', '4500'))   # MB; set per machine via the environment

# Some elements have no 6-311G-family basis at all (zinc is one: it exists neither in PySCF nor in the
# Basis Set Exchange).  A case may therefore name a different basis for individual elements; everything
# else stays at 6-311++G(d,p).
DEFAULT_BASIS_OVERRIDES = {}


def find_case_file(name):
    """Case JSONs live in molecules/<id>/cases/ (current layout) or cases/ (the original 3ATB layout)."""
    for c in sorted((ROOT / 'molecules').glob(f'*/cases/{name}.json')) + [ROOT / 'cases' / f'{name}.json']:
        if c.exists():
            return c
    raise FileNotFoundError(f'no case file for {name!r} under molecules/*/cases/ or cases/')


def load_case(name):
    """A case is a JSON file describing one molecule / ion pair / hydrogen-bonded cluster to study."""
    path = find_case_file(name)
    with open(path) as f:
        case = json.load(f)
    case['name'] = name
    case['case_file'] = path
    # results sit next to the case definition: molecules/<id>/results/<case>/, else results/<case>/
    if path.parent.name == 'cases' and path.parent.parent.parent.name == 'molecules':
        case['dir'] = path.parent.parent / 'results' / name
        case['molecule'] = path.parent.parent.name
    else:
        case['dir'] = ROOT / 'results' / name
        case['molecule'] = name
    case['dir'].mkdir(parents=True, exist_ok=True)
    # Level-of-theory knobs travel with the case file, not with the shell that launched the run, so every
    # process of a case (optimiser, finite-difference workers, analyses) uses the same definition.  Exported
    # to the environment because make_mf() does not receive the case.
    os.environ['DFT_DISP'] = str(case.get('disp') or 'none')
    os.environ['DFT_XC'] = str(case.get('xc') or XC)
    return case


def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0])
    syms, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        syms.append(p[0]); xyz.append([float(v) for v in p[1:4]])
    return syms, np.array(xyz)


def write_xyz(path, syms, xyz, comment=''):
    with open(path, 'w') as f:
        f.write(f'{len(syms)}\n{comment}\n')
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f'{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n')


def build_mol(xyz_path, charge=0, spin=0, basis=BASIS, verbose=3, basis_overrides=None):
    """basis_overrides maps an element symbol to its own basis, e.g. {'Zn': 'def2-tzvp'}."""
    lib.num_threads(NTHREADS)
    ov = dict(DEFAULT_BASIS_OVERRIDES)
    ov.update(basis_overrides or {})
    b = basis if not ov else dict({'default': basis}, **ov)
    return gto.M(atom=str(xyz_path), basis=b, charge=charge, spin=spin,
                 verbose=verbose, max_memory=MAX_MEMORY, symmetry=False)


def mol_for_case(case, xyz_path, verbose=3):
    """Build the PySCF Mole for a case, honouring a case-level 'basis' (e.g. "def2-tzvp" for the basis-set
    sensitivity twins *_tz; absent = the project basis 6-311++G(d,p)) and any per-element override."""
    return build_mol(xyz_path, case['charge'], case['spin'], basis=case.get('basis', BASIS), verbose=verbose,
                     basis_overrides=case.get('basis_overrides'))


# Printable names for the basis keywords PySCF accepts in lower case.
BASIS_LABELS = {'6-311++g(d,p)': '6-311++G(d,p)', '6-311++g**': '6-311++G(d,p)', 'def2-tzvp': 'def2-TZVP',
                'def2-tzvpd': 'def2-TZVPD', 'def2-svp': 'def2-SVP', '6-31g(d)': '6-31G(d)', 'sto-3g': 'STO-3G'}


def basis_label(basis):
    return BASIS_LABELS.get(str(basis).lower(), str(basis))


def basis_description(case):
    """Human-readable level of theory for the methods paragraph."""
    base = basis_label(case.get('basis', BASIS))
    ov = case.get('basis_overrides') or {}
    if not ov:
        return base
    rest = ', '.join(sorted(ov))
    return (base + ' for all atoms except ' + rest + ' ('
            + '; '.join(f'{k}: {v}' for k, v in sorted(ov.items())) + ')')


SOLVENT_EPS = {'water': 78.3553, 'dmso': 46.826, 'methanol': 32.613, 'ethanol': 24.852}

# Quality knobs.  The 3ATB run used grid level 3 (~Gaussian FineGrid) and PCM Lebedev order 17, both
# chosen to keep an 8 GB laptop usable.  On a machine with real memory, level 4 (~UltraFine) and PySCF's
# own default Lebedev order 29 are the better-converged settings.  Set GRID_LEVEL / PCM_LEBEDEV_ORDER.
GRID_LEVEL = int(os.environ.get('GRID_LEVEL', '4'))


def make_mf(mol, conv_tol=1e-9, grid_level=None, chkfile=None, solvent=None):
    """B3LYP/6-311++G(d,p) with density fitting; optional IEF-PCM continuum solvent (as Gaussian's SCRF=PCM)."""
    mf = (dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)).density_fit()
    mf.xc = os.environ.get('DFT_XC', XC)
    # Empirical dispersion (Grimme D3 with Becke-Johnson damping for B3LYP-D3(BJ)), from the case's 'disp'.
    # Set on the bare KS object, before any solvent wrapper, so energy and gradients both include it.
    disp = os.environ.get('DFT_DISP', 'none')
    if disp and disp.lower() != 'none':
        mf.disp = disp
    mf.grids.level = GRID_LEVEL if grid_level is None else grid_level
    mf.conv_tol = conv_tol
    mf.max_cycle = 150
    if chkfile:
        mf.chkfile = str(chkfile)
    if solvent:
        mf = mf.PCM()
        mf.with_solvent.method = 'IEF-PCM'
        mf.with_solvent.eps = SOLVENT_EPS[solvent.lower()]
        # cavity surface density: PySCF's default Lebedev order 29 (302 points/atom sphere) makes every SCF cycle ~10x
        # slower for 400 basis functions; order 17 (110 points/atom, ~2700 surface points for this pair) is comparable
        # to Gaussian's default tessellation density and is what we use.
        mf.with_solvent.lebedev_order = int(os.environ.get('PCM_LEBEDEV_ORDER', '29'))
    return mf


def stage_note(case, text):
    """Append a timestamped line to results/<case>/STAGES.md (the human-readable stage log)."""
    import datetime
    with open(case['dir'] / 'STAGES.md', 'a') as f:
        f.write(f"- {datetime.datetime.now():%Y-%m-%d %H:%M}  {text}\n")


def run_scf(mf, init_chk=None):
    """Run SCF, optionally starting from the density in a checkpoint file."""
    dm0 = None
    if init_chk and os.path.exists(init_chk):
        try:
            from pyscf import scf
            dm0 = scf.hf.init_guess_by_chkfile(mf.mol, str(init_chk))
        except Exception:
            dm0 = None
    t = time.time()
    mf.kernel(dm0=dm0)
    mf._wall = time.time() - t
    if not mf.converged:
        raise RuntimeError('SCF did not converge')
    return mf


def atom_labels(syms):
    """Element+index labels (C1, C2, ...) so tables are readable: numbering is per element, 1-based."""
    cnt, out = {}, []
    for s in syms:
        cnt[s] = cnt.get(s, 0) + 1
        out.append(f'{s}{cnt[s]}')
    return out


# Cordero et al. (2008) single-bond covalent radii, in Angstrom.  Keep every element that appears in any
# case here: a missing entry silently falls back to 0.8 A, which for zinc put the Zn-Cl cutoff at 2.28 A and
# made the connectivity check report the intact ZnCl4 unit as a broken bond.
COV_R = {'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66, 'Mg': 1.41,
         'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'K': 2.03, 'Ca': 1.76, 'Mn': 1.39,
         'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Br': 1.20, 'Cd': 1.44, 'I': 1.39}

def bonds_from_geometry(syms, xyz, scale=1.25):
    """Covalent bond list (i<j) from distance criterion; used for bond-length tables and rendering."""
    b = []
    for i in range(len(syms)):
        for j in range(i + 1, len(syms)):
            d = np.linalg.norm(xyz[i] - xyz[j])
            if d < scale * (COV_R.get(syms[i], 0.8) + COV_R.get(syms[j], 0.8)):
                b.append((i, j, d))
    return b


def save_json(path, obj):
    with open(path, 'w') as f:
        json.dump(obj, f, indent=2, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))


def hbond_geometry(xyz, D, H, A):
    """D-H, H...A, D...A distances (A) and D-H...A angle (deg) for atom indices D, H, A."""
    dh = xyz[H] - xyz[D]; ha = xyz[A] - xyz[H]; da = xyz[A] - xyz[D]
    v1 = xyz[D] - xyz[H]; v2 = xyz[A] - xyz[H]
    ang = np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
    return {'D-H_A': float(np.linalg.norm(dh)), 'H..A_A': float(np.linalg.norm(ha)), 'D..A_A': float(np.linalg.norm(da)), 'DHA_deg': float(ang)}
