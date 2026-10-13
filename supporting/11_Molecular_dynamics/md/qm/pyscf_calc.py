"""Adapter (a), second half: a PySCF calculation that stands for the Gaussian job of an MCPB.py input.

Runs with /root/cqc-sfa-tools/env/bin/python (PySCF 2.14.0).

What "B3LYP/6-31G* Integral=(Grid=UltraFine)" becomes here, and where it differs from Gaussian:

  functional   xc = 'b3lyp'.  In PySCF 2.14 this is libxc functional 402 (HYB_GGA_XC_B3LYP), the VWN-RPA
               ("VWN3") form, i.e. the B3LYP of Gaussian; PySCF before 2.3 used VWN5 under this name.
               `assert_b3lyp_is_gaussian()` checks it at run time.
  basis        '6-31g*' of PySCF's library: H 2s; C, N, O 3s2p1d; S 4s3p1d; Zn 5s4p2d1f (Rassolov et al.).
               PURE d and f functions (5D, 7F).  Gaussian takes Cartesian d functions (6D) for 6-31G*.
               Cartesian functions were tried first and dropped: PySCF then takes the fitting basis as
               Cartesian too, and its Coulomb metric is singular to working precision (Cholesky fails,
               5 eigenvalues below 1e-12 for the small model; md/qm/metric_conditioning.json).  PySCF
               falls back to an eigen-decomposition with a cut-off, which works for one geometry but is
               a threshold decision that can change between displaced geometries.  With pure functions
               the metric is positive definite (smallest eigenvalue 8e-7).  The price is a basis that
               differs from Gaussian's by one s-type combination per d shell.
  grid         (99, 590) points per atom, the size of Gaussian's UltraFine; pruning and atomic partition are
               PySCF's (NWChem pruning, Becke weights with Treutler-Ahlrichs radii), not Gaussian's.
  integrals    density fitting of Coulomb and exchange with def2-universal-JKFIT (Gaussian: exact
               four-centre integrals).  The error is measured in md/qm/test_small_chain.py.
  SCF          closed shell (RKS), conv_tol 1e-10 Hartree and conv_tol_grad 1e-6.

Everything a job needs is in `Settings`; the drivers store it in their JSON records.
"""
import dataclasses
import os
import time

import numpy as np

import gaucom

ANG = gaucom.BOHR


@dataclasses.dataclass
class Settings:
    xc: str = 'b3lyp'
    basis: str = '6-31g*'
    cart: bool = False
    grid: tuple = (99, 590)
    density_fit: bool = True
    auxbasis: str = 'def2-universal-jkfit'
    conv_tol: float = 1e-10
    conv_tol_grad: float = 1e-6
    max_cycle: int = 150
    max_memory_mb: int = 9000
    verbose: int = 4


def assert_b3lyp_is_gaussian():
    from pyscf.dft import libxc
    hyb, fn = libxc.parse_xc('b3lyp')
    if tuple(fn) != ((402, 1),):
        raise RuntimeError(f"xc 'b3lyp' is not libxc functional 402 in this PySCF: {fn}")
    return 'libxc 402 HYB_GGA_XC_B3LYP (VWN-RPA), the B3LYP of Gaussian'


def method_of_com(com, settings=None):
    """Settings for the route of a .com file; stops if the route asks for something else."""
    s = settings or Settings()
    if com.method.upper() != 'B3LYP':
        raise ValueError(f'method {com.method} in {com.path}: only B3LYP is set up')
    if com.basis.upper() not in ('6-31G*', '6-31G(D)'):
        raise ValueError(f'basis {com.basis} in {com.path}: only 6-31G* is set up')
    if com.keywords.get('grid') != 'ultrafine':
        raise ValueError(f'{com.path}: Integral=(Grid=UltraFine) expected')
    return s


def build_mol(symbols, coords_ang, charge, multiplicity, settings, output=None):
    from pyscf import gto
    mol = gto.Mole()
    mol.atom = [(s, tuple(float(x) for x in xyz)) for s, xyz in zip(symbols, coords_ang)]
    mol.unit = 'Angstrom'
    mol.charge = charge
    mol.spin = multiplicity - 1
    mol.basis = settings.basis
    mol.cart = settings.cart
    mol.symmetry = False
    mol.max_memory = settings.max_memory_mb
    mol.verbose = settings.verbose
    if output:
        mol.output = output
    mol.build()
    return mol


def build_mf(mol, settings, chkfile=None):
    from pyscf import dft
    if mol.spin != 0:
        raise ValueError('closed-shell set-up only')
    mf = dft.RKS(mol)
    mf.xc = settings.xc
    mf.grids.atom_grid = tuple(settings.grid)
    mf.conv_tol = settings.conv_tol
    mf.conv_tol_grad = settings.conv_tol_grad
    mf.max_cycle = settings.max_cycle
    if settings.density_fit:
        mf = mf.density_fit(auxbasis=settings.auxbasis)
    if chkfile:
        mf.chkfile = chkfile
    return mf


def run_scf(mf, dm0=None):
    t = time.time()
    e = mf.kernel(dm0=dm0)
    if not mf.converged:
        # second-order SCF from the last density; reported by the caller
        mf2 = mf.newton()
        e = mf2.kernel(dm0=mf.make_rdm1())
        if not mf2.converged:
            raise RuntimeError('SCF not converged, also not with the second-order solver')
        mf.mo_coeff, mf.mo_energy, mf.mo_occ, mf.e_tot, mf.converged = (
            mf2.mo_coeff, mf2.mo_energy, mf2.mo_occ, mf2.e_tot, True)
    return e, time.time() - t


def esp_nuclear(mol, points_bohr):
    z = mol.atom_charges()
    r = mol.atom_coords()
    d = np.linalg.norm(points_bohr[:, None, :] - r[None, :, :], axis=2)
    return (z[None, :] / d).sum(axis=1)


def esp_electronic(mol, dm, points_bohr, batch_mb=1500, method='int3c2e'):
    """- integral rho(r') / |r - r'| dr' at the points (atomic units; negative for electrons).

    method 'int3c2e': three-centre integrals with a unit point charge at every point (packed in mu <= nu);
    method 'rinv': one-electron integrals with the origin of 1/r moved to the point (slow; for checks).
    """
    from pyscf import df, gto
    npts = len(points_bohr)
    out = np.empty(npts)
    nao = mol.nao
    if method == 'rinv':
        for k, p in enumerate(points_bohr):
            with mol.with_rinv_origin(p):
                out[k] = -np.einsum('ij,ij->', mol.intor('int1e_rinv'), dm)
        return out
    npair = nao * (nao + 1) // 2
    # packed density: off-diagonal elements count twice
    idx = np.tril_indices(nao)
    dtri = dm[idx] * 2.0
    dtri[idx[0] == idx[1]] *= 0.5
    step = max(8, int(batch_mb * 1e6 / (8 * npair)))
    for k0 in range(0, npts, step):
        fake = gto.fakemol_for_charges(points_bohr[k0:k0 + step])
        j3c = df.incore.aux_e2(mol, fake, intor='int3c2e', aosym='s2ij')
        out[k0:k0 + step] = -dtri @ j3c
    return out


def esp_at_points(mol, dm, points_bohr, **kw):
    return esp_nuclear(mol, points_bohr) + esp_electronic(mol, dm, points_bohr, **kw)


def esp_at_nuclei(mol, dm):
    """Potential at each nucleus without the nucleus itself (what Gaussian prints on the 'Atom' lines)."""
    z = mol.atom_charges()
    r = mol.atom_coords()
    out = np.empty(mol.natm)
    for a in range(mol.natm):
        with mol.with_rinv_origin(r[a]):
            ve = -np.einsum('ij,ij->', mol.intor('int1e_rinv'), dm)
        d = np.linalg.norm(r - r[a], axis=1)
        d[a] = np.inf
        out[a] = ve + (z / d).sum()
    return out


def versions():
    import pyscf
    import scipy
    v = dict(pyscf=pyscf.__version__, numpy=np.__version__, scipy=scipy.__version__)
    try:
        import geometric
        v['geometric'] = geometric.__version__
    except Exception:
        pass
    from pyscf import lib
    v['threads'] = lib.num_threads()
    v['omp_num_threads_env'] = os.environ.get('OMP_NUM_THREADS')
    return v
