"""
Step 3 -- final single point at the optimised geometry: FMO / global reactivity descriptors,
Molden file for Multiwfn, and cube files (HOMO, LUMO, density, ESP) for 3-D rendering.
Usage: python dft/fmo.py <case>
"""
import sys, os, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from common import *
from pyscf.tools import molden, cubegen

case = load_case(sys.argv[1]); d = case['dir']
mol = mol_for_case(case, d / 'opt.xyz')
mf = make_mf(mol, conv_tol=1e-10, chkfile=d / 'final.chk', solvent=case.get('solvent'))
run_scf(mf, init_chk=d / 'final.chk' if (d / 'final.chk').exists() else d / 'opt.chk')

n = mol.nelectron // 2; e = mf.mo_energy * HARTREE_EV
H, L = e[n-1], e[n]
IP, EA = -H, -L                       # Koopmans-type approximations, as used in these papers
chi = (IP + EA) / 2; mu = -chi; eta = (IP - EA) / 2; S = 1 / (2 * eta); omega = mu**2 / (2 * eta)
dip = mf.dip_moment(unit='Debye', verbose=0)
desc = {'E_total_Ha': float(mf.e_tot), 'HOMO_eV': float(H), 'LUMO_eV': float(L), 'gap_eV': float(L - H),
        'HOMO-1_eV': float(e[n-2]), 'LUMO+1_eV': float(e[n+1]),
        'IP_eV': float(IP), 'EA_eV': float(EA), 'electronegativity_chi_eV': float(chi),
        'chemical_potential_mu_eV': float(mu), 'hardness_eta_eV': float(eta), 'softness_S_eV-1': float(S),
        'electrophilicity_omega_eV': float(omega), 'dipole_Debye': float(np.linalg.norm(dip)),
        'dipole_vector_Debye': [float(v) for v in dip], 'n_basis': int(mol.nao), 'n_electrons': int(mol.nelectron),
        'mo_energies_eV': [float(v) for v in e], 'n_occ': int(n)}
# where do HOMO/LUMO live? Mulliken-type fragment composition of each frontier orbital
S_ao = mol.intor('int1e_ovlp'); slices = mol.aoslice_by_atom()
def frag_share(mo, atoms):
    idx = np.concatenate([np.arange(slices[a][2], slices[a][3]) for a in atoms])
    return float((mo[idx] * (S_ao @ mo)[idx]).sum())
if case.get('fragments'):
    desc['fragment_composition'] = {name: {frag: round(100 * frag_share(mf.mo_coeff[:, k], atoms), 1) for frag, atoms in case['fragments'].items()}
                                    for name, k in [('HOMO-1', n-2), ('HOMO', n-1), ('LUMO', n), ('LUMO+1', n+1)]}
save_json(d / 'fmo.json', desc)
stage_note(case, f"FMO: HOMO {H:.3f} eV, LUMO {L:.3f} eV, gap {L-H:.3f} eV; fragment composition {desc.get('fragment_composition', {}).get('HOMO')} / {desc.get('fragment_composition', {}).get('LUMO')}")
with open(d / 'fmo_table.md', 'w') as f:
    f.write('| Descriptor | Formula | Value |\n|---|---|---|\n')
    for k, form, v, u in [('E(HOMO)', '', H, 'eV'), ('E(LUMO)', '', L, 'eV'), ('Energy gap', 'E(LUMO) - E(HOMO)', L-H, 'eV'),
                          ('Ionization potential (IP)', '-E(HOMO)', IP, 'eV'), ('Electron affinity (EA)', '-E(LUMO)', EA, 'eV'),
                          ('Electronegativity (chi)', '(IP+EA)/2', chi, 'eV'), ('Chemical potential (mu)', '-(IP+EA)/2', mu, 'eV'),
                          ('Chemical hardness (eta)', '(IP-EA)/2', eta, 'eV'), ('Chemical softness (S)', '1/(2 eta)', S, 'eV^-1'),
                          ('Electrophilicity index (omega)', 'mu^2/(2 eta)', omega, 'eV'), ('Dipole moment', '', np.linalg.norm(dip), 'Debye')]:
        f.write(f'| {k} | {form} | {v:.4f} {u} |\n')
print(open(d / 'fmo_table.md').read())

# Multiwfn segfaults while building the density matrix if the molden holds fewer MOs than basis
# functions.  PySCF drops near-linearly-dependent combinations (the diffuse ++ set on a 49-atom cluster
# loses 6 of 853), so mo_coeff is 853x847 and Multiwfn -- which assumes a square matrix -- walks off the
# end.  Pad the discarded vectors back with zeros: they carry no occupancy, so no physics changes.
nao, nmo = mf.mo_coeff.shape
if nmo < nao:
    pad = nao - nmo
    print(f'[fmo] molden: padding {pad} linearly-dependent vectors with zeros '
          f'({nmo} MOs -> {nao}) so Multiwfn sees a square coefficient matrix', flush=True)
    C = np.zeros((nao, nao)); C[:, :nmo] = mf.mo_coeff
    e = np.zeros(nao); e[:nmo] = mf.mo_energy; e[nmo:] = mf.mo_energy.max() + 1e3
    occ = np.zeros(nao); occ[:nmo] = mf.mo_occ
    molden.from_mo(mol, str(d / 'final.molden'), C, ene=e, occ=occ)
    stage_note(case, f'molden padded from {nmo} to {nao} MOs for Multiwfn compatibility')
else:
    molden.from_scf(mf, str(d / 'final.molden'))
print('[fmo] wrote final.molden (input for Multiwfn)')
# cube files: 80^3 grid with 4 A margin is enough for publication renders
kw = dict(nx=80, ny=80, nz=80, margin=4.5)
import time as _t
_t0 = _t.time()
cubegen.orbital(mol, str(d / 'homo.cube'), mf.mo_coeff[:, n-1], **kw)
cubegen.orbital(mol, str(d / 'lumo.cube'), mf.mo_coeff[:, n], **kw)
dm = mf.make_rdm1()
cubegen.density(mol, str(d / 'density.cube'), dm, **kw)
# ESP is ~100x costlier per point than the density: coarser grid, interpolated onto the density isosurface when rendering
print(f'[fmo] orbital/density cubes done ({_t.time()-_t0:.0f} s); ESP cube (the slow one, ~100x a density cube) ...', flush=True)
_t1 = _t.time()
# ESP is used only to colour the 3-D MEP render; Multiwfn computes its own surface ESP from the molden,
# so the reported MEP numbers do not depend on this grid.  ESP_CUBE_N trades render sharpness for time.
_n = int(os.environ.get('ESP_CUBE_N', '64'))
cubegen.mep(mol, str(d / 'esp.cube'), dm, nx=_n, ny=_n, nz=_n, margin=4.5)
print(f'[fmo] wrote homo/lumo/density cubes and the {_n}^3 ESP cube ({_t.time()-_t1:.0f} s)', flush=True)
stage_note(case, f'cubes written (ESP {_n}^3 in {(_t.time()-_t1)/60:.1f} min)')
