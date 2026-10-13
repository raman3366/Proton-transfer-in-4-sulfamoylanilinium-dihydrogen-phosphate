"""
Counterpoise-corrected interaction energy of the M1 ion pair, and the deformation energies of its fragments
against the separately relaxed fragments -- gas phase and IEF-PCM water.

Usage: python dft/interaction.py <case> [--partition auto|ionic|neutral] [--frag-cases A,B] [common switches]
Writes results/<case>/interaction[_<tag>].json and .md.

Definitions (Boys-Bernardi counterpoise; A, B = fragments; "cplx" = geometry inside the pair):
    E_int(raw) = E_AB - E_A(cplx; A basis)  - E_B(cplx; B basis)
    E_int(CP)  = E_AB - E_A(cplx; AB basis) - E_B(cplx; AB basis)        (partner atoms as ghosts)
    BSSE       = E_int(CP) - E_int(raw)  (>= 0)
    E_def(X)   = E_X(cplx; X basis) - E_X(relaxed; X basis)
    E_bind     = E_int(CP) + E_def(A) + E_def(B)   (pair relative to the relaxed, separated fragments)
Energies in kcal/mol and kJ/mol; the D3(BJ) part is listed separately.

Choices
  * Fragments: the atom lists of case['fragments'] (M1: 4-sulfamoylanilinium cation + H2PO4-).  If the
    optimisation moved the proton H1A to O5 (gas phase can), --partition auto (default) switches to the neutral
    partition (sulfanilamide + H3PO4: H1A moves to the acceptor's fragment), because a "cation" whose N-H proton
    sits on the phosphate is not a meaningful fragment.  The partition used is recorded.
  * Fragment charges and relaxed geometries come from the fragment case files (molecules/sfa/cases/): the case
    whose input has the same element sequence as the fragment (and, for the relaxed geometry, the same solvent,
    dispersion and basis as the pair) -- cation_<env>, h2po4_<env>, sulfanilamide_<env>, h3po4_<env>.  Override
    with --frag-cases.  Without a relaxed geometry the deformation terms are skipped and say so.
  * D3(BJ) is evaluated on real-atom molecules only (simple-dftd3 cannot take ghost atoms -- it tried to allocate
    25 GB on the Pi) and added to each term; dispersion carries no BSSE.  The SCF parts run without D3.
  * PCM (approximate, stated in the JSON): every species sits in its own cavity -- the pair in the pair cavity,
    each fragment in a cavity built from its real atoms only (ghost atoms get no sphere; PySCF would otherwise
    give a ghost atom the radius of element 0).  So E_int(PCM) is an interaction energy in solution, including
    the partial desolvation of the contact.  Ghost basis functions reach into the region the partner occupied,
    which in the fragment calculation is solvent, so part of the fragment density can leak outside its cavity:
    the counterpoise correction in PCM is therefore only approximate, and the gas-phase BSSE at the same
    geometry is reported next to it as a cross-check.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from props_common import *          # noqa: F401,F403 -- pins threads / scratch, so it must precede numpy and PySCF
from pyscf import gto, dft as _dft, lib
from pyscf.data.elements import charge as _z
import pyscf.solvent.pcm as _pcm

_gen_surface = _pcm.gen_surface


def _gen_surface_real_atoms(mol, *a, **k):
    """PCM cavity from the real atoms only (ghost atoms have zero nuclear charge)."""
    real = [i for i in range(mol.natm) if mol.atom_charge(i) > 0]
    if len(real) == mol.natm:
        return _gen_surface(mol, *a, **k)
    z = sum(int(mol.atom_charge(i)) for i in real)
    sub = gto.M(atom=[(mol.atom_pure_symbol(i), mol.atom_coord(i)) for i in real], unit='Bohr', basis='sto-3g',
                spin=z % 2, verbose=0)
    return _gen_surface(sub, *a, **k)


_pcm.gen_surface = _gen_surface_real_atoms

ap = parser(__doc__)
ap.add_argument('--partition', choices=['auto', 'ionic', 'neutral'], default='auto')
ap.add_argument('--frag-cases', help='comma-separated fragment case names, in case["fragments"] order')
ap.add_argument('--relaxed-xyz', help='comma-separated relaxed fragment geometries (else from the fragment cases)')
ap.add_argument('--no-gas-bsse', action='store_true', help='PCM cases: skip the gas-phase BSSE cross-check')
args = ap.parse_args()
t_start = time.time()
case, xyz, outdir, tag = setup(args)
if not case.get('fragments') or len(case['fragments']) != 2:
    raise SystemExit('[interaction] the case must define exactly two fragments')
syms, X = read_xyz(xyz)
basis = case.get('basis', BASIS)
solv = case.get('solvent')

# ---- partition -----------------------------------------------------------------------------------------------
fr = {k: list(v) for k, v in case['fragments'].items()}
transferred = None
if case.get('proton'):
    H, N, O = case['proton']['H'], case['proton']['donor_N'], case['proton']['acceptor_O']
    transferred = bool(np.linalg.norm(X[H] - X[O]) < np.linalg.norm(X[H] - X[N]))
partition = args.partition
if partition == 'auto':
    partition = 'neutral' if transferred else 'ionic'
if partition == 'neutral':
    if not case.get('proton'):
        raise SystemExit('[interaction] neutral partition needs case["proton"]')
    dn = next(k for k, v in fr.items() if H in v); ac = next(k for k, v in fr.items() if O in v)
    fr[dn] = [i for i in fr[dn] if i != H]; fr[ac] = fr[ac] + [H]
names = list(fr)


def all_cases():
    for p in sorted((ROOT / 'molecules').glob('*/cases/*.json')):
        c = json.load(open(p)); c['name'] = p.stem
        yield c


def frag_case(idx, forced=None):
    """(charge, relaxed-case name or None) for a fragment given by pair atom indices."""
    seq = [syms[i].capitalize() for i in idx]
    cands = []
    for c in all_cases():
        if forced and c['name'] != forced:
            continue
        s, _ = read_xyz(ROOT / c['xyz'])
        if [a.capitalize() for a in s] == seq:
            cands.append(c)
    if not cands:
        raise SystemExit(f'[interaction] no fragment case with element sequence {"".join(seq)}; give --frag-cases')
    same = [c for c in cands if c.get('solvent') == solv and (c.get('disp') or 'none') == (case.get('disp') or 'none')
            and c.get('basis', BASIS) == basis]
    return cands[0]['charge'], (same[0]['name'] if same else None), cands[0].get('chemical_name')


forced = args.frag_cases.split(',') if args.frag_cases else [None, None]
fmeta = {}
for k, f in zip(names, forced):
    q, rc, chem = frag_case(fr[k], f)
    fmeta[k] = {'atoms': fr[k], 'charge': q, 'relaxed_case': rc, 'chemical_name': chem}
if partition == 'neutral':
    # [review 28 Sep 2026] case['fragments'] names the IONS; after H1A has moved to the acceptor the fragments are
    # neutral sulfanilamide and H3PO4, so label them by the matched fragment case (its chemical_name).  Before this,
    # the JSON and .md reported e.g. "E_def(4-sulfamoylanilinium cation)" for neutral sulfanilamide.
    ren = {k: (fmeta[k]['chemical_name'] or k) for k in names}
    fr = {ren[k]: v for k, v in fr.items()}
    fmeta = {ren[k]: v for k, v in fmeta.items()}
    names = list(fr)
if sum(v['charge'] for v in fmeta.values()) != case['charge']:
    raise SystemExit(f'[interaction] fragment charges {[v["charge"] for v in fmeta.values()]} do not add up to '
                     f'{case["charge"]}')


def make_mol(idx_real, coords, charge, ghosts=()):
    atoms = [(syms[i], coords[i]) for i in idx_real] + [(f'GHOST-{syms[i]}', coords[i]) for i in ghosts]
    lib.num_threads(NTHREADS)
    return gto.M(atom=atoms, basis=basis, charge=charge, spin=0, unit='Angstrom', verbose=args.verbose,
                 max_memory=MAX_MEMORY, symmetry=False)


def d3(idx, coords):
    """D3(BJ) of the real atoms idx, from the case's functional (0 when the case has no dispersion)."""
    disp = case.get('disp') or 'none'
    if disp == 'none':
        return 0.0
    m = gto.M(atom=[(syms[i], coords[i]) for i in idx], basis='sto-3g', unit='Angstrom', verbose=0,
              spin=sum(_z(syms[i]) for i in idx) % 2)
    mf = _dft.RKS(m) if m.spin == 0 else _dft.UKS(m)
    mf.xc, mf.disp = case.get('xc') or XC, disp
    return float(mf.get_dispersion())


def energy(label, mol, solvent):
    t = time.time()
    mf = scf_for(case, mol, args, disp='none', solvent=solvent, conv_tol=1e-10)
    converge(mf, what=label)
    print(f'[interaction] {label:48s} E(SCF) = {mf.e_tot:.10f} Ha  ({time.time() - t:.0f} s)', flush=True)
    return float(mf.e_tot)


allidx = list(range(len(syms)))
A, B = names
res_env = {}
envs = [('pcm' if solv else 'gas', solv)]
if solv and not args.no_gas_bsse:
    envs.append(('gas_bsse_crosscheck', None))
for env, sv in envs:
    e = {}
    e['AB'] = energy(f'{env}: pair', make_mol(allidx, X, case['charge']), sv)
    for k in names:
        other = [i for i in allidx if i not in fr[k]]
        e[f'{k}|AB'] = energy(f'{env}: {k[:24]} in pair basis', make_mol(fr[k], X, fmeta[k]['charge'], other), sv)
        e[f'{k}|own'] = energy(f'{env}: {k[:24]} own basis', make_mol(fr[k], X, fmeta[k]['charge']), sv)
    res_env[env] = e
    if env == 'gas_bsse_crosscheck':
        break
    # relaxed fragments for the deformation energies
    given = [pathlib.Path(v).resolve() for v in args.relaxed_xyz.split(',')] if args.relaxed_xyz else [None, None]
    for k, pg in zip(names, given):
        rc = fmeta[k]['relaxed_case']
        p = pg or ((ROOT / 'molecules' / case['molecule'] / 'results' / rc / 'opt.xyz') if rc else None)
        if p is not None and p.exists():
            s, xr = read_xyz(p)
            if [a.capitalize() for a in s] != [syms[i].capitalize() for i in fr[k]]:
                raise SystemExit(f'[interaction] {p}: atoms do not match fragment {k}')
            full = np.zeros((len(syms), 3)); full[fr[k]] = xr
            e[f'{k}|relaxed'] = energy(f'{env}: {k[:24]} relaxed ({p.parent.name}/{p.name})',
                                       make_mol(fr[k], full, fmeta[k]['charge']), sv)
            fmeta[k]['relaxed_geometry'] = str(p.relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
            fmeta[k]['relaxed_sha256'] = sha256(p)
            e[f'{k}|relaxed_D3'] = d3(fr[k], full)
        else:
            fmeta[k]['relaxed_geometry'] = None
            print(f'[interaction] no relaxed geometry for {k} ({rc or "no matching case"}): deformation skipped', flush=True)

# dispersion (geometry only: same in every environment)
dAB = d3(allidx, X)
dX = {k: d3(fr[k], X) for k in names}
e = res_env['pcm' if solv else 'gas']
kc = KCAL_PER_HA
raw = e['AB'] - e[f'{A}|own'] - e[f'{B}|own']
cp = e['AB'] - e[f'{A}|AB'] - e[f'{B}|AB']
ddisp = dAB - dX[A] - dX[B]
out = {'E_int_raw_kcal': (raw + ddisp) * kc, 'E_int_CP_kcal': (cp + ddisp) * kc, 'BSSE_kcal': (cp - raw) * kc,
       'E_int_CP_electronic_kcal': cp * kc, 'E_int_D3_kcal': ddisp * kc}
defo = {}
for k in names:
    if f'{k}|relaxed' in e:
        defo[k] = ((e[f'{k}|own'] + dX[k]) - (e[f'{k}|relaxed'] + e[f'{k}|relaxed_D3'])) * kc
out['E_def_kcal'] = defo
if len(defo) == 2:
    out['E_bind_kcal'] = out['E_int_CP_kcal'] + sum(defo.values())
out_kj = {k.replace('_kcal', '_kJ'): (v * KJ_PER_HA / KCAL_PER_HA if not isinstance(v, dict)
                                      else {a: b * KJ_PER_HA / KCAL_PER_HA for a, b in v.items()})
          for k, v in out.items()}
cross = None
if 'gas_bsse_crosscheck' in res_env:
    g = res_env['gas_bsse_crosscheck']
    cross = {'BSSE_gas_same_geometry_kcal': ((g['AB'] - g[f'{A}|AB'] - g[f'{B}|AB'])
                                             - (g['AB'] - g[f'{A}|own'] - g[f'{B}|own'])) * kc,
             'E_int_CP_gas_same_geometry_kcal': (g['AB'] - g[f'{A}|AB'] - g[f'{B}|AB'] + ddisp) * kc}
res = {
    'provenance': provenance(case, args, xyz, t_start),
    'partition': partition, 'proton_transferred_in_geometry': transferred, 'fragments': fmeta,
    'energies_Ha': res_env, 'D3_Ha': {'AB': dAB, **dX},
    'results': out, 'results_kJ': out_kj, 'gas_crosscheck_for_pcm': cross,
    'definitions': 'E_int(CP)=E_AB-E_A(AB basis)-E_B(AB basis); BSSE=E_int(CP)-E_int(raw); '
                   'E_def(X)=E_X(cplx)-E_X(relaxed); E_bind=E_int(CP)+sum E_def; D3 added on real atoms (no BSSE)',
    'pcm_note': ('Approximate counterpoise in PCM: each species in its own cavity (fragments: real atoms only); '
                 'ghost functions extend into the solvent region where the partner was.') if solv else None,
}
name = out_name('interaction', tag)
save_json(outdir / name, res)
rows = [('E_int (raw)', f"{out['E_int_raw_kcal']:.2f}"), ('BSSE', f"{out['BSSE_kcal']:.2f}"),
        ('E_int (CP)', f"{out['E_int_CP_kcal']:.2f}"), ('  of which D3(BJ)', f"{out['E_int_D3_kcal']:.2f}")]
rows += [(f'E_def({k})', f'{v:.2f}') for k, v in defo.items()]
if 'E_bind_kcal' in out:
    rows.append(('E_bind = E_int(CP) + sum E_def', f"{out['E_bind_kcal']:.2f}"))
if cross:
    rows.append(('BSSE, gas phase at this geometry', f"{cross['BSSE_gas_same_geometry_kcal']:.2f}"))
(outdir / out_name('interaction', tag, '.md')).write_text(
    f"Interaction energy, {case['name']} ({environment_label(case)}), {partition} partition\n\n"
    + md_table(rows, ['Term', 'kcal/mol']))
if outdir == case['dir']:
    stage_note(case, f"interaction ({name}, {partition}): E_int(CP) {out['E_int_CP_kcal']:.2f} kcal/mol, BSSE "
                     f"{out['BSSE_kcal']:.2f}, E_bind {out.get('E_bind_kcal', float('nan')):.2f}")
print(open(outdir / out_name('interaction', tag, '.md')).read())
print(f"[interaction] {(time.time() - t_start) / 60:.1f} min -> {outdir / name}")
