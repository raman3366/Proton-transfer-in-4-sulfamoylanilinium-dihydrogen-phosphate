"""
Step 4b -- wavefunction analyses beyond dft/mwfn.py (which does the MEP surface, RDG/NCI and ELF/LOL planes),
run with Multiwfn 2026.9.20 as non-interactive stdin menu scripts on results/<case>/final.molden (dft/fmo.py).

  qtaim   QTAIM topology (main function 2): CPs from nuclei, atom-pair midpoints, triangles and pyramids, bond
          paths, Poincare-Hopf check, then every CP's properties (CPprop.txt).  For each non-covalent bond critical
          point: rho, Laplacian, G, V, H, -G/V, ellipticity, eta, and two H-bond energy estimates:
            E_HB = V/2                                  Espinosa, Molins, Lecomte, CPL 285 (1998) 170
            BE = -223.08 rho + 0.7423  (neutral fit)    Emamian, Lu, Kruse, Emamian, JCC 40 (2019) 2868,
            BE = -332.34 rho - 1.0661  (charged fit)    doi:10.1002/jcc.26068 (kcal/mol, rho in a.u.)
          The charged fit is used when either partner fragment carries a net charge.  NOTE: Emamian's "charged" set
          is a charged monomer bound to a neutral one; a cation...anion pair is outside both training sets, so
          both fits are reported and the numbers are indicative only.
          H-bonds are named from case['hbonds'] when the same H and donor/acceptor pair is found (either
          direction, so a transferred proton is still recognised), otherwise from the CIF labels.
  igmh    IGMH (20 -> 11; Lu & Chen, JCC 43 (2022) 539) between the molecular fragments: grid data
          (dg_inter and sign(lambda2)rho cubes for the isosurface figure), the scatter data, the integrated
          dg_inter, and the atom / atom-pair delta-g indices (atmdg.txt) mapped to atom labels.
  charges Hirshfeld, ADCH and CM5 (7 -> 1 / 11 / 16, built-in free-atom densities) from Multiwfn, plus NPA
          computed in-process from final.chk with PySCF's natural atomic orbitals (pyscf.lo.nao, Reed-Weinstock-
          Weinhold 1985).  Multiwfn has no NPA of its own without NBO, and we have no NBO licence.  PySCF's NAO
          partition is not identical to NBO 7's; differences of a few hundredths of e are expected.
  mep     MEP on the 0.001 a.u. density surface: runs dft/mwfn.py's MEP first if mep_surface.json is missing,
          then adds per-atom and per-fragment surface statistics (12 -> 0 -> 11 / 12) and collects the Vmin/Vmax
          list (with nearest atoms), sigma^2_tot, nu, Pi, MPI into one mep_summary.json.
  dos     TDOS and fragment PDOS (10; Hirshfeld orbital composition, because Mulliken/SCPA are unreliable with
          the diffuse ++ functions -- Multiwfn says so itself), Gaussian FWHM 0.3 eV; text data plus a per-MO
          fragment-composition table.

Fragments are NOT taken blindly from case['fragments']: they are the covalently connected components of the
final geometry, each mapped to the case fragment that owns its heavy atoms, with the charge corrected for any
H that moved (gas-phase sfa_pair transfers H1A to O5, turning the pair into neutral sulfanilamide + H3PO4).
Using the case lists there made IGMH report the H1A-O5 covalent bond as the dominant "interaction".

Output: results/<case>/mwfn_extra/ (inputs/*.in = the exact menu sequences, logs/*.out, CSV + JSON).
Usage: python dft/mwfn_extra.py <case> [qtaim igmh charges npa mep dos ...]      (default: all)
Environment: MWFN_BIN, MWFN_THREADS (default 4), IGMH_GRID (Multiwfn grid option, default 3 = ~1.7M points),
             DOS_FWHM_EV (0.3), DOS_RANGE_EV ("-30,5"), MWFN_KEEP_ALL=1 keeps dg.cub/dg_intra.cub/output.txt.
"""
import sys, os, re, subprocess, shutil, json, time, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *

MWFN = pathlib.Path(os.environ.get('MWFN_BIN') or ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI')
MWFN_HOME = MWFN.parent                                    # settings.ini lives next to the binary
NT = int(os.environ.get('MWFN_THREADS', '4'))
KCAL = 627.509474
EMAMIAN = {'neutral': (-223.08, 0.7423), 'charged': (-332.34, -1.0661)}      # BE = a*rho + b, kcal/mol
BONDI = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80, 'F': 1.47, 'Cl': 1.75}
ACCEPTORS = {'N', 'O', 'S', 'F', 'Cl', 'Br', 'P'}


class MwfnFailure(RuntimeError):
    pass


# ------------------------------------------------------------------ setup
def parse_molden(path):
    """Atoms (Angstrom) and MO energies/occupations from a PySCF molden file."""
    txt = path.read_text()
    m = re.search(r'\[Atoms\]\s*\(?(\w+)\)?\s*\n(.*?)\n\[', txt, re.S)
    unit = BOHR_A if m.group(1).lower().startswith('au') else 1.0
    syms, xyz = [], []
    for l in m.group(2).strip().splitlines():
        p = l.split()
        syms.append(re.sub(r'\d', '', p[0]).capitalize()); xyz.append([float(v) * unit for v in p[3:6]])
    ene = [float(v) for v in re.findall(r'Ene=\s*(\S+)', txt)]
    occ = [float(v) for v in re.findall(r'Occup=\s*(\S+)', txt)]
    return syms, np.array(xyz), np.array(ene), np.array(occ)


def case_fragment_charges(case):
    """Charge of each case fragment: explicit case['fragment_charges'], else from the name (cation/anion)."""
    fr = case.get('fragments') or {case['name']: list(range(len(case.get('cif_labels') or [])))}
    if case.get('fragment_charges'):
        q = dict(case['fragment_charges'])
    elif len(fr) == 1:
        q = {k: case['charge'] for k in fr}
    else:
        q = {k: (1 if re.search(r'\bcation\b', k.lower()) else -1 if re.search(r'\banion\b', k.lower()) else 0) for k in fr}
    if sum(q.values()) != case['charge']:
        raise SystemExit(f'fragment charges {q} do not add up to the case charge {case["charge"]}; '
                         f'add "fragment_charges" to the case file')
    return fr, q


def fragments_from_geometry(case, syms, xyz, lab):
    """Connected components of the final geometry, named and charged after the case fragments (see docstring)."""
    n = len(syms)
    bonds = bonds_from_geometry(syms, xyz)
    par = list(range(n))
    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]; i = par[i]
        return i
    for i, j, _ in bonds:
        par[find(i)] = find(j)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    comps = sorted(comps.values(), key=lambda c: min(c))
    cfr, cq = case_fragment_charges(case)
    if not case.get('fragments'):
        cfr = {case['name']: list(range(n))}
    owner = {i: k for k, atoms in cfr.items() for i in atoms}
    out = []
    for c in comps:
        heavy = [i for i in c if syms[i] != 'H']
        votes = {}
        for i in heavy or c:
            votes[owner.get(i)] = votes.get(owner.get(i), 0) + 1
        parent = max(votes, key=votes.get)
        orig = set(cfr.get(parent, []))
        gained = [i for i in c if i not in orig]; lost = [i for i in orig if i not in c]
        if any(syms[i] != 'H' for i in gained + lost) or parent is None:
            q, note = None, 'heavy-atom composition differs from the case fragment: charge unknown'
        else:
            q = cq[parent] + len(gained) - len(lost)
            note = ''.join([f' +{lab[i]}' for i in gained] + [f' -{lab[i]}' for i in lost]).strip()
        name = parent if not (gained or lost) else f'{parent} ({note})'
        out.append({'name': name, 'case_fragment': parent, 'atoms': sorted(c), 'charge': q,
                    'h_gained': [lab[i] for i in gained], 'h_lost': [lab[i] for i in lost]})
    return out, {(min(i, j), max(i, j)) for i, j, _ in bonds}


def ranges(idx):
    """0-based indices -> Multiwfn 1-based range string '1-15,17-20'."""
    idx = sorted(i + 1 for i in idx); out = []; s = p = idx[0]
    for i in idx[1:] + [None]:
        if i is not None and i == p + 1:
            p = i; continue
        out.append(f'{s}-{p}' if p > s else f'{s}'); s = p = i
    return ','.join(out)


def mwfn_home(esp_code=None):
    """Private Multiwfn settings dir (symlinks + own settings.ini) -- same scheme and reason as dft/mwfn.py."""
    home = OUT / 'mwfn_home'; home.mkdir(exist_ok=True)
    for f in MWFN_HOME.iterdir():
        if f.name != 'settings.ini' and not (home / f.name).exists():
            try:
                (home / f.name).symlink_to(f)
            except OSError:
                pass
    txt = (MWFN_HOME / 'settings.ini').read_text()
    # "-nt N" alone does not cap Multiwfn: it ran at 3.5-3.8 cores under -nt 2 (nested OpenMP regions; measured
    # 28 Sep 2026), which would break the burst thread budget.  nthreads is pinned here and OMP_THREAD_LIMIT in run().
    txt = re.sub(r'^(\s*nthreads=\s*)\d+', rf'\g<1>{NT}', txt, flags=re.M)
    if esp_code is not None:
        txt = re.sub(r'^(\s*iESPcode=\s*)\d+', rf'\g<1>{esp_code}', txt, flags=re.M)
    (home / 'settings.ini').write_text(txt)
    return home


def _unlimit_stack():
    # Multiwfn (Intel Fortran) keeps big arrays on the stack; 8 MB default segfaults for ~850 basis functions.
    import resource
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, OSError):
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))


def run(cmds, tag, expect=(), timeout=4 * 3600, esp_code=None):
    """Feed a menu sequence to Multiwfn in OUT/.  Success is judged on the expected output files, not on the exit
    code: Multiwfn exits non-zero after writing complete results and can segfault after writing nothing."""
    for f in expect:
        (OUT / f).unlink(missing_ok=True)                 # stale files must never count as success
    inp = OUT / 'inputs' / f'{tag}.in'; log = OUT / 'logs' / f'{tag}.out'
    inp.write_text('\n'.join(cmds) + '\n')
    env = dict(os.environ, Multiwfnpath=str(mwfn_home(esp_code)), OMP_STACKSIZE='200M', KMP_STACKSIZE='200M',
               OMP_THREAD_LIMIT=str(NT))            # hard cap incl. nested regions (see mwfn_home)
    t = time.time()
    with open(inp) as fi, open(log, 'w') as fo:
        rc = subprocess.run([str(MWFN), str(WFN), '-nt', str(NT), '-silent'], stdin=fi, stdout=fo,
                            stderr=subprocess.STDOUT, cwd=OUT, env=env, timeout=timeout,
                            preexec_fn=_unlimit_stack).returncode
    txt = log.read_text(errors='replace')
    wall = time.time() - t
    missing = [f for f in expect if not (OUT / f).exists() or (OUT / f).stat().st_size == 0]
    print(f'[mwfn_extra] {tag}: {wall:.0f} s (exit {rc})', flush=True)
    SUMMARY['runs'][tag] = {'menu_sequence': cmds, 'exit': rc, 'wall_s': round(wall, 1), 'threads': NT,
                            'esp_code': esp_code}
    v = re.search(r'Version (\S+)', txt)
    if v:
        SUMMARY['multiwfn_version'] = v.group(1)
    if missing:
        raise MwfnFailure(f'{tag}: exit {rc}' + (', SIGSEGV' if 'SIGSEGV' in txt else '') + f', missing/empty {missing}')
    return txt


def nearest_atoms(p, k=2):
    dist = np.linalg.norm(XYZ - p, axis=1); idx = np.argsort(dist)[:k]
    return ', '.join(f'{LAB[i]} ({dist[i]:.2f} A)' for i in idx)


def write_csv(path, rows, cols=None):
    if cols is None:
        cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, '')
                if isinstance(v, float):
                    v = f'{v:.6g}'
                v = str(v)
                vals.append(f'"{v}"' if ',' in v else v)
            f.write(','.join(vals) + '\n')


# ------------------------------------------------------------------ QTAIM
CP_FIELDS = {
    'rho': r'Density of all electrons:\s+(\S+)', 'G': r'Lagrangian kinetic energy G\(r\):\s+(\S+)',
    'K': r'Hamiltonian kinetic energy K\(r\):\s+(\S+)', 'V': r'Potential energy density V\(r\):\s+(\S+)',
    'H': r'Energy density E\(r\) or H\(r\):\s+(\S+)', 'lap': r'Laplacian of electron density:\s+(\S+)',
    'ELF': r'Electron localization function \(ELF\):\s+(\S+)', 'LOL': r'Localized orbital locator \(LOL\):\s+(\S+)',
    'RDG': r'Reduced density gradient \(RDG\):\s+(\S+)', 'sl2rho': r'Sign\(lambda2\)\*rho:\s+(\S+)',
    'dg_hirsh': r'Delta-g \(under Hirshfeld partition\):\s+(\S+)', 'ellipticity': r'Ellipticity of electron density:\s+(\S+)',
    'eta': r'eta index:\s+(\S+)', 'ESP_au': r'Total ESP:\s+(\S+) a\.u\.'}


def parse_cpprop(path):
    txt = path.read_text(errors='replace')
    heads = list(re.finditer(r'-{5,}\s+CP\s+(\d+),\s+Type\s+\((3,[+-]\d)\)\s+-{5,}', txt))
    cps = []
    for k, h in enumerate(heads):
        blk = txt[h.end(): heads[k + 1].start() if k + 1 < len(heads) else len(txt)]
        cp = {'index': int(h.group(1)), 'type': f'({h.group(2)})'}
        m = re.search(r'Position \(Angstrom\):\s+(\S+)\s+(\S+)\s+(\S+)', blk)
        cp['x'], cp['y'], cp['z'] = (float(v) for v in m.groups())
        m = re.search(r'Connected atoms:\s+(\d+)\(\s*\w+\s*\)\s+--\s+(\d+)\(', blk)
        if m:
            cp['atoms'] = (int(m.group(1)) - 1, int(m.group(2)) - 1)
        m = re.search(r'Corresponding nucleus:\s+(\d+)\(', blk)
        if m:
            cp['nucleus'] = int(m.group(1)) - 1
        for key, pat in CP_FIELDS.items():
            m = re.search(pat, blk)
            cp[key] = float(m.group(1)) if m else None
        m = re.search(r'Eigenvalues of Hessian:\s+(\S+)\s+(\S+)\s+(\S+)', blk)
        cp['lambda1'], cp['lambda2'], cp['lambda3'] = (float(v) for v in m.groups()) if m else (None,) * 3
        cps.append(cp)
    return cps


def bond_path_lengths(txt, cps):
    """Bond-path length (A) per BCP index, nucleus to nucleus.

    Multiwfn's path summary gives the two half-paths BCP -> (3,-3) CP.  With Gaussian basis sets the (3,-3) CP of an
    H atom sits 0.015-0.022 A off the nucleus, towards its bonding partner (measured on sfa_pair_gas_d3 in review,
    28 Sep 2026), so the plain sum of half-paths made every X-H bond path SHORTER than the X-H distance (N1-H1C
    0.994 vs 1.013 A, O4-H4A 0.942 vs 0.965 A), which is impossible for a path between nuclei.  The straight
    NCP -> nucleus segment is therefore added to each half-path (zero for heavy atoms, whose NCP is on the nucleus).
    """
    off = {cp['index']: float(np.linalg.norm(np.array([cp['x'], cp['y'], cp['z']]) - XYZ[cp['nucleus']]))
           for cp in cps if cp['type'] == '(3,-3)' and 'nucleus' in cp}
    plen = {}
    for m in re.finditer(r'CP:\s+(\d+) \(3,-1\) --> CP:\s+(\d+) \(3,-3\)\s+Length:\s+(\S+) Bohr', txt):
        b, n = int(m.group(1)), int(m.group(2))
        plen[b] = plen.get(b, 0.0) + float(m.group(3)) * BOHR_A + off.get(n, 0.0)
    return plen


def qtaim():
    # 2 topology; 2/3/4/5 search from nuclei, atom-pair midpoints, triangle and pyramid centres; 6 a sphere of
    # 200 starting points (radius 1.5 Bohr) around every nucleus in turn (10/11 set radius/points, -1 run, -9 back):
    # without it the midpoint search missed 3 of the 28 BCPs (short N-H bonds) on the X-ray ion pair and
    # Poincare-Hopf failed (27-21+2 = 8); with it all 28 are found (1), for ~3 s extra.
    # 8 bond paths; 0 summary (atom pairs of each BCP + Poincare-Hopf); 7 -> 0 all CP properties to CPprop.txt;
    # -4 -> 4 CPs.txt; -5 -> 4 paths.txt (for figures); -10 back to main menu; q.
    seq = ['2', '2', '3', '4', '5', '6', '10', '1.5', '11', '200', '-1', '-9',
           '8', '0', '7', '0', '-4', '4', '0', '-5', '4', '0', '-10', 'q']
    txt = run(seq, 'qtaim', expect=('CPprop.txt', 'CPs.txt'))
    cps = parse_cpprop(OUT / 'CPprop.txt')
    # Guard against a silently misaligned stdin script: every CP search prints "Totally find/found N new critical
    # points" (seen for options 2-5 in the builder's sandbox log; the sphere-search line was NOT inspected by the
    # reviewer).  Recorded so a burst log review can confirm that all 5 searches ran.
    searches = [int(v) for v in re.findall(r'Totally (?:find|found)\s+(\d+) new critical points', txt)]
    ph = re.search(r'Poincare-Hopf verification:\s*(\d+)\s*-\s*(\d+)\s*\+\s*(\d+)\s*-\s*(\d+)\s*=\s*(-?\d+)', txt)
    plen = bond_path_lengths(txt, cps)
    comp_of = {i: k for k, f in enumerate(FRAGS) for i in f['atoms']}
    hb_named = {tuple(v): k for k, v in (CASE.get('hbonds') or {}).items()}
    rows_all, rows_nc = [], []
    for cp in cps:
        row = {'cp': cp['index'], 'type': cp['type'], 'x_A': cp['x'], 'y_A': cp['y'], 'z_A': cp['z']}
        if 'nucleus' in cp:
            row['atoms'] = LAB[cp['nucleus']]
        if 'atoms' in cp:
            a, b = cp['atoms']
            row['atoms'] = f'{LAB[a]}-{LAB[b]}'
            row['d_AB_A'] = float(np.linalg.norm(XYZ[a] - XYZ[b]))
            row['bond_path_A'] = plen.get(cp['index'])
        row.update({k: cp[k] for k in ('rho', 'lap', 'G', 'V', 'H', 'ellipticity', 'eta', 'ELF', 'RDG', 'sl2rho',
                                       'dg_hirsh', 'lambda1', 'lambda2', 'lambda3')})
        rows_all.append(row)
        if cp['type'] != '(3,-1)' or 'atoms' not in cp:
            continue
        a, b = cp['atoms']
        if (min(a, b), max(a, b)) in BONDS:
            row['class'] = 'covalent'
            continue
        inter = comp_of[a] != comp_of[b]
        r = dict(row, cls='intermolecular' if inter else 'intramolecular',
                 fragments=f"{FRAGS[comp_of[a]]['name']} / {FRAGS[comp_of[b]]['name']}")
        row['class'] = r['cls']
        # H-bond geometry: H bonded to a donor D, the other BCP atom an acceptor
        h = a if SYMS[a] == 'H' else b if SYMS[b] == 'H' else None
        acc = (b if h == a else a) if h is not None else None
        name = f'{LAB[a]}...{LAB[b]}'; hb = False
        if h is not None and SYMS[acc] in ACCEPTORS:
            donors = [i if j == h else j for i, j in BONDS if h in (i, j)]
            if donors:
                D = min(donors, key=lambda j: np.linalg.norm(XYZ[j] - XYZ[h]))
                g = hbond_geometry(XYZ, D, h, acc)
                r.update({'donor': LAB[D], 'hydrogen': LAB[h], 'acceptor': LAB[acc], 'D-H_A': g['D-H_A'], 'H..A_A': g['H..A_A'],
                          'D..A_A': g['D..A_A'], 'DHA_deg': g['DHA_deg']})
                hb = g['H..A_A'] < BONDI['H'] + BONDI.get(SYMS[acc], 1.8) and g['DHA_deg'] >= 110
                name = f'{LAB[D]}-{LAB[h]}...{LAB[acc]}'
                for (cd, ch, ca), cname in hb_named.items():
                    if ch == h and {cd, ca} == {D, acc}:
                        name = cname if (cd, ca) == (D, acc) else f'{name} (reverse of case H-bond {cname}: proton transferred)'
                        r['case_hbond'] = cname
        qa, qb = FRAGS[comp_of[a]]['charge'], FRAGS[comp_of[b]]['charge']
        fit = 'charged' if (qa or qb) else 'neutral'
        r.update({'name': name, 'hbond_criteria_met': bool(hb), 'minusG_over_V': -cp['G'] / cp['V'],
                  'E_HB_Espinosa_kcal': cp['V'] / 2 * KCAL, 'E_HB_Espinosa_kJ': cp['V'] / 2 * KCAL * 4.184,
                  'BE_Emamian_neutral_kcal': EMAMIAN['neutral'][0] * cp['rho'] + EMAMIAN['neutral'][1],
                  'BE_Emamian_charged_kcal': EMAMIAN['charged'][0] * cp['rho'] + EMAMIAN['charged'][1],
                  'Emamian_fit_used': fit,
                  'nature_Rozas2000': ('weak (lap>0, H>0)' if cp['lap'] > 0 and cp['H'] > 0 else
                                       'medium (lap>0, H<0)' if cp['lap'] > 0 else 'strong (lap<0, H<0)')})
        r['BE_Emamian_kcal'] = r[f'BE_Emamian_{fit}_kcal']
        rows_nc.append(r)
    rows_nc.sort(key=lambda r: -r['rho'])
    write_csv(OUT / 'qtaim_cps.csv', rows_all, ['cp', 'type', 'atoms', 'class', 'd_AB_A', 'bond_path_A', 'rho', 'lap', 'G', 'V',
                                                'H', 'ellipticity', 'eta', 'ELF', 'RDG', 'sl2rho', 'dg_hirsh',
                                                'lambda1', 'lambda2', 'lambda3', 'x_A', 'y_A', 'z_A'])
    cols = ['name', 'cls', 'hbond_criteria_met', 'case_hbond', 'fragments', 'cp', 'donor', 'hydrogen', 'acceptor', 'D-H_A', 'H..A_A', 'D..A_A',
            'DHA_deg', 'd_AB_A', 'bond_path_A', 'rho', 'lap', 'G', 'V', 'H', 'minusG_over_V', 'ellipticity', 'eta',
            'E_HB_Espinosa_kcal', 'E_HB_Espinosa_kJ', 'BE_Emamian_neutral_kcal', 'BE_Emamian_charged_kcal',
            'Emamian_fit_used', 'BE_Emamian_kcal', 'nature_Rozas2000']
    write_csv(OUT / 'qtaim_noncovalent.csv', rows_nc, cols)
    missing = [k for k, v in (CASE.get('hbonds') or {}).items() if not any(r.get('case_hbond') == k for r in rows_nc)]
    counts = {t: sum(c['type'] == t for c in cps) for t in ('(3,-3)', '(3,-1)', '(3,+1)', '(3,+3)')}
    res = {'n_cp': len(cps), 'counts': counts,
           'cp_searches_logged': {'n': len(searches), 'new_cps_each': searches,
                                  'expected_n': 5, 'note': 'nuclei, midpoints, triangles, pyramids, spheres'},
           'poincare_hopf': {'n-b+r-c': int(ph.group(5)) if ph else None,
                             'satisfied': bool(ph and int(ph.group(5)) == 1)},
           'n_noncovalent_bcp': len(rows_nc), 'case_hbonds_without_bcp': missing,
           'n_bcp_without_bond_path': sum(c['type'] == '(3,-1)' and 'atoms' not in c for c in cps),
           'hbonds': [{k: r.get(k) for k in cols} for r in rows_nc if r['hbond_criteria_met']],
           'other_noncovalent': [{k: r.get(k) for k in cols} for r in rows_nc if not r['hbond_criteria_met']],
           'formulae': {'E_HB_Espinosa': 'V(r_BCP)/2 [Espinosa, Molins, Lecomte, CPL 285 (1998) 170, doi:10.1016/S0009-2614(98)00036-0]',
                        'BE_Emamian_neutral': 'BE/kcal mol-1 = -223.08 rho_BCP/a.u. + 0.7423 [doi:10.1002/jcc.26068]',
                        'BE_Emamian_charged': 'BE/kcal mol-1 = -332.34 rho_BCP/a.u. - 1.0661 [doi:10.1002/jcc.26068]',
                        'fit_choice': 'charged if either fragment has a net charge; ion pairs are outside both training sets',
                        'nature': 'Rozas, Alkorta, Elguero, JACS 122 (2000) 11154, doi:10.1021/ja0017864',
                        'hbond_criteria': 'BCP between H (bonded to D) and acceptor A, H...A < Bondi sum, D-H...A >= 110 deg'}}
    save_json(OUT / 'qtaim.json', res)
    SUMMARY['qtaim'] = {k: res[k] for k in ('n_cp', 'counts', 'cp_searches_logged', 'poincare_hopf', 'n_noncovalent_bcp',
                                            'case_hbonds_without_bcp', 'n_bcp_without_bond_path')}
    msg = (f"QTAIM: {len(cps)} CPs {counts}, Poincare-Hopf {'OK' if res['poincare_hopf']['satisfied'] else 'NOT satisfied'}; "
           + (f"{len(res['hbonds'])} H-bond BCP(s): " if res['hbonds'] else 'no H-bond BCPs')
           + '; '.join(f"{r['name']}: rho {r['rho']:.4f}, V/2 {r['E_HB_Espinosa_kcal']:.1f}, BE({r['Emamian_fit_used']}) "
                       f"{r['BE_Emamian_kcal']:.1f} kcal/mol" for r in rows_nc if r['hbond_criteria_met']))
    stage_note(CASE, msg); print('[mwfn_extra]', msg, flush=True)


# ------------------------------------------------------------------ IGMH
def igmh():
    if len(FRAGS) < 2:
        SUMMARY['igmh'] = 'skipped: one fragment only'; print('[mwfn_extra] IGMH skipped (one fragment)'); return
    grid = os.environ.get('IGMH_GRID', '3')
    seq = ['20', '11', str(len(FRAGS))] + [ranges(f['atoms']) for f in FRAGS] + [grid, '2', '3']
    two = len(FRAGS) == 2
    if two:
        seq += ['6', '2', 'y']                    # atom and atom-pair delta-g indices, high-quality integration grid
    seq += ['0', '0', 'q']
    expect = ('output.txt', 'sl2r.cub', 'dg_inter.cub') + (('atmdg.txt',) if two else ())
    txt = run(seq, 'igmh', expect=expect)
    g = lambda pat: (lambda m: float(m.group(1)) if m else None)(re.search(pat, txt))
    res = {'fragments': [{k: f[k] for k in ('name', 'charge')} | {'atoms': [LAB[i] for i in f['atoms']]} for f in FRAGS],
           'grid_option': grid, 'int_dg_au': g(r'Integral of delta-g over whole space:\s+(\S+)'),
           'int_dg_inter_au': g(r'Integral of delta-g_inter over whole space:\s+(\S+)'),
           'int_dg_intra_au': g(r'Integral of delta-g_intra over whole space:\s+(\S+)'),
           'figure_hint': 'starting point: isosurface of igmh_dg_inter.cube at ~0.005 a.u. coloured by igmh_sl2r.cube over -0.04..0.02 a.u.; tune per system'}
    m = re.search(r'Number of points in X,Y,Z is\s+(\d+)\s+(\d+)\s+(\d+)', txt)
    res['grid_points'] = [int(v) for v in m.groups()] if m else None
    # scatter: dg_inter, dg_intra, dg, sign(l2)rho per grid point -> compact npz (points carrying interfragment signal)
    dat = np.loadtxt(OUT / 'output.txt', dtype=np.float32)
    keep = dat[:, 0] > 1e-4
    np.savez_compressed(OUT / 'igmh_scatter.npz', dg_inter=dat[keep, 0], dg_intra=dat[keep, 1], dg=dat[keep, 2],
                        sl2rho=dat[keep, 3], n_total=len(dat))
    shutil.move(OUT / 'sl2r.cub', OUT / 'igmh_sl2r.cube'); shutil.move(OUT / 'dg_inter.cub', OUT / 'igmh_dg_inter.cube')
    if os.environ.get('MWFN_KEEP_ALL') != '1':
        for f in ('output.txt', 'dg.cub', 'dg_intra.cub'):
            (OUT / f).unlink(missing_ok=True)
    if two:
        t = (OUT / 'atmdg.txt').read_text()
        atoms, pairs = [], []
        for fr, blk in re.findall(r'fragment\s+(\d+) and percentage contributions\n(.*?)(?:\n\s*\n|\Z)', t, re.S):
            for i, v, p in re.findall(r'Atom\s+(\d+)\s*:\s+(\S+)\s+\(\s*(\S+) %', blk):
                atoms.append({'fragment': int(fr), 'atom': LAB[int(i) - 1], 'dg_index_au': float(v), 'pct': float(p)})
        blk = t.split('Atomic pair delta-g indices')[1]
        ibsiw = {}
        if (OUT / 'IBSIW.txt').exists():
            for i, j, v in re.findall(r'^\s+(\d+)\s+(\d+)\s*:\s+(\S+)', (OUT / 'IBSIW.txt').read_text().split('IBSIW index')[-1], re.M):
                ibsiw[(int(i), int(j))] = float(v)
        for i, j, v, p in re.findall(r'^\s+(\d+)\s+(\d+)\s*:\s+(\S+)\s+\(\s*(\S+) %', blk, re.M):
            a, b = int(i) - 1, int(j) - 1
            pairs.append({'pair': f'{LAB[a]}...{LAB[b]}', 'dist_A': float(np.linalg.norm(XYZ[a] - XYZ[b])),
                          'dg_pair_index_au': float(v), 'pct': float(p), 'IBSIW': ibsiw.get((int(i), int(j)))})
        write_csv(OUT / 'igmh_atoms.csv', atoms); write_csv(OUT / 'igmh_atom_pairs.csv', pairs)
        res['top_atom_pairs'] = pairs[:15]
    save_json(OUT / 'igmh.json', res)
    SUMMARY['igmh'] = {k: res[k] for k in ('int_dg_inter_au', 'grid_points')}
    top = ', '.join(f"{p['pair']} {p['pct']:.1f}%" for p in res.get('top_atom_pairs', [])[:4])
    stage_note(CASE, f"IGMH: integral dg_inter {res['int_dg_inter_au']} a.u.; top atom pairs {top}")
    print(f"[mwfn_extra] IGMH: int dg_inter {res['int_dg_inter_au']} a.u.; {top}", flush=True)


# ------------------------------------------------------------------ charges
def charges():
    # 7 population analysis; 1 Hirshfeld / 11 ADCH / 16 CM5, each "1" = built-in free-atom densities and
    # "y" = write the .chg file; 0 back; q.
    seq = ['7', '1', '1', 'y', '11', '1', 'y', '16', '1', 'y', '0', 'q']
    stem = WFN.stem
    txt = run(seq, 'charges', expect=(f'{stem}.chg',))
    blocks = re.findall(r'Final atomic charges:\s*\n((?:\s*Atom\s+\d+\(\s*\w+\s*\):\s+\S+\s*\n)+)', txt)
    if len(blocks) != 3:
        raise MwfnFailure(f'charges: expected 3 "Final atomic charges" blocks (Hirshfeld, ADCH, CM5), got {len(blocks)}')
    Q = {}
    for name, blk in zip(('Hirshfeld', 'ADCH', 'CM5'), blocks):
        Q[name] = [float(v) for v in re.findall(r'\):\s+(\S+)', blk)]
    (OUT / f'{stem}.chg').unlink(missing_ok=True)
    save_charges(Q)


def npa():
    """NPA charges from PySCF natural atomic orbitals (recipe of pyscf/examples/local_orb/02-pop_with_nao.py)."""
    from pyscf import scf
    from pyscf.lo import nao as nao_mod
    t = time.time()
    xyzfile = OUT / 'geom_for_npa.xyz'; write_xyz(xyzfile, SYMS, XYZ, 'from final.molden')
    mol = mol_for_case(CASE, xyzfile, verbose=0)
    chk = D / 'final.chk'
    if not chk.exists():
        raise MwfnFailure('npa: final.chk missing (written by dft/fmo.py)')
    # final.chk and final.molden are written by the same dft/fmo.py SCF, but a partial rerun (or a copy from
    # another run) could pair a chk with a different geometry; the electron-count test below cannot see that
    # (review 28 Sep 2026), so compare the chk's own molecule with the molden atoms first.
    from pyscf.lib import chkfile as _chk
    mchk = _chk.load_mol(str(chk))
    if mchk.nao != mol.nao or mchk.natm != mol.natm or \
            np.abs(mchk.atom_coords(unit='Angstrom') - XYZ).max() > 1e-4:
        raise MwfnFailure('npa: final.chk geometry/basis differs from final.molden; rerun dft/fmo.py')
    dat = scf.chkfile.load(str(chk), 'scf')
    mf = scf.RHF(mol) if np.asarray(dat['mo_occ']).ndim == 1 else scf.UHF(mol)
    mf.mo_coeff, mf.mo_occ = dat['mo_coeff'], dat['mo_occ']
    dm0 = mf.make_rdm1(); dm0 = dm0[0] + dm0[1] if np.asarray(dm0).ndim == 3 else dm0
    if abs(np.trace(dm0 @ mol.intor('int1e_ovlp')) - mol.nelectron) > 1e-6:
        raise MwfnFailure('npa: final.chk density does not integrate to the electron count (geometry/basis mismatch?)')
    C = nao_mod.nao(mol, mf)                        # AO -> NAO, orthonormal in the AO metric
    mo = np.linalg.solve(C, mf.mo_coeff) if mf.mo_occ.ndim == 1 else [np.linalg.solve(C, c) for c in mf.mo_coeff]
    dm = mf.make_rdm1(mo, mf.mo_occ)
    dm = dm[0] + dm[1] if np.asarray(dm).ndim == 3 else dm
    pop = np.diag(dm)
    q = [float(mol.atom_charge(a) - pop[s0:s1].sum()) for a, (_, _, s0, s1) in enumerate(mol.aoslice_by_atom())]
    os.remove(xyzfile)
    print(f'[mwfn_extra] NPA (PySCF NAO) in {time.time()-t:.0f} s; sum {sum(q):+.4f}', flush=True)
    save_charges({'NPA': q})


def save_charges(new):
    p = OUT / 'charges.json'
    old = json.loads(p.read_text()) if p.exists() else {'atomic': {}}
    old['atomic'].update(new)
    schemes = [s for s in ('NPA', 'Hirshfeld', 'ADCH', 'CM5') if s in old['atomic']]
    old['fragment_sums'] = {f['name']: {s: float(sum(old['atomic'][s][i] for i in f['atoms'])) for s in schemes}
                            | {'formal_charge': f['charge']} for f in FRAGS}
    old['sum'] = {s: float(sum(old['atomic'][s])) for s in schemes}
    old['notes'] = {'NPA': 'PySCF 2.14 pyscf.lo.nao (Reed, Weinstock, Weinhold 1985, doi:10.1063/1.449486); not NBO 7',
                    'Hirshfeld/ADCH/CM5': 'Multiwfn built-in sphericalized free-atom densities; ADCH doi:10.1142/S0219633612500113, '
                                          'CM5 doi:10.1021/ct200866d, Hirshfeld doi:10.1007/BF00549096'}
    save_json(p, old)
    rows = [{'index': i + 1, 'atom': LAB[i], 'element': SYMS[i],
             'fragment': next(f['name'] for f in FRAGS if i in f['atoms']), **{s: old['atomic'][s][i] for s in schemes}}
            for i in range(len(SYMS))]
    write_csv(OUT / 'charges.csv', rows)
    SUMMARY['charges'] = {'schemes': schemes, 'fragment_sums': old['fragment_sums']}
    stage_note(CASE, 'charges (' + '/'.join(schemes) + '): fragment sums ' +
               '; '.join(f"{k}: " + ', '.join(f"{s} {v[s]:+.3f}" for s in schemes) for k, v in old['fragment_sums'].items()))


# ------------------------------------------------------------------ MEP
def mep():
    if not (D / 'mep_surface.json').exists():
        print('[mwfn_extra] mep_surface.json missing -> running dft/mwfn.py <case> mep first', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'dft' / 'mwfn.py'), CASE['name'], 'mep'], check=True, cwd=ROOT,
                       env=dict(os.environ, MWFN_THREADS=str(NT), OMP_THREAD_LIMIT=str(NT)))
    # 12 surface analysis (defaults: rho = 0.001 a.u. isosurface, ESP mapped); 0 start; 11 per-atom table (n = no
    # locsurf.pqr); 12 <atoms> n per fragment; -1 -1 back to main; q.  (The vertices are already in the case
    # directory as vtx.txt, written by dft/mwfn.py.)
    seq = ['12', '0', '11', 'n']
    for f in FRAGS:
        seq += ['12', ranges(f['atoms']), 'n']
    seq += ['-1', '-1', 'q']
    last = None
    for code in (None, 3, 1):                    # LIBRETA fast ESP can segfault; same fallback as dft/mwfn.py
        try:
            txt = run(seq, 'mep_fragments', esp_code=code)
            nblk = txt.count('Properties on the surface of this fragment')
            if 'Atom#    All/Positive/Negative area' not in txt or nblk != len(FRAGS):
                raise MwfnFailure(f'mep_fragments: per-atom table or {len(FRAGS) - nblk} fragment block(s) missing')
            break
        except MwfnFailure as e:
            last = e
    else:
        raise last
    # per-atom tables (three blocks, in atom order)
    t1 = re.search(r'Atom#\s+All/Positive/Negative area \(Ang\^2\)\s+Minimal value\s+Maximal value\n(.*?)\n\s*Note', txt, re.S)
    t2 = re.search(r'Atom#\s+All/Positive/Negative average\s+All/Positive/Negative variance\n(.*?)\n\s*Note', txt, re.S)
    t3 = re.search(r'Atom#\s+Pi\s+nu\s+nu\*sigma\^2\n(.*?)\n\s*If outputting', txt, re.S)
    fl = lambda s: float('nan') if 'NaN' in s else float(s)
    rows = []
    for l1, l2, l3 in zip(t1.group(1).splitlines(), t2.group(1).splitlines(), t3.group(1).splitlines()):
        a, b, c = l1.split(), l2.split(), l3.split()
        i = int(a[0]) - 1
        rows.append({'atom': LAB[i], 'area_A2': fl(a[1]), 'pos_area_A2': fl(a[2]), 'neg_area_A2': fl(a[3]),
                     'min_kcal': fl(a[4]), 'max_kcal': fl(a[5]), 'avg_kcal': fl(b[1]), 'pos_avg_kcal': fl(b[2]),
                     'neg_avg_kcal': fl(b[3]), 'var_kcal2': fl(b[4]), 'Pi_kcal': fl(c[1]), 'nu': fl(c[2])})
    write_csv(OUT / 'mep_atoms.csv', rows)
    frs = []
    for f, blk in zip(FRAGS, re.findall(r'Properties on the surface of this fragment:(.*?)If outputting', txt, re.S)):
        g = lambda pat: (lambda m: fl(m.group(1)) if m else None)(re.search(pat, blk))
        frs.append({'fragment': f['name'], 'charge': f['charge'],
                    'min_kcal': g(r'Minimal value:\s+(\S+)'), 'max_kcal': g(r'Maximal value:\s+(\S+)'),
                    'area_A2': g(r'Overall surface area:.*?\(\s*(\S+) Angstrom'), 'pos_area_A2': g(r'Positive surface area:.*?\(\s*(\S+) Angstrom'),
                    'neg_area_A2': g(r'Negative surface area:.*?\(\s*(\S+) Angstrom'), 'avg_kcal': g(r'Overall average value:.*?\(\s*(\S+) kcal'),
                    'sigma2_tot_kcal2': g(r'Overall variance.*?\(\s*(\S+) \(kcal'), 'nu': g(r'Balance of charges \(nu\):\s+(\S+)'),
                    'Pi_kcal': g(r'Internal charge separation \(Pi\):.*?\(\s*(\S+) kcal')})
    ext = []
    if (D / 'mep_extrema.csv').exists():
        import csv
        ext = list(csv.DictReader(open(D / 'mep_extrema.csv')))
    summ = {'surface': 'electron density isosurface 0.001 a.u., ESP mapped (Multiwfn 12, defaults)',
            'whole': json.loads((D / 'mep_surface.json').read_text()),
            'extrema': [{k: (float(v) if k.endswith(('_au', '_eV', '_mol', '_A')) else v) for k, v in e.items()} for e in ext],
            'fragments': frs, 'per_atom_csv': 'mep_atoms.csv'}
    save_json(OUT / 'mep_summary.json', summ)
    SUMMARY['mep'] = {'n_extrema': len(ext), 'fragments': [{k: f[k] for k in ('fragment', 'min_kcal', 'max_kcal')} for f in frs]}
    print('[mwfn_extra] MEP per fragment: ' + '; '.join(f"{f['fragment']}: {f['min_kcal']:.1f}..{f['max_kcal']:.1f} kcal/mol"
                                                        for f in frs), flush=True)


# ------------------------------------------------------------------ DOS
def dos():
    fwhm = float(os.environ.get('DOS_FWHM_EV', '0.3'))
    lo_e, hi_e = (float(v) for v in os.environ.get('DOS_RANGE_EV', '-30,5').split(','))
    frags = FRAGS[:10]                             # Multiwfn allows 10 PDOS fragments
    # 10 DOS; 8 -> eV; 7 -> 3 Hirshfeld composition (set BEFORE fragments: changing the method clears them);
    # -1 fragments by atom list, 0 back; 2 energy range + label step; 3 FWHM; 0 draw TDOS+PDOS; 3 export
    # DOS_curve.txt / DOS_line.txt; 0 back; -10 main menu; q.
    seq = ['10', '8', '7', '3', '-1']
    for k, f in enumerate(frags):
        seq += [str(k + 1), ranges(f['atoms'])]
    seq += ['0', '2', f'{lo_e},{hi_e},5', '3', str(fwhm), '0', '3', '0', '-10', 'q']
    txt = run(seq, 'dos', expect=('DOS_curve.txt', 'DOS_line.txt'))
    names = [f['name'] for f in frags]
    cur = np.loadtxt(OUT / 'DOS_curve.txt')
    cols = ['E_eV', 'TDOS'] + [f'PDOS {n}' for n in names]
    with open(OUT / 'dos_curve.csv', 'w') as fo:
        fo.write(','.join(f'"{c}"' for c in cols) + '\n')
        for r in cur:
            fo.write(','.join(f'{v:.6g}' for v in r[:len(cols)]) + '\n')
    lines = np.loadtxt(OUT / 'DOS_line.txt')[1::3]               # middle row of each (0, height, 0) triple
    n_mo = len(ENE)
    real = np.abs(ENE) < 100                                    # fmo.py pads missing MOs at +1000 Ha
    if len(lines) != n_mo or np.abs(lines[real, 0] - ENE[real] * HARTREE_EV).max() > 0.01:
        raise MwfnFailure(f'dos: {len(lines)} level lines vs {n_mo} MOs in the molden, or energies do not match')
    rows = []
    for i in np.where(real)[0]:
        tot = lines[i, 1]
        rows.append({'mo': i + 1, 'E_eV': ENE[i] * HARTREE_EV, 'occ': OCC[i],
                     **{f'{n} %': 100 * lines[i, 2 + k] / tot if tot else float('nan') for k, n in enumerate(names)}})
    write_csv(OUT / 'mo_composition.csv', rows)
    homo = max(i for i in range(n_mo) if OCC[i] > 0.5); lumo = homo + 1
    centres = {names[int(k) - 1]: float(v) for k, v in re.findall(r'Center of PDOS\s+(\d+):\s+(\S+)', txt)}
    res = {'method': f'Hirshfeld orbital composition, Gaussian broadening FWHM {fwhm} eV, range {lo_e}..{hi_e} eV',
           'curve_units': 'Multiwfn DOS scale (arbitrary), eV axis', 'fragments': names,
           'tdos_centre_eV': float(re.search(r'Center of TDOS:\s+(\S+)', txt).group(1)) if 'Center of TDOS' in txt else None,
           'pdos_centres_eV': centres,
           'frontier': {lab: {'mo': i + 1, 'E_eV': float(ENE[i] * HARTREE_EV),
                              **{n: float(100 * lines[i, 2 + k] / lines[i, 1]) for k, n in enumerate(names)}}
                        for lab, i in (('HOMO-1', homo - 1), ('HOMO', homo), ('LUMO', lumo), ('LUMO+1', lumo + 1))}}
    save_json(OUT / 'dos.json', res)
    SUMMARY['dos'] = res['frontier']
    for f in ('DOS_curve.txt', 'DOS_line.txt'):         # converted above into dos_curve.csv / mo_composition.csv
        if os.environ.get('MWFN_KEEP_ALL') == '1':
            shutil.move(OUT / f, OUT / f'mwfn_{f}')
        else:
            (OUT / f).unlink()
    msg = 'DOS/PDOS (Hirshfeld): ' + '; '.join(f"{k} " + ', '.join(f'{n} {v[n]:.0f}%' for n in names) for k, v in res['frontier'].items() if k in ('HOMO', 'LUMO'))
    stage_note(CASE, msg); print('[mwfn_extra]', msg, flush=True)


# ------------------------------------------------------------------ main
if __name__ == '__main__':
    STEPS = {'qtaim': qtaim, 'igmh': igmh, 'charges': charges, 'npa': npa, 'mep': mep, 'dos': dos}
    CASE = load_case(sys.argv[1]); D = CASE['dir']
    what = [a for a in sys.argv[2:] if not a.startswith('-')] or list(STEPS)
    WFN = D / 'final.molden'
    if not WFN.exists():
        raise SystemExit(f'{WFN} missing: run the fmo stage first')
    OUT = D / 'mwfn_extra'
    for sub in ('', 'inputs', 'logs'):
        (OUT / sub).mkdir(exist_ok=True)
    SYMS, XYZ, ENE, OCC = parse_molden(WFN)
    if (D / 'opt.xyz').exists():
        s2, x2 = read_xyz(D / 'opt.xyz')
        if s2 != SYMS or np.abs(x2 - XYZ).max() > 1e-3:
            raise SystemExit('final.molden and opt.xyz disagree (atom order or geometry); refusing to mix them')
    LAB = CASE.get('cif_labels') or atom_labels(SYMS)
    FRAGS, BONDS = fragments_from_geometry(CASE, SYMS, XYZ, LAB)
    sp = OUT / 'summary.json'
    SUMMARY = json.loads(sp.read_text()) if sp.exists() else {}
    # summary.json and charges.json are MERGED across runs (e.g. 'npa' alone after a full run).  If final.molden
    # changed in between (re-optimised geometry, new fmo stage), the merge would silently mix two wavefunctions
    # in one table, so start both afresh whenever the wavefunction fingerprint differs (review 28 Sep 2026).
    import hashlib
    wfn_sha = hashlib.sha256(WFN.read_bytes()).hexdigest()
    if SUMMARY.get('wavefunction_sha256') != wfn_sha:     # also when summary.json is missing: charges.json unprovenanced
        if SUMMARY:
            print('[mwfn_extra] final.molden differs from the one in summary.json: discarding merged summary/charges',
                  flush=True)
        SUMMARY = {}
        (OUT / 'charges.json').unlink(missing_ok=True)
    SUMMARY['wavefunction_sha256'] = wfn_sha
    SUMMARY.setdefault('runs', {})
    SUMMARY.update({'case': CASE['name'], 'wavefunction': str(WFN.relative_to(ROOT)),
                    'fragments': [{k: f[k] for k in ('name', 'charge', 'h_gained', 'h_lost')} | {'atoms': ranges(f['atoms'])}
                                  for f in FRAGS]})
    needs_mwfn = [w for w in what if w != 'npa']
    if needs_mwfn and not MWFN.exists():
        raise SystemExit(f'Multiwfn binary not found at {MWFN} (x86 only); only "npa" runs without it')
    failed = []
    for w in what:
        t = time.time()
        try:
            STEPS[w]()
            SUMMARY.setdefault('status', {})[w] = {'ok': True, 'wall_s': round(time.time() - t, 1)}
        except Exception as e:                   # keep going: one failed analysis must not lose the others
            failed.append(w)
            SUMMARY.setdefault('status', {})[w] = {'ok': False, 'error': f'{type(e).__name__}: {e}'}
            stage_note(CASE, f'mwfn_extra {w} FAILED: {type(e).__name__}: {e}')
            print(f'[mwfn_extra] {w} FAILED: {type(e).__name__}: {e}', flush=True)
        save_json(sp, SUMMARY)
    sys.exit(1 if failed else 0)
