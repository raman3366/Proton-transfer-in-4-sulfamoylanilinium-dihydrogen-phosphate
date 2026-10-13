"""
Step 4 -- wavefunction analyses with Multiwfn on results/<case>/final.molden
  * MEP on the 0.001 a.u. electron-density surface (main function 12): surface extrema + statistics
      -> mep_extrema.csv, mep_surface.json, mep_range.json (colour range for rendering)
  * NCI / RDG (main function 20 -> 1, high-quality grid): rdg_scatter.txt, sl2r.cube (sign(lambda2)*rho), rdg.cube
  * ELF and LOL colour-filled plane data (main function 4 -> 9 / 10, plane through three ring atoms)
      -> elf_<plane>.npz, lol_<plane>.npz  (consumed by dft/plots.py)
Usage: python dft/mwfn.py <case> [mep|nci|elf|lol ...]
"""
import sys, os, re, subprocess, shutil, json, time, pathlib, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from common import *

MWFN = pathlib.Path(os.environ.get('MWFN_BIN') or next(
    (str(p) for p in (ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI',
                      ROOT / 'tools' / 'multiwfn' / 'Multiwfn',
                      ROOT / 'tools' / 'multiwfn-mac-build' / 'build' / 'multiwfn') if p.exists()),
    str(ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI')))
MWFN_HOME = ROOT / 'tools' / 'multiwfn'                     # holds settings.ini (nthreads, aug2D)
NT = int(os.environ.get('MWFN_THREADS', '8'))
KCAL = 627.509474

case = load_case(sys.argv[1]); d = case['dir']; (d / 'logs').mkdir(exist_ok=True)
what = [a for a in sys.argv[2:] if not a.startswith('-')] or ['mep', 'nci', 'elf', 'lol']
syms, xyz = read_xyz(d / 'opt.xyz'); lab = case.get('cif_labels') or atom_labels(syms)

def case_mwfn_home(esp_code=None):
    """Per-case Multiwfn settings directory.

    Multiwfn reads settings.ini from $Multiwfnpath.  Cases run concurrently, so editing the shared copy is
    unsafe; instead symlink everything and override only settings.ini.  esp_code sets iESPcode, which is the
    escape hatch for the LIBRETA segfault (see run()).
    """
    home = d / 'mwfn_home'
    home.mkdir(exist_ok=True)
    for f in MWFN_HOME.iterdir():
        link = home / f.name
        if f.name == 'settings.ini':
            continue
        if not link.exists():
            try:
                link.symlink_to(f)
            except OSError:
                pass
    txt = (MWFN_HOME / 'settings.ini').read_text()
    if esp_code is not None:
        txt = re.sub(r'^(\s*iESPcode=\s*)\d+', rf'\g<1>{esp_code}', txt, flags=re.M)
    (home / 'settings.ini').write_text(txt)
    return home


def run(cmds, tag, timeout=7200, esp_code=None, expect=()):
    """Feed a menu command sequence to Multiwfn (silent mode) with results/<case> as working directory.

    Multiwfn can die mid-analysis -- BTPC's MEP hit `forrtl: severe (174): SIGSEGV` inside the LIBRETA ESP
    routine after building the isosurface -- and the old version of this function ignored both the exit code
    and whether any output appeared, so the stage reported success while producing nothing.  Now a non-zero
    exit or a missing expected file raises.
    """
    inp = d / f'mwfn_{tag}.in'; out = d / 'logs' / f'mwfn_{tag}.out'
    inp.write_text('\n'.join(cmds) + '\n')
    home = case_mwfn_home(esp_code)
    # OMP_THREAD_LIMIT: "-nt N" alone lets Multiwfn's nested OpenMP regions run ~2N threads (measured 28 Sep 2026,
    # 3.5-3.8 cores under -nt 2), which breaks the burst thread budget.  A caller's own limit is kept.
    env = dict(os.environ, Multiwfnpath=str(home), OMP_STACKSIZE='200M', KMP_STACKSIZE='200M',
               OMP_THREAD_LIMIT=os.environ.get('OMP_THREAD_LIMIT', str(NT)))
    t = time.time()
    # Multiwfn is Fortran and puts large arrays on the STACK.  The default 8 MB process limit is fine up to
    # ~600 basis functions but segfaults on load for BTP (853 functions, 49 atoms) -- it dies in "Generating
    # density matrix" before any analysis starts.  OMP_STACKSIZE/KMP_STACKSIZE do not help: those size the
    # OpenMP worker stacks, not the main process.  Multiwfn's own manual asks for `ulimit -s unlimited`.
    def unlimit_stack():
        import resource
        try:
            resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        except (ValueError, OSError):
            soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))

    with open(inp) as fi, open(out, 'w') as fo:
        rc = subprocess.run([str(MWFN), 'final.molden', '-nt', str(NT), '-silent'], stdin=fi, stdout=fo,
                            stderr=subprocess.STDOUT, cwd=d, env=env, timeout=timeout,
                            preexec_fn=unlimit_stack).returncode
    txt = out.read_text(errors='replace')
    missing = [f for f in expect if not (d / f).exists() or (d / f).stat().st_size == 0]
    print(f'[mwfn] {tag}: {time.time()-t:.0f} s (exit {rc}'
          + (f', esp_code={esp_code}' if esp_code else '') + ')', flush=True)
    # Judge on OUTPUT, not exit code.  Multiwfn returns non-zero in situations where it has nevertheless
    # written complete results (the MEP for BTPC exited 59 after producing a valid surfanalysis.txt with
    # 10 minima and 17 maxima), and it can also die with SIGSEGV having produced nothing.  Requiring a zero
    # exit throws away good data; ignoring the exit entirely accepts crashes.  So: outputs present and
    # non-empty => success; otherwise fail with whatever diagnosis is available.
    if missing:
        raise MwfnFailure(f'{tag}: exit {rc}'
                          + (', SIGSEGV' if 'SIGSEGV' in txt else '')
                          + f', missing/empty {missing}')
    if rc != 0:
        print(f'[mwfn] {tag}: non-zero exit {rc} but all expected output present -- accepting', flush=True)
    return txt


class MwfnFailure(RuntimeError):
    pass

def nearest_atoms(p, k=2):
    dist = np.linalg.norm(xyz - p, axis=1); idx = np.argsort(dist)[:k]
    return ', '.join(f'{lab[i]} ({dist[i]:.2f} A)' for i in idx)

# ---------------- MEP on the molecular surface ----------------
def mep():
    # iESPcode 2 (LIBRETA fast) segfaults on some wavefunctions; fall back to 3 (LIBRETA slow) then
    # 1 (Multiwfn's original code), which are slower but robust.
    last = None
    for code in (None, 3, 1):
        try:
            txt = run(['12', '0', '1', '7', '-1', 'q'], 'mep', esp_code=code, expect=('surfanalysis.txt',))
            if code is not None:
                stage_note(case, f'MEP required iESPcode={code} (default ESP routine crashed)')
            return mep_post(txt)
        except MwfnFailure as e:
            last = e
            print(f'[mwfn] MEP attempt failed ({e}); retrying with a different ESP routine', flush=True)
    raise last

def mep_post(txt=None):
    if txt is None: txt = (d / 'logs' / 'mwfn_mep.out').read_text(errors='replace')
    surf = (d / 'surfanalysis.txt').read_text() if (d / 'surfanalysis.txt').exists() else txt
    ext = []
    for kind in ('minima', 'maxima'):
        m = re.search(rf'Number of surface {kind}:\s*(\d+)\s*\n(.*?)(?=\n\s*Number of surface|\n\s*=+|\Z)', surf, re.S)
        if not m: continue
        for line in m.group(2).splitlines():
            g = re.match(r'\s*\*?\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', line)
            if g:
                v = [float(x) for x in g.groups()[1:]]; p = np.array(v[3:6])
                ext.append({'type': kind[:-1] if kind != 'minima' else 'minimum', 'index': int(g.group(1)), 'au': v[0], 'eV': v[1], 'kcal_mol': v[2],
                            'x': v[3], 'y': v[4], 'z': v[5], 'nearest_atoms': nearest_atoms(p), 'global': '*' in line.split(str(g.group(1)))[0]})
    ext = [{**e, 'type': 'maximum' if e['type'] == 'maxim' else e['type']} for e in ext]
    with open(d / 'mep_extrema.csv', 'w') as f:
        f.write('type,index,global,ESP_au,ESP_eV,ESP_kcal_mol,x_A,y_A,z_A,nearest_atoms\n')
        for e in sorted(ext, key=lambda e: e['kcal_mol']):
            f.write(f"{e['type']},{e['index']},{'yes' if e['global'] else ''},{e['au']:.6f},{e['eV']:.4f},{e['kcal_mol']:.3f},{e['x']:.4f},{e['y']:.4f},{e['z']:.4f},\"{e['nearest_atoms']}\"\n")
    summ = {}
    for key, pat in {'global_min_au': r'Global surface minimum:\s*(-?\d+\.\d+)', 'global_max_au': r'Global surface maximum:\s*(-?\d+\.\d+)',
                     'volume_A3': r'Volume:.*?\(\s*(-?\d+\.\d+) Angstrom\^3', 'area_A2': r'Overall surface area:.*?\(\s*(-?\d+\.\d+) Angstrom\^2',
                     'positive_area_A2': r'Positive surface area:.*?\(\s*(-?\d+\.\d+) Angstrom', 'negative_area_A2': r'Negative surface area:.*?\(\s*(-?\d+\.\d+) Angstrom',
                     'overall_avg_kcal': r'Overall average value:.*?\(\s*(-?\d+\.\d+) kcal', 'positive_avg_kcal': r'Positive average value:.*?\(\s*(-?\d+\.\d+) kcal',
                     'negative_avg_kcal': r'Negative average value:.*?\(\s*(-?\d+\.\d+) kcal', 'sigma2_tot_kcal2': r'Overall variance.*?\(\s*(-?\d+\.\d+) \(kcal',
                     'nu': r'Balance of charges \(nu\):\s*(-?\d+\.\d+)', 'Pi_kcal': r'Internal charge separation \(Pi\):.*?\(\s*(-?\d+\.\d+) kcal',
                     'MPI_eV': r'Molecular polarity index \(MPI\):\s*(-?\d+\.\d+) eV', 'polar_area_pct': r'Polar surface area.*?\(\s*(-?\d+\.\d+) %'}.items():
        m = re.search(pat, txt); summ[key] = float(m.group(1)) if m else None
    summ['n_minima'] = sum(e['type'] == 'minimum' for e in ext); summ['n_maxima'] = sum(e['type'] == 'maximum' for e in ext)
    save_json(d / 'mep_surface.json', summ)
    if summ['global_min_au'] is not None:
        save_json(d / 'mep_range.json', {'vmin': summ['global_min_au'], 'vmax': summ['global_max_au']})
    stage_note(case, f"MEP surface: min {summ.get('global_min_au')} a.u., max {summ.get('global_max_au')} a.u.; {summ['n_minima']} minima, {summ['n_maxima']} maxima")
    print('[mwfn] MEP extrema:', len(ext), '| global min/max (a.u.):', summ['global_min_au'], summ['global_max_au'])

# ---------------- NCI / RDG ----------------
def nci():
    for f in ('output.txt', 'func1.cub', 'func2.cub'):
        (d / f).unlink(missing_ok=True)
    txt = run(['20', '1', '3', '2', '3', '0', '0', 'q'], 'nci', expect=('output.txt', 'func1.cub', 'func2.cub'))
    ok = (d / 'output.txt').exists() and (d / 'func1.cub').exists() and (d / 'func2.cub').exists()
    if not ok: raise RuntimeError('NCI outputs missing; see logs/mwfn_nci.out')
    shutil.move(d / 'output.txt', d / 'rdg_scatter.txt'); shutil.move(d / 'func1.cub', d / 'sl2r.cube'); shutil.move(d / 'func2.cub', d / 'rdg.cube')
    nci_post()

def nci_post():
    # columns of output.txt: x y z (Bohr), sign(lambda2)*rho, RDG.  Keep a compact copy for plotting.
    dat = np.loadtxt(d / 'rdg_scatter.txt', usecols=(3, 4)); x, y = dat[:, 0], dat[:, 1]
    keep = (x > -0.25) & (x < 0.06) & (y < 2.5); np.savez_compressed(d / 'rdg_scatter.npz', x=x[keep].astype(np.float32), y=y[keep].astype(np.float32), n_total=len(x))
    # spikes: local minima of RDG binned along sign(l2)rho -> table of interaction regions
    bins = np.arange(-0.20, 0.0501, 0.001); idx = np.digitize(x, bins); cand = []
    for b in range(1, len(bins)):
        sel = idx == b
        if sel.sum() > 20 and y[sel].min() < 0.3: cand.append(((bins[b-1] + bins[b]) / 2, float(y[sel].min()), int(sel.sum())))
    rows = []   # merge contiguous bins into one spike each (report the bin with the lowest RDG)
    for c, m, n in cand:
        if rows and abs(c - rows[-1][-1]) < 0.0015: 
            best = min([rows[-1][:3], (c, m, n)], key=lambda r: r[1]); rows[-1] = (*best, c)
        else: rows.append((c, m, n, c))
    rows = [r[:3] for r in rows]
    with open(d / 'nci_spikes.csv', 'w') as f:
        f.write('sign_l2_rho_au,min_RDG,n_points,type\n')
        for c, m, n in rows: f.write(f'{c:.4f},{m:.3f},{n},{"attractive (H-bond)" if c < -0.015 else ("van der Waals" if c < 0.01 else "steric repulsion")}\n')
    stage_note(case, f'NCI: {len(x)} grid points, {len(rows)} RDG spikes with RDG<0.3')
    print('[mwfn] NCI done; spikes:', rows[:12])

# ---------------- ELF / LOL plane maps ----------------
def plane(func_code, name):
    for pname, atoms in case.get('planes', {}).items():
        (d / 'plane.txt').unlink(missing_ok=True)
        a = ','.join(str(i + 1) for i in atoms)                    # Multiwfn atom indices are 1-based
        run(['4', str(func_code), '1', '300,300', '4', a, '-6', '-5', 'q'], f'{name}_{pname}', expect=('plane.txt',))
        if not (d / 'plane.txt').exists(): raise RuntimeError(f'plane.txt missing for {name} {pname}; see logs')
        dat = np.loadtxt(d / 'plane.txt'); P = dat[:, :3]; F = dat[:, 5]      # cols: X Y Z (Angstrom), graph-x, graph-y, value
        (d / 'plane.txt').rename(d / f'{name}_{pname}_plane.txt')
        # in-plane orthonormal axes from the three defining atoms
        p0, p1, p2 = xyz[atoms[0]], xyz[atoms[1]], xyz[atoms[2]]
        e1 = p1 - p0; e1 /= np.linalg.norm(e1); nrm = np.cross(e1, p2 - p0); nrm /= np.linalg.norm(nrm); e2 = np.cross(nrm, e1)
        c = P.mean(0); U = (P - c) @ e1; V = (P - c) @ e2
        n = int(round(np.sqrt(len(F)))); X = U.reshape(n, n); Y = V.reshape(n, n); Fz = F.reshape(n, n)
        if X[0, 0] == X[0, -1]: X, Y, Fz = X.T, Y.T, Fz.T
        atoms_xy = np.c_[(xyz - c) @ e1, (xyz - c) @ e2]; in_plane = np.abs((xyz - c) @ nrm) < 0.6
        np.savez(d / f'{name}_{pname}.npz', X=X, Y=Y, F=Fz, atoms_xy=atoms_xy, in_plane=in_plane, normal=nrm, center=c)
        print(f'[mwfn] {name} plane {pname}: {n}x{n} points, range {Fz.min():.3f}..{Fz.max():.3f}')
    stage_note(case, f'{name.upper()} plane data written for {list(case.get("planes", {}))}')

if __name__ == '__main__':
    if not MWFN.exists(): raise SystemExit(f'Multiwfn binary not found at {MWFN}')
    for w in what:
        {'mep': mep, 'meppost': mep_post, 'nci': nci, 'ncipost': nci_post, 'elf': lambda: plane(9, 'elf'), 'lol': lambda: plane(10, 'lol')}[w]()
