"""Gate G-R2, step 2: the single points of the new geometries (md/N2.md, section 3.5).  Runs on the rented machine.

    QM = /root/cqc-sfa-tools/env/bin/python      (PySCF 2.14.0; nothing of the MD environment is needed)
    QM md/qm/n2_batch.py run    [--workers 4] [--threads 4] [--deadline 18:30] [--list main|all] [--only ID,ID]
    QM md/qm/n2_batch.py point  ID            one single point (what `run` starts for every point)
    QM md/qm/n2_batch.py status               counts; no energy is printed
    QM md/qm/n2_batch.py d3check              the dispersion routine below against md/qm/gr/d3.json (38 geometries)
    QM md/qm/n2_batch.py checkx [--threads N] point X of gate G-R with THIS program; must give the energy of
                                              md/qm/gr/qm/X.json within 1e-7 Hartree -> md/qm/n2/check_X.json
    QM md/qm/n2_batch.py complete             exit code 0 when `run` has ended in state "finished"
    options of every command: --out DIR (default md/qm/n2: results qm/<id>.json), --bulk DIR (default
    /root/cqc-sfa-md/qm/n2: logs, checkpoints, status.json, batch.log, STOP), --override JSON (a test level,
    for the rehearsal only: for example '{"basis": "sto-3g", "grid": 1}'; such a result is marked and no
    program of the gate accepts it).

WHAT ONE POINT IS.  The settings of gate G-R, read from md/qm/config_prod.json as md/qm/gr_queue.py reads them
(Config.settings and run_point of the committed version, sha256 0fa1d56a...; the code is REPEATED here, not
imported, because gr_queue.py is being changed while this is written): B3LYP (the B3LYP of Gaussian, asserted
by pyscf_calc.assert_b3lyp_is_gaussian), 6-31G* with pure d functions, density fitting with
def2-universal-jkfit, grid level 3, convergence 1e-9 Hartree, gradient convergence max(1e-6, 0.1 sqrt(1e-9)),
at most 60 cycles and then the second-order solver, atomic guess (or the point's own checkpoint after a stop),
charge +1, singlet.  md/qm/pyscf_calc.py and md/qm/run_qm.py are imported (build_mol, Settings, set_grid,
has_scf), unchanged.  D3(BJ): class D3 below, the two-body term with the parameters of gate G-R, a verbatim
copy of the class of the committed md/qm/gr_queue.py, with md/qm/gr/d3_reference.json.

ONE JSON PER POINT, md/qm/n2/qm/<id>.json: e_tot (Hartree, B3LYP), d3_hartree, e_qm_hartree (the sum),
converged, scf_cycles, second_order_solver, wall_seconds, threads, the set ("fit" or "test") of
md/qm/n2/split.json, the checksum of the geometry file, settings, versions, host.  A file with
`test_level: true` (an --override) is never taken for a result.

THE RUN.  Points in the order of md/qm/n2/split.json (the main list of 144; with --list all the reserve of 48
after it), `workers` at a time, each a process of its own with `threads` threads.  Restartable: a point with a
converged result for the same geometry is skipped; a stopped SCF goes on from its own orbitals.  No new point
is started after --deadline (IST, today) or when <bulk>/STOP exists; what runs is left to finish (STOP: the SCF
stops at its next cycle and keeps its orbitals).  A point that does not converge is tried once more and then
left as failed, with its file.  The log <bulk>/batch.log has one line per point WITHOUT its energy, so that
nobody reads test energies by reading a log.
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n2_common as nc                                                     # noqa: E402

QM_CONFIG = HERE / 'config_prod.json'
D3_REF = HERE / 'gr' / 'd3_reference.json'
D3_Z = dict(H=1, C=6, N=7, O=8, S=16, Zn=30)
D3_FD_STEP = 1.0e-3
BOHR = 0.52917721092                                                       # as md/qm/gr/d3.json (angstrom_per_bohr)
COPIED_FROM = dict(file='md/qm/gr_queue.py', version='git HEAD when this was written',
                   sha256='0fa1d56a8e07b5e442262e239b81072c9a7b4213cef17db198a920c41fcf7ff8',
                   what='class D3 verbatim; Config.settings and run_point repeated')


def read_json(path, default=None):
    try:
        return json.loads(pathlib.Path(path).read_text())
    except Exception:                                                      # noqa: BLE001
        return default


class D3:
    """B3LYP-D3(BJ), two-body, for a fixed list of atomic numbers; coordinates in bohr."""

    def __init__(self, z, ref):
        import numpy as np
        self.np = np
        ref = read_json(ref) if not isinstance(ref, dict) else ref
        self.z = [int(v) for v in z]
        n = len(self.z)
        self.par = ref['b3lyp_bj']
        self.rc = 4.0 / 3.0 * np.array([ref['rcov_angstrom'][str(v)] for v in self.z]) / ref['angstrom_per_bohr']
        self.rc = self.rc[:, None] + self.rc[None, :]
        self.kcn, self.kw = ref['k_cn'], ref['k_weight']
        self.cnref = np.full((n, 7), 1.0e9)                                # a padded reference gets weight zero
        self.c6ref = np.zeros((n, n, 7, 7))
        for i, zi in enumerate(self.z):
            self.cnref[i, :ref['n_ref'][str(zi)]] = ref['cn_ref'][str(zi)]
            for j, zj in enumerate(self.z):
                m = np.array(ref['c6_ref'][f'{zi}-{zj}'])
                self.c6ref[i, j, :m.shape[0], :m.shape[1]] = m
        q = np.array([np.sqrt(0.5 * ref['r4_over_r2'][str(v)] * np.sqrt(v)) for v in self.z])
        self.qq = 3.0 * q[:, None] * q[None, :]
        self.r0 = self.par['a1'] * np.sqrt(self.qq) + self.par['a2']
        self.iu = np.triu_indices(n, 1)

    def c6(self, xyz):
        np = self.np
        d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
        np.fill_diagonal(d, 1.0)
        count = 1.0 / (1.0 + np.exp(-self.kcn * (self.rc / d - 1.0)))
        np.fill_diagonal(count, 0.0)
        cn = count.sum(axis=1)
        w = np.exp(-self.kw * (cn[:, None] - self.cnref) ** 2)
        norm = w.sum(axis=1)
        assert (norm > 1e-200).all(), 'a coordination number far from every reference'
        w /= norm[:, None]
        return np.einsum('ia,jb,ijab->ij', w, w, self.c6ref, optimize=True), d, cn

    def energy(self, xyz, c6=None):
        """Hartree.  With `c6` (a matrix from an earlier geometry) the coefficients are not recomputed."""
        np = self.np
        c, d, _ = self.c6(xyz)
        if c6 is not None:
            c = c6
        r = d[self.iu]
        p = self.par
        e = c[self.iu] * (p['s6'] / (r ** 6 + self.r0[self.iu] ** 6) + p['s8'] * self.qq[self.iu] / (r ** 8 + self.r0[self.iu] ** 8))
        return -float(e.sum())

    def gradient(self, xyz, atoms=None):
        """Hartree/bohr by central differences (step D3_FD_STEP); zero for atoms not in `atoms`."""
        np = self.np
        g = np.zeros_like(xyz)
        for i in (range(len(xyz)) if atoms is None else atoms):
            for k in range(3):
                x = xyz.copy()
                x[i, k] += D3_FD_STEP
                ep = self.energy(x)
                x[i, k] -= 2 * D3_FD_STEP
                g[i, k] = (ep - self.energy(x)) / (2 * D3_FD_STEP)
        return g


def settings(override=None):
    """As Config.settings of md/qm/gr_queue.py."""
    s2 = json.loads(QM_CONFIG.read_text())
    m = dict(s2['method'])
    st = dict(xc=m['xc'], basis=m['basis'], cart=m['cart'], density_fit=m['density_fit'],
              auxbasis=m.get('auxbasis', 'def2-universal-jkfit'), grid=s2['opt']['grid'],
              conv_tol=s2['opt']['conv_tol'], scf_max_cycle=int(s2.get('scf_max_cycle', 100)),
              max_memory_mb=int(s2['max_memory_mb']), source=str(QM_CONFIG), source_sha256=nc.sha(QM_CONFIG),
              overridden=sorted(override or {}))
    st.update(override or {})
    return st


def d3_energy(sym, xyz_ang):
    import numpy as np
    ref = read_json(D3_REF)
    d = D3([D3_Z[s] for s in sym], ref)
    return d.energy(np.asarray(xyz_ang) / ref['angstrom_per_bohr'])


class StopRequested(Exception):
    pass


def single_point(pid, xyz_file, out_file, bulk, st, extra):
    """As run_point of md/qm/gr_queue.py with the atomic guess; returns 0 (converged), 3 (stopped), 4 (not converged)."""
    import numpy as np
    import pyscf_calc as pc
    import run_qm as rq
    from pyscf import dft, lib
    t_start = time.time()
    for d in ('chk', 'logs'):
        (bulk / d).mkdir(parents=True, exist_ok=True)
    lib.param.TMPDIR = str(bulk)
    sym, xyz = nc.read_xyz(xyz_file)
    s = pc.Settings(xc=st['xc'], basis=st['basis'], cart=st['cart'], density_fit=st['density_fit'],
                    auxbasis=st['auxbasis'], conv_tol=st['conv_tol'], max_memory_mb=st['max_memory_mb'], verbose=4)
    functional = pc.assert_b3lyp_is_gaussian() if st['xc'].lower() == 'b3lyp' else f"test level {st['xc']}"
    log = bulk / 'logs' / f'{pid}.pyscf.log'
    mol = pc.build_mol(sym, xyz, 1, 1, s, output=str(log))
    mf = dft.RKS(mol)
    mf.xc = st['xc']
    rq.set_grid(mf, st['grid'] if isinstance(st['grid'], int) else tuple(st['grid']))
    if st['density_fit']:
        mf = mf.density_fit(auxbasis=st['auxbasis'])
    mf.conv_tol = st['conv_tol']
    mf.conv_tol_grad = max(1e-6, np.sqrt(st['conv_tol']) * 0.1)
    mf.max_cycle = st['scf_max_cycle']
    chk = bulk / 'chk' / f'{pid}.chk'
    own = chk.with_suffix('.part.chk')
    guess, dm0 = dict(kind='atomic (PySCF default)'), None
    if rq.has_scf(own):
        dm0 = mf.from_chk(str(own))
        guess = dict(kind='own checkpoint of a stopped run')
    mf.chkfile = str(own)
    cycles = []
    stop_file = bulk / 'STOP'

    def callback(env):
        cycles.append(dict(cycle=int(env['cycle']), e=float(env['e_tot']), t=round(time.time() - t_start, 1)))
        if stop_file.exists():
            raise StopRequested()
    mf.callback = callback
    started = nc.now()
    second_order = False
    try:
        mf.kernel(dm0=dm0)
        if not mf.converged:
            second_order = True
            mf2 = mf.newton()
            mf2.callback = callback
            mf2.kernel(dm0=mf.make_rdm1())
            if mf2.converged:
                mf.mo_coeff, mf.mo_energy, mf.mo_occ, mf.e_tot, mf.converged = (
                    mf2.mo_coeff, mf2.mo_energy, mf2.mo_occ, mf2.e_tot, True)
    except StopRequested:
        mol.stdout.flush()
        print(f'{pid}: stop requested after {len(cycles)} cycles; orbitals kept in {own.name}', flush=True)
        return 3
    mol.stdout.flush()
    e_d3 = d3_energy(sym, xyz)
    res = dict(id=pid, e_tot=float(mf.e_tot), d3_hartree=float(e_d3), e_qm_hartree=float(mf.e_tot) + float(e_d3),
               converged=bool(mf.converged), scf_cycles=len(cycles), second_order_solver=second_order,
               last_energy_change=(cycles[-1]['e'] - cycles[-2]['e']) if len(cycles) > 1 else None,
               homo_lumo_gap_hartree=round(float(mf.mo_energy[mf.mo_occ == 0].min() - mf.mo_energy[mf.mo_occ > 0].max()), 6),
               natm=mol.natm, nelectron=int(mol.nelectron), nao=int(mol.nao),
               naux=int(mf.with_df.get_naoaux()) if st['density_fit'] else None,
               grid_points=int(mf.grids.weights.size), guess=guess, settings=st, functional=functional,
               test_level=bool(st['overridden']), threads=lib.num_threads(), started=started, finished=nc.now(),
               wall_seconds=round(time.time() - t_start, 1), xyz=str(xyz_file), xyz_sha256=nc.sha(xyz_file),
               d3=dict(routine='class D3 of md/qm/n2_batch.py', reference_data_sha256=nc.sha(D3_REF), copied_from=COPIED_FROM),
               pyscf_log=str(log), versions=pc.versions(), host=os.uname().nodename, **extra)
    if mf.converged:
        os.replace(own, chk)
    nc.write_json(out_file, res)
    print(f"{pid}: converged {mf.converged}, {len(cycles)} cycles, {res['wall_seconds']:.0f} s", flush=True)
    return 0 if mf.converged else 4


# ------------------------------------------------------------------------------------------------- commands
def paths(args):
    out = pathlib.Path(args.out).resolve() if args.out else nc.N2
    bulk = pathlib.Path(args.bulk) if args.bulk else pathlib.Path('/root/cqc-sfa-md/qm/n2')
    return out, bulk


def split_of():
    return nc.split_record()


def finished(out, pid, split, test_level):
    r = read_json(out / 'qm' / f'{pid}.json')
    return bool(r and r.get('converged') and r.get('xyz_sha256') == split['xyz_sha256'][pid]
                and bool(r.get('test_level')) == test_level)


def cmd_point(args):
    out, bulk = paths(args)
    split = split_of()
    pid = args.id
    assert pid in split['set'], f'{pid} is not a point of md/qm/n2/split.json'
    xyz = nc.N2 / 'geom' / f'{pid}.xyz'
    assert nc.sha(xyz) == split['xyz_sha256'][pid], f'{xyz} is not the geometry the split was drawn for'
    override = json.loads(args.override) if args.override else None
    return single_point(pid, xyz, out / 'qm' / f'{pid}.json', bulk, settings(override),
                        dict(set=split['set'][pid], order=split['order'].index(pid) + 1))


def deadline_epoch(text):
    """'18:30' or '6:30 PM': that time of today in IST (UTC+5:30), as seconds since the epoch."""
    if not text:
        return None
    t = None
    for fmt in ('%H:%M', '%I:%M %p', '%I:%M%p'):
        try:
            t = time.strptime(text.strip().upper(), fmt)
            break
        except ValueError:
            pass
    assert t is not None, f'--deadline {text!r}: give 18:30 or "6:30 PM" (IST)'
    ist = time.gmtime(time.time() + 19800)
    midnight = time.time() + 19800 - (ist.tm_hour * 3600 + ist.tm_min * 60 + ist.tm_sec)
    return midnight + t.tm_hour * 3600 + t.tm_min * 60 - 19800


def cmd_run(args):
    out, bulk = paths(args)
    bulk.mkdir(parents=True, exist_ok=True)
    (out / 'qm').mkdir(parents=True, exist_ok=True)
    split = split_of()
    test_level = bool(args.override)
    order = split['order'] if args.list == 'all' else split['order'][:split['n_main']]
    if args.only:
        order = [p for p in order if p in args.only.split(',')]
    stop_file = bulk / 'STOP'
    if stop_file.exists():
        stop_file.unlink()                                                # as the tasks of gate G-R: a new start removes it
    dl = deadline_epoch(args.deadline)
    pending = [p for p in order if not finished(out, p, split, test_level)]
    tries, procs, failed, t0, done_now = {}, {}, [], time.time(), []
    logf = open(bulk / 'batch.log', 'a')

    def log(text):
        logf.write(f'{nc.now()}  {text}\n')
        logf.flush()
        print(text, flush=True)

    def save(state):
        n_done = sum(finished(out, p, split, test_level) for p in order)
        rate = (time.time() - t0) / max(1, len(done_now))
        nc.write_json(bulk / 'status.json', dict(
            state=state, host=os.uname().nodename, pid=os.getpid(), written=nc.now(), workers=args.workers,
            threads=args.threads, list=args.list, deadline=args.deadline, test_level=test_level,
            n_points=len(order), n_finished=n_done, running=sorted(procs), failed=failed,
            finished_in_this_run=len(done_now), seconds_per_point_at_this_rate=round(rate, 1) if done_now else None,
            minutes_left_estimate=round((len(pending) + len(procs)) * rate / 60.0, 1) if done_now else None,
            mean_wall_seconds=round(sum(done_now) / len(done_now), 1) if done_now else None))
    import signal

    def on_term(signum, frame):
        # the task runner ends a task with SIGTERM: the running SCFs stop at their next cycle and keep their orbitals
        stop_file.write_text(f'{nc.now()} written by n2_batch.py run on signal {signum}\n')
    signal.signal(signal.SIGTERM, on_term)
    signal.signal(signal.SIGINT, on_term)
    log(f'run: {len(order)} points, {len(pending)} to do, {args.workers} at a time with {args.threads} threads, '
        f'deadline {args.deadline or "none"} IST, list {args.list}' + (', TEST LEVEL' if test_level else ''))
    save('running')
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS=str(args.threads),
               MKL_NUM_THREADS=str(args.threads))
    while pending or procs:
        for pid, (p, ts) in list(procs.items()):
            rc = p.poll()
            if rc is None:
                continue
            del procs[pid]
            ok = finished(out, pid, split, test_level)
            if ok:
                done_now.append(time.time() - ts)
            elif rc != 3 and tries[pid] < 2:
                pending.insert(0, pid)
            elif rc != 3:
                failed.append(pid)
            log(f"{pid} ({split['set'][pid]}): exit code {rc}, {'converged' if ok else 'NOT finished'}, "
                f'{time.time() - ts:.0f} s, attempt {tries[pid]}')
            save('running')
        late = dl is not None and time.time() > dl
        halt = stop_file.exists() or late
        while pending and len(procs) < args.workers and not halt:
            pid = pending.pop(0)
            tries[pid] = tries.get(pid, 0) + 1
            cmd = [sys.executable, str(pathlib.Path(__file__).resolve()), 'point', pid, '--out', str(out), '--bulk', str(bulk)]
            if args.override:
                cmd += ['--override', args.override]
            lf = open(bulk / 'logs' / f'{pid}.out', 'a') if (bulk / 'logs').exists() or not (bulk / 'logs').mkdir(parents=True) else None
            procs[pid] = (subprocess.Popen(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT), time.time())
            save('running')
        if halt and not procs:
            break
        time.sleep(2 if test_level else 5)
    n_done = sum(finished(out, p, split, test_level) for p in order)
    state = 'finished' if (n_done + len(failed) == len(order) or (dl is not None and time.time() > dl)) else 'stopped'
    save(state)
    log(f'run ended: state {state}, {n_done} of {len(order)} points finished, failed {failed}')
    return 0 if state == 'finished' else 3


def cmd_status(args):
    out, bulk = paths(args)
    split = split_of()
    st = read_json(bulk / 'status.json', {})
    main = split['order'][:split['n_main']]
    f_main = [p for p in main if finished(out, p, split, False)]
    f_all = [p for p in split['order'] if finished(out, p, split, False)]
    print(json.dumps(dict(status=st, finished_main=len(f_main), of_main=len(main), finished_all=len(f_all),
                          of_all=len(split['order']),
                          finished_fit=sum(split['set'][p] == 'fit' for p in f_all),
                          finished_test=sum(split['set'][p] == 'test' for p in f_all)), indent=1))
    return 0


def cmd_complete(args):
    out, bulk = paths(args)
    st = read_json(bulk / 'status.json', {})
    return 0 if st.get('state') == 'finished' and not st.get('test_level') else 1


def cmd_d3check(args):
    import numpy as np
    out, bulk = paths(args)
    rec = read_json(HERE / 'gr' / 'd3.json')
    worst, n = 0.0, 0
    for pid, v in rec['points'].items():
        f = HERE / 'gr' / 'geom' / f'{pid}.xyz'
        if nc.sha(f) != v['xyz_sha256']:
            continue
        sym, xyz = nc.read_xyz(f)
        worst = max(worst, abs(d3_energy(sym, xyz) - v['e_d3bj']))
        n += 1
    ok = bool(n == len(rec['points']) and worst < 1e-10)
    nc.write_json(out / 'check_d3.json', dict(
        what='class D3 of md/qm/n2_batch.py against md/qm/gr/d3.json (the dispersion energies of gate G-R)',
        geometries=n, of=len(rec['points']), largest_abs_difference_hartree=worst, tolerance_hartree=1e-10,
        agree=ok, host=os.uname().nodename, written=nc.now(), numpy=np.__version__, copied_from=COPIED_FROM))
    print(f'D3: {n} geometries, largest difference {worst:.2e} Hartree, agree {ok}')
    return 0 if ok else 1


def cmd_checkx(args):
    out, bulk = paths(args)
    override = json.loads(args.override) if args.override else None
    xyz = HERE / 'gr' / 'geom' / 'X.xyz'
    tmp = bulk / 'check_X'
    rc = single_point('X_check', xyz, tmp / 'X_check.json', tmp, settings(override), dict(set='check'))
    r = read_json(tmp / 'X_check.json')
    gate = read_json(HERE / 'gr' / 'qm' / 'X.json')
    d3 = read_json(HERE / 'gr' / 'd3.json')['points']['X']['e_d3bj']
    d = r['e_tot'] - gate['e_tot']
    ok = bool(rc == 0 and abs(d) < 1e-7 and abs(r['d3_hartree'] - d3) < 1e-10 and r['xyz_sha256'] == gate['xyz_sha256']
              and not r['test_level'])
    nc.write_json(out / ('check_X.json' if not r['test_level'] else 'check_X_test_level.json'), dict(
        what='point X of gate G-R computed by md/qm/n2_batch.py against md/qm/gr/qm/X.json',
        this_program={k: r[k] for k in ('e_tot', 'd3_hartree', 'scf_cycles', 'wall_seconds', 'threads', 'host',
                                        'guess', 'grid_points', 'settings', 'test_level')},
        gate=dict(e_tot=gate['e_tot'], host=gate['host'], d3_hartree=d3),
        difference_hartree=d, tolerance_hartree=1e-7, agree=ok, written=nc.now()))
    print(f"X: this program {r['e_tot']:.10f}, gate {gate['e_tot']:.10f}, difference {d:.2e} Hartree, agree {ok}; "
          f"{r['scf_cycles']} cycles, {r['wall_seconds']:.0f} s with {r['threads']} threads")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'point', 'status', 'complete', 'd3check', 'checkx'])
    ap.add_argument('id', nargs='?')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--deadline', default=None)
    ap.add_argument('--list', choices=['main', 'all'], default='main')
    ap.add_argument('--only', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--bulk', default=None)
    ap.add_argument('--override', default=None)
    args = ap.parse_args()
    if args.cmd == 'checkx' and args.threads:
        from pyscf import lib
        lib.num_threads(args.threads)
    sys.exit(dict(run=cmd_run, point=cmd_point, status=cmd_status, complete=cmd_complete, d3check=cmd_d3check,
                  checkx=cmd_checkx)[args.cmd](args))


if __name__ == '__main__':
    main()
