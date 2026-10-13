"""Gate G-R, step 5: the queue.  Waits for the running quantum-chemical job, runs the charge step and the
force-field model, then the single points of md/qm/gr_points.py one at a time, and keeps a table of
quantum-chemical against force-field energies with the gate numbers up to date.  Since 6 Oct 2026, 9 AM IST
(md/PLAN.md, entry of 8:56 AM) also: the D3(BJ) energy of every geometry, the relaxed references of part B in
quantum chemistry, the binding check and both parts of the gate in log and status after every point.

    /root/cqc-sfa-tools/env/bin/python md/qm/gr_queue.py start  CONFIG.json     start detached (setsid, nohup, nice)
    /root/cqc-sfa-tools/env/bin/python md/qm/gr_queue.py status CONFIG.json     one screen of status
    ... gr_queue.py run CONFIG.json            the controller itself (what `start` launches)
    ... gr_queue.py point CONFIG.json ID       one single point (what the controller launches)
    ... gr_queue.py complete CONFIG.json       exit code 0 when every point of the list has a converged result
    ... gr_queue.py d3-extract SRCDIR OUT.json reference data of D3 for the six elements, from the s-dftd3 source
    ... gr_queue.py d3 CONFIG.json             B3LYP-D3(BJ) of every geometry of the list -> <small_out>/d3.json
    ... gr_queue.py relax CONFIG.json START    START = X or Rx: the ligand relaxed with zinc and imidazoles fixed,
                                               B3LYP-D3(BJ)/6-31G* with the settings of the single points, geomeTRIC;
                                               the signed volume at N2 is recorded at every step (see RELAX below)
    ... gr_queue.py relax-summary CONFIG.json  <small_out>/relax_qm.json from the two relaxations as far as they are
    ... gr_queue.py evaluate CONFIG.json [SRC] no single point is computed: results of the rented machine are taken
                                               over from SRC (a copy of its md/qm/gr), then force-field model, D3,
                                               force-field relaxation, binding check and gate, as after a point

Stop cleanly:  touch <bulk>/STOP      (the running SCF stops at its next cycle and keeps its orbitals; the
                                       controller ends with state "stopped"; remove STOP before a new start)
Stop at once:  kill $(cat <bulk>/queue.pid)   (the controller ends its single point and itself; by PID)
Restart:       the same `start` command.  Finished points are skipped, a stopped point continues from its
               own orbitals, the charge step is repeated only if its input changed.

STAGES
 (i)   wait: until the process named in `wait_for.pid_file` is gone and `wait_for.status` says a state out of
       `done_states`.  If the process is gone WITHOUT such a state (failed or killed), the queue waits
       `grace_minutes` for a restart and then goes on: the single points do not need that job.  With
       `wait_for.between_points` true the queue checks again before every single point and waits while such
       a process is alive (default false: a job of stage S2 started later shares the machine with the queue).
 (ii)  charges and force field (needs no quantum chemistry of its own): when <mcpb_source>/CA2_large_mk.log
       exists, a COPY of the MCPB.py folder is made in <bulk>/mcpb_<step2> and md/qm/after_qm.sh runs there
       (MCPB.py steps 2e, 3, 4, md/qm/resp_equiv.py, tleap, checks); md/qm/gr_mm.py then builds the MM model
       of the small model.  md/qm/mcpb itself is not touched.  Without the potential the configured
       placeholder parameters are used and every table says so.
 (iii) single points in the order of points.json, B3LYP/6-31G* with the settings of stage S2 for the small
       model (md/qm/config_prod.json: `method`, `opt.grid`, `opt.conv_tol`, `scf_max_cycle`, `max_memory_mb`;
       built with md/qm/pyscf_calc.py and the grid routine of md/qm/run_qm.py).  Starting orbitals: the
       point's own checkpoint when it was stopped; else the orbital coefficients of the finished point
       whose ligand lies nearest (heavy-atom r.m.s. distance), taken over at the new geometry (PySCF's
       init_guess_by_chkfile; the basis is the same, so nothing is projected: the coefficients move with
       their atoms); for the first point the checkpoint of the timing job of stage S2 when it exists (same
       geometry); else the atomic guess.  Every result file names its guess and its number of cycles.
       Each point writes <small_out>/qm/<id>.json as soon as it is done.
 (iv)  after every point: md/qm/gr_fit.py gate -> <small_out>/gr_table.csv, gr_gate.json (and the same for
       the other MM variants), status.json.  When all points are done and the gate has failed, the fit of
       the plan runs (md/qm/gr_fit.py fit), the topology is rebuilt with the fitted terms and the gate is
       evaluated again (gr_fit_table.csv, gr_fit_gate.json).

AFTER EVERY POINT the log line and `gate_so_far` of status.json carry: part A (terms at N2; r.m.s. with and
without D3, number of comparisons, largest deviation, verdict when complete), part B (contacts; X against Rx
with and without D3, the turns, the relaxed references), and the verdict of the binding check
(md/qm/gr_mm.py binding: pending until the optimum of the small model exists).

RELAX.  B3LYP/6-31G* with the settings of the single points plus the D3(BJ) energy and gradient of class D3
(NumPy; the gradient by central differences), optimised with geomeTRIC, atoms 1-46 frozen in Cartesian
coordinates.  Starting orbitals: atomic guess at the first step, the orbitals of the previous step after.
Every step is written to <bulk>/relax/<START>_steps.jsonl with the signed volume of Zn, H and S about N2
(+0.72 in X, -0.68 in Rx), the three angles at N2 and their sum, Zn-N2-S-C4, the largest shift of a frozen
atom, and the distance of the ligand (heavy atoms, r.m.s., no fitting: the frame is common) and the energy
from where the relaxation from the OTHER start has got to.
RULE FOR A CROSSING (md/PLAN.md, entry of 6 Oct 2026, 1:09 PM IST, point 1; until then the run was stopped at
the first change of sign, state "crossed"): a relaxation is NEVER stopped at a change of sign of the signed
volume.  It runs to convergence, and it counts as crossed (`crossed`: true) only when it HAS CONVERGED with a
signed volume of the sign opposite to that of its start.  The step of the first change of sign is recorded
(`sign_first_changed_at_step`).  A record with the old state "crossed" is not an end: the run continues from
its last geometry.  The key relax.on_crossing of older configurations is read and no longer acted on.
A run that was interrupted continues from its last geometry; at every (re)start the frozen atoms are put
back exactly where the start geometry has them (geomeTRIC holds `$freeze` only approximately unless
`enforce` is set: relax.convergence = {"enforce": 0.1} makes it exact once the constraints are met within
0.1 bohr; tested 6 Oct 2026, md/GR.md section 6.6).  Results: <small_out>/relax/qm_<START>.json (every
step), qm_<START>.xyz (the final geometry), and relax_qm.json (both starts, their difference).
Config key `relax`: {maxsteps (200), convergence (geomeTRIC keywords, default its own)}.

READINGS (md/PLAN.md, 1:09 PM, point 4).  Everything that depends on the force-field model (charges, mm_*.json,
gr_*gate.json, gr_*table.csv, fit/, relax_mm.json, relax/mm_*.xyz, binding_check.json) belongs to one READING
of the gate, named by `step2` (or by the config key `reading`).  The reading "2e" is the record of 6 Oct 2026
and writes into <small_out> itself, as before; every other reading writes into <small_out>/reading_<name>/
and its own scratch folders, and a reading refuses to write into a folder that holds another one.  Shared
by all readings, because they do not depend on the force field: points, geometries, qm/<id>.json, d3.json and
the quantum-chemical relaxations (relax/qm_*, relax_qm.json).

FAILURES.  A single point that ends without a result is tried once more with less memory (the fitting
integrals then go to disk); if it fails again it is marked failed and the queue goes on; three failed points
in a row stop the queue.  A failing charge or MM step is recorded and the single points run all the same.

CONFIG (JSON): label, bulk, small_out, points, qm_config, threads, nice, wait_for {pid_file, status,
done_states, grace_minutes, between_points}, mcpb_source, step2, placeholder_params (or null), guess ("nearest", the default
described above, or "atomic": always the atomic guess), guess_checkpoint (or null), timing_reference_status
(or null), only (list of ids, tests), method_override {xc, basis, grid, conv_tol} (tests), poll_seconds,
order_first (list of ids computed first, in that order; the rest follows in the order of the points file).
Relative paths are relative to the repository.
"""
import hashlib
import json
import os
import pathlib
import shutil
import signal
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
QM_PY = '/root/cqc-sfa-tools/env/bin/python'
MD_PY = '/root/cqc-sfa-tools/md/bin/python'
AMBERHOME = '/root/cqc-sfa-tools/md'
HARTREE = 627.5094740631
FIRST_READING = '2e'                                                       # the reading of 6 Oct 2026 (md/GR.md)
PLANAR_BELOW = 0.3                                                         # |signed volume|: as md/qm/gr_mm.py relax


def now(t=None):
    """IST (UTC+5:30), 12-hour, whatever time zone the machine is set to (the rented machine runs on UTC)."""
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime((time.time() if t is None else t) + 19800))


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def read_json(path, default=None):
    try:
        return json.loads(pathlib.Path(path).read_text())
    except Exception:                                                      # noqa: BLE001
        return default


def write_json(path, obj):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=1) + '\n')
    os.replace(tmp, path)


class Config:
    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()
        c = json.loads(self.path.read_text())
        self.raw = c

        def p(v):
            return None if v is None else (pathlib.Path(v) if os.path.isabs(v) else REPO / v)
        self.label = c['label']
        self.bulk = p(c['bulk'])
        self.small = p(c['small_out'])
        self.points_file = p(c['points'])
        self.qm_config = p(c['qm_config'])
        self.threads = int(c.get('threads', 4))
        self.nice = int(c.get('nice', 10))
        self.wait = c.get('wait_for') or {}
        self.mcpb_source = p(c['mcpb_source'])
        self.step2 = c.get('step2', '2e')
        self.placeholder = p(c.get('placeholder_params'))
        self.guess_chk = p(c.get('guess_checkpoint'))
        self.guess = c.get('guess') or 'nearest'
        assert self.guess in ('nearest', 'atomic')
        self.timing_ref = p(c.get('timing_reference_status'))
        self.only = c.get('only')
        self.override = c.get('method_override') or {}
        self.relax = dict(dict(maxsteps=200, convergence={}), **(c.get('relax') or {}))
        # `on_crossing` (older configurations): read, never acted on (plan entry of 6 Oct 2026, 1:09 PM IST)
        assert self.relax.get('on_crossing') in (None, 'stop', 'continue')
        # one reading of the gate per force-field model; "2e" is the record of 6 Oct 2026 and keeps its place
        self.reading = str(c.get('reading') or self.step2)
        assert self.reading.replace('_', '').replace('-', '').isalnum(), f'reading {self.reading!r}: letters, digits, _ and - only'
        first = self.reading == FIRST_READING
        self.ff_out = self.small if first else self.small / f'reading_{self.reading}'
        self.mm_work = self.bulk / ('mm' if first else f'mm_{self.reading}')
        self.evaluate_status = self.bulk / ('evaluate_status.json' if first else f'evaluate_status_{self.reading}.json')
        self.binding_fchk = [p(v) for v in (c.get('binding_fchk') or [
            c['mcpb_source'] + '/CA2_small_opt.fchk', '/root/cqc-sfa-md/from_gpu/repo_md/qm/mcpb/CA2_small_opt.fchk'])]
        self.binding_opt = [p(v) for v in (c.get('binding_opt_geometry') or ['/root/cqc-sfa-md/qm/prod/opt_final.npy'])]
        self.poll = float(c.get('poll_seconds', 30))
        pts = json.loads(self.points_file.read_text())
        self.points = pts
        # `order_first`: ids that are computed first, in this order; the others follow in the order of the
        # points file.  An id that the file does not know is an error (a typing mistake must not drop a point).
        first = list(c.get('order_first') or [])
        unknown = [i for i in first if i not in pts['order']]
        assert not unknown and len(set(first)) == len(first), f'order_first: unknown or double ids {unknown}'
        order = first + [i for i in pts['order'] if i not in first]
        self.order = [i for i in order if (not self.only or i in self.only)]
        self.geom = self.points_file.parent / 'geom'
        self.stop_file = self.bulk / 'STOP'
        self.pid_file = self.bulk / 'queue.pid'
        self.status_file = self.bulk / 'status.json'
        self.qm_dir = self.small / 'qm'

    def settings(self):
        s2 = json.loads(self.qm_config.read_text())
        m = dict(s2['method'])
        st = dict(xc=m['xc'], basis=m['basis'], cart=m['cart'], density_fit=m['density_fit'],
                  auxbasis=m.get('auxbasis', 'def2-universal-jkfit'), grid=s2['opt']['grid'],
                  conv_tol=s2['opt']['conv_tol'], scf_max_cycle=int(s2.get('scf_max_cycle', 100)),
                  max_memory_mb=int(s2['max_memory_mb']), source=str(self.qm_config), source_sha256=sha(self.qm_config),
                  overridden=sorted(self.override))
        st.update(self.override)
        return st


# ------------------------------------------------------------------------------------------------ one single point
class StopRequested(Exception):
    pass


def nearest_finished(cfg, pid):
    """The finished point whose ligand heavy atoms lie nearest to those of `pid`, with its checkpoint."""
    import numpy as np
    lig = {p['id']: p['ligand_coordinates'] for p in cfg.points['points']}
    heavy = [n for n in lig[pid] if not n.startswith('H')]
    best = None
    for other in cfg.order:
        chk = cfg.bulk / 'chk' / f'{other}.chk'
        res = read_json(cfg.qm_dir / f'{other}.json')
        if other == pid or not chk.exists() or not res or not res.get('converged'):
            continue
        d = float(np.sqrt(np.mean([sum((a - b) ** 2 for a, b in zip(lig[pid][n], lig[other][n])) for n in heavy])))
        if best is None or d < best[1]:
            best = (other, d, chk)
    return best


def run_point(cfg, pid, max_memory_mb=None):
    sys.path.insert(0, str(HERE))
    import numpy as np
    import pyscf_calc as pc
    import run_qm as rq
    from pyscf import dft, lib, scf
    t_start = time.time()
    st = cfg.settings()
    if max_memory_mb:
        st['max_memory_mb'] = int(max_memory_mb)
    for d in ('chk', 'logs'):
        (cfg.bulk / d).mkdir(parents=True, exist_ok=True)
    lib.param.TMPDIR = str(cfg.bulk)
    lines = (cfg.geom / f'{pid}.xyz').read_text().splitlines()
    n = int(lines[0])
    sym = [l.split()[0] for l in lines[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    model = cfg.points['model']
    s = pc.Settings(xc=st['xc'], basis=st['basis'], cart=st['cart'], density_fit=st['density_fit'],
                    auxbasis=st['auxbasis'], conv_tol=st['conv_tol'], max_memory_mb=st['max_memory_mb'], verbose=4)
    functional = pc.assert_b3lyp_is_gaussian() if st['xc'].lower() == 'b3lyp' else f"test level {st['xc']}"
    log = cfg.bulk / 'logs' / f'{pid}.pyscf.log'
    mol = pc.build_mol(sym, xyz, model['charge'], model['multiplicity'], s, output=str(log))
    # as Run.make_mf of md/qm/run_qm.py
    mf = dft.RKS(mol)
    mf.xc = st['xc']
    rq.set_grid(mf, st['grid'] if isinstance(st['grid'], int) else tuple(st['grid']))
    if st['density_fit']:
        mf = mf.density_fit(auxbasis=st['auxbasis'])
    mf.conv_tol = st['conv_tol']
    mf.conv_tol_grad = max(1e-6, np.sqrt(st['conv_tol']) * 0.1)
    mf.max_cycle = st['scf_max_cycle']
    chk = cfg.bulk / 'chk' / f'{pid}.chk'
    guess = dict(kind='atomic (PySCF default)')
    dm0 = None
    own = chk.with_suffix('.part.chk')
    if rq.has_scf(own):
        dm0 = mf.from_chk(str(own))
        guess = dict(kind='own checkpoint of a stopped run')
    elif cfg.guess == 'nearest':
        near = nearest_finished(cfg, pid)
        if near:
            dm0 = mf.from_chk(str(near[2]))
            guess = dict(kind='orbital coefficients of a finished point, taken over', point=near[0],
                         ligand_heavy_rmsd_A=round(near[1], 4))
        elif cfg.guess_chk and rq.has_scf(cfg.guess_chk):
            dm0 = mf.from_chk(str(cfg.guess_chk))
            guess = dict(kind='checkpoint named in the configuration', file=str(cfg.guess_chk))
    mf.chkfile = str(own)
    cycles = []

    def callback(env):
        cycles.append(dict(cycle=int(env['cycle']), e=float(env['e_tot']), t=round(time.time() - t_start, 1)))
        if cfg.stop_file.exists():
            raise StopRequested()
    mf.callback = callback
    load0 = list(os.getloadavg())
    started = now()
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
    homo = float(mf.mo_energy[mf.mo_occ > 0].max())
    lumo = float(mf.mo_energy[mf.mo_occ == 0].min())
    res = dict(id=pid, label=cfg.label, e_tot=float(mf.e_tot), converged=bool(mf.converged),
               scf_cycles=len(cycles), second_order_solver=second_order,
               last_energy_change=(cycles[-1]['e'] - cycles[-2]['e']) if len(cycles) > 1 else None,
               homo_lumo_gap_hartree=round(lumo - homo, 6),
               natm=mol.natm, nelectron=int(mol.nelectron), nao=int(mol.nao),
               naux=int(mf.with_df.get_naoaux()) if st['density_fit'] else None,
               grid_points=int(mf.grids.weights.size), guess=guess, settings=st, functional=functional,
               threads=lib.num_threads(), started=started, finished=now(),
               wall_seconds=round(time.time() - t_start, 1), load_average_start=load0,
               load_average_end=list(os.getloadavg()), xyz=str(cfg.geom / f'{pid}.xyz'),
               xyz_sha256=sha(cfg.geom / f'{pid}.xyz'), cycle_times=[c['t'] for c in cycles],
               pyscf_log=str(log), versions=pc.versions(), host=os.uname().nodename)
    if mf.converged:
        os.replace(own, chk)
    write_json(cfg.qm_dir / f'{pid}.json', res)
    print(f"{pid}: E = {mf.e_tot:.10f} Hartree, converged {mf.converged}, {len(cycles)} cycles, "
          f"{res['wall_seconds']:.0f} s", flush=True)
    return 0 if mf.converged else 4


# ------------------------------------------------------------------------------------- D3(BJ) dispersion
# B3LYP-D3(BJ), the two-body term only (no three-body term), in NumPy, so that the rented machine, which has no
# dispersion library, gives the same numbers as the Pi.  Formulas and reference data are those of the
# s-dftd3 library (source files of its tag v1.2.1; `d3-extract` reads them and writes the part for H, C, N, O,
# S and Zn to <small_out>/d3_reference.json, with the checksums of the source files):
#   CN_i  = sum_j 1 / (1 + exp(-16 (R_ij / r_ij - 1))),  R_ij = 4/3 (rcov_i + rcov_j)
#   w_ia  = exp(-4 (CN_i - CN_a)^2) / sum_a' ...          (the references a of the element of i)
#   C6_ij = sum_ab w_ia w_jb C6ref_ab
#   E     = - sum_{i<j} C6_ij [ s6 / (r^6 + R0^6) + s8 Q_ij / (r^8 + R0^8) ],  Q_ij = 3 q_i q_j,
#           R0 = a1 sqrt(Q_ij) + a2,  q = sqrt(0.5 sqrt(Z) <r4>/<r2>)
# with s6 = 1, s8 = 1.9889, a1 = 0.3981, a2 = 4.4211 bohr for B3LYP.  `d3` checks every value against the dftd3
# library (method 'b3lyp', which must give the same number as these four parameters) and one against
# PySCF's own dispersion interface, where the environment has them (the Pi has both, the machine neither).
D3_Z = dict(H=1, C=6, N=7, O=8, S=16, Zn=30)
D3_B3LYP_BJ = dict(s6=1.0, s8=1.9889, a1=0.3981, a2=4.4211)
D3_SOURCE = 'https://raw.githubusercontent.com/dftd3/simple-dftd3/v1.2.1/src/dftd3/'
D3_FD_STEP = 1.0e-3                                                       # bohr, central differences


def d3_extract(src, out):
    """Reference data of D3 for the six elements of the model, from the Fortran source of s-dftd3."""
    import re
    import numpy as np
    src = pathlib.Path(src)
    num = re.compile(r'([-+]?\d+\.\d*(?:[eEdD][-+]?\d+)?)_wp')

    def block(text, start, end):
        a = text.index(start)
        return text[a:text.index(end, a)]
    ref = (src / 'reference.f90').read_text()
    code = '\n'.join(l.split('!')[0] for l in ref.splitlines())
    nref = [int(v) for v in re.findall(r'\b\d+\b', block(code, 'number_of_references(max_elem) = [', ']').split('= [')[1])]
    assert len(nref) == 103, len(nref)
    cn = np.array([float(v) for v in num.findall(block(code, 'reference_cn(max_ref, max_elem) = reshape([', '[max_ref, max_elem])'))])
    assert cn.size == 7 * 103, cn.size
    cn = cn.reshape((7, 103), order='F')
    body = block(code, 'c6ab_view(1:2000) = [', 'end subroutine copy_c6')          # after the allocate statement
    c6 = np.array([float(v) for v in num.findall(body)])
    assert c6.size == 7 * 7 * 103 * 104 // 2, c6.size
    c6 = c6.reshape((7, 7, 103 * 104 // 2), order='F')

    def get_c6(i, j, zi, zj):                                              # get_c6 of reference.f90, 0-based refs
        if zi > zj:
            return c6[i, j, zj + zi * (zi - 1) // 2 - 1]
        return c6[j, i, zi + zj * (zj - 1) // 2 - 1]
    cov = '\n'.join(l.split('!')[0] for l in (src / 'data' / 'covrad.f90').read_text().splitlines())
    rcov = [float(v) for v in num.findall(block(cov, 'covalent_rad_2009(max_elem) = aatoau * [', ']'))]
    r42 = '\n'.join(l.split('!')[0] for l in (src / 'data' / 'r4r2.f90').read_text().splitlines())
    r4r2 = [float(v) for v in num.findall(block(r42, 'r4_over_r2(max_elem) = [', ']'))]
    assert len(rcov) == 118 and len(r4r2) == 118, (len(rcov), len(r4r2))
    zs = sorted(D3_Z.values())
    doc = dict(
        what='reference data of DFT-D3 for H, C, N, O, S, Zn, read from the source of the s-dftd3 library by '
             'md/qm/gr_queue.py d3-extract; used by its NumPy routine for B3LYP-D3(BJ)',
        source=D3_SOURCE, source_files={f: sha(src / f) for f in ('reference.f90', 'data/covrad.f90', 'data/r4r2.f90')},
        units='', angstrom_per_bohr=0.529177210903, b3lyp_bj=D3_B3LYP_BJ, k_cn=16.0, k_weight=4.0,
        n_ref={str(z): nref[z - 1] for z in zs},
        cn_ref={str(z): [float(v) for v in cn[:nref[z - 1], z - 1]] for z in zs},
        rcov_angstrom={str(z): rcov[z - 1] for z in zs}, r4_over_r2={str(z): r4r2[z - 1] for z in zs},
        c6_ref={f'{zi}-{zj}': [[float(get_c6(i, j, zi, zj)) for j in range(nref[zj - 1])] for i in range(nref[zi - 1])]
                for zi in zs for zj in zs})
    doc['units'] = ('C6 in Hartree bohr^6; covalent radii in Angstrom as tabulated (the routine multiplies them by '
                    '4/3 and converts with `angstrom_per_bohr`, the value of the library); <r4>/<r2> in atomic units')
    write_json(out, doc)
    print(f"{out}: references per element {doc['n_ref']}")


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


def read_xyz(path):
    import numpy as np
    lines = pathlib.Path(path).read_text().splitlines()
    n = int(lines[0])
    return [l.split()[0] for l in lines[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])


def d3_points(cfg, extra=()):
    """D3(BJ) of every geometry of the points file -> <small_out>/d3.json (no SCF; seconds)."""
    import numpy as np
    sys.path.insert(0, str(HERE))
    import gaucom
    ref_file = cfg.small / 'd3_reference.json'
    ref = read_json(ref_file)
    nfix = 46
    lib = psf = None
    try:
        from dftd3.interface import DispersionModel, RationalDampingParam, ZeroDampingParam
        import dftd3
        lib = dftd3.__version__
    except Exception:                                                      # noqa: BLE001
        pass
    doc = dict(what='B3LYP-D3(BJ) dispersion energy (two-body term, no three-body term) of every geometry of gate '
                    'G-R; added to the quantum-chemical energies in md/qm/gr_fit.py (md/PLAN.md, 6 Oct 2026, 8:56 AM IST, point 1)',
               routine='md/qm/gr_queue.py (class D3, NumPy)', reference_data=str(ref_file), reference_data_sha256=sha(ref_file),
               parameters=ref['b3lyp_bj'], angstrom_per_bohr=gaucom.BOHR, host=os.uname().nodename, written=now(),
               numpy=np.__version__, dftd3_library=lib, units='Hartree unless a key says kcal/mol',
               split='ligand: atoms 47-64 alone; frame: atoms 1-46 alone; interaction = total - ligand - frame', points={})
    worst = dict(library=0.0, explicit_parameters=0.0)
    gcheck = None
    for pid in cfg.points['order']:
        f = cfg.geom / f'{pid}.xyz'
        sym, xyz = read_xyz(f)
        z = [D3_Z[s] for s in sym]
        b = xyz / gaucom.BOHR
        tot = D3(z, ref).energy(b)
        lig = D3(z[nfix:], ref).energy(b[nfix:])
        frm = D3(z[:nfix], ref).energy(b[:nfix])
        rec = dict(xyz_sha256=sha(f), e_d3bj=tot, ligand=lig, frame=frm, interaction=tot - lig - frm)
        if lib:
            m = DispersionModel(np.array(z), b)
            e_lib = float(m.get_dispersion(RationalDampingParam(method='b3lyp'), grad=False)['energy'])
            e_par = float(m.get_dispersion(RationalDampingParam(s9=0.0, **ref['b3lyp_bj']), grad=False)['energy'])
            rec.update(library_d3bj=e_lib, numpy_minus_library=tot - e_lib,
                       library_d3bj_with_three_body=float(m.get_dispersion(RationalDampingParam(method='b3lyp', atm=True), grad=False)['energy']),
                       library_d3zero=float(m.get_dispersion(ZeroDampingParam(method='b3lyp'), grad=False)['energy']))
            worst['library'] = max(worst['library'], abs(tot - e_lib))
            worst['explicit_parameters'] = max(worst['explicit_parameters'], abs(e_par - e_lib))
            if pid in ('X', 'Rx'):
                g_lib = np.array(m.get_dispersion(RationalDampingParam(method='b3lyp'), grad=True)['gradient'])
                g_own = D3(z, ref).gradient(b)
                gcheck = max(gcheck or 0.0, float(np.abs(g_own - g_lib).max()))
        doc['points'][pid] = rec
    x = doc['points'].get('X')
    for pid, rec in doc['points'].items():
        if x:
            rec['rel_X_kcal_mol'] = round((rec['e_d3bj'] - x['e_d3bj']) * HARTREE, 4)
            rec['interaction_rel_X_kcal_mol'] = round((rec['interaction'] - x['interaction']) * HARTREE, 4)
    doc['checks'] = dict(
        library_available=bool(lib),
        largest_abs_numpy_minus_library_hartree=worst['library'] if lib else None,
        largest_abs_library_method_b3lyp_minus_explicit_parameters_hartree=worst['explicit_parameters'] if lib else None,
        largest_abs_gradient_numpy_minus_library_hartree_per_bohr_X_and_Rx=gcheck)
    try:                                                                   # PySCF's own interface, one value
        from pyscf import gto
        import pyscf.dispersion
        from pyscf.dispersion.dftd3 import DFTD3Dispersion
        sym, xyz = read_xyz(cfg.geom / 'X.xyz')
        mol = gto.M(atom=[(s, tuple(v)) for s, v in zip(sym, xyz)], basis='sto-3g', charge=cfg.points['model']['charge'],
                    unit='Angstrom', verbose=0)
        e_ps = float(DFTD3Dispersion(mol, xc='b3lyp', version='d3bj', atm=False).get_dispersion()['energy'])
        doc['checks']['pyscf_dispersion'] = dict(version=pyscf.dispersion.__version__, point='X', e_d3bj=e_ps,
                                                 numpy_minus_pyscf_dispersion_hartree=doc['points']['X']['e_d3bj'] - e_ps,
                                                 note="pyscf.dispersion.dftd3.DFTD3Dispersion(mol, xc='b3lyp', version='d3bj', atm=False); "
                                                      "PySCF's unit conversion for the coordinates")
    except Exception as e:                                                 # noqa: BLE001
        doc['checks']['pyscf_dispersion'] = dict(available=False, reason=f'{type(e).__name__}: {e}')
    c = doc['checks']
    ok = (not lib) or (c['largest_abs_numpy_minus_library_hartree'] < 1e-9 and
                       c['largest_abs_library_method_b3lyp_minus_explicit_parameters_hartree'] < 1e-12 and
                       (gcheck is None or gcheck < 1e-7))
    ps = c['pyscf_dispersion']
    if 'numpy_minus_pyscf_dispersion_hartree' in ps:
        ok = ok and abs(ps['numpy_minus_pyscf_dispersion_hartree']) < 1e-8
    c['passed'] = bool(ok)
    c['criteria'] = ('NumPy against the library within 1e-9 Hartree for every point; method b3lyp against the four '
                     'parameters within 1e-12; gradient within 1e-7 Hartree/bohr; PySCF interface within 1e-8 Hartree')
    write_json(cfg.small / 'd3.json', doc)
    print(f"D3(BJ) of {len(doc['points'])} geometries -> {cfg.small / 'd3.json'}; checks: {json.dumps({k: v for k, v in c.items() if k != 'criteria'})}")
    for pid, rec in doc['points'].items():
        print(f"  {pid:9s} {rec['e_d3bj']:.8f}  rel X {rec.get('rel_X_kcal_mol')}  interaction rel X {rec.get('interaction_rel_X_kcal_mol')}")
    return 0 if ok else 1


# ------------------------------------------------------------------- relaxed references of part B (QM)
def nitrogen_label(volume, v_start):
    """Words for the configuration at N2 from the signed volume (threshold PLANAR_BELOW, that of gr_mm.py relax)."""
    if abs(volume) < PLANAR_BELOW:
        return f'planar or nearly planar (|signed volume| below {PLANAR_BELOW})'
    same = volume * v_start > 0
    return 'pyramidal, configuration of the start' if same else 'pyramidal, the configuration opposite to the start'


def relax_point(cfg, start):
    """The ligand relaxed from geom/<start>.xyz with atoms 1-46 fixed; see RELAX in the header."""
    sys.path.insert(0, str(HERE))
    import numpy as np
    import gaucom
    import pyscf_calc as pc
    import run_qm as rq
    from pyscf import dft, lib
    from pyscf.geomopt import geometric_solver
    from pyscf.geomopt.addons import as_pyscf_method
    assert start in ('X', 'Rx'), start
    other = 'Rx' if start == 'X' else 'X'
    t_start = time.time()
    st = cfg.settings()
    work, small = cfg.bulk / 'relax', cfg.small / 'relax'
    for d in (work, small, cfg.bulk / 'logs'):
        d.mkdir(parents=True, exist_ok=True)
    lib.param.TMPDIR = str(cfg.bulk)
    sym, xyz0 = read_xyz(cfg.geom / f'{start}.xyz')
    model = cfg.points['model']
    names = model['ligand_names']
    nfix = 46
    i_zn = [k for k, a in enumerate(sym) if a == 'Zn'][0]
    i_n2, i_h, i_s, i_c4 = (nfix + names.index(n) for n in ('N2', 'HN21', 'S', 'C4'))
    heavy = [k for k in range(nfix, len(sym)) if sym[k] != 'H']
    steps_file, last, traj = work / f'{start}_steps.jsonl', work / f'{start}_last.npy', work / f'{start}_traj.npz'
    out_json = small / f'qm_{start}.json'
    prev = read_json(out_json, {}) or {}
    if prev.get('state') == 'done' and prev.get('xyz_start_sha256') == sha(cfg.geom / f'{start}.xyz'):
        print(f"relax {start}: state done already, nothing to do", flush=True)
        return 0
    steps = [json.loads(l) for l in steps_file.read_text().splitlines()] if steps_file.exists() else []
    geoms = list(np.load(traj)['xyz']) if traj.exists() else []
    xyz = np.load(last) if (last.exists() and steps) else xyz0.copy()
    restarted = len(steps)
    # the frozen atoms exactly where the start has them (a run made without `enforce` lets them drift by some 0.001 A)
    undone = float(np.abs(xyz[:nfix] - xyz0[:nfix]).max())
    xyz[:nfix] = xyz0[:nfix]
    # a record made under the rule in force until 6 Oct 2026, 1:09 PM IST: stopped at the first change of sign
    old_stop = prev.get('stopped_at_the_sign_change_by_the_old_rule')
    if prev.get('state') == 'crossed':
        old_stop = dict(step=prev.get('crossed_at_step'), time=prev.get('updated'),
                        note='stopped at the first change of sign of the signed volume at N2, unconverged, by the rule in '
                             'force until the plan entry of 6 Oct 2026, 1:09 PM IST; continued from that geometry')

    def geo(x):
        def ang(a, b, c):
            u, v = x[a] - x[b], x[c] - x[b]
            return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))
        z, h, u = (x[k] - x[i_n2] for k in (i_zn, i_h, i_s))
        z, h, u = z / np.linalg.norm(z), h / np.linalg.norm(h), u / np.linalg.norm(u)
        b1 = (x[i_s] - x[i_n2]) / np.linalg.norm(x[i_s] - x[i_n2])
        v = (x[i_zn] - x[i_n2]) - np.dot(x[i_zn] - x[i_n2], b1) * b1
        w = (x[i_c4] - x[i_s]) - np.dot(x[i_c4] - x[i_s], b1) * b1
        a1, a2, a3 = ang(i_zn, i_n2, i_s), ang(i_zn, i_n2, i_h), ang(i_s, i_n2, i_h)
        return dict(signed_volume=round(float(np.dot(z, np.cross(h, u))), 4),
                    Zn_N2=round(float(np.linalg.norm(x[i_zn] - x[i_n2])), 4), N2_S=round(float(np.linalg.norm(x[i_s] - x[i_n2])), 4),
                    Zn_N2_S=round(a1, 3), Zn_N2_H=round(a2, 3), S_N2_H=round(a3, 3), angle_sum_N2=round(a1 + a2 + a3, 3),
                    C4_S_N2=round(ang(i_c4, i_s, i_n2), 3),
                    Zn_N2_S_C4=round(float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))), 3),
                    rmsd_ligand_heavy_from_start=round(float(np.sqrt(((x[heavy] - xyz0[heavy]) ** 2).sum(1).mean())), 4),
                    fixed_atoms_max_shift_A=round(float(np.abs(x[:nfix] - xyz0[:nfix]).max()), 6))

    def other_now():
        """Where the relaxation from the other start has got to: its final geometry, else its last one."""
        f_fin, f_last, f_steps = work / f'{other}_final.npy', work / f'{other}_last.npy', work / f'{other}_steps.jsonl'
        try:
            so = [json.loads(l) for l in f_steps.read_text().splitlines()] if f_steps.exists() else []
            if not so:
                return None
            fin = f_fin.exists() and (read_json(small / f'qm_{other}.json', {}) or {}).get('state') == 'done'
            xo = np.load(f_fin if fin else f_last)
            return dict(x=xo, e_total=so[-1]['e_total'], step=so[-1]['step'], converged=bool(fin))
        except Exception:                                                  # noqa: BLE001  (a file being rewritten)
            return None

    def against_other(x, e_total):
        o = other_now()
        if o is None:
            return {}
        return {f'rmsd_ligand_heavy_from_{other}_run': round(float(np.sqrt(((x[heavy] - o['x'][heavy]) ** 2).sum(1).mean())), 4),
                f'e_total_minus_{other}_run_kcal_mol': round((e_total - o['e_total']) * HARTREE, 4),
                f'{other}_run_at': dict(step=o['step'], converged=o['converged'])}
    v_start = geo(xyz0)['signed_volume']
    s = pc.Settings(xc=st['xc'], basis=st['basis'], cart=st['cart'], density_fit=st['density_fit'],
                    auxbasis=st['auxbasis'], conv_tol=st['conv_tol'], max_memory_mb=st['max_memory_mb'], verbose=4)
    functional = pc.assert_b3lyp_is_gaussian() if st['xc'].lower() == 'b3lyp' else f"test level {st['xc']}"
    log = cfg.bulk / 'logs' / f'relax_{start}.pyscf.log'
    mol = pc.build_mol(sym, xyz, model['charge'], model['multiplicity'], s, output=str(log))
    mf = dft.RKS(mol)
    mf.xc = st['xc']
    rq.set_grid(mf, st['grid'] if isinstance(st['grid'], int) else tuple(st['grid']))
    if st['density_fit']:
        mf = mf.density_fit(auxbasis=st['auxbasis'])
    mf.conv_tol = st['conv_tol']
    mf.conv_tol_grad = max(1e-6, np.sqrt(st['conv_tol']) * 0.1)
    mf.max_cycle = st['scf_max_cycle']
    scanner = mf.nuc_grad_method().as_scanner()
    if st.get('surrogate'):
        # TEST OF THE MECHANICS ONLY (method_override.surrogate): no quantum chemistry; a harmonic pull of every
        # atom towards a target geometry stands for energy and gradient.  target "X" or "Rx": that point of the
        # list (from the other configuration the run must converge and be recorded as crossed); "self": the start,
        # shifted by seeded random numbers (the run must converge, not crossed).  Every result file carries the
        # label of the test.
        sg = st['surrogate']
        if sg['target'] == 'self':
            tgt = xyz0 + np.random.default_rng(int(sg.get('seed', 1))).normal(0.0, 0.05, xyz0.shape)
        else:
            tgt = read_xyz(cfg.geom / f"{sg['target']}.xyz")[1]
        tgt[:nfix] = xyz0[:nfix]

        class Surrogate:
            converged = True

            def __call__(self, mol_):
                d = mol_.atom_coords() - tgt / gaucom.BOHR
                return 0.5 * float(sg.get('k', 0.3)) * float((d ** 2).sum()), float(sg.get('k', 0.3)) * d
        scanner = Surrogate()
        functional = f'SURROGATE, NOT QUANTUM CHEMISTRY: harmonic pull towards {sg["target"]}'
    d3 = D3([D3_Z[a] for a in sym], cfg.small / 'd3_reference.json')
    t_last = [time.time()]
    restarts = list(prev.get('restarts') or [])
    if restarted:
        restarts.append(dict(at_steps=restarted, time=now(), frozen_atoms_put_back_by_at_most_A=round(undone, 6),
                             convergence=dict(cfg.relax.get('convergence') or {})))
    head = dict(what=f'part B of gate G-R: the ligand relaxed from {start} with zinc and imidazoles fixed, '
                     f"{st['xc']}-D3(BJ)/{st['basis']} (md/qm/gr_queue.py relax)", start=start, label=cfg.label,
                xyz_start=str(cfg.geom / f'{start}.xyz'), xyz_start_sha256=sha(cfg.geom / f'{start}.xyz'),
                frozen_atoms='1-46', settings=st, functional=functional, relax=cfg.relax,
                rule='md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 1: never stopped at a change of sign of the signed volume at N2; '
                     '`crossed` is true only for a run that has CONVERGED with a signed volume of the sign opposite to its start',
                d3=dict(routine='class D3 of md/qm/gr_queue.py', reference_data_sha256=sha(cfg.small / 'd3_reference.json'),
                        gradient=f'central differences, step {D3_FD_STEP} bohr'),
                signed_volume_start=v_start, host=os.uname().nodename, threads=lib.num_threads(), versions=pc.versions(),
                started=prev.get('started', now()), restarted_at_steps=(prev.get('restarted_at_steps', []) + [restarted]) if restarted else [],
                restarts=restarts, pyscf_log=str(log))
    if old_stop:
        head['stopped_at_the_sign_change_by_the_old_rule'] = old_stop

    def sign_change():
        for r in steps:
            if r['signed_volume'] * v_start <= 0:
                return r['step']
        return None

    def save(state, **kw):
        vols = [v_start] + [r['signed_volume'] for r in steps]
        doc = dict(head, state=state, updated=now(), n_steps=len(steps), steps=steps,
                   signed_volume_N2=dict(start=v_start, end=vols[-1], smallest=min(vols), largest=max(vols),
                                         changed_sign=bool(min(vols) * max(vols) <= 0)),
                   sign_first_changed_at_step=sign_change(), crossed=False, **kw)
        if steps:
            secs = [r['seconds'] for r in steps]
            doc['seconds_per_step_mean'] = round(sum(secs) / len(secs), 1)
            doc['wall_seconds_this_run'] = round(time.time() - t_start, 1)
            doc['largest_shift_of_a_fixed_atom_over_the_run_A'] = max(
                [r.get('fixed_atoms_max_shift_A', 0.0) for r in steps] + [float(np.abs(g[:nfix] - xyz0[:nfix]).max()) for g in geoms])
        write_json(out_json, doc)
        return doc

    def write_geom(path, x, text):
        path.write_text(f'{len(sym)}\n{text}\n' + '\n'.join(
            f'{a:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}' for a, v in zip(sym, x)) + '\n')

    def f(mol_):
        if cfg.stop_file.exists():
            raise StopRequested()
        e, g = scanner(mol_)
        if not scanner.converged:
            raise RuntimeError(f'relax {start}: the SCF did not converge at step {len(steps) + 1}')
        b = mol_.atom_coords()
        x = b * gaucom.BOHR
        ed, gd = d3.energy(b), d3.gradient(b)
        gt = np.array(g) + gd
        rec = dict(step=len(steps) + 1, time=now(), seconds=round(time.time() - t_last[0], 1), e_scf=float(e), e_d3=float(ed),
                   e_total=float(e + ed), grad_rms_ligand=float(np.sqrt((gt[nfix:] ** 2).mean())),
                   grad_max_ligand=float(np.abs(gt[nfix:]).max()), **geo(x))
        rec.update(against_other(x, rec['e_total']))
        t_last[0] = time.time()
        steps.append(rec)
        geoms.append(x.copy())
        with open(steps_file, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        np.savez(traj, xyz=np.array(geoms))
        np.save(last, x)
        write_geom(small / f'qm_{start}_last.xyz', x, f"relaxation from {start}, step {rec['step']}, E+D3 = {rec['e_total']:.10f} Hartree, NOT converged")
        save('running')
        line = (f"relax {start} step {rec['step']}: E = {rec['e_scf']:.8f} + D3 {rec['e_d3']:.8f} = {rec['e_total']:.8f}  rms grad(ligand) "
                f"{rec['grad_rms_ligand']:.2e}  volume at N2 {rec['signed_volume']:+.3f}  angles at N2 {rec['Zn_N2_S']:.1f} {rec['Zn_N2_H']:.1f} "
                f"{rec['S_N2_H']:.1f} (sum {rec['angle_sum_N2']:.1f})  Zn-N2-S-C4 {rec['Zn_N2_S_C4']:+.1f}  moved {rec['rmsd_ligand_heavy_from_start']:.3f} A")
        if f'rmsd_ligand_heavy_from_{other}_run' in rec:
            line += (f"  from the {other} run {rec[f'rmsd_ligand_heavy_from_{other}_run']:.3f} A, "
                     f"{rec[f'e_total_minus_{other}_run_kcal_mol']:+.3f} kcal/mol")
        print(line + f"  {rec['seconds']:.0f} s", flush=True)
        return e + ed, gt

    cons = work / f'{start}_constraints.txt'
    cons.write_text(f'$freeze\nxyz 1-{nfix}\n')
    fake = as_pyscf_method(mol, f)
    save('running')
    try:
        conv, mol_eq = geometric_solver.kernel(fake, constraints=str(cons), maxsteps=int(cfg.relax['maxsteps']),
                                               assert_convergence=False, **dict(cfg.relax.get('convergence') or {}))
    except StopRequested:
        save('stopped')
        print(f'relax {start}: stop requested after {len(steps)} steps', flush=True)
        return 3
    x = mol_eq.atom_coords() * gaucom.BOHR
    moved = float(np.abs(x[:nfix] - xyz0[:nfix]).max())
    if not conv:
        save('not converged', largest_move_of_a_fixed_atom_A=moved)
        print(f'relax {start}: not converged in {cfg.relax["maxsteps"]} steps; run again to continue', flush=True)
        return 4
    np.save(work / f'{start}_final.npy', x)
    write_geom(small / f'qm_{start}.xyz', x, f"quantum-chemical minimum of the ligand from {start}, zinc and imidazoles fixed; "
               f"E+D3 = {steps[-1]['e_total']:.10f} Hartree")
    (small / f'qm_{start}_last.xyz').unlink(missing_ok=True)
    fin = dict(e_scf=steps[-1]['e_scf'], e_d3=steps[-1]['e_d3'], e_total=steps[-1]['e_total'], **geo(x))
    fin.update(against_other(x, fin['e_total']))
    crossed = bool(fin['signed_volume'] * v_start < 0)
    doc = save('done', finished=now(), final=fin, largest_move_of_a_fixed_atom_A=moved, geometry=str(small / f'qm_{start}.xyz'))
    doc.update(crossed=crossed, nitrogen_at_the_end=nitrogen_label(fin['signed_volume'], v_start),
               crossed_rule='converged AND the signed volume at N2 has the sign opposite to that of the start')
    write_json(out_json, doc)
    print(f"relax {start}: converged in {len(steps)} steps, E+D3 = {steps[-1]['e_total']:.10f}; signed volume at N2 {v_start:+.3f} -> "
          f"{fin['signed_volume']:+.3f} ({doc['nitrogen_at_the_end']}); crossed: {crossed}", flush=True)
    return 0


def relax_summary(cfg):
    """<small_out>/relax_qm.json: both relaxations as far as they are, and their difference when both are done."""
    import numpy as np
    small = cfg.small / 'relax'
    doc = dict(what='part B of gate G-R, relaxed references in quantum chemistry (md/qm/gr_queue.py relax): the ligand '
                    'relaxed with zinc and imidazoles fixed, from X and from Rx; energies in Hartree, differences in kcal/mol',
               rule='md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 1: a run counts as crossed only when it has converged with a signed '
                    'volume at N2 of the sign opposite to its start; a run stopped at the change of sign (old rule) is not an end',
               written=now(), starts={})
    for start in ('X', 'Rx'):
        r = read_json(small / f'qm_{start}.json')
        if not r:
            doc['starts'][start] = dict(state='not started')
            continue
        steps = r.get('steps') or []
        state = r['state']
        if state == 'done' and (small / f'qm_{start}.xyz').exists():       # a copy of the last step, taken over while it ran
            (small / f'qm_{start}_last.xyz').unlink(missing_ok=True)
        if state == 'crossed':                                             # a file written under the old rule
            state = 'stopped at the change of sign by the old rule, NOT converged'
        v0 = r['signed_volume_N2']['start']
        fin = r.get('final')
        rec = dict(state=state, converged=bool(r['state'] == 'done'), n_steps=r['n_steps'], signed_volume_N2=r['signed_volume_N2'],
                   sign_first_changed_at_step=r.get('sign_first_changed_at_step', r.get('crossed_at_step')),
                   crossed=(bool(fin['signed_volume'] * v0 < 0) if (r['state'] == 'done' and fin) else None),
                   nitrogen_at_the_end=(nitrogen_label(fin['signed_volume'], v0) if (r['state'] == 'done' and fin) else None),
                   host=r.get('host'), seconds_per_step_mean=r.get('seconds_per_step_mean'), started=r.get('started'),
                   updated=r.get('updated'), restarts=r.get('restarts'),
                   stopped_at_the_sign_change_by_the_old_rule=(r.get('stopped_at_the_sign_change_by_the_old_rule') or (
                       dict(step=r.get('crossed_at_step'), time=r.get('updated')) if r['state'] == 'crossed' else None)),
                   largest_shift_of_a_fixed_atom_A=r.get('largest_shift_of_a_fixed_atom_over_the_run_A', r.get('largest_move_of_a_fixed_atom_A')),
                   final=fin, first=({k: steps[0][k] for k in ('e_scf', 'e_d3', 'e_total')} if steps else None),
                   last=(steps[-1] if steps else None))
        if steps:
            rec['relaxation_so_far_kcal_mol'] = round((steps[-1]['e_total'] - steps[0]['e_total']) * HARTREE, 4)
            ang = steps[-1].get('angle_sum_N2') or round(steps[-1]['Zn_N2_S'] + steps[-1]['Zn_N2_H'] + steps[-1]['S_N2_H'], 3)
            rec['angle_sum_N2_last'] = ang
        doc['starts'][start] = rec
    a, b = doc['starts']['X'], doc['starts']['Rx']
    both = bool(a.get('converged') and b.get('converged'))
    doc['both_converged'] = both
    doc['configurations_kept'] = bool(not a['crossed'] and not b['crossed']) if both else None
    # the two geometries as they are now (final, else last), compared atom by atom in the common frame
    geo = {}
    for start in ('X', 'Rx'):
        for name, conv in ((f'qm_{start}.xyz', True), (f'qm_{start}_last.xyz', False)):
            if (small / name).exists():
                sym, x = read_xyz(small / name)
                geo[start] = dict(file=name, converged=conv, sym=sym, x=x)
                break
    cmp_ = None
    if len(geo) == 2 and a.get('last') and b.get('last'):
        sym = geo['X']['sym']
        heavy = [k for k in range(46, len(sym)) if sym[k] != 'H']
        d = geo['Rx']['x'] - geo['X']['x']
        cmp_ = dict(geometries={s: geo[s]['file'] for s in geo}, converged={s: geo[s]['converged'] for s in geo},
                    rmsd_ligand_heavy_A=round(float(np.sqrt((d[heavy] ** 2).sum(1).mean())), 4),
                    largest_distance_of_a_ligand_heavy_atom_A=round(float(np.linalg.norm(d[heavy], axis=1).max()), 4),
                    Rx_minus_X_with_d3=round((b['last']['e_total'] - a['last']['e_total']) * HARTREE, 4),
                    Rx_minus_X_without_d3=round((b['last']['e_scf'] - a['last']['e_scf']) * HARTREE, 4),
                    Rx_minus_X_d3_part=round((b['last']['e_d3'] - a['last']['e_d3']) * HARTREE, 4),
                    note=('both runs converged: the difference between the two end points' if both else
                          'NOT both converged: the difference between the two runs as far as they are, not a result'))
    doc['end_points'] = cmp_
    if both:
        doc['Rx_minus_X_with_d3'] = round((b['final']['e_total'] - a['final']['e_total']) * HARTREE, 4)
        doc['Rx_minus_X_without_d3'] = round((b['final']['e_scf'] - a['final']['e_scf']) * HARTREE, 4)
        doc['Rx_minus_X_d3_part'] = round((b['final']['e_d3'] - a['final']['e_d3']) * HARTREE, 4)
        doc['note'] = ('without_d3: the B3LYP energies at the geometries relaxed WITH D3 (no relaxation without D3 was made)')
    else:
        doc['Rx_minus_X_with_d3'] = doc['Rx_minus_X_without_d3'] = None
    write_json(cfg.small / 'relax_qm.json', doc)
    print('relaxed references (QM):', json.dumps({k: (v.get('state'), v.get('n_steps'), (v.get('signed_volume_N2') or {}).get('end'),
                                                      'crossed', v.get('crossed')) for k, v in doc['starts'].items()}),
          'Rx - X with D3:', doc['Rx_minus_X_with_d3'], '| end points:', json.dumps({k: v for k, v in (cmp_ or {}).items() if k != 'note'}))
    return doc


# --------------------------------------------------------------------------------------------------- the controller
def gate_summary(g):
    """The numbers of gr_gate.json that go into status.json after every point."""
    v = g.get('verdict_part_A') or {}
    wd, wo = g.get('with_d3') or {}, g.get('without_d3') or {}
    b, b0 = wd.get('part_B') or {}, wo.get('part_B') or {}
    rel = b.get('relaxed') or {}
    bc = g.get('binding_check') or {}
    return dict(
        n_finished=g.get('n_finished'), n_points=g.get('n_points'), mm_label=g.get('mm_label'), plan=g.get('plan'),
        part_A=dict(what='terms at N2; the verdict of gate G-R', rms_with_d3=v.get('rms'), rms_without_d3=v.get('rms_without_d3'),
                    n_compared=v.get('n_compared'), n_expected=v.get('n_expected'), largest=v.get('largest'),
                    largest_at=v.get('largest_at'), complete=v.get('complete'), passed=v.get('passed'), criterion=v.get('criterion'),
                    rms_single_reference_with_d3=(wd.get('part_A') or {}).get('rms_single_reference'),
                    by_cap_scheme={k: x.get('rms') for k, x in (g.get('part_A_by_cap_scheme') or {}).items()}),
        part_B=dict(what='contacts; reported, no verdict, never fitted',
                    X_to_Rx_mm_minus_qm_with_d3=(b.get('X_to_Rx') or {}).get('mm_minus_qm'),
                    X_to_Rx_mm_minus_qm_without_d3=(b0.get('X_to_Rx') or {}).get('mm_minus_qm'),
                    X_to_Rx_mm_split=(b.get('X_to_Rx') or {}).get('mm_split'),
                    turns_rms_with_d3=b.get('turns_rms_mm_minus_qm_against_base'),
                    turns_rms_without_d3=b0.get('turns_rms_mm_minus_qm_against_base'), turns_n=b.get('turns_n'),
                    relaxed_Rx_minus_X_force_field=(rel.get('force_field') or {}).get('Rx_minus_X'),
                    relaxed_force_field_configurations_kept=(rel.get('force_field') or {}).get('configurations_kept'),
                    relaxed_Rx_minus_X_quantum_chemistry=(rel.get('quantum_chemistry') or {}).get('Rx_minus_X'),
                    relaxed_quantum_chemistry_state=(rel.get('quantum_chemistry') or {}).get('state'),
                    relaxed_quantum_chemistry_crossed=(rel.get('quantum_chemistry') or {}).get('crossed'),
                    relaxed_quantum_chemistry_end_points=(rel.get('quantum_chemistry') or {}).get('end_points'),
                    relaxed_force_field_at_the_quantum_chemical_geometries=(rel.get('force_field') or {}).get('at_the_quantum_chemical_geometries'),
                    relaxed_mismatch=rel.get('mismatch')),
        reported_only={k: (x or {}).get('mm_minus_qm') for k, x in (wd.get('reported_only') or {}).items()},
        binding_check=dict(verdict=bc.get('verdict'), comparison=bc.get('comparison'),
                           against_crystal_not_the_verdict=bc.get('against_crystal_not_the_verdict')))


def gate_line(gs):
    a, b, bc = gs.get('part_A') or {}, gs.get('part_B') or {}, gs.get('binding_check') or {}
    t = (f"part A (terms at N2): rms {a.get('rms_with_d3')} with D3 in two families ({a.get('rms_single_reference_with_d3')} against a "
         f"single reference, which does not decide), {a.get('rms_without_d3')} without D3, {a.get('n_compared')} "
         f"comparisons, largest {a.get('largest')} at {a.get('largest_at')}, complete {a.get('complete')}, passed {a.get('passed')}"
         f" | part B (contacts, no verdict): X->Rx {b.get('X_to_Rx_mm_minus_qm_with_d3')} with D3, {b.get('X_to_Rx_mm_minus_qm_without_d3')} "
         f"without; turns rms {b.get('turns_rms_with_d3')} (n {b.get('turns_n')}); relaxed Rx-X: force field "
         f"{b.get('relaxed_Rx_minus_X_force_field')} (configurations kept {b.get('relaxed_force_field_configurations_kept')}), "
         f"quantum chemistry {b.get('relaxed_Rx_minus_X_quantum_chemistry')} {b.get('relaxed_quantum_chemistry_state') or ''} "
         f"crossed {b.get('relaxed_quantum_chemistry_crossed')}"
         f" | binding check: {bc.get('verdict')}")
    if bc.get('comparison'):
        t += ' (' + ', '.join(f"{k} {c['difference']:+.3f} of {c['tolerance']}" for k, c in bc['comparison'].items()) + ')'
    return t


class Queue:
    def __init__(self, cfg):
        self.cfg = cfg
        cfg.bulk.mkdir(parents=True, exist_ok=True)
        cfg.small.mkdir(parents=True, exist_ok=True)
        self.logf = open(cfg.bulk / 'queue.log', 'a')
        self.child = None
        self.killed = False
        self.status = dict(label=cfg.label, state='starting', pid=os.getpid(), host=os.uname().nodename,
                           config=str(cfg.path), config_sha256=sha(cfg.path), started=now(), points_total=len(cfg.order),
                           order=cfg.order, points={}, stages={})
        signal.signal(signal.SIGTERM, self.on_signal)
        signal.signal(signal.SIGINT, self.on_signal)

    def on_signal(self, signum, frame):
        self.killed = True
        if self.child and self.child.poll() is None:
            self.child.terminate()

    def log(self, *a):
        line = now() + '  ' + ' '.join(str(x) for x in a)
        self.logf.write(line + '\n')
        self.logf.flush()
        print(line, flush=True)

    def save(self, **kw):
        self.status.update(kw)
        self.status['updated'] = now()
        write_json(self.cfg.status_file, self.status)

    def stop_wanted(self):
        return self.killed or self.cfg.stop_file.exists()

    def pause(self, seconds):
        t = time.time()
        while time.time() - t < seconds:
            if self.stop_wanted():
                return False
            time.sleep(min(2.0, seconds))
        return True

    def run_cmd(self, cmd, log_name, threads=1, env_extra=None, cwd=None):
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads),
                   MKL_NUM_THREADS=str(threads), OPENMM_CPU_THREADS=str(threads), AMBERHOME=AMBERHOME)
        env.update(env_extra or {})
        with open(self.cfg.bulk / log_name, 'a') as f:
            f.write(f'\n=== {now()}  {" ".join(str(c) for c in cmd)}\n')
            f.flush()
            # the controller runs at the configured niceness already; `nice -n` adds, so only the rest is added
            more = self.cfg.nice - os.nice(0)
            pre = ['nice', '-n', str(more)] if more > 0 else []
            self.child = subprocess.Popen(pre + [str(c) for c in cmd], stdout=f,
                                          stderr=subprocess.STDOUT, env=env, cwd=cwd or REPO)
            self.status['child_pid'] = self.child.pid
            rc = self.child.wait()
        self.child = None
        return rc

    # --------------------------------------------------------------------------------------------- (i) the wait
    def other_job(self):
        w = self.cfg.wait
        if not w:
            return dict(alive=False, state='(nothing to wait for)', done=True)
        pid = (pathlib.Path(w['pid_file']).read_text().strip() if pathlib.Path(w['pid_file']).exists() else '')
        st = (read_json(w['status'], {}) or {}).get('state')
        return dict(alive=bool(pid) and alive(pid), pid=pid, state=st, done=st in w.get('done_states', ['finished']))

    def wait_for_job(self, first=True):
        w = self.cfg.wait
        if not w:
            return True
        gone_since = None
        told = None
        while True:
            if self.stop_wanted():
                return False
            j = self.other_job()
            if not j['alive'] and (j['done'] or not first):
                return True
            if j['alive']:
                gone_since = None
                msg = f"waiting for PID {j['pid']} (state {j['state']})"
            else:
                gone_since = gone_since or time.time()
                left = float(w.get('grace_minutes', 10)) * 60 - (time.time() - gone_since)
                if left <= 0:
                    self.log(f"the awaited job is gone with state {j['state']!r}, not restarted within "
                             f"{w.get('grace_minutes', 10)} min: going on without it")
                    self.status['stages'].setdefault('wait', {}).update(ended=now(), awaited_state=j['state'],
                                                                        clean=False)
                    return True
                msg = f"awaited PID gone but state {j['state']!r}; grace {left / 60:.1f} min"
            if msg.split('(')[0] != told:
                self.log(msg)
                told = msg.split('(')[0]
            self.save(state='waiting', waiting_for=j)
            if not self.pause(self.cfg.poll):
                return False

    # ------------------------------------------------------------------------- (ii) charges and the force field
    def charge_step(self):
        cfg = self.cfg
        rec = dict(started=now())
        logf = cfg.mcpb_source / 'CA2_large_mk.log'
        work = cfg.bulk / f'mcpb_{cfg.step2}'
        params, real = None, False
        if logf.exists():
            key = sha(logf)
            done = read_json(work / 'gr_charge_step.json', {})
            if done.get('log_sha256') == key and done.get('ok') and (work / 'SN1.mol2').exists():
                rec.update(done, reused=True)
                params, real = work, True
            else:
                if work.exists():
                    shutil.rmtree(work)
                shutil.copytree(cfg.mcpb_source, work, ignore=shutil.ignore_patterns('*.chk', 'tmp*'))
                rep = cfg.bulk / f'after_qm_report_{cfg.step2}'
                rc = self.run_cmd(['bash', HERE / 'after_qm.sh', cfg.step2, work, rep], 'charge_step.log')
                have = all((work / f).exists() for f in ('SN1.mol2', 'HD1.mol2', 'HD2.mol2', 'HE1.mol2', 'ZN1.mol2',
                                                           'CA2_mcpbpy.frcmod', 'resp_equiv.json'))
                rec.update(after_qm_returncode=rc, log_sha256=key, ok=bool(have), reused=False,
                           all_hard_checks_passed=(read_json(rep / f'after_qm_check_{cfg.step2}.json', {}) or {}).get(
                               'all_hard_checks_passed'))
                if have:
                    params, real = work, True
                    out = cfg.ff_out / 'charges'
                    (out / 'unconstrained').mkdir(parents=True, exist_ok=True)
                    for f in ('SN1.mol2', 'HD1.mol2', 'HD2.mol2', 'HE1.mol2', 'ZN1.mol2'):
                        shutil.copyfile(work / f, out / f)
                        shutil.copyfile(work / 'charges_unconstrained' / f, out / 'unconstrained' / f)
                    shutil.copyfile(work / 'resp_equiv.json', out / 'resp_equiv.json')
                    shutil.copyfile(work / 'CA2_mcpbpy.frcmod', out / f'CA2_mcpbpy_{cfg.step2}.frcmod')
                    for f in (rep / f'after_qm_check_{cfg.step2}.json',):
                        if f.exists():
                            shutil.copyfile(f, out / f.name)
                    write_json(work / 'gr_charge_step.json', rec)
                else:
                    self.log('charge step FAILED, see', cfg.bulk / 'charge_step.log')
        else:
            rec.update(ok=False, reason=f'{logf} does not exist')
            self.log(f'no potential of the large model ({logf}): the charge step is skipped')
        if params is None and cfg.placeholder and (cfg.placeholder / 'SN1.mol2').exists():
            params = cfg.placeholder
        rec.update(params=str(params) if params else None, real_charges=real, ended=now())
        self.status['stages']['charges'] = rec
        self.save()
        return params, real

    def mm_step(self, params, real):
        cfg = self.cfg
        rec = dict(started=now(), variants={})
        if params is None:
            rec['ok'] = False
            self.status['stages']['mm'] = rec
            self.save()
            return {}
        tag = cfg.step2
        what = 'RESP charges with the ligand pairs equal' if real else 'PLACEHOLDER charges (not the RESP charges)'
        variants = [('main', params, ['--n2h', 'structure'], f'step {tag}, {what}, Zn-N2-H at its value in start X')]
        variants.append(('asis', params, ['--n2h', 'asis', '--no-min'], f'step {tag}, {what}, Zn-N2-H as MCPB.py wrote it'))
        free = pathlib.Path(params) / 'charges_unconstrained'
        if free.exists():
            alt = cfg.bulk / f'mcpb_{tag}_free'                             # carries the step
            if alt.exists():
                shutil.rmtree(alt)
            alt.mkdir(parents=True)
            for f in ('sfa_anion_gaff2.frcmod', 'sfa_anion_repair.frcmod', 'CA2_mcpbpy.frcmod'):
                shutil.copyfile(pathlib.Path(params) / f, alt / f)
            for f in free.glob('*.mol2'):
                shutil.copyfile(f, alt / f.name)
            variants.append(('free', alt, ['--n2h', 'structure', '--no-min'],
                             f'step {tag}, charges of MCPB.py WITHOUT the ligand equivalences, Zn-N2-H of start X'))
        out = {}
        for name, par, extra, label in variants:
            mmj = cfg.ff_out / f'mm_{tag}_{name}.json'
            rc = self.run_cmd([MD_PY, HERE / 'gr_mm.py', 'run', '--params', par, '--work', cfg.mm_work / name,
                               '--out', mmj, '--label', label, '--points', cfg.points_file] + extra, 'mm.log')
            rec['variants'][name] = dict(returncode=rc, file=str(mmj), label=label)
            if rc == 0 and mmj.exists():
                out[name] = mmj
            else:
                self.log(f'MM model {name} FAILED, see', cfg.bulk / 'mm.log')
        rec.update(ok=bool('main' in out), ended=now())
        self.status['stages']['mm'] = rec
        self.save()
        return out

    def side_cmd(self, what, cmd, log_name='gate.log', ok=(0,), **kw):
        """A step beside the single points; a return code outside `ok` goes into the log of the queue and the status."""
        rc = self.run_cmd(cmd, log_name, **kw)
        if rc not in ok:
            self.log(f'{what} ENDED WITH RETURN CODE {rc} (see {self.cfg.bulk / log_name}); its result file is not up to date')
            self.status.setdefault('side_failures', []).append(dict(step=what, returncode=rc, time=now(), log=str(self.cfg.bulk / log_name)))
        return rc

    def own_reading(self):
        """A reading never writes into the folder of another one (md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 4)."""
        cfg = self.cfg
        cfg.ff_out.mkdir(parents=True, exist_ok=True)
        mark = cfg.ff_out / 'reading.json'
        have = read_json(mark)
        if have is None:
            g = read_json(cfg.ff_out / 'gr_gate.json')
            if g is not None:                                              # the record of 6 Oct 2026 carries no mark
                have = dict(reading=g.get('reading', FIRST_READING), step2=g.get('step2', FIRST_READING))
        if have is not None and (have.get('reading') != cfg.reading or have.get('step2') != cfg.step2):
            raise SystemExit(f"{cfg.ff_out} holds the reading {have.get('reading')!r} (step {have.get('step2')}); the reading "
                             f"{cfg.reading!r} (step {cfg.step2}) must not overwrite it: nothing written")
        if cfg.ff_out == cfg.small and cfg.reading != FIRST_READING:
            raise SystemExit(f'the reading {cfg.reading!r} would write into {cfg.small}, the place of the reading {FIRST_READING!r}')
        if read_json(mark) is None:
            write_json(mark, dict(reading=cfg.reading, step2=cfg.step2, folder=str(cfg.ff_out), first_written=now(),
                                  what=('the reading of 6 Oct 2026 (md/GR.md): its files (charges/, mm_*.json, gr_*gate.json, '
                                        'gr_*table.csv, fit/, relax_mm.json, relax/mm_*.xyz, binding_check.json) lie in this folder beside '
                                        'the files shared by all readings' if cfg.reading == FIRST_READING else
                                        'one further reading of gate G-R: everything in this folder depends on the force-field model of '
                                        'this step') + '; shared, because independent of the force field: points, geometries, qm/, '
                                       'd3.json and the quantum-chemical relaxations',
                                  shared=str(cfg.small)))

    def side_results(self, mm):
        """What the gate reads beside the single points: the D3 energies, the binding check, the relaxed references."""
        cfg = self.cfg
        if 'main' not in mm:
            return
        d3 = read_json(cfg.small / 'd3.json', {}) or {}
        have = {k: v.get('xyz_sha256') for k, v in (d3.get('points') or {}).items()}
        want = {q['id']: q['xyz_sha256'] for q in cfg.points['points']}
        if have != want:                                                   # a geometry without its D3 energy
            self.side_cmd('D3 energies', [QM_PY, pathlib.Path(__file__).resolve(), 'd3', cfg.path])
        # binding check: 0 passed, 1 failed (a result), 2 pending (no optimum yet)
        self.side_cmd('binding check', [MD_PY, HERE / 'gr_mm.py', 'binding', '--mm', mm['main'], '--out', cfg.ff_out / 'binding_check.json',
                                        '--fchk'] + cfg.binding_fchk + ['--opt-geometry'] + cfg.binding_opt, ok=(0, 1, 2))
        # the summary of the quantum-chemical relaxations first: it also removes the copy of a last step once the run
        # has converged, so that the force field is not evaluated at a geometry that is no longer the last one
        if any((cfg.small / 'relax' / f'qm_{t}.json').exists() for t in ('X', 'Rx')):
            self.side_cmd('summary of the quantum-chemical relaxations', [QM_PY, pathlib.Path(__file__).resolve(), 'relax-summary', cfg.path])
        rm = cfg.ff_out / 'relax_mm.json'
        # quantum-chemical relaxed geometries: the final one of each start, and the last one while it is not converged
        qgeo = {f'qm_{t}{tag}': cfg.small / 'relax' / f'qm_{t}{tag}.xyz' for t in ('X', 'Rx') for tag in ('', '_last')}
        old = read_json(rm, {}) or {}
        seen = {k: (v or {}).get('sha256') for k, v in (old.get('at_other_geometries') or {}).items()}
        changed = [k for k, g in qgeo.items() if (sha(g) if g.exists() else None) != seen.get(k)]
        # the force-field relaxation itself (minutes on a loaded Pi) is repeated only when the model it was made with
        # has other inputs than the model of now (gr_mm.py relax --reuse decides); the single points at the
        # quantum-chemical geometries are made again whenever one of those geometries has changed
        mmj = read_json(mm['main'], {}) or {}
        now_in = dict(inputs_sha256=mmj.get('inputs_sha256'), frcmods=mmj.get('frcmods'), n2h=mmj.get('n2h'),
                      cap_scheme=mmj.get('used_scheme'), points_sha256=mmj.get('points_sha256'))
        if changed or old.get('mm_inputs') != now_in or 'fixed_atom_tolerance_A' not in old:
            self.side_cmd('force-field relaxation and single points at the quantum-chemical geometries',
                          [MD_PY, HERE / 'gr_mm.py', 'relax', '--mm', mm['main'], '--out', rm, '--geom-out', cfg.ff_out / 'relax',
                           '--reuse', '--qm-geom'] + [f'{k}={g}' for k, g in qgeo.items()])

    def gate(self, mm):
        cfg = self.cfg
        self.own_reading()
        self.side_results(mm)
        for name, mmj in mm.items():
            self.side_cmd(f'gate ({name})', [MD_PY, HERE / 'gr_fit.py', 'gate', '--mm', mmj, '--qm', cfg.qm_dir, '--out', cfg.ff_out,
                                            '--side', cfg.small, '--ff-side', cfg.ff_out, '--reading', cfg.reading, '--step2', cfg.step2,
                                            '--name', 'gr' if name == 'main' else f'gr_{name}'])
        g = read_json(cfg.ff_out / 'gr_gate.json', {}) or {}
        self.status['gate_so_far'] = gate_summary(g)
        return g

    # ------------------------------------------------------------------------------------------ (iii) the points
    def finished(self, pid):
        r = read_json(self.cfg.qm_dir / f'{pid}.json')
        return r if (r and r.get('converged')) else None

    def points_loop(self, mm):
        cfg = self.cfg
        fails_in_row = 0
        for pid in cfg.order:
            r = self.finished(pid)
            if r:
                self.status['points'][pid] = dict(state='done', e_tot=r['e_tot'], wall_seconds=r['wall_seconds'],
                                                  scf_cycles=r['scf_cycles'])
        if mm:
            self.gate(mm)
            self.log('before the first point | ' + gate_line(self.status['gate_so_far']))
        self.save(state='running')
        for k, pid in enumerate(cfg.order):
            if self.finished(pid):
                continue
            if cfg.wait.get('between_points') and not self.wait_for_job(first=False):
                return 'stopped'
            if self.stop_wanted():
                return 'stopped'
            pot = cfg.mcpb_source / 'CA2_large_mk.log'
            if not self.status['stages'].get('charges', {}).get('real_charges') and pot.exists() and \
                    sha(pot) != self.status['stages'].get('charges', {}).get('log_sha256'):
                # the potential of the large model has appeared (or changed) since the last attempt:
                # charges and MM model now; a failed attempt on the same file is not repeated
                self.log('the potential of the large model is there now: charge step')
                self.params, real = self.charge_step()
                mm = self.mm_step(self.params, real) or mm
                self.mm = mm
            self.save(state='running', current=pid, current_started=now())
            self.log(f'point {pid} ({k + 1}/{len(cfg.order)}) started')
            ok = False
            for attempt, mem in ((1, None), (2, 3000)):
                cmd = [QM_PY, pathlib.Path(__file__).resolve(), 'point', cfg.path, pid] + ([str(mem)] if mem else [])
                rc = self.run_cmd(cmd, 'points.log', threads=cfg.threads)
                if rc == 0 and self.finished(pid):
                    ok = True
                    break
                if rc == 3 or self.stop_wanted():
                    self.status['points'][pid] = dict(state='stopped')
                    return 'stopped'
                self.log(f'point {pid}: attempt {attempt} ended with code {rc}')
            if not ok:
                fails_in_row += 1
                self.status['points'][pid] = dict(state='failed')
                self.log(f'point {pid} FAILED; going on')
                if fails_in_row >= 3:
                    self.log('three failed points in a row: the queue stops')
                    return 'failed'
                self.save()
                continue
            fails_in_row = 0
            r = self.finished(pid)
            self.status['points'][pid] = dict(state='done', e_tot=r['e_tot'], wall_seconds=r['wall_seconds'],
                                              scf_cycles=r['scf_cycles'], guess=r['guess'].get('kind'))
            if pid == 'X' and cfg.timing_ref:
                t = ((read_json(cfg.timing_ref, {}) or {}).get('jobs', {}).get('timing', {}).get('scf') or {}).get('e_tot')
                if t is not None:
                    self.status['x_against_timing_job_of_S2'] = dict(
                        timing_job=t, here=r['e_tot'], difference_hartree=r['e_tot'] - t,
                        same_within_1e_7=bool(abs(r['e_tot'] - t) < 1e-7))
            done = [p for p in cfg.order if self.status['points'].get(p, {}).get('state') == 'done']
            walls = [self.status['points'][p]['wall_seconds'] for p in done]
            left = len(cfg.order) - len(done)
            mean = sum(walls) / len(walls)
            self.status['timing'] = dict(points_done=len(done), points_left=left, mean_minutes_per_point=round(mean / 60, 1),
                                         estimated_hours_left=round(left * mean / 3600, 1),
                                         estimated_end=now(time.time() + left * mean))
            if mm:
                self.gate(mm)
                self.log(f"point {pid} done: E = {r['e_tot']:.9f}, {r['scf_cycles']} cycles, {r['wall_seconds'] / 60:.1f} min | "
                         + gate_line(self.status['gate_so_far']))
            else:
                self.log(f"point {pid} done: E = {r['e_tot']:.9f}, {r['scf_cycles']} cycles, {r['wall_seconds'] / 60:.1f} min")
            self.save()
        return 'all_done' if all(self.finished(p) for p in cfg.order) else 'incomplete'

    # ------------------------------------------------------------------------------------------- (iv) the fit
    def fit_step(self, mm, params):
        """Only when part A is complete and has failed (md/PLAN.md, 6 Oct 2026, 8:56 AM IST, point 2).  Since the entry of
        1:09 PM IST, point 4: fitted terms are used only if part A WITH them, from the rebuilt topology, is at most 1.0;
        a fit that lowers part A without reaching 1.0 is reported as "terms improved, gate still failed" and no parameter
        file is left."""
        cfg = self.cfg
        g = read_json(cfg.ff_out / 'gr_gate.json', {}) or {}
        v = g.get('verdict_part_A') or {}
        crit = v.get('criterion', 1.0)
        rec = dict(part_A_complete=v.get('complete'), part_A_passed=v.get('passed'), part_A_rms=v.get('rms'), reading=cfg.reading)
        if v.get('complete') and v.get('passed') is False and 'main' in mm:
            out = cfg.ff_out / 'fit'
            rc = self.side_cmd('fit', [MD_PY, HERE / 'gr_fit.py', 'fit', '--mm', mm['main'], '--qm', cfg.qm_dir, '--out', out,
                                       '--side', cfg.small, '--ff-side', cfg.ff_out, '--reading', cfg.reading, '--step2', cfg.step2], 'fit.log')
            frc = out / 'CA2_n2_fit.frcmod'
            fit = read_json(out / 'fit.json', {}) or {}
            rec.update(fit_returncode=rc, chosen=fit.get('chosen'), reason=fit.get('reason'), fit_outcome=fit.get('outcome'),
                       frcmod=str(frc) if frc.exists() else None)
            if rc == 0 and frc.exists():
                mmj = cfg.ff_out / f'mm_{cfg.step2}_fit.json'
                rc2 = self.side_cmd('force-field model with the fitted terms',
                                    [MD_PY, HERE / 'gr_mm.py', 'run', '--params', params, '--work', cfg.mm_work / 'fit',
                                     '--out', mmj, '--label', f'step {cfg.step2}, terms at N2 fitted (gate G-R)',
                                     '--points', cfg.points_file, '--n2h', 'asis', '--no-min', '--extra', frc], 'mm.log')
                rec['mm_after_fit_returncode'] = rc2
                after = None
                if rc2 == 0:
                    self.side_cmd('gate with the fitted terms',
                                  [MD_PY, HERE / 'gr_fit.py', 'gate', '--mm', mmj, '--qm', cfg.qm_dir, '--out', cfg.ff_out, '--side', cfg.small,
                                   '--ff-side', cfg.ff_out, '--reading', cfg.reading, '--step2', cfg.step2, '--name', 'gr_fit'])
                    ga = gate_summary(read_json(cfg.ff_out / 'gr_fit_gate.json', {}) or {})
                    after = (ga.get('part_A') or {}).get('rms_with_d3')
                    rec.update(gate_after_fit=ga)
                falls = bool(after is not None and v.get('rms') is not None and after < v['rms'])
                passes = bool(after is not None and after <= crit)
                rec.update(part_A_rms_with_the_rebuilt_topology=after, part_A_number_falls=falls,
                           part_A_with_the_fitted_terms_at_most_1=passes)
                if not (falls and passes):                                 # the number that counts is that of the rebuilt topology
                    os.replace(frc, frc.with_suffix('.frcmod.rejected'))
                    why = ('the force-field model with the fitted terms could not be evaluated' if after is None else
                           'with the rebuilt topology the number of part A does not fall' if not falls else
                           f'terms improved, gate still failed: part A {v.get("rms")} -> {after} with the fitted terms, above {crit}')
                    rec.update(frcmod=None, rejected=why + ': no parameter file, the terms stay as built')
            if rec.get('frcmod'):
                rec['outcome'] = (f"gate G-R part A failed with the terms as built ({v.get('rms')}); with the fitted terms it is "
                                  f"{rec.get('part_A_rms_with_the_rebuilt_topology')}, at most {crit}: fitted terms accepted (fit/CA2_n2_fit.frcmod)")
            elif (rec.get('rejected') or '').startswith('terms improved') or fit.get('outcome') == 'terms improved, gate still failed':
                rec['outcome'] = ('gate G-R part A FAILED; terms improved, gate still failed: the fit is not used, no parameter file, '
                                  'the terms at N2 stay as built (fit/fit.json)')
            else:
                rec['outcome'] = 'gate G-R part A FAILED and no fit is accepted: the terms at N2 stay as built (split in fit/fit.json)'
            self.log(rec['outcome'])
        elif v.get('complete'):
            rec['outcome'] = 'gate G-R part A passed: nothing is fitted'
            self.log(rec['outcome'])
        self.status['stages']['fit'] = rec
        self.save()

    def evaluate(self, src=None):
        """No single point: take over results of the rented machine, then everything the queue does after a point."""
        cfg = self.cfg
        cfg.status_file = cfg.evaluate_status                              # the status of the queue itself is left alone
        self.own_reading()
        taken = collect(cfg, src) if src else []
        self.log(f"evaluation: {len(taken)} result file(s) taken over from {src}" if src else 'evaluation (nothing taken over)')
        self.params, real = self.charge_step()
        self.mm = self.mm_step(self.params, real)
        for pid in cfg.order:
            r = self.finished(pid)
            if r:
                self.status['points'][pid] = dict(state='done', e_tot=r['e_tot'], wall_seconds=r['wall_seconds'],
                                                  scf_cycles=r['scf_cycles'], host=r.get('host'))
        if self.mm:
            self.gate(self.mm)
            self.log(f"evaluation (reading {cfg.reading}), {len(self.status['points'])}/{len(cfg.order)} points | "
                     + gate_line(self.status['gate_so_far']))
            if all(self.finished(q) for q in cfg.order):
                self.fit_step(self.mm, self.params)
        self.save(state='evaluated', taken_over=taken)
        return 0

    def run(self):
        cfg = self.cfg
        cfg.pid_file.write_text(f'{os.getpid()}\n')
        self.log(f'queue {cfg.label} started, pid {os.getpid()}, {len(cfg.order)} points, threads {cfg.threads}')
        self.save(state='waiting')
        try:
            if cfg.stop_file.exists():
                self.log('STOP file present at the start: nothing done')
                self.save(state='stopped')
                return 0
            self.status['stages']['wait'] = dict(started=now())
            if not self.wait_for_job(first=True):
                self.save(state='stopped')
                self.log('stopped while waiting')
                return 0
            self.status['stages']['wait'].setdefault('ended', now())
            self.status['stages']['wait'].setdefault('clean', True)
            self.save(state='charges')
            self.log('the awaited job has ended; charge step')
            self.params, real = self.charge_step()
            c = self.status['stages']['charges']
            self.log(f"charge step: ok {c.get('ok')}, reused {c.get('reused')}, real charges {real}; MM model")
            self.mm = self.mm_step(self.params, real)
            self.log(f"MM model: {sorted(self.mm)} built; single points")
            end = self.points_loop(self.mm)
            if end == 'all_done' and self.mm:
                self.fit_step(self.mm, self.params)
            state = dict(all_done='finished', stopped='stopped', failed='failed', incomplete='finished_with_failed_points')[end]
            self.save(state=state, current=None, ended=now())
            self.log(f'queue ended: {state}')
        except Exception as e:                                             # noqa: BLE001
            import traceback
            self.log('queue FAILED:', traceback.format_exc())
            self.save(state='failed', error=f'{type(e).__name__}: {e}')
            return 1
        return 0


def collect(cfg, src):
    """Results of the rented machine -> <small_out>.  `src` is a copy of the machine's md/qm/gr (the watchdog's
    /root/cqc-sfa-md/from_gpu/repo_md/qm/gr, or a folder fetched by hand).  A single point is taken when it is
    converged, belongs to the geometry of the list (checksum) and the Pi has no result for it; a result the Pi
    has is never replaced (a different energy for the same point is reported).  Files of the relaxations and of
    the check of the machine are copied when they are newer."""
    src = pathlib.Path(src)
    want = {q['id']: q['xyz_sha256'] for q in cfg.points['points']}
    taken = []
    for f in sorted((src / 'qm').glob('*.json')):
        r = read_json(f)
        if not r or not r.get('converged') or want.get(f.stem) != r.get('xyz_sha256'):
            continue
        mine = read_json(cfg.qm_dir / f.name)
        if mine and mine.get('converged'):
            if abs(mine['e_tot'] - r['e_tot']) > 1e-7:
                print(f"collect: {f.stem} differs between {mine.get('host')} and {r.get('host')}: {mine['e_tot'] - r['e_tot']:.3e} Hartree")
            continue
        cfg.qm_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f, cfg.qm_dir / f.name)
        taken.append(f'qm/{f.name}')
    for sub in ('relax', 'check_vm', 'check_vm/qm'):
        for f in sorted((src / sub).glob('qm_*')) + sorted((src / sub).glob('*.json')) if (src / sub).is_dir() else []:
            if not f.is_file():
                continue
            tgt = cfg.small / sub / f.name
            if not tgt.exists() or f.stat().st_mtime > tgt.stat().st_mtime:
                tgt.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, tgt)
                taken.append(f'{sub}/{f.name}')
    return taken


def start(cfg):
    cfg.bulk.mkdir(parents=True, exist_ok=True)
    if cfg.pid_file.exists() and alive(cfg.pid_file.read_text().strip()):
        sys.exit(f'a queue is running already, PID {cfg.pid_file.read_text().strip()}')
    if cfg.stop_file.exists():
        sys.exit(f'{cfg.stop_file} exists: remove it first')
    cfg.pid_file.unlink(missing_ok=True)
    out = open(cfg.bulk / 'queue.stdout', 'a')
    subprocess.Popen(['setsid', 'nohup', 'nice', '-n', str(cfg.nice), QM_PY, str(pathlib.Path(__file__).resolve()),
                      'run', str(cfg.path)], stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                     cwd=REPO)
    t = time.time()
    while not cfg.pid_file.exists() and time.time() - t < 60:
        time.sleep(0.5)
    if not cfg.pid_file.exists():
        sys.exit(f'the queue did not start, see {cfg.bulk}/queue.stdout')
    print(f"started, PID {cfg.pid_file.read_text().strip()}, {now()}; status: {cfg.status_file}")


def show(cfg):
    s = read_json(cfg.status_file)
    if not s:
        sys.exit('no status file')
    pid = s.get('pid')
    print(f"state {s.get('state')}  (PID {pid}, {'alive' if alive(pid) else 'not running'}), updated {s.get('updated')}")
    if s.get('state') == 'waiting':
        print('  waiting for:', s.get('waiting_for'))
    done = [p for p, v in s.get('points', {}).items() if v.get('state') == 'done']
    print(f"  points done {len(done)}/{s.get('points_total')}; current {s.get('current')}; timing {s.get('timing')}")
    gs = s.get('gate_so_far') or {}
    print('  ' + (gate_line(gs) if 'part_A' in gs else 'gate so far (first version of the gate): ' + json.dumps(gs)))
    for k, v in (s.get('stages') or {}).items():
        print(f'  stage {k}:', json.dumps({a: b for a, b in v.items() if a != 'variants'}))


def complete(cfg):
    """Exit code 0 when every point of the list AS IT IS NOW has a converged result (the success test of the
    task on the rented machine: points added to the list after the controller started make it false)."""
    missing = [p for p in cfg.order if not (read_json(cfg.qm_dir / f'{p}.json') or {}).get('converged')]
    print(f"{len(cfg.order) - len(missing)} of {len(cfg.order)} points have a converged result"
          + (f"; missing: {' '.join(missing)}" if missing else ''))
    return 1 if missing else 0


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ('start', 'run', 'point', 'status', 'complete', 'd3', 'd3-extract', 'relax',
                                                 'relax-summary', 'evaluate'):
        sys.exit(__doc__)
    if sys.argv[1] == 'd3-extract':
        return d3_extract(sys.argv[2], sys.argv[3])
    cfg = Config(sys.argv[2])
    if sys.argv[1] == 'd3':
        sys.exit(d3_points(cfg))
    if sys.argv[1] == 'relax':
        rc = relax_point(cfg, sys.argv[3])
        relax_summary(cfg)
        sys.exit(rc)
    if sys.argv[1] == 'relax-summary':
        relax_summary(cfg)
        return
    if sys.argv[1] == 'evaluate':
        sys.exit(Queue(cfg).evaluate(sys.argv[3] if len(sys.argv) > 3 else None))
    if sys.argv[1] == 'start':
        start(cfg)
    elif sys.argv[1] == 'status':
        show(cfg)
    elif sys.argv[1] == 'complete':
        sys.exit(complete(cfg))
    elif sys.argv[1] == 'point':
        sys.exit(run_point(cfg, sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None))
    else:
        sys.exit(Queue(cfg).run())


if __name__ == '__main__':
    main()
