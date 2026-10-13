"""Stage S2: the quantum-chemical jobs of the zinc site, one after the other, restartable.

    OMP_NUM_THREADS=4 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/run_qm.py CONFIG.json [--jobs esp,write]

Input are the Gaussian files MCPB.py step 1 wrote (<prefix>_small_opt.com, <prefix>_large_mk.com); output
are the files MCPB.py steps 2 and 3 read (<prefix>_small_opt.fchk, <prefix>_large_mk.log).  The calculation
is done with PySCF (md/qm/pyscf_calc.py says what stands for which Gaussian keyword).

Jobs (each is skipped when its result exists; a stopped job resumes from its checkpoint):

  timing one SCF and one analytic gradient of the small model at the geometry of the input (settings of
         the optimisation); the wall times go to status.json, the checkpoint serves job opt.
  opt    geometry optimisation of the small model with geomeTRIC.  Every step is appended to
         <bulk>/opt_steps.jsonl and the geometry saved; a restart continues from the last geometry.
  hess   Cartesian force constants at the optimised geometry by central differences of analytic gradients,
         for the columns of the atoms named in the configuration (the zinc and the atoms bonded to it are
         what Seminario's method needs for the zinc terms; the sulfur adds the angles at sulfur).  Every
         displaced gradient is saved; a restart computes only the missing ones.  The matrix written to the
         fchk file is symmetric; blocks between two atoms that were both not displaced are ZERO and the file
         says so in its title line.
  esp    single point of the large model, Merz-Kollman fit points, electrostatic potential.
  write  the adapter files and a summary (distances, angles, Seminario constants by own code).

Status: <bulk>/status.json (rewritten after every step).  Log: <bulk>/run_qm.log; PySCF output of each job
in <bulk>/<job>.pyscf.log.  Small results go to the folder `outdir` of the configuration.

Configuration (JSON), all keys:
  prefix, indir, bulk, outdir      names and folders
  threads, max_memory_mb, scf_max_cycle (default 100; then the second-order solver takes over)
  guess_basis                      optional: a small basis (sto-3g) in which the SCF is converged first; its
                                   orbitals, projected, start the real SCF (see Run.guess)
  method   {xc, basis, cart, density_fit, auxbasis, allow_method_override}
           xc and basis must be what the .com route says unless allow_method_override is true (tests only)
  opt      {grid, conv_tol, maxsteps (0: no optimisation, input geometry taken), frozen (indices from 1,
           small model), convergence {...}}
  hess     {grid, conv_tol, step_bohr, grid_response (default true: gradients with the derivative of the grid
           weights), atoms ("zinc_site", "zinc_site+<residue>@<atom>" or a list of indices from 1)}
  esp      {grid_scf, grid_final, conv_tol, density (null: IOp(6/42) of the .com), layers}
  jobs     list, in order
`grid` is a PySCF level (integer) or [radial, angular].
"""
import argparse
import hashlib
import json
import os
import pathlib
import sys
import time
import traceback

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import gauout  # noqa: E402
import seminario  # noqa: E402

BOHR = gaucom.BOHR
IST = time.strftime('%Z')


def now():
    return time.strftime('%Y-%m-%d %H:%M:%S %Z')


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


class Run:
    def __init__(self, cfg_path, jobs=None):
        self.cfg_path = pathlib.Path(cfg_path).resolve()
        self.cfg = json.loads(self.cfg_path.read_text())
        c = self.cfg
        self.repo = HERE.parents[1]
        self.prefix = c['prefix']
        self.indir = (self.repo / c['indir']).resolve() if not os.path.isabs(c['indir']) else pathlib.Path(c['indir'])
        self.bulk = pathlib.Path(c['bulk'])
        self.outdir = (self.repo / c['outdir']).resolve() if not os.path.isabs(c['outdir']) else pathlib.Path(c['outdir'])
        self.bulk.mkdir(parents=True, exist_ok=True)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.jobs = jobs or c['jobs']
        self.status_path = self.bulk / 'status.json'
        self.status = json.loads(self.status_path.read_text()) if self.status_path.exists() else {'jobs': {}}
        self.logf = open(self.bulk / 'run_qm.log', 'a')
        self.small = gaucom.read_com(self.indir / f'{self.prefix}_small_opt.com')
        self.large = gaucom.read_com(self.indir / f'{self.prefix}_large_mk.com')

    # ------------------------------------------------------------------------------------------- bookkeeping
    def log(self, *a):
        line = now() + '  ' + ' '.join(str(x) for x in a)
        print(line, flush=True)
        self.logf.write(line + '\n')
        self.logf.flush()

    def save_status(self, **top):
        self.status.update(top)
        self.status['updated'] = now()
        tmp = self.status_path.with_suffix('.tmp')
        tmp.write_text(json.dumps(self.status, indent=1) + '\n')
        os.replace(tmp, self.status_path)

    def job(self, name):
        return self.status['jobs'].setdefault(name, {'state': 'pending'})

    # --------------------------------------------------------------------------------------------- PySCF side
    def settings(self, grid, conv_tol):
        import pyscf_calc as pc
        m = self.cfg['method']
        s = pc.Settings(xc=m['xc'], basis=m['basis'], cart=m['cart'], density_fit=m['density_fit'],
                        auxbasis=m.get('auxbasis', 'def2-universal-jkfit'), conv_tol=conv_tol,
                        max_memory_mb=self.cfg['max_memory_mb'], verbose=4)
        s.grid = grid
        return s

    def check_method(self):
        import pyscf_calc as pc
        m = self.cfg['method']
        for com in (self.small, self.large):
            same = (com.method.lower() == m['xc'].lower() and
                    com.basis.lower().replace('(d)', '*') == m['basis'].lower().replace('(d)', '*'))
            if not same and not m.get('allow_method_override'):
                raise SystemExit(f'{com.path} asks for {com.method}/{com.basis}, the configuration for '
                                 f"{m['xc']}/{m['basis']}; set allow_method_override only for tests")
        note = pc.assert_b3lyp_is_gaussian() if m['xc'].lower() == 'b3lyp' else f"test level {m['xc']}/{m['basis']}"
        return note

    def make_mf(self, com, coords_ang, grid, conv_tol, tag, chk_tag=None, basis=None):
        import pyscf_calc as pc
        from pyscf import dft, lib, scf
        lib.param.TMPDIR = str(self.bulk)
        s = self.settings(grid, conv_tol)
        if basis:
            s.basis = basis
        mol = pc.build_mol(com.symbols, coords_ang, com.charge, com.multiplicity, s,
                           output=str(self.bulk / f'{tag}.pyscf.log'))
        if s.xc.lower() == 'hf':
            mf = scf.RHF(mol)
            if s.density_fit:
                mf = mf.density_fit(auxbasis=s.auxbasis)
        else:
            mf = dft.RKS(mol)
            mf.xc = s.xc
            set_grid(mf, grid)
            if s.density_fit:
                mf = mf.density_fit(auxbasis=s.auxbasis)
        mf.conv_tol = conv_tol
        mf.conv_tol_grad = max(1e-6, np.sqrt(conv_tol) * 0.1)
        mf.max_cycle = int(self.cfg.get('scf_max_cycle', 100))
        mf.chkfile = str(self.bulk / f'{chk_tag or tag}.chk')
        return mol, mf

    def guess(self, com, coords_ang, tag):
        """Starting density from a converged calculation in a small basis (configuration: guess_basis).

        The superposition of atomic densities is a poor start for this zinc cation (the energy rises by
        Hartrees in the first cycles); one cycle of the real job costs minutes, one of this job seconds.
        The small-basis orbitals are projected on the real basis by PySCF when the checkpoint is read."""
        gb = self.cfg.get('guess_basis')
        target = self.bulk / f'{tag}.chk'
        if not gb or has_scf(target):
            return None
        mol, mf = self.make_mf(com, coords_ang, 1, 1e-6, f'{tag}_guess', basis=gb)
        r = self.scf(mf, f'{tag} guess in {gb}')
        mol.stdout.flush()
        import shutil
        shutil.copyfile(mf.chkfile, target)
        return dict(basis=gb, **r)

    def scf(self, mf, tag, dm0=None):
        chk = pathlib.Path(mf.chkfile)
        if dm0 is None and has_scf(chk):
            try:
                dm0 = mf.from_chk(str(chk))
                self.log(f'{tag}: SCF starts from the orbitals in {chk.name}')
            except Exception as e:                                         # noqa: BLE001
                self.log(f'{tag}: checkpoint {chk.name} not usable ({e}); guess from atoms')
                dm0 = None
        t = time.time()
        mf.kernel(dm0=dm0)
        second_order = False
        if not mf.converged:
            self.log(f'{tag}: SCF not converged in {mf.max_cycle} cycles; second-order solver')
            mf2 = mf.newton()
            mf2.kernel(dm0=mf.make_rdm1())
            if not mf2.converged:
                raise RuntimeError(f'{tag}: SCF not converged')
            mf.mo_coeff, mf.mo_energy, mf.mo_occ, mf.e_tot, mf.converged = (
                mf2.mo_coeff, mf2.mo_energy, mf2.mo_occ, mf2.e_tot, True)
            second_order = True
        dt = time.time() - t
        self.log(f'{tag}: SCF E = {mf.e_tot:.10f} Hartree, {dt:.0f} s')
        return dict(e_tot=float(mf.e_tot), seconds=round(dt, 1), second_order_solver=second_order)

    # ---------------------------------------------------------------------------------------------- timing
    def job_timing(self):
        """One SCF and one analytic gradient of the small model at the geometry of the input, with the
        settings of the optimisation.  The converged checkpoint is the one job `opt` starts from."""
        j = self.job('timing')
        if j.get('state') == 'done':
            self.log('timing: done before, skipped')
            return
        o = self.cfg['opt']
        j.update(state='running', started=now(), grid=o['grid'], load_average_start=list(os.getloadavg()))
        self.save_status()
        g = self.guess(self.small, np.array(self.small.coords), 'opt')
        if g:
            j.update(guess=g)
        mol, mf = self.make_mf(self.small, np.array(self.small.coords), o['grid'], o['conv_tol'], 'timing',
                               chk_tag='opt')
        r = self.scf(mf, 'timing')
        mol.stdout.flush()
        ncyc = sum(1 for l in open(self.bulk / 'timing.pyscf.log') if l.startswith('cycle='))
        j.update(natm=mol.natm, nao=int(mol.nao), scf=r, scf_cycles_in_log=ncyc,
                 naux=int(mf.with_df.get_naoaux()) if self.cfg['method']['density_fit'] else None,
                 grid_points=int(mf.grids.weights.size) if hasattr(mf, 'grids') else None)
        self.save_status()
        t = time.time()
        g = mf.nuc_grad_method().kernel()
        np.save(self.bulk / 'timing_gradient_input_geometry.npy', g)
        j.update(state='done', finished=now(), gradient_seconds=round(time.time() - t, 1),
                 gradient_rms=float(np.sqrt((g ** 2).mean())), gradient_max=float(np.abs(g).max()),
                 gradient_units='Hartree/Bohr', load_average_end=list(os.getloadavg()),
                 note='gradient of the small model at the geometry MCPB.py wrote (crystal heavy atoms)')
        self.save_status()
        self.log(f"timing: SCF {r['seconds']:.0f} s, gradient {j['gradient_seconds']:.0f} s, rms gradient "
                 f"{j['gradient_rms']:.2e} Hartree/Bohr")

    # ------------------------------------------------------------------------------------------------- opt
    def job_opt(self):
        j = self.job('opt')
        final = self.bulk / 'opt_final.npy'
        if j.get('state') == 'done' and final.exists():
            self.log('opt: done before, skipped')
            return
        o = self.cfg['opt']
        steps_file = self.bulk / 'opt_steps.jsonl'
        last = self.bulk / 'opt_last.npy'
        start = np.array(self.small.coords)
        nprev = 0
        if last.exists():
            start = np.load(last)
            nprev = sum(1 for _ in open(steps_file)) if steps_file.exists() else 0
            self.log(f'opt: restart from the geometry of step {nprev}')
        j.update(state='running', started=j.get('started', now()), restarted_at_step=nprev, grid=o['grid'],
                 frozen=o.get('frozen', []))
        self.save_status()
        if o['maxsteps'] == 0:
            # no optimisation: the force constants are taken at the geometry of the input (tests; or a
            # deliberate decision, which the summary then states)
            np.save(final, start)
            j.update(state='done', finished=now(), optimised=False,
                     note='maxsteps 0: geometry of the input taken as it is, NOT a stationary point')
            self.save_status()
            self.log('opt: maxsteps 0, geometry of the input taken as it is')
            write_xyz(self.outdir / f'{self.prefix}_small_opt.xyz', self.small.symbols, start,
                      'small model, geometry of the input, NOT optimised')
            return
        g = self.guess(self.small, start, 'opt')
        if g:
            j.update(guess=g)
        mol, mf = self.make_mf(self.small, start, o['grid'], o['conv_tol'], 'opt')
        j.update(natm=mol.natm, nao=int(mol.nao))
        from pyscf.geomopt import geometric_solver
        t_last = [time.time()]

        def callback(env):
            e, g, x = float(env['energy']), np.array(env['gradients']), np.array(env['coords']).reshape(-1, 3)
            n = nprev + int(env['self'].cycle)
            np.save(last, x * BOHR)
            rec = dict(step=n, time=now(), seconds=round(time.time() - t_last[0], 1), energy=e,
                       grad_rms=float(np.sqrt((g ** 2).mean())), grad_max=float(np.abs(g).max()),
                       scf_converged=bool(env['g_scanner'].converged))
            t_last[0] = time.time()
            with open(steps_file, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
            np.save(self.bulk / 'opt_last_gradient.npy', g)
            j.update(steps=n, last=rec)
            self.save_status()
            self.log(f"opt step {n}: E = {e:.8f}  rms grad {rec['grad_rms']:.2e}  max {rec['grad_max']:.2e}  "
                     f"{rec['seconds']:.0f} s")

        kw = dict(o.get('convergence', {}))
        cons = None
        if o.get('frozen'):
            cons = str(self.bulk / 'opt_constraints.txt')
            pathlib.Path(cons).write_text('$freeze\nxyz ' + ','.join(str(i) for i in o['frozen']) + '\n')
        self.scf(mf, 'opt')
        conv, mol_eq = geometric_solver.kernel(mf, callback=callback, maxsteps=o['maxsteps'], constraints=cons,
                                               assert_convergence=True, **kw)
        x = mol_eq.atom_coords() * BOHR
        np.save(last, x)
        if not conv:
            j.update(state='not converged', finished=now())
            self.save_status()
            raise RuntimeError(f"opt: not converged in {o['maxsteps']} steps; run again to continue")
        np.save(final, x)
        steps = [json.loads(l) for l in open(steps_file)]
        j.update(state='done', finished=now(), optimised=True, final_energy=steps[-1]['energy'], n_steps=len(steps),
                 energy_went_down=bool(steps[-1]['energy'] <= steps[0]['energy']))
        self.save_status()
        write_xyz(self.outdir / f'{self.prefix}_small_opt.xyz', self.small.symbols, x,
                  f"small model, optimised, E = {steps[-1]['energy']:.10f} Hartree")

    # ------------------------------------------------------------------------------------------------ hess
    def hess_atoms(self):
        spec = self.cfg['hess']['atoms']
        if isinstance(spec, list):
            return sorted(int(i) - 1 for i in spec)
        pdb = read_pdb(self.indir / f'{self.prefix}_small.pdb')
        serial_to_idx = {a['serial']: k for k, a in enumerate(pdb)}
        idx = set()
        for l in (self.indir / f'{self.prefix}_standard.fingerprint').read_text().splitlines():
            if l.startswith('LINK'):
                for t in l.split()[1:]:
                    idx.add(serial_to_idx[int(t.split('-')[0])])
        for extra in spec.split('+')[1:]:
            resn, name = extra.split('@')
            hit = [k for k, a in enumerate(pdb) if a['resn'] == resn and a['name'] == name]
            if len(hit) != 1:
                raise SystemExit(f'hess.atoms: {extra} matches {len(hit)} atoms of the small model')
            idx.add(hit[0])
        return sorted(idx)

    def job_hess(self):
        j = self.job('hess')
        hfile = self.bulk / 'hess_columns.npz'
        if j.get('state') == 'done' and hfile.exists():
            self.log('hess: done before, skipped')
            return
        final = self.bulk / 'opt_final.npy'
        if not final.exists():
            raise RuntimeError('hess: no optimised geometry (opt_final.npy)')
        h = self.cfg['hess']
        x0 = np.load(final) / BOHR
        atoms = self.hess_atoms()
        step = float(h['step_bohr'])
        ddir = self.bulk / 'hess_disp'
        ddir.mkdir(exist_ok=True)
        todo = [(a, k, s) for a in atoms for k in range(3) for s in (1, -1)]
        j.update(state='running', started=j.get('started', now()), atoms=[a + 1 for a in atoms], step_bohr=step,
                 n_gradients=len(todo) + 1, grid=h['grid'], grid_response=bool(h.get('grid_response', True)))
        self.save_status()
        mol, mf = self.make_mf(self.small, x0 * BOHR, h['grid'], h['conv_tol'], 'hess')
        scanner = None

        def gradient(x_bohr, tag):
            nonlocal scanner
            f = ddir / f'{tag}.npy'
            if f.exists():
                return np.load(f)
            t = time.time()
            if scanner is None:
                self.scf(mf, 'hess reference')
                gobj = mf.nuc_grad_method()
                if h.get('grid_response', True) and hasattr(mf, 'grids'):
                    # derivative of the grid weights: without it the forces do not sum to zero (1e-4
                    # Hartree/Bohr on grid level 2 or 3, md/qm/test_df_error.json) and the differences of
                    # gradients carry that error into the force constants
                    gobj.grid_response = True
                scanner = gobj.as_scanner()
            m = mol.copy()
            m.set_geom_(x_bohr, unit='Bohr')
            e, g = scanner(m)
            if not scanner.converged:
                raise RuntimeError(f'hess: SCF not converged at {tag}')
            np.save(f, g)
            self.log(f'hess {tag}: E = {e:.10f}  {time.time() - t:.0f} s')
            return g

        g0 = gradient(x0, 'reference')
        done = 0
        for a, k, s in todo:
            x = x0.copy()
            x[a, k] += s * step
            gradient(x, f'a{a + 1:03d}_{"xyz"[k]}_{"p" if s > 0 else "m"}')
            done += 1
            j.update(gradients_done=done + 1)
            self.save_status()
        n = len(x0)
        cols = {}
        for a in atoms:
            for k in range(3):
                gp = np.load(ddir / f'a{a + 1:03d}_{"xyz"[k]}_p.npy')
                gm = np.load(ddir / f'a{a + 1:03d}_{"xyz"[k]}_m.npy')
                cols[3 * a + k] = ((gp - gm) / (2 * step)).ravel()
        hess, quality = assemble_partial_hessian(n, atoms, cols)
        np.savez(hfile, hessian=hess, atoms=np.array(atoms), coords_bohr=x0, gradient=g0, step=step)
        j.update(state='done', finished=now(), reference_gradient_rms=float(np.sqrt((g0 ** 2).mean())),
                 reference_gradient_max=float(np.abs(g0).max()), **quality)
        self.save_status()

    # ------------------------------------------------------------------------------------------------- esp
    def job_esp(self):
        j = self.job('esp')
        efile = self.bulk / 'esp.npz'
        if j.get('state') == 'done' and efile.exists():
            self.log('esp: done before, skipped')
            return
        import pyscf_calc as pc
        e = self.cfg['esp']
        com = self.large
        xyz = np.array(com.coords)
        density = float(e['density']) if e.get('density') else float(com.keywords['iop'].get('6/42', 1))
        pts, lay, own = gauout.mk_points(com.symbols, xyz, com.radii, density=density, layers=tuple(e['layers']))
        j.update(state='running', started=j.get('started', now()), n_points=int(len(pts)), density_per_A2=density,
                 layers=list(e['layers']), radii={**{s: gauout.MK_RADII[s] for s in sorted(set(com.symbols)) if
                                                   s in gauout.MK_RADII}, **com.radii},
                 grid_scf=e['grid_scf'], grid_final=e['grid_final'])
        self.save_status()
        g = self.guess(com, xyz, 'esp')
        if g:
            j.update(guess=g)
            self.save_status()
        mol, mf = self.make_mf(com, xyz, e['grid_scf'], max(e['conv_tol'], 1e-7) if e['grid_scf'] != e['grid_final']
                               else e['conv_tol'], 'esp')
        j.update(natm=mol.natm, nao=int(mol.nao))
        r1 = self.scf(mf, 'esp (first grid)')
        j.update(scf_first=r1)
        self.save_status()
        if e['grid_scf'] != e['grid_final'] and self.cfg['method']['xc'].lower() != 'hf':
            set_grid(mf, e['grid_final'])
            mf.conv_tol = e['conv_tol']
            mf.conv_tol_grad = max(1e-6, np.sqrt(e['conv_tol']) * 0.1)
            r2 = self.scf(mf, 'esp (final grid)', dm0=mf.make_rdm1())
            j.update(scf_final=r2)
        dm = mf.make_rdm1()
        t = time.time()
        pb = pts / BOHR
        v = pc.esp_at_points(mol, dm, pb)
        va = pc.esp_at_nuclei(mol, dm)
        # check of the fast route on 20 points spread over the list, by the one-electron integrals
        pick = np.linspace(0, len(pb) - 1, 20).astype(int)
        v_check = pc.esp_nuclear(mol, pb[pick]) + pc.esp_electronic(mol, dm, pb[pick], method='rinv')
        dip = mf.dip_moment(mol, dm, unit='AU', verbose=0)
        q, rrms = esp_fit_charges(mol.atom_coords(), pb, v, com.charge)
        np.savez(efile, points_ang=pts, layer=lay, owner=own, esp=v, esp_atoms=va, coords_ang=xyz,
                 charges_unrestrained=q)
        j.update(state='done', finished=now(), e_tot=float(mf.e_tot), esp_seconds=round(time.time() - t, 1),
                 esp_min=float(v.min()), esp_max=float(v.max()),
                 check_rinv_max_abs_difference=float(np.abs(v_check - v[pick]).max()),
                 dipole_au=[float(x) for x in dip],
                 unrestrained_fit=dict(rrms=rrms, zinc_charge=float(q[[k for k, s in enumerate(com.symbols) if s == 'Zn'][0]])
                                       if 'Zn' in com.symbols else None, sum=float(q.sum())),
                 grid_points_final=int(mf.grids.weights.size) if hasattr(mf, 'grids') else None)
        self.save_status()

    # ----------------------------------------------------------------------------------------------- write
    def job_write(self):
        j = self.job('write')
        out = {}
        m = self.cfg['method']
        level = f"{m['xc'].upper()}/{m['basis']}" + (f" DF({m.get('auxbasis')})" if m['density_fit'] else '')
        hfile, efile = self.bulk / 'hess_columns.npz', self.bulk / 'esp.npz'
        if hfile.exists():
            d = np.load(hfile)
            atoms = [int(a) for a in d['atoms']]
            z = [gauout.Z_OF[s] for s in self.small.symbols]
            opt = self.status['jobs'].get('opt', {})
            title = (f'{self.prefix} small model, PySCF {level}; PARTIAL force constants: only blocks with atom '
                     + ','.join(str(a + 1) for a in atoms))[:72]
            f = self.indir / f'{self.prefix}_small_opt.fchk'
            gauout.write_fchk(f, z, d['coords_bohr'], d['hessian'], charge=self.small.charge,
                              multiplicity=self.small.multiplicity, energy=opt.get('final_energy'),
                              gradient=d['gradient'], title=title, method='R' + m['xc'].upper(), basis=m['basis'])
            out['fchk'] = dict(file=str(f), sha256=sha256(f), partial=True, atoms_with_columns=[a + 1 for a in atoms])
            out['small_model'] = self.summary_small(d)
        if efile.exists():
            d = np.load(efile)
            f = self.indir / f'{self.prefix}_large_mk.log'
            gauout.write_mk_log(f, self.large.symbols, d['coords_ang'], d['points_ang'], d['esp_atoms'], d['esp'],
                                header_lines=[f'{self.prefix} large model, PySCF {level}',
                                              f"{len(d['points_ang'])} points on the layers {self.cfg['esp']['layers']}"])
            out['mk_log'] = dict(file=str(f), sha256=sha256(f), n_points=int(len(d['points_ang'])))
            np.savetxt(self.outdir / f'{self.prefix}_large_unrestrained_esp_charges.txt',
                       np.c_[np.arange(1, len(d['charges_unrestrained']) + 1), d['charges_unrestrained']],
                       fmt=['%4d', '%10.6f'], header='atom of the large model, charge from an unrestrained fit '
                       '(own code, total charge fixed); for a first look, NOT the RESP charges')
        j.update(state='done', finished=now(), **out)
        self.save_status()
        summary = dict(prefix=self.prefix, level=level, config=self.cfg, status=self.status)
        (self.outdir / f'{self.prefix}_qm_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
        self.log('write: adapter files written:', ', '.join(k for k in out if k in ('fchk', 'mk_log')) or 'none')

    def summary_small(self, d):
        """Zinc geometry before and after, Seminario constants by own code for the terms the columns allow."""
        x0 = np.array(self.small.coords)
        x1 = d['coords_bohr'] * BOHR
        hess = d['hessian']
        atoms = [int(a) for a in d['atoms']]
        sym = self.small.symbols
        izn = [k for k, s in enumerate(sym) if s == 'Zn']
        res = dict(bonds=[], angles=[])
        bonded = bond_list(sym, x1)
        for zn in izn:
            lig = [a for a in atoms if a != zn and sym[a] in ('N', 'O', 'S') and np.linalg.norm(x1[a] - x1[zn]) < 2.8]
            for a in lig:
                r, k = seminario.bond(hess, d['coords_bohr'], a, zn)
                res['bonds'].append(dict(atoms=[a + 1, zn + 1], elements=[sym[a], sym[zn]],
                                         r_start=round(float(np.linalg.norm(x0[a] - x0[zn])), 4), r_opt=round(r, 4),
                                         k=round(k, 2)))
            for i in range(len(lig)):
                for l in range(i + 1, len(lig)):
                    th, k = seminario.angle(hess, d['coords_bohr'], lig[i], zn, lig[l])
                    res['angles'].append(dict(atoms=[lig[i] + 1, zn + 1, lig[l] + 1],
                                              theta_start=round(angle_deg(x0[lig[i]], x0[zn], x0[lig[l]]), 2),
                                              theta_opt=round(th, 2), k=round(k, 2)))
            for a in lig:
                for b in bonded[a]:
                    if b == zn:
                        continue
                    th, k = seminario.angle(hess, d['coords_bohr'], zn, a, b)
                    res['angles'].append(dict(atoms=[zn + 1, a + 1, b + 1], elements=[sym[zn], sym[a], sym[b]],
                                              theta_start=round(angle_deg(x0[zn], x0[a], x0[b]), 2),
                                              theta_opt=round(th, 2), k=round(k, 2)))
        # other centres whose columns were computed (for example the sulfur)
        for c in atoms:
            if c in izn:
                continue
            nb = bonded[c]
            for b in nb:
                if b > c and b in atoms:
                    continue                      # the pair is listed once, at its second atom
                r, k = seminario.bond(hess, d['coords_bohr'], b, c)
                res['bonds'].append(dict(atoms=[b + 1, c + 1], elements=[sym[b], sym[c]],
                                         r_start=round(float(np.linalg.norm(x0[b] - x0[c])), 4), r_opt=round(r, 4),
                                         k=round(k, 2)))
            for i in range(len(nb)):
                for l in range(i + 1, len(nb)):
                    th, k = seminario.angle(hess, d['coords_bohr'], nb[i], c, nb[l])
                    res['angles'].append(dict(atoms=[nb[i] + 1, c + 1, nb[l] + 1],
                                              elements=[sym[nb[i]], sym[c], sym[nb[l]]],
                                              theta_start=round(angle_deg(x0[nb[i]], x0[c], x0[nb[l]]), 2),
                                              theta_opt=round(th, 2), k=round(k, 2)))
        res['units'] = 'Angstrom, degrees, kcal/mol/A^2, kcal/mol/rad^2 (AMBER form k (q - q0)^2)'
        res['rmsd_heavy_atoms_start_to_opt'] = round(float(kabsch_rmsd(
            x0[[k for k, s in enumerate(sym) if s != 'H']], x1[[k for k, s in enumerate(sym) if s != 'H']])), 3)
        return res

    # ------------------------------------------------------------------------------------------------- run
    def run(self):
        import pyscf_calc as pc
        note = self.check_method()
        (self.bulk / 'driver.pid').write_text(f'{os.getpid()}\n')
        self.save_status(pid=os.getpid(), host=os.uname().nodename, config=str(self.cfg_path),
                         config_sha256=sha256(self.cfg_path), started=self.status.get('started', now()),
                         last_start=now(), versions=pc.versions(), functional=note, jobs_requested=self.jobs,
                         state='running')
        self.log(f'driver started, pid {os.getpid()}, jobs {self.jobs}, threads {pc.versions()["threads"]}; {note}')
        try:
            for name in self.jobs:
                self.save_status(current_job=name)
                getattr(self, f'job_{name}')()
            self.save_status(state='finished', current_job=None, finished=now())
            self.log('driver finished')
        except Exception as e:                                             # noqa: BLE001
            self.status['jobs'].setdefault(self.status.get('current_job', '?'), {})['error'] = f'{type(e).__name__}: {e}'
            self.save_status(state='failed', error=f'{type(e).__name__}: {e}', failed=now())
            self.log('driver FAILED:', traceback.format_exc())
            raise


# ------------------------------------------------------------------------------------------------------ helpers
def has_scf(chk):
    """True when the checkpoint file holds orbitals."""
    chk = pathlib.Path(chk)
    if not chk.exists():
        return False
    try:
        from pyscf.lib import chkfile
        return chkfile.load(str(chk), 'scf/mo_coeff') is not None
    except Exception:                                                      # noqa: BLE001
        return False


def set_grid(mf, grid):
    if isinstance(grid, int):
        mf.grids.atom_grid = {}
        mf.grids.level = grid
    else:
        mf.grids.atom_grid = tuple(grid)
    mf.grids.coords = None
    mf.grids.weights = None


def assemble_partial_hessian(n, atoms, cols):
    """Symmetric 3N x 3N matrix from the columns of the displaced atoms; quality numbers."""
    h = np.zeros((3 * n, 3 * n))
    have = np.zeros((3 * n, 3 * n), dtype=int)
    for c, col in cols.items():
        h[:, c] += col
        have[:, c] += 1
        h[c, :] += col
        have[c, :] += 1
    idx = sorted(cols)
    # blocks between two displaced atoms were filled twice: average; record their asymmetry first
    raw = np.array([cols[c][idx] for c in idx]).T          # raw[i, j] = dg_i / dx_j for i, j displaced
    asym = float(np.abs(raw - raw.T).max()) if len(idx) > 1 else 0.0
    h[have == 2] *= 0.5
    # sum rule: a rigid translation changes no force, so every column sums to zero over the atoms
    sr = max(float(np.abs(cols[c].reshape(n, 3).sum(axis=0)).max()) for c in idx)
    return h, dict(max_asymmetry_displaced_block=asym, max_translation_sum_rule=sr,
                   largest_element=float(np.abs(h).max()), units='Hartree/Bohr^2')


def esp_fit_charges(atoms_bohr, points_bohr, v, total):
    """Unrestrained least-squares charges with the total fixed (own code; a first look and a cross-check)."""
    a = 1.0 / np.linalg.norm(points_bohr[:, None, :] - atoms_bohr[None, :, :], axis=2)
    n = a.shape[1]
    m = np.zeros((n + 1, n + 1))
    m[:n, :n] = a.T @ a
    m[n, :n] = m[:n, n] = 1.0
    rhs = np.concatenate([a.T @ v, [total]])
    q = np.linalg.solve(m, rhs)[:n]
    rrms = float(np.sqrt(((a @ q - v) ** 2).sum() / (v ** 2).sum()))
    return q, rrms


def read_pdb(path):
    out = []
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            out.append(dict(serial=int(l[6:11]), name=l[12:16].strip(), resn=l[17:20].strip(), resi=int(l[22:26])))
    return out


COVALENT = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'Zn': 1.22}


def bond_list(symbols, x):
    n = len(symbols)
    nb = {i: [] for i in range(n)}
    for i in range(n):
        for k in range(i + 1, n):
            if 'Zn' in (symbols[i], symbols[k]):
                continue
            if np.linalg.norm(x[i] - x[k]) <= COVALENT[symbols[i]] + COVALENT[symbols[k]] + 0.4:
                nb[i].append(k)
                nb[k].append(i)
    return nb


def angle_deg(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def kabsch_rmsd(p, q):
    p0, q0 = p - p.mean(0), q - q.mean(0)
    u, s, vt = np.linalg.svd(p0.T @ q0)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return np.sqrt(((p0 @ r - q0) ** 2).sum() / len(p))


def write_xyz(path, symbols, x_ang, comment):
    lines = [str(len(symbols)), comment] + ['%-2s %14.8f %14.8f %14.8f' % (s, *r) for s, r in zip(symbols, x_ang)]
    pathlib.Path(path).write_text('\n'.join(lines) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('config')
    ap.add_argument('--jobs', help='comma-separated list out of esp,write,timing,opt,hess, run in the order given')
    a = ap.parse_args()
    Run(a.config, a.jobs.split(',') if a.jobs else None).run()


if __name__ == '__main__':
    main()
