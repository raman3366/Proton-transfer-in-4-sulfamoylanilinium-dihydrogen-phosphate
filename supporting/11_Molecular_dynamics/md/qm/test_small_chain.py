"""End-to-end test of the chain on a SMALL system, at the level of theory of the real jobs.

System: [Zn(imidazole)(NH3)(OH)]+, 16 atoms, charge +1, singlet.  Two parts, two environments:

    OMP_NUM_THREADS=2 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/test_small_chain.py qm
    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/test_small_chain.py readback

Part `qm` (PySCF): writes two Gaussian inputs in the form MCPB.py writes them (TINY_small_opt.com,
TINY_large_mk.com), runs the production driver md/qm/run_qm.py on them (optimisation, force constants by
finite differences for the zinc and its three ligating atoms, potential on Merz-Kollman points, adapter
files), and then checks the numerics against independent calculations:
  C1  the finite-difference columns against PySCF's analytic density-fitted Hessian (same grid);
  C2  all harmonic frequencies of the analytic Hessian real (the optimised geometry is a minimum);
  C3  Seminario constants from the partial matrix equal those from the full analytic matrix;
  C4  the error of density fitting: energy, potential at the fit points and Seminario constants against a
      calculation with exact four-centre integrals;
  C5  the potential by three-centre integrals against the one-electron route (in the driver, 20 points).
Part `readback` (AmberTools): reads the two adapter files with MCPB.py's own code and with espgen and resp:
  R1  geometry and force constants from the fchk through pymsmt = the arrays the driver saved;
  R2  Seminario bond and angle constants by pymsmt from the fchk = own code on the arrays;
  R3  the .esp file made by pymsmt and by espgen from the log = the arrays the driver saved;
  R4  resp without restraints on that .esp file gives the charges of an own least-squares fit of the arrays
      (this ties the units of the log, Angstrom, to the units of the fit, Bohr, through a physical result).

Scratch: the folder SCRATCH below.  Result: md/qm/test_small_chain.json (both parts merge into it).

    /root/cqc-sfa-tools/env/bin/python md/qm/test_small_chain.py finalize     writes the verdicts

The numerical part is slow: PySCF's analytic Hessian of this 16-atom complex took 21 minutes here (C1-C3),
and the one with exact integrals (C4) was stopped on 6 Oct 2026 to stay within the thread limit; the error
of density fitting on energy, gradient and potential is measured by md/qm/test_df_error.py instead.
`finalize` then marks C4 as not run.  As run on 6 Oct 2026, C1 and C2 FAIL: see md/ZINC.md, section 6.1.
"""
import json
import os
import pathlib
import subprocess
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import gauout  # noqa: E402
import seminario  # noqa: E402

SCRATCH = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                       '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdprep/zinc/chain_small')
OUT = HERE / 'test_small_chain.json'
BOHR = gaucom.BOHR

# imidazole in the xy plane with N3 towards +x, zinc on the x axis, NH3 and OH above and below the plane
ATOMS = [('N', 1.100, 0.000, 0.000), ('C', 0.340, 1.046, 0.000), ('N', -0.890, 0.647, 0.000),
         ('C', -0.890, -0.647, 0.000), ('C', 0.340, -1.046, 0.000), ('H', -1.707, 1.241, 0.000),
         ('H', 0.674, 2.073, 0.000), ('H', -1.764, -1.282, 0.000), ('H', 0.674, -2.073, 0.000),
         ('Zn', 3.100, 0.000, 0.000),
         ('N', 4.125, 0.000, 1.775), ('H', 4.295, 0.961, 2.072), ('H', 3.378, -0.481, 2.274),
         ('H', 5.043, -0.481, 1.312),
         ('O', 4.025, 0.000, -1.602), ('H', 4.803, 0.486, -1.923)]
ZN, N_IM, N_AM, O_OH = 10, 1, 11, 15          # indices from 1


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save(d):
    OUT.write_text(json.dumps(d, indent=1) + '\n')


def write_coms():
    body = '\n'.join('%-6s   %8.3f %8.3f %8.3f' % a for a in ATOMS)
    body_small = body
    seed = SCRATCH / 'seed.npy'
    if seed.exists():
        # to save time: the optimisation starts from the last geometry of an earlier run of this test
        # (40 steps with Cartesian d functions, stopped when the basis was changed to pure functions;
        # md/qm/test_small_chain_cart_steps.jsonl).  No criterion of the test depends on the start.
        x = np.load(seed)
        body_small = '\n'.join('%-6s   %9.4f %9.4f %9.4f' % (a[0], *r) for a, r in zip(ATOMS, x))
    (SCRATCH / 'TINY_small_opt.com').write_text(
        '%Chk=TINY_small_opt.chk\n%Mem=3000MB\n%NProcShared=2\n'
        '# B3LYP/6-31G* Geom=PrintInputOrient Integral=(Grid=UltraFine) Opt\n \nCLR\n \n1  1\n' + body_small + '\n \n')
    (SCRATCH / 'TINY_large_mk.com').write_text(
        '%Chk=TINY_large_mk.chk\n%Mem=3000MB\n%NProcShared=2\n'
        '# B3LYP/6-31G* Integral=(Grid=UltraFine) Pop(MK,ReadRadii)\nIOp(6/33=2,6/42=6)\n \nCLR\n \n1  1\n' + body +
        '\n \nZn 1.456 !VDW radius of Zn2+ from Smith et al. JCTC 2023, 19, 2064\n \n \n')


def after_driver(res):
    """What the read-back part needs, and check C5; written as soon as the driver has finished."""
    d = np.load(SCRATCH / 'bulk' / 'hess_columns.npz')
    e = np.load(SCRATCH / 'bulk' / 'esp.npz')
    np.savez(SCRATCH / 'reference_arrays.npz', coords_bohr=d['coords_bohr'], hessian=d['hessian'],
             points_ang=e['points_ang'], esp=e['esp'], coords_large_ang=e['coords_ang'], q_own=e['charges_unrestrained'])
    c5 = res['driver']['status']['jobs']['esp']['check_rinv_max_abs_difference']
    res['C5_potential_two_integral_routes'] = dict(max_abs_difference=c5, passed=bool(c5 < 1e-8))
    save(res)


def part_finalize():
    """Verdicts from what is in the result file (used when the numerical part was stopped before its end)."""
    res = load()
    if 'C5_potential_two_integral_routes' not in res:
        c5 = res['driver']['status']['jobs']['esp']['check_rinv_max_abs_difference']
        res['C5_potential_two_integral_routes'] = dict(max_abs_difference=c5, passed=bool(c5 < 1e-8))
    checks = ['C1_fd_against_analytic_hessian', 'C2_frequencies', 'C3_seminario_partial_against_full',
              'C4_density_fitting_error', 'C5_potential_two_integral_routes']
    res['numerical_checks'] = {k: (res[k]['passed'] if k in res else 'not run') for k in checks}
    res['qm_part_passed'] = all(v is True for v in res['numerical_checks'].values())
    jobs = res['driver']['status']['jobs']
    res['chain_passed'] = bool(all(jobs[j]['state'] == 'done' for j in ('opt', 'hess', 'esp', 'write')) and
                               res.get('readback_part_passed'))
    res['all_passed'] = bool(res['chain_passed'] and res['qm_part_passed'])
    save(res)
    print(json.dumps(dict(chain_passed=res['chain_passed'], numerical_checks=res['numerical_checks'],
                          all_passed=res['all_passed']), indent=1))


def part_qm():
    import pyscf_calc as pc
    import run_qm
    SCRATCH.mkdir(parents=True, exist_ok=True)
    write_coms()
    cfg = dict(prefix='TINY', indir=str(SCRATCH), bulk=str(SCRATCH / 'bulk'), outdir=str(SCRATCH / 'out'),
               threads=2, max_memory_mb=4000,
               method=dict(xc='b3lyp', basis='6-31g*', cart=False, density_fit=True, auxbasis='def2-universal-jkfit'),
               opt=dict(grid=2, conv_tol=1e-9, maxsteps=120, frozen=[]),
               hess=dict(grid=2, conv_tol=1e-10, step_bohr=0.005, atoms=[ZN, N_IM, N_AM, O_OH]),
               esp=dict(grid_scf=2, grid_final=[99, 590], conv_tol=1e-9, density=None, layers=[1.4, 1.6, 1.8, 2.0]),
               jobs=['opt', 'hess', 'esp', 'write'])
    (SCRATCH / 'config.json').write_text(json.dumps(cfg, indent=1) + '\n')
    t0 = time.time()
    run = run_qm.Run(SCRATCH / 'config.json')
    run.run()
    res = load()
    res['system'] = '[Zn(imidazole)(NH3)(OH)]+, 16 atoms, B3LYP/6-31G* (5D, 7F), density fitting def2-universal-JKFIT'
    res['start_of_optimisation'] = 'seed.npy (see write_coms)' if (SCRATCH / 'seed.npy').exists() else 'constructed'
    res['driver'] = dict(seconds=round(time.time() - t0, 1), status=run.status)
    save(res)
    after_driver(res)

    # ---- independent numerics at the optimised geometry -------------------------------------------------------
    d = np.load(SCRATCH / 'bulk' / 'hess_columns.npz')
    x_bohr, hfd, atoms = d['coords_bohr'], d['hessian'], [int(a) for a in d['atoms']]
    com = gaucom.read_com(SCRATCH / 'TINY_small_opt.com')
    s = pc.Settings(grid=2, conv_tol=1e-10, verbose=0, max_memory_mb=4000)
    from pyscf import dft
    mol = pc.build_mol(com.symbols, x_bohr * BOHR, 1, 1, s)

    def mf_for(density_fit):
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        mf.grids.level = 2
        mf.conv_tol = 1e-10
        if density_fit:
            mf = mf.density_fit(auxbasis='def2-universal-jkfit')
        mf.kernel()
        return mf

    t = time.time()
    mf = mf_for(True)
    h4 = mf.Hessian().kernel()
    n = mol.natm
    han = h4.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    t_an = time.time() - t
    np.save(SCRATCH / 'hessian_analytic_df.npy', han)
    cols = [3 * a + k for a in atoms for k in range(3)]
    c1 = float(np.abs(hfd[:, cols] - han[:, cols]).max())
    masses = [gauout.MASS[gauout.Z_OF[s_]] for s_ in com.symbols]
    freq = seminario.frequencies(han, masses, x_bohr)
    vib = np.sort(np.delete(freq, np.argsort(np.abs(freq))[:6]))     # without the six projected motions
    lig = [N_IM - 1, N_AM - 1, O_OH - 1]
    zn = ZN - 1

    def constants(h):
        out = {}
        for a in lig:
            out[f'bond {com.symbols[a]}{a + 1}-Zn'] = round(seminario.bond(h, x_bohr, a, zn)[1], 3)
        for i in range(3):
            for l in range(i + 1, 3):
                out[f'angle {lig[i] + 1}-Zn-{lig[l] + 1}'] = round(seminario.angle(h, x_bohr, lig[i], zn, lig[l])[1], 3)
        out['angle Zn-N1-C2'] = round(seminario.angle(h, x_bohr, zn, 0, 1)[1], 3)
        out['angle Zn-O-H'] = round(seminario.angle(h, x_bohr, zn, 14, 15)[1], 3)
        return out

    k_fd, k_an = constants(hfd), constants(han)
    res['C1_fd_against_analytic_hessian'] = dict(
        max_abs_difference=c1, largest_element=float(np.abs(han).max()), units='Hartree/Bohr^2',
        analytic_seconds=round(t_an, 1), passed=bool(c1 < 2e-4),
        quality_of_fd=res['driver']['status']['jobs']['hess'])
    res['C2_frequencies'] = dict(lowest_six_abs=[round(float(x), 1) for x in np.sort(np.abs(freq))[:6]],
                                 lowest_vibrations=[round(float(x), 1) for x in vib[:4]],
                                 n_imaginary=int((vib < -5).sum()), passed=bool((vib > -5).all()))
    rel = max(abs(k_fd[k] / k_an[k] - 1) for k in k_an)
    res['C3_seminario_partial_against_full'] = dict(partial=k_fd, full=k_an, largest_relative_difference=float(rel),
                                                    passed=bool(rel < 0.01))
    save(res)
    # exact integrals
    t = time.time()
    mfx = mf_for(False)
    hx4 = mfx.Hessian().kernel()
    hex_ = hx4.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    k_ex = constants(hex_)
    e = np.load(SCRATCH / 'bulk' / 'esp.npz')
    # potential on the fit points of the large-model input at ITS geometry, with and without density fitting,
    # both on grid level 2 (the driver's value is on (99,590): third column)
    coml = gaucom.read_com(SCRATCH / 'TINY_large_mk.com')
    moll = pc.build_mol(coml.symbols, coml.coords, 1, 1, s)
    pb = e['points_ang'] / BOHR
    vv = {}
    for name, dfit in (('df', True), ('exact', False)):
        m = dft.RKS(moll)
        m.xc = 'b3lyp'
        m.grids.level = 2
        m.conv_tol = 1e-10
        if dfit:
            m = m.density_fit(auxbasis='def2-universal-jkfit')
        m.kernel()
        vv[name] = (float(m.e_tot), pc.esp_at_points(moll, m.make_rdm1(), pb))
    rel = max(abs(k_an[k] / k_ex[k] - 1) for k in k_ex)
    res['C4_density_fitting_error'] = dict(
        energy_df_minus_exact=vv['df'][0] - vv['exact'][0],
        potential_max_abs_difference=float(np.abs(vv['df'][1] - vv['exact'][1]).max()),
        potential_rms=float(np.sqrt((vv['exact'][1] ** 2).mean())),
        potential_grid_level2_against_99_590=float(np.abs(vv['df'][1] - e['esp']).max()),
        hessian_max_abs_difference=float(np.abs(han - hex_).max()),
        seminario_df=k_an, seminario_exact=k_ex, largest_relative_difference_seminario=float(rel),
        exact_seconds=round(time.time() - t, 1),
        passed=bool(rel < 0.01 and np.abs(vv['df'][1] - vv['exact'][1]).max() < 1e-4))
    np.save(SCRATCH / 'hessian_analytic_exact.npy', hex_)
    res['qm_part_passed'] = all(res[k]['passed'] for k in res if k[:2] in ('C1', 'C2', 'C3', 'C4', 'C5'))
    res['versions_qm'] = pc.versions()
    save(res)
    print(json.dumps({k: v for k, v in res.items() if k[0] == 'C'}, indent=1))
    print('QM PART', 'PASSED' if res['qm_part_passed'] else 'FAILED')


def part_readback():
    from pymsmt.mol import gauio
    from pymsmt.mol.constants import B_TO_A
    from pymsmt.mcpb.gene_final_frcmod_file import get_bond_fc_with_sem, get_ang_fc_with_sem
    amber = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
    res = load()
    ref = np.load(SCRATCH / 'reference_arrays.npz')
    n = len(ATOMS)
    fchk, log = SCRATCH / 'TINY_small_opt.fchk', SCRATCH / 'TINY_large_mk.log'
    # R1
    crds = gauio.get_crds_from_fchk(str(fchk), n)
    fcm = gauio.get_matrix_from_fchk(str(fchk), 3 * n)
    e1 = float(np.abs(np.array(crds).reshape(n, 3) - ref['coords_bohr']).max())
    e2 = float(np.abs(fcm - ref['hessian']).max())
    res['R1_fchk_through_pymsmt'] = dict(max_abs_error_coordinates_bohr=e1, max_abs_error_force_constants=e2,
                                         title=fchk.read_text().splitlines()[0], passed=bool(e1 < 1e-7 and e2 < 1e-8))
    # R2
    lig, zn = [N_IM, N_AM, O_OH], ZN
    rows = []
    for a in lig:
        dis, fc = get_bond_fc_with_sem([], crds, fcm, a, zn, 1.0, 0)
        own = seminario.bond(ref['hessian'], ref['coords_bohr'], a - 1, zn - 1)
        rows.append(dict(term=f'bond {a}-{zn}', r_pymsmt=round(float(dis), 4), k_pymsmt=round(float(np.real(fc)), 3),
                         r_own=round(own[0], 4), k_own=round(own[1], 3)))
    for i in range(3):
        for l in range(i + 1, 3):
            ang, fc = get_ang_fc_with_sem([], crds, fcm, lig[i], zn, lig[l], 1.0, 0)
            own = seminario.angle(ref['hessian'], ref['coords_bohr'], lig[i] - 1, zn - 1, lig[l] - 1)
            rows.append(dict(term=f'angle {lig[i]}-{zn}-{lig[l]}', theta_pymsmt=round(float(ang), 3),
                             k_pymsmt=round(float(np.real(fc)), 3), theta_own=round(own[0], 3), k_own=round(own[1], 3)))
    ang, fc = get_ang_fc_with_sem([], crds, fcm, zn, 1, 2, 1.0, 0)
    own = seminario.angle(ref['hessian'], ref['coords_bohr'], zn - 1, 0, 1)
    rows.append(dict(term='angle Zn-N1-C2', theta_pymsmt=round(float(ang), 3), k_pymsmt=round(float(np.real(fc)), 3),
                     theta_own=round(own[0], 3), k_own=round(own[1], 3)))
    ok = all(abs(r['k_pymsmt'] / r['k_own'] - 1) < 1e-4 for r in rows) and all(r['k_pymsmt'] > 0 for r in rows)
    res['R2_seminario_pymsmt_from_file_against_own_code'] = dict(terms=rows, passed=bool(ok),
                                                                 units='Angstrom, degrees, kcal/mol/A^2 or /rad^2')
    # R3
    esp1, esp2 = SCRATCH / 'TINY_large_mk.pymsmt.esp', SCRATCH / 'TINY_large_mk.espgen.esp'
    gauio.get_esp_from_gau(str(log), str(esp1))
    p = subprocess.run([f'{amber}/bin/espgen', '-i', str(log), '-o', str(esp2)], capture_output=True, text=True,
                       env=dict(os.environ, AMBERHOME=amber))
    a1, p1, v1 = gauout.read_esp_file(esp1)
    same = esp1.read_text() == esp2.read_text() if esp2.exists() else False
    e_at = float(np.abs(a1 - ref['coords_large_ang'] / B_TO_A).max())
    e_pt = float(np.abs(p1 - ref['points_ang'] / B_TO_A).max())
    e_v = float(np.abs(v1 - ref['esp']).max())
    res['R3_esp_file_from_log'] = dict(n_atoms=int(len(a1)), n_points=int(len(p1)), max_abs_error_atoms_bohr=e_at,
                                       max_abs_error_points_bohr=e_pt, max_abs_error_potential=e_v,
                                       espgen_file_identical_to_pymsmt_file=bool(same),
                                       passed=bool(len(p1) == len(ref['esp']) and e_at < 2e-6 and e_pt < 2e-6
                                                   and e_v < 6e-7 and same))
    # R4: resp, one stage, no restraint (qwt = 0), all charges free
    z = [gauout.Z_OF[a[0]] for a in ATOMS]
    inp = ['tiny test', ' &cntrl', ' nmol = 1,', ' ihfree = 1,', ' ioutopt = 1,', ' qwt = 0.0,', ' &end', '    1.0',
           'tiny', '%5d%5d' % (1, n)] + ['%5d%5d' % (zz, 0) for zz in z] + ['', '']
    (SCRATCH / 'resp.in').write_text('\n'.join(inp) + '\n')
    p = subprocess.run([f'{amber}/bin/resp', '-O', '-i', 'resp.in', '-o', 'resp.out', '-p', 'resp.pch', '-t',
                        'resp.chg', '-e', esp1.name, '-s', 'resp_calc.esp'], cwd=SCRATCH, capture_output=True,
                       text=True, env=dict(os.environ, AMBERHOME=amber))
    q = np.array([float(x) for x in (SCRATCH / 'resp.chg').read_text().split()]) if (SCRATCH / 'resp.chg').exists() \
        else np.array([])
    ok = len(q) == n
    r4 = dict(resp_returncode=p.returncode)
    if ok:
        dq = float(np.abs(q - ref['q_own']).max())
        r4.update(charges_resp=[round(float(x), 4) for x in q], charges_own_fit=[round(float(x), 4) for x in ref['q_own']],
                  max_abs_difference=dq, sum_resp=round(float(q.sum()), 5),
                  zinc_charge=round(float(q[ZN - 1]), 4))
        ok = dq < 2e-3 and abs(q.sum() - 1.0) < 1e-4
    else:
        r4['error'] = (p.stdout + p.stderr)[-500:]
    r4['passed'] = bool(ok)
    res['R4_resp_against_own_fit'] = r4
    res['readback_part_passed'] = all(res[k]['passed'] for k in res if k[:2] in ('R1', 'R2', 'R3', 'R4'))
    res['all_passed'] = bool(res.get('qm_part_passed') and res['readback_part_passed'])
    save(res)
    print(json.dumps({k: v for k, v in res.items() if k[0] == 'R'}, indent=1))
    print('READBACK PART', 'PASSED' if res['readback_part_passed'] else 'FAILED')


if __name__ == '__main__':
    {'qm': part_qm, 'readback': part_readback, 'finalize': part_finalize}[sys.argv[1]]()
