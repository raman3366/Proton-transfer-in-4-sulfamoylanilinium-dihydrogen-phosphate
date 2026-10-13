"""Round-trip tests of the adapters with MCPB.py's own parsing code and with espgen.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/test_adapters.py

No quantum chemistry here: synthetic arrays with a fixed seed are written and read back.

 T1  .com reader: the three inputs of MCPB.py step 1 against the PDB files of the same step
     (number and order of atoms, elements, coordinates, charge, multiplicity, radius of zinc).
 T2  fchk writer: geometry and Hessian read back with pymsmt (get_crds_from_fchk, get_matrix_from_fchk)
     for 2, 3, 7, 64 and 101 atoms (line counts with and without a full last line), and with an
     independent reader; the order of the lower triangle is tested with a matrix whose elements encode
     their own indices.
 T3  units through MCPB.py's Seminario code (get_bond_fc_with_sem, get_ang_fc_with_sem) reading the fchk:
     (i) a diatomic with a known bond constant: exact for the method, must come back in kcal/mol/A^2 with
     the distance in Angstrom; (ii) a triatomic with ONE harmonic angle term and no bond term, Hessian
     K g g^T written analytically: Seminario's series formula then has the closed form
     k = k_in / (1/(1 + rho c) + 1/(1 + c/rho)), c = -cos(theta), rho = R_ab/R_cb (k_in/2 at 90 degrees),
     and pymsmt must give this number in kcal/mol/rad^2 with the angle in degrees; (iii) a bent triatomic
     with two bond terms and the angle term (Hessian by finite differences of the AMBER energy
     expression): bond constants within 2 % of the input, and pymsmt on the FILE equal to
     md/qm/seminario.py (own code) on the ARRAY in memory for bonds and angle.
 T4  MK log writer: atoms, fit points and potentials read back with pymsmt (get_esp_from_gau) and with
     the program espgen, compared with the arrays written (Angstrom -> Bohr with each reader's constant).
     Includes coordinates below -10 Angstrom (fields that touch), more than 9,999 points (the index field
     of the log overflows to asterisks) and 99,999 points, the most the header of the .esp file (2I5) and
     therefore the program resp can take; the writer refuses more.
 T5  MK points: counts, distances to the atoms, layers; compared with the area of the surface.

Writes md/qm/test_adapters.json; scratch files under the scratch folder given below.
"""
import json
import os
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import gauout  # noqa: E402
import seminario  # noqa: E402

from pymsmt.mol import gauio  # noqa: E402
from pymsmt.mol.constants import B_TO_A, HB2_TO_KCAL_MOL_A2, H_TO_KCAL_MOL  # noqa: E402
from pymsmt.mcpb.gene_final_frcmod_file import get_bond_fc_with_sem, get_ang_fc_with_sem  # noqa: E402

SCRATCH = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                       '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdprep/zinc/adapters')
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
SEED = 20261005
results = []


def record(name, passed, **kw):
    results.append(dict(test=name, passed=bool(passed), **kw))
    print(('PASS ' if passed else 'FAIL ') + name, json.dumps(kw))


def pdb_atoms(path):
    out = []
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            name, resn = l[12:16].strip(), l[17:20].strip()
            el = 'Zn' if resn == 'ZN' else name[0]
            out.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


# ------------------------------------------------------------------------------------------------------- T1
def t1():
    for tag, comf, pdbf, natm, chg in (('small', 'CA2_small_opt.com', 'CA2_small.pdb', 64, 1),
                                       ('large', 'CA2_large_mk.com', 'CA2_large.pdb', 101, 1)):
        c = gaucom.read_com(HERE / 'mcpb' / comf)
        p = pdb_atoms(HERE / 'mcpb' / pdbf)
        same_el = [a[0] for a in p] == c.symbols
        dmax = float(np.abs(np.array([a[1:] for a in p]) - np.array(c.coords)).max())
        ok = (c.natm == natm == len(p) and same_el and dmax < 1e-9 and c.charge == chg and c.multiplicity == 1
              and c.method == 'B3LYP' and c.basis == '6-31G*' and c.keywords['grid'] == 'ultrafine')
        extra = {}
        if tag == 'large':
            ok = ok and c.radii == {'Zn': 1.456} and c.keywords['pop_mk'] and c.keywords['iop'] == {'6/33': 2, '6/42': 6}
            extra = dict(radii=c.radii, iop=c.keywords['iop'], n_frozen=sum(c.frozen))
        record(f'T1 .com reader, {tag} model', ok, natm=c.natm, charge=c.charge, multiplicity=c.multiplicity,
               elements_match_pdb=same_el, max_coordinate_difference_to_pdb=dmax, **extra)
    c = gaucom.read_com(HERE / 'mcpb' / 'CA2_small_fc.com')
    record('T1 .com reader, force-constant input', c.geometry_from_checkpoint and c.natm == 0 and
           c.keywords['opt_calcall'], route=c.route)
    # the flagged form (large_opt 1) is not written with our settings; test the parser on a synthetic file
    f = SCRATCH / 'flag.com'
    f.write_text('%Chk=x.chk\n# B3LYP/6-31G* Integral=(Grid=UltraFine) Opt Pop(MK,ReadRadii)\nIOp(6/33=2,6/42=6)\n \n'
                 'CLR\n \n1  1\nH  0    1.000    2.000   -3.000\nC -1   -11.123  -12.456  -13.789\n \n'
                 'Zn 1.456 !comment\n \nZn 1.456 !comment\n \n \n')
    c = gaucom.read_com(f)
    record('T1 .com reader, flagged atoms (large_opt 1 form)', c.frozen == [False, True] and
           c.coords[1] == (-11.123, -12.456, -13.789) and c.radii == {'Zn': 1.456}, frozen=c.frozen)


# ------------------------------------------------------------------------------------------------------- T2
def t2():
    rng = np.random.default_rng(SEED)
    worst = dict(crd=0.0, hess=0.0, own=0.0)
    for n in (2, 3, 7, 64, 101):
        z = rng.choice([1, 6, 7, 8, 16, 30], size=n)
        xyz = rng.uniform(-40, 40, size=(n, 3))                      # Bohr
        a = rng.normal(0, 0.3, size=(3 * n, 3 * n))
        h = 0.5 * (a + a.T)
        f = SCRATCH / f'rt_{n}.fchk'
        gauout.write_fchk(f, z, xyz, h, charge=1, multiplicity=1, energy=-1234.5, title='round trip')
        crd = np.array(gauio.get_crds_from_fchk(str(f), n)).reshape(n, 3)
        hm = gauio.get_matrix_from_fchk(str(f), 3 * n)
        # nine significant digits are written
        e1 = float(np.abs(crd - xyz).max() / 40)
        e2 = float(np.abs(hm - h).max() / np.abs(h).max())
        own = gauout.read_fchk_array(f, 'Cartesian Force Constants')
        e3 = float(np.abs(own - h[np.tril_indices(3 * n)]).max() / np.abs(h).max())
        worst = dict(crd=max(worst['crd'], e1), hess=max(worst['hess'], e2), own=max(worst['own'], e3))
    record('T2 fchk round trip through pymsmt, random arrays, 2-101 atoms',
           worst['crd'] < 1e-8 and worst['hess'] < 1e-8 and worst['own'] < 1e-8,
           largest_relative_error_coordinates=worst['crd'], largest_relative_error_hessian=worst['hess'],
           largest_relative_error_own_reader=worst['own'])
    # order of the triangle: element (i, j), i >= j, holds 1000*i + j
    n = 4
    h = np.zeros((3 * n, 3 * n))
    for i in range(3 * n):
        for j in range(i + 1):
            h[i, j] = h[j, i] = 1000 * i + j
    f = SCRATCH / 'order.fchk'
    gauout.write_fchk(f, [1] * n, np.zeros((n, 3)), h, charge=0, multiplicity=1)
    hm = gauio.get_matrix_from_fchk(str(f), 3 * n)
    first = [float(x) for x in pathlib.Path(f).read_text().split('Cartesian Force Constants')[1].split()[3:9]]
    record('T2 fchk lower-triangle order', bool(np.array_equal(hm, h)) and first == [0.0, 1000.0, 1001.0, 2000.0, 2001.0, 2002.0],
           first_six_numbers=first)
    # a matrix that is not symmetric must be refused
    try:
        bad = h.copy()
        bad[0, 5] += 1.0
        gauout.write_fchk(SCRATCH / 'bad.fchk', [1] * n, np.zeros((n, 3)), bad, charge=0, multiplicity=1)
        record('T2 fchk refuses a non-symmetric matrix', False)
    except ValueError as e:
        record('T2 fchk refuses a non-symmetric matrix', True, message=str(e))


# ------------------------------------------------------------------------------------------------------- T3
def amber_energy(x, bonds, angles):
    """kcal/mol; x in Angstrom (N,3); bonds (i, j, k, r0); angles (i, j, k, kth, th0 in rad). AMBER: k (q - q0)^2."""
    e = 0.0
    for i, j, k, r0 in bonds:
        e += k * (np.linalg.norm(x[i] - x[j]) - r0) ** 2
    for i, j, l, k, t0 in angles:
        u, v = x[i] - x[j], x[l] - x[j]
        t = np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))
        e += k * (t - t0) ** 2
    return e


def numeric_hessian(x0, bonds, angles, h=1e-3):
    n = x0.size
    hess = np.zeros((n, n))
    f = lambda v: amber_energy(v.reshape(-1, 3), bonds, angles)  # noqa: E731
    x = x0.ravel().copy()
    for i in range(n):
        for j in range(i + 1):
            xpp, xpm, xmp, xmm = x.copy(), x.copy(), x.copy(), x.copy()
            xpp[i] += h; xpp[j] += h
            xpm[i] += h; xpm[j] -= h
            xmp[i] -= h; xmp[j] += h
            xmm[i] -= h; xmm[j] -= h
            hess[i, j] = hess[j, i] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h * h)
    return hess                                                     # kcal/mol/A^2


def t3():
    # diatomic, bond not along an axis
    kb, r0 = 123.4, 2.015
    d = np.array([1.0, 2.0, -2.0]) / 3.0
    x = np.array([[0.3, -0.2, 0.1], [0.3, -0.2, 0.1] + r0 * d])
    hk = numeric_hessian(x, [(0, 1, kb, r0)], [])
    h_au = hk / HB2_TO_KCAL_MOL_A2
    f = SCRATCH / 'diatomic.fchk'
    gauout.write_fchk(f, [30, 7], x / gauout.BOHR, h_au, charge=2, multiplicity=1)
    crds = gauio.get_crds_from_fchk(str(f), 2)
    fcm = gauio.get_matrix_from_fchk(str(f), 6)
    dis, fc = get_bond_fc_with_sem([], crds, fcm, 1, 2, 1.0, 0)
    fc = float(np.real(fc))
    record('T3 units, diatomic through fchk and MCPB.py Seminario', abs(fc - kb) < 0.01 and abs(dis - r0) < 1e-5,
           force_constant_in=kb, force_constant_out=round(float(fc), 4), distance_in=r0, distance_out=round(float(dis), 6))
    # bent triatomic: two bonds and the angle between them
    kb1, kb2, ka, th0 = 150.0, 80.0, 40.0, np.radians(109.5)
    r1, r2 = 2.02, 1.93
    x = np.array([[r1, 0, 0], [0, 0, 0], [r2 * np.cos(th0), r2 * np.sin(th0), 0]])
    rot = np.linalg.qr(np.random.default_rng(SEED).normal(size=(3, 3)))[0]
    x = x @ rot.T + np.array([1.0, -2.0, 0.5])
    hk = numeric_hessian(x, [(0, 1, kb1, r1), (2, 1, kb2, r2)], [(0, 1, 2, ka, th0)])
    f = SCRATCH / 'triatomic.fchk'
    gauout.write_fchk(f, [7, 30, 7], x / gauout.BOHR, hk / HB2_TO_KCAL_MOL_A2, charge=2, multiplicity=1)
    crds = gauio.get_crds_from_fchk(str(f), 3)
    fcm = gauio.get_matrix_from_fchk(str(f), 9)
    d1, f1 = get_bond_fc_with_sem([], crds, fcm, 1, 2, 1.0, 0)
    d2, f2 = get_bond_fc_with_sem([], crds, fcm, 3, 2, 1.0, 0)
    ang, fa = get_ang_fc_with_sem([], crds, fcm, 1, 2, 3, 1.0, 0)
    f1, f2, fa = float(np.real(f1)), float(np.real(f2)), float(np.real(fa))
    # own implementation on the array in memory (no file)
    h_au, x_bohr = hk / (seminario.HARTREE_KCAL / seminario.BOHR ** 2), x / seminario.BOHR
    o1 = seminario.bond(h_au, x_bohr, 0, 1)
    o2 = seminario.bond(h_au, x_bohr, 2, 1)
    oa = seminario.angle(h_au, x_bohr, 0, 1, 2)
    ok = (abs(f1 / kb1 - 1) < 0.02 and abs(f2 / kb2 - 1) < 0.02 and
          abs(d1 - r1) < 1e-5 and abs(d2 - r2) < 1e-5 and abs(ang - 109.5) < 1e-3 and
          abs(f1 / o1[1] - 1) < 1e-5 and abs(f2 / o2[1] - 1) < 1e-5 and abs(fa / oa[1] - 1) < 1e-5)
    record('T3 bent triatomic with bond and angle terms: file route (pymsmt) equals array route (own code)', ok,
           bonds_in=[kb1, kb2], bonds_out_pymsmt_from_file=[round(f1, 3), round(f2, 3)],
           bonds_out_own_code_from_array=[round(o1[1], 3), round(o2[1], 3)],
           angle_constant_in=ka, angle_constant_pymsmt_from_file=round(fa, 3),
           angle_constant_own_code_from_array=round(oa[1], 3),
           note='the angle constant of the method is an estimate; here it is 21 % above the input',
           distances_out=[round(float(d1), 5), round(float(d2), 5)], angle_out=round(float(ang), 3),
           constants=dict(B_TO_A=B_TO_A, H_TO_KCAL_MOL=H_TO_KCAL_MOL))
    # one harmonic angle term alone: analytic Hessian K g g^T, closed form of Seminario's formula
    out = []
    for th_deg in (90.0, 109.5, 120.0):
        th = np.radians(th_deg)
        xa = np.array([[r1, 0, 0], [0, 0, 0], [r2 * np.cos(th), r2 * np.sin(th), 0]]) @ rot.T + np.array([1.0, -2.0, 0.5])
        eab, ecb = (xa[0] - xa[1]) / r1, (xa[2] - xa[1]) / r2
        ga = (np.cos(th) * eab - ecb) / (r1 * np.sin(th))           # d(theta)/d(r_a)
        gc = (np.cos(th) * ecb - eab) / (r2 * np.sin(th))
        g = np.concatenate([ga, -(ga + gc), gc])
        hk = 2.0 * ka * np.outer(g, g)                              # kcal/mol/A^2, E = ka (theta - theta0)^2
        f = SCRATCH / f'angle_{int(th_deg)}.fchk'
        gauout.write_fchk(f, [7, 30, 7], xa / gauout.BOHR, hk / HB2_TO_KCAL_MOL_A2, charge=2, multiplicity=1)
        a_out, k_out = get_ang_fc_with_sem([], gauio.get_crds_from_fchk(str(f), 3),
                                           gauio.get_matrix_from_fchk(str(f), 9), 1, 2, 3, 1.0, 0)
        c, rho = -np.cos(th), r1 / r2
        expect = ka / (1.0 / (1.0 + rho * c) + 1.0 / (1.0 + c / rho))
        out.append(dict(angle_in=th_deg, angle_out=round(float(a_out), 4), k_in=ka,
                        k_closed_form=round(float(expect), 4), k_pymsmt_from_file=round(float(np.real(k_out)), 4)))
    ok = all(abs(o['k_pymsmt_from_file'] / o['k_closed_form'] - 1) < 1e-5 and abs(o['angle_out'] - o['angle_in']) < 1e-3
             for o in out)
    record('T3 units, one harmonic angle term through fchk and MCPB.py Seminario against the closed form', ok,
           cases=out)


# ------------------------------------------------------------------------------------------------------- T4
def t4():
    rng = np.random.default_rng(SEED + 1)
    for tag, nat, npt in (('small', 5, 37), ('real size', 101, 45000), ('99,999 points, the limit of the .esp header', 12, 99999)):
        xyz = rng.uniform(-30, 30, size=(nat, 3))
        xyz[0] = (-12.345678, -23.456789, -34.567891)                # three fields that touch
        pts = rng.uniform(-45, 45, size=(npt, 3))
        pts[0] = (-99.123456, -10.000001, -45.5)
        va = rng.uniform(-60, -1, size=nat)
        vp = rng.normal(0.1, 0.08, size=npt)
        log = SCRATCH / f'mk_{nat}_{npt}.log'
        gauout.write_mk_log(log, ['C'] * nat, xyz, pts, va, vp, header_lines=['synthetic test'])
        # (1) pymsmt
        esp1 = SCRATCH / f'mk_{nat}_{npt}.pymsmt.esp'
        gauio.get_esp_from_gau(str(log), str(esp1))
        a1, p1, v1 = gauout.read_esp_file(esp1)
        # (2) espgen (antechamber)
        esp2 = SCRATCH / f'mk_{nat}_{npt}.espgen.esp'
        p = subprocess.run([f'{AMBERHOME}/bin/espgen', '-i', str(log), '-o', str(esp2)], capture_output=True,
                           text=True, env=dict(os.environ, AMBERHOME=AMBERHOME))
        res = {}
        for name, (a, q, v) in (('pymsmt', (a1, p1, v1)),):
            res[name] = dict(n_atoms=len(a), n_points=len(q),
                             max_err_atoms_bohr=float(np.abs(a - xyz / B_TO_A).max()),
                             max_err_points_bohr=float(np.abs(q - pts / B_TO_A).max()),
                             max_err_potential=float(np.abs(v - vp).max()))
        ok = (res['pymsmt']['n_atoms'] == nat and res['pymsmt']['n_points'] == npt and
              res['pymsmt']['max_err_atoms_bohr'] < 2e-6 and res['pymsmt']['max_err_points_bohr'] < 2e-5 and
              res['pymsmt']['max_err_potential'] < 6e-7)
        if p.returncode == 0 and esp2.exists() and esp2.stat().st_size > 0:
            try:
                a2, p2, v2 = gauout.read_esp_file(esp2)
                res['espgen'] = dict(n_atoms=len(a2), n_points=len(p2),
                                     max_err_atoms_bohr=float(np.abs(a2 - xyz / gauout.BOHR).max()),
                                     max_err_points_bohr=float(np.abs(p2 - pts / gauout.BOHR).max()),
                                     max_err_potential=float(np.abs(v2 - vp).max()),
                                     max_difference_to_pymsmt_file=float(max(np.abs(a2 - a1).max(), np.abs(p2 - p1).max(),
                                                                             np.abs(v2 - v1).max())))
                ok = ok and res['espgen']['n_points'] == npt and res['espgen']['max_err_potential'] < 6e-7 and \
                    res['espgen']['max_err_points_bohr'] < 2e-5
            except Exception as e:                                      # noqa: BLE001
                res['espgen'] = dict(error=f'{type(e).__name__}: {e}')
                ok = False
        else:
            res['espgen'] = dict(error=(p.stdout + p.stderr)[-300:])
            ok = False
        record(f'T4 MK log round trip, {tag} ({nat} atoms, {npt} points)', ok, **res)


    try:
        gauout.write_mk_log(SCRATCH / 'toomany.log', ['C'], [[0, 0, 0]], np.zeros((100000, 3)) + 5.0, [0.0],
                            np.zeros(100000))
        record('T4 MK log writer refuses more than 99,999 points', False)
    except ValueError as e:
        record('T4 MK log writer refuses more than 99,999 points', True, message=str(e))


# ------------------------------------------------------------------------------------------------------- T5
def t5():
    c = gaucom.read_com(HERE / 'mcpb' / 'CA2_large_mk.com')
    xyz = np.array(c.coords)
    pts, lay, own = gauout.mk_points(c.symbols, xyz, c.radii, density=float(c.keywords['iop']['6/42']))
    rad = gauout.mk_radii(c.symbols, c.radii)
    d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2)
    ratio = d / rad[None, :]
    # every point lies on the scaled sphere of its own atom and outside all others
    on_own = float(np.abs(ratio[np.arange(len(pts)), own] - lay).max())
    inside = int(((ratio < lay[:, None] - 1e-9).sum(axis=1)).max())
    counts = {str(f): int((lay == f).sum()) for f in gauout.MK_LAYERS}
    # an isolated atom: number of points against density * area
    p1, _, _ = gauout.mk_points(['C'], [[0, 0, 0]], density=6.0)
    expect = sum(6.0 * 4 * np.pi * (1.5 * f) ** 2 for f in gauout.MK_LAYERS)
    record('T5 Merz-Kollman points of the large model', on_own < 1e-9 and inside == 0 and len(pts) < 99999 and
           abs(len(p1) / expect - 1) < 0.05,
           n_points=int(len(pts)), per_layer=counts, radii_used={s: float(r) for s, r in zip(c.symbols, rad)},
           density_per_A2=float(c.keywords['iop']['6/42']), max_deviation_from_own_sphere=on_own,
           nearest_point_to_any_atom_in_units_of_its_radius=float(ratio.min()),
           single_atom_points=int(len(p1)), single_atom_expected_from_area=round(float(expect), 1))


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    t1(); t2(); t3(); t4(); t5()
    out = dict(seed=SEED, all_passed=all(r['passed'] for r in results), results=results,
               readers=dict(pymsmt='pymsmt.mol.gauio of AmberTools 26 (MCPB.py version 2026)',
                            espgen=f'{AMBERHOME}/bin/espgen'))
    (HERE / 'test_adapters.json').write_text(json.dumps(out, indent=1) + '\n')
    print('ALL PASSED' if out['all_passed'] else 'FAILURES')


if __name__ == '__main__':
    main()
