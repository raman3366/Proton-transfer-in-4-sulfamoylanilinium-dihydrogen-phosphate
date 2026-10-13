"""Test of md/qm/after_qm.sh on the REAL models with adapter files whose content is known beforehand.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/test_after_qm.py

No quantum chemistry.  In a scratch copy of md/qm/mcpb two files are written with md/qm/gauout.py:

  CA2_large_mk.log   the potential of KNOWN point charges on the 101 atoms of the large model, at the
                     Merz-Kollman points the real job will use.  The charges are the ff14SB charges of the
                     residues and caps, the AM1-BCC charges of the ligand and +2 on the zinc (sum +1).  They
                     obey every constraint MCPB.py step 3b puts on the fit (caps neutral, backbone heavy
                     atoms at ff14SB values, equivalent hydrogens equal), so the two-stage RESP fit must give
                     them back, up to the pull of its hyperbolic restraint.
  CA2_small_opt.fchk the Hessian of a KNOWN harmonic model of the 64-atom small model at its input geometry,
                     written analytically, in two variants:
                     A  bond terms only: every covalent bond 300, the four Zn-N bonds 80, 75, 70 and 120
                        kcal/mol/A^2.  With bond terms only Seminario's bond constant is exact, so MCPB.py
                        step 2s must print these four numbers and the four distances of the structure.
                        (The angle constants are then zero; this variant runs step 2s only.)
                     B  the same bonds plus a harmonic term of 40 kcal/mol/rad^2 on each of the 14 angles
                        with zinc.  This variant runs the whole pipeline; the bond constants then come
                        back within a few per cent (the angle terms leak into the projection) and the
                        angle constants are positive (Seminario's angle value is an estimate, see
                        md/qm/test_adapters.py T3).

Then md/qm/after_qm.sh runs in that copy with variant B (2s) and without the fchk (2e): MCPB.py steps 2,
3, 4, tleap, the checks.  Passes when
  - both runs end with all hard checks passed,
  - variant A gives the four bond constants within 0.05 and the distances within 0.002 A,
  - resp, run on the files MCPB.py wrote but with the restraint weight set to zero, gives the known
    charges back within 0.005 e on all 101 atoms (this shows that order, units and constraints are right);
    the charges of MCPB.py's own run, with the restraint, are reported beside them: the restraint pulls the
    buried zinc and its four nitrogens by up to about 0.1 e, which is the method and not a fault.
The minimisation in the check uses the SYNTHETIC charges (zinc +1.9) and says nothing about the real model.

Result: md/qm/test_after_qm.json.  Scratch: SCRATCH below.  The real folder md/qm/mcpb is not touched.
"""
import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import parmed

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import gauout  # noqa: E402

SCRATCH = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                       '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/mdprep/zinc/after_qm_test')
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
HB2 = 627.5095 / 0.529177249 ** 2          # the constants MCPB.py uses
KZN = {'HID93@NE2': 80.0, 'HID95@NE2': 75.0, 'HIE118@ND1': 70.0, 'SAN261@N2': 120.0}
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'Zn': 1.22}


def pdb(path):
    out = []
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            out.append(dict(serial=int(l[6:11]), name=l[12:16].strip(), resn=l[17:20].strip(), resi=int(l[22:26]),
                            xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def known_charges(work):
    """ff14SB / AM1-BCC / +2 charges for the atoms of the large model, in its order."""
    lib = {}
    for a in parmed.load_file(f'{AMBERHOME}/dat/pymsmt/parm12.mol2', structure=True).atoms:
        lib[(a.residue.name, a.name)] = a.charge
    for a in parmed.load_file(str(work / 'SAN.mol2'), structure=True).atoms:
        lib[('SAN', a.name)] = a.charge
    lib[('ZN', 'ZN')] = 2.0
    q = []
    for a in pdb(work / 'CA2_large.pdb'):
        q.append(lib[(a['resn'], a['name'])])
    return np.array(q)


def bond_hessian(x_ang, bonds):
    """E = sum k (r - r0)^2 at r = r0: H = 2 k u u^T on the pair blocks.  kcal/mol/A^2."""
    n = len(x_ang)
    h = np.zeros((3 * n, 3 * n))
    for i, j, k in bonds:
        u = (x_ang[j] - x_ang[i]) / np.linalg.norm(x_ang[j] - x_ang[i])
        b = 2.0 * k * np.outer(u, u)
        h[3 * i:3 * i + 3, 3 * i:3 * i + 3] += b
        h[3 * j:3 * j + 3, 3 * j:3 * j + 3] += b
        h[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= b
        h[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= b
    return h


def angle_hessian(x_ang, angles):
    """E = sum k (theta - theta0)^2 at theta = theta0: H = 2 k g g^T, g = d(theta)/dx.  kcal/mol/A^2."""
    n = len(x_ang)
    h = np.zeros((3 * n, 3 * n))
    for a, b, c, k in angles:
        eab, ecb = x_ang[a] - x_ang[b], x_ang[c] - x_ang[b]
        r1, r2 = np.linalg.norm(eab), np.linalg.norm(ecb)
        eab, ecb = eab / r1, ecb / r2
        ct = float(np.dot(eab, ecb))
        st = np.sqrt(1.0 - ct * ct)
        ga = (ct * eab - ecb) / (r1 * st)
        gc = (ct * ecb - eab) / (r2 * st)
        g = np.zeros(3 * n)
        g[3 * a:3 * a + 3], g[3 * c:3 * c + 3], g[3 * b:3 * b + 3] = ga, gc, -(ga + gc)
        h += 2.0 * k * np.outer(g, g)
    return h


def frcmod_zn_bonds(path):
    out = {}
    sec = None
    for l in pathlib.Path(path).read_text().splitlines():
        if l.strip() in ('BOND', 'ANGL'):
            sec = l.strip()
            continue
        if sec == 'BOND' and 'M1' in l[:5]:
            t = l[:5].split('-')
            other = [x for x in t if x != 'M1'][0]
            out[other] = (float(l[5:].split()[0]), float(l[5:].split()[1]))
    return out


def main():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    work = SCRATCH / 'mcpb'
    shutil.copytree(HERE / 'mcpb', work)
    for f in list(work.glob('*.fchk')) + list(work.glob('CA2_large_mk.log')):
        f.unlink()
    res = {}
    # ---- the potential of known charges ----------------------------------------------------------------------
    com = gaucom.read_com(work / 'CA2_large_mk.com')
    xyz = np.array(com.coords)
    q = known_charges(work)
    pts, lay, own = gauout.mk_points(com.symbols, xyz, com.radii, density=float(com.keywords['iop']['6/42']))
    d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2) / gauout.BOHR
    v = (q[None, :] / d).sum(axis=1)
    gauout.write_mk_log(work / 'CA2_large_mk.log', com.symbols, xyz, pts, np.zeros(len(xyz)), v,
                        header_lines=['SYNTHETIC: potential of known point charges, test of after_qm.sh'])
    res['synthetic_esp'] = dict(n_points=int(len(pts)), sum_of_charges=round(float(q.sum()), 4))
    # ---- the Hessian of a known harmonic model ----------------------------------------------------------------
    sm = gaucom.read_com(work / 'CA2_small_opt.com')
    x = np.array(sm.coords)
    sp = pdb(work / 'CA2_small.pdb')
    izn = sm.symbols.index('Zn')
    bonds, zn_terms = [], {}
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            if izn in (i, j):
                continue
            if np.linalg.norm(x[i] - x[j]) <= COV[sm.symbols[i]] + COV[sm.symbols[j]] + 0.4:
                bonds.append((i, j, 300.0))
    for i, a in enumerate(sp):
        key = f"{a['resn']}{a['resi']}@{a['name']}"
        if key in KZN:
            bonds.append((i, izn, KZN[key]))
            zn_terms[key] = dict(k_in=KZN[key], r_in=round(float(np.linalg.norm(x[i] - x[izn])), 4))
    hA = bond_hessian(x, bonds) / HB2
    zatoms = [i for i, a in enumerate(sp) if f"{a['resn']}{a['resi']}@{a['name']}" in KZN]
    nb = {i: [j for j in range(len(x)) if j not in (i, izn) and
              np.linalg.norm(x[i] - x[j]) <= COV[sm.symbols[i]] + COV[sm.symbols[j]] + 0.4] for i in zatoms}
    angles = [(zatoms[i], izn, zatoms[j], 40.0) for i in range(4) for j in range(i + 1, 4)]
    angles += [(izn, i, j, 40.0) for i in zatoms for j in nb[i]]
    assert len(angles) == 14
    hB = hA + angle_hessian(x, angles) / HB2
    z = [gauout.Z_OF[s_] for s_ in sm.symbols]
    # variant A: step 2s alone, in its own copy
    workA = SCRATCH / 'mcpb_A'
    shutil.copytree(work, workA)
    gauout.write_fchk(workA / 'CA2_small_opt.fchk', z, x / 0.529177249, hA, charge=sm.charge,
                      multiplicity=sm.multiplicity, title='SYNTHETIC A: bond terms only, known constants')
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f"{AMBERHOME}/bin:" + os.environ.get('PATH', ''))
    p = subprocess.run([f'{AMBERHOME}/bin/python', str(HERE / 'mcpb_step2s.py'), 'CA2.in'], cwd=workA, env=env,
                       capture_output=True, text=True)
    (workA / 'step2s.log').write_text(p.stdout + p.stderr)
    types = {'Y1': 'HID93@NE2', 'Y2': 'HID95@NE2', 'Y3': 'HIE118@ND1', 'Y4': 'SAN261@N2'}
    got = frcmod_zn_bonds(workA / 'CA2_mcpbpy.frcmod') if (workA / 'CA2_mcpbpy.frcmod').exists() else {}
    cmpA = [dict(type=t, atom=types[t], k_in=zn_terms[types[t]]['k_in'], k_out=got.get(t, (None, None))[0],
                 r_in=zn_terms[types[t]]['r_in'], r_out=got.get(t, (None, None))[1]) for t in types]
    res['variant_A_bonds_only'] = dict(
        returncode=p.returncode, terms=cmpA, guard=(p.stdout.strip().splitlines() or [''])[-1],
        passed=bool(p.returncode == 0 and all(c['k_out'] is not None and abs(c['k_in'] - c['k_out']) < 0.05 and
                                              abs(c['r_in'] - c['r_out']) < 2e-3 for c in cmpA)))
    # the same step with plain MCPB.py: does it stop, as found?
    p0 = subprocess.run([f'{AMBERHOME}/bin/python', f'{AMBERHOME}/bin/MCPB.py', '-i', 'CA2.in', '-s', '2s'], cwd=workA,
                        env=env, capture_output=True, text=True)
    res['plain_mcpb_2s_on_variant_A'] = dict(returncode=p0.returncode,
                                             last_line=((p0.stdout + p0.stderr).strip().splitlines() or [''])[-1])
    # variant B for the pipeline
    gauout.write_fchk(work / 'CA2_small_opt.fchk', z, x / 0.529177249, hB, charge=sm.charge,
                      multiplicity=sm.multiplicity, title='SYNTHETIC B: bonds and the 14 angles with zinc')
    # ---- the pipeline, twice ------------------------------------------------------------------------------------
    for step2 in ('2s', '2e'):
        rep = SCRATCH / f'report_{step2}'
        p = subprocess.run(['bash', str(HERE / 'after_qm.sh'), step2, str(work), str(rep)], capture_output=True,
                           text=True)
        r = dict(returncode=p.returncode, tail=(p.stdout + p.stderr).splitlines()[-12:])
        f = rep / f'after_qm_check_{step2}.json'
        if f.exists():
            chk = json.loads(f.read_text())
            r['all_hard_checks_passed'] = chk['all_hard_checks_passed']
            r['hard'] = [dict(check=c['check'], passed=c['passed']) for c in chk['hard']]
            r['zinc_bonds'] = chk['reported']['zinc_bonds']
            r['zinc_charge'] = chk['reported']['zinc_charge']
            r['rules_of_thumb'] = chk['reported']['mcpb_rules_of_thumb']
            r['first_look_minimisation'] = chk['reported']['first_look_minimisation']
            if step2 == '2s':
                gotb = {b['atom']: b for b in chk['reported']['zinc_bonds']}
                names = {'HD1@NE2': 'HID93@NE2', 'HD2@NE2': 'HID95@NE2', 'HE1@ND1': 'HIE118@ND1', 'SN1@N2': 'SAN261@N2'}
                r['variant_B_bonds'] = [dict(atom=k, k_in=zn_terms[names[k]]['k_in'], k_out=gotb[k]['k']) for k in names]
                r['variant_B_angles'] = [dict(atoms=a['atoms'], k_in=40.0, k_out=a['k']) for a in chk['reported']['zinc_angles']]
                # RESP: (1) MCPB.py's own charges (with the restraint) against the known charges
                order = [(a['resi'], a['name']) for a in pdb(work / 'CA2_large.pdb')]
                q_mcpb = np.array([float(t) for t in (work / 'resp2.chg').read_text().split()])
                dev = np.abs(q_mcpb - q)
                worst = np.argsort(-dev)[:5]
                r['resp_with_restraint'] = dict(
                    largest_deviation=round(float(dev.max()), 4), rms_deviation=round(float(np.sqrt((dev ** 2).mean())), 4),
                    worst_five=[dict(residue=order[i][0], atom=order[i][1], resp=round(float(q_mcpb[i]), 4),
                                     known=round(float(q[i]), 4)) for i in worst])
                # (2) the same inputs with the restraint weight zero
                for f in ('resp1.in', 'resp2.in'):
                    t = (work / f).read_text()
                    t = t.replace('qwt = 0.00050,', 'qwt = 0.00000,').replace('qwt = 0.00100,', 'qwt = 0.00000,')
                    (work / f.replace('.in', '_q0.in')).write_text(t)
                env = dict(os.environ, AMBERHOME=AMBERHOME)
                subprocess.run([f'{AMBERHOME}/bin/resp', '-O', '-i', 'resp1_q0.in', '-o', 'resp1_q0.out', '-p',
                                'resp1_q0.pch', '-t', 'resp1_q0.chg', '-e', 'CA2_large_mk.esp', '-s', 'resp1_q0.esp'],
                               cwd=work, env=env, capture_output=True)
                subprocess.run([f'{AMBERHOME}/bin/resp', '-O', '-i', 'resp2_q0.in', '-o', 'resp2_q0.out', '-p',
                                'resp2_q0.pch', '-q', 'resp1_q0.chg', '-t', 'resp2_q0.chg', '-e', 'CA2_large_mk.esp',
                                '-s', 'resp2_q0.esp'], cwd=work, env=env, capture_output=True)
                q0 = np.array([float(t) for t in (work / 'resp2_q0.chg').read_text().split()])
                dev0 = np.abs(q0 - q)
                r['resp_without_restraint'] = dict(n_atoms=int(len(q0)), largest_deviation=round(float(dev0.max()), 5),
                                                   zinc=round(float(q0[com.symbols.index('Zn')]), 4))
                r['resp_known_answer_passed'] = bool(len(q0) == len(q) and dev0.max() < 0.005)
        res[step2] = r
    res['all_passed'] = bool(all(res[s].get('all_hard_checks_passed') for s in ('2s', '2e')) and
                             res['variant_A_bonds_only']['passed'] and res['2s'].get('resp_known_answer_passed'))
    (HERE / 'test_after_qm.json').write_text(json.dumps(res, indent=1) + '\n')
    for s in ('2s', '2e'):
        print(s, 'return code', res[s]['returncode'], 'hard checks', res[s].get('hard'))
        if res[s]['returncode'] != 0:
            print('\n'.join(res[s]['tail']))
    print('variant A', res['variant_A_bonds_only'])
    print('plain MCPB.py 2s', res['plain_mcpb_2s_on_variant_A'])
    print('variant B bonds', res['2s'].get('variant_B_bonds'))
    print('variant B angles', res['2s'].get('variant_B_angles'))
    print('resp with restraint', res['2s'].get('resp_with_restraint'))
    print('resp without restraint', res['2s'].get('resp_without_restraint'))
    print('ALL PASSED' if res['all_passed'] else 'FAILURES')


if __name__ == '__main__':
    main()
