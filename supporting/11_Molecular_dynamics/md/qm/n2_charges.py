"""Gate G-R2, the declared second step (md/PLAN.md, 6 Oct 2026, 1:34 PM IST, point 2): charges fitted to the
electrostatic potentials of TWO placements of the large model, X and the R-like one, instead of one.

    MD = AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python ;  QM = /root/cqc-sfa-tools/env/bin/python
    MD md/qm/n2_charges.py prepare-esp     the input of the potential job for the R-like placement (no energy):
                                           md/qm/n2/esp_Rx/in/CA2_large_mk.com, config.json, prepare.json
    QM md/qm/run_qm.py md/qm/n2/esp_Rx/config.json --jobs esp,write      (rented machine, 16 threads; task
                                           n2_esp_Rx of md/gpu/tasks_n2.json) -> md/qm/n2/esp_Rx/in/CA2_large_mk.log
    MD md/qm/n2_charges.py resp2 --work SCRATCH [--second X]   the two-stage RESP fit of MCPB.py with the equal
                                           ligand atoms (the inputs of md/qm/resp_equiv.py), by the program resp, on
                                           the two potentials together -> md/qm/n2/charges2/two_placements/.
                                           --second X: the same potential twice, which must give the charges of
                                           the one-placement fit (the proof of the two-molecule input)
    MD md/qm/n2_charges.py mol2 [--second X]   the five mol2 files with these charges, in a folder for
                                           md/build_topology.py --mcpb-dir -> md/qm/n2/charges2/mcpb_2e_two/
    MD md/qm/n2_charges.py ff --work SCRATCH [--second X] [--name ff2]   the force-field copy of the small model with
                                           them (cap scheme "methyl", as gate G-R) and the energies of every point
                                           -> md/qm/n2/ff2/; the fit and the gate use it with N2_FF=ff2

THE R-LIKE PLACEMENT is point Rx of gate G-R (md/qm/gr/points.json): the ligand of the crystal pose placed rigidly
as in the rotated start, Zn-N2 1.931 A, hydrogen by the construction of the starts.  The large model is the one
MCPB.py wrote (md/qm/mcpb/CA2_large_mk.com, 101 atoms, frame of 6RL9, the same frame as the small model); only
its last 18 atoms, the ligand, are replaced.  The potential job is md/qm/run_qm.py unchanged, with the settings of
md/qm/config_prod.json (the job that gave the potential of placement X): same level, same grids, same rule for the
Merz-Kollman points.

Computing the potential commits to nothing: whether the second step is taken is decided by the rule of md/N2.md,
section 1, from the fitting set alone.
"""
import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np

import n2_common as nc

MCPB = nc.HERE / 'mcpb'
OUT = nc.N2 / 'esp_Rx'
LIG_ORDER = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S', 'H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21']
BULK_VM = '/root/cqc-sfa-md/qm/n2/esp_Rx'


def read_large():
    lines = (MCPB / 'CA2_large_mk.com').read_text().splitlines()
    k = next(i for i, l in enumerate(lines) if l.split() == ['1', '1'])
    atoms, j = [], k + 1
    while lines[j].strip():
        t = lines[j].split()
        atoms.append((t[0], [float(v) for v in t[1:4]]))
        j += 1
    return lines, k + 1, j, atoms


def cmd_prepare_esp(args):
    lines, a0, a1, atoms = read_large()
    assert len(atoms) == 101 and atoms[-19][0] == 'Zn', 'the large model is not the one this was written for'
    pdb = [l for l in (MCPB / 'CA2_large.pdb').read_text().splitlines() if l.startswith(('ATOM', 'HETATM'))]
    names = [l[12:16].strip() for l in pdb]
    assert len(pdb) == 101 and names[-18:] == LIG_ORDER, 'order of the ligand atoms in CA2_large.pdb'
    pts = {p['id']: p for p in json.loads((nc.GR / 'points.json').read_text())['points']}
    x_lig = np.array([pts['X']['ligand_coordinates'][n] for n in LIG_ORDER])
    rx = np.array([pts['Rx']['ligand_coordinates'][n] for n in LIG_ORDER])
    old = np.array([a[1] for a in atoms[-18:]])
    # the ligand of the large model is the ligand of point X (the file has three decimals)
    dx = float(np.abs(old - x_lig).max())
    assert dx < 6e-4, f'the ligand of CA2_large_mk.com is not point X: {dx}'
    # the fixed atoms of the small model are where the large model has them (same frame)
    g = nc.Geometry()
    zn_large = np.array(atoms[-19][1])
    dzn = float(np.linalg.norm(zn_large - g.m.zn))
    assert dzn < 1e-3, f'zinc of the large model against the small model: {dzn}'
    el_lig = [a[0] for a in atoms[-18:]]
    assert el_lig == [n[0] for n in LIG_ORDER], el_lig
    new = lines[:a1 - 18] + [f'{e:<2s} {v[0]:15.8f} {v[1]:15.8f} {v[2]:15.8f}' for e, v in zip(el_lig, rx)] + lines[a1:]
    (OUT / 'in').mkdir(parents=True, exist_ok=True)
    (OUT / 'in' / 'CA2_large_mk.com').write_text('\n'.join(new) + '\n')
    shutil.copyfile(MCPB / 'CA2_small_opt.com', OUT / 'in' / 'CA2_small_opt.com')      # run_qm.py reads both inputs
    # contacts of the R-like ligand with the atoms of the large model
    rest = np.array([a[1] for a in atoms[:-19]])
    rest_el = [a[0] for a in atoms[:-19]]
    con = {}
    for tag, lig in (('X', old), ('Rx', rx)):
        best = dict(heavy_heavy=[9e9, ''], heavy_h=[9e9, ''], h_h=[9e9, ''])
        for i, n in enumerate(LIG_ORDER):
            d = np.linalg.norm(rest - lig[i], axis=1)
            for j in range(len(rest)):
                nh = (n[0] == 'H') + (rest_el[j] == 'H')
                t = best[('heavy_heavy', 'heavy_h', 'h_h')[nh]]
                if d[j] < t[0]:
                    t[0], t[1] = float(d[j]), f'{n}...{pdb[j][17:20].strip()}{pdb[j][22:26].strip()}:{names[j]}'
        con[tag] = {k: dict(d=round(v[0], 3), pair=v[1]) for k, v in best.items()}
    cfg = json.loads((nc.HERE / 'config_prod.json').read_text())
    cfg.update(indir='md/qm/n2/esp_Rx/in', bulk=BULK_VM, outdir='md/qm/n2/esp_Rx/results', threads=16, jobs=['esp', 'write'])
    nc.write_json(OUT / 'config.json', cfg)
    rec = dict(what='input of the potential job of the large model with the ligand in the R-like placement (point Rx of gate '
                    'G-R); written by md/qm/n2_charges.py prepare-esp; no energy', written=nc.now(),
               source=dict(com='md/qm/mcpb/CA2_large_mk.com', com_sha256=nc.sha(MCPB / 'CA2_large_mk.com'),
                           points='md/qm/gr/points.json', points_sha256=nc.sha(nc.GR / 'points.json'),
                           config='md/qm/config_prod.json', config_sha256=nc.sha(nc.HERE / 'config_prod.json')),
               natoms=101, ligand_atoms='84-101', ligand_order=LIG_ORDER,
               ligand_of_the_large_model_against_point_X_A=dx, zinc_large_against_small_A=dzn,
               ligand_moved_rms_A=round(float(np.sqrt(((rx - old) ** 2).sum(1).mean())), 3),
               closest_contacts_ligand_to_the_other_atoms_of_the_large_model=con,
               settings_changed_against_config_prod=dict(indir=cfg['indir'], bulk=cfg['bulk'], outdir=cfg['outdir'],
                                                         threads=16, jobs=cfg['jobs']),
               files={f.name: nc.sha(f) for f in sorted((OUT / 'in').iterdir()) if f.suffix == '.com'},
               config_sha256=nc.sha(OUT / 'config.json'))
    nc.write_json(OUT / 'prepare.json', rec)
    print(json.dumps({k: rec[k] for k in ('ligand_of_the_large_model_against_point_X_A', 'zinc_large_against_small_A',
                                          'ligand_moved_rms_A', 'closest_contacts_ligand_to_the_other_atoms_of_the_large_model')},
                     indent=1))
    return 0


# ------------------------------------------------------------------------------- the fit to two potentials
def two_molecule_input(src, dst):
    """A resp input of MCPB.py (one molecule) turned into the same fit on two molecules with the same atoms.

    Layout, found by trial with this build of resp and proven by the run on the same potential twice (which must
    and does give the charges of the one-molecule fit to the last digit): molecule 1 with its atom lines and its
    group constraints, a blank line; molecule 2 with the same atom lines, two blank lines; then the pairs
    (molecule 1, atom i) = (molecule 2, atom i) for every free centre that is not already tied to another atom
    inside its molecule, except ONE (its charge in molecule 2 follows from the total charge of that molecule;
    with it the two total-charge conditions would say the same thing twice)."""
    lines = pathlib.Path(src).read_text().splitlines()
    k = [i for i, l in enumerate(lines) if l.strip() == '&end'][0] + 3
    natm = int(lines[k].split()[1])
    head, atoms, cons = lines[:k + 1], lines[k + 1:k + 1 + natm], [l for l in lines[k + 1 + natm:] if l.strip()]
    assert sum('nmol = 1,' in l for l in head) == 1
    cntrl = [l.replace('nmol = 1,', 'nmol = 2,') for l in head[:-3]]
    mol = head[-3:]
    free = [i for i, l in enumerate(atoms, 1) if int(l.split()[1]) == 0]
    left_out = free[-1]
    pairs = []
    for i in free:
        if i != left_out:
            pairs += ['%5d' % 2, '%5d%5d%5d%5d' % (1, i, 2, i)]
    out = cntrl + mol + atoms + cons + [''] + mol + atoms + ['', ''] + pairs + ['', '', '']
    pathlib.Path(dst).write_text('\n'.join(out) + '\n')
    return dict(atoms=natm, free_centres=len(free), tied_between_the_molecules=len(free) - 1, left_out=left_out,
                group_constraints=len(cons) // 2)


def potential_rrms(q, esp):
    """Relative r.m.s. of the potential of point charges q against the potential on the fit points."""
    at, pts, v = esp
    d = np.linalg.norm(pts[:, None, :] - at[None, :, :], axis=2)
    vc = (q[None, :] / d).sum(1)
    return float(np.sqrt(((v - vc) ** 2).sum() / (v ** 2).sum()))


def cmd_resp2(args):
    import resp_equiv as rq
    src = pathlib.Path(args.mcpb_dir)
    work = pathlib.Path(args.work)
    tag = 'two_placements' if args.second == 'Rx' else 'check_same_potential_twice'
    out = nc.N2 / 'charges2' / tag
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    for f in ('resp1_equiv.in', 'resp2.in', 'CA2_large_mk.esp', 'resp1_equiv.chg', 'resp2_equiv.chg', 'CA2_large.pdb',
              'CA2_standard.pdb'):
        shutil.copyfile(src / f, work / f)
    env = dict(os.environ, AMBERHOME=nc.AMBERHOME)
    if args.second == 'Rx':
        log = OUT / 'in' / 'CA2_large_mk.log'
        if not log.exists():
            raise SystemExit(f'{log} is not there: the potential job of the R-like placement has not been copied from the machine')
        p = subprocess.run([f'{nc.AMBERHOME}/bin/espgen', '-i', str(log), '-o', str(work / 'second.esp')], env=env,
                           capture_output=True, text=True)
        if p.returncode != 0 or not (work / 'second.esp').exists():
            raise SystemExit(f'espgen failed: {p.stdout[-300:]} {p.stderr[-300:]}')
    else:
        shutil.copyfile(work / 'CA2_large_mk.esp', work / 'second.esp')
    e1, e2 = rq.read_esp(work / 'CA2_large_mk.esp'), rq.read_esp(work / 'second.esp')
    order = rq.pdb_atoms(work / 'CA2_large.pdb')
    assert len(e1[0]) == len(e2[0]) == len(order) == 101
    # the atoms of the two potentials: the same 83 fixed atoms; the ligand as in X and as in Rx
    pts = {p['id']: p for p in json.loads((nc.GR / 'points.json').read_text())['points']}
    bohr = 0.529177249
    lig = [k for k, (resi, resn, name) in enumerate(order) if resn == rq.LIG]
    names = [order[k][2] for k in lig]
    geo = {}
    for label, e, pid in (('first', e1, 'X'), ('second', e2, args.second)):
        want = np.array([pts[pid]['ligand_coordinates'][n] for n in names])
        geo[label] = dict(point=pid, ligand_against_the_point_A=float(np.abs(e[0][lig] * bohr - want).max()), n_fit_points=len(e[2]))
        assert geo[label]['ligand_against_the_point_A'] < 2e-3, geo
    fixed = [k for k in range(101) if k not in lig]
    geo['fixed_atoms_largest_difference_A'] = float(np.abs(e1[0][fixed] - e2[0][fixed]).max() * bohr)
    assert geo['fixed_atoms_largest_difference_A'] < 1e-5
    (work / 'two.esp').write_text((work / 'CA2_large_mk.esp').read_text() + (work / 'second.esp').read_text())
    lay1 = two_molecule_input(work / 'resp1_equiv.in', work / 'two1.in')
    lay2 = two_molecule_input(work / 'resp2.in', work / 'two2.in')
    s1 = rq.run_resp(work, 'two1.in', 'two1', 'two.esp')
    s2 = rq.run_resp(work, 'two2.in', 'two2', 'two.esp', qin='two1.chg')
    q_all = rq.read_chg(work / 'two2.chg')
    assert len(q_all) == 202
    q, qb = q_all[:101], q_all[101:]
    q_one = rq.read_chg(work / 'resp2_equiv.chg')
    q_free = rq.read_chg(src / 'resp2_free.chg') if (src / 'resp2_free.chg').exists() else None
    rec = dict(what='RESP charges of the large model fitted to the potentials of two placements of the ligand with the '
                    'program resp (two molecules with the same atoms); written by md/qm/n2_charges.py resp2',
               written=nc.now(), second_placement=args.second, mcpb_dir=str(src),
               inputs_sha256={f: nc.sha(work / f) for f in ('resp1_equiv.in', 'resp2.in', 'CA2_large_mk.esp', 'second.esp',
                                                            'two1.in', 'two2.in')},
               resp_equiv_py_sha256=nc.sha(nc.HERE / 'resp_equiv.py'), geometry=geo, layout=dict(stage1=lay1, stage2=lay2),
               fit=dict(stage1=s1, stage2=s2), checks={})
    rec['checks']['the_two_molecules_have_the_same_charges'] = dict(largest_difference=float(np.abs(q - qb).max()),
                                                                    passed=bool(np.abs(q - qb).max() < 2e-6))
    idx = {name: k for k, (resi, resn, name) in enumerate(order) if resn == rq.LIG}
    worst = max(abs(q[idx[a]] - q[idx[b]]) for a, b in rq.PAIRS)
    rec['checks']['ligand_pairs_equal'] = dict(largest_difference=float(worst), passed=bool(worst < 5e-7))
    rec['checks']['total_charge'] = dict(value=round(float(q.sum()), 6), one_placement=round(float(q_one.sum()), 6),
                                         passed=bool(abs(q.sum() - q_one.sum()) < 2e-5))
    std = rq.pdb_atoms(work / 'CA2_standard.pdb')
    ix2 = {(resi, name): k for k, (resi, resn, name) in enumerate(order)}
    site = [ix2[(resi, name)] for resi, resn, name in std]
    rec['checks']['site_charge'] = dict(value=round(float(q[site].sum()), 6), one_placement=round(float(q_one[site].sum()), 6),
                                        passed=bool(abs(q[site].sum() - q_one[site].sum()) < 2e-5))
    # atoms that MCPB.py holds at their ff14SB charges or in neutral groups keep them
    if args.second == 'X':
        d = float(np.abs(q - q_one).max())
        rec['checks']['same_potential_twice_gives_the_one_placement_charges'] = dict(largest_difference=d, passed=bool(d < 2e-6))
    zn = [k for k, (resi, resn, name) in enumerate(order) if resn == 'ZN'][0]
    rec['zinc'] = dict(two_placements=round(float(q[zn]), 6), one_placement=round(float(q_one[zn]), 6))
    rec['ligand'] = [dict(atom=order[k][2], two_placements=round(float(q[k]), 6), one_placement=round(float(q_one[k]), 6),
                          change=round(float(q[k] - q_one[k]), 6)) for k in lig]
    rec['ligand_total'] = dict(two_placements=round(float(q[lig].sum()), 6), one_placement=round(float(q_one[lig].sum()), 6))
    notlig = [k for k in range(101) if k not in lig]
    rec['largest_change'] = dict(ligand=round(float(np.abs(q[lig] - q_one[lig]).max()), 6),
                                 other_atoms=round(float(np.abs(q[notlig] - q_one[notlig]).max()), 6),
                                 at=f'{order[int(np.argmax(np.abs(q - q_one)))]}', rms_all=round(float(np.sqrt(((q - q_one) ** 2).mean())), 6))
    # how well each set gives each potential (own evaluation of sum q / r on the fit points)
    rr = {}
    sets = dict(one_placement_X=q_one, two_placements=q)
    if args.second == 'Rx':
        # for comparison only: the same two-stage fit on the potential of the R-like placement alone
        shutil.copyfile(work / 'second.esp', work / 'only_second.esp')
        rq.run_resp(work, 'resp1_equiv.in', 'only1', 'only_second.esp')
        rq.run_resp(work, 'resp2.in', 'only2', 'only_second.esp', qin='only1.chg')
        q_r = rq.read_chg(work / 'only2.chg')
        sets['one_placement_Rx'] = q_r
        rec['ligand_one_placement_Rx'] = {order[k][2]: round(float(q_r[k]), 6) for k in lig}
        rec['zinc']['one_placement_Rx'] = round(float(q_r[zn]), 6)
        shutil.copyfile(work / 'only2.chg', out / 'only_Rx_stage2.chg')
    for name, qq in sets.items():
        rr[name] = dict(potential_of_X=round(potential_rrms(qq, e1), 5), potential_of_the_second=round(potential_rrms(qq, e2), 5))
    rec['relative_rms_of_the_potential'] = rr
    rec['all_charges'] = dict(order=[f'{resi}-{resn}-{name}' for resi, resn, name in order],
                              two_placements=[round(float(v), 6) for v in q], one_placement=[round(float(v), 6) for v in q_one])
    # the two potential jobs (the large model in the two placements), from their status files
    jobs = {}
    for label, f in (('X', pathlib.Path('/root/cqc-sfa-md/qm/prod/status.json')), ('Rx', OUT / 'status.json')):
        if f.exists():
            e = json.loads(f.read_text())['jobs']['esp']
            jobs[label] = dict(status_file=str(f), e_tot_hartree=e['e_tot'], n_fit_points=e['n_points'], finished=e.get('finished'),
                               seconds_scf=[e['scf_first']['seconds'], e.get('scf_final', {}).get('seconds')],
                               unrestrained_fit=e.get('unrestrained_fit'), esp_min=e.get('esp_min'), esp_max=e.get('esp_max'),
                               check_of_the_potential_against_one_electron_integrals=e.get('check_rinv_max_abs_difference'))
    if args.second == 'Rx' and len(jobs) == 2:
        jobs['Rx_minus_X_kcal_B3LYP_without_D3'] = round((jobs['Rx']['e_tot_hartree'] - jobs['X']['e_tot_hartree']) * nc.HARTREE, 3)
        rec['potential_jobs_of_the_large_model'] = jobs
    rec['all_checks_passed'] = bool(all(c['passed'] for c in rec['checks'].values()))
    for f in ('two1.in', 'two2.in', 'two1.chg', 'two2.chg', 'two1.out', 'two2.out'):
        shutil.copyfile(work / f, out / f)
    nc.write_json(out / 'resp2.json', rec)
    print(json.dumps({k: rec[k] for k in ('second_placement', 'fit', 'checks', 'zinc', 'ligand_total', 'largest_change',
                                          'relative_rms_of_the_potential')}, indent=1))
    for r in rec['ligand']:
        print(f"  {r['atom']:5s} {r['one_placement']:9.4f} {r['two_placements']:9.4f} {r['change']:+8.4f}")
    return 0 if rec['all_checks_passed'] else 1


# ------------------------------------------------------------------- mol2 files and the small model with them
# read_mol2, write_mol2 and capped_residue are VERBATIM copies from md/qm/gr_mm.py (working copy of 6 Oct 2026,
# 3:00 PM IST, sha256 245fca5c...; the three functions are the same text in git HEAD): the capped histidines of
# the small model, cap scheme "methyl", exactly as gate G-R made them.
def read_mol2(path):
    atoms, bonds, sec = [], [], None
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith('@<TRIPOS>'):
            sec = l.strip()[9:]
            continue
        f = l.split()
        if sec == 'ATOM' and len(f) >= 9:
            atoms.append(dict(name=f[1], xyz=[float(f[2]), float(f[3]), float(f[4])], type=f[5], q=float(f[8])))
        elif sec == 'BOND' and len(f) >= 3:
            bonds.append((int(f[1]), int(f[2])))
    return atoms, bonds


def write_mol2(path, resname, atoms, bonds):
    o = ['@<TRIPOS>MOLECULE', resname, f'{len(atoms):5d}{len(bonds):6d}     1     0     0', 'SMALL',
         'capped for the small model, md/qm/gr_mm.py', ' ', ' ', '@<TRIPOS>ATOM']
    for k, a in enumerate(atoms, 1):
        o.append(f"{k:7d} {a['name']:<4s}    {a['xyz'][0]:10.4f}{a['xyz'][1]:10.4f}{a['xyz'][2]:10.4f} "
                 f"{a['type']:<4s} {1:6d} {resname:<4s} {a['q']:12.6f}")
    o.append('@<TRIPOS>BOND')
    for k, (i, j) in enumerate(bonds, 1):
        o.append(f'{k:6d}{i:5d}{j:5d} 1')
    o += ['@<TRIPOS>SUBSTRUCTURE', f'     1 {resname}         1 TEMP              0 ****  ****    0 ROOT']
    pathlib.Path(path).write_text('\n'.join(o) + '\n')


def capped_residue(resname, src_atoms, model_names, model_xyz, scheme):
    """The 15 atoms of the small model for one histidine, from the 17 atoms of its mol2 file."""
    s = {a['name']: a for a in src_atoms}
    removed = s['N']['q'] + s['H']['q'] + s['C']['q'] + s['O']['q']
    out = []
    for nm, xyz in zip(model_names, model_xyz):
        if nm == 'H1':
            a = dict(name='H1', type='H1', q=s['N']['q'] + s['H']['q'])
        elif nm == 'H3':
            a = dict(name='H3', type='H1', q=s['C']['q'] + s['O']['q'])
        elif nm == 'CH3':
            a = dict(name='CH3', type=s['CA']['type'], q=s['CA']['q'])
        elif nm == 'H2':
            a = dict(name='H2', type=s['HA']['type'], q=s['HA']['q'])
        else:
            a = dict(name=nm, type=s[nm]['type'], q=s[nm]['q'])
        a['xyz'] = [float(t) for t in xyz]
        out.append(a)
    d = {a['name']: a for a in out}
    if scheme == 'methyl':
        tot = removed + s['CA']['q'] + s['HA']['q']
        for h in ('H1', 'H2', 'H3'):
            d[h]['q'] = s['HA']['q']
        d['CH3']['q'] = tot - 3 * s['HA']['q']
    elif scheme == 'spread':
        d['H1']['q'] = d['H3']['q'] = 0.0
        for a in out:
            a['q'] += removed / len(out)
    elif scheme != 'group':
        raise ValueError(scheme)
    total_src = sum(a['q'] for a in src_atoms)
    assert abs(sum(a['q'] for a in out) - total_src) < 1e-9, 'the capped residue does not keep its charge'
    # bonds of the template: those of the source among kept atoms, plus the two caps on CA
    keep = {a['name'] for a in out}
    ren = {'CA': 'CH3', 'HA': 'H2'}
    return out, removed, ren, keep


MARK2 = 'RESP Charge from the potentials of two placements, ligand pairs equal (md/qm/n2_charges.py)'
RES = ('HD1', 'HD2', 'HE1')


def cmd_mol2(args):
    """The five mol2 files of the zinc site with the charges of a resp2 run, in a folder that md/build_topology.py
    can take as --mcpb-dir (a copy of the files of the gate's MCPB.py folder that the builder reads)."""
    import resp_equiv as rq
    src = pathlib.Path(args.mcpb_dir)
    tag = 'two_placements' if args.second == 'Rx' else 'check_same_potential_twice'
    rec2 = json.loads((nc.N2 / 'charges2' / tag / 'resp2.json').read_text())
    assert rec2['all_checks_passed'], 'the resp2 run did not pass its checks'
    q = np.array(rec2['all_charges']['two_placements'])
    dst = nc.N2 / 'charges2' / ('mcpb_2e_two' if args.second == 'Rx' else 'mcpb_2e_check')
    if dst.exists():
        shutil.rmtree(dst)
    (dst / 'charges_unconstrained').mkdir(parents=True)
    keep = ['CA2_tleap.in', 'CA2_mcpbpy.pdb', 'CA2_mcpbpy.frcmod', 'sfa_anion_gaff2.frcmod', 'sfa_anion_repair.frcmod',
            'CA2_large.pdb', 'CA2_standard.pdb'] + [f'{m}.mol2' for m in rq.MOL2]
    for f in keep:
        shutil.copyfile(src / f, dst / f)
    for m in rq.MOL2:
        shutil.copyfile(src / 'charges_unconstrained' / f'{m}.mol2', dst / 'charges_unconstrained' / f'{m}.mol2')
    order, std = rq.pdb_atoms(dst / 'CA2_large.pdb'), rq.pdb_atoms(dst / 'CA2_standard.pdb')
    qf = rq.read_chg(src / 'resp2_free.chg')
    rec = dict(checks={})
    rq.write_mol2_sets(dst, order, std, qf, q, rec)
    one, new, table = {}, {}, {}
    for m in rq.MOL2:
        t = (dst / f'{m}.mol2').read_text()
        assert rq.MARK in t
        (dst / f'{m}.mol2').write_text(t.replace(rq.MARK, MARK2))
        n1, q1 = rq.mol2_charges(src / f'{m}.mol2')
        n2, q2 = rq.mol2_charges(dst / f'{m}.mol2')
        assert n1 == n2
        one[m], new[m] = q1, q2
        table[m] = dict(total_one_placement=round(float(q1.sum()), 6), total_two_placements=round(float(q2.sum()), 6),
                        largest_change=round(float(np.abs(q2 - q1).max()), 6),
                        atoms={n: [round(float(a), 6), round(float(b), 6)] for n, a, b in zip(n1, q1, q2)})
    site = sum(float(v.sum()) for v in new.values())
    out = dict(what='mol2 files of the five residues of the zinc site with the charges fitted to two placements; a folder '
                    'for md/build_topology.py --mcpb-dir; written by md/qm/n2_charges.py mol2', written=nc.now(),
               second_placement=args.second, source=str(src), folder=str(dst), residues=table,
               site_charge=round(site, 6), mapping_check=rec['checks'],
               files_sha256={f.name: nc.sha(f) for f in sorted(dst.iterdir()) if f.is_file()},
               unchanged_files_same_as_source={f: bool(nc.sha(dst / f) == nc.sha(src / f)) for f in keep[:7]},
               passed=bool(all(c['passed'] for c in rec['checks'].values()) and abs(site - 1.0) < 5e-4))
    nc.write_json(dst / 'mol2_set.json', out)
    print(json.dumps({k: out[k] for k in ('folder', 'site_charge', 'mapping_check', 'passed')}, indent=1))
    print({m: (v['total_one_placement'], v['total_two_placements'], v['largest_change']) for m, v in table.items()})
    return 0 if out['passed'] else 1


def cmd_ff(args):
    """The force-field copy of the small model with the new charges: md/qm/n2/<name>/inputs (the tleap inputs of gate
    G-R with the five mol2 files replaced, cap scheme "methyl"), small_as_built.prmtop, and mm_points.json: total
    and rest of every point (old and new) from n2_common.FF on that topology."""
    import parmed
    name = args.name
    assert name != 'ff'
    mdir = nc.N2 / 'charges2' / ('mcpb_2e_two' if args.second == 'Rx' else 'mcpb_2e_check')
    assert json.loads((mdir / 'mol2_set.json').read_text())['passed']
    base = nc.N2 / 'ff'
    dst = nc.N2 / name
    if dst.exists():
        shutil.rmtree(dst)
    (dst / 'inputs').mkdir(parents=True)
    for f in (base / 'inputs').iterdir():
        shutil.copyfile(f, dst / 'inputs' / f.name)
    # the capped histidines, exactly as gr_mm.build makes them: names and places from geom/X.pdb and X.xyz
    lines = (nc.GR / 'geom' / 'X.pdb').read_text().splitlines()
    model = [(l[12:16].strip(), l[17:20].strip(), int(l[22:26])) for l in lines if l.startswith(('ATOM', 'HETATM'))]
    _, xyz = nc.read_xyz(nc.GR / 'geom' / 'X.xyz')
    caps = {}
    for k, rn in enumerate(RES, 1):
        src_atoms, src_bonds = read_mol2(mdir / f'{rn}.mol2')
        idx = [i for i, m in enumerate(model) if m[2] == k]
        names = [model[i][0] for i in idx]
        atoms, removed, ren, keep = capped_residue(rn, src_atoms, names, xyz[idx], 'methyl')
        pos = {a['name']: n for n, a in enumerate(atoms, 1)}
        sname = [a['name'] for a in src_atoms]
        bonds = []
        for i, j in src_bonds:
            a, b = sname[i - 1], sname[j - 1]
            a, b = ren.get(a, a), ren.get(b, b)
            if a in pos and b in pos:
                bonds.append((pos[a], pos[b]))
        bonds += [(pos['CH3'], pos['H1']), (pos['CH3'], pos['H3'])]
        assert len(bonds) == 15, (rn, len(bonds))
        write_mol2(dst / 'inputs' / f'{rn}.mol2', rn, atoms, bonds)
        caps[rn] = dict(removed_backbone_charge=round(removed, 6),
                        **{a['name']: round(a['q'], 6) for a in atoms if a['name'] in ('H1', 'CH3', 'H2', 'H3')},
                        residue_total=round(sum(a['q'] for a in atoms), 6))
    for f in ('ZN1.mol2', 'SN1.mol2'):
        shutil.copyfile(mdir / f, dst / 'inputs' / f)
    # the same route with the files of gate G-R must give the files of gate G-R (proof that the copy is faithful)
    env_ff = os.environ.get('N2_FF')
    work = pathlib.Path(args.work) / name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for f in (dst / 'inputs').iterdir():
        shutil.copyfile(f, work / f.name)
    env = dict(os.environ, AMBERHOME=nc.AMBERHOME, PATH=f'{nc.AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    p = subprocess.run(['tleap', '-s', '-f', 'leap.in'], cwd=work, env=env, capture_output=True, text=True)
    (work / 'leap.log').write_text(p.stdout + p.stderr)
    bad = [l for l in (p.stdout + p.stderr).splitlines() if 'Could not find' in l or 'FATAL' in l]
    if not (work / 'small.prmtop').exists() or bad:
        raise SystemExit(f'tleap failed: {bad[:5]}; see {work}/leap.log')
    shutil.copyfile(work / 'small.prmtop', dst / 'small_as_built.prmtop')
    # what differs from the copy as built: charges only
    a, b = parmed.load_file(str(base / 'small_as_built.prmtop')), parmed.load_file(str(dst / 'small_as_built.prmtop'))

    def strip(parm):
        nm = lambda x: f'{x.residue.name}:{x.name}'                       # noqa: E731
        return dict(atoms=[(nm(x), x.type, round(x.mass, 4), round(x.rmin, 5), round(x.epsilon, 6)) for x in parm.atoms],
                    bonds=sorted((nm(x.atom1), nm(x.atom2), round(x.type.k, 4), round(x.type.req, 5)) for x in parm.bonds),
                    angles=sorted((nm(x.atom1), nm(x.atom2), nm(x.atom3), round(x.type.k, 4), round(x.type.theteq, 4)) for x in parm.angles),
                    dihedrals=sorted((nm(x.atom1), nm(x.atom2), nm(x.atom3), nm(x.atom4), int(x.type.per), round(x.type.phi_k, 5),
                                      round(x.type.phase, 3), bool(x.improper), bool(x.ignore_end), round(x.type.scee, 4),
                                      round(x.type.scnb, 4)) for x in parm.dihedrals))
    sa, sb = strip(a), strip(b)
    dq = np.array([y.charge - x.charge for x, y in zip(a.atoms, b.atoms)])
    same_but_charges = {k: bool(sa[k] == sb[k]) for k in sa}
    # energies of every point
    ff = nc.FF(dst / 'small_as_built.prmtop')
    g = nc.Geometry()
    pj = json.loads((nc.N2 / 'points.json').read_text())
    pts = {}
    for pid in pj['old_rigid_points']:
        e = ff.energy(g.full(g.old[pid]))
        pts[pid] = dict(total=round(e['total'], 6), rest=round(e['rest'], 6))
    for q in pj['points']:
        _, x = nc.read_xyz(nc.N2 / 'geom' / f"{q['id']}.xyz")
        e = ff.energy(x)
        pts[q['id']] = dict(total=round(e['total'], 6), rest=round(e['rest'], 6))
    old = {k: v['mm_as_built'] for k, v in pj['old_rigid_points'].items()}
    old.update({q['id']: q['mm_as_built'] for q in pj['points']})
    d = np.array([(pts[k]['total'] - pts['X']['total']) - (old[k]['total'] - old['X']['total']) for k in pts])
    what = (f'small model with the charges fitted to the potentials of two placements (X and {args.second}), cap scheme '
            f'"methyl"; every other parameter as built; written by md/qm/n2_charges.py ff')
    nc.write_json(dst / 'mm_points.json', dict(
        what=what, written=nc.now(), prmtop_sha256=nc.sha(dst / 'small_as_built.prmtop'), points=pts,
        removed_for_rest=[[round(v, 5) if isinstance(v, float) else v for v in r] for r in ff.removed['rest']],
        change_of_the_energy_relative_to_X_against_the_copy_as_built=dict(
            rms=round(float(np.sqrt((d ** 2).mean())), 4), largest=round(float(d[np.argmax(np.abs(d))]), 4),
            Rx=round(float((pts['Rx']['total'] - pts['X']['total']) - (old['Rx']['total'] - old['X']['total'])), 4))))
    sums = ['  '.join((nc.sha(f), str(f.relative_to(dst)))) for f in [dst / 'small_as_built.prmtop'] + sorted((dst / 'inputs').iterdir())]
    (dst / 'SHA256SUMS').write_text('\n'.join(sums) + '\n')
    rec = dict(what=what, written=nc.now(), folder=str(dst), mol2_from=str(mdir), cap_charges=caps,
               everything_but_the_charges_equal_to_the_copy_as_built=same_but_charges,
               charge=dict(total=round(float(sum(x.charge for x in b.atoms)), 6), largest_change=round(float(np.abs(dq).max()), 6),
                           atoms_changed=int((np.abs(dq) > 1e-7).sum())),
               inputs_identical_to_ff={f.name: bool(nc.sha(f) == nc.sha(base / 'inputs' / f.name)) for f in sorted((dst / 'inputs').iterdir())},
               prmtop_identical_to_ff=bool(nc.sha(dst / 'small_as_built.prmtop') == nc.sha(base / 'small_as_built.prmtop')),
               energy_change_relative_to_X=dict(rms=round(float(np.sqrt((d ** 2).mean())), 4),
                                                largest=round(float(d[np.argmax(np.abs(d))]), 4), n_points=len(d)),
               gr_mm_py_copied_from_sha256='245fca5cf5245592f8a490dd1f3f05628b82383862d717c911a91f80dde68a07',
               passed=bool(all(same_but_charges.values()) and abs(sum(x.charge for x in b.atoms) - 1.0) < 5e-4))
    nc.write_json(dst / 'ff_copy.json', rec)
    print(json.dumps({k: rec[k] for k in ('folder', 'everything_but_the_charges_equal_to_the_copy_as_built', 'charge',
                                          'prmtop_identical_to_ff', 'energy_change_relative_to_X', 'passed')}, indent=1))
    del env_ff
    return 0 if rec['passed'] else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prepare-esp', 'resp2', 'mol2', 'ff'])
    ap.add_argument('--second', choices=['Rx', 'X'], default='Rx',
                    help='the second potential: Rx (the step itself) or X (the same potential twice: the check)')
    ap.add_argument('--mcpb-dir', default='/root/cqc-sfa-md/qm/gr/mcpb_2e')
    ap.add_argument('--work', default=None, help='scratch folder (resp2, ff)')
    ap.add_argument('--name', default='ff2', help='ff: name of the force-field copy under md/qm/n2/')
    args = ap.parse_args()
    if args.cmd in ('resp2', 'ff') and not args.work:
        ap.error('--work SCRATCH is needed')
    return dict([('prepare-esp', cmd_prepare_esp), ('resp2', cmd_resp2), ('mol2', cmd_mol2), ('ff', cmd_ff)])[args.cmd](args)


if __name__ == '__main__':
    sys.exit(main())
