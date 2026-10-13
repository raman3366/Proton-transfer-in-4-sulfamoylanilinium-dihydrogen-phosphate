"""RESP charges of the zinc site with the chemically equal ligand atoms made equal (md/PLAN.md, line of
6 Oct 2026, 3:02 AM IST, point 5).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/resp_equiv.py WORK [REPORT.json]
    ... md/qm/resp_equiv.py --selftest SCRATCH [REPORT.json]      (synthetic potentials; no quantum chemistry)

WORK is the MCPB.py folder after step 3 (md/qm/after_qm.sh calls this script right after that step).

WHY.  The potential of the large model is that of pose X.  MCPB.py puts no equivalence on ligand atoms (its
second resp stage only equalises hydrogens on CH2 and CH3 carbons), so O1 and O2, the ring pairs and the two
aniline hydrogens would get the unequal charges of that one pose and favour it.

WHAT IS DONE.  MCPB.py's two resp inputs are kept as they are, with one change: in the FIRST stage
(resp1.in, all atoms free, restraint weight 0.0005) the second atom of each pair is declared equivalent to
the first (the column `ivary` of resp):
      O2 = O1    C6 = C2    C5 = C3    H6 = H2    H5 = H3    HN12 = HN11
This is the stage in which the RESP procedure imposes equivalences of atoms other than methyl and
methylene hydrogens, so that all other charges adapt to them; the second stage (resp2.in, restraint weight
0.001) refits only the CH2 and CH3 groups of the protein residues and keeps every ligand atom at its
first-stage value, so the pairs stay equal.  All constraints MCPB.py wrote (total charge, neutral caps,
backbone atoms at their ff14SB charges) stay in force, hence the total charge is kept by the fit itself.
Averaging the pairs afterwards would also keep the total, but the other charges would not adapt; the averaged
set is computed and reported beside for comparison only.

FILES (in WORK)
  resp1_equiv.in, resp1_equiv.out/.chg, resp2_equiv.out/.chg      the constrained fit            (USED)
  resp1_free.out/.chg, resp2_free.out/.chg                        MCPB.py's own inputs, run again (BESIDE)
  HD1.mol2 HD2.mol2 HE1.mol2 ZN1.mol2 SN1.mol2                    rewritten with the constrained charges; the
                                                                  comment line says so
  charges_unconstrained/<the five mol2 files>                     the files MCPB.py step 3 wrote, copied
                                                                  unchanged before they are rewritten
  resp_equiv.json                                                 both sets, differences, checks
A second run gives the same files: everything is recomputed from resp1.in, resp2.in and the .esp file; the
mol2 files are rewritten from the atoms and bonds they hold, only the charge column and the comment change;
the copies in charges_unconstrained/ are renewed only when MCPB.py step 3 has written new files.

CHECKS made on every run (the script stops when one fails)
  - the resp order is the atom order of CA2_large.pdb; the six pairs are found by name in residue SAN;
  - the unconstrained charges computed here are those of MCPB.py's resp2.chg (when that file is there);
  - the charges MCPB.py wrote into its mol2 files are reproduced by the mapping used here;
  - in the constrained set the six pairs are equal, the five residues of the site sum to the same charge as
    in the unconstrained set, and resp ended normally.
"""
import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
GROUP = 'CA2'
PAIRS = (('O1', 'O2'), ('C2', 'C6'), ('C3', 'C5'), ('H2', 'H6'), ('H3', 'H5'), ('HN11', 'HN12'))
LIG = 'SAN'
MOL2 = ('HD1', 'HD2', 'HE1', 'ZN1', 'SN1')
MARK = 'RESP Charge, ligand pairs equal (md/qm/resp_equiv.py)'


def pdb_atoms(path):
    out = []
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            out.append((int(l[22:26]), l[17:20].strip(), l[12:16].strip()))
    return out


def read_chg(path):
    return np.array([float(t) for t in pathlib.Path(path).read_text().split()])


def split_resp_in(text):
    """Head (up to the line 'charge natoms'), the atom lines, the rest."""
    lines = text.splitlines()
    k = [i for i, l in enumerate(lines) if l.strip() == '&end'][0] + 3        # weight, title, charge/natoms
    natm = int(lines[k].split()[1])
    return lines[:k + 1], lines[k + 1:k + 1 + natm], lines[k + 1 + natm:]


def run_resp(work, inp, out_tag, esp, qin=None):
    cmd = [f'{AMBERHOME}/bin/resp', '-O', '-i', inp, '-o', f'{out_tag}.out', '-p', f'{out_tag}.pch', '-t',
           f'{out_tag}.chg', '-e', esp, '-s', f'{out_tag}_calc.esp']
    if qin:
        cmd += ['-q', qin]
    env = dict(os.environ, AMBERHOME=AMBERHOME)
    p = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True)
    out = (pathlib.Path(work) / f'{out_tag}.out')
    if p.returncode != 0 or not out.exists() or not (pathlib.Path(work) / f'{out_tag}.chg').exists():
        raise SystemExit(f'resp failed on {inp}: {p.stdout[-300:]} {p.stderr[-300:]}')
    stat = {}
    for l in out.read_text().splitlines():
        if 'ESP relative RMS' in l:
            stat['relative_rms'] = float(l.split()[-1])
        if 'std err of estimate' in l:
            stat['std_err'] = float(l.split()[-1])
    if 'relative_rms' not in stat:
        raise SystemExit(f'resp wrote no statistics for {inp}')
    for f in (f'{out_tag}_calc.esp', f'{out_tag}.pch'):
        (pathlib.Path(work) / f).unlink(missing_ok=True)
    return stat


def equivalence_input(work, order, qwt=None):
    """resp1_equiv.in from resp1.in; returns the pairs as indices (from 1)."""
    head, atoms, rest = split_resp_in((work / 'resp1.in').read_text())
    assert len(atoms) == len(order), 'resp1.in and CA2_large.pdb differ in the number of atoms'
    lig = {name: k for k, (resi, resn, name) in enumerate(order, 1) if resn == LIG}
    pairs = []
    for a, b in PAIRS:
        if a not in lig or b not in lig:
            raise SystemExit(f'ligand atoms {a}, {b} not found in the large model')
        i, j = lig[a], lig[b]
        assert i < j, 'the first atom of a pair must come first in the resp input'
        z, iv = atoms[j - 1].split()
        zi = atoms[i - 1].split()[0]
        assert z == zi and iv == '0', f'pair {a}/{b}: unexpected line {atoms[j - 1]!r}'
        atoms[j - 1] = '%5d %4d' % (int(z), i)
        pairs.append((a, b, i, j))
    text = '\n'.join(head + atoms + rest) + '\n'
    if qwt is not None:
        text = text.replace('qwt = 0.00050,', f'qwt = {qwt:.5f},')
    (work / 'resp1_equiv.in').write_text(text)
    return pairs


def zero_weight(work, name):
    t = (work / name).read_text()
    assert 'qwt = 0.00050,' in t or 'qwt = 0.00100,' in t
    t = t.replace('qwt = 0.00050,', 'qwt = 0.00000,').replace('qwt = 0.00100,', 'qwt = 0.00000,')
    q0 = name.replace('.in', '_q0.in')
    (work / q0).write_text(t)
    return q0


def rewrite_mol2(src, dst, charges, mark):
    """The mol2 file with other charges: only the charge column and the comment line change."""
    lines = pathlib.Path(src).read_text().splitlines()
    out, sec, k = [], None, 0
    for n, l in enumerate(lines):
        if l.startswith('@<TRIPOS>'):
            sec = l.strip()[9:]
            out.append(l)
            continue
        if sec == 'MOLECULE' and l.strip().startswith('RESP Charge'):
            out.append(mark)
            continue
        if sec == 'ATOM' and len(l.split()) >= 9:
            f = l.split()
            out.append('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f' % (
                int(f[0]), f[1], float(f[2]), float(f[3]), float(f[4]), f[5], int(f[6]), f[7], charges[k]))
            k += 1
            continue
        out.append(l)
    assert k == len(charges), (src, k, len(charges))
    pathlib.Path(dst).write_text('\n'.join(out) + '\n')


def mol2_charges(path):
    q, names, sec = [], [], None
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith('@<TRIPOS>'):
            sec = l.strip()[9:]
            continue
        f = l.split()
        if sec == 'ATOM' and len(f) >= 9:
            names.append(f[1])
            q.append(float(f[8]))
    return names, np.array(q)


def process(work, qwt_zero=False, tag=''):
    """Both fits in `work`; returns the record.  qwt_zero: restraint weight 0 in both stages (tests)."""
    work = pathlib.Path(work)
    esp = f'{GROUP}_large_mk.esp'
    for f in ('resp1.in', 'resp2.in', esp, f'{GROUP}_large.pdb', f'{GROUP}_standard.pdb'):
        if not (work / f).exists():
            raise SystemExit(f'{work}/{f} missing: run MCPB.py step 3 first')
    order = pdb_atoms(work / f'{GROUP}_large.pdb')
    std = pdb_atoms(work / f'{GROUP}_standard.pdb')
    pairs = equivalence_input(work, order)
    in1, in1e, in2 = 'resp1.in', 'resp1_equiv.in', 'resp2.in'
    if qwt_zero:
        in1, in1e, in2 = zero_weight(work, in1), zero_weight(work, in1e), zero_weight(work, in2)
    s1f = run_resp(work, in1, f'resp1_free{tag}', esp)
    s2f = run_resp(work, in2, f'resp2_free{tag}', esp, qin=f'resp1_free{tag}.chg')
    s1e = run_resp(work, in1e, f'resp1_equiv{tag}', esp)
    s2e = run_resp(work, in2, f'resp2_equiv{tag}', esp, qin=f'resp1_equiv{tag}.chg')
    qf, qe = read_chg(work / f'resp2_free{tag}.chg'), read_chg(work / f'resp2_equiv{tag}.chg')
    assert len(qf) == len(qe) == len(order)
    rec = dict(work=str(work), restraint_weight_zero=bool(qwt_zero), n_atoms_large_model=len(order),
               pairs=[dict(first=a, second=b, index_first=i, index_second=j) for a, b, i, j in pairs],
               fit=dict(free=dict(stage1=s1f, stage2=s2f), equal=dict(stage1=s1e, stage2=s2e)), checks={})
    # --- MCPB.py's own result, when it is there
    own = work / 'resp2.chg'
    if own.exists() and not qwt_zero:
        d = float(np.abs(read_chg(own) - qf).max())
        rec['checks']['free_fit_equals_resp2_chg_of_mcpb'] = dict(largest_difference=d, passed=bool(d < 5e-6))
    # --- pairs, totals
    idx = {(resi, name): k for k, (resi, resn, name) in enumerate(order)}
    worst_pair = max(abs(qe[i - 1] - qe[j - 1]) for _, _, i, j in pairs)
    rec['checks']['pairs_equal_in_constrained_set'] = dict(largest_difference=float(worst_pair),
                                                           passed=bool(worst_pair < 5e-7))
    site = [idx[(resi, name)] for resi, resn, name in std]
    tot_f, tot_e = float(qf[site].sum()), float(qe[site].sum())
    rec['checks']['site_charge_kept'] = dict(free=round(tot_f, 6), equal=round(tot_e, 6),
                                             whole_model_free=round(float(qf.sum()), 6),
                                             whole_model_equal=round(float(qe.sum()), 6),
                                             passed=bool(abs(tot_f - tot_e) < 2e-5 and abs(qf.sum() - qe.sum()) < 2e-5))
    # --- the averaged set, for comparison only
    qa = qf.copy()
    for _, _, i, j in pairs:
        qa[i - 1] = qa[j - 1] = 0.5 * (qf[i - 1] + qf[j - 1])
    lig = [k for k, (resi, resn, name) in enumerate(order) if resn == LIG]
    not_lig = [k for k in range(len(order)) if k not in lig]
    rec['ligand'] = [dict(atom=order[k][2], equal=round(float(qe[k]), 6), free=round(float(qf[k]), 6),
                          averaged=round(float(qa[k]), 6)) for k in lig]
    rec['ligand_total'] = dict(equal=round(float(qe[lig].sum()), 6), free=round(float(qf[lig].sum()), 6))
    zn = [k for k, (resi, resn, name) in enumerate(order) if resn == 'ZN'][0]
    rec['zinc'] = dict(equal=round(float(qe[zn]), 6), free=round(float(qf[zn]), 6))
    rec['pair_difference_in_free_set'] = {f'{a}-{b}': round(float(qf[i - 1] - qf[j - 1]), 6) for a, b, i, j in pairs}
    rec['largest_change_equal_minus_free'] = dict(
        ligand=round(float(np.abs(qe[lig] - qf[lig]).max()), 6),
        other_atoms=round(float(np.abs(qe[not_lig] - qf[not_lig]).max()), 6))
    rec['largest_difference_equal_minus_averaged'] = round(float(np.abs(qe - qa).max()), 6)
    rec['all_charges'] = dict(order=[f'{resi}-{resn}-{name}' for resi, resn, name in order],
                              equal=[round(float(t), 6) for t in qe], free=[round(float(t), 6) for t in qf])
    return rec, order, std, qf, qe


def write_mol2_sets(work, order, std, qf, qe, rec):
    work = pathlib.Path(work)
    idx = {(resi, name): k for k, (resi, resn, name) in enumerate(order)}
    resis = []
    for resi, resn, name in std:
        if resi not in resis:
            resis.append(resi)
    assert len(resis) == len(MOL2)
    keep = work / 'charges_unconstrained'
    keep.mkdir(exist_ok=True)
    worst = 0.0
    for resi, m in zip(resis, MOL2):
        if MARK not in (work / f'{m}.mol2').read_text():
            shutil.copyfile(work / f'{m}.mol2', keep / f'{m}.mol2')          # MCPB.py's own file, new
        if not (keep / f'{m}.mol2').exists():
            raise SystemExit(f'{m}.mol2 carries the constrained charges but {keep}/{m}.mol2 is missing')
        names_std = [name for r, resn, name in std if r == resi]
        names, q_file = mol2_charges(keep / f'{m}.mol2')
        assert names == names_std, f'{m}.mol2 and the standard model differ in atom names'
        k = [idx[(resi, n)] for n in names]
        # the charges MCPB.py wrote must be the unconstrained ones through this mapping
        worst = max(worst, float(np.abs(q_file - qf[k]).max()))
        rewrite_mol2(keep / f'{m}.mol2', work / f'{m}.mol2', qe[k], MARK)
    rec['checks']['mapping_reproduces_mol2_of_mcpb'] = dict(largest_difference=worst, passed=bool(worst < 2e-6))
    names, q = mol2_charges(work / 'SN1.mol2')
    d = dict(zip(names, q))
    rec['checks']['pairs_equal_in_SN1_mol2'] = dict(
        passed=bool(all(abs(d[a] - d[b]) < 1e-9 for a, b in PAIRS)))


def main_work(work, report=None):
    rec, order, std, qf, qe = process(work)
    write_mol2_sets(work, order, std, qf, qe, rec)
    rec['all_checks_passed'] = bool(all(c['passed'] for c in rec['checks'].values()))
    rec['used'] = 'equal (the mol2 files in the work folder); free: MCPB.py\'s own files, in charges_unconstrained/'
    text = json.dumps(rec, indent=1) + '\n'
    (pathlib.Path(work) / 'resp_equiv.json').write_text(text)
    if report:
        pathlib.Path(report).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(report).write_text(text)
    print('resp_equiv: pairs equal %s, site charge %s -> %s, largest change ligand %.4f, other atoms %.4f e; '
          'relative RMS %.5f (free) %.5f (equal)' % (
              rec['checks']['pairs_equal_in_constrained_set']['passed'], rec['checks']['site_charge_kept']['free'],
              rec['checks']['site_charge_kept']['equal'], rec['largest_change_equal_minus_free']['ligand'],
              rec['largest_change_equal_minus_free']['other_atoms'], rec['fit']['free']['stage2']['relative_rms'],
              rec['fit']['equal']['stage2']['relative_rms']))
    for name, c in rec['checks'].items():
        print(('PASS ' if c['passed'] else 'FAIL ') + name)
    if not rec['all_checks_passed']:
        raise SystemExit('resp_equiv: a check failed; see resp_equiv.json')
    return rec


# ---------------------------------------------------------------------------------------------------- self-test
def read_esp(path):
    lines = pathlib.Path(path).read_text().splitlines()
    natm, npts = int(lines[0][0:5]), int(lines[0][5:10])
    at = np.array([[float(l[16 + 16 * k:32 + 16 * k]) for k in range(3)] for l in lines[1:1 + natm]])
    body = np.array([[float(l[16 * k:16 + 16 * k]) for k in range(4)] for l in lines[1 + natm:1 + natm + npts]])
    return at, body[:, 1:], body[:, 0]


def constraints_of(text):
    """Stage-1 constraints of a resp input: equivalences, group charges; total charge."""
    head, atoms, rest = split_resp_in(text)
    total = float(head[-1].split()[0])
    equiv = [(k, int(l.split()[1])) for k, l in enumerate(atoms, 1) if int(l.split()[1]) > 0]
    groups, k = [], 0
    rest = [l for l in rest if l.strip()]
    while k < len(rest):
        n, q = int(rest[k].split()[0]), float(rest[k].split()[1])
        ids = [int(t) for t in rest[k + 1].split()[1::2]]
        assert len(ids) == n
        groups.append((ids, q))
        k += 2
    return len(atoms), total, equiv, groups


def constrained_lsq(at, pts, v, total, equiv, groups):
    """Least-squares charges under linear constraints (no restraint): own code, for the comparison."""
    a = 1.0 / np.linalg.norm(pts[:, None, :] - at[None, :, :], axis=2)
    n = at.shape[0]
    rows, rhs = [np.ones(n)], [total]
    for j, i in equiv:
        r = np.zeros(n)
        r[j - 1], r[i - 1] = 1.0, -1.0
        rows.append(r)
        rhs.append(0.0)
    for ids, q in groups:
        r = np.zeros(n)
        r[[t - 1 for t in ids]] = 1.0
        rows.append(r)
        rhs.append(q)
    c, d = np.array(rows), np.array(rhs)
    m = len(d)
    kkt = np.block([[a.T @ a, c.T], [c, np.zeros((m, m))]])
    sol = np.linalg.solve(kkt, np.concatenate([a.T @ v, d]))
    return sol[:n]


def selftest(scratch, report=None):
    """Synthetic potentials of known charges on the real large model, through MCPB.py step 3 and this script."""
    sys.path.insert(0, str(HERE))
    import gaucom
    import gauout
    import test_after_qm as zt
    scratch = pathlib.Path(scratch)
    if scratch.exists():
        shutil.rmtree(scratch)
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    res = dict(what='self-test of md/qm/resp_equiv.py on synthetic potentials (no quantum chemistry)', cases={})

    def make(case, q):
        work = scratch / case
        shutil.copytree(HERE / 'mcpb', work)
        for f in list(work.glob('*.fchk')) + list(work.glob(f'{GROUP}_large_mk.log')) + list(work.glob('resp*')):
            f.unlink()
        com = gaucom.read_com(work / f'{GROUP}_large_mk.com')
        xyz = np.array(com.coords)
        pts, _, _ = gauout.mk_points(com.symbols, xyz, com.radii, density=float(com.keywords['iop']['6/42']))
        d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2) / gauout.BOHR
        gauout.write_mk_log(work / f'{GROUP}_large_mk.log', com.symbols, xyz, pts, np.zeros(len(xyz)),
                            (q[None, :] / d).sum(axis=1),
                            header_lines=[f'SYNTHETIC ({case}): potential of known point charges, test of resp_equiv.py'])
        p = subprocess.run([f'{AMBERHOME}/bin/python', f'{AMBERHOME}/bin/MCPB.py', '-i', f'{GROUP}.in', '-s', '3'],
                           cwd=work, env=env, capture_output=True, text=True)
        (work / 'step3.log').write_text(p.stdout + p.stderr)
        if p.returncode != 0 or not (work / 'SN1.mol2').exists():
            raise SystemExit(f'MCPB.py step 3 failed in {work}')
        return work, len(pts)

    work0 = scratch / '_names'
    shutil.copytree(HERE / 'mcpb', work0)
    q_known = zt.known_charges(work0)
    order = pdb_atoms(work0 / f'{GROUP}_large.pdb')
    lig = {name: k for k, (resi, resn, name) in enumerate(order) if resn == LIG}
    q_sym = q_known.copy()
    for a, b in PAIRS:
        q_sym[lig[a]] = q_sym[lig[b]] = 0.5 * (q_known[lig[a]] + q_known[lig[b]])
    shutil.rmtree(work0)
    res['known_charges'] = dict(
        sum=round(float(q_known.sum()), 6), sum_symmetrised=round(float(q_sym.sum()), 6),
        pair_differences={f'{a}-{b}': round(float(q_known[lig[a]] - q_known[lig[b]]), 6) for a, b in PAIRS})

    # case A: the truth is NOT symmetric (O1 -0.736, O2 -0.780 and so on, the AM1-BCC charges of the anion)
    work, npts = make('A_unsymmetric_truth', q_known)
    rec0, order_, std, qf0, qe0 = process(work, qwt_zero=True, tag='_q0')
    at, pts, v = read_esp(work / f'{GROUP}_large_mk.esp')
    n, total, equiv, groups = constraints_of((work / 'resp1_equiv.in').read_text())
    q_own = constrained_lsq(at, pts, v, total, equiv, groups)
    q1e = read_chg(work / 'resp1_equiv_q0.chg')
    n_, total_, equiv_f, groups_f = constraints_of((work / 'resp1.in').read_text())
    q_own_free = constrained_lsq(at, pts, v, total_, equiv_f, groups_f)
    q1f = read_chg(work / 'resp1_free_q0.chg')
    lig_idx = sorted(lig.values())
    a = dict(
        n_points=npts, n_equivalences_stage1_free=len(equiv_f), n_equivalences_stage1_equal=len(equiv),
        n_group_constraints=len(groups),
        free_without_restraint_returns_known=dict(largest_deviation=float(np.abs(qf0 - q_known).max())),
        equal_without_restraint=dict(
            pairs_largest_difference=rec0['checks']['pairs_equal_in_constrained_set']['largest_difference'],
            largest_deviation_from_pair_average_of_truth_ligand=float(np.abs(qe0[lig_idx] - q_sym[lig_idx]).max()),
            largest_deviation_from_truth_other_atoms=float(np.abs(np.delete(qe0 - q_known, lig_idx)).max()),
            site_charge=rec0['checks']['site_charge_kept']),
        stage1_against_own_constrained_least_squares=dict(
            equal_largest_difference=float(np.abs(q1e - q_own).max()),
            free_largest_difference=float(np.abs(q1f - q_own_free).max())),
        relative_rms=dict(free=rec0['fit']['free']['stage2']['relative_rms'],
                          equal=rec0['fit']['equal']['stage2']['relative_rms']))
    a['passed'] = bool(a['free_without_restraint_returns_known']['largest_deviation'] < 5e-3 and
                       a['equal_without_restraint']['pairs_largest_difference'] < 5e-7 and
                       a['stage1_against_own_constrained_least_squares']['equal_largest_difference'] < 2e-4 and
                       a['stage1_against_own_constrained_least_squares']['free_largest_difference'] < 2e-4 and
                       rec0['checks']['site_charge_kept']['passed'])
    # the production path (with the restraint) on the same folder: this is what runs unattended
    recp = main_work(work, report=None)
    again = (work / 'SN1.mol2').read_text(), (work / 'charges_unconstrained' / 'SN1.mol2').read_text(), \
        (work / 'resp_equiv.json').read_text()
    main_work(work, report=None)
    same = again == ((work / 'SN1.mol2').read_text(), (work / 'charges_unconstrained' / 'SN1.mol2').read_text(),
                     (work / 'resp_equiv.json').read_text())
    a['production_path_with_restraint'] = dict(
        all_checks_passed=recp['all_checks_passed'], checks=recp['checks'], ligand=recp['ligand'], zinc=recp['zinc'],
        largest_change_equal_minus_free=recp['largest_change_equal_minus_free'],
        largest_difference_equal_minus_averaged=recp['largest_difference_equal_minus_averaged'],
        relative_rms=dict(free=recp['fit']['free']['stage2']['relative_rms'],
                          equal=recp['fit']['equal']['stage2']['relative_rms']),
        second_run_gives_the_same_files=bool(same))
    a['passed'] = bool(a['passed'] and recp['all_checks_passed'] and same)
    res['cases']['A_unsymmetric_truth'] = a

    # case B: the truth IS symmetric; the constrained fit must give it back
    work, npts = make('B_symmetric_truth', q_sym)
    rec0, order_, std, qf0, qe0 = process(work, qwt_zero=True, tag='_q0')
    b = dict(n_points=npts,
             equal_without_restraint_returns_known=dict(largest_deviation=float(np.abs(qe0 - q_sym).max())),
             free_without_restraint_returns_known=dict(largest_deviation=float(np.abs(qf0 - q_sym).max())),
             relative_rms=dict(free=rec0['fit']['free']['stage2']['relative_rms'],
                               equal=rec0['fit']['equal']['stage2']['relative_rms']))
    b['passed'] = bool(b['equal_without_restraint_returns_known']['largest_deviation'] < 5e-3 and
                       b['free_without_restraint_returns_known']['largest_deviation'] < 5e-3)
    res['cases']['B_symmetric_truth'] = b
    res['all_passed'] = bool(a['passed'] and b['passed'])
    text = json.dumps(res, indent=1) + '\n'
    if report:
        pathlib.Path(report).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(report).write_text(text)
    print(json.dumps({k: (v if k != 'cases' else {c: w['passed'] for c, w in v.items()}) for k, v in res.items()
                      if k != 'known_charges'}, indent=1))
    print('SELF-TEST PASSED' if res['all_passed'] else 'SELF-TEST FAILED')
    return 0 if res['all_passed'] else 1


if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == '--selftest':
        sys.exit(selftest(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None))
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main_work(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
