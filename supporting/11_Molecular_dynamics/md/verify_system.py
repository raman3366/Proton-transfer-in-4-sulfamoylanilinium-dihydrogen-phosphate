"""Independent check of the files written by md/build_system.py (stage S3, part 1).

    /root/cqc-sfa-tools/dock/bin/python md/verify_system.py

It does not import the builder or its helper module.  From the written PDB and SDF files alone it

  1. recomputes the numbers of the three starts (own geometry code) and compares them with starts.csv;
  2. checks that the heavy atoms of X are the deposited SAN 304 (through dock/T3_6RL9/ref_SAN.sdf, a file
     of the docking part), those of N the seed-3 record of the exported top poses
     (dock/results/poses/T3_6RL9/sfa_neutral__ad4zn__top_pose_per_seed.sdf) and those of R MODEL 3 of
     seed 14 of the production job;
  3. recomputes which waters lie within 2.2 A of each ligand and compares with the water files;
  4. looks for heavy-atom contacts below 2.2 A in each complete start (protein, zinc, ligand, water);
  5. checks the ligand files (18 atoms, charge -1, names, SDF and PDB coordinates equal);
  6. runs build_system.py and system_audit.py again and compares the checksums of every output.

Result: md/system/verify.json; exit status 0 only if every check passes.
"""
import csv
import hashlib
import json
import pathlib
import subprocess
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / 'md' / 'system'
DOCK = REPO / 'dock'
TOL = 2e-3          # A
TOL_ANG = 0.06      # deg
ANGLES = ('N2_H_OG1', 'C3_C4_S_N2', 'C5_C4_S_N2', 'tau_fold', 'C4_S_N2', 'S_N2_Zn')
NAMES = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S', 'H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21']
DEPOSITED_ORDER = ['C1', 'N1', 'C2', 'C3', 'C4', 'C5', 'C6', 'S', 'N2', 'O1', 'O2']     # SAN 304 in the entry


def read_pdb(path):
    out = []
    for ln in pathlib.Path(path).read_text().splitlines():
        if ln.startswith(('ATOM', 'HETATM')):
            out.append(dict(name=ln[12:16].strip(), res=ln[17:20].strip(), chain=ln[21], num=int(ln[22:26]),
                            xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]), el=ln[76:78].strip()))
    return out


def read_sdf(path, record=0):
    recs = pathlib.Path(path).read_text().split('$$$$\n')
    lines = recs[record].split('\n')
    na, nb = int(lines[3][:3]), int(lines[3][3:6])
    atoms = [(ln.split()[3], np.array([float(v) for v in ln.split()[:3]])) for ln in lines[4:4 + na]]
    chg = {}
    for ln in lines:
        if ln.startswith('M  CHG'):
            p = ln.split()
            for i in range(int(p[2])):
                chg[int(p[3 + 2 * i])] = int(p[4 + 2 * i])
    props = {}
    for i, ln in enumerate(lines):
        if ln.startswith('>  <'):
            props[ln[4:ln.index('>', 4)]] = lines[i + 1]
    return atoms, chg, props, nb


def d(a, b):
    return float(np.sqrt(((a - b) ** 2).sum()))


def ang(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v)))


def tors(p0, p1, p2, p3):
    # IUPAC sign, from the two plane normals: atan2(|b2| b1.(b2 x b3), (b1 x b2).(b2 x b3))
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return float(np.degrees(np.arctan2(np.linalg.norm(b2) * (b1 @ n2), n1 @ n2)))


def main():
    res = dict(written_by='md/verify_system.py', checks={})
    ok = True

    def check(name, passed, detail):
        nonlocal ok
        res['checks'][name] = dict(passed=bool(passed), **detail)
        ok = ok and bool(passed)
        print(('ok   ' if passed else 'FAIL ') + name)

    prot = read_pdb(OUT / 'protein_md.pdb')
    zn = read_pdb(OUT / 'zinc.pdb')[0]['xyz']
    P = {(a['num'], a['name']): a['xyz'] for a in prot}
    og1, cz, od1, tn = P[(199, 'OG1')], P[(131, 'CZ')], P[(67, 'OD1')], P[(199, 'N')]
    table = {r['start']: r for r in csv.DictReader(open(OUT / 'starts.csv'))}
    wall = read_pdb(OUT / 'waters_all.pdb')
    lig = {t: {a['name']: a['xyz'] for a in read_pdb(OUT / f'ligand_{t}.pdb')} for t in 'XNR'}

    # 1. numbers of the starts
    worst, worst_ang, recomputed = 0.0, 0.0, {}
    for t in 'XNR':
        L = lig[t]
        t1, t2 = tors(L['C3'], L['C4'], L['S'], L['N2']), tors(L['C5'], L['C4'], L['S'], L['N2'])
        mine = dict(Zn_N2=d(zn, L['N2']), Zn_O_min=min(d(zn, L['O1']), d(zn, L['O2'])), N2_OG1=d(L['N2'], og1),
                    H_OG1=d(L['HN21'], og1), N2_H_OG1=ang(L['N2'], L['HN21'], og1),
                    T199N_O_min=min(d(tn, L['O1']), d(tn, L['O2'])), C3_C4_S_N2=t1, C5_C4_S_N2=t2,
                    tau_fold=max(abs(t1), abs(t2)), C4_S_N2=ang(L['C4'], L['S'], L['N2']), S_N2=d(L['S'], L['N2']),
                    S_N2_Zn=ang(L['S'], L['N2'], zn), N1_Phe131CZ=d(L['N1'], cz), N1_Asn67OD1=d(L['N1'], od1))
        heavy = [n for n in NAMES if not n.startswith('H')]
        alts = [dict(), dict(C2='C6', C6='C2', C3='C5', C5='C3'), dict(O1='O2', O2='O1'),
                dict(C2='C6', C6='C2', C3='C5', C5='C3', O1='O2', O2='O1')]
        mine['rmsd_to_crystal_A'] = min(float(np.sqrt(np.mean([((L[n] - lig['X'][m.get(n, n)]) ** 2).sum() for n in heavy])))
                                        for m in alts)
        for k, v in mine.items():
            if k in ANGLES:
                worst_ang = max(worst_ang, abs(v - float(table[t][k])))
            else:
                worst = max(worst, abs(v - float(table[t][k])))
        recomputed[t] = {k: round(v, 3) for k, v in mine.items()}
    check('start numbers recomputed from the written files', worst < TOL and worst_ang < TOL_ANG,
          dict(largest_difference_distances_A=round(worst, 4), largest_difference_angles_deg=round(worst_ang, 4),
               tolerance=dict(distances_A=TOL, angles_deg=TOL_ANG), recomputed=recomputed,
               note='the builder computes from unrounded coordinates, this check from the PDB files with three '
                    'decimals; an angle at a hydrogen 1.02 A from its nitrogen can differ by 0.05 deg from that alone'))

    # 2. provenance of the heavy atoms
    ref, _, _, _ = read_sdf(DOCK / 'T3_6RL9' / 'ref_SAN.sdf')
    dx = max(d(lig['X'][n], p) for n, (_, p) in zip(DEPOSITED_ORDER, ref))
    top = (DOCK / 'results' / 'poses' / 'T3_6RL9' / 'sfa_neutral__ad4zn__top_pose_per_seed.sdf').read_text().split('$$$$\n')
    rec = [i for i, r in enumerate(top) if r.strip() and r.split('\n')[0].strip().endswith('seed 3')][0]
    a3, _, p3, _ = read_sdf(DOCK / 'results' / 'poses' / 'T3_6RL9' / 'sfa_neutral__ad4zn__top_pose_per_seed.sdf', rec)
    hv = [(e, p) for e, p in a3 if e != 'H']
    dn = max(min(d(lig['N'][n], p) for e, p in hv if e == n[0]) for n in NAMES[:11])
    pq = []
    cur = None
    for ln in (DOCK / 'T3_6RL9' / 'runs' / 'sfa_neutral__ad4zn__prod' / 'seed14' / 'out.pdbqt').read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = []
        elif ln.startswith('ATOM') and cur is not None:
            cur.append((ln[77:79].strip(), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
        elif ln.startswith('ENDMDL'):
            pq.append(cur)
    hv = [(e, p) for e, p in pq[2] if e != 'HD']
    dr = max(min(d(lig['R'][n], p) for e, p in hv) for n in NAMES[:11])
    # aniline hydrogens of N and R are the docked ones
    hd = [p for e, p in pq[2] if e == 'HD']
    dh = max(min(d(lig['R'][n], p) for p in hd) for n in ('HN11', 'HN12'))
    check('heavy atoms are those of the sources', max(dx, dn, dr, dh) < TOL,
          dict(X_vs_ref_SAN_sdf_A=round(dx, 4), N_vs_exported_seed3_top_pose_A=round(dn, 4),
               N_record_score=p3.get('score_kcal_mol'), R_vs_seed14_model3_A=round(dr, 4),
               R_aniline_H_vs_docked_A=round(dh, 4)))

    # 3. water of each start
    wdet, wok = {}, True
    for t in 'XNR':
        heavy = np.array([lig[t][n] for n in NAMES[:11]])
        clash = sorted(w['num'] for w in wall if np.min(np.linalg.norm(heavy - w['xyz'], axis=1)) < 2.2)
        have = {w['num'] for w in read_pdb(OUT / f'waters_{t}.pdb')}
        want = {w['num'] for w in wall} - set(clash)
        wok = wok and have == want and int(table[t]['n_waters']) == len(have)
        wdet[t] = dict(removed=clash, n=len(have))
    check('water within 2.2 A of the ligand removed, nothing else', wok, dict(per_start=wdet, kept_before_ligand=len(wall)))

    # 4. heavy-atom contacts in each complete start
    cdet, cok = {}, True
    for t in 'XNR':
        atoms = [(f'{a["res"]}{a["num"]}', a['num'], a['name'], a['xyz'], 'P') for a in prot]
        atoms.append(('ZN301', 301, 'ZN', zn, 'Z'))
        atoms += [('SAN304', 304, n, lig[t][n], 'L') for n in NAMES[:11]]
        atoms += [(f'HOH{w["num"]}', w['num'], 'O', w['xyz'], 'W') for w in read_pdb(OUT / f'waters_{t}.pdb')]
        X = np.array([a[3] for a in atoms])
        short = []
        for i in range(len(atoms)):
            dd = np.linalg.norm(X[i + 1:] - X[i], axis=1)
            for j in np.where(dd < 2.2)[0]:
                a, b = atoms[i], atoms[i + 1 + int(j)]
                if a[4] == 'P' and b[4] == 'P' and (a[1] == b[1] or (abs(a[1] - b[1]) <= 2 and dd[j] < 1.45
                                                                    and {a[2], b[2]} == {'C', 'N'})):
                    continue      # same residue, or the peptide bond (125-127 are neighbours)
                if a[4] == 'L' and b[4] == 'L':
                    continue
                short.append((a[0] + '.' + a[2], b[0] + '.' + b[2], round(float(dd[j]), 2)))
        allowed = [s for s in short if 'ZN301.ZN' in s[:2] and (s[0].split('.')[1] in ('NE2', 'ND1') or s[1] == 'SAN304.N2')]
        other = [s for s in short if s not in allowed]
        cok = cok and not other and len(allowed) == (4 if t == 'X' else 3)
        cdet[t] = dict(zinc_bonds_below_2p2A=allowed, other_contacts_below_2p2A=other)
    check('no heavy-atom contact below 2.2 A apart from the zinc bonds', cok,
          dict(per_start=cdet, note='X: four Zn-N bonds below 2.2 A; N and R: three, the docked Zn-N(sulfonamide) is 2.21-2.24 A'))

    # 5. ligand files
    ldet, lok = {}, True
    for t in 'XNR':
        atoms, chg, props, nb = read_sdf(OUT / f'ligand_{t}.sdf')
        names = props['atom_names'].split()
        dmax = max(d(p, lig[t][n]) for n, (_, p) in zip(names, atoms))
        els = ''.join(e for e, _ in atoms)
        good = (len(atoms) == 18 and nb == 18 and names == NAMES and chg == {NAMES.index('N2') + 1: -1}
                and dmax < 1e-3 and els == 'CCCCCCNNOOS' + 'H' * 7 and list(lig[t]) == NAMES)
        for h, parent in (('HN11', 'N1'), ('HN12', 'N1'), ('H2', 'C2'), ('H3', 'C3'), ('H5', 'C5'), ('H6', 'C6'), ('HN21', 'N2')):
            good = good and 0.95 < d(lig[t][h], lig[t][parent]) < 1.12
        lok = lok and good
        ldet[t] = dict(atoms=len(atoms), bonds=nb, charge=chg, sdf_vs_pdb_A=round(dmax, 4),
                       N2_HN21_A=round(d(lig[t]['N2'], lig[t]['HN21']), 3))
    check('ligand files: 18 atoms, charge -1 on N2, names, SDF equals PDB', lok, dict(per_start=ldet))

    # for information (no pass/fail): contacts of the ligand hydrogens, with the hydrogens tleap gave the protein
    ph = OUT / 'tleap' / 'protein_H.pdb'
    if ph.exists():
        pa = read_pdb(ph)
        info = {}
        for t in 'XNR':
            hs = [(n, lig[t][n]) for n in NAMES[11:]]
            wat = read_pdb(OUT / f'waters_{t}.pdb')
            hw = min((d(p, w['xyz']), n, w['num']) for n, p in hs for w in wat)
            lp = sorted((d(p, a['xyz']), n, f'{a["res"]}{a["num"]}.{a["name"]}') for n, p in lig[t].items() for a in pa)[:3]
            info[t] = dict(closest_ligand_H_to_water_O=dict(A=round(hw[0], 2), H=hw[1], water=hw[2]),
                           closest_ligand_protein_pairs=[dict(A=round(x[0], 2), ligand=x[1], protein=x[2]) for x in lp])
        res['ligand_contacts_with_tleap_hydrogens'] = dict(
            note='not a gate; the hydrogens of protein and water are as tleap places them, before any minimisation',
            per_start=info)

    # 6. the build is reproducible
    def hashes():
        fs = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'verify.json' and 'tleap' not in p.parts)
        return {str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in fs}
    before = hashes()
    for script in ('build_system.py', 'system_audit.py'):
        p = subprocess.run(['nice', '-n', '10', sys.executable, str(REPO / 'md' / script)], capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stderr[-800:])
    after = hashes()
    changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
    check('second run of build_system.py and system_audit.py gives identical files', not changed,
          dict(files=len(after), changed=changed))
    res['sha256'] = after
    res['all_passed'] = bool(ok)
    (OUT / 'verify.json').write_text(json.dumps(res, indent=1) + '\n')
    print('all checks passed' if ok else 'SOME CHECKS FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
