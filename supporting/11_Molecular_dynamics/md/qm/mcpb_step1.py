"""Stage S2, step 1: input of MCPB.py for the zinc site of hCA II with sulfanilamide, and MCPB.py step 1.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/mcpb_step1.py

Reads (nothing is changed there)
    md/system/tleap/protein_H.pdb   protein of 6RL9 with the hydrogens of tleap, deposited residue numbers
    md/system/zinc.pdb              the zinc ion
    md/system/ligand_X.pdb          the sulfonamidate anion in the crystal pose, 18 atoms, names of SAN
    md/system/residue_map.csv       deposited residue number -> index tleap assigns
    md/ligand/sfa_anion.mol2        GAFF2 types and AM1-BCC charges of the anion (sum -1)
    md/ligand/sfa_anion_gaff2.frcmod, md/ligand/sfa_anion_repair.frcmod   ligand terms of stage S1

Writes into md/qm/mcpb/ (the working folder of MCPB.py)
    CA2_H.pdb       protein + ZN + SAN, residues numbered as tleap numbers them (ACE = 1), atoms 1..N
    ZN.pdb, ZN.mol2 the ion (metalpdb2mol2.py, charge +2)
    SAN.mol2        copy of sfa_anion.mol2 (MCPB.py wants <residue name>.mol2)
    sfa_anion_gaff2.frcmod, sfa_anion_repair.frcmod   copies
    CA2.in          the MCPB.py input
    step1.log       screen output of MCPB.py -s 1
    CA2_small.pdb, CA2_small_opt.com, CA2_small_fc.com, CA2_small.res, CA2_standard.pdb,
    CA2_standard.fingerprint, CA2_large.pdb, CA2_large.fingerprint, CA2_large_mk.com (and the GAMESS inputs)
and md/qm/step1_report.json: what the models hold (atoms, caps, charge, multiplicity, bonded atoms, distances).

Nothing is random.  A second run gives the same files (checked on 6 Oct 2026 in a second folder, set with the
environment variable ZINC_STEP1_WORK; all files of step 1 identical).  The script refuses to run over a folder
that already holds results of the quantum-chemical jobs.
"""
import csv
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
MD = REPO / 'md'
QM = MD / 'qm'
# ZINC_STEP1_WORK: another work folder, for a reproducibility check that leaves md/qm/mcpb alone; the report
# then goes beside that folder
WORK = pathlib.Path(os.environ.get('ZINC_STEP1_WORK', QM / 'mcpb'))
REPORT = (QM if 'ZINC_STEP1_WORK' not in os.environ else WORK.parent) / 'step1_report.json'
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
BIN = pathlib.Path(AMBERHOME) / 'bin'
GROUP = 'CA2'
# large_opt 0: the hydrogens of the large model are not optimised (default of MCPB.py; see md/ZINC.md, section 4)
LARGE_OPT = 0
CUTOFF = 2.8

Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'ZN': 30, 'Zn': 30}


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def read_pdb_atoms(path):
    out = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            out.append(dict(rec=line[:6].strip(), name=line[12:16].strip(), resn=line[17:20].strip(),
                            chain=line[21], resi=int(line[22:26]),
                            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                            elem=line[76:78].strip()))
    return out


def pdb_line(rec, serial, name, resn, resi, xyz, elem):
    nm = name if len(name) == 4 else ' ' + name.ljust(3)
    if elem.upper() == 'ZN':
        nm = 'ZN  '
    return '%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s' % (
        rec, serial, nm, resn, 'A', resi, xyz[0], xyz[1], xyz[2], 1.0, 0.0, elem.upper().rjust(2))


def build_input_pdb():
    prot = read_pdb_atoms(MD / 'system' / 'tleap' / 'protein_H.pdb')
    zn = read_pdb_atoms(MD / 'system' / 'zinc.pdb')
    lig = read_pdb_atoms(MD / 'system' / 'ligand_X.pdb')
    assert len(zn) == 1 and len(lig) == 18
    resmap = {}
    with open(MD / 'system' / 'residue_map.csv') as fh:
        for row in csv.DictReader(fh):
            if row['tleap_index'] and row['name_md'] not in ('WAT', 'HOH'):
                resmap.setdefault(int(row['deposited_number']), int(row['tleap_index']))
    lines, serial, seq, last, table = [], 0, 0, None, {}
    for a in prot:
        if a['resi'] != last:
            seq += 1
            last = a['resi']
            # the index of tleap is the running number of the residue in the file
            assert resmap[a['resi']] == seq, (a['resi'], resmap[a['resi']], seq)
            table[a['resi']] = seq
        serial += 1
        lines.append(pdb_line('ATOM', serial, a['name'], a['resn'], seq, a['xyz'], a['elem']))
    nprot, nres = serial, seq
    lines.append('TER')
    serial += 1
    zn_serial, zn_res = serial, nres + 1
    lines.append(pdb_line('HETATM', serial, 'ZN', 'ZN', zn_res, zn[0]['xyz'], 'ZN'))
    lines.append('TER')
    lig_res = nres + 2
    lig_serials = {}
    for a in lig:
        serial += 1
        lig_serials[a['name']] = serial
        lines.append(pdb_line('HETATM', serial, a['name'], 'SAN', lig_res, a['xyz'], a['elem']))
    lines += ['TER', 'END']
    (WORK / f'{GROUP}_H.pdb').write_text('\n'.join(lines) + '\n')
    (WORK / 'ZN.pdb').write_text(pdb_line('HETATM', 1, 'ZN', 'ZN', 1, zn[0]['xyz'], 'ZN') + '\nEND\n')
    return dict(n_protein_atoms=nprot, n_protein_residues=nres, zn_serial=zn_serial, zn_residue=zn_res,
                ligand_residue=lig_res, ligand_serials=lig_serials, n_atoms=serial,
                seq_of={str(k): table[k] for k in (94, 96, 119, 199, 106)})


def run(cmd, logname):
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{BIN}:' + os.environ.get('PATH', ''),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    p = subprocess.run(cmd, cwd=WORK, env=env, capture_output=True, text=True)
    (WORK / logname).write_text(p.stdout + p.stderr)
    if p.returncode != 0:
        sys.exit(f'{cmd[0]} failed, see {WORK / logname}')
    return p.stdout + p.stderr


def read_com(path):
    """charge, multiplicity, [(element, x, y, z)] of a Gaussian input written by MCPB.py."""
    lines = pathlib.Path(path).read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip() == 'CLR')
    chg, mult = [int(x) for x in lines[i + 2].split()]
    atoms = []
    for l in lines[i + 3:]:
        t = l.split()
        if len(t) not in (4, 5):
            break
        atoms.append((t[0], *[float(x) for x in t[-3:]]))
    return chg, mult, atoms, lines


def dist(a, b):
    return float(np.linalg.norm(np.array(a) - np.array(b)))


def angle(a, b, c):
    u, v = np.array(a) - np.array(b), np.array(c) - np.array(b)
    return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    info = build_input_pdb()
    # the ion and the ligand
    run([str(BIN / 'python'), str(BIN / 'metalpdb2mol2.py'), '-i', 'ZN.pdb', '-o', 'ZN.mol2', '-c', '2'],
        'metalpdb2mol2.log')
    shutil.copyfile(MD / 'ligand' / 'sfa_anion.mol2', WORK / 'SAN.mol2')
    for f in ('sfa_anion_gaff2.frcmod', 'sfa_anion_repair.frcmod'):
        shutil.copyfile(MD / 'ligand' / f, WORK / f)
    inp = [f'original_pdb {GROUP}_H.pdb', f'group_name {GROUP}', f'cut_off {CUTOFF}',
           f'ion_ids {info["zn_serial"]}', 'ion_mol2files ZN.mol2', 'naa_mol2files SAN.mol2',
           'frcmod_files sfa_anion_gaff2.frcmod sfa_anion_repair.frcmod', 'force_field ff14SB', 'gaff 2',
           'water_model TIP3P', f'large_opt {LARGE_OPT}', 'software_version gau']
    (WORK / f'{GROUP}.in').write_text('\n'.join(inp) + '\n')
    # results of the quantum-chemical jobs live in the same folder: never run step 1 over them
    made = [WORK / f'{GROUP}_large_mk.log', WORK / f'{GROUP}_small_opt.fchk']
    if any(f.exists() for f in made) and '--force' not in sys.argv:
        sys.exit(f'{[f.name for f in made if f.exists()]} exist in {WORK}: step 1 was done and its models are in '
                 'use.  Run with --force only if the models are to be rebuilt (the adapter files stay, but they '
                 'then belong to the earlier models).')
    for name in ('small.pdb', 'small.res', 'small_opt.com', 'small_opt.inp', 'small_fc.com', 'small_fc.inp',
                 'standard.pdb', 'standard.fingerprint', 'large.pdb', 'large.fingerprint', 'large_mk.com',
                 'large_mk.inp'):
        f = WORK / f'{GROUP}_{name}'
        if f.exists():
            f.unlink()
    out = run([str(BIN / 'python'), str(BIN / 'MCPB.py'), '-i', f'{GROUP}.in', '-s', '1'], 'step1.log')

    # ---- inspection -------------------------------------------------------------------------------------------
    full = read_pdb_atoms(WORK / f'{GROUP}_H.pdb')
    znxyz = full[info['zn_serial'] - 1]['xyz']
    near = []
    for k, a in enumerate(full, 1):
        if k == info['zn_serial']:
            continue
        d = dist(a['xyz'], znxyz)
        if d <= 3.2 and a['name'][0] in 'NOS':
            near.append(dict(serial=k, residue=f"{a['resn']}{a['resi']}", atom=a['name'], distance=round(d, 3),
                             bonded_by_mcpb=bool(d <= CUTOFF)))
    near.sort(key=lambda r: r['distance'])
    links = [l.split()[1:] for l in (WORK / f'{GROUP}_standard.fingerprint').read_text().splitlines()
             if l.startswith('LINK')]
    types = [l.split() for l in (WORK / f'{GROUP}_standard.fingerprint').read_text().splitlines()
             if not l.startswith('LINK')]
    renamed = [dict(atom=t[0], serial=int(t[1]), type_before=t[2], type_after=t[4]) for t in types if t[2] != t[4]]

    report = dict(input=info, mcpb_input=inp, zinc_neighbours_within_3p2=near, links=links, renamed_types=renamed,
                  small_res=(WORK / f'{GROUP}_small.res').read_text().split('\n'))
    for tag, com in (('small', f'{GROUP}_small_opt.com'), ('large', f'{GROUP}_large_mk.com')):
        chg, mult, atoms, lines = read_com(WORK / com)
        pdb = read_pdb_atoms(WORK / f'{GROUP}_{tag}.pdb')
        assert len(pdb) == len(atoms)
        nel = sum(Z[a[0]] for a in atoms) - chg
        comp = {}
        for a in atoms:
            comp[a[0]] = comp.get(a[0], 0) + 1
        residues = []
        for a in pdb:
            key = f"{a['resn']}{a['resi']}"
            if not residues or residues[-1][0] != key:
                residues.append([key, 0])
            residues[-1][1] += 1
        izn = next(k for k, a in enumerate(atoms) if a[0].upper() == 'ZN')
        zn_n = []
        for k, (a, p) in enumerate(zip(atoms, pdb)):
            if a[0] == 'N' and dist(a[1:], atoms[izn][1:]) <= CUTOFF:
                zn_n.append(dict(index=k + 1, residue=f"{p['resn']}{p['resi']}", atom=p['name'],
                                 distance_com=round(dist(a[1:], atoms[izn][1:]), 3),
                                 distance_pdb=round(dist(p['xyz'], pdb[izn]['xyz']), 3)))
        angs = []
        for i in range(len(zn_n)):
            for j in range(i + 1, len(zn_n)):
                a, b = pdb[zn_n[i]['index'] - 1], pdb[zn_n[j]['index'] - 1]
                angs.append(dict(a=f"{zn_n[i]['residue']}@{zn_n[i]['atom']}", b=f"{zn_n[j]['residue']}@{zn_n[j]['atom']}",
                                 angle=round(angle(a['xyz'], pdb[izn]['xyz'], b['xyz']), 1)))
        report[tag] = dict(com=com, n_atoms=len(atoms), composition=comp, charge=chg, multiplicity=mult,
                           electrons=nel, residues=residues, route=[l for l in lines if l.startswith(('#', 'IOp'))],
                           zinc_index=izn + 1, zinc_nitrogens=zn_n, n_zn_n_angles=angs,
                           tail=[l for l in lines[-4:] if l.strip()])
    report['large']['mk_tail_lines'] = [l for l in (WORK / f'{GROUP}_large_mk.com').read_text().splitlines()[-6:]
                                        if l.strip()]
    report['software'] = dict(ambertools=subprocess.run([str(BIN / 'tleap'), '-h'], capture_output=True, text=True,
                                                        env=dict(os.environ, AMBERHOME=AMBERHOME)).stdout[:0] or
                              'AmberTools 26 (conda-forge), MCPB.py version 2026', python=sys.version.split()[0])
    report['sha256'] = {(str(p.relative_to(REPO)) if REPO in p.parents else f'<work>/{p.name}'): sha256(p) for p in sorted(
        [MD / 'system' / 'tleap' / 'protein_H.pdb', MD / 'system' / 'zinc.pdb', MD / 'system' / 'ligand_X.pdb',
         MD / 'ligand' / 'sfa_anion.mol2', MD / 'ligand' / 'sfa_anion_gaff2.frcmod',
         MD / 'ligand' / 'sfa_anion_repair.frcmod'] +
        [WORK / f for f in (f'{GROUP}_H.pdb', f'{GROUP}.in', 'ZN.mol2', 'SAN.mol2', f'{GROUP}_small.pdb',
                            f'{GROUP}_small_opt.com', f'{GROUP}_small_fc.com', f'{GROUP}_standard.pdb',
                            f'{GROUP}_standard.fingerprint', f'{GROUP}_large.pdb', f'{GROUP}_large.fingerprint',
                            f'{GROUP}_large_mk.com')])}
    REPORT.write_text(json.dumps(report, indent=1) + '\n')
    print(out[-1500:])
    print('small:', report['small']['n_atoms'], 'atoms, charge', report['small']['charge'], 'mult',
          report['small']['multiplicity'])
    print('large:', report['large']['n_atoms'], 'atoms, charge', report['large']['charge'], 'mult',
          report['large']['multiplicity'])


if __name__ == '__main__':
    main()
