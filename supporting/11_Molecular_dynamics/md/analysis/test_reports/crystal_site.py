# throwaway: the zinc site of the crystal start (6RL9 frame) from the files of md/system/, with plain formulas
# (acos of the normalised dot product; nothing imported from analyse_md.py or make_test_trajectory.py).
# Writes crystal_site.json: the six N-Zn-N angles (the reference values of gate G-Z / G-V in SPEC 1.1), the four
# Zn-N distances, the two ligand internal coordinates of GV4 and the two distances from the aniline N to Asn67.
import hashlib, json, math, sys
from pathlib import Path

SYS = Path(sys.argv[1])          # md/system
OUT = Path(sys.argv[2])


def atoms(path):
    out = {}
    for line in open(path):
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[12:16].strip(), line[22:26].strip())
            assert key not in out, key
            out[key] = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
    return out


def sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def norm(a):
    return math.sqrt(dot(a, a))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dist(a, b):
    return norm(sub(a, b))


def angle(a, b, c):
    u, v = sub(a, b), sub(c, b)
    return math.degrees(math.acos(dot(u, v) / (norm(u) * norm(v))))


def torsion(a, b, c, d):
    b1, b2, b3 = sub(b, a), sub(c, b), sub(d, c)
    n1, n2 = cross(b1, b2), cross(b2, b3)
    return math.degrees(math.atan2(norm(b2) * dot(b1, n2), dot(n1, n2)))


res = {}
for pfile in ('protein_md.pdb', 'tleap/protein_H.pdb'):
    P, Z, L = atoms(SYS / pfile), atoms(SYS / 'zinc.pdb'), atoms(SYS / 'ligand_X.pdb')
    zn = Z[('ZN', '301')]
    four = (('N2', 'N2 of SAN 304 (sulfonamide N)', L[('N2', '304')]), ('H94NE2', 'NE2 of His94', P[('NE2', '94')]),
            ('H96NE2', 'NE2 of His96', P[('NE2', '96')]), ('H119ND1', 'ND1 of His119', P[('ND1', '119')]))
    ang = []
    for i in range(4):
        for j in range(i + 1, 4):
            v = angle(four[i][2], zn, four[j][2])
            ang.append({'column': f'ang_{four[i][0]}_Zn_{four[j][0]}_deg', 'atoms': f'{four[i][1]} - Zn 301 - {four[j][1]}',
                        'value_deg': round(v, 4), 'one_decimal_deg': round(v, 1),
                        'distance_to_rounding_edge_deg': round(0.05 - abs(v - round(v, 1)), 4)})
    t3 = torsion(L[('C3', '304')], L[('C4', '304')], L[('S', '304')], L[('N2', '304')])
    t5 = torsion(L[('C5', '304')], L[('C4', '304')], L[('S', '304')], L[('N2', '304')])
    res[pfile] = {
        'angles': ang,
        'd_Zn_N_A': {k: round(dist(zn, x), 4) for k, _, x in four},
        'ang_C4_S_N2_deg': round(angle(L[('C4', '304')], L[('S', '304')], L[('N2', '304')]), 4),
        'tors_C3_C4_S_N2_deg': round(t3, 4), 'tors_C5_C4_S_N2_deg': round(t5, 4), 'tau_fold_deg': round(max(abs(t3), abs(t5)), 4),
        'd_Zn_O1_A': round(dist(zn, L[('O1', '304')]), 4), 'd_Zn_O2_A': round(dist(zn, L[('O2', '304')]), 4),
        'asn67': {s: {'d_N1_OD1_A': round(dist(atoms(SYS / f'ligand_{s}.pdb')[('N1', '304')], P[('OD1', '67')]), 4),
                      'd_N1_ND2_A': round(dist(atoms(SYS / f'ligand_{s}.pdb')[('N1', '304')], P[('ND2', '67')]), 4)} for s in 'XNR'}}
a, b = res['protein_md.pdb'], res['tleap/protein_H.pdb']
same = all(abs(x['value_deg'] - y['value_deg']) < 1e-9 for x, y in zip(a['angles'], b['angles']))
out = {'what': 'zinc site of the crystal start X (6RL9 frame), from md/system/zinc.pdb, ligand_X.pdb and the protein file named',
       'inputs_sha256': {f: hashlib.sha256((SYS / f).read_bytes()).hexdigest() for f in
                         ('zinc.pdb', 'ligand_X.pdb', 'ligand_N.pdb', 'ligand_R.pdb', 'protein_md.pdb', 'tleap/protein_H.pdb')},
       'the_two_protein_files_give_the_same_angles': same, 'from': res}
OUT.write_text(json.dumps(out, indent=1) + '\n')
for x in a['angles']:
    print(f"{x['column']:28s} {x['atoms']:62s} {x['value_deg']:9.4f}  -> {x['one_decimal_deg']:.1f}   (edge {x['distance_to_rounding_edge_deg']})")
print('same from both protein files:', same)
print({k: v for k, v in a.items() if k != 'angles'})
