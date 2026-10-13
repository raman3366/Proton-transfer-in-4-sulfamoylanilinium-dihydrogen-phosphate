"""Stage S1 of md/PLAN.md, side table: which internal coordinate of the ligand separates the two ring orientations.

  /root/cqc-sfa-tools/dock/bin/python md/ligand_bound_rotamer.py        (docking environment: gemmi and RDKit)

Reads the 1.0 A structure (dock/pdb/6RL9.cif.gz, SAN 304, zinc 301) and the docked poses of the paper
(dock/results/poses/T3_6RL9/: AutoDock4Zn with the DFT input geometry, with C-S-N 109.5 and with C-S-N 102) and
measures, for each, the rotation about N-S as seen from the zinc (Zn-N2-S-C4 and Zn-N2-S-O), the angles Zn-N2-S and
C4-S-N2, the ring torsion C3-C4-S-N2 and the two distances that md/PLAN.md fixes as the orientation coordinate.
Writes md/ligand/bound_rotamer.csv (one row per structure) and md/ligand/bound_rotamer.json (ranges per set).
Nothing is fitted here; the table shows that the two orientations differ in the rotation about N-S, a coordinate
that the tests of the free ligand do not reach.
"""
import json, math, pathlib
import numpy as np
import gemmi
from rdkit import Chem

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'md' / 'ligand'
POSES = ROOT / 'dock' / 'results' / 'poses' / 'T3_6RL9'
SETS = [('docked, AutoDock4Zn, input C-S-N 109.5 (crystal-like orientation)', 'sfa_neutral__ad4zn__csn109__top_pose_per_seed.sdf'),
        ('docked, AutoDock4Zn, input C-S-N 102 (rotated orientation)', 'sfa_neutral__ad4zn__csn102__top_pose_per_seed.sdf'),
        ('docked, AutoDock4Zn, production (DFT input geometry)', 'sfa_neutral__ad4zn__top_pose_per_seed.sdf')]


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return math.degrees(math.atan2(float(np.cross(b1, v) @ w), float(v @ w)))


def ang(a, b, c):
    u, v = a - b, c - b
    return math.degrees(math.acos(float(u @ v) / np.linalg.norm(u) / np.linalg.norm(v)))


def main():
    st = gemmi.read_structure(str(ROOT / 'dock' / 'pdb' / '6RL9.cif.gz'))
    model = st[0]

    def pos(res, seq, atom):
        hits = [np.array(a.pos.tolist()) for ch in model for r in ch if r.name == res and r.seqid.num == seq
                for a in r if a.name == atom and a.altloc in ('\0', 'A')]
        assert len(hits) == 1, (res, seq, atom, len(hits))
        return hits[0]
    zn, cz, od1, og1 = pos('ZN', 301, 'ZN'), pos('PHE', 131, 'CZ'), pos('ASN', 67, 'OD1'), pos('THR', 199, 'OG1')
    rows = []

    def row(label, k, X, hpos):
        r = dict(set=label, pose=k,
                 Zn_N2=float(np.linalg.norm(zn - X['N2'])), Zn_N2_S=ang(zn, X['N2'], X['S']),
                 Zn_N2_S_C4=dih(zn, X['N2'], X['S'], X['C4']),
                 Zn_N2_S_Oa=dih(zn, X['N2'], X['S'], X['O1']), Zn_N2_S_Ob=dih(zn, X['N2'], X['S'], X['O2']),
                 Zn_O_shortest=float(min(np.linalg.norm(zn - X['O1']), np.linalg.norm(zn - X['O2']))),
                 C4_S_N2=ang(X['C4'], X['S'], X['N2']), C3_C4_S_N2=dih(X['C3'], X['C4'], X['S'], X['N2']),
                 N1_CZ_Phe131=float(np.linalg.norm(X['N1'] - cz)), N1_OD1_Asn67=float(np.linalg.norm(X['N1'] - od1)),
                 OG1_Thr199_N2=float(np.linalg.norm(og1 - X['N2'])),
                 OG1_N2_S_C4=dih(og1, X['N2'], X['S'], X['C4']))
        r['tau_fold'] = max(abs(r['C3_C4_S_N2']), abs(dih(X['C5'], X['C4'], X['S'], X['N2'])))
        hs = sorted(dih(h, X['N2'], X['S'], X['C4']) for h in hpos)
        r['H_N2_S_C4'] = ' '.join(f'{h:.0f}' for h in hs)
        rows.append(r)

    # crystal ligand
    X = {a: pos('SAN', 304, a) for a in ('C1', 'N1', 'C2', 'C3', 'C4', 'C5', 'C6', 'S', 'N2', 'O1', 'O2')}
    row('crystal 6RL9, SAN 304', 0, X, [])
    patt = Chem.MolFromSmarts('Nc1ccc(cc1)S(=O)(=O)N')
    for label, fn in SETS:
        for k, mol in enumerate(Chem.SDMolSupplier(str(POSES / fn), removeHs=False)):
            mt = mol.GetSubstructMatch(patt)
            assert mt, fn
            P = mol.GetConformer().GetPositions()
            n1, c1, c2, c3, c4, c5, c6, s, o1, o2, n2 = mt
            X = dict(N1=P[n1], C1=P[c1], C2=P[c2], C3=P[c3], C4=P[c4], C5=P[c5], C6=P[c6], S=P[s], O1=P[o1], O2=P[o2],
                     N2=P[n2])
            hs = [P[a.GetIdx()] for a in mol.GetAtomWithIdx(n2).GetNeighbors() if a.GetAtomicNum() == 1]
            row(label, k, X, hs)
    cols = list(rows[0].keys())
    with open(OUT / 'bound_rotamer.csv', 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join(f'{r[c]:.3f}' if isinstance(r[c], float) else f'"{r[c]}"' if isinstance(r[c], str)
                             else str(r[c]) for c in cols) + '\n')
    summ = {}
    for label in [rows[0]['set']] + [s[0] for s in SETS]:
        sub = [r for r in rows if r['set'] == label]
        summ[label] = dict(n=len(sub), **{c: dict(min=min(r[c] for r in sub), max=max(r[c] for r in sub),
                                                  mean=float(np.mean([r[c] for r in sub])))
                                          for c in ('Zn_N2', 'Zn_N2_S', 'Zn_N2_S_C4', 'Zn_O_shortest', 'C4_S_N2',
                                                    'tau_fold', 'N1_CZ_Phe131', 'N1_OD1_Asn67', 'OG1_N2_S_C4')})
        # the oxygen that eclipses the zinc: smallest |Zn-N2-S-O| of the two
        summ[label]['Zn_N2_S_O_smallest_abs'] = dict(
            min=min(min(abs(r['Zn_N2_S_Oa']), abs(r['Zn_N2_S_Ob'])) for r in sub),
            max=max(min(abs(r['Zn_N2_S_Oa']), abs(r['Zn_N2_S_Ob'])) for r in sub))
    (OUT / 'bound_rotamer.json').write_text(json.dumps(dict(
        script='md/ligand_bound_rotamer.py', gemmi=gemmi.__version__, rdkit=Chem.rdBase.rdkitVersion,
        note='angles in degrees, distances in Angstrom; dihedral angles through the zinc are measured, no bond to '
             'the zinc is assumed', sets=summ), indent=1) + '\n')
    for label, v in summ.items():
        print(f"{label}: n={v['n']}  Zn-N2-S-C4 {v['Zn_N2_S_C4']['min']:.0f} to {v['Zn_N2_S_C4']['max']:.0f}; "
              f"smallest |Zn-N2-S-O| {v['Zn_N2_S_O_smallest_abs']['min']:.0f} to {v['Zn_N2_S_O_smallest_abs']['max']:.0f}; "
              f"Zn-N2-S {v['Zn_N2_S']['min']:.0f} to {v['Zn_N2_S']['max']:.0f}; tau_fold {v['tau_fold']['min']:.0f} to "
              f"{v['tau_fold']['max']:.0f}; N1..CZ {v['N1_CZ_Phe131']['min']:.1f} to {v['N1_CZ_Phe131']['max']:.1f}; "
              f"N1..OD1 {v['N1_OD1_Asn67']['min']:.1f} to {v['N1_OD1_Asn67']['max']:.1f}")


if __name__ == '__main__':
    main()
