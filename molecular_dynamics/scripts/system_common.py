"""Shared helpers of the system builder (stage S3 of md/PLAN.md): paths, reading of the deposited
model 6RL9, the rule for alternate conformers, small geometry functions and a PDB writer.

Used by md/system_audit.py (evidence), md/build_system.py (decisions -> files) and
md/check_tleap.py (test).  Run with the docking environment, which has gemmi and RDKit:

    /root/cqc-sfa-tools/dock/bin/python md/<script>.py

Nothing here writes a file.
"""
import hashlib
import pathlib

import gemmi
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
MD = REPO / 'md'
OUT = MD / 'system'
DOCK = REPO / 'dock'
CIF = DOCK / 'pdb' / '6RL9.cif.gz'
CHAIN = 'A'
PROD_JOB = DOCK / 'T3_6RL9' / 'runs' / 'sfa_neutral__ad4zn__prod'
DOCK_RECEPTOR = DOCK / 'T3_6RL9' / 'receptor.pdb'
DFT_ANION = REPO / 'molecules' / 'sfa' / 'results' / 'sulfonamidate_pcm_d3' / 'opt.xyz'

ZN_SEQ, LIG_SEQ = 301, 304
# heavy-atom names of SAN in 6RL9 (the order of the deposited residue), then the hydrogens of the anion
SAN_HEAVY = ['C1', 'N1', 'C2', 'C3', 'C4', 'C5', 'C6', 'S', 'N2', 'O1', 'O2']
SAN_BONDS = [('C1', 'N1', 1), ('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'C4', 2), ('C4', 'C5', 1),
             ('C5', 'C6', 2), ('C6', 'C1', 1), ('C4', 'S', 1), ('S', 'N2', 1), ('S', 'O1', 2), ('S', 'O2', 2)]
SAN_H = [('H2', 'C2'), ('H3', 'C3'), ('H5', 'C5'), ('H6', 'C6'), ('HN11', 'N1'), ('HN12', 'N1'),
         ('HN21', 'N2')]
# Atom order of the ligand files written here: the order of the paper's DFT files and of the ligand
# parameter files of stage S1 (md/ligand/sfa_anion.mol2), so that coordinates and parameters line up
# even for a tool that goes by position and not by name.
LIG_ORDER = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1', 'N2', 'O1', 'O2', 'S',
             'H2', 'H3', 'H5', 'H6', 'HN11', 'HN12', 'HN21']

# heavy atoms of the standard residues (PDB v3 names), to list what the deposited model lacks
HEAVY = {
    'ALA': 'N CA C O CB', 'ARG': 'N CA C O CB CG CD NE CZ NH1 NH2', 'ASN': 'N CA C O CB CG OD1 ND2',
    'ASP': 'N CA C O CB CG OD1 OD2', 'CYS': 'N CA C O CB SG', 'GLN': 'N CA C O CB CG CD OE1 NE2',
    'GLU': 'N CA C O CB CG CD OE1 OE2', 'GLY': 'N CA C O', 'HIS': 'N CA C O CB CG ND1 CD2 CE1 NE2',
    'ILE': 'N CA C O CB CG1 CG2 CD1', 'LEU': 'N CA C O CB CG CD1 CD2', 'LYS': 'N CA C O CB CG CD CE NZ',
    'MET': 'N CA C O CB CG SD CE', 'PHE': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'N CA C O CB CG CD',
    'SER': 'N CA C O CB OG', 'THR': 'N CA C O CB OG1 CG2', 'TRP': 'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'N CA C O CB CG1 CG2',
}
HEAVY = {k: v.split() for k, v in HEAVY.items()}
BACKBONE = ('N', 'CA', 'C', 'O')
COUPLING_CUT = 2.4      # A; atoms of two different labels closer than this cannot be present together


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def alt(atom):
    return atom.altloc if atom.altloc not in ('\x00', '') else '-'


def xyz(atom):
    return np.array(atom.pos.tolist())


def dist(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def angle(a, b, c):
    v1, v2 = np.asarray(a) - np.asarray(b), np.asarray(c) - np.asarray(b)
    return float(np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))


def dihedral(a, b, c, d):
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def place(a, b, c, bond, ang, dih):
    """Position d with |cd| = bond, angle b-c-d = ang and dihedral a-b-c-d = dih (degrees)."""
    a, b, c = (np.asarray(x, float) for x in (a, b, c))
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    t, p = np.radians(ang), np.radians(dih)
    d2 = np.array([-bond * np.cos(t), bond * np.sin(t) * np.cos(p), bond * np.sin(t) * np.sin(p)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def read_deposited(with_h=False):
    """The deposited model, chain A.  Returns (structure, chain)."""
    st = gemmi.read_structure(str(CIF))
    st.setup_entities()
    if not with_h:
        st.remove_hydrogens()
    return st, st[0][CHAIN]


def polymer_residues(chain):
    return [r for r in chain if r.het_flag == 'A' and r.name in HEAVY]


def waters(chain):
    return [r for r in chain if r.name == 'HOH']


def hetero(chain):
    return [r for r in chain if r.het_flag == 'H' and r.name != 'HOH']


def labels_of(res):
    return sorted({alt(a) for a in res} - {'-'})


def conformer_atoms(res, label):
    """Heavy atoms of one conformer: the unlabelled atoms and those of `label`."""
    return [a for a in res if alt(a) in ('-', label) and a.element.name not in ('H', 'D')]


def choose_conformers(chain):
    """One conformer per residue of the protein.

    Rule (the rule of dock/prep_receptor.py, with coupling added):
      1. residues with alternates form a group when they are neighbours in the sequence or when
         atoms of two different labels of the two residues come closer than COUPLING_CUT (they
         cannot be present together, so the labels are correlated); a group takes ONE label;
      2. the label of a group is the one with the most modelled heavy atoms, then the highest
         summed occupancy, then 'A'.
    Returns (choice {seqnum: label}, groups [list of dicts]).
    """
    res = [r for r in polymer_residues(chain) if labels_of(r)]
    parent = list(range(len(res)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    links = []
    order = [r.seqid.num for r in polymer_residues(chain)]
    pos = {n: i for i, n in enumerate(order)}
    for i, r in enumerate(res):
        for j in range(i + 1, len(res)):
            q = res[j]
            why = None
            if abs(pos[r.seqid.num] - pos[q.seqid.num]) == 1:
                why = 'neighbours in the sequence'
            else:
                best = None
                for a in r:
                    if alt(a) == '-' or a.element.name in ('H', 'D'):
                        continue
                    for b in q:
                        if alt(b) in ('-', alt(a)) or b.element.name in ('H', 'D'):
                            continue
                        d = a.pos.dist(b.pos)
                        if d < COUPLING_CUT and (best is None or d < best[0]):
                            best = (d, f'{r.name}{r.seqid.num}.{a.name}[{alt(a)}] - '
                                       f'{q.name}{q.seqid.num}.{b.name}[{alt(b)}]')
                if best:
                    why = f'{best[1]} {best[0]:.2f} A'
            if why:
                parent[find(i)] = find(j)
                links.append((r.seqid.num, q.seqid.num, why))
    groups = {}
    for i, r in enumerate(res):
        groups.setdefault(find(i), []).append(r)
    choice, out = {}, []
    for members in sorted(groups.values(), key=lambda m: m[0].seqid.num):
        labs = sorted({L for r in members for L in labels_of(r)})
        score = {}
        for L in labs:
            natoms = sum(len(conformer_atoms(r, L)) for r in members)
            occ = sum(next((a.occ for a in r if alt(a) == L), 0.0) for r in members)
            score[L] = (natoms, round(occ, 3), L == 'A')
        best = max(labs, key=lambda L: score[L])
        for r in members:
            choice[r.seqid.num] = best
        nums = [r.seqid.num for r in members]
        out.append(dict(
            residues=[f'{r.name}{r.seqid.num}' for r in members], label=best,
            occupancy={L: sorted({round(a.occ, 2) for r in members for a in r if alt(a) == L}) for L in labs},
            heavy_atoms={L: score[L][0] for L in labs},
            decided_by=('number of modelled atoms' if len({score[L][0] for L in labs}) > 1 and
                        score[best][0] > max(score[L][0] for L in labs if L != best)
                        else ('occupancy' if len({score[L][1] for L in labs}) > 1 else 'label A (tie)')),
            coupling=[w for a, b, w in links if a in nums and b in nums and
                      w != 'neighbours in the sequence'],
        ))
    return choice, out


def docking_conformers():
    """Labels the docking receptor took (dock/T3_6RL9/prep_report.json), {seqnum: label}."""
    import json
    import re
    rep = json.loads((DOCK / 'T3_6RL9' / 'prep_report.json').read_text())
    out = {}
    for k, v in rep['altloc_chosen'].items():
        m = re.match(r'A:([A-Z]{3})(\d+)$', k)
        if m and m.group(1) in HEAVY:
            out[int(m.group(2))] = v
    return out


def pdb_line(record, serial, name, resname, chain, resseq, pos, occ=1.0, b=0.0, element=None):
    element = element or ''.join(c for c in name if c.isalpha())[:1]
    field = f'{name:<4s}' if len(name) >= 4 or len(element) == 2 else f' {name:<3s}'
    return (f'{record:<6s}{serial % 100000:5d} {field} {resname:>3s} {chain}{resseq:4d}    '
            f'{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occ:6.2f}{b:6.2f}          {element.upper():>2s}')


def read_pdbqt_modes(path):
    """Modes of a Vina/AutoDock4Zn output: list of dicts(score, atoms={serial: (type, xyz)})."""
    modes, cur = [], None
    for ln in pathlib.Path(path).read_text().splitlines():
        if ln.startswith('MODEL'):
            cur = dict(score=None, atoms={})
        elif ln.startswith('REMARK VINA RESULT'):
            cur['score'] = float(ln.split()[3])
        elif ln.startswith(('ATOM', 'HETATM')):
            cur['atoms'][int(ln[6:11])] = (ln[77:79].strip(),
                                           np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
        elif ln.startswith('ENDMDL'):
            modes.append(cur)
    return modes


# The ligand PDBQT of the docking (dock/T3_6RL9/runs/*/seed*/ligand.pdbqt) numbers the atoms
# 1-6 ring (1 carries the aniline N, 4 the sulfur), 7 aniline N, 8-9 its H, 10 S, 11-12 O,
# 13 sulfonamide N, 14-15 its H.  The ring and the two sulfonyl oxygens are symmetric, so four
# assignments of the SAN names are possible; the one closest to the crystal pose is taken.
PDBQT_MAPS = []
for _ring in (('C1', 'C2', 'C3', 'C4', 'C5', 'C6'), ('C1', 'C6', 'C5', 'C4', 'C3', 'C2')):
    for _oo in (('O1', 'O2'), ('O2', 'O1')):
        PDBQT_MAPS.append({1: _ring[0], 2: _ring[1], 3: _ring[2], 4: _ring[3], 5: _ring[4], 6: _ring[5],
                           7: 'N1', 10: 'S', 11: _oo[0], 12: _oo[1], 13: 'N2'})
PDBQT_ANILINE_H, PDBQT_SULFONAMIDE_H = (8, 9), (14, 15)


def name_docked_atoms(mode, crystal):
    """SAN names for a docked mode: the symmetry-equivalent assignment with the lowest heavy-atom
    RMSD to the crystal pose (no superposition).  Returns (rmsd, {name: xyz}, map)."""
    best = None
    for mp in PDBQT_MAPS:
        M = {mp[s]: v[1] for s, v in mode['atoms'].items() if s in mp}
        r = float(np.sqrt(np.mean([np.sum((M[n] - crystal[n]) ** 2) for n in SAN_HEAVY])))
        if best is None or r < best[0] - 1e-12:
            best = (r, M, mp)
    return best


def site_atoms(chain, choice):
    """Reference atoms of the site in the deposited frame."""
    def P(num, name):
        for r in chain:
            if r.seqid.num == num:
                for a in r:
                    if a.name == name and alt(a) in ('-', choice.get(num, 'A')):
                        return xyz(a)
        raise KeyError((num, name))
    return dict(ZN=P(ZN_SEQ, 'ZN'), OG1=P(199, 'OG1'), T199N=P(199, 'N'), CZ131=P(131, 'CZ'),
                OD1_67=P(67, 'OD1'), ND2_67=P(67, 'ND2'), NE2_94=P(94, 'NE2'), NE2_96=P(96, 'NE2'), ND1_119=P(119, 'ND1'),
                OE1_106=P(106, 'OE1'), OE2_106=P(106, 'OE2'))


def pose_metrics(M, site, H=None):
    """Numbers reported for a ligand pose M {SAN name: xyz}; H = position of the sulfonamide hydrogen."""
    t1 = dihedral(M['C3'], M['C4'], M['S'], M['N2'])
    t2 = dihedral(M['C5'], M['C4'], M['S'], M['N2'])
    out = dict(
        Zn_N2=dist(M['N2'], site['ZN']),
        Zn_O1=dist(M['O1'], site['ZN']), Zn_O2=dist(M['O2'], site['ZN']),
        Zn_O_min=min(dist(M['O1'], site['ZN']), dist(M['O2'], site['ZN'])),
        N2_OG1=dist(M['N2'], site['OG1']),
        T199N_O_min=min(dist(M['O1'], site['T199N']), dist(M['O2'], site['T199N'])),
        C3_C4_S_N2=t1, C5_C4_S_N2=t2, tau_fold=max(abs(t1), abs(t2)),
        C4_S_N2=angle(M['C4'], M['S'], M['N2']), S_N2=dist(M['S'], M['N2']), C4_S=dist(M['C4'], M['S']),
        S_N2_Zn=angle(M['S'], M['N2'], site['ZN']),
        N1_Phe131CZ=dist(M['N1'], site['CZ131']), N1_Asn67OD1=dist(M['N1'], site['OD1_67']),
        N1_Asn67ND2=dist(M['N1'], site['ND2_67']),
    )
    if H is not None:
        out.update(H_OG1=dist(H, site['OG1']), N2_H_OG1=angle(M['N2'], H, site['OG1']),
                   H_Zn=dist(H, site['ZN']), Zn_N2_H=angle(site['ZN'], M['N2'], H),
                   S_N2_H=angle(M['S'], M['N2'], H))
    return out
