"""Gate G-R, step 2: force-field energies of the 64-atom small model at the geometries of md/qm/gr_points.py.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/gr_mm.py run \
        --params DIR --work WORK --out OUT.json --label TEXT [--extra FILE.frcmod ...] [--n2h structure|asis]
        [--points md/qm/gr/points.json] [--no-min]

  DIR    a folder as md/qm/after_qm.sh leaves it: HD1.mol2 HD2.mol2 HE1.mol2 ZN1.mol2 SN1.mol2 (the RESP
         charges and the atom types of the five residues of the site), CA2_mcpbpy.frcmod (step 2e or 2s),
         sfa_anion_gaff2.frcmod, sfa_anion_repair.frcmod; CA2_dry.prmtop if the comparison with the dry
         bonded complex is wanted (it is, whenever the file is there).
  WORK   scratch folder for tleap (templates, inputs, topologies); OUT.json the result (small).
  --extra  frcmod files loaded after CA2_mcpbpy.frcmod (for example the fitted terms of md/qm/gr_fit.py).
  --n2h    structure (default): the equilibrium value of Zn-N2-H is set to its value in start X (107.41 deg)
           by one more frcmod line, force constant unchanged; asis: left as MCPB.py wrote it (step 2e: 120).
           md/PLAN.md, 6 Oct 2026, 3:02 AM, point 2: "The 120 deg ... is not used."

THE MM MODEL OF THE SMALL MODEL.  MCPB.py writes no topology for its small model.  Built here with tleap from
the SAME templates and parameter files as the dry bonded complex, so that every term that exists in the real
system has the same value:
  - ligand SN1 and zinc ZN1: the mol2 files of MCPB.py step 3, unchanged;
  - each histidine: its mol2 file cut to the side chain, CA and HA; the backbone atoms N and C are replaced
    by the capping hydrogens of the quantum-chemical model (H1 at the place of N, H3 at the place of C; CA is
    the atom CH3 of that model, HA its H2), H, O are dropped.  Atom types: CA keeps CX, HA keeps H1, the two
    capping hydrogens get H1 (ff14SB has every bond and angle term for H1 on CX);
  - the four bonds to the zinc; parameters: leaprc.protein.ff14SB, leaprc.gaff2, the two ligand frcmod
    files, CA2_mcpbpy.frcmod, then the --extra files.
The atom order of the topology is the order of the quantum-chemical input (checked).
Checked against CA2_dry.prmtop when present: charge, type and Lennard-Jones parameters of every atom the two
have in common, and force constant and equilibrium value of every bond, angle and torsion among common atoms.

WHAT CANCELS.  Only energies RELATIVE to point X are used, zinc and imidazoles never move.  Every term among
the 46 fixed atoms is therefore the same number in all points: all bonded terms of the imidazoles and of
the caps (which is why the bonded parameters of the capping hydrogens do not matter), the Zn-N(His) bonds,
the N(His)-Zn-N(His) angles, all non-bonded terms among fixed atoms.  WHAT DOES NOT CANCEL: the ligand's own
terms where its geometry changes, the terms at N2, Zn-N2 and the N(His)-Zn-N2 angles, and the non-bonded
terms between the ligand and the fixed atoms, the caps among them.  The cap charges enter there and only
there.  The RESP charges come from the large model, which has the backbone; the small model has a methyl group
in its place, so the charge of the removed atoms N, H, C, O (about -0.11 e per histidine) must go somewhere.
Three ways are evaluated, every one keeps the total charge of the residue and of the model (+1):
  methyl  the three hydrogens of the methyl group equal, each with the charge of HA; the carbon takes the
          rest of q(N) + q(H) + q(CA) + q(HA) + q(C) + q(O).  A methyl group as the quantum-chemical model
          has one (carbon about -0.3, hydrogens about +0.1).  USED.
  group   each capping hydrogen carries the charge of the group it replaces: H1 <- q(N) + q(H),
          H3 <- q(C) + q(O); CA and HA keep their charges.  Nearest to the charges of the real system.
  spread  the capping hydrogens neutral, the charge of the removed atoms spread evenly over the 15 atoms.
          An extreme, to bracket.
The result holds the relative energies of all three and their largest differences (`cap_test`); md/qm/gr_fit.py
gives the gate numbers for each of them.

ENERGIES.  OpenMM (Reference platform, double precision), no cut-off, vacuum, kcal/mol.  Per point: total,
bond, angle, torsion and non-bonded parts, the "fit terms" (the three angle terms at N2 and all terms on
H-N2-S-C4) evaluated from the parameters in the topology, and `rest` = total - fit terms, which is also
computed by OpenMM with these terms switched off (the two must agree; checked).

SPLIT BY ENERGY PART (review of 6 Oct 2026, finding 3; md/PLAN.md, entry of 8:56 AM IST, point 2).  Every
point also carries `split`: bond; angle_N2 (the three angle terms at N2) and angle_other; torsion_N2S (every
torsion term about the bond N2-S: the H-N2-S-C4 series, H-N2-S-O, and the terms through the zinc, which are
zero) and torsion_other; nb_ligand (Lennard-Jones and Coulomb among the 18 ligand atoms, 1-4 pairs
included); lj_ligand_frame (with its r^-12 and r^-6 parts) and coulomb_ligand_frame (ligand against zinc and
imidazoles, the scaled 1-4 pairs through the zinc included); nb_frame (among the 46 fixed atoms: the same
number in every point).  The non-bonded parts are summed here from the particle and exception parameters of
OpenMM's own system; their sum must equal OpenMM's non-bonded energy, and the nine parts the total (checked).

BINDING CHECK (md/PLAN.md, 3:02 AM, point 6, and 8:56 AM, point 4).  The small model is minimised (a) with zinc
and imidazoles fixed, (b) with all atoms free; C4-S-N2, S-N2 and Zn-N2-S at the minimum (b) are judged against
the quantum-chemical optimum of the same model (tolerances 2 deg, 0.02 A, 3 deg); the values of 6RL9 (108.8 deg,
1.596 A, 117.6 deg) are given beside and do not decide (crystal against vacuum is not like for like).
WHAT A FAILURE MEANS (md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 3): the result is reported; the ligand terms are
NOT re-tuned to this optimum, because the optimum itself (B3LYP/6-31G*, vacuum, no dispersion term) is compared
with the 1.0 A structure in the same file (`quantum_chemistry_against_crystal`, `optimum_against_crystal`) and
misses it on S-N2 by more than the tolerance.  The like-for-like test of the ligand geometry is the validation
run in the protein against 6RL9 (gate G-V).
    gr_mm.py binding --mm MM.json --out OUT.json [--fchk FILE ...] [--opt-geometry FILE.npy ...]
  The first --fchk file that exists is read (md/qm/mcpb/CA2_small_opt.fchk, written by task qm_write of the
  rented machine; on the Pi it arrives under /root/cqc-sfa-md/from_gpu/repo_md/qm/mcpb/).  Before that file
  exists, the geometry of the finished optimisation (opt_final.npy of md/qm/run_qm.py, Angstrom) is taken if
  given and present, and the verdict is marked provisional.  With neither: verdict "pending", exit code 2.

RELAXED REFERENCES OF PART B (md/PLAN.md, 8:56 AM, point 2).
    gr_mm.py relax --mm MM.json --out OUT.json --geom-out DIR [--reuse] [--qm-geom NAME=FILE.xyz ...]
  From X and from Rx the ligand is minimised with zinc and imidazoles fixed (the topology of the used cap
  scheme named in MM.json, the terms as built).  The signed volume of Zn, H and S about N2 is recorded at
  every iteration of the minimiser, so that a crossing to the other configuration cannot go unnoticed.  The
  result holds both minima with their split, their difference, and single points with split at every geometry
  given with --qm-geom (the quantum-chemical relaxed geometries, when they exist: qm_X, qm_Rx the converged
  ones, qm_X_last, qm_Rx_last the last step of a run that has not converged).
  A quantum-chemical geometry is accepted when its fixed atoms lie within FIXED_TOL (0.01 A) of those of X:
  geomeTRIC holds frozen atoms only approximately unless `enforce` is set (0.005 A was measured, 6 Oct 2026).
  The largest shift is written into the result, and the single point is given at the geometry as it is (the
  one the quantum-chemical energy belongs to) and with the fixed atoms put back exactly.
  --reuse: the force-field minimisations of OUT.json are kept when that file was made with the same model
  (same parameter files, cap scheme and points); only the single points are made again.
"""
import argparse
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
KCAL = 4.184
FIXED_TOL = 0.01                                                           # A: fixed atoms of a quantum-chemical geometry
SCHEMES = ('methyl', 'group', 'spread')
USED = 'methyl'
RES = ('HD1', 'HD2', 'HE1')
ZN_BONDS = (('1', 'NE2'), ('2', 'NE2'), ('3', 'ND1'))
CRYSTAL = dict(C4_S_N2=(108.8, 2.0), S_N2=(1.596, 0.02), Zn_N2_S=(117.6, 3.0))


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


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


def read_xyz(path):
    l = pathlib.Path(path).read_text().splitlines()
    n = int(l[0])
    sym = [t.split()[0] for t in l[2:2 + n]]
    xyz = np.array([[float(v) for v in t.split()[1:4]] for t in l[2:2 + n]])
    return sym, xyz


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


def build(params, work, scheme, points_dir, extra, n2h):
    """tleap topology of the small model; returns the paths and a record."""
    params, work = pathlib.Path(params), pathlib.Path(work) / scheme
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    rec = dict(scheme=scheme, inputs={})
    lines = (points_dir / 'geom' / 'X.pdb').read_text().splitlines()
    model = [(l[12:16].strip(), l[17:20].strip(), int(l[22:26])) for l in lines if l.startswith(('ATOM', 'HETATM'))]
    sym, xyz = read_xyz(points_dir / 'geom' / 'X.xyz')
    cap_q = {}
    for k, rn in enumerate(RES, 1):
        src_atoms, src_bonds = read_mol2(params / f'{rn}.mol2')
        rec['inputs'][f'{rn}.mol2'] = sha(params / f'{rn}.mol2')
        idx = [i for i, m in enumerate(model) if m[2] == k]
        names = [model[i][0] for i in idx]
        atoms, removed, ren, keep = capped_residue(rn, src_atoms, names, xyz[idx], scheme)
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
        write_mol2(work / f'{rn}.mol2', rn, atoms, bonds)
        cap_q[rn] = dict(removed_backbone_charge=round(removed, 6),
                         **{a['name']: round(a['q'], 6) for a in atoms if a['name'] in ('H1', 'CH3', 'H2', 'H3')},
                         residue_total=round(sum(a['q'] for a in atoms), 6))
    for f in ('ZN1.mol2', 'SN1.mol2', 'sfa_anion_gaff2.frcmod', 'sfa_anion_repair.frcmod', 'CA2_mcpbpy.frcmod'):
        shutil.copyfile(params / f, work / f)
        rec['inputs'][f] = sha(params / f)
    frcmods = ['sfa_anion_gaff2.frcmod', 'sfa_anion_repair.frcmod', 'CA2_mcpbpy.frcmod']
    if n2h == 'structure':
        k_line = [l for l in (params / 'CA2_mcpbpy.frcmod').read_text().splitlines() if l.startswith('M1-Y4-hn')]
        assert len(k_line) == 1, 'angle M1-Y4-hn not found in CA2_mcpbpy.frcmod'
        kval = float(k_line[0][8:].split()[0])
        pts = json.loads((points_dir / 'points.json').read_text())
        th = [p for p in pts['points'] if p['id'] == 'X'][0]['internal']['Zn_N2_H']
        (work / 'n2h_structure.frcmod').write_text(
            'Zn-N2-H at its value in start X instead of the 120 deg of MCPB.py step 2e (md/qm/gr_mm.py)\n'
            'MASS\n\nBOND\n\nANGLE\n'
            f'M1-Y4-hn   {kval:7.2f}   {th:9.3f}    value of start X; force constant of CA2_mcpbpy.frcmod\n'
            '\nDIHE\n\nIMPROPER\n\nNONBON\n\n')
        frcmods.append('n2h_structure.frcmod')
        rec['n2h'] = dict(mode='structure', k=kval, theta0=th, replaced_line=k_line[0].strip())
    else:
        rec['n2h'] = dict(mode='asis')
    for n, f in enumerate(extra):
        tgt = f'extra{n + 1}.frcmod'
        shutil.copyfile(f, work / tgt)
        frcmods.append(tgt)
        rec['inputs'][f'extra{n + 1}:{pathlib.Path(f).name}'] = sha(f)
    shutil.copyfile(points_dir / 'geom' / 'X.pdb', work / 'X.pdb')
    leap = ['source leaprc.protein.ff14SB', 'source leaprc.gaff2',
            'addAtomTypes {', '        { "M1"  "Zn" "sp3" }', '        { "Y1"  "N" "sp3" }',
            '        { "Y2"  "N" "sp3" }', '        { "Y3"  "N" "sp3" }', '        { "Y4"  "N" "sp3" }', '}']
    leap += [f'{rn} = loadmol2 {rn}.mol2' for rn in RES + ('ZN1', 'SN1')]
    leap += [f'loadamberparams {f}' for f in frcmods]
    leap += ['mol = loadpdb X.pdb']
    leap += [f'bond mol.{r}.{a} mol.4.ZN' for r, a in ZN_BONDS] + ['bond mol.4.ZN mol.5.N2']
    leap += ['saveamberparm mol small.prmtop small.inpcrd', 'check mol', 'charge mol', 'quit']
    (work / 'leap.in').write_text('\n'.join(leap) + '\n')
    env = dict(os.environ, AMBERHOME=AMBERHOME, PATH=f'{AMBERHOME}/bin:' + os.environ.get('PATH', ''))
    p = subprocess.run(['tleap', '-s', '-f', 'leap.in'], cwd=work, env=env, capture_output=True, text=True)
    (work / 'leap.log').write_text(p.stdout + p.stderr)
    log = p.stdout + p.stderr
    bad = [l for l in log.splitlines() if 'Could not find' in l or 'FATAL' in l or 'missing' in l.lower()
           and 'parameter' in l.lower()]
    if not (work / 'small.prmtop').exists() or bad:
        raise SystemExit(f'tleap failed for the small model ({scheme}): {bad[:5]}; see {work}/leap.log')
    rec.update(frcmods=frcmods, cap_charges=cap_q, tleap_warnings=[l for l in log.splitlines() if 'arning' in l][:12],
               prmtop=str(work / 'small.prmtop'))
    return work / 'small.prmtop', rec


def check_topology(prmtop, points_dir, params):
    """Order, charge, and identity of the common terms with the dry bonded complex."""
    import parmed
    top = parmed.load_file(str(prmtop))
    sym, _ = read_xyz(points_dir / 'geom' / 'X.xyz')
    el = [parmed.periodic_table.Element[a.atomic_number] for a in top.atoms]
    rec = dict(natoms=len(top.atoms), order_is_qm_order=bool([e.capitalize() for e in el] == sym),
               total_charge=round(sum(a.charge for a in top.atoms), 6),
               n_bonds=len(top.bonds), n_angles=len(top.angles), n_dihedrals=len(top.dihedrals))
    zn = [a for a in top.atoms if a.atomic_number == 30][0]
    rec['zinc_bonded_to'] = sorted(f'{b.residue.name}:{b.name}' for b in zn.bond_partners)
    dry = pathlib.Path(params) / 'CA2_dry.prmtop'
    if not dry.exists():
        rec['dry_complex'] = None
        return rec, top
    big = parmed.load_file(str(dry))
    site = {}
    for r in big.residues:
        if r.name in RES + ('ZN1', 'SN1'):
            for a in r.atoms:
                site[(r.name, a.name)] = a
    ren = {'CH3': 'CA', 'H2': 'HA'}
    caps = {'H1', 'H3'}
    pair, cap_atoms = {}, set()
    for a in top.atoms:
        if a.residue.name in RES and a.name in caps:
            cap_atoms.add(a.idx)
            continue
        nm = ren.get(a.name, a.name) if a.residue.name in RES else a.name
        pair[a.idx] = site[(a.residue.name, nm)]
    dq = max(abs(a.charge - pair[a.idx].charge) for a in top.atoms if a.idx in pair)
    types = sum(a.type != pair[a.idx].type or abs(a.rmin - pair[a.idx].rmin) > 1e-6 or
                abs(a.epsilon - pair[a.idx].epsilon) > 1e-6 for a in top.atoms if a.idx in pair)

    def key(atoms):
        k = tuple(pair[a.idx].idx for a in atoms)
        return min(k, k[::-1])
    bb = {key((b.atom1, b.atom2)): (b.type.k, b.type.req) for b in big.bonds
          if all(x.residue.name in RES + ('ZN1', 'SN1') for x in (b.atom1, b.atom2))} if False else None
    big_b = {min((b.atom1.idx, b.atom2.idx), (b.atom2.idx, b.atom1.idx)): (b.type.k, b.type.req) for b in big.bonds}
    big_a = {min((a.atom1.idx, a.atom2.idx, a.atom3.idx), (a.atom3.idx, a.atom2.idx, a.atom1.idx)):
             (a.type.k, a.type.theteq) for a in big.angles}
    big_d = {}
    for d in big.dihedrals:
        k = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
        big_d.setdefault((min(k, k[::-1]), d.improper), []).append(
            (round(d.type.phi_k, 6), d.type.per, round(d.type.phase, 4), round(d.type.scee, 4), round(d.type.scnb, 4)))
    n = dict(bonds=0, angles=0, dihedrals=0)
    miss, diff = [], []
    for b in top.bonds:
        if b.atom1.idx in pair and b.atom2.idx in pair:
            n['bonds'] += 1
            v = big_b.get(key((b.atom1, b.atom2)))
            if v is None:
                miss.append(f'bond {b.atom1.name}-{b.atom2.name}')
            elif abs(v[0] - b.type.k) > 1e-6 or abs(v[1] - b.type.req) > 1e-6:
                diff.append(f'bond {b.atom1.name}-{b.atom2.name}')
    for a in top.angles:
        at = (a.atom1, a.atom2, a.atom3)
        if all(x.idx in pair for x in at):
            n['angles'] += 1
            v = big_a.get(key(at))
            if v is None:
                miss.append('angle ' + '-'.join(x.name for x in at))
            elif abs(v[0] - a.type.k) > 1e-6 or abs(v[1] - a.type.theteq) > 1e-6:
                diff.append('angle ' + '-'.join(x.name for x in at) + f' {a.type.k} {a.type.theteq} / {v}')
    small_d = {}
    for d in top.dihedrals:
        at = (d.atom1, d.atom2, d.atom3, d.atom4)
        if all(x.idx in pair for x in at):
            small_d.setdefault((key(at), d.improper, '-'.join(x.name for x in at)), []).append(
                (round(d.type.phi_k, 6), d.type.per, round(d.type.phase, 4), round(d.type.scee, 4), round(d.type.scnb, 4)))
    for (k, imp, label), terms in small_d.items():
        n['dihedrals'] += 1
        v = big_d.get((k, imp))
        if v is None:
            # an improper is stored with its central atom third; tleap may order the outer atoms differently
            cand = [t for (kk, ii), t in big_d.items() if ii == imp and sorted(kk) == sorted(k)]
            v = cand[0] if cand else None
        if v is None:
            miss.append(f"{'improper' if imp else 'torsion'} {label}")
        elif sorted(v) != sorted(terms):
            diff.append(f"{'improper' if imp else 'torsion'} {label}: {sorted(terms)} / {sorted(v)}")
    rec['dry_complex'] = dict(file=str(dry), sha256=sha(dry), common_atoms=len(pair), cap_atoms=len(cap_atoms),
                              largest_charge_difference_common_atoms=float(dq), atoms_with_other_type_or_lj=int(types),
                              compared=n, not_found_in_dry_complex=miss, different_parameters=diff)
    return rec, top


class Engine:
    def __init__(self, prmtop):
        import openmm
        from openmm import app, unit
        self.mm, self.unit = openmm, unit
        self.top = app.AmberPrmtopFile(str(prmtop))
        self.names = [(a.residue.name, a.name) for a in self.top.topology.atoms()]
        self.idx = {n: k for k, n in enumerate(self.names)}
        self.zn, self.n2 = self.idx[('ZN1', 'ZN')], self.idx[('SN1', 'N2')]
        self.s, self.h, self.c4 = self.idx[('SN1', 'S')], self.idx[('SN1', 'HN21')], self.idx[('SN1', 'C4')]
        self.system = self._system()
        self.fit_terms = self._fit_terms(self.system)
        self.rest_system = self._system(zero=True)
        self.ctx = self._context(self.system)
        self.rest_ctx = self._context(self.rest_system)

    def _system(self, zero=False):
        from openmm import app
        s = self.top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
                                  removeCMMotion=False)
        grp = dict(HarmonicBondForce=0, HarmonicAngleForce=1, PeriodicTorsionForce=2, NonbondedForce=3)
        for f in s.getForces():
            f.setForceGroup(grp.get(type(f).__name__, 4))
            if zero and type(f).__name__ == 'HarmonicAngleForce':
                for k in range(f.getNumAngles()):
                    a, b, c, th, kk = f.getAngleParameters(k)
                    if b == self.n2:
                        f.setAngleParameters(k, a, b, c, th, 0.0)
            if zero and type(f).__name__ == 'PeriodicTorsionForce':
                for k in range(f.getNumTorsions()):
                    a, b, c, d, per, ph, kk = f.getTorsionParameters(k)
                    if (a, b, c, d) in ((self.h, self.n2, self.s, self.c4), (self.c4, self.s, self.n2, self.h)):
                        f.setTorsionParameters(k, a, b, c, d, per, ph, 0.0)
        return s

    def _fit_terms(self, s):
        u = self.unit
        ang, tor, other_tor = {}, [], []
        for f in s.getForces():
            if type(f).__name__ == 'HarmonicAngleForce':
                for k in range(f.getNumAngles()):
                    a, b, c, th, kk = f.getAngleParameters(k)
                    if b == self.n2:
                        ends = {a, c}
                        name = ('Zn_N2_S' if ends == {self.zn, self.s} else 'Zn_N2_H' if ends == {self.zn, self.h}
                                else 'S_N2_H' if ends == {self.s, self.h} else None)
                        assert name, 'unexpected angle at N2'
                        # OpenMM: E = k/2 (theta - theta0)^2, kJ/mol/rad^2;  Amber: E = K (theta - theta0)^2
                        ang[name] = dict(k=kk.value_in_unit(u.kilojoule_per_mole / u.radian ** 2) / 2 / KCAL,
                                         theta0=th.value_in_unit(u.degree))
            if type(f).__name__ == 'PeriodicTorsionForce':
                for k in range(f.getNumTorsions()):
                    a, b, c, d, per, ph, kk = f.getTorsionParameters(k)
                    t = dict(n=int(per), phase=ph.value_in_unit(u.degree),
                             k=kk.value_in_unit(u.kilojoule_per_mole) / KCAL)
                    if (a, b, c, d) in ((self.h, self.n2, self.s, self.c4), (self.c4, self.s, self.n2, self.h)):
                        tor.append(t)
                    elif {b, c} == {self.n2, self.s}:
                        other_tor.append(dict(atoms='-'.join(self.names[x][1] for x in (a, b, c, d)), **t))
        assert len(ang) == 3, f'three angle terms at N2 expected, found {sorted(ang)}'
        return dict(angles=ang, torsion_H_N2_S_C4=sorted(tor, key=lambda t: t['n']),
                    other_torsions_about_N2_S=sorted(other_tor, key=lambda t: (t['atoms'], t['n'])))

    def _context(self, system):
        integ = self.mm.VerletIntegrator(0.001)
        return self.mm.Context(system, integ, self.mm.Platform.getPlatformByName('Reference'))

    def _tables(self):
        """Pair tables of the non-bonded force and the lists of angle and torsion terms, from OpenMM's system."""
        u = self.unit
        n = self.system.getNumParticles()
        for f in self.system.getForces():
            name = type(f).__name__
            if name == 'NonbondedForce':
                q, sg, ep = np.zeros(n), np.zeros(n), np.zeros(n)
                for k in range(n):
                    a, b, c = f.getParticleParameters(k)
                    q[k], sg[k], ep[k] = (a.value_in_unit(u.elementary_charge), b.value_in_unit(u.nanometer),
                                          c.value_in_unit(u.kilojoule_per_mole))
                qq, sig, eps = np.outer(q, q), 0.5 * (sg[:, None] + sg[None, :]), np.sqrt(np.outer(ep, ep))
                for k in range(f.getNumExceptions()):
                    i, j, a, b, c = f.getExceptionParameters(k)
                    qq[i, j] = qq[j, i] = a.value_in_unit(u.elementary_charge ** 2)
                    sig[i, j] = sig[j, i] = b.value_in_unit(u.nanometer)
                    eps[i, j] = eps[j, i] = c.value_in_unit(u.kilojoule_per_mole)
                self.nb = (qq, sig, eps)
            elif name == 'HarmonicAngleForce':
                self.ang_n2 = []
                for k in range(f.getNumAngles()):
                    a, b, c, th, kk = f.getAngleParameters(k)
                    if b == self.n2:
                        self.ang_n2.append((a, b, c, th.value_in_unit(u.radian),
                                            kk.value_in_unit(u.kilojoule_per_mole / u.radian ** 2)))
            elif name == 'PeriodicTorsionForce':
                self.tor_n2s = []
                for k in range(f.getNumTorsions()):
                    a, b, c, d, per, ph, kk = f.getTorsionParameters(k)
                    if {b, c} == {self.n2, self.s}:
                        self.tor_n2s.append((a, b, c, d, int(per), ph.value_in_unit(u.radian),
                                             kk.value_in_unit(u.kilojoule_per_mole)))
        self.lig = np.array([r == 'SN1' for r, _ in self.names])

    def split(self, xyz_ang, e):
        """The parts named in the header, kcal/mol; `e` is the result of energy() at the same geometry."""
        if not hasattr(self, 'nb'):
            self._tables()
        x = np.asarray(xyz_ang) * 0.1
        n = len(x)
        iu = np.triu_indices(n, 1)
        r = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)[iu]
        qq, sig, eps = (t[iu] for t in self.nb)
        coul = 138.93545764438198 * qq / r / KCAL
        s6 = (sig / r) ** 6
        l12, l6 = 4 * eps * s6 * s6 / KCAL, -4 * eps * s6 / KCAL
        li, lj = self.lig[iu[0]], self.lig[iu[1]]
        ll, lf, ff = li & lj, li ^ lj, ~li & ~lj
        ang = 0.0
        for a, b, c, th0, kk in self.ang_n2:
            v, w = x[a] - x[b], x[c] - x[b]
            th = np.arccos(np.dot(v, w) / np.linalg.norm(v) / np.linalg.norm(w))
            ang += 0.5 * kk * (th - th0) ** 2 / KCAL
        tor = 0.0
        for a, b, c, d, per, ph, kk in self.tor_n2s:
            b1 = (x[c] - x[b]) / np.linalg.norm(x[c] - x[b])
            v = (x[a] - x[b]) - np.dot(x[a] - x[b], b1) * b1
            w = (x[d] - x[c]) - np.dot(x[d] - x[c], b1) * b1
            phi = np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))
            tor += kk * (1.0 + np.cos(per * phi - ph)) / KCAL
        sp = dict(bond=e['bond'], angle_N2=float(ang), angle_other=e['angle'] - float(ang), torsion_N2S=float(tor),
                  torsion_other=e['torsion'] - float(tor),
                  nb_ligand=float((coul + l12 + l6)[ll].sum()), lj_ligand_frame=float((l12 + l6)[lf].sum()),
                  coulomb_ligand_frame=float(coul[lf].sum()), nb_frame=float((coul + l12 + l6)[ff].sum()))
        extra = dict(lj_ligand_frame_r12=float(l12[lf].sum()), lj_ligand_frame_r6=float(l6[lf].sum()),
                     nonbonded_own_sum_minus_openmm=float((coul + l12 + l6).sum() - e['nonbonded']),
                     parts_minus_total=float(sum(sp.values()) - e['total']))
        return sp, extra

    def energy(self, xyz_ang):
        u = self.unit
        self.ctx.setPositions(np.asarray(xyz_ang) * 0.1)
        self.rest_ctx.setPositions(np.asarray(xyz_ang) * 0.1)
        e = {}
        for name, g in (('bond', 0), ('angle', 1), ('torsion', 2), ('nonbonded', 3)):
            e[name] = self.ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(
                u.kilocalorie_per_mole)
        e['total'] = self.ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        e['rest_openmm'] = self.rest_ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            u.kilocalorie_per_mole)
        return e


def signed_volume(xyz, zn, n2, h, s):
    """Triple product of the unit vectors N2->Zn, N2->H, N2->S: +0.72 in start X, -0.68 in Rx, 0 for a planar N2
    (the quantity of md/analysis/SPEC.md, section 9A, and of md/qm/gr_points.py)."""
    z, hh, u = (np.asarray(xyz[k]) - np.asarray(xyz[n2]) for k in (zn, h, s))
    z, hh, u = z / np.linalg.norm(z), hh / np.linalg.norm(hh), u / np.linalg.norm(u)
    return float(np.dot(z, np.cross(hh, u)))


def internal(xyz, eng):
    def ang(a, b, c):
        u, v = xyz[a] - xyz[b], xyz[c] - xyz[b]
        return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))

    def dih(a, b, c, d):
        b1 = (xyz[c] - xyz[b]) / np.linalg.norm(xyz[c] - xyz[b])
        v = (xyz[a] - xyz[b]) - np.dot(xyz[a] - xyz[b], b1) * b1
        w = (xyz[d] - xyz[c]) - np.dot(xyz[d] - xyz[c], b1) * b1
        return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))
    i = eng.idx
    return dict(signed_volume=signed_volume(xyz, eng.zn, eng.n2, eng.h, eng.s),
                Zn_N2_S=ang(eng.zn, eng.n2, eng.s), Zn_N2_H=ang(eng.zn, eng.n2, eng.h), S_N2_H=ang(eng.s, eng.n2, eng.h),
                H_N2_S_C4=dih(eng.h, eng.n2, eng.s, eng.c4),
                H_N2_S_O1=dih(eng.h, eng.n2, eng.s, i[('SN1', 'O1')]), H_N2_S_O2=dih(eng.h, eng.n2, eng.s, i[('SN1', 'O2')]),
                C4_S_N2=ang(eng.c4, eng.s, eng.n2), S_N2=float(np.linalg.norm(xyz[eng.s] - xyz[eng.n2])),
                Zn_N2=float(np.linalg.norm(xyz[eng.zn] - xyz[eng.n2])),
                Zn_N2_S_C4=dih(eng.zn, eng.n2, eng.s, eng.c4))


def fit_term_energy(ic, terms):
    e = 0.0
    for name, p in terms['angles'].items():
        e += p['k'] * np.radians(ic[name] - p['theta0']) ** 2
    for t in terms['torsion_H_N2_S_C4']:
        e += t['k'] * (1.0 + np.cos(np.radians(t['n'] * ic['H_N2_S_C4'] - t['phase'])))
    return float(e)


def evaluate(prmtop, points_dir, order):
    eng = Engine(prmtop)
    out, worst = {}, 0.0
    for pid in order:
        _, xyz = read_xyz(points_dir / 'geom' / f'{pid}.xyz')
        e = eng.energy(xyz)
        ic = internal(xyz, eng)
        e['fit_terms'] = fit_term_energy(ic, eng.fit_terms)
        e['rest'] = e['total'] - e['fit_terms']
        worst = max(worst, abs(e['rest'] - e['rest_openmm']))
        sp, extra = eng.split(xyz, e)
        assert abs(extra['nonbonded_own_sum_minus_openmm']) < 1e-5 and abs(extra['parts_minus_total']) < 1e-5, (pid, extra)
        out[pid] = dict(energy={k: round(v, 6) for k, v in e.items()}, internal={k: round(v, 5) for k, v in ic.items()},
                        split={k: round(v, 6) for k, v in sp.items()}, split_extra={k: round(v, 8) for k, v in extra.items()})
    return eng, out, worst


def minimise(eng, xyz0, frozen):
    """Local minimum from the geometry of X (Reference platform; deterministic)."""
    mm, u = eng.mm, eng.unit
    s = eng._system()
    if frozen:
        for k in frozen:
            s.setParticleMass(k, 0.0)
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(np.asarray(xyz0) * 0.1)
    e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    vols = [signed_volume(np.asarray(xyz0), eng.zn, eng.n2, eng.h, eng.s)]

    class Watch(mm.MinimizationReporter):
        def report(self, iteration, x, grad, args):
            c = np.array(x).reshape(-1, 3)
            vols.append(signed_volume(c, eng.zn, eng.n2, eng.h, eng.s))
            return False
    watch = Watch()
    mm.LocalEnergyMinimizer.minimize(ctx, 1e-4, 20000, watch)
    st = ctx.getState(getEnergy=True, getPositions=True, getForces=True)
    x = st.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    f = st.getForces(asNumpy=True).value_in_unit(u.kilocalorie_per_mole / u.angstrom)
    free = [k for k in range(len(x)) if not frozen or k not in frozen]
    ic = internal(np.array(x), eng)
    zn = np.array(x)[eng.zn]
    his = [np.linalg.norm(np.array(x)[eng.idx[(r, a)]] - zn) for r, a in (('HD1', 'NE2'), ('HD2', 'NE2'), ('HE1', 'ND1'))]
    res = dict(energy_start=round(e0, 4),
               energy_min=round(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), 4),
               rms_force_free_atoms=round(float(np.sqrt((np.array(f)[free] ** 2).sum(1).mean())), 5),
               rmsd_ligand_heavy_from_X=None, Zn_N_His=[round(float(t), 4) for t in his],
               internal={k: round(v, 4) for k, v in ic.items()})
    lig_heavy = [k for k, (r, a) in enumerate(eng.names) if r == 'SN1' and not a.startswith('H')]
    res['rmsd_ligand_heavy_from_X'] = round(float(np.sqrt(((np.array(x)[lig_heavy] - np.asarray(xyz0)[lig_heavy]) ** 2
                                                           ).sum(1).mean())), 4)
    res['angle_sum_N2'] = round(ic['Zn_N2_S'] + ic['Zn_N2_H'] + ic['S_N2_H'], 3)
    res['signed_volume_N2'] = dict(start=round(vols[0], 4), end=round(ic['signed_volume'], 4), n_iterations=len(vols) - 1,
                                   smallest=round(min(vols), 4), largest=round(max(vols), 4),
                                   changed_sign=bool(min(vols) * max(vols) <= 0))
    chk = {}
    for k, (ref, tol) in CRYSTAL.items():
        chk[k] = dict(mm=round(ic[k], 4), crystal=ref, tolerance=tol, passed=bool(abs(ic[k] - ref) <= tol))
    res['binding_check'] = chk
    res['binding_check_passed'] = bool(all(c['passed'] for c in chk.values()))
    return res, np.array(x)


def run(args):
    points_dir = pathlib.Path(args.points).resolve().parent
    pts = json.loads(pathlib.Path(args.points).read_text())
    order = pts['order']
    work = pathlib.Path(args.work)
    res = dict(label=args.label, params=str(pathlib.Path(args.params).resolve()),
               points_file=str(pathlib.Path(args.points).resolve()),
               points_sha256=sha(args.points), schemes={}, used_scheme=USED,
               note='energies in kcal/mol; relative energies are relative to point X')
    energies = {}
    for scheme in SCHEMES:
        prmtop, rec = build(args.params, work, scheme, points_dir, args.extra or [], args.n2h)
        chk, _ = check_topology(prmtop, points_dir, args.params)
        eng, ev, worst = evaluate(prmtop, points_dir, order)
        rec.update(topology_check=chk, rest_two_ways_largest_difference=worst)
        if scheme == USED:
            res['fit_terms_in_topology'] = eng.fit_terms
            res['points'] = ev
            res['n2h'] = rec['n2h']
            res['frcmods'] = rec['frcmods']
            res['inputs_sha256'] = rec['inputs']
            if not args.no_min:
                _, x0 = read_xyz(points_dir / 'geom' / 'X.xyz')
                fixed = [k for k, (r, a) in enumerate(eng.names) if r != 'SN1']
                m1, x1 = minimise(eng, x0, fixed)
                m2, x2 = minimise(eng, x0, None)
                res['minimum_ligand_free'] = m1
                res['minimum_all_free'] = m2
                res['binding_check'] = dict(
                    rule='AGAINST THE CRYSTAL VALUES ONLY, which do not decide (md/PLAN.md, 6 Oct 2026, 8:56 AM IST, '
                         'point 4): at the force-field minimum of the small model C4-S-N2 within 2 deg, S-N2 within '
                         '0.02 A, Zn-N2-S within 3 deg of 6RL9.  The verdict is that of `gr_mm.py binding`, against '
                         'the quantum-chemical optimum of the small model',
                    verdict_on='minimum_all_free', passed=m2['binding_check_passed'],
                    passed_with_zinc_and_imidazoles_fixed=m1['binding_check_passed'])
                sym, _ = read_xyz(points_dir / 'geom' / 'X.xyz')
                for tag, x in (('min_ligand_free', x1), ('min_all_free', x2)):
                    (work / f'{tag}.xyz').write_text(
                        '64\nMM minimum of the small model, ' + tag + ', ' + args.label + '\n' +
                        '\n'.join(f'{s:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}' for s, v in zip(sym, x)) + '\n')
        rel = {p: ev[p]['energy']['total'] - ev['X']['energy']['total'] for p in order}
        energies[scheme] = rel
        rec['relative_energy'] = {p: round(v, 4) for p, v in rel.items()}
        rec['total'] = {p: ev[p]['energy']['total'] for p in order}
        rec['rest'] = {p: ev[p]['energy']['rest'] for p in order}
        res['schemes'][scheme] = rec
    cap = {}
    for a in range(len(SCHEMES)):
        for b in range(a + 1, len(SCHEMES)):
            d = {p: energies[SCHEMES[a]][p] - energies[SCHEMES[b]][p] for p in order}
            worst_p = max(d, key=lambda p: abs(d[p]))
            cap[f'{SCHEMES[a]}_minus_{SCHEMES[b]}'] = dict(
                largest_abs=round(abs(d[worst_p]), 4), at=worst_p,
                rms=round(float(np.sqrt(np.mean([v ** 2 for v in d.values()]))), 4),
                per_point={p: round(v, 4) for p, v in d.items()})
    res['cap_test'] = cap
    res['cap_test_largest_abs_two_reasonable'] = cap['methyl_minus_group']['largest_abs']
    res['cap_test_largest_abs_all'] = max(c['largest_abs'] for c in cap.values())
    import openmm
    import parmed
    res['versions'] = dict(openmm=openmm.__version__, parmed=parmed.__version__, numpy=np.__version__,
                           python=sys.version.split()[0], amberhome=AMBERHOME)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix('.tmp')
    tmp.write_text(json.dumps(res, indent=1) + '\n')
    os.replace(tmp, out)
    c = res['schemes'][USED]['topology_check']
    print(f"{args.label}: {c['natoms']} atoms, order as QM {c['order_is_qm_order']}, charge {c['total_charge']}, "
          f"zinc bonded to {c['zinc_bonded_to']}")
    if c['dry_complex']:
        d = c['dry_complex']
        print(f"  against the dry complex: {d['compared']}, not found {len(d['not_found_in_dry_complex'])}, "
              f"different {len(d['different_parameters'])}, charge difference {d['largest_charge_difference_common_atoms']:.1e}")
    print(f"  cap test: largest difference of a relative energy {res['cap_test_largest_abs_two_reasonable']} "
          f"(methyl against group), {res['cap_test_largest_abs_all']} (all three) kcal/mol")
    if not args.no_min:
        print('  binding check (all free):', res['minimum_all_free']['binding_check'])
    return res


# ------------------------------------------------------------------------------- relaxed references of part B
SPLIT_KEYS = ('bond', 'angle_N2', 'angle_other', 'torsion_N2S', 'torsion_other', 'nb_ligand', 'lj_ligand_frame',
              'coulomb_ligand_frame', 'nb_frame')


def single(eng, xyz):
    e = eng.energy(xyz)
    sp, extra = eng.split(xyz, e)
    assert abs(extra['nonbonded_own_sum_minus_openmm']) < 1e-5 and abs(extra['parts_minus_total']) < 1e-5, extra
    return dict(total=round(e['total'], 6), split={k: round(v, 6) for k, v in sp.items()},
                lj_ligand_frame_r12=round(extra['lj_ligand_frame_r12'], 6), lj_ligand_frame_r6=round(extra['lj_ligand_frame_r6'], 6),
                internal={k: round(v, 5) for k, v in internal(np.asarray(xyz), eng).items()})


def difference(a, b):
    """b minus a: total and every part."""
    return dict(total=round(b['total'] - a['total'], 4),
                split={k: round(b['split'][k] - a['split'][k], 4) for k in SPLIT_KEYS})


def careful_descent(eng, xyz0, frozen, max_step=0.02, max_iter=60000, f_tol=0.005):
    """Conjugate-gradient descent (Polak-Ribiere, restarted whenever a step is refused) in which no atom moves
    more than `max_step` A in a step and a step is taken only if the energy falls.  Unlike the L-BFGS minimiser
    it cannot step over a barrier: if the signed volume at N2 changes sign on this path, the start does not lie
    in the basin of a minimum of its own configuration."""
    mm, u = eng.mm, eng.unit
    s = eng._system()
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    free = np.array([k not in set(frozen) for k in range(len(xyz0))])

    def ef(x):
        ctx.setPositions(x * 0.1)
        st = ctx.getState(getEnergy=True, getForces=True)
        f = np.array(st.getForces(asNumpy=True).value_in_unit(u.kilocalorie_per_mole / u.angstrom))
        f[~free] = 0.0
        return st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), f
    x = np.array(xyz0, dtype=float)
    e, f = ef(x)
    d = f.copy()
    e0, step, v0 = e, max_step, signed_volume(x, eng.zn, eng.n2, eng.h, eng.s)
    vmin = vmax = v0
    crossed, n_acc, why, largest_move = None, 0, 'iteration limit', 0.0
    for it in range(1, max_iter + 1):
        if float(np.sqrt((f[free] ** 2).sum(1).mean())) < f_tol:
            why = f'rms force below {f_tol} kcal/mol/A'
            break
        if float((d * f).sum()) <= 0:
            d = f.copy()
        dmax = float(np.linalg.norm(d, axis=1).max())
        x2 = x + d * (step / dmax)
        e2, f2 = ef(x2)
        if e2 < e:
            beta = max(0.0, float((f2 * (f2 - f)).sum() / (f * f).sum()))
            largest_move = max(largest_move, step)
            x, e, d, f = x2, e2, f2 + beta * d, f2
            n_acc += 1
            step = min(max_step, step * 1.3)
            v = signed_volume(x, eng.zn, eng.n2, eng.h, eng.s)
            vmin, vmax = min(vmin, v), max(vmax, v)
            if crossed is None and v * v0 <= 0:
                crossed = dict(iteration=it, energy_below_start=round(e0 - e, 4))
        else:
            if np.allclose(d, f) and step < 1e-9:
                why = 'no lower energy along the force'
                break
            step = step * 0.5 if np.allclose(d, f) else step
            d = f.copy()
    lig_heavy = [k for k, (r, a) in enumerate(eng.names) if r == 'SN1' and not a.startswith('H')]
    return dict(method=f'conjugate-gradient descent, largest atomic move {max_step} A per step, a step only if the energy '
                       f'falls, at most {max_iter} steps',
                iterations=it, accepted_steps=n_acc, ended=why, energy_start=round(e0, 4), energy_end=round(e, 4),
                rms_force_end=round(float(np.sqrt((f[free] ** 2).sum(1).mean())), 4), largest_atomic_move_in_a_step_A=round(largest_move, 5),
                signed_volume_N2=dict(start=round(v0, 4), end=round(signed_volume(x, eng.zn, eng.n2, eng.h, eng.s), 4),
                                      smallest=round(vmin, 4), largest=round(vmax, 4), changed_sign=bool(vmin * vmax <= 0)),
                crossed=crossed,
                rmsd_ligand_heavy_from_start_A=round(float(np.sqrt(((x[lig_heavy] - np.asarray(xyz0)[lig_heavy]) ** 2).sum(1).mean())), 4))


def relax(args):
    mm = json.loads(pathlib.Path(args.mm).read_text())
    points_dir = pathlib.Path(mm['points_file']).parent
    prmtop = pathlib.Path(mm['schemes'][mm['used_scheme']]['prmtop'])
    eng = Engine(prmtop)
    fixed = [k for k, (r, a) in enumerate(eng.names) if r != 'SN1']
    lig_heavy = [k for k, (r, a) in enumerate(eng.names) if r == 'SN1' and not a.startswith('H')]
    sym, _ = read_xyz(points_dir / 'geom' / 'X.xyz')
    gout = pathlib.Path(args.geom_out)
    gout.mkdir(parents=True, exist_ok=True)
    import openmm
    res = dict(what='part B of gate G-R, relaxed references in the force field: the ligand minimised with zinc and '
                    'imidazoles fixed, from X and from Rx (md/qm/gr_mm.py relax); kcal/mol, vacuum, no cut-off',
               mm=str(args.mm), mm_label=mm['label'], cap_scheme=mm['used_scheme'], prmtop=str(prmtop),
               prmtop_sha256=sha(prmtop), fit_terms_in_topology=eng.fit_terms,
               minimiser='OpenMM LocalEnergyMinimizer (L-BFGS), Reference platform, tolerance 1e-4 kJ/mol/nm, '
                         'at most 20000 iterations; fixed atoms have mass zero',
               openmm=openmm.__version__, starts={},
               mm_inputs=dict(inputs_sha256=mm.get('inputs_sha256'), frcmods=mm.get('frcmods'), n2h=mm.get('n2h'),
                              cap_scheme=mm['used_scheme'], points_sha256=mm.get('points_sha256')))
    old = None
    if args.reuse and pathlib.Path(args.out).exists():
        try:
            old = json.loads(pathlib.Path(args.out).read_text())
        except ValueError:
            old = None
        ok = bool(old and old.get('mm_inputs') == res['mm_inputs'] and 'fixed_atom_tolerance_A' in old and
                  set(old.get('starts') or {}) == {'X', 'Rx'} and
                  all((gout / f'mm_{s}.xyz').exists() for s in ('X', 'Rx')))
        old = old if ok else None
    for start in ('X', 'Rx'):
        if old:                                                            # same model: the minimisations are not repeated
            res['starts'][start] = old['starts'][start]
            continue
        _, x0 = read_xyz(points_dir / 'geom' / f'{start}.xyz')
        m, x1 = minimise(eng, x0, fixed)
        rec = dict(start_geometry=str(points_dir / 'geom' / f'{start}.xyz'), minimisation=m, at_start=single(eng, x0),
                   at_minimum=single(eng, x1), relaxation=None,
                   largest_move_of_a_fixed_atom_A=float(np.abs(x1[fixed] - x0[fixed]).max()),
                   rmsd_ligand_heavy_from_start_A=round(float(np.sqrt(((x1[lig_heavy] - x0[lig_heavy]) ** 2).sum(1).mean())), 4))
        rec['relaxation'] = difference(rec['at_start'], rec['at_minimum'])
        (gout / f'mm_{start}.xyz').write_text(
            f'64\nforce-field minimum of the ligand from {start}, zinc and imidazoles fixed; {mm["label"]}\n' +
            '\n'.join(f'{a:2s} {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}' for a, v in zip(sym, x1)) + '\n')
        rec['geometry'] = str(gout / f'mm_{start}.xyz')
        rec['careful_descent'] = careful_descent(eng, x0, fixed)
        res['starts'][start] = rec
    res['minimisations'] = (dict(reused=True, from_file_written=old.get('written'), note='same parameter files, cap scheme and '
                                 'points as the model of now; the Reference platform is deterministic')
                            if old else dict(reused=False))
    res['written'] = time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))
    a, b = res['starts']['X'], res['starts']['Rx']
    res['Rx_minus_X'] = dict(unrelaxed=difference(a['at_start'], b['at_start']),
                             relaxed=difference(a['at_minimum'], b['at_minimum']))
    res['configurations_kept'] = bool(not a['minimisation']['signed_volume_N2']['changed_sign'] and
                                      not b['minimisation']['signed_volume_N2']['changed_sign'] and
                                      a['at_minimum']['internal']['signed_volume'] > 0.3 and
                                      b['at_minimum']['internal']['signed_volume'] < -0.3)
    res['careful_descent'] = {t: dict(crossed=res['starts'][t]['careful_descent']['crossed'],
                                      signed_volume_N2=res['starts'][t]['careful_descent']['signed_volume_N2'],
                                      energy_drop=round(res['starts'][t]['careful_descent']['energy_start'] -
                                                        res['starts'][t]['careful_descent']['energy_end'], 4),
                                      ended=res['starts'][t]['careful_descent']['ended']) for t in res['starts']}
    res['at_other_geometries'] = {}
    res['fixed_atom_tolerance_A'] = FIXED_TOL
    x_ref = read_xyz(points_dir / 'geom' / 'X.xyz')[1]
    mm_min = {s: read_xyz(gout / f'mm_{s}.xyz')[1] for s in ('X', 'Rx') if (gout / f'mm_{s}.xyz').exists()}
    for item in args.qm_geom or []:
        name, f = item.split('=', 1)
        if not pathlib.Path(f).exists():
            res['at_other_geometries'][name] = dict(file=f, exists=False)
            continue
        sy, x = read_xyz(f)
        if sy != sym:
            raise SystemExit(f'{f}: not the 64 atoms of the small model in its order')
        shift = float(np.abs(x[fixed] - x_ref[fixed]).max())
        if shift > FIXED_TOL:
            raise SystemExit(f'{f}: a fixed atom lies {shift:.4f} A from its place in X, more than {FIXED_TOL} A: not accepted')
        rec = dict(file=f, exists=True, sha256=sha(f), comment=pathlib.Path(f).read_text().splitlines()[1],
                   largest_shift_of_a_fixed_atom_A=round(shift, 6),
                   rms_shift_of_the_fixed_atoms_A=round(float(np.sqrt(((x[fixed] - x_ref[fixed]) ** 2).sum(1).mean())), 6),
                   **single(eng, x))
        xe = x.copy()
        xe[fixed] = x_ref[fixed]
        se = single(eng, xe)
        rec['with_the_fixed_atoms_put_back'] = dict(total=se['total'], split=se['split'],
                                                    total_minus_as_it_is=round(se['total'] - rec['total'], 6))
        rec['rmsd_ligand_heavy_from_X_start_A'] = round(float(np.sqrt(((x[lig_heavy] - x_ref[lig_heavy]) ** 2).sum(1).mean())), 4)
        for s, xm in mm_min.items():
            rec[f'rmsd_ligand_heavy_from_the_force_field_minimum_from_{s}_A'] = round(
                float(np.sqrt(((x[lig_heavy] - xm[lig_heavy]) ** 2).sum(1).mean())), 4)
            rec[f'above_the_force_field_minimum_from_{s}'] = round(rec['total'] - res['starts'][s]['at_minimum']['total'], 4)
        res['at_other_geometries'][name] = rec
    g = res['at_other_geometries']
    if g.get('qm_X', {}).get('exists') and g.get('qm_Rx', {}).get('exists'):
        res['Rx_minus_X']['force_field_at_the_quantum_chemical_minima'] = difference(g['qm_X'], g['qm_Rx'])
        res['Rx_minus_X']['force_field_at_the_quantum_chemical_minima_fixed_atoms_put_back'] = difference(
            g['qm_X']['with_the_fixed_atoms_put_back'], g['qm_Rx']['with_the_fixed_atoms_put_back'])
    # the two quantum-chemical geometries as they are now: the converged one of each start, else its last step
    pick = {s: next((k for k in (f'qm_{s}', f'qm_{s}_last') if g.get(k, {}).get('exists')), None) for s in ('X', 'Rx')}
    if pick['X'] and pick['Rx']:
        conv = {s: not pick[s].endswith('_last') for s in pick}
        xa, xb = (read_xyz(g[pick[s]]['file'])[1] for s in ('X', 'Rx'))
        res['Rx_minus_X']['force_field_at_the_quantum_chemical_geometries'] = dict(
            geometries=pick, converged=conv, both_converged=bool(all(conv.values())),
            note=('both quantum-chemical relaxations have converged' if all(conv.values()) else
                  'NOT both converged: a geometry named *_last is the last step of a run that is still going; not a result'),
            rmsd_ligand_heavy_between_the_two_A=round(float(np.sqrt(((xb[lig_heavy] - xa[lig_heavy]) ** 2).sum(1).mean())), 4),
            **difference(g[pick['X']], g[pick['Rx']]))
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix('.tmp')
    tmp.write_text(json.dumps(res, indent=1) + '\n')
    os.replace(tmp, out)
    for start in ('X', 'Rx'):
        r = res['starts'][start]
        v = r['minimisation']['signed_volume_N2']
        print(f"relaxed from {start}: E {r['at_start']['total']:.3f} -> {r['at_minimum']['total']:.3f} kcal/mol, ligand heavy atoms moved "
              f"{r['rmsd_ligand_heavy_from_start_A']} A, signed volume at N2 {v['start']} -> {v['end']} "
              f"(range {v['smallest']} to {v['largest']}, {v['n_iterations']} iterations), rms force {r['minimisation']['rms_force_free_atoms']}")
        c = r['careful_descent']
        print(f"   careful descent from {start}: {c['accepted_steps']} steps ({c['ended']}), E {c['energy_start']} -> {c['energy_end']}, "
              f"volume {c['signed_volume_N2']['start']} -> {c['signed_volume_N2']['end']}, crossed: {c['crossed']}")
    d = res['Rx_minus_X']
    print(f"Rx minus X: unrelaxed {d['unrelaxed']['total']}, relaxed {d['relaxed']['total']} kcal/mol; split of the relaxed "
          f"difference: {d['relaxed']['split']}; configurations kept: {res['configurations_kept']}")
    return res


# --------------------------------------------------------------------------------------------- the binding check
BOHR = 0.52917721092                                                       # as md/qm/gaucom.py, which wrote the fchk
TOL = dict(C4_S_N2=2.0, S_N2=0.02, Zn_N2_S=3.0)
UNIT = dict(C4_S_N2='deg', S_N2='A', Zn_N2_S='deg')


def read_fchk_geometry(path):
    lines = pathlib.Path(path).read_text().splitlines()

    def arr(label, cast):
        for k, l in enumerate(lines):
            if l.startswith(label):
                n = int(l.split('N=')[1])
                vals, j = [], k + 1
                while len(vals) < n:
                    vals += lines[j].split()
                    j += 1
                return [cast(v) for v in vals[:n]]
        raise KeyError(f'{label!r} not in {path}')
    return arr('Atomic numbers', int), np.array(arr('Current cartesian coordinates', float)).reshape(-1, 3) * BOHR


def site_geometry(xyz, names):
    i = {n: k for k, n in enumerate(names)}
    zn, n2, s, c4, h = i[('ZN1', 'ZN')], i[('SN1', 'N2')], i[('SN1', 'S')], i[('SN1', 'C4')], i[('SN1', 'HN21')]

    def ang(a, b, c):
        u, v = xyz[a] - xyz[b], xyz[c] - xyz[b]
        return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))
    his = [float(np.linalg.norm(xyz[i[k]] - xyz[zn])) for k in (('HD1', 'NE2'), ('HD2', 'NE2'), ('HE1', 'ND1'))]
    return dict(C4_S_N2=round(ang(c4, s, n2), 4), S_N2=round(float(np.linalg.norm(xyz[s] - xyz[n2])), 5),
                Zn_N2_S=round(ang(zn, n2, s), 4), Zn_N2=round(float(np.linalg.norm(xyz[zn] - xyz[n2])), 4),
                Zn_N_His=[round(v, 4) for v in his], Zn_N2_H=round(ang(zn, n2, h), 3), S_N2_H=round(ang(s, n2, h), 3),
                signed_volume=round(signed_volume(xyz, zn, n2, h, s), 4))


def binding(args):
    mm = json.loads(pathlib.Path(args.mm).read_text())
    points_dir = pathlib.Path(mm['points_file']).parent
    pdb = [l for l in (points_dir / 'geom' / 'X.pdb').read_text().splitlines() if l.startswith(('ATOM', 'HETATM'))]
    names = [(l[17:20].strip(), l[12:16].strip()) for l in pdb]
    sym, x_start = read_xyz(points_dir / 'geom' / 'X.xyz')
    zs = dict(H=1, C=6, N=7, O=8, S=16, Zn=30)
    res = dict(
        what='binding check of the ligand terms at sulfur and of Zn-N2-S (md/qm/gr_mm.py binding)',
        rule='md/PLAN.md, 6 Oct 2026, 3:02 AM IST, point 6, and 8:56 AM IST, point 4: at the force-field minimum of the '
             'small model (all atoms free, vacuum, final charges) C4-S-N2 within 2 deg, S-N2 within 0.02 A and Zn-N2-S '
             'within 3 deg of the quantum-chemical optimum of the same model (B3LYP/6-31G*, task qm_opt).  The values '
             'of 6RL9 are given beside and do not decide',
        consequence='md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 3 (which replaces the consequence named at 8:56 AM, a re-tuning of '
                    'the ligand terms at sulfur): the result is reported and the ligand terms are NOT re-tuned to this optimum, '
                    'because the reference itself, B3LYP/6-31G* in vacuum without a dispersion term, is outside the same '
                    'tolerances against the 1.0 A structure (`optimum_against_crystal`).  The like-for-like test of the ligand '
                    'geometry is the validation run in the protein against 6RL9 (gate G-V: C4-S-N2 within 3 deg, tau_fold '
                    'within 20 deg, ligand RMSD), which stands as defined',
        mm=str(args.mm), mm_label=mm['label'], tolerances=TOL, units=UNIT,
        written=time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800)),
        force_field=dict(minimum_all_free={k: mm['minimum_all_free']['internal'][k] for k in TOL},
                         minimum_with_zinc_and_imidazoles_fixed={k: mm['minimum_ligand_free']['internal'][k] for k in TOL},
                         ligand_heavy_rmsd_from_X_all_free=mm['minimum_all_free']['rmsd_ligand_heavy_from_X']),
        crystal_6RL9={k: v[0] for k, v in CRYSTAL.items()}, start_X=site_geometry(x_start, names))
    src = None
    for f in args.fchk or []:
        if pathlib.Path(f).exists():
            z, xyz = read_fchk_geometry(f)
            assert z == [zs[a] for a in sym], f'{f}: not the atoms of the small model in its order'
            src = dict(kind='fchk', file=str(f), sha256=sha(f), provisional=False)
            break
    if src is None:
        for f in args.opt_geometry or []:
            if pathlib.Path(f).exists():
                xyz = np.load(f)
                assert xyz.shape == (64, 3), f'{f}: shape {xyz.shape}'
                src = dict(kind='geometry of the finished optimisation (opt_final.npy); the fchk file is not written yet',
                           file=str(f), sha256=sha(f), provisional=True)
                break
    mmv = res['force_field']['minimum_all_free']
    res['against_crystal_not_the_verdict'] = {
        k: dict(force_field=mmv[k], crystal=CRYSTAL[k][0], difference=round(mmv[k] - CRYSTAL[k][0], 4), tolerance=TOL[k],
                within=bool(abs(mmv[k] - CRYSTAL[k][0]) <= TOL[k])) for k in TOL}
    if src is None:
        res.update(verdict='pending', quantum_chemical_optimum=None,
                   reason='no CA2_small_opt.fchk and no opt_final.npy at the paths given: ' + ', '.join(
                       list(args.fchk or []) + list(args.opt_geometry or [])))
        code = 2
    else:
        q = site_geometry(xyz, names)
        res['quantum_chemical_optimum'] = dict(source=src, **q)
        res['comparison'] = {
            k: dict(force_field=mmv[k], quantum_chemistry=q[k], difference=round(mmv[k] - q[k], 4), tolerance=TOL[k],
                    passed=bool(abs(mmv[k] - q[k]) <= TOL[k])) for k in TOL}
        res['quantum_chemistry_against_crystal'] = {k: round(q[k] - CRYSTAL[k][0], 4) for k in TOL}
        res['optimum_against_crystal'] = {
            k: dict(quantum_chemistry=q[k], crystal=CRYSTAL[k][0], difference=round(q[k] - CRYSTAL[k][0], 4), tolerance=TOL[k],
                    within=bool(abs(q[k] - CRYSTAL[k][0]) <= TOL[k])) for k in TOL}
        res['level_of_the_optimum'] = ('B3LYP/6-31G*, vacuum, all atoms free, WITHOUT a dispersion term (the energies of gate G-R '
                                       'carry D3(BJ); this geometry does not)')
        ok = all(c['passed'] for c in res['comparison'].values())
        res.update(verdict=('passed' if ok else 'failed') + (' (provisional)' if src['provisional'] else ''),
                   passed=bool(ok), provisional=src['provisional'])
        code = 0 if ok else 1
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix('.tmp')
    tmp.write_text(json.dumps(res, indent=1) + '\n')
    os.replace(tmp, out)
    line = f"binding check: {res['verdict']}"
    if src is not None:
        line += '; ' + ', '.join(f"{k} MM {c['force_field']} QM {c['quantum_chemistry']} ({c['difference']:+.3f}, tolerance {c['tolerance']})"
                                 for k, c in res['comparison'].items())
        line += '; the optimum against 6RL9: ' + ', '.join(
            f"{k} {c['difference']:+.3f} (tolerance {c['tolerance']}, within {c['within']})" for k, c in res['optimum_against_crystal'].items())
        line += '; by the plan (1:09 PM) the ligand terms are not re-tuned to this optimum'
    line += '; force field against 6RL9 (does not decide): ' + ', '.join(
        f"{k} {c['force_field']} / {c['crystal']}" for k, c in res['against_crystal_not_the_verdict'].items())
    print(line)
    return code


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    x = sub.add_parser('relax')
    x.add_argument('--mm', required=True)
    x.add_argument('--out', required=True)
    x.add_argument('--geom-out', required=True)
    x.add_argument('--qm-geom', nargs='*')
    x.add_argument('--reuse', action='store_true', help='keep the minimisations of OUT.json when it was made with the same model')
    b = sub.add_parser('binding')
    b.add_argument('--mm', required=True)
    b.add_argument('--out', required=True)
    b.add_argument('--fchk', nargs='*')
    b.add_argument('--opt-geometry', nargs='*')
    r = sub.add_parser('run')
    r.add_argument('--params', required=True)
    r.add_argument('--work', required=True)
    r.add_argument('--out', required=True)
    r.add_argument('--label', required=True)
    r.add_argument('--extra', nargs='*')
    r.add_argument('--n2h', choices=('structure', 'asis'), default='structure')
    r.add_argument('--points', default=str(HERE / 'gr' / 'points.json'))
    r.add_argument('--no-min', action='store_true')
    args = ap.parse_args()
    if args.cmd == 'relax':
        relax(args)
    elif args.cmd == 'binding':
        sys.exit(binding(args))
    else:
        run(args)


if __name__ == '__main__':
    main()
