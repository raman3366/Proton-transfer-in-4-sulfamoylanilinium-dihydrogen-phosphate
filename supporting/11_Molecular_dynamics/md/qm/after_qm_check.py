"""Checks after MCPB.py steps 2-4 and the tleap run of md/qm/after_qm.sh (called by that script).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/after_qm_check.py WORK REPORT STEP2

Reads in WORK: CA2_mcpbpy.frcmod, the five mol2 files of step 3, CA2_tleap_dry.log, CA2_dry.prmtop/.inpcrd,
CA2_standard.pdb, SAN.mol2.  Writes REPORT/after_qm_check_<step2>.json.  Exit code 1 if a hard check fails.

Hard checks (the script stops the pipeline when one fails)
  H1  every ligand term that carries the ligating nitrogen has, under the new type, the value of
      md/LIGAND.md section 6 (file: the table TERMS below, taken from that section)
  H2  tleap: no fatal error, no missing parameter; 4,094 atoms; net charge -1
  H3  the topology has the four Zn-N bonds, the 14 angles with zinc and only zero torsions through zinc
  H4  charges: the five residues of the site sum to +1; no charge was left at the value of the input
Reported, not hard (they are properties of the model, to be judged, and the gate G-Z comes later)
  the zinc terms (force constants, equilibrium values) beside the 1.0 A structure;
  MCPB.py's five rules of thumb (tutorial, section 5);
  the RESP charges of the ligand beside the AM1-BCC charges of stage S1;
  a minimisation in vacuum of everything within 12 A of the zinc (the rest fixed): zinc geometry before
  and after.  This is a first look, NOT gate G-Z (which needs the solvated system and 100 ps).
"""
import json
import pathlib
import re
import sys

import numpy as np
import parmed

WORK, REPORT, STEP2 = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), sys.argv[3]
GROUP = 'CA2'
SITE = ('HD1', 'HD2', 'HE1', 'ZN1', 'SN1')
# md/LIGAND.md, section 6: kind, types with n2 -> Y4, values that must arrive
TERMS = [('bond', ('Y4', 'sy'), (461.21, 1.5906)),
         ('angle', ('ca', 'sy', 'Y4'), (101.19, 117.133)),
         ('angle', ('hn', 'Y4', 'sy'), (61.55, 118.059)),
         ('angle', ('Y4', 'sy', 'o'), (132.73, 122.18)),
         ('bond', ('Y4', 'hn'), (454.59, 1.022)),
         ('dihedral', ('hn', 'Y4', 'sy', 'ca'), [(0.2624, 180.0, 1), (2.7030, 0.0, 2)]),
         ('dihedral', ('hn', 'Y4', 'sy', 'o'), [(0.0, None, None)])]
CRYSTAL = {'bonds': {'HD1@NE2': 2.010, 'HD2@NE2': 2.021, 'HE1@ND1': 2.018, 'SN1@N2': 1.931}}
out = dict(step2=STEP2, hard=[], reported={})
fail = False


def hard(name, ok, **kw):
    global fail
    out['hard'].append(dict(check=name, passed=bool(ok), **kw))
    print(('PASS ' if ok else 'FAIL ') + name)
    fail = fail or not ok


def read_frcmod(path):
    sec, d = None, dict(MASS=[], BOND=[], ANGL=[], DIHE=[], IMPR=[], NONB=[])
    for l in pathlib.Path(path).read_text().splitlines()[1:]:
        key = l.strip()[:4]
        if key in d and len(l.strip()) <= 8:
            sec = key
            continue
        if sec and l.strip():
            d[sec].append(l)
    return d


def types_of(line, n):
    head = line[:3 * n - 1]
    return tuple(t.strip() for t in head.split('-')), line[3 * n - 1:].split()


def same(a, b):
    return tuple(a) == tuple(b) or tuple(a) == tuple(b)[::-1]


# ----------------------------------------------------------------------------------------------------------- H1
frc = read_frcmod(WORK / f'{GROUP}_mcpbpy.frcmod')
rows = []
for kind, types, want in TERMS:
    if kind == 'bond':
        hit = [types_of(l, 2)[1] for l in frc['BOND'] if same(types_of(l, 2)[0], types)]
        got = (float(hit[0][0]), float(hit[0][1])) if len(hit) == 1 else None
        ok = got is not None and abs(got[0] - want[0]) < 5e-3 and abs(got[1] - want[1]) < 5e-4
    elif kind == 'angle':
        hit = [types_of(l, 3)[1] for l in frc['ANGL'] if same(types_of(l, 3)[0], types)]
        got = (float(hit[0][0]), float(hit[0][1])) if len(hit) == 1 else None
        ok = got is not None and abs(got[0] - want[0]) < 5e-3 and abs(got[1] - want[1]) < 5e-4
    else:
        hit = [types_of(l, 4)[1] for l in frc['DIHE'] if same(types_of(l, 4)[0], types)]
        got = [(float(h[1]) / float(h[0]), float(h[2]), abs(int(float(h[3])))) for h in hit]
        if want[0][1] is None:
            ok = len(got) >= 1 and all(abs(g[0]) < 1e-9 for g in got)
        else:
            ok = len(got) == len(want) and all(
                any(abs(g[0] - w[0]) < 5e-5 and abs(g[1] - w[1]) < 1e-3 and g[2] == w[2] for g in got) for w in want)
    rows.append(dict(kind=kind, types='-'.join(types), wanted=want, arrived=got, passed=bool(ok)))
hard('H1 ligand terms with the ligating nitrogen arrive under the new type', all(r['passed'] for r in rows), terms=rows)

# ----------------------------------------------------------------------------------------------------------- H2
log = (WORK / f'{GROUP}_tleap_dry.log').read_text()
missing = [l.strip() for l in log.splitlines() if re.search(r'Could not find|could not find|No torsion terms|'
                                                            r'FATAL|Fatal|\*\* No', l)]
top = parmed.load_file(str(WORK / f'{GROUP}_dry.prmtop'), str(WORK / f'{GROUP}_dry.inpcrd'))
qtot = float(sum(a.charge for a in top.atoms))
hard('H2 tleap: no error, 4,094 atoms, net charge -1', not missing and len(top.atoms) == 4094 and abs(qtot + 1.0) < 2e-3,
     messages=missing[:10], n_atoms=len(top.atoms), net_charge=round(qtot, 4),
     warnings=[l.strip() for l in log.splitlines() if 'Warning' in l or 'WARNING' in l][:6])

# ----------------------------------------------------------------------------------------------------------- H3
zn = [a for a in top.atoms if a.residue.name == 'ZN1'][0]
xyz = top.coordinates
zb = [b for b in top.bonds if zn in (b.atom1, b.atom2)]
za = [a for a in top.angles if zn in (a.atom1, a.atom2, a.atom3)]
zd = [d for d in top.dihedrals if zn in (d.atom1, d.atom2, d.atom3, d.atom4) and not d.improper]
bonds = []
for b in zb:
    o = b.atom2 if b.atom1 is zn else b.atom1
    key = f'{o.residue.name}@{o.name}'
    bonds.append(dict(atom=key, residue_index=o.residue.idx + 1, k=round(b.type.k, 2), r0=round(b.type.req, 4),
                      r_structure=round(float(np.linalg.norm(xyz[o.idx] - xyz[zn.idx])), 4),
                      r_6RL9=CRYSTAL['bonds'].get(key)))
angles = []
for a in za:
    at = (a.atom1, a.atom2, a.atom3)
    u, v = xyz[at[0].idx] - xyz[at[1].idx], xyz[at[2].idx] - xyz[at[1].idx]
    th = float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))
    angles.append(dict(atoms='-'.join(f'{x.residue.name}@{x.name}' for x in at), k=round(a.type.k, 2),
                       theta0=round(a.type.theteq, 2), theta_structure=round(th, 2)))
nonzero = [d for d in zd if abs(d.type.phi_k) > 1e-9]
hard('H3 four Zn-N bonds, 14 angles with zinc, torsions through zinc all zero',
     len(zb) == 4 and len(za) == 14 and not nonzero and all(b['k'] > 0 for b in bonds) and all(a['k'] > 0 for a in angles),
     n_bonds=len(zb), n_angles=len(za), n_torsions_with_zinc=len(zd), n_nonzero_torsions=len(nonzero))
out['reported']['zinc_bonds'] = bonds
out['reported']['zinc_angles'] = sorted(angles, key=lambda r: r['atoms'])

# ----------------------------------------------------------------------------------------------------------- H4
site = {}
for r in top.residues:
    if r.name in SITE:
        site[r.name] = round(float(sum(a.charge for a in r.atoms)), 4)
lig_in = parmed.load_file(str(WORK / 'SAN.mol2'))
lig = [r for r in top.residues if r.name == 'SN1'][0]
qin = {a.name: a.charge for a in lig_in.atoms}
ligq = [dict(atom=a.name, type=a.type, resp=round(a.charge, 4), am1bcc=round(qin[a.name], 4)) for a in lig.atoms]
unchanged = sum(abs(r['resp'] - r['am1bcc']) < 1e-4 for r in ligq)
hard('H4 charges of the site: five residues sum to +1, ligand charges refitted',
     abs(sum(site.values()) - 1.0) < 2e-3 and unchanged == 0, residue_charges=site,
     sum=round(sum(site.values()), 4), ligand_charges_left_at_input_value=unchanged)
out['reported']['ligand_charges'] = ligq
out['reported']['zinc_charge'] = round(zn.charge, 4)

# rules of thumb of the MCPB.py tutorial
out['reported']['mcpb_rules_of_thumb'] = dict(
    A_bond_k_below_200_and_r_below_2p8=all(b['k'] < 200 and b['r0'] < 2.8 for b in bonds),
    B_angle_k_below_100_and_theta_above_100=all(a['k'] < 100 and a['theta0'] > 100 for a in angles),
    B_exceptions=[a for a in angles if not (a['k'] < 100 and a['theta0'] > 100)],
    C_torsions_through_zinc_zero=not nonzero,
    D_zinc_charge_below_plus_one=bool(zn.charge < 1.0), zinc_charge=round(zn.charge, 4),
    E_zinc_lj_rmin_half_above_1=bool(zn.rmin > 1.0), zinc_rmin_half=round(zn.rmin, 4))

# ----------------------------------------------------------------- a first look: minimisation near the zinc
try:
    import openmm
    from openmm import app, unit
    system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    near = np.linalg.norm(xyz - xyz[zn.idx], axis=1) < 12.0
    for i in range(system.getNumParticles()):
        if not near[i]:
            system.setParticleMass(i, 0.0)
    integ = openmm.VerletIntegrator(0.001)
    ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName('Reference') if False else
                         openmm.Platform.getPlatformByName('CPU'), {'Threads': '1'})
    ctx.setPositions(top.positions)
    e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 5000)
    st = ctx.getState(getEnergy=True, getPositions=True)
    x1 = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    e1 = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    ligating = [b.atom2 if b.atom1 is zn else b.atom1 for b in zb]
    look = dict(note='vacuum, no cutoff, atoms beyond 12 A of the zinc fixed; a first look, NOT gate G-Z',
                n_free_atoms=int(near.sum()), energy_before=round(e0, 1), energy_after=round(e1, 1), bonds=[], angles=[])
    for o in ligating:
        look['bonds'].append(dict(atom=f'{o.residue.name}@{o.name}',
                                  before=round(float(np.linalg.norm(xyz[o.idx] - xyz[zn.idx])), 3),
                                  after=round(float(np.linalg.norm(x1[o.idx] - x1[zn.idx])), 3)))
    for i in range(4):
        for k in range(i + 1, 4):
            def ang(x):
                u, v = x[ligating[i].idx] - x[zn.idx], x[ligating[k].idx] - x[zn.idx]
                return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))
            look['angles'].append(dict(atoms=f'{ligating[i].residue.name}@{ligating[i].name}-Zn-'
                                             f'{ligating[k].residue.name}@{ligating[k].name}',
                                       before=round(ang(xyz), 1), after=round(ang(x1), 1)))
    lg = {a.name: a.idx for a in lig.atoms}

    def a3(x, i, j, k):
        u, v = x[i] - x[j], x[k] - x[j]
        return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))

    def dih(x, i, j, k, l):
        b0, b1, b2 = x[i] - x[j], x[k] - x[j], x[l] - x[k]
        b1 = b1 / np.linalg.norm(b1)
        v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
        return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

    look['ligand'] = dict(
        C4_S_N2=dict(before=round(a3(xyz, lg['C4'], lg['S'], lg['N2']), 2), after=round(a3(x1, lg['C4'], lg['S'], lg['N2']), 2)),
        S_N2=dict(before=round(float(np.linalg.norm(xyz[lg['S']] - xyz[lg['N2']])), 4),
                  after=round(float(np.linalg.norm(x1[lg['S']] - x1[lg['N2']])), 4)),
        Zn_N2_S=dict(before=round(a3(xyz, zn.idx, lg['N2'], lg['S']), 2), after=round(a3(x1, zn.idx, lg['N2'], lg['S']), 2)),
        Zn_N2_H=dict(before=round(a3(xyz, zn.idx, lg['N2'], lg['HN21']), 2), after=round(a3(x1, zn.idx, lg['N2'], lg['HN21']), 2)),
        Zn_N2_S_C4=dict(before=round(dih(xyz, zn.idx, lg['N2'], lg['S'], lg['C4']), 1),
                        after=round(dih(x1, zn.idx, lg['N2'], lg['S'], lg['C4']), 1)),
        C3_C4_S_N2=dict(before=round(dih(xyz, lg['C3'], lg['C4'], lg['S'], lg['N2']), 1),
                        after=round(dih(x1, lg['C3'], lg['C4'], lg['S'], lg['N2']), 1)),
        heavy_atom_rmsd_no_fit=round(float(np.sqrt(((x1 - xyz)[[a.idx for a in lig.atoms if a.atomic_number > 1]] ** 2)
                                                   .sum(axis=1).mean())), 3))
    out['reported']['first_look_minimisation'] = look
except Exception as e:                                                     # noqa: BLE001
    out['reported']['first_look_minimisation'] = dict(error=f'{type(e).__name__}: {e}')

out['all_hard_checks_passed'] = not fail
REPORT.mkdir(parents=True, exist_ok=True)
(REPORT / f'after_qm_check_{STEP2}.json').write_text(json.dumps(out, indent=1) + '\n')
sys.exit(1 if fail else 0)
