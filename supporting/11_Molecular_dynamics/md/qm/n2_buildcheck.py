"""Gate G-R2: a THROWAWAY build of start X with a file of terms at N2, and the comparison of every parameter with
the system as built (md/N2.md, section 6.6).  The built system is a test object and must never be simulated.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_buildcheck.py
        --frcmod md/qm/n2/n2_terms_new.frcmod --work SCRATCH [--reference /root/cqc-sfa-md/systems/bonded_2e/X]
        [--mcpb-dir DIR] [--out md/qm/n2/buildcheck.json] [--frames 3] [--charges-changed]
    --charges-changed (the declared second step: --mcpb-dir is then the folder of md/qm/n2_charges.py mol2): the
    charges of the five residues of the zinc site may differ from the reference build, and must be those of the mol2
    files of that folder; nothing else of an atom may differ, and no other atom at all.

 1  md/build_topology.py is IMPORTED (not edited) and its function build() is called for the bonded model, set 2e,
    start X, with --n2-frcmod = the file, into SCRATCH/systems and SCRATCH/records (so that neither
    /root/cqc-sfa-md/systems nor md/topology is touched).  The builder's own checks are reported as they come.
 2  The new system.prmtop against the one of the reference build, section by section of the file, and term by term
    (ParmEd): which sections are identical, which terms differ, and whether every difference is one the plan
    allows (the three angles at N2; torsions about N2-S with Zn or H at one end and C4, O1 or O2 at the other;
    torsions about Zn-N2 with a histidine nitrogen at one end and S or H at the other; one improper at N2).
 3  Every changed term in the topology against the number in the file.
 4  The energy: E(new topology) - E(reference topology), bonded terms only (OpenMM), at the coordinates as built
    and at frames of the rehearsal run TX0, against the same difference computed here from the numbers of the two
    frcmod files and the geometry.
"""
import argparse
import collections
import json
import pathlib
import sys

import numpy as np

import n2_common as nc

ALLOWED_SECTIONS = {'TITLE', 'POINTERS', 'ANGLE_FORCE_CONSTANT', 'ANGLE_EQUIL_VALUE', 'ANGLES_INC_HYDROGEN',
                    'ANGLES_WITHOUT_HYDROGEN', 'DIHEDRAL_FORCE_CONSTANT', 'DIHEDRAL_PERIODICITY', 'DIHEDRAL_PHASE',
                    'SCEE_SCALE_FACTOR', 'SCNB_SCALE_FACTOR', 'DIHEDRALS_INC_HYDROGEN', 'DIHEDRALS_WITHOUT_HYDROGEN'}


def sections(path):
    out, name = {}, None
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith('%FLAG'):
            name = line.split()[1]
            out[name] = []
        elif line.startswith('%VERSION') or line.startswith('%FORMAT') or line.startswith('%COMMENT'):
            continue
        elif name:
            out[name].append(line)
    return {k: '\n'.join(v) for k, v in out.items()}


def terms(parm):
    def canon(*i):
        return i if i <= i[::-1] else i[::-1]
    t = dict(atoms=[(a.name, a.type, round(a.charge, 6), round(a.mass, 4), round(a.rmin, 5), round(a.epsilon, 6),
                     a.residue.name) for a in parm.atoms],
             bonds=collections.Counter((canon(b.atom1.idx, b.atom2.idx), round(b.type.k, 4), round(b.type.req, 5)) for b in parm.bonds),
             angles={}, torsions=collections.defaultdict(list), impropers=collections.defaultdict(list), pairs14=set())
    for a in parm.angles:
        t['angles'][canon(a.atom1.idx, a.atom2.idx, a.atom3.idx)] = (round(a.type.k, 4), round(a.type.theteq, 4))
    for d in parm.dihedrals:
        idx = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
        term = (int(d.type.per), round(d.type.phi_k, 5), round(d.type.phase % 360.0, 3))
        if d.improper:
            t['impropers'][tuple(sorted(idx)) + (d.atom3.idx,)].append(term)
        else:
            t['torsions'][canon(*idx)].append(term)
            if not d.ignore_end:
                t['pairs14'].add((tuple(sorted((idx[0], idx[3]))), round(d.type.scee, 4), round(d.type.scnb, 4)))
    return t


def read_n2_file(path):
    """Angles, torsion series and the improper of a frcmod of terms at N2, by atom types."""
    sec, out = None, dict(ANGLE={}, DIHE=collections.defaultdict(list), IMPROPER={})
    for line in pathlib.Path(path).read_text().splitlines()[1:]:
        w = line.strip()
        if w[:4] in ('MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB'):
            sec = {'ANGL': 'ANGLE', 'IMPR': 'IMPROPER'}.get(w[:4], w[:4])
            continue
        if not w or sec not in out:
            continue
        n = dict(ANGLE=3, DIHE=4, IMPROPER=4)[sec]
        types = tuple(x.strip() for x in line[:3 * n - 1].split('-'))
        num = line[3 * n - 1:].split()
        if sec == 'ANGLE':
            out[sec][types] = (float(num[0]), float(num[1]))
        elif sec == 'DIHE':
            out[sec][types].append((abs(int(float(num[3]))), float(num[1]) / float(num[0]), float(num[2]) % 360.0))
        else:
            out[sec][types] = (int(float(num[2])), float(num[0]), float(num[1]) % 360.0)
    return out


def site(parm):
    ix = {(a.residue.name, a.name): a.idx for a in parm.atoms if a.residue.name in ('ZN1', 'SN1', 'HD1', 'HD2', 'HE1')}
    return dict(zn=ix[('ZN1', 'ZN')], n2=ix[('SN1', 'N2')], s=ix[('SN1', 'S')], h=ix[('SN1', 'HN21')], c4=ix[('SN1', 'C4')],
                o1=ix[('SN1', 'O1')], o2=ix[('SN1', 'O2')], n94=ix[('HD1', 'NE2')], n96=ix[('HD2', 'NE2')], n119=ix[('HE1', 'ND1')])


def n2_energy(frc, types, st, x):
    """Energy of the terms of a terms-at-N2 file at the geometry x (kcal/mol), in the AMBER forms."""
    def look(table, tt):
        return table.get(tt, table.get(tt[::-1]))
    e = 0.0
    for trio in (('zn', 'n2', 's'), ('zn', 'n2', 'h'), ('s', 'n2', 'h')):
        par = look(frc['ANGLE'], tuple(types[k] for k in trio))
        th = np.radians(nc.angle(*(x[st[k]] for k in trio)))
        e += par[0] * (th - np.radians(par[1])) ** 2
    quads = [('c4', 's', 'n2', 'h'), ('zn', 'n2', 's', 'c4'), ('zn', 'n2', 's', 'o1'), ('zn', 'n2', 's', 'o2'),
             ('h', 'n2', 's', 'o1'), ('h', 'n2', 's', 'o2')]
    quads += [(hn, 'zn', 'n2', end) for hn in ('n94', 'n96', 'n119') for end in ('s', 'h')]
    for q in quads:
        ser = look(frc['DIHE'], tuple(types[k] for k in q))
        if not ser:
            continue
        phi = np.radians(nc.dihedral(*(x[st[k]] for k in q)))
        e += sum(pk * (1.0 + np.cos(n * phi - np.radians(ph))) for n, pk, ph in ser)
    for tt, (n, pk, ph) in frc['IMPROPER'].items():
        om = np.radians(nc.dihedral(x[st['zn']], x[st['h']], x[st['n2']], x[st['s']]))
        e += pk * (1.0 + np.cos(n * om - np.radians(ph)))
    return float(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--frcmod', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--reference', default='/root/cqc-sfa-md/systems/bonded_2e/X')
    ap.add_argument('--mcpb-dir', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--frames', type=int, default=3)
    ap.add_argument('--reuse', action='store_true', help='compare a throwaway build that is there already')
    ap.add_argument('--charges-changed', action='store_true')
    args = ap.parse_args()
    import parmed
    work = pathlib.Path(args.work).resolve()
    if args.mcpb_dir:
        # the builder writes --mcpb-dir into the tleap input as given, and tleap runs in the build folder: a relative
        # path is not found there (seen 6 Oct 2026); it must be absolute
        args.mcpb_dir = str(pathlib.Path(args.mcpb_dir).resolve())
    assert nc.REPO not in work.parents and '/root/cqc-sfa-md/systems' not in str(work), 'the work folder must be scratch'
    ref = pathlib.Path(args.reference)
    sys.path.insert(0, str(nc.REPO / 'md'))
    import build_topology as bt
    new = work / 'systems' / 'bonded_2e' / 'X'
    builder = dict(script='md/build_topology.py', sha256=nc.sha(nc.REPO / 'md' / 'build_topology.py'))
    if not (args.reuse and (new / 'system.prmtop').exists()):
        try:
            bj = bt.build('bonded', 'X', str(work / 'systems'), mcpb_dir=args.mcpb_dir, records_root=str(work / 'records'),
                          pset='2e', n2_frcmod=str(pathlib.Path(args.frcmod).resolve()),
                          note='THROWAWAY build for the check of a file of terms at N2 (md/qm/n2_buildcheck.py); never to be simulated')
            builder.update(checks_failed=bj['checks_failed'], n_checks=len(bj['checks']),
                           checks_passed=sum(v is True for v in bj['checks'].values()),
                           terms_at_N2_record=bj.get('terms_at_N2', {}).get('in_the_topology'),
                           torsions_through_zinc=bj.get('zinc_site_bonded', {}).get('torsions_through_zinc'))
        except bt.BuildError as e:
            builder.update(build_error=str(e))
            print('the builder stopped:', e)
    else:
        bj = json.loads((new / 'build.json').read_text())
        builder.update(checks_failed=bj['checks_failed'], n_checks=len(bj['checks']), reused=True,
                       terms_at_N2_record=bj.get('terms_at_N2', {}).get('in_the_topology'),
                       torsions_through_zinc=bj.get('zinc_site_bonded', {}).get('torsions_through_zinc'))
    rec = dict(what='throwaway build of start X with a file of terms at N2, against the system as built; written by '
                    'md/qm/n2_buildcheck.py; the built system is a test object', written=nc.now(),
               frcmod=dict(file=str(args.frcmod), sha256=nc.sha(args.frcmod)), reference=str(ref),
               reference_n2_terms_sha256=nc.sha(ref / 'build' / 'n2_terms.frcmod'), work=str(work), builder=builder)
    if not (new / 'system.prmtop').exists():
        rec['passed'] = False
        nc.write_json(args.out or work / 'buildcheck.json', rec)
        return 1
    # ---- 2a: sections of the prmtop file
    s_old, s_new = sections(ref / 'system.prmtop'), sections(new / 'system.prmtop')
    same = sorted(k for k in s_old if s_new.get(k) == s_old[k])
    diff = sorted(set(s_old) ^ set(s_new)) + sorted(k for k in s_old if k in s_new and s_new[k] != s_old[k])
    allowed_sec = ALLOWED_SECTIONS | ({'CHARGE'} if args.charges_changed else set())
    rec['prmtop_sections'] = dict(n=len(s_old), identical=len(same), differ=diff,
                                  differ_outside_angles_and_torsions=[k for k in diff if k not in allowed_sec],
                                  identical_sections=same)
    rec['coordinates_identical'] = bool(nc.sha(ref / 'system.inpcrd') == nc.sha(new / 'system.inpcrd'))
    # ---- 2b: term by term
    p_old = parmed.load_file(str(ref / 'system.prmtop'), xyz=str(ref / 'system.inpcrd'))
    p_new = parmed.load_file(str(new / 'system.prmtop'), xyz=str(new / 'system.inpcrd'))
    t_old, t_new = terms(p_old), terms(p_new)
    st = site(p_new)
    assert st == site(p_old)
    name = {v: k for k, v in st.items()}
    lab = lambda idx: '-'.join(f'{p_new.atoms[i].residue.name}{p_new.atoms[i].residue.number}:{p_new.atoms[i].name}' for i in idx)  # noqa: E731
    other, allowed = [], dict(angles={}, torsions={}, improper={})
    if t_old['atoms'] != t_new['atoms']:
        if not args.charges_changed:
            other.append('atoms (names, types, charges, masses, Lennard-Jones)')
        else:
            # the declared second step: charges of the five site residues only, and as the mol2 files have them
            site_res = {'HD1', 'HD2', 'HE1', 'ZN1', 'SN1'}
            want = {}
            for rn in site_res:
                sec = None
                for l in (pathlib.Path(args.mcpb_dir) / f'{rn}.mol2').read_text().splitlines():
                    if l.startswith('@<TRIPOS>'):
                        sec = l.strip()[9:]
                        continue
                    f = l.split()
                    if sec == 'ATOM' and len(f) >= 9:
                        want[(rn, f[1])] = float(f[8])
            n_changed, worst_other, worst_file, bad_other = 0, 0.0, 0.0, 0
            for a, b in zip(t_old['atoms'], t_new['atoms']):
                if a[:2] + a[3:] != b[:2] + b[3:]:
                    bad_other += 1
                    continue
                if b[6] in site_res:
                    worst_file = max(worst_file, abs(b[2] - want[(b[6], b[0])]))
                    n_changed += abs(a[2] - b[2]) > 1e-7
                elif abs(a[2] - b[2]) > 1e-7:
                    worst_other = max(worst_other, abs(a[2] - b[2]))
                    bad_other += 1
            rec['charges'] = dict(site_atoms_with_another_charge=int(n_changed), atoms_outside_the_site_that_differ=int(bad_other),
                                  largest_difference_from_the_mol2_files=round(worst_file, 7),
                                  site_charge_new=round(sum(b[2] for b in t_new['atoms'] if b[6] in site_res), 5),
                                  site_charge_reference=round(sum(a[2] for a in t_old['atoms'] if a[6] in site_res), 5),
                                  total_charge_new=round(sum(b[2] for b in t_new['atoms']), 5))
            if bad_other or worst_file > 2e-6:
                other.append('atoms: a difference other than the charges of the five site residues, or charges not those of the mol2 files')
    if t_old['bonds'] != t_new['bonds']:
        other.append('bonds')
    if t_old['pairs14'] != t_new['pairs14']:
        other.append('the 1-4 pairs or their scaling')
    for k in set(t_old['angles']) | set(t_new['angles']):
        if t_old['angles'].get(k) != t_new['angles'].get(k):
            ok = k[1] == st['n2'] and {k[0], k[2]} in ({st['zn'], st['s']}, {st['zn'], st['h']}, {st['s'], st['h']}) \
                and k in t_old['angles'] and k in t_new['angles']
            if ok:
                allowed['angles'][lab(k)] = dict(as_built=t_old['angles'].get(k), new=t_new['angles'].get(k))
            else:
                other.append(f'angle {lab(k)}')
    his = {st['n94'], st['n96'], st['n119']}
    for k in set(t_old['torsions']) | set(t_new['torsions']):
        a, b = sorted(t_old['torsions'].get(k, [])), sorted(t_new['torsions'].get(k, []))
        if a == b:
            continue
        mid, ends = {k[1], k[2]}, {k[0], k[3]}
        ok = (mid == {st['n2'], st['s']} and len(ends & {st['h'], st['zn']}) == 1 and len(ends & {st['c4'], st['o1'], st['o2']}) == 1) or \
             (mid == {st['zn'], st['n2']} and len(ends & his) == 1 and len(ends & {st['s'], st['h']}) == 1)
        if ok and k in t_old['torsions'] and k in t_new['torsions']:
            allowed['torsions'][lab(k)] = dict(as_built=a, new=b, atoms=[name.get(i, i) for i in k])
        else:
            other.append(f'torsion {lab(k)}')
    for k in set(t_old['impropers']) | set(t_new['impropers']):
        a, b = sorted(t_old['impropers'].get(k, [])), sorted(t_new['impropers'].get(k, []))
        if a == b:
            continue
        if k[-1] == st['n2'] and set(k[:4]) == {st['zn'], st['n2'], st['s'], st['h']} and not a and len(b) == 1:
            allowed['improper'][lab(k[:4])] = b
        else:
            other.append(f'improper {lab(k[:4])}')
    rec['term_by_term'] = dict(atoms=len(t_new['atoms']), bonds=sum(t_new['bonds'].values()), angles=len(t_new['angles']),
                               torsion_quadruples=len(t_new['torsions']), improper_quadruples=len(t_new['impropers']),
                               pairs_1_4=len(t_new['pairs14']), pairs_1_4_identical=bool(t_old['pairs14'] == t_new['pairs14']),
                               differences_the_plan_allows=allowed, differences_the_plan_does_not_allow=other)
    # ---- 3: the changed terms against the file
    frc_new, frc_old = read_n2_file(args.frcmod), read_n2_file(ref / 'build' / 'n2_terms.frcmod')
    types = {k: p_new.atoms[i].type for k, i in st.items()}
    arrive = []

    def look(table, tt):
        return table.get(tt, table.get(tt[::-1]))
    for trio in (('zn', 'n2', 's'), ('zn', 'n2', 'h'), ('s', 'n2', 'h')):
        k = tuple(st[x] for x in trio)
        k = k if k <= k[::-1] else k[::-1]
        want, got = look(frc_new['ANGLE'], tuple(types[x] for x in trio)), t_new['angles'][k]
        arrive.append(dict(term='angle ' + '-'.join(trio), file=want, topology=got,
                           same=bool(want is not None and abs(want[0] - got[0]) < 5e-4 and abs(want[1] - got[1]) < 5e-4)))
    quads = [('c4', 's', 'n2', 'h'), ('zn', 'n2', 's', 'c4'), ('zn', 'n2', 's', 'o1'), ('zn', 'n2', 's', 'o2'),
             ('h', 'n2', 's', 'o1'), ('h', 'n2', 's', 'o2')] + [(hn, 'zn', 'n2', e) for hn in ('n94', 'n96', 'n119') for e in ('s', 'h')]
    for q in quads:
        k = tuple(st[x] for x in q)
        k = k if k <= k[::-1] else k[::-1]
        want = look(frc_new['DIHE'], tuple(types[x] for x in q))
        got = sorted(t for t in t_new['torsions'].get(k, []))
        if want is None:
            same = sorted(t_old['torsions'].get(k, [])) == got
            arrive.append(dict(term='torsion ' + '-'.join(q), file='not in the file: as built', topology=got, same=bool(same)))
            continue
        w = sorted((n, round(pk, 5), round(ph, 3)) for n, pk, ph in want)
        nz = lambda L: sorted((n, pk, (ph if pk else 0.0)) for n, pk, ph in L if abs(pk) > 1e-9)   # noqa: E731
        arrive.append(dict(term='torsion ' + '-'.join(q), file=w, topology=got,
                           same=bool(len(got) == len(w) and all(a[0] == b[0] and abs(a[1] - b[1]) < 5e-5 and (abs(a[1]) < 1e-9 or abs(a[2] - b[2]) < 1e-2)
                                                                for a, b in zip(sorted(got), w))),
                           non_zero_terms_same=bool(nz(got) == nz(w))))
    for tt, (n, pk, ph) in frc_new['IMPROPER'].items():
        got = [v for k, v in t_new['impropers'].items() if k[-1] == st['n2'] and set(k[:4]) == {st['zn'], st['n2'], st['s'], st['h']}]
        arrive.append(dict(term='improper at N2 ' + '-'.join(tt), file=(n, pk, ph), topology=got,
                           same=bool(len(got) == 1 and len(got[0]) == 1 and got[0][0][0] == n and abs(got[0][0][1] - pk) < 5e-5
                                     and abs(got[0][0][2] - ph) < 1e-2)))
    if not frc_new['IMPROPER']:
        got = [v for k, v in t_new['impropers'].items() if k[-1] == st['n2']]
        arrive.append(dict(term='improper at N2', file='none in the file', topology=got, same=not got))
    rec['terms_in_the_topology_against_the_file'] = arrive
    # ---- 4: energies of the bonded terms
    import openmm
    from openmm import app, unit as u

    def bonded_context(prmtop):
        top = app.AmberPrmtopFile(str(prmtop))
        sysm = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False, removeCMMotion=False)
        for k in reversed(range(sysm.getNumForces())):
            if type(sysm.getForce(k)).__name__ not in ('HarmonicAngleForce', 'PeriodicTorsionForce'):
                sysm.removeForce(k)
        return openmm.Context(sysm, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('Reference'))
    c_old, c_new = bonded_context(ref / 'system.prmtop'), bonded_context(new / 'system.prmtop')
    geoms = [('as built', np.array(p_new.coordinates, float))]
    if args.frames:
        import MDAnalysis as mda
        run = pathlib.Path('/root/cqc-sfa-md/runs_test/bonded_2e/TX0')
        uni = mda.Universe(str(run / 'full.pdb'), str(run / 'full_seg0001.dcd'))
        assert uni.atoms.n_atoms == len(p_new.atoms)
        nfr = uni.trajectory.n_frames
        for k in np.linspace(0, nfr - 1, args.frames).astype(int):
            uni.trajectory[int(k)]
            geoms.append((f'TX0 production, full frame {int(k)} of {nfr}', uni.atoms.positions.astype(float).copy()))
    en = []
    for label, x in geoms:
        e = []
        for c in (c_old, c_new):
            c.setPositions(x * 0.1)
            e.append(c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole))
        model = n2_energy(frc_new, types, st, x) - n2_energy(frc_old, types, st, x)
        en.append(dict(geometry=label, topology_new_minus_as_built=round(e[1] - e[0], 5), from_the_two_files=round(model, 5),
                       difference=round(e[1] - e[0] - model, 6)))
    rec['energy_of_the_changed_terms'] = dict(points=en, largest_difference_kcal=max(abs(v['difference']) for v in en),
                                              note='angles and torsions only, all 40,902 atoms; the difference of the two '
                                                   'topologies must be the difference of the terms at N2 and nothing else')
    rec['passed'] = bool(not other and not rec['prmtop_sections']['differ_outside_angles_and_torsions']
                         and all(a['same'] for a in arrive) and rec['energy_of_the_changed_terms']['largest_difference_kcal'] < 2e-3
                         and rec['coordinates_identical'])
    rec['builder_own_checks_all_passed'] = bool(builder.get('checks_failed') == [])
    nc.write_json(args.out or work / 'buildcheck.json', rec)
    print(json.dumps(dict(passed=rec['passed'], builder_checks_failed=builder.get('checks_failed'),
                          sections_that_differ=diff, other_differences=other,
                          allowed={k: list(v) for k, v in allowed.items()},
                          terms_arrived=all(a['same'] for a in arrive), energy=rec['energy_of_the_changed_terms']['points']), indent=1))
    return 0 if rec['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
