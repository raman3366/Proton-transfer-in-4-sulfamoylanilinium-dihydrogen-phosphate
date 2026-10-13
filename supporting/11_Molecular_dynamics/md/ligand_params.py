"""Stage S1 of md/PLAN.md, first half: the reference data and GAFF2 as assigned.

  /root/cqc-sfa-tools/md/bin/python md/ligand_params.py            (from the repository root; about 15 minutes)

What it does, for the sulfonamidate anion (charge -1, one H on the sulfonamide N) and for neutral sulfanilamide:
 1. Reference data of the paper: the DFT minima and relaxed scans (molecules/sfa/results/<case>/), the ligand of
    the 1.0 A structure (SAN 304 of dock/pdb/6RL9.cif.gz)                    -> ligand/reference_geometry.csv,
                                                                               ligand/dft_scan_<species>.csv
 2. antechamber (GAFF2 atom types, AM1-BCC charges) on the DFT minimum        -> ligand/sfa_<species>.mol2
    parmchk2 against the gaff2.dat of this install ("-s 2")                   -> ligand/sfa_<species>_gaff2.frcmod
    tleap (leaprc.gaff2 + that frcmod)                                         -> ligand/sfa_<species>_gaff2.prmtop
    every bond, angle and torsion term at sulfur with its source               -> ligand/terms_gaff2_<species>.csv
    check: no term of the topology differs from a native gaff2.dat term        -> ligand/gaff2_as_assigned.json
 3. Minima from 145 starts and the relaxed scan of C3-C4-S-N2 (ligand_eval.py), in vacuum and in implicit water
                                                -> ligand/geometry_before.csv, ligand/scan_before_<species>.csv
 4. Control: the same with the frcmod that "parmchk2 -s gaff2" writes.  In this AmberTools 26 build that spelling
    reads gaff.dat (GAFF 1.81), although the help text lists it for GAFF2; loaded after leaprc.gaff2 its terms
    replace native GAFF2 terms.  This is the route of the scoping test of 5 Oct 2026.

Nothing is random except the 0.02 A jitter of the starts (seed in ligand_eval.SEED).  Bulky scratch (antechamber
and sqm files) goes to /root/cqc-sfa-md/ligand/.  A second run overwrites the same files with the same content.
"""
import json, pathlib, shutil, sys
import numpy as np
import ligand_common as lc
from ligand_common import OUT, BULK, RESULTS, ROOT, SPECIES, HEAVY, GEOM_ROWS
import ligand_eval as le


def reference_data():
    ref = {}
    for sp, d in SPECIES.items():
        el, x, title = lc.read_xyz(RESULTS / d['case'] / 'opt.xyz')
        names = lc.names_for(sp, el, x)
        ref[sp] = dict(names=names, x=x, title=title, file=str((RESULTS / d['case'] / 'opt.xyz').relative_to(ROOT)),
                       m=lc.add_means(lc.measures(x, names)), scan=lc.dft_scan(sp))
    el, x, title = lc.read_xyz(RESULTS / 'sulfanilamide_gas_d3' / 'opt.xyz')
    names = lc.names_for('neutral', el, x)
    ref['neutral_gas'] = dict(names=names, x=x, title=title, m=lc.add_means(lc.measures(x, names)),
                              file='molecules/sfa/results/sulfanilamide_gas_d3/opt.xyz')
    names, x, atoms = lc.crystal_ligand(304)
    ref['crystal'] = dict(names=names, x=x, m=lc.add_means(lc.measures(x, names)), atoms=atoms,
                          file='dock/pdb/6RL9.cif.gz (SAN 304, chain A, occupancy 1.00, no alternate location)')
    # stored second derivatives?
    hess = sorted(str(p.relative_to(ROOT)) for p in RESULTS.glob('*/freq.json')
                  if p.parent.name in ('sulfanilamide_gas_d3', 'sulfanilamide_pcm_d3', 'sulfonamidate_pcm_d3'))
    ref['hessian_note'] = dict(
        freq_files=hess,
        note='Only molecules/sfa/results/sulfanilamide_gas_d3/freq.json exists for these species: harmonic '
             'frequencies and IR intensities of neutral sulfanilamide in the gas phase (finite-difference Hessian '
             'of 27 Sep 2026). The Cartesian Hessian itself was not stored, and nothing exists for the anion or '
             'for the PCM minima. Angle force constants therefore cannot be taken from a stored Hessian; the '
             'GAFF2 force constants are kept.')
    rows = []
    for key, label, unit in GEOM_ROWS:
        rows.append(dict(quantity=label, unit=unit,
                         dft_anion_pcm=ref['anion']['m'].get(key), dft_neutral_pcm=ref['neutral']['m'].get(key),
                         dft_neutral_gas=ref['neutral_gas']['m'].get(key), crystal_6RL9_SAN304=ref['crystal']['m'].get(key)))
    lc.write_csv(OUT / 'reference_geometry.csv', rows)
    for sp in SPECIES:
        t = []
        for p in ref[sp]['scan']:
            m = p['m']
            t.append(dict(file=p['file'], kind=p['kind'], branch=p['branch'], tau_C3_C4_S_N2=m['tau_C3'],
                          tau_C5_C4_S_N2=m['tau_C5'], E_rel_ref_kcal=p['E_rel'], converged=p['converged'],
                          a_C4_S_N2=m['a_C4_S_N2'], d_S_N2=m['d_S_N2'], d_C4_S=m['d_C4_S'], a_O1_S_O2=m['a_O1_S_O2'],
                          tau_C4_S_N2_HN21=m.get('tau_C4_S_N2_HN21'), tau_C4_S_N2_HN22=m.get('tau_C4_S_N2_HN22'),
                          N2_angle_sum=m.get('N2_angle_sum')))
        lc.write_csv(OUT / f'dft_scan_{sp}.csv', t)
    return ref


def write_input_mol2(path, names, el, x):
    ix = {n: i for i, n in enumerate(names)}
    bonds = [(ix[a], ix[b], 'ar' if a[0] == 'C' and b[0] == 'C' else ('2' if b[0] == 'O' else '1'))
             for a, b in lc.BONDS_HEAVY]
    bonds += [(ix[lc.H_PARENT[h]], ix[h], '1') for h in names[11:]]
    with open(path, 'w') as f:
        f.write(f'@<TRIPOS>MOLECULE\nSAN\n{len(names)} {len(bonds)} 1 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n')
        for i, n in enumerate(names):
            f.write(f'{i + 1:4d} {n:<5s} {x[i][0]:10.5f} {x[i][1]:10.5f} {x[i][2]:10.5f} {el[i]:<3s} 1 SAN 0.0000\n')
        f.write('@<TRIPOS>BOND\n')
        for k, (i, j, t) in enumerate(bonds):
            f.write(f'{k + 1:4d} {i + 1:4d} {j + 1:4d} {t}\n')


def read_mol2(path):
    atoms, bonds, sec = [], [], None
    for l in open(path):
        if l.startswith('@<TRIPOS>'):
            sec = l.strip()[9:]; continue
        s = l.split()
        if sec == 'ATOM' and len(s) >= 9:
            atoms.append(dict(name=s[1], x=[float(v) for v in s[2:5]], type=s[5], charge=float(s[8])))
        elif sec == 'BOND' and len(s) >= 4:
            bonds.append((int(s[1]), int(s[2]), s[3]))
    return atoms, bonds


def assign(sp, ref):
    """antechamber, parmchk2 and tleap for one species."""
    d = SPECIES[sp]
    work = BULK / f'antechamber_{sp}'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    el, x, _ = lc.read_xyz(RESULTS / d['case'] / 'opt.xyz')
    names = ref[sp]['names']
    inp = OUT / f'sfa_{sp}_input.mol2'
    write_input_mol2(inp, names, el, x)
    shutil.copy(inp, work / 'in.mol2')
    lc.run(['antechamber', '-i', 'in.mol2', '-fi', 'mol2', '-o', f'sfa_{sp}.mol2', '-fo', 'mol2', '-at', 'gaff2',
            '-c', 'bcc', '-nc', d['charge'], '-m', 1, '-rn', 'SAN', '-pf', 'y', '-s', 2], work,
           log=work / 'antechamber.log')
    resid = fix_total_charge(work / f'sfa_{sp}.mol2', OUT / f'sfa_{sp}.mol2', d['charge'])
    shutil.copy(OUT / f'sfa_{sp}.mol2', work / f'sfa_{sp}.mol2')      # the topology is built from the delivered file
    shutil.copy(work / 'sqm.out', OUT / f'sfa_{sp}_sqm.out')
    # "-s 2" reads gaff2.dat; "-s gaff2" reads gaff.dat in this build (see the control)
    lc.run(['parmchk2', '-i', f'sfa_{sp}.mol2', '-f', 'mol2', '-o', f'sfa_{sp}_gaff2.frcmod', '-s', 2], work)
    lc.run(['parmchk2', '-i', f'sfa_{sp}.mol2', '-f', 'mol2', '-o', f'sfa_{sp}_gaff2_all.frcmod', '-s', 2,
            '-a', 'Y'], work)
    lc.run(['parmchk2', '-i', f'sfa_{sp}.mol2', '-f', 'mol2', '-o', f'control_s_gaff2.frcmod', '-s', 'gaff2'], work)
    # which parameter file each spelling opens (strace, when available)
    opened = {}
    if shutil.which('strace'):
        for flag in ('2', 'gaff2'):
            lc.subprocess.run(['strace', '-f', '-e', 'trace=openat', '-o', f'strace_{flag}.txt', 'parmchk2', '-i',
                               f'sfa_{sp}.mol2', '-f', 'mol2', '-o', 'strace_tmp.frcmod', '-s', flag],
                              cwd=str(work), env=lc.env(), stdout=lc.subprocess.DEVNULL, stderr=lc.subprocess.DEVNULL)
            txt = (work / f'strace_{flag}.txt').read_text().split('\n')
            opened[flag] = sorted({l.split('"')[1].split('/')[-1] for l in txt
                                   if 'leap/parm/' in l and 'ENOENT' not in l})
    for f in (f'sfa_{sp}_gaff2.frcmod', f'sfa_{sp}_gaff2_all.frcmod'):
        shutil.copy(work / f, OUT / f)
    build_topology(work, f'sfa_{sp}.mol2', f'sfa_{sp}_gaff2.frcmod', f'sfa_{sp}_gaff2')
    build_topology(work, f'sfa_{sp}.mol2', 'control_s_gaff2.frcmod', 'control_s_gaff2')
    for ext in ('prmtop', 'inpcrd'):
        shutil.copy(work / f'sfa_{sp}_gaff2.{ext}', OUT / f'sfa_{sp}_gaff2.{ext}')
    atoms, bonds = read_mol2(OUT / f'sfa_{sp}.mol2')
    assert [a['name'] for a in atoms] == names
    return dict(work=str(work), parmchk2_opened=opened, charge_rounding_residual_removed=resid,
                atoms=[dict(name=a['name'], type=a['type'], charge=a['charge']) for a in atoms],
                bonds=[(names[i - 1], names[j - 1], t) for i, j, t in bonds],
                total_charge=round(sum(a['charge'] for a in atoms), 6))


def fix_total_charge(src, dst, target):
    """antechamber reads the AM1 Mulliken charges from sqm.out, where they are printed to three decimals; their sum
    can miss the net charge by a few 0.001 e.  The residual is spread evenly over all atoms so that the delivered
    mol2 has exactly the net charge.  Returns the residual that was removed."""
    lines = pathlib.Path(src).read_text().split('\n')
    a0 = lines.index('@<TRIPOS>ATOM') + 1
    a1 = lines.index('@<TRIPOS>BOND')
    q = [float(l.split()[8]) for l in lines[a0:a1]]
    resid = sum(q) - target
    qn = [round(v - resid / len(q), 6) for v in q]
    qn[-1] = round(qn[-1] - (sum(qn) - target), 6)
    assert abs(sum(qn) - target) < 5e-7 and abs(resid) < 0.01, (sum(qn), resid)
    for k, l in enumerate(lines[a0:a1]):
        head = l[:l.rindex(l.split()[8])]
        lines[a0 + k] = f'{head}{qn[k]:9.6f}'
    pathlib.Path(dst).write_text('\n'.join(lines))
    return resid


def build_topology(work, mol2, frcmod, stem):
    (work / f'{stem}.leap.in').write_text(
        f'source leaprc.gaff2\nloadamberparams {frcmod}\nSAN = loadmol2 {mol2}\ncheck SAN\n'
        f'saveamberparm SAN {stem}.prmtop {stem}.inpcrd\nquit\n')
    out = lc.run(['tleap', '-f', f'{stem}.leap.in'], work, log=work / f'{stem}.leap.out')
    import re
    m = re.search(r'Errors = (\d+); Warnings = (\d+)', out)
    bad = [l for l in out.split('\n') if 'Could not find' in l or 'No torsion terms' in l or 'missing' in l.lower()
           or 'FATAL' in l]
    assert m and int(m.group(1)) == 0 and not bad and (work / f'{stem}.prmtop').exists(), (m, bad)
    # the only warnings allowed are those about the net charge of the unit (the anion)
    nwarn = out.count('Warning!')
    assert nwarn == out.count('is not zero'), out
    return out


def term_tables(prmtop, frcmod_path, gaff2, repair_path=None):
    """Every valence term of the topology with its source; flags any term that differs from a native gaff2.dat
    term of the same atom types.  Terms of a repair frcmod (ligand_fit.py) are labelled as such and not flagged."""
    import parmed
    parm = parmed.load_file(str(prmtop))
    frc = lc.read_frcmod(frcmod_path)
    rep = lc.read_frcmod(repair_path) if repair_path else None
    rows, overridden = [], []

    def src(kind, types, value):
        if rep is not None:
            rkey, rrec = lc.lookup(rep[kind], types)
            if rrec is not None:
                c = rrec[0]['comment'] if kind == 'dihedral' else rrec['comment']
                return 'repair frcmod: ' + c, None, True
        key, rec = lc.lookup(gaff2[kind], types)
        fkey, frec = lc.lookup(frc[kind], types)
        if rec is not None:                      # a native term of exactly these types exists
            native = [(r['pn'], round(r['k'], 5), r['phase']) for r in rec] if kind == 'dihedral' else \
                (round(rec['k'], 3), round(rec['eq'], 3))
            line = rec[0]['line'] if kind == 'dihedral' else rec['line']
            ok = (sorted(native) == sorted(value)) if kind == 'dihedral' else \
                (abs(native[0] - value[0]) < 0.02 and abs(native[1] - value[1]) < 0.006)
            if not ok:
                overridden.append(dict(kind=kind, types='-'.join(types), topology=value, native=native))
            return f'gaff2.dat line {line}, native' + ('' if ok else ' BUT OVERRIDDEN'), native, ok
        if frec is not None:
            c = frec[-1]['comment'] if kind == 'dihedral' else frec['comment']
            return 'parmchk2 by analogy: ' + c, None, True
        if kind == 'dihedral':
            key, rec = lc.lookup(gaff2[kind], types, wildcard=True)
            if rec is not None:
                return f'gaff2.dat line {rec[0]["line"]}, generic {"-".join(key)}', None, True
        return 'not found', None, False

    for b in parm.bonds:
        t = (b.atom1.type, b.atom2.type); v = (round(b.type.k, 3), round(b.type.req, 4))
        s, native, ok = src('bond', t, v)
        rows.append(dict(kind='bond', atoms=f'{b.atom1.name}-{b.atom2.name}', types='-'.join(t), k=b.type.k,
                         eq=b.type.req, n=None, phase=None, source=s,
                         at_sulfur='S' in (b.atom1.name, b.atom2.name)))
    for a in parm.angles:
        t = (a.atom1.type, a.atom2.type, a.atom3.type); v = (round(a.type.k, 3), round(a.type.theteq, 3))
        s, native, ok = src('angle', t, v)
        rows.append(dict(kind='angle', atoms=f'{a.atom1.name}-{a.atom2.name}-{a.atom3.name}', types='-'.join(t),
                         k=a.type.k, eq=a.type.theteq, n=None, phase=None, source=s,
                         at_sulfur='S' in (a.atom1.name, a.atom2.name, a.atom3.name)))
    dih = {}
    for dd in parm.dihedrals:
        if dd.improper:
            continue
        key = (dd.atom1.idx, dd.atom2.idx, dd.atom3.idx, dd.atom4.idx)
        dih.setdefault(key, []).append((dd.type.per, round(dd.type.phi_k, 5), dd.type.phase))
    for key, terms in dih.items():
        at = [parm.atoms[i] for i in key]
        t = tuple(a.type for a in at)
        s, native, ok = src('dihedral', t, terms)
        for (n, k, ph) in terms:
            rows.append(dict(kind='torsion', atoms='-'.join(a.name for a in at), types='-'.join(t), k=k, eq=None,
                             n=n, phase=ph, source=s, at_sulfur=at[1].name == 'S' or at[2].name == 'S'))
    return rows, overridden


def unique_terms(rows):
    """One row per term type at sulfur: types, parameters, number of instances, atoms, source."""
    seen = {}
    for r in rows:
        if not r['at_sulfur']:
            continue
        t = r['types'].split('-')
        canon = '-'.join(min(t, t[::-1]))
        key = (r['kind'], canon, r['n'], r['phase'], round(r['k'], 5), r['eq'])
        e = seen.setdefault(key, dict(kind=r['kind'], types=canon, k=r['k'], eq=r['eq'], n=r['n'], phase=r['phase'],
                                      instances=0, atoms=[], source=r['source']))
        e['instances'] += 1; e['atoms'].append(r['atoms'])
    out = []
    for e in seen.values():
        e = dict(e); e['atoms'] = ' '.join(e['atoms']); out.append(e)
    order = {'bond': 0, 'angle': 1, 'torsion': 2}
    return sorted(out, key=lambda e: (order[e['kind']], e['types'], e['n'] or 0))


def geometry_table(sp, ref, res, extra=None):
    """Rows of the geometry comparison for one species (res = ligand_eval.evaluate result)."""
    rows = []
    for key, label, unit in GEOM_ROWS:
        dft = ref[sp]['m'].get(key)
        if dft is None:
            continue
        v = res['vacuum']['global_min']['m'].get(key); b = res['vacuum']['from_dft_min']['m'].get(key)
        o = res['OBC2']['global_min']['m'].get(key) if 'OBC2' in res else None
        r = dict(species=sp, quantity=label, unit=unit, dft_pcm=dft,
                 dft_gas=ref['neutral_gas']['m'].get(key) if sp == 'neutral' else None,
                 crystal_6RL9_SAN304=ref['crystal']['m'].get(key),
                 mm_vacuum_global_min=v, mm_vacuum_minus_dft=(v - dft) if v is not None else None,
                 mm_vacuum_min_from_dft_start=b, mm_obc2_global_min=o,
                 mm_obc2_minus_dft=(o - dft) if o is not None else None)
        rows.append(r)
    return rows


def minima_table(sp, res):
    rows = []
    for envn, r in res.items():
        for k, g in enumerate(r['minima']):
            m = g['m']
            rows.append(dict(species=sp, environment=envn, rank=k + 1, E_rel_kcal=g['E_rel'], n_starts=g['count'],
                             tau_fold=m['tau_fold'], tau_C3=m['tau_C3'], tau_C4_S_N2_HN21=m.get('tau_C4_S_N2_HN21'),
                             tau_C4_S_N2_HN22=m.get('tau_C4_S_N2_HN22'), a_C4_S_N2=m['a_C4_S_N2'],
                             d_S_N2=m['d_S_N2'], a_O1_S_O2=m['a_O1_S_O2'], N1_angle_sum=m.get('N1_angle_sum'),
                             N2_angle_sum=m.get('N2_angle_sum'), gmax=g['gmax']))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True); BULK.mkdir(parents=True, exist_ok=True)
    print(f'ligand_params.py started {lc.now_ist()}')
    ref = reference_data()
    gaff2 = lc.read_parm_dat(lc.GAFF2_DAT)
    summary = dict(script='md/ligand_params.py', started=lc.now_ist(), versions=lc.versions(), seed=le.SEED,
                   hessian_note=ref['hessian_note'], species={})
    geom_rows, min_rows = [], []
    for sp in SPECIES:
        print(f'== {sp}')
        a = assign(sp, ref)
        print('  types:', ' '.join(f"{x['name']}:{x['type']}" for x in a['atoms']), '| total charge', a['total_charge'])
        print('  parmchk2 opened:', a['parmchk2_opened'])
        rows, overridden = term_tables(OUT / f'sfa_{sp}_gaff2.prmtop', OUT / f'sfa_{sp}_gaff2.frcmod', gaff2)
        uniq = unique_terms(rows)
        lc.write_csv(OUT / f'terms_gaff2_{sp}.csv', uniq,
                     ['kind', 'types', 'k', 'eq', 'n', 'phase', 'instances', 'atoms', 'source'])
        crow, coverridden = term_tables(pathlib.Path(a['work']) / 'control_s_gaff2.prmtop',
                                        pathlib.Path(a['work']) / 'control_s_gaff2.frcmod', gaff2)
        print(f'  terms: {len(rows)} in the topology, {len(overridden)} differ from a native gaff2.dat term; '
              f'control ("-s gaff2"): {len(coverridden)} differ')
        res = le.evaluate(OUT / f'sfa_{sp}_gaff2.prmtop', sp, ref[sp]['x'], ref[sp]['names'], ref[sp]['m'])
        # control: the scoping route, geometry only
        cmm = lc.MM(pathlib.Path(a['work']) / 'control_s_gaff2.prmtop')
        cgroups, cfrom, _ = le.multistart(cmm, ref[sp]['x'], ref[sp]['names'])
        cm = cgroups[0]['m']
        print(f"  control (-s gaff2 frcmod over gaff2.dat): C-S-N {cm['a_C4_S_N2']:.2f}, S-N {cm['d_S_N2']:.4f}")
        geom_rows += geometry_table(sp, ref, res)
        min_rows += minima_table(sp, res)
        lc.write_csv(OUT / f'scan_before_{sp}.csv', le.scan_table(res['vacuum']['scan']))
        lc.write_csv(OUT / f'scan_before_{sp}_obc2.csv', le.scan_table(res['OBC2']['scan']))
        summary['species'][sp] = dict(
            assignment=a, n_terms=len(rows), native_terms_overridden=overridden,
            terms_at_sulfur=uniq,
            control_s_gaff2=dict(what='frcmod written by "parmchk2 -s gaff2" (reads gaff.dat, GAFF 1.81) loaded '
                                      'after leaprc.gaff2',
                                 native_terms_overridden=coverridden,
                                 global_min=dict(a_C4_S_N2=cm['a_C4_S_N2'], d_S_N2=cm['d_S_N2'],
                                                 a_O1_S_O2=cm['a_O1_S_O2'], tau_fold=cm['tau_fold'])),
            dft_min=dict(file=ref[sp]['file'], m=ref[sp]['m']),
            results={e: dict(n_starts=r['n_starts'], n_minima=len(r['minima']), global_min=r['global_min']['m'],
                             global_min_xyz=np.round(r['global_min']['x'], 5).tolist(),
                             from_dft_min=dict(E_rel=r['from_dft_min']['E_rel'], m=r['from_dft_min']['m'],
                                               xyz=np.round(r['from_dft_min']['x'], 5).tolist()),
                             hessian_lowest7=r['hessian_lowest7'], scan_summary=r['scan_summary'], gate=r['gate'],
                             gate_from_dft_min=r['gate_from_dft_min']) for e, r in res.items()})
    lc.write_csv(OUT / 'geometry_before.csv', geom_rows)
    lc.write_csv(OUT / 'minima_before.csv', min_rows)
    summary['crystal_6RL9_SAN304'] = dict(m=ref['crystal']['m'], file=ref['crystal']['file'],
                                          B_iso={n: ref['crystal']['atoms'][n]['b'] for n in HEAVY})
    summary['finished'] = lc.now_ist()
    (OUT / 'gaff2_as_assigned.json').write_text(json.dumps(le.strip(summary), indent=1) + '\n')
    print(f'ligand_params.py finished {lc.now_ist()}')


if __name__ == '__main__':
    main()
