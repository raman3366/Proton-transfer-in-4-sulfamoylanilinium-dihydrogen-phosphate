import json, hashlib
from pathlib import Path
R = Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/md')
BULK = Path('/root/cqc-sfa-md')
B = {s: json.load(open(R/'topology'/('bonded_'+s)/'X'/'build.json')) for s in ('2e','2s')}
sb = json.load(open(R/'production'/'seminario_blocks.json'))
fz = sb['frcmod_zinc_lines']
def key(t):
    p = t.split('-'); return '-'.join(p if p[0] <= p[-1] else p[::-1])
F = {s: {key(k): v for k, v in fz[s].items()} for s in ('2e','2s','2ms')}
out = []; P = out.append
z = {s: B[s]['zinc_site_bonded'] for s in B}
P('| Bond | Types | In 6RL9 (Å) | 2e: k | 2e: r0 (Å) | 2s: k | 2s: r0 (Å) | 2ms: k | 2ms: r0 (Å) |')
P('|---|---|---|---|---|---|---|---|---|')
be = {b['bond']: b for b in z['2e']['bonds_to_zinc_against_the_parameter_files']}
order = ['Zn-SN1304:N2', 'Zn-HD194:NE2', 'Zn-HD296:NE2', 'Zn-HE1119:ND1']
bs = {b['bond']: b for b in z['2s']['bonds_to_zinc_against_the_parameter_files']}
nm = lambda s: s.replace('HD194', 'His94').replace('HD296', 'His96').replace('HE1119', 'His119').replace('SN1304:', '').replace(':', ' ')
for k in order:
    b, e = bs[k], be[k]; m = F['2ms'][key(b['types'])]
    assert abs(F['2s'][key(b['types'])]['k'] - b['k_kcal_per_mol_A2']) < 1e-9 and abs(F['2e'][key(e['types'])]['k'] - e['k_kcal_per_mol_A2']) < 1e-9
    P(f"| {nm(k)} | {b['types']} | {b['distance_in_6RL9_A']:.4f} | {e['k_kcal_per_mol_A2']:.1f} | {e['r0_A']:.4f} | {b['k_kcal_per_mol_A2']:.1f} | {b['r0_A']:.4f} | {m['k']:.1f} | {m['equilibrium']:.4f} |")
P('')
P('| Angle | Types | In 6RL9 (°) | 2e: k | 2e: θ0 (°) | 2s: k | 2s: θ0 (°) | 2s θ0 − 6RL9 (°) | 2ms: k |')
P('|---|---|---|---|---|---|---|---|---|')
ae = {a['angle']: a for a in z['2e']['angles_with_zinc_against_the_parameter_files']}
for a in z['2s']['angles_with_zinc_against_the_parameter_files']:
    if a['kind'] not in ('N-Zn-N', 'Zn-N-C'): continue
    e = ae[a['angle']]; x = a['angle_in_this_start_deg']; m = F['2ms'][key(a['types'])]
    assert abs(F['2s'][key(a['types'])]['k'] - a['k_kcal_per_mol_rad2']) < 1e-9 and abs(F['2e'][key(e['types'])]['equilibrium'] - e['theta0_deg']) < 1e-9
    assert abs(m['equilibrium'] - a['theta0_deg']) < 1e-9
    P(f"| {nm(a['angle'])} | {a['types']} | {x:.2f} | {e['k_kcal_per_mol_rad2']:.2f} | {e['theta0_deg']:.2f} | {a['k_kcal_per_mol_rad2']:.2f} | {a['theta0_deg']:.2f} | {a['theta0_deg']-x:+.2f} | {m['k']:.2f} |")
P('')
P('Lines of the MCPB.py files at N2 that are NOT in force (replaced by the file of the terms at N2):\n')
P('| Types | 2e: k, θ0 | 2s: k, θ0 | 2ms: k, θ0 |'); P('|---|---|---|---|')
for t in ('M1-Y4-hn', 'M1-Y4-sy'):
    P(f"| {t} | {F['2e'][key(t)]['k']:.2f}, {F['2e'][key(t)]['equilibrium']:.2f} | {F['2s'][key(t)]['k']:.2f}, {F['2s'][key(t)]['equilibrium']:.2f} | {F['2ms'][key(t)]['k']:.2f}, {F['2ms'][key(t)]['equilibrium']:.2f} |")
P('')
# torsions through zinc, improper
for s in ('2e','2s'):
    t = z[s]['torsions_through_zinc']
    P(f"set {s}: torsions through the zinc {t['n']}, largest |k| {t['largest_abs_k_kcal_per_mol']}, impropers with the zinc {len(t['impropers_with_zinc'])}; angles by kind {z[s]['n_angles_by_kind']}")
P('')
# N2 terms
n2 = B['2s']['terms_at_N2']
P('N2 file sha256 ' + n2['sha256'] + ' | 2e: ' + B['2e']['terms_at_N2']['sha256'])
for e in n2['in_the_topology']['entries_of_the_file']:
    P('  ' + json.dumps(e))
P('')
# charges
ch = z['2s']['charges']; che = z['2e']['charges']
P('| Residue | Sum of charges (e) | Largest difference from the unconstrained fit (e) |'); P('|---|---|---|')
for rn, rv in ch['residues'].items():
    assert che['residues'][rn]['atoms_e'] == rv['atoms_e']
    P(f"| {rv['residue']} | {rv['sum_e']:+.6f} | {rv['largest_difference_from_the_unconstrained_fit_e']:.6f} |")
P(f"| all five | {ch['site_sum_e']:+.4f} | |")
P('')
r = ch['residues']
P('| Atom | Charge (e) |'); P('|---|---|')
P(f"| Zn | {r['ZN1']['atoms_e']['ZN']:+.6f} |")
P(f"| NE2 of His94 | {r['HD1']['atoms_e']['NE2']:+.6f} |")
P(f"| NE2 of His96 | {r['HD2']['atoms_e']['NE2']:+.6f} |")
P(f"| ND1 of His119 | {r['HE1']['atoms_e']['ND1']:+.6f} |")
for a, q in r['SN1']['atoms_e'].items():
    P(f"| {a} of the ligand | {q:+.6f} |")
P('')
P('charges of 2e and 2s systems identical: ' + str(all(che['residues'][k]['atoms_e'] == ch['residues'][k]['atoms_e'] for k in ch['residues'])) + '; resp json sha ' + ch['resp_equiv_json_sha256'])
P('')
# systems
P('| System | Start | Atoms | Water | Na+ / Cl− | Checks true / n.a. / false | system.prmtop sha256 | system.inpcrd sha256 | builder sha256 |')
P('|---|---|---|---|---|---|---|---|---|')
for model in ('bonded_2e', 'bonded_2s', 'bonded_2e_n2fit', 'bonded_2s_n2fit'):
    for st in ('X', 'N', 'R'):
        p = R/'topology'/model/st/'build.json'
        if not p.exists(): continue
        b = json.load(open(p)); c = b['checks']
        tr = sum(1 for v in c.values() if v is True); na = sum(1 for v in c.values() if v is None); fa = sum(1 for v in c.values() if v is False)
        o = b['outputs_sha256']
        sysd = BULK/'systems'/model/st
        for fn in ('system.prmtop', 'system.inpcrd'):
            h = hashlib.sha256(open(sysd/fn, 'rb').read()).hexdigest()
            ov = o[fn] if isinstance(o[fn], str) else o[fn].get('sha256')
            assert h == ov, (model, st, fn)
        g = lambda fn: (o[fn] if isinstance(o[fn], str) else o[fn]['sha256'])
        P(f"| {model} | {st} | {b['composition']['atoms']:,} | {b['composition']['water_molecules']:,} | {b['composition']['Na+']} / {b['composition']['Cl-']} | {tr} / {na} / {fa} | `{g('system.prmtop')}` | `{g('system.inpcrd')}` | `{b['script_sha256'][:16]}…` |")
P('')
for model in ('bonded_2s', 'bonded_2e_n2fit', 'bonded_2s_n2fit'):
    for st in ('X', 'N', 'R'):
        p = R/'topology'/model/st/'build.json'
        if not p.exists(): continue
        c = json.load(open(p))['compared_with']
        P(f"{model}/{st} against {c['other_model']}: differences by class {c['differences_by_class']}; classes that differ but may not {c['classes_that_differ_but_may_not']}; 1-4 pairs {c['pairs_1_4']}; atoms differ {sum(c['atoms_that_differ_in'].values())}; coordinates {c['largest_coordinate_difference_A']}; prmtop sections that differ {c['prmtop_sections']['that_differ']}; only allowed {c['only_the_allowed_terms_differ']}")
print('\n'.join(out))
