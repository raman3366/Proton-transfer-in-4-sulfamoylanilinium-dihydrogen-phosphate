#!/usr/bin/env python
"""Gate G-R2, for the record (no weight in the verdict): THROWAWAY builds of the protein systems of start X and
start R with a candidate file of terms at N2, into the examiner's scratch folder, and a comparison of each built
topology with the system as built (/root/cqc-sfa-md/systems/bonded_2e/<start>) by a parser written here.

    AMBERHOME=/root/cqc-sfa-tools/md nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_build.py \
        --candidate second_step|charges_as_built [--starts X R] [--reuse]

md/build_topology.py is imported unchanged and its function build() is called with out_root and records_root in
scratch (it removes the files of a system before it builds: the folders of the real systems are never passed).
The checksums of the real systems are compared with md/topology/bonded_2e/<start>/SHA256SUMS before and after.
Writes md/qm/n2/final_test/built_system_<candidate>.json.  The built systems are NOT ACCEPTED terms: never to be simulated
for a result.
"""
import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
N2 = REPO / 'md' / 'qm' / 'n2'
OUT = N2 / 'final_test'
SCR = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                   '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/n2/final')
REAL = pathlib.Path('/root/cqc-sfa-md/systems/bonded_2e')
CAND = dict(second_step=dict(frcmod=N2 / 'n2_terms_new.frcmod', mcpb_dir=N2 / 'charges2' / 'mcpb_2e_two'),
            charges_as_built=dict(frcmod=N2 / 'n2_terms_new_charges_as_built.frcmod', mcpb_dir=None))


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def real_checksums(start):
    p = subprocess.run(['sha256sum', '-c', str(REPO / 'md' / 'topology' / 'bonded_2e' / start / 'SHA256SUMS')],
                       cwd=REAL / start, capture_output=True, text=True)
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    return dict(ok=sum(l.endswith(': OK') for l in lines), of=len(lines), exit=p.returncode)


# ----------------------------------------------------------------------------- an AMBER topology, read here
def prmtop(path):
    sec, name, fmt = {}, None, None
    for l in pathlib.Path(path).read_text().splitlines():
        if l.startswith('%FLAG'):
            name = l.split()[1]
            sec[name] = []
        elif l.startswith('%FORMAT'):
            fmt = l[l.index('(') + 1:l.index(')')]
            sec[name] = dict(fmt=fmt, text=[])
        elif name and isinstance(sec.get(name), dict):
            sec[name]['text'].append(l)
    out = {}
    for k, v in sec.items():
        if not isinstance(v, dict):
            continue
        f = v['fmt'].upper()
        w = int(''.join(ch for ch in f.split('A' if 'A' in f else 'I' if 'I' in f else 'E')[1].split('.')[0] if ch.isdigit()))
        items = []
        for l in v['text']:
            items += [l[i:i + w] for i in range(0, len(l.rstrip('\n')), w)]
        if 'A' in f:
            out[k] = [s.strip() for s in items]
        elif 'I' in f:
            out[k] = [int(s) for s in items if s.strip()]
        else:
            out[k] = [float(s) for s in items if s.strip()]
        out['_text_' + k] = '\n'.join(v['text'])
    return out


def terms(t):
    """Atoms, bonds, angles, torsions (proper, improper), 1-4 pairs of a topology, keyed by atom index."""
    n = t['POINTERS'][0]
    res = [0] * n
    rp = t['RESIDUE_POINTER'] + [n + 1]
    for r in range(len(rp) - 1):
        for a in range(rp[r] - 1, rp[r + 1] - 1):
            res[a] = r
    name = [f"{t['RESIDUE_LABEL'][res[a]]}{res[a] + 1}:{t['ATOM_NAME'][a]}" for a in range(n)]
    out = dict(n=n, name=name, charge=[round(q / 18.2223, 6) for q in t['CHARGE']], type=t['AMBER_ATOM_TYPE'],
               lj_type=t['ATOM_TYPE_INDEX'], mass=t['MASS'], bonds={}, angles={}, proper={}, improper={}, pairs14={})
    for key in ('BONDS_INC_HYDROGEN', 'BONDS_WITHOUT_HYDROGEN'):
        v = t[key]
        for i in range(0, len(v), 3):
            a, b, k = v[i] // 3, v[i + 1] // 3, v[i + 2] - 1
            out['bonds'][tuple(sorted((a, b)))] = (t['BOND_FORCE_CONSTANT'][k], t['BOND_EQUIL_VALUE'][k])
    for key in ('ANGLES_INC_HYDROGEN', 'ANGLES_WITHOUT_HYDROGEN'):
        v = t[key]
        for i in range(0, len(v), 4):
            a, b, c, k = v[i] // 3, v[i + 1] // 3, v[i + 2] // 3, v[i + 3] - 1
            out['angles'][(a, b, c) if a < c else (c, b, a)] = (t['ANGLE_FORCE_CONSTANT'][k], round(float(np.degrees(t['ANGLE_EQUIL_VALUE'][k])), 4))
    for key in ('DIHEDRALS_INC_HYDROGEN', 'DIHEDRALS_WITHOUT_HYDROGEN'):
        v = t[key]
        for i in range(0, len(v), 5):
            a, b, c, d, k = v[i] // 3, v[i + 1] // 3, abs(v[i + 2]) // 3, abs(v[i + 3]) // 3, v[i + 4] - 1
            term = (int(round(t['DIHEDRAL_PERIODICITY'][k])), round(t['DIHEDRAL_FORCE_CONSTANT'][k], 6),
                    round(float(np.degrees(t['DIHEDRAL_PHASE'][k])) % 360.0, 3))
            if v[i + 3] < 0:
                out['improper'].setdefault((a, b, c, d), []).append(term)
            else:
                out['proper'].setdefault((a, b, c, d) if a < d else (d, c, b, a), []).append(term)
                if v[i + 2] >= 0:
                    out['pairs14'][tuple(sorted((a, d)))] = (round(t['SCEE_SCALE_FACTOR'][k], 4), round(t['SCNB_SCALE_FACTOR'][k], 4))
    for d in (out['proper'], out['improper']):
        for v in d.values():
            v.sort()
    return out


def compare(ref, new):
    nm = lambda k: '-'.join(new['name'][i] for i in k)                    # noqa: E731
    d = dict(atoms=ref['n'] == new['n'] and ref['name'] == new['name'] and ref['type'] == new['type'] and
             ref['lj_type'] == new['lj_type'] and ref['mass'] == new['mass'])
    ch = [i for i in range(new['n']) if ref['charge'][i] != new['charge'][i]]
    d['charges_changed'] = dict(n=len(ch), residues=sorted({new['name'][i].split(':')[0] for i in ch}),
                                largest=max([abs(ref['charge'][i] - new['charge'][i]) for i in ch], default=0.0),
                                total_before=round(sum(ref['charge']), 5), total_after=round(sum(new['charge']), 5))
    for sec in ('bonds', 'angles', 'proper', 'improper', 'pairs14'):
        keys = sorted(set(ref[sec]) | set(new[sec]))
        diff = [k for k in keys if ref[sec].get(k) != new[sec].get(k)]
        d[sec] = dict(n_before=len(ref[sec]), n_after=len(new[sec]), n_different=len(diff),
                      different=[dict(atoms=nm(k), before=ref[sec].get(k), after=new[sec].get(k)) for k in diff[:40]])
    return d


def allowed(d, new):
    """Which differences are terms at N2 as point 2 of the plan allows; everything else is listed."""
    bad = []
    lig = [n for n in new['name'] if n.endswith(':N2')][0].split(':')[0]
    zn = [n for n in new['name'] if n.endswith(':ZN')][0]
    n2, s, h = f'{lig}:N2', f'{lig}:S', f'{lig}:HN21'
    if not d['atoms']:
        bad.append('atoms, types, Lennard-Jones types or masses differ')
    for sec in ('bonds', 'pairs14'):
        if d[sec]['n_different']:
            bad.append(f"{d[sec]['n_different']} {sec} differ")
    for x in d['angles']['different']:
        a = x['atoms'].split('-')
        if not (a[1] == n2 and {a[0], a[2]} in ({zn, s}, {zn, h}, {s, h})):
            bad.append('angle ' + x['atoms'])
    for x in d['proper']['different']:
        a = x['atoms'].split('-')
        mid, ends = {a[1], a[2]}, {a[0], a[3]}
        ok = (mid == {n2, s} and ends & {h, zn}) or (mid == {zn, n2} and ends & {s, h})
        if not ok:
            bad.append('torsion ' + x['atoms'])
    for x in d['improper']['different']:
        if set(x['atoms'].split('-')) != {zn, n2, s, h}:
            bad.append('improper ' + x['atoms'])
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate', required=True, choices=sorted(CAND))
    ap.add_argument('--starts', nargs='+', default=['X', 'R'])
    ap.add_argument('--reuse', action='store_true')
    args = ap.parse_args()
    c = CAND[args.candidate]
    work = SCR / f'build_{args.candidate}'
    assert str(work).startswith('/tmp/claude-0/') and '/root/cqc-sfa-md' not in str(work) and REPO not in work.parents
    sys.path.insert(0, str(REPO / 'md'))
    import build_topology as bt
    doc = dict(what='THROWAWAY builds of the protein systems with a candidate file of terms at N2 (NOT ACCEPTED), and their comparison '
                    'with the systems as built; for the record of gate G-R2, no weight in the verdict', written=now(), candidate=args.candidate,
               frcmod=str(c['frcmod']), frcmod_sha256=sha(c['frcmod']), mcpb_dir=str(c['mcpb_dir']) if c['mcpb_dir'] else 'default of the builder (charges as built)',
               builder=dict(script='md/build_topology.py', sha256=sha(REPO / 'md' / 'build_topology.py')), work=str(work), starts={})
    for start in args.starts:
        before = real_checksums(start)
        new = work / 'systems' / 'bonded_2e' / start
        rec = dict(real_system_checksums_before=before)
        if not (args.reuse and (new / 'system.prmtop').exists()):
            try:
                bj = bt.build('bonded', start, str(work / 'systems'), mcpb_dir=str(c['mcpb_dir'].resolve()) if c['mcpb_dir'] else None,
                              records_root=str(work / 'records'), pset='2e', n2_frcmod=str(c['frcmod'].resolve()),
                              note='THROWAWAY build of the examiner of gate G-R2 (md/qm/n2/final_test/n2_exam_build.py); terms NOT ACCEPTED; never to be simulated for a result')
                rec.update(checks_failed=bj['checks_failed'], n_checks=len(bj['checks']), checks_passed=sum(v is True for v in bj['checks'].values()))
            except bt.BuildError as e:
                rec.update(build_error=str(e))
                doc['starts'][start] = rec
                continue
        rec['real_system_checksums_after'] = real_checksums(start)
        t_ref, t_new = prmtop(REAL / start / 'system.prmtop'), prmtop(new / 'system.prmtop')
        flags = sorted(k for k in t_new if not k.startswith('_'))
        rec['sections'] = dict(n=len(flags), identical_text=sum(t_ref.get('_text_' + k) == t_new['_text_' + k] for k in flags),
                               different=[k for k in flags if t_ref.get('_text_' + k) != t_new['_text_' + k]])
        a, b = terms(t_ref), terms(t_new)
        d = compare(a, b)
        rec['comparison'] = d
        rec['differences_the_plan_does_not_allow'] = allowed(d, b)
        rec['coordinates_identical'] = sha(REAL / start / 'system.inpcrd') == sha(new / 'system.inpcrd')
        rec['prmtop_sha256'] = sha(new / 'system.prmtop')
        rec['system'] = str(new)
        doc['starts'][start] = rec
        print(start, 'atoms', b['n'], '| sections identical', rec['sections']['identical_text'], 'of', rec['sections']['n'], '| different:', rec['sections']['different'])
        print('   charges changed:', d['charges_changed'], '| angles', d['angles']['n_different'], '| proper', d['proper']['n_different'],
              '| improper', d['improper']['n_different'], '| bonds', d['bonds']['n_different'], '| 1-4 pairs', d['pairs14']['n_different'])
        print('   not allowed:', rec['differences_the_plan_does_not_allow'], '| builder checks failed:', rec.get('checks_failed'),
              '| real system', before, rec['real_system_checksums_after'], '| coordinates identical', rec['coordinates_identical'])
    (OUT / f'built_system_{args.candidate}.json').write_text(json.dumps(doc, indent=1) + '\n')


if __name__ == '__main__':
    main()
