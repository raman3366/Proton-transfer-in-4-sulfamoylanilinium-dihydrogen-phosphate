"""Test that tleap loads the protein and the kept crystal water of each start (stage S3, part 1).

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/check_tleap.py

For each start (X, N, R) tleap (AmberTools 26; leaprc.protein.ff14SB, leaprc.water.tip3p) loads
md/system/protein_md.pdb and md/system/waters_<start>.pdb.  The zinc and the ligand are left out: their
parameters are made in stages S1 and S2.  The test passes when

  - tleap ends without an error and writes a topology,
  - the heavy atoms tleap adds are exactly those listed in decisions.json as not modelled in the entry
    (plus nothing: OXT and the ACE cap are in the PDB file),
  - no atom name of the PDB file is unknown to the residue templates,
  - the heavy atoms of the input keep their coordinates,
  - the residues have the hydrogens decided (His tautomers, Glu106 without a carboxyl hydrogen,
    Thr199 with HG1, Cys206 with HG, no charged N-terminus, C-terminal carboxylate).

Small results go to md/system/tleap/ (leap inputs and logs, protein_H.pdb = the protein with the
hydrogens tleap adds, tleap_check.json); the test topologies go to /root/cqc-sfa-md/system/tleap/.
Nothing is random: tleap builds missing atoms from the internal coordinates of its templates.
"""
import json
import os
import pathlib
import re
import subprocess
import sys

import numpy as np
import parmed

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / 'md' / 'system'
SMALL = OUT / 'tleap'
BULK = pathlib.Path('/root/cqc-sfa-md/system/tleap')
AMBERHOME = os.environ.get('AMBERHOME', '/root/cqc-sfa-tools/md')
TLEAP = pathlib.Path(AMBERHOME) / 'bin' / 'tleap'
THREADS = '2'


def run_leap(tag, waters):
    """waters: file name or None (protein only)."""
    inp = ['source leaprc.protein.ff14SB', 'source leaprc.water.tip3p',
           f'prot = loadpdb {OUT / "protein_md.pdb"}']
    if waters:
        inp += [f'wat = loadpdb {OUT / waters}', 'sys = combine { prot wat }']
    else:
        inp += ['sys = prot']
    inp += ['check sys', 'charge sys',
            f'saveamberparm sys {BULK / (tag + ".prmtop")} {BULK / (tag + ".inpcrd")}',
            f'savepdb sys {BULK / (tag + ".pdb")}', 'quit']
    (SMALL / f'leap_{tag}.in').write_text('\n'.join(inp) + '\n')
    env = dict(os.environ, AMBERHOME=AMBERHOME, OMP_NUM_THREADS=THREADS, OPENBLAS_NUM_THREADS=THREADS,
               MKL_NUM_THREADS=THREADS)
    for f in (BULK / f'{tag}.prmtop', BULK / f'{tag}.inpcrd', BULK / f'{tag}.pdb'):
        if f.exists():
            f.unlink()
    p = subprocess.run(['nice', '-n', '10', str(TLEAP), '-f', str(SMALL / f'leap_{tag}.in')], cwd=BULK, env=env,
                       capture_output=True, text=True)
    log = p.stdout + p.stderr
    # absolute paths of this machine are of no use in the stored log
    (SMALL / f'leap_{tag}.log').write_text(log.replace(str(REPO) + '/', '').replace(str(BULK) + '/', ''))
    (BULK / 'leap.log').unlink(missing_ok=True)
    return p.returncode, log


def parse_log(log):
    added = re.findall(r'Added missing heavy atom: \.R<(\S+) (\d+)>\.A<(\S+) \d+>', log)
    created = re.findall(r'Created a new atom named: (\S+) within residue: \.R<(\S+) (\d+)>', log)
    summary = re.search(r'Leap added (\d+) missing atoms? according to residue templates:\s*(\d+) Heavy\s*(\d+) H', log)
    total = re.search(r'Total unperturbed charge:\s*(-?[\d.]+)', log)
    exiting = re.search(r'Exiting LEaP: Errors = (\d+); Warnings = (\d+); Notes = (\d+)', log)
    close = re.findall(r'Close contact of ([\d.]+) angstroms between nonbonded atoms (\S+) and (\S+)', log)
    heavy_close = [(d, a, b) for d, a, b in close if not (a.startswith('H') and b.startswith('H'))]
    return dict(added_heavy=[(n, int(i), a) for n, i, a in added], created_unknown=created,
                leap_summary=dict(total=int(summary.group(1)), heavy=int(summary.group(2)), hydrogen=int(summary.group(3))) if summary else None,
                total_charge=float(total.group(1)) if total else None,
                errors=int(exiting.group(1)) if exiting else None, warnings=int(exiting.group(2)) if exiting else None,
                notes=int(exiting.group(3)) if exiting else None,
                unknown_residue='Unknown residue' in log, fatal=('FATAL' in log or 'Fatal' in log),
                close_contacts_reported=len(close),
                close_contacts_between_hydrogens=len(close) - len(heavy_close),
                close_contacts_with_a_heavy_atom=[dict(dist_A=float(d), atoms=[a, b]) for d, a, b in heavy_close],
                shortest_close_contact_A=min((float(d) for d, _, _ in close), default=None),
                other_warnings=sorted({ln.strip() for ln in log.splitlines()
                                       if 'unperturbed charge' in ln or 'Unknown' in ln or 'Could not' in ln}))


def input_heavy():
    out = []
    for ln in (OUT / 'protein_md.pdb').read_text().splitlines():
        if ln.startswith('ATOM'):
            out.append((ln[17:20].strip(), int(ln[22:26]), ln[12:16].strip(),
                        np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return out


def angle(a, b, c):
    v1, v2 = a - b, c - b
    return float(np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2))))


def analyse(tag, dec, rmap):
    st = parmed.load_file(str(BULK / f'{tag}.prmtop'), str(BULK / f'{tag}.inpcrd'))
    xyz = st.coordinates
    by_index = {r.idx + 1: r for r in st.residues}
    idx_of = {int(row['deposited_number']): int(row['tleap_index']) for row in rmap}
    dep_of = {v: k for k, v in idx_of.items()}

    def res(dep):
        return by_index[idx_of[dep]]

    def pos(dep, name):
        return xyz[[a for a in res(dep).atoms if a.name == name][0].idx]

    prot_res = [r for r in st.residues if r.name not in ('WAT', 'HOH')]
    wat_res = [r for r in st.residues if r.name in ('WAT', 'HOH')]
    names = {}
    for r in st.residues:
        names[r.name] = names.get(r.name, 0) + 1
    # heavy atoms of the input must be where they were
    inp = input_heavy()
    key = {(dep_of[r.idx + 1], a.name): a.idx for r in prot_res for a in r.atoms if a.atomic_number > 1}
    disp = [float(np.linalg.norm(xyz[key[(n, nm)]] - p)) for _, n, nm, p in inp]
    in_keys = {(n, nm) for _, n, nm, _ in inp}
    added = sorted(k for k in key if k not in in_keys)
    # tleap-built heavy atoms: nearest other heavy atom that is not in the same residue
    heavy_idx = np.array([a.idx for a in st.atoms if a.atomic_number > 1])
    heavy_res = np.array([st.atoms[i].residue.idx for i in heavy_idx])
    built = []
    for n, nm in added:
        i = key[(n, nm)]
        d = np.linalg.norm(xyz[heavy_idx] - xyz[i], axis=1)
        d[heavy_res == st.atoms[i].residue.idx] = 99
        j = int(d.argmin())
        o = st.atoms[int(heavy_idx[j])]
        on = 'HOH' if o.residue.name == 'WAT' else f'{o.residue.name}{dep_of.get(o.residue.idx + 1, "")}'
        built.append(dict(atom=f'{res(n).name}{n}.{nm}', nearest=f'{on}.{o.name}' + (f' (residue index {o.residue.idx + 1})' if on == 'HOH' else ''),
                          dist_A=round(float(d[j]), 2)))
    # hydrogens of the residues decided
    hyd = {}
    for dep, (state, _d, _r) in ((int(k[3:]), (v['state'], None, None)) for k, v in dec['protein']['histidines'].items()):
        r = res(dep)
        hn = sorted(a.name for a in r.atoms if a.name in ('HD1', 'HE2'))
        want = {'HID': ['HD1'], 'HIE': ['HE2'], 'HIP': ['HD1', 'HE2']}[state]
        hyd[f'HIS{dep}'] = dict(residue_name=r.name, ring_NH=hn, as_decided=bool(r.name == state and hn == want))
    g = res(106)
    hyd['GLU106'] = dict(residue_name=g.name, carboxyl_H=[a.name for a in g.atoms if a.name in ('HE1', 'HE2')],
                         as_decided=bool(g.name == 'GLU' and not [a for a in g.atoms if a.name in ('HE1', 'HE2')]))
    g = res(117)
    hyd['GLU117'] = dict(residue_name=g.name, as_decided=bool(g.name == 'GLU' and not [a for a in g.atoms if a.name in ('HE1', 'HE2')]))
    t = res(199)
    hg1, og1, oe1, oe2 = pos(199, 'HG1'), pos(199, 'OG1'), pos(106, 'OE1'), pos(106, 'OE2')
    hyd['THR199'] = dict(residue_name=t.name, hydrogens=sorted(a.name for a in t.atoms if a.atomic_number == 1),
                         as_decided=bool('HG1' in [a.name for a in t.atoms]),
                         HG1_as_placed_by_tleap=dict(HG1_OE1_A=round(float(np.linalg.norm(hg1 - oe1)), 2),
                                                     OG1_HG1_OE1_deg=round(angle(og1, hg1, oe1), 1),
                                                     HG1_OE2_A=round(float(np.linalg.norm(hg1 - oe2)), 2),
                                                     note='tleap places hydroxyl hydrogens from the template, not by hydrogen bonds; the entry has HG1 1.81 A from OE1 (162 deg)'))
    c = res(206)
    hyd['CYS206'] = dict(residue_name=c.name, as_decided=bool(c.name == 'CYS' and 'HG' in [a.name for a in c.atoms]))
    y = res(7)
    hyd['TYR7'] = dict(residue_name=y.name, as_decided=bool('HH' in [a.name for a in y.atoms]))
    first, h3, last = by_index[1], res(3), res(261)
    hyd['N_terminus'] = dict(first_residue=first.name, atoms=sorted(a.name for a in first.atoms),
                             His3_N_hydrogens=sorted(a.name for a in h3.atoms if a.name in ('H', 'H1', 'H2', 'H3')),
                             as_decided=bool(first.name == 'ACE' and sorted(a.name for a in h3.atoms if a.name in ('H', 'H1', 'H2', 'H3')) == ['H']))
    hyd['C_terminus'] = dict(last_residue=last.name, has_OXT=bool('OXT' in [a.name for a in last.atoms]),
                             as_decided=bool('OXT' in [a.name for a in last.atoms]))
    q = np.array([a.charge for a in st.atoms])
    pq = float(sum(a.charge for r in prot_res for a in r.atoms))
    return dict(
        atoms=len(st.atoms), residues=len(st.residues), protein_residues=len(prot_res), waters=len(wat_res),
        protein_atoms=sum(len(r.atoms) for r in prot_res), protein_heavy_atoms=sum(1 for r in prot_res for a in r.atoms if a.atomic_number > 1),
        water_atoms=sum(len(r.atoms) for r in wat_res),
        net_charge=round(float(q.sum()), 4), protein_charge=round(pq, 4),
        residue_names=dict(sorted(names.items())),
        input_heavy_atoms=len(inp), max_displacement_of_input_heavy_atoms_A=round(max(disp), 4),
        heavy_atoms_added_by_tleap=[f'{res(n).name}{n}.{nm}' for n, nm in added],
        built_atom_contacts=built, hydrogens=hyd,
    )


def main():
    SMALL.mkdir(parents=True, exist_ok=True)
    BULK.mkdir(parents=True, exist_ok=True)
    dec = json.loads((OUT / 'decisions.json').read_text())
    import csv
    rmap = list(csv.DictReader(open(OUT / 'residue_map.csv')))
    expected = sorted(dec['protein']['expected_tleap_added_heavy_atoms'])
    meta = sorted(p.name[:-5] for p in (pathlib.Path(AMBERHOME) / 'conda-meta').glob('ambertools-*.json'))
    result = dict(written_by='md/check_tleap.py', tleap=str(TLEAP), amberhome=AMBERHOME,
                  force_field=['leaprc.protein.ff14SB', 'leaprc.water.tip3p'],
                  software=dict(ambertools=meta[0] if meta else None, parmed=parmed.__version__, numpy=np.__version__,
                                python=sys.version.split()[0]),
                  not_in_this_test='zinc and ligand (no parameters yet); no solvent box, no ions',
                  expected_added_heavy_atoms=expected, runs={})
    ok_all = True
    for tag, wat in (('protein', None), ('X', 'waters_X.pdb'), ('N', 'waters_N.pdb'), ('R', 'waters_R.pdb')):
        rc, log = run_leap(tag, wat)
        info = parse_log(log)
        wrote = (BULK / f'{tag}.prmtop').exists() and (BULK / f'{tag}.prmtop').stat().st_size > 0
        entry = dict(returncode=rc, topology_written=bool(wrote), log=info)
        if wrote:
            a = analyse(tag, dec, rmap)
            entry.update(a)
            got = sorted(x.replace('CLYS', 'LYS').replace('NLYS', 'LYS') for x in a['heavy_atoms_added_by_tleap'])
            got = sorted(re.sub(r'^(HID|HIE|HIP)', 'HIS', x) for x in got)
            entry['added_heavy_atoms_as_expected'] = bool(got == expected and info['leap_summary'] and info['leap_summary']['heavy'] == len(expected))
            entry['hydrogens_as_decided'] = bool(all(v['as_decided'] for v in a['hydrogens'].values()))
            passed = (rc == 0 and info['errors'] == 0 and not info['unknown_residue'] and not info['fatal']
                      and not info['created_unknown'] and not info['close_contacts_with_a_heavy_atom']
                      and entry['added_heavy_atoms_as_expected']
                      and entry['hydrogens_as_decided'] and a['max_displacement_of_input_heavy_atoms_A'] < 0.002)
        else:
            passed = False
        entry['passed'] = bool(passed)
        ok_all = ok_all and passed
        result['runs'][tag] = entry
        print(f'{tag:8s} tleap rc {rc}, errors {info["errors"]}, warnings {info["warnings"]}; '
              + (f'{entry["atoms"]} atoms, {entry["waters"]} waters, charge {entry["net_charge"]:+.3f}; '
                 f'added heavy {len(entry["heavy_atoms_added_by_tleap"])} (expected {len(expected)}); '
                 f'passed {passed}' if wrote else 'NO TOPOLOGY'))
    result['gate_tleap_load_ok'] = bool(ok_all)
    # the protein with the hydrogens tleap adds, deposited residue numbers restored, for the next stages
    if (BULK / 'protein.pdb').exists():
        idx2dep = {int(r['tleap_index']): int(r['deposited_number']) for r in rmap}
        out = ['REMARK   protein_md.pdb after tleap (ff14SB): hydrogens and the side-chain atoms without density added',
               'REMARK   residue numbers of 6RL9 restored; hydroxyl and water hydrogens are NOT hydrogen-bond optimised']
        for ln in (BULK / 'protein.pdb').read_text().splitlines():
            if ln.startswith(('ATOM', 'HETATM')):
                out.append(f'{ln[:21]}A{idx2dep[int(ln[22:26])]:4d}{ln[26:]}')
            elif ln.startswith('TER'):
                out.append('TER')
        out.append('END')
        (SMALL / 'protein_H.pdb').write_text('\n'.join(out) + '\n')
    (SMALL / 'tleap_check.json').write_text(json.dumps(result, indent=1) + '\n')
    print('gate "tleap load ok":', 'PASSED' if ok_all else 'FAILED')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
