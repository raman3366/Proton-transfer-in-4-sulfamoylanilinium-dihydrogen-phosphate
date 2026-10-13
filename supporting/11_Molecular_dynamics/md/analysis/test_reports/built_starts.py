# throwaway (versions 1.2 and 1.3): the three built systems of the restrained-zinc model as runs of one frame, analysed by
# analyse_md.py, against the same numbers from plain formulas in this file (own reader of the Amber files, minimum
# image by a search over the 27 neighbouring cells).  Checks the reading of the real selections.json and
# reference.pdb (40,900 atoms) and the new quantities on real coordinates.  No simulation: the one frame of each
# "run" is system.inpcrd as built.   usage: built_starts.py SCRATCH ANALYSIS_DIR
import csv, hashlib, json, math, subprocess, sys
from pathlib import Path
import numpy as np
S, A = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(A.parent))
import run_md                                              # the driver's own writers of DCD and PDB files
import openmm.app as app
import openmm.unit as u
SYS = Path('/root/cqc-sfa-md/systems/restrained')
PY = sys.executable
out = {'what': 'the built systems (system.inpcrd, restrained-zinc model) as one-frame runs: analyse_md.py against plain formulas',
       'systems': str(SYS), 'analysis_script_sha256': hashlib.sha256((A / 'analyse_md.py').read_bytes()).hexdigest(),
       'driver_script_sha256': hashlib.sha256((A.parent / 'run_md.py').read_bytes()).hexdigest(), 'starts': {}}
worst = {'A': 0.0, 'vol': 0.0}
ok_all = True
for s in 'XNR':
    sysd = SYS / s
    top = app.AmberPrmtopFile(str(sysd / 'system.prmtop'))
    crd = app.AmberInpcrdFile(str(sysd / 'system.inpcrd'))
    pos = np.array(crd.positions.value_in_unit(u.nanometer))
    box = np.array([[v.x, v.y, v.z] for v in crd.boxVectors.value_in_unit(u.nanometer)])
    sel = json.loads((sysd / 'selections.json').read_text())
    an, rn, roa = run_md.prmtop_names(str(sysd / 'system.prmtop'))      # names as tleap wrote them (HID, HIE, WAT), as the driver reads them
    el = [a.element.symbol if a.element is not None else 'X' for a in top.topology.atoms()]
    atoms = [(i, an[i], rn[roa[i]], int(roa[i]), el[i]) for i in range(len(an))]
    rd = S / 'built' / 'runs' / 'restrained' / f'{s}1'
    rd.mkdir(parents=True, exist_ok=True)
    solute = sel['solute']
    (rd / 'full.pdb').write_text(run_md.pdb_text(atoms, pos, box, ['built start as a one-frame run (test)']))
    (rd / 'solute.pdb').write_text(run_md.pdb_text([atoms[i] for i in solute], pos[solute], box, ['built start as a one-frame run (test)']))
    for name, idx, every in (('solute_seg0001.dcd', solute, 2500), ('full_seg0001.dcd', list(range(len(atoms))), 2500)):
        w = run_md.DcdWriter(str(rd / name), len(idx), every, every, 0.004, ['built start, one frame', 'test'])
        w.write(pos[idx], box)
        w.close()
    (rd / 'manifest.json').write_text(json.dumps({'model': 'restrained', 'start': s, 'replica': 1, 'planned_ns': 0.01,
                                                  'system_dir': str(sysd), 'full_interval_ps': 10.0,
                                                  'note': 'one frame: system.inpcrd as built; not a simulation'}) + '\n')
    r = subprocess.run(['nice', '-n', '10', PY, str(A / 'analyse_md.py'), str(rd), '--no-figures'], capture_output=True, text=True)
    summ = json.loads((rd / 'analysis' / 'summary.json').read_text())
    ts = list(csv.DictReader(open(rd / 'analysis' / 'timeseries_solute.csv')))[0]
    tf = list(csv.DictReader(open(rd / 'analysis' / 'timeseries_full.csv')))[0]

    # ---- plain formulas (Angstrom)
    x = pos * 10.0
    B = box * 10.0
    cells = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float) @ B

    def mi(a, b):                                          # minimum-image distances from the points a [m, 3] to b [n, 3]
        d = b[None, :, None, :] - a[:, None, None, :] + cells[None, None, :, :]
        return np.sqrt((d ** 2).sum(-1)).min(-1)
    L = {k: x[i] for k, i in sel['lig_heavy'].items()}
    zn, h = x[sel['zn']], x[sel['lig_HN2']]
    e = [(v - L['N2']) / math.sqrt(float((v - L['N2']) @ (v - L['N2']))) for v in (zn, h, L['S'])]
    vol = float(np.linalg.det(np.array(e)))
    ref = [l for l in open(sysd / 'reference.pdb') if l.startswith(('ATOM', 'HETATM'))]
    rxyz = {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in ref if l.startswith('HETATM')}
    rca = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ref if l[12:16].strip() == 'CA'])
    ca = x[sel['protein_CA']]
    cx, cr = ca.mean(0), rca.mean(0)
    U, _, Vt = np.linalg.svd((ca - cx).T @ (rca - cr))
    Q = U @ np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))]) @ Vt        # frame -> reference: (y - cx) Q + cr
    wat = x[sel['water_O']]
    # version 1.3: the count that the driver acts on (six ring carbons, 2.5 A; these one-frame manifests name no
    # settings, so the analysis uses its defaults) and the wider count of version 1.2 (seven positions, 3.5 A)
    site = np.array([(rxyz[a] - cr) @ Q.T + cx for a in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6')])
    dsite = mi(site, wat).min(0)
    wide = np.array([(rxyz[a] - cr) @ Q.T + cx for a in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1')])
    dwide = mi(wide, wat).min(0)
    fit_lig = np.array([(L[a] - cx) @ Q + cr for a in sorted(L)])
    rmsd = math.sqrt(float(((fit_lig - np.array([rxyz[a] for a in sorted(L)])) ** 2).sum(1).mean()))
    plain = {'n2_volume': vol, 'n2_sign': int(np.sign(vol)),
             'd_Zn_E106OE1_A': float(np.linalg.norm(zn - x[sel['glu106_OE1']])), 'd_Zn_E106OE2_A': float(np.linalg.norm(zn - x[sel['glu106_OE2']])),
             'n_water_ring_place': int((dsite <= 2.5).sum()), 'd_water_ring_place_min_A': float(dsite.min()),
             'n_water_ring_wide_3p5A': int((dwide <= 3.5).sum()), 'd_water_ring_wide_min_A': float(dwide.min()),
             'n_water_N1_3p5A': int((mi(L['N1'][None], wat)[0] <= 3.5).sum()), 'd_Zn_Owater_min_A': float(mi(zn[None], wat)[0].min()),
             'lig_rmsd_identity_A': rmsd, 'd1_N1_F131CZ_A': float(np.linalg.norm(L['N1'] - x[sel['phe131_CZ']])),
             'd2_N1_N67OD1_A': float(np.linalg.norm(L['N1'] - x[sel['asn67_OD1']]))}
    got = {k: (float(ts[k]) if k in ts else float(tf[k])) for k in plain}
    rows = []
    for k, v in plain.items():
        kind = 'count' if k.startswith('n_') or k == 'n2_sign' else 'vol' if k == 'n2_volume' else 'A'
        d = abs(v - got[k])
        ok = d == 0 if kind == 'count' else d <= 0.001
        if kind != 'count':
            worst[kind] = max(worst[kind], d)
        ok_all &= ok
        rows.append({'quantity': k, 'plain_formulas': round(v, 6), 'analyse_md': got[k], 'abs_diff': round(d, 6), 'agrees': bool(ok)})
    out['starts'][s] = {'n_atoms': len(atoms), 'n_solute': len(solute), 'n_water_O': len(sel['water_O']), 'analysis_exit_code': r.returncode,
                        'warnings': summ['warnings'], 'state': ts['state'], 'role': summ['run']['role'],
                        'box_A_deg': [round(float(v), 4) for v in run_md.box_lengths_angles(box)[:3]] and
                                     [round(10 * float(run_md.box_lengths_angles(box)[i]), 4) for i in range(3)] + [round(float(run_md.box_lengths_angles(box)[i]), 4) for i in range(3, 6)],
                        'rows': rows,
                        'water_ring_place': summ['water_ring_place'],
                        'three_nearest_water_O_to_the_crystal_positions_of_the_six_ring_carbons_A': [round(float(v), 3) for v in np.sort(dsite)[:3]],
                        'three_nearest_water_O_to_the_seven_positions_of_the_wider_count_A': [round(float(v), 3) for v in np.sort(dwide)[:3]],
                        'water_O_within_2p5_3p0_3p5_A_of_the_six_ring_carbon_positions': [int((dsite <= c).sum()) for c in (2.5, 3.0, 3.5)]}
    ok_all &= r.returncode == 0
    print(s, 'exit', r.returncode, 'warnings', len(summ['warnings']), {k: round(v, 4) for k, v in got.items()})
out.update(tolerance_A=0.001, largest_abs_diff_A=round(worst['A'], 6), largest_abs_diff_signed_volume=round(worst['vol'], 6), all_agree=bool(ok_all))
(S / 'built' / 'built_starts.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'starts'}))
