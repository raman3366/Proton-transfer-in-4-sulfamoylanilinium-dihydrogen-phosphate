"""d04: proton scan, functional check (from the log), dielectric ladder, ZPE against the barrier.

Energies are re-read from the per-point JSON files and the fsp log and converted here; distances and delta are
recomputed from each point's .xyz.  The functional-check JSON of burst 1b was lost to the DF-MP2 memory crash, so its
four DFT results are parsed from the archived log (burst1b_cont2/fsp_scan_pcm.log); the burst-2 rerun
(sfa_pair_pcm_d3/fsp_scan.json, with DF-MP2) reproduces them within 3.9e-11 Ha (SI Section S1).
"""
import re
import glob
from common import *

SC = RES / 'sfa_pair_pcm_d3/scan'
pts = []
for f in sorted(glob.glob(str(SC / '[AB]_0*.json'))):
    d = json.loads(Path(f).read_text())
    stem = Path(f).stem
    xyz = SC / f'{stem}.xyz'
    s, X, com = read_xyz(xyz)
    P = {l: X[i] for i, l in enumerate(PAIR_LABELS)}
    dNH, dHO = dist(P['N1'], P['H1A']), dist(P['H1A'], P['O5'])
    h = d.get('history') or []
    last = h[-1] if h else {}
    pts.append(dict(file=rel(f), chain=d['chain'], index=d['index'], kind=d['kind'], target=d.get('target_d_A'),
                    E_Ha=d['E_Ha'], converged=d['converged'], steps=d['steps'], accepted_on=d.get('accepted_on'),
                    d_N1_H1A=dNH, d_H1A_O5=dHO, delta=dNH - dHO, d_N1_O5=dist(P['N1'], P['O5']),
                    ang=angle(P['N1'], P['H1A'], P['O5']), d_P1_O5=dist(P['P1'], P['O5']),
                    last_gmax_conv=last.get('gmax_conv'), last_grms_conv=last.get('grms_conv'),
                    last_dE_Ha=(h[-1]['E_Ha'] - h[-2]['E_Ha']) if len(h) > 1 else None))
ref = [p for p in pts if p['chain'] == 'A' and p['index'] == 0][0]['E_Ha']
for p in pts:
    p['E_rel_kcal'] = (p['E_Ha'] - ref) * HA_KCAL
    p['E_rel_kJ'] = (p['E_Ha'] - ref) * HA_KJ
# lower envelope on the grid of d(N1-H1A) (rounded to 0.001 A)
grid = {}
for p in pts:
    k = round(p['d_N1_H1A'], 3)
    if k not in grid or p['E_Ha'] < grid[k]['E_Ha']:
        grid[k] = p
env = [grid[k] for k in sorted(grid)]
hyst = {}
for k in sorted(set(round(p['d_N1_H1A'], 3) for p in pts)):
    same = [p for p in pts if round(p['d_N1_H1A'], 3) == k]
    if len(same) == 2:
        a = [p for p in same if p['chain'] == 'A'][0]; b = [p for p in same if p['chain'] == 'B'][0]
        hyst[f'{k:.3f}'] = dict(A_minus_B_kcal=(a['E_Ha'] - b['E_Ha']) * HA_KCAL,
                                A_minus_B_kJ=(a['E_Ha'] - b['E_Ha']) * HA_KJ,
                                dNH_A=a['d_N1_H1A'], dNH_B=b['d_N1_H1A'],
                                delta_A=a['delta'], delta_B=b['delta'], dNO_A=a['d_N1_O5'], dNO_B=b['d_N1_O5'])
top = max(env, key=lambda p: p['E_rel_kcal'])
omin = min(env, key=lambda p: p['E_rel_kcal'])
out = dict(points=pts, reference_E_Ha=ref, envelope=[dict(chain=p['chain'], index=p['index'],
                                                          d=p['d_N1_H1A'], E_rel_kcal=p['E_rel_kcal']) for p in env],
           barrier_from_N_kcal=top['E_rel_kcal'], barrier_from_N_kJ=top['E_rel_kJ'], barrier_d=top['d_N1_H1A'],
           barrier_point=f"{top['chain']}_{top['index']:02d}", barrier_point_converged=top['converged'],
           O_min_rel_kcal=omin['E_rel_kcal'], O_min_rel_kJ=omin['E_rel_kJ'], O_min_d=omin['d_N1_H1A'],
           reverse_barrier_kcal=top['E_rel_kcal'] - omin['E_rel_kcal'],
           reverse_barrier_kJ=top['E_rel_kJ'] - omin['E_rel_kJ'],
           hysteresis=hyst, n_points=len(pts), n_unconverged=sum(1 for p in pts if not p['converged']))

# ---- functional check, parsed from the archived log
LOG = ARCH / 'logs/burst/burst1b_cont2/fsp_scan_pcm.log'
fsp = {}
txt = LOG.read_text()
last_scf = None
for line in txt.splitlines():
    m = re.match(r'converged SCF energy = (-?\d+\.\d+)', line)
    if m:
        last_scf = float(m.group(1))
    m = re.match(r'\[fsp\]\s+(\S+)\s+(\S+)\s+E = (-?\d+\.\d+) Ha\s+\((\d+) s\)', line)
    if m:
        meth, geom, E, t = m.group(1), m.group(2), float(m.group(3)), int(m.group(4))
        fsp.setdefault(meth, {})[geom] = dict(E_Ha=E, scf_E_Ha=last_scf, extra_disp_Ha=E - last_scf, wall_s=t,
                                              log_line=line.strip())
mp2_crash = 'MemoryError' in txt and 'Insufficient memory for holding t2 incore' in txt
fres = {}
for meth, g in fsp.items():
    eN, eT, eO = g['end_N_d1.098']['E_Ha'], g['top1_d1.250']['E_Ha'], g['end_O_d1.631']['E_Ha']
    fres[meth] = dict(dE_O_minus_N_kcal=(eO - eN) * HA_KCAL, dE_O_minus_N_kJ=(eO - eN) * HA_KJ,
                      top_minus_N_kcal=(eT - eN) * HA_KCAL, top_minus_N_kJ=(eT - eN) * HA_KJ,
                      top_minus_O_kcal=(eT - eO) * HA_KCAL, top_minus_O_kJ=(eT - eO) * HA_KJ,
                      top_is_max=bool(eT > max(eN, eO)))
out['functional_check'] = dict(source=str(LOG), status='from the burst-1b log (its JSON was lost to the DF-MP2 MemoryError); the burst-2 rerun '
                               '(fsp_scan.json) reproduces these energies within 3.9e-11 Ha',
                               geometries={'end_N_d1.098': 'scan/A_00_anchor.xyz', 'top1_d1.250': 'scan/A_02_d1.250.xyz',
                                           'end_O_d1.631': 'scan/B_00_anchor.xyz'},
                               raw=fsp, relative=fres, mp2_crashed=mp2_crash,
                               b3lyp_log_vs_scan_json_Ha=dict(
                                   N=fsp['B3LYP-D3(BJ)']['end_N_d1.098']['E_Ha'] - ref,
                                   O=fsp['B3LYP-D3(BJ)']['end_O_d1.631']['E_Ha'] - [p for p in pts if p['chain'] == 'B' and p['index'] == 0][0]['E_Ha']))

# ---- dielectric ladder (energies are not comparable across eps; only the proton position is)
lad = json.loads((RES / 'sfa_pair_pcm_d3/eps_ladder/eps_ladder.json').read_text())
out['eps_ladder'] = dict(bracket=lad['bracket_eps'], crossover_text=lad['crossover'],
                         rows=[dict(eps=r['eps'], steps=r['steps'], converged=r['converged'], E_Ha=r['E_Ha'],
                                    proton_on=r['metrics']['proton_on'], source=r['source'],
                                    started_from=r['started_from']) for r in lad['rows']])

# ---- ZPE of the N1-H1A stretch against the barrier
fr = json.loads((RES / 'sfa_pair_pcm_d3/freq.json').read_text())
zb = json.loads((RES / 'sfa_pair_pcm_d3/scan/zpe_barrier.json').read_text())
nu = zb['nh_stretch']['normal_mode']['freq_cm-1']
out['zpe'] = dict(nu_NH_cm1=nu, zpe_half_hnu_kcal=0.5 * nu * CM1_KCAL, zpe_half_hnu_kJ=0.5 * nu * CM1_KCAL * 4.184,
                  ratio_zpe_over_barrier=0.5 * nu * CM1_KCAL / top['E_rel_kcal'],
                  local_mode_cm1=zb['nh_stretch']['local_mode']['freq_cm-1'],
                  local_mode_zpe_kcal=0.5 * zb['nh_stretch']['local_mode']['freq_cm-1'] * CM1_KCAL,
                  scan_curvature_cm1=zb['scan_curvature']['freq_cm-1'],
                  one_d_E0_kcal=zb['one_d_model']['E0_above_reference_kcal'],
                  one_d_E0_minus_top_kcal=zb['one_d_model']['E0_minus_barrier_top_kcal'],
                  one_d_P_O_side=zb['one_d_model']['levels'][0]['P_O_side'],
                  barrier_spline_kcal=zb['profile']['barrier_cubic_spline_kcal'],
                  barrier_spline_d=zb['profile']['d_barrier_cubic_spline_A'],
                  effective_mass_amu=zb['effective_mass_amu'],
                  nh_character=zb['nh_stretch']['normal_mode']['nh_stretch_character'],
                  mode_index=zb['nh_stretch']['normal_mode']['index'],
                  freq_json_has_mode=bool(min(abs(np.array(fr['frequencies_cm-1']) - nu)) < 1e-6))
dump('d04_scan.json', out)
