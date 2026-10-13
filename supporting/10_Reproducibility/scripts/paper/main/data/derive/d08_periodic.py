"""d08: periodic PBE-D3(BJ) model -- energies, proton position, all nine crystal H-bonds, RMSD to the X-ray.

Reads the relaxed cells (periodic_qe/relax_{i,ii,iifix,iirfix}.extxyz) and the X-ray start (start_i_xray.extxyz)
with our own code (relax_iirfix, the constrained neutral with H1A and O5 held in space, added 29 Sep 2026).
H-bond images are fixed on the X-ray start (the acceptor image that reproduces the CIF H...A) and reused for the
relaxed cells.  Heavy-atom RMSD is computed after removing the mean shift (Pc has a floating
origin in the a-c plane), and also raw.  Energies from periodic_qe.json (Ry) converted here.
"""
import itertools
from common import *

PD = RES / 'periodic_qe'
x01 = json.loads((OUT / 'd01_xray.json').read_text())
pj = json.loads((PD / 'periodic_qe.json').read_text())


def load(name):
    s, X, lab, cell, head = read_extxyz(PD / f'{name}.extxyz')
    return dict(zip(lab, X)), cell, lab, s, head


x0, cell, labels, syms0, _ = load('start_i_xray')
A = cell  # rows are lattice vectors
imgs = [np.array(v) @ A for v in itertools.product(range(-2, 3), repeat=3)]


def acc_label(A_, sym):
    return A_ if sym.startswith('1') else A_ + '_2'


# fix the acceptor image on the X-ray start
hbdef = []
for h in x01['hbonds']:
    D, H_, A_, sym = h['D'], h['H'], h['A'], h['sym']
    al = acc_label(A_, sym)
    cda, cha = float(h['cif']['DA'].split('(')[0]), float(h['cif']['HA'].split('(')[0])
    tH = min(imgs, key=lambda t: np.linalg.norm(x0[H_] + t - x0[D]))          # H next to its donor
    tA = min(imgs, key=lambda t: abs(np.linalg.norm(x0[al] + t - x0[D]) - cda)
             + abs(np.linalg.norm(x0[al] + t - x0[H_] - tH) - cha))            # the acceptor image of the CIF
    hbdef.append(dict(name=f"{D}-{H_}...{A_}", D=D, H=H_, A=al, sym=sym, shiftH=tH, shift=tA,
                      cif_DA=cda, cif_HA=cha, cif_DHA=float(h['cif']['DHA'].split('(')[0])))

res = {}
for name in ['start_i_xray', 'relax_i', 'relax_ii', 'relax_iifix', 'relax_iirfix', 'start_ii_OH']:
    P, cell_, lab, s, head = load(name)
    maxdisp = max(np.linalg.norm(P[l] - x0[l]) for l in lab)
    rows = []
    for d in hbdef:
        dD, dH, dA = P[d['D']], P[d['H']] + d['shiftH'], P[d['A']] + d['shift']
        rows.append(dict(name=d['name'], sym=d['sym'], DH=dist(dD, dH), HA=dist(dH, dA), DA=dist(dD, dA),
                         DHA=angle(dD, dH, dA), cif_DA=d['cif_DA'], dDA=dist(dD, dA) - d['cif_DA']))
    heavy = [l for l, e in zip(lab, s) if e != 'H']
    D = np.array([P[l] - x0[l] for l in heavy])
    shift = D.mean(0)
    rmsd_raw = float(np.sqrt((D ** 2).sum(1).mean()))
    rmsd = float(np.sqrt(((D - shift) ** 2).sum(1).mean()))
    DA_ = np.array([P[l] - x0[l] for l in lab])
    sh_all = DA_.mean(0)
    rmsd_all = float(np.sqrt(((DA_ - sh_all) ** 2).sum(1).mean()))
    kmax = int(np.argmax(np.linalg.norm(D - shift, axis=1)))
    res[name] = dict(max_atom_displacement_from_xray_A=float(maxdisp), hbonds=rows,
                     rmsd_heavy_shift_removed=rmsd, rmsd_heavy_raw=rmsd_raw, mean_shift=shift.tolist(),
                     max_dev_heavy=float(np.linalg.norm(D - shift, axis=1).max()), max_dev_atom=heavy[kmax],
                     rmsd_all_shift_removed=rmsd_all, n_heavy=len(heavy), n_all=len(lab),
                     H1A_position_change_from_start_ii=None)
    # proton metric for both formula units
    for fu, suf in [(1, ''), (2, '_2')]:
        n1, h1, o5 = P['N1' + suf], P['H1A' + suf], None
        # O5 partner: the image nearest to H1A
        o5 = min((P['O5' + suf] + t for t in imgs), key=lambda q: np.linalg.norm(q - h1))
        res[name][f'fu{fu}'] = dict(d_N1_H1A=dist(n1, h1), d_H1A_O5=dist(h1, o5), delta=dist(n1, h1) - dist(h1, o5),
                                    d_N1_O5=dist(n1, o5), ang=angle(n1, h1, o5))
        o4, h4 = P['O4' + suf], P['H4A' + suf]
        o5b = min((P['O5' + suf] + t for t in imgs), key=lambda q: np.linalg.norm(q - h4))
        res[name][f'fu{fu}_O4H4AO5'] = dict(OH=dist(o4, h4), HO=dist(h4, o5b), OO=dist(o4, o5b),
                                            ang=angle(o4, h4, o5b), frac=float((h4 - o4) @ (o5b - o4) / dist(o4, o5b) ** 2))

# the constrained 'neutral' start: where did H1A end up relative to where it was fixed?
Pii, *_ = load('start_ii_OH')
Pfx, *_ = load('relax_iifix')
res['relax_iifix']['H1A_moved_from_fixed_position_A'] = dist(Pii['H1A'], Pfx['H1A'])
res['relax_iifix']['N1_moved_A'] = dist(Pii['N1'], Pfx['N1'])
res['relax_iifix']['O5_moved_A'] = dist(Pii['O5'], Pfx['O5'])
heavy = [l for l, e in zip(labels, syms0) if e != 'H']
res['relax_iifix']['mean_heavy_shift_vs_start_ii'] = np.mean([Pfx[l] - Pii[l] for l in heavy], axis=0).tolist()

# the constrained neutral that held (burst 2): relax_i with H1A moved 0.98 A from O5, then H1A and O5 held in space
Pri, *_ = load('relax_i')
Prr, *_ = load('relax_iirfix')
sh_r = np.mean([Prr[l] - Pri[l] for l in heavy], axis=0)
res['relax_iirfix']['mean_heavy_shift_vs_relax_i'] = sh_r.tolist()
res['relax_iirfix']['mean_heavy_shift_vs_relax_i_norm'] = float(np.linalg.norm(sh_r))
res['relax_iirfix']['N1_moved_vs_relax_i_A'] = dist(Pri['N1'], Prr['N1'])
res['relax_iirfix']['rmsd_heavy_vs_relax_i_shift_removed'] = float(np.sqrt(np.mean(
    [np.sum((Prr[l] - Pri[l] - sh_r) ** 2) for l in heavy])))
for fu, suf in [(1, ''), (2, '_2')]:
    o5 = Prr['O5' + suf]
    res['relax_iirfix'][f'fu{fu}']['P1_O5'] = min(dist(Prr['P1' + suf], o5 + t) for t in imgs)
    res['relax_iirfix'][f'fu{fu}']['C1_N1'] = min(dist(Prr['C1' + suf], Prr['N1' + suf] + t) for t in imgs)

# energies (Ry) -> per formula unit
E = {k: pj['relax'][k]['energy_Ry'] for k in ['i', 'ii', 'iifix', 'iirfix']}
Z = 2
en = dict(E_Ry=E,
          dE_ii_minus_i_meV_cell=(E['ii'] - E['i']) * RY_EV * 1000,
          dE_ii_minus_i_kJ_fu=(E['ii'] - E['i']) * RY_EV * EV_KJ / Z,
          dE_iifix_minus_i_meV_cell=(E['iifix'] - E['i']) * RY_EV * 1000,
          dE_iifix_minus_i_kJ_fu=(E['iifix'] - E['i']) * RY_EV * EV_KJ / Z,
          dE_iifix_minus_i_kcal_fu=(E['iifix'] - E['i']) * RY_EV * EV_KJ / Z / 4.184,
          dE_iirfix_minus_i_meV_cell=(E['iirfix'] - E['i']) * RY_EV * 1000,
          dE_iirfix_minus_i_kJ_fu=(E['iirfix'] - E['i']) * RY_EV * EV_KJ / Z,
          json_comparison=pj['comparison'],
          max_force=dict((k, pj['relax'][k]['max_force_Ry_bohr']) for k in E),
          pressure_kbar=dict((k, pj['relax'][k]['pressure_kbar_final']) for k in E),
          bfgs_steps=dict((k, pj['relax'][k]['bfgs_steps']) for k in E),
          ecut_Ry=pj['relax']['i']['ecut_Ry'], k=pj['relax']['i']['k'], qe_version=pj['relax']['i']['qe_version'],
          json_rmsd_heavy=pj['relax']['i']['rmsd_heavy_vs_xray'])
conv = json.loads((PD / 'convergence.json').read_text())
# hybrid_estimate.json is the record of the cost estimate made before the hybrid runs; its own status string calls it an
# unverified order-of-magnitude estimate. The hybrid single points replaced it (periodic_qe.json -> hybrid, profile.json
# -> hybrid), so it is labelled here as what it is (ACS integration review, 30 Sep 2026; values unchanged).
hyb = json.loads((PD / 'hybrid_estimate.json').read_text())
hyb['status'] = ('pre-run order-of-magnitude estimate of the cost of the hybrid single points; replaced by the '
                 'hybrid single points themselves (periodic_qe.json hybrid, profile.json hybrid)')
out = dict(structures=res, energies=en, convergence=conv, method=pj['method'], hybrid=hyb)
dump('d08_periodic.json', out)
