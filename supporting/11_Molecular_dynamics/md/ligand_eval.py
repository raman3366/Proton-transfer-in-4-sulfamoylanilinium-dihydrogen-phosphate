"""Evaluation of one ligand topology against the DFT reference (stage S1 of md/PLAN.md): the minima found from many
starts, and the relaxed scan of C3-C4-S-N2 with the protocol of the DFT scan.  Used unchanged before and after the
repair, by ligand_params.py and ligand_fit.py.

Protocol of the scan (the DFT one, dft/torsion_scan.py): every point is independent; it starts from the species'
unconstrained minimum with the sulfur side (S, O1, O2, N2 and its hydrogens) rotated rigidly about C4-S to the
target; C3-C4-S-N2 is held, every other coordinate relaxes.  Here the start is (a) the MM minimum reached from the
DFT minimum and (b) the DFT geometry of the same point; the lower energy is kept and the difference is recorded.
A third set of starts (all rotations about S-N and both pyramids of the aniline NH2) looks for a lower basin that
the DFT scan, which follows one basin, would not have seen; it is reported beside the scan, not used for the gate.
"""
import math
import numpy as np
from ligand_common import (MM, SPECIES, dft_scan, measures, add_means, rotate_about, set_scan_dihedral, wrap,
                           dihedral)

SEED = 20261005


def _aniline_flip(x, ix):
    """Mirror the two aniline hydrogens through the ring plane (the other pyramid at N1)."""
    y = np.array(x, float)
    n = np.cross(x[ix['C2']] - x[ix['C1']], x[ix['C6']] - x[ix['C1']]); n /= np.linalg.norm(n)
    for h in ('HN11', 'HN12'):
        v = y[ix[h]] - x[ix['C1']]
        y[ix[h]] = y[ix[h]] - 2.0 * (v @ n) * n
    return y


def rotor_starts(x, names, sn_steps=(0, 60, 120, 180, 240, 300), flips=(False, True)):
    """Starts that differ in the rotation of the hydrogens of N2 about S-N and in the pyramid at N1."""
    ix = {n: i for i, n in enumerate(names)}
    hs = [ix[h] for h in ('HN21', 'HN22') if h in ix]
    out = []
    for fl in flips:
        y = _aniline_flip(x, ix) if fl else np.array(x, float)
        for a in sn_steps:
            out.append(rotate_about(y, ix['S'], ix['N2'], hs, a))
    return out


def multistart(mm, x_dft, names, seed=SEED, jitter=0.02):
    """Unconstrained minimisations from the DFT minimum and from 144 rotated starts (12 x C-S, 6 x S-N, 2 x N1
    pyramid, Gaussian jitter of 0.02 A, fixed seed).  Returns the list of distinct minima, lowest first, and the
    minimum reached from the DFT geometry itself."""
    rng = np.random.default_rng(seed)
    ix = {n: i for i, n in enumerate(names)}
    tau0 = dihedral(x_dft, ix['C3'], ix['C4'], ix['S'], ix['N2'])
    runs = []
    r0 = mm.minimise(x_dft)
    runs.append(dict(start='DFT minimum', **r0))
    for dcs in range(0, 360, 30):
        y = set_scan_dihedral(x_dft, names, wrap(tau0 + dcs))
        for k, z in enumerate(rotor_starts(y, names)):
            z = z + rng.normal(0.0, jitter, z.shape)
            r = mm.minimise(z)
            runs.append(dict(start=f'C-S +{dcs}, rotor start {k}', **r))
    for r in runs:
        r['m'] = add_means(measures(r['x'], names))
    emin = min(r['E'] for r in runs)
    groups = []
    for r in sorted(runs, key=lambda r: r['E']):
        for g in groups:
            if abs(r['E'] - g['E']) < 2e-3 and abs(r['m']['tau_fold'] - g['m']['tau_fold']) < 3.0 \
                    and abs(r['m']['a_C4_S_N2'] - g['m']['a_C4_S_N2']) < 0.3:
                g['count'] += 1
                break
        else:
            groups.append(dict(E=r['E'], E_rel=r['E'] - emin, count=1, m=r['m'], x=r['x'], gmax=r['gmax'],
                               start=r['start']))
    return groups, runs[0], len(runs)


def scan(mm, species, x_min_basin, names, global_search=True):
    """Relaxed scan at the torsions of the DFT scan (C3 branch).  Returns one dict per point."""
    pts = [p for p in dft_scan(species) if p['kind'] == 'constrained' and p['branch'] == 'C3']
    rows = []
    for p in pts:
        tau = p['m']['tau_C3']                         # the value the DFT point really has (180 is stored as -180)
        a = mm.minimise(set_scan_dihedral(x_min_basin, names, tau), tau=tau)
        b = mm.minimise(p['x'], tau=tau)
        best = a if a['E'] <= b['E'] else b
        row = dict(tau=abs(p['target']) if p['target'] else tau, tau_C3=tau, E_dft_rel_ref=p['E_rel'],
                   E_mm=best['E'], E_mm_from_min=a['E'], E_mm_from_dft_geom=b['E'],
                   start_diff=abs(a['E'] - b['E']), tau_dev=best['tau_dev'], gmax=best['gmax'], x=best['x'],
                   m=add_means(measures(best['x'], names)), m_dft=add_means(p['m']))
        if global_search:
            lo = best
            for z in rotor_starts(set_scan_dihedral(x_min_basin, names, tau), names):
                c = mm.minimise(z, tau=tau)
                if c['E'] < lo['E']:
                    lo = c
            row['E_mm_global'] = lo['E']
            row['m_global'] = add_means(measures(lo['x'], names))
        rows.append(row)
    e_mm0 = min(r['E_mm'] for r in rows)
    e_dft0 = min(r['E_dft_rel_ref'] for r in rows)
    for r in rows:
        r['E_mm_rel'] = r['E_mm'] - e_mm0
        r['E_dft_rel'] = r['E_dft_rel_ref'] - e_dft0
        r['dev'] = r['E_mm_rel'] - r['E_dft_rel']
        if global_search:
            r['E_mm_global_rel'] = r['E_mm_global'] - e_mm0
    return rows


def scan_summary(rows):
    dev = np.array([r['dev'] for r in rows])
    i_mm = int(np.argmin([r['E_mm_rel'] for r in rows])); i_dft = int(np.argmin([r['E_dft_rel'] for r in rows]))
    j_mm = int(np.argmax([r['E_mm_rel'] for r in rows])); j_dft = int(np.argmax([r['E_dft_rel'] for r in rows]))
    s = dict(n_points=len(rows), max_abs_dev=float(np.abs(dev).max()), tau_of_max_dev=rows[int(np.abs(dev).argmax())]['tau'],
             rms_dev=float(np.sqrt((dev ** 2).mean())),
             mm_min_at=rows[i_mm]['tau'], dft_min_at=rows[i_dft]['tau'],
             mm_max_at=rows[j_mm]['tau'], dft_max_at=rows[j_dft]['tau'],
             mm_barrier=rows[j_mm]['E_mm_rel'], dft_barrier=rows[j_dft]['E_dft_rel'],
             max_start_diff=float(max(r['start_diff'] for r in rows)),
             max_tau_dev=float(max(abs(r['tau_dev']) for r in rows)))
    if 'E_mm_global' in rows[0]:
        s['max_gain_other_basin'] = float(max(r['E_mm'] - r['E_mm_global'] for r in rows))
    return s


def scan_table(rows):
    out = []
    for r in rows:
        o = dict(tau_C3_C4_S_N2=round(r['tau'], 2), E_dft_rel=r['E_dft_rel'], E_mm_rel=r['E_mm_rel'], dev=r['dev'],
                 E_mm_global_rel=r.get('E_mm_global_rel'),
                 csn_dft=r['m_dft']['a_C4_S_N2'], csn_mm=r['m']['a_C4_S_N2'],
                 sn_dft=r['m_dft']['d_S_N2'], sn_mm=r['m']['d_S_N2'],
                 tauH21_dft=r['m_dft'].get('tau_C4_S_N2_HN21'), tauH21_mm=r['m'].get('tau_C4_S_N2_HN21'),
                 tauH22_dft=r['m_dft'].get('tau_C4_S_N2_HN22'), tauH22_mm=r['m'].get('tau_C4_S_N2_HN22'),
                 tauH21_mm_global=(r.get('m_global') or {}).get('tau_C4_S_N2_HN21'),
                 start_diff=r['start_diff'], tau_dev=r['tau_dev'], gmax=r['gmax'])
        out.append(o)
    return out


def gate(geom_min, m_dft, ssum):
    """Gate G-L of md/PLAN.md: at the MM minimum C-S-N within 2 deg and S-N within 0.02 A of DFT; torsion profile
    within 0.5 kcal/mol of the DFT scan (largest deviation over the scan, each profile relative to its minimum)."""
    d_csn = geom_min['a_C4_S_N2'] - m_dft['a_C4_S_N2']
    d_sn = geom_min['d_S_N2'] - m_dft['d_S_N2']
    g = dict(csn_mm=geom_min['a_C4_S_N2'], csn_dft=m_dft['a_C4_S_N2'], d_csn=d_csn, csn_pass=bool(abs(d_csn) <= 2.0),
             sn_mm=geom_min['d_S_N2'], sn_dft=m_dft['d_S_N2'], d_sn=d_sn, sn_pass=bool(abs(d_sn) <= 0.02),
             scan_max_abs_dev=ssum['max_abs_dev'], scan_pass=bool(ssum['max_abs_dev'] <= 0.5))
    g['pass'] = bool(g['csn_pass'] and g['sn_pass'] and g['scan_pass'])
    return g


def evaluate(prmtop, species, x_dft, names, m_dft, environments=('vacuum', 'OBC2'), log=print):
    """Minima, scan and gate of one topology in vacuum (the gate) and in implicit water (OBC2 generalised Born,
    dielectric 78.3553 as the IEF-PCM of the DFT; reported as a check of the sensitivity to the environment)."""
    out = {}
    for envn in environments:
        mm = MM(prmtop, implicit=None if envn == 'vacuum' else envn)
        assert mm.names == names, (mm.names, names)
        groups, from_dft, nruns = multistart(mm, x_dft, names)
        gmin = groups[0]
        w = mm.hessian_min_eig(gmin['x'])
        # the minimum of the basin the DFT scan follows: reached from the DFT minimum
        basin = dict(E=from_dft['E'], E_rel=from_dft['E'] - gmin['E'], m=from_dft['m'], x=from_dft['x'],
                     gmax=from_dft['gmax'])
        rows = scan(mm, species, basin['x'], names, global_search=(envn == 'vacuum'))
        ssum = scan_summary(rows)
        res = dict(n_starts=nruns, minima=groups, global_min=gmin, from_dft_min=basin,
                   hessian_lowest7=[float(v) for v in w[:7]], scan=rows, scan_summary=ssum,
                   gate=gate(gmin['m'], m_dft, ssum), gate_from_dft_min=gate(basin['m'], m_dft, ssum))
        out[envn] = res
        g = res['gate']
        log(f"  [{species} {envn}] {nruns} starts, {len(groups)} distinct minima; global minimum: "
            f"C-S-N {g['csn_mm']:.2f} (DFT {g['csn_dft']:.2f}), S-N {g['sn_mm']:.4f} (DFT {g['sn_dft']:.4f}), "
            f"tau_fold {gmin['m']['tau_fold']:.1f}; scan max|dev| {ssum['max_abs_dev']:.3f} kcal/mol at "
            f"{ssum['tau_of_max_dev']:.1f}; gate {'PASS' if g['pass'] else 'FAIL'}")
    return out


def strip(res):
    """JSON-ready copy of an evaluate() result (coordinates rounded, arrays to lists)."""
    def conv(o):
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        if isinstance(o, np.ndarray):
            return np.round(o, 5).tolist()
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o
    return conv(res)
