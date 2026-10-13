"""Stage S1 of md/PLAN.md, side table: the rotation about S-N of the free ligand in the two force fields.

  /root/cqc-sfa-tools/md/bin/python md/ligand_sn_rotor.py          (about 3 minutes)

Relaxed MM scan of C4-S-N2-HN21 in steps of 15 degrees (implicit water, OBC2; every other coordinate relaxes; each
point from three starts: the neighbouring point, the DFT minimum with the hydrogens turned rigidly, and the same
with the ring at the other side), for GAFF2 as assigned and for the repaired terms.  No DFT curve exists for this
coordinate; the DFT evidence is the minimum itself and the single points of md/ligand/rigid_single_points.csv.
The table shows how far the two parameter sets differ along the coordinate that, in the zinc-bound form, turns the
ring from one orientation into the other.          -> md/ligand/sn_rotor_profile.csv
"""
import numpy as np
import ligand_common as lc
from ligand_common import OUT, SPECIES


def main():
    rows = []
    for sp in SPECIES:
        el, x0, _ = lc.read_xyz(lc.RESULTS / SPECIES[sp]['case'] / 'opt.xyz')
        names = lc.names_for(sp, el, x0)
        ix = {n: i for i, n in enumerate(names)}
        hs = [ix[h] for h in ('HN21', 'HN22') if h in ix]
        phi0 = lc.dihedral(x0, ix['C4'], ix['S'], ix['N2'], ix['HN21'])
        prof = {}
        for tag in ('gaff2', 'repaired'):
            mm = lc.MM(OUT / f'sfa_{sp}_{tag}.prmtop', scan_atoms=('C4', 'S', 'N2', 'HN21'), implicit='OBC2')
            prev = None
            res = {}
            for phi in list(range(-180, 180, 15)):
                starts = [lc.rotate_about(x0, ix['S'], ix['N2'], hs, phi - phi0)]
                starts.append(lc.set_scan_dihedral(starts[0], names, -lc.dihedral(x0, ix['C3'], ix['C4'], ix['S'], ix['N2'])))
                if prev is not None:
                    starts.append(lc.rotate_about(prev, ix['S'], ix['N2'], hs,
                                                  phi - lc.dihedral(prev, ix['C4'], ix['S'], ix['N2'], ix['HN21'])))
                best = min((mm.minimise(s, tau=float(phi)) for s in starts), key=lambda r: r['E'])
                prev = best['x']
                res[phi] = best
            e0 = min(r['E'] for r in res.values())
            prof[tag] = {phi: (r['E'] - e0, lc.measures(r['x'], names)) for phi, r in res.items()}
        for phi in sorted(prof['gaff2']):
            g, r = prof['gaff2'][phi], prof['repaired'][phi]
            rows.append(dict(species=sp, tau_C4_S_N2_HN21=float(phi), E_gaff2_as_assigned=g[0], E_repaired=r[0],
                             tau_fold_gaff2=g[1]['tau_fold'], tau_fold_repaired=r[1]['tau_fold'],
                             csn_gaff2=g[1]['a_C4_S_N2'], csn_repaired=r[1]['a_C4_S_N2']))
    lc.write_csv(OUT / 'sn_rotor_profile.csv', rows)
    for sp in SPECIES:
        sub = [r for r in rows if r['species'] == sp]
        print(sp, 'phi: as assigned / repaired (kcal/mol above each minimum)')
        print('  ' + '  '.join(f"{r['tau_C4_S_N2_HN21']:.0f}: {r['E_gaff2_as_assigned']:.1f}/{r['E_repaired']:.1f}" for r in sub))


if __name__ == '__main__':
    main()
