"""Stage S1 of md/PLAN.md, side table: how the fit of the S-N rotor depends on the environment and on the number
of Fourier terms.

  /root/cqc-sfa-tools/md/bin/python md/ligand_rotor_models.py          (seconds)

Same data and same least squares as ligand_fit.py (rigid DFT single points, sigma equal for all; coefficients of
cos(n phi) on H-N-S-C4, one coefficient for the aniline rotor), but for vacuum and implicit water and for n = (1, 2)
and n = (1, 2, 3).  The C4-S torsion is GAFF2's here, as in the first cycle of ligand_fit.py.  Documents two
choices of ligand_fit.py: the environment (the n = 1 coefficient changes by several kcal/mol between vacuum and
implicit water) and the two-term series (a third term does not lower the leave-one-out error for both species).
                                                                              -> md/ligand/rotor_models.json
"""
import json
import numpy as np
import ligand_common as lc
import ligand_fit as lf
from ligand_common import OUT, SPECIES


def main():
    out = dict(script='md/ligand_rotor_models.py', versions=lc.versions(), rows=[])
    for sp in SPECIES:
        names, rigid, _ = lf.rigid_points(sp)
        xs = [p['x'] for p in rigid]
        d = np.array([p['E'] for p in rigid])
        for envn in ('vacuum', 'OBC2'):
            model = lf.Model(sp, OUT / f'sfa_{sp}_gaff2.prmtop', None if envn == 'vacuum' else envn)
            model.zero_original('rotor')
            er = model.e_rest(xs, ['rotor'])
            y = d[1:] - (er[1:] - er[0])
            for ns in ((1, 2), (1, 2, 3)):
                B = np.array([[model.basis(x, 'rotor', n) for n in ns] + [model.aniline_basis(x)] for x in xs])
                A = (B - B[0])[1:]
                c = np.linalg.lstsq(A, y, rcond=None)[0]
                loo = []
                for k in range(len(y)):
                    m = np.ones(len(y), bool); m[k] = False
                    ck = np.linalg.lstsq(A[m], y[m], rcond=None)[0]
                    loo.append(float(A[k] @ ck - y[k]))
                row = dict(species=sp, environment=envn, n_terms=list(ns),
                           coefficients_of_cos={str(n): float(v) for n, v in zip(ns, c)},
                           aniline_cos2=float(c[-1]), rms=float(np.sqrt(((A @ c - y) ** 2).mean())),
                           largest=float(np.abs(A @ c - y).max()), leave_one_out_rms=float(np.sqrt(np.mean(np.square(loo)))),
                           n_points=len(y))
                out['rows'].append(row)
                print(f"{sp:8s} {envn:7s} n={ns}: coefficients {[round(float(v), 2) for v in c[:-1]]}, aniline "
                      f"{c[-1]:.2f}, rms {row['rms']:.2f}, leave-one-out {row['leave_one_out_rms']:.2f}")
    (OUT / 'rotor_models.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
