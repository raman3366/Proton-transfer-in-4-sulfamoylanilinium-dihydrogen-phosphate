"""Why the calculations use pure (5D, 7F) functions: conditioning of the density-fitting metric.

    OMP_NUM_THREADS=1 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/metric_conditioning.py

For the small and the large model, with Cartesian and with pure functions: number of orbital and fitting
functions (def2-universal-JKFIT), whether the Cholesky factorisation of the Coulomb metric (aux|aux)
succeeds, and for the small model the extreme eigenvalues.  Result: md/qm/metric_conditioning.json.
"""
import json
import pathlib
import sys

import numpy as np
import scipy.linalg

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import pyscf_calc as pc  # noqa: E402
from pyscf import df  # noqa: E402

out = []
for tag, f in (('small', 'CA2_small_opt.com'), ('large', 'CA2_large_mk.com')):
    com = gaucom.read_com(HERE / 'mcpb' / f)
    for cart in (True, False):
        s = pc.Settings(verbose=0, cart=cart)
        mol = pc.build_mol(com.symbols, com.coords, com.charge, com.multiplicity, s)
        aux = df.addons.make_auxmol(mol, s.auxbasis)
        j2c = aux.intor('int2c2e', hermi=1)
        try:
            scipy.linalg.cholesky(j2c, lower=True)
            chol = True
        except scipy.linalg.LinAlgError:
            chol = False
        r = dict(model=tag, functions='Cartesian (6D, 10F)' if cart else 'pure (5D, 7F)', nao=int(mol.nao),
                 naux=int(aux.nao), cholesky_succeeds=chol)
        if tag == 'small':
            w = np.linalg.eigvalsh(j2c)
            r.update(smallest_eigenvalue=float(w.min()), largest_eigenvalue=float(w.max()),
                     eigenvalues_below_1e_12=int((w < 1e-12).sum()), eigenvalues_below_1e_9=int((w < 1e-9).sum()))
        out.append(r)
        print(r, flush=True)
(HERE / 'metric_conditioning.json').write_text(json.dumps(dict(auxbasis='def2-universal-jkfit', rows=out,
                                                               versions=pc.versions()), indent=1) + '\n')
