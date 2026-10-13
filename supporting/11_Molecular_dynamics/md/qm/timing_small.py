"""Stage S2, step 3: first timing of one SCF of the real small model, at the settings tried first.

    OMP_NUM_THREADS=4 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/timing_small.py [small|large]

Settings of this script (fixed here, NOT those of the later jobs): B3LYP/6-31G* with Cartesian d functions,
density fitting, grid (99,590).  Run on 5 Oct 2026 from 10:25 PM IST and stopped by its PID after two SCF
cycles (3.3 and 3.6 minutes each), because that already showed the cost; md/qm/timing_ultrafine_cycles.txt
holds the time stamps.  No timing_small.json exists for that reason.  The full timing of one SCF and one
gradient at the settings chosen afterwards is the job `timing` of md/qm/run_qm.py.
Log and checkpoint: /root/cqc-sfa-md/qm/timing/.
"""
import json
import os
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import pyscf_calc as pc  # noqa: E402

BULK = pathlib.Path('/root/cqc-sfa-md/qm/timing')


def main():
    model = sys.argv[1] if len(sys.argv) > 1 else 'small'
    com = gaucom.read_com(HERE / 'mcpb' / ('CA2_small_opt.com' if model == 'small' else 'CA2_large_mk.com'))
    BULK.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('PYSCF_TMPDIR', str(BULK))
    from pyscf import lib
    lib.param.TMPDIR = str(BULK)
    s = pc.method_of_com(com, pc.Settings(cart=True, grid=(99, 590)))
    rec = dict(model=model, com=os.path.relpath(com.path, HERE.parents[1]), settings=s.__dict__ | {'grid': list(s.grid)},
               b3lyp=pc.assert_b3lyp_is_gaussian(), versions=pc.versions(), started=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
               load_average_start=os.getloadavg())
    mol = pc.build_mol(com.symbols, com.coords, com.charge, com.multiplicity, s,
                       output=str(BULK / f'{model}.log'))
    mf = pc.build_mf(mol, s, chkfile=str(BULK / f'{model}_xtal.chk'))
    rec.update(natm=mol.natm, nao=int(mol.nao), nelectron=int(mol.nelectron))
    t0 = time.time()
    e, t_scf = pc.run_scf(mf)
    rec.update(e_tot=float(e), scf_seconds=round(t_scf, 1), scf_cycles=int(getattr(mf, 'cycles', -1)),
               naux=int(mf.with_df.get_naoaux()), grid_points=int(mf.grids.weights.size))
    (HERE / f'timing_{model}.json').write_text(json.dumps(rec, indent=1) + '\n')
    if model == 'small':
        t1 = time.time()
        g = mf.nuc_grad_method().kernel()
        rec.update(gradient_seconds=round(time.time() - t1, 1), gradient_rms=float(np.sqrt((g ** 2).mean())),
                   gradient_max=float(np.abs(g).max()))
        np.save(BULK / 'small_xtal_gradient.npy', g)
    rec.update(total_seconds=round(time.time() - t0, 1), load_average_end=os.getloadavg(),
               finished=time.strftime('%Y-%m-%d %H:%M:%S %Z'))
    import resource
    rec['max_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 0)
    (HERE / f'timing_{model}.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps(rec, indent=1))


if __name__ == '__main__':
    main()
