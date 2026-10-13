"""Stage S2, step 3: where the time of one SCF cycle of the small model goes (64 atoms, 585 functions).

    OMP_NUM_THREADS=4 nice -n 10 /root/cqc-sfa-tools/env/bin/python md/qm/bench_components.py

Times, with the guess density, (1) the build of the density-fitting integrals, (2) one Coulomb-and-exchange
build, (3) one exchange-correlation build on several grids.  Result: md/qm/bench_components.json.
Run on 5 Oct 2026, 10:45 to 11:05 PM IST, with Cartesian d functions (585 functions, 3,837 fitting functions
after PySCF dropped the linearly dependent ones); the machine was shared (load average 8 to 9.5).
"""
import json
import os
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gaucom  # noqa: E402
import pyscf_calc as pc  # noqa: E402

BULK = pathlib.Path('/root/cqc-sfa-md/qm/timing')


def main():
    from pyscf import dft, lib
    BULK.mkdir(parents=True, exist_ok=True)
    lib.param.TMPDIR = str(BULK)
    com = gaucom.read_com(HERE / 'mcpb' / 'CA2_small_opt.com')
    s = pc.Settings(verbose=0, cart=True)      # as run on 5 Oct 2026 (585 functions); the jobs use pure functions
    mol = pc.build_mol(com.symbols, com.coords, com.charge, com.multiplicity, s)
    rec = dict(natm=mol.natm, nao=int(mol.nao), versions=pc.versions(), load_average_start=os.getloadavg(),
               started=time.strftime('%Y-%m-%d %H:%M:%S %Z'))
    mf = pc.build_mf(mol, s)
    dm = mf.get_init_guess()
    t = time.time()
    mf.with_df.build()
    rec['df_build_seconds'] = round(time.time() - t, 1)
    rec['naux'] = int(mf.with_df.get_naoaux())
    t = time.time()
    vj, vk = mf.with_df.get_jk(dm, hermi=1)
    rec['jk_seconds'] = round(time.time() - t, 1)
    t = time.time()
    vj = mf.with_df.get_jk(dm, hermi=1, with_k=False)[0]
    rec['j_only_seconds'] = round(time.time() - t, 1)
    rec['grids'] = []
    ni = dft.numint.NumInt()
    for label, kw in (('level 1', dict(level=1)), ('level 2', dict(level=2)), ('level 3 (PySCF default)', dict(level=3)),
                      ('(75,302)', dict(atom_grid=(75, 302))), ('(99,590)', dict(atom_grid=(99, 590)))):
        g = dft.gen_grid.Grids(mol)
        for k, v in kw.items():
            setattr(g, k, v)
        t = time.time()
        g.build()
        tb = time.time() - t
        t = time.time()
        n, exc, vxc = ni.nr_rks(mol, g, 'b3lyp', dm)
        rec['grids'].append(dict(grid=label, points=int(g.weights.size), build_seconds=round(tb, 1),
                                 xc_seconds=round(time.time() - t, 1), electrons=float(n), exc=float(exc)))
        (HERE / 'bench_components.json').write_text(json.dumps(rec, indent=1) + '\n')
    rec['load_average_end'] = os.getloadavg()
    (HERE / 'bench_components.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps(rec, indent=1))


if __name__ == '__main__':
    main()
