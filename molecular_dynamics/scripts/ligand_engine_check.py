"""Stage S1 of md/PLAN.md, check of the engine: the energies that the fit used (OpenMM through ParmEd, Reference
platform) against sander of AmberTools on the same topology and coordinates, in vacuum and in implicit water.

  /root/cqc-sfa-tools/md/bin/python md/ligand_engine_check.py          (a few seconds)

For each delivered topology (as assigned and repaired, both species) the potential energy is computed at two
geometries (the DFT minimum of the inpcrd file and the repaired minimum) with both programs.  Compared is the
energy DIFFERENCE between the two geometries, which is what fits and scans use, and the absolute energy.
sander: imin = 1, maxcyc = 1 (energy of the input geometry), no cut-off; igb = 0 for vacuum, igb = 5 (OBC2) with
extdiel = 78.3553, intdiel = 1, no salt, no surface term for implicit water.
Writes md/ligand/engine_check.json.
"""
import json, re, shutil
import numpy as np
import ligand_common as lc
from ligand_common import OUT, BULK, SPECIES

MDIN = {'vacuum': " single point\n &cntrl\n  imin=1, maxcyc=1, ncyc=1, ntb=0, igb=0, cut=999.0, ntpr=1,\n /\n",
        'OBC2': " single point\n &cntrl\n  imin=1, maxcyc=1, ncyc=1, ntb=0, igb=5, extdiel=78.3553, intdiel=1.0,\n"
                "  saltcon=0.0, gbsa=0, cut=999.0, rgbmax=999.0, ntpr=1,\n /\n"}


def write_rst7(path, x):
    v = np.asarray(x).ravel()
    with open(path, 'w') as f:
        f.write('single point\n%5d\n' % (len(v) // 3))
        for i in range(0, len(v), 6):
            f.write(''.join('%12.7f' % a for a in v[i:i + 6]) + '\n')


def sander_energy(work, prmtop, x, envn):
    (work / 'mdin').write_text(MDIN[envn])
    write_rst7(work / 'in.rst7', x)
    lc.run(['sander', '-O', '-i', 'mdin', '-o', 'mdout', '-p', prmtop, '-c', 'in.rst7', '-r', 'out.rst7'], work)
    txt = (work / 'mdout').read_text()
    m = re.search(r'NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER\s*\n\s*1\s+(\S+)', txt)
    return float(m.group(1))


def main():
    work = BULK / 'engine_check'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    out = dict(script='md/ligand_engine_check.py', versions=lc.versions(), rows=[])
    worst = 0.0
    for sp in SPECIES:
        el, x_min, _ = lc.read_xyz(OUT / f'sfa_{sp}_repaired_min.xyz')
        el, x_dft, _ = lc.read_xyz(lc.RESULTS / SPECIES[sp]['case'] / 'opt.xyz')
        for tag in ('gaff2', 'repaired'):
            prm = OUT / f'sfa_{sp}_{tag}.prmtop'
            shutil.copy(prm, work / 'x.prmtop')
            for envn in ('vacuum', 'OBC2'):
                mm = lc.MM(prm, implicit=None if envn == 'vacuum' else envn)
                eo = [mm.energy(x_dft), mm.energy(x_min)]
                es = [sander_energy(work, 'x.prmtop', x_dft, envn), sander_energy(work, 'x.prmtop', x_min, envn)]
                row = dict(species=sp, topology=f'sfa_{sp}_{tag}.prmtop', environment=envn,
                           openmm_E_dft_geom=eo[0], sander_E_dft_geom=es[0],
                           openmm_dE=eo[1] - eo[0], sander_dE=es[1] - es[0],
                           diff_absolute=eo[0] - es[0], diff_dE=(eo[1] - eo[0]) - (es[1] - es[0]))
                worst = max(worst, abs(row['diff_dE']))
                out['rows'].append(row)
                print(f"{sp:8s} {tag:9s} {envn:7s} E(DFT geometry) OpenMM {eo[0]:11.4f} sander {es[0]:11.4f}; "
                      f"dE OpenMM {eo[1] - eo[0]:9.4f} sander {es[1] - es[0]:9.4f}; difference of dE {row['diff_dE']:+.4f}")
    out['largest_difference_of_dE_kcal'] = worst
    out['note'] = ('sander prints energies to four decimals. Implicit water: sander igb=5 and the OBC2 force that ParmEd '
                   'builds for OpenMM use the same radii (those of the topology) and the same formula; differences of '
                   'absolute energies come from the constants of the two programs.')
    (OUT / 'engine_check.json').write_text(json.dumps(out, indent=1) + '\n')
    print('largest difference of relative energies:', round(worst, 4), 'kcal/mol')


if __name__ == '__main__':
    main()
