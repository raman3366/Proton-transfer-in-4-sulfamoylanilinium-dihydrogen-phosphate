"""Gate G-R2, declared second step: what a change of the ligand charges does to the three numbers of gate G-L.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_gl_check.py --work SCRATCH
        [--out md/qm/n2/gl_check.json]

Gate G-L (md/LIGAND.md, section 5; md/ligand_fit.py) calibrated the ligand terms at sulfur on the FREE anion with
AM1-BCC charges, in implicit water: at the force-field minimum C4-S-N2 within 2 deg and S-N2 within 0.02 A of the
DFT minimum, and the relaxed scan about C4-S within 0.5 kcal/mol.  The bound ligand of the simulation never
carries those charges: it carries RESP charges of the zinc-site model.  This script repeats the three numbers of
G-L with the tools of gate G-L (md/ligand_common.py, md/ligand_eval.py, imported and not changed; the procedure of
ligand_fit.charge_sensitivity) for three charge sets on the same free anion, every other term as delivered
(md/ligand/sfa_anion_gaff2.frcmod, sfa_anion_repair.frcmod):

    am1bcc          md/ligand/sfa_anion.mol2 as it is: must return the numbers of record (a check of this script)
    resp_as_built   the ligand charges of the systems as built      (md/qm/n2/ff/inputs/SN1.mol2)
    resp_two        the ligand charges fitted to two placements     (md/qm/n2/ff2/inputs/SN1.mol2)

WHAT IT CAN AND CANNOT SAY.  The RESP charges of the bound ligand sum to about -0.5 e, not -1 (the rest sits on
the zinc and the histidines), so a free molecule with these charges is not the anion in water: the numbers are a
SENSITIVITY of the ligand's own geometry and C4-S profile to the charge set, like for like between the two RESP
sets, and not a second reading of gate G-L.  Nothing is written outside --work and --out.  No energy of a point
of gate G-R2 is read.
"""
import argparse
import json
import pathlib
import shutil
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))                                       # md/: the ligand tools of stage S1
sys.path.insert(0, str(HERE))
import ligand_common as lc                                                 # noqa: E402
import ligand_eval as le                                                   # noqa: E402
import n2_common as nc                                                     # noqa: E402

SP = 'anion'
SETS = (('am1bcc', None), ('resp_as_built', nc.N2 / 'ff' / 'inputs' / 'SN1.mol2'),
        ('resp_two', nc.N2 / 'ff2' / 'inputs' / 'SN1.mol2'))


def charges_of(mol2):
    lines = pathlib.Path(mol2).read_text().split('\n')
    a0, a1 = lines.index('@<TRIPOS>ATOM') + 1, lines.index('@<TRIPOS>BOND')
    return {l.split()[1]: float(l.split()[8]) for l in lines[a0:a1]}


def with_charges(src, q, dst):
    lines = pathlib.Path(src).read_text().split('\n')
    a0, a1 = lines.index('@<TRIPOS>ATOM') + 1, lines.index('@<TRIPOS>BOND')
    for k in range(a0, a1):
        f = lines[k].split()
        lines[k] = lines[k][:lines[k].rindex(f[8])] + f'{q[f[1]]:9.6f}'.strip().rjust(len(f[8]))
    pathlib.Path(dst).write_text('\n'.join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--out', default=str(nc.N2 / 'gl_check.json'))
    args = ap.parse_args()
    work = pathlib.Path(args.work)
    el, x_dft, _ = lc.read_xyz(lc.RESULTS / lc.SPECIES[SP]['case'] / 'opt.xyz')
    names = lc.names_for(SP, el, x_dft)
    m_dft = lc.add_means(lc.measures(x_dft, names))
    base = lc.OUT / f'sfa_{SP}.mol2'
    q0 = charges_of(base)
    doc = dict(what='the three numbers of gate G-L on the free anion with three charge sets; a sensitivity, not a '
                    'second reading of the gate; written by md/qm/n2_gl_check.py', written=nc.now(),
               environment='implicit water, OBC2, dielectric 78.3553 (as md/ligand_fit.py)',
               dft=dict(a_C4_S_N2=m_dft['a_C4_S_N2'], d_S_N2=m_dft['d_S_N2']),
               tolerances_of_G_L=dict(a_C4_S_N2_deg=2.0, d_S_N2_A=0.02, scan_kcal=0.5),
               inputs={str(p.relative_to(nc.REPO)): nc.sha(p) for p in (
                   base, lc.OUT / f'sfa_{SP}_gaff2.frcmod', lc.OUT / f'sfa_{SP}_repair.frcmod',
                   nc.N2 / 'ff' / 'inputs' / 'SN1.mol2', nc.N2 / 'ff2' / 'inputs' / 'SN1.mol2')},
               sets={})
    for tag, src in SETS:
        w = work / tag
        if w.exists():
            shutil.rmtree(w)
        w.mkdir(parents=True)
        q = dict(q0) if src is None else charges_of(src)
        assert set(q) == set(q0), (tag, sorted(set(q) ^ set(q0)))
        with_charges(base, q, w / 'lig.mol2')
        for f in (f'sfa_{SP}_gaff2.frcmod', f'sfa_{SP}_repair.frcmod'):
            shutil.copy(lc.OUT / f, w / f)
        (w / 'leap.in').write_text(
            f'source leaprc.gaff2\nloadamberparams sfa_{SP}_gaff2.frcmod\nloadamberparams sfa_{SP}_repair.frcmod\n'
            'SAN = loadmol2 lig.mol2\nsaveamberparm SAN lig.prmtop lig.inpcrd\nquit\n')
        lc.run(['tleap', '-f', 'leap.in'], w, log=w / 'leap.out')
        assert (w / 'lig.prmtop').exists(), (w / 'leap.out').read_text()[-1500:]
        mm = lc.MM(w / 'lig.prmtop', implicit='OBC2')
        r = mm.minimise(x_dft)
        m = lc.add_means(lc.measures(r['x'], names))
        ss = le.scan_summary(le.scan(mm, SP, r['x'], names, global_search=False))
        g = le.gate(m, m_dft, ss)
        doc['sets'][tag] = dict(
            charges_from='md/ligand/sfa_anion.mol2 (AM1-BCC, the charges of the calibration)' if src is None
            else str(src.relative_to(nc.REPO)), total_charge=round(sum(q.values()), 6),
            charges={k: q[k] for k in ('N2', 'S', 'O1', 'O2', 'C4', 'HN21')},
            a_C4_S_N2=round(m['a_C4_S_N2'], 3), a_C4_S_N2_minus_DFT=round(g['d_csn'], 3),
            d_S_N2=round(m['d_S_N2'], 4), d_S_N2_minus_DFT=round(g['d_sn'], 4),
            scan_largest_deviation=round(ss['max_abs_dev'], 4), scan_mm_barrier=round(ss['mm_barrier'], 4),
            scan_dft_barrier=round(ss['dft_barrier'], 4),
            within_the_tolerances_of_G_L=dict(a_C4_S_N2=g['csn_pass'], d_S_N2=g['sn_pass'], scan=g['scan_pass']))
        print(tag, json.dumps({k: v for k, v in doc['sets'][tag].items() if k not in ('charges', 'charges_from')}), flush=True)
    a, b = doc['sets']['resp_as_built'], doc['sets']['resp_two']
    doc['two_placements_minus_as_built'] = dict(
        a_C4_S_N2=round(b['a_C4_S_N2'] - a['a_C4_S_N2'], 3), d_S_N2=round(b['d_S_N2'] - a['d_S_N2'], 4),
        scan_largest_deviation=round(b['scan_largest_deviation'] - a['scan_largest_deviation'], 4))
    nc.write_json(pathlib.Path(args.out), doc)
    print('written', args.out)


if __name__ == '__main__':
    main()
