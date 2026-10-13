"""Gate G-R2: does tleap put a non-zero torsion through the zinc, and an improper at N2, on the intended atoms
only, and does the list of 1-4 pairs stay the same?  (md/N2.md, section 2.3.)  A throwaway build, no model.

    AMBERHOME=/root/cqc-sfa-tools/md /root/cqc-sfa-tools/md/bin/python md/qm/n2_leaptest.py WORKDIR

Builds the small model twice with the inputs of md/qm/n2/ff/inputs (those of gate G-R): as built, and with a
test frcmod loaded last that holds one term of every class the plan allows (arbitrary numbers).  Compares the
two topologies term by term (n2_gate.topology_terms, compare_topologies) and the exception lists of OpenMM.
Writes md/qm/n2/leap_test.json.
"""
import json
import pathlib
import sys

import n2_common as nc
import n2_gate as ng

TEST = {'K_ZnN2S': 40.0, 'B_ZnN2S': -2 * 40.0 * 1.95, 'K_ZnN2H': 30.0, 'B_ZnN2H': -2 * 30.0 * 2.05,
        'K_SN2H': 55.0, 'B_SN2H': -2 * 55.0 * 1.98,
        'H_NS_C4_1': -0.9, 'H_NS_C4_2': 2.1, 'H_NS_C4_3': 0.4, 'Zn_NS_C4_1': -0.5, 'Zn_NS_C4_2': 0.3, 'Zn_NS_C4_3': 0.2,
        'His_ZnN_S_3c': 0.5, 'His_ZnN_S_3s': -0.5, 'His_ZnN_H_3c': 0.0, 'His_ZnN_H_3s': 0.4, 'improper': 1.1}


def main():
    work = pathlib.Path(sys.argv[1])
    frc = work / 'test.frcmod'
    work.mkdir(parents=True, exist_ok=True)
    terms = nc.write_frcmod(TEST, frc, 'throwaway test of md/qm/n2_leaptest.py: one term of every allowed class')
    a = nc.build_topology(None, work / 'as_built')
    b = nc.build_topology(frc, work / 'test')
    ta, tb = ng.topology_terms(a), ng.topology_terms(b)
    allowed, other = ng.compare_topologies(ta, tb)
    from openmm import app

    def exceptions(f):
        s = app.AmberPrmtopFile(str(f)).createSystem(nonbondedMethod=app.NoCutoff)
        nb = [x for x in s.getForces() if type(x).__name__ == 'NonbondedForce'][0]
        return sorted((nb.getExceptionParameters(k)[0], nb.getExceptionParameters(k)[1],
                       round(nb.getExceptionParameters(k)[2]._value, 8), round(nb.getExceptionParameters(k)[4]._value, 8))
                      for k in range(nb.getNumExceptions()))
    ea, eb = exceptions(a), exceptions(b)
    zn14 = sorted(p[0][0] if p[0][1] == 'ZN1:ZN' else p[0][1] for p in ta['pairs14'] if 'ZN1:ZN' in p[0])
    changed = {k: v for k, v in allowed['torsions'].items() if v['changed']}
    doc = dict(
        what='throwaway tleap build: one term of every class gate G-R2 may change; written by md/qm/n2_leaptest.py',
        written=nc.now(), test_parameters=TEST, frcmod_text=frc.read_text().splitlines(), amber_numbers_written=terms,
        arrived=dict(angles=allowed['angles'], torsions_changed=changed, improper=allowed['improper'],
                     torsions_of_the_allowed_classes_unchanged=sorted(k for k, v in allowed['torsions'].items() if not v['changed'])),
        differences_outside_the_allowed_classes=other,
        pairs14=dict(n_as_built=len(ta['pairs14']), n_test=len(tb['pairs14']), identical_with_scaling=ta['pairs14'] == tb['pairs14'],
                     scaling=sorted({(p[1], p[2]) for p in ta['pairs14']}), partners_of_the_zinc=zn14),
        openmm_exceptions=dict(n=len(ea), identical=ea == eb),
        rebuilt_as_built_equals_kept_topology=ng.topology_terms(nc.FFDIR / 'small_as_built.prmtop') == ta,
        ok=bool(not other and ta['pairs14'] == tb['pairs14'] and ea == eb and len(allowed['improper']) == 1))
    nc.write_json(nc.N2 / 'leap_test.json', doc)
    print(json.dumps({k: doc[k] for k in ('arrived', 'differences_outside_the_allowed_classes', 'pairs14', 'openmm_exceptions',
                                           'rebuilt_as_built_equals_kept_topology', 'ok')}, indent=1))


if __name__ == '__main__':
    main()
