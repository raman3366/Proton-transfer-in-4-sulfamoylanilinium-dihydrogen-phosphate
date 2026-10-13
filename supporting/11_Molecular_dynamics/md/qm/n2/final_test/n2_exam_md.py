#!/usr/bin/env python
"""Gate G-R2, for the record (no weight in the verdict): a 20 ps SMOKE TEST of a protein system built with a
candidate file of terms at N2 (NOT ACCEPTED terms; a throwaway build of n2_exam_build.py).  Not a result.

    OPENMM_CPU_THREADS=1 nice -n 10 /root/cqc-sfa-tools/md/bin/python md/qm/n2/final_test/n2_exam_md.py run \
        --candidate second_step --start X|R [--ps 20]
    ... n2_exam_md.py collect --candidate second_step      -> md/qm/n2/final_test/stability_<candidate>.json

Why not md/run_md.py: its minimisation alone took 29 minutes on this machine for the system as built (300 iterations
per pass, 2 threads; /root/cqc-sfa-md/systems/_checks/bonded_2e/R/min_300/run_md.log), and the shortened
equilibration as long again.  This driver is the examiner's own and is much shorter:
  - topology: the throwaway system (candidate terms and charges); positions and cell: the system AS BUILT after
    that minimisation (/root/cqc-sfa-md/systems/_checks/bonded_2e/<start>/min_300/state_min.xml; same atoms);
  - settings of md/protocol.json: PME, 1.0 nm cut-off, constraints on bonds to hydrogen, rigid water, hydrogen
    mass 3.024 amu, Langevin middle integrator at 300 K with 1 ps-1, Monte Carlo barostat at 1 bar every 25 steps;
  - 100 L-BFGS iterations with the new topology, NO restraint; velocities at 300 K from seed 20261006; 2 ps with a
    2 fs step; then --ps (20) with the 4 fs step of the production; nothing is held, the ligand is free from the
    first step;
  - every 0.2 ps: energies, temperature, the four Zn-N distances, the three angles at N2, the signed volume at
    N2, Zn-N2-S-C4, the two orientation distances of md/PLAN.md (aniline N to Phe131 CZ and to Asn67 OD1), N2 to
    Thr199 OG1, the nearest sulfonyl oxygen to the zinc.
  - --candidate as_built: the CONTROL, the system as built (/root/cqc-sfa-md/systems/bonded_2e/<start>, read only),
    same driver, same start coordinates, same seed.
The box as built has a density of 0.91 g/cm3 and no equilibration was made: this tests whether the system runs and
the zinc site holds with the new terms, nothing else.
"""
import argparse
import hashlib
import json
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path('/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate')
OUT = REPO / 'md' / 'qm' / 'n2' / 'final_test'
SCR = pathlib.Path('/tmp/claude-0/-root-workspace-cqc-4-sulfamoylanilinium-dihydrogen-phosphate/'
                   '82cc65d7-a11f-4d60-91ee-436dc9534941/scratchpad/n2/final')
SEED = 20261006


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def now():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def ang(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def measure(x, sel, lig):
    zn, n2, s, h, c4, n1 = x[sel['zn']], x[lig['N2']], x[lig['S']], x[sel['lig_HN2']], x[lig['C4']], x[lig['N1']]
    z, hh, uu = [(v - n2) / np.linalg.norm(v - n2) for v in (zn, h, s)]
    d = lambda a, b: float(np.linalg.norm(a - b))                        # noqa: E731
    return dict(Zn_N2=d(zn, n2), Zn_H94=d(zn, x[sel['his94_NE2']]), Zn_H96=d(zn, x[sel['his96_NE2']]), Zn_H119=d(zn, x[sel['his119_ND1']]),
                Zn_N2_S=ang(zn, n2, s), Zn_N2_H=ang(zn, n2, h), S_N2_H=ang(s, n2, h), signed_volume=float(np.dot(z, np.cross(hh, uu))),
                Zn_N2_S_C4=dih(zn, n2, s, c4), d1_N1_Phe131CZ=d(n1, x[sel['phe131_CZ']]), d2_N1_Asn67OD1=d(n1, x[sel['asn67_OD1']]),
                N2_Thr199OG1=d(n2, x[sel['thr199_OG1']]), Zn_nearest_O=min(d(zn, x[lig['O1']]), d(zn, x[lig['O2']])),
                C4_S_N2=ang(c4, s, n2))


def run(args):
    import openmm
    from openmm import app, unit as u
    # 'as_built': the CONTROL, the system as built (read only), with the same driver and the same start coordinates
    sysdir = (pathlib.Path('/root/cqc-sfa-md/systems/bonded_2e') / args.start if args.candidate == 'as_built' else
              SCR / f'build_{args.candidate}' / 'systems' / 'bonded_2e' / args.start)
    state0 = pathlib.Path(f'/root/cqc-sfa-md/systems/_checks/bonded_2e/{args.start}/min_300/state_min.xml')
    out = SCR / 'md' / args.candidate / args.start
    out.mkdir(parents=True, exist_ok=True)
    sel = json.loads((sysdir / 'selections.json').read_text())
    lig = dict(sel['lig_heavy'])
    prm = app.AmberPrmtopFile(str(sysdir / 'system.prmtop'))
    names = [a.name for a in prm.topology.atoms()]
    for k, i in lig.items():
        assert names[i] == k
    st0 = openmm.XmlSerializer.deserialize(state0.read_text())
    system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometer, constraints=app.HBonds, rigidWater=True,
                              ewaldErrorTolerance=0.0005, hydrogenMass=3.024 * u.amu, removeCMMotion=True)
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            f.setUseDispersionCorrection(True)
            f.setUseSwitchingFunction(False)
    system.addForce(openmm.MonteCarloBarostat(1.0 * u.bar, 300.0 * u.kelvin, 25))
    system.getForce(system.getNumForces() - 1).setRandomNumberSeed(SEED)
    integ = openmm.LangevinMiddleIntegrator(300.0 * u.kelvin, 1.0 / u.picosecond, 0.002 * u.picosecond)
    integ.setConstraintTolerance(1e-6)
    integ.setRandomNumberSeed(SEED)
    plat = openmm.Platform.getPlatformByName('CPU')
    sim = app.Simulation(prm.topology, system, integ, plat, {'Threads': str(args.threads)})
    sim.context.setPeriodicBoxVectors(*st0.getPeriodicBoxVectors())
    sim.context.setPositions(st0.getPositions())
    dof = sum(3 for i in range(system.getNumParticles()) if system.getParticleMass(i).value_in_unit(u.amu) > 0) - system.getNumConstraints() - 3
    kB = 0.0083144626

    def state():
        st = sim.context.getState(getEnergy=True, getPositions=True)
        x = np.array(st.getPositions(asNumpy=True).value_in_unit(u.angstrom))
        m = measure(x, sel, lig)
        m.update(potential_kJ=st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole),
                 temperature_K=2.0 * st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole) / (dof * kB),
                 volume_nm3=st.getPeriodicBoxVolume().value_in_unit(u.nanometer ** 3))
        return m, st
    log = open(out / 'progress.log', 'a')

    def say(t):
        log.write(f'{now()}  {t}\n')
        log.flush()
    rec = dict(what=('20 ps smoke test, CONTROL with the system as built' if args.candidate == 'as_built' else
                     '20 ps smoke test with candidate terms at N2 (NOT ACCEPTED)') + '; not a result', candidate=args.candidate, start=args.start,
               system=str(sysdir), prmtop_sha256=sha(sysdir / 'system.prmtop'), start_coordinates=str(state0), start_coordinates_sha256=sha(state0),
               seed=SEED, threads=args.threads, platform='CPU', openmm=openmm.__version__, started=now(), frames=[])
    m, _ = state()
    rec['as_started'] = m
    say(f'start: potential {m["potential_kJ"]:.1f} kJ/mol, signed volume {m["signed_volume"]:.3f}')
    t0 = time.time()
    sim.minimizeEnergy(maxIterations=100)
    m, _ = state()
    rec['after_100_iterations'] = m
    say(f'after 100 L-BFGS iterations ({time.time() - t0:.0f} s): potential {m["potential_kJ"]:.1f}')
    sim.context.setVelocitiesToTemperature(300.0 * u.kelvin, SEED)
    t = 0.0
    finite = True
    for phase, dt, n_rec in (('settle_2fs', 0.002, 10), ('run_4fs', 0.004, int(round(args.ps / 0.2)))):
        integ.setStepSize(dt * u.picosecond)
        for k in range(n_rec):
            t1 = time.time()
            sim.step(int(round(0.2 / dt)))
            t += 0.2
            m, st = state()
            m.update(phase=phase, time_ps=round(t, 3))
            rec['frames'].append(m)
            if not np.isfinite(m['potential_kJ']):
                finite = False
                break
            if k % 5 == 4:
                say(f'{phase} t {t:.1f} ps: T {m["temperature_K"]:.0f} K, Zn-N2 {m["Zn_N2"]:.3f}, Zn-N2-S {m["Zn_N2_S"]:.1f}, signed volume {m["signed_volume"]:.3f}, '
                    f'{(time.time() - t1) / 0.2:.0f} s per ps')
                rec.update(finite=finite, written=now(), wall_seconds=round(time.time() - t0))
                (out / 'record.json').write_text(json.dumps(rec, indent=1) + '\n')
        if not finite:
            break
    rec.update(finite=finite, finished=now(), wall_seconds=round(time.time() - t0), complete=bool(finite and t >= 2.0 + args.ps - 1e-6))
    (out / 'record.json').write_text(json.dumps(rec, indent=1) + '\n')
    (out / 'final_state.xml').write_text(openmm.XmlSerializer.serialize(sim.context.getState(getPositions=True, getVelocities=True)))
    say('done' if rec['complete'] else 'stopped: energy not finite')


def collect(args):
    doc = dict(what=('20 ps smoke tests, CONTROL: the systems as built' if args.candidate == 'as_built' else
                     '20 ps smoke tests of throwaway protein systems with candidate terms at N2 (NOT ACCEPTED)') +
                    '; for the record of gate G-R2, no weight in the verdict; not an equilibrated run and not a result',
               written=now(), candidate=args.candidate, this_program=sha(__file__), starts={})
    for start in ('X', 'R'):
        f = SCR / 'md' / args.candidate / start / 'record.json'
        if not f.exists():
            doc['starts'][start] = dict(done=False)
            continue
        r = json.loads(f.read_text())
        fr = [x for x in r['frames'] if x['phase'] == 'run_4fs']
        keys = ['Zn_N2', 'Zn_H94', 'Zn_H96', 'Zn_H119', 'Zn_N2_S', 'Zn_N2_H', 'S_N2_H', 'signed_volume', 'd1_N1_Phe131CZ', 'd2_N1_Asn67OD1',
                'N2_Thr199OG1', 'Zn_nearest_O', 'C4_S_N2', 'temperature_K']
        summ = {k: dict(mean=round(float(np.mean([x[k] for x in fr])), 3), min=round(float(np.min([x[k] for x in fr])), 3),
                        max=round(float(np.max([x[k] for x in fr])), 3)) for k in keys} if fr else {}
        doc['starts'][start] = dict(
            done=bool(r.get('complete')), finite=r.get('finite'), ps_at_4_fs=round(len(fr) * 0.2, 1), frames=len(fr), wall_seconds=r.get('wall_seconds'),
            started=r.get('started'), finished=r.get('finished'), prmtop_sha256=r['prmtop_sha256'], start_coordinates=r['start_coordinates'],
            as_started={k: round(v, 3) for k, v in r['as_started'].items() if isinstance(v, float)},
            summary_of_the_4_fs_part=summ,
            frames_with_signed_volume_below_0_3=int(sum(x['signed_volume'] < 0.3 for x in fr)),
            frames_with_signed_volume_above_minus_0_3=int(sum(x['signed_volume'] > -0.3 for x in fr)),
            sign_changes_of_the_signed_volume=int(sum(np.sign(a['signed_volume']) != np.sign(b['signed_volume']) for a, b in zip(fr[:-1], fr[1:]))),
            frames_d1_above_6_7_A=int(sum(x['d1_N1_Phe131CZ'] > 6.7 for x in fr)), frames_d2_below_5_2_A=int(sum(x['d2_N1_Asn67OD1'] < 5.2 for x in fr)),
            signed_volume_series=[round(x['signed_volume'], 2) for x in fr],
            sign_of_the_configuration_at_N2_changed=bool(fr and (np.sign(fr[0]['signed_volume']) != np.sign(fr[-1]['signed_volume'])
                                                               or len({np.sign(x['signed_volume']) for x in fr}) > 1)),
            last_frame={k: round(v, 3) for k, v in fr[-1].items() if isinstance(v, float)} if fr else None)
    (OUT / f'stability_{args.candidate}.json').write_text(json.dumps(doc, indent=1) + '\n')
    for s, v in doc['starts'].items():
        print(s, {k: v.get(k) for k in ('done', 'finite', 'ps_at_4_fs', 'wall_seconds')})
        for k, m in (v.get('summary_of_the_4_fs_part') or {}).items():
            print('   ', k, m)
        print('    signed volume:', v.get('signed_volume_series'))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--candidate', required=True)
    r.add_argument('--start', required=True, choices=['X', 'R'])
    r.add_argument('--ps', type=float, default=20.0)
    r.add_argument('--threads', type=int, default=1)
    c = sub.add_parser('collect')
    c.add_argument('--candidate', required=True)
    args = ap.parse_args()
    (run if args.cmd == 'run' else collect)(args)


if __name__ == '__main__':
    sys.exit(main())
