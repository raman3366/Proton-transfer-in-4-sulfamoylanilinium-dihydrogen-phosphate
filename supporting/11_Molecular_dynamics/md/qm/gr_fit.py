"""Gate G-R, steps 3 and 4: quantum chemistry against the force field on the small zinc model, read in two
parts (md/PLAN.md, entry of 6 Oct 2026, 8:56 AM IST), and the fit of the terms at N2 if part A fails.

    /root/cqc-sfa-tools/md/bin/python md/qm/gr_fit.py gate --mm MM.json --qm QMDIR --out OUTDIR [--name gr]
    /root/cqc-sfa-tools/md/bin/python md/qm/gr_fit.py fit  --mm MM.json --qm QMDIR --out OUTDIR
    /root/cqc-sfa-tools/md/bin/python md/qm/gr_fit.py selftest --mm MM.json --out OUTDIR

  MM.json  result of md/qm/gr_mm.py run (energies and their split for all points and three cap schemes)
  QMDIR    one <id>.json per finished single point (md/qm/gr_queue.py): e_tot in Hartree, converged
  Points and sets come from the points.json named in MM.json; the D3(BJ) energy of every geometry from
  d3.json beside it (md/qm/gr_queue.py d3); the relaxed references from relax_mm.json and relax_qm.json beside
  it and the binding check from binding_check.json, when these files exist (--side DIR names another folder
  for these four files; --ff-side DIR the folder of the two that depend on the force-field model,
  relax_mm.json and binding_check.json, when a reading of the gate has its own folder: md/qm/gr_queue.py,
  READINGS; --reading and --step2 are written into the result).

ENERGIES.  E_QM = B3LYP/6-31G* energy + B3LYP-D3(BJ) energy of the same geometry (plan, point 1: the force field
has dispersion in its Lennard-Jones terms, the functional has none).  Every number is also written without the
D3 term (`without_d3`); the verdict uses the numbers with it.

PART A, the terms at N2.  Points: X, the inner points of the inversion path (a2 to a6), the six points in which
only the hydrogen on N2 moves (d), the turns of the ligand about Zn-N2 (e), the bend of Zn-N2-S (f), and Rx with
the six hydrogen points made on it (g).  The points form two FAMILIES by the geometry they were made from: X
(with a, d, e, f) and Rx (with g).  In each family the reference r is its point of lowest E_QM, and
      d_i = [E_MM(i) - E_MM(r)] - [E_QM(i) - E_QM(r)]
      rms = sqrt( mean of d_i^2 over all points of part A that are not a reference )     criterion: <= 1.0
  Why two references: the points (g) differ from X by the whole difference X to Rx, which is contact energy and
  belongs to part B; against their own base Rx they test the terms at N2, as (d) does against X.  The number
  with ONE reference for all points (the wording of the plan before the families existed) is written beside
  as `rms_single_reference` and does not decide.  Also beside: without the points flagged protein_contact,
  and for the other two cap schemes.

PART B, contacts.  X against Rx, the turns about N2-S (c1 from X, c2 from Rx), the N-like point Nx, and the two
relaxed references.  Reported with the split of the force-field energy by part (md/qm/gr_mm.py), never
fitted, no verdict.  For the relaxed references: E(relaxed from Rx) - E(relaxed from X) in the force field and
in quantum chemistry; a mismatch above 3.0 kcal/mol is flagged as a limit of the model.

REPORTED ONLY: R and N as built (docked geometries of the neutral molecule).

FIT, only if part A is complete and has failed.  Everything but the terms at N2 is fixed:
      y_i = E_QM(i) - rest(i),   rest = force-field energy without the terms fitted
      y_i ~ c_family + sum_a k_a (theta_a(i) - theta0_a)^2 + c1 cos(phi_i) + c2 cos(2 phi_i)
  a = Zn-N2-S, Zn-N2-H, S-N2-H; phi = H-N2-S-C4; c_n > 0 is an Amber term with phase 0, c_n < 0 with phase
  180 deg and constant |c_n|; one free constant per family.  Both models are linear least-squares problems:
      M1  force constants kept, three equilibrium values and two torsion coefficients: 5 parameters
      M2  force constants free as well: 8 parameters
  Points: those of part A.  Leave-one-out: each point is predicted from a fit without it.
  A model is ACCEPTABLE when all of this holds (rule of md/GR.md of 5 AM with the additions of the plan entry of
  8:56 AM, fixed before any energy but those of X, R and Rx was read; the last two lines as amended by the plan
  entry of 6 Oct 2026, 1:09 PM IST, point 4, for every reading of the gate from then on):
      every k_a in 10-300 kcal/mol/rad^2; every standard error of a theta0 at most 5 deg and (M2) of a k_a at
      most 30 % of k_a;
      residual r.m.s. <= 1.0 and leave-one-out r.m.s. <= 1.0 kcal/mol;
      |c1| < 5 and |c2| < 5 kcal/mol;
      every theta0 within 10 deg of its reference: Zn-N2-S of the crystal value (start X, 117.55 deg); Zn-N2-H
      and S-N2-H of their values in the QUANTUM-CHEMICAL OPTIMUM of the small model with all atoms free (the
      geometry task qm_opt delivers, read from binding_check.json), not of the constructed hydrogen of start X
      (107.41 deg).  Until that optimum exists the condition cannot be met: it says so and refuses;
      the number of part A with the fitted terms lower than before AND at most 1.0.
  M2 is taken when it is acceptable and its leave-one-out r.m.s. is at least 20 % below that of M1; otherwise
  M1 when it is acceptable; otherwise NO parameter file is written, the terms stay as built and the gate is
  reported as failed, with the split by energy part (`fit.json`: failed, part_A_split).  A model that meets
  every condition but the last half of the last one (part A falls, but stays above 1.0) is reported as
  "terms improved, gate still failed" (`fit.json`: outcome) and is not used either.
  Output: fit.json, fit_residuals.csv and, only for an accepted model, CA2_n2_fit.frcmod.  md/qm/gr_queue.py
  then rebuilds the topology with the frcmod and runs `gate` again; that number is the one that counts.

SELF-TEST.  Synthetic energies on the real geometries, seed 20261006: (T0) the force field itself, with and
without noise: the gate must pass and nothing may be fitted; (T0b) the reviewer's test, force field minus
D3(BJ) read without the D3 term: the right answer is "change nothing"; (T0c) the force field plus a contact
error on everything made from Rx; (T1, T2) known other parameters inside the rule: they must come back and be
accepted; (T3) known parameters outside the rule: no parameter file; (T4) the rule without the
quantum-chemical optimum: refused, with the reason; (T5) a fit that lowers part A but leaves it above 1.0:
"terms improved, gate still failed", no parameter file; (T6) Zn-N2-H 11 deg from the optimum: refused; (T7)
Zn-N2-H and S-N2-H at the optimum itself (15 deg from the constructed hydrogen): accepted.
"""
import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np

HARTREE = 627.5094740631
ANGLES = ('Zn_N2_S', 'Zn_N2_H', 'S_N2_H')
FRC_NAMES = dict(Zn_N2_S='M1-Y4-sy', Zn_N2_H='M1-Y4-hn', S_N2_H='hn-Y4-sy')
K_RANGE = (10.0, 300.0)
SE_THETA_MAX, SE_K_REL_MAX = 5.0, 0.30
GATE = 1.0                      # kcal/mol: part A, residual and leave-one-out of a fit
C_MAX = 5.0                     # kcal/mol: torsion coefficients of a fit
THETA_WINDOW = 10.0             # deg: equilibrium angles of a fit around the structure
RELAX_LIMIT = 3.0               # kcal/mol: mismatch of the relaxed references stated as a limit of the model
SEED = 20261006
SPLIT = ('bond', 'angle_N2', 'angle_other', 'torsion_N2S', 'torsion_other', 'nb_ligand', 'lj_ligand_frame',
         'coulomb_ligand_frame')
PLAN = 'md/PLAN.md, entries of 6 Oct 2026, 8:56 AM IST and 1:09 PM IST'
IMPROVED_ONLY = 'terms improved, gate still failed'


def now_ist():
    return time.strftime('%Y-%m-%d %I:%M:%S %p IST', time.gmtime(time.time() + 19800))


def read_json(path, default=None):
    try:
        return json.loads(pathlib.Path(path).read_text())
    except Exception:                                                      # noqa: BLE001
        return default


def load(mm_path, qm_dir, side=None, ff_side=None):
    mm = json.loads(pathlib.Path(mm_path).read_text())
    pfile = pathlib.Path(mm['points_file'])
    side = pathlib.Path(side) if side else pfile.parent
    ff_side = pathlib.Path(ff_side) if ff_side else side
    pts = json.loads(pfile.read_text())
    info = {p['id']: p for p in pts['points']}
    for p in info.values():                                                # a points file of the first version
        s = p['sets']
        s.setdefault('base', 'Rx' if (p['id'].startswith('c2_') or p['id'] == 'Rx') else 'X')
        s.setdefault('part_A', p['id'] in ('X', 'Rx') or (p['group'] in ('a', 'd')))
        s.setdefault('part_B', p['id'] in ('X', 'Rx', 'Nx') or p['group'] == 'c')
        s.setdefault('report_only', p['id'] in ('R', 'N'))
    qm = {}
    if qm_dir:
        for pid in pts['order']:
            r = read_json(pathlib.Path(qm_dir) / f'{pid}.json')
            if r and r.get('converged') and r.get('e_tot') is not None and r.get('xyz_sha256') == info[pid]['xyz_sha256']:
                qm[pid] = r
    d3doc = read_json(side / 'd3.json') or {}
    d3 = {}
    for pid, rec in (d3doc.get('points') or {}).items():
        if pid in info and rec.get('xyz_sha256') == info[pid]['xyz_sha256']:
            d3[pid] = rec['e_d3bj']
    aux = dict(dir=side, ff_dir=ff_side, d3_file=str(side / 'd3.json'), d3_checks=(d3doc.get('checks') or {}).get('passed'),
               d3_routine=d3doc.get('routine'))
    return mm, pts, info, qm, d3, aux


FIRST_READING = '2e'                                                       # the reading of 6 Oct 2026 (md/GR.md)


def reading_of(args, mm):
    """Name of the reading and its step: as given, else from the label of the force-field model ("step 2e, ...")."""
    lab = str(mm.get('label', ''))
    step = args.step2 or (lab.split(',')[0].split()[1] if lab.startswith('step ') else FIRST_READING)
    return (args.reading or step), step


def own_file(path, reading):
    """A reading never overwrites the record of another one (md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 4): a file
    that exists and names another reading (no name: the reading of 6 Oct 2026, "2e") stops the program."""
    old = read_json(path)
    if old is not None and old.get('reading', FIRST_READING) != reading:
        raise SystemExit(f"{path} is the record of the reading {old.get('reading', FIRST_READING)!r}; the reading {reading!r} "
                         f"must write under its own names (md/qm/gr_queue.py, READINGS): nothing written")


# ------------------------------------------------------------------------------------------------ the two parts
def families(info):
    fam = {}
    for p, i in info.items():
        if i['sets']['part_A']:
            fam.setdefault(i['sets']['base'], []).append(p)
    return fam


def part_a(e_qm, e_mm, info, drop_flag=None):
    """e_qm, e_mm: kcal/mol for the finished points (any common zero)."""
    fam = families(info)
    need = [p for ids in fam.values() for p in ids if not (drop_flag and drop_flag in info[p]['flags'])]
    done = [p for p in need if p in e_qm and p in e_mm]
    out = dict(n_points=len(done), n_expected=len(need), complete=bool(len(done) == len(need)), families={}, deviations={})
    dev = {}
    for base, ids in fam.items():
        ids = [p for p in ids if p in done]
        if len(ids) < 2:
            out['families'][base] = dict(points=ids, reference=None)
            continue
        ref = min(ids, key=lambda p: e_qm[p])
        d = {p: (e_mm[p] - e_mm[ref]) - (e_qm[p] - e_qm[ref]) for p in ids}
        v = np.array([d[p] for p in ids if p != ref])
        allv = np.array([d[p] for p in ids])
        out['families'][base] = dict(points=ids, reference=ref, n_compared=len(ids) - 1,
                                     rms=round(float(np.sqrt((v ** 2).mean())), 4),
                                     rms_after_best_shift=round(float(np.sqrt(((allv - allv.mean()) ** 2).mean())), 4))
        dev.update({p: d[p] for p in ids if p != ref})
    out['deviations'] = {p: round(v, 4) for p, v in dev.items()}
    if dev:
        v = np.array(list(dev.values()))
        worst = max(dev, key=lambda p: abs(dev[p]))
        out.update(n_compared=len(dev), rms=round(float(np.sqrt((v ** 2).mean())), 4), largest=round(dev[worst], 4),
                   largest_at=worst, passed=bool(np.sqrt((v ** 2).mean()) <= GATE))
        ref = min(done, key=lambda p: e_qm[p])
        one = np.array([(e_mm[p] - e_mm[ref]) - (e_qm[p] - e_qm[ref]) for p in done if p != ref])
        if len(one):
            out.update(rms_single_reference=round(float(np.sqrt((one ** 2).mean())), 4), single_reference=ref)
    return out


def subsets(e_qm, e_mm, info):
    """Part A over parts of its set (reported only; check of 6 Oct 2026, finding 7).  Each subset is read like part
    A itself: two families, each against its own lowest quantum-chemical point WITHIN the subset."""
    fam = families(info)
    done = {b: [p for p in ids if p in e_qm and p in e_mm] for b, ids in fam.items()}

    def rms_of(keep, own_base=False):
        dv = []
        for b, ids in done.items():
            ids = [p for p in ids if keep(p)]
            if len(ids) < 2:
                continue
            r = b if (own_base and b in ids) else min(ids, key=lambda p: e_qm[p])
            dv += [(e_mm[p] - e_mm[r]) - (e_qm[p] - e_qm[r]) for p in ids if p != r]
        return dict(n_compared=len(dv), rms=(round(float(np.sqrt(np.mean(np.square(dv)))), 4) if dv else None))
    g = lambda p: info[p]['group']
    return dict(
        how='two families, each against its lowest quantum-chemical point within the subset; reported only, no verdict',
        without_the_inversion_path=rms_of(lambda p: g(p) != 'a'),
        without_the_turns_about_Zn_N2=rms_of(lambda p: g(p) != 'e'),
        hydrogen_points_and_bends_only=rms_of(lambda p: g(p) not in ('a', 'e')),
        X_and_Rx_themselves_as_references=rms_of(lambda p: True, own_base=True))


def composition(e_qm, mm, info):
    """What moves the force-field energy of each point of part A against the geometry it was made from (X or Rx):
    the fitted terms (three angles at N2 and H-N2-S-C4), the non-bonded energy between ligand and fixed atoms, the
    non-bonded energy inside the ligand, the rest (check of 6 Oct 2026, finding 7).  kcal/mol."""
    out = {}
    for p, i in info.items():
        b = i['sets']['base']
        if not i['sets']['part_A'] or p == b or p not in mm['points'] or b not in mm['points']:
            continue
        sp, sb = mm['points'][p].get('split'), mm['points'][b].get('split')
        if not sp or not sb:
            continue
        d = {k: sp[k] - sb[k] for k in SPLIT}
        tot = mm['points'][p]['energy']['total'] - mm['points'][b]['energy']['total']
        fit = mm['points'][p]['energy']['fit_terms'] - mm['points'][b]['energy']['fit_terms']
        out[p] = dict(base=b, group=i['group'], mm=round(tot, 3),
                      qm=(round(e_qm[p] - e_qm[b], 3) if (p in e_qm and b in e_qm) else None),
                      fitted_terms=round(fit, 3), angles_at_N2=round(d['angle_N2'], 3),
                      ligand_frame_nonbonded=round(d['lj_ligand_frame'] + d['coulomb_ligand_frame'], 3),
                      ligand_frame_electrostatic=round(d['coulomb_ligand_frame'], 3),
                      nonbonded_inside_ligand=round(d['nb_ligand'], 3),
                      rest=round(tot - fit - d['lj_ligand_frame'] - d['coulomb_ligand_frame'] - d['nb_ligand'], 3))
    return out


def rel_split(mm, pid, base):
    a, b = mm['points'][base].get('split'), mm['points'][pid].get('split')
    if not a or not b:
        return None
    return {k: round(b[k] - a[k], 3) for k in SPLIT}


def compare(e_qm, e_mm, mm, pid, base):
    """One point against another: quantum chemistry, force field, their difference, and the split."""
    if pid not in e_qm or base not in e_qm or pid not in e_mm:
        return None
    q, m = e_qm[pid] - e_qm[base], e_mm[pid] - e_mm[base]
    return dict(against=base, qm=round(q, 4), mm=round(m, 4), mm_minus_qm=round(m - q, 4), mm_split=rel_split(mm, pid, base))


def part_b(e_qm, e_mm, mm, info):
    out = dict(X_to_Rx=compare(e_qm, e_mm, mm, 'Rx', 'X'), Nx=compare(e_qm, e_mm, mm, 'Nx', 'X'), turns={})
    dev = []
    for p, i in info.items():
        if i['group'] == 'c':
            c = compare(e_qm, e_mm, mm, p, i['sets']['base'])
            if c:
                c['against_X'] = {k: v for k, v in (compare(e_qm, e_mm, mm, p, 'X') or {}).items() if k != 'mm_split'} or None
                c['flags'] = i['flags']
                out['turns'][p] = c
                dev.append(c['mm_minus_qm'])
    if dev:
        out['turns_rms_mm_minus_qm_against_base'] = round(float(np.sqrt(np.mean(np.square(dev)))), 4)
        out['turns_n'] = len(dev)
    return out


def relaxed_block(aux, with_d3):
    """E(relaxed from Rx) - E(relaxed from X) in both descriptions, from the files of the relaxations.  A
    quantum-chemical run counts as crossed only when it has converged with a signed volume at N2 of the sign
    opposite to its start (md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 1)."""
    rm, rq = read_json(aux['ff_dir'] / 'relax_mm.json'), read_json(aux['dir'] / 'relax_qm.json')
    out = dict(force_field=None, quantum_chemistry=None, mismatch=None)
    if rm:
        rx = rm['Rx_minus_X']
        out['force_field'] = dict(
            Rx_minus_X=rx['relaxed']['total'], split=rx['relaxed']['split'],
            configurations_kept=rm.get('configurations_kept'),
            signed_volume_N2={s: rm['starts'][s]['minimisation']['signed_volume_N2'] for s in rm['starts']},
            careful_descent=rm.get('careful_descent'),
            at_the_quantum_chemical_minima=rx.get('force_field_at_the_quantum_chemical_minima'),
            at_the_quantum_chemical_geometries=rx.get('force_field_at_the_quantum_chemical_geometries'),
            quantum_chemical_geometries={k: ({a: v.get(a) for a in (
                'file', 'sha256', 'largest_shift_of_a_fixed_atom_A', 'total', 'rmsd_ligand_heavy_from_X_start_A',
                'rmsd_ligand_heavy_from_the_force_field_minimum_from_X_A', 'above_the_force_field_minimum_from_X')}
                if v.get('exists') else None) for k, v in (rm.get('at_other_geometries') or {}).items()},
            file=str(aux['ff_dir'] / 'relax_mm.json'))
    if rq:
        key = 'Rx_minus_X_with_d3' if with_d3 else 'Rx_minus_X_without_d3'
        st = rq.get('starts', {})
        ep = rq.get('end_points')
        out['quantum_chemistry'] = dict(
            Rx_minus_X=rq.get(key), state={s: v.get('state') for s, v in st.items()},
            converged={s: v.get('converged', v.get('state') == 'done') for s, v in st.items()},
            crossed={s: v.get('crossed') for s, v in st.items()},
            nitrogen_at_the_end={s: v.get('nitrogen_at_the_end') for s, v in st.items()},
            configurations_kept=rq.get('configurations_kept'), both_converged=rq.get('both_converged'),
            signed_volume_N2={s: v.get('signed_volume_N2') for s, v in st.items()},
            n_steps={s: v.get('n_steps') for s, v in st.items()},
            end_points=(dict(ep, Rx_minus_X=ep.get(key)) if ep else None),
            rule=rq.get('rule'), file=str(aux['dir'] / 'relax_qm.json'))
    f, q = out['force_field'], out['quantum_chemistry']
    if f and q and q.get('Rx_minus_X') is not None and f.get('configurations_kept') and q.get('configurations_kept'):
        mis = f['Rx_minus_X'] - q['Rx_minus_X']
        out['mismatch'] = dict(mm_minus_qm=round(mis, 4), limit=RELAX_LIMIT, above_limit=bool(abs(mis) > RELAX_LIMIT))
    elif f and f.get('configurations_kept') is False:
        out['mismatch'] = dict(mm_minus_qm=None, note='in the force field the relaxation from Rx does not keep the R-like '
                               'configuration at N2, so there is no force-field difference between two relaxed configurations')
        if q and q.get('both_converged'):
            g = f.get('at_the_quantum_chemical_minima')
            out['mismatch'].update(
                quantum_chemistry_configurations_kept=q.get('configurations_kept'),
                at_the_quantum_chemical_end_points=(dict(
                    qm=q['Rx_minus_X'], mm=g['total'], mm_minus_qm=round(g['total'] - q['Rx_minus_X'], 4), mm_split=g['split'],
                    limit=RELAX_LIMIT, above_limit=bool(abs(g['total'] - q['Rx_minus_X']) > RELAX_LIMIT),
                    what='E(end point from Rx) - E(end point from X): quantum chemistry, and the force field as single points '
                         'at the two quantum-chemical end geometries') if (g and q.get('Rx_minus_X') is not None) else None))
    return out


def reading(qm, d3, mm, info, aux, scheme, with_d3):
    """Everything for one way of reading the quantum-chemical energies (with or without D3)."""
    ids = [p for p in qm if (not with_d3 or p in d3)]
    e_qm = {p: (qm[p]['e_tot'] + (d3[p] if with_d3 else 0.0)) * HARTREE for p in ids}
    e_mm = {p: mm['schemes'][scheme]['total'][p] for p in ids if p in mm['schemes'][scheme]['total']}
    a = part_a(e_qm, e_mm, info)
    a['without_protein_contact'] = {k: v for k, v in part_a(e_qm, e_mm, info, 'protein_contact').items()
                                    if k in ('n_points', 'n_expected', 'n_compared', 'rms', 'largest', 'largest_at')}
    b = part_b(e_qm, e_mm, mm, info)
    b['relaxed'] = relaxed_block(aux, with_d3)
    rep = {p: compare(e_qm, e_mm, mm, p, 'X') for p, i in info.items() if i['sets']['report_only']}
    return dict(part_A=a, part_B=b, reported_only=rep), e_qm, e_mm


def first_version(e_qm, e_mm, info):
    """The numbers of the gate as it stood until 8:56 AM (plan entries of 3:02 and 8:01 AM), kept for the record."""
    out = {}
    done = set(e_qm) & set(e_mm)

    def one(ids, need):
        ids = [p for p in ids if p in done]
        r = dict(n_points=len(ids), n_expected=len(need), complete=bool(set(need) <= done))
        if len(ids) < 2:
            return r
        ref = min(ids, key=lambda p: e_qm[p])
        d = {p: (e_mm[p] - e_mm[ref]) - (e_qm[p] - e_qm[ref]) for p in ids}
        v = np.array([d[p] for p in ids if p != ref])
        worst = max((p for p in ids if p != ref), key=lambda p: abs(d[p]))
        r.update(reference=ref, rms=round(float(np.sqrt((v ** 2).mean())), 4), largest=round(float(d[worst]), 4),
                 largest_at=worst, passed=bool(np.sqrt((v ** 2).mean()) <= GATE))
        return r
    for name in ('gate_plan', 'gate_rigid'):
        need = [p for p, i in info.items() if i['sets'].get(name)]
        out[name] = one(need, need)
    for a, b, key in (('X', 'R', 'dd_X_R_as_built'), ('X', 'Rx', 'dd_X_Rx')):
        if a in done and b in done:
            dd = (e_mm[b] - e_mm[a]) - (e_qm[b] - e_qm[a])
            out[key] = dict(qm=round(e_qm[b] - e_qm[a], 4), mm=round(e_mm[b] - e_mm[a], 4), dd=round(dd, 4),
                            passed=bool(abs(dd) <= GATE))
        else:
            out[key] = None
    for name, g, dd in (('verdict_plan', out['gate_plan'], out['dd_X_R_as_built']),
                        ('verdict_rigid', out['gate_rigid'], out['dd_X_Rx'])):
        out[name] = dict(complete=bool(g.get('complete') and dd is not None),
                         passed=(bool(g['passed'] and dd['passed']) if (g.get('rms') is not None and dd is not None) else None),
                         text='first version of the gate, quantum chemistry WITHOUT D3; superseded by the two parts')
    return out


def binding_block(mm, aux):
    b = read_json(aux['ff_dir'] / 'binding_check.json')
    if b:
        return dict(verdict=b.get('verdict'), comparison=b.get('comparison'), force_field=b.get('force_field'),
                    quantum_chemical_optimum=b.get('quantum_chemical_optimum'),
                    quantum_chemistry_against_crystal=b.get('quantum_chemistry_against_crystal'),
                    against_crystal_not_the_verdict=b.get('against_crystal_not_the_verdict'),
                    consequence=b.get('consequence'),
                    file=str(aux['ff_dir'] / 'binding_check.json'))
    if 'minimum_all_free' in mm:
        return dict(verdict='pending', note='binding_check.json not written yet (md/qm/gr_mm.py binding)',
                    against_crystal_not_the_verdict=mm['minimum_all_free'].get('binding_check'))
    return None


def cmd_gate(args):
    mm, pts, info, qm, d3, aux = load(args.mm, args.qm, args.side, args.ff_side)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    args.reading, args.step2 = reading_of(args, mm)
    own_file(out / f'{args.name}_gate.json', args.reading)
    used = mm['used_scheme']
    covered = [p for p in pts['order'] if p in mm['points']]
    missing_mm = [p for p in pts['order'] if p not in mm['points']]
    missing_d3 = [p for p in qm if p not in d3]
    res = dict(what='gate G-R: quantum chemistry against the force field on the small zinc model, read in two parts',
               plan=PLAN, reading=args.reading, step2=args.step2, written=now_ist(),
               mm=str(args.mm), mm_label=mm['label'], qm_dir=str(args.qm), used_scheme=used,
               d3=dict(file=aux['d3_file'], routine=aux['d3_routine'], checks_passed=aux['d3_checks'],
                       finished_points_without_d3=missing_d3),
               n_points=len(pts['order']), n_finished=len(qm), finished=[p for p in pts['order'] if p in qm],
               pending=[p for p in pts['order'] if p not in qm], points_not_in_the_mm_file=missing_mm,
               hosts=sorted({r.get('host', '?') for r in qm.values()}))
    wd, e_qm_d, e_mm = reading(qm, d3, mm, info, aux, used, True)
    wo, e_qm_0, _ = reading(qm, d3, mm, info, aux, used, False)
    res['with_d3'], res['without_d3'] = wd, wo
    a = wd['part_A']
    res['verdict_part_A'] = dict(
        complete=a['complete'], passed=(a.get('passed') if a['complete'] else None), rms=a.get('rms'),
        rms_without_d3=wo['part_A'].get('rms'), n_compared=a.get('n_compared'), n_expected=a['n_expected'],
        largest=a.get('largest'), largest_at=a.get('largest_at'), criterion=GATE,
        rms_single_reference=a.get('rms_single_reference'), rms_single_reference_without_d3=wo['part_A'].get('rms_single_reference'),
        subsets=subsets(e_qm_d, e_mm, info),
        text='part A (terms at N2): r.m.s. of force field minus quantum chemistry with D3(BJ), each family against its '
             'lowest point, at most 1.0 kcal/mol; so far only, while `complete` is false.  `rms_single_reference` (one '
             'reference for all points) is quoted beside and does not decide; `subsets`: the same number over parts of '
             'the set, reported only')
    res['part_B_is_not_a_gate'] = 'part B (contacts) is reported with the split by energy part; it has no verdict and is never fitted'
    res['part_A_composition'] = dict(
        what='force-field energy of every point of part A against the geometry it was made from (X or Rx), by what moves it; '
             'in the inner points of the inversion path (group a) and the turns about Zn-N2 (group e) more than the atoms at N2 '
             'moves, so part A tests the whole description of the nitrogen (fixed charges and 1-4 terms included), not the '
             'three angle terms alone',
        points=composition(e_qm_d, mm, info))
    res['part_A_by_cap_scheme'] = {}
    for s in mm['schemes']:
        r, _, _ = reading(qm, d3, mm, info, aux, s, True)
        res['part_A_by_cap_scheme'][s] = dict(rms=r['part_A'].get('rms'), rms_single_reference=r['part_A'].get('rms_single_reference'),
                                              largest=r['part_A'].get('largest'),
                                              X_to_Rx_mm_minus_qm=(r['part_B']['X_to_Rx'] or {}).get('mm_minus_qm'))
    res['cap_test'] = dict(largest_abs_two_reasonable=mm['cap_test_largest_abs_two_reasonable'],
                           largest_abs_all=mm['cap_test_largest_abs_all'])
    res['binding_check'] = binding_block(mm, aux)
    # the numbers of the first version of the gate (until 8:56 AM, without D3).  This block was also written under the
    # key `gate` for controllers started before 9 AM; none runs any more, and the key is gone since 6 Oct 2026, 1:30 PM
    # IST, so that nobody takes its pass flags for the verdict (the verdict is `verdict_part_A`)
    res['first_version_not_a_verdict'] = first_version(e_qm_0, e_mm, info)
    # ---- the table
    cols = ['order', 'id', 'group', 'part', 'base', 'flags', 'E_QM', 'E_MM', 'MM_minus_QM', 'vs_base_QM', 'vs_base_MM',
            'vs_base_MM_minus_QM', 'E_QM_no_D3', 'MM_minus_QM_no_D3', 'D3_rel_X', 'E_MM_group', 'E_MM_spread'] + \
           [f'MM_{k}' for k in SPLIT] + ['QM_hartree', 'D3_hartree', 'scf_cycles', 'wall_min', 'host']
    rows = [','.join(cols)]
    x_ok = 'X' in e_qm_d and 'X' in mm['points']
    for p in pts['points']:
        pid = p['id']
        s = p['sets']
        part = 'A' if s['part_A'] and pid != 'Rx' else ''
        part += 'B' if s['part_B'] else ''
        part = 'A+B' if pid in ('X', 'Rx') else (part or ('reported' if s['report_only'] else ''))
        r = {c: '' for c in cols}
        r.update(order=p['order'], id=pid, group=p['group'], part=part, base=s['base'], flags=' '.join(p['flags']))
        if pid in mm['points']:
            r.update(E_MM=round(mm['schemes'][used]['relative_energy'][pid], 3),
                     E_MM_group=round(mm['schemes']['group']['relative_energy'][pid], 3),
                     E_MM_spread=round(mm['schemes']['spread']['relative_energy'][pid], 3))
            sp = rel_split(mm, pid, 'X')
            if sp:
                r.update({f'MM_{k}': sp[k] for k in SPLIT})
        if pid in d3:
            r['D3_hartree'] = f'{d3[pid]:.9f}'
            if 'X' in d3:
                r['D3_rel_X'] = round((d3[pid] - d3['X']) * HARTREE, 3)
        if pid in qm:
            r.update(QM_hartree=f"{qm[pid]['e_tot']:.9f}", scf_cycles=qm[pid].get('scf_cycles', ''),
                     wall_min=round(qm[pid].get('wall_seconds', 0) / 60.0, 1), host=qm[pid].get('host', ''))
            if x_ok and pid in e_qm_d and r['E_MM'] != '':
                eq = e_qm_d[pid] - e_qm_d['X']
                r.update(E_QM=round(eq, 3), MM_minus_QM=round(r['E_MM'] - eq, 3))
                b = s['base']
                if b in e_qm_d and b in e_mm:
                    q, m = e_qm_d[pid] - e_qm_d[b], e_mm[pid] - e_mm[b]
                    r.update(vs_base_QM=round(q, 3), vs_base_MM=round(m, 3), vs_base_MM_minus_QM=round(m - q, 3))
            if 'X' in e_qm_0 and r['E_MM'] != '':
                e0 = e_qm_0[pid] - e_qm_0['X']
                r.update(E_QM_no_D3=round(e0, 3), MM_minus_QM_no_D3=round(r['E_MM'] - e0, 3))
        rows.append(','.join(str(r[c]) for c in cols))
    (out / f'{args.name}_table.csv').write_text('\n'.join(rows) + '\n')
    tmp = out / f'{args.name}_gate.json.tmp'
    tmp.write_text(json.dumps(res, indent=1) + '\n')
    os.replace(tmp, out / f'{args.name}_gate.json')
    print(summary_line(res))
    return res


def summary_line(res):
    v = res['verdict_part_A']
    b = res['with_d3']['part_B']
    xr, xr0 = b['X_to_Rx'] or {}, res['without_d3']['part_B']['X_to_Rx'] or {}
    rel = b['relaxed']
    bc = res.get('binding_check') or {}
    t = (f"G-R ({res['mm_label']}): {res['n_finished']}/{res['n_points']} points | part A (terms at N2): rms "
         f"{v['rms']} with D3 in two families ({v.get('rms_single_reference')} against a single reference, which does not decide), "
         f"{v['rms_without_d3']} without D3, {v['n_compared']} of {v['n_expected'] - len(res['with_d3']['part_A']['families'])} "
         f"comparisons, largest {v['largest']} at {v['largest_at']}, complete {v['complete']}, passed {v['passed']} | "
         f"part B (contacts, no verdict): X->Rx MM-QM {xr.get('mm_minus_qm')} with D3, {xr0.get('mm_minus_qm')} without; "
         f"turns rms {b.get('turns_rms_mm_minus_qm_against_base')} (n {b.get('turns_n')}); relaxed Rx-X MM "
         f"{(rel['force_field'] or {}).get('Rx_minus_X')} (configurations kept {(rel['force_field'] or {}).get('configurations_kept')}) "
         f"QM {(rel['quantum_chemistry'] or {}).get('Rx_minus_X')} | binding check: {bc.get('verdict')}")
    return t


# ------------------------------------------------------------------------------------------------------- the fit
def design(ic, fam, model, k_fixed):
    """Columns of the linear problem; `fam` gives the family index (0, 1, ...) of every point."""
    th = np.radians(np.array([[c[a] for a in ANGLES] for c in ic]))
    phi = np.radians(np.array([c['H_N2_S_C4'] for c in ic]))
    nf = max(fam) + 1
    cols = [np.array([1.0 if f == j else 0.0 for f in fam]) for j in range(nf)]
    if model == 'M2':
        cols += [th[:, j] ** 2 for j in range(3)]
    cols += [th[:, j] for j in range(3)] + [np.cos(phi), np.cos(2 * phi)]
    shift = np.zeros(len(ic))
    if model == 'M1':
        shift = (np.array([k_fixed[a] for a in ANGLES])[None, :] * th ** 2).sum(axis=1)
    return np.array(cols).T, shift, nf


def unpack(coef, nf, model, k_fixed):
    if model == 'M2':
        k = {a: float(coef[nf + j]) for j, a in enumerate(ANGLES)}
        b = coef[nf + 3:nf + 6]
        c1, c2 = float(coef[nf + 6]), float(coef[nf + 7])
    else:
        k = {a: float(k_fixed[a]) for a in ANGLES}
        b = coef[nf:nf + 3]
        c1, c2 = float(coef[nf + 3]), float(coef[nf + 4])
    th0 = {a: (float(np.degrees(-b[j] / (2.0 * k[a]))) if k[a] != 0 else float('nan')) for j, a in enumerate(ANGLES)}
    return dict(k=k, theta0=th0, c1=c1, c2=c2)


def solve(ic, fam, y, model, k_fixed):
    a, shift, nf = design(ic, fam, model, k_fixed)
    scale = np.abs(a).max(axis=0)
    scale[scale == 0] = 1.0
    coef, _, rank, sv = np.linalg.lstsq(a / scale, y - shift, rcond=None)
    coef = coef / scale
    return coef, a @ coef + shift, int(rank), sv, nf


def standard_errors(ic, fam, y, coef, pred, model, k_fixed):
    a, _, nf = design(ic, fam, model, k_fixed)
    n, p = a.shape
    if n <= p:
        return None
    scale = np.abs(a).max(axis=0)
    scale[scale == 0] = 1.0
    cov = np.linalg.pinv((a / scale).T @ (a / scale)) / np.outer(scale, scale) * float(((y - pred) ** 2).sum() / (n - p))
    se = dict(k={}, theta0={})
    for j, name in enumerate(ANGLES):
        jac = np.zeros(p)
        if model == 'M2':
            ik, ib = nf + j, nf + 3 + j
            k, b = coef[ik], coef[ib]
            jac[ik], jac[ib] = b / (2.0 * k * k), -1.0 / (2.0 * k)
            se['k'][name] = float(np.sqrt(max(cov[ik, ik], 0.0)))
        else:
            jac[nf + j] = -1.0 / (2.0 * k_fixed[name])
            se['k'][name] = 0.0
        se['theta0'][name] = float(np.degrees(np.sqrt(max(jac @ cov @ jac, 0.0))))
    i1 = nf + (6 if model == 'M2' else 3)
    se['c1'], se['c2'] = float(np.sqrt(max(cov[i1, i1], 0.0))), float(np.sqrt(max(cov[i1 + 1, i1 + 1], 0.0)))
    return se


def model_energy(ic, par):
    e = 0.0
    for a in ANGLES:
        e += par['k'][a] * np.radians(ic[a] - par['theta0'][a]) ** 2
    e += par['c1'] * np.cos(np.radians(ic['H_N2_S_C4'])) + par['c2'] * np.cos(2 * np.radians(ic['H_N2_S_C4']))
    return float(e)


def checks(f, structure, rms_before, rms_after):
    """Every condition of the acceptance rule, one by one.  `structure`: the centres of the 10-degree window of the
    equilibrium angles (window_centres(): Zn-N2-S of the crystal, Zn-N2-H and S-N2-H of the quantum-chemical optimum
    of the small model); a centre that is None (no optimum yet) cannot be met."""
    p, se = f['parameters'], f['standard_errors']
    loo = f['leave_one_out']['rms']
    known = all(structure.get(a) is not None for a in ANGLES)
    c = dict(
        force_constants_in_range=bool(all(K_RANGE[0] <= p['k'][a] <= K_RANGE[1] for a in ANGLES)),
        standard_errors_small=bool(se is not None and all(se['theta0'][a] <= SE_THETA_MAX for a in ANGLES) and
                                   all(se['k'][a] <= SE_K_REL_MAX * abs(p['k'][a]) for a in ANGLES)),
        residual_rms_at_most_1=bool(f['rms_residual'] <= GATE),
        leave_one_out_rms_at_most_1=bool(loo is not None and loo <= GATE),
        torsion_coefficients_below_5=bool(abs(p['c1']) < C_MAX and abs(p['c2']) < C_MAX),
        equilibrium_angles_within_10_deg_of_their_reference=bool(
            known and all(abs(p['theta0'][a] - structure[a]) <= THETA_WINDOW for a in ANGLES)),
        part_A_number_falls=bool(rms_before is not None and rms_after is not None and rms_after < rms_before),
        part_A_after_at_most_1=bool(rms_after is not None and rms_after <= GATE))
    return c, bool(all(c.values()))


def window_centres(info, aux=None, optimum='read'):
    """Centres of the 10-degree window of the fitted equilibrium angles (md/PLAN.md, 6 Oct 2026, 1:09 PM IST, point 4).
    Zn-N2-S: the crystal value (start X).  Zn-N2-H and S-N2-H: their values in the quantum-chemical optimum of the small
    model with all atoms free (binding_check.json, key quantum_chemical_optimum, written by md/qm/gr_mm.py binding from
    CA2_small_opt.fchk or, before that file exists, from the geometry of the finished optimisation), NOT the
    constructed hydrogen of start X.  Without the optimum the two centres are None and the condition refuses.
    `optimum`: 'read' (from aux), None (as if it did not exist), or a dict with Zn_N2_H and S_N2_H."""
    src = dict(Zn_N2_S='start X (crystal structure 6RL9)')
    val = dict(Zn_N2_S=info['X']['internal']['Zn_N2_S'], Zn_N2_H=None, S_N2_H=None)
    q, rec = None, None
    if optimum == 'read':
        b = read_json(aux['ff_dir'] / 'binding_check.json') if aux else None
        q = (b or {}).get('quantum_chemical_optimum')
        rec = dict(file=str(aux['ff_dir'] / 'binding_check.json') if aux else None, source=(q or {}).get('source'))
    elif optimum is not None:
        q, rec = optimum, dict(source='given to the function (a test)')
    for a in ('Zn_N2_H', 'S_N2_H'):
        if q and q.get(a) is not None:
            val[a] = float(q[a])
            src[a] = 'quantum-chemical optimum of the small model, all atoms free'
        else:
            src[a] = ('NOT AVAILABLE: the quantum-chemical optimum of the small model does not exist yet; the condition on the '
                      'equilibrium angles cannot be met and no fit is accepted until it does')
    return val, dict(sources=src, optimum=rec, constructed_hydrogen_of_start_X={a: info['X']['internal'][a] for a in ('Zn_N2_H', 'S_N2_H')},
                     complete=bool(all(v is not None for v in val.values())))


def fit_model(ids, fam_of, ic_all, y_all, model, k_fixed):
    ic = [ic_all[p] for p in ids]
    fam = [fam_of[p] for p in ids]
    y = np.array([y_all[p] for p in ids])
    coef, pred, rank, sv, nf = solve(ic, fam, y, model, k_fixed)
    par = unpack(coef, nf, model, k_fixed)
    res = y - pred
    se = standard_errors(ic, fam, y, coef, pred, model, k_fixed)
    nfree = 8 if model == 'M2' else 5
    loo = {}
    if len(ids) > nfree + nf + 1:
        for n, p in enumerate(ids):
            keep = [m for m in range(len(ids)) if m != n]
            if len({fam[m] for m in keep}) < nf:
                continue
            c2 = solve([ic[m] for m in keep], [fam[m] for m in keep], y[keep], model, k_fixed)[0]
            a1, s1, _ = design([ic[n]], [fam[n]], model, k_fixed)
            if a1.shape[1] != len(c2):                                     # the family column of a lone family
                a1 = np.hstack([np.eye(nf)[fam[n]][None, :], a1[:, 1:]])
            loo[p] = float(y[n] - (a1 @ c2 + s1)[0])
    return dict(model=model, n_points=len(ids), n_free_parameters=nfree, n_family_constants=nf, rank=rank,
                condition_number=float(sv.max() / sv.min()) if sv.min() > 0 else float('inf'),
                parameters=par, standard_errors=se,
                rms_residual=float(np.sqrt((res ** 2).mean())), largest_residual=float(np.abs(res).max()),
                residuals={p: round(float(r), 4) for p, r in zip(ids, res)},
                leave_one_out=dict(rms=(float(np.sqrt(np.mean([v ** 2 for v in loo.values()]))) if loo else None),
                                   largest=(float(max(abs(v) for v in loo.values())) if loo else None),
                                   per_point={p: round(v, 4) for p, v in loo.items()}))


def choose(m1, m2):
    l1, l2 = m1['leave_one_out']['rms'], (m2.get('leave_one_out') or {}).get('rms')
    if m2.get('acceptable') and l1 is not None and l2 is not None and l2 <= 0.8 * l1:
        return 'M2', 'M2 acceptable and its leave-one-out r.m.s. at least 20 % below that of M1'
    if m1.get('acceptable'):
        return 'M1', 'M1 acceptable; M2 not acceptable or not 20 % better in leave-one-out'
    imp = [m['model'] for m in (m1, m2) if m.get('improved_only')]
    if imp:
        return None, (f"{IMPROVED_ONLY}: {' and '.join(imp)} meet(s) every other condition and lower(s) part A, but part A with "
                      f"the fitted terms stays above {GATE}; not used, no parameter file, the terms stay as built")
    return None, 'neither model acceptable: no parameter file, the terms stay as built, the gate is reported as failed'


def amber_torsion(c, n):
    return dict(n=n, k=abs(c), phase=0.0 if c >= 0 else 180.0)


def write_frcmod(path, par, label):
    t1, t2 = amber_torsion(par['c1'], 1), amber_torsion(par['c2'], 2)
    lines = [f'Terms at N2 of the zinc-bound sulfonamidate fitted on B3LYP-D3(BJ)/6-31G* single points of the small zinc '
             f'model, gate G-R part A ({label}; md/qm/gr_fit.py); load after CA2_mcpbpy.frcmod',
             'MASS', '', 'BOND', '', 'ANGLE']
    for a in ANGLES:
        lines.append(f"{FRC_NAMES[a]}   {par['k'][a]:8.3f}   {par['theta0'][a]:9.3f}    fitted, gate G-R")
    lines += ['', 'DIHE',
              f"ca-sy-Y4-hn   1   {t1['k']:8.4f}   {t1['phase']:8.3f}   -1.0      fitted, gate G-R",
              f"ca-sy-Y4-hn   1   {t2['k']:8.4f}   {t2['phase']:8.3f}    2.0",
              '', 'IMPROPER', '', 'NONBON', '', '']
    pathlib.Path(path).write_text('\n'.join(lines))


def run_fit(e_qm, rest, total, ic, info, k_fixed, structure):
    """The whole decision for one set of energies (kcal/mol).  Returns the record and the chosen parameters."""
    fam = families(info)
    order = sorted(fam)
    fam_of = {p: order.index(b) for b, ids in fam.items() for p in ids}
    ids = [p for b in order for p in fam[b] if p in e_qm and p in rest]
    before = part_a(e_qm, {p: total[p] for p in e_qm if p in total}, info)
    rec = dict(fit_points=ids, families={b: [p for p in fam[b] if p in ids] for b in order}, part_A_before=before)
    if len(ids) < 12:
        rec.update(chosen=None, reason=f'only {len(ids)} points of part A are finished', outcome='no fit accepted')
        return rec, None
    y = {p: e_qm[p] - rest[p] for p in ids}
    models = {}
    for m in ('M1', 'M2'):
        if m == 'M2' and len(ids) < 16:
            models[m] = dict(model='M2', acceptable=False, leave_one_out=dict(rms=None), note='fewer than 16 points')
            continue
        f = fit_model(ids, fam_of, ic, y, m, k_fixed)
        e_after = {p: rest[p] + model_energy(ic[p], f['parameters']) for p in e_qm if p in rest}
        after = part_a(e_qm, e_after, info)
        f['part_A_after'] = dict(rms=after.get('rms'), rms_single_reference=after.get('rms_single_reference'),
                                 largest=after.get('largest'), largest_at=after.get('largest_at'),
                                 deviations=after.get('deviations'))
        f['checks'], f['acceptable'] = checks(f, structure, before.get('rms'), after.get('rms'))
        f['improved_only'] = bool(not f['acceptable'] and all(v for k, v in f['checks'].items() if k != 'part_A_after_at_most_1'))
        models[m] = f
    pick, why = choose(models['M1'], models['M2'])
    rec.update(M1=models['M1'], M2=models['M2'], chosen=pick, reason=why,
               outcome=('accepted' if pick else IMPROVED_ONLY if why.startswith(IMPROVED_ONLY) else 'no fit accepted'))
    par = models[pick]['parameters'] if pick else None
    if pick:
        rec['amber_terms'] = dict(
            angles={FRC_NAMES[a]: dict(k=round(par['k'][a], 3), theta0=round(par['theta0'][a], 3)) for a in ANGLES},
            torsion_ca_sy_Y4_hn=[amber_torsion(par['c1'], 1), amber_torsion(par['c2'], 2)])
    return rec, par


def rule_text():
    return dict(k_range=K_RANGE, largest_standard_error_theta0=SE_THETA_MAX, largest_relative_standard_error_k=SE_K_REL_MAX,
                residual_rms_max=GATE, leave_one_out_rms_max=GATE, torsion_coefficient_max=C_MAX,
                theta0_window_around_its_reference=THETA_WINDOW, part_A_with_the_fitted_terms_max=GATE,
                window_reference='Zn-N2-S: crystal value (start X); Zn-N2-H and S-N2-H: quantum-chemical optimum of the small '
                                 'model, all atoms free; without that optimum the condition refuses (plan, 6 Oct 2026, 1:09 PM IST, point 4)',
                text='acceptable: all conditions of `checks`; M2 when acceptable and leave-one-out r.m.s. <= 0.8 x that '
                     'of M1; else M1 when acceptable; else no parameter file and the gate is reported as failed.  A model that '
                     f'fails only on part A with the fitted terms being above {GATE} is reported as "{IMPROVED_ONLY}" and not used')


def structure_angles(info):
    return {a: info['X']['internal'][a] for a in ANGLES}


def cmd_fit(args):
    mm, pts, info, qm, d3, aux = load(args.mm, args.qm, args.side, args.ff_side)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    args.reading, args.step2 = reading_of(args, mm)
    own_file(out / 'fit.json', args.reading)
    used = mm['used_scheme']
    k_fixed = {a: mm['fit_terms_in_topology']['angles'][a]['k'] for a in ANGLES}
    structure, centres = window_centres(info, aux)
    ids = [p for p in qm if p in d3]
    e_qm = {p: (qm[p]['e_tot'] + d3[p]) * HARTREE for p in ids}
    ic = {p: mm['points'][p]['internal'] for p in mm['points']}
    res = dict(what='fit of the terms at N2 on the points of part A, quantum chemistry with D3(BJ) (md/qm/gr_fit.py)',
               plan=PLAN, reading=args.reading, step2=args.step2, written=now_ist(),
               mm=str(args.mm), mm_label=mm['label'], used_scheme=used, terms_before=mm['fit_terms_in_topology'],
               window_centres=structure, window_centres_record=centres, angles_in_start_X=structure_angles(info),
               rule=rule_text(), n_finished=len(qm),
               pending=[p for p in pts['order'] if p not in qm], schemes={})
    pars = {}
    for s in mm['schemes']:
        rec, par = run_fit(e_qm, mm['schemes'][s]['rest'], mm['schemes'][s]['total'], ic, info, k_fixed, structure)
        res['schemes'][s] = rec
        pars[s] = par
    u = res['schemes'][used]
    res['chosen'] = u['chosen']
    res['reason'] = u['reason']
    res['outcome'] = u.get('outcome')
    frc = out / 'CA2_n2_fit.frcmod'
    if u['chosen']:
        par = pars[used]
        write_frcmod(frc, par, f"model {u['chosen']}, cap scheme {used}, {len(u['fit_points'])} points")
        res.update(frcmod=str(frc), failed=False)
        res['spread_over_cap_schemes_same_model'] = {}
        for key in ('theta0', 'k'):
            for a in ANGLES:
                vals = [res['schemes'][s][u['chosen']]['parameters'][key][a] for s in res['schemes'] if u['chosen'] in res['schemes'][s]]
                res['spread_over_cap_schemes_same_model'][f'{key}_{a}'] = round(max(vals) - min(vals), 4)
        for key in ('c1', 'c2'):
            vals = [res['schemes'][s][u['chosen']]['parameters'][key] for s in res['schemes'] if u['chosen'] in res['schemes'][s]]
            res['spread_over_cap_schemes_same_model'][key] = round(max(vals) - min(vals), 4)
    else:
        frc.unlink(missing_ok=True)
        res.update(frcmod=None, failed=True,
                   consequence=('the terms at N2 stay as built; gate G-R part A is reported as failed, split by energy part' +
                                (f'; reported as "{IMPROVED_ONLY}"' if u.get('outcome') == IMPROVED_ONLY else '')))
        # the split of the deviations of part A: which part of the force-field energy carries each of them
        split = {}
        for p, dv in (u.get('part_A_before') or {}).get('deviations', {}).items():
            b = info[p]['sets']['base']
            ref = u['part_A_before']['families'][b]['reference']
            split[p] = dict(against=ref, mm_minus_qm=dv, mm=round(mm['schemes'][used]['total'][p] - mm['schemes'][used]['total'][ref], 4),
                            qm=round(e_qm[p] - e_qm[ref], 4), mm_split=rel_split(mm, p, ref))
        res['part_A_split'] = split
    if 'M1' in u:
        rows = ['id,group,base,residual_M1,residual_M2,leave_one_out_M1,leave_one_out_M2,part_A_before,part_A_after_M1,part_A_after_M2']
        for p in u['fit_points']:
            rows.append(','.join(str(t) for t in (
                p, info[p]['group'], info[p]['sets']['base'], u['M1']['residuals'].get(p, ''),
                (u['M2'].get('residuals') or {}).get(p, ''), u['M1']['leave_one_out']['per_point'].get(p, ''),
                (u['M2'].get('leave_one_out') or {}).get('per_point', {}).get(p, ''),
                u['part_A_before']['deviations'].get(p, ''), u['M1']['part_A_after']['deviations'].get(p, ''),
                (u['M2'].get('part_A_after') or {}).get('deviations', {}).get(p, ''))))
        (out / 'fit_residuals.csv').write_text('\n'.join(rows) + '\n')
    tmp = out / 'fit.json.tmp'
    tmp.write_text(json.dumps(res, indent=1) + '\n')
    os.replace(tmp, out / 'fit.json')
    print(f"fit ({mm['label']}, scheme {used}): chosen {u['chosen']} ({u['reason']}); outcome: {u.get('outcome')}; window centres "
          f"{ {a: (round(v, 2) if v is not None else None) for a, v in structure.items()} }")
    for m in ('M1', 'M2'):
        r = u.get(m) or {}
        if 'parameters' in r:
            p = r['parameters']
            print(f"  {m}: residual {r['rms_residual']:.3f} loo {r['leave_one_out']['rms']} part A {u['part_A_before'].get('rms')} -> "
                  f"{r['part_A_after']['rms']} acceptable {r['acceptable']} theta0 {[round(p['theta0'][a], 2) for a in ANGLES]} "
                  f"k {[round(p['k'][a], 1) for a in ANGLES]} c1 {p['c1']:.3f} c2 {p['c2']:.3f}; failed checks "
                  f"{[k for k, v in r['checks'].items() if not v]}")
    return res


# ------------------------------------------------------------------------------------------------------ self-test
def cmd_selftest(args):
    mm, pts, info, _, d3, aux = load(args.mm, None, args.side, args.ff_side)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    used = mm['used_scheme']
    rest, total = mm['schemes'][used]['rest'], mm['schemes'][used]['total']
    ic = {p: mm['points'][p]['internal'] for p in rest}
    k_top = {a: mm['fit_terms_in_topology']['angles'][a]['k'] for a in ANGLES}
    structure, centres = window_centres(info, aux)
    stand_in = not centres['complete']
    if stand_in:                                                           # no optimum yet: the test needs centres all the same
        structure, centres = window_centres(info, optimum={a: info['X']['internal'][a] for a in ('Zn_N2_H', 'S_N2_H')})
        centres['stand_in'] = ('the quantum-chemical optimum does not exist yet: for this self-test ONLY the constructed hydrogen '
                               'angles of start X stand in for it')
    no_optimum, _ = window_centres(info, optimum=None)
    order = [p for p in pts['order'] if p in rest]
    top = mm['fit_terms_in_topology']
    c_top = {1: 0.0, 2: 0.0}
    for t in top['torsion_H_N2_S_C4']:
        c_top[t['n']] = t['k'] if abs(t['phase']) < 1.0 else -t['k']
    truth_top = dict(k=dict(k_top), theta0={a: top['angles'][a]['theta0'] for a in ANGLES}, c1=c_top[1], c2=c_top[2])
    # the force field's own terms in the form of the model: must reproduce total - rest up to a constant
    const = [total[p] - rest[p] - model_energy(ic[p], truth_top) for p in order]
    form_ok = bool(max(const) - min(const) < 1e-3)
    rng = np.random.default_rng(SEED)
    res = dict(what='self-test of gate and fit on synthetic energies (md/qm/gr_fit.py selftest)', plan=PLAN, mm=str(args.mm),
               mm_label=mm['label'], written=now_ist(), seed=SEED, noise_sigma=0.3, n_draws=200, rule=rule_text(),
               window_centres=structure, window_centres_record=centres, angles_in_start_X=structure_angles(info),
               terms_of_the_topology=truth_top,
               model_form_reproduces_the_topology_terms=dict(spread_of_the_constant=max(const) - min(const), passed=form_ok),
               cases={})
    ok = form_ok

    def decide(e, centres_=None):
        rec, par = run_fit(e, rest, total, ic, info, k_top, structure if centres_ is None else centres_)
        b = rec['part_A_before']
        fitted = bool(b.get('passed') is False)                            # the fit is run only when part A fails
        return rec, par, b, fitted

    def noise_runs(e0, expect):
        counts = dict(gate_passed_no_fit=0, M1=0, M2=0, none=0)
        pars = []
        for _ in range(200):
            e = {p: v + rng.normal(0.0, 0.3) for p, v in e0.items()}
            rec, par, b, fitted = decide(e)
            if not fitted:
                counts['gate_passed_no_fit'] += 1
            else:
                counts[rec['chosen'] or 'none'] += 1
                if par:
                    pars.append([par['theta0'][a] for a in ANGLES] + [par['k'][a] for a in ANGLES] + [par['c1'], par['c2']])
        blk = dict(outcomes_of_200=counts, expected=expect)
        if pars:
            arr = np.array(pars)
            names = [f'theta0_{a}' for a in ANGLES] + [f'k_{a}' for a in ANGLES] + ['c1', 'c2']
            blk['released_parameters'] = dict(mean={n: round(float(v), 3) for n, v in zip(names, arr.mean(axis=0))},
                                              sd={n: round(float(v), 3) for n, v in zip(names, arr.std(axis=0, ddof=1) if len(arr) > 1 else arr.std(axis=0))})
        return blk

    def brief(rec, par, b, fitted):
        d = dict(part_A_before=dict(rms=b.get('rms'), passed=b.get('passed'), rms_single_reference=b.get('rms_single_reference')),
                 fit_run=fitted, chosen=(rec.get('chosen') if fitted else None), outcome=(rec.get('outcome') if fitted else None))
        if fitted:
            for m in ('M1', 'M2'):
                if 'parameters' in (rec.get(m) or {}):
                    d[m] = dict(parameters=rec[m]['parameters'], rms_residual=rec[m]['rms_residual'],
                                leave_one_out_rms=rec[m]['leave_one_out']['rms'], part_A_after=rec[m]['part_A_after']['rms'],
                                acceptable=rec[m]['acceptable'], improved_only=rec[m].get('improved_only'),
                                failed_checks=[k for k, v in rec[m]['checks'].items() if not v])
        return d

    def offset(d_s, d_h1, d_h2):
        return dict(Zn_N2_S=structure['Zn_N2_S'] + d_s, Zn_N2_H=structure['Zn_N2_H'] + d_h1, S_N2_H=structure['S_N2_H'] + d_h2)

    def energies(truth, extra=None):
        e = {p: rest[p] + model_energy(ic[p], truth) + 7.0 for p in order}
        for p, v in (extra or {}).items():
            e[p] += v
        return e

    # T0: the force field itself
    e0 = {p: total[p] + 7.0 for p in order}
    r = decide(e0)
    case = dict(what='energies = the force field itself; right answer: the gate passes, nothing is fitted',
                without_noise=brief(*r), with_noise=noise_runs(e0, 'gate passes in (nearly) every draw; no parameter file'))
    case['passed'] = bool(r[2].get('passed') is True and r[2].get('rms') < 1e-6 and
                          case['with_noise']['outcomes_of_200']['M1'] + case['with_noise']['outcomes_of_200']['M2'] == 0)
    res['cases']['T0_change_nothing'] = case
    ok = ok and case['passed']
    # T0b: the reviewer's test: force field minus D3(BJ), read WITHOUT the D3 term
    if all(p in d3 for p in order):
        eb = {p: total[p] - d3[p] * HARTREE for p in order}
        r = decide(eb)
        r_with = decide({p: eb[p] + d3[p] * HARTREE for p in order})
        case = dict(what="the reviewer's test (review of 6 Oct 2026, fit_sensitivity case B): energies = force field minus "
                         'D3(BJ); read without the D3 term the right answer is "change nothing"; read with it the gate '
                         'must pass exactly',
                    read_without_d3=brief(*r), read_with_d3=brief(*r_with))
        case['passed'] = bool((not r[3] or r[0].get('chosen') is None) and r_with[2].get('passed') is True and
                              r_with[2].get('rms') < 1e-6)
        case['nothing_changed'] = bool(not r[3] or r[0].get('chosen') is None)
        res['cases']['T0b_reviewer_force_field_minus_D3'] = case
        ok = ok and case['passed']
    # T0c: a contact error on everything made from Rx, and on the turns
    ec = {p: total[p] + (3.4 if info[p]['sets']['base'] == 'Rx' else 0.0) + (2.0 if info[p]['group'] == 'c' else 0.0) for p in order}
    r = decide(ec)
    case = dict(what='energies = force field + 3.4 kcal/mol on every point made from Rx + 2.0 on every turn about N2-S '
                     '(errors of part B); right answer: part A passes, nothing is fitted',
                without_noise=brief(*r))
    case['passed'] = bool(r[2].get('passed') is True and r[2].get('rms') < 1e-6)
    res['cases']['T0c_contact_error_only'] = case
    ok = ok and case['passed']
    # T1, T2: other parameters inside the rule; T3: outside it.  The equilibrium angles are given as offsets from the
    # centres of the window (3.45, 4.59, 2.59 deg: with the centres of the rule before 1:09 PM these were 121, 112, 110)
    truths = dict(
        T1_other_angles_and_torsion_constants_of_the_topology=dict(k=dict(k_top), theta0=offset(3.45, 4.59, 2.59), c1=0.8, c2=-1.5),
        T2_other_constants_as_well=dict(k=dict(Zn_N2_S=38.0, Zn_N2_H=28.0, S_N2_H=45.0), theta0=offset(3.45, 4.59, 2.59), c1=0.8, c2=-1.5),
        T3_outside_the_rule=dict(k=dict(k_top), theta0=offset(13.45, 4.59, 2.59), c1=0.8, c2=-6.5))
    for name, truth in truths.items():
        et = energies(truth)
        rec, par, b, fitted = decide(et)
        inside = not name.startswith('T3')
        case = dict(truth=truth, inside_the_rule=inside, without_noise=brief(rec, par, b, fitted))
        if inside:
            exact = {}
            for m in ('M1', 'M2'):
                p = (rec.get(m) or {}).get('parameters')
                if p:
                    exact[m] = dict(theta0=max(abs(p['theta0'][a] - truth['theta0'][a]) for a in ANGLES),
                                    k=max(abs(p['k'][a] - truth['k'][a]) for a in ANGLES),
                                    c=max(abs(p['c1'] - truth['c1']), abs(p['c2'] - truth['c2'])))
            case['largest_errors_without_noise'] = exact
            can = ['M2'] if name.startswith('T2') else ['M1', 'M2']
            good = bool(fitted and rec.get('chosen') in can and
                        all(exact[m]['theta0'] < 1e-5 and exact[m]['k'] < 1e-5 and exact[m]['c'] < 1e-7 for m in can) and
                        rec[rec['chosen']]['part_A_after']['rms'] < 1e-6)
            case['passed'] = good
            case['with_noise'] = noise_runs(et, 'a parameter set is released in most draws and scatters about the truth')
        else:
            case['passed'] = bool(fitted and rec.get('chosen') is None)
            case['with_noise'] = noise_runs(et, 'no parameter file in (nearly) every draw')
        res['cases'][name] = case
        ok = ok and case['passed']
    # ---- the two additions of the plan entry of 6 Oct 2026, 1:09 PM IST, point 4
    window = 'equilibrium_angles_within_10_deg_of_their_reference'
    t1 = truths['T1_other_angles_and_torsion_constants_of_the_topology']
    # T4: the same data as T1, read while the quantum-chemical optimum does not exist: refused, on that condition alone
    r = decide(energies(t1), no_optimum)
    w = brief(*r)
    case = dict(what='the data of T1 (inside the rule), read WITHOUT the quantum-chemical optimum of the small model: the window of '
                     'the two hydrogen angles has no centre; right answer: no parameter file, refused on that condition alone',
                window_centres=no_optimum, without_noise=w)
    case['passed'] = bool(r[3] and r[0].get('chosen') is None and all(w[m]['failed_checks'] == [window] for m in ('M1', 'M2') if m in w))
    res['cases']['T4_without_the_quantum_chemical_optimum'] = case
    ok = ok and case['passed']
    # T5: a fit that lowers part A but leaves it above 1.0 (the case of the check of 6 Oct 2026, finding 4): X and the four
    # turns about Zn-N2 have the same fitted coordinates, so no fitted term can change their mutual deviations
    far = dict(k=dict(k_top), theta0=offset(6.45, 5.59, 4.59), c1=1.5, c2=-2.5)
    same = [p for p in ('X', 'e_p015', 'e_p030', 'e_m120', 'e_p120') if p in order]
    e5 = energies(far)
    low = min(same, key=lambda p: e5[p])
    r = decide(energies(far, {low: -2.4}))
    w = brief(*r)
    case = dict(what=f'parameters inside the rule; one of the five points with identical fitted coordinates ({low}) lowered by 2.4 '
                     f'kcal/mol; right answer: "{IMPROVED_ONLY}", no parameter file',
                truth=far, lowered=low, without_noise=w, reason=r[0].get('reason'))
    case['passed'] = bool(r[3] and r[0].get('chosen') is None and r[0].get('outcome') == IMPROVED_ONLY and
                          any(w[m]['failed_checks'] == ['part_A_after_at_most_1'] for m in ('M1', 'M2') if m in w))
    res['cases']['T5_terms_improved_gate_still_failed'] = case
    ok = ok and case['passed']
    # T6: Zn-N2-H 11 deg from the optimum (exact data): refused on the window alone
    t6 = dict(k=dict(k_top), theta0=offset(3.45, 11.0, 2.59), c1=0.8, c2=-1.5)
    r = decide(energies(t6))
    w = brief(*r)
    case = dict(what='exact data; equilibrium value of Zn-N2-H 11 deg from the centre of its window; right answer: no parameter file, '
                     'refused on that condition alone', truth=t6, without_noise=w)
    case['passed'] = bool(r[3] and r[0].get('chosen') is None and w['M1']['failed_checks'] == [window])
    res['cases']['T6_Zn_N2_H_11_deg_from_its_centre'] = case
    ok = ok and case['passed']
    # T7: both hydrogen angles at the centres themselves (with the real optimum: 15 deg from the constructed hydrogen)
    t7 = dict(k=dict(k_top), theta0=offset(3.45, 0.0, 0.0), c1=0.8, c2=-1.5)
    r = decide(energies(t7))
    w = brief(*r)
    case = dict(what='exact data; equilibrium values of Zn-N2-H and S-N2-H at the centres of their windows (the quantum-chemical '
                     'optimum); right answer: accepted, parameters back',
                truth=t7, distance_from_the_constructed_hydrogen={a: round(structure[a] - info['X']['internal'][a], 3) for a in ('Zn_N2_H', 'S_N2_H')},
                without_noise=w)
    case['passed'] = bool(r[3] and r[0].get('chosen') in ('M1', 'M2') and
                          max(abs(r[1]['theta0'][a] - t7['theta0'][a]) for a in ANGLES) < 1e-5)
    res['cases']['T7_hydrogen_angles_at_the_optimum'] = case
    ok = ok and case['passed']
    res['passed'] = bool(ok)
    (out / 'fit_selftest.json').write_text(json.dumps(res, indent=1) + '\n')
    print('window centres:', structure, '| stand-in for the optimum:', stand_in)
    print('the model form reproduces the terms of the topology:', res['model_form_reproduces_the_topology_terms'])
    for name, case in res['cases'].items():
        print(name, 'PASSED' if case['passed'] else 'FAILED')
        for key in ('without_noise', 'read_without_d3', 'read_with_d3'):
            if key in case:
                w = case[key]
                line = f"  {key}: part A before {w['part_A_before']}, fit run {w['fit_run']}, chosen {w['chosen']}, outcome {w.get('outcome')}"
                for m in ('M1', 'M2'):
                    if m in w:
                        line += (f"; {m} residual {w[m]['rms_residual']:.3f} loo {w[m]['leave_one_out_rms']} part A after "
                                 f"{w[m]['part_A_after']} acceptable {w[m]['acceptable']} failed {w[m]['failed_checks']}")
                print(line)
        if 'with_noise' in case:
            print('  with noise 0.3:', case['with_noise']['outcomes_of_200'], '| sd of released parameters:',
                  (case['with_noise'].get('released_parameters') or {}).get('sd'))
    print('SELF-TEST PASSED' if ok else 'SELF-TEST FAILED')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('gate', 'fit', 'selftest'):
        s = sub.add_parser(name)
        s.add_argument('--mm', required=True)
        s.add_argument('--out', required=True)
        s.add_argument('--side', default=None, help='folder with d3.json, relax_mm.json, relax_qm.json, binding_check.json (default: that of the points file)')
        s.add_argument('--ff-side', default=None, help='folder with relax_mm.json and binding_check.json of this reading (default: --side)')
        s.add_argument('--reading', default=None, help='name of the reading (default: the step named in the label of MM.json)')
        s.add_argument('--step2', default=None)
        if name != 'selftest':
            s.add_argument('--qm', required=True)
        if name == 'gate':
            s.add_argument('--name', default='gr')
    args = ap.parse_args()
    if args.cmd == 'gate':
        cmd_gate(args)
    elif args.cmd == 'fit':
        cmd_fit(args)
    else:
        sys.exit(cmd_selftest(args))


if __name__ == '__main__':
    main()
