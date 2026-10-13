"""Gate G-R2: the shortest contact between the ligand and the three zinc histidines that the rehearsal runs
visited (md/N2.md, sections 5.1 and 5.2: the contact limit of the new points).

    /root/cqc-sfa-tools/md/bin/python md/qm/n2_contacts.py      (MDAnalysis, NumPy; one thread; under a minute)

The same 2,100 solute frames as md/qm/n2_ranges.py (equilibration stages 1 to 6 and the 5 ns of production of
TX0, TR0, TN0; /root/cqc-sfa-md/runs_test/bonded_2e/).  For every frame, in the coordinates of the run itself:
the shortest distance between a heavy atom of the ligand and a heavy atom of His94, His96 or His119 that the
small model contains (CA, CB, CG, ND1, CD2, CE1, NE2), the quantity `heavy_heavy` of gr_points.Model.contacts
that the exclusion rule of md/qm/n2_points.py uses.  No energy is read.
Added 6 Oct 2026, 3:10 PM IST (after the first fitting energies showed where the fit misses): the same for the
contacts with hydrogens, heavy atom to hydrogen and hydrogen to hydrogen (all atoms of the three side chains from
CA outwards against all ligand atoms), and the number of new points below what the runs visited and below the
limits gate G-R used for its own turns (2.1 and 1.9 A).  Reported; no point is removed after its energy exists.

Writes md/qm/n2/contacts_visited.json.
"""
import json
import pathlib

import numpy as np

import n2_common as nc

RUNS = pathlib.Path('/root/cqc-sfa-md/runs_test/bonded_2e')
STAGES = ['eq01', 'eq02', 'eq03', 'eq04', 'eq05', 'eq06', 'seg0001']
SIDE = ('CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2')


def main():
    import MDAnalysis as mda
    out = dict(what='shortest heavy-atom contact ligand to zinc histidines (atoms of the small model) in the frames of '
                    'the three rehearsal runs; written by md/qm/n2_contacts.py', written=nc.now(), runs={},
               stages=STAGES, histidine_atoms=SIDE, mdanalysis=mda.__version__)
    allmin, all_hh, all_h2 = [], [], []
    for run in ('TX0', 'TR0', 'TN0'):
        d = RUNS / run
        u = mda.Universe(str(d / 'solute.pdb'), [str(d / f'solute_{s}.dcd') for s in STAGES])
        lig = u.select_atoms('resname SN1 and not name H*')
        his = u.select_atoms('resid 93 95 118 and name ' + ' '.join(SIDE))
        assert len(lig) == 11 and len(his) == 21, (len(lig), len(his))
        lig_all = u.select_atoms('resname SN1')
        his_all = u.select_atoms('resid 93 95 118 and not name N H C O')      # the atoms the small model has (CA, HA outwards)
        lh = np.array([a.name.startswith('H') for a in lig_all])
        hh = np.array([a.name.startswith('H') for a in his_all])
        assert len(lig_all) == 18 and int(hh.sum()) + 21 == len(his_all), (len(lig_all), len(his_all))
        mins, pair, m_hh, m_h2 = [], {}, [], []
        for ts in u.trajectory:
            dm = np.linalg.norm(lig.positions[:, None, :].astype(float) - his.positions[None, :, :].astype(float), axis=2)
            i, j = np.unravel_index(dm.argmin(), dm.shape)
            mins.append(float(dm[i, j]))
            k = f'{lig[i].name}...{his[j].resname}{his[j].resid}:{his[j].name}'
            pair[k] = min(pair.get(k, 9e9), float(dm[i, j]))
            da = np.linalg.norm(lig_all.positions[:, None, :].astype(float) - his_all.positions[None, :, :].astype(float), axis=2)
            nh = lh[:, None].astype(int) + hh[None, :].astype(int)
            m_hh.append(float(da[nh == 1].min()))
            m_h2.append(float(da[nh == 2].min()))
        mins = np.array(mins)
        allmin += list(mins)
        all_hh += m_hh
        all_h2 += m_h2
        out['runs'][run] = dict(frames=len(mins), shortest=round(float(mins.min()), 3),
                                first_percentile=round(float(np.percentile(mins, 1)), 3),
                                median=round(float(np.median(mins)), 3),
                                closest_pairs={k: round(v, 3) for k, v in sorted(pair.items(), key=lambda kv: kv[1])[:4]},
                                heavy_to_hydrogen=dict(shortest=round(min(m_hh), 3), first_percentile=round(float(np.percentile(m_hh, 1)), 3)),
                                hydrogen_to_hydrogen=dict(shortest=round(min(m_h2), 3), first_percentile=round(float(np.percentile(m_h2, 1)), 3)))
    allmin = np.array(allmin)
    pj = json.loads((nc.N2 / 'points.json').read_text())
    hh = [p['contacts']['model']['heavy_heavy']['d'] for p in pj['points']]
    out['all_runs'] = dict(frames=len(allmin), shortest=round(float(allmin.min()), 3),
                           first_percentile=round(float(np.percentile(allmin, 1)), 3))
    out['limit_in_n2_points'] = pj['limits']['contact_heavy_heavy_A']
    out['new_points'] = dict(n=len(hh), shortest=min(hh), below_the_shortest_visited=int(sum(h < allmin.min() for h in hh)))
    split = nc.split_record()
    vis = dict(heavy_h=min(all_hh), h_h=min(all_h2))
    out['contacts_with_hydrogens'] = dict(
        visited=dict(heavy_to_hydrogen_shortest=round(vis['heavy_h'], 3), heavy_to_hydrogen_first_percentile=round(float(np.percentile(all_hh, 1)), 3),
                     hydrogen_to_hydrogen_shortest=round(vis['h_h'], 3), hydrogen_to_hydrogen_first_percentile=round(float(np.percentile(all_h2, 1)), 3)),
        limits_of_gate_GR_for_its_turns=dict(heavy_to_hydrogen=2.1, hydrogen_to_hydrogen=1.9), new_points={})
    for label, sel in (('main_list', [p for p in pj['points'] if p['list'] == 'main']), ('main_list_fit', [p for p in pj['points'] if p['list'] == 'main' and split['set'][p['id']] == 'fit']),
                       ('main_list_test', [p for p in pj['points'] if p['list'] == 'main' and split['set'][p['id']] == 'test']), ('all', pj['points']),
                       ('all_fit', [p for p in pj['points'] if split['set'][p['id']] == 'fit']), ('all_test', [p for p in pj['points'] if split['set'][p['id']] == 'test'])):
        a = [p['contacts']['model']['heavy_h']['d'] for p in sel]
        b = [p['contacts']['model']['h_h']['d'] for p in sel]
        out['contacts_with_hydrogens']['new_points'][label] = dict(
            n=len(sel), heavy_to_hydrogen_shortest=min(a), hydrogen_to_hydrogen_shortest=min(b),
            below_the_limits_of_gate_GR=int(sum(x < 2.1 or y < 1.9 for x, y in zip(a, b))),
            below_what_the_runs_visited=int(sum(x < vis['heavy_h'] or y < vis['h_h'] for x, y in zip(a, b))))
    ho = [p for p in json.loads((nc.GR / 'points.json').read_text())['points'] if p['id'] not in ('R', 'N')]
    out['contacts_with_hydrogens']['old_rigid_points'] = dict(
        heavy_to_hydrogen_shortest=min(p['contacts']['model']['heavy_h']['d'] for p in ho),
        hydrogen_to_hydrogen_shortest=min(p['contacts']['model']['h_h']['d'] for p in ho))
    old = json.loads((nc.GR / 'points.json').read_text())['points']
    hh_old = {p['id']: p['contacts']['model']['heavy_heavy']['d'] for p in old if p['id'] not in ('R', 'N')}
    out['old_rigid_points'] = dict(n=len(hh_old), shortest=min(hh_old.values()), at=min(hh_old, key=hh_old.get))
    nc.write_json(nc.N2 / 'contacts_visited.json', out)
    print(json.dumps({k: v for k, v in out.items() if k in ('all_runs', 'limit_in_n2_points', 'new_points',
                                                           'old_rigid_points', 'contacts_with_hydrogens')}, indent=1))


if __name__ == '__main__':
    main()
