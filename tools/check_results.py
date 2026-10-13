#!/usr/bin/env python3
"""Cross-check selected stored results against raw outputs and printed tables.

Uses only the standard library. Does not run a scientific calculation.
"""
import csv
from decimal import Decimal
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DFT = ROOT / 'supporting/06_DFT_results'
TABLES = ROOT / 'supporting/05_Tables/SI'
HA_KCAL = Decimal('627.5094740631')
RY_KJ = Decimal('13.605693122994') * Decimal('96.48533212331002')


def read(path):
    return json.loads(path.read_text(), parse_float=Decimal)


def close(a, b, tolerance, label):
    if abs(Decimal(str(a)) - Decimal(str(b))) > Decimal(str(tolerance)):
        raise ValueError('Stored result disagreement: ' + label)


def table(name):
    with (TABLES / name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.reader(stream))[1:]


def qe_energy(run):
    p = DFT / 'periodic_QE/runs' / run / (run + '.out')
    matches = re.findall(r'(?m)^\s*!\s+total energy\s*=\s*([-+\d.Ee]+)\s+Ry', p.read_text())
    if not matches: raise ValueError('No energy: ' + run)
    return Decimal(matches[-1])


def main():
    scan = {}
    for p in (DFT / 'proton_scan').glob('[AB]_*.json'):
        d = read(p)
        if 'E_Ha' in d: scan[(d['chain'], str(d['index']))] = d
    reference = scan[('A', '0')]['E_Ha']
    for row in table('Table_S07_scan.csv'):
        d = scan[(row[0], row[1])]
        close((d['E_Ha'] - reference) * HA_KCAL, row[8], '.0051', 'Table S7 energy')
        close(d['metrics']['d_N1_H1A_A'], row[4], '.00051', 'Table S7 distance')
        if d['converged'] != row[9].startswith('yes'):
            raise ValueError('Table S7 convergence flag differs')
    if scan[('A', '2')]['converged']:
        raise ValueError('Unconverged scan point lost its flag')
    methods = read(DFT / 'functional_check/fsp_scan.json')
    keys = {'B3LYP-D3(BJ)': 'b3lyp-d3bj', 'ωB97X-D': 'wb97x-d', 'PBE0-D3(BJ)': 'pbe0-d3bj', 'M06-2X': 'm06-2x', 'DF-MP2': 'df-mp2'}
    for row in table('Table_S09_fcheck.csv'):
        e = methods['energies'][keys[row[0]]]
        n, t, o = (e[g]['E_Ha'] for g in ['end_N_d1.098', 'top1_d1.250', 'end_O_d1.631'])
        for actual, printed in zip([(t - n) * HA_KCAL, (o - n) * HA_KCAL, (t - o) * HA_KCAL], row[1:4]):
            close(actual, printed, '.0051', 'Table S9 energy difference')
    ladder = read(DFT / 'eps_ladder/eps_ladder.json')['rows']
    for data, printed in zip(ladder, table('Table_S08_eps.csv')):
        close(data['metrics']['d_N1_H1A_A'], printed[2], '.00051', 'Table S8 distance')
        close(data['metrics']['d_H1A_O5_A'], printed[3], '.00051', 'Table S8 distance')
    profile = read(DFT / 'periodic_QE/profile.json')
    salt = qe_energy('relax_i')
    close(salt, profile['relax_i']['energy_Ry'], '1e-8', 'Salt QE output')
    for row in profile['profile']:
        point = row['point']
        e = qe_energy('relax_' + point)
        close(e, profile['points'][point]['energy_Ry'], '1e-8', 'Profile QE output')
        pbe = (e - salt) * RY_KJ / 2
        close(pbe, row['E_PBE_D3BJ_kJ_mol_per_fu'], '.0006', 'PBE per-formula-unit conversion')
        hybrid_delta = qe_energy('hybrid_' + point + '_sg15_e100') - qe_energy('hybrid_i_sg15_e100')
        pbe_delta = qe_energy('pbe_' + point + '_sg15_e100') - qe_energy('pbe_i_sg15_e100')
        best = pbe + (hybrid_delta - pbe_delta) * RY_KJ / 2
        close(best, row['E_PBE0_D3BJ_kJ_mol_per_fu'], '.001', 'PBE0 hybrid correction')
    with (ROOT / 'data/derived/md_run_index.csv').open() as stream:
        runs = list(csv.DictReader(stream))
    roles = {r: sum(x['role'] == r for x in runs) for r in ['primary', 'sensitivity', 'validation']}
    if roles != {'primary': 9, 'sensitivity': 2, 'validation': 2}:
        raise ValueError('MD roles differ from the manuscript design')
    for row in runs:
        for key in ['series_path', 'final_frame_path']:
            if row[key] and not (ROOT / row[key]).is_file(): raise ValueError('Missing indexed MD file')
    relax = ROOT / 'molecular_dynamics/quantum_chemistry/gate_G-R/relax'
    old = read(relax / 'qm_Rx.json'); new = read(relax / 'qm_Rx_converged_step66.json')
    if old['state'] != 'running' or new['state'] != 'done' or new['n_steps'] != 66:
        raise ValueError('Provisional/completed Rx distinction lost')
    historical = ROOT / 'supporting/11_Molecular_dynamics/trajectories_not_included/trajectories.csv'
    with historical.open() as stream: old_rows = list(csv.DictReader(stream))
    with (ROOT / 'metadata/external_md_inventory.csv').open() as stream:
        external = {r['relative_path']: r for r in csv.DictReader(stream)}
    verified = 0
    for row in old_rows:
        # The exact historical header is asserted rather than guessed.
        name = row['file']
        rel = 'runs/' + row['run'] + '/' + name
        actual = external[rel]
        if actual['sha256'] != row['sha256'] or int(actual['size_bytes']) != int(row['bytes']):
            raise ValueError('External trajectory differs from historical manifest: ' + rel)
        verified += 1
    print(json.dumps({'passed': True, 'checks': ['S7 raw energies/geometries/convergence flags',
          'S8 stored geometries', 'S9 five-method energy differences',
          'S15 raw QE PBE energies and hybrid corrections with two formula units',
          'MD primary/sensitivity/validation roles and file references',
          'provisional versus completed Rx records'], 'external_trajectory_hashes_matching_historical_manifest': verified,
          'scientific_calculations_reproduced': False}, indent=2))


if __name__ == '__main__':
    main()
