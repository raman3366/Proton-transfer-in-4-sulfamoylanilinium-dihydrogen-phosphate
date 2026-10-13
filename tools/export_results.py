#!/usr/bin/env python3
"""Export source-referenced numerical tables from stored calculation records.

No calculations are rerun. JSON numbers are read as Decimal to preserve the
recorded decimal representation. Blank cells mean absent/not applicable.
"""
import csv
from decimal import Decimal
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DFT = Path('supporting/06_DFT_results')
OUT = ROOT / 'data/derived'


def read(rel):
    return json.loads((ROOT / rel).read_text(), parse_float=Decimal)


def write(name, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(f'{name}: {len(rows)} rows')


def main():
    rows = []
    for path in sorted((ROOT / DFT / 'proton_scan').glob('[AB]_*.json')):
        rel = path.relative_to(ROOT); d = read(rel)
        if 'E_Ha' not in d: continue
        m = d['metrics']
        rows.append(dict(chain=d['chain'], point=d['index'], kind=d['kind'],
                         target_d_A=d.get('target_d_A'), E_Ha=d['E_Ha'],
                         d_N1_H1A_A=m['d_N1_H1A_A'], d_H1A_O5_A=m['d_H1A_O5_A'],
                         delta_A=m['delta_A'], d_N1_O5_A=m['d_N1_O5_A'],
                         converged=d['converged'], steps=d['steps'], accepted_on=d.get('accepted_on'),
                         source_json=rel.as_posix(), source_key='E_Ha; metrics; converged; steps', si_table='S7'))
    write('proton_scan.csv', rows)
    rel = DFT / 'eps_ladder/eps_ladder.json'; d = read(rel); rows = []
    for i, row in enumerate(d['rows']):
        m = row['metrics']
        rows.append(dict(epsilon=row['eps'], start=row['started_from'], E_Ha=row['E_Ha'],
                         d_N1_H1A_A=m['d_N1_H1A_A'], d_H1A_O5_A=m['d_H1A_O5_A'], delta_A=m['delta_A'],
                         d_N1_O5_A=m['d_N1_O5_A'], d_P_O5_A=m['d_P_O5_A'], proton_on=m['proton_on'],
                         converged=row['converged'], accepted_on=row.get('accepted_on'), steps=row['steps'],
                         source_json=rel.as_posix(), source_key=f'rows[{i}]', si_table='S8'))
    write('dielectric_ladder.csv', rows)
    rel = DFT / 'functional_check/fsp_scan.json'; d = read(rel); rows = []
    for method, geometries in d['energies'].items():
        for geometry, vals in geometries.items():
            rows.append(dict(method=method, geometry=geometry, E_Ha=vals['E_Ha'],
                             relative_to_salt_kcal_mol=d['relative'][method]['vs_first_kcal'][geometry],
                             source_json=rel.as_posix(), source_key=f'energies.{method}.{geometry}; relative.{method}.vs_first_kcal.{geometry}', si_table='S9'))
    write('method_comparison.csv', rows)
    rel = DFT / 'periodic_QE/profile.json'; d = read(rel); rows = []
    for i, vals in enumerate(d['profile']):
        row = dict(vals)
        row.update(source_json=rel.as_posix(), source_key=f'profile[{i}]', si_table='S15')
        rows.append(row)
    write('periodic_profile.csv', rows)
    rows = []
    for path in sorted((ROOT / DFT / 'periodic_QE/runs').glob('*/*.out')):
        lines = path.read_text().splitlines(); energies = []
        for i, line in enumerate(lines, 1):
            m = re.match(r'\s*!\s+total energy\s*=\s*([-+\d.Ee]+)\s+Ry', line)
            if m: energies.append((i, m[1]))
        if not energies: raise ValueError('No final total-energy record: ' + str(path))
        line, energy = energies[-1]
        rows.append(dict(run=path.stem, final_printed_total_energy_Ry=energy,
                         total_energy_records=len(energies), job_done_marker=any('JOB DONE.' in x for x in lines),
                         convergence_interpretation='See periodic_qe.json/profile.json; job completion alone is not ionic convergence',
                         source_file=path.relative_to(ROOT).as_posix(), source_line=line))
    write('qe_output_energies.csv', rows)
    rows = []
    for filename in ['qm_X.json', 'qm_Rx_converged_step66.json']:
        rel = Path('molecular_dynamics/quantum_chemistry/gate_G-R/relax') / filename
        d = read(rel); f = d['final']
        rows.append(dict(start=d['start'], state=d['state'], steps=d['n_steps'],
                         E_SCF_Ha=f['e_scf'], E_D3_Ha=f['e_d3'], E_total_Ha=f['e_total'],
                         signed_volume_N2_dimensionless=f['signed_volume'], Zn_N2_A=f['Zn_N2'],
                         source_json=rel.as_posix(), source_key='final', si_table='S40'))
    write('converged_zinc_model_relaxations.csv', rows)
    health = read(Path('molecular_dynamics/results/run_health.json'))
    conversions = {r['label']: r for r in read(Path('molecular_dynamics/results/conversion_times.json'))['runs']}
    campaign_path = Path('molecular_dynamics/results/campaign/first_program/campaign_table.csv')
    with (ROOT / campaign_path).open() as stream:
        campaign = {r['run']: r for r in csv.DictReader(stream)}
    rows = []
    for record in health['runs']:
        run = record['run']; model, name = run.split('/')
        role = 'validation' if name == 'X0' else 'primary' if model == 'bonded_2e' else 'sensitivity'
        c = conversions.get(run, {}).get('conversion') or {}
        row = dict(run=run, role=role, start=name[0], duration_ns=5 if role == 'validation' else 100,
                   verdict=campaign.get(run, {}).get('verdict', 'see validation gate'),
                   time_of_first_long_opposite_visit_from_release_ps=c.get('time_from_release_ps'),
                   visit_stage=c.get('stage'),
                   note='MIXED sensitivity visit is not a conversion' if model == 'bonded_2e_n2fit' else '',
                   series_path=f'molecular_dynamics/results/series/{model}_{name}.csv.gz' if role != 'validation' else '',
                   final_frame_path=f'molecular_dynamics/results/final_frames/{model}_{name}_final_solute.pdb.gz' if role != 'validation' else '',
                   source_json='molecular_dynamics/results/run_health.json; molecular_dynamics/results/conversion_times.json',
                   source_key=f'runs[run={run}]; runs[label={run}]', si_table='S40' if role == 'validation' else 'S39')
        rows.append(row)
    write('md_run_index.csv', rows)


if __name__ == '__main__':
    main()
