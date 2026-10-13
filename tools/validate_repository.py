#!/usr/bin/env python3
"""Validate preservation, checksums, file structure and selected references.

This is a packaging/consistency check, not a reproduction of DFT, docking or MD.
The core checks use only the Python standard library.
"""
import argparse
import ast
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
import zipfile

from make_inventory import ROOT, INVENTORY, publication_paths, sha256


def require(condition, message):
    if not condition:
        raise ValueError(message)


def numeric(values):
    for v in values:
        require(math.isfinite(float(v)), 'Nonfinite coordinate')


def xyz(text):
    lines = text.splitlines()
    i = frames = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i]); require(n > 0, 'Invalid XYZ atom count')
        require(i + n + 2 <= len(lines), 'Truncated XYZ frame')
        for line in lines[i + 2:i + n + 2]:
            fields = line.split()
            require(len(fields) >= 4, 'Short XYZ coordinate record')
            numeric(fields[1:4])
        i += n + 2; frames += 1
    require(frames, 'No XYZ frames')


def pdb(text):
    count = 0
    for line in text.splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            numeric([line[30:38], line[38:46], line[46:54]])
            count += 1
    require(count, 'No PDB coordinate records')


def sdf(text):
    require(text.rstrip().endswith('$$$$'), 'Missing SDF record terminator')
    for block in text.split('$$$$')[:-1]:
        if block.startswith('\n'):
            block = block[1:]
        lines = block.splitlines()
        require(len(lines) >= 4, 'Short SDF record')
        counts = lines[3]
        require('V2000' in counts, 'Only V2000 SDF checked by this parser')
        n, nb = int(counts[:3]), int(counts[3:6])
        require(len(lines) >= 4 + n + nb, 'Truncated SDF')
        for line in lines[4:4 + n]:
            numeric([line[:10], line[10:20], line[20:30]])
        for line in lines[4 + n:4 + n + nb]:
            require(1 <= int(line[:3]) <= n and 1 <= int(line[3:6]) <= n, 'SDF bond index out of range')
        require('M  END' in lines, 'Missing SDF M END')


def check_format(path, data, counters, extended):
    suffix = path.suffix.lower()
    if suffix == '.gz':
        data = gzip.decompress(data)
        counters['gzip_crc'] += 1
        path = path.with_suffix(''); suffix = path.suffix.lower()
    if suffix in {'.json', '.jsonl', '.csv', '.xyz', '.extxyz', '.pdb', '.pdbqt', '.sdf', '.py', '.cif', '.yml', '.cff'}:
        text = data.decode('utf-8')
        if suffix == '.json':
            constants = []
            json.loads(text, parse_constant=lambda value: constants.append(value))
            if constants:
                counters['json_with_nonstandard_constants'] += 1
        elif suffix == '.jsonl':
            for line in text.splitlines():
                if line.strip(): json.loads(line)
        elif suffix == '.csv':
            rows = list(csv.reader(io.StringIO(text), strict=True))
            require(rows, 'Empty CSV')
            # CSV is a record format and need not be a rectangular table. Ragged
            # tables are counted explicitly; critical tables have separate checks.
            if len({len(r) for r in rows if r}) > 1:
                counters['csv_variable_width'] += 1
        elif suffix in {'.xyz', '.extxyz'}:
            if path.name.endswith('.maps.xyz'):
                rows = [line.split() for line in text.splitlines() if line.strip()]
                require(len(rows) == 3 and all(len(row) == 2 for row in rows), 'AutoGrid extent format')
                for row in rows:
                    numeric(row)
                    require(float(row[0]) < float(row[1]), 'AutoGrid extent ordering')
                counters['autogrid_extents'] += 1
            else:
                xyz(text)
        elif suffix in {'.pdb', '.pdbqt'}:
            if path.name == 'cofactors.pdbqt' and not text:
                counters['empty_cofactor_placeholders'] += 1
            else:
                pdb(text)
        elif suffix == '.sdf':
            sdf(text)
        elif suffix == '.py':
            ast.parse(text)
        elif suffix == '.cif':
            require(re.search(r'(?m)^data_', text), 'No CIF data block')
            if extended:
                import gemmi
                gemmi.cif.read_string(text)
                counters['cif_gemmi'] += 1
        elif suffix in {'.yml', '.cff'} and extended:
            import yaml
            yaml.safe_load(text)
            counters['yaml_parsed'] += 1
        counters[suffix[1:]] += 1
    elif suffix in {'.zip', '.npz'}:
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            require(archive.testzip() is None, 'ZIP/NPZ CRC failure')
        counters['zip_npz_crc'] += 1


def archive_preservation():
    with zipfile.ZipFile(ROOT / '17_SAPH_Supplementary_structures.zip') as archive:
        require(archive.testzip() is None, 'Original ZIP CRC failure')
        count = 0
        for entry in archive.infolist():
            if entry.is_dir(): continue
            path = PurePosixPath(entry.filename)
            require(not path.is_absolute() and '..' not in path.parts, 'Unsafe ZIP member')
            rel = Path(*path.parts[1:])
            require((ROOT / rel).read_bytes() == archive.read(entry), f'Original file changed: {rel}')
            count += 1
        require(count == 749, 'Unexpected original archive file count')
    return count


def imported_preservation():
    rows = list(csv.DictReader((ROOT / 'metadata/handover_import.csv').open()))
    for row in rows:
        data = (ROOT / row['repository_path']).read_bytes()
        if row['transformation'] == 'gzip_lossless': data = gzip.decompress(data)
        require(len(data) == int(row['source_size_bytes']), 'Imported size mismatch')
        require(hashlib.sha256(data).hexdigest() == row['source_sha256'], 'Imported source checksum mismatch: ' + row['repository_path'])
    return len(rows)


def references():
    for path in (ROOT / 'inputs/pyscf/case_definitions').glob('*.json'):
        case = json.loads(path.read_text())
        start = ROOT / 'inputs/pyscf/start_geometries' / Path(case['xyz']).name
        require(start.is_file(), f'Missing start geometry: {path.name}')
        n = int(start.read_text().splitlines()[0])
        if 'cif_labels' in case:
            require(n == len(case['cif_labels']), f'Label/atom count mismatch: {path.name}')
    for path in (ROOT / 'docking/receptors').glob('*/vina_config.txt'):
        receptor = re.search(r'(?m)^receptor\s*=\s*(.+)$', path.read_text())[1].strip()
        require((path.parent / receptor).is_file(), f'Missing docking receptor: {path}')
    expected = {f['filename'] for f in json.loads((ROOT / 'metadata/pseudopotentials.json').read_text())['files']}
    for path in (ROOT / 'inputs/quantum_espresso').glob('*.in'):
        named = re.findall(r'(?m)^\s*[A-Z][a-z]?\s+[\d.]+\s+(\S+\.(?:UPF|upf))\s*$', path.read_text())
        require(len(named) == 6 and set(named) <= expected, f'Unregistered QE potential: {path.name}')
    # Check newly authored Markdown links; archived historical paths are preserved
    # and explicitly documented as requiring reconstruction.
    for path in [ROOT / 'README.md', *(ROOT / 'docs').glob('*.md')]:
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            if '://' in target or target.startswith('#'): continue
            clean = target.split('#')[0]
            require((path.parent / clean).exists(), f'Broken documentation link: {path.name}: {target}')


def md_series():
    files = sorted((ROOT / 'molecular_dynamics/results/series').glob('*.csv.gz'))
    require(len(files) == 11, 'Expected eleven compact MD series')
    primary_rotated = 0
    for path in files:
        with gzip.open(path, 'rt') as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == 10130, f'MD series row count: {path.name}')
        prod = [r for r in rows if r['part'] == 'prod']
        require(len(prod) == 10000, 'Expected 10000 production frames')
        require(abs(float(rows[0]['time_ns']) + 1.29) < 1e-9 and float(rows[-1]['time_ns']) == 100, 'MD time endpoints')
        for a, b in zip(rows, rows[1:]):
            require(abs(float(b['time_ns']) - float(a['time_ns']) - .01) < 1e-8, 'MD frame interval')
        for row in rows:
            d1, d2 = float(row['d1_A']), float(row['d2_A'])
            state = 'C' if d1 < 6.7 and d2 > 5.2 else 'R' if d1 > 6.7 and d2 < 5.2 else 'I'
            require(row['state'] == state, f'MD state inconsistent with stored distances: {path.name}')
        if re.fullmatch(r'bonded_2e_[XNR][123]\.csv\.gz', path.name):
            primary_rotated += sum(r['state'] == 'R' for r in prod)
        frame = ROOT / 'molecular_dynamics/results/final_frames' / path.name.replace('.csv.gz', '_final_solute.pdb.gz')
        lines = gzip.decompress(frame.read_bytes()).decode().splitlines()
        atoms = {}
        for line in lines:
            if line.startswith(('ATOM  ', 'HETATM')):
                atoms[(int(line[22:26]), line[12:16].strip())] = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        require(len(atoms) == 4094, 'Final solute atom count')
        n1 = atoms[(261, 'N1')]
        for atom, column in [((129, 'CZ'), 'd1_A'), ((66, 'OD1'), 'd2_A')]:
            require(abs(math.dist(n1, atoms[atom]) - float(rows[-1][column])) < .002, 'Final frame/series distance mismatch')
    require(primary_rotated == 0, 'Primary production contains rotated frames')
    return {'series': len(files), 'rows_per_series': 10130, 'primary_production_rotated_frames': primary_rotated}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-inventory', action='store_true', help='Preparation only; not a release validation')
    parser.add_argument('--extended', action='store_true', help='Also parse CIF using gemmi and YAML using PyYAML')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    errors, counts = [], Counter()
    report = {'scientific_calculations_reproduced': False}
    def check(name, fn):
        try: report[name] = fn()
        except Exception as exc: errors.append({'check': name, 'error': str(exc)})
    check('original_files_byte_identical', archive_preservation)
    check('handover_files_source_verified', imported_preservation)
    paths = publication_paths()
    def inventory():
        rows = list(csv.DictReader((ROOT / INVENTORY).open()))
        require(len(rows) == len({r['relative_path'] for r in rows}), 'Duplicate inventory path')
        require(set(paths) == {r['relative_path'] for r in rows}, 'Inventory path coverage differs')
        for row in rows:
            path = ROOT / row['relative_path']
            require(path.stat().st_size == int(row['size_bytes']) and sha256(path) == row['sha256'], 'Inventory mismatch: ' + row['relative_path'])
        return len(rows)
    if not args.skip_inventory: check('inventory_files_verified', inventory)
    for rel in paths:
        path = ROOT / rel
        try:
            require(path.is_file() and not path.is_symlink(), 'Not a regular file')
            require(path.stat().st_size < 50 * 1024 * 1024, 'Publication file exceeds 50 MiB')
            require(path.suffix.lower() not in {'.dcd', '.xtc', '.trr', '.nc', '.ncdf'}, 'Trajectory in publication files')
            check_format(path, path.read_bytes(), counts, args.extended)
        except Exception as exc:
            errors.append({'check': 'file_structure', 'path': rel, 'error': str(exc)})
    check('important_references', references)
    check('md_series_consistency', md_series)
    report.update(file_checks=dict(sorted(counts.items())), errors=errors, passed=not errors,
                  inventory_checked=not args.skip_inventory)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
