#!/usr/bin/env python3
"""Create an isolated calculation layout without running a calculation."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]


def molecular(dest):
    # Real copies are essential: common.py resolves __file__, so symlinks fail.
    shutil.copytree(ROOT / 'inputs/pyscf/scripts', dest / 'dft')
    base = dest / 'molecules/sfa'
    shutil.copytree(ROOT / 'inputs/pyscf/case_definitions', base / 'cases')
    shutil.copytree(ROOT / 'inputs/pyscf/start_geometries', base / 'input')
    for case in (base / 'cases').glob('*.json'):
        data = json.loads(case.read_text())
        target = (dest / data['xyz']).resolve()
        if not target.is_relative_to(dest) or not target.is_file():
            raise ValueError(f'Unresolved start geometry in {case.name}')
    print(f'Prepared {len(list((base / "cases").glob("*.json")))} cases in {dest}')
    print('No scientific program was imported or run. See docs/REPRODUCING.md.')


def project(dest):
    """Restore original relative paths using the handover's provenance map."""
    source_roots = [
        '/root/workspace/cqc-4-sulfamoylanilinium-dihydrogen-phosphate/',
        '/root/cqc-sfa-archive/vm/cqc-sfa/',
    ]
    imports = {}
    with (ROOT / 'metadata/handover_import.csv').open() as stream:
        for row in csv.DictReader(stream):
            package_path = row['archive_member'].split('/', 1)[1]
            imports[package_path] = row
    planned = {}
    skipped = 0
    with (ROOT / 'supporting/10_Reproducibility/source_map.csv').open() as stream:
        for row in csv.DictReader(stream):
            entry = imports.get(row['package_path'])
            if entry is None:
                continue
            prefix = next((r for r in source_roots if row['source'].startswith(r)), None)
            if prefix is None:
                skipped += 1; continue
            relative = Path(row['source'][len(prefix):])
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('Unsafe source-map destination')
            target = dest / relative
            if target in planned and planned[target]['source_sha256'] != entry['source_sha256']:
                raise ValueError('Conflicting source copies: ' + str(relative))
            planned[target] = entry
    # Preflight all destinations before making a directory or copying files.
    for target in planned:
        if not target.resolve().is_relative_to(dest):
            raise ValueError('Source-map path escapes destination')
    dest.mkdir(parents=True)
    for target, entry in sorted(planned.items()):
        data = (ROOT / entry['repository_path']).read_bytes()
        if entry['transformation'] == 'gzip_lossless': data = gzip.decompress(data)
        if hashlib.sha256(data).hexdigest() != entry['source_sha256']:
            raise ValueError('Source checksum mismatch: ' + entry['repository_path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    (dest / 'RESTORATION.json').write_text(json.dumps({
        'copied_paths': len(planned), 'unmapped_generated_or_external_files': skipped,
        'note': 'Only verified original relative paths restored. Absolute machine paths, external software, grids, trajectories and generated-only assets require separate configuration; no calculation was run.'
    }, indent=2) + '\n')
    print(f'Restored {len(planned)} original relative paths in {dest}; {skipped} generated/external records not mapped.')


def qe(dest, source, pseudo, scratch):
    original = source.read_text()
    if not source.resolve().is_relative_to(ROOT / 'inputs/quantum_espresso'):
        raise ValueError('Select an archived Quantum ESPRESSO input')
    if any(x in str(pseudo) + str(scratch) for x in ["'", '\n', '\r']):
        raise ValueError('Paths cannot contain quotes or line breaks')
    text = original
    for key, value in [('pseudo_dir', pseudo), ('outdir', scratch)]:
        text, count = re.subn(r"(?m)^(\s*" + key + r"\s*=\s*)'[^']*'", lambda m: m[1] + "'" + str(value) + "'", text)
        if count != 1:
            raise ValueError(f'Expected one {key} in {source.name}')
    species = original.split('ATOMIC_SPECIES', 1)[1].strip().splitlines()
    names = []
    for line in species:
        bits = line.split()
        if len(bits) != 3 or not bits[2].lower().endswith('.upf'):
            break
        names.append(bits[2])
    records = json.loads((ROOT / 'metadata/pseudopotentials.json').read_text())['files']
    expected = {r['filename']: r['recorded_md5'] for r in records}
    for name in names:
        data = (pseudo / name).read_bytes()
        if hashlib.md5(data).hexdigest() != expected[name]:
            raise ValueError('Pseudopotential MD5 differs from the archived record: ' + name)
    dest.mkdir(parents=True)
    (dest / source.name).write_text(text)
    print(f'Prepared {dest / source.name}; verified {len(names)} pseudopotentials. No calculation run.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('molecular')
    p.add_argument('--dest', type=Path, default=ROOT / 'work/molecular')
    p = sub.add_parser('project')
    p.add_argument('--dest', type=Path, default=ROOT / 'work/project')
    p = sub.add_parser('qe')
    p.add_argument('input', type=Path)
    p.add_argument('--dest', type=Path, default=ROOT / 'work/qe')
    p.add_argument('--pseudo-dir', required=True, type=Path)
    p.add_argument('--scratch-dir', required=True, type=Path)
    args = parser.parse_args()
    dest = args.dest.resolve()
    if dest.exists():
        parser.error('Destination exists; choose a new directory (no overwrites).')
    if args.command == 'molecular':
        molecular(dest)
    elif args.command == 'project':
        project(dest)
    else:
        qe(dest, args.input.resolve(), args.pseudo_dir.resolve(), args.scratch_dir.resolve())


if __name__ == '__main__':
    main()
