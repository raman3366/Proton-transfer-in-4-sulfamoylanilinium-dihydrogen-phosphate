#!/usr/bin/env python3
"""Inventory every non-ignored publication file except the inventory itself."""
import csv
import hashlib
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = 'metadata/inventory.csv'


def publication_paths():
    output = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=ROOT)
    return sorted({p.decode('utf-8') for p in output.split(b'\0') if p and p.decode('utf-8') != INVENTORY})


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    rows = []
    for rel in publication_paths():
        path = ROOT / rel
        if path.is_symlink() or not path.is_file():
            raise ValueError('Publication path must be a regular file: ' + rel)
        rows.append([rel, path.stat().st_size, sha256(path)])
    dest = ROOT / INVENTORY
    dest.parent.mkdir(exist_ok=True)
    with dest.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['relative_path', 'size_bytes', 'sha256'])
        writer.writerows(rows)
    print(f'Wrote {len(rows)} file records to {INVENTORY}')


if __name__ == '__main__':
    main()
