#!/usr/bin/env python3
"""Stream SHA-256 inventory of a separate data-deposit directory; never copy it."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys


def role(path):
    parts = path.parts
    if len(parts) > 2 and parts[0] == 'runs':
        model, run = parts[1:3]
        if run == 'X0':
            return 'validation'
        if model == 'bonded_2e' and run in {s + str(n) for s in 'XNR' for n in (1, 2, 3)}:
            return 'primary'
        if (model, run) in {('bonded_2s', 'R1'), ('bonded_2e_n2fit', 'R1')}:
            return 'sensitivity'
    return 'other_or_support'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--subdirectories', nargs='+', default=['runs', 'systems'])
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['relative_path', 'size_bytes', 'sha256', 'role', 'payload_type'])
    count = 0
    for name in args.subdirectories:
        base = (root / name).resolve(strict=True)
        if not base.is_relative_to(root):
            raise ValueError('Subdirectory escapes inventory root')
        for path in sorted(base.rglob('*')):
            if path.is_symlink() or not path.is_file():
                continue
            before = path.stat()
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
                    digest.update(block)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise RuntimeError('File changed during hashing: ' + str(path.relative_to(root)))
            rel = path.relative_to(root)
            payload = 'trajectory' if path.suffix.lower() in {'.dcd', '.xtc', '.trr', '.nc', '.ncdf'} else 'support'
            writer.writerow([rel.as_posix(), after.st_size, digest.hexdigest(), role(rel), payload])
            count += 1
    print(f'Inventoried {count} external files; no data copied.', file=sys.stderr)


if __name__ == '__main__':
    main()
