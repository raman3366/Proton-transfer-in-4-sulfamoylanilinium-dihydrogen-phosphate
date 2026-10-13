#!/usr/bin/env python3
"""Safely import selected research records from the verified 12 Oct handover ZIP."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path, PurePosixPath
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SHA256 = 'bc474c9a09bbf78e9f7e4768445e75035deae9e10f5f4125963fe8c22417c29b'
PREFIX = 'SAPH_deliverables_2026-10-12/'


def selected(rel):
    return (rel.startswith(('05_Tables/', '06_DFT_results/', '07_Docking/',
                            '08_ADMET_and_Fukui/', '11_Molecular_dynamics/'))
            or rel.startswith(('10_Reproducibility/inputs/',
                               '10_Reproducibility/scripts/dft/',
                               '10_Reproducibility/scripts/dock/',
                               '10_Reproducibility/scripts/paper/',
                               '10_Reproducibility/scripts/structures/'))
            or rel in {'10_Reproducibility/README.md',
                       '10_Reproducibility/SOFTWARE_VERSIONS.md',
                       '10_Reproducibility/source_map.csv',
                       '10_Reproducibility/excluded_files.csv',
                       '10_Reproducibility/git_commit.txt',
                       '10_Reproducibility/vm_archive_MANIFEST.sha256'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    h = hashlib.sha256()
    with args.archive.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    if h.hexdigest() != EXPECTED_SHA256:
        raise ValueError('Handover ZIP SHA-256 differs from the source on the Pi')
    if (ROOT / 'supporting').exists():
        raise FileExistsError('supporting/ already exists; refusing to overwrite')
    mappings = []
    with zipfile.ZipFile(args.archive) as archive:
        if archive.testzip() is not None:
            raise ValueError('ZIP CRC failure')
        manifest = {}
        for line in archive.read(PREFIX + 'MANIFEST.sha256').decode().splitlines():
            digest, name = line.split(None, 1)
            manifest[name.lstrip('*')] = digest
        planned = []
        names = set()
        for member in archive.infolist():
            p = PurePosixPath(member.filename)
            if p.is_absolute() or '..' in p.parts or stat.S_ISLNK(member.external_attr >> 16):
                raise ValueError('Unsafe ZIP member')
            if not member.filename.startswith(PREFIX):
                raise ValueError('Unexpected archive root')
            if member.is_dir():
                continue
            rel = member.filename[len(PREFIX):]
            if not selected(rel):
                continue
            # No full trajectories enter Git. Large wavefunctions are losslessly
            # compressed; their source hashes remain verifiable after decompression.
            if p.suffix.lower() in {'.dcd', '.xtc', '.trr', '.nc', '.ncdf'}:
                raise ValueError('Unexpected trajectory in handover')
            target = 'supporting/' + rel
            compress = member.file_size > 40 * 1024 * 1024
            if compress:
                target += '.gz'
            if target in names:
                raise ValueError('Duplicate destination')
            names.add(target)
            planned.append((member, rel, target, compress))
        for member, rel, target, compress in planned:
            data = archive.read(member)
            digest = hashlib.sha256(data).hexdigest()
            if digest != manifest[rel]:
                raise ValueError('Source manifest mismatch: ' + rel)
            dest = ROOT / target
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(gzip.compress(data, mtime=0) if compress else data)
            mappings.append([target, member.filename, len(data), digest,
                             'gzip_lossless' if compress else 'byte_copy'])
    with (ROOT / 'metadata/handover_import.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['repository_path', 'archive_member', 'source_size_bytes', 'source_sha256', 'transformation'])
        writer.writerows(mappings)
    print(f'Imported {len(mappings)} records; every record matched the handover manifest.')


if __name__ == '__main__':
    main()
