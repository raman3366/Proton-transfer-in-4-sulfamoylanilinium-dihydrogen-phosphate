#!/usr/bin/env python3
"""Heuristic credential/privacy screen; output locations and categories, never secrets.

Review findings manually. Absence of matches is not a guarantee of absence.
"""
import argparse
import gzip
import io
import ipaddress
import json
from pathlib import Path
import re
import zipfile

from make_inventory import ROOT, publication_paths

PATTERNS = {
    'private_key': re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
    'github_token': re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,})'),
    'cloud_access_key': re.compile(rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
    'slack_token': re.compile(rb'xox[baprs]-[A-Za-z0-9-]{20,}'),
    'service_key': re.compile(rb'\bsk-(?:proj-|ant-)?[A-Za-z0-9_-]{35,}'),
    'credential_assignment': re.compile(rb'''(?i)(?:api[_-]?key|access[_-]?token|password|client[_-]?secret)\s*[:=]\s*["'][A-Za-z0-9_+/=-]{16,}["']'''),
    'credential_url': re.compile(rb'https?://[^\s/:]+:[^\s/@]{3,}@'),
}
IP = re.compile(rb'(?<![\d.])(?:\d{1,3}\.){3}\d{1,3}(?![\d.])')


def inspect(name, data, findings, containers=True):
    for category, pattern in PATTERNS.items():
        for match in pattern.finditer(data):
            findings.append({'path': name, 'line': data[:match.start()].count(b'\n') + 1, 'category': category})
    for match in IP.finditer(data):
        try: addr = ipaddress.ip_address(match[0].decode())
        except ValueError: continue
        if addr.is_loopback or addr.is_unspecified: continue
        # Package/version strings can also look like addresses; report for review.
        findings.append({'path': name, 'line': data[:match.start()].count(b'\n') + 1, 'category': 'address_like_text'})
    if containers and name.endswith('.gz'):
        inspect(name + ' [decompressed]', gzip.decompress(data), findings, False)
    elif containers and name.endswith(('.zip', '.npz', '.docx')):
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for member in archive.infolist():
                if not member.is_dir(): inspect(name + '::' + member.filename, archive.read(member), findings, False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    findings = []
    paths = publication_paths()
    for rel in paths:
        if rel == 'tools/scan_publication.py': continue
        inspect(rel, (ROOT / rel).read_bytes(), findings)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps({'files_screened': len(paths), 'findings': findings}, indent=2) + '\n')
    categories = {}
    for finding in findings:
        categories[finding['category']] = categories.get(finding['category'], 0) + 1
    print(json.dumps({'files_screened': len(paths), 'categories': categories, 'report': str(args.report)}, indent=2))


if __name__ == '__main__':
    main()
