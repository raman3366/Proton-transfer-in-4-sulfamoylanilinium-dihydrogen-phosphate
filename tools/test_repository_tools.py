#!/usr/bin/env python3
"""Regression checks for preservation and safe input rewriting; no scientific jobs."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import prepare_workspace as prepare
import validate_repository as validate


class PreservationTests(unittest.TestCase):
    def test_corruption_of_extracted_member_is_detected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder).resolve()
            archive_path = validate.ROOT / '17_SAPH_Supplementary_structures.zip'
            (root / archive_path.name).write_bytes(archive_path.read_bytes())
            with zipfile.ZipFile(archive_path) as archive:
                for member in archive.infolist():
                    if member.is_dir(): continue
                    path = root / member.filename.split('/', 1)[1]
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(archive.read(member))
            with patch.object(validate, 'ROOT', root):
                self.assertEqual(validate.archive_preservation(), 749)
                path = root / 'molecular/sfa_pair_pcm_d3.xyz'
                path.write_bytes(path.read_bytes() + b'\n')
                with self.assertRaisesRegex(ValueError, 'Original file changed'):
                    validate.archive_preservation()

    def test_qe_rewrite_changes_only_two_paths_and_rejects_bad_hash(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder).resolve()
            source = root / 'inputs/quantum_espresso/relax_i.in'
            source.parent.mkdir(parents=True)
            original = (prepare.ROOT / 'inputs/quantum_espresso/relax_i.in').read_text()
            source.write_text(original)
            metadata = json.loads((prepare.ROOT / 'metadata/pseudopotentials.json').read_text())
            pseudo = root / 'potentials'; pseudo.mkdir()
            # Non-scientific placeholder bytes exist only inside this temporary
            # test fixture; the production metadata and inputs are never altered.
            for record in metadata['files']:
                data = ('TEST FIXTURE ' + record['filename']).encode()
                (pseudo / record['filename']).write_bytes(data)
                record['recorded_md5'] = hashlib.md5(data).hexdigest()
            (root / 'metadata').mkdir()
            (root / 'metadata/pseudopotentials.json').write_text(json.dumps(metadata))
            with patch.object(prepare, 'ROOT', root):
                dest = root / 'prepared'; scratch = root / 'scratch'
                prepare.qe(dest, source, pseudo, scratch)
                result = (dest / source.name).read_text()
                normalized = result.replace(str(pseudo), '/home/ubuntu/cqc-sfa/tools/qe/pseudo').replace(str(scratch), '/mnt/scratch/periodic_qe/relax_i')
                self.assertEqual(original, normalized)
                (pseudo / 'H.pbe-rrkjus_psl.1.0.0.UPF').write_bytes(b'wrong')
                with self.assertRaisesRegex(ValueError, 'MD5 differs'):
                    prepare.qe(root / 'rejected', source, pseudo, scratch)
                self.assertFalse((root / 'rejected').exists())
            self.assertEqual(source.read_text(), original)

    def test_truncated_structure_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Truncated'):
            validate.xyz('2\nfixture\nH 0 0 0\n')


if __name__ == '__main__':
    unittest.main()
