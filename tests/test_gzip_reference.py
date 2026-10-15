"""Reference-reader regressions using the compiled legacy CLI."""

import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


BINARY = Path(os.environ.get("DRAGEN_OS_BINARY", "build/release/dragen-os")).resolve()


class GzipReferenceTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        rng = random.Random(13)
        self.fasta = b">chr1\n" + b"".join(
            ("".join(rng.choices("ACGT", k=80)) + "\n").encode()
            for _ in range(1500)
        )
        self.compressed = gzip.compress(self.fasta, mtime=0)
        self.plain = self.root / "reference.fasta"
        self.plain.write_bytes(self.fasta)
        self.gzipped = self.root / "reference.fasta.gz"
        self.gzipped.write_bytes(self.compressed)
        self.decoy = self.root / "decoy.fasta"
        self.decoy.write_bytes(b">decoy\n" + b"ACGT" * 100 + b"\n")

    def read_reference(self, reference, output_name, decoy=None):
        output = self.root / output_name
        output.mkdir()
        command = [
            str(BINARY),
            "--build-hash-table", "true",
            "--ht-reference", str(reference),
            "--output-directory", str(output),
            "--output-file-prefix", "test",
            "--ht-num-threads", "1",
            "--ht-size", "1MB",
            "--ht-mem-limit", "64MB",
            "--ht-test-only", "true",
        ]
        if decoy is not None:
            command.extend(["--ht-decoys", str(decoy)])
        result = subprocess.run(command, capture_output=True, text=True, timeout=120)
        return result, output

    def test_gzip_matches_plain_reference_with_and_without_decoys(self):
        for with_decoy in (False, True):
            with self.subTest(with_decoy=with_decoy):
                decoy = self.decoy if with_decoy else None
                plain_result, plain_output = self.read_reference(
                    self.plain, f"plain-{with_decoy}", decoy
                )
                gzip_result, gzip_output = self.read_reference(
                    self.gzipped, f"gzip-{with_decoy}", decoy
                )
                for result in (plain_result, gzip_result):
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    count = "2 sequences" if with_decoy else "1 sequence"
                    self.assertIn(f"Total: {count}, 120000 bases", result.stdout)
                plain_reference = (plain_output / "reference.bin").read_bytes()
                self.assertTrue(plain_reference)
                self.assertEqual(
                    plain_reference, (gzip_output / "reference.bin").read_bytes()
                )

    def test_rejects_incomplete_or_corrupt_gzip_before_writing_reference(self):
        corrupt_crc = bytearray(self.compressed)
        corrupt_crc[-8] ^= 1
        cases = {
            "truncated-data": self.compressed[: len(self.compressed) // 2],
            "missing-trailer": self.compressed[:-8],
            "invalid-crc": bytes(corrupt_crc),
        }
        for name, contents in cases.items():
            with self.subTest(case=name):
                reference = self.root / f"{name}.fasta.gz"
                reference.write_bytes(contents)
                result, output = self.read_reference(reference, name)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("compressed reference input", result.stderr)
                self.assertFalse((output / "reference.bin").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
