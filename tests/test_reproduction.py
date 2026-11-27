"""Regression checks for the verifier entry point and evidence reporting."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import reproduce
from research.single_target_recovery import check_single_target


class ReproductionTests(unittest.TestCase):
    def test_optimized_verifier_rejects_before_printing_evidence(self):
        result = subprocess.run(
            [sys.executable, "-O", "-B", "-m", "research.single_target_recovery.check_single_target"],
            cwd=reproduce.ROOT, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertTrue(result.stderr.rstrip().endswith("RuntimeError: This verifier requires Python assertions."))

    def test_exceptional_target_substitution_preserves_count_but_fails(self):
        path = reproduce.ROOT / "research" / "sparse_target_runs" / "witnesses.json"
        case = json.loads(path.read_text(encoding="utf-8"))[0]
        self.assertEqual(case["exceptional_logs"], [0, 17, 120, 137])
        case["exceptional_logs"][1] = 18
        with TemporaryDirectory() as tmp:
            witness_file = Path(tmp) / "witnesses.json"
            witness_file.write_text(json.dumps([case]), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "exceptional targets differ"):
                check_single_target.check_sparse_witnesses(witness_file)

    def test_unreadable_evidence_is_reported_and_later_checks_run(self):
        names = ("different", "missing", "bad_utf8", "directory", "last")
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            scripts = root / "checks"
            evidence = root / "evidence"
            scripts.mkdir()
            evidence.mkdir()
            for name in names:
                (scripts / f"{name}.py").write_text("print('observed')", encoding="utf-8")
            (evidence / "different.txt").write_text("other", encoding="utf-8")
            (evidence / "bad_utf8.txt").write_bytes(bytes([0xff]))
            (evidence / "directory.txt").mkdir()
            (evidence / "last.txt").write_text("observed\n", encoding="utf-8")

            with (patch.object(reproduce, "ROOT", root), patch.object(reproduce, "CHECKS", names),
                  patch.object(sys, "argv", ["reproduce.py"]), contextlib.redirect_stdout(io.StringIO()) as out,
                  self.assertRaises(SystemExit) as stopped):
                reproduce.main()

        output = out.getvalue()
        for name in ("missing", "bad_utf8", "directory"):
            self.assertIn(f"{name}: FAILED to read {name}.txt", output)
        self.assertIn("last: matches recorded evidence.", output)
        self.assertEqual(str(stopped.exception), "4 of 5 checks failed or differed: different, missing, bad_utf8, directory")
