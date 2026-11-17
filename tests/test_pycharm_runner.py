"""Exercise the plain-Python settings launcher without requiring a notebook IDE."""

from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_fastq_reader as launcher


class PyCharmRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reads = self.root / "reads"
        self.reads.mkdir()
        self.reference = self.root / "chosen reference.fa"
        rng = random.Random(42)
        sequence = "".join(rng.choices("ACGT", k=3000))
        self.selected = self.reads / "chosen sample.fastq"
        self.selected.write_text(f"@r\n{sequence}\n+\n" + "I" * len(sequence) + "\n")
        self.reference.write_text(">ref\n" + sequence + "\n")

    def test_pycharm_folder_run_produces_qc_without_jupyter(self):
        with patch.multiple(launcher, FASTQ_DIR=self.reads, FASTQ_FILE=None,
                            REFERENCE_FILE=None, REFERENCE_DIR=None, OUTPUT_DIR=None, QC_ONLY=True):
            report = launcher.run_from_settings()
        self.assertEqual(report.iloc[0]["status"], "qc_only")
        self.assertTrue((self.reads / "results" / "chosen sample_summary.csv").is_file())

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_single_file_pairing_outputs_bam_and_igv_session(self):
        (self.reads / "unselected.fastq").write_text("bad data must not be read")
        with patch.multiple(launcher, FASTQ_DIR=None, FASTQ_FILE=self.selected,
                            REFERENCE_FILE=self.reference, OUTPUT_DIR=None, QC_ONLY=False):
            report = launcher.run_from_settings()
        self.assertEqual(len(report), 1)
        self.assertEqual(report.iloc[0]["status"], "aligned")
        subprocess.run(["samtools", "quickcheck", report.iloc[0]["bam"]], check=True)
        self.assertTrue(Path(report.iloc[0]["igv_session"]).is_file())

    def test_main_reports_configuration_errors(self):
        with patch.multiple(launcher, FASTQ_DIR=None, FASTQ_FILE=None, OUTPUT_DIR=None):
            self.assertEqual(launcher.main(), 2)

    def test_import_does_not_start_analysis(self):
        result = subprocess.run([sys.executable, "-c", "import run_fastq_reader"],
                                cwd=Path(launcher.__file__).parent, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
