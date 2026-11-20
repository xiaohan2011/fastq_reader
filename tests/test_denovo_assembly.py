import gzip
import json
import os
from pathlib import Path
import shutil
import struct
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from Bio import SeqIO

from denovo_assembly import read_qc, run_assembly


class AssemblyValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reads = self.root / "sample.fastq"
        self.reads.write_text("@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n")

    def test_fastq_qc_plain_and_gzip(self):
        compressed = self.root / "sample.fastq.gz"
        with gzip.open(compressed, "wt") as handle:
            handle.write(self.reads.read_text())
        self.assertEqual(read_qc(self.reads), read_qc(compressed))
        qc = read_qc(self.reads)
        self.assertEqual(qc["read_count"], 2)
        self.assertEqual(qc["total_bases"], 6)
        self.assertEqual(qc["read_n50_bp"], 4)
        self.assertEqual(qc["mean_phred_score"], 40)

    def test_empty_truncated_and_invalid_fastq_fail(self):
        for text in ("", "@r\nACGT\n+\n", "@r\nACGT\n+\nII\n", "r\nAC\n+\nII\n"):
            with self.subTest(text=text):
                self.reads.write_text(text)
                with self.assertRaises(ValueError):
                    read_qc(self.reads)

    def test_inputs_and_existing_files_are_protected(self):
        with self.assertRaisesRegex(ValueError, "input files"):
            run_assembly(self.reads, self.root)
        output = self.root / "existing"
        output.mkdir()
        keep = output / "keep.txt"
        keep.write_text("preserve")
        with patch("denovo_assembly.shutil.which", return_value="/tool"), \
             patch("denovo_assembly.subprocess.check_output", return_value="version"), \
             patch("denovo_assembly.run_command", side_effect=RuntimeError("stop after preparing output")):
            with self.assertRaisesRegex(RuntimeError, "stop after preparing output"):
                run_assembly(self.reads, output)
        siblings = list(self.root.glob("existing_*"))
        self.assertEqual(len(siblings), 1)
        self.assertRegex(siblings[0].name, r"^existing_\d{8}_\d{6}_\d{6}$")
        self.assertTrue((siblings[0] / "run_state.json").is_file())
        self.assertEqual(keep.read_text(), "preserve")
        for kwargs in ({"threads": 0}, {"min_overlap": 100}, {"read_type": "unknown"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                run_assembly(self.reads, self.root / "new", **kwargs)

    def test_high_coverage_empty_assembly_retries_without_manual_settings(self):
        output = self.root / "new"
        attempts = []

        def flye(command, log, commands):
            attempts.append(command)
            destination = Path(command[command.index("--out-dir") + 1])
            if len(attempts) == 1:
                destination.mkdir(exist_ok=True)
                (destination / "flye.log").write_text(
                    "Overlap-based coverage: 500\nNo disjointigs were assembled\n"
                )
                raise RuntimeError("Flye produced no disjointigs")
            raise RuntimeError("stop after verifying retry")

        with patch("denovo_assembly.shutil.which", return_value="/tool"), \
             patch("denovo_assembly.subprocess.check_output", return_value="version"), \
             patch("denovo_assembly.read_qc", return_value={"total_bases": 50_000_000}), \
             patch("denovo_assembly.run_command", side_effect=flye):
            with self.assertRaisesRegex(RuntimeError, "stop after verifying retry"):
                run_assembly(self.reads, output)
        initial, retry = attempts
        self.assertNotIn("--genome-size", initial)
        self.assertNotIn("--asm-coverage", initial)
        self.assertEqual(retry[retry.index("--genome-size") + 1], "100000")
        self.assertEqual(retry[retry.index("--asm-coverage") + 1], "100")
        self.assertTrue((output / "flye" / "flye.log").is_file())

    def test_missing_tools_are_actionable(self):
        with patch("denovo_assembly.shutil.which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "flye is missing"):
                run_assembly(self.reads, self.root / "new")

    def test_resume_requires_same_input_fingerprint(self):
        output = self.root / "existing"
        output.mkdir()
        (output / "run_state.json").write_text(json.dumps({"signature": {}, "status": "failed"}))
        with self.assertRaisesRegex(ValueError, "Inputs/settings changed"):
            run_assembly(self.reads, output, resume=True)


DATA = Path(os.environ.get("FLYE_TEST_DATA", "/nonexistent"))
TOOLS = all(shutil.which(tool) for tool in ("flye", "minimap2", "samtools"))


@unittest.skipUnless(TOOLS and (DATA / "ecoli_500kb_reads_hifi.fastq.gz").is_file(),
                     "Set FLYE_TEST_DATA to Flye's bundled tests/data and install assembly tools")
class RealAssemblyTests(unittest.TestCase):
    def test_flye_snapgene_comparison_coverage_igv_and_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (DATA / "ecoli_500kb.fasta").open() as handle:
                record = next(SeqIO.parse(handle, "fasta"))
            def packet(kind, payload):
                return struct.pack(">BI", kind, len(payload)) + payload
            reference = root / "optional comparison.dna"
            reference.write_bytes(packet(9, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1)) +
                                  packet(0, b"\x00" + str(record.seq).encode("ascii")))
            reads = DATA / "ecoli_500kb_reads_hifi.fastq.gz"
            output = root / "assembly"
            settings = dict(read_type="pacbio-corr", threads=4, genome_size="500k",
                            min_overlap=1000, polish_iterations=1)
            summary = run_assembly(reads, output, reference, **settings)
            fasta = Path(summary["assembly_fasta"])
            self.assertEqual(fasta.read_bytes(), (output / "flye" / "assembly.fasta").read_bytes())
            self.assertGreater(sum(c["length_bp"] for c in summary["contigs"]), 400000)
            self.assertTrue(all(c["read_back_mean_depth"] > 0 for c in summary["contigs"]))
            self.assertGreater(summary["reference_alignment_blocks"], 0)
            session = ET.parse(summary["igv_session"]).getroot()
            self.assertEqual(session.get("genome"), fasta.name)
            self.assertTrue((output / session.find("Resources/Resource").get("path")).is_file())
            state = json.loads((output / "run_state.json").read_text())
            flye_command = next(c for c in state["commands"] if c[0] == "flye")
            self.assertNotIn(str(reference), flye_command)
            self.assertEqual(state["status"], "complete")
            before = fasta.stat().st_mtime_ns
            self.assertEqual(run_assembly(reads, output, reference, resume=True, **settings), summary)
            self.assertEqual(fasta.stat().st_mtime_ns, before)


if __name__ == "__main__":
    unittest.main()
