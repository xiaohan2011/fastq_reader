"""Functional checks using deterministic FASTQ, FASTA and SnapGene fixtures."""

import gzip
import json
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

import pandas as pd
from PIL import Image

from batch_fastq_reader import create_igv_sessions, main, run_batch, sample_id


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.reads = self.root / "reads"
        self.refs = self.root / "references"
        self.output = self.root / "output"
        self.reads.mkdir()
        self.refs.mkdir()
        rng = random.Random(42)
        self.sequence = "".join(rng.choices("ACGT", k=15000))

    def fastq(self, name, lengths=(1000, 12000, 5000)):
        path = self.reads / name
        text = "".join(f"@read{i}\n{self.sequence[:length]}\n+\n{'I' * length}\n"
                       for i, length in enumerate(lengths))
        if name.lower().endswith(".gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return path

    def fasta(self, name):
        path = self.refs / name
        path.write_text(">plasmid\n" + self.sequence + "\n")
        return path

    def snapgene(self, name):
        def packet(kind, data):
            return struct.pack(">BI", kind, len(data)) + data
        path = self.refs / name
        path.write_bytes(packet(9, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1)) +
                         packet(0, b"\x01" + self.sequence.encode("ascii")))
        return path

    def run_batch(self, **kwargs):
        return run_batch(self.reads, self.output, self.refs, threads=2, **kwargs)

    def test_missing_reference_gzip_qc_and_skip(self):
        self.fastq("1132.fastq.gz")
        first = self.run_batch()
        self.assertEqual(first.iloc[0]["status"], "qc_only_missing_reference")
        df = pd.read_csv(self.output / "1132_summary.csv")
        self.assertEqual(df["Length"].tolist(), [12000, 5000, 1000])
        self.assertEqual(df["Average Quality"].tolist(), [40, 40, 40])
        with Image.open(self.output / "1132_length_hist.png") as image:
            image.verify()
        before = (self.output / "1132_summary.csv").stat().st_mtime_ns
        second = self.run_batch()
        self.assertEqual(second.iloc[0]["action"], "skipped")
        self.assertEqual((self.output / "1132_summary.csv").stat().st_mtime_ns, before)
        self.assertTrue((self.output / "batch_report.csv").is_file())
        self.assertFalse(any(path.is_dir() for path in self.output.iterdir()))

    def test_matching_is_exact_and_search_is_top_level(self):
        self.fastq("1132.fastq")
        self.fasta("11320 other.fa")
        nested = self.refs / "nested"
        nested.mkdir()
        (nested / "1132.fa").write_text(">hidden\nACGT\n")
        nested_reads = self.reads / "nested"
        nested_reads.mkdir()
        (nested_reads / "999.fastq").write_text("bad")
        result = self.run_batch()
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["status"], "qc_only_missing_reference")
        self.assertEqual(sample_id(Path("1132 LPPC EE R-S.dna")), "1132")

    def test_ambiguous_references_get_qc_and_error_exit(self):
        self.fastq("1132.fastq")
        self.fasta("1132 A.fa")
        self.snapgene("1132 B.dna")
        result = self.run_batch()
        self.assertEqual(result.iloc[0]["status"], "ambiguous_reference")
        self.assertTrue((self.output / "1132_summary.csv").is_file())
        self.assertFalse((self.output / "1132_sorted.bam").exists())
        self.assertEqual(main([str(self.reads), "--reference-dir", str(self.refs),
                               "--output-dir", str(self.output), "--threads", "2"]), 1)

    def test_malformed_sample_does_not_stop_batch(self):
        (self.reads / "1132.fastq").write_text("@broken\nACGT\n+\nI\n")
        self.fastq("1133.FQ.GZ")
        result = self.run_batch()
        self.assertEqual(result["status"].tolist(), ["failed", "qc_only_missing_reference"])
        self.assertIn("length mismatch", result.iloc[0]["message"])
        self.assertTrue((self.output / "1133_summary.csv").is_file())

    def test_empty_fastq_and_incomplete_record_are_failures(self):
        (self.reads / "1.fastq").write_text("")
        (self.reads / "2.fastq").write_text("@r\nACGT\n+\n")
        result = self.run_batch()
        self.assertEqual(result["status"].tolist(), ["failed", "failed"])
        self.assertIn("No FASTQ records", result.iloc[0]["message"])
        self.assertIn("Incomplete", result.iloc[1]["message"])

    def test_changed_input_and_missing_artifact_require_overwrite(self):
        path = self.fastq("1132.fastq")
        self.run_batch()
        path.write_text(path.read_text() + "@new\nACGT\n+\nIIII\n")
        result = self.run_batch()
        self.assertEqual(result.iloc[0]["status"], "failed")
        self.assertIn("--overwrite", result.iloc[0]["message"])
        self.assertEqual(self.run_batch(overwrite=True).iloc[0]["reads"], 4)
        (self.output / "1132_length_hist.png").unlink()
        self.assertEqual(self.run_batch().iloc[0]["status"], "failed")
        self.assertEqual(self.run_batch(overwrite=True).iloc[0]["status"], "qc_only_missing_reference")

    def test_new_reference_requires_explicit_rebuild(self):
        self.fastq("1132.fastq")
        self.run_batch()
        self.fasta("1132.fa")
        self.assertEqual(self.run_batch().iloc[0]["status"], "failed")

    def test_duplicate_sample_id_and_same_stem_do_not_clobber(self):
        for name in ("1132.fastq", "1132.fastq.gz", "1132 extra.fq"):
            self.fastq(name)
        self.assertEqual(len(self.run_batch()), 3)
        for name in ("1132.fastq_summary.csv", "1132.fastq.gz_summary.csv", "1132 extra_summary.csv"):
            self.assertTrue((self.output / name).is_file())
        self.assertTrue(self.run_batch()["action"].eq("skipped").all())

    def test_qc_only_and_filename_without_id(self):
        self.fastq("sample.fastq")
        self.assertEqual(self.run_batch().iloc[0]["status"], "qc_only_no_sample_id")
        self.assertEqual(self.run_batch(qc_only=True, overwrite=True).iloc[0]["status"], "qc_only")

    def test_invalid_reference_preserves_qc(self):
        self.fastq("1132.fastq")
        (self.refs / "1132.dna").write_bytes(b"not a SnapGene file")
        result = self.run_batch()
        self.assertEqual(result.iloc[0]["status"], "failed")
        self.assertTrue((self.output / "1132_summary.csv").is_file())
        self.assertFalse((self.output / ".1132_complete.json").exists())

    def test_missing_tools_reports_failure_after_qc(self):
        self.fastq("1132.fastq")
        self.fasta("1132.fa")
        with patch("batch_fastq_reader.shutil.which", return_value=None):
            result = self.run_batch()
        self.assertIn("minimap2 is missing", result.iloc[0]["message"])
        self.assertTrue((self.output / "1132_summary.csv").is_file())

    @unittest.skipUnless(shutil.which("samtools"), "samtools not installed")
    def test_mapper_failure_keeps_diagnostics_and_continues(self):
        self.fastq("1132.fastq")
        self.fastq("1133.fastq")
        self.fasta("1132.fa")
        tools = self.root / "tools"
        tools.mkdir()
        mapper = tools / "minimap2"
        mapper.write_text("#!/bin/sh\necho 'deliberate mapper failure' >&2\nexit 23\n")
        mapper.chmod(0o755)
        with patch.dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"]):
            result = self.run_batch()
        self.assertEqual(result["status"].tolist(), ["failed", "qc_only_missing_reference"])
        self.assertIn("minimap2=23", result.iloc[0]["message"])
        self.assertIn("deliberate mapper failure", (self.output / "1132_alignment.log").read_text())
        self.assertTrue((self.output / "1132_summary.csv").is_file())
        self.assertFalse((self.output / "1132_sorted.bam").exists())
        self.assertFalse((self.output / ".1132_complete.json").exists())

    def test_invalid_folders_and_settings(self):
        with self.assertRaisesRegex(ValueError, "No supported FASTQ"):
            self.run_batch()
        self.fastq("1132.fastq")
        with self.assertRaisesRegex(ValueError, "separate output"):
            run_batch(self.reads, self.reads)
        with self.assertRaisesRegex(ValueError, "positive"):
            self.run_batch(bin_size=0)

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_fasta_snapgene_alignment_and_overwrite(self):
        for identifier in ("1132", "1133"):
            self.fastq(identifier + ".fastq.gz")
        self.fasta("1132 LPPC EE R-S.fa")
        self.snapgene("1133 LPPC EE R-T.dna")
        result = self.run_batch()
        self.assertEqual(result["status"].tolist(), ["aligned", "aligned"])
        for identifier in ("1132", "1133"):
            bam = self.output / f"{identifier}_sorted.bam"
            subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
            count = subprocess.check_output(["samtools", "view", "-c", "-F", "4", str(bam)])
            self.assertEqual(int(count), 3)
            self.assertTrue(Path(str(bam) + ".bai").is_file())
            session = ET.parse(self.output / f"{identifier}_igv.xml").getroot()
            self.assertEqual(session.tag, "Session")
            self.assertEqual(session.get("genome"), f"{identifier}_reference.fasta")
            self.assertEqual(session.find("Resources/Resource").get("path"), bam.name)
            self.assertTrue((self.output / f"{identifier}_reference.fasta.fai").is_file())
        converted = (self.output / "1133_reference.fasta").read_text()
        self.assertEqual("".join(converted.splitlines()[1:]), self.sequence)
        self.assertTrue(self.run_batch()["action"].eq("skipped").all())
        (self.output / "unrelated.txt").write_text("keep me")
        rerun = self.run_batch(overwrite=True, qc_only=True)
        self.assertTrue(rerun["status"].eq("qc_only").all())
        self.assertFalse((self.output / "1132_sorted.bam").exists())
        self.assertFalse((self.output / "1132_igv.xml").exists())
        self.assertFalse((self.output / "1132_reference.fasta.fai").exists())
        self.assertEqual((self.output / "unrelated.txt").read_text(), "keep me")

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_igv_sessions_upgrade_old_results_and_move_with_folder(self):
        self.fastq("1132 & sample.fastq")
        self.fasta("1132 reference.fa")
        self.run_batch()
        bam = self.output / "1132 & sample_sorted.bam"
        before = bam.stat().st_mtime_ns
        session = self.output / "1132 & sample_igv.xml"
        marker = self.output / ".1132 & sample_complete.json"
        old = json.loads(marker.read_text())
        old["artifacts"] = [name for name in old["artifacts"] if not name.endswith(("_igv.xml", ".fai"))]
        old["result"].pop("igv_session")
        marker.write_text(json.dumps(old))
        session.unlink()
        (self.output / "1132 & sample_reference.fasta.fai").unlink()
        self.assertEqual(self.run_batch().iloc[0]["action"], "skipped")
        self.assertEqual(bam.stat().st_mtime_ns, before)
        self.assertTrue(session.is_file())
        self.assertIn("igv_session", json.loads(marker.read_text())["result"])
        moved = self.root / "moved results"
        shutil.move(self.output, moved)
        root = ET.parse(moved / session.name).getroot()
        self.assertTrue((moved / root.get("genome")).is_file())
        self.assertTrue((moved / root.find("Resources/Resource").get("path")).is_file())
        (moved / session.name).unlink()
        report = create_igv_sessions(moved)
        self.assertEqual(report.iloc[0]["status"], "ready")
        self.assertEqual((moved / bam.name).stat().st_mtime_ns, before)
        saved = (moved / session.name).read_bytes()
        self.assertEqual(create_igv_sessions(moved).iloc[0]["status"], "ready")
        self.assertEqual((moved / session.name).read_bytes(), saved)
        self.assertTrue((moved / "igv_session_report.csv").is_file())

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_igv_only_continues_after_missing_index_and_cli(self):
        for identifier in ("1132", "1133"):
            self.fastq(identifier + ".fastq")
            self.fasta(identifier + ".fa")
        self.run_batch()
        for identifier in ("1132", "1133"):
            (self.output / f"{identifier}_igv.xml").unlink()
        (self.output / "1132_sorted.bam.bai").unlink()
        result = create_igv_sessions(self.output)
        self.assertEqual(result["status"].tolist(), ["failed", "ready"])
        self.assertFalse((self.output / "1132_igv.xml").exists())
        self.assertEqual(main([str(self.reads), "--output-dir", str(self.output), "--igv-only"]), 1)

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_shared_snapgene_reference_combines_all_reads(self):
        self.output = self.reads / "results"
        for name in ("1132.fastq.gz", "1133.fastq", "replicate.fq"):
            self.fastq(name)
        reference = self.snapgene("shared plasmid & construct.dna")
        result = self.run_batch(reference_file=reference)
        self.assertTrue(result["status"].eq("aligned").all())
        session = self.output / "shared plasmid & construct_all_reads_igv.xml"
        self.assertTrue(result["combined_igv_session"].eq(str(session)).all())
        root = ET.parse(session).getroot()
        tracks = root.findall("Resources/Resource")
        self.assertEqual(len(tracks), 3)
        self.assertEqual({track.get("name") for track in tracks}, {"1132", "1133", "replicate"})
        reference_path = self.output / root.get("genome")
        self.assertTrue(reference_path.is_file())
        for track in tracks:
            bam = self.output / track.get("path")
            header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
            self.assertIn("SN:shared_plasmid_&_construct", header)
        before = {path.name: path.stat().st_mtime_ns for path in self.output.glob("*_sorted.bam")}
        self.assertTrue(self.run_batch(reference_file=reference)["action"].eq("skipped").all())
        self.assertEqual(before, {path.name: path.stat().st_mtime_ns for path in self.output.glob("*_sorted.bam")})
        session.unlink()
        self.assertTrue(create_igv_sessions(self.output)["combined_igv_session"].eq(str(session)).all())

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "Alignment tools not installed")
    def test_shared_fasta_combined_session_includes_new_fastq(self):
        reference = self.fasta("shared.fa")
        for name in ("1132.fastq", "1133.fastq"):
            self.fastq(name)
        self.run_batch(reference_file=reference)
        session = self.output / "shared_all_reads_igv.xml"
        self.assertEqual(len(ET.parse(session).findall("Resources/Resource")), 2)
        self.fastq("1134.fastq")
        self.run_batch(reference_file=reference)
        self.assertEqual(len(ET.parse(session).findall("Resources/Resource")), 3)
        self.assertEqual(main([str(self.reads), "--output-dir", str(self.output),
                               "--reference-file", str(reference), "--threads", "2"]), 0)


if __name__ == "__main__":
    unittest.main()
