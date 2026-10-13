#!/usr/bin/env python3
"""Batch FASTQ QC and optional ONT reference alignment (top-level folders only)."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
REFERENCE_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".dna")
WORKFLOW_VERSION = 1
REPORT_FIELDS = (
    "sample_id", "fastq", "reference", "action", "status", "reads",
    "summary_csv", "histogram_png", "bam", "bai", "alignment_log", "message",
)


def sample_id(path: Path) -> str:
    """Match the complete leading digit sequence, never a numeric prefix."""
    match = re.match(r"\d+", path.name)
    return match.group() if match else ""


def fastq_stem(path: Path) -> str:
    for extension in FASTQ_EXTENSIONS:
        if path.name.lower().endswith(extension):
            return path.name[:-len(extension)]
    raise ValueError(f"Not a supported FASTQ filename: {path}")


def discover(folder: Path, extensions: tuple[str, ...]) -> list[Path]:
    return sorted(
        (p for p in folder.iterdir() if p.is_file() and p.name.lower().endswith(extensions)),
        key=lambda p: p.name,
    )


def read_fastq(path: Path) -> pd.DataFrame:
    """Read strict four-line FASTQ records, including compressed input."""
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="ascii") as handle:
        record = 0
        while True:
            header = handle.readline()
            if not header:
                break
            record += 1
            sequence, plus, quality = (handle.readline() for _ in range(3))
            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record {record} in {path.name}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ header/separator at record {record}")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch at record {record}")
            if any(ord(c) < 33 or ord(c) > 126 for c in quality):
                raise ValueError(f"Invalid quality character at record {record}")
            rows.append((sequence, sum(ord(c) - 33 for c in quality) / len(quality), len(sequence)))
    if not rows:
        raise ValueError(f"No FASTQ records found in {path.name}")
    return pd.DataFrame(rows, columns=["Sequence", "Average Quality", "Length"]).sort_values(
        "Length", ascending=False, kind="stable",
    )


def write_qc(df: pd.DataFrame, summary: Path, histogram: Path, bin_size: int, log_y: bool) -> None:
    df.to_csv(summary, index=False)
    max_length = int(df["Length"].max())
    # Always include an upper edge beyond the longest read.
    edges = range(0, (max_length // bin_size + 1) * bin_size + 1, bin_size)
    fig, ax = plt.subplots()
    try:
        ax.hist(df["Length"] / 1000, bins=[edge / 1000 for edge in edges])
        ax.set(xlabel="Read length (kb)", ylabel="Number of reads",
               title=f"Read-length distribution (bin = {bin_size / 1000:g} kb)")
        if log_y:
            ax.set_yscale("log")
        fig.tight_layout()
        fig.savefig(histogram, dpi=300)
    finally:
        plt.close(fig)


def prepare_reference(source: Path, destination: Path, identifier: str) -> None:
    """Validate FASTA or extract the DNA sequence from a binary SnapGene file."""
    from Bio import SeqIO

    if source.suffix.lower() == ".dna":
        with source.open("rb") as handle:
            record = SeqIO.read(handle, "snapgene")
        record.id = identifier
        record.name = identifier
        record.description = "sequence extracted from SnapGene"
        records = [record]
    else:
        with source.open("r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    if not records or any(not len(record.seq) for record in records):
        raise ValueError(f"Reference contains no usable sequence: {source.name}")
    if len({record.id for record in records}) != len(records):
        raise ValueError(f"Reference has duplicate sequence IDs: {source.name}")
    SeqIO.write(records, destination, "fasta")


def align(reference: Path, fastq: Path, bam: Path, log: Path, threads: int) -> None:
    for tool in ("minimap2", "samtools"):
        if not shutil.which(tool):
            raise RuntimeError(f"{tool} is missing from PATH; install it to enable alignment")
    with log.open("w") as errors:
        # Write stderr directly to disk to avoid a blocked stderr pipe on large jobs.
        mapper = subprocess.Popen(
            ["minimap2", "-t", str(threads), "-ax", "map-ont", str(reference), str(fastq)],
            stdout=subprocess.PIPE, stderr=errors,
        )
        try:
            sorter = subprocess.Popen(
                ["samtools", "sort", "-@", str(max(1, threads - 1)), "-o", str(bam)],
                stdin=mapper.stdout, stderr=errors,
            )
            mapper.stdout.close()
            sort_code = sorter.wait()
            map_code = mapper.wait()
        finally:
            if mapper.stdout:
                mapper.stdout.close()
            if mapper.poll() is None:
                mapper.terminate()
                mapper.wait()
        if map_code or sort_code:
            raise RuntimeError(f"Alignment failed (minimap2={map_code}, samtools={sort_code}); see {log.name}")
        subprocess.run(["samtools", "index", str(bam)], stderr=errors, check=True)
        subprocess.run(["samtools", "quickcheck", str(bam)], stderr=errors, check=True)


def fingerprint(path: Path) -> dict:
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def process_sample(
    fastq: Path, candidates: list[Path], output: Path, prefix: str, *,
    bin_size: int, log_y: bool, threads: int, overwrite: bool, qc_only: bool,
) -> dict:
    identifier = sample_id(fastq)
    reference = candidates[0] if len(candidates) == 1 else None
    row = dict.fromkeys(REPORT_FIELDS, "")
    row.update(sample_id=identifier, fastq=str(fastq), reference=str(reference) if reference else "",
               action="processed", status="failed")
    names = {
        "summary_csv": f"{prefix}_summary.csv",
        "histogram_png": f"{prefix}_length_hist.png",
        "reference_fasta": f"{prefix}_reference.fasta",
        "bam": f"{prefix}_sorted.bam",
        "bai": f"{prefix}_sorted.bam.bai",
        "alignment_log": f"{prefix}_alignment.log",
    }
    marker = output / f".{prefix}_complete.json"
    signature = {
        "version": WORKFLOW_VERSION, "fastq": fingerprint(fastq),
        "references": [fingerprint(p) for p in candidates],
        "bin_size": bin_size, "log_y": log_y, "threads": threads, "qc_only": qc_only,
    }
    if not overwrite:
        if marker.is_file():
            try:
                saved = json.loads(marker.read_text())
                valid = all((output / name).is_file() and (output / name).stat().st_size > 0
                            for name in saved["artifacts"])
                if saved["signature"] == signature and valid:
                    return dict(saved["result"], action="skipped")
            except (ValueError, KeyError, TypeError):
                pass
        if marker.exists() or any((output / name).exists() for name in names.values()):
            row["message"] = "Existing outputs are incomplete or inputs/settings changed; use --overwrite to rebuild"
            return row
    # Keep each sample's work isolated until an artifact is complete.
    with tempfile.TemporaryDirectory(prefix=".fastq-work-", dir=output) as temporary:
        work = Path(temporary)
        artifacts = []

        def publish(key: str) -> None:
            name = names[key]
            (work / name).replace(output / name)
            artifacts.append(name)
            if key in row:
                row[key] = str(output / name)

        if overwrite:
            marker.unlink(missing_ok=True)
            # Remove only this sample's known generated files, never a directory tree.
            for name in names.values():
                (output / name).unlink(missing_ok=True)
        try:
            df = read_fastq(fastq)
            row["reads"] = len(df)
            write_qc(df, work / names["summary_csv"], work / names["histogram_png"], bin_size, log_y)
            publish("summary_csv")
            publish("histogram_png")
            if qc_only:
                row.update(status="qc_only", message="Alignment disabled")
            elif not identifier:
                row.update(status="qc_only_no_sample_id", message="Filename has no leading number; no reference matched")
            elif len(candidates) > 1:
                row.update(status="ambiguous_reference", message="QC complete; alignment skipped. Matches: " +
                           "; ".join(p.name for p in candidates))
            elif reference is None:
                row.update(status="qc_only_missing_reference", message="No matching reference; QC complete")
            else:
                prepare_reference(reference, work / names["reference_fasta"], identifier)
                align(work / names["reference_fasta"], fastq, work / names["bam"],
                      work / names["alignment_log"], threads)
                for key in ("reference_fasta", "bam", "bai", "alignment_log"):
                    publish(key)
                row.update(status="aligned", message="QC and alignment complete")
        except Exception as error:
            row.update(status="failed", message=f"{type(error).__name__}: {error}")
            # Preserve diagnostics and successful QC even when alignment fails.
            if (work / names["alignment_log"]).exists():
                publish("alignment_log")
            return row
        atomic_json(marker, {"signature": signature, "artifacts": artifacts, "result": row})
    return row


def run_batch(
    fastq_dir: str | Path, output_dir: str | Path, reference_dir: str | Path | None = None,
    *, overwrite: bool = False, qc_only: bool = False, bin_size: int = 2000,
    log_y: bool = False, threads: int = 4,
) -> pd.DataFrame:
    """Process all top-level FASTQs; return and save the per-file batch report."""
    fastq_dir = Path(fastq_dir).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    reference_dir = Path(reference_dir).expanduser().resolve() if reference_dir else fastq_dir
    if not fastq_dir.is_dir():
        raise ValueError(f"FASTQ folder does not exist: {fastq_dir}")
    if not qc_only and not reference_dir.is_dir():
        raise ValueError(f"Reference folder does not exist: {reference_dir}")
    if output == fastq_dir or (not qc_only and output == reference_dir):
        raise ValueError("Use a separate output folder, not an input/reference folder")
    if bin_size < 1 or threads < 1:
        raise ValueError("bin_size and threads must be positive integers")
    fastqs = discover(fastq_dir, FASTQ_EXTENSIONS)
    if not fastqs:
        raise ValueError(f"No supported FASTQ files in {fastq_dir} (subfolders are not scanned)")
    references = defaultdict(list)
    if not qc_only:
        for path in discover(reference_dir, REFERENCE_EXTENSIONS):
            if sample_id(path):
                references[sample_id(path)].append(path)
    stems = Counter(fastq_stem(path) for path in fastqs)
    # Distinguish same-stem inputs such as 1132.fastq and 1132.fastq.gz.
    prefixes = [fastq_stem(path) if stems[fastq_stem(path)] == 1 else path.name for path in fastqs]
    if len(set(prefixes)) != len(prefixes):
        raise ValueError("FASTQ filenames cause output name collisions; rename the conflicting inputs")
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for path, prefix in zip(fastqs, prefixes):
        try:
            result = process_sample(
                path, references[sample_id(path)], output, prefix, bin_size=bin_size,
                log_y=log_y, threads=threads, overwrite=overwrite, qc_only=qc_only,
            )
        except Exception as error:
            result = dict.fromkeys(REPORT_FIELDS, "")
            result.update(sample_id=sample_id(path), fastq=str(path), action="processed", status="failed",
                          message=f"{type(error).__name__}: {error}")
        results.append(result)
        print(f"[{result['action']}] {path.name}: {result['status']} — {result['message']}")
    report = pd.DataFrame(results, columns=REPORT_FIELDS)
    temporary_report = output / ".batch_report.csv.tmp"
    report.to_csv(temporary_report, index=False, quoting=csv.QUOTE_MINIMAL)
    temporary_report.replace(output / "batch_report.csv")
    print(f"Batch report: {output / 'batch_report.csv'}")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_dir", type=Path, help="Folder containing FASTQ files; no recursive search")
    parser.add_argument("--reference-dir", type=Path, help="FASTA/SnapGene folder (default: FASTQ folder)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Separate folder for all generated files")
    parser.add_argument("--overwrite", action="store_true", help="Rebuild this batch's outputs")
    parser.add_argument("--qc-only", action="store_true", help="Skip reference search and alignment")
    parser.add_argument("--bin-size", type=int, default=2000, help="Histogram bin size in bp (default: 2000)")
    parser.add_argument("--log-y", action="store_true", help="Use a logarithmic histogram y axis")
    parser.add_argument("--threads", type=int, default=min(4, os.cpu_count() or 1), help="Alignment threads per sample")
    args = parser.parse_args(argv)
    try:
        report = run_batch(args.fastq_dir, args.output_dir, args.reference_dir,
                           overwrite=args.overwrite, qc_only=args.qc_only, bin_size=args.bin_size,
                           log_y=args.log_y, threads=args.threads)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    return int(report["status"].isin(["failed", "ambiguous_reference"]).any())


if __name__ == "__main__":
    raise SystemExit(main())
