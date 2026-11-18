#!/usr/bin/env python3
"""Independent long-read assembly with Flye, read-back QC and optional comparison."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import sys

from Bio import SeqIO

from batch_fastq_reader import align, fastq_stem, fingerprint, prepare_reference, write_igv_session

READ_TYPES = ("nano-hq", "nano-raw", "nano-corr", "pacbio-hifi", "pacbio-corr", "pacbio-raw")


def read_qc(path: Path) -> dict:
    """Validate four-line FASTQ while collecting read stats without retaining sequences."""
    lengths, quality_sum, bases = [], 0, 0
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, separator, quality = (handle.readline() for _ in range(3))
            if not sequence or not separator or not quality:
                raise ValueError(f"Truncated FASTQ record {len(lengths) + 1}")
            sequence, quality = sequence.rstrip("\r\n"), quality.rstrip("\r\n")
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError(f"Invalid FASTQ record {len(lengths) + 1}")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch at record {len(lengths) + 1}")
            scores = [ord(c) - 33 for c in quality]
            if any(score < 0 or score > 93 for score in scores):
                raise ValueError("FASTQ quality scores must use Phred+33 encoding")
            lengths.append(len(sequence))
            bases += len(sequence)
            quality_sum += sum(scores)
    if not lengths:
        raise ValueError("FASTQ contains no reads")
    cumulative = 0
    for n50 in sorted(lengths, reverse=True):
        cumulative += n50
        if cumulative >= bases / 2:
            break
    return {"read_count": len(lengths), "total_bases": bases, "read_n50_bp": n50,
            "median_read_length_bp": statistics.median(lengths),
            "mean_phred_score": quality_sum / bases, "longest_read_bp": max(lengths)}


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def run_command(command: list[str], log: Path, commands: list, *, stdout: Path | None = None) -> None:
    commands.append(command)
    with log.open("a") as errors:
        if stdout is None:
            result = subprocess.run(command, stdout=errors, stderr=subprocess.STDOUT)
        else:
            with stdout.open("w") as output:
                result = subprocess.run(command, stdout=output, stderr=errors)
    if result.returncode:
        raise RuntimeError(f"{Path(command[0]).name} exited {result.returncode}; see {log}")


def compare_reference(reference: Path, assembly: Path, output: Path, commands: list) -> list[dict]:
    """PAF alignment blocks, including split/repeat mappings; not a variant caller."""
    paf = output / "assembly_vs_reference.paf"
    run_command(["minimap2", "-x", "asm5", "-c", "--cs=long", str(reference), str(assembly)],
                output / "comparison.log", commands, stdout=paf)
    rows = []
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise ValueError("Malformed minimap2 PAF output")
            rows.append({"assembly_contig": fields[0], "contig_length": int(fields[1]),
                         "query_start_0based": int(fields[2]), "query_end_0based": int(fields[3]),
                         "strand": fields[4], "reference_contig": fields[5],
                         "reference_length": int(fields[6]), "ref_start_0based": int(fields[7]),
                         "ref_end_0based": int(fields[8]), "matching_bases": int(fields[9]),
                         "alignment_block_length": int(fields[10]), "mapq": int(fields[11]),
                         "identity_percent": 100 * int(fields[9]) / int(fields[10]) if int(fields[10]) else 0})
    fieldnames = list(rows[0]) if rows else ["assembly_contig", "reference_contig", "identity_percent"]
    with (output / "reference_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run_assembly(fastq_file: str | Path, output_dir: str | Path | None = None,
                 reference_file: str | Path | None = None, *, read_type: str = "nano-hq",
                 threads: int = 4, genome_size: str | None = None, min_overlap: int | None = None,
                 polish_iterations: int = 2, resume: bool = False) -> dict:
    """Assemble reads independently; the optional reference is used only afterward."""
    fastq = Path(fastq_file).expanduser().resolve()
    prefix = fastq_stem(fastq)
    if not fastq.is_file():
        raise ValueError(f"FASTQ does not exist: {fastq}")
    reference = Path(reference_file).expanduser().resolve() if reference_file else None
    if reference is not None and (not reference.is_file() or reference.suffix.lower() not in (".dna", ".fa", ".fasta", ".fna", ".fas")):
        raise ValueError(f"Reference must be an existing SnapGene or FASTA file: {reference}")
    if read_type not in READ_TYPES or threads < 1 or polish_iterations < 1:
        raise ValueError("Choose a supported read_type and positive threads/polish_iterations")
    if min_overlap is not None and min_overlap < 1000:
        raise ValueError("Flye min_overlap must be at least 1000 bp")
    output = Path(output_dir).expanduser().resolve() if output_dir else fastq.parent / f"{prefix}_denovo"
    if any(path.is_relative_to(output) for path in (fastq, reference) if path is not None):
        raise ValueError("Use a dedicated output folder that does not contain your input files")
    signature = {"workflow_version": 1, "fastq": fingerprint(fastq),
                 "reference": fingerprint(reference) if reference else None, "read_type": read_type,
                 "threads": threads, "genome_size": genome_size, "min_overlap": min_overlap,
                 "polish_iterations": polish_iterations}
    state_path = output / "run_state.json"
    if output.exists() and any(output.iterdir()):
        if not resume or not state_path.is_file():
            raise ValueError("Output folder is not empty; choose a new folder or use resume=True for this run")
        saved = json.loads(state_path.read_text())
        if saved.get("signature") != signature:
            raise ValueError("Inputs/settings changed; choose a new output folder instead of resuming")
        if saved.get("status") == "complete" and all((output / p).is_file() for p in saved.get("artifacts", [])):
            return json.loads((output / "summary.json").read_text())
    for tool in ("flye", "minimap2", "samtools"):
        if not shutil.which(tool):
            raise RuntimeError(f"{tool} is missing from PATH; see ASSEMBLY.md for installation")
    qc = read_qc(fastq)
    output.mkdir(parents=True, exist_ok=True)
    commands = []
    state = {"signature": signature, "status": "running", "commands": commands}
    write_json(state_path, state)
    try:
        write_json(output / "read_qc.json", qc)
        versions = {}
        for tool in ("flye", "minimap2", "samtools"):
            option = "--version-only" if tool == "samtools" else "--version"
            versions[tool] = subprocess.check_output([tool, option], text=True).splitlines()[0]
        flye_dir = output / "flye"
        command = ["flye", f"--{read_type}", str(fastq), "--out-dir", str(flye_dir),
                   "--threads", str(threads), "--iterations", str(polish_iterations)]
        if genome_size:
            command += ["--genome-size", str(genome_size)]
        if min_overlap is not None:
            command += ["--min-overlap", str(min_overlap)]
        if resume and (flye_dir / "params.json").is_file():
            command += ["--resume"]
        # No reference sequence/path is supplied to Flye.
        run_command(command, output / "assembly.log", commands)
        source = flye_dir / "assembly.fasta"
        with source.open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        if not records or any(not len(record.seq) for record in records):
            raise RuntimeError("Flye did not produce a nonempty assembly FASTA")
        assembly = output / f"{prefix}_assembly.fasta"
        shutil.copyfile(source, assembly)
        info = {}
        with (flye_dir / "assembly_info.txt").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                info[row["#seq_name"]] = row
        bam = output / f"{prefix}_reads_to_assembly.bam"
        preset = "map-hifi" if read_type == "pacbio-hifi" else "map-pb" if read_type.startswith("pacbio") else "map-ont"
        commands.append(["minimap2", "-t", str(threads), "-ax", preset, str(assembly), str(fastq)])
        align(assembly, fastq, bam, output / "read_back_alignment.log", threads, preset=preset)
        coverage_file = output / "read_back_coverage.tsv"
        run_command(["samtools", "coverage", str(bam)], output / "coverage.log", commands, stdout=coverage_file)
        coverage = {}
        with coverage_file.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                coverage[row["#rname"]] = row
        contigs = []
        for record in records:
            cov = coverage.get(record.id, {})
            contigs.append({"contig": record.id, "length_bp": len(record.seq),
                            "ambiguous_bases": sum(base not in "ACGT" for base in str(record.seq).upper()),
                            "flye_reports_circular": info.get(record.id, {}).get("circ.") == "Y",
                            "flye_coverage": float(info.get(record.id, {}).get("cov.", 0)),
                            "read_back_mean_depth": float(cov.get("meandepth", 0)),
                            "bases_covered_percent": float(cov.get("coverage", 0))})
        igv = output / f"{prefix}_assembly_igv.xml"
        write_igv_session(assembly, bam, igv)
        comparisons = []
        artifacts = [assembly.name, bam.name, bam.name + ".bai", assembly.name + ".fai", igv.name,
                     "read_qc.json", coverage_file.name, "summary.json", "REPORT.md"]
        if reference is not None:
            converted = output / "comparison_reference.fasta"
            prepare_reference(reference, converted, "comparison_reference")
            comparisons = compare_reference(converted, assembly, output, commands)
            artifacts += [converted.name, "assembly_vs_reference.paf", "reference_comparison.csv"]
        notes = ["All Flye contigs are retained; no reference bases were inserted into the assembly.",
                 "Circularity is Flye's graph classification, not independent proof of a complete plasmid.",
                 "Read-back coverage does not establish per-base accuracy or exclude mixed plasmids."]
        if len(contigs) != 1:
            notes.append("Multiple contigs: do not treat the output as a single finished plasmid.")
        if any(c["read_back_mean_depth"] < 10 or c["bases_covered_percent"] < 100 for c in contigs):
            notes.append("At least one contig has low mean depth (<10x) or uncovered bases; inspect coverage before use.")
        if reference is not None and not comparisons:
            notes.append("No assembly-to-reference alignment was found; inspect read type, sample identity and assembly quality.")
        summary = {"assembly_fasta": str(assembly), "igv_session": str(igv), "read_qc": qc,
                   "contigs": contigs, "reference_alignment_blocks": len(comparisons),
                   "tool_versions": versions, "notes": notes}
        write_json(output / "summary.json", summary)
        lines = ["# De novo assembly report", "", f"Assembly: `{assembly.name}`", "",
                 f"Reads: {qc['read_count']:,}; bases: {qc['total_bases']:,}; read N50: {qc['read_n50_bp']:,} bp", "",
                 "| Contig | Length (bp) | Flye circular | Mean read depth | Bases covered (%) |",
                 "| --- | ---: | --- | ---: | ---: |"]
        for contig in contigs:
            lines.append(f"| {contig['contig']} | {contig['length_bp']} | {contig['flye_reports_circular']} | "
                         f"{contig['read_back_mean_depth']:.2f} | {contig['bases_covered_percent']:.2f} |")
        lines += ["", *[f"- {note}" for note in notes]]
        if reference is not None:
            lines += ["", "Reference comparison: `reference_comparison.csv` and `assembly_vs_reference.paf`.",
                      "Coordinates are zero-based, half-open. Rows are local alignment blocks, not independent variants.",
                      "Circular sequences may align in split blocks because their starting positions differ."]
        (output / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        state.update(status="complete", artifacts=artifacts, tool_versions=versions)
        write_json(state_path, state)
        return summary
    except Exception as error:
        state.update(status="failed", error=f"{type(error).__name__}: {error}")
        write_json(state_path, state)
        raise


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--reference-file", type=Path, help="Comparison only; not supplied to the assembler")
    parser.add_argument("--read-type", choices=READ_TYPES, default="nano-hq")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--genome-size", help="Optional size estimate, e.g. 105k; does not supply reference sequence")
    parser.add_argument("--min-overlap", type=int)
    parser.add_argument("--polish-iterations", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary = run_assembly(**vars(args))
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"Assembly failed: {error}", file=sys.stderr)
        return 1
    print("Assembly FASTA:", summary["assembly_fasta"])
    print("Read-back IGV session:", summary["igv_session"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
