#!/usr/bin/env python3
"""Independent long-read assembly with Flye, read-back QC and optional comparison."""

from __future__ import annotations

import argparse
import csv
from contextlib import ExitStack, contextmanager
from datetime import datetime
import gzip
import json
import re
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Event, Thread
from time import monotonic

from Bio import SeqIO

from local_settings import find_tool

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


def progress(message: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {message}", flush=True)


@contextmanager
def show_progress(label: str):
    progress(label)
    started = monotonic()
    finished = Event()

    def report_while_running():
        while not finished.wait(60):
            elapsed = int(monotonic() - started)
            progress(f"{label} still running ({elapsed // 60}m {elapsed % 60:02d}s elapsed)")

    reporter = Thread(target=report_while_running, daemon=True)
    reporter.start()
    try:
        yield
    finally:
        finished.set()
        reporter.join()


def run_command(command: list[str], log: Path, commands: list, *, stdout: Path | None = None) -> None:
    commands.append(command)
    tool = Path(command[0]).name
    executable = command[0] if Path(command[0]).is_absolute() else find_tool(tool)
    if not executable:
        raise RuntimeError(f"{tool} is missing from PATH or local tool_dirs")
    label = f"Running {tool} (log: {log.name})"
    with show_progress(label), ExitStack() as files:
        errors = files.enter_context(log.open("a", encoding="utf-8"))
        if stdout is None:
            output = subprocess.PIPE
            diagnostics = subprocess.STDOUT
        else:
            output = files.enter_context(stdout.open("wb"))
            diagnostics = subprocess.PIPE
        with subprocess.Popen([executable, *command[1:]], stdout=output, stderr=diagnostics, text=True,
                              encoding="utf-8", errors="replace") as process:
            stream = process.stdout if stdout is None else process.stderr
            assert stream is not None
            for line in stream:
                errors.write(line)
                errors.flush()
                # PyCharm gives stderr its error color. Keep ordinary tool
                # diagnostics on stdout even when the tool wrote them to stderr.
                is_error = re.match(r"^(?:\[[^\]]+\]\s*)*(?:ERROR|FATAL)\b", line.lstrip(), re.I)
                is_error = is_error or re.match(r"^\[E::", line.lstrip())
                if is_error and sys.stderr.isatty():
                    print(f"\033[31m{line.rstrip(chr(10))}\033[0m", file=sys.stderr, flush=True)
                else:
                    print(line, end="", file=sys.stderr if is_error else sys.stdout, flush=True)
            return_code = process.wait()
    if return_code:
        raise RuntimeError(f"{tool} exited {return_code}; see {log}")


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
                 polish_iterations: int = 2, resume: bool = False,
                 asm_coverage: int | None = None) -> dict:
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
    if asm_coverage is not None and (asm_coverage < 1 or not genome_size):
        raise ValueError("Flye asm_coverage requires a positive coverage and a genome size estimate")
    output = Path(output_dir).expanduser().resolve() if output_dir else fastq.parent / f"{prefix}_denovo"
    if any(path.is_relative_to(output) for path in (fastq, reference) if path is not None):
        raise ValueError("Use a dedicated output folder that does not contain your input files")
    signature = {"workflow_version": 1, "fastq": fingerprint(fastq),
                 "reference": fingerprint(reference) if reference else None, "read_type": read_type,
                 "threads": threads, "genome_size": genome_size, "min_overlap": min_overlap,
                 "polish_iterations": polish_iterations}
    if asm_coverage is not None:
        signature["asm_coverage"] = asm_coverage
    if not resume and output.is_dir() and any(output.iterdir()):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        output = output.with_name(f"{output.name}_{stamp}")
        print(f"Output folder already contains files. Writing to {output}", flush=True)
    state_path = output / "run_state.json"
    saved = {}
    if resume and output.is_dir() and any(output.iterdir()):
        if not state_path.is_file():
            raise ValueError("Output folder is not empty; choose a new folder or use resume=True for this run")
        saved = json.loads(state_path.read_text())
        if saved.get("signature") != signature:
            raise ValueError("Inputs/settings changed; choose a new output folder instead of resuming")
        if saved.get("status") == "complete" and all((output / p).is_file() for p in saved.get("artifacts", [])):
            return json.loads((output / "summary.json").read_text())
    tools = {name: find_tool(name) for name in ("flye", "minimap2", "samtools")}
    for tool, executable in tools.items():
        if not executable:
            raise RuntimeError(f"{tool} is missing from PATH or local tool_dirs; see ASSEMBLY.md for installation")
    progress(f"Output folder: {output}")
    with show_progress(f"Checking FASTQ: {fastq.name}"):
        qc = read_qc(fastq)
    progress(f"FASTQ check complete: {qc['total_bases']:,} bases")
    output.mkdir(parents=True, exist_ok=True)
    commands = []
    state = {"signature": signature, "status": "running", "commands": commands}
    write_json(state_path, state)
    try:
        write_json(output / "read_qc.json", qc)
        versions = {}
        for tool, executable in tools.items():
            option = "--version-only" if tool == "samtools" else "--version"
            versions[tool] = subprocess.check_output([executable, option], text=True).splitlines()[0]
        flye_dir = output / "flye"
        retry_settings = saved.get("flye_retry")
        if retry_settings:
            genome_size = retry_settings["genome_size"]
            asm_coverage = retry_settings["asm_coverage"]
            min_overlap = retry_settings["min_overlap"]
            flye_dir = output / "flye_retry"
            state["flye_retry"] = retry_settings
            write_json(state_path, state)
        while True:
            try:
                # Flye rejects spaces in read paths and passes paths to native tools.
                # Keep the actual inputs and checkpoints in their original locations.
                with ExitStack() as cleanup:
                    reads_arg, output_arg = fastq, flye_dir
                    if " " in str(fastq) or " " in str(flye_dir):
                        alias_dir = Path(cleanup.enter_context(TemporaryDirectory(prefix="fastq-flye-")))
                        if " " in str(alias_dir):
                            raise ValueError("Flye requires a temporary directory with no spaces in its path")
                        reads_arg = alias_dir / ("reads.fastq.gz" if fastq.name.lower().endswith(".gz") else "reads.fastq")
                        reads_arg.symlink_to(fastq)
                        flye_dir.mkdir(exist_ok=True)
                        output_arg = alias_dir / "flye"
                        output_arg.symlink_to(flye_dir, target_is_directory=True)
                    command = ["flye", f"--{read_type}", str(reads_arg), "--out-dir", str(output_arg),
                               "--threads", str(threads), "--iterations", str(polish_iterations)]
                    if genome_size:
                        command += ["--genome-size", str(genome_size)]
                    if asm_coverage is not None:
                        command += ["--asm-coverage", str(asm_coverage)]
                    if min_overlap is not None:
                        command += ["--min-overlap", str(min_overlap)]
                    if resume and (flye_dir / "params.json").is_file():
                        command += ["--resume"]
                    run_command(command, output / "assembly.log", commands)
                break
            except RuntimeError:
                flye_log = flye_dir / "flye.log"
                if asm_coverage is not None or flye_dir.name != "flye" or not flye_log.is_file():
                    raise
                log_text = flye_log.read_text(errors="replace")
                if "No disjointigs were assembled" not in log_text:
                    raise
                coverage = re.findall(r"Overlap-based coverage:\s*(\d+)", log_text)
                if not coverage or int(coverage[-1]) <= 100:
                    raise
                if not genome_size:
                    genome_size = str(max(1, round(qc["total_bases"] / int(coverage[-1]))))
                asm_coverage = 100
                min_overlap = 1000 if min_overlap is None else min_overlap
                flye_dir = output / "flye_retry"
                if flye_dir.exists():
                    raise RuntimeError(f"Automatic retry directory already exists: {flye_dir}")
                print(f"Flye produced no disjointigs. Retrying with a {genome_size} bp size "
                      f"estimate and a {asm_coverage}x initial coverage limit.", flush=True)
                state["flye_retry"] = {"genome_size": genome_size, "asm_coverage": asm_coverage,
                                       "min_overlap": min_overlap}
                write_json(state_path, state)
        progress("Flye finished. Reading assembled contigs.")
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
        with show_progress(f"Aligning reads back to {len(records)} assembled contig(s)"):
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
            progress("Preparing the optional reference comparison")
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
        progress(f"Assembly complete. Report: {output / 'REPORT.md'}")
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
    parser.add_argument("--asm-coverage", type=int, help="Limit initial disjointig coverage; requires --genome-size")
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
