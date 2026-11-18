"""Edit these settings and click Run in PyCharm. See ASSEMBLY.md for tool setup."""

from pathlib import Path
import os

from denovo_assembly import run_assembly

# %% Settings
FASTQ_FILE = Path("/path/to/reads.fastq.gz")
REFERENCE_FILE = None  # optional Path("/path/to/reference.dna") or FASTA, for comparison only
OUTPUT_DIR = None  # default: <FASTQ filename>_denovo beside the FASTQ
READ_TYPE = "nano-hq"  # modern ONT; use nano-raw for older/noisier ONT
THREADS = 4
GENOME_SIZE = None  # optional size estimate, e.g. "105k"; not a reference sequence
MIN_OVERLAP = None  # Flye default; use 1000 for short plasmids/reads if needed
POLISH_ITERATIONS = 2
RESUME = False  # True reuses an existing run only if inputs/settings match
# If the IDE cannot see your tools, add their bin directories here.
TOOL_DIRS = []  # e.g. [Path("/opt/homebrew/bin"), Path("/path/to/Flye/bin")]


# %% Run
def main():
    if TOOL_DIRS:
        os.environ["PATH"] = os.pathsep.join(str(Path(p).expanduser()) for p in TOOL_DIRS) + os.pathsep + os.environ.get("PATH", "")
    summary = run_assembly(FASTQ_FILE, OUTPUT_DIR, REFERENCE_FILE, read_type=READ_TYPE,
                           threads=THREADS, genome_size=GENOME_SIZE, min_overlap=MIN_OVERLAP,
                           polish_iterations=POLISH_ITERATIONS, resume=RESUME)
    print("Assembly FASTA:", summary["assembly_fasta"])
    print("IGV session:", summary["igv_session"])
    return summary


if __name__ == "__main__":
    main()
