"""Edit these settings and click Run in PyCharm. See ASSEMBLY.md for tool setup."""

from pathlib import Path

from denovo_assembly import run_assembly
from local_settings import local_path

# %% Settings
FASTQ_FILE = local_path("assembly", "fastq_file")
REFERENCE_FILE = local_path("assembly", "reference_file")
OUTPUT_DIR = local_path("assembly", "output_dir")
READ_TYPE = "nano-hq"  # modern ONT; use nano-raw for older/noisier ONT
THREADS = 4
GENOME_SIZE = "100000"  # Optional size estimate
ASM_COVERAGE = 100  # Optional initial coverage limit (requires GENOME_SIZE)
MIN_OVERLAP = None  # Let Flye choose the initial minimum overlap
POLISH_ITERATIONS = 2
RESUME = False  # True reuses an existing run only if inputs/settings match
# %% Run
def main():
    if FASTQ_FILE is None:
        raise ValueError("Set assembly.fastq_file in .fastq-reader.local.json, or set FASTQ_FILE here")
    summary = run_assembly(FASTQ_FILE, OUTPUT_DIR, REFERENCE_FILE, read_type=READ_TYPE,
                           threads=THREADS, genome_size=GENOME_SIZE, min_overlap=MIN_OVERLAP,
                           polish_iterations=POLISH_ITERATIONS, resume=RESUME, asm_coverage=ASM_COVERAGE)
    print("Assembly FASTA:", summary["assembly_fasta"])
    print("IGV session:", summary["igv_session"])
    return summary


if __name__ == "__main__":
    main()
