"""Edit these settings and click Run in PyCharm. See ASSEMBLY.md for tool setup."""

from pathlib import Path
import os

from denovo_assembly import run_assembly

# %% Settings
FASTQ_FILE = Path('/Users/wangx21/Desktop/1147/J792H6_1_1147-F2-2.fastq')
REFERENCE_FILE = Path('/Users/wangx21/Library/CloudStorage/Dropbox-NYULangoneHealth/Xiaohan Wang/XW_reagents/1_DNA&Protein/1_pXW_record/1147_Hba_EF1a_CT3_SEQinv.dna')  # optional Path("/path/to/reference.dna") or FASTA, for comparison only
OUTPUT_DIR = None  # Default: <FASTQ stem>_denovo beside the FASTQ
READ_TYPE = "nano-hq"  # modern ONT; use nano-raw for older/noisier ONT
THREADS = 4
GENOME_SIZE = "100000"  # Optional size estimate
ASM_COVERAGE = 100  # Optional initial coverage limit (requires GENOME_SIZE)
MIN_OVERLAP = None  # Let Flye choose the initial minimum overlap
POLISH_ITERATIONS = 2
RESUME = False  # True reuses an existing run only if inputs/settings match
# If the IDE cannot see your tools, add their bin directories here.
TOOL_DIRS = [Path.home() / "PycharmProjects" / ".fastq-flye" / "bin", Path("/usr/local/bin")]


# %% Run
def main():
    if TOOL_DIRS:
        os.environ["PATH"] = os.pathsep.join(str(Path(p).expanduser()) for p in TOOL_DIRS) + os.pathsep + os.environ.get("PATH", "")
    summary = run_assembly(FASTQ_FILE, OUTPUT_DIR, REFERENCE_FILE, read_type=READ_TYPE,
                           threads=THREADS, genome_size=GENOME_SIZE, min_overlap=MIN_OVERLAP,
                           polish_iterations=POLISH_ITERATIONS, resume=RESUME, asm_coverage=ASM_COVERAGE)
    print("Assembly FASTA:", summary["assembly_fasta"])
    print("IGV session:", summary["igv_session"])
    return summary


if __name__ == "__main__":
    main()
