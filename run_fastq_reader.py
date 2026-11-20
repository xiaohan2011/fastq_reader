"""Edit the settings below, then run this file in PyCharm or with Python.

Keep batch_fastq_reader.py in the same folder. This is ordinary Python, not
notebook JSON; it does not require Jupyter or a special notebook editor.
"""

from __future__ import annotations

import importlib
from pathlib import Path
import sys

from local_settings import local_path


# %% Settings — edit these values, then run the entire file.
# Folder mode: set FASTQ_DIR and leave FASTQ_FILE = None.
FASTQ_DIR = local_path("batch", "fastq_dir")
FASTQ_FILE = local_path("batch", "fastq_file")  # or Path("/path/to/selected.fastq.gz") for just one file

# Automatic matching: set REFERENCE_DIR, or leave None to search the FASTQ folder.
# Explicit pairing: set FASTQ_FILE above and REFERENCE_FILE below.
REFERENCE_DIR = local_path("batch", "reference_dir")  # or Path("/path/to/reference_folder")
REFERENCE_FILE = local_path("batch", "reference_file")  # or Path("/path/to/reference.dna") / Path("/path/to/reference.fa")

# None creates a results subfolder beside the selected FASTQ, or inside FASTQ_DIR.
OUTPUT_DIR = local_path("batch", "output_dir")  # or Path("/path/to/custom/results")
OVERWRITE = False
QC_ONLY = False
BIN_SIZE_BP = 2000
LOG_Y = False
THREADS = 4

# Windows paths: use raw strings, e.g. Path(r"C:\Users\yourname\data").
# macOS example: Path("/Volumes/your_share/fastqs").


# %% Run analysis — normally no changes are needed below this line.
def run_from_settings():
    """Run the selected folder or file and return the same report as the notebook."""
    import batch_fastq_reader

    # Also pick up updates when this section is reused in a Python console.
    importlib.invalidate_caches()
    reader = importlib.reload(batch_fastq_reader)
    print("Using batch code:", reader.__file__)
    if OUTPUT_DIR is not None:
        output = Path(OUTPUT_DIR)
    elif FASTQ_FILE is not None:
        output = Path(FASTQ_FILE).expanduser().parent / "results"
    elif FASTQ_DIR is not None:
        output = Path(FASTQ_DIR).expanduser() / "results"
    else:
        raise ValueError("Set FASTQ_DIR or FASTQ_FILE in the Settings section")

    report = reader.run_batch(
        fastq_dir=FASTQ_DIR,
        fastq_file=FASTQ_FILE,
        reference_dir=REFERENCE_DIR,
        reference_file=REFERENCE_FILE,
        output_dir=output,
        overwrite=OVERWRITE,
        qc_only=QC_ONLY,
        bin_size=BIN_SIZE_BP,
        log_y=LOG_Y,
        threads=THREADS,
    )
    # Ordinary Python does not provide Jupyter's display() function.
    print(report.to_string(index=False))
    return report


def main() -> int:
    try:
        report = run_from_settings()
    except (ValueError, OSError) as error:
        print(f"Configuration or file error: {error}", file=sys.stderr)
        return 2
    return int(report["status"].isin(["failed", "ambiguous_reference"]).any())


if __name__ == "__main__":
    raise SystemExit(main())
