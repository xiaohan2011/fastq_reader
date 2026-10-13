# FASTQ reader

Analyze Nanopore FASTQ files with a length-sorted CSV and read-length histogram,
and optionally align reads to a reference with minimap2 and samtools.
The original single-file notebook (`1_fastq_reader.ipynb`) and interactive shell
helper (`run_minimap2_ont.sh`) remain available.

## Batch processing

`batch_fastq_reader.py` scans **only files directly inside the selected folders**.
It accepts `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files (case insensitive).
Reference formats are `.fa`, `.fasta`, `.fna`, `.fas`, and binary SnapGene `.dna`.

References match by the **complete leading number** in each filename:

| FASTQ | Reference | Match |
| --- | --- | --- |
| `1132.fastq.gz` | `1132 LPPC EE R-S.dna` | Yes |
| `1133.fastq.gz` | `1133 LPPC EE R-T.fa` | Yes |
| `1132.fastq.gz` | `11320 plasmid.fa` | No |

Numbers are compared as written (`001132` and `1132` are different IDs). Every
FASTQ gets its own outputs, including files sharing a sample ID. Missing
references or filenames without a leading number still get QC. Multiple matching
references produce QC and an `ambiguous_reference` result; alignment is skipped
until the reference folder has exactly one match for that ID.

### Installation

Use Python 3.10 or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Alignment additionally requires `minimap2` and `samtools` on `PATH`. For example,
install them with `brew install minimap2 samtools` on macOS, or
`sudo apt-get install minimap2 samtools` on Debian/Ubuntu. QC alone does not
require those tools. SnapGene extraction uses Biopython and does not require the
SnapGene application.

In the prepared cloud environment, dependencies are already installed:

```bash
source /workspace/fastq-reader-local/activate.sh
cd /workspace/fastq_reader
```

### Command line

```bash
python batch_fastq_reader.py "/path/to/fastq folder" \
  --reference-dir "/path/to/reference folder" \
  --output-dir "/path/to/results"
```

Omit `--reference-dir` to search the FASTQ folder for references. Use `--qc-only`
to disable reference search and alignment entirely. Other options:

```bash
python batch_fastq_reader.py --help
# --overwrite     Rebuild outputs, including previously completed samples
# --threads 4     Alignment threads per sample (samples run sequentially)
# --bin-size 2000 Histogram bin size in base pairs
# --log-y         Logarithmic histogram y axis
```

### Notebook

Install `requirements-notebook.txt` instead of `requirements.txt` if you need
Jupyter, then open `2_batch_fastq_reader.ipynb` from the repository folder. Set
`FASTQ_DIR`, `REFERENCE_DIR`, and `OUTPUT_DIR` in the settings cell and run the
notebook. It calls the same implementation as the command line.

### Outputs and reruns

All generated results are written **directly into one separate output folder**:

```text
results/
  1132_summary.csv
  1132_length_hist.png
  1132_reference.fasta
  1132_sorted.bam
  1132_sorted.bam.bai
  1132_alignment.log
  1133_summary.csv
  1133_length_hist.png
  batch_report.csv
  .1132_complete.json
  .1133_complete.json
```

QC columns match the original notebook: `Sequence`, `Average Quality` (mean
Phred+33 score), and `Length`, sorted by descending read length. Reference FASTA,
BAM, BAM index, and alignment log are created for aligned samples. SnapGene
sequences are extracted into FASTA; annotations and circular topology are not
used by the `map-ont` alignment, matching the original linear FASTA workflow.

Output prefixes normally use the FASTQ stem. If both `1132.fastq` and
`1132.fastq.gz` are present, their full filenames are used as prefixes, such as
`1132.fastq.gz_summary.csv`. Keep input filenames stable across reruns.

Completion manifests track input/reference paths, sizes, modification times,
settings, and expected outputs. A rerun skips a completed sample when these
match and its artifacts still exist and are nonempty. This is a practical
resume check, not a content-hash or BAM integrity check on every rerun.
Incomplete outputs, changed inputs/settings, or new reference matches require
`--overwrite`. Overwrite replaces only the selected samples' known generated
files; unrelated files remain. Do not run multiple batches concurrently into
the same output folder.

Failures are reported per file and do not stop the other samples. If alignment
fails, successful CSV/histogram outputs and any alignment log are kept. Empty,
truncated, or malformed four-line FASTQ records are errors rather than silently
discarded reads. `batch_report.csv` is replaced with the latest run's report;
its `action` column distinguishes processed from skipped samples. The CLI exits
with status 1 for failures or ambiguous references, 0 for a successful batch
(including missing-reference QC), and 2 for invalid arguments/folders.

Like the original notebook, one FASTQ's sequences are held in memory for sorting
and CSV export. Samples run one at a time; very large individual FASTQs need
enough RAM for their sequences and dataframe.

### Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover compressed/plain FASTQ QC, leading-number matching, top-level-only
scanning, missing/ambiguous references, failures, output collisions, resume and
overwrite, and binary SnapGene extraction plus real alignment. The alignment
test is skipped if minimap2 or samtools is absent.

