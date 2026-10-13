# FASTQ reader

Analyze Nanopore FASTQ files with a length-sorted CSV and read-length histogram,
and optionally align reads to a reference with minimap2 and samtools.
The original single-file notebook (`1_fastq_reader.ipynb`) and interactive shell
helper (`run_minimap2_ont.sh`) remain available.

## Batch processing

`batch_fastq_reader.py` scans **only files directly inside the selected folders**.
It accepts `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files (case insensitive).
Reference formats are `.fa`, `.fasta`, `.fna`, `.fas`, and binary SnapGene `.dna`.

References match by the **sample ID**, either a complete leading number or the
third field in a `RUN_BARCODE_ID` FASTQ filename:

| FASTQ | Reference | Match |
| --- | --- | --- |
| `1132.fastq.gz` | `1132 LPPC EE R-S.dna` | Yes |
| `1133.fastq.gz` | `1133 LPPC EE R-T.fa` | Yes |
| `1132.fastq.gz` | `11320 plasmid.fa` | No |
| `9M53XY_1_1108-G1.fastq` | `1108 plasmid.dna` | Yes |

For `9M53XY_1_1108-G1.fastq`, `9M53XY` is the sequencing-run token, `1` is the
barcode, and **`1108` is the sample ID**. Numeric or `barcode01`-style barcode
fields are supported. A sample ID ends at a separator; digits inside an
alphanumeric run token are not used as IDs. Multiple run/barcode files with ID
`1108` automatically share its reference and combined IGV session. Output prefixes
and track labels retain the complete FASTQ stem to distinguish those read files.

Numbers are compared as written (`001132` and `1132` are different IDs). Every
FASTQ gets its own outputs, including files sharing a sample ID. Missing
references or filenames without a recognized sample ID still get QC. Multiple matching
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
  1132_reference.fasta.fai
  1132_sorted.bam
  1132_sorted.bam.bai
  1132_alignment.log
  1132_igv.xml
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

Batches also create an **IGV desktop session** for every successfully aligned
sample. For `1132.fastq.gz`, open `1132_igv.xml` in IGV with **File → Open Session**.
The session loads the exact generated reference FASTA used for alignment and
the sorted BAM, with its BAM index. It opens at the first reference sequence.
There is no need to separately choose a genome, load a BAM, or save an initial
session. You can adjust the view in IGV and use **File → Save Session** afterward.

### IGV sessions for existing results

To create sessions for BAMs you already generated, without rerunning QC or alignment:

```bash
python batch_fastq_reader.py "/path/to/fastq_folder" \
  --output-dir "/path/to/fastq_folder/results" --igv-only
```

Or in a notebook cell:

```python
from batch_fastq_reader import create_igv_sessions
igv_report = create_igv_sessions(OUTPUT_DIR)
display(igv_report)
```

This scans the output folder for `*_sorted.bam`, uses each corresponding
`*_reference.fasta`, checks that the `.bam.bai` exists, creates a FASTA index
with `samtools faidx` if needed, and writes `*_igv.xml`. `igv_session_report.csv`
lists successes and failures. Existing nonempty sessions are preserved unless
`--overwrite` (or `overwrite=True`) is supplied.

Sessions use paths relative to the XML file, so keep the XML, reference FASTA,
`.fasta.fai`, BAM, and `.bam.bai` together when moving/downloading the results
folder. IGV must have access to those files on the machine where it runs. Session
generation does not require IGV to be installed or running and does not click or
open the desktop application automatically. Missing-reference QC samples do not
get an IGV session.

Completed older batch outputs are upgraded with sessions on the next matching
batch run, without rerunning alignment; deleted sessions or FASTA indexes are
also regenerated. Existing customized sessions are preserved on skipped runs.

### One shared reference, multiple FASTQs

**Automatic matching is the default:** FASTQs sharing the same detected sample
ID use the reference with that ID. For example, `1183_2-3.fastq.gz`,
`1183_4-6.fastq.gz`, `1183_12-17.fastq.gz`, and `1183_19-21.fastq.gz` all use
`1183 ref.dna` (or an `1183` FASTA). The suffix ranges label read files and are
ignored when choosing the reference. `11830` is a different ID and will not
match `1183`. Multiple matching references remain ambiguous.

Put the reference in the FASTQ folder or select its folder with `REFERENCE_DIR` /
`--reference-dir`. Leave `REFERENCE_FILE = None`; there is no need to select a
shared file. The batch reports the detected shared-reference groups and creates
one `1183 ref_all_reads_igv.xml` with four separately labeled BAM tracks, in
addition to the individual sessions. Multiple sample IDs in one folder are
grouped independently.

**Optional override for different IDs:**
Use `--reference-file` to align every FASTQ to a single reference instead of
matching references by sample ID:

```bash
python batch_fastq_reader.py "/path/to/fastq_folder" \
  --reference-file "/path/to/plasmid.dna" \
  --output-dir "/path/to/fastq_folder/results"
```

In the notebook, set `REFERENCE_FILE = Path("/path/to/plasmid.dna")`; FASTA files
also work. `REFERENCE_FILE` overrides `REFERENCE_DIR` and ID-based matching.
Each FASTQ still produces its own QC outputs and BAM. The shared reference uses
the same sequence identifiers across BAMs, including when converted from SnapGene.

The code creates **`plasmid_all_reads_igv.xml`**, containing every BAM aligned to
that reference as a separately labeled track (`1132`, `1133`, etc.), as well as
the individual sample sessions. Open the combined session in IGV to compare the
samples together. Combined sessions are also created when normal folder matching
maps several FASTQs to an identical generated reference, and by `--igv-only` for
compatible existing results. References with different sequence IDs/content are
kept separate so BAMs are not displayed against an incompatible genome.

Combined sessions are regenerated to include the current batch's successful
alignments on each run. Save customized overview sessions under another filename
to retain your layout. `combined_igv_session` in the batch/session report points
to the relevant combined file. If older SnapGene BAMs used differing sample IDs
as their reference sequence names, rebuild with `--reference-file` and
`--overwrite` to get compatible contig names before combining them.

### Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover compressed/plain FASTQ QC, leading-number matching, top-level-only
scanning, missing/ambiguous references, failures, output collisions, resume and
overwrite, and binary SnapGene extraction plus real alignment. The alignment
test is skipped if minimap2 or samtools is absent.

