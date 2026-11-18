# Local de novo plasmid assembly

This branch adds `denovo_assembly.py` and a plain-Python PyCharm launcher,
`run_denovo_assembly.py`. The existing alignment/QC workflow is still available.

Flye assembles the reads independently and polishes its contigs using the reads.
An optional SnapGene/FASTA reference is compared **after** assembly; it never
supplies bases, contig ordering, gap fills, or a forced plasmid sequence.
All assembled contigs are retained in the output FASTA.

## Install once

Use Linux, macOS, or WSL Ubuntu on Windows. Flye is not a native Windows program.
One option is a dedicated Conda/Miniforge environment:

```bash
conda create -n plasmid-assembly -c conda-forge -c bioconda --strict-channel-priority \
  python=3.12 flye=2.9.6 minimap2 samtools pandas matplotlib biopython
conda activate plasmid-assembly
```

Availability of Bioconda binaries depends on your operating system/architecture.
If Flye is unavailable there, follow its supported
[source installation instructions](https://github.com/fenderglass/Flye/blob/2.9.6/docs/INSTALL.md).
The source build requires a C++ compiler, make, and zlib development files. The
tested source version is Flye 2.9.6. `pip install -r requirements.txt` installs
this project's Python dependencies, **not** Flye/minimap2/samtools.

Verify these commands in the same environment that runs Python:

```bash
flye --version
minimap2 --version
samtools --version-only
```

In PyCharm, select this environment's Python interpreter. If tools are missing
from the IDE's PATH, set `TOOL_DIRS` in the launcher to their `bin` directories.
Keep the downloaded branch files together in one project folder. No data needs
to be uploaded to GitHub or the cloud environment.

## PyCharm: edit settings and click Run

Open `run_denovo_assembly.py` and set:

```python
FASTQ_FILE = Path("/path/to/9M53XY_1_1108-G1.fastq.gz")
REFERENCE_FILE = Path("/path/to/1108_reference.dna")  # optional; None skips comparison
OUTPUT_DIR = None  # creates 9M53XY_1_1108-G1_denovo beside the FASTQ
READ_TYPE = "nano-hq"
THREADS = 4
RESUME = False
```

Then run the entire file. `.fastq`, `.fq`, and gzipped variants are supported.
For Windows paths, use raw strings, e.g. `Path(r"C:\Users\you\data\reads.fastq.gz")`;
when running in WSL, use Linux paths such as `/mnt/c/Users/you/data/reads.fastq.gz`.
WSL network shares need to be mounted separately.

`nano-hq` is intended for modern ONT R10/Q20 reads. Use `nano-raw` for older/noisier
ONT reads; do not choose the mode purely from the file extension. The read QC
report includes average encoded Phred score, which is descriptive rather than a
calibration guarantee. Other supported modes include `pacbio-hifi` and `pacbio-corr`.
`GENOME_SIZE` is an optional estimate such as `"105k"`; it does not constrain the
assembly to a reference. `MIN_OVERLAP = None` uses Flye's default. For short
plasmids/reads you may need `MIN_OVERLAP = 1000`; reads that do not overlap enough
cannot support a complete assembly. Insufficient coverage or a fragmented result
is reported rather than replaced with reference sequence.

## Command line

```bash
python denovo_assembly.py "/path/to/reads.fastq.gz" \
  --reference-file "/path/to/reference.dna" \
  --output-dir "/path/to/new_assembly_folder" \
  --read-type nano-hq --threads 4
```

Omit `--reference-file` to assemble without any reference. The workflow runs one
selected FASTQ at a time; filenames do not need matching IDs. Multiple FASTQ
files are not automatically pooled, since they may represent different plasmids.

## Outputs and interpretation

- `<FASTQ stem>_assembly.fasta`: the Flye assembly, with **all contigs** retained.
- `REPORT.md` / `summary.json`: contig lengths, assembler circularity flags,
  ambiguous-base counts, read-back depth/coverage, notes and tool versions.
- `read_qc.json`: number/bases of reads, median length, N50, longest read, mean Phred.
- `<stem>_reads_to_assembly.bam` and `.bai`: reads mapped to the assembly.
- `<stem>_assembly_igv.xml`: open in IGV with **File → Open Session** to inspect the
  reads against the **new assembly**, not the original reference.
- `read_back_coverage.tsv`: samtools per-contig coverage summary (not per-base accuracy).
- With a reference: `comparison_reference.fasta`, `reference_comparison.csv`,
  `assembly_vs_reference.paf` (including CIGAR and `cs` difference tags).
- `flye/`: full original Flye outputs, assembly graph, contig table, polishing data.
- `run_state.json` and logs: input fingerprints, settings, executed commands,
  tool versions, completion or failure information.

**An assembly is not automatically a verified complete plasmid.** Multiple contigs
are not merged artificially. Flye's circular flag is a graph classification;
read support and repeat structure still need inspection. Read-back coverage alone
does not exclude mixtures or establish a consensus quality score. Flye polishing
does not guarantee every homopolymer or low-support base is correct.

Reference comparison rows are **local alignment blocks**, not a final variant
list. Split alignments can arise from circular starting-position differences or
repeats. Do not sum overlapping blocks as unique coverage or infer a deletion
from an unaligned region without checking the assembly graph/read support.

Keep the assembly FASTA, its `.fai`, BAM, `.bai`, and XML together when moving
results so the relative paths in IGV continue to work.

## Existing runs and failures

The output directory must be new/empty for a new run. Existing input/output files
are never deleted automatically. For a failed or completed run with identical
inputs/settings, set `RESUME = True` or add `--resume`. A completed run is reused
when its expected artifacts still exist; an incomplete Flye run resumes using
Flye's own checkpoints. Inputs are checked by path, size and modification time,
not full content hashes. If settings or inputs change, select a new output folder.
Inspect `assembly.log`, `read_back_alignment.log`, and `run_state.json` for failures.

## Validation

```bash
python -m unittest discover -s tests -v
```

For the real end-to-end assembly integration test, install the tools and set
`FLYE_TEST_DATA` to Flye's bundled `flye/tests/data` directory. It uses the bundled
500 kb read dataset, validates FASTA/BAM/IGV output and optional SnapGene comparison,
and checks that the reference is absent from the assembler command. Without those
fixtures/tools that integration test is explicitly skipped. No uploaded user reads
or references are included in this branch.
