#!/usr/bin/env bash
set -euo pipefail

echo "=== Nanopore → plasmid alignment helper ==="
echo

#############################################
# Helper: resolve filename with or without extension
#############################################
resolve_file() {
  local folder="$1"
  local name="$2"
  local kind="$3"   # "ref" or "ont"

  # Strip surrounding quotes
  name="${name%\"}"
  name="${name#\"}"
  name="${name%\'}"
  name="${name#\'}"

  # 1) Try as-is
  if [ -f "${folder}/${name}" ]; then
    echo "${folder}/${name}"
    return 0
  fi

  # 2) Normalize base name (strip common extensions if user half-typed)
  local base="$name"
  base="${base%.gz}"
  base="${base%.fa}"
  base="${base%.fasta}"
  base="${base%.fna}"
  base="${base%.fas}"
  base="${base%.fastq}"
  base="${base%.fq}"

  local exts=()
  if [ "$kind" = "ref" ]; then
    exts=(.fa .fasta .fna .fas)
  else
    exts=(.fastq .fq .fastq.gz .fq.gz)
  fi

  local cand
  for ext in "${exts[@]}"; do
    cand="${folder}/${base}${ext}"
    if [ -f "$cand" ]; then
      echo "$cand"
      return 0
    fi
  done

  # If nothing found, return empty string + non-zero status
  echo ""
  return 1
}

#############################################
# 1) Ask for folder
#############################################

read -r -p "Paste the folder where ref and ONT files are stored: " seq_folder
if [ -z "$seq_folder" ]; then
  echo "ERROR: Folder path cannot be empty." >&2
  exit 1
fi

# Strip surrounding single/double quotes if present (Finder “Copy as Pathname”)
seq_folder="${seq_folder%\"}"
seq_folder="${seq_folder#\"}"
seq_folder="${seq_folder%\'}"
seq_folder="${seq_folder#\'}"

# Expand leading ~ to $HOME
seq_folder="${seq_folder/#\~/$HOME}"

# Remove trailing slash, if any
seq_folder="${seq_folder%/}"

if [ ! -d "$seq_folder" ]; then
  echo "ERROR: Folder does not exist: $seq_folder" >&2
  exit 1
fi

echo
echo "Folder set to: $seq_folder"
echo "Files in folder (preview):"
ls "$seq_folder" || true
echo

#############################################
# 2) Ask for filenames (with or without extension)
#############################################

read -r -p "Paste the reference name (with or without extension, e.g. plasmid or plasmid.fa): " ref_input
if [ -z "$ref_input" ]; then
  echo "ERROR: Reference name cannot be empty." >&2
  exit 1
fi

read -r -p "Paste the ONT reads name (with or without extension, e.g. sample1 or sample1.fastq.gz): " ont_input
if [ -z "$ont_input" ]; then
  echo "ERROR: ONT name cannot be empty." >&2
  exit 1
fi

# Resolve to actual files on disk
ref_path="$(resolve_file "$seq_folder" "$ref_input" "ref" || true)"
ont_path="$(resolve_file "$seq_folder" "$ont_input" "ont" || true)"

if [ -z "$ref_path" ]; then
  echo "ERROR: Could not find reference file for input '$ref_input' in folder: $seq_folder" >&2
  echo "       Tried common FASTA extensions (.fa, .fasta, .fna, .fas)." >&2
  exit 1
fi

if [ -z "$ont_path" ]; then
  echo "ERROR: Could not find ONT file for input '$ont_input' in folder: $seq_folder" >&2
  echo "       Tried common FASTQ extensions (.fastq, .fq, .fastq.gz, .fq.gz)." >&2
  exit 1
fi

#############################################
# 3) Check dependencies
#############################################

if ! command -v minimap2 >/dev/null 2>&1; then
  echo "ERROR: minimap2 not found in PATH. Install via e.g. 'brew install minimap2'." >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH. Install via e.g. 'brew install samtools'." >&2
  exit 1
fi

#############################################
# 4) Build output names (fastaRoot_vs_fastqRoot)
#############################################

ref_base=$(basename "$ref_path")
ont_base=$(basename "$ont_path")

# Strip possible extensions
ref_root="$ref_base"
ref_root="${ref_root%.fa}"
ref_root="${ref_root%.fasta}"
ref_root="${ref_root%.fna}"
ref_root="${ref_root%.fas}"

ont_root="$ont_base"
ont_root="${ont_root%.gz}"
ont_root="${ont_root%.fastq}"
ont_root="${ont_root%.fq}"

out_sam="${seq_folder}/${ref_root}_vs_${ont_root}.sam"
out_bam="${seq_folder}/${ref_root}_vs_${ont_root}.sorted.bam"

echo
echo "Resolved files:"
echo "  Reference: $ref_path"
echo "  Reads:     $ont_path"
echo
echo "Planned outputs:"
echo "  SAM : $out_sam"
echo "  BAM : $out_bam"
echo

read -r -p "Proceed with alignment? [y/N]: " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
  echo "Aborted by user."
  exit 0
fi

#############################################
# 5) Run minimap2 + samtools
#############################################

echo
echo "Running minimap2 + samtools..."
echo

minimap2 -ax map-ont --secondary=no "$ref_path" "$ont_path" \
  | tee "$out_sam" \
  | samtools view -bS - \
  | samtools sort -o "$out_bam"

samtools index "$out_bam"

echo
echo "Done."
echo "  SAM : $out_sam"
echo "  BAM : $out_bam"
echo "  BAI : ${out_bam}.bai"
echo
