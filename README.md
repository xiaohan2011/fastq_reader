A short script to read fastq files from Nanopore long-read sequencing data. Sort based on length.

Input: fastq file path, optional: fasta file path for reference.

Output: 
  if fastq only: csv file of reads sorted based on length and a histogram.
  if fastq and fasta: in addition to above, generate bam file in the same folder of fastq.

