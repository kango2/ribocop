#!/bin/bash

# Define reference genome file and directory with FASTQ files
REF_GENOME=path/to/reference/genome.fasta
FASTQ_DIR=path/to/directory/with/fastq/files

# Loop through all FASTQ files in the directory
for FASTQ in $FASTQ_DIR/*.fastq
do
  # Extract base file name without extension
  BASE=$(basename $FASTQ .fastq)
  
  # Run minimap2 alignment
  minimap2 -ax map-ont $REF_GENOME $FASTQ > $BASE.sam
  
  # Convert SAM to BAM
  samtools view -bS $BASE.sam > $BASE.bam
  
  # Sort BAM file
  samtools sort $BASE.bam -o $BASE.sorted.bam
  
  # Index BAM file
  samtools index $BASE.sorted.bam
  
  # Remove SAM file
  rm $BASE.sam
done
