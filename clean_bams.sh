#!/bin/bash
#
# Usage:
#   clean_bams.sh <regions.bed> <indir> <outdir> <njobs>
#
# Description:
#   Remove unmapped (4), secondary (256), QC failed (512), duplicate (1024), and
#   supplementary (2048) reads from indel-realigned BAMs, and keep only reads
#   mapped in a proper pair (2) to regions in a BED file (non-repetitive regions
#   in scaffolds >= 1 Mb) using SAMtools.
#
# Requirements:
#   samtools
#   parallel
#   'no_repeats_1mb_scaffolds.bed' file generated with 'process_assembly.sh'
#
# Important:
#   Each job uses 4 CPU threads.

# clean indel-realigned BAMs
parallel -j $4 --plus \
  "samtools view -@ 4 -b -F 3844 -f 2 -L $1 -o $3/{/..}.clean.bam {}" \
  ::: $2/*.realigned.bam

# index clean BAMs
parallel -j $4 samtools index -b {} ::: $3/*.clean.bam
