#!/bin/bash
#
# Usage:
#   annotate_coding.sh <assembly.fa> <repeats.bed> <proteins.faa> <ncpus>
#
# Description:
#   Soft-mask assembly and perform the automated prediction and annotation of
#   protein-coding genes with BRAKER. Protein sequences derived from cow (Bos
#   taurus; GCA_002263795.2) are used as extrinsic evidence for homolog-based
#   predictions.
#
# Requirements:
#   bedtools
#   braker (and its many dependencies)
#
# Important:
#   <repeats.bed> is the 'repeats_merged_adj.bed' file generated with
#   'process_assembly.sh'.
#   <proteins.faa> is the 'GCF_002263795.1_ARS-UCD1.2_protein.faa' file
#   downloaded from NCBI.

# soft-mask assembly for annotation of coding regions
bedtools maskfasta -soft -fi $1 -bed $2 -fo ${1/.fa/.softmasked.fa}

# perform automated prediction and annotation of protein-coding genes
braker.pl \
  --cores $4 \
  --species=giraffe \
  --genome=${1/.fa/.softmasked.fa} \
  --prot_seq=$3 \
  --softmasking \
  --UTR=off \
  --prg=gth \
  --ALIGNMENT_TOOL_PATH=/opt/gth-1.7.1-Linux_x86_64-64bit/bin \
  --AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
  --AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
  --AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
  --BLAST_PATH=/opt/ncbi-blast-2.2.31+/bin
