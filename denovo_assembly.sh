#!/bin/bash
#
# Usage:
#   denovo_assembly.sh <fastqs_dir> <outdir> <ncpus>
#
# Description:
#   Generates a phased whole-genome de novo assembly from a Chromium-prepared
#   library using Supernova.
#
# Requirements:
#   supernova
#
# Note:
#   FASTQ files must be named 'PLA01_S1_L003_R1_001.fastq.gz' and
#   'PLA01_S1_L003_R2_001.fastq.gz'

# generate whole genome de novo assembly
supernova run --id='PLA01' --fastqs=$1 --maxreads=all --localcores=$3

# generate pseudohaploid assembly FASTA
supernova mkoutput \
  --style=pseudohap \
  --asmdir=$2/PLA01/outs/assembly \
  --outprefix=$2/PLA01/outs/PLA01_assembly_pseudohap
