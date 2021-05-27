#!/bin/bash
#
# Usage:
#   assembly_qc.sh <asm_fasta> <ref_fasta> <ncpus>
#
# Description:
#   Calculate assembly statistics with QUAST, assess assembly completeness with
#   BUSCO, and generate sinteny plot with JupiterPlot.
#
# Requirements:
#   busco
#   jupiterplot
#   quast

# calculate assembly statistics
quast \
  -o $(dirname $1)/quast_results \
  --min-contig 1 \
  --threads $3 \
  --split-scaffolds \
  --labels 'Kordofan giraffe' \
  --large \
  --contig-thresholds 1000,10000,100000,1000000,10000000 \
  $1

# assess assembly completeness
run_BUSCO.py \
  --input $1 \
  --output $(dirname $1)/busco_results \
  --mode genome \
  --lineage_dataset mammalia_odb9 \
  --long \
  --cpu $3

# generate sinteny plot by aligning the assembly to a reference genome
jupiter name='PLA01_vs_Agaba-et-al-2016' ref=$2 fa=$1 ng=95 m=100000 t=$3
