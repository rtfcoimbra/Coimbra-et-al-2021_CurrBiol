#!/bin/bash
#
# Usage:
#   mapping_qc.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate mapping statistics and generate plots of various mapping quality
#   measurements with QualiMap.
#
# Requirements:
#   qualimap

# iterate over duplicate marked BAM files
for dedup_bam in $1/*.dedup.bam; do
  # output directory
  outdir=$2/$(basename ${dedup_bam//bam/qualimap})
  # run QualiMap
  qualimap bamqc -bam ${dedup_bam} -ip -nt $3 -outdir ${outdir} --java-mem-size=16G &> ${outdir}.log
done
