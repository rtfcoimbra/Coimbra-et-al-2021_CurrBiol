#!/bin/bash
#
# Usage:
#   merge_bams.sh <indir> <outdir> <njobs>
#
# Description:
#   Merge same sample BAMs from different lanes with SAMtools.
#
# Requirements:
#   BAMs must be named 'sample_L?.sorted.bam'
#   samtools
#   parallel

# set shell nullglob
shopt -s nullglob

# create an array of BAMs to merge
lane_bams=($1/*_L*.sorted.bam)

# unset shell nullglob
shopt -u nullglob

# iterate over elements in the array
for ((i=0; i<${#lane_bams[@]}; ++i)); do
  # output name
  sample_bam=${lane_bams[i]/_L*./.}
  # create list of SAMtools commands for GNU parallel
  echo "samtools merge ${sample_bam} ${lane_bams[i]} ${lane_bams[((++i))]}" >> $2/samtools-merge.jobs
done

# run SAMtools in parallel
cat $2/samtools-merge.jobs | parallel -j $3
