#!/bin/bash
#
# Usage:
#   trim_reads.sh <indir> <outdir> <njobs>
#
# Description:
#   Trim paired-end reads with Trimmomatic.
#
# Requirements:
#   reads must be named 'sample_{1,2}.fq.gz' or 'sample_lane_{1,2}.fq.gz'
#   parallel
#   trimmomatic
#
# Important:
#    Each job uses 4 CPU threads.
#    Adjust the path to the file 'TruSeq3-PE-2.fa' in the variable ${adapters}.

# trimmomatic adapters list
adapters=~/software/trimmomatic-0.38-1/adapters/TruSeq3-PE-2.fa

# iterate over read 1 files
for r1 in $1/*_1.fq.gz; do
  # read 2
  r2=${r1/_1.fq.gz/_2.fq.gz}
  # output directory and basename
  sample=$2/$(basename ${r1%_1.fq.gz})
  # echo trimmomatic commands to a file
  if [[ ${r1} = $1/ZNP01_1.fq.gz ]]; then # ZNP01 reads are phred+64
    echo "trimmomatic PE -threads 4 -phred64 ${r1} ${r2} -baseout ${sample}.fq.gz ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40 TOPHRED33 2> ${sample}.trimmomatic.log" >> $2/trimmomatic.jobs
  else # all other reads are phred+33 encoded
    echo "trimmomatic PE -threads 4 -phred33 ${r1} ${r2} -baseout ${sample}.fq.gz ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40 2> ${sample}.trimmomatic.log" >> $2/trimmomatic.jobs
  fi
done

# run trimmomatic in parallel
cat $2/trimmomatic.jobs | parallel -j $3
