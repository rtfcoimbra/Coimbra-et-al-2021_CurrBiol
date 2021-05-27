#!/bin/bash
#
# Usage:
#   indel_realigner.sh <ref.fa> <indir_bams> <indir_intervals> <outdir_bams> <outdir_logs> <njobs>
#
# Description:
#   Performs local realignment on reads coincident with a list of target
#   intervals using GATK's IndelRealigner.
#
# Requirements:
#   gatk 3.8-1
#   java 1.8
#   parallel
#
# Important:
#   Each job will use 4 CPUs.
#   Adjust path to GATK in the variable ${gatk}.

# path to GATK v3.8-1
gatk=~/software/GenomeAnalysisTK.jar

# iterate over deduplicated BAMs
for dedup_bam in $2/*.dedup.bam; do
  # set the species
  if [[ ${dedup_bam} == *WOAK* ]]; then
    target_intervals=$(find $3 -name 'okapi.realigner.intervals')
  else
    target_intervals=$(find $3 -name 'giraffe.realigner.intervals')
  fi
  # output BAM
  realigned_bam=$(basename ${dedup_bam//dedup/realigned})
  # logfile
  log=$(basename ${dedup_bam/dedup.bam/indel-realigner.log})
  # echo command to a file
  echo \
    "java -Xmx16G -XX:ParallelGCThreads=4 -jar ${gatk} \
      -T IndelRealigner \
      -R $1 \
      -I ${dedup_bam} \
      -targetIntervals ${target_intervals} \
      -o $4/${realigned_bam} \
      &> $5/${log}" \
    >> $5/indel-realigner.jobs
done

# run IndelRealigner in parallel
cat $5/indel-realigner.jobs | parallel -j $6

# search for warn messages
n_warns=$(grep -L 'Done. There were no warn messages.' $5/*.indel-realigner.log | wc -l)

# if there are no warnings, print message
if [[ ${n_warns} -eq 0 ]]; then
  echo 'IndelRealigner completed successfully!'
# if there are warnings, print a list of files in which they occurred
else
  echo 'These files had warn messages:'
  grep -L 'Done. There were no warn messages.' $5/*.indel-realigner.log
fi
