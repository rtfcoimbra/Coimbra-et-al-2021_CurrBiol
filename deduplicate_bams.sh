#!/bin/bash
#
# Usage:
#   deduplicate_bams.sh <indir> <outdir> <njobs>
#
# Description:
#   Mark PCR/optical duplicate reads with Picard's MarkDuplicates.
#
# Requirements:
#   BAMs must be named 'sample.sorted.bam'
#   parallel
#   picard
#   samtools
#
# Important:
#   Each MarkDuplicates' job uses up to 2 CPU threads.
#   Adjust the path to Picard in the variable ${picard}.

# path to Picard
picard=~/software/picard.jar

# create list of Picard's MarkDuplicates commands for GNU parallel
for sorted_bam in $1/*.sorted.bam; do
  # output filenames
  dedup_bam=$2/$(basename ${sorted_bam/sorted/dedup})
  metrics=$2/${dedup_bam/bam/metrics.txt}
  log=$2/${dedup_bam/dedup.bam/markduplicates.log}

  # run MarkDuplicates without optical duplicate detection
  # samples MA1 and WOAK do not have proper read headers
  if [[ ${sorted_bam} == $1/@(MA1|WOAK).sorted.bam ]]; then
    echo \
      "java -XX:ParallelGCThreads=2 -Xmx8G -jar ${picard} \
        MarkDuplicates \
          MAX_FILE_HANDLES=2048 \
          I=${sorted_bam} \
          O=${dedup_bam} \
          M=${metrics} \
          READ_NAME_REGEX='null' \
          &> ${log}" >> $2/markduplicates.jobs

  # run MarkDuplicates with optical duplicate pixel distance recommended for
  # reads sequenced in unpatterned flowcells
  elif [[ ${sorted_bam} == $1/ZNP01.sorted.bam ]]; then
    echo \
      "java -XX:ParallelGCThreads=2 -Xmx8G -jar ${picard} \
        MarkDuplicates \
          MAX_FILE_HANDLES=2048 \
          I=${sorted_bam} \
          O=${dedup_bam} \
          M=${metrics} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
          &> ${log}" >> $2/markduplicates.jobs

  # run MarkDuplicates with optical duplicate pixel distance recommended for
  # reads sequenced in patterned flowcells
  else
    echo \
      "java -XX:ParallelGCThreads=2 -Xmx8G -jar ${picard} \
        MarkDuplicates \
          MAX_FILE_HANDLES=2048 \
          I=${sorted_bam} \
          O=${dedup_bam} \
          M=${metrics} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          &> ${log}" >> $2/markduplicates.jobs
  fi
done

# run picard in parallel
cat $2/markduplicates.jobs | parallel -j $3

# index new BAMs
parallel -j $3 samtools index -b {} ::: $2/*.dedup.bam
