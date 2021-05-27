#!/bin/bash
#
# Usage:
#   realigner_target_creator.sh <giraffe|okapi> <ref.fa> <indir> <outdir> <ndata_threads> <njobs>
#
# Description:
#   Identifies regions where alignments may potentially be improved using
#   GATK RealignerTargetCreator.
#
# Requirements:
#   gatk 3.8-1
#   java 1.8
#   parallel
#   picard
#   samtools
#
# Important:
#   Each job will use 'N data threads' + 4 CPUs ('-nt' + ParallelGCThreads).
#   Each data thread will use 20GB RAM.
#   Adjust path to Picard and GATK in the variables ${picard} and ${gatk}.

# path to Picard
picard=~/software/picard.jar

# path to GATK v3.8-1
gatk=~/software/GenomeAnalysisTK.jar

# index reference genome FASTA
samtools faidx $2

# create a sequence dictionary for the reference genome FASTA
java -jar ${picard} CreateSequenceDictionary R=$2 O=${2/.fa/.dict}

# get formatted string with input BAMs
if [[ $1 == 'giraffe' ]]; then
  bams=$(find $3 -name '*.dedup.bam' ! -name 'WOAK*' -printf '-I %p ')
elif [[ $1 == 'okapi' ]]; then
  bams=$(find $3 -name 'WOAK.dedup.bam' -printf '-I %p ')
else
  echo "Unrecognized argument: $1 ... expected 'giraffe' or 'okapi'"
  exit 1
fi

# get scaffold IDs and spans formatted as 'scaff:start-end'
cat $2.fai | awk '{ print $1":1-"$2}' \
  | while read scaff; do # iterate over scaffolds
    # echo command to a file
    echo \
      "java -Xmx20G -XX:ParallelGCThreads=4 -jar ${gatk} \
        -T RealignerTargetCreator \
        -nt $5 \
        -R $2 \
        -L ${scaff} \
        ${bams} \
        -o $4/${scaff%:*}.intervals \
        &> $4/${scaff%:*}.rtc.log" \
      >> $4/rtc.jobs
    done

# run RealignerTargetCreator in parallel
cat $4/rtc.jobs | parallel -j $6

# search for warn messages
n_warns=$(grep -L 'Done. There were no warn messages.' $4/*.rtc.log | wc -l)

# if there are no warnings, concatenate intervals and remove intermediate files
if [[ ${n_warns} -eq 0 ]]; then
  ls -v $4/*.intervals > $4/scaff.intervals.list &&
  cat $(cat $4/scaff.intervals.list) > $4/$1.realigner.intervals &&
  echo 'RealignmentTargetCreator completed successfully!' &&
  rm $4/*.{jobs,list}
# if there are warnings, print a list of files in which they occurred
else
  echo 'These files had warn messages:'
  grep -L 'Done. There were no warn messages.' $4/*.rtc.log
fi
