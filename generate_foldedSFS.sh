#!/bin/bash
#
# Usage:
#   generate_foldedSFS.sh <ref.fa> <indir_bam> <indir_depth> <outdir> <njobs>
#
# Description:
#   Generate the folded site frequency espectrum per sample with ANGSD.
#
# Requirements:
#   angsd (incl. the realSFS companion program)
#   parallel
#
# Important:
#   Each job uses up to 4 CPU threads.
#   <ref.fa> must be the masked Kordofan giraffe assembly with scaffolds < 1 Mbp
#   removed, i.e. the 'kordofan-giraffe_1mb_scaffolds.masked.fa' file generated
#   with 'process_assembly.sh'.
#   <indir_depth> is the directory containing the '*.depth.stats' files for each
#   sample previously generated with 'generate_fasta.sh'.

for bam in $2/*.clean.bam; do
  # sample name
  sample=$(basename ${bam%.clean.bam})
  # get the 95th percentile of the site depth distribution
  max_dp=$(grep -Po '95th percentile: \K\d+' $3/${sample}.depth.stats)
  # set filters
  filters="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 30 -setMaxDepth ${max_dp}"
  # set tasks
  todo="-doCounts 1 -doSaf 1"
  # estimate the folded sample allele frequency (SAF) likelihoods for each site
  echo "angsd -i ${bam} -ref $1 -anc $1 ${filters} ${todo} -fold 1 -GL 1 -P 4 -out $4/${sample} &> $4/${sample}.angsd.log" >> $4/angsd_het.jobs
  # estimate the folded site frequency spectrum (SFS)
  echo "realSFS $4/${sample}.saf.idx -bootstrap 200 > $4/${sample}.sfs 2> $4/${sample}.realSFS.log" >> $4/realSFS_het.jobs
done

# run ANGSD in parallel
cat $4/angsd_het.jobs | parallel -j $5
# run realSFS in parallel
cat $4/realSFS_het.jobs | parallel -j $5
