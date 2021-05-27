#!/bin/bash
#
# Usage:
#   ngsadmix.sh <indir> <outdir> <k_min> <k_max> <ncpus>
#
# Description:
#   Estimate admixture proportions assuming varying K values in NGSadmix with
#   100 replicates per K, create a list of run likelihoods per K, and identify
#   the run with the highest likelihood for each K.
#
# Requirements:
#   ngsadmix

# find genotype likelihoods file
gl=$(find $1 -name '*.beagle.gz')

# set output name and directory
fout=$2/$(basename ${gl%.beagle.gz})

# iterate over a range of K values
for k in $(seq $3 $4); do
  # run 100 replicates for each K
  for run in {1..100}; do
    NGSadmix -likes ${gl} -K ${k} -minMaf 0.05 -o ${fout}.k${k}_r${run} -P $5
  done
done

# create a list of run likelihoods in descending order for each K
for k in $(seq $3 $4); do
  ( for log in $(ls ${fout}.k${k}_r*.log); do grep -Po 'like=\K[^ ]+' ${log}; done ) | sort -gr > $2/likelihoods_k${k}
done

# concatenate run likelihoods lists
cat $(ls -v $2/likelihoods_k*) > $2/likelihoods.list

# create a list with the highest likelihood run of each K
for k in $(seq $3 $4); do
  # list all log files for each K
  k_log=$(ls -v ${fout}.k${k}_r*.log)
  # find the highest likelihood run for each K
  best_run=$(head -1 $2/likelihoods_k${k} | grep -f - ${k_log} | grep -Po '_r\K[0-9]+')
  # print the highest likelihood run for each K
  echo ${fout}.k${k}_r${best_run}.qopt >> $2/best_runs.list
done
