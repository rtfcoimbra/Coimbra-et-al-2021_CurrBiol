#!/bin/bash
#
# Usage:
#   pcangsd_hwe.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate a covariance matrix for SNPs filtering out those that strongly
#   deviate from the Hardy-Weinberg equilibrium using PCAngsd.
#
# Requirements:
#   pcangsd
#   python3
#
# Note:
#   The input directory here should be the output directory of 'ld_pruning.sh'.
#   Adjust the path to PCAngsd in the variable ${pcangsd}.

# path to PCAngsd
pcangsd=~/software/pcangsd/pcangsd.py

# iterate over genotype likelihood files
for gl in $1/*.beagle.gz; do
  # output file directory and name
  out=$2/$(basename ${gl%.beagle.gz})

  # calculate covariance matrix, estimate per-site inbreeding coefficients, and
  # perform likehood ratio tests for HWE
  python ${pcangsd} \
    -beagle ${gl} \
    -minMaf 0.05 \
    -o ${out} \
    -inbreedSites \
    -threads $3 \
    &> ${out}.pcangsd.log

  # add delimiter between logs to 'pcangsd.log'
  echo -e "\n---\n" >> ${out}.pcangsd.log

  # calculate covariance matrix with HWE filter
  python ${pcangsd} \
    -beagle ${gl} \
    -minMaf 0.05 \
    -o ${out}.hwe_filter \
    -hwe ${out}.lrt.sites.npy \
    -hwe_tole 1e-6 \
    -threads $3 \
    &>> ${out}.pcangsd.log
done
