#!/bin/bash
#
# Usage:
#   gendist_tree.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate bootstrapped genetic distance matrices with ngsDist and generate
#   a neighbor-joining tree with FastME.
#
# Requirements:
#   'gendist_tree.sh' and 'make_pops_label.py' must be in the same directory
#   fastme
#   iqtree
#   ngsdist
#   python3
#
# Note:
#   The input directory here should be the output directory of 'snp_calling_gendist.sh'.

# find MAF file generated with ANGSD '-doMaf'
maf=$(find $1 -name '*.mafs.gz')
# count the number of SNP sites
n_sites=$(zcat ${maf} | tail -n +2 | wc -l)
# find genotype likelihoods generated with ANGSD '-doGlf 2'
gl=$(find $1 -name '*.beagle.gz')
# number of samples
n_ind=$(wc -l $1/bamlist)

# create a file with population labels for each sample
python3 ./make_pops_label.py $1/bamlist $2/pops.label

# compute pairwise genetic distances with bootstrap replicates
ngsDist \
  --geno ${gl} \
  --probs \
  --n_ind ${n_ind} \
  --n_sites ${n_sites} \
  --labels $2/pops.label \
  --n_boot_rep 1000 \
  --boot_block_size 500 \
  --out $2/all.boot.dist \
  --n_threads $3 \
  &> $2/ngsdist.log

# infer trees based on bootstrapped distance matrices
fastme -i $2/all.boot.dist -o $2/all.boot.tree -D 1001 -s -T $3

# split the main tree from the bootstraped ones
head -n 1 $2/all.boot.tree > $2/main.tree
tail -n +2 $2/all.boot.tree | awk 'NF' > $2/boot.tree

# assign branch support values onto the main tree
iqtree -sup $2/main.tree $2/boot.tree
