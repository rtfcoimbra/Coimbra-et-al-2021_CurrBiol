#!/bin/bash
#
# Usage:
#   snp_calling.sh <ref.fa> <indir> <outdir>
#
# Description:
#   Call and filter SNPs with ANGSD. Uses giraffe BAMs.
#
# Requirements:
#   angsd
#   python3
#   sambamba
#   'site_depth_stats.py' (must be in same directory as 'snp_calling.sh')
#
# Important:
#   <ref.fa> must be the masked Kordofan giraffe assembly with scaffolds < 1 Mbp
#   removed, i.e. the 'kordofan-giraffe_1mb_scaffolds.masked.fa' file generated
#   with 'process_assembly.sh'.

# find giraffe BAMs
bams=$(find $2 -name '*.clean.bam' ! -name 'WOAK.clean.bam' -printf '%p ')

# calculate global site depth
sambamba depth base -t 8 --combined ${bams} | awk 'NR>1 { print $3 }' > $3/giraffe.depth

# calculate site depth statistics
python3 ./site_depth_stats.py $3/giraffe.depth > $3/giraffe.depth.stats

# create an array of population names
pops=(WA GNP ZNP SNR ETH MF RET ISC LVNP MA SGR MTNP BNP SUN KKR V23 ENP HNB)
# create list of input BAMs sorted by population
for pop in ${pops[@]}; do
  ls -1 $2/${pop}*.clean.bam >> $3/bamlist
done

# count the number of individuals
n_ind=$(cat $3/bamlist | wc -l)
# set minimum number of individuals per site
min_ind=$(python -c "print(f'{round(${n_ind} * 0.9)}')")
# set minimum and maximum depth cutoff
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K-*\d+' $3/giraffe.depth.stats)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3/giraffe.depth.stats)

# SNP calling with filters
angsd -b $3/bamlist -ref $1 -GL 1 -P 4 -out $3/snps \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1 \
  -minMapQ 30 -minQ 30 -minInd ${min_ind} \
  -doSnpStat 1 -doHWE 1 -sb_pval 1e-6 -hwe_pval 1e-6 -hetbias_pval 1e-6 \
  -doCounts 1 -setMinDepth ${min_dp} -setMaxDepth ${max_dp} \
  -doMajorMinor 1 -skipTriallelic 1 \
  -doMaf 2 -doPost 1 -minMaf 0.05 -SNP_pval 1e-6 \
  -doGeno 8 -geno_minDepth 3 -doGlf 2 \
  &> $3/angsd.log
