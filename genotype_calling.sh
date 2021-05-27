#!/bin/bash
#
# Usage:
#   genotype_calling.sh <ref.fa> <bams_dir> <depth.stats> <sites_file> <outdir>
#
# Description:
#   Call genotypes with bcftools using the SNPs called with ANGSD and convert
#   BCF to Oxford GEN format.
#
# Requirements:
#   bcftools
#
# Important:
#   <ref.fa> must be the masked Kordofan giraffe assembly with scaffolds < 1 Mbp
#   removed, i.e. the 'kordofan-giraffe_1mb_scaffolds.masked.fa' file generated
#   with 'process_assembly.sh'.
#   <depth.stats> is the 'giraffe.depth.stats' file generated with 'snp_calling.sh'.
#   <sites_file> is the file containing the SNPs called by ANGSD.

# create an array of population IDs
pops=(WA GNP ZNP SNR ETH MF RET ISC LVNP MA SGR MTNP BNP SUN KKR V23 ENP HNB)
# create a list of BAM files sorted by population
for pop in ${pops[@]}; do
  ls -1 -v $2/${pop}*.clean.bam >> $5/bamlist
done

# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $3)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3)
# if minimum depth threshold is negative make it 1
if [[ ${min_dp} -lt 0 ]]; then min_dp=1; fi

# joint genotype calling
bcftools mpileup -b $5/bamlist -C 50 -f $1 -T $4 -Ou \
  | bcftools call -cv -Ou - \
  | bcftools filter -e "DP<${min_dp} || DP>${max_dp} || MQ<30 || QUAL<30 || F_MISSING>0.1" -s 'FAIL' -Ou - \
  | bcftools view -v snps -m 2 -M 2 -c 1:minor -i 'FILTER="PASS"' -Ob -o $5/snps.filtered.bcf -

# index BCF
bcftools index $5/snps.filtered.bcf

# convert BCF to Oxford GEN format
bcftools convert -t ^chrX,chrY -g $5/snps.filtered --chrom --tag PL $5/snps.filtered.bcf
