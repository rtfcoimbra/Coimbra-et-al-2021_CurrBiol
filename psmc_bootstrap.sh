#!/bin/bash
#
# Usage:
#   psmc_bootstrap.sh <ref.fa> <bams_dir> <outdir> <ncpus>
#
# Description:
#   Calculate per-base depth from cleaned BAMs with sambamba, generate consensus
#   genome sequences in FASTQ format with BCFtools, and infer giraffe effective
#   population size changes through time with the PSMC model.
#
# Requirements:
#   bcftools
#   parallel
#   perl
#   psmc
#   python3
#   sambamba
#   'site_depth_stats.py' (must be in same directory as 'psmc.sh')
#
# Important:
#   <ref.fa> must be the masked Kordofan giraffe assembly with scaffolds < 1 Mbp
#   removed, i.e. the 'kordofan-giraffe_1mb_scaffolds.masked.fa' file generated
#   with 'process_assembly.sh'.

# create an array of sample names
samples=(
  WA720 WA806 WA808
  GNP04 GNP05 ZNP01
  ETH3 MF22 MF24
  ISC04 ISC08 RET4
  MA1 SGR14 SGR05
  KKR01 SUN3 V23
  ENP16 ENP19 HNB110
)

# iterate over samples in the array
for sample in ${samples[@]}; do
  # calculate site depth per sample
  sambamba depth base -t $4 $2/${sample}.clean.bam \
    | awk 'NR>1 { print $3 }' > $3/${sample}.depth
done

# calculate site depth statistics
parallel -j $4 "python3 site_depth_stats.py {} > {}.stats" ::: $3/*.depth

# iterate over samples in the array
for sample in ${samples[@]}; do
  # extract the median of the read depth distribution
  median_dp=$(grep -Po 'Median: \K\d+' $3/${sample}.depth.stats)
  # calculate maximum depth filter
  max_dp=$(python -c "print(f'{int(${median_dp}) * 2}')")
  # generate consensus genome sequence in FASTQ format
  echo "bcftools mpileup -C 50 -f $1 -Ou $2/${sample}.clean.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D ${max_dp} -Q 30 | gzip > $3/${sample}.consensus.fq.gz" >> $3/consensus_call.jobs
  # transform consensus sequence into a FASTA-like format where the i-th
  # character in the output sequence indicates whether there is at least one
  # heterozygote in the bin
  echo "fq2psmcfa -q 30 $3/${sample}.consensus.fq.gz > $3/${sample}.psmcfa" >> $3/fq2psmcfa.jobs
  # split long scaffold sequences to shorter segments
  echo "splitfa $3/${sample}.psmcfa 10000 > $3/${sample}.split.psmcfa" >> $3/splitfa.jobs
  # infer population size history
  echo "psmc -N 25 -t 15 -r 5 -p '4+25*2+4+6' -o $3/${sample}.psmc $3/${sample}.psmcfa" >> $3/psmc.jobs
  # perform 100 rounds of bootstrapping
  for round in {1..100}; do
    echo "psmc -N 25 -t 15 -r 5 -b -p '4+25*2+4+6' -o $3/${sample}.round-$(printf '%03d' ${round}).psmc $3/${sample}.split.psmcfa" >> $3/psmc_boot.jobs
  done
done

# run jobs in parallel
cat $3/consensus_call.jobs | parallel -j $4 &&
cat $3/fq2psmcfa.jobs | parallel -j $4 &&
cat $3/splitfa.jobs | parallel -j $4 &&
cat $3/psmc.jobs | parallel -j $4 &&
cat $3/psmc_boot.jobs | parallel -j $4

# multiline plot for one sample per subspecies
psmc_plot.pl \
  -u 2.12e-08 \
  -g 10 \
  -X 1e07 \
  -Y 4 \
  -M 'West African','Kordofan','Nubian','Reticulated','Masai s. str.','South African','Angolan' \
  -P 'right bottom' \
  -p \
  $3/multiline_plot \
  $3/{WA720,ZNP01,MF22,ISC08,SGR14,KKR01,ENP19}.psmc

# iterate over samples in the array
for sample in ${samples[@]}; do
  # concatenate bootstrapped psmc outputs
  cat $3/${sample}.psmc $3/${sample}.round-*.psmc > $3/${sample}.boot.psmc
  # plot bootstrapped individual results
  psmc_plot.pl -u 2.12e-08 -g 10 -X 1e07 -Y 4 -T "${sample}" -p $3/${sample}.boot $3/${sample}.boot.psmc
done

# organize directory
mkdir $3/psmc_boot; mv $3/*.round-*.psmc $3/psmc_boot
