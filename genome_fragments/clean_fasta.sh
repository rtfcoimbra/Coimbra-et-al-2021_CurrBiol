#!/bin/bash
#
# Usage:
#   clean_fasta.sh <indir> <no_repeats.bed> <outdir> <njobs>
#
# Description:
#   Remove repetitive regions from the consensus sequences and calculate the
#   percentage of ambiguous sites in each sequence.
#
# Requirements:
#   bedtools
#   python3
#   parallel
#   samtools
#   'concatenate_fasta.py' (must be in same directory as 'clean_masked_fasta.sh')
#
# Important:
#   <no_repeats.bed> is the BED file containing non-repetitive regions in
#   scaffolds >= 1 Mbp (i.e. the file 'no_repeats_1mb_scaffolds.bed' generated
#   with 'process_assembly.sh').

# extract non-repetitive regions from masked FASTA
parallel -j $4 bedtools getfasta -fi {} -bed $2 -fo $3/{/.}.clean.fa ::: $1/*.fa

# concatenate FASTA regions per sample
parallel -j $4 python3 ./concatenate_fasta.py {} $3/{/.}.concat.fa ::: $3/*.clean.fa

# index FASTAs
parallel -j $4 samtools faidx {} ::: $3/*.concat.fa

# calculate proportion of N's per FASTA
for fasta in $3/*.concat.fa; do
  sample=$(basename ${fasta%.clean.concat.fa})
  seq_len=$(grep -v '^>' ${fasta} | tr -d '\n' | wc -c)
  n_count=$(grep -v '^>' ${fasta} | tr -cd N | wc -c)
  n_percent=$(python -c "print(f'{round((${n_count} / ${seq_len}) * 100, 2)}')")
  echo -e "${sample}\t${seq_len}\t${n_count}\t${n_percent}" >> $3/n_percent.tmp
done
cat <(echo -e "file\tlength_bp\t#N\t%N") <(sort -grk 4,4 $3/n_percent.tmp) > $3/n_percent.tbl

# remove intermediate files
rm $3/*.clean.fa $3/n_percent.tmp
