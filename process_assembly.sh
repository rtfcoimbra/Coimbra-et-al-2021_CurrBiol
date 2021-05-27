#!/bin/bash
#
# Usage:
#   process_assembly.sh <workdir>
#
# Description:
#   Index assembly, annotate repeats with RepeatMasker and RepeatModeler, mask
#   FASTA, extract scaffolds >= 1 Mbp from masked FASTA, and create a BED file
#   of non-repetitive regions in scaffolds >= 1 Mbp.
#
# Requirements:
#   bedtools
#   bwa
#   picard
#   repeatmasker (incl. buildSummary.pl, queryRepeatDatabase.pl)
#   repeatmodeler
#   samtools
#
# Note:
#   https://github.com/rmhubley/RepeatMasker/issues/5

# working directory
dir=$1

# decompress FASTA and change scaffold headers
zcat ${dir}/PLA01_assembly_pseudohap.fasta.gz \
  | sed -r 's/^>([0-9]+).*/>Scaffold_\1/' > ${dir}/kordofan-giraffe_assembly.fa

# genome assembly
asm=${dir}/kordofan-giraffe_assembly.fa

# index assembly
bwa index -a bwtsw ${asm}
samtools faidx ${asm}
picard CreateSequenceDictionary R=${asm} O=${asm/.fa/.dict}

# hard mask known repeats based on a repository query
RepeatMasker -s -engine ncbi -pa 32 -species 'cetartiodactyla' -dir ${dir} -gff ${asm} > ${dir}/repeatmasker_step1.log

# organize output files
mkdir ${dir}/repeatmasker_step1_results
mv ${dir}/*.{nhr,nin,nnd,nni,nog,nsq,translation} ${dir}/*.{cat.gz,gff,log,masked,out,tbl} ${dir}/repeatmasker_step1_results

# masked assembly
masked_asm=${dir}/repeatmasker_step1_results/kordofan-giraffe_assembly.fa.masked

# create a de novo repeat library for the masked assembly
BuildDatabase -engine ncbi -name ${dir}/kordofan-giraffe_assembly.masked.db ${masked_asm}
RepeatModeler -engine ncbi -pa 32 -database ${dir}/kordofan-giraffe_assembly.masked.db > ${dir}/repeatmodeler.log

# rename RepeatModeler output directory and organize output files
mv ${dir}/RM_* ${dir}/repeatmodeler_results
mkdir ${dir}/trfResults && mv ${dir}/trfResults* ${dir}/trfResults
mv ${dir}/{*.{db*,log},trfResults,unaligned.fa} ${dir}/repeatmodeler_results

# exclude repeats with unknown classification from de novo repeat library
awk '/#[^Unknown]/{flag=1} /#Unknown/{flag=0} flag' ${dir}/repeatmodeler_results/kordofan-giraffe_assembly.masked.db-families.fa > ${dir}/repeatmodeler_results/repmod_classified.lib

# mask predicted repeats based on a custom repeat library
RepeatMasker -s -engine ncbi -pa 32 -lib ${dir}/repeatmodeler_results/repmod_classified.lib -dir ${dir} -gff ${masked_asm} > ${dir}/repeatmasker_step2.log

# organize output files
mkdir ${dir}/repeatmasker_step2_results
mv ${dir}/*.{nhr,nin,nnd,nni,nog,nsq,translation} ${dir}/*.{cat.gz,gff,log,masked,out,tbl} ${dir}/repeatmasker_step2_results

# concatenate RepeatMasker output
tail -n +4 ${dir}/repeatmasker_step2_results/*.out \
  | cat ${dir}/repeatmasker_step1_results/*.out - > ${dir}/repeatmasker.combined.out

# generate detailed summary of RepeatMasker output
cut -f 1,2 ${asm}.fai > ${asm/.fa/.tsv}
buildSummary.pl -species 'cetartiodactyla' -genome ${asm/.fa/.tsv} -useAbsoluteGenomeSize ${dir}/repeatmasker.combined.out > ${dir}/repeatmasker.summary.tbl

# create a sorted BED file of identified repeat regions
tail -n +4 ${dir}/repeatmasker.combined.out \
  | sed 's/^\s*//' \
  | sed -E 's/\s+/\t/g' \
  | cut -f 5-7 \
  | awk 'BEGIN { OFS="\t" } { print $1, $2-1, $3 }' \
  | sort -V > ${dir}/repeats.sorted

# merge adjacent repeats
bedtools merge -i ${dir}/repeats.sorted > ${dir}/repeats_merged_adj.bed

# merge repeats within 1 Kbp from each other
bedtools merge -i ${dir}/repeats.sorted -d 1000 > ${dir}/repeats_merged_1kb.bed

# hard mask original assembly FASTA based on the previously generated BED file
bedtools maskfasta -fi ${asm} -bed ${dir}/repeats_merged_1kb.bed -fo ${asm/.fa/.masked.fa}

# index final hard masked assembly
samtools faidx ${asm/.fa/.masked.fa}

# create BED and BED4 files containing assembly scaffolds >= 1 Mbp
scaffs_1mb=${dir}/kordofan-giraffe_1mb_scaffolds.bed
awk '$2>=1000000 { print $0 }' ${asm/.fa/.tsv} \
  | awk -F'\t' '$1 = $1' 'OFS=\t0\t' > ${scaffs_1mb}
awk 'BEGIN { OFS="\t" } { $4 = $1 }1' ${scaffs_1mb} > ${scaffs_1mb/.bed/.bed4}

# create and index a new FASTA containing only scaffolds >= 1 Mbp from the final masked assembly
bedtools getfasta -name -fi ${asm/.fa/.masked.fa} -bed ${scaffs_1mb/.bed/.bed4} \
  | fold > ${dir}/kordofan-giraffe_1mb_scaffolds.masked.fa
samtools faidx ${dir}/kordofan-giraffe_1mb_scaffolds.masked.fa

# create a BED file of repeats present in scaffolds >= 1 Mbp
awk '{ print $1 }' ${scaffs_1mb} > ${dir}/scaffolds_1mb.list
while read scaff; do
  grep -P "${scaff}\t" ${dir}/repeats_merged_1kb.bed >> ${dir}/repeats_1mb_scaffolds.bed
done < ${dir}/scaffolds_1mb.list

# create a BED file of non-repetitive regions of scaffolds >= 1 Mbp
awk 'BEGIN { OFS="\t" } { print $1, $3 }' ${scaffs_1mb} > ${dir}/scaffolds_1mb.sizes
bedtools complement -i ${dir}/repeats_1mb_scaffolds.bed -g ${dir}/scaffolds_1mb.sizes > ${dir}/no_repeats_1mb_scaffolds.bed
