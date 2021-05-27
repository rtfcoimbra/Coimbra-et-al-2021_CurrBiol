#!/bin/bash
#
# Usage:
#   generate_fasta.sh <ref.fa> <indir> <outdir> <ncpus>
#
# Description:
#   Calculate per-base depth from cleaned BAMs using sambamba and generate
#   consensus sequences in FASTA format with ANGSD.
#
# Requirements:
#   angsd
#   parallel
#   python3
#   sambamba
#   'site_depth_stats.py' (must be in same directory as 'generate_fasta.sh')
#
# Important:
#   <ref.fa> must be the masked Kordofan giraffe assembly with scaffolds < 1 Mbp
#   removed, i.e. the 'kordofan-giraffe_1mb_scaffolds.masked.fa' file generated
#   with 'process_assembly.sh'.

# iterate over BAM files
for bam in $2/*.clean.bam; do
  # get sample name
  sample=$(basename ${bam%.clean.bam})
  # calculate site depth per sample
  sambamba depth base -t $4 ${bam} | awk 'NR>1 { print $3 }' > $3/${sample}.depth
done

# calculate site depth statistics
parallel -j $4 "python3 site_depth_stats.py {} > {}.stats" ::: $3/*.depth

# iterate over BAM files
for bam in $2/*.clean.bam; do
  # get sample name
  sample=$(basename ${bam%.clean.bam})
  # get the 95th percentile of the site depth distribution
  max_dp=$(grep -Po '95th percentile: \K\d+' $3/${sample}.depth.stats)
  # set read filters
  filter_reads="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1"
  # set site filters
  filter_sites="-minMapQ 30 -minQ 30 -setMinDepthInd 4 -setMaxDepthInd ${max_dp}"
  # set angsd tasks
  todo="-doCounts 1 -doFasta 4 -iupacRatio 0.33 -basesPerLine 80"
  # command to generate FASTA file with IUPAC ambiguity codes
  echo "angsd -i ${bam} -ref $1 ${filter_reads} ${filter_sites} ${todo} -P 4 -out $3/${sample} &> $3/${sample}.log" >> $3/doFasta.jobs
done

# run ANGSD in parallel
cat $3/doFasta.jobs | parallel -j $(python -c "print(f'{round($4 / 4)}')")
