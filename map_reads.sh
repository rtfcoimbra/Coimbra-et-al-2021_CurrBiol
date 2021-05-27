#!/bin/bash
#
# Usage:
#   map_reads.sh <ref.fa> <indir> <outdir> <ncpus>
#
# Description:
#   Map reads against the Kordofan giraffe assembly using BWA and sort the
#   output BAMs with SAMtools.
#
# Requirements:
#   bwa
#   samtools

# reference assembly
ref=$1

# create an array of name patterns
patterns=(WA GNP ZNP SNR ETH MF RET ISC LVNP MA SGR MTNP BNP SUN KKR V23 ENP HNB WOAK)
# iterate over patterns in the array
for pattern in ${patterns[@]}; do
  # create sorted list of trimmed reads 1
  ls -1 $2/${pattern}*_1P.fq.gz >> $3/r1.list
  # create sorted list of trimmed reads 2
  ls -1 $2/${pattern}*_2P.fq.gz >> $3/r2.list
  # create a sorted template for a read group information list
  grep -Po "${pattern}.*" $3/r1.list >> $3/rg.list.tmp

  # prepare ID and SM tags in read group information list
  case ${pattern} in
    WA)
      id='peralta'
      ;;
    GNP|ZNP|SNR)
      id='antiquorum'
      ;;
    ETH|MF)
      id='camelopardalis'
      ;;
    RET|ISC)
      id='reticulata'
      ;;
    LVNP|MA|SGR)
      id='tippelskirchi'
      ;;
    MTNP|BNP|SUN|KKR|V23)
      id='giraffa'
      ;;
    ENP|HNB)
      id='angolensis'
      ;;
    *)
      id='okapia'
      ;;
  esac
  sed -ri "s/(${pattern}[Rot]*[0-9]*-*[0-9]*)([_L0-9]*)_1P.+/${id}\2\t\1/" $3/rg.list.tmp
done

# reformat read group information list and add PL, and LB tags
awk '{ print "@RG\\tID:"$1"\\tSM:"$2"\\tPL:ILLUMINA\\tLB:nebnext" }' $3/rg.list.tmp > $3/rg.list
sed -ie '/ENP11/ s/nebnext/chromium/' $3/rg.list
sed -ie '/\(MA1\|WOAK\)/ s/nebnext/truseq/' $3/rg.list
sed -ie '/ZNP01/ s/nebnext/bgi-1/' $3/rg.list
sed -ie '/\(WA746\|ZNP01\|ETH1\|MF06\|RET3\|RET6\|ISC08\|LVNP8-04\|KKR08\)/ s/nebnext/bgi-2/' $3/rg.list

# combine RG, R1, and R2 lists
paste -d" " $3/rg.list $3/r1.list $3/r2.list > $3/bwa.args

# read command arguments from 'bwa.args'
while IFS=" " read -r rg trimmed_r1 trimmed_r2; do
  # output directory and name
  sorted_bam=$3/$(basename ${trimmed_r1/_1P.fq.gz/.sorted.bam})
  # map reads against reference and sort BAMs
  bwa mem -M -t $4 -R ${rg} ${ref} ${trimmed_r1} ${trimmed_r2} 2>> $3/bwa-mem.log \
    | samtools sort -@ 4 -o ${sorted_bam} - 2>> $3/samtools-sort.log
done < $3/bwa.args
