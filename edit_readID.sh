#!/bin/bash
#
# Usage:
#   edit_readID.sh
#
# Description:
#   Remove 10X Chromium barcode sequences annotated by 'process_10xReads.py'
#   from the read ID in the BAM file of sample ENP11.
#
# Requirements:
#   samtools

# remove barcode sequences from read ID
samtools view -h ENP11.sorted.bam \
  | sed -r 's/([A-Z]{16})://' \
  | samtools view -b -o ENP11.sorted.bam.edited -

# remove old BAM
rm ENP11.sorted.bam

# rename new BAM
mv ENP11.sorted.bam.edited ENP11.sorted.bam
