#!/bin/bash
#
# Usage:
#   remove_10Xbarcodes.sh <proc10xG_dir> <indir> <outdir>
#
# Description:
#   Remove 10X Chromium barcodes from sample ENP11 with proc10xG.
#
# Requirements:
#   proc10xG
#   python2

# extract GEM barcodes and trim primer from read 1; compare barcode sequence to
# a white list; append status, barcode and trimmed sequence to read ID; only
# output reads with status 'MATCH' or 'MISMATCH1'
$1/process_10xReads.py \
  -1 $2/FCHGVVFBBXX_L6_CWHPEI17040130_1.fq.gz \
  -2 $2/FCHGVVFBBXX_L6_CWHPEI17040130_2.fq.gz \
  --bctrim 16 \
  --trim 7 \
  --output $3/ENP11 \
  2> $3/process_10xReads.log
