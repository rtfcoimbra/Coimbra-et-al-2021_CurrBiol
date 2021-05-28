#!/bin/bash
#
# Usage:
#   iqtree_au_test.sh <infile>
#
# Description:
#   Perform AU test with IQ-TREE and extract p-AU values from the log files.
#
# Requirements:
#   iqtree
#   parallel
#   python3
#   scripts 'genome_fragments.py', 'iqtree_parallel.sh', 'iqtree_au_test.sh',
#   and 'parse_iqlog.py' must be in the same directory.
#
# Important:
#   If no file containing tree topologies is given, the AU test is skipped.

# path to file of tree topologies
alt_trees=$(find .. -name 'topologies.tree' -printf '%p')

echo "Tree reconstruction for $(basename $1) in progress..."
echo "Approximately unbiased (AU) tree topology test for $(basename $1) in progress..."

# perform AU tree topology test with 10,000 replicates using ultrafast model selection
iqtree -s $1 -o 'WOAK' -n 0 -z ${alt_trees} -zb 10000 -au &&

# extract p-AU values from IQ-TREE's log file
python3 parse_iqlog.py ${1}.log > ${1}.au
