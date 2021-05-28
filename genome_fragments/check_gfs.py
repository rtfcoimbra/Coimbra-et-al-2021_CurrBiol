import sys
from Bio import SeqIO

filename = sys.argv[1]  # input file name
frag_size = sys.argv[2]  # alignment length
n_seqs = sys.argv[3]  # number of sequences in the alignments
acceptable_seqs = 0
bad_seqs = 0

# calculate the %N per sequence, categorize them into 'bad' or 'acceptable',
# and count the number of 'bad' and 'acceptable' sequences
for seq_record in SeqIO.parse(sys.stdin, 'fasta'):
    sequence = str(seq_record.seq).upper()
    percent_n = round((sequence.count('N') / int(frag_size)) * 100, 2)
    print(f"{seq_record.id}\t{percent_n}%")
    # acceptable: %N is between 20% and 50%
    if 20 < percent_n <= 50:
        acceptable_seqs += 1
    # bad: %N is higher than 50%
    elif percent_n > 50:
        bad_seqs += 1

# write output files
with open('good.gfs', 'a') as fout1, open('bad.gfs', 'a') as fout2:
    # good alignments: do not contain 'bad' sequences and can contain up to 20%
    # of 'acceptable' sequences
    if bad_seqs == 0 and acceptable_seqs <= int(n_seqs) * 0.2:
        fout1.write(f"{filename}\t{acceptable_seqs}\t{bad_seqs}\n")
    # bad alignments: do not fit good GF criteria
    else:
        fout2.write(f"{filename}\t{acceptable_seqs}\t{bad_seqs}\n")

    fout1.close()
    fout2.close()
