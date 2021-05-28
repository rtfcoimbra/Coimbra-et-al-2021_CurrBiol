from sys import argv

infile = argv[1]
outfile = argv[2]
last_chrom = ''

with open(infile) as fin, open(outfile, 'w') as fout:
    for line in fin.readlines():
        if line.startswith('>'):  # FASTA header line
            chrom = line.strip().split(':')[0]  # extract chrom. id from header

            print(f"Concatenating chromosome '{chrom.replace('>', '')}'...")

            if not last_chrom:  # first chromosome
                fout.write(f'{chrom}\n')  # write header
            elif chrom != last_chrom:  # second chromosome onwards
                fout.write(f'\n{chrom}\n')  # write header
            else:
                pass

        else:  # sequence line
            fout.write(line.strip())  # write sequence

        print("Done!")

        last_chrom = chrom  # remember last concatenated chromosome
