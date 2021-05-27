from sys import argv

infile = argv[1]
outfile = argv[2]

with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    for line in fin.readlines():
        filename = line.strip().split('/')[-1]
        sample = filename.strip().split('.')[0]
        fout.write(f'{sample}\n')

    fin.close()
    fout.close()
