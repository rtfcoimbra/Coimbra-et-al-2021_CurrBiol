import argparse
import logging
import random
import sys
import time
from Bio import Align, AlignIO, SeqIO
from os.path import basename


def read_fasta(fasta_files):
    """
    Read input files in FASTA format.

    Read multiple DNA consensus sequences in FASTA format and return a list of
    filenames, a dictionary of sequences, and a list of scaffolds.
    """
    filenames = []
    sequences = {}
    scaffolds = []
    for filepath in fasta_files:
        filename = basename(filepath)
        filenames.append(filename)
        sequences[filename] = SeqIO.index(filepath, 'fasta')
        scaffolds_per_fasta = sequences.get(filename).keys()
        scaffolds.extend(scaffolds_per_fasta)
        logging.info(f"{filename} successfully loaded!")

    scaffolds = set(scaffolds)

    return filenames, sequences, scaffolds


def read_bed(bed_file):
    """
    Read input BED file.

    Read a BED file to generate a query list from scaffold names. Start and end
    positions (columns 2 and 3) are not used!
    """
    queries = []
    with open(bed_file) as bed:
        for line in bed.readlines():
            queries.append(line.split('\t')[0])

    logging.info(f"BED file successfully loaded!")

    return queries


def build_alignment(filenames, sequences, scaffold):
    """
    Build scaffold alignment.

    Extract sequences of a given scaffold from each input FASTA and build a
    scaffold alignment containing all subjects in the input FASTA list.
    """
    seqs = []
    for filename in filenames:
        seq_to_add = sequences.get(filename).get(scaffold)
        seq_to_add.id = filename.split('.')[0]
        seqs.append(seq_to_add)

    logging.info(f"Building alignment for '{scaffold}'")

    min_len_alignment = min([len(seq) for seq in seqs])
    alignment = Align.MultipleSeqAlignment(
        [seq[:min_len_alignment] for seq in seqs])

    return alignment


def clean_alignment(alignment):
    """
    Remove ambiguities from alignment.

    Iterate over sites in the alignment and build a new alignment containing
    either only pure ATGC sites (-c) or sites with up to a specified proportion
    of N's (-c -n FLOAT).
    """
    site_length = len(alignment[:, 0])
    cleaned_alignment = Align.MultipleSeqAlignment(
        [seq[:0] for seq in alignment])

    if args.n_ratio:
        logging.info(f"Removing sites with > {int(args.n_ratio * 100)}% of "
                     + f"N's from '{alignment[0].name}'")
        for pos in range(0, len(alignment[0])):
            site = alignment[:, pos: pos + 1]
            site_nucleotides = alignment[:, pos]
            n_count = site_nucleotides.upper().count('N')
            n_ratio = n_count / site_length
            if n_ratio <= args.n_ratio:
                cleaned_alignment += site

    else:
        logging.info("Removing sites with ambiguities from "
                     + f"'{alignment[0].name}'")
        iupac = ['N', 'Y', 'R', 'K', 'M', 'W', 'S', 'B', 'D', 'H', 'V', '-']
        iupac_length = len(iupac)
        for pos in range(0, len(alignment[0])):
            site = alignment[:, pos: pos + 1]
            site_nucleotides = alignment[:, pos]
            bad_char = False
            if site_length > iupac_length:
                if any([char in site_nucleotides.upper() for char in iupac]):
                    bad_char = True
                    break
            else:
                for char in site:
                    if str(char.seq).upper() in 'NYRKMWSBDHV-':
                        bad_char = True
                        break
            if not bad_char:
                cleaned_alignment += site

    return cleaned_alignment


def chop_alignment(alignment, fragment_size):
    """
    Chop alignment into fragments of a given size.

    Chop alignment into fragments of a given size. If a fragment becomes
    shorter than the requested size it is discarded.
    """
    fragment_count = 0
    for pos in range(0, len(alignment[0]), fragment_size):
        fragment = alignment[:, pos: pos + fragment_size]
        if len(fragment[0]) < fragment_size:
            continue

        with open(f"GF{fragment_count + 1:05}_{fragment_size}bp_"
                  + f"{alignment[0].name}.fa", 'w') as fout:
            AlignIO.write(fragment, fout, 'fasta')
            fragment_count += 1

    logging.info(f"{fragment_count} genome fragments of {fragment_size} bp "
                 + f"were successfully generated from '{alignment[0].name}'!")


def random_fragments(alignments, fragment_size):
    """
    Randomly sample fragments from scaffold alignments.

    This function randomly samples an alignment from a list of scaffold
    alignments and then a genome fragment from it in each iteration. Genome
    fragments are saved as new FASTAs.
    """
    size_buffer = fragment_size * 2
    fragment_count = 0
    while fragment_count < args.random:
        random_alignment = random.choice(alignments)
        len_random_alignment = len(random_alignment[0])
        if len_random_alignment > size_buffer:
            random_start = random.randint(
                0, len_random_alignment - fragment_size)
            if (args.clean and
                    len_random_alignment - random_start <= size_buffer):
                continue

            elif (args.clean and
                    len_random_alignment - random_start > size_buffer):
                logging.info(f"Sampling a random {size_buffer} bp chunk "
                             + f"from '{random_alignment[0].name}' "
                             + "alignment")
                chunk = random_alignment[
                    :, random_start: random_start + size_buffer]
                processed_alignment = clean_alignment(chunk)
                if processed_alignment:
                    len_processed_alignment = len(processed_alignment[0])
                    if len_processed_alignment < fragment_size:
                        logging.info("Cleaned alignment length shorter "
                                     + "than fragment size... skipping")
                        continue

                else:
                    logging.warning(f"Failed to clean alignment of "
                                    + f"'{chunk[0].name}'!")
                    continue

                logging.info(f"Sampling a {fragment_size} bp fragment from "
                             + "the cleaned alignment chunk of "
                             + f"'{chunk[0].name}'")
                fragment = processed_alignment[:, 0:fragment_size]

            elif not args.clean:
                logging.info(f"Sampling a random {fragment_size} bp fragment "
                             + f"from the '{random_alignment[0].name}' "
                             + "alignment")
                fragment = random_alignment[
                    :, random_start: random_start + fragment_size]

        elif fragment_size < len_random_alignment <= size_buffer:
            if args.clean:
                processed_alignment = clean_alignment(random_alignment)
                if processed_alignment:
                    len_processed_alignment = len(processed_alignment[0])
                    if len_processed_alignment < fragment_size:
                        logging.info("Cleaned alignment length shorter "
                                     + "than fragment size... skipping")
                        continue

                    else:
                        logging.info(f"Samplig a random {fragment_size} bp "
                                     + "fragment from the cleaned alignment "
                                     + f" of '{random_alignment[0].name}'")

                else:
                    logging.warning(f"Failed to clean alignment of "
                                    + f"'{random_alignment[0].name}'!")
                    continue

            else:
                processed_alignment = random_alignment
                len_processed_alignment = len_random_alignment
                logging.info(f"Samplig a random {fragment_size} bp fragment "
                             + f"from the '{random_alignment[0].name}' "
                             + "alignment")

            random_start = random.randint(
                0, len_processed_alignment - fragment_size)
            fragment = processed_alignment[
                :, random_start: random_start + fragment_size]

        else:
            logging.info(f"Alignment of '{random_alignment[0].name}' shorter "
                         + "than fragment size... skipping")
            continue

        with open(f"GF{fragment_count:04}_{fragment_size}bp_"
                  + f"{fragment[0].name}_{random_start}.fa",
                  'w') as fout:
            AlignIO.write(fragment, fout, 'fasta')
            fragment_count += 1

    logging.info(f"{fragment_count} random genome fragments of "
                 + f"{fragment_size} bp were successfully generated!")


def main():
    """Set the workflow and call all other functions."""
    # read FASTA file(s) and extract filenames, sequences, scaffold names
    filenames, sequences, scaffolds = read_fasta(args.fasta)

    # read BED file and create scaffold query list
    queries = read_bed(args.bed) if args.bed else []

    # generate alignments list
    alignments = []
    for scaffold in scaffolds:
        # skip scaffolds not found in the input BED
        if args.bed and scaffold not in queries:
            logging.info(f"'{scaffold}' not in query list... skipping")
            continue

        else:
            # build alignment for each scaffold
            alignment = build_alignment(filenames, sequences, scaffold)
            alignments.append(alignment)

    if args.random:
        # randomly sample fragments of a given size from the alignments
        random_fragments(alignments, args.fragment_size)

    else:
        processed_alignments = []
        # remove ambiguities from alignments
        if args.clean:
            for alignment in alignments:
                cleaned_alignment = clean_alignment(alignment)
                if cleaned_alignment:
                    logging.info(f"Cleaned '{cleaned_alignment[0].name}' "
                                 + f"length: {len(cleaned_alignment[0])} bp")
                    processed_alignments.append(cleaned_alignment)

                else:
                    logging.warning(f"Failed to clean '{alignment[0].name}'!")
                    continue
        # take raw alignments
        else:
            for alignment in alignments:
                logging.info(f"'{alignment[0].name}' length: "
                             + f"{len(alignment[0])} bp")
                processed_alignments.append(alignment)

        # split alignments in fragments of a given size
        for processed_alignment in processed_alignments:
            chop_alignment(processed_alignment, args.fragment_size)


if __name__ == '__main__':
    print("\n### GENOME FRAGMENTS GENERATOR ###\n"
          + "                                  \n"
          + "(c) Fritjof Lammers, 2015-2019    \n"
          + "    Raphael Coimbra, 2019-2021    \n")

    # set up argument parser
    parser = argparse.ArgumentParser(
        description="""This program loads several genomic consensus sequences
            in FASTA format, aligns them by identical headers, and saves
            alignments for fragments of a given size. Output files are written
            to working directory.""")

    # add arguments to parser
    parser.add_argument(
        '-f',
        '--fasta',
        action='append',
        required=True,
        help="""Input file in FASTA format. Give multiple files with separate
            '-f' parameters.""",
        metavar='FILE')
    parser.add_argument(
        '-b',
        '--bed',
        required=False,
        help="BED file with scaffolds to be considered.",
        metavar='FILE')
    parser.add_argument(
        '-c',
        '--clean',
        action='store_true',
        required=False,
        help="""Remove sites containing any ambiguity characters from scaffold
            alignments.""")
    parser.add_argument(
        '-n',
        '--n_ratio',
        type=float,
        required=False,
        help="""Proportion of N's allowed per site. This will ignore any other
            ambiguity characters in the alignments. Requires '-c'.""",
        metavar='FLOAT')
    parser.add_argument(
        '-s',
        '--fragment_size',
        type=int,
        required=True,
        help="Length of the genome fragments to be generated.",
        metavar='INT')
    parser.add_argument(
        '-r',
        '--random',
        type=int,
        required=False,
        help="Number of fragments to be randomly sampled from the genome.",
        metavar='INT')

    args = parser.parse_args()

    # set error messages for dependent arguments
    if args.n_ratio and not args.clean:
        parser.error("argument -n/--n_ratio: requires -c/--clean.")
    # set error message for incorrect float value
    if args.n_ratio is not None:
        if not 0 < args.n_ratio < 1:
            parser.error("argument -n/--n_ratio: value should be a float "
                         + "between, but not including, 0 and 1.")

    # set up logging to file
    logging.basicConfig(
        filename='./genome_fragments_%sbp.log' % args.fragment_size,
        filemode='w',
        format='%(asctime)-16s %(levelname)-7s %(message)s',
        datefmt='%y-%m-%d  %H:%M',
        level=logging.DEBUG,)
    # define a handler which prints INFO messages or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set format for console use
    formatter = logging.Formatter(
        fmt='%(asctime)-16s %(levelname)-7s %(message)s',
        datefmt='%y-%m-%d  %H:%M')
    console.setFormatter(formatter)
    # add handler to root logger
    logging.getLogger('').addHandler(console)

    # log command
    logging.info("./genome_fragments.py "
                 + f"--fasta {args.fasta} "
                 + f"--bed {args.bed} "
                 + f"--clean {args.clean} "
                 + f"--n_ratio {args.n_ratio} "
                 + f"--fragment_size {args.fragment_size} "
                 + f"--random {args.random}")

    # record execution starting time
    start_time = time.time()

    # set seed for pseudo-random number generator
    if args.random:
        seed = random.randrange(sys.maxsize)
        random.Random(seed)
        logging.info(f"Random seed: {seed}")

    # call main function
    main()

    # log running time
    logging.info(f"Running time: {round(time.time() - start_time, 2)} seconds")
