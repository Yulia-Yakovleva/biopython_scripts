# USAGE: python3 extract_unique_sequences.py -i /PATH_TO_YOUR_FILES/example.fa

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Process arguments from user')
parser.add_argument('-i', '--infile', nargs='*',
                    help='Input file of files. WARNING: your input file(s) will be OVERWRITTEN!')
args = parser.parse_args()


def drop_nonunique_sequences(in_file):
    print(str(in_file) + ' is in process...')
    out_file = in_file + '_unique.faa'
    with open(in_file, 'a') as outFile:
        record_ids = list()
        print('Initial number of sequences in file: ', len(list(SeqIO.parse(in_file, 'fasta'))))
        for record in SeqIO.parse(in_file, 'fasta'):
            if record.id not in record_ids:
                record_ids.append(record.id)
                SeqIO.write(record, outFile, 'fasta')
        print('Final number of sequences after non-unique removal: ', len(list(SeqIO.parse(in_file, 'fasta'))))


for infile in args.infile:
    drop_nonunique_sequences(infile)
