from sys import argv
import argparse
import glob
import subprocess

from Bio import SeqIO
import pandas as pd

def parse_args(argv):
    desc = ''
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Table prodiced by "indentify_horizontal_transfers.R".')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output directory.')
    
    arguments = parser.parse_args()

    return arguments

def split_ids(table, column, row):    
    df = pd.DataFrame({'id': table.loc[row, column].split(';')})
    df[['id', 'species']] = df['id'].str.split('-', expand=True).iloc[0,:2]
    df = df.drop('species', axis=1)

    return df

def extract_sequences(header_subset, fasta_subset, header_complete, fasta_complete):
    with open(header_subset, 'r') as file:
        h1 = file.readlines()

    with open(header_complete, 'r') as file:
        h2 = file.readlines()

    with open(fasta_complete) as original_fasta, open(fasta_subset, 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        i = 0
        for record in records:
            if h2[i] not in h1:
                i += 1
                continue
            else:
                SeqIO.write(record, corrected_fasta, 'fasta')
                i += 1
                continue

def main():    
    arguments = parse_args(argv)
    table = pd.read_csv(arguments.input, sep="\t")
    dir = arguments.output.split('/')[0]
        
    subprocess.run('mkdir %s' %(arguments.output), shell=True)
        
    for i in range(len(table.iloc[:,0])):
        ls = glob.glob('%s/separated_groups/*' %(dir))

        if len(list(filter(lambda x:table.iloc[i,0] in x, ls))) == 0:
            j = 1

        else: 
            j += 1

        subdir1 = '-'.join([table.iloc[i,0], table.iloc[i,1], str(j), 'grp1'])
        subprocess.run('mkdir %s/%s' %(arguments.output, subdir1), shell=True)
    
        ids1 = split_ids(table, "Species1", i)
        ids1.to_csv('%s/%s/accession_ids.txt' %(arguments.output, subdir1), index=False, header=False)

        subprocess.run('while read line; do grep "$line" data/%s.hmm/fasta_headers.txt >> %s/%s/header_subset.txt; done<%s/%s/accession_ids.txt' %(table.iloc[i,1], arguments.output, subdir1, arguments.output, subdir1), shell=True)

        extract_sequences('%s/%s/header_subset.txt' %(arguments.output, subdir1), '%s/%s/nucleotides.fna' %(arguments.output, subdir1), 'data/%s.hmm/fasta_headers.txt' %(table.iloc[i,1]), 'data/%s.hmm/predicted-orfs.fasta' %(table.iloc[i,1]))

        subdir2 = '-'.join([table.iloc[i,0], table.iloc[i,1], str(j), 'grp2'])
        subprocess.run('mkdir %s/%s' %(arguments.output, subdir2), shell=True)
    
        ids2 = split_ids(table, "Species2", i)
        ids2.to_csv('%s/%s/accession_ids.txt' %(arguments.output, subdir2), index=False, header=False)

        subprocess.run('while read line; do grep "$line" data/%s.hmm/fasta_headers.txt >> %s/%s/header_subset.txt; done<%s/%s/accession_ids.txt' %(table.iloc[i,1], arguments.output, subdir2, arguments.output, subdir2), shell=True)

        extract_sequences('%s/%s/header_subset.txt' %(arguments.output, subdir2), '%s/%s/nucleotides.fna' %(arguments.output, subdir2), 'data/%s.hmm/fasta_headers.txt' %(table.iloc[i,1]), 'data/%s.hmm/predicted-orfs.fasta' %(table.iloc[i,1]))

if __name__ == '__main__':
    main()