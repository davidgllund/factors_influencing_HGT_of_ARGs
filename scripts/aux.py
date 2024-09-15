#!/usr/bin/env python
from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
import progressbar

def setup_progressbar(length):
    widgets = [' [',
            progressbar.Timer(format='elapsed time: %(elapsed)s'),
            '] ',
            progressbar.Bar('*'),' (',
            progressbar.ETA(), ') ',
            ]
    bar = progressbar.ProgressBar(max_value=length, widgets=widgets).start()

    return bar

def export_output(results, output_filename):
    strToExport = ''

    for line in results:
        strToExport += str(line) + '\n'

    with open(output_filename, 'w') as output:
        output.write(strToExport)

def get_kmers(seq, k):
    kmers = []
    
    for i in range(len(seq)-(k-1)):
        kmers.append(str(seq[i:i+k]))
        
    return kmers

def generate_possible_kmers(k):
    comb = list(product('ACGT', repeat=k))
    possible_kmers = []

    for i in range(len(comb)):
        possible_kmers.append(''.join(comb[i]))

    return possible_kmers

def extract_accession_ids(dataframe):
    ids = list(dataframe['AssemblyAccessionID1'])
    ids.extend(list(dataframe['AssemblyAccessionID2']))

    return ids

def filepath_to_dict(path):
    dictionary = {}
    with open(path) as f:
        for line in f:
            items = line.split('/')
            dictionary['_'.join(items[-1].split("_")[0:2])] = line.rstrip()

    return dictionary

def read_file(filename):
    with open(filename, 'r') as file:
        output = [x.strip() for x in file.readlines()]

    return output

def fasta_reader(filename):
    fasta_df = pd.read_csv(filename, sep='>', lineterminator='>', header=None)
    fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', n=1, expand=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
    fasta_df.index = list(['_'.join(x.split('_')[-2:]) for x in list(fasta_df['Accession'])])
    
    return fasta_df

def subset_fasta(fasta, headerfile):
    with open(headerfile, 'r') as f:
        header_subset = [x.replace('>', '').strip() for x in f.readlines()]
    
    header_subset2 = ['_'.join(x.split('_')[-2:]) for x in header_subset]

    fasta_subset = fasta.loc[header_subset2]

    return fasta_subset

def get_kmer_distribution(kmers, possible_kmers):
    kmers_counted = Counter(kmers)
    dist_list = []

    for kmer in possible_kmers:
        dist_list.append(kmers_counted[kmer]/len(kmers))

    distribution = pd.DataFrame({'fraction': dist_list}, index=possible_kmers)

    return distribution
