#!/usr/bin/env python
from collections import ChainMap
from sys import argv
import argparse
import math
import time
import gzip

from Bio import SeqIO
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd

import aux

def parse_args(argv):
    desc = 'Calculates the mean 5mer distributions of genomes involved in horizontal gene transfer and the distance beween them'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Table prodiced by "indentify_horizontal_transfers.R".')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    parser.add_argument('-k', required=True, type=int,
                        help='k to use for splitting sequence into kmers.')
    parser.add_argument('--num_cores', required=False, default=2, type=int,
                        help='Number of cores to use for parallel processing')
    
    arguments = parser.parse_args()
    
    return arguments

def calc_kmer_distribution(genomes, k, possible_kmers, results, paths):
    genomes_split = genomes.split(';')
    kmer_distribution = pd.DataFrame(columns = possible_kmers, index = genomes_split)

    for j in range(len(genomes_split)):
        if genomes_split[j] in paths.keys():
            host_genome = SeqIO.to_dict(SeqIO.parse(gzip.open(paths[genomes_split[j]], 'rt'), 'fasta'))

            genome_kmers = []
            for l in range(len(host_genome)):
                seq_entry = host_genome[list(host_genome.keys())[l]]
                genome_kmers.extend(aux.get_kmers(seq_entry.seq, k))

            distribution = aux.get_kmer_distribution(genome_kmers, possible_kmers)
            kmer_distribution.loc[genomes_split[j]] = list(distribution['fraction'])

    results[genomes] = list(np.mean(kmer_distribution, axis=0))

    return results

def calc_kmer_distance(input_data, kmer_distributions, bar):
    kmer_distribution1 = []
    kmer_distribution2 = []
    genome_kmer_distance = []

    for i in range(input_data.shape[0]):
        kmer_distribution1.append(kmer_distributions[input_data.loc[i, 'AssemblyAccessionID1']])
        kmer_distribution2.append(kmer_distributions[input_data.loc[i, 'AssemblyAccessionID2']])
        genome_kmer_distance.append(math.dist(kmer_distributions[input_data.loc[i, 'AssemblyAccessionID1']], kmer_distributions[input_data.loc[i, 'AssemblyAccessionID2']]))

        time.sleep(0.1)
        bar.update(i)
    
    results = pd.DataFrame({'Kmer_distribution1': kmer_distribution1, 'Kmer_distribution2': kmer_distribution2, 'Genome_kmer_distance': genome_kmer_distance})
    event_ids = []
    for i in range(input_data.shape[0]):
        if '-'.join([input_data.loc[i, 'Node'], input_data.loc[i, 'Gene class'], '1']) not in event_ids:
            j = 1
        else:
            j += 1

        event_ids.append('-'.join([input_data.loc[i, 'Node'], input_data.loc[i, 'Gene class'], str(j)]))

    results.index = event_ids

    return results

def main():
    arguments = parse_args(argv)
    input_data = pd.read_csv(arguments.input, delimiter='\t')
    genome_paths = aux.filepath_to_dict('auxiliary_files/genome_filepaths.txt')
    accession_ids = aux.extract_accession_ids(input_data)
    all_kmers = aux.generate_possible_kmers(arguments.k)

    print('Generating mean kmer distributions')
    function_inputs = tqdm(set(accession_ids))
    mean_kmer_distributions = {}
    mean_kmer_distributions = Parallel(n_jobs = arguments.num_cores)(delayed(calc_kmer_distribution)(i, arguments.k, all_kmers, mean_kmer_distributions, genome_paths) for i in function_inputs)
    mean_kmer_distributions = dict(ChainMap(*mean_kmer_distributions))

    print('Calculating kmer distance')
    bar = aux.setup_progressbar(input_data.shape[0])
    results = calc_kmer_distance(input_data, mean_kmer_distributions, bar)
    results.to_csv(arguments.output, sep = '\t', header = True, index = True)

if __name__ == '__main__':
    main()