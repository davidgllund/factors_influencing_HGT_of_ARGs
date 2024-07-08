#!/usr/bin/env python
from itertools import islice
from operator import itemgetter
from sys import argv
import argparse
import glob
import pickle
import random
import time

from Bio import SeqIO
import numpy as np
import pandas as pd

import aux

def parse_args(argv):
    desc = 'Creates a dictionary containing information about predicted ARG clusters'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--taxonomy', required=True,
                        help='Table containing host taxonomy.')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def setup_dictionaries(taxpath, genomepaths, gene_class):
    taxonomy = {}
    
    with open(taxpath) as f:
        for line in islice(f, 1, None):
            items = line.split('\t')
            taxonomy[items[0]] = {'phylum': items[2], 'class': items[3], 'order': items[4], 'family': items[5], 'genus': items[6], 'species': items[7].rstrip()}

    genome_file_paths = {}
    gene_file_paths = {}

    with open(genomepaths) as f:
        for line in f:
            items = line.split('/')
            genome_file_paths['_'.join(items[8].split('_')[0:2])] = line.rstrip()
            gene_file_paths['_'.join(items[8].split('_')[0:2])] = '/' + '/'.join(line.split('/')[:-1]) + '/' + gene_class + '.hmm/predictedGenes/predicted-orfs.fasta'

    return taxonomy, genome_file_paths, gene_file_paths

def subset_data(cluster_data):
    if len(set(list(cluster_data['order']))) > 1:
        selected_order = random.sample(list(cluster_data['order']), 1)[0]
        cluster_data = cluster_data.loc[cluster_data['order'] == selected_order]

    if len(list(cluster_data['assembly_accession'])) > 1000:
        selected_ids = random.sample(list(cluster_data['assembly_accession']), 1000)
        cluster_data = cluster_data.loc[selected_ids,:]

    return cluster_data

def calc_gene_kmer_distributions(cluster_data, k, file_path_dict, possible_kmers):
    predicted_args = aux.fasta_reader(itemgetter(*tuple(random.sample(list(cluster_data['assembly_accession']),1)))(file_path_dict))
    predicted_args = predicted_args.drop(labels=list(set(predicted_args.index).difference(list(cluster_data['acc']))))

    gene_kmers = aux.get_kmers(predicted_args.iloc[0, 1], k)
    distribution = aux.get_kmer_distribution(gene_kmers, possible_kmers)

    results = pd.DataFrame(data=list(distribution['fraction']))

    return results

def calc_genome_kmer_distributions(cluster_data, k, file_path_dict, possible_kmers):
    kmer_distributions_genome = []

    for j in range(len(list(cluster_data['assembly_accession']))):
        host_genome = SeqIO.to_dict(SeqIO.parse(file_path_dict[list(cluster_data['assembly_accession'])[j]], 'fasta'))

        genome_kmers = []
        for k in range(len(host_genome)):
            seq_entry = host_genome[list(host_genome.keys())[k]]
            genome_kmers.extend(aux.get_kmers(seq_entry.seq, 5))

        distribution = aux.get_kmer_distribution(genome_kmers, possible_kmers)
        kmer_distributions_genome.append(list(distribution['fraction']))

    results = pd.DataFrame(data=kmer_distributions_genome)

    return results

def main():
    arguments = parse_args(argv)
    genomepaths = '/storage/dlund/HGT_inference_project/paths_taxonomy_check_ok_no_contam.txt'
    gene_class = arguments.taxonomy.split('/')[1]
    taxonomy, genome_file_paths, gene_file_paths = setup_dictionaries(arguments.taxonomy, genomepaths, gene_class)
    clusters = glob.glob('example_data/%s/clusters/*' %(gene_class))
    all_kmers = aux.generate_possible_kmers(5)

    sample_dict = {}
    bar = aux.setup_progressbar(len(clusters))

    for i in range(len(clusters)):
        key = clusters[i].split('/')[-1]

        if key not in taxonomy.keys():
            time.sleep(0.1)
            bar.update(i)
            continue

        else:
            species = aux.read_file('%s/hidden.txt' %(clusters[i]))
            assembly_accession = aux.read_file('%s/assembly_accessions.txt' %(clusters[i]))

            if len(assembly_accession) == 0:
                time.sleep(0.1)
                bar.update(i)
                
                continue

            cluster_df = pd.DataFrame({'assembly_accession': assembly_accession, 'acc': [x.split('-')[0] for x in species], 'order': [taxonomy.get(x).get('order') for x in species], 'species': species}, index=assembly_accession)
            cluster_df = subset_data(cluster_df)

            gene_kmer_distr = calc_gene_kmer_distributions(cluster_df, 5, gene_file_paths, all_kmers)
            genome_kmer_distr = calc_genome_kmer_distributions(cluster_df, 5, genome_file_paths, all_kmers)

            sample_dict[key] = {'gene_class': gene_class, 'order': list(cluster_df['order']), 'assembly_accessions': list(cluster_df['assembly_accession']), 'species': list(cluster_df['species']), '5mer_distribution_genome': list(np.mean(genome_kmer_distr, axis=0)), '5mer_distribution_gene': gene_kmer_distr}

            time.sleep(0.1)
            bar.update(i)

    with open(arguments.output, 'wb') as fp:
        pickle.dump(sample_dict, fp)

if __name__ == '__main__':
    main()
