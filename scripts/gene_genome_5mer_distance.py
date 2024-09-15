#!/usr/bin/env python
from sys import argv
import argparse
import math
import os
import random
import time

import aux

def parse_args(argv):
    desc = 'Calculates the maximum euclidean distance between the 5mer distributions of genes and genomes'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='File containing 5mer distributions of bacterial genomes.')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def read_genome_5mer_dist(filename):
    genome_5mer_distr = {}
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:
            items = line.replace('[', '').replace(']', '').split('\t')
            genome_5mer_distr[items[0]] = {'distr1': [float(x) for x in items[1].split(', ')], 'distr2': [float(x) for x in items[2].split(', ')]}

    return genome_5mer_distr

def select_gene(grps):
    if os.path.isfile('/'.join(['example_data', 'separated_groups', grps[0], '5mer_distributions_gene.txt'])) and os.path.isfile('/'.join(['example_data', 'separated_groups', grps[1], '5mer_distributions_gene.txt'])):
        sample = random.sample(grps, 1)
        selection = sample[0]

    elif os.path.isfile('/'.join(['example_data','separated_groups', grps[0], '5mer_distributions_gene.txt'])) and not os.path.isfile('/'.join(['example_data', 'separated_groups', grps[1], '5mer_distributions_gene.txt'])):
        selection = grps[0]

    elif not os.path.isfile('/'.join(['example_data', 'separated_groups', grps[0], '5mer_distributions_gene.txt'])) and os.path.isfile('/'.join(['example_data', 'separated_groups', grps[1], '5mer_distributions_gene.txt'])):
        selection = grps[1]
    
    elif not os.path.isfile('/'.join(['example_data', 'separated_groups', grps[0], '5mer_distributions_gene.txt'])) and not os.path.isfile('/'.join(['example_data', 'separated_groups', grps[1], '5mer_distributions_gene.txt'])):
        selection = 'NA'

    return selection

def calc_gene_genome_distance(selection, genome_5mer_distr, key):
    with open('/'.join(['example_data', 'separated_groups', selection, '5mer_distributions_gene.txt']), 'r') as file:
        gene_list = [x.strip() for x in file.readlines()]

    gene_distr = [float(x) for x in gene_list[0].split("\t")]

    distance1 = math.dist(gene_distr, genome_5mer_distr[key]['distr1'])
    distance2 = math.dist(gene_distr, genome_5mer_distr[key]['distr2'])

    return max([distance1, distance2])


def main():
    arguments = parse_args(argv)
    genome_5mer_distr = read_genome_5mer_dist(arguments.input)
    gene_genome_dist = ['Gene_genome_5mer_distance']

    bar = aux.setup_progressbar(len(genome_5mer_distr.keys()))

    for i in range(len(genome_5mer_distr.keys())):
        key = list(genome_5mer_distr.keys())[i]
        grps = ['-'.join([key, 'grp1']), '-'.join([key, 'grp2'])]
        selection = select_gene(grps)

        if selection == 'NA':
            gene_genome_dist.append('NA')
            
            time.sleep(0.1)
            bar.update(i)
            continue

        else:
            max_distance = calc_gene_genome_distance(selection, genome_5mer_distr, key)
            gene_genome_dist.append(max_distance)

            time.sleep(0.1)
            bar.update(i)
    
    aux.export_output(gene_genome_dist, arguments.output)

if __name__ == '__main__':
    main()
