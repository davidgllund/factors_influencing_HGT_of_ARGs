#!/usr/bin/env python
from collections import ChainMap
from sys import argv
import argparse
import subprocess
import time

from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
import progressbar

import aux

def parse_args(argv):
    desc = 'Estimates the proportional mean difference in size between two sets of genomes'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Table prodiced by "indentify_horizontal_transfers.R".')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    parser.add_argument('--num_cores', required=False, default=2, type=int,
                        help='Number of cores to use for parallel processing.')
    
    arguments = parser.parse_args()

    return arguments

def calc_mean_size(genomes, results, paths):
    genomes_split = genomes.split(';')
    genome_size = []

    for j in range(len(genomes_split)):
        size = subprocess.check_output(' '.join(['cat', paths[genomes_split[j]], '| grep -v ">" | tr -d "\n" | wc -c']), shell=True, text=True)
        genome_size.append(int(size.strip()))

    results[genomes] = np.mean(genome_size)

    return results

def calc_genome_size_difference(input_data, mean_size, bar):
    results = pd.DataFrame(columns = ['Genome_size_difference'])

    for i in range(input_data.shape[0]):
        results.loc[i,'Genome_size_difference'] = abs(mean_size[input_data.loc[i, 'AssemblyAccessionID1']] - mean_size[input_data.loc[i, 'AssemblyAccessionID2']]) / (mean_size[input_data.loc[i, 'AssemblyAccessionID1']] + mean_size[input_data.loc[i, 'AssemblyAccessionID2']])

        time.sleep(0.1)
        bar.update(i)

    return results

def main():
    arguments = parse_args(argv)
    input_data = pd.read_csv(arguments.input, sep='\t')
    genome_paths = aux.filepath_to_dict('/home/dlund/HGT_2.0/analysis_updated/paths_tax_check_ok_no_contam.txt')
    accession_ids = aux.extract_accession_ids(input_data)

    mean_size = {}
    function_inputs = tqdm(set(accession_ids))
    mean_size = Parallel(n_jobs = arguments.num_cores)(delayed(calc_mean_size)(i, mean_size, genome_paths) for i in function_inputs)
    mean_size = dict(ChainMap(*mean_size))

    bar = aux.setup_progressbar(input_data.shape[0])
    results = calc_genome_size_difference(input_data, mean_size, bar)
    results.to_csv(arguments.output, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    main()
