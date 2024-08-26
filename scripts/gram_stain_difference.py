#!/usr/bin/env python
import argparse
from sys import argv

import pandas as pd

def parse_args(argv):
    desc = 'Checks the gram staining of two bacterial taxa'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Table prodiced by "indentify_horizontal_transfers.R".')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def setup_dictionaries():
    gram_stain = {}
    with open('/home/dlund/HGT_inference_project/analysis/analysis_files/phylum_stain.txt') as f:
        for line in f:
            items = line.split('\t')
            gram_stain[items[0]] = items[1].rstrip()

    taxonomy = {}
    with open('/home/dlund/HGT_inference_project/analysis/analysis_files/taxonomy_table.txt') as f:
        for line in f:
            items = line.split('\t')
            taxonomy[items[2].strip()] = {'phylum': items[0], 'class': items[1]}

    return gram_stain, taxonomy

def measure_gram_stain(order1, order2, gram_stain, taxonomy):
    if taxonomy[order1]['phylum'] in gram_stain.keys() and taxonomy[order2]['phylum'] in gram_stain.keys():
        stain1 = gram_stain[taxonomy[order1]['phylum']]
        stain2 = gram_stain[taxonomy[order2]['phylum']]

        if stain1 == 'NA' or stain2 == 'NA':
            stain_diff = 'NA'

        else:
            if stain1 == 'Negative' and stain2 == 'Negative':
                stain_diff = 'NN'
    
            elif stain1 == 'Positive' and stain2 == 'Positive':
                stain_diff = 'PP'
        
            elif stain1 == 'Negative' and stain2 == 'Positive':
                stain_diff = 'NP'

            elif stain1 == 'Positive' and stain2 == 'Negative':
                stain_diff ='NP'

    else:
        stain_diff ='NA'

    return stain_diff

def main():
    arguments = parse_args(argv)
    input_data  = pd.read_csv(arguments.input, sep = '\t')
    gram_stain, taxonomy = setup_dictionaries()

    order1 = input_data['Order1']
    order2 = input_data['Order2']

    stain_diff = []
    for i in range(len(order1)):
        stain_diff.append(measure_gram_stain(order1[i], order2[i], gram_stain, taxonomy))

    gram_stain_diff = pd.DataFrame(data=stain_diff)
    gram_stain_diff.columns = ['Gram_stain_difference'] 
    gram_stain_diff.to_csv(arguments.output, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    main()