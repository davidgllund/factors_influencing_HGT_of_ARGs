#!/usr/bin/env python
from sys import argv
import argparse
import time

import pandas as pd

import aux

def parse_args(argv):
    desc = 'Translates assembly accession IDs to representative OTU IDs'
    copyright = 'Copyright (c) David Lund 2024.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', '-i', required=True,
                        help='Table produced by "identify_horizontal_transfers.R.')
    parser.add_argument('--database', '-d', required=True,
                        help='Database to use ("emp" or "gwmc")')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def lookup(item, dictionary):
    lines = item.split(';')

    otus = []
    for line in lines:
        if line in dictionary.keys():
            otus.append(dictionary[line])

    if len(otus) == 0:
        otus.append("NA")

    merged_otus = ""
    for j in range(len(otus)):
        merged_otus += str(otus[j])
            
        if j is not len(otus)-1:
            merged_otus += ";"

    return merged_otus.strip()

def read_mapping_file(db):
    otu_mapping = {}

    if db == 'emp':
        with open('auxiliary_files/otu_mapping_file_emp.txt') as f:
            for line in f:
                (key, val) = line.split('\t')
                otu_mapping[key] = val.strip()

    elif db == 'gwmc':
        with open('auxiliary_files/otu_mapping_file_gwmc.txt') as f:
            for line in f:
                (key, val) = line.split('\t')
                otu_mapping[key] = val.strip()

    return otu_mapping

def generate_header(db):
    if db == 'emp':
        headerline = 'Matched_OTU_EMP1' + '\t' + 'Matched_OTU_EMP2' + '\n'

    elif db == 'gwmc':
        headerline = 'Matched_OTU_GWMC1' + '\t' + 'Matched_OTU_GWMC2' + '\n'

    return headerline

def main():
    arguments = parse_args(argv)
    df = pd.read_csv(arguments.input, delimiter='\t')

    otu_mapping = read_mapping_file(arguments.database)
    matched_otus = generate_header(arguments.database)

    bar = aux.setup_progressbar(df.shape[0])

    for i in range(df.shape[0]):
        accno1 = df.loc[i,'AssemblyAccessionID1']
        accno2 = df.loc[i,'AssemblyAccessionID2']

        if accno1 != accno1 and accno2 != accno2:
            matched_otus += 'NA' + '\t' + 'NA' + '\n'
    
        elif accno1 != accno1 and accno2 == accno2:
            matched_otus += 'NA' + '\t' + lookup(accno2, otu_mapping) + '\n'
    
        elif accno1 == accno1 and accno2 != accno2:
            matched_otus += lookup(accno1, otu_mapping) + '\t' + 'NA' + '\n'

        else:
            matched_otus += lookup(accno1, otu_mapping) + '\t' + lookup(accno2, otu_mapping) + '\n'

        time.sleep(0.1)
        bar.update(i)

    output=open(arguments.output, 'w')
    output.write(matched_otus)
    output.close()

if __name__ == '__main__':
    main()
