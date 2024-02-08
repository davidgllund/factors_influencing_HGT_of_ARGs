#!/usr/bin/env python
from sys import argv
import argparse
import time

import pandas as pd
import progressbar

import aux

def parse_args(argv):
    desc = 'Translate NCBI Nucleotide accession IDs to Assembly accession IDs'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--input', required=True,
                        help='Input file.')
    parser.add_argument('--index', required=True,
                        help='ID index file.')
    parser.add_argument('--output', required=True,
                        help='Output file.')
    
    arguments = parser.parse_args()

    return arguments

def translate_ids(species, accession_index):
    bar = aux.setup_progressbar(len(species))

    accession_ids = []
    for i in range(len(species)):
        all_ids = []
        
        if species[i] != species[i]:
            accession_ids.append("NA")
            continue

        for item in species[i].split(";"):
            if item == "NA":
                continue
            else:
                if item.split('_')[0] in accession_index.keys():
                    all_ids.append(accession_index[item.split('_')[0]])
                else:
                    all_ids.append("NA")
    
        unique_ids = list(set(all_ids))
        id_line = ""
        for j in range(len(unique_ids)):
                id_line += str(unique_ids[j])
            
                if j != len(unique_ids)-1:
                    id_line += ";"
            
        accession_ids.append(id_line.strip())

        time.sleep(0.1)
        bar.update(i)
    
    return accession_ids

def read_files(arguments):
    accession_index = {}
    with open(arguments.index) as f:
        for line in f:
            (val, key) = line.split("\t")
            accession_index[key.strip()] = val

    df_total = pd.read_csv(arguments.input, delimiter='\t')

    species1 = list(df_total["Species1"])
    species2 = list(df_total["Species2"])

    return accession_index, species1, species2

def main():
    arguments = parse_args(argv)

    accession_index, species1, species2 = read_files(arguments)

    accession_id1=translate_ids(species1, accession_index)
    accession_id2=translate_ids(species2, accession_index)

    id_table = "AssemblyAccessionID1" + "\t" + "AssemblyAccessionID2" + "\n"
    
    for i in range(len(accession_id1)):
        id_table += str(accession_id1[i]) + "\t" + str(accession_id2[i]) + "\n"

    output=open(arguments.output, 'w')
    output.write(id_table)
    output.close()

if __name__ == '__main__':
    main()