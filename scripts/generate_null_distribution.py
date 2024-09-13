#!/usr/bin/env python
from sys import argv
import argparse
import glob
import math
import os
import pickle
import random
import time

import pandas as pd

import aux

def parse_args(argv):
    desc = 'Generates pairs of gene clusters carried by hosts from distantly related taxa'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--max_number', required=False, default=100000,
                        help='Maximum number of pairs that will be generated (if not set, the script will continue until no more pairs can be generated).')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def read_otu_mapping(path):
    otu_mapping = {}
    with open(path) as f:
        for line in f:
            (key, val) = line.split("\t")
            otu_mapping[key] = val.strip()

    return otu_mapping

def read_taxonomy(path):
    taxonomy = {}
    with open(path) as f:
        for line in f:
            items = line.split('\t')
            taxonomy[items[2].strip()] = {'phylum': items[0], 'class': items[1]}

    return taxonomy

def lookup_otus(item, dictionary):
    lines = item.split(';')

    otus = []
    for line in lines:
        if line in dictionary.keys():
            otus.append(dictionary[line])

    if len(otus) == 0:
        otus.append('NA')

    merged_otus = ''
    for j in range(len(otus)):
        merged_otus += str(otus[j])
            
        if j != len(otus)-1:
            merged_otus += ';'

    return merged_otus.strip()

def measure_taxonomic_distance(order1, order2, taxonomy):
    if order1 not in taxonomy.keys() or order2 not in taxonomy.keys():
        diff = 'NA'
    elif taxonomy[order1]['phylum'] != taxonomy[order2]['phylum']:
        diff = 'phylum'
    elif taxonomy[order1]['phylum'] == taxonomy[order2]['phylum'] and taxonomy[order1]['class'] != taxonomy[order2]['class']:
        diff = 'class'
    else:
        diff = 'order'
    
    return diff

def find_antibiotic_class(gene_class):
    if gene_class in ['aac2p', 'aac3_class1', 'aac3_class2', 'aac6p_class1', 'aac6p_class2', 'aac6p_class3', 'aph2b', 'aph3p', 'aph6']:
        antibiotic_class = 'aminoglycoside'
    
    elif gene_class in ['class_A', 'class_B_1_2', 'class_B_3', 'class_C', 'class_D_1', 'class_D_2']:
        antibiotic_class = 'beta_lactam'
    
    elif gene_class in ['erm_typeA', 'erm_typeF', 'mph']:
        antibiotic_class = 'macrolide'
    
    elif gene_class in ['qnr']:
        antibiotic_class = 'quinolone'
    
    elif gene_class in ['tet_efflux', 'tet_enzyme', 'tet_rpg']:
        antibiotic_class = 'tetracycline'

    else:
        antibiotic_class = 'NA'
    
    return antibiotic_class

def combine_sample_dicts():
    combined_dict = {}
    categories = [dir for dir in glob.glob("example_data/*") if os.path.isdir( dir)]

    for item in categories:
        gene_class = item.split('.')[0]
        path = '%s/sample_dict.pkl' % (item)

        if os.path.isfile(path):
            with open(path, 'rb') as handle:
                sample_dict = pickle.load(handle)

            for key in sample_dict.keys():
                if sample_dict[key]['order'][0] != 'NA':
                    combined_dict['%s-%s' % (key, gene_class)] = sample_dict[key]
    
    return combined_dict

def setup_dataframe(dictionary):
    leaves = list(dictionary.keys())
    orders = []

    for leaf in leaves:
        orders.append(dictionary[leaf]['order'][0])

    df = pd.DataFrame(data = orders)
    df.index = leaves

    return df

def extract_information(key, dictionary):
    data = dictionary[key[0]]
    order = data['order'][0]
    species = ';'.join(data['species'])
    assembly_accession = ';'.join(data['assembly_accessions'])

    return data, order, species, assembly_accession

def calc_gene_genome_5mer_distance(dictionary, selected, data1, data2):
    if len(dictionary[selected[0]]['5mer_distribution_gene']) == 1024:
        gene = list(dictionary[selected[0]]['5mer_distribution_gene'].iloc[:,0])

    else:
        gene = list(dictionary[selected[0]]['5mer_distribution_gene'].iloc[0,:])

    gene_genome_5mer_dist = max([math.dist(data1['5mer_distribution_genome'], gene), math.dist(data2['5mer_distribution_genome'], gene)])

    return gene_genome_5mer_dist

def main():
    arguments = parse_args(argv)
    otu_mapping_emp = read_otu_mapping('auxiliary_files/otu_mapping_file_emp.txt')
    otu_mapping_gwmc = read_otu_mapping('auxiliary_files/otu_mapping_file_gwmc.txt')
    taxonomy = read_taxonomy('auxiliary_files/order_taxonomy.txt')
    combined_dict = combine_sample_dicts()
    df = setup_dataframe(combined_dict)
    bar = aux.setup_progressbar(int(arguments.max_number))
    
    null_data = []

    for i in range(int(arguments.max_number)):
        if df.shape[0] == 0:
            break
        else:
            key1 = random.sample(list(df.index), 1)
            data1, order1, species1, assembly_accession1 = extract_information(key1, combined_dict)

            if len(list(df.index[df.loc[:,0] != order1])) == 0:
                break

            else:
                key2 = random.sample(list(df.index[df.loc[:,0] != order1]), 1)
                data2, order2, species2, assembly_accession2 = extract_information(key2, combined_dict)

                tax_diff = measure_taxonomic_distance(order1, order2, taxonomy)
                otus_emp1 = lookup_otus(assembly_accession1, otu_mapping_emp)
                otus_emp2 = lookup_otus(assembly_accession2, otu_mapping_emp)
                otus_gwmc1 = lookup_otus(assembly_accession1, otu_mapping_gwmc)
                otus_gwmc2 = lookup_otus(assembly_accession2, otu_mapping_gwmc)
                genome_5mer_dist = math.dist(data1['5mer_distribution_genome'], data2['5mer_distribution_genome'])

                selected = random.sample([key1, key2], 1)[0]
                gene_genome_5mer_dist = calc_gene_genome_5mer_distance(combined_dict, selected, data1, data2)
                gene_class = combined_dict[selected[0]]['gene_class']
                antibiotic_class = find_antibiotic_class(gene_class)

                null_data.append(['Null_%s' %(i), gene_class, antibiotic_class, tax_diff, order1, order2, species1, species2, assembly_accession1, assembly_accession2, otus_emp1, otus_emp2, otus_gwmc1, otus_gwmc2, genome_5mer_dist, gene_genome_5mer_dist])
                df = df.drop(key1[0])
                df = df.drop(key2[0])

                time.sleep(0.1)
                bar.update(i)

    null_table = pd.DataFrame(data=null_data)
    null_table.columns = ['Node', 
                          'Gene class', 
                          'Antibiotic class', 
                          'Taxonomic difference', 
                          'Order1', 
                          'Order2', 
                          'Species1', 
                          'Species2', 
                          'AssemblyAccessionID1', 
                          'AssemblyAccessionID2', 
                          'Matched_OTU_EMP1', 
                          'Matched_OTU_EMP2', 
                          'Matched_OTU_GWMC1', 
                          'Matched_OTU_GWMC2', 
                          'Genome_5mer_distance', 
                          'Gene_genome_5mer_distance']
    
    null_table.to_csv(arguments.output, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    main()
