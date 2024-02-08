from sys import argv
import argparse
import os
import subprocess
import time

from Bio import SeqIO
import pandas as pd

import aux

def parse_args(argv):
    desc = 'Generates files required to run "identify_horizontal_transfers.R"'
    copyright = 'Copyright (c) David Lund 2023.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--nucleotides', required=True,
                        help='Nucleotide multifasta produced by fARGene [INPUT].')
    parser.add_argument('--proteins', required=True,
                        help='Protein multifasta produced by fARGene [INPUT].')
    parser.add_argument('--taxonomy', required=True,
                        help='File containing taxonomy of bacterial hosts [OUTPUT].')
    parser.add_argument('--fasta_w_species', required=True,
                        help='Protein multifasta with host species included in the headers [OUTPUT].')
    parser.add_argument('--clusters', required=True,
                        help='Cluster directory [OUTPUT].')
    parser.add_argument('--centroids', required=True,
                        help='Centroid multifasta file [OUTPUT].')
    parser.add_argument('--alignment', required=True,
                        help='Multiple sequence alignment file [OUTPUT].')
    parser.add_argument('--tree', required=True,
                        help='Phylogenetic tree file [OUTPUT].')
    parser.add_argument('--blastout', required=True,
                        help='Blast output file [OUTPUT].')
    
    arguments = parser.parse_args()

    return arguments

def restructure_headerlines(taxonomy, arguments):
    s = os.popen('grep ">" %s | tr -d ">"' %(arguments.nucleotides))
    headers = s.readlines()

    new_headers = []

    bar = aux.setup_progressbar(len(headers))

    for i in range(len(headers)):
        assembly_accno = '_'.join(headers[i].split('_')[0:2])
        contig_accno = '_'.join(headers[i].split('-')[-1].split('_')[1:]).rstrip()

        new_head = '-'.join([contig_accno, taxonomy[assembly_accno]['species'].replace(' ', '_')])
        new_headers.append(new_head)

        time.sleep(0.1)
        bar.update(i)

    return headers, new_headers

def compile_taxonomy(headers, taxonomy, arguments):
    print('Compiling taxonomy')

    df = pd.DataFrame(columns=['ID', 'superkingdom', 'phylum',  'class', 'order', 'family', 'genus', 'species'])
    bar = aux.setup_progressbar(len(headers))

    for i in range(len(headers)):
        assembly_accno = '_'.join(headers[i].split('_')[0:2])
        contig_accno = '_'.join(headers[i].split('-')[-1].split('_')[1:]).rstrip()

        new_head = '-'.join([contig_accno, taxonomy[assembly_accno]['species'].replace(' ', '_')])

        df.loc[len(df)] = {'ID': new_head,
                           'superkingdom': taxonomy[assembly_accno]['superkingdom'], 
                           'phylum': taxonomy[assembly_accno]['phylum'],  
                           'class': taxonomy[assembly_accno]['class'], 
                           'order': taxonomy[assembly_accno]['order'], 
                           'family': taxonomy[assembly_accno]['family'], 
                           'genus': taxonomy[assembly_accno]['genus'], 
                           'species': taxonomy[assembly_accno]['species']}
    
    time.sleep(0.1)
    bar.update(i)

    df.to_csv(arguments.taxonomy, sep = '\t', header = True, index = False)

def phylogenetic_analysis(new_headers, arguments):
    print('Performing phylogenetic analysis')

    i = 0

    with open('%s' %(arguments.proteins)) as original_fasta, open('%s' %(arguments.fasta_w_species), 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            record.id = new_headers[i]
            record.description = ""
            i += 1
            SeqIO.write(record, corrected_fasta, 'fasta')

    subprocess.run('usearch -cluster_fast %s -id 1 -centroids %s' %(arguments.fasta_w_species, arguments.centroids), shell=True)
    subprocess.run('mafft %s > %s' %(arguments.centroids, arguments.alignment), shell=True, check=False, capture_output=True)
    subprocess.run('FastTree %s > %s' %(arguments.alignment, arguments.tree), shell=True)

def generate_cluster_directory(arguments):
    print('Generating cluster directory')

    subprocess.run('mkdir %s' %(arguments.clusters), shell=True)
    subprocess.run('usearch -cluster_fast %s -id 1 -clusters %s/c_' %(arguments.fasta_w_species, arguments.clusters), shell=True)
    subprocess.run('for f in %s/c*; do name=$(grep ">" $f | head -1 | tr -d ">"); mkdir %s/$name; grep ">" $f | tr -d ">" > %s/$name/hidden.txt; mv $f %s/$name/$name.fna; done' %(arguments.clusters, arguments.clusters, arguments.clusters, arguments.clusters), shell=True)

def calc_sequence_similarity(arguments):
    print('Calculating sequence similarities')

    subprocess.run('blastp -query %s -db /home/dlund/index_files/card_db/card_proteins.fasta -out %s -max_target_seqs 1 -outfmt 6' %(arguments.centroids, arguments.blastout), shell=True)

def main():
    arguments = parse_args(argv)
    
    taxonomy_index = {}
    with open('/home/dlund/index_files/assembly_id_w_full_taxonomy.txt') as f:
        for line in f:
            adj = line.split('\t')
            taxonomy_index[adj[0]] = {'superkingdom': adj[1], 'phylum': adj[2],  'class': adj[3], 'order': adj[4], 'family': adj[5], 'genus': adj[6], 'species': adj[7].rstrip()}

    headers, headers_adj = restructure_headerlines(taxonomy_index, arguments)
    
    compile_taxonomy(headers, taxonomy_index, arguments)

    phylogenetic_analysis(headers_adj, arguments)
    
    generate_cluster_directory(arguments)

    calc_sequence_similarity(arguments)

if __name__ == '__main__':
    main()