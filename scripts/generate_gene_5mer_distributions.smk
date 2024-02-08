# This pipeline is used to calculate the 5mer distribution of resistance
# genes predicted by fARGene.  
#-------------------------------------------------------------------------------
# 0 IMPORT LIBRARIES
#-------------------------------------------------------------------------------
from collections import Counter
import glob

from Bio import SeqIO
import numpy as np
import pandas as pd

DIR, = glob_wildcards('separated_groups/{dir}/nucleotides.fna')

#-------------------------------------------------------------------------------
# 1 FUNCTIONS
#-------------------------------------------------------------------------------
def get_kmers(seq, k):
    kmers = []
    
    for i in range(len(seq)-(k-1)):
        kmers.append(str(seq[i:i+k]))
        
    return(kmers)

def get_possible_kmers(k):
    comb = list(product('ACGT', repeat=k))
    possible_kmers = []
    for i in range(len(comb)):
        possible_kmers.append(''.join(comb[i]))

    return possible_kmers

#-------------------------------------------------------------------------------
# 2 RULES
#-------------------------------------------------------------------------------
rule all:
    input:
        expand('separated_groups/{dir}/5mer_distributions_gene.txt', dir=DIR)

rule cluster_genes:
    input:
        'separated_groups/{dir}/nucleotides.fna'
    output:
        temp('separated_groups/{dir}/nucleotides_unique.fna')
    shell:
        '''
        usearch -cluster_fast {input} -id 1 -centroids {output}
        '''

rule generate_kmer_distributions:
    input:
        'separated_groups/{dir}/nucleotides_unique.fna'
    output:
        'separated_groups/{dir}/5mer_distributions_gene.txt'
    run:
        genes = SeqIO.to_dict(SeqIO.parse(input[0], 'fasta'))
        k = 5
        all_kmers = get_possible_kmers(k)

        kmer_distributions_gene = []
        gene_labels = []

        for i in range(len(genes)):
            seq_entry = genes[list(genes.keys())[i]]
            kmers_in_gene = get_kmers(seq_entry.seq, k)
            num_kmers = Counter(kmers_in_gene)
    
            d = np.zeros(shape=len(all_kmers))
            distr = pd.DataFrame(data=d)
            distr.index = all_kmers
            distr.columns = ['fraction']

            for kmer in all_kmers:
                distr.loc[kmer] = num_kmers[kmer]/len(kmers_in_gene)
        
            kmer_distributions_gene.append(list(distr['fraction']))
            gene_labels.append('_'.join(seq_entry.id.split('_')[0:3]))
            
        results_gene = pd.DataFrame(data=kmer_distributions_gene)
        results_gene.index = gene_labels
        results_gene.to_csv(output[0], sep = '\t', header = False, index = False)
        