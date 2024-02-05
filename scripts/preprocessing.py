import os
import subprocess
import time

from Bio import SeqIO
import pandas as pd

import aux

def restructure_headerlines(taxonomy):
    s = os.popen('grep ">" predicted-orfs.fasta | tr -d ">"')
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

def compile_taxonomy(headers, taxonomy):
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

    df.to_csv('host_taxonomy.txt', sep = '\t', header = True, index = False)

def phylogenetic_analysis(new_headers):
    print('Performing phylogenetic analysis')

    i = 0

    with open('predicted-orfs-amino.fasta') as original_fasta, open('args_w_species.fasta', 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            record.id = new_headers[i]
            record.description = ""
            i += 1
            SeqIO.write(record, corrected_fasta, 'fasta')

    subprocess.run('usearch -cluster_fast args_w_species.fasta -id 1 -centroids args_clustered.fasta', shell=True)
    subprocess.run('mafft args_clustered.fasta > alignment_clustered.aln', shell=True, check=False, capture_output=True)
    subprocess.run('FastTree alignment_clustered.aln > tree_clustered.txt', shell=True)

def generate_cluster_directory():
    print('Generating cluster directory')

    subprocess.run('mkdir clusters/', shell=True)
    subprocess.run('usearch -cluster_fast args_w_species.fasta -id 1 -clusters clusters/c_', shell=True)
    subprocess.run('for f in clusters/c*; do name=$(grep ">" $f | head -1 | tr -d ">"); mkdir clusters/$name; grep ">" $f | tr -d ">" > clusters/$name/hidden.txt; mv $f clusters/$name/$name.fna; done', shell=True)

def calc_sequence_similarity():
    print('Calculating sequence similarities')

    subprocess.run('blastp -query args_clustered.fasta -db /home/dlund/index_files/card_db/card_proteins.fasta -out blastout.txt -max_target_seqs 1 -outfmt 6', shell=True)

def main():
    
    taxonomy_index = {}
    with open('/home/dlund/index_files/assembly_id_w_full_taxonomy.txt') as f:
        for line in f:
            adj = line.split('\t')
            taxonomy_index[adj[0]] = {'superkingdom': adj[1], 'phylum': adj[2],  'class': adj[3], 'order': adj[4], 'family': adj[5], 'genus': adj[6], 'species': adj[7].rstrip()}

    headers, headers_adj = restructure_headerlines(taxonomy_index)
    
    compile_taxonomy(headers, taxonomy_index)

    phylogenetic_analysis(headers_adj)
    
    generate_cluster_directory()

    calc_sequence_similarity()

if __name__ == '__main__':
    main()