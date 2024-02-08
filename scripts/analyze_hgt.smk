#-------------------------------------------------------------------------------
# 0 IMPORT LIBRARIES
#-------------------------------------------------------------------------------
import glob
import subprocess

import pandas as pd

DIR, = glob_wildcards('{dir}/predicted-orfs-amino.fasta')

#-------------------------------------------------------------------------------
# 1 FUNCTIONS
#-------------------------------------------------------------------------------
def split_ids(table, column):    
    df = pd.DataFrame({'id': table.loc[1,column].split(';')})
    df[['id', 'species']] = df['id'].str.split('-', expand=True)
    df = df.drop('species', axis=1)

    return df

def extract_sequences(header_subset, fasta_subset, header_complete, fasta_complete):
    with open(fasta_complete) as original_fasta, open(fasta_subset, 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        i = 0
        for record in records:
            if header_complete[i] not in header_subset:
                i += 1
                continue
            else:
                SeqIO.write(record, corrected_fasta, 'fasta')
                i += 1
                continue

#-------------------------------------------------------------------------------
# 2 RULES
#-------------------------------------------------------------------------------
rule all:
    input:
        #expand('{dir}/tree_clustered.txt', dir=DIR)
        'observed_horizontal_transfers.txt',
        'null_distribution.txt'

rule preprocessing:
    input:
        nucleotides = '{dir}/predicted-orfs.fasta',
        proteins = '{dir}/predicted-orfs-amino.fasta'
    output:
        taxonomy = '{dir}/host_taxonomy.txt',
        fasta_w_species = '{dir}/args_w_species.fasta',
        clusters = directory('{dir}/clusters'),
        centroids = '{dir}/args_clustered.fasta',
        alignment = '{dir}/alignment_clustered.aln',
        tree = '{dir}/tree_clustered.txt',
        blastout = '{dir}/blastout.txt',
        headerlines = '{dir}/fasta_headers.txt'
    shell:
        '''
        python ../scripts/preprocessing.py --nucleotides {input.nucleotides} --proteins {input.proteins} --taxonomy {output.taxonomy} --fasta_w_species {output.fasta_w_species} --clusters {output.clusters} --centroids {output.centroids} --alignment {output.alignment} --tree {output.tree} --blastout {output.blastout}
        grep '>' {input.nucleotides} > {output.headerlines}
        '''

rule identify_horizontal_transfers:
    input:
        tree = '{dir}/tree_clustered.txt',
        taxonomy = '{dir}/host_taxonomy.txt',
        clusters = '{dir}/clusters'
    output:
        '{dir}/observed_hgt.txt'
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript ../scripts/identify_horizontal_transfers.R --input {input.tree} --taxonomy {input.taxonomy} --clusters {input.clusters} --output {output}
        '''

rule add_labels:
    input:
        table = '{dir}/observed_hgt.txt',
        taxonomy = '{dir}/host_taxonomy.txt',
        cl = '{dir}'
    output:
        tdist = temp('{dir}/taxonomic_distance.txt'),
        gclass = temp('{dir}/gene_class.txt'),
        split1 = temp('{dir}/s1.txt'),
        split2 = temp('{dir}/s2.txt'),
        table = temp('{dir}/table_w_labels.txt')
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript ../scripts/taxonomic_difference.R --taxonomy {input.taxonomy} --input {input.table} --output {output.tdist}
        python ../scripts/make_labels.py --list {input.table} --word {input.cl} --header "Gene class" --output {output.gclass}
        cat {input.table} | cut -f 1 > {output.split1}
        cat {input.table} | cut -f 2 > {output.split2}
        paste {output.split1} {output.gclass} {output.tdist} {output.split2} | tail -n +2 > {output.table}
        '''

rule extract_header:
    output:
        temp("header.txt")
    run:
        head = ""
        head += 'Node' + '\t' + 'Gene class' + '\t' + 'Taxonomic difference' + '\t' + 'Order1' + '\t' + 'Order2' + '\t' + 'Species1' + '\t' + 'Species2' + '\n'

        with open(output[0], 'w') as file:
            file.write(head)

rule combine_hgt_table:
    input:
        table = expand("{dir}/table_w_labels.txt", dir=DIR),
        header = "header.txt"
    output:
        temp("hgt_table1.txt")
    run:
        '''
        cat {input.header} */table_w_labels.txt > {output}
        '''

rule translate_accession_ids:
    input:
        table = 'hgt_table1.txt',
        index = '/home/dlund/index_files/ID_index.txt'
    output:
        ids = temp('assembly_ids.txt'),
        table = temp('hgt_table2.txt')
    shell:
        '''
        python ../scripts/translate_accession_ids.py --input {input.table} --index {input.index} --output {output}
        paste {input.table} {output.ids} > {output.table}
        '''

rule lookup_otus:
    input:
        table = 'hgt_table2.txt'
    output:
        emp = temp('otus_emp.txt'),
        gwmc = temp('otus_gwmc.txt'),
        table = temp('hgt_table3.txt')
    shell:
        '''
        python ../scripts/lookup_otus.py --input {input.table} --database emp --output {output.emp}
        python ../scripts/lookup_otus.py --input {input.table} --database gwmc --output {output.gwmc}
        paste {input.table} {output.emp} {output.gwmc} > {output.table}
        '''

rule genome_5mer_distance:
    input:
        'hgt_table3.txt'
    output:
        temp('genome_5mer_distance.txt')
    shell:
        '''
        python ../scripts/genome_5mer_distance.py --input {input} -k 5 --num_cores 1 --output {output}
        '''

rule separate_genome_groups:
    input:
        'hgt_table3.txt'
    output:
        dname = directory('separated_groups')
    run:
        table = pd.read_csv(input[0], sep="\t")
        
        subprocess.run('mkdir %s' %(output[0]), shell=True)
        
        for i in range(len(table.iloc[:,0])):
            subdir1 = '-'.join([table.iloc[i,0], table.iloc[i,1], 'grp1'])
            subprocess.run('mkdir %s/%s' %(output[0], subdir1), shell=True)
    
            ids1 = split_ids(table, "Species1")
            ids1.to_csv('%s/´%s/accession_ids.txt' %(output[0], subdir1), index=False, header=False)

            subdir2 = '-'.join([table.iloc[i,0], table.iloc[i,1], 'grp2'])
            subprocess.run('mkdir %s/%s' %(output[0], subdir2), shell=True)
    
            ids2 = split_ids(table, "Species2")
            ids1.to_csv('%s/´%s/accession_ids.txt' %(output[0], subdir2), index=False, header=False)

rule gene_genome_5mer_distance:
    input:
        dname = 'separated_groups',
        distr = 'genome_5mer_distance.txt'
    output:
        temp('gene_genome_5mer_distance.txt')
    run:
        subdir = glob.glob(input[0] + '/*')
        
        for d in subdir:
            gclass = d.split('-')[1]
            subprocess.run('while read line; do grep ">" %s/predicted-orfs.fasta >> %s/header_subset.txt' %(gclass, d), shell=True)

            with open('%s/header_subset.txt' %(d), 'r') as file:
                header_subset = file.readlines()

            with open('%s/fasta_headers.txt' %(gclass), 'r') as file:
                header_complete = file.readlines()

            extract_sequences(header_subset, '%s/nucleotides.fna' %(d), header_complete, '%s/predicted-orfs.fasta')

        subprocess.run('snakemake -s ../scripts/generate_gene_5mer_distributions.smk --cores 3 all', shell=True)
        subprocess.run('python ../scripts/gene_genome_5mer_distance.py --input %s --output %s' %(input[1], output[0]))

rule genome_size_difference:
    input:
        true = 'hgt_table3.txt',
        null = 'null_table1.txt'
    output:
        true = temp('genome_size_diff.txt'),
        null = temp('genome_size_diff_null.txt')
    shell:
        '''
        python ../scripts/genome_size_difference.py --input {input.true} --output {output.true} --num_cores 1
        python ../scripts/genome_size_difference.py --input {input.null} --output {output.null} --num_cores 1
        '''

rule cooccurrence:
    input:
        true = 'hgt_table3.txt',
        null = 'null_table1.txt'
    output:
        emp_true = temp('cooccurrence_emp.txt'),
        gwmc_true = temp('cooccurrence_gwmc.txt'),
        emp_null = temp('cooccurrence_emp_null.txt'),
        gwmc_null = temp('cooccurrence_gwmc_null.txt')
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript ../scripts/cooccurrence.R --database emp --input {input.true} --output {output.emp_true} --num_cores 1
        Rscript ../scripts/cooccurrence.R --database gwmc --input {input.true} --output {output.gwmc_true} --num_cores 1
        Rscript ../scripts/cooccurrence.R --database emp --input {input.null} --output {output.emp_null} --num_cores 1
        Rscript ../scripts/cooccurrence.R --database gwmc --input {input.null} --output {output.gwmc_null} --num_cores 1
        '''

rule gram_stain_difference:
    input:
        true = 'hgt_table3.txt',
        null = 'null_table1.txt'
    output:
        true = temp('gram_stain_diff.txt'),
        null = temp('gram_stain_diff_null.txt')
    shell:
        '''
        python ../scripts/genome_size_difference.py --input {input.true} --output {output.true} --num_cores 1
        python ../scripts/genome_size_difference.py --input {input.null} --output {output.null} --num_cores 1
        '''

rule complete_hgt_table:
    input:
        table = 'hgt_table3.txt',
        genome_dist = 'genome_5mer_distance.txt',
        gene_genome_dist = 'gene_genome_5mer_distance.txt',
        genome_size = 'genome_size_diff.txt',
        emp = 'cooccurrence_emp.txt',
        gwmc = 'cooccurrence_gwmc.txt',
        gram_stain = 'gram_stain_diff.txt'
    output:
        'observed_horizontal_transfers.txt'
    shell:
        '''
        paste {input.table} {input.genome_dist} {input.gene_genome_dist} {input.genome_size} {input.gram_stain} {input.emp} {input.gwmc} > {output}
        '''

rule create_sample_dict:
    input:
        '{dir}'
    output:
        '{dir}/sample_dict.pkl'
    shell:
        '''
        python ../scripts/create_sample_dict.py --gene_class {input} --output {output}
        '''

rule generate_null_distribution:
    input:
        expand("{dir}/sample_dict.pkl", dir=DIR)
    output:
        temp('null_table1.txt')
    shell:
        '''
        python ../scripts/generate_null_distribution.py --max_number 100000 --output {output}
        '''

rule complete_null_distribution:
    input:
        table = 'null_table1.txt',
        genome_size = 'genome_size_diff_null.txt',
        emp = 'cooccurrence_emp_null.txt',
        gwmc = 'cooccurrence_gwmc_null.txt',
        gram_stain = 'gram_stain_diff_null.txt'
    output:
        'null_distribution.txt'
    shell:
        '''
        paste {input.table} {input.genome_size} {input.gram_stain} {input.emp} {input.gwmc} > {output}
        '''
