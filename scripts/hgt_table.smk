#-------------------------------------------------------------------------------
# 0 IMPORT LIBRARIES
#-------------------------------------------------------------------------------
import glob
import subprocess

import pandas as pd

DIR, = glob_wildcards('example_data/{dir}/predicted-orfs-amino.fasta')

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
        'observed_horizontal_transfers.txt',
        'null_distribution.txt'

rule preprocessing:
    input:
        nucleotides = 'example_data/{dir}/predicted-orfs.fasta',
        proteins = 'example_data/{dir}/predicted-orfs-amino.fasta'
    output:
        taxonomy = 'example_data/{dir}/host_taxonomy.txt',
        fasta_w_species = 'example_data/{dir}/args_w_species.fasta',
        clusters = directory('example_data/{dir}/clusters'),
        centroids = 'example_data/{dir}/args_clustered.fasta',
        alignment = 'example_data/{dir}/alignment_clustered.aln',
        tree = 'example_data/{dir}/tree_clustered.txt',
        blastout = 'example_data/{dir}/blastout.txt',
        headerlines = 'example_data/{dir}/fasta_headers.txt'
    shell:
        '''
        python scripts/preprocessing.py --nucleotides {input.nucleotides} --proteins {input.proteins} --taxonomy {output.taxonomy} --fasta_w_species {output.fasta_w_species} --clusters {output.clusters} --centroids {output.centroids} --alignment {output.alignment} --tree {output.tree} --blastout {output.blastout}
        grep '>' {input.nucleotides} > {output.headerlines}
        '''

rule identify_horizontal_transfers:
    input:
        tree = 'example_data/{dir}/tree_clustered.txt',
        taxonomy = 'example_data/{dir}/host_taxonomy.txt',
        clusters = 'example_data/{dir}/clusters',
        blast = 'example_data/{dir}/blastout.txt'
    output:
        table = 'example_data/{dir}/observed_hgt.txt',
        pdf = 'example_data/{dir}/tree_annotated.pdf'
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript scripts/identify_horizontal_transfers.R --input {input.tree} --taxonomy {input.taxonomy} --clusters {input.clusters} --blast {input.blast} --pdf {output.pdf} --output {output.table}
        '''

rule calc_taxonomic_distance:
    input:
        table = 'example_data/{dir}/observed_hgt.txt',
        taxonomy = 'example_data/{dir}/host_taxonomy.txt'
    output:
        tdist = temp('example_data/{dir}/taxonomic_distance.txt')
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript scripts/taxonomic_difference.R --taxonomy {input.taxonomy} --input {input.table} --output {output.tdist}
        '''

rule add_labels:
    input:
        table = 'example_data/{dir}/observed_hgt.txt',
        tdist = 'example_data/{dir}/taxonomic_distance.txt'
    output:
        gclass = temp('example_data/{dir}/gene_class.txt'),
        split1 = temp('example_data/{dir}/s1.txt'),
        split2 = temp('example_data/{dir}/s2.txt'),
        table = temp('example_data/{dir}/table_w_labels.txt')
    shell:
        '''
        name=$(echo {input.table} | cut -d "/" -f 2)
        python scripts/make_labels.py --list {input.table} --word $name --header "Gene class" --output {output.gclass}
        cat {input.table} | cut -f 1 > {output.split1}
        cat {input.table} | cut -f 2- > {output.split2}
        paste {output.split1} {output.gclass} {input.tdist} {output.split2} | tail -n +2 > {output.table}
        '''

rule extract_header:
    output:
        temp("example_data/header.txt")
    run:
        head = ""
        head += 'Node' + '\t' + 'Gene class' + '\t' + 'Taxonomic difference' + '\t' + 'Order1' + '\t' + 'Order2' + '\t' + 'Species1' + '\t' + 'Species2' + '\n'

        with open(output[0], 'w') as file:
            file.write(head)

rule combine_hgt_table:
    input:
        table = expand("example_data/{dir}/table_w_labels.txt", dir=DIR),
        header = "example_data/header.txt"
    output:
        table = temp("example_data/hgt_table1.txt")
    shell:
        '''
        cat {input.header} example_data/*/table_w_labels.txt > {output.table}
        '''

rule translate_accession_ids:
    input:
        table = 'example_data/hgt_table1.txt',
        index = '/home/dlund/index_files/ID_index.txt'
    output:
        ids = temp('example_data/assembly_ids.txt'),
        table = temp('example_data/hgt_table2.txt')
    shell:
        '''
        python scripts/translate_accession_ids.py --input {input.table} --index {input.index} --output {output.ids}
        paste {input.table} {output.ids} > {output.table}
        '''

rule lookup_otus:
    input:
        table = 'example_data/hgt_table2.txt'
    output:
        emp = temp('example_data/otus_emp.txt'),
        gwmc = temp('example_data/otus_gwmc.txt'),
        table = temp('example_data/hgt_table3.txt')
    shell:
        '''
        python scripts/lookup_otus.py --input {input.table} --database emp --output {output.emp}
        python scripts/lookup_otus.py --input {input.table} --database gwmc --output {output.gwmc}
        paste {input.table} {output.emp} {output.gwmc} > {output.table}
        '''

rule genome_5mer_distance:
    input:
        'example_data/hgt_table3.txt'
    output:
        temp('example_data/genome_5mer_distance.txt')
    shell:
        '''
        python scripts/genome_5mer_distance.py --input {input} -k 5 --num_cores 1 --output {output}
        '''

rule separate_genome_groups:
    input:
        'example_data/hgt_table3.txt'
    output:
        dname = directory('example_data/separated_groups')
    run:
        table = pd.read_csv(input[0], sep="\t")
        
        subprocess.run('mkdir %s' %(output[0]), shell=True)
        
        for i in range(len(table.iloc[:,0])):
            subdir1 = '-'.join([table.iloc[i,0], table.iloc[i,1], 'grp1'])
            subprocess.run('mkdir %s/%s' %(output[0], subdir1), shell=True)
    
            ids1 = split_ids(table, "Species1")
            ids1.to_csv('%s/%s/accession_ids.txt' %(output[0], subdir1), index=False, header=False)

            subdir2 = '-'.join([table.iloc[i,0], table.iloc[i,1], 'grp2'])
            subprocess.run('mkdir %s/%s' %(output[0], subdir2), shell=True)
    
            ids2 = split_ids(table, "Species2")
            ids1.to_csv('%s/%s/accession_ids.txt' %(output[0], subdir2), index=False, header=False)

rule gene_genome_5mer_distance:
    input:
        dname = 'example_data/separated_groups',
        distr = 'example_data/genome_5mer_distance.txt'
    output:
        temp('example_data/gene_genome_5mer_distance.txt')
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

        subprocess.run('snakemake -s scripts/generate_gene_5mer_distributions.smk --cores 3 all', shell=True)
        subprocess.run('python scripts/gene_genome_5mer_distance.py --input %s --output %s' %(input[1], output[0]))

rule genome_size_difference:
    input:
        true = 'example_data/hgt_table3.txt',
        null = 'example_data/null_table1.txt'
    output:
        true = temp('example_data/genome_size_diff.txt'),
        null = temp('example_data/genome_size_diff_null.txt')
    shell:
        '''
        python scripts/genome_size_difference.py --input {input.true} --output {output.true} --num_cores 1
        python scripts/genome_size_difference.py --input {input.null} --output {output.null} --num_cores 1
        '''

rule cooccurrence:
    input:
        true = 'example_data/hgt_table3.txt',
        null = 'example_data/null_table1.txt'
    output:
        emp_true = temp('example_data/cooccurrence_emp.txt'),
        gwmc_true = temp('example_data/cooccurrence_gwmc.txt'),
        emp_null = temp('example_data/cooccurrence_emp_null.txt'),
        gwmc_null = temp('example_data/cooccurrence_gwmc_null.txt')
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript scripts/cooccurrence.R --database emp --input {input.true} --output {output.emp_true} --num_cores 1
        Rscript scripts/cooccurrence.R --database gwmc --input {input.true} --output {output.gwmc_true} --num_cores 1
        Rscript scripts/cooccurrence.R --database emp --input {input.null} --output {output.emp_null} --num_cores 1
        Rscript scripts/cooccurrence.R --database gwmc --input {input.null} --output {output.gwmc_null} --num_cores 1
        '''

rule gram_stain_difference:
    input:
        true = 'example_data/hgt_table3.txt',
        null = 'example_data/null_table1.txt'
    output:
        true = temp('example_data/gram_stain_diff.txt'),
        null = temp('example_data/gram_stain_diff_null.txt')
    shell:
        '''
        python scripts/genome_size_difference.py --input {input.true} --output {output.true} --num_cores 1
        python scripts/genome_size_difference.py --input {input.null} --output {output.null} --num_cores 1
        '''

rule complete_hgt_table:
    input:
        table = 'example_data/hgt_table3.txt',
        genome_dist = 'example_data/genome_5mer_distance.txt',
        gene_genome_dist = 'example_data/gene_genome_5mer_distance.txt',
        genome_size = 'example_data/genome_size_diff.txt',
        emp = 'example_data/cooccurrence_emp.txt',
        gwmc = 'example_data/cooccurrence_gwmc.txt',
        gram_stain = 'example_data/gram_stain_diff.txt'
    output:
        'observed_horizontal_transfers.txt'
    shell:
        '''
        paste {input.table} {input.genome_dist} {input.gene_genome_dist} {input.genome_size} {input.gram_stain} {input.emp} {input.gwmc} > {output}
        '''

rule create_sample_dict:
    input:
        rules.preprocessing.output.taxonomy
    output:
        'example_data/{dir}/sample_dict.pkl'
    shell:
        '''
        python scripts/create_sample_dict.py --taxonomy {input} --output {output}
        '''

rule generate_null_distribution:
    input:
        expand('example_data/{dir}/sample_dict.pkl', dir=DIR)
    output:
        temp('example_data/null_table1.txt')
    shell:
        '''
        python scripts/generate_null_distribution.py --max_number 1000 --output {output}
        '''

rule complete_null_distribution:
    input:
        table = 'example_data/null_table1.txt',
        genome_size = 'example_data/genome_size_diff_null.txt',
        emp = 'example_data/cooccurrence_emp_null.txt',
        gwmc = 'example_data/cooccurrence_gwmc_null.txt',
        gram_stain = 'example_data/gram_stain_diff_null.txt'
    output:
        'null_distribution.txt'
    shell:
        '''
        paste {input.table} {input.genome_size} {input.gram_stain} {input.emp} {input.gwmc} > {output}
        '''
