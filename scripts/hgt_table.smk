#-------------------------------------------------------------------------------
# 0 IMPORT LIBRARIES
#-------------------------------------------------------------------------------
import glob

DIR, = glob_wildcards('example_data/{dir}/predicted-orfs-amino.fasta')

#-------------------------------------------------------------------------------
# 1 RULES
#-------------------------------------------------------------------------------
rule all:
    input:
        'example_data/hgt_table.txt'

rule preprocessing:
    input:
        nucleotides = 'example_data/{dir}/predicted-orfs.fasta',
        proteins = 'example_data/{dir}/predicted-orfs-amino.fasta',
        taxonomy = 'index_files/assembly_taxonomy.txt'
    output:
        taxonomy = 'example_data/{dir}/host_taxonomy.txt',
        fasta_w_species = temp('example_data/{dir}/args_w_species.fasta'),
        clusters = directory('example_data/{dir}/clusters'),
        centroids = temp('example_data/{dir}/args_clustered.fasta'),
        alignment = temp('example_data/{dir}/alignment_clustered.aln'),
        tree = temp('example_data/{dir}/tree_clustered.txt'),
        blastout = temp('example_data/{dir}/blastout.txt'),
        headerlines = 'example_data/{dir}/fasta_headers.txt'
    shell:
        '''
        python scripts/preprocessing.py --nucleotides {input.nucleotides} --proteins {input.proteins} --taxonomy_index {input.taxonomy} --taxonomy {output.taxonomy} --fasta_w_species {output.fasta_w_species} --clusters {output.clusters} --centroids {output.centroids} --alignment {output.alignment} --tree {output.tree} --blastout {output.blastout}
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
        table = "example_data/hgt_table.txt"
    shell:
        '''
        cat {input.header} example_data/*/table_w_labels.txt > {output.table}
        '''
