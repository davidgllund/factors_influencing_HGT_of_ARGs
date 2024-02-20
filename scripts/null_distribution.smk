#-------------------------------------------------------------------------------
# 0 IMPORT LIBRARIES
#-------------------------------------------------------------------------------
import glob

DIR, = glob_wildcards('example_data/{dir}/host_taxonomy.txt')

#-------------------------------------------------------------------------------
# 1 RULES
#-------------------------------------------------------------------------------
rule all:
    input:
        'example_data/null_table.txt'

rule create_sample_dict:
    input:
        'example_data/{dir}/host_taxonomy.txt'
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
        temp('example_data/null_table.txt')
    shell:
        '''
        python scripts/generate_null_distribution.py --max_number 1000 --output {output}
        '''
