rule all:
    input:
        config["output_emp"],
        config["output_gwmc"]

rule cooccurrence_emp:
    input:
        config["input"]
    output:
        config["output_emp"]
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript scripts/cooccurrence.R --database emp --input {input} --output {output} --num_cores 1
        '''

rule cooccurrence_gwmc:
    input:
        config["input"]
    output:
        config["output_gwmc"]
    conda:
        '../envs/arg_hgt_R.yml'
    shell:
        '''
        Rscript scripts/cooccurrence.R --database gwmc --input {input} --output {output} --num_cores 1
        '''