#!/bin/bash
while getopts "p:" opt
do
   case "$opt" in
      p ) p="$OPTARG" ;;
   esac
done

# GENERATE POSITIVE DATASET
# Collect data on horizontal gene transfers
snakemake -s scripts/hgt_table.smk --cores $p --use-conda --conda-frontend conda all

python scripts/translate_accession_ids.py --input example_data/hgt_table.txt --index /home/dlund/index_files/ID_index.txt --output example_data/assembly_ids.txt
paste example_data/hgt_table.txt example_data/assembly_ids.txt > example_data/hgt_table2.txt
rm example_data/assembly_ids.txt

python scripts/lookup_otus.py --input example_data/hgt_table2.txt --database emp --output example_data/otus_emp.txt
python scripts/lookup_otus.py --input example_data/hgt_table2.txt --database gwmc --output example_data/otus_gwmc.txt
paste example_data/hgt_table2.txt example_data/otus_emp.txt example_data/otus_gwmc.txt > example_data/hgt_table3.txt
rm example_data/otus_* example_data/hgt_table2.txt

# Calculate features of observed transfers
python scripts/genome_5mer_distance.py --input example_data/hgt_table3.txt -k 5 --num_cores $p --output example_data/genome_5mer_data.txt
cat example_data/genome_5mer_data.txt | cut -f 3 > example_data/genome_5mer_distance.txt

python scripts/separate_genome_groups.py --input example_data/hgt_table3.txt --output example_data/separated_groups

snakemake -s scripts/generate_gene_5mer_distributions.smk --cores $p all
python scripts/gene_genome_5mer_distance.py --input example_data/genome_5mer_data.txt --output example_data/gene_genome_5mer_distance.txt

python scripts/genome_size_difference.py --input example_data/hgt_table3.txt --output example_data/genome_size_diff.txt --num_cores $p

Rscript scripts/cooccurrence.R --database emp --input example_data/hgt_table3.txt --output example_data/cooccurrence_emp.txt --num_cores $p
Rscript scripts/cooccurrence.R --database gwmc --input example_data/hgt_table3.txt --output example_data/cooccurrence_gwmc.txt --num_cores $p

python scripts/gram_stain_difference.py --input example_data/hgt_table3.txt --output example_data/gram_stain_diff.txt

paste example_data/hgt_table3.txt example_data/genome_5mer_distance.txt example_data/gene_genome_5mer_distance.txt example_data/genome_size_diff.txt example_data/gram_stain_diff.txt example_data/cooccurrence_emp.txt example_data/cooccurrence_gwmc.txt > observed_horizontal_transfers.txt
rm example_data/hgt_table3.txt example_data/genome_5mer_distance.txt example_data/gene_genome_5mer_distance.txt example_data/genome_size_diff.txt example_data/gram_stain_diff.txt example_data/cooccurrence_emp.txt example_data/cooccurrence_gwmc.txt example_data/hgt_table.txt

# GENERATE NEGATIVE DATASET
snakemake -s scripts/null_distribution.smk --cores $p --use-conda --conda-frontend conda all

python scripts/genome_size_difference.py --input example_data/null_table.txt --output example_data/genome_size_diff_null.txt --num_cores $p

Rscript scripts/cooccurrence.R --database emp --input example_data/null_table.txt --output example_data/cooccurrence_emp_null.txt --num_cores $p
Rscript scripts/cooccurrence.R --database gwmc --input example_data/null_table.txt --output example_data/cooccurrence_gwmc_null.txt --num_cores $p

python scripts/gram_stain_difference.py --input example_data/null_table.txt --output example_data/gram_stain_diff_null.txt --num_cores $p

paste example_data/null_table.txt example_data/genome_size_diff_null.txt example_data/gram_stain_diff_null.txt example_data/cooccurrence_emp_null.txt example_data/cooccurrence_gwmc_null.txt > null_distribution.txt
rm example_data/null_table.txt example_data/genome_size_diff_null.txt example_data/gram_stain_diff_null.txt example_data/cooccurrence_emp_null.txt example_data/cooccurrence_gwmc_null.txt