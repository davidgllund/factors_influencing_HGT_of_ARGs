### Introduction
This repository contains scripts and files needed to generate the input data used to train random forest models as described in the paper "Genetic compatibility and ecological connectivity drive the dissemination of antibiotic resistance genes" by Lund et al. 2024. Alongside the scripts, example data representing predicted ARGs from two gene classes is also provided, alongside .yml files with which to create conda environments with the required software. 

### Dependencies
To run the scripts in this repository, the following software is required:
- Python >= 3.12.1
    - biopython >= 1.83
    - joblib >= 1.2.0
    - numpy >= 1.26.3
    - pandas >= 2.1.1
    - progressbar2 >= 4.3.2
    - tqdm >= 4.66.1
- Snakemake >= 5.32.0
- R >= 4.2.0
    - ape >= 5.7.1
    - caret >= 6.0_94
    - data.table >= 1.14.10
    - ggtree >= 3.10.0
    - gtools >= 3.9.4
    - optparse >= 1.7.4
    - pbapply >= 1.7_2
    - phangorn >= 2.11.1
    - taxonomizr >= 0.10.2
    - tidyverse >= 2.0.0
    - viridis >= 0.6.5
- biom-format >= 2.1.6
- BLAST >= 2.10.1
- FastTree >= 2.1.9
- mafft >= 7.3.10

Files that can be used to create conda environments with the required software installed can be found in /envs. The installation and setup of the conda environments is expected to take around 30 minutes on a "normal" computer.

### Tutorial
Below is a step-by-step guide on how to generate the data. Please note that the pipeline is designed to run on Unix based servers. Additionally, please note that running the example code will result in the downloading and generation of several large files. As such, please make sure that you have at least 12GB of free disk space when executing the code to avoid any issues. Finally, please note that this is a large and demanding analysis, which is expected to take around four hours to complete, depending on the number of cores used.

1. Clone the repository using
    ```
    git clone https://github.com/davidgllund/factors_influencing_HGT_of_ARGs.git
    ```

2. Setup conda environments
    ```
    cd factors_influencing_HGT_of_ARGs
    conda env create -f envs/arg_hgt_setup.yml
    conda env create -f envs/arg_hgt_python.yml
    ```
3. Download bacterial genomes and metagenomes, and setup auxiliary files
   ```
   conda activate arg_hgt_setup
   python scripts/setup_and_download.py
   ```
    
4. Identify horizontal transfers, generate null distribution and calculate input features using
    ```
    conda activate arg_hgt_python
    bash scripts/analyze_hgt.sh -p [number of cores]
    ```

#### Output files
After running the example, the results of the analysis can be found in two files:
- **observed_horizontal_transfers.txt**
- **randomized_transfers.txt**

These two tables represent the positive and negative datasets respectively, where each row represents an observation. For each observation, the folowing features are recorded:
- **Node:** Position in the phylogenetic tree at which the transfer was observed.
- **Gene class:** Class of resistance gene which the observation represents (nomenclature taken from fARGene [1]).
- **Antibitoic class:** Class of antibiotics which the corresponding ARG confers resistance to.
- **Taxonomic difference:** Level of taxonomic difference between the two orders included in the observation (order, class, or phylum).
- **Order1:** The first bacterial order involved in the observation.
- **Order2:** The second bacterial order involved in the observation.
- **Species1:** Species from the first bacterial order involved in the observation.
- **Species2:** Species from the second bacterial order involved in the observation.
- **AssemblyAccession1:** Assembly accessions of the geomes from the first bacterial order involved in the observation.
- **AssemblyAccession2:** Assembly accessions of the geomes from the second bacterial order involved in the observation.
- **Matched_OTU_EMP1:** OTU(s) from the Earth Microbiome Project mapping to the genomes from order 1.
- **Matched_OTU_EMP2:** OTU(s) from the Earth Microbiome Project mapping to the genomes from order 2.
- **Matched_OTU_GWMC1:** OTU(s) from the Global Water Microbiome Consortium matching to the genomes from order 1.
- **Matched_OTU_GWMC2:** OTU(s) from the Global Water Microbiome Consortium matching to the genomes from order 2.
- **Genome_5mer_distance:** Euclidean distance between the mean 5mer distributions of the genomes from different orders.
- **Gene_genome_5mer_distance:** Maximial observed euclidean distance between a the 5mer distribution of a random ARG included in the observation and the mean 5mer distribution of the genomes from order 1 and order 2.
- **Genome_size_difference:** Proportional difference in mean size between the genomes from different orders.
- **Gram_stain_difference:** Gram staining properties of order 1 and order 2. Encoded as NN (both gram negative), PP (both gram positive), or NP (different gram staining).
- **Animal:** Estimated co-occurrence of the genomes from different orders in animal microbiomes.
- **Human:** Estimated co-occurrence of the genomes from different orders in human microbiomes.
- **Soil:** Estimated co-occurrence of the genomes from different orders in soil microbiomes.
- **Water:** Estimated co-occurrence of the genomes from different orders in water microbiomes.
- **Wastewater:** Estimated co-occurrence of the genomes from different orders in wastewater microbiomes.

#### References
1. Berglund F, Österlund T, Boulund F, Marathe NP, Larsson DJ, Kristiansson E. Identification and reconstruction of novel antibiotic resistance genes from metagenomes. Microbiome. 2019 Dec;7:1-4.
