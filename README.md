### Introduction

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

Additionally, the following software should be located in your $PATH:
- biom-format >= 2.1.6
- BLAST >= 2.10.1
- FastTree >= 2.1.9
- mafft >= 7.3.10

### Tutorial
Below is a step-by-step guide on how to generate the input data used to train random forest models as described in the paper "Genetic compatibility and ecological connectivity drive the dissemination of antibiotic resistance genes" by Lund et al. 2024.

1. Clone the repository using
    ```
    git clone https://github.com/davidgllund/factors_influencing_HGT_of_ARGs.git
    ```

2. Setup conda environments
    ```
    conda env create -f envs/arg_hgt_setup.yml
    conda env create -f envs/arg_hgt_python.yml
    conda env create -f envs/arg_hgt_R.yml
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

The results of the analysis can then be found in two files:
- observed_horizontal_transfers.txt
- randomized_transfers.txt
