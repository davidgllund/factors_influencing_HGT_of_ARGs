### Introduction

### Dependencies
To run the scripts in this repository, the following software is required:
- Python >= 3.12.1
    - pandas >= 2.1.1
    - Biopython >= 1.83
    - numpy >= 1.26.3
    - progressbar2 >= 4.3.2
    - joblib >= 1.2.0
    - tqdm >= 4.66.1
- Snakemake >= 5.32.0
- R >= 4.2.0
    - Taxonomizr >= 0.10.2
    - tidyverse >=
    - gtools >=
    - phangorn >=
    - Dict >=
    - ggtree >=
    - ape >=
    - data.table >=
    - adephylo >=
    - viridis >=
    - pbapply >=
    - optparse >=
    - caret >=

Additionally, the following software should be located in your $PATH:
- mafft >= 7.3.10
- FastTree >= 2.1.9
- BLAST >= 2.10.1

### Tutorial
Below is a step-by-step guide on how to generate the input data used to train random forest models in ...

1. Clone the repository using
    ```
    git clone https://github.com/davidgllund/factors_influencing_HGT_of_ARGs.git
    ```

2. Enter the folder "example data", where the scripts will be executed
    ```
    cd example_data
    ```
    
3. Identify horizontal transfers, generate null distribution and calculate input features using
    ```
    snakemake -s ../scripts/analyze_hgt.smk --cores [number of cores] all
    ```

The results of the analysis can then be found in ...
