
#!/usr/bin/env Rscript

# This script converts co-occurrence data from the Earth Microbiome project into
# a more useable format. The initial data, consisting of a large count matrix
# where each row represents an OTU and each column represents a sample, is
# converted into five files, each representing a type of environment (Human,
# Animal, Soil, Water, Other). Each file contains all available pairs of
# OTUs from the original data, and the level to which those OTUs co-occurr in
# that kind of environment.

# Copyright (c) David Lund 2023.

#-------------------------------------------------------------------------------
# 0 LOAD LIBRARIES
#-------------------------------------------------------------------------------
library(data.table)
suppressMessages(library(tidyverse))
library(utils)
library(pbapply)
library(optparse)
library(dplyr)

pbo <- pboptions(type = "timer")

#-------------------------------------------------------------------------------
# 1 INPUT ARUGMENTS
#-------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--database"), type = "character", default = NULL,
              help = "Database to use ('emp' or 'gwmc')", metavar = "character"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file name", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file name", metavar = "character"),
  make_option(c("--num_cores"), type = "integer", default = 2,
              help = "Number of cores to use for parallel processing", metavar = "number")
)
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$database) | is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("ERROR: Missing input argument(s)", call. = FALSE)
}

#-------------------------------------------------------------------------------
# 2 READ AND FILTER DATA
#-------------------------------------------------------------------------------
if (opt$database == "emp") {
  source("/home/dlund/HGT_inference_project/code_for_publication/scripts/estimate_cooccurrence_emp.R")

  count_table <- data.frame(fread("/home/dlund/HGT_inference_project/analysis/coocurrence_data_emp/otus_gg_13_8.txt"))

  metadata <- data.frame(fread("/home/dlund/HGT_inference_project/analysis/coocurrence_data_emp/emp_qiime_mapping_qc_filtered.tsv"))
  metadata$sample_scientific_name[metadata$sample_scientific_name == "skin metagenome" & metadata$host_common_name_provided == "human"] <- "human skin metagenome"
  metadata$sample_scientific_name[metadata$sample_scientific_name == "upper respitatory tract metagenome" & metadata$host_common_name_provided == "human"] <- "human upper respiratory tract metagenome"

  classification_list <- list(
  "activated carbon metagenome" = "Others", "air metagenome" = "Others",
  "algae metagenome" = "Water", "bat metagenome" = "Animal",
  "biofilm metagenome" = "Others", "bovine gut metagenome" = "Animal",
  "coral metagenome" = "Water", "ecological metagenomes" = "Others",
  "estuary metagenome" = "Water", "fish metagenome" = "Animal",
  "freshwater metagenome" = "Water", "freshwater sediment metagenome" = "Water",
  "groundwater metagenome" = "Water", "gut metagenome" = "Animal",
  "human gut metagenome" = "Human", "human metagenome" = "Human",
  "human nasal/pharyngeal metagenome" = "Human",
  "human oral metagenome" = "Human", "human skin metagenome" = "Human",
  "hydrocarbon metagenome" = "Others", "hydrothermal vent metagenome" = "Water",
  "hypersaline lake metagenome" = "Water", "indoor metagenome" = "Others",
  "insect metagenome" = "Animal", "marine metagenome" = "Water",
  "marine sediment metagenome" = "Water", "metagenomes" = "Others",
  "microbial mat metagenome" = "Others", "mine drainage metagenome" = "Soil",
  "organismal metagenomes" = "Others", "permafrost metagenome" = "Soil",
  "plant metagenome" = "Soil", "primate metagenome" = "Animal",
  "rat metagenome" = "Animal", "rhizosphere metagenome" = "Soil",
  "root metagenome" = "Soil", "seawater metagenome" = "Water",
  "sediment metagenome" = "Soil", "skin metagenome" = "Animal",
  "soil metagenome" = "Soil", "Thalassia testudinum" = "Water",
  "upper respiratory tract metagenome" = "Animal",
  "human upper respiratory tract metagenome" = "Human")

  environment_categories <- c("Animal", "Human", "Water", "Soil")

  samples_from_env <- list()

  for (i in seq_len(length(environment_categories))) {
    samples_from_env[[environment_categories[i]]] <- metadata[metadata$sample_scientific_name %in% names(classification_list[classification_list == environment_categories[i]]), 1]
  }

} else if (opt$database == "gwmc") {
  source("/home/dlund/HGT_inference_project/code_for_publication/scripts/estimate_cooccurrence_gwmc.R")

  count_table <- data.frame(fread("/home/dlund/HGT_inference_project/analysis/coocurrence_data_gwmc/GWMC_16S_otutab.txt"))
}

horizontal_transfers <- data.frame(fread(opt$input))

rownames(count_table) <- count_table[, 1]
count_table <- count_table[, c(-1, -ncol(count_table))]
count_table <- count_table[, !(apply(count_table, 2, sum) < 10000)]

count_table[count_table < 3] <- 0
count_table[count_table >= 3] <- 1

#-------------------------------------------------------------------------------
# 3 ESTIMATE COOCCURRENCE
#-------------------------------------------------------------------------------
print("Calculating co-occurrence")
start_time <- Sys.time()

if (opt$database == "emp") {
  otus <- horizontal_transfers[,c("Matched_OTU_EMP1", "Matched_OTU_EMP2")]
  estimated_cooccurrence <- pbapply(otus, 1, analyze.cooccurrence.emp, cl = opt$num_cores) %>% t() %>% as.data.frame()
  colnames(estimated_cooccurrence) <- c("Animal", "Human", "Soil", "Water")

} else if (opt$database == "gwmc") {
  otus <- horizontal_transfers[,c("Matched_OTU_GWMC1", "Matched_OTU_GWMC2")]
  estimated_cooccurrence <- pbapply(otus, 1, analyze.cooccurrence.gwmc, cl = opt$num_cores) %>% as.data.frame()
  colnames(estimated_cooccurrence) <- c("Wastewater")
}

end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "hours")
print(paste("Finished measusing cooccurrence, time elapased:", time_diff, "hours"))

write.table(estimated_cooccurrence, opt$output, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")