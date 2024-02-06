#!/usr/bin/env Rscript

# This is a version of the script 'Antibiotic_resistance_genes_analyzer.R' from
# https://github.com/marcosparrasmolto/EventsAnalysis/blob/main/Antibiotic_resistance_genes_analyzer.R
# adapted to deal with larger sets of resistance genes.

# Compared with the original script, this version handles only the centroid
# sequences, which allows for the creation of larger phylogenetic trees, and the
# redundancy is instead implemented as a dictionary containing the complete set
# of isolates that carried each unique gene.

# Copyright (c) David Lund 2023.

#-------------------------------------------------------------------------------
# 0 LOAD LIBRARIES
#-------------------------------------------------------------------------------
library(taxonomizr)
library(ggplot2)
library(stringr)
library(gtools)
library(phangorn)
library(Dict)
library(utils)
library(ggnewscale)
library(optparse)
suppressMessages(library(ggtree))
suppressMessages(library(ape))
suppressMessages(library(data.table))
suppressMessages(library(adephylo))
suppressMessages(library(viridis))
suppressMessages(library(tidyr))

#-------------------------------------------------------------------------------
# 1 INPUT ARUGMENTS
#-------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file name", metavar = "character"),
  make_option(c("-t", "--taxonomy"), type = "character", default = NULL,
              help = "Host taxonomy table", metavar = "character"),
  make_option(c("-d", "--directory"), type = "character", default = NULL,
              help = "Directory name", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file name", metavar = "character")
)
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$taxonomy) | is.null(opt$directory) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("ERROR: Missing input argument(s)", call. = FALSE)
}

#-------------------------------------------------------------------------------
# 2 DEFINE FUNCTIONS
#-------------------------------------------------------------------------------
select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree) + 1, tree$tip.label[element], tree$node.label[element - Ntip(tree)])
}

resolve.multilevel.event <- function(leaf_ids, leaf_dict, taxa) {
  leaf_taxonomy <- taxonomy_table[leaf_dict[leaf_ids], ]
  leaf_taxonomy <- na.omit(leaf_taxonomy)

  order_combinations <- expand.grid(x = taxa, y = taxa)
  order_combinations <- order_combinations[order_combinations$x != order_combinations$y, ]
  order_combinations$x <- as.character(order_combinations$x)
  order_combinations$y <- as.character(order_combinations$y)

  included <- c()
  output <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(output) <- c("Node", "Order1", "Order2", "Species1", "Species2")

  for (l in 1:nrow(order_combinations)) {
    if (paste(order_combinations$x[l], order_combinations$y[l]) %in% included || paste(order_combinations$y[l], order_combinations$x[l]) %in% included) {
      next
    }
    else{
      output[nrow(output) + 1, ] <- c(paste(c("Leaf", grep(leaf_ids, leaf_dict$keys)), collapse = ""), order_combinations$x[l], order_combinations$y[l], paste(rownames(leaf_taxonomy[leaf_taxonomy$order == order_combinations$x[l], ]), collapse = ";"), paste(rownames(leaf_taxonomy[leaf_taxonomy$order == order_combinations$y[l], ]), collapse = ";"))
      included <- c(included, paste(order_combinations$x[l], order_combinations$y[l]))
    }
  }

  return(output)
}

get.phylum <- function(query) {
  included <- c("Actinomycetota", "Bacteroidota", "Bacillota", "Campylobacterota", "Pseudomonadota")
  retrieved <- unique(na.omit(taxonomy_table[leaves[query], "phylum"]))

  if (length(retrieved) == 0) {
    phylum <- NA
  }

  else if (length(retrieved) > 1) {
    phylum <- "Multiple"
  }

  else if (length(retrieved) == 1 & retrieved %in% included) {
    phylum <- retrieved
  }

  else if (length(retrieved) == 1 & !(retrieved %in% included)) {
    phylum <- "Other"
  }

  else {
    phylum <- NA
  }

  return(phylum)
}

#-------------------------------------------------------------------------------
# 3 PREPROCESSING
#-------------------------------------------------------------------------------
setwd(opt$directory)

taxonomy_table <- read.table(opt$taxonomy, sep = "\t", stringsAsFactors = FALSE)

tree <- read.tree(opt$input)
tree <- makeNodeLabel(tree)

leaves <- Dict$new(
  token = 0,
  .overwite = TRUE
)

keys <- data.frame(V1 = system(paste("ls clusters"), intern = TRUE))
leaves[as.character(keys$V1[1])] = system(paste('cat clusters/', as.character(keys$V1[1]), '/hidden.txt', sep = ''), intern = TRUE)

print("Preparing leaf dictionary")
pb <- txtProgressBar(min = 0, max = length(keys$V1), initial = 0, style = 3)

for (i in 1:length(keys$V1)) {
  leaves[as.character(keys$V1[i])] = system(paste('cat clusters/', as.character(keys$V1[i]), '/hidden.txt', sep = ''), intern = TRUE)
  setTxtProgressBar(pb, i)
}

close(pb)

edge_table <- data.frame(
  "parent_ID" = tree$edge[, 1],
  "parent_label" = sapply(tree$edge[, 1], select.tip.or.node, tree = tree),
  "child_ID" = tree$edge[, 2],
  "child_label" = sapply(tree$edge[, 2], select.tip.or.node, tree = tree)
)

edge_table$parent_ID <- as.numeric(edge_table$parent_ID)
edge_table$child_ID <- as.numeric(edge_table$child_ID)
edge_table$parent_label <- as.character(edge_table$parent_label)
edge_table$child_label <- as.character(edge_table$child_label)

node_list <- unique(edge_table$parent_label)

results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) <- c("Node", "Order1", "Order2", "Species1", "Species2")
event_nodes <- c()
nodes_to_color <- c()

#-------------------------------------------------------------------------------
# 4 SEARCH FOR HORIZONTAL TRANSFERS
#-------------------------------------------------------------------------------
# This loop will traverse the phylogenetic tree, examining each position in turn
# looking for horizontally transferred genes. To classify a horizontal transfer,
# the following scenarios are possible:

# 1) A single leaf is comprised of multiple genes (100% amino acid identity)
#    carried by species from multiple different taxonomic orders.

# 2) The descendents of a node in the tree represent genes carried by species
#    from two different taxonomic orders.

print("Searching for horizontal transfers")
start_time <- Sys.time()
pb <- txtProgressBar(min = 0, max = length(node_list), initial = 0, style = 3)

for (i in 1:length(node_list)) {
  leaf_ids <- edge_table$child_label[edge_table$child_ID %in% Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]]]

  all_leaves <- c()
  for (k in 1:length(leaf_ids)) {
    if (leaf_ids[k] %in% keys$V1) {
      all_leaves <- c(all_leaves, leaves[leaf_ids[k]])
    }
    else {
      all_leaves <- c(all_leaves, leaf_ids[k])
    }
  }

  # Check if node is furthest down in the tree (only two leaves as descendants)
  if (!sort(str_detect(all_leaves, "Node"), decreasing = TRUE)[1]) {
    child_node_levels <- levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% all_leaves])))

    # Find the orders represented in the two descendant leaves
    orders1 <- levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% leaves[leaf_ids[1]]])))
    orders2 <- levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% leaves[leaf_ids[2]]])))

    # Check if either leaf includes more than one order
    if (length(orders1) > 1 || length(orders2) > 1) {

      # If the first leaf includes more than one order, generate all possible
      # transfer events
      if (length(orders1) > 1) {
        possible_events <- resolve.multilevel.event(leaf_ids[1], leaves, orders1)

        for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j, ]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][1])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
      }

      # If the second leaf includes more than one order, generate all possible
      # transfer events
      if (length(orders2) > 1) {
        possible_events <- resolve.multilevel.event(leaf_ids[2], leaves, orders2)

        for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j, ]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][2])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
      }
    }

    # If each leaf represents a different order (a single order each), record
    # as horizontal transfer
    else if((length(child_node_levels) == 2 && !('Unknown' %in% child_node_levels)) && length(orders1) == 1 && length(orders2) == 1) {
      results[nrow(results) +1, ] <- c(node_list[length(node_list)-i+1], paste(levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% leaves[leaf_ids[1]]]))), collapse = ";"), paste(levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% leaves[leaf_ids[2]]]))), collapse = ";"), paste(leaves[leaf_ids[1]], collapse = ";"),paste(leaves[leaf_ids[2]], collapse = ";"))
      event_nodes <- c(event_nodes, unique(edge_table[edge_table$parent_label %in% results$Node[nrow(results)], 1]))
    }
  }

  else {
    descendant_nodes = Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]]

    # Given that all direct descendants are nodes, skip if a recorded event
    # exists further down in the tree (to avoid nested events)
    if ((sum(event_nodes %in% Descendants(tree, descendant_nodes[1], type = "all")) > 0 || sum(event_nodes %in% Descendants(tree, descendant_nodes[2], type = "all")) > 0) && sum(str_detect(leaf_ids, "Node")) == 2) {
      next
    }

    else {
      leaf_ids1 <- edge_table[edge_table$child_ID %in% Descendants(tree, descendant_nodes[1])[[1]], 4]
      leaf_ids2 <- edge_table[edge_table$child_ID %in% Descendants(tree, descendant_nodes[2])[[1]], 4]

      vector1 <- NULL
      vector2 <- NULL
    
      for(o in 1:length(leaf_ids1)) {
        vector1 <- c(vector1, leaves[leaf_ids1[o]])
      }

      for(o in 1:length(leaf_ids2)) {
        vector2 <- c(vector2, leaves[leaf_ids2[o]])
      }

      vector1 <- sort(vector1)
      vector2 <- sort(vector2)

      orders1 <- levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% vector1])))
      orders2 <- levels(factor(as.character(taxonomy_table$order[rownames(taxonomy_table) %in% vector2])))

      # Skip if node is tripartite
      if (length(descendant_nodes) > 2) {
        next
      }

      # Skip if the two descendant groups represent the same order (one
      # order each)
      else if (length(unique(c(orders1, orders2))) == 1) {
        next
      }
      
      # If the first descendant is a single leaf that includes more than one
      # order while the second descendant is a node representing a single order,
      # generate all possible transfer events from the leaf
      else if (length(orders1) >= 2 && length(orders2) <= 1 && length(Descendants(tree, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][1], type = "tips")[[1]]) == 1) {
        possible_events <- resolve.multilevel.event(leaf_ids[1], leaves, orders1)

        for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j, ]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][1])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
      }

      # If the second descendant is a single leaf that includes more than one
      # orderwhile the first descendant is a node representing a single order,
      # generate all possible HGT events from this leaf
      else if (length(orders2) >= 2 && length(orders1) <= 1 && length(Descendants(tree, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][2], type = "tips")[[1]]) == 1) {
        possible_events <- resolve.multilevel.event(leaf_ids[2], leaves, orders2)

        for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j,]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][2])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
      }

      # Check if both descendant groups represent more than one order
      else if (length(orders1) >= 2 && length(orders2) >= 2) {
        # If the first descendant is a single leaf, generate all possible
        # HGT events from this leaf
        if (length(Descendants(tree, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][1], type = "tips")[[1]]) == 1) {
          possible_events <- resolve.multilevel.event(leaf_ids[1], leaves, orders1)

          for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j, ]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][1])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
        }

        # If the second descendant is a single leaf, generate all possible
        # HGT events from this leaf
        else if (length(Descendants(tree, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[2]][2], type = "tips")[[1]]) == 1) {
          possible_events <- resolve.multilevel.event(leaf_ids[2], leaves, orders2)

          for (j in 1:nrow(possible_events)) {
            results[nrow(results) + 1, ] <- possible_events[j, ]
            event_nodes <- c(event_nodes, Descendants(tree, edge_table[edge_table$parent_label == node_list[length(node_list) - i + 1], 1], type = "child")[[1]][2])
            nodes_to_color <- c(nodes_to_color, event_nodes[length(event_nodes)])
          }
        }
      }

      # If both descendants are nodes representing different orders (one
      # order each), record as HGT event
      else if (length(descendant_nodes) == 2 && orders1 != "Unknown" && orders2 != "Unknown" && length(orders1) == 1 && length(orders2) == 1 && orders1 != orders2) {
        results[nrow(results) + 1, ] <- c(node_list[length(node_list) - i + 1], paste(orders1, collapse = ";"), paste(orders2, collapse = ";"), paste(vector1, collapse = ";"), paste(vector2, collapse = ";"))
        event_nodes <- c(event_nodes, unique(edge_table[edge_table$parent_label %in% results$Node[nrow(results)], 1]))
      }
    }
  }
  setTxtProgressBar(pb, i)
}

end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, unit = "mins")
close(pb)
print(paste("Finished searching, time elapsed:", time_diff, "min"))

#-------------------------------------------------------------------------------
# 5 POST-PROCESSING AND EXPORTING RESULTS
#-------------------------------------------------------------------------------
# Remove any nested transfers that were previously missed
results_adj <- data.frame(matrix(ncol = 5))
colnames(results_adj) <- c("Node", "Order1", "Order2", "Species1", "Species2")
k <- 1

for (i in 1:length(event_nodes)) {
  if (length(grep("Leaf", results$Node[i])) != 0) {
    results_adj[k, ] <- results[i, ]
    k <- k + 1
    next
  }
  else if (sum(unique(event_nodes) %in% Descendants(tree, event_nodes[i], type = "all")) > 0) {
    next
  }
  else{
    results_adj[k, ] <- results[i, ]
    k <- k + 1
    next
  }
}

print(paste("A total of", nrow(results_adj), "HGT events were identified"))
write.table(results_adj, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------------------------------------------------------------
# 6 PLOT TREE
#-------------------------------------------------------------------------------
nodes_to_color <- unique(nodes_to_color)
nodes_to_color <- c(nodes_to_color, unique(edge_table[edge_table$parent_label %in% results_adj$Node,1]))

headers <- as.character(tree$tip.label)
phylum_list <- sapply(headers, get.phylum)
phylum_list <- as.data.frame(phylum_list)

blastout <- read.table("blastout.txt", sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = "")
blast_df <- data.frame(perc_id = blastout$V3)
rownames(blast_df) <- blastout$V1

blast_df <- as.data.frame(blast_df[headers, ])
rownames(blast_df) <- headers

p1 <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tiplab(size = 0) +
  geom_point2(aes(subset = (node %in% nodes_to_color)), color = "green", size = 1.5) +
  geom_nodelab2(aes(subset = (node %in% nodes_to_color)), color = "green", size = 1.5)

p2 <- gheatmap(p1, phylum_list, offset = -9.5, width = 0.2, color = NULL) +
  scale_fill_manual(values=c("#4daf4a", "#377eb8", "#ff7f00", "#984ea3", "#999999", "#e41a1c"), na.value="white", name = "Phylum")

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, blast_df, offset = 14, width = 0.2, color = NULL) +
  scale_fill_viridis_c(direction = -1)

pdf("tree_annotated.pdf", height = 10, width = 10)
plot(p3)
dev.off()

#-------------------------------------------------------------------------------
