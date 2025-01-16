#----------------- Supplementary Fig 2 ----------------------------------------
# Load packages
library(ggplot2)

# Plotting settings
fig_theme <-  theme(axis.text.x=element_blank(),
                    axis.text.y=element_text(size = rel(1.3)),
                    axis.title=element_text(size=rel(1.3)),
                    legend.title=element_text(size=rel(1.2)), 
                    legend.text=element_text(size=rel(1)),
                    plot.title = element_text(size = rel(3)))

mechanism_palette <- c(
  "AAC"='#a50026',
  "APH" ="#d73027",
  "Class A, C, D"='#fdae61',
  "Class B"="#fee090",
  "Erm" = "#ffffff",
  "Mph" = "#e0f3f8",
  "Tet efflux" = "#74add1",
  "Tet enzyme" = "#4575b4",
  "Tet RPG" = "#313695",
  "Qnr" = "#000000"
)

# Compile matrices
aac <- matrix(c(65094, 1101566, 571094, 928248), nrow=2, ncol=2, byrow = TRUE)
aph <- matrix(c(387203, 779457, 67231, 1432111), nrow=2, ncol=2, byrow = TRUE)
class_A_C_D <- matrix(c(287169, 879491, 261245, 1238097), nrow=2, ncol=2, byrow = TRUE)
class_B <- matrix(c(5926, 1160734, 468248, 1031094), nrow=2, ncol=2, byrow = TRUE)
erm <- matrix(c(43686, 1122974, 14453, 1484889), nrow=2, ncol=2, byrow = TRUE)
mph <- matrix(c(56203, 1110457, 6721, 1492621), nrow=2, ncol=2, byrow = TRUE)
tet_efflux <- matrix(c(234174, 932486, 28278, 1471064), nrow=2, ncol=2, byrow = TRUE)
tet_enzyme <- matrix(c(1671, 1164989, 541, 1498801), nrow=2, ncol=2, byrow = TRUE)
tet_rpg <- matrix(c(47694, 1118966, 65391, 1433951), nrow=2, ncol=2, byrow = TRUE)
qnr <- matrix(c(37840, 1128820, 16140, 1483202), nrow=2, ncol=2, byrow = TRUE)

# Analysis
ratio <- c(fisher.test(aac)$estimate,
           fisher.test(aph)$estimate,
           fisher.test(class_A_C_D)$estimate,
           fisher.test(class_B)$estimate,
           fisher.test(erm)$estimate,
           fisher.test(mph)$estimate,
           fisher.test(tet_efflux)$estimate,
           fisher.test(tet_enzyme)$estimate,
           fisher.test(tet_rpg)$estimate,
           fisher.test(qnr)$estimate)

pval <- c(fisher.test(aac)$p.value,
          fisher.test(aph)$p.value,
          fisher.test(class_A_C_D)$p.value,
          fisher.test(class_B)$p.value,
          fisher.test(erm)$p.value,
          fisher.test(mph)$p.value,
          fisher.test(tet_efflux)$p.value,
          fisher.test(tet_enzyme)$p.value,
          fisher.test(tet_rpg)$p.value,
          fisher.test(qnr)$p.value)

# Compile results
significance <- rep('', length(pval))

for (i in 1:length(pval)){
  if (pval[i] < 0.01){
    significance[i] <- '*'
  }
}

mechanism <- c("AAC", "APH", "Class A, C, D", "Class B", "Erm", "Mph", "Tet efflux", "Tet enzyme", "Tet RPG", "Qnr")
data <- data.frame(log(ratio), pval, mechanism, significance)

fig_s2 <- ggplot(data, aes(y=log.ratio., 
                           x=factor(mechanism, levels = mechanism),
                           fill = factor(mechanism, levels = mechanism))) + 
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_text(aes(label=significance), position=position_dodge(width=0.9), vjust=-0.25, size=7) +
  xlab("") +
  ylab("log(Odds ratio)") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_palette, name = "Resistance mechanism") +
  fig_theme

# Grouped
pdf("sup_fig_2.pdf", height = 6, width=7)
plot(fig_s2)
dev.off()

#----------------- Supplementary Fig 3 ----------------------------------------
# Load packages
library(ggplot2)

# Plotting settings
fig_theme <- theme(axis.text.x=element_text(size = rel(1.8)),
                   axis.text.y=element_text(size = rel(1.5)),
                   axis.title=element_text(size=rel(1.5)),
                   legend.title=element_text(size=rel(1.4)), 
                   legend.text=element_text(size=rel(1.2)))

mechanism_palette <- c(
  "AAC"='#a50026',
  "APH" ="#d73027",
  "Class A, C, D"='#fdae61',
  "Class B"="#fee090",
  "Erm" = "#ffffff",
  "Mph" = "#e0f3f8",
  "Tet efflux" = "#74add1",
  "Tet enzyme" = "#4575b4",
  "Tet RPG" = "#313695",
  "Qnr" = "#000000"
)

# Compile data
df <- read.table("data_for_figures/sequence_similarity.txt", sep="\t")
colnames(df) <- c("class", "similarity")
df$mechanism <- rep(NA, rep=nrow(df))

mech_names <- c("aac", "aph", "class_A_C_D", "class_B", "erm", "mph", "tet_efflux", "tet_enzyme", "tet_rpg", "qnr")
mech_names_adj <- c("AAC", "APH", "Class A, C, D", "Class B", "Erm", "Mph", "Tet efflux", "Tet enzyme", "Tet RPG", "Qnr")

for (i in 1:length(mech_names)) {
  if (mech_names[i] == "class_A_C_D") {
    df$mechanism[df$class %in% c("class_A", "class_C", "class_D_1", "class_D_2")] <- mech_names_adj[i]
  }
  else {
    df$mechanism[grep(mech_names[i], df$class)] <- mech_names_adj[i]
  }
}

bins <- c(-Inf, 50, 60, 70, 80, 90, 95, 99, Inf)
labels <- c("<50", "50-60", "60-70", "70-80", "80-90", "90-95", "95-99", "99<")

df <- df %>%
  mutate(category = cut(similarity, breaks = bins, labels = labels, right = FALSE))

df2 <- df %>%
  group_by(mechanism, category) %>%
  summarise(init = n(), .groups = "drop") %>%
  group_by(mechanism) %>%
  mutate(init = init / sum(init))

# Aggregate data
df2 <- df %>%
  count(mechanism, category) %>%
  complete(mechanism, category = categories, fill = list(n = 0)) %>%
  group_by(mechanism) %>%
  mutate(init = n / sum(n)) %>%
  ungroup() %>%
  select(init, categories = category, mechanism)

fig_s3 <- ggplot(df2, aes(fill=factor(mechanism, levels = mech_names_adj),
                     y=init, x=categories)) + 
  geom_bar(position="dodge", stat="identity", color="black") +
  xlab("% Amino acid identity") +
  ylab("Propotion") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_palette, name = "Resistance mechanism",
  labels = c("AAC (n = 325)",
             "APH (n = 451)",
             "Class A, C, D (n = 521)",
             "Class B (n = 155)",
             "Erm (n = 148)",
             "Mph (n = 52)",
             "Tet efflux (n = 103)",
             "Tet enzyme (n = 36)",
             "Tet RPG (n = 364)",
             "Qnr (n = 27)")) +
  fig_theme

pdf("sup_fig_3.pdf", width = 9, height = 6)
plot(fig_s3)
dev.off()

#----------------- Supplementary Fig 4 ----------------------------------------
# Load packages
library(caret)
library(dplyr)
library(forcats)
library(data.table)
library(randomForest)
library(pROC)
library(data.table)
library(patchwork)

# Functions
remove_redundant_events <- function(events) {
  uniq_events <- unique(events[, c("Node", "Gene.class")])
  subtable <- data.frame(matrix(nrow = 0, ncol = ncol(events)))
  colnames(subtable) <- colnames(events)
  
  for (i in seq_len(nrow(uniq_events))) {
    redundant_events <- events[events$Node == uniq_events$Node[i] & events$Gene.class == uniq_events$Gene.class[i], ]
    if (nrow(redundant_events) > 10) {
      subtable <- rbind(subtable, redundant_events[sample(seq_len(nrow(redundant_events)), 10), ])
    } else {
      subtable <- rbind(subtable, redundant_events)
    }
  }
  return(subtable)
}

subsample_null_events <- function(true, null) {
  gene_class <- unique(true$Gene.class)
  downsampled_events <- data.frame()
  tmp <- data.frame()
  
  for (class in gene_class) {
    n <- sum(true$Gene.class == class)
    subset <- null[null$Gene.class == class, ]
    if (nrow(subset) > n) {
      selection <- sample(seq_len(nrow(subset)), n)
      downsampled_events <- rbind(downsampled_events, subset[selection, ])
      tmp <- rbind(tmp, subset[-selection, ])
    } else {
      downsampled_events <- rbind(downsampled_events, subset)
    }
  }
  
  n_left <- nrow(true) - nrow(downsampled_events)
  downsampled_events <- rbind(downsampled_events, tmp[sample(seq_len(nrow(tmp)), n_left), ])
  
  return(downsampled_events)
}

format_input_data <- function(true_data, null_data, mechanism) {
  true_data$Transfer <- 1
  null_data$Transfer <- 0
  
  if (mechanism != "all") {
    if (mechanism == "class_A_C_D") {
      true_data <- true_data[true_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
      null_data <- null_data[null_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
    }
    
    else {
      true_data <- true_data[grep(mechanism, true_data$Gene.class), ]
      null_data <- null_data[grep(mechanism, null_data$Gene.class), ]
    }
  }
  
  true_data <- remove_redundant_events(true_data)
  if (nrow(null_data) > nrow(true_data)) {
    null_data <- subsample_null_events(true_data, null_data)
  }
  
  input_data <- rbind(true_data, null_data)
  input_data$Transfer <- as.factor(input_data$Transfer)
  input_data$Gene.class <- as.factor(input_data$Gene.class)
  
  input_data <- input_data %>%
    mutate(
      NN = as.integer(Gram_stain_difference == "NN"),
      PP = as.integer(Gram_stain_difference == "PP"),
      NP = as.integer(Gram_stain_difference == "NP")
    )
  
  return(input_data)
}

get_distance_top_left <- function(sens_spec) {
  d <- sqrt(sum((sens_spec - c(1,1))^2))
  return(d)
}

# Plotting settings
fig_theme <- theme(axis.text.x=element_text(size = rel(1.3)),
                   axis.text.y=element_text(size = rel(1.3)),
                   axis.title=element_text(size=rel(1.3)),
                   legend.title=element_text(size=rel(1.2)), 
                   legend.text=element_text(size=rel(1)),
                   plot.title = element_text(size = rel(2)),
                   plot.tag = element_text(size = rel(3), face="bold"),
                   plot.margin = ggplot2:::margin(t = 10,
                                                  r = 10,
                                                  b = 10,
                                                  l = 10))

# Analysis
observed_transfers <- data.frame(fread("data_for_figures/observed_horizontal_transfers.txt")) %>% 
  subset(select = -c(Header1, Header2)) %>%
  na.omit()

genes <- c("all", "aac", "aph", "class_A_C_D", "class_B", "erm", "tet_efflux", "tet_rpg")
titles <- c("All resistance mechanisms", "AAC", "APH", "Class A, C, D", "Class B", "Erm", "Tet efflux", "Tet RPG")
tags <- c("a", "b", "c", "d", "e", "f", "g", "h")
j <- 1

for (class in genes) {
  
  print(paste(c("Generating", class, "models"), collapse=" "))
  
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  for (i in 1:10){
    randomized_transfers <- read.delim(paste(c("data_for_figures/randomized_transfers", i, ".txt"), collapse = "")) %>%
      na.omit()
    
    input_data <- format_input_data(observed_transfers, randomized_transfers, class) %>% 
      select(Gene.class, NN, PP, NP, Genome_5mer_distance,
             Gene_genome_5mer_distance, Genome_size_difference,
             Animal, Human, Soil, Water, Wastewater, Transfer)
    
    set.seed(1)
    
    train_index <- createDataPartition(input_data$Transfer, times = 1, p = 0.7, list = FALSE)
    train_set <- input_data %>% slice(train_index)
    test_set <- input_data %>% slice(-train_index)
    
    rf_model <- randomForest(Transfer ~ .,
                             data = train_set, 
                             ntree = 500, 
                             type = "classification", 
                             na.action=na.omit, 
                             importance = TRUE)
    
    p_hat_test <- as.data.frame(predict(rf_model, newdata = test_set, type = "prob"))
    suppressMessages(assign(paste(c("roc", i), collapse = ""), roc(test_set$Transfer ~ p_hat_test$`1`, plot=FALSE, smooth = FALSE, quiet = TRUE)))
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Compile results
  sensitivities <- list(roc1$sensitivities, roc2$sensitivities,
                        roc3$sensitivities, roc4$sensitivities,
                        roc5$sensitivities, roc6$sensitivities,
                        roc7$sensitivities, roc8$sensitivities,
                        roc9$sensitivities, roc10$sensitivities)
  
  specificities <- list(roc1$specificities, roc2$specificities,
                        roc3$specificities, roc4$specificities,
                        roc5$specificities, roc6$specificities,
                        roc7$specificities, roc8$specificities,
                        roc9$specificities, roc10$specificities)
  
  roc_data <- data.frame(Sensitivity = tapply(unlist(sensitivities), sequence(lengths(sensitivities)), mean), 
                         Specificity = tapply(unlist(specificities), sequence(lengths(specificities)), mean))
  
  euc_dist <- apply(roc_data, 1, get_distance_top_left)
  
  # Plot and save subfigure
  p <- ggplot() +
    geom_line(aes(x = roc1$specificities, y = roc1$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc2$specificities, y = roc2$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc3$specificities, y = roc3$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc4$specificities, y = roc4$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc5$specificities, y = roc5$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc6$specificities, y = roc6$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc7$specificities, y = roc7$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc8$specificities, y = roc8$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc9$specificities, y = roc9$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc10$specificities, y = roc10$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
    geom_line(aes(x = roc_data$Specificity, y = roc_data$Sensitivity), col = "#000000", linewidth=1.1) +
    geom_point(aes(x = roc_data$Specificity[euc_dist==min(euc_dist)], y = roc_data$Sensitivity[euc_dist==min(euc_dist)]), size = 4.5, color= "#000000") +
    scale_x_reverse(labels = function(x) 1 - x) +
    theme_minimal() +
    labs(title = paste(c(titles[j], " (n = ", sum(input_data$Transfer == 1), ")"), collapse = ""), tag = tags[j]) +
    xlab("1-Specificity") +
    ylab("Sensitivity") +
    fig_theme
  
  pdf(paste(c("sup_fig_4", tags[j], ".pdf"), collapse=""),  width=7, height=6)
  plot(p)
  dev.off()
  
  j <- j+1
  
}

#----------------- Supplementary Fig 5 ----------------------------------------
# Load packages
library(caret)
library(dplyr)
library(forcats)
library(data.table)
library(pROC)
library(boot)
library(rfPermute)
library(rfUtilities)
library(parallel)
library(patchwork)

# Functions
remove_redundant_events <- function(events) {
  uniq_events <- unique(events[, c("Node", "Gene.class")])
  subtable <- data.frame(matrix(nrow = 0, ncol = ncol(events)))
  colnames(subtable) <- colnames(events)
  
  for (i in seq_len(nrow(uniq_events))) {
    redundant_events <- events[events$Node == uniq_events$Node[i] & events$Gene.class == uniq_events$Gene.class[i], ]
    if (nrow(redundant_events) > 10) {
      subtable <- rbind(subtable, redundant_events[sample(seq_len(nrow(redundant_events)), 10), ])
    } else {
      subtable <- rbind(subtable, redundant_events)
    }
  }
  return(subtable)
}

subsample_null_events <- function(true, null) {
  gene_class <- unique(true$Gene.class)
  downsampled_events <- data.frame()
  tmp <- data.frame()
  
  for (class in gene_class) {
    n <- sum(true$Gene.class == class)
    subset <- null[null$Gene.class == class, ]
    if (nrow(subset) > n) {
      selection <- sample(seq_len(nrow(subset)), n)
      downsampled_events <- rbind(downsampled_events, subset[selection, ])
      tmp <- rbind(tmp, subset[-selection, ])
    } else {
      downsampled_events <- rbind(downsampled_events, subset)
    }
  }
  
  n_left <- nrow(true) - nrow(downsampled_events)
  downsampled_events <- rbind(downsampled_events, tmp[sample(seq_len(nrow(tmp)), n_left), ])
  
  return(downsampled_events)
}

format_input_data <- function(true_data, null_data, mechanism) {
  true_data$Transfer <- 1
  null_data$Transfer <- 0
  
  if (mechanism != "all") {
    if (mechanism == "class_A_C_D") {
      true_data <- true_data[true_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
      null_data <- null_data[null_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
    }
    
    else {
      true_data <- true_data[grep(mechanism, true_data$Gene.class), ]
      null_data <- null_data[grep(mechanism, null_data$Gene.class), ]
    }
  }
  
  true_data <- remove_redundant_events(true_data)
  if (nrow(null_data) > nrow(true_data)) {
    null_data <- subsample_null_events(true_data, null_data)
  }
  
  input_data <- rbind(true_data, null_data)
  input_data$Transfer <- as.factor(input_data$Transfer)
  input_data$Gene.class <- as.factor(input_data$Gene.class)
  
  input_data <- input_data %>%
    mutate(
      NN = as.integer(Gram_stain_difference == "NN"),
      PP = as.integer(Gram_stain_difference == "PP"),
      NP = as.integer(Gram_stain_difference == "NP")
    )
  
  return(input_data)
}

# Plotting settings
fig_theme <- theme(axis.text.x=element_text(size = rel(2.4), angle = 45, vjust = 1, hjust = 1),
                   axis.text.y=element_text(size = rel(2.4)),
                   axis.title=element_text(size=rel(2.4)),
                   legend.title=element_text(size=rel(2.25)), 
                   legend.text=element_text(size=rel(1.95)),
                   plot.title = element_text(size = rel(3)),
                   plot.tag = element_text(size = rel(4.5), face="bold"),
                   plot.margin = ggplot2:::margin(t = 10,
                                                  r = 10,
                                                  b = 10,
                                                  l = 40))

feature_category_palette <- c(
  "Genetic incompatibility" ="#225ea8",
  "Environmental co-occurrence"='#41b6c4',
  "Gram staining" ="#a1dab4",
  "Other"="#ffffcc"
)

# Analysis
observed_transfers <- data.frame(fread("data_for_figures/observed_horizontal_transfers.txt")) %>% 
  subset(select = -c(Header1, Header2)) %>%
  na.omit()

genes <- c("aac", "aph", "class_A_C_D", "class_B", "erm", "tet_efflux", "tet_rpg")
titles <- c("AAC", "APH", "Class A, C, D", "Class B", "Erm", "Tet efflux", "Tet RPG")
tags <- c("a", "b", "c", "d", "e", "f", "g")
k <- 1

for (class in genes) {
  imp_data <- data.frame()
  
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  for (j in 1:10) {
    randomized_transfers <- read.delim(paste(c("data_for_figures/randomized_transfers", j, ".txt"), collapse = "")) %>%
      na.omit()
    
    input_data <- format_input_data(observed_transfers, randomized_transfers, class) %>% 
      select(Gene.class, NN, PP, NP, Genome_5mer_distance,
             Gene_genome_5mer_distance, Genome_size_difference,
             Animal, Human, Soil, Water, Wastewater, Transfer)
    
    set.seed(1)
    
    train_index <- createDataPartition(input_data$Transfer, times = 1, p = 0.7, list = FALSE)
    train_set <- input_data %>% slice(train_index)
    test_set <- input_data %>% slice(-train_index)
    
    rf_permute <- rfPermute(Transfer ~ .,
                            data = train_set,
                            na.action = na.omit,
                            num.cores = detectCores()-1)
    
    imp <- importance(rf_permute, scale = FALSE)[,5]
    pval <-importance(rf_permute, scale = FALSE)[,6]
    
    imp <- data.frame(MeanDecreaseAccuracy = imp, p.value = pval)
    imp$Feature <- rownames(imp)
    rownames(imp) <- NULL
    
    imp_data <- rbind(imp_data, imp)
    setTxtProgressBar(pb, j)
  }
  
  close(pb)
  
  # Compile results
  plot_data <- data.frame(Feature = c("NN", "PP", "NP", "Genome_5mer_distance", 
                                    "Gene_genome_5mer_distance", "Genome_size_difference", 
                                    "Animal", "Human", "Soil", "Water", "Wastewater", 
                                    "Gene.class"), 
                        mean_imp = rep(0, times = 12),
                        sd_imp = rep(0, times = 12), 
                        significance = rep("", times = 12))
  
  for (f in plot_data$Feature){
    plot_data$mean_imp[plot_data$Feature == f] <- mean(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f])
    plot_data$sd_imp[plot_data$Feature == f] <- sd(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f])
    
    if(sum(imp_data$p.value[imp_data$Feature == f] < 0.01) == 10) {
      plot_data$significance[plot_data$Feature == f] <- "*"
    }
  }
  
  plot_data$var_categ <- rep(NA, times = nrow(plot_data))
  
  plot_data$Feature[plot_data$Feature == "Genome_5mer_distance"] <- "Genome 5mer distance"
  plot_data$Feature[plot_data$Feature == "Gene_genome_5mer_distance"] <- "Genome-gene 5mer distance"
  plot_data$Feature[plot_data$Feature == "Genome_size_difference"] <- "Difference in genome size"
  plot_data$Feature[plot_data$Feature == "Gene.class"] <- "Gene class"
  plot_data$Feature[plot_data$Feature == "NN"] <- "G-/G-"
  plot_data$Feature[plot_data$Feature == "PP"] <- "G+/G+"
  plot_data$Feature[plot_data$Feature == "NP"] <- "G-/G+"
  
  genomic <- c("Genome 5mer distance",
               "Genome-gene 5mer distance",
               "Difference in genome size")
  
  cooccurrence <- c("Animal",
                    "Human",
                    "Soil",
                    "Water",
                    "Wastewater")
  
  gram_diff <- c("G-/G-", "G+/G+", "G-/G+")
  
  for (i in 1:nrow(plot_data)){
    if (plot_data$Feature[i] %in% genomic){
      plot_data$var_categ[i] <- "Genetic incompatibility"
    }
    else if (plot_data$Feature[i] %in% cooccurrence){
      plot_data$var_categ[i] <- "Environmental co-occurrence"
    } 
    else if (plot_data$Feature[i] %in% gram_diff){
      plot_data$var_categ[i] <- "Gram staining"
    }
    else {
      plot_data$var_categ[i] <- "Other"
    }
  }
  
  # Plot subfigure
  p <- ggplot(plot_data, aes(x=factor(Feature, levels=Feature[order(mean_imp, decreasing = TRUE)]), y=mean_imp, fill=as.factor(var_categ))) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    geom_errorbar(aes(ymin=mean_imp-sd_imp, ymax=mean_imp+sd_imp), width=.3, position=position_dodge(.9)) +
    xlab("") +
    ylab("MeanDecreaseAccuracy") +
    labs(title = titles[k], tag = tags[k]) +
    scale_fill_manual(values = feature_category_palette, 
    name = "Feature group",
    breaks=c("Genetic incompatibility",
             "Environmental co-occurrence",
             "Gram staining",
             "Other")) +
    theme_minimal() +
    geom_text(aes(label=significance),position=position_dodge(width=0.75), hjust=0.6, vjust = 0.3, size=10) +
    fig_theme
  
  assign(paste(c("p", k), collapse = ""), p)
  
  k <- k+1
  
}

# Save Sup. Fig 5
pdf("sup_fig_5.pdf",  width=36, height=30)
(p1 | p2 | p3) /
  (p4 | p5 | p6) /
  (p7 | plot_spacer()  | plot_spacer() | plot_spacer())
dev.off()

#----------------- Supplementary Fig 6, 7 -------------------------------------
# Load packages
library(caret)
library(dplyr)
library(forcats)
library(data.table)
library(randomForest)
library(pROC)
library(boot)
library(rfPermute)
library(ggpattern)
library(rfUtilities)
library(pdp)
library(aplot)

# Functions
remove_redundant_events <- function(events) {
  uniq_events <- unique(events[, c("Node", "Gene.class")])
  subtable <- data.frame(matrix(nrow = 0, ncol = ncol(events)))
  colnames(subtable) <- colnames(events)
  
  for (i in seq_len(nrow(uniq_events))) {
    redundant_events <- events[events$Node == uniq_events$Node[i] & events$Gene.class == uniq_events$Gene.class[i], ]
    if (nrow(redundant_events) > 10) {
      subtable <- rbind(subtable, redundant_events[sample(seq_len(nrow(redundant_events)), 10), ])
    } else {
      subtable <- rbind(subtable, redundant_events)
    }
  }
  return(subtable)
}

subsample_null_events <- function(true, null) {
  gene_class <- unique(true$Gene.class)
  downsampled_events <- data.frame()
  tmp <- data.frame()
  
  for (class in gene_class) {
    n <- sum(true$Gene.class == class)
    subset <- null[null$Gene.class == class, ]
    if (nrow(subset) > n) {
      selection <- sample(seq_len(nrow(subset)), n)
      downsampled_events <- rbind(downsampled_events, subset[selection, ])
      tmp <- rbind(tmp, subset[-selection, ])
    } else {
      downsampled_events <- rbind(downsampled_events, subset)
    }
  }
  
  n_left <- nrow(true) - nrow(downsampled_events)
  downsampled_events <- rbind(downsampled_events, tmp[sample(seq_len(nrow(tmp)), n_left), ])
  
  return(downsampled_events)
}

format_input_data <- function(true_data, null_data, mechanism) {
  true_data$Transfer <- 1
  null_data$Transfer <- 0
  
  if (mechanism != "all") {
    if (mechanism == "class_A_C_D") {
      true_data <- true_data[true_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
      null_data <- null_data[null_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
    }
    
    else {
      true_data <- true_data[grep(mechanism, true_data$Gene.class), ]
      null_data <- null_data[grep(mechanism, null_data$Gene.class), ]
    }
  }
  
  true_data <- remove_redundant_events(true_data)
  if (nrow(null_data) > nrow(true_data)) {
    null_data <- subsample_null_events(true_data, null_data)
  }
  
  input_data <- rbind(true_data, null_data)
  input_data$Transfer <- as.factor(input_data$Transfer)
  input_data$Gene.class <- as.factor(input_data$Gene.class)
  
  input_data <- input_data %>%
    mutate(
      NN = as.integer(Gram_stain_difference == "NN"),
      PP = as.integer(Gram_stain_difference == "PP"),
      NP = as.integer(Gram_stain_difference == "NP")
    )
  
  return(input_data)
}

compile_distance_data <- function(input_data, variable) {
  distance_data <- input_data[, c("Taxonomic.difference", variable, "Transfer")]
  distance_data <- distance_data %>%
    mutate(
      Taxonomic.difference = recode(
        Taxonomic.difference,
        order = "Order",
        class = "Class",
        phylum = "Phylum"
      ) %>%
        as.factor() %>%
        relevel("Order")
    ) %>%
    mutate(
      Transfer = recode(
        Transfer,
        "1" = "Observed",
        "0" = "Random"
      ) %>%
        as.factor() %>%
        relevel("Random")
    )
}

create_partial_plot <- function(model, input_data, pred_var, class_idx, xlab, ylab, title, tag) {
  part <- partial(model, pred.var = pred_var, type = "classification", which.class = class_idx)
  df <- data.frame(x = part[[pred_var]], y = part$yhat)
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    geom_smooth(size = 1.3, color = "#4575b4", method = "gam") +
    theme_minimal() +
    labs(title = title,
         tag = tag,
         x = xlab,
         y = ylab) +
    ylim(-2, 2) +
    theme(
      axis.text.x = element_text(size = rel(1.3)),
      axis.text.y = element_text(size = rel(1.3)),
      axis.title = element_text(size = rel(1.3)),
      plot.title = element_text(size = rel(2)),
      plot.tag = element_text(size = rel(3), face = "bold")
    )
  
  return(p)
}

combine_plots <- function(main_plot, subplots, heights) {
  combined <- main_plot
  for (i in seq_along(subplots)) {
    combined <- insert_bottom(combined, subplots[[i]], height = heights[i])
  }
  return(combined)
}

# Analysis
observed_transfers <- data.frame(fread("data_for_figures/observed_horizontal_transfers.txt")) %>% 
  subset(select = -c(Header1, Header2)) %>%
  na.omit()

randomized_transfers <- read.delim(paste(c("data_for_figures/randomized_transfers1.txt"), collapse = "")) %>% 
  na.omit()

genes <- c("aac", "aph", "class_A_C_D", "class_B", "erm", "tet_efflux", "tet_rpg")
titles <- c("AAC", "APH", "Class A, C, D", "Class B", "Erm", "Tet efflux", "Tet RPG")
tags <- c("a", "b", "c", "d", "e", "f", "g")

for (i in 1:length(genes)) {
  input_data <- format_input_data(observed_transfers, randomized_transfers, genes[i]) %>% 
    subset(select = c(Taxonomic.difference, Gene.class, NN, PP, NP, Genome_5mer_distance,
                      Gene_genome_5mer_distance, Genome_size_difference,
                      Animal, Human, Soil, Water, Wastewater, Transfer))
  
  set.seed(1)

  train_index <- createDataPartition(input_data$Transfer, times = 1, p = 0.7, list = FALSE)
  train_set <- input_data %>% slice(train_index)
  test_set <- input_data %>% slice(-train_index)
  
  rf_model <- randomForest(Transfer ~ 
                             Gene.class +
                             NN +
                             PP +
                             NP +
                             Genome_5mer_distance + 
                             Gene_genome_5mer_distance +
                             Genome_size_difference + 
                             Animal +
                             Human +
                             Soil +
                             Water +
                             Wastewater,
                           data = train_set, 
                           ntree = 500, 
                           type = "classification", 
                           na.action=na.omit, 
                           importance = TRUE)
  
  # Compile results for Sup. Fig 6
  distance_data <- compile_distance_data(input_data, "Genome_5mer_distance") %>%
    mutate(Taxonomic.difference = factor(
      Taxonomic.difference,
      levels = c("Phylum", "Class", "Order") 
    ))
  
  # Plot and save subfigure
  p_main <- create_partial_plot(
    rf_model, input_data, "Genome_5mer_distance", 2,
    xlab = "Genome 5mer distance",
    ylab = "Partial dependence",
    title = titles[i],
    tag = tags[i]
  )

  box_plots <- lapply(levels(distance_data[["Taxonomic.difference"]]), function(level) {
    ggplot(distance_data %>% filter(!!sym("Taxonomic.difference") == level), 
           aes(x = Genome_5mer_distance, y = Transfer, fill = Transfer)) +
      geom_boxplot() +
      theme_minimal() +
      labs(x = NULL, y = level) +
      scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
      theme(
        axis.text.y = element_blank(),
        legend.position = ifelse(level == "Phylum", "right", "none"),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.3)),
        legend.key.size = unit(0.75, 'cm')
      )
  })
  
  combined_plot <- combine_plots(p_main, box_plots, heights = c(0.2, 0.2, 0.2))
  
  pdf(paste(c("sup_fig_6_", tags[i], ".pdf"), collapse=""), width = 7, height = 5)
  print(combined_plot)
  dev.off()
  
  # Compile results for Sup. Fig 7
  distance_data <- compile_distance_data(input_data, "Gene_genome_5mer_distance") %>%
    mutate(Taxonomic.difference = factor(
      Taxonomic.difference,
      levels = c("Phylum", "Class", "Order") 
    ))
  
  # Plot and save subfigure
  p_main <- create_partial_plot(
    rf_model, input_data, "Gene_genome_5mer_distance", 2,
    xlab = "Gene-genome 5mer distance",
    ylab = "Partial dependence",
    title = titles[i],
    tag = tags[i]
  )

  box_plots <- lapply(levels(distance_data[["Taxonomic.difference"]]), function(level) {
    ggplot(distance_data %>% filter(!!sym("Taxonomic.difference") == level), 
           aes(x = Gene_genome_5mer_distance, y = Transfer, fill = Transfer)) +
      geom_boxplot() +
      theme_minimal() +
      labs(x = NULL, y = level) +
      scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
      theme(
        axis.text.y = element_blank(),
        legend.position = ifelse(level == "Phylum", "right", "none"),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.3)),
        legend.key.size = unit(0.75, 'cm')
      )
  })
  
  combined_plot <- combine_plots(p_main, box_plots, heights = c(0.2, 0.2, 0.2))
  
  pdf(paste(c("sup_fig_7_", tags[i], ".pdf"), collapse=""), width = 7, height = 5)
  print(combined_plot)
  dev.off()
}

#----------------- Supplementary Fig 8 ----------------------------------------
# Load packages
library(igraph)
library(stringr)
library(GGally)
library(ggplot2)
library(tidyr)
library(data.table)
library(dplyr)

# Functions
plot_cooccurrence_network <- function(graph, environment) {
  size <- c()
  for (name in V(graph)$name) {
    size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
  }
  
  genus <- c()
  for (name in V(graph)$name) {
    genus <- c(genus, str_split(name, " ")[[1]][1])
  }
  
  V(graph)$label <- taxonomy[genus,3]
  
  for (i in 1:length(V(graph)$label)) {
    if (is.na(V(graph)$label[i])) {
      V(graph)$label[i] <- taxonomy[taxonomy$family == genus[i],3][1]
    }
  }
  
  env_list <- c("Animal", "Human", "Soil", "Water", "Wastewater")
  
  imp <- edge_attr(graph)[[grep(environment, env_list)]]
  
  imp[is.na(imp)] <- 0
  imp[imp < 0.01] <- 0
  texture <- c()
  alph <- c()
  
  for (i in 1:length(imp)) {
    if (is.na(imp[i])) {
      texture <- c(texture, 0)
      alph <- c(alph, 0)
    }
    else if (imp[i] == 0) {
      texture <- c(texture, 3)
      alph <- c(alph, 0)
    }
    else if (imp[i] > 0 && imp[i] <= 0.01) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.01 && imp[i] <= 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
  }
  
  E(graph)$lty <- texture
  E(graph)$color <- "#bdbdbd"
  
  node_color <- c("#41ab5d", "#005a32", "#9ecae1", "#4292c6", "#084594", 
                  "#feb24c", "#dd3497", "#fc9272", "#ef3b2c", "#99000d")
  
  names(node_color) <-c("Actinomycetes", "Coriobacteriia", "Bacilli", "Clostridia",
                        "Erysipelotrichia", "Bacteroidia", "Epsilonproteobacteria",
                        "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")
  
  esize <- sqrt(imp*20)
  esize[esize == 0] <- 0.5
  
  shape <- rep(15, length(V(graph)$name))
  shape[grep(" ", V(graph)$name)] <- 19
  
  set.seed(456)
  
  p <- ggnet2(graph, 
              mode = "kamadakawai", 
              color=V(graph)$label, 
              palette = node_color, 
              edge.alpha = alph, 
              edge.size = 2*esize, 
              alpha = 1, size = 10*size, 
              edge.color = E(graph)$color, 
              edge.lty = E(graph)$lty, 
              legend.size = 9, 
              layout.exp = 0.1,
              node.shape = shape) +
    geom_text(aes(label = V(graph)$name), color = "black", size=1.5, fontface = "bold") +
    guides(size = FALSE) +
    coord_fixed()
}

# Load data
taxonomy <- data.frame(fread("data_for_figures/assembly_full_lineage.txt", header=FALSE)) %>%
  select(-c(1,ncol(.))) %>%
  na.omit() %>%
  unique() %>%
  `colnames<-`(c("kingdom", "phylum", "class", "order", "family", "genus")) %>%
  `rownames<-`(.$genus)

co_occurrence <- read.table("data_for_figures/species_cooccurrence.txt", header = TRUE, sep="\t") %>%
  .[!(is.na(.$Species2) | is.na(.$Species1)),]

# Filter non-frequent pairs
df <- co_occurrence[co_occurrence$n.true >= 5, ]

# Generate network
g <- graph.data.frame(d = df, directed = FALSE)

# Plot and save subfigures
fig_s8a <- plot_cooccurrence_network(g, "Animal")

pdf("sup_fig_8a.pdf", width = 70, height = 35)
plot(fig_s8a)
dev.off()

fig_s8b <- plot_cooccurrence_network(g, "Human")

pdf("sup_fig_8b.pdf", width = 70, height = 35)
plot(fig_s8b)
dev.off()

fig_s8c <- plot_cooccurrence_network(g, "Soil")

pdf("sup_fig_8c.pdf", width = 70, height = 35)
plot(fig_s8c)
dev.off()

fig_s8d <- plot_cooccurrence_network(g, "Water")

pdf("sup_fig_8d.pdf", width = 70, height = 35)
plot(fig_s8d)
dev.off()

fig_s8e <- plot_cooccurrence_network(g, "Wastewater")

pdf("sup_fig_8e.pdf", width = 70, height = 35)
plot(fig_s8e)
dev.off()

#----------------- Supplementary Fig 9 ----------------------------------------
# Load packages
library(caret)
library(dplyr)
library(forcats)
library(data.table)
library(randomForest)
library(pROC)
library(boot)
library(rfPermute)
library(rfUtilities)
library(parallel)
library(patchwork)


# Functions
remove_redundant_events <- function(events) {
  uniq_events <- unique(events[, c("Node", "Gene.class")])
  subtable <- data.frame(matrix(nrow = 0, ncol = ncol(events)))
  colnames(subtable) <- colnames(events)
  
  for (i in seq_len(nrow(uniq_events))) {
    redundant_events <- events[events$Node == uniq_events$Node[i] & events$Gene.class == uniq_events$Gene.class[i], ]
    if (nrow(redundant_events) > 10) {
      subtable <- rbind(subtable, redundant_events[sample(seq_len(nrow(redundant_events)), 10), ])
    } else {
      subtable <- rbind(subtable, redundant_events)
    }
  }
  return(subtable)
}

subsample_null_events <- function(true, null) {
  gene_class <- unique(true$Gene.class)
  downsampled_events <- data.frame()
  tmp <- data.frame()
  
  for (class in gene_class) {
    n <- sum(true$Gene.class == class)
    subset <- null[null$Gene.class == class, ]
    if (nrow(subset) > n) {
      selection <- sample(seq_len(nrow(subset)), n)
      downsampled_events <- rbind(downsampled_events, subset[selection, ])
      tmp <- rbind(tmp, subset[-selection, ])
    } else {
      downsampled_events <- rbind(downsampled_events, subset)
    }
  }
  
  n_left <- nrow(true) - nrow(downsampled_events)
  downsampled_events <- rbind(downsampled_events, tmp[sample(seq_len(nrow(tmp)), n_left), ])
  
  return(downsampled_events)
}

format_input_data <- function(true_data, null_data, mechanism) {
  true_data$Transfer <- 1
  null_data$Transfer <- 0
  
  if (mechanism != "all") {
    if (mechanism == "class_A_C_D") {
      true_data <- true_data[true_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
      null_data <- null_data[null_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
    }
    
    else {
      true_data <- true_data[grep(mechanism, true_data$Gene.class), ]
      null_data <- null_data[grep(mechanism, null_data$Gene.class), ]
    }
  }
  
  true_data <- remove_redundant_events(true_data)
  if (nrow(null_data) > nrow(true_data)) {
    null_data <- subsample_null_events(true_data, null_data)
  }
  
  input_data <- rbind(true_data, null_data)
  input_data$Transfer <- as.factor(input_data$Transfer)
  input_data$Gene.class <- as.factor(input_data$Gene.class)
  
  input_data <- input_data %>%
    mutate(
      NN = as.integer(Gram_stain_difference == "NN"),
      PP = as.integer(Gram_stain_difference == "PP"),
      NP = as.integer(Gram_stain_difference == "NP")
    )
  
  return(input_data)
}

get_distance_top_left <- function(sens_spec) {
  d <- sqrt(sum((sens_spec - c(1,1))^2))
  return(d)
}

calc_effect_size <- function(model) {
  es <- data.frame(matrix(nrow=1, ncol=15))
  colnames(es) <- c("NN", "PP", "NP", "Genome_5mer_distance", "Gene_genome_5mer_distance", "Genome_size_difference", "Animal", "Human.gut", "Human.skin", "Human.oral", "Freshwater", "Marine", "Soil", "Wastewater", "Gene.class")
  
  nn <- try(rf.effectSize(model$rf, y = test_set$NN, pred.data = test_set, x.var = NN), silent = TRUE)
  if (inherits(nn, "try-error") | is.na(nn)) {
    nn <- rf.effectSize(model$rf, y = train_set$NN, pred.data = train_set, x.var = NN)
  }
  es[1, "NN"] <- nn
  
  pp <- try(rf.effectSize(model$rf, y = test_set$PP, pred.data = test_set, x.var = PP), silent = TRUE)
  if (inherits(pp, "try-error") | is.na(pp)) {
    pp <- rf.effectSize(model$rf, y = train_set$PP, pred.data = train_set, x.var = PP)
  }
  es[1, "PP"] <- pp
  
  np <- try(rf.effectSize(model$rf, y = test_set$NP, pred.data = test_set, x.var = NP), silent = TRUE)
  if (inherits(np, "try-error") | is.na(np)) {
    np <- rf.effectSize(model$rf, y = train_set$NP, pred.data = train_set, x.var = NP)
  }
  es[1, "NP"] <- np
  
  
  gdist <- try(rf.effectSize(model$rf, y = test_set$Genome_5mer_distance, pred.data = test_set, x.var = Genome_5mer_distance), silent = TRUE)
  if (inherits(gdist, "try-error") | is.na(gdist)) {
    gdist <- rf.effectSize(model$rf, y = train_set$Genome_5mer_distance, pred.data = train_set, x.var = Genome_5mer_distance)
  }
  es[1, "Genome_5mer_distance"] <- gdist
  
  ggdist <- try(rf.effectSize(model$rf, y = test_set$Gene_genome_5mer_distance, pred.data = test_set, x.var = Gene_genome_5mer_distance), silent = TRUE)
  if (inherits(ggdist, "try-error") | is.na(ggdist)) {
    ggdist <- rf.effectSize(model$rf, y = train_set$Gene_genome_5mer_distance, pred.data = train_set, x.var = Gene_genome_5mer_distance)
  }
  es[1, "Gene_genome_5mer_distance"] <- ggdist
  
  gsize <- try(rf.effectSize(model$rf, y = test_set$Genome_size_difference, pred.data = test_set, x.var = Genome_size_difference), silent = TRUE)
  if (inherits(gsize, "try-error") | is.na(gsize)) {
    gsize <- rf.effectSize(model$rf, y = train_set$Genome_size_difference, pred.data = train_set, x.var = Genome_size_difference)
  }
  es[1, "Genome_size_difference"] <- gsize
  
  animal <- try(rf.effectSize(model$rf, y = test_set$Animal, pred.data = test_set, x.var = Animal), silent = TRUE)
  if (inherits(animal, "try-error") | is.na(animal)) {
    animal <- rf.effectSize(model$rf, y = train_set$Animal, pred.data = train_set, x.var = Animal)
  }
  es[1, "Animal"] <- animal
  
  hgut <- try(rf.effectSize(model$rf, y = test_set$Human.gut, pred.data = test_set, x.var = Human.gut), silent = TRUE)
  if (inherits(hgut, "try-error") | is.na(hgut)) {
    hgut <- rf.effectSize(model$rf, y = train_set$Human.gut, pred.data = train_set, x.var = Human.gut)
  }
  es[1, "Human.gut"] <- hgut
  
  hskin <- try(rf.effectSize(model$rf, y = test_set$Human.skin, pred.data = test_set, x.var = Human.skin), silent = TRUE)
  if (inherits(hskin, "try-error") | is.na(hskin)) {
    hskin <- rf.effectSize(model$rf, y = train_set$Human.skin, pred.data = train_set, x.var = Human.skin)
  }
  es[1, "Human.skin"] <- hskin
  
  horal <- try(rf.effectSize(model$rf, y = test_set$Human.oral, pred.data = test_set, x.var = Human.oral), silent = TRUE)
  if (inherits(horal, "try-error") | is.na(horal)) {
    horal <- rf.effectSize(model$rf, y = train_set$Human.oral, pred.data = train_set, x.var = Human.oral)
  }
  es[1, "Human.oral"] <- horal
  
  fwater <- try(rf.effectSize(model$rf, y = test_set$Freshwater, pred.data = test_set, x.var = Freshwater), silent = TRUE)
  if (inherits(fwater, "try-error") | is.na(fwater)) {
    fwater <- rf.effectSize(model$rf, y = train_set$Freshwater, pred.data = train_set, x.var = Freshwater)
  }
  es[1, "Freshwater"] <- fwater
  
  marine <- try(rf.effectSize(model$rf, y = test_set$Marine, pred.data = test_set, x.var = Marine), silent = TRUE)
  if (inherits(marine, "try-error") | is.na(marine)) {
    marine <- rf.effectSize(model$rf, y = train_set$Marine, pred.data = train_set, x.var = Marine)
  }
  es[1, "Marine"] <- marine
  
  soil <- try(rf.effectSize(model$rf, y = test_set$Soil, pred.data = test_set, x.var = Soil), silent = TRUE)
  if (inherits(soil, "try-error") | is.na(soil)) {
    soil <- rf.effectSize(model$rf, y = train_set$Soil, pred.data = train_set, x.var = Soil)
  }
  es[1, "Soil"] <- soil
  
  wwater <- try(rf.effectSize(model$rf, y = test_set$Wastewater, pred.data = test_set, x.var = Wastewater), silent = TRUE)
  if (inherits(wwater, "try-error") | is.na(wwater)) {
    wwater <- rf.effectSize(model$rf, y = train_set$Wastewater, pred.data = train_set, x.var = Wastewater)
  }
  es[1, "Wastewater"] <- wwater
  
  es[1, "Gene.class"] <- 1
  
  return(es)
}

# Plotting settings
fig_a_theme <- theme(axis.text.x=element_text(size = rel(1.95), angle = 45, vjust = 1, hjust = 1),
                     axis.text.y=element_text(size = rel(1.95)),
                     axis.title=element_text(size=rel(1.95)),
                     legend.title=element_text(size=rel(1.8)), 
                     legend.text=element_text(size=rel(1.5)),
                     plot.title = element_text(size = rel(2.5)),
                     plot.tag = element_text(size = rel(5), face = "bold"))

fig_bcd_theme <- theme(axis.text.x=element_text(size = rel(1.95)),
                       axis.text.y=element_text(size = rel(1.95)),
                       axis.title=element_text(size=rel(1.95)),
                       legend.title=element_text(size=rel(1.8)), 
                       legend.text=element_text(size=rel(1.5)),
                       plot.title = element_text(size = rel(2.5)),
                       plot.tag = element_text(size = rel(5), face = "bold"))

feature_category_palette <- c(
  "Genetic incompatibility" ="#225ea8",
  "Environmental co-occurrence"='#41b6c4',
  "Gram staining" ="#a1dab4",
  "Other"="#ffffcc"
)

mechanism_palette <- c(
  "All" ="#000000",
  "AAC"='#a50026',
  "APH" ="#d73027",
  "Class A, C, D"='#fdae61',
  "Class B"="#fee090",
  "Erm" = "#ffffff",
  "Tet efflux" = "#74add1",
  "Tet RPG" = "#313695"
)

# Analysis
observed_transfers <- data.frame(fread("data_for_figures/observed_transfers_cooc_divided.txt")) %>% 
  subset(select = -c(Header1, Header2)) %>%
  na.omit()

df_total <- data.frame()
mechanisms <- c("all", "aac", "aph", "class_A_C_D", "class_B", "erm", "tet_efflux", "tet_rpg")

for (class in mechanisms){
  print(class)
  imp_data <- data.frame()
  df <- data.frame(Feature = c("NN", "PP", "NP", "Genome_5mer_distance", 
                               "Gene_genome_5mer_distance", "Genome_size_difference", 
                               "Animal", "Human.gut", "Human.skin", "Human.oral", 
                               "Freshwater", "Marine", "Soil", "Wastewater", 
                               "Gene.class"), 
                   mean_imp = rep(0, times = 15),
                   sd_imp = rep(0, times = 15), 
                   significance = rep("", times = 15))
  
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  for (j in 1:10) {
    randomized_transfers <- read.delim(paste(c("data_for_figures/randomized_transfers", j, "_cooc_divided", ".txt"), collapse = "")) %>% 
      na.omit()
    
    input_data <- format_input_data(observed_transfers, randomized_transfers, class) %>% 
      subset(select = c(Gene.class, NN, PP, NP, Genome_5mer_distance,
                        Gene_genome_5mer_distance, Genome_size_difference,
                        Animal, Human.gut, Human.skin, Human.oral, Freshwater, 
                        Marine, Soil, Wastewater, Transfer))
    
    set.seed(1)
    
    train_index <- createDataPartition(input_data$Transfer, times = 1, p = 0.7, list = FALSE)
    train_set <- input_data %>% slice(train_index)
    test_set <- input_data %>% slice(-train_index)
    
    rf_permute <- rfPermute(Transfer ~ ., 
                            data = train_set, 
                            na.action = na.omit,
                            num.cores = detectCores()-1)
    
    imp <- importance(rf_permute, scale = FALSE)[,5]
    pval <- importance(rf_permute, scale = FALSE)[,6]
    effect_size <- calc_effect_size(rf_permute)
    
    imp <- data.frame(MeanDecreaseAccuracy = imp, p.value = pval)
    imp$Feature <- rownames(imp)
    rownames(imp) <- NULL
    imp$MeanDecreaseAccuracy <- as.numeric(sign(effect_size[imp$Feature])*imp$MeanDecreaseAccuracy)
    
    imp_data <- rbind(imp_data, imp)
    
    setTxtProgressBar(pb, j)
  }
  
  for (f in df$Feature){
    df$mean_imp[df$Feature == f] <- mean(na.omit(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f]))
    df$sd_imp[df$Feature == f] <- sd(na.omit(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f]))
    
    if(sum(imp_data$p.value[imp_data$Feature == f] < 0.01) == 10) {
      df$significance[df$Feature == f] <- "*"
    }
  }
  
  df$class <- rep(class, times = 15)
  df_total <- rbind(df_total, df)
  
  if (class == "all") {
    # Compile results for Fig 3a
    plot_data_a <- data.frame(Feature = c("NN", "PP", "NP", "Genome_5mer_distance", 
                                          "Gene_genome_5mer_distance", "Genome_size_difference", 
                                          "Animal", "Human.gut", "Human.skin", "Human.oral", 
                                          "Freshwater", "Marine", "Soil", "Wastewater", 
                                          "Gene.class"), 
                              mean_imp = rep(0, times = 15),
                              sd_imp = rep(0, times = 15), 
                              significance = rep("", times = 15))
    
    for (f in plot_data_a$Feature){
      plot_data_a$mean_imp[plot_data_a$Feature == f] <- abs(mean(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f]))
      plot_data_a$sd_imp[plot_data_a$Feature == f] <- sd(imp_data$MeanDecreaseAccuracy[imp_data$Feature == f])
      
      if(sum(imp_data$p.value[imp_data$Feature == f] < 0.01) == 10) {
        plot_data_a$significance[plot_data_a$Feature == f] <- "*"
      }
    }
    
    plot_data_a$var_categ <- rep(NA, times = nrow(plot_data_a))
    
    feat_names_to_adjust <- c("Genome_5mer_distance", "Gene_genome_5mer_distance", 
                              "Genome_size_difference", "Human.gut", "Human.skin", 
                              "Human.oral", "Gene.class", "NN", "PP", "NP")
    feat_names_adj <- c("Genome 5mer distance", "Genome-gene 5mer distance", 
                        "Difference in genome size", "Human gut", "Human skin", 
                        "Human oral", "Gene class", "G-/G-","G+/G+", "G-/G+")
    
    for (i in 1:length(feat_names_to_adjust)) {
      plot_data_a$Feature[plot_data_a$Feature == feat_names_to_adjust[i]] <- feat_names_adj[i]
    }
    
    genomic <- feat_names_adj[1:3]
    gram_diff <- feat_names_adj[8:10]
    cooccurrence <- c("Animal", "Human gut", "Human skin", "Human oral",
                      "Freshwater", "Marine", "Soil", "Wastewater")
    
    for (i in 1:nrow(plot_data_a)){
      if (plot_data_a$Feature[i] %in% genomic){
        plot_data_a$var_categ[i] <- "Genetic incompatibility"
      }
      else if (plot_data_a$Feature[i] %in% cooccurrence){
        plot_data_a$var_categ[i] <- "Environmental co-occurrence"
      } 
      else if (plot_data_a$Feature[i] %in% gram_diff){
        plot_data_a$var_categ[i] <- "Gram staining"
      }
      else {
        plot_data_a$var_categ[i] <- "Other"
      }
    }
    
    # Plot Fig 3a
    fig_s9a <- ggplot(plot_data_a, aes(x=factor(Feature, levels=Feature[order(mean_imp, decreasing = TRUE)]), 
                                      y=mean_imp, fill=as.factor(var_categ))) + 
      geom_bar(position="dodge", stat="identity", color="black") +
      geom_errorbar(aes(ymin=mean_imp-sd_imp, ymax=mean_imp+sd_imp), 
                    width=.3, position=position_dodge(.9)) +
      labs(title = "", x = "", y = "MeanDecreaseAccuracy", tag = "a") +
      scale_fill_manual(values=feature_category_palette, 
                        name = "Feature group",
                        breaks=c("Genetic incompatibility",
                                 "Environmental co-occurrence",
                                 "Gram staining",
                                 "Other")) +
      theme_minimal() +
      geom_text(aes(label=significance),
                position=position_dodge(width=0.75), 
                hjust=0.6, vjust = 0.3, size=10) +
      fig_a_theme
    
  }
  
  close(pb)
}

# Compile results for Fig 3b,c,d
plot_data <- df_total[df_total$Feature != "Gene.class",]

mech_names_adj <- c("All", "AAC", "APH", "Class A, C, D", "Class B", "Erm", "Tet efflux", "Tet RPG")

for (i in 1:length(mechanisms)) {
  plot_data$class[plot_data$class == mechanisms[i]] <- mech_names_adj[i]
}

for (i in 1:length(feat_names_to_adjust)) {
  plot_data$Feature[plot_data$Feature == feat_names_to_adjust[i]] <- feat_names_adj[i]
}

plot_data$class <- factor(plot_data$class, levels=mech_names_adj)
plot_data$class <- fct_rev(factor(plot_data$class))

genetic_compatibility <- plot_data[plot_data$Feature %in% genomic,]
env_cooccurrence <- plot_data[plot_data$Feature %in% cooccurrence,]
gram <- plot_data[plot_data$Feature %in% gram_diff,]

# Plot Fig 3b,c,d
fig_s9b <- ggplot(genetic_compatibility, aes(fill=class, y=mean_imp, x=reorder(Feature, -mean_imp))) + 
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_imp-sd_imp, ymax=mean_imp+sd_imp), 
                width=.3, position=position_dodge(.9)) +
  xlab(NULL) +
  ylab("MeanDecreaseAccuracy") +
  labs(title = "Genetic incompatibility", tag = "b") +
  geom_text(aes(label=significance),
            position=position_dodge(width=0.9), 
            hjust=-5*sign(genetic_compatibility$mean_imp), 
            vjust = 0.75, size=10) +
  scale_fill_manual(values=mechanism_palette, 
                    name = "Resistance mechanism",
                    breaks = rev(levels(genetic_compatibility$class))) + 
  coord_flip() +
  theme_minimal() +
  fig_bcd_theme

fig_s9c <- ggplot(env_cooccurrence, aes(fill=class, y=mean_imp, x=reorder(Feature, -mean_imp))) + 
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_imp-sd_imp, ymax=mean_imp+sd_imp), 
                width=.3, position=position_dodge(.9)) +
  xlab(NULL) +
  ylab("MeanDecreaseAccuracy") +
  labs(title = "Environmental co-occurrence", tag = "c") +
  geom_text(aes(label=significance),
            position=position_dodge(width=0.9), 
            hjust=-5*sign(env_cooccurrence$mean_imp), 
            vjust = 0.75, size=10) +
  scale_fill_manual(values=mechanism_palette, 
                    name = "Resistance mechanism",
                    breaks = rev(levels(env_cooccurrence$class))) + 
  coord_flip() +
  theme_minimal() +
  fig_bcd_theme

fig_s9d <- ggplot(gram, aes(fill=class, y=mean_imp, x=reorder(Feature, -mean_imp))) + 
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_imp-sd_imp, ymax=mean_imp+sd_imp), 
                width=.3, position=position_dodge(.9)) +
  xlab(NULL) +
  ylab("MeanDecreaseAccuracy") +
  labs(title = "Gram staining", tag = "d") +
  geom_text(aes(label=significance),
            position=position_dodge(width=0.9), 
            hjust=-5*sign(gram$mean_imp), 
            vjust = 0.75, size=10) +
  scale_fill_manual(values=mechanism_palette, 
                    name = "Resistance mechanism",
                    breaks = rev(levels(gram$class))) + 
  coord_flip() +
  theme_minimal() +
  fig_bcd_theme

pdf("sup_fig_9.pdf", width = 22.5, height = 18.75)
fig_s9a + fig_s9b + fig_s9c + fig_s9d
dev.off()

#----------------- Supplementary Fig 10 ---------------------------------------
# Load packages
library(igraph)
library(stringr)
library(GGally)
library(ggplot2)
library(tidyr)
library(data.table)
library(dplyr)

# Functions
plot_cooccurrence_network <- function(graph, environment) {
  size <- c()
  for (name in V(graph)$name) {
    size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
  }
  
  genus <- c()
  for (name in V(graph)$name) {
    genus <- c(genus, str_split(name, " ")[[1]][1])
  }
  
  V(graph)$label <- taxonomy[genus,3]
  
  for (i in 1:length(V(graph)$label)) {
    if (is.na(V(graph)$label[i])) {
      V(graph)$label[i] <- taxonomy[taxonomy$family == genus[i],3][1]
    }
  }
  
  env_list <- c("Animal", "Human", "Soil", "Water", "Wastewater")
  
  imp <- edge_attr(graph)[[grep(environment, env_list)]]
  
  imp[is.na(imp)] <- 0
  imp[imp < 0.01] <- 0
  texture <- c()
  alph <- c()
  
  for (i in 1:length(imp)) {
    if (is.na(imp[i])) {
      texture <- c(texture, 0)
      alph <- c(alph, 0)
    }
    else if (imp[i] == 0) {
      texture <- c(texture, 3)
      alph <- c(alph, 0)
    }
    else if (imp[i] > 0 && imp[i] <= 0.01) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.01 && imp[i] <= 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
  }
  
  E(graph)$lty <- texture
  E(graph)$color <- "#bdbdbd"
  
  node_color <- c("#41ab5d", "#005a32", "#9ecae1", "#4292c6", "#084594", 
                  "#feb24c", "#dd3497", "#fc9272", "#ef3b2c", "#99000d")
  
  names(node_color) <-c("Actinomycetes", "Coriobacteriia", "Bacilli", "Clostridia",
                        "Erysipelotrichia", "Bacteroidia", "Epsilonproteobacteria",
                        "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")
  
  esize <- sqrt(imp*20)
  esize[esize == 0] <- 0.5
  
  shape <- rep(15, length(V(graph)$name))
  shape[grep(" ", V(graph)$name)] <- 19
  
  set.seed(444)
  
  p <- ggnet2(graph, 
              mode = "kamadakawai", 
              color=V(graph)$label, 
              palette = node_color, 
              edge.alpha = alph, 
              edge.size = esize, 
              alpha = 1, size = size, 
              edge.color = E(graph)$color, 
              edge.lty = E(graph)$lty, 
              legend.size = 9, 
              layout.exp = 0.1,
              node.shape = shape) +
    geom_text(aes(label = V(graph)$name), color = "black", size=3, fontface = "bold") +
    guides(size = FALSE) +
    coord_fixed()
}

plot_all_combinations <- function(graph) {
  size <- c()
  for (name in V(graph)$name) {
    size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
  }
  
  genus <- c()
  for (name in V(graph)$name) {
    genus <- c(genus, str_split(name, " ")[[1]][1])
  }
  
  V(graph)$label <- taxonomy[genus,3]
  
  for (i in 1:length(V(graph)$label)) {
    if (is.na(V(graph)$label[i])) {
      V(graph)$label[i] <- taxonomy[taxonomy$family == genus[i],3][1]
    }
  }
  
  node_color <- c("#41ab5d", "#005a32", "#9ecae1", "#4292c6", "#084594", 
                  "#feb24c", "#dd3497", "#fc9272", "#ef3b2c", "#99000d")
  
  names(node_color) <-c("Actinomycetes", "Coriobacteriia", "Bacilli", "Clostridia",
                        "Erysipelotrichia", "Bacteroidia", "Epsilonproteobacteria",
                        "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")
  
  shape <- rep(15, length(V(graph)$name))
  shape[grep(" ", V(graph)$name)] <- 19
  
  set.seed(444)
  
  p <- ggnet2(graph, 
              mode = "kamadakawai", 
              color=V(graph)$label, 
              palette = node_color, 
              edge.alpha = 0.6, 
              edge.size = 0.75, 
              alpha = 1, size = size, 
              edge.color = "#4575b4", 
              legend.size = 9, 
              layout.exp = 0.1,
              node.shape = shape) +
    geom_text(aes(label = V(graph)$name), color = "black", size=3, fontface = "bold") +
    guides(size = FALSE) +
    coord_fixed()
}


# Load data
taxonomy <- data.frame(fread("data_for_figures/assembly_full_lineage.txt", header=FALSE)) %>%
  select(-c(1,ncol(.))) %>%
  na.omit() %>%
  unique() %>%
  `colnames<-`(c("kingdom", "phylum", "class", "order", "family", "genus")) %>%
  `rownames<-`(.$genus)

co_occurrence <- read.table("data_for_figures/species_cooccurrence.txt", header = TRUE, sep="\t") %>%
  .[!(is.na(.$Species2) | is.na(.$Species1)),]

co_occurrence <- co_occurrence[order(co_occurrence$n.true, decreasing = TRUE), ]

prio_list <- c("Mycobacterium tuberculosis", "Staphylococcus aureus", "Escherichia coli", 
               "Streptococcus pneumoniae",  "Klebsiella pneumoniae", "Pseudomonas aeruginosa", 
               "Acinetobacter baumannii", "Enterobacter cloacae", "Streptococcus agalactiae", 
               "Enterococcus faecalis", "Enterococcus faecium", "Salmonella enterica", 
               "Streptococcus pyogenes", "Neisseria meningitidis","Campylobacter jejuni", 
               "Vibrio cholerae", "Clostridioides difficile")

# Aggregate species to genus-level
sp <- unique(c(co_occurrence$Species1, co_occurrence$Species2))
for (name in sp) {
  if (name %in% prio_list) {
    next
  }
  else {
    co_occurrence$Species1[co_occurrence$Species1 == name] <- str_split(name, " ")[[1]][1]
    co_occurrence$Species2[co_occurrence$Species2 == name] <- str_split(name, " ")[[1]][1]
  }
}

# Filter non-frequent pairs
df <- co_occurrence[co_occurrence$n.true >= 5, ]

sp <- unique(c(df$Species1, df$Species2))
df2 <- data.frame(matrix(nrow = 0, ncol = 8))
included <- c()

for (name in sp) {
  subset <- df[df$Species1 == name | df$Species2 == name,]
  
  for (i in 1:nrow(subset)) {
    if (paste(subset$Species1[i], subset$Species2[i]) %in% included || paste(subset$Species2[i], subset$Species1[i]) %in% included) {
      next
    }
    
    else {
      included <- c(included, paste(subset$Species1[i], subset$Species2[i]))
      subset2 <- subset[subset$Species1 == subset$Species1[i] & subset$Species2 == subset$Species2[i] | subset$Species1 == subset$Species2[i] & subset$Species2 == subset$Species1[i],]
      df2 <- rbind(df2, c(subset$Species1[i], subset$Species2[i], max(na.omit(subset2$Animal)), max(na.omit(subset2$Human)), max(na.omit(subset2$Soil)), max(na.omit(subset2$Water)), max(na.omit(subset2$Wastewater)), sum(na.omit(subset2$n.true))))
    }
  }
}

colnames(df2) <- colnames(df)

for (i in 3:8){
  df2[,i] <- as.numeric(df2[,i])
}

df2[df2 == -Inf] <- 0

# Generate initial graph
g <- graph.data.frame(d = df2, directed = FALSE)

# Aggregate uncommon genera to family level
size <- c()
for (name in V(g)$name) {
  size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
}

prio_list <- data.frame(V1 = c(prio_list, head(V(g)$name[order(size, decreasing = TRUE)], 15)))

sp <- unique(c(df2$Species1, df2$Species2))
for (name in sp) {
  if (name %in% prio_list$V1) {
    next
  }
  else {
    co_occurrence$Species1[co_occurrence$Species1 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    co_occurrence$Species2[co_occurrence$Species2 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    
    df2$Species1[df2$Species1 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    df2$Species2[df2$Species2 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
  }
}

sp <- unique(c(df2$Species1, df2$Species2))
df3 <- data.frame(matrix(nrow = 0, ncol = 8))
included <- c()

for (name in sp) {
  subset <- df2[df2$Species1 == name | df2$Species2 == name,]
  
  for (i in 1:nrow(subset)) {
    if (paste(subset$Species1[i], subset$Species2[i]) %in% included || paste(subset$Species2[i], subset$Species1[i]) %in% included) {
      next
    }
    
    else {
      included <- c(included, paste(subset$Species1[i], subset$Species2[i]))
      subset2 <- subset[subset$Species1 == subset$Species1[i] & subset$Species2 == subset$Species2[i] | subset$Species1 == subset$Species2[i] & subset$Species2 == subset$Species1[i],]
      df3 <- rbind(df3, c(subset$Species1[i], subset$Species2[i], max(na.omit(subset2$Animal)), max(na.omit(subset2$Human)), max(na.omit(subset2$Soil)), max(na.omit(subset2$Water)), max(na.omit(subset2$Wastewater)), sum(na.omit(subset2$n.true))))
    }
  }
}

colnames(df3) <- colnames(df2)

for (i in 3:8){
  df3[,i] <- as.numeric(df3[,i])
}

# Generate final network
g <- graph.data.frame(d = df3, directed = FALSE)

# Plot and save subfigures
fig_s10a <- plot_cooccurrence_network(g, "Animal")

pdf("sup_fig_10a.pdf", width = 15, height = 10)
plot(fig_s10a)
dev.off()

fig_s10b <- plot_cooccurrence_network(g, "Soil")

pdf("sup_fig_10b.pdf", width = 15, height = 10)
plot(fig_s10b)
dev.off()

fig_s10c <- plot_cooccurrence_network(g, "Water")

pdf("sup_fig_10c.pdf", width = 15, height = 10)
plot(fig_s10c)
dev.off()

fig_s10d <- plot_all_combinations(g)

pdf("sup_fig_10d.pdf", width = 15, height = 10)
plot(fig_s10d)
dev.off()

#----------------- Supplementary Fig 11 ---------------------------------------
# Load packages
library(igraph)
library(randomcoloR)
library(stringr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyr)
library(intergraph)
library(taxizedb)
library(data.table)
library(dplyr)

# Functions
plot_cooccurrence_network <- function(graph, environment) {
  size <- c()
  for (name in V(graph)$name) {
    size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
  }
  
  genus <- c()
  for (name in V(graph)$name) {
    genus <- c(genus, str_split(name, " ")[[1]][1])
  }
  
  V(graph)$label <- taxonomy[genus,3]
  
  for (i in 1:length(V(graph)$label)) {
    if (is.na(V(graph)$label[i])) {
      V(graph)$label[i] <- taxonomy[taxonomy$family == genus[i],3][1]
    }
  }
  
  env_list <- c("Animal", "Human", "Soil", "Water", "Wastewater")
  
  imp <- edge_attr(graph)[[grep(environment, env_list)]]
  
  imp[is.na(imp)] <- 0
  imp[imp < 0.01] <- 0
  texture <- c()
  alph <- c()
  
  for (i in 1:length(imp)) {
    if (is.na(imp[i])) {
      texture <- c(texture, 0)
      alph <- c(alph, 0)
    }
    else if (imp[i] == 0) {
      texture <- c(texture, 3)
      alph <- c(alph, 0)
    }
    else if (imp[i] > 0 && imp[i] <= 0.01) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.01 && imp[i] <= 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
    else if (imp[i] > 0.03) {
      texture <- c(texture, 1)
      alph <- c(alph, 0.6)
    }
  }
  
  E(graph)$lty <- texture
  E(graph)$color <- "#bdbdbd"
  
  node_color <- c("#41ab5d", "#005a32", "#9ecae1", "#4292c6", "#084594", 
                  "#feb24c", "#dd3497", "#fc9272", "#ef3b2c", "#99000d")
  
  names(node_color) <-c("Actinomycetes", "Coriobacteriia", "Bacilli", "Clostridia",
                        "Erysipelotrichia", "Bacteroidia", "Epsilonproteobacteria",
                        "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria")
  
  esize <- sqrt(imp*20)
  esize[esize == 0] <- 0.5
  
  shape <- rep(15, length(V(graph)$name))
  shape[grep(" ", V(graph)$name)] <- 19
  
  set.seed(444)
  
  p <- ggnet2(graph, 
              mode = "kamadakawai", 
              color=V(graph)$label, 
              palette = node_color, 
              edge.alpha = alph, 
              edge.size = esize, 
              alpha = 1, size = size, 
              edge.color = E(graph)$color, 
              edge.lty = E(graph)$lty, 
              legend.size = 9, 
              layout.exp = 0.1,
              node.shape = shape) +
    geom_text(aes(label = V(graph)$name), color = "black", size=3, fontface = "bold") +
    guides(size = FALSE) +
    coord_fixed()
}

# Load data
taxonomy <- data.frame(fread("data_for_figures/assembly_full_lineage.txt", header=FALSE)) %>%
  select(-c(1,ncol(.))) %>%
  na.omit() %>%
  unique() %>%
  `colnames<-`(c("kingdom", "phylum", "class", "order", "family", "genus")) %>%
  `rownames<-`(.$genus)

co_occurrence <- read.table("data_for_figures/species_cooccurrence.txt", header = TRUE, sep="\t") %>%
  .[!(is.na(.$Species2) | is.na(.$Species1)),]

co_occurrence <- co_occurrence[order(co_occurrence$n.true, decreasing = TRUE), ]

prio_list <- c("Acinetobacter guillouiae", "Acinetobacter gyllenbergii",
               "Aeromonas caviae", "Aeromonas allosaccharophila", 
               "Aeromonas sanarelli", "Aeromonas media", "Citrobacter freundii",
               "Morganella morganii", "Enterobacter asburiae", "Enterobacter cloacae",
               "Hafnia alvei", "Hafnia paralvei", "Klebsiella pneumoniae",
               "Acinetobacter radioresistens", "Shewanella xiamenensis",
               "Acinetobacter baumannii", "Kluyvera georgina", "Kluyvera ascorbata",
               "Rheinheimera pacifica", "Moraxella pluranimalium", "Shewanella frigidimarina",
               "Shewanella algae", "Vibrio splendidus", "Kluyvera georgiana", "Leclercia adecarboxylata")

# Aggregate species to genus-level
sp <- unique(c(co_occurrence$Species1, co_occurrence$Species2))
for (name in sp) {
  if (name %in% prio_list) {
    next
  }
  else {
    co_occurrence$Species1[co_occurrence$Species1 == name] <- str_split(name, " ")[[1]][1]
    co_occurrence$Species2[co_occurrence$Species2 == name] <- str_split(name, " ")[[1]][1]
  }
}

# Filter non-frequent pairs
df <- co_occurrence[co_occurrence$n.true >= 5, ]

sp <- unique(c(df$Species1, df$Species2))
df2 <- data.frame(matrix(nrow = 0, ncol = 8))
included <- c()

for (name in sp) {
  subset <- df[df$Species1 == name | df$Species2 == name,]
  
  for (i in 1:nrow(subset)) {
    if (paste(subset$Species1[i], subset$Species2[i]) %in% included || paste(subset$Species2[i], subset$Species1[i]) %in% included) {
      next
    }
    
    else {
      included <- c(included, paste(subset$Species1[i], subset$Species2[i]))
      subset2 <- subset[subset$Species1 == subset$Species1[i] & subset$Species2 == subset$Species2[i] | subset$Species1 == subset$Species2[i] & subset$Species2 == subset$Species1[i],]
      df2 <- rbind(df2, c(subset$Species1[i], subset$Species2[i], max(na.omit(subset2$Animal)), max(na.omit(subset2$Human)), max(na.omit(subset2$Soil)), max(na.omit(subset2$Water)), max(na.omit(subset2$Wastewater)), sum(na.omit(subset2$n.true))))
    }
  }
}

colnames(df2) <- colnames(df)

for (i in 3:8){
  df2[,i] <- as.numeric(df2[,i])
}

df2[df2 == -Inf] <- 0

# Generate initial graph
g <- graph.data.frame(d = df2, directed = FALSE)

# Aggregate uncommon genera to family level
size <- c()
for (name in V(g)$name) {
  size <- c(size, sum(co_occurrence[co_occurrence$Species1 == name | co_occurrence$Species2 == name, "n.true"]))
}

prio_list <- data.frame(V1 = c(prio_list, head(V(g)$name[order(size, decreasing = TRUE)], 15)))

sp <- unique(c(df2$Species1, df2$Species2))
for (name in sp) {
  if (name %in% prio_list$V1) {
    next
  }
  else {
    co_occurrence$Species1[co_occurrence$Species1 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    co_occurrence$Species2[co_occurrence$Species2 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    
    df2$Species1[df2$Species1 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
    df2$Species2[df2$Species2 == name] <- taxonomy[str_split(name, " ")[[1]][1], 5]
  }
}

sp <- unique(c(df2$Species1, df2$Species2))
df3 <- data.frame(matrix(nrow = 0, ncol = 8))
included <- c()

for (name in sp) {
  subset <- df2[df2$Species1 == name | df2$Species2 == name,]
  
  for (i in 1:nrow(subset)) {
    if (paste(subset$Species1[i], subset$Species2[i]) %in% included || paste(subset$Species2[i], subset$Species1[i]) %in% included) {
      next
    }
    
    else {
      included <- c(included, paste(subset$Species1[i], subset$Species2[i]))
      subset2 <- subset[subset$Species1 == subset$Species1[i] & subset$Species2 == subset$Species2[i] | subset$Species1 == subset$Species2[i] & subset$Species2 == subset$Species1[i],]
      df3 <- rbind(df3, c(subset$Species1[i], subset$Species2[i], max(na.omit(subset2$Animal)), max(na.omit(subset2$Human)), max(na.omit(subset2$Soil)), max(na.omit(subset2$Water)), max(na.omit(subset2$Wastewater)), sum(na.omit(subset2$n.true))))
    }
  }
}

colnames(df3) <- colnames(df2)

for (i in 3:8){
  df3[,i] <- as.numeric(df3[,i])
}

# Generate final network
g <- graph.data.frame(d = df3, directed = FALSE)

# Plot and save subfigures
fig_s10a <- plot_cooccurrence_network(g, "Animal")

pdf("sup_fig_11a.pdf", width = 15, height = 10)
plot(fig_s10a)
dev.off()

fig_s10b <- plot_cooccurrence_network(g, "Human")

pdf("sup_fig_11b.pdf", width = 15, height = 10)
plot(fig_s10b)
dev.off()

fig_s10c <- plot_cooccurrence_network(g, "Soil")

pdf("sup_fig_11c.pdf", width = 15, height = 10)
plot(fig_s10c)
dev.off()

fig_s10d <- plot_cooccurrence_network(g, "Human")

pdf("sup_fig_11d.pdf", width = 15, height = 10)
plot(fig_s10d)
dev.off()

fig_s10e <- plot_cooccurrence_network(g, "Soil")

pdf("sup_fig_11e.pdf", width = 15, height = 10)
plot(fig_s10e)
dev.off()

#----------------- Supplementary Fig 12 ---------------------------------------
# Load packages
library(ggplot2)
library(corrplot)
library(stringr)
library(data.table)

count_otus <- function(otus) {
  otus_split <- strsplit(otus, ";")
  n <- length(unique(otus_split[[1]]))
  
  return(as.numeric(n))
}

calc_overall_coccurrence <- function(intra1, intra2, inter){
  if (intra1 == 0 & intra2 == 0) {
    overall_cooccurrence <- inter/1e-10
  }
  else {
    overall_cooccurrence <- inter/((intra1 + intra2)/2)
  }
  
  return(as.numeric(overall_cooccurrence))
}

prep_for_correlation_analysis <- function(filepath) {
  table <- data.frame(fread(filepath)) %>%
    na.omit() %>%
    select(., Genome_5mer_distance, Gene_genome_5mer_distance,
           Genome_size_difference, Gram_stain_difference,
           Animal, Human, Soil, Water, Wastewater)
  
  for (i in 1:nrow(table)){
    if (table$Gram_stain_difference[i] == "NN"){
      table$Gram_stain_difference[i] <- 0
    }
    else if (table$Gram_stain_difference[i] == "PP"){
      table$Gram_stain_difference[i] <- 1
    }
    else if (table$Gram_stain_difference[i] == "NP"){
      table$Gram_stain_difference[i] <- 2
    }
  }
  
  table$Gram_stain_difference <- as.numeric(table$Gram_stain_difference)
  
  colnames(table) <- c("Genome 5mer distance", "Genome-gene 5mer distance",
                                    "Difference in genome size", "Gram stain", "Animal",
                                    "Human", "Soil", "Water", "Wastewater")
  
  return(table)
}

# Analysis for Sup. Fig 12a
observed_transfers <- prep_for_correlation_analysis("data_for_figures/observed_horizontal_transfers.txt")

M <- cor(observed_transfers, method = "spearman")
testRes = cor.mtest(observed_transfers, conf.level = 0.99)

col <- colorRampPalette(c("#313695","#4575b4","#74add1", "#abd9e9","#e0f3f8",
                          "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"))

# Plot and save subfigure
pdf("sup_fig_12a.pdf", width=8.5, height=7)
corrplot(M, p.mat = testRes$p, sig.level = 1, order = 'hclust', 
         method="color", col=col(200), addCoef.col = "black", tl.col="black", 
         tl.srt=55, insig = "blank", diag = FALSE, type="upper", tl.cex = 1)
dev.off()

# Analysis for Sup. Fig 12a
observed_transfers <- prep_for_correlation_analysis("data_for_figures/randomized_transfers1.txt")

M <- cor(observed_transfers, method = "spearman")
testRes = cor.mtest(observed_transfers, conf.level = 0.99)

col <- colorRampPalette(c("#313695","#4575b4","#74add1", "#abd9e9","#e0f3f8",
                          "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"))

# Plot and save subfigure
pdf("sup_fig_12b.pdf", width=8.5, height=7)
corrplot(M, p.mat = testRes$p, sig.level = 1, order = 'hclust', 
         method="color", col=col(200), addCoef.col = "black", tl.col="black", 
         tl.srt=55, insig = "blank", diag = FALSE, type="upper", tl.cex = 1)
dev.off()

#----------------- Supplementary Fig 13 ---------------------------------------
# Load packages
library(caret)
library(dplyr)
library(forcats)
library(data.table)
library(randomForest)
library(pROC)
library(data.table)
library(utils)
library(patchwork)

# Functions
remove_redundant_events <- function(events) {
  uniq_events <- unique(events[, c("Node", "Gene.class")])
  subtable <- data.frame(matrix(nrow = 0, ncol = ncol(events)))
  colnames(subtable) <- colnames(events)
  
  for (i in seq_len(nrow(uniq_events))) {
    redundant_events <- events[events$Node == uniq_events$Node[i] & events$Gene.class == uniq_events$Gene.class[i], ]
    if (nrow(redundant_events) > 10) {
      subtable <- rbind(subtable, redundant_events[sample(seq_len(nrow(redundant_events)), 10), ])
    } else {
      subtable <- rbind(subtable, redundant_events)
    }
  }
  return(subtable)
}

subsample_null_events <- function(true, null) {
  gene_class <- unique(true$Gene.class)
  downsampled_events <- data.frame()
  tmp <- data.frame()
  
  for (class in gene_class) {
    n <- sum(true$Gene.class == class)
    subset <- null[null$Gene.class == class, ]
    if (nrow(subset) > n) {
      selection <- sample(seq_len(nrow(subset)), n)
      downsampled_events <- rbind(downsampled_events, subset[selection, ])
      tmp <- rbind(tmp, subset[-selection, ])
    } else {
      downsampled_events <- rbind(downsampled_events, subset)
    }
  }
  
  n_left <- nrow(true) - nrow(downsampled_events)
  downsampled_events <- rbind(downsampled_events, tmp[sample(seq_len(nrow(tmp)), n_left), ])
  
  return(downsampled_events)
}

format_input_data <- function(true_data, null_data, mechanism) {
  true_data$Transfer <- 1
  null_data$Transfer <- 0
  
  if (mechanism != "all") {
    if (mechanism == "class_A_C_D") {
      true_data <- true_data[true_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
      null_data <- null_data[null_data$Gene.class %in% c("class_A", "class_C", "class_D_1", "class_D_2"), ]
    }
    
    else {
      true_data <- true_data[grep(mechanism, true_data$Gene.class), ]
      null_data <- null_data[grep(mechanism, null_data$Gene.class), ]
    }
  }
  
  true_data <- remove_redundant_events(true_data)
  if (nrow(null_data) > nrow(true_data)) {
    null_data <- subsample_null_events(true_data, null_data)
  }
  
  input_data <- rbind(true_data, null_data)
  input_data$Transfer <- as.factor(input_data$Transfer)
  input_data$Gene.class <- as.factor(input_data$Gene.class)
  
  input_data <- input_data %>%
    mutate(
      NN = as.integer(Gram_stain_difference == "NN"),
      PP = as.integer(Gram_stain_difference == "PP"),
      NP = as.integer(Gram_stain_difference == "NP")
    )
  
  return(input_data)
}

get_distance_top_left <- function(sens_spec) {
  d <- sqrt(sum((sens_spec - c(1,1))^2))
  return(d)
}

# Plotting settings
fig_theme <- theme(
  axis.text.x = element_text(size = rel(1.3)),
  axis.text.y = element_text(size = rel(1.3)),
  axis.title = element_text(size = rel(1.3)),
  legend.title = element_text(size = rel(1.2)),
  legend.text = element_text(size = rel(1)),
  plot.title = element_text(size = rel(3), face = "bold")
)

aggregated_palette <- c(
  "All (n = 10)" ="#000000",
  "AAC"='#a50026',
  "APH" ="#d73027",
  "Class A, C, D"='#fdae61',
  "Class B"="#fee090",
  "Erm" = "#ffffff",
  "Tet efflux" = "#74add1",
  "Tet RPG" = "#313695"
)

observed_transfers <- data.frame(fread("data_for_figures/observed_transfers_cooc_divided.txt")) %>% 
  subset(select = -c(Header1, Header2)) %>%
  na.omit()
mechanisms <- c("all", "aac", "aph", "class_A_C_D", "class_B", "erm", "tet_efflux", "tet_rpg")
performance <- data.frame(class = mechanisms, 
                          auc_mean = rep(0, 8), auc_sd = rep(0, 8), 
                          sens_mean = rep(0, 8), sens_sd = rep(0, 8), 
                          spec_mean = rep(0, 8), spec_sd = rep(0, 8))

# Analysis
for (class in mechanisms) {
  perf <- data.frame(matrix(ncol=4, nrow=10))
  colnames(perf) <- c("auc", "acc", "sens", "spec")
  
  print(paste(c("Generating", class, "models"), collapse=" "))
  
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  for (i in 1:10){
    randomized_transfers <- read.delim(paste(c("data_for_figures/randomized_transfers", i, "_cooc_divided", ".txt"), collapse = "")) %>%
      na.omit()
    
    input_data <- format_input_data(observed_transfers, randomized_transfers, class) %>% 
      select(Gene.class, NN, PP, NP, Genome_5mer_distance,
             Gene_genome_5mer_distance, Genome_size_difference,
             Animal, Human.gut, Human.skin, Human.oral, Freshwater, Marine, 
             Soil, Wastewater, Transfer)
    
    set.seed(1)
    
    train_index <- createDataPartition(input_data$Transfer, times = 1, p = 0.7, list = FALSE)
    train_set <- input_data %>% slice(train_index)
    test_set <- input_data %>% slice(-train_index)
    
    rf_model <- randomForest(Transfer ~ .,
                             data = train_set, 
                             ntree = 500, 
                             type = "classification", 
                             na.action=na.omit, 
                             importance = TRUE)
    
    p_hat_test <- as.data.frame(predict(rf_model, newdata = test_set, type = "prob"))
    suppressMessages(roc_test <- roc(test_set$Transfer ~ p_hat_test$`1`, plot=FALSE, smooth = FALSE))
    best_point <- coords(roc_test, "best", best.method = "closest.topleft")
    
    y_hat_test <- predict(rf_model, newdata = test_set, cutoff = c(1 - as.numeric(best_point["threshold"][[1]])[1], as.numeric(best_point["threshold"][[1]])[1]))
    
    conf <- caret::confusionMatrix(y_hat_test, test_set$Transfer, positive='1')
    
    perf$auc[i] <- auc(roc_test)
    perf$acc[i] <- conf$overall[[1]]
    perf$sens[i] <- conf$byClass[[1]]
    perf$spec[i] <- conf$byClass[[2]]
    
    if (class == "all") {
      suppressMessages(assign(paste(c("roc", i), collapse = ""), roc_test))
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  performance$auc_mean[performance$class == class] <- mean(perf$auc)
  performance$auc_sd[performance$class == class] <- sd(perf$auc)
  performance$sens_mean[performance$class == class] <- mean(perf$sens)
  performance$sens_sd[performance$class == class] <- sd(perf$sens)
  performance$spec_mean[performance$class == class] <- mean(perf$spec)
  performance$spec_sd[performance$class == class] <- sd(perf$spec)
  
}

# Compile results
sensitivities <- list(roc1$sensitivities, roc2$sensitivities,
                      roc3$sensitivities, roc4$sensitivities,
                      roc5$sensitivities, roc6$sensitivities,
                      roc7$sensitivities, roc8$sensitivities,
                      roc9$sensitivities, roc10$sensitivities)

specificities <- list(roc1$specificities, roc2$specificities,
                      roc3$specificities, roc4$specificities,
                      roc5$specificities, roc6$specificities,
                      roc7$specificities, roc8$specificities,
                      roc9$specificities, roc10$specificities)

roc_data <- data.frame(Sensitivity = tapply(unlist(sensitivities), sequence(lengths(sensitivities)), mean), 
                       Specificity = tapply(unlist(specificities), sequence(lengths(specificities)), mean))

euc_dist <- apply(roc_data, 1, get_distance_top_left)

plot_data <- data.frame(Performance = c(performance$auc_mean, performance$sens_mean, performance$spec_mean),
                        SD = c(performance$auc_sd, performance$sens_sd, performance$spec_sd),
                        Mechanism = rep(performance$class, 3),
                        Metric = c(rep("AUROC", 8), rep("Sensitivity", 8), rep("Specificity", 8)))

mech_names_adj <- c("All (n = 10)", "AAC", "APH", "Class A, C, D", "Class B", "Erm", "Tet efflux", "Tet RPG")

for (i in 1:length(mechanisms)) {
  plot_data$Mechanism[plot_data$Mechanism == mechanisms[i]] <- mech_names_adj[i]
}

plot_data$Mechanism <- as.factor(plot_data$Mechanism)

# Plot figures
fig_s13a <- ggplot() +
  geom_line(aes(x = roc1$specificities, y = roc1$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc2$specificities, y = roc2$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc3$specificities, y = roc3$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc4$specificities, y = roc4$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc5$specificities, y = roc5$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc6$specificities, y = roc6$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc7$specificities, y = roc7$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc8$specificities, y = roc8$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc9$specificities, y = roc9$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc10$specificities, y = roc10$sensitivities), col = "#4575b4", linewidth=0.8, alpha = 0.4) +
  geom_line(aes(x = roc_data$Specificity, y = roc_data$Sensitivity), col = "#000000", linewidth=1.1) +
  geom_point(aes(x = roc_data$Specificity[euc_dist==min(euc_dist)], y = roc_data$Sensitivity[euc_dist==min(euc_dist)]), size = 4.5, color= "#000000") +
  scale_x_reverse(labels = function(x) 1 - x) +
  theme_minimal() +
  ggtitle("a") +
  ylab("Sensitivity") +
  xlab("1-Specificity") +
  fig_theme

fig_s13b <- ggplot(data = plot_data, aes(x = Metric, y = Performance, 
                                       fill = factor(Mechanism, levels = mech_names_adj))) +
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_errorbar(aes(ymin=Performance-SD, ymax=Performance+SD), width=.3, position=position_dodge(.9)) +
  theme_minimal() +
  scale_fill_manual(values=aggregated_palette, name = "Resistance mechanism") + 
  ggtitle("b") +
  xlab("") +
  ylab("Performance") +
  ylim(0,1) +
  fig_theme

# Compile and save Fig 2
pdf("sup_fig_13.pdf", width = 15, height = 10.5)
fig_s13a + fig_s13b
dev.off()