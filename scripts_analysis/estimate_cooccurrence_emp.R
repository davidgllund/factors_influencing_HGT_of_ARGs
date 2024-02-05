#!/usr/bin/env Rscript
calc.presence.proportion <- function(key, presence, otus) {
  sample_subset <- presence[colnames(otus) %in% c("X.OTU.ID", paste("X", samples_from_env[[key]], sep = ""))]

  proportion <- sum(sample_subset == 2) / length(sample_subset)

  return(proportion)
}

compute.cooccurrence.emp <- function(otus) {
  otu1 <- otus[1]
  otu2 <- otus[2]

  if (otu1 == otu2) {
    cooccurrence <- rep(0, times = 4)
  }
  else {
    cooccurrence <- rep(NA, times = 4)
    extracted_otus <- count_table[c(as.character(otu1), as.character(otu2)), ]
    sum_counts <- colSums(extracted_otus)

    if (sum(sum_counts) == 0) {
      cooccurrence <- rep(0, times = 4)
    }

    else {
      cooccurrence <- sapply(environment_categories, calc.presence.proportion, presence = sum_counts, otus = extracted_otus)
      names(cooccurrence) <- NULL
    }
  }

  output <- data.frame(OTU1 = otu1, OTU2 = otu2, Animal = cooccurrence[1], Human = cooccurrence[2], Soil = cooccurrence[3], Water = cooccurrence[4])
  rownames(output) <- NULL
  return(output)
}

compile.cooccurrence.emp <- function(otu1, otu2) {
  otus_crossed <- tidyr::crossing(x = otu1, y = otu2)
  otus_crossed <- data.frame(OTU1 = otus_crossed$x, OTU2 = otus_crossed$y)

  otus_crossed$OTU1 <- as.character(otus_crossed$OTU1)
  otus_crossed$OTU2 <- as.character(otus_crossed$OTU2)

  otus_crossed <- otus_crossed[otus_crossed$OTU1 != otus_crossed$OTU2,]

  included <- c()
  otu_combinations <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(otu_combinations) <- c("OTU1", "OTU2")

  for (l in 1:nrow(otus_crossed)) {
    if (paste(otus_crossed$OTU1[l], otus_crossed$OTU2[l]) %in% included || paste(otus_crossed$OTU2[l], otus_crossed$OTU1[l]) %in% included) {
      next
    }
    else{
      included <- c(included, paste(otus_crossed$OTU1[l], otus_crossed$OTU2[l]))
      otu_combinations <- rbind(otu_combinations, otus_crossed[l,])
    }
  }

  otu_combinations <- apply(otu_combinations, 1, compute.cooccurrence.emp)
  tmp <- as.data.frame(matrix(nrow = 0, ncol = 6))

  for (l in 1:length(otu_combinations)) {
    tmp <- rbind(tmp, otu_combinations[[l]])
  }

  results <- colSums(tmp[, 3:6]) / nrow(tmp)
  names(results) <- NULL

  return(results)
}

analyze.cooccurrence.emp <- function(input_data) {
    otu1 <- unique(str_split(as.character(input_data[1]), ";")[[1]])
    otu2 <- unique(str_split(as.character(input_data[2]), ";")[[1]])

    if (is.na(otu1) || is.na(otu2)) {
      results <- rep(NA, times = 4)
    }

    else if (length(otu1) == length(otu2) && otu1 == otu2) {
      results <- rep(NA, times = 4)
    }

    else {
      results <- compile.cooccurrence.emp(otu1, otu2)
    }
    return(results)
}

