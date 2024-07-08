#!/usr/bin/env Rscript
compute.cooccurrence.gwmc <- function(otus) {
  otu1 <- otus[1]
  otu2 <- otus[2]

  if (otu1 == otu2) {
    cooccurrence <- 0
  }
  else {
    cooccurrence <- NA
    extracted_otus <- count_table[c(otu1, otu2), ]
    sum_counts <- colSums(extracted_otus)

    if (sum(sum_counts) == 0) {
      cooccurrence <- 0
    }

    else {
      cooccurrence <- sum(sum_counts == 2) / length(sum_counts)
      names(cooccurrence) <- NULL
    }
  }

  output <- data.frame(OTU1 = otu1, OTU2 = otu2, Wastewater = cooccurrence[1])
  rownames(output) <- NULL
  return(output)
}

compile.cooccurrence.gwmc <- function(otu1, otu2) {
  otus_crossed <- tidyr::crossing(x = otu1, y = otu2)
  otus_crossed <- data.frame(OTU1 = otus_crossed$x, OTU2 = otus_crossed$y)

  otus_crossed$OTU1 <- as.character(otus_crossed$OTU1)
  otus_crossed$OTU2 <- as.character(otus_crossed$OTU2)

  otus_crossed <- otus_crossed[otus_crossed$OTU1 != otus_crossed$OTU2, ]

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

  otu_combinations <- apply(otu_combinations, 1, compute.cooccurrence.gwmc)
  tmp <- as.data.frame(matrix(nrow = 0, ncol = 6))

  for (l in 1:length(otu_combinations)) {
    tmp <- rbind(tmp, otu_combinations[[l]])
  }

  results <- sum(tmp[, 3]) / nrow(tmp)
  names(results) <- NULL

  return(results)
}

analyze.cooccurrence.gwmc <- function(input_data) {
    otu1 <- unique(str_split(as.character(input_data[1]), ";")[[1]]) %>% .[. != ""]
    otu2 <- unique(str_split(as.character(input_data[2]), ";")[[1]]) %>% .[. != ""]

    if (is.na(input_data[1]) || is.na(input_data[2])) {
      results <- NA
    }

    else if (length(otu1) == length(otu2) && input_data[1] == input_data[2]) {
      results <- NA
    }

    else {
      results <- compile.cooccurrence.gwmc(otu1, otu2)
    }
    return(results)
}
