# Copyright (C) 2024 Troels Mouritzen
# 
# Author: Troels Mouritzen
# Affiliation: Aarhus University
# Email: twm@mbg.au.dk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Description: This script calculates HAPk scores
# on a SNP table in non-overlapping windows along the genome
# up to a certain amount of clusters
# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

main <- function() {
  library(tidyverse)
  
  args <- commandArgs(T) #snpfile, window size (in kb), jobname and max number of clusters
  
  set.seed(1234)
  
  print_args(args)
  window_size <- as.integer(args[3])
  max_clusters <- as.integer(args[5])
  
  snpfile <- read.csv(args[1])
  
  GM <- snpfile %>%
    select(1:2) %>%
    rename(Chromosome = CHR, Position = POS) %>%
    mutate(SNP = str_c("SNP", Chromosome, Position, sep = "_")) %>%
    as_tibble()
  
  GD <- snpfile %>%
    select(-c(1:2)) %>%
    t() %>%
    as_tibble()
  
  colnames(GD) <- GM$SNP
  
  GM %>%
    make_kmeans_matrix(GD, window_size, max_clusters) %>%
    mutate(Start = Position, End = Position + window_size*10^3,
           Position = as.integer(Position + 0.5*window_size*10^3)) %>%
    write_csv(str_c(args[4], "_cluster_result.csv"))
}

print_args <- function(args) {
  print(str_c("Snpfile = ", args[2]))
  print(str_c("Window size = ", args[3], "kb"))
  print(str_c("Jobname = ", args[4]))
  print(str_c("Number of clusters = ", args[5]))
}

myKmeans <- function(data, n_k) {
  if (n_k < 2) {
    return()
  }
  save <- rep(NA, n_k - 1)
  clu <- 2:n_k
  for (k in 2:n_k) {
    try_R2 <- rep(NA, 5)
    for (l in 1:50) {
      km <- kmeans(data, centers = k, iter.max = 1000)
      try_R2[l] <- km$betweenss/km$totss
    }
    save[k - 1] <- max(try_R2)
  }
  return(save)
}

find_locations <- function(map, ranges, windows, current, i, j) {
  locations <- c()
  while(T) {
    if (map$Chromosome[current] > ranges$Chromosome[i]) {
      break
    }
    if (map$Chromosome[current] < ranges$Chromosome[i]) {
      current = current + 1
      next
    }
    
    if (map$Position[current] < windows[j]) {
      locations <- c(locations, current)
      current = current + 1
    } else {
      break
    }
  }
  return(list(locations, current))
}



make_kmeans_matrix <- function(map, genotypes, window_size, n_k) {
  map %>%
    group_by(Chromosome) %>%
    summarise(min = min(Position), max = max(Position)) -> ranges
  
  current <- 1
  chrom <- c()
  window_start <- c()
  clusters <- c()
  R2 <- c()
  n_SNPs <- c()
  for (i in 1:nrow(ranges)) {
    cat("\n")
    print(str_c("Chromosome ", as.character(ranges$Chromosome[i])))
    windows <- seq(from = ranges$min[i], to = ranges$max[i], by = as.integer(window_size)*10^3)
    print(str_c("Number of windows = ", as.character(length(windows))))
    progress <- seq.int(from = 1, to = length(windows), length.out = 20) %>%
      as.integer()
    progress_print <- seq(from = 0, to = 100, by = 5) %>%
      as.character() %>%
      str_c("%..")
    progress_print[20] <- "100%"
    for (j in 2:length(windows)) {
      if (is.na(windows[j])) {
        break
      }
      if (j %in% progress) {
        cat(progress_print[match(j, progress)])
      }
      l_c <- find_locations(map, ranges, windows, current, i, j)
      locations <- l_c[[1]]
      current <- l_c[[2]]
      if (length(locations) == 0) {
        next
      }
      my_geno <- genotypes[,locations]
      n_clusters <- min(n_k, as.integer(length(locations)))
      if (n_clusters < 2) {
        next
      }
      R2_temp <- myKmeans(my_geno, n_clusters)
      
      chrom <- c(chrom, rep(ranges$Chromosome[i], n_clusters - 1))
      window_start <- c(window_start, rep(windows[j - 1], n_clusters - 1))
      n_SNPs <- c(n_SNPs, rep(length(locations), n_clusters - 1))
      clusters <- c(clusters, 2:n_clusters)
      R2 <- c(R2, R2_temp)
    }
  }
  
  tibble(Chromosome = chrom, Position = window_start, n_SNPs = n_SNPs, clusters = clusters, R2 = R2) %>%
    mutate(clusters = str_c("V", clusters)) %>%
    pivot_wider(id_cols = c("Chromosome", "Position", "n_SNPs"), names_from = "clusters", values_from = "R2") %>%
    return()
}

main()
