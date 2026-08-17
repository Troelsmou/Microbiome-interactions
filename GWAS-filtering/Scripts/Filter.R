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

# Description: This script filters raw GWAS results
# according to the criteria used in the associated publication.

# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]
pacman::p_load(tidyverse, data.table)

################################################################################
# Functions
################################################################################

get_window_starts_n <- function(chromosomes, n) {
  
  start_indices <- 1L
  current_chr <- chromosomes[1]
  count <- 1L
  
  for (i in seq_along(chromosomes)[-1]) {
    if (chromosomes[i] != current_chr) {
      # New chromosome always starts a window
      start_indices <- c(start_indices, i)
      current_chr <- chromosomes[i]
      count <- 1L
    } else {
      count <- count + 1L
      if (count > n) {
        # New window when exceeding n SNPs
        start_indices <- c(start_indices, i)
        count <- 1L
      }
    }
  }
  
  return(start_indices)
}

Meff.simpleM <- function(cor_r, pca_cutoff_perc) {
  num_of_snps <- ncol(cor_r)
  if(num_of_snps == 0) return(0)
  
  eigen_values <- eigen(cor_r, only.values = TRUE)$values
  sum_eigen_values <- sum(eigen_values)
  eigen_values <- sort(eigen_values, decreasing = TRUE)
  
  M_eff_G <- 1
  for(k in 1:num_of_snps) {
    temp <- sum(eigen_values[1:k])/sum_eigen_values
    if(temp >= pca_cutoff_perc) {
      M_eff_G <- k
      break
    }
  }
  return(M_eff_G)
}

compute_window_meff <- function(genotype, window_starts, pca_cutoff_perc) {
  if(length(window_starts) == 0) return(numeric(0))
  
  n_windows <- length(window_starts)
  meff_results <- numeric(n_windows)
  total_snps <- ncol(genotype)
  window_ends = c(window_starts[-1] - 1, total_snps)
  
  for(i in seq_along(window_starts)) {
    start_idx <- window_starts[i]
    end_idx <- window_ends[i]
    
    if (i %% 50 == 0) {
      print(paste("Processing window", i, "of", n_windows))
    }
    
    # Handle empty windows gracefully
    if(start_idx > end_idx) {
      meff_results[i] <- 0
      next
    }
    
    window_data <- genotype[, start_idx:end_idx, drop = FALSE]
    cor_matrix <- cor(window_data, use = "pairwise.complete.obs")
    
    meff_results[i] <- Meff.simpleM(cor_matrix, pca_cutoff_perc)
  }
  
  return(meff_results)
}

determine_LD_blocks_snpstats <- function(order, snp_mat, chromosomes, positions, pos_threshold, R2_threshold = 0.8) {
  require(snpStats)
  require(data.table)
  LD_block = integer(ncol(snp_mat))
  ld_cur = 1L
  snp_mat <- as(snp_mat, "SnpMatrix")
  
  if (is.unsorted(chromosomes)) {
    stop("SNPs must be sorted by chromosome and position.")
  }
  
  dt <- data.table(chr=chromosomes, pos=positions, idx=seq_along(positions))
  chr_ranges <- dt[, .(start = min(idx), end = max(idx)), by = chr]
  
  for (i in order) {
    if (LD_block[i] != 0L) {
      next
    }
    # Fast chromosomal boundary check
    current_chr <- chromosomes[i]
    chr_start <- chr_ranges$start[chr_ranges$chr == current_chr]
    chr_end <- chr_ranges$end[chr_ranges$chr == current_chr]
    
    # Binary search for position window
    current_pos <- positions[i]
    current_positions <- positions[chr_start:chr_end]
    lower_bound = current_pos - pos_threshold
    upper_bound = current_pos + pos_threshold
    window_start <- chr_start + findInterval(lower_bound, current_positions, checkSorted = F, checkNA = F)
    window_end <- chr_start + findInterval(upper_bound, current_positions, checkSorted = F, checkNA = F) - 1
    
    candidate_idx <- window_start:window_end
    candidate_idx <- candidate_idx[LD_block[candidate_idx] == 0L]
    mysnp <- which(i == candidate_idx)
    if (length(candidate_idx) != 0L) {
      ld_mat <- ld(snp_mat[, candidate_idx], stats = "R.squared", depth = length(candidate_idx) - 1, symmetric = T)
      in_ld <- which(ld_mat[mysnp, ] > R2_threshold)
      LD_block[candidate_idx[in_ld]] <- ld_cur
    }
    LD_block[i] <- ld_cur
    ld_cur <- ld_cur + 1L
  }
  tibble(SNP = colnames(snp_mat), LD_block = LD_block)
}

# Efficiently checks if a SNP is within or near any gene (using binary search)
snp_in_gene_binary <- function(pos, gene_starts, gene_ends, tol) {
  left <- 1
  right <- length(gene_starts)
  while (left <= right) {
    mid <- floor((left + right) / 2)
    if (pos < gene_starts[mid] - tol) {
      right <- mid - 1
    } else if (pos > gene_ends[mid] + tol) {
      left <- mid + 1
    } else {
      return(TRUE)
    }
  }
  return(FALSE)
}

snp_in_gene_vectorized_binary <- Vectorize(snp_in_gene_binary, vectorize.args = "pos")

snp_in_gene_allchromosomes_binary <- function(snps, genes, tol) {
  log_vector <- logical(nrow(snps))
  for (chrom in unique(snps$Chromosome)) {
    print(chrom)
    snp_idx <- which(snps$Chromosome == chrom)
    gene_idx <- which(genes$Chromosome == chrom)
    gene_starts <- genes$Start[gene_idx]
    gene_ends <- genes$End[gene_idx]
    log_vector[snp_idx] <- snp_in_gene_vectorized_binary(
      snps$Position[snp_idx],
      gene_starts,
      gene_ends,
      tol
    )
  }
  snps$in_gene <- log_vector
  return(snps)
}

variance_of_fitted <- function(names_vec, slopes_vec, data_df) {
  # Check input lengths
  if(length(names_vec) != length(slopes_vec)) {
    stop("names_vec and slopes_vec must be the same length.")
  }
  # Calculate variance for each (name, slope) pair
  sapply(seq_along(names_vec), function(i) {
    col_name <- names_vec[i]
    slope <- slopes_vec[i]
    fitted_values <- slope * data_df[[col_name]]
    var(fitted_values, na.rm = TRUE)
  })
}

my_variance_explained <- function(pheno, SNP, trait, effect, GD) {
  pheno_var_exp <- pheno %>%
    group_by(DescriptionOfTrait) %>%
    summarize(var = var(Score))
  
  sig_var_exp <- variance_of_fitted(SNP, effect, GD)
  
  sig_var_exp <- tibble(SNP = SNP, Trait = trait, Var_explained = sig_var_exp)
  
  sig_var_exp <- sig_var_exp %>%
    left_join(pheno_var_exp, by = c("Trait" = "DescriptionOfTrait")) %>%
    mutate(Var_explained = Var_explained / var) %>%
    pull(Var_explained)
  return(sig_var_exp)
}


################################################################################
# Preprocessing
################################################################################

GD <- fread("../Data/Lotus_GD.csv") %>%
  as_tibble()
GM <- fread("../Data/Lotus_GM.csv") %>%
  as_tibble()

GM$SNP <- str_c("SNP_", GM$SNP)
colnames(GD) <- c("taxa", GM$SNP)

genes <- read_delim("20210713_Lj_Gifu_v1.3_predictedGenes.gff3", delim = "\t",
                    col_names = F, skip = 10)
genes %>%
  filter(X3 == "gene") %>%
  dplyr::select(X1, X4, X5) %>%
  rename(Chromosome = X1, Start = X4, End = X5) %>%
  mutate(Chromosome = as.numeric(str_remove(Chromosome, "LjG1.1_chr"))) %>%
  na.omit() %>%
  arrange(Chromosome, Start) -> genes

permuted <- read_csv("20250408_Lotus_permutation_results.csv")
sig <- read_csv("20240508_Lotus_rarefied_GWAS_Results_sig.csv") %>%
  filter(P.value < 10^-6 & Chromosome %in% 1:6) %>%
  arrange(Chromosome, Position) %>%
  snp_in_gene_allchromosomes_binary(genes, 10^3) %>%
  filter(in_gene) %>%
  select(!in_gene)

permuted %>%
  right_join(sig) %>%
  select(!c(permutation_count, exceeding_permutations_count, Trial)) %>%
  mutate(Replicate = str_extract(Trait, "(?<=_)([1-3]|ave)$")) -> full

full %>%
  filter(is.na(Replicate)) -> non_bact
full %>%
  filter(!is.na(Replicate)) %>%
  mutate(Trait = str_remove(Trait, "_([1-3]|ave)$")) %>%
  mutate(Trait_type = case_when(
    str_detect(Trait, "(MDS|pca)") ~ "Dimensionality reduction",
    str_detect(Trait, "ceae") ~ "Bacterial families",
    str_detect(Trait, "Lj") ~ "Isolates",
    TRUE ~ "Contaminants")) -> bact

################################################################################
# Simple M
################################################################################

new_GD <- as.matrix(GD[-1])
maf <- colSums(new_GD) / (nrow(new_GD) * 2)
maf[maf > 0.5] <- 1 - maf[maf > 0.5]
new_GD <- new_GD[, maf > 0.05]
new_GM <- GM[maf > 0.05,] %>%
  arrange(Chromosome, Position)
gene_vector <- snp_in_gene_allchromosomes_binary(new_GM, genes, 0)$in_gene
new_GM <- new_GM[gene_vector, ]
new_GD <- new_GD[, gene_vector]

windows <- get_window_starts_n(new_GM$Chromosome, 159)
simplem <- compute_window_meff(new_GD, windows, 0.995)

meff <- sum(simplem)
cutoff <- 1 - (1 - 0.05)^(1 / meff)

################################################################################
# Bacterial GWAS filtering
################################################################################

bact %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> bact_snps

snp_names <- bact_snps$SNP
chromosomes <- bact_snps$Chromosome
positions <- bact_snps$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- bact %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(bact) -> bact

bact %>%
  group_by(LD_block, Trait, Replicate) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(bact) %>%
  distinct(SNP, Trait, Replicate, Method) %>%
  mutate(Significance = "Trait and Replicate")-> trait_rep_sig

bact %>%
  group_by(LD_block, Trait) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(bact) %>%
  distinct(SNP, Trait, Replicate, Method) %>%
  anti_join(trait_rep_sig) %>%
  mutate(Significance = "Trait") -> trait_sig

bact %>%
  group_by(LD_block) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(bact) %>%
  distinct(SNP, Trait, Replicate, Method) %>%
  anti_join(trait_rep_sig) %>%
  anti_join(trait_sig) %>%
  mutate(Significance = "General") %>%
  bind_rows(trait_rep_sig, trait_sig) %>%
  left_join(bact) -> all_sig

################################################################################
# Permuted Bacterial GWAS
################################################################################

permuted_permuted <- read_csv("20250808_Lotus_nonsense_GWAS_permutation_results.csv")
sig_permuted <- read_csv("20250714_Lotus_permutation_GWAS_Results_sig.csv") %>%
  filter(P.value < 10^-6 & Chromosome %in% 1:6) %>%
  arrange(Chromosome, Position) %>%
  snp_in_gene_allchromosomes_binary(genes, 0) %>%
  filter(in_gene) %>%
  select(!in_gene)

permuted_permuted %>%
  right_join(sig_permuted) %>%
  select(!c(permutation_count, exceeding_permutations_count, Trial)) %>%
  mutate(Permutation = str_extract(Trait, "(?<=Permuted_)[0-9]+$"),
         Replicate = str_extract(Trait, "(?<=_)([1-3]|ave)(?=_Permuted_[0-9]+$)")) -> full_permuted

bact_permuted <- full_permuted %>%
  filter(!is.na(Replicate) & !is.na(permuted_p_value))

bact_permuted %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> bact_snps_permuted

snp_names <- bact_snps_permuted$SNP
chromosomes <- bact_snps_permuted$Chromosome
positions <- bact_snps_permuted$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- bact_permuted %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(bact_permuted) -> bact_permuted

bact_permuted %>%
  group_by(LD_block) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(bact_permuted) -> permuted_sig

permuted_sig %>%
  rename(GEMMA_p.value = original_p_value, Permuted_p.value = permuted_p_value,
         Effect = effect, QTL = LD_block) %>%
  select(QTL, SNP, Chromosome, Position, Trait, Replicate, Permutation,
         Method, P.value, GEMMA_p.value, Permuted_p.value, Effect, MAF) %>%
  write_csv("20250825_Permutation_filtered_nonsense_significants.csv")

abiotic_permuted <- full_permuted %>%
  filter(is.na(Replicate) & !is.na(permuted_p_value))

abiotic_permuted %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> abiotic_snps_permuted

snp_names <- abiotic_snps_permuted$SNP
chromosomes <- abiotic_snps_permuted$Chromosome
positions <- abiotic_snps_permuted$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- abiotic_permuted %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(abiotic_permuted) -> abiotic_permuted

abiotic_permuted %>%
  group_by(LD_block) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(abiotic_permuted) -> permuted_sig

permuted_sig %>%
  rename(GEMMA_p.value = original_p_value, Permuted_p.value = permuted_p_value,
         Effect = effect, QTL = LD_block) %>%
  select(QTL, SNP, Chromosome, Position, Trait, Permutation,
         Method, P.value, GEMMA_p.value, Permuted_p.value, Effect, MAF) %>%
  write_csv("20250902_Permutation_filtered_nonsense_abiotic_significants.csv")

################################################################################
# Flowering GWAS filtering
################################################################################

non_bact %>%
  filter(Trait %in% c("1st_ft", "2015_16_2nd_year_FT", "2017_FT", "2018_FT",
                      "2nd_ft", "2nd_ft_period", "FP_2014", "FT_DK")) -> flowering

flowering %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> flowering_snps

snp_names <- flowering_snps$SNP
chromosomes <- flowering_snps$Chromosome
positions <- flowering_snps$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- flowering %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>% 
  left_join(flowering) -> flowering

flowering %>%
  group_by(LD_block, Trait) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(flowering) %>%
  mutate(Significance = "Exact Trait", LD_block = str_c(LD_block, "_flowering")) -> flowering_sig

################################################################################
# Temperature GWAS filtering
################################################################################

non_bact %>%
  filter(Trait %in% c("Min_temp", "Altitude", "Mean_temp", "OW_2014", 
                      "OW_2015", "OW_2016", "OW_2017")) -> temperature

temperature %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> temperature_snps

snp_names <- temperature_snps$SNP
chromosomes <- temperature_snps$Chromosome
positions <- temperature_snps$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- temperature %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(temperature) -> temperature

temperature %>%
  group_by(LD_block, Trait) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(temperature) %>%
  mutate(Significance = "Exact Trait", LD_block = str_c(LD_block, "_temperature")) -> temperature_sig

################################################################################
# Seed traits
################################################################################

non_bact %>%
  filter(Trait %in% c("Seed_perimeter", "Seed_weight", "Seed_size", "Seed_length",
                      "Seed_width", "Seed_LWR", "Seed_circularity")) -> seed

seed %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> seed_snps

snp_names <- seed_snps$SNP
chromosomes <- seed_snps$Chromosome
positions <- seed_snps$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- seed %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(seed) -> seed

seed %>%
  group_by(LD_block, Trait) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(seed) %>%
  distinct(SNP, Trait) %>%
  mutate(Significance = "Exact Trait") -> seed_sig

seed %>%
  group_by(LD_block) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(seed) %>%
  distinct(SNP, Trait) %>%
  anti_join(seed_sig) %>%
  mutate(Significance = "Other Seed Trait") %>%
  bind_rows(seed_sig) %>%
  left_join(seed) %>%
  mutate(LD_block = str_c(LD_block, "_seed"))-> seed_sig

################################################################################
# Salt
################################################################################

non_bact %>%
  filter(Trait %in% c("K_salt", "N_salt", "K_contr", "N_contr")) -> salt

salt %>%
  distinct(SNP, Chromosome, Position) %>%
  arrange(Chromosome, Position) -> salt_snps

snp_names <- salt_snps$SNP
chromosomes <- salt_snps$Chromosome
positions <- salt_snps$Position

GD %>%
  select(all_of(snp_names)) %>%
  as.matrix() -> GD_matrix_gwas

order <- salt %>%
  count(SNP) %>%
  arrange(-n) %>%
  pull(SNP) %>%
  match(colnames(GD_matrix_gwas))

determine_LD_blocks_snpstats(order, GD_matrix_gwas, chromosomes, positions, 5*10^4, R2_threshold = 0.8) %>%
  left_join(salt) -> salt

salt %>%
  group_by(LD_block, Trait) %>%
  summarize(sig = any(permuted_p_value < cutoff), .groups = "drop") %>%
  filter(sig) %>%
  select(!sig) %>%
  left_join(salt) %>%
  mutate(Significance = "Exact Trait", LD_block = str_c(LD_block, "_salt")) -> salt_sig

################################################################################
# Non_bacterial together
################################################################################

bind_rows(salt_sig, flowering_sig, seed_sig, temperature_sig) %>%
  rename(GEMMA_p.value = original_p_value, Permuted_p.value = permuted_p_value,
         Effect = effect) -> non_bact_sig

################################################################################
# Var explained and write
################################################################################

pheno <- read_csv("Lotus_Johan_ave_rep_family_mds_rarefied.csv")

all_sig$VarianceExplained <- my_variance_explained(pheno, all_sig$SNP, str_c(all_sig$Trait, "_", all_sig$Replicate),
                      all_sig$effect, GD)

all_sig %>%
  group_by(LD_block) %>%
  summarize(MethodsPerQTL = length(unique(Method))) %>%
  left_join(all_sig) %>%
  rename(GEMMA_p.value = original_p_value, Permuted_p.value = permuted_p_value,
         Effect = effect, QTL = LD_block, TraitType = Trait_type) %>%
  select(QTL, SNP, Chromosome, Position, Trait, TraitType, Replicate, Significance,
         Method, MethodsPerQTL, P.value, GEMMA_p.value, Permuted_p.value, Effect, MAF, VarianceExplained) %>%
  write_csv("20250514_Permutation_filtered_significants.csv")

non_bact_sig$VarianceExplained <- my_variance_explained(pheno, non_bact_sig$SNP, non_bact_sig$Trait,
                                                          non_bact_sig$Effect, GD)
non_bact_sig %>%
  group_by(LD_block) %>%
  summarize(MethodsPerQTL = length(unique(Method))) %>%
  left_join(non_bact_sig) %>%
  rename(QTL = LD_block) %>%
  select(QTL, SNP, Chromosome, Position, Trait, Significance,
         Method, MethodsPerQTL, P.value, GEMMA_p.value, Permuted_p.value, Effect, MAF, VarianceExplained) %>%
  write_csv("20250514_Permutation_filtered_significants_non_bact.csv")
