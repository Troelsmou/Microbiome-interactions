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
################################################################################
# Functions
################################################################################
#Helper function
snp_in_gene <- function(pos, gene_starts, gene_ends, tol) {
  return(any(pos >= gene_starts - tol & pos <= gene_ends + tol))
}
#Helper function
snp_in_gene_vectorized <- Vectorize(snp_in_gene, vectorize.args = "pos")

# Function for determining if a SNP is in a gene
# Or generally if all rows in a dataframe with columns Chromosome and position is within
# a window in a dataframe defined by columns Chromosome, start and end +- tol
# Where tol is the length outside the window that is included in the window
snp_in_gene_allchromosomes <- function(snps, genes, tol) {
  log_vector <- c()
  snps %>%
    arrange(Chromosome, Position) -> snps
  for (chrom in unique(snps$Chromosome)) {
    gene_starts <- genes$start[genes$Chromosome == chrom]
    gene_ends <- genes$end[genes$Chromosome == chrom]
    log_vector <- c(log_vector, snp_in_gene_vectorized(snps$Position[snps$Chromosome == chrom],
                                                       gene_starts, gene_ends, tol))
  }
  snps %>%
    mutate(in_gene = log_vector) %>%
    return()
}

################################################################################
# Load packages
################################################################################
library(tidyverse)
################################################################################
# Load files
################################################################################

genes <- read_delim("20210713_Lj_Gifu_v1.3_predictedGenes.gff3", delim = "\t",
                    col_names = F, skip = 10) %>%
  filter(X3 == "gene") %>%
  select(X1, X4, X5) %>%
  rename(Chromosome = X1, start = X4, end = X5) %>%
  mutate(Chromosome = as.numeric(str_remove(Chromosome, "LjG1.1_chr"))) %>%
  na.omit()

sig <- read_csv("20240508_Lotus_rarefied_GWAS_Results_sig.csv.gz")

################################################################################
# GWAS Filters
################################################################################

Methods <- c("FarmCPU", "Gemma", "ISIS EM-BLASSO") #Which methods to include
p_value <- 10^-6 #P.value filter
filter_genes <- TRUE #Whether to only include SNPs within x bp of gene
gene_range = 1000 #Range to include genes (in both directions)
snps_per_trait <- 200 #Max number of significant unique SNPs per trait
round_position <- -4 #How much to round position of snps of later filters.
# 0, -1, -2, -3, -4 corresponds to rounding to nearest 1, 10, 100, 1000, 10000 bp
reps_per_trait <- 2 #How many methods/trials/snps per rounded position per trait to keep
Trials_plus_traits_across_traits <- 3 #How many trials plus traits per rounded position to keep (minimum is 2)

# Remove snps below p value of 10^-6 and BLINK rows
sig <- sig %>%
  filter(Method %in% Methods) %>%
  filter(P.value < p_value)

# Remove SNPs not within 1000 bp of a gene
sig %>%
  snp_in_gene_allchromosomes(genes, gene_range) %>%
  filter(!filter_genes & in_gene) %>%
  select(!in_gene) -> sig_gene

# Split SNPs into bacterial gwas and non_bacterial_gwas
sig_gene %>%
  mutate(Trial = str_extract(Trait, "(?<=_)(ave)|[1-3]$"),
         Trait = str_remove(Trait, "_((ave)|[1-3])$")) -> sig_gene
sig_gene %>%
  filter(is.na(Trial) & Trait != "OW_2017") -> other_gwas
sig_gene %>%
  filter(!is.na(Trial)) -> sig_bact

# Remove traits with more than 200 significant SNPs
sig_bact %>%
  count(Trait, SNP) %>%
  count(Trait) %>%
  filter(n < snps_per_trait) %>%
  select(!n) %>%
  left_join(sig_bact) -> sig_bact_traits

# Rounds positions to nearest 10000
sig_bact_traits %>%
  mutate(new_position = round(Position, round_positions)) -> try
# Remove SNPs with no significant SNPs at the same rounded position for the same trait
# (i.e. remove SNPs that are not significant for another replicate or method or closeby snp)
try %>%
  count(Chromosome, new_position, Trait) %>%
  filter(n >= reps_per_trait) %>%
  select(!n) %>%
  left_join(try) -> sig_bact_traits_snps

# Remove SNPs that are not significant at the same rounded position for all traits
# for a different trial or method
sig_bact_traits_snps %>%
  group_by(Chromosome, new_position) %>%
  summarize(n = length(unique(Trait)) + length(unique(Trial))) %>%
  filter(n >= reps_across_traits) %>%
  left_join(sig_bact_traits_snps) %>%
  ungroup() -> sig_bact_traits_snps_consistent

# Rounds position to nearest 10000 and groups traits
other_gwas %>%
  mutate(new_pos = round(Position, -4),
         old_trait = Trait, 
         Trait = case_when(Trait %like% "temp" | Trait %like% "OW" ~ "Temperature",
                           Trait %like% "Seed" ~ "Seed",
                           Trait %like% "ft" | Trait %like% "FT" | Trait %like% "FP" ~ "Flowering",
                           Trait %like% "contr" | Trait %like% "salt" ~ "Salt",
                           .default = Trait)) -> other_gwas

# Filter SNPs that are not significant twice for the same trait in the same rounded position
# Can be 2 methods, 2 replicates or nearby snp
other_gwas %>%
  count(Chromosome, new_pos, Trait) %>%
  filter(n > 1) %>%
  left_join(other_gwas) %>%
  mutate(Trait = old_trait) %>%
  select(!c(in_gene, old_trait, Trial)) -> other_gwas_filtered

# write files
write_csv(sig_bact_traits_snps_consistent, "20240508_Lotus_rarefied_GWAS_Results_sig_filtered_bacteria.csv")
write_csv(other_gwas_filtered, "20240508_Lotus_rarefied_GWAS_Results_sig_filtered_other.csv")
