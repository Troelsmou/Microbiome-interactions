########### Pick a microbiome GWA SNP from each LD block ##########

######### This part is from TWM ##########
# Load packages
pacman::p_load(tidyverse, data.table)

# Load files
GD <- fread("Lotus_GD.csv", sep = ",") %>% as_tibble()  # Unzip the data first.
gwas_snps <- read_csv("microbiomeSNPsExtended.csv")  # 78 SNPs

# Transform data
GD_matrix <- GD[, -1] %>%
  as.matrix()
gwas_snps %>%
  pull(SNP) %>%
  str_replace("SNP_", "") %>%
  unique() -> snp_names
GD_matrix_gwas <- GD_matrix[, snp_names]

# Define function to determine LD blocks
determine_LD_blocks <- function(snp_mat, R2_threshold = 0.8) {
  # Calculate LD
  LD <- snp_mat %>%
    cor() %>%
    abs()
  LD[lower.tri(LD, diag = T)] <- 0
  which(LD > R2_threshold, arr.ind = T) -> LD_pairs
  # Determine LD blocks
  indices = 1:ncol(snp_mat)
  found = rep(FALSE, length(indices))
  LD_blocks <- list()
  for (i in indices) {
    if (found[i]) {
      next
    }
    relevant_cells = which(LD_pairs[, 1] == i | LD_pairs[, 2] == i)
    c(LD_pairs[relevant_cells, 1], LD_pairs[relevant_cells, 2])-> relevant_indices
    ld_block <- unique(c(i, relevant_indices))
    found[ld_block] <- TRUE
    LD_blocks <- c(LD_blocks, list(ld_block))
  }
  return(LD_blocks)
}

# So at this point, GD_matrix_gwas is a matrix with a column for each unique SNP
# It might work with unsorted positions, but it is probably best to sort them first
LD_blocks <- determine_LD_blocks(GD_matrix_gwas)
# LD blocks is a list with each element being a vector of column indices in GD_matrix_gwas
# corresponding to the SNPs in that LD block
n <- length(LD_blocks) # The number of LD blocks / the number of random SNPs in the whole genome sampling
# We can then sample 1 from each LD block into a vector of column indices

#################################################
# Pick 3 sets of microbiome GWA SNPs with size of 36
set.seed(57)
popo = list()
for(i in 1:100) {
  popo[[i]] = sapply(LD_blocks, function(x) if(length(x) > 1) sample(x, 1) else x)
}
sapply(popo, length)  # all 36

popo1 = popo[[1]]
popo2 = popo[[2]]
popo3 = popo[[3]]

popo1 = paste0('SNP_', colnames(GD_matrix_gwas)[popo1])
popo2 = paste0('SNP_', colnames(GD_matrix_gwas)[popo2])
popo3 = paste0('SNP_', colnames(GD_matrix_gwas)[popo3])

saveRDS(popo1, file = 'popo1.rds')
saveRDS(popo2, file = 'popo2.rds')
saveRDS(popo3, file = 'popo3.rds')
