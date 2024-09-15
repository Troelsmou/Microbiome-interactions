# Copyright (C) 2024 Turgut Akyol
# 
# Author: Turgut Akyol
# Affiliation: Aarhus University
# Email: tya@mbg.au.dk
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

# Description: This script performs mantel tests on
# a set of microbiome GWA SNPs and compares it to
# a distance matrix from the geographical locations
# or population designations of the samples.
# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

########## Mantel test for microbiome GWA SNPs ##########
library(vegan)
library(ggplot2)

# GWAS SNPs =====================================
ssig = read.delim('microbiomeSNPs.csv', sep = ',')
ssig = unique(ssig$significantSNPs)

# Other SNPs from Troels ========================
unsig = read.delim('SNPs_to_generate_random_sets.txt', sep = '\t', header = F)
unsig = paste0('SNP_', unsig$V1)
unsig = setdiff(unsig, ssig)  # remove 78

# All SNP data ==================================
allSNPs = read.delim('./20220330_lotus_snps.csv', sep = ',')  # Unzip the data first.
row.names(allSNPs) = paste('SNP', allSNPs$CHR, allSNPs$POS, sep = '_')
allSNPs = allSNPs[, c(-1, -2)]
ssigDF = allSNPs[ssig, ]
set.seed(44); unsigDF = allSNPs[sample(unsig, 36), ]

# Significant SNPs with LD ======================
popo1 = readRDS('popo1.rds')
popo2 = readRDS('popo2.rds')
popo3 = readRDS('popo3.rds')

# Distance matrices of GWAS SNPs
dp1 = dist(t(ssigDF[popo1,]))
dp2 = dist(t(ssigDF[popo2,]))
dp3 = dist(t(ssigDF[popo3,]))

# PCA of GWAS SNPs
pca1 = prcomp(dp1)
pca2 = prcomp(dp2)
pca3 = prcomp(dp3)

# Nothern-southern information ==================
pop = read.delim('./pop.csv', sep = ';', row.names = 1)
row.names(pop)[1] = 'Gifu_'
sss = intersect(row.names(pop), row.names(pca1$x))
pop = pop[row.names(pca1$x), ,drop = FALSE]

# Mantel test with nothern-southern
mm = model.matrix( ~ pop$Population)[,-1, drop = FALSE]
d4 = dist(mm)

mantel(dp1, d4, permutations = 0)
mantel(dp2, d4, permutations = 0)
mantel(dp3, d4, permutations = 0)

# Exact geographical locations ==================
lloc2 = read.delim('./Lotus_accessions_location_2.csv', sep = ',', row.names = 1)
identical(row.names(lloc2), row.names(pop))  # TRUE
identical(row.names(lloc2), colnames(ssigDF))  # TRUE

# Mantel test
lloc3 = lloc2[complete.cases(lloc2), c('lon', 'lat')]
mantel(dist(t(ssigDF[popo1, row.names(lloc3)])), dist(lloc3), permutations = 0)
mantel(dist(t(ssigDF[popo2, row.names(lloc3)])), dist(lloc3), permutations = 0)
mantel(dist(t(ssigDF[popo3, row.names(lloc3)])), dist(lloc3), permutations = 0)
mantel(dist(t(unsigDF[, row.names(lloc3)])), dist(lloc3), permutations = 0)
