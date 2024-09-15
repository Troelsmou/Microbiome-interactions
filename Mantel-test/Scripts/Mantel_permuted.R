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
# a set of random snps of size equal to the number of
# LD blocks in the significant SNPs and to make a null
# distribution of test scores.
# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

########## Mantel test for 1,000 sets ##########
library(vegan)

popo1 = readRDS('popo1.rds')

ssig = read.delim('microbiomeSNPs.csv', sep = ',')
ssig = unique(ssig$significantSNPs)

unsig = read.delim('SNPs_to_generate_random_sets.txt', sep = '\t', header = F)
unsig = paste0('SNP_', unsig$V1)
unsig = setdiff(unsig, ssig)

allSNPs = read.delim('20220330_lotus_snps.csv', sep = ',')  # Unzip the data first.
row.names(allSNPs) = paste('SNP', allSNPs$CHR, allSNPs$POS, sep = '_')
allSNPs = allSNPs[, c(-1, -2)]

mm__ = read.delim('./Lotus_accessions_location_2.csv', sep = ',', row.names = 1)
mm__ = mm__[complete.cases(mm__),]
pipi = dist(mm__[, c('lat', 'lon')])

allSNPs = allSNPs[, row.names(mm__)]
unsigDF = list()
set.seed(35)
for(i in 1:1000) {
  unsigDF[[i]] = allSNPs[sample(unsig, length(popo1), replace = F), ]
}
ddd = lapply(unsigDF, function(x) dist(t(x)))

Mres = lapply(ddd, function(x) mantel(x, pipi))
Mres__ = sapply(Mres, function(x) x$statistic)

# pepe = mget(ls())
saveRDS(Mres__, file = 'mantelNull.rds')
