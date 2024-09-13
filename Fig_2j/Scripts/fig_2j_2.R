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
