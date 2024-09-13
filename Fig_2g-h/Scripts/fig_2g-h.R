library(vegan)
library(ggplot2)

# GWAS SNPs =====================================
ssig = read.delim('./microbiomeSNPs.csv', sep = ',')
ssig = unique(ssig$significantSNPs)

# Other SNPs from Troels ========================
unsig = read.delim('./SNPs_to_generate_random_sets.txt', sep = '\t', header = F)
unsig = paste0('SNP_', unsig$V1)
unsig = setdiff(unsig, ssig)  # remove 78

# All SNP data ==================================
allSNPs = read.delim('./20220330_lotus_snps.csv', sep = ',')  # Unzip the data first.
row.names(allSNPs) = paste('SNP', allSNPs$CHR, allSNPs$POS, sep = '_')
allSNPs = allSNPs[, c(-1, -2)]
ssigDF = allSNPs[ssig, ]

# This part is from cluster.
snpDist = readRDS('./snpDist.rds')

# Distance matrices of GWAS SNPs
dp1 = dist(t(ssigDF))

# PCAs
pca1 = prcomp(dp1)
pca2 = prcomp(snpDist)

# Nothern-southern information ==================
pop = read.delim('./lotusPop.csv', sep = ';', row.names = 1)
row.names(pop)[1] = 'Gifu_'
sss = intersect(row.names(pop), row.names(pca1$x))
pop = pop[row.names(pca1$x), ,drop = FALSE]

# PCA plots with nothern-southern colors
qplot(pca1$x[,1], pca1$x[,2], geom = 'blank') + 
  geom_point(aes(color = pop$Population), size = 6) + 
  scale_color_manual(values = c('royalblue3', 'firebrick3')) +
  xlab('PC1') + ylab('PC2') + 
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = 'bottom',
        axis.ticks = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
ggsave('./fig_2h.pdf', height = 5, width = 6)

qplot(pca2$x[,1], pca2$x[,2], geom = 'blank') + 
  geom_point(aes(color = pop$Population), size = 6) + 
  scale_color_manual(values = c('royalblue3', 'firebrick3')) +
  xlab('PC1') + ylab('PC2') + 
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = 'bottom',
        axis.ticks = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
ggsave('./fig_2g.pdf', height = 5, width = 6)
