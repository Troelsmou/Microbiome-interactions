library(magrittr)
library(stringr)
library(vegan)
library(pheatmap)
library(ggplot2)

# Colors ========================================
colors_ = list(black = 'gray10', 
               white = 'antiquewhite1', 
               gray = 'darkgrey', 
               pink = 'deeppink', 
               green = 'springgreen3', 
               yellow = 'gold', 
               blue = 'dodgerblue3', 
               turk = 'darkturquoise', 
               coffee = 'sienna4', 
               red = 'firebrick2', 
               purple = 'purple4', 
               orange = 'darkorange1', 
               lightred = 'coral1', 
               lightblue = 'lightblue1', 
               lightpurple = 'lavenderblush2', 
               lightgreen = 'palegreen2',
               bluegreen = 'aquamarine1',
               anotherpurple = 'darkorchid1',
               yetanotherpurple = 'mediumvioletred',
               laci = 'midnightblue',
               anotherred = 'palevioletred1')

# Data preparation ==============================
# Count table
asv = read.delim('./isolate210722HQ.csv', sep = ',', header = TRUE, row.names = 1)
asv = asv[, 1:495]
accessions = strsplit(colnames(asv), split = '_') %>% sapply(function(x) x[[1]])
numbers = gsub(pattern = '[A-Z]', replacement = '', x = accessions, ignore.case = TRUE)
numbers = str_pad(string = numbers, width = 3, side = 'left', pad = '0')[4:495]
accessions = c('Gifu_', 'Gifu_', 'Gifu_', paste0('MG', numbers))
asv = asv[grepl(pattern = '^Lj', x = row.names(asv)), ]
asv = asv[rowSums(asv) > 0, ]
lj = sapply(strsplit(row.names(asv), split = '='), function(x) x[[1]])
row.names(asv) = lj
colnames(asv) = paste(accessions, 1:3, sep = '_')

# Sample table
st = data.frame(genotype = accessions)
row.names(st) = colnames(asv)
st$Batch = strsplit(row.names(st), split = '_|__') %>% sapply(function(x) x[[2]]) %>% as.character

# Taxonomy
taxa = read.delim('./taxonomy.csv', sep = ';', row.names = 1)
taxa = taxa[lj, ]

# Relative abundance
raav = vegan::decostand(asv, method = 'total', MARGIN = 2) * 100000
raav = round(raav)
raav2 = log2(raav + 1)
raav3 = raav2[order(apply(raav, 1, mean), decreasing = TRUE),]

# Heatmap colors
colForHM = colorRampPalette(c(colors_$gray, colors_$lightblue, colors_$laci))(n = 100)
colForAnn = list(Batch = c('1' = colors_$coffee,
                           '2' = colors_$red,
                           '3' = colors_$orange),
                 Cl    = c('1' = colors_$white,
                           '2' = colors_$green,
                           '3' = colors_$yellow,
                           '4' = colors_$black)
)

pheatmap(raav3, 
         border_color = NA, scale = 'none', 
         cluster_rows = F, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE, 
         clustering_distance_cols = vegdist(t(raav), 'cao'), clustering_method = 'ward.D2',
         clustering_distance_rows = vegdist(raav, 'cao'),
         annotation_colors = colForAnn,
         annotation_col = st[, 'Batch', drop = FALSE],
         color = colForHM, 
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         filename = './fig_1c.pdf', width = 8, height = 7,
         annotation_legend = TRUE, legend = TRUE
         )
