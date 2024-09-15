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

# Description: This script plots the histogram of the Mantel
# test results from the null distribution and the observed
# mantel statistic from GWAS snps.
# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

########## Generate the histogram from calculated values ##########
library(ggplot2)

results1 = readRDS('mantelNull.rds')
summary(results1)

qplot(x = results1, geom = 'blank') + 
  scale_x_continuous(limits = c(0, max(results1) + 0.1)) +
  geom_histogram(fill = 'mediumvioletred', bins = 100) +
  xlab("Mantel statistic") + 
  ylab('Frequency') +
  geom_vline(xintercept = mean(c(0.1015, 0.1016, 0.09884)), 
             color = 'black', linetype = 'dashed', size = 1.25) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = 'bottom',
        axis.ticks = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
ggsave('./fig_2j.pdf', height = 5, width = 6)
