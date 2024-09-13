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
