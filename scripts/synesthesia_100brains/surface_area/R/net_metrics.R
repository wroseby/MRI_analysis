# SCRIPT TO PLOT NETWORK METRICS #

library(ggplot2)
library(reshape2)

clust_wrd60_net_metrics = read.csv(paste0(savepath_data, '/hclust/clust_wrd60_net_metrics.csv'))
clust_wrd60_net_metrics = melt(clust_wrd60_net_metrics)

clust_wrd60_net_metrics = clust_wrd60_net_metrics %>%
  group_by(variable) %>%
  mutate(value_norm = normalise(value))

ggplot(data = clust_wrd60_net_metrics, aes(x = variable, y = log(value), fill = Group))+
  geom_col(position = position_dodge())+
  scale_fill_manual(values = c('blue', 'red'))+
  labs(x = 'Network metric', y = 'log(value)')+
  theme_minimal()
ggsave(filename = 'clust_wrd60_SA_abs_positions.png', path = paste0(savepath_outputs), # save as output
       width = 2200, height = 1080, units = 'px')
