library(Spectrum)

Spectrum(pcor(SA_abs_lr[SA_abs_lr$Group == "Control", 8:ncol(SA_abs_lr)])$estimate,
         runrange = T, krangemax = 34)

Spectrum_out = Spectrum(pcor(SA_abs_lr[SA_abs_lr$Group == "Control", 8:ncol(SA_abs_lr)])$estimate,
                        runrange = T, krangemax = 34)[[30]]

clust_spec_SA_abs = data.frame(Area = colnames(Spectrum_out$similarity_matrix), cluster = Spectrum_out$assignments)

# Plot cluster anatomical positions
clust_spec_SA_abs_pos = merge(parcel_positions, clust_spec_SA_abs) # merge clusters and positions
clust_spec_SA_abs_pos$cluster = factor(clust_spec_SA_abs_pos$cluster) # factor clusters for plotting
ggplot(data = clust_spec_SA_abs_pos, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 0.5, shape = 21, stroke = 0.5)+
  geom_text(aes(label = gsub(x = cluster, pattern = 'cluster ', replacement = '')), size = 1.5, vjust = -0.5)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_spec_SA_abs_positions.png', path = savepath_outputs, # save as output
       width = 1500, height = 1080, units = 'px')

# Summarise within clusters
SA_abs_spec90_avg = summarise_clusters(SA_abs_lr, clust_spec_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_spec90_avg[[1]][SA_abs_spec90_avg[[1]]$Group == 'Syn', 8:ncol(SA_abs_spec90_avg[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_spec90_avg[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec90_avg[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_avg.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec90_avg[[3]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls to match syn sample size and repeat cluster averaging, multiple times
set.seed(6475) # set seed with rn
for (i in seq(1:10)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),] # reform data with n = 102
  SA_abs_spec90_avg_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_spec_SA_abs_pos)
  write.csv(SA_abs_spec90_avg_102[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_spec90_avg[[1]][SA_abs_spec90_avg[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum_syn.csv'), row.names = F) # write just the syn