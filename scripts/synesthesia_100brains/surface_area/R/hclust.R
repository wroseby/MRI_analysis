# SCRIPT TO PERORM HIERARCHICAL CLUSTERING OF PARCELS AND PLOT RESULTS #

# Load packages
library(ppcor)
library(dendextend)
library(dplyr)
library(viridis)
library(ggplot2)
library(reshape2)
library(mclust)
library(cluster)
library(clValid)


set.seed(6780) # for consistency

# Function to perform hierarchical cluster and plot dendrogram
hclust.plot = function(data, group, clust_method, cut = 2, clustn){
 
  # generate correlation matrix - partial pearson
  cormat = pcor(x = data[data$Group == group, 8:ncol(data)],
           method = 'pearson')$estimate
      
  distance = dist(cormat, method = 'euclidean') # euclidean distance: sum of squared differences; positive and negative correlations can be considered separately
  clust = hclust(distance, method = clust_method) # single method to do nearest neighbour (paired clusters); Ward's method to minimise within-cluster variance
  
  # plot dendrogram
  dend = as.dendrogram(clust) # convert to dendrogram
  dend %>% 
    set('branches_k_color', viridis(clustn), k = clustn) %>% # colour branches by cluster
    set('branches_lwd', 1.8) %>% # set linewidth
    set('labels_cex', 0.8) %>% # set label size
    plot(main = paste0(clustn,' clusters of 180 parcels')) # plot with title
  
  # reorder original cormat to match clustering order and re-generate clustering
  cormat = cormat[clust$order, clust$order] # reorder
  distance = dist(cormat, method = 'euclidean') # euclidean distance: sum of squared differences; positive and negative correlations can be considered separately
  clust = hclust(distance, method = clust_method) # single method to do nearest neighbour (paired clusters); Ward's method to minimise within-cluster variance
  
  # acquire cluster identities
  if (clust_method == 'single'){
    clusters = stats::cutree(clust, k = clustn) # desired number of clusters if nearest neighbour
  }
  if (clust_method == 'ward.D'){
    clusters = stats::cutree(clust, h = cut) # level to cut at if Ward's method
  }
  
  clusters = data.frame(Area = names(clusters), cluster = clusters) # create dataframe
  
  return(list(clusters, distance))
} 

# Hierarchical clustering - Nearest neighbour
# Perform nearest neighbour (aka single linkage) clustering on full control data and retrieve pairwise clusters
png(paste0(savepath_outputs,"clust_nn90_dend.png"), width = 1600, height = 720) # dendrogram plot .png
clust_nn_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', clust_method = 'single', clustn = 90)
dev.off() # save .png

# Plot cluster anatomical positions
clust_nn_SA_abs_pos = merge(parcel_positions, clust_nn_SA_abs[[1]]) # merge clusters and positions
clust_nn_SA_abs_pos$cluster = factor(clust_nn_SA_abs_pos$cluster) # factor clusters for plotting
ggplot(data = clust_nn_SA_abs_pos, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 1, shape = 21, stroke = 0.5)+
  geom_text(aes(label = gsub(x = cluster, pattern = 'cluster ', replacement = '')), size = 1.5, vjust = -0.8)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_nn_SA_abs_positions.png', path = paste0(savepath_outputs), # save as output
       width = 1500, height = 1080, units = 'px')

# Diagnostics - nearest neighbour
# Silhouette score (internal validity)
sil_scores = silhouette(x = as.integer(clust_nn_SA_abs[[1]]$cluster), dist = clust_nn_SA_abs[[2]]) # silhouette score using cluster assingments and similarity matrix
mean(sil_scores[,"sil_width"]) # mean silhouette width = 0.031

# Repeat separately for left and right hemispheres (consistency)
SA_abs_l = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^L_')]) # select left parcels
clust_nn_SA_abs_l = hclust.plot(SA_abs_l, group = 'Control',  clust_method = 'single', clustn = 90) # cluster left parcels
SA_abs_r = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^R_')]) # select right parcels
clust_nn_SA_abs_r = hclust.plot(SA_abs_r, group = 'Control', clust_method = 'single', clustn = 90)  # cluster right parcels
adjustedRandIndex(clust_nn_SA_abs_l[[1]]$cluster, clust_nn_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.426

# Repeat for 10 resamples of control dataset (consistency)
clust_nn_resamp = list() # initialise list of resample clusters
set.seed(6780)# set seed with rn
for (i in seq(1:10)){ # for 10 resamples
  SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
  clust_nn_resamp[[i]] = hclust.plot(SA_abs_lr_resamp, group = 'Control', clust_method = 'single', clustn = 90) # nn clustering
}
clust_nn_resamp_ari = vector() # initialise vector of adjusted rand indexes
for (i in 1:ncol(combn(10,2))){ # for every pair of clusterings
  clust_nn_resamp_ari[i] = adjustedRandIndex(clust_nn_resamp[[combn(10,2)[,i][1]]][[1]]$cluster, # get the adjusted rand index
                                               clust_nn_resamp[[combn(10,2)[,i][2]]][[1]]$cluster)
}
mean(clust_nn_resamp_ari) # mean ari of resample comparison = 0.059

# Summarise data within 90 nearest neighbour clusters for both controls and synesthetes
summarise_clusters = function(data, clusters){ # function to average across clusters
  data_long = melt(data[,c(1,8:ncol(data))]) # convert parcel data to long form for merging
  colnames(data_long) = c('ID', 'Area', 'value') # rename columns for merging
  data_long = merge(data_long, clusters, by = 'Area') # merge data with cluster assignments
  
  # summarise within clusters 
  data_summarised = data_long %>% 
    group_by(ID, cluster) %>%
    summarise(sum_value = sum(value), # sum the SA values
              mean_value = mean(value), # average the SA values
              x = mean(x), # average the x coordinate
              y = mean(y),
              Region = first(Region)) # average the y coordinate
  
  # convert back to wide form and bind with metadata
  data_wide_sum = dcast(data_summarised, ID ~ cluster, value.var = 'sum_value')
  data_wide_sum = merge(data[,1:7], data_wide_sum, by = 'ID')
  data_wide_mean = dcast(data_summarised, ID ~ cluster, value.var = 'mean_value')
  data_wide_mean = merge(data[,1:7], data_wide_mean, by = 'ID')
  
  # get centroids of x and y positions for clusters
  data_pos = unique(dplyr::select(ungroup(data_summarised), 2,5,6,7)) # select cluster column, mean x and mean y columns
  
  return(list(data_wide_sum, data_wide_mean, data_pos))
}

# Average the clusters and export results
SA_abs_nn90_sum = summarise_clusters(SA_abs_lr, clust_nn_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_nn90_sum[[1]][SA_abs_nn90_sum[[1]]$Group == 'Syn', 8:ncol(SA_abs_nn90_sum[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_nn90_sum[[1]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum.csv'), row.names = F) # write summarised clusters to .csv
write.csv(SA_abs_nn90_sum[[3]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls to match syn sample size and repeat cluster averaging
set.seed(6780)# set seed with rn
SA_abs_lr_resamp = rbind(SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),], # reform data with n = 102
                         SA_abs_lr[SA_abs_lr$Group == 'Syn',])
SA_abs_nn90_sum_102 = summarise_clusters(SA_abs_lr_resamp, clust_nn_SA_abs_pos)
write.csv(SA_abs_nn90_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_resamp.csv'), row.names = F) # write to .csv

# Resample the controls 6 times without replacement for unique groups
set.seed(6780)
sample_numbers = sample(1:650, size = 6 * 102, replace = F) # generate random samples
sample_numbers = split(sample_numbers, ceiling(seq_along(sample_numbers)/102)) # split into 6 groups
for (i in seq(1:6)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample_numbers[[i]],]
  SA_abs_nn90_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_nn_SA_abs_pos)
  write.csv(SA_abs_nn90_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_nn90_sum[[1]][SA_abs_nn90_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_syn.csv'), row.names = F) # write just the syn

# Hierarhical cluster - Ward's method
# Perform clustering using Ward's method
png(paste0(savepath_outputs,"clust_wrd60_dend.png"), width = 1600, height = 720) # dendrogram plot .png
clust_wrd_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', clust_method = 'ward.D', cut = 2, clustn = 60)
dev.off() # save png

# Plot cluster anatomical positions
clust_wrd_SA_abs_pos = merge(parcel_positions, clust_wrd_SA_abs[[1]]) # merge clusters and positions
clust_wrd_SA_abs_pos$cluster = factor(clust_wrd_SA_abs_pos$cluster) # factor clusters for plotting
ggplot(data = clust_wrd_SA_abs_pos, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 0.6, shape = 21, stroke = 0.5)+
  geom_text(aes(label = gsub(x = cluster, pattern = 'cluster ', replacement = '')), size = 1.5, vjust = -0.5)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_wrd_SA_abs_positions.png', path = paste0(savepath_outputs), # save as output
       width = 1500, height = 1080, units = 'px')

# Diagnostics - nearest neighbour
# Silhouette score (internal validity)
sil_scores = silhouette(x = as.integer(clust_wrd_SA_abs[[1]]$cluster), dist = clust_wrd_SA_abs[[2]]) # silhouette score using cluster assingments and similarity matrix
mean(sil_scores[,"sil_width"]) # mean silhouette width = 0.070

# Repeat separately for left and right hemispheres (consistency)
SA_abs_l = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^L_')]) # select left parcels
clust_wrd_SA_abs_l = hclust.plot(SA_abs_l, group = 'Control',  clust_method = 'ward.D', clustn = 60, cut = 2) # cluster left parcels
SA_abs_r = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^R_')]) # select right parcels
clust_wrd_SA_abs_r = hclust.plot(SA_abs_r, group = 'Control', clust_method = 'ward.D', clustn = 60, cut = 2)  # cluster right parcels
adjustedRandIndex(clust_wrd_SA_abs_l[[1]][order(clust_wrd_SA_abs_l[[1]]$Area),]$cluster, clust_wrd_SA_abs_r[[1]][order(clust_wrd_SA_abs_l[[1]]$Area),]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.479

# Repeat for resamples of control dataset (consistency)
clust_wrd_resamp = list() # initialise list of resample clusters
set.seed(6780)# set seed with rn
for (i in seq(1:10)){ # for 10 resamples
  SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
  clust_wrd_resamp[[i]] = hclust.plot(SA_abs_lr_resamp, group = 'Control', clust_method = 'ward.D', clustn = 60, cut = 2) # wrd clustering
}
clust_wrd_resamp_ari = vector() # initialise vector of adjusted rand indexes
for (i in 1:ncol(combn(10,2))){ # for every pair of clusterings
  resamp1 = clust_wrd_resamp[[combn(10,2)[,i][1]]][[1]]
  resamp2 = clust_wrd_resamp[[combn(10,2)[,i][2]]][[1]]
  clust_wrd_resamp_ari[i] = adjustedRandIndex(resamp1[order(resamp1$Area),]$cluster,
                                              resamp2[order(resamp2$Area),]$cluster)
}
mean(clust_wrd_resamp_ari) # mean ari of resample comparison = 0.024

# Summarise data within 60 ward clusters for both controls and synesthetes
# Average the clusters and export results
SA_abs_wrd60_sum = summarise_clusters(SA_abs_lr, clust_wrd_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_wrd60_sum[[1]][SA_abs_wrd60_sum[[1]]$Group == 'Syn', 8:ncol(SA_abs_wrd60_sum[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_wrd60_sum[[1]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum.csv'), row.names = F) # write summarised clusters to .csv
write.csv(SA_abs_wrd60_sum[[3]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls to match syn sample size and repeat cluster averaging
set.seed(6780)# set seed with rn
SA_abs_lr_resamp = rbind(SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),], # reform data with n = 102
                         SA_abs_lr[SA_abs_lr$Group == 'Syn',])
SA_abs_wrd60_sum_102 = summarise_clusters(SA_abs_lr_resamp, clust_wrd_SA_abs_pos)
write.csv(SA_abs_wrd60_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_resamp.csv'), row.names = F) # write to .csv

# Resample the controls without replacement for 6 unique groups
set.seed(6780)
sample_numbers = sample(1:650, size = 6 * 102, replace = F) # generate random samples
sample_numbers = split(sample_numbers, ceiling(seq_along(sample_numbers)/102)) # split into 6 groups
for (i in seq(1:6)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample_numbers[[i]],]
  SA_abs_wrd60_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_wrd_SA_abs_pos)
  write.csv(SA_abs_wrd60_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_wrd60_sum[[1]][SA_abs_wrd60_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_syn.csv'), row.names = F) # write just the syn

# Plot the correlation distributions
SA_abs_wrd60_resamps = list()
for (i in 1:6){ # re-read the written csvs for consistency
  SA_abs_wrd60_resamps[[i]] = read.csv(paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_cntrl_resamp_', i, '.csv'))
}
SA_abs_wrd60_resamps[[7]] = read.csv(paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_syn.csv')) # read syn csv
SA_abs_wrd60_resamps = lapply(SA_abs_wrd60_resamps, FUN = function(df){df %>% mutate(ID = as.character(ID))}) # convert all IDs to character
SA_abs_wrd60_resamps_pcors = lapply(SA_abs_wrd60_resamps, FUN = function(df){pcor(df[,8:67])$estimate}) # generate pcor matrix
SA_abs_wrd60_resamps_pcors = lapply(SA_abs_wrd60_resamps_pcors, FUN = function(df){df = melt(df)}) # convert matrix to long form
SA_abs_wrd60_resamps_pcors = bind_rows(SA_abs_wrd60_resamps_pcors, .id = 'Sample') # bind data into one by resample
SA_abs_wrd60_resamps_pcors$Group = ifelse(SA_abs_wrd60_resamps_pcors$Sample == 7, 'Syn', 'Control') # 

ggplot(data = SA_abs_wrd60_resamps_pcors, aes(x = value, color = Group, group = Sample))+ # plot parcel locations with fill colour for cluster
  geom_density(lwd = 1.05)+
  scale_color_manual(values = c('blue', 'red'))+
  labs(x = 'Partial correlation', y = 'Density')
ggsave(filename = 'clust_wrd_SA_abs_pcor_distributions.png', path = paste0(savepath_outputs), # save as output
       width = 1800, height = 1080, units = 'px')




# Iterate over cut values to determine best clusterings - Ward's method
mean_silhouettes = vector()
dunn_indices = vector()
lr_aris = vector()
resample_mean_aris = vector()

cut_values = seq(from = 1.40, to = 2.50, by = 0.01)

for (i in seq_along(cut_values)){
  clust_wrd_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', clust_method = 'ward.D', cut = cut_values[i], clustn = 10)
  
  sil_scores = silhouette(x = as.integer(clust_wrd_SA_abs[[1]]$cluster), dist = clust_wrd_SA_abs[[2]])
  mean_silhouettes[i] = mean(sil_scores[,"sil_width"]) 
  dunn_indices[i] = dunn(distance = clust_wrd_SA_abs[[2]], clusters = as.integer(clust_wrd_SA_abs[[1]]$cluster))
  
  set.seed(6780)
  clust_wrd_SA_abs_l = hclust.plot(SA_abs_l, group = 'Control',  clust_method = 'ward.D', clustn = 10, cut = cut_values[i]) # cluster left parcels
  clust_wrd_SA_abs_r = hclust.plot(SA_abs_r, group = 'Control', clust_method = 'ward.D', clustn = 10, cut = cut_values[i])  # cluster right parcels
  lr_aris[i] = adjustedRandIndex(clust_wrd_SA_abs_l[[1]]$cluster, clust_wrd_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.495
  
  set.seed(6780)
  for (j in seq(1:10)){ # for 10 resamples
    SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
    clust_wrd_resamp[[j]] = hclust.plot(SA_abs_lr_resamp, group = 'Control', clust_method = 'ward.D', clustn = 10, cut = cut_values[i]) # wrd clustering
  }
  clust_wrd_resamp_ari = vector() # initialise vector of adjusted rand indexes
  for (j in 1:ncol(combn(10,2))){ # for every pair of clusterings
    resamp1 = clust_wrd_resamp[[combn(10,2)[,j][1]]][[1]]
    resamp2 = clust_wrd_resamp[[combn(10,2)[,j][2]]][[1]]
    clust_wrd_resamp_ari[j] = adjustedRandIndex(resamp1[order(resamp1$Area),]$cluster,
                                                resamp2[order(resamp2$Area),]$cluster)
  }
  resample_mean_aris[i] = mean(clust_wrd_resamp_ari) # mean ari of resample comparison = 0.024
}

cut_values[which(mean_silhouettes == max(mean_silhouettes))] # best silhouette cuts = 1.81, 1.82
cut_values[which(dunn_indices == max(dunn_indices))] # best silhouette cuts = 1.42, 1.43; decreases with increasing cut (not useful)
cut_values[which(lr_aris == max(lr_aris))] # best lr consistency = 2.41
cut_values[which(resample_mean_aris == max(na.omit(resample_mean_aris)))] # 2.06
# find the maximum of all; normalise then aggregate
normalise = function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
mean_silhouettes_n = normalise(mean_silhouettes)
lr_aris_n = normalise(lr_aris)
resample_mean_aris_n = normalise(resample_mean_aris)
agg_score = (mean_silhouettes_n + lr_aris_n + resample_mean_aris_n)
cut_values[which(agg_score == max(agg_score, na.rm = T))] # 2.03

