# SCRIPT TO PERORM HIERARCHICAL CLUSTERING OF PARCELS AND PLOT RESULTS #

# Load packages
library(ppcor)
library(dplyr)
library(ggplot2)
library(reshape2)

set.seed(154) # for consistency

# Function to perform hierarchical cluster and plot dendrogram
hclust.plot <- function(data, group, clust_method, cut){
 
  # generate correlation matrix - partial pearson
  cormat = pcor(x = data[data$Group == group, 8:ncol(data)],
           method = 'pearson')$estimate
      
  distance = dist(cormat, method = 'euclidean') # euclidean distance: sum of squared differences; positive and negative correlations can be considered separately
  clust = hclust(distance, method = clust_method) # single method to do nearest neighbour (paired clusters); Ward's method to minimise within-cluster variance
  
  plot(clust, main = paste('Clustering:', group), sub = paste0("Euclidean distance, method = ", clust_method), # plot dendrogram
       cex = 1.3, cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.5, lwd = 2.0)  # adjust overall text size
  
  # acquire cluster identities
  if (clust_method == 'single'){
    clusters = cutree(clust, k = cut) # desired number of clusters if nearest neighbour
  }
  if (clust_method == 'ward.D'){
    clusters = cutree(clust, h = cut) # level to cut at if Ward's method
  }
  
  clusters = data.frame(Area = names(clusters), cluster = clusters)
  return(list(clusters, distance))
} 

# Hierarchical clustering - Nearest neighbour
# Perform nearest neighbour (aka single linkage) clustering on full control data and retrieve pairwise clusters
clust_nn_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', clust_method = 'single', cut = 90)

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
ggsave(filename = 'clust_nn_SA_abs_positions.png', path = savepath_outputs, # save as output
       width = 1500, height = 1080, units = 'px')

# Diagnostics - nearest neighbour
# Silhouette score (internal validity)
sil_scores = silhouette(x = as.integer(clust_nn_SA_abs[[1]]$cluster), dist = clust_nn_SA_abs[[2]]) # silhouette score using cluster assingments and similarity matrix
mean(sil_scores[,"sil_width"]) # mean silhouette width = 0.031

# Repeat separately for left and right hemispheres (consistency)
SA_abs_l = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^L_')]) # select left parcels
clust_nn_SA_abs_l = hclust.plot(SA_abs_l, group = 'Control',  clust_method = 'single', cut = 90) # cluster left parcels
SA_abs_r = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^R_')]) # select right parcels
clust_nn_SA_abs_r = hclust.plot(SA_abs_r, group = 'Control', clust_method = 'single', cut = 90)  # cluster right parcels
adjustedRandIndex(clust_nn_SA_abs_l[[1]]$cluster, clust_nn_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.426

# Repeat for resamples of control dataset (consistency)
clust_nn_resamp = list() # initialise list of resample clusters
set.seed(6475) # set seed with rn
for (i in seq(1:10)){ # for 10 resamples
  SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
  clust_nn_resamp[[i]] = hclust.plot(SA_abs_lr_resamp, group = 'Control', clust_method = 'single', cut = 90) # nn clustering
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
              mean_x = mean(x), # average the x coordinate
              mean_y = mean(y)) # average the y coordinate
  
  # convert back to wide form and bind with metadata
  data_wide_sum = dcast(data_summarised, ID ~ cluster, value.var = 'sum_value')
  data_wide_sum = merge(data[,1:7], data_wide_sum, by = 'ID')
  data_wide_mean = dcast(data_summarised, ID ~ cluster, value.var = 'mean_value')
  data_wide_mean = merge(data[,1:7], data_wide_mean, by = 'ID')
  
  # get centroids of x and y positions for clusters
  data_pos = unique(dplyr::select(ungroup(data_summarised), 2,4,5)) # select cluster column, mean x and mean y columns
  
  return(list(data_wide_sum, data_wide_mean, data_pos))
}

# Average the clusters and export results
SA_abs_nn90_sum = summarise_clusters(SA_abs_lr, clust_nn_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_nn90_sum[[1]][SA_abs_nn90_sum[[1]]$Group == 'Syn', 8:ncol(SA_abs_nn90_sum[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_nn90_sum[[1]], file = paste0(savepath_data, 'SA_abs_nn90_sum.csv'), row.names = F) # write summarised clusters to .csv
write.csv(SA_abs_nn90_sum[[3]], file = paste0(savepath_data, 'SA_abs_nn90_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls to match syn sample size and repeat cluster averaging
set.seed(6475) # set seed with rn
SA_abs_lr_resamp = rbind(SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),], # reform data with n = 102
                         SA_abs_lr[SA_abs_lr$Group == 'Syn',])
SA_abs_nn90_sum_102 = summarise_clusters(SA_abs_lr_resamp, clust_nn_SA_abs_pos)
write.csv(SA_abs_nn90_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_resamp.csv'), row.names = F) # write to .csv

# Resample the controls to match syn sample size and repeat cluster averaging, multiple times
set.seed(6475) # set seed with rn
for (i in seq(1:10)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),] # reform data with n = 102
  SA_abs_nn90_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_nn_SA_abs_pos)
  write.csv(SA_abs_nn90_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_nn90_sum[[1]][SA_abs_nn90_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'hclust/SA_abs_nn90_sum_syn.csv'), row.names = F) # write just the syn




# Hierarhical cluster - Ward's method
# Perform clustering using Ward's method
clust_wrd_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', clust_method = 'ward.D', cut = 2)

# Plot cluster anatomical positions
clust_wrd_SA_abs_pos = merge(parcel_positions, clust_wrd_SA_abs[[1]]) # merge clusters and positions
clust_wrd_SA_abs_pos$cluster = factor(clust_wrd_SA_abs_pos$cluster) # factor clusters for plotting
ggplot(data = clust_wrd_SA_abs_pos, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 1, shape = 21, stroke = 0.5)+
  geom_text(aes(label = gsub(x = cluster, pattern = 'cluster ', replacement = '')), size = 1.5, vjust = -0.8)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_wrd_SA_abs_positions.png', path = savepath_outputs, # save as output
       width = 1500, height = 1080, units = 'px')

# Diagnostics - nearest neighbour
# Silhouette score (internal validity)
sil_scores = silhouette(x = as.integer(clust_wrd_SA_abs[[1]]$cluster), dist = clust_wrd_SA_abs[[2]]) # silhouette score using cluster assingments and similarity matrix
mean(sil_scores[,"sil_width"]) # mean silhouette width = 0.070

# Repeat separately for left and right hemispheres (consistency)
SA_abs_l = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^L_')]) # select left parcels
clust_wrd_SA_abs_l = hclust.plot(SA_abs_l, group = 'Control',  clust_method = 'ward.D', cut = 2) # cluster left parcels
SA_abs_r = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^R_')]) # select right parcels
clust_wrd_SA_abs_r = hclust.plot(SA_abs_r, group = 'Control', clust_method = 'ward.D', cut = 2)  # cluster right parcels
adjustedRandIndex(clust_wrd_SA_abs_l[[1]]$cluster, clust_wrd_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.495

# Repeat for resamples of control dataset (consistency)
clust_wrd_resamp = list() # initialise list of resample clusters
set.seed(6475) # set seed with rn
for (i in seq(1:10)){ # for 10 resamples
  SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
  clust_wrd_resamp[[i]] = hclust.plot(SA_abs_lr_resamp, group = 'Control', clust_method = 'ward.D', cut = 2) # wrd clustering
}
clust_wrd_resamp_ari = vector() # initialise vector of adjusted rand indexes
for (i in 1:ncol(combn(10,2))){ # for every pair of clusterings
  clust_wrd_resamp_ari[i] = adjustedRandIndex(clust_wrd_resamp[[combn(10,2)[,i][1]]][[1]]$cluster, # get the adjusted rand index
                                             clust_wrd_resamp[[combn(10,2)[,i][2]]][[1]]$cluster)
}
mean(clust_wrd_resamp_ari) # mean ari of resample comparison = 0.024

# Summarise data within 60 ward clusters for both controls and synesthetes
# Average the clusters and export results
SA_abs_wrd60_sum = summarise_clusters(SA_abs_lr, clust_wrd_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_wrd60_sum[[1]][SA_abs_wrd60_sum[[1]]$Group == 'Syn', 8:ncol(SA_abs_wrd60_sum[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_wrd60_sum[[1]], file = paste0(savepath_data, 'SA_abs_wrd60_sum.csv'), row.names = F) # write summarised clusters to .csv
write.csv(SA_abs_wrd60_sum[[3]], file = paste0(savepath_data, 'SA_abs_wrd60_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls to match syn sample size and repeat cluster averaging
set.seed(6475) # set seed with rn
SA_abs_lr_resamp = rbind(SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),], # reform data with n = 102
                         SA_abs_lr[SA_abs_lr$Group == 'Syn',])
SA_abs_wrd60_sum_102 = summarise_clusters(SA_abs_lr_resamp, clust_wrd_SA_abs_pos)
write.csv(SA_abs_wrd60_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_resamp.csv'), row.names = F) # write to .csv

# Resample the controls to match syn sample size and repeat cluster averaging, multiple times
set.seed(6475) # set seed with rn
for (i in seq(1:10)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(seq(1:650), size = 102, replace = F),] # reform data with n = 102
  SA_abs_wrd60_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_wrd_SA_abs_pos)
  write.csv(SA_abs_wrd60_sum_102[[1]], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_wrd60_sum[[1]][SA_abs_wrd60_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'hclust/SA_abs_wrd60_sum_syn.csv'), row.names = F) # write just the syn