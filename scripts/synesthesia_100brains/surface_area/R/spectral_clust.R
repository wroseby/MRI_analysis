# SCRIPT TO PERFORM SPECTRAL CLUSTERING OF PARCELS #

# Load packages
library(cluster)
library(Matrix)
library(ppcor)
library(MASS)
library(ggplot2)
library(reshape2)
library(dplyr)
library(mclust)

# Function to calculate eta squared
eta_squared = function(x, y) {
  1 - sum((x - y)^2) / (sum((x - mean(x))^2) + sum((y - mean(y))^2))
}

# Function to calculate gaussian kernel distance
gaussian_sim = function(x){as.matrix(exp(-dist(x)^2 / (2 * 1^2)))}

# Function to perform spectral clustering
spectral_clust = function(data, group, n_clusters, normalise = F){
  
  # Get partial correlation in desired group
  cormat = pcor(x = data[data$Group == group, 8:ncol(data)],
                method = 'pearson')$estimate
  
  # create similarity matrix using eta-squared from the input matrix
  n = ncol(cormat)
  S = matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      S[i, j] <- eta_squared(cormat[, i], cormat[, j])
    }
  }
  diag(S) <- 1 # set diagonal to 1
  D = diag(rowSums(S)) # degree matrix
  L = D - S # Laplacian
  
  if (normalise == T){ 
    # symmetrically normalised laplacian
    #L = sqrt(ginv(D)) %*% L %*% sqrt(ginv(D)) # produces NaNs due to small weights of some nodes
    # random walk laplacian
    L = ginv(D) %*% L
  }
  
  eigen_decomp = eigen(L) # eigendecomposition of the Laplacian

  U = eigen_decomp$vectors[, (ncol(eigen_decomp$vectors) - n_clusters + 1):ncol(eigen_decomp$vectors)] # get smallest eigenvectors
  row_norms = sqrt(rowSums(U^2))
  U = sweep(U, 1, row_norms, "/") # normalise U
  clusters = kmeans(U, centers = n_clusters, nstart = 10)$cluster # k-means clustering
  
  # reorder the clusters to match data order
  unique_clusters = unique(clusters) # get unique clusters
  cluster_order = match(clusters, unique_clusters) # map the original clusters to unique clusters
  ordered_clusters = rep(NA, length(clusters)) # initialise ordered clusters
  for (i in seq_along(unique_clusters)) {
    ordered_clusters[clusters == unique_clusters[i]] <- i # reset order
  }
  
  spectral_clusters = data.frame(Area = colnames(cormat), cluster = ordered_clusters) # build the output data
  
  plot(U[, 1], U[, 2], col = clusters, pch = 19,
       xlab = "Eigenvector 1", ylab = "Eigenvector 2", main = "Spectral Clustering") # plot the clusters on the first two eigenvectors
  
  return(list(spectral_clusters, S)) # return the clusters and similarity matrix
}

# Spectral clustering of partial correlation matrix (unormalised laplacian)
set.seed(6780) # set seed with rn
clust_spec_SA_abs = spectral_clust(SA_abs_lr, group = 'Control', n_clusters = 90)

# Plot cluster anatomical positions
clust_spec_SA_abs_pos = merge(parcel_positions, clust_spec_SA_abs[[1]]) # merge clusters and positions
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

# Diagnostics
# Silhouette score (intenal validity)
sil_scores = silhouette(x = as.integer(clust_spec_SA_abs[[1]][,2]), dist = as.dist(1 - clust_spec_SA_abs[[2]])) # silhouette score using cluster assingments and similarity matrix
mean(sil_scores[,"sil_width"]) # mean silhouette width = 0.127

# Repeat separately for left and right hemispheres (consistency)
set.seed(6780) # set seed with rn
SA_abs_l = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^L_')]) # select left parcels
clust_spec_SA_abs_l = spectral_clust(SA_abs_l, group = 'Control', n_clusters = 90, normalise = F) # cluster left parcels
SA_abs_r = cbind(SA_abs[,1:7],SA_abs[,grep(x = colnames(SA_abs), pattern = '^R_')]) # select right parcels
clust_spec_SA_abs_r = spectral_clust(SA_abs_r, group = 'Control', n_clusters = 90, normalise = F)  # cluster pcors to 90
adjustedRandIndex(clust_spec_SA_abs_l[[1]]$cluster, clust_spec_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.485

# Repeat for resamples of control dataset (consistency)
clust_spec_resamp = list() # initialise list of resample clusters
set.seed(6780) # set seed with rn
for (i in seq(1:10)){
  SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
  clust_spec_resamp[[i]] = spectral_clust(SA_abs_lr_resamp, group = 'Control', n_clusters = 90, normalise = F) # spectral clustering
}
clust_spec_resamp_ari = vector() # initialise vector of adjusted rand indexes
for (i in 1:ncol(combn(10,2))){ # for every pair of clusterings
  clust_spec_resamp_ari[i] = adjustedRandIndex(clust_spec_resamp[[combn(10,2)[,i][1]]][[1]]$cluster, # get the adjusted rand index
                                               clust_spec_resamp[[combn(10,2)[,i][2]]][[1]]$cluster)
}
mean(clust_spec_resamp_ari) # mean ari of resample comparison = 0.267

# Summarise data within 90 spectral clusters for both controls and synesthetes
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

# Summarise the clusters and export results
SA_abs_spec90_avg = summarise_clusters(SA_abs_lr, clust_spec_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_spec90_avg[[1]][SA_abs_spec90_avg[[1]]$Group == 'Syn', 8:ncol(SA_abs_spec90_avg[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_spec90_avg[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec90_avg[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_avg.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec90_avg[[3]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls 6 times without replacement for unique groups
set.seed(6780)
sample_numbers = sample(1:650, size = 6 * 102, replace = F) # generate random samples
sample_numbers = split(sample_numbers, ceiling(seq_along(sample_numbers)/102)) # split into 6 groups
for (i in seq(1:6)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample_numbers[[i]],]
  SA_abs_spec90_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_spec_SA_abs_pos)
  write.csv(SA_abs_spec90_sum_102[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_spec90_sum[[1]][SA_abs_spec90_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec90_sum_syn.csv'), row.names = F) # write just the syn



# Iterate over cut values to determine best clusterings
mean_silhouettes = vector()
dunn_indices = vector()
lr_aris = vector()
resample_mean_aris = vector()

cluster_n = seq(from = 10, to = 100, by = 1)

for (i in seq_along(cluster_n)){
  clust_spec_SA_abs = spectral_clust(SA_abs_lr, group = 'Control', n_clusters = cluster_n[i])
  
  sil_scores = silhouette(x = as.integer(clust_spec_SA_abs[[1]][,2]), dist = as.dist(1 - clust_spec_SA_abs[[2]]))
  mean_silhouettes[i] = mean(sil_scores[,"sil_width"]) 
  dunn_indices[i] = dunn(distance = as.dist(1 - clust_spec_SA_abs[[2]]), clusters = as.integer(clust_spec_SA_abs[[1]][,2]))
  
  set.seed(6780)
  clust_spec_SA_abs_l = spectral_clust(SA_abs_l, group = 'Control', n_clusters = cluster_n[i]) # cluster left parcels
  clust_spec_SA_abs_r = spectral_clust(SA_abs_r, group = 'Control', n_clusters = cluster_n[i])  # cluster right parcels
  lr_aris[i] = adjustedRandIndex(clust_spec_SA_abs_l[[1]]$cluster, clust_spec_SA_abs_r[[1]]$cluster) # adjusted Rand index for l vs r cluster similarity = 0.495
  
  clust_spec_resamp = list()
  set.seed(6780)
  for (j in seq(1:10)){ # for 10 resamples
    SA_abs_lr_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample(1:650, size = 300, replace = F),] # take samples of 300 controls
    clust_spec_resamp[[j]] = spectral_clust(SA_abs_lr_resamp, group = 'Control', n_clusters = cluster_n[i])
  }
  clust_spec_resamp_ari = vector() # initialise vector of adjusted rand indexes
  for (j in 1:ncol(combn(10,2))){ # for every pair of clusterings
    resamp1 = clust_spec_resamp[[combn(10,2)[,j][1]]][[1]]
    resamp2 = clust_spec_resamp[[combn(10,2)[,j][2]]][[1]]
    clust_spec_resamp_ari[j] = adjustedRandIndex(resamp1[order(resamp1$Area),]$cluster,
                                                resamp2[order(resamp2$Area),]$cluster)
  }
  resample_mean_aris[i] = mean(clust_spec_resamp_ari) # mean ari of resample comparison = 0.024
}

cluster_n[which(mean_silhouettes == max(mean_silhouettes))] # best silhouette k = 69
cluster_n[which(dunn_indices == max(dunn_indices))] # best dunn k = 89
cluster_n[which(lr_aris == max(lr_aris))] # best lr consistency = 72
cluster_n[which(resample_mean_aris == max(na.omit(resample_mean_aris)))] # best resample k = 93
# find the maximum of all; normalise then aggregate
normalise = function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
mean_silhouettes_n = normalise(mean_silhouettes)
dunn_indices_n = normalise(dunn_indices)
lr_aris_n = normalise(lr_aris)
resample_mean_aris_n = normalise(resample_mean_aris)
agg_score = (mean_silhouettes_n + dunn_indices_n + lr_aris_n + resample_mean_aris_n)
cluster_n[which(agg_score == max(agg_score, na.rm = T))] # 79, 72 not including Dunn, also 72 not including Dunn or silhouette



# Redo clustering for 79 clusters
# Spectral clustering of partial correlation matrix (unormalised laplacian)
set.seed(7842) # set seed with rn
clust_spec_SA_abs = spectral_clust(SA_abs_lr, group = 'Control', n_clusters = 79)

# Plot cluster anatomical positions
clust_spec_SA_abs_pos = merge(parcel_positions, clust_spec_SA_abs[[1]]) # merge clusters and positions
clust_spec_SA_abs_pos$cluster = factor(clust_spec_SA_abs_pos$cluster) # factor clusters for plotting
ggplot(data = clust_spec_SA_abs_pos, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 1, shape = 21, stroke = 0.5)+
  geom_text(aes(label = gsub(x = cluster, pattern = 'cluster ', replacement = '')), size = 1.5, vjust = -0.8)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_spec_SA_abs_positions.png', path = savepath_outputs, # save as output
       width = 1500, height = 1080, units = 'px')

# Summarise data within 79 spectral clusters for both controls and synesthetes
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

# Summarise the clusters and export results
SA_abs_spec79_sum = summarise_clusters(SA_abs_lr, clust_spec_SA_abs_pos) # average within clusters for absolute SA
pcor(SA_abs_spec79_sum[[1]][SA_abs_spec79_sum[[1]]$Group == 'Syn', 8:ncol(SA_abs_spec79_sum[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_spec79_sum[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec79_sum.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec79_sum[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec79_sum.csv'), row.names = F) # write to .csv
write.csv(SA_abs_spec79_sum[[3]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec79_pos.csv'), row.names = F) # write cluster centroids to .csv

# Resample the controls 6 times without replacement for unique groups
set.seed(7842)
sample_numbers = sample(1:650, size = 6 * 102, replace = F) # generate random samples
sample_numbers = split(sample_numbers, ceiling(seq_along(sample_numbers)/102)) # split into 6 groups
for (i in seq(1:6)){
  SA_abs_lr_cntrl_resamp = SA_abs_lr[SA_abs_lr$Group == 'Control',][sample_numbers[[i]],]
  SA_abs_spec79_sum_102 = summarise_clusters(SA_abs_lr_cntrl_resamp, clust_spec_SA_abs_pos)
  write.csv(SA_abs_spec79_sum_102[[1]], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec79_sum_cntrl_resamp_', i, '.csv'), row.names = F) # write to .csv
}
write.csv(SA_abs_spec79_sum[[1]][SA_abs_spec79_sum[[1]]$Group == 'Syn',], file = paste0(savepath_data, 'spectral_clust/SA_abs_spec79_sum_syn.csv'), row.names = F) # write just the syn