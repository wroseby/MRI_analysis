# SCRIPT TO PERORM HIERARCHICAL CLUSTERING OF PARCELS AND PLOT RESULTS #

# Load packages
library(ppcor)
library(dplyr)
library(ggplot2)
library(reshape2)

# Function to perform hierarchical cluster and plot dendrogram
hclust.plot <- function(data, group, cortype = c('pearson', 'partial'), clust_method, cut){
 
  # generate correlation matrix - standard pearson or partial pearson
  
  if (cortype == 'pearson'){
    cormat = cor(x = data[data$Group == group, 8:ncol(data)], 
               y = data[data$Group == group, 8:ncol(data)]) 
  }
    
  if (cortype == 'partial'){
      cormat = pcor(x = data[data$Group == group, 8:ncol(data)],
           method = 'pearson')$estimate
  }
      
  distance = dist(cormat, method = 'euclidean') # euclidean distance: sum of squared differences; positive and negative correlations can be considered separately
  clust = hclust(distance, method = clust_method) # single method to do nearest neighbour (paired clusters); Ward's method to minimise within-cluster variance
  
  plot(clust, main = paste('Clustering:', group), sub = paste0("Euclidean distance, method = ", clust_method), # plot dendrogram
       cex = 1.3, cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.5, lwd = 2.0)  # adjust overall text size
  
  # acquire cluster identities
  if (clust_method == 'single'){
    clusters = cutree(clust, k = cut) # desired number of clusters if nearest neighbour
  }
  if (clust_method == 'Ward.D'){
    clusters = cutree(clust, h = cut) # level to cut at if Ward's method
  }
  
  clusters = data.frame(Area = names(clusters), cluster = paste0(rep('cluster ', length(clusters)), clusters))
  return(clusters)
} 

# Perform nearest neighbour (aka single linkage) clustering on control data and retrieve pairwise clusters
clust_nn_SA_abs = hclust.plot(SA_abs_lr, group = 'Control', cortype = 'partial', clust_method = 'single', cut = 90)

# Plot cluster anatomical positions
clust_nn_SA_abs = merge(parcel_positions, clust_nn_SA_abs) # merge clusters and positions
clust_nn_SA_abs$cluster = factor(clust_nn_SA_abs$cluster) # factor clusters for plotting

ggplot(data = clust_nn_SA_abs, aes(x = x, y = max(y) - y, fill = cluster))+ # plot parcel locations with fill colour for cluster
  geom_point(size = 1, shape = 21, stroke = 0.5)+
  geom_text(aes(label = cluster), size = 1.5, vjust = -0.8)+
  scale_fill_viridis_d(option = "C")+
  labs(x = 'x', y = 'y')+
  theme(legend.position = 'none',
        axis.text = element_blank())
ggsave(filename = 'clust_nn_SA_abs_positions.png', path = savepath_outputs, # save as output
       width = 1500, height = 1080, units = 'px')

# Average data within 90 nearest neighbour clusters for both controls and synesthetes
average_clusters = function(data, clusters){ # function to average across clusters
  data_long = melt(data[,c(1,8:ncol(data))]) # convert parcel data to long form for merging
  colnames(data_long) = c('ID', 'Area', 'value') # rename columns for merging
  data_long = merge(data_long, clusters, by = 'Area') # merge data with cluster assignments
  
  # average within clusters 
  data_averaged = data_long %>% 
    group_by(ID, cluster) %>%
    summarise(mean_value = mean(value), # average the SA value
              mean_x = mean(x), # average the x coordinate
              mean_y = mean(y)) # average the y coordinate
  
  # convert back to wide form and bind with metadata
  data_wide = dcast(data_averaged, ID ~ cluster, value.var = 'mean_value')
  data_wide = cbind(data[,2:7], data_wide)
  
  # get centroids of x and y positions for clusters
  data_pos = unique(dplyr::select(ungroup(data_averaged), 2,4,5)) # select cluster column, mean x and mean y columns
  
  return(list(data_wide, data_pos))
}

# Average the clusters and export results
SA_abs_nn90_avg = average_clusters(SA_abs_lr, clust_nn_SA_abs) # average within clusters for absolute SA
pcor(SA_abs_nn90_avg[[1]][SA_abs_nn90_avg[[1]]$Group == 'Syn', 8:ncol(SA_abs_nn90_avg[[1]])])$estimate # quick check of synesthete partial correlation
write.csv(SA_abs_nn90_avg[[1]], file = paste0(savepath_data, 'SA_abs_nn90_avg.csv'), row.names = F) # write to .csv
write.csv(SA_abs_nn90_avg[[2]], file = paste0(savepath_data, 'SA_abs_nn90_pos.csv'), row.names = F) # write cluster centroids to .csv
