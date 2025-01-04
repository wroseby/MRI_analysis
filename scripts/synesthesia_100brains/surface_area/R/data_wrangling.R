# SCRIPT TO LOAD AND ORGAISE SURFACE AREA DATA FOR SYNESTHETES AND CONTROLS #

library(dplyr)

# Load data
datapath = 'D:/Documents/Academia/projects/ward_lab/MRI_analysis/datasets/synesthesia_100brains/' # path to data
SA_abs = read.csv(paste0(datapath, '/surface_area/biomarker_S1C_harm.csv'), row.names = 1) # read the csv for absolute SA
SA_rel = read.csv(paste0(datapath, '/surface_area/biomarker_S1A_harm.csv'), row.names = 1) # read the csv for relative SA
parcel_positions = read.csv(paste0(datapath, '/common/parcel_positions_flat.csv'))

# Set the path for saving any intermediate data (shared) or plots (outputs)
savepath_data = 'D:/Documents/Academia/projects/ward_lab/MRI_analysis/shared/synesthesia_100brains/surface_area/'
savepath_outputs = 'D:/Documents/Academia/projects/ward_lab/MRI_analysis/outputs/synesthesia_100brains/surface_area/R/'

# Average data across left and right parcels
lr_average = function(data){ # function to perform lr averaging
  parcel_data = data[,grepl(pattern='^[L,R]_', x=colnames(data))] # get data for parcels only (L_ and R_ columns)
  parcels_LR = unique(sub("^[LR]_", "", colnames(parcel_data))) # get unique parcels
  data_lr = list() # initialise list
  for (parcel in parcels_LR) {
    l_col = paste0("L_", parcel)
    r_col = paste0("R_", parcel)
    data_lr[[parcel]] = rowMeans(parcel_data[, c(l_col, r_col)], na.rm = TRUE) # get means of left and right parcels
  }
  data_lr = as.data.frame(data_lr) # convert list to dataframe
  data_lr = cbind(data[,!grepl(pattern='^[L,R]_', x=colnames(data))], data_lr) # bind averaged parcel data with metadata
  colnames(data_lr) = gsub(x = colnames(data_lr), pattern = '^X', replacement = '') # remove leading Xs in column names with numbers
  
  return(data_lr)
}

SA_abs_lr = lr_average(SA_abs)
SA_rel_lr = lr_average(SA_rel)