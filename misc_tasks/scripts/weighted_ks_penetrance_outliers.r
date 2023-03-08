# Calculate weighted penetrances
# Purpose: this is a new approach we're taking for identifying high-confidence
# outlier cells

# Alex Daiejavad
# March 03, 2023
#-------------------------------------------------------------------------------
library(dplyr)

input_dir = "/home/morphology/mpg6/mojcamu/Marker_Project/OD_results/"
output_dir = "/home/morphology/shared/morphology/Users/Alex/weighted_pen_outliers/"

markers = c('Cdc11', 'Dad2', 'Heh2', 'Hta2', 'Nop10', 'Nuf2', 'Om45', 'Pil1', 
            'Psr1', 'Sac6', 'Sec7', 'Sec21', 'Snf7', 'Spf1', 'Vph1')

# num_cells and ks_pen are vectors; calculate weighted_penetrance using 
# vectorization
weighted_penetrance = function(num_cells, ks_pen) {
  wp = sum(num_cells * ks_pen) / sum(num_cells)
  
  return(wp)
}


get_wp_outliers = function(marker) {
  
  # Step 1: Calculate weighted penetrances and attach to per_strain data ---------
  setwd(paste(input_dir, marker, "_all_PCA", sep = ""))
  
  per_well = read.csv(paste(marker, "_all_PCA_GMM_OD_results_well.csv", sep = ""))
  per_strain = read.csv(paste(marker, "_all_PCA_GMM_OD_results_strain.csv", sep = ""))
  strains = unique(per_strain$Strain.ID)
  
  weighted_penetrances = c()
  for (strain in strains) {
    strain_subset = dplyr::filter(per_well, per_well$Strain.ID == strain)
    weighted_penetrances = c(weighted_penetrances, 
                             weighted_penetrance(strain_subset$Num_cells, strain_subset$KS_Penetrance))
  }
  
  per_strain = cbind(per_strain, weighted_penetrance = weighted_penetrances)
  
  
  # Step 2: Compare each strain's weighted penetrance to its overall KS_Penetrance
  
  ks_below_weighted = per_strain$KS_Penetrance < per_strain$weighted_penetrance
  per_strain = cbind(per_strain, below_threshold = ks_below_weighted)
  
  strains_to_be_kept = dplyr::filter(per_strain, per_strain$below_threshold == FALSE)
  strains_to_be_kept = strains_to_be_kept$Strain.ID
  
  # Step 3: From the list of outlier cells, remove those that come from strains
  #         filtered out
  
  od_results = read.csv(paste(marker, "_all_PCA_GMM_OD_results_score_strain.csv", sep = ""))
  od_results_kept = dplyr::filter(od_results, od_results$Strain.ID %in% strains_to_be_kept)
  
  kept_outliers = dplyr::filter(od_results_kept, od_results_kept$Is_outlier == "True")
  
  
  # Step 4: Export outlier cells
  write.csv(kept_outliers, paste(output_dir, marker, 
                                 "_weighted_penetrance_outliers.csv", sep = ""), 
            row.names = FALSE)
}


for (m in markers){
  get_wp_outliers(m)
}