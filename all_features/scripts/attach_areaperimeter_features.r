# Attach Area and Perimeter Features
# Purpose: prepare files needed to calculate cell size distributions

# Alex Daiejavad
# Date: January 24, 2023
#-------------------------------------------------------------------------------
library(naturalsort)
library(dplyr)

feature_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features"
marker_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/cell_cycle_files/cc_statistics"
output_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/attached_feature_outputs"

#feature_dir = "/home/morphology/shared/morphology/Users/Alex/areaperimeter/features"
#marker_dir = "/home/morphology/shared/morphology/Users/Alex//cc_statistics"
#output_dir = "/home/morphology/shared/morphology/Users/Alex/areaperimeter/outputs"

attach_features = function(marker, AS) { # AS = c('with_AreaShape', 'without_AreaShape')
  feature_file = paste(marker, '_AreaPerimeter.csv', sep = '')
  features = read.csv(paste(feature_dir, feature_file, sep = '/'))
  
  setwd(paste(marker_dir, AS, marker, sep = '/'))
  marker_clusters = naturalsort::naturalsort(list.files())
  
  for (mc in marker_clusters) {
    mc_table = read.csv(mc)
    area_features = c()
    perimeter_features = c()
    
    for (r in 1:nrow(mc_table)) {
      cellid = mc_table[r, 1]
      filter_features = dplyr::filter(features, features$cell_id == cellid)
      area_features = c(area_features, filter_features$Cells_AreaShape_Area)
      perimeter_features = c(perimeter_features, filter_features$Cells_AreaShape_Perimeter)
    }
    mc_table = cbind(mc_table, Area = area_features, Perimeter = perimeter_features)
    
    marker_output_dir = paste(output_dir, AS, marker, sep = "/")
    write.csv(mc_table, file = paste(marker_output_dir, mc, sep = "/"), row.names = FALSE)
    print(mc)
  }
}


attach_features_vph1 = function(marker, AS) { # AS = c('with_AreaShape', 'without_AreaShape')
  feature_file = paste(marker, '_AreaPerimeter.csv', sep = '')
  features = read.csv(paste(feature_dir, feature_file, sep = '/'))
  
  setwd(paste(marker_dir, AS, marker, sep = '/'))
  marker_clusters = c('normal_clusters_n12_per_cell_cluster-0.csv', 
                      'normal_clusters_n12_per_cell_cluster-1.csv', 
                      'normal_clusters_n12_per_cell_cluster-2.csv')
  
  for (mc in marker_clusters) {
    mc_table = read.csv(mc)
    area_features = c()
    perimeter_features = c()
    
    for (r in 1:nrow(mc_table)) {
      cellid = mc_table[r, 1]
      filter_features = dplyr::filter(features, features$cell_id == cellid)
      area_features = c(area_features, filter_features$Cells_AreaShape_Area)
      perimeter_features = c(perimeter_features, filter_features$Cells_AreaShape_Perimeter)
    }
    mc_table = cbind(mc_table, Area = area_features, Perimeter = perimeter_features)
    
    marker_output_dir = paste(output_dir, AS, marker, sep = "/")
    write.csv(mc_table, file = paste(marker_output_dir, mc, sep = "/"), row.names = FALSE)
    print(mc)
  }
}

attach_features("Cdc11", "with_AreaShape")
attach_features("Cdc11", "without_AreaShape")

attach_features("Om45", "with_AreaShape")
attach_features("Om45", "without_AreaShape")

attach_features("Sec21", "with_AreaShape")
attach_features("Sec21", "without_AreaShape")

attach_features("Spf1", "with_AreaShape")
attach_features("Spf1", "without_AreaShape")

#attach_features("Vph1", "with_AreaShape")
#attach_features("Vph1", "without_AreaShape")

#attach_features_vph1("Vph1", "with_AreaShape")