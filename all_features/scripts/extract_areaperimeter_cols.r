feature_dir = "/home/morphology/mpg6/mojcamu/Marker_Project/PCA_features"
output_dir = "/home/morphology/shared/morphology/Users/Alex/areaperimeter/features2"

feature_files = c("Cdc11_all_CP.csv", "Om45_all_CP.csv", "Sec21_all_CP.csv", 
                  "Spf1_all_CP.csv", "Vph1_all_CP.csv", )

for (ff in feature_files) {
  marker = stringr::str_split(ff, "_")[[1]][1]
  marker_features = read.csv(ff)
  feature_subset = data.frame(marker_features$cell_id, 
                              marker_features$Cells_AreaShape_Area,
                              marker_features$Cells_AreaShape_Perimeter)
  colnames(feature_subset) = c("cell_id", "Cells_AreaShape_Area", "Cells_AreaShape_Perimeter")
  setwd(output_dir)
  feature_subset_name = paste(marker, "AreaPerimeter.csv", sep = "_")
  write.csv(feature_subset, feature_subset_name)
}