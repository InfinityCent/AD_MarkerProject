library(dplyr)
library(naturalsort)

input_f = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/attached_feature_outputs"
output_f = "C:/Users/peree/OneDrive/Desktop/combined_cc_as"

as = c("with_AreaShape", "without_AreaShape")
markers = c("Cdc11", "Om45", "Sec21", "Spf1", "Vph1")

for (as_dir in as) {
  for (m in markers) {
    setwd(paste(input_f, as_dir, m, sep = "/"))
    
    combined_df = data.frame(matrix(nrow = 0, ncol = 13))
    colnames(combined_df) = c("cell_id", "ORF", "Name", "Allele", "Strain.ID", 
                              "Plate", "Row", "Column", "ImageNumber", 
                              "ObjectNumber", "Predictions", "Area", "Perimeter")
    
    m_files = naturalsort::naturalsort(list.files())
    
    for (m_file in m_files) {
      cluster_df = read.csv(m_file)
      combined_df = rbind(combined_df, cluster_df)
    }
    
    setwd(output_f)
    combined_df_file = paste(m, as_dir, "cc_as.csv", sep = "_")
    write.csv(combined_df, combined_df_file, row.names = FALSE)
  }
}


sec21 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cc_as/Sec21_without_AreaShape_cc_as.csv")
sec21_pca = read.csv("C:/Users/peree/OneDrive/Desktop/Sec21_all_PCA.csv")
sec21_pca_cell = dplyr::filter(sec21_pca, sec21_pca$cell_id == 13843)
sec21_pca_plate = dplyr::filter(sec21_pca, sec21_pca$Plate == "Sec21_TS2_37C_Plate02")
sec21_pca_plate_orf = dplyr::filter(sec21_pca_plate, sec21_pca_plate$ORF == "YER069W")
