# Attach cell cycle predictions to each marker-cluster
# Purpose: knowing the distribution of CC predictions can determine if there's 
# CC bias in clustering data

# Alex Daiejavad
# Date: January 20, 2023
#-------------------------------------------------------------------------------
library(dplyr)

cc_predictions_dir = "C:/Users/peree/OneDrive/Desktop/combined_cell_cycle"
input_dir = "C:/Users/peree/OneDrive/Desktop/clustering_profiles_recluster_after_applying_thresh1"
output_dir = "C:/Users/peree/OneDrive/Desktop/cc_statistics"

attach_cc_predictions = function(marker, AS) { #AS = c("with_AreaShape", "without_AreaShape")
  all_cc_preds_file = paste(marker, "_all_predictions.csv", sep = "")
  all_cc_preds = read.csv(paste(cc_predictions_dir, all_cc_preds_file, sep = "/"))
  
  AS_dir = NA
  if (AS == "with_AreaShape") {AS_dir = "_all_PCA_combinedMA"}
  if (AS == "without_AreaShape") {AS_dir = "_all_PCA_combinedMA_noAreaShape"}
  AS_dir = paste(marker, AS_dir, sep = "")
  
  setwd(paste(input_dir, AS, AS_dir, "normal_clusters", "Clusters", sep = "/"))
  clusters = list.files()
  
  
  for (c in clusters) {
    cc_predictions = c()
    cluster_cells = read.csv(c)
    
    for (r in 1:nrow(cluster_cells)) {
      cell_id = cluster_cells[r, 1]
      cell_cc_pred = dplyr::filter(all_cc_preds, all_cc_preds$Cell_ID == cell_id)$Prediction
      
      if (length(cell_cc_pred) == 0) {cc_predictions = c(cc_predictions, NA)}
      if (length(cell_cc_pred) > 0) {cc_predictions = c(cc_predictions, cell_cc_pred)}
      #cc_predictions = c(cc_predictions, cell_cc_pred)
    }
    
    marker_output_dir = paste(output_dir, AS, marker, sep = "/")
    cluster_cells = cbind(cluster_cells, Predictions = cc_predictions)
    write.csv(cluster_cells, file = paste(marker_output_dir, c, sep = "/"), row.names = FALSE)
    print(c)
  }
}

#attach_cc_predictions("Cdc11", "without_AreaShape")
#attach_cc_predictions("Om45", "without_AreaShape")
#attach_cc_predictions("Spf1", "without_AreaShape")
#attach_cc_predictions("Sec21", "without_AreaShape")
attach_cc_predictions("Vph1", "without_AreaShape")

#attach_cc_predictions("Cdc11", "with_AreaShape")
#attach_cc_predictions("Om45", "with_AreaShape")
#attach_cc_predictions("Spf1", "with_AreaShape")
#attach_cc_predictions("Sec21", "with_AreaShape")
attach_cc_predictions("Vph1", "with_AreaShape")