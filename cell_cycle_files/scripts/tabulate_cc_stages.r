# Create a table of wAS and woAS marker-clusters showing proportion of cells 
# falling into each CC stage

# Alex Daiejavad
# Date: January 21, 2023
#-------------------------------------------------------------------------------
library(naturalsort)
dir = "C:/Users/peree/OneDrive/Desktop/cc_statistics"

AS_splits = c("with_AreaShape", "without_AreaShape")
markers = c("Cdc11", "Om45", "Sec21", "Spf1", "Vph1")

cc_tabulate = function(AS, marker, cc_stages) {
  setwd(paste(dir, AS, marker, sep = "/"))
  marker_clusters = naturalsort::naturalsort(list.files())
  summary_cc = data.frame(matrix(nrow = length(cc_stages), ncol = 0))
  
  for (mc in marker_clusters) {
    mc_table = read.csv(mc, na.strings = "")
    stage_counts = c()
    
    for (stage in cc_stages) {
      stage_counts = c(stage_counts, sum((mc_table$Predictions == stage), na.rm = TRUE))
    }
    
    names(stage_counts) = cc_stages
    summary_cc = cbind(summary_cc, stage_counts)
  }
  colnames(summary_cc) = paste("Cluster", seq(from = 0, to = length(marker_clusters) - 1, by = 1))
  return(summary_cc)
}


cc_tabulate_overall = function(all_preds, cc_labels) {
  cc_counts = c()
  for (cc_label in cc_labels) {
    cc_counts = c(cc_counts, sum(all_preds$Prediction == cc_label))
  }
  names(cc_counts) = cc_labels
  return(cc_counts)
}

cdc11 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cell_cycle/Cdc11_all_predictions.csv")
om45 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cell_cycle/Om45_all_predictions.csv")
sec21 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cell_cycle/Sec21_all_predictions.csv")
spf1 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cell_cycle/Spf1_all_predictions.csv")
vph1 = read.csv("C:/Users/peree/OneDrive/Desktop/combined_cell_cycle/Vph1_all_mpg3_predictions.csv")

#cdc11_cc = cc_tabulate("without_AreaShape", "Cdc11", unique(cdc11$Prediction))
#om45_cc = cc_tabulate("without_AreaShape", "Om45", unique(om45$Prediction))
#sec21_cc = cc_tabulate("without_AreaShape", "Sec21", unique(sec21$Prediction))
#spf1_cc = cc_tabulate("without_AreaShape", "Spf1", unique(spf1$Prediction))
vph1_cc = cc_tabulate("without_AreaShape", "Vph1", unique(vph1$Prediction))

#cdc11_cc_wAS = cc_tabulate("with_AreaShape", "Cdc11", unique(cdc11$Prediction))
#om45_cc_wAS = cc_tabulate("with_AreaShape", "Om45", unique(om45$Prediction))
#sec21_cc = cc_tabulate("with_AreaShape", "Sec21", unique(sec21$Prediction))
#spf1_cc_wAS = cc_tabulate("with_AreaShape", "Spf1", unique(spf1$Prediction))
vph1_cc_wAS = cc_tabulate("with_AreaShape", "Vph1", unique(vph1$Prediction))


setwd(dir)
#write.csv(cdc11_cc, file = "Cdc11_CellCycle.csv")
#write.csv(om45_cc, file = "Om45_CellCycle.csv")
#write.csv(sec21_cc, file = "Sec21_CellCycle.csv")
#write.csv(spf1_cc, file = "Spf1_CellCycle.csv")
write.csv(vph1_cc, file = "Vph1_CellCycle_woAS.csv")

#write.csv(cdc11_cc_wAS, file = "Cdc11_CellCycle_wAS.csv")
#write.csv(om45_cc_wAS, file = "Om45_CellCycle_wAS.csv")
#write.csv(sec21_cc, file = "Sec21_CellCycle_wAS.csv")
#write.csv(spf1_cc_wAS, file = "Spf1_CellCycle_wAS.csv")
write.csv(vph1_cc_wAS, file = "Vph1_CellCycle_wAS.csv")

cdc11_tab = cc_tabulate_overall(cdc11, unique(cdc11$Prediction))
om45_tab = cc_tabulate_overall(om45, unique(om45$Prediction))
sec21_tab = cc_tabulate_overall(sec21, unique(sec21$Prediction))
spf1_tab = cc_tabulate_overall(spf1, unique(spf1$Prediction))
vph1_tab = cc_tabulate_overall(vph1, unique(vph1$Prediction))
