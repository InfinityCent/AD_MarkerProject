# Marker-cluster size analysis
# Purpose: get the average cell area/perimeter for each marker-cluster to see 
# if there are 'big cell' or 'small cell' clusters between wAS and woAS datasets

# Date: February 13, 2023
# Alex Daiejavad
#-------------------------------------------------------------------------------
library(naturalsort)
library(reshape)
library(ggplot2)

input_files = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/attached_feature_outputs"
markers = c("Cdc11", "Om45", "Sec21", "Spf1", "Vph1") # only those with wAS/woAS data

avg_area_perimeter = function(AS_data) {
  all_markers_df = data.frame(matrix(nrow = 0, ncol = 13))
  
  for (m in markers) {
    setwd(paste(input_files, AS_data, m, sep = "/"))
    clusters = naturalsort::naturalsort(list.files())
    
    avg_area = c(paste(m, "Area", sep = " "))
    avg_perimeter = c(paste(m, "Perimeter", sep = " "))
    for (c in clusters) {
      cluster_df = read.csv(c)
      avg_area = c(avg_area, mean(cluster_df$Area))
      avg_perimeter = c(avg_perimeter, mean(cluster_df$Perimeter))
    }
    # In case length of vectors < 13
    avg_area = c(avg_area, rep(NA, 13 - length(avg_area)))
    avg_perimeter = c(avg_perimeter, rep(NA, 13 - length(avg_perimeter)))
    
    all_markers_df = rbind(all_markers_df, avg_area, avg_perimeter)
  }
  df_cols = c("Marker-Feature", paste("Cluster", 0:11, sep = " "))
  colnames(all_markers_df) = df_cols
  
  return(all_markers_df)
}


med_area_perimeter = function(AS_data) {
  all_markers_df = data.frame(matrix(nrow = 0, ncol = 13))
  
  for (m in markers) {
    setwd(paste(input_files, AS_data, m, sep = "/"))
    clusters = naturalsort::naturalsort(list.files())
    
    med_area = c(paste(m, "Area", sep = " "))
    med_perimeter = c(paste(m, "Perimeter", sep = " "))
    for (c in clusters) {
      cluster_df = read.csv(c)
      med_area = c(med_area, median(cluster_df$Area))
      med_perimeter = c(med_perimeter, median(cluster_df$Perimeter))
    }
    # In case length of vectors < 13
    med_area = c(med_area, rep(NA, 13 - length(med_area)))
    med_perimeter = c(med_perimeter, rep(NA, 13 - length(med_perimeter)))
    
    all_markers_df = rbind(all_markers_df, med_area, med_perimeter)
  }
  df_cols = c("Marker-Feature", paste("Cluster", 0:11, sep = " "))
  colnames(all_markers_df) = df_cols
  
  return(all_markers_df)
}



wAS_mean = avg_area_perimeter("with_AreaShape")
woAS_mean = avg_area_perimeter("without_AreaShape")

wAS_med = med_area_perimeter("with_AreaShape")
woAS_med = med_area_perimeter("without_AreaShape")

setwd(input_files)
write.csv(wAS_mean, "wAS_avg_area_perimeter.csv", row.names = FALSE)
write.csv(woAS_mean, "woAS_avg_area_perimeter.csv", row.names = FALSE)

#-------------------------------------------------------------------------------

create_barplots_mean = function(marker, r, wAS, woAS) {
  wAS = as.numeric(wAS[r, 2:13])
  wAS = sort(replace(wAS, is.na(wAS), 0))
  
  woAS = as.numeric(woAS[r, 2:13])
  woAS = sort(replace(woAS, is.na(woAS), 0))
  
  positions = c(paste("Position", 1:12, sep = " "))
  
  to_plot = data.frame(positions, wAS, woAS)
  melted = melt(to_plot, id = "positions")
  
  setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/areaperimeter_plots")
  if (r %in% c(1, 3, 5, 7, 9)) {
    fn = 'MC_Average_Area'
    feature = 'Area'}
  
  if (r %in% c(2, 4, 6, 8, 10)) {
    fn = 'MC_Average_Perimeter'
    feature = 'Perimeter'}
  
  file_name = paste(marker, "_", fn, ".png", sep = "")
  
  print(ggplot(melted, aes(x = positions, y = value, fill = variable)) + 
          geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
          scale_x_discrete(limits = positions) +
          scale_y_continuous(breaks = seq(round(min(melted$value - 0.1)), 
                                          round(max(melted$value + 0.1)), 
                                          by = 0.2)) +
          scale_fill_discrete(labels = c("With AreaShape", "Without AreaShape")) +
          labs(title = paste(marker, "Clusters: Average", feature, sep = " "), 
               y = paste("Average ", feature, sep = ""), 
               x = "Positions",
               fill = "Dataset") +
          theme(plot.title = element_text(hjust = 0.5)))
  
  
  ggsave(file_name, width = 26, height = 10, units = 'cm')
  dev.off()
  
}


create_barplots_median = function(marker, r, wAS, woAS) {
  wAS = as.numeric(wAS[r, 2:13])
  wAS = sort(replace(wAS, is.na(wAS), 0))
  
  woAS = as.numeric(woAS[r, 2:13])
  woAS = sort(replace(woAS, is.na(woAS), 0))
  
  positions = c(paste("Position", 1:12, sep = " "))
  
  to_plot = data.frame(positions, wAS, woAS)
  melted = melt(to_plot, id = "positions")
  
  setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/areaperimeter_plots")
  if (r %in% c(1, 3, 5, 7, 9)) {
    fn = 'MC_Median_Area'
    feature = 'Area'}
  
  if (r %in% c(2, 4, 6, 8, 10)) {
    fn = 'MC_Median_Perimeter'
    feature = 'Perimeter'}
  
  file_name = paste(marker, "_", fn, ".png", sep = "")
  
  print(ggplot(melted, aes(x = positions, y = value, fill = variable)) + 
          geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
          scale_x_discrete(limits = positions) +
          scale_y_continuous(breaks = seq(round(min(melted$value - 0.1)), 
                                          round(max(melted$value + 0.1)), 
                                          by = 0.2)) +
          scale_fill_discrete(labels = c("With AreaShape", "Without AreaShape")) +
          labs(title = paste(marker, "Clusters: Median", feature, sep = " "), 
               y = paste("Median ", feature, sep = ""), 
               x = "Positions",
               fill = "Dataset") +
          theme(plot.title = element_text(hjust = 0.5)))
  
  
  ggsave(file_name, width = 26, height = 10, units = 'cm')
  dev.off()
  
}


create_barplots_mean("Cdc11", 1, wAS_mean, woAS_mean)
create_barplots_mean("Cdc11", 2, wAS_mean, woAS_mean)
create_barplots_mean("Om45", 3, wAS_mean, woAS_mean)
create_barplots_mean("Om45", 4, wAS_mean, woAS_mean)
create_barplots_mean("Sec21", 5, wAS_mean, woAS_mean)
create_barplots_mean("Sec21", 6, wAS_mean, woAS_mean)
create_barplots_mean("Spf1", 7, wAS_mean, woAS_mean)
create_barplots_mean("Spf1", 8, wAS_mean, woAS_mean)
create_barplots_mean("Vph1", 9, wAS_mean, woAS_mean)
create_barplots_mean("Vph1", 10, wAS_mean, woAS_mean)

create_barplots_median("Cdc11", 1, wAS_med, woAS_med)
create_barplots_median("Cdc11", 2, wAS_med, woAS_med)
create_barplots_median("Om45", 3, wAS_med, woAS_med)
create_barplots_median("Om45", 4, wAS_med, woAS_med)
create_barplots_median("Sec21", 5, wAS_med, woAS_med)
create_barplots_median("Sec21", 6, wAS_med, woAS_med)
create_barplots_median("Spf1", 7, wAS_med, woAS_med)
create_barplots_median("Spf1", 8, wAS_med, woAS_med)
create_barplots_median("Vph1", 9, wAS_med, woAS_med)
create_barplots_median("Vph1", 10, wAS_med, woAS_med)