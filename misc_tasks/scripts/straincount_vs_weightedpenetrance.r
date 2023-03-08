# Create # strains vs. weighted penetrance after Threshold 0
# Purpose: when doing wt percentile penetrance, there's a huge peak in the 0-10%
#          bin. We want to see if doing weighted_penetrance eliminates this issue

# Date: February 22, 2023
# Alex Daiejavad
#--------------------------------------------------------------------
library(dbplyr)
library(ggplot2)

od_outputs = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/misc_tasks/strain_penetrance_hist/OD_results"
high_conf_strains = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/misc_tasks/strain_penetrance_hist/weighted_pen_outliers"
plot_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/misc_tasks/strain_penetrance_hist/plots_weighted_pen"

od_files = list.files(od_outputs)
strain_files = list.files(high_conf_strains)

for (i in 1:length(od_files)) {
  marker_od = read.csv(paste(od_outputs, od_files[i], sep = "/"))
  marker_strain = read.csv(paste(high_conf_strains, strain_files[i], sep = "/"))
  included_strains = unique(marker_strain$Strain.ID)
  
  marker_od_filtered = dplyr::filter(marker_od, marker_od$Strain.ID %in% included_strains)
  
  plot_title = paste(strsplit(od_files[i], split = "_")[[1]][1],
                     "Strain Count vs. Weighted Penetrance",
                     sep = " ")
  
  my_plot = ggplot(marker_od_filtered, aes(x = KS_Penetrance)) +
    geom_histogram(color = "white", binwidth = 5, bins = 20, 
                   center = 2.5) +
    scale_x_continuous(breaks = round(seq(0, 100, by = 5),1)) +
    ylab("Strain Count") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 16))
  
  file_name = paste(strsplit(od_files[i], split = "_")[[1]][1],
                    "weighted_pen.png",
                    sep = "_")
  
  file_name = paste(plot_dir, file_name, sep = "/")
  
  ggsave(file_name, plot = my_plot, width = 250, height = 250, units = "mm")
}
