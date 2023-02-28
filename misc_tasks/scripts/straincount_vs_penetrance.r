# Create # strains vs. penetrance after Threshold 0
# Purpose: this is something Thanasis asked for after my lab meeting

# Date: February 22, 2023
# Alex Daiejavad
#--------------------------------------------------------------------
library(dbplyr)
library(ggplot2)

od_outputs = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/myra_tasks/strain_penetrance_hist/OD_results"
strains_removed = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/myra_tasks/strain_penetrance_hist/thr0_low_strain_cellcount"
plot_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/myra_tasks/strain_penetrance_hist/plots"

od_files = list.files(od_outputs)
strain_files = list.files(strains_removed)

for (i in 1:length(od_files)) {
  marker_od = read.csv(paste(od_outputs, od_files[i], sep = "/"))
  marker_strain = read.csv(paste(strains_removed, strain_files[i], sep = "/"))
  
  marker_od_filtered = dplyr::filter(marker_od, !(marker_od$ORF %in% marker_strain$ORF))
  
  plot_title = paste(strsplit(od_files[i], split = "_")[[1]][1],
                     "Strain Count vs. Penetrance",
                     sep = " ")
  
  my_plot = ggplot(marker_od_filtered, aes(x = Penetrance)) +
    geom_histogram(color = "white", binwidth = 5, bins = 20, 
                   center = 2.5) +
    scale_x_continuous(breaks = round(seq(0, 100, by = 5),1)) +
    scale_y_continuous(breaks = round(seq(0, 7000, by = 500),1)) +
    ylab("Strain Count") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 16))
  
  file_name = paste(strsplit(od_files[i], split = "_")[[1]][1],
                    "histogram.png",
                    sep = "_")
  
  file_name = paste(plot_dir, file_name, sep = "/")
  
  ggsave(file_name, plot = my_plot, width = 250, height = 250, units = "mm")
}
