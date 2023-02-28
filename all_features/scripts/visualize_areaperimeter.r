# Create Distributions of Area and Perimeter Features
# Purpose: Distributions are used to see if there's any difference between wAS and 
# woAS

# Alex Daiejavad
# Date: January 25, 2023
#-------------------------------------------------------------------------------
library(naturalsort)
library(ggplot2)

input_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/attached_feature_outputs"

# This function reads all of the individual marker-cluster csv files with area 
# and perimeter feature info, and combines them into one giant df
combine_clusters = function(marker, AS) {
  setwd(paste(input_dir, AS, marker, sep = "/"))
  marker_clusters = naturalsort::naturalsort(list.files())
  
  combined_df = data.frame(matrix(ncol = 2, nrow = 0))
  
  for (mc in marker_clusters) {
    mc_table = read.csv(mc)
    combined_df = rbind(combined_df, mc_table[ , 12:13])
  }
  
  colnames(combined_df) = c("Area", "Perimeter")
  return(combined_df)
}

#-------------------------------------------------------------------------------

#cdc11_wAS = combine_clusters("Cdc11", "with_AreaShape")
#om45_wAS = combine_clusters("Om45", "with_AreaShape")
#sec21_wAS = combine_clusters("Sec21", "with_AreaShape")
#spf1_wAS = combine_clusters("Spf1", "with_AreaShape")
vph1_wAS = combine_clusters("Vph1", "with_AreaShape")

#cdc11_woAS = combine_clusters("Cdc11", "without_AreaShape")
#om45_woAS = combine_clusters("Om45", "without_AreaShape")
#sec21_woAS = combine_clusters("Sec21", "without_AreaShape")
#spf1_woAS = combine_clusters("Spf1", "without_AreaShape")
vph1_woAS = combine_clusters("Vph1", "without_AreaShape")

#-------------------------------------------------------------------------------

# Since wAS and woAS have different numbers of observations, historgram is not a 
# good idea. Density plots are much better for this purpose. They're also better
# than scatterplots because there are so many individual datapoints for scatterplots
# to actually show much info.
# https://towardsdatascience.com/how-to-compare-two-or-more-distributions-9b06ee4d30bf
create_densityplots = function(marker, wAS, woAS, feature_num) { # feature = 1 is Area, 2 is Perimeter
  setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/areaperimeter_features/areaperimeter_plots")
  
  if (feature_num == 1) {feature = "Area"}
  if(feature_num == 2) {feature = "Perimeter"}
  file_name = paste(marker, "_", feature, ".png", sep = "")
  
  print(ggplot() +
          geom_density(wAS, mapping = aes(x = wAS[, feature_num], fill = "With AreaShape"), alpha = 0.25) +
          geom_density(woAS, mapping = aes(x = woAS[, feature_num], fill = "Without AreaShape"), alpha = 0.25) + 
          labs(x = paste(feature, "Feature Values", sep = " "), y = "Density", fill = "Dataset") +
          ggtitle(paste(marker, feature, "Measurements", sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5)))
  
  ggsave(file_name, width = 15, height = 10, units = 'cm')
  dev.off()
}

#create_densityplots("Cdc11", cdc11_wAS, cdc11_woAS, 1)
#create_densityplots("Om45", om45_wAS, om45_woAS, 1)
#create_densityplots("Sec21", sec21_wAS, sec21_woAS, 1)
#create_densityplots("Spf1", spf1_wAS, spf1_woAS, 1)
create_densityplots("Vph1", vph1_wAS, vph1_woAS, 1)

#create_densityplots("Cdc11", cdc11_wAS, cdc11_woAS, 2)
#create_densityplots("Om45", om45_wAS, om45_woAS, 2)
#create_densityplots("Sec21", sec21_wAS, sec21_woAS, 2)
#create_densityplots("Spf1", spf1_wAS, spf1_woAS, 2)
create_densityplots("Vph1", vph1_wAS, vph1_woAS, 2)

#-------------------------------------------------------------------------------

# I can't use t-tests to compare the distributions of wAS and woAS for each marker
# because t-test assumes that distributions are normal. For our datasets, the
# distributions are right-skewed.
# https://www.statology.org/determine-equal-or-unequal-variance/

# Our data is non-parametric (not normal) so it seems that Mann-Whitney U Test
# (aka Wilcoxon Rank Sum Test) is the best option here.
# https://www.technologynetworks.com/informatics/articles/mann-whitney-u-test-assumptions-and-example-363425

wilcox.test(cdc11_wAS$Area, cdc11_woAS$Area)
wilcox.test(cdc11_wAS$Perimeter, cdc11_woAS$Perimeter)

wilcox.test(om45_wAS$Area, om45_woAS$Area)
wilcox.test(om45_wAS$Perimeter, om45_woAS$Perimeter)

wilcox.test(sec21_wAS$Area, sec21_woAS$Area)
wilcox.test(sec21_wAS$Perimeter, sec21_woAS$Perimeter)

wilcox.test(spf1_wAS$Area, spf1_woAS$Area)
wilcox.test(spf1_wAS$Perimeter, spf1_woAS$Perimeter)
