# Feature PCC
# Purpose: create a PCC matrix for each batch of features for each marker

# Alex Daiejavad
# Date: January 10, 2023
#-------------------------------------------------------------------------------

library(naturalsort)
library(Hmisc)
library(corrplot)

setwd("C:/Users/peree/OneDrive/Desktop/Om45_wAS_splits")  # Windows
setwd("~/Documents/split_feature_files/wAS")  # Linux
setwd("E:/alex_feature_files/split_features") # Seagate

directory = "E:/alex_feature_files/split_features"

feature_pcc = function(marker, feature_group) {
  filename = paste(marker, feature_group, "features.csv", sep = "_")
  cell_features = read.csv(filename)
  
  cf_trunc = cell_features[ , 11:ncol(cell_features)]  # remove cell info columns
  ordered_colname = naturalsort::naturalsort(colnames(cf_trunc))
  cf_trunc = cf_trunc[ , ordered_colname]  # rearrange columns so they're sorted
  
  pcc_matrix = Hmisc::rcorr(as.matrix(cf_trunc), type = c("pearson"))$r
  
  return(pcc_matrix)
}

markers = c("Cdc11")
# feature_groups = c("AreaShape", "Granularity", "Intensity", "RadialDistribution", "Texture")
feature_groups = c("Texture1", "Texture2")

for (m in markers) {
  setwd(paste(directory, m, sep = "/"))
  
  for (fg in feature_groups) {
    pcc_matrix = feature_pcc(marker = m, feature_group = fg)
    pcc_matrix_file = paste(m, fg, "PCC.csv", sep = "_")
    pcc_matrix_file = paste(directory, pcc_matrix_file, sep = "/")
    
    write.csv(pcc_matrix, file = pcc_matrix_file)
   }
}

output_dir = "E:/alex_feature_files/split_features/Outputs"
corrplotter = function(output_dir, marker) {
  setwd(output_dir)
  all_files = list.files()
  
  for (f in all_files) {
    pcc_matrix = read.csv(f, row.names = 1)
    fg = strsplit(f, split = "_")[[1]][2]
    plot_title = paste(marker, " | Pearson Correlation Coefficients (", fg, ")" , sep = "")
    plot_file = paste(marker, fg, "plot.png", sep = "_")
    
    png(file = plot_file, width = 948, height = 874)
    
    print(corrplot::corrplot(as.matrix(pcc_matrix), 
                             method = c("color"),
                             title = plot_title,
                             mar=c(0,0,2,0),
                             outline = TRUE, cl.pos = 'b', 
                             tl.pos = 'n'))
    
    dev.off()
  }
}

corrplotter(output_dir, "Cdc11")



setwd(directory)

plot_title = paste(marker, " | Pearson Correlation Coefficients (", feature_group, ")" , sep = "")
corrplot::corrplot(pcc_matrix, 
                   method = c("color"),
                   title = plot_title,
                   mar=c(0,0,2,0),
                   outline = TRUE)

