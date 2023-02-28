# Create a heatmap showing included and excluded features for each marker
# Purpose: get an overall understanding of how many features were included for
# analysis in each marker

# Date: February 10, 2023
# Alex Daiejavad
#-------------------------------------------------------------------------------
library(naturalsort)
library(pheatmap)
library(RColorBrewer)

all_features = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/marker_features_all/wAS"
included_features = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/marker_features_original_clustering/wAS_features_for_clustering"
excluded_features = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/all_features/marker_features_excluded_from_clustering/wAS"

# Features used in OD
marker_files = naturalsort::naturalsort(list.files(path = excluded_features))
marker_files = marker_files[!(marker_files == 'Heh2.txt')]
marker_files = marker_files[!(marker_files == 'Hta2.txt')]

all_marker_features_df = data.frame(matrix(ncol = 0, nrow = 184))

for (mf in marker_files) {
  included = read.delim(paste(included_features, "/", mf, sep = ""), header = FALSE)
  excluded = read.delim(paste(excluded_features, "/", mf, sep = ""), header = FALSE)
  inc_exc = c(c(rep('included', nrow(included))), c(rep('excluded', nrow(excluded))))
  
  all_features_df = data.frame(cbind(inc_exc, c(included[ , 1], excluded[ , 1])))
  colnames(all_features_df) = c("Included/Excluded", "Feature")
  all_features_df = all_features_df[naturalorder(all_features_df$Feature), ]
  
  inc_exc = all_features_df$`Included/Excluded`
  inc_exc = as.integer(inc_exc == 'included')
  
  all_marker_features_df = cbind(all_marker_features_df, inc_exc)
}

rownames(all_marker_features_df) = rbind(all_features_df$Feature)
colnames(all_marker_features_df) = c('Cdc11', 'Dad2', 'Nop10', 'Nuf2', 
                                     'Om45', 'Pil1', 'Psr1', 'Rad52', 'Sac6', 
                                     'Sec7', 'Sec21', 'Spf1', 'Vph1')


setwd("C:/Users/peree/OneDrive/Desktop")

tiff("C:/Users/peree/OneDrive/Desktop/Plot3.png", width = 10, height = 30, units = 'in', res = 300)
png("C:/Users/peree/OneDrive/Desktop/Plot4.png", width = 10, height = 30, units = 'in', res = 300)

pheatmap(as.matrix(all_marker_features_df),
         border_color = "white", 
         cellwidth = 8,
         cellheight = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = c(1:13),
         gaps_row = c(44, 54, 84, 132),
         color = colorRampPalette(brewer.pal(n = 5, name = "RdYlBu"))(100),
         main = "Included and Excluded Features per Marker")

dev.off()

#-------------------------------------------------------------------------------
marker_files = c('Heh2.txt', 'Hta2.txt')

all_marker_features_df = data.frame(matrix(ncol = 0, nrow = 246))

for (mf in marker_files) {
  included = read.delim(paste(included_features, "/", mf, sep = ""), header = FALSE)
  excluded = read.delim(paste(excluded_features, "/", mf, sep = ""), header = FALSE)
  inc_exc = c(c(rep('included', nrow(included))), c(rep('excluded', nrow(excluded))))
  
  all_features_df = data.frame(cbind(inc_exc, c(included[ , 1], excluded[ , 1])))
  colnames(all_features_df) = c("Included/Excluded", "Feature")
  all_features_df = all_features_df[naturalorder(all_features_df$Feature), ]
  
  inc_exc = all_features_df$`Included/Excluded`
  inc_exc = as.integer(inc_exc == 'included')
  
  all_marker_features_df = cbind(all_marker_features_df, inc_exc)
}

rownames(all_marker_features_df) = rbind(all_features_df$Feature)
colnames(all_marker_features_df) = c('Heh2', 'Hta2')


setwd("C:/Users/peree/OneDrive/Desktop")

png("C:/Users/peree/OneDrive/Desktop/Plot4.png", width = 10, height = 30, units = 'in', res = 300)

pheatmap(as.matrix(all_marker_features_df),
         border_color = "white", 
         cellwidth = 8,
         cellheight = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = c(1:2),
         gaps_row = c(44, 64, 94, 142),
         color = colorRampPalette(brewer.pal(n = 5, name = "RdYlBu"))(100),
         main = "Included and Excluded Features per Marker")

dev.off()
