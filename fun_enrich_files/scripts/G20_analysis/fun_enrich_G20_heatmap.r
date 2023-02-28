# G20 Heatmapper
# Alex Daiejavad
# December 16, 2022
#------------------------------------------------------------------------------
setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_G20")
# setwd("/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim")

library(readxl)
library(openxlsx)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

wAS_G20 = readxl::read_xlsx('fun_enrich_per_marker_ordered_g20_custom_bg_wAS.xlsx')
woAS_G20 = readxl::read_xlsx('fun_enrich_per_marker_ordered_g20_custom_bg_woAS.xlsx')

# remove Vphh1WPC rows
wAS_G20 = wAS_G20[1:128, ]
woAS_G20 = woAS_G20[1:41, ]

# rename Vph1NPC to just Vph1
wAS_G20_clusters = wAS_G20$Cluster
woAS_G20_clusters = woAS_G20$Cluster

wAS_G20_clusters = wAS_G20_clusters[1:122]
woAS_G20_clusters = woAS_G20_clusters[1:33]

wAS_G20_clusters = c(wAS_G20_clusters, "Vph1: Normal 3", "Vph1: Normal 4", 
                    "Vph1: Normal 5", "Vph1: Normal 7", "Vph1: Normal 8",
                    "Vph1: Normal 9")
woAS_G20_clusters = c(woAS_G20_clusters, "Vph1: Normal 1", "Vph1: Normal 2", 
                      "Vph1: Normal 4", "Vph1: Normal 5", "Vph1: Normal 7",  
                      "Vph1: Normal 8", "Vph1: Normal 9", "Vph1: Normal 10")

wAS_G20$Cluster = wAS_G20_clusters
woAS_G20$Cluster = woAS_G20_clusters

terms = c("G1/S and G2/M cell cycle progression/meiosis",
          "polarity/morphogenesis/cytokenesis/endocytosis/exocytosis",
          
          "chromatin/transcription",
          "chromosome segregation/kinetochore/spindle/microtubule",
          "DNA replication/repair/HR/cohesion",
          "RNA processing",
          
          "ribosome/translation and tRNA processing",
          "amino acid biosynth&transport/nitrogen utilization",
          "proteosome/protein degradation/ER-dependent protein folding/degradation",
          
          #"drug/ion transport/plasma membrane pumps",
          "ER<->Golgi traffic/ER translocation",
          "Golgi<->endosome<->vacuole<->plasma membrane/intra-golgi transport/vatpase",
          "lipid/sterol/fatty acid biosynth & transport",
          "nuclear-cytoplasic transport",
          
          #"autophagy/CVT",
          "glycosylation/GPI anchor/cell wall",
          #"metabolism/mitochondria",
          #"peroxisome",
          "signaling/stress response"
          
          #"multi-functional could not assign to specific terms"
          )


# get terms that have 0 value in both wAS and woAS; they can just be removed entirely
wAS_zeroes = names(wAS_G20)[sapply(wAS_G20, function(x) all(x == 0))]
woAS_zeroes = names(woAS_G20)[sapply(woAS_G20, function(x) all(x == 0))]
both_zeroes = intersect(wAS_zeroes, woAS_zeroes)

# change the counts tables so only terms with non-zero values are kept
wAS_G20 = wAS_G20[,!(names(wAS_G20) %in% both_zeroes)]
woAS_G20 = woAS_G20[,!(names(woAS_G20) %in% both_zeroes)]

# turn all placeholder 0s to NA
wAS_G20[wAS_G20 == 0] = NA
woAS_G20[woAS_G20 == 0] = NA

# reorder columns based on manually grouped processes
wAS_G20 = wAS_G20[, c("Cluster", terms)]
woAS_G20 = woAS_G20[, c("Cluster", terms)]

# remove markers column
wAS_G20 = wAS_G20[, -1]
woAS_G20 = woAS_G20[, -1]

# essentially, flip the rows and columns (needed to do this for the plotting function below)
wAS_G20_test = data.frame(t(wAS_G20)) # t() flips rows and columns; rownames and colnames are preserved
wAS_G20_test = cbind(terms = rownames(wAS_G20_test), wAS_G20_test) # add a terms column at beginning of df

woAS_G20_test = data.frame(t(woAS_G20))
woAS_G20_test = cbind(terms = rownames(woAS_G20_test), woAS_G20_test)

colnames(wAS_G20_test) = c('Terms', wAS_G20_clusters)
colnames(woAS_G20_test) = c('Terms', woAS_G20_clusters)

#------------------------------------------------------------------------------

# Create the heatmaps
wAS_matrix = as.matrix(wAS_G20_test[,-1])
woAS_matrix = as.matrix(woAS_G20_test[,-1])

legend_divs = c(0, 0.002, 0.004, 0.006, 0.008, 0.01)

col_divs_wAS = c(11, 22, 32, 42, 53, 61, 71, 81, 88, 98, 104, 112, 122)
col_divs_woAS = c(11, 21, 28, 33)

row_divs = c(2, 6, 9, 13)


pheatmap(wAS_matrix,
         border_color = "white", 
         cellwidth = 7.5,
         cellheight = 11,
         legend_breaks = legend_divs,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = col_divs_wAS,
         gaps_row = row_divs,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
         main = "With AreaShape | Group 19 Annotations")


pheatmap(woAS_matrix,
         border_color = "white", 
         cellwidth = 15,
         cellheight = 15,
         legend_breaks = legend_divs,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = col_divs_woAS,
         gaps_row = row_divs,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
         main = "Without AreaShape | Group 19 Annotations")

# plot_zoom_png?width=1920&height=1017