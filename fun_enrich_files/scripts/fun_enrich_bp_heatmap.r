# BP Heatmapper
# Alex Daiejavad
# October 26, 2022
#------------------------------------------------------------------------------
setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim")
# setwd("/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim")

library(readxl)
library(openxlsx)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

wAS_bp = readxl::read_xlsx('fun_enrich_per_marker_ordered_p_custom_bg_wAS.xlsx')
woAS_bp = readxl::read_xlsx('fun_enrich_per_marker_ordered_p_custom_bg_woAS.xlsx')

# remove Vphh1WPC rows
wAS_bp = wAS_bp[1:129, ] 
woAS_bp = woAS_bp[1:43, ]

# rename Vph1NPC to just Vph1
wAS_bp_clusters = wAS_bp$Cluster
woAS_bp_clusters = woAS_bp$Cluster

wAS_bp_clusters = wAS_bp_clusters[1:124]
woAS_bp_clusters = woAS_bp_clusters[1:36]

wAS_bp_clusters = c(wAS_bp_clusters, "Vph1: Normal 3", "Vph1: Normal 4", 
                    "Vph1: Normal 5", "Vph1: Normal 7", "Vph1: Normal 9")
woAS_bp_clusters = c(woAS_bp_clusters, "Vph1: Normal 1",  "Vph1: Normal 4",
                     "Vph1: Normal 5", "Vph1: Normal 7",  "Vph1: Normal 8",  
                     "Vph1: Normal 9", "Vph1: Normal 10")

wAS_bp$Cluster = wAS_bp_clusters
woAS_bp$Cluster = woAS_bp_clusters

# processes = colnames(wAS_bp)

processes = c("cell budding",
              "cell morphogenesis",
              "cell wall organization or biogenesis",
              "chromosome segregation",
              "conjugation",
              "cytokinesis",
              "cytoskeleton organization",
              "membrane fusion",
              "meiotic cell cycle",
              "mitotic cell cycle",
              "regulation of cell cycle",
              "sporulation",
              
              
              "cellular response to DNA damage stimulus",
              "chromatin organization",
              "DNA recombination",
              "DNA repair",
              "DNA replication",
              "DNA-templated transcription, elongation",
              "DNA-templated transcription, initiation",
              "DNA-templated transcription, termination",
              "histone modification",
              "nucleus organization",
              "telomere organization",
              "transposition",
              
              
              "mRNA processing",
              "RNA catabolic process",
              "RNA splicing",
              "rRNA processing",
              "snoRNA processing",
              "tRNA aminoacylation for protein translation",
              
              
              "transcription from RNA polymerase I promoter",
              "transcription from RNA polymerase II promoter",
              "transcription from RNA polymerase III promoter",
              "cytoplasmic translation",
              "translational elongation",
              "translational initiation",
              
              
              "peptidyl-amino acid modification",
              "protein acylation",
              "protein alkylation",
              "protein complex biogenesis",
              "protein dephosphorylation",
              "protein glycosylation",
              "protein lipidation",
              "protein maturation",
              "protein modification by small protein conjugation or removal",
              "protein phosphorylation",
              "protein targeting",
              "proteolysis involved in cellular protein catabolic process",
              "ribosomal large subunit biogenesis",
              "regulation of protein modification process",
              
              
              "endocytosis",
              "endosomal transport",
              "exocytosis",
              "Golgi vesicle transport",
              "nuclear transport",
              "nucleobase-containing compound transport",
              "regulation of transport",
              "ribosomal subunit export from nucleus",
              
              
              "cofactor metabolic process",
              "invasive growth in response to glucose limitation",
              "lipid metabolic process",
              "oligosaccharide metabolic process",
              "regulation of DNA metabolic process",
              
              
              
              "mitochondrion organization",
              "organelle assembly",
              "organelle fission",
              "organelle fusion",
              "organelle inheritance",
              "peroxisome organization",
              "regulation of organelle organization",
              "vacuole organization",
              "vesicle organization",
              
              
              "response to osmotic stress",
              "signaling")


# get terms that have 0 value in both wAS and woAS; they can just be removed entirely
wAS_zeroes = names(wAS_bp)[sapply(wAS_bp, function(x) all(x == 0))]
woAS_zeroes = names(woAS_bp)[sapply(woAS_bp, function(x) all(x == 0))]
both_zeroes = intersect(wAS_zeroes, woAS_zeroes)

# change the counts tables so only terms with non-zero values are kept
wAS_bp = wAS_bp[,!(names(wAS_bp) %in% both_zeroes)]
woAS_bp = woAS_bp[,!(names(woAS_bp) %in% both_zeroes)]

# turn all placeholder 0s to NA
wAS_bp[wAS_bp == 0] = NA
woAS_bp[woAS_bp == 0] = NA

# reorder columns based on manually grouped processes
wAS_bp = wAS_bp[, c("Cluster", processes)]
woAS_bp = woAS_bp[, c("Cluster", processes)]

# remove markers column
wAS_bp = wAS_bp[, -1]
woAS_bp = woAS_bp[, -1]

# essentially, flip the rows and columns (needed to do this for the plotting function below)
wAS_bp_test = data.frame(t(wAS_bp)) # t() flips rows and columns; rownames and colnames are preserved
wAS_bp_test = cbind(terms = rownames(wAS_bp_test), wAS_bp_test) # add a terms column at beginning of df

woAS_bp_test = data.frame(t(woAS_bp))
woAS_bp_test = cbind(terms = rownames(woAS_bp_test), woAS_bp_test)

colnames(wAS_bp_test) = c('Terms', wAS_bp_clusters)
colnames(woAS_bp_test) = c('Terms', woAS_bp_clusters)

#------------------------------------------------------------------------------

# Create the heatmaps
wAS_matrix = as.matrix(wAS_bp_test[,-1])
woAS_matrix = as.matrix(woAS_bp_test[,-1])

legend_divs = c(0, 0.002, 0.004, 0.006, 0.008, 0.01)

col_divs_wAS = c(11, 22, 33, 44, 54, 60, 70, 80, 88, 99, 107, 115, 124)
col_divs_woAS = c(10, 22, 30, 36)

row_divs = c(12, 24, 30, 36, 50, 58, 63, 72)


pheatmap(wAS_matrix,
         border_color = "white", 
         cellwidth = 8,
         legend_breaks = legend_divs,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = col_divs_wAS,
         gaps_row = row_divs,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
         main = "With AreaShape | GO Slim Biological Processes")


pheatmap(woAS_matrix,
         border_color = "white", 
         legend_breaks = legend_divs,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         gaps_col = col_divs_woAS,
         gaps_row = row_divs,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
         main = "Without AreaShape | GO Slim Biological Processes")

# plot_zoom_png?width=1920&height=1017