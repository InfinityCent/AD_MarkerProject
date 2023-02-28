# Functional enrichment values compiler
# Alex Daiejavad
# August 23, 2022

# What does this script do?
## Plots cluster vs. P-value for each CC term in ordered per_marker datasets
# -----------------------------------------------------------------------------

# Install and load packages
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(stringr)
library(naturalsort)
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)

# define directory with enrichment data, data splits, and cc terms
setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim')
dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim'
splitAS = c('with_areashape', 'without_areashape')

compartments = c("cell cortex", "cell wall", "cellular bud", "cellular_component", "chromosome", 
                 "cytoplasm", "cytoplasmic vesicle", "cytoskeleton", 
                 "endomembrane system", "endoplasmic reticulum", "extracellular region", "Golgi apparatus", 
                 "membrane", "microtubule organizing center", 
                 "mitochondrial envelope", "mitochondrion", "nucleolus", "nucleus", "other", 
                 "peroxisome", "plasma membrane", "ribosome", "site of polarized growth",
                 "vacuole")

# get human-readable marker-cluster identifier
f_identifier = function(file_name) {
  marker = stringr::str_split(file_name, "_")[[1]][1]
  cluster = stringr::str_split(file_name, "_")[[1]][2]
  
  cluster_num = stringr::str_split(cluster, "-")[[1]][2]
  cluster_num = paste("Normal ", cluster_num, sep = "")
  
  id = paste(marker, ": ", cluster_num, sep = "")
  
  return(id)
}

# create a vector containing P-values of all terms for each marker-cluster
fun_enrich_values = function(comps, fun_enrich_df) {
  values = c()
  i = 1
  
  for (c in comps) { # for each cc term...
    if (c %in% fun_enrich_df[, 1]) { 
      # if the term has a P-value, record it and increase i by 1
      values = append(values, as.numeric(fun_enrich_df[i, 2]))
      i = i + 1
    } else { # otherwise, record NA
      values = append(values, NA)
    }
  }
  return(values)
}

# create df
# columns: terms // rows: marker-clusters // cells: P-values
fun_enrich_table = function(directory, AS) {
  comp_df = data.frame()
  
  # get a list of markers (different between wAS and woAS)
  marker_dir = paste(directory, "/", AS, sep = "")
  markers = list.files(marker_dir)
  
  for (marker in markers) {
    # get a vector containing only cc enrichment files for each marker
    file_dir = paste(marker_dir, "/", marker, sep = "")
    all_files = naturalsort::naturalsort(list.files(file_dir))
    c_files = c()
    
    for (f in all_files) {
      if (stringr::str_split(f, "_")[[1]][5] == "cc.xlsx") {
        c_files = append(c_files, f)
      }
    }
    
    for (f in c_files) {
      # load each enrichment file, order alphabetically by term, filter by P-value < 0.01
      f_path = paste(file_dir, "/", f, sep = "")
      f_df = readxl::read_excel(path = f_path)
      f_df = f_df[order(f_df$'term_name'), ]
      f_df = dplyr::filter(f_df, f_df$'p_value' < 0.01)
      
      if (nrow(f_df) == 0) { # if there are no enrichments with P-value < 0.01, move on
        next
      }
      
      # get human-readable marker-cluster identifier, create subset of enrichment
      # data containing term names and corresponding P-values only
      file_id = f_identifier(f)
      f_df_fun_enrich = data.frame(f_df$'term_name', f_df$'p_value')
      colnames(f_df_fun_enrich) = c('Term', 'P Value')
      
      # get vector containing P-values (or NAs) for all terms, add it to master df
      fun_enrich_row = c(file_id, fun_enrich_values(compartments, f_df_fun_enrich))
      comp_df = rbind(comp_df, fun_enrich_row)
    }
  }
  columns = c("Cluster", compartments)
  colnames(comp_df) = columns
  
  return(comp_df)
}

# create marker-clusters vs. P-value plots for each term
fun_enrich_plots = function(fe_table, AS, plots_path) {
  # get subset of fe-table without marker-cluster column
  fe_table_noclusters = fe_table[ , 2:ncol(fe_table)]
  terms = colnames(fe_table_noclusters) # get all terms
  
  for (i in 1:ncol(fe_table_noclusters)) { # for each term...
    term = terms[i] # get actual term name
    # create a subset of fe_table containing only marker-clusters and corresponding 
    # P-values for each marker-cluster
    term_df = data.frame(fe_table[, 1], fe_table_noclusters[, i])
    colnames(term_df) = c('Cluster', term)
    
    # order P-values from smallest to largest, ensure they're actually numeric
    term_df = term_df[order(as.numeric(term_df[,2]), na.last = NA), ]
    term_df[ , 2] = as.numeric(term_df[ , 2])
    
    if (nrow(term_df) == 0) {next} # if no marker-clusters for that term, move on
    # if over 30 marker-clusters for that term, only graph top 30
    if (nrow(term_df) > 20) {term_df = term_df[1:20, ]}
    
    # get plot file name and path
    plot_file = paste(term, "_", AS, ".png", sep = "")
    plot_path = paste(plots_path, "/", plot_file, sep = "")
    
    # get plot subtitle
    plot_subtitle = ''
    if (AS == 'wAS') {plot_subtitle = 'With AreaShape | Ordered Gene Set'}
    if (AS == 'woAS') {plot_subtitle = 'Without AreaShape | Ordered Gene Set'}
    
    # prime R to create png image with hard-coded dimensions
    png(file = plot_path, width = 948, height = 874)
    
    # create marker-cluster vs. P-values plot for each term
    print(ggplot(term_df, aes(x = term_df[ , 2], y = fct_rev(fct_inorder(term_df[ , 1])))) + 
            geom_bar(stat='identity', fill='#FFB516', width = 0.5) +
            labs(x='P Value', y='Marker: Cluster', title=term, 
                 subtitle=plot_subtitle) + 
            theme(axis.text = element_text(size=10), 
                  axis.title = element_text(face="bold", size=13), 
                  plot.title = element_text(face="bold", hjust = 0.5, size=15), 
                  plot.subtitle = element_text(hjust = 0.5)))
    
    
    dev.off()
  }
}

#-------------------create tables and plots------------------------------------
wAS = fun_enrich_table(dir, "with_areashape")
woAS = fun_enrich_table(dir, "without_areashape")

sheet_names = list('with_AreaShape' = wAS, 'without_AreaShape' = woAS)
write.xlsx(sheet_names, file = 'fun_enrich_per_marker_ordered_c_custom_bg.xlsx')


wAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/plots/histograms_per_term/with_areashape_cc'
woAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/plots/histograms_per_term/without_areashape_cc'
fun_enrich_plots(wAS, "wAS", wAS_plots_directory)
fun_enrich_plots(woAS, "woAS", woAS_plots_directory)



