# Visualizing functional enrichment outputs
# Alex Daiejavad
# August 03, 2022

# What does this script do?
## Create term vs. proportion bar plots for combined cluster data

#==========================================================

# Install and load packages
if (! requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (! requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (! requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
if (! requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
#-------------------------------------------------------------------------------

# Set location of data that needs to be plotted and where plots should be output to
setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted')
directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/combined_cluster_outputs'
plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/plots/'

# Data splits
split_CC = c('no_cellcycle', 'with_cellcycle')
split_AS = c('with_AreaShape', 'without_AreaShape')
split_genes = c('all_genes', 'ess_only', 'no_vesicle', 'noness_only')

# https://www.geeksforgeeks.org/how-to-read-a-xlsx-file-with-multiple-sheets-in-r/
# Never ended up using this but still helpful
multiplesheets = function(path_name) {
  sheets = readxl::excel_sheets(path_name)
  tibble = lapply(sheets, function(x) readxl::read_excel(path_name, sheet = x))
  data_frame = lapply(tibble, as.data.frame)
  
  names(data_frame) = sheets
  
  return(data_frame)
}

# given a filename, extract relevant data (normal/cluster, number, marker, aspect)
clean_file_name = function(file_name) {
  break_up = stringr::str_split(file_name, "-")[[1]]
  
  transformed = c()
  for (s in break_up) {
    if (startsWith(s, "normal")) {
      nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "normal")[[1]][2]
      nc_num = paste("Normal", nc_num, "of", stringr::str_split(s, "of")[[1]][2])
      transformed = c(transformed, nc_num)
    }
    if (startsWith(s, "Cluster")) {
      nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "Cluster")[[1]][2]
      nc_num = paste("Cluster", nc_num, "of", stringr::str_split(s, "of")[[1]][2])
      transformed = c(transformed, nc_num)
    }
    if (s == "go_slim_p.xlsx") {
      transformed = c(transformed, "Biological Process")
    }
    if (s == "go_slim_c.xlsx") {
      transformed = c(transformed, "Cellular Compartment")
    }
  }
  return(transformed)
}

# simple function that just makes given CC/AS/gene_dir human-readable
plot_info = function(CC, AS, gene_dir) {
  transformed = c()
  if (CC == "with_cellcycle") {
    transformed = c(transformed, "With Cell Cycle")
  }
  if (CC == "no_cellcycle") {
    transformed = c(transformed, "Without Cell Cycle")
  }
  
  if (AS == "with_AreaShape") {
    transformed = c(transformed, "With AreaShape")
  }
  if (AS == "without_AreaShape") {
    transformed = c(transformed, "Without AreaShape")
  }
  if (gene_dir == "all_genes") {
    transformed = c(transformed, "All Genes")
  }
  if (gene_dir == "ess_only") {
    transformed = c(transformed, "Essential Genes")
  }
  if (gene_dir == "noness_only") {
    transformed = c(transformed, "Nonessential Genes")
  }
  if (gene_dir == "no_vesicle") {
    transformed = c(transformed, "No Vesicle Genes")
  }
  
  clean_name = paste(transformed, collapse = " | ")
  
  return(clean_name)
}

for (CC in split_CC) {
  for (AS in split_AS) {
    for (gene_dir in split_genes) {
      # generate necessary file path to reach each terminal directory, load all files
      files_path = paste(directory, "/", CC, "/", AS, "/", gene_dir, sep = "")
      files = list.files(files_path)
      
      # get sheet path, load it and set filter so only terms with P-value < 0.01 are used
      # also create new column with proportion values (never ended up using this in 
      # the future -- not the best approach)
      for (f in files) {
        sheet_path = paste(files_path, "/", f, sep = "")
        sheet_df = readxl::read_excel(path = sheet_path)
        Proportion = sheet_df$'Hits in term' / sheet_df$'Term size'
        sheet_df = cbind(sheet_df, Proportion)
        sheet_df = dplyr::filter(sheet_df, sheet_df$'P-value' < 0.01)
        
        # some sheets have no terms with P-value < 0.01 and thus they're empty; 
        # ignore them and move on
        if (nrow(sheet_df) == 0) {
          next
        }
        
        # Extract just the Ontology and Fold enrichment columns for analysis,
        # break up terms that are too long on a single line
        # factor() used as a way to prevent ggplot from ordering terms alphabetically
        sheet_df_subset = data.frame(sheet_df$'Ontology', sheet_df$'Proportion')
        colnames(sheet_df_subset) = c('Term Name', 'Proportion')
        sheet_df_subset$'Term Name' = factor(stringr::str_wrap(sheet_df_subset$'Term Name', 30))
        
#------ just stuff that needs to be included in each generated plot------------
        cluster_num = clean_file_name(f)
        
        plot_title = plot_info(CC, AS, gene_dir)
        
        plot_subtitle = cluster_num[1]
        
        plot_color = c()
        if (cluster_num[2] == "Biological Process") {plot_color = '#6586CA'}
        if (cluster_num[2] == "Cellular Compartment") {plot_color = '#FFB516'}
        
        size = c()
        if (cluster_num[2] == "Biological Process") {size = NULL}
        if (cluster_num[2] == "Cellular Compartment") {size = 13}
# ------------------------------------------------------------------------------
        # create plot name and its path
        plot_file = paste(CC, "-", AS, "-", gene_dir, "-", stringr::str_split(f, "\\.")[[1]][1], ".png", sep = "")
        plot_path = paste(plots_directory, plot_file, sep = "")
        
        # prime R to create a png image with hard-coded dimensions
        png(file = plot_path, width = 948, height = 874)
        
        # create the actual barplot (term vs proportion for each datasplit)
        print(ggplot(sheet_df_subset, aes(sheet_df_subset[ , 2], sheet_df_subset[ , 1])) + 
                geom_bar(stat='identity', fill=plot_color, width = 0.5) +
                labs(x='Proportion', y=cluster_num[2], title=plot_title, subtitle=plot_subtitle) + 
                theme(axis.text = element_text(size=10), 
                      axis.title = element_text(face="bold", size=13), 
                      plot.title = element_text(face="bold", hjust = 0.5, size=15), 
                      plot.subtitle = element_text(hjust = 0.5),
                      axis.text.x = element_text(size = size),
                      axis.text.y = element_text(size = size)))
        
        
        dev.off()
      }
    }
  }
}
