# FUNCTIONAL ENCRICHMENT SCRIPT
# JULY 18, 2022
# ALEX DAIEJAVAD

# What does this script do?
## Do functional enrichment for cc and bp using GO Slims and gProfiler2 for
## unordered gene lists per marker-cluster

# -----------------------------------------------------------------------------

# Install and load packages
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

if (! requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (! requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}

if (! requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (! requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (! requireNamespace("naturalsort", quietly = TRUE)) {
  install.packages("naturalsort")
}

library(openxlsx)
library(readxl)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(stringr)
library(naturalsort)
#------------------------------------------------------------------------------

# Define directories and data splits
# Windows
setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_unordered_custom_bg")
INPUT_DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/genes_list"
OUTPUT_DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_unordered_custom_bg"
BACKGROUND_DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/genes_list/background_genes"

# Linux
setwd("/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_unordered_custom_bg")
INPUT_DIRECTORY = "/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/genes_list"
OUTPUT_DIRECTORY = "/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_unordered_custom_bg"
BACKGROUND_DIRECTORY = "/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/genes_list/background_genes"

split_AS = c('with_areashape', 'without_areashape')

# Create standardized output xlsx file name
file_name = function(marker_file, aspect) {
  split_name = stringr::str_split(marker_file, "_")[[1]]
  file_name = paste(split_name[1], split_name[3], 'unordered_enrichment', aspect, sep = "_")
  file_name = paste(file_name, '.xlsx', sep = "")
  
  return(file_name)
}

# Upload GMT files containing GO Slim mappings
#cc_gmt = gprofiler2::upload_GMT_file('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_c_20210430.gmt')
#bp_gmt = gprofiler2::upload_GMT_file('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_p_20210430.gmt')
#cc_gmt = gprofiler2::upload_GMT_file('/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_c_20210430.gmt') # gp__nb33_fXYE_RDA
#bp_gmt = gprofiler2::upload_GMT_file('/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_p_20210430.gmt') # gp__F3En_WJYo_Rao

# Create df with paths to background genes for custom_bg
m = c('Cdc11', 'Dad2', 'Heh2', 'Hta2', 'Nop10', 'Nuf2', 'Om45', 'Pil1', 'Psr1', 
      'Rad52', 'Sec7', 'Sec21', 'Spf1', 'Vph1')
f = list.files(BACKGROUND_DIRECTORY)
f = paste(INPUT_DIRECTORY, "/background_genes/", f, sep = "")
background_files_df = data.frame(marker = m, filepath = f)

# Determine correct file path for custom_bg in gprofiler2::gost
bg_path = function(m, bg_df) {
  path = c()
  if (m == "Vph1NoPosCtrl" | m == "Vph1WithPosCtrl") {
    path = dplyr::filter(bg_df, bg_df$marker == "Vph1")[1, 2]
  }
  else {
    path = dplyr::filter(bg_df, bg_df$marker == m)[1, 2]
  }
  
  return(path)
}

# Do functional enrichment
# go_standard is the gmt file, aspect refers to cc or bp enrichment
gost_enrich = function(go_standard, aspect) {
  for (AS in split_AS) {
    AS_path = paste(INPUT_DIRECTORY, "/", AS, sep = "") # go into AS-specific folder and load markers
    markers = list.files(AS_path) # different markers for wAS and woAS
    setwd(paste(OUTPUT_DIRECTORY, "/outputs_slim/", AS, sep = "")) # go to output folder
    
    for (marker in markers) {
      marker_path = paste(AS_path, "/", marker, sep = "") # go into marker file with the gene lists, load files
      marker_files = naturalsort::naturalsort(list.files(marker_path))
      bg_genes = read.csv(bg_path(marker, background_files_df)) # get path of bg file and read it
      bg_genes = bg_genes[ , 1] # turn it into a vector
      
      for (marker_file in marker_files) {
        xlsx_file = file_name(marker_file, aspect) # create output file name
        file_path = paste(marker_path, "/", marker_file, sep = "") # get full gene list file path, load it
        marker_data = read.csv(file_path)
        
        # do actual enrichment (unordered)
        # gost output is a list, only save the actual results (the other
        # components are unnecessary)
        enrichment = gprofiler2::gost(query = marker_data$ORF,
                                      organism = go_standard,
                                      ordered_query = FALSE,
                                      significant = TRUE,
                                      user_threshold = 0.01,
                                      correction_method = "fdr",
                                      sources = c(),
                                      custom_bg = bg_genes)$result
        
        # if there are no enriched terms, just move on (otherwise it throws out
        # an error and interrupts the loop)
        if (is.null(enrichment)) {next}
        
        # save enrichments
        write.xlsx(enrichment, file = xlsx_file)
      }
    }
  }
}

cc_gmt = "gp__nb33_fXYE_RDA"
bp_gmt = "gp__F3En_WJYo_Rao"

# separate enrichments for cc and bp
gost_enrich(cc_gmt, 'cc')
gost_enrich(bp_gmt, 'bp')
