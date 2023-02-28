# Functional enrichment values compiler
# Alex Daiejavad
# August 23, 2022

# What does this script do?
## Plots cluster vs. P-value for each BP term in ordered per_marker datasets

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

# define directory with enrichment data, data splits, and bp terms
setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim')
dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim'
splitAS = c('with_areashape', 'without_areashape')

processes = sort(c('RNA splicing', 'mRNA processing', 'biological_process', 
                   'cellular respiration', 'generation of precursor metabolites and energy', 
                   'nucleobase-containing small molecule metabolic process', 'DNA recombination', 
                   'ion transport', 'transmembrane transport', 'mitochondrial translation', 
                   'mitochondrion organization', 'other', 'protein targeting', 
                   'transcription from RNA polymerase III promoter', 'endosomal transport', 
                   'cytoskeleton organization', 'regulation of organelle organization', 
                   'regulation of translation', 'translational elongation', 'cytoplasmic translation', 
                   'nuclear transport', 'protein folding', 
                   'protein modification by small protein conjugation or removal', 
                   'proteolysis involved in cellular protein catabolic process', 'Golgi vesicle transport', 
                   'lipid metabolic process', 'nucleus organization', 'protein dephosphorylation', 
                   'lipid transport', 'protein complex biogenesis', 'chromatin organization', 
                   'cellular amino acid metabolic process', 'DNA replication', 
                   'carbohydrate metabolic process', 'histone modification', 
                   'regulation of DNA metabolic process', 'response to heat', 
                   'transcription from RNA polymerase II promoter', 'DNA repair', 
                   'cellular response to DNA damage stimulus', 'response to chemical', 
                   'response to oxidative stress', 'chromosome segregation', 'mitotic cell cycle', 
                   'organelle fission', 'regulation of cell cycle', 'protein phosphorylation', 
                   'cell wall organization or biogenesis', 'meiotic cell cycle', 'sporulation', 
                   'RNA modification', 'cell budding', 'tRNA processing', 
                   'DNA-templated transcription, elongation', 'RNA catabolic process', 
                   'nucleobase-containing compound transport', 'protein glycosylation', 'signaling', 
                   'rRNA processing', 'ribosomal large subunit biogenesis', 'conjugation', 
                   'endocytosis', 'organelle assembly', 'ribosomal small subunit biogenesis', 
                   'ribosome assembly', 'response to osmotic stress', 'organelle inheritance', 
                   'exocytosis', 'membrane fusion', 'organelle fusion', 'vesicle organization', 
                   'regulation of protein modification process', 'snoRNA processing', 
                   'translational initiation', 'response to starvation', 'cofactor metabolic process', 
                   'monocarboxylic acid metabolic process', 'vacuole organization', 'cytokinesis', 
                   'DNA-templated transcription, termination', 'protein maturation', 
                   'peptidyl-amino acid modification', 'protein acylation', 'peroxisome organization', 
                   'invasive growth in response to glucose limitation', 'pseudohyphal growth', 
                   'ribosomal subunit export from nucleus', 'protein alkylation', 'telomere organization', 
                   'cell morphogenesis', 'regulation of transport', 'DNA-templated transcription, initiation', 
                   'transcription from RNA polymerase I promoter', 'transposition', 
                   'vitamin metabolic process', 'tRNA aminoacylation for protein translation', 
                   'oligosaccharide metabolic process', 'protein lipidation', 
                   'cellular ion homeostasis', 'amino acid transport', 'carbohydrate transport', 
                   'not_yet_annotated'))

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
fun_enrich_values = function(procs, fun_enrich_df) {
  values = c()
  i = 1
  
  for (p in procs) { # for each bp term...
    if (p %in% fun_enrich_df[, 1]) {
      # if the term has a P-value, record it and increase i by 1
      values = append(values, as.numeric(fun_enrich_df[i, 2]))
      i = i + 1
    } else {
      values = append(values, NA) # otherwise, record NA
    }
  }
  return(values)
}

# create df
# columns: terms // rows: marker-clusters // cells: P-values
fun_enrich_table = function(directory, AS) {
  proc_df = data.frame()
  
  # get a list of markers (different between wAS and woAS)
  marker_dir = paste(directory, "/", AS, sep = "")
  markers = list.files(marker_dir)
  
  for (marker in markers) {
    # get a vector containing only bp enrichment files for each marker
    file_dir = paste(marker_dir, "/", marker, sep = "")
    all_files = naturalsort::naturalsort(list.files(file_dir))
    p_files = c()
    
    for (f in all_files) {
      if (stringr::str_split(f, "_")[[1]][5] == "bp.xlsx") {
        p_files = append(p_files, f)
      }
    }
    
    for (f in p_files) {
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
      fun_enrich_row = c(file_id, fun_enrich_values(processes, f_df_fun_enrich))
      proc_df = rbind(proc_df, fun_enrich_row)
    }
  }
  columns = c("Cluster", processes)
  colnames(proc_df) = columns
  
  return(proc_df)
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
    
    if (nrow(term_df) == 0) {next} # if no marker-clsuters for that term, move on
    # if over 30 marker-clusters for that term, only graph top 30
    if (nrow(term_df) > 30) {term_df = term_df[1:30, ]}
    
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
            geom_bar(stat='identity', fill='#6586CA', width = 0.5) +
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
wAS = fun_enrich_table(dir, "with_AreaShape")
woAS = fun_enrich_table(dir, "without_AreaShape")

sheet_names = list('with_AreaShape' = wAS, 'without_AreaShape' = woAS)
write.xlsx(sheet_names, file = 'fun_enrich_per_marker_ordered_p_custom_bg.xlsx')


wAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/plots/histograms_per_term/with_areashape_bp'
woAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/plots/histograms_per_term/without_areashape_bp'
fun_enrich_plots(wAS, "wAS", wAS_plots_directory)
fun_enrich_plots(woAS, "woAS", woAS_plots_directory)



