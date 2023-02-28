# Cumulative Compartment Bar Chart Grapher
# Alex Daiejavad
# October 20, 2022

# What does this script do?
## Combines all cc terms in all clusters for each marker, then counts how many
## times a particular cc term showed up for each marker. Visuallizes the data
## using term vs. num. occurrences bar plots

#------------------------------------------------------------------------------
setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim")
plots_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/plots/per_marker_combined_bargraphs_custom_bg"
#setwd("/home/alex/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim")

# Load packages
library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(viridis)

# Load combined functional enrichment data (normally each marker-cluster has
# its own file but I combined them all into one file)
# columns: terms // rows: marker-clusters // cells: p-values
wAS_cc = readxl::read_xlsx('fun_enrich_per_marker_ordered_c_custom_bg_wAS.xlsx')
woAS_cc = readxl::read_xlsx('fun_enrich_per_marker_ordered_c_custom_bg_woAS.xlsx')

# get list of compartment terms (just take header of funenrich data minus marker)
compartments = c('marker', colnames(wAS_cc)[-1])

# wAS and woAS markers
wAS_markers = c('Cdc11', 'Dad2', 'Heh2', 'Hta2', 'Nop10', 'Nuf2', 'Om45', 
                'Pil1', 'Psr1', 'Rad52', 'Sec21', 'Sec7', 'Spf1', 
                'Vph1NoPosCtrl')
woAS_markers = c('Cdc11', 'Om45', 'Sec21', 'Spf1', 'Vph1NoPosCtrl')

# create df with all clusters thrown together and counts of how many times 
# each term showed up per marker
# columns: terms // rows: markers // cells: counts
count_compartments_table = function(markers, fun_enrich_table) {
  
  # create df skeleton
  compartment_counts = data.frame(matrix(nrow = 0, ncol = length(colnames(fun_enrich_table))))
  
  # get all clusters from first row of combined functional enrichment data
  clusters = dplyr::pull(fun_enrich_table[,1])
  
  for (marker in markers) {
    # the pattern ^marker is just a regex expression saying "anything that begins with
    # the string marker..."
    pattern = paste("^", marker, sep = "")
    marker_rows = grep(pattern, clusters) # get the row indices that begin with <marker>
    # create subset of the complete table that only contains specific marker information
    marker_table = fun_enrich_table[marker_rows,]
    
    counts = c(marker)
    for (i in 2:ncol(marker_table)) { # first column is just marker-cluster names
      # use vectorization to get a count of how many cells have values larger than 0 (actual p-values;
      # indicate that a term was present in this marker-cluster)
      counts = c(counts, sum(dplyr::pull(marker_table[,i]) > 0))
    }
    compartment_counts = rbind(compartment_counts, counts) # add counts data to master df
  }
  colnames(compartment_counts) = compartments
  return(compartment_counts)
}

# get colour palette, create df mapping each marker to a certain colour
# 14 markers -> 14 colours
inferno = viridis::inferno(14)
# inferno = c(inferno, inferno[length(inferno)]) # handle the two Vph1s
colours = data.frame(marker = wAS_markers, colour = inferno) 


#--------------------------------------------------------------------------------
# get counts table for wAS and woAS
wAS_cc_counts = count_compartments_table(wAS_markers, wAS_cc)
woAS_cc_counts = count_compartments_table(woAS_markers, woAS_cc)

# change Vph1NPC to Vph1
wAS_cc_counts$marker = c(wAS_markers[1:13], "Vph1")
woAS_cc_counts$marker = c(woAS_markers[1:4], "Vph1")

# get terms that have 0 value in both wAS and woAS; they can just be removed entirely
# I did this so I'd have to generate fewer colours for the colour palette
wAS_zeroes = names(wAS_cc_counts)[sapply(wAS_cc_counts, function(x) all(x == 0))]
woAS_zeroes = names(woAS_cc_counts)[sapply(woAS_cc_counts, function(x) all(x == 0))]
both_zeroes = intersect(wAS_zeroes, woAS_zeroes)

# change the counts tables so only terms with non-zero values are kept
wAS_cc_counts = wAS_cc_counts[,!(names(wAS_cc_counts) %in% both_zeroes)]
woAS_cc_counts = woAS_cc_counts[,!(names(woAS_cc_counts) %in% both_zeroes)]

# not a necessary step; just renamed row names with corresponding marker
rownames(wAS_cc_counts) = wAS_cc_counts$marker
rownames(woAS_cc_counts) = woAS_cc_counts$marker

# remove markers column
wAS_cc_counts = wAS_cc_counts[, -1]
woAS_cc_counts = woAS_cc_counts[, -1]

# essentially, flip the rows and columns (needed to do this for the plotting function below)
wAS_cc_counts_test = data.frame(t(wAS_cc_counts)) # t() flips rows and columns; rownames and colnames are preserved
wAS_cc_counts_test = cbind(terms = rownames(wAS_cc_counts_test), wAS_cc_counts_test) # add a terms column at beginning of df

woAS_cc_counts_test = data.frame(t(woAS_cc_counts))
woAS_cc_counts_test = cbind(terms = rownames(woAS_cc_counts_test), woAS_cc_counts_test)
#--------------------------------------------------------------------------------


# create actual plots
# counts_table must be:
## columns: markers // rows: terms // cells: counts
bar_chart_grapher = function(counts_table, AS) {
  for (i in 2:ncol(counts_table)) { # first column is just markers
    bar_data = counts_table[, c(1, i)] # create subset of counts_table that only contains terms + marker counts
    bar_data[, 2] = as.numeric(bar_data[, 2]) # make counts numeric
    bar_data = bar_data[order(-bar_data[, 2]), ] # order df from greatest counts to least
    bar_data = dplyr::filter(bar_data, bar_data[, 2] > 0) # remove 0s
    
    Compartments = bar_data[, 1]
    Values = bar_data[, 2]
    
    # get marker name and its corresponding colour
    m = colnames(counts_table)[i]
    c = dplyr::filter(colours, colours$marker == m)[1, 2]
    
    # set plot subtitle
    plot_sub = ""
    if (AS == 'wAS') {plot_sub = 'With AreaShape'}
    else if (AS == 'woAS') {plot_sub = 'Without AreaShape'}
    
    # create name of plot file, determine its path and directory
    plot_file = paste(m, "_", AS, sep = "")
    plots_directory = plots_dir
    plot_path = paste(plots_directory, "/", plot_file, ".png", sep = "")
    
    # prime R to create a png image with hard-coded dimensions
    png(file = plot_path, width = 948, height = 874)
    
    # create the actual plot (terms vs. counts)
    print(ggplot(bar_data, aes(Compartments, Values)) + 
            geom_bar(stat = 'identity', fill = c, colour = 'black') +
            labs(y = 'Number of Occurrences in Marker', x = 'Compartments', 
                 title = colnames(counts_table)[i], subtitle = plot_sub) + 
            scale_x_discrete(limits = bar_data[, 1]) + 
            scale_y_discrete(limits = 0:12) + 
            coord_flip() + 
            theme(axis.title = element_text(face="bold", size=20), 
                  plot.title = element_text(face="bold", hjust = 0.5, size=25), 
                  plot.subtitle = element_text(hjust = 0.5, size=20),
                  axis.text.x = element_text(size = 17),
                  axis.text.y = element_text(size = 17)))
    
    dev.off()
    
  }
}


# generate plots
bar_chart_grapher(wAS_cc_counts_test, 'wAS')
bar_chart_grapher(woAS_cc_counts_test, 'woAS')

# https://commentpicker.com/random-color-generator.php

# testing corner; figuring out colour palettes and visualizing them using pie chart
palette.pals() 
x = viridis::inferno(14)
x
pie(c(1:14), col = x)