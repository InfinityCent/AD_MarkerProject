# Comparing unordered and ordered functional enrichment outputs
#September 26, 2022
# Alex Daiejavad

# What does this script do?
## Creates xlsx files containing marker-cluster columns put next to each other
## containing enriched terms for cc and bp between ordered and unordered 
## functional enrichment data

## The excel files don't look nice until further processing is done in excel
# (add colours, etc...)

#-----------------------------------------------------------------------------

# define directories containing ordered and unordered functional enrichment data
setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code')
ordered_dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_ordered_custom_bg/outputs_slim/'
unordered_dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/per_marker_unordered_custom_bg/outputs_slim/'

# data splits
split_AS = c('with_areashape', 'without_areashape')

# obtain files that are for cc or bp enrichment data only
aspect_files = function(files, aspect) {
  af = c()
  for (f in files) {
    a = stringr::str_split(f, "_")[[1]][5]
    a = stringr::str_split(a, "\\.")[[1]][1]
    
    if (a == aspect) {af = append(af, f)}
  }
  return(af)
}

# ensure that vector 1 and vector 2 are of the same length by adding the necessary
# number of NAs to the shorter vector
equal_length = function(vec1, vec2, max_length) {
  # one or both of these will be 0
  diff_vec1 = max_length - length(vec1)
  diff_vec2 = max_length - length(vec2)
  
  vec1 = append(vec1, rep(NA, diff_vec1))
  vec2 = append(vec2, rep(NA, diff_vec2))
  
  return(data.frame(vec1, vec2)) # return vec1 and vec2 next to each other as df
}

# given vector of ordered and unordered enrichment files for the same marker, 
# return the cluster *numbers* that are present in both ordered and unordered 
# enrichments
common_clusters = function(vec1, vec2) {
  vec1_cluster_nums = c()
  vec2_cluster_nums = c()
  
  for (f in vec1) {
    cluster = stringr::str_split(f, "_")[[1]][2]
    cluster_num = stringr::str_split(cluster, "-")[[1]][2]
    vec1_cluster_nums = append(vec1_cluster_nums, cluster_num)
  }
  
  for (f in vec2) {
    cluster = stringr::str_split(f, "_")[[1]][2]
    cluster_num = stringr::str_split(cluster, "-")[[1]][2]
    vec2_cluster_nums = append(vec2_cluster_nums, cluster_num)
  }

  # intersect() returns the common elements between two vectors
  return(intersect(vec1_cluster_nums, vec2_cluster_nums))
}

# given the vector containing cluster numbers seen in both ordered and unordered 
# enrichment data, get a vector containing the relevant files pertaining to these
# cluster numbers (run function twice -- once with all ordered files, once with 
# all unordered files)
common_files = function(file_names, common_cluster_nums) {
  vec_clusters = c()
  for (f in file_names) {
    cluster_num = stringr::str_split(stringr::str_split(f, "_")[[1]][2], "-")[[1]][2]
    if (cluster_num %in% common_cluster_nums) {vec_clusters = append(vec_clusters, f)}
  }
  return(vec_clusters)
}

# get human-readable marker-cluster identifier given file name and whether it
# contains ordered/unordered data
f_identifier = function(file_name, order_type) {
  marker = stringr::str_split(file_name, "_")[[1]][1]
  cluster = stringr::str_split(file_name, "_")[[1]][2]
  
  cluster_num = stringr::str_split(cluster, "-")[[1]][2]
  
  id = paste(marker, cluster_num, order_type, sep = "_")
  
  return(id)
}

# create the excel file combining ordered and unordered enrichment terms for
# each marker-cluster
direct_comparison = function(aspect, AS, col_length) {
  comparison_df = data.frame(misc = c(rep(NA, col_length))) # df skeleton
  
  # directories containing enrichment data
  ordered_as = paste(ordered_dir, AS, sep = "")
  unordered_as = paste(unordered_dir, AS, sep = "")
  
  # get all markers (different between wAS and woAS)
  markers = list.files(ordered_as)
  
  for (marker in markers) {
    # directories leading to enrichment data for each marker
    ordered_marker = paste(ordered_as, "/", marker, sep = "")
    unordered_marker = paste(unordered_as, "/", marker, sep = "")
    
    # get all file names
    ordered_marker_files = naturalsort::naturalsort(list.files(ordered_marker))
    unordered_marker_files = naturalsort::naturalsort(list.files(unordered_marker))
    
    # get file names containing only bp or cc data
    ordered_aspect_files = aspect_files(ordered_marker_files, aspect)
    unordered_aspect_files = aspect_files(unordered_marker_files, aspect)
    
    # get file names only containing data for clusters seen in both ordered and 
    # unordered datasets
    common_cluster_nums = common_clusters(ordered_aspect_files, unordered_aspect_files)
    common_ordered_files = common_files(ordered_aspect_files, common_cluster_nums)
    common_unordered_files = common_files(unordered_aspect_files, common_cluster_nums)
    
    for (i in 1:length(common_cluster_nums)) { # for each cluster...
      # load its ordered and unordered enrichment data
      ordered_f = paste(ordered_marker, "/", common_ordered_files[i], sep = "")
      unordered_f = paste(unordered_marker, "/", common_unordered_files[i], sep = "")
      
      # only save term_name column
      ordered_terms = sort(readxl::read_excel(path = ordered_f)$term_name)
      unordered_terms = sort(readxl::read_excel(path = unordered_f)$term_name)
      
      # create a df containing the ordered and unordered enrichments together,
      # made to have the same length using NA values
      mini_df = equal_length(ordered_terms, unordered_terms, col_length)
      col_names = c(f_identifier(common_ordered_files[i], 'ordered'), f_identifier(common_unordered_files[i], 'unordered'))
      colnames(mini_df) = col_names
      
      comparison_df = cbind(comparison_df, mini_df) # add df to master df
    }
  }
  return(comparison_df[ ,2:(ncol(comparison_df))]) # ignore first column; just placeholder
}

#----------------------gemerate tables and save them----------------------------
cc_wAS = direct_comparison('cc', 'with_areashape', 30)
cc_woAS = direct_comparison('cc', 'without_areashape', 30)

bp_wAS = direct_comparison('bp', 'with_areashape', 102)
bp_woAS = direct_comparison('bp', 'without_areashape', 102)

sheet_names_cc = list('With_AreaShape' = cc_wAS, 'Without_AreaShape' = cc_woAS)
openxlsx::write.xlsx(sheet_names_cc, file = 'fun_enrich_comparisons_cc_custom_bg.xlsx')

sheet_names_bp = list('With_AreaShape' = bp_wAS, 'Without_AreaShape' = bp_woAS)
openxlsx::write.xlsx(sheet_names_bp, file = 'fun_enrich_comparisons_bp_custom_bg.xlsx')
