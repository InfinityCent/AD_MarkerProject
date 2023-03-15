# Apply thresholds 1-2 to new clustering data from March 2023
# Purpose: we re-did clustering again and I want to do functional enrichment on
#          our new clusters, but I have to process the input files first

# Date: March 14, 2023
# Alex Daiejavad
# -----------------------------------------------------------------------------

library(naturalsort)
library(dplyr)

input_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/clustering_profiles"
setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files")

# steps:
# 1. find strains with <10 cells in cluster
# 2. remove strains with <10% cluster penetrance
# 3. remove strains from step 1
# 4. sort gene lists from greatest to least penetrance
#------------------------------------------------------------------------------

# Step 1: find strains with <10 cells for each marker-cluster
under_ten_cells = function(marker, total_clusters, cluster_num) {
  # input marker: str
  # total_clusters: str
  # cluster_num: str
  
  # return: vector of strain IDs with less than ten cells in that cluster
  
  setwd(sprintf("%s/%s_all_PCA/normal_clusters/Clusters", input_dir, marker))
  
  cluster_file = read.csv(sprintf("normal_clusters_n%s_per_cell_cluster-%s.csv", total_clusters, cluster_num))
  num_strain_cells_tabulate = table(cluster_file$Strain.ID)
  strains_under_ten = names(num_strain_cells_tabulate[num_strain_cells_tabulate < 10])
  
  return(strains_under_ten)
}


# Step 2/3/4: apply thresholds 2/3 to each cluster, sort, save lists
ordered_gene_lists = function(marker) {
  setwd(sprintf("%s/%s_all_PCA", input_dir, marker))
  strain_profile = read.csv("Strain_profile.csv")
  
  # Apply threshold 1: remove strains with overall penetrance < 10%
  strain_profile = dplyr::filter(strain_profile, strain_profile$Penetrance > 0.1)
  
  for (cluster in 9:ncol(strain_profile)) {
    cluster_profile = strain_profile[ , c(1, 2, 4, cluster)]
    
    # Apply threshold 2: remove strains with cluster penetrance < 10%
    cluster_profile = dplyr::filter(cluster_profile, cluster_profile[ , 4] > 0.1)
    
    # Apply threshold 3: remove strains with cluster cell count < 10
    cluster_num = strsplit(colnames(cluster_profile)[4], "\\.")[[1]][2]
    total_clusters = ncol(strain_profile) - 8
    
    remove_strains = under_ten_cells(marker, total_clusters, cluster_num)
    cluster_profile = dplyr::filter(cluster_profile, !(cluster_profile$Strain.ID %in% remove_strains))
    
    # Sort cluster penetrances from greatest to least
    cluster_profile = dplyr::arrange(cluster_profile, desc(cluster_profile[ , 4]))
    
    # Export list
    file_name = sprintf("%s/%s_all_PCA/gene_lists/genelist_cluster-%s.csv", input_dir, marker, cluster_num)
    write.csv(cluster_profile, file_name, row.names = FALSE)
  }
}

ordered_gene_lists("Cdc11")
ordered_gene_lists("Dad2")
ordered_gene_lists("Heh2")
ordered_gene_lists("Nuf2")
ordered_gene_lists("Snf7")
