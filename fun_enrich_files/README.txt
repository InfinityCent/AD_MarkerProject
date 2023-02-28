genes_list
  > contains list of all genes in the deletion and ts arrays

mapping_files
  > contains GO Slim JSON files used in functional enrichment

per_marker_unordered
  > contains subfolders with/without AS; markers
  > each marker subfolder contains combined excel sheet of all subclusters
    with all ORFs and corresponding penetrance values, as well as simple gene
    lists of ORFs found in each cluster
  > NOTE: combined_cluster files used in functional enrichment are found in
    CompBio_Code/mojca_clustering

per_marker_sorted
  > similar to per_marker but gene lists are sorted in order of greatest to
    least penetrance

plots
  > plots_proportion
    >>  plots where the x-axis is term proportion
    >> outdated, do not use
  > plots_fold_enrich
    >> plots where the x-axis is fold enrichment value of each term

plotting_scripts
  > functional_enrichment_visualization_[marker/cluster].r: generate
      functional enrichemnt plots for per_marker or combined_cluster
      enrichment outputs
  > functional_enrichment.r: the original functional enrichment script that
      uses gProfiler2
  > fun_enrich_values_[permarker/cluster]_[C/BP].r: script that combines the
      functional enrichment value of a C/BP term across different marker
      clusters or combined clusters
  > functional_enrichment_[marker/cluster].py: functional enrichment using
      GO Slim files with per_marker or combined_cluster gene lists
