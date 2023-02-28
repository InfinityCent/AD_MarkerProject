Scripts are run in the following order:

1. functional_enrichment_[ordered/unordered].r  <- do actual functional enrichment analysis
2. compare_ordered_unordered_data.r  <- put ordered and unordered terms for each marker-cluster together for comparison
3. fun_enrich_values_permarker_[ordered/unordered]_[BP/CC].r  <- create table with all functional enrichment terms/p-values and per-term histograms
4. fun_enrich_cumulative_cc_bargraph.r  <- get quick overview of most enriched terms per marker (combine clusters together)
5. fun_enrich_bp_heatmap.r <- display BP results in a pretty heatmap
