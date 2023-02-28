# GMT (Gene Matrix Transposed) File Writer
# Alex Daiejavad
# September 25, 2022

#------------------------------------------------------------------------------

if (!requireNamespace("ActivePathways", quietly = TRUE)) {
  install.packages("ActivePathways")
}

if (!requireNamespace("rWikiPathways", quietly = TRUE)) {
  install.packages("rWikiPathways")
}

if (!requireNamespace("rjson", quietly = TRUE)) {
  install.packages("rjson")
}

library(ActivePathways)
library(rWikiPathways)
library(rjson)

setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files")
DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files"

# ----------------------------- Load JSON Files -------------------------------
yeast_cc = rjson::fromJSON(file = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_c_20210430.json")[[1]]
yeast_bp = rjson::fromJSON(file = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_p_20210430.json")[[1]]
yeast_mf = rjson::fromJSON(file = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/go_slim_f_20210430.json")[[1]]

json_files = list(yeast_cc, yeast_bp, yeast_mf)
# ----------------------------- Write GMT Files -------------------------------
# Requirements for a GMT File:
## List of terms, where each term is a list with the character vectors:
### id: the term id
### name: the full name of the term
### genes: a character vector of genes annotated to this term

json_to_GMT_AP = function(json_file) {
  gmt_object = list()
  list_names = c()
  for (i in 1:length(json_file)) {
    id = json_file[[i]][[1]][1]
    name = json_file[[i]][[1]][2]
    genes = json_file[[i]][[2]]
    
    term_list = list(id = id, name = name, genes = genes)
    gmt_object = append(gmt_object, list(term_list))
    list_names = append(list_names, name)
    
  }
  names(gmt_object) = list_names
  return(gmt_object)
}

json_to_GMT_df = function(json_file) {
  gmt_object = data.frame(id = NA, name = NA, genes = NA)
  for (i in 1:length(json_file)) {
    id = json_file[[i]][[1]][1]
    name = json_file[[i]][[1]][2]
    genes = paste(json_file[[i]][[2]], collapse = "\t")
    
    gmt_object = rbind(gmt_object, c(id, name, genes))
  }

  return(gmt_object[2:nrow(gmt_object), ])
}

cc_gmt = json_to_GMT_df(yeast_cc)

bp_gmt = json_to_GMT_df(yeast_bp)
bp_gmt = bp_gmt[1:101, ]

mf_gmt = json_to_GMT_df(yeast_mf)

write.table(cc_gmt, file = 'yeast_GOSlim_cc.gmt', append = FALSE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bp_gmt, file = 'yeast_GOSlim_bp.gmt', append = FALSE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mf_gmt, file = 'yeast_GOSlim_mf.gmt', append = FALSE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


my_gmt = ActivePathways::read.GMT(file = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/yeast_GOSlim_bp.gmt')
ActivePathways::is.GMT(my_gmt)
