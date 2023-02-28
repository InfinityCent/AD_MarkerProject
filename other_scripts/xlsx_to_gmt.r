# GMT (Gene Matrix Transposed) File Writer
# Alex Daiejavad
# December 16, 2022

#------------------------------------------------------------------------------

if (!requireNamespace("ActivePathways", quietly = TRUE)) {
  install.packages("ActivePathways")
}

if (!requireNamespace("rWikiPathways", quietly = TRUE)) {
  install.packages("rWikiPathways")
}


library(ActivePathways)
library(openxlsx)
library(dplyr)

setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files")
DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files"

# ----------------------------- Load xlsx Files -------------------------------
yeast_g20 = openxlsx::read.xlsx("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/G20_Costanzo_AD.xlsx")

# ----------------------------- Write GMT Files -------------------------------
# Requirements for a GMT File:
## List of terms, where each term is a list with the character vectors:
### id: the term id
### name: the full name of the term
### genes: a character vector of genes annotated to this term

xlsx_to_GMT_df = function(xlsx_file_path) {
  yeast_g20 = openxlsx::read.xlsx(xlsx_file_path)
  
  gmt_object = data.frame(id = NA, name = NA, genes = NA)
  terms = unique(yeast_g20$Term.name)
  
  for (term in terms) {
    sub_df = dplyr::filter(yeast_g20, yeast_g20$Term.name == term)
    
    id = sub_df[1, 1]
    name = term
    genes = paste(sub_df$ORF, collapse = "\t")
    
    gmt_object = rbind(gmt_object, c(id, name, genes))
  }

  return(gmt_object[2:nrow(gmt_object), ])
}

gmt_df = xlsx_to_GMT_df("C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/G20_Costanzo_AD.xlsx")

write.table(gmt_df, file = 'G20_annotations.gmt', append = FALSE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


my_gmt = ActivePathways::read.GMT(file = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/mapping_files/G20_annotations.gmt')
ActivePathways::is.GMT(my_gmt)
