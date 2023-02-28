# Scraping SGD
# Purpose: accesses SGD HTML scripts to extract short gene descriptions

# Alex Daiejavad
# January 5, 2023
#------------------------------------------------------------------------------

library(readxl)
library(writexl)
library(rvest)

setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code")


all_orfs = readxl::read_xlsx("ORFs.xlsx")[[1]]

descriptions = c()
for (orf in all_orfs) {
  site_url = paste("https://www.yeastgenome.org/locus", orf, "go", sep = "/")
  site = rvest::read_html(site_url)
  go_summary = rvest::html_text(rvest::html_nodes(site, "dd"))[[1]]
  descriptions = c(descriptions, go_summary)
}

annotated_orfs = data.frame(ORF = all_orfs, "GO Summary" = descriptions)

writexl::write_xlsx(annotated_orfs, "annotated_orfs.xlsx")
