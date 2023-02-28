# % of wildtype cells thrown out by garbage collector
# Alex Daiejavad
# December 2, 2022

# Purpose: one of Myra's remaining taks that she gave me to complete. The aim
# is to determine the percentage of wildtype cells thrown out by the garbage 
# collector.
#-------------------------------------------------------------------------------

setwd("//metax.ccbr.utoronto.ca/morphology/mpg6/mojcamu/Marker_Project/cell_cycle")
root_dir = "//metax.ccbr.utoronto.ca/morphology/mpg6/mojcamu/Marker_Project/cell_cycle"

library(dplyr)

# All marker screens
screens = list.files("//metax.ccbr.utoronto.ca/morphology/mpg6/mojcamu/Marker_Project/cell_cycle")
# df skeleton
cellcounts_df = data.frame(matrix(nrow = 0, 
                                  ncol = 4, 
                                  dimnames = list(c(), c("Marker", "Bad", "Normal", "Total"))))

for (screen in screens) { # for each screen...
  screen_dir = paste(root_dir, screen, sep = "/")
  plates = list.files(screen_dir) # get all plates
  
  marker = strsplit(screen, split = "_")[[1]][2] # get marker name
  bad = c()
  normal = c()
  
  for (plate in plates) { # for each plate...
    plate_cellcount_path = paste(screen_dir, plate, "Predictions_normal_bad_cellcounts.csv", sep = "/")
    plate_cellcounts = read.csv(file = plate_cellcount_path) # read file with normal/bad cellcounts
    
    # record cell counts
    bad = c(bad, plate_cellcounts[1, 2])
    normal = c(normal, plate_cellcounts[2, 2])
  }
  # for each screen, record the marker name, total number of normal cells, bad cells, and their sum
  cellcounts_df = rbind(cellcounts_df, data.frame(Marker = marker, 
                                                  Bad = sum(bad), 
                                                  Normal = sum(normal), 
                                                  Total = sum(sum(bad), sum(normal))))
}

# save results in case I mess up the original df
save_df = cellcounts_df
# add total cell numbers taken from CP feature files
cellcounts_df = cbind(cellcounts_df, CP_Cell_Counts = c(6081015,
                                                        3526704,
                                                        5512994,
                                                        16351241,
                                                        6064198,
                                                        3365404,
                                                        5790962,
                                                        5875768,
                                                        9230089,
                                                        12951472,
                                                        2454811,
                                                        5368374,
                                                        11175048,
                                                        4400124,
                                                        27959666))

# Source of confusion #1: the 'total's taken from sum of normal and bad is larger
# then the 'total' taken from taking the number of rows from the feature files. Why?
## ANSWER: there are actually 3 classifications by the garbage collector: "normal", 
## "bad", and "no_normalbad". It appears that no_normalbad cells were not included in
## the cellcount files I used for creating the df (are no_normalbad cells kept or
## discarded?)

# Source of confusion #2: the objective is to find the number of WILDTYPEs tossed out;
# but here we're looking at the TOTAL number of cells tossed out. This is prior to OD
# step so the data hasn't been separated into WT and mutant yet.

## I need to find the number of wildtype cells for each plate somewhere and divide that 
## by the total number of cells
## OR I need to find the total number of non-wildtype cells, subtract that from total
## number of cells, and divide by total number of cells

#-------------------------------------------------------------------------------
start.time = Sys.time()

root_dir = "C:/Users/peree/OneDrive/Desktop/cell_cycle"
screens = list.files("C:/Users/peree/OneDrive/Desktop/cell_cycle")
# df skeleton
wt_cellcounts_df = data.frame(matrix(nrow = 0,
                                     ncol = 12,
                                     dimnames = list(c(), c("Marker",
                                                            "Plate",
                                                            "Total Num. Perimeter Cells",
                                                            "Total Num. HIS3 Cells",
                                                            "Bad Cells (Perimeter)",
                                                            "Bad Cells (HIS3)",
                                                            "No_NormalBad Cells (Perimeter)",
                                                            "No_NormalBad Cells (HIS3)",
                                                            "% Bad Cells (Perimeter)",
                                                            "% Bad Cells (HIS3)",
                                                            "% Bad + No_NormalBad Cells (Perimeter)",
                                                            "% Bad + No_NormalBad Cells (HIS3)"))))
# 76 wt coordinates
for (screen in screens) {
  screen_dir = paste(root_dir, screen, sep = "/")
  plates = list.files(screen_dir)
  
  marker = strsplit(screen, split = "_")[[1]][2] # get marker name
  total_wt = c()
  bad_wt = c()
  no_normalbad_wt = c()
  
  for (plate in plates) {
    classifications_path = paste(screen_dir, plate, "Predictions_normal_bad.csv", sep = "/")
    classes = read.csv(file = classifications_path)
    classes_his3_only = dplyr::filter(classes, classes$Name == "HIS3")
    
    wt_top = dplyr::filter(classes, classes$Row == 1)
    wt_bottom = dplyr::filter(classes, classes$Row == 16)
    
    wt_left = dplyr::filter(classes, classes$Column == 1)
    wt_left = dplyr::filter(wt_left, wt_left$Row != 1)
    wt_left = dplyr::filter(wt_left, wt_left$Row != 16)
    
    wt_right = dplyr::filter(classes, classes$Column == 24)
    wt_right = dplyr::filter(wt_right, wt_right$Row != 1)
    wt_right = dplyr::filter(wt_right, wt_right$Row != 16)
    
    all_wt = data.frame(matrix(nrow = 0, ncol = 13, dimnames = list(c(), colnames(wt_top))))
    all_wt = rbind(all_wt, wt_top, wt_bottom, wt_left, wt_right)
    
    num_perimeter = nrow(all_wt)
    num_his3 = nrow(classes_his3_only)
    
    num_bad_perimeter = nrow(dplyr::filter(all_wt, all_wt$Prediction == 'bad'))
    num_bad_his3 = nrow(dplyr::filter(classes_his3_only, classes_his3_only$Prediction == 'bad'))
    
    num_no_normalbad_perimeter = nrow(dplyr::filter(all_wt, all_wt$Prediction == 'no_normalbad'))
    num_no_normalbad_his3 = nrow(dplyr::filter(classes_his3_only, classes_his3_only$Prediction == 'no_normalbad'))
    
    wt_cellcounts_df = rbind(wt_cellcounts_df, data.frame("Marker" = marker,
                                                          "Plate" = plate, 
                                                          "Total Num. Perimeter Cells" = num_perimeter,
                                                          "Total Num. HIS3 Cells" = num_his3,
                                                          "Bad Cells (Perimeter)" = num_bad_perimeter,
                                                          "Bad Cells (HIS3)" = num_bad_his3,
                                                          "No_NormalBad Cells (Perimeter)" = num_no_normalbad_perimeter,
                                                          "No_NormalBad Cells (HIS3)" = num_no_normalbad_his3,
                                                          "% Bad Cells (Perimeter)" = (num_bad_perimeter/num_perimeter)*100,
                                                          "% Bad Cells (HIS3)" = (num_bad_his3/num_his3)*100,
                                                          "% Bad + No_NormalBad Cells (Perimeter)" = ((num_no_normalbad_perimeter + num_bad_perimeter)/num_perimeter)*100,
                                                          "% Bad + No_NormalBad Cells (HIS3)" = ((num_no_normalbad_his3 + num_bad_his3)/num_his3)*100))
    
    print(plate)
    
  }
}

colnames(wt_cellcounts_df) = c("Marker",
                               "Plate",
                               "Total Num. Perimeter Cells",
                               "Total Num. HIS3 Cells",
                               "Bad Cells (Perimeter)",
                               "Bad Cells (HIS3)",
                               "No_NormalBad Cells (Perimeter)",
                               "No_NormalBad Cells (HIS3)",
                               "% Bad Cells (Perimeter)",
                               "% Bad Cells (HIS3)",
                               "% Bad + No_NormalBad Cells (Perimeter)",
                               "% Bad + No_NormalBad Cells (HIS3)")

setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/myra_tasks")
openxlsx::write.xlsx(wt_cellcounts_df, file = "wt_cellcounts_his3.xlsx")

end.time = Sys.time()
time.taken = end.time - start.time
time.taken


#------------------------------------------------------------------------------
wt_cellcounts_df = data.frame(matrix(nrow = 0,
                                     ncol = 12,
                                     dimnames = list(c(), c("Marker",
                                                            "Plate",
                                                            "Total Num. Perimeter Cells",
                                                            "Total Num. URA10 Cells",
                                                            "Bad Cells (Perimeter)",
                                                            "Bad Cells (URA10)",
                                                            "No_NormalBad Cells (Perimeter)",
                                                            "No_NormalBad Cells (URA10)",
                                                            "% Bad Cells (Perimeter)",
                                                            "% Bad Cells (URA10)",
                                                            "% Bad + No_NormalBad Cells (Perimeter)",
                                                            "% Bad + No_NormalBad Cells (URA10)"))))


screens = c("Screens_Sac6", "Screens_Snf7", "Screens_Vph1")
for (screen in screens) {
  screen_dir = paste(root_dir, screen, sep = "/")
  plates = list.files(screen_dir)
  
  marker = strsplit(screen, split = "_")[[1]][2] # get marker name
  total_wt = c()
  bad_wt = c()
  no_normalbad_wt = c()
  
  for (plate in plates) {
    classifications_path = paste(screen_dir, plate, "Predictions_normal_bad.csv", sep = "/")
    classes = read.csv(file = classifications_path)
    classes_ura10_only = dplyr::filter(classes, classes$Name == "URA10")
    
    wt_top = dplyr::filter(classes, classes$Row == 1)
    wt_bottom = dplyr::filter(classes, classes$Row == 16)
    
    wt_left = dplyr::filter(classes, classes$Column == 1)
    wt_left = dplyr::filter(wt_left, wt_left$Row != 1)
    wt_left = dplyr::filter(wt_left, wt_left$Row != 16)
    
    wt_right = dplyr::filter(classes, classes$Column == 24)
    wt_right = dplyr::filter(wt_right, wt_right$Row != 1)
    wt_right = dplyr::filter(wt_right, wt_right$Row != 16)
    
    all_wt = data.frame(matrix(nrow = 0, ncol = 13, dimnames = list(c(), colnames(wt_top))))
    all_wt = rbind(all_wt, wt_top, wt_bottom, wt_left, wt_right)
    
    num_perimeter = nrow(all_wt)
    num_ura10 = nrow(classes_ura10_only)
    
    num_bad_perimeter = nrow(dplyr::filter(all_wt, all_wt$Prediction == 'bad'))
    num_bad_ura10 = nrow(dplyr::filter(classes_ura10_only, classes_ura10_only$Prediction == 'bad'))
    
    num_no_normalbad_perimeter = nrow(dplyr::filter(all_wt, all_wt$Prediction == 'no_normalbad'))
    num_no_normalbad_ura10 = nrow(dplyr::filter(classes_ura10_only, classes_ura10_only$Prediction == 'no_normalbad'))
    
    wt_cellcounts_df = rbind(wt_cellcounts_df, data.frame("Marker" = marker,
                                                          "Plate" = plate, 
                                                          "Total Num. Perimeter Cells" = num_perimeter,
                                                          "Total Num. URA10 Cells" = num_ura10,
                                                          "Bad Cells (Perimeter)" = num_bad_perimeter,
                                                          "Bad Cells (URA10)" = num_bad_ura10,
                                                          "No_NormalBad Cells (Perimeter)" = num_no_normalbad_perimeter,
                                                          "No_NormalBad Cells (URA10)" = num_no_normalbad_ura10,
                                                          "% Bad Cells (Perimeter)" = (num_bad_perimeter/num_perimeter)*100,
                                                          "% Bad Cells (URA10)" = (num_bad_ura10/num_ura10)*100,
                                                          "% Bad + No_NormalBad Cells (Perimeter)" = ((num_no_normalbad_perimeter + num_bad_perimeter)/num_perimeter)*100,
                                                          "% Bad + No_NormalBad Cells (URA10)" = ((num_no_normalbad_ura10 + num_bad_ura10)/num_ura10)*100))
    
    print(plate)
    
  }
}

colnames(wt_cellcounts_df) = c("Marker",
                               "Plate",
                               "Total Num. Perimeter Cells",
                               "Total Num. URA10 Cells",
                               "Bad Cells (Perimeter)",
                               "Bad Cells (URA10)",
                               "No_NormalBad Cells (Perimeter)",
                               "No_NormalBad Cells (URA10)",
                               "% Bad Cells (Perimeter)",
                               "% Bad Cells (URA10)",
                               "% Bad + No_NormalBad Cells (Perimeter)",
                               "% Bad + No_NormalBad Cells (URA10)")

setwd("C:/Users/peree/OneDrive/Desktop/CompBio_Code/myra_tasks")
openxlsx::write.xlsx(wt_cellcounts_df, file = "wt_cellcounts_ura10.xlsx")