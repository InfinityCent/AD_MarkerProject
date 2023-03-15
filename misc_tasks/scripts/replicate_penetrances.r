# Get number of strains with outlier cells in only 1-2 reps
# Done to get justification for threshold 0

# Alex Daiejavad
# March 13, 2023
#-------------------------------------------------------------------------------
library(naturalsort)
library(dplyr)

input_dir = "C:/Users/peree/OneDrive/Desktop/CompBio_Code/misc_tasks/OD penetrance comparison - plate level"
markers = list.files(input_dir)

master_df = data.frame(matrix(nrow = 0, ncol = 6))

for (marker in markers) {
  setwd(paste(input_dir, marker, sep = "/"))
  plate_files = naturalsort::naturalsort(list.files())
  
  r1r2 = plate_files[1:14]
  r1r3 = plate_files[15:28]
  r2r3 = plate_files[29:42]
  
  # R!_R2 info
  r1r2_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(r1r2_df) = c("ORF", "Name", "Strain ID", "Penetrance_R1", "Penetrance_R2")
  
  for (plate_file in r1r2) {
    plate_csv = read.csv(plate_file)
    r1r2_df = rbind(r1r2_df, plate_csv)
  }
  
  # R1_R3 info
  r1r3_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(r1r3_df) = c("ORF", "Name", "Strain ID", "Penetrance_R1", "Penetrance_R2")
  
  for (plate_file in r1r3) {
    plate_csv = read.csv(plate_file)
    r1r3_df = rbind(r1r3_df, plate_csv)
  }
  
  # R2_R3 info
  r2r3_df = data.frame(matrix(nrow = 0, ncol = 5))
  colnames(r2r3_df) = c("ORF", "Name", "Strain ID", "Penetrance_R1", "Penetrance_R2")
  
  for (plate_file in r2r3) {
    plate_csv = read.csv(plate_file)
    r2r3_df = rbind(r2r3_df, plate_csv)
  }
  
  # Confirmed that these are the same set of strains
  r1_strains = intersect(r1r2_df$Strain.ID, r1r3_df$Strain.ID)
  r2_strains = intersect(r1r2_df$Strain.ID, r2r3_df$Strain.ID)
  r3_strains = intersect(r1r3_df$Strain.ID, r2r3_df$Strain.ID)
  core_strains = unique(c(r1_strains, r2_strains, r3_strains))
  
  no_rep = 0
  one_rep = 0
  two_reps = 0
  three_reps = 0
  for (strain in core_strains) {
    # max() is to account for strains that have multiple wells within a rep 
    # (just take the biggest penetrance, whether that's 0 or larger)
    rep1 = max(dplyr::filter(r1r2_df, r1r2_df$Strain.ID == strain)$Penetrance_R1)
    rep2 = max(dplyr::filter(r1r2_df, r1r2_df$Strain.ID == strain)$Penetrance_R2)
    rep3 = max(dplyr::filter(r1r3_df, r1r3_df$Strain.ID == strain)$Penetrance_R3)
    
    rep_pens = c(rep1, rep2, rep3)
    rep_pens_larger_zero = sum(rep_pens > 0)
    
    if (rep_pens_larger_zero == 0) {no_rep = no_rep + 1}
    if (rep_pens_larger_zero == 1) {one_rep = one_rep + 1}
    if (rep_pens_larger_zero == 2) {two_reps = two_reps + 1}
    if (rep_pens_larger_zero == 3) {three_reps = three_reps + 1}
  }
  
  marker_info = c(marker, length(core_strains), no_rep, one_rep, two_reps, three_reps)
  master_df = rbind(master_df, marker_info)
  print(marker)
}

colnames(master_df) = c("Marker", "Strains Common to All Replicates", 
                        "Strains With Cells From No Reps",
                        "Strains With Cells From One Rep",
                        "Strains With Cells From Two Reps",
                        "Strains With Cells From Three Reps")

write.csv(master_df, file = paste(input_dir, "rep_info.csv", sep = "/"), row.names = FALSE)
