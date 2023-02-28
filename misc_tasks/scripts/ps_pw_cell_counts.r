library(dplyr)

setwd("C:/Users/peree/OneDrive/Desktop/Cdc11_all_PCA")

per_well = read.csv("Cdc11_all_PCA_GMM_OD_results_well.csv")
per_strain = read.csv("Cdc11_all_PCA_GMM_OD_results_strain.csv")
outliers = read.csv("Cdc11_all_PCA_GMM_outlier_cells.csv")

ps_pw = per_strain[,1:5]
strains = per_strain$Strain.ID
# Test 1: making sure total cell counts for each strain is the same between PW 
# and PS
pw_counts = c()
for (strain in strains) {
  pw_filter = dplyr::filter(per_well, per_well$Strain.ID == strain)
  pw_counts = c(pw_counts, sum(pw_filter$Num_cells))
}
colnames(ps_pw) = c('ORF', 'Name', 'Allele', 'StrainID', 'Num_Cells_PS')
ps_pw = cbind(ps_pw, Num_Cells_PW = pw_counts)
sum(ps_pw$Num_Cells_PS == ps_pw$Num_Cells_PW)  # 6243 -> same number as all strains so no discrepancies


# Test 2: making sure total outliers for each strain < Num_Cells
outlier_counts = c()
for (strain in strains) {
  outlier_filter = dplyr::filter(outliers, outliers$Strain.ID == strain)
  outlier_counts = c(outlier_counts, nrow(outlier_filter))
}
ps_pw = cbind(ps_pw, HighConf_Outliers = outlier_counts)
sum(ps_pw$Num_Cells_PS > ps_pw$HighConf_Outlier)  # 6243 -> no discrepancies


# Test 3: per-strain outlier counts using Num_Cells * KS_Penetrance
ks_outliers_strain = c()
for (strain in strains) {
  ps_filter = dplyr::filter(per_strain, per_strain$Strain.ID == strain)
  ks_outliers_strain = c(ks_outliers_strain, ps_filter$Num_cells * (ps_filter$KS_Penetrance / 100))
}
ps_pw = cbind(ps_pw, KS_Outliers_Strain = ks_outliers_strain)
sum(ps_pw$KS_Outliers_Strain, na.rm = TRUE)  # 260903


# Test 4: per-well outlier counts using Num_Cells * KS_Penetrance
ks_outliers_well = c()
for (strain in strains) {
  pw_filter = dplyr::filter(per_well, per_well$Strain.ID == strain)
  ks_outliers_well = c(ks_outliers_well, sum(pw_filter$Num_cells * (pw_filter$KS_Penetrance / 100)))
}
ps_pw = cbind(ps_pw, KS_Outliers_Well = ks_outliers_well)
sum(ps_pw$KS_Outliers_Well, na.rm = TRUE)  # 333899

write.csv(ps_pw, "strain_level_cell_counts.csv", row.names = FALSE)
#-------------------------------------------------------------------------------
#pw = dplyr::filter(per_well, per_well$Strain.ID == "tsa31-37C")
#ps = dplyr::filter(per_strain, per_strain$Strain.ID == "tsa31-37C")
#out = dplyr::filter(outliers, outliers$Strain.ID == "tsa31-37C")
pspw = dplyr::filter(ps_pw, ps_pw$StrainID == 'dma3345')

new_per_well = cbind(per_well, KS_Outliers = per_well$Num_cells * (per_well$KS_Penetrance / 100))
new_per_well = new_per_well[order(new_per_well$Strain.ID),]
strains = sort(unique(per_well$Strain.ID))

cell_counts = c()

for (strain in strains) {
  fil = dplyr::filter(per_well, per_well$Strain.ID == strain)
  cell_counts = c(cell_counts, sum(fil$KS_Outliers))
}
new_pspw = ps_pw[order(ps_pw$StrainID),]
new_pspw = cbind(new_pspw, PS_Outliers = cell_counts)
