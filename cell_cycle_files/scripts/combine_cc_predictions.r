# Combine cell cycle predictions for each marker
# Purpose: CC predictions are separated by plate, pool them together first to 
# make future analysis easier

# Alex Daiejavad
# Date: January 20, 2023
#-------------------------------------------------------------------------------
input_f = "C:/Users/peree/OneDrive/Desktop/cell_cycle"
output_f = "C:/Users/peree/OneDrive/Desktop/combined_cell_cycle"

#markers = c("Cdc11", "Dad2", "Heh2", "Nop10", "Nuf2", "Om45", "Pil1", "Psr1", 
#            "Rad52", "Sac6", "Sec7", "Sec21", "Snf7", "Spf1", "Vph1")

markers = c("Vph1")

for (m in markers) { # all markers except vph1
  m_screen = paste("Screens_", m, sep = "")
  setwd(paste(input_f, m_screen, sep = "/"))
  all_screens = list.files()
  
  combined_predictions = data.frame(matrix(nrow = 0, ncol = 2))
  colnames(combined_predictions) = c("Cell_ID", "Prediction")
  for (screen in all_screens) {
    screen_predictions = read.csv(file = paste(input_f, m_screen, screen, 
                                               "Predictions_cellcycle_stage.csv", sep = "/"))
    screen_predictions = data.frame(Cell_ID = screen_predictions$cell_id, 
                                    Prediction = screen_predictions$Prediction)
    
    combined_predictions = rbind(combined_predictions, screen_predictions)
  }
  file_name = paste(m, "all_mpg3_predictions.csv", sep = "_")
  write.csv(combined_predictions, 
            file = paste(output_f, file_name, sep = "/"),
            row.names = FALSE)
}


for (m in markers) { # vph1 only
  m_screen = paste("Screens_", m, sep = "")
  setwd(paste(input_f, m_screen, sep = "/"))
  all_screens = list.files()
  relevant_screens = c(all_screens[grepl("Vph1_R1", all_screens)], all_screens[grepl("Vph1_TS1", all_screens)])
  
  combined_predictions = data.frame(matrix(nrow = 0, ncol = 2))
  colnames(combined_predictions) = c("Cell_ID", "Prediction")
  for (screen in relevant_screens) {
    screen_predictions = read.csv(file = paste(input_f, m_screen, screen, 
                                               "Predictions_cellcycle_stage.csv", sep = "/"))
    screen_predictions = data.frame(Cell_ID = screen_predictions$cell_id, 
                                    Prediction = screen_predictions$Prediction)
    
    combined_predictions = rbind(combined_predictions, screen_predictions)
  }
  file_name = paste(m, "R1_predictions.csv", sep = "_")
  write.csv(combined_predictions, 
            file = paste(output_f, file_name, sep = "/"),
            row.names = FALSE)
}
