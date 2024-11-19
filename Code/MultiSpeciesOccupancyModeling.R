# MultiSpeciesOccu

# Clear environment
rm(list = ls())
gc()

# Disable scientific notation
options(scipen = 999) 

# Set seed for reproducibility
set.seed(111)

# Load packages
library(unmarked)
library(tidyr)
library(dplyr)
# library(AICcmodavg)

# Import Data
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Data-Prep.RData")

# Initialize an empty list to store detection matrices for each species
detection_matrices <- list()

# Loop through each species
for (species in species_of_interest) {
  
  # Filter data for the current species
  species_data <- filter(det_hist_full, scientificName == species)
  
  # Remove site-related info columns to focus on the detection history
  detection_data <- species_data %>%
    select(starts_with("2024")) %>%  # All columns starting with a date
    as.data.frame()  # Convert to a data frame
  
  # Convert detection data to a matrix and assign it to the species list
  detection_matrices[[species]] <- as.matrix(detection_data)
  
  # Optionally, assign row and column names for easier understanding
  rownames(detection_matrices[[species]]) <- species_data$locationName
  colnames(detection_matrices[[species]]) <- colnames(detection_data)
}

rm(detection_data, species_data, species)


umf <- unmarkedFrameOccuMulti(y = detection_matrices, siteCovs=sitecovs)

colnames(umf@fDesign)

stateformulas1 <- c("~ PatchID",rep("~1", 126))
stateformulas2 <- c("~ PatchID", "~ Matrix", rep("~1", 125))
stateformulas3 <- c("~ PatchID", "~ log_Area", rep("~1", 125))
stateformulas4 <- c("~ PatchID", "~ Matrix", "~ log_Area", rep("~1", 124))
stateformulas5 <- c("~ PatchID", "~ Matrix", "~ log_Area", rep("~1", 124))
stateformulas6 <- c("~ PatchID", "~ Matrix", "~ log_Area", "~ min_distance_to_next_patch_km", rep("~1", 123))

stateformulas <- c(rep("~1", 127))  # Sets all formulas to intercept-only
detformulas <- c("~1","~1","~1","~1","~1","~1","~1")

  Model1.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas1, data = umf, control = list(maxit = 1000))
  Model2.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas2, data = umf, control = list(maxit = 1000))
  Model3.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas3, data = umf, control = list(maxit = 1000))
  Model4.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas4, data = umf, control = list(maxit = 1000))
  Model5.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas5, data = umf, control = list(maxit = 1000))
  Model6.4 <- occuMulti(detformulas=detformulas, stateformulas=stateformulas6, data = umf, control = list(maxit = 1000))

