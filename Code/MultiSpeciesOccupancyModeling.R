# MultiSpeciesOccu: Multi-species occupancy modeling workflow

# Clear environment to ensure a clean workspace
rm(list = ls())  # Remove all objects from the environment
gc()  # Trigger garbage collection to free memory

# Disable scientific notation for better readability of outputs
options(scipen = 999) 

# Set seed for reproducibility of results
set.seed(111)

# Load required packages
library(unmarked)  # For occupancy modeling
library(tidyr)     # For data manipulation
library(dplyr)     # For data wrangling
library(AICcmodavg)  # For model comparison 

# Load pre-processed data
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Data-Prep.RData")

# Initialize an empty list to store detection matrices for each species
detection_matrices <- list()

# Loop through each species of interest
for (species in species_of_interest) {
  
  # Filter detection history for the current species
  species_data <- filter(det_hist_patch, scientificName == species)
  
  # Extract only detection history columns (e.g., date columns) and convert to a data frame
  detection_data <- species_data %>%
    select(starts_with("2024")) %>%  # Extract detection history
    as.data.frame()
  
  # Convert detection data to a matrix format for use in unmarkedFrame
  detection_matrices[[species]] <- as.matrix(detection_data)
  
  # Assign meaningful row and column names for the matrix
  rownames(detection_matrices[[species]]) <- species_data$PatchID
  colnames(detection_matrices[[species]]) <- colnames(detection_data)
}

# Clean up intermediate objects
rm(detection_data, species_data, species)

# Create an unmarkedFrame for multi-species occupancy modeling
umf <- unmarkedFrameOccuMulti(
  y = detection_matrices,  # Detection matrices for each species
  siteCovs = p_sitecovs      # Site-specific covariates
)

# Check the design matrix for state formulas
colnames(umf@fDesign)  # Lists covariates used for state formulas

# Define state formulas (occupancy) for testing different covariate combinations
# stateformulas <- c(rep("~1", 7))  # Intercept-only model for all species
# stateformulas0 <- c(rep("~ PatchID", 7))  # Include PatchID as a covariate
# stateformulas1 <- c(rep("~ PatchID + Matrix", 7))  # Add habitat matrix
# stateformulas2 <- c(rep("~ PatchID + log_Area", 7))  # Add log-transformed area
# stateformulas3 <- c(rep("~ PatchID + Matrix + log_Area", 7))  # Combine covariates
# stateformulas4 <- c(rep("~ PatchID + Matrix * log_Area", 7))  # Add interaction term
# stateformulas5 <- c(rep("~ PatchID + Matrix + log_Area + min_distance_to_next_patch_km", 7))  # Add Connectivity
# stateformulas6 <- c(rep("~ PatchID + Matrix * log_Area + min_distance_to_next_patch_km", 7))  # Add interaction term

# Define state formulas (occupancy) for testing different covariate combinations
stateformulas0 <- c(rep("~ (1 | PatchID)", 7))  # Intercept-only model for all species
stateformulas1 <- c(rep("~ Matrix", 7))  # Add habitat matrix
stateformulas2 <- c(rep("~ log_Area", 7))  # Add log-transformed area
stateformulas3 <- c(rep("~ Matrix + log_Area", 7))  # Combine covariates
stateformulas4 <- c(rep("~ Matrix * log_Area", 7))  # Add interaction term
stateformulas5 <- c(rep("~ Matrix + log_Area + min_distance_to_next_patch_km", 7))  # Add Connectivity
stateformulas6 <- c(rep("~ Matrix * log_Area + min_distance_to_next_patch_km", 7))  # Add interaction term


# Detection formulas (shared for all models)
# detformulas <- c(rep("~ 1", 7))
detformulas <- c("~n_cameras", "~n_cameras", "~n_cameras", "~n_cameras", "~n_cameras", "~n_cameras", "~n_cameras")  # number of camera traps per patch
# detformulas <- c(rep("~ VegetationCover", 7)) # VegetationCover as Detection Covariate for all Species
# detformulas <- c(rep("~ TreeDensity", 7)) # TreeDensity as Detection Covariate for all Species
# detformulas <- c(rep("~ VegetationCover + TreeDensity", 7)) # VegetationCover and TreeDensity as Detection Covariates for all Species

# Fit multi-species occupancy models with different covariates
#Model  <- occuMulti(detformulas = detformulas, stateformulas = stateformulas , data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model0 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas0, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model1 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas1, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model2 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas2, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model3 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas3, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model4 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas4, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model5 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas5, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)
Model6 <- occuMulti(detformulas = detformulas, stateformulas = stateformulas6, data = umf, control = list(maxit = 1000, trace = 1), maxOrder = 1)

# Save Workspace 
#save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/MultiOcc_results.RData")
#load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/MultiOcc_results.RData")

# Store models in a list
models <- list(
# Model  = Model,
  Model0 = Model0,
  Model1 = Model1,
  Model2 = Model2,
  Model3 = Model3#,
#  Model4 = Model4,
#  Model5 = Model5 #,
#  Model6 = Model6
)

# Perform AIC model comparison
aictab(cand.set = models, modnames = names(models))

summary(Model2)
