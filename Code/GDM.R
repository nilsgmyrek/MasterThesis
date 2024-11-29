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
