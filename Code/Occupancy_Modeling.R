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
library(AICcmodavg)

# Import Data
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Data-Prep.RData")

## Occupancy Modeling ####
# test Loop
# species <- species_of_interest[3]

# Helper function to check for NA in model summary components 
valid_model_check <- function(model) {
  # Check convergence first
  if (model@opt$convergence == 0) {
    # Extract the summary of the model
    model_summary <- summary(model)
    
    # Check for NAs in fixed effects (Occupancy) and detection (logit-scale) estimates, SE, z, and P values
    fixed_effects <- model_summary$state$coefficients
    detection_effects <- model_summary$det$coefficients
    
    # Check if any NA values are present in coefficients, SE, z, or P
    if (any(is.na(fixed_effects)) || any(is.na(detection_effects))) {
      return(FALSE)
    }
    
    # Check if there are NAs in SE, z, or P in both occupancy and detection
    fixed_effects_full <- cbind(fixed_effects, model_summary$state$SE, model_summary$state$z, model_summary$state$p)
    detection_effects_full <- cbind(detection_effects, model_summary$det$SE, model_summary$det$z, model_summary$det$p)
    
    if (any(is.na(fixed_effects_full)) || any(is.na(detection_effects_full))) {
      return(FALSE)
    }
    
    return(TRUE)
  }
  return(FALSE)
}

# Create empty lists to store results for each model type
occupancy_results <- list()

# Initialize list to store best models per species
best_models <- list()

# Loop through each species
for (species in species_of_interest) {
  
  # Filter for species data
  speciesData <- det_hist_full %>%
    filter(scientificName %in% species)
 
  # Extract detection history (columns with dates)
  detection_history <- speciesData[, grep("^2024", names(speciesData))]
  
  # Assign row names using the locationName column
  row.names(detection_history) <- speciesData$locationName
  
  # Convert the detection history into the correct format for unmarked
  detection_matrix <- as.matrix(detection_history)
  
  # Create the unmarked frame
  umf <- unmarkedFrameOccu(y = detection_matrix, siteCovs = sitecovs, obsCovs = obscovs)
  
  # Scale continuous variables 
  umf@siteCovs$log_Area <- scale(umf@siteCovs$Area)
  umf@siteCovs$min_distance_to_next_patch_km <- scale(umf@siteCovs$min_distance_to_next_patch_km)
  
  # Fit both models
  occ_model <- occu(~ TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_random_model <- occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_model_pred <- occu(~ TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))
  occ_random_model_pred <- occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))
  occ_random <- occu(~ (1 | weeks) ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_random_pred <- occu(~ (1 | weeks) ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))
  
  # Initialize the model list
  model_list <- list(
    Original = list(model = occ_model, formula = "TreeDensity + VegetationCover ~ (1 | PatchID)"),
    Random = list(model = occ_random_model, formula = "(1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID)"),
    Predictors_Occ = list(model = occ_model_pred, formula = "TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km"),
    Predictors_Random = list(model = occ_random_model_pred, formula = "(1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km"),
    RandomOnly = list(model = occ_random, formula = "(1 | weeks) ~ (1 | PatchID)"),
    Random_Predictors = list(model = occ_random_pred, formula = "(1 | weeks) ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km")
  )
  
  # Track the best model for the species
  best_aic <- Inf
  best_model_info <- NULL
  
  # Loop through each model type
  for (model_type in names(model_list)) {
    current_model <- model_list[[model_type]]$model
    current_formula <- model_list[[model_type]]$formula
    
    # Validate the model before proceeding
    if (valid_model_check(current_model)) {
      
      # Get Occupancy Estimates
      occ_estimate <- predict(current_model, newdata = sitecovs, type = "state")
      row.names(occ_estimate) <- sitecovs$locationName
      occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      
      # Get Detection Estimates
      det_estimate <- predict(current_model, newdata = umf, type = "det")
      det_estimate <- det_estimate[seq(9, nrow(det_estimate), by = 9), ]
      row.names(det_estimate) <- sitecovs$locationName
      det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      
      # Goodness of Fit Test
      Gof <- mb.gof.test(current_model, nsim = 1000)
      
      # Combine estimates into final results
      results <- merge(occ_estimate, det_estimate, by = "Row.names", suffixes = c("_occ", "_det"))
      
      # Save the results for this model and species
      occupancy_results[[paste(species, model_type, sep = "_")]] <- data.frame(
        Patch_ID = results$PatchID_occ,
        Camera = results$Row.names,
        occupancy_estimate = results$Predicted_occ,
        occupancy_SE = results$SE_occ,
        occupancy_upper = results$upper_occ,
        occupancy_lower = results$lower_occ,
        detection_estimate = results$Predicted_det,
        detection_SE = results$SE_det,
        detection_upper = results$upper_det,
        detection_lower = results$lower_det,
        formula = current_formula,
        AIC = current_model@AIC,
        GoodnessOfFit_chi_square = Gof$chi.square,
        GoodnessOfFit_p_value = Gof$p.value,
        GoodnessOfFit_c_hat = Gof$c.hat.est 
      )
      
      # Check if this model has the lowest AIC
      if (current_model@AIC < best_aic) {
        best_aic <- current_model@AIC
        best_model_info <- list(
          species = species,
          model_type = model_type,
          model = current_model,
          formula = current_formula
        )
      }
    } else {
      cat("Model for", species, model_type, "did not pass validation and was skipped.\n")
    }
  }
  
  # Store the best model information for this species
  best_models[[species]] <- best_model_info
}

# Clean environment - remove unnecessary data
rm(species, current_model, model_list, occ_model_pred, occ_random, occ_random_pred, best_aic, best_model_info, occ_random_model_pred, current_formula, model_type, occ_model, occ_random_model, occ_estimate, det_estimate, detection_history, detection_matrix, Gof, results, sitecovs, speciesData, umf)

# Save Workspace 
save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy.RData")
