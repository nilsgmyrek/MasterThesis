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
# species <- species_of_interest[4]

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

# Define model formulas to fit
model_definitions <- list(
  Original = "occu(~ TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
  Random = "occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
  Predictors_Occ = "occu(~ TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
  Predictors_Random = "occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
  RandomOnly = "occu(~ (1 | weeks) ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
  Random_Predictors = "occu(~ (1 | weeks) ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))"
)

# Initialize lists to store results
occupancy_results <- list()
best_models_aic <- list()
best_models_gof <- list()

# Loop through each species
for (species in species_of_interest) {
  
  # Initialize list to store models for each species
  fitted_models <- list()
  
  # Filter for species data
  speciesData <- det_hist_full %>%
    filter(scientificName %in% species)
  
  # Extract detection history
  detection_history <- speciesData[, grep("^2024", names(speciesData))]
  row.names(detection_history) <- speciesData$locationName
  detection_matrix <- as.matrix(detection_history)
  
  # Create unmarked frame
  umf <- unmarkedFrameOccu(y = detection_matrix, siteCovs = sitecovs, obsCovs = obscovs)
  
  # Scale continuous variables
  umf@siteCovs$Area <- scale(umf@siteCovs$Area)
  umf@siteCovs$log_Area <- scale(umf@siteCovs$log_Area)
  umf@siteCovs$min_distance_to_next_patch_km <- scale(umf@siteCovs$min_distance_to_next_patch_km)
  
  
  # Fit models
  for (model_name in names(model_definitions)) {
    model_formula <- model_definitions[[model_name]]
    
    model <- tryCatch(
      {
        eval(parse(text = model_formula))
      },
      error = function(e) {
        message(paste("Error in model", model_name, "for species", species, ":", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(model)) {
      fitted_models[[model_name]] <- model
    } else {
      cat("Skipping model", model_name, "due to fitting error.\n")
    }
  }
  
  # Find best models by AIC and GOF
  best_aic <- Inf
  best_gof <- list(pvalue = -Inf, chat_diff = Inf) # Initialize best GOF criteria
  best_model_aic <- NULL
  best_model_gof <- NULL
  
  for (model_type in names(fitted_models)) {
    current_model <- fitted_models[[model_type]]
    current_formula <- model_definitions[[model_type]]
    
    if (valid_model_check(current_model)) {
      occ_estimate <- predict(current_model, newdata = sitecovs, type = "state")
      det_estimate <- predict(current_model, newdata = umf, type = "det")
      det_estimate <- det_estimate[seq(11, nrow(det_estimate), by = 11), ]
      row.names(occ_estimate) <- sitecovs$locationName
      row.names(det_estimate) <- sitecovs$locationName
      occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      
      # Keep track 
      print(species)
      summary(current_model)
      Gof <- mb.gof.test(current_model, nsim = 10)
      Gof
      
      # Store results
      occupancy_results[[paste(species, model_type, sep = "_")]] <- data.frame(
        Patch_ID = occ_estimate$PatchID,
        Camera = occ_estimate$Row.names,
        occupancy_estimate = occ_estimate$Predicted,
        occupancy_SE = occ_estimate$SE,
        occupancy_upper = occ_estimate$upper,
        occupancy_lower = occ_estimate$lower,
        detection_estimate = det_estimate$Predicted,
        detection_SE = det_estimate$SE,
        detection_upper = det_estimate$upper,
        detection_lower = det_estimate$lower,
        formula = current_formula,
        AIC = current_model@AIC,
        GOF_chi_square = Gof$chi.square,
        GOF_p_value = Gof$p.value,
        GOF_c_hat = Gof$c.hat.est 
      )
      
      # AIC comparison
      if (current_model@AIC < best_aic) {
        best_aic <- current_model@AIC
        best_model_aic <- list(
          species = species,
          model_type = model_type,
          model = current_model,
          formula = current_formula
        )
      }
      
      # GOF comparison with criteria: p-value > 0.05 and `c-hat` closest to 1
      chat_diff <- abs(Gof$c.hat.est - 1)
      if (Gof$p.value > 0.05 && (Gof$p.value > best_gof$pvalue || (Gof$p.value == best_gof$pvalue && chat_diff < best_gof$chat_diff))) {
        best_gof <- list(pvalue = Gof$p.value, chat_diff = chat_diff)
        best_model_gof <- list(
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
  
  # Store best models for this species
  best_models_aic[[species]] <- best_model_aic
  best_models_gof[[species]] <- best_model_gof
}

# Clean environment - remove unnecessary data
rm()

# Save Workspace
#save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy1.RData")
