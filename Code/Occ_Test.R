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

# Define model formulas to fit
model_definitions <- list(
  Model1 = "occu(~ (1 | weeks) ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
  Model2 = "occu(~ (1 | weeks) ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
  Model2 = "occu(~ (1 | weeks) ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
  Model3 = "occu(~ (1 | weeks) ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
  Model4 = "occu(~ (1 | weeks) ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
  Model5 = "occu(~ (1 | weeks) ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
  Model6 = "occu(~ (1 | weeks) ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))"
)

# Initialize lists to store results
aic_results <- list()
occupancy_results <- list()
best_model_aic_results <- list()  # Stores the best models by AIC for each species

# Loop through each species
for (species in species_of_interest) {
  
  # Initialize list to store models and AIC for each species
  fitted_models <- list()
  species_aic <- data.frame(model = character(), AIC = numeric(), stringsAsFactors = FALSE)
  best_aic <- Inf  # Track the best (lowest) AIC for comparison
  best_model_aic <- NULL  # Placeholder for the best model
  
  # Filter for species data
  speciesData <- det_hist_full %>%
    filter(scientificName == species)
  
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
  
  # Fit models and collect AIC and occupancy results
  for (model_name in names(model_definitions)) {
    model_formula <- model_definitions[[model_name]]
    current_formula <- model_formula  # Store current formula for reference in output
    
    model <- tryCatch(
      {
        eval(parse(text = model_formula))
      },
      error = function(e) {
        message(paste("Error in model", model_name, "for species", species, ":", e$message))
        return(NULL)
      }
    )
    
    # If model fitting succeeds, store AIC and Model
    if (!is.null(model)) {
      fitted_models[[model_name]] <- model
      current_aic <- model@AIC
      species_aic <- rbind(species_aic, data.frame(model = model_name, AIC = current_aic))
      
      # Check if this model has the lowest AIC so far
      if (current_aic < best_aic) {
        best_aic <- current_aic
        best_model_aic <- list(
          species = species,
          model_type = model_name,
          model = model,
          formula = current_formula
        )
      }
      
      # Predict occupancy and detection estimates
      occ_estimate <- predict(model, newdata = sitecovs, type = "state")
      det_estimate <- predict(model, newdata = umf, type = "det")
      
      # Restructure detection estimates to match rows with `sitecovs`
      det_estimate <- det_estimate[seq(11, nrow(det_estimate), by = 11), ]
      row.names(occ_estimate) <- sitecovs$locationName
      row.names(det_estimate) <- sitecovs$locationName
      occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
      
      # Goodness-of-fit test
#     Gof <- mb.gof.test(model, nsim = 10)
      
      # Store results
      occupancy_results[[paste(species, model_name, sep = "_")]] <- data.frame(
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
        AIC = current_aic #,
        # GOF_chi_square = Gof$chi.square,
        # GOF_p_value = Gof$p.value,
        # GOF_c_hat = Gof$c.hat.est
      )
    } else {
      cat("Skipping model", model_name, "due to fitting error.\n")
    }
  }
  
  # Store AIC results for the species
  aic_results[[species]] <- species_aic
  
  # Store the best model by AIC for the current species
  if (!is.null(best_model_aic)) {
    best_model_aic_results[[species]] <- best_model_aic
  }
}

# Clean environment - remove unnecessary data
rm(obscovs, sitecovs)

rm(best_model_aic , det_estimate, detection_history, detection_matrix, fitted_models, model, model_definitions, occ_estimate, species_aic, speciesData, umf, best_aic, current_aic, current_formula, model_formula, model_name, species)

# Save Workspace
save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy.RData")

