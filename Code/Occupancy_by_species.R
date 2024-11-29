# Models for species occupancy

# Clean environment
rm(list = ls())

# Disable scientific notation
options(scipen = 999) 

# Set seed for reproducibility
set.seed(111)

# Load packages
library(unmarked)
library(tidyr)
library(dplyr)
# library(AICcmodavg)


### Model selection function ###

# model.sel funciton. Inputs:
# -- det.hist: data frame with the species detection history, including all sp.
# ----- the column for the species is named 'scientificName'
# ----- The column for locations is named 'locationName'
# -- model.list: list of models to be compared
# siteCovs, obsCovs: matrix of site and observation covariates, respectively

model.sel <- function(det.hist, model.list, siteCovs, obsCovs) {
  
  sp.list <- unique(det.hist$scientificName)
  
  # Initialize containers for results
  Model_container <- list()         # Store all fitted models
  aic_results <- list()             # Store AIC results per species
  best_model_aic_results <- list()  # Store best model per species
  occupancy_results <- list()       # Store occupancy predictions
  
  for (species in sp.list) {
    # Initialize lists to store models and AIC for each species
    fitted_models <- list()
    species_aic <- data.frame(model = character(), AIC = numeric(), stringsAsFactors = FALSE)
    best_aic <- Inf
    best_model_aic <- NULL
    
    # Filter for species data
    speciesData <- det.hist %>% # Detection history
      filter(scientificName == species)
    
    # Extract detection history for one species
  # locationName <- speciesData$locationName
    speciesData <- as.matrix(speciesData[, grep("^2024", names(speciesData))])
  # row.names(speciesData) <- locationName
    
    # Create unmarked frame
    umf <- unmarkedFrameOccu(y = speciesData, siteCovs = siteCovs, obsCovs = obsCovs)  
    
    # Fit models and collect AIC results
    for (model_name in names(model.list)) {
      model_formula <- model.list[[model_name]]
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
        
        # Check for best AIC model
        if (current_aic < best_aic) {
          best_aic <- current_aic
          best_model_aic <- list(
            species = species,
            model_type = model_name,
            model = model,
            formula = current_formula
          )
        }
        
        #Predict occupancy and detection estimates
        occ_estimate <- predict(model, newdata = sitecovs, type = "state")
        det_estimate <- predict(model, newdata = umf, type = "det")
        
        # Restructure detection estimates to match rows with `sitecovs`
        det_estimate <- det_estimate[seq(11, nrow(det_estimate), by = 11), ]
        row.names(occ_estimate) <- sitecovs$locationName
        row.names(det_estimate) <- sitecovs$locationName
        occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
        det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)

        # # Store results
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
        AIC = current_aic
        )
        # Store Model
        Model_container[[paste(species, model_name, sep = "_")]] <- model  
      } 
    }
    
    # Store AIC results for the species
    aic_results[[species]] <- species_aic
    
    # Store the best model by AIC for the current species
    if (!is.null(best_model_aic)) {
      best_model_aic_results[[species]] <- best_model_aic
    }
  }
  
  return(list(all_models = Model_container, aic_results = aic_results, best_models = best_model_aic_results, occupancy = occupancy_results))
}


## Occupancy model selection analyses

# 1. Read data
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Data-Prep.RData")


# 2. Define models to be tested

# 2.1 
# models_all <- list(
#   Model1.1 = "occu(~ 1 ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
#   Model1.2 = "occu(~ VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
#   Model1.3 = "occu(~ TreeDensity ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
#   Model1.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
#   Model2.1 = "occu(~ 1 ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
#   Model2.2 = "occu(~ VegetationCover ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
#   Model2.3 = "occu(~ TreeDensity ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
#   Model2.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
#   Model3.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
#   Model3.2 = "occu(~ VegetationCover ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
#   Model3.3 = "occu(~ TreeDensity ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
#   Model3.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
#   Model4.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
#   Model4.2 = "occu(~ VegetationCover ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
#   Model4.3 = "occu(~ TreeDensity ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
#   Model4.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
#   Model5.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
#   Model5.2 = "occu(~ VegetationCover ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
#   Model5.3 = "occu(~ TreeDensity ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
#   Model5.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
#   Model6.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model6.2 = "occu(~ VegetationCover ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model6.3 = "occu(~ TreeDensity ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model6.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model7.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model7.2 = "occu(~ VegetationCover~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model7.3 = "occu(~ TreeDensity ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#   Model7.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))"
#   )

# compare.detect <- model.sel(det_hist_full, models_all, sitecovs, NULL)


# 2.2 

# models_detect_occu <- list(
#  Model1.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
#  Model2.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
#  Model3.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
#  Model4.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
#  Model5.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
#  Model6.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
#  Model7.4 = "occu(~ VegetationCover + TreeDensity ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))"
# )

# Run function
# compare.detect.occu <- model.sel(det_hist_full, models_detect_occu, sitecovs, NULL)


# 2.3

models_occu <- list(
  Model1.1 = "occu(~ 1 ~ (1 | PatchID), data = umf, control = list(maxit = 1000))",
  Model2.1 = "occu(~ 1 ~ (1 | PatchID) + Matrix, data = umf, control = list(maxit = 1000))",
  Model3.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area, data = umf, control = list(maxit = 1000))",
  Model4.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area + Matrix, data = umf, control = list(maxit = 1000))",
  Model5.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area * Matrix, data = umf, control = list(maxit = 1000))",
  Model6.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))",
  Model7.1 = "occu(~ 1 ~ (1 | PatchID) + log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))"
)

compare.occu <- model.sel(det_hist_full, models_occu, sitecovs, NULL)


# # 2.4
# models_patch <- list(
#   Model1.1 = "occu(~ 1 ~ 1, data = umf, control = list(maxit = 5000))",
#   Model2.1 = "occu(~ 1 ~ Matrix, data = umf, control = list(maxit = 5000))",
#   Model3.1 = "occu(~ 1 ~ log_Area, data = umf, control = list(maxit = 5000))",
#   Model4.1 = "occu(~ 1 ~ log_Area + Matrix, data = umf, control = list(maxit = 5000))",
#   Model5.1 = "occu(~ 1 ~ log_Area * Matrix, data = umf, control = list(maxit = 5000))",
#   Model6.1 = "occu(~ 1 ~ log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 5000))",
#   Model7.1 = "occu(~ 1 ~ log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 5000))"
# )
# 
# 
# models_patch <- list(
#   Model1.1 = "occu(~ n_cameras ~ 1, data = umf, control = list(maxit = 5000))",
#   Model2.1 = "occu(~ n_cameras ~ Matrix, data = umf, control = list(maxit = 5000))",
#   Model3.1 = "occu(~ n_cameras ~ log_Area, data = umf, control = list(maxit = 5000))",
#   Model4.1 = "occu(~ n_cameras ~ log_Area + Matrix, data = umf, control = list(maxit = 5000))",
#   Model5.1 = "occu(~ n_cameras ~ log_Area * Matrix, data = umf, control = list(maxit = 5000))",
#   Model6.1 = "occu(~ n_cameras ~ log_Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 5000))",
#   Model7.1 = "occu(~ n_cameras ~ log_Area * Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 5000))"
# )
# 
# compare.patch <- model.sel(det_hist_patch, models_patch, p_sitecovs, NULL)


# # Save Workspace 
save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occ_results.RData")

# -----------------------------
# PLAYGROUND
