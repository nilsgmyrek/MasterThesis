# Clear environment
rm(list = ls())
gc()

# Disable scientific notation
options(scipen = 999) 
# options(unmarked.maxit = 10000)  # number of iterations for unmarked

# Load packages
library(lubridate)
library(unmarked)
library(tidyr)
library(dplyr)
library(sf)
library(AICcmodavg)

## Import and prepare Data #### 
# Load Presence-Absence Data
# presence_absence_df <- st_read('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/presence_absence_by_patch_new.gpkg') %>%
#   st_drop_geometry() %>%
#   write.csv('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/presence_absence_by_patch_new.csv')

# Load Presence-Absence Data
presence_absence_df <- read.csv('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/presence_absence_by_patch_new.csv')

# Select relevant columns from presence_absence_df
presence_absence_df <- presence_absence_df[, c("Patch_ID", "Perimeter", "Area", "CT", "SU_trunc", "min_distance_to_next_patch_km")]

# Load Covariates
covs <- read.csv('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/covariate_estimation.csv')

# Load Deployment-Data
deployments <- read.csv('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Agouti Output Aug/deployments.csv')

# Ensure deploymentStart and deploymentEnd are in POSIXct format
deployments <- deployments %>%
  mutate(
    deploymentStart = as.POSIXct(deploymentStart),
    deploymentEnd = as.POSIXct(deploymentEnd)
  )

# Correct the year to 2024 and fix the "LE39" issue
deployments <- deployments %>%
  mutate(
    # Change year to 2024 if it is not already 2024
    deploymentStart = if_else(year(deploymentStart) != 2024, 
                              update(deploymentStart, year = 2024), 
                              deploymentStart),
    deploymentEnd = if_else(year(deploymentEnd) != 2024, 
                            update(deploymentEnd, year = 2024), 
                            deploymentEnd)
  ) %>%
  mutate(
    # Correct the deploymentStart for location "LE39"
    deploymentStart = if_else(locationName == "LE39", 
                              as.POSIXct("2024-03-04 00:00:00"), 
                              deploymentStart)
  )

# Merge Patch characteristics with deployments
deployments <- merge(deployments, presence_absence_df, by.x = "deploymentTags", by.y = "Patch_ID", all.x = TRUE)

# Merge deployments with Covariates
deployments <- merge(deployments, covs, by.x = "locationName", by.y = "locationName", all.x = T)

# Check Structure of Deployments-Table
str(deployments)

# Load Observations-Data
observations <- read.csv('/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Agouti Output Aug/observations.csv')

# Update all years to 2024 except for 2024
observations <- observations %>%
  mutate(
    # Parse the eventStart and eventEnd columns into POSIXct objects
    eventStart = ymd_hms(eventStart, tz = "Europe/Berlin"),
    eventEnd = ymd_hms(eventEnd, tz = "Europe/Berlin"),
    
    # Change year to 2024 if it is not already 2024
    eventStart = if_else(year(eventStart) != 2024, 
                         update(eventStart, year = 2024), 
                         eventStart),
    eventEnd = if_else(year(eventEnd) != 2024, 
                       update(eventEnd, year = 2024), 
                       eventEnd),
    
    # Convert back to ISO 8601 datetime strings
    eventStart = format(eventStart, "%Y-%m-%dT%H:%M:%S%z"),
    eventEnd = format(eventEnd, "%Y-%m-%dT%H:%M:%S%z")
  )

# Convert dates to proper format
observations$eventStart <- as.POSIXct(observations$eventStart)
observations$eventEnd <- as.POSIXct(observations$eventEnd)

# Merge locationName from deployments to observations based on deploymentID
observations <- observations %>%
  left_join(deployments %>% select(deploymentID, locationName), by = "deploymentID")

# Check structure of Observations-Table
str(observations)

# Set Survey Start-Date
# start_date <- as.Date(min(deployments$deploymentStart))

# Calculate the week number
# observations$week <- as.numeric(difftime(observations$eventStart, start_date, units = "weeks")) + 1


## Create detection History ####
# Define Species of Interest
species_of_interest <- c("Capreolus capreolus", "Martes", "Felis silvestris", "Sus scrofa", "Procyon lotor", "Meles meles", "Vulpes vulpes")

# Filter observations for species of interest
obs_filtered <- observations %>%
  filter(scientificName %in% species_of_interest)

# Define weekly time bins
obs_filtered <- obs_filtered %>%
  mutate(week = floor_date(eventStart, "week")) %>%
  mutate(day = floor_date(eventStart, "day"))

# Create detection history matrix - Weeks
det_hist <- obs_filtered %>%
  group_by(deploymentID, week, scientificName) %>%
  summarize(detection = 1) %>%
  spread(key = week, value = detection, fill = 0)

# remove first and last week of deployment to make things more even
det_hist <-  det_hist[,-c(3,13)]

# # Create detection history matrix - Days
# det_hist <- obs_filtered %>%
#   group_by(deploymentID, day, scientificName) %>%
#   summarize(detection = 1) %>%
#   spread(key = day, value = detection, fill = 0)
# 
# # remove all days were not all cameras were switched on
# det_hist <-  det_hist[,-c(3:15, 52:72)]

# Add missing species
all_combinations <- expand.grid(
  deploymentID = unique(deployments$deploymentID),
  scientificName = unique(det_hist$scientificName)
)

existing_combinations <- det_hist %>%
  select(deploymentID, scientificName) %>%
  distinct()

missing_combinations <- anti_join(all_combinations, existing_combinations, by = c("deploymentID", "scientificName"))
missing_combinations

# Combine with Existing Data
det_hist_full <- bind_rows(det_hist, missing_combinations) %>%
  distinct() %>%
  mutate(across(everything(), ~replace_na(., 0)))

# Add deployments that didn't detect the species (fill missing deployments)
det_hist_full <- deployments %>%
  left_join(det_hist_full, by = "deploymentID") %>%
  select(locationName, scientificName, everything())

# View the completed data
det_hist_full <-  det_hist_full[,-c(4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,25)]
str(det_hist_full)

# Remove unnecessary Data
rm(all_combinations, covs, deployments, det_hist, existing_combinations, missing_combinations, observations, start_date, presence_absence_df, obs_filtered)

# Filter Patches with wrong number of camera traps
det_hist_full <- det_hist_full %>%
  filter(CT == SU_trunc) %>%  # Keep only rows where SU matches SU_trunc
  select(-CT,-SU_trunc) %>%
  rename(Matrix = deploymentGroups) %>%
  rename(PatchID = deploymentTags)



## Create Observation Covariates ####
# Prepare Observation-Covariate Data - weeks
obscovs <- list()

sdraft <- det_hist_full %>%
  filter(scientificName %in% 'Capreolus capreolus')

draft <-  sdraft[, grep("^2024", names(sdraft))]
#row.names(draft) <- sdraft$locationName

for (i in 1:ncol(draft)) {
  draft[,i] <- i
}

# Create Week Observation Matrix
weeks <- as.matrix(draft)

# As Characters - categorical vs. continuous variables
weeks <- apply(weeks, 2, as.factor)
row.names(weeks) <- sdraft$locationName

# Add 'weeks' matrix to the existing 'obscovs' list
obscovs$weeks <- weeks

# remove unnecessary Data
rm(draft, sdraft, weeks, i)


## Simple Occupancy Model ####
# Occupancy List
occupancy_results <- list()

# test Loop
# species <- species_of_interest[1]

# Loop through each species
for (species in species_of_interest) {

  # filter for Species
  speciesData <- det_hist_full %>%
    filter(scientificName %in% species)

  # create Site specific Covariate-Data
  sitecovs <- speciesData %>%
    select(locationName, PatchID, TreeDensity, VegetationCover, cameraHeight, Matrix, min_distance_to_next_patch_km, Area) %>%
    mutate_if(is.character, as.factor)
  
  # Extract the detection history (columns with dates)
  detection_history <- speciesData[, grep("^2024", names(speciesData))]

  # Assign row names using the locationName column
  row.names(detection_history) <- speciesData$locationName

  # Convert the detection history into the correct format for unmarked
  detection_matrix <- as.matrix(detection_history)

  # Check structure of Covariate Table
  str(sitecovs)
  
  # Create the unmarked frame for the current row
  umf <- unmarkedFrameOccu(y = detection_matrix, siteCovs = sitecovs, obsCovs = obscovs) 

  # Fit the models
  occ_model <- occu( ~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID), data = umf)
  occ_random_model <- occu( ~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf)

  # Summary
  summary(occ_model)
  summary(occ_random_model)

  # Access AIC values directly from the models
  aic_occ_model <- occ_model@AIC
  aic_occ_random_model <- occ_random_model@AIC

  print(species)

  # Compare the AIC values and store the model with the lower AIC
  if (aic_occ_model < aic_occ_random_model) {
    occ_model <- occ_model  # Keep the original model
    cat("Original Model selected with AIC:", aic_occ_model, "\n")
  } else {
    occ_model <- occ_random_model  # Assign the random model
    cat("Random Model selected with AIC:", aic_occ_random_model, "\n")
  }

  # Get the Occupancy Estimates
  occ_estimate <-  predict(occ_model,
                           newdata = sitecovs,
                           type = "state")

  # Add Identifiers to Occupancy estimates
  row.names(occ_estimate) <- sitecovs$locationName
  occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)

  # Get Detection Estimates
  det_estimate <- predict(occ_model,
                          newdata = umf,
                          type = "det")
  
  # remove duplicates
  det_estimate <- det_estimate[seq(9, nrow(det_estimate), by = 9), ]
  
  # Add Identifiers to Occupancy estimates
  row.names(det_estimate) <- sitecovs$locationName
  det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)

  # Godness of Fit-Test
  Gof <- mb.gof.test(occ_model, nsim = 100)

  # Ensure that both estimates are correctly combined
  results <- merge(occ_estimate, det_estimate, by = "Row.names", suffixes = c("_occ", "_det"))

  # Create the final data frame
  occupancy_results[[species]] <- data.frame(
    Patch_ID = results$PatchID_occ,                # PatchID 
    Camera = results$Row.names,                    # Camera ID (locationName)
    occupancy_estimate = results$Predicted_occ,    # Occupancy Predictions
    occupancy_SE = results$SE_occ,                 # Occupancy Standard Error
    occupancy_upper = results$upper_occ,           # Occupancy Upper
    occupancy_lower = results$lower_occ,           # Occupancy Lower
    detection_estimate = results$Predicted_det,    # Detection Predictions
    detection_SE = results$SE_det,                 # Occupancy Standard Error
    detection_upper = results$upper_det,           # Occupancy Upper
    detection_lower = results$lower_det,           # Occupancy Lower
    formula = as.character(occ_model@formula)[3],  # Model Formula (Original Model vs. Random Model)
    AIC = occ_model@AIC,                           # Akaike Information Criterion
    GodnessOfFit_chi.square = Gof$chi.square,      # Godness Of Fit-Statistic
    GodnessOfFit_p.value = Gof$p.value,            # Godness Of Fit-Statistic
    GodnessOfFit_c.hat = Gof$c.hat.est             # Godness Of Fit-Statistic
  )
}

# Clean environment - remove unnecessary data
rm(species, occ_model, occ_random_model, occ_estimate, det_estimate, detection_history, detection_matrix, Gof, results, sitecovs, speciesData, umf, aic_occ_model, aic_occ_random_model )

 ## Random-Effects perform better for All species of Interest
 ## Models without Area, Matrix Type and Connectiovity perform better for all species of Interest 

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

## Occupancy List
occupancy_det_results <- list()

# Loop through each species with Covariates - TreeDensity and Vegetation Cover
for (species in species_of_interest) {
  
  # filter for Species
  speciesData <- det_hist_full %>%
    filter(scientificName %in% species)
  
  # create Site specific Covariate data
  sitecovs <- speciesData %>%
    select(locationName, PatchID, TreeDensity, VegetationCover, cameraHeight, Matrixs, Area, min_distance_to_next_patch_km)
  
  # Extract the detection history (columns with dates) 
  detection_history <- speciesData[, grep("^2024", names(speciesData))]
  
  # Assign row names using the locationName column
  row.names(detection_history) <- speciesData$locationName
  
  # Convert the detection history into the correct format for unmarked
  detection_matrix <- as.matrix(detection_history)
  
  # Convert character columns to factors
  sitecovs <- sitecovs %>%
    mutate_if(is.character, as.factor)
  
  # Check structure of Covariate Table
  str(sitecovs)
  
  # Create the unmarked frame for the current row
  umf <- unmarkedFrameOccu(y = detection_matrix, siteCovs = sitecovs, obsCovs = obscovs)
  
  # Fit the models 
  occ_model <- occu(~ (1 | weeks) ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_cov_TD <- occu(~ (1 | weeks) + TreeDensity ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_cov_VC <- occu(~ (1 | weeks) + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_cov_TD_VC <- occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  
  # Access AIC values only for valid models (converged and no NA in coefficients or summary)
  aic_occ_model <- ifelse(valid_model_check(occ_model), occ_model@AIC, NA)
  aic_occ_cov_TD <- ifelse(valid_model_check(occ_cov_TD), occ_cov_TD@AIC, NA)
  aic_occ_cov_VC <- ifelse(valid_model_check(occ_cov_VC), occ_cov_VC@AIC, NA)
  aic_occ_cov_TD_VC <- ifelse(valid_model_check(occ_cov_TD_VC), occ_cov_TD_VC@AIC, NA)
  
  # Print species name
  print(species)
  
  # Compare the AIC values and store the model with the lowest AIC (ignore models with NA AIC)
  aic_values <- c(aic_occ_model, aic_occ_cov_TD, aic_occ_cov_VC, aic_occ_cov_TD_VC)
  model_names <- c("Original Model", "TreeDensity Model", "VegetationCover Model", "TreeDensity + VegetationCover Model")
  
  # Find the model with the lowest valid AIC (excluding NA values)
  if (all(is.na(aic_values))) {
    cat("None of the models are valid or converged properly.\n")
  } else {
    # Get the lowest AIC among valid models
    lowest_aic <- min(aic_values, na.rm = TRUE)
    selected_model_index <- which(aic_values == lowest_aic)
    selected_model <- list(occ_model, occ_cov_TD, occ_cov_VC, occ_cov_TD_VC)[[selected_model_index]]
    cat(model_names[selected_model_index], "selected due to lowest AIC:", lowest_aic, "\n")
  }
  
  # Summary of the selected model
  summary(selected_model)
  
  # random Effects of selected model
  # ranef(selected_model)
  
  # Get the Occupancy Estimates
  occ_estimate <-  predict(selected_model, 
                           newdata = umf,
                           type = "state")
  
  # Add Identifiers to Occupancy estimates
  row.names(occ_estimate) <- sitecovs$locationName
  occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
  
  # Get Detection Estimates
  det_estimate <- predict(selected_model, 
                          newdata = umf,
                          type = "det")
  
  # remove duplicates
  det_estimate <- det_estimate[seq(11, nrow(det_estimate), by = 11), ] ## correct??
  
  # Add Identifiers to Occupancy estimates
  row.names(det_estimate) <- sitecovs$locationName
  det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
  
  # Godness of Fit-Test
  Gof <- mb.gof.test(selected_model, nsim = 100)
  Gof
  
  # Ensure that both estimates are correctly combined
  #results <- occ_estimate
  results <- merge(occ_estimate, det_estimate, by = "Row.names", suffixes = c("_occ", "_det"))
  
  # Create the final data frame
  occupancy_det_results[[species]] <- data.frame(
    Patch_ID = results$PatchID_occ,              # PatchID (PatchID)
    Camera = results$Row.names,                         # Camera ID (locationName)
    occupancy_estimate = results$Predicted_occ,         # Occupancy Predictions
    occupancy_SE = results$SE_occ,                      # Occupancy Standard Error
    occupancy_upper = results$upper_occ,                # Occupancy Upper
    occupancy_lower = results$lower_occ,                # Occupancy Lower
    detection_estimate = results$Predicted_det,         # Detection Predictions
    detection_SE = results$SE_det,                      # Occupancy Standard Error
    detection_upper = results$upper_det,                # Occupancy Upper
    detection_lower = results$lower_det,                # Occupancy Lower
    formula = as.character(selected_model@formula)[2],  # Covariates
    AIC = selected_model@AIC,                           # Akaike Information Criterion
    GodnessOfFit_chi.square = Gof$chi.square,           # Chi-Square  Godness Of Fit-Statistic
    GodnessOfFit_p.value = Gof$p.value,                 # P-Value     Godness Of Fit-Statistic
    GodnessOfFit_c.hat = Gof$c.hat.est                  # C-Hat       Godness Of Fit-Statistic
  )
}

# Clean environment - remove unnecessary data
rm(aic_occ_cov_TD, aic_occ_cov_TD_VC, aic_occ_cov_VC, aic_occ_model, lowest_aic, species, aic_values, model_names, selected_model_index)
rm(det_estimate, detection_history, detection_matrix, Gof, occ_cov_TD, occ_cov_TD_VC, occ_cov_VC, occ_estimate, occ_model, results, selected_model, sitecovs, speciesData, umf)


## Check for weeks - yes or no. And effect sizes of coefficients 
## Occupancy List
occupancy_det_results <- list()

# Loop through each species with Covariates - TreeDensity and Vegetation Cover
for (species in species_of_interest) {
  
  # filter for Species
  speciesData <- det_hist_full %>%
    filter(scientificName %in% species)
  
  # create Site specific Covariate data
  sitecovs <- speciesData %>%
    select(locationName, PatchID, TreeDensity, VegetationCover, cameraHeight, Matrix, Area, min_distance_to_next_patch_km)
  
  # Extract the detection history (columns with dates) 
  detection_history <- speciesData[, grep("^2024", names(speciesData))]
  
  # Assign row names using the locationName column
  row.names(detection_history) <- speciesData$locationName
  
  # Convert the detection history into the correct format for unmarked
  detection_matrix <- as.matrix(detection_history)
  
  # Convert character columns to factors
  sitecovs <- sitecovs %>%
    mutate_if(is.character, as.factor)
  
  # Check structure of Covariate Table
  str(sitecovs)
  
  # Create the unmarked frame for the current row
  umf <- unmarkedFrameOccu(y = detection_matrix, siteCovs = sitecovs, obsCovs = obscovs)
  
  # Fit the models 
  occ_model <- occu(~ TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_cov_TD <- occu(~ TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))
  occ_cov_VC <- occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID), data = umf, control = list(maxit = 1000))
  occ_cov_TD_VC <- occu(~ (1 | weeks) + TreeDensity + VegetationCover ~ (1 | PatchID) + Area + Matrix + min_distance_to_next_patch_km, data = umf, control = list(maxit = 1000))
  
  # Access AIC values only for valid models (converged and no NA in coefficients or summary)
  aic_occ_model <- ifelse(valid_model_check(occ_model), occ_model@AIC, NA)
  aic_occ_cov_TD <- ifelse(valid_model_check(occ_cov_TD), occ_cov_TD@AIC, NA)
  aic_occ_cov_VC <- ifelse(valid_model_check(occ_cov_VC), occ_cov_VC@AIC, NA)
  aic_occ_cov_TD_VC <- ifelse(valid_model_check(occ_cov_TD_VC), occ_cov_TD_VC@AIC, NA)
  
  # Print species name
  print(species)
  
  # Compare the AIC values and store the model with the lowest AIC (ignore models with NA AIC)
  aic_values <- c(aic_occ_model, aic_occ_cov_TD, aic_occ_cov_VC, aic_occ_cov_TD_VC)
  model_names <- c("Original Model", "TreeDensity Model", "VegetationCover Model", "TreeDensity + VegetationCover Model")
  
  # Find the model with the lowest valid AIC (excluding NA values)
  if (all(is.na(aic_values))) {
    cat("None of the models are valid or converged properly.\n")
  } else {
    # Get the lowest AIC among valid models
    lowest_aic <- min(aic_values, na.rm = TRUE)
    selected_model_index <- which(aic_values == lowest_aic)
    selected_model <- list(occ_model, occ_cov_TD, occ_cov_VC, occ_cov_TD_VC)[[selected_model_index]]
    cat(model_names[selected_model_index], "selected due to lowest AIC:", lowest_aic, "\n")
  }
  
  # Summary of the selected model
  summary(selected_model)
  
  # random Effects of selected model
  # ranef(selected_model)
  
  # Get the Occupancy Estimates
  occ_estimate <-  predict(selected_model, 
                           newdata = umf,
                           type = "state")
  
  # Add Identifiers to Occupancy estimates
  row.names(occ_estimate) <- sitecovs$locationName
  occ_estimate <- merge(occ_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
  
  # Get Detection Estimates
  det_estimate <- predict(selected_model, 
                          newdata = umf,
                          type = "det")
  
  # remove duplicates
  det_estimate <- det_estimate[seq(9, nrow(det_estimate), by = 9), ] ## correct??
  
  # Add Identifiers to Occupancy estimates
  row.names(det_estimate) <- sitecovs$locationName
  det_estimate <- merge(det_estimate, sitecovs[, c("locationName", "PatchID")], by.x = "row.names", by.y = "locationName", all.x = TRUE)
  
  # Godness of Fit-Test
  Gof <- mb.gof.test(selected_model, nsim = 100)
  Gof
  
  # Ensure that both estimates are correctly combined
  #results <- occ_estimate
  results <- merge(occ_estimate, det_estimate, by = "Row.names", suffixes = c("_occ", "_det"))
  
  # Create the final data frame
  occupancy_det_results[[species]] <- data.frame(
    Patch_ID = results$PatchID_occ,              # PatchID (PatchID)
    Camera = results$Row.names,                         # Camera ID (locationName)
    occupancy_estimate = results$Predicted_occ,         # Occupancy Predictions
    occupancy_SE = results$SE_occ,                      # Occupancy Standard Error
    occupancy_upper = results$upper_occ,                # Occupancy Upper
    occupancy_lower = results$lower_occ,                # Occupancy Lower
    detection_estimate = results$Predicted_det,         # Detection Predictions
    detection_SE = results$SE_det,                      # Occupancy Standard Error
    detection_upper = results$upper_det,                # Occupancy Upper
    detection_lower = results$lower_det,                # Occupancy Lower
    formula = as.character(selected_model@formula)[4],  # Covariates
    AIC = selected_model@AIC,                           # Akaike Information Criterion
    GodnessOfFit_chi.square = Gof$chi.square,           # Chi-Square  Godness Of Fit-Statistic
    GodnessOfFit_p.value = Gof$p.value,                 # P-Value     Godness Of Fit-Statistic
    GodnessOfFit_c.hat = Gof$c.hat.est                  # C-Hat       Godness Of Fit-Statistic
  )
}

# Clean environment - remove unnecessary data
rm(aic_occ_cov_TD, aic_occ_cov_TD_VC, aic_occ_cov_VC, aic_occ_model, lowest_aic, species, aic_values, model_names, selected_model_index)
rm(det_estimate, detection_history, detection_matrix, Gof, occ_cov_TD, occ_cov_TD_VC, occ_cov_VC, occ_estimate, occ_model, results, selected_model, sitecovs, speciesData, umf)


# Save Workspace 
# save.image("~/Desktop/Master Thesis/Occupancy.RData")

# Clean Environment
rm(list = ls())

# Reaload/Load required R-packages 
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load Workspace
load("~/Desktop/Master Thesis/Occupancy.RData")


# Unpack Results (List)
# List of species in occupancy_det_results

# Initialize an empty dataframe
combined_df <- data.frame()

# Loop through each species' dataframe and add a species column
for (species in species_of_interest) {
  # Extract the dataframe for each species
  df <- occupancy_det_results[[species]]
  
  # Add a new column 'Species' with the current species name
  df$Species <- species
  
  # Combine each species' dataframe into a single dataframe
  combined_df <- rbind(combined_df, df)
}

# Filter for relevant columns
filtered_df <- combined_df %>%
  select(Patch_ID, Species, occupancy_estimate, detection_estimate) 

# Group by Patch_ID and Species, summarizing occupancy estimates
grouped_df <- filtered_df %>%
  group_by(Patch_ID, Species) %>%
  summarise(occupancy_estimate = mean(occupancy_estimate), .groups = 'drop')

# Spread to wide format
wide_df <- grouped_df %>%
  pivot_wider(names_from = Species, values_from = occupancy_estimate, values_fill = 0)

# Display the resulting wide dataframe
head(wide_df)

# remove unnecessarily Data
rm(combined_df, df, filtered_df, grouped_df)

# Calculate community metrics
wide_df$C_capreolus <- round(wide_df$`Capreolus capreolus` * 100)
wide_df$F_silvestris <- round(wide_df$`Felis silvestris` * 100)
wide_df$Martes <- round(wide_df$Martes * 100)
wide_df$M_meles <- round(wide_df$`Meles meles` * 100)
wide_df$P_lotor <- round(wide_df$`Procyon lotor` * 100)
wide_df$S_scrufa <- round(wide_df$`Sus scrofa` * 100)
wide_df$V_vuples <- round(wide_df$`Vulpes vulpes` * 100)


# create subset for merging
det_hist_selected <- det_hist_full %>%
# mutate(PatchID = as.factor(PatchID)) %>%
  select(PatchID, Area, Perimeter, min_distance_to_next_patch_km, Matrixs) %>%
  distinct() # Ensure only unique values are joined

# Join the selected columns to the community table using Patch_ID
community <- wide_df %>%
  left_join(det_hist_selected, by = c("Patch_ID" = "PatchID"))

# remove unnecessarily Data
rm(det_hist_selected)


## Species Analysis ## 
# Binomial Logistic Regression (glm) on Capreolus Capreolus:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Cc.0 <- glm(cbind(community$C_capreolus, 100 - community$C_capreolus) ~ Matrix, family= binomial, data = community)
glm.Cc.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Cc.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Cc.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Cc.0, glm.Cc.1, glm.Cc.2, glm.Cc.3)


# Binomial Logistic Regression (glm) on Felis sylvestris:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Fs.0 <- glm(cbind(community$F_silvestris, 100 - community$F_silvestris) ~ Matrix, family= binomial, data = community)
glm.Fs.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Fs.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Fs.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Fs.0, glm.Fs.1, glm.Fs.2, glm.Fs.3)


# Binomial Logistic Regression (glm) on Martes:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.M.0 <- glm(cbind(community$Martes, 100 - community$Martes) ~ Matrix, family= binomial, data = community)
glm.M.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.M.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.M.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.M.0, glm.M.1, glm.M.2, glm.M.3)


# Binomial Logistic Regression (glm) on Meles meles:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Mm.0 <- glm(cbind(community$M_meles, 100 - community$M_meles) ~ Matrix, family= binomial, data = community)
glm.Mm.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Mm.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Mm.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Mm.0, glm.Mm.1, glm.Mm.2, glm.Mm.3)


# Binomial Logistic Regression (glm) on Procyon lotor:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Pl.0 <- glm(cbind(community$P_lotor, 100 - community$P_lotor) ~ Matrix, family= binomial, data = community)
glm.Pl.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Pl.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Pl.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Pl.0, glm.Pl.1, glm.Pl.2, glm.Pl.3)


# Binomial Logistic Regression (glm) on Sus scrufa:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Ss.0 <- glm(cbind(community$S_scrufa, 100 - community$S_scrufa) ~ Matrix, family= binomial, data = community)
glm.Ss.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Ss.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Ss.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.n.sp.0, glm.n.sp.1, glm.n.sp.2, glm.n.sp.3)


# Binomial Logistic Regression (glm) on Vulpes vulpes:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Vv.0 <- glm(cbind(community$V_vuples, 100 - community$V_vuples) ~ Matrix, family= binomial, data = community)
glm.Vv.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.Vv.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.Vv.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Vv.0, glm.Vv.1, glm.Vv.2, glm.Vv.3)


## Communty Analysis ## 
# Remove the 'Patch_ID' column and convert data to a matrix
occupancy_prob <- as.matrix(ifelse(wide_df[,-1] >= 0.8, 1, 0))
# occupancy_prob <- as.matrix(det_hist_full[,c(12:22)])

# Calculate Sørensen-like similarity using occupancy probabilities
# soerensen_sim_prob <- 1 - vegdist(as.matrix(wide_df[,-1]), method = "bray")

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~ log(Area) + Matrixs + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000)
adonis_result2 <- adonis2(occupancy_prob ~ log(Area) * Matrixs + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000)

# View the results
adonis_result
adonis_result2

# Adding significance codes based on p-values
adonis_result$significance <- cut(adonis_result$`Pr(>F)`,
                                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                   labels = c("***", "**", "*", ""))

# Plot barplot of R² values
ggplot(adonis_result, aes(x = rownames(adonis_result) , y = R2)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = significance), vjust = -0.5, size = 5) + 
  theme_minimal() +
  labs(title = "Effect Size (R²) from PERMANOVA",
       x = "Terms", y = "R²") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Calculate mean occupancy for each patch and categorize
community$occupancy_category <- cut(rowMeans(occupancy_prob),
                                    breaks = c(-Inf, 0.3, 0.7, Inf),
                                    labels = c("Low", "Mid", "High"))

# PERMANOVA for each occupancy category
adonis_high <- adonis2(occupancy_prob[community$occupancy_category == "High", ] ~ log(Area) + deploymentGroups + min_distance_to_next_patch_km,
                       data = community[community$occupancy_category == "High", ], method = "bray", permutations = 5000)

adonis_mid <- adonis2(occupancy_prob[community$occupancy_category == "Mid", ] ~ log(Area) + deploymentGroups + min_distance_to_next_patch_km,
                      data = community[community$occupancy_category == "Mid", ], method = "bray", permutations = 5000)

adonis_low <- adonis2(occupancy_prob[community$occupancy_category == "Low", ] ~ log(Area) + deploymentGroups + min_distance_to_next_patch_km,
                      data = community[community$occupancy_category == "Low", ], method = "bray", permutations = 5000)

# Compile results into a data frame for plotting
adonis_results <- data.frame(
  Term = rep(c("log(Area)", "deploymentGroups", "min_distance_to_next_patch_km"), 3),
  R2 = c(adonis_high$R2, adonis_mid$R2, adonis_low$R2),
  Category = rep(c("High", "Mid", "Low"), each = 3)
)

# Plotting R² values by category
ggplot(adonis_results, aes(x = Term, y = R2, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Effect Size (R²) from PERMANOVA by Occupancy Category",
       x = "Terms", y = "R²") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Pairwise Comparisons Using pairwise.adonis
# Custom pairwise.adonis function
pairwise.adonis <- function(x, factors, sim.method = 'bray', p.adjust.m = 'bonferroni') {
  library(vegan)
  
  co = combn(unique(factors), 2)
  results = data.frame()
  
  for (elem in 1:ncol(co)) {
    fac1 = co[1, elem]
    fac2 = co[2, elem]
    
    subset_x = x[factors %in% c(fac1, fac2), ]
    subset_factors = factors[factors %in% c(fac1, fac2)]
    
    ad = adonis(subset_x ~ subset_factors, method = sim.method)
    
    p_value = ad$aov.tab$`Pr(>F)`[1]
    
    results = rbind(results, data.frame(Comparison = paste(fac1, "vs", fac2),
                                        P.value = p_value))
  }
  
  # Adjust p-values
  results$P.adjusted = p.adjust(results$P.value, method = p.adjust.m)
  return(results)
}

# Usage Example:
# Assuming 'occupancy_prob' is your species matrix and 'community$Matrixs' are the factors
pairwise_results <- pairwise.adonis(occupancy_prob, community$Matrixs, sim.method = "bray", p.adjust.m = "bonferroni")

# View the pairwise results
pairwise_results

# Testing for Homogeneity of Group Dispersions Using betadisper
# Calculate the Sørensen similarity distance (Bray-Curtis)
dist_matrix <- vegdist(occupancy_prob, method = "bray")

# Test for homogeneity of group dispersions
betadisp_result <- betadisper(dist_matrix, community$Matrixs)

# View the betadisper result
summary(betadisp_result)

# Perform ANOVA to check for significance
anova_result <- anova(betadisp_result)

# View the ANOVA result for group dispersions
anova_result

# Plot the result
plot(betadisp_result)

# Perform Tukey HSD to identify which groups differ in their dispersions
tukey_result <- TukeyHSD(betadisp_result)
tukey_result



### END ###
