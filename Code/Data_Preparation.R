# Clear environment
rm(list = ls())
gc()

# Disable scientific notation
options(scipen = 999) 

# Load packages
library(lubridate)
library(tidyr)
library(dplyr)
library(sf)

## Import and prepare Data #### 
# Load Presence-Absenc Data
presence_absence_df <- st_read('/Users/nr72kini/Desktop/Master Thesis/R/Output/presence_absence_by_patch_new.gpkg') %>%
  st_drop_geometry()

# Select relevant columns from presence_absence_df
presence_absence_df <- presence_absence_df[, c("Patch_ID", "Perimeter", "Area", "CT", "SU_trunc", "min_distance_to_next_patch_km")]

# Load Covariates
covs <- read.csv('/Users/nr72kini/Desktop/Master Thesis/R/Data/covariate_estimation.csv')

# Load Deployment-Data
deployments <- read.csv('/Users/nr72kini/Desktop/Master Thesis/R/Data/Agouti Output Aug/deployments.csv')

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
observations <- read.csv('/Users/nr72kini/Desktop/Master Thesis/R/Data/Agouti Output Aug/observations.csv')

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

# Create all possible combinations to assess missing data
all_combinations <- expand.grid(
  deploymentID = unique(deployments$deploymentID),
  scientificName = unique(det_hist$scientificName)
)

# Subsample existing combinations
existing_combinations <- det_hist %>%
  select(deploymentID, scientificName) %>%
  distinct()

# Calculate missing combinations
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
rm(all_combinations, covs, deployments, det_hist, existing_combinations, missing_combinations, observations, presence_absence_df, obs_filtered)

# Filter Patches with wrong number of camera traps
det_hist_full <- det_hist_full %>%
  filter(CT == SU_trunc) %>%  # Keep only rows where SU matches SU_trunc
  select(-CT,-SU_trunc) %>%
  rename(Matrix = deploymentGroups) %>%
  rename(PatchID = deploymentTags) %>%
  mutate(log_Area = log(Area))


## Create Site-Covariates ####
# Create site-specific covariate data (outside the loop)
sitecovs <- det_hist_full %>%
  select(locationName, PatchID, TreeDensity, VegetationCover, cameraHeight, Matrix, 
         min_distance_to_next_patch_km, Area, log_Area) %>%
  mutate_if(is.character, as.factor) %>%
  distinct()  # Keep unique rows only

# Scale continuous variables in sitecovs
sitecovs$log_Area <- scale(sitecovs$log_Area)
sitecovs$min_distance_to_next_patch_km <- scale(sitecovs$min_distance_to_next_patch_km)

# Check sitecovs-Data
str(sitecovs)


## Create Observation Covariates ####

# Initialize a list to store observation covariates
obscovs <- list()

# Extract week names (columns 12 to 20 of det_hist_full)
week_names <- as.factor(names(det_hist_full)[12:20])

# Initialize a matrix with the correct number of rows (locationNames) and columns (week names)
weeks_matrix <- matrix(NA, nrow = length(unique(det_hist_full$locationName)), ncol = length(week_names))

# Assign row names to the matrix based on locationName
rownames(weeks_matrix) <- unique(det_hist_full$locationName)

# Assign column names to the matrix based on week names
colnames(weeks_matrix) <- week_names

# Populate the matrix with week numbers (1, 2, 3, ..., for each week)
for (i in 1:ncol(weeks_matrix)) {
  weeks_matrix[, i] <- i  # Each column will have the corresponding week number
}

# Convert the matrix to a factor matrix (each element will be a factor representing a week)
weeks_matrix <- apply(weeks_matrix, c(1, 2), factor)

# Store in the 'obscovs' list as the 'weeks' matrix
obscovs$weeks <- weeks_matrix

# Check the structure of the resulting factor matrix
str(obscovs$weeks)
head(obscovs$weeks)

# remove unnecessary Data
rm(week_names, weeks_matrix, i) 

# Save Workspace 
save.image("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Data-Prep.RData")
