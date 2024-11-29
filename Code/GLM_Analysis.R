rm(list =ls())
gc()

# Analysis
library(sf)
library(dplyr)
library(car)

# Disable scientific notation
options(scipen = 999)

# Load Presence-Absenc Data
presence_absence_sf <- st_read('/Users/nr72kini/Desktop/Master Thesis/R/Output/presence_absence_by_patch_new.gpkg')
str(presence_absence_sf)

# Create a non-spatial version of the data
presence_absence_df <- st_drop_geometry(presence_absence_sf)
str(presence_absence_df)

# Check for missing values
sum(is.na(presence_absence_df))

# Filter out the rows with missing data: CT-Stolen
presence_absence_df <- na.omit(presence_absence_df)

# Filter Patches with wrong number of camera traps
presence_absence_df <- presence_absence_df %>%
  filter(CT == SU_trunc)             # Keep only rows where SU matches SU_trunc

# Binomial Logistic Regression (glm) on Species Richness:
presence_absence_df$miss_species <- (7 - presence_absence_df$Species_Richness) 
glm.n.sp.0 <- glm(cbind(presence_absence_df$Species_Richness, presence_absence_df$miss_species) ~ Matrix, family=binomial, data = presence_absence_df)
glm.n.sp.1 <- update(glm.n.sp.0, .~. + log(Area))
glm.n.sp.2 <- update(glm.n.sp.1, .~. + Matrix : log(Area))
glm.n.sp.3 <- update(glm.n.sp.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.n.sp.0, glm.n.sp.1, glm.n.sp.2, glm.n.sp.3)

coef(glm.n.sp.1)

# Summary of best performing Model
summary(glm.n.sp.1)

vif(glm.n.sp.3)

# Plot best Model
par(mfrow=c(2,2))
plot(glm.n.sp.3)

# Count non-zero species richness
non_zero_species <- sum(presence_absence_df$Species_Richness > 0)

# Number of predictors in glm.n.sp.1
num_predictors <- length(coef(glm.n.sp.1)) - 1  # Exclude the intercept

# Check ratio
events_per_predictor <- non_zero_species / num_predictors

# Output result
cat("Events per predictor:", events_per_predictor, "\n")
if (events_per_predictor < 10) {
  cat("Warning: Sample size may be insufficient for reliable estimates.\n")
} else {
  cat("Sample size is adequate.\n")
}


## Occ
occ.dat <- read.csv("/Users/nr72kini/Desktop/Master Thesis/R/Output/Occupancy_Species_vs_Patch.csv")
