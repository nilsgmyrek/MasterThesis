# Clean Environment
rm(list = ls())

# Reaload/Load required R-packages 
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load Workspace
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occ_results.RData")

# Unpack Occupancy-Results
# List of species in occupancy_results - AIC
models_of_interest <- c("Capreolus capreolus_Model1", 
                        "Martes_Model1", 
                        "Felis silvestris_Model1", 
                        "Sus scrofa_Model4", 
                        "Procyon lotor_Model4", 
                        "Meles meles_Model1", 
                        "Vulpes vulpes_Model1")

# List of species in occupancy_results - GOF
# models_of_interest <- c("Capreolus capreolus_Model2", "Martes_Model1", "Felis silvestris_Model4", "Sus scrofa_Model3", "Procyon lotor_Model7", "Meles meles_Model3", "Vulpes vulpes_Model4" )

# Initialize an empty dataframe to store combined results
combined_df <- data.frame()

# Loop through each species' dataframe and add a species column
for (model in models_of_interest) {
  # Extract the dataframe for each species
  df <- occupancy_results[[model]]
  
  # Add a new column 'Species' with the current species name
  df$Model <- model
  
  # Combine each species' dataframe into a single dataframe
  combined_df <- rbind(combined_df, df)
}

# Filter for relevant columns
filtered_df <- combined_df %>%
  select(Patch_ID, Model, occupancy_estimate, detection_estimate) 

# Group by Patch_ID and Species, summarizing occupancy estimates
grouped_df <- filtered_df %>%
  group_by(Patch_ID, Model) %>%
  summarise(occupancy_estimate = mean(occupancy_estimate), .groups = 'drop')

# Convert to wide format for occupancy estimates by species
wide_df <- grouped_df %>%
  pivot_wider(names_from = Model, values_from = occupancy_estimate)

# View the final wide-format dataframe
print(wide_df)

# remove unnecessarily Data
rm(combined_df, df, filtered_df, grouped_df)

# create subset for merging
det_hist_selected <- det_hist_full %>%
  select(PatchID, Area, min_distance_to_next_patch_km, Matrix, long, lat) %>%
  distinct() %>% # Ensure only unique values are joined
  rename(Patch_ID = PatchID)

# Join the selected columns to the community table using Patch_ID
community <- wide_df %>%
  left_join(det_hist_selected, by = "Patch_ID")

# remove unnecessarily Data
rm(det_hist_selected)
  
# # Communty Analysis ##
# Remove the 'Patch_ID' column and convert data to a matrix
# occupancy_prob <- as.matrix(ifelse(wide_df[,c(2:8)] >= 0.7, 1, 0))
occupancy_prob <- as.matrix(wide_df[,c(2:8)])

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~ log(Area) + Matrix + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000)

# View the results
adonis_result

# Adding significance codes based on p-values
adonis_result$significance <- cut(adonis_result$`Pr(>F)`,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("***", "**", "*", ""))

# Plot barplot of R² values
ggplot(adonis_result, aes(x = rownames(adonis_result) , y = R2)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = significance), vjust = -0.5, size = 5) +
  theme_minimal() +
  labs(
    title = "Effect Size (R²) from PERMANOVA",
    x = "Classifications",  # Rename x-axis to "Classifications"
    y = "R²"
  ) +
  labs(title = "Effect Size (R²) from PERMANOVA",
       x = "Terms", y = "R²") +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1))



### END ###
#beta <- vegdist(as.matrix(wide_df[,c(2:8)]), method = "jaccard")
library(sf)
# Load Presence-Absenc Data
presence_absence_sf <- st_read('/Users/nr72kini/Desktop/Master Thesis/R/Output/presence_absence_by_patch_new.gpkg')
str(presence_absence_sf)

# Create a non-spatial version of the data
presence_absence_df <- st_drop_geometry(presence_absence_sf)
str(presence_absence_df)

# Filter Patches with wrong number of camera traps
presence_absence_df <- presence_absence_df %>%
  filter(CT == SU_trunc) %>%
  na.omit()

# Create the species data matrix
species_data <- presence_absence_df %>%
  select("Patch_ID", 
         "Capreolus.capreolus", 
         "Martes", 
         "Felis.silvestris", 
         "Sus.scrofa", 
         "Procyon.lotor", 
         "Meles.meles", 
         "Vulpes.vulpes")


# # Generalized dissimilarity model #
# Load the gdm package
library(gdm)

# Create the species data matrix
species_data <- as.data.frame(community[, c("Patch_ID", "Capreolus capreolus_Model1", 
                              "Felis silvestris_Model1", 
                              "Martes_Model1", 
                              "Meles meles_Model1", 
                              "Procyon lotor_Model4", 
                              "Sus scrofa_Model4", 
                              "Vulpes vulpes_Model1")])

str(species_data)

# convert to presence absence data
species_data[,c(2:8)] <- as.matrix(ifelse(species_data[,c(2:8)] >= 0.8, 1, 0))

# Create a binary variable for land use: 0 for agriculture, 1 for urban
community <- community %>%
  mutate(land_use_binary = case_when(
    Matrix == "Agricultural Matrix" ~ 0,  # Assuming 'Agricultural' represents agriculture
    Matrix == "Urban Matrix" ~ 1,         # Assuming 'Urban' represents urban
    TRUE ~ NA_real_                # Any other values will be assigned NA
  ))


# Create the environmental data frame
env_data <- as.data.frame(community[, c("Patch_ID", "Area", "min_distance_to_next_patch_km", "land_use_binary", "long", "lat")])

# Check the structures of both
str(species_data)
str(env_data)

# Format the data using formatsitepair()
gdm_data <- formatsitepair(
  bioData = species_data,         # Your species data
  bioFormat = 1,                  # Site-by-species format
  dist = "jaccard",               # Dissimilarity measure (can be "bray", "jaccard", etc.)
  abundance = FALSE,              # Species presence/absence data
  siteColumn = "Patch_ID",         # Unique site identifier column
  predData = env_data,             # Your environmental predictor data
  XColumn = "long",
  YColumn = "lat"
  )

# gdm.crossvalidation(gdm_data, train.proportion=0.9, n.crossvalid.tests=1,
#                    geo=FALSE, splines=NULL, knots=NULL)

gdm.varImp(gdm_data, geo = F, splines = NULL, knots = NULL,
           predSelect = F, nPerm = 500, pValue=0.05, parallel = F, outFile = NULL)

# Fit the Generalized Dissimilarity Model (GDM)
gdm_model <- gdm(gdm_data, geo = F)

# Summarize the GDM
summary(gdm_model)

# plot
plot(gdm_model, type = "spline", xlab = "Predictor Value", ylab = "Partial Dissimilarity")

plot(gdm_model)


  
  