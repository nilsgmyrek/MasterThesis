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
rm(model, models_of_interest)

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
rm(p_sitecovs, detection_columns)


## Rarefaction Curve ##
# occupancy_prob <- as.matrix(ifelse(wide_df[,c(2:8)] >= 0.7, 1, 0))
occupancy_prob <- as.matrix(wide_df[,c(2:8)])

# create rownames 
rownames(occupancy_prob) <- wide_df$Patch_ID

# Convert probabilities to estimated occurrences
expected_abundance <- round(occupancy_prob * 10)  # Adjust scaling factor as needed

# Create a rarefaction curve
rarecurve(expected_abundance, step = 1, col = "black", label = F)


# # Communty Analysis ##

# distance-based Redundancy Analysis (dbRDA)
# Calculate beta diversity (Bray-Curtis)
beta_dist <- vegdist(occupancy_prob, method = "bray")

# Visualize the distance matrix 
heatmap(as.matrix(beta_dist))


# dbRDA
dbRDA_result <- capscale(beta_dist ~ Matrix + log(Area) + min_distance_to_next_patch_km,
                         data = community)

# Summary and significance testing
summary(dbRDA_result)
anova(dbRDA_result)

anova_result <- anova(dbRDA_result, by = "margin")
print(anova_result)

plot(dbRDA_result, scaling = 2, display = "sites")

# Q-Q plot
qqnorm(residuals(dbRDA_result))
qqline(residuals(dbRDA_result))

# Residual plot to check for homoscedasticity
plot(fitted(dbRDA_result), residuals(dbRDA_result))


  
# Non-metric Multidimensional Scaling
# Run NMDS
#nmds_result <- metaMDS(sqrt(occupancy_prob), distance = "bray", k = 2)  # k=2 for 2D NMDS
nmds_result <- metaMDS(occupancy_prob, distance = "bray", k = 2)  # k=2 for 2D NMDS

# Extract scores
nmds_result 

# change categorical variable Matrix in binary variable
community <-  community %>%
  mutate(Matrix = ifelse(Matrix == "Urban Matrix", 0, 1))

# Create the environmental data frame
env_data <- as.data.frame(community[, c("Area", 
                                        "Matrix",
                                        "min_distance_to_next_patch_km")])

# Fit environmental variables using envfit
env_fit <- envfit(nmds_result, env_data, perm = 999)

# Print environmental fit results to check the significance of variables
print(env_fit)


## PERMANOVA ##
# betadisperser check
occ_dist <- vegdist(occupancy_prob, method = "bray")

Matrix_disp <- betadisper(occ_dist, group =community$Matrix, type = "centroid", sqrt.dist = T)
anova(Matrix_disp)
boxplot(Matrix_disp, xlab = "Matrix Type", ylab = "Distances to Centroids", names =c("Urban","Agriculture"))  # Boxplot of distances to group centroids

#plot(Matrix_disp, main = "Homogeneity of Multivariate Dispersion by Matrix")

community$Area_cat <- cut(log(community$Area), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, Inf), labels = c("0-1", "1-2", "2-3","3-4","4-5","5-6","6-7","7-8"))
Area_disp <- betadisper(occ_dist, group = community$Area_cat, type = "centroid")
#Area_disp <- betadisper(occ_dist, group = log(community$Area), type = "centroid")
anova(Area_disp)
boxplot(Area_disp)  # Boxplot of distances to group centroids

# Connectivity #
# Create the categorical variable for min_distance_to_next_patch_km
community$min_dist_cat <- cut(community$min_distance_to_next_patch_km,  
                              breaks = c(-Inf, 0, 0.25, 0.5, 0.75, Inf), 
                              labels = c("0", "0-0.25", "0.25-0.5", "0.5-0.75", ">=0.75"))

# Run betadisper with the newly created categorical variable
min_dist_disp <- betadisper(occ_dist, group = community$min_dist_cat, type = "centroid")

# Perform ANOVA to test for dispersion differences
anova_result <- anova(min_dist_disp)
print(anova_result)

# Visualize the dispersion with a boxplot
boxplot(min_dist_disp, main = "Dispersion of Min Distance to Next Patch Categories", 
        xlab = "Distance Categories", ylab = "Distances to Centroid", 
        col = "lightblue")


community$Matrix <- as.factor(community$Matrix)

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~ log(Area) + Matrix + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000, by = "margin")
# View the results
adonis_result
# Check AIC
AICc_permanova2(adonis_result)

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~  Matrix + min_distance_to_next_patch_km + log(Area), data = community, method = "bray", permutations = 5000, by = "margin")
# View the results
adonis_result
# Check AIC
AICc_permanova2(adonis_result)

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~  min_distance_to_next_patch_km + log(Area) + Matrix, data = community, method = "bray", permutations = 5000, by = "margin")
# View the results
adonis_result
# Check AIC
AICc_permanova2(adonis_result)


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
    x = "",  # Rename x-axis to "Classifications"
    y = "R²"
  ) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1))


# # Generalized dissimilarity model #
# Load the gdm package
library(gdm)

# Create the species data matrix
species_data <- as.data.frame(community[, c("Capreolus capreolus_Model1", 
                              "Felis silvestris_Model1", 
                              "Martes_Model1", 
                              "Meles meles_Model1", 
                              "Procyon lotor_Model4", 
                              "Sus scrofa_Model4", 
                              "Vulpes vulpes_Model1")])


species_data <- as.matrix(vegdist(species_data, method = "bray"))


site <- as.numeric(as.factor(community$Patch_ID))

species_data <- as.data.frame(cbind(site, species_data))

colnames(species_data)[1] <- ''

str(species_data)

# convert to presence absence data
#species_data[,c(2:8)] <- as.matrix(ifelse(species_data[,c(2:8)] >= 0.8, 1, 0))

# Create the environmental data frame
env_data <- as.data.frame(community[, c("Area", 
                                        "Matrix",
                                        "min_distance_to_next_patch_km",
                                        "long",
                                        "lat")])
 
env_data <- cbind(site, env_data)

# Format the data using formatsitepair()
gdm_data <- formatsitepair(
  bioData = species_data,              # Your species data
  bioFormat = 3,               # Site-by-species format
#  dist = "bray",               # Dissimilarity measure (can be "bray", "jaccard", etc.)
#  abundance = FALSE,           # Species presence/absence data
  siteColumn = "Patch_ID",     # Unique site identifier column
  predData = env_data,         # Your environmental predictor data
  XColumn = "long",
  YColumn = "lat"
  )

# Calc Variable importance
gdm.varImp(gdm_data, geo = F, splines = NULL, knots = NULL,
           predSelect = F, nPerm = 500, pValue=0.05, parallel = F, outFile = NULL)

# Fit the Generalized Dissimilarity Model (GDM)
gdm_model <- gdm(gdm_data, geo = F)

# Summarize the GDM
summary(gdm_model)

# plot
plot(gdm_model)

# beta <- vegdist(wide_df[,c(2:8)], method = "bray")
# beta <- as.matrix(beta)
# beta <- as.data.frame(beta)
# range(beta)  # Should be between 0 and 1
# #rownames(beta) <-  as.numeric(as.factor(rownames(beta)))
# rownames(beta) <-  species_data$Patch_ID
# #colnames(beta) <- species_data$Patch_ID
# beta$Patch_ID <- as.numeric(as.factor(rownames(beta)))
# beta <- as.matrix(beta)
# #beta <- vegdist(species_data[,c(2:8)])
# str(beta)
