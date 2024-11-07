# Clean Environment
rm(list = ls())

# Reaload/Load required R-packages 
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load Workspace
# load("~/Desktop/Master Thesis/Occupancy.RData")
# load("~/Desktop/Master Thesis/Occupancy.RData_new")
load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy.RData")

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
  select(PatchID, Area, Perimeter, min_distance_to_next_patch_km, Matrix) %>%
  distinct() %>% # Ensure only unique values are joined %>%
  rename(Patch_ID = PatchID)

# Join the selected columns to the community table using Patch_ID
community <- wide_df %>%
  left_join(det_hist_selected, by = "Patch_ID")

# remove unnecessarily Data
rm(det_hist_selected)


## Species Analysis ## 
# Binomial Logistic Regression (glm) on Capreolus Capreolus:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Cc.0 <- glm(cbind(community$C_capreolus, 100) ~ Matrix, family= binomial, data = community)
glm.Cc.1 <- update(glm.Cc.0, .~. + log(Area))
glm.Cc.2 <- update(glm.Cc.1, .~. + Matrix : log(Area))
glm.Cc.3 <- update(glm.Cc.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Cc.0, glm.Cc.1, glm.Cc.2, glm.Cc.3)


# Binomial Logistic Regression (glm) on Felis sylvestris:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Fs.0 <- glm(cbind(community$F_silvestris, 100) ~ Matrix, family= binomial, data = community)
glm.Fs.1 <- update(glm.Fs.0, .~. + log(Area))
glm.Fs.2 <- update(glm.Fs.1, .~. + Matrix : log(Area))
glm.Fs.3 <- update(glm.Fs.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Fs.0, glm.Fs.1, glm.Fs.2, glm.Fs.3)


# Binomial Logistic Regression (glm) on Martes:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.M.0 <- glm(cbind(community$Martes, 100) ~ Matrix, family= binomial, data = community)
glm.M.1 <- update(glm.M.0, .~. + log(Area))
glm.M.2 <- update(glm.M.1, .~. + Matrix : log(Area))
glm.M.3 <- update(glm.M.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.M.0, glm.M.1, glm.M.2, glm.M.3)


# Binomial Logistic Regression (glm) on Meles meles:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Mm.0 <- glm(cbind(community$M_meles, 100) ~ Matrix, family= binomial, data = community)
glm.Mm.1 <- update(glm.Mm.0, .~. + log(Area))
glm.Mm.2 <- update(glm.Mm.1, .~. + Matrix : log(Area))
glm.Mm.3 <- update(glm.Mm.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Mm.0, glm.Mm.1, glm.Mm.2, glm.Mm.3)


# Binomial Logistic Regression (glm) on Procyon lotor:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Pl.0 <- glm(cbind(community$P_lotor, 100) ~ Matrix, family= binomial, data = community)
glm.Pl.1 <- update(glm.Pl.0, .~. + log(Area))
glm.Pl.2 <- update(glm.Pl.1, .~. + Matrix : log(Area))
glm.Pl.3 <- update(glm.Pl.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Pl.0, glm.Pl.1, glm.Pl.2, glm.Pl.3)


# Binomial Logistic Regression (glm) on Sus scrufa:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Ss.0 <- glm(cbind(community$S_scrufa, 100) ~ Matrix, family= binomial, data = community)
glm.Ss.1 <- update(glm.Ss.0, .~. + log(Area))
glm.Ss.2 <- update(glm.Ss.1, .~. + Matrix : log(Area))
glm.Ss.3 <- update(glm.Ss.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Ss.0, glm.Ss.1, glm.Ss.2, glm.Ss.3)


# Binomial Logistic Regression (glm) on Vulpes vulpes:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Vv.0 <- glm(cbind(community$V_vuples, 100) ~ Matrix, family= binomial, data = community)
glm.Vv.1 <- update(glm.Vv.0, .~. + log(Area))
glm.Vv.2 <- update(glm.Vv.1, .~. + Matrix : log(Area))
glm.Vv.3 <- update(glm.Vv.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Vv.0, glm.Vv.1, glm.Vv.2, glm.Vv.3)


## Communty Analysis ## 
# Remove the 'Patch_ID' column and convert data to a matrix
occupancy_prob <- as.matrix(ifelse(wide_df[,-1] >= 0.8, 1, 0))
# occupancy_prob <- as.matrix(det_hist_full[,c(12:22)])

# Calculate Sørensen-like similarity using occupancy probabilities
# soerensen_sim_prob <- 1 - vegdist(as.matrix(wide_df[,-1]), method = "bray")

# Run PERMANOVA on the species composition data
adonis_result <- adonis2(occupancy_prob ~ log(Area) + Matrix + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000)
adonis_result2 <- adonis2(occupancy_prob ~ log(Area) * Matrix + min_distance_to_next_patch_km, data = community, method = "bray", permutations = 5000)

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
community$occupancy_category <- cut(rowMeans(wide_df[,-1]),
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
# Assuming 'occupancy_prob' is your species matrix and 'community$Matrix' are the factors
pairwise_results <- pairwise.adonis(occupancy_prob, community$Matrix, sim.method = "bray", p.adjust.m = "bonferroni")

# View the pairwise results
pairwise_results

# Testing for Homogeneity of Group Dispersions Using betadisper
# Calculate the Sørensen similarity distance (Bray-Curtis)
dist_matrix <- vegdist(occupancy_prob, method = "bray")

# Test for homogeneity of group dispersions
betadisp_result <- betadisper(dist_matrix, community$Matrix)

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