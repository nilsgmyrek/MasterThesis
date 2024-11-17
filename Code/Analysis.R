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
#load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy1.RData")
#load("/Users/nr72kini/Desktop/Master Thesis/Github/MasterThesis/Data/Occupancy_all.RData")

# Create empty data frame to store results
results <- data.frame(
  species = character(),
  model = character(),
  AIC = numeric(),
  p_value = numeric(),
  c_hat = numeric(),
  stringsAsFactors = FALSE
)

# Combine AIC and gof results
for (species in names(aic_results)) {
  aic_data <- aic_results[[species]]
  gof_data <- gof_results[[species]]
  
  # Merge AIC and GOF data by model
  species_results <- merge(aic_data, gof_data, by = "model")
  species_results$species <- species  # Add species name
  
  # Append to the combined results
  results <- rbind(results, species_results)
}

# View combined results
print(results)


# Unpack Occupancy-Results

# List of species in occupancy_results - AIC
# models_of_interest <- c("Capreolus capreolus_Model1", "Martes_Model1", "Felis silvestris_Model1", "Sus scrofa_Model4", "Procyon lotor_Model4", "Meles meles_Model1", "Vulpes vulpes_Model1" )

# List of species in occupancy_results - GOF
models_of_interest <- c("Capreolus capreolus_Model2", "Martes_Model1", "Felis silvestris_Model4", "Sus scrofa_Model3", "Procyon lotor_Model7", "Meles meles_Model3", "Vulpes vulpes_Model4" )

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

# Calculate community metrics
wide_df <- wide_df %>%
  mutate(C_capreolus = round(wide_df[[2]] * 100)) %>%
  mutate(F_silvestris = round(wide_df[[3]] * 100)) %>%
  mutate(Martes = round(wide_df[[4]] * 100)) %>%
  mutate(M_meles = round(wide_df[[5]] * 100)) %>%
  mutate(P_lotor = round(wide_df[[6]] * 100)) %>%
  mutate(S_scrufa = round(wide_df[[7]] * 100)) %>%
  mutate(V_vuples = round(wide_df[[8]] * 100))

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
glm.Cc.0 <- glm(cbind(community$C_capreolus, 100 - community$C_capreolus) ~ Matrix, family= binomial, data = community)
glm.Cc.1 <- update(glm.Cc.0, .~. + log(Area))
glm.Cc.2 <- update(glm.Cc.1, .~. + Matrix : log(Area))
glm.Cc.3 <- update(glm.Cc.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Cc.0, glm.Cc.1, glm.Cc.2, glm.Cc.3)

coef(glm.Cc.2)

# Binomial Logistic Regression (glm) on Felis sylvestris:
glm.Fs.0 <- glm(cbind(community$F_silvestris, 100 - community$F_silvestris) ~ Matrix, family= binomial, data = community)
glm.Fs.1 <- update(glm.Fs.0, .~. + log(Area))
glm.Fs.2 <- update(glm.Fs.1, .~. + Matrix : log(Area))
glm.Fs.3 <- update(glm.Fs.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Fs.0, glm.Fs.1, glm.Fs.2, glm.Fs.3)

coef(glm.Fs.3)

# Binomial Logistic Regression (glm) on Martes:
glm.M.0 <- glm(cbind(community$Martes, 100 - community$Martes) ~ Matrix, family= binomial, data = community)
glm.M.1 <- update(glm.M.0, .~. + log(Area))
glm.M.2 <- update(glm.M.1, .~. + Matrix : log(Area))
glm.M.3 <- update(glm.M.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.M.0, glm.M.1, glm.M.2, glm.M.3)

coef(glm.M.1)

# Binomial Logistic Regression (glm) on Meles meles:
glm.Mm.0 <- glm(cbind(community$M_meles, 100 - community$M_meles) ~ Matrix, family= binomial, data = community)
glm.Mm.1 <- update(glm.Mm.0, .~. + log(Area))
glm.Mm.2 <- update(glm.Mm.1, .~. + Matrix : log(Area))
glm.Mm.3 <- update(glm.Mm.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Mm.0, glm.Mm.1, glm.Mm.2, glm.Mm.3)

coef(glm.Mm.1)

# Binomial Logistic Regression (glm) on Procyon lotor:
glm.Pl.0 <- glm(cbind(community$P_lotor, 100 - community$P_lotor) ~ Matrix, family= binomial, data = community)
glm.Pl.1 <- update(glm.Pl.0, .~. + log(Area))
glm.Pl.2 <- update(glm.Pl.1, .~. + Matrix : log(Area))
glm.Pl.3 <- update(glm.Pl.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Pl.0, glm.Pl.1, glm.Pl.2, glm.Pl.3)

coef(glm.Pl.3)

# Binomial Logistic Regression (glm) on Sus scrufa:
#community$miss_percentage <- (100 - community$C_capreolus) 
glm.Ss.0 <- glm(cbind(community$S_scrufa, 100- community$S_scrufa) ~ Matrix, family= binomial, data = community)
glm.Ss.1 <- update(glm.Ss.0, .~. + log(Area))
glm.Ss.2 <- update(glm.Ss.1, .~. + Matrix : log(Area))
glm.Ss.3 <- update(glm.Ss.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Ss.0, glm.Ss.1, glm.Ss.2, glm.Ss.3)

coef(glm.Ss.1)

# Binomial Logistic Regression (glm) on Vulpes vulpes:
glm.Vv.0 <- glm(cbind(community$V_vuples, 100- community$V_vuples) ~ Matrix, family= binomial, data = community)
glm.Vv.1 <- update(glm.Vv.0, .~. + log(Area))
glm.Vv.2 <- update(glm.Vv.1, .~. + Matrix : log(Area))
glm.Vv.3 <- update(glm.Vv.1, .~. + min_distance_to_next_patch_km)

# Check for Model Performance
AIC(glm.Vv.0, glm.Vv.1, glm.Vv.2, glm.Vv.3)

coef(glm.Vv.1)

## Communty Analysis ## 
# Remove the 'Patch_ID' column and convert data to a matrix
occupancy_prob <- as.matrix(ifelse(wide_df[,c(2:8)] >= 0.7, 1, 0))

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
  labs(title = "Effect Size (R²) from PERMANOVA",
       x = "Terms", y = "R²") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### END ###


























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

