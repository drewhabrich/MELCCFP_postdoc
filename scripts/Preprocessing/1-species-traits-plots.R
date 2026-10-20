## Species parameters visualizations
##
## Written by Carina Rauen Firkowski
## 14 April 2026
## 
## Description:
## This code loads the species parameters collected in the first round of the 
## literature review to visualize species traits.



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Load data --------------------------------------------------------------------

# Species parameters
spp_param <- read.csv(file.path(tabular_dir, "species-parameters.csv"))

# Mammals only
spp_mammals <- spp_param[spp_param$Class == "Mammalia",]



# Plot -------------------------------------------------------------------------

# Dispersal vs. Longevity
ggplot(data = spp_param, 
       aes(x = log(Dispersion_km), y = log(Longevite_ans), label = Abreviation)) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, size = 2) +
  labs(x = "log(Dispersion)", y = "log(Longévité)") +
  theme_classic(base_size = 14)

# Dispersal vs. Body size
ggplot(data = spp_param, 
       aes(x = log(Dispersion_km), y = log(Longueur_cm), label = Abreviation)) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, size = 2) +
  labs(x = "log(Dispersion)", y = "log(Longueur)") +
  theme_classic(base_size = 14)

# Dispersal vs. Home range
ggplot(data = spp_param, 
       aes(x = log10(Domaine_vital_ha), y = log10(Dispersion_km), label = Abreviation)) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, size = 2) +
  labs(x = "log10(Domaine vital)", y = "log10(Dispersion)") +
  theme_classic(base_size = 14)



# Validate ---------------------------------------------------------------------

# Calculate Pearson correlation
x = log10(spp_mammals$Domaine_vital_ha)
y = log10(spp_mammals$Dispersion_km)
r <- cor(x, y)
print(r)

# Fit Linear Model
model <- lm(y ~ x)
summary_model <- summary(model)

# Get intercept and slope
pvalue <- summary_model$coefficients[2,4]
slope <- coef(model)[2]
print(paste("Intercept:", intercept, "Slope:", slope))

# Get R-squared
r_squared <- summary_model$r.squared
print(paste("R-squared:", r_squared))

# Dispersal vs. Home range
ggplot(data = spp_mammals, 
       aes(x = log10(Domaine_vital_ha), y = log10(Dispersion_km), label = Abreviation)) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, size = 2) +
  labs(x = "log10(Domaine vital)", y = "log10(Dispersion)") +
  geom_smooth(method = "lm") +
  annotate("text", x = c(0.31, 0.5, 0.65, 0.6), y = c(2.75, 2.5, 2.25, 2), 
           label = c(
             paste0("r = ", round(r, 2)),
             paste0("Slope = ", round(slope, 2)),
             paste0("R-squared = ", round(r_squared, 2)),
             paste0("p-value = ", round(summary_model$coefficients[2,4], 3)))) +
  theme_classic(base_size = 14)



# K-means clustering -----------------------------------------------------------

# Identify clusters
result_BS <- kmeans(c(log(spp_param$Longueur_cm)), centers = 5, nstart = 25)
result_DV <- kmeans(log(spp_param$Domaine_vital_ha), centers = 5, nstart = 25)
result_D <- kmeans(log(spp_param$Dispersion_km), centers = 5, nstart = 25)

# Add results to species parameter table
spp_param$D <- result_D$cluster
spp_param$DV <- result_DV$cluster
spp_param$BS <- result_BS$cluster

# Save to file
write.csv(spp_param, file = file.path(tabular_dir, "spp_parameters_2.csv"))
