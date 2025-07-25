

# ===============================
# Load required libraries
library(ropls)
library(ggplot2)
library(dplyr)

df <- read.csv("E:/DATA/HSI_Metrics.csv")

# Filter out rows with missing SpeciesID (if present)
df <- df %>% filter(!is.na(SpeciesID))

# Ensure SpeciesID is a factor
df$SpeciesID <- as.factor(df$SpeciesID)

# Step 2: Select predictor columns (e.g., Max_*.nm)
spectra_cols <- grep("^Mean_", names(df), value = TRUE)

# Step 3: Create predictor matrix X and response Y
X <- df[, spectra_cols]
Y <- df$SpeciesID


# Run PLS-DA instead of O-PLS-DA
plsda_model <- opls(
  X,
  Y,
  predI = 2,           # You can increase predictive components if needed
  orthoI = 0,          # 0 orthogonal components = PLS-DA
  scaleC = "standard"
)


# Step 5: Extract VIP Scores
vip_scores <- getVipVn(plsda_model)

# Step 6: Extract Wavelengths from column names
wavelengths <- as.numeric(gsub("Mean_|\\.nm", "", colnames(X)))

# Step 7: Prepare Data for Plotting
vip_df <- data.frame(
  Wavelength = wavelengths,
  VIP = vip_scores
)

# Step 8: Plot VIP Scores
ggplot(vip_df, aes(x = Wavelength, y = VIP)) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "VIP Scores Across Spectral Range (PLS-DA)",
    x = "Wavelength (nm)",
    y = "VIP Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
