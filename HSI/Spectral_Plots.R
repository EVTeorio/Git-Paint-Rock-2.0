

library(tidyverse)


sp_spectra <- read.csv("E:/DATA/HSI_Metrics.csv")


# ---- 1. Prepare spectral data ----
# Extract spectral column names and wavelengths
spectral_cols <- grep("^Mean_", names(sp_spectra), value = TRUE)
wavelengths <- as.numeric(gsub("Mean_|\\.nm", "", spectral_cols))

# Reshape to long format
long_df <- sp_spectra %>%
  select(TreeID, SpeciesID, all_of(spectral_cols)) %>%
  pivot_longer(
    cols = all_of(spectral_cols),
    names_to = "Wavelength",
    values_to = "Reflectance"
  ) %>%
  mutate(Wavelength = as.numeric(gsub("Mean_|\\.nm", "", Wavelength)))

# ---- 2. Create visible spectrum bands ----
visible_bands <- tibble(
  xmin = c(400, 495, 570),
  xmax = c(495, 570, 700),
  fill = c("blue", "green", "red")
)

# ---- 3. Loop and plot by species ----
species_list <- unique(long_df$SpeciesID)


for (sp in species_list) {
  sp_data <- long_df %>% filter(SpeciesID == sp)
  n_indiv <- sp_data %>% distinct(TreeID) %>% nrow()
  
  summary_df <- sp_data %>%
    group_by(Wavelength) %>%
    summarise(
      mean_reflectance = mean(Reflectance, na.rm = TRUE),
      p25 = quantile(Reflectance, 0.25, na.rm = TRUE),
      p75 = quantile(Reflectance, 0.75, na.rm = TRUE),
      p12.5 = quantile(Reflectance, 0.125, na.rm = TRUE),
      p87.5 = quantile(Reflectance, 0.875, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot() +
    # Visible bands
    geom_rect(data = visible_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.35, inherit.aes = FALSE) +
    scale_fill_manual(values = c("blue" = "#ADD8E6", "green" = "#90EE90", "red" = "#FFC0CB"), guide = "none") +
    
    # Light gray: 75% band
    geom_ribbon(data = summary_df, aes(x = Wavelength, ymin = p12.5, ymax = p87.5), fill = "lightgray", alpha = 0.5) +
    # Dark gray: 50% band
    geom_ribbon(data = summary_df, aes(x = Wavelength, ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.6) +
    
    # Mean line
    geom_line(data = summary_df, aes(x = Wavelength, y = mean_reflectance), color = "black", size = 1.1) +
    
    # Plot formatting
    labs(
      title = paste0("Spectral Signature: ", sp, " (n=", n_indiv, ")"),
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Save plot
  ggsave(filename = paste0("E:/DATA/Spectral_Plots/", sp, "_spectrum.png"),
         plot = p, width = 10, height = 6)
}
beep()
