


library(dplyr)
library(tidyr)
library(entropy)

spectral_df <- read.csv("E:/DATA/VI_Metrics.csv")

# ---- Step 0: Define entropy function ----
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# ---- Step 1: Identify VI columns ----
id_cols <- c("TileNumber", "SpeciesID", "TreeID")  # Adjust if other non-VI columns exist
vi_cols <- setdiff(names(spectral_df), id_cols)

# ---- Step 2: Convert to long format ----
spectral_long <- spectral_df %>%
  pivot_longer(
    cols = all_of(vi_cols),
    names_to = "Layer",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# ---- Step 3: Summarize by TreeID and Layer (VI) ----
vi_summary <- spectral_long %>%
  group_by(TreeID, SpeciesID, Layer) %>%
  summarise(
    Max = max(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Entropy = compute_entropy(Value),
    .groups = "drop"
  )

# Step 5: Pivot to wide format
vi_wide <- vi_summary %>%
  pivot_wider(
    names_from = Layer,
    values_from = c(Max, Min, Mean, SD, Entropy),
    names_glue = "{.value}_{Layer}"
  )

# ---- Optional: Write to CSV ----
write.csv(vi_wide, "E:/DATA/VI_Metrics_CanopySummary.csv", row.names = FALSE)
