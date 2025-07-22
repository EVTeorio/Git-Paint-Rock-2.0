

library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)

# Load data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# Identify VI columns (exclude ID columns)
vi_cols <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

# Remove constant or all-NA columns
vi_data_clean <- combined_df %>%
  select(all_of(vi_cols)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# Correlation matrix
vi_corr_matrix <- cor(vi_data_clean, use = "pairwise.complete.obs", method = "pearson")

# View result
print(round(vi_corr_matrix, 3))


# ---- Optional: Heatmap using corrplot ----
corrplot::corrplot(
  vi_corr_matrix,
  method = "color",
  type = "upper",
  tl.cex = 0.7,
  number.cex = 0.6,
  tl.col = "black",
  order = "hclust",
  addCoef.col = "black"
)

beep()
