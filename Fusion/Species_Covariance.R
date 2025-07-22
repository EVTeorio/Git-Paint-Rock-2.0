

library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# ---- Load data ----
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# ---- Identify VI columns ----
id_cols <- c("TreeID", "SpeciesID")
vi_cols <- setdiff(names(combined_df), id_cols)

# ---- Function to compute covariance matrix for each species ----
species_cov_matrices <- combined_df %>%
  group_by(SpeciesID) %>%
  group_split() %>%
  map(function(df) {
    species <- unique(df$SpeciesID)
    data <- dplyr::select(df, all_of(vi_cols))  # full namespace
    cov_matrix <- cov(data, use = "pairwise.complete.obs")
    list(SpeciesID = species, CovMatrix = cov_matrix)
  })

# ---- Example: View covariance matrix for one species ----
example_cov <- species_cov_matrices[[1]]


# ---- Function to flatten matrix into vector (upper triangle only) ----
flatten_cov <- function(mat) {
  mat[upper.tri(mat, diag = TRUE)]
}

# ---- Create a named list of flattened matrices ----
cov_list <- map(species_cov_matrices, ~ flatten_cov(.x$CovMatrix))
names(cov_list) <- map_chr(species_cov_matrices, "SpeciesID")

# ---- Compute pairwise matrix distances ----
species_pairs <- combn(names(cov_list), 2, simplify = FALSE)

matrix_distances <- map_df(species_pairs, function(pair) {
  mat1 <- cov_list[[pair[1]]]
  mat2 <- cov_list[[pair[2]]]
  dist <- sqrt(sum((mat1 - mat2)^2))  # Frobenius norm
  tibble(Species1 = pair[1], Species2 = pair[2], Frobenius_Distance = dist)
})

# ---- View matrix distances between species ----
print(matrix_distances)

# Convert covariance matrix to long format for heatmap
cov_df <- as.data.frame(as.table(example_cov$CovMatrix))
colnames(cov_df) <- c("VI1", "VI2", "Covariance")

ggplot(cov_df, aes(VI1, VI2, fill = Covariance)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_minimal() +
  ggtitle(paste("Covariance Matrix for", example_cov$SpeciesID))
