

library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(tidyverse)
library(corrr)
library(igraph)
library(vegan)

# Load data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# Load and average importance results across iterations
importance_list <- readRDS("E:/DATA/Metric_Importance_Selection/Metric_Importance_Values.rds")

# Convert list of named vectors into a data frame and average by variable
importance_df <- as.data.frame(do.call(rbind, importance_list))
mean_importance <- colMeans(importance_df, na.rm = TRUE)

# Identify VI columns (exclude ID columns)
vi_cols <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

# Filter out known problematic TreeIDs
Canopies_noVIs <- c("080736","090800","090830","090848","090863",
                    "100987","110563","110602","110618","110752","110819",
                    "110883","110913","120866","140806","150720","80917")

samples_clean <- combined_df[!combined_df$TreeID %in% Canopies_noVIs, ]

# Remove constant or all-NA columns
vi_data_clean <- samples_clean %>%
  select(all_of(vi_cols)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# Compute correlation matrix
vi_corr_matrix <- cor(vi_data_clean, use = "pairwise.complete.obs", method = "pearson")

# Set correlation threshold
corr_threshold <- 0.7

# Keep only upper triangle
corr_matrix_upper <- vi_corr_matrix
corr_matrix_upper[lower.tri(corr_matrix_upper, diag = TRUE)] <- 0

# Get correlated variable pairs (edges)
edges <- which(abs(corr_matrix_upper) >= corr_threshold, arr.ind = TRUE)
edge_list <- data.frame(
  from = rownames(corr_matrix_upper)[edges[, 1]],
  to   = colnames(corr_matrix_upper)[edges[, 2]],
  corr = corr_matrix_upper[edges]
)

# Build correlation graph and find components
g <- graph_from_data_frame(edge_list, directed = FALSE)
components <- components(g)
groups <- split(names(components$membership), components$membership)

# Also include unconnected (uncorrelated) variables
all_vars <- colnames(vi_data_clean)
connected_vars <- unlist(groups)
unconnected_vars <- setdiff(all_vars, connected_vars)

# Select top variable per group using averaged importance
select_top_by_importance <- function(vars, importance_scores) {
  vars_with_scores <- importance_scores[names(importance_scores) %in% vars]
  if (length(vars_with_scores) == 0) return(NA)
  return(names(which.max(vars_with_scores)))
}

selected_vars <- sapply(groups, select_top_by_importance, importance_scores = mean_importance)
print(selected_vars)
