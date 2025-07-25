library(party)

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
importance_list <- readRDS("E:/DATA/Perfomance/Performance_Model_All_Importance.rds")

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
##############################################################################
library(ggplot2)

# Convert to data frame and sort by importance
importance_df_plot <- data.frame(
  Variable = names(importance_df),
  Importance = mean_importance
) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 50)  # keep only top 50

# Plot top 50
ggplot(importance_df_plot, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top Variables by Average Importance",
    x = "609 Metrics",
    y = "Mean Importance Score"
  ) +
  theme_minimal(base_size = 13)


row.names(importance_df_plot)

########################################################################################

# Compute correlation matrix
vi_corr_matrix <- cor(vi_data_clean, use = "pairwise.complete.obs", method = "pearson")

# Set correlation threshold
corr_threshold <- 0.75
importance_threshold <- 0

# Build correlation graph
corr_matrix_upper <- vi_corr_matrix
corr_matrix_upper[lower.tri(corr_matrix_upper, diag = TRUE)] <- 0
edges <- which(abs(corr_matrix_upper) >= corr_threshold, arr.ind = TRUE)
edge_list <- data.frame(
  from = rownames(corr_matrix_upper)[edges[, 1]],
  to   = colnames(corr_matrix_upper)[edges[, 2]],
  corr = corr_matrix_upper[edges]
)

g <- graph_from_data_frame(edge_list, directed = FALSE)
components <- components(g)
groups <- split(names(components$membership), components$membership)

# Get all correlated and uncorrelated variables
all_vars <- colnames(vi_data_clean)
connected_vars <- unlist(groups)
unconnected_vars <- setdiff(all_vars, connected_vars)

# Select top variable in group
select_top_if_above_threshold <- function(vars, importance_scores, threshold) {
  vars_with_scores <- importance_scores[names(importance_scores) %in% vars]
  vars_above_threshold <- vars_with_scores[vars_with_scores >= threshold]
  if (length(vars_above_threshold) == 0) return(NA)
  return(names(which.max(vars_above_threshold)))
}

selected_group_vars <- sapply(groups, select_top_if_above_threshold,
                              importance_scores = mean_importance,
                              threshold = importance_threshold)

# Remove NA entries (groups without any variable ≥ threshold)
selected_group_vars <- selected_group_vars[!is.na(selected_group_vars)]
selected_metrics <- na.omit(selected_group_vars)
print(selected_metrics)





