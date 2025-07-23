
# Load required packages
library(readr)
library(dplyr)
library(partykit)
library(VIM)      # For KNN imputation
library(tibble)   # For rownames_to_column if needed

# Step 1: Load data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# Step 2: Define TreeIDs to exclude
Canopies_noVIs <- c("080736","090800","090830","090848","090863",
                    "100987","110563","110602","110618","110752","110819",
                    "110883","110913","120866","140806","150720","80917")

# Step 3: Remove excluded TreeIDs
samples_clean <- combined_df[!combined_df$TreeID %in% Canopies_noVIs, ]

# Step 4: Identify metric columns (excluding ID columns)
all_metrics <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

# Step 5: Filter constant/all-NA columns
vi_filtered <- samples_clean %>%
  select(all_of(all_metrics)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

low_var_cols <- sapply(vi_filtered, function(x) length(unique(na.omit(x))) < 2)
vi_filtered <- vi_filtered[, !low_var_cols]

# Step 6: Recombine with IDs
canopy_means <- samples_clean %>%
  select(TreeID, SpeciesID) %>%
  bind_cols(vi_filtered)

# Step 7: Ensure SpeciesID is a factor (target variable)
canopy_means$SpeciesID <- as.factor(canopy_means$SpeciesID)

# Step 8: KNN Imputation (excluding TreeID)
# Note: kNN automatically adds suffixes like "_imp" — we remove them afterward
imputed_data <- kNN(canopy_means %>% select(-TreeID), 
                    variable = names(vi_filtered), 
                    k = 5, 
                    imp_var = FALSE)
beep()

# Step 9: Recombine TreeID (if needed)
imputed_data$TreeID <- canopy_means$TreeID

# Convert to formula
predictors <- imputed_data %>%
  select(-TreeID) %>%
  select(-SpeciesID) %>%
  colnames()

# Step 10: Fit cforest model
set.seed(123)
model <- party::cforest(SpeciesID ~ .,
                 predictors, collapse = "+",
                 data = imputed_data %>% select(-TreeID),
                 control = ctree_control(mtry = 4, mincriterion = 0.95),
                 ntree = 25)

# Step 11: Compute variable importance (conditional permutation importance)
importance_scores <- varimp(model, conditional = TRUE)

# Step 12: Sort and display top 20 most important variables
importance_df <- as.data.frame(sort(importance_scores, decreasing = TRUE))
colnames(importance_df) <- "Importance"
importance_df$Variable <- rownames(importance_df)
rownames(importance_df) <- NULL

# View top 20
print(head(importance_df, 20))

# Step 13: Optional — Save to CSV
write_csv(importance_df, "E:/DATA/Variable_Importance_Scores.csv")
