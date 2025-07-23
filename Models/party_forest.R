
library(VIM)      # For KNN imputation
library(tibble) 
library(party)

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

df <- imputed_data

# Convert character columns to factors
df[] <- lapply(df, function(x) if (is.character(x)) as.factor(x) else x)

# Define response and predictors
# Set the target variable (e.g., "Ozone") and ensure it's not all NA
df <- df[!is.na(df$SpeciesID), ]  # Remove rows where target is NA
response_var <- "SpeciesID"

# Convert to formula
predictors <- setdiff(names(df), response_var)
formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))

# Train conditional inference forest
set.seed(123)  # for reproducibility
cforest_model <- party::cforest(
  formula,
  data = df,
  controls = cforest_unbiased(mtry = 2, ntree = 1000)
)
beep()

# Calculate variable importance (conditional = TRUE for conditional importance)
varimp_result <- party::varimp(cforest_model, conditional = TRUE)
beep(3)
# Sort and display importance
varimp_sorted <- sort(varimp_result, decreasing = TRUE)
print(varimp_sorted)

varimp_frame <- as.data.frame(varimp_sorted)


write.csv(varimp_frame, "E:/DATA/Metric_Importance_Selection/Metric_Importance_Values.csv")
