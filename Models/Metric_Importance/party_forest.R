
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

all_metrics <- c(
  "Mean_NPCI", "Mean_SIPI", "Mean_mNDVI",
  "Mean_PSRI", "Mean_PSND", "Mean_NDVI3",
  "Min_CRI3", "Entropy_mSR", "Mean_MPRI",
  "SD_Carter3", "Entropy_PAD_30_35_on", "Max_Height_LeafOff",
  "Mean_DD", "Max_Carter6", "Max_TVI",
  "SD_TGI", "Entropy_PAD_35_40_on", "Max_EVI",
  "Mean_GreenNDVI", "Entropy_Vogelmann3", "Mean_NDVI",
  "SD_mNDVI", "Mean_SR5", "SD_TVI",
  "Max_TCARI2"
)

# Step 4: Identify metric columns (excluding ID columns)
#all_metrics <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

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
  controls = cforest_unbiased(mtry = 3, ntree = 2000)
)
beep()

# Calculate variable importance (conditional = TRUE for conditional importance)
varimp_result <- party::varimp(cforest_model, conditional = TRUE)
beep()

# Sort and display importance
varimp_sorted <- sort(varimp_result, decreasing = TRUE)
print(varimp_sorted)

varimp_frame <- as.data.frame(varimp_sorted)

# Convert to data frame and sort by importance
importance_df_plot <- data.frame(
  Variable = rownames(varimp_frame),
  Importance = varimp_sorted
) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 50)  # keep only top 50

# Plot top 50
ggplot(importance_df_plot, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top Variables by Average Importance",
    x = "Metrics",
    y = "Mean Importance Score"
  ) +
  theme_minimal(base_size = 13)




print(rownames(varimp_frame))

write.csv(varimp_frame, "E:/DATA/Metric_Importance_Selection/cforest50Metric_Importance_Values.csv")
varimp_df <- read.csv("E:/DATA/Metric_Importance_Selection/Metric_Importance_Values.csv")

