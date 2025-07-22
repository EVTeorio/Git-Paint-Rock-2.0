

library(tidyverse)
library(ranger)
library(caret)
library(beepr)

# Load data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

Canopies_noVIs <- c("080736","090800","090830","090848","090863",
                    "100987","110563","110602","110618","110752","110819",
                    "110883","110913","120866","140806","150720","80917")

# Identify VI columns
all_metrics <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

samples_clean <- combined_df[!combined_df$TreeID %in% Canopies_noVIs, ]

# Remove constant or all-NA columns
vi_filtered <- samples_clean %>%
  select(all_of(all_metrics)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# Identify VI columns (exclude ID columns)
metrics <- names(vi_filtered)[!(names(vi_filtered) %in% c("TreeID", "SpeciesID"))]

# Recombine with ID columns
canopy_means <- samples_clean %>%
  select(TreeID, SpeciesID) %>%
  bind_cols(vi_filtered)

# Count species
species_counts <- canopy_means %>% count(SpeciesID)
All_species <- species_counts %>% filter(n > 0) %>% pull(SpeciesID)
  
all_results <- list()
importance_results <- list()

for (i in 1:50) {
  cat("\n--- Object-Based Sample", i, "---\n")
  set.seed(50 + i)
  
  # Split common species into train/test
  train_list <- list()
  test_list <- list()
  
  for (sp in All_species) {
    sp_data <- canopy_means %>% filter(SpeciesID == sp)
    sp_test <- sp_data %>% slice_sample(prop = 0.25)
    sp_train <- anti_join(sp_data, sp_test, by = "TreeID")
    train_list[[sp]] <- sp_train
    test_list[[sp]] <- sp_test
  }
  
  train_canopies <- bind_rows(train_list)
  test_canopies <- bind_rows(test_list)
  
  # Ensure levels match and no NAs
  train_data <- train_canopies %>%
    select(SpeciesID, all_of(metrics)) #%>%
    #drop_na()
  
  test_data <- test_canopies %>%
    select(SpeciesID, all_of(metrics)) #%>%
    #drop_na()
  
  train_data$SpeciesID <- factor(train_data$SpeciesID)
  test_data$SpeciesID <- factor(test_data$SpeciesID, levels = levels(train_data$SpeciesID))
  
  # Train using ranger
  rf_model <- ranger(
    formula = SpeciesID ~ .,
    data = train_data,
    importance = "impurity",
    num.trees = 3000,
    classification = TRUE,
    probability = FALSE
  )
  
  # Predict
  preds <- predict(rf_model, data = test_data)$predictions
  preds_factor <- factor(preds, levels = levels(train_data$SpeciesID))
  
  # Confusion matrix
  cm <- confusionMatrix(preds_factor, test_data$SpeciesID)
  
  precision <- cm$byClass[, "Precision"]
  recall <- cm$byClass[, "Recall"]
  f1_by_class <- 2 * (precision * recall) / (precision + recall)
  f1_macro <- mean(f1_by_class, na.rm = TRUE)
  
  # Store results
  all_results[[paste0("Sample_", i)]] <- list(
    model = rf_model,
    confusion = cm,
    accuracy = cm$overall["Accuracy"],
    f1_by_class = f1_by_class,
    f1_macro = f1_macro
  )
  
  importance_results[[paste0("Sample_", i)]] <- rf_model$variable.importance
}

beep()

help("importance")

saveRDS(all_results, "E:/DATA/Metric_Importance_Selection/Metric_Importance_Model.rds")
saveRDS(importance_results, "E:/DATA/Metric_Importance_Selection/Metric_Importance_Values.rds")

