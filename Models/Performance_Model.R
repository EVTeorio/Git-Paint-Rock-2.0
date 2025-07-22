

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

# Final metrics
metrics <- c("Mean_NPCI", "Mean_Height_LeafOn")

# Recombine with ID columns
canopy_means <- samples_clean %>%
  select(TreeID, SpeciesID) %>%
  bind_cols(vi_filtered)

# Species to keep
species_counts <- canopy_means %>% count(SpeciesID)
All_species <- species_counts %>% filter(n > 0) %>% pull(SpeciesID)

# Output lists
all_results <- list()
importance_results <- list()

for (i in 1:50) {
  cat("\n--- Object-Based Sample", i, "---\n")
  set.seed(50 + i)
  
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
  
  train_data <- train_canopies %>% select(SpeciesID, all_of(metrics))
  test_data  <- test_canopies %>% select(SpeciesID, all_of(metrics))
  
  train_data$SpeciesID <- factor(train_data$SpeciesID)
  test_data$SpeciesID <- factor(test_data$SpeciesID, levels = levels(train_data$SpeciesID))
  
  # Train model
  rf_mod <- ranger(
    formula = SpeciesID ~ .,
    data = train_data,
    importance = "impurity",
    num.trees = 3000,
    classification = TRUE,
    probability = TRUE
  )
  
  # Prediction (get probabilities)
  rf_preds <- predict(rf_mod, data = test_data)
  prob_matrix <- rf_preds$predictions
  pred_labels <- colnames(prob_matrix)[apply(prob_matrix, 1, which.max)]
  preds_factor <- factor(pred_labels, levels = levels(train_data$SpeciesID))
  
  # Evaluation
  cm <- confusionMatrix(preds_factor, test_data$SpeciesID)
  
  precision <- cm$byClass[, "Precision"]
  recall <- cm$byClass[, "Recall"]
  f1_scores <- 2 * (precision * recall) / (precision + recall)
  
  sensitivity <- cm$byClass[, "Sensitivity"]
  specificity <- cm$byClass[, "Specificity"]
  balanced_accuracy_scores <- (sensitivity + specificity) / 2
  
  macro_f1 <- mean(f1_scores, na.rm = TRUE)
  macro_balanced_accuracy <- mean(balanced_accuracy_scores, na.rm = TRUE)
  accuracy <- as.numeric(cm$overall["Accuracy"])
  kappa_value <- cm$overall["Kappa"]
  
  # Confidence scores
  confidence <- apply(prob_matrix, 1, max)
  
  # Results DF
  results_df <- tibble(
    TreeID = test_canopies$TreeID,
    TrueLabel = test_data$SpeciesID,
    PredictedLabel = preds_factor,
    Confidence = confidence
  )
  
  # Store all
  all_results[[paste0("Sample_", i)]] <- list(
    model = rf_mod,
    accuracy = accuracy,
    kappa = as.numeric(kappa_value),
    macro_f1 = macro_f1,
    macro_balanced_accuracy = macro_balanced_accuracy,
    f1_per_class = f1_scores,
    balanced_accuracy_per_class = balanced_accuracy_scores,
    results = results_df,
    confidence = confidence,
    confusion_matrix = cm
  )
  
  importance_results[[paste0("Sample_", i)]] <- rf_mod$variable.importance
}

# Done sound
beep()

# Save
saveRDS(all_results, "E:/DATA/Perfomance/Performance_Model.rds")


