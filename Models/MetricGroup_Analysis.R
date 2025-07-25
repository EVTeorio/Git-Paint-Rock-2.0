

library(dplyr)
library(ranger)
library(caret)
library(beepr)
library(readr)

# Load pre-processed canopy metrics data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# Remove specific canopies (no valid VI metrics)
Canopies_noVIs <- c("080736","090800","090830","090848","090863",
                    "100987","110563","110602","110618","110752","110819",
                    "110883","110913","120866","140806","150720","80917")

samples_clean <- combined_df[!combined_df$TreeID %in% Canopies_noVIs, ]

# Identify metric columns (exclude TreeID and SpeciesID)
all_metrics <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

# Remove columns with constant or all NA values
vi_filtered <- samples_clean %>%
  select(all_of(all_metrics)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# Recombine with IDs
canopy_means <- samples_clean %>%
  select(TreeID, SpeciesID) %>%
  bind_cols(vi_filtered)

# Count species
species_counts <- canopy_means %>%
  count(SpeciesID)

# Define metric groups
vi_vars <- c("Mean_NPCI", "Mean_SIPI", "Entropy_mSR",
             "Entropy_PAD_35_40_on", "SD_PAD_35_40_on", "Mean_Height_LeafOn",
             "Mean_PSRI", "SD_Carter2", "SD_PAD_30_35_on",
             "Max_Height_LeafOff", "Max_Carter6", "Mean_PSND",
             "Min_CRI3", "Entropy_PAD_30_35_on", "Mean_NDVI",
             "Max_Height_LeafOn", "Mean_GreenNDVI", "Mean_mNDVI",
             "Mean_PAD_35_40_on", "Max_DD", "Mean_OSAVI",
             "Max_PAD_30_35_off", "Mean_NDVI3"
)

leafon_vars <- c("Mean_NPCI", "Mean_PSRI" , "Mean_mNDVI", "Mean_SIPI",
                 "Mean_NDVI3", "Mean_Datt5", "Mean_MCARI2OSAVI2", "Mean_MCARI2",
                 "Mean_Height_LeafOn", "Max_Height_LeafOn", "Entropy_mSR", "Mean_PAD_35_40_on"
)

leafoff_vars <- c("Mean_PSRI", "Mean_Datt5", "Mean_MCARI2OSAVI2",
                  "Mean_Height_LeafOn", "Mean_PAD_35_40_on"
)

seasonal_var <- c(
  "Mean_Rumble_LeafOn", "Mean_Rumble_LeafOff", "Mean_Height_LeafOn", "Mean_Intensity_LeafOn",
  "Mean_Height_LeafOff", "Mean_PAD_40_45_on", "Mean_PAD_25_30_off", "Min_PAD_5_10_off",
  "Min_PAD_30_35_off", "Mean_PAD_35_40_on", "Entropy_PAD_40_45_off", "Max_PAD_0_5_off",
  "Mean_PAD_0_5_on", "Max_PAD_10_15_off", "Entropy_PAD_10_15_on", "Max_PAD_15_20_off",
  "SD_PAD_15_20_on", "SD_PAD_20_25_off", "SD_PAD_5_10_off", "Entropy_PAD_5_10_on",
  "Mean_PAD_15_20_on", "Mean_PAD_20_25_on", "Mean_NPCI", "Mean_CRI4",
  "Mean_MCARI", "Mean_PSRI", "Mean_DWSI4", "Max_DPI",
  "Min_MTCI", "Min_SRPI", "Mean_Vogelmann2", "Mean_CI",
  "Mean_Carter", "Min_Datt3", "Mean_PSSR", "Mean_REPLE",
  "Min_mSR705", "Max_EVI", "SD_OSAVI2", "SD_mSR",
  "Entropy_Datt6", "Entropy_GDVI3", "Entropy_SPVI", "Entropy_SR3",
  "Entropy_mNDVI"
)


model_inputs <- list(
  cforest23 = vi_vars,
  rforest5 = leafoff_vars,
  rforest10 = leafon_vars,
  rforest45 = seasonal_var
)

Final_grouped_results <- list()
importance_results <- list()

for (group_name in names(model_inputs)) {
  cat("\n============== Modeling Group:", group_name, "===============\n")
  
  metrics_in_use <- model_inputs[[group_name]]
  all_results <- list()
  group_importance <- list()
  
  for (i in 1:20) {
    cat("\n--- Iteration", i, "---\n")
    set.seed(50 + i)
    
    # Sampling logic (same as before)
    species_counts <- canopy_means %>% count(SpeciesID)
    species_12_plus <- species_counts %>% filter(n >= 25) %>% pull(SpeciesID)
    species_6_11 <- species_counts %>% filter(n >= 10 & n < 20) %>% pull(SpeciesID)
    rare_species <- species_counts %>% filter(n < 10) %>% pull(SpeciesID)
    
    sampled_12_plus <- canopy_means %>%
      filter(SpeciesID %in% species_12_plus) %>%
      group_by(SpeciesID) %>%
      slice_sample(n = 25) %>%
      ungroup()
    train_12_plus <- sampled_12_plus %>%
      group_by(SpeciesID) %>%
      slice_sample(n = 20) %>%
      ungroup()
    test_12_plus <- sampled_12_plus %>%
      filter(!TreeID %in% train_12_plus$TreeID)
    
    sampled_6_11 <- canopy_means %>%
      filter(SpeciesID %in% species_6_11) %>%
      group_by(SpeciesID) %>%
      slice_sample(n = 10) %>%
      ungroup()
    train_6_11 <- sampled_6_11 %>%
      group_by(SpeciesID) %>%
      slice_sample(n = 6) %>%
      ungroup()
    test_6_11 <- sampled_6_11 %>%
      filter(!TreeID %in% train_6_11$TreeID)
    
    rare_canopies <- canopy_means %>%
      filter(SpeciesID %in% rare_species)
    rare_train <- rare_canopies %>%
      group_by(SpeciesID) %>%
      slice_sample(n = 3) %>%
      ungroup() %>%
      mutate(SpeciesID = "others")
    rare_test <- rare_canopies %>%
      filter(!TreeID %in% rare_train$TreeID) %>%
      mutate(SpeciesID = "others")
    
    train_canopies <- bind_rows(train_12_plus, train_6_11, rare_train)
    test_canopies <- bind_rows(test_12_plus, test_6_11, rare_test)
    
    train_data <- train_canopies %>%
      select(SpeciesID, all_of(metrics_in_use)) %>%
      drop_na()
    test_data <- test_canopies %>%
      select(SpeciesID, all_of(metrics_in_use)) %>%
      drop_na()
    
    train_data$SpeciesID <- factor(train_data$SpeciesID)
    test_data$SpeciesID <- factor(test_data$SpeciesID, levels = levels(train_data$SpeciesID))
    
    # ==== MODEL TRAINING ====
    rf_mod <- ranger(
      formula = SpeciesID ~ .,
      data = train_data,
      importance = "impurity",
      num.trees = 3000,
      classification = TRUE,
      probability = TRUE
    )
    
    # ==== PREDICTION ====
    rf_preds <- predict(rf_mod, data = test_data)
    prob_matrix <- rf_preds$predictions
    pred_labels <- colnames(prob_matrix)[apply(prob_matrix, 1, which.max)]
    preds_factor <- factor(pred_labels, levels = levels(train_data$SpeciesID))
    
    # ==== EVALUATION ====
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
    
    # ==== CONFIDENCE AND RESULTS ====
    confidence <- apply(prob_matrix, 1, max)
    
    # results_df <- tibble(
    #   TreeID = test_canopies$TreeID,
    #   TrueLabel = test_data$SpeciesID,
    #   PredictedLabel = preds_factor,
    #   Confidence = confidence
    # )
    results_df <- tibble(
      TreeID = test_data$TreeID,
      TrueLabel = test_data$SpeciesID,
      PredictedLabel = preds_factor,
      Confidence = confidence
    )
    
    # ==== STORE RESULTS ====
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
    
    # Variable importance
    group_importance[[paste0("Sample_", i)]] <- rf_mod$variable.importance
  }
  
  Final_grouped_results[[group_name]] <- all_results
  importance_results[[group_name]] <- group_importance
}

beep()

# Save final results and importance metrics
saveRDS(Final_grouped_results, file = "E:/DATA/1Final_Model.rds")
saveRDS(importance_results, file = "E:/DATA/1Final_Importance.rds")
