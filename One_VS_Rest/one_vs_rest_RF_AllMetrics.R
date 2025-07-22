


library(tidyverse)
library(ranger)
library(caret)
library(beepr)

# Load data
spec_chem_canopy <- read.csv("E:/Thesis_Final_Data/ALLmetrics_clean_sunlit_no_nas.csv")
beep()
str(spec_chem_canopy)


# Get canopy and species info
canopies <- spec_chem_canopy %>%
  group_by(TreeID) %>%
  slice(1) %>%
  ungroup() %>%
  select(TreeID, SpeciesID)

# Define vegetation indices (VIs)
vi_vars <- c(
  "Boochs", "Boochs2", "CARI", "Carter", "Carter2", "Carter3", "Carter4", "Carter5", "Carter6",
  "CI", "CI2", "ClAInt", "CRI1", "CRI2", "CRI3", "CRI4", "D1", "D2", "Datt", "Datt2", "Datt3",
  "Datt4", "Datt5", "Datt6", "DD", "DDn", "DPI", "DWSI4", "EGFN", "EGFR", "EVI", "GDVI2",
  "GDVI3", "GDVI4", "GI", "Gitelson", "Gitelson2", "GMI1", "GMI2", "GreenNDVI", "Maccioni",
  "MCARI", "MCARIOSAVI", "MCARI2", "MCARI2OSAVI2", "mND705", "mNDVI", "MPRI", "MSAVI", "mSR",
  "mSR2", "mSR705", "MTCI", "MTVI", "NDVI", "NDVI2", "NDVI3", "NPCI", "OSAVI", "OSAVI2",
  "PARS", "PRI", "PRICI2", "PRInorm", "PSND", "PSRI", "PSSR", "RDVI", "REPLE", "REPLi",
  "SAVI", "SIPI", "SPVI", "SR", "SR1", "SR2", "SR3", "SR4", "SR5", "SR6", "SR7", "SR8", "SRPI",
  "SumDr1", "SumDr2", "TCARI", "TCARIOSAVI", "TCARI2", "TCARI2OSAVI2", "TGI", "TVI",
  "Vogelmann", "Vogelmann2", "Vogelmann3", "Vogelmann4"
)

# Identify LiDAR variables
leafon_vars <- c("PAD_0_5_on", "PAD_10_15_on", "PAD_15_20_on", "PAD_20_25_on", "PAD_25_30_on",
                 "PAD_30_35_on", "PAD_35_40_on","PAD_5_10_on")

leafoff_vars <- c("PAD_0_5_off", "PAD_10_15_off", "PAD_15_20_off", "PAD_20_25_off", "PAD_25_30_off",
                  "PAD_30_35_off", "PAD_35_40_off", "PAD_5_10_off")

seasonal_var <- "Seasonal_Occupancy_20_35m"

model_inputs <- list(
  VIs_only = vi_vars,
  VIs_LiDARleafon = c(vi_vars, leafon_vars),
  VIs_LiDARleafoff = c(vi_vars, leafoff_vars),
  VIs_allLiDAR = c(vi_vars, leafon_vars, leafoff_vars, seasonal_var)
)

metrics <- c(vi_vars, leafon_vars, leafoff_vars, seasonal_var)

# Average metrics by canopy
canopy_means <- spec_chem_canopy %>%
  group_by(TreeID, SpeciesID) %>%
  summarise(across(all_of(metrics), mean, na.rm = TRUE), .groups = "drop")


# Containers for results
all_results <- list()

for (i in 1:50) {
  cat("\n--- Sample", i, "---\n")
  set.seed(50 + i)
  
  # Species grouping
  species_counts <- canopy_means %>% count(SpeciesID)
  rare_species <- species_counts %>% filter(n < 6) %>% pull(SpeciesID)
  common_species <- species_counts %>% filter(n > 5) %>% pull(SpeciesID)
  
  # Split common species into train/test
  train_list <- list()
  test_list <- list()
  
  for (sp in common_species) {
    sp_data <- canopy_means %>% filter(SpeciesID == sp)
    sp_train <- sp_data %>% slice_sample(prop = 0.6)
    sp_test <- anti_join(sp_data, sp_train, by = "TreeID")
    train_list[[sp]] <- sp_train
    test_list[[sp]] <- sp_test
  }
  
  train_common <- bind_rows(train_list)
  test_common <- bind_rows(test_list)
  
  # Handle rare species: 1 sample for train, rest for test, recode as "others"
  rare_data <- canopy_means %>% filter(SpeciesID %in% rare_species)
  rare_train <- rare_data %>% group_by(SpeciesID) %>% slice_sample(n = 1) %>% ungroup() %>%
    mutate(SpeciesID = "others")
  rare_test <- anti_join(rare_data, rare_train, by = "TreeID") %>%
    mutate(SpeciesID = "others")
  
  # Final train/test sets
  train_canopies <- bind_rows(train_common, rare_train)
  test_canopies <- bind_rows(test_common, rare_test)
  
  # Select metrics and drop NA
  train_data <- train_canopies %>%
    select(SpeciesID, all_of(metrics)) #%>%
  #drop_na()
  test_data <- test_canopies %>%
    select(SpeciesID, all_of(metrics)) #%>%
  #drop_na()
  
  # Ensure factor consistency
  all_classes <- union(unique(train_data$SpeciesID), unique(test_data$SpeciesID))
  train_data$SpeciesID <- factor(train_data$SpeciesID, levels = all_classes)
  test_data$SpeciesID <- factor(test_data$SpeciesID, levels = all_classes)
  
  # One-vs-Rest loop
  one_vs_rest_results <- list()
  
  for (cls in levels(train_data$SpeciesID)) {
    cat("→ Training one-vs-rest model for:", cls, "\n")
    
    train_binary <- train_data %>%
      mutate(target = if_else(SpeciesID == cls, 1L, 0L))
    
    test_binary <- test_data %>%
      mutate(target = if_else(SpeciesID == cls, 1L, 0L))
    
    # Skip if only one class is present in training or testing
    if (length(unique(train_binary$target)) < 2 || length(unique(test_binary$target)) < 2) {
      cat("Skipping", cls, "- insufficient class diversity\n")
      next
    }
    
    # Train model
    rf_model <- ranger(
      formula = target ~ . -SpeciesID,
      data = train_binary,
      num.trees = 3000,
      #probability = TRUE,
      importance = "impurity"
    )
    
    # Predict safely
    pred_out <- predict(rf_model, data = test_binary)$predictions
    
    if (is.matrix(pred_out) && "1" %in% colnames(pred_out)) {
      preds_prob <- pred_out[, "1"]
    } else if (is.vector(pred_out)) {
      preds_prob <- pred_out
    } else if (is.matrix(pred_out) && ncol(pred_out) == 1) {
      preds_prob <- pred_out[, 1]
    } else {
      stop("Unknown format in ranger prediction output.")
    }
    
    preds_class <- factor(if_else(pred_out >= 0.5, 1L, 0L), levels = c(0, 1))
    actual_class <- factor(test_binary$target, levels = c(0, 1))
    
    if (length(preds_class) != length(actual_class)) {
      stop("Prediction and target lengths do not match.")
    }
    
    cm <- confusionMatrix(preds_class, actual_class)
    
    precision <- cm$byClass["Pos Pred Value"]
    recall <- cm$byClass["Sensitivity"]
    f1 <- if (!is.na(precision) && !is.na(recall)) {
      2 * (precision * recall) / (precision + recall)
    } else {
      NA
    }
    
    one_vs_rest_results[[cls]] <- list(
      model = rf_model,
      confusion = cm,
      accuracy = cm$overall["Accuracy"],
      precision = precision,
      recall = recall,
      f1 = f1,
      importance = rf_model$variable.importance
    )
  }
  
  all_results[[paste0("Sample_", i)]] <- one_vs_rest_results
}
beep()

# Save results
saveRDS(all_results, "E:/Results/One_vs_Rest_Metric Importance/OneVsRest_Full_RF_results_all_metrics.rds")
