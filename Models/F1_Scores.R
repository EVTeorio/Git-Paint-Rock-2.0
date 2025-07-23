




library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

F1_results <- readRDS("E:/DATA/Perfomance/Performance_Model_ALL_Metric.rds")

# Initialize a list to hold each iteration's summary
summary_list <- list()

# Loop through each iteration in F1_results
for (iter_name in names(F1_results)) {
  res <- F1_results[[iter_name]]
  
  # Extract overall metrics
  accuracy <- res$accuracy %||% NA
  macro_f1 <- res$macro_f1 %||% NA
  kappa <- res$kappa %||% NA
  
  # Extract and clean per-class F1 scores
  class_f1 <- if (!is.null(res$f1_per_class) && length(res$f1_per_class) > 0) {
    vals <- as.numeric(res$f1_per_class)
    names(vals) <- gsub("^Class: ", "F1_", names(res$f1_per_class))
    vals
  } else {
    numeric(0)
  }
  
  # Assemble a data frame for this iteration
  df_row <- tibble(
    Iteration = iter_name,
    Accuracy = accuracy,
    F1_Macro = macro_f1,
    Kappa = kappa
  ) %>% bind_cols(as_tibble(as.list(class_f1)))
  
  summary_list[[iter_name]] <- df_row
}

# Combine all into one summary data frame
summary_df <- bind_rows(summary_list)

# View the compiled results
print(summary_df)

# Compute summary stats across iterations
model_stats <- summary_df %>%
  summarise(
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    SD_Accuracy   = sd(Accuracy, na.rm = TRUE),
    Mean_F1_Macro = mean(F1_Macro, na.rm = TRUE),
    SD_F1_Macro   = sd(F1_Macro, na.rm = TRUE),
    Mean_Kappa    = mean(Kappa, na.rm = TRUE),
    SD_Kappa      = sd(Kappa, na.rm = TRUE)
  )

print(model_stats)

# Optional: Plot Accuracy and Macro F1 over iterations
summary_df %>%
  select(Iteration, Accuracy, F1_Macro) %>%
  pivot_longer(-Iteration, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Iteration, y = Value, color = Metric, group = Metric)) +
  geom_line() +
  geom_point() +
  labs(title = "Model Performance 609 Metrics",
       x = "Iteration", y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(summary_df, aes(y = F1_Macro)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Distribution of Macro F1 Scores 609 Metrics",
       y = "Macro F1 Score") +
  theme_minimal()


# Pivot to long format for per-class F1
f1_long <- summary_df %>%
  pivot_longer(cols = starts_with("F1_"), 
               names_to = "Class", 
               values_to = "F1_Score") %>%
  filter(Class != "F1_Macro")  # exclude macro if needed

# Plot boxplot by class
ggplot(f1_long, aes(x = Class, y = F1_Score)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Per-Class F1 Score 609 metrics",
       x = "Class", y = "F1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




##############################################################################

# Ensure numeric columns are treated as such
numeric_cols <- c("Accuracy", "F1_Macro", "F1_ACSA3C", "F1_CACA38", "F1_CAOV2",
                  "F1_FAGR", "F1_FRAMCO", "F1_LIST2", "F1_LITU", "F1_PIEC2",
                  "F1_QUAL", "F1_QURU", "F1_QUSH", "F1_TIAM", "F1_others")
summary_df[numeric_cols] <- lapply(summary_df[numeric_cols], as.numeric)

# Summarize mean and standard deviation by Group
summary_stats <- summary_df %>%
  group_by(Group) %>%
  summarise(across(all_of(numeric_cols),
                   list(mean = ~mean(. , na.rm = TRUE),
                        sd = ~sd(. , na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

# View result
print(summary_stats, n = Inf, width = Inf)


##############################################################################
# box plot per species/class F1
# Pivot longer for F1 per class, exclude Macro
f1_long <- summary_df %>%
  pivot_longer(cols = starts_with("F1_") & !starts_with("F1_Macro"),
               names_to = "Class", values_to = "F1") %>%
  filter(!is.na(F1))

ggplot(f1_long, aes(x = Class, y = F1, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "F1 Score by Class and Group", x = "Class", y = "F1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
########################################################################

# Pivot longer for F1 per class, exclude Macro, and remove specific species
f1_long <- summary_df %>%
  pivot_longer(
    cols = starts_with("F1_") & !starts_with("F1_Macro"),
    names_to = "Class", values_to = "F1"
  ) %>%
  filter(!is.na(F1)) %>%
  # Remove "F1_" prefix and filter out unwanted species
  mutate(Class = str_remove(Class, "^F1_")) %>%
  filter(!Class %in% c("ACSA3C", "FAGR", "QURU", "QUSH", "TIAM"))

# Plot
ggplot(f1_long, aes(x = Class, y = F1, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "F1 Score by Class and Group", x = "Class", y = "F1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

