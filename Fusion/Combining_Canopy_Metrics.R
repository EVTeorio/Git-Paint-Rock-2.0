

library(dplyr)
library(purrr)
library(readr)
library(stringr)

# ---- Step 1: Set folder path ----
csv_folder <- "E:/DATA/Canopy_values_dfs/"  # Change this to your folder
output_csv <- "E:/DATA/All_CanopyMetrics.csv"

# ---- Step 2: List all CSV files in the folder ----
csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

# ---- Step 3: Read and merge all CSVs by TreeID ----
# This assumes each CSV has a column "TreeID" to join on

combined_df <- csv_files %>%
  map(read_csv) %>%
  reduce(full_join, by = "TreeID")

# ---- Step 4: Save combined CSV (optional) ----
write_csv(combined_df, output_csv)

# ---- View the result ----
print(head(combined_df))
