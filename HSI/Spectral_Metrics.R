

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(exactextractr)
library(entropy)
library(beepr)

# ---- File paths ----
hsi_folder <- "E:/HSI_Files_Parsing/"
shapefile_folder <- "E:/Git Paint Rock 1.0/Updated Canopy Polygons/Updated/"
output_csv <- "E:/DATA/HSI_Metrics.csv"

# ---- Step 1: Load and filter filenames ----
shapefiles <- list.files(shapefile_folder, pattern = "\\.shp$", full.names = TRUE)
shapefile_numbers <- gsub("\\D", "", basename(shapefiles))  # Extract numbers from shapefile names

allfiles <- list.files(hsi_folder)
hsi_files_all <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]  # Filter only .tif
img_numbers <- gsub("\\D", "", hsi_files_all)  # Extract numbers from HSI filenames

# ---- Step 2: Match shapefiles and HSI files by image number ----
matched_files <- intersect(img_numbers, shapefile_numbers)  # Only process where both exist

# Function to compute entropy
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# ---- Step 4: Loop and extract ----
extracted_list <- lapply(matched_files, function(img_number) {
  cat("Processing image number:", img_number, "\n")
  
  # Match paths
  img_idx <- which(img_numbers == img_number)
  shp_idx <- which(shapefile_numbers == img_number)
  
  if (length(img_idx) == 0 || length(shp_idx) == 0) return(NULL)  # Skip if missing
  
  img_path <- file.path(hsi_folder, hsi_files_all[img_idx])
  shp_path <- shapefiles[shp_idx]
  
  # Load data
  hsi_raster <- rast(img_path)
  band_names <- names(hsi_raster)
  canopy_sf <- st_read(shp_path, quiet = TRUE)
  
  # Extract values
  extracted <- exact_extract(
    hsi_raster,
    canopy_sf,
    include_cols = "Canopies",
    progress = FALSE
  )
  
  polygon_pixels_df <- map_df(extracted, function(df) {
    df <- df %>% filter(!is.na(df[, 1]))
    df$Canopies <- unique(df$Canopies)
    df
  })
  
  polygon_pixels_df$ImageID <- img_number
  return(polygon_pixels_df)
})
beep()


# ---- Step 5: Combine and reshape ----
hsi_pixels_df <- bind_rows(extracted_list) %>%
  separate(Canopies, into = c("TreeID", "SpeciesID"), sep = "_", remove = TRUE)
print(hsi_wide)
# Long format
hsi_long <- hsi_pixels_df %>%
  pivot_longer(
    cols = matches("^\\d+\\.\\d+\\s*nm$"),
    names_to = "Wavelength",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Summarize per tree
hsi_summary <- hsi_long %>%
  group_by(TreeID, Wavelength) %>%
  summarise(
    Max = max(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Entropy = compute_entropy(Value),
    .groups = "drop"
  )
beep()
# Wide format
hsi_wide <- hsi_summary %>%
  pivot_wider(
    names_from = Wavelength,
    values_from = c(Max, Min, Mean, SD, Entropy),
    names_glue = "{.value}_{Wavelength}"
  )

# Export
write.csv(hsi_wide, output_csv, row.names = FALSE)
beep()
