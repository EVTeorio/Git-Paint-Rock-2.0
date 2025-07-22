

# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(entropy)
library(purrr)
library(beepr)

# Load PAD raster
PAD_raster <- rast("E:/LiDAR_Metrics/ALL_Metrics.tif")
pad_names <- names(PAD_raster)

# Load canopy polygons
canopies_path <- "E:/Git Paint Rock 1.0/Updated Canopy Polygons/Updated/"
shapefiles <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)

crs_raster_path <- "E:/Updated LiDAR/PRFPD_CHM_leafOff.tiff"
justforcrs <- rast(crs_raster_path)

# Apply CRS to LiDAR
crs(PAD_raster) <- crs(justforcrs)

# Function to extract pixel values for each PAD layer
extract_PAD_pixels <- function(pad_raster, shapefiles) {
  pixels_list <- list()
  
  for (shp_path in shapefiles) {
    cat("Processing:", basename(shp_path), "\n")
    canopy_sf <- st_read(shp_path, quiet = TRUE)
    canopy_sf <- st_transform(canopy_sf, crs = st_crs(pad_raster))
    
    extracted <- exact_extract(
      pad_raster,
      canopy_sf,
      include_cols = "Canopies",
      progress = FALSE
    )
    
    # Flatten each record: long format with all layer values
    polygon_pixels_df <- map_df(extracted, function(df) {
      df <- df %>% filter(!is.na(df[,1]))  # filter out fully NA rows
      df$Canopies <- unique(df$Canopies)
      df
    })
    
    pixels_list[[basename(shp_path)]] <- polygon_pixels_df
  }
  
  all_pixels <- bind_rows(pixels_list, .id = "Shapefile_Source") %>%
    separate(Canopies, into = c("TreeID", "SpeciesID"), sep = "_", remove = TRUE)
  
  return(all_pixels)
}

# Step 1: Extract all pixel values
pad_pixels_df <- extract_PAD_pixels(PAD_raster, shapefiles)

# Step 2: Calculate stats for each TreeID and layer
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# Step 3: Convert to long format for grouped summary
pad_long <- pad_pixels_df %>%
  pivot_longer(
    cols = all_of(pad_names),
    names_to = "Layer",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Step 4: Summarize for each TreeID and Layer
pad_summary <- pad_long %>%
  group_by(TreeID, Layer) %>%
  summarise(
    Max = max(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Entropy = compute_entropy(Value),
    .groups = "drop"
  )

# Step 5: Pivot to wide format
pad_wide <- pad_summary %>%
  pivot_wider(
    names_from = Layer,
    values_from = c(Max, Min, Mean, SD, Entropy),
    names_glue = "{.value}_{Layer}"
  )

# Step 6: Export
write.csv(pad_wide, "E:/DATA/PAD_Metrics.csv", row.names = FALSE)
beep()