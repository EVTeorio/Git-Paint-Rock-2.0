

# Load necessary libraries
library(lidR)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(exactextractr)
library(entropy)
library(beepr)


# Define file paths
canopies_path <- "E:/Git Paint Rock 1.0/Updated Canopy Polygons/Updated/"
LAZon_path <- "E:/Updated LiDAR/PaintRock_20ha_leafOn_class.laz"
LAZoff_path <- "E:/Updated LiDAR/PaintRock_20ha_leafOff_class.laz"

# Function to generate an intensity raster from a LAZ file
create_intensity_raster <- function(laz_path, res = 1) {
  las <- readLAS(laz_path)
  if (is.empty(las)) stop("LAS file is empty or invalid.")
  grid_metrics(las, ~mean(Intensity, na.rm = TRUE), res = res)
}

# Create rasters for both LAZ files
intensity_raster_on <- create_intensity_raster(LAZon_path)
intensity_raster_off <- create_intensity_raster(LAZoff_path)
beep()
plot(intensity_raster_on)

# List shapefiles
shapefiles <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)

# Function to extract intensity pixels and return data
extract_intensity_pixels <- function(intensity_raster, shapefiles, label) {
  chm_crs <- st_crs(intensity_raster)
  pixels_list <- list()
  
  for (shp_path in shapefiles) {
    cat("Processing:", basename(shp_path), "-", label, "\n")
    canopy_sf <- st_read(shp_path, quiet = TRUE)
    canopy_sf <- st_transform(canopy_sf, crs = chm_crs)
    
    extracted <- exact_extract(
      intensity_raster,
      canopy_sf,
      include_cols = "Canopies",
      progress = FALSE
    )
    
    df <- map_df(extracted, ~ filter(.x, !is.na(value)) %>% select(Canopies, value))
    pixels_list[[basename(shp_path)]] <- df
  }
  
  all_pixels <- bind_rows(pixels_list, .id = "Shapefile_Source") %>%
    separate(Canopies, into = c("TreeID", "SpeciesID"), sep = "_", remove = TRUE) %>%
    rename(Intensity = value)
  
  all_pixels$Source <- label
  return(all_pixels)
}

# Extract pixel values for intensity rasters
pixels_int_on <- extract_intensity_pixels(intensity_raster_on, shapefiles, "LeafOn")
pixels_int_off <- extract_intensity_pixels(intensity_raster_off, shapefiles, "LeafOff")

# Combine both into one dataframe
all_intensity_pixels <- bind_rows(pixels_int_on, pixels_int_off)

# Entropy calculation function
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# Summarise for Leaf-On
intensity_stats_on <- pixels_int_on %>%
  group_by(TreeID) %>%
  summarise(
    Max_Intensity_LeafOn = max(Intensity, na.rm = TRUE),
    Min_Intensity_LeafOn = min(Intensity, na.rm = TRUE),
    Mean_Intensity_LeafOn = mean(Intensity, na.rm = TRUE),
    SD_Intensity_LeafOn = sd(Intensity, na.rm = TRUE),
    Entropy_Intensity_LeafOn = compute_entropy(Intensity)
  )

# Summarise for Leaf-Off
intensity_stats_off <- pixels_int_off %>%
  group_by(TreeID) %>%
  summarise(
    Max_Intensity_LeafOff = max(Intensity, na.rm = TRUE),
    Min_Intensity_LeafOff = min(Intensity, na.rm = TRUE),
    Mean_Intensity_LeafOff = mean(Intensity, na.rm = TRUE),
    SD_Intensity_LeafOff = sd(Intensity, na.rm = TRUE),
    Entropy_Intensity_LeafOff = compute_entropy(Intensity)
  )

# Join both summaries
canopy_intensity_stats <- full_join(intensity_stats_on, intensity_stats_off, by = "TreeID")

# Write result
write.csv(canopy_intensity_stats, "E:/DATA/Intensity_Metrics.csv", row.names = FALSE)

beep()
