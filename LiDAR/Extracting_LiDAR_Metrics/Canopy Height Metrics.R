

# Load required libraries
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
chm_on_path <- "E:/Updated LiDAR/PRFPD_CHM_leafOn.tiff"
chm_off_path <- "E:/Updated LiDAR/PRFPD_CHM_leafOff.tiff"

# Load rasters
chm_on <- rast(chm_on_path)
chm_off <- rast(chm_off_path)


# List shapefiles
shapefiles <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)

# Function to extract pixel values and return a combined dataframe
extract_pixels <- function(chm_raster, shapefiles, chm_label) {
  chm_crs <- terra::crs(chm_raster)
  pixels_list <- list()
  
  for (shp_path in shapefiles) {
    cat("Processing:", basename(shp_path), "-", chm_label, "\n")
    
    canopy_sf <- st_read(shp_path, quiet = TRUE)
    canopy_sf <- st_transform(canopy_sf, crs = chm_crs)
    
    extracted_pixels <- exact_extract(
      chm_raster, 
      canopy_sf,
      include_cols = "Canopies",
      progress = FALSE
    )
    
    polygon_pixels_df <- map_df(extracted_pixels, ~ {
      df <- filter(.x, !is.na(value))
      df %>% select(Canopies, value)
    })
    
    pixels_list[[basename(shp_path)]] <- polygon_pixels_df
  }
  
  all_pixels_df <- bind_rows(pixels_list, .id = "Shapefile_Source") %>%
    separate(Canopies, into = c("TreeID", "SpeciesID"), sep = "_", remove = TRUE) %>%
    rename(Height = value)
  
  all_pixels_df$CHM_Type <- chm_label
  return(all_pixels_df)
}

# Extract pixel data for both CHMs
pixels_on <- extract_pixels(chm_on, shapefiles, "LeafOn")
pixels_off <- extract_pixels(chm_off, shapefiles, "LeafOff")
beep()

# Combine both datasets
all_pixels_combined <- bind_rows(pixels_on, pixels_off)

# Function to compute entropy
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# Compute stats for both CHMs separately
canopy_stats_on <- pixels_on %>%
  group_by(TreeID) %>%
  summarise(
    Max_Height_LeafOn = max(Height, na.rm = TRUE),
    Min_Height_LeafOn = min(Height, na.rm = TRUE),
    Mean_Height_LeafOn = mean(Height, na.rm = TRUE),
    SD_Height_LeafOn = sd(Height, na.rm = TRUE),
    Entropy_Height_LeafOn = compute_entropy(Height)
  )

canopy_stats_off <- pixels_off %>%
  group_by(TreeID) %>%
  summarise(
    Max_Height_LeafOff = max(Height, na.rm = TRUE),
    Min_Height_LeafOff = min(Height, na.rm = TRUE),
    Mean_Height_LeafOff = mean(Height, na.rm = TRUE),
    SD_Height_LeafOff = sd(Height, na.rm = TRUE),
    Entropy_Height_LeafOff = compute_entropy(Height)
  )

# Merge both CHM stats by TreeID
canopy_stats <- full_join(canopy_stats_on, canopy_stats_off, by = "TreeID")

# Optional: write output
write.csv(canopy_stats, "E:/DATA/CHM_Metrics.csv", row.names = FALSE)

beep()
