

# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(entropy)
library(purrr)
library(beepr)


# ---- USER PATHS ----
canopies_path <- "E:/Git Paint Rock 1.0/Updated Canopy Polygons/Updated/"
chm_on_path <- "E:/Updated LiDAR/PRFPD_CHM_leafOn.tiff"
chm_off_path <- "E:/Updated LiDAR/PRFPD_CHM_leafOff.tiff"

# ---- LOAD RASTERS ----
chm_on <- rast(chm_on_path)
chm_off <- rast(chm_off_path)

# ---- RUMBLE INDEX FUNCTION ----
compute_rumble <- function(chm, window_size = 2) {
  # Apply focal standard deviation
  focal(chm, w = matrix(1, nrow = window_size, ncol = window_size), fun = sd, na.policy = "omit", na.rm = TRUE)
}

# ---- CREATE RUMBLE RASTERS ----
rumble_on <- compute_rumble(chm_on, window_size = 3)
rumble_off <- compute_rumble(chm_off, window_size = 3)
plot(rumble_off)

names(rumble_on) <- "Rumble_LeafOn"
names(rumble_off) <- "Rumble_LeafOff"

# ---- LOAD CANOPY SHAPEFILES ----
shapefiles <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)

# ---- EXTRACT RUMBLE PIXELS FROM EACH CANOPY ----
extract_rumble_pixels <- function(rumble_raster, shapefiles, label) {
  chm_crs <- st_crs(rumble_raster)
  pixels_list <- list()
  
  for (shp_path in shapefiles) {
    cat("Processing:", basename(shp_path), "-", label, "\n")
    
    canopy_sf <- st_read(shp_path, quiet = TRUE)
    canopy_sf <- st_transform(canopy_sf, crs = chm_crs)
    
    extracted <- exact_extract(
      rumble_raster,
      canopy_sf,
      include_cols = "Canopies",
      progress = FALSE
    )
    
    df <- map_df(extracted, function(x) {
      df <- filter(x, !is.na(value)) %>% select(Canopies, value)
      return(df)
    })
    
    pixels_list[[basename(shp_path)]] <- df
  }
  
  all_pixels <- bind_rows(pixels_list, .id = "Shapefile_Source") %>%
    separate(Canopies, into = c("TreeID", "SpeciesID"), sep = "_", remove = TRUE) %>%
    rename(Rumble = value)
  
  all_pixels$Source <- label
  return(all_pixels)
}

# ---- EXTRACT PIXELS FOR BOTH RUMBLES ----
rumble_pixels_on <- extract_rumble_pixels(rumble_on, shapefiles, "LeafOn")
rumble_pixels_off <- extract_rumble_pixels(rumble_off, shapefiles, "LeafOff")

# ---- ENTROPY FUNCTION ----
compute_entropy <- function(x, bins = 10) {
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}

# ---- SUMMARY FOR EACH ----
rumble_stats_on <- rumble_pixels_on %>%
  group_by(TreeID) %>%
  summarise(
    Max_Rumble_LeafOn = max(Rumble, na.rm = TRUE),
    Min_Rumble_LeafOn = min(Rumble, na.rm = TRUE),
    Mean_Rumble_LeafOn = mean(Rumble, na.rm = TRUE),
    SD_Rumble_LeafOn = sd(Rumble, na.rm = TRUE),
    Entropy_Rumble_LeafOn = compute_entropy(Rumble)
  )

rumble_stats_off <- rumble_pixels_off %>%
  group_by(TreeID) %>%
  summarise(
    Max_Rumble_LeafOff = max(Rumble, na.rm = TRUE),
    Min_Rumble_LeafOff = min(Rumble, na.rm = TRUE),
    Mean_Rumble_LeafOff = mean(Rumble, na.rm = TRUE),
    SD_Rumble_LeafOff = sd(Rumble, na.rm = TRUE),
    Entropy_Rumble_LeafOff = compute_entropy(Rumble)
  )

# ---- MERGE FINAL TABLE ----
canopy_rumble_stats <- full_join(rumble_stats_on, rumble_stats_off, by = "TreeID")

# ---- EXPORT ----
write.csv(canopy_rumble_stats, "E:/DATA/canopy_rumble_stats.csv", row.names = FALSE)
beep()
