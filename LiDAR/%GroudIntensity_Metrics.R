

library(lidR)
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(entropy)
library(beepr)


# ---- Paths ----
canopies_path <- "E:/Git Paint Rock 1.0/Updated Canopy Polygons/Updated/"
LAZon_path <- "E:/Updated LiDAR/PaintRock_20ha_leafOn_class.laz"
LAZoff_path <- "E:/Updated LiDAR/PaintRock_20ha_leafOff_class.laz"

# ---- Function: Compute Entropy ----
compute_entropy <- function(x, bins = 10) {
  if (length(x) == 0 || all(is.na(x))) return(NA)
  counts <- hist(x, breaks = bins, plot = FALSE)$counts
  if (sum(counts) == 0) return(NA)
  entropy.empirical(counts, unit = "log2")
}


# ---- Function: Create % Ground Intensity Raster ----
create_ground_intensity_pct_raster <- function(laz_path, res = 1) {
  las <- readLAS(laz_path)
  if (is.empty(las)) stop("LAS file is empty.")
  
  cat("Rasterizing total intensity...\n")
  all_intensity <- grid_metrics(las, ~sum(Intensity, na.rm = TRUE), res = res)
  
  cat("Rasterizing ground intensity...\n")
  ground_las <- filter_poi(las, Classification == 2)
  ground_intensity <- grid_metrics(ground_las, ~sum(Intensity, na.rm = TRUE), res = res)
  
  # Align rasters
  all_intensity <- resample(all_intensity, ground_intensity, method = "bilinear")
  
  # Compute % ground intensity
  pct_ground <- (ground_intensity / all_intensity) * 100
  names(pct_ground) <- "Pct_Ground_Intensity"
  
  return(pct_ground)
}

# ---- Create rasters ----
pct_ground_raster_on <- create_ground_intensity_pct_raster(LAZon_path)
pct_ground_raster_off <- create_ground_intensity_pct_raster(LAZoff_path)
beep()

# ---- Load canopy polygons ----
shapefiles <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)
canopies_list <- lapply(shapefiles, st_read, quiet = TRUE)
canopies_sf <- do.call(rbind, canopies_list)

# ---- Ensure CRS match ----
canopies_sf <- st_transform(canopies_sf, crs = st_crs(pct_ground_raster_on))

# ---- Extract stats from LeafOn ----
extract_on <- exact_extract(pct_ground_raster_on, canopies_sf, fun = function(values, cov) {
  values <- values[!is.na(values)]
  if(length(values) == 0) return(c(Mean = NA, Max = NA, SD = NA, Entropy = NA))
  
  list(
    Mean = mean(values),
    Max = max(values),
    SD = sd(values),
    Entropy = compute_entropy(values)
  )
}, progress = FALSE)

# ---- Extract stats from LeafOff ----
extract_off <- exact_extract(pct_ground_raster_off, canopies_sf, fun = function(values, cov) {
  values <- values[!is.na(values)]
  if(length(values) == 0) return(c(Mean = NA, Max = NA, SD = NA, Entropy = NA))
  
  list(
    Mean = mean(values),
    Max = max(values),
    SD = sd(values),
    Entropy = compute_entropy(values)
  )
}, progress = FALSE)

table(sapply(extract_on, length))
table(sapply(extract_off, length))


# ---- Combine with TreeID ----
canopy_ids <- canopies_sf$Canopies
canopy_stats <- data.frame(
  TreeID = gsub("_.*", "", canopy_ids),
  SpeciesID = gsub(".*_", "", canopy_ids),
  
  Mean_PctGround_LeafOn = sapply(extract_on, `[[`, 1),
  Max_PctGround_LeafOn = sapply(extract_on, `[[`, 2),
  SD_PctGround_LeafOn = sapply(extract_on, `[[`, 3),
  Entropy_PctGround_LeafOn = sapply(extract_on, `[[`, 4),
  
  Mean_PctGround_LeafOff = sapply(extract_off, `[[`, 1),
  Max_PctGround_LeafOff = sapply(extract_off, `[[`, 2),
  SD_PctGround_LeafOff = sapply(extract_off, `[[`, 3),
  Entropy_PctGround_LeafOff = sapply(extract_off, `[[`, 4)
)

# ---- Export ----
write.csv(canopy_stats, "canopy_ground_intensity_pct_stats.csv", row.names = FALSE)
