
library(terra)
library(ggplot2)
library(dplyr)


seasonal_Occ <- rast("D:/LiDAR_Metrics/ALL_Metrics.tif")
plot(seasonal_Occ[[21]])

r <- seasonal_Occ[[21]]

# Convert raster to a data frame for ggplot
r_df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)

# Rename the value column to "value" for clarity
names(r_df)[3] <- "value"

# Create a custom color scale
color_palette <- scale_fill_gradientn(
  colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
  na.value = "white",
  name = "Value"
)

# Replace 0s with NA (they will appear white)
r_df <- r_df %>%
  mutate(value = ifelse(value == 0, NA, value))

# Plot with ggplot
ggplot(r_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  coord_equal() +
  color_palette +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")
