# -------------------------------
# Side-by-side plot of Snow Precip and Basal Area + Dog Man points
# -------------------------------

library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(tigris)
library(dplyr)
library(patchwork)  # for arranging plots side by side

set.wd("C:/Users/gmjon/Documents/Git Local/michigan-dogman/")

# 0️⃣ File paths
snow_file <- "demo_snowPrecip.tif"
ba_file   <- "demo_basalAreaValue.tif"
csv_path  <- "demo_dogMan2.csv"

# 1️⃣ Load rasters
snow_rast <- rast(snow_file)
ba_rast   <- rast(ba_file)

# 2️⃣ Load Michigan boundary
mi <- tigris::states(cb = TRUE) %>%
  filter(NAME == "Michigan")

# 3️⃣ Ensure all layers use same CRS (EPSG:3078)
mi_crs <- "EPSG:3078"

# Reproject Michigan boundary
mi <- st_transform(mi, crs = mi_crs)
mi_vect <- vect(mi)

# Assign CRS to rasters if missing
if(is.na(crs(snow_rast)) || crs(snow_rast) == ""){
  crs(snow_rast) <- mi_crs
}
if(is.na(crs(ba_rast)) || crs(ba_rast) == ""){
  crs(ba_rast) <- mi_crs
}

# 4️⃣ Resample BA raster to snow raster grid
ba_rast <- resample(ba_rast, snow_rast, method = "bilinear")

# 5️⃣ Mask both rasters to Michigan
snow_rast_mi <- mask(snow_rast, mi_vect)
ba_rast_mi   <- mask(ba_rast, mi_vect)

# 6️⃣ Convert rasters to data.frames for ggplot
snow_df <- as.data.frame(snow_rast_mi, xy = TRUE)
names(snow_df)[3] <- "SnowPrecip"

ba_df <- as.data.frame(ba_rast_mi, xy = TRUE)
names(ba_df)[3] <- "BasalArea"

# 7️⃣ Load and project Dog Man points
points_df <- read.csv(csv_path, sep = ",")
points_sf <- st_as_sf(points_df, coords = c("x_coord", "y_coord"), crs = 3078)
# (optional) make sure in same CRS
points_sf <- st_transform(points_sf, crs = mi_crs)

# 8️⃣ Create plots with point overlay
p1 <- ggplot() +
  geom_tile(data = snow_df, aes(x = x, y = y, fill = SnowPrecip)) +
  geom_sf(data = mi, fill = NA, color = "black", size = 0.7) +
  geom_sf(data = points_sf, shape = 21, fill = "white", color = "black", size = 2, alpha = 0.8) +
  scale_fill_viridis(option = "viridis", na.value = "grey90") +
  coord_sf(crs = st_crs(mi)) +
  theme_minimal() +
  labs(fill = "mm")

p2 <- ggplot() +
  geom_tile(data = ba_df, aes(x = x, y = y, fill = BasalArea)) +
  geom_sf(data = mi, fill = NA, color = "black", size = 0.7) +
  geom_sf(data = points_sf, shape = 21, fill = "white", color = "black", size = 2, alpha = 0.8) +
  scale_fill_viridis(option = "plasma", na.value = "grey90", trans = "sqrt") +
  coord_sf(crs = st_crs(mi)) +
  theme_minimal() +
  labs(fill = "ft²/acre")

# 9️⃣ Arrange side by side
tiff("demo_variables_points.tif", res=300,
     width=13.3, height=7.5, units="in", compression="lzw")
p1 + p2
dev.off()
