# Install these if you don't have them
# install.packages(c("sf", "maptiles", "tigris", "ggplot2", "ggspatial"))

library(sf)
library(ggplot2)
library(tigris)
library(maptiles)
library(ggspatial)

set.wd("C:/Users/gmjon/Documents/Git Local/michigan-dogman/")

# 1️⃣ Load your CSV
csv_path <- "demo_dogMan.csv"
points_df <- read.csv(csv_path, sep = ",")

# 2️⃣ Convert to sf (EPSG:3078)
points_sf <- st_as_sf(points_df, coords = c("x_coord", "y_coord"), crs = 3078)

# 3️⃣ Get Michigan boundary and match CRS
mi <- states(cb = TRUE) |>
  filter(NAME == "Michigan") |>
  st_transform(3078)

# 4️⃣ Get satellite basemap tiles (ESRI)
mi_bbox <- st_bbox(mi)
basemap <- get_tiles(mi_bbox, provider = "Esri.WorldImagery", zoom = 7, crop = TRUE)

# 5️⃣ Plot
a <- ggplot() +
  layer_spatial(basemap) +
  geom_sf(data = mi, fill = NA, color = "white", linewidth = 0.6) +
  geom_sf(data = points_sf, color = "#FFCD00", size = 2) +
  theme_minimal()

tiff("demo_pointMap.tif", res=300,
     width=8, height=7.5, units="in", compression="lzw")
a
dev.off()
