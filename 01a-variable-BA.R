# ============================================================
# Basal Area Rasterization + Gaussian Focal Smoothing (Michigan)
# ============================================================

library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(tigris)
library(dplyr)

set.wd("C:/Users/gmjon/Documents/Git Local/michigan-dogman/")
# -----------------------------
# 1️⃣ Load forest stands shapefile
# -----------------------------
shp <- st_read("AdminMiStateForest_4050457808582882637/Stands.shp")
shp_vect <- vect(shp)

# -----------------------------
# 2️⃣ Michigan boundary
# -----------------------------
mi <- states(cb = TRUE) |>
  filter(NAME == "Michigan") |>
  st_transform(st_crs(shp))
mi_vect <- vect(mi)

# -----------------------------
# 3️⃣ Raster template over Michigan
# -----------------------------
r_template <- rast(ext(mi_vect), resolution = 1000)
crs(r_template) <- crs(shp_vect)
values(r_template) <- NA_real_

# -----------------------------
# 4️⃣ Convert BARange categories to numeric
# -----------------------------
range_lookup <- data.frame(
  BARange = c("No Data","Unspecified","Immature","51-80","81-110",
              "111-140","141-170","171-200","201+"),
  min_val = c(1,1,1,51,81,111,141,171,200),
  max_val = c(50,50,50,80,110,140,170,200,300)
)

# Compute mean value per polygon for rasterization
shp$BA_numeric <- sapply(shp$BARange, function(x) {
  idx <- which(range_lookup$BARange == x)
  mean(c(range_lookup$min_val[idx], range_lookup$max_val[idx]))
})

# -----------------------------
# 5️⃣ Rasterize numeric basal area
# -----------------------------
r_ba <- rasterize(shp, r_template, field = "BA_numeric", fun = "mean")

# -----------------------------
# 6️⃣ Add variability within ranges
# -----------------------------
vals <- values(r_ba)
ba_numeric <- vals

for (i in seq_len(nrow(range_lookup))) {
  idx <- which(vals == mean(c(range_lookup$min_val[i], range_lookup$max_val[i])))
  if (length(idx) > 0) {
    ba_numeric[idx] <- runif(length(idx),
                             range_lookup$min_val[i],
                             range_lookup$max_val[i])
  }
}

# Strongly bias "1-50" range toward lower end
idx_small <- which(ba_numeric <= 50)
ba_numeric[idx_small] <- 1 + (50 - 1) * (runif(length(idx_small))^12)

values(r_ba) <- ba_numeric

# -----------------------------
# 6️⃣b Impute low basal area (1–50) for originally missing areas
# -----------------------------
vals <- values(r_ba)
na_idx <- which(is.na(vals))
vals[na_idx] <- 1 + (50 - 1) * (runif(length(na_idx))^12)
values(r_ba) <- vals

# -----------------------------
# 7️⃣ Mask raster to Michigan
# -----------------------------
r_ba <- mask(r_ba, mi_vect)

# -----------------------------
# 8️⃣ Gaussian focal smoothing (buffer-based edge correction)
# -----------------------------
# Create a 10 km buffer around Michigan boundary
mi_vect_buffer <- terra::buffer(mi_vect, width = 10000)

# Crop raster to buffered extent (less memory, preserves edges)
r_ba_crop <- terra::crop(r_ba, mi_vect_buffer)

# Create Gaussian weight matrix with 5 km radius
w <- terra::focalMat(r_ba_crop, d = 10000, type = "Gauss")

# Apply focal smoothing
r_ba_smooth_temp <- terra::focal(
  r_ba_crop,
  w = w,
  fun = mean,
  na.policy = "omit",
  na.rm = TRUE
)

# Mask final smoothed raster back to Michigan boundary
r_ba_smooth <- terra::mask(r_ba_smooth_temp, mi_vect)

# -----------------------------
# 9️⃣ Export GeoTIFF
# -----------------------------
writeRaster(
  r_ba_smooth,
  "demo_basalAreaValue_smooth.tif",
  overwrite = TRUE
)

# -----------------------------
# 🔟 Plot smoothed raster
# -----------------------------
r_df <- as.data.frame(r_ba_smooth, xy = TRUE, na.rm = TRUE)
colnames(r_df) <- c("x", "y", "BA")

ggplot(r_df) +
  geom_tile(aes(x = x, y = y, fill = BA)) +
  geom_sf(data = mi, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90",
                       name = "Basal Area (ft²/acre)",
                       trans = "sqrt") +
  coord_sf() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Smoothed Basal Area (1 km resolution, Gaussian kernel)",
    subtitle = "Full-state coverage with edge-corrected smoothing"
  )
