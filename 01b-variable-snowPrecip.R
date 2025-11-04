# -------------------------------
# Full Workflow: Cumulative Snow Precipitation (SWE) with Export
# -------------------------------

library(terra)
library(sf)
library(tigris)
library(ggplot2)
library(viridis)
library(dplyr)
library(nohrsc)

set.wd("C:/Users/gmjon/Documents/Git Local/michigan-dogman/")

# 0️⃣ Specify date range
start_date <- "2022-01-01"
end_date   <- "2022-04-01"

# 1️⃣ Local folder to save SNODAS files
snodas_folder <- "SNODAS"

# 2️⃣ Download snow precipitation for the period
snow_swe_list <- nsa_get_snodas(
  start.date = start_date,
  end.date   = end_date,
  product    = "snow precipitation",
  path       = snodas_folder
)

# 3️⃣ Find all .bil files for the period
bil_files <- list.files(snodas_folder, pattern = "\\.bil$", full.names = TRUE, recursive = TRUE)
if(length(bil_files) == 0) stop("No .bil files found in SNODAS folder.")

# 4️⃣ Load all rasters and stack them
rasters <- lapply(bil_files, rast)
swe_stack <- rast(rasters)

# 5️⃣ Compute cumulative snow precipitation (sum across days)
cum_swe <- app(swe_stack, fun = sum, na.rm = TRUE)

# 6️⃣ Reproject to Michigan Oblique Mercator (EPSG:3078)
cum_swe <- project(cum_swe, "EPSG:3078")

# 7️⃣ Get Michigan boundary in same CRS
mi <- tigris::states(cb = TRUE) %>%
  filter(NAME == "Michigan") %>%
  st_transform(crs(cum_swe))
mi_vect <- vect(mi)

# 8️⃣ Mask cumulative raster to Michigan only
cum_swe_mi <- mask(cum_swe, mi_vect)

# 9️⃣ Export the raster to disk
# Write out as GeoTIFF
writeRaster(cum_swe_mi, 
            filename = "demo_snowPrecip.tif", 
            overwrite = TRUE)

# 🔟 Convert to data.frame for ggplot
swe_df <- as.data.frame(cum_swe_mi, xy = TRUE)
names(swe_df)[3] <- "SWE"

# 1️⃣1️⃣ Plot
ggplot() +
  geom_tile(data = swe_df, aes(x = x, y = y, fill = SWE)) +
  geom_sf(data = mi, fill = NA, color = "black", size = 0.7) +
  scale_fill_viridis(option = "plasma", na.value = "grey90") +
  coord_sf(xlim = c(st_bbox(mi)$xmin, st_bbox(mi)$xmax),
           ylim = c(st_bbox(mi)$ymin, st_bbox(mi)$ymax)) +
  theme_minimal() +
  labs(title = paste0("Cumulative Snow Precipitation in Michigan: ", start_date, " to ", end_date),
       fill = "SWE (mm)")
