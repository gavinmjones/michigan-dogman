# ================================================================
# Logistic Regression with Quadratic Terms + Confidence Intervals
# ================================================================

library(sf)
library(terra)
library(tigris)
library(dplyr)
library(ggplot2)
library(patchwork)

set.wd("C:/Users/gmjon/Documents/Git Local/michigan-dogman/")

# -------------------------------
# 1️⃣ File paths
# -------------------------------
csv_path          <- "demo_dogMan3.csv"
ba_raster_path    <- "demo_basalAreaValue_smooth.tif"
snow_raster_path  <- "demo_snowPrecip.tif"

# -------------------------------
# 2️⃣ Load rasters
# -------------------------------
ba_rast   <- rast(ba_raster_path)
snow_rast <- rast(snow_raster_path)

# -------------------------------
# 3️⃣ Load observation points
# -------------------------------
points_df <- read.csv(csv_path, sep = ",")
points_sf <- st_as_sf(points_df, coords = c("x_coord", "y_coord"), crs = 3078)

# -------------------------------
# 4️⃣ Extract raster values to points
# -------------------------------
points_sf$BA   <- terra::extract(ba_rast, vect(points_sf))[, 2]
points_sf$Snow <- terra::extract(snow_rast, vect(points_sf))[, 2]

# -------------------------------
# 5️⃣ Generate random background (pseudo-absence) points
# -------------------------------
mi <- states(cb = TRUE) %>%
  filter(NAME == "Michigan") %>%
  st_transform(3078)

set.seed(123)
random_pts <- st_sample(mi, size = 600)
random_sf  <- st_sf(geometry = random_pts)

random_sf$BA   <- terra::extract(ba_rast, vect(random_sf))[, 2]
random_sf$Snow <- terra::extract(snow_rast, vect(random_sf))[, 2]

# Label presences (1) and absences (0)
points_sf$type <- 1
random_sf$type <- 0

# -------------------------------
# 6️⃣ Combine datasets
# -------------------------------
combined_sf <- bind_rows(
  points_sf %>% st_set_geometry(NULL),
  random_sf %>% st_set_geometry(NULL)
)

# Remove rows with missing data
combined_df_clean <- combined_sf %>%
  filter(!is.na(Snow), !is.na(BA))

# -------------------------------
# 7️⃣ Fit logistic regression with quadratic effects
# -------------------------------
logit_model <- glm(
  type ~ Snow + I(Snow^2) + BA,
  data = combined_df_clean,
  family = binomial
)

summary(logit_model)

# -------------------------------
# 8️⃣ Generate prediction data for plotting
# -------------------------------

# Predictor sequences
snow_seq <- seq(min(combined_df_clean$Snow, na.rm = TRUE),
                max(combined_df_clean$Snow, na.rm = TRUE),
                length.out = 200)
ba_seq <- seq(min(combined_df_clean$BA, na.rm = TRUE),
              max(combined_df_clean$BA, na.rm = TRUE),
              length.out = 200)

# Hold the other variable constant at its mean
snow_mean <- mean(combined_df_clean$Snow, na.rm = TRUE)
ba_mean   <- mean(combined_df_clean$BA, na.rm = TRUE)

# New data frames for predictions
new_snow <- data.frame(Snow = snow_seq, BA = ba_mean)
new_ba   <- data.frame(Snow = snow_mean, BA = ba_seq)

# Predict fitted probabilities and 95% confidence intervals
pred_snow <- predict(logit_model, newdata = new_snow, type = "link", se.fit = TRUE)
new_snow$fit  <- plogis(pred_snow$fit)
new_snow$lwr  <- plogis(pred_snow$fit - 1.96 * pred_snow$se.fit)
new_snow$upr  <- plogis(pred_snow$fit + 1.96 * pred_snow$se.fit)

pred_ba <- predict(logit_model, newdata = new_ba, type = "link", se.fit = TRUE)
new_ba$fit  <- plogis(pred_ba$fit)
new_ba$lwr  <- plogis(pred_ba$fit - 1.96 * pred_ba$se.fit)
new_ba$upr  <- plogis(pred_ba$fit + 1.96 * pred_ba$se.fit)

# -------------------------------
# 9️⃣ Plot with confidence intervals
# -------------------------------
library(ggplot2)
library(patchwork)
library(viridis)

set.seed(1223)
# --- Snow plot ---
a <- ggplot(combined_df_clean, aes(x = Snow, y = type)) +
  geom_jitter(aes(color = Snow), height = 0.02, alpha = 0.7) +
  geom_ribbon(
    data = new_snow,
    aes(x = Snow, ymin = lwr, ymax = upr),
    inherit.aes = FALSE, fill = "black", alpha = 0.2
  ) +
  geom_line(
    data = new_snow,
    aes(x = Snow, y = fit),
    inherit.aes = FALSE, color = "black", linewidth = 1.2
  ) +
  scale_color_viridis_c(option = "viridis", name = "Snow (mm)") +
  labs(
    y = "Predicted Probability of Michigan Dogman",
    x = "Snow Precipitation (mm)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# --- Basal Area plot ---
b <- ggplot(combined_df_clean, aes(x = BA, y = type)) +
  geom_jitter(aes(color = BA), height = 0.02, alpha = 0.7) +
  geom_ribbon(
    data = new_ba,
    aes(x = BA, ymin = lwr, ymax = upr),
    inherit.aes = FALSE, fill = "black", alpha = 0.2
  ) +
  geom_line(
    data = new_ba,
    aes(x = BA, y = fit),
    inherit.aes = FALSE, color = "black", linewidth = 1.2
  ) +
  scale_color_viridis_c(option = "plasma", name = "Basal Area (ft²/acre)") +
  labs(
    y = "Predicted Probability of Michigan Dogman",
    x = "Basal Area (ft²/acre)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# --- Combine plots ---
a + b





tiff("demo_logisticReg_additional.tif", res=300,
     width=13.3, height=7.5, units="in", compression="lzw")
a + b
dev.off()




# Ensure same extent and CRS
snow_rast <- resample(snow_rast, ba_rast, method = "bilinear")

# -----------------------------
# 2️⃣ Stack predictors
# -----------------------------
pred_stack <- c(ba_rast, snow_rast)
names(pred_stack) <- c("BA", "Snow")

# -----------------------------
# 3️⃣ Predict across raster
# -----------------------------
# logistic regression model already fitted as 'logit_model'

# terra::predict expects a SpatRaster with predictor names matching model
sdm_pred <- terra::predict(pred_stack, logit_model, type = "response")

# -----------------------------
# 4️⃣ Mask to Michigan state boundary
# -----------------------------
mi <- states(cb = TRUE) %>%
  filter(NAME == "Michigan") %>%
  st_transform(crs(ba_rast))
mi_vect <- vect(mi)

sdm_pred <- mask(sdm_pred, mi_vect)

# -----------------------------
# 5️⃣ Export SDM as GeoTIFF
# -----------------------------
writeRaster(sdm_pred, "demo_dogMan_SDM.tif",
            overwrite = TRUE)

# -----------------------------
# 6️⃣ Quick plot
# -----------------------------
sdm_df <- as.data.frame(sdm_pred, xy = TRUE, na.rm = TRUE)
colnames(sdm_df) <- c("x", "y", "prob")

tiff("C:/Users/gmjon/Desktop/MTU/demo_dogmanMap_v2.tif", res=300,
     width=8, height=7.5, units="in", compression="lzw")
ggplot() +
  # 1️⃣ Raster SDM
  geom_tile(data = sdm_df, aes(x = x, y = y, fill = prob)) +
  
  # 2️⃣ Michigan boundary
  geom_sf(data = mi, fill = NA, color = "black", size = 0.7) +
  
  
  # 4️⃣ Color scale for SDM
  scale_fill_viridis_c(option = "cividis", name = "P(Dogman)", trans = "sqrt") +
  
  coord_sf() +
  labs(x = "Easting", y = "Northing") +
  theme_minimal(base_size = 12)
dev.off()


tiff("demo_dogmanMap_points_v2.tif", res=300,
     width=8, height=7.5, units="in", compression="lzw")
ggplot() +
  # 1️⃣ Raster SDM
  geom_tile(data = sdm_df, aes(x = x, y = y, fill = prob)) +
  
  # 2️⃣ Michigan boundary
  geom_sf(data = mi, fill = NA, color = "black", size = 0.7) +
  
  # 3️⃣ Original point locations
  geom_sf(data = points_sf, size = 2, shape = 19, fill = "black") +
  
  # 4️⃣ Color scale for SDM
  scale_fill_viridis_c(option = "cividis", name = "P(Dogman)", trans = "sqrt") +
  
  coord_sf() +
  labs(x = "Easting", y = "Northing") +
  theme_minimal(base_size = 12)
dev.off()




# -----------------------------
# Install/load package
# -----------------------------
# install.packages("pROC")
library(pROC)

# -----------------------------
# 1️⃣ Prepare data
# -----------------------------
# 'type' is 1 for original points, 0 for random points
# 'prob_fit' contains model-predicted probability
eval_df <- combined_df_clean  # already has prob_fit

# -----------------------------
# 2️⃣ Compute ROC / AUC
# -----------------------------
roc_obj <- roc(response = eval_df$type,
               predictor = eval_df$prob_fit)

# AUC value
auc_val <- auc(roc_obj)
print(paste("AUC =", round(auc_val, 3)))

# -----------------------------
# 3️⃣ Plot ROC curve
# -----------------------------
plot(roc_obj, col = "blue", lwd = 2,
     main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))
