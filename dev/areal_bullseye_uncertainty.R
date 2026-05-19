# Areal uncertainty bullseye check around observed plots.
#
# This script uses the example data, fits the COS model, then predicts the
# observed plot supports and merged concentric distance bands around all
# observed plots. The goal is to check whether support-level uncertainty is
# smallest on/near the observed support and changes as prediction support moves
# away.

library(sf)
library(raster)
library(ggplot2)
library(cosTapered)

set.seed(1)

# ------------------------------------------------
# Settings
# ------------------------------------------------

n_threads <- as.integer(Sys.getenv("CTV_N_THREADS", "25"))
n_iter <- as.integer(Sys.getenv("CTV_N_ITER", "10000"))
burn_in <- as.integer(Sys.getenv("CTV_BURN_IN", "1000"))
n_pred_samples <- as.integer(Sys.getenv("CTV_N_PRED_SAMPLES", "25"))

plot_ids_env <- Sys.getenv("CTV_PLOT_IDS", "all")

ring_width <- as.numeric(Sys.getenv("CTV_RING_WIDTH", "10"))
max_radius <- as.numeric(Sys.getenv("CTV_MAX_RADIUS", "500"))

batch_length <- 25L
n_batch <- n_iter / batch_length
if (n_batch != as.integer(n_batch)) {
  stop("n_iter must be divisible by batch_length.")
}
n_batch <- as.integer(n_batch)

out_dir <- file.path("dev", "areal_bullseye_uncertainty")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Read example data
# ------------------------------------------------

chm <- raster(system.file("extdata", "example_chm.tif",
                          package = "cosTapered"))
names(chm) <- "chm"

B_sf <- st_read(system.file("extdata", "example_B.gpkg",
                            package = "cosTapered"), quiet = TRUE)
forest_sf <- st_read(system.file("extdata", "example_forest.gpkg",
                                 package = "cosTapered"), quiet = TRUE)
truth <- readRDS(system.file("extdata", "example_truth.rds",
                             package = "cosTapered"))

if (plot_ids_env == "all") {
  plot_ids <- seq_len(nrow(B_sf))
} else {
  plot_ids <- strsplit(plot_ids_env, ",")[[1]]
  plot_ids <- as.integer(trimws(plot_ids))
  plot_ids <- plot_ids[is.finite(plot_ids)]
  if (length(plot_ids) == 0L || any(plot_ids < 1L | plot_ids > nrow(B_sf))) {
    stop("CTV_PLOT_IDS must be 'all' or contain valid B_sf row numbers.")
  }
}

intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"

X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

# ------------------------------------------------
# Prepare, fit, and recover posterior samples
# ------------------------------------------------

prep <- cos_prepare(
  X_rast = X_rast,
  B_sf = B_sf,
  response_col = "y_B",
  gamma = truth$gamma,
  taper_code = truth$taper_code,
  n_threads = n_threads,
  verbose = TRUE
)

priors <- cos_default_priors(
  prep = prep,
  beta_mu = c(intercept = 40, chm = 1),
  beta_sd = c(intercept = 20, chm = 2),
  tau_B_shape = 2,
  tau_B_scale = truth$tau_B_sq_true,
  sigma_shape = 2,
  sigma_scale = truth$sigma_sq_true,
  phi_lower = 3 / 600,
  phi_upper = 3 / 80
)

fit <- cos_fit(
  prep = prep,
  priors = priors,
  n_chains = 1L,
  starting = data.frame(
    tau_B_sq = truth$tau_B_sq_true,
    sigma_sq = truth$sigma_sq_true,
    phi = truth$phi_true
  ),
  tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
  n_batch = n_batch,
  batch_length = batch_length,
  seed = 1,
  report = 100,
  verbose = TRUE
)

rec_B <- cos_recover_B(
  fit = fit,
  burn_in = burn_in,
  n_samples = n_pred_samples,
  seed = 2,
  verbose = TRUE
)

saveRDS(fit, file.path(out_dir, "fit.rds"))
saveRDS(rec_B, file.path(out_dir, "rec_B.rds"))

# ------------------------------------------------
# Build observed supports and merged concentric bands
# ------------------------------------------------

domain <- st_union(st_geometry(forest_sf))
plot_geom <- st_union(st_geometry(B_sf[plot_ids, ]))
plot_centers <- st_centroid(st_geometry(B_sf[plot_ids, ]))

supports <- list()

obs_geom <- suppressWarnings(st_intersection(plot_geom, domain))
supports[[1L]] <- st_sf(
  band_id = 0L,
  distance_min = 0,
  distance_max = NA_real_,
  distance_mid = 0,
  support = "observed_plots",
  geometry = obs_geom,
  crs = st_crs(B_sf)
)

previous <- plot_geom
radius_max <- seq(ring_width, max_radius, by = ring_width)

for (ii in seq_along(radius_max)) {
  outer <- st_union(st_buffer(plot_centers, dist = radius_max[ii]))
  ring <- suppressWarnings(st_difference(outer, previous))
  ring <- suppressWarnings(st_intersection(ring, domain))
  ring <- st_make_valid(ring)
  ring <- suppressWarnings(st_collection_extract(ring, "POLYGON"))
  ring <- st_union(ring)

  supports[[ii + 1L]] <- st_sf(
    band_id = ii,
    distance_min = radius_max[ii] - ring_width,
    distance_max = radius_max[ii],
    distance_mid = radius_max[ii] - ring_width / 2,
    support = paste0("band_", radius_max[ii] - ring_width, "_", radius_max[ii]),
    geometry = ring,
    crs = st_crs(B_sf)
  )

  previous <- st_union(outer)
}

U_sf <- do.call(rbind, supports)
U_sf <- U_sf[!st_is_empty(U_sf), , drop = FALSE]
U_sf$bullseye_id <- seq_len(nrow(U_sf))

message("Using ", length(plot_ids), " observed plot location(s).")
message("Built ", nrow(U_sf), " merged/clipped support band(s).")

U_blocks <- cos_make_areal_blocks(
  U_sf = U_sf,
  X_rast = X_rast,
  id_col = "bullseye_id",
  verbose = TRUE
)

pred_U <- cos_predict_areal(
  fit = fit,
  rec_B = rec_B,
  U_blocks = U_blocks,
  method = "sample",
  keep_samples = TRUE,
  verbose = TRUE
)

pred_df <- cbind(
  st_drop_geometry(U_blocks$U_sf),
  pred_U$summary[, c("omega_mean", "omega_sd", "eta_mean", "eta_sd",
                     "omega_q025", "omega_q975", "eta_q025", "eta_q975")]
)
pred_sf <- cbind(U_blocks$U_sf, pred_U$summary[, c("omega_sd", "eta_sd")])

saveRDS(pred_U, file.path(out_dir, "areal_bullseye_prediction.rds"))
write.csv(pred_df, file.path(out_dir, "areal_bullseye_summary.csv"),
          row.names = FALSE)
st_write(pred_sf, file.path(out_dir, "areal_bullseye_supports.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# ------------------------------------------------
# Plot uncertainty against distance and as maps
# ------------------------------------------------

p_distance <- ggplot(pred_df, aes(x = distance_mid)) +
  geom_line(aes(y = omega_sd, color = "omega"), linewidth = 0.7) +
  geom_point(aes(y = omega_sd, color = "omega"), size = 1.8) +
  geom_line(aes(y = eta_sd, color = "eta"), linewidth = 0.7) +
  geom_point(aes(y = eta_sd, color = "eta"), size = 1.8) +
  scale_color_manual(values = c(omega = "#0072B2", eta = "#D55E00"),
                     name = NULL) +
  labs(
    title = "Areal latent uncertainty by nearest-plot distance band",
    x = "Distance from observed plot centers",
    y = "Posterior SD"
  ) +
  theme_bw()

p_map_omega <- ggplot() +
  geom_sf(data = pred_sf, aes(fill = omega_sd), color = "white", linewidth = 0.15) +
  geom_sf(data = B_sf[plot_ids, ], fill = NA, color = "black", linewidth = 0.2) +
  geom_sf(data = forest_sf, fill = NA, color = "grey20", linewidth = 0.35) +
  scale_fill_viridis_c(option = "cividis", name = "omega SD") +
  coord_sf(expand = FALSE) +
  labs(title = "Spatial-effect uncertainty by support") +
  theme_bw()

p_map_eta <- ggplot() +
  geom_sf(data = pred_sf, aes(fill = eta_sd), color = "white", linewidth = 0.15) +
  geom_sf(data = B_sf[plot_ids, ], fill = NA, color = "black", linewidth = 0.2) +
  geom_sf(data = forest_sf, fill = NA, color = "grey20", linewidth = 0.35) +
  scale_fill_viridis_c(option = "magma", name = "eta SD") +
  coord_sf(expand = FALSE) +
  labs(title = "Latent response uncertainty by support") +
  theme_bw()

ggsave(file.path(out_dir, "areal_bullseye_sd_by_distance.png"),
       p_distance, width = 8, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "areal_bullseye_omega_sd_map.png"),
       p_map_omega, width = 7, height = 7, dpi = 300)
ggsave(file.path(out_dir, "areal_bullseye_eta_sd_map.png"),
       p_map_eta, width = 7, height = 7, dpi = 300)

saveRDS(
  list(distance = p_distance, omega_map = p_map_omega, eta_map = p_map_eta),
  file.path(out_dir, "areal_bullseye_plots.rds")
)

message("Wrote outputs to: ", normalizePath(out_dir))
