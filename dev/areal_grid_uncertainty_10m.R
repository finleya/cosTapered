# Areal prediction uncertainty on a 10 m x 10 m grid.
#
# This script fits the example COS model and predicts latent areal quantities
# for a regular grid clipped to the forest/domain polygon. It is intended as a
# diagnostic for how areal prediction uncertainty changes near observed plots.

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
grid_size <- as.numeric(Sys.getenv("CTV_GRID_SIZE", "10"))

batch_length <- 25L
n_batch <- n_iter / batch_length
if (n_batch != as.integer(n_batch)) {
  stop("n_iter must be divisible by batch_length.")
}
n_batch <- as.integer(n_batch)

if (!is.finite(grid_size) || grid_size <= 0) {
  stop("CTV_GRID_SIZE must be positive.")
}

out_dir <- file.path("dev", "areal_grid_uncertainty_10m")
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
# Build 10 m grid clipped to domain
# ------------------------------------------------

domain <- st_union(st_geometry(forest_sf))
grid_geom <- st_make_grid(
  domain,
  cellsize = c(grid_size, grid_size),
  square = TRUE
)
grid_sf <- st_sf(grid_id = seq_along(grid_geom), geometry = grid_geom,
                 crs = st_crs(forest_sf))
grid_sf <- suppressWarnings(st_intersection(grid_sf, st_sf(geometry = domain)))
grid_sf <- st_make_valid(grid_sf)
grid_sf <- suppressWarnings(st_collection_extract(grid_sf, "POLYGON"))
grid_sf <- grid_sf[!st_is_empty(grid_sf), , drop = FALSE]
grid_sf$grid_id <- seq_len(nrow(grid_sf))

message("Prediction grid cells after clipping: ", nrow(grid_sf))

U_blocks <- cos_make_areal_blocks(
  U_sf = grid_sf,
  X_rast = X_rast,
  id_col = "grid_id",
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

pred_sf <- cbind(
  U_blocks$U_sf,
  pred_U$summary[, c("omega_mean", "omega_sd", "eta_mean", "eta_sd",
                     "omega_q025", "omega_q975", "eta_q025", "eta_q975")]
)
pred_df <- st_drop_geometry(pred_sf)

saveRDS(pred_U, file.path(out_dir, "areal_grid_prediction.rds"))
write.csv(pred_df, file.path(out_dir, "areal_grid_summary.csv"),
          row.names = FALSE)
st_write(pred_sf, file.path(out_dir, "areal_grid_supports.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# ------------------------------------------------
# Plot areal uncertainty and means
# ------------------------------------------------

base_map_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title.position = "plot"
  )

p_omega_sd <- ggplot() +
  geom_sf(data = pred_sf, aes(fill = omega_sd), color = NA) +
  geom_sf(data = B_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "grey10", linewidth = 0.35) +
  scale_fill_viridis_c(option = "cividis", name = "omega SD") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal spatial-effect uncertainty on 10 m grid") +
  base_map_theme

p_eta_sd <- ggplot() +
  geom_sf(data = pred_sf, aes(fill = eta_sd), color = NA) +
  geom_sf(data = B_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "grey10", linewidth = 0.35) +
  scale_fill_viridis_c(option = "magma", name = "eta SD") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal latent-response uncertainty on 10 m grid") +
  base_map_theme

p_eta_mean <- ggplot() +
  geom_sf(data = pred_sf, aes(fill = eta_mean), color = NA) +
  geom_sf(data = B_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "grey10", linewidth = 0.35) +
  scale_fill_viridis_c(option = "magma", name = "eta mean") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal latent-response posterior mean on 10 m grid") +
  base_map_theme

ggsave(file.path(out_dir, "areal_grid_omega_sd.png"),
       p_omega_sd, width = 7, height = 7, dpi = 300)
ggsave(file.path(out_dir, "areal_grid_eta_sd.png"),
       p_eta_sd, width = 7, height = 7, dpi = 300)
ggsave(file.path(out_dir, "areal_grid_eta_mean.png"),
       p_eta_mean, width = 7, height = 7, dpi = 300)

saveRDS(
  list(omega_sd = p_omega_sd, eta_sd = p_eta_sd, eta_mean = p_eta_mean),
  file.path(out_dir, "areal_grid_plots.rds")
)

message("Wrote outputs to: ", normalizePath(out_dir))
