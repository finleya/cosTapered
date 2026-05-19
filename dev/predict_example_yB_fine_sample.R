# Predict fine-support observed y over the example study area.
#
# This script uses the installed cosTapered package, fits the example data for
# 10000 MCMC iterations, recovers 25 posterior samples after 1000 burn-in
# iterations, and makes pixel-level posterior predictive maps for y over a
# dense local subset around one or two observed plots.

library(sf)
library(raster)
library(ggplot2)
library(cosTapered)

set.seed(1)

# ------------------------------------------------
# Settings
# ------------------------------------------------

n_threads <- 25L
n_iter <- as.integer(Sys.getenv("CTV_N_ITER", "10000"))
burn_in <- as.integer(Sys.getenv("CTV_BURN_IN", "1000"))
n_pred_samples <- as.integer(Sys.getenv("CTV_N_PRED_SAMPLES", "25"))
pred_n_cells <- as.integer(Sys.getenv("CTV_PRED_N_CELLS", "12000"))
pred_plot_ids <- strsplit(Sys.getenv("CTV_PRED_PLOT_IDS", "1,2"), ",")[[1]]
pred_plot_ids <- as.integer(trimws(pred_plot_ids))
pred_plot_ids <- pred_plot_ids[is.finite(pred_plot_ids)]

batch_length <- 25L
n_batch <- n_iter / batch_length
if (n_batch != as.integer(n_batch)) {
  stop("n_iter must be divisible by batch_length.")
}
n_batch <- as.integer(n_batch)

if (!is.finite(pred_n_cells) || pred_n_cells < 1L) {
  stop("CTV_PRED_N_CELLS must be a positive integer.")
}

out_dir <- file.path("dev", "example_fine_yB_prediction_sample")
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

saveRDS(fit, file.path(out_dir, "fit.rds"))

rec_B <- cos_recover_B(
  fit = fit,
  burn_in = burn_in,
  n_samples = n_pred_samples,
  seed = 2,
  verbose = TRUE
)

saveRDS(rec_B, file.path(out_dir, "rec_B.rds"))

# ------------------------------------------------
# Predict fine-support y for a dense local subset
# ------------------------------------------------

X_values <- raster::getValues(X_rast)
if (is.null(dim(X_values))) {
  X_values <- matrix(X_values, ncol = 1)
}
X_values <- as.matrix(X_values)
colnames(X_values) <- names(X_rast)

pred_cells <- which(stats::complete.cases(X_values))
if (length(pred_plot_ids) == 0L ||
    any(pred_plot_ids < 1L | pred_plot_ids > nrow(B_sf))) {
  stop("CTV_PRED_PLOT_IDS must contain valid B_sf row numbers.")
}
pred_coords <- raster::xyFromCell(chm, pred_cells)
plot_centers <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(B_sf[pred_plot_ids, ])))

dist_sq <- matrix(NA_real_, nrow = nrow(pred_coords), ncol = nrow(plot_centers))
for (jj in seq_len(nrow(plot_centers))) {
  dist_sq[, jj] <- (pred_coords[, 1] - plot_centers[jj, 1])^2 +
    (pred_coords[, 2] - plot_centers[jj, 2])^2
}
nearest_dist_sq <- apply(dist_sq, 1, min)
keep_ord <- order(nearest_dist_sq)
n_keep <- min(pred_n_cells, length(keep_ord))
pred_cells <- pred_cells[keep_ord[seq_len(n_keep)]]
pred_coords <- pred_coords[keep_ord[seq_len(n_keep)], , drop = FALSE]
X_pred <- X_values[pred_cells, , drop = FALSE]

n_pred <- nrow(X_pred)
dense_gb <- 8 * n_pred^2 / 1024^3
message("Prediction cells: ", n_pred)
message("Prediction subset: nearest cells to B_sf row(s) ",
        paste(pred_plot_ids, collapse = ", "))
message("Dense prediction covariance is about ",
        round(dense_gb, 1), " GB per matrix before intermediates.")

pred_fine <- cos_predict_fine(
  fit = fit,
  rec_B = rec_B,
  pred_coords = pred_coords,
  X_pred = X_pred,
  method = "sample",
  keep_samples = TRUE,
  verbose = TRUE
)

eta_samples <- pred_fine$eta_samples
nugget_samples <- matrix(
  stats::rnorm(nrow(eta_samples) * ncol(eta_samples)),
  nrow = nrow(eta_samples),
  ncol = ncol(eta_samples)
)
nugget_samples <- sweep(
  nugget_samples,
  1,
  sqrt(rec_B$theta_samples$tau_sq),
  FUN = "*"
)
y_samples <- eta_samples + nugget_samples

pred_fine_df <- pred_fine$summary
pred_fine_df$cell_id <- pred_cells
pred_fine_df$y_mean <- colMeans(y_samples)
pred_fine_df$y_q025 <- apply(
  y_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE
)
pred_fine_df$y_q50 <- apply(
  y_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE
)
pred_fine_df$y_q975 <- apply(
  y_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE
)
pred_fine_df$y_ci_width_95 <- pred_fine_df$y_q975 - pred_fine_df$y_q025

saveRDS(pred_fine_df, file.path(out_dir, "fine_yB_prediction_summary.rds"))
write.csv(pred_fine_df,
          file.path(out_dir, "fine_yB_prediction_summary.csv"),
          row.names = FALSE)

# ------------------------------------------------
# Plot posterior predictive mean and 95 percent interval width
# ------------------------------------------------

base_map_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title.position = "plot",
    legend.position = "right"
  )

p_mean <- ggplot() +
  geom_tile(
    data = pred_fine_df,
    aes(x = x, y = y, fill = y_mean),
    width = xres(chm),
    height = yres(chm)
  ) +
  geom_sf(data = forest_sf, fill = NA, color = "grey10", linewidth = 0.45) +
  geom_sf(data = B_sf, fill = NA, color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "magma", name = "Mean") +
  coord_sf(expand = FALSE) +
  labs(
    title = "Fine-resolution posterior predictive mean",
    subtitle = "Observed fine-support target y, method = sample",
    x = NULL,
    y = NULL
  ) +
  base_map_theme

p_ci_width <- ggplot() +
  geom_tile(
    data = pred_fine_df,
    aes(x = x, y = y, fill = ometa_sd),
    width = xres(chm),
    height = yres(chm)
  ) +
  geom_sf(data = forest_sf, fill = NA, color = "grey10", linewidth = 0.45) +
  geom_sf(data = B_sf, fill = NA, color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "cividis", name = "95% CI width") +
  coord_sf(expand = FALSE) +
  labs(
    title = "Fine-resolution posterior predictive uncertainty",
    subtitle = "Width of the pixel-level 95% posterior predictive interval",
    x = NULL,
    y = NULL
  ) +
  base_map_theme


p_ci_width

ggsave(file.path(out_dir, "fine_yB_posterior_predictive_mean.png"),
       p_mean, width = 7, height = 7, dpi = 300)
ggsave(file.path(out_dir, "fine_yB_posterior_predictive_ci_width.png"),
       p_ci_width, width = 7, height = 7, dpi = 300)

saveRDS(
  list(mean = p_mean, ci_width = p_ci_width),
  file.path(out_dir, "fine_yB_prediction_plots.rds")
)

message("Wrote outputs to: ", normalizePath(out_dir))
