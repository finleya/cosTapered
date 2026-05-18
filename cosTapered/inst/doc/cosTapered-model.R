## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----range-bounds-------------------------------------------------------------
range_min <- 80
range_max <- 600
c(phi_lower = 3 / range_max, phi_upper = 3 / range_min)

## ----packages-----------------------------------------------------------------
library(sf)
library(raster)
library(ggplot2)
library(cosTapered)

## ----read-data----------------------------------------------------------------
chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"

B_sf <- st_read(system.file("extdata", "example_B.gpkg",
                            package = "cosTapered"), quiet = TRUE)
U_sf <- st_read(system.file("extdata", "example_U.gpkg",
                            package = "cosTapered"), quiet = TRUE)
forest_sf <- st_read(system.file("extdata", "example_forest.gpkg",
                                 package = "cosTapered"), quiet = TRUE)
truth <- readRDS(system.file("extdata", "example_truth.rds",
                             package = "cosTapered"))

truth[c("beta_true", "tau_B_sq_true", "sigma_sq_true",
        "eff_range_true", "phi_true", "gamma", "n_B", "n_U")]

## ----build-raster-design------------------------------------------------------
intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"

X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

## ----input-plots--------------------------------------------------------------
chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)
names(chm_df) <- c("x", "y", "chm")

ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = chm),
            width = xres(chm), height = yres(chm)) +
  geom_sf(data = U_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "black", linewidth = 0.55) +
  geom_sf(data = B_sf, fill = NA, color = "#c1121f", linewidth = 0.45) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Fine-support covariate and support polygons",
       fill = "CHM", x = NULL, y = NULL) +
  theme_bw()

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
          linewidth = 0.35) +
  geom_sf(data = B_sf, aes(fill = y_B), color = "black", linewidth = 0.25) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf(expand = FALSE) +
  labs(title = "Observed-support response", fill = "y_B",
       x = NULL, y = NULL) +
  theme_bw()

## ----prepare, eval=FALSE------------------------------------------------------
# prep <- cos_prepare(
#   X_rast = X_rast,
#   B_sf = B_sf,
#   response_col = "y_B",
#   gamma = truth$gamma,
#   taper_code = truth$taper_code,
#   n_threads = 1,
#   verbose = FALSE
# )
# 
# summary(prep)

## ----priors, eval=FALSE-------------------------------------------------------
# priors <- cos_default_priors(
#   prep = prep,
#   beta_mu = c(intercept = 40, chm = 1),
#   beta_sd = c(intercept = 20, chm = 2),
#   tau_B_shape = 2,
#   tau_B_scale = truth$tau_B_sq_true,
#   sigma_shape = 2,
#   sigma_scale = truth$sigma_sq_true,
#   phi_lower = 3 / 600,
#   phi_upper = 3 / 80
# )
# 
# priors

## ----fit, eval=FALSE----------------------------------------------------------
# fit <- cos_fit(
#   prep = prep,
#   priors = priors,
#   n_chains = 3,
#   starting = data.frame(
#     tau_B_sq = truth$tau_B_sq_true,
#     sigma_sq = truth$sigma_sq_true,
#     phi = truth$phi_true
#   ),
#   tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
#   n_batch = 50,
#   batch_length = 10,
#   seed = 1,
#   verbose = TRUE
# )
# 
# summary(fit, burn_in = 25)
# plot(fit, burn_in = 25)

## ----fit-multichain, eval=FALSE-----------------------------------------------
# fit <- cos_fit(
#   prep = prep,
#   priors = priors,
#   n_chains = 3,
#   starting = data.frame(
#     tau_B_sq = truth$tau_B_sq_true * c(0.5, 1, 2),
#     sigma_sq = truth$sigma_sq_true * c(0.5, 1, 2),
#     phi = truth$phi_true * c(0.75, 1, 1.25)
#   ),
#   #tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
#   n_batch = 100,
#   batch_length = 25,
#   seed = 1,
#   verbose = TRUE
# )
# 
# summary(fit, burn_in = 50)
# plot(fit, burn_in = 50)

## ----parameter-truth-plot, eval=FALSE-----------------------------------------
# theta <- fit$theta_samples
# theta <- theta[theta$iter >= 50, ]
# 
# theta_long <- reshape(
#   theta[, c("chain", "iter", "tau_B_sq", "sigma_sq", "phi", "eff_range")],
#   varying = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
#   v.names = "value",
#   timevar = "parameter",
#   times = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
#   direction = "long"
# )
# 
# truth_df <- data.frame(
#   parameter = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
#   truth = c(truth$tau_B_sq_true, truth$sigma_sq_true,
#             truth$phi_true, truth$eff_range_true)
# )
# 
# ggplot(theta_long, aes(x = value)) +
#   geom_density(fill = "grey80", color = "grey25", linewidth = 0.35) +
#   geom_vline(data = truth_df, aes(xintercept = truth),
#              color = "#c1121f", linewidth = 0.7) +
#   facet_wrap(~ parameter, scales = "free", ncol = 2) +
#   labs(title = "Posterior distributions and known simulation values",
#        x = NULL, y = NULL) +
#   theme_bw()

## ----recover, eval=FALSE------------------------------------------------------
# rec_B <- cos_recover_B(
#   fit = fit,
#   burn_in = 50,
#   n_samples = 100,
#   seed = 1,
#   verbose = FALSE
# )
# 
# rec_sum <- summary(rec_B, include_eta_B = TRUE)
# rec_sum

## ----recover-plots, eval=FALSE------------------------------------------------
# B_plot <- B_sf
# B_plot$eta_B_mean <- rec_sum$eta_B_summary$mean[
#   match(B_plot$B_id, rec_sum$eta_B_summary$B_id)
# ]
# B_plot$eta_B_error <- B_plot$eta_B_mean - B_plot$eta_B_true
# 
# ggplot(B_plot, aes(x = eta_B_true, y = eta_B_mean)) +
#   geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.4) +
#   geom_point(size = 2) +
#   coord_equal() +
#   labs(title = "Recovered observed-support latent means",
#        x = "Simulated eta_B truth",
#        y = "Posterior mean") +
#   theme_bw()
# 
# ggplot() +
#   geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
#           linewidth = 0.35) +
#   geom_sf(data = B_plot, aes(fill = eta_B_error),
#           color = "black", linewidth = 0.25) +
#   scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1121f") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Fitted eta_B posterior mean minus truth",
#        fill = "error", x = NULL, y = NULL) +
#   theme_bw()

## ----predict-fine, eval=FALSE-------------------------------------------------
# pred_fine <- cos_predict_fine(
#   fit = fit,
#   rec_B = rec_B,
#   method = "mean",
#   verbose = FALSE
# )
# 
# summary(pred_fine)

## ----predict-fine-plots, eval=FALSE-------------------------------------------
# fine_df <- pred_fine$summary[, c("x", "y", "eta_mean", "eta_sd")]
# 
# ggplot() +
#   geom_tile(data = fine_df, aes(x = x, y = y, fill = eta_mean),
#             width = xres(chm), height = yres(chm)) +
#   geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.25) +
#   geom_sf(data = forest_sf, fill = NA, color = "white", linewidth = 0.35) +
#   scale_fill_viridis_c(option = "magma") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Fine-support posterior mean",
#        fill = "eta mean", x = NULL, y = NULL) +
#   theme_bw()
# 
# ggplot() +
#   geom_tile(data = fine_df, aes(x = x, y = y, fill = eta_sd),
#             width = xres(chm), height = yres(chm)) +
#   geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.25) +
#   geom_sf(data = forest_sf, fill = NA, color = "white", linewidth = 0.35) +
#   scale_fill_viridis_c(option = "cividis") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Fine-support posterior SD",
#        fill = "eta SD", x = NULL, y = NULL) +
#   theme_bw()

## ----areal-blocks, eval=FALSE-------------------------------------------------
# U_blocks <- cos_make_areal_blocks(
#   U_sf = U_sf,
#   X_rast = X_rast,
#   id_col = "U_id",
#   verbose = FALSE
# )
# 
# summary(U_blocks)

## ----predict-areal, eval=FALSE------------------------------------------------
# pred_U <- cos_predict_areal(
#   fit = fit,
#   rec_B = rec_B,
#   U_blocks = U_blocks,
#   method = "mean",
#   verbose = TRUE
# )
# 
# summary(pred_U)

## ----predict-areal-plots, eval=FALSE------------------------------------------
# U_plot <- pred_U$U_blocks$U_sf
# U_plot$eta_mean <- pred_U$summary$eta_mean[
#   match(U_plot$U_id, pred_U$summary$U_id)
# ]
# U_plot$eta_sd <- pred_U$summary$eta_sd[
#   match(U_plot$U_id, pred_U$summary$U_id)
# ]
# 
# ggplot() +
#   geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
#           linewidth = 0.35) +
#   geom_sf(data = U_plot, aes(fill = eta_mean),
#           color = "white", linewidth = 0.25) +
#   geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.2) +
#   scale_fill_viridis_c(option = "magma") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Areal posterior mean",
#        fill = "eta mean", x = NULL, y = NULL) +
#   theme_bw()
# 
# ggplot() +
#   geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
#           linewidth = 0.35) +
#   geom_sf(data = U_plot, aes(fill = eta_sd),
#           color = "white", linewidth = 0.25) +
#   geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.2) +
#   scale_fill_viridis_c(option = "cividis") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Areal posterior SD",
#        fill = "eta SD", x = NULL, y = NULL) +
#   theme_bw()

## ----cv, eval=FALSE-----------------------------------------------------------
# set.seed(1)
# fold_id <- sample(rep(seq_len(10), length.out = nrow(B_sf)))
# 
# cv <- cos_cv_observed(
#   fit = fit,
#   X_rast = X_rast,
#   B_sf = B_sf,
#   target = c("latent", "observed"),
#   latent_col = "eta_B_true",
#   k = 10,
#   fold_id = fold_id,
#   fit_args = list(n_batch = 50, batch_length = 8),
#   recover_args = list(burn_in = 100, n_samples = 60),
#   seed = 1,
#   verbose = FALSE
# )
# 
# summary(cv)
# cv$fold_summary

## ----cv-plot, eval=FALSE------------------------------------------------------
# ggplot(cv$predictions, aes(x = truth, y = pred_mean)) +
#   geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.4) +
#   geom_errorbar(aes(ymin = q025, ymax = q975),
#                 width = 0, alpha = 0.55) +
#   geom_point(aes(color = covered_95), size = 2) +
#   facet_wrap(~ target, scales = "free") +
#   labs(title = "Held-out observed-support predictions",
#        x = "Truth", y = "Posterior predictive mean",
#        color = "Covered") +
#   theme_bw()

## ----nonspatial, eval=FALSE---------------------------------------------------
# prep_ns <- cos_prepare(
#   X_rast = X_rast,
#   B_sf = B_sf,
#   response_col = "y_B",
#   spatial = FALSE,
#   verbose = FALSE
# )
# 
# priors_ns <- cos_default_priors(
#   prep = prep_ns,
#   beta_mu = c(intercept = 40, chm = 1),
#   beta_sd = c(intercept = 20, chm = 2),
#   tau_B_shape = 2,
#   tau_B_scale = truth$tau_B_sq_true
# )
# 
# fit_ns <- cos_fit(
#   prep = prep_ns,
#   priors = priors_ns,
#   starting = c(tau_B_sq = truth$tau_B_sq_true),
#   tuning = c(log_tau_B_sq = 0.25),
#   n_batch = 50,
#   batch_length = 8,
#   seed = 1,
#   verbose = FALSE
# )
# 
# summary(fit_ns, burn_in = 100)
# plot(fit_ns, burn_in = 100)

