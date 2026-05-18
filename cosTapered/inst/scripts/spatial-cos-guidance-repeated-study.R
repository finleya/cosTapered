set.seed(20260517)

library(sf)
library(raster)

pkg_dir <- "cosTapered"
use_pkgload <- identical(Sys.getenv("CTV_USE_PKGLOAD", unset = "FALSE"), "TRUE")
if (use_pkgload && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_dir, quiet = TRUE, recompile = TRUE)
} else {
  library(cosTapered)
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) return(default)
  as.integer(value)
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) return(default)
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n", sep = "")
}

summarize_metric <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q10 = unname(stats::quantile(x, 0.10, na.rm = TRUE)),
    q50 = unname(stats::quantile(x, 0.50, na.rm = TRUE)),
    q90 = unname(stats::quantile(x, 0.90, na.rm = TRUE))
  )
}

summarize_groups <- function(dat, by, metrics) {
  groups <- split(dat, dat[by], drop = TRUE)
  rows <- lapply(groups, function(z) {
    out <- z[1, by, drop = FALSE]
    for (metric in metrics) {
      vals <- summarize_metric(z[[metric]])
      for (stat in names(vals)) {
        out[[paste(metric, stat, sep = "_")]] <- vals[[stat]]
      }
    }
    out$n_reps <- length(unique(z$rep_id))
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

make_improvement <- function(results) {
  wide <- reshape(
    results[, c(
      "rep_id", "target", "eff_range", "model", "RMSPE", "MAE", "CRPS",
      "coverage_95", "mean_interval_width"
    )],
    idvar = c("rep_id", "target", "eff_range"),
    timevar = "model",
    direction = "wide"
  )
  names(wide) <- sub("^RMSPE\\.", "RMSPE_", names(wide))
  names(wide) <- sub("^MAE\\.", "MAE_", names(wide))
  names(wide) <- sub("^CRPS\\.", "CRPS_", names(wide))
  names(wide) <- sub("^coverage_95\\.", "coverage_95_", names(wide))
  names(wide) <- sub("^mean_interval_width\\.", "mean_interval_width_", names(wide))
  wide$RMSPE_improvement_pct <-
    100 * (wide$RMSPE_nonspatial - wide$RMSPE_spatial) / wide$RMSPE_nonspatial
  wide$CRPS_improvement_pct <-
    100 * (wide$CRPS_nonspatial - wide$CRPS_spatial) / wide$CRPS_nonspatial
  wide$spatial_wins_RMSPE <- wide$RMSPE_spatial < wide$RMSPE_nonspatial
  wide$spatial_wins_CRPS <- wide$CRPS_spatial < wide$CRPS_nonspatial
  wide
}

write_outputs <- function(results, predictions, out_dir) {
  improvement <- make_improvement(results)

  metric_cols <- c(
    "RMSPE", "MAE", "bias", "cor", "slope", "CRPS", "coverage_95",
    "mean_interval_width", "tau_B_sq_median", "sigma_sq_median",
    "eff_range_median"
  )
  improvement_cols <- c(
    "RMSPE_spatial", "RMSPE_nonspatial", "RMSPE_improvement_pct",
    "CRPS_spatial", "CRPS_nonspatial", "CRPS_improvement_pct",
    "coverage_95_spatial", "coverage_95_nonspatial",
    "mean_interval_width_spatial", "mean_interval_width_nonspatial",
    "spatial_wins_RMSPE", "spatial_wins_CRPS"
  )

  mean_summary <- summarize_groups(
    results,
    by = c("target", "model", "eff_range"),
    metrics = metric_cols
  )
  mean_improvement <- summarize_groups(
    improvement,
    by = c("target", "eff_range"),
    metrics = improvement_cols
  )

  write.csv(results, file.path(out_dir, "cv_range_targets_replicate_summary.csv"),
            row.names = FALSE)
  write.csv(predictions, file.path(out_dir, "cv_range_targets_predictions.csv"),
            row.names = FALSE)
  write.csv(improvement, file.path(out_dir, "cv_range_targets_replicate_improvement.csv"),
            row.names = FALSE)
  write.csv(mean_summary, file.path(out_dir, "cv_range_targets_mean_summary.csv"),
            row.names = FALSE)
  write.csv(mean_improvement, file.path(out_dir, "cv_range_targets_mean_improvement.csv"),
            row.names = FALSE)

  target_names <- c("latent_eta_B", "observed_y_B")
  target_labels <- c("Latent eta_B", "Observed y_B")
  model_cols <- c(spatial = "#2563eb", nonspatial = "#dc2626")

  png(file.path(out_dir, "cv_range_targets_mean_rmspe_crps.png"),
      width = 1000, height = 700, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (target_i in seq_along(target_names)) {
    target <- target_names[target_i]
    dat <- mean_summary[mean_summary$target == target, ]
    ylim <- range(dat$RMSPE_mean, finite = TRUE)
    plot(NA, xlim = range(dat$eff_range), ylim = ylim,
         xlab = "True effective range", ylab = "Mean RMSPE",
         main = target_labels[target_i])
    for (model in c("spatial", "nonspatial")) {
      sub <- dat[dat$model == model, ]
      lines(sub$eff_range, sub$RMSPE_mean, type = "b", pch = 19, lwd = 2,
            col = model_cols[[model]])
    }
    if (target_i == 1L) {
      legend("topleft", legend = c("spatial COS", "non-spatial COS"),
             col = model_cols, pch = 19, lwd = 2, bty = "n")
    }
  }
  for (target_i in seq_along(target_names)) {
    target <- target_names[target_i]
    dat <- mean_summary[mean_summary$target == target, ]
    ylim <- range(dat$CRPS_mean, finite = TRUE)
    plot(NA, xlim = range(dat$eff_range), ylim = ylim,
         xlab = "True effective range", ylab = "Mean CRPS",
         main = target_labels[target_i])
    for (model in c("spatial", "nonspatial")) {
      sub <- dat[dat$model == model, ]
      lines(sub$eff_range, sub$CRPS_mean, type = "b", pch = 19, lwd = 2,
            col = model_cols[[model]])
    }
  }
  par(op)
  dev.off()

  png(file.path(out_dir, "cv_range_targets_mean_improvement.png"),
      width = 1000, height = 700, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (metric in c("RMSPE_improvement_pct", "CRPS_improvement_pct",
                   "spatial_wins_RMSPE", "spatial_wins_CRPS")) {
    ylim <- range(mean_improvement[[paste0(metric, "_mean")]], finite = TRUE)
    if (grepl("spatial_wins", metric)) ylim <- c(0, 1)
    plot(NA, xlim = range(mean_improvement$eff_range), ylim = ylim,
         xlab = "True effective range", ylab = metric, main = metric)
    abline(h = if (grepl("spatial_wins", metric)) 0.5 else 0,
           lty = 2, col = "grey45")
    for (target_i in seq_along(target_names)) {
      target <- target_names[target_i]
      sub <- mean_improvement[mean_improvement$target == target, ]
      lines(sub$eff_range, sub[[paste0(metric, "_mean")]], type = "b",
            pch = 19, lwd = 2, col = c("#111827", "#059669")[target_i])
    }
    if (metric == "RMSPE_improvement_pct") {
      legend("topleft", legend = target_labels,
             col = c("#111827", "#059669"), pch = 19, lwd = 2, bty = "n")
    }
  }
  par(op)
  dev.off()

  list(
    improvement = improvement,
    mean_summary = mean_summary,
    mean_improvement = mean_improvement
  )
}

out_dir <- Sys.getenv(
  "CTV_OUT_DIR",
  unset = "spatial-cos-guidance-output"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(out_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

msg("Reading example raster and forest boundary")
chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"
forest_sf <- st_read(
  system.file("extdata", "example_forest.gpkg", package = "cosTapered"),
  quiet = TRUE
)
forest_sf <- st_transform(forest_sf, st_crs(crs(chm)))

intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"
X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

beta_true <- c(intercept = -10, chm = 15.87132)
tau_B_sq_true <- 250
sigma_sq_true <- 750
eff_range_grid <- env_num_vec(
  "CTV_EFF_RANGES",
  c(25, 50, 100, 150, 250, 350, 500, 750, 1000)
)
n_reps <- env_int("CTV_N_REPS", 50L)
n_threads <- env_int("CTV_N_THREADS", 25L)
k_folds <- env_int("CTV_K_FOLDS", 10L)

n_chains <- env_int("CTV_N_CHAINS", 1L)
n_batch <- env_int("CTV_N_BATCH", 60L)
batch_length <- env_int("CTV_BATCH_LENGTH", 10L)
burn_in <- env_int("CTV_BURN_IN", 150L)
thin <- env_int("CTV_THIN", 2L)
n_samples <- env_int("CTV_N_SAMPLES", 100L)
n_B <- env_int("CTV_N_B", 49L)
plot_radius <- env_int("CTV_PLOT_RADIUS", 16L)
design_seed <- env_int("CTV_DESIGN_SEED", 1101L)
fold_seed <- env_int("CTV_FOLD_SEED", 3303L)
sim_seed_base <- env_int("CTV_SIM_SEED_BASE", 220200L)

run_settings <- list(
  beta_true = beta_true,
  tau_B_sq_true = tau_B_sq_true,
  sigma_sq_true = sigma_sq_true,
  eff_range_grid = eff_range_grid,
  gamma_rule = "effective range + 100",
  n_reps = n_reps,
  n_threads = n_threads,
  k_folds = k_folds,
  n_chains = n_chains,
  n_batch = n_batch,
  batch_length = batch_length,
  burn_in = burn_in,
  thin = thin,
  n_samples = n_samples,
  n_B = n_B,
  plot_radius = plot_radius,
  design_seed = design_seed,
  fold_seed = fold_seed,
  sim_seed_base = sim_seed_base
)
saveRDS(run_settings, file.path(out_dir, "cv_range_targets_settings.rds"))

results <- list()
predictions <- list()
row_id <- 1L
fold_id_fixed <- NULL

for (rep_id in seq_len(n_reps)) {
  for (rr in seq_along(eff_range_grid)) {
    eff_range <- eff_range_grid[rr]
    phi <- 3 / eff_range
    gamma <- eff_range + 100
    sim_seed <- sim_seed_base + rep_id
    checkpoint_file <- file.path(
      checkpoint_dir,
      sprintf("rep_%03d_range_%04d_checkpoint.rds", rep_id, eff_range)
    )

    msg("Replicate ", rep_id, " of ", n_reps,
        "; effective range = ", eff_range,
        "; gamma = ", gamma,
        "; sim_seed = ", sim_seed)

    sim <- cos_simulate_data(
      X_rast = X_rast,
      forest_sf = forest_sf,
      n_B = n_B,
      plot_radius = plot_radius,
      gamma = gamma,
      beta = beta_true,
      tau_B_sq = tau_B_sq_true,
      sigma_sq = sigma_sq_true,
      phi = phi,
      taper_code = 1L,
      design_seed = design_seed,
      sim_seed = sim_seed,
      n_threads = n_threads,
      verbose = FALSE
    )
    B_sf <- sim$B_sf

    if (is.null(fold_id_fixed)) {
      set.seed(fold_seed)
      fold_id_fixed <- sample(rep(seq_len(k_folds), length.out = nrow(B_sf)))
      saveRDS(fold_id_fixed, file.path(out_dir, "cv_range_targets_fold_id.rds"))
    }
    if (length(fold_id_fixed) != nrow(B_sf)) {
      stop(
        "The simulated number of observed-support polygons changed from ",
        length(fold_id_fixed), " to ", nrow(B_sf),
        ". Use a fixed design setup before running this repeated study."
      )
    }

    prep_sp <- cos_prepare(
      X_rast = X_rast,
      B_sf = B_sf,
      response_col = "y_B",
      gamma = gamma,
      spatial = TRUE,
      n_threads = n_threads,
      verbose = FALSE
    )
    priors_sp <- cos_default_priors(
      prep = prep_sp,
      beta_mu = beta_true,
      beta_sd = c(intercept = 60, chm = 4),
      tau_B_shape = 2,
      tau_B_scale = tau_B_sq_true,
      sigma_shape = 2,
      sigma_scale = sigma_sq_true,
      phi_lower = 3 / 1500,
      phi_upper = 3 / 10
    )
    fit_sp <- cos_fit(
      prep = prep_sp,
      priors = priors_sp,
      n_chains = n_chains,
      starting = data.frame(
        tau_B_sq = tau_B_sq_true,
        sigma_sq = sigma_sq_true,
        phi = phi
      ),
      tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.35),
      n_batch = n_batch,
      batch_length = batch_length,
      seed = 100000 + 1000 * rep_id + rr,
      report = 999,
      verbose = FALSE
    )

    prep_ns <- cos_prepare(
      X_rast = X_rast,
      B_sf = B_sf,
      response_col = "y_B",
      spatial = FALSE,
      n_threads = n_threads,
      verbose = FALSE
    )
    priors_ns <- cos_default_priors(
      prep = prep_ns,
      beta_mu = beta_true,
      beta_sd = c(intercept = 60, chm = 4),
      tau_B_shape = 2,
      tau_B_scale = stats::var(B_sf$y_B)
    )
    fit_ns <- cos_fit(
      prep = prep_ns,
      priors = priors_ns,
      n_chains = n_chains,
      starting = c(tau_B_sq = stats::var(B_sf$y_B)),
      tuning = c(log_tau_B_sq = 0.25),
      n_batch = n_batch,
      batch_length = batch_length,
      seed = 200000 + 1000 * rep_id + rr,
      report = 999,
      verbose = FALSE
    )

    fit_args_sp <- list(
      n_chains = n_chains,
      n_batch = n_batch,
      batch_length = batch_length,
      tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.35),
      report = 999
    )
    fit_args_ns <- list(
      n_chains = n_chains,
      n_batch = n_batch,
      batch_length = batch_length,
      tuning = c(log_tau_B_sq = 0.25),
      report = 999
    )
    recover_args <- list(burn_in = burn_in, thin = thin, n_samples = n_samples)

    cv_sp <- cos_cv_observed(
      fit = fit_sp,
      X_rast = X_rast,
      B_sf = B_sf,
      response_col = "y_B",
      target = c("latent", "observed"),
      latent_col = "eta_B_true",
      k = k_folds,
      fold_id = fold_id_fixed,
      fit_args = fit_args_sp,
      recover_args = recover_args,
      seed = 300000 + 1000 * rep_id + rr,
      verbose = FALSE
    )
    cv_ns <- cos_cv_observed(
      fit = fit_ns,
      X_rast = X_rast,
      B_sf = B_sf,
      response_col = "y_B",
      target = c("latent", "observed"),
      latent_col = "eta_B_true",
      k = k_folds,
      fold_id = fold_id_fixed,
      fit_args = fit_args_ns,
      recover_args = recover_args,
      seed = 400000 + 1000 * rep_id + rr,
      verbose = FALSE
    )

    for (model in c("spatial", "nonspatial")) {
      cv <- if (model == "spatial") cv_sp else cv_ns
      for (target in cv$summary$target) {
        target_label <- if (target == "latent") "latent_eta_B" else "observed_y_B"
        summary_row <- cv$summary[cv$summary$target == target, , drop = FALSE]
        theta <- if (model == "spatial") fit_sp$theta_samples else fit_ns$theta_samples
        theta <- theta[theta$iter >= burn_in, , drop = FALSE]
        results[[row_id]] <- data.frame(
          rep_id = rep_id,
          target = target_label,
          model = model,
          eff_range = eff_range,
          summary_row[, setdiff(names(summary_row), "target"), drop = FALSE],
          tau_B_sq_median = stats::median(theta$tau_B_sq),
          sigma_sq_median = if (all(is.na(theta$sigma_sq))) {
            NA_real_
          } else {
            stats::median(theta$sigma_sq, na.rm = TRUE)
          },
          eff_range_median = if (all(is.na(theta$eff_range))) {
            NA_real_
          } else {
            stats::median(theta$eff_range, na.rm = TRUE)
          },
          stringsAsFactors = FALSE
        )

        pred <- cv$predictions[cv$predictions$target == target, , drop = FALSE]
        pred$target <- target_label
        pred$model <- model
        pred$eff_range <- eff_range
        pred$rep_id <- rep_id
        predictions[[row_id]] <- pred
        row_id <- row_id + 1L
      }
    }

    completed_results <- do.call(rbind, results)
    completed_predictions <- do.call(rbind, predictions)
    saveRDS(
      list(
        settings = run_settings,
        rep_id = rep_id,
        eff_range = eff_range,
        gamma = gamma,
        phi = phi,
        sim_seed = sim_seed,
        truth = sim$truth,
        B_sf = B_sf,
        fold_id = fold_id_fixed,
        fit_spatial = fit_sp,
        fit_nonspatial = fit_ns,
        cv_spatial = cv_sp,
        cv_nonspatial = cv_ns,
        completed_results = completed_results,
        completed_predictions = completed_predictions
      ),
      checkpoint_file
    )
    write_outputs(completed_results, completed_predictions, out_dir)
    msg("Checkpointed replicate ", rep_id, "; effective range = ", eff_range)
  }
}

results <- do.call(rbind, results)
predictions <- do.call(rbind, predictions)
summaries <- write_outputs(results, predictions, out_dir)

sink(file.path(out_dir, "assessment.txt"))
cat("Repeated observed-support CV range target check\n")
cat("================================================\n\n")
cat("Settings:\n")
print(run_settings)
cat("\nMean summary:\n")
print(summaries$mean_summary)
cat("\nMean improvement:\n")
print(summaries$mean_improvement)
cat("\nInterpretation notes:\n")
cat("* Domain, plot design, and fold assignment are fixed across replicates and ranges.\n")
cat("* Replicates vary only the simulated spatial process and nugget realization.\n")
cat("* Each replicate uses the same simulation seed across ranges to couple process draws where possible.\n")
cat("* latent_eta_B scores posterior eta samples against simulated eta_B_true.\n")
cat("* observed_y_B scores nugget-added posterior samples against simulated y_B.\n")
sink()

msg("Wrote outputs to ", out_dir)
print(summaries$mean_improvement)
