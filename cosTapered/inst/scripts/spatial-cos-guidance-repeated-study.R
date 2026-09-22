# Reproduce the spatial-guidance article with the installed package.
# source(system.file("scripts", "spatial-cos-guidance-repeated-study.R",
#                    package = "cosTapered"))
# Settings are recorded with every run. Checkpoints can resume the same settings.

study_env_number <- function(name, default, integer = FALSE) {
  value <- Sys.getenv(name, unset = "")
  out <- if (value == "") default else as.numeric(value)
  if (length(out) != 1L || !is.finite(out) || out <= 0 ||
      (integer && out != floor(out))) stop("Invalid ", name)
  out
}
study_env_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  out <- if (value == "") default else as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (!length(out) || any(!is.finite(out)) || anyDuplicated(out)) stop("Invalid ", name)
  out
}
study_message <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n", sep = "")
}
study_compatible_settings <- function(saved, current) {
  # Increasing the retry budget can reuse successful cases; their actual chain
  # lengths remain recorded in diagnostics.csv and attempts.csv.
  budget_ok <- current$max_attempts >= saved$max_attempts
  saved$max_attempts <- current$max_attempts <- NULL
  budget_ok && identical(saved, current)
}

# These priors depend only on declared study settings, never on responses or
# generating parameter values. The same beta and nugget priors enter both fits.
study_priors <- function(prep, prior_scale) {
  cosTapered::cos_default_priors(
    prep, beta_mu = c(intercept = 0, chm = 0),
    beta_sd = c(intercept = 100, chm = 10),
    tau_B_shape = 2, tau_B_scale = prior_scale,
    sigma_shape = 2, sigma_scale = prior_scale,
    phi_lower = 3 / 1500, phi_upper = 3 / 10
  )
}
study_starts <- function(spatial, prior_scale, n_chains) {
  factors <- exp(seq(log(0.25), log(4), length.out = n_chains))
  out <- data.frame(tau_B_sq = prior_scale * factors)
  if (spatial) {
    out$sigma_sq <- prior_scale * rev(factors)
    out$phi <- 3 / exp(seq(log(25), log(1200), length.out = n_chains))
  }
  out
}
study_diagnostics <- function(fits, burn_in, settings) {
  rows <- lapply(seq_along(fits), function(fold) {
    z <- fits[[fold]]$theta_samples
    z <- z[z$iter >= burn_in, ]
    variables <- if (isTRUE(fits[[fold]]$spatial)) {
      c("tau_B_sq", "sigma_sq", "eff_range", "lp")
    } else c("tau_B_sq", "lp")
    do.call(rbind, lapply(variables, function(variable) {
      chains <- split(z[[variable]], z$chain)
      stopifnot(length(chains) == settings$n_chains,
                length(unique(lengths(chains))) == 1L)
      draws <- do.call(cbind, chains)
      rh <- posterior::rhat(draws)
      eb <- posterior::ess_bulk(draws)
      et <- posterior::ess_tail(draws)
      data.frame(fold = fold, variable = variable, rhat = rh,
                 ess_bulk = eb, ess_tail = et,
                 passed = is.finite(rh) && is.finite(eb) && is.finite(et) &&
                   rh < settings$rhat_limit && eb >= settings$ess_min &&
                   et >= settings$ess_min)
    }))
  })
  do.call(rbind, rows)
}

study_case <- function(task, X_rast, forest_sf, settings, out_dir) {
  tag <- sprintf("rep_%03d_range_%04d_prior_%04d", task$rep_id,
                 task$eff_range, task$prior_scale)
  checkpoint <- file.path(out_dir, "checkpoints", paste0(tag, ".rds"))
  previous <- NULL
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (!study_compatible_settings(saved$settings, settings)) {
      stop("Checkpoint settings differ; choose a new CTV_OUT_DIR: ", checkpoint)
    }
    if (all(saved$diagnostics$passed)) {
      study_message("Resuming ", tag)
      return(saved)
    }
    previous <- saved
  }
  study_message("Starting ", tag)
  sim <- cosTapered::cos_simulate_data(
    X_rast = X_rast, forest_sf = forest_sf, n_B = 49L, plot_radius = 16,
    gamma = settings$gamma, beta = c(intercept = -10, chm = 15.87132),
    tau_B_sq = 250, sigma_sq = 750,
    phi = 3 / if (task$eff_range == 0) 150 else task$eff_range,
    design_seed = 1101L, sim_seed = 220200L + task$rep_id,
    n_threads = settings$n_threads, verbose = FALSE
  )
  B_sf <- sim$B_sf
  # The simulator requires sigma_sq > 0. Removing its known spatial component
  # produces the exact nonspatial truth case, retaining the same Gaussian nugget.
  if (task$eff_range == 0) {
    B_sf$eta_B_true <- as.numeric(sim$prep$X_B %*% c(-10, 15.87132))
    B_sf$omega_B_true <- 0
    B_sf$y_B <- B_sf$eta_B_true + B_sf$eps_B_true
  }
  set.seed(3303L)
  fold_id <- sample(rep(seq_len(settings$k_folds), length.out = nrow(B_sf)))
  scores <- predictions <- diagnostics <- attempts <- list()
  for (model in c("spatial", "nonspatial")) {
    first_attempt <- 1L
    if (!is.null(previous)) {
      old_diag <- previous$diagnostics[previous$diagnostics$model == model, ]
      old_attempts <- previous$attempts[previous$attempts$model == model, ]
      attempts[[length(attempts) + 1L]] <- old_attempts[, setdiff(names(old_attempts), names(task))]
      if (all(old_diag$passed)) {
        diagnostics[[model]] <- old_diag[, setdiff(names(old_diag), names(task))]
        scores[[model]] <- previous$scores[previous$scores$model == model, setdiff(names(previous$scores), names(task))]
        predictions[[model]] <- previous$predictions[previous$predictions$model == model, setdiff(names(previous$predictions), names(task))]
        next
      }
      first_attempt <- max(old_attempts$attempt) + 1L
      if (first_attempt > settings$max_attempts) {
        stop("Increase CTV_MAX_ATTEMPTS to retry failed fits in ", tag)
      }
    }
    spatial <- model == "spatial"
    prep <- cosTapered::cos_prepare(
      X_rast, B_sf, response_col = "y_B", spatial = spatial,
      gamma = settings$gamma, n_threads = settings$n_threads, verbose = FALSE
    )
    priors <- study_priors(prep, task$prior_scale)
    starts <- study_starts(spatial, task$prior_scale, settings$n_chains)
    tuning <- if (spatial) c(log_tau_B_sq = 0.5, log_sigma_sq = 0.5, z_phi = 0.5) else c(log_tau_B_sq = 0.5)
    # This one-draw template only supplies geometry and the fixed priors to CV.
    # Every scored prediction below comes from a fresh training-fold refit.
    template <- cosTapered::cos_fit(
      prep, priors, n_chains = 1L, starting = starts[1, , drop = FALSE],
      tuning = tuning, n_batch = 1L, batch_length = 1L, seed = 1L, verbose = FALSE
    )
    seed <- 300000L + 10000L * task$rep_id + 10L * task$eff_range + as.integer(spatial)
    for (attempt in seq.int(first_attempt, settings$max_attempts)) {
      n_batch <- settings$n_batch * 2L^(attempt - 1L)
      iterations <- n_batch * settings$batch_length
      burn_in <- floor(iterations / 2) + 1L
      thin <- max(1L, ceiling((iterations - burn_in + 1L) / settings$draws_per_chain))
      study_message(tag, " ", model, ": ", iterations, " iterations per chain")
      cv <- cosTapered::cos_cv_observed(
        template, X_rast, B_sf, response_col = "y_B",
        target = c("latent", "observed"), latent_col = "eta_B_true",
        k = settings$k_folds, fold_id = fold_id,
        fit_args = list(n_chains = settings$n_chains, starting = starts,
                        tuning = tuning, n_batch = n_batch,
                        batch_length = settings$batch_length, report = n_batch + 1L),
        recover_args = list(burn_in = burn_in, thin = thin),
        seed = seed, keep_fits = TRUE, verbose = FALSE
      )
      diag <- study_diagnostics(cv$fits, burn_in, settings)
      passed <- all(diag$passed)
      attempts[[length(attempts) + 1L]] <- data.frame(
        model = model, attempt = attempt, iterations = iterations,
        max_rhat = max(diag$rhat), min_ess_bulk = min(diag$ess_bulk),
        min_ess_tail = min(diag$ess_tail), passed = passed,
        elapsed_seconds = cv$timing[["elapsed"]]
      )
      if (passed) break
      study_message(tag, " ", model, " attempt ", attempt,
                    ": Rhat=", round(max(diag$rhat), 3),
                    ", min ESS=", round(min(diag$ess_bulk, diag$ess_tail)))
    }
    diag$model <- model
    diag$iterations <- iterations
    diag$burn_in <- burn_in
    diag$thin <- thin
    diag$prediction_draws <- settings$n_chains * length(seq.int(burn_in, iterations, by = thin))
    diagnostics[[model]] <- diag
    z <- cv$summary
    z$model <- model
    z$diagnostics_passed <- passed
    scores[[model]] <- z
    z <- cv$predictions
    z$model <- model
    predictions[[model]] <- z
    if (!passed) {
      saveRDS(lapply(cv$fits, function(fit) fit$theta_samples),
              file.path(out_dir, "checkpoints", paste0(tag, "_", model, "_failed_traces.rds")))
    }
  }
  add_task <- function(z) {
    z <- do.call(rbind, z)
    for (name in names(task)) z[[name]] <- task[[name]]
    rownames(z) <- NULL
    z
  }
  result <- list(settings = settings, scores = add_task(scores),
                 predictions = add_task(predictions), diagnostics = add_task(diagnostics),
                 attempts = add_task(attempts),
                 design = data.frame(B_id = B_sf$B_id,
                                     sf::st_coordinates(sf::st_centroid(sf::st_geometry(B_sf))),
                                     area_m2 = B_sf$area_m2, fold = fold_id))
  saveRDS(result, paste0(checkpoint, ".tmp"))
  stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
  study_message("Completed ", tag)
  result
}

study_summarize <- function(data, groups, metrics) {
  parts <- split(data, interaction(data[groups], drop = TRUE))
  out <- lapply(parts, function(z) {
    row <- z[1, groups, drop = FALSE]
    row$n_reps <- nrow(z)
    for (metric in metrics) {
      x <- z[[metric]]
      stopifnot(all(is.finite(x)))
      row[[paste0(metric, "_mean")]] <- mean(x)
      row[[paste0(metric, "_se")]] <- stats::sd(x) / sqrt(length(x))
      row[[paste0(metric, "_sd")]] <- stats::sd(x)
    }
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[do.call(order, out[groups]), ]
}
study_improvement <- function(scores) {
  keys <- c("rep_id", "eff_range", "prior_scale", "target")
  z <- merge(scores[scores$model == "spatial", ],
             scores[scores$model == "nonspatial", ], by = keys,
             suffixes = c("_spatial", "_nonspatial"))
  stopifnot(nrow(z) * 2L == nrow(scores))
  for (metric in c("RMSPE", "CRPS")) {
    z[[paste0(metric, "_improvement_pct")]] <-
      100 * (z[[paste0(metric, "_nonspatial")]] - z[[paste0(metric, "_spatial")]]) /
      z[[paste0(metric, "_nonspatial")]]
    z[[paste0(metric, "_win")]] <- as.numeric(z[[paste0(metric, "_spatial")]] < z[[paste0(metric, "_nonspatial")]])
  }
  z
}

run_spatial_guidance_study <- function() {
  for (pkg in c("cosTapered", "sf", "raster", "posterior")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package ", pkg, " first.")
  }
  settings <- list(
    study_version = 2L,
    n_reps = study_env_number("CTV_N_REPS", 10L, TRUE),
    eff_ranges = study_env_vector("CTV_EFF_RANGES", c(0, 50, 150, 350, 750, 1000)),
    prior_scales = study_env_vector("CTV_PRIOR_SCALES", c(500, 250, 1000)),
    aggregate_factor = study_env_number("CTV_AGG_FACTOR", 4L, TRUE),
    gamma = study_env_number("CTV_GAMMA", 1500),
    k_folds = study_env_number("CTV_K_FOLDS", 5L, TRUE),
    n_chains = study_env_number("CTV_N_CHAINS", 4L, TRUE),
    n_batch = study_env_number("CTV_N_BATCH", 400L, TRUE),
    batch_length = study_env_number("CTV_BATCH_LENGTH", 25L, TRUE),
    max_attempts = study_env_number("CTV_MAX_ATTEMPTS", 3L, TRUE),
    draws_per_chain = study_env_number("CTV_DRAWS_PER_CHAIN", 500L, TRUE),
    rhat_limit = 1.01, ess_min = 400,
    n_threads = study_env_number("CTV_N_THREADS", 1L, TRUE)
  )
  if (any(settings$eff_ranges < 0) || any(settings$prior_scales <= 0) ||
      settings$n_chains < 4 || settings$k_folds < 2 ||
      settings$n_batch * settings$batch_length < 1000) stop("Invalid study settings.")
  workers <- study_env_number("CTV_WORKERS", 1L, TRUE)
  if (.Platform$OS.type != "unix" && workers > 1L) stop("Use CTV_WORKERS=1 on Windows.")
  out_dir <- Sys.getenv("CTV_OUT_DIR", "spatial-cos-guidance-output")
  dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
  existing <- file.path(out_dir, "settings.rds")
  if (file.exists(existing) && !study_compatible_settings(readRDS(existing), settings)) {
    stop("Settings differ from existing run. Choose a new CTV_OUT_DIR.")
  }
  saveRDS(settings, existing)
  dput(settings, file = file.path(out_dir, "settings.R"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session-info.txt"))
  inputs <- system.file("extdata", c("example_chm.tif", "example_forest.gpkg"), package = "cosTapered")
  write.csv(data.frame(file = basename(inputs), md5 = unname(tools::md5sum(inputs))),
            file.path(out_dir, "input-checksums.csv"), row.names = FALSE)
  chm <- raster::raster(inputs[1])
  if (settings$aggregate_factor > 1) {
    chm <- raster::aggregate(chm, fact = settings$aggregate_factor, fun = mean, na.rm = TRUE)
  }
  names(chm) <- "chm"
  intercept <- chm
  intercept[] <- ifelse(is.na(raster::getValues(chm)), NA_real_, 1)
  X_rast <- raster::stack(intercept, chm)
  names(X_rast) <- c("intercept", "chm")
  forest_sf <- sf::st_read(inputs[2], quiet = TRUE)
  forest_sf <- sf::st_transform(forest_sf, sf::st_crs(raster::crs(chm)))
  tasks <- expand.grid(rep_id = seq_len(settings$n_reps),
                       eff_range = settings$eff_ranges, prior_scale = settings$prior_scales)
  jobs <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))
  started <- Sys.time()
  results <- parallel::mclapply(jobs, study_case, X_rast = X_rast,
                               forest_sf = forest_sf, settings = settings,
                               out_dir = out_dir, mc.cores = workers,
                               mc.preschedule = FALSE, mc.set.seed = FALSE)
  failed <- vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) stop("Study task failed; check the log. Completed checkpoints are retained.")
  tables <- list()
  for (name in c("scores", "predictions", "diagnostics", "attempts")) {
    tables[[name]] <- do.call(rbind, lapply(results, `[[`, name))
    write.csv(tables[[name]], file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  stopifnot(all(vapply(results, function(x) identical(x$design, results[[1]]$design), logical(1))))
  write.csv(results[[1]]$design, file.path(out_dir, "design.csv"), row.names = FALSE)
  improvements <- study_improvement(tables$scores)
  write.csv(improvements, file.path(out_dir, "paired-scores.csv"), row.names = FALSE)
  summary <- study_summarize(tables$scores, c("prior_scale", "target", "model", "eff_range"),
                             c("RMSPE", "CRPS", "coverage_95", "mean_interval_width"))
  gains <- study_summarize(improvements, c("prior_scale", "target", "eff_range"),
                           c("RMSPE_improvement_pct", "CRPS_improvement_pct", "RMSPE_win", "CRPS_win"))
  write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(gains, file.path(out_dir, "improvement.csv"), row.names = FALSE)
  writeLines(sprintf("Elapsed wall time: %.1f minutes", as.numeric(difftime(Sys.time(), started, units = "mins"))),
              file.path(out_dir, "runtime.txt"))
  if (!all(tables$diagnostics$passed)) {
    stop("Some fits failed the diagnostic thresholds. Scores are saved, but do not publish them as validated results.")
  }
  study_message("All fits passed diagnostics. Results saved to ", out_dir)
  invisible(tables)
}

if (!identical(Sys.getenv("CTV_DEFINE_ONLY"), "TRUE")) run_spatial_guidance_study()
