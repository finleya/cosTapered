cos_cv_observed <- function(fit,
                            X_rast,
                            B_sf = NULL,
                            response_col = NULL,
                            y_B = NULL,
                            target = "observed",
                            latent_col = NULL,
                            k = 10L,
                            fold_id = NULL,
                            row_sum_tol = 1e-5,
                            fit_args = list(),
                            recover_args = list(),
                            seed = NULL,
                            keep_fits = FALSE,
                            verbose = TRUE) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(fit, "cos_fit")) {
    stop("fit must be a cos_fit object from cos_fit().")
  }
  if (!inherits(X_rast, "Raster")) {
    stop("X_rast must be a raster RasterLayer, RasterStack, or RasterBrick.")
  }
  target <- match.arg(target, choices = c("observed", "latent"), several.ok = TRUE)
  if (anyDuplicated(target)) {
    target <- unique(target)
  }

  if (is.null(B_sf)) {
    B_sf <- fit$prep$B_sf
  }
  if (!inherits(B_sf, "sf")) {
    stop("B_sf must be an sf object containing observed-support polygons.")
  }

  k <- as.integer(k)
  if (!is.finite(k) || k < 2L) {
    stop("k must be an integer greater than or equal to 2.")
  }
  if (k > nrow(B_sf)) {
    stop("k cannot exceed nrow(B_sf).")
  }

  if (!is.list(fit_args)) stop("fit_args must be a list.")
  if (!is.list(recover_args)) stop("recover_args must be a list.")

  if (is.null(y_B)) {
    if (!is.null(response_col)) {
      if (!response_col %in% names(B_sf)) {
        stop("response_col not found in B_sf.")
      }
      y_B <- B_sf[[response_col]]
    } else {
      y_B <- fit$prep$y_B
    }
  }
  y_B <- as.numeric(y_B)
  if (length(y_B) != nrow(B_sf)) {
    stop("length(y_B) must equal nrow(B_sf).")
  }
  if (any(!is.finite(y_B))) {
    stop("y_B contains non-finite values.")
  }

  eta_B <- NULL
  if ("latent" %in% target) {
    if (is.null(latent_col) || length(latent_col) != 1L || is.na(latent_col)) {
      stop("latent_col must be supplied when target includes 'latent'.")
    }
    if (!latent_col %in% names(B_sf)) {
      stop("latent_col not found in B_sf.")
    }
    eta_B <- as.numeric(B_sf[[latent_col]])
    if (length(eta_B) != nrow(B_sf)) {
      stop("length of latent_col must equal nrow(B_sf).")
    }
    if (any(!is.finite(eta_B))) {
      stop("latent_col contains non-finite values.")
    }
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(fold_id)) {
    fold_id <- sample(rep(seq_len(k), length.out = nrow(B_sf)))
  } else {
    fold_id <- as.integer(fold_id)
    if (length(fold_id) != nrow(B_sf)) {
      stop("length(fold_id) must equal nrow(B_sf).")
    }
    if (any(is.na(fold_id)) || any(fold_id < 1L) || any(fold_id > k)) {
      stop("fold_id values must be integers in 1:k.")
    }
    if (!all(seq_len(k) %in% fold_id)) {
      stop("fold_id must assign at least one observation to every fold.")
    }
  }

  # ------------------------------------------------
  # Pull fixed settings into local objects
  # ------------------------------------------------

  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE

  gamma <- fit$prep$gamma
  taper_code <- fit$prep$taper_code
  n_threads <- fit$prep$n_threads
  priors <- fit$priors

  if (spatial && (!is.finite(gamma) || gamma <= 0)) {
    stop("fit$prep$gamma must be positive and finite for spatial CV.")
  }
  if (!spatial) {
    gamma <- NA_real_
  }
  if (is.null(taper_code) || length(taper_code) != 1L) {
    taper_code <- 1L
  }
  if (is.null(n_threads) || length(n_threads) != 1L || !is.finite(n_threads)) {
    n_threads <- 1L
  }

  B_work <- B_sf
  B_work$.cos_cv_id <- seq_len(nrow(B_work))
  B_work$.cos_cv_y <- y_B
  if (!is.null(eta_B)) {
    B_work$.cos_cv_eta <- eta_B
  }

  predictions <- vector("list", k)
  fold_summary <- vector("list", k)
  fits <- if (keep_fits) vector("list", k) else NULL

  fit_args_user <- fit_args
  recover_args_user <- recover_args

  fit_defaults <- list(
    n_chains = fit$sampler$n_chains,
    tuning = fit$tuning,
    n_batch = fit$sampler$n_batch,
    batch_length = fit$sampler$batch_length,
    accept_rate = fit$sampler$accept_rate
  )
  fit_defaults <- fit_defaults[!vapply(fit_defaults, is.null, logical(1))]
  fit_args_user <- utils::modifyList(fit_defaults, fit_args_user)

  if (is.null(recover_args_user$burn_in)) {
    recover_args_user$burn_in <- 1L
  }

  # ------------------------------------------------
  # Run K-fold validation
  # ------------------------------------------------

  for (fold in seq_len(k)) {
    if (verbose) {
      message("Running observed-support CV fold ", fold, " of ", k)
    }

    test_ind <- which(fold_id == fold)
    train_ind <- which(fold_id != fold)

    B_train <- B_work[train_ind, , drop = FALSE]
    B_test <- B_work[test_ind, , drop = FALSE]

    prep_args <- list(
      X_rast = X_rast,
      B_sf = B_train,
      response_col = ".cos_cv_y",
      gamma = gamma,
      spatial = spatial,
      taper_code = taper_code,
      n_threads = n_threads,
      row_sum_tol = row_sum_tol,
      verbose = FALSE
    )
    if (!spatial) {
      prep_args$gamma <- NULL
    }

    prep <- do.call(cos_prepare, prep_args)

    fit_call_args <- c(
      list(
        prep = prep,
        priors = priors,
        verbose = FALSE
      ),
      fit_args_user
    )
    if (is.null(fit_call_args$seed) && !is.null(seed)) {
      fit_call_args$seed <- seed + 1000L * fold
    }

    fit <- do.call(cos_fit, fit_call_args)

    recover_call_args <- c(
      list(
        fit = fit,
        verbose = FALSE
      ),
      recover_args_user
    )
    if (is.null(recover_call_args$seed) && !is.null(seed)) {
      recover_call_args$seed <- seed + 2000L * fold
    }

    rec <- do.call(cos_recover_B, recover_call_args)

    holdout_blocks <- cos_make_areal_blocks(
      U_sf = B_test,
      X_rast = X_rast,
      id_col = ".cos_cv_id",
      row_sum_tol = row_sum_tol,
      verbose = FALSE
    )

    pred <- cos_predict_areal(
      fit = fit,
      rec_B = rec,
      U_blocks = holdout_blocks,
      method = "sample",
      target = "latent",
      keep_samples = TRUE,
      verbose = FALSE
    )

    eta_samples <- as.matrix(pred$eta_samples)
    theta_samples <- as.data.frame(rec$theta_samples)
    tau_sq <- theta_samples$tau_sq

    d_h <- vapply(holdout_blocks$blocks, function(z) sum(z$h^2), numeric(1))
    holdout_id <- as.integer(pred$summary$U_label)
    observed <- B_work$.cos_cv_y[match(holdout_id, B_work$.cos_cv_id)]
    latent <- if (!is.null(eta_B)) {
      B_work$.cos_cv_eta[match(holdout_id, B_work$.cos_cv_id)]
    } else {
      rep(NA_real_, length(holdout_id))
    }

    sample_sets <- list()
    truth_sets <- list()

    if ("latent" %in% target) {
      sample_sets$latent <- eta_samples
      truth_sets$latent <- latent
    }

    if ("observed" %in% target) {
      y_samples <- eta_samples
      for (jj in seq_len(ncol(y_samples))) {
        y_samples[, jj] <- y_samples[, jj] +
          stats::rnorm(nrow(y_samples), mean = 0, sd = sqrt(tau_sq * d_h[jj]))
      }
      sample_sets$observed <- y_samples
      truth_sets$observed <- observed
    }

    pred_fold_list <- lapply(names(sample_sets), function(target_i) {
      samples <- sample_sets[[target_i]]
      truth <- truth_sets[[target_i]]
      out <- data.frame(
        target = target_i,
        fold = fold,
        B_id = holdout_id,
        truth = truth,
        observed = observed,
        latent = latent,
        pred_mean = colMeans(samples),
        q025 = apply(samples, 2, stats::quantile, probs = 0.025),
        q50 = apply(samples, 2, stats::quantile, probs = 0.5),
        q975 = apply(samples, 2, stats::quantile, probs = 0.975),
        crps = vapply(seq_len(ncol(samples)), function(jj) {
          cos_crps_sample(samples[, jj], truth[jj])
        }, numeric(1)),
        stringsAsFactors = FALSE
      )
      out$interval_width <- out$q975 - out$q025
      out$error <- out$pred_mean - out$truth
      out$abs_error <- abs(out$error)
      out$covered_95 <- out$truth >= out$q025 & out$truth <= out$q975
      out
    })
    pred_fold <- do.call(rbind, pred_fold_list)

    theta_keep <- as.data.frame(fit$theta_samples)
    theta_keep <- theta_keep[
      theta_keep$iter >= as.integer(recover_args_user$burn_in),
      ,
      drop = FALSE
    ]

    fold_summary[[fold]] <- do.call(rbind, lapply(split(pred_fold, pred_fold$target), function(z) {
      data.frame(
        target = z$target[1],
        fold = fold,
        n_train = length(train_ind),
        n_test = length(test_ind),
        tau_B_sq_median = stats::median(theta_keep$tau_B_sq),
        sigma_sq_median = if (all(is.na(theta_keep$sigma_sq))) {
          NA_real_
        } else {
          stats::median(theta_keep$sigma_sq, na.rm = TRUE)
        },
        eff_range_median = if (all(is.na(theta_keep$eff_range))) {
          NA_real_
        } else {
          stats::median(theta_keep$eff_range, na.rm = TRUE)
        },
        mean_interval_width = mean(z$interval_width),
        mean_error = mean(z$error),
        MAE = mean(z$abs_error),
        bias = mean(z$error),
        cor = if (stats::sd(z$truth) > 0 && stats::sd(z$pred_mean) > 0) {
          stats::cor(z$truth, z$pred_mean)
        } else {
          NA_real_
        },
        slope = if (stats::sd(z$pred_mean) > 0) {
          unname(stats::coef(stats::lm(truth ~ pred_mean, data = z))[2])
        } else {
          NA_real_
        },
        RMSPE = sqrt(mean(z$error^2)),
        CRPS = mean(z$crps),
        coverage_95 = mean(z$covered_95),
        stringsAsFactors = FALSE
      )
    }))

    predictions[[fold]] <- pred_fold
    if (keep_fits) {
      fits[[fold]] <- fit
    }
  }

  predictions <- do.call(rbind, predictions)
  rownames(predictions) <- NULL

  fold_summary <- do.call(rbind, fold_summary)
  rownames(fold_summary) <- NULL

  # ------------------------------------------------
  # Overall summaries
  # ------------------------------------------------

  overall <- do.call(rbind, lapply(split(predictions, predictions$target), function(z) {
    data.frame(
      target = z$target[1],
      n = nrow(z),
      k = k,
      RMSPE = sqrt(mean(z$error^2)),
      CRPS = mean(z$crps),
      coverage_95 = mean(z$covered_95),
      mean_interval_width = mean(z$interval_width),
      mean_error = mean(z$error),
      MAE = mean(z$abs_error),
      bias = mean(z$error),
      cor = if (stats::sd(z$truth) > 0 && stats::sd(z$pred_mean) > 0) {
        stats::cor(z$truth, z$pred_mean)
      } else {
        NA_real_
      },
      slope = if (stats::sd(z$pred_mean) > 0) {
        unname(stats::coef(stats::lm(truth ~ pred_mean, data = z))[2])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall) <- NULL

  out <- list(
    summary = overall,
    fold_summary = fold_summary,
    predictions = predictions,
    fold_id = fold_id,
    spatial = spatial,
    target = target,
    k = k,
    fit_args = fit_args_user,
    recover_args = recover_args_user,
    fits = fits,
    call = match.call()
  )

  class(out) <- "cos_cv_observed"
  out
}

cos_crps_sample <- function(draws, y) {
  draws <- sort(as.numeric(draws))
  draws <- draws[is.finite(draws)]
  y <- as.numeric(y)

  if (length(y) != 1L || !is.finite(y)) {
    return(NA_real_)
  }
  n <- length(draws)
  if (n == 0L) {
    return(NA_real_)
  }

  mean(abs(draws - y)) -
    sum((2 * seq_len(n) - n - 1) * draws) / (n^2)
}
