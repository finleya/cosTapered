print.cos_prep <- function(x, ...) {
  spatial <- isTRUE(x$spatial)
  if (is.null(x$spatial)) spatial <- TRUE

  if (spatial) {
    cat("Tapered COS preparation object\n")
  } else {
    cat("Non-spatial COS preparation object\n")
  }
  cat("  observed supports B: ", nrow(x$B_sf), "\n", sep = "")
  cat("  fine cells A:        ", nrow(x$A_df), "\n", sep = "")
  cat("  covariates:          ", paste(colnames(x$X), collapse = ", "), "\n", sep = "")
  if (spatial) {
    cat("  gamma:               ", x$gamma, "\n", sep = "")
  }
  cat("Use summary() for support-weight diagnostics.\n")

  invisible(x)
}

print.cos_priors <- function(x, ...) {
  p <- length(x$beta$mu)
  spatial <- is.finite(x$phi$lower) && is.finite(x$phi$upper)

  if (spatial) {
    cat("Tapered COS prior object\n")
  } else {
    cat("Non-spatial COS prior object\n")
  }
  cat("  beta coefficients: ", p, "\n", sep = "")
  cat("  tau_B_sq IG:       shape = ", x$tau_B_sq$shape,
      ", scale = ", x$tau_B_sq$scale, "\n", sep = "")
  if (spatial) {
    cat("  sigma_sq IG:       shape = ", x$sigma_sq$shape,
        ", scale = ", x$sigma_sq$scale, "\n", sep = "")
    cat("  phi bounds:        [", x$phi$lower, ", ", x$phi$upper, "]\n", sep = "")
  }

  invisible(x)
}

print.cos_fit <- function(x, ...) {
  spatial <- isTRUE(x$spatial)
  if (is.null(x$spatial)) spatial <- TRUE

  if (spatial) {
    cat("Tapered COS fit\n")
  } else {
    cat("Non-spatial COS fit\n")
  }
  cat("  chains:       ", x$sampler$n_chains, "\n", sep = "")
  cat("  batches:      ", x$sampler$n_batch, "\n", sep = "")
  cat("  batch length: ", x$sampler$batch_length, "\n", sep = "")
  cat("  saved draws:  ", nrow(x$theta_samples), "\n", sep = "")
  cat("Use summary() for posterior summaries and plot() for trace plots.\n")

  invisible(x)
}

print.cos_recovery_B <- function(x, ...) {
  spatial <- isTRUE(x$spatial)
  if (is.null(x$spatial)) spatial <- TRUE

  if (spatial) {
    cat("Tapered COS observed-support recovery\n")
  } else {
    cat("Non-spatial COS observed-support recovery\n")
  }
  cat("  posterior draws: ", nrow(x$beta_samples), "\n", sep = "")
  cat("  observed units:  ", ncol(x$eta_B_samples), "\n", sep = "")
  cat("  beta terms:      ", ncol(x$beta_samples), "\n", sep = "")
  cat("Use summary() for beta, omega_B, and eta_B summaries.\n")

  invisible(x)
}

print.cos_prediction_fine <- function(x, ...) {
  spatial <- isTRUE(x$spatial)
  if (is.null(x$spatial)) spatial <- TRUE

  if (spatial) {
    cat("Tapered COS fine prediction\n")
  } else {
    cat("Non-spatial COS fine prediction\n")
  }
  cat("  method:          ", x$method, "\n", sep = "")
  cat("  target:          ", x$target, "\n", sep = "")
  cat("  prediction cells:", x$n_pred, "\n", sep = "")
  cat("  posterior draws: ", x$n_samples, "\n", sep = "")
  cat("  stored samples:  ", ifelse(x$keep_samples, "yes", "no"), "\n", sep = "")
  cat("Use summary() for prediction summaries.\n")

  invisible(x)
}

print.cos_areal_blocks <- function(x, ...) {
  cat("Tapered COS areal prediction blocks\n")
  cat("  prediction supports U: ", length(x$blocks), "\n", sep = "")
  cat("  covariates:            ", paste(x$x_names, collapse = ", "), "\n", sep = "")
  cat("Use summary() for support-weight diagnostics.\n")

  invisible(x)
}

print.cos_prediction_areal <- function(x, ...) {
  spatial <- isTRUE(x$spatial)
  if (is.null(x$spatial)) spatial <- TRUE

  if (spatial) {
    cat("Tapered COS areal prediction\n")
  } else {
    cat("Non-spatial COS areal prediction\n")
  }
  cat("  method:          ", x$method, "\n", sep = "")
  cat("  target:          ", x$target, "\n", sep = "")
  cat("  prediction units:", x$n_U, "\n", sep = "")
  cat("  posterior draws: ", x$n_samples, "\n", sep = "")
  cat("  stored samples:  ", ifelse(x$keep_samples, "yes", "no"), "\n", sep = "")
  cat("Use summary() for prediction summaries.\n")

  invisible(x)
}

print.cos_cv_observed <- function(x, ...) {
  spatial <- isTRUE(x$spatial)

  if (spatial) {
    cat("Tapered COS observed-support K-fold validation\n")
  } else {
    cat("Non-spatial COS observed-support K-fold validation\n")
  }
  cat("  folds:          ", x$k, "\n", sep = "")
  cat("  targets:        ", paste(x$target, collapse = ", "), "\n", sep = "")
  cat("  held-out units: ", length(unique(x$predictions$B_id)), "\n", sep = "")
  if (nrow(x$summary) == 1L) {
    cat("  RMSPE:          ", signif(x$summary$RMSPE, 4), "\n", sep = "")
    cat("  CRPS:           ", signif(x$summary$CRPS, 4), "\n", sep = "")
    cat("  coverage_95:    ", signif(x$summary$coverage_95, 4), "\n", sep = "")
  } else {
    cat("  scores:         use summary() for target-specific scores\n", sep = "")
  }
  cat("Use summary() for fold summaries and held-out predictions.\n")

  invisible(x)
}

print.cos_simulation <- function(x, ...) {
  cat("COS simulated dataset\n")
  cat("  observed supports B: ", nrow(x$B_sf), "\n", sep = "")
  cat("  fine cells A:        ", x$truth$n_A, "\n", sep = "")
  cat("  gamma:               ", x$truth$gamma, "\n", sep = "")
  cat("  effective range:     ", x$truth$eff_range_true, "\n", sep = "")
  cat("Known parameters are available in $truth.\n")

  invisible(x)
}
