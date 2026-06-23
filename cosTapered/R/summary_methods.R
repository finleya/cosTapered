summary.cos_fit <- function(object,
                            burn_in = 1,
                            pars = c("tau_B_sq", "tau_sq", "sigma_sq", "phi", "eff_range", "lp"),
                            probs = c(0.5, 0.025, 0.975),
                            digits = 3,
                            ...) {
  # ------------------------------------------------
  # Check inputs and select retained samples
  # ------------------------------------------------

  fit <- object
  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE
  family <- fit$family
  if (is.null(family)) family <- "gaussian"

  burn_in <- as.integer(burn_in)
  if (!is.finite(burn_in) || burn_in < 1L) {
    stop("burn_in must be a positive integer.")
  }

  theta <- as.data.frame(fit$theta_samples)
  theta <- theta[theta$iter >= burn_in, , drop = FALSE]
  if (nrow(theta) == 0L) stop("No theta samples remain after applying burn_in.")

  pars <- intersect(pars, names(theta))
  pars <- pars[vapply(pars, function(nm) any(is.finite(theta[[nm]])), logical(1))]

  summary_source <- theta
  summary_label <- "Parameter posterior summaries"
  if (length(pars) == 0L && family %in% c("binomial", "negative_binomial") && !is.null(fit$beta_samples)) {
    beta <- as.data.frame(fit$beta_samples)
    beta <- beta[beta$iter >= burn_in, , drop = FALSE]
    beta_names <- setdiff(names(beta), c("chain", "iter"))
    beta_names <- beta_names[vapply(beta_names, function(nm) any(is.finite(beta[[nm]])), logical(1))]
    if (length(beta_names) > 0L) {
      summary_source <- beta
      pars <- beta_names
      summary_label <- "Beta posterior summaries"
    }
  }
  if (length(pars) == 0L) stop("No requested parameters have finite retained samples.")

  # ------------------------------------------------
  # Parameter summaries
  # ------------------------------------------------

  theta_summary <- t(vapply(pars, function(nm) {
    stats::quantile(summary_source[[nm]], probs = probs, na.rm = TRUE)
  }, numeric(length(probs))))

  colnames(theta_summary) <- paste0(100 * probs, "%")
  theta_summary <- round(theta_summary, digits)

  # ------------------------------------------------
  # Chain and sampler summaries
  # ------------------------------------------------

  chain_summary <- do.call(rbind, lapply(fit$chain_fits, function(z) {
    c(
      mean_accept = mean(z$batch.accept.rate, na.rm = TRUE),
      final_accept = tail(z$batch.accept.rate, 1),
      final_sigma_sq_m = z$sigma_sq_m
    )
  }))
  chain_summary <- as.data.frame(chain_summary)
  chain_summary$chain <- seq_len(nrow(chain_summary))
  chain_summary <- chain_summary[, c("chain", "mean_accept", "final_accept", "final_sigma_sq_m")]

  out <- list(
    call = fit$call,
    burn_in = burn_in,
    n_chains = fit$sampler$n_chains,
    n_batch = fit$sampler$n_batch,
    batch_length = fit$sampler$batch_length,
    n_retained = nrow(theta),
    theta_summary = theta_summary,
    summary_label = summary_label,
    chain_summary = chain_summary,
    spatial = spatial,
    family = family,
    priors = fit$priors
  )

  class(out) <- "summary.cos_fit"
  out
}

print.summary.cos_fit <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS fit\n")
  } else {
    cat("Non-spatial COS fit\n")
  }
  cat("  family:        ", x$family, "\n", sep = "")
  cat("  chains:        ", x$n_chains, "\n", sep = "")
  cat("  batches:       ", x$n_batch, "\n", sep = "")
  cat("  batch length:  ", x$batch_length, "\n", sep = "")
  cat("  burn-in start: ", x$burn_in, "\n", sep = "")
  cat("  retained draws:", x$n_retained, "\n\n", sep = "")

  cat(x$summary_label, ":\n", sep = "")
  print(x$theta_summary)

  cat("\nChain summaries:\n")
  print(round(x$chain_summary, 4), row.names = FALSE)

  invisible(x)
}

summary.cos_recovery_B <- function(object,
                                   probs = c(0.5, 0.025, 0.975),
                                   digits = 3,
                                   include_eta_B = TRUE,
                                   ...) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  rec <- object
  spatial <- isTRUE(rec$spatial)
  if (is.null(rec$spatial)) spatial <- TRUE

  beta_samples <- as.matrix(rec$beta_samples)
  omega_B_samples <- as.matrix(rec$omega_B_samples)
  eta_B_samples <- as.matrix(rec$eta_B_samples)

  if (nrow(beta_samples) == 0L) stop("No beta samples found in recovery object.")

  # ------------------------------------------------
  # Parameter summaries
  # ------------------------------------------------

  beta_summary <- t(apply(beta_samples, 2, function(z) {
    stats::quantile(z, probs = probs, na.rm = TRUE)
  }))
  colnames(beta_summary) <- paste0(100 * probs, "%")
  beta_summary <- round(beta_summary, digits)

  omega_B_summary <- data.frame(
    B_id = seq_len(ncol(omega_B_samples)),
    mean = colMeans(omega_B_samples),
    sd = apply(omega_B_samples, 2, stats::sd),
    q025 = apply(omega_B_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    q50 = apply(omega_B_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
    q975 = apply(omega_B_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
  omega_B_summary[, -1] <- round(omega_B_summary[, -1, drop = FALSE], digits)

  eta_B_summary <- NULL
  if (include_eta_B) {
    eta_B_summary <- data.frame(
      B_id = seq_len(ncol(eta_B_samples)),
      mean = colMeans(eta_B_samples),
      sd = apply(eta_B_samples, 2, stats::sd),
      q025 = apply(eta_B_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
      q50 = apply(eta_B_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
      q975 = apply(eta_B_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    )
    eta_B_summary[, -1] <- round(eta_B_summary[, -1, drop = FALSE], digits)
  }

  out <- list(
    call = rec$call,
    n_samples = nrow(beta_samples),
    n_B = ncol(omega_B_samples),
    spatial = spatial,
    beta_summary = beta_summary,
    omega_B_summary = omega_B_summary,
    eta_B_summary = eta_B_summary
  )

  class(out) <- "summary.cos_recovery_B"
  out
}

print.summary.cos_recovery_B <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS observed-support recovery\n")
  } else {
    cat("Non-spatial COS observed-support recovery\n")
  }
  cat("  retained draws: ", x$n_samples, "\n", sep = "")
  cat("  observed units: ", x$n_B, "\n\n", sep = "")

  cat("Beta posterior summaries:\n")
  print(x$beta_summary)

  if (x$spatial) {
    cat("\nOmega_B summaries are available in $omega_B_summary.")
  } else {
    cat("\nOmega_B is fixed at zero; summaries are available in $omega_B_summary.")
  }
  if (!is.null(x$eta_B_summary)) {
    cat("\nEta_B summaries are available in $eta_B_summary.")
  }
  cat("\n")

  invisible(x)
}

summary.cos_prep <- function(object, digits = 3, ...) {
  prep <- object
  spatial <- isTRUE(prep$spatial)
  if (is.null(prep$spatial)) spatial <- TRUE

  row_sums <- prep$row_sums
  row_sum_summary <- summary(row_sums)
  row_sum_summary <- round(row_sum_summary, digits)

  out <- list(
    call = prep$call,
    n_B = nrow(prep$B_sf),
    n_A = nrow(prep$A_df),
    p = ncol(prep$X),
    x_names = colnames(prep$X),
    spatial = spatial,
    gamma = prep$gamma,
    taper_code = prep$taper_code,
    missing = prep$missing,
    row_sum_summary = row_sum_summary,
    row_sums = row_sums
  )

  class(out) <- "summary.cos_prep"
  out
}

print.summary.cos_prep <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS preparation object\n")
  } else {
    cat("Non-spatial COS preparation object\n")
  }
  cat("  observed supports B: ", x$n_B, "\n", sep = "")
  cat("  fine cells A:        ", x$n_A, "\n", sep = "")
  cat("  covariates:          ", paste(x$x_names, collapse = ", "), "\n", sep = "")
  if (x$spatial) {
    cat("  gamma:               ", x$gamma, "\n", sep = "")
    cat("  taper_code:          ", x$taper_code, "\n", sep = "")
  }
  if (!is.null(x$missing)) {
    cat("  missing:             ", x$missing, "\n", sep = "")
  }
  cat("\n")

  cat("H_BA row-sum summary:\n")
  print(x$row_sum_summary)

  invisible(x)
}

summary.cos_prediction_fine <- function(object, digits = 3, ...) {
  pred <- object
  spatial <- isTRUE(pred$spatial)
  if (is.null(pred$spatial)) spatial <- TRUE
  s <- pred$summary

  out <- list(
    call = pred$call,
    spatial_uncertainty = pred$spatial_uncertainty,
    uncertainty = pred$uncertainty,
    target = pred$target,
    n_pred = pred$n_pred,
    n_samples = pred$n_samples,
    keep_samples = pred$keep_samples,
    spatial = spatial,
    eta_mean_summary = round(summary(s$eta_mean), digits),
    eta_sd_summary = round(summary(s$eta_sd), digits),
    omega_mean_summary = round(summary(s$omega_mean), digits),
    omega_sd_summary = round(summary(s$omega_sd), digits),
    prediction_summary = s
  )

  if (!is.null(s$y_mean)) {
    out$y_mean_summary <- round(summary(s$y_mean), digits)
    out$y_sd_summary <- round(summary(s$y_sd), digits)
  }

  class(out) <- "summary.cos_prediction_fine"
  out
}

print.summary.cos_prediction_fine <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS fine prediction/recovery\n")
  } else {
    cat("Non-spatial COS fine prediction/recovery\n")
  }
  cat("  spatial uncertainty: ", x$spatial_uncertainty, "\n", sep = "")
  if (!is.null(x$uncertainty$note)) {
    cat("  uncertainty:     ", x$uncertainty$note, "\n", sep = "")
  }
  cat("  target:          ", x$target, "\n", sep = "")
  cat("  prediction cells:", x$n_pred, "\n", sep = "")
  cat("  posterior draws: ", x$n_samples, "\n", sep = "")
  cat("  stored samples:  ", ifelse(x$keep_samples, "yes", "no"), "\n\n", sep = "")

  cat("Eta mean summary:\n")
  print(x$eta_mean_summary)
  cat("\nEta SD summary:\n")
  print(x$eta_sd_summary)

  cat("\nOmega mean summary:\n")
  print(x$omega_mean_summary)
  cat("\nOmega SD summary:\n")
  print(x$omega_sd_summary)

  if (!is.null(x$y_mean_summary)) {
    cat("\nY mean summary:\n")
    print(x$y_mean_summary)
    cat("\nY SD summary:\n")
    print(x$y_sd_summary)
  }

  cat("\nPrediction summaries are available in $prediction_summary.\n")

  invisible(x)
}

summary.cos_areal_blocks <- function(object, digits = 3, ...) {
  U_blocks <- object

  row_sums <- U_blocks$row_sums
  row_sum_summary <- round(summary(row_sums), digits)

  n_cell <- vapply(U_blocks$blocks, `[[`, numeric(1), "n_cell")
  n_cell_summary <- round(summary(n_cell), digits)

  out <- list(
    call = U_blocks$call,
    n_U = length(U_blocks$blocks),
    x_names = U_blocks$x_names,
    missing = U_blocks$missing,
    row_sum_summary = row_sum_summary,
    n_cell_summary = n_cell_summary,
    row_sums = row_sums
  )

  class(out) <- "summary.cos_areal_blocks"
  out
}

print.summary.cos_areal_blocks <- function(x, ...) {
  cat("Tapered COS areal prediction blocks\n")
  cat("  prediction supports U: ", x$n_U, "\n", sep = "")
  cat("  covariates:            ", paste(x$x_names, collapse = ", "), "\n\n", sep = "")
  if (!is.null(x$missing)) {
    cat("  missing:              ", x$missing, "\n\n", sep = "")
  }

  cat("H_UA row-sum summary:\n")
  print(x$row_sum_summary)

  cat("\nFine cells per areal unit:\n")
  print(x$n_cell_summary)

  invisible(x)
}

summary.cos_prediction_areal <- function(object, digits = 3, ...) {
  pred <- object
  spatial <- isTRUE(pred$spatial)
  if (is.null(pred$spatial)) spatial <- TRUE

  eta_mean_summary <- round(summary(pred$summary$eta_mean), digits)
  eta_sd_summary <- round(summary(pred$summary$eta_sd), digits)
  omega_mean_summary <- round(summary(pred$summary$omega_mean), digits)
  omega_sd_summary <- round(summary(pred$summary$omega_sd), digits)
  row_sum_summary <- round(summary(pred$summary$row_sum), digits)

  out <- list(
    call = pred$call,
    spatial_uncertainty = pred$spatial_uncertainty,
    uncertainty = pred$uncertainty,
    target = pred$target,
    n_U = pred$n_U,
    n_samples = pred$n_samples,
    keep_samples = pred$keep_samples,
    spatial = spatial,
    row_sum_summary = row_sum_summary,
    eta_mean_summary = eta_mean_summary,
    eta_sd_summary = eta_sd_summary,
    omega_mean_summary = omega_mean_summary,
    omega_sd_summary = omega_sd_summary,
    prediction_summary = pred$summary
  )

  if (!is.null(pred$summary$y_mean)) {
    out$y_mean_summary <- round(summary(pred$summary$y_mean), digits)
    out$y_sd_summary <- round(summary(pred$summary$y_sd), digits)
  }

  class(out) <- "summary.cos_prediction_areal"
  out
}

print.summary.cos_prediction_areal <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS areal prediction\n")
  } else {
    cat("Non-spatial COS areal prediction\n")
  }
  cat("  spatial uncertainty: ", x$spatial_uncertainty, "\n", sep = "")
  if (!is.null(x$uncertainty$note)) {
    cat("  uncertainty:     ", x$uncertainty$note, "\n", sep = "")
  }
  cat("  target:          ", x$target, "\n", sep = "")
  cat("  prediction units:", x$n_U, "\n", sep = "")
  cat("  posterior draws: ", x$n_samples, "\n", sep = "")
  cat("  stored samples:  ", ifelse(x$keep_samples, "yes", "no"), "\n\n", sep = "")

  cat("H_UA row-sum summary:\n")
  print(x$row_sum_summary)

  cat("\nEta_U mean summary:\n")
  print(x$eta_mean_summary)

  cat("\nEta_U SD summary:\n")
  print(x$eta_sd_summary)

  cat("\nOmega_U mean summary:\n")
  print(x$omega_mean_summary)

  cat("\nOmega_U SD summary:\n")
  print(x$omega_sd_summary)

  if (!is.null(x$y_mean_summary)) {
    cat("\nY_U mean summary:\n")
    print(x$y_mean_summary)
    cat("\nY_U SD summary:\n")
    print(x$y_sd_summary)
  }

  cat("\nPrediction summaries are available in $prediction_summary.\n")

  invisible(x)
}

summary.cos_cv_observed <- function(object, digits = 3, ...) {
  cv <- object
  spatial <- isTRUE(cv$spatial)

  out <- list(
    call = cv$call,
    k = cv$k,
    target = cv$target,
    spatial = spatial,
    n = length(unique(cv$predictions$B_id)),
    summary = cv$summary,
    fold_summary = cv$fold_summary,
    predictions = cv$predictions
  )

  out$summary[, c("RMSPE", "CRPS", "coverage_95",
                  "mean_interval_width", "mean_error", "MAE",
                  "bias", "cor", "slope")] <-
    round(out$summary[, c("RMSPE", "CRPS", "coverage_95",
                          "mean_interval_width", "mean_error", "MAE",
                          "bias", "cor", "slope"),
                      drop = FALSE], digits)

  num_cols <- vapply(out$fold_summary, is.numeric, logical(1))
  out$fold_summary[, num_cols] <- round(out$fold_summary[, num_cols, drop = FALSE], digits)

  class(out) <- "summary.cos_cv_observed"
  out
}

print.summary.cos_cv_observed <- function(x, ...) {
  if (x$spatial) {
    cat("Tapered COS observed-support K-fold validation\n")
  } else {
    cat("Non-spatial COS observed-support K-fold validation\n")
  }
  cat("  folds:            ", x$k, "\n", sep = "")
  cat("  targets:          ", paste(x$target, collapse = ", "), "\n", sep = "")
  cat("  held-out units:   ", x$n, "\n\n", sep = "")

  cat("Overall prediction scores:\n")
  print(x$summary, row.names = FALSE)

  cat("\nFold summaries are available in $fold_summary.")
  cat("\nHeld-out predictions are available in $predictions.\n")

  invisible(x)
}
