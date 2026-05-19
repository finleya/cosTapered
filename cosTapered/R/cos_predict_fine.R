cos_predict_fine <- function(fit,
                             rec_B,
                             pred_coords = NULL,
                             X_pred = NULL,
                             method = c("mean", "sample"),
                             target = c("latent", "observed"),
                             keep_samples = FALSE,
                             verbose = TRUE) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(fit, "cos_fit")) {
    stop("fit must be a cos_fit object from cos_fit().")
  }
  if (!inherits(rec_B, "cos_recovery_B")) {
    stop("rec_B must be a cos_recovery_B object from cos_recover_B().")
  }

  method <- match.arg(method)
  target <- match.arg(target)
  keep_samples <- isTRUE(keep_samples)

  if (is.null(pred_coords) && is.null(X_pred)) {
    pred_coords <- fit$prep$A_coords
    X_pred <- fit$prep$X
  } else if (is.null(pred_coords) || is.null(X_pred)) {
    stop("Provide both pred_coords and X_pred, or neither to use the original fine support.")
  }

  pred_coords <- as.matrix(pred_coords)
  X_pred <- as.matrix(X_pred)
  storage.mode(pred_coords) <- "double"
  storage.mode(X_pred) <- "double"

  if (ncol(pred_coords) != 2L) {
    stop("pred_coords must have two columns containing x and y coordinates.")
  }
  if (nrow(pred_coords) != nrow(X_pred)) {
    stop("nrow(pred_coords) must equal nrow(X_pred).")
  }

  x_names <- colnames(fit$prep$X)
  if (is.null(colnames(X_pred))) {
    stop("X_pred must have column names matching the fitted covariates.")
  }
  if (!all(x_names %in% colnames(X_pred))) {
    stop("X_pred must contain columns: ", paste(x_names, collapse = ", "))
  }
  X_pred <- X_pred[, x_names, drop = FALSE]

  # ------------------------------------------------
  # Pull objects into local variables
  # ------------------------------------------------

  theta_samples <- as.data.frame(rec_B$theta_samples)
  beta_samples <- as.matrix(rec_B$beta_samples)
  omega_B_samples <- as.matrix(rec_B$omega_B_samples)
  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE

  if (verbose && method == "sample") {
    if (spatial) {
      message("method = 'sample' forms dense n_pred x n_pred covariance matrices; ensure sufficient memory is available.")
    } else {
      message("method = 'sample' has no spatial effect to sample for non-spatial fits; using beta draws only.")
    }
  }

  C_B_pairs <- fit$prep$C_B_pairs
  H_BA_comp <- fit$prep$H_BA_comp
  gamma <- fit$prep$gamma
  taper_code <- fit$prep$taper_code
  n_threads <- fit$prep$n_threads

  n_save <- nrow(beta_samples)
  n_pred <- nrow(X_pred)
  n_b <- ncol(omega_B_samples)

  if (nrow(theta_samples) != n_save || nrow(omega_B_samples) != n_save) {
    stop("rec_B theta, beta, and omega_B samples must have the same number of rows.")
  }
  if (!all(c("sigma_sq", "phi") %in% names(theta_samples))) {
    stop("rec_B$theta_samples must contain sigma_sq and phi.")
  }
  if (ncol(beta_samples) != ncol(X_pred)) {
    stop("ncol(X_pred) must match the number of recovered beta coefficients.")
  }
  if (method == "sample" && as.double(n_pred) * as.double(n_pred) > .Machine$integer.max) {
    stop("method = 'sample' requires a dense n_pred x n_pred covariance matrix; use fewer prediction cells or method = 'mean'.")
  }

  # ------------------------------------------------
  # Allocate output storage
  # ------------------------------------------------

  omega_mean <- rep(0, n_pred)
  eta_mean <- rep(0, n_pred)
  y_mean <- rep(0, n_pred)
  omega_M2 <- rep(0, n_pred)
  eta_M2 <- rep(0, n_pred)
  y_M2 <- rep(0, n_pred)

  omega_samples <- NULL
  eta_samples <- NULL
  y_samples <- NULL

  if (keep_samples) {
    omega_samples <- matrix(NA_real_, n_save, n_pred)
    eta_samples <- matrix(NA_real_, n_save, n_pred)
    colnames(omega_samples) <- paste0("pred_", seq_len(n_pred))
    colnames(eta_samples) <- paste0("pred_", seq_len(n_pred))
    if (target == "observed") {
      y_samples <- matrix(NA_real_, n_save, n_pred)
      colnames(y_samples) <- colnames(eta_samples)
    }
  }

  psd_warning_given <- FALSE

  # ------------------------------------------------
  # Loop over retained recovery samples
  # ------------------------------------------------

  for (ii in seq_len(n_save)) {
    if (verbose && (ii == 1L || ii %% 100L == 0L || ii == n_save)) {
      message("Predicting fine surface: ", ii, " of ", n_save)
    }

    beta <- as.numeric(beta_samples[ii, ])
    tau_sq <- theta_samples$tau_sq[ii]

    if (spatial) {
      sigma_sq <- theta_samples$sigma_sq[ii]
      phi <- theta_samples$phi[ii]
      omega_B <- as.numeric(omega_B_samples[ii, ])

      C_B <- make_C_B_tapered_from_pairs_Rcall(
        C_B_pairs = C_B_pairs,
        phi = phi,
        n_threads = n_threads
      )
      R_C_B <- chol(C_B)

      C_pred_B <- make_C_pred_B_tapered_Rcall(
        pred_coords = pred_coords,
        H_BA_comp = H_BA_comp,
        phi = phi,
        gamma = gamma,
        taper_code = taper_code,
        n_threads = n_threads
      )

      C_B_inv_omega_B <- backsolve(R_C_B, forwardsolve(t(R_C_B), omega_B))
      omega_pred_mean <- as.numeric(C_pred_B %*% C_B_inv_omega_B)

      if (method == "mean") {
        omega_pred <- omega_pred_mean
      } else {
        C_pred <- make_C_pred_tapered_dense_Rcall(
          pred_coords = pred_coords,
          phi = phi,
          gamma = gamma,
          taper_code = taper_code,
          n_threads = n_threads
        )
        C_B_inv_C_B_pred <- backsolve(R_C_B, forwardsolve(t(R_C_B), t(C_pred_B)))
        cond_cor <- C_pred - C_pred_B %*% C_B_inv_C_B_pred

        R_cond <- tryCatch(
          chol(sigma_sq * cond_cor),
          error = function(e) NULL
        )

        if (!is.null(R_cond)) {
          omega_pred <- omega_pred_mean + as.numeric(t(R_cond) %*% stats::rnorm(n_pred))
        } else {
          cond_cor_sym <- (cond_cor + t(cond_cor)) / 2
          eig <- eigen(cond_cor_sym, symmetric = TRUE)
          eig_tol <- 1e-8 * max(abs(eig$values), 1)

          if (min(eig$values) < -eig_tol) {
            stop("Fine-support conditional covariance has materially negative eigenvalues.")
          }

          if (!psd_warning_given) {
            warning(
              "Fine-support conditional covariance is positive semidefinite, not positive definite; using eigen-based sampling. ",
              "This can happen when method = 'sample' predicts fine cells constrained by observed-support effects.",
              call. = FALSE
            )
            psd_warning_given <- TRUE
          }

          eig_values <- pmax(eig$values, 0)
          omega_pred <- omega_pred_mean +
            sqrt(sigma_sq) * as.numeric(eig$vectors %*% (sqrt(eig_values) * stats::rnorm(n_pred)))
        }
      }
    } else {
      omega_pred <- rep(0, n_pred)
    }

    eta_pred <- as.numeric(X_pred %*% beta + omega_pred)
    if (target == "observed") {
      y_pred <- eta_pred + stats::rnorm(n_pred, mean = 0, sd = sqrt(tau_sq))
    }

    if (keep_samples) {
      omega_samples[ii, ] <- omega_pred
      eta_samples[ii, ] <- eta_pred
      if (target == "observed") {
        y_samples[ii, ] <- y_pred
      }
    }

    delta_omega <- omega_pred - omega_mean
    omega_mean <- omega_mean + delta_omega / ii
    omega_M2 <- omega_M2 + delta_omega * (omega_pred - omega_mean)

    delta_eta <- eta_pred - eta_mean
    eta_mean <- eta_mean + delta_eta / ii
    eta_M2 <- eta_M2 + delta_eta * (eta_pred - eta_mean)

    if (target == "observed") {
      delta_y <- y_pred - y_mean
      y_mean <- y_mean + delta_y / ii
      y_M2 <- y_M2 + delta_y * (y_pred - y_mean)
    }
  }

  # ------------------------------------------------
  # Prediction summaries
  # ------------------------------------------------

  omega_sd <- sqrt(omega_M2 / pmax(n_save - 1L, 1L))
  eta_sd <- sqrt(eta_M2 / pmax(n_save - 1L, 1L))
  y_sd <- sqrt(y_M2 / pmax(n_save - 1L, 1L))

  pred_summary <- data.frame(
    pred_id = seq_len(n_pred),
    x = pred_coords[, 1],
    y = pred_coords[, 2],
    omega_mean = omega_mean,
    omega_sd = omega_sd,
    eta_mean = eta_mean,
    eta_sd = eta_sd
  )

  if (target == "observed") {
    pred_summary$y_mean <- y_mean
    pred_summary$y_sd <- y_sd
  }

  if (keep_samples) {
    pred_summary$omega_q025 <- apply(omega_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    pred_summary$omega_q50 <- apply(omega_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
    pred_summary$omega_q975 <- apply(omega_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    pred_summary$eta_q025 <- apply(eta_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    pred_summary$eta_q50 <- apply(eta_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
    pred_summary$eta_q975 <- apply(eta_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    if (target == "observed") {
      pred_summary$y_q025 <- apply(y_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      pred_summary$y_q50 <- apply(y_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
      pred_summary$y_q975 <- apply(y_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    }
  }

  # ------------------------------------------------
  # Return prediction object
  # ------------------------------------------------

  out <- list(
    fit = fit,
    rec_B = rec_B,
    summary = pred_summary,
    coords = pred_coords,
    X_pred = X_pred,
    method = method,
    target = target,
    keep_samples = keep_samples,
    n_pred = n_pred,
    n_samples = n_save,
    omega_samples = omega_samples,
    eta_samples = eta_samples,
    y_samples = y_samples,
    spatial = spatial,
    call = match.call()
  )

  class(out) <- "cos_prediction_fine"
  out
}
