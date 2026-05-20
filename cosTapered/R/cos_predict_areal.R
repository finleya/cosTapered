cos_predict_areal <- function(fit,
                              rec_B,
                              U_blocks,
                              method = c("mean", "sample"),
                              target = c("latent", "observed"),
                              keep_samples = FALSE,
                              n_threads = NULL,
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
  if (!inherits(U_blocks, "cos_areal_blocks")) {
    stop("U_blocks must be a cos_areal_blocks object from cos_make_areal_blocks().")
  }

  method <- match.arg(method)
  target <- match.arg(target)
  keep_samples <- isTRUE(keep_samples)

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
      message("method = 'sample' draws from the one-dimensional conditional Gaussian for each areal unit and posterior sample.")
    } else {
      message("method = 'sample' has no spatial effect to sample for non-spatial fits; using beta draws only.")
    }
  }

  blocks <- U_blocks$blocks
  C_B_pairs <- fit$prep$C_B_pairs
  H_BA_comp <- fit$prep$H_BA_comp
  gamma <- fit$prep$gamma
  taper_code <- fit$prep$taper_code
  if (is.null(n_threads)) {
    n_threads <- fit$prep$n_threads
  } else {
    n_threads <- as.integer(n_threads)
    if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L
  }

  x_names <- colnames(fit$prep$X)
  n_save <- nrow(beta_samples)
  n_u <- length(blocks)
  n_b <- ncol(omega_B_samples)

  if (n_u == 0L) {
    stop("U_blocks contains no usable areal prediction units.")
  }
  if (nrow(theta_samples) != n_save || nrow(omega_B_samples) != n_save) {
    stop("rec_B theta, beta, and omega_B samples must have the same number of rows.")
  }
  if (!all(c("sigma_sq", "phi") %in% names(theta_samples))) {
    stop("rec_B$theta_samples must contain sigma_sq and phi.")
  }
  if (ncol(beta_samples) != length(x_names)) {
    stop("Recovered beta samples are not compatible with fit$prep$X.")
  }
  if (!all(x_names %in% U_blocks$x_names)) {
    stop("U_blocks covariates are not compatible with the fitted covariates.")
  }

  # ------------------------------------------------
  # Allocate output storage
  # ------------------------------------------------

  omega_mean <- rep(0, n_u)
  eta_mean <- rep(0, n_u)
  y_mean <- rep(0, n_u)
  omega_M2 <- rep(0, n_u)
  eta_M2 <- rep(0, n_u)
  y_M2 <- rep(0, n_u)

  omega_samples <- NULL
  eta_samples <- NULL
  y_samples <- NULL

  if (keep_samples) {
    omega_samples <- matrix(NA_real_, n_save, n_u)
    eta_samples <- matrix(NA_real_, n_save, n_u)
    colnames(omega_samples) <- paste0("U_", vapply(blocks, `[[`, integer(1), "U_id"))
    colnames(eta_samples) <- colnames(omega_samples)
    if (target == "observed") {
      y_samples <- matrix(NA_real_, n_save, n_u)
      colnames(y_samples) <- colnames(omega_samples)
    }
  }

  # ------------------------------------------------
  # Loop over retained recovery samples
  # ------------------------------------------------

  for (ii in seq_len(n_save)) {
    if (verbose && (ii == 1L || ii %% 100L == 0L || ii == n_save)) {
      message("Predicting areal units: posterior sample ", ii, " of ", n_save)
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

      C_B_inv_omega_B <- backsolve(R_C_B, forwardsolve(t(R_C_B), omega_B))
    }

    omega_draw <- numeric(n_u)
    eta_draw <- numeric(n_u)
    y_draw <- numeric(n_u)

    for (kk in seq_len(n_u)) {
      U_k <- blocks[[kk]]

      if (verbose && (ii == 1L || ii == n_save)) {
        message(
          "  areal unit ", kk, " of ", n_u,
          " (U_label = ", U_k$U_label, ", n_cell = ", U_k$n_cell, ")"
        )
      }

      if (spatial) {
        C_pixel_B <- make_C_pred_B_tapered_Rcall(
          pred_coords = U_k$coords,
          H_BA_comp = H_BA_comp,
          phi = phi,
          gamma = gamma,
          taper_code = taper_code,
          n_threads = n_threads
        )

        C_UB <- as.numeric(crossprod(U_k$h, C_pixel_B))
        omega_U_mean <- as.numeric(crossprod(C_UB, C_B_inv_omega_B))
        omega_U <- omega_U_mean

        if (method == "sample") {
          plot_start <- c(1L, length(U_k$h) + 1L)

          C_UU <- make_C_B_tapered_Rcall(
            plot_start = plot_start,
            h = U_k$h,
            x = U_k$coords[, 1],
            y = U_k$coords[, 2],
            phi = phi,
            gamma = gamma,
            taper_code = taper_code,
            n_threads = n_threads
          )
          C_UU <- as.numeric(C_UU[1, 1])

          C_B_inv_C_BU <- backsolve(R_C_B, forwardsolve(t(R_C_B), C_UB))
          cond_var_cor <- C_UU - as.numeric(crossprod(C_UB, C_B_inv_C_BU))
          cond_var <- sigma_sq * max(0, cond_var_cor)

          omega_U <- omega_U_mean + stats::rnorm(1L, mean = 0, sd = sqrt(cond_var))
        }
      } else {
        omega_U <- 0
      }

      omega_draw[kk] <- omega_U
      eta_draw[kk] <- as.numeric(crossprod(U_k$X_U[x_names], beta)) + omega_U

      if (target == "observed") {
        d_U <- sum(U_k$h^2)
        y_draw[kk] <- eta_draw[kk] +
          stats::rnorm(1L, mean = 0, sd = sqrt(tau_sq * d_U))
      }
    }

    if (keep_samples) {
      omega_samples[ii, ] <- omega_draw
      eta_samples[ii, ] <- eta_draw
      if (target == "observed") {
        y_samples[ii, ] <- y_draw
      }
    }

    delta_omega <- omega_draw - omega_mean
    omega_mean <- omega_mean + delta_omega / ii
    omega_M2 <- omega_M2 + delta_omega * (omega_draw - omega_mean)

    delta_eta <- eta_draw - eta_mean
    eta_mean <- eta_mean + delta_eta / ii
    eta_M2 <- eta_M2 + delta_eta * (eta_draw - eta_mean)

    if (target == "observed") {
      delta_y <- y_draw - y_mean
      y_mean <- y_mean + delta_y / ii
      y_M2 <- y_M2 + delta_y * (y_draw - y_mean)
    }
  }

  # ------------------------------------------------
  # Prediction summaries
  # ------------------------------------------------

  omega_sd <- sqrt(omega_M2 / pmax(n_save - 1L, 1L))
  eta_sd <- sqrt(eta_M2 / pmax(n_save - 1L, 1L))
  y_sd <- sqrt(y_M2 / pmax(n_save - 1L, 1L))

  pred_summary <- data.frame(
    U_id = vapply(blocks, `[[`, integer(1), "U_id"),
    U_label = vapply(blocks, `[[`, character(1), "U_label"),
    n_cell = vapply(blocks, `[[`, integer(1), "n_cell"),
    row_sum = vapply(blocks, `[[`, numeric(1), "row_sum"),
    omega_mean = omega_mean,
    omega_sd = omega_sd,
    eta_mean = eta_mean,
    eta_sd = eta_sd,
    stringsAsFactors = FALSE
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
    U_blocks = U_blocks,
    summary = pred_summary,
    method = method,
    keep_samples = keep_samples,
    n_U = n_u,
    n_samples = n_save,
    omega_samples = omega_samples,
    eta_samples = eta_samples,
    y_samples = y_samples,
    spatial = spatial,
    target = target,
    call = match.call()
  )

  class(out) <- "cos_prediction_areal"
  out
}
