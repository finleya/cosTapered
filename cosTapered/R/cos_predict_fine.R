cos_predict_fine <- function(fit,
                             rec_B,
                             pred_coords = NULL,
                             X_pred = NULL,
                             offset_pred = NULL,
                             spatial_uncertainty = c("conditional_mean", "marginal", "joint"),
                             target = c("latent", "observed"),
                             keep_samples = FALSE,
                             n_threads = NULL,
                             verbose = TRUE) {
  t_start <- proc.time()

  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(fit, "cos_fit")) {
    stop("fit must be a cos_fit object from cos_fit().")
  }
  if (!inherits(rec_B, "cos_recovery_B")) {
    stop("rec_B must be a cos_recovery_B object from cos_recover_B().")
  }

  spatial_uncertainty <- match.arg(spatial_uncertainty)
  target <- match.arg(target)
  keep_samples <- isTRUE(keep_samples)
  if (!fit$family %in% c("gaussian", "negative_binomial") && target == "observed") {
    stop("target = 'observed' is currently only implemented for Gaussian and negative-binomial fits.")
  }
  offset_supplied <- !is.null(offset_pred)

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

  if (is.null(offset_pred)) {
    offset_pred <- rep(0, nrow(X_pred))
  }
  if (length(offset_pred) == 1L) {
    offset_pred <- rep(offset_pred, nrow(X_pred))
  }
  offset_pred <- as.numeric(offset_pred)
  if (length(offset_pred) != nrow(X_pred) || any(!is.finite(offset_pred))) {
    stop("offset_pred must be NULL or a finite numeric vector of length one or nrow(X_pred).")
  }
  if (identical(fit$family, "gaussian") && offset_supplied) {
    stop("offset_pred is used only with Polya-Gamma response families.")
  }

  # ------------------------------------------------
  # Pull objects into local variables
  # ------------------------------------------------

  theta_samples <- as.data.frame(rec_B$theta_samples)
  beta_samples <- as.matrix(rec_B$beta_samples)
  omega_B_samples <- as.matrix(rec_B$omega_B_samples)
  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE

  if (verbose && spatial_uncertainty %in% c("marginal", "joint")) {
    if (spatial) {
      if (spatial_uncertainty == "joint") {
        message("spatial_uncertainty = 'joint' forms dense n_pred x n_pred covariance matrices; ensure sufficient memory is available.")
      } else {
        message("spatial_uncertainty = 'marginal' draws independent marginal conditional values for each prediction cell.")
      }
    } else {
      message("No spatial effect is present for non-spatial fits; using beta draws only.")
    }
  }
  if (spatial && spatial_uncertainty == "conditional_mean") {
    warning(
      "spatial_uncertainty = 'conditional_mean' uses the conditional mean of unobserved spatial effects; ",
      "reported eta_sd omits conditional spatial prediction uncertainty. ",
      "Use spatial_uncertainty = 'marginal' or 'joint' when posterior predictive uncertainty is needed.",
      call. = FALSE
    )
  }

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

  n_save <- nrow(beta_samples)
  n_pred <- nrow(X_pred)
  n_b <- ncol(omega_B_samples)

  if (nrow(theta_samples) != n_save || nrow(omega_B_samples) != n_save) {
    stop("rec_B theta, beta, and omega_B samples must have the same number of rows.")
  }
  if (spatial && !all(c("sigma_sq", "phi") %in% names(theta_samples))) {
    stop("rec_B$theta_samples must contain sigma_sq and phi.")
  }
  if (ncol(beta_samples) != ncol(X_pred)) {
    stop("ncol(X_pred) must match the number of recovered beta coefficients.")
  }
  if (spatial_uncertainty == "joint" && as.double(n_pred) * as.double(n_pred) > .Machine$integer.max) {
    if (identical(fit$family, "gaussian")) {
      stop("spatial_uncertainty = 'joint' requires a dense n_pred x n_pred covariance matrix; use fewer prediction cells or spatial_uncertainty = 'marginal'.")
    }
  }

  if (fit$family %in% c("binomial", "negative_binomial")) {
    is_nb <- identical(fit$family, "negative_binomial")
    if (!is_nb && target != "latent") {
      stop("Binomial fine prediction currently uses target = 'latent' and reports link-scale eta and response probability p.")
    }
    if (spatial_uncertainty == "joint") {
      stop("Polya-Gamma fine prediction currently supports spatial_uncertainty = 'conditional_mean' or 'marginal'.")
    }
    if (is_nb && (is.null(fit$size) || length(fit$size) != 1L || !is.finite(fit$size) || fit$size <= 0)) {
      stop("negative-binomial fit is missing a finite positive size.")
    }

    eta_mean <- rep(0, n_pred)
    eta_M2 <- rep(0, n_pred)
    omega_mean <- rep(0, n_pred)
    omega_M2 <- rep(0, n_pred)
    p_mean <- rep(0, n_pred)
    p_M2 <- rep(0, n_pred)
    mu_mean <- rep(0, n_pred)
    mu_M2 <- rep(0, n_pred)
    y_mean <- rep(0, n_pred)
    y_M2 <- rep(0, n_pred)
    cond_omega_sd_mean <- rep(0, n_pred)

    omega_samples <- NULL
    eta_samples <- NULL
    p_samples <- NULL
    mu_samples <- NULL
    y_samples <- NULL
    if (keep_samples) {
      omega_samples <- matrix(NA_real_, n_save, n_pred)
      eta_samples <- matrix(NA_real_, n_save, n_pred)
      colnames(omega_samples) <- paste0("pred_", seq_len(n_pred))
      colnames(eta_samples) <- colnames(omega_samples)
      if (is_nb) {
        mu_samples <- matrix(NA_real_, n_save, n_pred)
        colnames(mu_samples) <- colnames(omega_samples)
        if (target == "observed") {
          y_samples <- matrix(NA_real_, n_save, n_pred)
          colnames(y_samples) <- colnames(omega_samples)
        }
      } else {
        p_samples <- matrix(NA_real_, n_save, n_pred)
        colnames(p_samples) <- colnames(omega_samples)
      }
    }

    if (verbose && spatial && spatial_uncertainty == "marginal") {
      message("Polya-Gamma fine prediction uses independent marginal conditional draws for each prediction cell; no joint fine-grid covariance is formed.")
    }

    for (ii in seq_len(n_save)) {
      if (verbose && (ii == 1L || ii %% 100L == 0L || ii == n_save)) {
        if (is_nb) {
          message("Predicting negative-binomial fine means: ", ii, " of ", n_save)
        } else {
          message("Predicting binomial fine probabilities: ", ii, " of ", n_save)
        }
      }

      beta <- as.numeric(beta_samples[ii, ])

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

        C_B_inv_C_B_pred <- backsolve(R_C_B, forwardsolve(t(R_C_B), t(C_pred_B)))
        cond_var_cor <- pmax(0, 1 - rowSums(C_pred_B * t(C_B_inv_C_B_pred)))
        cond_omega_sd <- sqrt(sigma_sq * cond_var_cor)
        cond_omega_sd_mean <- cond_omega_sd_mean + cond_omega_sd / n_save

        if (spatial_uncertainty == "marginal") {
          omega_pred <- omega_pred_mean + cond_omega_sd * stats::rnorm(n_pred)
        } else {
          omega_pred <- omega_pred_mean
        }
      } else {
        omega_pred <- rep(0, n_pred)
        cond_omega_sd <- rep(0, n_pred)
      }

      eta_pred <- as.numeric(offset_pred + X_pred %*% beta + omega_pred)
      if (is_nb) {
        mu_pred <- exp(eta_pred)
        if (target == "observed") {
          y_pred <- stats::rnbinom(n_pred, size = fit$size, mu = mu_pred)
        }
      } else {
        p_pred <- stats::plogis(eta_pred)
      }

      if (keep_samples) {
        omega_samples[ii, ] <- omega_pred
        eta_samples[ii, ] <- eta_pred
        if (is_nb) {
          mu_samples[ii, ] <- mu_pred
          if (target == "observed") y_samples[ii, ] <- y_pred
        } else {
          p_samples[ii, ] <- p_pred
        }
      }

      delta_omega <- omega_pred - omega_mean
      omega_mean <- omega_mean + delta_omega / ii
      omega_M2 <- omega_M2 + delta_omega * (omega_pred - omega_mean)

      delta_eta <- eta_pred - eta_mean
      eta_mean <- eta_mean + delta_eta / ii
      eta_M2 <- eta_M2 + delta_eta * (eta_pred - eta_mean)

      if (is_nb) {
        delta_mu <- mu_pred - mu_mean
        mu_mean <- mu_mean + delta_mu / ii
        mu_M2 <- mu_M2 + delta_mu * (mu_pred - mu_mean)

        if (target == "observed") {
          delta_y <- y_pred - y_mean
          y_mean <- y_mean + delta_y / ii
          y_M2 <- y_M2 + delta_y * (y_pred - y_mean)
        }
      } else {
        delta_p <- p_pred - p_mean
        p_mean <- p_mean + delta_p / ii
        p_M2 <- p_M2 + delta_p * (p_pred - p_mean)
      }
    }

    pred_summary <- data.frame(
      pred_id = seq_len(n_pred),
      x = pred_coords[, 1],
      y = pred_coords[, 2],
      offset = offset_pred,
      omega_mean = omega_mean,
      omega_sd = sqrt(omega_M2 / pmax(n_save - 1L, 1L)),
      eta_mean = eta_mean,
      eta_sd = sqrt(eta_M2 / pmax(n_save - 1L, 1L)),
      conditional_omega_sd_mean = cond_omega_sd_mean
    )
    if (is_nb) {
      pred_summary$mu_mean <- mu_mean
      pred_summary$mu_sd <- sqrt(mu_M2 / pmax(n_save - 1L, 1L))
      if (target == "observed") {
        pred_summary$y_mean <- y_mean
        pred_summary$y_sd <- sqrt(y_M2 / pmax(n_save - 1L, 1L))
      }
    } else {
      pred_summary$p_mean <- p_mean
      pred_summary$p_sd <- sqrt(p_M2 / pmax(n_save - 1L, 1L))
    }

    if (keep_samples) {
      if (is_nb) {
        pred_summary$mu_q025 <- apply(mu_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
        pred_summary$mu_q50 <- apply(mu_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
        pred_summary$mu_q975 <- apply(mu_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
        if (target == "observed") {
          pred_summary$y_q025 <- apply(y_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
          pred_summary$y_q50 <- apply(y_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
          pred_summary$y_q975 <- apply(y_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
        }
      } else {
        pred_summary$p_q025 <- apply(p_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
        pred_summary$p_q50 <- apply(p_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
        pred_summary$p_q975 <- apply(p_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
      }
      pred_summary$eta_q025 <- apply(eta_samples, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      pred_summary$eta_q50 <- apply(eta_samples, 2, stats::quantile, probs = 0.5, na.rm = TRUE)
      pred_summary$eta_q975 <- apply(eta_samples, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    }

    dt <- proc.time() - t_start
    timing <- c(
      user = unname(dt[["user.self"]]),
      system = unname(dt[["sys.self"]]),
      elapsed = unname(dt[["elapsed"]])
    )

    out <- list(
      fit = fit,
      rec_B = rec_B,
      summary = pred_summary,
      coords = pred_coords,
      X_pred = X_pred,
      offset_pred = offset_pred,
      spatial_uncertainty = spatial_uncertainty,
      uncertainty = list(
        spatial_effect = if (spatial && spatial_uncertainty == "marginal") {
          "marginal_conditional_sample"
        } else if (spatial) {
          "conditional_mean_only"
        } else {
          "none_nonspatial"
        },
        joint_samples = FALSE,
        note = if (spatial && spatial_uncertainty == "marginal") {
          if (is_nb) {
            "Mean count summaries use independent marginal conditional draws for each prediction cell; no joint fine-grid covariance is formed."
          } else {
            "Response probability summaries use independent marginal conditional draws for each prediction cell; no joint fine-grid covariance is formed."
          }
        } else if (spatial) {
          if (is_nb) {
            "Mean count summaries use the conditional mean of unobserved spatial effects."
          } else {
            "Response probability summaries use the conditional mean of unobserved spatial effects."
          }
        } else {
          "No residual spatial effect is present in the nonspatial model."
        }
      ),
      target = target,
      keep_samples = keep_samples,
      n_pred = n_pred,
      n_samples = n_save,
      omega_samples = omega_samples,
      eta_samples = eta_samples,
      p_samples = p_samples,
      mu_samples = mu_samples,
      y_samples = y_samples,
      spatial = spatial,
      timing = timing,
      call = match.call()
    )

    class(out) <- "cos_prediction_fine"
    return(out)
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

      if (spatial_uncertainty == "conditional_mean") {
        omega_pred <- omega_pred_mean
      } else if (spatial_uncertainty == "joint") {
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
              "This can happen when spatial_uncertainty = 'joint' predicts fine cells constrained by observed-support effects.",
              call. = FALSE
            )
            psd_warning_given <- TRUE
          }

          eig_values <- pmax(eig$values, 0)
          omega_pred <- omega_pred_mean +
            sqrt(sigma_sq) * as.numeric(eig$vectors %*% (sqrt(eig_values) * stats::rnorm(n_pred)))
        }
      } else {
        C_B_inv_C_B_pred <- backsolve(R_C_B, forwardsolve(t(R_C_B), t(C_pred_B)))
        cond_var_cor <- pmax(0, 1 - rowSums(C_pred_B * t(C_B_inv_C_B_pred)))
        omega_pred <- omega_pred_mean + sqrt(sigma_sq * cond_var_cor) * stats::rnorm(n_pred)
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

  dt <- proc.time() - t_start
  timing <- c(
    user = unname(dt[["user.self"]]),
    system = unname(dt[["sys.self"]]),
    elapsed = unname(dt[["elapsed"]])
  )

  out <- list(
    fit = fit,
    rec_B = rec_B,
    summary = pred_summary,
    coords = pred_coords,
    X_pred = X_pred,
    spatial_uncertainty = spatial_uncertainty,
    uncertainty = list(
      spatial_effect = if (spatial && spatial_uncertainty == "conditional_mean") {
        "conditional_mean_only"
      } else if (spatial && spatial_uncertainty == "marginal") {
        "marginal_conditional_sample"
      } else if (spatial) {
        "joint_conditional_sample"
      } else {
        "none_nonspatial"
      },
      joint_samples = spatial && spatial_uncertainty == "joint",
      note = if (spatial && spatial_uncertainty == "conditional_mean") {
        "eta_sd omits conditional spatial prediction uncertainty."
      } else if (spatial && spatial_uncertainty == "marginal") {
        "Samples include independent marginal conditional spatial prediction uncertainty for each fine prediction cell; they are not joint samples across cells."
      } else if (spatial) {
        "Samples include joint conditional spatial prediction uncertainty across fine prediction cells."
      } else {
        "No residual spatial effect is present in the nonspatial model."
      }
    ),
    target = target,
    keep_samples = keep_samples,
    n_pred = n_pred,
    n_samples = n_save,
    omega_samples = omega_samples,
    eta_samples = eta_samples,
    y_samples = y_samples,
    spatial = spatial,
    timing = timing,
    call = match.call()
  )

  class(out) <- "cos_prediction_fine"
  out
}
