cos_predict_areal <- function(fit,
                              rec_B,
                              U_blocks,
                              offset_U = NULL,
                              spatial_uncertainty = c("conditional_mean", "marginal"),
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
  if (!inherits(U_blocks, "cos_areal_blocks")) {
    stop("U_blocks must be a cos_areal_blocks object from cos_make_areal_blocks().")
  }

  spatial_uncertainty <- match.arg(spatial_uncertainty)
  target <- match.arg(target)
  keep_samples <- isTRUE(keep_samples)
  family <- fit$family
  if (is.null(family)) family <- "gaussian"
  if (!family %in% c("gaussian", "negative_binomial") && target == "observed") {
    stop("target = 'observed' is currently only implemented for Gaussian and negative-binomial fits.")
  }
  offset_supplied <- !is.null(offset_U)

  # ------------------------------------------------
  # Pull objects into local variables
  # ------------------------------------------------

  theta_samples <- as.data.frame(rec_B$theta_samples)
  beta_samples <- as.matrix(rec_B$beta_samples)
  omega_B_samples <- as.matrix(rec_B$omega_B_samples)
  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE

  if (verbose && spatial_uncertainty == "marginal") {
    if (spatial) {
      message("spatial_uncertainty = 'marginal' draws independent one-dimensional conditional Gaussians for each areal unit and posterior sample.")
    } else {
      message("spatial_uncertainty = 'marginal' has no spatial effect to sample for non-spatial fits; using beta draws only.")
    }
  }
  if (spatial && spatial_uncertainty == "conditional_mean") {
    warning(
      "spatial_uncertainty = 'conditional_mean' uses the conditional mean of unobserved spatial effects; ",
      "reported eta_sd omits conditional spatial prediction uncertainty. ",
      "Use spatial_uncertainty = 'marginal' for marginal areal posterior predictive uncertainty.",
      call. = FALSE
    )
  }
  if (verbose && spatial && spatial_uncertainty == "marginal" && keep_samples) {
    warning(
      "Areal samples are marginal by prediction unit, not joint samples across units; ",
      "do not use eta_samples or y_samples for totals, contrasts, or rankings that require cross-unit covariance.",
      call. = FALSE
    )
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

  if (!is.null(offset_U) && is.character(offset_U) && length(offset_U) == 1L) {
    if (is.null(U_blocks$U_sf) || !offset_U %in% names(U_blocks$U_sf)) {
      stop("offset_U column not found in U_blocks$U_sf.")
    }
    offset_U <- U_blocks$U_sf[[offset_U]]
  }
  if (is.null(offset_U)) {
    offset_U <- rep(0, n_u)
  }
  if (length(offset_U) == 1L) {
    offset_U <- rep(offset_U, n_u)
  }
  offset_U <- as.numeric(offset_U)
  if (length(offset_U) != n_u || any(!is.finite(offset_U))) {
    stop("offset_U must be NULL, a column name in U_blocks$U_sf, or a finite numeric vector of length one or length(U_blocks$blocks).")
  }
  if (identical(family, "gaussian") && offset_supplied) {
    stop("offset_U is used only with Polya-Gamma response families.")
  }

  if (n_u == 0L) {
    stop("U_blocks contains no usable areal prediction units.")
  }
  if (nrow(theta_samples) != n_save || nrow(omega_B_samples) != n_save) {
    stop("rec_B theta, beta, and omega_B samples must have the same number of rows.")
  }
  if (spatial && !all(c("sigma_sq", "phi") %in% names(theta_samples))) {
    stop("rec_B$theta_samples must contain sigma_sq and phi.")
  }
  if (ncol(beta_samples) != length(x_names)) {
    stop("Recovered beta samples are not compatible with fit$prep$X.")
  }
  if (!all(x_names %in% U_blocks$x_names)) {
    stop("U_blocks covariates are not compatible with the fitted covariates.")
  }

  if (family %in% c("binomial", "negative_binomial")) {
    is_nb <- identical(family, "negative_binomial")
    if (!is_nb && target != "latent") {
      stop("Binomial areal prediction currently uses target = 'latent' and reports link-scale eta and response probability p.")
    }
    if (is_nb && (is.null(fit$size) || length(fit$size) != 1L || !is.finite(fit$size) || fit$size <= 0)) {
      stop("negative-binomial fit is missing a finite positive size.")
    }

    omega_mean <- rep(0, n_u)
    eta_mean <- rep(0, n_u)
    p_mean <- rep(0, n_u)
    mu_mean <- rep(0, n_u)
    y_mean <- rep(0, n_u)
    omega_M2 <- rep(0, n_u)
    eta_M2 <- rep(0, n_u)
    p_M2 <- rep(0, n_u)
    mu_M2 <- rep(0, n_u)
    y_M2 <- rep(0, n_u)
    cond_omega_sd_mean <- rep(0, n_u)

    omega_samples <- NULL
    eta_samples <- NULL
    p_samples <- NULL
    mu_samples <- NULL
    y_samples <- NULL
    if (keep_samples) {
      omega_samples <- matrix(NA_real_, n_save, n_u)
      eta_samples <- matrix(NA_real_, n_save, n_u)
      colnames(omega_samples) <- paste0("U_", vapply(blocks, `[[`, integer(1), "U_id"))
      colnames(eta_samples) <- colnames(omega_samples)
      if (is_nb) {
        mu_samples <- matrix(NA_real_, n_save, n_u)
        colnames(mu_samples) <- colnames(omega_samples)
        if (target == "observed") {
          y_samples <- matrix(NA_real_, n_save, n_u)
          colnames(y_samples) <- colnames(omega_samples)
        }
      } else {
        p_samples <- matrix(NA_real_, n_save, n_u)
        colnames(p_samples) <- colnames(omega_samples)
      }
    }

    if (verbose && spatial && spatial_uncertainty == "marginal") {
      message("Polya-Gamma areal prediction uses independent marginal conditional draws for each areal unit; no joint areal covariance is formed.")
    }

    for (ii in seq_len(n_save)) {
      if (verbose && (ii == 1L || ii %% 100L == 0L || ii == n_save)) {
        if (is_nb) {
          message("Predicting negative-binomial areal means: posterior sample ", ii, " of ", n_save)
        } else {
          message("Predicting binomial areal probabilities: posterior sample ", ii, " of ", n_save)
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
        C_B_inv_omega_B <- backsolve(R_C_B, forwardsolve(t(R_C_B), omega_B))
      }

      omega_draw <- numeric(n_u)
      eta_draw <- numeric(n_u)
      p_draw <- numeric(n_u)
      mu_draw <- numeric(n_u)
      y_draw <- numeric(n_u)

      for (kk in seq_len(n_u)) {
        U_k <- blocks[[kk]]

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
          cond_var_cor <- max(0, C_UU - as.numeric(crossprod(C_UB, C_B_inv_C_BU)))
          cond_omega_sd <- sqrt(sigma_sq * cond_var_cor)
          cond_omega_sd_mean[kk] <- cond_omega_sd_mean[kk] + cond_omega_sd / n_save

          if (spatial_uncertainty == "marginal") {
            omega_U <- omega_U_mean + stats::rnorm(1L, mean = 0, sd = cond_omega_sd)
          } else {
            omega_U <- omega_U_mean
          }
        } else {
          omega_U <- 0
        }

        omega_draw[kk] <- omega_U
        eta_draw[kk] <- offset_U[kk] + as.numeric(crossprod(U_k$X_U[x_names], beta)) + omega_U
        if (is_nb) {
          mu_draw[kk] <- exp(eta_draw[kk])
          if (target == "observed") {
            y_draw[kk] <- stats::rnbinom(1L, size = fit$size, mu = mu_draw[kk])
          }
        } else {
          p_draw[kk] <- stats::plogis(eta_draw[kk])
        }
      }

      if (keep_samples) {
        omega_samples[ii, ] <- omega_draw
        eta_samples[ii, ] <- eta_draw
        if (is_nb) {
          mu_samples[ii, ] <- mu_draw
          if (target == "observed") y_samples[ii, ] <- y_draw
        } else {
          p_samples[ii, ] <- p_draw
        }
      }

      delta_omega <- omega_draw - omega_mean
      omega_mean <- omega_mean + delta_omega / ii
      omega_M2 <- omega_M2 + delta_omega * (omega_draw - omega_mean)

      delta_eta <- eta_draw - eta_mean
      eta_mean <- eta_mean + delta_eta / ii
      eta_M2 <- eta_M2 + delta_eta * (eta_draw - eta_mean)

      if (is_nb) {
        delta_mu <- mu_draw - mu_mean
        mu_mean <- mu_mean + delta_mu / ii
        mu_M2 <- mu_M2 + delta_mu * (mu_draw - mu_mean)

        if (target == "observed") {
          delta_y <- y_draw - y_mean
          y_mean <- y_mean + delta_y / ii
          y_M2 <- y_M2 + delta_y * (y_draw - y_mean)
        }
      } else {
        delta_p <- p_draw - p_mean
        p_mean <- p_mean + delta_p / ii
        p_M2 <- p_M2 + delta_p * (p_draw - p_mean)
      }
    }

    pred_summary <- data.frame(
      U_id = vapply(blocks, `[[`, integer(1), "U_id"),
      U_label = vapply(blocks, `[[`, character(1), "U_label"),
      n_cell = vapply(blocks, `[[`, integer(1), "n_cell"),
      row_sum = vapply(blocks, `[[`, numeric(1), "row_sum"),
      offset = offset_U,
      omega_mean = omega_mean,
      omega_sd = sqrt(omega_M2 / pmax(n_save - 1L, 1L)),
      eta_mean = eta_mean,
      eta_sd = sqrt(eta_M2 / pmax(n_save - 1L, 1L)),
      conditional_omega_sd_mean = cond_omega_sd_mean,
      stringsAsFactors = FALSE
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
      U_blocks = U_blocks,
      offset_U = offset_U,
      summary = pred_summary,
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
            "Mean count summaries use independent marginal conditional draws for each areal unit; no joint areal covariance is formed."
          } else {
            "Response probability summaries use independent marginal conditional draws for each areal unit; no joint areal covariance is formed."
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
      keep_samples = keep_samples,
      n_U = n_u,
      n_samples = n_save,
      omega_samples = omega_samples,
      eta_samples = eta_samples,
      p_samples = p_samples,
      mu_samples = mu_samples,
      y_samples = y_samples,
      spatial = spatial,
      target = target,
      timing = timing,
      call = match.call()
    )

    class(out) <- "cos_prediction_areal"
    return(out)
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

        if (spatial_uncertainty == "marginal") {
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

  dt <- proc.time() - t_start
  timing <- c(
    user = unname(dt[["user.self"]]),
    system = unname(dt[["sys.self"]]),
    elapsed = unname(dt[["elapsed"]])
  )

  out <- list(
    fit = fit,
    rec_B = rec_B,
    U_blocks = U_blocks,
    summary = pred_summary,
    spatial_uncertainty = spatial_uncertainty,
    uncertainty = list(
      spatial_effect = if (spatial && spatial_uncertainty == "conditional_mean") {
        "conditional_mean_only"
      } else if (spatial) {
        "marginal_conditional_sample"
      } else {
        "none_nonspatial"
      },
      joint_samples = FALSE,
      note = if (spatial && spatial_uncertainty == "conditional_mean") {
        "eta_sd omits conditional spatial prediction uncertainty."
      } else if (spatial) {
        "Areal samples include marginal conditional uncertainty for each unit but are not joint samples across units."
      } else {
        "No residual spatial effect is present in the nonspatial model."
      }
    ),
    keep_samples = keep_samples,
    n_U = n_u,
    n_samples = n_save,
    omega_samples = omega_samples,
    eta_samples = eta_samples,
    y_samples = y_samples,
    spatial = spatial,
    target = target,
    timing = timing,
    call = match.call()
  )

  class(out) <- "cos_prediction_areal"
  out
}
