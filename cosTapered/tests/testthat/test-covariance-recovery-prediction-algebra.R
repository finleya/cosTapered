test_that("full and pair-based tapered observed-support covariance agree", {
  plot_start <- c(1L, 3L, 6L, 7L)
  h <- c(0.6, 0.4, 0.25, 0.5, 0.25, 1)
  x <- c(0, 1, 1, 0, 2, 2)
  y <- c(0, 0, 0, 1, 1, 1)
  gamma <- 2.5
  phi <- 0.4

  for (taper_code in c(1L, 2L)) {
    C_full <- make_C_B_tapered_Rcall(
      plot_start = plot_start,
      h = h,
      x = x,
      y = y,
      phi = phi,
      gamma = gamma,
      taper_code = taper_code,
      n_threads = 1L
    )

    C_pairs <- prep_C_B_tapered_pairs_Rcall(
      plot_start = plot_start,
      h = h,
      x = x,
      y = y,
      gamma = gamma,
      taper_code = taper_code,
      n_threads = 1L
    )
    C_fast <- make_C_B_tapered_from_pairs_Rcall(
      C_B_pairs = C_pairs,
      phi = phi,
      n_threads = 1L
    )

    taper_fun <- function(d) {
      r <- d / gamma
      if (taper_code == 1L) {
        out <- (1 - r)^4 * (1 + 4 * r)
      } else {
        out <- 1 - 1.5 * r + 0.5 * r^3
      }
      ifelse(r < 1, out, 0)
    }

    C_ref <- matrix(0, 3, 3)
    for (ll in seq_len(3)) {
      l_ind <- plot_start[ll]:(plot_start[ll + 1L] - 1L)
      for (kk in seq_len(3)) {
        k_ind <- plot_start[kk]:(plot_start[kk + 1L] - 1L)
        for (ii in l_ind) {
          for (jj in k_ind) {
            d <- sqrt((x[ii] - x[jj])^2 + (y[ii] - y[jj])^2)
            C_ref[ll, kk] <- C_ref[ll, kk] +
              h[ii] * h[jj] * exp(-phi * d) * taper_fun(d)
          }
        }
      }
    }

    expect_equal(C_fast, C_full, tolerance = 1e-12)
    expect_equal(C_full, C_ref, tolerance = 1e-12)
  }
})

test_that("Gaussian recovery means match direct precision calculation", {
  plot_start <- c(1L, 3L, 6L, 7L)
  h <- c(0.6, 0.4, 0.25, 0.5, 0.25, 1)
  x <- c(0, 1, 1, 0, 2, 2)
  y <- c(0, 0, 0, 1, 1, 1)
  gamma <- 2.5
  phi <- 0.35

  C_B_pairs <- prep_C_B_tapered_pairs_Rcall(
    plot_start = plot_start,
    h = h,
    x = x,
    y = y,
    gamma = gamma,
    taper_code = 1L,
    n_threads = 1L
  )

  y_B <- c(1.2, -0.5, 0.8)
  X_B <- cbind(intercept = 1, x = c(-0.3, 0.4, 1.1))
  D_h <- matrix(c(
    0.52, 0.10, 0.00,
    0.10, 0.375, 0.25,
    0.00, 0.25, 1.00
  ), 3, 3, byrow = TRUE)
  mu_beta <- c(intercept = 0.2, x = -0.1)
  V_beta <- matrix(c(1.5, 0.2, 0.2, 1.0), 2, 2)
  tau_sq <- 0.7
  sigma_sq <- 1.3

  fit <- list(
    prep = list(
      y_B = y_B,
      X_B = X_B,
      D_h = D_h,
      C_B_pairs = C_B_pairs,
      n_threads = 1L
    ),
    priors = list(beta = list(mu = mu_beta, V = V_beta)),
    theta_samples = data.frame(
      chain = 1L,
      iter = 1L,
      tau_sq = tau_sq,
      sigma_sq = sigma_sq,
      phi = phi
    ),
    spatial = TRUE,
    family = "gaussian"
  )
  class(fit) <- "cos_fit"

  rec <- cos_recover_B(
    fit = fit,
    burn_in = 1,
    n_samples = 1,
    seed = 1,
    verbose = FALSE
  )

  C_B <- make_C_B_tapered_from_pairs_Rcall(C_B_pairs, phi = phi, n_threads = 1L)
  Z <- cbind(X_B, diag(3))
  D_h_inv <- chol2inv(chol(D_h))
  V_beta_inv <- chol2inv(chol(V_beta))
  Q_alpha <- crossprod(Z, D_h_inv %*% Z) / tau_sq
  Q_alpha[1:2, 1:2] <- Q_alpha[1:2, 1:2] + V_beta_inv
  Q_alpha[3:5, 3:5] <- Q_alpha[3:5, 3:5] + chol2inv(chol(C_B)) / sigma_sq
  b_alpha <- as.numeric(crossprod(Z, D_h_inv %*% y_B) / tau_sq)
  b_alpha[1:2] <- b_alpha[1:2] + as.numeric(V_beta_inv %*% mu_beta)
  alpha_mean <- as.numeric(solve(Q_alpha, b_alpha))

  expect_equal(as.numeric(rec$beta_mean_samples[1, ]), alpha_mean[1:2], tolerance = 1e-12)
  expect_equal(as.numeric(rec$omega_B_mean_samples[1, ]), alpha_mean[3:5], tolerance = 1e-12)
})

test_that("Gaussian prediction conditional means and variances match direct algebra", {
  plot_start <- c(1L, 3L, 6L, 7L)
  h <- c(0.6, 0.4, 0.25, 0.5, 0.25, 1)
  x <- c(0, 1, 1, 0, 2, 2)
  y <- c(0, 0, 0, 1, 1, 1)
  gamma <- 2.5
  phi <- 0.35
  sigma_sq <- 1.3
  beta <- c(intercept = 0.4, x = -0.2)
  omega_B <- c(0.3, -0.1, 0.2)

  H_BA_comp <- list(plot_start = plot_start, h = h, x = x, y = y)
  C_B_pairs <- prep_C_B_tapered_pairs_Rcall(
    plot_start = plot_start,
    h = h,
    x = x,
    y = y,
    gamma = gamma,
    taper_code = 1L,
    n_threads = 1L
  )
  C_B <- make_C_B_tapered_from_pairs_Rcall(C_B_pairs, phi = phi, n_threads = 1L)
  R_C_B <- chol(C_B)
  C_B_inv_omega_B <- backsolve(R_C_B, forwardsolve(t(R_C_B), omega_B))

  fit <- list(
    prep = list(
      X = matrix(0, 1, 2, dimnames = list(NULL, names(beta))),
      C_B_pairs = C_B_pairs,
      H_BA_comp = H_BA_comp,
      gamma = gamma,
      taper_code = 1L,
      n_threads = 1L
    ),
    spatial = TRUE,
    family = "gaussian"
  )
  class(fit) <- "cos_fit"

  rec <- list(
    fit = fit,
    theta_samples = data.frame(
      chain = 1L,
      iter = 1L,
      tau_sq = 0.5,
      sigma_sq = sigma_sq,
      phi = phi
    ),
    beta_samples = matrix(beta, nrow = 1L, dimnames = list(NULL, names(beta))),
    omega_B_samples = matrix(omega_B, nrow = 1L),
    spatial = TRUE
  )
  class(rec) <- "cos_recovery_B"

  pred_coords <- matrix(c(
    0.5, 0.2,
    1.5, 0.8
  ), ncol = 2, byrow = TRUE)
  X_pred <- cbind(intercept = 1, x = c(0.1, 0.9))

  C_pred_B <- make_C_pred_B_tapered_Rcall(
    pred_coords = pred_coords,
    H_BA_comp = H_BA_comp,
    phi = phi,
    gamma = gamma,
    taper_code = 1L,
    n_threads = 1L
  )
  omega_pred_mean <- as.numeric(C_pred_B %*% C_B_inv_omega_B)
  eta_pred_mean <- as.numeric(X_pred %*% beta + omega_pred_mean)
  C_B_inv_C_B_pred <- backsolve(R_C_B, forwardsolve(t(R_C_B), t(C_pred_B)))
  fine_cond_var <- sigma_sq *
    pmax(0, 1 - rowSums(C_pred_B * t(C_B_inv_C_B_pred)))

  expect_warning(
    pred_fine <- cos_predict_fine(
      fit = fit,
      rec_B = rec,
      pred_coords = pred_coords,
      X_pred = X_pred,
      spatial_uncertainty = "conditional_mean",
      keep_samples = TRUE,
      verbose = FALSE
    ),
    "conditional_mean"
  )

  expect_equal(pred_fine$summary$omega_mean, omega_pred_mean, tolerance = 1e-12)
  expect_equal(pred_fine$summary$eta_mean, eta_pred_mean, tolerance = 1e-12)
  expect_equal(as.numeric(pred_fine$eta_samples[1, ]), eta_pred_mean, tolerance = 1e-12)
  expect_true(all(fine_cond_var > 0))

  U_blocks <- list(
    U_sf = NULL,
    x_names = names(beta),
    blocks = list(
      list(
        U_id = 1L,
        U_label = "U1",
        n_cell = 2L,
        row_sum = 1,
        h = c(0.7, 0.3),
        coords = pred_coords,
        X_UA = X_pred,
        X_U = as.numeric(crossprod(c(0.7, 0.3), X_pred)),
        cell_id = 1:2
      )
    )
  )
  names(U_blocks$blocks[[1]]$X_U) <- names(beta)
  class(U_blocks) <- "cos_areal_blocks"

  C_UB <- as.numeric(crossprod(U_blocks$blocks[[1]]$h, C_pred_B))
  omega_U_mean <- as.numeric(crossprod(C_UB, C_B_inv_omega_B))
  C_UU <- make_C_B_tapered_Rcall(
    plot_start = c(1L, 3L),
    h = U_blocks$blocks[[1]]$h,
    x = pred_coords[, 1],
    y = pred_coords[, 2],
    phi = phi,
    gamma = gamma,
    taper_code = 1L,
    n_threads = 1L
  )
  C_B_inv_C_BU <- backsolve(R_C_B, forwardsolve(t(R_C_B), C_UB))
  areal_cond_var <- sigma_sq *
    max(0, as.numeric(C_UU[1, 1]) - as.numeric(crossprod(C_UB, C_B_inv_C_BU)))
  eta_U_mean <- as.numeric(crossprod(U_blocks$blocks[[1]]$X_U[names(beta)], beta)) +
    omega_U_mean

  expect_warning(
    pred_areal <- cos_predict_areal(
      fit = fit,
      rec_B = rec,
      U_blocks = U_blocks,
      spatial_uncertainty = "conditional_mean",
      keep_samples = TRUE,
      verbose = FALSE
    ),
    "conditional_mean"
  )

  expect_equal(pred_areal$summary$omega_mean, omega_U_mean, tolerance = 1e-12)
  expect_equal(pred_areal$summary$eta_mean, eta_U_mean, tolerance = 1e-12)
  expect_equal(as.numeric(pred_areal$eta_samples[1, ]), eta_U_mean, tolerance = 1e-12)
  expect_true(areal_cond_var > 0)
})
