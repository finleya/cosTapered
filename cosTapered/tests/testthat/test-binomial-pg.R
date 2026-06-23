make_binomial_test_prep <- function(spatial = FALSE) {
  X_B <- cbind(intercept = 1, x = c(-1.2, -0.4, 0.2, 0.9, 1.4))
  beta <- c(intercept = -0.2, x = 0.9)
  eta <- as.numeric(X_B %*% beta)
  trials <- c(5L, 7L, 6L, 8L, 5L)
  y_B <- c(1, 2, 3, 6, 4)

  prep <- list(
    y_B = y_B,
    X_B = X_B,
    D_h = diag(length(y_B)),
    spatial = spatial,
    n_threads = 1L,
    B_sf = data.frame(y_B = y_B, trials = trials),
    eta_true = eta
  )

  if (spatial) {
    coords <- cbind(c(0, 1, 2, 3, 4), c(0, 0.2, 0.1, 0.3, 0))
    dist_B <- as.matrix(stats::dist(coords))
    pair_l <- integer()
    pair_k <- integer()
    pair_d <- numeric()
    pair_wtap <- numeric()
    for (l in seq_len(nrow(dist_B))) {
      for (k in seq_len(l)) {
        pair_l <- c(pair_l, l)
        pair_k <- c(pair_k, k)
        pair_d <- c(pair_d, dist_B[l, k])
        pair_wtap <- c(pair_wtap, 1)
      }
    }
    prep$C_B_pairs <- list(
      pair_l = pair_l,
      pair_k = pair_k,
      pair_d = pair_d,
      pair_wtap = pair_wtap,
      n_b = length(y_B)
    )
    prep$H_BA_comp <- list(
      plot_start = seq_len(length(y_B) + 1L),
      h = rep(1, length(y_B)),
      x = coords[, 1],
      y = coords[, 2]
    )
    prep$gamma <- 10
    prep$taper_code <- 1L
  } else {
    prep$C_B_pairs <- NULL
  }

  class(prep) <- "cos_prep"
  prep
}

test_that("binomial PG nonspatial fit and recovery run", {
  prep <- make_binomial_test_prep(spatial = FALSE)
  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 0, x = 0),
    beta_sd = c(intercept = 2, x = 2)
  )

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    family = "binomial",
    trials = prep$B_sf$trials,
    n_batch = 3,
    batch_length = 4,
    seed = 11,
    verbose = FALSE
  )

  expect_s3_class(fit, "cos_fit")
  expect_equal(fit$family, "binomial")
  expect_equal(nrow(fit$beta_samples), 12)
  expect_equal(nrow(fit$omega_B_samples), 12)
  expect_true(all(is.finite(as.matrix(fit$beta_samples[, c("intercept", "x")]))))
  expect_true(all(fit$omega_B_samples[, paste0("B_", seq_along(prep$y_B))] == 0))
  expect_s3_class(summary(fit), "summary.cos_fit")

  rec <- cos_recover_B(fit, burn_in = 2, verbose = FALSE)
  expect_s3_class(rec, "cos_recovery_B")
  expect_equal(nrow(rec$beta_samples), 11)
  expect_equal(ncol(rec$eta_B_samples), length(prep$y_B))
})

test_that("binomial PG spatial fit runs", {
  prep <- make_binomial_test_prep(spatial = TRUE)
  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 0, x = 0),
    beta_sd = c(intercept = 2, x = 2),
    sigma_shape = 2,
    sigma_scale = 1,
    phi_lower = 0.1,
    phi_upper = 3
  )

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    family = "binomial",
    trials = "trials",
    n_batch = 3,
    batch_length = 3,
    starting = data.frame(sigma_sq = 0.8, phi = 0.8),
    tuning = c(log_sigma_sq = 0.15, z_phi = 0.15),
    seed = 12,
    verbose = FALSE
  )

  expect_s3_class(fit, "cos_fit")
  expect_true(all(c("sigma_sq", "phi", "eff_range", "lp") %in% names(fit$theta_samples)))
  expect_true(all(is.finite(fit$theta_samples$sigma_sq)))
  expect_true(all(is.finite(fit$theta_samples$phi)))
  expect_true(all(is.finite(as.matrix(fit$omega_B_samples[, paste0("B_", seq_along(prep$y_B))]))))
  expect_s3_class(summary(fit), "summary.cos_fit")

  rec <- cos_recover_B(fit, burn_in = 1, n_samples = 4, seed = 1, verbose = FALSE)
  expect_equal(nrow(rec$eta_B_samples), 4)

  pred <- cos_predict_fine(
    fit = fit,
    rec_B = rec,
    pred_coords = cbind(c(0.25, 1.5, 3.5), c(0.1, 0.2, 0.2)),
    X_pred = cbind(intercept = 1, x = c(-1, 0.1, 1)),
    spatial_uncertainty = "marginal",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred, "cos_prediction_fine")
  expect_true(all(c("p_mean", "p_sd", "conditional_omega_sd_mean") %in% names(pred$summary)))
  expect_true(all(pred$summary$p_mean >= 0 & pred$summary$p_mean <= 1))
  expect_equal(dim(pred$p_samples), c(4, 3))
  expect_false(pred$uncertainty$joint_samples)

  U_blocks <- list(
    U_sf = data.frame(U_id = 1L),
    blocks = list(list(
      U_id = 1L,
      U_label = "1",
      n_cell = 2L,
      row_sum = 1,
      h = c(0.5, 0.5),
      coords = cbind(c(0.25, 1.5), c(0.1, 0.2)),
      X_UA = cbind(intercept = c(1, 1), x = c(-1, 0.1)),
      X_U = c(intercept = 1, x = -0.45),
      cell_id = 1:2
    )),
    row_sums = 1,
    x_names = c("intercept", "x"),
    missing = "renormalize"
  )
  class(U_blocks) <- "cos_areal_blocks"

  pred_U <- cos_predict_areal(
    fit = fit,
    rec_B = rec,
    U_blocks = U_blocks,
    spatial_uncertainty = "marginal",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred_U, "cos_prediction_areal")
  expect_true(all(c("p_mean", "p_sd", "conditional_omega_sd_mean") %in% names(pred_U$summary)))
  expect_true(pred_U$summary$p_mean >= 0 && pred_U$summary$p_mean <= 1)
  expect_equal(dim(pred_U$p_samples), c(4, 1))
  expect_false(pred_U$uncertainty$joint_samples)
})

test_that("binomial PG validates trials and responses", {
  prep <- make_binomial_test_prep(spatial = FALSE)
  priors <- cos_default_priors(prep = prep, beta_mu = c(intercept = 0, x = 0))

  expect_error(
    cos_fit(prep, priors, family = "binomial", n_batch = 1, batch_length = 1, verbose = FALSE),
    "trials must be supplied"
  )
  expect_error(
    cos_fit(prep, priors, family = "binomial", trials = c(1, 1), n_batch = 1, batch_length = 1, verbose = FALSE),
    "trials must have length"
  )

  prep_bad <- prep
  prep_bad$y_B[1] <- 9
  expect_error(
    cos_fit(prep_bad, priors, family = "binomial", trials = prep$B_sf$trials, n_batch = 1, batch_length = 1, verbose = FALSE),
    "cannot exceed trials"
  )
})
