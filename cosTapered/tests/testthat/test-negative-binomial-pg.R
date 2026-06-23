make_nb_test_prep <- function(spatial = FALSE, n = 6L) {
  x <- seq(-1, 1, length.out = n)
  X_B <- cbind(intercept = 1, x = x)
  offset <- log(seq(0.8, 1.3, length.out = n))
  beta <- c(intercept = 0.2, x = -0.35)
  eta <- as.numeric(offset + X_B %*% beta)
  size <- 4
  y_B <- c(1, 2, 1, 3, 2, 4)[seq_len(n)]

  prep <- list(
    y_B = y_B,
    X_B = X_B,
    X = X_B,
    A_coords = cbind(seq_len(n), rep(0, n)),
    D_h = diag(n),
    spatial = spatial,
    n_threads = 1L,
    B_sf = data.frame(y_B = y_B, offset = offset),
    eta_true = eta,
    size_true = size
  )

  if (spatial) {
    coords <- cbind(seq_len(n) - 1, rep(0, n))
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
      n_b = n
    )
    prep$H_BA_comp <- list(
      plot_start = seq_len(n + 1L),
      h = rep(1, n),
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

test_that("negative-binomial PG nonspatial fit and recovery run", {
  prep <- make_nb_test_prep(spatial = FALSE)
  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 0, x = 0),
    beta_sd = c(intercept = 3, x = 3)
  )

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    family = "nb",
    size = prep$size_true,
    offset = "offset",
    n_batch = 4,
    batch_length = 5,
    seed = 21,
    verbose = FALSE
  )

  expect_s3_class(fit, "cos_fit")
  expect_equal(fit$family, "negative_binomial")
  expect_equal(fit$size, prep$size_true)
  expect_equal(fit$offset, prep$B_sf$offset)
  expect_equal(nrow(fit$beta_samples), 20)
  expect_true(all(is.finite(as.matrix(fit$beta_samples[, c("intercept", "x")]))))
  expect_true(all(fit$omega_B_samples[, paste0("B_", seq_along(prep$y_B))] == 0))
  expect_s3_class(summary(fit), "summary.cos_fit")

  rec <- cos_recover_B(fit, burn_in = 2, verbose = FALSE)
  expect_s3_class(rec, "cos_recovery_B")
  expect_equal(nrow(rec$beta_samples), 19)
  expect_equal(ncol(rec$eta_B_samples), length(prep$y_B))
})

test_that("negative-binomial PG spatial prediction returns mean and count summaries", {
  prep <- make_nb_test_prep(spatial = TRUE)
  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 0, x = 0),
    beta_sd = c(intercept = 3, x = 3),
    sigma_shape = 2,
    sigma_scale = 1,
    phi_lower = 0.1,
    phi_upper = 3
  )

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    family = "negative_binomial",
    size = 5,
    offset = prep$B_sf$offset,
    n_batch = 3,
    batch_length = 4,
    starting = data.frame(sigma_sq = 0.8, phi = 0.8),
    tuning = c(log_sigma_sq = 0.15, z_phi = 0.15),
    seed = 22,
    verbose = FALSE
  )

  rec <- cos_recover_B(fit, burn_in = 1, n_samples = 4, seed = 1, verbose = FALSE)

  pred <- cos_predict_fine(
    fit = fit,
    rec_B = rec,
    pred_coords = cbind(c(0.25, 1.5, 3.5), c(0, 0, 0)),
    X_pred = cbind(intercept = 1, x = c(-1, 0.1, 1)),
    offset_pred = log(c(1, 1.2, 0.9)),
    spatial_uncertainty = "marginal",
    target = "observed",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred, "cos_prediction_fine")
  expect_true(all(c("mu_mean", "mu_sd", "y_mean", "y_sd", "conditional_omega_sd_mean") %in% names(pred$summary)))
  expect_true(all(pred$summary$mu_mean > 0))
  expect_equal(dim(pred$mu_samples), c(4, 3))
  expect_equal(dim(pred$y_samples), c(4, 3))
  expect_true(all(pred$y_samples >= 0))
  expect_true(all(abs(pred$y_samples - round(pred$y_samples)) < sqrt(.Machine$double.eps)))

  U_blocks <- list(
    U_sf = data.frame(U_id = 1L, log_exposure = log(2)),
    blocks = list(list(
      U_id = 1L,
      U_label = "1",
      n_cell = 2L,
      row_sum = 1,
      h = c(0.5, 0.5),
      coords = cbind(c(0.25, 1.5), c(0, 0)),
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
    offset_U = "log_exposure",
    spatial_uncertainty = "marginal",
    target = "observed",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred_U, "cos_prediction_areal")
  expect_true(all(c("mu_mean", "mu_sd", "y_mean", "y_sd", "conditional_omega_sd_mean") %in% names(pred_U$summary)))
  expect_true(pred_U$summary$mu_mean > 0)
  expect_equal(dim(pred_U$mu_samples), c(4, 1))
  expect_equal(dim(pred_U$y_samples), c(4, 1))
})

test_that("negative-binomial PG validates size, counts, and incompatible trials", {
  prep <- make_nb_test_prep(spatial = FALSE)
  priors <- cos_default_priors(prep = prep, beta_mu = c(intercept = 0, x = 0))

  expect_error(
    cos_fit(prep, priors, family = "negative_binomial", n_batch = 1, batch_length = 1, verbose = FALSE),
    "size must be supplied"
  )
  expect_error(
    cos_fit(prep, priors, family = "negative_binomial", size = 0, n_batch = 1, batch_length = 1, verbose = FALSE),
    "size must be a finite positive scalar"
  )
  expect_error(
    cos_fit(prep, priors, family = "negative_binomial", size = 2, trials = 1, n_batch = 1, batch_length = 1, verbose = FALSE),
    "trials is used only"
  )

  prep_bad <- prep
  prep_bad$y_B[1] <- 1.5
  expect_error(
    cos_fit(prep_bad, priors, family = "negative_binomial", size = 2, n_batch = 1, batch_length = 1, verbose = FALSE),
    "nonnegative integers"
  )
})
