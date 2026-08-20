test_that("observed-support CV returns scores", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("exactextractr")

  chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
  B_file <- system.file("extdata", "example_B.gpkg", package = "cosTapered")
  if (chm_file == "") {
    chm_file <- testthat::test_path("..", "..", "inst", "extdata", "example_chm.tif")
  }
  if (B_file == "") {
    B_file <- testthat::test_path("..", "..", "inst", "extdata", "example_B.gpkg")
  }
  skip_if(chm_file == "" || B_file == "", "example data not installed")
  skip_if(!file.exists(chm_file) || !file.exists(B_file), "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"

  intercept <- chm
  intercept[] <- ifelse(is.na(raster::getValues(chm)), NA_real_, 1)
  names(intercept) <- "intercept"

  X_rast <- raster::stack(intercept, chm)
  names(X_rast) <- c("intercept", "chm")

  B_sf <- sf::st_read(B_file, quiet = TRUE)
  B_sf$eta_test <- B_sf$y_B

  prep <- cos_prepare(
    X_rast = X_rast,
    B_sf = B_sf,
    response_col = "y_B",
    spatial = FALSE,
    verbose = FALSE
  )
  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 40, chm = 1),
    beta_sd = c(intercept = 20, chm = 2),
    tau_B_shape = 2,
    tau_B_scale = 8
  )
  fit <- cos_fit(
    prep = prep,
    priors = priors,
    n_batch = 2,
    batch_length = 2,
    tuning = c(log_tau_B_sq = 0.25),
    seed = 1,
    verbose = FALSE
  )

  fold_id <- rep(1:2, length.out = nrow(B_sf))
  cv <- cos_cv_observed(
    fit = fit,
    X_rast = X_rast,
    k = 2,
    fold_id = fold_id,
    fit_args = list(n_batch = 2, batch_length = 2, report = 99),
    recover_args = list(burn_in = 1, n_samples = 2),
    seed = 1,
    verbose = FALSE
  )

  expect_s3_class(cv, "cos_cv_observed")
  expect_named(cv$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(cv$timing)))
  expect_true(cv$timing[["elapsed"]] >= 0)
  expect_s3_class(summary(cv), "summary.cos_cv_observed")
  expect_equal(nrow(cv$predictions), nrow(B_sf))
  expect_true(all(c("RMSPE", "CRPS", "coverage_95") %in% names(cv$summary)))
  expect_true(is.finite(cv$summary$RMSPE))
  expect_true(is.finite(cv$summary$CRPS))

  cv_both <- cos_cv_observed(
    fit = fit,
    X_rast = X_rast,
    B_sf = B_sf,
    response_col = "y_B",
    target = c("observed", "latent"),
    latent_col = "eta_test",
    k = 2,
    fold_id = fold_id,
    fit_args = list(n_batch = 2, batch_length = 2, report = 99),
    recover_args = list(burn_in = 1, n_samples = 2),
    seed = 1,
    verbose = FALSE
  )

  expect_s3_class(cv_both, "cos_cv_observed")
  expect_equal(sort(unique(cv_both$predictions$target)), c("latent", "observed"))
  expect_equal(nrow(cv_both$predictions), 2L * nrow(B_sf))
  expect_equal(nrow(cv_both$summary), 2L)
  expect_true(all(c("target", "truth", "observed", "latent") %in% names(cv_both$predictions)))
  expect_error(
    cos_cv_observed(
      fit = fit,
      X_rast = X_rast,
      B_sf = B_sf,
      target = "latent",
      k = 2,
      fold_id = fold_id,
      fit_args = list(n_batch = 1, batch_length = 1),
      recover_args = list(burn_in = 1, n_samples = 1),
      verbose = FALSE
    ),
    "latent_col"
  )
})

test_that("cos_cv_observed rejects target = 'observed' for non-Gaussian families", {
  X_B <- cbind(intercept = 1, x = c(-1.2, -0.4, 0.2, 0.9, 1.4, 1.8))
  beta <- c(intercept = -0.2, x = 0.9)
  eta <- as.numeric(X_B %*% beta)
  trials <- c(5L, 7L, 6L, 8L, 5L, 6L)
  y_B <- c(1, 2, 3, 6, 4, 5)

  prep <- list(
    y_B = y_B,
    X_B = X_B,
    D_h = diag(length(y_B)),
    spatial = FALSE,
    n_threads = 1L,
    B_sf = data.frame(y_B = y_B, trials = trials),
    eta_true = eta,
    C_B_pairs = NULL
  )
  class(prep) <- "cos_prep"

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
    n_batch = 2,
    batch_length = 2,
    seed = 21,
    verbose = FALSE
  )

  B_sf <- data.frame(y_B = y_B)

  expect_error(
    cos_cv_observed(
      fit = fit,
      X_rast = raster::raster(matrix(1)),
      B_sf = B_sf,
      target = "observed",
      k = 2,
      verbose = FALSE
    ),
    "cos_cv_observed\\(\\) is currently only implemented for Gaussian fits"
  )
})
