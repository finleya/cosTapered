test_that("fine prediction samples semidefinite conditional covariance", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("exactextractr")

  chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
  B_file <- system.file("extdata", "example_B.gpkg", package = "cosTapered")
  truth_file <- system.file("extdata", "example_truth.rds", package = "cosTapered")
  if (chm_file == "") {
    chm_file <- testthat::test_path("..", "..", "inst", "extdata", "example_chm.tif")
  }
  if (B_file == "") {
    B_file <- testthat::test_path("..", "..", "inst", "extdata", "example_B.gpkg")
  }
  if (truth_file == "") {
    truth_file <- testthat::test_path("..", "..", "inst", "extdata", "example_truth.rds")
  }
  skip_if(!file.exists(chm_file) || !file.exists(B_file) || !file.exists(truth_file),
          "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"

  intercept <- chm
  intercept[] <- ifelse(is.na(raster::getValues(chm)), NA_real_, 1)
  names(intercept) <- "intercept"

  X_rast <- raster::stack(intercept, chm)
  names(X_rast) <- c("intercept", "chm")

  B_sf <- sf::st_read(B_file, quiet = TRUE)
  truth <- readRDS(truth_file)

  prep <- cos_prepare(
    X_rast = X_rast,
    B_sf = B_sf,
    response_col = "y_B",
    gamma = truth$gamma,
    taper_code = truth$taper_code,
    verbose = FALSE
  )

  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 40, chm = 1),
    beta_sd = c(intercept = 20, chm = 2),
    tau_B_shape = 2,
    tau_B_scale = truth$tau_B_sq_true,
    sigma_shape = 2,
    sigma_scale = truth$sigma_sq_true,
    phi_lower = 3 / 600,
    phi_upper = 3 / 80
  )

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    starting = data.frame(
      tau_B_sq = truth$tau_B_sq_true,
      sigma_sq = truth$sigma_sq_true,
      phi = truth$phi_true
    ),
    tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
    n_batch = 2,
    batch_length = 2,
    seed = 1,
    verbose = FALSE
  )

  rec <- cos_recover_B(
    fit = fit,
    burn_in = 1,
    n_samples = 1,
    seed = 2,
    verbose = FALSE
  )

  rec$theta_samples$phi <- 3 / 600
  pred_ind <- seq_len(750)

  expect_warning(
    pred <- cos_predict_fine(
      fit = fit,
      rec_B = rec,
      pred_coords = prep$A_coords[pred_ind, , drop = FALSE],
      X_pred = prep$X[pred_ind, , drop = FALSE],
      spatial_uncertainty = "joint",
      keep_samples = TRUE,
      verbose = FALSE
    ),
    "positive semidefinite"
  )

  expect_s3_class(pred, "cos_prediction_fine")
  expect_equal(ncol(pred$eta_samples), length(pred_ind))
  expect_true(all(is.finite(pred$summary$eta_mean)))
})
