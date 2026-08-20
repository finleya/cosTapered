test_that("cos_recover_B, cos_predict_areal, and cos_predict_fine default to gaussian when fit$family is NULL", {
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
  skip_if(!file.exists(chm_file) || !file.exists(B_file), "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"
  intercept <- chm
  intercept[] <- ifelse(is.na(raster::getValues(chm)), NA_real_, 1)
  names(intercept) <- "intercept"
  X_rast <- raster::stack(intercept, chm)
  names(X_rast) <- c("intercept", "chm")
  B_sf <- sf::st_read(B_file, quiet = TRUE)

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

  # Simulate an older or hand-built fit object that never recorded $family.
  fit$family <- NULL

  rec <- cos_recover_B(fit = fit, burn_in = 1, n_samples = 2, seed = 2, verbose = FALSE)
  expect_s3_class(rec, "cos_recovery_B")

  pred_fine <- cos_predict_fine(fit = fit, rec_B = rec, verbose = FALSE)
  expect_s3_class(pred_fine, "cos_prediction_fine")

  U_blocks <- cos_make_areal_blocks(
    U_sf = B_sf[seq_len(2), ],
    X_rast = X_rast,
    verbose = FALSE
  )
  pred_areal <- cos_predict_areal(fit = fit, rec_B = rec, U_blocks = U_blocks, verbose = FALSE)
  expect_s3_class(pred_areal, "cos_prediction_areal")
})
