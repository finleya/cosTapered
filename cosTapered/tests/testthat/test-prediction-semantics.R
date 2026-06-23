test_that("prediction objects record uncertainty semantics", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("exactextractr")
  skip_if_not_installed("ggplot2")

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
  rec <- cos_recover_B(
    fit = fit,
    burn_in = 1,
    n_samples = 2,
    seed = 2,
    verbose = FALSE
  )

  pred_fine <- cos_predict_fine(fit = fit, rec_B = rec, verbose = FALSE)
  expect_equal(pred_fine$uncertainty$spatial_effect, "none_nonspatial")
  expect_false(pred_fine$uncertainty$joint_samples)
  expect_match(summary(pred_fine)$uncertainty$note, "No residual spatial effect")

  U_blocks <- cos_make_areal_blocks(
    U_sf = B_sf[seq_len(2), ],
    X_rast = X_rast,
    verbose = FALSE
  )
  pred_areal <- cos_predict_areal(
    fit = fit,
    rec_B = rec,
    U_blocks = U_blocks,
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_equal(pred_areal$uncertainty$spatial_effect, "none_nonspatial")
  expect_false(pred_areal$uncertainty$joint_samples)
  expect_match(summary(pred_areal)$uncertainty$note, "No residual spatial effect")
})

test_that("cos_simulate_data requires explicit intercept raster", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")

  chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
  forest_file <- system.file("extdata", "example_forest.gpkg", package = "cosTapered")
  if (chm_file == "") {
    chm_file <- testthat::test_path("..", "..", "inst", "extdata", "example_chm.tif")
  }
  if (forest_file == "") {
    forest_file <- testthat::test_path("..", "..", "inst", "extdata", "example_forest.gpkg")
  }
  skip_if(!file.exists(chm_file) || !file.exists(forest_file), "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"
  forest_sf <- sf::st_read(forest_file, quiet = TRUE)

  expect_error(
    cos_simulate_data(
      X_rast = chm,
      forest_sf = forest_sf,
      n_B = 3,
      plot_radius = 10,
      gamma = 100,
      beta = c(intercept = 1, chm = 2),
      tau_B_sq = 1,
      sigma_sq = 1,
      phi = 3 / 100,
      verbose = FALSE
    ),
    "Missing: intercept"
  )
})
