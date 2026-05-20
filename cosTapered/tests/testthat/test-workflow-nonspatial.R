test_that("non-spatial workflow and S3 methods run on example data", {
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

  prep <- cos_prepare(
    X_rast = X_rast,
    B_sf = B_sf,
    response_col = "y_B",
    spatial = FALSE,
    verbose = FALSE
  )
  expect_s3_class(prep, "cos_prep")
  expect_named(prep$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(prep$timing)))
  expect_true(prep$timing[["elapsed"]] >= 0)
  expect_s3_class(summary(prep), "summary.cos_prep")

  priors <- cos_default_priors(
    prep = prep,
    beta_mu = c(intercept = 40, chm = 1),
    beta_sd = c(intercept = 20, chm = 2),
    tau_B_shape = 2,
    tau_B_scale = 8
  )
  expect_s3_class(priors, "cos_priors")

  fit <- cos_fit(
    prep = prep,
    priors = priors,
    n_batch = 2,
    batch_length = 2,
    tuning = c(log_tau_B_sq = 0.25),
    seed = 1,
    verbose = FALSE
  )
  expect_s3_class(fit, "cos_fit")
  expect_named(fit$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(fit$timing)))
  expect_true(fit$timing[["elapsed"]] >= 0)
  expect_s3_class(summary(fit), "summary.cos_fit")
  expect_s3_class(plot(fit), "ggplot")
  expect_true(is.matrix(cos_summary_theta(fit)))
  expect_s3_class(cos_plot_trace(fit), "ggplot")

  rec <- cos_recover_B(
    fit = fit,
    burn_in = 1,
    n_samples = 2,
    seed = 2,
    verbose = FALSE
  )
  expect_s3_class(rec, "cos_recovery_B")
  expect_named(rec$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(rec$timing)))
  expect_true(rec$timing[["elapsed"]] >= 0)
  expect_s3_class(summary(rec), "summary.cos_recovery_B")

  pred <- cos_predict_fine(
    fit = fit,
    rec_B = rec,
    verbose = FALSE
  )
  expect_s3_class(pred, "cos_prediction_fine")
  expect_named(pred$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(pred$timing)))
  expect_true(pred$timing[["elapsed"]] >= 0)
  expect_s3_class(summary(pred), "summary.cos_prediction_fine")
  expect_equal(nrow(pred$summary), nrow(prep$X))

  pred_y <- cos_predict_fine(
    fit = fit,
    rec_B = rec,
    target = "observed",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred_y, "cos_prediction_fine")
  expect_true(all(c("y_mean", "y_sd", "y_q025", "y_q975") %in% names(pred_y$summary)))
  expect_equal(dim(pred_y$y_samples), dim(pred_y$eta_samples))
  expect_true(mean(pred_y$summary$y_sd, na.rm = TRUE) > 0)

  B_sub <- B_sf[seq_len(2), ]
  B_sub$.test_id <- seq_len(nrow(B_sub))
  U_blocks <- cos_make_areal_blocks(
    U_sf = B_sub,
    X_rast = X_rast,
    id_col = ".test_id",
    verbose = FALSE
  )
  expect_named(U_blocks$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(U_blocks$timing)))
  expect_true(U_blocks$timing[["elapsed"]] >= 0)
  pred_U <- cos_predict_areal(
    fit = fit,
    rec_B = rec,
    U_blocks = U_blocks,
    target = "observed",
    keep_samples = TRUE,
    verbose = FALSE
  )
  expect_s3_class(pred_U, "cos_prediction_areal")
  expect_named(pred_U$timing, c("user", "system", "elapsed"))
  expect_true(all(is.finite(pred_U$timing)))
  expect_true(pred_U$timing[["elapsed"]] >= 0)
  expect_true(all(c("y_mean", "y_sd", "y_q025", "y_q975") %in% names(pred_U$summary)))
  expect_equal(dim(pred_U$y_samples), dim(pred_U$eta_samples))
  expect_equal(ncol(pred_U$y_samples), length(U_blocks$blocks))
})
