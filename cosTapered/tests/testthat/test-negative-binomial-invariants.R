make_nb_prediction_fixture <- function(n_save = 4L) {
  X_fit <- cbind(intercept = 1, x = c(-1, 0, 1))
  fit <- list(
    family = "negative_binomial",
    size = 5,
    spatial = FALSE,
    prep = list(
      X = X_fit,
      A_coords = cbind(seq_len(nrow(X_fit)), rep(0, nrow(X_fit))),
      n_threads = 1L
    )
  )
  class(fit) <- "cos_fit"

  beta_draw <- c(intercept = log(10), x = 0.25)
  beta_samples <- matrix(rep(beta_draw, each = n_save), nrow = n_save)
  colnames(beta_samples) <- names(beta_draw)

  rec_B <- list(
    theta_samples = data.frame(lp = rep(0, n_save)),
    beta_samples = beta_samples,
    omega_B_samples = matrix(0, n_save, nrow(X_fit))
  )
  class(rec_B) <- "cos_recovery_B"

  list(fit = fit, rec_B = rec_B)
}

make_nb_areal_blocks_fixture <- function() {
  U_blocks <- list(
    U_sf = data.frame(
      U_id = 1:2,
      log_exposure = log(c(1.5, 0.75))
    ),
    blocks = list(
      list(
        U_id = 1L,
        U_label = "1",
        n_cell = 1L,
        row_sum = 1,
        h = 1,
        coords = cbind(0, 0),
        X_UA = cbind(intercept = 1, x = -0.5),
        X_U = c(intercept = 1, x = -0.5),
        cell_id = 1L
      ),
      list(
        U_id = 2L,
        U_label = "2",
        n_cell = 1L,
        row_sum = 1,
        h = 1,
        coords = cbind(1, 0),
        X_UA = cbind(intercept = 1, x = 0.5),
        X_U = c(intercept = 1, x = 0.5),
        cell_id = 2L
      )
    ),
    row_sums = c(1, 1),
    x_names = c("intercept", "x"),
    missing = "renormalize"
  )
  class(U_blocks) <- "cos_areal_blocks"
  U_blocks
}

test_that("negative-binomial fine prediction applies offsets on the mean-count scale", {
  fx <- make_nb_prediction_fixture(n_save = 4L)
  X_pred <- cbind(intercept = 1, x = c(-1, 0.5))
  pred_coords <- cbind(c(0, 1), c(0, 0))

  pred_unit <- cos_predict_fine(
    fit = fx$fit,
    rec_B = fx$rec_B,
    pred_coords = pred_coords,
    X_pred = X_pred,
    offset_pred = 0,
    target = "latent",
    keep_samples = TRUE,
    verbose = FALSE
  )

  pred_double <- cos_predict_fine(
    fit = fx$fit,
    rec_B = fx$rec_B,
    pred_coords = pred_coords,
    X_pred = X_pred,
    offset_pred = log(2),
    target = "latent",
    keep_samples = TRUE,
    verbose = FALSE
  )

  expect_equal(pred_unit$summary$eta_mean + log(2), pred_double$summary$eta_mean)
  expect_equal(2 * pred_unit$summary$mu_mean, pred_double$summary$mu_mean)
  expect_equal(2 * pred_unit$mu_samples, pred_double$mu_samples)
  expect_null(pred_unit$y_samples)
  expect_false("y_mean" %in% names(pred_unit$summary))
})

test_that("negative-binomial observed target adds count draws without changing latent means", {
  fx <- make_nb_prediction_fixture(n_save = 200L)
  X_pred <- cbind(intercept = 1, x = c(0, 1))
  pred_coords <- cbind(c(0, 1), c(0, 0))

  pred_latent <- cos_predict_fine(
    fit = fx$fit,
    rec_B = fx$rec_B,
    pred_coords = pred_coords,
    X_pred = X_pred,
    target = "latent",
    keep_samples = TRUE,
    verbose = FALSE
  )

  set.seed(1)
  pred_observed <- cos_predict_fine(
    fit = fx$fit,
    rec_B = fx$rec_B,
    pred_coords = pred_coords,
    X_pred = X_pred,
    target = "observed",
    keep_samples = TRUE,
    verbose = FALSE
  )

  expect_equal(pred_latent$summary$eta_mean, pred_observed$summary$eta_mean)
  expect_equal(pred_latent$summary$mu_mean, pred_observed$summary$mu_mean)
  expect_true(all(c("y_mean", "y_sd", "y_q025", "y_q975") %in% names(pred_observed$summary)))
  expect_equal(dim(pred_observed$y_samples), dim(pred_observed$mu_samples))
  expect_true(all(pred_observed$y_samples >= 0))
  expect_true(all(abs(pred_observed$y_samples - round(pred_observed$y_samples)) < sqrt(.Machine$double.eps)))
  expect_true(all(pred_latent$summary$mu_sd == 0))
  expect_true(all(pred_observed$summary$y_sd > pred_observed$summary$mu_sd))
})

test_that("negative-binomial areal offset column and vector paths agree", {
  fx <- make_nb_prediction_fixture(n_save = 4L)
  U_blocks <- make_nb_areal_blocks_fixture()

  pred_col <- cos_predict_areal(
    fit = fx$fit,
    rec_B = fx$rec_B,
    U_blocks = U_blocks,
    offset_U = "log_exposure",
    target = "latent",
    keep_samples = TRUE,
    verbose = FALSE
  )

  pred_vec <- cos_predict_areal(
    fit = fx$fit,
    rec_B = fx$rec_B,
    U_blocks = U_blocks,
    offset_U = U_blocks$U_sf$log_exposure,
    target = "latent",
    keep_samples = TRUE,
    verbose = FALSE
  )

  expect_equal(pred_col$summary, pred_vec$summary)
  expect_equal(pred_col$eta_samples, pred_vec$eta_samples)
  expect_equal(pred_col$mu_samples, pred_vec$mu_samples)
})

test_that("prediction guardrails enforce family-specific response options", {
  fx <- make_nb_prediction_fixture(n_save = 4L)
  U_blocks <- make_nb_areal_blocks_fixture()
  X_pred <- cbind(intercept = 1, x = 0)
  pred_coords <- cbind(0, 0)

  expect_error(
    cos_predict_fine(
      fit = fx$fit,
      rec_B = fx$rec_B,
      pred_coords = pred_coords,
      X_pred = X_pred,
      spatial_uncertainty = "joint",
      verbose = FALSE
    ),
    "Polya-Gamma fine prediction currently supports"
  )

  bin_fit <- fx$fit
  bin_fit$family <- "binomial"
  expect_error(
    cos_predict_fine(
      fit = bin_fit,
      rec_B = fx$rec_B,
      pred_coords = pred_coords,
      X_pred = X_pred,
      target = "observed",
      verbose = FALSE
    ),
    "target = 'observed'"
  )
  expect_error(
    cos_predict_areal(
      fit = bin_fit,
      rec_B = fx$rec_B,
      U_blocks = U_blocks,
      target = "observed",
      verbose = FALSE
    ),
    "target = 'observed'"
  )

  gauss_fit <- fx$fit
  gauss_fit$family <- "gaussian"
  expect_error(
    cos_predict_fine(
      fit = gauss_fit,
      rec_B = fx$rec_B,
      pred_coords = pred_coords,
      X_pred = X_pred,
      offset_pred = 0,
      verbose = FALSE
    ),
    "offset_pred is used only"
  )
  expect_error(
    cos_predict_areal(
      fit = gauss_fit,
      rec_B = fx$rec_B,
      U_blocks = U_blocks,
      offset_U = 0,
      verbose = FALSE
    ),
    "offset_U is used only"
  )
})
