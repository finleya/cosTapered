grid_covariance_fixture <- function(nx = 37L, ny = 29L) {
  grid <- list(res = c(3, 3), origin = c(680001.5, 5960001.5))
  index <- as.matrix(expand.grid(seq_len(nx) - 1, seq_len(ny) - 1))
  # Holes, disconnected cells and fractional polygon boundaries.
  index <- index[(index[, 1] + 2 * index[, 2]) %% 7 != 0, , drop = FALSE]
  coords <- sweep(index * 3, 2, grid$origin, "+")
  h <- (seq_len(nrow(index)) %% 11 + 1)
  h <- 0.83 * h / sum(h) # Retain unnormalized weights, as with missing='drop'.
  bcoords <- coords[c(1, 2, 200, 201, 700, 701), , drop = FALSE]
  H <- list(plot_start = c(1L, 3L, 5L, 7L),
            h = c(.25, .75, .6, .4, .9, .1), x = bcoords[, 1], y = bcoords[, 2])
  blocks <- list(blocks = list(list(h = h, coords = coords)), grid = grid)
  list(U = blocks, H = H)
}

test_that("cached grid covariance agrees with direct pair sums", {
  f <- grid_covariance_fixture()
  b <- f$U$blocks[[1]]
  for (gamma in c(1, 12, 18.7, 150)) {
    for (taper in c(1L, 2L)) {
      cache <- cos_prepare_areal_grid_covariance(
        f$U, f$H, gamma, taper, TRUE, n_threads = 2L,
        min_cells = 1L, min_cross_cells = 1L)[[1]]
      expect_false(is.null(cache$self))
      expect_false(is.null(cache$cross))
      for (phi in c(.003, .04, .7)) {
        got <- cos_eval_areal_grid_covariance(cache, phi)
        self <- make_C_B_tapered_Rcall(c(1L, length(b$h) + 1L), b$h,
          b$coords[, 1], b$coords[, 2], phi, gamma, taper, 1L)[1, 1]
        cross <- as.numeric(crossprod(b$h, make_C_pred_B_tapered_Rcall(
          b$coords, f$H, phi, gamma, taper, 1L)))
        expect_equal(got$C_UU, self, tolerance = 1e-11)
        expect_equal(got$C_UB, cross, tolerance = 1e-11)
      }
    }
  }
})

test_that("grid covariance has bounded fallbacks and handles duplicate cells", {
  f <- grid_covariance_fixture()
  args <- list(H_BA_comp = f$H, gamma = 18, taper_code = 1L,
               need_self = TRUE, min_cells = 1L, min_cross_cells = 1L)
  prep <- function(u, ...) do.call(cos_prepare_areal_grid_covariance,
                                  c(list(U_blocks = u), args, list(...)))
  expect_null(prep(within(f$U, rm(grid)))[[1]])
  bad <- f$U
  bad$grid$res[2] <- 6
  expect_null(prep(bad)[[1]])
  bad <- f$U
  bad$blocks[[1]]$coords[1, 1] <- bad$blocks[[1]]$coords[1, 1] + .1
  expect_null(prep(bad)[[1]])
  expect_null(prep(f$U, cache_bytes = 1)[[1]])
  # An off-grid observed support prevents cross caching only.
  args$H_BA_comp$x[1] <- args$H_BA_comp$x[1] + .1
  out <- prep(f$U)[[1]]
  expect_null(out$cross)
  expect_false(is.null(out$self))
  index <- matrix(c(0, 0, 0, 0, 1, 0), ncol = 2, byrow = TRUE)
  weights <- cos_areal_grid_self_weights(index, c(.2, .3, .5), 4)
  expect_equal(weights, c(.5, .5, 0, 0), tolerance = 1e-12)
  expect_equal(cos_areal_grid_self_weights(index, c(.2, .3, .5), 4, 8L),
               weights, tolerance = 1e-12)
  # Widely separated cells need no dense raster allocation.
  sparse <- matrix(c(0L, 0L, 10000000L, 10000000L), ncol = 2, byrow = TRUE)
  expect_equal(cos_areal_grid_self_weights(sparse, c(.3, .7), 4, 2L),
               c(.58, 0, 0, 0), tolerance = 1e-12)
})

test_that("areal prediction preserves samples for all supported families", {
  f <- grid_covariance_fixture(112L, 108L)
  b <- f$U$blocks[[1]]
  b$U_id <- 1L; b$U_label <- "forest"; b$n_cell <- length(b$h)
  b$row_sum <- sum(b$h); b$X_U <- c(intercept = sum(b$h))
  U <- f$U; U$blocks <- list(b); U$x_names <- "intercept"
  class(U) <- "cos_areal_blocks"
  pairs <- prep_C_B_tapered_pairs_Rcall(f$H$plot_start, f$H$h,
                                       f$H$x, f$H$y, 18, 1L, 1L)
  fit <- structure(list(prep = list(X = matrix(1, 1, 1,
    dimnames = list(NULL, "intercept")), H_BA_comp = f$H,
    C_B_pairs = pairs, gamma = 18, taper_code = 1L, n_threads = 1L),
    spatial = TRUE, size = 4), class = "cos_fit")
  rec <- structure(list(theta_samples = data.frame(sigma_sq = c(1, 2, 1),
    phi = c(.03, .1, .05), tau_sq = c(.1, .2, .1)),
    beta_samples = matrix(c(.2, .3, .1), ncol = 1),
    omega_B_samples = matrix(seq_len(9)/10, nrow = 3)), class = "cos_recovery_B")
  legacy <- U; legacy$grid <- NULL
  for (family in c("gaussian", "binomial", "negative_binomial")) {
    fit$family <- family
    for (uncertainty in c("conditional_mean", "marginal")) {
      for (target in if (family == "binomial") "latent" else c("latent", "observed")) {
        set.seed(81)
        fast <- suppressWarnings(cos_predict_areal(fit, rec, U,
          spatial_uncertainty = uncertainty, target = target,
          keep_samples = TRUE, verbose = FALSE))
        set.seed(81)
        slow <- suppressWarnings(cos_predict_areal(fit, rec, legacy,
          spatial_uncertainty = uncertainty, target = target,
          keep_samples = TRUE, verbose = FALSE))
        expect_equal(fast$summary, slow$summary, tolerance = 1e-10)
        for (field in c("eta_samples", "omega_samples", "mu_samples", "p_samples", "y_samples"))
          expect_equal(fast[[field]], slow[[field]], tolerance = 1e-10)
      }
    }
  }
})


test_that("direct self-cache matches unordered pair accumulation across threads", {
  index <- matrix(c(2, -3, -1, 4, 2, -3, 0, 0, 5, 1, -1, -2),
                  ncol = 2, byrow = TRUE)
  h <- c(.1, .2, .15, .05, .3, .2)
  n_bin <- 26L
  ref <- numeric(n_bin)
  for (i in seq_along(h)) {
    for (j in seq_along(h)) {
      d2 <- sum((index[i, ] - index[j, ])^2)
      if (d2 < n_bin) ref[d2 + 1L] <- ref[d2 + 1L] + h[i] * h[j]
    }
  }
  order <- c(6, 2, 4, 1, 5, 3)
  for (threads in c(1L, 2L, 8L)) {
    expect_equal(cos_areal_grid_self_weights(index, h, n_bin, threads),
                 ref, tolerance = 1e-12)
    expect_equal(cos_areal_grid_self_weights(index[order, ], h[order], n_bin, threads),
                 ref, tolerance = 1e-12)
  }
})
