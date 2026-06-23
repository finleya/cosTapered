test_that("C++ Gaussian Metropolis sampler matches R reference path", {
  y_B <- c(1.1, -0.4, 0.7)
  X_B <- cbind(intercept = 1, x = c(-1, 0.5, 1.2))
  mu_beta <- c(intercept = 0.2, x = 0.1)
  V_beta <- diag(c(2, 1.5))
  D_h <- diag(c(0.7, 1.1, 0.9))
  X_B_V_beta_X_B <- X_B %*% V_beta %*% t(X_B)
  mu_y <- as.numeric(X_B %*% mu_beta)

  d_h_bar <- mean(diag(D_h))
  tau_B_shape <- 2.1
  tau_B_scale <- 1.3
  starting <- c(log_tau_B_sq = log(0.8))
  tuning <- c(log_tau_B_sq = 0.35)

  ltd <- function(theta_z) {
    tau_B_sq <- exp(theta_z[1])
    tau_sq <- tau_B_sq / d_h_bar
    V_y <- X_B_V_beta_X_B + tau_sq * D_h
    R <- chol(V_y)
    r <- y_B - mu_y
    V_inv_r <- backsolve(R, forwardsolve(t(R), r))

    log_lik <- -0.5 * (
      length(y_B) * log(2 * pi) +
        2 * sum(log(diag(R))) +
        sum(r * V_inv_r)
    )
    log_prior_tau_B <-
      tau_B_shape * log(tau_B_scale) -
      lgamma(tau_B_shape) -
      (tau_B_shape + 1) * log(tau_B_sq) -
      tau_B_scale / tau_B_sq

    log_lik + log_prior_tau_B + log(tau_B_sq)
  }

  set.seed(44)
  r_fit <- lapMetropBlock(
    ltd = ltd,
    starting = starting,
    tuning = tuning,
    batch = 5,
    batch.length = 4,
    accept.rate = 0.234,
    c0 = 1,
    c1 = 0.8,
    report = 100,
    verbose = FALSE
  )

  set.seed(44)
  cpp_fit <- gaussian_metrop_sampler_Rcall(
    y_B = y_B,
    mu_y = mu_y,
    X_B_V_beta_X_B = X_B_V_beta_X_B,
    D_h = D_h,
    spatial = FALSE,
    starting = starting,
    tuning = tuning,
    n_batch = 5,
    batch_length = 4,
    accept_rate = 0.234,
    c0 = 1,
    c1 = 0.8,
    report = 100,
    verbose = FALSE,
    d_h_bar = d_h_bar,
    tau_B_shape = tau_B_shape,
    tau_B_scale = tau_B_scale
  )

  expect_equal(cpp_fit$p.theta.samples, r_fit$p.theta.samples, tolerance = 1e-12)
  expect_equal(cpp_fit$p.lp.samples, r_fit$p.lp.samples, tolerance = 1e-12)
  expect_equal(cpp_fit$accept, r_fit$accept)
  expect_equal(cpp_fit$batch.accept.rate, r_fit$batch.accept.rate)
  expect_equal(unname(cpp_fit$Sigma0), unname(r_fit$Sigma0), tolerance = 1e-12)
  expect_equal(unname(cpp_fit$proposal.cov), unname(r_fit$proposal.cov), tolerance = 1e-12)
})

test_that("C++ Gaussian Metropolis sampler matches R spatial reference path", {
  y_B <- c(0.2, -0.8, 1.3)
  X_B <- cbind(intercept = 1, x = c(-0.5, 0.3, 1.0))
  mu_beta <- c(intercept = 0.1, x = -0.2)
  V_beta <- diag(c(1.4, 0.8))
  D_h <- diag(c(0.9, 0.8, 1.2))
  X_B_V_beta_X_B <- X_B %*% V_beta %*% t(X_B)
  mu_y <- as.numeric(X_B %*% mu_beta)

  dist_B <- as.matrix(stats::dist(cbind(c(0, 1, 2.5), c(0, 0.4, 0.9))))
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
  C_B_pairs <- list(
    pair_l = pair_l,
    pair_k = pair_k,
    pair_d = pair_d,
    pair_wtap = pair_wtap,
    n_b = length(y_B)
  )

  d_h_bar <- mean(diag(D_h))
  tau_B_shape <- 2.2
  tau_B_scale <- 1.4
  sigma_shape <- 2.5
  sigma_scale <- 1.1
  phi_lower <- 0.05
  phi_upper <- 2.5

  phi_to_z <- function(phi) {
    log((phi - phi_lower) / (phi_upper - phi))
  }
  z_to_phi <- function(z) {
    phi_upper - (phi_upper - phi_lower) / (1 + exp(z))
  }

  starting <- c(
    log_tau_B_sq = log(0.7),
    log_sigma_sq = log(0.9),
    z_phi = phi_to_z(0.6)
  )
  tuning <- c(log_tau_B_sq = 0.25, log_sigma_sq = 0.3, z_phi = 0.2)

  ltd <- function(theta_z) {
    tau_B_sq <- exp(theta_z[1])
    sigma_sq <- exp(theta_z[2])
    phi <- z_to_phi(theta_z[3])
    tau_sq <- tau_B_sq / d_h_bar

    C_B <- exp(-phi * dist_B)
    V_y <- sigma_sq * C_B + X_B_V_beta_X_B + tau_sq * D_h
    R <- chol(V_y)
    r <- y_B - mu_y
    V_inv_r <- backsolve(R, forwardsolve(t(R), r))

    log_lik <- -0.5 * (
      length(y_B) * log(2 * pi) +
        2 * sum(log(diag(R))) +
        sum(r * V_inv_r)
    )
    log_prior_tau_B <-
      tau_B_shape * log(tau_B_scale) -
      lgamma(tau_B_shape) -
      (tau_B_shape + 1) * log(tau_B_sq) -
      tau_B_scale / tau_B_sq
    log_prior_sigma <-
      sigma_shape * log(sigma_scale) -
      lgamma(sigma_shape) -
      (sigma_shape + 1) * log(sigma_sq) -
      sigma_scale / sigma_sq
    log_prior_phi <- -log(phi_upper - phi_lower)
    log_jac <- log(tau_B_sq) + log(sigma_sq) +
      log(phi - phi_lower) + log(phi_upper - phi) -
      log(phi_upper - phi_lower)

    log_lik + log_prior_tau_B + log_prior_sigma + log_prior_phi + log_jac
  }

  set.seed(45)
  r_fit <- lapMetropBlock(
    ltd = ltd,
    starting = starting,
    tuning = tuning,
    batch = 4,
    batch.length = 3,
    accept.rate = 0.234,
    c0 = 1,
    c1 = 0.8,
    report = 100,
    verbose = FALSE
  )

  set.seed(45)
  cpp_fit <- gaussian_metrop_sampler_Rcall(
    y_B = y_B,
    mu_y = mu_y,
    X_B_V_beta_X_B = X_B_V_beta_X_B,
    D_h = D_h,
    C_B_pairs = C_B_pairs,
    spatial = TRUE,
    starting = starting,
    tuning = tuning,
    n_batch = 4,
    batch_length = 3,
    accept_rate = 0.234,
    c0 = 1,
    c1 = 0.8,
    report = 100,
    verbose = FALSE,
    d_h_bar = d_h_bar,
    tau_B_shape = tau_B_shape,
    tau_B_scale = tau_B_scale,
    sigma_shape = sigma_shape,
    sigma_scale = sigma_scale,
    phi_lower = phi_lower,
    phi_upper = phi_upper
  )

  expect_equal(cpp_fit$p.theta.samples, r_fit$p.theta.samples, tolerance = 1e-12)
  expect_equal(cpp_fit$p.lp.samples, r_fit$p.lp.samples, tolerance = 1e-12)
  expect_equal(cpp_fit$accept, r_fit$accept)
  expect_equal(cpp_fit$batch.accept.rate, r_fit$batch.accept.rate)
})
