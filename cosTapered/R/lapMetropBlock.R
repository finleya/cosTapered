lapMetropBlock <- function(ltd, starting, tuning = NULL, batch = 100, batch.length = 50,
                           accept.rate = 0.234, c0 = 1, c1 = 0.8,
                           report = 10, verbose = TRUE) {
  theta <- as.numeric(starting)
  theta_names <- names(starting)
  d <- length(theta)

  if (is.null(tuning)) {
    tuning <- rep(1, d)
  }
  tuning <- as.numeric(tuning)

  n_save <- batch * batch.length
  samples <- matrix(NA_real_, n_save, d)
  colnames(samples) <- theta_names
  lp_samples <- rep(NA_real_, n_save)

  accept <- matrix(0L, batch, batch.length)
  batch_accept_rate <- rep(NA_real_, batch)
  sigma_sq_m <- 2.4^2 / d
  Sigma0 <- diag(tuning^2, d)

  log_post <- function(x) {
    out <- ltd(x)
    if (!is.finite(out)) return(-Inf)
    out
  }

  lp <- log_post(theta)
  if (!is.finite(lp)) stop("Starting value has non-finite log target density.")

  make_chol <- function(S) {
    chol(S)
  }

  save_i <- 0L

  for (b in seq_len(batch)) {
    batch_samples <- matrix(NA_real_, batch.length, d)
    R_prop <- make_chol(sigma_sq_m * Sigma0)

    for (j in seq_len(batch.length)) {
      theta_prop <- as.numeric(theta + drop(t(R_prop) %*% rnorm(d)))
      lp_prop <- log_post(theta_prop)

      if (log(runif(1)) < lp_prop - lp) {
        theta <- theta_prop
        lp <- lp_prop
        accept[b, j] <- 1L
      }

      save_i <- save_i + 1L
      samples[save_i, ] <- theta
      lp_samples[save_i] <- lp
      batch_samples[j, ] <- theta
    }

    rhat <- mean(accept[b, ])
    batch_accept_rate[b] <- rhat

    gamma1 <- 1 / (b + 1)^c1
    gamma2 <- c0 * gamma1

    sigma_sq_m <- exp(log(sigma_sq_m) + gamma2 * (rhat - accept.rate))

    if (batch.length > 1) {
      Sigma_hat <- stats::cov(batch_samples)
      if (all(is.finite(Sigma_hat))) {
        Sigma0 <- Sigma0 + gamma1 * (Sigma_hat - Sigma0)
      }
    }

    if (verbose && (b == 1 || b %% report == 0)) {
      cat("batch", b,
          "accept =", round(rhat, 3),
          "sigma_sq_m =", signif(sigma_sq_m, 4),
          "lp =", signif(lp, 6), "\n")
    }
  }

  list(
    p.theta.samples = samples,
    p.lp.samples = lp_samples,
    accept = accept,
    batch.accept.rate = batch_accept_rate,
    sigma_sq_m = sigma_sq_m,
    Sigma0 = Sigma0,
    proposal.cov = sigma_sq_m * Sigma0,
    call = match.call()
  )
}
