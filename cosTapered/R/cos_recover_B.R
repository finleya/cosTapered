cos_recover_B <- function(fit,
                          burn_in = 1,
                          thin = 1L,
                          n_samples = NULL,
                          seed = NULL,
                          verbose = TRUE) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(fit, "cos_fit")) {
    stop("fit must be a cos_fit object from cos_fit().")
  }

  burn_in <- as.integer(burn_in)
  thin <- as.integer(thin)

  if (!is.finite(burn_in) || burn_in < 1L) {
    stop("burn_in must be a positive integer.")
  }
  if (!is.finite(thin) || thin < 1L) {
    stop("thin must be a positive integer.")
  }
  if (!is.null(n_samples)) {
    n_samples <- as.integer(n_samples)
    if (!is.finite(n_samples) || n_samples < 1L) {
      stop("n_samples must be NULL or a positive integer.")
    }
  }
  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------
  # Pull objects into local variables
  # ------------------------------------------------

  prep <- fit$prep
  priors <- fit$priors
  spatial <- isTRUE(fit$spatial)
  if (is.null(fit$spatial)) spatial <- TRUE

  y_B <- as.numeric(prep$y_B)
  X_B <- as.matrix(prep$X_B)
  D_h <- as.matrix(prep$D_h)
  C_B_pairs <- prep$C_B_pairs
  n_threads <- prep$n_threads

  theta_samples <- as.data.frame(fit$theta_samples)

  mu_beta <- as.numeric(priors$beta$mu)
  names(mu_beta) <- colnames(X_B)
  V_beta <- as.matrix(priors$beta$V)

  n_b <- length(y_B)
  p <- ncol(X_B)

  if (!all(c("tau_sq", "sigma_sq", "phi") %in% names(theta_samples))) {
    stop("fit$theta_samples must contain tau_sq, sigma_sq, and phi.")
  }
  if (length(mu_beta) != p || nrow(V_beta) != p || ncol(V_beta) != p) {
    stop("beta prior dimensions are not compatible with prep$X_B.")
  }

  # ------------------------------------------------
  # Select retained theta samples
  # ------------------------------------------------

  theta_keep <- theta_samples[theta_samples$iter >= burn_in, , drop = FALSE]

  if (nrow(theta_keep) == 0L) {
    stop("No theta samples remain after applying burn_in.")
  }

  theta_keep <- theta_keep[order(theta_keep$chain, theta_keep$iter), , drop = FALSE]

  if (thin > 1L) {
    keep_by_chain <- ave(theta_keep$iter, theta_keep$chain, FUN = function(z) {
      seq_along(z)
    })
    theta_keep <- theta_keep[(keep_by_chain - 1L) %% thin == 0L, , drop = FALSE]
  }

  if (nrow(theta_keep) == 0L) {
    stop("No theta samples remain after applying thinning.")
  }

  if (!is.null(n_samples) && nrow(theta_keep) > n_samples) {
    theta_keep <- theta_keep[sample(seq_len(nrow(theta_keep)), n_samples), , drop = FALSE]
    theta_keep <- theta_keep[order(theta_keep$chain, theta_keep$iter), , drop = FALSE]
  }

  n_save <- nrow(theta_keep)

  # ------------------------------------------------
  # Fixed matrices for observed-support recovery
  # ------------------------------------------------

  R_D <- chol(D_h)

  D_h_inv_y <- as.numeric(backsolve(R_D, forwardsolve(t(R_D), y_B)))

  R_V_beta <- chol(V_beta)

  V_beta_inv <- chol2inv(R_V_beta)
  V_beta_inv_mu_beta <- as.numeric(V_beta_inv %*% mu_beta)

  if (spatial) {
    Z <- cbind(X_B, diag(n_b))
    D_h_inv_Z <- backsolve(R_D, forwardsolve(t(R_D), Z))
    omega_ind <- p + seq_len(n_b)
  } else {
    D_h_inv_X_B <- backsolve(R_D, forwardsolve(t(R_D), X_B))
  }

  beta_samples <- matrix(NA_real_, n_save, p)
  omega_B_samples <- matrix(NA_real_, n_save, n_b)
  beta_mean_samples <- matrix(NA_real_, n_save, p)
  omega_B_mean_samples <- matrix(NA_real_, n_save, n_b)

  colnames(beta_samples) <- colnames(X_B)
  colnames(beta_mean_samples) <- colnames(X_B)
  colnames(omega_B_samples) <- paste0("B_", seq_len(n_b))
  colnames(omega_B_mean_samples) <- paste0("B_", seq_len(n_b))

  # ------------------------------------------------
  # Loop over retained theta samples
  # ------------------------------------------------

  for (ii in seq_len(n_save)) {
    if (verbose && (ii == 1L || ii %% 100L == 0L || ii == n_save)) {
      if (spatial) {
        message("Recovering beta and omega_B: ", ii, " of ", n_save)
      } else {
        message("Recovering beta: ", ii, " of ", n_save)
      }
    }

    tau_sq <- theta_keep$tau_sq[ii]

    if (spatial) {
      sigma_sq <- theta_keep$sigma_sq[ii]
      phi <- theta_keep$phi[ii]

      C_B <- make_C_B_tapered_from_pairs_Rcall(
        C_B_pairs = C_B_pairs,
        phi = phi,
        n_threads = n_threads
      )

      R_C_B <- chol(C_B)

      C_B_inv <- chol2inv(R_C_B)

      Q_alpha <- crossprod(Z, D_h_inv_Z) / tau_sq
      Q_alpha[seq_len(p), seq_len(p)] <-
        Q_alpha[seq_len(p), seq_len(p)] + V_beta_inv
      Q_alpha[omega_ind, omega_ind] <-
        Q_alpha[omega_ind, omega_ind] + C_B_inv / sigma_sq
      b_alpha <- as.numeric(crossprod(Z, D_h_inv_y) / tau_sq)
      b_alpha[seq_len(p)] <- b_alpha[seq_len(p)] + V_beta_inv_mu_beta
    } else {
      Q_alpha <- crossprod(X_B, D_h_inv_X_B) / tau_sq + V_beta_inv

      b_alpha <- as.numeric(crossprod(X_B, D_h_inv_y) / tau_sq)
      b_alpha <- b_alpha + V_beta_inv_mu_beta
    }

    R_Q <- chol(Q_alpha)

    alpha_mean <- as.numeric(backsolve(R_Q, forwardsolve(t(R_Q), b_alpha)))
    alpha_draw <- alpha_mean + as.numeric(backsolve(R_Q, stats::rnorm(length(b_alpha))))

    beta_samples[ii, ] <- alpha_draw[seq_len(p)]
    beta_mean_samples[ii, ] <- alpha_mean[seq_len(p)]

    if (spatial) {
      omega_B_samples[ii, ] <- alpha_draw[omega_ind]
      omega_B_mean_samples[ii, ] <- alpha_mean[omega_ind]
    } else {
      omega_B_samples[ii, ] <- rep(0, n_b)
      omega_B_mean_samples[ii, ] <- rep(0, n_b)
    }
  }

  # ------------------------------------------------
  # Observed-support latent response summaries
  # ------------------------------------------------

  eta_B_samples <- beta_samples %*% t(X_B) + omega_B_samples
  colnames(eta_B_samples) <- paste0("B_", seq_len(n_b))

  # ------------------------------------------------
  # Return recovery object
  # ------------------------------------------------

  out <- list(
    fit = fit,
    theta_samples = theta_keep,
    beta_samples = beta_samples,
    omega_B_samples = omega_B_samples,
    eta_B_samples = eta_B_samples,
    beta_mean_samples = beta_mean_samples,
    omega_B_mean_samples = omega_B_mean_samples,
    spatial = spatial,
    call = match.call()
  )

  class(out) <- "cos_recovery_B"
  out
}
