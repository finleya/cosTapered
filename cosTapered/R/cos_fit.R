cos_fit <- function(prep,
                    priors,
                    family = c("gaussian", "binomial", "negative_binomial", "negbin", "nb"),
                    trials = NULL,
                    size = NULL,
                    offset = NULL,
                    n_chains = 1L,
                    starting = NULL,
                    tuning = c(log_tau_B_sq = 1, log_sigma_sq = 1, z_phi = 1),
                    n_batch = 1000,
                    batch_length = 25,
                    accept_rate = 0.234,
                    seed = NULL,
                    report = 100,
                    n_threads = NULL,
                    verbose = TRUE) {
  t_start <- proc.time()

  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(prep, "cos_prep")) {
    stop("prep must be a cos_prep object from cos_prepare().")
  }
  if (!inherits(priors, "cos_priors")) {
    stop("priors must be a cos_priors object from cos_default_priors().")
  }
  family <- match.arg(family)
  if (family %in% c("negbin", "nb")) {
    family <- "negative_binomial"
  }

  n_chains <- as.integer(n_chains)
  n_batch <- as.integer(n_batch)
  batch_length <- as.integer(batch_length)
  if (is.null(n_threads)) {
    n_threads <- prep$n_threads
  } else {
    n_threads <- as.integer(n_threads)
    if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L
  }

  if (!is.finite(n_chains) || n_chains < 1L) stop("n_chains must be positive.")
  if (!is.finite(n_batch) || n_batch < 1L) stop("n_batch must be positive.")
  if (!is.finite(batch_length) || batch_length < 1L) stop("batch_length must be positive.")
  if (!is.finite(accept_rate) || accept_rate <= 0 || accept_rate >= 1) {
    stop("accept_rate must be in (0, 1).")
  }

  spatial <- isTRUE(prep$spatial)
  if (is.null(prep$spatial)) spatial <- TRUE

  if (spatial) {
    theta_names <- c("log_tau_B_sq", "log_sigma_sq", "z_phi")
  } else {
    theta_names <- "log_tau_B_sq"
  }

  if (family %in% c("binomial", "negative_binomial")) {
    theta_names <- if (spatial) c("log_sigma_sq", "z_phi") else character(0)
  }

  if (length(theta_names) > 0L) {
    if (is.null(names(tuning)) || !all(theta_names %in% names(tuning))) {
      stop("tuning must be a named vector with names: ", paste(theta_names, collapse = ", "))
    }
    tuning <- as.numeric(tuning[theta_names])
    names(tuning) <- theta_names
    if (any(!is.finite(tuning) | tuning <= 0)) stop("tuning values must be positive and finite.")
  } else {
    tuning <- numeric()
  }

  if (family == "gaussian") {
    if (!is.null(trials)) stop("trials is used only with family = 'binomial'.")
    if (!is.null(size)) stop("size is used only with family = 'negative_binomial'.")
    if (!is.null(offset)) stop("offset is used only with Polya-Gamma response families.")
  }

  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------
  # Pull prep and prior objects into local variables
  # ------------------------------------------------

  y_B <- as.numeric(prep$y_B)
  X_B <- as.matrix(prep$X_B)
  D_h <- as.matrix(prep$D_h)
  C_B_pairs <- prep$C_B_pairs
  n_b <- length(y_B)

  mu_beta <- as.numeric(priors$beta$mu)
  names(mu_beta) <- colnames(X_B)
  V_beta <- as.matrix(priors$beta$V)

  tau_B_shape <- priors$tau_B_sq$shape
  tau_B_scale <- priors$tau_B_sq$scale
  sigma_shape <- priors$sigma_sq$shape
  sigma_scale <- priors$sigma_sq$scale
  phi_lower <- priors$phi$lower
  phi_upper <- priors$phi$upper

  if (spatial && is.null(C_B_pairs)) {
    stop("prep does not contain tapered covariance pairs; rerun cos_prepare() with spatial = TRUE.")
  }
  if (spatial && (!is.finite(phi_lower) || !is.finite(phi_upper) ||
      phi_lower <= 0 || phi_upper <= phi_lower)) {
    stop("Spatial fits require finite prior bounds 0 < phi_lower < phi_upper.")
  }

  d_h_bar <- mean(diag(D_h))
  if (!is.finite(d_h_bar) || d_h_bar <= 0) {
    stop("mean(diag(D_h)) must be positive and finite.")
  }

  X_B_V_beta_X_B <- X_B %*% V_beta %*% t(X_B)
  mu_y <- as.numeric(X_B %*% mu_beta)

  phi_to_z <- function(phi) {
    log((phi - phi_lower) / (phi_upper - phi))
  }

  z_to_phi <- function(z) {
    phi_upper - (phi_upper - phi_lower) / (1 + exp(z))
  }

  if (family %in% c("binomial", "negative_binomial")) {
    family_code <- if (family == "binomial") 1L else 2L
    family_label <- if (family == "binomial") "binomial" else "negative-binomial"

    if (!is.null(offset) && is.character(offset) && length(offset) == 1L) {
      if (!offset %in% names(prep$B_sf)) {
        stop("offset column not found in prep$B_sf.")
      }
      offset <- prep$B_sf[[offset]]
    }
    if (is.null(offset)) {
      offset <- rep(0, n_b)
    }
    if (length(offset) == 1L) {
      offset <- rep(offset, n_b)
    }
    offset <- as.numeric(offset)
    if (length(offset) != n_b || any(!is.finite(offset))) {
      stop("offset must be NULL, a column name in prep$B_sf, or a finite numeric vector of length one or length(prep$y_B).")
    }

    if (family == "binomial") {
      if (!is.null(size)) {
        stop("size is used only with family = 'negative_binomial'.")
      }
      if (is.null(trials)) {
        stop("trials must be supplied for family = 'binomial'.")
      }
      if (is.character(trials) && length(trials) == 1L) {
        if (!trials %in% names(prep$B_sf)) {
          stop("trials column not found in prep$B_sf.")
        }
        trials <- prep$B_sf[[trials]]
      }
      if (length(trials) == 1L) {
        trials <- rep(trials, n_b)
      }
      if (length(trials) != n_b) {
        stop("trials must have length one or length(prep$y_B).")
      }
      trials <- as.integer(trials)
      if (any(!is.finite(trials) | trials < 1L)) {
        stop("trials must contain positive integers.")
      }
      if (any(!is.finite(y_B) | y_B < 0 | y_B != floor(y_B))) {
        stop("Binomial responses in prep$y_B must be nonnegative integers.")
      }
      if (any(y_B > trials)) {
        stop("Binomial responses in prep$y_B cannot exceed trials.")
      }
      nb_size <- NA_real_
    } else {
      if (!is.null(trials)) {
        stop("trials is used only with family = 'binomial'.")
      }
      if (is.null(size)) {
        stop("size must be supplied for family = 'negative_binomial'.")
      }
      size <- as.numeric(size)
      if (length(size) != 1L || !is.finite(size) || size <= 0) {
        stop("size must be a finite positive scalar for family = 'negative_binomial'.")
      }
      nb_size <- size
      trials <- rep(1L, n_b)
      if (any(!is.finite(y_B) | y_B < 0 | y_B != floor(y_B))) {
        stop("Negative-binomial responses in prep$y_B must be nonnegative integers.")
      }
    }

    # ------------------------------------------------
    # Starting values for explicit PG state
    # ------------------------------------------------

    beta_start <- matrix(rep(mu_beta, each = n_chains), nrow = n_chains)
    colnames(beta_start) <- colnames(X_B)

    if (is.null(starting)) {
      if (spatial) {
        if (n_chains == 1L) {
          phi_grid <- exp(mean(log(c(phi_lower, phi_upper))))
        } else {
          phi_grid <- exp(seq(log(phi_lower), log(phi_upper), length.out = n_chains + 2L))
          phi_grid <- phi_grid[-c(1L, length(phi_grid))]
        }
        starting <- data.frame(
          sigma_sq = rep(1, n_chains),
          phi = phi_grid
        )
      } else {
        starting <- data.frame(.chain = seq_len(n_chains))
      }
    } else {
      starting <- as.data.frame(starting)
      if (nrow(starting) == 1L && n_chains > 1L) {
        starting <- starting[rep(1L, n_chains), , drop = FALSE]
      }
      if (nrow(starting) != n_chains) {
        stop("starting must have either one row or n_chains rows.")
      }
      beta_cols <- intersect(colnames(X_B), names(starting))
      if (length(beta_cols) > 0L) {
        if (length(beta_cols) != ncol(X_B)) {
          stop("For Polya-Gamma fits, starting beta values must include all covariates or none.")
        }
        beta_start <- as.matrix(starting[, colnames(X_B), drop = FALSE])
      }
    }

    if (spatial) {
      if (!all(c("sigma_sq", "phi") %in% names(starting))) {
        stop("Spatial Polya-Gamma starting values must include sigma_sq and phi.")
      }
      if (any(!is.finite(starting$sigma_sq) | starting$sigma_sq <= 0)) {
        stop("starting$sigma_sq must be positive and finite.")
      }
      if (any(!is.finite(starting$phi) | starting$phi <= phi_lower | starting$phi >= phi_upper)) {
        stop("starting$phi must be finite and inside (phi_lower, phi_upper).")
      }
    }

    # ------------------------------------------------
    # Run PG chains
    # ------------------------------------------------

    chain_fits <- vector("list", n_chains)
    theta_z_list <- vector("list", n_chains)
    theta_list <- vector("list", n_chains)
    lp_list <- vector("list", n_chains)
    beta_list <- vector("list", n_chains)
    omega_list <- vector("list", n_chains)
    eta_list <- vector("list", n_chains)

    for (ch in seq_len(n_chains)) {
      if (verbose) {
        cat("\nStarting ", family_label, " PG chain ", ch, " of ", n_chains, "\n", sep = "")
      }

      chain_fits[[ch]] <- binomial_pg_sampler_Rcall(
        y_B = y_B,
        trials = trials,
        offset = offset,
        family_code = family_code,
        nb_size = nb_size,
        X_B = X_B,
        C_B_pairs = C_B_pairs,
        spatial = spatial,
        beta_mu = mu_beta,
        V_beta = V_beta,
        starting_beta = beta_start[ch, ],
        starting_sigma_sq = if (spatial) starting$sigma_sq[ch] else NA_real_,
        starting_phi = if (spatial) starting$phi[ch] else NA_real_,
        tuning = tuning,
        n_batch = n_batch,
        batch_length = batch_length,
        accept_rate = accept_rate,
        c0 = 1,
        c1 = 0.8,
        report = report,
        verbose = verbose,
        sigma_shape = sigma_shape,
        sigma_scale = sigma_scale,
        phi_lower = phi_lower,
        phi_upper = phi_upper,
        n_threads = n_threads
      )

      theta_z <- as.data.frame(chain_fits[[ch]]$theta_z_samples)
      if (spatial) names(theta_z) <- theta_names
      theta_z$chain <- ch
      theta_z$iter <- seq_len(nrow(theta_z))
      theta_z_list[[ch]] <- theta_z

      theta <- theta_z[, c("chain", "iter"), drop = FALSE]
      if (spatial) {
        theta$sigma_sq <- exp(theta_z$log_sigma_sq)
        theta$phi <- z_to_phi(theta_z$z_phi)
        theta$eff_range <- 3 / theta$phi
      } else {
        theta$sigma_sq <- NA_real_
        theta$phi <- NA_real_
        theta$eff_range <- NA_real_
      }
      theta$lp <- as.numeric(chain_fits[[ch]]$p.lp.samples)
      theta_list[[ch]] <- theta

      lp_list[[ch]] <- data.frame(
        chain = ch,
        iter = seq_along(chain_fits[[ch]]$p.lp.samples),
        lp = as.numeric(chain_fits[[ch]]$p.lp.samples)
      )

      beta_dat <- as.data.frame(chain_fits[[ch]]$beta_samples)
      beta_dat$chain <- ch
      beta_dat$iter <- seq_len(nrow(beta_dat))
      beta_list[[ch]] <- beta_dat

      omega_dat <- as.data.frame(chain_fits[[ch]]$omega_B_samples)
      omega_dat$chain <- ch
      omega_dat$iter <- seq_len(nrow(omega_dat))
      omega_list[[ch]] <- omega_dat

      eta_dat <- as.data.frame(chain_fits[[ch]]$eta_B_samples)
      eta_dat$chain <- ch
      eta_dat$iter <- seq_len(nrow(eta_dat))
      eta_list[[ch]] <- eta_dat
    }

    theta_z_samples <- do.call(rbind, theta_z_list)
    rownames(theta_z_samples) <- NULL
    theta_samples <- do.call(rbind, theta_list)
    rownames(theta_samples) <- NULL
    lp_samples <- do.call(rbind, lp_list)
    rownames(lp_samples) <- NULL

    beta_samples <- do.call(rbind, beta_list)
    beta_samples <- beta_samples[order(beta_samples$chain, beta_samples$iter), , drop = FALSE]
    rownames(beta_samples) <- NULL

    omega_B_samples <- do.call(rbind, omega_list)
    omega_B_samples <- omega_B_samples[order(omega_B_samples$chain, omega_B_samples$iter), , drop = FALSE]
    rownames(omega_B_samples) <- NULL

    eta_B_samples <- do.call(rbind, eta_list)
    eta_B_samples <- eta_B_samples[order(eta_B_samples$chain, eta_B_samples$iter), , drop = FALSE]
    rownames(eta_B_samples) <- NULL

    dt <- proc.time() - t_start
    timing <- c(
      user = unname(dt[["user.self"]]),
      system = unname(dt[["sys.self"]]),
      elapsed = unname(dt[["elapsed"]])
    )

    out <- list(
      prep = prep,
      priors = priors,
      family = family,
      trials = trials,
      size = if (family == "negative_binomial") nb_size else NULL,
      offset = offset,
      theta_samples = theta_samples,
      theta_z_samples = theta_z_samples,
      lp_samples = lp_samples,
      beta_samples = beta_samples,
      omega_B_samples = omega_B_samples,
      eta_B_samples = eta_B_samples,
      chain_fits = chain_fits,
      starting = starting,
      tuning = tuning,
      sampler = list(
        n_chains = n_chains,
        n_batch = n_batch,
        batch_length = batch_length,
        accept_rate = accept_rate
      ),
      spatial = spatial,
      d_h_bar = d_h_bar,
      timing = timing,
      call = match.call()
    )

    class(out) <- "cos_fit"
    return(out)
  }

  # ------------------------------------------------
  # Starting values on natural scale
  # ------------------------------------------------

  if (is.null(starting)) {
    y_var <- stats::var(y_B)
    if (!is.finite(y_var) || y_var <= 0) y_var <- 1

    if (spatial) {
      if (n_chains == 1L) {
        phi_grid <- exp(mean(log(c(phi_lower, phi_upper))))
      } else {
        # Use interior points on the log-phi scale. The endpoints are excluded
        # because phi must be strictly inside (phi_lower, phi_upper).
        phi_grid <- exp(seq(log(phi_lower), log(phi_upper), length.out = n_chains + 2L))
        phi_grid <- phi_grid[-c(1L, length(phi_grid))]
      }
    }

    if (spatial) {
      starting <- data.frame(
        tau_B_sq = y_var * exp(stats::rnorm(n_chains, log(0.25), 0.5)),
        sigma_sq = y_var * exp(stats::rnorm(n_chains, log(0.50), 0.5)),
        phi = phi_grid
      )
    } else {
      starting <- data.frame(
        tau_B_sq = y_var * exp(stats::rnorm(n_chains, log(0.25), 0.5))
      )
    }
  } else {
    if (spatial) {
      req <- c("tau_B_sq", "sigma_sq", "phi")
    } else {
      req <- "tau_B_sq"
    }

    if (is.numeric(starting) && is.null(dim(starting))) {
      if (is.null(names(starting)) || !all(req %in% names(starting))) {
        stop("starting must include named values: tau_B_sq, sigma_sq, phi.")
      }
      starting <- as.data.frame(as.list(starting[req]))
    } else {
      starting <- as.data.frame(starting)
      if (!all(req %in% names(starting))) {
        stop("starting must include columns: tau_B_sq, sigma_sq, phi.")
      }
      starting <- starting[, req, drop = FALSE]
    }

    if (nrow(starting) == 1L && n_chains > 1L) {
      starting <- starting[rep(1L, n_chains), , drop = FALSE]
    }
    if (nrow(starting) != n_chains) {
      stop("starting must have either one row or n_chains rows.")
    }
  }

  if (any(!is.finite(starting$tau_B_sq) | starting$tau_B_sq <= 0)) {
    stop("starting$tau_B_sq must be positive and finite.")
  }
  if (spatial && any(!is.finite(starting$sigma_sq) | starting$sigma_sq <= 0)) {
    stop("starting$sigma_sq must be positive and finite.")
  }
  if (spatial && any(!is.finite(starting$phi) | starting$phi <= phi_lower | starting$phi >= phi_upper)) {
    stop("starting$phi must be finite and inside (phi_lower, phi_upper).")
  }

  # ------------------------------------------------
  # Collapsed log target for covariance parameters
  # ------------------------------------------------

  theta_log_target <- function(theta_z) {
    names(theta_z) <- theta_names

    tau_B_sq <- exp(theta_z["log_tau_B_sq"])

    if (tau_B_sq <= 0 || !is.finite(tau_B_sq)) {
      return(-Inf)
    }

    tau_sq <- tau_B_sq / d_h_bar

    if (spatial) {
      sigma_sq <- exp(theta_z["log_sigma_sq"])

      phi <- z_to_phi(theta_z["z_phi"])

      if (sigma_sq <= 0 ||
          phi <= phi_lower || phi >= phi_upper ||
          !is.finite(sigma_sq) || !is.finite(phi)) {
        return(-Inf)
      }

      C_B <- make_C_B_tapered_from_pairs_Rcall(
        C_B_pairs = C_B_pairs,
        phi = phi,
        n_threads = n_threads
      )

      V_y <- sigma_sq * C_B + X_B_V_beta_X_B + tau_sq * D_h
    } else {
      V_y <- X_B_V_beta_X_B + tau_sq * D_h
    }
    R <- chol(V_y)

    r <- y_B - mu_y
    V_inv_r <- backsolve(R, forwardsolve(t(R), r))
    log_det_V <- 2 * sum(log(diag(R)))

    log_lik <- -0.5 * (
      n_b * log(2 * pi) +
        log_det_V +
        sum(r * V_inv_r)
    )

    log_prior_tau_B <-
      tau_B_shape * log(tau_B_scale) -
      lgamma(tau_B_shape) -
      (tau_B_shape + 1) * log(tau_B_sq) -
      tau_B_scale / tau_B_sq

    if (spatial) {
      log_prior_sigma <-
        sigma_shape * log(sigma_scale) -
        lgamma(sigma_shape) -
        (sigma_shape + 1) * log(sigma_sq) -
        sigma_scale / sigma_sq

      log_prior_phi <- -log(phi_upper - phi_lower)

      log_jac <- log(tau_B_sq) + log(sigma_sq) +
        log(phi - phi_lower) + log(phi_upper - phi) -
        log(phi_upper - phi_lower)

      out <- log_lik + log_prior_tau_B + log_prior_sigma + log_prior_phi + log_jac
    } else {
      log_jac <- log(tau_B_sq)
      out <- log_lik + log_prior_tau_B + log_jac
    }
    if (!is.finite(out)) return(-Inf)
    unname(out)
  }

  # ------------------------------------------------
  # Run chains
  # ------------------------------------------------

  chain_fits <- vector("list", n_chains)
  theta_z_list <- vector("list", n_chains)
  lp_list <- vector("list", n_chains)

  for (ch in seq_len(n_chains)) {
    if (spatial) {
      theta_start_z <- c(
        log_tau_B_sq = log(starting$tau_B_sq[ch]),
        log_sigma_sq = log(starting$sigma_sq[ch]),
        z_phi = phi_to_z(starting$phi[ch])
      )
    } else {
      theta_start_z <- c(
        log_tau_B_sq = log(starting$tau_B_sq[ch])
      )
    }

    if (verbose) {
      cat("\nStarting chain ", ch, " of ", n_chains,
          "; log posterior = ", signif(theta_log_target(theta_start_z), 6), "\n", sep = "")
    }

    chain_fits[[ch]] <- gaussian_metrop_sampler_Rcall(
      y_B = y_B,
      mu_y = mu_y,
      X_B_V_beta_X_B = X_B_V_beta_X_B,
      D_h = D_h,
      C_B_pairs = C_B_pairs,
      spatial = spatial,
      starting = theta_start_z,
      tuning = tuning,
      n_batch = n_batch,
      batch_length = batch_length,
      accept_rate = accept_rate,
      c0 = 1,
      c1 = 0.8,
      report = report,
      verbose = verbose,
      d_h_bar = d_h_bar,
      tau_B_shape = tau_B_shape,
      tau_B_scale = tau_B_scale,
      sigma_shape = sigma_shape,
      sigma_scale = sigma_scale,
      phi_lower = phi_lower,
      phi_upper = phi_upper,
      n_threads = n_threads
    )

    theta_z <- as.data.frame(chain_fits[[ch]]$p.theta.samples)
    names(theta_z) <- theta_names
    theta_z$chain <- ch
    theta_z$iter <- seq_len(nrow(theta_z))
    theta_z_list[[ch]] <- theta_z

    lp_dat <- data.frame(
      chain = ch,
      iter = seq_along(chain_fits[[ch]]$p.lp.samples),
      lp = as.numeric(chain_fits[[ch]]$p.lp.samples)
    )
    lp_list[[ch]] <- lp_dat
  }

  theta_z_samples <- do.call(rbind, theta_z_list)
  rownames(theta_z_samples) <- NULL

  lp_samples <- do.call(rbind, lp_list)
  rownames(lp_samples) <- NULL

  # ------------------------------------------------
  # Transform and collect samples
  # ------------------------------------------------

  theta_samples <- theta_z_samples[, c("chain", "iter"), drop = FALSE]
  theta_samples$tau_B_sq <- exp(theta_z_samples$log_tau_B_sq)
  theta_samples$tau_sq <- theta_samples$tau_B_sq / d_h_bar

  if (spatial) {
    theta_samples$sigma_sq <- exp(theta_z_samples$log_sigma_sq)
    theta_samples$phi <- z_to_phi(theta_z_samples$z_phi)
    theta_samples$eff_range <- 3 / theta_samples$phi
  } else {
    theta_samples$sigma_sq <- NA_real_
    theta_samples$phi <- NA_real_
    theta_samples$eff_range <- NA_real_
  }
  theta_samples$lp <- lp_samples$lp

  # ------------------------------------------------
  # Return fit object
  # ------------------------------------------------

  dt <- proc.time() - t_start
  timing <- c(
    user = unname(dt[["user.self"]]),
    system = unname(dt[["sys.self"]]),
    elapsed = unname(dt[["elapsed"]])
  )

  out <- list(
    prep = prep,
    priors = priors,
    family = family,
    theta_samples = theta_samples,
    theta_z_samples = theta_z_samples,
    lp_samples = lp_samples,
    chain_fits = chain_fits,
    starting = starting,
    tuning = tuning,
    sampler = list(
      n_chains = n_chains,
      n_batch = n_batch,
      batch_length = batch_length,
      accept_rate = accept_rate
    ),
    spatial = spatial,
    d_h_bar = d_h_bar,
    timing = timing,
    call = match.call()
  )

  class(out) <- "cos_fit"
  out
}

cos_summary_theta <- function(fit, burn_in = 1,
                              pars = c("tau_B_sq", "tau_sq", "sigma_sq", "phi", "eff_range"),
                              probs = c(0.5, 0.025, 0.975),
                              digits = 3) {
  if (!inherits(fit, "cos_fit")) stop("fit must be a cos_fit object.")
  summary(fit, burn_in = burn_in, pars = pars, probs = probs, digits = digits)$theta_summary
}

cos_plot_trace <- function(fit, burn_in = 1,
                           pars = c("tau_B_sq", "sigma_sq", "phi", "eff_range")) {
  if (!inherits(fit, "cos_fit")) stop("fit must be a cos_fit object.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")

  theta <- fit$theta_samples
  theta <- theta[theta$iter >= burn_in, , drop = FALSE]
  pars <- intersect(pars, names(theta))
  pars <- pars[vapply(pars, function(nm) any(is.finite(theta[[nm]])), logical(1))]
  if (length(pars) == 0L) stop("No requested parameters have finite retained samples.")

  dat <- data.frame(
    chain = theta$chain,
    iter = theta$iter,
    theta[, pars, drop = FALSE]
  )

  dat <- stats::reshape(
    dat,
    varying = pars,
    v.names = "value",
    timevar = "parameter",
    times = pars,
    direction = "long"
  )

  ggplot2::ggplot(dat, ggplot2::aes(iter, value, group = chain)) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.8) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "MCMC iteration", y = NULL, title = "Collapsed theta chains") +
    ggplot2::theme_bw()
}

plot.cos_fit <- function(x, burn_in = 1,
                         pars = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
                         ...) {
  cos_plot_trace(fit = x, burn_in = burn_in, pars = pars)
}
