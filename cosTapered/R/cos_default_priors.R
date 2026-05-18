cos_default_priors <- function(prep,
                               beta_mu = NULL,
                               beta_sd = 10,
                               tau_B_shape = 2,
                               tau_B_scale = 50^2,
                               sigma_shape = 2,
                               sigma_scale = 10,
                               phi_lower,
                               phi_upper) {
  if (!inherits(prep, "cos_prep")) {
    stop("prep must be a cos_prep object from cos_prepare().")
  }
  spatial <- isTRUE(prep$spatial)
  if (is.null(prep$spatial)) spatial <- TRUE

  if (spatial && (missing(phi_lower) || missing(phi_upper))) {
    stop("phi_lower and phi_upper must be supplied.")
  }
  if (spatial && (!is.finite(phi_lower) || !is.finite(phi_upper) ||
      phi_lower <= 0 || phi_upper <= phi_lower)) {
    stop("Require 0 < phi_lower < phi_upper.")
  }
  if (!spatial) {
    phi_lower <- NA_real_
    phi_upper <- NA_real_
  }

  x_names <- colnames(prep$X_B)
  p <- length(x_names)

  if (is.null(beta_mu)) {
    beta_mu <- rep(0, p)
    names(beta_mu) <- x_names
  } else {
    if (is.null(names(beta_mu))) {
      if (length(beta_mu) != p) stop("beta_mu must have length ncol(prep$X_B).")
      names(beta_mu) <- x_names
    }
    if (!all(x_names %in% names(beta_mu))) {
      stop("beta_mu must contain names: ", paste(x_names, collapse = ", "))
    }
    beta_mu <- as.numeric(beta_mu[x_names])
    names(beta_mu) <- x_names
  }

  if (length(beta_sd) == 1) {
    beta_sd <- rep(beta_sd, p)
    names(beta_sd) <- x_names
  } else {
    if (is.null(names(beta_sd))) {
      if (length(beta_sd) != p) stop("beta_sd must have length ncol(prep$X_B).")
      names(beta_sd) <- x_names
    }
    if (!all(x_names %in% names(beta_sd))) {
      stop("beta_sd must contain names: ", paste(x_names, collapse = ", "))
    }
    beta_sd <- as.numeric(beta_sd[x_names])
    names(beta_sd) <- x_names
  }

  if (any(!is.finite(beta_mu))) stop("beta_mu contains non-finite values.")
  if (any(!is.finite(beta_sd) | beta_sd <= 0)) stop("beta_sd must be positive and finite.")

  V_beta <- diag(beta_sd^2, p)
  dimnames(V_beta) <- list(x_names, x_names)

  if (!is.finite(tau_B_shape) || !is.finite(tau_B_scale) ||
      tau_B_shape <= 0 || tau_B_scale <= 0) {
    stop("tau_B_shape and tau_B_scale must be positive and finite.")
  }
  if (!is.finite(sigma_shape) || !is.finite(sigma_scale) ||
      sigma_shape <= 0 || sigma_scale <= 0) {
    stop("sigma_shape and sigma_scale must be positive and finite.")
  }

  out <- list(
    beta = list(mu = beta_mu, V = V_beta),
    tau_B_sq = list(shape = tau_B_shape, scale = tau_B_scale),
    sigma_sq = list(shape = sigma_shape, scale = sigma_scale),
    phi = list(lower = phi_lower, upper = phi_upper)
  )

  class(out) <- "cos_priors"
  out
}
