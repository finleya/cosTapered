make_C_B_tapered_Rcall <- function(plot_start, h, x, y, phi, gamma,
                                   taper_code = 1L, n_threads = 1L) {
  .Call("make_C_B_tapered",
        as.integer(plot_start),
        as.numeric(h),
        as.numeric(x),
        as.numeric(y),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

prep_C_B_tapered_pairs_Rcall <- function(plot_start, h, x, y, gamma,
                                          taper_code = 1L, n_threads = 1L) {
  .Call("prep_C_B_tapered_pairs",
        as.integer(plot_start),
        as.numeric(h),
        as.numeric(x),
        as.numeric(y),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_B_tapered_from_pairs_Rcall <- function(C_B_pairs, phi,
                                              n_threads = 1L) {
  .Call("make_C_B_tapered_from_pairs",
        as.integer(C_B_pairs$pair_l),
        as.integer(C_B_pairs$pair_k),
        as.numeric(C_B_pairs$pair_d),
        as.numeric(C_B_pairs$pair_wtap),
        as.integer(C_B_pairs$n_b),
        as.numeric(phi),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_pred_B_tapered_Rcall <- function(pred_coords, H_BA_comp, phi, gamma,
                                        taper_code = 1L, n_threads = 1L) {
  pred_coords <- as.matrix(pred_coords)

  .Call("make_C_AB_tapered",
        as.numeric(pred_coords[, 1]),
        as.numeric(pred_coords[, 2]),
        as.integer(H_BA_comp$plot_start),
        as.numeric(H_BA_comp$h),
        as.numeric(H_BA_comp$x),
        as.numeric(H_BA_comp$y),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_pred_tapered_dense_Rcall <- function(pred_coords, phi, gamma,
                                            taper_code = 1L, n_threads = 1L) {
  pred_coords <- as.matrix(pred_coords)

  .Call("make_C_A_tapered_dense",
        as.numeric(pred_coords[, 1]),
        as.numeric(pred_coords[, 2]),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

gaussian_metrop_sampler_Rcall <- function(y_B, mu_y, X_B_V_beta_X_B, D_h,
                                          C_B_pairs = NULL, spatial,
                                          starting, tuning,
                                          n_batch, batch_length,
                                          accept_rate, c0 = 1, c1 = 0.8,
                                          report = 100, verbose = TRUE,
                                          d_h_bar,
                                          tau_B_shape, tau_B_scale,
                                          sigma_shape = NA_real_,
                                          sigma_scale = NA_real_,
                                          phi_lower = NA_real_,
                                          phi_upper = NA_real_,
                                          n_threads = 1L) {
  if (isTRUE(spatial)) {
    pair_l <- C_B_pairs$pair_l
    pair_k <- C_B_pairs$pair_k
    pair_d <- C_B_pairs$pair_d
    pair_wtap <- C_B_pairs$pair_wtap
    n_b_pairs <- C_B_pairs$n_b
  } else {
    pair_l <- integer()
    pair_k <- integer()
    pair_d <- numeric()
    pair_wtap <- numeric()
    n_b_pairs <- length(y_B)
  }

  out <- .Call(
    "gaussian_metrop_sampler",
    as.numeric(y_B),
    as.numeric(mu_y),
    as.matrix(X_B_V_beta_X_B),
    as.matrix(D_h),
    as.integer(pair_l),
    as.integer(pair_k),
    as.numeric(pair_d),
    as.numeric(pair_wtap),
    as.integer(n_b_pairs),
    as.logical(spatial),
    as.numeric(starting),
    as.numeric(tuning),
    as.integer(n_batch),
    as.integer(batch_length),
    as.numeric(accept_rate),
    as.numeric(c0),
    as.numeric(c1),
    as.integer(report),
    as.logical(verbose),
    as.numeric(d_h_bar),
    as.numeric(tau_B_shape),
    as.numeric(tau_B_scale),
    as.numeric(sigma_shape),
    as.numeric(sigma_scale),
    as.numeric(phi_lower),
    as.numeric(phi_upper),
    as.integer(n_threads),
    PACKAGE = "cosTapered"
  )

  colnames(out$p.theta.samples) <- names(starting)
  colnames(out$Sigma0) <- names(starting)
  rownames(out$Sigma0) <- names(starting)
  colnames(out$proposal.cov) <- names(starting)
  rownames(out$proposal.cov) <- names(starting)
  out$call <- match.call()
  out
}

binomial_pg_sampler_Rcall <- function(y_B, trials, X_B, C_B_pairs = NULL,
                                      offset = NULL,
                                      family_code = 1L,
                                      nb_size = NA_real_,
                                      spatial, beta_mu, V_beta,
                                      starting_beta, starting_sigma_sq = NA_real_,
                                      starting_phi = NA_real_,
                                      tuning = c(log_sigma_sq = 1, z_phi = 1),
                                      n_batch, batch_length,
                                      accept_rate, c0 = 1, c1 = 0.8,
                                      report = 100, verbose = TRUE,
                                      sigma_shape = NA_real_,
                                      sigma_scale = NA_real_,
                                      phi_lower = NA_real_,
                                      phi_upper = NA_real_,
                                      n_threads = 1L) {
  if (isTRUE(spatial)) {
    pair_l <- C_B_pairs$pair_l
    pair_k <- C_B_pairs$pair_k
    pair_d <- C_B_pairs$pair_d
    pair_wtap <- C_B_pairs$pair_wtap
    n_b_pairs <- C_B_pairs$n_b
  } else {
    pair_l <- integer()
    pair_k <- integer()
    pair_d <- numeric()
    pair_wtap <- numeric()
    n_b_pairs <- length(y_B)
    tuning <- numeric()
  }
  if (is.null(offset)) {
    offset <- rep(0, length(y_B))
  }

  out <- .Call(
    "binomial_pg_sampler",
    as.numeric(y_B),
    as.integer(trials),
    as.numeric(offset),
    as.integer(family_code),
    as.numeric(nb_size),
    as.matrix(X_B),
    as.integer(pair_l),
    as.integer(pair_k),
    as.numeric(pair_d),
    as.numeric(pair_wtap),
    as.integer(n_b_pairs),
    as.logical(spatial),
    as.numeric(beta_mu),
    as.matrix(V_beta),
    as.numeric(starting_beta),
    as.numeric(starting_sigma_sq),
    as.numeric(starting_phi),
    as.numeric(tuning),
    as.integer(n_batch),
    as.integer(batch_length),
    as.numeric(accept_rate),
    as.numeric(c0),
    as.numeric(c1),
    as.integer(report),
    as.logical(verbose),
    as.numeric(sigma_shape),
    as.numeric(sigma_scale),
    as.numeric(phi_lower),
    as.numeric(phi_upper),
    as.integer(n_threads),
    PACKAGE = "cosTapered"
  )

  colnames(out$beta_samples) <- colnames(X_B)
  colnames(out$omega_B_samples) <- paste0("B_", seq_along(y_B))
  colnames(out$eta_B_samples) <- paste0("B_", seq_along(y_B))

  if (isTRUE(spatial)) {
    colnames(out$theta_z_samples) <- c("log_sigma_sq", "z_phi")
    colnames(out$Sigma0) <- c("log_sigma_sq", "z_phi")
    rownames(out$Sigma0) <- c("log_sigma_sq", "z_phi")
    colnames(out$proposal.cov) <- c("log_sigma_sq", "z_phi")
    rownames(out$proposal.cov) <- c("log_sigma_sq", "z_phi")
  }

  out$call <- match.call()
  out
}
