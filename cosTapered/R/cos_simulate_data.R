cos_simulate_data <- function(X_rast,
                              forest_sf,
                              n_B = 36L,
                              plot_radius,
                              gamma,
                              beta,
                              tau_B_sq,
                              sigma_sq,
                              phi,
                              taper_code = 1L,
                              design_seed = NULL,
                              sim_seed = NULL,
                              n_threads = 1L,
                              verbose = FALSE) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(X_rast, "Raster")) {
    stop("X_rast must be a raster RasterLayer, RasterStack, or RasterBrick.")
  }
  if (!inherits(forest_sf, "sf")) {
    stop("forest_sf must be an sf object containing the simulation boundary.")
  }
  if (sf::st_is_longlat(forest_sf)) {
    stop("forest_sf must use a projected CRS with linear units.")
  }
  if (!cos_same_crs(forest_sf, raster::crs(X_rast))) {
    stop("forest_sf and X_rast must have the same CRS.")
  }
  if (!isTRUE(sf::st_crs(forest_sf) == sf::st_crs(raster::crs(X_rast)))) {
    suppressWarnings(sf::st_crs(forest_sf) <- sf::st_crs(raster::crs(X_rast)))
  }

  n_B <- as.integer(n_B)
  n_threads <- as.integer(n_threads)

  if (!is.finite(n_B) || n_B < 3L) stop("n_B must be an integer greater than or equal to 3.")
  if (missing(plot_radius) || length(plot_radius) != 1L ||
      !is.finite(plot_radius) || plot_radius <= 0) {
    stop("plot_radius must be a positive distance in the CRS units.")
  }
  if (missing(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a positive taper distance in the CRS units.")
  }
  if (!taper_code %in% c(1L, 2L)) {
    stop("taper_code must be 1L (Wendland) or 2L (spherical).")
  }
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L

  beta_names <- names(beta)
  beta <- as.numeric(beta)
  names(beta) <- beta_names
  if (is.null(names(beta)) || any(names(beta) == "")) {
    stop("beta must be a named numeric vector.")
  }
  if (any(!is.finite(beta))) stop("beta contains non-finite values.")

  if (!is.finite(tau_B_sq) || tau_B_sq <= 0) stop("tau_B_sq must be positive and finite.")
  if (!is.finite(sigma_sq) || sigma_sq <= 0) stop("sigma_sq must be positive and finite.")
  if (!is.finite(phi) || phi <= 0) stop("phi must be positive and finite.")

  x_names <- names(X_rast)
  if (is.null(x_names) || anyNA(x_names) || any(x_names == "")) {
    stop("X_rast must have valid layer names.")
  }
  if (anyDuplicated(x_names)) {
    stop("X_rast layer names must be unique.")
  }

  # ------------------------------------------------
  # Add intercept layer if requested by beta
  # ------------------------------------------------

  if ("intercept" %in% names(beta) && !"intercept" %in% x_names) {
    intercept <- X_rast[[1]]
    intercept[] <- ifelse(is.na(raster::getValues(intercept)), NA_real_, 1)
    names(intercept) <- "intercept"
    X_rast <- raster::stack(intercept, X_rast)
    x_names <- names(X_rast)
  }

  if (!all(names(beta) %in% x_names)) {
    stop("All beta names must match X_rast layer names. Missing: ",
         paste(setdiff(names(beta), x_names), collapse = ", "))
  }
  if (!all(x_names %in% names(beta))) {
    stop("beta must include coefficients for every X_rast layer. Missing: ",
         paste(setdiff(x_names, names(beta)), collapse = ", "))
  }

  beta <- beta[x_names]

  # ------------------------------------------------
  # Build systematic observed-support polygons
  # ------------------------------------------------

  if (!is.null(design_seed)) set.seed(design_seed)

  forest_sf <- sf::st_make_valid(forest_sf)
  forest_sf <- suppressWarnings(sf::st_collection_extract(forest_sf, "POLYGON"))
  forest_sf <- suppressWarnings(sf::st_cast(forest_sf, "MULTIPOLYGON"))
  forest_sf <- forest_sf[!sf::st_is_empty(forest_sf), , drop = FALSE]
  if (nrow(forest_sf) < 1L) stop("No usable forest polygon found.")

  forest_union <- sf::st_union(forest_sf)

  B_centers_geom <- sf::st_sample(
    forest_union,
    size = n_B,
    type = "regular"
  )

  B_centers <- sf::st_sf(
    center_id = seq_along(B_centers_geom),
    geometry = B_centers_geom,
    crs = sf::st_crs(raster::crs(X_rast))
  )

  if (nrow(B_centers) < 3L) {
    stop("Fewer than 3 observed-support centers were generated inside forest_sf.")
  }

  B_sf <- sf::st_buffer(B_centers, dist = plot_radius)
  B_sf <- suppressWarnings(sf::st_intersection(B_sf, forest_union))
  B_sf <- sf::st_make_valid(B_sf)
  B_sf <- suppressWarnings(sf::st_collection_extract(B_sf, "POLYGON"))
  B_sf <- suppressWarnings(sf::st_cast(B_sf, "MULTIPOLYGON"))
  B_sf <- B_sf[!sf::st_is_empty(B_sf), , drop = FALSE]

  if (nrow(B_sf) < 3L) {
    stop("Fewer than 3 observed-support polygons remain after clipping to forest_sf.")
  }

  B_sf$B_id <- seq_len(nrow(B_sf))
  B_sf$area_m2 <- as.numeric(sf::st_area(B_sf))
  B_sf$y_B <- 0

  full_plot_area <- pi * plot_radius^2
  n_boundary_clipped <- sum(B_sf$area_m2 < 0.99 * full_plot_area)

  # ------------------------------------------------
  # Build COS objects and simulate observed response
  # ------------------------------------------------

  prep <- cos_prepare(
    X_rast = X_rast,
    B_sf = B_sf,
    response_col = "y_B",
    gamma = gamma,
    taper_code = taper_code,
    n_threads = n_threads,
    verbose = verbose
  )

  d_h_bar <- mean(diag(prep$D_h))
  tau_sq <- tau_B_sq / d_h_bar

  C_B <- make_C_B_tapered_from_pairs_Rcall(
    C_B_pairs = prep$C_B_pairs,
    phi = phi,
    n_threads = n_threads
  )

  if (!is.null(sim_seed)) set.seed(sim_seed)

  n_b <- nrow(prep$X_B)
  R_omega <- chol(sigma_sq * C_B)
  R_eps <- chol(tau_sq * as.matrix(prep$D_h))

  omega_B <- as.numeric(t(R_omega) %*% stats::rnorm(n_b))
  eps_B <- as.numeric(t(R_eps) %*% stats::rnorm(n_b))
  eta_B <- as.numeric(prep$X_B %*% beta + omega_B)
  y_B <- eta_B + eps_B

  B_sf$y_B <- y_B
  B_sf$eta_B_true <- eta_B
  B_sf$omega_B_true <- omega_B
  B_sf$eps_B_true <- eps_B

  # ------------------------------------------------
  # Return simulated data object
  # ------------------------------------------------

  truth <- list(
    beta_true = beta,
    tau_B_sq_true = tau_B_sq,
    tau_sq_true = tau_sq,
    d_h_bar = d_h_bar,
    sigma_sq_true = sigma_sq,
    eff_range_true = 3 / phi,
    phi_true = phi,
    gamma = gamma,
    taper_code = taper_code,
    plot_radius_m = plot_radius,
    n_plot_target = n_B,
    n_boundary_clipped = n_boundary_clipped,
    design_seed = design_seed,
    sim_seed = sim_seed,
    seed = sim_seed,
    n_B = nrow(B_sf),
    n_A = nrow(prep$A_df),
    row_sum_summary = summary(prep$row_sums)
  )

  out <- list(
    X_rast = X_rast,
    forest_sf = forest_sf,
    B_sf = B_sf,
    truth = truth,
    prep = prep
  )

  class(out) <- "cos_simulation"
  out
}
