# Exact (up to floating-point roundoff) covariance sums on a square raster.
# Group pair weights by squared lattice distance once, then reuse them for
# every phi. No distance binning or change of prediction support is involved.
cos_areal_grid_indices <- function(coords, grid) {
  if (is.null(grid) || length(grid$res) != 2L ||
      length(grid$origin) != 2L || any(!is.finite(c(grid$res, grid$origin))) ||
      any(grid$res <= 0) ||
      abs(grid$res[1] - grid$res[2]) > 1e-12 * max(grid$res)) return(NULL)
  coords <- as.matrix(coords)
  if (ncol(coords) != 2L || any(!is.finite(coords))) return(NULL)
  z <- sweep(sweep(coords, 2L, grid$origin, "-"), 2L, grid$res, "/")
  # Accept only machine-level coordinate noise, not off-grid prediction sites.
  if (any(abs(z - round(z)) > 1e-8) ||
      any(abs(z) > .Machine$integer.max / 4)) return(NULL)
  z <- round(z)
  storage.mode(z) <- "integer"
  z
}

cos_areal_grid_self_weights <- function(index, h, n_bin, n_threads = 1L) {
  .Call("cos_areal_grid_self_weights_cpp",
        as.integer(index[, 1]), as.integer(index[, 2]), as.numeric(h),
        as.integer(n_bin), as.integer(n_threads), PACKAGE = "cosTapered")
}

cos_prepare_areal_grid_covariance <- function(U_blocks, H_BA_comp, gamma,
                                               taper_code, need_self,
                                               n_threads = 1L,
                                               verbose = FALSE,
                                               min_cells = 2048L,
                                               min_cross_cells = 10000L,
                                               cache_bytes = 256 * 1024^2) {
  blocks <- U_blocks$blocks
  cache <- vector("list", length(blocks))
  grid <- U_blocks$grid
  if (is.null(grid)) return(cache) # Previously saved block objects still work.
  b_index <- cos_areal_grid_indices(cbind(H_BA_comp$x, H_BA_comp$y), grid)
  n_b <- length(H_BA_comp$plot_start) - 1L
  for (k in order(vapply(blocks, function(b) length(b$h), integer(1)),
                  decreasing = TRUE)) {
    b <- blocks[[k]]
    if (length(b$h) < min_cells || any(!is.finite(b$h)) || any(b$h < 0)) next
    index <- cos_areal_grid_indices(b$coords, grid)
    if (is.null(index)) next
    # Integer squared distances identify exact lags on a square lattice.
    n_bin <- ceiling((gamma / grid$res[1])^2)
    if (!is.finite(n_bin) || n_bin < 1 || n_bin > 1e6) next
    do_cross <- !is.null(b_index) && length(b$h) >= min_cross_cells &&
      all(is.finite(H_BA_comp$h)) && all(H_BA_comp$h >= 0) &&
      8 * n_bin * (n_b + 1) <= cache_bytes
    if (!need_self && !do_cross) next
    if (8 * n_bin > cache_bytes) next
    if (verbose) message("Preparing grid covariance weights: areal unit ", k,
                         " (", length(b$h), " cells)")
    self <- if (need_self) cos_areal_grid_self_weights(index, b$h, n_bin, n_threads) else NULL
    cross <- if (do_cross) {
      .Call("cos_areal_grid_cross_weights",
            as.integer(index[, 1]), as.integer(index[, 2]), as.numeric(b$h),
            as.integer(b_index[, 1]), as.integer(b_index[, 2]),
            as.numeric(H_BA_comp$h), as.integer(H_BA_comp$plot_start),
            as.integer(n_bin), as.integer(n_threads), PACKAGE = "cosTapered")
    } else NULL
    if (is.null(self) && is.null(cross)) next
    keep <- rep(FALSE, n_bin)
    if (!is.null(self)) keep <- keep | self != 0
    if (!is.null(cross)) keep <- keep | rowSums(cross) != 0
    d <- sqrt(which(keep) - 1) * grid$res[1]
    r <- d / gamma
    taper <- if (taper_code == 1L) (1 - r)^4 * (1 + 4 * r) else
      1 - 1.5 * r + 0.5 * r^3
    taper[d >= gamma] <- 0
    value <- list(distance = d, taper = taper,
                  self = if (!is.null(self)) self[keep] else NULL,
                  cross = if (!is.null(cross)) cross[keep, , drop = FALSE] else NULL)
    cache[k] <- list(value)
    cache_bytes <- cache_bytes - as.numeric(object.size(value))
  }
  cache
}

cos_eval_areal_grid_covariance <- function(cache, phi) {
  if (is.null(cache)) return(list(C_UU = NULL, C_UB = NULL))
  kernel <- exp(-phi * cache$distance) * cache$taper
  list(C_UU = if (!is.null(cache$self)) sum(kernel * cache$self) else NULL,
       C_UB = if (!is.null(cache$cross)) as.numeric(crossprod(kernel, cache$cross)) else NULL)
}
