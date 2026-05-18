cos_prepare <- function(X_rast,
                        B_sf,
                        response_col = NULL,
                        y_B = NULL,
                        gamma,
                        spatial = TRUE,
                        taper_code = 1L,
                        n_threads = 1L,
                        row_sum_tol = 1e-5,
                        verbose = TRUE) {
  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(X_rast, "Raster")) {
    stop("X_rast must be a raster RasterLayer, RasterStack, or RasterBrick.")
  }
  if (!inherits(B_sf, "sf")) {
    stop("B_sf must be an sf object containing observed-support polygons.")
  }
  if (sf::st_is_longlat(B_sf)) {
    stop("B_sf must use a projected CRS with linear units.")
  }
  if (!cos_same_crs(B_sf, raster::crs(X_rast))) {
    stop("B_sf and X_rast must have the same CRS.")
  }
  if (!isTRUE(sf::st_crs(B_sf) == sf::st_crs(raster::crs(X_rast)))) {
    suppressWarnings(sf::st_crs(B_sf) <- sf::st_crs(raster::crs(X_rast)))
  }
  spatial <- isTRUE(spatial)
  if (spatial && (missing(gamma) || length(gamma) != 1 || !is.finite(gamma) || gamma <= 0)) {
    stop("gamma must be supplied as a positive taper distance in the CRS units.")
  }
  if (!spatial) {
    gamma <- NA_real_
  }
  if (!taper_code %in% c(1L, 2L)) {
    stop("taper_code must be 1L (Wendland) or 2L (spherical).")
  }
  n_threads <- as.integer(n_threads)
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L

  x_names <- names(X_rast)
  if (is.null(x_names) || anyNA(x_names) || any(x_names == "")) {
    stop("X_rast must have valid layer names.")
  }
  if (anyDuplicated(x_names)) {
    stop("X_rast layer names must be unique.")
  }

  if (is.null(y_B)) {
    if (is.null(response_col)) {
      stop("Provide y_B or response_col.")
    }
    if (!response_col %in% names(B_sf)) {
      stop("response_col not found in B_sf.")
    }
    y_B <- B_sf[[response_col]]
  }
  y_B <- as.numeric(y_B)
  if (length(y_B) != nrow(B_sf)) {
    stop("length(y_B) must equal nrow(B_sf).")
  }
  if (any(!is.finite(y_B))) {
    stop("y_B contains non-finite values.")
  }

  # ------------------------------------------------
  # Prepare observed supports B
  # ------------------------------------------------

  B_sf <- sf::st_make_valid(B_sf)
  B_sf$B_id <- seq_len(nrow(B_sf))
  n_b <- nrow(B_sf)

  B_area <- as.numeric(sf::st_area(B_sf))
  if (any(!is.finite(B_area) | B_area <= 0)) {
    stop("Some B_sf polygons have non-positive or non-finite area.")
  }

  # ------------------------------------------------
  # Exact raster-cell support weights H_BA
  # ------------------------------------------------

  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop("Package exactextractr is required.")
  }

  r_template <- X_rast[[1]]
  cell_area <- prod(raster::res(r_template))

  ex <- exactextractr::exact_extract(
    r_template,
    B_sf,
    include_cell = TRUE,
    include_xy = TRUE,
    progress = verbose
  )

  trip_list <- vector("list", length(ex))

  for (l in seq_along(ex)) {
    z <- ex[[l]]
    if (is.null(z) || nrow(z) == 0) next

    trip_list[[l]] <- data.frame(
      B_id = l,
      cell_id = z$cell,
      h = z$coverage_fraction * cell_area / B_area[l],
      x = z$x,
      y = z$y
    )
  }

  trip <- do.call(rbind, trip_list)
  if (is.null(trip) || nrow(trip) == 0) {
    stop("No raster cells intersected B_sf.")
  }

  trip <- trip[is.finite(trip$h) & trip$h > 0, , drop = FALSE]
  if (nrow(trip) == 0) {
    stop("No positive support weights were found.")
  }

  trip <- stats::aggregate(
    h ~ B_id + cell_id,
    data = trip,
    FUN = sum
  )
  trip <- trip[order(trip$B_id, trip$cell_id), , drop = FALSE]

  unique_cells <- sort(unique(trip$cell_id))
  A_df <- data.frame(
    cell_id = unique_cells,
    stringsAsFactors = FALSE
  )

  xy <- raster::xyFromCell(r_template, unique_cells)
  A_df$x <- xy[, 1]
  A_df$y <- xy[, 2]

  X <- raster::extract(X_rast, unique_cells)

  # raster::extract() returns a vector for a single-layer raster and a
  # matrix/data frame for multi-layer rasters. Force a standard numeric
  # matrix so downstream code is unchanged.
  if (is.null(dim(X))) {
    X <- matrix(X, ncol = 1)
  } else {
    X <- as.matrix(X)
  }

  storage.mode(X) <- "double"

  if (ncol(X) != length(x_names)) {
    stop("Number of extracted X columns does not match number of X_rast layers.")
  }

  colnames(X) <- x_names

  keep_A <- stats::complete.cases(X)
  if (!all(keep_A)) {
    if (verbose) {
      message("Dropping ", sum(!keep_A), " fine cell(s) with incomplete X_rast values.")
    }
    A_df <- A_df[keep_A, , drop = FALSE]
    X <- X[keep_A, , drop = FALSE]
    unique_cells <- A_df$cell_id
    trip <- trip[trip$cell_id %in% unique_cells, , drop = FALSE]
  }

  if (nrow(A_df) == 0) {
    stop("No complete fine-support covariate cells remain after filtering.")
  }

  A_df$A_id <- seq_len(nrow(A_df))
  trip <- merge(trip, A_df[, c("cell_id", "A_id")], by = "cell_id", all.x = TRUE)
  trip <- trip[order(trip$B_id, trip$A_id), , drop = FALSE]

  H_BA <- Matrix::sparseMatrix(
    i = trip$B_id,
    j = trip$A_id,
    x = trip$h,
    dims = c(n_b, nrow(A_df))
  )
  H_BA <- methods::as(H_BA, "dgCMatrix")

  row_sums <- as.numeric(Matrix::rowSums(H_BA))

  if (any(row_sums < 1 - row_sum_tol)) {
    warning("Some H_BA rows sum to less than 1. Weights are left as constructed; inspect summary(prep)$row_sum_summary.")
  }
  if (any(row_sums > 1 + row_sum_tol)) {
    warning("Some H_BA rows sum to greater than 1. Weights are left as constructed; inspect summary(prep)$row_sum_summary.")
  }

  # ------------------------------------------------
  # Observed-support design and nugget matrices
  # ------------------------------------------------

  A_coords <- as.matrix(A_df[, c("x", "y")])
  X_B <- as.matrix(H_BA %*% X)
  colnames(X_B) <- colnames(X)
  D_h <- as.matrix(H_BA %*% Matrix::t(H_BA))

  # ------------------------------------------------
  # Compressed H_BA representation for C++ covariance code
  # ------------------------------------------------

  Hs <- Matrix::summary(H_BA)
  Hs <- as.data.frame(Hs)

  if (ncol(Hs) < 3) {
    stop("Unexpected summary(H_BA) structure.")
  }

  names(Hs)[1:3] <- c("i", "j", "x")
  Hs <- Hs[order(Hs$i, Hs$j), , drop = FALSE]
  counts <- tabulate(Hs$i, nbins = n_b)

  H_BA_comp <- list(
    plot_start = c(1L, as.integer(cumsum(counts) + 1L)),
    h = as.numeric(Hs$x),
    x = as.numeric(A_df$x[Hs$j]),
    y = as.numeric(A_df$y[Hs$j]),
    B_id = as.integer(Hs$i),
    A_id = as.integer(Hs$j)
  )

  # ------------------------------------------------
  # Tapered covariance pair precomputation
  # ------------------------------------------------

  C_B_pairs <- NULL
  if (spatial) {
    C_B_pairs <- prep_C_B_tapered_pairs_Rcall(
      plot_start = H_BA_comp$plot_start,
      h = H_BA_comp$h,
      x = H_BA_comp$x,
      y = H_BA_comp$y,
      gamma = gamma,
      taper_code = taper_code,
      n_threads = n_threads
    )
  }

  # ------------------------------------------------
  # Return prep object
  # ------------------------------------------------

  out <- list(
    B_sf = B_sf,
    y_B = y_B,
    A_df = A_df,
    A_coords = A_coords,
    X = X,
    X_B = X_B,
    H_BA = H_BA,
    H_BA_comp = H_BA_comp,
    D_h = D_h,
    row_sums = row_sums,
    B_area = B_area,
    C_B_pairs = C_B_pairs,
    gamma = gamma,
    spatial = spatial,
    taper_code = taper_code,
    n_threads = n_threads,
    call = match.call()
  )

  class(out) <- "cos_prep"
  out
}
