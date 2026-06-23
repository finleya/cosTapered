cos_make_areal_blocks <- function(U_sf,
                                  X_rast,
                                  id_col = NULL,
                                  missing = c("error", "drop", "renormalize"),
                                  row_sum_tol = 1e-5,
                                  verbose = TRUE) {
  t_start <- proc.time()

  # ------------------------------------------------
  # Check inputs
  # ------------------------------------------------

  if (!inherits(U_sf, "sf")) {
    stop("U_sf must be an sf object containing areal prediction polygons.")
  }
  if (!inherits(X_rast, "Raster")) {
    stop("X_rast must be a raster RasterLayer, RasterStack, or RasterBrick.")
  }
  if (sf::st_is_longlat(U_sf)) {
    stop("U_sf must use a projected CRS with linear units.")
  }
  if (!cos_same_crs(U_sf, raster::crs(X_rast))) {
    stop("U_sf and X_rast must have the same CRS.")
  }
  if (!isTRUE(sf::st_crs(U_sf) == sf::st_crs(raster::crs(X_rast)))) {
    suppressWarnings(sf::st_crs(U_sf) <- sf::st_crs(raster::crs(X_rast)))
  }
  if (!is.null(id_col) && !id_col %in% names(U_sf)) {
    stop("id_col not found in U_sf.")
  }
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop("Package exactextractr is required.")
  }
  missing <- match.arg(missing)

  x_names <- names(X_rast)
  if (is.null(x_names) || anyNA(x_names) || any(x_names == "")) {
    stop("X_rast must have valid layer names.")
  }
  if (anyDuplicated(x_names)) {
    stop("X_rast layer names must be unique.")
  }

  # ------------------------------------------------
  # Prepare areal prediction supports U
  # ------------------------------------------------

  U_sf <- sf::st_make_valid(U_sf)
  U_sf$U_id <- seq_len(nrow(U_sf))

  if (is.null(id_col)) {
    U_sf$U_label <- as.character(U_sf$U_id)
  } else {
    U_sf$U_label <- as.character(U_sf[[id_col]])
  }

  U_area <- as.numeric(sf::st_area(U_sf))
  if (any(!is.finite(U_area) | U_area <= 0)) {
    stop("Some U_sf polygons have non-positive or non-finite area.")
  }

  # ------------------------------------------------
  # Extract raster cells and area weights for each U
  # ------------------------------------------------

  cell_area <- prod(raster::res(X_rast))

  ex <- exactextractr::exact_extract(
    X_rast,
    U_sf,
    include_cell = TRUE,
    include_xy = TRUE,
    progress = verbose
  )

  blocks <- vector("list", length(ex))

  for (k in seq_along(ex)) {
    z <- ex[[k]]

    if (is.null(z) || nrow(z) == 0L) {
      blocks[[k]] <- NULL
      next
    }

    z <- as.data.frame(z)

    # exactextractr may call a single unnamed raster layer "value".
    if (!all(x_names %in% names(z)) && length(x_names) == 1L && "value" %in% names(z)) {
      names(z)[names(z) == "value"] <- x_names
    }

    if (!all(c("coverage_fraction", "cell", "x", "y") %in% names(z))) {
      stop("Unexpected exactextractr output when building areal blocks.")
    }
    if (!all(x_names %in% names(z))) {
      stop("Could not find all X_rast layer names in exactextractr output.")
    }

    df <- data.frame(
      cell_id = z$cell,
      h = z$coverage_fraction * cell_area / U_area[k],
      x = z$x,
      y = z$y,
      z[, x_names, drop = FALSE],
      check.names = FALSE
    )

    df <- df[is.finite(df$h) & df$h > 0, , drop = FALSE]
    if (nrow(df) == 0L) {
      blocks[[k]] <- NULL
      next
    }

    keep <- stats::complete.cases(df[, x_names, drop = FALSE])
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0L) {
      if (missing == "error") {
        stop("Areal prediction unit ", k, " has no complete raster covariate cells.")
      } else {
        blocks[[k]] <- NULL
        next
      }
    }

    # Aggregate in case exact extraction returns duplicate cell records.
    h_sum <- stats::aggregate(h ~ cell_id, data = df, FUN = sum)
    first_ind <- match(h_sum$cell_id, df$cell_id)

    df2 <- data.frame(
      cell_id = h_sum$cell_id,
      h = h_sum$h,
      x = df$x[first_ind],
      y = df$y[first_ind],
      df[first_ind, x_names, drop = FALSE],
      check.names = FALSE
    )
    df2 <- df2[order(df2$cell_id), , drop = FALSE]

    X_UA <- as.matrix(df2[, x_names, drop = FALSE])
    storage.mode(X_UA) <- "double"
    colnames(X_UA) <- x_names

    coords <- as.matrix(df2[, c("x", "y")])
    storage.mode(coords) <- "double"

    h <- as.numeric(df2$h)
    row_sum <- sum(h)

    if (row_sum < 1 - row_sum_tol && missing == "error") {
      stop(
        "Incomplete X_rast values leave areal prediction support weights below 1 for unit ",
        k, ". Use missing = 'drop' to keep the covered-cell weights or ",
        "missing = 'renormalize' to average over covered cells only."
      )
    }
    if (missing == "renormalize") {
      if (!is.finite(row_sum) || row_sum <= 0) {
        stop("Cannot renormalize support weights because areal prediction unit ", k, " has no complete raster cells.")
      }
      h <- h / row_sum
      row_sum <- sum(h)
    }

    X_U <- as.numeric(crossprod(h, X_UA))
    names(X_U) <- x_names

    if (row_sum < 1 - row_sum_tol) {
      warning("Some H_UA rows sum to less than 1. Weights are left as constructed; inspect summary(U_blocks)$row_sum_summary.")
    }
    if (row_sum > 1 + row_sum_tol) {
      warning("Some H_UA rows sum to greater than 1. Weights are left as constructed; inspect summary(U_blocks)$row_sum_summary.")
    }

    blocks[[k]] <- list(
      U_id = U_sf$U_id[k],
      U_label = U_sf$U_label[k],
      n_cell = nrow(df2),
      row_sum = row_sum,
      h = h,
      coords = coords,
      X_UA = X_UA,
      X_U = X_U,
      cell_id = df2$cell_id
    )
  }

  keep <- !vapply(blocks, is.null, logical(1))
  if (!all(keep)) {
    warning(sum(!keep), " areal prediction unit(s) had no usable raster cells and were omitted.")
  }

  blocks <- blocks[keep]
  U_sf_kept <- U_sf[keep, , drop = FALSE]
  row_sums <- vapply(blocks, `[[`, numeric(1), "row_sum")

  # ------------------------------------------------
  # Return areal-block object
  # ------------------------------------------------

  dt <- proc.time() - t_start
  timing <- c(
    user = unname(dt[["user.self"]]),
    system = unname(dt[["sys.self"]]),
    elapsed = unname(dt[["elapsed"]])
  )

  out <- list(
    U_sf = U_sf_kept,
    blocks = blocks,
    row_sums = row_sums,
    x_names = x_names,
    missing = missing,
    timing = timing,
    call = match.call()
  )

  class(out) <- "cos_areal_blocks"
  out
}
