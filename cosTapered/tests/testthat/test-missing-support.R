test_that("cos_prepare handles incomplete raster support explicitly", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("exactextractr")

  chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
  B_file <- system.file("extdata", "example_B.gpkg", package = "cosTapered")
  if (chm_file == "") {
    chm_file <- testthat::test_path("..", "..", "inst", "extdata", "example_chm.tif")
  }
  if (B_file == "") {
    B_file <- testthat::test_path("..", "..", "inst", "extdata", "example_B.gpkg")
  }
  skip_if(!file.exists(chm_file) || !file.exists(B_file), "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"
  B_sf <- sf::st_read(B_file, quiet = TRUE)[seq_len(2), ]

  ex <- suppressWarnings(
    exactextractr::exact_extract(chm, B_sf[1, ], include_cell = TRUE)
  )[[1]]
  ex <- ex[is.finite(ex$coverage_fraction) & ex$coverage_fraction > 0, , drop = FALSE]
  skip_if(nrow(ex) == 0L, "test polygon has no raster overlap")

  chm_bad <- chm
  chm_bad[ex$cell[1]] <- NA_real_
  names(chm_bad) <- "chm"

  intercept <- chm_bad
  intercept[] <- ifelse(is.na(raster::getValues(chm_bad)), NA_real_, 1)
  names(intercept) <- "intercept"
  X_bad <- raster::stack(intercept, chm_bad)
  names(X_bad) <- c("intercept", "chm")

  expect_error(
    cos_prepare(
      X_rast = X_bad,
      B_sf = B_sf,
      response_col = "y_B",
      spatial = FALSE,
      verbose = FALSE
    ),
    "Incomplete X_rast values|No raster cells intersected|No positive support weights|No complete fine-support"
  )

  expect_warning(
    prep_drop <- cos_prepare(
      X_rast = X_bad,
      B_sf = B_sf,
      response_col = "y_B",
      spatial = FALSE,
      missing = "drop",
      verbose = FALSE
    ),
    "rows sum to less than 1"
  )
  expect_s3_class(prep_drop, "cos_prep")
  expect_true(any(prep_drop$row_sums < 1))
  expect_equal(prep_drop$missing, "drop")

  prep_renorm <- cos_prepare(
    X_rast = X_bad,
    B_sf = B_sf,
    response_col = "y_B",
    spatial = FALSE,
    missing = "renormalize",
    verbose = FALSE
  )
  expect_s3_class(prep_renorm, "cos_prep")
  expect_true(all(abs(prep_renorm$row_sums - 1) < 1e-8))
  expect_equal(prep_renorm$missing, "renormalize")

  chm_all_bad <- chm
  chm_all_bad[ex$cell] <- NA_real_
  names(chm_all_bad) <- "chm"
  intercept_all_bad <- chm_all_bad
  intercept_all_bad[] <- ifelse(is.na(raster::getValues(chm_all_bad)), NA_real_, 1)
  names(intercept_all_bad) <- "intercept"
  X_all_bad <- raster::stack(intercept_all_bad, chm_all_bad)
  names(X_all_bad) <- c("intercept", "chm")

  expect_error(
    cos_prepare(
      X_rast = X_all_bad,
      B_sf = B_sf[1, ],
      response_col = "y_B",
      spatial = FALSE,
      verbose = FALSE
    ),
    "Incomplete X_rast values|No raster cells intersected|No positive support weights|No complete fine-support"
  )
  expect_error(
    cos_prepare(
      X_rast = X_all_bad,
      B_sf = B_sf[1, ],
      response_col = "y_B",
      spatial = FALSE,
      missing = "renormalize",
      verbose = FALSE
    ),
    "No raster cells intersected|Cannot renormalize|No positive support weights|No complete fine-support"
  )
})

test_that("cos_make_areal_blocks handles incomplete raster support explicitly", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("exactextractr")

  chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
  U_file <- system.file("extdata", "example_U.gpkg", package = "cosTapered")
  if (chm_file == "") {
    chm_file <- testthat::test_path("..", "..", "inst", "extdata", "example_chm.tif")
  }
  if (U_file == "") {
    U_file <- testthat::test_path("..", "..", "inst", "extdata", "example_U.gpkg")
  }
  skip_if(!file.exists(chm_file) || !file.exists(U_file), "example data not installed")

  chm <- raster::raster(chm_file)
  names(chm) <- "chm"
  U_sf <- sf::st_read(U_file, quiet = TRUE)[1, ]

  ex <- suppressWarnings(
    exactextractr::exact_extract(chm, U_sf, include_cell = TRUE)
  )[[1]]
  ex <- ex[is.finite(ex$coverage_fraction) & ex$coverage_fraction > 0, , drop = FALSE]
  skip_if(nrow(ex) == 0L, "test polygon has no raster overlap")

  chm_bad <- chm
  chm_bad[ex$cell[1]] <- NA_real_
  names(chm_bad) <- "chm"

  intercept <- chm_bad
  intercept[] <- ifelse(is.na(raster::getValues(chm_bad)), NA_real_, 1)
  names(intercept) <- "intercept"
  X_bad <- raster::stack(intercept, chm_bad)
  names(X_bad) <- c("intercept", "chm")

  expect_error(
    cos_make_areal_blocks(
      U_sf = U_sf,
      X_rast = X_bad,
      id_col = "U_id",
      verbose = FALSE
    ),
    "Incomplete X_rast values"
  )

  expect_warning(
    blocks_drop <- cos_make_areal_blocks(
      U_sf = U_sf,
      X_rast = X_bad,
      id_col = "U_id",
      missing = "drop",
      verbose = FALSE
    ),
    "rows sum to less than 1"
  )
  expect_s3_class(blocks_drop, "cos_areal_blocks")
  expect_true(any(blocks_drop$row_sums < 1))
  expect_equal(blocks_drop$missing, "drop")

  blocks_renorm <- cos_make_areal_blocks(
    U_sf = U_sf,
    X_rast = X_bad,
    id_col = "U_id",
    missing = "renormalize",
    verbose = FALSE
  )
  expect_s3_class(blocks_renorm, "cos_areal_blocks")
  expect_true(all(abs(blocks_renorm$row_sums - 1) < 1e-8))
  expect_equal(blocks_renorm$missing, "renormalize")
})
