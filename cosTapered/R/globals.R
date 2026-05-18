if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("chain", "iter", "value"))
}

cos_same_crs <- function(x, y) {
  crs_x <- sf::st_crs(x)
  crs_y <- sf::st_crs(y)

  if (is.na(crs_x) || is.na(crs_y)) {
    return(FALSE)
  }
  if (isTRUE(crs_x == crs_y)) {
    return(TRUE)
  }

  proj_x <- crs_x$proj4string
  proj_y <- crs_y$proj4string
  if (is.null(proj_x) || is.null(proj_y) || is.na(proj_x) || is.na(proj_y)) {
    return(FALSE)
  }

  identical(proj_x, proj_y)
}
