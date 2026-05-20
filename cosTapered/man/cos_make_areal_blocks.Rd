\name{cos_make_areal_blocks}
\alias{cos_make_areal_blocks}
\title{Build Areal Prediction Blocks}
\usage{
cos_make_areal_blocks(U_sf, X_rast, id_col = NULL,
  row_sum_tol = 1e-5, verbose = TRUE)
}
\arguments{
\item{U_sf}{An \pkg{sf} object containing areal prediction polygons in the
same projected CRS as \code{X_rast}.}
\item{X_rast}{A \pkg{raster} object with prepared fine-support covariates.}
\item{id_col}{Optional column in \code{U_sf} used as prediction-unit label.}
\item{row_sum_tol}{Tolerance used when warning about support-weight row
sums.}
\item{verbose}{Logical; if \code{TRUE}, print exact extraction progress.}
}
\description{
Computes fine-cell weights, covariate averages, coordinates, and metadata
for each areal prediction polygon.
}
\value{
A list with class \code{"cos_areal_blocks"} containing areal support blocks
and \code{timing}, a named numeric vector with user, system, and elapsed
seconds.
}
\examples{
\dontrun{
library(sf)
library(raster)
library(ggplot2)

chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"
U_sf <- st_read(system.file("extdata", "example_U.gpkg",
                            package = "cosTapered"), quiet = TRUE)
forest_sf <- st_read(system.file("extdata", "example_forest.gpkg",
                                 package = "cosTapered"), quiet = TRUE)

intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"
X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

U_blocks <- cos_make_areal_blocks(U_sf, X_rast, id_col = "U_id",
                                  verbose = FALSE)
summary(U_blocks)

chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)
names(chm_df) <- c("x", "y", "chm")

ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = chm),
            width = xres(chm), height = yres(chm)) +
  geom_sf(data = U_sf, fill = NA, color = "white", linewidth = 0.35) +
  geom_sf(data = forest_sf, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal prediction supports over fine-support covariate",
       fill = "CHM", x = NULL, y = NULL) +
  theme_bw()
}
}
