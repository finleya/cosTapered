\name{cos_prepare}
\alias{cos_prepare}
\title{Prepare Observed-Support Objects for a Tapered COS Model}
\usage{
cos_prepare(X_rast, B_sf, response_col = NULL, y_B = NULL, gamma,
  spatial = TRUE, taper_code = 1L, n_threads = 1L,
  missing = c("error", "drop", "renormalize"), row_sum_tol = 1e-5,
  verbose = TRUE)
}
\arguments{
\item{X_rast}{A \pkg{raster} \code{RasterLayer}, \code{RasterStack}, or
\code{RasterBrick} containing prepared fine-support covariates. Layer names
must be valid and unique.}
\item{B_sf}{An \pkg{sf} object containing observed-support polygons in the
same projected CRS as \code{X_rast}.}
\item{response_col}{Optional column name in \code{B_sf} containing observed
responses.}
\item{y_B}{Optional numeric vector of observed responses. If supplied,
\code{response_col} is not used.}
\item{gamma}{Positive taper distance in CRS units. Required when
\code{spatial = TRUE}.}
\item{spatial}{Logical; if \code{FALSE}, prepare the non-spatial COS
submodel without covariance-pair precomputation.}
\item{taper_code}{Integer taper code: \code{1L} for Wendland or \code{2L}
for spherical.}
\item{n_threads}{Number of OpenMP threads used by compiled covariance code.}
\item{missing}{How to handle incomplete raster covariate cells that overlap
observed supports. \code{"error"} stops if retained support weights fall below
one. \code{"drop"} keeps the covered-cell weights without rescaling.
\code{"renormalize"} rescales retained weights to sum to one within each
observed support.}
\item{row_sum_tol}{Tolerance used when warning about support-weight row
sums.}
\item{verbose}{Logical; if \code{TRUE}, print progress messages from exact
extraction.}
}
\description{
Constructs the observed-support design matrix, support-weight matrix,
nugget matrix, and, when \code{spatial = TRUE}, tapered covariance
precomputations needed by \code{\link{cos_fit}}.
}
\details{
By default, incomplete raster covariate cells are treated as an error when
they would change the observed-support averaging weights. This protects the
model target: support averages are defined over the supplied polygons, not
only over the subset of cells with complete covariates. Use
\code{missing = "drop"} or \code{missing = "renormalize"} only when that
behavior matches the intended analysis.
}
\value{
A list with class \code{"cos_prep"} containing observed supports,
fine-support covariates, \code{H_BA}, \code{D_h}, compressed support
weights, optional precomputed covariance pairs, model settings, and
\code{timing}, a named numeric vector with user, system, and elapsed seconds.
}
\references{
Zhang, L., Finley, A. O., Nothdurft, A., and Banerjee, S. (2024).
Bayesian modeling of incompatible spatial data: A case study involving
Post-Adrian storm forest damage assessment. \emph{International Journal
of Applied Earth Observation and Geoinformation}, 135, 104224.
\doi{10.1016/j.jag.2024.104224}
}
\examples{
\dontrun{
library(sf)
library(raster)
library(ggplot2)

chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"
B_sf <- st_read(system.file("extdata", "example_B.gpkg",
                            package = "cosTapered"), quiet = TRUE)
U_sf <- st_read(system.file("extdata", "example_U.gpkg",
                            package = "cosTapered"), quiet = TRUE)
forest_sf <- st_read(system.file("extdata", "example_forest.gpkg",
                                 package = "cosTapered"), quiet = TRUE)
truth <- readRDS(system.file("extdata", "example_truth.rds",
                             package = "cosTapered"))

intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"
X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

prep <- cos_prepare(
  X_rast = X_rast,
  B_sf = B_sf,
  response_col = "y_B",
  gamma = truth$gamma,
  taper_code = truth$taper_code,
  verbose = FALSE
)
summary(prep)

# Plot the support geometry and the simulated observed response.
chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)
names(chm_df) <- c("x", "y", "chm")

ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = chm),
            width = xres(chm), height = yres(chm)) +
  geom_sf(data = U_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = B_sf, fill = NA, color = "#c1121f", linewidth = 0.45) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Fine-support covariate and support polygons",
       fill = "CHM", x = NULL, y = NULL) +
  theme_bw()

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey55",
          linewidth = 0.4) +
  geom_sf(data = B_sf, aes(fill = y_B), color = "black", linewidth = 0.25) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf(expand = FALSE) +
  labs(title = "Observed-support response", fill = "y_B",
       x = NULL, y = NULL) +
  theme_bw()
}
}
