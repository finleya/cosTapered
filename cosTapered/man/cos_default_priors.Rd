\name{cos_default_priors}
\alias{cos_default_priors}
\title{Construct Default Prior Objects}
\usage{
cos_default_priors(prep, beta_mu = NULL, beta_sd = 10,
  tau_B_shape = 2, tau_B_scale = 50^2, sigma_shape = 2,
  sigma_scale = 10, phi_lower, phi_upper)
}
\arguments{
\item{prep}{A \code{"cos_prep"} object from \code{\link{cos_prepare}}.}
\item{beta_mu}{Prior mean vector for regression coefficients. If unnamed,
values are matched to covariates by position. If \code{NULL}, zeros are used.}
\item{beta_sd}{Prior standard deviation vector for regression coefficients.}
\item{tau_B_shape, tau_B_scale}{Shape and scale for the inverse-gamma prior
on \code{tau_B_sq}, the average observed-support nugget variance.}
\item{sigma_shape, sigma_scale}{Shape and scale for the inverse-gamma prior
on \code{sigma_sq}. Ignored by \code{\link{cos_fit}} for non-spatial
preparation objects.}
\item{phi_lower, phi_upper}{Positive lower and upper bounds for the uniform
prior on \code{phi}. Required for spatial preparation objects and ignored for
non-spatial preparation objects.}
}
\description{
Builds a simple prior specification for the supported model: proper normal
regression prior, inverse-gamma variance priors, and, for spatial fits,
bounded uniform \code{phi} prior.
}
\value{
A list with class \code{"cos_priors"}.
}
\examples{
\dontrun{
library(sf)
library(raster)

chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"
B_sf <- st_read(system.file("extdata", "example_B.gpkg",
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

priors <- cos_default_priors(
  prep = prep,
  beta_mu = c(intercept = 40, chm = 1),
  beta_sd = c(intercept = 20, chm = 2),
  tau_B_shape = 2,
  tau_B_scale = truth$tau_B_sq_true,
  sigma_shape = 2,
  sigma_scale = truth$sigma_sq_true,
  phi_lower = 3 / 600,
  phi_upper = 3 / 80
)

priors
}
}
