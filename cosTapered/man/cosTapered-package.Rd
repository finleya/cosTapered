\name{cosTapered-package}
\alias{cosTapered-package}
\alias{cosTapered}
\docType{package}
\title{Tapered Change-of-Support Gaussian Process Model}
\description{
The \pkg{cosTapered} package implements a Bayesian change-of-support
Gaussian process model for areal observations, fine-support raster
covariates, and areal prediction units. The spatial covariance is a
tapered exponential covariance. Observed-support and prediction-support
weights are computed from polygon and raster intersections.
}
\details{
The main workflow is:
\enumerate{
  \item Build observed-support objects with \code{\link{cos_prepare}}.
  \item Specify priors with \code{\link{cos_default_priors}}.
  \item Fit the collapsed posterior with \code{\link{cos_fit}}.
  \item Recover observed-support latent quantities with
    \code{\link{cos_recover_B}}.
  \item Predict fine-support or areal-support quantities with
    \code{\link{cos_predict_fine}} and \code{\link{cos_predict_areal}}.
}
Users provide covariate rasters as prepared inputs. The package does not
add intercepts, scale covariates, or buffer points automatically.
}
\author{
Andrew O. Finley, Lu Zhang, Arne Nothdurft, and Sudipto Banerjee.
}
\references{
Zhang, L., Finley, A. O., Nothdurft, A., and Banerjee, S. (2024).
Bayesian modeling of incompatible spatial data: A case study involving
Post-Adrian storm forest damage assessment. \emph{International Journal
of Applied Earth Observation and Geoinformation}, 135, 104224.
\doi{10.1016/j.jag.2024.104224}
}
\keyword{package}
