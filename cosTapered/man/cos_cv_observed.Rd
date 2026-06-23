\name{cos_cv_observed}
\alias{cos_cv_observed}
\title{K-Fold Observed-Support Prediction Diagnostics}
\usage{
cos_cv_observed(fit, X_rast, B_sf = NULL, response_col = NULL,
  y_B = NULL, target = "observed", latent_col = NULL, k = 10L,
  fold_id = NULL, missing = c("error", "drop", "renormalize"),
  row_sum_tol = 1e-5, fit_args = list(), recover_args = list(),
  seed = NULL, keep_fits = FALSE, verbose = TRUE)
}
\arguments{
\item{fit}{A \code{"cos_fit"} object. The fitted object is used as the model
template: spatial/non-spatial status, taper settings, priors, and sampler
defaults are inherited from it.}
\item{X_rast}{A \pkg{raster} \code{RasterLayer}, \code{RasterStack}, or
\code{RasterBrick} containing fine-support covariates.}
\item{B_sf}{Optional \pkg{sf} polygon object containing observed-support
polygons. If \code{NULL}, \code{fit$prep$B_sf} is used.}
\item{response_col}{Optional name of the response column in \code{B_sf}.}
\item{y_B}{Optional numeric response vector. If supplied, this is used instead
of \code{response_col}. If neither \code{y_B} nor \code{response_col} is
supplied, \code{fit$prep$y_B} is used.}
\item{target}{Prediction target or targets to score. \code{"observed"} scores
support responses with nugget variation included against held-out \code{y_B}.
\code{"latent"} scores latent means against \code{B_sf[[latent_col]]}. Supply
\code{c("observed", "latent")} to score both targets from the same fold fits.}
\item{latent_col}{Name of the latent truth column in \code{B_sf}. Required
when \code{target} includes \code{"latent"}. This is mainly for simulation
studies where latent \code{eta_B} is known.}
\item{k}{Number of folds.}
\item{fold_id}{Optional integer vector assigning each observed unit to a fold
in \code{1:k}. If \code{NULL}, folds are generated at random.}
\item{missing}{How to handle incomplete raster covariate cells in training or
held-out supports. Passed to \code{\link{cos_prepare}} and
\code{\link{cos_make_areal_blocks}}.}
\item{row_sum_tol}{Tolerance used when checking support-weight row sums.}
\item{fit_args}{Named list of additional arguments passed to
\code{\link{cos_fit}}. Values supplied here override sampler defaults
inherited from \code{fit}.}
\item{recover_args}{Named list of additional arguments passed to
\code{\link{cos_recover_B}}.}
\item{seed}{Optional seed used for random fold assignment and default
per-fold fit/recovery seeds.}
\item{keep_fits}{Logical; if \code{TRUE}, keep fitted fold models in the
returned object.}
\item{verbose}{Logical; if \code{TRUE}, print fold progress.}
}
\description{
Runs K-fold validation for observed-support prediction targets. Each fold
refits the model to the training observed supports and predicts the held-out
supports using \code{\link{cos_predict_areal}}. By default, it scores
observed-support responses with nugget variation included. In simulation
studies, it can also score latent held-out means when latent truth is supplied.
}
\details{
Held-out observed supports are treated as areal prediction units. The function
uses \code{cos_predict_areal(spatial_uncertainty = "marginal",
keep_samples = TRUE)} to draw latent areal means. For
\code{target = "observed"}, it then adds independent nugget draws with
variance \code{tau_sq * sum(h^2)} for each held-out support, where \code{h}
contains the held-out support's fine-cell area weights.

For \code{target = "latent"}, the latent areal mean samples are scored against
\code{B_sf[[latent_col]]}. This target is useful in simulation studies for
separating recovery of the underlying support mean from prediction of the
support response with nugget variation included.

The CRPS is computed from posterior predictive samples using the empirical
scoring rule
\deqn{CRPS(F, y) = E_F |Y - y| - \frac{1}{2} E_F |Y - Y'|.}

This function is intended for model assessment and simulation studies rather
than lightweight examples. For simulation-based checks, generate data with
\code{\link{cos_simulate_data}} and pass the simulated \code{B_sf} and
\code{X_rast} to this function after fitting a model to the simulated data.
}
\value{
A list with class \code{"cos_cv_observed"} containing:
\itemize{
  \item \code{summary}: overall RMSPE, CRPS, coverage, interval width, and
    error summaries. A \code{target} column is included.
  \item \code{fold_summary}: fold-level diagnostics and selected parameter
    medians, by target.
  \item \code{predictions}: held-out predictions, intervals, CRPS, errors,
    and coverage indicators, by target.
  \item \code{fold_id}: fold assignments.
  \item \code{fits}: optional fitted fold models when \code{keep_fits = TRUE}.
  \item \code{timing}: named numeric vector with user, system, and elapsed
    seconds.
}
}
\examples{
\dontrun{
cv <- cos_cv_observed(
  fit = fit,
  X_rast = X_rast,
  k = 10,
  fit_args = list(n_batch = 50, batch_length = 8),
  recover_args = list(burn_in = 100, n_samples = 60),
  seed = 1
)

summary(cv)
cv$predictions
}
}
