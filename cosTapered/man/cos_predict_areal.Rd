\name{cos_predict_areal}
\alias{cos_predict_areal}
\title{Predict Areal-Support Quantities}
\usage{
cos_predict_areal(fit, rec_B, U_blocks, offset_U = NULL,
  spatial_uncertainty = c("conditional_mean", "marginal"),
  target = c("latent", "observed"), keep_samples = FALSE,
  n_threads = NULL, verbose = TRUE)
}
\arguments{
\item{fit}{A \code{"cos_fit"} object.}
\item{rec_B}{A \code{"cos_recovery_B"} object.}
\item{U_blocks}{A \code{"cos_areal_blocks"} object from
\code{\link{cos_make_areal_blocks}}.}
\item{offset_U}{Optional prediction offset for Polya-Gamma response
families. Supply \code{NULL}, a finite numeric vector of length one or
\code{length(U_blocks$blocks)}, or the name of a column in
\code{U_blocks$U_sf}. For negative-binomial areal prediction this is commonly
a log exposure for each prediction unit.}
\item{spatial_uncertainty}{How to handle unobserved areal-support spatial
effects. \code{"conditional_mean"} is faster but omits conditional spatial
prediction uncertainty. \code{"marginal"} draws independent marginal
conditional effects for each areal unit, but samples are not joint across
units.}
\item{target}{Prediction target. For Gaussian fits, \code{"latent"} returns
latent areal means \eqn{\eta(U_k)=\sum_i g_{ki}\eta(A_i)}, while
\code{"observed"} also adds independent support-averaged nugget draws with
variance \eqn{\tau^2\sum_i g_{ki}^2} for prediction support \eqn{U_k}. For
negative-binomial fits, \code{"latent"} summarizes the log mean \code{eta}
and mean count \code{mu = exp(eta)}, while \code{"observed"} also simulates
count draws from the fitted negative-binomial likelihood.}
\item{keep_samples}{Logical; if \code{TRUE}, store posterior sample matrices.}
\item{n_threads}{Optional number of OpenMP threads used by compiled
covariance code. If \code{NULL}, uses \code{fit$prep$n_threads}.}
\item{verbose}{Logical; if \code{TRUE}, print progress messages.}
}
\description{
Predicts areal-support effects and latent means directly at the
prediction-unit support. For spatial fits this uses tapered covariance
cross-covariances. For non-spatial fits the spatial effect is fixed at zero
and predictions use recovered regression coefficients only.
}
\details{
By default the function reports the latent areal mean target. This is the
usual target for stand means and inventory summaries. Use
\code{target = "observed"} when the desired posterior predictive target is an
areal-support response with support-averaged nugget variation included rather
than the underlying latent mean.

With \code{spatial_uncertainty = "conditional_mean"}, posterior summaries
reflect uncertainty in recovered parameters and observed-support effects, but
they omit conditional uncertainty in unobserved areal spatial effects. With
\code{spatial_uncertainty = "marginal"}, samples include marginal conditional
spatial uncertainty for each prediction unit, but they are not joint samples
across units. Do not use those samples for totals, contrasts, rankings, or
other summaries that require cross-unit posterior covariance.

For binomial fits, the function reports link-scale \code{eta_*} summaries and
response probability \code{p_*} summaries. For negative-binomial fits, it
reports link-scale \code{eta_*} summaries and mean-count \code{mu_*}
summaries. When \code{target = "observed"} for a negative-binomial fit, it
also reports simulated count \code{y_*} summaries. The \code{offset_U}
argument defines the prediction-unit exposure or other known offset; the
package does not infer count exposure from polygon area.
}
\value{
A list with class \code{"cos_prediction_areal"} containing areal prediction
summaries and, optionally, posterior samples. Existing \code{eta_*} summaries
are always returned. When \code{target = "observed"}, \code{y_*} summaries and,
if requested, \code{y_samples} are also returned. Binomial fits return
\code{p_*} summaries and optional \code{p_samples}; negative-binomial fits
return \code{mu_*} summaries and optional \code{mu_samples}. The object also
contains \code{timing}, a named numeric vector with user, system, and elapsed
seconds.
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

prep <- cos_prepare(X_rast, B_sf, response_col = "y_B",
                    gamma = truth$gamma, taper_code = truth$taper_code,
                    verbose = FALSE)
priors <- cos_default_priors(
  prep,
  beta_mu = c(intercept = 40, chm = 1),
  beta_sd = c(intercept = 20, chm = 2),
  tau_B_shape = 2,
  tau_B_scale = truth$tau_B_sq_true,
  sigma_shape = 2,
  sigma_scale = truth$sigma_sq_true,
  phi_lower = 3 / 600,
  phi_upper = 3 / 80
)
fit <- cos_fit(
  prep, priors,
  starting = data.frame(tau_B_sq = truth$tau_B_sq_true,
                        sigma_sq = truth$sigma_sq_true,
                        phi = truth$phi_true),
  tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
  n_batch = 50,
  batch_length = 8,
  seed = 1,
  verbose = FALSE
)
rec_B <- cos_recover_B(fit, burn_in = 100, n_samples = 100,
                       seed = 1, verbose = FALSE)

U_blocks <- cos_make_areal_blocks(U_sf, X_rast, id_col = "U_id",
                                  verbose = FALSE)
pred_U <- cos_predict_areal(fit, rec_B, U_blocks,
                            spatial_uncertainty = "conditional_mean",
                            verbose = FALSE)
summary(pred_U)

U_plot <- pred_U$U_blocks$U_sf
U_plot$eta_mean <- pred_U$summary$eta_mean[
  match(U_plot$U_id, pred_U$summary$U_id)
]
U_plot$eta_sd <- pred_U$summary$eta_sd[
  match(U_plot$U_id, pred_U$summary$U_id)
]

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
          linewidth = 0.35) +
  geom_sf(data = U_plot, aes(fill = eta_mean),
          color = "white", linewidth = 0.25) +
  geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal posterior mean",
       fill = "eta mean", x = NULL, y = NULL) +
  theme_bw()

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
          linewidth = 0.35) +
  geom_sf(data = U_plot, aes(fill = eta_sd),
          color = "white", linewidth = 0.25) +
  geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.2) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Areal posterior SD",
       fill = "eta SD", x = NULL, y = NULL) +
  theme_bw()
}
}
