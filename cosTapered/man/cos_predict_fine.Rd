\name{cos_predict_fine}
\alias{cos_predict_fine}
\title{Predict Fine-Support Latent Quantities}
\usage{
cos_predict_fine(fit, rec_B, pred_coords = NULL, X_pred = NULL,
  offset_pred = NULL,
  spatial_uncertainty = c("conditional_mean", "marginal", "joint"),
  target = c("latent", "observed"), keep_samples = FALSE,
  n_threads = NULL, verbose = TRUE)
}
\arguments{
\item{fit}{A \code{"cos_fit"} object.}
\item{rec_B}{A \code{"cos_recovery_B"} object.}
\item{pred_coords}{Optional matrix of prediction coordinates with columns
\code{x} and \code{y}. If omitted, the original fine support in \code{fit}
is used.}
\item{X_pred}{Optional matrix of prediction covariates with column names
matching the fitted covariates.}
\item{offset_pred}{Optional prediction offset for Polya-Gamma response
families. Supply \code{NULL} or a finite numeric vector of length one or
\code{nrow(X_pred)}. For negative-binomial prediction this is commonly a log
exposure.}
\item{spatial_uncertainty}{How to handle unobserved fine-support spatial
effects. \code{"conditional_mean"} is faster but omits conditional spatial
prediction uncertainty. \code{"marginal"} draws independent marginal
conditional effects for each prediction cell. \code{"joint"} draws the joint
conditional spatial effect across fine prediction cells and can require a
large dense covariance matrix.}
\item{target}{Prediction target. For Gaussian fits, \code{"latent"} returns
latent fine-scale means \code{eta_A = X_A beta + omega_A}, while
\code{"observed"} also adds independent fine-scale nugget draws with variance
\code{tau_sq}. For negative-binomial fits, \code{"latent"} summarizes the log
mean \code{eta} and mean count \code{mu = exp(eta)}, while
\code{"observed"} also simulates count draws from the fitted
negative-binomial likelihood.}
\item{keep_samples}{Logical; if \code{TRUE}, store posterior sample matrices.}
\item{n_threads}{Optional number of OpenMP threads used by compiled
covariance code. If \code{NULL}, uses \code{fit$prep$n_threads}.}
\item{verbose}{Logical; if \code{TRUE}, print progress messages.}
}
\description{
Predicts fine-support effects and latent means. For spatial fits this uses
the conditional Gaussian distribution of fine-support effects given
observed-support effects. For non-spatial fits the spatial effect is fixed at
zero and predictions use recovered regression coefficients only.
}
\details{
By default the function reports the latent fine-scale mean target. This is the
usual target for mapping the underlying forest condition. Use
\code{target = "observed"} when the desired posterior predictive target is a
fine-support response with finite-cell nugget variation included.

With \code{spatial_uncertainty = "conditional_mean"}, posterior summaries
reflect uncertainty in recovered parameters and observed-support effects, but
they omit conditional uncertainty in unobserved fine-support spatial effects.
Use \code{spatial_uncertainty = "marginal"} for independent per-cell marginal
posterior predictive uncertainty or \code{spatial_uncertainty = "joint"} when
joint fine-surface samples are needed.

With \code{spatial_uncertainty = "joint"}, the conditional covariance can be
positive semidefinite when predicting fine cells constrained by
observed-support effects. In that case the function uses eigen-based sampling
for the valid degenerate Gaussian distribution.

For binomial fits, the function reports link-scale \code{eta_*} summaries and
response probability \code{p_*} summaries. For negative-binomial fits, it
reports link-scale \code{eta_*} summaries and mean-count \code{mu_*}
summaries. When \code{target = "observed"} for a negative-binomial fit, it
also reports simulated count \code{y_*} summaries.
}
\value{
A list with class \code{"cos_prediction_fine"} containing prediction
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

pred <- cos_predict_fine(fit, rec_B,
                         spatial_uncertainty = "conditional_mean",
                         verbose = FALSE)
summary(pred)

fine_df <- pred$summary[, c("x", "y", "eta_mean", "eta_sd")]

ggplot() +
  geom_tile(data = fine_df, aes(x = x, y = y, fill = eta_mean),
            width = xres(chm), height = yres(chm)) +
  geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "white", linewidth = 0.35) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf(expand = FALSE) +
  labs(title = "Fine-support posterior mean",
       fill = "eta mean", x = NULL, y = NULL) +
  theme_bw()

ggplot() +
  geom_tile(data = fine_df, aes(x = x, y = y, fill = eta_sd),
            width = xres(chm), height = yres(chm)) +
  geom_sf(data = B_sf, fill = NA, color = "black", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "white", linewidth = 0.35) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Fine-support posterior SD",
       fill = "eta SD", x = NULL, y = NULL) +
  theme_bw()
}
}
