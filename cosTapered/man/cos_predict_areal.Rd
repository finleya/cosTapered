\name{cos_predict_areal}
\alias{cos_predict_areal}
\title{Predict Areal-Support Latent Quantities}
\usage{
cos_predict_areal(fit, rec_B, U_blocks, method = c("mean", "sample"),
  target = c("latent", "observed"), keep_samples = FALSE,
  n_threads = NULL, verbose = TRUE)
}
\arguments{
\item{fit}{A \code{"cos_fit"} object.}
\item{rec_B}{A \code{"cos_recovery_B"} object.}
\item{U_blocks}{A \code{"cos_areal_blocks"} object from
\code{\link{cos_make_areal_blocks}}.}
\item{method}{Use conditional means or sample one-dimensional conditional
Gaussian effects for each areal unit.}
\item{target}{Prediction target. \code{"latent"} returns latent areal means
\code{eta_U = X_U beta + omega_U}. \code{"observed"} also adds independent
support-level nugget draws with variance \code{tau_sq * sum(h^2)} for each
prediction support.}
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
observed/support response rather than the underlying latent mean.
}
\value{
A list with class \code{"cos_prediction_areal"} containing areal prediction
summaries and, optionally, posterior samples. Existing \code{eta_*} summaries
are always returned. When \code{target = "observed"}, \code{y_*} summaries and,
if requested, \code{y_samples} are also returned.
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
pred_U <- cos_predict_areal(fit, rec_B, U_blocks, method = "mean",
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
