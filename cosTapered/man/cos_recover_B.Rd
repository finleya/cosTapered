\name{cos_recover_B}
\alias{cos_recover_B}
\title{Recover Observed-Support Latent Quantities}
\usage{
cos_recover_B(fit, burn_in = 1, thin = 1L, n_samples = NULL,
  seed = NULL, n_threads = NULL, verbose = TRUE)
}
\arguments{
\item{fit}{A \code{"cos_fit"} object from \code{\link{cos_fit}}.}
\item{burn_in}{First retained iteration index.}
\item{thin}{Within-chain thinning interval.}
\item{n_samples}{Optional maximum number of retained posterior samples.}
\item{seed}{Optional random seed used when subsampling retained draws.}
\item{n_threads}{Optional number of OpenMP threads used by compiled
covariance code. If \code{NULL}, uses \code{fit$prep$n_threads}.}
\item{verbose}{Logical; if \code{TRUE}, print progress messages.}
}
\description{
Given retained posterior samples of covariance parameters, draws
regression coefficients and observed-support spatial effects from their
Gaussian full conditional distribution.
}
\value{
A list with class \code{"cos_recovery_B"} containing posterior samples of
\code{beta}, \code{omega_B}, \code{eta_B}, and \code{timing}, a named numeric
vector with user, system, and elapsed seconds.
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
rec_sum <- summary(rec_B, include_eta_B = TRUE)
rec_sum

# Compare fitted observed-support latent means with known simulated truth.
B_plot <- B_sf
B_plot$eta_B_mean <- rec_sum$eta_B_summary$mean[
  match(B_plot$B_id, rec_sum$eta_B_summary$B_id)
]
B_plot$eta_B_error <- B_plot$eta_B_mean - B_plot$eta_B_true

ggplot(B_plot, aes(x = eta_B_true, y = eta_B_mean)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.4) +
  geom_point(size = 2) +
  coord_equal() +
  labs(title = "Recovered observed-support latent means",
       x = "Simulated eta_B truth",
       y = "Posterior mean") +
  theme_bw()

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey60",
          linewidth = 0.35) +
  geom_sf(data = B_plot, aes(fill = eta_B_error),
          color = "black", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c1121f") +
  coord_sf(expand = FALSE) +
  labs(title = "Fitted eta_B posterior mean minus truth",
       fill = "error", x = NULL, y = NULL) +
  theme_bw()
}
}
