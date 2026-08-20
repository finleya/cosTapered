\name{cos_fit}
\alias{cos_fit}
\alias{cos_summary_theta}
\alias{cos_plot_trace}
\title{Fit a COS Model}
\usage{
cos_fit(prep, priors,
  family = c("gaussian", "binomial", "negative_binomial", "negbin", "nb"),
  trials = NULL, size = NULL, offset = NULL, n_chains = 1L, starting = NULL,
  tuning = c(log_tau_B_sq = 1, log_sigma_sq = 1, z_phi = 1),
  n_batch = 1000, batch_length = 25, accept_rate = 0.234,
  seed = NULL, report = 100, n_threads = NULL, verbose = TRUE)

cos_summary_theta(fit, burn_in = 1,
  pars = c("tau_B_sq", "tau_sq", "sigma_sq", "phi", "eff_range"),
  probs = c(0.5, 0.025, 0.975), digits = 3)

cos_plot_trace(fit, burn_in = 1,
  pars = c("tau_B_sq", "sigma_sq", "phi", "eff_range"))
}
\arguments{
\item{prep}{A \code{"cos_prep"} object.}
\item{priors}{A \code{"cos_priors"} object.}
\item{family}{Response family. The default \code{"gaussian"} uses the
collapsed Gaussian COS sampler. \code{"binomial"} uses a separate
Polya-Gamma sampler with explicit \code{beta} and observed-support spatial
effect draws. \code{"negative_binomial"} uses the same Polya-Gamma path for
fixed-size count models; aliases \code{"negbin"} and \code{"nb"} are
accepted.}
\item{trials}{Binomial trial counts for \code{family = "binomial"}. Supply
a positive integer vector of length one or \code{length(prep$y_B)}, or the
name of a column in \code{prep$B_sf}.}
\item{size}{Fixed negative-binomial size parameter for
\code{family = "negative_binomial"}.}
\item{offset}{Optional linear-predictor offset for Polya-Gamma response
families. Supply \code{NULL}, a finite numeric vector of length one or
\code{length(prep$y_B)}, or the name of a column in \code{prep$B_sf}. For
count models this is commonly a log exposure.}
\item{n_chains}{Number of independent MCMC chains.}
\item{starting}{Optional starting values. Spatial fits use
\code{tau_B_sq}, \code{sigma_sq}, and \code{phi} for Gaussian fits. Spatial
Polya-Gamma fits use \code{sigma_sq} and \code{phi}; optional beta starting
values can be supplied in columns named like the fitted covariates.}
\item{tuning}{Named proposal tuning vector on the transformed scale.
Spatial fits use \code{log_tau_B_sq}, \code{log_sigma_sq}, and
\code{z_phi} for Gaussian fits. Spatial Polya-Gamma fits use
\code{log_sigma_sq} and \code{z_phi}. Non-spatial Polya-Gamma fits do not use
a Metropolis tuning vector.}
\item{n_batch}{Number of adaptive Metropolis batches.}
\item{batch_length}{Iterations per batch.}
\item{accept_rate}{Target acceptance rate.}
\item{seed}{Optional random seed.}
\item{report}{Batch reporting interval.}
\item{n_threads}{Optional number of OpenMP threads used by compiled
covariance code. If \code{NULL}, uses \code{prep$n_threads}.}
\item{verbose}{Logical; if \code{TRUE}, print sampler progress.}
\item{fit}{A \code{"cos_fit"} object.}
\item{burn_in}{First retained iteration index.}
\item{pars}{Parameter names to summarize or plot.}
\item{probs}{Quantile probabilities.}
\item{digits}{Number of digits for rounded summaries.}
}
\description{
Fits a COS model. For \code{family = "gaussian"}, this fits the collapsed
posterior. Spatial preparation objects fit
\code{tau_B_sq}, \code{sigma_sq}, and \code{phi}. Non-spatial preparation
objects fit the COS submodel without spatial random effects and sample only
\code{tau_B_sq}. The fine-support nugget variance \code{tau_sq} is recovered
from \code{tau_B_sq / mean(diag(D_h))}.

For \code{family = "binomial"}, the sampler uses Polya-Gamma augmentation on
the observed support, with linear predictor
\eqn{\eta_l = offset_l + x_B(B_l)^\top\beta + \omega_B(B_l)} for observed
support \eqn{B_l}. Stacked over supports, this is
\code{eta_B = offset + X_B beta + omega_B}. This first implementation keeps
\code{beta} and \code{omega_B} explicit in the Markov state. There is no
Gaussian nugget variance \code{tau_sq} in the binomial likelihood.

For \code{family = "negative_binomial"}, the sampler uses a fixed-size
negative-binomial likelihood with mean \eqn{\mu_l = \exp(\eta_l)} on observed
support \eqn{B_l}. Stacked over supports,
\code{eta_B = offset + X_B beta + omega_B}. The size parameter is fixed by
the \code{size} argument. There is no Gaussian nugget variance
\code{tau_sq} in this likelihood.
}
\value{
\code{cos_fit()} returns a list with class \code{"cos_fit"} containing
posterior samples and \code{timing}, a named numeric vector with user,
system, and elapsed seconds. The convenience function
\code{cos_summary_theta()} returns a parameter summary table. \code{cos_plot_trace()}
returns a \pkg{ggplot2} object.

For Polya-Gamma fits, the returned object also contains
\code{beta_samples}, \code{omega_B_samples}, and \code{eta_B_samples}. These
are used directly by \code{cos_recover_B()}.
}
\examples{
\dontrun{
library(sf)
library(raster)
library(ggplot2)

# Read the bundled example data.
chm <- raster(system.file("extdata", "example_chm.tif", package = "cosTapered"))
names(chm) <- "chm"
B_sf <- st_read(system.file("extdata", "example_B.gpkg",
                            package = "cosTapered"), quiet = TRUE)
truth <- readRDS(system.file("extdata", "example_truth.rds",
                             package = "cosTapered"))

# Users prepare the raster design matrix. Here we add an intercept layer.
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

fit <- cos_fit(
  prep = prep,
  priors = priors,
  n_chains = 1,
  starting = data.frame(
    tau_B_sq = truth$tau_B_sq_true,
    sigma_sq = truth$sigma_sq_true,
    phi = truth$phi_true
  ),
  tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
  n_batch = 50,
  batch_length = 8,
  seed = 1,
  verbose = FALSE
)

summary(fit, burn_in = 100)
plot(fit, burn_in = 100)

# Compare posterior parameter samples with the known simulation values.
theta <- fit$theta_samples
theta <- theta[theta$iter >= 100, ]
theta_long <- reshape(
  theta[, c("chain", "iter", "tau_B_sq", "sigma_sq", "phi", "eff_range")],
  varying = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
  v.names = "value",
  timevar = "parameter",
  times = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
  direction = "long"
)
truth_df <- data.frame(
  parameter = c("tau_B_sq", "sigma_sq", "phi", "eff_range"),
  truth = c(truth$tau_B_sq_true, truth$sigma_sq_true,
            truth$phi_true, truth$eff_range_true)
)

ggplot(theta_long, aes(x = value)) +
  geom_density(fill = "grey80", color = "grey25", linewidth = 0.35) +
  geom_vline(data = truth_df, aes(xintercept = truth),
             color = "#c1121f", linewidth = 0.7) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +
  labs(title = "Posterior distributions and known simulation values",
       x = NULL, y = NULL) +
  theme_bw()
}
}
