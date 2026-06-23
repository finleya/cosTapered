# cosTapered

`cosTapered` fits Bayesian change-of-support spatial models for forest
inventory and remote-sensing problems where responses, predictors, and
prediction targets live on different spatial supports.

The motivating case is common in forest applications: LiDAR-derived predictors
are available on a fine raster grid, field responses are observed on larger
fixed-area plots, and final summaries are needed for irregular stands or other
management units. Rather than forcing everything onto one arbitrary support
before modeling, `cosTapered` defines a latent response surface on the fine
predictor grid and links plots and stands to that surface through
area-weighted averages.

![Example CHM raster, stand polygons, observed plots, and simulated plot response](README_files/example-supports.png)

Plots and stands are treated symmetrically: both are averages of the same
fine-scale process over their footprints. The main Gaussian response model
uses a tapered exponential covariance for spatial residual variation so that
large raster and polygon-support problems remain computationally manageable.
An experimental binomial logit response path is also available for binary or
binomial observed-support responses, along with a fixed-size
negative-binomial path for count responses.

![Example fine-support and stand-support latent predictions](README_files/fine-prediction.png)

## Features

- Prepare raster covariates and observed-support polygons with
  `cos_prepare()`.
- Fit Gaussian, experimental binomial, or fixed-size negative-binomial
  change-of-support models with `cos_fit()`.
- Recover observed-support latent quantities with `cos_recover_B()`.
- Predict fine-support or areal-support targets with `cos_predict_fine()` and
  `cos_predict_areal()`.
- Choose latent or observed prediction targets with `target = "latent"` or
  `target = "observed"`.
- Choose how unobserved spatial prediction uncertainty is handled with
  `spatial_uncertainty = "conditional_mean"`, `"marginal"`, or `"joint"`.
- Run observed-support K-fold validation with `cos_cv_observed()`.

## Installation

`cosTapered` requires R 4.1 or later. The main package dependencies are
installed by `remotes::install_github()`, including `sf`, `raster`,
`exactextractr`, `Matrix`, and `BayesLogit`. `BayesLogit` is required because
the experimental binomial and negative-binomial response samplers use
Polya-Gamma random variates through its C/C++ interface. `exactextractr` is
required to compute the exact fraction of each raster cell covered by plot and
areal prediction polygons, which defines the change-of-support averaging
weights used by the model.

On Linux, `sf` and `exactextractr` may also require geospatial system
libraries such as GDAL, GEOS, PROJ, and udunits before their R packages can be
compiled.

```r
install.packages("remotes")
remotes::install_github("finleya/cosTapered", subdir = "cosTapered")
```

To build the package vignettes from source, install the suggested vignette
tools and the Quarto command-line tool:

```r
install.packages(c("ggplot2", "knitr", "quarto", "rmarkdown"))
remotes::install_github(
  "finleya/cosTapered",
  subdir = "cosTapered",
  build_vignettes = TRUE
)
```

Quarto itself must also be available on the system path; see
<https://quarto.org/docs/get-started/>. Pre-built vignette HTML/PDF files in
the repository can be viewed without rebuilding them.

## Example

```r
library(cosTapered)
library(raster)
library(sf)

chm <- raster(system.file("extdata", "example_chm.tif",
                          package = "cosTapered"))
names(chm) <- "chm"

B_sf <- st_read(system.file("extdata", "example_B.gpkg",
                            package = "cosTapered"), quiet = TRUE)

intercept <- chm
intercept[] <- ifelse(is.na(getValues(chm)), NA_real_, 1)
names(intercept) <- "intercept"

X_rast <- stack(intercept, chm)
names(X_rast) <- c("intercept", "chm")

prep <- cos_prepare(
  X_rast = X_rast,
  B_sf = B_sf,
  response_col = "y_B",
  gamma = 350,
  n_threads = 1,
  verbose = FALSE
)

summary(prep)

truth <- readRDS(system.file("extdata", "example_truth.rds",
                             package = "cosTapered"))

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
  starting = c(
    tau_B_sq = truth$tau_B_sq_true,
    sigma_sq = truth$sigma_sq_true,
    phi = truth$phi_true
  ),
  tuning = c(log_tau_B_sq = 0.25, log_sigma_sq = 0.25, z_phi = 0.25),
  n_batch = 50,
  batch_length = 10,
  seed = 1,
  verbose = FALSE
)

summary(fit, burn_in = 25)
```

See the introductory package vignette for a complete Gaussian workflow:

```r
vignette("cosTapered-introduction", package = "cosTapered")
```

The package includes complementary documents for different levels of detail:

- `cosTapered-introduction`: basic Gaussian workflow.
- `cosTapered-model`: model parameterization, computation, support targets,
  and prediction uncertainty.
- `cosTapered-spatial-guidance`: practical guidance on when spatial COS helps.
- `cosTapered-binomial`: experimental Polya-Gamma binomial response workflow.
- `cosTapered-negative-binomial`: experimental fixed-size negative-binomial
  count response workflow with offsets and areal count summaries.
- `cosTapered.pdf`: TeX/PDF model and software details with posterior
  calculations.

```r
vignette("cosTapered-model", package = "cosTapered")
vignette("cosTapered-spatial-guidance", package = "cosTapered")
vignette("cosTapered-binomial", package = "cosTapered")
vignette("cosTapered-negative-binomial", package = "cosTapered")
system.file("doc", "cosTapered.pdf", package = "cosTapered")
```

## Model Scope

The current package supports the standard tapered exponential covariance and a
proper normal prior for regression coefficients. Users supply the raster
design matrix they want the model to use; `cosTapered` does not automatically
add intercepts, scale covariates, or buffer points.

Prediction functions can return either the underlying support average or the
observed response on that support. Use `target = "latent"` for the underlying
mean response, denoted eta in the model, and `target = "observed"` when the
prediction should include support-level nugget variation.

The main prediction choices are:

```mermaid
flowchart TD
  A[Start with fitted model and recovered B effects] --> B{Prediction support}
  B --> C[Fine cells: cos_predict_fine]
  B --> D[Areal units: cos_predict_areal]

  C --> E{Response family}
  D --> F{Response family}

  E --> G[Gaussian]
  E --> H[Binomial]
  E --> NB[Negative-binomial]
  F --> I[Gaussian]
  F --> J[Binomial]
  F --> NBU[Negative-binomial]

  G --> K{target}
  I --> L{target}
  H --> M[latent link eta and probability p]
  J --> N[latent link eta and probability p]
  NB --> NBF[latent link eta and mean count mu]
  NBU --> NBA[latent link eta and mean count mu]

  K --> O[latent eta]
  K --> P[observed y = eta + nugget]
  L --> Q[latent eta]
  L --> R[observed y = eta + support nugget]
  NBF --> NBY[observed count y via NB sampling]
  NBA --> NBU_Y[observed count y via NB sampling]

  O --> S{spatial_uncertainty}
  P --> S
  Q --> T{spatial_uncertainty}
  R --> T
  M --> U{spatial_uncertainty}
  N --> V{spatial_uncertainty}

  S --> W[conditional_mean: plug-in surface]
  S --> X[marginal: per-cell uncertainty]
  S --> Y[joint: joint fine-surface samples]
  T --> Z[conditional_mean: plug-in unit summaries]
  T --> AA[marginal: per-unit uncertainty]
  U --> AB[conditional_mean: plug-in probability surface]
  U --> AD[marginal: per-cell probability uncertainty]
  V --> AC[conditional_mean: plug-in areal probabilities]
  V --> AE[marginal: per-unit probability uncertainty]
  NBF --> NBFM[marginal: per-cell mean-count uncertainty]
  NBA --> NBAM[marginal: per-unit mean-count uncertainty]
```

For Gaussian fine-support prediction, `spatial_uncertainty = "joint"` draws
from the joint conditional Gaussian surface and can be expensive for large
rasters. `spatial_uncertainty = "marginal"` draws independent marginal
conditional values for each prediction cell or areal unit. Use the marginal
option for per-cell, per-stand, or per-polygon posterior means, SDs, and
intervals when units are interpreted separately. For areal prediction, samples
are marginal by unit and should not be used for totals, contrasts, or rankings
that require cross-unit covariance. `spatial_uncertainty = "conditional_mean"`
is a plug-in fitted-surface summary and does not represent full prediction
uncertainty.

The binomial and fixed-size negative-binomial response models use Polya-Gamma
augmentation through `BayesLogit`. Binomial prediction reports link-scale
summaries and response probability summaries. Negative-binomial prediction
uses a fixed `size` parameter, supports offsets such as log exposure, reports
link-scale and mean-count summaries with `target = "latent"`, and can
simulate observed count predictions with `target = "observed"`.

## Package Source

The R package source is in `cosTapered/`.

## Reference

Zhang, L., Finley, A. O., Nothdurft, A., and Banerjee, S. (2024). Bayesian
modeling of incompatible spatial data: A case study involving Post-Adrian storm
forest damage assessment. *International Journal of Applied Earth Observation
and Geoinformation*, 135, 104224. <https://doi.org/10.1016/j.jag.2024.104224>
