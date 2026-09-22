# Documentation review

Reviewed against the R and C++ implementation of cosTapered 0.0.2 on
2026-09-22. Package algorithms and exported interfaces were not changed.

## Corrections

- Fine prediction without explicit coordinates uses the complete cells
  intersecting observed polygons. The introductory workflow now explicitly
  supplies all 65,792 cells for its full-raster prediction.
- Prediction quantiles require `keep_samples = TRUE`. The README, introduction,
  model article, PDF, and prediction reference pages now distinguish those
  outputs from the means and SDs returned by default.
- Both prediction functions default to `spatial_uncertainty = "conditional_mean"`.
  The documentation now states that default and explains when to request
  marginal uncertainty.
- `burn_in` is the first retained iteration within each chain. The README and
  Gaussian workflow now explain this convention.
- Non-Gaussian areal predictions transform a support-averaged linear predictor.
  They are not averages of fine-cell probabilities or sums of fine-cell counts.
  Prediction offsets default to zero and must be supplied explicitly.
- The documented cross-validation workflow is Gaussian. The current function
  does not automatically pass a non-Gaussian fit's family, trials, size, and
  offsets to fold refits; the reference page now makes that limitation explicit.
- Removed references to two unavailable `dev/` figure-generation scripts and
  identified saved figures and tables as outputs from earlier runs. Corrected
  the binomial article's claim that its displayed probability SDs increase
  monotonically.
- The spatial study changes the taper distance along with exponential range:
  `gamma = effective_range + 100`. This is now included in the study description.
- Completed the model article's Gaussian recovery formula by defining its
  information vector, and rebuilt the matching model PDF.
- Corrected installation instructions: the repository has a prebuilt PDF;
  HTML vignettes and installed `vignette()` entries require a vignette build.

## Website

The root README is the single source for the homepage. The build script stages
the nested package and homepage assets in a temporary directory, builds all
five Quarto articles and the reference index, includes the model PDF, and
writes `docs/`. Native Quarto image links replace R chunks that could not
locate saved figures when pkgdown changed their execution directory.

The GitHub Actions workflow builds pull requests and deploys `main` through
GitHub Pages. Repository Pages settings must use GitHub Actions as the source;
see [maintenance instructions](README.md).

## Validation

- Package installation and the existing test suite passed.
- `R CMD build` and `R CMD check --no-manual` passed with no errors or warnings.
  The sole note concerns installed package size, primarily vignette assets.
- The README example and the binomial and negative-binomial workflow chunks
  ran, including chunks marked `eval: false`. The simulations reproduced the
  documented 15 binary successes and total count of 147.
- The full introductory Gaussian workflow also ran successfully: three
  chains, recovery, prediction over all 65,792 raster cells and 49 polygons,
  and ten-fold cross-validation for both latent and observed targets.
- A reduced spatial study completed with one replicate, one range setting,
  two folds, and short chains. The full 50-replicate study was not rerun, so its
  saved numerical results and original figures were not independently regenerated.
- The pkgdown build and configuration checks passed. All 52 generated HTML
  pages were scanned for local links, anchors, images, scripts, stylesheets,
  and CSS resources; none were missing.
- Chrome checks covered the homepage, all five articles, reference and article
  indexes, equations, search, and mobile navigation. No broken images,
  equation errors, horizontal page overflow, or browser errors were found.

The GitHub Actions workflow has been checked locally but has not run on GitHub.
