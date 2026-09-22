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
- Replaced the spatial-guidance study with a controlled comparison. Its priors
  are independent of held-out responses and shared parameters use identical
  priors across models. The taper stays fixed as exponential range changes.
  Added a no-spatial-effect case, variance-prior sensitivity, multiple chains,
  diagnostic gates, paired simulation uncertainty, and archived results.
  The article distinguishes latent and observed targets and reports interval
  undercoverage as well as prediction gains.
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
  The sole note concerns installed package size, from vignette assets and
  archived study data.
- The README example and the binomial and negative-binomial workflow chunks
  ran, including chunks marked `eval: false`. The simulations reproduced the
  documented 15 binary successes and total count of 147.
- The full introductory Gaussian workflow also ran successfully: three
  chains, recovery, prediction over all 65,792 raster cells and 49 polygons,
  and ten-fold cross-validation for both latent and observed targets.
- The redesigned spatial study ran in full: ten replicates, six generating
  cases, three prior scales, and five folds for both models. All 1,800 final
  fold fits passed the four-chain diagnostic gates, with 10,000–40,000
  iterations per chain. Maximum R-hat was 1.009859; minimum bulk and tail ESS
  were 517.1 and 402.8. Failed shorter attempts were rerun and recorded.
- An independent audit reconstructed dataset-level scores from all 35,280
  held-out prediction records and checked common observations, latent truth,
  and folds across models and priors. Every comparison passed. New tables and
  figures use these archived results; the previous figures were removed.
- Regression tests check response-independent shared priors, paired score
  aggregation, compatible checkpoint settings, disagreement or stuck chains,
  and refusal to plot results with failed diagnostics.
- The pkgdown build and configuration checks passed. All 52 generated HTML
  pages were scanned for local links, anchors, images, scripts, stylesheets,
  and CSS resources; none were missing.
- Chrome checks covered the homepage, all five articles, reference and article
  indexes, equations, search, and mobile navigation. No broken images,
  equation errors, horizontal page overflow, or browser errors were found.

The GitHub Actions workflow has successfully built and deployed the site.
