# Spatial-guidance simulation results

Generated on 2026-09-22 with cosTapered 0.0.2. These files replace the older
spatial-guidance results; they are a new study, not a reproduction of its
previous 50-replicate design. The package algorithms used here are unchanged
from repository commit `5ed01f9`.

The installed scripts `spatial-cos-guidance-repeated-study.R` and
`spatial-cos-guidance-plots.R` reproduce the calculations and figures. Their
SHA-256 checksums are in `code-checksums.csv`; the input data MD5 checksums
are in `input-checksums.csv`. The article describes the design and limitations.

The full run used the scripts' default statistical settings, 24 parallel
workers, and one thread per worker. `settings.R` records the configurable
study settings; fixed generating parameters and seeds are in the checksummed
study script. `session-info.txt` records the software environment and `runtime.txt`
records elapsed wall time. All 1,800 final fold fits passed the declared
R-hat and effective-sample-size thresholds. Retry attempts are retained.

## Files

- `scores.csv`: scores over all 49 held-out plots for each simulated dataset,
  fitted model, generating case, prior scale, and prediction target.
- `predictions.csv.gz`: one row per held-out plot, dataset, model, prior scale,
  and target, including truth, prediction means and quantiles, and CRPS.
  R reads it directly with `read.csv("predictions.csv.gz")`.
- `paired-scores.csv`: spatial and nonspatial scores paired within dataset,
  with percentage improvements and win indicators.
- `summary.csv`: dataset-level score means, simulation SEs, and SDs.
- `improvement.csv`: means, simulation SEs, and SDs of paired improvements
  and win indicators. Means of paired percentages are not percentages of
  grand mean scores.
- `diagnostics.csv`: rank-normalized split R-hat, bulk ESS, and tail ESS for
  each final fold fit and sampled covariance parameter or log posterior,
  plus actual iteration counts, first retained iteration, thinning, and
  the number of prediction draws.
- `attempts.csv`: diagnostics and elapsed time for each model's five-fold
  run, including attempts that were rerun with longer chains.
- `design.csv`: plot IDs, centroids in the raster's projected coordinate
  system (meters), clipped areas, and fixed validation fold assignments.

`rep_id` identifies independent simulated datasets within a generating case.
The same dataset is reused for both models and every prior scale. Cases also
reuse simulation seeds, so cases are not independent of one another.
`eff_range` is the nominal exponential range in meters; zero denotes the
no-spatial-effect control, not a zero correlation-range parameter.
`prior_scale` is the inverse-gamma scale shared by the nugget priors and, for
the spatial model, its additional spatial-variance prior.

The archive omits full chain traces and per-case RDS checkpoints to keep the
package small. Running the study creates checkpoints for resuming it. It does
not rerun during package or website builds.
