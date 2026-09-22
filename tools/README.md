# Documentation maintenance

Run commands from the repository root. The R package is in `cosTapered/`.

Install the package's dependencies, its suggested vignette packages, and
`pkgdown >= 2.2.1`. Quarto CLI 1.5 or later must be available. Build the site:

```sh
Rscript tools/build-site.R
```

The script installs the current package in a temporary library, builds from a
temporary source copy, and writes the site to `docs/`. It uses the root
`README.md` and `README_files/` for the homepage, renders all five Quarto
vignettes as articles, includes the model PDF, and checks the finished site.
The temporary library is passed to Quarto's separate R sessions, so a prior
installation of `cosTapered` is not needed. Rendering output is included in
the build log to make failures easier to diagnose.
Generated HTML is ignored by Git. Preview with:

```sh
python3 -m http.server 8000 --directory docs
```

Open <http://localhost:8000/>. Edit the root README, the files in
`cosTapered/vignettes/`, or `cosTapered/_pkgdown.yml` and rebuild.

The pkgdown workflow builds pull requests and publishes pushes to `main`.
For initial publication, set the repository's **Settings > Pages > Build and
deployment > Source** to **GitHub Actions**, then run the `pkgdown` workflow
or push to `main`. The site URL is <https://finleya.github.io/cosTapered/>.

The articles retain the vignettes' execution settings. Most fitting examples
and the long spatial simulation study are not run during site builds. The
spatial-guidance article uses archived, validated study results; other saved
figures are illustrative outputs from earlier runs. Validate displayed
workflow code separately when the API changes; rendering alone does not check
chunks marked `eval: false`.

Build and check the package, including its standalone vignettes, separately:

```sh
R CMD build cosTapered
R CMD check --no-manual cosTapered_0.0.2.tar.gz
```

## Spatial-guidance study

The spatial-guidance article reads its tables from the archived CSV files in
`cosTapered/inst/extdata/spatial-guidance/`. Its figures come from the same
results. Regenerate the study with the installed current package and the
`posterior` and `ggplot2` packages:

```sh
CTV_WORKERS=4 CTV_OUT_DIR=spatial-cos-guidance-output Rscript cosTapered/inst/scripts/spatial-cos-guidance-repeated-study.R
CTV_OUT_DIR=spatial-cos-guidance-output Rscript cosTapered/inst/scripts/spatial-cos-guidance-plots.R
```

This runs ten replicates of six generating cases at three prior scales,
with five-fold validation of both models. Each fit uses four chains and
must meet the recorded R-hat and effective-sample-size thresholds. The script
retains checkpoints and retries with longer chains. An unresolved diagnostic
failure stops figure generation. Parallel workers each use one thread by
default; `CTV_WORKERS=1` also works on Windows.

For a trial, set `CTV_N_REPS=2`, `CTV_EFF_RANGES=0,350`, and
`CTV_PRIOR_SCALES=500` with a separate output directory. A trial is not a
replacement for the archived article results. When refreshing the article,
copy validated scores, summaries, diagnostics, settings, and provenance into
its results directory and the generated figures into `vignettes/figures/`.
Check the written interpretation against those results before publishing.
