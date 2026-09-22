study_env <- new.env(parent = globalenv())
study_script <- testthat::test_path("..", "..", "inst", "scripts",
                                    "spatial-cos-guidance-repeated-study.R")
if (!file.exists(study_script)) {
  study_script <- system.file("scripts", "spatial-cos-guidance-repeated-study.R",
                              package = "cosTapered", mustWork = TRUE)
}
local({
  old <- Sys.getenv("CTV_DEFINE_ONLY", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("CTV_DEFINE_ONLY") else
    Sys.setenv(CTV_DEFINE_ONLY = old))
  Sys.setenv(CTV_DEFINE_ONLY = "TRUE")
  sys.source(study_script, envir = study_env)
  sys.source(file.path(dirname(study_script), "spatial-cos-guidance-plots.R"),
             envir = study_env)
})

test_that("study priors do not learn from validation responses", {
  prep <- structure(list(X_B = cbind(intercept = 1, chm = 1:6),
                         y_B = 1:6, spatial = TRUE), class = "cos_prep")
  before <- study_env$study_priors(prep, 500)
  prep$y_B[6] <- 1e9
  expect_identical(study_env$study_priors(prep, 500), before)
  prep$spatial <- FALSE
  nonspatial <- study_env$study_priors(prep, 500)
  expect_identical(nonspatial$beta, before$beta)
  expect_identical(nonspatial$tau_B_sq, before$tau_B_sq)
})

test_that("study improvements pair replicates before averaging", {
  scores <- data.frame(rep_id = rep(1:2, 2), eff_range = 350,
                       prior_scale = 500, target = "observed",
                       model = rep(c("spatial", "nonspatial"), each = 2),
                       RMSPE = c(1, 12, 2, 10), CRPS = c(1, 12, 2, 10))
  paired <- study_env$study_improvement(scores[c(4, 1, 3, 2), ])
  expect_equal(paired$RMSPE_improvement_pct, c(50, -20))
  summary <- study_env$study_summarize(
    paired, c("target", "eff_range", "prior_scale"), "RMSPE_improvement_pct"
  )
  expect_equal(summary$RMSPE_improvement_pct_mean, 15)
  expect_equal(summary$RMSPE_improvement_pct_se, 35)
  expect_equal(summary$n_reps, 2)
})

test_that("resuming cannot mix study designs or shorten a retry budget", {
  settings <- list(gamma = 1500, prior_scales = c(500, 250, 1000), max_attempts = 3)
  longer <- settings
  longer$max_attempts <- 4
  expect_true(study_env$study_compatible_settings(settings, longer))
  expect_false(study_env$study_compatible_settings(longer, settings))
  changed <- longer
  changed$gamma <- 1000
  expect_false(study_env$study_compatible_settings(settings, changed))
  changed <- longer
  changed$prior_scales <- 250
  expect_false(study_env$study_compatible_settings(settings, changed))
})

test_that("the study diagnostic gate detects chains that disagree or get stuck", {
  skip_if_not_installed("posterior")
  set.seed(29)
  draws <- data.frame(chain = rep(1:4, each = 2000), iter = rep(1:2000, 4),
                      tau_B_sq = exp(rnorm(8000)), lp = rnorm(8000))
  fit <- list(theta_samples = draws, spatial = FALSE)
  settings <- list(n_chains = 4, rhat_limit = 1.01, ess_min = 400)
  good <- study_env$study_diagnostics(list(fit), 1, settings)
  expect_true(all(good$passed))
  fit$theta_samples$tau_B_sq <- exp(rnorm(8000) + 3 * draws$chain)
  bad <- study_env$study_diagnostics(list(fit), 1, settings)
  expect_false(bad$passed[bad$variable == "tau_B_sq"])
  fit$theta_samples$tau_B_sq <- 1
  stuck <- study_env$study_diagnostics(list(fit), 1, settings)
  expect_false(stuck$passed[stuck$variable == "tau_B_sq"])
})

test_that("figures cannot be generated from fits with failed diagnostics", {
  skip_if_not_installed("ggplot2")
  results <- tempfile("study-failed-")
  dir.create(results)
  on.exit(unlink(results, recursive = TRUE))
  write.csv(data.frame(passed = FALSE), file.path(results, "diagnostics.csv"), row.names = FALSE)
  write.csv(data.frame(diagnostics_passed = FALSE), file.path(results, "scores.csv"), row.names = FALSE)
  expect_error(study_env$plot_spatial_guidance_study(results), "Resolve failed diagnostics")
  expect_false(dir.exists(file.path(results, "figures")))
})
