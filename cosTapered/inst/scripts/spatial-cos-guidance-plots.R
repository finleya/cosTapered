# Regenerate article figures from a completed, diagnostically accepted study.
# Set CTV_OUT_DIR to the directory written by the repeated-study script.
plot_spatial_guidance_study <- function(out_dir = Sys.getenv("CTV_OUT_DIR", "spatial-cos-guidance-output")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  read <- function(name) utils::read.csv(file.path(out_dir, paste0(name, ".csv")))
  diagnostics <- read("diagnostics")
  scores <- read("scores")
  if (!all(diagnostics$passed) || !all(scores$diagnostics_passed)) {
    stop("Resolve failed diagnostics before generating publication figures.")
  }
  summary <- read("summary")
  improvement <- read("improvement")
  if (any(summary$n_reps < 2L) || any(improvement$n_reps < 2L)) {
    stop("At least two replicates per case are needed for Monte Carlo intervals.")
  }
  design <- read("design")
  ranges <- sort(unique(summary$eff_range))
  label <- function(x) factor(x, levels = ranges,
                              labels = ifelse(ranges == 0, "None", as.character(ranges)))
  prepare <- function(data) {
    data$case <- label(data$eff_range)
    data$target <- factor(data$target, c("latent", "observed"),
                           c("Latent support mean", "Observed response"))
    if ("model" %in% names(data)) {
      data$model <- factor(data$model, c("spatial", "nonspatial"),
                            c("Spatial COS", "Nonspatial COS"))
    }
    data
  }
  summary <- prepare(summary)
  improvement <- prepare(improvement)
  primary <- summary[summary$prior_scale == 500, ]
  if (!nrow(primary)) stop("Include prior scale 500 for the article's primary comparison.")
  longer <- function(data, metrics, labels) {
    do.call(rbind, lapply(seq_along(metrics), function(i) {
      out <- data
      out$value <- data[[paste0(metrics[i], "_mean")]]
      out$se <- data[[paste0(metrics[i], "_se")]]
      out$metric <- factor(labels[i], levels = labels)
      out$lower <- out$value - stats::qt(0.975, out$n_reps - 1L) * out$se
      out$upper <- out$value + stats::qt(0.975, out$n_reps - 1L) * out$se
      out
    }))
  }
  theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "#f0f2f4"),
                   plot.margin = ggplot2::margin(10, 15, 10, 10))
  model_colors <- c("Spatial COS" = "#0072B2", "Nonspatial COS" = "#D55E00")
  x_label <- "Nominal exponential range (m); None = no spatial effect"
  figures <- file.path(out_dir, "figures")
  dir.create(figures, showWarnings = FALSE)
  save <- function(plot, name, height = 6) {
    ggplot2::ggsave(file.path(figures, name), plot, width = 10, height = height,
                    dpi = 150, bg = "white")
  }
  accuracy <- longer(primary, c("RMSPE", "CRPS"), c("RMSPE", "CRPS"))
  p <- ggplot2::ggplot(accuracy, ggplot2::aes(case, value, color = model, group = model)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15,
                           position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::facet_grid(metric ~ target, scales = "free_y") +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::labs(x = x_label, y = "Mean score (lower is better)", color = NULL,
                   caption = "Bars: 95% Monte Carlo intervals for the mean across independent simulated datasets.") + theme
  save(p, "spatial-guidance-accuracy.png")
  calibration <- longer(primary, c("coverage_95", "mean_interval_width"),
                          c("95% interval coverage", "Mean interval width"))
  p <- ggplot2::ggplot(calibration, ggplot2::aes(case, value, color = model, group = model)) +
    ggplot2::geom_hline(data = data.frame(metric = "95% interval coverage", value = 0.95),
                         ggplot2::aes(yintercept = value), inherit.aes = FALSE,
                         linetype = 2, color = "grey40") +
    ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::facet_grid(metric ~ target, scales = "free_y") +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::labs(x = x_label, y = NULL, color = NULL,
                   caption = "Coverage is evaluated against known latent truth or held-out responses, as labeled.") + theme
  save(p, "spatial-guidance-calibration.png")
  gains <- longer(improvement, c("RMSPE_improvement_pct", "CRPS_improvement_pct"),
                    c("RMSPE improvement", "CRPS improvement"))
  p <- ggplot2::ggplot(gains[gains$prior_scale == 500, ], ggplot2::aes(case, value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15) +
    ggplot2::geom_point(size = 2.4, color = "#0072B2") +
    ggplot2::facet_wrap(~ metric + target, ncol = 2, scales = "free_y") +
    ggplot2::labs(x = x_label, y = "Mean paired improvement (%)",
                   caption = "Positive favors spatial COS. Bars: 95% Monte Carlo intervals\nfor mean paired improvement across simulated datasets.") + theme
  save(p, "spatial-guidance-improvement.png")
  gains$prior_scale <- factor(gains$prior_scale, c(250, 500, 1000))
  p <- ggplot2::ggplot(gains, ggplot2::aes(case, value, color = prior_scale)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.12,
                           position = ggplot2::position_dodge(width = 0.45), alpha = 0.65) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.45)) +
    ggplot2::facet_wrap(~ metric + target, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("250" = "#009E73", "500" = "#0072B2", "1000" = "#D55E00")) +
    ggplot2::labs(x = x_label, y = "Mean paired improvement (%)", color = "Variance-prior scale",
                   caption = "Same simulated datasets and folds at every prior scale.\nPositive favors spatial COS.") + theme
  save(p, "spatial-guidance-prior-sensitivity.png")
  design$east <- design$X - min(design$X)
  design$north <- design$Y - min(design$Y)
  p <- ggplot2::ggplot(design, ggplot2::aes(east, north, color = factor(fold))) +
    ggplot2::geom_point(size = 3) + ggplot2::coord_equal() +
    ggplot2::labs(x = "Relative easting (m)", y = "Relative northing (m)", color = "CV fold",
                   caption = "Points show plot centroids. All scenarios use this layout and fold assignment.") + theme
  save(p, "spatial-guidance-design.png", height = 6.5)
  invisible(figures)
}
if (!identical(Sys.getenv("CTV_DEFINE_ONLY"), "TRUE")) plot_spatial_guidance_study()
