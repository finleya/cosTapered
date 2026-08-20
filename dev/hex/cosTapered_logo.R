rm(list = ls())

library(ggplot2)
library(scico)

set.seed(1)

# ---- Logo controls ---------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if(length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}

panel_file <- file.path(script_dir, "cosTapered_gp_pixels_panel.png")
sticker_file <- normalizePath(
  file.path(script_dir, "..", "..", "cosTapered", "man", "figures", "logo.png"),
  mustWork = FALSE
)
vignette_sticker_file <- normalizePath(
  file.path(script_dir, "..", "..", "cosTapered", "vignettes", "figures", "logo.png"),
  mustWork = FALSE
)

img_width <- 2600
img_height <- 1900
gp_grid_n <- 132
gp_range <- 0.13
gp_nugget <- 1e-8

plot_grid_n_col <- 4
plot_grid_n_row <- 3
plot_grid_spacing <- 0.22
plot_grid_center_x <- 0.4
plot_grid_center_y <- 0.34
n_plots <- plot_grid_n_col * plot_grid_n_row
plot_radius <- rep(0.046, n_plots)
plot_response_noise_sd <- 0.00

# Edit these controls directly and rerun this script.
# x and y shifts are in surface-coordinate units; z shift is in rendered
# height units for the flat plot plane.
# Useful nudges are small, usually between -0.06 and 0.06.
plot_grid_x_shift <- 0.00
plot_grid_y_shift <- 0.00
plot_grid_z_shift <- 0.10
plot_plane_z_base <- 0.36
plot_grid_x <- plot_grid_center_x +
  (seq_len(plot_grid_n_col) - mean(seq_len(plot_grid_n_col))) *
  plot_grid_spacing
plot_grid_y <- plot_grid_center_y +
  (seq_len(plot_grid_n_row) - mean(seq_len(plot_grid_n_row))) *
  plot_grid_spacing * 4 / 3
plot_grid <- expand.grid(
  x = plot_grid_x,
  y = plot_grid_y,
  KEEP.OUT.ATTRS = FALSE
)
plot_centers <- data.frame(
  id = seq_along(plot_radius),
  x = plot_grid$x,
  y = plot_grid$y,
  radius = plot_radius
)

plot_shadow_color <- grDevices::adjustcolor("#000000", alpha.f = 0.45)
plot_edge_color <- grDevices::adjustcolor("#031017", alpha.f = 0.78)
plot_plane_z <- function() plot_plane_z_base + plot_grid_z_shift
plot_x_range <- range(plot_centers$x + plot_grid_x_shift)
plot_y_range <- range(plot_centers$y + plot_grid_y_shift)
if(any(plot_x_range < 0 | plot_x_range > 1) ||
   any(plot_y_range < 0 | plot_y_range > 1)) {
  warning(
    "The shifted plot grid extends outside the 0 to 1 surface domain; ",
    "reduce plot_grid_x_shift or plot_grid_y_shift.",
    call. = FALSE
  )
}

#surface_palette <- c(
#  "#1F4D5F", "#1F8A86", "#35B7A4", "#BFE8A7",
#  "#FFF5B7", "#F9B24A", "#DD5A2A", "#8E2740"
#)
surface_palette <- scico(30, palette = 'bam')

sticker_hex_fill <- "#050505"
sticker_hex_border <- "#FFFFFF"
sticker_hex_border_size <- 0.25

sticker_image_x <- 1.00
sticker_image_y <- 1.22
sticker_image_width <- 1.36
sticker_image_height <- 1.18

sticker_text_left <- "cos"
sticker_text_right <- "Tapered"
sticker_text_x <- 1.08
sticker_text_y <- 0.405
sticker_text_width <- 1.30
sticker_text_family <- "Latin Modern Roman"
sticker_text_face <- "plain"
sticker_text_size <- 28
sticker_text_shadow_offset <- 0.012
sticker_panel_pad <- 0.008
sticker_padding_px <- 5

font_files <- c(
  "Latin Modern Roman" = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf",
  "Liberation Sans" = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
)

if(requireNamespace("sysfonts", quietly = TRUE) &&
   requireNamespace("showtext", quietly = TRUE) &&
   sticker_text_family %in% names(font_files) &&
   file.exists(font_files[[sticker_text_family]])) {
  sysfonts::font_add(sticker_text_family, regular = font_files[[sticker_text_family]])
  showtext::showtext_auto(TRUE)
}

# ---- Surface and plot-support data ----------------------------------------

cov_1d_sqexp <- function(u, range, sigma_sq = 1) {
  D <- as.matrix(dist(u))
  sigma_sq * exp(-0.5 * (D / range)^2)
}

simulate_gp_surface <- function(n = gp_grid_n,
                                range = gp_range,
                                nugget = gp_nugget) {
  x <- seq(0, 1, length.out = n)
  y <- seq(0, 1, length.out = n)

  Lx <- chol(cov_1d_sqexp(x, range) + diag(nugget, n))
  Ly <- chol(cov_1d_sqexp(y, range * 1.10) + diag(nugget, n))

  z <- t(Lx) %*% matrix(rnorm(n * n), n, n) %*% Ly

  # Add a weak large-scale trend so the field reads as a spatial covariate
  # surface rather than pure texture.
  xx <- matrix(rep(x, n), n, n)
  yy <- matrix(rep(y, each = n), n, n)
  z <- z + 0.35 * sin(2 * pi * xx) - 0.25 * cos(2 * pi * yy)

  z <- z - min(z)
  z <- z / max(z)
  z <- z^1.05

  list(x = x, y = y, z = z)
}

surface <- simulate_gp_surface()

surface_dat <- expand.grid(
  x = surface$x,
  y = surface$y,
  KEEP.OUT.ATTRS = FALSE
)
surface_dat$z <- as.vector(surface$z)

plot_centers$response <- vapply(seq_len(nrow(plot_centers)), function(i) {
  plot_x <- plot_centers$x[i] + plot_grid_x_shift
  plot_y <- plot_centers$y[i] + plot_grid_y_shift
  d <- sqrt(
    (surface_dat$x - plot_x)^2 +
      (surface_dat$y - plot_y)^2
  )
  mean(surface_dat$z[d <= plot_centers$radius[i]])
}, numeric(1))
if(plot_response_noise_sd > 0) {
  plot_centers$response <- plot_centers$response +
    rnorm(nrow(plot_centers), 0, plot_response_noise_sd)
}
plot_centers$response <- pmax(0, pmin(1, plot_centers$response))

plot_response_palette <- grDevices::colorRampPalette(surface_palette)(256)
plot_response_idx <- pmax(1L, pmin(256L, round(plot_centers$response * 255) + 1L))
sticker_text_left_color <- plot_response_palette[[min(plot_response_idx)]]
sticker_text_right_color <- plot_response_palette[[max(plot_response_idx)]]

# ---- Render 3D panel -------------------------------------------------------

face_colors <- function(z, pal = surface_palette, n = 256) {
  ramp <- grDevices::colorRampPalette(pal)(n)
  z_face <- (z[-nrow(z), -ncol(z)] +
               z[-1L, -ncol(z)] +
               z[-nrow(z), -1L] +
               z[-1L, -1L]) / 4
  idx <- pmax(1L, pmin(n, round(z_face * (n - 1L)) + 1L))
  matrix(ramp[idx], nrow = nrow(z_face), ncol = ncol(z_face))
}

draw_circle3d <- function(pm, center_x, center_y, radius, z_base,
                          fill, edge = "#031017") {
  theta <- seq(0, 2 * pi, length.out = 181L)
  x <- center_x + radius * cos(theta)
  y <- center_y + radius * sin(theta)
  z <- rep(z_base, length(theta))
  xy <- grDevices::trans3d(x, y, z, pm)

  shadow_xy <- grDevices::trans3d(
    center_x + 0.030 + radius * cos(theta),
    center_y - 0.030 + radius * sin(theta),
    z_base - 0.035,
    pm
  )
  polygon(shadow_xy$x, shadow_xy$y, col = plot_shadow_color, border = NA)
  polygon(xy$x, xy$y, col = fill, border = edge, lwd = 1.0)
}

draw_3d_panel <- function(surface, file) {
  if(requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      file,
      width = img_width,
      height = img_height,
      units = "px",
      background = "transparent",
      res = 300
    )
  } else {
    png(file, width = img_width, height = img_height,
        bg = "transparent", res = 300)
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = rep(0, 4), bg = NA)

  x <- surface$x
  y <- surface$y
  z <- surface$z
  z3 <- 0.46 * z - 0.06

  cols <- face_colors(z)

  pm <- persp(
    x = x,
    y = y,
    z = z3,
    theta = -38,
    phi = 40,
    expand = 0.70,
    scale = FALSE,
    col = cols,
    border = grDevices::adjustcolor("#171717", alpha.f = 0.20),
    ltheta = -35,
    lphi = 55,
    shade = 0.18,
    axes = FALSE,
    box = FALSE,
    xlab = "",
    ylab = "",
    zlab = "",
    xlim = c(-0.07, 1.08),
    ylim = c(-0.06, 1.08),
    zlim = c(-0.08, 0.68)
  )

  # Draw back-to-front so overlapping disks and shadows obey the camera view.
  plot_x <- plot_centers$x + plot_grid_x_shift
  plot_y <- plot_centers$y + plot_grid_y_shift
  draw_order <- order(plot_y, decreasing = TRUE)
  for(ii in draw_order) {
    draw_circle3d(
      pm = pm,
      center_x = plot_x[ii],
      center_y = plot_y[ii],
      radius = plot_centers$radius[ii],
      z_base = plot_plane_z(),
      fill = plot_response_palette[plot_response_idx[ii]],
      edge = plot_edge_color
    )
  }

  invisible(TRUE)
}

draw_3d_panel(surface, panel_file)

message("Wrote pixel surface panel PNG: ", normalizePath(panel_file, mustWork = FALSE))

# ---- Hex sticker -----------------------------------------------------------

save_sticker_plot <- function(sticker_plot, output_file, panel_pad) {
  center <- 1
  radius <- 1
  half_width <- sqrt(3) / 2 * radius
  border_room <- sticker_hex_border_size * 0.04
  x_pad <- half_width * panel_pad + border_room
  y_pad <- radius * panel_pad + border_room

  sticker_plot <- suppressMessages(
    sticker_plot +
      coord_fixed(clip = "off") +
      scale_x_continuous(
        expand = c(0, 0),
        limits = c(center - half_width - x_pad, center + half_width + x_pad)
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(center - radius - y_pad, center + radius + y_pad)
      ) +
      theme(plot.margin = margin(0, 0, 0, 0, unit = "lines"))
  )

  ggsave(
    filename = output_file,
    plot = sticker_plot,
    width = 43.9,
    height = 50.8,
    units = "mm",
    bg = "transparent",
    dpi = 600
  )
}

add_sticker_padding <- function(input_file, output_file, padding_px) {
  if(!requireNamespace("magick", quietly = TRUE) || padding_px <= 0) {
    if(!file.copy(input_file, output_file, overwrite = TRUE))
      stop("Could not copy sticker PNG: ", output_file)
    return(invisible(output_file))
  }

  img <- magick::image_read(input_file)
  info <- magick::image_info(img)
  geometry <- sprintf("%dx%d", info$width + 2 * padding_px, info$height + 2 * padding_px)
  img <- magick::image_extent(img, geometry = geometry, gravity = "center", color = "none")
  magick::image_write(img, path = output_file, format = "png")

  invisible(output_file)
}

add_sticker_text <- function(sticker_plot) {
  if(requireNamespace("systemfonts", quietly = TRUE)) {
    widths <- systemfonts::string_width(
      c(sticker_text_left, sticker_text_right),
      family = sticker_text_family,
      size = sticker_text_size
    )
    split <- widths[1] / sum(widths)
  } else {
    split <- nchar(sticker_text_left) /
      nchar(paste0(sticker_text_left, sticker_text_right))
  }

  join_x <- sticker_text_x + (split - 0.5) * sticker_text_width

  sticker_plot +
    annotate(
      "text",
      x = join_x - sticker_text_shadow_offset,
      y = sticker_text_y - sticker_text_shadow_offset,
      label = sticker_text_left,
      hjust = 1,
      color = "#000000",
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    ) +
    annotate(
      "text",
      x = join_x + sticker_text_shadow_offset,
      y = sticker_text_y - sticker_text_shadow_offset,
      label = sticker_text_right,
      hjust = 0,
      color = "#000000",
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    ) +
    annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_left,
      hjust = 1,
      color = sticker_text_left_color,
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    ) +
    annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_right,
      hjust = 0,
      color = sticker_text_right_color,
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    )
}

if(!requireNamespace("hexSticker", quietly = TRUE))
  stop("Package hexSticker is required to build the package logo.")

dir.create(dirname(sticker_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(vignette_sticker_file), recursive = TRUE, showWarnings = FALSE)

sticker_tmp <- tempfile(fileext = ".png")
sticker_panel_tmp <- tempfile(fileext = ".png")

sticker_plot <- hexSticker::sticker(
  subplot = panel_file,
  s_x = sticker_image_x,
  s_y = sticker_image_y,
  s_width = sticker_image_width,
  s_height = sticker_image_height,
  package = "",
  p_x = sticker_text_x,
  p_y = sticker_text_y,
  p_color = "transparent",
  p_family = sticker_text_family,
  p_fontface = sticker_text_face,
  p_size = 0,
  h_fill = sticker_hex_fill,
  h_color = sticker_hex_border,
  h_size = sticker_hex_border_size,
  white_around_sticker = FALSE,
  filename = sticker_tmp,
  dpi = 600
)

sticker_plot <- add_sticker_text(sticker_plot)
save_sticker_plot(sticker_plot, sticker_panel_tmp, sticker_panel_pad)
add_sticker_padding(sticker_panel_tmp, sticker_file, sticker_padding_px)
if(!file.copy(sticker_file, vignette_sticker_file, overwrite = TRUE))
  stop("Could not copy the package logo PNG to the vignette figures directory.")

message("Wrote package logo PNG: ", normalizePath(sticker_file, mustWork = FALSE))
message("Wrote vignette logo PNG: ", normalizePath(vignette_sticker_file, mustWork = FALSE))
