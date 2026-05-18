\name{cosTapered-example-data}
\alias{cosTapered-example-data}
\title{Small Example Dataset for cosTapered}
\description{
A small set of external files included with \pkg{cosTapered} for examples,
manual development checks, and simulation tests.
}
\details{
The files are installed under \file{extdata} and can be located with
\code{system.file("extdata", ..., package = "cosTapered")}. They include:
\itemize{
  \item \file{example_chm.tif}: canopy height raster used as a fine-support
    covariate.
  \item \file{example_B.gpkg}: observed-support polygons with simulated
    response \code{y_B} and simulated latent truth columns.
  \item \file{example_U.gpkg}: areal prediction polygons.
  \item \file{example_forest.gpkg}: forest boundary for plotting.
  \item \file{example_truth.rds}: known simulation parameter values.
}
The example data are a small subset of a larger dataset collected by
Arthur Spitta. The subset was prepared for package examples and tests.

The observed responses in \file{example_B.gpkg} were generated with
\code{\link{cos_simulate_data}} using the included canopy-height raster,
forest boundary, and known parameter values stored in
\file{example_truth.rds}. The development script
\file{dev/make_example_data.R} contains the full reproducible data-generation
workflow: it prepares the source raster and vector files, calls
\code{\link{cos_simulate_data}}, and writes the resulting files into
\file{cosTapered/inst/extdata}.
}
\references{
Zhang, L., Finley, A. O., Nothdurft, A., and Banerjee, S. (2024).
Bayesian modeling of incompatible spatial data: A case study involving
Post-Adrian storm forest damage assessment. \emph{International Journal
of Applied Earth Observation and Geoinformation}, 135, 104224.
\doi{10.1016/j.jag.2024.104224}
}
\examples{
chm_file <- system.file("extdata", "example_chm.tif", package = "cosTapered")
B_file <- system.file("extdata", "example_B.gpkg", package = "cosTapered")
truth_file <- system.file("extdata", "example_truth.rds", package = "cosTapered")

file.exists(chm_file)
file.exists(B_file)
readRDS(truth_file)

\dontrun{
library(sf)
library(raster)
library(ggplot2)

chm <- raster(chm_file)
names(chm) <- "chm"

B_sf <- st_read(B_file, quiet = TRUE)
U_sf <- st_read(system.file("extdata", "example_U.gpkg",
                            package = "cosTapered"), quiet = TRUE)
forest_sf <- st_read(system.file("extdata", "example_forest.gpkg",
                                 package = "cosTapered"), quiet = TRUE)

chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)
names(chm_df) <- c("x", "y", "chm")
tile_width <- xres(chm)
tile_height <- yres(chm)

ggplot() +
  geom_tile(data = chm_df, aes(x = x, y = y, fill = chm),
            width = tile_width, height = tile_height) +
  geom_sf(data = U_sf, fill = NA, color = "white", linewidth = 0.25) +
  geom_sf(data = forest_sf, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = B_sf, fill = NA, color = "#c1121f", linewidth = 0.45) +
  scale_fill_viridis_c(option = "cividis") +
  coord_sf(expand = FALSE) +
  labs(title = "Example supports", fill = "CHM", x = NULL, y = NULL) +
  theme_bw()

ggplot() +
  geom_sf(data = forest_sf, fill = "grey95", color = "grey55",
          linewidth = 0.4) +
  geom_sf(data = B_sf, aes(fill = y_B), color = "black", linewidth = 0.25) +
  scale_fill_viridis_c(option = "magma") +
  coord_sf(expand = FALSE) +
  labs(title = "Simulated observed response", fill = "y_B",
       x = NULL, y = NULL) +
  theme_bw()
}
}
\keyword{datasets}
