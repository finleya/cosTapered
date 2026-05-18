\name{cos_fit-methods}
\alias{cos_fit-methods}
\alias{summary.cos_fit}
\alias{print.summary.cos_fit}
\alias{plot.cos_fit}
\title{Summarize and Plot Fitted COS Models}
\usage{
\method{summary}{cos_fit}(object, burn_in = 1,
  pars = c("tau_B_sq", "tau_sq", "sigma_sq", "phi", "eff_range", "lp"),
  probs = c(0.5, 0.025, 0.975), digits = 3, ...)

\method{print}{summary.cos_fit}(x, ...)

\method{plot}{cos_fit}(x, burn_in = 1,
  pars = c("tau_B_sq", "sigma_sq", "phi", "eff_range"), ...)
}
\arguments{
\item{object}{A \code{"cos_fit"} object from \code{\link{cos_fit}}.}
\item{burn_in}{First retained iteration index.}
\item{pars}{Parameter names to summarize or plot.}
\item{probs}{Quantile probabilities.}
\item{digits}{Number of digits for rounded summaries.}
\item{x}{A \code{"summary.cos_fit"} object for \code{print()}, or a
\code{"cos_fit"} object for \code{plot()}.}
\item{...}{Additional arguments passed to methods. Currently unused.}
}
\description{
S3 methods for inspecting fitted COS model objects. These methods do not fit
the model; they summarize retained posterior samples or draw MCMC trace plots
from an existing \code{"cos_fit"} object.
}
\details{
For non-spatial fits, spatial covariance parameters are stored as
\code{NA_real_}. The summary and plot methods omit requested parameters that
have no finite retained samples.
}
\value{
\code{summary.cos_fit()} returns a list with class
\code{"summary.cos_fit"}. \code{print.summary.cos_fit()} invisibly returns
its input. \code{plot.cos_fit()} returns a \pkg{ggplot2} object.
}
\examples{
\dontrun{
summary(fit, burn_in = 100)
plot(fit, burn_in = 100)
}
}
