\name{cos_print_methods}
\alias{cos_print_methods}
\alias{print.cos_prep}
\alias{print.cos_priors}
\alias{print.cos_fit}
\alias{print.cos_recovery_B}
\alias{print.cos_prediction_fine}
\alias{print.cos_areal_blocks}
\alias{print.cos_prediction_areal}
\alias{print.cos_cv_observed}
\alias{print.cos_simulation}
\title{Print COS Package Objects}
\usage{
\method{print}{cos_prep}(x, ...)

\method{print}{cos_priors}(x, ...)

\method{print}{cos_fit}(x, ...)

\method{print}{cos_recovery_B}(x, ...)

\method{print}{cos_prediction_fine}(x, ...)

\method{print}{cos_areal_blocks}(x, ...)

\method{print}{cos_prediction_areal}(x, ...)

\method{print}{cos_cv_observed}(x, ...)

\method{print}{cos_simulation}(x, ...)
}
\arguments{
\item{x}{A COS package object.}
\item{...}{Additional arguments passed to methods. Currently unused.}
}
\description{
S3 \code{print()} methods for the main objects returned by \pkg{cosTapered}.
These methods provide compact object overviews when users type an object name
at the console. Use \code{summary()} for fuller diagnostics and posterior or
prediction summaries.
}
\value{
Each method invisibly returns its input.
}
\examples{
\dontrun{
prep
priors
fit
summary(fit)
}
}
