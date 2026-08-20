\name{cos_summary_methods}
\alias{cos_summary_methods}
\alias{summary.cos_prep}
\alias{print.summary.cos_prep}
\alias{summary.cos_recovery_B}
\alias{print.summary.cos_recovery_B}
\alias{summary.cos_prediction_fine}
\alias{print.summary.cos_prediction_fine}
\alias{summary.cos_areal_blocks}
\alias{print.summary.cos_areal_blocks}
\alias{summary.cos_prediction_areal}
\alias{print.summary.cos_prediction_areal}
\alias{summary.cos_cv_observed}
\alias{print.summary.cos_cv_observed}
\title{Summarize COS Package Objects}
\usage{
\method{summary}{cos_prep}(object, digits = 3, ...)

\method{print}{summary.cos_prep}(x, ...)

\method{summary}{cos_recovery_B}(object, probs = c(0.5, 0.025, 0.975),
  digits = 3, include_eta_B = TRUE, ...)

\method{print}{summary.cos_recovery_B}(x, ...)

\method{summary}{cos_prediction_fine}(object, digits = 3, ...)

\method{print}{summary.cos_prediction_fine}(x, ...)

\method{summary}{cos_areal_blocks}(object, digits = 3, ...)

\method{print}{summary.cos_areal_blocks}(x, ...)

\method{summary}{cos_prediction_areal}(object, digits = 3, ...)

\method{print}{summary.cos_prediction_areal}(x, ...)

\method{summary}{cos_cv_observed}(object, digits = 3, ...)

\method{print}{summary.cos_cv_observed}(x, ...)
}
\arguments{
\item{object}{A \code{"cos_prep"} object from \code{\link{cos_prepare}}, a
\code{"cos_recovery_B"} object from \code{\link{cos_recover_B}}, a
\code{"cos_prediction_fine"} object from \code{\link{cos_predict_fine}}, a
\code{"cos_areal_blocks"} object from \code{\link{cos_make_areal_blocks}}, a
\code{"cos_prediction_areal"} object from \code{\link{cos_predict_areal}}, or
a \code{"cos_cv_observed"} object from \code{\link{cos_cv_observed}},
matching the method being dispatched.}
\item{x}{The matching \code{"summary.*"} object returned by one of the
\code{summary()} methods documented here, for use by the corresponding
\code{print()} method.}
\item{digits}{Number of digits used when rounding numeric summaries.}
\item{probs}{Quantile probabilities used for the beta and observed-support
effect summaries in \code{summary.cos_recovery_B()}.}
\item{include_eta_B}{Logical; if \code{TRUE}, \code{summary.cos_recovery_B()}
also summarizes recovered observed-support latent means \code{eta_B}.}
\item{...}{Additional arguments passed to methods. Currently unused.}
}
\description{
S3 \code{summary()} methods for the main objects returned by
\pkg{cosTapered}, together with matching \code{print()} methods for the
returned summary objects. These methods report fuller diagnostics and
posterior or prediction summaries than the compact overviews printed by
\code{\link{cos_print_methods}}.
}
\value{
Each \code{summary()} method documented here returns a list with class
\code{"summary.<x>"} matching the input object's class, and each matching
\code{print()} method prints a compact overview and invisibly returns its
input.
\itemize{
  \item \code{summary.cos_prep()} returns \code{n_B}, \code{n_A}, \code{p},
    \code{x_names}, \code{spatial}, \code{gamma}, \code{taper_code},
    \code{missing}, and \code{row_sum_summary}.
  \item \code{summary.cos_recovery_B()} returns \code{n_samples},
    \code{n_B}, \code{spatial}, \code{beta_summary}, and
    \code{omega_B_summary}. When \code{include_eta_B = TRUE}, it also
    returns \code{eta_B_summary}.
  \item \code{summary.cos_prediction_fine()} returns
    \code{spatial_uncertainty}, \code{uncertainty}, \code{target},
    \code{n_pred}, \code{n_samples}, \code{keep_samples}, \code{spatial},
    rounded five-number summaries for \code{eta} and \code{omega} means and
    standard deviations, and the full \code{prediction_summary} data frame.
    When \code{target = "observed"}, it also returns rounded five-number
    summaries for \code{y}.
  \item \code{summary.cos_areal_blocks()} returns \code{n_U},
    \code{x_names}, \code{missing}, \code{row_sum_summary}, and
    \code{n_cell_summary}.
  \item \code{summary.cos_prediction_areal()} returns the same fields as
    \code{summary.cos_prediction_fine()}, describing areal prediction units
    in place of fine prediction cells, plus \code{row_sum_summary} and
    \code{n_U} in place of \code{n_pred}.
  \item \code{summary.cos_cv_observed()} returns \code{k}, \code{target},
    \code{spatial}, \code{n}, and rounded \code{summary}, \code{fold_summary},
    and \code{predictions} data frames.
}
}
\examples{
\dontrun{
summary(prep)
summary(rec_B, include_eta_B = TRUE)
summary(pred)
summary(U_blocks)
summary(pred_U)
summary(cv)
}
}
