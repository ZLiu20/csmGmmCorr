\name{csm_test}
\alias{csm_test}
\alias{print.csm_test}
\title{Run a complete composite-null analysis}
\usage{
csm_test(fit, z, target = c("pleiotropy", "replication", "mediation"),
  t = 1L, q = 0.1, pairs = NULL, adjust = FALSE,
  extrapolation = TRUE, allow_unconverged = FALSE)
\method{print}{csm_test}(x, ...)
}
\arguments{
\item{fit}{A fitted model.}
\item{z}{Matrix of z-statistics.}
\item{target, t, pairs}{Composite-null definition.}
\item{q}{Input false discovery rate level.}
\item{adjust}{Apply the optional monotone adjustment.}
\item{extrapolation}{Compute the winsorization diagnostic.}
\item{allow_unconverged}{Permit explicitly labeled inference from an unconverged fit.}
\item{x}{A test result.}
\item{...}{Additional print arguments.}
}
\value{A \code{csm_test} result.}
\description{Compute target-specific local false discovery rates, apply the
step-up decision rule, and optionally add diagnostics.}
\details{
The raw posterior score is the primary model-based result.  With
\code{adjust = TRUE}, the function replaces it by the majorant returned by
\code{csm_adjust()} before applying the step-up rule.  The majorant enforces
the specified finite-sample ordering on the observed evaluation set; it does
not refit the mixture and is not a general frequentist FDR guarantee under
model misspecification.

The input \code{q} is the target average local false discovery rate among the
selected prefix.  It defaults to 0.1 but remains a user-controlled analysis
level.  An unconverged fit returns no discoveries unless
\code{allow_unconverged = TRUE}; such inference is labeled in the result.
}
\examples{
\dontrun{
result <- csm_test(fit, Z, target = "replication", t = 1,
                   q = 0.1, adjust = FALSE)
result$n_rejected
}
}
