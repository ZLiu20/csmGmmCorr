\name{csm_test}
\alias{csm_test}
\alias{print.csm_test}
\title{Discoveries for one target}
\description{
Computes local false discovery rates on the test statistics, applies the
step-up rule, flags extrapolation-driven discoveries, and optionally applies
the monotone adjustment. A failed fit makes no discoveries.
}
\usage{
csm_test(fit, Z, target = c("pleiotropy", "replication", "mediation"), t = 1L,
         q = 0.1, pairs = NULL, adjust = FALSE, extrapolation = TRUE,
         allow_unconverged = FALSE)

\method{print}{csm_test}(x, ...)
}
\arguments{
\item{fit}{a \code{\link{csm_fit}} object.}
\item{Z}{test statistics.}
\item{target, t, pairs}{see \code{\link{csm_lfdr}}.}
\item{q}{level.}
\item{adjust}{also report the monotone adjustment (\code{\link{csm_adjust}}),
computed on all rows of \code{Z}.}
\item{extrapolation}{apply the extrapolation diagnostic
(\code{\link{csm_extrapolation}}).}
\item{allow_unconverged}{use an unconverged fit instead of recording a
failure.}
\item{x}{a \code{csm_test} object.}
\item{...}{ignored.}
}
\value{
An object of class \code{csm_test} with \code{status}, \code{lfdr},
\code{rejected}, \code{n_rejected}, \code{cutoff}, and, if requested,
\code{extrapolation}, \code{lfdr_winsorized} and \code{adjusted}.
}
