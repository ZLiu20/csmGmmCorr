\name{csm_extrapolation}
\alias{csm_extrapolation}
\title{Extrapolation diagnostic}
\description{
Recomputes the local false discovery rate of each discovery after winsorizing
\eqn{|z_k|} at the largest fitted magnitude of coordinate \eqn{k} plus
\code{margin}, and flags discoveries whose recomputed score exceeds the largest
unadjusted score among the discoveries. This is a descriptive check, not an
inferential procedure.
}
\usage{
csm_extrapolation(fit, Z, target = c("pleiotropy", "replication", "mediation"),
                  t = 1L, rejected, cutoff, pairs = NULL, margin = 3)
}
\arguments{
\item{fit}{a \code{\link{csm_fit}} object.}
\item{Z}{test statistics.}
\item{target, t, pairs}{see \code{\link{csm_lfdr}}.}
\item{rejected}{logical vector of discoveries.}
\item{cutoff}{largest unadjusted score among the discoveries.}
\item{margin}{winsorization margin.}
}
\value{A list with \code{flag} (\code{NA} for non-discoveries),
\code{lfdr_winsorized} and \code{caps}.}
