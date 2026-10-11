\name{csm_extrapolation}
\alias{csm_extrapolation}
\title{Diagnose magnitude extrapolation}
\usage{
csm_extrapolation(fit, z,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, rejected, cutoff, pairs = NULL, margin = 3)
}
\arguments{
\item{fit}{A fitted model.}
\item{z}{Matrix of z-statistics.}
\item{target, t, pairs}{Composite-null definition.}
\item{rejected}{Logical rejection vector.}
\item{cutoff}{Local-fdr cutoff from the step-up rule.}
\item{margin}{Amount added to the largest fitted magnitude before winsorization.}
}
\value{Flags, winsorized local false discovery rates, and coordinate caps.}
\description{Flags discoveries that disappear after extreme coordinates are
winsorized to a range supported by the fitted magnitude components.}
