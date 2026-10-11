\name{csmGmmCorr-package}
\alias{csmGmmCorr-package}
\alias{csmGmmCorr}
\docType{package}
\title{Correlated csmGmm for composite null testing}
\description{
The correlated conditionally symmetric multidimensional Gaussian mixture extends
csmGmm (Sun, McCaw and Lin, 2025) to correlated test statistics by a common
within-feature correlation matrix. Support-specific weights, support-specific
magnitude components and the magnitude-ordering constraint of csmGmm are kept,
and all \eqn{3^K} signed configurations are represented. One fit serves the
pleiotropy, replication and mediation targets \eqn{K(h) > t} for every
threshold \eqn{t}.

Typical use: \code{\link{csm_spec}} to specify the model,
\code{\link{csm_fit}} to estimate it on training statistics, and
\code{\link{csm_test}} to compute local false discovery rates and discoveries
on test statistics. \code{\link{csm_adjust}} gives the optional monotone
adjustment.
}
\references{
Sun, R., McCaw, Z. R., and Lin, X. (2025). Testing a large number of composite
null hypotheses using conditionally symmetric multidimensional Gaussian
mixtures in genome-wide studies. \emph{Journal of the American Statistical
Association} 120, 605--617.

Sun, W., and Cai, T. T. (2007). Oracle and adaptive compound decision rules for
false discovery rate control. \emph{Journal of the American Statistical
Association} 102, 901--912.

Meng, X.-L., and Rubin, D. B. (1993). Maximum likelihood estimation via the
ECM algorithm: a general framework. \emph{Biometrika} 80, 267--278.
}
\keyword{package}
