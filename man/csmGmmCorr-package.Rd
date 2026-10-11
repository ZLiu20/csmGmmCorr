\name{csmGmmCorr-package}
\alias{csmGmmCorr-package}
\alias{csmGmmCorr}
\title{Correlated conditionally symmetric Gaussian mixtures}
\description{
The package fits a conditionally symmetric multidimensional Gaussian mixture
with a common correlation matrix.  Rows are features and columns are the
possibly correlated statistics observed for each feature.  One fitted model
can be scored under pleiotropy, same-sign replication, or paired-mediation
composite null hypotheses.
}
\details{
For every retained binary support, the model assigns one or more magnitude
components and distributes each component's probability uniformly over its
signed configurations.  All components share a positive-definite correlation
matrix.  Active magnitudes are estimated subject to bounds and the specified
cross-support ordering restriction.

The main workflow is \code{csm_spec()}, \code{csm_fit()}, \code{csm_lfdr()},
and \code{csm_stepup()}.  The convenience function \code{csm_test()} combines
the final two steps.  \code{csm_adjust()} is an optional order-enforcing
post-processing operation; it is not another mixture fit.  Convergence status
and diagnostics should be inspected before scientific interpretation.

This software extends the csmGmm framework of Sun, McCaw, and Lin (2025),
doi:10.1080/01621459.2024.2422124.
}
\docType{package}
