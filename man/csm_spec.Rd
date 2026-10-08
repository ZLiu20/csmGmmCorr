\name{csm_spec}
\alias{csm_spec}
\alias{print.csm_spec}
\title{Specify a correlated csmGmm}
\description{
Fixes the model structure and the parameter region before fitting: the number
of statistics, the magnitude components of each support, the largest support
size, the magnitude box and the eigenvalue bound on the correlation matrix.
Every support has its own weights and magnitudes, as in csmGmm.
}
\usage{
csm_spec(K, M = csm_structure(K, 2L, 1L), S = K, a_mu = 0.5, b_mu = NULL, kappa = 0.05)

\method{print}{csm_spec}(x, ...)
}
\arguments{
\item{K}{number of statistics per feature (2 to 12).}
\item{M}{number of magnitude components per support: a single integer, or a
vector of length \eqn{2^K} whose element \eqn{b+1} refers to the support with
binary label \eqn{b} (coordinate 1 is the most significant bit, as in csmGmm).
The null support always has one component. The default,
\code{csm_structure(K, 2, 1)}, gives two components to each support of size one
and one component to each larger support (the structure of the example in Sun,
McCaw and Lin, 2025): supports of size one carry most non-null features and their
effect sizes are the most heterogeneous. \code{M = 1} gives one component per support.}
\item{S}{largest support size represented. \code{S = K} (default) is the
correlated csmGmm; \code{S < K} truncates it to supports of size at most
\code{S}, keeping support-specific weights, magnitudes and the ordering
constraint for the supports retained. The ordering constraint involves the
full support and is therefore vacuous when \code{S < K}.}
\item{a_mu, b_mu}{lower and upper bounds on the magnitudes. If \code{b_mu} is
\code{NULL}, \code{\link{csm_fit}} sets it to \code{ceiling(max(abs(X))) + 1}.}
\item{kappa}{lower bound on the eigenvalues of the correlation matrix; it
should not exceed the smallest eigenvalue of the null-statistic correlation.}
\item{x}{a \code{csm_spec} object.}
\item{...}{ignored.}
}
\value{An object of class \code{csm_spec}.}
\details{
The ordering constraint requires every full-support magnitude component to be at
least as large, coordinatewise, as every magnitude component of every other
support.
}
\examples{
spec <- csm_spec(4)
spec
}
