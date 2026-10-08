\name{csm_params}
\alias{csm_params}
\alias{as_sun_lists}
\title{Parameters of a correlated csmGmm}
\description{
\code{csm_params} builds and checks a parameter for a given specification.
\code{as_sun_lists} converts a parameter to the list format of the csmGmm
package.
}
\usage{
csm_params(spec, varpi = NULL, mu, Sigma = diag(spec$K), pi = NULL)

as_sun_lists(object, spec = NULL)
}
\arguments{
\item{spec}{a \code{\link{csm_spec}} object.}
\item{varpi}{total probability of each support and magnitude component, either
a vector with one entry per group (support and component, in the order of
\code{spec$grp_b} and \code{spec$grp_m}) or a list of length \eqn{2^K} whose
element \eqn{b+1} is the vector for support \eqn{b}.}
\item{pi}{alternatively, the probability of one signed configuration and
component, \eqn{\pi_{bm} = \varpi_{bm}/2^{|A_b|}}, in either format.}
\item{mu}{magnitudes: a group-by-\eqn{K} matrix with zeros outside each
support; a list of length \eqn{2^K} whose element \eqn{b+1} is a
\eqn{K \times M_b} matrix (the csmGmm \code{initMuList} format); or a vector of
length \eqn{K}, used for every support.}
\item{Sigma}{correlation matrix.}
\item{object}{a \code{csm_fit} object or a parameter.}
}
\value{
\code{csm_params} returns an object of class \code{csm_params}, a list with
\code{varpi}, \code{mu} and \code{Sigma}. It stops if the parameter violates the
simplex, the magnitude box, the ordering constraint, the restrictions of
\code{spec}, or \eqn{\Sigma \succeq \kappa I}. \code{as_sun_lists} returns a
list with \code{muList}, \code{piList} and \code{Sigma}.
}
\examples{
spec <- csm_spec(3, b_mu = 10)
lam <- c(0.8, 0.12, 0.05, 0.03)   # probabilities of support sizes 0..3
p <- csm_params(spec, varpi = lam[spec$grp_size + 1] / choose(3, spec$grp_size),
                mu = c(2, 2.5, 3), Sigma = 0.8 * diag(3) + 0.2)
str(as_sun_lists(p, spec)$piList)
}
