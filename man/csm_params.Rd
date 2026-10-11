\name{csm_params}
\alias{csm_params}
\alias{as_sun_lists}
\title{Construct or convert csmGmm parameters}
\usage{
csm_params(spec, varpi = NULL, mu, Sigma = diag(spec$K), pi = NULL)
as_sun_lists(object, spec = NULL)
}
\arguments{
\item{spec}{A model specification from \code{csm_spec}.}
\item{varpi}{Support-component probabilities, as a vector or support list.}
\item{mu}{Magnitude matrix, coordinate vector, or csmGmm magnitude list.}
\item{Sigma}{Common conditional correlation matrix.}
\item{pi}{Alternative signed-component probabilities in csmGmm list format.}
\item{object}{A parameter object or fitted model.}
}
\value{\code{csm_params()} returns validated parameters;
\code{as_sun_lists()} returns \code{muList}, \code{piList}, and \code{Sigma}.}
\description{Create validated model parameters or convert them to the list
representation used by the original csmGmm implementation.}
\details{
\code{varpi[g]} is the total probability of support-magnitude group \code{g}.
For a support of size \eqn{s}, that mass is divided uniformly among the
\eqn{2^s} signs during likelihood evaluation.  The elements of a csmGmm
\code{piList} use the same block-mass convention; they must not be divided by
\eqn{2^s} before conversion.

Inactive entries of \code{mu} must be zero.  Active entries must lie between
the bounds in \code{spec}, obey the magnitude ordering, and be paired with a
correlation matrix whose smallest eigenvalue is at least \code{kappa}.
}
\examples{
spec <- csm_spec(2, M = 1, b_mu = 10)
weights <- c(0.7, 0.1, 0.1, 0.1)
params <- csm_params(spec, varpi = weights,
                     mu = 2 * spec$grp_A, Sigma = diag(2))
lists <- as_sun_lists(params, spec)
sum(unlist(lists$piList))
}
