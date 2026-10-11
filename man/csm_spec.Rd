\name{csm_spec}
\alias{csm_spec}
\alias{csm_structure}
\title{Specify a correlated csmGmm}
\usage{
csm_spec(K, M = csm_structure(K, 2L, 1L), S = K,
  a_mu = 0.5, b_mu = NULL, kappa = 0.05)
csm_structure(K, single = 2L, multiple = 1L)
}
\arguments{
\item{K}{Number of statistic coordinates.}
\item{M}{Number of magnitude components for every binary support.}
\item{S}{Largest support size retained in the model.}
\item{a_mu, b_mu}{Lower and optional upper bounds for active magnitudes.}
\item{kappa}{Lower bound on eigenvalues of the common correlation matrix.}
\item{single, multiple}{Component counts for size-one and larger supports.}
}
\value{\code{csm_spec()} returns a model specification. \code{csm_structure()}
returns an integer vector of length \eqn{2^K}.}
\description{Build the finite state space and magnitude-component structure
used by a correlated csmGmm fit.}
\details{
Supports are represented by integers from zero through \eqn{2^K-1}.  The
first coordinate is the most significant bit.  A support of size \eqn{s}
contributes \eqn{2^s} signed configurations, but the implementation evaluates
opposite signs as a pair without deleting either state.

\code{M[b + 1]} is the number of magnitude components for support \code{b}.
\code{S} is a modeling restriction: supports larger than \code{S} are absent,
so using \code{S < K} changes the fitted probability model.  In contrast,
\code{block_cells} and \code{cores} in other functions change computation but
not the retained state space.

For the full state space use \code{S = K}.  The default structure gives two
components to size-one supports and one component to every larger support.
}
\examples{
structure <- csm_structure(4, single = 2, multiple = 1)
spec <- csm_spec(4, M = structure, S = 4, b_mu = 15)
spec
table(spec$grp_size)
}
