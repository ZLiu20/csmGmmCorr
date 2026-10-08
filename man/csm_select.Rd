\name{csm_select}
\alias{csm_select}
\alias{csm_structure}
\title{Choose the number of magnitude components per support by BIC}
\description{
Fits the correlated csmGmm with several prespecified numbers of magnitude
components per support and returns the fit with the smallest BIC.
}
\usage{
csm_select(X, K, structures = list(A = c(single = 1L, multiple = 1L),
           B = c(single = 2L, multiple = 1L), C = c(single = 2L, multiple = 2L)), bic_tol = 0,
           starts = c("null", "shared"), control = csm_control(),
           S = K, a_mu = 0.5, b_mu = NULL, kappa = 0.05, seed = NULL)
csm_structure(K, single = 1L, multiple = 1L)
}
\arguments{
  \item{X}{J-by-K matrix of test statistics.}
  \item{K}{number of columns of \code{X}.}
  \item{structures}{named list; each element gives the number of components for supports of
    size one (\code{single}) and of size two or more (\code{multiple}). The null support has one.}
  \item{bic_tol}{among structures whose BIC is within \code{bic_tol} of the smallest, the one with
    the most parameters is chosen (0: plain BIC).}
  \item{starts, control, seed}{passed to \code{\link{csm_fit}}.}
  \item{S, a_mu, b_mu, kappa}{passed to \code{\link{csm_spec}}.}
  \item{single, multiple}{numbers of components for supports of size one and of size two or more.}
}
\details{
The alternative distribution of a support is a finite Gaussian mixture; when effect sizes vary,
too few components can let one support absorb the spread of another and bias the local false
discovery rates, whereas additional components mainly increase computation. The BIC uses the
log-likelihood at the fitted (by default penalized) parameters and the number of free weights,
magnitudes and correlations (see \code{logLik.csm_fit}).
}
\value{
The chosen \code{csm_fit} object, with elements \code{selection} (one row per structure: numbers
of components, degrees of freedom, log-likelihood, BIC, difference from the smallest BIC, whether it
converged, time, and which was chosen) and \code{all_fits}. \code{csm_structure} returns the vector
\code{M} of length \eqn{2^K} for \code{\link{csm_spec}}.
}
\examples{
\donttest{
sp <- csm_spec(2, b_mu = 10)
p <- csm_params(sp, varpi = c(0.9, 0.04, 0.04, 0.02), mu = sp$grp_A * 3)
z <- csm_simulate(3000, p, sp)$Z
fs <- csm_select(z, 2, structures = list(A = c(single = 1, multiple = 1),
                                         B = c(single = 2, multiple = 1)),
                 control = csm_control(max_iter = 200))
fs$selection
}
}
