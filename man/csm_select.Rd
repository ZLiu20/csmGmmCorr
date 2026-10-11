\name{csm_select}
\alias{csm_select}
\title{Select a magnitude-component structure}
\usage{
csm_select(x, K, structures = .default_structures(), bic_tol = 0,
  starts = c("null", "shared"), control = csm_control(), S = K,
  a_mu = 0.5, b_mu = NULL, kappa = 0.05, seed = NULL)
}
\arguments{
\item{x}{Matrix of z-statistics.}
\item{K}{Number of coordinates.}
\item{structures}{Named component-count candidates.}
\item{bic_tol}{BIC tolerance for choosing the richer tied candidate.}
\item{starts, control}{Fitting starts and control options.}
\item{S, a_mu, b_mu, kappa}{Specification arguments.}
\item{seed}{Optional random seed.}
}
\value{The selected fit with a BIC table and all candidate fits attached.}
\description{Fits candidate component structures and chooses by BIC.}
