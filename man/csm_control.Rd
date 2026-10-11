\name{csm_control}
\alias{csm_control}
\title{Control correlated csmGmm fitting}
\usage{
csm_control(max_iter = 5000L, tol_loglik = NULL, tol_param = Inf,
  accelerate = TRUE, warm_n = 0L, penalty = 30, weight_penalty = 10,
  block_cells = 5e6, sigma_steps = 50L, fix_sigma = FALSE,
  cores = 1L, verbose = FALSE)
}
\arguments{
\item{max_iter}{Maximum number of ECM maps.}
\item{tol_loglik, tol_param}{Likelihood and parameter stopping tolerances.}
\item{accelerate}{Use SQUAREM acceleration.}
\item{warm_n}{Deterministic subset size for an optional warm start.}
\item{penalty, weight_penalty}{Magnitude and within-size weight penalties.}
\item{block_cells}{Maximum working-matrix size in the R reference E-step.}
\item{sigma_steps}{Maximum projected-gradient correlation updates per map.}
\item{fix_sigma}{\code{FALSE}, \code{TRUE} for identity, or a fixed matrix.}
\item{cores}{Requested native-kernel worker count.}
\item{verbose}{Report fitting progress.}
}
\value{A fitting-control list.}
\description{Sets numerical, penalty, acceleration, and covariance options.}
