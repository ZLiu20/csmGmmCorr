\name{csm_fit}
\alias{csm_fit}
\alias{print.csm_fit}
\alias{predict.csm_fit}
\alias{logLik.csm_fit}
\title{Fit a correlated csmGmm}
\usage{
csm_fit(x, spec, starts = c("null", "equal", "shared"), n_random = 0L,
  init = NULL, control = csm_control(), seed = NULL)
\method{print}{csm_fit}(x, ...)
\method{predict}{csm_fit}(object, newdata,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, pairs = NULL, ...)
\method{logLik}{csm_fit}(object, ...)
}
\arguments{
\item{x, newdata}{Numeric matrices with features in rows and statistics in columns.}
\item{spec}{A model specification.}
\item{starts}{Named deterministic starting values.}
\item{n_random}{Number of random additional starts.}
\item{init}{One parameter object or a list of parameter objects.}
\item{control}{Options from \code{csm_control}.}
\item{seed}{Optional random seed.}
\item{object}{A fitted model.}
\item{target, t, pairs}{Composite-null definition.}
\item{...}{Additional method arguments.}
}
\value{A \code{csm_fit} object; prediction returns local false discovery rates.}
\description{Fit the constrained correlated mixture from one or more starting
values and retain the eligible run with the largest objective.}
\details{
Each ECM map updates support-component weights, constrained active
magnitudes, and, unless fixed, the common correlation matrix.  Likelihood
terms are evaluated with log-sum-exp stabilization.  Native code pairs global
sign reversals and may use several worker threads; this is an exact
computational reorganization of the retained state space.

The returned \code{runs} table records the status of every initialization.
\code{converged} is true only when the selected run satisfies the stopping
rule.  If no run converges but a finite run exists, its parameters are kept
for diagnosis and the fit is explicitly marked unconverged.  By default,
\code{csm_test()} makes no discoveries from such a fit.

\code{fix_sigma = TRUE} fixes the covariance to the identity.  Supplying a
matrix fixes it to that correlation matrix.  These are different statistical
models from estimating \code{Sigma}.
}
\examples{
\dontrun{
spec <- csm_spec(ncol(Z), b_mu = 15)
fit <- csm_fit(Z, spec, starts = c("null", "shared"),
               control = csm_control(cores = 4))
fit
predict(fit, Z, target = "pleiotropy", t = 1)
}
}
