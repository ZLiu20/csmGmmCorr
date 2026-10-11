\name{csm_fit}
\alias{csm_fit}
\alias{print.csm_fit}
\alias{logLik.csm_fit}
\alias{predict.csm_fit}
\title{Fit a correlated csmGmm}
\description{
Maximizes the likelihood of the correlated csmGmm over the prespecified region
by a constrained expectation conditional maximization (ECM) algorithm, run from
several starting values. The converged fit with the largest training
likelihood is selected. If no start converges, the fit is recorded as a
failure, and \code{\link{csm_test}} makes no discoveries.
}
\usage{
csm_fit(X, spec, starts = c("null", "equal", "shared"), n_random = 0L,
        init = NULL, control = csm_control(), seed = NULL)

\method{print}{csm_fit}(x, ...)

\method{logLik}{csm_fit}(object, ...)

\method{predict}{csm_fit}(object, newdata,
        target = c("pleiotropy", "replication", "mediation"),
        t = 1L, pairs = NULL, ...)
}
\arguments{
\item{X}{training statistics, an \eqn{N \times K} matrix. They should be
independent of the test statistics.}
\item{spec}{a \code{\link{csm_spec}} object.}
\item{starts}{built-in starting values: \code{"null"} (half of the mass on
the null support, the rest spread equally), \code{"equal"} (equal weights on
all supports and components), and \code{"shared"} (a fit of the
model with one magnitude per coordinate shared by all supports, embedded in
the correlated csmGmm; it is used only as a starting value). All use magnitudes
satisfying the ordering constraint; the first two use \eqn{\Sigma = I}.}
\item{n_random}{number of additional random perturbations of the
\code{"null"} start.}
\item{init}{optional list of \code{\link{csm_params}} objects used as
additional starting values.}
\item{control}{see \code{\link{csm_control}}.}
\item{seed}{optional seed for the random starts.}
\item{x, object}{a \code{csm_fit} object.}
\item{newdata}{test statistics.}
\item{target, t, pairs}{see \code{\link{csm_lfdr}}.}
\item{...}{ignored.}
}
\details{
The fit maximizes the penalized log-likelihood set by \code{penalty} and
\code{weight_penalty} in \code{\link{csm_control}} (the log-likelihood when
both are 0). Each iteration computes the responsibilities and then updates, in
turn, the weights (closed form; without the weight penalty fitted weights may
be exactly zero), the magnitudes (a convex quadratic program under the box and ordering
constraints, solved blockwise and by ADMM when the ordering constraint binds),
and the correlation matrix (projected gradient steps on
\eqn{\{\Sigma: \mathrm{diag}(\Sigma) = 1, \Sigma \succeq \kappa I\}}). Every
block is accepted only if it does not decrease its conditional objective, so
the objective does not decrease. A start converges when the Aitken
estimate of the remaining increase of the objective falls below
\code{control$tol_loglik} (see \code{\link{csm_control}}).
}
\value{
An object of class \code{csm_fit} with components \code{params} (the selected
parameter), \code{loglik} (average log-likelihood), \code{objective}
(penalized average log-likelihood, the criterion used to select among
starts), \code{selected},
\code{converged}, \code{runs} (one row per start), \code{traces}
(objective paths), \code{spec} (with \code{b_mu} filled in), \code{n}
and \code{control}. \code{predict} returns local false discovery rates.
}
\examples{
spec <- csm_spec(3, b_mu = 10)
lam <- c(0.8, 0.12, 0.05, 0.03)   # probabilities of support sizes 0..3
p <- csm_params(spec, varpi = lam[spec$grp_size + 1] / choose(3, spec$grp_size),
                mu = c(2, 2.5, 3), Sigma = 0.8 * diag(3) + 0.2)
set.seed(1)
train <- csm_simulate(3000, p, spec)$Z
fit <- csm_fit(train, csm_spec(3, b_mu = 10), starts = "null")
fit
}
