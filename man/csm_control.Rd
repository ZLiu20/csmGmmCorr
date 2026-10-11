\name{csm_control}
\alias{csm_control}
\title{Control parameters for fitting}
\description{Tolerances and limits of the constrained ECM algorithm.}
\usage{
csm_control(max_iter = 5000L, tol_loglik = NULL, tol_param = Inf,
            accelerate = TRUE, warm_n = 0L, penalty = 30, weight_penalty = 10,
            block_cells = 5e6,
            sigma_steps = 50L, fix_sigma = FALSE, cores = 1L, verbose = FALSE)
}
\arguments{
\item{max_iter}{maximum number of ECM steps per start; a start that reaches
it is recorded as unconverged.}
\item{tol_loglik}{a start converges when the Aitken estimate of the remaining
increase of the average log-likelihood, \eqn{\Delta_t/(1-a_t)} with
\eqn{\Delta_t} the last increment and \eqn{a_t} the ratio of the last two
increments, is at most \code{tol_loglik}. The default \code{NULL} uses
\eqn{0.01/N}, small relative to the \eqn{O(1/N)} deficit allowed by the
approximate-maximization condition.}
\item{tol_param}{optionally, also require the largest absolute parameter
change in the last iteration to be at most \code{tol_param}.}
\item{accelerate}{use SQUAREM extrapolation (Varadhan and Roland, 2008):
after two ECM steps, extrapolate, project onto the parameter region and take
one more ECM step, keeping the result only if it does not lower the
log-likelihood. It changes the path to the maximum, not the model or the
estimator. Like plain ECM it finds a local maximum; on multimodal
likelihoods different paths can end at different local maxima, so compare
log-likelihoods across starts.}
\item{warm_n}{when the number of training features exceeds \code{3 * warm_n},
each start is first fitted on a fixed, evenly spread subset of \code{warm_n}
features and then continued on all features from that point, with the
full-data stopping rule. This changes only the starting point of the full-data
run, not the estimator. Off by default (\code{0}): in our K = 10 benchmark the
subset fit led the full-data run to a slower path and a lower local maximum.}
\item{penalty}{strength \eqn{\lambda \ge 0} of the magnitude penalty. The
fit maximizes \eqn{F(\zeta) - (\lambda/N) P(\mu)}, where \eqn{F} is the
average log-likelihood and
\eqn{P(\mu) = \frac12 \sum_{b,m} \sum_{k \in A_b} (\mu_{b,m,k} - \bar\mu_{|A_b|,k})^2}{P(mu) = 0.5 * sum (mu[b,m,k] - mbar[|A_b|, k])^2},
with \eqn{\bar\mu_{s,k}}{mbar[s,k]} the mean magnitude of coordinate
\eqn{k} over supports of size \eqn{s}. Supports with few training features
are shrunk toward supports of the same size; supports with many features are
essentially unaffected, and the effect vanishes as \eqn{N} grows. The model
is unchanged, and because \eqn{P} is bounded on the parameter region the
penalized maximizer is within \eqn{O(\lambda/N)} of the maximum likelihood.
\code{0} gives maximum likelihood. Roughly, \eqn{\lambda} is the number of
pseudo-features of pooled information each support borrows. The default
\code{30} was chosen by simulation (K = 10, N = 1e5), where it brought the
false discovery proportion close to that of the oracle.}
\item{weight_penalty}{number \eqn{\alpha \ge 0} of pseudo-features given to
each support (and magnitude component) when its weight is divided within its
support size. The fit adds \eqn{(\alpha/N)\sum_g \log(\varpi_g/W_{s(g)})} to
the objective, where \eqn{W_s} is the total weight of supports of size
\eqn{s}; the size totals are unaffected and each share within a size becomes
\eqn{(n_g + \alpha)/(n_s + \alpha J_s)} in feature counts, with \eqn{J_s} the
number of supports of that size. The effect vanishes as \eqn{N} grows.
\code{0} gives maximum likelihood. The default \code{10} was chosen by
simulation together with \code{penalty = 30}.}
\item{block_cells}{memory limit of the plain-R reference implementation;
not used by the compiled E-step.}
\item{sigma_steps}{maximum projected-gradient steps in each correlation
update.}
\item{fix_sigma}{\code{FALSE} (default) estimates the correlation matrix.
\code{TRUE} fixes it at the identity: the independence working model of the csmGmm
of Sun, McCaw and Lin (2025). A correlation matrix (eigenvalues at least
\code{kappa}) fixes it at that matrix, e.g. a plug-in estimate as in their
\code{symm_fit_cor_EM}. The fixed value replaces the correlation matrix of every
starting value. Used for comparisons.}
\item{cores}{number of threads used by the compiled E-step and the local false
discovery rates; features are split into blocks processed in parallel.}
\item{verbose}{print progress messages.}
}
\value{A list.}
\references{
Varadhan, R., and Roland, C. (2008). Simple and globally convergent methods for
accelerating the convergence of any EM algorithm. \emph{Scandinavian Journal
of Statistics} 35, 335--353.

Bohning, D., Dietz, E., Schaub, R., Schlattmann, P., and Lindsay, B. G. (1994).
The distribution of the likelihood ratio for mixtures of densities from the
one-parameter exponential family. \emph{Annals of the Institute of
Statistical Mathematics} 46, 373--388.
}
