\name{csm_adjust}
\alias{csm_adjust}
\title{Monotone adjustment of local false discovery rates}
\description{
Replaces each score by the largest score among the features that are at least
as strong in a prespecified order,
\eqn{(M s)_i = \max\{s_j : j \mbox{ at least as strong as } i\}}{(M s)_i = max{s_j : j at least as strong as i}}, and applies the
step-up rule to the adjusted scores. This is the smallest adjustment that
dominates the scores and is nonincreasing in the order; it rejects no more
features than the unadjusted scores.
}
\usage{
csm_adjust(Z, lfdr, target = c("pleiotropy", "replication", "mediation"),
           t = 1L, q = 0.1, all = FALSE)
}
\arguments{
\item{Z}{statistics of the entire test set.}
\item{lfdr}{their local false discovery rates.}
\item{target}{determines the order: coordinatewise absolute magnitudes within
the same open sign orthant for pleiotropy and mediation (vectors with a zero
coordinate are compared only with themselves), and the directional branches
of Supplementary Section S4 for replication.}
\item{t}{threshold; used by the replication order.}
\item{q}{level.}
\item{all}{if \code{TRUE}, compute every adjusted score. If \code{FALSE}
(default), compute adjusted scores exactly only for the features with the
smallest unadjusted scores, enlarging this set until the decisions are
certified to equal those of the full adjustment; other entries of
\code{adjusted} are \code{NA}.}
}
\details{
Because adjusted scores are never below the unadjusted ones, the decisions
computed with exact adjusted scores on a set \eqn{C} and unadjusted scores
elsewhere equal the full decisions whenever every rejected feature lies in
\eqn{C} and its adjusted score is smaller than every unadjusted score outside
\eqn{C}.
}
\value{A list with \code{adjusted}, \code{rejected}, \code{n_rejected} and
\code{exact} (features whose adjusted score was computed).}
