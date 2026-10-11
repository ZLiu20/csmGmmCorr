\name{csm_incongruity}
\alias{csm_incongruity}
\title{Incongruity diagnostic}
\description{
Counts comparable pairs in the prespecified order of \code{\link{csm_adjust}}
and the pairs in which the stronger feature has the larger score.
}
\usage{
csm_incongruity(Z, scores, target = c("pleiotropy", "replication", "mediation"),
                t = 1L, max_n = 5000L, tolerance = 0, seed = NULL)
}
\arguments{
\item{Z}{statistics.}
\item{scores}{scores to check.}
\item{target, t}{determine the order.}
\item{max_n}{if there are more rows, a random subset of this size is used.}
\item{tolerance}{score differences up to \code{tolerance} are not counted.}
\item{seed}{optional seed for the subset.}
}
\value{A list with \code{n}, \code{comparable_pairs}, \code{violations} and
\code{violation_rate} (\code{NA} if no pair is comparable).}
