\name{csm_adjust}
\alias{csm_adjust}
\title{Monotone local-fdr adjustment}
\usage{
csm_adjust(z, lfdr,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, q = 0.1, all = FALSE)
}
\arguments{
\item{z}{Matrix of z-statistics.}
\item{lfdr}{Unadjusted local false discovery rates.}
\item{target, t}{Ordering associated with the scientific target.}
\item{q}{Decision level used in adaptive mode.}
\item{all}{Adjust every row rather than only the decision-relevant subset.}
}
\value{Adjusted scores, rejection indicators, and exactness indicators.}
\description{Construct the least observed-set majorant under the selected
scientific evidence order.}
\details{
The adjusted value for an observation is the largest raw score among the
observed rows that it dominates in the target-specific partial order.  The
adjusted score is therefore never smaller than the raw score.  With
\code{all = TRUE}, every row is adjusted and the resulting observed set has
zero score-order violations up to numerical tolerance.

With \code{all = FALSE}, the implementation evaluates only the subset needed
to determine the decision at level \code{q}.  This can reduce work while
preserving the rejection set; entries outside that subset are not a complete
global majorant.  The operation changes posterior scores and should be
reported separately from raw model inference.
}
