\name{csm_incongruity}
\alias{csm_incongruity}
\title{Count violations of a target ordering}
\usage{
csm_incongruity(z, scores,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, max_n = 5000L, tolerance = 0, seed = NULL)
}
\arguments{
\item{z}{Matrix of z-statistics.}
\item{scores}{Scores ordered so that smaller values are stronger evidence.}
\item{target, t}{Scientific ordering.}
\item{max_n}{Maximum number of rows included in the diagnostic.}
\item{tolerance}{Numerical score tolerance.}
\item{seed}{Seed used when subsampling.}
}
\value{Numbers of sampled rows, comparable pairs, violations, and their rate.}
\description{Counts comparable ordered pairs whose scores move in the wrong direction.}
