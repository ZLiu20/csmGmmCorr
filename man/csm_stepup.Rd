\name{csm_stepup}
\alias{csm_stepup}
\title{Step-up rule on the average local false discovery rate}
\description{
Rejects the largest number \eqn{r} of features with the smallest scores whose
average score is at most \eqn{q} (Sun and Cai, 2007). Ties are broken by the
feature index.
}
\usage{csm_stepup(lfdr, q = 0.1)}
\arguments{
\item{lfdr}{scores in \eqn{[0,1]}.}
\item{q}{level.}
}
\value{A list with \code{rejected} (logical), \code{n_rejected},
\code{cutoff} (largest rejected score) and \code{average}.}
\examples{
csm_stepup(c(0.01, 0.5, 0.02, 0.2, 0.9), q = 0.1)$rejected
}
