\name{csm_stepup}
\alias{csm_stepup}
\title{Local-fdr step-up rule}
\usage{csm_stepup(lfdr, q = 0.1)}
\arguments{
\item{lfdr}{Local false discovery rates in \eqn{[0,1]}.}
\item{q}{Target Bayesian false discovery rate.}
}
\value{A list containing rejection indicators, count, cutoff, and selected average.}
\description{Selects the largest ordered prefix whose average local false
discovery rate does not exceed \code{q}.}
