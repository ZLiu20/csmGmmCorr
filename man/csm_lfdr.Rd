\name{csm_lfdr}
\alias{csm_lfdr}
\alias{csm_null_components}
\title{Local false discovery rates for composite nulls}
\description{
Computes the posterior probability that a feature's signed configuration lies
in the null set of the chosen target. Only the null set depends on the
target, so one fit serves all targets and thresholds.
}
\usage{
csm_lfdr(object, Z, target = c("pleiotropy", "replication", "mediation"),
         t = 1L, pairs = NULL, spec = NULL, block_cells = 5e6, cores = 1L)

csm_null_components(H, target = c("pleiotropy", "replication", "mediation"),
                    t = 1L, pairs = NULL)
}
\arguments{
\item{object}{a \code{\link{csm_fit}} object or a \code{\link{csm_params}}
object (then \code{spec} is required).}
\item{Z}{test statistics, a matrix with \eqn{K} columns or one vector.}
\item{target}{\code{"pleiotropy"}: more than \code{t} nonzero effects;
\code{"replication"}: more than \code{t} nonzero effects, all with the same
sign; \code{"mediation"}: at least one complete pair and more than \code{t}
nonzero effects.}
\item{t}{threshold, \eqn{1 \le t < K}.}
\item{pairs}{for mediation, a two-column matrix of the exposure--mediator and
mediator--outcome coordinates of each mediator; by default
\eqn{(1,2), (3,4), \ldots}.}
\item{spec}{a \code{\link{csm_spec}} object.}
\item{block_cells, cores}{see \code{\link{csm_control}}.}
\item{H}{a matrix of signed configurations with entries in \eqn{\{-1,0,1\}}.}
}
\value{
\code{csm_lfdr} returns a numeric vector. \code{csm_null_components} returns
a logical vector that is \code{TRUE} for null configurations.
}
