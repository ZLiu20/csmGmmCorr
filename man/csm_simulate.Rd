\name{csm_simulate}
\alias{csm_simulate}
\title{Simulate from a correlated csmGmm}
\description{Draws configurations, components and statistics.}
\usage{csm_simulate(n, params, spec)}
\arguments{
\item{n}{number of features.}
\item{params}{a \code{\link{csm_params}} object.}
\item{spec}{the matching \code{\link{csm_spec}}.}
}
\value{A list with \code{Z} (statistics), \code{H} (signed configurations) and
\code{group} (support and component index).}
