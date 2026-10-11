\name{csm_simulate}
\alias{csm_simulate}
\title{Simulate correlated csmGmm observations}
\usage{csm_simulate(n, params, spec)}
\arguments{
\item{n}{Number of observations.}
\item{params}{Validated mixture parameters.}
\item{spec}{Model specification.}
}
\value{A list containing simulated statistics, signed states, and component groups.}
\description{Draws signed configurations and Gaussian observations from a specified model.}
