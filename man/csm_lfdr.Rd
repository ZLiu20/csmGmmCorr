\name{csm_lfdr}
\alias{csm_lfdr}
\alias{csm_null_components}
\title{Target-specific local false discovery rates}
\usage{
csm_lfdr(object, z,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, pairs = NULL, spec = NULL, block_cells = 5e6, cores = 1L)
csm_null_components(H,
  target = c("pleiotropy", "replication", "mediation"),
  t = 1L, pairs = NULL)
}
\arguments{
\item{object}{A fitted model or parameter object.}
\item{z}{Matrix of z-statistics.}
\item{target}{Scientific composite-null target.}
\item{t}{Required support threshold.}
\item{pairs}{Two-column coordinate-pair matrix for mediation.}
\item{spec}{Required when \code{object} contains parameters only.}
\item{block_cells, cores}{Evaluation controls.}
\item{H}{Matrix of signed configurations.}
}
\value{A local-fdr vector, or a logical null-component indicator.}
\description{Sum posterior mixture mass over the scientific null states.}
\details{
For pleiotropy, the alternative contains configurations with more than
\code{t} active coordinates.  For replication, those coordinates must also
have a common sign.  For mediation, the alternative requires more than
\code{t} active coordinates and at least one complete pair listed in
\code{pairs}.  Consequently, the mediation target can depend on \code{t}; a
complete path alone is not sufficient when its support does not exceed the
threshold.

When \code{pairs} is omitted, even-dimensional statistics are paired as
(1,2), (3,4), and so forth.  Odd dimensions require an explicit two-column
pair matrix.  \code{block_cells} and \code{cores} affect evaluation resources,
not the fitted model or target definition.
}
\examples{
H <- rbind(c(0, 0, 0, 0), c(1, 1, 0, 0), c(1, -1, 1, 0))
csm_null_components(H, "pleiotropy", t = 1)
csm_null_components(H, "replication", t = 1)
csm_null_components(H, "mediation", t = 2)
}
