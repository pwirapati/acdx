\name{aggr}
\alias{aggr}
\title{Aggregate Single-Cell Expression Profiles}
\description{Produce aggregate matrix from raw counts}
\usage{
  aggr( x, cellgroups, v0 = 1/12 )
}
\arguments{
  \item{x}{Raw counts in sparse dgCMatrix format, with genes as rows
  and cells as columns}
  \item{cellgroups}{Labels for grouping the cells into the same aggregate,
    a vector with the same length as the number of columns of \code{x}}
  \item{v0}{Baseline variance, with the default 1/12 equals to that of
    the uniform distribution in an interval of width one}
}
\details{
  The summary for each gene in an aggregate is the mean of raw counts,
  accompanied by an estimate of the variance of the mean (i.e., squared
  standard error). The aggregate data is a 3d array with the pair
  of mean and variance in the first dimension, the aggregate (as
  identified by \code{cellgroups} in the second, and the genes
  as the third.
}
\value{
  A list of:
  \item{annot}{Aggregate annotation which contains the aggregate
    identifier (\code{id}) and the number of cells per aggregate (\code{N}).}
  \item{y}{The aggregate array}
}
\seealso{
  \code{\link{acdx}}
}

\examples{}
