\name{acplot}
\alias{acplot}
\title{Aggregate Cells Plot}
\description{Forest plot of the aggregated cells summaries for a given
  gene}
\usage{
  acplot( ac , gene, o=NULL, xf=function(x)log(x+.5),
    col.key = NULL, col.map = NULL,
    pch=20, cex=NULL,
    xt.key=NULL,xt.key2=NULL,bar=FALSE,
    ...)
}
\arguments{
  \item{ac}{Aggregated cell data produced by \code{\link{aggr}}}
  \item{gene}{Name of the gene to be displayed}
  \item{o}{Ordering of the aggregates in the plot}
  \item{xf}{Function to transform the data and axis}
  \item{col.key}{Color key: the name of a column in the \code{ac}
    object annotation data frame}
  \item{col.map}{Color map: a vector that maps the \code{col.key}
    string into actual R color}
  \item{xt.key}{Horizontal axis label (a column name in \code{ac}
  object annotation)}
  \item{xt.key2}{2nd horizontal axis label}
  \item{...}{Other options to \code{plot}}
}
\details{
}
\seealso{
}
\examples{
}
