\name{acplot}
\alias{acplot}
\title{Aggregate Cells Plot}
\description{Forest plot of the aggregated cells summaries for a given
  gene}
\usage{
  acplot( ac , gene, o=NULL, xf=function(x)log(x+.5),
    col.key = NULL, col.map = NULL,
    pch=20, cex=NULL,
    bar=FALSE,
    xt.key=NULL,xt.key2=NULL,
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
  \item{pch}{Point type, see \code{\link{plot.default}}}
  \item{cex}{The size of the points, defaults to semi-quantitative
    function of the number of cells per aggregate}
  \item{bar}{Show 95-percent confidence bar of the aggregate
    data, based on normal approximation
    in the log scale}
  \item{...}{Other options to \code{plot}}
}
\details{
}
\seealso{
}
\examples{
}
