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
    xt.key=NULL,
    xtv=0,xtl=0,
    xt.key2=NULL,
    xtv2=0,xtl2=2,
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
  \item{pch}{Point type, see \code{\link{plot.default}}}
  \item{cex}{The size of the points, defaults to semi-quantitative
    function of the number of cells per aggregate}
  \item{bar}{Show 95-percent confidence bar of the aggregate
    data, based on normal approximation
    in the log scale}
  \item{xt.key}{Horizontal axis label (a column name in \code{ac}
  object annotation)}
  \item{xtv}{If zero, print in horizontal direction If greater than
   zero it correspond to the number of lines allocated for the label}
  \item{xtl}{vertical lines separating the groups}
  \item{xt.key2}{2nd horizontal axis label}
  \item{xtv2}{As above but for the second labels}
  \item{xtl2}{As above for for the second labels}
  \item{...}{Other options to \code{plot}}
}
\details{
}
\seealso{
}
\examples{
data(KUL3)
acplot(KUL3,"EPCAM",
  xt.key="region",xt.key2="ctype",
  col.key="region",col.map=c(core="purple",border="violet",normal="darkgray"))
}
