\name{acgplot}
\alias{acgplot}
\title{Aggregated Cells Gene Plot}
\description{Forest plot of the aggregated cells summaries for a given
  gene}
\usage{
acgplot( ac, gene, 
  o = NULL,
  sample_tags = NULL,
  fit=NULL,
  adjusted=FALSE,
  conf.bar=FALSE,
  s=.5,
  pch=20,
  col.dot=NULL,
  N.dot=TRUE,
  cex.dot=0.8,
  lwd.conf=1,
  lty.conf=1,
  las=c(0,0),
  legend.key=NULL,
  pbulk=F,
  ... )    # options to 'plot.default'
}
\arguments{
  \item{ac}{Aggregated cell data object produced by \code{\link{greg}}}
  \item{gene}{Name of the gene to be displayed}
  \item{o}{Ordering of the sample. Default to `1:n_sample`.}
  \item{sample_tags}{A data frame of sample tags to be displayed
    below the x-axis. The rows must match the samples in the \code{ac}
    object.}
  \item{fit}{A fitted model object from \code{\link{acdx}}. This
    allows displaying adjusted/normalized data.}
  \item{adjusted}{If \code{fit=} is defined, show adjusted counts
    instead of raw counts.}
  \item{conf.bar}{Show confidence bars. (Note: this is
    only approximate 95\% confidence bars, not for inference purpose.)}
  \item{s}{A constant shift for the log-transformation for plotting.
    This is unrelated to any analysis, and only to allow the
    zero counts to be displayed on the plot, with y-axis still
    matching the original value before shifting.}
  \item{pch}{Dot type (default is filled circles)}
  \item{col.dot}{Colors for the dots. The colors can
    be made sample specific and has to be supplied in the
    same order as in the \code{ac} object. If
    \code{o=} is supplied, it will be used.}
  \item{N.dot}{Dot size is proportional to log number of cells}
  \item{cex.dot}{Dot size scale factor}
  \item{lwd.conf}{Line width of confidence bars.}
  \item{lty.conf}{Line type of confidence bars.}
  \item{las}{Group label directionality}
  \item{legend.key}{A vector of dot colors with names corresponding
    to the labels to be shown in the legend. Used in conjuction
    with \code{col.dot=}}
  \item{pbulk}{Show the sum (mean multiplied by the number of cells)}
  \item{...}{Other options to \code{plot}}
}
\examples{
  ## see the 'quick tutorial'
}
