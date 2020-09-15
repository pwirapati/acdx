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
  s=.5,
  col.dot=1,
  ... )    # options to 'plot.default'
}
\arguments{
  \item{ac}{Aggregated cell data object produced by \code{\link{greg}}}
  \item{gene}{Name of the gene to be displayed}
  \item{sample_tags}{A data frame of sample tags to be displayed
    below the x-axis. The rows must match the samples in the \code{ac}
    object.}
  \item{o}{Ordering of the sample. Default to `1:n_sample`.}
  \item{fit}{A fitted model object from \code{\link{acdx}}. This
    allows displaying adjusted/normalized data.}
  \item{adjusted}{If \code{fit=} is defined, show adjusted counts
    instead of raw counts.}
  \item{s}{A constant shift for the log-transformation for plotting.
    This is unrelated to any analysis, and only to allow the
    zero counts to be displayed on the plot, with y-axis still
    matching the original value before shifting.}
  \item{col.dot}{Colors for the dots. The colors can
    be made sample specific and has to be supplied in the
    same order as in the \code{ac} object. If
  \code{o=} is supplied, it will be used.}
  \item{...}{Other options to \code{plot}}
}
\examples{
  ## see the 'quick tutorial'
}