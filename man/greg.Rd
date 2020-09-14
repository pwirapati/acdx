\name{greg}
\alias{greg}
\title{Aggregate Single-Cell Expression Profiles}
\description{Produce aggregate array from raw counts}
\usage{
  greg( x, cell_sample, cell_type, v0 = 1/12, sep="\t", verbose=1 )
}
\arguments{
  \item{x}{Raw counts in sparse dgCMatrix format, with genes as rows
  and cells as columns. Alternatively, if it is a single string, it
  is interpreted as local file name containing text table (which
  can be compressed).}
  \item{cell_sample}{Sample labels of the cells,
    a vector with the same length as the number of columns of \code{x}}
  \item{cell_type}{Cell type labels of the cells,
    a vector with the same length as the number of columns of \code{x}}
  \item{v0}{Baseline variance, with the default 1/12 equals to that of
    the uniform distribution in an interval of width one.}
  \item{sep}{Separator character; relevant only when \code{x} is
    a file of text table}
  \item{verbose}{Output progress messages when performing the aggregation}
}
\details{
  The summary for each gene in an aggregate is the mean of raw counts,
  accompanied by an estimate of the variance of the mean (i.e., squared
  standard error). This variance is estimated by the sample variance
  plus \code{v0}, divided by the number of cells minus one.

  If the cell count is zero or one, the variance is infinite (the
  mean is set to zero or the single-cell value, respectively).
  This aggregate will be ignored in subsequent analyses, although
  they are still formally represented by the data structure.
}
 
\value{
  A list of:
  \item{N}{A matrix of samples (rows) by cell types (columns), containing
    the number of cells in each aggregate.}
  \item{y}{The aggregate array with four dimensions, corresponding
  to the summaries, samples, genes and cell types.}
}
\examples{
  ## see the 'quick tutorial'
}
