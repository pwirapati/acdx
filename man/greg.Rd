\name{greg}
\alias{greg}
\title{Aggregate single-cell expression profiles}
\description{Produce aggregate expression profiles from raw counts}
\usage{
  greg( x, cell_sample, cell_type, cell_names=NULL,
    s2_0 = 1/12, sep="\t", verbose=1 )
}
\arguments{
  \item{x}{Raw counts in sparse dgCMatrix format, with genes as rows
  and cells as columns. Alternatively, if it is a single string, it
  is interpreted as local file name containing text table (which
  can be compressed).}
  \item{cell_sample}{Sample identifiers of the cells.}
  \item{cell_type}{Cell type labels of the cells.}
  \item{cell_names}{Cell identifiers to be matched by the names in \code{x}
    (should have exactly the same length as \code{cell_sample}).
    If \code{NULL}, then the data and annotations are assumed to
    have the same length and in the same order. If the names
    is missing in \code{x} but \code{cell_names} is supplied,
    then it is considered an error. The both names are
    present, only the intersection are used, silently dropping
    the mismatches, and the ordering is not important.}
  \item{s2_0}{Baseline variance, with the default 1/12 equals to that of
    the uniform distribution in an interval of width one.}
  \item{sep}{Column separator character; relevant only when \code{x} is
    a file of text table.}
  \item{verbose}{Output progress messages when performing the aggregation}
}
\details{
  The summary for each gene in an aggregate is the mean of raw counts,
  accompanied by an estimate of the variance of the mean (i.e., squared
  standard error). This variance is estimated by the sample variance
  plus \code{s2_0}, divided by the number of cells minus one.

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
