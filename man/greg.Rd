\name{greg}
\alias{greg}
\title{Aggregate single-cell expression profiles}
\description{Produce aggregate expression profiles from raw counts}
\usage{
  greg( x, cell_sample, cell_type, cell_names=NULL,
    u_0 = NULL, s2_0 = 1/12, 
    sep="\t", verbose=1 )
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
  \item{u_0}{Constant shift from the origin. The default is 0.5
    divided by the average number of cells per aggregate.}
  \item{s2_0}{Baseline constant variance of the data. 
    The default is 1/12 corresponding
    to variance of uniform distribution in an interval with a width
    equals to one.}
  \item{sep}{Column separator character; relevant only when \code{x} is
    a file of text table.}
  \item{verbose}{Output progress messages when performing the aggregation}
}
\details{
  The expression of each gene in an aggregate is summarized
  by a pair of statistics: the mean of raw counts, and its variance
  (i.e., squared standard error, which is the sample variance
  divided by the number of cells). These pairs are the first
  dimension of the array \code{y} in the output object.

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
  \item{source}{A string describing the data source}
  \item{min_nz}{The minimum of non-zero mean expression}
  \item{gene_cv}{Gene-specific coefficient of variations}
  \item{gene_scale}{Gene-specific scale factors}
  \item{aggr_scale}{Aggregate-speicfic scale factors}
}
\examples{
  ## see the 'quick tutorial'
}
