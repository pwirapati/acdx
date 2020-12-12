\name{top}
\alias{top}
\title{Show top-ranking genes}
\description{Sort, show and select top-ranking genes from
  a significance testing summary matrix.}
\usage{
  top( 
    mtp,
    sort.by="T",
    aux=NULL,
    n=Inf,
    decreasing=T,
    ...    # other options for sort.list
    )
}
\arguments{
  \item{mtp}{Multiple-testing objects, such as the result
    of \code{\link{bremt}}}
  \item{sort.by}{Sorting rules. Can be a vector of selected column names
    in \code{sigsum} and \code{aux}. Ordering is prioritized
    according to the given order (left first). A prefix \code{"-"} to
    the name reverses the direction of the sort.}
  \item{aux}{Additional columns to be appending the results. It has
    to be a dataframe or matrix in the same gene order as \code{sigsum}.
    It can contain non-numeric values. In addition to be displayed
  together, the columns of \code{aux} can be used in the sorting keys.}
  \item{n}{The number of rows shown. Any values outside \code{1:nrow(sigsum)}
   (including non-numeric or NA) is interpreted as all rows.}
  \item{decreasing}{The default sorting direction, if the key is not prefixed
    by \code{"-"}.}
  \item{...}{Other options for \code{\link{order}}}
}
\examples{
  ## see the 'quick tutorial'
}
