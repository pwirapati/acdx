\name{tab2spar}
\alias{tab2spar}
\title{Read a tab-delimited table into a sparse matrix}
\description{Read a tab-delimited table of densely formatted numerical
  data with many of zeroes
  into a sparse matrix object \code{dgCMatrix}}
\usage{
  tab2spar(con, transpose=FALSE)
}
\arguments{
  \item{con}{Input as R connection object, which can be a name to a local (compressed)
    file}
  \item{transpose}{Transpose the matrix}
}
\value{
  An object of \code{dgCMatrix} class (see the \code{Matrix} package)
}

