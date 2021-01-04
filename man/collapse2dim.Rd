\name{collapse2dim}
\alias{collapse2dim}
\title{Collapse 2 dimensions in a multi-dimensional array}
\description{Combine 2 dimensiones in a multi-dimension array into one, replacing it
  with the flattened Cartesian product}
\usage{
  collapse2dim(a,cid,pos=1,sep="|")
}
\arguments{
  \item{a}{Array to be collapsed}
  \item{cid}{An ordered pair of indices to collapse}
  \item{pos}{The position of the new dimension in the resulting array}
  \item{sep}{Separator string for the collpsed dimension names}
}
