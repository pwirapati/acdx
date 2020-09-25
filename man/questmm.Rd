\name{questmm}
\alias{questmm}
\title{Question differential expression patterns using min-max model.}
\description{Rank genes according to difference between the minimum
  of a set of parameter and the maximum of another sets.}
\usage{
  questmm(
    param,
    showcomb = F,
    default="L",
    H=NULL,
    N=NULL,
    L=NULL,
    sep=",",
    na.rm=FALSE,
    include=NULL,
    exclude=NULL,
    multest_funct=bremt
    )
}
\arguments{
  \item{param}{A parameter object that is a part of \code{\link{acdx}}
   output, typically \code{fit$beta}, but can be other parameters
   or transformed values thereof.}
 \item{showcomb}{If \code{TRUE}, just return the list of possible
   combinations of coefficients and cell types, without doing any
 analysis.}
 \item{default}{The default class for comparions, which is \code{"L"}
   (the "low" group}
 \item{H}{A set of regular expression patterns for the "high" group}
 \item{N}{A set of regular expression patterns for the "neutral" group}
 \item{L}{A set of regular expression patterns for the "low" group}
 \item{sep}{String separator for coeffcient-celltype combinations.
   Can be changed if comma is already used as part of the names.}
 \item{na.rm}{Skip coefficients if they are \code{NA}. The default
   is to make the interest score \code{NA} for any missing coefficient.}
  \item{include}{Vector of gene names to be included. If \code{NULL} (default),
    all are included. If specified, it excludes all but those
  in the set.}
  \item{exclude}{A vector of gene names to be excluded. If \code{NULL} (default),
    none is excluded (from the result of applying \code{include=}).}
  \item{multest_func}{Multiple-testing procedure callback function. It
    should take a matrix of interest scores with
    genes as rows and bootstraps as column.}
}
\details{
  \code{questmm} calculates, for each gene and each bootstrap
 replicate, a score defined as
  the difference (log-ratio) between the minimum values of the "high"
  group and maximum value of the "low" group.

  The groups of coefficients can span across multiple cell types, and
  they can be arbitrarily chosen.

  The resulting scores are passed on a multiple-testing correction
  method \code{\link{bremt}}, which collapse the bootstrap
  replicates into significance test statistics.
}
\seealso{
  \code{\link{bremt}}
}
