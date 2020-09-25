\name{bremt}
\alias{bremt}
\title{Bootstrap Resampling Multiple-Testing}
\description{Calculate various significance testing statistics
  from a genome-wide scores with bootstrap replicates.}
\usage{
  bremt(theta, T0=0.05)
}
\arguments{
  \item{theta}{A matrix of genes (rows) versus bootstraps (columns)
    containing the score of interest (bigger value is better).
  The original (non-bootstrap) values are in the first column.}
  \item{alpha}{Significance level threshold.}
}
\details{
  This function computes per-comparison error rate (PCER),
  family-wise error rate (FWER)
  and false-discovery rate (FDR) from a set of gene-wise
  scores with bootstrap replicates.

  FWER and FDR are calculated using algorithm in box 2 and 5,
  respectively, of
  Ge, Dudoit, Speed (2003) Test 12:1.
}
\value{
  A numeric matrix with genes as the row (in the original input
  order) and columns:
  \item{theta}{The first column of the original input, typically the
    log ratios or any functional of the fitted parameters.}
  \item{SE}{Bootstrap standard error, after centering using
    the original value (not the mean).}
  \item{T}{\code{theta/SE}}
  \item{PCER}{Bootstrap one-sided p-value for each gene. The null
    distribution is from the centered bootstrap replicates.}
  \item{FWER}{Multiple testing p-value based on the step-down
    MaxT method}
  \item{FDR}{False discovery rate.}
}
\examples{
  ## see the 'quick tutorial'
}
