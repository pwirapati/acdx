\name{bremt}
\alias{bremt}
\title{Bootstrap Resampling Multiple-Testing}
\description{Calculate various significance testing statistics
  from a genome-wide scores with bootstrap replicates.}
\usage{
  bremt(theta, alpha=0.05)
}
\arguments{
  \item{theta}{A matrix of genes (rows) versus bootstraps (columns)
    containing the score of interest (bigger value is better).
  The original (non-bootstrap) values are in the first column.}
  \item{alpha}{Significance level threshold.}
}
\details{
  This function computes per-comparison error rate (PCER), family-wise error rate (FWER)
  and false-discovery proportion (FDP) from a set of gene-wise
  scores with bootstrap replicates.
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
  \item{FWER}{Multiple testing p-value based on the MaxT method
    of Westfall-Young (1993)}
  \item{FDP}{False discovery proportion, similar to FDR, but based
  on the ranked bootstrap values instead of theoretical distribution.}
}
\examples{
  ## see the 'quick tutorial'
}
