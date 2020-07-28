\name{quest}
\alias{quest}
\title{Question differentially expressed genes}
\description{Rank genes according to an interest functional of the
  model parameters}
\usage{
  quest(fit, ifun, siglevel=0.05, sort.by=NULL)
}
\arguments{
  \item{fit}{\code{acdx} object of fitted model, with bootstrap data}
  \item{ifun}{A function that map the parameters of each gene
    into a scalar value. The function will be called with arguments
  \code{alpha=},\code{beta=} and \code{phi=} (some of which
  can be \code{null} if not present in the fitted object)}
  \item{siglevel}{The cutoffs for constructing FDP quantiles}
  \item{sort.by}{Sorting of the \code{top} table in the output. Default
    to bootstrap t-statistics, can be changed to \code{'fc'} for
    fold change}
}
\details{
}
\seealso{
}
\examples{
}
