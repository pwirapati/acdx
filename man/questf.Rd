\name{questf}
\alias{questf}
\title{Question differential expression patterns.}
\description{Rank genes according to an interest functional of the
  model parameters}
\usage{
  questf(fit, ifun)
}
\arguments{
  \item{fit}{\code{acdx} object of fitted model, with bootstrap data}
  \item{ifun}{A function that map the parameters of each gene
    into a scalar value. The function will be called with arguments
  \code{alpha=},\code{beta=} and \code{phi=} (some of which
  can be \code{null} if not present in the fitted object)}
}
