\name{acdx}
\alias{acdx}
\title{Aggregated Cell Differential Expression}
\description{Gene-by-gene meta-regression analysis of aggregated single-cell profiles}
\usage{
  acdx( ac, X, Gi=NULL, n_boot=1,
    id_boot=NULL, seed_boot=0,
    pbulk=FALSE, mean2sum = TRUE,
    norm_method=1,
    u_0=NULL, s2_0 = 1/12, verbose=0 )
  }
\arguments{
  \item{ac}{Aggregated cell object (output of \code{\link{greg}})}
  \item{X}{Design matrix}
  \item{Gi}{Group indicators of dispersion parameters. Default: one
  parameter for all samples.}
  \item{n_boot}{Number of bootstrap resamples}
  \item{id_boot}{Identifier of bootstrap blocks}
  \item{seed_boot}{Random number seed for bootstrapping}
  \item{pbulk}{Use pseudo-bulk approach (negative binomial model)}
  \item{mean2sum}{\code{y} is converted from means to sums, if \code{pbulk=TRUE}.}
  \item{norm_method}{normalization method: 0 = none, 1 = per cell-type, 2 = global}
  \item{u_0}{Small constant added to the data (default: half of smallest
  nonzero value in the entire dataset)}
  \item{s2_0}{Small constant variance.}
  \item{verbose}{Print \code{#} for each
    bootstrap replicate when `verbose=1`. The default is silent.}
}
\details{
Fit a gamma-GLM-like regression with variance components
that include per data point standard errors from cell-level
aggregation.
The model
is analogous to random-effect meta-regression models used
in meta-analysis, but with the variance as quadratic
function of the mean.
The dispersion parameters can be interpreted
as the square of coefficient of variations, modelling relative
(multiplicative) error due to between aggregate heterogeneity,
after subtractin the cell-level variability. 

There are two models fitted step-wise. The first one fits
simultaneous row (sample)
and column (gene) intercepts, without considering the design matrix.
The sample intercepts (which correspond to normalization or size 
factors) are then used as offset in the second model. The
gene intercepts and dispersions are not used for the second model,
but are still kept.

The second model is fitted separately for all combination
of gene, cell type and bootstrap replicates.

Missing aggregates due to zero or one counts are skipped.
The parameters can be \code{NA} due to missing aggregates
or collinearity (which favor the term that appear first in
the design matrix).
}
\value{
\item{gamma}{Gene-specific intercepts from first-step model, in natural
  logarithmic scale.}
\item{alpha}{Sample-specific intercepts from first-step model, in natural
  logarithmic scale.}
\item{psi}{Gene-specific gamma dispersion parameter from the first
  model. Not in logarithmic scale.}
\item{beta}{Coefficients from the second-step model, corresponding
  to columns of \code{X}. In natural
  logarithmic scale.}
\item{phi}{Gene- and group-specific gamma dispersion parameter from the second
  model. Not in logarithmic scale.}
}
\examples{
  ## See `quick tutorial`
}
