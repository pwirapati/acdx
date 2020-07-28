\name{acdx}
\alias{acdx}
\title{Aggregated Cells Differential Expression}
\description{Gene-by-gene meta-regression analysis of aggregated single-cell profiles}
\usage{
  acdx(Y,Zf,X1,G1,Xf=NULL,Xc=NULL,w=NULL,
    iter_tol=1e-5,bt_tol=-0.1,inv_tol=1e-12,
    iter_max=30,bt_max=20,wssize_min=4096/8,
    y_eps=NULL, init=TRUE, keepY=FALSe,keepE=FALSE,
    bootid=NULL,nboot=99,bootpar="beta" )
  }
\arguments{
  \item{Y}{Aggregated cells array}
  \item{Zf}{Scaling-factor groups}
  \item{X1}{Intercept groups}
  \item{G1}{Dispersion groups}
  \item{Xf}{Explanatory factors}
  \item{Xc}{Continuous-valued explanatory variables}
  \item{w}{Weights per aggregate}
  \item{iter_tol}{Convergence tolerance for the overall iterations}
  \item{bt_tol}{Backtracking tolerance}
  \item{inv_tol}{Inversion tolerance}
  \item{iter_max}{Maximum number of overall iterations}
  \item{bt_max}{Maximum number of backtracking iterations}
  \item{wssize_min}{Minimum per-thread temporary workspace size}
  \item{y_eps}{Data shift, default to half of the smallest non-zero value}
  \item{verbose}{Verbosity: quiet (0), one line per fit (1) one line per
    iterations}
  \item{init}{Use default initializer}
  \item{keepY}{Return a copy of \code{Y} in the output object}
  \item{keepE}{Return a copy of the expected mean and variance
   per data point, in the same format as \code{Y}}
  \item{bootid}{labels for bootstrap resampling units}
  \item{nboot}{number of bootstrap resampling}
  \item{bootpar}{parameters whose bootstrap values are to be kept; any
    subset of \code{c("alpha","beta","phi")}}
  }
\details{
  Fit a regression model to aggregated cells data. This a gamma meta-GLM
  model where the error is the sum of given aggregate summary variance
  and the multiplicative gamma error between aggregates. Logarithmic
  link is used to related the mean parameters.
}
\value{
  \item{alpha}{Coefficients for scaling factors, in log scale}
  \item{beta}{Gene-wise effects corresponding to the matrices \code{X1},
    \code{Xf} and \code{Xc}}
  \item{phi}{Gene-wise dispersion parameters}
  \item{status}{A vector of iteration status, number of iterations
    and number of backtracking}
}
\seealso{
  \code{\link{aggr}}
}
\examples{
  acdx()
}
