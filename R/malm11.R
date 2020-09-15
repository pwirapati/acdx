malm11 <- function( y, w = NULL,
  iter_max=30, epsilon=1e-4,y0=0,verbose=0 )
{
  n <- dim(y)[2]
  m <- dim(y)[3]
  if(is.null(w)) w <- rep(1,n)
  r <- .C("malm11",
    dim=as.integer(c(n,m)),
    y = as.double(y + c(y0,0)),
    w=as.double(w),
    E = as.double(rep(1,2*n*m)),
    alpha=as.double(ifelse(w==0,NA,1)),
    gamma=as.double(rep(1,m)),
    psi=as.double(rep(1,m)),
    L=double(1),
    iopt=as.integer(c(verbose,iter_max)),
    dopt=as.double(c(epsilon)),
    NAOK=TRUE,
    PACKAGE="acdx"
    )
  r$y <- NULL
  r$E <- NULL
  names(r$alpha) <- dimnames(y)[[2]]
  names(r$gamma) <- dimnames(y)[[3]]
  names(r$psi) <- dimnames(y)[[3]]
  r
}