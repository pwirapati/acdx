# bootstrap resampling multiple testing
bremt <- function(theta,T0=.2)
{
  M <- nrow(theta)
  B <- ncol(theta)
  r <- .C("bremt",
    dims=as.integer(c(M,B)),
    theta=as.double(theta),
    T0=as.double(T0),
    SE=double(M),
    PCER=double(M),
    FWER=double(M),
    FDR=double(M),
    NAOK=TRUE,
    PACKAGE="acdx"
    )
  dim(r$theta) <- dim(theta)
  out <- cbind(
    theta=theta[,1],SE=r$SE,'T'=r$theta[,1],
    PCER=r$PCER,
    FWER=r$FWER,
    FDR=r$FDR)
  rownames(out) <- rownames(theta)
  out
}
