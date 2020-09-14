# bootstrap resampling multiple testing
bremt <- function(theta,alpha=.05)
{
  M <- nrow(theta)
  B <- ncol(theta)
  r <- .C("bremt",
    dims=as.integer(c(M,B)),
    theta=as.double(theta),
    alpha=as.double(alpha),
    SE=double(M),
    PCER=double(M),
    FWER=double(M),
    FDP=double(M),
    NAOK=TRUE,
    PACKAGE="acdx"
    )
  dim(r$theta) <- dim(theta)
  out <- cbind(
    theta=theta[,1],SE=r$SE,'T'=r$theta[,1],
    PCER=r$PCER,
    FWER=r$FWER,
    FDP=r$FDP)
  rownames(out) <- rownames(theta)
  out
}
