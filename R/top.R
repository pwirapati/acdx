top <- function(
  mtp,
  sort.by="T",
  aux=NULL,
  n=10,
  decreasing=T,
  ...    # other options for order
  )
{
  M <- attr(mtp,"M")
  mtp <- data.frame(mtp)
  if(!is.null(aux)) mtp <- cbind(mtp,aux)
  
  revdir <- ifelse(grepl("^-",sort.by),T,F)
  sort.by <- sub("^-","",sort.by)

  notfound <- setdiff(sort.by,colnames(mtp))
  if(length(notfound) > 0)
    stop("sort key(s) not found: ",notfound)

  decreasing <- xor( rep(decreasing,length(sort.by)),revdir)

  o <- do.call("order",
    c(mtp[,sort.by,drop=FALSE],list(decreasing=decreasing,...)))
  if(!is.numeric(n) || n < 1 || n > nrow(mtp))
    n <- nrow(mtp)

  structure(
    mtp[o[1:n],],
    M=M,
    class=unique(c("top",class(mtp)))
    )
}

plot.top <- function( ttab, alpha=0.05,
  xlim=NULL, xlim.frac=0.25,
  lwd=1,col=c("red3","blue3",gray(.333)), ... )
{
  nPCER <- sum(ttab$PCER <= alpha)
  nFWER <- sum(ttab$FWER <= alpha)
  nFDR <- sum(ttab$FDR <= alpha)
  n <- length(ttab$FDR)

  if( is.null(xlim))
    {
    if( nFDR < xlim.frac * n)
      xlim <- c(1,min(n,round(nFDR/xlim.frac)))
    else xlim <- c(1,n)
    }
  else if(is.na(xlim[1]))
    xlim <- c(1,n)
  
  if( xlim[2] < 5 ) xlim[2] <- 5
  if( xlim[2] < 1.1*nFDR ) xlim[2] <- 1.1*nFDR

  plot( ttab$PCER, pch=20,
    xlim=xlim, col=col[3], ylim=c(0,1),
    xlab=paste("rank (of",attr(ttab,"M"),"tests)"),ylab="p-value",
    lwd=lwd,...)
  lines( ttab$FWER, col=col[1],lwd=lwd )
  lines( ttab$FDR, col=col[2],lwd=lwd )
  abline(h=alpha,lty=2,lwd=lwd)
  abline(v=c(nFWER,nFDR),lty=2,col=col[-3],lwd=lwd)

  if(nFWER < nFDR) tadj <- c(1,0,0)
  else tadj <- c(0,1,0)
  tat <- c(nFWER,nFDR,xlim[2])
  mtext(text=c(nFWER,nFDR,nPCER),side=3,at=tat,adj=tadj,col=col)
  mtext(text=paste("\u2264", alpha),side=3,line=1,at=tat,adj=tadj,col=col)
  mtext(text=c("FWER","FDR","PCER"),side=3,line=2,at=tat,adj=tadj,col=col)
}
