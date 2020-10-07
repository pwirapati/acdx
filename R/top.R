top <- function(
  mtp,
  sort.by="T",
  aux=NULL,
  n=Inf,
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

plot_top_rank <- function( ttab, alpha=0.05,
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


local_max <- function( x, y, nbins )
{
  rx <- range(x)
  dx <- diff(rx)/nbins
  b <- .bincode( x, seq(rx[1]-dx,rx[2]+dx,dx ))
  names(y) <- 1:length(y)
  as.integer(tapply(y,b,function(u)names(which.max(u))))
}


plot_top_stat <- function( ttab, n_names=40, alpha=0.05, onesided=TRUE, 
  cex=1,pch=20,lwd=2,
  ... )
{
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  mar <- par('mar')
  if(mar[4] < 4) { mar[4] <- 4; par(mar=mar) }

  if( onesided )
    {
    ttab <- ttab[ ttab$theta >= 0, ]
    xmax <- 1.2 * max(ttab$theta)
    ymax <- 1.2 * max(ttab$T)
    }
  else
    stop("two sided not implemented")

  plot( ttab$theta, ttab$T, cex=cex, pch=pch, axes=FALSE,
    frame=TRUE,xlab="exp(theta)",ylab="T",
    xlim=c(0,xmax),ylim=c(0,ymax),
    ... )
  fc <- sapply(-1:3,function(p) c(1,2,5)*10^p,simplify=TRUE)
  axis(side=1, at=log(fc),labels=as.character(fc))
  axis(side=2)
  FWERcutoff <- ttab[ sum(ttab$FWER <= alpha),"T"]
  abline(h=FWERcutoff,col=2,lty=2,lwd=lwd )
  mtext(side=4,at=FWERcutoff,text=paste0("FWER\n",alpha),col=2,las=2,line=1)

  FDRcutoff <- ttab[ sum(ttab$FDR <= alpha),"T"]
  abline(h=FDRcutoff,col=4,lty=2,lwd=lwd )
  mtext(side=4,at=FDRcutoff,text=paste0("FDR\n",alpha),col=4,las=2,line=1)

  gh <- local_max( ttab$T, ttab$theta, n_names)
  gv <- local_max( ttab$theta, ttab$T, n_names)

  g <- union(gh,gv) 
  g <- g[ttab[g,"T"] >= FWERcutoff | ttab[g,"T"] >=FDRcutoff]
  text( labels=rownames(ttab)[g], ttab[g,"theta"],ttab[g,"T"],adj=c(-.1,.5))

}

plot.top <- function( ttab, type="rank", ... )
{
  typeid <- pmatch(type,c("rank","pvalue","stat","volcano"))
  if(typeid ==1 || typeid == 2)
    plot_top_rank( ttab, ... )
  else if( typeid == 3 || typeid == 4 )
    plot_top_stat( ttab, ... )
} 
