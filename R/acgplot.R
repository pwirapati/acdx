# aggregated-cell gene-wise plot
#
# TODO
# - ordering for samples, ordering for ctypes
# - multiple genes on separate panels (with common x-axis labels)
# - '%' or '#' in gene names to plot logit(percentage) or log(counts)
# - mean lines if 'fit' is supplied
# - auto choice of variance stabilizing parameter (for geometric placement)
#
acgplot <- function( ac, gene,
  o=NULL,
  sample_tags = NULL,
  fit=NULL,
  adjusted=FALSE,
  conf.bar=FALSE,
  s=.5,
  pch=20,
  col.dot=NULL,
  N.dot=TRUE,
  cex.dot=0.5,
  lwd.conf=1,
  lty.conf=1,
  legend.key=NULL,
  pbulk=F,
  ...     # options to 'plot.default'
)
{
  # pardefault <- par(no.readonly=TRUE)
  # on.exit(par(pardefault))

  if(is.null(o)) o <- 1:nrow(ac$N)
  else if(!is.numeric(o) || min(o) < 1 || max(o) > nrow(ac$N))
    stop("o= is not properly formatted.")

  g <- gene[1]
  if( g == '#' )
    {
    y <- c(ac$N[o,])
    v <- y
    }
  else if( g == '%' )
    {
    rs <- rowSums(ac$N)
    y <- c(ac$N[o,])
    y <- y/(rep(rs,ncol(ac$N))-y)
    v <- 0
    }
  else
    {
    v <- c(ac$y[2,o,g,])
    y <- ifelse(is.finite(v),c(ac$y[1,o,g,]),NA)
    if(pbulk)
      {
      y <- y * c(ac$N)
      s <- s * mean(ac$N)
      }
    }

  if(adjusted && gene != '#' && gene != '%')
    {
    if(is.null(fit)) stop("fit object not supplied")
    y <- y/exp(c(fit$alpha[o,,1]))
    }
  n_sample <- length(o)
  n_ctype <- dim(ac$y)[4]
  if( !is.null(sample_tags))
    {
    sample_tags <- cbind(sample_tags)
    n_tags <- ncol(sample_tags)
    }
  else
    n_tags <- 0

  if(is.null(col.dot)) col.dot <- rep(1,length(y))
  else
    {
    if(length(col.dot)==1) col.dot <- rep(col.dot[1],length(y))
    col.dot <- col.dot[o]
    }

  if(!(N.dot)) cex <- cex.dot * rep(1,length(y))
  else cex <- cex.dot * log10(1+c(ac$N[o,]))

  if(!is.null(legend.key))
    {
    m4 <- 1+max(strwidth(names(legend.key),"inches"))/par("cin")[2]
    mar <- par("mar")
    if( mar[4] < m4 )
      {
      mar[4] <- m4
      par(mar=mar)
      }
    }

  if(N.dot)
    {
    m4 <- 3+max(strwidth("N=1000","inches"))/par("cin")[2]
    mar <- par("mar")
    if( mar[4] < m4 )
      {
      mar[4] <- m4
      par(mar=mar)
      }
    }

  plot( log(y+s), axes=FALSE,frame=TRUE, xlab="",
        ylab=paste(ifelse(adjusted,"adjusted","raw"),"counts"),main=g,
        xaxs="i",xlim=c(.5,length(y)+.5),
        col=col.dot,cex=cex,pch=pch,... )

  if(conf.bar)
    {
    col.bar <- rep(col.dot,n_ctype)
    cv <- sqrt(v)/(y+s)
    ub <- exp( log(y+s) + 1.96*cv )
    lb <- exp( log(y+s) - 1.97*cv )
    for(i in 1:length(y))
      lines( c(i,i), log(c(ub[i],lb[i])), col=col.bar[i],lty=lty.conf,
      lwd=lwd.conf)
    }

  p <- ceiling(log10(max(y,na.rm=T)))
  axt <- c(0,c(sapply((p-3):(p+1),function(r) c(.2,.5,1)*10^r)))
  axis(side=2,at=log(axt+s),labels=axt,las=2)

  # cell type labels and partition
  abline(v=n_sample * (1:(n_ctype-1)) + .5,lty=2)
  mtext(side=1,line=n_tags+1,at=n_sample*seq(0.5,n_ctype,1),text=dimnames(ac$y)[[4]])

  if(!is.null(sample_tags))
    for(k in 1:n_ctype)
      for(a in 1:n_tags)
        gmtext(sample_tags[o,a],line=a,off=(k-1)*n_sample)

  if(!is.null(legend.key))
    {
    legend(par('usr')[2],par('usr')[4], xpd=T,
        legend=names(legend.key),xjust=0,yjust=1,
        col=legend.key,pch=20,pt.cex=2,bty="n")
    }

  if(N.dot)
    {
    legend(par('usr')[2],par('usr')[3], xpd=T,
        legend=c("N=1000","N=100","N=10"),xjust=0,yjust=0,
        col=1,pch=20,pt.cex=cex.dot*log10(1+c(1000,100,10)),bty="n")
    }
}
