print.ac <- function(x,...)
{
  cat("Aggregated cells from", x$source,"\n")
  cat("number of cells = ", sum(x$N),"\n",sep="")
  cat("number of samples = ",nrow(x$N),"\n",sep="")
  cat("number of cell types = ",ncol(x$N),"\n",sep="")
  cat("number of aggregates = ",length(x$N),", (non-empty = ",sum(x$N > 1),")\n",sep="")
  cat("number of genes = ",dim(x$y)[3]," = ", 
      sum(x$gene_mean > 0)," (some positive) + ",
      sum(x$gene_mean==0)," (all zeroes))\n",sep="")
}

plot_ac_genesum <- function(ac,cex=.25,pch=20,n_names=20,...)
{
  j <- (ac$gene_mean0 > 0)
  gene_mean <- ac$gene_mean0[j]
  gene_cv <- ac$gene_sd0[j]/gene_mean
  rcv <- range(gene_cv)
  wcv <- (rcv[2]-rcv[1])/n_names
  plot( log10(gene_mean), gene_cv,
   xlab="gene mean", ylab="gene CV", axes=FALSE,
   cex=cex,pch=pch, ylim=c(0,rcv[2]),
    ,...)
  axis(side=2)
  axis(side=1,at=-5:5,labels=10^{-5:5})
  
  gbin <- .bincode(gene_cv,seq(rcv[1]-wcv,rcv[2]+wcv, wcv))
  g <- tapply( gene_mean, gbin, function(u) names(which.max(u)))
  text( log10(gene_mean[g]),gene_cv[g],labels=g, adj=c(-.2,.5),xpd=TRUE)
}

plot_ac_aggrsum <- function( ac, o=NULL, s=1e-3, pad=2, ... )
{
  if(is.null(o)) o <- 1:nrow(ac$N)
  xlim <- c(1-pad,length(ac$N)+pad)

  on.exit(layout(1))
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar),add=TRUE,after=FALSE)
  oldmar <- par('mar')

  layout(cbind(c(1,2,3)),widths=c(1),heights=c(1,1,1.5))
  par(cex=1)

  par(mar=oldmar+c(-5,0,-2,0))
  plot( log(.5+c(ac$N[o,])), xlim=xlim, xaxs="i", ylab="cell counts",
    axes=FALSE, frame=TRUE, xlab=NA )
  axis(side=2,at=log(.5+c(0,1,10,100,1000)),labels=c(0,1,10,100,1000),las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)

  Ns <- rep(rowSums(ac$N[o,]),ncol(ac$N))
  plot( log((.01+c(ac$N[o,]))/(.01+Ns-c(ac$N[o,]))), xlim=xlim, xaxs="i", ylab="cell proportion",
    axes=FALSE, frame=TRUE, xlab=NA )
  axis(side=2,at=log((+c(.001,.01,.1,.25,.5,.75,.9))/(c(.999,.99,.9,.75,.5,.25,.1))),
    labels=c("0.1%","1%","10%","25%","50%","75%","90%"),las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)


  par(mar=oldmar+c(4,0,-2,0))
  plot.new()
  plot.window(xlim=xlim,xaxs="i",ylim=log(s+c(0,1000)))
  boxplot.matrix( log(s+apply(ac$y[1,o,ac$gene_mean > 0,],2,c)),
    border=gray(.5), 
    use.cols=FALSE,outline=F,add=TRUE,axes=FALSE,frame=TRUE
    )
  title(ylab="mean raw counts")
  ti <- c(0,unlist(sapply(-3:4,function(p) c(1)*10^p)))
  axis(side=2,at=log(s+ti),labels=ti,las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)

  L <- list( rep( rownames(ac$N)[o],ncol(ac$N)), 
    c( sapply(colnames(ac$N),function(u) rep(u,nrow(ac$N)))))
  dtext( L, side=1,line=0, bspace=1,las=c(2,0))

}

