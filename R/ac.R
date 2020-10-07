print.ac <- function(ac)
{
  cat("Aggregated cells from", ac$source,"\n")
  cat("number of cells = ", sum(ac$N),"\n",sep="")
  cat("number of samples = ",nrow(ac$N),"\n",sep="")
  cat("number of cell types = ",ncol(ac$N),"\n",sep="")
  cat("number of aggregates = ",length(ac$N),", (non-empty = ",sum(ac$N > 1),")\n",sep="")
  cat("number of genes = ",dim(ac$y)[3]," = ", 
      sum(ac$gene_mean > 0)," (some positive) + ",
      sum(ac$gene_mean==0)," (all zeroes))\n",sep="")
}

plot.ac <- function( ac, genes=NULL, ... )
{
  if(is.null(genes))
    plot_ac_genesum( ac, ... )
  else if(is.na(genes))
    plot_ac_aggrsum( ac, ... )
}

plot_ac_genesum <- function(ac,cex=.25,pch=20,n_names=20,...)
{
  rcv <- range(ac$gene_cv[is.finite(ac$gene_cv)])
  wcv <- (rcv[2]-rcv[1])/n_names
  plot( log10(ac$gene_mean), ac$gene_cv,
   xlab="gene mean", ylab="gene CV", axes=FALSE,
   cex=cex,pch=pch, ylim=c(0,rcv[2]),
    ,...)
  axis(side=2)
  axis(side=1,at=-5:5,labels=10^{-5:5})
  
  gbin <- .bincode(ac$gene_cv,seq(rcv[1]-wcv,rcv[2]+wcv, wcv))
  g <- tapply( ac$gene_mean, gbin, function(u) names(which.max(u)))
  text( log10(ac$gene_mean[g]),ac$gene_cv[g],labels=g, adj=c(-.2,.5),xpd=TRUE)
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
  plot( log((.01+c(ac$N[o,]))/(.01+Ns-c(ac$N[o,]))), xlim=xlim, xaxs="i", ylab="cell counts",
    axes=FALSE, frame=TRUE, xlab=NA )
  axis(side=2,at=log((+c(.001,.01,.1,.25,.5,.75,.9))/(c(.999,.99,.9,.75,.5,.25,.1))),
    labels=c("0.1%","1%","10%","25%","50%","75%","90%"),las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)


  par(mar=oldmar+c(4,0,-2,0))
  plot.new()
  plot.window(xlim=xlim,xaxs="i",ylim=log(s+c(0,1000)))
  boxplot.matrix( log(s+apply(ac$y[1,o,KUL3$gene_mean > 0,],2,c)),
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

