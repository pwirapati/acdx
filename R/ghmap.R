ghmap <- function(...) UseMethod("ghmap")

ghmap.ac <- function(ac,genes,
  o=NULL,
  normalize=c(TRUE,TRUE),modnorm=0,
  cluster_genes = TRUE,
  cluster_samples = NULL, cstrat = NULL,
  col.sep="white",lwd.sep=2,
  saturation=.1,
  col=gray.colors(20,start=.95,end=.05),
  tags=NULL,las.tags=NULL,
  ...
  )
{
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  par(mar=c(8,4,4,8))

  n_sample <- dim(ac$y)[2]
  n_ctype <- dim(ac$y)[4]

  if(is.null(o)) o <- 1:nrow(ac$N)
  genes <- rev(intersect(genes,dimnames(ac$y)[[3]]))
  z <- collapse2dim( ac$y[1,o,genes,], c(1,3) )

  if(normalize[1])
    {
    s <- c(ac$aggr_scale[o,])
    s <- ifelse(!is.finite(s),1,s)
    z <- z/s
    }

  if(normalize[2]==TRUE)
    z <- apply(z,2,function(u) u/sqrt(sum(modnorm+u^2,na.rm=TRUE)))
    #z <- apply(z,2,function(u) u/(s+mean(u,na.rm=TRUE)))

  if(cluster_genes)
    {
    cz <- 1-cor(z,use="pairwise.complete.obs")
    cz <- ifelse(!is.finite(cz),0,cz)
    og <- hclust(as.dist(cz),method="average")$order
    z <- z[,og]
    }

  if(!is.null(cluster_samples))
    {
    k <- intersect(cluster_samples,1:n_ctype)
    x <- ac$y[1,,genes,k,drop=F]
    dim(x) <- c(dim(x)[2],dim(x)[3]*dim(x)[4])
    x <- apply(x,1,function(u) u/mean(u,na.rm=TRUE))
    cx <- 1-cor(x,use="pairwise.complete.obs")
    cx <- ifelse(!is.finite(cx),0,cx)
    os <- hclust(as.dist(cx),method="average")$order
    if(!is.null(cstrat))
      os <- order(cstrat,order(os))
    oo <- c(sapply(0:(n_ctype-1),function(i) os+i*n_sample))
    z <- z[oo,]
    }


  image(
    tanh(saturation*z),
    zlim=c(0,1),
    y=seq(1,ncol(z)),
    ylim=c(0,ncol(z))+0.5,
    x=seq(1,nrow(z)),
    xlim=c(0,nrow(z))+0.5,
    col=col,
    useRaster=TRUE,axes=FALSE,xlab=NA,ylab=NA,...)
  rowticks <- pretty(1:ncol(z))
  rowticks[1] <- 1
  rowticks <- c(rowticks,ncol(z))
  axis(side=2,las=2,at=rev(rowticks),labels=rowticks)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + 0.5,col=col.sep,lwd=lwd.sep)

  d <- ceiling(2*strheight("","u"))
  ig <- seq(1,ncol(z),d)
  mtext(side=4,line=.5,at=1:ncol(z),text=colnames(z),las=2)

  L <-  list( c( sapply(colnames(ac$N),function(u) rep(u,nrow(ac$N)))) )
  if( !is.null(tags) )
    {
    if(!is.list(tags)) tags <- list(tags)
    if( any(lapply(tags,length) != nrow(ac$N)))
       stop("tags: wrong length")
    L <- c( lapply(tags,function(u)rep(u,ncol(ac$N))), L )
    if( is.null(las.tags))
      las.tags <- rep(0,length(L))
    else if(length(las.tags)==length(L)-1) 
      las.tags <- c(las.tags,0)
    else
      stop("las.tags: wrong length")
    }
  else
    las.tags <- 0

  acdx:::dtext( L, side=1,line=0, bspace=1,las=las.tags)

  invisible(z)
}


