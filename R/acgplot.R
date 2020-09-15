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
  s=.5,
  col.dot=NULL,
  ...     # options to 'plot.default'
)
{
  g <- gene[1]
  if(is.null(o)) o <- 1:nrow(ac$N)
  else if(!is.numeric(o) || min(o) < 1 || max(o) > nrow(ac$N))
    stop("o= is not properly formatted.")

  v <- c(ac$y[2,o,g,])
  y <- ifelse(is.finite(v),c(ac$y[1,o,g,]),NA)
  if(adjusted)
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

  plot( log(y+s), axes=FALSE,frame=TRUE, xlab="",
        ylab=paste(ifelse(adjusted,"adjusted","raw"),"counts"),main=g,
        xaxs="i",xlim=c(.5,length(y)+.5),
        col=col.dot[o],... )

  p <- ceiling(log10(max(y,na.rm=T)))
  axt <- c(0,c(sapply((p-2):(p+1),function(r) c(.2,.5,1)*10^r)))
  axis(side=2,at=log(axt+s),labels=axt,las=2)

  # cell type labels and partition
  abline(v=n_sample * (1:(n_ctype-1)) + .5,lty=2)
  mtext(side=1,line=n_tags+1,at=n_sample*seq(0.5,n_ctype,1),text=dimnames(ac$y)[[4]])

  if(!is.null(sample_tags))
    for(k in 1:n_ctype)
      for(a in 1:n_tags)
        gmtext(sample_tags[o,a],line=a,off=(k-1)*n_sample)

}
