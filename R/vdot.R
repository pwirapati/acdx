vdot <- function(sc,gene,ctype=NULL,cell_ctype=NULL,sampleid=NULL,
  sample_tags = NULL,
  bw=.5,s=0.5,width=.95,
  xlim=NULL,ylim=NULL,xaxs="i",...)
{
  if(length(cell_ctype) == 1 && !is.null(sc$cell[,cell_ctype]))
    cell_type <- sc$cell[,cell_ctype]

  if(is.null(cell_ctype))
    if(is.null(sc$cell[,2])) 
      stop("cell type not specified.")
    else
      cell_ctype <- sc$cell[,2]

  if(is.null(ctype))
    ctype <- unique(cell_ctype)

  if(is.null(sampleid))
    sampleid <- rownames(sc$sample)

  Y <- split( log(s+unname(sc$y[gene[1],])),
              list( factor(cell_ctype,levels=ctype), 
                    factor(sc$cell$sample,levels=sampleid) ), lex.order=T )

  Y <- lapply( Y, function(y) 
    {
    if( length(y) > 2 )
      {
      d <- density(y,bw=bw)
      fy <- d$y[ findInterval(y,d$x) ]
      cbind(runif( d$n, -.5*fy, .5*fy ),y)
      }
    else if( length(y) > 0 )
      {
      cbind(0,y)
      }
    else
      {
      cbind(NA,NA)
      }
    })
  maxw <- max(0,sapply(Y,function(u)max(0,abs(u[,1]),na.rm=T)),na.rm=T)
  n <- length(Y)
  if(is.null(xlim))
    xlim <- c(-.5,n+.5)
  if(is.null(ylim))
    {
    r <- sapply(Y,function(u) suppressWarnings(range(u[,2],na.rm=T)))
    ylim <- c(min(Inf,r[1,],na.rm=T),max(-Inf,r[2,],na.rm=T))
    }
  
  olmar <- par(mar=c(10,5,4,5))
  plot(NULL,xlim=xlim,ylim=ylim,xaxs=xaxs,xlab=NA,ylab="raw counts",
    axes=FALSE,frame=TRUE,main=gene)
  axt <- c(0,c(sapply(0:7,function(p) c(1,2,5)*10^p)))
  axis(side=2,at=log(s+axt),labels=axt,las=2)

  s <- 0.5*width/maxw
  for(i in 1:n)
    points(i+Y[[i]][,1]*s,Y[[i]][,2],...)

  if(is.null(sample_tags)) 
    L <- data.frame( rep(ctype,each=length(sampleid)),stringsAsFactors=FALSE)
  else
    L <- expand.grid(sample_tags,ctype)
  acdx:::dtext( L, side=1,line=0,bspace=1,las=c(2,0))
  abline(v=acdx:::sep(L[,ncol(L)]),lty=2,col="red")

  par(olmar)
  invisible(Y)
}
