
log_axt <- c("0.001","0.002","0.005","0.01","0.02","0.05","0.1","0.2","0.5",
  "1","2","5","10","20","50","100","200","500","1000")

logit_axt <- c("0.01","0.05","0.1","0.25","0.5","0.75","0.9","0.95","0.99")

acplot <- function( ac , gene, o=NULL, xf=function(x)log(x+.5),
  col.key = NULL, col.map = NULL,
  pch=20, cex=NULL,bar=FALSE,
  xt.key=NULL,
  xtv=0,xtl=0,
  xt.key2=NULL,
  xtv2=0,xtl2=2,
  ...
)
{
  gene <- gene[1]
  nac <- nrow(ac$a)
  if(is.null(o)) o <- 1:nac
  else nac <- length(o)

  if(!is.null(col.key) && !is.null(col.map) )
    col <- col.map[ ac$a[o,col.key] ]
  else
    col <- rep(1,nac)
  
  if(is.null(cex))
    cex <- 1 + sqrt(ac$a$N/max(ac$a$N))*2

  if(!is.null(xt.key))
    xt <- rle(ac$a[o,xt.key])
  if(!is.null(xt.key2))
    xt2 <- rle(ac$a[o,xt.key2])
  
  par(mar=c(5+xtv+xtv2,4,4,2))
  plot( xf(ac$y[1,o,gene]),
    axes=FALSE, xaxs="i",xlim=c(0,nac+1),
    xlab="",ylab="ave. counts", main=gene, 
    col=col,pch=pch,cex=cex,
    ... )
  axis(side=2,at=xf(as.numeric(log_axt)),labels=log_axt,las=2)
  if(!is.null(xt.key))
    {
    mtext( xt$values, at=cumsum(xt$lengths)-xt$lengths/2,side=1,
      line=1.5,las=ifelse(xtv > 0,2,1))
    if(xtl) abline(v=0.5+cumsum(xt$lengths),lty=xtl)
    }
  if(!is.null(xt.key2))
    {
    mtext( xt2$values, at=cumsum(xt2$lengths)-xt2$lengths/2,side=1,
      line=3+xtv,las=ifelse(xtv2 > 0,2,1))
    if(xtl2) abline(v=0.5+cumsum(xt2$lengths),lty=xtl2)
    }
  if(!is.null(ac$a$dataset) && length(unique(ac$a$dataset))==1 )
    mtext( side=3, at=0,adj=0,line=2,
      text=paste("dataset: ",ac$a$dataset[1]))
  if(bar==TRUE)
  {
    for(i in 1:nac)
      {
      u <- ac$y[1,o[i],gene]
      cv <- ac$y[2,o[i],gene]/(.001+u^2)
      d <- exp(sqrt(cv)*2)
      lines( c(i,i), xf(c(u*d,u/d)),col=col[i])
      }
  }
}

