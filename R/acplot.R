
log_axt <- c("0.001","0.002","0.005","0.01","0.02","0.05","0.1","0.2","0.5",
  "1","2","5","10","20","50","100","200","500","1000")

logit_axt <- c("0.01","0.05","0.1","0.25","0.5","0.75","0.9","0.95","0.99")

acplot <- function( ac , gene, o=NULL, xf=function(x)log(x+.5),
  col.key = NULL, col.map = NULL,
  pch=20, cex=NULL,
  xt.key=NULL,xt.key2=NULL,bar=FALSE,
  ...
)
{
  gene <- gene[1]
  nac <- nrow(ac$annot)
  if(is.null(o)) o <- 1:nac
  if(!is.null(col.key) && !is.null(col.map) )
    col <- col.map[ ac$annot[o,col.key] ]
  else
    col <- 1
  
  if(is.null(cex))
    cex <- 1 + sqrt(ac$annot$N/max(ac$annot$N))*2

  if(!is.null(xt.key))
    xt <- rle(ac$annot[o,xt.key])
  if(!is.null(xt.key2))
    xt2 <- rle(ac$annot[o,xt.key2])


  plot( xf(ac$y[1,o,gene]),
    axes=FALSE, xaxs="i",xlim=c(0,nac+1),
    xlab="",ylab="ave. counts", main=gene, 
    col=col,pch=pch,cex=cex,
    ... )
  axis(side=2,at=xf(as.numeric(log_axt)),lab=log_axt,las=2)
  if(!is.null(xt.key))
    {
    mtext( xt$values, at=cumsum(xt$lengths)-xt$lengths/2,side=1,line=1.5)
    #abline(v=0.5+cumsum(xt$lengths),lty=3)
    }
  if(!is.null(xt.key2))
    {
    mtext( xt2$values, at=cumsum(xt2$lengths)-xt2$lengths/2,side=1,line=3)
    abline(v=0.5+cumsum(xt2$lengths),lty=2)
    }
  mtext( side=3, at=0,adj=0,text=paste("dataset: ",ac$annot$dataset[1]),line=2)
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

