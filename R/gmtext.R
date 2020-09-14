# grouped label for mtext on categorical axis
#
gmtext <- function(text,side=1,las=1,line=1,off=0, col.tag=NULL)
{
  r <- rle(text)
  b <- cumsum(r$lengths)
  m <- length(b)
  a <- ((c(0,b)+c(b,0))/2)[1:m]
  col <- 1
  if(!is.null(col.tag)) col <- col.tag[ r$values ]
  mtext(side=side,at=off+a,text=r$values,las=las,line=line,col=col)
}

