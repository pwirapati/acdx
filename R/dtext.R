dtext <- function(
  L,  # list of labels, from inner- to outer-most
  side=1,
  off=0,
  line=0,
  bspace=1,
  bpos=c(.2,.8),
  las=NULL
  )
{
  if(!is.list(L))
    L <- as.list(cbind(L))

  if(is.null(las))
    las <- rep(0,length(L))

  if(length(las) == 1)
    las <- rep(las,length(L))

  R <- sapply( L, function(labels)
    {
    r <- rle(as.character(labels))
    wmax <- max(sapply(r$values,strwidth,"i"),na.rm=T)
    b <- cumsum(r$lengths)
    m <- length(b)
    a <- c(0,b[-m])+1
    inc <- !is.na(r$values)
    list(labels=r$values[inc],a=a[inc],b=b[inc],wmax=wmax)
    },simplify=FALSE)

  l <- line
  for(i in 1:length(R))
    {
    R[[i]]$line <- l
    mtext(side=side,at=off+(R[[i]]$a+R[[i]]$b)/2,
      text=R[[i]]$labels,las=las[i],
      line=l+bspace)

    if(las[i]==0 || (las[i]==1 && (side==1||side==3)))
      l <- l + 1 + bspace
    else 
      l <- l + R[[i]]$wmax/par('csi') + bspace
    }
  if( l > par('mar')[side] )
    warning("dtext: margin ",side," too small (",par('mar')[side]," < ",l,")",
      call.=FALSE)

  oxpd <- par("xpd")
  par(xpd=TRUE)
  recordGraphics(
    {
    if(side==1||side==3)
      s <- par("cxy")[2]
    else
      s <- xinch(par("csi"))

    for(i in 1:length(R))
      {
      y <- par('usr')[ c(3,1,4,2)[side] ] + 
            c(-1,-1,1,1)[side]*(R[[i]]$line + bspace - bpos) * s

      if(side==1||side==3)
        for(j in 1:length(R[[i]]$a))
          lines( off+c(R[[i]]$a[j],R[[i]]$a[j],R[[i]]$b[j],R[[i]]$b[j]),
            c(y[2],y[1],y[1],y[2]))
      else
        for(j in 1:length(R[[i]]$a))
          lines( c(y[2],y[1],y[1],y[2]),
            off+c(R[[i]]$a[j],R[[i]]$a[j],R[[i]]$b[j],R[[i]]$b[j]))
      }
    },
    list(R=R,side=side,off=off,bspace=bspace,bpos=bpos),
    getNamespace("graphics"))
  par(xpd=oxpd)

  invisible(l)
}
