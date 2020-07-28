aggr <- function( x, cellgroups, v0 = 1/12 )
{
  cellgroups <- as.character(cellgroups)
  uid <- unique(cellgroups)
  naggr <- length(uid)
  ngene <- nrow(x)
  if( length(cellgroups) != ncol(x) ) stop("length(cellgroups) != ncol(x)")
  y <- array(0,dim=c(2,naggr,ngene),
    dimnames=list(c("u","v"),uid,rownames(x)))
  a <- data.frame(id=uid,N=0)
  rownames(a) <- a$id

  for(i in uid)
    {
    print(i)
    j <- grep(i,cellgroups,fixed=T)
    Ni <- length(j)
    a[i,"N"] <- Ni

    if(Ni == 1 )
      {
      y[1,i,] <- x[,j]
      y[2,i,] <- Inf
      }
    else
      {
      Sx <- Matrix::rowSums( x[,j] )
      Sxx <- Matrix::rowSums( x[,j]^2 )
      y[1,i,] <- Sx/Ni
      y[2,i,] <- (v0 + (Sxx-Sx*Sx/Ni)/(Ni-1))/Ni
      }
    }
  list(annot=a,y=y)
}
