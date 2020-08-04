aggr <- function( x, cellgroups, v0 = 1/12 )
{
  #require(Matrix)

  cellgroups <- as.character(cellgroups)
  uid <- na.omit(unique(cellgroups))
  n_ac <- length(uid)
  a <- data.frame(id=uid,N=0)
  rownames(a) <- a$id

  if(is.matrix(x) || is(x,"sparseMatrix"))
    {
    n_gene <- nrow(x)
    if(length(cellgroups) != ncol(x))
      stop("length(cellgroups) != ncol(x)")
    y <- array(0,dim=c(2,n_ac,n_gene),
      dimnames=list(c("u","v"),uid,rownames(x)))

    for(i in uid) # loop over aggregates defined by unique id
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
        Sx <- rowSums( x[,j] )
        Sxx <- rowSums( x[,j]^2 )
        y[1,i,] <- Sx/Ni
        y[2,i,] <- (v0 + (Sxx-Sx*Sx/Ni)/(Ni-1))/Ni
        }
      }
    }
  list(annot=a,y=y)
}
