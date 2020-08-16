aggr <- function( x, cellgroups, v0 = 1/12,
  sep="\t",
  verbose=1 )
{
  cellgroups <- as.character(cellgroups)
  ncells <- length(cellgroups)
  uid <- na.omit(unique(cellgroups))
  n_ac <- length(uid)
  a <- data.frame(id=uid,N=0)
  rownames(a) <- a$id

  if(is.matrix(x) || is(x,"sparseMatrix"))
    {
    n_gene <- nrow(x)
    if( ncells != ncol(x))
      stop("length(cellgroups) != ncol(x)")
    y <- array(0,dim=c(2,n_ac,n_gene),
      dimnames=list(c("u","v"),uid,rownames(x)))

    for(i in uid) # loop over aggregates defined by unique id
      {
      if(verbose) message(i)
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
  else if( is.character(x) && length(x)==1 ) # interprete as filename
    {
    f <- gzfile(x,open="rt")

    h <- strsplit( readLines(f,1), sep, fixed=TRUE )[[1]]
    nh <- length(h)
    if(verbose)
      message(sprintf("nh = %d ncells = %d",nh,ncells))
    if( nh != ncells + 1 )
      stop("number of columns in header does not match cellgroups")

    chunk <- 1200
    n_gene <- 0
    y <- array(0,dim=c(2,n_ac, chunk),
                dimnames=list(c("u","v"),uid,1:chunk))

    j <- sapply( uid, function(i) grep(i,cellgroups,fixed=T))
    N <- sapply(j,length)
    a$N <- N

    repeat
      {
      gene <- scan(f,"",1,quiet=TRUE)
      if(length(gene)==0) break
      
      counts <- scan(f,0.0,ncells,quiet=TRUE)
      n_gene <- n_gene + 1
      dimnames(y)[[3]][n_gene] <- gene
      for(i in 1:length(uid)) 
        {
        yc <- counts[j[[i]]]
        y[1,i,n_gene] <- mean(yc)
        SE2 <- (var(yc) + v0)/N[i]
        if(is.na(SE2)) SE2 <- Inf
        y[2,i,n_gene] <- SE2
        }
      if( verbose && ( n_gene %% 50 == 0) )
        message(sprintf("\rprocessed: %6d genes (alloc. chunk %7d)",n_gene,chunk),appendLF=FALSE)

      if( n_gene == chunk )
        {
        chunk <- 2*chunk
        y <- abind( y, array(0,dim=c(2,n_ac,chunk)),along=3)
        gc()
        }
      }
    if(verbose)
      message(sprintf("\rprocessed: %6d genes",n_gene),appendLF=TRUE)
    close(f)
    y <- y[,,1:n_gene]
    gc()
    }
  list(annot=a,y=y)
}
