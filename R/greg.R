greg <- function( 
  x, cell_sample, cell_type, v0 = 1/12, sep="\t", verbose = 1 )
{
  if(length(cell_sample) != length(cell_type))
    stop("length(cell_sample) != length(cell_type)")
  n_cells <- length(cell_sample)

  u_sample <- unique(as.character(cell_sample))
  u_ctype <- unique(as.character(cell_type))
  
  n_sample <- length(u_sample)
  n_ctype <- length(u_ctype)

  # cell
  N <- array(0,dim=c(n_sample,n_ctype),dimnames=list(u_sample,u_ctype))
  
  if( is.matrix(x) || is(x,"sparseMatrix" ))
    {
    n_gene <- nrow(x)
    n_cells <- ncol(x)
    if( n_cells != length(cell_sample))
      stop("length of cell labels doesn't match the number of cells in data")
    
    y <- array( c(0,Inf), dim=c(2,n_sample,n_gene,n_ctype),
                dimnames=list(c("u","v"),u_sample,rownames(x),u_ctype))

    for(k in u_ctype )
      {
      for(i in u_sample )
        {
        message(i," ",k)

        l <- intersect( 
          grep(i, cell_sample, fixed=T),
          grep(k, cell_type, fixed=T))
        N[i,k] <- Nik <- length(l)
        
        if( Nik == 1 )
          {
          y[1,i,,k] <- x[,l]
          }
        else
          {
          Sx <- rowSums( x[,l] )
          Sxx <- rowSums( x[,l]^2 )
          y[1,i,,k] <- Sx/Nik
          y[2,i,,k] <- (v0 + (Sxx-Sx*Sx/Nik)/(Nik-1))/Nik
          }
        } # k
      } # i
    }
  else if( is.character(x) && length(x)==1 ) # interpret as filename
    {
    f <- gzfile(x,open="rt")

    h <- strsplit( readLines(f,1), sep, fixed=TRUE )[[1]]
    nh <- length(h)
    if(verbose)
      message(sprintf("nh = %d n_cells = %d",nh,n_cells))
    if( nh != n_cells + 1 )
      stop("number of columns in header does not match cellgroups")

    chunk <- 1200
    n_gene <- 0

    y <- array(c(0,Inf),dim=c(2,n_sample, chunk, n_ctype),
                dimnames=list(c("u","v"),u_sample,1:chunk,u_ctype))

    cidx <- sapply( u_ctype,function(k)
             sapply( u_sample,function(i)
               intersect(grep(k,cell_type,fixed=T),
                 grep(i,cell_sample,fixed=T))))

    N <- sapply( u_ctype, function(k)
            sapply( u_sample, function(i) length(cidx[i,k][[1]])))

    repeat
      {
      gene <- scan(f,"",1,quiet=TRUE)
      if(length(gene)==0) break  # EOF
      
      counts <- scan(f,0.0,n_cells,quiet=TRUE)
      n_gene <- n_gene + 1
      dimnames(y)[[3]][n_gene] <- gene

      for(k in u_ctype )
        for(i in u_sample) 
          {
          yc <- counts[cidx[i,k][[1]]]
          y[1,i,n_gene,k] <- mean(yc)
          SE2 <- (stats::var(yc) + v0)/N[i,k]
          if(is.na(SE2)) SE2 <- Inf
          y[2,i,n_gene,k] <- SE2
          }

      if( verbose && ( n_gene %% 50 == 0) )
        message(sprintf("\rprocessed: %6d genes (alloc. chunk %7d)",n_gene,chunk),appendLF=FALSE)

      if( n_gene == chunk )
        {
        chunk <- 2*chunk
        y <- abind( y, array(c(0,Inf),dim=c(2,n_sample,chunk,n_ctype)),along=3)
        gc()
        }
      }
    if(verbose)
      message(sprintf("\rprocessed: %6d genes",n_gene),appendLF=TRUE)
    close(f)
    y <- y[,,1:n_gene,]
    gc()
    }
  
  structure( list(N=N,y=y), class="greg")
}

rm_ctype <- function(ac,cid=NULL)
{
  ac$N <- ac$N[,-cid]
  ac$y <- ac$y[,,,-cid]
  ac
}
