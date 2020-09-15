greg <- function( 
  x,                # matrix object (dgCMatrix or dense), or a file name
  cell_sample,
  cell_type, 
  cell_names = NULL, # cell names, to be matched with names in x
  v0 = 1/12, 
  sep="\t", verbose = 1 )
{
  if(length(cell_sample) != length(cell_type))
    stop("length(cell_sample) != length(cell_type)")
  n_cells <- length(cell_sample)

  if(!is.null(cell_names) && length(cell_names) != n_cells )
    stop("length(cell_names) != length(cell_sample)")

  u_sample <- unique(as.character(cell_sample))
  u_ctype <- unique(as.character(cell_type))
  
  n_sample <- length(u_sample)
  n_ctype <- length(u_ctype)

  # cell
  N <- array(0,dim=c(n_sample,n_ctype),dimnames=list(u_sample,u_ctype))
  
  if( is.matrix(x) || is(x,"sparseMatrix" ))
    {
    # cell_names && no colnames: error
    # cell_names && colnames: intersect
    # no cell_nmaes && no colnames: assume integer, length must match
    # no cell_names && colnames: ditto

    if( is.null(cell_names))
      {
      if( ncol(x) != n_cells )
        stop("number of cells annotated doesn't match number of cells in x")
      la2lx <- 1:n_cells
      }
    else 
      {
      la2lx <- pmatch(cell_names, colnames(x))
      n_nomatch <- sum(is.na(la2lx))
      if(n_nomatch != 0)
        warning(n_nomatch," cells are not matched.")
      }

    n_gene <- nrow(x)

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

        l <- na.omit(la2lx[l])
        
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
  else if( is.character(x) && length(x)==1 ) # interpret x as a filename
    {
    f <- gzfile(x,open="rt")

    h <- strsplit( readLines(f,1), sep, fixed=TRUE )[[1]]
    n_headers <- length(h)
    if(verbose)
      message(sprintf("n_headers = %d, n_cells = %d",n_headers,n_cells))

    if( is.null(cell_names))
      {
      if( n_headers != n_cells + 1 )
        stop("number of cells annotated doesn't match number of cells in x")
      la2lx <- 1:n_cells
      }
    else
      {
      la2lx <- pmatch( cell_names, h[-1])
      n_nomatch <- sum(is.na(la2lx))
      if(n_nomatch != 0)
        warning(n_nomatch," cells are not matched.")
      }

    chunk <- 1200
    n_gene <- 0

    y <- array(c(0,Inf),dim=c(2,n_sample, chunk, n_ctype),
                dimnames=list(c("u","v"),u_sample,1:chunk,u_ctype))

    cidx <- sapply( u_ctype,function(k)
             sapply( u_sample,function(i)
               na.omit(la2lx[intersect(grep(k,cell_type,fixed=T),
                 grep(i,cell_sample,fixed=T))]),simplify=FALSE))

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
