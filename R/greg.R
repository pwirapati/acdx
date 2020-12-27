greg <- function( 
  x,                # matrix object (dgCMatrix or dense), or a file name
  cell_sample,
  cell_type, 
  cell_names = NULL, # cell names, to be matched with names in x
  u_0 = NULL,
  s2_0 = 1/12, 
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
    src <- paste("Matrix",deparse(substitute(x)))
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
          y[2,i,,k] <- ((Sxx-Sx*Sx/Nik)/(Nik-1))/Nik
          }
        } # k
      } # i
    }
  else if( is.character(x) && length(x)==1 ) # interpret x as a filename
    {
    src <- paste("file",x)
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
          SE2 <- (stats::var(yc))/N[i,k]
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
  else if( class(x) == "ac" || class(x) == "greg") # to convert old format
    {
    src <- paste(class(x), deparse(substitute(x)))
    y <- x$y
    N <- x$N
    n_sample <- dim(y)[2]
    n_gene <- dim(y)[3]
    n_ctype <- dim(y)[4]
    }

  if( is.null(u_0)) u_0 <- 0.5/mean(N)

  gene_mean0 <- apply(y[1,,,],2,mean,na.rm=T)
  gene_sd0 <- apply(y[1,,,],2,sd,na.rm=T)

  for(k in 1:ncol(N))
    y[,,,k] <- y[,,,k] + c( rbind(u_0,s2_0/N[,k]))

  w <- apply(1/y[2,,,],2,sum,na.rm=T)
  gene_wmean_u <- apply( y[1,,,]/y[2,,,],2,sum,na.rm=T)/w
  gene_wmean_v <- apply( y[2,,,]/y[2,,,],2,sum,na.rm=T)/w

  # aggregate and gene scale factors
  gene_scale <- array(dim=c(n_gene,n_ctype))
  gene_phi0 <- array(dim=c(n_gene,n_ctype))
  aggr_scale <- array(dim=c(n_sample,n_ctype))

  for(k in 1:n_ctype)
    {
    fit <- malm11( y[,,,k], w=ifelse(N[,k] > 1,1,0))
    gene_scale[,k] <- exp(fit$beta0)
    gene_phi0[,k] <- exp(fit$phi0)
    aggr_scale[,k] <- exp(fit$alpha)
    }

  structure(
    list(
      N=N,
      y=y,
      source=src,
      u_0 = u_0,
      s2_0 = s2_0,
      gene_mean0 = gene_mean0,
      gene_sd0 = gene_sd0,
      gene_wmean_u = gene_wmean_u,
      gene_wmean_v = gene_wmean_v,
      gene_scale = gene_scale,
      gene_phi0 = gene_phi0,
      aggr_scale = aggr_scale
      ),
   class="ac")
}

