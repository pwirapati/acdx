greg <- function( 
  x,                # matrix object (dgCMatrix or dense), or a file name
  cell_sample,
  cell_type, 
  cell_names = NULL, # cell names, to be matched with names in x
  s2_0 = 0, 
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
          y[2,i,,k] <- (s2_0 + (Sxx-Sx*Sx/Nik)/(Nik-1))/Nik
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
          SE2 <- (stats::var(yc) + s2_0)/N[i,k]
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
    }
  
  gene_mean <- apply(y[1,,,],2,mean,na.rm=T)
  gene_cv <- apply(y[1,,,],2,sd,na.rm=T)/gene_mean

  structure(
    list(
      N=N,
      y=y,
      source=src,
      gene_mean = gene_mean,
      gene_cv = gene_cv
      ),
   class="ac")
}

print.ac <- function(ac)
{
  cat("Aggregated cells from", ac$source,"\n")
  cat("number of cells = ", sum(ac$N),"\n",sep="")
  cat("number of samples = ",nrow(ac$N),"\n",sep="")
  cat("number of cell types = ",ncol(ac$N),"\n",sep="")
  cat("number of aggregates = ",length(ac$N),", (non-empty = ",sum(ac$N > 1),")\n",sep="")
  cat("number of genes = ",dim(ac$y)[3]," = ", 
      sum(ac$gene_mean > 0)," (some positive) + ",
      sum(ac$gene_mean==0)," (all zeroes))\n",sep="")
}

plot_ac_genesum <- function(ac,cex=.25,pch=20,n_names=20,...)
{
  rcv <- range(ac$gene_cv[is.finite(ac$gene_cv)])
  wcv <- (rcv[2]-rcv[1])/n_names
  plot( log10(ac$gene_mean), ac$gene_cv,
   xlab="gene mean", ylab="gene CV", axes=FALSE,
   cex=cex,pch=pch, ylim=c(0,rcv[2]),
    ,...)
  axis(side=2)
  axis(side=1,at=-5:5,labels=10^{-5:5})
  
  gbin <- .bincode(ac$gene_cv,seq(rcv[1]-wcv,rcv[2]+wcv, wcv))
  g <- tapply( ac$gene_mean, gbin, function(u) names(which.max(u)))
  text( log10(ac$gene_mean[g]),ac$gene_cv[g],labels=g, adj=c(-.2,.5),xpd=TRUE)
}

#
#plot cell counts, cell_type % / sample (in logit), box plot
# of average counts, boxplot of sums (pseudo-bulk)
#
plot_ac_aggrsum <- function( ac, o=NULL, s=1e-3, pad=2, ... )
{
  if(is.null(o)) o <- 1:nrow(ac$N)
  xlim <- c(1-pad,length(ac$N)+pad)

  on.exit(layout(1))
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar),add=TRUE,after=FALSE)
  oldmar <- par('mar')

  layout(cbind(c(1,2,3)),widths=c(1),heights=c(1,1,1.5))
  par(cex=1)

  par(mar=oldmar+c(-5,0,-2,0))
  plot( log(.5+c(ac$N[o,])), xlim=xlim, xaxs="i", ylab="cell counts",
    axes=FALSE, frame=TRUE, xlab=NA )
  axis(side=2,at=log(.5+c(0,1,10,100,1000)),labels=c(0,1,10,100,1000),las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)

  Ns <- rep(rowSums(ac$N[o,]),ncol(ac$N))
  plot( log((.01+c(ac$N[o,]))/(.01+Ns-c(ac$N[o,]))), xlim=xlim, xaxs="i", ylab="cell counts",
    axes=FALSE, frame=TRUE, xlab=NA )
  axis(side=2,at=log((+c(.001,.01,.1,.25,.5,.75,.9))/(c(.999,.99,.9,.75,.5,.25,.1))),
    labels=c("0.1%","1%","10%","25%","50%","75%","90%"),las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)


  par(mar=oldmar+c(4,0,-2,0))
  plot.new()
  plot.window(xlim=xlim,xaxs="i",ylim=log(s+c(0,1000)))
  boxplot.matrix( log(s+apply(ac$y[1,o,KUL3$gene_mean > 0,],2,c)),
    border=gray(.5), 
    use.cols=FALSE,outline=F,add=TRUE,axes=FALSE,frame=TRUE
    )
  title(ylab="mean raw counts")
  ti <- c(0,unlist(sapply(-3:4,function(p) c(1)*10^p)))
  axis(side=2,at=log(s+ti),labels=ti,las=2)
  abline(v=(1:(dim(ac$y)[4]-1))*dim(ac$y)[2] + .5)

  L <- list( rep( rownames(ac$N)[o],ncol(ac$N)), 
    c( sapply(colnames(ac$N),function(u) rep(u,nrow(ac$N)))))
  dtext( L, side=1,line=0, bspace=1,las=c(2,0))

}

plot.ac <- function( ac, genes=NULL, ... )
{
  if(is.null(genes))
    plot_ac_genesum( ac, ... )
  else if(is.na(genes))
    plot_ac_aggrsum( ac, ... )
}

rm_ctype <- function(ac,cid=NULL)
{
  ac$N <- ac$N[,-cid]
  ac$y <- ac$y[,,,-cid]
  ac
}
