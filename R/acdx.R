# iteration control options
acdx.ctrl <- list(
    iter_max11 = 30, epsilon11 = 1e-4,
    iter_max = 30, epsilon = 1e-4,
    bt_iter_max = 4,
    bt_tol = 1e-4,
    bt_step = 0.5)

acdx <- function(
  ac,           # aggregated cell object
  X,            # design matrix for mean parameters
  Gi=NULL,      # group intercepts for dispersion parameters

  n_boot=1,     # number of bootstrap resampling
  id_boot=NULL, # block bootstrap id (default: sample id)
  seed_boot=0,  # random number seed
  
  y0=NULL,      # 'pseudo-count'
  verbose=0
  )
{
  n_sample <- nrow(ac$N)
  n_ctype <- ncol(ac$N)
  n_gene <- dim(ac$y)[3]
  n_coef <- ncol(X)

  if(is.null(Gi)) Gi <- rep(0,n_sample)
  Gi_names <- unique(Gi)
  Gi <- as.integer(factor(Gi,levels=Gi_names))-1
  n_disp <- length(unique(Gi))

  if(is.null(id_boot)) id_boot <- 0:(n_sample-1)
  id_boot <- as.integer(factor(id_boot,levels=unique(id_boot)))-1
  n_bid <- length(unique(id_boot))
 
  set.seed(seed_boot)
  B <- sapply(1:n_boot,function(u)sample.int(n_bid,n_bid,replace=TRUE)-1)
  B[,1] <- 0:(n_bid-1)

  if(is.null(y0)) y0 <- 0.5*min(ac$y[1,,,][ ac$y[1,,,] > 0],na.rm=T)

  r <- .C("malm_acdx",
    dim=as.integer(c(n_sample,n_gene,n_coef,n_disp,n_ctype,n_bid,n_boot)),
    Y=as.double(ac$y),
    N=as.double(ac$N),
    Xt=as.double(t(X)),
    Gi=as.integer(Gi),
    id_boot=as.integer(id_boot),
    B=as.integer(B),
    gamma=double(n_gene*n_ctype*n_boot),
    psi=double(n_gene*n_ctype*n_boot),
    alpha=double(n_sample*n_ctype*n_boot),
    beta=double(n_coef*n_gene*n_ctype*n_boot),
    phi=double(n_disp*n_gene*n_ctype*n_boot),
    y0=as.double(y0),
    iopt11=as.integer(c(verbose,acdx.ctrl$iter_max11)),
    dopt11=as.double(c(acdx.ctrl$epsilon11)),
    iopt=as.integer(c(verbose,acdx.ctrl$iter_max,acdx.ctrl$bt_iter_max)),
    dopt=as.double(c(acdx.ctrl$epsilon,acdx.ctrl$bt_tol,acdx.ctrl$bt_step)),
    NAOK=TRUE,
    PACKAGE="acdx"
    )
  r$Y <- NULL
  r$N <- NULL
  sample_names <- rownames(ac$N)
  ctype_names <- colnames(ac$N)
  gene_names <- dimnames(ac$y)[[3]]

  dim(r$gamma) <- c(n_gene,n_ctype,n_boot)
  dim(r$psi) <- c(n_gene,n_ctype,n_boot)
  dimnames(r$gamma) <- dimnames(r$psi) <- list(gene_names,ctype_names,NULL)
  dim(r$alpha) <- c(n_sample,n_ctype,n_boot)
  dimnames(r$alpha) <- list(sample_names,ctype_names,NULL)
  dim(r$beta) <- c(n_coef,n_gene,n_ctype,n_boot)
  dimnames(r$beta) <- list(colnames(X),gene_names,ctype_names,NULL)
  dim(r$phi) <- c(n_disp,n_gene,n_ctype,n_boot)
  dimnames(r$phi) <- list( Gi_names,gene_names,ctype_names,NULL)
  r$seed_boot <- seed_boot
  dim(r$B) <- c(n_bid,n_boot)
  structure(r,class="acdx_fit")
}
