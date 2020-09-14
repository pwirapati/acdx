malm <- function(
  Y, X, w=NULL, offset=NULL, G=NULL,
  iter_max=40,epsilon=1e-4,y0=1e-5,
  bt_iter_max=4,bt_step=0.5,bt_tol=1,
  verbose=0  
)
{
  if(length(dim(Y)) < 3 ) stop("length(dim(Y)) < 3")
  n <- dim(Y)[2]
  m <- dim(Y)[3]
  if(!is.matrix(X)) X <- cbind(X)
  if(nrow(X)!= n) stop("nrow(x) != n")
  p <- ncol(X)
  if(is.null(w)) w <- rep(1,n)
  if(is.null(offset)) offset <- rep(0,n)
  if(is.null(G)) G <- rep("intercept",n)
  Gi <- as.integer(factor(G,levels=unique(G)))-1
  q <- length(unique(Gi))
  r <- .C("malm",
    dim=as.integer(c(n,m,p,q)),
    Y=as.double(Y + c(y0,0)),
    w=as.double(w),
    Xt=as.double(t(X)),
    Gi=as.integer(Gi),
    E=double(2*n*m),
    offset=as.double(offset),
    beta=double(p*m),
    phi=double(q*m),
    L=double(m),
    iopt=as.integer(c(verbose,iter_max,bt_iter_max)),
    dopt=as.double(c(epsilon,bt_tol,bt_step)),
    n_iter=integer(1*m),
    NAOK=TRUE,
    PACKAGE="acdx"
    )
  r$Y <- NULL
  dim(r$E) <- c(2,n,m)
  dimnames(r$E) <- dimnames(Y)
  dim(r$beta) <- c(p,m)
  dim(r$phi) <- c(q,m)
  rownames(r$phi) <- unique(G)
  colnames(r$phi) <- colnames(r$beta) <- dimnames(Y)[[3]]
  r
}
