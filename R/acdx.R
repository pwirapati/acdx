#
# bootstrap resampling by modifying the weights
#
bootweight <-
function(x) {
  x <- as.character(x)
  ux <- unique(x)
  b <- table(sample(ux,replace=T))
  w <- rep(0,length(x))
  names(w) <- x
  z <- setdiff(ux,names(b))
  wz <- rep(0,length(z))
  names(wz) <- z
  b <- c(b,wz)
  b[x]
}


# turn data frame into integer-coded factors for C routine
cfactor <- function(x)
{
  if(!is.data.frame(x)) x <- as.data.frame(x)

  x <- lapply(x,function(u) 
    if(is.factor(u)) u
    else factor(u,levels=unique(u)))
  level_names <-  lapply(x,levels)
  level_k <- lapply(level_names,length)
  x <- sapply(x,function(u) as.integer(u)-1)
  attr(x,"level_names") <- level_names
  attr(x,"level_k") <- unlist(level_k)
  return(x)
}

acdx <- function(
  Y,
  Zf,
  X1,
  G1,
  Xf=NULL,
  Xc=NULL,
  w=NULL,
  iter_tol=1e-5,
  bt_tol=-0.1,
  inv_tol=1e-12,
  iter_max=30,
  bt_max=20,
  wssize_min=4096/8,
  y_eps=NULL,
  verbose=1,
  init=TRUE,
  keepY=FALSE,
  keepE=FALSE,
  bootid=NULL,
  nboot=99,
  bootpar="beta"
)
{
  if(length(dim(Y))!=3) stop("Y is not a 3d array ")
  if(dim(Y)[1] != 2) stop("first dim of Y is not 2")
  n_ac <- dim(Y)[2]
  n_gene <- dim(Y)[3]

  if(is.null(y_eps))
    y_eps <- 0.5*min( Y[1,,][Y[1,,] > 0] )
  if(!is.finite(y_eps))
    stop("y_eps not finite")
  Y[1,,] <- Y[1,,] + y_eps

  if(is.null(w)) w <- rep(1,n_ac)
  if( length(w) != n_ac ) stop("length(w) != number of aggregates")
  w <- ifelse(w < 0, 0, w)

  if(is.null(Zf)) Zf <- rep("raw",n_ac)
  if(length(Zf) != n_ac) stop("length(Zf) != number of aggregates")
  if(any(is.na(Zf))) stop("NA in Zf")
  Zf <- cfactor(Zf)
  k_Zf <- attr(Zf,"level_k")[1]

  if(is.null(X1)) X1 <- rep("Intercept",n_ac)
  if(length(X1) != n_ac) stop("length(X1) != number of aggregates")
  if(any(is.na(X1))) stop("NA in X1")
  X1 <- cfactor(X1)
  k_X1 <- attr(X1,"level_k")[1]

  if(!is.null(Xf))
    {
    if(any(is.na(Xf))) stop("NA in Xf")
    Xf <- cfactor(Xf)
    if(nrow(Xf) != n_ac) stop("nrow(Xf) != number of aggregates")

    n_Xf <- ncol(Xf)
    k_Xf <- attr(Xf,"level_k")
    }
  else
    {
    n_Xf <- k_Xf <- 0 
    Xf <- integer(0)
    }

  if(!is.null(Xc))
    {
    if(!is.matrix(Xc)) Xc <- cbind(Xc=Xc)
    if( nrow(Xc) != n_ac) stop("nrow(Xc) != number of aggregates")
    if( !is.numeric(Xc)) stop("Xc is not numeric")
    if(any(is.na(Xc))) stop("NA in Xc")
    n_Xc <- ncol(Xc)
    }
  else
    {
    n_Xc <- 0
    Xc <- double(0)
    }

  if(is.null(G1)) rep("Intercept",n_ac)
  if(length(G1) != n_ac) stop("length(G1) != number of aggregates")
  if(any(is.na(G1))) stop("NA in G1")
  G1 <- cfactor(G1)
  k_G1 <- attr(G1,"level_k")[1]

  r <- .C("malm_acdx",
    dims=as.integer( c(n_ac, n_gene,
        1,1,n_Xf,n_Xc,1) ),
    w=as.double(w),
    Y=as.double(Y),
    X1=as.integer(X1),
    Xf=as.integer(Xf),
    Xc=as.double(Xc),
    Zf=as.integer(Zf),
    G1=as.integer(G1),

    beta=double(n_gene * (k_X1 + sum(k_Xf) + n_Xc)),
    alpha=double(k_Zf),
    phi=double(n_gene * k_G1),

    E=double(2*n_ac*n_gene),

    tol=as.double(c(iter_tol,bt_tol,inv_tol)),
    limits=as.integer(c(iter_max,bt_max,wssize_min)),
    options=as.integer(c(init,verbose)),
    status=integer(3),

    NAOK=T,
    PACKAGE="acdx")

  names(r$dims) <- c("n_ac","n_gene","n_Zf","n_X1","n_Xf","n_Xc","n_G1")
  names(r$alpha) <- attr(Zf,"level_names")[[1]]

  dim(r$phi) <- c(k_G1,n_gene) 
  rownames(r$phi) <- attr(G1,"level_names")[[1]]
  colnames(r$phi) <- dimnames(Y)[[3]]

  dim(r$beta) <- c(k_X1+sum(k_Xf),n_gene)
  colnames(r$beta) <- dimnames(Y)[[3]]
  rownames(r$beta) <- unname(c(
                        attr(X1,"level_names")[[1]],
                        unlist(attr(Xf,"level_names")),
                        colnames(Xc)))

  names(r$status) <- c("status","n_iter","n_backtrack")
  names(r$options) <- c("init","verbose")
  names(r$limits) <- c("iter_max","bt_max","wssize_min")
  names(r$tol) <- c("iter_tol","bt_tol","inv_tol")

  if(keepY) 
    attributes(r$Y) <- attributes(Y)
  else
    r$Y <- NULL
  if(keepE) 
    attributes(r$E) <- attributes(Y)
  else
    r$E <- NULL

  if( !is.null(bootid) && nboot > 0 )
    {
    if("alpha" %in% bootpar)
      r$Alpha <- array(0,dim=c(length(r$alpha),nboot))
    if("beta" %in% bootpar)
      r$Beta <- array(0,dim=c(nrow(r$beta),ncol(r$beta),nboot))
    if("phi" %in% bootpar)
      r$Phi <- array(0,dim=c(nrow(r$phi),ncol(r$phi),nboot))

    for(b in 1:nboot)
      {
      set.seed(b)
      if(verbose > 0 ) cat(b,":")
      rb <- .C("malm_acdx",
        dims=as.integer(r$dims),
        w=as.double( w * bootweight(bootid)),
        Y=as.double(Y),
        X1=as.integer(X1),
        Xf=as.integer(Xf),
        Xc=as.double(Xc),
        Zf=as.integer(Zf),
        G1=as.integer(G1),
        beta=double(n_gene * (k_X1 + sum(k_Xf) + n_Xc)),
        alpha=double(k_Zf),
        phi=double(n_gene *k_G1),
        E=double(2*n_ac*n_gene),
        tol=as.double(r$tol),
        limits=as.integer(r$limits),
        options=as.integer(r$options),
        status=integer(3),
        NAOK=T,
        PACKAGE="acdx"
        )
      if(!is.null(r$Alpha)) r$Alpha[,b] <- rb$alpha
      if(!is.null(r$Beta)) r$Beta[,,b] <- rb$beta
      if(!is.null(r$Phi)) r$Phi[,,b] <- rb$phi
      }
    if(!is.null(r$Alpha))
      dimnames(r$Alpha) <- list(names(r$alpha),1:nboot)
    if(!is.null(r$Beta))
      dimnames(r$Beta) <- c(dimnames(r$beta),list(1:nboot))
    if(!is.null(r$Phi))
      dimnames(r$Phi) <- c(dimnames(r$phi),list(1:nboot))
    r$bootid <- bootid
    r$nboot <- nboot
    r$bootpar <- bootpar
    }
  
  return(r)
}
