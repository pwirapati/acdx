questf <- function(fit, ifun)
{
  n_boot <- ncol(fit$B)

  cat("calculating interest functional...")
  theta_boot <- sapply( 1:n_boot,function(b)
                   ifun(beta=fit$beta[,,,b],phi=fit$phi[,,,b]))
  cat("DONE.\n")

  bremt(theta_boot)
}

questmm <- function(...)
{
  quest(method="minmax",...)
}

quest <- function(
  param,
  showcomb = F,
  default="L",
  H=NULL,
  N=NULL,
  L=NULL,
  sep=",",
  na.rm=FALSE,
  include=NULL,
  exclude=NULL,
  method="minmax",
  multest_func=bremt,
  ...   # options to multest_func
  )
{
  if(class(param)=="acdx_fit")
    param <- param$beta

  coef_names <- dimnames(param)[[1]]
  n_coef <- length(coef_names)
  ctype_names <- dimnames(param)[[3]]
  n_ctype <- length(ctype_names)
  ptab <- expand.grid(coef=coef_names,ctype=ctype_names)
  rownames(ptab) <- paste(ptab[,1],ptab[,2],sep=sep)

  if(showcomb) { return(cbind(rownames(ptab))) }

  if(!grepl("^[HNL]$",default))
    stop("default= should be either H, N or L")
  ptab[["set"]] <- default
  if(!is.null(H))
    {
    if(default=="H") stop("H is the default")
    for(pat in H)
      {
      i <- grep(pat,rownames(ptab))
      if(length(i)==0) warning("H: ",pat," not matched");
      ptab[i,"set"] <- "H"
      }
    }
  if(!is.null(N))
    {
    if(default=="N") stop("N is the default")
    for(pat in N)
      {
      i <- grep(pat,rownames(ptab))
      if(length(i)==0) warning("N: ",pat," not matched");
      ptab[i,"set"] <- "N"
      }
    }
  if(!is.null(L))
    {
    if(default=="L") stop("L is the default")
    for(pat in L)
      {
      i <- grep(pat,rownames(ptab))
      if(length(i)==0) warning("L: ",pat," not matched");
      ptab[i,"set"] <- "L"
      }
    }

  idH <- grep("H",ptab$set)
  idHcoef <- match( ptab[idH,1], coef_names)
  idHctype <- match( ptab[idH,2], ctype_names)
  idL <- grep("L",ptab$set)
  idLcoef <- match( ptab[idL,1], coef_names)
  idLctype <- match( ptab[idL,2], ctype_names)

  g <- dimnames(param)[[2]]
  if(!is.null(include))
    g <- intersect(g,include)
  if(!is.null(exclude))
    g <- setdiff(g,exclude)

  if(length(g)==0)
    stop("No genes are pre-selected")

  param <- param[,g,,]

  if(method=="ave") cfunc <- "avediff"
  else cfunc <- "mmdiff"

  r <- .C(cfunc,
    dims=as.integer(dim(param)),
    param=as.double(param),
    nH=as.integer(length(idH)),
    Hctype=as.integer(idHctype-1),
    Hcoef=as.integer(idHcoef-1),
    nL=as.integer(length(idL)),
    Lctype=as.integer(idLctype-1),
    Lcoef=as.integer(idLcoef-1),
    theta=double(dim(param)[2]*dim(param)[4]),
    na_rm=as.integer(na.rm),
    NAOK=TRUE,
    PACKAGE="acdx"
  )
  r$param <- NULL
  dim(r$theta) <- c(dim(param)[2],dim(param)[4])
  rownames(r$theta) <- dimnames(param)[[2]]
  multest_func(r$theta,...)
}
