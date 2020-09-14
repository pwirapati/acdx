top <- function(
  sigsum,
  sort.by="T",
  aux=NULL,
  n=10,
  digits=4,
  decreasing=T,
  ...    # other options for order
  )
{
  sigsum <- data.frame(sigsum)
  if(!is.null(aux)) sigsum <- cbind(sigsum,aux)
  
  revdir <- ifelse(grepl("^-",sort.by),T,F)
  sort.by <- sub("^-","",sort.by)

  notfound <- setdiff(sort.by,colnames(sigsum))
  if(length(notfound) > 0)
    stop("sort key(s) not found: ",notfound)

  decreasing <- xor( rep(decreasing,length(sort.by)),revdir)

  o <- do.call("order",
    c(sigsum[,sort.by,drop=FALSE],list(decreasing=decreasing,...)))
  if(!is.numeric(n) || n < 1 || n > nrow(sigsum))
    n <- nrow(sigsum)
  round(sigsum[o[1:n],],digits)
}
