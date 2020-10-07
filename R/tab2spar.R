tab2spar <- function(con, transpose=FALSE)
{
  if(!inherits(con,"connection"))
    {
    if( is.character(con) && file.exists(con) )
      con <- gzfile(con,"rb")
    else
      stop("first arg not a connection nor a file")
    }

  r <- .Call("tab2spar", con, as.logical(transpose))
  close(con)
  
  sparseMatrix(
    i=r[[3]],p=r[[5]],x=r[[4]],
    dims=c(length(r[[1]]),length(r[[2]])),
    dimnames=list(r[[1]],r[[2]]),
    index1=FALSE
  )
}
