#
# utilities for multidimensional array
#

# collapse two dimensions of an array
# - the new combined dimensions are now at position `pos`
# - the length is the product of the old one
# - the names is the concatenated labels separated by "|"
collapse2dim <- function(
  a, # array with dimensions >= 2
  cid,  # an ordered pair of dimension indices to be collapsed
  pos = 1,
  sep = "|"
  )
{
  if(!is.array(a) || length(dim(a)) < 2)
    stop("a is not array with 2 or more dimensions")

  nda <- length(dim(a))
  if( !is.numeric(cid) || length(cid) != 2 || max(cid) > nda
    || min(cid) < 1 )
    stop("wrong format of collapsed indices")
  
  dperm <- append( setdiff(1:nda,cid), cid, after=pos-1)
  a <- aperm(a, dperm)
  dnames <- dimnames(a)
  da <- dim(a)
  da[pos] <- da[pos]*da[pos+1]
  dim(a) <- da[-(pos+1)]
  dnames[[pos]] <- c(outer(dnames[[pos]],dnames[[pos+1]],paste,sep=sep))
  dnames <- dnames[-(pos+1)]
  dimnames(a) <- dnames
  a
}
