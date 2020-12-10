# position of separators
sep <- function(x)
{
  y <- cumsum(rle(as.character(x))$lengths)
  y[-length(y)] + 0.5
}
