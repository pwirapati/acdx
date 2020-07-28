quest <- function(fit, ifun, siglevel=0.05, sort.by=NULL)
{
  if( is.null(fit$nboot) || fit$nboot < 2 ) stop("no bootstrap data")

  cat("calculating interest functional...")
  theta_boot <- sapply( 1:fit$nboot,function(b)
                   ifun(beta=fit$Beta[,,b],phi=fit$Phi[,,b],alpha=fit$Alpha[,b]))
  cat("DONE.\n")

  theta <-  ifun(beta=fit$beta,phi=fit$phi,alpha=fit$alpha)
  se <- apply(theta_boot,1,sd)
  
  cat("calculating bootstrap distribution...")
  maxT <- apply( scale(t(theta_boot),center=T,scale=T), 1, sort, decr=T )
  tnull=apply(maxT,1,quantile,siglevel)
  cat("DONE.\n")

  top=data.frame(theta=theta,se=se,t=theta/se,fc=exp(theta),
    row.names=names(theta))
  top <- top[order(top$t,decreasing=T),]
  top$FWER <- sapply(top$t,function(tt) (1+sum(tt < maxT[1,]))/(1+fit$nboot))
  top$FDP <- sapply( 1:nrow(top),function(r) (sum(tnull > top$t[r])+1)/(r+1) )
  for(r in (nrow(top)-1):1) top$FDP[r] <- min(top$FDP[r],top$FDP[r+1])
  top$FDP <- round(top$FDP,5)

  if(!is.null(sort.by) && sort.by=="fc")
    top <- top[order(top$fc,decreasing=T),]

  list(
    ifun=ifun,
    top=top,
    tnull=tnull,
    maxT.1=maxT[1,]
  )

}

