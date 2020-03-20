### R program for computing Thei's L individual decomposision
### Version 1.0, February 2020
### Tim Liao, University of Illinois

iTheilL <- function(x,g,w=rep(1,length(x)))
{
  ng <- length(unique(g))
  n <- length(x)
  N1 <- sum(w)
  N0 <- sum(w[x>0])
  sx <- sum(w*as.numeric(x))
  c <- 2*N1*sx
  f1 <- w/N1              # weighting factor
  f0 <- w[x>0]/N0
  xk <- aggregate(w*x,by=list(G=g),FUN="sum")
  yk <- xk[["x"]]/sx

  if (ng<=1)
    stop("data must have 2 or more groups to use this function")

    output <- matrix(NA,nrow=n,ncol=3)
    Tl.i=Tl.ib=Tl.iw=nk0 <- NULL
    # code for vectorization
    for (j in 1:ng) {
      nk0 <- sum(w[g==j])/N0
      ni0 <- w[g==j]/sum(w[g==j])
      Tl.i <- cbind(Tl.i,t(na.omit({w[g==j]/N0}*log({w[g==j]/N0}/{{w[g==j]*x[g==j]}/sx}))))
      Tl.ib <- cbind(Tl.ib,t(na.omit(ni0*nk0*log(nk0/yk[j]))))
      Tl.iw <- cbind(Tl.iw,t(na.omit(ni0*nk0*log(ni0/{{w[g==j]*x[g==j]}/xk$x[j]}))))
    }
    output[,1] <- Tl.i
    output[,2] <- Tl.ib
    output[,3] <- Tl.iw
    colnames(output) <- c("Tl.i","Tl.ib","Tl.iw")
    return(output)
}
