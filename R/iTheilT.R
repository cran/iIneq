### R program for computing Theil'sT individual decompositon
### Version 1.0.2, January 2021
### Tim Liao, University of Illinois

iTheilT <- function(x,g,w=rep(1,length(x)))
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
    Tt.i=Tt.ib=Tt.iw <- NULL
    # code for vectorization
    for (j in 1:ng) {
      nk1 <- sum(w[g==j])/N1
      yi <- {w[g==j]*x[g==j]}/xk[["x"]][j]
      Tt.i <- cbind(Tt.i,t(na.omit({{w[g==j]*x[g==j]}/sx}*log({{w[g==j]*x[g==j]}/sx}/{w[g==j]/N1}))))
      Tt.ib <- cbind(Tt.ib,t(na.omit(yi*yk[j]*log(yk[j]/nk1))))
      Tt.iw <- cbind(Tt.iw,t(na.omit(yi*yk[j]*log({{w[g==j]*x[g==j]}/xk$x[j]}/{w[g==j]/sum(w[g==j])}))))
    }
    output[,1] <- Tt.i
    output[,2] <- Tt.ib
    output[,3] <- Tt.iw
    colnames(output) <- c("Tt.i","Tt.b","Tt.iw")
    return(output)
}

