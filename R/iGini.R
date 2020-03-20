### R program for computing iGini decomposision
### Version 1.0, February 2020
### Tim Liao, University of Illinois

iGini <- function(x,g,w=rep(1,length(x)),core=1)
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
    cl <- parallel::makeCluster(core)
    doParallel::registerDoParallel(cl)
    i <- NULL

    res <- foreach (i=1:n, .combine='rbind') %dopar% {
      g.ikb=g.ikw <- 0
      # new code for vectorization
      g.ikb <- sum(w[i]*w[1:n][g[i]!=g[1:n]]*abs(x[i]-x[1:n][g[i]!=g[1:n]]))
      g.ikw <- sum(w[i]*w[1:n][g[i]==g[1:n]]*abs(x[i]-x[1:n][g[i]==g[1:n]]))
      c({g.ikb+g.ikw}/c,g.ikb/c,g.ikw/c)
    }
    parallel::stopCluster(cl)
    output[,1] <- res[,1]
    output[,2] <- res[,2]
    output[,3] <- res[,3]
    colnames(output) <- c("g.i","g.ikb","g.ikw")
    return(output)
}
