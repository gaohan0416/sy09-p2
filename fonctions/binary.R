library(Rlab)
library(MASS)

binaryNBCfit <- function(Xapp, zapp)
{
  #nb ligne n
  n <- nrow(Xapp)
  #nb col p
  p <- ncol(Xapp)
  #nb classes
  g <- max(unique(zapp))
  
  param <- NULL
  pkj <- array(0, c(g,p))
  pik <- rep(0, g)
  
  for (k in 1:g)
  {
    #indices ou z = k
    indk <- which(zapp==k)
    #individus de la classe k
    Xk <- Xapp[indk,]
    #nb d'individus dans la classe k
    nk <- length(indk)
    
    #moyenne des individus de la classe k
    pkj[k,] <- apply(Xk, 2, mean)
    #probabilité à priori d'être dans la classe k
    pik[k] <- nk / n
    
  }
  
  param$pkj <- pkj
  param$pik <- pik
  
  return(param)
}

binaryNBCval <- function(param, Xtst)
{
  #nb ligne n
  n <- nrow(Xtst)
  #nb col p
  p <- ncol(Xtst)
  #nb classes
  g <- length(param$pik)
  
  res <- list()
  prob <- matrix(0, nrow=n, ncol=g)
  
  for(k in 1:g)
  {
    pik <- param$pik[k]
    for(i in 1:n)
    {
      prob.ik <- 1
      for(j in 1:p)
      {
        pkj <- param$pkj[k,j]
        prob.ij <- pkj ** Xtst[i,j] * (1-pkj) ** (1 - Xtst[i,j])
        prob.ik <- prob.ik * prob.ij
      }
      prob[i,k] <- prob.ik * pik
    }
  }
  
  prob <- prob / apply(prob,1,sum) 
  pred <- max.col(prob)
  
  res$prob <- prob
  res$pred <- pred
  
  return(res)
  
}


calculTauxErreurNBC <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  fit <- binaryNBCfit(Xapp, zapp)
  val <- binaryNBCval(fit, Xtst)
  z.res <- val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}


scriptTauxErreurNBC <- function(donn, nbIter){
  errVec <- rep(0, nbIter)
  for(i in 1:nbIter){
    errVec[i] <- calculTauxErreurNBC(donn)
  }
  res <- list()
  res$errVec <- errVec
  res$moy <- mean(errVec)
  res$var <- var(errVec)
  
  return(res)
}

