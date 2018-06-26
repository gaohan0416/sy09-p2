adq.app <- function(Xapp, zapp)
{
  n <- dim(Xapp)[1]
  p <- dim(Xapp)[2]
  g <- max(unique(zapp))
  
  param <- NULL
  param$MCov <- array(0, c(p,p,g))
  param$mean <- array(0, c(g,p))
  param$prop <- rep(0, g)
  
  for (k in 1:g)
  {
    indk <- which(zapp==k)
    #individus de la classe k
    Xk <- Xapp[indk,]
    #nb d'individus dans la classe k
    nk <- length(indk)
    
    param$MCov[,,k] <- cov(Xk)
    param$mean[k,] <- apply(Xk, 2, mean)
    param$prop[k] <- nk / n
  }
  
  return(param)
}

adl.app <- function(Xapp, zapp)
{
  n <- dim(Xapp)[1]
  p <- dim(Xapp)[2]
  g <- max(unique(zapp))
  
  param <- NULL
  MCov <- array(0, c(p,p))
  param$MCov <- array(0, c(p,p,g))
  param$mean <- array(0, c(g,p))
  param$prop <- rep(0, g)
  
  for (k in 1:g)
  {
    indk <- which(zapp==k)
    #individus de la classe k
    Xk <- Xapp[indk,]
    #nb d'individus dans la classe k
    nk <- length(indk)
    
    MCov <- MCov + nk * cov(Xk)
    param$mean[k,] <- apply(Xk, 2, mean)
    param$prop[k] <- nk / n
  }
  MCov <- MCov / n
  
  for (k in 1:g)
  {
    param$MCov[,,k] <- MCov
  }
  
  return(param)
}

nba.app <- function(Xapp, zapp)
{
  n <- dim(Xapp)[1]
  p <- dim(Xapp)[2]
  g <- max(unique(zapp))
  
  param <- NULL
  
  param$MCov <- array(0, c(p,p,g))
  param$mean <- array(0, c(g,p))
  param$prop <- rep(0, g)
  
  for (k in 1:g)
  {
    indk <- which(zapp==k)
    #individus de la classe k
    Xk <- Xapp[indk,]
    #nb d'individus dans la classe k
    nk <- length(indk)
    
    MCov <- array(0, c(p,p))
    diag(MCov) <- diag(cov(Xk))
    param$MCov[,,k] <- MCov
    param$mean[k,] <- apply(Xk, 2, mean)
    param$prop[k] <- nk / n
  }
  
  return(param)
}

ad.val <- function(param, Xtst)
{
  n <- dim(Xtst)[1]
  p <- dim(Xtst)[2]
  g <- length(param$prop)
  
  out <- NULL
  
  prob <- matrix(0, nrow=n, ncol=g)
  
  for (k in 1:g)
  {
    # loi de densité pour la classe k
    fk <- mvdnorm(Xtst, param$mean[k,], param$MCov[,,k])
    #on multiplie par la probabilité à priori pik
    prob[,k] <- param$prop[k] * fk
  }
  
  #pour obtenir la probabilité à postériori on divise par la somme des nk * fk
  prob <- prob / apply(prob,1,sum) 
  pred <- max.col(prob)
  
  out$prob <- prob
  out$pred <- pred
  
  return(out)
}
