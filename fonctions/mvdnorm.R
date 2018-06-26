library(MASS)
mvdnorm <- function(X, mu, Sigma)
{
  X <- as.matrix(X)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  dens <- tryCatch(
    {
      
      
      B <- chol(Sigma)
      U <- (X-matrix(rep(mu,n),nrow=n,byrow=T))%*%ginv(B)
      
      dens <- exp(-rowSums(U*U)/2) * (2*pi)^(-p/2) / det(B)
      
      
    },
    error=function(cond) {
      #message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      #message(cond)
      # Choose a return value in case of warning
      return(NA)
    },
    finally={
      
    }
  )
  return(dens) 
  
}

