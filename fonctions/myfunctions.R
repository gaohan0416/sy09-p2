library(Rlab)
library(MASS)
library(rpart)

#calcul du taux d'erreur pour l'analyse discriminante quadratique (ADQ)
calculTauxErreurADQnoSep <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  param.adq <- adq.app(Xapp, zapp)
  adq.val <- ad.val(param.adq, Xtst)
  z.res <- adq.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

calculTauxErreurADQ <- function(Xapp, zapp, Xtst, ztst){
  param.adq <- adq.app(Xapp, zapp)
  adq.val <- ad.val(param.adq, Xtst)
  z.res <- adq.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

#calcul du taux d'erreur pour l'analyse discriminante lineaire (ADL)
calculTauxErreurADLnoSep <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  param.adl <- adl.app(Xapp, zapp)
  adl.val <- ad.val(param.adl, Xtst)
  z.res <- adl.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

calculTauxErreurADL <- function(Xapp, zapp, Xtst, ztst){
  param.adl <- adl.app(Xapp, zapp)
  adl.val <- ad.val(param.adl, Xtst)
  z.res <- adl.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

#calcul du taux d'erreur pour le classifieur bayésien naif (NBA)
calculTauxErreurNBAnoSep <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  param.nba <- nba.app(Xapp, zapp)
  nba.val <- ad.val(param.nba, Xtst)
  z.res <- nba.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

calculTauxErreurNBA <- function(Xapp, zapp, Xtst, ztst){
  param.nba <- nba.app(Xapp, zapp)
  nba.val <- ad.val(param.nba, Xtst)
  z.res <- nba.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}



#calcul du taux d'erreur pour la regression logistique (RL)
calculTauxErreurRLnoSep <- function(donn, intercept, epsi = 1*(10^-5)){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  res.app <- log.app(Xapp, zapp, intercept, epsi)
  
  res.val <- log.val(res.app$beta, Xtst)
  z.res <- res.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

calculTauxErreurRL <- function(Xapp, zapp, Xtst, ztst, intercept){
  res.app <- log.app(Xapp, zapp, intercept, 1*(10^-5))
  
  res.val <- log.val(res.app$beta, Xtst)
  z.res <- res.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

scriptCalculTauxErreurRL <- function(donn, intercept, nbIter){
  tauxErreur <- rep(0, nbIter)
  for(i in 1:nbIter){
    tauxErreur[i] <- calculTauxErreurRLnoSep(donn, intercept)
  }
  moy <- sum(tauxErreur) / nbIter
  res <- list(tauxErreur, moy)
  names(res) <- c("Vector", "Mean")
  
  return(res)
}

#calcul du taux erreur avec la méthode des arbres binaires de décison
calculTauxErreurABD <- function(Xapp, zapp, Xtst, ztst){
  
  #calcul de l'arbre de décision
  tree <- rpart(zapp~.,data=Xapp)
  #élagage de l'abre pour obtenir l'arbre optimal
  treeOptimal <- prune(tree,cp=tree$cptable[which.min(tree$cptable[,4]),1])
  
  #prediction de l'arbre
  pred <- predict(treeOptimal, Xtst)
  z.res <- round(pred)
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
  
}

calculTauxErreurRLQ <- function(Xapp, zapp, Xtst, ztst, intercept){
  donn <- transformDataForRLQ(Xapp, Xtst)
  Xapp2 <- donn$Xapp2
  Xtst2 <- donn$Xtst2
  
  rlq.app <- log.app(Xapp2, zapp, intercept, 1*(10^-5))
  rlq.val <- log.val(test.app$beta, Xtst2)
  z.res <- rlq.val$pred
  
  nbSim.tst <- length(which(z.res==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  return(tauxErreur.tst)
}

compareModelsError <- function(donn, nbIter, intercept = T){
  matErr <- matrix(0, nrow = nbIter, ncol = 5)
  i <- 0
  nbError <- 0
  while(i <= nbIter){
    n <- ncol(donn)
    X <- donn[,1:(n-1)]
    z <- donn[,n]
    donn.sep <- separ1(X, z)
    Xapp <- donn.sep$Xapp
    zapp <- donn.sep$zapp
    Xtst <- donn.sep$Xtst
    ztst <- donn.sep$ztst
    
    errADL = calculTauxErreurADL(Xapp, zapp, Xtst, ztst)
    errADQ = calculTauxErreurADQ(Xapp, zapp, Xtst, ztst)
    errNBA = calculTauxErreurNBA(Xapp, zapp, Xtst, ztst)
    errRL = calculTauxErreurRL(Xapp, zapp, Xtst, ztst, intercept)
    errABD = calculTauxErreurABD(Xapp, zapp, Xtst, ztst)
    
    if(errADL != 1 && errADQ != 1 && errNBA != 1 && errRL !=1 && errABD!=1){
      matErr[i,] <- c(errADL, errADQ, errNBA, errRL, errABD)
      i <- i + 1
    }
    else {
      #erreur
      nbError <- nbError + 1
    }
    
    
  }
  
  colnames(matErr) <- c("errADL", "errADQ", "errNBA", "errRL", "errABD")
  moy <- colSums(matErr) / nbIter
  var <- apply(matErr, 2, var)
  res <- NULL
  res$matErr <- matErr
  res$moy <- moy
  res$var <- var
  
  return(res)
}

compareModelsErrorOnSonar <- function(donn, nbIter){
  matErr <- matrix(0, nrow = nbIter, ncol = 4)
  i <- 0
  nbError <- 0
  while(i <= nbIter){
    n <- ncol(donn)
    X <- donn[,1:(n-1)]
    z <- donn[,n]
    donn.sep <- separ1(X, z)
    Xapp <- donn.sep$Xapp
    zapp <- donn.sep$zapp
    Xtst <- donn.sep$Xtst
    ztst <- donn.sep$ztst
    
    errADL = calculTauxErreurADL(Xapp, zapp, Xtst, ztst)
    errADQ = calculTauxErreurADQ(Xapp, zapp, Xtst, ztst)
    errNBA = calculTauxErreurNBA(Xapp, zapp, Xtst, ztst)
    errABD = calculTauxErreurABD(Xapp, zapp, Xtst, ztst)
    
    if(errADL != 1 && errADQ != 1 && errNBA != 1 && errABD!=1){
      matErr[i,] <- c(errADL, errADQ, errNBA, errABD)
      i <- i + 1
    }
    else {
      #erreur
      nbError <- nbError + 1
    }
    
    
  }
  
  colnames(matErr) <- c("errADL", "errADQ", "errNBA", "errABD")
  moy <- colSums(matErr) / nbIter
  var <- apply(matErr, 2, var)
  res <- NULL
  res$matErr <- matErr
  res$moy <- moy
  res$var <- var
  
  return(res)
}

#script pour tester les résultats de l'ADL, ADQ et le NBA un nombre nbIter de fois
compareModelsErrorClassDiscri <- function(donn, nbIter){
  matErr <- matrix(0, nrow = nbIter, ncol = 3)
  for(i in 1:nbIter){
    n <- ncol(donn)
    X <- donn[,1:(n-1)]
    z <- donn[,n]
    donn.sep <- separ1(X, z)
    Xapp <- donn.sep$Xapp
    zapp <- donn.sep$zapp
    Xtst <- donn.sep$Xtst
    ztst <- donn.sep$ztst
    
    errADL = calculTauxErreurADL(Xapp, zapp, Xtst, ztst)
    errADQ = calculTauxErreurADQ(Xapp, zapp, Xtst, ztst)
    errNBA = calculTauxErreurNBA(Xapp, zapp, Xtst, ztst)
    
    matErr[i,] <- c(errADL, errADQ, errNBA)
    
  }
  colnames(matErr) <- c("errADL", "errADQ", "errNBA")
  moy <- colSums(matErr) / nbIter
  
  res <- NULL
  res$matErr <- matErr
  res$moy <- moy
  
  return(res)
}
transformDataForRLQ <- function(Xapp, Xtst){
  Xapp2 <- Xapp
  Xtst2 <- Xtst
  for (p in 1:(dim(Xapp)[2]-1))
  {
    for (q in (p+1):dim(Xapp)[2])
    {
      Xapp2 <- cbind(Xapp2, Xapp[,p]*Xapp[,q])
      Xtst2 <- cbind(Xtst2, Xtst[,p]*Xtst[,q])
    }
  }
  for (p in 1:dim(Xapp)[2])
  {
    Xapp2 <- cbind(Xapp2, Xapp[,p]^2)
    Xtst2 <- cbind(Xtst2, Xtst[,p]^2)
  }
  res <- NULL
  res$Xapp2 <- Xapp2
  res$Xtst2 <- Xtst2
  return(res)
}


#fonction pour calculer LDA avec utilisation de fct de R
ldaR <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  lda.model <- lda(Xapp, grouping = zapp)
  
  predmodel.lda.app <- predict(lda.model, Xapp)
  predmodel.lda.tst <- predict(lda.model, Xtst)
  
  table.lda.app <- table(Predicted=predmodel.lda.app$class, z=zapp)
  table.lda.tst <- table(Predicted=predmodel.lda.tst$class, z=ztst)
  accuracy.app <- (table.lda.app[1,1] + table.lda.app[2,2])/sum(colSums(table.lda.app))
  accuracy.tst <-  (table.lda.tst[1,1] + table.lda.tst[2,2])/sum(colSums(table.lda.tst))
  
  z.res.app <- predmodel.lda.app$class
  nbSim.app <- length(which(z.res.app==zapp))
  tauxErreur.app <- 1 - nbSim.app/length(zapp)
  z.res.tst <- predmodel.lda.tst$class
  nbSim.tst <- length(which(z.res.tst==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  err <- list(tauxErreur.app, tauxErreur.tst)
  names(err) <- c("Training","Test")
  
  acc <- list(accuracy.app, accuracy.tst)
  names(acc) <- c("Training","Test")
  
  pred <- list(predmodel.lda.app,predmodel.lda.tst)
  names(pred) <- c("Training", "Test")
  
  res <- list(lda.model, pred, acc, err)
  names(res) <- c("Modele", "Prediction", "Accuracy", "Error")
  
  return (res)
}

#fonction pour calculer LDA avec utilisation de fct de R
qdaR <- function(donn){
  n <- ncol(donn)
  X <- donn[,1:(n-1)]
  z <- donn[,n]
  donn.sep <- separ1(X, z)
  Xapp <- donn.sep$Xapp
  zapp <- donn.sep$zapp
  Xtst <- donn.sep$Xtst
  ztst <- donn.sep$ztst
  
  qda.model <- qda(Xapp, grouping = zapp)
  
  predmodel.qda.app <- predict(qda.model, Xapp)
  predmodel.qda.tst <- predict(qda.model, Xtst)
  
  table.qda.app <- table(Predicted=predmodel.qda.app$class, z=zapp)
  table.qda.tst <- table(Predicted=predmodel.qda.tst$class, z=ztst)
  accuracy.app <- (table.qda.app[1,1] + table.qda.app[2,2])/sum(colSums(table.qda.app))
  accuracy.tst <-  (table.qda.tst[1,1] + table.qda.tst[2,2])/sum(colSums(table.qda.tst))
  
  z.res.app <- predmodel.qda.app$class
  nbSim.app <- length(which(z.res.app==zapp))
  tauxErreur.app <- 1 - nbSim.app/length(zapp)
  z.res.tst <- predmodel.qda.tst$class
  nbSim.tst <- length(which(z.res.tst==ztst))
  tauxErreur.tst <- 1 - nbSim.tst/length(ztst)
  
  err <- list(tauxErreur.app, tauxErreur.tst)
  names(err) <- c("Training","Test")
  
  acc <- list(accuracy.app, accuracy.tst)
  names(acc) <- c("Training","Test")
  
  pred <- list(predmodel.qda.app,predmodel.qda.tst)
  names(pred) <- c("Training", "Test")
  
  res <- list(qda.model, pred, acc, err)
  names(res) <- c("Modele", "Prediction", "Accuracy", "Error")
  
  return (res)
}
