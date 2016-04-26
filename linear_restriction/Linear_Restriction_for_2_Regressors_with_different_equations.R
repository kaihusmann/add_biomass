#------------------------------------------------------------------------------#
# Estimation with linear Restriction for 2 Regressors with different equations #
#------------------------------------------------------------------------------#

rlm_dh2 <- function(komp, BHD, Height, Vollbaum, n_boot = 10000){
  # Eingangsvariablen:
  # komp: Vektor mit log Kompartimentsgewichten. Die Kompartimente muessen in 
  # der Reihenfolge: Ast, Reisig, Derbholz, Derbholzrinde untereinander
  #
  # bhd: Vektor mit log BHD
  # hoehe: Vektor mit log hoehe
  # Vollbaum: Vektor mit log addierten Kompartimentsgewichten
  
  ## Laden zusaetzlicher Pakete und Variablendeklaration ##
  
  nlength <- length(komp) / 4
  
  x.11 <- kronecker(diag(4), rep(1, nlength))
  x.22 <- kronecker(diag(4), BHD)
  x.33 <- kronecker(diag(4), Height)
  XX <- cbind(x.11, x.22, x.33[,c(3,4)])
  
  # Y matrix
  YY <- komp
  
  # Regression unrestringiert
  beta_ur <- solve(t(XX)%*%XX) %*% t(XX) %*% YY
  
  # Parameter Bootstrap matrx
  Parameter <- matrix(NA, nrow = n_boot, ncol = 10)
  
  # Covariance Bootstrap array
  Covariance.Matrix <- array(NA, c(10,10, n_boot))
  
  #set.seed(1234)
  for(i in 1 : n_boot){
    Z <- sample(c(1 : nlength), 3, replace = FALSE)
    v <- Vollbaum[c(Z)]
    a.1 <- matrix(1, nrow = 3, ncol = 4)
    a.2 <- matrix(BHD[c(Z)], nrow = 3, ncol = 4)
    a.3 <- matrix(Height[c(Z)], nrow = 3, ncol = 2)
    A <- cbind(a.1, a.2, a.3)
    #Regression mit Restriktion
    Beta_R <- NA
    try(
      Beta_R <- beta_ur - solve(t(XX)%*%XX) %*% t(A) %*% solve(A %*% solve(t(XX)%*%XX) %*% t(A)) %*% (A %*% beta_ur - v),
      silent = TRUE
    )
    if(any(!is.na(Beta_R))){
      if(any( abs(Beta_R[5 : 10]) > abs(beta_ur[5 : 10]) * 2) |  any(Beta_R[5:10] < 0) ){
        Beta_R <- NA
      }
    }
    if(!any(is.na(Beta_R))){
      Parameter[i,] <- Beta_R 
      Covariance.Matrix[,,i] <- solve(t(XX)%*%XX) -  solve(t(XX)%*%XX) %*% t(A) %*% solve(A %*% solve(t(XX)%*%XX) %*% t(A)) %*% A %*% solve(t(XX)%*%XX)
      RSS <- t(YY - XX %*% Beta_R) %*% (YY - XX %*% Beta_R)
      S.2 <- as.numeric(RSS/(nrow(XX) - ncol(XX)))
      Covariance.Matrix[,,i] <-  S.2 * Covariance.Matrix[,,i]
    }
    
  }
  # Parameter und Covariance Matrix reduzieren = NA Zeilen entfernen
  Covariance.Matrix <- Covariance.Matrix[, , !(is.na(Parameter[, 1]))]
  Parameter <- Parameter[ !(is.na(Parameter[, 1])), ]
  r <- n_boot/nrow(Parameter)
  
  # Mediane ermitteln (robuster gegen?ber ausrei?ern zum mean)
  a1 <- median(Parameter[,1])
  a2 <- median(Parameter[,2])
  a3 <- median(Parameter[,3])
  a4 <- median(Parameter[,4])
  
  b1 <- median(Parameter[,5])
  b2 <- median(Parameter[,6])
  b3 <- median(Parameter[,7])
  b4 <- median(Parameter[,8])
  
  g1 <- median(Parameter[,9])
  g2 <- median(Parameter[,10])
  
  # Fuer Cov
  Covariance_return <- apply(Covariance.Matrix, c(1,2), median)
  Covariance_return <- round(Covariance_return, digits = 6)
  
  # Rueckgabe
  param_return <- c(a1, a2, a3, a4, b1, b2, b3, b4, g1, g2)
  
  # Output erstellen
  output <- list(point_estimates = param_return,
                 covariance      = Covariance_return,
                 success         = r                 )
  return(output)
  
}