#-------------------------------#
# restricted linear estimation  #
#-------------------------------#

rlm <- function(komp, bhd){
  # Eingangsvariablen:
  # komp: Data Frame mit Kompartimentsgewichten. Die Kompartimente muessen in 
  # der Reihenfolge: Ast, Reisig, Derbholz, Derbholzrinde uebergeben werden
  
  # bhd: Vektor mit BHD
  
  ## Laden zusaetzlicher Pakete und Variablendeklaration ##
  require(reshape)
  
  ## Variablendeklaration und Berechnungen
  # # Speichern von Rohdaten Vollbaum
  # VB <- apply(komp, 1, sum)
  nlength <- dim(komp)[1]
  # 
  # #d[,10] <- d[, 6] * d[, 7] * d[, 8] * d[, 9]
  # 
  # # log linearisieren der Holzmessungen
  # Holz <- log(komp)
  # 
  # # transormation der Vollbaum variable: Y = A^a*R^r...
  Baum <- komp[, 1] * komp[, 2] * komp[ ,3] * komp[ ,4]
  # 
  # restriction 
  # reisig <- Holz[,2]
  # Baum <- log(d[,10])
  # reisig <- reisig - Baum
  
  # Abhaengige Variable Y erstellen
  # Holz <- Holz[,-2]
  # Holz[,4] <- reisig
  # Y <- melt(Holz)
  # Y <- Y[,2]
  # 
  # Design matrix
  #x.1 <- rep(1, nlength)
  #X <- kronecker(diag(3), x.1)
  #x.2 <- matrix(-1, nrow = nlength, ncol = 3)
  #X <- rbind(X, x.2)
  
  # Regression
  #beta <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # berechnen von reisig intercept
  #beta_r <- mean(Baum) - beta[1] - beta[2] - beta[3] 
  
  # zusammenf?gen der regressionsparameter
  #beta <- rbind(beta, beta_r)
  
  # testen ob Summe der Parameter den Vollbaum ergeben
  #sum(beta) == mean(Baum)
  
  # Ruecktransformation der Parameter
  #beta <- exp(beta)
  #rownames(beta) <- c('?ste', 'Derbholz', 'Rinde', 'Reisig')
  
  
  # Design matrix erstellen
  BHD <- log(bhd)
  x.11 <- kronecker(diag(4), rep(1, nlength))
  x.22 <- kronecker(diag(4), BHD)
  XX <- cbind(x.11, x.22) 
  
  # Y matrix
  Holz2 <- log(d[,6:9])
  YY <- melt(Holz2)
  YY <- YY[,2]
  # 
  # # erstellen der Resktriktion
  # Vollbaum <- log(d[,10])
  # set.seed(121684)
  # Z <- sample(c(1:35), 2)
  # Vollbaum <- Vollbaum[c(Z)]
  # a.1 <- matrix(1, nrow = 2, ncol = 4)
  # a.2 <- matrix(BHD[c(Z)], nrow = 2, ncol = 4)
  # A <- cbind(a.1, a.2)
  # 
  # # Regression unrestringiert
  # beta_ur <- solve(t(XX)%*%XX) %*% t(XX) %*% YY
  # 
  # # Regression mit Restriktion
  # Beta_R <- beta_ur - solve(t(XX)%*%XX) %*% t(A) %*% solve(A %*% solve(t(XX)%*%XX) %*% t(A)) %*% (A %*% beta_ur - Vollbaum)
  # 
  # summe_geschaetzt <- Beta_R[3,] + Beta_R[7,] * log(d$BHD) +
  #   Beta_R[1,] + Beta_R[5,] * log(d$BHD) +
  #   Beta_R[2,] + Beta_R[6,] * log(d$BHD) +
  #   Beta_R[4,] + Beta_R[8,] * log(d$BHD)
  # 
  # plot(summe_geschaetzt ~ log(d$BHD), ylim = c(-5, 40))
  # points(log(d$vollbaum) ~ log(d$BHD), pch = 2)
  
  #---------------#
  # Bootstrapping #
  #---------------#
  
  # Regression unrestringiert
  beta_ur <- solve(t(XX)%*%XX) %*% t(XX) %*% YY
  Fitted <- matrix(0, nrow = 10000, ncol = nlength)
  Parameter <- matrix(0, nrow = 10000, ncol = 8)
  #set.seed(1)
  for(i in 1:10000){
    Vollbaum <- log(komp[, 1] * komp[, 2] * komp[ ,3] * komp[ ,4])
    Z <- sample(c(1:35), 2, replace = FALSE)
    Vollbaum <- Vollbaum[c(Z)]
    a.1 <- matrix(1, nrow = 2, ncol = 4)
    a.2 <- matrix(BHD[c(Z)], nrow = 2, ncol = 4)
    A <- cbind(a.1, a.2)
    #Regression mit Restriktion
    try(
      Beta_R <- beta_ur - solve(t(XX)%*%XX) %*% t(A) %*% solve(A %*% solve(t(XX)%*%XX) %*% t(A)) %*% (A %*% beta_ur - Vollbaum)
    )
    
    Fitted[i,] <- Beta_R[3,] + Beta_R[7,] * log(d$BHD) +
      Beta_R[1,] + Beta_R[5,] * log(d$BHD) +
      Beta_R[2,] + Beta_R[6,] * log(d$BHD) +
      Beta_R[4,] + Beta_R[8,] * log(d$BHD)
    Parameter[i,] <- Beta_R
  }
  
  return(apply(Parameter, 2, median))
}