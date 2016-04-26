rm(list = ls())

# GITHUB

## Pakete laden und Daten vorbereiten ##
library(ggplot2)
library(reshape2)

setwd('H:\\Lange\\linear_restriction')
#setwd('/home/khusmann/mnt/H/additive_Biomassefunktionen/')

# Funktionen zur Linearen Restringierten Regression laden:
# 1.: Schaetzung mit BHD
source('Linear_Restriction_with_DBH_Function.R')
# 2.: Schaetzung mit BHD und Hoehe
source('Linear_Restriction_with_DBH_and_Height_Function.R')
# 3.: Schaetzung mit BHD und Hoehe in Derbholz und Rinde
source('Linear_Restriction_for_2_Regressors_with_different_equations.R')

# Datensaete laden
d <- read.csv2("kompartimentsgewichte.csv")

# logarithmieren
bhd <- log(d[, 3])
h <- log(d[, 4])

 # centering and scaling
# mb <- mean(bhd)
# bhd <- (bhd - mb)/sd(bhd)
# mh <- mean(h)
# h <- (h - mh)/sd(h)

# Vektor der abhängigen Variablen (Biomasse)
komp <- log(d[,c(6 : 9)])
Y <- melt(komp)
Y <- Y[,2]


## Aufruf der Regression ##
#para <- rlm_dh(komp = Y, BHD = bhd, Height = h, Vollbaum = rowSums(log(d[,c(6 : 9)])), n_boot = 10000)
para2 <- rlm_dh2(komp = Y, BHD = bhd, Height = h, Vollbaum = rowSums(log(d[,c(6 : 9)])), n_boot = 100000)
# Anzeige der Parameter
#para
para2
RLM.Parameters <- rbind(para2[3], para2[7], para2[9], para2[4], para2[8], para2[10], para2[1], para2[5]
                        , para2[2], para2[6])


## Datenausgabe ##
ret <- data.frame(component = c(
  'Derbholz.b0', 'Derbholz.b1', 'Derbholz.b2', 
  'Derbholzrinde.b0', 'Derbholzrinde.b1', 'Derbholzrinde.b2', 
  'Ast.b0', 'Ast.b1', 
  'Reisig.b0', 'Reisig.b1'
),
parameter.linrest = RLM.Parameters
)

write.csv2(ret, 'Coef_RLM.csv')


# LM zur Schaetzung der Signifikanz der Parameter
summary(lm(log(d$ast) ~ log(d$BHD) + log(d$Hoehe)))
summary(lm(log(d$reisig) ~ log(d$BHD) + log(d$Hoehe)))
summary(lm(log(d$derbholz) ~ log(d$BHD) + log(d$Hoehe)))
summary(lm(log(d$derbholzrinde) ~ log(d$BHD) + log(d$Hoehe)))

summary(lm(log(d$derbholzrinde) ~ log(d$BHD) ))
summary(lm(log(d$derbholzrinde) ~ log(d$BHD) ))


## Erstellen von Abb. und Testen der Modellergebnisse ##

# Illsutration der Sch?tzung der Kompartimente
# Derbholz
D.plot.boot <- function(x,y){
  z <-   para[3] + para[7]*x + para[11]*y   
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), D.plot.boot)
Dplot <- persp(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$derbholz), pmat = Dplot), col = 'black', pch = 16)

# Derbholz, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para[3]) * x^para[7] * y^para[11] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), d$derbholz, pmat = Dplot), col = 'black', pch = 16)

ges_derbholz <- z.werte

# Derbholzrinde
DR.plot.boot <- function(x,y){
  z <-   para[4] + para[8]*x + para[12]*y
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), DR.plot.boot)
DRplot <- persp(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$derbholzrinde), pmat = DRplot), col = 'black', pch = 16)



# Derbholzrinde, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para[4]) * x^para[8] * y^para[12] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), d$derbholzrinde, pmat = Dplot), col = 'black', pch = 16)

ges_derbrinde <- z.werte

# ?ste
A.plot.boot <- function(x,y){
  z <-   para[1] + para[5]*x + para[9]*y
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), A.plot.boot)
Aplot <- persp(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$ast), pmat = Aplot), col = 'black', pch = 16)


# Aeste, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para[1]) * x^para[5] * y^para[9] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), d$ast, pmat = Dplot), col = 'black', pch = 16)

ges_aeste <- z.werte

# Reisig
R.plot.boot <- function(x,y){
  z <-   para[2] + para[6]*x + para[10]*y 
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), R.plot.boot)
Rplot <- persp(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$reisig), pmat = Rplot), col = 'black', pch = 16)

# Reisig, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para[2]) * x^para[6] * y^para[10] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 90, phi = 0)
points(trans3d(exp(bhd), exp(h), d$reisig, pmat = Dplot), col = 'black', pch = 16)

ges_reisig <- z.werte

# Vollbaumschtzung
VB.plot.boot <- function(x,y){
  z <-  para[1] + para[2] + para[3] + para[4] + para[5]*x + para[6]*x + para[7]*x + para[8]*x + para[9]*y + para[10]*y + para[11]*y + para[12]*y  
  return(z)
}

z.werte <- outer(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), VB.plot.boot)
Baumplot <- persp(seq(min(bhd),max(bhd), 0.1), seq(min(h), max(h), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$vollbaum), pmat = Baumplot), col = 'black', pch = 16)


# Vollbaumschtzung, ruecktransformiert
VB.plot.boot <- function(x,y){
  #z <-  exp(para[1] + para[2] + para[3] + para[4] + para[5]*log(x) + para[6]*log(x) + para[7]*log(x) + para[8]*log(x) + para[9]*log(y) + para[10]*log(y) + para[11]*log(y) + para[12]*log(y))
  z <- (exp(para[1]) * x^para[5] * y^para[9]) + (exp(para[2]) * x^para[6] * y^para[10]) + (exp(para[3]) * x^para[7] * y^para[11]) + (exp(para[4]) * x^para[8] * y^para[12])
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), VB.plot.boot)
Baumplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), 0.1), seq(min(exp(h)), max(exp(h)), 0.1), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 45, phi = 0)
points(trans3d(exp(bhd), exp(h), (d$vollbaum), pmat = Baumplot), col = 'black', pch = 16)

hist(d$vollbaum - VB.plot.boot(exp(bhd), exp(h)))

#-------------------#
# Plots f?r rlm_dh2 #
#-------------------#

# Derbholz
D.plot.boot <- function(x,y){
  z <-   para2[3] + para2[7]*x + para2[9]*y   
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), D.plot.boot)
Dplot <- persp(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$derbholz), pmat = Dplot), col = 'black', pch = 16)

# Derbholz, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para2[3]) * x^para2[7] * y^para2[9] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), d$derbholz, pmat = Dplot), col = 'black', pch = 16)

# Derbholzrinde
DR.plot.boot <- function(x,y){
  z <-   para2[4] + para2[8]*x + para2[10]*y
  return(z)
}
z.werte <- outer(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), DR.plot.boot)
DRplot <- persp(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$derbholzrinde), pmat = DRplot), col = 'black', pch = 16)



# Derbholzrinde, ruecktransformiert
D.plot.boot <- function(x, y){
  z <-   exp(para2[4]) * x^para2[8] * y^para2[10] 
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), D.plot.boot)
Dplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), d$derbholzrinde, pmat = Dplot), col = 'black', pch = 16)

# Aeste
theme_set(theme_bw(15))
Ast <- para2[1] + para2[5]*log(d$BHD)
AA <- as.data.frame(cbind(Ast, log(d$ast), log(d$BHD)))
names(AA) <- c('fit', 'true', 'BHD')
ggplot(AA, aes(x = BHD, y = fit)) + geom_line(size = 1.2) + geom_point(aes(x = BHD, y = true), size = 2) + ylab('ln Biomasse') +
  xlab('ln BHD') + ggtitle('?ste')

# Aeste, ruecktransformiert
Ast <- exp(para2[1]) * d$BHD^para2[5]
AA <- as.data.frame(cbind(Ast, d$ast, d$BHD))
names(AA) <- c('fit', 'true', 'BHD')
ggplot(AA, aes(x = BHD, y = fit)) + geom_line(size = 1.2) + geom_point(aes(x = BHD, y = true), size = 2) + ylab('Biomasse') +
  xlab('BHD') + ggtitle('?ste')


# Reisig
Reisig <- para2[2] + para2[6]*log(d$BHD)
RR <- as.data.frame(cbind(Reisig, log(d$reisig), log(d$BHD)))
names(RR) <- c('fit', 'true', 'BHD')
ggplot(RR, aes(x = BHD, y = fit)) + geom_line(size = 1.2) + geom_point(aes(x = BHD, y = true), size = 2) + ylab('ln Biomasse') +
  xlab('ln BHD') + ggtitle('Reisig')

# Reisig, ruecktransformiert
Reisig <- exp(para2[2]) * d$BHD^para2[6]
RR <- as.data.frame(cbind(Reisig, d$reisig, d$BHD))
names(RR) <- c('fit', 'true', 'BHD')
ggplot(RR, aes(x = BHD, y = fit)) + geom_line(size = 1.2) + geom_point(aes(x = BHD, y = true), size = 2) + ylab('Biomasse') +
  xlab('BHD') + ggtitle('Reisig')

# Vollbaumschtzung
VB.plot.boot <- function(x,y){
  z <-  para2[1] + para2[2] + para2[3] + para2[4] + para2[5]*x + para2[6]*x + para2[7]*x + para2[8]*x + para2[9]*y + para2[10]*y 
  return(z)
}

z.werte <- outer(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), VB.plot.boot)
Baumplot <- persp(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(bhd, h, log(d$vollbaum), pmat = Baumplot), col = 'black', pch = 16)


# Vollbaumschtzung, ruecktransformiert
VB.plot.boot <- function(x,y){
  z <- exp(para2[3]) * x^para2[7] * y^para2[9] + exp(para2[4]) * x^para2[8] * y^para2[10] + exp(para2[1]) * x^para2[5] + exp(para2[2]) * x^para2[6]
  return(z)
}

z.werte <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)),length.out = 30), VB.plot.boot)
Baumplot <- persp(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), z.werte, xlab = 'DBH', ylab = 'Height', zlab = 'Biomass', expand = 0.45, theta = 0, phi = 0)
points(trans3d(exp(bhd), exp(h), (d$vollbaum), pmat = Baumplot), col = 'black', pch = 16)

hist(d$vollbaum - VB.plot.boot(exp(bhd), exp(h)))

