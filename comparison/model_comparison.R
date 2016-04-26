

#--------------------#
#### Vorbereitung ####
#--------------------#

rm(list = ls())
setwd('/home/khusmann/mnt/H/additive_Biomassefunktionen')

## Daten einlesen und vorbereiten ## 

imp1 <- read.csv2('./linear_unrestricted/coef_ols.csv')
imp2 <- read.csv2('./linear_restriction/Coef_RLM.csv')
imp3 <- read.csv2('SUR/coef_SUR_sepGNLS.csv')

b <- cbind(imp1, imp2$parameter.linrest, imp3$parameter.sur, imp3$parameter.separateGNLS)
names(b) <- c('component', 'parameter.ols', 'parameter.linrest', 'parameter.sur', 'parameter.gnls')

rm(imp1, imp2, imp3)

d <- read.csv2('kompartimentsgewichte.csv')
d <- d[, c(2 : dim(d)[2])]
names(d)[5 : 8] <- c('as', 're', 'ho', 'ri')

load('./linear_unrestricted/lin_unrestr_models.RData')
load('./linear_restriction/lin_restr_model.RData')
load('./SUR/nonlin_models.RData')

# Speichern aller Koeffizienten in einer csv Datei
write.csv2(b, 'coef.csv')

## Funktionen deklarieren ## 
# Funktion zur Anwendung der Biomassefunktionen, Berechnung der Fitted Values
fit <- function(b0, b1, b2 = NA, dbh, h = NA){
  if(any(is.na(c(b2, h)))){
    exp(b0) * (dbh  ^ b1)
  }else{
    exp(b0) * (dbh  ^ b1) * h ^ b2
  }
}

## Berechnungen ##

# Fitted Values
d$fit.ho.ols <- fit(b0 = b$parameter.ols[1],     b1 = b$parameter.ols[2],     b2 = b$parameter.ols[3],     dbh = d$BHD, h = d$Hoehe)
d$fit.ho.rls <- fit(b0 = b$parameter.linrest[1], b1 = b$parameter.linrest[2], b2 = b$parameter.linrest[3], dbh = d$BHD, h = d$Hoehe)
d$fit.ho.sur <- fit(b0 = b$parameter.sur[1],     b1 = b$parameter.sur[2],     b2 = b$parameter.sur[3],     dbh = d$BHD, h = d$Hoehe)
d$fit.ho.gnl <- fit(b0 = b$parameter.gnls[1],    b1 = b$parameter.gnls[2],    b2 = b$parameter.gnls[3],    dbh = d$BHD, h = d$Hoehe)

d$fit.ri.ols <- fit(b0 = b$parameter.ols[4],     b1 = b$parameter.ols[5],     b2 = b$parameter.ols[6],     dbh = d$BHD, h = d$Hoehe)
d$fit.ri.rls <- fit(b0 = b$parameter.linrest[4], b1 = b$parameter.linrest[5], b2 = b$parameter.linrest[6], dbh = d$BHD, h = d$Hoehe)
d$fit.ri.sur <- fit(b0 = b$parameter.sur[4],     b1 = b$parameter.sur[5],     b2 = b$parameter.sur[6],     dbh = d$BHD, h = d$Hoehe)
d$fit.ri.gnl <- fit(b0 = b$parameter.gnls[4],    b1 = b$parameter.gnls[5],    b2 = b$parameter.gnls[6],    dbh = d$BHD, h = d$Hoehe)

d$fit.as.ols <- fit(b0 = b$parameter.ols[7],     b1 = b$parameter.ols[8],     dbh = d$BHD)
d$fit.as.rls <- fit(b0 = b$parameter.linrest[7], b1 = b$parameter.linrest[8], dbh = d$BHD)
d$fit.as.sur <- fit(b0 = b$parameter.sur[7],     b1 = b$parameter.sur[8],     dbh = d$BHD)
d$fit.as.gnl <- fit(b0 = b$parameter.gnls[7],    b1 = b$parameter.gnls[8],    dbh = d$BHD)

d$fit.re.ols <- fit(b0 = b$parameter.ols[9],     b1 = b$parameter.ols[10],     dbh = d$BHD)
d$fit.re.rls <- fit(b0 = b$parameter.linrest[9], b1 = b$parameter.linrest[10], dbh = d$BHD)
d$fit.re.sur <- fit(b0 = b$parameter.sur[9],     b1 = b$parameter.sur[10],     dbh = d$BHD)
d$fit.re.gnl <- fit(b0 = b$parameter.gnls[9],    b1 = b$parameter.gnls[10],    dbh = d$BHD)


## Histogramme: Observed - Fitted ##
setwd('/home/khusmann/mnt/H/additive_Biomassefunktionen/comparison/observed_fitted_hist/')
for(i in c('ho', 'ri', 'as', 're')){
  pdf(paste(i, 'pdf', sep= '.'))
  par(mfrow = c(2,2))
  for(j in c('ols', 'rls', 'gnl', 'sur')){
    hist(d[, i] - d[, paste('fit', i, j, sep = '.')], main = paste(i, j, sep = ' - '))
    mtext(side = 3, line = -0.5, paste('bias:', round(mean(d[, i] - d[, paste('fit', i, j, sep = '.')]), digits = 2)))
    mtext(side = 3, line = 0.5, paste('rmse:', sqrt( mean((d[, i] - d[, paste('fit', i, j, sep = '.')])^2)) ))
  }
  dev.off()
}
pdf('vb.pdf')
par(mfrow = c(2,2))
for(j in c('ols', 'rls', 'gnl', 'sur')){
  
  hist(apply(d[, 5 : 8], 1, sum) - apply(cbind(d[, paste('fit.ho', j, sep = '.')], d[, paste('fit.ri', j, sep = '.')], d[, paste('fit.as', j, sep = '.')], d[, paste('fit.re', j, sep = '.')]), 1, sum), main = paste(i, j, sep = ' - '))
  i_bias <- mean(apply(d[, 5 : 8], 1, sum) - apply(cbind(d[, paste('fit.ho', j, sep = '.')], d[, paste('fit.ri', j, sep = '.')], d[, paste('fit.as', j, sep = '.')], d[, paste('fit.re', j, sep = '.')]), 1, sum), main = paste(i, j, sep = ' - '))
  i_rmse <- sqrt(mean((apply(d[, 5 : 8], 1, sum) - apply(cbind(d[, paste('fit.ho', j, sep = '.')], d[, paste('fit.ri', j, sep = '.')], d[, paste('fit.as', j, sep = '.')], d[, paste('fit.re', j, sep = '.')]), 1, sum))^2))
  mtext(side = 3, line = -0.5, paste('bias:',round(i_bias, digits = 2)))
  mtext(side = 3, line =  0.5, paste('rmse:',round(i_rmse, digits = 2)))
}
dev.off()
